// Steady solves: the state a transient from the present one ends in, found without waiting it
// out, by Newton on the steady equations where it reaches, by continuation along the drives
// where it doesn't, and by pseudo-transient time steps where nothing else does.
//
// CONTRACT
//
// Methods of the Solver (src/solver.js; its layout is in docs/internals.md), installed on its
// prototype like the assembly's. Device.solve and Device.impedance (src/index.js) call
// solveSteady, from a snapshot that Device.solve restores where it fails; time steps
// (src/timeSteps.js) read fastestTime; Solution and the checks read stuckReactions. What they
// use of the rest:
//   this.step(dt, opts)  one time step (src/timeSteps.js), whose contract this leans on: at
//                        dt = ∞ with this.steady, the steady equations with the laws' rows in
//                        place (src/lawRows.js), sources at this.steadyAt; accepted, the time
//                        moves on by dt and the history gains the step; failed, all as before
//                        it.
//   captureLaws, captureLocal (src/lawRows.js)  the amounts and charges the laws hold, taken
//                        from the state when a solve begins; laws.js's findLaws (through
//                        this._laws, on a change of drive) says which laws there are, and
//                        this.unkept which no row can keep.
//   this.redrive(), this._refreshSources(), this.sourceOverride  a terminal's drive, changed
//                        and put back; the held voltages read at this.sourceTime, or the
//                        override where one is set. A redrive() that changes how a terminal is
//                        driven (a voltage, a current) finds the laws again (this._laws): new
//                        constraints, local constraints with their references at 0, islands
//                        off, no pins; this.atSteady false, the history empty.
//   this.initFromComposition() (src/start.js), this._snapshot()/_restore(), this.amount(st),
//   this.computeConcentrations(), maxPotentialStep (src/newton.js).
//
// solveSteady(opts) → { converged, steps, iterations, history, residual?, ramped? }
//   opts: maxSteps (80: time steps per pseudo-transient ramp), tol (1e-11: the largest potential
//   change, in thermal units, that ends a ramp of giant steps), continuation (false: Newton from
//   here and the pseudo-transient ramp only). Sets this.violent (violentStart, where nothing has
//   been stepped or solved), clears this.unreached, this.charging and this.drifted, captures the
//   laws, and
//   solves with this.steady = true and this.steadyAt = this.time (false again on return).
//   Stops at once, converged false, steps 0, where a floating terminal is driven by a nonzero
//   current that only charges a capacitance or feeds a conserved combination (a constraint):
//   this.charging says which terminal and what it fills. Else _steadyFromHere, and
//   this.atSteady = converged. Not converged: this.solvedV and this.stepped as before the call.
//
// _steadyFromHere(opts): in stages until one converges, each from the state the solve began in
// (u, uLo and termV, put back between stages, the concentrations recomputed) unless it says
// otherwise. A continuation applies with continuation not off; canContinue: both contacts held
// at voltages (not behind a resistance) that differ.
//   1. Newton from here (_steadySteps): a single step (maxSteps 1) where canContinue or there's
//      generation, and the device can be solved directly (_directSteady); else the
//      pseudo-transient ramp with opts' maxSteps. Failing on a state reached by time steps
//      (this.stepped), on a direct device with continuation not off, the same again from a cold
//      start (initFromComposition), whose u and uLo become the stages' start (termV stays the
//      start's: a floating terminal keeps the voltage the transient left it at).
//   2. Continuation: generation ramped up (_generationContinuation); a terminal driven by a
//      current, not bound (_bound), held at voltages instead (_floatingContinuation); with
//      canContinue, the right contact's voltage ramped to its target (_continuation): from the
//      voltages of the last converged solve (this.solvedV), the left contact's first where it
//      moved, else from level contacts, a cold start.
//   3. The pseudo-transient ramp: after the contact's ramp fails, from the start; otherwise only
//      where stage 1 was a single step, from the dimmer light's solution where the generation
//      ramp got past its first solve (u, uLo and termV as that ramp left them), else from the
//      start.
//   Converged (but by _floatingContinuation): this.solvedV = [termV[0], termV[1]], the contacts'
//   own voltages at the solution (driven by a current or behind a resistance too), this.stepped =
//   false. (By _floatingContinuation: as its last held solve set them.)
//
// _steadySteps({ maxSteps, tol }): up to maxSteps time steps (tau = slowestTime), dt from giant
//   (∞ where the device can be solved directly; else 1e6·tau): a first step at ∞ that fails is
//   retried once with maxStep 3. A step that fails sends dt down (to tau·1e-6 from giant, then
//   ÷4, stopping below min(1e-15·tau, 1e-2·fastestTime)); one that converges sends it back up
//   (×2 after a failure, growing by ×1.5 to ×10 a step; past 1e6·tau, giant); on a direct device,
//   once dt ≥ tau, a step at ∞ is tried once a decade of dt, the state kept where it fails. Past
//   giant (not direct), dt grows ×10 a step to 1e6·giant. Converged where a step at ∞ converged,
//   or, not direct, where a step of dt ≥ giant moved the potentials under tol (finished by one
//   step of exactly giant where dt had grown past it, which must converge too). Not direct,
//   before each step of dt ≥ giant and the finishing one, each spectator's amount restored
//   (_renormalizeSpectators). Not direct, with unkept laws: after each step of dt ≥ giant, each
//   unkept combination of stretches' and surfaces' amounts (_unkeptAmounts) must be within 1e-6
//   of its size (Σ |w·amount|) of what it was when the call began, else the step is put back
//   (_restore) and counts as failed (`drifted` in its history entry; this.drifted names the
//   combination's parts); the finishing step, failing that, ends the call unconverged. (A
//   combination through an electrode's or a contact's charge isn't checked.) On return, whatever
//   the outcome: this.time as it was, the history
//   empty, this.dtNext undefined; converged, every connected stretch's referenceAmounts its
//   amount now and its boundaryIntake 0.
//
// The continuations (each warm from the last converged solve; a solve that fails puts back
// u, uLo and termV and takes a smaller step):
//   _generationContinuation  generationScale from 1e-30 up to 1, by a factor of 100 at most
//                            (squared after a success, its square root after a failure, giving
//                            up below 1.5); fails if the first, at 1e-30, does; this.dimmer once
//                            it converges. generationScale is 1 again on return.
//   _continuation(k, from, to)  contact k's voltage through sourceOverride (the other's held
//                            too, where given), in steps of an eighth of the span, ×1.5 on
//                            success, ÷4 on a failure, giving up below 1e-6 of the span. Cold
//                            (not warm), first solved at `from`: failing, again from a cold start
//                            laid out there, then with the light ramped where there's
//                            generation. The overrides are cleared on return.
//   _floatingContinuation    the terminal held at voltages (its drive replaced, then put back):
//                            marched from its voltage (or 0) in steps from 20 mV, ×1.5, toward
//                            the target current (both ways, to ±20 V), the crossing bisected to
//                            0.1 mV, then floated from the nearer end. Floated and failing, the
//                            held solve whose current is within 1e-10 of the target (by regula
//                            falsi, Illinois, from the bracket; the scale is the largest current
//                            seen) is returned as converged, its drive put back. No crossing:
//                            this.unreached says what the held voltages passed. Each time its
//                            drive is put back, the local constraints' references are taken from
//                            the state again (captureLocal), which redrive() zeroed.
//
// Also: slowestTime() = L²/D_min over the species present with D > 0 in the materials some
// region is made of (1 where none), fastestTime() = h_min²/D_max over the same (slowestTime where
// there's no such cell); _directSteady(): no law unkept and every
// stretch connected or conserved; stuckReactions(): after a steady solve, each reaction that
// alone changes a combination nothing feeds and whose affinity keeps one sign, over 1e-6, wherever
// it acts (see its comment).
//
// Invariants (each says what a reviewer should be able to check here):
//   1. Converged: the state solves the steady equations as step(∞) solved them; or, not direct,
//      a giant step's (whose storage term is 1e-6/tau of the rest at most); or, by
//      _floatingContinuation's regula falsi, the held equations at a voltage whose current is
//      within 1e-10 of the target's scale.
//   2. On return: this.time, every terminal's drive, sourceOverride, generationScale and
//      this.steady as before the call; the history empty, but where the call stopped at once
//      (this.charging) before any step.
//   3. A stage that fails leaves no trace in the next one's start: u, uLo and termV are put back
//      (the concentrations recomputed); only stage 1's cold start, on a stepped state, moves the
//      start, and stage 3 may begin from the generation ramp's solution.
//   4. Deterministic: the same state and options give the same result, bit for bit. The state
//      here is more than u, uLo, termV and time: this.solvedV, this.stepped, whether the history
//      is empty (this.violent) and the drives, all of which choose the stages.
//   5. Reach: time steps reach a mode as slow as ~1e12·tau (1e6·tau a step, past the giant; on a
//      direct device no finite step is longer than 1e6·tau). A reaction that alone settles a
//      closed combination, slower than that, isn't waited out: the solve fails or converges where
//      it began, and stuckReactions names it (Solution's warnings, both ways).
// Tests: every Device.solve goes through here; test/steady.test.js (the stages, the laws, the
// continuations) and test/generation.test.js; the devices that needed each stage are named in its
// comments.

import { sourceAt } from './device.js';
import { chargeOnly, conservedMoieties, feeds, noCurrent } from './laws.js';
import { captureLaws, captureLocal } from './lawRows.js';
import { maxPotentialStep } from './newton.js';
import { violentStart } from './timeSteps.js';

export class SteadySolves {
  // Slowest diffusion time across the device, used to size "giant" steps.
  // (Of the materials some region is made of: one defined but unused, its D = 1e-22, made a
  // giant step 1e7 times longer, too long to keep a fast species' amount: see _steadySteps.)
  slowestTime() {
    const { model } = this;
    let dMin = Infinity;
    for (const mat of this._usedMaterials()) {
      for (let i = 0; i < this.n; i++) if (mat.present[i] && mat.D[i] > 0) dMin = Math.min(dMin, mat.D[i]);
    }
    const L = model.grid.length;
    return Number.isFinite(dMin) ? (L * L) / dMin : 1;
  }

  _usedMaterials() {
    return [...new Set(this.model.regions.map((reg) => reg.material))].map((m) => this.model.materials[m]);
  }

  // The fastest diffusion time: across the smallest cell, at the largest D.
  fastestTime() {
    let dMax = 0;
    for (const mat of this._usedMaterials()) {
      for (let i = 0; i < this.n; i++) if (mat.present[i]) dMax = Math.max(dMax, mat.D[i]);
    }
    let h = Infinity;
    for (const s of this.model.grid.segLength) if (s > 0) h = Math.min(h, s);
    return dMax > 0 && Number.isFinite(h) ? (h * h) / dMax : this.slowestTime();
  }

  /**
   * Steady state: where a transient from the present state ends, found as the solution of the
   * steady equations (dt = ∞, no storage term), with no slow transport to wait out (e.g.
   * exponentially scarce minority carriers); but a slow reaction that alone settles a closed
   * combination is waited out by time steps, to ~1e12 of the slowest diffusion time (invariant
   * 5). What the storage term held is held by
   * the steady system's own rows instead (see src/laws.js): each conserved amount in place of a balance
   * row that says nothing new, levels that carry nothing held flat, islands' summed balances.
   * Where Newton from here fails: continuation along the drives, then the pseudo-transient ramp.
   * The clock is not advanced, and open-system conservation bookkeeping restarts here.
   */
  // Steady state from the present one; `atSteady` says whether the state is one (until a step or
  // a change of drive).
  solveSteady(opts = {}) {
    this.violent = this.history.length === 0 && !this.stepped ? violentStart(this) : []; // (for a failure to say: see src/timeSteps.js)
    this.unreached = null; // (a driven current no held voltage reaches: see _floatingContinuation)
    this.charging = null; // (a driven current that only charges: see _solveSteadyAll)
    captureLaws(this);
    this.steady = true; // (the steady system's rows; sources at the present time)
    this.steadyAt = this.time; // (which a pseudo-transient's steps don't move on)
    this.drifted = null; // (an amount no row keeps, lost by a giant step: see _steadySteps)
    // A solve that fails leaves where the state was solved, and whether it was stepped, as they
    // were (Device.solve restores the state): a floating continuation's held solves had set them,
    // and the next held solve ramped warm from a voltage the state was never at, and failed.
    const was = [this.solvedV, this.stepped];
    try {
      const r = this._solveSteadyAll(opts);
      if (!r.converged) [this.solvedV, this.stepped] = was;
      return r;
    } finally {
      this.steady = false;
    }
  }

  _solveSteadyAll(opts) {
    // A terminal driven by a current that can only store what comes in (a battery's cathode
    // behind it, closed but for that terminal; a gate, or a capacitance with no species through
    // it) has no steady state but at no current, where it keeps what it holds.
    const { species, regions } = this.model;
    for (const k of this.floating) {
      const I = sourceAt(this.terms[k].drive.src, this.time);
      if (this.terms[k].drive.kind !== 'I' || I === 0) continue;
      const laws = this.constraints.filter((cs) => cs.terminals.includes(k));
      if (laws.length === 0 && !chargeOnly(this, k)) continue;
      const what = laws.flatMap((cs) => cs.parts.map(({ stretch }) => {
        const st = this.stretches[stretch], [r0, r1] = st.regions;
        return `${species[st.species].name} in ${regions[r0].name}${r1 > r0 ? `…${regions[r1].name}` : ''}`;
      }));
      this.charging = { terminal: this.terms[k].name, current: I, what };
      return { converged: false, steps: 0, iterations: 0, history: [] };
    }
    const r = this._steadyFromHere(opts);
    this.atSteady = r.converged;
    return r;
  }

  // Steady state from the present one, in three stages: Newton from here; failing that, from a
  // solved state along the drives (the light ramped up, a current's terminal held at the voltages
  // that pass it, the voltages ramped from where they were last solved); failing that, the
  // pseudo-transient ramp, through time.
  _steadyFromHere(opts) {
    this.conditioning = null;
    this.computeConcentrations();
    captureLocal(this); // immobile combinations keep, node by node, what they hold now
    this.sourceTime = this.time;
    // Continuation applies where the right terminal is held at a voltage, against a held left one.
    const [dl, dr] = [this.terms[0].drive, this.terms[1].drive];
    const held = (d) => d.kind === 'V' && !(d.R > 0);
    const target = sourceAt(dr.src, this.time), level = sourceAt(dl.src, this.time);
    const allowed = opts.continuation !== false;
    const canContinue = allowed && held(dl) && held(dr) && target !== level;
    const direct = this._directSteady();
    const u0 = Float64Array.from(this.u), u0Lo = Float64Array.from(this.uLo), v0 = Float64Array.from(this.termV);
    const restart = () => {
      this.u.set(u0);
      this.uLo.set(u0Lo);
      this.termV.set(v0);
      this.computeConcentrations();
    };
    const solved = (q) => {
      if (q.converged) {
        this.solvedV = [this.termV[0], this.termV[1]]; // (a contact's own: driven by a current, not the current)
        this.stepped = false;
      }
      return q;
    };
    const after = (q, r) => ({ ...q, steps: q.steps + r.steps, iterations: q.iterations + r.iterations });

    // 1. Newton from here. Where a continuation could follow, one direct attempt only, rather
    // than the pseudo-transient ramp straight away.
    const quick = (canContinue || this.hasGeneration) && direct && allowed;
    let r = this._steadySteps(quick ? { ...opts, maxSteps: 1 } : opts);
    if (r.converged) return solved(r);
    // A transient's state can be a worse start than a cold one: 1 ns after a lit MOS capacitor's
    // gate jumps, its surface is in deep depletion, minority carriers ten orders below where they
    // end, and Newton fails from there where it converges from the layout for the drives (and
    // the stages below took 9000 iterations from there, 2000 from the layout). What a steady
    // solve keeps it took from the state already, so from here on it goes as a cold device would.
    if (this.stepped && direct && allowed) {
      this._refreshSources();
      this.initFromComposition();
      u0.set(this.u);
      u0Lo.set(this.uLo);
      const q = this._steadySteps(quick ? { ...opts, maxSteps: 1 } : opts);
      if (q.converged) return solved(after(q, r));
      restart();
      r = after(q, r);
    }

    // 2. Continuation. Generation (e.g. light) holding the device far from equilibrium: ramp it up
    // from nearly nothing, each solve warm from the last.
    let dimmer = null;
    if (this.hasGeneration && allowed) {
      restart();
      this.dimmer = false;
      const g = this._generationContinuation(opts, r);
      if (this.dimmer) dimmer = [Float64Array.from(this.u), Float64Array.from(this.uLo), Float64Array.from(this.termV)];
      if (g.converged) return solved(g);
      r = g;
    }
    // A terminal driven by a current (open circuit, say): hold it at a voltage instead, march the
    // voltage until the current crosses its target, and float it from there.
    if (allowed && this.floating.some((k) => this.terms[k].drive.kind === 'I' && !this._bound(k))) {
      restart();
      const q = this._floatingContinuation(opts, r);
      if (q.converged) return q;
      r = q;
    }
    if (canContinue) {
      // The right terminal's voltage ramped to its target: from the contacts' voltages at the last
      // converged solve when the state is that solution (a warm start), the left one first where
      // it moved, the right held where it was; otherwise from level terminals, consistent with a
      // cold start.
      restart();
      const was = this.solvedV;
      let c;
      if (was) {
        c = was[0] === level ? { ...r, converged: true, ramped: false } : this._continuation(opts, 0, was[0], level, r, true, was[1]);
        // (Ramping the right from level terminals, solve there first: a sweep's start.)
        if (c.converged) c = this._continuation(opts, 1, was[1], target, c, c.ramped || was[1] !== level);
      } else c = this._continuation(opts, 1, level, target, r, false);
      if (c.converged) return solved(c);
      // 3. The pseudo-transient ramp, from the start.
      restart();
      return solved(after(this._steadySteps(opts), c));
    }
    // 3. The pseudo-transient ramp (where the first attempt wasn't it already): from the dimmer
    // light's solution where the light's continuation got anywhere (a lit floating base charges
    // from there as it would in time), else from the start.
    if (!quick) return r;
    if (!dimmer) restart();
    else {
      // (where it was left, not where a floating continuation since left the state)
      this.u.set(dimmer[0]);
      this.uLo.set(dimmer[1]);
      this.termV.set(dimmer[2]);
      this.computeConcentrations();
    }
    return after(this._steadySteps(opts), r);
  }

  // Whether terminal k, driven by a current, feeds a combination that it alone can change, or
  // only charges a capacitance (so that its steady states are one per amount or charge, not one
  // per voltage).
  _bound(k) {
    return this.constraints.some((cs) => cs.terminals.includes(k)) || chargeOnly(this, k);
  }

  // A cold steady solve with a terminal driven by a current, which can fail where the same state
  // is easy warm (a solar cell's open circuit, from cold on a fine grid): the terminal held at a
  // voltage instead, whose steady solves have their own continuation, the voltage marched from
  // the start's until the current crosses its target, the crossing narrowed by bisection, and
  // the terminal floated again from beside it.
  _floatingContinuation(opts, r) {
    const k = this.floating.find((j) => this.terms[j].drive.kind === 'I' && !this._bound(j));
    const drive = this.terms[k].drive, target = sourceAt(drive.src, this.time);
    // Its drive put back: redrive() finds the laws again for it, with the local constraints'
    // references at zero, so they're taken from the state again (an immobile species held at 0
    // at every node: the floated solve failed, 80 steps, and the regula falsi rescued it at 6×).
    const floatAgain = () => {
      this.terms[k].drive = drive;
      this.redrive();
      captureLocal(this);
    };
    let steps = r.steps ?? 0, iterations = r.iterations ?? 0;
    const history = (r.history ?? []).slice();
    const hold = (V) => {
      this.terms[k].drive = { kind: 'V', src: { value: V }, R: 0 };
      this.redrive();
      const q = this._steadyFromHere(opts);
      steps += q.steps ?? 0;
      iterations += q.iterations ?? 0;
      history.push({ held: V, converged: q.converged });
      if (q.converged) seen.push([V, this.termI[k]]);
      return q.converged ? this.termI[k] - target : NaN;
    };
    const fail = () => ({ converged: false, steps, iterations, history });
    const seen = []; // [V, I] of each held solve that converged
    let Va = this.termV[k], fa, Vb, fb;
    try {
      fa = hold(Va);
      if (!Number.isFinite(fa)) {
        Va = 0;
        fa = hold(Va);
      }
      if (!Number.isFinite(fa)) return fail();
      // March, the step growing: up where the current is below its target, since raising a
      // terminal's voltage raises the current into a passive device (the other way, from the
      // start again, if 20 V brings no crossing).
      const start = Va, f0 = fa;
      [Vb, fb] = [Va, fa];
      for (const dir of [fa < 0 ? 1 : -1, fa < 0 ? -1 : 1]) {
        [Va, fa] = [start, f0];
        if (dir !== (f0 < 0 ? 1 : -1)) hold(start); // (back to the start, for a warm march)
        let dV = 0.02;
        Vb = Va + dir * dV;
        fb = hold(Vb);
        while (Number.isFinite(fb) && Math.sign(fb) === Math.sign(fa) && fa !== 0 && Math.abs(Vb - start) < 20) {
          [Va, fa] = [Vb, fb];
          dV *= 1.5;
          Vb = Va + dir * dV;
          fb = hold(Vb);
        }
        if (Number.isFinite(fb) && (Math.sign(fb) !== Math.sign(fa) || fa === 0)) break;
      }
      if (!Number.isFinite(fb) || (Math.sign(fb) === Math.sign(fa) && fa !== 0)) {
        // What the held voltages passed, at the two extremes reached.
        if (seen.length > 1) {
          seen.sort((a, b) => a[0] - b[0]);
          this.unreached = { terminal: this.terms[k].name, target, low: seen[0], high: seen.at(-1) };
        }
        return fail();
      }
      // Bisect (each held solve warm from the last) to a tenth of a millivolt.
      for (let it = 0; it < 40 && Math.abs(Vb - Va) > 1e-4; it++) {
        const Vm = (Va + Vb) / 2, fm = hold(Vm);
        if (!Number.isFinite(fm)) return fail();
        if (Math.sign(fm) === Math.sign(fa)) [Va, fa] = [Vm, fm];
        else [Vb, fb] = [Vm, fm];
      }
      hold(Math.abs(fa) < Math.abs(fb) ? Va : Vb);
    } finally {
      floatAgain();
    }
    const q = this._steadySteps(opts);
    steps += q.steps ?? 0;
    iterations += q.iterations ?? 0;
    history.push(...(q.history ?? []));
    if (q.converged) return { ...q, steps, iterations, history };
    // Floated, the system can lose what held it doesn't (a closed redox cell driven near its
    // limit: its bordered row 33 digits short), but a current-driven steady state is the held
    // one at the voltage that passes the target: found by held solves, by regula falsi
    // (Illinois) within the bracket, to 1e-10 of the target, or at open circuit of the largest
    // current the march saw (a target of 0 is never met exactly).
    const scale = Math.max(Math.abs(target), ...seen.map(([, I]) => Math.abs(I)));
    try {
      let side = 0;
      for (let it = 0; it < 60; it++) {
        const Vm = Va - (fa * (Vb - Va)) / (fb - fa), fm = hold(Vm);
        if (!Number.isFinite(fm)) return fail();
        if (Math.abs(fm) <= 1e-10 * scale) return { converged: true, steps, iterations, history, residual: Math.abs(fm) };
        if (Math.sign(fm) === Math.sign(fb)) {
          [Vb, fb] = [Vm, fm];
          if (side === -1) fa /= 2;
          side = -1;
        } else {
          [Va, fa] = [Vm, fm];
          if (side === 1) fb /= 2;
          side = 1;
        }
      }
      return fail();
    } finally {
      floatAgain();
    }
  }

  // Steady solves with the generation reactions' rates scaled from 1e-30 up to 1, ×100 a step
  // while they converge, smaller steps where they don't.
  _generationContinuation(opts, r) {
    const sub = this._directSteady() ? { ...opts, maxSteps: 1 } : opts;
    let s = 1e-30, factor = 100, steps = r.steps, iterations = r.iterations;
    const history = r.history.slice();
    try {
      this.generationScale = s;
      // (nearly dark: a cold start can need the full ramp, as a dark device's does)
      let q = this._steadySteps(sub);
      steps += q.steps;
      iterations += q.iterations;
      if (!q.converged) return { converged: false, steps, iterations, history };
      this.dimmer = true; // (the state is now a solution under dimmer light: a better start than cold)
      while (s < 1) {
        const next = Math.min(1, s * factor);
        const u1 = Float64Array.from(this.u), u1Lo = Float64Array.from(this.uLo), v1 = Float64Array.from(this.termV);
        this.generationScale = next;
        q = this._steadySteps(sub);
        steps += q.steps;
        iterations += q.iterations;
        history.push({ generation: next, converged: q.converged });
        if (q.converged) {
          s = next;
          factor = Math.min(100, factor * factor);
        } else {
          this.u.set(u1);
          this.uLo.set(u1Lo);
          this.termV.set(v1);
          this.computeConcentrations();
          factor = Math.sqrt(factor);
          if (factor < 1.5) return { converged: false, steps, iterations, history }; // (crawling: time steps do better)
        }
      }
      return { converged: true, steps, iterations, history };
    } finally {
      this.generationScale = 1;
    }
  }

  // Terminal k (a contact) ramped from `level` to `target` through steady solves, each warm from
  // the last (`hold`: the other contact's voltage meanwhile, if not its own). A warm ramp that has
  // nowhere to go converges only if what it started from did (`ramped`): the state it carries
  // must be some ramp's solution, never just the one a failed solve left.
  _continuation(opts, k, level, target, r, warm, hold) {
    const sub = this._directSteady() ? { ...opts, maxSteps: 1 } : opts;
    let V = level, dV = (target - level) / 8, steps = r.steps, iterations = r.iterations;
    const history = r.history.slice();
    if (warm && V === target) return { ...r, converged: r.ramped === true, steps, iterations, history, ramped: r.ramped };
    try {
      this.sourceOverride.set(k, V);
      if (hold !== undefined) this.sourceOverride.set(1 - k, hold);
      let q;
      if (!warm) {
        q = this._steadySteps(opts);
        if (!q.converged) {
          // A cold device's start was laid out for the target (a set() before the first solve
          // builds it there): lay it out again, level.
          steps += q.steps;
          iterations += q.iterations;
          this._refreshSources();
          this.initFromComposition();
          q = this._steadySteps(opts);
        }
        // (lit, the level start may need the light ramped up too)
        if (!q.converged && this.hasGeneration) q = this._generationContinuation(opts, q);
        steps += q.steps;
        iterations += q.iterations;
        if (!q.converged) return { ...r, steps, iterations };
      }
      while (V !== target) {
        const next = Math.abs(target - V) <= Math.abs(dV) ? target : V + dV;
        const u1 = Float64Array.from(this.u), u1Lo = Float64Array.from(this.uLo), v1 = Float64Array.from(this.termV);
        this.sourceOverride.set(k, next);
        q = this._steadySteps(sub);
        steps += q.steps;
        iterations += q.iterations;
        history.push({ continuation: next, converged: q.converged });
        if (q.converged) {
          V = next;
          dV *= 1.5;
        } else {
          this.u.set(u1);
          this.uLo.set(u1Lo);
          this.termV.set(v1);
          this.computeConcentrations();
          dV /= 4;
          if (Math.abs(dV) < 1e-6 * Math.abs(target - level)) return { converged: false, steps, iterations, history };
        }
      }
      return { converged: true, steps, iterations, history, ramped: true };
    } finally {
      this.sourceOverride.delete(k);
      if (hold !== undefined) this.sourceOverride.delete(1 - k);
      this._refreshSources();
    }
  }

  _steadySteps({ maxSteps = 80, tol = 1e-11 } = {}) {
    const time = this.time;
    const tau = this.slowestTime();
    const direct = this._directSteady();
    const giant = direct ? Infinity : 1e6 * tau;
    let dt = giant, grow = 10; // (after a failure, dt grows by less, then back up to ×10)
    let triedDirect = 0;
    let totalIter = 0, steps = 0, converged = false;
    const history = [];
    // What no row keeps (laws.js's unkept) a giant step keeps only by its storage term, which can
    // fall below round-off against the fast species' conduction, and Newton then converges
    // anywhere along it: a closed Y + Z, Y still in part of its stretch, converged 1e8× too large.
    // So each such amount is checked after every step of dt ≥ giant, and a step that moved it
    // fails.
    const unkept = !direct && this.unkept.length > 0 ? this._unkeptAmounts() : null;
    const drifted = () => {
      const now = this._unkeptAmounts(), lost = now.findIndex(([a], j) => !(Math.abs(a - unkept[j][0]) <= 1e-6 * unkept[j][1]));
      if (lost >= 0) this.drifted = now[lost][2];
      return lost >= 0;
    };
    while (steps < maxSteps) {
      steps++;
      // Restore conserved amounts exactly before each huge step; the step then re-solves, so
      // the final state satisfies every equation.
      if (!direct && dt >= giant) this._renormalizeSpectators();
      const snap = unkept && dt >= giant ? this._snapshot() : null;
      let r = this.step(dt);
      if (snap && r.converged && drifted()) {
        this._restore(snap);
        r = { ...r, converged: false, drifted: true };
      }
      totalIter += r.iterations;
      history.push({ dt, converged: r.converged, iterations: r.iterations, ...(r.drifted ? { drifted: true } : {}) });
      if (!r.converged && dt === Infinity && steps === 1) {
        // A large change can overshoot Newton's usual damping: retry with a tighter limit.
        r = this.step(dt, { maxStep: 3, maxIter: 80 });
        totalIter += r.iterations;
        history.push({ dt, converged: r.converged, iterations: r.iterations, maxStep: 3 });
      }
      if (!r.converged) {
        // Pseudo-transient continuation: down from a slow time scale until a step converges,
        // as far as the fastest (slow ions mustn't stop it short of what a cold start needs).
        dt = dt >= giant ? tau * 1e-6 : dt / 4;
        grow = 2;
        if (dt < Math.min(tau * 1e-15, this.fastestTime() * 1e-2)) break;
        continue;
      }
      if (dt === Infinity) {
        converged = true; // the steady equations themselves were solved
        break;
      }
      // Past the slowest diffusion time, the direct solve may already reach from here: try it
      // once a decade (the state kept if it fails), rather than step on to the giant dt through
      // steps that can fail for slow ions' sake.
      if (direct && dt < giant && dt >= tau && dt >= 10 * triedDirect) {
        triedDirect = dt;
        const snap = this._snapshot();
        const d = this.step(Infinity);
        totalIter += d.iterations;
        history.push({ dt: Infinity, converged: d.converged, iterations: d.iterations });
        if (d.converged) {
          converged = true;
          break;
        }
        this._restore(snap);
      }
      if (dt >= giant && maxPotentialStep(this, this._diff()) < tol) {
        converged = true;
        if (dt > giant) {
          // Finish at the base giant step, where the storage term pins amounts tightly.
          this._renormalizeSpectators();
          const f = this.step(giant);
          totalIter += f.iterations;
          converged = f.converged && !(unkept && drifted());
        }
        break;
      }
      if (dt < giant) {
        dt *= grow; // ramping up after a failure
        grow = Math.min(10, grow * 1.5);
        if (dt > 1e6 * tau) dt = giant; // then the direct steady solve (or the giant step)
      } else if (dt < 1e6 * giant) {
        dt *= 10; // conserved amounts present and still moving: let dt keep growing (capped)
      }
    }
    this.time = time;
    this.history = []; // a steady solve isn't a trajectory: transients restart from here
    this.dtNext = undefined;
    if (converged) {
      this.stretches.forEach((st, k) => {
        if (!st.connected) return;
        this.referenceAmounts[k] = this.amount(st);
        this.boundaryIntake[k] = 0;
      });
    }
    return { converged, steps, iterations: totalIter, history };
  }

  // The amounts no row keeps (unkept) that are of stretches and surfaces alone: each as [its value,
  // the sum of its parts' sizes, the names of its parts].
  _unkeptAmounts() {
    const S = this.stretches.length, X = S + this.surfCols.length, { species, ports } = this.model;
    return this.unkept.filter((w) => w.every((x, k) => k < X || x === 0)).map((w) => {
      let a = 0, size = 0;
      const names = [];
      w.forEach((x, k) => {
        if (x === 0) return;
        const v = x * (k < S ? this.amount(this.stretches[k]) : this.surfaceAmount(k - S));
        a += v;
        size += Math.abs(v);
        names.push(k < S ? species[this.stretches[k].species].name : ports[this.surfCols[k - S][0]].surface[this.surfCols[k - S][1]].name);
      });
      return [a, size, names];
    });
  }

  // Huge steps pin each conserved amount only through a tiny storage term, so round-off can let
  // it creep. Shift each spectator's level uniformly to restore its amount exactly (in one step
  // for ideal statistics, where c ∝ e^η; by Newton on the shift otherwise).
  _renormalizeSpectators() {
    const { n, M, u } = this;
    let changed = false;
    this.stretches.forEach((st, k) => {
      // A floating conductor holds its charge on its faces, not as an amount: the huge steps keep
      // it, and there's nothing here to shift.
      if (!st.spectator || this.model.materials[this.model.regions[st.regions[0]].material].conductor) return;
      const want = this.referenceAmounts[k];
      const ideal = this.nodeIdeal.subarray(st.nodes[0], st.nodes[1] + 1).every((v) => v === 1);
      for (let it = 0; it < (ideal ? 1 : 30); it++) {
        const now = this.amount(st);
        if (!(now > 0) || now === want) return;
        let shift;
        if (ideal) shift = Math.log(want / now);
        else {
          let dA = 0;
          for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
            const dc = this.nodeIdeal[g] ? this.c[g * n + st.species] : this.K[g * n * n + st.species * (n + 1)];
            dA += this.model.grid.vol[g] * dc;
          }
          shift = (want - now) / dA;
          if (!(Math.abs(shift) > 1e-16)) return;
        }
        for (let g = st.nodes[0]; g <= st.nodes[1]; g++) u[this.blockOfNode[g] * M + 1 + st.species] += shift;
        changed = true;
        if (!ideal) this.computeConcentrations();
      }
    });
    if (changed) this.computeConcentrations();
  }

  _diff() {
    const d = this.delta; // reuse as scratch
    const { fullOf } = this;
    for (let j = 0; j < fullOf.length; j++) d[j] = this.u[fullOf[j]] - this.uPrev[fullOf[j]];
    return d;
  }

  // Every stretch either reaches a contact or is kept by a law, and no law is left unkept: the
  // steady equations can be solved directly (dt = ∞), with the laws in place of their rows. Else
  // the steady solve goes by time steps, which keep every amount themselves.
  _directSteady() {
    return this.unkept.length === 0 && this.stretches.every((st) => st.connected || st.conserved);
  }

  // After a steady solve: each reaction that alone changes some combination nothing else feeds
  // (see src/laws.js), whose total rate a steady state therefore sets to zero, but whose affinity has one
  // sign wherever it acts, so that its total can't be zero. Too slow to show, its rate lost to
  // round-off, it was left where the solve began: a slow complexation 43 mV off, an electrode
  // port's slow couple at its starting composition, each converged without a word. (A face
  // reaction stands on its affinity, so it's only caught here if that couldn't be arranged.)
  stuckReactions() {
    const { n, M, u, uLo, model } = this, grid = model.grid, S = this.stretches.length, X = S + this.surfCols.length, out = [];
    const passes = (k) => !noCurrent(this.terms[k]);
    const fed = (w) => w.some((x, k) => x !== 0 && ((k >= X && passes(this.carriers[k - X])) || (k < S && this.terms.some((t, j) => t.drive.kind === 'I' && passes(j) && feeds(j, this.stretches[k])))));
    const needed = (rx) => {
      const basis = conservedMoieties(this._structure(), { drives: true, skip: rx });
      return basis.some((w) => !fed(w) && basis.skipped.some((row) => Math.abs(row.reduce((t, x, k) => t + x * w[k], 0)) > 1e-9));
    };
    const eta = (g, slot) => u[this.blockOfNode[g] * M + slot] + uLo[this.blockOfNode[g] * M + slot];
    const report = (where, as) => {
      if (as.length === 0 || !(as.every((a) => a > 1e-6) || as.every((a) => a < -1e-6))) return;
      out.push({ where, affinity: as.reduce((m, a) => (Math.abs(a) < Math.abs(m) ? a : m)) });
    };
    model.interfaces.forEach((itf, f) => {
      const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
      itf.reactions.forEach((rx, q) => {
        if (this.affinities.some((af) => af.rx === rx) || !needed(rx)) return;
        report(`interfaces[${f}].reactions[${q}]`, [rx.part.reduce((a, p) => a - p.nu * eta(p.side ? gR : gL, 1 + p.i), rx.fixedA)]);
      });
    });
    model.reactions.forEach((mrx, x) => {
      const rx = this.rxs[x];
      if (rx.generation || !needed(mrx)) return;
      const as = [];
      for (let g = 0; g < this.nNodes; g++) if (rx.kfNode[g] > 0 || rx.srh?.[this.nodeMaterial[g]]) as.push(rx.sp.reduce((a, i, p) => a + rx.nu[p] * eta(g, 1 + i), rx.fixedA));
      report(`bulkReactions[${x}]`, as);
    });
    model.ports.forEach((port, k) => {
      port.reactions.forEach((rx, q) => {
        if (!needed(rx)) return;
        let aV = rx.fixedA;
        for (const p of rx.part) if (p.side === 1) aV -= p.nu * this.portEta(port, p.i);
        report(`ports[${k}].reactions[${q}]`, port.nodes.map((g) => rx.part.reduce((a, p) => (p.side === 1 ? a : a - p.nu * eta(g, p.side === 2 ? 1 + n + p.s : 1 + p.i)), aV)));
      });
    });
    return out;
  }
}
