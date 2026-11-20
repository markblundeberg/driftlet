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
//   here and the pseudo-transient ramp only). Starts the call's diagnostics fresh
//   (this._freshDiagnostics(): this.violent, violentStart where nothing has been stepped or
//   solved; this.unreached, charging, drifted, settling, modelError and conditioning cleared),
//   finds the laws again where how a terminal is driven, read at this.time, has changed since
//   they were found (this._lawsNow: a waveform current zero now and not then, or the reverse),
//   captures the laws and the local constraints' sums (keepLocal: each stage holds those the
//   solve began with, a cold layout's or a redrive()'s notwithstanding), and
//   solves with this.steady = true and this.steadyAt = this.time (false again on return).
//   Stops at once, converged false, steps 0, where a floating terminal is driven by a nonzero
//   current that only charges a capacitance or feeds a conserved combination (a constraint):
//   this.charging says which terminal and what it fills. Else _steadyFromHere, and
//   this.atSteady = converged. Not converged: this.solvedV and this.stepped as before the call.
//
// _steadyFromHere(opts): in stages until one converges, each from the state the solve began in
// (the whole of STATE, put back between stages by _restore, and this.carried with it) unless it
// says otherwise. A continuation applies with continuation not off; canContinue: both contacts held
// at voltages (not behind a resistance) that differ.
//   1. Newton from here (_steadySteps): a single step (maxSteps 1) where canContinue or there's
//      generation, and the device can be solved directly (_directSteady); else the
//      pseudo-transient ramp with opts' maxSteps. Failing on a state reached by time steps
//      (this.stepped), on a direct device with continuation not off, the same again from a cold
//      start (initFromComposition), whose u and uLo become the stages' start (termV stays the
//      start's: a floating terminal keeps the voltage the transient left it at).
//   2. Continuation: generation ramped up (_generationContinuation); terminals behind a
//      resistance (_resistanceContinuation): one whose current Kirchhoff doesn't give, held at
//      voltages marched to its load line's crossing (_floatingContinuation with a line); else, or
//      where that fails, each held (at its source voltage, or where every other terminal is
//      driven by a current at V_src − I R, Kirchhoff's I; from a cold start not stepped, the start
//      laid out again for the held drive, initFromComposition), then the resistances put back and
//      the steady system solved from there (and where that fails for one terminal held by
//      Kirchhoff, its load line); a terminal driven by a current, not bound (_bound), held at
//      voltages instead (_floatingContinuation); a port held at a voltage across a capacitance,
//      ramped from its zeroCharge (_continuation, cold there, with the contacts held level where
//      canContinue, then the right one ramped warm to its target); with
//      canContinue, the right contact's voltage ramped to its target (_continuation): from the
//      voltages of the last converged solve (this.solvedV), the left contact's first where it
//      moved, else from level contacts, a cold start.
//   3. The pseudo-transient ramp: after the contact's ramp fails, from the start; otherwise only
//      where stage 1 was a single step, from the dimmer light's solution where the generation
//      ramp got past its first solve (the state as that ramp left it), else from the start.
//   Converged (but by _floatingContinuation): this.solvedV = [termV[0], termV[1]], the contacts'
//   own voltages at the solution (driven by a current or behind a resistance too), this.stepped =
//   false. (By _floatingContinuation: as its last held solve set them.)
//
// _steadySteps({ maxSteps, tol }): up to maxSteps time steps (tau = slowestTime), dt from giant
//   (∞ where the device can be solved directly; else 1e6·tau): a first step at ∞ that fails is
//   retried once with maxStep 3. A step that fails sends dt down (to tau·1e-6 from giant, then
//   ÷4, stopping below min(1e-15·tau, 1e-2·fastestTime)); one that converges sends it back up
//   (×2 after a failure, growing by ×1.5 to ×10 a step; past 1e6·tau, giant, or on a direct
//   device past 1e12·tau); on a direct device,
//   once dt ≥ tau, a step at ∞ is tried once a decade of dt, the state kept where it fails. Past
//   giant (not direct), dt grows ×10 a step to 1e6·giant (Number.MAX_VALUE at most). Converged where a step at ∞ converged,
//   or, not direct, where a step of dt ≥ giant moved the potentials and the floating terminals' voltages (termV − termVPrev) under tol (finished by one
//   step of exactly giant where dt had grown past it, which must converge too). Not direct,
//   before each step of dt ≥ giant and the finishing one, each spectator's amount restored
//   (_renormalizeSpectators: in turn, again while a shift by Newton, a coupled statistics', left
//   one off by over 1e-12, fifty passes at most). Not direct, with unkept laws: after each step of dt ≥ giant, each
//   unkept combination of stretches' and surfaces' amounts (_unkeptAmounts) must be within 1e-6
//   of its size (Σ |w·amount|) of what it was when solveSteady began (this.solveAmounts; when
//   the call began, called alone), else the step is put back
//   (_restore) and counts as failed (`drifted` in its history entry; this.drifted names the
//   combination's parts); the finishing step, failing that, ends the call unconverged. (A
//   combination through an electrode's or a contact's charge isn't checked.) On return, whatever
//   the outcome: this.time as it was, the history
//   empty, this.dtNext undefined; converged, every connected stretch's referenceAmounts its
//   amount now and its boundaryIntake 0. Out of steps with the last three giant steps' motions
//   each smaller than the one before, by ratios under 0.95 within 0.2 of each other: this.settling = { maxSteps, ratio, motion, tol }, which the
//   Solution's failure says (with the maxSteps that would finish).
//
// The continuations (each warm from the last converged solve; a solve that fails puts back the
// state, as between stages, and takes a smaller step):
//   _generationContinuation  generationScale from 1e-30 up to 1, by a factor of 100 at most
//                            (squared after a success, its square root after a failure, giving
//                            up below 1.5); fails if the first, at 1e-30, does; this.dimmer once
//                            it converges. generationScale is 1 again on return.
//   _continuation(k, from, to)  contact k's voltage through sourceOverride (the other's held
//                            too, where given), in steps of an eighth of the span, ×1.5 on
//                            success, ÷4 on a failure, giving up below 1e-6 of the span (or
//                            at a stride of 0, where that millionth underflows). Cold
//                            (not warm), first solved at `from`: failing, with no law unkept,
//                            again from a cold start laid out there, then with the light ramped where there's
//                            generation. The overrides are cleared on return.
//   _floatingContinuation    the terminal held at voltages (its drive replaced, then put back):
//                            marched from its voltage (or 0) in steps from 20 mV, ×1.5, toward
//                            the target current (not at all where the start passes it exactly;
//                            both ways, to ±20 V; a held solve that fails
//                            retried at a quarter of the stride from the last converged one, to
//                            1 mV), the crossing bisected to 0.1 mV, then floated from the
//                            nearer end. Floated and failing, the
//                            held solve whose current is within 1e-10 of the target (by regula
//                            falsi, Illinois, from the bracket; the scale is the larger of the
//                            target and the bracket's ends' currents; an end that passes it
//                            exactly, that end) is returned as converged, its drive put back. No
//                            crossing: this.unreached = { terminal, target, V: [lo, hi], I: [lo,
//                            hi], nearing }, the ranges of the voltages held and the currents
//                            they passed, where each way showed a limit (reached ±20 V no
//                            nearer the target, to a part in 1e6, than before, or stopped by
//                            failing held solves after the current had turned away from it) or
//                            reached ±20 V still nearing it at a pace that, kept up, passes it
//                            within ten times as far again (nearing: { at: that voltage, guess:
//                            where that pace passes it }, else null; one way at most). Each time its
//                            drive is put back, the local constraints' references are put back
//                            as the solve began (captureLocal, from this.localKept), which
//                            redrive() zeroed.
//
// Also: slowestTime() = L²/D_min over the species present with D > 0 in each region (a graded
// D's values over its region in place of its material's; 1 where none), fastestTime() = h_min²/D_max over the same (slowestTime where
// there's no such cell); _directSteady(): no law unkept and every
// stretch connected or conserved; stuckReactions(): after a steady solve, each reaction that
// alone changes a combination nothing feeds and whose affinity keeps one sign, over 1e-6, wherever
// it acts, { where, affinity } (the smallest); each neutral species fed only by contacts' rate
// links (exchange, velocity), unreacted, with no port or flow, whose links all push it one way
// by over 1e-6 (η outside − η at the contact), { where: the links, affinity }, or whose links'
// summed rates over their slope put it over 1e-6 from where they balance, { where, affinity,
// links: how many }; each bulk reaction out of equilibrium, over 1e-6, at a node where it alone
// changes a participant that can't move, { where, affinity (the largest there), immobile: how
// many nodes } (see its comment); each bulk or face reaction out of equilibrium, over 1e-6, in a
// part no contact, port, flow or fixed participant (fixedA ≠ 0, or a generation) reaches,
// { where, affinity (the largest), closed: how many nodes }; each face's permeability or
// conductance in a stretch nothing makes or uses, with no port and at most one end open, whose
// levels differ by over 1e-6, { where, affinity, face: true }; and each terminal behind a
// resistance that passes current only by charging, over 1e-6 V_T from its source, { where: the
// terminal, affinity: (V − V_src)/V_T, charging: true }. With { all: false } (a failed solve's
// warning), the one-sign and immobile kinds only.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. Converged: the state solves the steady equations as step(∞) solved them; or, not direct,
//      a giant step's (whose storage term is 1e-6/tau of the rest at most); or, by
//      _floatingContinuation's regula falsi, the held equations at a voltage whose current is
//      within 1e-10 of the target's scale.
//   2. On return: this.time, every terminal's drive, sourceOverride, generationScale and
//      this.steady as before the call; the history empty, but where the call stopped at once
//      (this.charging) before any step.
//   3. A stage that fails leaves no trace in the next one's start: the whole of STATE is put back
//      (the capacitors' charges with the levels; the concentrations recomputed); only stage 1's
//      cold start, on a stepped state, moves the start (its u and uLo), and stage 3 may begin
//      from the generation ramp's solution.
//   4. Deterministic: the same state and options give the same result, bit for bit. The state
//      here is more than u, uLo, termV and time: this.solvedV, this.stepped, whether the history
//      is empty (this.violent) and the drives, all of which choose the stages.
//   5. Reach: time steps reach a mode as slow as ~1e12·tau (1e6·tau a step, past the giant; on a
//      direct device by finite steps to 1e12·tau). A reaction that alone settles a
//      closed combination, slower than that, isn't waited out: the solve fails or converges where
//      it began, and stuckReactions names it (Solution's warnings, both ways).
// Tests: every Device.solve goes through here; test/steady.test.js (the stages, the laws, the
// continuations), test/generation.test.js, test/round7-solve.test.js (behind R, cold) and
// test/round8-deep.test.js (a stage's rollback, the march's claims, absurd spans); the
// devices that needed each stage are named in its comments.

import { profileAt, sourceAt } from './device.js';
import { solveFlat } from './dense.js';
import { DeviceError } from './errors.js';
import { carrierZ, chargeOnly, conservedMoieties, feeds, noCurrent } from './laws.js';
import { captureLaws, captureLocal, keepLocal, planeSlot } from './lawRows.js';
import { maxPotentialStep } from './newton.js';

export class SteadySolves {
  // Slowest diffusion time across the device, used to size "giant" steps.
  // (Of the materials some region is made of: one defined but unused, its D = 1e-22, made a
  // giant step 1e7 times longer, too long to keep a fast species' amount: see _steadySteps.)
  // What a steady solve's stages save and put back: the whole of STATE, by the one rollback
  // (_snapshot/_restore, the concentrations recomputed), and with it this.carried, the call's own
  // say of whether the state is the solved one it began from. (As u, uLo and termV alone, a
  // stage that failed after finite steps left its capacitors' charges, contactDEnd and portQEnd,
  // beside the levels put back, and the next stage's steps stored against them: a floating
  // plate's Fe³⁺ + Q/F 0.58% off, 0.16 mV, unwarned; Astra's A16, round 8. Written out by hand
  // five times before that; a code-smell review.)
  _levels() {
    return { snap: this._snapshot(), carried: this.carried };
  }

  _putLevels({ snap, carried }) {
    this._restore(snap);
    this.carried = carried;
  }

  slowestTime() {
    const [dMin] = this._dRange();
    const L = this.model.grid.length;
    return Number.isFinite(dMin) ? (L * L) / dMin : 1;
  }

  // The smallest and largest D > 0 of a species present in some region: a graded D's values over
  // its region in place of its material's, which the grading replaces there (counted too, a
  // material's 1e-22 the grade replaced everywhere made the slowest time 1e12 s for 0.1 s;
  // Fable's re-check).
  _dRange() {
    const { model, n } = this, grid = model.grid;
    let dMin = Infinity, dMax = 0;
    const see = (d) => d > 0 && ((dMin = Math.min(dMin, d)), (dMax = Math.max(dMax, d)));
    model.regions.forEach((reg, r) => {
      const mat = model.materials[reg.material], a = grid.x[grid.regionStart[r]], b = grid.x[grid.regionEnd[r]];
      for (let i = 0; i < n; i++) {
        if (!mat.present[i]) continue;
        const p = reg.grade?.[i]?.D;
        if (!p) see(mat.D[i]);
        else [profileAt(p, a), profileAt(p, b), ...Array.from(p.values).filter((_, k) => p.x[k] > a && p.x[k] < b)].forEach(see);
      }
    });
    return [dMin, dMax];
  }

  // The fastest diffusion time: across the smallest cell, at the largest D.
  fastestTime() {
    const [, dMax] = this._dRange();
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
    this._freshDiagnostics(); // (violent, unreached, charging, drifted, settling, modelError, conditioning: see src/solver.js)
    this._lawsNow(); // (a waveform's current, zero or not at this time: see noCurrent in src/laws.js)
    captureLaws(this);
    this.computeConcentrations();
    keepLocal(this); // (what each node's immobile combinations hold now, through every stage)
    this.solveAmounts = this._columnAmounts(); // (what the giant steps keep of an unkept amount: see _steadySteps)
    this.steady = true; // (the steady system's rows; sources at the present time)
    this.steadyAt = this.time; // (which a pseudo-transient's steps don't move on)
    // A solve that fails leaves where the state was solved, and whether it was stepped, as they
    // were (Device.solve restores the state): a floating continuation's held solves had set them,
    // and the next held solve ramped warm from a voltage the state was never at, and failed.
    // (a solved state this call began from, which a nested stage's own solves don't count as)
    this.carried = this.solvedV != null;
    const was = [this.solvedV, this.stepped];
    try {
      const r = this._solveSteadyAll(opts);
      if (!r.converged) [this.solvedV, this.stepped] = was;
      return r;
    } finally {
      this.steady = false;
      this.solveAmounts = null;
      this.localKept = null;
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
      // (A law fed by several terminals charges only if their currents' net feed into it isn't
      // zero, Σ λ_t I_t, λ_t = w/z on what t feeds: +1 A/m² in at one terminal and −1 out at
      // another, each passing K⁺ into one pool, was refused as charging, where the transient
      // settles; Astra's fourth pass. An electrode's current enters through its charge column,
      // at λ = w/z there: two ports running Fe²⁺ = Fe³⁺ + e⁻ at ±1 mA/m², what one makes the
      // other takes, were refused at zero steps as charging the closed Fe³⁺; Astra's re-read of
      // the steady solve, round 8. A capacitance alone still charges: its column is all it feeds.)
      const net = (cs) => {
        let sum = 0, size = 0;
        for (const t of cs.terminals) {
          if (this.terms[t].drive.kind !== 'I') continue;
          const cr = cs.carriers?.find((q) => q.k === t);
          if (cr) {
            const x = (cr.w / carrierZ(this, t)) * sourceAt(this.terms[t].drive.src, this.time);
            sum += x;
            size += Math.abs(x);
            continue;
          }
          const p = cs.parts.find((q) => feeds(t, this.stretches[q.stretch]) && this.z[this.stretches[q.stretch].species] !== 0);
          if (!p) return true;
          const x = (p.w / this.z[this.stretches[p.stretch].species]) * sourceAt(this.terms[t].drive.src, this.time);
          sum += x;
          size += Math.abs(x);
        }
        return Math.abs(sum) > 1e-12 * size;
      };
      const laws = this.constraints.filter((cs) => cs.terminals.includes(k) && net(cs));
      if (laws.length === 0 && !chargeOnly(this, k)) continue;
      const what = laws.flatMap((cs) => [
        ...cs.parts.map(({ stretch }) => {
          const st = this.stretches[stretch], [r0, r1] = st.regions;
          return `${species[st.species].name} in ${regions[r0].name}${r1 > r0 ? `…${regions[r1].name}` : ''}`;
        }),
        // (an electrode surface's species: a tethered couple fed by the current, not a capacitance)
        ...(cs.surface ?? []).map(({ col }) => {
          const [k, slot] = this.surfCols[col], sf = this.model.surfaces[k];
          return `${sf.species[slot].name} on ${sf.face >= 0 ? `interfaces[${sf.face}]` : `ports ${sf.name}`}`;
        }),
      ]);
      // (and the other terminals driven by currents into the same pool, whose net feed it is)
      const others = [...new Set(laws.flatMap((cs) => cs.terminals))].filter((t) => t !== k && this.terms[t].drive.kind === 'I').map((t) => this.terms[t].name);
      this.charging = { terminal: this.terms[k].name, current: I, what, others };
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
    captureLocal(this); // immobile combinations keep, node by node, what they held as the solve began
    this.sourceTime = this.time;
    // Continuation applies where the right terminal is held at a voltage, against a held left one.
    const [dl, dr] = [this.terms[0].drive, this.terms[1].drive];
    const held = (d) => d.kind === 'V' && !(d.R > 0);
    const target = sourceAt(dr.src, this.time), level = sourceAt(dl.src, this.time);
    const allowed = opts.continuation !== false;
    // (Level contacts too, where the state was last solved at other voltages: back to level by
    // the same ramps. A Schottky diode at 0 V from −0.5 V failed, which a sweep reached; Fable
    // round 2's cold user.)
    const was = this.solvedV;
    const canContinue = allowed && held(dl) && held(dr) && (target !== level || (was != null && (was[0] !== level || was[1] !== target)));
    const direct = this._directSteady();
    const start0 = this._levels();
    const restart = () => {
      this._putLevels(start0);
    };
    const solved = (q) => {
      if (q.converged) {
        this.solvedV = [this.termV[0], this.termV[1]]; // (a contact's own: driven by a current, not the current)
        this.stepped = false;
        this.carried = false;
        this.coldStart = false; // (solved: no t = 0⁺ left for a transient to start from)
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
    // So can a solved state carried through set(): a Si diode at −5 V taken from 150 K to 300 K
    // failed warm (n_i moved 20 orders), and converges from the layout (Fable round 3's hunter).
    if ((this.stepped || this.carried) && direct && allowed) {
      this._refreshSources();
      this.initFromComposition();
      // (A transient's state is replaced by the layout for the stages below; a carried solved
      // state stays their start, which a continuation from its solved voltages wants: replaced,
      // GaAs taken warm through T at 0.3 V failed 17 of 56 jumps against 1.9.0's 4; Fable round
      // 3's hunter, on its re-check.)
      if (this.stepped) {
        start0.snap.u.set(this.u);
        start0.snap.uLo.set(this.uLo);
      }
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
      if (this.dimmer) dimmer = this._levels();
      if (g.converged) return solved(g);
      r = g;
    }
    // A terminal behind a resistance, cold: held at its source voltage first (with that its own
    // continuations), then the resistance put back and solved from there (a diode behind 7 Ω·m²
    // failed cold, and converged so; round 14's fuzzer, 4 of 9 such re-solves).
    if (allowed && this.floating.some((k) => this.terms[k].drive.kind === 'V' && this.terms[k].drive.R > 0)) {
      restart();
      const q = this._resistanceContinuation(opts, r);
      if (q.converged) return solved(q);
      r = q;
    }
    // A terminal driven by a current (open circuit, say): hold it at a voltage instead, march the
    // voltage until the current crosses its target, and float it from there.
    if (allowed && this.floating.some((k) => this.terms[k].drive.kind === 'I' && !this._bound(k))) {
      restart();
      const q = this._floatingContinuation(opts, r);
      if (q.converged) return q;
      r = q;
    }
    // A port held at a voltage across a capacitance (a transistor's gate): ramped there from its
    // zeroCharge, where the layout's neutral start is right. Cold, a thin-film transistor 0.8 V
    // below threshold failed, its channel laid out e^{−31} above where it ends, which a sweep of
    // the gate in 0.1 V steps reached (Fable round 3's hunter).
    // Several such gates: each starts at its zeroCharge, then they're ramped in turn, warm (two
    // gates both below threshold, each ramped alone with the other at its target, failed: the
    // hunter's re-check).
    const gates = allowed
      ? this.model.ports.flatMap((port, p) => {
        const d = this.terms[2 + p].drive;
        if (!port.capacitance || d.kind !== 'V' || d.R > 0) return [];
        const V = sourceAt(d.src, this.time), from = port.capacitance.zeroCharge;
        return V === from ? [] : [{ k: 2 + p, from, V }];
      })
      : [];
    if (gates.length > 0) {
      restart();
      // (with the contacts level while they ramp, where they can be, then the right one ramped
      // to its target from there: cold at the drain's voltage fails too)
      if (canContinue) this.sourceOverride.set(1, level);
      for (const g of gates) this.sourceOverride.set(g.k, g.from);
      let q = r;
      try {
        for (const [j, g] of gates.entries()) {
          q = this._continuation(opts, g.k, g.from, g.V, q, j > 0);
          // (_continuation clears its own override on return; this one is now held at its V)
          this.sourceOverride.set(g.k, g.V);
          if (!q.converged) break;
        }
      } finally {
        this.sourceOverride.delete(1);
        for (const g of gates) this.sourceOverride.delete(g.k);
        this._refreshSources();
      }
      if (q.converged && canContinue) q = this._continuation(opts, 1, level, target, q, true);
      if (q.converged) return solved(q);
      r = q;
    }
    if (canContinue) {
      // The right terminal's voltage ramped to its target: from the contacts' voltages at the last
      // converged solve when the state is that solution (a warm start), the left one first where
      // it moved, the right held where it was; otherwise from level terminals, consistent with a
      // cold start.
      restart();
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
      this._putLevels(dimmer);
    }
    // (solved() too: without it, a stepped state solved here stayed `stepped`, and its solution read
    // a waveform's voltage from before the step it had stepped onto, 0.1 V for 0; Astra's fourth pass)
    return solved(after(this._steadySteps(opts), r));
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
  // (With `line`, terminal k behind a resistance: the current its load line passes at V,
  // (V_src − V)/R, the target, the march from the first held terminal's level.)
  _floatingContinuation(opts, r, k = this.floating.find((j) => this.terms[j].drive.kind === 'I' && !this._bound(j)), line = null) {
    const drive = this.terms[k].drive, target = line ? NaN : sourceAt(drive.src, this.time), want = (V) => (line ? line(V) : target);
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
      return q.converged ? this.termI[k] - want(V) : NaN;
    };
    const fail = () => ({ converged: false, steps, iterations, history });
    const seen = []; // [V, I] of each held solve that converged
    const heldK = this.terms.findIndex((t) => t.drive.kind === 'V' && !(t.drive.R > 0));
    let Va = line ? (heldK >= 0 ? this.termV[heldK] : 0) : this.termV[k], fa, Vb, fb;
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
      // (A held solve that fails is retried at a quarter of the stride, warm from the last that
      // converged, to 1 mV: ending the march there, one failed sample of an iron couple at
      // 3.12 V sent it the other way to −23 V, and the solve failed claiming a limiting current
      // where 3.38 V passes it; Astra's fourth pass.)
      const start = Va, f0 = fa;
      // (Each way that shows a limit: marched to ±20 V with its current no nearer the target
      // than before, or stopped by its held solves after its current had turned away from the
      // target, past an extremum. Stopped while still nearing it, a march says nothing of a
      // limit; nor does one at ±20 V still nearing it, which says only that the voltage that
      // passes it, if any, is further out: an ohmic resistor driven at 100 V's current was
      // refused as "a limiting current"; Astra's A25, round 8.)
      let limited = 0, nearing = null;
      [Vb, fb] = [Va, fa];
      // (Passed exactly at the start, no march: one whose neighbours failed was refused, where
      // the start was the answer; Astra's A25, round 8.)
      for (const dir of f0 === 0 ? [] : [fa < 0 ? 1 : -1, fa < 0 ? -1 : 1]) {
        [Va, fa] = [start, f0];
        if (dir !== (f0 < 0 ? 1 : -1)) hold(start); // (back to the start, for a warm march)
        let dV = 0.02, nearest = Math.abs(fa);
        for (;;) {
          const before = this._levels();
          Vb = Va + dir * dV;
          fb = hold(Vb);
          if (!Number.isFinite(fb)) {
            this._putLevels(before);
            dV /= 4;
            if (dV < 1e-3) {
              if (Math.abs(fa) > nearest) limited++;
              break;
            }
            continue;
          }
          if (Math.sign(fb) !== Math.sign(fa) || fa === 0) break;
          if (Math.abs(Vb - start) >= 20) {
            // (Still nearing: nearer than before, to a part in a million, at a pace that, kept up,
            // passes it within ten times as far again. A reverse-biased diode's current, still
            // growing by 1e-5 of the gap a stride at 26 V, would pass 100 A/m² only past 1e8 V:
            // that's a limit.)
            const ahead = (Math.abs(fb) * Math.abs(Vb - Va)) / (Math.abs(fa) - Math.abs(fb));
            if (Math.abs(fb) < (1 - 1e-6) * nearest && ahead <= 10 * Math.abs(Vb - start)) nearing = { at: Vb, guess: Vb + dir * ahead };
            else limited++;
            break;
          }
          [Va, fa] = [Vb, fb];
          nearest = Math.min(nearest, Math.abs(fa));
          dV *= 1.5;
        }
        if (Number.isFinite(fb) && (Math.sign(fb) !== Math.sign(fa) || fa === 0)) break;
      }
      if (!Number.isFinite(fb) || (Math.sign(fb) === Math.sign(fa) && fa !== 0)) {
        // What the held voltages passed, the range of voltages and the range of currents (not
        // the currents at the two extreme voltages: a peak of 10 A/m² between them read "0 to 0"):
        // only where both ways showed a limit, or reached ±20 V still nearing the target (nearing:
        // the voltage it was still nearing at, and where its last stride's pace would pass it).
        if (seen.length > 1 && !line && limited + (nearing === null ? 0 : 1) === 2) {
          const Vs = seen.map(([V]) => V), Is = seen.map(([, I]) => I);
          this.unreached = { terminal: this.terms[k].name, target, V: [Math.min(...Vs), Math.max(...Vs)], I: [Math.min(...Is), Math.max(...Is)], nearing };
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
    // (Illinois) within the bracket, to 1e-10 of the target, or at open circuit of the larger
    // current at the bracket's ends (a target of 0 is never met exactly). (Of every current the
    // march saw, 7e9 A/m² at its far end, −5.3e-7 was accepted for a target of +1e-12; Astra's
    // fourth pass.)
    const scale = Math.max(line ? 0 : Math.abs(target), Math.abs(fa + want(Va)), Math.abs(fb + want(Vb)));
    try {
      let side = 0;
      for (let it = 0; it < 60; it++) {
        // (both ends passing the target exactly, fb − fa = 0: the end itself, not NaN; Astra's
        // A25, round 8)
        const Vm = fb === fa ? Va : Va - (fa * (Vb - Va)) / (fb - fa), fm = hold(Vm);
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

  // One terminal behind a resistance: its load line's crossing by held solves. Otherwise (or where
  // that fails), each terminal behind a resistance held at its source voltage (R taken off; at
  // V_src − I R where Kirchhoff gives I), solved, then the resistances put back and the steady
  // system solved from that state.
  _resistanceContinuation(opts, r) {
    const ks = this.floating.filter((k) => this.terms[k].drive.kind === 'V' && this.terms[k].drive.R > 0);
    const drives = ks.map((k) => this.terms[k].drive);
    let steps = r.steps ?? 0, iterations = r.iterations ?? 0;
    const history = (r.history ?? []).slice();
    // Where every other terminal is driven by a current, Kirchhoff gives the one behind R its
    // current, and so the voltage it ends at, V_src − I R: held there, putting R back moves
    // nothing (held at V_src, it had to slide the whole device by IR, and failed past volts).
    const others = this.terms.every((t, j) => ks.includes(j) || t.drive.kind === 'I');
    const at = (k) => {
      const d = this.terms[k].drive, V = sourceAt(d.src, this.time);
      if (!(others && ks.length === 1)) return V;
      const I = -this.terms.reduce((t, tm, j) => (j === k ? t : t + sourceAt(tm.drive.src, this.time)), 0);
      return V - I * d.R;
    };
    const start = this._levels();
    const restore = () => {
      this._putLevels(start);
    };
    // One terminal behind R, the rest not all driven by currents (so its current isn't known):
    // its crossing with the load line, (V_src − V)/R, by held solves, as an open circuit's is
    // with its target, first. Held at its source instead, the state is I R from the answer and
    // the pseudo-transient relaxes it over R C at Newton's pace through that swing (a lit Schottky
    // diode passing 0.9 A/m²: 0.05 s at 1 Ω·m², 10 s at 1e4, 105 s at 1e5, steps ∝ R).
    const d = drives[0], Vs = sourceAt(d.src, this.time), line = (V) => (Vs - V) / d.R;
    const byLine = ks.length === 1 && !others;
    if (byLine) {
      const q = this._floatingContinuation(opts, r, ks[0], line);
      if (q.converged) return q;
      restore(); // (from where this stage began, not where the failed one left it)
      steps = q.steps ?? steps;
      iterations = q.iterations ?? iterations;
      history.splice(0, history.length, ...(q.history ?? []));
    }
    let heldOk = false;
    try {
      const held = ks.map(at);
      ks.forEach((k, j) => (this.terms[k].drive = { ...this.terms[k].drive, src: { value: held[j] }, R: 0 }));
      this.redrive();
      // From the layout, the layout at the held drive: the one a cold start lays out behind R
      // has the terminal at its held neighbour's voltage, its metal's carriers at that level (a
      // closed Cu | Zn cell held at its open circuit, −1.126 V, behind 1e-4 Ω·m², with Stern
      // layers: held there from the zinc laid out at 0 V, the solve failed, where the same cell
      // without R converged cold; the echem combos' case 38).
      if (this.coldStart && !this.stepped) {
        this._refreshSources();
        this.initFromComposition();
      }
      const q = this._steadyFromHere(opts);
      steps += q.steps ?? 0;
      iterations += q.iterations ?? 0;
      history.push({ resistancesOff: true, converged: q.converged });
      heldOk = q.converged;
    } finally {
      ks.forEach((k, j) => (this.terms[k].drive = drives[j]));
      this.redrive();
      captureLocal(this);
    }
    const q = heldOk ? this._steadySteps(opts) : { converged: false, steps: 0, iterations: 0, history: [] };
    const done = { ...q, steps: steps + q.steps, iterations: iterations + q.iterations, history: [...history, ...q.history] };
    if (q.converged || ks.length !== 1 || byLine) return done;
    // Where Kirchhoff's hold failed too, the load line.
    restore();
    return this._floatingContinuation(opts, done, ks[0], line);
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
        const before = this._levels();
        this.generationScale = next;
        q = this._steadySteps(sub);
        steps += q.steps;
        iterations += q.iterations;
        history.push({ generation: next, converged: q.converged });
        if (q.converged) {
          s = next;
          factor = Math.min(100, factor * factor);
        } else {
          this._putLevels(before);
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
        // (Not where an amount no row keeps is there to lose: laid out from c0, a closed Y + Z went
        // from the state's 1e-5 to c0's 1.03e-6; Astra's fourth pass.)
        if (!q.converged && this.unkept.length === 0) {
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
        if (!q.converged) return { ...r, converged: false, steps, iterations }; // (said, not inherited from r: a caller may pass one that converged)
      }
      while (V !== target) {
        // (and the target where the stride no longer moves V: a span of 5e-324 V, its eighth 0,
        // ramped without end; Astra's fourth pass)
        const next = Math.abs(target - V) <= Math.abs(dV) || V + dV === V ? target : V + dV;
        const before = this._levels();
        this.sourceOverride.set(k, next);
        q = this._steadySteps(sub);
        steps += q.steps;
        iterations += q.iterations;
        history.push({ continuation: next, converged: q.converged });
        if (q.converged) {
          V = next;
          dV *= 1.5;
        } else {
          this._putLevels(before);
          dV /= 4;
          // (and at 0: a span of 1e-320 V, its millionth 0, its stride quartered to 0 and held
          // there without end; Astra's A29, round 8)
          if (Math.abs(dV) < 1e-6 * Math.abs(target - level) || dV === 0) return { converged: false, steps, iterations, history };
        }
      }
      return { converged: true, steps, iterations, history, ramped: true };
    } finally {
      this.sourceOverride.delete(k);
      if (hold !== undefined) this.sourceOverride.delete(1 - k);
      this._refreshSources();
    }
  }

  // (Its clock, history and step put back on a throw too, as on a return: a throw from a custom
  // statistics inside a step left the clock at 44000 s, and the next advance() did nothing;
  // Astra's fourth pass.)
  _steadySteps(opts) {
    const time = this.time;
    try {
      return this._stepsToSteady(opts);
    } catch (e) {
      this.time = time;
      this.history = [];
      this.dtNext = undefined;
      throw e;
    }
  }

  _stepsToSteady({ maxSteps = 80, tol = 1e-11 } = {}) {
    const time = this.time;
    const tau = this.slowestTime();
    const direct = this._directSteady();
    const giant = direct ? Infinity : 1e6 * tau;
    // (A giant step that overflows is a steady solve without the storage term that keeps what no
    // row does: a species at D = 1e-315 made it ∞, and the solve converged with 18% of Y + Z lost;
    // Astra's fourth pass.)
    if (!direct && !(giant < Infinity)) {
      const [dMin] = this._dRange(), { species, materials } = this.model;
      const slow = materials.flatMap((mat) => species.flatMap((sp, i) => (mat.present[i] && mat.D[i] === dMin ? [`materials.${mat.name}.species['${sp.name}'].D`] : [])));
      throw new DeviceError(`${slow[0] ?? 'a graded D'}: ${dMin} m²/s makes the steady solve's giant time step (1e6 of the slowest diffusion time, L²/D) overflow, and without it nothing keeps the amounts no row of the steady system holds. Make that species immobile (D: 0), or give it a D the device can wait out (above ~1e-290 m²/s).`);
    }
    let dt = giant, grow = 10; // (after a failure, dt grows by less, then back up to ×10)
    let triedDirect = 0;
    let totalIter = 0, steps = 0, converged = false;
    const history = [];
    // What no row keeps (laws.js's unkept) a giant step keeps only by its storage term, which can
    // fall below round-off against the fast species' conduction, and Newton then converges
    // anywhere along it: a closed Y + Z, Y still in part of its stretch, converged 1e8× too large.
    // So each such amount is checked after every step of dt ≥ giant, and a step that moved it
    // fails.
    // (From when the public solve began, not this call: a continuation's cold retry laid Y + Z out
    // again from c0, and the next calls each took that as theirs, keeping 1.03e-6 for 1e-5; Astra's
    // fourth pass.)
    const unkept = !direct && this.unkept.length > 0 ? this._unkeptAmounts(this.solveAmounts ?? undefined) : null;
    const drifted = () => {
      const now = this._unkeptAmounts(), lost = now.findIndex(([a], j) => !(Math.abs(a - unkept[j][0]) <= 1e-6 * unkept[j][1]));
      if (lost >= 0) this.drifted = now[lost][2];
      return lost >= 0;
    };
    const motions = []; // (each giant step's largest potential change: see settling, below)
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
      // (the floating terminals' motion over the step too, not Newton's last update to them: a
      // gate behind R ~1e18 Ω·m², its RC 6e5 s, moved 0.26 VT in its one giant step, its last
      // update 3e-15, and was accepted far from its source's 1 V; Astra's fourth pass)
      const motion = dt >= giant ? maxPotentialStep(this, this._diff(), Float64Array.from(this.floating, (k) => this.termV[k] - this.termVPrev[k])) : NaN;
      if (dt >= giant) motions.push(motion);
      if (dt >= giant && motion < tol) {
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
        // Then the giant step; or, solved directly, on to 1e12 τ, the direct solve tried once a
        // decade (above): jumping to dt = ∞ at 1e6 τ, its failure sent dt back to 1e-6 τ, round
        // and round below 70 s, where traps filled by 1e9 cm⁻³ electrons took 1e3 s (τ 7e-5 s).
        if (dt > (direct ? 1e12 : 1e6) * tau) dt = giant;
      } else if (dt < 1e6 * giant) {
        // conserved amounts present and still moving: let dt keep growing (capped; and finite,
        // where 1e6·giant overflows: a giant of 1e303 s grew to ∞, a step without the storage
        // term that keeps what no row does; Astra's A29, round 8)
        dt = Math.min(dt * 10, Number.MAX_VALUE);
      }
    }
    // Out of steps while each one still moved the state less than the last (a mode slower than
    // the capped step, shrinking by 1/(1 + k dt) a step): say so, and that more steps finish it.
    // (It took the generic "nothing specific showed why"; Fable round 3's review of src/steady.js.)
    const [m0, m1, m2] = motions.slice(-3);
    // (a steady ratio, below 0.95: a round-off floor or a cycle can show three falling values)
    if (!converged && steps >= maxSteps && m2 < m1 && m1 < m0 && m2 / m1 < 0.95 && Math.abs(m2 / m1 - m1 / m0) < 0.2) this.settling = { maxSteps, ratio: m2 / m1, motion: m2, tol };
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
  // the sum of its parts' sizes, the names of its parts], from the columns' amounts (now, or as
  // _columnAmounts took them before).
  _unkeptAmounts(cols = this._columnAmounts()) {
    const S = this.stretches.length, X = S + this.surfCols.length, { species } = this.model;
    return this.unkept.filter((w) => w.every((x, k) => k < X || x === 0)).map((w) => {
      let a = 0, size = 0;
      const names = [];
      w.forEach((x, k) => {
        if (x === 0) return;
        const v = x * cols[k];
        a += v;
        size += Math.abs(v);
        names.push(k < S ? species[this.stretches[k].species].name : this.model.surfaces[this.surfCols[k - S][0]].species[this.surfCols[k - S][1]].name);
      });
      return [a, size, names];
    });
  }

  // Each stretch's amount and each surface column's, in column order (laws.js's [0, X)).
  _columnAmounts() {
    const S = this.stretches.length;
    return Float64Array.from({ length: S + this.surfCols.length }, (_, k) => (k < S ? this.amount(this.stretches[k]) : this.surfaceAmount(k - S)));
  }

  // Huge steps pin each conserved amount only through a tiny storage term, so round-off can let
  // it creep. Shift each spectator's level uniformly to restore its amount exactly (in one step
  // for ideal statistics, where c ∝ e^η; by Newton on the shift otherwise).
  // (Swept again while a shift by Newton left another off: species sharing a lattice's sites
  // move each other's amounts, and in one pass B, listed before A, was left 4.5% off after A's
  // shift; Astra's fourth pass. An ideal one's shift moves only its own: one pass, as before.)
  _renormalizeSpectators() {
    const { n, M, u } = this;
    let changed = false, coupled = false;
    // (A floating conductor holds its charge on its faces, not as an amount: the huge steps keep
    // it, and there's nothing here to shift.)
    const spectators = this.stretches.flatMap((st, k) => (st.spectator && !this.model.materials[this.model.regions[st.regions[0]].material].conductor ? [k] : []));
    const sweep = (st, k) => {
      const want = this.referenceAmounts[k];
      // (one shift of η scales an ideal stretch's amount exactly, but not what a diffuse layer
      // at its edge holds, whose excess isn't proportional to c: iterate there too)
      const ideal = this.nodeIdeal.subarray(st.nodes[0], st.nodes[1] + 1).every((v) => v === 1) && this._diffuseHeld(st, this.gam) === 0;
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
          // (a layer's excess Γ_s by every concentration at its edge node, each moving with η_s:
          // Σ_j ∂Γ_s/∂c_j ∂c_j/∂η_s. On a lattice the others move too, and ∂Γ_s/∂c_s c_s alone
          // had the wrong sign there; Argon's review of 1.6.1.)
          const t = this.dV2 ?? (this.dV2 = new Float64Array(M));
          this.diffList.forEach((dl) => {
            if (!(dl.g >= st.nodes[0] && dl.g <= st.nodes[1])) return;
            let dG = 0;
            for (let j = 0; j < n; j++) {
              const gj = dl.dl.dgdc[st.species * n + j];
              if (gj === 0 || !this.present[dl.g * n + j]) continue;
              this._dc(dl.g, j, t);
              dG += gj * t[1 + st.species];
            }
            dA += this.model.grid.area[dl.g] * dG;
          });
          shift = (want - now) / dA;
          if (!(Math.abs(shift) > 1e-16)) return;
        }
        for (let g = st.nodes[0]; g <= st.nodes[1]; g++) u[this.blockOfNode[g] * M + 1 + st.species] += shift;
        changed = true;
        if (!ideal) {
          coupled = true;
          this.computeConcentrations();
        }
      }
    };
    const kept = () => spectators.every((k) => !(Math.abs(this.amount(this.stretches[k]) - this.referenceAmounts[k]) > 1e-12 * Math.abs(this.referenceAmounts[k])));
    for (let pass = 0; pass < 50; pass++) {
      coupled = false;
      for (const k of spectators) sweep(this.stretches[k], k);
      if (!coupled) break;
      if (kept()) break;
    }
    // Coupled tightly, the sweeps creep: P and Q filling a lattice to 4.999 of its 10 sites each
    // move each other almost one for one, and after the 50 passes Q, swept first, was still
    // 7.0e-5 high (P first: exact). Then a joint Newton on every spectator's shift at once
    // (Astra's re-read of the steady solve, round 8).
    if (coupled && !kept()) {
      const m = spectators.length, J = new Float64Array(m * m), r = new Float64Array(m), x = new Float64Array(m);
      const t = this.dV2 ?? (this.dV2 = new Float64Array(M)), sts = spectators.map((k) => this.stretches[k]);
      const within = (st, g) => g >= st.nodes[0] && g <= st.nodes[1];
      for (let it = 0; it < 30 && !kept(); it++) {
        J.fill(0);
        sts.forEach((st, a) => {
          r[a] = this.referenceAmounts[spectators[a]] - this.amount(st);
          // (∂ amount_a / ∂ shift_b: Σ vol ∂c_a/∂η_b over the nodes they share, and a diffuse
          // layer's Σ_j ∂Γ_a/∂c_j ∂c_j/∂η_b at its edge node, as the sweep reads its own)
          for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
            this._dc(g, st.species, t);
            sts.forEach((sb, b) => {
              if (within(sb, g)) J[a * m + b] += this.model.grid.vol[g] * t[1 + sb.species];
            });
          }
          this.diffList.forEach((dl) => {
            if (!within(st, dl.g)) return;
            for (let j = 0; j < n; j++) {
              const gj = dl.dl.dgdc[st.species * n + j];
              if (gj === 0 || !this.present[dl.g * n + j]) continue;
              this._dc(dl.g, j, t);
              sts.forEach((sb, b) => {
                if (within(sb, dl.g)) J[a * m + b] += this.model.grid.area[dl.g] * gj * t[1 + sb.species];
              });
            }
          });
        });
        solveFlat(J, r, x, m);
        if (!x.subarray(0, m).every(Number.isFinite)) break;
        sts.forEach((st, b) => {
          for (let g = st.nodes[0]; g <= st.nodes[1]; g++) u[this.blockOfNode[g] * M + 1 + st.species] += x[b];
        });
        this.computeConcentrations();
      }
    }
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
  // reaction stands on its affinity where the steady system is solved directly, so there it's
  // only caught here if that couldn't be arranged; by giant steps, it's checked like the rest.)
  stuckReactions({ all = true } = {}) {
    const { n, M, u, uLo, model } = this, grid = model.grid, S = this.stretches.length, X = S + this.surfCols.length, out = [], direct = this._directSteady();
    // Where no reservoir reaches (no contact, port or flow brings any of it in, and no reaction
    // trades with a fixed participant: fixedA = 0), a steady state is equilibrium: Σ η_i (∇·J_i
    // − Σ ν r) = 0 over such a part, by parts, is its dissipation, Σ J·(−∇η) + Σ r a, each term
    // ≥ 0, so each is 0. So every reaction in it is at equilibrium wherever it acts, and every
    // face's link passes nothing. The checks of a reaction alone (needed, below) missed two: a
    // cycle U = V = W = U, any two spanning what the third changes, at kf 1e-32 left at its c0
    // (1, 2, 9; its equilibrium 4, 4, 4); and P = Q in two pools across a blocked face, one at
    // equilibrium hiding the other at 1.9 RT (Astra's re-read of the steady solve, round 8).
    // (Groups: the stretches, joined by what reacts them, a face's surface species too.)
    const group = this.stretches.map((_, k) => k), root = (k) => (group[k] === k ? k : (group[k] = root(group[k])));
    const join = (ks) => ks.forEach((k) => (group[root(k)] = root(ks[0])));
    const surfNode = new Map(), surfOf = (f, q) => surfNode.get(`${f}:${q}`) ?? surfNode.set(`${f}:${q}`, group.push(group.length) - 1).get(`${f}:${q}`);
    const opens = [];
    this.stretches.forEach((st, k) => {
      if (st.leftOpen || st.rightOpen || st.ports.length > 0 || model.regions.slice(st.regions[0], st.regions[1] + 1).some((reg) => reg.velocity !== 0)) opens.push(k);
    });
    model.regions.forEach((reg, q) => this.rxs.forEach((rx, x) => {
      if (!this.runsIn[q][x]) return;
      const ks = Array.from(rx.sp, (i) => this.stretchOf[q * n + i]);
      join(ks);
      if (rx.fixedA !== 0 || model.reactions[x].generation) opens.push(ks[0]);
    }));
    const faceNodes = (f, rx) => [...rx.part.map((p) => this.stretchOf[(f + p.side) * n + p.i]), ...rx.surf.map((p) => surfOf(f, p.s))];
    model.interfaces.forEach((itf, f) => itf.reactions.forEach((rx) => {
      const ks = faceNodes(f, rx).filter((k) => k >= 0);
      if (ks.length === 0) return;
      join(ks);
      if (rx.fixedA !== 0 || ks.length < rx.part.length + rx.surf.length) opens.push(ks[0]);
    }));
    const open = new Set(opens.map(root));
    const closed = (k) => all && k >= 0 && !open.has(root(k));
    const shut = (where, as) => {
      const off = as.filter((a) => Math.abs(a) > 1e-6);
      if (off.length > 0) out.push({ where, affinity: off.reduce((m, a) => (Math.abs(a) > Math.abs(m) ? a : m)), closed: off.length });
    };
    const passes = (k) => !noCurrent(this.terms[k], this.time);
    const fed = (w) => w.some((x, k) => x !== 0 && ((k >= X && passes(this.carriers[k - X])) || (k < S && this.terms.some((t, j) => t.drive.kind === 'I' && passes(j) && feeds(j, this.stretches[k])))));
    const needed = (rx) => {
      const basis = conservedMoieties(this._structure(), { drives: true, skip: rx });
      return basis.some((w) => !fed(w) && basis.skipped.some((row) => Math.abs(row.reduce((t, x, k) => t + x * w[k], 0)) > 1e-9));
    };
    const eta = (g, slot) => u[this.blockOfNode[g] * M + slot] + uLo[this.blockOfNode[g] * M + slot];
    const report = (where, as) => {
      if (as.length === 0 || !(as.every((a) => a > 1e-6) || as.every((a) => a < -1e-6))) return false;
      out.push({ where, affinity: as.reduce((m, a) => (Math.abs(a) < Math.abs(m) ? a : m)) });
      return true;
    };
    model.interfaces.forEach((itf, f) => {
      const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
      itf.reactions.forEach((rx, q) => {
        // (an affinity row holds it at equilibrium only where the steady system is solved
        // directly; by giant steps, it's never enforced: a face reaction at 1e-30 left its c0,
        // converged, unwarned; Fable round 3's review of src/steady.js)
        const gS = itf.surfaceSide ? gR : gL;
        const bf = this.blockOfFace[f], w = (p) => { const s = planeSlot(this, f, p); return s < 0 ? 0 : u[bf * M + s] + uLo[bf * M + s]; }; // (at a layer's plane)
        const a = rx.surf.reduce((t, p) => t - p.nu * eta(gS, 1 + n + p.s), rx.part.reduce((t, p) => t - p.nu * (eta(p.side ? gR : gL, 1 + p.i) + w(p)), rx.fixedA));
        const where = `interfaces[${f}].reactions[${q}]`;
        if (!(direct && this.affinities.some((af) => af.rx === rx)) && needed(rx) && report(where, [a])) return;
        if (closed(faceNodes(f, rx)[0] ?? -1)) shut(where, [a]);
      });
    });
    // Where a participant can't move (D = 0, inside a region without flow or mixing, away from its
    // faces and any port) and no other bulk reaction makes or uses it, each node holds its own
    // amount, so the reaction must be at equilibrium at that node by itself. Out of equilibrium one
    // way in one region and the other way in the next (Y immobile in the second), the sum of a
    // frozen reaction's rates over the device was balanced and its affinities of both signs: it
    // passed for steady (Fable round 3's review of src/steady.js).
    // (A port moves a species at its nodes only if it links it or reacts it: an electron-only port
    // over the region switched the check off for every species there, and an immobile Y left 50×
    // off converged unwarned; Astra's fourth pass.)
    const inPort = new Uint8Array(this.nNodes * n);
    for (const port of model.ports) {
      for (let i = 0; i < n; i++) {
        if ((port.species[i]?.type ?? 'blocked') === 'blocked' && !port.reactions.some((r) => r.part.some((p) => p.side === 0 && p.i === i))) continue;
        for (const g of port.nodes) inPort[g * n + i] = 1;
      }
    }
    const held = (x, g) => {
      const reg = model.regions[grid.nodeRegion[g]], mat = model.materials[reg.material];
      if (reg.velocity !== 0 || reg.mixing !== 0 || g === grid.regionStart[grid.nodeRegion[g]] || g === grid.regionEnd[grid.nodeRegion[g]]) return false;
      const rx = this.rxs[x];
      return Array.from(rx.sp).some((i) => !inPort[g * n + i] && !(mat.D[i] > 0) && this.rxs.every((o, y) => y === x || !(o.kfNode[g] > 0) || !o.sp.includes(i)));
    };
    model.reactions.forEach((mrx, x) => {
      // (Generations too: what's classed as one, made from or turned into fixed reservoirs only,
      // includes water formed from H⁺ and OH⁻, which a closed solution needs at equilibrium; a
      // frozen one passed for steady at 4.6 RT of affinity, check() saying "equilibrium" (Fable
      // round 5). Light that's balanced by recombination isn't alone in what it changes.)
      const rx = this.rxs[x];
      const as = [], own = [], inClosed = [];
      for (let g = 0; g < this.nNodes; g++) {
        if (!(rx.kfNode[g] > 0 || rx.srh?.[this.nodeMaterial[g]])) continue;
        const a = rx.sp.reduce((t, i, p) => t + rx.nu[p] * eta(g, 1 + i), rx.fixedA);
        as.push(a);
        if (held(x, g)) own.push(a);
        if (closed(this.stretchOf[grid.nodeRegion[g] * n + rx.sp[0]])) inClosed.push(a);
      }
      if (needed(mrx) && report(`bulkReactions[${x}]`, as)) return;
      const off = own.filter((a) => Math.abs(a) > 1e-6);
      if (off.length > 0) out.push({ where: `bulkReactions[${x}]`, affinity: off.reduce((m, a) => (Math.abs(a) > Math.abs(m) ? a : m)), immobile: off.length });
      else shut(`bulkReactions[${x}]`, inClosed);
    });
    model.ports.forEach((port, k) => {
      port.reactions.forEach((rx, q) => {
        if (!needed(rx)) return;
        let aV = rx.fixedA;
        for (const p of rx.part) if (p.side === 1) aV -= p.nu * this.portEta(port, p.i);
        report(`ports[${k}].reactions[${q}]`, Array.from(port.nodes, (g) => rx.part.reduce((a, p) => (p.side === 1 ? a : a - p.nu * eta(g, p.side === 2 ? 1 + n + p.s : 1 + p.i)), aV)));
      });
    });
    // A neutral species that nothing makes or uses, in regions without flow, fed only by contacts'
    // rate links (exchange, velocity): a steady state has no net flow through them, so not all of
    // them can push it the same way. A link too slow for the giant steps (an exchange at
    // k ≤ 1e-23 mol/(m²·s)) left it at its c0, 0.69 RT from the reservoir, converged, unwarned
    // (Astra's fourth pass). Its affinity: η outside − η at the contact's node.
    for (const st of this.stretches) {
      const i = st.species, [r0, r1] = st.regions;
      if (this.z[i] !== 0 || st.reactive || st.linked.length > 0 || !(st.leftOpen || st.rightOpen)) continue;
      if (model.regions.slice(r0, r1 + 1).some((reg) => reg.velocity !== 0)) continue;
      const sides = [...(st.leftOpen ? ['left'] : []), ...(st.rightOpen ? ['right'] : [])], links = sides.map((side) => model.contacts[side].species[i]);
      if (!links.every((l) => l.type === 'exchange' || l.type === 'velocity')) continue;
      const where = sides.map((side) => `contacts.${side}.species['${model.species[i].name}']`).join(', ');
      const as = sides.map((side, j) => (links[j].type === 'exchange' ? links[j].mu / model.RT : this.contactEta(side, i)) - eta(side === 'left' ? 0 : this.nNodes - 1, 1 + i));
      if (report(where, as) || !all) continue;
      // Pushing both ways, their rates summed, against their slope (as rateHeld does where the
      // steady system is solved directly): two exchange contacts at k = 1e-32 toward 0.5 and 8,
      // on giant steps, left W at its c0 of 1 or 3 unwarned, the two through-fluxes 3× apart
      // (Astra's re-read of the steady solve, round 8).
      let q = 0, slope = 0;
      sides.forEach((side, j) => {
        const g = side === 'left' ? 0 : this.nNodes - 1, A = grid.area[g], l = links[j], x = Math.min(as[j], 700);
        const [r, d] = l.type === 'exchange' ? [l.k * x, l.k] : [l.v * this.c[g * n + i] * Math.expm1(x), l.v * this.c[g * n + i] * Math.exp(x)];
        q += A * r;
        slope += A * d;
      });
      if (Math.abs(q / slope) > 1e-6) out.push({ where, affinity: q / slope, links: sides.length });
    }
    // A face's rate link (a permeability, a conductance) inside a stretch that nothing makes or
    // uses, with no port and at most one end open: its flux is zero everywhere in a steady
    // state, so each such face's levels meet. W at 1 and 7 in two closed halves behind
    // P = 1e-35 m/s was taken as steady from its first giant step (4 and 4 at any P > 0;
    // Astra's re-read of the steady solve, round 8).
    for (const st of this.stretches) {
      if (!all || st.reactive || st.ports.length > 0 || (st.leftOpen && st.rightOpen)) continue;
      const i = st.species;
      for (let f = st.regions[0]; f < st.regions[1]; f++) {
        const type = model.interfaces[f].links[i].type;
        if (type !== 'permeability' && type !== 'conductance') continue;
        const d = eta(grid.regionStart[f + 1], 1 + i) - eta(grid.regionEnd[f], 1 + i);
        if (Math.abs(d) > 1e-6) out.push({ where: `interfaces[${f}].species['${model.species[i].name}']`, affinity: d, face: true });
      }
    }
    // A terminal behind a resistance that passes current only by charging (a gate, a capacitive
    // contact with no ions through it): a steady state passes none, so it sits at its source. An
    // RC too long for the giant steps (R = 1e30 Ω·m², C = 1 µF/m²) moved it less than their
    // tolerance and was taken as steady at −0.18 V for a 1 V source (Astra's re-read of the steady
    // solve, round 8). Its affinity: (V − V_src)/V_T.
    for (const k of this.floating) {
      const d = this.terms[k].drive;
      if (!all || d.kind !== 'V' || !(d.R > 0) || !chargeOnly(this, k)) continue;
      const off = (this.termV[k] - sourceAt(d.src, this.time)) / this.VT;
      if (Math.abs(off) > 1e-6) out.push({ where: this.terms[k].name, affinity: off, charging: true });
    }
    return out;
  }
}
