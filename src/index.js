// driftlet: 1D drift–diffusion–reaction solver.

import { normalizeDevice, DeviceError } from './device.js';
import { Solver, SolverError } from './solver.js';
import { makeSolution } from './solution.js';
import { deepEqual, markWhole, merge, WHOLE } from './merge.js';
import { carryState } from './carry.js';
import { checkOptions } from './errors.js';

export { DeviceError } from './device.js';
import { normalizeDrives } from './device.js';

// (a contact the kit wrote stays whole through set()'s own copies of it: see merge.js)
const keepWhole = (from, to) => (from[WHOLE] ? markWhole(to) : to);
export { SolverError } from './errors.js';
export * from './constants.js';

// Whether two definitions are the same (or, with `terminal`, the same but for a terminal's
// drive). Functions (custom statistics) compare by identity.
const DRIVE = new Set(['V', 'I', 'R']);
const eq = (x, y, terminal) => {
  if (x === y) return true;
  if (typeof x !== 'object' || typeof y !== 'object' || x === null || y === null) return false;
  if (Array.isArray(x) !== Array.isArray(y)) return false;
  const keys = new Set([...Object.keys(x), ...Object.keys(y)]);
  for (const k of keys) {
    if (terminal && DRIVE.has(k)) continue;
    if (!eq(x[k], y[k], false)) return false;
  }
  return true;
};
// Whether two definitions differ only in their terminals' drives (V, I, R on the contacts and
// ports).
function sameExceptDrives(a, b) {
  const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
  for (const k of keys) {
    if (k === 'contacts') {
      const ca = a.contacts ?? {}, cb = b.contacts ?? {};
      for (const side of new Set([...Object.keys(ca), ...Object.keys(cb)])) {
        if (side !== 'left' && side !== 'right') return false; // a misspelt side: rebuild, which says so
        if (!eq(ca[side] ?? {}, cb[side] ?? {}, true)) return false;
      }
    } else if (k === 'ports') {
      const pa = a.ports ?? [], pb = b.ports ?? [];
      if (pa.length !== pb.length || pa.some((p, i) => !eq(p, pb[i], true))) return false;
    } else if (!eq(a[k], b[k], false)) return false;
  }
  return true;
}

// A definition's contacts as the drive-only path can take them: absent, or a plain object (an
// array or a number goes through a rebuild, which refuses it; Astra's review).
const plainContacts = (def) => def.contacts === undefined || (typeof def.contacts === 'object' && def.contacts !== null && !Array.isArray(def.contacts));

// A state that can't carry over (the species changed): the device restarts from its regions' c0,
// at the sources then, the clock carrying on, its terminals' running charges by name; said in
// every solution up to the next steady solve.
function restartFrom(old, s) {
  s.restarted = `set(): the state couldn't carry over (the species changed), so the device restarted from its regions' c0 at t = ${old.time} (the clock carries on): what follows is from there, not from the state before set(). (Said in every solution up to the next steady solve, which keeps a closed region's restarted amounts.)`;
  // (The clock carries on, so a run's next advance(t) is still at t: it read 0, and an
  // advance(20.001) after a change at t = 20 ran 20 s from c0. And its start is at the sources
  // then: built at t = 0's, a ramped contact restarted 50 mV off; Astra's review.)
  s.time = old.time;
  s.sourceTime = old.sourceTime;
  s._refreshSources();
  s.initFromComposition();
  s.referenceAmounts = s.stretches.map((st) => s.amount(st));
  s.surfaceRef = Float64Array.from(s.surfCols, (_, j) => s.surfaceAmount(j));
  s.terms.forEach((t, k) => {
    const o = old.terms.findIndex((u) => u.name === t.name);
    if (o >= 0) s.termCharge[k] = old.termCharge[o];
  });
}

// Each terminal's drive (the contacts, then the ports), set in the model in place.
function applyDrives(model, drives) {
  model.terminals.forEach((t, k) => (t.drive = drives[k]));
  model.contacts.left.drive = drives[0];
  model.contacts.right.drive = drives[1];
  model.ports.forEach((p, k) => (p.drive = drives[2 + k]));
}

// A method's options: a name it doesn't take (a typo, dtmax for dtMax) is an error, not ignored.
const options = (opts, method, allowed) => checkOptions(method, opts, allowed, SolverError);

/** @typedef {import('./types.js').DeviceDefinition} DeviceDefinition */
/** @typedef {import('./types.js').Solution} Solution */

export class Device {
  /** @param {DeviceDefinition} def plain, serialisable device definition (see docs/device.md) */
  constructor(def) {
    this.def = def;
    this.model = normalizeDevice(def);
    this._solver = null;
    this._touched = false; // (solved or stepped: a state of its own to carry through set())
  }

  /** The solver, built on first use (so construction only validates the definition). */
  get solver() {
    if (!this._solver) this._solver = new Solver(this.model);
    return this._solver;
  }

  get grid() {
    return this.model.grid;
  }

  /**
   * Change part of the definition (deep-merged). A contact given V or I drops the other (and
   * given I, any series R), so `{ contacts: { right: { I: 0 } } }` switches it to open circuit. Ports are patched by name,
   * `{ ports: { gate: { V: 0.2 } } }` (an array replaces them all). A change to the terminals'
   * drives alone (V, I, R) is applied in place, cheaply. Otherwise the device is rebuilt, carrying the current state
   * over by what each value means (src/carry.js): onto a new grid, through new statistics, μ° or T; a change of the
   * species restarts it from the regions' c0.
   * @param {object} patch a partial device definition, merged into the current one
   * @returns {this}
   */
  set(patch) {
    // A terminal given a new kind of drive drops the old one: { I: 0 } replaces a held V.
    // (A current drives a terminal alone: given I, a contact drops V and any series R too.)
    const redrive = (c) => (c && typeof c === 'object' && !Array.isArray(c) ? ('I' in c && !('V' in c) ? keepWhole(c, { ...c, V: undefined, ...('R' in c ? {} : { R: undefined }) }) : 'V' in c && !('I' in c) ? keepWhole(c, { ...c, I: undefined }) : c) : c);
    for (const side of ['left', 'right']) {
      const c = patch?.contacts?.[side];
      if (c && typeof c === 'object') patch = { ...patch, contacts: { ...patch.contacts, [side]: redrive(c) } };
    }
    // Ports: an array replaces them all; an object patches them by name, { gate: { V: 0.2 } }.
    const ports = patch?.ports;
    if (ports !== undefined && !Array.isArray(ports)) {
      if (ports === null || typeof ports !== 'object') throw new DeviceError('set: ports must be an array (replacing them all) or an object of patches by port name, e.g. { gate: { V: 0.2 } }');
      const old = this.def.ports ?? [];
      for (const name of Object.keys(ports)) {
        if (!old.some((p, k) => (p.name ?? `ports[${k}]`) === name)) throw new DeviceError(`set: ports.${name}: no port named '${name}'${old.length ? ` (the ports: ${old.map((p, k) => p.name ?? `ports[${k}]`).join(', ')})` : ''}`);
      }
      patch = { ...patch, ports: old.map((p, k) => {
        const q = ports[p.name ?? `ports[${k}]`];
        if (q === undefined) return p;
        const m = merge(p, redrive(q));
        for (const key of ['V', 'I']) if (q[key] !== undefined) m[key] = q[key]; // (whole: see below)
        // (and a link given with its type, as at a contact: merged, a conductance link kept its
        // offset over one given without, and the port's current halved; Astra's review)
        for (const [name, link] of Object.entries(q.species ?? {})) if (link && typeof link === 'object' && 'type' in link) m.species = { ...m.species, [name]: link };
        return m;
      }) };
    }
    const def = merge(this.def, patch);
    // A drive is replaced whole, not merged: a waveform's keys are its own (merged, a ramp set
    // over a triangle kept the triangle's repeat; Fable's review of the front door).
    // Likewise a link or a φ law given with its type: a velocity link set over an equilibrium one
    // kept its offset, a pinned φ over a capacitive one its C (Fable's re-check). (A contact the kit wrote replaces the old one whole: merge's WHOLE.)
    for (const side of ['left', 'right']) {
      const c = patch?.contacts?.[side];
      if (c && typeof c === 'object') for (const key of ['V', 'I']) if (c[key] !== undefined) def.contacts[side][key] = c[key];
      for (const [name, link] of Object.entries(c?.species ?? {})) if (link && typeof link === 'object' && 'type' in link) def.contacts[side].species[name] = link;
      if (c?.phi && typeof c.phi === 'object' && 'type' in c.phi) def.contacts[side].phi = c.phi;
    }
    // Nothing changed (a page resending its drives every frame): nothing to do. (Redriven, a
    // transient's BDF2 history was cleared and a steady state forgot it was one; Astra's review.)
    if (deepEqual(def, this.def)) return this;
    // A regridded steady state not yet solved on its grid stands for the steady state at the old
    // definition: solve it there first (see Device.advance).
    this._settleRegrid();
    // A solver nothing has solved or stepped holds only its start, which the new definition's
    // replaces. (Kept, a drive patch after a bare solution() kept the old drive's cold start, and
    // the first step read −73846 A/m² where 26154 is right; Astra's review.)
    const live = this._solver && (this._touched || this._solver.restarted) ? this._solver : null;
    // Only the terminals' drives changed: update them in place, keeping the solver and its state
    // (a step in a source restarts the time stepping's order, as at any discontinuity).
    if (live && sameExceptDrives(this.def, def) && plainContacts(def)) {
      this._redefineDrives(def);
      return this;
    }
    // A rebuild, carrying the state by what each value means (src/carry.js). It changes the device
    // only once it has all succeeded: a change that fails leaves it as it was.
    const was = { def: this.def, model: this.model, solver: this._solver };
    try {
      // A new grid and anything else, on a steady state: the new grid first, the steady state
      // solved again there (the state it stands for), then the rest, as two set()s are. (In one,
      // the regrid was settled at the new definition, and the change's transient lost: a bias
      // stepped, −859 A/m² where −2.0e7 is right, Sol's review; a contact's level, a midpoint at 5
      // where 1 is right, Astra's.)
      const mid = { ...this.def, grid: def.grid };
      if (live?.atSteady && !deepEqual(this.def.grid, def.grid) && !deepEqual(mid, def)) {
        let ok = true;
        try {
          normalizeDevice(mid);
        } catch (err) {
          if (!(err instanceof DeviceError)) throw err;
          ok = false; // (the new grid doesn't fit the old device: in one)
        }
        if (ok) {
          this._rebuild(mid, live);
          this._settleRegrid();
          if (sameExceptDrives(mid, def) && plainContacts(def)) this._redefineDrives(def);
          else this._rebuild(def, this._solver);
          return this;
        }
      }
      this._rebuild(def, live);
    } catch (err) {
      this.def = was.def;
      this.model = was.model;
      this._solver = was.solver;
      throw err;
    }
    return this;
  }

  // Rebuild at def, carrying old's state (src/carry.js); a device with no state of its own just
  // takes the new definition.
  _rebuild(def, old) {
    const model = normalizeDevice(def);
    const solver = this._solver ? new Solver(model) : null;
    if (old) {
      // The same species (by name and charge, in any order) carry over; other species restart
      // the device from its regions' c0.
      const key = (m) => m.species.map((sp) => `${sp.name}:${sp.z}`).sort().join('\u0000');
      if (key(old.model) === key(model)) carryState(old, solver, { oldDef: this.def, newDef: def, pure: deepEqual({ ...this.def, grid: undefined }, { ...def, grid: undefined }) });
      else restartFrom(old, solver);
    }
    this.def = def;
    this.model = model;
    this._solver = solver;
  }

  // The definition def, which differs from the current one only in the terminals' drives: set
  // them in place. (Validated before anything changes.)
  _redefineDrives(def) {
    const drives = normalizeDrives(def, this.model);
    // A regridded steady state stands for the old drives' (see Device.advance): solve it now, before
    // they change. (Left for advance(), it was solved at the new drive, and the step response
    // was lost; Fable round 3's steady.js review.)
    this._settleRegrid();
    applyDrives(this.model, drives);
    this.def = def;
    this._solver.redrive();
  }

  /**
   * Steady state (or equilibrium), warm-started from the current state. Doesn't advance time.
   * @param {{ maxSteps?: number, tol?: number, continuation?: boolean }} [opts]
   * @returns {Solution}
   */
  solve(opts) {
    options(opts, 'solve', ['maxSteps', 'tol', 'continuation']);
    // A solve that fails reports the state it failed in, and leaves the device as it was: a
    // device with no steady state (a host filling at a current) can still advance() from there.
    const snap = this.solver._snapshot();
    const result = this.solver.solveSteady(opts);
    // (a failed solve's state isn't a solution: its accuracy warnings are left out)
    const sol = makeSolution(this.solver, result.converged ? result : { ...result, unsolved: true });
    if (!sol.converged) this.solver._restore(snap);
    else {
      this._touched = true;
      this.solver.regriddedSteady = false; // (solved on the new grid already: see Device.advance)
      this.solver.jumpedT = undefined; // (solved at the new temperature: see integrate in src/timeSteps.js)
    }
    return sol;
  }

  /**
   * Advance the transient by dt (s) with backward Euler. If Newton fails, the interval is
   * split into smaller steps automatically; `substeps` in the result says how many.
   * With `{ method: 'bdf2' }`, steps after the first use second-order BDF2.
   * @param {number} dt s
   * @param {{ method?: 'be' | 'bdf2' }} [opts]
   * @returns {Solution}
   */
  step(dt, opts) {
    options(opts, 'step', ['method']);
    // dt is a duration. Unchecked, step(-1) stepped back to t = −1 and converged, and step('1')
    // made the clock the string '-11' (the state audit). As advance() does with its tEnd, a NaN
    // (an animation's first frame) or 0 does nothing, and says so.
    const length = `dt is the step's length in s, a positive number (the device is at t = ${this.solver.time})`;
    if (typeof dt !== 'number' || dt < 0 || dt === Infinity) throw new SolverError(`step(${typeof dt === 'string' ? `'${dt}'` : dt}): ${length}${dt === Infinity ? '. For the long-time limit, solve() the steady state' : ''}`);
    if (!(dt > 0)) return makeSolution(this.solver, { converged: true, iterations: 0, substeps: 0, stopped: `step(${dt}) did nothing: ${length}` });
    this._settleRegrid(); // (as advance() does: stepped unsolved, 27 A/m² where 68.8 is right)
    this.solver._freshDiagnostics();
    return this._touching(() => makeSolution(this.solver, this.solver.advance(dt, opts)));
  }

  /**
   * Adaptive transient to time tEnd (s): variable-step BDF2 with local error control. With
   * `budgetMs`, it returns after that much wall time even if tEnd isn't reached (`done` says
   * which), so an animation can call it once per frame; the step size carries over. `probes`
   * read a species (or φ) at points inside the device after every accepted step, into
   * `trace.probes` (one array per probe, alongside `trace.t`): what a detector at x sees. At a
   * face, `{ interface, species }` reads a species' flux through it (mol/(m²·s), toward +x), and
   * `{ interface, gate }` a gate's fraction open.
   * @param {number} tEnd s
   * @param {{ tol?: number, dt0?: number, dtMax?: number, budgetMs?: number, maxSteps?: number, method?: 'bdf2' | 'be',
   *   probes?: ({ x: number, species?: string, quantity?: 'c' | 'V' | 'phi', region?: string | number }
   *     | { interface: number, species: string } | { interface: number, gate: string })[] }} [opts]
   * @returns {Solution}
   */
  advance(tEnd, opts) {
    options(opts, 'advance', ['tol', 'dt0', 'dtMax', 'budgetMs', 'maxSteps', 'method', 'probes']);
    // A steady state interpolated onto a new grid isn't quite one there (a strictly neutral
    // region's interpolated levels miss neutrality), and a transient's first short steps can't
    // close that gap: storage resists it (Fable round 2's cold user: set({ grid }) on a solved
    // cell, then advance(), failed down to 1e-17 s). Its steady state on the new grid is what it
    // stands for: solve that first. (Mid-run, an interpolated transient state steps on as it is.)
    this._settleRegrid();
    this.solver._freshDiagnostics();
    return this._touching(() => makeSolution(this.solver, this.solver.integrate(tEnd, opts)));
  }

  // A call that moved the clock has a state of its own to carry through set(). (Marked before it
  // ran, a step(0), a refused step or a failed first solve made c0 a "state" that a c0 patch
  // then couldn't replace; Astra's review.)
  _touching(call) {
    const t0 = this.solver.time;
    const out = call();
    if (this.solver.time !== t0) this._touched = true;
    return out;
  }

  // A steady state carried onto a new grid, solved there (once; see Device.advance).
  // (A settle that fails leaves the interpolated state, and says so in the next solutions: it
  // was dropped silently, the flag cleared before the snapshot that restored it; Astra's review.)
  _settleRegrid() {
    const s = this._solver;
    if (!s?.regriddedSteady) return;
    const snap = s._snapshot();
    if (s.solveSteady().converged) s.regriddedSteady = false;
    else {
      s._restore(snap);
      s.regriddedSteady = false;
      s.restarted = 'set(): the steady state carried onto the new grid could not be solved again there, so what follows goes on from the interpolated state. (Said in every solution up to the next steady solve.)';
    }
  }

  /**
   * Small-signal impedance Z(f) about the steady state (solved first), at one terminal: its
   * held voltage is perturbed (or its driven current), the others keep their drives.
   * Z = δV/δI with I into the device (Ω·m²), with the steady state's `warnings`.
   * @param {ArrayLike<number>} frequencies Hz
   * @param {{ terminal?: string, profiles?: boolean }} [opts] the terminal (default 'right'), and
   *   whether to return complex profiles per frequency
   * @returns {import('./types.js').ImpedanceResult}
   */
  impedance(frequencies, opts) {
    options(opts, 'impedance', ['terminal', 'profiles']);
    // The steady state, solved again even where the device holds one: the solve from there also
    // polishes it, and taken as it stood after a solve(), a stress electrode's Z at 1 µHz came
    // out non-passive, Re Z −2.3e3 of 1.6e6 (stress echem 133). (A steady solve that fails leaves
    // the device as it was, as solve() does: a dark cell driven past saturation was left at a held
    // stage's 26.9 V, and the advance() the error suggests went on from there; Fable round 3's
    // steady.js review.)
    const snap = this.solver._snapshot();
    const steady = this.solver.solveSteady();
    if (!steady.converged) {
      // (with why, as solve() would say it: a terminal driven past what any steady state passes,
      // say, where advance() is the way; Fable round 2's cold user, charging a battery)
      const why = makeSolution(this.solver, { ...steady, unsolved: true }).warnings;
      this.solver._restore(snap);
      throw new SolverError(`impedance: the steady state did not converge${why.length > 0 ? `: ${why[0].replace(/\.$/, '')}` : ''}. An impedance is about a steady state; about a transient one, drive a small sine with advance() and fit it`, { warnings: why });
    }
    this._touched = true; // (solved: a state of its own)
    // The steady state's warnings (an unresolved double layer, say) hold for the impedance too.
    const { warnings } = makeSolution(this.solver, steady);
    const out = this.solver.impedance(frequencies, opts);
    // Where the current a frequency carries is lost to round-off (far below the device's slowest
    // relaxation, or read as a small difference of large terms), Z is said to be unresolved: an
    // electrolyte between blocking electrodes at 1e-4 Hz read Re Z 1e4 times too large, and a
    // closed host's held terminal at 1e-7 Hz a negative one, each returned without a word; Re Z
    // as well, where its error estimate is over 10 times it (the same electrodes at 1e-3 Hz: 153×
    // too large, its estimate 170 times it; at 100 Hz, right, 2 times).
    const g = (x) => String(Number(x.toPrecision(3))), lost = [], rough = [];
    out.f.forEach((f, q) => {
      const e = out.relError[q], Z = Math.hypot(out.Z.re[q], out.Z.im[q]);
      if (!(e <= 1e-2)) lost.push(f);
      else if (!(10 * Math.abs(out.Z.re[q]) > e * Z)) rough.push(f); // (the estimate is cautious: right ones sat within 3× of it)
    });
    if (lost.length > 0) {
      warnings.push(
        `impedance: Z isn't resolved at ${g(lost[0])}${lost.length > 1 ? `…${g(lost.at(-1))}` : ''} Hz (${lost.length} of ${out.f.length} frequencies): its estimated error (relError) is over 1%, the current there lost to round-off ` +
          '(far below a capacitive device\'s relaxation, where almost no current flows, or where the current is a small difference of large terms). ' +
          'Measuring at another terminal, or with this one driven by a current (at its steady current, e.g. I: 0), can avoid it: a closed host at open circuit resolved so.',
      );
    }
    if (rough.length > 0) {
      warnings.push(`impedance: Re Z isn't resolved at ${g(rough[0])}${rough.length > 1 ? `…${g(rough.at(-1))}` : ''} Hz (${rough.length} of ${out.f.length} frequencies): it's under a tenth of its estimated error, relError·|Z| (a cautious estimate), so there read it only as small next to Im Z.`);
    }
    return { ...out, warnings };
  }

  /**
   * A checkpoint of the definition and solver state, for `_rollback()` (the live wrapper keeps
   * its last good state this way).
   * @internal
   */
  _checkpoint() {
    const solver = this.solver;
    // (the drives too: a change of drive alone is made in place, in the model this checkpoint keeps)
    return { def: this.def, model: this.model, solver, touched: this._touched, drives: this.model.terminals.map((t) => t.drive), snap: solver._snapshot() };
  }

  /**
   * Back to a checkpoint taken on this device.
   * @param {ReturnType<Device['_checkpoint']>} cp
   * @internal
   */
  _rollback(cp) {
    this.def = cp.def;
    this.model = cp.model;
    this._solver = cp.solver;
    this._touched = cp.touched;
    // The drives as they were (a drive changed since, then a rebuild, left the checkpoint's model
    // driven at the new value while its definition said the old one), and the solver's own
    // reading of them; then its state.
    applyDrives(this.model, cp.drives);
    cp.solver.redrive();
    cp.solver._restore(cp.snap);
  }

  /**
   * Snapshot of the current state.
   * @returns {Solution}
   */
  solution() {
    return makeSolution(this.solver);
  }
}
