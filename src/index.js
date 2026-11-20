// driftlet: 1D drift–diffusion–reaction solver.

import { normalizeDevice, DeviceError } from './device.js';
import { Solver, SolverError } from './solver.js';
import { makeSolution } from './solution.js';
import { deepEqual, markWhole, merge, ownCopy, WHOLE } from './merge.js';
import { carryState, carryTermCharges } from './carry.js';
import { checkOptions } from './errors.js';

export { DeviceError } from './device.js';
import { nextBreakpoint, normalizeDrives, sourceAt } from './device.js';

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
  if (Array.isArray(x) && x.length !== y.length) return false; // (holes and a trailing undefined count)
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

// A definition's grid: the device's, and each region's own (grid, cells); def with `from`'s, where
// they have as many regions; and def without any. (Read off def.grid alone, a grid given per region
// missed the regrid's settle and the two-step patch: a bias step's first current the wrong sign;
// Astra's second pass.)
const REGION_GRID = ['grid', 'cells'];
const gridOf = (d) => ({ grid: d.grid, regions: (d.regions ?? []).map((r) => REGION_GRID.map((k) => r?.[k])) });
function withGridOf(d, from) {
  const out = { ...d, grid: from.grid };
  if (out.grid === undefined) delete out.grid;
  if (!Array.isArray(d.regions) || !Array.isArray(from.regions) || d.regions.length !== from.regions.length) return null;
  out.regions = d.regions.map((r, k) => {
    const o = { ...r };
    for (const key of REGION_GRID) {
      if (from.regions[k]?.[key] === undefined) delete o[key];
      else o[key] = from.regions[k][key];
    }
    return o;
  });
  return out;
}
const withoutGrid = (d) => ({ ...d, grid: undefined, regions: (d.regions ?? []).map((r) => (r && typeof r === 'object' ? { ...r, grid: undefined, cells: undefined } : r)) });

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
  carryTermCharges(old, s); // (by name, scaled as the carry scales them: one rule; Astra's second pass)
}

// Each terminal's drive (the contacts, then the ports), set in the model in place.
function applyDrives(model, drives) {
  model.terminals.forEach((t, k) => (t.drive = drives[k]));
  model.contacts.left.drive = drives[0];
  model.contacts.right.drive = drives[1];
  model.ports.forEach((p, k) => (p.drive = drives[2 + k]));
}

// Which state of which device a solution is of, for check() and polarization(): every set(),
// and every solve or step that moves the state, gives the device a new mark, and each solution
// carries the mark of the state it was made from (non-enumerable, so a copy carries none and is
// compared by its numbers). Compared by its numbers alone, a solution from before a set() of a
// material passed as the present state: φ is NaN in a metal, which masked every node after it,
// and a current 60× the device's passed check() (Astra's fifth pass). The marks are counted across
// devices, so another device's solution is told apart (a twin's, of the same definition, is
// compared by its numbers, as a copy is).
const MARK = Symbol('driftlet.state');
let marks = 0;

// A method's options: a name it doesn't take (a typo, dtmax for dtMax) is an error, not ignored.
const options = (opts, method, allowed) => checkOptions(method, opts, allowed, SolverError);

/** @typedef {import('./types.js').DeviceDefinition} DeviceDefinition */
/** @typedef {import('./types.js').Solution} Solution */

export class Device {
  /** @param {DeviceDefinition} def plain, serialisable device definition (see docs/device.md) */
  constructor(def) {
    // (its own, copied before it's read: normalized from the caller's, a statistics model kept
    // the caller's arrays, and an edit to one after construction made a hybrid of two models;
    // Astra's third pass)
    this.def = ownCopy(def); // (a caller's later edits to def reach it only through set())
    this.model = normalizeDevice(this.def);
    this._solver = null;
    this._touched = false; // (solved or stepped: a state of its own to carry through set())
    // (def: the definition's mark, which a set() that rebuilds changes, and source the definition
    // itself, for another device's solution; state: the state's, which any set() and every solve
    // or step changes; solved: whether that state is a solution, not a layout or a state set() has
    // since changed the equations of; by: what last moved it)
    this._mark = { device: ++marks, def: ++marks, source: this.def, state: ++marks, solved: false, by: 'built' };
  }

  // The state moved (by a set(), a solve or a step): a new mark. (A rebuild's new model is a new
  // definition too.)
  _moved(by, solved, rebuilt = false) {
    this._mark = { device: this._mark.device, def: rebuilt ? ++marks : this._mark.def, source: this.def, state: ++marks, solved, by };
  }

  // A solution marked as of the present state (of: 0 for a state the device didn't keep).
  _marked(sol, of = this._mark) {
    Object.defineProperty(sol, MARK, { value: of, configurable: true });
    return sol;
  }

  /**
   * Which state of this device a solution is of, by its mark: undefined for a solution without
   * one (a copy, `{ ...sol }`, or one sent from a worker) or a twin's (another device's, of the
   * same definition), for the caller to compare by its numbers; otherwise `moved`, why it isn't
   * the present state (null if it is), `unsolved`, why the present state isn't a solution (null
   * if it is), and `definition`, whether it is of the present definition (a later solve or step
   * leaves that, a set() that rebuilds doesn't).
   * @param {object} sol
   * @returns {{ moved: string | null, unsolved: string | null, definition: boolean } | undefined}
   * @internal
   */
  _stateOf(sol) {
    const m = sol?.[MARK], now = this._mark;
    if (m === undefined) return undefined;
    const unsolved = now.solved ? null : now.by === 'built' ? 'the device holds its starting layout, never solved or advanced' : `the device's state hasn't been solved since a ${now.by} changed its equations`;
    if (m === 0) return { moved: "it's the state a solve gave up in, which the device didn't keep", unsolved, definition: false };
    // (another device's, of the same definition: a twin's state, compared by its numbers)
    if (m.device !== now.device) return deepEqual(m.source, this.def) ? undefined : { moved: "it is another device's, of another definition: check a solution with the device that returned it", unsolved, definition: false };
    const definition = m.def === now.def;
    const moved = m.state === now.state ? null : definition ? `the device has moved on since it was made, most recently by a ${now.by}` : 'the device has been set() to a new definition since it was made';
    return { moved, unsolved, definition };
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
    const def = ownCopy(merge(this.def, patch)); // (its own: see ownCopy)
    // A drive is replaced whole, not merged: a waveform's keys are its own (merged, a ramp set
    // over a triangle kept the triangle's repeat; Fable's review of the front door).
    // Likewise a link or a φ law given with its type: a velocity link set over an equilibrium one
    // kept its offset, a pinned φ over a capacitive one its C (Fable's re-check). (A contact the kit wrote replaces the old one whole: merge's WHOLE.)
    for (const side of ['left', 'right']) {
      const c = patch?.contacts?.[side];
      // (copies: the patch's own objects put in, a caller's later edit to one made the next set()
      // of the same value a no-op and leaked in at an unrelated rebuild; Astra's third pass)
      if (c && typeof c === 'object') for (const key of ['V', 'I']) if (c[key] !== undefined) def.contacts[side][key] = ownCopy(c[key]);
      for (const [name, link] of Object.entries(c?.species ?? {})) if (link && typeof link === 'object' && 'type' in link) def.contacts[side].species[name] = ownCopy(link);
      if (c?.phi && typeof c.phi === 'object' && 'type' in c.phi) def.contacts[side].phi = ownCopy(c.phi);
    }
    // Nothing changed (a page resending its drives every frame): nothing to do. (Redriven, a
    // transient's BDF2 history was cleared and a steady state forgot it was one; Astra's review.)
    if (deepEqual(def, this.def)) return this;
    // A regridded steady state not yet solved on its grid stands for the steady state at the old
    // definition: solve it there first (see Device.advance); once the patch is known to be valid
    // (refused, it settled the earlier change and left the device not as it was; Astra's second pass).
    if (this._solver?.regriddedSteady) normalizeDevice(def);
    this._settleRegrid();
    // A solver nothing has solved or stepped holds only its start, which the new definition's
    // replaces. (Kept, a drive patch after a bare solution() kept the old drive's cold start, and
    // the first step read −73846 A/m² where 26154 is right; Astra's review.)
    const live = this._solver && (this._touched || this._solver.restarted) ? this._solver : null;
    // Only the terminals' drives changed: update them in place, keeping the solver and its state
    // (a step in a source restarts the time stepping's order, as at any discontinuity).
    if (live && sameExceptDrives(this.def, def) && plainContacts(def)) {
      this._redefineDrives(def);
      this._moved('set()', false);
      return this;
    }
    // A rebuild, carrying the state by what each value means (src/carry.js). It changes the device
    // only once it has all succeeded: a change that fails leaves it as it was.
    const was = { def: this.def, model: this.model, solver: this._solver, mark: this._mark };
    try {
      // A new grid and anything else, on a steady state: the new grid first, the steady state
      // solved again there (the state it stands for), then the rest, as two set()s are. (In one,
      // the regrid was settled at the new definition, and the change's transient lost: a bias
      // stepped, −859 A/m² where −2.0e7 is right, Sol's review; a contact's level, a midpoint at 5
      // where 1 is right, Astra's.)
      const mid = withGridOf(this.def, def);
      if (live?.atSteady && mid && !deepEqual(gridOf(this.def), gridOf(def)) && !deepEqual(mid, def)) {
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
          this._moved('set()', false, true);
          return this;
        }
      }
      this._rebuild(def, live);
    } catch (err) {
      this.def = was.def;
      this.model = was.model;
      this._solver = was.solver;
      this._mark = was.mark;
      throw err;
    }
    this._moved('set()', false, true);
    return this;
  }

  // Rebuild at def, carrying old's state (src/carry.js); a device with no state of its own just
  // takes the new definition.
  _rebuild(def, old) {
    const model = normalizeDevice(def);
    const solver = this._solver ? new Solver(model) : null;
    if (old) {
      // Its own start at the sources then, not at t = 0's: what nothing carries starts from there
      // (a floating port added mid-ramp started at the ramp's first voltage, and kept it; Astra's
      // second pass).
      solver.time = old.time;
      solver.sourceTime = old.sourceTime;
      solver._refreshSources();
      solver.initFromComposition();
      // The same species (by name and charge, in any order) carry over; other species restart
      // the device from its regions' c0.
      const key = (m) => JSON.stringify(m.species.map((sp) => [sp.name, sp.z]).sort());
      if (key(old.model) === key(model)) carryState(old, solver, { oldDef: this.def, newDef: def, pure: deepEqual(withoutGrid(this.def), withoutGrid(def)) });
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
    // (a held terminal set floating starts at the voltage the state was solved at, the source as
    // the last step saw it: at a waveform's landing, read after the jump, it started 1 V off;
    // Astra's second pass)
    const s = this._solver, seen = s ? s.terms.map((t) => (t.drive.kind === 'V' && !(t.drive.R > 0) ? sourceAt(t.drive.src, s.time, s.stepped) : NaN)) : [];
    // A regridded steady state stands for the old drives' (see Device.advance): solve it now, before
    // they change. (Left for advance(), it was solved at the new drive, and the step response
    // was lost; Fable round 3's steady.js review.)
    this._settleRegrid();
    // (and is read at the voltages it was solved at, until a step or a solve moves it: read at
    // the new ones, a patch from 0 to 0.05 V read a plate charge of 1e-2 C/m² and 1e4 A/m² that
    // never flowed; Sol's third review. A step from it is the step response.) A floating terminal
    // likewise passes what its old drive said: read at the new one, a capacitive plate set from
    // I: 0 to 2 read 2 A/m² flowing at no step, its terminals no longer summing to zero, without a
    // word (Sol's fourth pass).
    if (!s.readV) s.readI = Float64Array.from(s.terms, (t, k) => (s.floating.includes(k) ? s.passedCurrent(k) : NaN));
    s.readV = s.readV ?? Float64Array.from(s.termV, (v, k) => (Number.isFinite(seen[k]) ? seen[k] : v)); // (a held one as the last step saw its source)
    applyDrives(this.model, drives);
    this.def = def;
    this._solver.redrive();
    s.terms.forEach((t, k) => {
      if (Number.isFinite(seen[k]) && s.floating.includes(k)) s.termV[k] = seen[k];
    });
  }

  /**
   * Steady state (or equilibrium), warm-started from the current state. Doesn't advance time.
   * @param {{ maxSteps?: number, tol?: number, continuation?: boolean }} [opts]
   * @returns {Solution}
   */
  solve(opts) {
    options(opts, 'solve', ['maxSteps', 'tol', 'continuation']);
    // (as advance() checks its own: tol NaN with maxSteps Infinity stepped without end, a NaN
    // motion never under tol; Astra's fourth pass)
    const { tol, maxSteps } = opts ?? {};
    if (tol !== undefined && !(tol > 0 && tol < Infinity)) throw new SolverError(`solve(): tol is the largest potential change, in thermal units, that ends a ramp of giant steps, a positive number; got ${tol}`);
    if (maxSteps !== undefined && !(Number.isInteger(maxSteps) && maxSteps >= 1)) throw new SolverError(`solve(): maxSteps is the time steps per pseudo-transient ramp, a whole number ≥ 1; got ${maxSteps}`);
    // A solve that fails reports the state it failed in, and leaves the device as it was: a
    // device with no steady state (a host filling at a current) can still advance() from there.
    // (And one that throws, as _touching does for the transients: a custom statistics that threw
    // mid-solve left the state moved and the clock at 44000 s, and solution() threw after it;
    // Astra's fourth pass.)
    const snap = this.solver._snapshot();
    let result, sol;
    try {
      result = this.solver.solveSteady(opts);
      if (result.converged) this.solver.readV = undefined; // (solved at its drives, and read so: see _redefineDrives)
      // (a failed solve's state isn't a solution: its accuracy warnings are left out)
      sol = makeSolution(this.solver, result.converged ? result : { ...result, unsolved: true });
    } catch (e) {
      this.solver._restore(snap);
      throw e;
    }
    if (!sol.converged) {
      this.solver._restore(snap);
      return this._marked(sol, 0);
    } else {
      this._moved('solve()', true);
      this._touched = true;
      this.solver.regriddedSteady = false; // (solved on the new grid already: see Device.advance)
      this.solver.jumpedT = undefined; // (solved at the new temperature: see integrate in src/timeSteps.js)
    }
    return this._marked(sol);
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
    if (!(dt > 0)) return this._marked(makeSolution(this.solver, { converged: true, iterations: 0, substeps: 0, stopped: `step(${dt}) did nothing: ${length}` }));
    // A step too short to move the clock (1e-17 s at t = 1 s) can't be taken: it posted charge
    // while the time stood still (Sol's third review).
    if (this.solver.time + dt === this.solver.time) return this._marked(makeSolution(this.solver, { converged: true, iterations: 0, substeps: 0, stopped: `step(${dt}) did nothing: it's below the time's resolution at t = ${this.solver.time} s` }));
    this._settleRegrid(); // (as advance() does: stepped unsolved, 27 A/m² where 68.8 is right)
    this.solver._freshDiagnostics();
    const crossed = this._crossed(dt);
    const sol = this._touching(() => makeSolution(this.solver, this.solver.advance(dt, opts)), 'step()');
    if (crossed) sol.warnings.push(crossed);
    return sol;
  }

  // A drive's breakpoint strictly inside (t, t + dt), which a step crosses unlanded, reading the
  // drive at its end alone: what it says, or null. (One step of 3 ms over a 1 ms pulse of current
  // passed no charge at all, silently; Astra's fifth pass.)
  _crossed(dt) {
    const s = this.solver, t0 = s.time, t1 = t0 + dt;
    let first = null;
    for (const term of s.terms) {
      let tb;
      try {
        tb = nextBreakpoint(term.drive.src, t0);
      } catch {
        continue; // (a period below the time's resolution: the step reads it as it always has)
      }
      if (tb < t1 - 1e-12 * Math.abs(t1) && (first === null || tb < first.tb)) first = { tb, name: term.name };
    }
    if (first === null) return null;
    return `step(${dt}): the drive of terminal '${first.name}' has a breakpoint at t = ${first.tb} s, inside this step, which reads the drive at its end alone and so crosses it unlanded (a pulse inside one step is missed whole, its charge lost): advance(tEnd) lands on every breakpoint, or step to it first`;
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
    // (an advance to now takes no step: nothing to settle for; nor does a tEnd that isn't a number,
    // which integrate refuses: checked first, advance('1') settled the regrid before it threw;
    // Sol's third review)
    if (typeof tEnd === 'number' && tEnd > this.solver.time) this._settleRegrid();
    this.solver._freshDiagnostics();
    return this._touching(() => makeSolution(this.solver, this.solver.integrate(tEnd, opts)), 'advance()');
  }

  // A call that moved the clock has a state of its own to carry through set(). (Marked before it
  // ran, a step(0), a refused step or a failed first solve made c0 a "state" that a c0 patch
  // then couldn't replace; Astra's review.)
  // A call that throws (a statistics model of the user's that throws an Error, not a refusal)
  // leaves the state as it found it, then throws: left at the refused trial, a later call
  // converged 0.7% off in its amounts, or the device couldn't be read or set() at all (Sol's
  // third review).
  // (Its steps, accepted, are solutions of the device as it stands: its state is solved, and its
  // solution marked so. A call that took none left the state as it was, and its solution is of it.)
  _touching(call, by) {
    const t0 = this.solver.time, snap = this.solver._snapshot();
    let out;
    try {
      out = call();
    } catch (e) {
      this.solver._restore(snap);
      this.solver.sourceBefore = false;
      this.solver.landing = undefined;
      this.solver.combining = false;
      throw e;
    }
    if (this.solver.time !== t0) {
      this._touched = true;
      this._moved(by, true);
    }
    return this._marked(out);
  }

  // A steady state carried onto a new grid, solved there (once; see Device.advance).
  // (A settle that fails leaves the interpolated state, and says so in the next solutions: it
  // was dropped silently, the flag cleared before the snapshot that restored it; Astra's review.)
  _settleRegrid() {
    const s = this._solver;
    if (!s?.regriddedSteady) return;
    const snap = s._snapshot();
    if (s.solveSteady().converged) {
      s.regriddedSteady = false;
      this._moved('solve()', true); // (the steady state there: see Device.advance)
    } else {
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
    // steady.js review.) What the arguments alone refuse is refused first, and any later refusal
    // (a frequency where the linearised system is singular) leaves it as it was too: a refused
    // call that had solved the steady state took a pending step response with it, 3420 A/m² of
    // a plate set to 0.1 V read 0 at the next step (Sol's fourth pass).
    this.solver.impedanceArgs(frequencies, opts?.terminal);
    const snap = this.solver._snapshot();
    const steady = this.solver.solveSteady();
    if (!steady.converged) {
      // (with why, as solve() would say it: a terminal driven past what any steady state passes,
      // say, where advance() is the way; Fable round 2's cold user, charging a battery)
      const why = makeSolution(this.solver, { ...steady, unsolved: true }).warnings;
      this.solver._restore(snap);
      throw new SolverError(`impedance: the steady state did not converge${why.length > 0 ? `: ${why[0].replace(/\.$/, '')}` : ''}. An impedance is about a steady state; about a transient one, drive a small sine with advance() and fit it`, { warnings: why });
    }
    this.solver.regriddedSteady = false; // (solved on its grid: nothing left to settle; Astra's second pass)
    // The steady state's warnings (an unresolved double layer, say) hold for the impedance too.
    const { warnings } = makeSolution(this.solver, steady);
    let out, unconverged;
    try {
      ({ unconverged, ...out } = this.solver.impedance(frequencies, opts));
    } catch (err) {
      this.solver._restore(snap);
      throw err;
    }
    this._touched = true; // (solved: a state of its own)
    this._moved('impedance()', true); // (its steady state solved again)
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
    if (unconverged.length > 0) {
      warnings.push(`impedance: the iterative solve stopped short of its tolerance at ${g(unconverged[0])}${unconverged.length > 1 ? `…${g(unconverged.at(-1))}` : ''} Hz (${unconverged.length} of ${out.f.length} frequencies), by more than relError there allows for: Z there may be less accurate than relError says.`);
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
    return { def: this.def, model: this.model, solver, touched: this._touched, mark: this._mark, drives: this.model.terminals.map((t) => t.drive), snap: solver._snapshot() };
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
    this._mark = cp.mark; // (the state as it was, and the solutions of it are of the present state again)
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
    return this._marked(makeSolution(this.solver));
  }
}
