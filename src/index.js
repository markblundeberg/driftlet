// driftlet: 1D drift–diffusion–reaction solver.

import { normalizeDevice } from './device.js';
import { Solver, SolverError } from './solver.js';
import { makeSolution } from './solution.js';
import { merge } from './merge.js';

export { DeviceError, normalizeDevice } from './device.js';
import { normalizeDrives } from './device.js';
export { SolverError } from './solver.js';
export { buildGrid, gradedCells } from './grid.js';
export { BlockTridiagonal } from './blockTridiagonal.js';
export { bernoulli, bernoulliDerivative } from './bernoulli.js';
export * from './constants.js';

// Whether two definitions differ only in their terminals' drives (V, I, R on the contacts and
// ports). Functions (custom statistics) compare by identity.
function sameExceptDrives(a, b) {
  const drive = new Set(['V', 'I', 'R']);
  const eq = (x, y, terminal) => {
    if (x === y) return true;
    if (typeof x !== 'object' || typeof y !== 'object' || x === null || y === null) return false;
    if (Array.isArray(x) !== Array.isArray(y)) return false;
    const keys = new Set([...Object.keys(x), ...Object.keys(y)]);
    for (const k of keys) {
      if (terminal && drive.has(k)) continue;
      if (!eq(x[k], y[k], false)) return false;
    }
    return true;
  };
  const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
  for (const k of keys) {
    if (k === 'contacts') {
      const ca = a.contacts ?? {}, cb = b.contacts ?? {};
      for (const side of new Set([...Object.keys(ca), ...Object.keys(cb)])) if (!eq(ca[side] ?? {}, cb[side] ?? {}, true)) return false;
    } else if (k === 'ports') {
      const pa = a.ports ?? [], pb = b.ports ?? [];
      if (pa.length !== pb.length || pa.some((p, i) => !eq(p, pb[i], true))) return false;
    } else if (!eq(a[k], b[k], false)) return false;
  }
  return true;
}

/** @typedef {import('./types.js').DeviceDefinition} DeviceDefinition */
/** @typedef {import('./types.js').Solution} Solution */

export class Device {
  /** @param {DeviceDefinition} def plain, serialisable device definition (see docs/device.md) */
  constructor(def) {
    this.def = def;
    this.model = normalizeDevice(def);
    this._solver = null;
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
   * Change part of the definition (deep-merged). A contact given V or I drops the other, so
   * `{ contacts: { right: { I: 0 } } }` switches it to open circuit. A change to the terminals'
   * drives alone (V, I, R) is applied in place, cheaply. Otherwise the device is rebuilt, keeping the current state as
   * the warm start when the grid and species are unchanged (else restarting from the regions' c0).
   * @param {object} patch a partial device definition, merged into the current one
   * @returns {this}
   */
  set(patch) {
    // A contact given a new kind of drive drops the old one: { I: 0 } replaces a held V.
    for (const side of ['left', 'right']) {
      const c = patch?.contacts?.[side];
      if (c && typeof c === 'object') {
        if ('I' in c && !('V' in c)) patch = { ...patch, contacts: { ...patch.contacts, [side]: { ...c, V: undefined } } };
        else if ('V' in c && !('I' in c)) patch = { ...patch, contacts: { ...patch.contacts, [side]: { ...c, I: undefined } } };
      }
    }
    const def = merge(this.def, patch);
    // Only the terminals' drives changed: update them in place, keeping the solver and its state
    // (a step in a source restarts the time stepping's order, as at any discontinuity).
    if (this._solver && sameExceptDrives(this.def, def)) {
      this._redefineDrives(def);
      return this;
    }
    // A rebuild. Everything that can fail (validation, building the solver) happens before the
    // device changes, so a change that fails leaves it as it was.
    const model = normalizeDevice(def);
    const old = this._solver;
    const solver = old ? new Solver(model) : null;
    this.def = def;
    this.model = model;
    this._solver = solver;
    if (old && old.u.length === solver.u.length && old.n === solver.n) {
      solver.u.set(old.u);
      solver.uLo.set(old.uLo);
      solver.computeConcentrations();
      solver.time = old.time;
      solver.contactDEnd = old.contactDEnd;
      if (old.portQEnd?.length === solver.portQ.length) solver.portQEnd = old.portQEnd;
      solver.solvedV = old.solvedV; // where the carried-over state was solved (for continuation)
      // Floating terminals keep their voltages (a good start), where the terminals match.
      if (old.terms.length === solver.terms.length) {
        for (const k of solver.floating) if (old.terms[k].name === solver.terms[k].name) solver.termV[k] = old.termV[k];
      }
      // What each stretch conserves from here: a stretch that was closed and still is (same
      // species, same regions, same boxes) keeps its amount; any other (newly closed, split,
      // merged, open, or with its boxes resized by a new grid or geometry) starts from what the
      // carried state holds.
      const key = (st) => `${st.species}:${st.regions[0]}-${st.regions[1]}`;
      const before = new Map(old.stretches.map((st, k) => [key(st), k]));
      const va = old.model.grid.vol, vb = model.grid.vol, sameBoxes = va.length === vb.length && va.every((v, g) => v === vb[g]);
      solver.referenceAmounts = solver.stretches.map((st) => {
        const k = before.get(key(st));
        return sameBoxes && k !== undefined && !old.stretches[k].connected && !st.connected ? old.referenceAmounts[k] + old.boundaryIntake[k] : solver.amount(st);
      });
      solver.boundaryIntake.fill(0);
      // Each surface's reference likewise: kept where its window and the boxes are unchanged.
      const sameSurface = (j) => {
        const [k, q] = solver.surfCols[j], o = old.surfCols.findIndex(([kk, qq]) => old.model.ports[kk].name === model.ports[k].name && qq === q);
        const a = old.model.ports[old.surfCols[o]?.[0]], b = model.ports[k];
        return o >= 0 && a.surface[q].name === b.surface[q].name && a.nodes.length === b.nodes.length && a.nodes.every((g, w) => g === b.nodes[w]) ? o : -1;
      };
      solver.surfaceRef = Float64Array.from(solver.surfCols, (_, j) => {
        const o = sameBoxes ? sameSurface(j) : -1;
        return o >= 0 ? old.surfaceRef[o] : solver.surfaceAmount(j);
      });
    }
    return this;
  }

  // The definition def, which differs from the current one only in the terminals' drives: set
  // them in place. (Validated before anything changes.)
  _redefineDrives(def) {
    const drives = normalizeDrives(def, this.model);
    this.model.terminals.forEach((t, k) => (t.drive = drives[k]));
    this.model.contacts.left.drive = drives[0];
    this.model.contacts.right.drive = drives[1];
    this.model.ports.forEach((p, k) => (p.drive = drives[2 + k]));
    this.def = def;
    this._solver.redrive();
  }

  /**
   * Steady state (or equilibrium), warm-started from the current state. Doesn't advance time.
   * @param {{ maxSteps?: number, tol?: number, continuation?: boolean }} [opts]
   * @returns {Solution}
   */
  solve(opts) {
    return makeSolution(this.solver, this.solver.solveSteady(opts));
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
    return makeSolution(this.solver, this.solver.advance(dt, opts));
  }

  /**
   * Adaptive transient to time tEnd (s): variable-step BDF2 with local error control. With
   * `budgetMs`, it returns after that much wall time even if tEnd isn't reached (`done` says
   * which), so an animation can call it once per frame; the step size carries over. `probes`
   * read a species at points inside the device after every accepted step, into `trace.probes`
   * (one array per probe, alongside `trace.t`): what a detector at x sees.
   * @param {number} tEnd s
   * @param {{ tol?: number, dt0?: number, dtMax?: number, budgetMs?: number, maxSteps?: number, method?: 'bdf2' | 'be',
   *   probes?: { x: number, species?: string, quantity?: 'c' | 'V' | 'phi', region?: string | number }[] }} [opts]
   * @returns {Solution}
   */
  advance(tEnd, opts) {
    return makeSolution(this.solver, this.solver.integrate(tEnd, opts));
  }

  /**
   * Small-signal impedance Z(f) about the steady state (solved first), at one terminal: its
   * held voltage is perturbed (or its driven current), the others keep their drives.
   * Z = δV/δI with I into the device (Ω·m²).
   * @param {ArrayLike<number>} frequencies Hz
   * @param {{ terminal?: string, profiles?: boolean }} [opts] the terminal (default 'right'), and
   *   whether to return complex profiles per frequency
   * @returns {import('./types.js').ImpedanceResult}
   */
  impedance(frequencies, opts) {
    const steady = this.solver.solveSteady();
    if (!steady.converged) throw new SolverError('impedance: the steady state did not converge');
    return this.solver.impedance(frequencies, opts);
  }

  /**
   * A checkpoint of the definition and solver state, for `_rollback()` (the live wrapper keeps
   * its last good state this way).
   * @internal
   */
  _checkpoint() {
    const solver = this.solver;
    return { def: this.def, model: this.model, solver, snap: solver._snapshot(), solvedV: solver.solvedV, referenceAmounts: solver.referenceAmounts.slice(), surfaceRef: solver.surfaceRef.slice() };
  }

  /**
   * Back to a checkpoint taken on this device.
   * @param {ReturnType<Device['_checkpoint']>} cp
   * @internal
   */
  _rollback(cp) {
    if (this._solver === cp.solver) this._redefineDrives(cp.def); // only drives changed since: restore them exactly
    else {
      this.def = cp.def;
      this.model = cp.model;
      this._solver = cp.solver;
    }
    cp.solver._restore(cp.snap);
    cp.solver.solvedV = cp.solvedV;
    cp.solver.referenceAmounts = cp.referenceAmounts.slice();
    cp.solver.surfaceRef = cp.surfaceRef.slice();
  }

  /**
   * Snapshot of the current state.
   * @returns {Solution}
   */
  solution() {
    return makeSolution(this.solver);
  }
}
