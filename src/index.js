// driftlet: 1D drift–diffusion–reaction solver.

import { normalizeDevice } from './device.js';
import { Solver, SolverError } from './solver.js';
import { makeSolution } from './solution.js';

export { DeviceError, normalizeDevice } from './device.js';
export { SolverError } from './solver.js';
export { buildGrid, gradedCells } from './grid.js';
export { BlockTridiagonal } from './blockTridiagonal.js';
export { bernoulli, bernoulliDerivative } from './bernoulli.js';
export * from './constants.js';

const isPlainObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v) && !ArrayBuffer.isView(v);

// Deep-merge plain objects; arrays and other values are replaced wholesale.
function merge(base, patch) {
  if (!isPlainObject(base) || !isPlainObject(patch)) return patch;
  const out = { ...base };
  for (const [k, v] of Object.entries(patch)) out[k] = isPlainObject(v) && isPlainObject(base[k]) ? merge(base[k], v) : v;
  return out;
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
   * Change part of the definition (deep-merged). The current state is kept as the warm start
   * when the grid and species are unchanged; otherwise it restarts from the regions' c0.
   * @param {object} patch a partial device definition, merged into the current one
   * @returns {this}
   */
  set(patch) {
    const def = merge(this.def, patch);
    const model = normalizeDevice(def);
    const old = this._solver;
    this.def = def;
    this.model = model;
    this._solver = null;
    if (old) {
      const solver = this.solver;
      if (old.u.length === solver.u.length && old.n === solver.n) {
        solver.u.set(old.u);
        solver.uLo.set(old.uLo);
        solver.computeConcentrations();
        solver.time = old.time;
        solver.contactDEnd = old.contactDEnd;
        solver.solvedV = old.solvedV; // where the carried-over state was solved (for continuation)
        solver.referenceAmounts = old.referenceAmounts.slice();
        if (old.stretches.length === solver.stretches.length) solver.boundaryIntake.set(old.boundaryIntake);
      }
    }
    return this;
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
   * which), so an animation can call it once per frame; the step size carries over.
   * @param {number} tEnd s
   * @param {{ tol?: number, dt0?: number, dtMax?: number, budgetMs?: number, maxSteps?: number, method?: 'bdf2' | 'be' }} [opts]
   * @returns {Solution}
   */
  advance(tEnd, opts) {
    return makeSolution(this.solver, this.solver.integrate(tEnd, opts));
  }

  /**
   * Small-signal impedance Z(f) about the steady state (solved first). In voltage mode the
   * right terminal's voltage is perturbed; in current mode, the circuit current.
   * Z = −δV/δI, the impedance seen at the terminals (Ω·m²).
   * @param {ArrayLike<number>} frequencies Hz
   * @param {{ profiles?: boolean }} [opts] also return complex profiles per frequency
   * @returns {import('./types.js').ImpedanceResult}
   */
  impedance(frequencies, opts) {
    const steady = this.solver.solveSteady();
    if (!steady.converged) throw new SolverError('impedance: the steady state did not converge');
    return this.solver.impedance(frequencies, opts);
  }

  /**
   * Snapshot of the current state.
   * @returns {Solution}
   */
  solution() {
    return makeSolution(this.solver);
  }
}
