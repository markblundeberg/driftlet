// driftlet: 1D drift–diffusion–reaction solver.

import { normalizeDevice, DeviceError } from './device.js';
import { Solver, SolverError } from './solver.js';
import { makeSolution } from './solution.js';
import { merge } from './merge.js';
import { checkOptions } from './errors.js';

export { DeviceError } from './device.js';
import { normalizeDrives } from './device.js';
export { SolverError } from './solver.js';
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
   * drives alone (V, I, R) is applied in place, cheaply. Otherwise the device is rebuilt, keeping the current state as
   * the warm start when the grid and species are unchanged (else restarting from the regions' c0).
   * @param {object} patch a partial device definition, merged into the current one
   * @returns {this}
   */
  set(patch) {
    // A terminal given a new kind of drive drops the old one: { I: 0 } replaces a held V.
    // (A current drives a terminal alone: given I, a contact drops V and any series R too.)
    const redrive = (c) => (c && typeof c === 'object' && !Array.isArray(c) ? ('I' in c && !('V' in c) ? { ...c, V: undefined, ...('R' in c ? {} : { R: undefined }) } : 'V' in c && !('I' in c) ? { ...c, I: undefined } : c) : c);
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
        return q === undefined ? p : merge(p, redrive(q));
      }) };
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
    // The state carries over where the unknowns line up; on a new grid (the same regions, faces
    // and terminals), interpolated onto it. (Started cold instead, a floating gate lost its
    // charge to a refined grid.)
    const sameUnknowns = old && old.u.length === solver.u.length && old.n === solver.n;
    const regridded = old && !sameUnknowns && old.M === solver.M && old.n === solver.n && old.nFaces === solver.nFaces && old.terms.length === solver.terms.length && old.model.grid.regionStart.length === model.grid.regionStart.length;
    if (sameUnknowns || regridded) {
      if (sameUnknowns) {
        solver.u.set(old.u);
        solver.uLo.set(old.uLo);
        // A node newly on an electrode's surface (a window moved, a surface added) starts at its
        // starting coverages, not at whatever its slots held.
        const surfaceAt = (sv, g) => (sv.surfAt[g] < 0 ? '' : sv.model.surfaces[sv.surfAt[g]].species.map((sp) => sp.name).join());
        for (let g = 0; g < solver.nNodes; g++) if (solver.surfAt[g] >= 0 && surfaceAt(solver, g) !== surfaceAt(old, g)) solver._surfaceStart(g);
      } else solver._warmFrom(old);
      solver.computeConcentrations();
      solver.time = old.time;
      // Held sources read at the same time as before, and an adaptive transient goes on at the
      // step size it had reached (its order restarts, as at any change of the equations).
      solver.sourceTime = old.sourceTime;
      solver._refreshSources();
      solver.dtNext = old.dtNext;
      solver.coldStart = old.coldStart; // (a carried state is as cold as it was: see equilibrateFaces)
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
        const [k, q] = solver.surfCols[j], o = old.surfCols.findIndex(([kk, qq]) => old.model.surfaces[kk].name === model.surfaces[k].name && qq === q);
        const a = old.model.surfaces[old.surfCols[o]?.[0]], b = model.surfaces[k];
        return o >= 0 && a.species[q].name === b.species[q].name && a.nodes.length === b.nodes.length && a.nodes.every((g, w) => g === b.nodes[w]) ? o : -1;
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
    options(opts, 'solve', ['maxSteps', 'tol', 'continuation']);
    // A solve that fails reports the state it failed in, and leaves the device as it was: a
    // device with no steady state (a host filling at a current) can still advance() from there.
    const snap = this.solver._snapshot();
    const result = this.solver.solveSteady(opts);
    // (a failed solve's state isn't a solution: its accuracy warnings are left out)
    const sol = makeSolution(this.solver, result.converged ? result : { ...result, unsolved: true });
    if (!sol.converged) this.solver._restore(snap);
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
    return makeSolution(this.solver, this.solver.advance(dt, opts));
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
    return makeSolution(this.solver, this.solver.integrate(tEnd, opts));
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
    const steady = this.solver.solveSteady();
    if (!steady.converged) throw new SolverError('impedance: the steady state did not converge');
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
          'Measuring at another terminal, or with this one driven by a current (at its steady current, e.g. I: 0), can avoid it: a closed host at open circuit resolved so, blocking electrodes did not.',
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
