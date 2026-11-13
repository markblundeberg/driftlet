// driftlet: 1D drift–diffusion–reaction solver.

import { normalizeDevice, DeviceError } from './device.js';
import { Solver, SolverError } from './solver.js';
import { makeSolution } from './solution.js';
import { markWhole, merge, WHOLE } from './merge.js';
import { checkOptions } from './errors.js';

export { DeviceError } from './device.js';
import { normalizeDrives } from './device.js';

// (a contact the kit wrote stays whole through set()'s own copies of it: see merge.js)
const keepWhole = (from, to) => (from[WHOLE] ? markWhole(to) : to);
export { SolverError } from './errors.js';
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
   * drives alone (V, I, R) is applied in place, cheaply. Otherwise the device is rebuilt, keeping the current state as
   * the warm start when the grid and species are unchanged (else restarting from the regions' c0).
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
    // Only the terminals' drives changed: update them in place, keeping the solver and its state
    // (a step in a source restarts the time stepping's order, as at any discontinuity).
    if (this._solver && sameExceptDrives(this.def, def)) {
      this._redefineDrives(def);
      return this;
    }
    // A rebuild. Everything that can fail (validation, building the solver) happens before the
    // device changes, so a change that fails leaves it as it was.
    const model = normalizeDevice(def);
    const solver = this._solver ? new Solver(model) : null;
    // A regridded steady state not yet solved on its grid stands for the steady state at the old
    // definition: solve it there first (see advance), or the flag that says so is lost with the
    // old solver, and the next transient stepped from the interpolated state (a second set()
    // after set({ grid }): 23.7 A/m² where 68.8 is right; Fable's re-check).
    this._settleRegrid();
    // A solver nothing has solved or stepped holds only its c0, which the new definition's c0
    // replaces (a c0 patch after a bare solution() was otherwise ignored).
    const old = this._solver && (this._touched || this._solver.restarted) ? this._solver : null;
    this.def = def;
    this.model = model;
    this._solver = solver;
    // The state carries over where the unknowns line up; on a new grid (the same regions, faces
    // and terminals), interpolated onto it. What happens to each STATE field is named in CARRY
    // (src/solver.js), and test/state.test.js holds this code to it: a field added to STATE needs
    // a line there, and a 'kept' one a line here. (Started cold instead, a floating gate lost its
    // charge to a refined grid.)
    // (The same species in the same order, not just as many: tested by count, a reordered list
    // swapped two species' amounts; Fable's review.)
    const sameSpecies = old && old.n === solver.n && old.model.species.every((sp, i) => sp.name === model.species[i].name && sp.z === model.species[i].z);
    // (and regions not reordered: the same names at other places swapped two closed
    // compartments' contents; a renamed region carries by place)
    const names = (m) => m.regions.map((rg) => rg.name), was = old ? names(old.model) : [], now = names(model);
    const reordered = old && was.length === now.length && was.some((nm, k) => nm !== now[k]) && [...was].sort().join('\u0000') === [...now].sort().join('\u0000');
    const sameUnknowns = sameSpecies && !reordered && old.u.length === solver.u.length;
    const regridded = old && !sameUnknowns && sameSpecies && !reordered && old.M === solver.M && old.nFaces === solver.nFaces && old.terms.length === solver.terms.length && old.model.grid.regionStart.length === model.grid.regionStart.length;
    // A state that can't carry over restarts from c0. Say so in the next solution, where there
    // was a state to lose: a pump taken off a membrane's reactions restarted a closed cell from
    // its c0, converged and unwarned (Fable round 2's cold user).
    if (old && !sameUnknowns && !regridded) {
      const why = !sameSpecies ? 'the species changed'
        : reordered ? 'the regions were reordered'
        : old.model.grid.regionStart.length !== model.grid.regionStart.length || old.nFaces !== solver.nFaces ? 'the regions or faces changed'
        : old.terms.length !== solver.terms.length ? 'the terminals changed'
        : "a face's or port's unknowns changed: its reactions, gates, surface or diffuse layers";
      solver.restarted = `set(): the state couldn't carry over (${why}), so the device restarted from its regions' c0 at t = ${old.time} (the clock carries on): what follows is from there, not from the state before set(). ` +
        'To switch a reaction off and keep the state, keep it listed and make it negligible (k0 or vmax tiny). (Said in every solution up to the next steady solve, which keeps a closed region\'s restarted amounts.)';
      // The clock carries on, so a run's next advance(t) is still at t (it read 0, and an
      // advance(20.001) after a change at t = 20 ran 20 s from c0).
      solver.time = old.time;
      solver.sourceTime = old.sourceTime;
      solver._refreshSources();
    } else if (old?.restarted) solver.restarted = old.restarted; // (said until a steady solve, through further set()s)
    // Each terminal's running charge carries on with the clock, by name (zeroed by a rebuild, a
    // charge read across a set() that changed anything but a drive missed all that had passed:
    // Fable's review; termCharge is in STATE, which this hand-written carry predates).
    if (old) for (let k = 0; k < solver.terms.length; k++) {
      const o = old.terms.findIndex((t) => t.name === solver.terms[k].name);
      if (o >= 0) solver.termCharge[k] = old.termCharge[o];
    }
    if (sameUnknowns || regridded) {
      if (sameUnknowns) {
        const fresh = Float64Array.from(solver.u); // (its own start, from the new c0)
        solver.u.set(old.u);
        solver.uLo.set(old.uLo);
        // (a new T, or materials' μ° and c_ref: carry the concentrations and φ, not their scaled
        // potentials)
        if (old.model.RT !== model.RT || old.mu0hat.some((m, k) => m !== solver.mu0hat[k]) || old.cRef.some((c, k) => c !== solver.cRef[k])) solver._carryAcross(old);
        // (a jump in T, for a transient's first step to name if it fails: see advance)
        if (old.model.T !== model.T || old.jumpedT) solver.jumpedT = [old.jumpedT?.[0] ?? old.model.T, model.T];
        if (solver.jumpedT?.[0] === solver.jumpedT?.[1]) solver.jumpedT = undefined; // (back where it was solved: no jump to blame)
        // A species newly present at a node (a region's material changed to one that holds it)
        // starts from the new c0: its slot held nothing the old solver solved for, and a closed
        // stretch then kept that stale amount for good (X at 133 mol/m³ where c0 said 7; Fable
        // round 5's review of src/solver.js).
        for (let g = 0; g < solver.nNodes; g++) {
          for (let i = 0; i < solver.n; i++) {
            if (!solver.present[g * solver.n + i] || old.present[g * old.n + i]) continue;
            const o = solver.blockOfNode[g] * solver.M + 1 + i;
            solver.u[o] = fresh[o];
            solver.uLo[o] = 0;
          }
        }
        // A node newly on an electrode's surface (a window moved, a surface added) starts at its
        // starting coverages, not at whatever its slots held.
        const surfaceAt = (sv, g) => (sv.surfAt[g] < 0 ? '' : sv.model.surfaces[sv.surfAt[g]].species.map((sp) => sp.name).join());
        for (let g = 0; g < solver.nNodes; g++) if (solver.surfAt[g] >= 0 && surfaceAt(solver, g) !== surfaceAt(old, g)) solver._surfaceStart(g);
      } else {
        solver._warmFrom(old);
        solver.regriddedSteady = old.atSteady; // (see advance)
      }
      solver.computeConcentrations();
      solver.time = old.time;
      // Held sources read at the same time as before, and an adaptive transient goes on at the
      // step size it had reached (its order restarts, as at any change of the equations).
      solver.sourceTime = old.sourceTime;
      solver._refreshSources();
      solver.dtNext = old.dtNext;
      solver.coldStart = old.coldStart; // (a carried state is as cold as it was: see equilibrateFaces)
      solver.stepped = old.stepped; // (and as stepped: a transient's state retries from the layout; see src/steady.js)
      solver.contactDEnd = old.contactDEnd;
      // (each port's by name: swapped by index, two capacitive ports read each other's charge)
      const portOf = model.ports.map((p) => old.model.ports.findIndex((q) => q.name === p.name));
      if (old.portQEnd && portOf.every((o) => o >= 0)) solver.portQEnd = Float64Array.from(portOf, (o) => old.portQEnd[o]);
      solver.solvedV = old.solvedV; // where the carried-over state was solved (for continuation)
      // Floating terminals keep their voltages (a good start), by name.
      for (const k of solver.floating) {
        const o = old.terms.findIndex((t) => t.name === solver.terms[k].name);
        if (o >= 0) solver.termV[k] = old.termV[o];
      }
      // What each stretch conserves from here: a stretch that was closed and still is (same
      // species, same regions, same boxes) keeps its amount; any other (newly closed, split,
      // merged, open, or with its boxes resized by a new grid or geometry) starts from what the
      // carried state holds.
      const key = (st) => `${st.species}:${st.regions[0]}-${st.regions[1]}`;
      const before = new Map(old.stretches.map((st, k) => [key(st), k]));
      const va = old.model.grid.vol, vb = model.grid.vol, sameBoxes = va.length === vb.length && va.every((v, g) => v === vb[g]);
      // (a new grid over the same regions and geometry holds the same amounts: read off the
      // interpolated state instead, a closed one ratcheted down at each regrid, −0.47% over 60;
      // Fable's re-check)
      const sameSpace = sameBoxes || (old.model.regions.every((rg, k) => rg.length === model.regions[k].length) && JSON.stringify(old.model.geometry) === JSON.stringify(model.geometry));
      solver.referenceAmounts = solver.stretches.map((st) => {
        const k = before.get(key(st));
        return sameSpace && k !== undefined && !old.stretches[k].connected && !st.connected ? old.referenceAmounts[k] + old.boundaryIntake[k] : solver.amount(st);
      });
      solver.boundaryIntake.fill(0);
      // ...and on a new grid, a closed stretch's state holds exactly that amount: the
      // interpolation keeps it only to its own error, which a transient would keep for good
      // (−0.02% at a refinement 4× mid-run; a steady solve restores it).
      if (regridded && sameSpace) solver.stretches.forEach((st, j) => {
        const k = before.get(key(st));
        if (k === undefined || old.stretches[k].connected || st.connected) return;
        const have = solver.amount(st), want = solver.referenceAmounts[j], shift = Math.log(want / have);
        if (!Number.isFinite(shift) || shift === 0) return;
        for (let g = st.nodes[0]; g <= st.nodes[1]; g++) if (solver.nodeConductor[g] < 0) solver.u[solver.blockOfNode[g] * solver.M + 1 + st.species] += shift;
        solver.computeConcentrations();
      });
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
    // A regridded steady state stands for the old drives' (see advance): solve it now, before
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
    this._touched = true;
    options(opts, 'solve', ['maxSteps', 'tol', 'continuation']);
    // A solve that fails reports the state it failed in, and leaves the device as it was: a
    // device with no steady state (a host filling at a current) can still advance() from there.
    const snap = this.solver._snapshot();
    const result = this.solver.solveSteady(opts);
    // (a failed solve's state isn't a solution: its accuracy warnings are left out)
    const sol = makeSolution(this.solver, result.converged ? result : { ...result, unsolved: true });
    if (!sol.converged) this.solver._restore(snap);
    else {
      this.solver.regriddedSteady = false; // (solved on the new grid already: see advance)
      this.solver.jumpedT = undefined; // (solved at the new temperature: see advance)
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
    this._touched = true;
    options(opts, 'step', ['method']);
    // dt is a duration. Unchecked, step(-1) stepped back to t = −1 and converged, and step('1')
    // made the clock the string '-11' (the state audit). As advance() does with its tEnd, a NaN
    // (an animation's first frame) or 0 does nothing, and says so.
    const length = `dt is the step's length in s, a positive number (the device is at t = ${this.solver.time})`;
    if (typeof dt !== 'number' || dt < 0 || dt === Infinity) throw new SolverError(`step(${typeof dt === 'string' ? `'${dt}'` : dt}): ${length}${dt === Infinity ? '. For the long-time limit, solve() the steady state' : ''}`);
    if (!(dt > 0)) return makeSolution(this.solver, { converged: true, iterations: 0, substeps: 0, stopped: `step(${dt}) did nothing: ${length}` });
    this._settleRegrid(); // (as advance() does: stepped unsolved, 27 A/m² where 68.8 is right)
    this.solver._freshDiagnostics();
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
    this._touched = true;
    options(opts, 'advance', ['tol', 'dt0', 'dtMax', 'budgetMs', 'maxSteps', 'method', 'probes']);
    // A steady state interpolated onto a new grid isn't quite one there (a strictly neutral
    // region's interpolated levels miss neutrality), and a transient's first short steps can't
    // close that gap: storage resists it (Fable round 2's cold user: set({ grid }) on a solved
    // cell, then advance(), failed down to 1e-17 s). Its steady state on the new grid is what it
    // stands for: solve that first. (Mid-run, an interpolated transient state steps on as it is.)
    this._settleRegrid();
    this.solver._freshDiagnostics();
    return makeSolution(this.solver, this.solver.integrate(tEnd, opts));
  }

  // A steady state carried onto a new grid, solved there (once; see advance).
  _settleRegrid() {
    if (!this._solver?.regriddedSteady) return;
    this._solver.regriddedSteady = false;
    const snap = this._solver._snapshot();
    if (!this._solver.solveSteady().converged) this._solver._restore(snap);
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
    this._touched = true;
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
    return { def: this.def, model: this.model, solver, drives: this.model.terminals.map((t) => t.drive), snap: solver._snapshot() };
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
