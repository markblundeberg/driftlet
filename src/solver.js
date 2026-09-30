// Discretisation and nonlinear solver.
//
// Unknowns, per solver block (block size M = 1 + nSpecies):
//   grid node: [φ̂, η_1 … η_n], with φ̂ = Fφ/RT and η_i = μ̄_i/RT
//   flux node: [D, N_1 … N_n], the displacement and particle fluxes through an interface
// Every region boundary is a doubled grid node (one per side) with a zero-volume flux node
// between them in the linear system, so the Jacobian stays block-tridiagonal.
//
// Balance rows (node g, box volume v per unit area):
//   φ:   D_out − D_in − v·(F Σ z_i c_i + ρ_fixed) = 0
//   i:   v·(c_i − c_i,old)/dt + N_out − N_in = 0
// Fluxes along segments are Scharfetter–Gummel. Each flux is computed once per segment and
// added with opposite signs to both neighbours, so sums over boxes telescope exactly.

import { BlockTridiagonal } from './blockTridiagonal.js';
import { bernoulli, bernoulliDerivative } from './bernoulli.js';
import { EPS0, FARADAY } from './constants.js';

export class SolverError extends Error {
  constructor(message, details) {
    super(message);
    this.name = 'SolverError';
    this.details = details;
  }
}

export class Solver {
  constructor(model) {
    this.model = model;
    const { grid, species, materials, regions, contacts } = model;
    const n = species.length;
    const M = n + 1;
    const nNodes = grid.nNodes;
    const nFaces = regions.length - 1;
    const nB = nNodes + nFaces;
    this.n = n;
    this.M = M;
    this.nNodes = nNodes;
    this.nFaces = nFaces;
    this.nB = nB;
    this.VT = model.RT / FARADAY; // thermal voltage, V

    for (const side of ['left', 'right']) {
      contacts[side].species.forEach((link, i) => {
        if (link.type !== 'blocked') {
          throw new SolverError(`contacts.${side}.species.${species[i].name}: '${link.type}' links are not implemented yet`);
        }
      });
      if (contacts[side].phi.type === 'neutral') {
        throw new SolverError(`contacts.${side}.phi: 'neutral' links are not implemented yet`);
      }
    }

    // Solver block of each grid node and of each interface flux node.
    this.blockOfNode = new Int32Array(nNodes);
    for (let g = 0; g < nNodes; g++) this.blockOfNode[g] = g + grid.nodeRegion[g];
    this.blockOfFace = new Int32Array(nFaces);
    for (let f = 0; f < nFaces; f++) this.blockOfFace[f] = grid.regionEnd[f] + f + 1;

    // Per-node material data, flattened [g·n + i].
    this.present = new Uint8Array(nNodes * n);
    this.cRef = new Float64Array(nNodes * n);
    this.mu0hat = new Float64Array(nNodes * n); // μ°/RT
    this.rhoFixed = new Float64Array(nNodes);
    this.z = Int32Array.from(species, (s) => s.z);
    for (let g = 0; g < nNodes; g++) {
      const reg = regions[grid.nodeRegion[g]];
      const mat = materials[reg.material];
      this.rhoFixed[g] = reg.fixedCharge;
      for (let i = 0; i < n; i++) {
        this.present[g * n + i] = mat.present[i];
        this.cRef[g * n + i] = mat.cRef[i];
        this.mu0hat[g * n + i] = mat.mu0[i] / model.RT;
      }
    }

    this.sys = new BlockTridiagonal(nB, M);
    this.u = new Float64Array(nB * M);
    this.res = new Float64Array(nB * M);
    this.delta = new Float64Array(nB * M);
    this.uPrev = new Float64Array(nB * M);
    this.c = new Float64Array(nNodes * n);
    this.cOld = new Float64Array(nNodes * n);
    this.time = 0;

    this._findStretches();
    this.initFromComposition();
    this.referenceAmounts = this.stretches.map((st) => this.amount(st));
  }

  // Connected stretches of regions where a species is present. A stretch not connected to
  // a contact (and, later, with no reaction) is a conserved spectator.
  _findStretches() {
    const { model, n } = this;
    const { regions, materials, contacts, grid } = model;
    const last = regions.length - 1;
    this.stretches = [];
    for (let i = 0; i < n; i++) {
      let r = 0;
      while (r <= last) {
        if (!materials[regions[r].material].present[i]) { r++; continue; }
        const r0 = r;
        while (r + 1 <= last && materials[regions[r + 1].material].present[i]) r++;
        const leftOpen = r0 === 0 && contacts.left.species[i].type !== 'blocked';
        const rightOpen = r === last && contacts.right.species[i].type !== 'blocked';
        this.stretches.push({
          species: i,
          regions: [r0, r],
          nodes: [grid.regionStart[r0], grid.regionEnd[r]],
          spectator: !(leftOpen || rightOpen),
        });
        r++;
      }
    }
  }

  /** Total amount (mol per unit area) of a stretch's species in the current state. */
  amount(stretch) {
    const { n, c, model } = this;
    const vol = model.grid.vol;
    let s = 0;
    for (let g = stretch.nodes[0]; g <= stretch.nodes[1]; g++) s += vol[g] * c[g * n + stretch.species];
    return s;
  }

  /** Set the state from each region's c0, with φ stepping by the interface dipoles. */
  initFromComposition() {
    const { model, n, M, u } = this;
    const grid = model.grid;
    const { regions, species, interfaces } = model;
    u.fill(0);
    let phiHat = 0;
    for (let r = 0; r < regions.length; r++) {
      if (r > 0) phiHat += interfaces[r - 1].dipole / this.VT;
      const reg = regions[r];
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const b = this.blockOfNode[g];
        u[b * M] = phiHat;
        for (let i = 0; i < n; i++) {
          if (!this.present[g * n + i]) continue;
          const c0 = reg.c0[i];
          if (!(c0 > 0)) {
            throw new SolverError(`regions[${r}].c0.${species[i].name}: an initial concentration is needed`);
          }
          u[b * M + 1 + i] = Math.log(c0 / this.cRef[g * n + i]) + this.mu0hat[g * n + i] + species[i].z * phiHat;
        }
      }
    }
    this.computeConcentrations();
  }

  computeConcentrations() {
    const { n, M, u, c, z } = this;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      const phiHat = u[b * M];
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        c[k] = this.present[k] ? this.cRef[k] * Math.exp(u[b * M + 1 + i] - this.mu0hat[k] - z[i] * phiHat) : 0;
      }
    }
  }

  // Add ∂(row rb, slot rs)/∂(block cb, slot cs) to the Jacobian.
  _j(rb, rs, cb, cs, v) {
    const M = this.M, o = rb * M * M + rs * M + cs;
    if (cb === rb) this.sys.B[o] += v;
    else if (cb === rb - 1) this.sys.A[o] += v;
    else if (cb === rb + 1) this.sys.C[o] += v;
    else throw new Error(`internal: non-tridiagonal coupling ${rb}→${cb}`);
  }

  /** Assemble residual and Jacobian for a backward-Euler step of size dt. */
  assemble(dt) {
    const { model, n, M, u, res, c, cOld, z, VT, sys } = this;
    const { grid, materials, regions, interfaces, contacts } = model;
    const F = FARADAY;
    sys.clear();
    res.fill(0);
    this.computeConcentrations();

    // Node terms: storage, space charge, identity rows for absent species.
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g], v = grid.vol[g];
      let q = this.rhoFixed[g], dq = 0;
      for (let i = 0; i < n; i++) {
        const k = g * n + i, r = 1 + i;
        if (!this.present[k]) {
          this._j(b, r, b, r, 1);
          continue;
        }
        const ck = c[k];
        res[b * M + r] += (v * (ck - cOld[k])) / dt;
        this._j(b, r, b, r, (v * ck) / dt);
        this._j(b, r, b, 0, (-v * z[i] * ck) / dt);
        q += F * z[i] * ck;
        dq += F * z[i] * z[i] * ck;
        if (z[i] !== 0) this._j(b, 0, b, r, -v * F * z[i] * ck);
      }
      res[b * M] -= v * q;
      this._j(b, 0, b, 0, v * dq);
    }

    // Ordinary segments: displacement and Scharfetter–Gummel fluxes.
    for (let s = 0; s < this.nNodes - 1; s++) {
      const reg = grid.segRegion[s];
      if (reg < 0) continue;
      const mat = materials[regions[reg].material];
      const h = grid.segLength[s];
      const bL = this.blockOfNode[s], bR = bL + 1;
      const phiL = u[bL * M], phiR = u[bR * M];

      const k = (mat.epsr * EPS0 * VT) / h;
      const D = -k * (phiR - phiL);
      res[bL * M] += D;
      res[bR * M] -= D;
      this._j(bL, 0, bL, 0, k);
      this._j(bL, 0, bR, 0, -k);
      this._j(bR, 0, bL, 0, -k);
      this._j(bR, 0, bR, 0, k);

      for (let i = 0; i < n; i++) {
        if (!mat.present[i] || mat.D[i] === 0) continue;
        const r = 1 + i, zi = z[i];
        const cL = c[s * n + i], cR = c[(s + 1) * n + i];
        const g = mat.D[i] / h;
        const d = zi * (phiR - phiL);
        const Bp = bernoulli(d), Bm = bernoulli(-d);
        const N = g * (Bp * cL - Bm * cR);
        const dNdd = g * (bernoulliDerivative(d) * cL + bernoulliDerivative(-d) * cR);
        const dNdEtaL = g * Bp * cL;
        const dNdEtaR = -g * Bm * cR;
        const dNdPhiL = -zi * dNdEtaL - zi * dNdd;
        const dNdPhiR = -zi * dNdEtaR + zi * dNdd;
        res[bL * M + r] += N;
        res[bR * M + r] -= N;
        this._j(bL, r, bL, r, dNdEtaL);
        this._j(bL, r, bR, r, dNdEtaR);
        this._j(bL, r, bL, 0, dNdPhiL);
        this._j(bL, r, bR, 0, dNdPhiR);
        this._j(bR, r, bL, r, -dNdEtaL);
        this._j(bR, r, bR, r, -dNdEtaR);
        this._j(bR, r, bL, 0, -dNdPhiL);
        this._j(bR, r, bR, 0, -dNdPhiR);
      }
    }

    // Interfaces: the flux node carries D and N_i; its rows are the interface laws.
    for (let f = 0; f < this.nFaces; f++) {
      const bf = this.blockOfFace[f], bL = bf - 1, bR = bf + 1;
      const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
      // φ: fixed offset (the dipole); displacement passes through.
      res[bf * M] = u[bR * M] - u[bL * M] - interfaces[f].dipole / VT;
      this._j(bf, 0, bR, 0, 1);
      this._j(bf, 0, bL, 0, -1);
      res[bL * M] += u[bf * M];
      this._j(bL, 0, bf, 0, 1);
      res[bR * M] -= u[bf * M] + interfaces[f].sheetCharge;
      this._j(bR, 0, bf, 0, -1);
      for (let i = 0; i < n; i++) {
        const r = 1 + i;
        const pL = this.present[gL * n + i], pR = this.present[gR * n + i];
        if (pL && pR) {
          // Local equilibrium: μ̄ continuous.
          res[bf * M + r] = u[bR * M + r] - u[bL * M + r];
          this._j(bf, r, bR, r, 1);
          this._j(bf, r, bL, r, -1);
          res[bL * M + r] += u[bf * M + r];
          this._j(bL, r, bf, r, 1);
          res[bR * M + r] -= u[bf * M + r];
          this._j(bR, r, bf, r, -1);
        } else {
          res[bf * M + r] = u[bf * M + r];
          this._j(bf, r, bf, r, 1);
        }
      }
    }

    // Contacts: electrostatic links (species links are all blocked for now).
    const first = 0, last = this.nB - 1;
    const pl = contacts.left.phi, pr = contacts.right.phi;
    if (pl.type === 'capacitive') {
      const phiG = pl.V - pl.zeroCharge;
      res[first * M] -= pl.C * (phiG - VT * u[first * M]); // −D_in
      this._j(first, 0, first, 0, pl.C * VT);
    }
    if (pr.type === 'capacitive') {
      const phiG = pr.V - pr.zeroCharge;
      res[last * M] += pr.C * (VT * u[last * M] - phiG); // +D_out
      this._j(last, 0, last, 0, pr.C * VT);
    }
  }

  // Scale every row by its largest Jacobian entry (in place, residual too).
  _equilibrate() {
    const { M, nB, sys, res } = this;
    const { A, B, C } = sys;
    const mm = M * M;
    for (let b = 0; b < nB; b++) {
      for (let r = 0; r < M; r++) {
        const o = b * mm + r * M;
        let mx = 0;
        for (let k = 0; k < M; k++) {
          mx = Math.max(mx, Math.abs(A[o + k]), Math.abs(B[o + k]), Math.abs(C[o + k]));
        }
        if (mx === 0) continue;
        const s = 1 / mx;
        for (let k = 0; k < M; k++) {
          A[o + k] *= s;
          B[o + k] *= s;
          C[o + k] *= s;
        }
        res[b * M + r] *= s;
      }
    }
  }

  // Largest update among the potential-like unknowns (φ̂ and η at grid nodes).
  _maxPotentialStep(delta) {
    const { n, M } = this;
    let mx = 0;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      mx = Math.max(mx, Math.abs(delta[b * M]));
      for (let i = 0; i < n; i++) {
        if (this.present[g * n + i]) mx = Math.max(mx, Math.abs(delta[b * M + 1 + i]));
      }
    }
    return mx;
  }

  /**
   * Newton iteration for one backward-Euler step from the current cOld.
   * @returns {{converged: boolean, iterations: number, history: number[]}}
   */
  newton(dt, { maxIter = 60, tol = 1e-10, maxStep = 10 } = {}) {
    const { u, delta, res } = this;
    const history = [];
    for (let it = 1; it <= maxIter; it++) {
      this.assemble(dt);
      this._equilibrate();
      let rmax = 0;
      for (let k = 0; k < res.length; k++) rmax = Math.max(rmax, Math.abs(res[k]));
      try {
        this.sys.factor();
      } catch (err) {
        return { converged: false, iterations: it, history, error: err.message };
      }
      this.sys.solve(res, delta);
      const step = this._maxPotentialStep(delta);
      history.push(step);
      if (!Number.isFinite(step)) return { converged: false, iterations: it, history, error: 'non-finite update' };
      const alpha = step > maxStep ? maxStep / step : 1;
      for (let k = 0; k < u.length; k++) u[k] -= alpha * delta[k];
      if (alpha === 1 && step < tol) {
        this.computeConcentrations();
        return { converged: true, iterations: it, history, residual: rmax };
      }
    }
    this.computeConcentrations();
    return { converged: false, iterations: maxIter, history };
  }

  /** One backward-Euler step. On failure the state is restored. */
  step(dt, opts) {
    this.uPrev.set(this.u);
    this.computeConcentrations();
    this.cOld.set(this.c);
    const result = this.newton(dt, opts);
    if (result.converged) {
      this.time += dt;
    } else {
      this.u.set(this.uPrev);
      this.computeConcentrations();
    }
    return result;
  }

  // Slowest diffusion time across the device, used to size "giant" steps.
  slowestTime() {
    const { model } = this;
    let dMin = Infinity;
    for (const mat of model.materials) {
      for (let i = 0; i < this.n; i++) if (mat.present[i] && mat.D[i] > 0) dMin = Math.min(dMin, mat.D[i]);
    }
    const L = model.grid.length;
    return Number.isFinite(dMin) ? (L * L) / dMin : 1;
  }

  /**
   * Steady state by backward-Euler steps at a dt far beyond the slowest time constant.
   * The storage term keeps spectator amounts exactly conserved. If Newton fails at the
   * giant dt, ramp dt up from a small value instead (pseudo-transient continuation).
   */
  solveSteady({ maxSteps = 60, tol = 1e-10 } = {}) {
    const tau = this.slowestTime();
    const giant = 1e6 * tau;
    let dt = giant;
    let totalIter = 0, steps = 0;
    const history = [];
    while (steps < maxSteps) {
      steps++;
      const r = this.step(dt);
      totalIter += r.iterations;
      history.push({ dt, converged: r.converged, iterations: r.iterations });
      if (!r.converged) {
        dt = dt === giant ? tau * 1e-6 : dt / 4;
        if (dt < tau * 1e-15) break;
        continue;
      }
      if (dt >= giant && this._maxPotentialStep(this._diff()) < tol) {
        return { converged: true, steps, iterations: totalIter, history };
      }
      dt = Math.min(dt * 10, giant);
    }
    return { converged: false, steps, iterations: totalIter, history };
  }

  _diff() {
    const d = this.delta; // reuse as scratch
    for (let k = 0; k < d.length; k++) d[k] = this.u[k] - this.uPrev[k];
    return d;
  }
}
