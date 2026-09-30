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
        if (link.type !== 'blocked' && link.type !== 'fixed') {
          throw new SolverError(`contacts.${side}.species.${species[i].name}: '${link.type}' links are not implemented yet`);
        }
      });
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
    this.nodeMaterial = new Int32Array(nNodes);
    this.z = Int32Array.from(species, (s) => s.z);
    for (let g = 0; g < nNodes; g++) {
      const reg = regions[grid.nodeRegion[g]];
      const mat = materials[reg.material];
      this.rhoFixed[g] = reg.fixedCharge;
      this.nodeMaterial[g] = reg.material;
      for (let i = 0; i < n; i++) {
        this.present[g * n + i] = mat.present[i];
        this.cRef[g * n + i] = mat.cRef[i];
        this.mu0hat[g * n + i] = mat.mu0[i] / model.RT;
      }
    }

    this.sys = new BlockTridiagonal(nB, M);
    // Unknowns as compensated double-doubles, u + uLo. Only differences of η need the extra
    // precision: a majority carrier carrying a small current has a quasi-Fermi step between
    // nodes far below the ulp of η itself (e.g. 1e-19 vs 3.5e-15), and would otherwise carry
    // exactly zero current there. The Jacobian and linear solve stay in plain doubles.
    this.u = new Float64Array(nB * M);
    this.uLo = new Float64Array(nB * M);
    this.res = new Float64Array(nB * M);
    this.delta = new Float64Array(nB * M);
    this.uPrev = new Float64Array(nB * M);
    this.uPrevLo = new Float64Array(nB * M);
    this.c = new Float64Array(nNodes * n);
    this.cOld = new Float64Array(nNodes * n);
    this.time = 0;
    this.lastDt = Infinity;
    // Contact bookkeeping, filled by assemble(): particle flux toward +x through each contact,
    // and the displacement there (the metal's surface charge for a neutral link).
    this.contactFlux = { left: new Float64Array(n), right: new Float64Array(n) };
    this.contactD = { left: 0, right: 0 };
    this.contactDOld = { left: 0, right: 0 };
    // Total current through the last segment and its derivatives (for current/load circuits).
    this.segI = 0;
    this.segIJac = new Float64Array(2 * M); // [∂/∂(block last−1 slots), ∂/∂(block last slots)]
    this.segDOld = 0;

    this._findStretches();
    this.initFromComposition();
    this.referenceAmounts = this.stretches.map((st) => this.amount(st));
    // ∫ (flux in − flux out) dt through the contacts, per stretch, since the reference.
    this.boundaryIntake = new Float64Array(this.stretches.length);
  }

  // Connected stretches of regions where a species is present. A stretch not connected to
  // a contact (and, later, with no reaction) is a conserved spectator.
  _findStretches() {
    const { model, n } = this;
    const { regions, materials, contacts, grid } = model;
    const last = regions.length - 1;
    this.stretches = [];
    this.stretchOf = new Int32Array(regions.length * n).fill(-1);
    for (let i = 0; i < n; i++) {
      let r = 0;
      while (r <= last) {
        if (!materials[regions[r].material].present[i]) { r++; continue; }
        const r0 = r;
        while (r + 1 <= last && materials[regions[r + 1].material].present[i]) r++;
        const leftOpen = r0 === 0 && contacts.left.species[i].type !== 'blocked';
        const rightOpen = r === last && contacts.right.species[i].type !== 'blocked';
        let reactive = false;
        for (let q = r0; q <= r; q++) {
          this.stretchOf[q * n + i] = this.stretches.length;
          const m = regions[q].material;
          for (const rx of model.reactions) {
            if (rx.kf[m] > 0 && [...rx.reactants, ...rx.products].some((p) => p.i === i)) reactive = true;
          }
        }
        const connected = leftOpen || rightOpen;
        this.stretches.push({
          species: i,
          regions: [r0, r],
          nodes: [grid.regionStart[r0], grid.regionEnd[r]],
          connected,
          reactive,
          // Conserved on its own: not fed by a contact, not made or consumed by a reaction.
          spectator: !connected && !reactive,
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

  /** η_i/RT that a fixed contact link imposes (or NaN if the link isn't fixed). */
  contactEta(side, i) {
    const ct = this.model.contacts[side], link = ct.species[i];
    if (link.type !== 'fixed') return NaN;
    const z = this.z[i];
    return z === 0 ? link.mu / this.model.RT : (z * (ct.V + link.offset)) / this.VT;
  }

  /**
   * Cold start. Species connected to a contact take that contact's level; spectators take
   * their region's c0 (which fixes their conserved amount). φ in each region is then chosen
   * for local neutrality, or continued across the interface dipole if nothing there responds.
   */
  initFromComposition() {
    const { model, n, M, u, z } = this;
    const { regions, species, interfaces, materials } = model;
    const grid = model.grid;
    u.fill(0);
    this.uLo.fill(0);
    const eta = new Float64Array(n), cFix = new Float64Array(n), mode = new Int8Array(n); // 1 level, 2 amount
    let phiHat = 0;
    for (let r = 0; r < regions.length; r++) {
      if (r > 0) phiHat += interfaces[r - 1].dipole / this.VT;
      const reg = regions[r], mat = materials[reg.material];
      mode.fill(0);
      for (let i = 0; i < n; i++) {
        if (!mat.present[i]) continue;
        const st = this.stretches[this.stretchOf[r * n + i]];
        if (!st.connected) {
          if (!(reg.c0[i] > 0)) {
            throw new SolverError(
              `regions[${r}].c0.${species[i].name}: a species not connected to a contact needs its initial concentration`,
            );
          }
          mode[i] = 2;
          cFix[i] = reg.c0[i];
        } else {
          const left = st.regions[0] === 0 ? this.contactEta('left', i) : NaN;
          eta[i] = Number.isFinite(left) ? left : this.contactEta('right', i);
          mode[i] = 1;
        }
      }
      // Net charge (mol/m³) at trial φ̂; decreasing in φ̂ wherever a level-fixed ion responds.
      const charge = (ph) => {
        let q = reg.fixedCharge / FARADAY;
        for (let i = 0; i < n; i++) {
          if (mode[i] === 2) q += z[i] * cFix[i];
          else if (mode[i] === 1 && z[i] !== 0) {
            const ex = Math.min(700, Math.max(-700, eta[i] - mat.mu0[i] / model.RT - z[i] * ph));
            q += z[i] * mat.cRef[i] * Math.exp(ex);
          }
        }
        return q;
      };
      if (mode.some((m, i) => m === 1 && z[i] !== 0)) {
        let lo = phiHat - 1, hi = phiHat + 1;
        while (charge(lo) < 0 && lo > -1e4) lo -= 2 * (hi - lo);
        while (charge(hi) > 0 && hi < 1e4) hi += 2 * (hi - lo);
        for (let it = 0; it < 200 && hi - lo > 1e-12; it++) {
          const m = 0.5 * (lo + hi);
          if (charge(m) > 0) lo = m;
          else hi = m;
        }
        phiHat = 0.5 * (lo + hi);
      }
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const b = this.blockOfNode[g];
        u[b * M] = phiHat;
        for (let i = 0; i < n; i++) {
          if (mode[i] === 1) u[b * M + 1 + i] = eta[i];
          else if (mode[i] === 2) u[b * M + 1 + i] = Math.log(cFix[i] / mat.cRef[i]) + mat.mu0[i] / model.RT + z[i] * phiHat;
        }
      }
    }
    this.computeConcentrations();
  }

  computeConcentrations() {
    const { n, M, u, uLo, c, z } = this;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      const phiHat = u[b * M];
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        c[k] = this.present[k] ? this.cRef[k] * Math.exp(u[b * M + 1 + i] + uLo[b * M + 1 + i] - this.mu0hat[k] - z[i] * phiHat) : 0;
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
    const { model, n, M, u, uLo, res, c, cOld, z, VT, sys } = this;
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

      // Bulk reactions: r = k_f Π c_R^ν · (−expm1(−a)), with a = A/RT from the (compensated) η.
      const m = this.nodeMaterial[g];
      for (const rx of model.reactions) {
        const kf = rx.kf[m];
        if (!(kf > 0)) continue;
        let P = kf, zsum = 0, aHi = rx.fixedA, aLo = 0;
        for (const { i, nu } of rx.reactants) {
          P *= c[g * n + i] ** nu;
          zsum += nu * z[i];
          aHi += nu * u[b * M + 1 + i];
          aLo += nu * uLo[b * M + 1 + i];
        }
        for (const { i, nu } of rx.products) {
          aHi -= nu * u[b * M + 1 + i];
          aLo -= nu * uLo[b * M + 1 + i];
        }
        const a = aHi + aLo;
        const f = -Math.expm1(-a); // 1 − e^{−a}
        const rate = P * f;
        const dfda = 1 - f; // e^{−a}
        // Row contributions: reactants consumed (+v·ν·r in their balance), products made (−).
        const add = (list, sign) => {
          for (const { i, nu } of list) {
            const row = 1 + i;
            res[b * M + row] += sign * v * nu * rate;
            this._j(b, row, b, 0, sign * v * nu * (-zsum * rate));
            for (const { i: j, nu: nj } of rx.reactants) this._j(b, row, b, 1 + j, sign * v * nu * nj * (rate + P * dfda));
            for (const { i: j, nu: nj } of rx.products) this._j(b, row, b, 1 + j, sign * v * nu * (-nj * P * dfda));
          }
        };
        add(rx.reactants, 1);
        add(rx.products, -1);
      }
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
      const lastSeg = s === this.nNodes - 2;
      if (lastSeg) {
        // Displacement current through this segment; conduction is added per species below.
        this.segD = D;
        this.segI = Number.isFinite(dt) ? (D - this.segDOld) / dt : 0;
        this.segIJac.fill(0);
        if (Number.isFinite(dt)) {
          this.segIJac[0] = k / dt;
          this.segIJac[M] = -k / dt;
        }
      }
      res[bL * M] += D;
      res[bR * M] -= D;
      this._j(bL, 0, bL, 0, k);
      this._j(bL, 0, bR, 0, -k);
      this._j(bR, 0, bL, 0, -k);
      this._j(bR, 0, bR, 0, k);

      for (let i = 0; i < n; i++) {
        if (!mat.present[i] || mat.D[i] === 0) continue;
        const r = 1 + i, zi = z[i];
        // SG flux N = g[B(Δ)c_L − B(−Δ)c_R], rewritten with B(−Δ) = B(Δ)e^Δ and
        // c_R e^Δ = c_L e^{Δη} as N = −g·B(Δ)·c_L·expm1(Δη). This is precise relative to the
        // quasi-Fermi difference Δη, so tiny fluxes (e.g. majority carriers carrying a small
        // current) don't vanish in the cancellation of two huge drift and diffusion terms.
        const cL = c[s * n + i];
        const g = mat.D[i] / h;
        const d = zi * (phiR - phiL);
        const deta = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]);
        const E = Math.expm1(deta);
        const gBc = g * bernoulli(d) * cL;
        const N = -gBc * E;
        const dNdd = -g * bernoulliDerivative(d) * cL * E;
        const dNdEtaL = gBc;
        const dNdEtaR = -gBc * (E + 1);
        const dNdPhiL = -zi * dNdd + zi * gBc * E; // via Δ, and via c_L ∝ e^{−zφ̂_L}
        const dNdPhiR = zi * dNdd;
        res[bL * M + r] += N;
        res[bR * M + r] -= N;
        if (lastSeg) {
          const q = F * zi;
          this.segI += q * N;
          this.segIJac[r] += q * dNdEtaL;
          this.segIJac[M + r] += q * dNdEtaR;
          this.segIJac[0] += q * dNdPhiL;
          this.segIJac[M] += q * dNdPhiR;
        }
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
      res[bf * M] = u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) - interfaces[f].dipole / VT;
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
          res[bf * M + r] = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]);
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

    // Contacts. The boundary node's balance rows are complete except for the flux through the
    // contact, so before a row is replaced by a contact condition its residual *is* that flux:
    // entering at the left (+res), leaving at the right (−res), both counted toward +x.
    const circuit = model.circuit;
    for (const side of ['left', 'right']) {
      const ct = contacts[side];
      const g = side === 'left' ? 0 : this.nNodes - 1;
      const b = this.blockOfNode[g];
      const sgn = side === 'left' ? 1 : -1;
      const flux = this.contactFlux[side];
      // In current/load mode the right terminal floats: its voltage is read off the terminal
      // species' own level at the contact node, and its row becomes the circuit law instead.
      const floating = side === 'right' && circuit.mode !== 'voltage';
      const t = ct.terminal;
      let Vt = ct.V;
      if (floating) {
        const ot = b * M + 1 + t;
        Vt = (VT * (u[ot] + uLo[ot])) / z[t] - ct.species[t].offset;
        this.terminalV = Vt;
      }
      for (let i = 0; i < n; i++) {
        flux[i] = sgn * res[b * M + 1 + i];
        if (ct.species[i].type !== 'fixed') continue;
        const o = b * M + 1 + i;
        this._replaceRow(b, 1 + i);
        if (floating && i === t) {
          // I_segment − I_circuit(V_t) = 0. The last segment's total current equals the
          // terminal current exactly (box balance plus Poisson, differenced in time).
          const I = circuit.mode === 'current' ? circuit.I : (Vt - contacts.left.V - circuit.V) / circuit.R;
          res[o] = this.segI - I;
          for (let r = 0; r < M; r++) {
            this._j(b, 1 + i, b - 1, r, this.segIJac[r]);
            this._j(b, 1 + i, b, r, this.segIJac[M + r]);
          }
          if (circuit.mode === 'load') this._j(b, 1 + i, b, 1 + t, -VT / (z[t] * circuit.R));
        } else if (floating && ct.species[i].mu === undefined) {
          // Charged species tied to the floating terminal: η_i = z_i (V_t + offset_i)/V_T.
          res[o] = u[o] + uLo[o] - (z[i] * (Vt + ct.species[i].offset)) / VT;
          this._j(b, 1 + i, b, 1 + i, 1);
          this._j(b, 1 + i, b, 1 + t, -z[i] / z[t]);
        } else {
          this._j(b, 1 + i, b, 1 + i, 1);
          res[o] = u[o] - this.contactEta(side, i) + uLo[o];
        }
      }
      const link = ct.phi;
      if (link.type === 'capacitive') {
        // Gate (or metal across a Stern layer) at φ_g = V − zeroCharge; D = C·(φ_g − φ) inward.
        const phiG = ct.V - link.zeroCharge;
        const D = link.C * (phiG - VT * u[b * M]) * sgn;
        res[b * M] -= sgn * D;
        this._j(b, 0, b, 0, link.C * VT);
        this.contactD[side] = D;
      } else if (link.type === 'neutral') {
        this.contactD[side] = sgn * res[b * M];
        this._replaceRow(b, 0);
        let q = this.rhoFixed[g];
        let dq = 0;
        for (let i = 0; i < n; i++) {
          const k = g * n + i;
          if (!this.present[k] || z[i] === 0) continue;
          q += F * z[i] * c[k];
          dq -= F * z[i] * z[i] * c[k];
          this._j(b, 0, b, 1 + i, F * z[i] * c[k]);
        }
        this._j(b, 0, b, 0, dq);
        res[b * M] = q;
      } else {
        this.contactD[side] = 0;
      }
    }
  }

  // Displacement toward +x through the last segment, from the current state.
  _lastSegmentD() {
    const { model, M, u, VT } = this;
    const grid = model.grid, s = this.nNodes - 2;
    const mat = model.materials[model.regions[grid.segRegion[s]].material];
    const b = this.blockOfNode[s];
    return (-(mat.epsr * EPS0 * VT) / grid.segLength[s]) * (u[(b + 1) * M] - u[b * M]);
  }

  // Zero one row of the Jacobian (all three blocks), ready to be replaced.
  _replaceRow(b, r) {
    const M = this.M, o = b * M * M + r * M;
    for (let k = 0; k < M; k++) {
      this.sys.A[o + k] = 0;
      this.sys.B[o + k] = 0;
      this.sys.C[o + k] = 0;
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
      // Give up early on clear divergence; the caller will take a smaller step instead.
      if (step > 1e4 || (it > 6 && step > 10 * history[0])) {
        this.computeConcentrations();
        return { converged: false, iterations: it, history, error: 'diverging' };
      }
      const alpha = step > maxStep ? maxStep / step : 1;
      this._addToState(delta, -alpha);
      if (alpha === 1 && step < tol) {
        this.computeConcentrations();
        return { converged: true, iterations: it, history, residual: rmax };
      }
    }
    this.computeConcentrations();
    return { converged: false, iterations: maxIter, history };
  }

  /**
   * Advance by dt, splitting the interval into halves (recursively) wherever Newton fails.
   * The end time is always t + dt; the result reports the substeps taken.
   */
  advance(dt, opts, depth = 0) {
    const r = this.step(dt, opts);
    if (r.converged || depth >= 30) return { ...r, substeps: 1 };
    const a = this.advance(dt / 2, opts, depth + 1);
    if (!a.converged) return a;
    const b = this.advance(dt / 2, opts, depth + 1);
    return {
      converged: b.converged,
      iterations: r.iterations + a.iterations + b.iterations,
      history: b.history,
      substeps: a.substeps + b.substeps,
    };
  }

  // u += scale·d in compensated arithmetic (Knuth two-sum, then renormalise hi/lo).
  _addToState(d, scale) {
    const { u, uLo } = this;
    for (let k = 0; k < u.length; k++) {
      const a = u[k], b = scale * d[k];
      const s = a + b, bb = s - a;
      const err = a - (s - bb) + (b - bb);
      const lo = uLo[k] + err;
      const hi = s + lo;
      u[k] = hi;
      uLo[k] = lo - (hi - s);
    }
  }

  /** One backward-Euler step. On failure the state is restored. */
  step(dt, opts) {
    this.uPrev.set(this.u);
    this.uPrevLo.set(this.uLo);
    this.computeConcentrations();
    this.cOld.set(this.c);
    this.segDOld = this._lastSegmentD();
    // Contact displacement before the step: as it was at the end of the previous step (under
    // the parameters then), so a gate-voltage change shows up as displacement current.
    if (!this.contactDEnd) {
      this.assemble(dt);
      this.contactDEnd = { ...this.contactD };
    }
    const DOld = { ...this.contactDEnd };
    const result = this.newton(dt, opts);
    if (result.converged) {
      this.time += dt;
      this.lastDt = dt;
      this.contactDOld = DOld;
      this._accumulateBoundaryIntake(dt);
    } else {
      this.u.set(this.uPrev);
      this.uLo.set(this.uPrevLo);
      this.computeConcentrations();
    }
    return result;
  }

  // Add this step's contact fluxes (at the converged state) to each stretch's intake.
  _accumulateBoundaryIntake(dt) {
    this.assemble(dt);
    this.contactDEnd = { ...this.contactD };
    const last = this.model.regions.length - 1;
    this.stretches.forEach((st, k) => {
      if (!st.connected) return;
      let q = 0;
      if (st.regions[0] === 0) q += this.contactFlux.left[st.species];
      if (st.regions[1] === last) q -= this.contactFlux.right[st.species];
      this.boundaryIntake[k] += q * dt;
    });
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
  solveSteady({ maxSteps = 60, tol = 1e-13 } = {}) {
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
