// Butler–Volmer factor g(a) = e^{αa} − e^{−(1−α)a} = e^{−(1−α)a}·expm1(a) (precise near a = 0),
// and its derivative g′(a) = α e^{αa} + (1−α) e^{−(1−α)a}.
function bvFactor(a, alpha) {
  const e = Math.exp(-(1 - alpha) * a);
  return { g: e * Math.expm1(a), gp: alpha * Math.exp(alpha * a) + (1 - alpha) * e };
}

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
    // A floating terminal with no fixed species to read its voltage from gets its own block.
    const terminalUnknown = model.circuit.mode !== 'voltage' && model.circuit.terminalUnknown;
    const nB = nNodes + nFaces + (terminalUnknown ? 1 : 0);
    this.terminalBlock = terminalUnknown ? nB - 1 : -1;
    this.n = n;
    this.M = M;
    this.nNodes = nNodes;
    this.nFaces = nFaces;
    this.nB = nB;
    this.VT = model.RT / FARADAY; // thermal voltage, V


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

    // Regions where φ is undefined: no charged species, and either ε = 0 or nothing couples the
    // region electrostatically (neutral faces, no gate). φ there gets an identity row.
    this.phiUndefined = new Uint8Array(nNodes);
    regions.forEach((reg, r) => {
      const mat = materials[reg.material];
      if (species.some((sp, i) => mat.present[i] && sp.z !== 0)) return;
      const leftOpen = r === 0 ? contacts.left.phi.type === 'capacitive' : model.interfaces[r - 1].phi.type !== 'neutral';
      const rightOpen =
        r === regions.length - 1 ? contacts.right.phi.type === 'capacitive' : model.interfaces[r].phi.type !== 'neutral';
      if (mat.epsr === 0 || !(leftOpen || rightOpen)) {
        for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) this.phiUndefined[g] = 1;
      }
    });

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
    this._termI = 0;
    this._termIJac = new Float64Array(M + 1);
    this.contactDStart = { left: 0, right: 0 };

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
        while (r + 1 <= last && materials[regions[r + 1].material].present[i] && model.interfaces[r].links[i].type !== 'blocked') r++;
        const touches = (ct) => ct.species[i].type !== 'blocked' || ct.reactions.some((rx) => [...rx.reactants, ...rx.products].some((p) => p.i === i));
        const leftOpen = r0 === 0 && touches(contacts.left);
        const rightOpen = r === last && touches(contacts.right);
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
    if (this.terminalBlock >= 0) u[this.terminalBlock * M] = model.contacts.right.V / this.VT;
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
      if (this.phiUndefined[g]) {
        this._j(b, 0, b, 0, 1); // no charge responds and no field reaches: φ is not defined here
      } else {
        res[b * M] -= v * q;
        this._j(b, 0, b, 0, v * dq);
      }

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

      const k = this.phiUndefined[s] ? 0 : (mat.epsr * EPS0 * VT) / h; // ε = 0: no displacement
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
      const itf = interfaces[f];
      // φ law: pinned jump (dipole), Helmholtz capacitor, or no charge at all (neutral).
      const law = itf.phi.type;
      if (law === 'neutral') {
        res[bf * M] = u[bf * M]; // D = 0; the jump is whatever each side's neutrality needs
        this._j(bf, 0, bf, 0, 1);
      } else {
        const jump = u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) - itf.dipole / VT;
        if (law === 'dipole') {
          res[bf * M] = jump;
          this._j(bf, 0, bR, 0, 1);
          this._j(bf, 0, bL, 0, -1);
        } else {
          // D = −C (φ_R − φ_L − dipole): displacement toward +x drops across the layer.
          const kC = itf.phi.C * VT;
          res[bf * M] = u[bf * M] + kC * jump;
          this._j(bf, 0, bf, 0, 1);
          this._j(bf, 0, bR, 0, kC);
          this._j(bf, 0, bL, 0, -kC);
        }
      }
      res[bL * M] += u[bf * M];
      this._j(bL, 0, bf, 0, 1);
      res[bR * M] -= u[bf * M] + itf.sheetCharge;
      this._j(bR, 0, bf, 0, -1);
      for (let i = 0; i < n; i++) {
        const r = 1 + i, o = bf * M + r;
        const type = itf.links[i].type;
        const deta = u[bL * M + r] - u[bR * M + r] + (uLo[bL * M + r] - uLo[bR * M + r]); // η_L − η_R
        if (type === 'blocked') {
          res[o] = u[o];
          this._j(bf, r, bf, r, 1);
          continue;
        }
        if (type === 'equilibrium') {
          res[o] = -deta; // μ̄ continuous
          this._j(bf, r, bR, r, 1);
          this._j(bf, r, bL, r, -1);
        } else if (type === 'conductance') {
          // J = G (V_L − V_R), V = V_T η / z  ⇒  N = G V_T (η_L − η_R) / (z² F)
          const gG = (itf.links[i].G * VT) / (z[i] * z[i] * F);
          res[o] = u[o] - gG * deta;
          this._j(bf, r, bf, r, 1);
          this._j(bf, r, bL, r, -gG);
          this._j(bf, r, bR, r, gG);
        } else {
          res[o] = u[o]; // kinetic: N − Σ ν r, the rates are subtracted below
          this._j(bf, r, bf, r, 1);
        }
        res[bL * M + r] += u[o];
        this._j(bL, r, bf, r, 1);
        res[bR * M + r] -= u[o];
        this._j(bR, r, bf, r, -1);
      }
      for (const tr of itf.transfers) this._transfer(tr, f, bf, bL, bR);
    }

    // Contacts.
    const circuit = model.circuit;
    const tb = this.terminalBlock; // extra block carrying a floating terminal voltage, or −1
    this._termI = 0;
    this._termIJac.fill(0); // [∂/∂(last node slots), ∂/∂V_t]
    for (const side of ['left', 'right']) this._contact(side, dt);
    if (tb >= 0) {
      // Circuit law for the floating terminal: I_contact(V_t, last node) − I_circuit(V_t) = 0.
      const Vt = this.terminalV;
      let I = circuit.I, dIdV = 0;
      if (circuit.mode === 'load') {
        I = (Vt - contacts.left.V - circuit.V) / circuit.R;
        dIdV = 1 / circuit.R;
      }
      res[tb * M] = this._termI - I;
      for (let r = 0; r < M; r++) this._j(tb, 0, tb - 1, r, this._termIJac[r]);
      this._j(tb, 0, tb, 0, (this._termIJac[M] - dIdV) * VT);
      for (let r = 1; r < M; r++) {
        res[tb * M + r] = u[tb * M + r];
        this._j(tb, r, tb, r, 1);
      }
    }
  }

  // Kinetic transfer across interface f (Butler–Volmer form, forward = left to right):
  //   r = k0 Π [(c_L/c_ref,L)^{ν(1−α)} (c_R/c_ref,R)^{να}] (e^{αa} − e^{−(1−α)a}),
  //   a = Σ ν (η_L − η_R).  Each species' flux-node row gets −ν·r.
  _transfer(tr, f, bf, bL, bR) {
    const { n, M, u, uLo, c, z, res } = this;
    const gL = this.model.grid.regionEnd[f], gR = this.model.grid.regionStart[f + 1];
    const al = tr.alpha;
    let pref = tr.k0, aHi = 0, aLo = 0, zL = 0, zR = 0;
    for (const { i, nu } of tr.species) {
      pref *= (c[gL * n + i] / this.cRef[gL * n + i]) ** (nu * (1 - al)) * (c[gR * n + i] / this.cRef[gR * n + i]) ** (nu * al);
      aHi += nu * (u[bL * M + 1 + i] - u[bR * M + 1 + i]);
      aLo += nu * (uLo[bL * M + 1 + i] - uLo[bR * M + 1 + i]);
      zL += nu * (1 - al) * z[i];
      zR += nu * al * z[i];
    }
    const { g, gp } = bvFactor(aHi + aLo, al);
    const rate = pref * g;
    for (const { i, nu } of tr.species) {
      const row = 1 + i;
      res[bf * M + row] -= nu * rate;
      for (const { i: j, nu: nj } of tr.species) {
        this._j(bf, row, bL, 1 + j, -nu * (nj * (1 - al) * rate + pref * gp * nj));
        this._j(bf, row, bR, 1 + j, -nu * (nj * al * rate - pref * gp * nj));
      }
      this._j(bf, row, bL, 0, -nu * (-zL * rate));
      this._j(bf, row, bR, 0, -nu * (-zR * rate));
    }
  }

  // One contact: record the flux through its outer face, then add its exchange terms (electrode
  // reactions, conductance links), then apply fixed links and the φ link.
  _contact(side, dt) {
    const { model, n, M, u, uLo, res, c, z, VT } = this;
    const F = FARADAY;
    const contacts = model.contacts, circuit = model.circuit;
    const ct = contacts[side];
    const g = side === 'left' ? 0 : this.nNodes - 1;
    const b = this.blockOfNode[g];
    const sgn = side === 'left' ? 1 : -1;
    const flux = this.contactFlux[side];
    const floating = side === 'right' && circuit.mode !== 'voltage';
    const readout = floating && !circuit.terminalUnknown;
    const tb = this.terminalBlock;
    const t = ct.terminal;

    // Terminal voltage, and how it depends on the unknowns.
    let Vt = ct.V;
    if (readout) {
      const ot = b * M + 1 + t;
      Vt = (VT * (u[ot] + uLo[ot])) / z[t] - ct.species[t].offset;
    } else if (floating) {
      Vt = VT * u[tb * M];
    }
    if (floating) this.terminalV = Vt;
    const addVt = (rs, val) => {
      if (readout) this._j(b, rs, b, 1 + t, (val * VT) / z[t]);
      else if (floating) this._j(b, rs, tb, 0, val * VT);
    };
    const toTerminal = floating && !readout; // accumulate the terminal-current row
    const termJ = this._termIJac;

    // Before anything is added, each balance residual is the flux through this face.
    for (let i = 0; i < n; i++) flux[i] = sgn * res[b * M + 1 + i];

    // Electrode reactions: Σ ν_R R + n e⁻(metal, μ̄ = −F V_t) ⇌ Σ ν_P P.
    for (const rx of ct.reactions) {
      const al = rx.alpha;
      let pref = rx.k0, zs = 0, aHi = rx.fixedA - (rx.electrons * Vt) / VT, aLo = 0;
      for (const { i, nu } of rx.reactants) {
        pref *= (c[g * n + i] / this.cRef[g * n + i]) ** (nu * (1 - al));
        zs += nu * (1 - al) * z[i];
        aHi += nu * u[b * M + 1 + i];
        aLo += nu * uLo[b * M + 1 + i];
      }
      for (const { i, nu } of rx.products) {
        pref *= (c[g * n + i] / this.cRef[g * n + i]) ** (nu * al);
        zs += nu * al * z[i];
        aHi -= nu * u[b * M + 1 + i];
        aLo -= nu * uLo[b * M + 1 + i];
      }
      const { g: gf, gp } = bvFactor(aHi + aLo, al);
      const rate = pref * gf;
      const dVt = (pref * gp * -rx.electrons) / VT;
      // ∂rate/∂η_j for participant j: prefactor exponent·rate ± pref·g′·ν
      const dEta = (j, e, sign, nu) => e * rate + sign * pref * gp * nu;
      const terms = [
        ...rx.reactants.map(({ i, nu }) => ({ i, nu, row: 1, d: dEta(i, nu * (1 - al), 1, nu) })),
        ...rx.products.map(({ i, nu }) => ({ i, nu, row: -1, d: dEta(i, nu * al, -1, nu) })),
      ];
      for (const tA of terms) {
        const rs = 1 + tA.i, w = tA.row * tA.nu; // consumed (+) or produced (−) at this node
        res[b * M + rs] += w * rate;
        for (const tB of terms) this._j(b, rs, b, 1 + tB.i, w * tB.d);
        this._j(b, rs, b, 0, w * -zs * rate);
        addVt(rs, w * dVt);
      }
      if (toTerminal) {
        // Electrons taken from the metal at the right: current toward +x of n F r.
        const q = rx.electrons * F;
        this._termI += q * rate;
        for (const tB of terms) termJ[1 + tB.i] += q * tB.d;
        termJ[0] += q * -zs * rate;
        termJ[M] += q * dVt;
      }
    }

    // Conductance links: J (toward the device) = G (V_out − V_i), V_out = V_t + offset.
    for (let i = 0; i < n; i++) {
      const link = ct.species[i];
      if (link.type !== 'conductance') continue;
      const o = b * M + 1 + i;
      const Vi = (VT * (u[o] + uLo[o])) / z[i];
      const Nin = (link.G * (Vt + link.offset - Vi)) / (z[i] * F); // particles entering
      const k = (link.G * VT) / (z[i] * z[i] * F);
      res[o] -= Nin;
      this._j(b, 1 + i, b, 1 + i, k);
      addVt(1 + i, -link.G / (z[i] * F));
      if (toTerminal) {
        // Current toward +x leaving through the right face: −z F N_in.
        this._termI += -z[i] * F * Nin;
        termJ[1 + i] += link.G * VT / z[i];
        termJ[M] += -link.G;
      }
    }

    // Fixed links.
    for (let i = 0; i < n; i++) {
      if (ct.species[i].type !== 'fixed') continue;
      const o = b * M + 1 + i;
      this._replaceRow(b, 1 + i);
      if (readout && i === t) {
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
        addVt(1 + i, -z[i] / VT);
      } else {
        this._j(b, 1 + i, b, 1 + i, 1);
        res[o] = u[o] - this.contactEta(side, i) + uLo[o];
      }
    }

    // φ link.
    const link = ct.phi;
    if (link.type === 'capacitive') {
      // Gate, or metal across a Stern layer, at φ_g = V_t − zeroCharge. D toward +x.
      const phiG = Vt - link.zeroCharge;
      const D = link.C * (phiG - VT * u[b * M]) * sgn;
      res[b * M] -= sgn * D;
      this._j(b, 0, b, 0, link.C * VT);
      addVt(0, -sgn * sgn * link.C); // ∂(−sgn·D)/∂V_t = −C
      this.contactD[side] = D;
      if (toTerminal && Number.isFinite(dt)) {
        this._termI += (D - this.contactDStart.right) / dt;
        termJ[0] += (link.C * VT) / dt;
        termJ[M] += -link.C / dt;
      }
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
    let mx = this.terminalBlock >= 0 ? Math.abs(delta[this.terminalBlock * M]) : 0;
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
    this.contactDStart = DOld;
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
    if (!Number.isFinite(dt)) return;
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
   * Steady state.
   * - If every species stretch is fed by a contact, nothing is conserved on its own, so the
   *   true steady equations (dt = ∞, no storage term) are solved directly: no slow modes to
   *   wait out, however slow the physics (e.g. exponentially scarce minority carriers).
   * - Otherwise backward-Euler steps at a huge dt, whose storage term pins each conserved
   *   amount exactly; dt keeps growing ×10 while the state still moves.
   * If Newton fails, dt ramps up from a small value (pseudo-transient continuation) instead.
   * The clock is not advanced, and open-system conservation bookkeeping restarts here.
   */
  solveSteady({ maxSteps = 80, tol = 1e-11 } = {}) {
    const time = this.time;
    const tau = this.slowestTime();
    const direct = this.stretches.every((st) => st.connected);
    const giant = direct ? Infinity : 1e6 * tau;
    let dt = giant;
    let totalIter = 0, steps = 0, converged = false;
    const history = [];
    while (steps < maxSteps) {
      steps++;
      // Restore conserved amounts exactly before each huge step; the step then re-solves, so
      // the final state satisfies every equation.
      if (!direct && dt >= giant) this._renormalizeSpectators();
      const r = this.step(dt);
      totalIter += r.iterations;
      history.push({ dt, converged: r.converged, iterations: r.iterations });
      if (!r.converged) {
        dt = dt >= giant ? tau * 1e-6 : dt / 4;
        if (dt < tau * 1e-15) break;
        continue;
      }
      if (dt === Infinity) {
        converged = true; // the steady equations themselves were solved
        break;
      }
      if (dt >= giant && this._maxPotentialStep(this._diff()) < tol) {
        converged = true;
        if (dt > giant) {
          // Finish at the base giant step, where the storage term pins amounts tightly.
          this._renormalizeSpectators();
          const f = this.step(giant);
          totalIter += f.iterations;
          converged = f.converged;
        }
        break;
      }
      if (dt < giant) {
        dt *= 10; // ramping up after a failure
        if (dt > 1e6 * tau) dt = giant; // then the direct steady solve (or the giant step)
      } else if (dt < 1e6 * giant) {
        dt *= 10; // conserved amounts present and still moving: let dt keep growing (capped)
      }
    }
    this.time = time;
    if (converged) {
      this.stretches.forEach((st, k) => {
        if (!st.connected) return;
        this.referenceAmounts[k] = this.amount(st);
        this.boundaryIntake[k] = 0;
      });
    }
    return { converged, steps, iterations: totalIter, history };
  }

  // Huge steps pin each conserved amount only through a tiny storage term, so round-off can let
  // it creep. Shift each spectator's level uniformly to restore its amount exactly (exact for
  // ideal statistics, where c ∝ e^η).
  _renormalizeSpectators() {
    const { n, M, u } = this;
    let changed = false;
    this.stretches.forEach((st, k) => {
      if (!st.spectator) return;
      const now = this.amount(st), want = this.referenceAmounts[k];
      if (!(now > 0) || now === want) return;
      const shift = Math.log(want / now);
      for (let g = st.nodes[0]; g <= st.nodes[1]; g++) u[this.blockOfNode[g] * M + 1 + st.species] += shift;
      changed = true;
    });
    if (changed) this.computeConcentrations();
  }

  _diff() {
    const d = this.delta; // reuse as scratch
    for (let k = 0; k < d.length; k++) d[k] = this.u[k] - this.uPrev[k];
    return d;
  }
}
