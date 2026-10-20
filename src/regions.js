// The equations inside the regions: each node's statistics (c from the levels), its storage and
// space charge, the fluxes along each segment, and the bulk reactions; and the strictly neutral
// nodes' rows in better-conditioned terms.
//
// CONTRACT
//
// Methods of the Solver, called by assemble() (src/assemble.js, whose contract sets the rows'
// conventions: a balance row is storage + flux out − flux in + consumed − made, the Gauss row
// Σ D out − vol·ρ, totals over the box) and its bookkeeping. They read the state u + uLo (φ̂, η),
// the per-node data the Solver's constructor set (cRef, mu0hat, rhoFixed, present,
// nodeMaterial, nodeIdeal, nodeConductor, sheetFace, phiUndefined, chargeNode, lateStorage, the
// bulk reactions rxs), cOld, and write res and sys through _j or by local index (the dilute fast
// paths), or their terms into this.lin where it's set.
//
// Statistics, computeConcentrations(): at every node, from ζ_i = η_i − μ°_i/RT − z_i φ̂ (η
// compensated, u + uLo; φ̂ as u alone):
//   ideal: c_i = cRef_i e^ζ_i, and K = diag(c) (not stored);
//   non-ideal: each material's statistics models (src/statistics.js) give c and K = ∂c/∂ζ
//     (n×n at each node, symmetric), and ex_i = ζ_i − ln(c_i/cRef_i) the excess chemical
//     potential; a species no model holds is ideal there;
//   a conductor's node: no bulk concentration; its carrier's c is the excess charge of a charged
//     face's sheet, sheetSign·D·A/(zF·vol);
//   absent species: c = 0; then the coverages (_coverages: Langmuir, θ_s = e^{ζ_s}/(1 + Σe^ζ)).
//   ∂c_i/∂η_j = K_ij, ∂c_i/∂φ̂ = −(Kz)_i (_Kz; _dlnc adds these as ∂ln c, _dc sets them as ∂c).
//   At a non-ideal node of a region with eddy mixing (nodeDK), K's slope dK = ∂K/∂ζ as well.
//
// The terms:
//   node (dilute, concentrated): storage vol·(c − c_old)/dt in each present species' row, and
//     −vol·ρ in the Gauss row (none where φ is undefined, or, where the node's charge rows will
//     carry them, at a chargeNode on a time step: there only a neutral species' storage); at a
//     lateStorage node on a time step the ions' storage is left to the contact (_contact); a
//     species a port holds at a device end node keeps its storage here.
//   node (conductor): its carrier's storage (a charged face's sheet, through that face's D).
//   displacement: D = −ε V_T (φ̂_R − φ̂_L)/h toward +x, out of the left box and into the right;
//     none where ε = 0 or φ is undefined at the segment's left node.
//   Scharfetter–Gummel (dilute): N = −(D_i/h) B(Δ) c_L expm1(Δη), Δ = z Δφ̂ − Pe, Δη the
//     compensated η_R − η_L − Pe, Pe = v h/D_i (advection at the region's velocity; flow is
//     planar only, so h is the segment's length): exactly zero where μ̄ is level and v = 0, and
//     exact in the cancellation of drift against diffusion. Past a drift Δ of SG_FAR (300), the
//     same flux from the downstream end, g B(−Δ) c_R expm1(−Δη) (sgFlux in src/bernoulli.js):
//     B(Δ) and expm1(Δη) underflow and overflow there, toward −x at a large Péclet number.
//   concentrated: the same with Δ = zΔφ̂ + Δex − Pe, and c_L, ex through K at both ends.
//   conductor: the flux J in slot 0 of the segment's left node, its row η_R − η_L + J/g = 0,
//     g = σRT/(z²F²h) (Ohm's law), and J out of the left node's carrier balance, into the right's.
//   eddy mixing: N = −(D_mix/h) P̄ Δη, P̄ = K̄ − (K̄z)(K̄z)ᵀ/(zᵀK̄z) over the mobile species: no
//     current, and zero where μ̄ is level. Between ideal ends K̄ = C̄, the logarithmic means of
//     the two ends' c, and a neutral species' mixing is exactly −D_mix Δc/h; elsewhere
//     (_mixCrowded) K̄ = ½(R_L K_L R_L + R_R K_R R_R), R = diag(√(c̄/c)), and in a neutral
//     solution every species' mixing is −D_mix Δc/h to second order (eddies mix concentration).
//   bulk reaction at a node: rate r = P(1 − e^{−a}), a = fixedA + Σ ν η (reactants ν > 0,
//     products ν < 0; compensated), P = k_f Π c^ν over reactants (times generationScale for a
//     generation); or SRH (_srhRate) where the material has an SRH law, its derivatives through
//     ln c and the excess potentials ex = ζ − ln c (_srhInto: exact for any statistics, the
//     carriers' coupled or not). vol·ν·r in each
//     participant's row. bulkRate, bulkForward and bulkOneWay read the same rate, its forward
//     one-way part, and the larger one-way part, at the present state.
//   charge rows (_chargeRows, after the laws, on a time step where combining): at each
//     chargeNode, the balance of the charged species with the largest z²c becomes Σ (z_i/z_k)
//     × the charged balances (current continuity, without storage); the columns become φ̂' and
//     η'_i = η_i − z_i φ̂' (in every block row and termC that reads them; _untransform maps an
//     update back); then the other charged species' storage and the neutrality row, Σ z F c +
//     ρ_fixed = 0 (with a port's capacitance there, its charge, and the ions' storage that
//     follows from it in the combined row).
//
// Invariants (each says what a reviewer should be able to check here):
//   1. Every term's Jacobian entries are its residual's derivatives (via c(ζ) and K); in
//      difference form (this.lin), the same terms, each on differences of δη and δφ̂. One
//      exception: crowded mixing's K̄ needs K's own slope, which no statistics model gives, so it
//      comes from central differences of the model (_statSlopes: near 1e-10 of K).
//   2. A segment's flux goes into its two nodes' rows with opposite signs, the same number.
//   3. At a level μ̄ along a segment (Δη = 0, no velocity), its flux and mixing are exactly 0;
//      at a = 0 a bulk reaction's rate is exactly 0.
//   4. computeConcentrations writes c, K, ex, th, th0, dK at the nodeDK nodes (and scratch:
//      zeta, the statistics models' work arrays) and nothing else; the kernels write res,
//      sys (or lin), termB/termC (the charge rows), and scratch kept on the Solver (dA, dB, jL,
//      jR, sg, mixWork, mixCrowd) or cached on a material (mobile).
//   5. The charge rows are an invertible transformation of the equations at the node (the same
//      solution), and _untransform inverts the column change exactly.
// Tests: test/jacobian.test.js (1), test/layout.test.js (difference form), test/srh.test.js,
// test/statistics.test.js, the validation suite (analytic profiles: test/equilibrium.test.js,
// test/junction.test.js, test/neutral.test.js and others).

import { sgFlux } from './bernoulli.js';
import { EPS0, FARADAY } from './constants.js';
import { powi } from './pow.js';

// The logarithmic mean of a and b = a·e^l, a·f(l) with f = expm1(l)/l, into cb[k], with
// ∂/∂a = f − f′ and ∂/∂b = a f′/b into dLa[k] and dLb[k].
function logMean(a, b, l, cb, dLa, dLb, k) {
  let f, fp;
  if (Math.abs(l) < 1e-3) {
    f = 1 + l / 2 + (l * l) / 6 + (l * l * l) / 24;
    fp = 0.5 + l / 3 + (l * l) / 8 + (l * l * l) / 30;
  } else {
    const em = Math.expm1(l);
    f = em / l;
    fp = (l * (em + 1) - em) / (l * l);
  }
  cb[k] = a * f;
  dLa[k] = f - fp;
  dLb[k] = b > 0 ? (a * fp) / b : 0.5;
}

// _statSlopes' step in ζ.
const STAT_STEP = 1e-5;

// _mixCrowded's work arrays.
function mixWork(n, M) {
  const W = {};
  for (const k of ['a', 'b', 'cb', 'dLa', 'dLb', 'rL', 'rR', 'de', 'zc', 'g', 'kz', 'N', 'tg', 'tkz', 'tr', 'tcbr']) W[k] = new Float64Array(n);
  for (const k of ['KL', 'KR', 'Kb', 'tK']) W[k] = new Float64Array(n * n);
  W.Jl = new Float64Array(n * (n + 1));
  W.Jr = new Float64Array(n * (n + 1));
  W.dcs = new Float64Array(n * M);
  return W;
}

export class RegionTerms {
  // Each electrode surface's coverages from its η: Langmuir on shared sites, θ_s = e^{ζ_s}/(1 + Σ e^{ζ}),
  // ζ = η − μ°/RT (so θ_s/θ₀ = e^{ζ_s}), at every node of a surfaced window.
  _coverages() {
    const { n, M, u, uLo, nSurf, th } = this, RT = this.model.RT;
    for (let g = 0; g < this.nNodes; g++) {
      const k = this.surfPort[g];
      if (k < 0) continue;
      const surf = this.model.ports[k].surface, b = this.blockOfNode[g];
      let m = 0;
      for (let q = 0; q < surf.length; q++) m = Math.max(m, u[b * M + 1 + n + q] + uLo[b * M + 1 + n + q] - surf[q].mu0 / RT);
      let sum = Math.exp(-m);
      for (let q = 0; q < surf.length; q++) sum += (th[g * nSurf + q] = Math.exp(u[b * M + 1 + n + q] + uLo[b * M + 1 + n + q] - surf[q].mu0 / RT - m));
      for (let q = 0; q < surf.length; q++) th[g * nSurf + q] /= sum;
      this.th0[g] = Math.exp(-m) / sum;
    }
  }

  computeConcentrations() {
    const { n, M, u, uLo, c, z } = this;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      const phiHat = u[b * M];
      if (!this.nodeIdeal[g]) {
        this._nodeStatistics(g, b, phiHat);
        continue;
      }
      const im = this.nodeConductor[g];
      if (im >= 0) {
        // Excess carriers: zero in the bulk, the surface sheet at a charged face.
        for (let i = 0; i < n; i++) c[g * n + i] = 0;
        const f = this.sheetFace[g];
        if (f >= 0) c[g * n + im] = (this.sheetSign[g] * u[this.blockOfFace[f] * M] * this.model.grid.area[g]) / (z[im] * FARADAY * this.model.grid.vol[g]);
        continue;
      }
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        c[k] = this.present[k] ? this.cRef[k] * Math.exp(u[b * M + 1 + i] + uLo[b * M + 1 + i] - this.mu0hat[k] - z[i] * phiHat) : 0;
      }
    }
    this._coverages();
  }

  _work(md) {
    let w = this.scratch.get(md);
    if (!w) {
      const k = md.idx.length;
      w = { z: new Float64Array(k), c: new Float64Array(k), K: new Float64Array(k * k), fixed: new Uint8Array(k), t: new Float64Array(k) };
      this.scratch.set(md, w);
    }
    return w;
  }

  // c, K and the excess at a node of a non-ideal material.
  _nodeStatistics(g, b, phiHat) {
    const { n, M, u, uLo, c, z, K, ex, zeta } = this;
    const mat = this.model.materials[this.nodeMaterial[g]];
    const Kg = g * n * n;
    K.fill(0, Kg, Kg + n * n);
    for (let i = 0; i < n; i++) {
      const k = g * n + i;
      ex[k] = 0;
      if (!this.present[k]) {
        c[k] = 0;
        continue;
      }
      zeta[i] = u[b * M + 1 + i] + uLo[b * M + 1 + i] - this.mu0hat[k] - z[i] * phiHat;
      if (mat.modelOf[i] < 0) {
        c[k] = this.cRef[k] * Math.exp(zeta[i]);
        K[Kg + i * n + i] = c[k];
      }
    }
    for (const md of mat.models) {
      const w = this._work(md), idx = md.idx, k = idx.length;
      for (let a = 0; a < k; a++) w.z[a] = zeta[idx[a]];
      md.evaluate(w.z, w.c, w.K, this.background[g]);
      for (let a = 0; a < k; a++) {
        const i = idx[a], ci = w.c[a];
        c[g * n + i] = ci;
        ex[g * n + i] = ci > 0 ? zeta[i] - Math.log(ci / this.cRef[g * n + i]) : 0;
        for (let q = 0; q < k; q++) K[Kg + i * n + idx[q]] = w.K[a * k + q];
      }
    }
    if (this.nodeDK?.[g]) this._statSlopes(g, mat);
  }

  // K's own slope at node g, dK[g·n³ + (i·n + j)·n + r] = ∂K_ij/∂ζ_r, which no statistics model
  // gives: by central differences of each model's evaluate, a step of 1e-5 in ζ (truncation
  // and round-off near 1e-10 of K). Only where crowded eddy mixing needs it (nodeDK).
  _statSlopes(g, mat) {
    const { n, c, dK, zeta } = this;
    const o = g * n * n * n;
    dK.fill(0, o, o + n * n * n);
    for (let i = 0; i < n; i++) if (this.present[g * n + i] && mat.modelOf[i] < 0) dK[o + (i * n + i) * n + i] = c[g * n + i];
    for (const md of mat.models) {
      const w = this._work(md), idx = md.idx, k = idx.length;
      w.zp ??= new Float64Array(k);
      w.cp ??= new Float64Array(k);
      w.Kp ??= new Float64Array(k * k);
      w.Km ??= new Float64Array(k * k);
      for (let r = 0; r < k; r++) {
        for (const [sgn, Kx] of [[1, w.Kp], [-1, w.Km]]) {
          for (let a = 0; a < k; a++) w.zp[a] = zeta[idx[a]];
          w.zp[r] += sgn * STAT_STEP;
          md.evaluate(w.zp, w.cp, Kx, this.background[g]);
        }
        for (let a = 0; a < k; a++) {
          for (let q = 0; q < k; q++) dK[o + (idx[a] * n + idx[q]) * n + idx[r]] = (w.Kp[a * k + q] - w.Km[a * k + q]) / (2 * STAT_STEP);
        }
      }
    }
  }

  // Σ_j K_ij z_j at node g (= z_i c_i for ideal statistics).
  _Kz(g, i) {
    const { n, z } = this;
    if (this.nodeIdeal[g]) return z[i] * this.c[g * n + i];
    let s = 0;
    const o = g * n * n + i * n;
    for (let j = 0; j < n; j++) s += this.K[o + j] * z[j];
    return s;
  }

  /** zᵀKz at node g (per RT): the charge response that sets the screening length. */
  screening(g) {
    let s = 0;
    for (let i = 0; i < this.n; i++) s += this.z[i] * this._Kz(g, i);
    return s;
  }

  // Add w·∂(ln c_i)/∂(slot) at node g into d (slot 0 = φ̂, 1 + j = η_j).
  _dlnc(g, i, w, d) {
    const { n, z } = this;
    if (this.nodeIdeal[g]) {
      d[1 + i] += w;
      d[0] -= w * z[i];
      return;
    }
    const ci = this.c[g * n + i];
    if (!(ci > 0)) return;
    const o = g * n * n + i * n, f = w / ci;
    for (let j = 0; j < n; j++) d[1 + j] += f * this.K[o + j];
    d[0] -= f * this._Kz(g, i);
  }

  // A conductor node: only its carrier's balance (the segment flux J has its own row, in
  // _segmentConductor). Its bulk stores nothing; a charged face's sheet sits in the edge node.
  _nodeConductor(g, dt) {
    const { n, M, res, c, cOld, z } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g], im = this.nodeConductor[g];
    const r = 1 + im, k = g * n + im;
    res[R[b * M + r]] += (v * (c[k] - cOld[k])) / dt;
    const f = this.sheetFace[g];
    if (f >= 0) this._j(b, r, this.blockOfFace[f], 0, (this.sheetSign[g] * this.model.grid.area[g]) / (z[im] * FARADAY * dt));
  }

  // Whether a port holds species i at node g, a device end node (see portHeldEnd).
  _portHeldAt(g, i) {
    return g === 0 ? this.portHeldEnd[i] === 1 : g === this.nNodes - 1 && this.portHeldEnd[this.n + i] === 1;
  }

  // Dilute nodes g0…g1 (ideal statistics, c = c_ref e^ζ): storage of each species and the space
  // charge in the Gauss row, written straight into the diagonal blocks by local index. A
  // dielectric is the case with no species.
  _nodesDilute(g0, g1, dt) {
    const { n, M, res, c, cOld, z, sys, loc, present, rhoFixed, blockOfNode, lin } = this, R = this.rix, F = FARADAY;
    const JB = sys.B, sizes = sys.sizes, offB = sys.offB, vol = this.model.grid.vol;
    for (let g = g0; g <= g1; g++) {
      const b = blockOfNode[g], v = vol[g];
      const m = sizes[b], oB = offB[b], lb = b * M, p = loc[lb]; // p: φ's row, −1 where undefined
      if (this.combining && this.chargeNode[g] === 1 && dt !== Infinity) {
        // Storage and neutrality go in later, in _chargeRows; a neutral species' storage here.
        for (let i = 0; i < n; i++) if (present[g * n + i] && z[i] === 0) this._storage(g, i, dt);
        continue;
      }
      let q = rhoFixed[g], dq = 0;
      const late = this.combining && this.lateStorage[g] === 1 && dt !== Infinity; // (the ions' storage: in _contact)
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        if (!present[k]) continue;
        const ck = c[k], l = loc[lb + 1 + i], stores = !(late && z[i] !== 0 && !this._portHeldAt(g, i));
        if (stores) res[R[lb + 1 + i]] += (v * (ck - cOld[k])) / dt;
        q += F * z[i] * ck;
        dq += F * z[i] * z[i] * ck;
        if (lin) {
          lin.node(R[lb + 1 + i], R[lb], z[i], stores ? v * ck : 0, -v * F * z[i] * ck);
          continue;
        }
        if (stores) {
          JB[oB + l * m + l] += (v * ck) / dt;
          if (p >= 0 && z[i] !== 0) JB[oB + l * m + p] += (-v * z[i] * ck) / dt;
        }
        if (p >= 0 && z[i] !== 0) JB[oB + p * m + l] += -v * F * z[i] * ck;
      }
      if (p >= 0) {
        res[R[lb]] -= v * q;
        if (!lin) JB[oB + p * m + p] += v * dq;
      }
    }
  }

  // A concentrated node (any statistics): ∂c_i/∂η_j = K_ij and ∂c_i/∂φ̂ = −(Kz)_i.
  _nodeConcentrated(g, dt) {
    const { n, M, res, c, z } = this, R = this.rix, F = FARADAY;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g];
    if (this.combining && this.chargeNode[g] === 1 && dt !== Infinity) {
      // Storage and neutrality go in later, in _chargeRows; a neutral species' storage here.
      for (let i = 0; i < n; i++) if (this.present[g * n + i] && z[i] === 0) this._storage(g, i, dt);
      return;
    }
    let q = this.rhoFixed[g], dq = 0;
    const lin = this.lin;
    for (let i = 0; i < n; i++) {
      const k = g * n + i, r = 1 + i;
      if (!this.present[k]) continue;
      const Kz = this._Kz(g, i);
      if (lin) {
        // In difference form (see DifferenceTerms): storage of i and its charge, in each
        // δη_j − z_j δφ̂ through K_ij.
        res[R[b * M + r]] += (v * (c[k] - this.cOld[k])) / dt;
        q += F * z[i] * c[k];
        const Kg = g * n * n + i * n;
        for (let j = 0; j < n; j++) {
          const Kij = this.nodeIdeal[g] ? (i === j ? c[k] : 0) : this.K[Kg + j];
          if (Kij !== 0) lin.node(R[b * M + 1 + j], R[b * M], z[j], v * Kij, -v * F * z[i] * Kij, R[b * M + r]);
        }
        continue;
      }
      this._storage(g, i, dt, Kz);
      q += F * z[i] * c[k];
      dq += F * z[i] * Kz;
      if (Kz !== 0) this._j(b, 0, b, r, -v * F * Kz); // K symmetric: ∂(Σ z c)/∂η_i = (Kz)_i
    }
    // (where no charge responds and no field reaches, φ isn't defined, nor an unknown)
    if (!this.phiUndefined[g]) {
      res[R[b * M]] -= v * q;
      if (!lin) this._j(b, 0, b, 0, v * dq);
    }
  }

  // Species i's storage at node g, v(c − c_old)/dt, in its balance row and the Jacobian
  // (∂c_i/∂η_j = K_ij, ∂c_i/∂φ̂ = −(Kz)_i; for ideal statistics K is diagonal, c_i).
  _storage(g, i, dt, Kz = this._Kz(g, i)) {
    const { n, M, res, c, cOld } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g], k = g * n + i, r = 1 + i;
    res[R[b * M + r]] += (v * (c[k] - cOld[k])) / dt;
    if (this.nodeIdeal[g]) this._j(b, r, b, r, (v * c[k]) / dt);
    else {
      const Kg = g * n * n;
      for (let j = 0; j < n; j++) {
        const Kij = this.K[Kg + i * n + j];
        if (Kij !== 0) this._j(b, r, b, 1 + j, (v * Kij) / dt);
      }
    }
    if (Kz !== 0 && !this.phiUndefined[g]) this._j(b, r, b, 0, (-v * Kz) / dt);
  }

  // Strictly neutral nodes on a transient step, solved in better-conditioned terms. There, a
  // change of φ̂ with every η_i shifted by z_i times it leaves every concentration as it was, so
  // storage and neutrality don't see it; but in (φ̂, η) they see it as pairs of huge entries that
  // cancel, and round-off in that cancellation (storage/flux ~ h²/(D·dt), times the range of
  // concentrations: a trace ion beside 3 M KCl) swamps the fluxes that do fix it. So, at each node
  // that stays neutral (see chargeNode), with storage and neutrality left out of assembly:
  // - rows: the balance of the most abundant charged species (by z²c) becomes Σ (z_i/z_k) × each
  //   balance, which without storage is current continuity;
  // - columns: the unknowns become φ̂' and η'_i = η_i − z_i φ̂ (newton() maps the update back),
  //   in which storage and neutrality have no φ̂' term at all;
  // - then storage (on the other balances) and neutrality go in, exactly, in those terms.
  // A start-of-step net charge (round-off in a solved state) isn't carried over: the neutrality
  // row holds the new state neutral.
  _chargeRows(dt) {
    const { n, M, res, c, cOld, z, sys, loc, present, termB, termC } = this, R = this.rix, F = FARADAY, live = this._liveTerms();
    const { A, B, C, sizes, offA, offB, offC } = sys, vol = this.model.grid.vol, last = this.nB - 1;
    for (const g of this.chargeNodes) {
      const b = this.blockOfNode[g], lb = b * M, m = sizes[b], mp = b > 0 ? sizes[b - 1] : 0, mn = b < last ? sizes[b + 1] : 0;
      const p = loc[lb], gn = g * n;
      let k = -1, best = -1;
      for (let i = 0; i < n; i++) {
        if (!present[gn + i] || z[i] === 0) continue;
        const w = z[i] * z[i] * c[gn + i];
        if (w > best) {
          k = i;
          best = w;
        }
      }
      if (k < 0 || p < 0) continue;
      // Rows: the pivot's balance becomes the z-weighted sum (all still without storage).
      const lk = loc[lb + 1 + k], rk = R[lb + 1 + k];
      for (let i = 0; i < n; i++) {
        if (i === k || !present[gn + i] || z[i] === 0) continue;
        const li = loc[lb + 1 + i], w = z[i] / z[k], ri = R[lb + 1 + i];
        for (let q = 0, o = offA[b] + lk * mp, s = offA[b] + li * mp; q < mp; q++) A[o + q] += w * A[s + q];
        for (let q = 0, o = offB[b] + lk * m, s = offB[b] + li * m; q < m; q++) B[o + q] += w * B[s + q];
        for (let q = 0, o = offC[b] + lk * mn, s = offC[b] + li * mn; q < mn; q++) C[o + q] += w * C[s + q];
        res[rk] += w * res[ri];
        for (const t of live) termB[t][rk] += w * termB[t][ri];
        if (this.lin) this.lin.combine(rk, ri, w);
      }
      // Columns: φ̂' = φ̂ with η fixed becomes φ̂ with η_i shifted by z_i, wherever node g's
      // unknowns appear (its own rows, its neighbours', a terminal's current).
      for (let i = 0; i < n; i++) {
        if (!present[gn + i] || z[i] === 0) continue;
        const li = loc[lb + 1 + i], zi = z[i];
        for (let r = 0, o = offB[b]; r < m; r++) B[o + r * m + p] += zi * B[o + r * m + li];
        if (b > 0) for (let r = 0, mr = mp, o = offC[b - 1]; r < mr; r++) C[o + r * m + p] += zi * C[o + r * m + li];
        if (b < last) for (let r = 0, mr = mn, o = offA[b + 1]; r < mr; r++) A[o + r * m + p] += zi * A[o + r * m + li];
        for (const t of live) termC[t][R[lb]] += zi * termC[t][R[lb + 1 + i]];
      }
      // Storage on the other charged balances, and neutrality: no φ̂' terms, exactly.
      const v = vol[g], oB = offB[b];
      let q = this.rhoFixed[g];
      if (this.nodeIdeal[g]) {
        for (let i = 0; i < n; i++) {
          if (!present[gn + i] || z[i] === 0) continue;
          const ci = c[gn + i], li = loc[lb + 1 + i];
          if (i !== k) {
            res[R[lb + 1 + i]] += (v * (ci - cOld[gn + i])) / dt;
            B[oB + li * m + li] += (v * ci) / dt;
          }
          q += F * z[i] * ci;
          B[oB + p * m + li] += -v * F * z[i] * ci;
        }
      } else {
        for (let i = 0; i < n; i++) {
          if (!present[gn + i]) continue;
          if (i !== k && z[i] !== 0) this._storage(g, i, dt, 0);
          q += F * z[i] * c[gn + i];
          const Kz = this._Kz(g, i); // (nonzero for a neutral species only on a shared lattice)
          if (Kz !== 0) this._j(b, 0, b, 1 + i, -v * F * Kz);
        }
      }
      res[R[lb]] -= v * q;
      // The pivot's row (charge continuity) left storage out, the ions' net charge being constant
      // under neutrality. Against a capacitance it isn't: neutrality makes it −(aσ + ρ_fixed)/F,
      // so their storage is that less what the step started with, Σ z c_old, over dt (in units
      // of the pivot's balance).
      if (this.qgTerms[g].length > 0) {
        const f = -v / (FARADAY * z[k] * dt);
        let qOld = this.rhoFixed[g];
        for (let i = 0; i < n; i++) if (present[gn + i]) qOld += F * z[i] * cOld[gn + i];
        res[rk] += f * (this.qg[g] + qOld);
        B[offB[b] + lk * m + p] += -f * this.qgSlope[g] * this.VT;
        for (const [kq, w] of this.qgTerms[g]) termB[kq][rk] += f * this.portArea[kq - 2][w] * this.model.ports[kq - 2].capacitance.C;
      }
    }
  }

  // Newton's update in the transformed unknowns (see _chargeRows) back in η: η_i = η'_i + z_i φ̂'.
  _untransform(delta) {
    const { n, M, z, loc, present } = this, R = this.rix;
    for (const g of this.chargeNodes) {
      const lb = this.blockOfNode[g] * M;
      if (loc[lb] < 0) continue;
      const dphi = delta[R[lb]];
      for (let i = 0; i < n; i++) if (present[g * n + i] && z[i] !== 0) delta[R[lb + 1 + i]] += z[i] * dphi;
    }
  }

  // Bulk reactions at node g: r = k_f Π c_R^ν · (−expm1(−a)), with a = A/RT from the
  // (compensated) η. Reactants are consumed (+v·ν·r in their balance), products made (−v·ν·r).
  _bulkReactions(g, rxs) {
    const { n, M, u, uLo, res, c } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g];
    const dr = this.dA;
    for (let x = 0; x < rxs.length; x++) {
      const rx = rxs[x];
      const law = rx.srh?.[this.nodeMaterial[g]];
      if (law) {
        // SRH: dr is ∂rate/∂slot directly, from its sensitivities to ln n, ln p, a and the excesses.
        let aHi = rx.fixedA, aLo = 0;
        for (let p = 0; p < rx.sp.length; p++) {
          aHi += rx.nu[p] * u[b * M + 1 + rx.sp[p]];
          aLo += rx.nu[p] * uLo[b * M + 1 + rx.sp[p]];
        }
        const { rate, Cn, Cp, Ca, Cxn, Cxp } = this._srhRate(law, c[g * n + law.n.i], c[g * n + law.p.i], aHi + aLo, ...this._srhDegen(g, g, law));
        dr.fill(0);
        this._srhInto(g, law.n.i, Cn, Cxn, dr);
        this._srhInto(g, law.p.i, Cp, Cxp, dr);
        for (let p = 0; p < rx.sp.length; p++) dr[1 + rx.sp[p]] += Ca * rx.nu[p];
        for (let p = 0; p < rx.sp.length; p++) {
          const row = 1 + rx.sp[p], w = v * rx.nu[p];
          res[R[b * M + row]] += w * rate;
          for (let k = 0; k < M; k++) if (dr[k] !== 0) this._j(b, row, b, k, w * dr[k]);
        }
        continue;
      }
      let P = rx.generation ? rx.kfNode[g] * this.generationScale : rx.kfNode[g], aHi = rx.fixedA, aLo = 0;
      dr.fill(0); // ∂ ln P / ∂slot
      // Participants: reactants with ν > 0, then products with ν < 0.
      for (let p = 0; p < rx.sp.length; p++) {
        const i = rx.sp[p], nu = rx.nu[p];
        if (nu > 0) {
          const ci = c[g * n + i];
          P *= nu === 1 ? ci : powi(ci, nu);
          this._dlnc(g, i, nu, dr);
        }
        aHi += nu * u[b * M + 1 + i];
        aLo += nu * uLo[b * M + 1 + i];
      }
      const a = aHi + aLo;
      const f = -Math.expm1(-a); // 1 − e^{−a}
      const rate = P * f;
      const Pd = P * (1 - f); // P·df/da, df/da = e^{−a}
      // ∂rate/∂slot = rate·∂lnP + P·f′·∂a
      for (let k = 0; k < M; k++) dr[k] *= rate;
      for (let p = 0; p < rx.sp.length; p++) dr[1 + rx.sp[p]] += Pd * rx.nu[p];
      for (let p = 0; p < rx.sp.length; p++) {
        const row = 1 + rx.sp[p], w = v * rx.nu[p];
        res[R[b * M + row]] += w * rate;
        for (let k = 0; k < M; k++) if (dr[k] !== 0) this._j(b, row, b, k, w * dr[k]);
      }
    }
  }

  // Displacement along a segment where φ is defined: D = −ε(φ_R − φ_L)/h into the two Gauss rows.
  _segmentDisplacement(s, bL, bR, mat, h) {
    const { M, u, res, sys, loc } = this, R = this.rix;
    const k = this.phiUndefined[s] ? 0 : (mat.epsr * EPS0 * this.VT) / h; // ε = 0: no displacement
    if (k === 0) return;
    const D = -k * (u[bR * M] - u[bL * M]);
    const mL = sys.sizes[bL], mR = sys.sizes[bR], pL = loc[bL * M], pR = loc[bR * M];
    res[R[bL * M]] += D;
    res[R[bR * M]] -= D;
    if (this.lin) return this.lin.displacement(R[bL * M], R[bR * M], k);
    sys.B[sys.offB[bL] + pL * mL + pL] += k;
    sys.C[sys.offC[bL] + pL * mR + pR] -= k;
    sys.A[sys.offA[bR] + pR * mL + pL] -= k;
    sys.B[sys.offB[bR] + pR * mR + pR] += k;
  }

  // Scharfetter–Gummel fluxes along a dilute region's segments (nodes g0…g1), written
  // by local index. The flux
  // N = g[B(Δ)c_L − B(−Δ)c_R] is rewritten with B(−Δ) = B(Δ)e^Δ and c_R e^Δ = c_L e^{Δη} as
  // N = −g·B(Δ)·c_L·expm1(Δη). This is precise relative to the quasi-Fermi difference Δη, so tiny
  // fluxes (e.g. majority carriers carrying a small current) don't vanish in the cancellation of
  // two huge drift and diffusion terms.
  _segmentsDilute(g0, g1, region, mat, vel) {
    const { n, M, u, uLo, res, c, z, sys, loc, lin } = this, R = this.rix, sg = this.sg;
    const { A: JA, B: JB, C: JC, sizes, offA, offB, offC } = sys;
    const { segLength, segArea } = this.model.grid;
    for (let s = g0; s < g1; s++) {
      const bL = s + region, bR = bL + 1, h = segLength[s] / segArea[s]; // (fluxes are totals through the cross-section)
      const mL = sizes[bL], mR = sizes[bR], pL = loc[bL * M], pR = loc[bR * M];
      const phiL = u[bL * M], phiR = u[bR * M];
      for (let i = 0; i < n; i++) {
        if (!mat.present[i] || mat.D[i] === 0) continue;
        const r = 1 + i, zi = z[i];
        const cL = c[s * n + i];
        const g = mat.D[i] / h;
        const pe = (vel * h) / mat.D[i]; // advection: a Péclet shift of the drift potential
        const d = zi * (phiR - phiL) - pe;
        const deta = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]) - pe;
        sgFlux(g, d, deta, cL, c[(s + 1) * n + i], sg);
        const N = sg[0], dNdd = sg[1];
        const dNdEtaL = sg[2];
        const dNdEtaR = -sg[3];
        const dNdPhiL = -zi * dNdd - zi * N; // via Δ, and via c_L ∝ e^{−zφ̂_L}
        const dNdPhiR = zi * dNdd;
        res[R[bL * M + r]] += N;
        res[R[bR * M + r]] -= N;
        if (lin) {
          lin.segment(R[bL * M + r], R[bR * M + r], R[bL * M], R[bR * M], zi, N, dNdEtaR, zi * dNdd);
          continue;
        }
        // Row r of block bL (couplings to itself in B, to bR in C) and of bR (to bL in A).
        const lL = loc[bL * M + r], lR = loc[bR * M + r];
        const oBL = offB[bL] + lL * mL, oCL = offC[bL] + lL * mR, oAR = offA[bR] + lR * mL, oBR = offB[bR] + lR * mR;
        JB[oBL + lL] += dNdEtaL;
        JC[oCL + lR] += dNdEtaR;
        JA[oAR + lL] -= dNdEtaL;
        JB[oBR + lR] -= dNdEtaR;
        if (zi !== 0 && pL >= 0) {
          JB[oBL + pL] += dNdPhiL;
          JC[oCL + pR] += dNdPhiR;
          JA[oAR + pL] -= dNdPhiL;
          JB[oBR + pR] -= dNdPhiR;
        }
      }
    }
  }

  // Scharfetter–Gummel with non-ideal statistics. The excess ex = ζ − ln(c/c_ref) acts as an
  // extra potential, linear along the segment like φ, so Δ = zΔφ̂ + Δex and
  //   N = −(D/h)·B(Δ)·c_L·expm1(η_R − η_L),
  // still exactly zero at equilibrium. c_L and ex depend on every ζ at their node through K.
  _segmentConcentrated(s, bL, bR, mat, h, vel) {
    const R = this.rix;
    const { n, M, res, z, K } = this;
    const KL = s * n * n, KR = (s + 1) * n * n;
    for (let i = 0; i < n; i++) {
      if (!mat.present[i] || mat.D[i] === 0) continue;
      const r = 1 + i, zi = z[i];
      const { N, gBcE1, dNdd, fL, cR } = this._fluxConcentrated(s, bL, bR, mat, h, vel, i);
      res[R[bL * M + r]] += N;
      res[R[bR * M + r]] -= N;
      if (this.lin) {
        // In difference form (see DifferenceTerms): b on Δη_i and z_i·∂N/∂d on Δφ̂, and the rest
        // as side terms in each δη_j − z_j δφ̂ at either end (c_L and the excess potentials).
        const rL = R[bL * M + r], rR = R[bR * M + r], pL = R[bL * M], pR = R[bR * M];
        this.lin.segment(rL, rR, pL, pR, zi, 0, -gBcE1, zi * dNdd);
        for (let j = 0; j < n; j++) {
          const wL = fL * K[KL + i * n + j] - (j === i ? dNdd : 0);
          const wR = cR > 0 ? dNdd * ((j === i ? 1 : 0) - K[KR + i * n + j] / cR) : 0;
          if (wL !== 0) this.lin.side(rL, rR, R[bL * M + 1 + j], pL, z[j], wL);
          if (wR !== 0) this.lin.side(rL, rR, R[bR * M + 1 + j], pR, z[j], wR);
        }
        continue;
      }
      const { jL, jR } = this;
      for (let k = 0; k < M; k++) {
        if (jL[k] !== 0) {
          this._j(bL, r, bL, k, jL[k]);
          this._j(bR, r, bL, k, -jL[k]);
        }
        if (jR[k] !== 0) {
          this._j(bL, r, bR, k, jR[k]);
          this._j(bR, r, bR, k, -jR[k]);
        }
      }
    }
  }

  // Species i's flux across segment s of a non-ideal region, toward +x, and its derivatives in
  // this.jL and this.jR (over the unknowns of the blocks bL and bR).
  _fluxConcentrated(s, bL, bR, mat, h, vel, i) {
    const { n, M, u, uLo, c, z, K, ex, jL, jR } = this;
    const gL = s, gR = s + 1, KL = gL * n * n, KR = gR * n * n, r = 1 + i, zi = z[i];
    const cL = c[gL * n + i], cR = c[gR * n + i];
    const g = mat.D[i] / h;
    const pe = (vel * h) / mat.D[i];
    const d = zi * (u[bR * M] - u[bL * M]) + ex[gR * n + i] - ex[gL * n + i] - pe;
    const deta = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]) - pe;
    const sg = sgFlux(g, d, deta, cL, cR, this.sg);
    const N = sg[0], dNdd = sg[1], gBcE1 = sg[3]; // (gBcE1 = g·B·c_L·e^{Δη})
    // Left node: ∂N/∂ζ_Lj = −gE(B + B′)K_L,ij − dNdd·δ_ij, plus the direct η_L,i term.
    const fL = sg[4];
    jL[0] = -fL * this._Kz(gL, i);
    for (let j = 0; j < n; j++) jL[1 + j] = fL * K[KL + i * n + j];
    jL[r] += -dNdd + gBcE1;
    // Right node: ∂N/∂ζ_Rj = dNdd·(δ_ij − K_R,ij/c_R), plus the direct η_R,i term.
    jR.fill(0);
    if (cR > 0) {
      const fR = -dNdd / cR;
      jR[0] = -fR * this._Kz(gR, i);
      for (let j = 0; j < n; j++) jR[1 + j] = fR * K[KR + i * n + j];
    } else {
      jR[0] = dNdd * zi; // K/c → δ as c → 0
      jR[r] -= dNdd;
    }
    jR[r] += dNdd - gBcE1;
    return { N, gBcE1, dNdd, fL, cR };
  }

  // Ohmic conduction in a metal, in mixed form. The flux J (slot 0 of the left node) obeys
  //   η_R − η_L + J/g = 0,  g = σ RT/(z²F² h)   (Ohm's law, J = −g Δη),
  // and enters the two nodes' carrier balances as outflow and inflow.
  _segmentConductor(s, bL, bR, mat, h) {
    const R = this.rix;
    const { M, u, uLo, res, z } = this;
    const i = mat.conductor.i, r = 1 + i;
    const g = (mat.conductor.sigma * this.model.RT) / (z[i] * z[i] * FARADAY * FARADAY * h);
    const J = u[bL * M];
    res[R[bL * M]] = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]) + J / g;
    this._j(bL, 0, bR, r, 1);
    this._j(bL, 0, bL, r, -1);
    this._j(bL, 0, bL, 0, 1 / g);
    res[R[bL * M + r]] += J;
    res[R[bR * M + r]] -= J;
    this._j(bL, r, bL, 0, 1);
    this._j(bR, r, bL, 0, -1);
  }

  // Eddy mixing: N = −(D_mix/RT) P ∇μ̄ with P = K − (Kz)(Kz)ᵀ/(zᵀKz), K = ∂c/∂ζ. It mixes
  // composition without carrying current (zᵀP = 0) and vanishes exactly at equilibrium; in a
  // neutral solution it's −D_mix ∇c for every species, under any statistics (eddies carry parcels
  // of fluid, so they mix concentration). On a segment, N_i = −(D_mix/h) Σ_j P̄_ij Δη_j, with
  // P̄ from K̄ = C̄, the logarithmic mean of each c, where both ends are ideal (so K̄ Δζ = Δc
  // exactly), and from the secant K̄ of _mixCrowded otherwise. Only mobile species (D > 0) take part.
  _segmentMixing(s, bL, bR, Dm, h, mat) {
    const R = this.rix;
    const { n, M, u, uLo, res, c, z } = this;
    const gL = s, gR = s + 1, k0 = Dm / h;
    const on = mat.mobile ?? (mat.mobile = [...Array(n).keys()].filter((i) => mat.present[i] && mat.D[i] > 0));
    if (on.length === 0) return;
    if (!this.nodeIdeal[gL] || !this.nodeIdeal[gR]) return this._mixCrowded(gL, gR, bL, bR, on, k0);
    const w = this.mixWork ?? (this.mixWork = [0, 1, 2, 3].map(() => new Float64Array(n)));
    const [cb, dLa, dLb, de] = w;
    let S = 0, Q = 0;
    for (const k of on) {
      const a = c[gL * n + k], b = c[gR * n + k];
      // l = ln(c_R/c_L), from the compensated potentials (statistics ideal at both ends)
      const l = u[bR * M + 1 + k] - u[bL * M + 1 + k] + (uLo[bR * M + 1 + k] - uLo[bL * M + 1 + k]) - z[k] * (u[bR * M] - u[bL * M]);
      logMean(a, b, l, cb, dLa, dLb, k);
      de[k] = u[bR * M + 1 + k] - u[bL * M + 1 + k] + (uLo[bR * M + 1 + k] - uLo[bL * M + 1 + k]);
      S += z[k] * z[k] * cb[k];
      Q += z[k] * cb[k] * de[k];
    }
    const q = S > 0 ? Q / S : 0;
    const { jL, jR } = this;
    const dcL = this.dA, dcR = this.dB; // ∂c_k/∂slot at each end
    for (const i of on) {
      const r = 1 + i;
      const N = -k0 * cb[i] * (de[i] - z[i] * q);
      jL.fill(0);
      jR.fill(0);
      for (const j of on) {
        // direct dependence on Δη_j
        const P = cb[i] * ((i === j ? 1 : 0) - (S > 0 ? (z[i] * z[j] * cb[j]) / S : 0));
        jL[1 + j] += k0 * P;
        jR[1 + j] -= k0 * P;
      }
      // through c̄_k: ∂a_i/∂c̄_k = δ_ik (Δη_i − z_i q) − c̄_i z_i (z_k Δη_k − q z_k²)/S
      for (const k of on) {
        let da = (i === k ? de[i] - z[i] * q : 0) - (S > 0 ? (cb[i] * z[i] * (z[k] * de[k] - q * z[k] * z[k])) / S : 0);
        if (da === 0) continue;
        da *= -k0;
        this._dc(gL, k, dcL);
        this._dc(gR, k, dcR);
        for (let t = 0; t < M; t++) {
          jL[t] += da * dLa[k] * dcL[t];
          jR[t] += da * dLb[k] * dcR[t];
        }
      }
      res[R[bL * M + r]] += N;
      res[R[bR * M + r]] -= N;
      for (let t = 0; t < M; t++) {
        if (jL[t] !== 0) {
          this._j(bL, r, bL, t, jL[t]);
          this._j(bR, r, bL, t, -jL[t]);
        }
        if (jR[t] !== 0) {
          this._j(bL, r, bR, t, jR[t]);
          this._j(bR, r, bR, t, -jR[t]);
        }
      }
    }
  }

  // Eddy mixing where either end is crowded (non-ideal): N = −k0 P̄ Δη as in the dilute case,
  // with C̄ replaced by K̄ = ½(R_L K_L R_L + R_R K_R R_R), R = diag(√(c̄/c)) at each end (c̄ the
  // logarithmic mean): each end's K rescaled to the segment's mean concentrations, and averaged.
  // It's C̄ where statistics are ideal, symmetric and positive definite, and smooth in the state,
  // so P̄ keeps everything the dilute term has (no current, exactly zero at a level μ̄); and in
  // a neutral solution it gives −k0 Δc to second order in the segment's Δζ. Its Jacobian runs
  // through K's own slope, dK (_statSlopes).
  _mixCrowded(gL, gR, bL, bR, on, k0) {
    const R = this.rix;
    const { n, M, res } = this;
    const m = on.length, W = this._mixBase(gL, gR, bL, bR, on);
    const { N, a, b, cb, dLa, dLb, rL, rR, KL, KR, Kb, de, zc, kz, tK, tg, tkz, tr, dcs } = W;
    const { S, Q } = W;
    const Jl = W.Jl, Jr = W.Jr; // ∂N_p/∂slot t at either end, m × (n + 1)
    for (const end of [0, 1]) {
      const gX = end ? gR : gL, cX = end ? b : a, rX = end ? rR : rL, KX = end ? KR : KL, J = end ? Jr : Jl;
      const dLX = end ? dLb : dLa, sgn = end ? 1 : -1;
      for (let q = 0; q < m; q++) this._dc(gX, on[q], dcs.subarray(q * M, (q + 1) * M));
      for (let t = 0; t <= n; t++) {
        // the tangent of K̄ along slot t: through c at this end (its own ρ and c̄) and K there
        for (let p = 0; p < m; p++) {
          const tc = dcs[p * M + t], tcb = dLX[p] * tc;
          tr[p] = 0.5 * rX[p] * (tcb / cb[p] - tc / cX[p]); // tρ at this end
          W.tcbr[p] = tcb; // and c̄'s, which moves the other end's ρ too
        }
        for (let p = 0; p < m; p++) {
          for (let q = 0; q < m; q++) {
            const dk = this._dKslot(gX, on[p], on[q], t);
            const rO = end ? rL : rR, KO = end ? KL : KR;
            const trO_p = 0.5 * rO[p] * (W.tcbr[p] / cb[p]), trO_q = 0.5 * rO[q] * (W.tcbr[q] / cb[q]);
            tK[p * m + q] =
              0.5 * ((tr[p] * rX[q] + rX[p] * tr[q]) * KX[p * m + q] + rX[p] * rX[q] * dk) +
              0.5 * (trO_p * rO[q] + rO[p] * trO_q) * KO[p * m + q];
          }
        }
        for (let p = 0; p < m; p++) {
          let x = 0, y = 0;
          for (let q = 0; q < m; q++) {
            x += tK[p * m + q] * de[q] + (t === 1 + on[q] ? sgn * Kb[p * m + q] : 0);
            y += tK[p * m + q] * zc[q];
          }
          tg[p] = x;
          tkz[p] = y;
        }
        let tS = 0, tQ = 0;
        for (let p = 0; p < m; p++) {
          tS += zc[p] * tkz[p];
          tQ += zc[p] * tg[p];
        }
        for (let p = 0; p < m; p++) {
          J[p * (n + 1) + t] = -(S > 0 ? tg[p] - (tkz[p] * Q + kz[p] * tQ) / S + (kz[p] * Q * tS) / (S * S) : tg[p]);
        }
      }
    }
    for (let p = 0; p < m; p++) {
      const r = 1 + on[p], Ni = k0 * N[p];
      res[R[bL * M + r]] += Ni;
      res[R[bR * M + r]] -= Ni;
      for (let t = 0; t <= n; t++) {
        const vl = k0 * Jl[p * (n + 1) + t], vr = k0 * Jr[p * (n + 1) + t];
        if (vl !== 0) {
          this._j(bL, r, bL, t, vl);
          this._j(bR, r, bL, t, -vl);
        }
        if (vr !== 0) {
          this._j(bL, r, bR, t, vr);
          this._j(bR, r, bR, t, -vr);
        }
      }
    }
  }

  // ∂K_ij/∂(slot t) at node g: through ζ, slot 1 + r moving ζ_r and slot 0 (φ̂) every ζ_r by −z_r.
  _dKslot(g, i, j, t) {
    const { n, z } = this;
    if (this.nodeIdeal[g]) {
      if (i !== j) return 0;
      const ci = this.c[g * n + i];
      return t === 0 ? -z[i] * ci : t === 1 + i ? ci : 0;
    }
    const o = g * n * n * n + (i * n + j) * n;
    if (t > 0) return this.dK[o + t - 1];
    let v = 0;
    for (let r = 0; r < n; r++) v -= this.dK[o + r] * z[r];
    return v;
  }

  // _mixCrowded's flux per k0 on a segment, into W.N (indexed as on), and what its Jacobian
  // needs. Also the solution's readout (segmentFluxes).
  _mixBase(gL, gR, bL, bR, on) {
    const { n, M, u, uLo, c, z } = this;
    const m = on.length;
    const W = this.mixCrowd ?? (this.mixCrowd = mixWork(n, M));
    const { a, b, cb, dLa, dLb, rL, rR, KL, KR, Kb, de, zc, g, kz, N } = W;
    for (let p = 0; p < m; p++) {
      const k = on[p], oL = bL * M + 1 + k, oR = bR * M + 1 + k;
      a[p] = c[gL * n + k];
      b[p] = c[gR * n + k];
      zc[p] = z[k];
      logMean(a[p], b[p], Math.log(b[p] / a[p]), cb, dLa, dLb, p);
      rL[p] = Math.sqrt(cb[p] / a[p]);
      rR[p] = Math.sqrt(cb[p] / b[p]);
      de[p] = u[oR] - u[oL] + (uLo[oR] - uLo[oL]);
    }
    for (let p = 0; p < m; p++) {
      for (let q = 0; q < m; q++) {
        KL[p * m + q] = this._Kat(gL, on[p], on[q]);
        KR[p * m + q] = this._Kat(gR, on[p], on[q]);
        Kb[p * m + q] = 0.5 * (rL[p] * rL[q] * KL[p * m + q] + rR[p] * rR[q] * KR[p * m + q]);
      }
    }
    let S = 0, Q = 0;
    for (let p = 0; p < m; p++) {
      let x = 0, y = 0;
      for (let q = 0; q < m; q++) {
        x += Kb[p * m + q] * de[q];
        y += Kb[p * m + q] * zc[q];
      }
      g[p] = x;
      kz[p] = y;
    }
    for (let p = 0; p < m; p++) {
      S += zc[p] * kz[p];
      Q += zc[p] * g[p];
    }
    W.S = S;
    W.Q = Q;
    for (let p = 0; p < m; p++) N[p] = -(S > 0 ? g[p] - (kz[p] * Q) / S : g[p]);
    return W;
  }

  // K_ij at node g.
  _Kat(g, i, j) {
    const { n } = this;
    if (this.nodeIdeal[g]) return i === j ? this.c[g * n + i] : 0;
    return this.K[g * n * n + i * n + j];
  }

  // ∂c_k/∂slot at node g into d (slot 0 = φ̂, 1 + j = η_j).
  _dc(g, k, d) {
    const { n } = this;
    d.fill(0);
    if (this.nodeIdeal[g]) {
      const ck = this.c[g * n + k];
      d[1 + k] = ck;
      d[0] = -this.z[k] * ck;
      return;
    }
    const o = g * n * n + k * n;
    for (let j = 0; j < n; j++) d[1 + j] = this.K[o + j];
    d[0] = -this._Kz(g, k);
  }

  // A bulk reaction's rate at node g (mol/(m³·s), forward), as assembled: for the solution.
  bulkRate(rx, g) {
    const { n, M, u, uLo, c } = this, b = this.blockOfNode[g];
    let aHi = rx.fixedA, aLo = 0;
    for (let p = 0; p < rx.sp.length; p++) {
      aHi += rx.nu[p] * u[b * M + 1 + rx.sp[p]];
      aLo += rx.nu[p] * uLo[b * M + 1 + rx.sp[p]];
    }
    const law = rx.srh?.[this.nodeMaterial[g]];
    if (law) return this._srhRate(law, c[g * n + law.n.i], c[g * n + law.p.i], aHi + aLo, ...this._srhDegen(g, g, law)).rate;
    let P = rx.generation ? rx.kfNode[g] * this.generationScale : rx.kfNode[g];
    for (let p = 0; p < rx.sp.length; p++) if (rx.nu[p] > 0) P *= rx.nu[p] === 1 ? c[g * n + rx.sp[p]] : powi(c[g * n + rx.sp[p]], rx.nu[p]);
    return P * -Math.expm1(-(aHi + aLo));
  }

  // A bulk reaction's forward one-way rate at node g (mol/(m³·s)); its backward one is that
  // times e^{−a}, so their difference is bulkRate's.
  bulkForward(rx, g) {
    const { n, M, u, uLo, c } = this, b = this.blockOfNode[g];
    const law = rx.srh?.[this.nodeMaterial[g]];
    if (law) {
      let a = rx.fixedA;
      for (let p = 0; p < rx.sp.length; p++) a += rx.nu[p] * (u[b * M + 1 + rx.sp[p]] + uLo[b * M + 1 + rx.sp[p]]);
      return this._srhRate(law, c[g * n + law.n.i], c[g * n + law.p.i], a, ...this._srhDegen(g, g, law)).forward;
    }
    let P = rx.generation ? rx.kfNode[g] * this.generationScale : rx.kfNode[g];
    for (let p = 0; p < rx.sp.length; p++) if (rx.nu[p] > 0) P *= rx.nu[p] === 1 ? c[g * n + rx.sp[p]] : powi(c[g * n + rx.sp[p]], rx.nu[p]);
    return P;
  }

  // The larger of a bulk reaction's two one-way rates at node g (mol/(m³·s)): the scale its net
  // rate is read against (in equilibrium the net rate is their round-off).
  bulkOneWay(rx, g) {
    const { n, M, u, uLo, c } = this, b = this.blockOfNode[g];
    let a = rx.fixedA;
    for (let p = 0; p < rx.sp.length; p++) a += rx.nu[p] * (u[b * M + 1 + rx.sp[p]] + uLo[b * M + 1 + rx.sp[p]]);
    const law = rx.srh?.[this.nodeMaterial[g]];
    if (law) return this._srhRate(law, c[g * n + law.n.i], c[g * n + law.p.i], a, ...this._srhDegen(g, g, law)).oneWay;
    let P = rx.generation ? rx.kfNode[g] * this.generationScale : rx.kfNode[g];
    for (let p = 0; p < rx.sp.length; p++) if (rx.nu[p] > 0) P *= rx.nu[p] === 1 ? c[g * n + rx.sp[p]] : powi(c[g * n + rx.sp[p]], rx.nu[p]);
    return P * Math.max(1, Math.exp(-a));
  }

  // SRH kinetics (see srhLaw in device.js): the rate, and its sensitivities to ln n, ln p and
  // SRH kinetics (see srhLaw in device.js) at electron and hole concentrations cn, cp and affinity
  // a = A/RT. With no n₁ given, the trap is midgap: n₁ = p₁ = n_i = √(n p e^{−a}). n₁p₁ = n p e^{−a}
  // (n_i² for Boltzmann carriers), and under degeneracy each carrier's emission scales with its own
  // γ = e^{−ex} (xn, xp: the carriers' excess chemical potentials): n₁ = γ_n n₁₀ and p₁ = γ_p p₁₀,
  // n₁₀ the law's n1, or n_i at midgap. (γ_n γ_p split evenly onto both, as it was, put a degenerate
  // n⁺ region's γ_n on its holes' emission too.) The rate's partial derivatives in ln n, ln p, a,
  // xn and xp, each with the others held (_srhInto takes them to the unknowns). For large negative
  // a (strong extraction), np·(1 − e^{−a}) as −n p e^{−a}.
  _srhRate(law, cn, cp, a, xn = 0, xp = 0) {
    const { tn, tp } = law, lnnp = Math.log(cn) + Math.log(cp);
    let lnn1, kN = 0, kP = 0, kA = 0, kX = 0; // ln n₁ = kN ln n + kP ln p + kA a + kX (xn − xp) + const
    if (Number.isNaN(law.n1)) {
      lnn1 = (lnnp - a - xn + xp) / 2; // γ_n n_i, with n_i² = n p e^{−a}/(γ_n γ_p)
      kN = 0.5;
      kP = 0.5;
      kA = -0.5;
      kX = -0.5;
    } else {
      lnn1 = Math.log(law.n1) - xn;
      kX = -1; // (and nothing in xp)
    }
    const n1 = Math.exp(lnn1), p1 = Math.exp(lnnp - a - lnn1), np = cn * cp;
    const den = tp * (cn + n1) + tn * (cp + p1);
    // np(1 − e^{−a}), with e^{−a} folded into np where it would overflow
    const f = -Math.expm1(-a), num = a > -700 ? np * f : np - Math.exp(lnnp - a);
    const rate = num / den;
    const dDn = tp * cn + tp * n1 * kN + tn * p1 * (1 - kN);
    const dDp = tp * n1 * kP + tn * cp + tn * p1 * (1 - kP);
    const dDa = tp * n1 * kA + tn * p1 * (-1 - kA);
    const dDxn = (tp * n1 - tn * p1) * kX, dDxp = Number.isNaN(law.n1) ? -dDxn : 0;
    return {
      rate,
      forward: np / den,
      oneWay: a > -700 ? (np / den) * Math.max(1, Math.exp(-a)) : Math.exp(lnnp - a) / den,
      Cn: rate * (1 - dDn / den),
      Cp: rate * (1 - dDp / den),
      Ca: (a > -700 ? (np / den) * (1 - f) : Math.exp(lnnp - a) / den) - (rate * dDa) / den,
      Cxn: (-rate * dDxn) / den,
      Cxp: (-rate * dDxp) / den,
    };
  }

  // _srhRate's excess chemical potentials: the electrons' at node gn, the holes' at gq.
  _srhDegen(gn, gq, law) {
    return [this.nodeIdeal[gn] ? 0 : this.ex[gn * this.n + law.n.i], this.nodeIdeal[gq] ? 0 : this.ex[gq * this.n + law.p.i]];
  }

  // Add C·∂(ln c_i) + Cx·∂(ex_i) at node g into d, over its slots: ex_i = ζ_i − ln c_i (0 where
  // the statistics are ideal), ζ_i = η_i − μ°_i/RT − z_i φ̂.
  _srhInto(g, i, C, Cx, d) {
    if (this.nodeIdeal[g] || Cx === 0) return this._dlnc(g, i, C, d);
    this._dlnc(g, i, C - Cx, d);
    d[1 + i] += Cx;
    d[0] -= this.z[i] * Cx;
  }

}
