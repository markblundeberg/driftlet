// Small signal: the impedance at one terminal about a steady state, by linearising the steady
// equations and adding the storage term's iω, and the complex profiles of the response.
//
// CONTRACT
//
// Methods of the Solver (src/solver.js; its layout is in docs/internals.md), installed on its
// prototype like the assembly's. Device.impedance (src/index.js) calls impedance() right after a
// steady solve that converged (src/steady.js), at the state it left. What they use of the rest:
//   this.assemble(dt) (src/assemble.js)  at dt = ∞ the steady Jacobian in this.sys (blocks A, B,
//                        C), and for every terminal while allTerminals is set, termB[k] = ∂res/∂V_k,
//                        termC[k] = ∂I_k/∂x, termDI[k] = ∂I_k/∂V_k, termDX[k][k2] = ∂I_k/∂V_k2,
//                        each compact (sys order). At a finite dt the same plus the storage
//                        term's derivatives over dt.
//   this._assembleDifference(dt) (src/difference.js)  the same Jacobian with the dilute kernels'
//                        terms kept apart in difference form (`lin`, applied exactly) and the rest
//                        in this.sysRest: J·v without the cancellation of an assembled matrix.
//   this._fluxConcentrated (src/regions.js), this.jL / this.jR, sgFlux (src/bernoulli.js)  a
//                        segment's flux and its derivatives, for the current cuts.
//   ComplexBlockTridiagonal (src/blockTridiagonal.js)  a complex block-tridiagonal factorisation.
//
// impedance(frequencies, { terminal = 'right', profiles = false })
//   → { f, Z: { re, im }, profiles? }, Z in Ω·m² (Ω for a non-planar device) at each frequency.
//   Throws a SolverError where the terminal doesn't exist, has a series resistance, or every
//   other terminal is driven by a current (no return path: Z infinite), or where a frequency
//   isn't positive.
//   Linearised about the present state: (J + iωS)·δx + Σ_k B_k δV_k = 0 for the grid's rows and,
//   for each floating terminal k (driven by a current), (C_k + iωC′_k)·δx + (DI_k + iωDI′_k) δV_k
//   + Σ_k2 (DX + iωDX′)_k,k2 δV_k2 = δI_k, with J the steady Jacobian (dt = ∞) and S, C′, DI′,
//   DX′ the storage parts, (X(dt) − X(∞))·dt at dt = 1e-30. The measured terminal: held at a
//   voltage, δV = 1 and Z = 1/δI, δI its current into the device; driven by a current, δI = 1
//   and Z = δV. Every other terminal keeps its drive: held, δV = 0 (AC ground); driven, δI = 0
//   (AC open).
//   Each frequency: J + iωS with each row divided by its largest entry, factorised, preconditions
//   GMRES (m 12, 8 restarts, tol 1e-8; its best iterate kept) on the exact J·v (difference form).
//   Solved: y, the response with the floating terminals held (to the measured one's δV, or none),
//   and X_k, to each floating terminal's unit δV_k; then the floating terminals' rows, a small
//   complex system, for their δV; δx = y − Σ X_k δV_k.
//   The current read (held terminal): DI + iωDI′ + Σ (DX + iωDX′) δV + (C + iωC′)·δx, its error
//   estimated from its terms' round-off and the solves' residuals; for a contact, also across
//   each cut (_currentCuts: a segment's linearised total current, conduction and displacement,
//   between that contact and the nearest port's window where there are ports), the reading with
//   the least error taken. With no floating terminal, cuts, and that error over 1% of the current:
//   the response split into each region's uniform shift (_regionShift) and a remainder solved to
//   GMRES's floor, the shift's effect taken only at the rows it can touch (_regionEdgeRows), and
//   that reading kept where its error is less.
//   profiles: per frequency, _smallSignalProfiles(δx): δφ (V), δμ̄_i (J/mol) and δc_i (mol/m³)
//   at each node, complex, per volt (held) or per unit current (driven); NaN where φ is undefined,
//   where a species is absent, and for a metal's concentration.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. u, uLo, termV and the time are as before the call; allTerminals as it was.
//   2. Z is the linear response of the equations the steady solve solved: as f → 0 it tends to
//      the steady state's differential dV/dI at the measured terminal (the other terminals'
//      drives held), and the readout doesn't depend on which cut or the contact it was read at
//      beyond its error estimate.
//   3. A term the steady system holds by a law in place of a balance row (src/lawRows.js) is in
//      the linearisation only as the assembly's rows are: steady rows replaced by laws are not
//      in J (the assembly here runs with this.steady false).
//   4. Deterministic: the same state, frequencies and options give the same Z, bit for bit.
// Also written, left as the last assembly at dt = ∞ left them: this.sourceTime (the present
// time), cOld, thOld, th0Old, gateOld (the present state's), contactDStart, portQStart, the
// residual, the Jacobian and the terminal terms, sysRest and the concentrations.
// Tests: test/impedance.test.js (Macdonald's blocking electrodes, a junction's depletion
// capacitance, a MOS capacitor at both terminals, ports, profiles, Fermi–Dirac against
// Boltzmann); the demos' EIS.

import { ComplexBlockTridiagonal } from './blockTridiagonal.js';
import { sgFlux } from './bernoulli.js';
import { EPS0, FARADAY } from './constants.js';
import { SolverError } from './solver.js';

// Solve the K×K complex system (Ar + i Ai) x = (last column), partial pivoting by modulus.
function complexSolve(Ar, Ai, K) {
  for (let col = 0; col < K; col++) {
    let p = col;
    for (let r = col + 1; r < K; r++) if (Math.hypot(Ar[r][col], Ai[r][col]) > Math.hypot(Ar[p][col], Ai[p][col])) p = r;
    [Ar[col], Ar[p]] = [Ar[p], Ar[col]];
    [Ai[col], Ai[p]] = [Ai[p], Ai[col]];
    const dr = Ar[col][col], di = Ai[col][col], d2 = dr * dr + di * di;
    for (let r = col + 1; r < K; r++) {
      const fr = (Ar[r][col] * dr + Ai[r][col] * di) / d2, fi = (Ai[r][col] * dr - Ar[r][col] * di) / d2;
      for (let k = col; k <= K; k++) {
        const ur = Ar[col][k], ui = Ai[col][k];
        Ar[r][k] -= fr * ur - fi * ui;
        Ai[r][k] -= fr * ui + fi * ur;
      }
    }
  }
  const re = new Float64Array(K), im = new Float64Array(K);
  for (let r = K - 1; r >= 0; r--) {
    let vr = Ar[r][K], vi = Ai[r][K];
    for (let k = r + 1; k < K; k++) {
      vr -= Ar[r][k] * re[k] - Ai[r][k] * im[k];
      vi -= Ar[r][k] * im[k] + Ai[r][k] * re[k];
    }
    const dr = Ar[r][r], di = Ai[r][r], d2 = dr * dr + di * di;
    re[r] = (vr * dr + vi * di) / d2;
    im[r] = (vi * dr - vr * di) / d2;
  }
  return { re, im };
}

// GMRES(m), restarted, for a complex system op(x) = b, with vectors as [re, im] pairs and op
// already preconditioned; x holds the starting guess and receives the answer. A cycle's own
// residual estimate can drift far from the true one (a preconditioner nearly singular at low
// frequency, or an op whose round-off the preconditioner amplifies), so each cycle is judged by
// the true residual, at the next one's start: it returns once that's below tol, and gives up when
// two cycles in a row fail to lower it, leaving x at the best iterate. Returns the true residuals (relative
// to b) at the start and at that iterate.
function gmres(op, b, x, { m = 12, restarts = 8, tol = 1e-8 } = {}) {
  const N = b[0].length;
  const cdot = (a, c) => {
    let re = 0, im = 0;
    for (let j = 0; j < N; j++) {
      re += a[0][j] * c[0][j] + a[1][j] * c[1][j];
      im += a[0][j] * c[1][j] - a[1][j] * c[0][j];
    }
    return [re, im]; // conj(a)·c
  };
  const norm = (a) => Math.sqrt(cdot(a, a)[0]);
  const bn = norm(b) || 1;
  const best = [Float64Array.from(x[0]), Float64Array.from(x[1])];
  let initial = Infinity, least = Infinity, idle = 0;
  for (let cycle = 0; ; cycle++) {
    const ax = op(x), r = [Float64Array.from(b[0], (v, j) => v - ax[0][j]), Float64Array.from(b[1], (v, j) => v - ax[1][j])];
    const beta = norm(r) / bn;
    if (cycle === 0) initial = beta;
    idle = beta < least ? 0 : idle + 1;
    if (idle === 0) {
      least = beta;
      best[0].set(x[0]);
      best[1].set(x[1]);
    }
    if (!(beta > tol) || cycle === restarts || idle === 2) break;
    const V = [[r[0].map((v) => v / (beta * bn)), r[1].map((v) => v / (beta * bn))]];
    const H = Array.from({ length: m + 1 }, () => Array.from({ length: m }, () => [0, 0]));
    const cs = [], sn = [], g = Array.from({ length: m + 1 }, () => [0, 0]);
    g[0] = [beta * bn, 0];
    let k = 0;
    for (; k < m; k++) {
      const w = op(V[k]);
      for (let i = 0; i <= k; i++) {
        const [hr, hi] = cdot(V[i], w);
        H[i][k] = [hr, hi];
        for (let j = 0; j < N; j++) {
          w[0][j] -= hr * V[i][0][j] - hi * V[i][1][j];
          w[1][j] -= hr * V[i][1][j] + hi * V[i][0][j];
        }
      }
      const hn = norm(w);
      H[k + 1][k] = [hn, 0];
      V.push([w[0].map((v) => v / (hn || 1)), w[1].map((v) => v / (hn || 1))]);
      // The earlier rotations, then a new one zeroing H[k+1][k]: x' = c x + s y, y' = −s̄ x + c y.
      for (let i = 0; i < k; i++) {
        const [xr, xi] = H[i][k], [yr, yi] = H[i + 1][k], c = cs[i], [sr, si] = sn[i];
        H[i][k] = [c * xr + sr * yr - si * yi, c * xi + sr * yi + si * yr];
        H[i + 1][k] = [-(sr * xr + si * xi) + c * yr, -(sr * xi - si * xr) + c * yi];
      }
      const [ar, ai] = H[k][k], am = Math.hypot(ar, ai), t = Math.hypot(am, hn);
      const c = t > 0 ? am / t : 1, s = am > 0 ? [(ar / am) * (hn / t), (ai / am) * (hn / t)] : [1, 0];
      cs.push(c);
      sn.push(s);
      H[k][k] = am > 0 ? [(ar / am) * t, (ai / am) * t] : [hn, 0];
      H[k + 1][k] = [0, 0];
      const [gr, gi] = g[k];
      g[k] = [c * gr, c * gi];
      g[k + 1] = [-(s[0] * gr + s[1] * gi), -(s[0] * gi - s[1] * gr)];
      if (!(Math.hypot(g[k + 1][0], g[k + 1][1]) > tol * bn)) {
        k++;
        break;
      }
    }
    // Back-substitute H y = g over the first k columns, then x += V y.
    const y = Array.from({ length: k }, () => [0, 0]);
    for (let i = k - 1; i >= 0; i--) {
      let [vr, vi] = g[i];
      for (let j = i + 1; j < k; j++) {
        const [hr, hi] = H[i][j], [yr, yi] = y[j];
        vr -= hr * yr - hi * yi;
        vi -= hr * yi + hi * yr;
      }
      const [dr, di] = H[i][i], d2 = dr * dr + di * di;
      y[i] = [(vr * dr + vi * di) / d2, (vi * dr - vr * di) / d2];
    }
    for (let i = 0; i < k; i++) {
      const [yr, yi] = y[i];
      for (let j = 0; j < N; j++) {
        x[0][j] += yr * V[i][0][j] - yi * V[i][1][j];
        x[1][j] += yr * V[i][1][j] + yi * V[i][0][j];
      }
    }
  }
  x[0].set(best[0]);
  x[1].set(best[1]);
  return { initial, residual: least };
}

export class SmallSignal {
  /**
   * Small-signal impedance about the current (steady) state, seen at one terminal. Linearising
   * gives (J + iωM)·δx + Σ B_k δV_k = rhs for the grid's rows, and (C_k + iωC′_k)·δx +
   * (∂I_k/∂V_k)·δV_k = δI_set,k for each floating terminal's circuit row, with J the steady
   * Jacobian and M the storage part (from J(dt) = J + M/dt), B_k = ∂res/∂V_k and C_k = ∂I_k/∂x.
   * At the measured terminal, a held voltage is perturbed (rhs = −B δV) and its current read
   * from C, or a driven current is perturbed and its voltage solved for. The other terminals
   * keep their drives: held ones at AC ground, driven ones at AC open circuit.
   * Z = δV/δI with I into the device (positive real part for a passive device).
   * @param {ArrayLike<number>} frequencies Hz
   */
  impedance(frequencies, opts) {
    const was = this.allTerminals;
    this.allTerminals = true; // (it reads its terminal's ∂res/∂V and ∂I/∂x)
    try {
      return this._impedance(frequencies, opts);
    } finally {
      this.allTerminals = was;
    }
  }

  _impedance(frequencies, { profiles = false, terminal = 'right' } = {}) {
    const { M, nB, sys, VT } = this;
    const kT = this.terms.findIndex((t) => t.name === terminal);
    if (kT < 0) throw new SolverError(`impedance: no terminal named '${terminal}' (${this.terms.map((t) => t.name).join(', ')})`);
    const dT = this.terms[kT].drive;
    if (dT.R > 0) throw new SolverError(`impedance: terminal '${terminal}' has a series resistance; that belongs to the external circuit`);
    // A small signal returns through the terminals held at a voltage: with every other one driven
    // by a current (held open), none comes back, and the impedance is infinite (else round-off).
    const others = this.terms.filter((t, k) => k !== kT);
    if (others.length && others.every((t) => t.drive.kind === 'I')) {
      throw new SolverError(
        `impedance: every other terminal (${others.map((t) => t.name).join(', ')}) is driven by a current, so it is open to a small signal and no current returns through '${terminal}': ` +
          `its impedance is infinite. Measure at one of those instead (impedance(f, { terminal: '${others[0].name}' }), its voltage against '${terminal}'), or hold one of them at a voltage`,
      );
    }
    const N = nB * M, fl = this.floating, K = fl.length;
    this.sourceTime = this.time;
    this.computeConcentrations();
    this.cOld.set(this.c);
    this.thOld.set(this.th);
    this.th0Old.set(this.th0);
    this.gateList.forEach((gt, j) => (this.gateOld[j] = this.u[gt.o]));
    this.assemble(Infinity);
    this.contactDStart = { ...this.contactD };
    this.portQStart = Float64Array.from(this.portQ);
    this.assemble(Infinity);
    // The steady parts, then the storage parts from an assembly at a tiny dt.
    const J = { A: sys.A.slice(), B: sys.B.slice(), C: sys.C.slice() };
    const Bt = this.termB.map((v) => v.slice()), Ct = this.termC.map((v) => v.slice()), DI = Float64Array.from(this.termDI), DX = this.termDX.map((v) => v.slice());
    const dts = 1e-30;
    this.assemble(dts);
    const S = {};
    for (const X of ['A', 'B', 'C']) S[X] = sys[X].map((v, k) => (v - J[X][k]) * dts);
    const Cs = this.termC.map((v, k) => v.map((x, j) => (x - Ct[k][j]) * dts));
    const DIs = this.termDI.map((x, k) => (x - DI[k]) * dts);
    const DXs = this.termDX.map((v, k) => v.map((x, j) => (x - DX[k][j]) * dts));
    // Where the measured contact's current may be read across a segment: anywhere without
    // ports; with them, only between it and the nearest window (a port's exchange takes current
    // in or out along it). Never for a port's own terminal.
    let cuts = kT < 2 ? this._currentCuts() : null;
    if (cuts && this.model.ports.length > 0) {
      const windows = this.model.ports.flatMap((p) => Array.from(p.nodes));
      const lo = Math.min(...windows), hi = Math.max(...windows);
      cuts = cuts.filter((cut) => (kT === 0 ? cut.s + 1 <= lo : cut.s >= hi));
    }
    // J·v for GMRES, exact: the assembled J holds a flux's dependence on η_L and η_R as two
    // entries, and where both are huge (an inversion layer) and v nearly uniform, J·v loses the
    // flux to round-off. So the dilute kernels' terms are kept in difference form (see
    // DifferenceTerms), with the rest from a matrix of their own, steady (Jr) and storage (Sr).
    this._assembleDifference(Infinity);
    const Jr = { A: this.sysRest.A.slice(), B: this.sysRest.B.slice(), C: this.sysRest.C.slice() };
    const { lin } = this._assembleDifference(dts);
    const Sr = {};
    for (const X of ['A', 'B', 'C']) Sr[X] = this.sysRest[X].map((v, k) => (v - Jr[X][k]) * dts);
    this.assemble(Infinity);

    const { sizes, offA, offB, offC, offX } = sys;
    const csys = new ComplexBlockTridiagonal(nB, sizes);
    const Nc = csys.size;
    const rr = new Float64Array(Nc), ri = new Float64Array(Nc);
    const yr = new Float64Array(Nc), yi = new Float64Array(Nc);
    const Xr = fl.map(() => new Float64Array(Nc)), Xi = fl.map(() => new Float64Array(Nc));
    const xr = new Float64Array(N), xi = new Float64Array(N);
    const rRes = new Float64Array(Nc), iRes = new Float64Array(Nc);
    const jr = new Float64Array(Nc), ji = new Float64Array(Nc), sr = new Float64Array(Nc), si = new Float64Array(Nc);
    const out = { f: Float64Array.from(frequencies), Z: { re: new Float64Array(frequencies.length), im: new Float64Array(frequencies.length) } };
    if (profiles) out.profiles = [];
    // (C + iωC′)·(ar + i ai) for terminal k
    const dot = (k, w, ar, ai) => {
      let re = 0, im = 0;
      const c = Ct[k], cs = Cs[k];
      for (let j = 0; j < Nc; j++) {
        re += c[j] * ar[j] - w * cs[j] * ai[j];
        im += c[j] * ai[j] + w * cs[j] * ar[j];
      }
      return [re, im];
    };
    Array.from(frequencies).forEach((f, q) => {
      const w = 2 * Math.PI * f;
      if (!(w > 0)) throw new SolverError('impedance: frequencies must be positive');
      // Rows scaled by their largest entry (and the terminal columns with them).
      const scale = new Float64Array(Nc);
      for (let blk = 0; blk < nB; blk++) {
        const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
        for (let r = 0; r < m; r++) {
          const rows = [['A', offA[blk] + r * mp, mp], ['B', offB[blk] + r * m, m], ['C', offC[blk] + r * mn, mn]];
          let mx = 0;
          for (const [X, o, len] of rows) for (let c = 0; c < len; c++) mx = Math.max(mx, Math.abs(J[X][o + c]), w * Math.abs(S[X][o + c]));
          const sc = mx > 0 ? 1 / mx : 1;
          scale[offX[blk] + r] = sc;
          for (const [X, o, len] of rows) {
            const JX = J[X], SX = S[X], Xre = csys[X + 'r'], Xim = csys[X + 'i'];
            for (let c = 0; c < len; c++) {
              Xre[o + c] = JX[o + c] * sc;
              Xim[o + c] = w * SX[o + c] * sc;
            }
          }
        }
      }
      csys.factor();
      const sx = (X, vr, out) => {
        for (let blk = 0; blk < nB; blk++) {
          const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
          for (let r = 0; r < m; r++) {
            let acc = 0;
            for (let c = 0; c < mp; c++) acc += X.A[offA[blk] + r * mp + c] * vr[offX[blk - 1] + c];
            for (let c = 0; c < m; c++) acc += X.B[offB[blk] + r * m + c] * vr[offX[blk] + c];
            for (let c = 0; c < mn; c++) acc += X.C[offC[blk] + r * mn + c] * vr[offX[blk + 1] + c];
            out[offX[blk] + r] = acc;
          }
        }
      };
      // A x = −σ ∂res/∂V_k, x = outr + i outi (σ = 1 for the response to the measured terminal,
      // −1 for a floating one's unit δV): GMRES on the system preconditioned by the factorised one
      // (whose error is concentrated in a few slow modes, which it finds in a few iterations).
      const zero = new Float64Array(Nc);
      const precondition = (ar, ai) => {
        const outR = new Float64Array(Nc), outI = new Float64Array(Nc);
        for (let j = 0; j < Nc; j++) {
          rRes[j] = ar[j] * scale[j];
          iRes[j] = ai[j] * scale[j];
        }
        csys.solve(rRes, iRes, outR, outI);
        return [outR, outI];
      };
      // (J + iωS)·v, exactly (see DifferenceTerms), into (ar, ai).
      const jv = (vr, vi, ar, ai) => {
        sx(Jr, vr, jr);
        sx(Jr, vi, ji);
        sx(Sr, vi, si);
        sx(Sr, vr, sr);
        for (let j = 0; j < Nc; j++) {
          ar[j] = jr[j] - w * si[j];
          ai[j] = ji[j] + w * sr[j];
        }
        lin.apply(vr, ar, 1, 0);
        lin.apply(vi, ar, 0, -w);
        lin.apply(vi, ai, 1, 0);
        lin.apply(vr, ai, 0, w);
      };
      const op = ([vr, vi]) => {
        const ar = new Float64Array(Nc), ai = new Float64Array(Nc);
        jv(vr, vi, ar, ai);
        return precondition(ar, ai);
      };
      const solve = (br, bi, outr, outi, k, sigma) => {
        const rhs = k >= 0 ? precondition(Bt[k].map((v) => -sigma * v), zero) : [new Float64Array(Nc), new Float64Array(Nc)];
        outr.set(rhs[0]);
        outi.set(rhs[1]);
        // (GMRES keeps its best iterate, which is never worse than the factorised solve.)
        const { residual } = gmres(op, rhs, [outr, outi]);
        accuracy = Math.max(accuracy, residual);
      };
      let accuracy = 1e-14; // (relative, of the solves: the readouts' error is weighed with it)
      // y: the response with the floating terminals held; X_k: to a unit δV_k.
      const heldT = dT.kind === 'V';
      for (let j = 0; j < Nc; j++) {
        rr[j] = heldT ? -Bt[kT][j] * scale[j] : 0;
        ri[j] = 0;
      }
      solve(rr, ri, yr, yi, heldT ? kT : -1, heldT ? 1 : 0);
      fl.forEach((k, a) => {
        for (let j = 0; j < Nc; j++) rr[j] = Bt[k][j] * scale[j];
        ri.fill(0);
        solve(rr, ri, Xr[a], Xi[a], k, -1);
      });
      // The floating terminals' circuit rows, a small complex system for their δV.
      const Ar = Array.from({ length: K }, () => new Float64Array(K + 1)), Ai = Array.from({ length: K }, () => new Float64Array(K + 1));
      fl.forEach((k, a) => {
        for (let b = 0; b < K; b++) {
          const [re, im] = dot(k, w, Xr[b], Xi[b]);
          Ar[a][b] = -re;
          Ai[a][b] = -im;
        }
        Ar[a][a] += DI[k];
        Ai[a][a] += w * DIs[k];
        // (another terminal's voltage, where k's current reads rows it reaches: a port's window
        // in a contact's end box)
        fl.forEach((k2, b) => {
          if (b === a) return;
          Ar[a][b] += DX[k][k2];
          Ai[a][b] += w * DXs[k][k2];
        });
        const [yre, yim] = dot(k, w, yr, yi);
        Ar[a][K] = (k === kT ? 1 : 0) - yre - (heldT ? DX[k][kT] : 0);
        Ai[a][K] = -yim - (heldT ? w * DXs[k][kT] : 0);
      });
      const dV = complexSolve(Ar, Ai, K);
      // δx = y − Σ X_k δV_k
      for (let j = 0; j < Nc; j++) {
        let re = yr[j], im = yi[j];
        for (let a = 0; a < K; a++) {
          re -= Xr[a][j] * dV.re[a] - Xi[a][j] * dV.im[a];
          im -= Xr[a][j] * dV.im[a] + Xi[a][j] * dV.re[a];
        }
        yr[j] = re;
        yi[j] = im;
        xr[this.fullOf[j]] = re;
        xi[this.fullOf[j]] = im;
      }
      // The current into terminal kT per volt, (C + iωC′)·δx + (∂I/∂V)(1 + iω′), with its error.
      // At the contact it's a sum whose terms can be far larger than itself (in a junction's
      // neutral ends, a conductance times a δη that cancels its neighbours' to 1e-17), and its
      // round-off reads as a constant conductance. Between the two contacts of a device without
      // ports, the total current (conduction and displacement) is the same through every cut, so
      // it's also read across each segment, and the reading with the least error taken. Each
      // reading's error is its terms' round-off and the solves' error in δx (to `accuracy` of its
      // largest entry) through its coefficients: a contact's or a cut's flux is a difference of
      // levels that the solved rows don't pin.
      // `parts`: the response as a sum, each part compact (r, i) and full-index (fr, fi), with the
      // absolute error of its entries.
      const read = (parts) => {
        const eps = 1.1e-16, c = Ct[kT], cs = Cs[kT];
        let Ir = DI[kT], Ii = w * DIs[kT], err = eps * (Math.abs(DI[kT]) + w * Math.abs(DIs[kT]));
        // (the floating terminals' voltages, where kT's current reads rows they reach)
        fl.forEach((k2, b) => {
          const xr = DX[kT][k2], xi = w * DXs[kT][k2];
          Ir += xr * dV.re[b] - xi * dV.im[b];
          Ii += xr * dV.im[b] + xi * dV.re[b];
        });
        for (const p of parts) {
          const [re, im] = dot(kT, w, p.r, p.i);
          Ir += re;
          Ii += im;
          for (let j = 0; j < Nc; j++) err += (Math.abs(c[j]) + w * Math.abs(cs[j])) * (eps * Math.hypot(p.r[j], p.i[j]) + p.floor);
        }
        if (!cuts) return { Ir, Ii, err };
        const sgn = kT === 0 ? 1 : -1;
        for (const cut of cuts) {
          let re = 0, im = 0, e = 0;
          for (const p of parts) {
            for (let a = 0; a < cut.idx.length; a++) {
              const x = cut.idx[a], gk = cut.g[a], dk = w * cut.d[a];
              // (g + iωd)(xr + i xi)
              re += gk * p.fr[x] - dk * p.fi[x];
              im += gk * p.fi[x] + dk * p.fr[x];
              e += (Math.abs(gk) + Math.abs(dk)) * (eps * Math.hypot(p.fr[x], p.fi[x]) + p.floor);
            }
          }
          if (e < err) [err, Ir, Ii] = [e, sgn * re, sgn * im];
        }
        return { Ir, Ii, err };
      };
      const largest = (r, i) => {
        let m = 0;
        for (let j = 0; j < r.length; j++) m = Math.max(m, Math.hypot(r[j], i[j]));
        return m;
      };
      const full = (r, i) => {
        const fr = new Float64Array(N), fi = new Float64Array(N);
        for (let j = 0; j < Nc; j++) {
          fr[this.fullOf[j]] = r[j];
          fi[this.fullOf[j]] = i[j];
        }
        return [fr, fi];
      };
      let Zr, Zi;
      if (heldT) {
        // δI per volt: (C + iωC′)·δx + (∂I/∂V)(1 + iω′) at the terminal…
        let { Ir, Ii, err } = read([{ r: yr, i: yi, fr: xr, fi: xi, floor: accuracy * largest(yr, yi) }]);
        // …and where the solves leave that unresolved (the estimate past 1% of the current), the
        // response is split: a uniform shift of each region (z_i s on every η_i, s on φ̂), which
        // changes no flux and no concentration, so that J·e vanishes inside each region and is
        // taken at its edges; and what's left, small, solved for to GMRES's floor, its error now
        // relative to itself. An electrolyte at its open circuit needs it: the response to its
        // terminal is nearly such a shift, and its current the slope of levels uniform to 1e-14,
        // more than a double can hold. (A port's exchange sees the shift: its window's rows are
        // kept with the edges'.) So does a MOS capacitor's back contact, its channel held by a port.
        if (K === 0 && cuts && err > 1e-2 * Math.hypot(Ir, Ii)) {
          const [er, ei] = this._regionShift(yr, yi);
          const tr = new Float64Array(Nc), ti = new Float64Array(Nc);
          jv(er, ei, tr, ti);
          const edge = this._regionEdgeRows();
          const ar = new Float64Array(Nc), ai = new Float64Array(Nc);
          for (let j = 0; j < Nc; j++) {
            ar[j] = -Bt[kT][j] - (edge[j] ? tr[j] : 0);
            ai[j] = edge[j] ? -ti[j] : 0;
          }
          const rhs = precondition(ar, ai);
          const dr = yr.map((v, j) => v - er[j]), di = yi.map((v, j) => v - ei[j]);
          const { residual } = gmres(op, rhs, [dr, di], { tol: 1e-11 });
          const [efr, efi] = full(er, ei), [dfr, dfi] = full(dr, di);
          const r2 = read([
            { r: er, i: ei, fr: efr, fi: efi, floor: 0 },
            { r: dr, i: di, fr: dfr, fi: dfi, floor: Math.max(1e-14, residual) * largest(dr, di) },
          ]);
          if (r2.err < err) ({ Ir, Ii, err } = r2);
          if (profiles) {
            for (let j = 0; j < N; j++) {
              xr[j] = efr[j] + dfr[j];
              xi[j] = efi[j] + dfi[j];
            }
          }
        }
        const d2 = Ir * Ir + Ii * Ii; // Z = 1/δI
        Zr = Ir / d2;
        Zi = -Ii / d2;
      } else {
        // δV per unit current
        const a = fl.indexOf(kT);
        Zr = dV.re[a];
        Zi = dV.im[a];
      }
      out.Z.re[q] = Zr;
      out.Z.im[q] = Zi;
      if (profiles) out.profiles.push(this._smallSignalProfiles(xr, xi));
    });
    return out;
  }

  // A uniform shift of each region, fitted to a small-signal response (compact, complex): φ̂ by
  // s_r and each species' η by z_i s_r, with s_r the region's mean δφ̂ (or, where φ̂ isn't an
  // unknown, its carriers' mean δη/z). Coverages and flux unknowns aren't shifted.
  _regionShift(yr, yi) {
    const { n, M, z, loc, model } = this, R = this.rix, { grid } = model;
    const er = new Float64Array(yr.length), ei = new Float64Array(yr.length);
    for (let r = 0; r < model.regions.length; r++) {
      let sr = 0, si = 0, k = 0;
      const g0 = grid.regionStart[r], g1 = grid.regionEnd[r];
      for (let g = g0; g <= g1; g++) {
        const o = this.blockOfNode[g] * M;
        if (loc[o] >= 0) {
          sr += yr[R[o]];
          si += yi[R[o]];
          k++;
        }
      }
      if (k === 0) {
        for (let g = g0; g <= g1; g++) {
          const o = this.blockOfNode[g] * M;
          for (let i = 0; i < n; i++) {
            if (z[i] === 0 || loc[o + 1 + i] < 0) continue;
            sr += yr[R[o + 1 + i]] / z[i];
            si += yi[R[o + 1 + i]] / z[i];
            k++;
          }
        }
      }
      if (k === 0) continue;
      sr /= k;
      si /= k;
      for (let g = g0; g <= g1; g++) {
        const o = this.blockOfNode[g] * M;
        if (loc[o] >= 0) {
          er[R[o]] = sr;
          ei[R[o]] = si;
        }
        for (let i = 0; i < n; i++) {
          if (loc[o + 1 + i] < 0) continue;
          er[R[o + 1 + i]] = z[i] * sr;
          ei[R[o + 1 + i]] = z[i] * si;
        }
      }
    }
    return [er, ei];
  }

  // Compact rows a region's uniform shift can touch: faces, the nodes beside them, the two ends
  // and the ports' windows (whatever holds a level from outside). Elsewhere inside a region,
  // every row is unchanged by it.
  _regionEdgeRows() {
    const { M, loc, nB } = this, R = this.rix;
    const edge = new Uint8Array(this.sys.size);
    const mark = (b) => {
      if (b < 0 || b >= nB) return;
      for (let r = 0; r < M; r++) if (loc[b * M + r] >= 0) edge[R[b * M + r]] = 1;
    };
    mark(0);
    mark(nB - 1);
    for (let f = 0; f < this.nFaces; f++) {
      const b = this.blockOfFace[f];
      mark(b - 1);
      mark(b);
      mark(b + 1);
    }
    for (const port of this.model.ports) for (const g of port.nodes) mark(this.blockOfNode[g]);
    return edge;
  }

  // The cuts where the impedance may read a small-signal current. Each is the linearised total
  // current toward +x through one segment s of a region, Σ (g + iω d)·δu over the unknowns idx:
  // the Scharfetter–Gummel fluxes' derivatives (as in _segmentsDilute, or _segmentConcentrated
  // for non-ideal statistics) and the displacement's. Without ports it's the same through every cut; with them, only between a
  // contact and the nearest port's window (see _impedance).
  _currentCuts() {
    const { model, n, M, u, uLo, c, z } = this, F = FARADAY;
    const { grid, materials, regions } = model;
    const cuts = [];
    for (let r = 0; r < regions.length; r++) {
      const reg = regions[r], mat = materials[reg.material];
      if (mat.conductor || reg.mixing > 0) continue;
      const vel = reg.velocity;
      for (let s = grid.regionStart[r]; s < grid.regionEnd[r]; s++) {
        const bL = s + r, bR = bL + 1, h = grid.segLength[s] / grid.segArea[s];
        const terms = new Map();
        const add = (x, g, d) => {
          const t = terms.get(x) || [0, 0];
          t[0] += g;
          t[1] += d;
          terms.set(x, t);
        };
        const phiL = u[bL * M], phiR = u[bR * M], hasPhi = !this.phiUndefined[s];
        for (let i = 0; i < n; i++) {
          if (z[i] === 0 || !mat.present[i] || mat.D[i] === 0) continue;
          if (!mat.ideal) {
            this._fluxConcentrated(s, bL, bR, mat, h, vel, i);
            for (let k = 0; k < M; k++) {
              if (k === 0 && !hasPhi) continue;
              if (this.jL[k] !== 0) add(bL * M + k, z[i] * F * this.jL[k], 0);
              if (this.jR[k] !== 0) add(bR * M + k, z[i] * F * this.jR[k], 0);
            }
            continue;
          }
          const zi = z[i], r1 = 1 + i, cL = c[s * n + i], g = mat.D[i] / h;
          const pe = (vel * h) / mat.D[i];
          const d = zi * (phiR - phiL) - pe;
          const sg = sgFlux(g, d, u[bR * M + r1] - u[bL * M + r1] + (uLo[bR * M + r1] - uLo[bL * M + r1]) - pe, cL, c[(s + 1) * n + i], this.sg);
          const dNdd = sg[1];
          add(bL * M + r1, zi * F * sg[2], 0);
          add(bR * M + r1, -zi * F * sg[3], 0);
          if (hasPhi) {
            add(bL * M, zi * F * (-zi * dNdd - zi * sg[0]), 0);
            add(bR * M, zi * F * zi * dNdd, 0);
          }
        }
        // The displacement toward +x, ε(φ_L − φ_R)/h (in φ̂ = φ/V_T).
        const k = hasPhi ? (mat.epsr * EPS0 * this.VT) / h : 0;
        if (k > 0) {
          add(bL * M, 0, k);
          add(bR * M, 0, -k);
        }
        if (terms.size === 0) continue;
        const idx = [...terms.keys()];
        cuts.push({ s, idx, g: idx.map((x) => terms.get(x)[0]), d: idx.map((x) => terms.get(x)[1]) });
      }
    }
    return cuts;
  }

  // Complex profiles of δφ (V), δμ̄ (J/mol) and δc (mol/m³) from a small-signal solution.
  _smallSignalProfiles(xr, xi) {
    const { n, M, z, VT, model } = this;
    const nN = this.nNodes, RT = model.RT;
    const pair = () => ({ re: new Float64Array(nN), im: new Float64Array(nN) });
    const phi = pair(), mu = {}, c = {};
    for (const sp of model.species) {
      mu[sp.name] = pair();
      c[sp.name] = pair();
    }
    for (let g = 0; g < nN; g++) {
      const b = this.blockOfNode[g];
      const undef = this.phiUndefined[g];
      phi.re[g] = undef ? NaN : VT * xr[b * M];
      phi.im[g] = undef ? NaN : VT * xi[b * M];
      for (let i = 0; i < n; i++) {
        const name = model.species[i].name;
        if (!this.present[g * n + i]) {
          mu[name].re[g] = mu[name].im[g] = c[name].re[g] = c[name].im[g] = NaN;
          continue;
        }
        mu[name].re[g] = RT * xr[b * M + 1 + i];
        mu[name].im[g] = RT * xi[b * M + 1 + i];
        if (this.nodeConductor[g] >= 0) {
          // a metal has a Fermi level but no concentration (as a solution reports it)
          c[name].re[g] = c[name].im[g] = NaN;
          continue;
        }
        // δc_i = Σ_j K_ij (δη_j − z_j δφ̂)
        let sr = 0, si = 0;
        for (let j = 0; j < n; j++) {
          const Kij = this.nodeIdeal[g] ? (i === j ? this.c[g * n + i] : 0) : this.K[g * n * n + i * n + j];
          if (Kij === 0) continue;
          sr += Kij * (xr[b * M + 1 + j] - z[j] * xr[b * M]);
          si += Kij * (xi[b * M + 1 + j] - z[j] * xi[b * M]);
        }
        c[name].re[g] = sr;
        c[name].im[g] = si;
      }
    }
    return { phi, mu, c };
  }
}
