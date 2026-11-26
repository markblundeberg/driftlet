// Linear stability of a converged steady state: the rightmost (slowest-decaying or growing)
// eigenvalues of the dynamic linearised DAE pencil −J v = λ M v about the device's current state,
// with conserved amounts deflated by projection, left-eigenvector error bars, difference-form
// residual verification, and an honest windowed verdict.
//
// CONTRACT
//
// Methods of the Solver (src/solver.js; its layout is in docs/internals.md), installed on its
// prototype like SmallSignal. Device.stability (src/index.js) calls stability() on the device's
// current, converged steady state without re-solving (preserving the branch).
//
// What it uses of the rest:
//   this._assembleSplit() (src/assemble.js)            steady Jacobian J (this.sys) and storage
//                                                      Jacobian M (this.sysMass), plus terminal
//                                                      border terms (termB, termC/Cs, termDI/DIs,
//                                                      termDX/DXs) on the dynamic rows (steady
//                                                      false).
//   this._assembleSplitDifference() (src/difference.js) cancellation-resistant J·v and M·v for
//                                                      held-vector residual checks.
//   this._currentCuts() (src/impedance.js)             cut-current readings for floating contacts
//                                                      so extreme face capacitances or metal G
//                                                      don't cancel the series capacitance in the
//                                                      Schur complement.
//   this._smallSignalProfiles() (src/impedance.js)     complex mode profiles (δφ, δμ̄, δc).
//   lawAmount, sharedChange (src/lawRows.js), chargeOnly (src/laws.js)
//                                                      exact gradients c = Mᵀw of conserved
//                                                      combinations, held-charge terminals and
//                                                      node-local invariants for deflation.
//   BlockTridiagonal, ComplexBlockTridiagonal (src/blockTridiagonal.js)
//                                                      forward and transposed (solveT / solveH)
//                                                      solves for shift-invert Arnoldi, inverse
//                                                      iteration, and left eigenvectors.
//   cdiv (src/impedance.js)                            overflow-safe complex division.
//
// stability(opts) → StabilityResult
//   opts:
//     shifts?: number[]         real positive shifts σ (1/s); default log-spaced over the device's
//                               own rates, never 0.
//     frequencies?: number[]    optional frequencies f (Hz) for imaginary-axis / complex shifts
//                               σ = σ_r + 2πi f (for Hopf pairs far from the real axis).
//     modes?: number            number of rightmost modes to return (default 8).
//     profiles?: boolean        attach complex mode shapes { phi, mu, c, termV } (default false).
//     kBar?: number             significance multiple of the error bar (default 10).
//     debug?: boolean           include raw { J, M, deflation, labels } as plain arrays for
//                               independent QZ cross-checks (default false).
//
// Precondition: the solver must hold a converged, un-regridded steady state (atSteady &&
// !regriddedSteady), else SolverError. Before running Arnoldi, one un-kicked dynamic step(h) is
// taken on a snapshot (immediately restored): if it fails or moves the state by more than 1e-5
// thermal units, stability() returns verdict 'undecided' with the measured drift rather than
// returning eigenvalues of a linearisation about a non-steady state.
//
// Deflation: every exact conserved combination in constraints (whose terminals are all driven by
// current), every current-driven chargeOnly terminal, and every localConstraint (per node) gives a
// gradient c_j = Mᵀw_j with w_jᵀJ = 0. The hyperplane c_j·x = 0 is invariant under T_σ =
// (J + σM)⁻¹M (since c_j·T_σ x = c_j·x / σ) and contains every finite eigenvector with λ ≠ 0.
// Arnoldi starts in this hyperplane (purified twice through T_σ to suppress index-1/2 DAE infinity
// directions) and re-projects at every step. Dependent gradients are dropped automatically.
//
// Refinement, residuals and error bars: each candidate eigenpair λ = σ − 1/μ is refined at a
// nearby shift λ_ref via two-sided inverse iteration (right v, left w via solveT/solveH, two-sided
// Rayleigh quotient). Held-vector residual r = (J + λM)v is evaluated in difference form; modes
// with row-scaled residual > 1e-5 or non-finite error bar (catching spurious infinite DAE roots)
// are rejected. Each accepted mode reports its first-order componentwise error bar:
//   error = ε · |w|ᴴ (|J| + |λ| |M|) |v| / |wᴴ M v|.
//
// Verdict (never a bare 'stable'):
//   'unstable'                 a verified mode has Re λ > kBar · error;
//   'neutral within resolution' no growing mode, and a verified non-deflated mode has
//                              |Re λ| ≤ kBar · error;
//   'no growing mode found'    all verified modes have Re λ < −kBar · error, qualified by the
//                              searched window { rLo, rHi, omegaHi };
//   'undecided'                precondition failed (drifting state) or no shift factored cleanly.
//
// Invariants:
//   1. Completely non-destructive: STATE (u, uLo, termV, time, history, atSteady, stepped, solvedV,
//      ...), LEARNED, allTerminals, steady and lawForm are identical before and after the call,
//      even on error or undecided.
//   2. Dynamic rows only: assembly always runs with this.steady = false; steady law pins, flat
//      rows, affinities and island sums never enter J or M.
//   3. Deterministic: fixed start-vector PRNG seed; same state and options give bit-identical
//      results.
//   4. Pure JS, zero dependencies, block-tridiagonal + tiny K×K Schur border only, no ** or
//      Math.pow.
//
// Tests: test/stability.test.js (Hodgkin–Huxley Hopf vs advance(), RC + Sol's C_face=1e16
// cancellation cut, blocked diffusion discrete cosine spectrum with n=0 deflation, immobile/trap
// local constraints, Fe³⁺/Fe²⁺ d26594 falling vs rising branch vs kicked advance(), non-steady
// refusal, state preservation, and BlockTridiagonal solveT/solveH).

import { BlockTridiagonal, ComplexBlockTridiagonal } from './blockTridiagonal.js';
import { FARADAY } from './constants.js';
import { SolverError } from './errors.js';
import { cdiv } from './impedance.js';
import { chargeOnly } from './laws.js';
import { lawAmount, sharedChange } from './lawRows.js';

const EPS = Number.EPSILON; // 2.220446049250313e-16

// --- Small dense real/complex linear solvers (for K×K bordered Schur complements) ---

function realLuFactor(A, K) {
  const lu = A.map((row) => Float64Array.from(row));
  const piv = new Int32Array(K);
  for (let k = 0; k < K; k++) {
    let p = k, best = Math.abs(lu[k][k]);
    for (let r = k + 1; r < K; r++) {
      const v = Math.abs(lu[r][k]);
      if (v > best) { best = v; p = r; }
    }
    piv[k] = p;
    if (best === 0 || !Number.isFinite(best)) return null;
    if (p !== k) { const t = lu[k]; lu[k] = lu[p]; lu[p] = t; }
    const inv = 1 / lu[k][k];
    for (let r = k + 1; r < K; r++) {
      const f = (lu[r][k] *= inv);
      if (f !== 0) for (let c = k + 1; c < K; c++) lu[r][c] -= f * lu[k][c];
    }
  }
  return { lu, piv, K };
}

function realLuSolve(fac, rhs) {
  const { lu, piv, K } = fac, x = Float64Array.from(rhs);
  for (let k = 0; k < K; k++) {
    const p = piv[k];
    if (p !== k) { const t = x[k]; x[k] = x[p]; x[p] = t; }
  }
  for (let r = 1; r < K; r++) {
    let s = x[r];
    for (let c = 0; c < r; c++) s -= lu[r][c] * x[c];
    x[r] = s;
  }
  for (let r = K - 1; r >= 0; r--) {
    let s = x[r];
    for (let c = r + 1; c < K; c++) s -= lu[r][c] * x[c];
    x[r] = s / lu[r][r];
  }
  return x;
}

function realLuSolveT(fac, rhs) {
  const { lu, piv, K } = fac, x = Float64Array.from(rhs);
  for (let r = 0; r < K; r++) {
    let s = x[r];
    for (let c = 0; c < r; c++) s -= lu[c][r] * x[c];
    x[r] = s / lu[r][r];
  }
  for (let r = K - 2; r >= 0; r--) {
    let s = x[r];
    for (let c = r + 1; c < K; c++) s -= lu[c][r] * x[c];
    x[r] = s;
  }
  for (let k = K - 1; k >= 0; k--) {
    const p = piv[k];
    if (p !== k) { const t = x[k]; x[k] = x[p]; x[p] = t; }
  }
  return x;
}

function complexLuFactor(Ar, Ai, K) {
  const lr = Ar.map((row) => Float64Array.from(row));
  const li = Ai.map((row) => Float64Array.from(row));
  const piv = new Int32Array(K);
  for (let k = 0; k < K; k++) {
    let p = k, best = Math.hypot(lr[k][k], li[k][k]);
    for (let r = k + 1; r < K; r++) {
      const v = Math.hypot(lr[r][k], li[r][k]);
      if (v > best) { best = v; p = r; }
    }
    piv[k] = p;
    if (best === 0 || !Number.isFinite(best)) return null;
    if (p !== k) {
      let t = lr[k]; lr[k] = lr[p]; lr[p] = t;
      t = li[k]; li[k] = li[p]; li[p] = t;
    }
    const dr = lr[k][k], di = li[k][k];
    for (let r = k + 1; r < K; r++) {
      const [fr, fi] = cdiv(lr[r][k], li[r][k], dr, di);
      lr[r][k] = fr;
      li[r][k] = fi;
      if (fr === 0 && fi === 0) continue;
      for (let c = k + 1; c < K; c++) {
        const ur = lr[k][c], ui = li[k][c];
        lr[r][c] -= fr * ur - fi * ui;
        li[r][c] -= fr * ui + fi * ur;
      }
    }
  }
  return { lr, li, piv, K };
}

function complexLuSolve(fac, br, bi) {
  const { lr, li, piv, K } = fac;
  const xr = Float64Array.from(br), xi = Float64Array.from(bi);
  for (let k = 0; k < K; k++) {
    const p = piv[k];
    if (p !== k) {
      let t = xr[k]; xr[k] = xr[p]; xr[p] = t;
      t = xi[k]; xi[k] = xi[p]; xi[p] = t;
    }
  }
  for (let r = 1; r < K; r++) {
    let vr = xr[r], vi = xi[r];
    for (let c = 0; c < r; c++) {
      vr -= lr[r][c] * xr[c] - li[r][c] * xi[c];
      vi -= lr[r][c] * xi[c] + li[r][c] * xr[c];
    }
    xr[r] = vr;
    xi[r] = vi;
  }
  for (let r = K - 1; r >= 0; r--) {
    let vr = xr[r], vi = xi[r];
    for (let c = r + 1; c < K; c++) {
      vr -= lr[r][c] * xr[c] - li[r][c] * xi[c];
      vi -= lr[r][c] * xi[c] + li[r][c] * xr[c];
    }
    [xr[r], xi[r]] = cdiv(vr, vi, lr[r][r], li[r][r]);
  }
  return [xr, xi];
}

function complexLuSolveH(fac, br, bi) {
  const { lr, li, piv, K } = fac;
  // Solve Aᴴ x = b ⇔ Aᵀ conj(x) = conj(b)
  const xr = Float64Array.from(br), xi = Float64Array.from(bi, (v) => -v);
  for (let r = 0; r < K; r++) {
    let vr = xr[r], vi = xi[r];
    for (let c = 0; c < r; c++) {
      vr -= lr[c][r] * xr[c] - li[c][r] * xi[c];
      vi -= lr[c][r] * xi[c] + li[c][r] * xr[c];
    }
    [xr[r], xi[r]] = cdiv(vr, vi, lr[r][r], li[r][r]);
  }
  for (let r = K - 2; r >= 0; r--) {
    let vr = xr[r], vi = xi[r];
    for (let c = r + 1; c < K; c++) {
      vr -= lr[c][r] * xr[c] - li[c][r] * xi[c];
      vi -= lr[c][r] * xi[c] + li[c][r] * xr[c];
    }
    xr[r] = vr;
    xi[r] = vi;
  }
  for (let k = K - 1; k >= 0; k--) {
    const p = piv[k];
    if (p !== k) {
      let t = xr[k]; xr[k] = xr[p]; xr[p] = t;
      t = xi[k]; xi[k] = xi[p]; xi[p] = t;
    }
  }
  for (let k = 0; k < K; k++) xi[k] = -xi[k];
  return [xr, xi];
}

// --- Small dense real upper-Hessenberg Francis double-shift QR eigensolver (m ≤ 40) ---

function realHessenbergEig(Hin, m) {
  const H = Array.from({ length: m }, (_, i) => Float64Array.from(Hin[i].subarray(0, m)));
  const Q = Array.from({ length: m }, (_, i) => {
    const row = new Float64Array(m);
    row[i] = 1;
    return row;
  });
  let normH = 0;
  for (let i = 0; i < m; i++) for (let j = Math.max(0, i - 1); j < m; j++) normH += Math.abs(H[i][j]);

  let nn = m - 1, iter = 0, totalIter = 0;
  const maxTotal = 60 * m;
  while (nn >= 0 && totalIter < maxTotal) {
    let l = nn;
    while (l > 0) {
      const s = Math.abs(H[l - 1][l - 1]) + Math.abs(H[l][l]);
      const tol = EPS * (s > 0 ? s : normH || 1);
      if (Math.abs(H[l][l - 1]) <= tol) { H[l][l - 1] = 0; break; }
      l--;
    }
    if (l === nn) {
      nn--;
      iter = 0;
      continue;
    }
    if (l === nn - 1) {
      // Standardize 2×2 block: if real eigenvalues, reduce to triangular via Givens rotation.
      const a = H[l][l], b = H[l][nn], c = H[nn][l], d = H[nn][nn];
      const tr = 0.5 * (a + d), df = 0.5 * (a - d), disc = df * df + b * c;
      if (disc >= 0) {
        const z = df >= 0 ? Math.sqrt(disc) : -Math.sqrt(disc);
        const p = df + z;
        // Eigenvector of (A - (tr - z)I) is [p, c] or [b, -p]
        let cs, sn;
        if (Math.abs(p) + Math.abs(c) >= Math.abs(b) + Math.abs(p)) {
          const r = Math.hypot(p, c);
          cs = r > 0 ? p / r : 1;
          sn = r > 0 ? c / r : 0;
        } else {
          const r = Math.hypot(b, p);
          cs = r > 0 ? b / r : 1;
          sn = r > 0 ? -p / r : 0;
        }
        // Apply similarity Gᵀ H G with G = [[cs, -sn], [sn, cs]]
        for (let j = l; j < m; j++) {
          const t1 = H[l][j], t2 = H[nn][j];
          H[l][j] = cs * t1 + sn * t2;
          H[nn][j] = -sn * t1 + cs * t2;
        }
        for (let i = 0; i <= nn; i++) {
          const t1 = H[i][l], t2 = H[i][nn];
          H[i][l] = cs * t1 + sn * t2;
          H[i][nn] = -sn * t1 + cs * t2;
        }
        for (let i = 0; i < m; i++) {
          const t1 = Q[i][l], t2 = Q[i][nn];
          Q[i][l] = cs * t1 + sn * t2;
          Q[i][nn] = -sn * t1 + cs * t2;
        }
        H[nn][l] = 0;
      }
      nn -= 2;
      iter = 0;
      continue;
    }
    iter++;
    totalIter++;
    // Francis double shift (or exceptional shift every 10 iters).
    let s, t;
    if (iter % 10 === 0) {
      const w = Math.abs(H[nn][nn - 1]) + Math.abs(H[nn - 1][nn - 2]);
      s = 1.5 * w;
      t = 0.5625 * w * w;
    } else {
      const a = H[nn - 1][nn - 1], b = H[nn - 1][nn], c = H[nn][nn - 1], d = H[nn][nn];
      s = a + d;
      t = a * d - b * c;
    }
    // Look for two consecutive small subdiagonal elements to start bulge chase.
    let k0 = nn - 2;
    let x = 0, y = 0, z = 0;
    for (; k0 >= l; k0--) {
      const h00 = H[k0][k0], h11 = H[k0 + 1][k0 + 1], h10 = H[k0 + 1][k0], h01 = H[k0][k0 + 1];
      x = h00 * h00 - s * h00 + t + h01 * h10;
      y = h10 * (h00 + h11 - s);
      z = h10 * H[k0 + 2][k0 + 1];
      if (k0 === l) break;
      const lhs = Math.abs(H[k0][k0 - 1]) * (Math.abs(y) + Math.abs(z));
      const rhs = EPS * Math.abs(x) * (Math.abs(H[k0 - 1][k0 - 1]) + Math.abs(h00) + Math.abs(h11));
      if (lhs <= rhs) break;
    }
    for (let k = k0; k <= nn - 1; k++) {
      const notLast = k < nn - 1;
      if (k > k0) {
        x = H[k][k - 1];
        y = H[k + 1][k - 1];
        z = notLast ? H[k + 2][k - 1] : 0;
      }
      const sc = Math.abs(x) + Math.abs(y) + Math.abs(z);
      if (sc === 0) continue;
      x /= sc; y /= sc; z /= sc;
      let r = Math.sqrt(x * x + y * y + z * z);
      if (x < 0) r = -r;
      if (r === 0) continue;
      if (k > k0) H[k][k - 1] = -r * sc;
      else if (k0 > l) H[k][k - 1] = -H[k][k - 1];
      x += r;
      const p = x / r, q = y / r, u = z / r;
      const wy = y / x, wz = z / x;
      // Row modification: H[k..k+2, k..m-1]
      for (let j = k; j < m; j++) {
        let v = H[k][j] + wy * H[k + 1][j];
        if (notLast) {
          v += wz * H[k + 2][j];
          H[k + 2][j] -= v * u;
        }
        H[k + 1][j] -= v * q;
        H[k][j] -= v * p;
      }
      // Column modification: H[0..min(nn, k+3), k..k+2]
      const iMax = Math.min(nn, k + 3);
      for (let i = 0; i <= iMax; i++) {
        let v = p * H[i][k] + q * H[i][k + 1];
        if (notLast) {
          v += u * H[i][k + 2];
          H[i][k + 2] -= v * wz;
        }
        H[i][k + 1] -= v * wy;
        H[i][k] -= v;
      }
      // Accumulate Q
      for (let i = 0; i < m; i++) {
        let v = p * Q[i][k] + q * Q[i][k + 1];
        if (notLast) {
          v += u * Q[i][k + 2];
          Q[i][k + 2] -= v * wz;
        }
        Q[i][k + 1] -= v * wy;
        Q[i][k] -= v;
      }
      if (k > k0) {
        H[k + 1][k - 1] = 0;
        if (notLast) H[k + 2][k - 1] = 0;
      }
    }
  }

  // Extract eigenvalues and eigenvectors of real Schur matrix H, then lift by Q.
  const eigs = [];
  for (let i = 0; i < m; ) {
    if (i + 1 < m && H[i + 1][i] !== 0) {
      const a = H[i][i], b = H[i][i + 1], c = H[i + 1][i], d = H[i + 1][i + 1];
      const re = 0.5 * (a + d), df = 0.5 * (a - d);
      const disc = -(df * df + b * c);
      const im = disc > 0 ? Math.sqrt(disc) : 0;
      // Solve (H - (re + i*im)I) y = 0 by complex back-substitution for upper-quasi-triangular H.
      const yr = new Float64Array(m), yi = new Float64Array(m);
      if (Math.abs(b) >= Math.abs(c)) {
        yr[i] = b; yi[i] = 0;
        yr[i + 1] = re - a; yi[i + 1] = im;
      } else {
        yr[i] = re - d; yi[i] = im;
        yr[i + 1] = c; yi[i + 1] = 0;
      }
      for (let r = i - 1; r >= 0; ) {
        if (r > 0 && H[r][r - 1] !== 0) {
          // 2×2 block at (r-1, r)
          let rhs1r = 0, rhs1i = 0, rhs2r = 0, rhs2i = 0;
          for (let j = r + 1; j <= i + 1; j++) {
            rhs1r -= H[r - 1][j] * yr[j]; rhs1i -= H[r - 1][j] * yi[j];
            rhs2r -= H[r][j] * yr[j]; rhs2i -= H[r][j] * yi[j];
          }
          const a11r = H[r - 1][r - 1] - re, a11i = -im, a12 = H[r - 1][r];
          const a21 = H[r][r - 1], a22r = H[r][r] - re, a22i = -im;
          const detr = a11r * a22r - a11i * a22i - a12 * a21;
          const deti = a11r * a22i + a11i * a22r;
          const n1r = rhs1r * a22r - rhs1i * a22i - a12 * rhs2r;
          const n1i = rhs1r * a22i + rhs1i * a22r - a12 * rhs2i;
          const n2r = a11r * rhs2r - a11i * rhs2i - a21 * rhs1r;
          const n2i = a11r * rhs2i + a11i * rhs2r - a21 * rhs1i;
          const dm = Math.hypot(detr, deti);
          if (dm > 0) {
            [yr[r - 1], yi[r - 1]] = cdiv(n1r, n1i, detr, deti);
            [yr[r], yi[r]] = cdiv(n2r, n2i, detr, deti);
          }
          r -= 2;
        } else {
          let sr = 0, si = 0;
          for (let j = r + 1; j <= i + 1; j++) {
            sr -= H[r][j] * yr[j];
            si -= H[r][j] * yi[j];
          }
          const dr = H[r][r] - re, di = -im;
          if (Math.hypot(dr, di) > 0) [yr[r], yi[r]] = cdiv(sr, si, dr, di);
          r--;
        }
      }
      const zr = new Float64Array(m), zi = new Float64Array(m);
      let n2 = 0;
      for (let r = 0; r < m; r++) {
        let ar = 0, ai = 0;
        for (let j = 0; j <= i + 1; j++) { ar += Q[r][j] * yr[j]; ai += Q[r][j] * yi[j]; }
        zr[r] = ar; zi[r] = ai;
        n2 += ar * ar + ai * ai;
      }
      const inv = n2 > 0 ? 1 / Math.sqrt(n2) : 1;
      for (let r = 0; r < m; r++) { zr[r] *= inv; zi[r] *= inv; }
      eigs.push({ re, im, zr, zi });
      const zic = Float64Array.from(zi, (v) => -v);
      eigs.push({ re, im: -im, zr: Float64Array.from(zr), zi: zic });
      i += 2;
    } else {
      const re = H[i][i];
      const y = new Float64Array(m);
      y[i] = 1;
      for (let r = i - 1; r >= 0; ) {
        if (r > 0 && H[r][r - 1] !== 0) {
          let rhs1 = 0, rhs2 = 0;
          for (let j = r + 1; j <= i; j++) { rhs1 -= H[r - 1][j] * y[j]; rhs2 -= H[r][j] * y[j]; }
          const a11 = H[r - 1][r - 1] - re, a12 = H[r - 1][r], a21 = H[r][r - 1], a22 = H[r][r] - re;
          const det = a11 * a22 - a12 * a21;
          if (Math.abs(det) > 0) {
            y[r - 1] = (rhs1 * a22 - a12 * rhs2) / det;
            y[r] = (a11 * rhs2 - a21 * rhs1) / det;
          }
          r -= 2;
        } else {
          let s = 0;
          for (let j = r + 1; j <= i; j++) s -= H[r][j] * y[j];
          const d = H[r][r] - re;
          y[r] = Math.abs(d) > 0 ? s / d : 0;
          r--;
        }
      }
      const zr = new Float64Array(m), zi = new Float64Array(m);
      let n2 = 0;
      for (let r = 0; r < m; r++) {
        let a = 0;
        for (let j = 0; j <= i; j++) a += Q[r][j] * y[j];
        zr[r] = a;
        n2 += a * a;
      }
      const inv = n2 > 0 ? 1 / Math.sqrt(n2) : 1;
      for (let r = 0; r < m; r++) zr[r] *= inv;
      eigs.push({ re, im: 0, zr, zi });
      i++;
    }
  }
  return eigs;
}

// --- Small dense complex upper-Hessenberg Wilkinson QR eigensolver (for complex shifts) ---

function complexHessenbergEig(Hrin, Hiin, m) {
  const Hr = Array.from({ length: m }, (_, i) => Float64Array.from(Hrin[i].subarray(0, m)));
  const Hi = Array.from({ length: m }, (_, i) => Float64Array.from(Hiin[i].subarray(0, m)));
  const Qr = Array.from({ length: m }, (_, i) => { const r = new Float64Array(m); r[i] = 1; return r; });
  const Qi = Array.from({ length: m }, () => new Float64Array(m));
  let normH = 0;
  for (let i = 0; i < m; i++) for (let j = Math.max(0, i - 1); j < m; j++) normH += Math.hypot(Hr[i][j], Hi[i][j]);

  let nn = m - 1, iter = 0, totalIter = 0;
  const maxTotal = 60 * m;
  while (nn >= 0 && totalIter < maxTotal) {
    let l = nn;
    while (l > 0) {
      const s = Math.hypot(Hr[l - 1][l - 1], Hi[l - 1][l - 1]) + Math.hypot(Hr[l][l], Hi[l][l]);
      if (Math.hypot(Hr[l][l - 1], Hi[l][l - 1]) <= EPS * (s > 0 ? s : normH || 1)) {
        Hr[l][l - 1] = 0; Hi[l][l - 1] = 0;
        break;
      }
      l--;
    }
    if (l === nn) { nn--; iter = 0; continue; }
    iter++;
    totalIter++;
    // Wilkinson shift from trailing 2×2 block
    let sr, si;
    if (iter % 10 === 0) {
      const w = Math.hypot(Hr[nn][nn - 1], Hi[nn][nn - 1]);
      sr = Hr[nn][nn] + 0.75 * w;
      si = Hi[nn][nn] + 0.4375 * w;
    } else {
      const a00r = Hr[nn - 1][nn - 1], a00i = Hi[nn - 1][nn - 1];
      const a01r = Hr[nn - 1][nn], a01i = Hi[nn - 1][nn];
      const a10r = Hr[nn][nn - 1], a10i = Hi[nn][nn - 1];
      const a11r = Hr[nn][nn], a11i = Hi[nn][nn];
      const trr = 0.5 * (a00r + a11r), tri = 0.5 * (a00i + a11i);
      const dfr = 0.5 * (a00r - a11r), dfi = 0.5 * (a00i - a11i);
      const discr = dfr * dfr - dfi * dfi + (a01r * a10r - a01i * a10i);
      const disci = 2 * dfr * dfi + (a01r * a10i + a01i * a10r);
      const dm = Math.hypot(discr, disci);
      const sqr = Math.sqrt(0.5 * (dm + Math.abs(discr)));
      const sqi = sqr > 0 ? Math.abs(disci) / (2 * sqr) : 0;
      let rr, ri;
      if (discr >= 0) { rr = sqr; ri = disci >= 0 ? sqi : -sqi; }
      else { rr = sqi; ri = disci >= 0 ? sqr : -sqr; }
      const l1r = trr + rr, l1i = tri + ri, l2r = trr - rr, l2i = tri - ri;
      if (Math.hypot(l1r - a11r, l1i - a11i) <= Math.hypot(l2r - a11r, l2i - a11i)) { sr = l1r; si = l1i; }
      else { sr = l2r; si = l2i; }
    }
    // Single-shift QR step via Givens rotations on H[l..nn, l..nn]
    for (let i = l; i <= nn; i++) { Hr[i][i] -= sr; Hi[i][i] -= si; }
    const cs = new Float64Array(nn - l), snr = new Float64Array(nn - l), sni = new Float64Array(nn - l);
    for (let k = l; k < nn; k++) {
      const ar = Hr[k][k], ai = Hi[k][k], br = Hr[k + 1][k], bi = Hi[k + 1][k];
      const am = Math.hypot(ar, ai), bm = Math.hypot(br, bi), r = Math.hypot(am, bm);
      if (r === 0) { cs[k - l] = 1; snr[k - l] = 0; sni[k - l] = 0; continue; }
      const c = am / r;
      let svr, svi;
      if (am > 0) {
        // s = (a / |a|) * conj(b) / r so that c * a + s * b = (a/|a|) * r and -conj(s)*a + c*b = 0
        const er = ar / am, ei = ai / am;
        svr = (er * br + ei * bi) / r;
        svi = (ei * br - er * bi) / r;
      } else {
        svr = 1; svi = 0;
      }
      cs[k - l] = c; snr[k - l] = svr; sni[k - l] = svi;
      // Left multiply rows k, k+1 by G = [[c, s], [-conj(s), c]]
      for (let j = k; j < m; j++) {
        const xr = Hr[k][j], xi = Hi[k][j], yr = Hr[k + 1][j], yi = Hi[k + 1][j];
        Hr[k][j] = c * xr + (svr * yr - svi * yi);
        Hi[k][j] = c * xi + (svr * yi + svi * yr);
        Hr[k + 1][j] = -(svr * xr + svi * xi) + c * yr;
        Hi[k + 1][j] = -(svr * xi - svi * xr) + c * yi;
      }
      Hr[k + 1][k] = 0; Hi[k + 1][k] = 0;
    }
    // Right multiply columns k, k+1 by Gᴴ = [[c, -s], [conj(s), c]]
    for (let k = l; k < nn; k++) {
      const c = cs[k - l], svr = snr[k - l], svi = sni[k - l];
      const iMax = Math.min(nn, k + 1);
      for (let i = 0; i <= iMax; i++) {
        const xr = Hr[i][k], xi = Hi[i][k], yr = Hr[i][k + 1], yi = Hi[i][k + 1];
        Hr[i][k] = c * xr + (svr * yr + svi * yi);
        Hi[i][k] = c * xi + (svr * yi - svi * yr);
        Hr[i][k + 1] = -(svr * xr - svi * xi) + c * yr;
        Hi[i][k + 1] = -(svr * xi + svi * xr) + c * yi;
      }
      for (let i = 0; i < m; i++) {
        const xr = Qr[i][k], xi = Qi[i][k], yr = Qr[i][k + 1], yi = Qi[i][k + 1];
        Qr[i][k] = c * xr + (svr * yr + svi * yi);
        Qi[i][k] = c * xi + (svr * yi - svi * yr);
        Qr[i][k + 1] = -(svr * xr - svi * xi) + c * yr;
        Qi[i][k + 1] = -(svr * xi + svi * xr) + c * yi;
      }
    }
    for (let i = l; i <= nn; i++) { Hr[i][i] += sr; Hi[i][i] += si; }
  }

  const eigs = [];
  for (let i = 0; i < m; i++) {
    const re = Hr[i][i], im = Hi[i][i];
    const yr = new Float64Array(m), yi = new Float64Array(m);
    yr[i] = 1;
    for (let r = i - 1; r >= 0; r--) {
      let sr = 0, si = 0;
      for (let j = r + 1; j <= i; j++) {
        sr -= Hr[r][j] * yr[j] - Hi[r][j] * yi[j];
        si -= Hr[r][j] * yi[j] + Hi[r][j] * yr[j];
      }
      const dr = Hr[r][r] - re, di = Hi[r][r] - im;
      if (Math.hypot(dr, di) > 0) [yr[r], yi[r]] = cdiv(sr, si, dr, di);
    }
    const zr = new Float64Array(m), zi = new Float64Array(m);
    let n2 = 0;
    for (let r = 0; r < m; r++) {
      let ar = 0, ai = 0;
      for (let j = 0; j <= i; j++) {
        ar += Qr[r][j] * yr[j] - Qi[r][j] * yi[j];
        ai += Qr[r][j] * yi[j] + Qi[r][j] * yr[j];
      }
      zr[r] = ar; zi[r] = ai;
      n2 += ar * ar + ai * ai;
    }
    const inv = n2 > 0 ? 1 / Math.sqrt(n2) : 1;
    for (let r = 0; r < m; r++) { zr[r] *= inv; zi[r] *= inv; }
    eigs.push({ re, im, zr, zi });
  }
  return eigs;
}

// --- Dynamic bordered pencil + Sol's cut-current readings ---

function btMul(bt, x, out) {
  const { n, A, B, C, sizes, offA, offB, offC, offX } = bt;
  for (let blk = 0; blk < n; blk++) {
    const m = sizes[blk], v = offX[blk];
    const mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < n - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      let s = 0;
      for (let c = 0; c < mp; c++) s += A[offA[blk] + r * mp + c] * x[v - mp + c];
      for (let c = 0; c < m; c++) s += B[offB[blk] + r * m + c] * x[v + c];
      for (let c = 0; c < mn; c++) s += C[offC[blk] + r * mn + c] * x[v + m + c];
      out[v + r] = s;
    }
  }
  return out;
}

function btMulT(bt, x, out) {
  out.fill(0);
  const { n, A, B, C, sizes, offA, offB, offC, offX } = bt;
  for (let blk = 0; blk < n; blk++) {
    const m = sizes[blk], v = offX[blk];
    const mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < n - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      const xr = x[v + r];
      if (xr === 0) continue;
      for (let c = 0; c < mp; c++) out[v - mp + c] += A[offA[blk] + r * mp + c] * xr;
      for (let c = 0; c < m; c++) out[v + c] += B[offB[blk] + r * m + c] * xr;
      for (let c = 0; c < mn; c++) out[v + m + c] += C[offC[blk] + r * mn + c] * xr;
    }
  }
  return out;
}

function btAbsMul(bt, x, out, factor = 1) {
  const { n, A, B, C, sizes, offA, offB, offC, offX } = bt;
  for (let blk = 0; blk < n; blk++) {
    const m = sizes[blk], v = offX[blk];
    const mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < n - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      let s = 0;
      for (let c = 0; c < mp; c++) s += Math.abs(A[offA[blk] + r * mp + c]) * x[v - mp + c];
      for (let c = 0; c < m; c++) s += Math.abs(B[offB[blk] + r * m + c]) * x[v + c];
      for (let c = 0; c < mn; c++) s += Math.abs(C[offC[blk] + r * mn + c]) * x[v + m + c];
      out[v + r] += factor * s;
    }
  }
}

function denseOfBt(bt, K) {
  const { n, A, B, C, sizes, offA, offB, offC, offX } = bt, Nc = bt.size, N = Nc + K;
  const out = Array.from({ length: N }, () => new Array(N).fill(0));
  for (let blk = 0; blk < n; blk++) {
    const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < n - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      const row = out[offX[blk] + r];
      for (let c = 0; c < mp; c++) row[offX[blk - 1] + c] = A[offA[blk] + r * mp + c];
      for (let c = 0; c < m; c++) row[offX[blk] + c] = B[offB[blk] + r * m + c];
      for (let c = 0; c < mn; c++) row[offX[blk + 1] + c] = C[offC[blk] + r * mn + c];
    }
  }
  return out;
}

function buildPencil(s) {
  s.sourceTime = s.time;
  s.computeConcentrations();
  s.cOld.set(s.c);
  s.thOld.set(s.th);
  s.th0Old.set(s.th0);
  s.gamOld.set(s.gam);
  s.gateList.forEach((gt, j) => (s.gateOld[j] = s.u[gt.o] + s.uLo[gt.o]));

  const byChange = (cs) => cs.circuit >= 0 && !(cs.contacts?.length > 0) && Number.isFinite(cs.lambda) && cs.lambda !== 0 && !sharedChange(s, cs);
  s.lawForm = new Map([...s.constraints, ...s.timeLaws].map((cs) => [cs, byChange(cs)]));
  const lawRead = new Set([...s.constraints, ...s.timeLaws].filter(byChange).map((cs) => cs.circuit));

  const { rest: Jr, massRest: Sr, lin } = s._assembleSplitDifference();
  const { sys: Jg, mass: Mg } = s._assembleSplit();
  const Nc = Jg.size, fl = s.floating.slice(), K = fl.length, N = Nc + K;

  // Keep raw contact storage rows before cut replacement (for contact-charge deflation terms).
  const rawCs = s.termCs.map((v) => v.slice(0, Nc));
  const rawDIs = Float64Array.from(s.termDIs);
  const rawDXs = s.termDXs.map((v) => Float64Array.from(v));

  const B = fl.map((k) => s.termB[k].slice(0, Nc));
  const C = fl.map((k) => s.termC[k].slice(0, Nc));
  const Cs = fl.map((k) => s.termCs[k].slice(0, Nc));
  const D = fl.map((k, a) => Float64Array.from(fl, (k2, b) => (a === b ? s.termDI[k] : s.termDX[k][k2])));
  const Ds = fl.map((k, a) => Float64Array.from(fl, (k2, b) => (a === b ? s.termDIs[k] : s.termDXs[k][k2])));

  // Sol's cut-current reading for floating contacts not read by a circuit law's change:
  // reading total current at an interior dielectric/electrolyte segment avoids Schur cancellation
  // when C_face or metal conductance G is huge.
  const { species } = s.model;
  const charged = (p) => p.capacitance || p.reactions.length > 0 || species.some((_, i) => s.z[i] !== 0 && p.species[i].type !== 'blocked');
  const allCuts = s._currentCuts();
  const windows = s.model.ports.filter(charged).flatMap((p) => Array.from(p.nodes));
  const lo = windows.length ? Math.min(...windows) : Infinity;
  const hi = windows.length ? Math.max(...windows) : -Infinity;
  const cutsOf = (k) => allCuts.filter((cut) => (k === 0 ? cut.s + 1 <= lo : cut.s >= hi));

  const rowMax = (c, cs, d, ds) => {
    let mx = 0;
    for (let j = 0; j < c.length; j++) { const v = Math.abs(c[j]); if (v > mx) mx = v; }
    for (let j = 0; j < cs.length; j++) { const v = Math.abs(cs[j]); if (v > mx) mx = v; }
    for (let j = 0; j < d.length; j++) { const v = Math.abs(d[j]); if (v > mx) mx = v; }
    for (let j = 0; j < ds.length; j++) { const v = Math.abs(ds[j]); if (v > mx) mx = v; }
    return mx;
  };
  const cutUsed = new Array(K).fill(-1);
  fl.forEach((k, a) => {
    if (k > 1 || lawRead.has(k)) return;
    const ownMax = rowMax(C[a], Cs[a], D[a], Ds[a]);
    const sgn = k === 0 ? 1 : -1, R = s.terms[k].drive.R;
    let bestMax = ownMax, best = null;
    for (const cut of cutsOf(k)) {
      const cc = new Float64Array(Nc), ccs = new Float64Array(Nc);
      let hasTerm = false;
      cut.idx.forEach((x, idxA) => {
        const r = s.rix[x];
        if (r >= Nc) return;
        const g = sgn * cut.g[idxA], d = sgn * cut.d[idxA];
        if (g !== 0 || d !== 0) hasTerm = true;
        cc[r] += g;
        ccs[r] += d;
      });
      if (!hasTerm) continue;
      const dd = new Float64Array(K), dds = new Float64Array(K);
      if (R > 0) dd[a] = 1 / R;
      const cm = rowMax(cc, ccs, dd, dds);
      if (cm < bestMax) {
        bestMax = cm;
        best = { s: cut.s, c: cc, cs: ccs, d: dd, ds: dds };
      }
    }
    if (best && bestMax < 0.1 * ownMax) {
      C[a] = best.c;
      Cs[a] = best.cs;
      D[a] = best.d;
      Ds[a] = best.ds;
      cutUsed[a] = best.s;
    }
  });

  // Row scales of J and M separately (for row-scaled held-vector residual).
  const jRowMax = new Float64Array(N), mRowMax = new Float64Array(N);
  const { sizes, offA, offB, offC, offX } = Jg, nB = Jg.n;
  for (let blk = 0; blk < nB; blk++) {
    const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      let jm = 0, mm = 0;
      const idx = offX[blk] + r;
      for (let c = 0; c < mp; c++) {
        const jv = Math.abs(Jg.A[offA[blk] + r * mp + c]), mv = Math.abs(Mg.A[offA[blk] + r * mp + c]);
        if (jv > jm) jm = jv; if (mv > mm) mm = mv;
      }
      for (let c = 0; c < m; c++) {
        const jv = Math.abs(Jg.B[offB[blk] + r * m + c]), mv = Math.abs(Mg.B[offB[blk] + r * m + c]);
        if (jv > jm) jm = jv; if (mv > mm) mm = mv;
      }
      for (let c = 0; c < mn; c++) {
        const jv = Math.abs(Jg.C[offC[blk] + r * mn + c]), mv = Math.abs(Mg.C[offC[blk] + r * mn + c]);
        if (jv > jm) jm = jv; if (mv > mm) mm = mv;
      }
      for (let a = 0; a < K; a++) {
        const bv = Math.abs(B[a][idx]);
        if (bv > jm) jm = bv;
      }
      jRowMax[idx] = jm;
      mRowMax[idx] = mm;
    }
  }
  for (let a = 0; a < K; a++) {
    let jm = 0, mm = 0;
    for (let j = 0; j < Nc; j++) {
      const jv = Math.abs(C[a][j]), mv = Math.abs(Cs[a][j]);
      if (jv > jm) jm = jv; if (mv > mm) mm = mv;
    }
    for (let b = 0; b < K; b++) {
      const jv = Math.abs(D[a][b]), mv = Math.abs(Ds[a][b]);
      if (jv > jm) jm = jv; if (mv > mm) mm = mv;
    }
    jRowMax[Nc + a] = jm;
    mRowMax[Nc + a] = mm;
  }

  // Scratch vectors for matvecs.
  const tmp1 = new Float64Array(Nc), tmp2 = new Float64Array(Nc);
  const mulJDiff = (v, out) => {
    btMul(Jr, v, out);
    lin.apply(v, out, 1, 0);
    for (let a = 0; a < K; a++) {
      const va = v[Nc + a];
      if (va !== 0) for (let r = 0; r < Nc; r++) out[r] += B[a][r] * va;
    }
    for (let a = 0; a < K; a++) {
      let acc = 0;
      const Ca = C[a], Da = D[a];
      for (let j = 0; j < Nc; j++) acc += Ca[j] * v[j];
      for (let b = 0; b < K; b++) acc += Da[b] * v[Nc + b];
      out[Nc + a] = acc;
    }
    return out;
  };
  const mulMDiff = (v, out) => {
    btMul(Sr, v, out);
    lin.apply(v, out, 0, 1);
    for (let a = 0; a < K; a++) {
      let acc = 0;
      const Csa = Cs[a], Dsa = Ds[a];
      for (let j = 0; j < Nc; j++) acc += Csa[j] * v[j];
      for (let b = 0; b < K; b++) acc += Dsa[b] * v[Nc + b];
      out[Nc + a] = acc;
    }
    return out;
  };
  const mulM = (v, out) => {
    btMul(Mg, v, out);
    for (let a = 0; a < K; a++) {
      let acc = 0;
      const Csa = Cs[a], Dsa = Ds[a];
      for (let j = 0; j < Nc; j++) acc += Csa[j] * v[j];
      for (let b = 0; b < K; b++) acc += Dsa[b] * v[Nc + b];
      out[Nc + a] = acc;
    }
    return out;
  };
  const mulMT = (w, out) => {
    btMulT(Mg, w, out);
    for (let j = Nc; j < N; j++) out[j] = 0;
    for (let a = 0; a < K; a++) {
      const wa = w[Nc + a];
      if (wa === 0) continue;
      const Csa = Cs[a], Dsa = Ds[a];
      for (let j = 0; j < Nc; j++) out[j] += Csa[j] * wa;
      for (let b = 0; b < K; b++) out[Nc + b] += Dsa[b] * wa;
    }
    return out;
  };
  const absMulJM = (absV, lamAbs, out) => {
    out.fill(0);
    btAbsMul(Jg, absV, out, 1);
    if (lamAbs > 0) btAbsMul(Mg, absV, out, lamAbs);
    for (let a = 0; a < K; a++) {
      const va = absV[Nc + a];
      if (va !== 0) for (let r = 0; r < Nc; r++) out[r] += Math.abs(B[a][r]) * va;
    }
    for (let a = 0; a < K; a++) {
      let acc = 0;
      const Ca = C[a], Csa = Cs[a], Da = D[a], Dsa = Ds[a];
      for (let j = 0; j < Nc; j++) acc += (Math.abs(Ca[j]) + lamAbs * Math.abs(Csa[j])) * absV[j];
      for (let b = 0; b < K; b++) acc += (Math.abs(Da[b]) + lamAbs * Math.abs(Dsa[b])) * absV[Nc + b];
      out[Nc + a] = acc;
    }
    return out;
  };

  return {
    Jg, Mg, Jr, Sr, lin, Nc, K, N, fl, B, C, Cs, D, Ds,
    rawCs, rawDIs, rawDXs, lawRead, cutUsed, jRowMax, mRowMax,
    mulJDiff, mulMDiff, mulM, mulMT, absMulJM, tmp1, tmp2,
  };
}

// --- Deflation basis: exact conserved quantities, held-charge terminals, node-local invariants ---

function buildDeflation(s, pencil, coordScale) {
  const { Nc, K, N, fl, Cs, Ds, rawCs, rawDIs, rawDXs } = pencil;
  const { model, n, M, th, nSurf } = s;
  const vol = model.grid.vol;
  const invScale = Float64Array.from(coordScale, (v) => 1 / v);
  const rawVecs = [];
  const meta = [];

  // 1. Global conserved combinations from s.constraints (all feeding terminals driven by current).
  const handledCircuit = new Set();
  for (const cs of s.constraints) {
    if (cs.terminals.some((t) => s.terms[t].drive.kind !== 'I')) continue;
    const { dV } = lawAmount(s, cs);
    const vec = new Float64Array(N);
    for (let j = 0; j < cs.len; j++) vec[cs.idx[j]] += cs.w[j];
    fl.forEach((k, a) => { vec[Nc + a] += dV[k] ?? 0; });
    for (const { t, lambda } of cs.contacts ?? []) {
      const f = lambda / FARADAY;
      const rc = rawCs[t];
      for (let j = 0; j < Nc; j++) vec[j] += f * rc[j];
      fl.forEach((k, a) => {
        vec[Nc + a] += f * (k === t ? rawDIs[t] : rawDXs[t][k]);
      });
    }
    if (cs.circuit >= 0) handledCircuit.add(cs.circuit);
    const names = cs.parts.map(({ stretch }) => {
      const st = s.stretches[stretch], [r0, r1] = st.regions;
      return `${model.species[st.species].name} in ${model.regions[r0].name}${r1 > r0 ? `…${model.regions[r1].name}` : ''}`;
    });
    rawVecs.push({ kind: 'constraint', vec, label: `conserved: ${names.join(' + ') || cs.key}` });
  }

  // 2. Floating current-driven terminals that pass current only by charging (chargeOnly).
  fl.forEach((k, a) => {
    if (s.terms[k].drive.kind !== 'I' || !chargeOnly(s, k) || handledCircuit.has(k)) return;
    const vec = new Float64Array(N);
    vec.set(Cs[a], 0);
    for (let b = 0; b < K; b++) vec[Nc + b] = Ds[a][b];
    rawVecs.push({ kind: 'charge', vec, label: `held charge: terminal ${s.terms[k].name}` });
  });

  // 3. Node-local invariants (immobile species / trap pools), stored sparsely by block for O(N) projection.
  const localByBlock = new Map(); // blk → [{ vec: Float64Array(m), raw: ... }]
  const dA = s.dA;
  let localCount = 0;
  for (const lc of s.localConstraints) {
    for (const g of lc.nodes) {
      const blk = s.blockOfNode[g], off = pencil.Jg.offX[blk], m = pencil.Jg.sizes[blk];
      if (m === 0) continue;
      const locVec = new Float64Array(m);
      const vg = vol[g];
      for (const { i, q, k, w } of lc.parts) {
        if (i !== undefined) {
          s._dc(g, i, dA);
          for (let r = 0; r < M; r++) {
            const loc = s.loc[blk * M + r];
            if (loc >= 0 && dA[r] !== 0) locVec[loc] += w * vg * dA[r];
          }
        } else if (s.surfAt[g] === k) {
          const sf = model.surfaces[k], ns = sf.species.length;
          const tq = th[g * nSurf + q];
          const ws = w * vg * sf.species[0].capacity * s.surfArea[k][s.surfW[g]];
          for (let x = 0; x < ns; x++) {
            const loc = s.loc[blk * M + 1 + n + x];
            if (loc >= 0) locVec[loc] += ws * tq * ((q === x ? 1 : 0) - th[g * nSurf + x]);
          }
        }
      }
      // Scale to y-space (y = D x ⇒ cᵀ x = (D⁻¹ c)ᵀ y)
      let n0 = 0;
      for (let r = 0; r < m; r++) { locVec[r] *= invScale[off + r]; n0 += locVec[r] * locVec[r]; }
      n0 = Math.sqrt(n0);
      if (!(n0 > 0)) continue;
      const list = localByBlock.get(blk) ?? (localByBlock.set(blk, []), localByBlock.get(blk));
      // Twice-MGS against existing local vectors on this same block.
      for (let pass = 0; pass < 2; pass++) {
        for (const qv of list) {
          let dot = 0;
          for (let r = 0; r < m; r++) dot += qv[r] * locVec[r];
          for (let r = 0; r < m; r++) locVec[r] -= dot * qv[r];
        }
      }
      let n1 = 0;
      for (let r = 0; r < m; r++) n1 += locVec[r] * locVec[r];
      n1 = Math.sqrt(n1);
      if (n1 <= 1e-12 * n0) continue;
      for (let r = 0; r < m; r++) locVec[r] /= n1;
      list.push(locVec);
      localCount++;
    }
  }
  if (localCount > 0) meta.push({ kind: 'local', count: localCount, label: `${localCount} node-local invariant${localCount > 1 ? 's' : ''}` });

  // Orthonormalize global vectors in y-space (against local blocks first, then against each other).
  const globalQ = [];
  const debugVecs = [];
  for (const item of rawVecs) {
    const v = new Float64Array(N);
    let n0 = 0;
    for (let j = 0; j < N; j++) { v[j] = item.vec[j] * invScale[j]; n0 += v[j] * v[j]; }
    n0 = Math.sqrt(n0);
    if (!(n0 > 0)) continue;
    for (let pass = 0; pass < 2; pass++) {
      for (const [blk, list] of localByBlock) {
        const off = pencil.Jg.offX[blk], m = pencil.Jg.sizes[blk];
        for (const qv of list) {
          let dot = 0;
          for (let r = 0; r < m; r++) dot += qv[r] * v[off + r];
          for (let r = 0; r < m; r++) v[off + r] -= dot * qv[r];
        }
      }
      for (const qv of globalQ) {
        let dot = 0;
        for (let j = 0; j < N; j++) dot += qv[j] * v[j];
        for (let j = 0; j < N; j++) v[j] -= dot * qv[j];
      }
    }
    let n1 = 0;
    for (let j = 0; j < N; j++) n1 += v[j] * v[j];
    n1 = Math.sqrt(n1);
    if (n1 <= 1e-12 * n0) continue;
    for (let j = 0; j < N; j++) v[j] /= n1;
    globalQ.push(v);
    debugVecs.push({ label: item.label, vector: Array.from(item.vec) });
    meta.push({ kind: item.kind, label: item.label });
  }

  const totalDeflated = globalQ.length + localCount;
  const projectReal = (y) => {
    if (totalDeflated === 0) return y;
    for (let pass = 0; pass < 2; pass++) {
      for (const [blk, list] of localByBlock) {
        const off = pencil.Jg.offX[blk], m = pencil.Jg.sizes[blk];
        for (const qv of list) {
          let dot = 0;
          for (let r = 0; r < m; r++) dot += qv[r] * y[off + r];
          for (let r = 0; r < m; r++) y[off + r] -= dot * qv[r];
        }
      }
      for (const qv of globalQ) {
        let dot = 0;
        for (let j = 0; j < N; j++) dot += qv[j] * y[j];
        for (let j = 0; j < N; j++) y[j] -= dot * qv[j];
      }
    }
    return y;
  };
  const projectComplex = (yr, yi) => {
    if (totalDeflated === 0) return;
    projectReal(yr);
    projectReal(yi);
  };

  return { totalDeflated, localCount, meta, debugVecs, localByBlock, projectReal, projectComplex };
}

// --- Bordered shifted factorisations (real & complex) with forward + transposed solves ---

function factorRealShift(pencil, sigma) {
  const { Jg, Mg, Nc, K, N, B, C, Cs, D, Ds } = pencil;
  const { sizes, offA, offB, offC, offX } = Jg, nB = Jg.n;
  const bt = new BlockTridiagonal(nB, sizes);
  bt.bothWays = true;
  const sc = new Float64Array(Nc);
  for (let blk = 0; blk < nB; blk++) {
    const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      const idx = offX[blk] + r;
      let mx = 0;
      const oA = offA[blk] + r * mp, oB = offB[blk] + r * m, oC = offC[blk] + r * mn;
      for (let c = 0; c < mp; c++) { const v = Math.abs(Jg.A[oA + c] + sigma * Mg.A[oA + c]); if (v > mx) mx = v; }
      for (let c = 0; c < m; c++) { const v = Math.abs(Jg.B[oB + c] + sigma * Mg.B[oB + c]); if (v > mx) mx = v; }
      for (let c = 0; c < mn; c++) { const v = Math.abs(Jg.C[oC + c] + sigma * Mg.C[oC + c]); if (v > mx) mx = v; }
      for (let a = 0; a < K; a++) { const v = Math.abs(B[a][idx]); if (v > mx) mx = v; }
      const s = mx > 0 ? 1 / mx : 1;
      sc[idx] = s;
      for (let c = 0; c < mp; c++) bt.A[oA + c] = (Jg.A[oA + c] + sigma * Mg.A[oA + c]) * s;
      for (let c = 0; c < m; c++) bt.B[oB + c] = (Jg.B[oB + c] + sigma * Mg.B[oB + c]) * s;
      for (let c = 0; c < mn; c++) bt.C[oC + c] = (Jg.C[oC + c] + sigma * Mg.C[oC + c]) * s;
    }
  }
  try {
    bt.factor();
  } catch {
    return null;
  }
  const Csig = C.map((Ca, a) => {
    const out = new Float64Array(Nc), Csa = Cs[a];
    for (let j = 0; j < Nc; j++) out[j] = Ca[j] + sigma * Csa[j];
    return out;
  });
  const Dsig = D.map((Da, a) => {
    const out = new Float64Array(K), Dsa = Ds[a];
    for (let b = 0; b < K; b++) out[b] = Da[b] + sigma * Dsa[b];
    return out;
  });
  const scB = new Float64Array(K);
  for (let a = 0; a < K; a++) {
    let mx = 0;
    for (let j = 0; j < Nc; j++) { const v = Math.abs(Csig[a][j]); if (v > mx) mx = v; }
    for (let b = 0; b < K; b++) { const v = Math.abs(Dsig[a][b]); if (v > mx) mx = v; }
    scB[a] = mx > 0 ? 1 / mx : 1;
  }
  const X = [];
  let schurFac = null;
  if (K > 0) {
    const rhsTmp = new Float64Array(Nc);
    for (let a = 0; a < K; a++) {
      for (let j = 0; j < Nc; j++) rhsTmp[j] = B[a][j] * sc[j];
      X.push(bt.solve(rhsTmp));
    }
    const S = Array.from({ length: K }, () => new Float64Array(K));
    for (let a = 0; a < K; a++) {
      const Ca = Csig[a];
      for (let b = 0; b < K; b++) {
        let dot = 0;
        const Xb = X[b];
        for (let j = 0; j < Nc; j++) dot += Ca[j] * Xb[j];
        S[a][b] = (Dsig[a][b] - dot) * scB[a];
      }
    }
    schurFac = realLuFactor(S, K);
    if (!schurFac) return null;
  }
  const rhsG = new Float64Array(Nc), yG = new Float64Array(Nc);
  const solve = (rhs, out = new Float64Array(N)) => {
    for (let j = 0; j < Nc; j++) rhsG[j] = rhs[j] * sc[j];
    bt.solve(rhsG, yG);
    if (K === 0) {
      out.set(yG, 0);
      return out;
    }
    const rV = new Float64Array(K);
    for (let a = 0; a < K; a++) {
      let dot = 0;
      const Ca = Csig[a];
      for (let j = 0; j < Nc; j++) dot += Ca[j] * yG[j];
      rV[a] = (rhs[Nc + a] - dot) * scB[a];
    }
    const xV = realLuSolve(schurFac, rV);
    out.set(yG, 0);
    for (let a = 0; a < K; a++) {
      const xa = xV[a], Xa = X[a];
      out[Nc + a] = xa;
      if (xa !== 0) for (let j = 0; j < Nc; j++) out[j] -= Xa[j] * xa;
    }
    return out;
  };
  const solveT = (rhs, out = new Float64Array(N)) => {
    if (K === 0) {
      bt.solveT(rhs.subarray(0, Nc), yG);
      for (let j = 0; j < Nc; j++) out[j] = yG[j] * sc[j];
      return out;
    }
    const rV = new Float64Array(K);
    for (let a = 0; a < K; a++) {
      let dot = 0;
      const Xa = X[a];
      for (let j = 0; j < Nc; j++) dot += Xa[j] * rhs[j];
      rV[a] = rhs[Nc + a] - dot;
    }
    const wVt = realLuSolveT(schurFac, rV);
    for (let j = 0; j < Nc; j++) rhsG[j] = rhs[j];
    for (let a = 0; a < K; a++) {
      const wa = wVt[a] * scB[a];
      out[Nc + a] = wa;
      const Ca = Csig[a];
      if (wa !== 0) for (let j = 0; j < Nc; j++) rhsG[j] -= Ca[j] * wa;
    }
    bt.solveT(rhsG, yG);
    for (let j = 0; j < Nc; j++) out[j] = yG[j] * sc[j];
    return out;
  };
  return { solve, solveT, minPivotRatio: bt.minPivotRatio };
}

function factorComplexShift(pencil, sr, si) {
  const { Jg, Mg, Nc, K, N, B, C, Cs, D, Ds } = pencil;
  const { sizes, offA, offB, offC, offX } = Jg, nB = Jg.n;
  const cbt = new ComplexBlockTridiagonal(nB, sizes);
  const sc = new Float64Array(Nc);
  for (let blk = 0; blk < nB; blk++) {
    const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      const idx = offX[blk] + r;
      let mx = 0;
      const oA = offA[blk] + r * mp, oB = offB[blk] + r * m, oC = offC[blk] + r * mn;
      for (let c = 0; c < mp; c++) { const v = Math.hypot(Jg.A[oA + c] + sr * Mg.A[oA + c], si * Mg.A[oA + c]); if (v > mx) mx = v; }
      for (let c = 0; c < m; c++) { const v = Math.hypot(Jg.B[oB + c] + sr * Mg.B[oB + c], si * Mg.B[oB + c]); if (v > mx) mx = v; }
      for (let c = 0; c < mn; c++) { const v = Math.hypot(Jg.C[oC + c] + sr * Mg.C[oC + c], si * Mg.C[oC + c]); if (v > mx) mx = v; }
      for (let a = 0; a < K; a++) { const v = Math.abs(B[a][idx]); if (v > mx) mx = v; }
      const s = mx > 0 ? 1 / mx : 1;
      sc[idx] = s;
      for (let c = 0; c < mp; c++) { cbt.Ar[oA + c] = (Jg.A[oA + c] + sr * Mg.A[oA + c]) * s; cbt.Ai[oA + c] = si * Mg.A[oA + c] * s; }
      for (let c = 0; c < m; c++) { cbt.Br[oB + c] = (Jg.B[oB + c] + sr * Mg.B[oB + c]) * s; cbt.Bi[oB + c] = si * Mg.B[oB + c] * s; }
      for (let c = 0; c < mn; c++) { cbt.Cr[oC + c] = (Jg.C[oC + c] + sr * Mg.C[oC + c]) * s; cbt.Ci[oC + c] = si * Mg.C[oC + c] * s; }
    }
  }
  try {
    cbt.factor();
  } catch {
    return null;
  }
  const Cr = C.map((Ca, a) => { const o = new Float64Array(Nc); for (let j = 0; j < Nc; j++) o[j] = Ca[j] + sr * Cs[a][j]; return o; });
  const Ci = Cs.map((Csa) => { const o = new Float64Array(Nc); for (let j = 0; j < Nc; j++) o[j] = si * Csa[j]; return o; });
  const Dr = D.map((Da, a) => Float64Array.from(Da, (v, b) => v + sr * Ds[a][b]));
  const Di = Ds.map((Dsa) => Float64Array.from(Dsa, (v) => si * v));
  const scB = new Float64Array(K);
  for (let a = 0; a < K; a++) {
    let mx = 0;
    for (let j = 0; j < Nc; j++) { const v = Math.hypot(Cr[a][j], Ci[a][j]); if (v > mx) mx = v; }
    for (let b = 0; b < K; b++) { const v = Math.hypot(Dr[a][b], Di[a][b]); if (v > mx) mx = v; }
    scB[a] = mx > 0 ? 1 / mx : 1;
  }
  const Xr = [], Xi = [];
  let schurFac = null;
  if (K > 0) {
    const br = new Float64Array(Nc), bi = new Float64Array(Nc);
    for (let a = 0; a < K; a++) {
      for (let j = 0; j < Nc; j++) br[j] = B[a][j] * sc[j];
      const xr = new Float64Array(Nc), xi = new Float64Array(Nc);
      cbt.solve(br, bi, xr, xi);
      Xr.push(xr);
      Xi.push(xi);
    }
    const Sr = Array.from({ length: K }, () => new Float64Array(K));
    const Si = Array.from({ length: K }, () => new Float64Array(K));
    for (let a = 0; a < K; a++) {
      const cra = Cr[a], cia = Ci[a];
      for (let b = 0; b < K; b++) {
        let dr = 0, di = 0;
        const xrb = Xr[b], xib = Xi[b];
        for (let j = 0; j < Nc; j++) {
          dr += cra[j] * xrb[j] - cia[j] * xib[j];
          di += cra[j] * xib[j] + cia[j] * xrb[j];
        }
        Sr[a][b] = (Dr[a][b] - dr) * scB[a];
        Si[a][b] = (Di[a][b] - di) * scB[a];
      }
    }
    schurFac = complexLuFactor(Sr, Si, K);
    if (!schurFac) return null;
  }
  const rhsGr = new Float64Array(Nc), rhsGi = new Float64Array(Nc);
  const yGr = new Float64Array(Nc), yGi = new Float64Array(Nc);
  const solve = (br, bi, outR = new Float64Array(N), outI = new Float64Array(N)) => {
    for (let j = 0; j < Nc; j++) { rhsGr[j] = br[j] * sc[j]; rhsGi[j] = bi[j] * sc[j]; }
    cbt.solve(rhsGr, rhsGi, yGr, yGi);
    if (K === 0) { outR.set(yGr, 0); outI.set(yGi, 0); return [outR, outI]; }
    const rVr = new Float64Array(K), rVi = new Float64Array(K);
    for (let a = 0; a < K; a++) {
      let dr = 0, di = 0;
      const cra = Cr[a], cia = Ci[a];
      for (let j = 0; j < Nc; j++) {
        dr += cra[j] * yGr[j] - cia[j] * yGi[j];
        di += cra[j] * yGi[j] + cia[j] * yGr[j];
      }
      rVr[a] = (br[Nc + a] - dr) * scB[a];
      rVi[a] = (bi[Nc + a] - di) * scB[a];
    }
    const [xVr, xVi] = complexLuSolve(schurFac, rVr, rVi);
    outR.set(yGr, 0); outI.set(yGi, 0);
    for (let a = 0; a < K; a++) {
      const ar = xVr[a], ai = xVi[a], xra = Xr[a], xia = Xi[a];
      outR[Nc + a] = ar; outI[Nc + a] = ai;
      for (let j = 0; j < Nc; j++) {
        outR[j] -= xra[j] * ar - xia[j] * ai;
        outI[j] -= xra[j] * ai + xia[j] * ar;
      }
    }
    return [outR, outI];
  };
  const solveH = (br, bi, outR = new Float64Array(N), outI = new Float64Array(N)) => {
    if (K === 0) {
      cbt.solveH(br.subarray(0, Nc), bi.subarray(0, Nc), yGr, yGi);
      for (let j = 0; j < Nc; j++) { outR[j] = yGr[j] * sc[j]; outI[j] = yGi[j] * sc[j]; }
      return [outR, outI];
    }
    // (Pᴴ) [w_g; w_V] = [r_g; r_V] with P_gg = diag(sc)⁻¹ T, S = diag(scB)(D - C X):
    // Sᴴ (w_V / scB) = r_V − Xᴴ r_g, then Tᴴ (w_g / sc) = r_g − Cᴴ w_V.
    const rVr = new Float64Array(K), rVi = new Float64Array(K);
    for (let a = 0; a < K; a++) {
      let dr = 0, di = 0;
      const xra = Xr[a], xia = Xi[a];
      for (let j = 0; j < Nc; j++) {
        dr += xra[j] * br[j] + xia[j] * bi[j];
        di += xra[j] * bi[j] - xia[j] * br[j];
      }
      rVr[a] = br[Nc + a] - dr;
      rVi[a] = bi[Nc + a] - di;
    }
    const [wVtr, wVti] = complexLuSolveH(schurFac, rVr, rVi);
    for (let j = 0; j < Nc; j++) { rhsGr[j] = br[j]; rhsGi[j] = bi[j]; }
    for (let a = 0; a < K; a++) {
      const war = wVtr[a] * scB[a], wai = wVti[a] * scB[a];
      outR[Nc + a] = war; outI[Nc + a] = wai;
      const cra = Cr[a], cia = Ci[a];
      for (let j = 0; j < Nc; j++) {
        // Cᴴ w_V: conj(C[a,j]) * w_V[a]
        rhsGr[j] -= cra[j] * war + cia[j] * wai;
        rhsGi[j] -= cra[j] * wai - cia[j] * war;
      }
    }
    cbt.solveH(rhsGr, rhsGi, yGr, yGi);
    for (let j = 0; j < Nc; j++) { outR[j] = yGr[j] * sc[j]; outI[j] = yGi[j] * sc[j]; }
    return [outR, outI];
  };
  return { solve, solveH };
}

// --- Deterministic LCG vector for reproducible Arnoldi starts ---

function seededVector(N, seed = 0x51a7b11) {
  const v = new Float64Array(N);
  let state = seed >>> 0;
  for (let j = 0; j < N; j++) {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
    v[j] = state / 2147483648 - 1;
  }
  return v;
}

// --- Shift-invert Arnoldi (real & complex) with twice-MGS reorthogonalization ---

function arnoldiReal(op, N, mMax, startVec, project) {
  const V = [];
  const v0 = Float64Array.from(startVec);
  project(v0);
  // Two purification passes through T_σ suppress index-1/2 DAE infinity directions before Arnoldi.
  for (let p = 0; p < 2; p++) {
    const w = op(v0);
    project(w);
    let nrm = 0;
    for (let j = 0; j < N; j++) nrm += w[j] * w[j];
    nrm = Math.sqrt(nrm);
    if (!(nrm > 0) || !Number.isFinite(nrm)) return [];
    for (let j = 0; j < N; j++) v0[j] = w[j] / nrm;
  }
  V.push(v0);
  const m = Math.min(mMax, N);
  const H = Array.from({ length: m + 1 }, () => new Float64Array(m));
  let k = 0;
  for (; k < m; k++) {
    const w = op(V[k]);
    project(w);
    let normW0 = 0;
    for (let j = 0; j < N; j++) normW0 += w[j] * w[j];
    normW0 = Math.sqrt(normW0);
    if (!Number.isFinite(normW0)) break;
    // Twice-MGS reorthogonalization
    for (let pass = 0; pass < 2; pass++) {
      for (let i = 0; i <= k; i++) {
        const Vi = V[i];
        let dot = 0;
        for (let j = 0; j < N; j++) dot += Vi[j] * w[j];
        H[i][k] += dot;
        for (let j = 0; j < N; j++) w[j] -= dot * Vi[j];
      }
      project(w);
    }
    let hn = 0;
    for (let j = 0; j < N; j++) hn += w[j] * w[j];
    hn = Math.sqrt(hn);
    H[k + 1][k] = hn;
    if (!(hn > 1e-14 * Math.max(1, normW0))) {
      k++;
      break;
    }
    const vNext = new Float64Array(N);
    for (let j = 0; j < N; j++) vNext[j] = w[j] / hn;
    V.push(vNext);
  }
  if (k === 0) return [];
  const eigs = realHessenbergEig(H, k);
  const out = [];
  for (const e of eigs) {
    const vr = new Float64Array(N), vi = new Float64Array(N);
    for (let i = 0; i < k; i++) {
      const Vi = V[i], zr = e.zr[i], zi = e.zi[i];
      for (let j = 0; j < N; j++) {
        vr[j] += zr * Vi[j];
        vi[j] += zi * Vi[j];
      }
    }
    const ritzEst = Math.abs(H[k][k - 1]) * Math.hypot(e.zr[k - 1], e.zi[k - 1]);
    out.push({ muRe: e.re, muIm: e.im, vr, vi, ritzEst });
  }
  return out;
}

function arnoldiComplex(op, N, mMax, startVec, projectComplex) {
  const Vr = [], Vi = [];
  let v0r = Float64Array.from(startVec), v0i = new Float64Array(N);
  projectComplex(v0r, v0i);
  for (let p = 0; p < 2; p++) {
    const [wr, wi] = op(v0r, v0i);
    projectComplex(wr, wi);
    let n2 = 0;
    for (let j = 0; j < N; j++) n2 += wr[j] * wr[j] + wi[j] * wi[j];
    const nrm = Math.sqrt(n2);
    if (!(nrm > 0) || !Number.isFinite(nrm)) return [];
    for (let j = 0; j < N; j++) { v0r[j] = wr[j] / nrm; v0i[j] = wi[j] / nrm; }
  }
  Vr.push(v0r); Vi.push(v0i);
  const m = Math.min(mMax, N);
  const Hr = Array.from({ length: m + 1 }, () => new Float64Array(m));
  const Hi = Array.from({ length: m + 1 }, () => new Float64Array(m));
  let k = 0;
  for (; k < m; k++) {
    const [wr, wi] = op(Vr[k], Vi[k]);
    projectComplex(wr, wi);
    let n0 = 0;
    for (let j = 0; j < N; j++) n0 += wr[j] * wr[j] + wi[j] * wi[j];
    n0 = Math.sqrt(n0);
    if (!Number.isFinite(n0)) break;
    for (let pass = 0; pass < 2; pass++) {
      for (let i = 0; i <= k; i++) {
        const vr = Vr[i], vi = Vi[i];
        let dr = 0, di = 0;
        for (let j = 0; j < N; j++) {
          dr += vr[j] * wr[j] + vi[j] * wi[j];
          di += vr[j] * wi[j] - vi[j] * wr[j];
        }
        Hr[i][k] += dr; Hi[i][k] += di;
        for (let j = 0; j < N; j++) {
          wr[j] -= dr * vr[j] - di * vi[j];
          wi[j] -= dr * vi[j] + di * vr[j];
        }
      }
      projectComplex(wr, wi);
    }
    let n2 = 0;
    for (let j = 0; j < N; j++) n2 += wr[j] * wr[j] + wi[j] * wi[j];
    const hn = Math.sqrt(n2);
    Hr[k + 1][k] = hn;
    if (!(hn > 1e-14 * Math.max(1, n0))) { k++; break; }
    const nr = new Float64Array(N), ni = new Float64Array(N);
    for (let j = 0; j < N; j++) { nr[j] = wr[j] / hn; ni[j] = wi[j] / hn; }
    Vr.push(nr); Vi.push(ni);
  }
  if (k === 0) return [];
  const eigs = complexHessenbergEig(Hr, Hi, k);
  const out = [];
  for (const e of eigs) {
    const vr = new Float64Array(N), vi = new Float64Array(N);
    for (let i = 0; i < k; i++) {
      const vri = Vr[i], vii = Vi[i], zr = e.zr[i], zi = e.zi[i];
      for (let j = 0; j < N; j++) {
        vr[j] += zr * vri[j] - zi * vii[j];
        vi[j] += zr * vii[j] + zi * vri[j];
      }
    }
    const ritzEst = Math.hypot(Hr[k][k - 1], Hi[k][k - 1]) * Math.hypot(e.zr[k - 1], e.zi[k - 1]);
    out.push({ muRe: e.re, muIm: e.im, vr, vi, ritzEst });
  }
  return out;
}

// --- Refinement, left eigenvector, held-vector difference-form residual, and error bar ---

function normalizeComplex(vr, vi) {
  let n2 = 0;
  for (let j = 0; j < vr.length; j++) n2 += vr[j] * vr[j] + vi[j] * vi[j];
  const nrm = Math.sqrt(n2);
  if (!(nrm > 0) || !Number.isFinite(nrm)) return 0;
  for (let j = 0; j < vr.length; j++) { vr[j] /= nrm; vi[j] /= nrm; }
  return nrm;
}

function evalResidual(pencil, vr, vi, lamRe, lamIm) {
  const { N, jRowMax, mRowMax, mulJDiff, mulMDiff } = pencil;
  const jvr = new Float64Array(N), jvi = new Float64Array(N);
  const mvr = new Float64Array(N), mvi = new Float64Array(N);
  mulJDiff(vr, jvr);
  mulMDiff(vr, mvr);
  const isComplex = lamIm !== 0 || vi.some((x) => x !== 0);
  if (isComplex) {
    mulJDiff(vi, jvi);
    mulMDiff(vi, mvi);
  }
  const lamAbs = Math.hypot(lamRe, lamIm);
  let vInf = 0;
  for (let j = 0; j < N; j++) { const a = Math.hypot(vr[j], vi[j]); if (a > vInf) vInf = a; }
  if (!(vInf > 0)) return { residual: Infinity, jvr, jvi, mvr, mvi };
  let resInf = 0;
  for (let i = 0; i < N; i++) {
    const rr = jvr[i] + lamRe * mvr[i] - lamIm * mvi[i];
    const ri = jvi[i] + lamRe * mvi[i] + lamIm * mvr[i];
    const sc = jRowMax[i] + lamAbs * mRowMax[i];
    if (sc > 0) {
      const r = Math.hypot(rr, ri) / (sc * vInf);
      if (r > resInf) resInf = r;
    }
  }
  return { residual: resInf, jvr, jvi, mvr, mvi };
}

function vecOverlap(a, b) {
  let dr = 0, di = 0;
  for (let j = 0; j < a.vr.length; j++) {
    dr += a.vr[j] * b.vr[j] + a.vi[j] * b.vi[j];
    di += a.vr[j] * b.vi[j] - a.vi[j] * b.vr[j];
  }
  return Math.hypot(dr, di);
}

function refineAndBound(pencil, deflation, coordScale, cand) {
  const { N, mulM, mulMT, absMulJM } = pencil;
  const invScale = Float64Array.from(coordScale, (v) => 1 / v);
  // Convert Ritz vector from y-space (y = D x) to x-space (x = D⁻¹ y).
  let vr = new Float64Array(N), vi = new Float64Array(N);
  for (let j = 0; j < N; j++) { vr[j] = cand.vr[j] * invScale[j]; vi[j] = cand.vi[j] * invScale[j]; }
  let lamRe = cand.re, lamIm = cand.im;
  const isReal = Math.abs(lamIm) <= 1e-10 * Math.max(1, Math.abs(lamRe));
  if (isReal) { lamIm = 0; vi.fill(0); }
  if (!normalizeComplex(vr, vi)) return null;

  const wr = new Float64Array(N), wi = new Float64Array(N);

  const twoSidedPass = (sRe, sIm) => {
    const bump = Math.max(1e-7, 1e-7 * Math.hypot(sRe, sIm));
    const dir = sRe >= 0 ? 1 : -1;
    let haveLeft = false;
    if (isReal) {
      const fac = factorRealShift(pencil, sRe + dir * bump) || factorRealShift(pencil, sRe - dir * bump);
      if (!fac) return null;
      const mv = new Float64Array(N), y = new Float64Array(N);
      for (let it = 0; it < 2; it++) {
        mulM(vr, mv);
        fac.solve(mv, vr);
        for (let j = 0; j < N; j++) y[j] = vr[j] * coordScale[j];
        deflation.projectReal(y);
        for (let j = 0; j < N; j++) vr[j] = y[j] * invScale[j];
        if (!normalizeComplex(vr, vi)) return null;
      }
      mulM(vr, wr);
      if (!normalizeComplex(wr, wi)) { wr.set(seededVector(N, 0x1e57)); normalizeComplex(wr, wi); }
      const mtw = new Float64Array(N);
      for (let it = 0; it < 3; it++) {
        mulMT(wr, mtw);
        fac.solveT(mtw, wr);
        if (!normalizeComplex(wr, wi)) break;
      }
      haveLeft = normalizeComplex(wr, wi) > 0;
    } else {
      const fac = factorComplexShift(pencil, sRe + dir * bump, sIm) || factorComplexShift(pencil, sRe - dir * bump, sIm);
      if (!fac) return null;
      const mvr = new Float64Array(N), mvi = new Float64Array(N);
      const yr = new Float64Array(N), yi = new Float64Array(N);
      for (let it = 0; it < 2; it++) {
        mulM(vr, mvr); mulM(vi, mvi);
        fac.solve(mvr, mvi, vr, vi);
        for (let j = 0; j < N; j++) { yr[j] = vr[j] * coordScale[j]; yi[j] = vi[j] * coordScale[j]; }
        deflation.projectComplex(yr, yi);
        for (let j = 0; j < N; j++) { vr[j] = yr[j] * invScale[j]; vi[j] = yi[j] * invScale[j]; }
        if (!normalizeComplex(vr, vi)) return null;
      }
      mulM(vr, wr); mulM(vi, wi);
      if (!normalizeComplex(wr, wi)) { wr.set(seededVector(N, 0x1e57)); wi.fill(0); normalizeComplex(wr, wi); }
      const mtwr = new Float64Array(N), mtwi = new Float64Array(N);
      for (let it = 0; it < 3; it++) {
        mulMT(wr, mtwr); mulMT(wi, mtwi);
        fac.solveH(mtwr, mtwi, wr, wi);
        if (!normalizeComplex(wr, wi)) break;
      }
      haveLeft = normalizeComplex(wr, wi) > 0;
    }
    if (!haveLeft) return null;
    let check = evalResidual(pencil, vr, vi, sRe, sIm);
    let numR = 0, numI = 0, denR = 0, denI = 0;
    for (let j = 0; j < N; j++) {
      numR -= wr[j] * check.jvr[j] + wi[j] * check.jvi[j];
      numI -= wr[j] * check.jvi[j] - wi[j] * check.jvr[j];
      denR += wr[j] * check.mvr[j] + wi[j] * check.mvi[j];
      denI += wr[j] * check.mvi[j] - wi[j] * check.mvr[j];
    }
    const denAbs = Math.hypot(denR, denI);
    if (!(denAbs > 0) || !Number.isFinite(denAbs)) return null;
    const [rqRe, rqImRaw] = cdiv(numR, numI, denR, denI);
    const rqIm = isReal ? 0 : rqImRaw;
    let outRe = sRe, outIm = sIm;
    if (Number.isFinite(rqRe) && Number.isFinite(rqIm)) {
      const rqCheck = evalResidual(pencil, vr, vi, rqRe, rqIm);
      if (rqCheck.residual <= check.residual * 1.1 + 1e-15) {
        outRe = rqRe;
        outIm = rqIm;
        check = rqCheck;
      }
    }
    return { re: outRe, im: outIm, check, denAbs };
  };

  let best = twoSidedPass(lamRe, lamIm);
  if (!best) return null;
  // Re-center once if Rayleigh quotient moved noticeably or residual is not yet tight.
  const moved = Math.hypot(best.re - lamRe, best.im - lamIm) > 1e-5 * Math.max(1, Math.hypot(best.re, best.im));
  if (moved || best.check.residual > 1e-9) {
    const vrKeep = vr.slice(), viKeep = vi.slice(), wrKeep = wr.slice(), wiKeep = wi.slice();
    const second = twoSidedPass(best.re, best.im);
    if (second && second.check.residual <= best.check.residual) {
      best = second;
    } else {
      vr.set(vrKeep); vi.set(viKeep); wr.set(wrKeep); wi.set(wiKeep);
    }
  }
  lamRe = best.re;
  lamIm = best.im;
  const { check, denAbs } = best;

  // First-order componentwise error bar: ε · |w|ᴴ (|J| + |λ| |M|) |v| / |wᴴ M v|.
  const absV = new Float64Array(N), absJMv = new Float64Array(N);
  for (let j = 0; j < N; j++) absV[j] = Math.hypot(vr[j], vi[j]);
  absMulJM(absV, Math.hypot(lamRe, lamIm), absJMv);
  let numBar = 0;
  for (let j = 0; j < N; j++) numBar += Math.hypot(wr[j], wi[j]) * absJMv[j];
  const error = (EPS * numBar) / denAbs;
  if (!Number.isFinite(error) || !Number.isFinite(check.residual)) return null;

  return { re: lamRe, im: lamIm, residual: check.residual, error, vr, vi };
}

// --- Deterministic-phase profile extraction ---

function modeProfiles(s, pencil, vr, vi) {
  const { M, nB, VT } = s, { Nc, K, N, fl } = pencil;
  // Rotate complex eigenvector so its largest component is real and positive.
  let idxMax = 0, best = 0;
  for (let j = 0; j < N; j++) { const a = Math.hypot(vr[j], vi[j]); if (a > best) { best = a; idxMax = j; } }
  const cr = best > 0 ? vr[idxMax] / best : 1, ci = best > 0 ? -vi[idxMax] / best : 0;
  const scale = best > 0 ? 1 / best : 1;
  const xr = new Float64Array(nB * M), xi = new Float64Array(nB * M);
  for (let o = 0; o < nB * M; o++) {
    const r = s.rix[o];
    if (r >= Nc) continue;
    xr[o] = (vr[r] * cr - vi[r] * ci) * scale;
    xi[o] = (vr[r] * ci + vi[r] * cr) * scale;
  }
  const prof = s._smallSignalProfiles(xr, xi);
  const termV = {};
  fl.forEach((k, a) => {
    const re = (vr[Nc + a] * cr - vi[Nc + a] * ci) * scale;
    const im = (vr[Nc + a] * ci + vi[Nc + a] * cr) * scale;
    termV[s.terms[k].name] = { re, im };
  });
  return { ...prof, termV, scaleV: VT };
}

// --- Default real shift grid from the device's own timescales (never σ = 0) ---

function defaultShifts(s) {
  const tauSlow = s.slowestTime(), tauFast = s.fastestTime();
  const rSlow = Number.isFinite(tauSlow) && tauSlow > 0 ? 1 / tauSlow : 1;
  const rFast = Number.isFinite(tauFast) && tauFast > 0 ? 1 / tauFast : 1e6;
  const lo = Math.max(1e-4, Math.min(1, rSlow * 0.1));
  const hi = Math.min(1e8, Math.max(1e3, rFast * 10));
  const count = 6;
  const shifts = [];
  const logLo = Math.log(lo), logHi = Math.log(hi);
  for (let i = 0; i < count; i++) {
    shifts.push(Math.exp(logLo + ((logHi - logLo) * i) / (count - 1)));
  }
  return shifts;
}

// --- Un-kicked dynamic steadiness check (Fable §5, Sol §6) ---

function checkDynamicSteady(s) {
  const tau = s.slowestTime();
  const baseH = Math.min(1e-1, Math.max(1e-9, Number.isFinite(tau) && tau > 0 ? 1e-2 * tau : 1e-3));
  const u0 = s.u.slice(), uLo0 = s.uLo.slice(), v0 = s.termV.slice();
  const snap = s._snapshot();
  let lastDrift = Infinity, stepOk = false;
  try {
    for (const h of [baseH, baseH * 0.01]) {
      s._restore(snap);
      const r = s.step(h);
      if (!r.converged) continue;
      let drift = 0;
      for (let j = 0; j < s.u.length; j++) {
        const d = Math.abs((s.u[j] - u0[j]) + (s.uLo[j] - uLo0[j]));
        if (d > drift) drift = d;
      }
      for (let k = 0; k < s.termV.length; k++) {
        const d = Math.abs(s.termV[k] - v0[k]) / s.VT;
        if (d > drift) drift = d;
      }
      lastDrift = drift;
      stepOk = true;
      break;
    }
  } finally {
    s._restore(snap);
  }
  return { ok: stepOk && lastDrift < 1e-5, converged: stepOk, drift: lastDrift, h: baseH };
}

export class Stability {
  /**
   * Linear stability of the current converged steady state on the dynamic DAE pencil −J v = λ M v.
   * @param {{ shifts?: number[], frequencies?: number[], modes?: number, profiles?: boolean, kBar?: number, debug?: boolean }} [opts]
   */
  stability(opts = {}) {
    if (!this.atSteady || this.regriddedSteady) {
      throw new SolverError('stability: solve() first — stability() analyses the current converged steady state (preserving its branch)');
    }
    const { shifts: userShifts, frequencies = [], modes: nWant = 8, profiles = false, kBar = 10, debug = false } = opts;
    if (userShifts !== undefined) {
      for (const sv of userShifts) {
        if (!(sv > 0 && Number.isFinite(sv))) throw new SolverError(`stability: real shifts must be positive and finite (never 0); got ${sv}`);
      }
    }
    for (const f of frequencies) {
      if (!(f > 0 && Number.isFinite(f))) throw new SolverError(`stability: frequencies must be positive and finite (Hz); got ${f}`);
    }
    if (!(Number.isInteger(nWant) && nWant >= 1)) throw new SolverError(`stability: modes must be a positive integer; got ${nWant}`);

    const snap = this._snapshot();
    const wasAll = this.allTerminals, wasSteady = this.steady, wasLawForm = this.lawForm;
    const learnedKeep = { lostMode: this.lostMode, accepted: this.accepted, firstErrs: this.firstErrs?.slice() };
    const warnings = [];
    if (this.unkept.length > 0) {
      warnings.push(`stability: ${this.unkept.length} conserved combination(s) are in unkept; neutral directions may remain.`);
    }
    const nonAuto = this.terms.some((t) => Boolean(t.drive.src?.t));
    if (nonAuto) {
      warnings.push(`stability: non-autonomous waveform drive frozen at t = ${this.time} s; instantaneous spectrum does not establish trajectory stability.`);
    }

    try {
      // 1. Precondition: verify state is actually steady for the time-dependent equations.
      const steadyCheck = checkDynamicSteady(this);

      // 2. Assemble dynamic bordered pencil (steady = false) with Sol's cut-current readings.
      //    (Also supplies the dynamic steady residual when an extreme C_face/G makes the un-cut
      //    time stepper's Newton solve fail to converge despite zero dynamic residual.)
      this.allTerminals = true;
      this.steady = false;
      const pencil = buildPencil(this);
      const { Nc, K, N } = pencil;

      let dynamicResInf = 0;
      for (let i = 0; i < Nc; i++) {
        const sc = pencil.jRowMax[i];
        const r = sc > 0 ? Math.abs(this.res[i]) / sc : Math.abs(this.res[i]);
        if (r > dynamicResInf) dynamicResInf = r;
      }
      for (let a = 0; a < K; a++) {
        const k = pencil.fl[a];
        let mx = 0;
        const rc = this.termC[k];
        for (let j = 0; j < Nc; j++) { const v = Math.abs(rc[j]); if (v > mx) mx = v; }
        mx = Math.max(mx, Math.abs(this.termDI[k]));
        for (const k2 of pencil.fl) if (k2 !== k) mx = Math.max(mx, Math.abs(this.termDX[k][k2]));
        const r = mx > 0 ? Math.abs(this.termRes[k]) / mx : Math.abs(this.termRes[k]);
        if (r > dynamicResInf) dynamicResInf = r;
      }

      let steadyOk = steadyCheck.ok;
      let reportedDrift = steadyCheck.drift;
      if (!steadyCheck.converged && dynamicResInf < 1e-10) {
        steadyOk = true;
        reportedDrift = 0;
        warnings.push(`stability: un-kicked time step did not converge, but row-scaled dynamic residual is ${dynamicResInf.toExponential(2)}; proceeding on cut-current pencil.`);
      }
      if (!steadyOk) {
        const reason = !steadyCheck.converged
          ? `not a steady state of the time-dependent equations (an un-kicked time step did not converge; dynamic residual ${dynamicResInf.toExponential(2)})`
          : `not a steady state of the time-dependent equations (an un-kicked step of ${steadyCheck.h.toPrecision(2)} s moved the state by ${steadyCheck.drift.toPrecision(3)} thermal units)`;
        return {
          verdict: 'undecided',
          reason,
          drift: steadyCheck.drift,
          modes: [],
          deflated: [],
          shifts: { attempted: [], succeeded: [] },
          window: null,
          warnings,
        };
      }

      // Coordinate weights D for Arnoldi metric: thermal units everywhere (terminal V scaled by 1/VT).
      const coordScale = new Float64Array(N).fill(1);
      for (let a = 0; a < K; a++) coordScale[Nc + a] = 1 / this.VT;
      const invScale = Float64Array.from(coordScale, (v) => 1 / v);

      // 3. Exact deflation basis for conserved amounts, held-charge terminals, and node-local invariants.
      const deflation = buildDeflation(this, pencil, coordScale);
      const startVec = seededVector(N);

      const debugData = debug
        ? {
            J: (() => {
              const Mj = denseOfBt(pencil.Jg, K);
              for (let r = 0; r < Nc; r++) for (let a = 0; a < K; a++) Mj[r][Nc + a] = pencil.B[a][r];
              for (let a = 0; a < K; a++) {
                for (let c = 0; c < Nc; c++) Mj[Nc + a][c] = pencil.C[a][c];
                for (let b = 0; b < K; b++) Mj[Nc + a][Nc + b] = pencil.D[a][b];
              }
              return Mj;
            })(),
            M: (() => {
              const Mm = denseOfBt(pencil.Mg, K);
              for (let a = 0; a < K; a++) {
                for (let c = 0; c < Nc; c++) Mm[Nc + a][c] = pencil.Cs[a][c];
                for (let b = 0; b < K; b++) Mm[Nc + a][Nc + b] = pencil.Ds[a][b];
              }
              return Mm;
            })(),
            deflation: deflation.debugVecs,
            localDeflated: deflation.localCount,
            cutUsed: pencil.cutUsed,
            labels: (() => {
              const lbl = new Array(N);
              for (let o = 0; o < this.nB * this.M; o++) if (this.rix[o] < Nc) lbl[this.rix[o]] = [Math.floor(o / this.M), o % this.M];
              pencil.fl.forEach((k, a) => (lbl[Nc + a] = ['V', this.terms[k].name]));
              return lbl;
            })(),
          }
        : undefined;

      // 4. Shift-invert Arnoldi over real positive shifts (+ optional complex frequency shifts).
      const realShifts = userShifts ? Array.from(userShifts) : defaultShifts(this);
      const attempted = [];
      const succeeded = [];
      const candidates = [];
      const mKrylov = Math.min(32, Math.max(4, N - deflation.totalDeflated));
      const rCeiling = Math.max(1e9, ...realShifts.map((x) => x * 1e3));

      const xTmp = new Float64Array(N), mxTmp = new Float64Array(N), solTmp = new Float64Array(N);
      for (const sigma of realShifts) {
        attempted.push({ re: sigma, im: 0 });
        const fac = factorRealShift(pencil, sigma);
        if (!fac) continue;
        succeeded.push({ re: sigma, im: 0 });
        const op = (y) => {
          for (let j = 0; j < N; j++) xTmp[j] = y[j] * invScale[j];
          pencil.mulM(xTmp, mxTmp);
          fac.solve(mxTmp, solTmp);
          const out = new Float64Array(N);
          for (let j = 0; j < N; j++) out[j] = solTmp[j] * coordScale[j];
          return out;
        };
        const ritz = arnoldiReal(op, N, mKrylov, startVec, deflation.projectReal);
        for (const r of ritz) {
          const muMag = Math.hypot(r.muRe, r.muIm);
          if (!(muMag > 1e-11 / Math.max(1, sigma))) continue;
          const [invRe, invIm] = cdiv(1, 0, r.muRe, r.muIm);
          const re = sigma - invRe, im = -invIm;
          const lamMag = Math.hypot(re, im);
          if (!Number.isFinite(re) || !Number.isFinite(im) || lamMag > rCeiling) continue;
          // First-order Arnoldi error estimate in λ-space (used to rank duplicate candidates, not to drop non-normal modes).
          const deltaLam = r.ritzEst / (muMag * muMag);
          // Skip lower-half duplicate of conjugate pair; we'll emit both members cleanly after refinement.
          if (im < -1e-10 * Math.max(1, Math.abs(re))) continue;
          candidates.push({ re, im, vr: r.vr, vi: r.vi, ritzEst: r.ritzEst, deltaLam, shift: sigma });
        }
      }

      const xTmpR = new Float64Array(N), xTmpI = new Float64Array(N);
      const mxTmpR = new Float64Array(N), mxTmpI = new Float64Array(N);
      const solTmpR = new Float64Array(N), solTmpI = new Float64Array(N);
      let maxOmega = 0;
      for (const f of frequencies) {
        const w = 2 * Math.PI * f;
        const sr = Math.max(1e-2, 1e-3 * w);
        attempted.push({ re: sr, im: w });
        const fac = factorComplexShift(pencil, sr, w);
        if (!fac) continue;
        succeeded.push({ re: sr, im: w });
        if (w > maxOmega) maxOmega = w;
        const op = (yr, yi) => {
          for (let j = 0; j < N; j++) { xTmpR[j] = yr[j] * invScale[j]; xTmpI[j] = yi[j] * invScale[j]; }
          pencil.mulM(xTmpR, mxTmpR);
          pencil.mulM(xTmpI, mxTmpI);
          fac.solve(mxTmpR, mxTmpI, solTmpR, solTmpI);
          const outR = new Float64Array(N), outI = new Float64Array(N);
          for (let j = 0; j < N; j++) { outR[j] = solTmpR[j] * coordScale[j]; outI[j] = solTmpI[j] * coordScale[j]; }
          return [outR, outI];
        };
        const ritz = arnoldiComplex(op, N, mKrylov, startVec, deflation.projectComplex);
        for (const r of ritz) {
          const muMag = Math.hypot(r.muRe, r.muIm);
          if (!(muMag > 1e-11 / Math.max(1, Math.hypot(sr, w)))) continue;
          const [invRe, invIm] = cdiv(1, 0, r.muRe, r.muIm);
          const re = sr - invRe, im = w - invIm;
          const lamMag = Math.hypot(re, im);
          if (!Number.isFinite(re) || !Number.isFinite(im) || lamMag > rCeiling) continue;
          const deltaLam = r.ritzEst / (muMag * muMag);
          if (im < -1e-10 * Math.max(1, Math.abs(re))) continue;
          candidates.push({ re, im, vr: r.vr, vi: r.vi, ritzEst: r.ritzEst, deltaLam, shift: w });
        }
      }

      if (succeeded.length === 0) {
        return {
          verdict: 'undecided',
          reason: 'no shifted linearisation factored without singularity',
          drift: reportedDrift,
          modes: [],
          deflated: deflation.meta,
          shifts: { attempted, succeeded },
          window: null,
          warnings,
          ...(debug ? { debug: debugData } : {}),
        };
      }

      // 5. Pre-deduplicate Ritz candidates, refine rightmost ones, compute left-vector error bars.
      candidates.sort((a, b) => b.re - a.re || Math.abs(a.im) - Math.abs(b.im));
      const uniqueCands = [];
      for (const c of candidates) {
        const dup = uniqueCands.find((u) => Math.hypot(u.re - c.re, u.im - c.im) <= 1e-3 * Math.max(1, Math.hypot(u.re, u.im), Math.hypot(c.re, c.im)));
        if (!dup) uniqueCands.push(c);
        else if (c.deltaLam < dup.deltaLam) {
          dup.re = c.re; dup.im = c.im; dup.vr = c.vr; dup.vi = c.vi; dup.ritzEst = c.ritzEst; dup.deltaLam = c.deltaLam;
        }
      }

      const refined = [];
      let acceptedSlots = 0;
      const targetSlots = Math.max(nWant + 4, 12);
      const maxAttempts = Math.min(uniqueCands.length, Math.max(nWant * 12, 160));
      let rejectedSpurious = 0;
      for (let i = 0; i < maxAttempts; i++) {
        if (acceptedSlots >= targetSlots) {
          // Since uniqueCands is sorted by descending Re λ, stop once remaining candidates lie
          // strictly left of our current nWant-th accepted mode.
          let worst = Infinity;
          for (const u of refined) if (u.re < worst) worst = u.re;
          if (uniqueCands[i].re < worst - Math.max(1, Math.abs(worst)) * 1e-2) break;
        }
        const res = refineAndBound(pencil, deflation, coordScale, uniqueCands[i]);
        if (!res || res.residual > 1e-5) {
          rejectedSpurious++;
          continue;
        }
        // Final deduplication after refinement (eigenvalues close, or close-ish with collinear eigenvectors).
        const dup = refined.find((u) => {
          const dist = Math.hypot(u.re - res.re, u.im - res.im);
          const scale = Math.max(1, Math.hypot(u.re, u.im), Math.hypot(res.re, res.im));
          return dist <= 1e-4 * scale + 10 * (u.error + res.error) || (dist <= 2e-2 * scale && vecOverlap(u, res) > 0.99);
        });
        if (!dup) {
          refined.push(res);
          acceptedSlots += Math.abs(res.im) > 1e-10 * Math.max(1, Math.abs(res.re)) ? 2 : 1;
        } else if (res.residual < dup.residual) {
          Object.assign(dup, res);
        }
      }

      // Expand complex conjugate pairs (re ± i im) for reporting, sorted rightmost first.
      refined.sort((a, b) => b.re - a.re || Math.abs(a.im) - Math.abs(b.im));
      const allModes = [];
      for (const m of refined) {
        const base = {
          re: m.re,
          im: m.im,
          residual: m.residual,
          error: m.error,
          ...(profiles ? { profiles: modeProfiles(this, pencil, m.vr, m.vi) } : {}),
        };
        allModes.push(base);
        if (Math.abs(m.im) > 1e-10 * Math.max(1, Math.abs(m.re))) {
          const vic = Float64Array.from(m.vi, (v) => -v);
          allModes.push({
            re: m.re,
            im: -m.im,
            residual: m.residual,
            error: m.error,
            ...(profiles ? { profiles: modeProfiles(this, pencil, m.vr, vic) } : {}),
          });
        }
        if (allModes.length >= nWant) break;
      }
      const modesOut = allModes.slice(0, nWant);

      // 6. Compute honest window and verdict (never a bare 'stable').
      const rHi = Math.max(...succeeded.filter((x) => x.im === 0).map((x) => x.re), 0);
      const rLo = modesOut.length > 0 ? Math.max(...modesOut.map((m) => kBar * m.error)) : EPS * Math.max(1, rHi);
      const window = { rLo, rHi, omegaHi: maxOmega };

      const growing = modesOut.find((m) => m.re > kBar * m.error);
      const neutral = modesOut.find((m) => Math.abs(m.re) <= kBar * m.error);
      let verdict, reason;
      if (growing) {
        verdict = 'unstable';
        reason = `growing mode at Re λ = ${growing.re.toPrecision(4)} ± ${growing.error.toPrecision(2)} /s` +
          (Math.abs(growing.im) > 0 ? ` (ω = ${Math.abs(growing.im).toPrecision(4)} rad/s)` : '');
      } else if (neutral) {
        verdict = 'neutral within resolution';
        reason = `mode at Re λ = ${neutral.re.toPrecision(3)} /s lies within ${kBar}× its error bar (${neutral.error.toPrecision(2)} /s)`;
      } else {
        verdict = 'no growing mode found';
        reason = `no growing mode found in [${rLo.toPrecision(2)}, ${rHi.toPrecision(2)}] /s` +
          (maxOmega > 0 ? ` (|Im λ| ≤ ${maxOmega.toPrecision(3)} rad/s)` : ' (real shifts only)');
      }

      return {
        verdict,
        reason,
        drift: reportedDrift,
        modes: modesOut,
        deflated: deflation.meta,
        rejectedSpurious,
        shifts: { attempted, succeeded },
        window,
        warnings,
        ...(debug ? { debug: debugData } : {}),
      };
    } finally {
      this.allTerminals = wasAll;
      this.steady = wasSteady;
      this.lawForm = wasLawForm;
      this.lostMode = learnedKeep.lostMode;
      this.accepted = learnedKeep.accepted;
      this.firstErrs = learnedKeep.firstErrs;
      this._restore(snap);
    }
  }
}
