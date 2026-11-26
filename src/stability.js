// Linear stability of a converged steady state: the rightmost (slowest-decaying or growing)
// eigenvalues of the dynamic linearised DAE pencil −J v = λ M v about the device's current state,
// with matrix-verified deflation of conserved amounts, scale-invariant iterative refinement,
// pole verification against infinite DAE ghosts, left-eigenvector + residual error bars, and an
// honest windowed verdict.
//
// CONTRACT
//
// Methods of the Solver (src/solver.js; its layout is in docs/internals.md), installed on its
// prototype like SmallSignal. Device.stability (src/index.js) calls stability() on the device's
// current, converged steady state without re-solving (preserving the branch).
//
// What it uses of the rest:
//   this._assembleSplit() (src/assemble.js)            steady Jacobian J (this.sys), storage
//                                                      Jacobian M (this.sysMass), dynamic steady
//                                                      residuals (this.res, this.termRes), and
//                                                      terminal border terms on the dynamic rows
//                                                      (steady false, frozen sourceTime = time).
//   this._assembleSplitDifference() (src/difference.js) cancellation-resistant J·v and M·v for
//                                                      held-vector residual checks.
//   this._currentCuts() (src/impedance.js)             cut-current readings for floating contacts
//                                                      so extreme face capacitances or metal G
//                                                      don't cancel the series capacitance in the
//                                                      Schur complement.
//   this._sgSegment(), this._fluxConcentrated() (src/regions.js)
//                                                      steady conduction current across a cut
//                                                      segment (for cancellation-free cut residual).
//   this._smallSignalProfiles() (src/impedance.js)     complex mode profiles (δφ, δμ̄, δc).
//   maxPotentialStep (src/newton.js)                   thermal-unit state-motion norm for the
//                                                      steadiness check.
//   lawAmount, sharedChange (src/lawRows.js), chargeOnly (src/laws.js)
//                                                      seeds for exact left null vectors w of the
//                                                      assembled J (wᵀ J = 0, c = Mᵀ w).
//   BlockTridiagonal, ComplexBlockTridiagonal (src/blockTridiagonal.js)
//                                                      forward and transposed (solveT / solveH)
//                                                      solves for shift-invert Arnoldi, left-null
//                                                      polish, inverse iteration, and left vectors.
//   cdiv (src/impedance.js)                            overflow-safe complex division.
//   sourceAt (src/sources.js)                          frozen terminal source at this.time.
//
// stability(opts) → StabilityResult
//   opts:
//     shifts?: number[]         real positive shifts σ (1/s); default log-spaced over the device's
//                               own rates (scale-invariant), never 0.
//     frequencies?: number[]    optional frequencies f (Hz) for complex shifts σ = σ_r + 2πi f
//                               (for Hopf pairs far from the real axis).
//     modes?: number            number of rightmost modes to return (default 8).
//     profiles?: boolean        attach complex mode shapes { phi, mu, c, termV } (default false).
//     kBar?: number             positive finite significance multiple of the error bar (default 10).
//     debug?: boolean           include raw { J, M, deflation, localDeflated, cutUsed, labels }
//                               as plain arrays for independent cross-checks (default false).
//
// Precondition: the solver must hold a converged, un-regridded steady state (atSteady &&
// !regriddedSteady), else SolverError. On the assembled cut-current dynamic pencil (frozen at
// this.time, never mutating time-stepper state), we check (a) row-scaled dynamic residual ≤ 1e-6
// and (b) the linearised backward-Euler step δ = (J + σ M)⁻¹ F (projected onto the conserved
// hyperplane, measured in thermal units via maxPotentialStep) at solve()'s own giant-step rate
// σ = 1/(1e6·τ_slow) (tested at two non-resonant rates σ and σ/√10 to avoid an accidental σ = λ
// pole): if drift ≥ 1e-5 thermal units, stability() returns verdict 'undecided'.
//
// Deflation: every global conserved combination in constraints (whose terminals are I-driven) and
// every I-driven chargeOnly terminal (seeded from raw contact storage before cut replacement) is
// polished via transposed shift-invert (J + σ_z M)⁻ᵀ Mᵀ on the assembled pencil until it satisfies
// wᵀ J = 0 to round-off (checked directly with mulJT against |J|ᵀ |w|), then c = Mᵀ w is taken
// from the assembled M. Node-local invariants (localConstraints) are projected per block in O(N).
// The hyperplane c·x = 0 is invariant under T_σ = (J + σ M)⁻¹ M (c·T_σ x = c·x / σ).
//
// Refinement, anti-ghost verification, residuals, error bars, and deduplication:
//   - Each Ritz candidate must lie within a Meerbergen–Roose trust region of the shift that found
//     it (|λ − σ| ≤ 25 |σ|).
//   - Scale-invariant iterative two-sided inverse iteration (relative shift offset 1e-6 |λ|, up to
//     5 re-centered passes) must converge (|Δλ| ≤ 1e-5 |λ|) with difference-form row-scaled
//     residual ≤ 1e-8.
//   - Every growing, far-from-shift, or looser-residual candidate also passes a 0.1%-offset pole
//     check at s_p = λ (1 + 1e-3): one solve u = (J + s_p M)⁻¹ M v must reproduce v / (s_p − λ),
//     separating genuine finite poles from infinite DAE eigenvalue ghosts.
//   - Each accepted mode reports a first-order error bar covering both rounding and residual:
//       error = (ε · |w|ᴴ (|J| + |λ| |M|) |v| + |wᴴ r| + ρ_L · |w|ᴴ |r|) / |wᴴ M v|
//     where ρ_L is the relative left-eigenvector residual.
//   - Duplicates are merged only by eigenvalue proximity within bars AND biorthogonality
//     (w₁ᴴ M v₂ ≠ 0), never by Euclidean eigenvector overlap and never merging growing with
//     decaying modes.
//
// Verdict (evaluated over all verified modes inside the searched window; never bare 'stable'):
//   'unstable'                 a verified mode inside the searched window has Re λ > kBar · error;
//   'undecided'                precondition failed, no shift factored, unkept conserved amounts
//                              remain undeflated (when no growing mode was found), or rLo ≥ rHi;
//   'neutral within resolution' no growing mode, and a verified mode has |Re λ| ≤ kBar · error;
//   'no growing mode found'    all verified modes have Re λ < −kBar · error, qualified by
//                              { rLo, rHi, omegaHi }.
//
// Invariants:
//   1. Completely non-destructive: STATE, LEARNED, allTerminals, steady, lawForm, sourceBefore,
//      and offNeutral are identical before and after the call, even on error or undecided.
//   2. Dynamic rows only: assembly always runs with this.steady = false; steady law pins, flat
//      rows, affinities and island sums never enter J or M.
//   3. Deterministic & scale-invariant: fixed PRNG seed; scaling all rates by s scales shifts and
//      eigenvalues by s with the same verdict.
//   4. Pure JS, zero dependencies, block-tridiagonal + tiny K×K Schur border only, no ** or
//      Math.pow.
//
// Tests: test/stability.test.js.

import { BlockTridiagonal, ComplexBlockTridiagonal } from './blockTridiagonal.js';
import { FARADAY } from './constants.js';
import { SolverError } from './errors.js';
import { cdiv } from './impedance.js';
import { chargeOnly } from './laws.js';
import { lawAmount, sharedChange } from './lawRows.js';
import { maxPotentialStep } from './newton.js';
import { sourceAt } from './device.js';

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
      const a = H[l][l], b = H[l][nn], c = H[nn][l], d = H[nn][nn];
      const df = 0.5 * (a - d), disc = df * df + b * c;
      if (disc >= 0) {
        const z = df >= 0 ? Math.sqrt(disc) : -Math.sqrt(disc);
        const p = df + z;
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
      for (let j = k; j < m; j++) {
        let v = H[k][j] + wy * H[k + 1][j];
        if (notLast) {
          v += wz * H[k + 2][j];
          H[k + 2][j] -= v * u;
        }
        H[k + 1][j] -= v * q;
        H[k][j] -= v * p;
      }
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

  // Only extract converged blocks (indices > nn if the iteration cap was hit).
  const startIdx = nn >= 0 ? nn + 1 : 0;
  const eigs = [];
  for (let i = startIdx; i < m; ) {
    if (i + 1 < m && H[i + 1][i] !== 0) {
      const a = H[i][i], b = H[i][i + 1], c = H[i + 1][i], d = H[i + 1][i + 1];
      const re = 0.5 * (a + d), df = 0.5 * (a - d);
      const disc = -(df * df + b * c);
      const im = disc > 0 ? Math.sqrt(disc) : 0;
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
    for (let i = l; i <= nn; i++) { Hr[i][i] -= sr; Hi[i][i] -= si; }
    const cs = new Float64Array(nn - l), snr = new Float64Array(nn - l), sni = new Float64Array(nn - l);
    for (let k = l; k < nn; k++) {
      const ar = Hr[k][k], ai = Hi[k][k], br = Hr[k + 1][k], bi = Hi[k + 1][k];
      const am = Math.hypot(ar, ai), bm = Math.hypot(br, bi), r = Math.hypot(am, bm);
      if (r === 0) { cs[k - l] = 1; snr[k - l] = 0; sni[k - l] = 0; continue; }
      const c = am / r;
      let svr, svi;
      if (am > 0) {
        const er = ar / am, ei = ai / am;
        svr = (er * br + ei * bi) / r;
        svi = (ei * br - er * bi) / r;
      } else {
        svr = 1; svi = 0;
      }
      cs[k - l] = c; snr[k - l] = svr; sni[k - l] = svi;
      for (let j = k; j < m; j++) {
        const xr = Hr[k][j], xi = Hi[k][j], yr = Hr[k + 1][j], yi = Hi[k + 1][j];
        Hr[k][j] = c * xr + (svr * yr - svi * yi);
        Hi[k][j] = c * xi + (svr * yi + svi * yr);
        Hr[k + 1][j] = -(svr * xr + svi * xi) + c * yr;
        Hi[k + 1][j] = -(svr * xi - svi * xr) + c * yi;
      }
      Hr[k + 1][k] = 0; Hi[k + 1][k] = 0;
    }
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

  const startIdx = nn >= 0 ? nn + 1 : 0;
  const eigs = [];
  for (let i = startIdx; i < m; i++) {
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

function btAbsMulT(bt, absW, out) {
  out.fill(0);
  const { n, A, B, C, sizes, offA, offB, offC, offX } = bt;
  for (let blk = 0; blk < n; blk++) {
    const m = sizes[blk], v = offX[blk];
    const mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < n - 1 ? sizes[blk + 1] : 0;
    for (let r = 0; r < m; r++) {
      const wr = absW[v + r];
      if (wr === 0) continue;
      for (let c = 0; c < mp; c++) out[v - mp + c] += Math.abs(A[offA[blk] + r * mp + c]) * wr;
      for (let c = 0; c < m; c++) out[v + c] += Math.abs(B[offB[blk] + r * m + c]) * wr;
      for (let c = 0; c < mn; c++) out[v + m + c] += Math.abs(C[offC[blk] + r * mn + c]) * wr;
    }
  }
  return out;
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

// Steady conduction current across interior segment cutSeg (toward +x, A or A/m²).
// Displacement ∂D/∂t is identically zero at steady state.
function steadyCutCurrent(s, cutSeg) {
  const { model, n, z } = s;
  const { grid, materials, regions } = model;
  let rIdx = 0;
  for (let r = 0; r < regions.length; r++) {
    if (cutSeg >= grid.regionStart[r] && cutSeg < grid.regionEnd[r]) { rIdx = r; break; }
  }
  const reg = regions[rIdx], mat = materials[reg.material];
  if (mat.conductor) return 0;
  const bL = cutSeg + rIdx, bR = bL + 1, h = grid.segLength[cutSeg] / grid.segArea[cutSeg];
  const vel = reg.velocity;
  let Icut = 0;
  for (let i = 0; i < n; i++) {
    if (z[i] === 0 || !mat.present[i] || mat.D[i] === 0) continue;
    if (!mat.ideal) {
      const { N } = s._fluxConcentrated(cutSeg, bL, bR, mat, h, vel, i);
      Icut += z[i] * FARADAY * N;
    } else {
      const sg = s._sgSegment(cutSeg, bL, i, h, vel);
      Icut += z[i] * FARADAY * sg[0];
    }
  }
  return Icut;
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

  // Keep raw contact storage rows before cut replacement (for held-charge deflation seeding).
  const rawCs = s.termCs.map((v) => v.slice(0, Nc));
  const rawDIs = Float64Array.from(s.termDIs);
  const rawDXs = s.termDXs.map((v) => Float64Array.from(v));

  const B = fl.map((k) => s.termB[k].slice(0, Nc));
  const C = fl.map((k) => s.termC[k].slice(0, Nc));
  const Cs = fl.map((k) => s.termCs[k].slice(0, Nc));
  const D = fl.map((k, a) => Float64Array.from(fl, (k2, b) => (a === b ? s.termDI[k] : s.termDX[k][k2])));
  const Ds = fl.map((k, a) => Float64Array.from(fl, (k2, b) => (a === b ? s.termDIs[k] : s.termDXs[k][k2])));

  // Sol's cut-current reading for floating contacts not read by a circuit law's change.
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

  // Dynamic steady residual vector F (frozen at s.time, with cancellation-free cut / lawForm rows).
  const F = new Float64Array(N);
  F.set(s.res.subarray(0, Nc), 0);
  fl.forEach((k, a) => {
    if (cutUsed[a] >= 0) {
      const sgn = k === 0 ? 1 : -1;
      const Icut = steadyCutCurrent(s, cutUsed[a]);
      const d = s.terms[k].drive;
      const src = s.sourceOverride.has(k) ? s.sourceOverride.get(k) : sourceAt(d.src, s.sourceTime, s.sourceBefore);
      F[Nc + a] = d.kind === 'I' ? sgn * Icut - src : sgn * Icut - (src - s.termV[k]) / d.R;
    } else if (lawRead.has(k)) {
      // circuitLaws in split mode replaces C/D/Cs/Ds with the by-change law without overwriting
      // termRes; at u = uOld (amount === lawStart), that row's residual is purely −I_drive.
      const d = s.terms[k].drive;
      const src = s.sourceOverride.has(k) ? s.sourceOverride.get(k) : sourceAt(d.src, s.sourceTime, s.sourceBefore);
      F[Nc + a] = -(d.kind === 'I' ? src : (src - s.termV[k]) / d.R);
    } else {
      F[Nc + a] = s.termRes[k];
    }
  });

  // Row scales of J and M separately, plus dynamic rate samples jRowMax/mRowMax.
  const jRowMax = new Float64Array(N), mRowMax = new Float64Array(N);
  let jInf = 0;
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
      if (jm > jInf) jInf = jm;
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
    if (jm > jInf) jInf = jm;
  }

  // Row-rate range from differential rows (mRowMax > 0 && jRowMax > 0).
  const rowRates = [];
  for (let i = 0; i < N; i++) {
    if (mRowMax[i] > 0 && jRowMax[i] > 0) {
      const r = jRowMax[i] / mRowMax[i];
      if (Number.isFinite(r) && r > 0) rowRates.push(r);
    }
  }
  rowRates.sort((a, b) => a - b);
  const rMinDyn = rowRates.length ? rowRates[0] : 1;
  const rMaxDyn = rowRates.length ? rowRates[rowRates.length - 1] : 1e6;

  // Matvec helpers.
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
  const mulJT = (w, out) => {
    btMulT(Jg, w, out);
    for (let j = Nc; j < N; j++) out[j] = 0;
    for (let a = 0; a < K; a++) {
      let dotB = 0;
      const Ba = B[a];
      for (let r = 0; r < Nc; r++) dotB += Ba[r] * w[r];
      out[Nc + a] += dotB;
    }
    for (let a = 0; a < K; a++) {
      const wa = w[Nc + a];
      if (wa === 0) continue;
      const Ca = C[a], Da = D[a];
      for (let j = 0; j < Nc; j++) out[j] += Ca[j] * wa;
      for (let b = 0; b < K; b++) out[Nc + b] += Da[b] * wa;
    }
    return out;
  };
  const absMulJT = (absW, out) => {
    btAbsMulT(Jg, absW, out);
    for (let j = Nc; j < N; j++) out[j] = 0;
    for (let a = 0; a < K; a++) {
      let dotB = 0;
      const Ba = B[a];
      for (let r = 0; r < Nc; r++) dotB += Math.abs(Ba[r]) * absW[r];
      out[Nc + a] += dotB;
    }
    for (let a = 0; a < K; a++) {
      const wa = absW[Nc + a];
      if (wa === 0) continue;
      const Ca = C[a], Da = D[a];
      for (let j = 0; j < Nc; j++) out[j] += Math.abs(Ca[j]) * wa;
      for (let b = 0; b < K; b++) out[Nc + b] += Math.abs(Da[b]) * wa;
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

  const absMulJMT = (absW, lamAbs, out) => {
    btAbsMulT(Jg, absW, out);
    for (let j = Nc; j < N; j++) out[j] = 0;
    if (lamAbs > 0) {
      const mTmp = new Float64Array(N);
      btAbsMulT(Mg, absW, mTmp);
      for (let j = 0; j < Nc; j++) out[j] += lamAbs * mTmp[j];
    }
    for (let a = 0; a < K; a++) {
      let dotB = 0;
      const Ba = B[a];
      for (let r = 0; r < Nc; r++) dotB += Math.abs(Ba[r]) * absW[r];
      out[Nc + a] += dotB;
    }
    for (let a = 0; a < K; a++) {
      const wa = absW[Nc + a];
      if (wa === 0) continue;
      const Ca = C[a], Csa = Cs[a], Da = D[a], Dsa = Ds[a];
      for (let j = 0; j < Nc; j++) out[j] += (Math.abs(Ca[j]) + lamAbs * Math.abs(Csa[j])) * wa;
      for (let b = 0; b < K; b++) out[Nc + b] += (Math.abs(Da[b]) + lamAbs * Math.abs(Dsa[b])) * wa;
    }
    return out;
  };

  return {
    Jg, Mg, Jr, Sr, lin, Nc, K, N, fl, B, C, Cs, D, Ds, F,
    rawCs, rawDIs, rawDXs, lawRead, cutUsed, jRowMax, mRowMax, jInf, rMinDyn, rMaxDyn,
    mulJDiff, mulMDiff, mulM, mulMT, mulJT, absMulJT, absMulJM, absMulJMT,
  };
}

// --- Bordered shifted factorisations (real & complex) with forward + transposed solves ---

function factorRealShift(pencil, sigma, strict = false) {
  const { Jg, Mg, Nc, K, N, B, C, Cs, D, Ds } = pencil;
  const { sizes, offA, offB, offC, offX } = Jg, nB = Jg.n;
  const bt = new BlockTridiagonal(nB, sizes);
  bt.bothWays = true;
  bt.staticPivots = true;
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
  if (strict && bt.lostColumns > 0) return null;
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

// --- Matrix-verified deflation: exact wᵀ J = 0, c = Mᵀ w on the assembled pencil ---

function leftNullResidual(pencil, w) {
  const { N, jInf } = pencil;
  const jt = new Float64Array(N), absW = new Float64Array(N), absJT = new Float64Array(N);
  let wInf = 0;
  for (let j = 0; j < N; j++) { const a = Math.abs(w[j]); absW[j] = a; if (a > wInf) wInf = a; }
  if (!(wInf > 0)) return Infinity;
  pencil.mulJT(w, jt);
  pencil.absMulJT(absW, absJT);
  let maxJT = 0, maxAbs = 0;
  for (let j = 0; j < N; j++) {
    const a = Math.abs(jt[j]);
    if (a > maxJT) maxJT = a;
    if (absJT[j] > maxAbs) maxAbs = absJT[j];
  }
  const denom = Math.max(maxAbs, 1e-12 * jInf * wInf);
  return denom > 0 ? maxJT / denom : 0;
}

function buildDeflation(s, pencil, coordScale, warnings) {
  const { Nc, K, N, fl, rawCs, rawDIs, rawDXs, rMinDyn } = pencil;
  const { model, n, M, th, nSurf } = s;
  const vol = model.grid.vol;
  const invScale = Float64Array.from(coordScale, (v) => 1 / v);
  const seeds = [];
  const meta = [];

  // 1. Global conserved combinations from s.constraints (all feeding terminals I-driven).
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
    seeds.push({ kind: 'constraint', vec, label: `conserved: ${names.join(' + ') || cs.key}` });
  }

  // 2. Floating I-driven chargeOnly terminals: seed from raw (pre-cut) storage rows.
  fl.forEach((k, a) => {
    if (s.terms[k].drive.kind !== 'I' || !chargeOnly(s, k) || handledCircuit.has(k)) return;
    const vec = new Float64Array(N);
    vec.set(rawCs[k], 0);
    fl.forEach((k2, b) => { vec[Nc + b] = k2 === k ? rawDIs[k] : rawDXs[k][k2]; });
    seeds.push({ kind: 'charge', vec, label: `held charge: terminal ${s.terms[k].name}` });
  });

  // 3. Node-local invariants (immobile species / trap pools), stored sparsely by block for O(N) projection.
  const localByBlock = new Map();
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
      let n0 = 0;
      for (let r = 0; r < m; r++) { locVec[r] *= invScale[off + r]; n0 += locVec[r] * locVec[r]; }
      n0 = Math.sqrt(n0);
      if (!(n0 > 0)) continue;
      const list = localByBlock.get(blk) ?? (localByBlock.set(blk, []), localByBlock.get(blk));
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

  // Helper: project a y-space vector orthogonal to localByBlock and a list of orthonormal global vectors.
  const projectY = (v, qList) => {
    for (let pass = 0; pass < 2; pass++) {
      for (const [blk, list] of localByBlock) {
        const off = pencil.Jg.offX[blk], m = pencil.Jg.sizes[blk];
        for (const qv of list) {
          let dot = 0;
          for (let r = 0; r < m; r++) dot += qv[r] * v[off + r];
          for (let r = 0; r < m; r++) v[off + r] -= dot * qv[r];
        }
      }
      for (const qv of qList) {
        let dot = 0;
        for (let j = 0; j < N; j++) dot += qv[j] * v[j];
        for (let j = 0; j < N; j++) v[j] -= dot * qv[j];
      }
    }
  };

  // Polish each global seed via transposed shift-invert (J + σ_z M)⁻ᵀ Mᵀ to obtain an exact
  // left null vector w of the ASSEMBLED J (wᵀ J = 0 to round-off) and set c = Mᵀ w from M.
  const globalQ = [];
  const debugVecs = [];
  if (seeds.length > 0) {
    const tauSlow = s.slowestTime();
    const rBase = Math.min(Number.isFinite(tauSlow) && tauSlow > 0 ? 1 / tauSlow : rMinDyn, rMinDyn);
    const trySigs = [rBase * 1e-6, rBase * 1e-9, rBase * 1e-3, rBase * 1e-12];
    const facs = [];
    for (const sig of trySigs) {
      if (!(sig > 0 && Number.isFinite(sig))) continue;
      const f = factorRealShift(pencil, sig);
      if (f) { facs.push({ sig, fac: f }); if (facs.length >= 2) break; }
    }
    const cTmp = new Float64Array(N), wTmp = new Float64Array(N), yTmp = new Float64Array(N);
    for (const item of seeds) {
      let best = null;
      for (const { fac } of facs) {
        // Seed w from (J + σ M)⁻ᵀ c_seed.
        let nSeed = 0;
        for (let j = 0; j < N; j++) nSeed += item.vec[j] * item.vec[j];
        if (!(nSeed > 0)) break;
        fac.solveT(item.vec, wTmp);
        let nw = 0;
        for (let j = 0; j < N; j++) nw += wTmp[j] * wTmp[j];
        nw = Math.sqrt(nw);
        if (!(nw > 0) || !Number.isFinite(nw)) continue;
        for (let j = 0; j < N; j++) wTmp[j] /= nw;

        let relJ = leftNullResidual(pencil, wTmp);
        for (let it = 0; it < 6 && relJ > 1e-13; it++) {
          pencil.mulMT(wTmp, cTmp);
          // Orthogonalize cTmp in y-space against already-accepted directions so independent
          // constraints converge to distinct left null vectors.
          for (let j = 0; j < N; j++) yTmp[j] = cTmp[j] * invScale[j];
          projectY(yTmp, globalQ);
          let ny = 0;
          for (let j = 0; j < N; j++) ny += yTmp[j] * yTmp[j];
          ny = Math.sqrt(ny);
          if (!(ny > 0) || !Number.isFinite(ny)) { relJ = Infinity; break; }
          for (let j = 0; j < N; j++) cTmp[j] = (yTmp[j] / ny) * coordScale[j];
          fac.solveT(cTmp, wTmp);
          nw = 0;
          for (let j = 0; j < N; j++) nw += wTmp[j] * wTmp[j];
          nw = Math.sqrt(nw);
          if (!(nw > 0) || !Number.isFinite(nw)) { relJ = Infinity; break; }
          for (let j = 0; j < N; j++) wTmp[j] /= nw;
          relJ = leftNullResidual(pencil, wTmp);
        }
        if (relJ < 1e-9 && (!best || relJ < best.relJ)) {
          const cExact = pencil.mulMT(wTmp, new Float64Array(N));
          best = { relJ, c: cExact, w: wTmp.slice() };
          if (relJ <= 1e-13) break;
        }
      }
      if (!best) {
        warnings.push(`stability: could not verify exact left null vector for "${item.label}"; leaving it undeflated.`);
        continue;
      }
      const v = new Float64Array(N);
      let n0 = 0;
      for (let j = 0; j < N; j++) { v[j] = best.c[j] * invScale[j]; n0 += v[j] * v[j]; }
      n0 = Math.sqrt(n0);
      if (!(n0 > 0)) continue;
      projectY(v, globalQ);
      let n1 = 0;
      for (let j = 0; j < N; j++) n1 += v[j] * v[j];
      n1 = Math.sqrt(n1);
      if (n1 <= 1e-12 * n0) continue;
      for (let j = 0; j < N; j++) v[j] /= n1;
      globalQ.push(v);
      debugVecs.push({ label: item.label, vector: Array.from(best.c), wResidual: best.relJ });
      meta.push({ kind: item.kind, label: item.label });
    }
  }

  const totalDeflated = globalQ.length + localCount;
  const projectReal = (y) => {
    if (totalDeflated === 0) return y;
    projectY(y, globalQ);
    return y;
  };
  const projectComplex = (yr, yi) => {
    if (totalDeflated === 0) return;
    projectReal(yr);
    projectReal(yi);
  };

  return { totalDeflated, localCount, meta, debugVecs, localByBlock, projectReal, projectComplex };
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
    if (!(normW0 > 0) || !Number.isFinite(normW0)) break;
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
    if (!(hn > 1e-14 * normW0)) {
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
  const v0r = Float64Array.from(startVec), v0i = new Float64Array(N);
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
    if (!(n0 > 0) || !Number.isFinite(n0)) break;
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
    if (!(hn > 1e-14 * n0)) { k++; break; }
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

// --- Refinement, anti-ghost pole check, held-vector difference-form residual, and error bar ---

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
  const rrArr = new Float64Array(N), riArr = new Float64Array(N);
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
  if (!(vInf > 0)) return { residual: Infinity, jvr, jvi, mvr, mvi, rrArr, riArr };
  let resInf = 0;
  for (let i = 0; i < N; i++) {
    const rr = jvr[i] + lamRe * mvr[i] - lamIm * mvi[i];
    const ri = jvi[i] + lamRe * mvi[i] + lamIm * mvr[i];
    rrArr[i] = rr;
    riArr[i] = ri;
    const sc = jRowMax[i] + lamAbs * mRowMax[i];
    if (sc > 0) {
      const r = Math.hypot(rr, ri) / (sc * vInf);
      if (r > resInf) resInf = r;
    }
  }
  return { residual: resInf, jvr, jvi, mvr, mvi, rrArr, riArr };
}

// Biorthogonal overlap between two refined modes:
//   bio(a, b) = sqrt(|a.wᴴ M b.v| · |b.wᴴ M a.v| / (|a.wᴴ M a.v| · |b.wᴴ M b.v|)).
// Equals ~1 for two copies of the same simple mode, and ~0 for distinct modes even in strongly
// non-normal pencils where unweighted Euclidean |v₁ᴴ v₂| can exceed 0.999.
function biorthOverlap(a, b) {
  const N = a.wr.length;
  let d1r = 0, d1i = 0, d2r = 0, d2i = 0;
  for (let j = 0; j < N; j++) {
    d1r += a.wr[j] * b.mvr[j] + a.wi[j] * b.mvi[j];
    d1i += a.wr[j] * b.mvi[j] - a.wi[j] * b.mvr[j];
    d2r += b.wr[j] * a.mvr[j] + b.wi[j] * a.mvi[j];
    d2i += b.wr[j] * a.mvi[j] - b.wi[j] * a.mvr[j];
  }
  const prod = Math.hypot(d1r, d1i) * Math.hypot(d2r, d2i);
  const denom = a.denAbs * b.denAbs;
  return denom > 0 ? Math.sqrt(prod / denom) : 0;
}

// 0.1%-offset pole check: for a true finite eigenvalue λ with eigenvector v,
// (J + s_p M)⁻¹ M v = v / (s_p − λ) at s_p = λ (1 + 1e-3). An infinite-eigenvalue ghost has no
// pole at λ, so (s_p − λ) u differs from v by O(1) and its Rayleigh quotient leaves λ completely.
function verifyPole(pencil, isReal, lamRe, lamIm, vr, vi, mvr, mvi) {
  const lamMag = Math.hypot(lamRe, lamIm);
  if (!(lamMag > 0)) return true;
  const epsP = 1e-3;
  const sRe = lamRe * (1 + epsP), sIm = lamIm * (1 + epsP);
  const dRe = sRe - lamRe, dIm = sIm - lamIm;
  const { N } = pencil;
  let ur, ui;
  if (isReal) {
    const fac = factorRealShift(pencil, sRe);
    if (!fac) return false;
    ur = fac.solve(mvr);
    ui = new Float64Array(N);
  } else {
    const fac = factorComplexShift(pencil, sRe, sIm);
    if (!fac) return false;
    [ur, ui] = fac.solve(mvr, mvi);
  }
  // Check (s_p − λ) u ≈ v and Rayleigh estimate λ_probe = s_p − (vᴴ v)/(vᴴ u) ≈ λ.
  let diff2 = 0, v2 = 0, dotR = 0, dotI = 0;
  for (let j = 0; j < N; j++) {
    const suR = dRe * ur[j] - dIm * ui[j] - vr[j];
    const suI = dRe * ui[j] + dIm * ur[j] - vi[j];
    diff2 += suR * suR + suI * suI;
    v2 += vr[j] * vr[j] + vi[j] * vi[j];
    dotR += vr[j] * ur[j] + vi[j] * ui[j];
    dotI += vr[j] * ui[j] - vi[j] * ur[j];
  }
  if (!(v2 > 0)) return false;
  const vecErr = Math.sqrt(diff2 / v2);
  if (!(vecErr <= 0.1)) return false;
  const muMag = Math.hypot(dotR, dotI) / v2;
  if (!(muMag > 0)) return false;
  const [invR, invI] = cdiv(v2, 0, dotR, dotI);
  const jump = Math.hypot(sRe - invR - lamRe, sIm - invI - lamIm);
  return jump <= 1e-2 * lamMag;
}

function refineAndBound(pencil, deflation, coordScale, cand, shiftScale) {
  const { N, mulM, mulMT, absMulJM } = pencil;
  const invScale = Float64Array.from(coordScale, (v) => 1 / v);
  const vr = new Float64Array(N), vi = new Float64Array(N);
  for (let j = 0; j < N; j++) { vr[j] = cand.vr[j] * invScale[j]; vi[j] = cand.vi[j] * invScale[j]; }
  let lamRe = cand.re, lamIm = cand.im;
  const isReal = Math.abs(lamIm) <= 1e-10 * Math.abs(lamRe);
  if (isReal) { lamIm = 0; vi.fill(0); }
  if (!normalizeComplex(vr, vi)) return null;

  const wr = new Float64Array(N), wi = new Float64Array(N);

  const twoSidedPass = (sRe, sIm) => {
    const mag = Math.hypot(sRe, sIm);
    const bump = mag > 0 ? 1e-6 * mag : 1e-12 * shiftScale;
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
      for (let it = 0; it < 4; it++) {
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
      for (let it = 0; it < 4; it++) {
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

  // Iterative re-centered two-sided refinement until λ stops moving (up to 5 passes).
  let best = twoSidedPass(lamRe, lamIm);
  if (!best) return null;
  let lastStep = Math.hypot(best.re - lamRe, best.im - lamIm);
  let bestVr = vr.slice(), bestVi = vi.slice(), bestWr = wr.slice(), bestWi = wi.slice();
  for (let pass = 1; pass < 5; pass++) {
    const mag = Math.hypot(best.re, best.im);
    if (lastStep <= 1e-8 * Math.max(mag, 1e-30 * shiftScale) && best.check.residual <= 1e-11) break;
    const prevRe = best.re, prevIm = best.im;
    const next = twoSidedPass(prevRe, prevIm);
    if (!next) break;
    lastStep = Math.hypot(next.re - prevRe, next.im - prevIm);
    if (next.check.residual <= best.check.residual || lastStep <= 1e-6 * Math.max(Math.hypot(next.re, next.im), 1e-30 * shiftScale)) {
      best = next;
      bestVr = vr.slice(); bestVi = vi.slice(); bestWr = wr.slice(); bestWi = wi.slice();
    } else {
      break;
    }
  }
  vr.set(bestVr); vi.set(bestVi); wr.set(bestWr); wi.set(bestWi);
  lamRe = best.re;
  lamIm = best.im;
  const { check, denAbs } = best;
  const lamMag = Math.hypot(lamRe, lamIm);

  // Convergence & tight residual gate.
  if (!Number.isFinite(check.residual) || check.residual > 1e-8) return null;
  if (lastStep > 1e-5 * Math.max(lamMag, 1e-30 * shiftScale)) return null;

  // Anti-ghost 0.1%-offset pole verification on any growing, far-from-shift, or looser-residual candidate.
  if (lamRe > 0 || lamMag > 2 * cand.shiftMag || check.residual > 1e-10) {
    if (!verifyPole(pencil, isReal, lamRe, lamIm, vr, vi, check.mvr, check.mvi)) return null;
  }

  // First-order + second-order Rayleigh-quotient error bar:
  //   error = (ε · |w|ᴴ (|J| + |λ| |M|) |v| + |wᴴ r| + ρ_L · |w|ᴴ |r|) / |wᴴ M v|
  // where ρ_L = ‖(J + λ M)ᴴ w‖_∞ / ‖(|J| + |λ| |M|)ᵀ |w|‖_∞ measures left-vector convergence.
  const absV = new Float64Array(N), absW = new Float64Array(N), absJMv = new Float64Array(N);
  for (let j = 0; j < N; j++) { absV[j] = Math.hypot(vr[j], vi[j]); absW[j] = Math.hypot(wr[j], wi[j]); }
  absMulJM(absV, lamMag, absJMv);
  // Left residual (Jᵀ + conj(λ) Mᵀ) w
  const jtwr = pencil.mulJT(wr, new Float64Array(N));
  const mtwr = pencil.mulMT(wr, new Float64Array(N));
  const jtwi = isReal ? new Float64Array(N) : pencil.mulJT(wi, new Float64Array(N));
  const mtwi = isReal ? new Float64Array(N) : pencil.mulMT(wi, new Float64Array(N));
  const absJMTw = pencil.absMulJMT(absW, lamMag, new Float64Array(N));
  let maxRL = 0, maxAbsL = 0;
  for (let j = 0; j < N; j++) {
    // (J + λ M)ᴴ w = (Jᵀ + conj(λ) Mᵀ) w; conj(λ) = lamRe − i lamIm
    const lrr = jtwr[j] + lamRe * mtwr[j] + lamIm * mtwi[j];
    const lri = jtwi[j] + lamRe * mtwi[j] - lamIm * mtwr[j];
    const a = Math.hypot(lrr, lri);
    if (a > maxRL) maxRL = a;
    if (absJMTw[j] > maxAbsL) maxAbsL = absJMTw[j];
  }
  const rhoL = maxAbsL > 0 ? Math.min(1, maxRL / maxAbsL) : 1;

  let numRound = 0, numRes = 0, wrDotR = 0, wiDotR = 0;
  for (let j = 0; j < N; j++) {
    const wj = absW[j];
    numRound += wj * absJMv[j];
    numRes += wj * Math.hypot(check.rrArr[j], check.riArr[j]);
    wrDotR += wr[j] * check.rrArr[j] + wi[j] * check.riArr[j];
    wiDotR += wr[j] * check.riArr[j] - wi[j] * check.rrArr[j];
  }
  const wDotR = Math.hypot(wrDotR, wiDotR);
  const error = (EPS * numRound + wDotR + rhoL * numRes) / denAbs;
  if (!Number.isFinite(error)) return null;
  // A fast mode far from 0 whose error bar exceeds half its own |λ| has collapsed |wᴴ M v|
  // (a near-infinite DAE mode, not a resolved finite eigenvalue); slow/neutral modes near 0 keep
  // their error bars so genuinely neutral directions report 'neutral within resolution'.
  if (lamMag > 100 * shiftScale && error > 0.5 * lamMag) return null;

  return {
    re: lamRe, im: lamIm, residual: check.residual, error,
    vr, vi, wr, wi, mvr: check.mvr, mvi: check.mvi, denAbs,
  };
}

// --- Deterministic-phase profile extraction ---

function modeProfiles(s, pencil, vr, vi) {
  const { M, nB, VT } = s, { Nc, N, fl } = pencil;
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

// --- Scale-invariant default real shift grid from the device's own timescales (never σ = 0) ---

function hasAnyDiffusion(s) {
  for (const reg of s.model.regions) {
    const mat = s.model.materials[reg.material];
    if (mat.conductor) continue;
    for (let i = 0; i < s.n; i++) if (mat.present[i] && mat.D[i] > 0) return true;
  }
  return false;
}

function defaultShifts(s, pencil) {
  const { rMinDyn, rMaxDyn } = pencil;
  const hasDiff = hasAnyDiffusion(s);
  const tauSlow = s.slowestTime(), tauFast = s.fastestTime();
  let lo, hi;
  if (hasDiff && Number.isFinite(tauSlow) && tauSlow > 0 && Number.isFinite(tauFast) && tauFast > 0) {
    const rSlow = 1 / tauSlow, rFast = 1 / tauFast;
    // Widen toward actual pencil dynamic rates (RC, reactions, gates), bounded relative to diffusion
    // so an ultra-slow frozen reaction (1e-25 /s) doesn't stretch the grid over 30 decades.
    const dynLo = Math.max(rMinDyn, rSlow * 1e-3);
    const dynHi = Math.min(rMaxDyn, rFast * 1e3);
    lo = Math.min(rSlow, dynLo) * 0.1;
    hi = Math.max(rFast, dynHi) * 10;
  } else {
    lo = rMinDyn * 0.1;
    hi = Math.max(rMaxDyn * 10, lo * 1e4);
  }
  if (!(hi > lo)) hi = lo * 1e4;
  // Keep span sane (at most 14 decades from lo to hi) and floating-point finite.
  if (hi / lo > 1e14) lo = hi * 1e-14;
  const count = 8;
  const shifts = [];
  const logLo = Math.log(lo), logHi = Math.log(hi);
  for (let i = 0; i < count; i++) {
    shifts.push(Math.exp(logLo + ((logHi - logLo) * i) / (count - 1)));
  }
  return shifts;
}

// --- Linearised cut-current giant-step steadiness check (frozen at s.time, zero stepper leaks) ---

function checkDynamicSteady(s, pencil, deflation, coordScale) {
  const { Nc, K, N, F, jRowMax, rMinDyn } = pencil;
  const invScale = Float64Array.from(coordScale, (v) => 1 / v);

  // 1. Row-scaled dynamic steady residual (catches large bulk-reaction / terminal imbalances).
  let dynamicResInf = 0;
  for (let i = 0; i < N; i++) {
    const sc = jRowMax[i];
    const r = sc > 0 ? Math.abs(F[i]) / sc : Math.abs(F[i]);
    if (r > dynamicResInf) dynamicResInf = r;
  }
  if (dynamicResInf > 1e-6) {
    return {
      ok: false,
      drift: dynamicResInf,
      reason: `not a steady state of the time-dependent equations (row-scaled dynamic residual ${dynamicResInf.toExponential(2)})`,
    };
  }

  // 2. Linearised backward-Euler step δ = (J + σ M)⁻¹ F at solve()'s own giant-step rate
  //    σ = 1 / (1e6 · τ_slow), projected onto the conserved hyperplane and measured via
  //    maxPotentialStep (thermal units). Two non-resonant rates σ and 1.7 σ guard against an
  //    accidental σ = λ pole of a real growing mode.
  const hasDiff = hasAnyDiffusion(s);
  const tau = hasDiff ? s.slowestTime() : 1 / Math.max(rMinDyn, 1e-30);
  const giant = Number.isFinite(tau) && tau > 0 ? 1e6 * tau : 1e6;
  const sig1 = 1 / giant;
  const sig2 = sig1 * 1.7;

  const dx = new Float64Array(Nc + 1); // entry Nc is the sink (0) for maxPotentialStep
  const dV = new Float64Array(K);
  const y = new Float64Array(N);
  const evalStep = (sig) => {
    const fac = factorRealShift(pencil, sig);
    if (!fac) return null;
    const sol = fac.solve(F);
    // Project out conserved directions so O(ε/σ) round-off along exact zeros doesn't inflate drift.
    for (let j = 0; j < N; j++) y[j] = sol[j] * coordScale[j];
    deflation.projectReal(y);
    for (let j = 0; j < Nc; j++) dx[j] = y[j] * invScale[j];
    dx[Nc] = 0;
    for (let a = 0; a < K; a++) dV[a] = y[Nc + a] * invScale[Nc + a];
    const d = maxPotentialStep(s, dx, dV);
    return Number.isFinite(d) ? { drift: d, h: 1 / sig } : null;
  };
  let tries = [];
  for (const [a, b] of [[sig1, sig2], [sig1 * 100, sig2 * 100], [sig1 * 1e4, sig2 * 1e4]]) {
    const t1 = evalStep(a), t2 = evalStep(b);
    if (t1) tries.push(t1);
    if (t2) tries.push(t2);
    if (tries.length > 0) break;
  }
  if (tries.length === 0) {
    return {
      ok: false,
      drift: Infinity,
      reason: 'could not factor the dynamic linearisation for the steadiness check',
    };
  }
  // Taking the smaller of the two non-resonant steps avoids a single-point σ = λ resonance while
  // still seeing any genuine drift (which is O(1) at both σ and 1.7 σ).
  const chosen = tries.length === 2 && tries[1].drift < tries[0].drift ? tries[1] : tries[0];
  if (!(chosen.drift < 1e-5)) {
    return {
      ok: false,
      drift: chosen.drift,
      reason: `not a steady state of the time-dependent equations (a linearised step of ${chosen.h.toPrecision(2)} s moved the state by ${chosen.drift.toPrecision(3)} thermal units)`,
    };
  }
  return { ok: true, drift: chosen.drift, h: chosen.h };
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
      if (userShifts.length === 0 && frequencies.length === 0) {
        throw new SolverError('stability: provide at least one positive shift or frequency');
      }
    }
    for (const f of frequencies) {
      if (!(f > 0 && Number.isFinite(f))) throw new SolverError(`stability: frequencies must be positive and finite (Hz); got ${f}`);
    }
    if (!(Number.isInteger(nWant) && nWant >= 1)) throw new SolverError(`stability: modes must be a positive integer; got ${nWant}`);
    if (!(kBar > 0 && Number.isFinite(kBar))) throw new SolverError(`stability: kBar must be positive and finite; got ${kBar}`);

    const snap = this._snapshot();
    const wasAll = this.allTerminals, wasSteady = this.steady, wasLawForm = this.lawForm;
    const wasSourceBefore = this.sourceBefore, wasOffNeutral = this.offNeutral;
    const learnedKeep = { lostMode: this.lostMode, accepted: this.accepted, firstErrs: this.firstErrs?.slice() };
    const warnings = [];
    if (this.unkept.length > 0) {
      warnings.push(`stability: ${this.unkept.length} conserved combination(s) are in unkept; undeflated neutral directions remain.`);
    }
    const nonAuto = this.terms.some((t) => Boolean(t.drive.src?.t));
    if (nonAuto) {
      warnings.push(`stability: non-autonomous waveform drive frozen at t = ${this.time} s; instantaneous spectrum does not establish trajectory stability.`);
    }

    try {
      // 1. Assemble dynamic bordered pencil (steady = false, frozen at this.time) with cut-current rows.
      this.allTerminals = true;
      this.steady = false;
      const pencil = buildPencil(this);
      const { Nc, K, N } = pencil;

      // Coordinate weights D for Arnoldi metric: thermal units everywhere (terminal V scaled by 1/VT).
      const coordScale = new Float64Array(N).fill(1);
      for (let a = 0; a < K; a++) coordScale[Nc + a] = 1 / this.VT;
      const invScale = Float64Array.from(coordScale, (v) => 1 / v);

      // 2. Matrix-verified deflation basis (wᵀ J = 0 to round-off, c = Mᵀ w from assembled M).
      const deflation = buildDeflation(this, pencil, coordScale, warnings);

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

      // 3. Precondition: verify state is actually steady for the frozen dynamic equations.
      const steadyCheck = checkDynamicSteady(this, pencil, deflation, coordScale);
      if (!steadyCheck.ok) {
        return {
          verdict: 'undecided',
          reason: steadyCheck.reason,
          drift: steadyCheck.drift,
          modes: [],
          deflated: deflation.meta,
          shifts: { attempted: [], succeeded: [] },
          window: null,
          warnings,
          ...(debug ? { debug: debugData } : {}),
        };
      }

      // 4. Shift-invert Arnoldi over real positive shifts (+ optional complex frequency shifts).
      const realShifts = userShifts ? Array.from(userShifts) : defaultShifts(this, pencil);
      const startVec = seededVector(N);
      const attempted = [];
      const succeeded = [];
      const candidates = [];
      const mKrylov = Math.min(32, Math.max(4, N - deflation.totalDeflated));

      // Meerbergen–Roose trust-region multiplier around each shift σ: |λ − σ| ≤ TRUST · |σ|.
      const TRUST = 25;

      const xTmp = new Float64Array(N), mxTmp = new Float64Array(N), solTmp = new Float64Array(N);
      for (const sigma of realShifts) {
        attempted.push({ re: sigma, im: 0 });
        const fac = factorRealShift(pencil, sigma, true);
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
          if (!(muMag > 1e-10 / sigma)) continue;
          const [invRe, invIm] = cdiv(1, 0, r.muRe, r.muIm);
          const re = sigma - invRe, im = -invIm;
          if (!Number.isFinite(re) || !Number.isFinite(im)) continue;
          if (Math.hypot(re - sigma, im) > TRUST * sigma) continue;
          if (im < -1e-10 * Math.abs(re)) continue;
          const deltaLam = r.ritzEst / (muMag * muMag);
          candidates.push({ re, im, vr: r.vr, vi: r.vi, ritzEst: r.ritzEst, deltaLam, shiftMag: sigma });
        }
      }

      const xTmpR = new Float64Array(N), xTmpI = new Float64Array(N);
      const mxTmpR = new Float64Array(N), mxTmpI = new Float64Array(N);
      const solTmpR = new Float64Array(N), solTmpI = new Float64Array(N);
      let maxOmega = 0;
      for (const f of frequencies) {
        const w = 2 * Math.PI * f;
        const sr = 1e-3 * w;
        const sMag = Math.hypot(sr, w);
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
          if (!(muMag > 1e-10 / sMag)) continue;
          const [invRe, invIm] = cdiv(1, 0, r.muRe, r.muIm);
          const re = sr - invRe, im = w - invIm;
          if (!Number.isFinite(re) || !Number.isFinite(im)) continue;
          if (Math.hypot(re - sr, im - w) > TRUST * sMag) continue;
          if (im < -1e-10 * Math.abs(re)) continue;
          const deltaLam = r.ritzEst / (muMag * muMag);
          candidates.push({ re, im, vr: r.vr, vi: r.vi, ritzEst: r.ritzEst, deltaLam, shiftMag: sMag });
        }
      }

      if (succeeded.length === 0) {
        return {
          verdict: 'undecided',
          reason: 'no shifted linearisation factored without singularity',
          drift: steadyCheck.drift,
          modes: [],
          deflated: deflation.meta,
          shifts: { attempted, succeeded },
          window: null,
          warnings,
          ...(debug ? { debug: debugData } : {}),
        };
      }

      const rHi = Math.max(...succeeded.filter((x) => x.im === 0).map((x) => x.re), 0);
      const rSearch = Math.max(rHi, maxOmega);
      const shiftScale = Math.min(...succeeded.map((x) => Math.hypot(x.re, x.im)));

      // 5. Scale-invariant pre-dedup (never merging opposite-sign Re λ), then iterative refinement
      //    and biorthogonal post-deduplication.
      candidates.sort((a, b) => b.re - a.re || Math.abs(a.im) - Math.abs(b.im));
      const uniqueCands = [];
      for (const c of candidates) {
        const cMag = Math.hypot(c.re, c.im);
        const dup = uniqueCands.find((u) => {
          if ((u.re > 0) !== (c.re > 0) && u.re !== 0 && c.re !== 0) return false;
          const scale = Math.max(Math.hypot(u.re, u.im), cMag);
          return Math.hypot(u.re - c.re, u.im - c.im) <= 1e-4 * scale;
        });
        if (!dup) uniqueCands.push(c);
        else if (c.deltaLam < dup.deltaLam) {
          dup.re = c.re; dup.im = c.im; dup.vr = c.vr; dup.vi = c.vi;
          dup.ritzEst = c.ritzEst; dup.deltaLam = c.deltaLam; dup.shiftMag = c.shiftMag;
        }
      }

      const refined = [];
      let acceptedSlots = 0;
      const targetSlots = Math.max(nWant + 4, 12);
      const maxAttempts = Math.min(uniqueCands.length, Math.max(nWant * 12, 160));
      let rejectedSpurious = 0;
      for (let i = 0; i < maxAttempts; i++) {
        if (acceptedSlots >= targetSlots) {
          let worst = Infinity;
          for (const u of refined) if (u.re < worst) worst = u.re;
          if (uniqueCands[i].re < worst - Math.abs(worst) * 1e-2) break;
        }
        const res = refineAndBound(pencil, deflation, coordScale, uniqueCands[i], shiftScale);
        if (!res) {
          rejectedSpurious++;
          continue;
        }
        // Biorthogonal + bar-based post-dedup (never merges growing into decaying).
        const resMag = Math.hypot(res.re, res.im);
        const dup = refined.find((u) => {
          const uPos = u.re > kBar * u.error, uNeg = u.re < -kBar * u.error;
          const rPos = res.re > kBar * res.error, rNeg = res.re < -kBar * res.error;
          if ((uPos && rNeg) || (uNeg && rPos)) return false;
          const dist = Math.hypot(u.re - res.re, u.im - res.im);
          const scale = Math.max(Math.hypot(u.re, u.im), resMag);
          const barTol = 5 * (u.error + res.error) + 1e-9 * scale;
          if (dist <= barTol) return biorthOverlap(u, res) > 0.25;
          if (dist <= 1e-3 * scale) return biorthOverlap(u, res) > 0.5;
          return false;
        });
        if (!dup) {
          refined.push(res);
          acceptedSlots += Math.abs(res.im) > 1e-10 * Math.abs(res.re) ? 2 : 1;
        } else if (res.error < dup.error || (res.error === dup.error && res.residual < dup.residual)) {
          Object.assign(dup, res);
        }
      }

      // Expand complex conjugate pairs (re ± i im), sorted rightmost first.
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
        if (Math.abs(m.im) > 1e-10 * Math.abs(m.re)) {
          const vic = Float64Array.from(m.vi, (v) => -v);
          allModes.push({
            re: m.re,
            im: -m.im,
            residual: m.residual,
            error: m.error,
            ...(profiles ? { profiles: modeProfiles(this, pencil, m.vr, vic) } : {}),
          });
        }
      }
      const modesOut = allModes.slice(0, nWant);

      // 6. Compute honest window and verdict over ALL verified modes inside the searched window.
      //    No mode outside the searched window (Re λ > rSearch) may decide the verdict.
      const inWindow = allModes.filter((m) => m.re <= rSearch);
      let rLo;
      if (inWindow.length > 0) {
        // Resolution floor near the imaginary axis (rightmost / smallest-|Re λ| verified modes).
        let minAbsRe = Infinity;
        for (const m of inWindow) { const a = Math.abs(m.re); if (a < minAbsRe) minAbsRe = a; }
        const near = inWindow.filter((m) => Math.abs(m.re) <= 10 * minAbsRe);
        const pool = near.length > 0 ? near : [inWindow[0]];
        rLo = Math.max(...pool.map((m) => kBar * m.error));
      } else {
        rLo = EPS * Math.max(rSearch, shiftScale);
      }
      const window = { rLo, rHi, omegaHi: maxOmega };

      const growing = inWindow.find((m) => m.re > kBar * m.error);
      const neutral = inWindow.find((m) => Math.abs(m.re) <= kBar * m.error);
      let verdict, reason;
      if (growing) {
        verdict = 'unstable';
        reason = `growing mode at Re λ = ${growing.re.toPrecision(4)} ± ${growing.error.toPrecision(2)} /s` +
          (Math.abs(growing.im) > 0 ? ` (ω = ${Math.abs(growing.im).toPrecision(4)} rad/s)` : '');
      } else if (this.unkept.length > 0) {
        verdict = 'undecided';
        reason = `${this.unkept.length} conserved combination(s) are in unkept; undeflated neutral directions remain in the pencil`;
      } else if (rHi > 0 && rLo >= rHi) {
        verdict = 'undecided';
        reason = `resolution floor (${rLo.toPrecision(2)} /s) exceeds the searched window rHi (${rHi.toPrecision(2)} /s)`;
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
        drift: steadyCheck.drift,
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
      this.sourceBefore = wasSourceBefore;
      this.offNeutral = wasOffNeutral;
      this.lostMode = learnedKeep.lostMode;
      this.accepted = learnedKeep.accepted;
      this.firstErrs = learnedKeep.firstErrs;
      this._restore(snap);
    }
  }
}
