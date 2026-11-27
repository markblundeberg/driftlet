// Pseudo-arclength branch continuation with linear stability marking and bifurcation location.
//
// CONTRACT
//
// Method `branch(opts)` on Solver (installed on Solver.prototype in src/solver.js) and exposed
// as `Device.prototype.branch(opts)` in src/index.js.
//
//   device.branch({ terminal = 'right', to, step, maxSteps = 200, voltageBound = 20,
//                   budgetMs = Infinity, stability = true, solutions = false })
//     → { terminal, kind, points, events, end, warnings }
//
// Precondition:
//   - The device must hold a converged, un-regridded steady state (`atSteady && !regriddedSteady`),
//     else throws SolverError. Invalid options throw SolverError.
//   - If the device cannot be continued directly (`!_directSteady()`, i.e. unkept conservation
//     laws), or the chosen terminal is current-driven and fixed by Kirchhoff (`_bound(k)`) or only
//     charges a capacitance/conserved pool at `to !== p0`, returns cleanly with the start point and
//     `end: 'refused: ...'`.
//
// Continuation & step control:
//   - Traces the 1D steady solution manifold F(x, V_fl, p) = 0 in (state, p), where p is the
//     chosen terminal's own drive value (V for a held or R-backed voltage terminal, I for a
//     current-driven terminal), keeping all conserved amounts (lawRefs), floating charges
//     (chargeRefs) and immobile sums (localKept) fixed at their starting values.
//   - Dimensionless span coordinate s = (p − p0) / (to − p0) (0 at start, 1 at `to`), oriented
//     initially with ds/d(arc) > 0. The arclength metric combines bulk potentials/levels (RMS over
//     bulk nodes, mesh-independent), localized unknowns (gates weighted 10×, surfaces, floating
//     terminal voltages in thermal units) without bulk dilution, and the drive coordinate scaled
//     dynamically so the state carries at least half the tangent weight (`wsEff·|ts| ≤ 0.85` when
//     the state moves). All drive tolerances are relative to the local drive/bracket scale — never
//     absolute in drive units. Bordered Woodbury solves reuse s.sys (pure JS, zero dependencies,
//     no ** or Math.pow).
//   - Every output point is converged to `tol = 1e-10` thermal units. On devices with floating
//     islands (`this.islands.length > 0`), if the plain corrector leaves Kirchhoff unbalanced, a
//     second-phase islands round (`islandsOn = true`, fixed s) polishes the island levels (as
//     newton() does).
//   - Step acceptance guards (applied uniformly to ordinary and landing steps; never switched off
//     at the minimum stride):
//       1. Corrector contraction: no Newton damping (`alpha === 1`) and classical contraction
//          `δ_2 / δ_1 ≤ 0.5` (for `δ_1 > 1e-6`), plus bounded pointwise state motion.
//       2. Tangent continuity: `cos θ ≥ 0.75` (≤ 41°) in the balanced metric and forward alignment
//          `t_prev · Δz > 0`.
//       3. Hard `step` cap: when `opts.step` is given, `|p_next − p_cur| ≤ step` on every step.
//       4. All-mode nearer-end guard (when `stability !== false`): unique modes (im ≥ 0) are matched
//          1-to-1 by biorthogonal eigenvector overlap. On the same side of Re λ = 0, every matched
//          active mode must satisfy `|Re λ_B − Re λ_A| ≤ 0.75 · max(min(|Re λ_A|, |Re λ_B|), kBar · error)`.
//          Modes within `kBar · error` of zero are not resolved below their error bar. A genuine
//          stability verdict flip with no fold or Hopf (e.g. a real eigenvalue passing through
//          infinity) is audited as a `stability mismatch` event rather than stopping blind.
//   - When pseudo-arclength cannot advance at the minimum stride (corrector failure, singular
//     Jacobian, or model wall), continuation stops honestly and notes that `solve()` may continue
//     from the last accepted drive rather than taking an unguarded natural-parameter jump.
//
// One span rule:
//   - The trace stays inside the requested drive span `s ∈ [0, 1]` (between start drive `p0` and
//     `to`). It ends where it first leaves that span, landing cleanly on the crossed edge (`to` or
//     `p0`) and naming which edge in `end`. When a step turns a fold, the fold's drive is located
//     first: if the fold lies past `to` (sFold ≥ 1), the trace lands on `to` on the pre-fold leg
//     with no fold event; if a returning leg crosses back through `p0` (even on the same step as
//     the fold), it records the fold and lands on `p0` on the far branch.
//
// Stability marks & events (and what a mark does not promise):
//   - `points` always contains at least the starting point `points[0]`. If `budgetMs` runs out
//     during the start point's `stability()` call, `points[0]` is still kept with
//     `stability: { verdict: 'undecided', checked: false, ..., reason }`.
//   - With `stability !== false`, every accepted point is evaluated by a full `this.stability(...)`
//     call at that point's state and drive, agreeing with a standalone `dev.stability()` call at
//     the same state. Marks are for the drive as held (a curve can be stable under V-drive and
//     unstable under I-drive).
//   - What a mark does NOT promise: a stability mark holds at the sampled point. The nearer-end
//     guard zooms steps proportionally to every matched mode's distance from Re λ = 0, catching
//     narrow windows approached by listed modes, but modes within their error bar (`kBar · error`)
//     of zero are not resolved and an unstable window narrower than the step that is unseen by
//     listed modes at both endpoints can still be missed.
//   - Fold (`dp/ds` changes sign): located by safeguarded secant/bisection along the exact arc to
//     `|ts| ≤ 1e-8 · max(|tsA|, |tsB|)` (or relative chord `< 1e-9`), reporting
//     `{ kind: 'fold', drive, terminals, between }`.
//   - Hopf (`Re λ` of a complex conjugate pair `Im λ > 0` crosses 0): matched by eigenvector
//     overlap and located by safeguarded secant/bisection along the exact arc to
//     `|Re λ| ≤ max(kBar·error, 1e-5·ω)` or `|pB − pA| ≤ 1e-6 · max(|pA|, |pB|, |pB0 − pA0|)`,
//     reporting `{ kind: 'hopf', drive, terminals, omega, between }`.
//   - Fallbacks are always reported: if a fold or Hopf location fails to converge or is cut short
//     by `budgetMs`, the event is still emitted with `located: false` and a `note` saying why.
//   - Consistency audit: runs on both normal completion and `budgetMs` stops. At a fold in the
//     held parameter, one real eigenvalue must pass through 0. If stability marks do not change
//     across a fold (without walking across any neighbouring event), or flip between decided states
//     with no fold or Hopf located, a `{ kind: 'stability mismatch', ... }` event is emitted.
//
// Stop conditions (every stop is named in `end`):
//   - Reached `to` (lands cleanly at p = to).
//   - Returned to starting drive `p0` after a fold (lands cleanly at p = p0 on the far branch).
//   - Left `voltageBound` (any floating terminal moved more than `voltageBound` V from its start).
//   - `maxSteps` reached.
//   - `budgetMs` spent (unlimited `budgetMs = Infinity` never reads the wall clock; finite
//     `budgetMs` is checked inside every Newton iteration, shift factorisation, mode refinement,
//     and event bisection so overruns stay within roughly one Newton round or shift factorisation).
//   - Continuation step failed (stride underflow, singular bordered system, model wall, Kirchhoff).
//
// Invariants:
//   1. Completely non-destructive: state, drives, `sourceOverride`, `drivenKey`, laws, and LEARNED
//      counters are restored bit-for-bit in `finally` (even on throw or budget stop).
//   2. Deterministic when `budgetMs = Infinity` (never reads the wall clock).
//   3. Zero change to `solve()`, `advance()`, `impedance()`, `stability()` (when called directly),
//      or `homotopy.js`.

import { sourceAt } from './device.js';
import { SolverError } from './errors.js';
import { carrierZ, chargeOnly, feeds } from './laws.js';
import { captureLaws, keepLocal } from './lawRows.js';
import { addToState, clipLevels, equilibrate, maxPotentialStep } from './newton.js';
import { makeSolution } from './solution.js';

const GATE_WEIGHT = 10; // weight on [0, 1] gate fractions (Δgate = 0.1 ↔ 1 thermal unit)
const BRANCH_TOL = 1e-10; // corrector Newton tolerance in thermal units (~2.5e-12 V)
const MIN_ARC = 1e-6;
const MAX_ARC = 1e5;
const COS_MIN = 0.75; // max tangent turn ≈ 41° per step
const KAPPA_MAX = 0.5; // classical Newton contraction factor bound δ_{k} / δ_{k-1}
const ALPHA_RE = 0.75; // nearer-end mode rule: |ΔRe λ| ≤ ALPHA_RE · max(min(|a|, |b|), kBar · err)
const BUDGET_SPENT = Symbol('branch budget spent');
const now = () => (globalThis.performance ? globalThis.performance.now() : Date.now());

function branchBudget(ms) {
  const start = now();
  let attempted = false;
  return {
    check() {
      if (attempted && now() - start >= ms) throw BUDGET_SPENT;
    },
    checkWork() {
      if (attempted && now() - start >= ms) throw BUDGET_SPENT;
      attempted = true;
    },
    markAttempted() {
      attempted = true;
    },
    record() {
      attempted = true;
    },
    remaining() {
      return Math.max(0, ms - (now() - start));
    },
    elapsed() {
      return now() - start;
    },
  };
}

// Balanced, mesh-independent squared state weights: bulk node unknowns are averaged (1 / nBulk),
// while localized 0D unknowns (gates, surfaces, diffusion face vars) are NOT diluted by nBulk.
function stateMetricEntries(s) {
  const { n, M, rix: R } = s;
  const bulkIdx = [];
  const surfIdx = [];
  const gateIdx = [];
  const diffIdx = [];
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (s.nodeConductor[g] < 0 && s.loc[b * M] >= 0) bulkIdx.push(R[b * M]);
    for (let i = 0; i < n; i++) {
      if (s.present[g * n + i] && s.loc[b * M + 1 + i] >= 0) bulkIdx.push(R[b * M + 1 + i]);
    }
    const k = s.surfAt[g];
    if (k >= 0) {
      for (let q = 0; q < s.model.surfaces[k].species.length; q++) {
        if (s.loc[b * M + 1 + n + q] >= 0) surfIdx.push(R[b * M + 1 + n + q]);
      }
    }
  }
  for (const gt of s.gateList) {
    if (s.loc[gt.o] >= 0) gateIdx.push(R[gt.o]);
  }
  for (const d of s.diffList) {
    const bf = s.blockOfFace[d.f];
    if (s.loc[d.o] >= 0) diffIdx.push(R[d.o]);
    for (const q of d.drawSlot) if (q >= 0 && s.loc[bf * M + q] >= 0) diffIdx.push(R[bf * M + q]);
  }
  const total = bulkIdx.length + surfIdx.length + gateIdx.length + diffIdx.length;
  const potIdx = new Int32Array(total);
  const potWt2 = new Float64Array(total);
  let p = 0;
  const wBulk = bulkIdx.length > 0 ? 1 / bulkIdx.length : 1;
  for (let i = 0; i < bulkIdx.length; i++) { potIdx[p] = bulkIdx[i]; potWt2[p++] = wBulk; }
  const wSurf = surfIdx.length > 0 ? 1 / surfIdx.length : 1;
  for (let i = 0; i < surfIdx.length; i++) { potIdx[p] = surfIdx[i]; potWt2[p++] = wSurf; }
  const wGate = GATE_WEIGHT * GATE_WEIGHT;
  for (let i = 0; i < gateIdx.length; i++) { potIdx[p] = gateIdx[i]; potWt2[p++] = wGate; }
  const wDiff = diffIdx.length > 0 ? 1 / diffIdx.length : 1;
  for (let i = 0; i < diffIdx.length; i++) { potIdx[p] = diffIdx[i]; potWt2[p++] = wDiff; }
  return { potIdx, potWt2 };
}

function maxStateStep(s, dx, dV) {
  let mx = maxPotentialStep(s, dx, dV);
  const R = s.rix;
  for (const gt of s.gateList) {
    const v = GATE_WEIGHT * Math.abs(dx[R[gt.o]]);
    if (v > mx) mx = v;
  }
  return mx;
}

function maxStateDiff(s, snapA, uB, uLoB, termVB) {
  const { n, M, floating: fl, VT } = s;
  let mx = 0;
  for (let a = 0; a < fl.length; a++) {
    const k = fl[a];
    mx = Math.max(mx, Math.abs(termVB[k] - snapA.termV[k]) / VT);
  }
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (s.nodeConductor[g] < 0) {
      const o = b * M;
      mx = Math.max(mx, Math.abs((uB[o] - snapA.u[o]) + (uLoB[o] - snapA.uLo[o])));
    }
    for (let i = 0; i < n; i++) {
      if (s.present[g * n + i]) {
        const o = b * M + 1 + i;
        mx = Math.max(mx, Math.abs((uB[o] - snapA.u[o]) + (uLoB[o] - snapA.uLo[o])));
      }
    }
    const k = s.surfAt[g];
    if (k >= 0) {
      for (let q = 0; q < s.model.surfaces[k].species.length; q++) {
        const o = b * M + 1 + n + q;
        mx = Math.max(mx, Math.abs((uB[o] - snapA.u[o]) + (uLoB[o] - snapA.uLo[o])));
      }
    }
  }
  for (const gt of s.gateList) {
    const o = gt.o;
    mx = Math.max(mx, GATE_WEIGHT * Math.abs((uB[o] - snapA.u[o]) + (uLoB[o] - snapA.uLo[o])));
  }
  for (const d of s.diffList) {
    const bf = s.blockOfFace[d.f];
    mx = Math.max(mx, Math.abs((uB[d.o] - snapA.u[d.o]) + (uLoB[d.o] - snapA.uLo[d.o])));
    for (const q of d.drawSlot) {
      if (q >= 0) {
        const o = bf * M + q;
        mx = Math.max(mx, Math.abs((uB[o] - snapA.u[o]) + (uLoB[o] - snapA.uLo[o])));
      }
    }
  }
  return mx;
}

function stateOnlyNorm(potIdx, potWt2, VT, tx, tV) {
  let sum = 0;
  for (let i = 0; i < potIdx.length; i++) {
    const v = tx[potIdx[i]];
    sum += potWt2[i] * v * v;
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tV.length; a++) sum += invVT2 * tV[a] * tV[a];
  return Math.sqrt(sum);
}

function tangentNorm(potIdx, potWt2, wsEff, VT, tx, tV, ts) {
  let sum = wsEff * wsEff * ts * ts;
  for (let i = 0; i < potIdx.length; i++) {
    const v = tx[potIdx[i]];
    sum += potWt2[i] * v * v;
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tV.length; a++) sum += invVT2 * tV[a] * tV[a];
  return Math.sqrt(sum);
}

function buildWeights(potIdx, potWt2, wsEff, VT, tx, tV, ts, wx, wV) {
  wx.fill(0);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i];
    wx[j] = potWt2[i] * tx[j];
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tV.length; a++) wV[a] = invVT2 * tV[a];
  return wsEff * wsEff * ts;
}

function stateDistW(s, ctx, uA, uLoA, termVA, sA, uB, uLoB, termVB, sB) {
  const { potIdx, potWt2, wsEff } = ctx;
  const fl = s.floating, VT = s.VT;
  let sum = wsEff * wsEff * (sB - sA) * (sB - sA);
  for (let i = 0; i < potIdx.length; i++) {
    const o = s.fullOf[potIdx[i]];
    const d = (uB[o] - uA[o]) + (uLoB[o] - uLoA[o]);
    sum += potWt2[i] * d * d;
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < fl.length; a++) {
    const k = fl[a];
    const d = termVB[k] - termVA[k];
    sum += invVT2 * d * d;
  }
  return Math.sqrt(sum);
}

function stepAlignment(s, ctx, refTx, refTV, refTs, snapA, sA, uB, uLoB, termVB, sB) {
  const { potIdx, potWt2, wsEff } = ctx;
  const fl = s.floating, invVT2 = 1 / (s.VT * s.VT);
  let align = wsEff * wsEff * refTs * (sB - sA);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i], o = s.fullOf[j];
    align += potWt2[i] * refTx[j] * ((uB[o] - snapA.u[o]) + (uLoB[o] - snapA.uLo[o]));
  }
  for (let a = 0; a < fl.length; a++) {
    const k = fl[a];
    align += invVT2 * refTV[a] * (termVB[k] - snapA.termV[k]);
  }
  return align;
}

// Exact ∂F/∂s for single-terminal continuation p(s) = p0 + s * dpSpan after assemble(Infinity)
// and equilibrate(s).
function computeParamDs(s, termIdx, dpSpan, Bs, pinDs, circDs) {
  const N = s.sys.size, fl = s.floating, pins = s.pins;
  Bs.fill(0);
  pinDs.fill(0);
  circDs.fill(0);
  const d = s.terms[termIdx].drive;
  if (d.kind === 'V' && !(d.R > 0)) {
    const Bk = s.termB[termIdx];
    for (let j = 0; j < N; j++) if (Bk[j] !== 0) Bs[j] = dpSpan * Bk[j];
    for (let q = 0; q < pins.length; q++) if (pins[q].dV) pinDs[q] = dpSpan * pins[q].dV[termIdx];
    for (let a = 0; a < fl.length; a++) circDs[a] = dpSpan * s.termDX[fl[a]][termIdx];
  } else {
    const a = fl.indexOf(termIdx);
    if (a >= 0) circDs[a] = d.kind === 'I' ? -dpSpan : -dpSpan / d.R;
  }
}

// Solve the (P + K + 1)-bordered system reusing the factored s.sys.
function solveBordered(s, work, rhs, pinRhs, circRhs, gRhs, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws) {
  const N = s.sys.size, fl = s.floating, pins = s.pins, P = pins.length, K = fl.length, S = P + K + 1;
  s._solveLinear(rhs, delta);
  const e = work.e, cols = work.cols;
  for (let a = 0; a < S; a++) {
    if (!cols[a]) cols[a] = new Float64Array(N + 1);
    if (a < P) {
      const col = pins[a].col ?? [pins[a].row];
      for (const o of col) e[s.rix[o]] = 1;
      try {
        s._solveLinear(e, cols[a]);
      } finally {
        for (const o of col) e[s.rix[o]] = 0;
      }
    } else if (a < P + K) {
      s._solveLinear(s.termB[fl[a - P]], cols[a]);
    } else {
      s._solveLinear(Bs, cols[a]);
    }
  }
  const dot = (a, v) => {
    let t = 0;
    if (a < P) {
      const p = pins[a];
      for (let j = 0; j < p.len; j++) t += p.w[j] * v[p.idx[j]];
    } else if (a < P + K) {
      const C = s.termC[fl[a - P]];
      for (let j = 0; j < N; j++) t += C[j] * v[j];
    } else {
      for (let j = 0; j < N; j++) if (wx[j] !== 0) t += wx[j] * v[j];
    }
    return t;
  };
  const A = Array.from({ length: S }, () => new Float64Array(S + 1));
  for (let a = 0; a < S; a++) {
    for (let b = 0; b < S; b++) A[a][b] = (b < P ? 1 : -1) * dot(a, cols[b]);
    if (a < P) {
      if (pins[a].termSlope) for (let b = P; b < P + K; b++) A[a][b] += pins[a].termSlope.get(fl[b - P]) ?? 0;
      A[a][P + K] += pinDs[a];
      A[a][S] = pinRhs[a] - dot(a, delta);
    } else if (a < P + K) {
      const fk = fl[a - P];
      A[a][a] += s.termDI[fk];
      for (let b = P; b < P + K; b++) if (b !== a) A[a][b] += s.termDX[fk][fl[b - P]];
      A[a][P + K] += circDs[a - P];
      A[a][S] = circRhs[a - P] - dot(a, delta);
    } else {
      for (let b = 0; b < K; b++) A[a][P + b] += wV[b];
      A[a][P + K] += ws;
      A[a][S] = gRhs - dot(a, delta);
    }
    let mx = 0;
    for (let b = 0; b < S; b++) mx = Math.max(mx, Math.abs(A[a][b]));
    if (mx > 0) for (let b = 0; b <= S; b++) A[a][b] /= mx;
  }
  for (let col = 0; col < S; col++) {
    let p = col;
    for (let r = col + 1; r < S; r++) if (Math.abs(A[r][col]) > Math.abs(A[p][col])) p = r;
    if (p !== col) [A[col], A[p]] = [A[p], A[col]];
    const piv = A[col][col];
    if (!(piv !== 0 && Number.isFinite(piv))) throw new SolverError('the bordered arclength system is singular');
    for (let r = col + 1; r < S; r++) {
      const f = A[r][col] / piv;
      for (let k = col; k <= S; k++) A[r][k] -= f * A[col][k];
    }
  }
  const x = new Float64Array(S);
  for (let r = S - 1; r >= 0; r--) {
    let v = A[r][S];
    for (let k = r + 1; k < S; k++) v -= A[r][k] * x[k];
    x[r] = v / A[r][r];
  }
  for (let a = 0; a < S; a++) {
    const sa = a < P ? x[a] : -x[a];
    const col = cols[a];
    for (let j = 0; j < N; j++) delta[j] += sa * col[j];
    if (a >= P && a < P + K) deltaV[a - P] = x[a];
  }
  return { ds: x[P + K] };
}

// Compute unit tangent (tx, tV, ts) at current converged state + sVal, oriented to agree with
// (prevTx, prevTV, prevTs), or with ts > 0 when prevTx is null.
function computeTangent(s, ctx, sVal, prevTx, prevTV, prevTs, tx, tV, reuseFactor = false) {
  const { termIdx, p0, dpSpan, potIdx, potWt2, wsEff, Bs, pinDs, circDs, wx, wV, work } = ctx;
  const N = s.sys.size, K = s.floating.length;
  if (!reuseFactor) {
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    s.islandsOn = false;
    s.assemble(Infinity);
    equilibrate(s);
    computeParamDs(s, termIdx, dpSpan, Bs, pinDs, circDs);
    s._factor();
  }
  let ws;
  if (prevTx) {
    ws = buildWeights(potIdx, potWt2, wsEff, s.VT, prevTx, prevTV, prevTs, wx, wV);
  } else {
    wx.fill(0);
    wV.fill(0);
    ws = 1;
  }
  const res = solveBordered(s, work, work.zero, work.zeroPin, work.zeroCirc, 1, tx, tV, Bs, pinDs, circDs, wx, wV, ws);
  const tsRaw = res.ds;
  if (!Number.isFinite(tsRaw)) return { ok: false };
  const nrm = tangentNorm(potIdx, potWt2, wsEff, s.VT, tx, tV, tsRaw);
  if (!(nrm > 0 && Number.isFinite(nrm))) return { ok: false };
  for (let j = 0; j < N; j++) tx[j] /= nrm;
  for (let a = 0; a < K; a++) tV[a] /= nrm;
  let ts = tsRaw / nrm;
  let cosTheta = 1;
  if (!prevTx) {
    if (ts < 0) {
      for (let j = 0; j < N; j++) tx[j] = -tx[j];
      for (let a = 0; a < K; a++) tV[a] = -tV[a];
      ts = -ts;
    }
  } else {
    cosTheta = 1 / nrm;
  }
  return { ok: true, ts, cosTheta };
}

function kirchhoffSum(s) {
  let sum = 0;
  for (const I of s.termI) sum += I;
  return Math.abs(sum);
}

// Corrector Newton on F(x, V_fl, s) = 0 and wx·(x − xPred) + wV·(V − vPred) + ws·(s − sPred) = 0.
// Tracks classical Newton contraction factor `kappa = max(δ_k / δ_{k-1})` and damping flag.
function correctorNewton(s, ctx, sStart, uPred, uLoPred, vPred, sPred, wx, wV, ws, { maxIter = 40, tol = BRANCH_TOL, maxStep = 10, fixedS = false } = {}) {
  const { termIdx, p0, dpSpan, potIdx, potWt2, wsEff, Bs, pinDs, circDs, work, budget } = ctx;
  budget?.check();
  const fl = s.floating, K = fl.length;
  const delta = work.delta, deltaV = work.deltaV;
  const pinRhs = work.pinRhs, circRhs = work.circRhs;
  let sVal = fixedS ? sPred : sStart;

  const runLoop = (islandsOn, isFixed, limit) => {
    let d1 = 0;
    let d2 = 0;
    let damped = false;
    for (let it = 1; it <= limit; it++) {
      budget?.check();
      s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
      s.islandsOn = islandsOn;
      try {
        s.assemble(Infinity);
      } catch (err) {
        if (!(err instanceof SolverError)) throw err;
        return { converged: false, iterations: it, s: sVal, modelError: err.message };
      }
      equilibrate(s);
      computeParamDs(s, termIdx, dpSpan, Bs, pinDs, circDs);
      try {
        s._factor();
      } catch (err) {
        if (!(err instanceof Error)) throw err;
        return { converged: false, iterations: it, s: sVal, error: err.message };
      }
      let g = isFixed ? 0 : ws * (sVal - sPred);
      if (!isFixed) {
        for (let i = 0; i < potIdx.length; i++) {
          const j = potIdx[i], o = s.fullOf[j];
          g += wx[j] * ((s.u[o] - uPred[o]) + (s.uLo[o] - uLoPred[o]));
        }
        for (let a = 0; a < K; a++) g += wV[a] * (s.termV[fl[a]] - vPred[a]);
      }
      const pins = s.pins;
      for (let q = 0; q < pins.length; q++) pinRhs[q] = pins[q].res;
      for (let a = 0; a < K; a++) circRhs[a] = s.termRes[fl[a]];
      let solved;
      try {
        deltaV.fill(0);
        if (isFixed) {
          work.zeroWx.fill(0);
          work.zeroWv.fill(0);
          solved = solveBordered(s, work, s.res, pinRhs, circRhs, 0, delta, deltaV, Bs, pinDs, circDs, work.zeroWx, work.zeroWv, 1);
        } else {
          solved = solveBordered(s, work, s.res, pinRhs, circRhs, g, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws);
        }
      } catch (err) {
        if (!(err instanceof SolverError)) throw err;
        return { converged: false, iterations: it, s: sVal, error: err.message };
      }
      const ds = isFixed ? 0 : solved.ds;
      clipLevels(s, delta, maxStep);
      const step = Math.max(maxPotentialStep(s, delta, deltaV), wsEff * Math.abs(ds));
      const own = Math.max(maxPotentialStep(s, delta, deltaV, false), wsEff * Math.abs(ds));
      if (!Number.isFinite(step)) return { converged: false, iterations: it, s: sVal, error: 'non-finite update' };
      if (it > 1 && own > 1e5) return { converged: false, iterations: it, s: sVal, error: 'diverging' };
      const dNorm = tangentNorm(potIdx, potWt2, wsEff, s.VT, delta, deltaV, ds);
      if (it === 1) d1 = dNorm;
      else if (it === 2) d2 = dNorm;
      const alpha = own > maxStep ? maxStep / own : 1;
      if (alpha < 1) damped = true;
      addToState(s, delta, -alpha);
      for (let a = 0; a < K; a++) s.termV[fl[a]] -= alpha * deltaV[a];
      if (!isFixed) sVal -= alpha * ds;
      budget?.markAttempted();
      if (alpha === 1 && step < tol) {
        try {
          s.computeConcentrations();
        } catch (err) {
          if (!(err instanceof SolverError)) throw err;
          return { converged: false, iterations: it, s: sVal, modelError: err.message };
        }
        return {
          converged: true,
          iterations: it,
          s: sVal,
          d1,
          kappa: d1 > 1e-6 ? d2 / d1 : 0,
          damped,
          reusable: !islandsOn && s.sys.perturbedPivots === 0 && s.sys.lostColumns === 0,
        };
      }
    }
    return { converged: false, iterations: limit, s: sVal, error: 'Newton did not converge' };
  };

  const r1 = runLoop(false, fixedS, maxIter);
  if (!r1.converged) return r1;

  // Phase 2: on floating-island devices, if plain convergence left Kirchhoff unbalanced, polish
  // islands at fixed s (mirroring newton()'s islands round).
  if (s.islands.length > 0 && s._missesKirchhoff()) {
    const preIslSnap = s._snapshot();
    for (const isl of s.islands) isl.cutAt = null;
    const r2 = runLoop(true, true, 25);
    s.islandsOn = false;
    const okKirchhoff = () => !s._missesKirchhoff() || (ctx.startMissesKirchhoff && kirchhoffSum(s) <= 10 * ctx.startKirchhoffSum + 1e-13);
    if (r2.converged && okKirchhoff()) {
      return { converged: true, iterations: Math.max(r1.iterations, 6), s: sVal, d1: r1.d1, kappa: r1.kappa, damped: r1.damped, reusable: false };
    }
    s._restore(preIslSnap);
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    const qPol = s._steadySteps(budget ? { maxSteps: 1, budget } : { maxSteps: 1 });
    if (qPol.converged && okKirchhoff()) {
      return { converged: true, iterations: Math.max(r1.iterations, 6), s: sVal, d1: r1.d1, kappa: r1.kappa, damped: r1.damped, reusable: false };
    }
    s._restore(preIslSnap);
    return { converged: false, iterations: r1.iterations + r2.iterations, s: sVal, error: 'floating-island current balance lost' };
  }
  return r1;
}

function chargingReason(s, k, I) {
  if (s._bound(k)) {
    return `terminal '${s.terms[k].name}' is bound by Kirchhoff (all other terminals are current-driven)`;
  }
  if (chargeOnly(s, k)) {
    return `terminal '${s.terms[k].name}' only charges a capacitance (no steady state at I ≠ 0)`;
  }
  const net = (cs) => {
    let sum = 0, size = 0;
    for (const t of cs.terminals) {
      if (s.terms[t].drive.kind !== 'I') continue;
      const it = t === k ? I : sourceAt(s.terms[t].drive.src, s.time);
      const cr = cs.carriers?.find((q) => q.k === t);
      if (cr) {
        const x = (cr.w / carrierZ(s, t)) * it;
        sum += x;
        size += Math.abs(x);
        continue;
      }
      const p = cs.parts.find((q) => feeds(t, s.stretches[q.stretch]) && s.z[s.stretches[q.stretch].species] !== 0);
      if (!p) return true;
      const x = (p.w / s.z[s.stretches[p.stretch].species]) * it;
      sum += x;
      size += Math.abs(x);
    }
    return Math.abs(sum) > 1e-12 * size;
  };
  const laws = s.constraints.filter((cs) => cs.terminals.includes(k) && net(cs));
  if (laws.length > 0) {
    return `terminal '${s.terms[k].name}' feeds a conserved pool with no exit (no steady state at I ≠ 0)`;
  }
  return null;
}

function readPointData(s, termIdx, p, wantSolution) {
  const term = s.terms[termIdx];
  const origSrc = term.drive.src;
  const hadOverride = s.sourceOverride.has(termIdx);
  const oldOverride = hadOverride ? s.sourceOverride.get(termIdx) : undefined;
  const wasSteady = s.steady, wasAtSteady = s.atSteady, wasLastDt = s.lastDt, wasReadV = s.readV;
  const wasRestarted = s.restarted;
  term.drive.src = { value: p };
  s.sourceOverride.delete(termIdx);
  s.steady = false;
  s.atSteady = true;
  s.lastDt = Infinity;
  s.readV = undefined;
  try {
    s._refreshSources();
    s.computeConcentrations();
    s._assembleBookkeeping(Infinity);
    if (wantSolution) {
      const sol = makeSolution(s, { converged: true });
      const terminals = {};
      for (const [name, rec] of Object.entries(sol.terminals)) terminals[name] = { V: rec.V, current: rec.current };
      return { terminals, solution: sol };
    }
    const terminals = {};
    for (let k = 0; k < s.terms.length; k++) {
      terminals[s.terms[k].name] = { V: s.termV[k], current: s._currentInto(k) };
    }
    return { terminals };
  } finally {
    term.drive.src = origSrc;
    if (hadOverride) s.sourceOverride.set(termIdx, oldOverride);
    s.steady = wasSteady;
    s.atSteady = wasAtSteady;
    s.lastDt = wasLastDt;
    s.readV = wasReadV;
    s.restarted = wasRestarted;
    s._refreshSources();
  }
}

function evaluateStability(s, termIdx, p, stabOpts, isCurrentDrive, budget) {
  const term = s.terms[termIdx];
  const origSrc = term.drive.src;
  const hadOverride = s.sourceOverride.has(termIdx);
  const oldOverride = hadOverride ? s.sourceOverride.get(termIdx) : undefined;
  const wasSteady = s.steady, wasAtSteady = s.atSteady, wasRegrid = s.regriddedSteady;
  const wasDrivenKey = s.drivenKey;
  const wasLawRefs = s.lawRefs, wasChargeRefs = s.chargeRefs, wasLocalKept = s.localKept;
  term.drive.src = { value: p };
  s.sourceOverride.delete(termIdx);
  s.steady = false;
  s.atSteady = true;
  s.regriddedSteady = false;
  let relawed = false;
  const callOpts = budget
    ? { ...stabOpts, _rawModes: true, _budgetCheck: () => budget.checkWork() }
    : { ...stabOpts, _rawModes: true };
  try {
    s._refreshSources();
    if (isCurrentDrive && s._drivenKey() !== wasDrivenKey) {
      s._lawsNow();
      relawed = true;
    }
    return s.stability(callOpts);
  } finally {
    term.drive.src = origSrc;
    if (hadOverride) s.sourceOverride.set(termIdx, oldOverride);
    if (relawed) {
      s._lawsNow();
      s.lawRefs = wasLawRefs;
      s.chargeRefs = wasChargeRefs;
      s.localKept = wasLocalKept;
    }
    s.steady = wasSteady;
    s.atSteady = wasAtSteady;
    s.regriddedSteady = wasRegrid;
    s._refreshSources();
  }
}

function summarizeStability(st) {
  const m0 = st.modes.length > 0 ? st.modes[0] : null;
  return {
    verdict: st.verdict,
    rate: m0 ? m0.re : NaN,
    im: m0 ? Math.abs(m0.im) : NaN,
    error: m0 ? m0.error : NaN,
    reason: st.reason,
  };
}

function vecOverlap(a, b) {
  const ar = a._vr, ai = a._vi, br = b._vr, bi = b._vi;
  if (!ar || !br || ar.length !== br.length) return -1;
  let dr = 0, di = 0;
  for (let j = 0; j < ar.length; j++) {
    dr += ar[j] * br[j] + ai[j] * bi[j];
    di += ar[j] * bi[j] - ai[j] * br[j];
  }
  return Math.hypot(dr, di);
}

function uniqueModes(st) {
  if (!st || !st.modes) return [];
  const out = [];
  for (const m of st.modes) if (m.im >= 0) out.push(m);
  return out;
}

// Biorthogonal mode overlap using dual vectors u = Mᴴ w / |wᴴ M v|, falling back to right-vector
// cosine when duals are absent.
function biorthOverlap(a, b) {
  const ar = a._vr, ai = a._vi, br = b._vr, bi = b._vi;
  if (!ar || !br || ar.length !== br.length) return -1;
  const aur = a._ur, aui = a._ui, bur = b._ur, bui = b._ui;
  if (aur && bur) {
    let dr1 = 0, di1 = 0, dr2 = 0, di2 = 0;
    for (let j = 0; j < ar.length; j++) {
      dr1 += aur[j] * br[j] + aui[j] * bi[j];
      di1 += aur[j] * bi[j] - aui[j] * br[j];
      dr2 += bur[j] * ar[j] + bui[j] * ai[j];
      di2 += bur[j] * ai[j] - bui[j] * ar[j];
    }
    const ov = Math.sqrt(Math.hypot(dr1, di1) * Math.hypot(dr2, di2));
    if (Number.isFinite(ov)) return Math.min(1, ov);
  }
  return vecOverlap(a, b);
}

// 1-to-1 mode matching by biorthogonal eigenvector overlap.
function matchModes(stA, stB) {
  const A = uniqueModes(stA);
  const B = uniqueModes(stB);
  const nA = A.length, nB = B.length;
  const cands = [];
  for (let i = 0; i < nA; i++) {
    for (let j = 0; j < nB; j++) {
      let s = biorthOverlap(A[i], B[j]);
      if (s < 0) {
        const sc = Math.max(Math.hypot(A[i].re, A[i].im), Math.hypot(B[j].re, B[j].im), 1e-12);
        s = 1 / (1 + Math.hypot(A[i].re - B[j].re, A[i].im - B[j].im) / sc);
      }
      cands.push({ i, j, s });
    }
  }
  cands.sort((x, y) => y.s - x.s);
  const matchedA = new Int32Array(nA).fill(-1);
  const matchedB = new Int32Array(nB).fill(-1);
  const pairs = [];
  for (const c of cands) {
    if (c.s < 0.45) break;
    if (matchedA[c.i] >= 0 || matchedB[c.j] >= 0) continue;
    matchedA[c.i] = c.j;
    matchedB[c.j] = c.i;
    pairs.push({ i: c.i, j: c.j, a: A[c.i], b: B[c.j], overlap: c.s });
  }
  return { A, B, pairs };
}

// Single principled mode-tracking guard across all matched active modes:
//   |Re λ_B − Re λ_A| ≤ ALPHA_RE · max(min(|Re λ_A|, |Re λ_B|), kBar · err)
// on the same side of Re λ = 0.
function modeStepSmooth(stPrev, stNew, kBar, nWant, stats) {
  if (!stPrev || !stNew) return { ok: true };
  const { A, B, pairs } = matchModes(stPrev, stNew);
  if (A.length === 0 || B.length === 0) return { ok: true };
  const isOsc = (m) => m.im > 0 && (m.re >= 0 || m.im > 0.25 * Math.abs(m.re));
  if (stPrev.verdict === 'undecided' || stNew.verdict === 'undecided') {
    if (!A.some((m) => m.re >= 0 || isOsc(m)) && !B.some((m) => m.re >= 0 || isOsc(m))) return { ok: true };
  }
  const cutA = stPrev.modes.length >= nWant && A.length > 3 ? A[A.length - 1].re : -Infinity;
  const cutB = stNew.modes.length >= nWant && B.length > 3 ? B[B.length - 1].re : -Infinity;
  const isActive = (arr, st, cut, idx) => {
    const m = arr[idx];
    if (idx === 0 || m.re >= 0 || arr.length <= 2) return true;
    if (isOsc(m)) return !(st.modes.length >= nWant && idx === arr.length - 1);
    return !(st.modes.length >= nWant && (idx >= arr.length - 2 || m.re <= 0.7 * cut));
  };
  for (const { i, j, a, b } of pairs) {
    if (!isActive(A, stPrev, cutA, i) && !isActive(B, stNew, cutB, j)) continue;
    if ((a.re < 0) !== (b.re < 0)) continue;
    const err = Math.max(a.error, b.error);
    const dNear = Math.max(Math.min(Math.abs(a.re), Math.abs(b.re)), kBar * err);
    if (Math.abs(b.re - a.re) > ALPHA_RE * dNear) {
      if (stats) stats.modeNearEnd++;
      return { ok: false, reason: 'mode Re λ moved past nearer-end scale' };
    }
  }
  return { ok: true };
}

// Chord-hyperplane solve at fraction theta in (0, 1) between two converged states A and B.
function solveOnChord(s, ctx, A, B, theta) {
  const { potIdx, potWt2, wsEff, wx, wV, uPred, uLoPred, vPred } = ctx;
  const N = s.sys.size, fl = s.floating, K = fl.length;
  const chX = ctx.work.chordX, chV = ctx.work.chordV;
  chX.fill(0);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i], o = s.fullOf[j];
    chX[j] = (B.snap.u[o] - A.snap.u[o]) + (B.snap.uLo[o] - A.snap.uLo[o]);
  }
  for (let a = 0; a < K; a++) chV[a] = B.snap.termV[fl[a]] - A.snap.termV[fl[a]];
  const chS = B.s - A.s;
  const nrm = tangentNorm(potIdx, potWt2, wsEff, s.VT, chX, chV, chS);
  if (!(nrm > 0 && Number.isFinite(nrm))) return { ok: false, chordLen: 0 };
  for (let j = 0; j < N; j++) chX[j] /= nrm;
  for (let a = 0; a < K; a++) chV[a] /= nrm;
  const chSu = chS / nrm;

  s._restore(A.snap);
  for (let o = 0; o < s.u.length; o++) {
    uPred[o] = A.snap.u[o] + theta * (B.snap.u[o] - A.snap.u[o]);
    uLoPred[o] = A.snap.uLo[o] + theta * (B.snap.uLo[o] - A.snap.uLo[o]);
  }
  s.u.set(uPred);
  s.uLo.set(uLoPred);
  for (let a = 0; a < K; a++) {
    const v = A.snap.termV[fl[a]] + theta * (B.snap.termV[fl[a]] - A.snap.termV[fl[a]]);
    vPred[a] = v;
    s.termV[fl[a]] = v;
  }
  const sPred = A.s + theta * (B.s - A.s);
  const ws = buildWeights(potIdx, potWt2, wsEff, s.VT, chX, chV, chSu, wx, wV);
  const cor = correctorNewton(s, ctx, sPred, uPred, uLoPred, vPred, sPred, wx, wV, ws);
  if (!cor.converged) return { ok: false, chordLen: nrm };
  const tan = computeTangent(s, ctx, cor.s, A.tx, A.tV, A.ts, ctx.work.locTx, ctx.work.locTV, !!cor.reusable);
  if (!tan.ok) return { ok: false, chordLen: nrm };
  return {
    ok: true,
    chordLen: nrm,
    state: {
      snap: s._snapshot(),
      s: cor.s,
      tx: Float64Array.from(ctx.work.locTx),
      tV: Float64Array.from(ctx.work.locTV),
      ts: tan.ts,
    },
  };
}

// Locate a fold (ts = 0) between A0 and B0 (tsA * tsB < 0) by safeguarded secant/bisection,
// using purely relative tolerances scaled to the bracket's own tangent and chord length.
function locateFold(s, ctx, A0, B0) {
  let A = A0, B = B0;
  let best = Math.abs(A.ts) <= Math.abs(B.ts) ? A : B;
  const tsScale = Math.max(Math.abs(A0.ts), Math.abs(B0.ts));
  const chordLen0 = stateDistW(s, ctx, A0.snap.u, A0.snap.uLo, A0.snap.termV, A0.s, B0.snap.u, B0.snap.uLo, B0.snap.termV, B0.s);
  let converged = false;
  let solvedAny = false;
  let failNote = null;
  for (let it = 0; it < 32; it++) {
    ctx.budget?.check();
    const denom = Math.abs(A.ts) + Math.abs(B.ts);
    let theta = it % 3 === 2 || !(denom > 0) ? 0.5 : Math.abs(A.ts) / denom;
    if (theta < 0.15) theta = 0.15;
    if (theta > 0.85) theta = 0.85;
    let m = solveOnChord(s, ctx, A, B, theta);
    if (!m.ok && theta !== 0.5) m = solveOnChord(s, ctx, A, B, 0.5);
    if (!m.ok) {
      failNote = 'fold location chord solve failed (using best bracket point)';
      break;
    }
    solvedAny = true;
    const M = m.state;
    if (Math.abs(M.ts) < Math.abs(best.ts)) best = M;
    if (Math.abs(M.ts) <= 1e-8 * tsScale || (chordLen0 > 0 && m.chordLen <= 1e-9 * chordLen0)) {
      best = M;
      converged = true;
      break;
    }
    if ((M.ts > 0) === (A.ts > 0)) A = M;
    else B = M;
  }
  s._restore(best.snap);
  const drive = ctx.p0 + best.s * ctx.dpSpan;
  const { terminals } = readPointData(s, ctx.termIdx, drive, false);
  const located = converged && solvedAny;
  const note = located ? undefined : (failNote ?? 'fold location did not reach tolerance (using best bracket point)');
  return { drive, terminals, s: best.s, state: best, located, note };
}

function findComplexPair(st, omegaEst, refMode = null) {
  let bestVec = null, bestOv = 0.45;
  if (refMode && refMode._vr) {
    for (const m of st.modes) {
      if (!(m.im > 0)) continue;
      const ov = vecOverlap(refMode, m);
      if (ov > bestOv) {
        bestOv = ov;
        bestVec = m;
      }
    }
    if (bestVec) return bestVec;
  }
  let best = null, bestDist = Infinity;
  for (const m of st.modes) {
    if (!(m.im > 0)) continue;
    const dist = Math.abs(m.im - omegaEst);
    if (dist < bestDist) {
      bestDist = dist;
      best = m;
    }
  }
  return best;
}

function evalHopfMode(s, ctx, M, omegaEst, refMode = null, extraOmega = 0) {
  const drive = ctx.p0 + M.s * ctx.dpSpan;
  s._restore(M.snap);
  const fHz = omegaEst / (2 * Math.PI);
  const baseOpts = ctx.stabOpts;
  const freqs = baseOpts.frequencies ? [...baseOpts.frequencies, fHz] : [fHz];
  if (extraOmega > 0 && Math.abs(extraOmega - omegaEst) > 0.1 * Math.max(extraOmega, omegaEst)) {
    freqs.push(extraOmega / (2 * Math.PI));
  }
  const fastOpts = { ...baseOpts, shifts: [Math.max(1, 0.25 * omegaEst)], frequencies: freqs };
  let st = evaluateStability(s, ctx.termIdx, drive, fastOpts, ctx.isCurrentDrive, ctx.budget);
  let pair = findComplexPair(st, omegaEst, refMode);
  if (!pair) {
    st = evaluateStability(s, ctx.termIdx, drive, baseOpts, ctx.isCurrentDrive, ctx.budget);
    pair = findComplexPair(st, omegaEst, refMode);
  }
  return pair;
}

// Locate a Hopf bifurcation (Re λ = 0 for complex conjugate pair with Im λ > 0) between A0 and B0,
// using purely relative drive and chord tolerances (never absolute in drive units).
function locateHopf(s, ctx, A0, B0, pairA0, pairB0, omegaInit) {
  let A = A0, B = B0;
  let omega = omegaInit;
  let mA = pairA0 ?? evalHopfMode(s, ctx, A, omega, pairB0);
  let mB = pairB0 ?? evalHopfMode(s, ctx, B, omega, pairA0);
  if (!mA || !mB) return null;
  if ((mA.re > 0) === (mB.re > 0) && mA.re !== 0 && mB.re !== 0) return null;
  let best = Math.abs(mA.re) <= Math.abs(mB.re) ? { state: A, mode: mA } : { state: B, mode: mB };
  const kBar = ctx.stabOpts.kBar ?? 10;
  const dpBracket0 = Math.abs((B0.s - A0.s) * ctx.dpSpan);
  const chordLen0 = stateDistW(s, ctx, A0.snap.u, A0.snap.uLo, A0.snap.termV, A0.s, B0.snap.u, B0.snap.uLo, B0.snap.termV, B0.s);
  let converged = false;
  let solvedAny = false;
  let failNote = null;
  for (let it = 0; it < 20; it++) {
    ctx.budget?.check();
    const denom = Math.abs(mA.re) + Math.abs(mB.re);
    let theta = it % 3 === 2 || !(denom > 0) ? 0.5 : Math.abs(mA.re) / denom;
    if (theta < 0.15) theta = 0.15;
    if (theta > 0.85) theta = 0.85;
    let m = solveOnChord(s, ctx, A, B, theta);
    if (!m.ok && theta !== 0.5) m = solveOnChord(s, ctx, A, B, 0.5);
    if (!m.ok) {
      failNote = 'Hopf location chord solve failed (using best bracket point)';
      break;
    }
    const M = m.state;
    const refM = theta <= 0.5 ? mA : mB;
    const omegaTry = (1 - theta) * mA.im + theta * mB.im;
    const mM = evalHopfMode(s, ctx, M, omegaTry > 0 ? omegaTry : omega, refM, omega);
    if (!mM) {
      failNote = 'Hopf pair lost during bisection (using best bracket point)';
      break;
    }
    solvedAny = true;
    omega = mM.im;
    if (Math.abs(mM.re) < Math.abs(best.mode.re)) best = { state: M, mode: mM };
    const pA = ctx.p0 + A.s * ctx.dpSpan;
    const pB = ctx.p0 + B.s * ctx.dpSpan;
    const pScale = Math.max(Math.abs(pA), Math.abs(pB), dpBracket0);
    if (
      Math.abs(mM.re) <= Math.max(kBar * mM.error, 1e-5 * omega) ||
      (pScale > 0 && Math.abs(pB - pA) <= 1e-6 * pScale) ||
      (chordLen0 > 0 && m.chordLen <= 1e-9 * chordLen0)
    ) {
      best = { state: M, mode: mM };
      converged = true;
      break;
    }
    if ((mM.re > 0) === (mA.re > 0)) { A = M; mA = mM; }
    else { B = M; mB = mM; }
  }
  s._restore(best.state.snap);
  const drive = ctx.p0 + best.state.s * ctx.dpSpan;
  const { terminals } = readPointData(s, ctx.termIdx, drive, false);
  const located = converged && solvedAny;
  const note = located ? undefined : (failNote ?? 'Hopf location did not reach tolerance (using best bracket point)');
  return { drive, terminals, omega: best.mode.im, located, note };
}

function detectHopfBracket(stA, stB, kBar, stats = null) {
  const sig = (m) => (m.re > kBar * m.error ? 1 : m.re < -kBar * m.error ? -1 : 0);
  const hasIm = (m) => m.im > kBar * m.error;
  const matched = matchModes(stA, stB);
  for (const p of matched.pairs) {
    const { a: ma, b: mb } = p;
    if (!(hasIm(ma) && hasIm(mb))) continue;
    const sa = sig(ma), sb = sig(mb);
    if (sa !== 0 && sb !== 0 && sa !== sb) {
      if (stats) stats.hopfByVecMatch++;
      return { pairA: ma, pairB: mb, omega: 0.5 * (ma.im + mb.im) };
    }
  }
  if (stA.verdict === 'no growing mode found' && stB.verdict === 'unstable') {
    const grow = stB.modes.find((m) => m.re > kBar * m.error && hasIm(m));
    if (grow) {
      if (stats) stats.hopfByVecMatch++;
      return { pairA: findComplexPair(stA, grow.im, grow), pairB: grow, omega: grow.im };
    }
  }
  if (stA.verdict === 'unstable' && stB.verdict === 'no growing mode found') {
    const grow = stA.modes.find((m) => m.re > kBar * m.error && hasIm(m));
    if (grow) {
      if (stats) stats.hopfByVecMatch++;
      return { pairA: grow, pairB: findComplexPair(stB, grow.im, grow), omega: grow.im };
    }
  }
  return null;
}

function countPositiveReal(st, kBar) {
  let count = 0;
  for (const m of st.modes) if (m.im === 0 && m.re > kBar * m.error) count++;
  return count;
}

function formatDrive(val, kind, hasGeometry) {
  const unit = kind === 'V' ? 'V' : hasGeometry ? 'A' : 'A/m²';
  return `${Number(val.toPrecision(6))} ${unit}`;
}

export class Branch {
  /**
   * Trace a steady-state branch in one terminal's drive parameter using pseudo-arclength
   * continuation, marking each point with its linear stability and locating folds and Hopf
   * bifurcations.
   * @param {{ terminal?: string, to: number, step?: number, maxSteps?: number, voltageBound?: number, budgetMs?: number, stability?: boolean | object, solutions?: boolean }} opts
   */
  branch(opts = {}) {
    if (!this.atSteady || this.regriddedSteady) {
      throw new SolverError('branch: solve() first — branch() traces from the current converged steady state');
    }
    const {
      terminal = 'right',
      to,
      step,
      maxSteps = 200,
      voltageBound = 20,
      budgetMs = Infinity,
      stability: stabArg = true,
      solutions = false,
      _stats = null,
    } = opts;
    if (_stats) {
      Object.assign(_stats, {
        metricRebalanced: 0,
        rejectContraction: 0,
        rejectStateCap: 0,
        rejectAngle: 0,
        modeNearEnd: 0,
        hopfByVecMatch: 0,
        startUncheckedBudget: 0,
        auditOnBudget: 0,
      });
    }

    const termIdx = this.terms.findIndex((t) => t.name === terminal);
    if (termIdx < 0) {
      throw new SolverError(`branch: unknown terminal '${terminal}' (terminals: ${this.terms.map((t) => t.name).join(', ')})`);
    }
    if (typeof to !== 'number' || !Number.isFinite(to)) {
      throw new SolverError(`branch: 'to' must be a finite target drive value; got ${to}`);
    }
    if (step !== undefined && !(typeof step === 'number' && step > 0 && Number.isFinite(step))) {
      throw new SolverError(`branch: 'step' must be a positive finite number; got ${step}`);
    }
    if (!(Number.isInteger(maxSteps) && maxSteps >= 1)) {
      throw new SolverError(`branch: 'maxSteps' must be a positive integer; got ${maxSteps}`);
    }
    if (typeof voltageBound !== 'number' || !(voltageBound > 0)) {
      throw new SolverError(`branch: 'voltageBound' must be a positive number or Infinity; got ${voltageBound}`);
    }
    if (typeof budgetMs !== 'number' || !(budgetMs >= 0)) {
      throw new SolverError(`branch: 'budgetMs' must be a non-negative number or Infinity; got ${budgetMs}`);
    }
    if (typeof stabArg !== 'boolean' && (typeof stabArg !== 'object' || stabArg === null || Array.isArray(stabArg))) {
      throw new SolverError('branch: \'stability\' must be a boolean or stability() options object');
    }
    if (typeof solutions !== 'boolean') {
      throw new SolverError('branch: \'solutions\' must be a boolean');
    }

    const doStab = stabArg !== false;
    const stabOpts = typeof stabArg === 'object' ? { ...stabArg } : {};
    const kBar = stabOpts.kBar ?? 10;
    const nWant = stabOpts.modes ?? 8;
    const term = this.terms[termIdx];
    const kind = term.drive.kind === 'I' ? 'I' : 'V';
    const isCurrentDrive = kind === 'I';
    const hasGeometry = this.model.geometry.type !== 'planar';
    const p0 = sourceAt(term.drive.src, this.time);
    const warnings = [];
    if (this.terms.some((t) => Boolean(t.drive.src?.t))) {
      warnings.push(`branch: non-autonomous waveform drive frozen at t = ${this.time} s; drives are held at their values there.`);
    }

    // Snapshot everything so branch() is completely non-destructive (even on throw or budget stop).
    const snap0 = this._snapshot();
    const origSrc = term.drive.src;
    const origOverrides = new Map(this.sourceOverride);
    const origDrivenKey = this.drivenKey;
    const wasAll = this.allTerminals;
    const wasSteady = this.steady;
    const wasSteadyAt = this.steadyAt;
    const wasLawForm = this.lawForm;
    const wasSourceBefore = this.sourceBefore;
    const wasOffNeutral = this.offNeutral;
    const wasIslandsOn = this.islandsOn;
    const wasSolveAmounts = this.solveAmounts;
    const wasLocalKept = this.localKept;
    const wasCarried = this.carried;
    const wasLawRefs = this.lawRefs;
    const wasChargeRefs = this.chargeRefs;
    const wasSolvedV = this.solvedV ? this.solvedV.slice() : this.solvedV;
    const wasStepped = this.stepped;
    const wasColdStart = this.coldStart;
    const learnedKeep = { lostMode: this.lostMode, accepted: this.accepted, firstErrs: this.firstErrs?.slice() };

    const budget = Number.isFinite(budgetMs) ? branchBudget(budgetMs) : undefined;
    const points = [];
    const events = [];
    const stFullList = [];
    const states = [];
    let end = '';
    let ctx = null;

    // Post-trace consistency audit between geometric events and stability marks.
    const runAudit = (allowLocate) => {
      if (!doStab || points.length < 2) return;
      const isNeutral = (v) => v === 'neutral within resolution';
      const isDecided = (v) => v === 'no growing mode found' || v === 'unstable';
      const geomEvents = events.filter((e) => e.kind === 'fold' || e.kind === 'hopf');

      for (const ev of geomEvents) {
        if (ev.kind !== 'fold') continue;
        const [kL, kR] = ev.between;
        if (kL < 0 || kR >= points.length) continue;
        let minL = 0, maxR = points.length - 1;
        for (const other of geomEvents) {
          if (other === ev) continue;
          if (other.between[1] <= kL && other.between[1] > minL) minL = other.between[1];
          if (other.between[0] >= kR && other.between[0] < maxR) maxR = other.between[0];
        }
        let iL = kL;
        while (iL > minL && isNeutral(points[iL].stability.verdict) && kL - iL < 2) iL--;
        let iR = kR;
        while (iR < maxR && isNeutral(points[iR].stability.verdict) && iR - kR < 2) iR++;
        const stL = stFullList[iL], stR = stFullList[iR];
        if (!stL || !stR || !isDecided(stL.verdict) || !isDecided(stR.verdict)) continue;
        if (stL.verdict === stR.verdict && countPositiveReal(stL, kBar) === countPositiveReal(stR, kBar)) {
          events.push({
            kind: 'stability mismatch',
            drive: ev.drive,
            terminals: ev.terminals,
            between: ev.between,
            note: `stability marks did not change across fold (${stL.verdict} on both sides)`,
          });
        }
      }

      let prevDec = -1;
      for (let i = 0; i < points.length; i++) {
        const v = points[i].stability?.verdict;
        if (!isDecided(v)) continue;
        if (prevDec >= 0 && points[prevDec].stability.verdict !== v) {
          const covered = events.some(
            (ev) => (ev.kind === 'fold' || ev.kind === 'hopf') && ev.between[0] >= prevDec && ev.between[1] <= i,
          );
          if (!covered) {
            let resolved = false;
            const hopfBracket = stFullList[prevDec] && stFullList[i]
              ? detectHopfBracket(stFullList[prevDec], stFullList[i], kBar, _stats)
              : null;
            if (hopfBracket) {
              if (allowLocate && ctx && states[prevDec] && states[i]) {
                try {
                  const hopf = locateHopf(this, ctx, states[prevDec], states[i], hopfBracket.pairA, hopfBracket.pairB, hopfBracket.omega);
                  if (hopf) {
                    events.push({
                      kind: 'hopf',
                      drive: hopf.drive,
                      terminals: hopf.terminals,
                      omega: hopf.omega,
                      between: [prevDec, i],
                      ...(hopf.located ? {} : { located: false, note: hopf.note }),
                    });
                    resolved = true;
                  }
                } catch (e) {
                  if (e === BUDGET_SPENT) {
                    events.push({
                      kind: 'hopf',
                      drive: 0.5 * (points[prevDec].drive + points[i].drive),
                      terminals: points[i].terminals,
                      omega: hopfBracket.omega,
                      between: [prevDec, i],
                      located: false,
                      note: 'Hopf location cut short by budgetMs',
                    });
                    resolved = true;
                  } else {
                    throw e;
                  }
                }
              } else if (!allowLocate) {
                events.push({
                  kind: 'hopf',
                  drive: 0.5 * (points[prevDec].drive + points[i].drive),
                  terminals: points[i].terminals,
                  omega: hopfBracket.omega,
                  between: [prevDec, i],
                  located: false,
                  note: 'Hopf location cut short by budgetMs',
                });
                resolved = true;
              }
            }
            if (!resolved) {
              const mid = Math.floor((prevDec + i) / 2);
              events.push({
                kind: 'stability mismatch',
                drive: points[mid].drive,
                terminals: points[mid].terminals,
                between: [prevDec, i],
                note: `stability changed from '${points[prevDec].stability.verdict}' to '${v}' with no fold or Hopf located`,
              });
            }
          }
        }
        prevDec = i;
      }
      events.sort((a, b) => a.between[0] - b.between[0] || a.between[1] - b.between[1]);
    };

    try {
      // Push point 0 BEFORE stability() so start point is always kept even on budgetMs = 0.
      const pt0Data = readPointData(this, termIdx, p0, solutions);
      const pt0 = {
        drive: p0,
        terminals: pt0Data.terminals,
        ...(doStab ? {
          stability: {
            verdict: 'undecided',
            checked: false,
            rate: NaN,
            im: NaN,
            error: NaN,
            reason: 'not checked: budgetMs spent during start-point stability()',
          },
        } : {}),
        ...(solutions ? { solution: pt0Data.solution } : {}),
      };
      points.push(pt0);
      stFullList.push(null);
      if (doStab) {
        try {
          const st0 = evaluateStability(this, termIdx, p0, stabOpts, isCurrentDrive, budget);
          pt0.stability = summarizeStability(st0);
          stFullList[0] = st0;
        } catch (err) {
          if (err === BUDGET_SPENT && _stats) _stats.startUncheckedBudget++;
          throw err;
        }
      }

      if (to === p0) {
        return {
          terminal,
          kind,
          points,
          events,
          end: `reached to (${formatDrive(to, kind, hasGeometry)}; already at target)`,
          warnings,
        };
      }

      if (isCurrentDrive && p0 === 0) {
        term.drive.src = { value: to };
        this._lawsNow();
        term.drive.src = origSrc;
      } else {
        this._lawsNow();
      }
      if (!this._directSteady()) {
        return {
          terminal,
          kind,
          points,
          events,
          end: `refused: device has ${this.unkept.length} unkept conservation law(s) or an unconserved disconnected stretch (cannot continue steady equations directly)`,
          warnings,
        };
      }
      if (isCurrentDrive) {
        const why = chargingReason(this, termIdx, to);
        if (why) return { terminal, kind, points, events, end: `refused: ${why}`, warnings };
      }

      captureLaws(this);
      this.computeConcentrations();
      keepLocal(this);
      this.steady = true;
      this.steadyAt = this.time;
      this.allTerminals = true;

      const dpSpan = to - p0;
      const absSpan = Math.abs(dpSpan);
      const N = this.sys.size, fl = this.floating, K = fl.length;
      const { potIdx, potWt2 } = stateMetricEntries(this);
      const maxPins = Math.max(1, this.constraints.length + this.affinities.length + this.islands.length + 4);
      const work = {
        e: new Float64Array(N + 1),
        cols: [],
        zero: new Float64Array(N + 1),
        zeroPin: new Float64Array(maxPins),
        zeroCirc: new Float64Array(K),
        zeroWx: new Float64Array(N + 1),
        zeroWv: new Float64Array(K),
        delta: new Float64Array(N + 1),
        deltaV: new Float64Array(K),
        pinRhs: new Float64Array(maxPins),
        circRhs: new Float64Array(K),
        chordX: new Float64Array(N + 1),
        chordV: new Float64Array(K),
        locTx: new Float64Array(N + 1),
        locTV: new Float64Array(K),
        landTx: new Float64Array(N + 1),
        landTV: new Float64Array(K),
      };
      ctx = {
        termIdx,
        p0,
        dpSpan,
        potIdx,
        potWt2,
        wsEff: 1,
        Bs: new Float64Array(N + 1),
        pinDs: new Float64Array(maxPins),
        circDs: new Float64Array(K),
        wx: new Float64Array(N + 1),
        wV: new Float64Array(K),
        uPred: new Float64Array(this.u.length),
        uLoPred: new Float64Array(this.uLo.length),
        vPred: new Float64Array(K),
        work,
        budget,
        stabOpts,
        isCurrentDrive,
        startMissesKirchhoff: this.islands.length > 0 && this._missesKirchhoff(),
        startKirchhoffSum: this.islands.length > 0 ? kirchhoffSum(this) : 0,
      };

      const tx = new Float64Array(N + 1), tV = new Float64Array(K);
      const prevTx = new Float64Array(N + 1), prevTV = new Float64Array(K);
      const vStart = Float64Array.from(fl, (k) => this.termV[k]);

      // 1. Initial unit tangent at s = 0, then scale `wsEff` using the EXACT same state norm so
      //    state and drive start at an exact 50/50 split (`wsEff·|ts| = 1/√2`), 100% invariant to
      //    drive units and `to`.
      let tan0;
      try {
        tan0 = computeTangent(this, ctx, 0, null, null, null, tx, tV, false);
      } catch (err) {
        if (!(err instanceof Error)) throw err;
        return { terminal, kind, points, events, end: `refused: initial tangent failed (${err.message})`, warnings };
      }
      if (!tan0.ok) {
        return { terminal, kind, points, events, end: 'refused: initial tangent is non-finite or singular', warnings };
      }
      if (Math.abs(tan0.ts) > 1e-12) {
        const stateSlope = stateOnlyNorm(potIdx, potWt2, this.VT, tx, tV) / Math.abs(tan0.ts);
        if (stateSlope > 1e-12 && Number.isFinite(stateSlope)) {
          ctx.wsEff = stateSlope;
          tan0 = computeTangent(this, ctx, 0, null, null, null, tx, tV, true);
        }
      }

      let ts = tan0.ts;
      let sCur = 0;

      // Rebalance `wsEff` whenever the drive component drifts above 0.85 of the unit tangent while
      // the state moves, restoring 50/50 balance so opposite-state branches can never look parallel.
      const rebalanceMetric = () => {
        if (!(ctx.wsEff * Math.abs(ts) > 0.85)) return 1;
        const stNorm = stateOnlyNorm(potIdx, potWt2, this.VT, tx, tV);
        if (!(stNorm > 1e-12 && Math.abs(ts) > 1e-12)) return 1;
        ctx.wsEff = stNorm / Math.abs(ts);
        const newNrm = tangentNorm(potIdx, potWt2, ctx.wsEff, this.VT, tx, tV, ts);
        for (let j = 0; j < N; j++) tx[j] /= newNrm;
        for (let a = 0; a < K; a++) tV[a] /= newNrm;
        ts /= newNrm;
        // Keep the latest recorded state tangent consistent with the updated metric.
        const last = states[states.length - 1];
        if (last) {
          last.tx = Float64Array.from(tx);
          last.tV = Float64Array.from(tV);
          last.ts = ts;
        }
        if (_stats) _stats.metricRebalanced++;
        return newNrm;
      };

      states.push({ snap: this._snapshot(), s: 0, tx: Float64Array.from(tx), tV: Float64Array.from(tV), ts });

      // Initial stride in the balanced metric (1 unit ≈ 1 thermal/gate unit of motion).
      const maxPotPerStep = 6;
      const maxDVFloat = 0.25;
      const maxSPerStep = step !== undefined ? step / absSpan : 0.2;
      let dsArc = 2;
      if (Math.abs(ts) > 1e-12) dsArc = Math.min(dsArc, maxSPerStep / Math.abs(ts));
      if (isCurrentDrive) {
        const flIdx = fl.indexOf(termIdx);
        if (flIdx >= 0 && Math.abs(tV[flIdx]) > 1e-12) dsArc = Math.min(dsArc, 0.1 / Math.abs(tV[flIdx]));
      }
      const initStateRate = maxStateStep(this, tx, tV);
      if (initStateRate > 1e-12) dsArc = Math.min(dsArc, 5 / initStateRate);
      dsArc = Math.max(MIN_ARC, Math.min(MAX_ARC, dsArc));

      let geomFail = 0;
      let lastFailWasModeOnly = false;
      let stepCount = 0;

      const acceptPoint = (sNew, newTx, newTV, newTs, preSt = null, preFold = null) => {
        const kPrev = points.length - 1;
        const prevState = states[kPrev];
        const curSnap = this._snapshot();
        const curState = {
          snap: curSnap,
          s: sNew,
          tx: Float64Array.from(newTx),
          tV: Float64Array.from(newTV),
          ts: newTs,
        };
        const driveNew = p0 + sNew * dpSpan;
        const ptData = readPointData(this, termIdx, driveNew, solutions);
        let stNew = preSt;
        if (doStab && !stNew) {
          budget?.check();
          stNew = evaluateStability(this, termIdx, driveNew, stabOpts, isCurrentDrive, budget);
        }
        points.push({
          drive: driveNew,
          terminals: ptData.terminals,
          ...(doStab ? { stability: summarizeStability(stNew) } : {}),
          ...(solutions ? { solution: ptData.solution } : {}),
        });
        stFullList.push(stNew);
        states.push(curState);
        const kNew = points.length - 1;

        const hasFold = (prevState.ts > 0 && newTs < 0) || (prevState.ts < 0 && newTs > 0) || (newTs === 0 && prevState.ts !== 0);
        if (hasFold) {
          if (preFold) {
            events.push({
              kind: 'fold',
              drive: preFold.drive,
              terminals: preFold.terminals,
              between: [kPrev, kNew],
              ...(preFold.located ? {} : { located: false, note: preFold.note }),
            });
          } else {
            try {
              const fold = locateFold(this, ctx, prevState, curState);
              events.push({
                kind: 'fold',
                drive: fold.drive,
                terminals: fold.terminals,
                between: [kPrev, kNew],
                ...(fold.located ? {} : { located: false, note: fold.note }),
              });
            } catch (err) {
              if (err === BUDGET_SPENT) {
                const bestEnd = Math.abs(prevState.ts) <= Math.abs(newTs) ? points[kPrev] : points[kNew];
                events.push({
                  kind: 'fold',
                  drive: bestEnd.drive,
                  terminals: bestEnd.terminals,
                  between: [kPrev, kNew],
                  located: false,
                  note: 'fold location cut short by budgetMs',
                });
              }
              throw err;
            }
          }
        }

        if (doStab && stFullList[kPrev] && stNew) {
          const hopfBracket = detectHopfBracket(stFullList[kPrev], stNew, kBar, _stats);
          if (hopfBracket) {
            try {
              const hopf = locateHopf(this, ctx, prevState, curState, hopfBracket.pairA, hopfBracket.pairB, hopfBracket.omega);
              if (hopf) {
                events.push({
                  kind: 'hopf',
                  drive: hopf.drive,
                  terminals: hopf.terminals,
                  omega: hopf.omega,
                  between: [kPrev, kNew],
                  ...(hopf.located ? {} : { located: false, note: hopf.note }),
                });
              }
            } catch (err) {
              if (err === BUDGET_SPENT) {
                events.push({
                  kind: 'hopf',
                  drive: 0.5 * (points[kPrev].drive + points[kNew].drive),
                  terminals: points[kNew].terminals,
                  omega: hopfBracket.omega,
                  between: [kPrev, kNew],
                  located: false,
                  note: 'Hopf location cut short by budgetMs',
                });
              }
              throw err;
            }
          }
        }

        this._restore(curSnap);
        this.sourceOverride.set(termIdx, driveNew);
      };

      const landExact = (A, B, sTarget, refTx, refTV, refTs) => {
        const frac = B.s !== A.s ? (sTarget - A.s) / (B.s - A.s) : 1;
        const fClamped = Math.max(0, Math.min(1, frac));
        this._restore(A.snap);
        for (let o = 0; o < this.u.length; o++) {
          this.u[o] = A.snap.u[o] + fClamped * (B.snap.u[o] - A.snap.u[o]);
          this.uLo[o] = A.snap.uLo[o] + fClamped * (B.snap.uLo[o] - A.snap.uLo[o]);
        }
        for (let a = 0; a < K; a++) {
          this.termV[fl[a]] = A.snap.termV[fl[a]] + fClamped * (B.snap.termV[fl[a]] - A.snap.termV[fl[a]]);
        }
        ctx.uPred.set(this.u);
        ctx.uLoPred.set(this.uLo);
        for (let a = 0; a < K; a++) ctx.vPred[a] = this.termV[fl[a]];
        const predTermV = Float64Array.from(this.termV);
        ctx.wx.fill(0);
        ctx.wV.fill(0);
        const corLand = correctorNewton(this, ctx, sTarget, ctx.uPred, ctx.uLoPred, ctx.vPred, sTarget, ctx.wx, ctx.wV, 1, { fixedS: true });
        if (!corLand.converged || corLand.damped) return null;

        for (let a = 0; a < K; a++) {
          if (Math.abs(this.termV[fl[a]] - vStart[a]) > voltageBound) return null;
        }
        const lastPt = states[points.length - 1];
        if (step !== undefined && Math.abs((sTarget - lastPt.s) * dpSpan) > step * (1 + 1e-6)) return null;

        const bracketLen = stateDistW(this, ctx, A.snap.u, A.snap.uLo, A.snap.termV, A.s, B.snap.u, B.snap.uLo, B.snap.termV, B.s);
        const dCor = stateDistW(this, ctx, ctx.uPred, ctx.uLoPred, predTermV, sTarget, this.u, this.uLo, this.termV, sTarget);
        if (bracketLen > 1e-12 && dCor > (bracketLen > 1e-3 ? 0.75 : 2.0) * bracketLen) return null;
        if (maxStateDiff(this, lastPt.snap, this.u, this.uLo, this.termV) > maxPotPerStep * 1.5) return null;

        const tanLand = computeTangent(this, ctx, sTarget, refTx, refTV, refTs, work.landTx, work.landTV, !!corLand.reusable);
        if (!tanLand.ok || tanLand.cosTheta < COS_MIN) return null;
        if (!(stepAlignment(this, ctx, refTx, refTV, refTs, A.snap, A.s, this.u, this.uLo, this.termV, sTarget) > 0)) return null;

        return { s: sTarget, tx: Float64Array.from(work.landTx), tV: Float64Array.from(work.landTV), tan: tanLand };
      };

      // 2. Main pseudo-arclength continuation loop.
      while (stepCount < maxSteps) {
        budget?.check();

        let dsStep = dsArc;
        if (Math.abs(ts) > 1e-12) dsStep = Math.min(dsStep, maxSPerStep / Math.abs(ts));
        const stateRate = maxStateStep(this, tx, tV);
        if (stateRate > 1e-12) dsStep = Math.min(dsStep, maxPotPerStep / stateRate);
        for (let a = 0; a < K; a++) {
          const vt = Math.abs(tV[a]);
          if (vt > 1e-12) dsStep = Math.min(dsStep, maxDVFloat / vt);
        }
        if (dsStep < MIN_ARC) dsStep = MIN_ARC;

        if (sCur < 1 && ts > 1e-12 && sCur + dsStep * ts > 1.02) {
          dsStep = Math.max(MIN_ARC, (1.01 - sCur) / ts);
        }

        const beforeSnap = this._snapshot();
        const prevStateObj = states[points.length - 1];
        const sPred = sCur + dsStep * ts;
        const maxPt = stateRate * dsStep;
        const predAlpha = maxPt > 10 ? 10 / maxPt : 1;

        addToState(this, tx, dsStep);
        for (let a = 0; a < K; a++) this.termV[fl[a]] += dsStep * tV[a];
        ctx.uPred.set(this.u);
        ctx.uLoPred.set(this.uLo);
        for (let a = 0; a < K; a++) ctx.vPred[a] = this.termV[fl[a]];
        if (predAlpha < 1) {
          this._restore(beforeSnap);
          addToState(this, tx, dsStep * predAlpha);
          for (let a = 0; a < K; a++) this.termV[fl[a]] += dsStep * predAlpha * tV[a];
        }
        const sInit = sCur + dsStep * predAlpha * ts;
        const ws = buildWeights(potIdx, potWt2, ctx.wsEff, this.VT, tx, tV, ts, ctx.wx, ctx.wV);

        const cor = correctorNewton(this, ctx, sInit, ctx.uPred, ctx.uLoPred, ctx.vPred, sPred, ctx.wx, ctx.wV, ws, {
          maxIter: geomFail > 0 ? 20 : 35,
        });

        const rejectStep = (factor, reasonIfMin, isModeOnly = false) => {
          this._restore(beforeSnap);
          if (isModeOnly) {
            lastFailWasModeOnly = true;
          } else {
            geomFail = 2;
            lastFailWasModeOnly = false;
          }
          dsArc = dsStep * factor;
          if (dsArc < MIN_ARC) {
            end = `${reasonIfMin} (solve() may continue from here)`;
            return false;
          }
          return true;
        };

        if (!cor.converged) {
          const msg = cor.modelError
            ? `continuation stopped: model refused (${cor.modelError})`
            : `continuation stopped: corrector failed at ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)} (${cor.error ?? 'step shrunk below minimum'})`;
          if (!rejectStep(0.5, msg)) break;
          continue;
        }

        const sNext = cor.s;

        // Guard 1: Corrector contraction (Deuflhard / Allgower-Georg) + state & hard step cap.
        if (cor.damped || (cor.d1 > 1e-6 && cor.kappa > KAPPA_MAX)) {
          if (_stats) _stats.rejectContraction++;
          if (!rejectStep(0.5, `continuation stopped: corrector did not contract near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`)) break;
          continue;
        }
        const actualStateMove = maxStateDiff(this, beforeSnap, this.u, this.uLo, this.termV);
        const actualDp = Math.abs((sNext - sCur) * dpSpan);
        const stepExceeded = step !== undefined && actualDp > step * (1 + 1e-6);
        if (stepExceeded || actualStateMove > maxPotPerStep * 1.25) {
          if (_stats) _stats.rejectStateCap++;
          const factor = stepExceeded && actualDp > 0 ? Math.min(0.5, 0.8 * (step / actualDp)) : 0.5;
          if (!rejectStep(factor, `continuation stopped: step control could not bound state step near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`)) break;
          continue;
        }

        // Guard 2: Tangent continuity in the balanced metric + forward alignment.
        prevTx.set(tx);
        prevTV.set(tV);
        const prevTs = ts;
        let tanNext;
        try {
          tanNext = computeTangent(this, ctx, sNext, prevTx, prevTV, prevTs, tx, tV, !!cor.reusable);
        } catch (err) {
          if (!(err instanceof Error)) throw err;
          tanNext = { ok: false };
        }
        if (!tanNext.ok || tanNext.cosTheta < COS_MIN) {
          if (_stats) _stats.rejectAngle++;
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          if (!rejectStep(0.5, `continuation stopped: tangent turned sharply or became singular near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`)) break;
          continue;
        }
        if (!(stepAlignment(this, ctx, prevTx, prevTV, prevTs, beforeSnap, sCur, this.u, this.uLo, this.termV, sNext) > 0)) {
          if (_stats) _stats.rejectAngle++;
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          if (!rejectStep(0.5, `continuation stopped: corrector jumped backward near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`)) break;
          continue;
        }

        // Guard 3: Voltage bound on floating terminals.
        let hitVBound = null;
        for (let a = 0; a < K; a++) {
          const dv = Math.abs(this.termV[fl[a]] - vStart[a]);
          if (dv > voltageBound) {
            hitVBound = { name: this.terms[fl[a]].name, V: this.termV[fl[a]], dv };
            break;
          }
        }
        if (hitVBound) {
          this._restore(beforeSnap);
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          end = `left voltage bound (±${voltageBound} V from start at terminal '${hitVBound.name}', reached ${hitVBound.V.toPrecision(4)} V)`;
          break;
        }

        // Check for fold FIRST so the One Span Rule handles folds past `to` and immediate return
        // past `p0` uniformly.
        const candState = {
          snap: this._snapshot(),
          s: sNext,
          tx: Float64Array.from(tx),
          tV: Float64Array.from(tV),
          ts: tanNext.ts,
        };
        const stepFold = (prevTs > 0 && tanNext.ts < 0) || (prevTs < 0 && tanNext.ts > 0);
        let preFold = null;
        if (stepFold) {
          try {
            preFold = locateFold(this, ctx, prevStateObj, candState);
          } catch (err) {
            if (err === BUDGET_SPENT && !doStab) {
              const kPrev = points.length - 1;
              const candDrive = p0 + candState.s * dpSpan;
              this._restore(candState.snap);
              const ptData = readPointData(this, termIdx, candDrive, solutions);
              points.push({ drive: candDrive, terminals: ptData.terminals, ...(solutions ? { solution: ptData.solution } : {}) });
              const bestEnd = Math.abs(prevTs) <= Math.abs(tanNext.ts) ? points[kPrev] : points[kPrev + 1];
              events.push({
                kind: 'fold',
                drive: bestEnd.drive,
                terminals: bestEnd.terminals,
                between: [kPrev, kPrev + 1],
                located: false,
                note: 'fold location cut short by budgetMs',
              });
            }
            throw err;
          }
          this._restore(candState.snap);

          if (prevTs > 0 && sCur < 1 && preFold.s >= 1) {
            const landed = landExact(prevStateObj, preFold.state, 1, prevTx, prevTV, prevTs);
            if (landed) {
              acceptPoint(1, landed.tx, landed.tV, landed.tan.ts);
              end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
              break;
            }
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            if (!rejectStep(0.5, `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`)) break;
            continue;
          }
          if (prevTs < 0 && sCur > 0 && preFold.s <= 0) {
            const landed = landExact(prevStateObj, preFold.state, 0, prevTx, prevTV, prevTs);
            if (landed) {
              acceptPoint(0, landed.tx, landed.tV, landed.tan.ts);
              const vNote = isCurrentDrive ? `, V = ${points[points.length - 1].terminals[terminal].V.toPrecision(5)} V` : '';
              end = `returned to starting drive (${formatDrive(p0, kind, hasGeometry)}${vNote}) after fold`;
              break;
            }
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            if (!rejectStep(0.5, `continuation stopped: could not land on starting drive ${formatDrive(p0, kind, hasGeometry)}`)) break;
            continue;
          }
        }

        // One Span Rule — edge s = 1 (`to`).
        if (sCur < 1 && sNext >= 1) {
          const landed = landExact(prevStateObj, candState, 1, prevTx, prevTV, prevTs);
          if (landed) {
            let stLand = null;
            if (doStab) {
              budget?.check();
              stLand = evaluateStability(this, termIdx, to, stabOpts, isCurrentDrive, budget);
              const sm = modeStepSmooth(stFullList[points.length - 1], stLand, kBar, nWant, _stats);
              if (!sm.ok) {
                tx.set(prevTx);
                tV.set(prevTV);
                ts = prevTs;
                if (!rejectStep(0.5, `continuation stopped: mode guard (${sm.reason}) near ${formatDrive(to, kind, hasGeometry)}`, true)) break;
                continue;
              }
            }
            acceptPoint(1, landed.tx, landed.tV, landed.tan.ts, stLand);
            end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
            break;
          }
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          if (!rejectStep(0.25, `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`)) break;
          continue;
        }

        // One Span Rule — edge s = 0 (`p0`, returning after a fold).
        if (sNext <= 0 && (sCur > 0 || (stepFold && preFold && preFold.s > 0))) {
          const bracketA = stepFold && preFold ? preFold.state : prevStateObj;
          const landed = landExact(bracketA, candState, 0, bracketA.tx, bracketA.tV, bracketA.ts);
          if (landed) {
            let stLand = null;
            if (doStab) {
              budget?.check();
              stLand = evaluateStability(this, termIdx, p0, stabOpts, isCurrentDrive, budget);
              const sm = modeStepSmooth(stFullList[points.length - 1], stLand, kBar, nWant, _stats);
              if (!sm.ok) {
                tx.set(prevTx);
                tV.set(prevTV);
                ts = prevTs;
                if (!rejectStep(0.5, `continuation stopped: mode guard (${sm.reason}) near ${formatDrive(p0, kind, hasGeometry)}`, true)) break;
                continue;
              }
            }
            acceptPoint(0, landed.tx, landed.tV, landed.tan.ts, stLand, preFold);
            const vNote = isCurrentDrive ? `, V = ${points[points.length - 1].terminals[terminal].V.toPrecision(5)} V` : '';
            end = `returned to starting drive (${formatDrive(p0, kind, hasGeometry)}${vNote}) after fold`;
            break;
          }
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          if (!rejectStep(0.25, `returned past starting drive (${formatDrive(p0, kind, hasGeometry)}) after fold`)) break;
          continue;
        }

        // Guard 4: Single principled nearer-end mode guard across all matched active modes.
        let stCand = null;
        if (doStab) {
          budget?.check();
          const driveCand = p0 + sNext * dpSpan;
          stCand = evaluateStability(this, termIdx, driveCand, stabOpts, isCurrentDrive, budget);
          const sm = modeStepSmooth(stFullList[points.length - 1], stCand, kBar, nWant, _stats);
          if (!sm.ok) {
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            if (!rejectStep(0.5, `continuation stopped: mode guard (${sm.reason}) near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`, true)) break;
            continue;
          }
        }

        const wasModeHalved = lastFailWasModeOnly;
        lastFailWasModeOnly = false;

        // Accept this point!
        ts = tanNext.ts;
        sCur = sNext;
        stepCount++;
        acceptPoint(sNext, tx, tV, ts, stCand, preFold);

        // Adapt arclength stride from Newton contraction factor and tangent smoothness.
        if (geomFail > 0) {
          geomFail--;
          dsArc = dsStep;
        } else if (cor.kappa < 0.08 && cor.iterations <= 4 && tanNext.cosTheta > 0.94) {
          dsArc = Math.min(MAX_ARC, dsStep * (wasModeHalved ? 2.0 : 1.6));
        } else if (cor.kappa < 0.25 && cor.iterations <= 6 && tanNext.cosTheta > 0.88) {
          dsArc = Math.min(MAX_ARC, dsStep * (wasModeHalved ? 1.8 : 1.3));
        } else if (cor.kappa > 0.35 || cor.iterations >= 10 || tanNext.cosTheta < 0.82) {
          dsArc = Math.max(MIN_ARC, dsStep * 0.7);
        } else {
          dsArc = dsStep;
        }

        // Keep state and drive balanced if the drive component has drifted above 0.85.
        const rebalScale = rebalanceMetric();
        if (rebalScale !== 1) dsArc = Math.max(MIN_ARC, Math.min(MAX_ARC, dsArc * rebalScale));
      }

      if (!end && stepCount >= maxSteps) {
        end = `maxSteps reached (${maxSteps} steps, stopped at ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)})`;
      }

      // 3. Post-trace consistency audit between geometric events and stability marks.
      runAudit(true);

      return { terminal, kind, points, events, end, warnings };
    } catch (err) {
      if (err === BUDGET_SPENT) {
        if (doStab && points.length >= 2) {
          if (_stats) _stats.auditOnBudget++;
          runAudit(false);
        }
        const el = Number(budget.elapsed().toPrecision(3));
        return {
          terminal,
          kind,
          points,
          events,
          end: `budget spent after ${el} ms (budgetMs: ${budgetMs})`,
          warnings,
        };
      }
      throw err;
    } finally {
      term.drive.src = origSrc;
      this.sourceOverride.clear();
      for (const [k, v] of origOverrides) this.sourceOverride.set(k, v);
      this.allTerminals = wasAll;
      this.steady = wasSteady;
      this.steadyAt = wasSteadyAt;
      this.lawForm = wasLawForm;
      this.sourceBefore = wasSourceBefore;
      this.offNeutral = wasOffNeutral;
      this.islandsOn = wasIslandsOn;
      this.solveAmounts = wasSolveAmounts;
      this.localKept = wasLocalKept;
      this.carried = wasCarried;
      this.lawRefs = wasLawRefs;
      this.chargeRefs = wasChargeRefs;
      this.solvedV = wasSolvedV;
      this.stepped = wasStepped;
      this.coldStart = wasColdStart;
      this.lostMode = learnedKeep.lostMode;
      this.accepted = learnedKeep.accepted;
      this.firstErrs = learnedKeep.firstErrs;
      if (this.drivenKey !== origDrivenKey) this._lawsNow();
      this._restore(snap0);
      this._refreshSources();
    }
  }
}
