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
//     initially with ds/d(arc) > 0. The arclength metric scales s by the device's own physical
//     response scale at p0 (so a far `to` does not inflate the initial step) and weights gate
//     variables (in [0, 1]) by 10× alongside potentials in thermal units. Bordered Woodbury solves
//     reuse s.sys (pure JS, zero dependencies, no ** or Math.pow).
//   - Every output point is converged to `tol = 1e-10` thermal units. On devices with floating
//     islands (`this.islands.length > 0`), if the plain corrector leaves Kirchhoff unbalanced, a
//     second-phase islands round (`islandsOn = true`, fixed s) polishes the island levels (as
//     newton() does).
//   - Heuristics guarding step continuity (applied to ordinary steps and landing steps alike; they
//     prevent gross leaps and missed narrow windows, not a mathematical proof against jumping onto
//     an adjacent parallel branch):
//       1. Corrector displacement: `||z_next − z_pred||_W ≤ 0.5 · dsStep` and
//          `||z_next − z_cur||_W ≤ 1.5 · dsStep`, plus pointwise state motion bounded.
//       2. Tangent turning angle: `cos θ ≥ 0.75` (≤ 41°) and positive forward alignment.
//       3. Hard `step` cap: when `opts.step` is given, `|p_next − p_cur| ≤ step` on every step.
//       4. Rightmost-mode smoothness (when `stability !== false`): arclength-extrapolated Re λ
//          must not predict a missed zero crossing or disagree wildly with observed Re λ, and
//          single-step mode magnitude / complex-pair frequency jumps are capped.
//
// One span rule:
//   - The trace stays inside the requested drive span `s ∈ [0, 1]` (between start drive `p0` and
//     `to`). It ends where it first leaves that span, landing cleanly on the crossed edge (`to` or
//     `p0`) and naming which edge in `end`. When a step turns a fold, the fold's drive is located
//     first: if the fold lies past `to` (sFold ≥ 1), the trace lands on `to` on the pre-fold leg
//     with no fold event; if a returning leg crosses back through `p0` (even on the same step as
//     the fold), it records the fold and lands on `p0` on the far branch.
//
// Stability marks & events:
//   - With `stability !== false`, every accepted point is evaluated by a full `this.stability(...)`
//     call at that point's state and drive, agreeing with a standalone `dev.stability()` call at
//     the same state. Marks are for the drive as held (a curve can be stable under V-drive and
//     unstable under I-drive).
//   - Fold (`dp/ds` changes sign): located by safeguarded secant/bisection along the exact arc to
//     `|ts| < 1e-8` (or chord < 1e-9), reporting `{ kind: 'fold', drive, terminals, between }`.
//   - Hopf (`Re λ` of a complex conjugate pair `Im λ > 0` crosses 0): located by safeguarded
//     secant/bisection along the exact arc to `|Re λ| ≤ max(kBar·error, 1e-5·ω)` or
//     `|Δp| ≤ 1e-6·max(1, |p|)`, reporting `{ kind: 'hopf', drive, terminals, omega, between }`.
//   - Fallbacks are always reported: if a fold or Hopf location fails to converge or is cut short
//     by `budgetMs`, the event is still emitted with `located: false` and a `note` saying why.
//   - Consistency audit: at a fold in the held parameter, one real eigenvalue must pass through 0.
//     If stability marks do not change across a fold (without walking across any neighbouring
//     event), or flip between decided states with no fold or Hopf located, a
//     `{ kind: 'stability mismatch', ... }` event is emitted (never silently ignored).
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

const WS_BASE = 10; // base arclength weight on normalized drive motion
const GATE_WEIGHT = 10; // weight on [0, 1] gate fractions (Δgate = 0.1 ↔ 1 thermal unit)
const BRANCH_TOL = 1e-10; // corrector Newton tolerance in thermal units (~2.5e-12 V)
const MIN_ARC = 1e-6;
const MAX_ARC = 1e5;
const COS_MIN = 0.75; // max tangent turn ≈ 41° per step (heuristic continuity guard)
const BUDGET_SPENT = Symbol('branch budget spent');
const now = () => (globalThis.performance ? globalThis.performance.now() : Date.now());

function branchBudget(ms) {
  const start = now();
  let attempted = false;
  let elapsed = 0;
  return {
    check() {
      if (attempted && (elapsed = now() - start) >= ms) throw BUDGET_SPENT;
    },
    // Check inside stability() shifts/refinements: allows at least one shift when ms = 0.
    checkWork() {
      if (attempted && (elapsed = now() - start) >= ms) throw BUDGET_SPENT;
      attempted = true;
    },
    markAttempted() {
      attempted = true;
    },
    elapsed() {
      return elapsed || now() - start;
    },
  };
}

// Dimensionless potential/level unknowns + per-entry weight (1 for potentials/levels, GATE_WEIGHT
// for gate fractions in [0, 1]).
function stateMetricEntries(s) {
  const { n, M, rix: R } = s;
  const idx = [];
  const wt = [];
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (s.nodeConductor[g] < 0 && s.loc[b * M] >= 0) { idx.push(R[b * M]); wt.push(1); }
    for (let i = 0; i < n; i++) {
      if (s.present[g * n + i] && s.loc[b * M + 1 + i] >= 0) { idx.push(R[b * M + 1 + i]); wt.push(1); }
    }
    const k = s.surfAt[g];
    if (k >= 0) {
      for (let q = 0; q < s.model.surfaces[k].species.length; q++) {
        if (s.loc[b * M + 1 + n + q] >= 0) { idx.push(R[b * M + 1 + n + q]); wt.push(1); }
      }
    }
  }
  for (const gt of s.gateList) {
    if (s.loc[gt.o] >= 0) { idx.push(R[gt.o]); wt.push(GATE_WEIGHT); }
  }
  for (const d of s.diffList) {
    const bf = s.blockOfFace[d.f];
    if (s.loc[d.o] >= 0) { idx.push(R[d.o]); wt.push(1); }
    for (const q of d.drawSlot) if (q >= 0 && s.loc[bf * M + q] >= 0) { idx.push(R[bf * M + q]); wt.push(1); }
  }
  return { potIdx: Int32Array.from(idx), potWt: Float64Array.from(wt) };
}

// Pointwise max state change: maxPotentialStep with gate variables weighted by GATE_WEIGHT.
function maxStateStep(s, dx, dV) {
  let mx = maxPotentialStep(s, dx, dV);
  const R = s.rix;
  for (const gt of s.gateList) {
    const v = GATE_WEIGHT * Math.abs(dx[R[gt.o]]);
    if (v > mx) mx = v;
  }
  return mx;
}

// Pointwise max state difference between two snapshots (or snapshot vs current solver state).
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

// Weighted arclength norm of difference between (uB, uLoB, termVB, sB) and (uA, uLoA, termVA, sA).
function stateDistW(s, ctx, uA, uLoA, termVA, sA, uB, uLoB, termVB, sB) {
  const { potIdx, potWt, wsEff } = ctx;
  const fl = s.floating, VT = s.VT;
  let sum = wsEff * wsEff * (sB - sA) * (sB - sA);
  const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i], o = s.fullOf[j], w = potWt[i];
    const d = (uB[o] - uA[o]) + (uLoB[o] - uLoA[o]);
    sum += invN * w * w * d * d;
  }
  for (let a = 0; a < fl.length; a++) {
    const k = fl[a];
    const d = (termVB[k] - termVA[k]) / VT;
    sum += d * d;
  }
  return Math.sqrt(sum);
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
    if (a >= 0) {
      if (d.kind === 'I') circDs[a] = -dpSpan;
      else circDs[a] = -dpSpan / d.R;
    }
  }
}

// Solve the (P + K + 1)-bordered system reusing the factored s.sys, returning
// { ds, minSchurPivot } and writing δx into `delta`, δV into `deltaV`.
function solveBordered(s, work, rhs, pinRhs, circRhs, gRhs, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws) {
  const N = s.sys.size, fl = s.floating, pins = s.pins, P = pins.length, K = fl.length, S = P + K + 1;
  s._solveLinear(rhs, delta);
  const e = work.e;
  const cols = work.cols;
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
  let minSchurPivot = Infinity;
  for (let col = 0; col < S; col++) {
    let p = col;
    for (let r = col + 1; r < S; r++) if (Math.abs(A[r][col]) > Math.abs(A[p][col])) p = r;
    if (p !== col) [A[col], A[p]] = [A[p], A[col]];
    const piv = A[col][col];
    if (!(piv !== 0 && Number.isFinite(piv))) throw new SolverError('the bordered arclength system is singular');
    const apiv = Math.abs(piv);
    if (apiv < minSchurPivot) minSchurPivot = apiv;
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
  return { ds: x[P + K], minSchurPivot };
}

function buildWeights(potIdx, potWt, wsEff, VT, tx, tV, ts, wx, wV) {
  const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
  wx.fill(0);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i], w = potWt[i];
    wx[j] = invN * w * w * tx[j];
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tV.length; a++) wV[a] = invVT2 * tV[a];
  return wsEff * wsEff * ts;
}

function tangentNorm(potIdx, potWt, wsEff, VT, tx, tV, ts) {
  let sum = wsEff * wsEff * ts * ts;
  const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
  for (let i = 0; i < potIdx.length; i++) {
    const w = potWt[i];
    const v = w * tx[potIdx[i]];
    sum += invN * v * v;
  }
  for (let a = 0; a < tV.length; a++) {
    const v = tV[a] / VT;
    sum += v * v;
  }
  return Math.sqrt(sum);
}

// Compute unit tangent (tx, tV, ts) at current converged state + sVal, oriented to agree with
// (prevTx, prevTV, prevTs), or with ts > 0 when prevTx is null.
function computeTangent(s, ctx, sVal, prevTx, prevTV, prevTs, tx, tV, reuseFactor = false) {
  const { termIdx, p0, dpSpan, potIdx, potWt, wsEff, Bs, pinDs, circDs, wx, wV, work } = ctx;
  const N = s.sys.size, K = s.floating.length;
  if (!reuseFactor) {
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    s.islandsOn = false;
    s.assemble(Infinity);
    equilibrate(s);
    computeParamDs(s, termIdx, dpSpan, Bs, pinDs, circDs);
    s._factor();
  }
  const zero = work.zero;
  const zeroPin = work.zeroPin;
  const zeroCirc = work.zeroCirc;
  let ws;
  if (prevTx) {
    ws = buildWeights(potIdx, potWt, wsEff, s.VT, prevTx, prevTV, prevTs, wx, wV);
  } else {
    wx.fill(0);
    wV.fill(0);
    ws = 1;
  }
  const res = solveBordered(s, work, zero, zeroPin, zeroCirc, 1, tx, tV, Bs, pinDs, circDs, wx, wV, ws);
  const tsRaw = res.ds;
  if (!Number.isFinite(tsRaw)) return { ok: false };
  const nrm = tangentNorm(potIdx, potWt, wsEff, s.VT, tx, tV, tsRaw);
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
  return { ok: true, ts, cosTheta, minSchurPivot: res.minSchurPivot };
}

// Corrector Newton on F(x, V_fl, s) = 0 and wx·(x − xPred) + wV·(V − vPred) + ws·(s − sPred) = 0.
// On floating-island devices, if plain convergence leaves Kirchhoff unbalanced, runs a second-phase
// fixed-s islands polish (`islandsOn = true`) mirroring newton()'s islands round.
function correctorNewton(s, ctx, sStart, uPred, uLoPred, vPred, sPred, wx, wV, ws, { maxIter = 50, tol = BRANCH_TOL, maxStep = 10, fixedS = false } = {}) {
  const { termIdx, p0, dpSpan, potIdx, wsEff, Bs, pinDs, circDs, work, budget } = ctx;
  budget?.check();
  const N = s.sys.size, fl = s.floating, K = fl.length;
  const delta = work.delta, deltaV = work.deltaV;
  const pinRhs = work.pinRhs, circRhs = work.circRhs;
  let sVal = fixedS ? sPred : sStart;

  const runLoop = (islandsOn, isFixed, limit) => {
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
      const alpha = own > maxStep ? maxStep / own : 1;
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
          reusable: !islandsOn && s.sys.perturbedPivots === 0 && s.sys.lostColumns === 0,
          minSchurPivot: solved.minSchurPivot,
        };
      }
    }
    return { converged: false, iterations: limit, s: sVal, error: 'Newton did not converge' };
  };

  const r1 = runLoop(false, fixedS, maxIter);
  if (!r1.converged) return r1;

  // Phase 2: on floating-island devices, if plain convergence left Kirchhoff unbalanced, run an
  // islands round at fixed s (matching newton()'s islands round) to settle weakly held islands.
  if (s.islands.length > 0 && s._missesKirchhoff()) {
    const preIslSnap = s._snapshot();
    for (const isl of s.islands) isl.cutAt = null;
    const r2 = runLoop(true, true, 25);
    s.islandsOn = false;
    if (r2.converged && !s._missesKirchhoff()) {
      return {
        converged: true,
        iterations: r1.iterations + r2.iterations,
        s: sVal,
        reusable: false,
        minSchurPivot: r2.minSchurPivot,
      };
    }
    // If direct islands round failed or still missed Kirchhoff, try one full _steadySteps polish
    // (which includes GMRES refinement) at fixed sVal.
    s._restore(preIslSnap);
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    const qPol = s._steadySteps({ maxSteps: 1 });
    if (qPol.converged && !s._missesKirchhoff()) {
      return {
        converged: true,
        iterations: r1.iterations + qPol.iterations,
        s: sVal,
        reusable: false,
        minSchurPivot: r1.minSchurPivot,
      };
    }
    s._restore(preIslSnap);
    return { converged: false, iterations: r1.iterations + r2.iterations, s: sVal, error: 'floating-island current balance lost' };
  }
  return r1;
}

// Check whether terminal `k` (current-driven) is forced to charge a capacitance or a conserved
// pool at non-zero current `I`, meaning no steady branch exists as `I` varies.
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

// Read terminal voltages and currents (into the device, matching sol.terminals) and optionally
// build a full Solution at drive value `p`.
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

// Run full stability() at the current state and drive `p`, passing any active budget check so
// overruns stay bounded to roughly one shift factorisation / mode refinement.
function evaluateStability(s, termIdx, p, stabOpts, isCurrentDrive, budget) {
  const term = s.terms[termIdx];
  const origSrc = term.drive.src;
  const hadOverride = s.sourceOverride.has(termIdx);
  const oldOverride = hadOverride ? s.sourceOverride.get(termIdx) : undefined;
  const wasSteady = s.steady, wasAtSteady = s.atSteady, wasRegrid = s.regriddedSteady;
  const wasDrivenKey = s.drivenKey;
  const wasLawRefs = s.lawRefs, wasChargeRefs = s.chargeRefs;
  const wasLocalKept = s.localKept;
  term.drive.src = { value: p };
  s.sourceOverride.delete(termIdx);
  s.steady = false;
  s.atSteady = true;
  s.regriddedSteady = false;
  let relawed = false;
  const callOpts = budget ? { ...stabOpts, _budgetCheck: () => budget.checkWork() } : stabOpts;
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

// Heuristic smoothness check on rightmost modes between consecutive states (prev -> next, with
// optional prev2 for arclength linear extrapolation). Returns true if modes connect smoothly.
function modeStepSmooth(stPrev2, stPrev, stNew, dPrevArc, dActual, kBar) {
  if (!stPrev || !stNew) return true;
  if (stPrev.verdict === 'undecided' || stNew.verdict === 'undecided') return true;
  if (stPrev.modes.length === 0 || stNew.modes.length === 0) return true;
  const m1 = stPrev.modes[0], m2 = stNew.modes[0];
  const err = Math.max(m1.error, m2.error);

  // 1. Arclength linear extrapolation of Re λ from (prev2, prev) to next.
  if (stPrev2 && stPrev2.verdict !== 'undecided' && stPrev2.modes.length > 0 && dPrevArc > 1e-12) {
    const m0 = stPrev2.modes[0];
    const rho = dActual / dPrevArc;
    const rPred = m1.re + (m1.re - m0.re) * rho;
    // Predicted zero crossing that did not appear at m2 -> possible skipped narrow unstable window
    // or narrow double fold.
    if (m1.re < -kBar * m1.error && m2.re < -kBar * m2.error && rPred > kBar * err) return false;
    if (m1.re > kBar * m1.error && m2.re > kBar * m2.error && rPred < -kBar * err) return false;
    // Large deviation from predicted Re λ (captures strong curvature near folds / Hopf windows).
    const rScale = Math.max(Math.abs(m1.re), Math.abs(m2.re), 0.25 * Math.max(Math.abs(m1.im), Math.abs(m2.im)), 20 * err, 1e-6);
    if (Math.abs(m2.re - rPred) > 0.75 * rScale) return false;
  }

  // 2. Single-step mode jump / frequency continuity (protects the very first step too).
  const im2 = Math.abs(m2.im);
  if (im2 > 0) {
    const mc = findComplexPair(stPrev, im2, 0.9);
    if (mc) {
      const imc = Math.abs(mc.im);
      if (Math.abs(im2 - imc) > 0.35 * Math.max(im2, imc)) return false;
      const minR = Math.min(Math.abs(mc.re), Math.abs(m2.re));
      const maxR = Math.max(Math.abs(mc.re), Math.abs(m2.re));
      if ((mc.re < 0) === (m2.re < 0) && minR > 20 * err && maxR > 4 * minR + 0.3 * Math.max(im2, imc)) return false;
    } else if (im2 > 2 * Math.max(Math.abs(m1.re), Math.abs(m2.re), 20 * err)) {
      return false;
    }
  }
  const minR = Math.min(Math.abs(m1.re), Math.abs(m2.re));
  const maxR = Math.max(Math.abs(m1.re), Math.abs(m2.re));
  if ((m1.re < 0) === (m2.re < 0) && minR > 20 * err && maxR > 4 * minR + 0.5 * Math.max(Math.abs(m1.im), im2)) {
    return false;
  }
  return true;
}

// Chord-hyperplane solve at fraction theta in (0, 1) between two converged states A and B.
function solveOnChord(s, ctx, A, B, theta) {
  const { potIdx, potWt, wsEff, wx, wV, uPred, uLoPred, vPred } = ctx;
  const N = s.sys.size, fl = s.floating, K = fl.length;
  const chX = ctx.work.chordX;
  const chV = ctx.work.chordV;
  chX.fill(0);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i], o = s.fullOf[j];
    chX[j] = (B.snap.u[o] - A.snap.u[o]) + (B.snap.uLo[o] - A.snap.uLo[o]);
  }
  for (let a = 0; a < K; a++) chV[a] = B.snap.termV[fl[a]] - A.snap.termV[fl[a]];
  const chS = B.s - A.s;
  const nrm = tangentNorm(potIdx, potWt, wsEff, s.VT, chX, chV, chS);
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
  const ws = buildWeights(potIdx, potWt, wsEff, s.VT, chX, chV, chSu, wx, wV);
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

// Locate a fold (ts = 0) between A0 and B0 (tsA * tsB < 0) by safeguarded secant/bisection.
// Always returns { drive, terminals, s, state, located, note? }; never silently falls back.
function locateFold(s, ctx, A0, B0) {
  let A = A0, B = B0;
  let best = Math.abs(A.ts) <= Math.abs(B.ts) ? A : B;
  let converged = false;
  let solvedAny = false;
  let failNote = null;
  for (let it = 0; it < 24; it++) {
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
    if (Math.abs(M.ts) < 1e-8 || m.chordLen < 1e-9) {
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

// Find the complex mode (im > 0) in `st.modes` closest in frequency to `omegaEst`.
function findComplexPair(st, omegaEst, maxRelDist = 0.5) {
  let best = null, bestDist = Infinity;
  for (const m of st.modes) {
    if (!(m.im > 0)) continue;
    const dist = Math.abs(m.im - omegaEst);
    if (dist < bestDist) {
      bestDist = dist;
      best = m;
    }
  }
  if (!best) return null;
  if (omegaEst > 0 && bestDist > maxRelDist * Math.max(omegaEst, best.im)) return null;
  return best;
}

// Evaluate the tracked complex pair at state M during Hopf bisection.
function evalHopfMode(s, ctx, M, omegaEst) {
  const drive = ctx.p0 + M.s * ctx.dpSpan;
  s._restore(M.snap);
  const fHz = omegaEst / (2 * Math.PI);
  const baseOpts = ctx.stabOpts;
  const fastOpts = {
    ...baseOpts,
    shifts: [Math.max(1, 0.25 * omegaEst)],
    frequencies: baseOpts.frequencies ? [...baseOpts.frequencies, fHz] : [fHz],
  };
  let st = evaluateStability(s, ctx.termIdx, drive, fastOpts, ctx.isCurrentDrive, ctx.budget);
  let pair = findComplexPair(st, omegaEst);
  if (!pair) {
    st = evaluateStability(s, ctx.termIdx, drive, baseOpts, ctx.isCurrentDrive, ctx.budget);
    pair = findComplexPair(st, omegaEst);
  }
  return pair;
}

// Locate a Hopf bifurcation (Re λ = 0 for complex conjugate pair with Im λ > 0) between A0 and B0.
// Always reports `located: false` and `note` when location falls back instead of converging.
function locateHopf(s, ctx, A0, B0, pairA0, pairB0, omegaInit) {
  let A = A0, B = B0;
  let omega = omegaInit;
  let mA = pairA0 ?? evalHopfMode(s, ctx, A, omega);
  let mB = pairB0 ?? evalHopfMode(s, ctx, B, omega);
  if (!mA || !mB) return null;
  if ((mA.re > 0) === (mB.re > 0) && mA.re !== 0 && mB.re !== 0) return null;
  let best = Math.abs(mA.re) <= Math.abs(mB.re) ? { state: A, mode: mA } : { state: B, mode: mB };
  const kBar = ctx.stabOpts.kBar ?? 10;
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
    const mM = evalHopfMode(s, ctx, M, omega);
    if (!mM) {
      failNote = 'Hopf pair lost during bisection (using best bracket point)';
      break;
    }
    solvedAny = true;
    omega = mM.im;
    if (Math.abs(mM.re) < Math.abs(best.mode.re)) best = { state: M, mode: mM };
    const pA = ctx.p0 + A.s * ctx.dpSpan;
    const pB = ctx.p0 + B.s * ctx.dpSpan;
    const pM = ctx.p0 + M.s * ctx.dpSpan;
    if (Math.abs(mM.re) <= Math.max(kBar * mM.error, 1e-5 * omega) || Math.abs(pB - pA) <= 1e-6 * Math.max(1, Math.abs(pM))) {
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

// Detect whether any complex conjugate pair crossed Re λ = 0 between two stability results.
function detectHopfBracket(stA, stB, kBar) {
  const sig = (m) => (m.re > kBar * m.error ? 1 : m.re < -kBar * m.error ? -1 : 0);
  for (const ma of stA.modes) {
    if (!(ma.im > 0)) continue;
    const sa = sig(ma);
    if (sa === 0) continue;
    const mb = findComplexPair(stB, ma.im);
    if (!mb) continue;
    const sb = sig(mb);
    if (sb !== 0 && sa !== sb) return { pairA: ma, pairB: mb, omega: 0.5 * (ma.im + mb.im) };
  }
  if (stA.verdict === 'no growing mode found' && stB.verdict === 'unstable') {
    const grow = stB.modes.find((m) => m.re > kBar * m.error && m.im > 0);
    if (grow) return { pairA: findComplexPair(stA, grow.im), pairB: grow, omega: grow.im };
  }
  if (stA.verdict === 'unstable' && stB.verdict === 'no growing mode found') {
    const grow = stA.modes.find((m) => m.re > kBar * m.error && m.im > 0);
    if (grow) return { pairA: grow, pairB: findComplexPair(stB, grow.im), omega: grow.im };
  }
  return null;
}

// Count verified positive real modes (im === 0, re > kBar * error).
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
    } = opts;

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
    const term = this.terms[termIdx];
    const kind = term.drive.kind === 'I' ? 'I' : 'V';
    const isCurrentDrive = kind === 'I';
    const hasGeometry = this.model.geometry.type !== 'planar';
    const p0 = sourceAt(term.drive.src, this.time);
    const warnings = [];
    const nonAuto = this.terms.some((t) => Boolean(t.drive.src?.t));
    if (nonAuto) {
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
    const learnedKeep = { lostMode: this.lostMode, accepted: this.accepted, firstErrs: this.firstErrs?.slice() };

    const budget = Number.isFinite(budgetMs) ? branchBudget(budgetMs) : undefined;
    const points = [];
    const events = [];
    const stFullList = [];
    const states = [];
    const stepArcs = [0]; // arclength of step ending at points[k]
    let end = '';

    try {
      // Point 0 from current converged state.
      const pt0Data = readPointData(this, termIdx, p0, solutions);
      let st0 = null;
      if (doStab) {
        st0 = evaluateStability(this, termIdx, p0, stabOpts, isCurrentDrive, budget);
      }
      points.push({
        drive: p0,
        terminals: pt0Data.terminals,
        ...(doStab ? { stability: summarizeStability(st0) } : {}),
        ...(solutions ? { solution: pt0Data.solution } : {}),
      });
      stFullList.push(st0);

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

      // Prepare steady equations & laws for the branch.
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
        if (why) {
          return { terminal, kind, points, events, end: `refused: ${why}`, warnings };
        }
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
      const { potIdx, potWt } = stateMetricEntries(this);
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
      };
      const ctx = {
        termIdx,
        p0,
        dpSpan,
        potIdx,
        potWt,
        wsEff: WS_BASE,
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
      };

      const tx = new Float64Array(N + 1), tV = new Float64Array(K);
      const prevTx = new Float64Array(N + 1), prevTV = new Float64Array(K);
      const vStart = Float64Array.from(fl, (k) => this.termV[k]);

      // 1. Initial unit tangent at s = 0 and physical drive scaling so a far `to` never inflates
      //    the first step beyond the device's own response scale.
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

      // Compute physical drive scale at p0: Δp that moves the state by ~1 thermal/gate unit.
      if (Math.abs(tan0.ts) > 1e-12) {
        const statePerS = maxStateStep(this, tx, tV) / Math.abs(tan0.ts);
        if (statePerS > 1e-12) {
          const pScalePhys = absSpan / statePerS;
          const pScaleNorm = Math.min(absSpan, Math.max(1e-12, 10 * pScalePhys));
          if (pScaleNorm < absSpan) {
            ctx.wsEff = Math.min(WS_BASE * 1e6, WS_BASE * (absSpan / pScaleNorm));
            // Re-normalise initial tangent under the updated physical metric (reusing factor).
            tan0 = computeTangent(this, ctx, 0, null, null, null, tx, tV, true);
          }
        }
      }

      let ts = tan0.ts;
      let sCur = 0;

      states.push({
        snap: this._snapshot(),
        s: 0,
        tx: Float64Array.from(tx),
        tV: Float64Array.from(tV),
        ts,
      });

      // Step-size parameters. `step` is strictly enforced as a hard |Δp| cap; never loosen state
      // motion caps when `step` is provided. `wsEff` only scales the initial stride so a far `to`
      // never jumps over features at p0 while still allowing steps to grow smoothly in linear zones.
      const maxPotPerStep = 6; // thermal units (~0.155 V, or Δgate ≤ 0.6)
      const maxDVFloat = 0.25; // V per step on floating terminals
      const maxSPerStep = step !== undefined ? step / absSpan : 0.2;
      let dsArc = Math.abs(ts) > 1e-9 ? Math.min(maxSPerStep, (0.08 * WS_BASE) / ctx.wsEff) / Math.abs(ts) : 2;
      if (isCurrentDrive) {
        const flIdx = fl.indexOf(termIdx);
        if (flIdx >= 0 && Math.abs(tV[flIdx]) > 1e-9) dsArc = Math.min(dsArc, 0.1 / Math.abs(tV[flIdx]));
      }
      const initStateRate = maxStateStep(this, tx, tV);
      if (initStateRate > 1e-12) dsArc = Math.min(dsArc, 3 / initStateRate);
      dsArc = Math.max(MIN_ARC, Math.min(MAX_ARC, dsArc));

      let consecFails = 0;
      let consecModeFails = 0;
      let recentFail = 0;
      let stepCount = 0;

      // Record an accepted point (with optional precomputed stNew and precomputed foldInfo) and
      // locate any fold or Hopf on [kPrev, kNew]. If budget expires during event location, still
      // emit the event with `located: false` before unwinding.
      const acceptPoint = (sNew, newTx, newTV, newTs, dActual, preSt = null, preFold = null) => {
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
        stepArcs.push(dActual);

        const kNew = points.length - 1;

        // (a) Fold: dp/ds (ts) changes sign.
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

        // (b) Hopf: complex conjugate pair crosses Re λ = 0.
        if (doStab && stFullList[kPrev] && stNew) {
          const hopfBracket = detectHopfBracket(stFullList[kPrev], stNew, kBar);
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

        // Restore solver to curState so continuation proceeds cleanly.
        this._restore(curSnap);
        this.sourceOverride.set(termIdx, driveNew);
      };

      // Land cleanly on an exact span edge `sTarget` (0 or 1) from bracket `[A, B]`, enforcing
      // the SAME continuity guards (alignment, angle, step cap, state motion, voltage bound) as
      // any ordinary step.
      const landExact = (A, B, sTarget, refTx, refTV, refTs, enforceHardStep = true) => {
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
        if (!corLand.converged) return null;

        // Check voltage bound on the landed point.
        for (let a = 0; a < K; a++) {
          if (Math.abs(this.termV[fl[a]] - vStart[a]) > voltageBound) return null;
        }
        // Check hard `step` cap from the last accepted point (states[points.length - 1]).
        const lastPt = states[points.length - 1];
        const pLast = p0 + lastPt.s * dpSpan;
        const pLand = p0 + sTarget * dpSpan;
        if (enforceHardStep && step !== undefined && Math.abs(pLand - pLast) > step * (1 + 1e-6)) return null;

        // Continuity guards: forward alignment, corrector displacement, state step, and angle.
        const bracketLen = stateDistW(this, ctx, A.snap.u, A.snap.uLo, A.snap.termV, A.s, B.snap.u, B.snap.uLo, B.snap.termV, B.s);
        const dCor = stateDistW(this, ctx, ctx.uPred, ctx.uLoPred, predTermV, sTarget, this.u, this.uLo, this.termV, sTarget);
        const maxCorRatio = bracketLen > 1e-3 ? 0.75 : 2.0;
        if (bracketLen > 1e-12 && dCor > maxCorRatio * bracketLen) return null;
        const dActual = stateDistW(this, ctx, lastPt.snap.u, lastPt.snap.uLo, lastPt.snap.termV, lastPt.s, this.u, this.uLo, this.termV, sTarget);
        const stateMove = maxStateDiff(this, lastPt.snap, this.u, this.uLo, this.termV);
        if (stateMove > maxPotPerStep * 1.5 && bracketLen > MIN_ARC * 4) return null;

        const tanLand = computeTangent(this, ctx, sTarget, refTx, refTV, refTs, prevTx, prevTV, !!corLand.reusable);
        if (!tanLand.ok) return null;
        if (tanLand.cosTheta < COS_MIN && bracketLen > MIN_ARC * 4) return null;

        let align = ctx.wsEff * ctx.wsEff * refTs * (sTarget - A.s);
        const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
        for (let i = 0; i < potIdx.length; i++) {
          const j = potIdx[i], o = this.fullOf[j], w = potWt[i];
          align += invN * w * w * refTx[j] * ((this.u[o] - A.snap.u[o]) + (this.uLo[o] - A.snap.uLo[o]));
        }
        const invVT2 = 1 / (this.VT * this.VT);
        for (let a = 0; a < K; a++) align += invVT2 * refTV[a] * (this.termV[fl[a]] - A.snap.termV[fl[a]]);
        if (!(align > 0)) return null;

        return { s: sTarget, tx: Float64Array.from(prevTx), tV: Float64Array.from(prevTV), tan: tanLand, dActual };
      };

      // 2. Main pseudo-arclength continuation loop.
      while (stepCount < maxSteps) {
        budget?.check();

        // Cap stride by drive step, pointwise state step (potentials + 10× gates), and floating V.
        let dsStep = dsArc;
        if (Math.abs(ts) > 1e-9) dsStep = Math.min(dsStep, maxSPerStep / Math.abs(ts));
        const stateRate = maxStateStep(this, tx, tV);
        if (stateRate > 1e-12) dsStep = Math.min(dsStep, maxPotPerStep / stateRate);
        for (let a = 0; a < K; a++) {
          const vt = Math.abs(tV[a]);
          if (vt > 1e-9) dsStep = Math.min(dsStep, maxDVFloat / vt);
        }
        if (dsStep < MIN_ARC) dsStep = MIN_ARC;

        // If approaching s = 1 smoothly, cap predictor slightly past 1 so landing is tight.
        if (sCur < 1 && ts > 1e-9 && sCur + dsStep * ts > 1.02) {
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
        const predTermV = Float64Array.from(this.termV);
        if (predAlpha < 1) {
          this._restore(beforeSnap);
          addToState(this, tx, dsStep * predAlpha);
          for (let a = 0; a < K; a++) this.termV[fl[a]] += dsStep * predAlpha * tV[a];
        }
        const sInit = sCur + dsStep * predAlpha * ts;
        const ws = buildWeights(potIdx, potWt, ctx.wsEff, this.VT, tx, tV, ts, ctx.wx, ctx.wV);

        const cor = correctorNewton(this, ctx, sInit, ctx.uPred, ctx.uLoPred, ctx.vPred, sPred, ctx.wx, ctx.wV, ws);

        const rejectStep = (factor, reasonIfMin) => {
          this._restore(beforeSnap);
          consecFails++;
          recentFail = 3;
          dsArc = dsStep * factor;
          if (dsArc < MIN_ARC) {
            end = reasonIfMin;
            return true;
          }
          return false;
        };

        if (!cor.converged) {
          const stopped = rejectStep(
            consecFails > 0 ? 0.25 : 0.5,
            cor.modelError
              ? `continuation stopped: model refused (${cor.modelError})`
              : `continuation stopped: corrector failed at ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)} (${cor.error ?? 'step shrunk below minimum'})`,
          );
          if (stopped) break;
          continue;
        }

        const sNext = cor.s;

        // Guard 1: Corrector displacement, total arclength, pointwise state motion, and hard `step`.
        const dCor = stateDistW(this, ctx, ctx.uPred, ctx.uLoPred, predTermV, sPred, this.u, this.uLo, this.termV, sNext);
        const dActual = stateDistW(this, ctx, beforeSnap.u, beforeSnap.uLo, beforeSnap.termV, sCur, this.u, this.uLo, this.termV, sNext);
        const actualStateMove = maxStateDiff(this, beforeSnap, this.u, this.uLo, this.termV);
        const actualDp = Math.abs((sNext - sCur) * dpSpan);
        const stepExceeded = step !== undefined && actualDp > step * (1 + 1e-6);
        const dispTooLarge = (dCor > 0.5 * dsStep || dActual > 1.5 * dsStep || actualStateMove > maxPotPerStep * 1.25) && dsStep > MIN_ARC * 4;
        if (stepExceeded || dispTooLarge) {
          let factor = 0.5;
          if (stepExceeded && actualDp > 0) factor = Math.min(0.5, 0.8 * (step / actualDp));
          else if (dCor > 0.5 * dsStep) factor = Math.max(0.25, Math.min(0.5, 0.4 * dsStep / dCor));
          const stopped = rejectStep(
            factor,
            `continuation stopped: step control could not bound displacement near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
          );
          if (stopped) break;
          continue;
        }

        // Guard 2: Compute new unit tangent, check turning angle (cosTheta ≥ 0.75) and forward alignment.
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
        if (!tanNext.ok || (tanNext.cosTheta < COS_MIN && dsStep > MIN_ARC * 4)) {
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          const stopped = rejectStep(
            0.5,
            `continuation stopped: tangent turned sharply or became singular near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
          );
          if (stopped) break;
          continue;
        }

        let align = ctx.wsEff * ctx.wsEff * prevTs * (sNext - sCur);
        const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
        for (let i = 0; i < potIdx.length; i++) {
          const j = potIdx[i], o = this.fullOf[j], w = potWt[i];
          align += invN * w * w * prevTx[j] * ((this.u[o] - beforeSnap.u[o]) + (this.uLo[o] - beforeSnap.uLo[o]));
        }
        const invVT2 = 1 / (this.VT * this.VT);
        for (let a = 0; a < K; a++) align += invVT2 * prevTV[a] * (this.termV[fl[a]] - beforeSnap.termV[fl[a]]);
        if (!(align > 0)) {
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          const stopped = rejectStep(
            0.5,
            `continuation stopped: corrector jumped backward near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
          );
          if (stopped) break;
          continue;
        }

        // Guard 3: Check voltage bound on floating terminals BEFORE accepting a point outside the bound.
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

        // Check for fold on this step FIRST so the One Span Rule handles folds past `to` (Finding 2)
        // and immediate fold+return past `p0` (Finding 3) uniformly!
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
              // When stability is off, record the converged post-fold point and the unlocated fold
              // so `between: [kPrev, kNew]` references valid indices in `points`.
              const kPrev = points.length - 1;
              const candDrive = p0 + candState.s * dpSpan;
              this._restore(candState.snap);
              const ptData = readPointData(this, termIdx, candDrive, solutions);
              points.push({
                drive: candDrive,
                terminals: ptData.terminals,
                ...(solutions ? { solution: ptData.solution } : {}),
              });
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
          // Restore candState after fold location.
          this._restore(candState.snap);

          // Finding 2: if the fold itself lies past `to` (sFold >= 1) while sCur < 1, the curve
          // exited the span at `to` BEFORE reaching the fold! Land on s = 1 on the pre-fold leg.
          if (prevTs > 0 && sCur < 1 && preFold.s >= 1) {
            const landed = landExact(prevStateObj, preFold.state, 1, prevTx, prevTV, prevTs);
            if (landed) {
              // Run mode smoothness check on landed point if stability is on? Only if large step.
              acceptPoint(1, landed.tx, landed.tV, landed.tan.ts, landed.dActual);
              end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
              break;
            }
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            const stopped = rejectStep(0.5, `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`);
            if (stopped) break;
            continue;
          }
          // Symmetric edge: turning past s = 0 while heading downward from sCur > 0.
          if (prevTs < 0 && sCur > 0 && preFold.s <= 0) {
            const landed = landExact(prevStateObj, preFold.state, 0, prevTx, prevTV, prevTs);
            if (landed) {
              acceptPoint(0, landed.tx, landed.tV, landed.tan.ts, landed.dActual);
              const vNote = isCurrentDrive ? `, V = ${points[points.length - 1].terminals[terminal].V.toPrecision(5)} V` : '';
              end = `returned to starting drive (${formatDrive(p0, kind, hasGeometry)}${vNote}) after fold`;
              break;
            }
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            const stopped = rejectStep(0.5, `continuation stopped: could not land on starting drive ${formatDrive(p0, kind, hasGeometry)}`);
            if (stopped) break;
            continue;
          }
        }

        // One Span Rule — edge s = 1 (`to`): crossed from inside [0, 1] to >= 1.
        if (sCur < 1 && sNext >= 1) {
          const landed = landExact(prevStateObj, candState, 1, prevTx, prevTV, prevTs);
          if (landed) {
            let stLand = null;
            if (doStab) {
              budget?.check();
              stLand = evaluateStability(this, termIdx, to, stabOpts, isCurrentDrive, budget);
              const kPrev = points.length - 1;
              const smooth = modeStepSmooth(
                kPrev >= 1 ? stFullList[kPrev - 1] : null,
                stFullList[kPrev],
                stLand,
                stepArcs[kPrev],
                landed.dActual,
                kBar,
              );
              if (!smooth && dsStep > MIN_ARC * 4 && consecModeFails < 6) {
                consecModeFails++;
                tx.set(prevTx);
                tV.set(prevTV);
                ts = prevTs;
                const stopped = rejectStep(0.5, `continuation stopped: mode jump near ${formatDrive(to, kind, hasGeometry)}`);
                if (stopped) break;
                continue;
              }
            }
            consecModeFails = 0;
            acceptPoint(1, landed.tx, landed.tV, landed.tan.ts, landed.dActual, stLand);
            end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
            break;
          }
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          const stopped = rejectStep(0.25, `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`);
          if (stopped) break;
          continue;
        }

        // One Span Rule — edge s = 0 (`p0`): crossed back to <= 0 (including on the same step as
        // a fold, Finding 3!).
        if (sNext <= 0 && (sCur > 0 || (stepFold && preFold && preFold.s > 0))) {
          const bracketA = stepFold && preFold ? preFold.state : prevStateObj;
          const refT = stepFold && preFold ? preFold.state : prevStateObj;
          const landed = landExact(bracketA, candState, 0, refT.tx, refT.tV, refT.ts);
          if (landed) {
            acceptPoint(0, landed.tx, landed.tV, landed.tan.ts, landed.dActual, null, preFold);
            const vNote = isCurrentDrive ? `, V = ${points[points.length - 1].terminals[terminal].V.toPrecision(5)} V` : '';
            end = `returned to starting drive (${formatDrive(p0, kind, hasGeometry)}${vNote}) after fold`;
            break;
          }
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          const stopped = rejectStep(0.25, `returned past starting drive (${formatDrive(p0, kind, hasGeometry)}) after fold`);
          if (stopped) break;
          continue;
        }

        // Guard 4: Rightmost-mode smoothness check (when stability is enabled) BEFORE accepting!
        let stCand = null;
        if (doStab) {
          budget?.check();
          const driveCand = p0 + sNext * dpSpan;
          stCand = evaluateStability(this, termIdx, driveCand, stabOpts, isCurrentDrive, budget);
          const kPrev = points.length - 1;
          const smooth = modeStepSmooth(
            kPrev >= 1 ? stFullList[kPrev - 1] : null,
            stFullList[kPrev],
            stCand,
            stepArcs[kPrev],
            dActual,
            kBar,
          );
          if (!smooth && dsStep > MIN_ARC * 4 && consecModeFails < 6) {
            consecModeFails++;
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            const stopped = rejectStep(
              0.5,
              `continuation stopped: mode jumped sharply near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
            );
            if (stopped) break;
            continue;
          }
        }

        consecFails = 0;
        consecModeFails = 0;

        // Accept this point!
        ts = tanNext.ts;
        sCur = sNext;
        stepCount++;
        acceptPoint(sNext, tx, tV, ts, dActual, stCand, preFold);

        // Adapt arclength stride.
        if (recentFail > 0) {
          recentFail--;
          if (cor.iterations <= 8 && tanNext.cosTheta > 0.92 && dCor < 0.2 * dsStep) dsArc = Math.min(MAX_ARC, dsStep * 1.25);
          else dsArc = dsStep;
        } else if (cor.iterations <= 5 && tanNext.cosTheta > 0.95 && dCor < 0.15 * dsStep) {
          dsArc = Math.min(MAX_ARC, dsStep * 1.5);
        } else if (cor.iterations <= 8 && tanNext.cosTheta > 0.88 && dCor < 0.25 * dsStep) {
          dsArc = Math.min(MAX_ARC, dsStep * 1.2);
        } else if (cor.iterations >= 15 || tanNext.cosTheta < 0.82 || dCor > 0.35 * dsStep) {
          dsArc = Math.max(MIN_ARC, dsStep * 0.7);
        } else {
          dsArc = dsStep;
        }
      }

      if (!end && stepCount >= maxSteps) {
        end = `maxSteps reached (${maxSteps} steps, stopped at ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)})`;
      }

      // 3. Post-trace consistency audit between geometric events and stability marks.
      if (doStab && points.length >= 2) {
        const isNeutral = (v) => v === 'neutral within resolution';
        const isDecided = (v) => v === 'no growing mode found' || v === 'unstable';
        const geomEvents = events.filter((e) => e.kind === 'fold' || e.kind === 'hopf');

        // (a) Every fold must correspond to a stability change (or positive real mode count change),
        //     NEVER walking across another fold/Hopf event (prevents self-induced artefacts).
        for (const ev of geomEvents) {
          if (ev.kind !== 'fold') continue;
          const [kL, kR] = ev.between;
          let minL = 0;
          let maxR = points.length - 1;
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
          if (!stL || !stR) continue;
          if (!isDecided(stL.verdict) || !isDecided(stR.verdict)) continue;
          const flipped = stL.verdict !== stR.verdict;
          const realCountChanged = countPositiveReal(stL, kBar) !== countPositiveReal(stR, kBar);
          if (!flipped && !realCountChanged) {
            events.push({
              kind: 'stability mismatch',
              drive: ev.drive,
              terminals: ev.terminals,
              between: ev.between,
              note: `stability marks did not change across fold (${stL.verdict} on both sides)`,
            });
          }
        }

        // (b) Every flip between decided non-neutral states must be explained by a fold or Hopf on
        //     the spanning interval.
        let prevDec = -1;
        for (let i = 0; i < points.length; i++) {
          const v = points[i].stability.verdict;
          if (!isDecided(v)) continue;
          if (prevDec >= 0 && points[prevDec].stability.verdict !== v) {
            const covered = events.some(
              (ev) => (ev.kind === 'fold' || ev.kind === 'hopf') && ev.between[0] >= prevDec && ev.between[1] <= i,
            );
            if (!covered) {
              let resolved = false;
              const hopfBracket = detectHopfBracket(stFullList[prevDec], stFullList[i], kBar);
              if (hopfBracket && states[prevDec] && states[i]) {
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
      }

      return { terminal, kind, points, events, end, warnings };
    } catch (err) {
      if (err === BUDGET_SPENT) {
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
      this.lostMode = learnedKeep.lostMode;
      this.accepted = learnedKeep.accepted;
      this.firstErrs = learnedKeep.firstErrs;
      if (this.drivenKey !== origDrivenKey) this._lawsNow();
      this._restore(snap0);
      this._refreshSources();
    }
  }
}
