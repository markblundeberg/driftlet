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
//   - Heuristics guarding step continuity (applied to ordinary steps, landing steps, and
//     natural-parameter fallback steps alike; never switched off at the minimum stride):
//       1. Corrector displacement: `||z_next − z_pred||_W ≤ 0.5 · dsStep` and
//          `||z_next − z_cur||_W ≤ 1.5 · dsStep`, plus pointwise state motion bounded.
//       2. Tangent turning angle: `cos θ ≥ 0.75` (≤ 41°) and positive forward alignment.
//       3. Hard `step` cap: when `opts.step` is given, `|p_next − p_cur| ≤ step` on every step.
//       4. All-mode tracking (when `stability !== false`): every listed mode (complex pairs and
//          real modes alike, not only the rightmost) is matched from point to point by eigenvector
//          overlap. Steps are halved whenever matching is ambiguous, whenever any matched mode's
//          Re λ (or a complex pair's Im λ) moves by more than a fraction of its own scale (with
//          smooth Boltzmann exponential decay of strictly negative real modes recognised in log
//          space), or whenever arclength extrapolation predicts a missed crossing or sharp bend.
//   - Natural-parameter fallback: when pseudo-arclength stalls (corrector fails, guards shrink the
//     stride below minimum or exhaust halvings, or a rank-deficient Jacobian produces a degenerate
//     tangent), branch() tries a plain steady solve at the next drive value from the last accepted
//     state and accepts it only if it passes the same continuity and mode-tracking guards; failing
//     that, it stops and states why in `end`.
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
//   - What a mark does NOT promise: a stability mark holds at the sampled point, and all-mode
//     tracking catches narrow windows that any listed mode approaches, but an unstable window
//     narrower than the step that is unseen by all listed modes at both endpoints can still be
//     missed. (In the midpoint audit across random multi-feature devices and 72 downward HH sweeps
//     including narrowed g.Na = 850, zero windows were missed.)
//   - Fold (`dp/ds` changes sign): located by safeguarded secant/bisection along the exact arc to
//     `|ts| < 1e-8` (or chord < 1e-9), and confirmed geometrically as a genuine drive extremum
//     with a real change in `p` beyond tangent noise, reporting `{ kind: 'fold', drive, terminals,
//     between }`.
//   - Hopf (`Re λ` of a complex conjugate pair `Im λ > 0` crosses 0): matched by eigenvector
//     overlap (not a narrow frequency window, so large frequency shifts across a step are still
//     matched) and located by safeguarded secant/bisection along the exact arc to
//     `|Re λ| ≤ max(kBar·error, 1e-5·ω)` or `|Δp| ≤ 1e-6·max(1, |p|)`, reporting
//     `{ kind: 'hopf', drive, terminals, omega, between }`.
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
    // Compatible with steady.js's _steadySteps({ budget }) interface.
    record() {
      attempted = true;
    },
    remaining() {
      return Math.max(0, ms - (now() - start));
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

function kirchhoffSum(s) {
  let sum = 0;
  for (const I of s.termI) sum += I;
  return Math.abs(sum);
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
    const okKirchhoff = () => !s._missesKirchhoff() || (ctx.startMissesKirchhoff && kirchhoffSum(s) <= 10 * ctx.startKirchhoffSum + 1e-13);
    if (r2.converged && okKirchhoff()) {
      return {
        converged: true,
        iterations: Math.max(r1.iterations, 6),
        s: sVal,
        reusable: false,
        minSchurPivot: r2.minSchurPivot,
      };
    }
    // If direct islands round failed or still missed Kirchhoff, try one budgeted _steadySteps
    // polish (which includes GMRES refinement) at fixed sVal.
    s._restore(preIslSnap);
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    const qPol = s._steadySteps(budget ? { maxSteps: 1, budget } : { maxSteps: 1 });
    if (qPol.converged && okKirchhoff()) {
      return {
        converged: true,
        iterations: Math.max(r1.iterations, 6),
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

// Unique modes (im >= 0): each real mode once, each complex conjugate pair once (with im > 0).
function uniqueModes(st) {
  if (!st || !st.modes) return [];
  const out = [];
  for (const m of st.modes) {
    if (m.im >= 0) out.push(m);
  }
  return out;
}

// Complex eigenvector overlap |v_aᴴ v_b| in [0, 1] (normalized coordScale vectors), or -1 if absent.
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

// Biorthogonal mode overlap between two modes using the dual vectors u = Mᴴ w / |wᴴ M v| (which
// satisfy |u_iᴴ v_j| = δ_{ij} even on strongly non-normal pencils, independent of M's units),
// falling back to right-vector cosine |v_aᴴ v_b| only when duals are absent or non-finite.
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
    if (Number.isFinite(ov)) {
      const clamped = Math.min(1, ov);
      if (clamped >= 0.5 || !(a._cond > 1e6 || b._cond > 1e6)) return clamped;
      // When a mode has a hugely ill-conditioned left vector (e.g. depleted floating island),
      // fall back to right-vector cosine if eigenvalues match tightly.
      const sc = Math.max(Math.hypot(a.re, a.im), Math.hypot(b.re, b.im), 1e-12);
      if (Math.hypot(a.re - b.re, a.im - b.im) <= 0.1 * sc) {
        return Math.max(clamped, vecOverlap(a, b));
      }
      return clamped;
    }
  }
  return vecOverlap(a, b);
}

// Match unique modes between two stability results 1-to-1 by biorthogonal eigenvector overlap
// (Stage 1), λ-proximity / Boltzmann depletion for localized moving-edge real modes (Stage 2),
// and weak eigenvector overlap (Stage 3). Returns { A, B, pairs, matchedA, matchedB, O }.
function matchModes(stA, stB, actualStateMove = 0) {
  const A = uniqueModes(stA);
  const B = uniqueModes(stB);
  const nA = A.length, nB = B.length;
  const O = Array.from({ length: nA }, () => new Float64Array(nB));
  const cands = [];
  for (let i = 0; i < nA; i++) {
    for (let j = 0; j < nB; j++) {
      let s = biorthOverlap(A[i], B[j]);
      if (s < 0) {
        const scale = Math.max(Math.hypot(A[i].re, A[i].im), Math.hypot(B[j].re, B[j].im), 1e-12);
        const dist = Math.hypot(A[i].re - B[j].re, A[i].im - B[j].im);
        s = 1 / (1 + dist / scale);
      }
      O[i][j] = s;
      cands.push({ i, j, s });
    }
  }
  cands.sort((x, y) => y.s - x.s);
  const matchedA = new Int32Array(nA).fill(-1);
  const matchedB = new Int32Array(nB).fill(-1);
  const pairs = [];
  // Stage 1: 1-to-1 matching by biorthogonal eigenvector overlap.
  for (const c of cands) {
    if (c.s < 0.45) break;
    if (matchedA[c.i] >= 0 || matchedB[c.j] >= 0) continue;
    matchedA[c.i] = c.j;
    matchedB[c.j] = c.i;
    // Compute runner-up overlap against non-nearly-degenerate modes on either side.
    const ai = A[c.i], bj = B[c.j];
    let runnerUp = 0;
    for (let jj = 0; jj < nB; jj++) {
      if (jj === c.j) continue;
      const sc = Math.max(Math.hypot(bj.re, bj.im), Math.hypot(B[jj].re, B[jj].im), 1e-12);
      if (Math.hypot(bj.re - B[jj].re, bj.im - B[jj].im) > 0.1 * sc && O[c.i][jj] > runnerUp) {
        runnerUp = O[c.i][jj];
      }
    }
    for (let ii = 0; ii < nA; ii++) {
      if (ii === c.i) continue;
      const sc = Math.max(Math.hypot(ai.re, ai.im), Math.hypot(A[ii].re, A[ii].im), 1e-12);
      if (Math.hypot(ai.re - A[ii].re, ai.im - A[ii].im) > 0.1 * sc && O[ii][c.j] > runnerUp) {
        runnerUp = O[ii][c.j];
      }
    }
    pairs.push({ i: c.i, j: c.j, a: ai, b: bj, overlap: c.s, runnerUp, byEig: false });
  }

  // Stage 2: λ-proximity fallback for remaining unmatched same-type, same-sign modes (handles
  // localized depletion-edge / trap modes whose eigenvector shifts by 1 grid cell so overlap drops
  // to ~0 while λ barely moves, or leading strictly negative real modes following smooth Boltzmann
  // depletion across moving depletion edges). Never matches across Re λ = 0 or real <-> oscillatory.
  const eigCands = [];
  for (let i = 0; i < nA; i++) {
    if (matchedA[i] >= 0) continue;
    for (let j = 0; j < nB; j++) {
      if (matchedB[j] >= 0) continue;
      const ai = A[i], bj = B[j];
      if ((ai.re < 0) !== (bj.re < 0)) continue;
      if (isOscillatory(ai) !== isOscillatory(bj)) continue;
      const err = Math.max(ai.error, bj.error);
      const sc = Math.max(Math.hypot(ai.re, ai.im), Math.hypot(bj.re, bj.im), 20 * err, 1e-12);
      const rel = Math.hypot(bj.re - ai.re, bj.im - ai.im) / sc;
      if (rel <= 0.35) {
        eigCands.push({ i, j, rel });
      } else if (i === 0 && j === 0 && !isOscillatory(ai) && ai.re < -20 * err && bj.re < -20 * err && actualStateMove >= 0.5) {
        const dLn = Math.abs(Math.log(-bj.re / -ai.re));
        if (dLn <= 3.2 * actualStateMove + 0.5) eigCands.push({ i, j, rel: 0.35 + 0.01 * dLn });
      }
    }
  }
  eigCands.sort((x, y) => x.rel - y.rel);
  for (const c of eigCands) {
    if (matchedA[c.i] >= 0 || matchedB[c.j] >= 0) continue;
    matchedA[c.i] = c.j;
    matchedB[c.j] = c.i;
    pairs.push({ i: c.i, j: c.j, a: A[c.i], b: B[c.j], overlap: O[c.i][c.j], runnerUp: 0, byEig: true });
  }

  // Stage 3: Weak eigenvector matches (0.30..0.45) for any still-unmatched pair (e.g. a mode that
  // moved rapidly across a coarse step so Check 2 can inspect and halve it).
  for (const c of cands) {
    if (c.s < 0.3) break;
    if (matchedA[c.i] >= 0 || matchedB[c.j] >= 0) continue;
    matchedA[c.i] = c.j;
    matchedB[c.j] = c.i;
    pairs.push({ i: c.i, j: c.j, a: A[c.i], b: B[c.j], overlap: c.s, runnerUp: 0, byEig: false });
  }
  return { A, B, pairs, matchedA, matchedB, O };
}

// Check if an unmatched real mode `r` on one side is the coalescing partner of a matched complex
// pair `c` on the other side (2 real modes <-> 1 complex pair across an exceptional point).
function isExceptionalPartner(r, otherList, otherMatched) {
  if (r.im > 0) return false;
  for (let j = 0; j < otherList.length; j++) {
    const c = otherList[j];
    if (!(c.im > 0) || otherMatched[j] < 0) continue;
    const ovB = biorthOverlap(r, c);
    const ovV = vecOverlap(r, c);
    if (Math.max(ovB, ovV) >= 0.45) return true;
    const sc = Math.max(Math.abs(r.re), Math.abs(c.re), 1e-12);
    if (c.im <= 0.55 * sc && Math.abs(r.re - c.re) <= 0.55 * sc) return true;
  }
  return false;
}

// True for a genuinely oscillatory complex pair (not a deep-LHP real mode with small numerical Im).
function isOscillatory(m) {
  return m.im > 0 && (m.re >= 0 || m.im > 0.25 * Math.abs(m.re));
}

// All-mode tracking & smoothness check between consecutive states (prev -> next, with optional
// prev2 for arclength extrapolation). Tracks EVERY listed mode (complex pairs and real modes
// alike) by biorthogonal eigenvector overlap — including when a point's verdict is 'undecided'
// so long as it has an oscillatory or growing mode! Returns { ok, reason }.
function modeStepSmooth(stPrev2, stPrev, stNew, dPrevArc, dActual, actualStateMove, kBar, nWant, stats) {
  if (!stPrev || !stNew) return { ok: true };
  const m12 = matchModes(stPrev, stNew, actualStateMove);
  const { A, B, pairs, matchedA, matchedB } = m12;
  if (A.length === 0 || B.length === 0) return { ok: true };
  // If either endpoint is 'undecided' and neither side has any oscillatory or growing mode (e.g.
  // broken Arnoldi shift grid on a depleted floating island), do not stall on LHP real jitter.
  // (When an oscillatory or growing mode IS present — e.g. HH right beside a Hopf — still track it!)
  if (stPrev.verdict === 'undecided' || stNew.verdict === 'undecided') {
    const hasKey = (arr) => arr.some((m) => m.re >= 0 || isOscillatory(m));
    if (!hasKey(A) && !hasKey(B)) return { ok: true };
  }

  // Modes strictly inside the truncation window (excluding modes near the deepest cutoff tail when
  // the list hit nWant slots). Leading mode, any growing mode, and any non-tail oscillatory pair
  // are always active.
  const cutA = stPrev.modes.length >= nWant && A.length > 3 ? A[A.length - 1].re : -Infinity;
  const cutB = stNew.modes.length >= nWant && B.length > 3 ? B[B.length - 1].re : -Infinity;
  const leadA = A[0].re;
  const leadB = B[0].re;
  const isInteriorA = (i) => {
    const m = A[i];
    if (i === 0 || m.re >= 0 || A.length <= 2) return true;
    if (isOscillatory(m)) return !(stPrev.modes.length >= nWant && i === A.length - 1);
    if (stPrev.modes.length >= nWant && (i >= A.length - 2 || m.re <= 0.7 * cutA)) return false;
    return true;
  };
  const isInteriorB = (j) => {
    const m = B[j];
    if (j === 0 || m.re >= 0 || B.length <= 2) return true;
    if (isOscillatory(m)) return !(stNew.modes.length >= nWant && j === B.length - 1);
    if (stNew.modes.length >= nWant && (j >= B.length - 2 || m.re <= 0.7 * cutB)) return false;
    return true;
  };

  // 1. Unmatched check for critical modes: leading mode (index 0), any growing mode, and any
  //    interior oscillatory complex pair (secondary deep-LHP real modes may enter/exit Arnoldi's
  //    shift-invert subset without moving near 0).
  for (let i = 0; i < A.length; i++) {
    if (!isInteriorA(i) || matchedA[i] >= 0) continue;
    const m = A[i];
    if (isExceptionalPartner(m, B, matchedB)) continue;
    // Strictly negative real leading mode entering/leaving by smooth Boltzmann depletion from the
    // bulk continuum (or Arnoldi LHP shift flicker on a small step) while the adjacent bulk mode
    // matches cleanly on both decaying endpoints.
    if (i === 0 && !isOscillatory(m) && m.re < -kBar * m.error && matchedB[0] >= 0 && stPrev.verdict === 'no growing mode found' && stNew.verdict === 'no growing mode found') {
      const b0 = B[0];
      if (actualStateMove < 1.0 || (b0.re < -20 * b0.error && Math.abs(Math.log(-m.re / -b0.re)) <= 3.2 * actualStateMove + 0.5)) continue;
    }
    if (i > 0 && !isOscillatory(m) && m.re < 0 && (leadA < 0 ? m.re < 1.5 * leadA : false)) continue;
    if (i > 0 && !isOscillatory(m) && m.re < 0 && (m._cond > 1e3 || stPrev.modes.length >= nWant - 1)) continue;
    if (stats) stats.modeMatchAmbiguous++;
    return { ok: false, reason: 'unmatched mode' };
  }
  for (let j = 0; j < B.length; j++) {
    if (!isInteriorB(j) || matchedB[j] >= 0) continue;
    const m = B[j];
    if (isExceptionalPartner(m, A, matchedA)) continue;
    if (j === 0 && !isOscillatory(m) && m.re < -kBar * m.error && matchedA[0] >= 0 && stPrev.verdict === 'no growing mode found' && stNew.verdict === 'no growing mode found') {
      const a0 = A[0];
      if (actualStateMove < 1.0 || (a0.re < -20 * a0.error && Math.abs(Math.log(-m.re / -a0.re)) <= 3.2 * actualStateMove + 0.5)) continue;
    }
    if (j > 0 && !isOscillatory(m) && m.re < 0 && (leadB < 0 ? m.re < 1.5 * leadB : false)) continue;
    if (j > 0 && !isOscillatory(m) && m.re < 0 && (m._cond > 1e3 || stNew.modes.length >= nWant - 1)) continue;
    if (stats) stats.modeMatchAmbiguous++;
    return { ok: false, reason: 'new unmatched mode' };
  }

  // Optional predecessor match (prev2 -> prev) for per-mode arclength extrapolation.
  const m01 = stPrev2 && stPrev2.modes.length > 0 && dPrevArc > 1e-12 ? matchModes(stPrev2, stPrev) : null;
  const prev2OfA = new Map();
  if (m01) {
    for (const p of m01.pairs) {
      if (p.overlap >= 0.6 || p.byEig) prev2OfA.set(p.j, p.a);
    }
  }

  // 2. Per-mode match quality, Re λ scale movement, Im λ movement, and arclength extrapolation
  //    across ALL matched interior modes (complex pairs and real modes alike).
  const rho = dPrevArc > 1e-12 ? dActual / dPrevArc : 0;
  for (const p of pairs) {
    const { i, j, a, b, overlap, runnerUp, byEig } = p;
    const active = isInteriorA(i) || isInteriorB(j);
    if (!active) continue;

    const err = Math.max(a.error, b.error);
    const sameSign = (a.re < 0) === (b.re < 0);
    const maxR = Math.max(Math.abs(a.re), Math.abs(b.re));
    const oscA = isOscillatory(a), oscB = isOscillatory(b);
    const isComplex = oscA || oscB;
    const isLeading = i === 0 || j === 0;
    const realCoal = oscA !== oscB;

    // Eigenvector match ambiguity: strict on leading modes & oscillatory complex pairs unless
    // matched cleanly by Stage-2 λ-proximity (localized moving-edge modes whose eigenvector shifted
    // by a grid cell while λ barely moved).
    const eigTight = Math.hypot(b.re - a.re, b.im - a.im) <= 0.15 * Math.max(Math.hypot(a.re, a.im), Math.hypot(b.re, b.im), 20 * err);
    const deepRealTight = !isComplex && !isLeading && sameSign && a.re < 0 && Math.abs(b.re - a.re) <= 0.25 * Math.max(maxR, 20 * err);
    const minOv = isLeading || (isComplex && !realCoal) ? 0.6 : 0.5;
    if (!byEig && !deepRealTight && !(eigTight && overlap >= 0.5)) {
      if (overlap < minOv || (!realCoal && runnerUp > 0.88 * overlap && overlap < 0.90)) {
        if (stats) stats.modeMatchAmbiguous++;
        return { ok: false, reason: 'ambiguous mode match' };
      }
    }

    const a0 = prev2OfA.get(i) ?? null;

    if (isComplex) {
      // Never cross BOTH an exceptional point (real <-> complex) AND Re λ = 0 in a single step:
      // halving separates them so both sides of a Hopf hold the complex pair.
      if (realCoal && !sameSign) {
        if (stats) stats.modeImMove++;
        return { ok: false, reason: 'real-complex split crossed Re λ = 0' };
      }
      // Frequency continuity (scaled by max(a.im, b.im, maxR) so emerging exceptional points deep
      // in the LHP where Im λ ~ sqrt(Δs) << |Re λ| do not halve forever, while near-axis oscillatory
      // pairs like HH still guard large frequency shifts).
      if (oscA && oscB) {
        const imScale = Math.max(a.im, b.im, maxR);
        if (Math.abs(b.im - a.im) > 0.45 * imScale) {
          if (stats) stats.modeImMove++;
          return { ok: false, reason: 'complex frequency jumped' };
        }
      } else {
        const imC = Math.max(a.im, b.im);
        const deepBehind = a.re < 0 && b.re < 0 && leadA < 0 && leadB < 0 && maxR > 10 * Math.max(Math.abs(leadA), Math.abs(leadB));
        if (!deepBehind && imC > 0.55 * Math.max(maxR, 20 * err)) {
          if (stats) stats.modeImMove++;
          return { ok: false, reason: 'unresolved real-complex split' };
        }
      }
      if (sameSign) {
        const maxIm = Math.max(a.im, b.im);
        // Scale for Re λ movement: include a small Im floor so a point landing right beside a Hopf
        // (Re λ ≈ ±0.1 with ω ≈ 1000) does not stall with |ΔRe|/maxR ≈ 1, while keeping the stable-side
        // floor tiny (0.5% of ω) so narrow unstable windows (like HH g.Na = 850) are never skipped.
        const reScale = a.re > 0
          ? Math.max(maxR, 0.05 * maxIm, 20 * err)
          : Math.max(maxR, 0.005 * maxIm, 20 * err);
        if (maxR > 20 * err && Math.abs(b.re - a.re) > 0.55 * reScale) {
          if (stats) {
            stats.modeReMove++;
            if (i > 0 || j > 0) stats.modeReMoveSecondary++;
            else stats.modeReMoveLeading++;
          }
          return { ok: false, reason: 'complex Re λ moved past scale' };
        }
        if (a0) {
          const rPred = a.re + (a.re - a0.re) * rho;
          if ((a.re < -kBar * a.error && b.re < -kBar * b.error && rPred > kBar * err) ||
              (a.re > kBar * a.error && b.re > kBar * b.error && rPred < -kBar * err)) {
            if (stats) stats.modeExtrap++;
            return { ok: false, reason: 'extrapolated complex zero crossing' };
          }
          if (Math.abs(b.re - rPred) > 0.6 * reScale) {
            if (stats) stats.modeExtrap++;
            return { ok: false, reason: 'complex Re λ curved sharply' };
          }
        }
      }
    } else if (sameSign) {
      // Real mode on the same side of 0: recognise smooth Boltzmann exponential depletion
      // λ(s) ~ -A exp(-k u(s)) (covering |z| = 1..3 ions and multi-carrier recombination) in log
      // space so reverse-bias depletion sweeps do not crawl, while catching genuine linear /
      // non-exponential approaches to zero or fast real-mode moves.
      let boltzSmooth = false;
      if (a.re < -20 * err && b.re < -20 * err) {
        const dLn = Math.log(-b.re / -a.re);
        if (a0 && a0.re < -20 * a0.error) {
          const lnPred = Math.log(-a.re) + Math.log(-a.re / -a0.re) * rho;
          if (Math.abs(Math.log(-b.re) - lnPred) <= 0.75 && Math.abs(dLn) <= 3.2 * actualStateMove + 0.5) {
            boltzSmooth = true;
          }
        } else if (actualStateMove >= 1.0 && Math.abs(dLn) <= 3.0 * actualStateMove + 0.5) {
          boltzSmooth = true;
        }
      }
      if (boltzSmooth) {
        if (stats) stats.modeBoltzExempt++;
        continue;
      }
      if (maxR > 20 * err && Math.abs(b.re - a.re) > 0.65 * maxR) {
        if (stats) {
          stats.modeReMove++;
          if (i > 0 || j > 0) stats.modeReMoveSecondary++;
          else stats.modeReMoveLeading++;
        }
        return { ok: false, reason: 'real Re λ moved past scale' };
      }
      if (a0) {
        // For strictly negative real modes already decaying exponentially (|a| < 0.5 |a0|), linear
        // extrapolation always overshoots 0; check log-space consistency instead (using kBar * error
        // so rapidly depleting modes approaching resolution do not trigger false zero-crossing stops).
        const expDecay = a0.re < -kBar * a0.error && a.re < -kBar * a.error && b.re < 0 && Math.abs(a.re) < 0.5 * Math.abs(a0.re);
        if (expDecay) {
          if (b.re < -kBar * b.error) {
            const lnPred = Math.log(-a.re) + Math.log(-a.re / -a0.re) * rho;
            if (Math.abs(Math.log(-b.re) - lnPred) > 1.5) {
              if (stats) stats.modeExtrap++;
              return { ok: false, reason: 'exponential real decay bent sharply' };
            }
          }
        } else {
          const rPred = a.re + (a.re - a0.re) * rho;
          if ((a.re < -kBar * a.error && b.re < -kBar * b.error && rPred > kBar * err) ||
              (a.re > kBar * a.error && b.re > kBar * b.error && rPred < -kBar * err)) {
            if (stats) stats.modeExtrap++;
            return { ok: false, reason: 'extrapolated real zero crossing' };
          }
          if (Math.abs(b.re - rPred) > 0.75 * Math.max(maxR, 20 * err)) {
            if (stats) stats.modeExtrap++;
            return { ok: false, reason: 'real Re λ curved sharply' };
          }
        }
      }
    }
  }
  return { ok: true };
}

// Verify that a located fold is a genuine drive extremum with a real change in p beyond tangent
// noise (not a numerical sign wobble on a monotonically advancing curve).
function isConfirmedFold(A, B, fold, absSpan) {
  if (!fold || !fold.located) return false;
  const sA = A.s, sB = B.s, sF = fold.s;
  const sTol = 1e-8 * Math.max(1, Math.abs(sF));
  // Smooth turn: ts > 0 -> ts < 0 implies local maximum sF >= max(sA, sB); opposite implies local minimum.
  if (A.ts > 0 && B.ts < 0 && sF < Math.max(sA, sB) - sTol) return false;
  if (A.ts < 0 && B.ts > 0 && sF > Math.min(sA, sB) + sTol) return false;
  // Require real drive turn beyond tangent/round-off noise.
  const dpTurn = Math.max(Math.abs(sF - sA), Math.abs(sF - sB)) * absSpan;
  const pScale = Math.max(1, Math.abs(fold.drive), absSpan);
  if (!(dpTurn > 1e-8 * pScale)) return false;
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

// Find the complex mode (im > 0) in `st.modes` best matching `refMode` by eigenvector overlap,
// falling back to closest frequency when eigenvectors are absent. Never imposes a narrow 50% ω cap.
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

// Evaluate the tracked complex pair at state M during Hopf bisection.
function evalHopfMode(s, ctx, M, omegaEst, refMode = null, extraOmega = 0) {
  const drive = ctx.p0 + M.s * ctx.dpSpan;
  s._restore(M.snap);
  const fHz = omegaEst / (2 * Math.PI);
  const baseOpts = ctx.stabOpts;
  const freqs = baseOpts.frequencies ? [...baseOpts.frequencies, fHz] : [fHz];
  if (extraOmega > 0 && Math.abs(extraOmega - omegaEst) > 0.1 * Math.max(extraOmega, omegaEst)) {
    freqs.push(extraOmega / (2 * Math.PI));
  }
  const fastOpts = {
    ...baseOpts,
    shifts: [Math.max(1, 0.25 * omegaEst)],
    frequencies: freqs,
  };
  let st = evaluateStability(s, ctx.termIdx, drive, fastOpts, ctx.isCurrentDrive, ctx.budget);
  let pair = findComplexPair(st, omegaEst, refMode);
  if (!pair) {
    st = evaluateStability(s, ctx.termIdx, drive, baseOpts, ctx.isCurrentDrive, ctx.budget);
    pair = findComplexPair(st, omegaEst, refMode);
  }
  return pair;
}

// Locate a Hopf bifurcation (Re λ = 0 for complex conjugate pair with Im λ > 0) between A0 and B0.
// Always reports `located: false` and `note` when location falls back instead of converging.
function locateHopf(s, ctx, A0, B0, pairA0, pairB0, omegaInit) {
  let A = A0, B = B0;
  let omega = omegaInit;
  let mA = pairA0 ?? evalHopfMode(s, ctx, A, omega, pairB0);
  let mB = pairB0 ?? evalHopfMode(s, ctx, B, omega, pairA0);
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

// Detect whether any complex conjugate pair crossed Re λ = 0 between two stability results using
// eigenvector-overlap mode matching (no narrow frequency cutoff).
function detectHopfBracket(stA, stB, kBar, stats = null) {
  const sig = (m) => (m.re > kBar * m.error ? 1 : m.re < -kBar * m.error ? -1 : 0);
  const matched = matchModes(stA, stB);
  for (const p of matched.pairs) {
    const { a: ma, b: mb } = p;
    if (!(ma.im > 0 && mb.im > 0)) continue;
    const sa = sig(ma), sb = sig(mb);
    if (sa !== 0 && sb !== 0 && sa !== sb) {
      if (stats) {
        stats.hopfByVecMatch++;
        if (Math.abs(ma.im - mb.im) > 0.5 * Math.max(ma.im, mb.im)) stats.hopfLargeFreqShift++;
      }
      return { pairA: ma, pairB: mb, omega: 0.5 * (ma.im + mb.im) };
    }
  }
  if (stA.verdict === 'no growing mode found' && stB.verdict === 'unstable') {
    const grow = stB.modes.find((m) => m.re > kBar * m.error && m.im > 0);
    if (grow) {
      const pairA = findComplexPair(stA, grow.im, grow);
      if (stats) stats.hopfByVecMatch++;
      return { pairA, pairB: grow, omega: grow.im };
    }
  }
  if (stA.verdict === 'unstable' && stB.verdict === 'no growing mode found') {
    const grow = stA.modes.find((m) => m.re > kBar * m.error && m.im > 0);
    if (grow) {
      const pairB = findComplexPair(stB, grow.im, grow);
      if (stats) stats.hopfByVecMatch++;
      return { pairA: grow, pairB, omega: grow.im };
    }
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
      _stats = null,
    } = opts;
    if (_stats) {
      Object.assign(_stats, {
        modeMatchAmbiguous: 0,
        modeReMove: 0,
        modeReMoveLeading: 0,
        modeReMoveSecondary: 0,
        modeImMove: 0,
        modeExtrap: 0,
        modeBoltzExempt: 0,
        falseFoldRejected: 0,
        guardAtMinArc: 0,
        naturalFallbackTried: 0,
        naturalFallbackAccepted: 0,
        naturalFallbackRejected: 0,
        hopfByVecMatch: 0,
        hopfLargeFreqShift: 0,
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
    const wasSolvedV = this.solvedV ? this.solvedV.slice() : this.solvedV;
    const wasStepped = this.stepped;
    const wasColdStart = this.coldStart;
    const learnedKeep = { lostMode: this.lostMode, accepted: this.accepted, firstErrs: this.firstErrs?.slice() };

    const budget = Number.isFinite(budgetMs) ? branchBudget(budgetMs) : undefined;
    const points = [];
    const events = [];
    const stFullList = [];
    const states = [];
    const stepArcs = [0]; // arclength of step ending at points[k]
    let end = '';
    let ctx = null;

    // Post-trace consistency audit between geometric events and stability marks. Runs on both
    // normal completion (`allowLocate = true`) and budget stops (`allowLocate = false`, no extra
    // stability evaluations).
    const runAudit = (allowLocate) => {
      if (!doStab || points.length < 2) return;
      const isNeutral = (v) => v === 'neutral within resolution';
      const isDecided = (v) => v === 'no growing mode found' || v === 'unstable';
      const geomEvents = events.filter((e) => e.kind === 'fold' || e.kind === 'hopf');

      // (a) Every fold must correspond to a stability change (or positive real mode count change),
      //     NEVER walking across another fold/Hopf event (prevents self-induced artefacts).
      for (const ev of geomEvents) {
        if (ev.kind !== 'fold') continue;
        const [kL, kR] = ev.between;
        if (kL < 0 || kR >= points.length) continue;
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
      // Point 0 from current converged state — always pushed BEFORE stability() runs so `points`
      // always contains at least the start point even when `budgetMs` runs out during point 0's
      // stability check.
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
      const branchLawRefs = this.lawRefs;
      const branchChargeRefs = this.chargeRefs;
      const branchLocalKept = this.localKept;
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
        landTx: new Float64Array(N + 1),
        landTV: new Float64Array(K),
      };
      ctx = {
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
        startMissesKirchhoff: this.islands.length > 0 && this._missesKirchhoff(),
        startKirchhoffSum: this.islands.length > 0 ? kirchhoffSum(this) : 0,
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
      let lastRankDeficient = this.sys.lostColumns > 0 || this.sys.perturbedPivots > 0;
      let inNaturalMode = this.sys.lostColumns >= 10;

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
      let dsArc = Math.abs(ts) > 1e-9 ? Math.min(maxSPerStep, (0.5 * WS_BASE) / ctx.wsEff) / Math.abs(ts) : 2;
      if (isCurrentDrive) {
        const flIdx = fl.indexOf(termIdx);
        if (flIdx >= 0 && Math.abs(tV[flIdx]) > 1e-9) dsArc = Math.min(dsArc, 0.1 / Math.abs(tV[flIdx]));
      }
      const initStateRate = maxStateStep(this, tx, tV);
      if (initStateRate > 1e-12) dsArc = Math.min(dsArc, 5 / initStateRate);
      dsArc = Math.max(MIN_ARC, Math.min(MAX_ARC, dsArc));

      let consecFails = 0;
      let consecModeFails = 0;
      let recentFail = 0;
      let naturalStreak = 0;
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

        // (a) Fold: dp/ds (ts) changes sign — only emit when confirmed as a genuine turn.
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
              if (isConfirmedFold(prevState, curState, fold, absSpan)) {
                events.push({
                  kind: 'fold',
                  drive: fold.drive,
                  terminals: fold.terminals,
                  between: [kPrev, kNew],
                });
              } else if (_stats) {
                _stats.falseFoldRejected++;
              }
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

        // (b) Hopf: complex conjugate pair crosses Re λ = 0 (matched by eigenvector overlap).
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

        // Restore solver to curState so continuation proceeds cleanly.
        this._restore(curSnap);
        this.sourceOverride.set(termIdx, driveNew);
      };

      // Land cleanly on an exact span edge `sTarget` (0 or 1) from bracket `[A, B]`, enforcing
      // the SAME continuity guards (alignment, angle, step cap, state motion, voltage bound) as
      // any ordinary step — never switching off at small stride and never aliasing `refTx`.
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
        if (stateMove > maxPotPerStep * 1.5) return null;

        // Write landed tangent into dedicated work.landTx/landTV buffers so refTx/prevTx is NEVER aliased.
        const tanLand = computeTangent(this, ctx, sTarget, refTx, refTV, refTs, work.landTx, work.landTV, !!corLand.reusable);
        if (!tanLand.ok) return null;
        if (tanLand.cosTheta < COS_MIN) return null;

        let align = ctx.wsEff * ctx.wsEff * refTs * (sTarget - A.s);
        const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
        for (let i = 0; i < potIdx.length; i++) {
          const j = potIdx[i], o = this.fullOf[j], w = potWt[i];
          align += invN * w * w * refTx[j] * ((this.u[o] - A.snap.u[o]) + (this.uLo[o] - A.snap.uLo[o]));
        }
        const invVT2 = 1 / (this.VT * this.VT);
        for (let a = 0; a < K; a++) align += invVT2 * refTV[a] * (this.termV[fl[a]] - A.snap.termV[fl[a]]);
        if (!(align > 0)) return null;

        return { s: sTarget, tx: Float64Array.from(work.landTx), tV: Float64Array.from(work.landTV), tan: tanLand, dActual, stateMove };
      };

      // Item 3: Natural-parameter fallback — plain steady solve at the next drive value from the
      // last accepted state, accepted only if it passes the same continuity and mode-tracking
      // guards. Returns 'accepted' (loop continues) or 'done' (end is set, loop breaks).
      const tryNaturalStep = (triggerReason) => {
        if (_stats) _stats.naturalFallbackTried++;
        const prevStateObj = states[points.length - 1];
        const dirS = sCur <= 0 ? 1 : (ts >= 0 ? 1 : -1);
        const maxNatS = Math.min(maxSPerStep, 0.15);
        let dsTry = maxNatS;
        let lastFail = 'steady solve did not converge';
        let hitVBoundInfo = null;

        for (let attempt = 0; attempt < 8; attempt++) {
          budget?.check();
          let sTry = sCur + dirS * dsTry;
          if (dirS > 0 && sTry >= 1) sTry = 1;
          if (dirS < 0 && sTry <= 0) sTry = 0;
          const pTry = p0 + sTry * dpSpan;
          if (step !== undefined && Math.abs((sTry - sCur) * dpSpan) > step * (1 + 1e-6)) {
            dsTry = (step / absSpan) * 0.99;
            continue;
          }

          this._restore(prevStateObj.snap);
          term.drive.src = { value: pTry };
          this.sourceOverride.delete(termIdx);
          this.steady = false;
          this.allTerminals = false;
          this.solvedV = [this.termV[0], this.termV[1]];
          this.stepped = false;
          this.carried = false;
          this.coldStart = false;
          let res;
          try {
            this._refreshSources();
            res = this.solveSteady(budget ? { budgetMs: budget.remaining() } : {});
          } finally {
            term.drive.src = origSrc;
            this.steady = true;
            this.steadyAt = this.time;
            this.allTerminals = true;
            this.lawRefs = branchLawRefs;
            this.chargeRefs = branchChargeRefs;
            this.localKept = branchLocalKept;
            this.sourceOverride.set(termIdx, pTry);
            this._refreshSources();
          }

          if (res.budgetSpent) throw BUDGET_SPENT;
          if (!res.converged) {
            lastFail = res.stopped ?? 'steady solve did not converge';
            dsTry *= 0.5;
            if (dsTry < 1e-5 * Math.min(1, maxSPerStep)) break;
            continue;
          }

          // Check voltage bound.
          let vb = null;
          for (let a = 0; a < K; a++) {
            const dv = Math.abs(this.termV[fl[a]] - vStart[a]);
            if (dv > voltageBound) {
              vb = { name: this.terms[fl[a]].name, V: this.termV[fl[a]], dv };
              break;
            }
          }
          if (vb) {
            hitVBoundInfo = vb;
            dsTry *= 0.5;
            if (dsTry < 1e-4 * Math.min(1, maxSPerStep)) break;
            continue;
          }

          // Continuity guard: same state-move cap as ordinary steps.
          const actualStateMove = maxStateDiff(this, prevStateObj.snap, this.u, this.uLo, this.termV);
          if (actualStateMove > maxPotPerStep * 1.5) {
            lastFail = `state jumped by ${actualStateMove.toPrecision(3)} thermal units`;
            dsTry *= 0.5;
            if (dsTry < 1e-5 * Math.min(1, maxSPerStep)) break;
            continue;
          }

          const dActual = stateDistW(this, ctx, prevStateObj.snap.u, prevStateObj.snap.uLo, prevStateObj.snap.termV, sCur, this.u, this.uLo, this.termV, sTry);

          // All-mode continuity guard (when stability is enabled).
          let stCand = null;
          if (doStab) {
            budget?.check();
            stCand = evaluateStability(this, termIdx, pTry, stabOpts, isCurrentDrive, budget);
            const kPrev = points.length - 1;
            const sm = modeStepSmooth(
              kPrev >= 1 ? stFullList[kPrev - 1] : null,
              stFullList[kPrev],
              stCand,
              stepArcs[kPrev],
              dActual,
              actualStateMove,
              kBar,
              nWant,
              _stats,
            );
            if (!sm.ok) {
              lastFail = `mode guard (${sm.reason})`;
              dsTry *= 0.5;
              if (dsTry < 1e-5 * Math.min(1, maxSPerStep)) break;
              continue;
            }
          }

          // Natural-parameter step passed all guards! Avoid ping-ponging back into doomed
          // pseudo-arclength iterations when pseudo-arclength has repeatedly failed (stay in
          // natural mode once naturalStreak >= 2, re-probing arclength only every 5th step).
          if (_stats) _stats.naturalFallbackAccepted++;
          naturalStreak++;
          const acceptedSnap = this._snapshot();
          prevTx.set(tx);
          prevTV.set(tV);
          const prevTsLocal = ts;
          let resumedArc = false;
          const tryResume = naturalStreak < 2 || (naturalStreak % 5 === 0);
          if (tryResume) {
            try {
              this.sourceOverride.set(termIdx, pTry);
              this.islandsOn = false;
              this.assemble(Infinity);
              equilibrate(this);
              computeParamDs(this, termIdx, dpSpan, ctx.Bs, ctx.pinDs, ctx.circDs);
              this._factor();
              const rankBad = this.sys.lostColumns > 0 || this.sys.perturbedPivots > 0;
              if (!rankBad) {
                const tanTry = computeTangent(this, ctx, sTry, prevTx, prevTV, prevTsLocal, work.landTx, work.landTV, true);
                if (tanTry.ok && tanTry.cosTheta >= COS_MIN && (tanTry.ts * dirS > 0)) {
                  tx.set(work.landTx);
                  tV.set(work.landTV);
                  ts = tanTry.ts;
                  lastRankDeficient = false;
                  inNaturalMode = false;
                  dsArc = Math.max(MIN_ARC, Math.min(MAX_ARC, dActual));
                  resumedArc = true;
                }
              } else {
                lastRankDeficient = true;
              }
            } catch {
              lastRankDeficient = true;
            }
          }
          if (!resumedArc) {
            tx.fill(0);
            tV.fill(0);
            ts = dirS / ctx.wsEff;
            inNaturalMode = true;
          }
          this._restore(acceptedSnap);
          this.sourceOverride.set(termIdx, pTry);

          sCur = sTry;
          stepCount++;
          consecFails = 0;
          consecModeFails = 0;
          acceptPoint(sTry, tx, tV, ts, dActual, stCand, null);

          if (dirS > 0 && sTry >= 1) {
            end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
            return 'done';
          }
          if (dirS < 0 && sTry <= 0) {
            const vNote = isCurrentDrive ? `, V = ${points[points.length - 1].terminals[terminal].V.toPrecision(5)} V` : '';
            end = `returned to starting drive (${formatDrive(p0, kind, hasGeometry)}${vNote}) after fold`;
            return 'done';
          }
          return 'accepted';
        }

        // Fallback failed: restore last accepted state cleanly and say why.
        this._restore(prevStateObj.snap);
        tx.set(prevStateObj.tx);
        tV.set(prevStateObj.tV);
        ts = prevStateObj.ts;
        if (_stats) _stats.naturalFallbackRejected++;
        if (hitVBoundInfo) {
          end = `left voltage bound (±${voltageBound} V from start at terminal '${hitVBoundInfo.name}', reached ${hitVBoundInfo.V.toPrecision(4)} V)`;
        } else {
          end = `${triggerReason}; natural-parameter step failed (${lastFail})`;
        }
        return 'done';
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

        // If already in natural-parameter mode or a rank-deficient Jacobian collapsed the tangent
        // stride to micro-steps (e.g. depleted floating semiconductor island), take a continuous
        // natural-parameter step rather than crawling for maxSteps.
        const rankCollapsed = lastRankDeficient && Math.abs(dsStep * ts) < 1e-2 * Math.min(1, maxSPerStep);
        if (inNaturalMode || rankCollapsed) {
          const nat = tryNaturalStep(
            `continuation stopped: rank-deficient Jacobian collapsed arclength stride near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
          );
          if (nat === 'done') break;
          continue;
        }

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

        const cor = correctorNewton(this, ctx, sInit, ctx.uPred, ctx.uLoPred, ctx.vPred, sPred, ctx.wx, ctx.wV, ws, {
          maxIter: consecFails > 0 ? 20 : 35,
        });

        // Reject step and halve stride; when stride drops below MIN_ARC (or a hard halving cap is
        // hit), try Item 3's natural-parameter step instead of switching guards off or stopping cold.
        const rejectStep = (factor, reasonIfMin, forceFallback = false) => {
          this._restore(beforeSnap);
          consecFails++;
          recentFail = 3;
          dsArc = dsStep * factor;
          if (forceFallback || dsArc < MIN_ARC) {
            if (_stats && dsStep <= MIN_ARC * 4) _stats.guardAtMinArc++;
            const nat = tryNaturalStep(reasonIfMin);
            return nat === 'done' ? 'done' : 'natural';
          }
          return 'retry';
        };

        if (!cor.converged) {
          const act = rejectStep(
            consecFails > 0 ? 0.25 : 0.5,
            cor.modelError
              ? `continuation stopped: model refused (${cor.modelError})`
              : `continuation stopped: corrector failed at ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)} (${cor.error ?? 'step shrunk below minimum'})`,
            Boolean(cor.modelError) || lastRankDeficient || consecFails >= 4,
          );
          if (act === 'done') break;
          continue;
        }

        const sNext = cor.s;

        // Guard 1: Corrector displacement, total arclength, pointwise state motion, and hard `step`.
        // NEVER switches off at MIN_ARC.
        const dCor = stateDistW(this, ctx, ctx.uPred, ctx.uLoPred, predTermV, sPred, this.u, this.uLo, this.termV, sNext);
        const dActual = stateDistW(this, ctx, beforeSnap.u, beforeSnap.uLo, beforeSnap.termV, sCur, this.u, this.uLo, this.termV, sNext);
        const actualStateMove = maxStateDiff(this, beforeSnap, this.u, this.uLo, this.termV);
        const actualDp = Math.abs((sNext - sCur) * dpSpan);
        const stepExceeded = step !== undefined && actualDp > step * (1 + 1e-6);
        const dispTooLarge = dCor > 0.5 * dsStep || dActual > 1.5 * dsStep || actualStateMove > maxPotPerStep * 1.25;
        if (stepExceeded || dispTooLarge) {
          if (_stats && dsStep <= MIN_ARC * 4) _stats.guardAtMinArc++;
          let factor = 0.5;
          if (stepExceeded && actualDp > 0) factor = Math.min(0.5, 0.8 * (step / actualDp));
          else if (dCor > 0.5 * dsStep) factor = Math.max(0.25, Math.min(0.5, 0.4 * dsStep / dCor));
          const act = rejectStep(
            factor,
            `continuation stopped: step control could not bound displacement near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
          );
          if (act === 'done') break;
          continue;
        }

        // Guard 2: Compute new unit tangent, check turning angle (cosTheta ≥ 0.75) and forward
        // alignment. NEVER switches off at MIN_ARC.
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
        const candRankDeficient = this.sys.lostColumns > 0 || this.sys.perturbedPivots > 0;
        if (!tanNext.ok || tanNext.cosTheta < COS_MIN) {
          if (_stats && dsStep <= MIN_ARC * 4) _stats.guardAtMinArc++;
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          const act = rejectStep(
            0.5,
            `continuation stopped: tangent turned sharply or became singular near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
            candRankDeficient && consecFails >= 3,
          );
          if (act === 'done') break;
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
          const act = rejectStep(
            0.5,
            `continuation stopped: corrector jumped backward near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
          );
          if (act === 'done') break;
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

        // Check for fold on this step FIRST so the One Span Rule handles folds past `to` and
        // immediate fold+return past `p0` uniformly — and confirm every fold geometrically so
        // tangent noise on a rank-deficient Jacobian never emits a false fold.
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
          this._restore(candState.snap);

          // Confirm the fold geometrically: dp/ds through 0 with a real extremum and drive change
          // beyond tangent noise. If not confirmed, reject and halve (falling back to natural step).
          if (!isConfirmedFold(prevStateObj, candState, preFold, absSpan)) {
            if (_stats) _stats.falseFoldRejected++;
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            const act = rejectStep(
              0.5,
              `continuation stopped: unconfirmed tangent sign flip near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
              candRankDeficient || consecFails >= 3,
            );
            if (act === 'done') break;
            continue;
          }

          // If the fold itself lies past `to` (sFold >= 1) while sCur < 1, land on s = 1 on the
          // pre-fold leg with no fold event.
          if (prevTs > 0 && sCur < 1 && preFold.s >= 1) {
            const landed = landExact(prevStateObj, preFold.state, 1, prevTx, prevTV, prevTs);
            if (landed) {
              acceptPoint(1, landed.tx, landed.tV, landed.tan.ts, landed.dActual);
              end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
              break;
            }
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            const act = rejectStep(0.5, `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`);
            if (act === 'done') break;
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
            const act = rejectStep(0.5, `continuation stopped: could not land on starting drive ${formatDrive(p0, kind, hasGeometry)}`);
            if (act === 'done') break;
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
              const sm = modeStepSmooth(
                kPrev >= 1 ? stFullList[kPrev - 1] : null,
                stFullList[kPrev],
                stLand,
                stepArcs[kPrev],
                landed.dActual,
                landed.stateMove,
                kBar,
                nWant,
                _stats,
              );
              if (!sm.ok) {
                consecModeFails++;
                tx.set(prevTx);
                tV.set(prevTV);
                ts = prevTs;
                const act = rejectStep(
                  0.5,
                  `continuation stopped: mode guard (${sm.reason}) near ${formatDrive(to, kind, hasGeometry)}`,
                  consecModeFails >= 8,
                );
                if (act === 'done') break;
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
          const act = rejectStep(0.25, `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`);
          if (act === 'done') break;
          continue;
        }

        // One Span Rule — edge s = 0 (`p0`): crossed back to <= 0 (including on the same step as
        // a fold).
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
          const act = rejectStep(0.25, `returned past starting drive (${formatDrive(p0, kind, hasGeometry)}) after fold`);
          if (act === 'done') break;
          continue;
        }

        // Guard 4: All-mode tracking & smoothness check (when stability is enabled) BEFORE accepting!
        // NEVER switches off at MIN_ARC and NEVER gives up silently after a fixed halving count.
        let stCand = null;
        if (doStab) {
          budget?.check();
          const driveCand = p0 + sNext * dpSpan;
          stCand = evaluateStability(this, termIdx, driveCand, stabOpts, isCurrentDrive, budget);
          const kPrev = points.length - 1;
          const sm = modeStepSmooth(
            kPrev >= 1 ? stFullList[kPrev - 1] : null,
            stFullList[kPrev],
            stCand,
            stepArcs[kPrev],
            dActual,
            actualStateMove,
            kBar,
            nWant,
            _stats,
          );
          if (!sm.ok) {
            consecModeFails++;
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            const act = rejectStep(
              0.5,
              `continuation stopped: mode guard (${sm.reason}) near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`,
              consecModeFails >= 8,
            );
            if (act === 'done') break;
            continue;
          }
        }

        consecFails = 0;
        consecModeFails = 0;
        naturalStreak = 0;
        lastRankDeficient = candRankDeficient;

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
