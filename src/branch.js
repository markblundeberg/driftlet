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
// Continuation:
//   - Traces the 1D steady solution manifold F(x, V_fl, p) = 0 in (state, p), where p is the
//     chosen terminal's own drive value (V for a held or R-backed voltage terminal, I for a
//     current-driven terminal), keeping all conserved amounts (lawRefs), floating charges
//     (chargeRefs) and immobile sums (localKept) fixed at their starting values.
//   - Dimensionless parameter s = (p − p0) / (to − p0) (0 at start, 1 at `to`), oriented initially
//     with ds/d(arc) > 0, using grid-independent potential arclength weights matching homotopy.js
//     and bordered Woodbury linear solves (no extra matrix storage; pure JS, zero dependencies,
//     no ** or Math.pow).
//   - Every output point is converged to `tol = 1e-10` thermal units (and Kirchhoff-checked on
//     floating-island devices), with a tangent-angle guard (`cos θ ≥ 0.45`) that shrinks the
//     arclength stride around sharp folds so no step jumps onto an unrelated branch.
//
// Stability marks & events:
//   - With `stability !== false`, every accepted point is evaluated by a full `this.stability(...)`
//     call at that point's drive (guaranteeing agreement with a standalone `dev.stability()` call).
//     Marks are for the drive as held (a curve can be stable under V-drive and unstable under
//     I-drive).
//   - Fold (`dp/ds` changes sign): located by safeguarded secant/bisection along the exact arc to
//     `|tp| < 1e-8` (or chord < 1e-9), reporting `{ kind: 'fold', drive, terminals, between }`.
//   - Branch point (`sign(det A)` changes sign without a turn in p, where A is the tangent-bordered
//     Jacobian): verified by bisection along the exact arc; reported as
//     `{ kind: 'branch point', drive, terminals, between }` without switching branches.
//   - Hopf (`Re λ` of a complex conjugate pair `Im λ > 0` crosses 0): located by safeguarded
//     secant/bisection along the exact arc, reporting
//     `{ kind: 'hopf', drive, terminals, omega, between }`.
//   - Consistency audit: at a fold in the held parameter, a real eigenvalue must pass through 0.
//     If stability marks do not change across a fold, or flip between decided states with no fold,
//     Hopf, or branch point located, a `{ kind: 'stability mismatch', ... }` event is emitted
//     (never silently ignored).
//
// Stop conditions (every stop is named in `end`):
//   - Reached `to` (lands cleanly at p = to).
//   - Returned to starting drive `p0` after a fold (lands cleanly at p = p0 on the far branch).
//   - Left `voltageBound` (any floating terminal moved more than `voltageBound` V from its start).
//   - `maxSteps` reached.
//   - `budgetMs` spent (unlimited `budgetMs = Infinity` never reads the wall clock; finite
//     `budgetMs` allows at least one step attempt and stops cleanly when expired).
//   - Continuation step failed (stride underflow, singular bordered system, model wall, Kirchhoff).
//
// Invariants:
//   1. Completely non-destructive: state, drives, `sourceOverride`, `drivenKey`, laws, and LEARNED
//      counters are restored bit-for-bit in `finally` (even on throw or budget stop).
//   2. Deterministic when `budgetMs = Infinity` (never reads the wall clock).
//   3. Zero change to `solve()`, `advance()`, `impedance()`, `stability()`, or `homotopy.js`.

import { sourceAt } from './device.js';
import { SolverError } from './errors.js';
import { carrierZ, chargeOnly, feeds } from './laws.js';
import { captureLaws, keepLocal } from './lawRows.js';
import { addToState, clipLevels, equilibrate, maxPotentialStep } from './newton.js';
import { makeSolution } from './solution.js';

const WS = 10; // arclength weight on dimensionless s (Δs = 0.1 ↔ ~1 thermal unit RMS motion)
const BRANCH_TOL = 1e-10; // corrector Newton tolerance in thermal units (~2.5e-12 V)
const MIN_ARC = 1e-5;
const MAX_ARC = 1e5;
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
    markAttempted() {
      attempted = true;
    },
    elapsed() {
      return elapsed || now() - start;
    },
  };
}

// Dimensionless potential/level unknowns (matching potentialIndices in homotopy.js).
function potentialIndices(s) {
  const { n, M, rix: R } = s;
  const idx = [];
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (s.nodeConductor[g] < 0 && s.loc[b * M] >= 0) idx.push(R[b * M]);
    for (let i = 0; i < n; i++) if (s.present[g * n + i] && s.loc[b * M + 1 + i] >= 0) idx.push(R[b * M + 1 + i]);
    const k = s.surfAt[g];
    if (k >= 0) for (let q = 0; q < s.model.surfaces[k].species.length; q++) if (s.loc[b * M + 1 + n + q] >= 0) idx.push(R[b * M + 1 + n + q]);
  }
  for (const gt of s.gateList) if (s.loc[gt.o] >= 0) idx.push(R[gt.o]);
  for (const d of s.diffList) {
    const bf = s.blockOfFace[d.f];
    if (s.loc[d.o] >= 0) idx.push(R[d.o]);
    for (const q of d.drawSlot) if (q >= 0 && s.loc[bf * M + q] >= 0) idx.push(R[bf * M + q]);
  }
  return Int32Array.from(idx);
}

// Sign of det(T) from the already-factored BlockTridiagonal s.sys (returns +1, -1, or 0).
// Block reversal Q T Qᵀ (fromRight === true) is a similarity transform, so det sign is identical.
function blockDetSign(sys) {
  if (!sys.factored) return 0;
  const t = sys.fromRight ? sys._rev : sys;
  const { n, sizes, offB, offX, _lu: lu, _piv: piv } = t;
  let sign = 1;
  for (let i = 0; i < n; i++) {
    const m = sizes[i];
    if (m === 0) continue;
    const o = offB[i];
    const p = offX[i];
    for (let k = 0; k < m; k++) {
      if (piv[p + k] !== k) sign = -sign;
      const d = lu[o + k * m + k];
      if (!(d !== 0 && Number.isFinite(d))) return 0;
      if (d < 0) sign = -sign;
    }
  }
  return sign;
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
// { ds, detSignA, minSchurPivot } and writing δx into `delta`, δV into `deltaV`.
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
      s._solveLinear(e, cols[a]);
      for (const o of col) e[s.rix[o]] = 0;
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
  let schurSign = 1;
  let minSchurPivot = Infinity;
  for (let col = 0; col < S; col++) {
    let p = col;
    for (let r = col + 1; r < S; r++) if (Math.abs(A[r][col]) > Math.abs(A[p][col])) p = r;
    if (p !== col) {
      [A[col], A[p]] = [A[p], A[col]];
      schurSign = -schurSign;
    }
    const piv = A[col][col];
    if (!(piv !== 0 && Number.isFinite(piv))) throw new SolverError('the bordered arclength system is singular');
    const apiv = Math.abs(piv);
    if (apiv < minSchurPivot) minSchurPivot = apiv;
    if (piv < 0) schurSign = -schurSign;
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
  const detSignA = blockDetSign(s.sys) * schurSign;
  return { ds: x[P + K], detSignA, minSchurPivot };
}

function buildWeights(potIdx, VT, tx, tV, ts, wx, wV) {
  const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
  wx.fill(0);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i];
    wx[j] = invN * tx[j];
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tV.length; a++) wV[a] = invVT2 * tV[a];
  return WS * WS * ts;
}

function tangentNorm(potIdx, VT, tx, tV, ts) {
  let sum = WS * WS * ts * ts;
  const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
  for (let i = 0; i < potIdx.length; i++) {
    const v = tx[potIdx[i]];
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
// Also returns detSignA (sign of det of the tangent-bordered Jacobian A with row t^T W),
// cosTheta (alignment with prevT), and cleanFactor.
function computeTangent(s, ctx, sVal, prevTx, prevTV, prevTs, tx, tV, reuseFactor = false) {
  const { termIdx, p0, dpSpan, potIdx, Bs, pinDs, circDs, wx, wV, work } = ctx;
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
    ws = buildWeights(potIdx, s.VT, prevTx, prevTV, prevTs, wx, wV);
  } else {
    wx.fill(0);
    wV.fill(0);
    ws = 1;
  }
  let res = solveBordered(s, work, zero, zeroPin, zeroCirc, 1, tx, tV, Bs, pinDs, circDs, wx, wV, ws);
  let tsRaw = res.ds;
  if (!Number.isFinite(tsRaw)) return { ok: false };
  const nrm = tangentNorm(potIdx, s.VT, tx, tV, tsRaw);
  if (!(nrm > 0 && Number.isFinite(nrm))) return { ok: false };
  for (let j = 0; j < N; j++) tx[j] /= nrm;
  for (let a = 0; a < K; a++) tV[a] /= nrm;
  let ts = tsRaw / nrm;
  let detSignA = res.detSignA;
  let minSchurPivot = res.minSchurPivot;
  let cosTheta = 1;
  if (!prevTx) {
    if (ts < 0) {
      for (let j = 0; j < N; j++) tx[j] = -tx[j];
      for (let a = 0; a < K; a++) tV[a] = -tV[a];
      ts = -ts;
    }
    // Evaluate detSignA with the canonical tangent border w = t^T W so k = 0 and k >= 1 share
    // the exact same augmented matrix convention.
    const ws0 = buildWeights(potIdx, s.VT, tx, tV, ts, wx, wV);
    const r0 = solveBordered(s, work, zero, zeroPin, zeroCirc, 1, work.tmpX, work.tmpV, Bs, pinDs, circDs, wx, wV, ws0);
    detSignA = r0.detSignA;
    minSchurPivot = r0.minSchurPivot;
  } else {
    // Since prevT is unit-length in W and prevT^T W t_raw = 1, prevT^T W t = 1 / nrm.
    cosTheta = 1 / nrm;
  }
  const cleanFactor = s.sys.perturbedPivots === 0 && s.sys.lostColumns === 0 && minSchurPivot > 1e-12;
  return { ok: true, ts, detSignA, minSchurPivot, cosTheta, cleanFactor };
}

// Corrector Newton on F(x, V_fl, s) = 0 and wx·(x − xPred) + wV·(V − vPred) + ws·(s − sPred) = 0.
function correctorNewton(s, ctx, sStart, uPred, uLoPred, vPred, sPred, wx, wV, ws, { maxIter = 50, tol = BRANCH_TOL, maxStep = 10, fixedS = false } = {}) {
  const { termIdx, p0, dpSpan, potIdx, Bs, pinDs, circDs, work, budget } = ctx;
  budget?.check();
  const N = s.sys.size, fl = s.floating, K = fl.length;
  const delta = work.delta, deltaV = work.deltaV;
  const pinRhs = work.pinRhs, circRhs = work.circRhs;
  let sVal = fixedS ? sPred : sStart;
  for (let it = 1; it <= maxIter; it++) {
    budget?.check();
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    s.islandsOn = false;
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
    let g = ws * (sVal - sPred);
    if (!fixedS) {
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
      solved = solveBordered(s, work, s.res, pinRhs, circRhs, g, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws);
    } catch (err) {
      if (!(err instanceof SolverError)) throw err;
      return { converged: false, iterations: it, s: sVal, error: err.message };
    }
    const ds = fixedS ? 0 : solved.ds;
    clipLevels(s, delta, maxStep);
    const step = Math.max(maxPotentialStep(s, delta, deltaV), WS * Math.abs(ds));
    const own = Math.max(maxPotentialStep(s, delta, deltaV, false), WS * Math.abs(ds));
    if (!Number.isFinite(step)) return { converged: false, iterations: it, s: sVal, error: 'non-finite update' };
    if (it > 1 && own > 1e5) return { converged: false, iterations: it, s: sVal, error: 'diverging' };
    const alpha = own > maxStep ? maxStep / own : 1;
    addToState(s, delta, -alpha);
    for (let a = 0; a < K; a++) s.termV[fl[a]] -= alpha * deltaV[a];
    if (!fixedS) sVal -= alpha * ds;
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
        reusable: s.sys.perturbedPivots === 0 && s.sys.lostColumns === 0,
        minSchurPivot: solved.minSchurPivot,
      };
    }
  }
  return { converged: false, iterations: maxIter, s: sVal, error: 'Newton did not converge' };
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

// Run full stability() at the current state and drive `p`, guaranteeing agreement with a
// standalone dev.stability() call there.
function evaluateStability(s, termIdx, p, stabOpts, isCurrentDrive) {
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
  try {
    s._refreshSources();
    if (isCurrentDrive && s._drivenKey() !== wasDrivenKey) {
      s._lawsNow();
      relawed = true;
    }
    return s.stability(stabOpts);
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

// Chord-hyperplane solve at fraction theta in (0, 1) between two converged states A and B.
// Restores A's state if convergence fails.
function solveOnChord(s, ctx, A, B, theta) {
  const { potIdx, wx, wV, uPred, uLoPred, vPred } = ctx;
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
  const nrm = tangentNorm(potIdx, s.VT, chX, chV, chS);
  if (!(nrm > 0 && Number.isFinite(nrm))) return { ok: false, chordLen: 0 };
  for (let j = 0; j < N; j++) chX[j] /= nrm;
  for (let a = 0; a < K; a++) chV[a] /= nrm;
  const chSu = chS / nrm;

  // Interpolated starting guess & hyperplane point z_pred = (1 - theta) A + theta B.
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
  const ws = buildWeights(potIdx, s.VT, chX, chV, chSu, wx, wV);
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
      detSignA: tan.detSignA,
      cleanFactor: tan.cleanFactor,
      minSchurPivot: tan.minSchurPivot,
    },
  };
}

// Locate a fold (ts = 0) between A and B (tsA * tsB < 0) by safeguarded secant/bisection.
function locateFold(s, ctx, A0, B0) {
  let A = A0, B = B0;
  let best = Math.abs(A.ts) <= Math.abs(B.ts) ? A : B;
  for (let it = 0; it < 24; it++) {
    ctx.budget?.check();
    const denom = Math.abs(A.ts) + Math.abs(B.ts);
    let theta = it % 3 === 2 || !(denom > 0) ? 0.5 : Math.abs(A.ts) / denom;
    if (theta < 0.15) theta = 0.15;
    if (theta > 0.85) theta = 0.85;
    let m = solveOnChord(s, ctx, A, B, theta);
    if (!m.ok && theta !== 0.5) m = solveOnChord(s, ctx, A, B, 0.5);
    if (!m.ok) break;
    const M = m.state;
    if (Math.abs(M.ts) < Math.abs(best.ts)) best = M;
    if (Math.abs(M.ts) < 1e-8 || m.chordLen < 1e-9) {
      best = M;
      break;
    }
    if ((M.ts > 0) === (A.ts > 0)) A = M;
    else B = M;
  }
  s._restore(best.snap);
  const drive = ctx.p0 + best.s * ctx.dpSpan;
  const { terminals } = readPointData(s, ctx.termIdx, drive, false);
  return { drive, terminals, ts: best.ts };
}

// Verify and locate a branch point (detSignA changes sign without a turn in p) between A and B.
function locateBranchPoint(s, ctx, A0, B0) {
  if (!A0.cleanFactor || !B0.cleanFactor) return null;
  let A = A0, B = B0;
  let lastGood = null;
  let cleanSteps = 0;
  for (let it = 0; it < 18; it++) {
    ctx.budget?.check();
    const m = solveOnChord(s, ctx, A, B, 0.5);
    if (!m.ok) break;
    const M = m.state;
    if (M.detSignA === 0) {
      lastGood = M;
      break;
    }
    if (M.cleanFactor) cleanSteps++;
    lastGood = M;
    if (M.detSignA === A.detSignA) A = M;
    else if (M.detSignA === B.detSignA) B = M;
    else break;
    if (m.chordLen < 1e-7) break;
  }
  if (!lastGood || cleanSteps < 2 || A.detSignA === B.detSignA) return null;
  s._restore(lastGood.snap);
  const drive = ctx.p0 + lastGood.s * ctx.dpSpan;
  const { terminals } = readPointData(s, ctx.termIdx, drive, false);
  return { drive, terminals };
}

// Find the complex mode (im > 0) in `st.modes` closest in frequency to `omegaEst`.
function findComplexPair(st, omegaEst) {
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
  if (omegaEst > 0 && bestDist > 0.5 * Math.max(omegaEst, best.im)) return null;
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
  let st = evaluateStability(s, ctx.termIdx, drive, fastOpts, ctx.isCurrentDrive);
  let pair = findComplexPair(st, omegaEst);
  if (!pair) {
    st = evaluateStability(s, ctx.termIdx, drive, baseOpts, ctx.isCurrentDrive);
    pair = findComplexPair(st, omegaEst);
  }
  return pair;
}

// Locate a Hopf bifurcation (Re λ = 0 for complex conjugate pair with Im λ > 0) between A and B.
function locateHopf(s, ctx, A0, B0, pairA0, pairB0, omegaInit) {
  let A = A0, B = B0;
  let omega = omegaInit;
  let mA = pairA0 ?? evalHopfMode(s, ctx, A, omega);
  let mB = pairB0 ?? evalHopfMode(s, ctx, B, omega);
  if (!mA || !mB) return null;
  if ((mA.re > 0) === (mB.re > 0) && mA.re !== 0 && mB.re !== 0) return null;
  let best = Math.abs(mA.re) <= Math.abs(mB.re) ? { state: A, mode: mA } : { state: B, mode: mB };
  const kBar = ctx.stabOpts.kBar ?? 10;
  for (let it = 0; it < 18; it++) {
    ctx.budget?.check();
    const denom = Math.abs(mA.re) + Math.abs(mB.re);
    let theta = it % 3 === 2 || !(denom > 0) ? 0.5 : Math.abs(mA.re) / denom;
    if (theta < 0.15) theta = 0.15;
    if (theta > 0.85) theta = 0.85;
    let m = solveOnChord(s, ctx, A, B, theta);
    if (!m.ok && theta !== 0.5) m = solveOnChord(s, ctx, A, B, 0.5);
    if (!m.ok) break;
    const M = m.state;
    const mM = evalHopfMode(s, ctx, M, omega);
    if (!mM) break;
    omega = mM.im;
    if (Math.abs(mM.re) < Math.abs(best.mode.re)) best = { state: M, mode: mM };
    const pA = ctx.p0 + A.s * ctx.dpSpan;
    const pB = ctx.p0 + B.s * ctx.dpSpan;
    const pM = ctx.p0 + M.s * ctx.dpSpan;
    if (Math.abs(mM.re) <= Math.max(kBar * mM.error, 1e-4 * omega) || Math.abs(pB - pA) <= 1e-6 * Math.max(1, Math.abs(pM))) {
      best = { state: M, mode: mM };
      break;
    }
    if ((mM.re > 0) === (mA.re > 0)) { A = M; mA = mM; }
    else { B = M; mB = mM; }
  }
  s._restore(best.state.snap);
  const drive = ctx.p0 + best.state.s * ctx.dpSpan;
  const { terminals } = readPointData(s, ctx.termIdx, drive, false);
  return { drive, terminals, omega: best.mode.im };
}

// Detect whether any complex conjugate pair crossed Re λ = 0 between two stability results.
function detectHopfBracket(stA, stB, kBar) {
  const sig = (m) => (m.re > kBar * m.error ? 1 : m.re < -kBar * m.error ? -1 : 0);
  // First: matched complex pair present on both sides with opposite clear signs.
  for (const ma of stA.modes) {
    if (!(ma.im > 0)) continue;
    const sa = sig(ma);
    if (sa === 0) continue;
    const mb = findComplexPair(stB, ma.im);
    if (!mb) continue;
    const sb = sig(mb);
    if (sb !== 0 && sa !== sb) return { pairA: ma, pairB: mb, omega: 0.5 * (ma.im + mb.im) };
  }
  // Second: verdict flipped between 'no growing mode found' and 'unstable' where the unstable
  // side's growing mode is a complex pair (the stable side may have had slower real modes ahead
  // of the complex pair).
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
   * continuation, marking each point with its linear stability and locating folds, Hopf
   * bifurcations, and branch points.
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
    const stFullList = []; // parallel to points: full StabilityResult (when doStab)
    const states = []; // parallel to points: { snap, s, tx, tV, ts, detSignA, cleanFactor }
    let end = '';

    try {
      // Build point 0 from the current converged state.
      const pt0Data = readPointData(this, termIdx, p0, solutions);
      let st0 = null;
      if (doStab) {
        st0 = evaluateStability(this, termIdx, p0, stabOpts, isCurrentDrive);
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

      // Prepare steady equations & laws for the branch. If starting from I = 0 on a current-driven
      // terminal, evaluate laws at a non-zero current (`to`) so noCurrent reflects the active branch.
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
      const N = this.sys.size, fl = this.floating, K = fl.length;
      const potIdx = potentialIndices(this);
      const maxPins = Math.max(1, this.constraints.length + this.affinities.length + this.islands.length + 4);
      const work = {
        e: new Float64Array(N + 1),
        cols: [],
        zero: new Float64Array(N + 1),
        zeroPin: new Float64Array(maxPins),
        zeroCirc: new Float64Array(K),
        delta: new Float64Array(N + 1),
        deltaV: new Float64Array(K),
        pinRhs: new Float64Array(maxPins),
        circRhs: new Float64Array(K),
        tmpX: new Float64Array(N + 1),
        tmpV: new Float64Array(K),
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

      // 1. Initial unit tangent at s = 0 (oriented toward `to`, i.e. ts >= 0).
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
      let ts = tan0.ts;
      let detSignA = tan0.detSignA;
      let cleanFactor = tan0.cleanFactor;
      let sCur = 0;

      states.push({
        snap: this._snapshot(),
        s: 0,
        tx: Float64Array.from(tx),
        tV: Float64Array.from(tV),
        ts,
        detSignA,
        cleanFactor,
      });

      // Step-size parameters.
      const absSpan = Math.abs(dpSpan);
      const maxSPerStep = step !== undefined ? step / absSpan : 0.12;
      const maxPotPerStep = step !== undefined ? 10 : 6; // thermal units (~0.155 V default)
      const maxDVFloat = 0.25; // V per step on floating terminals
      let dsArc = Math.abs(ts) > 1e-6 ? Math.min(maxSPerStep, 0.08) / Math.abs(ts) : 2;
      if (isCurrentDrive) {
        const flIdx = fl.indexOf(termIdx);
        if (flIdx >= 0 && Math.abs(tV[flIdx]) > 1e-9) dsArc = Math.min(dsArc, 0.1 / Math.abs(tV[flIdx]));
      }
      dsArc = Math.max(MIN_ARC, Math.min(MAX_ARC, dsArc));

      let foldSeen = false;
      let consecFails = 0;
      let recentFail = 0;
      let stepCount = 0;

      // Record an accepted point (and detect/locate fold, branch point, and Hopf on [k, k+1]).
      const acceptPoint = (sNew, newTx, newTV, newTs, newDetSignA, newCleanFactor) => {
        const kPrev = points.length - 1;
        const prevState = states[kPrev];
        const curSnap = this._snapshot();
        const curState = {
          snap: curSnap,
          s: sNew,
          tx: Float64Array.from(newTx),
          tV: Float64Array.from(newTV),
          ts: newTs,
          detSignA: newDetSignA,
          cleanFactor: newCleanFactor,
        };
        const driveNew = p0 + sNew * dpSpan;
        const ptData = readPointData(this, termIdx, driveNew, solutions);
        let stNew = null;
        if (doStab) {
          budget?.check();
          stNew = evaluateStability(this, termIdx, driveNew, stabOpts, isCurrentDrive);
        }
        points.push({
          drive: driveNew,
          terminals: ptData.terminals,
          ...(doStab ? { stability: summarizeStability(stNew) } : {}),
          ...(solutions ? { solution: ptData.solution } : {}),
        });
        stFullList.push(stNew);
        states.push(curState);

        // Event detection on segment [kPrev, kPrev + 1]:
        const kNew = points.length - 1;
        let segmentHadFold = false;
        let segmentHadBP = false;

        // (a) Fold: dp/ds (ts) changes sign.
        if ((prevState.ts > 0 && newTs < 0) || (prevState.ts < 0 && newTs > 0) || (newTs === 0 && prevState.ts !== 0)) {
          foldSeen = true;
          segmentHadFold = true;
          const fold = locateFold(this, ctx, prevState, curState);
          events.push({
            kind: 'fold',
            drive: fold.drive,
            terminals: fold.terminals,
            between: [kPrev, kNew],
          });
        } else if (
          prevState.detSignA !== 0 &&
          newDetSignA !== 0 &&
          prevState.detSignA !== newDetSignA &&
          prevState.cleanFactor &&
          newCleanFactor
        ) {
          // (b) Determinant sign change without a turn in p -> verify branch point by bisection.
          const bp = locateBranchPoint(this, ctx, prevState, curState);
          if (bp) {
            segmentHadBP = true;
            events.push({
              kind: 'branch point',
              drive: bp.drive,
              terminals: bp.terminals,
              between: [kPrev, kNew],
            });
          }
        }

        // (c) Hopf: complex conjugate pair crosses Re λ = 0.
        if (doStab && stFullList[kPrev] && stNew) {
          const hopfBracket = detectHopfBracket(stFullList[kPrev], stNew, kBar);
          if (hopfBracket) {
            const hopf = locateHopf(this, ctx, prevState, curState, hopfBracket.pairA, hopfBracket.pairB, hopfBracket.omega);
            if (hopf) {
              events.push({
                kind: 'hopf',
                drive: hopf.drive,
                terminals: hopf.terminals,
                omega: hopf.omega,
                between: [kPrev, kNew],
              });
            }
          }
        }

        // Restore solver to curState so continuation proceeds from exact curSnap.
        this._restore(curSnap);
        this.sourceOverride.set(termIdx, driveNew);
        return { segmentHadFold, segmentHadBP };
      };

      // Helper to land cleanly on an exact target sTarget (either s = 1 for `to` or s = 0 for return).
      const landExact = (beforeSnap, sBefore, sAfter, sTarget) => {
        const afterSnap = this._snapshot();
        const frac = sAfter !== sBefore ? (sTarget - sBefore) / (sAfter - sBefore) : 1;
        const fClamped = Math.max(0, Math.min(1, frac));
        this._restore(beforeSnap);
        for (let o = 0; o < this.u.length; o++) {
          this.u[o] = beforeSnap.u[o] + fClamped * (afterSnap.u[o] - beforeSnap.u[o]);
          this.uLo[o] = beforeSnap.uLo[o] + fClamped * (afterSnap.uLo[o] - beforeSnap.uLo[o]);
        }
        for (let a = 0; a < K; a++) {
          this.termV[fl[a]] = beforeSnap.termV[fl[a]] + fClamped * (afterSnap.termV[fl[a]] - beforeSnap.termV[fl[a]]);
        }
        ctx.wx.fill(0);
        ctx.wV.fill(0);
        const corLand = correctorNewton(this, ctx, sTarget, ctx.uPred, ctx.uLoPred, ctx.vPred, sTarget, ctx.wx, ctx.wV, 1, { fixedS: true });
        if (!corLand.converged) {
          this._restore(afterSnap);
          return null;
        }
        if (this.islands.length > 0 && this._missesKirchhoff()) {
          this._restore(afterSnap);
          return null;
        }
        const tanLand = computeTangent(this, ctx, sTarget, tx, tV, ts, prevTx, prevTV, !!corLand.reusable);
        if (!tanLand.ok) {
          this._restore(afterSnap);
          return null;
        }
        return { s: sTarget, tx: prevTx, tV: prevTV, tan: tanLand };
      };

      // 2. Main pseudo-arclength continuation loop.
      while (stepCount < maxSteps) {
        budget?.check();

        // Cap stride by drive step, pointwise potential step, and floating terminal voltage change.
        let dsStep = dsArc;
        if (Math.abs(ts) > 1e-9) dsStep = Math.min(dsStep, maxSPerStep / Math.abs(ts));
        const potRate = maxPotentialStep(this, tx, tV);
        if (potRate > 1e-12) dsStep = Math.min(dsStep, maxPotPerStep / potRate);
        for (let a = 0; a < K; a++) {
          const vt = Math.abs(tV[a]);
          if (vt > 1e-9) dsStep = Math.min(dsStep, maxDVFloat / vt);
        }
        if (dsStep < MIN_ARC) dsStep = MIN_ARC;

        // If the predictor is about to step cleanly past s = 1 without having folded, cap to land near 1.
        if (!foldSeen && sCur < 1 && ts > 1e-9 && sCur + dsStep * ts > 1.05) {
          dsStep = Math.max(MIN_ARC, (1.02 - sCur) / ts);
        }

        const beforeSnap = this._snapshot();
        const sPred = sCur + dsStep * ts;
        const maxPt = potRate * dsStep;
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
        const ws = buildWeights(potIdx, this.VT, tx, tV, ts, ctx.wx, ctx.wV);

        const cor = correctorNewton(this, ctx, sInit, ctx.uPred, ctx.uLoPred, ctx.vPred, sPred, ctx.wx, ctx.wV, ws);
        budget?.markAttempted();

        if (!cor.converged) {
          this._restore(beforeSnap);
          consecFails++;
          recentFail = 3;
          dsArc = dsStep * (consecFails > 1 ? 0.25 : 0.5);
          if (dsArc < MIN_ARC) {
            end = cor.modelError
              ? `continuation stopped: model refused (${cor.modelError})`
              : `continuation stopped: corrector failed at ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)} (${cor.error ?? 'step shrunk below minimum'})`;
            break;
          }
          continue;
        }

        // Verify Kirchhoff on devices with floating islands; reject step if a coarse step lost balance.
        if (this.islands.length > 0 && this._missesKirchhoff()) {
          this._restore(beforeSnap);
          consecFails++;
          recentFail = 3;
          dsArc = dsStep * 0.5;
          if (dsArc < MIN_ARC) {
            end = `continuation stopped: floating-island current balance lost at ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`;
            break;
          }
          continue;
        }

        const sNext = cor.s;

        // Compute new unit tangent and verify turning angle (no sharp jump onto another branch).
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
        if (!tanNext.ok || (tanNext.cosTheta < 0.45 && dsStep > MIN_ARC * 4)) {
          this._restore(beforeSnap);
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          consecFails++;
          recentFail = 3;
          dsArc = dsStep * 0.5;
          if (dsArc < MIN_ARC) {
            end = `continuation stopped: tangent turned sharply or became singular near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`;
            break;
          }
          continue;
        }

        // Guard against backward-jumping corrector (step moved opposite to tangent).
        let align = WS * WS * prevTs * (sNext - sCur);
        const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
        for (let i = 0; i < potIdx.length; i++) {
          const j = potIdx[i], o = this.fullOf[j];
          align += invN * prevTx[j] * ((this.u[o] - beforeSnap.u[o]) + (this.uLo[o] - beforeSnap.uLo[o]));
        }
        const invVT2 = 1 / (this.VT * this.VT);
        for (let a = 0; a < K; a++) align += invVT2 * prevTV[a] * (this.termV[fl[a]] - beforeSnap.termV[fl[a]]);
        if (!(align > 0)) {
          this._restore(beforeSnap);
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          consecFails++;
          recentFail = 3;
          dsArc = dsStep * 0.5;
          if (dsArc < MIN_ARC) {
            end = `continuation stopped: corrector jumped backward near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`;
            break;
          }
          continue;
        }

        consecFails = 0;

        // Check voltage bound on floating terminals BEFORE accepting a point outside the bound.
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
          end = `left voltage bound (±${voltageBound} V from start at terminal '${hitVBound.name}', reached ${hitVBound.V.toPrecision(4)} V)`;
          break;
        }

        // Check if this step crossed `to` (s = 1) -> land cleanly at s = 1 and stop!
        if ((sCur < 1 && sNext >= 1) || (sCur > 1 && sNext <= 1)) {
          // Restore prev tangent for landExact's orientation, then land at s = 1.
          const savedTx = Float64Array.from(tx), savedTV = Float64Array.from(tV);
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          const landed = landExact(beforeSnap, sCur, sNext, 1);
          if (landed) {
            acceptPoint(1, landed.tx, landed.tV, landed.tan.ts, landed.tan.detSignA, landed.tan.cleanFactor);
            end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
            break;
          }
          // If exact landing failed, shrink step and try approaching closer.
          tx.set(savedTx);
          tV.set(savedTV);
          this._restore(beforeSnap);
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          dsArc = dsStep * 0.25;
          if (dsArc < MIN_ARC) {
            end = `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`;
            break;
          }
          continue;
        }

        // Check if a post-fold returning branch crossed back through the starting drive (s = 0) ->
        // land cleanly at s = 0 on the far branch and stop!
        if (foldSeen && sCur > 0 && sNext <= 0) {
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          const landed = landExact(beforeSnap, sCur, sNext, 0);
          if (landed) {
            acceptPoint(0, landed.tx, landed.tV, landed.tan.ts, landed.tan.detSignA, landed.tan.cleanFactor);
            const vNote = isCurrentDrive ? `, V = ${points[points.length - 1].terminals[terminal].V.toPrecision(5)} V` : '';
            end = `returned to starting drive (${formatDrive(p0, kind, hasGeometry)}${vNote}) after fold`;
            break;
          }
          this._restore(beforeSnap);
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          dsArc = dsStep * 0.25;
          if (dsArc < MIN_ARC) {
            end = `returned past starting drive (${formatDrive(p0, kind, hasGeometry)}) after fold`;
            break;
          }
          continue;
        }

        // Accept this point.
        ts = tanNext.ts;
        detSignA = tanNext.detSignA;
        cleanFactor = tanNext.cleanFactor;
        sCur = sNext;
        stepCount++;
        acceptPoint(sNext, tx, tV, ts, detSignA, cleanFactor);

        // Adapt arclength stride.
        if (recentFail > 0) {
          recentFail--;
          if (cor.iterations <= 8) dsArc = Math.min(MAX_ARC, dsStep * 1.25);
          else dsArc = dsStep;
        } else if (cor.iterations <= 5 && tanNext.cosTheta > 0.92) {
          dsArc = Math.min(MAX_ARC, dsStep * 1.5);
        } else if (cor.iterations <= 8 && tanNext.cosTheta > 0.8) {
          dsArc = Math.min(MAX_ARC, dsStep * 1.2);
        } else if (cor.iterations >= 15 || tanNext.cosTheta < 0.7) {
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

        // (a) Every fold must correspond to a stability change (or positive real mode count change).
        for (const ev of events) {
          if (ev.kind !== 'fold') continue;
          const [kL, kR] = ev.between;
          let iL = kL;
          while (iL > 0 && isNeutral(points[iL].stability.verdict) && kL - iL < 2) iL--;
          let iR = kR;
          while (iR < points.length - 1 && isNeutral(points[iR].stability.verdict) && iR - kR < 2) iR++;
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

        // (b) Every flip between decided non-neutral states must be explained by a fold, Hopf, or
        //     branch point on the spanning interval.
        let prevDec = -1;
        for (let i = 0; i < points.length; i++) {
          const v = points[i].stability.verdict;
          if (!isDecided(v)) continue;
          if (prevDec >= 0 && points[prevDec].stability.verdict !== v) {
            const covered = events.some(
              (ev) =>
                (ev.kind === 'fold' || ev.kind === 'hopf' || ev.kind === 'branch point') &&
                ev.between[0] >= prevDec &&
                ev.between[1] <= i,
            );
            if (!covered) {
              // If intermediate neutral point(s) split a Hopf crossing across [prevDec, i], locate it now.
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
                    });
                    resolved = true;
                  }
                } catch (e) {
                  if (e !== BUDGET_SPENT) throw e;
                }
              }
              if (!resolved) {
                const mid = Math.floor((prevDec + i) / 2);
                events.push({
                  kind: 'stability mismatch',
                  drive: points[mid].drive,
                  terminals: points[mid].terminals,
                  between: [prevDec, i],
                  note: `stability changed from '${points[prevDec].stability.verdict}' to '${v}' with no fold, Hopf, or branch point located`,
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
