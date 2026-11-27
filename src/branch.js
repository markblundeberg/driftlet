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
//     initially with ds/d(arc) > 0. The arclength metric combines the RMS of state unknowns with
//     floating terminal voltages in thermal units and a dynamically rebalanced drive weight
//     `wsEff ≥ 1` so the state carries at least half the tangent weight (`wsEff·|ts| ≤ 0.85` when
//     the state responds) while insensitive terminals retain dimensionless span scaling. All drive
//     tolerances are relative to the local drive/bracket scale — never absolute in drive units.
//     Bordered Woodbury solves reuse s.sys (pure JS, zero dependencies, no ** or Math.pow), using
//     a single locked interior current cut (_currentCuts) per corrector or tangent solve when a
//     floating contact's own flux derivative cancels >3 digits (e.g. heavily doped semiconductor
//     contacts), and cut-aware difference-form GMRES refinement (_assembleDifference) on lost-mode
//     devices (depleted pools or islands).
//   - Every output point is converged to `tol = 1e-10` thermal units in state AND `1e-9` relative
//     to the local drive scale (`max(|p|, |∂I/∂V|·V_T)` for I-drive, `max(|p|, V_T)` for V-drive —
//     independent of `to`). Before accepting a floating-terminal point, the driven terminal's own
//     reported equation (contact flux and, when active, cut current) and two-terminal Kirchhoff
//     balance are verified directly on the converged state. On devices with floating islands
//     (`this.islands.length > 0`), if the plain corrector leaves Kirchhoff unbalanced, a second-
//     phase islands round (`islandsOn = true`, fixed s) polishes the island levels (as newton() does).
//   - Step acceptance guards (applied uniformly; never switched off at the minimum stride):
//       1. Hard `step` cap (`|p_next − p_cur| ≤ step` when `opts.step` is given) and pointwise
//          state-step bound (`maxStateDiff ≤ 7.5` thermal units, with `10×` weight on `[0, 1]` gates).
//       2. Chord-tangent, tangent-turn & Hermite continuity (`stepHasHiddenS`): unit chord agrees
//          with both endpoint tangents (`min(c·t_A, c·t_B) ≥ 0.85`), successive tangents satisfy
//          `t_A·t_B ≥ 0.90` (≤ 26° per step), and the Hermite cubic along the arc does not reverse
//          inside `(0, 1)` (catches narrow S-curves inside one step).
//       3. Oscillatory nearer-end guard (when `stability !== false`): unique modes (im ≥ 0) are
//          matched 1-to-1 by biorthogonal eigenvector overlap. On the same side of Re λ = 0, every
//          matched oscillatory pair (`im > kBar · error`) satisfies
//          `|Re λ_B − Re λ_A| ≤ 0.75 · max(min(|Re λ_A|, |Re λ_B|), kBar · error, 1e-3 · ω)`.
//       4. Morse-index & undecided-mark refinement (when `stability !== false`): every change in the
//          count of strictly growing modes `N_+` across a step must be accounted for by a located
//          fold (`|ΔN_+| = 1` and `N_+(fold) = min(N_+(A), N_+(B))`, not counting the fold's own
//          critical mode: the real mode nearest zero at the located fold, nearer zero than the
//          nearest real mode at either end of the step) or Hopf (located cleanly with
//          no hidden fold or pair loss inside the bracket). An `undecided` mark has unknown `N_+`
//          and triggers refinement across fold steps and count-changing transitions (up to
//          `MAX_REFINE` halvings); if still unresolved, `runAudit` does not count a fold with an
//          undecided side as covering a stability flip and emits a `stability mismatch` event instead.
//   - When pseudo-arclength cannot advance at the minimum stride (corrector failure, singular
//     Jacobian, or model wall), continuation stops honestly and notes that `solve()` may continue
//     from the last accepted drive.
//
// One span rule:
//   - The trace stays inside the requested drive span `s ∈ [0, 1]` (between start drive `p0` and
//     `to`). Initial orientation is confirmed by the first accepted step's actual drive (flipping
//     if noisy `ts` at `s = 0` pointed out of the span), and steps from `s ≤ 0` cannot wander into
//     `s < 0`. Subsequent tangents are oriented as a whole by the sign of their inner product with
//     the previous tangent in the trace's metric (never overwriting `ts` alone). The trace ends
//     where it first leaves the span, landing cleanly on the crossed edge (`to` or `p0`) and naming
//     which edge in `end`. When a step turns a fold, the fold's drive is located first: if the fold
//     lies past `to` (sFold ≥ 1), the trace lands on `to` on the pre-fold leg with no fold event;
//     if a returning leg crosses back through `p0` (even on the same step as the fold), it records
//     the fold and lands on `p0` on the far branch.
//
// Stability marks & events (and what a mark does not promise):
//   - `points` always contains at least the starting point `points[0]`. If `budgetMs` runs out
//     during the start point's `stability()` call, `points[0]` is still kept with
//     `stability: { verdict: 'undecided', checked: false, ..., reason }`.
//   - With `stability !== false`, every accepted point is evaluated by a full `this.stability(...)`
//     call at that point's state and drive, agreeing with a standalone `dev.stability()` call at
//     the same state. Marks are for the drive as held (a curve can be stable under V-drive and
//     unstable under I-drive).
//   - What a mark does NOT promise: with `stability: false`, folds narrower than a step (e.g. a
//     narrow S stepped over cleanly) can be missed because there is no stability signal to trigger
//     refinement. With `stability !== false`, modes within `kBar · error` of zero are not resolved
//     and a narrow unstable window unseen by all listed modes at both endpoints can still be missed.
//   - Fold (`dp/ds` turns against the current travel direction AND accepted drives reverse): located
//     by safeguarded secant/bisection along the exact arc to `|ts| ≤ 1e-8 · max(|tsA|, |tsB|)` (or
//     relative chord `< 1e-9`). If the fold step itself retreats past the previous point by more
//     than the drive round-off floor, the fold is recorded immediately; otherwise it is held as a
//     provisional fold and confirmed only when a subsequent accepted point's drive retreats beyond
//     that floor. A provisional fold whose subsequent accepted point advances forward instead is
//     dropped. One still pending when the trace ends is emitted, with a `note` saying no later point
//     confirmed it, only if its located drive lies beyond both of its step's accepted drives by
//     more than the floor; otherwise it is dropped. On a parameter
//     plateau, continuation walks on with no invented folds until `to`, `voltageBound`, `maxSteps`,
//     or `budgetMs`.
//   - Hopf (`Re λ` of a complex conjugate pair `Im λ > kBar · error` crosses 0): matched by
//     eigenvector overlap and located by safeguarded secant/bisection along the exact arc to
//     `|pB − pA| ≤ 1e-6 · max(|pA|, |pB|, |pB0 − pA0|)` or relative chord `≤ 1e-9`, reporting
//     `{ kind: 'hopf', drive, terminals, omega, between }`.
//   - Fallbacks are always reported: if a fold or Hopf location fails to converge or is cut short
//     by `budgetMs`, the event is still emitted with `located: false` and a `note` saying why.
//   - Consistency audit: runs on both normal completion and `budgetMs` stops. If stability marks
//     flip between decided states with no covering fold (with non-undecided endpoints) or Hopf
//     located after refinement hits resolution or budget, a `{ kind: 'stability mismatch', ... }`
//     event is emitted.
//
// Stop conditions (every stop is named in `end`):
//   - Reached `to` (lands cleanly at p = to).
//   - Returned to starting drive `p0` after a fold (lands cleanly at p = p0 on the far branch).
//   - Left `voltageBound` (any floating terminal moved more than `voltageBound` V from its start).
//   - `maxSteps` reached.
//   - `budgetMs` spent (unlimited `budgetMs = Infinity` never reads the wall clock; finite
//     `budgetMs` is checked inside every Newton iteration, shift factorisation, mode refinement,
//     and event bisection).
//   - Continuation step failed (stride underflow, singular bordered system, model wall, Kirchhoff).
//
// Invariants:
//   1. Completely non-destructive: state, drives, `sourceOverride`, `drivenKey`, laws, and LEARNED
//      counters are restored bit-for-bit in `finally` (even on throw or budget stop).
//   2. Deterministic when `budgetMs = Infinity` (never reads the wall clock).
//   3. Zero change to `solve()`, `advance()`, `impedance()`, `stability()` (when called directly),
//      or `homotopy.js`.

import { FARADAY } from './constants.js';
import { sourceAt } from './device.js';
import { SolverError } from './errors.js';
import { carrierZ, chargeOnly, feeds } from './laws.js';
import { captureLaws, keepLocal } from './lawRows.js';
import { addToState, clipLevels, equilibrate, maxPotentialStep } from './newton.js';
import { makeSolution } from './solution.js';

const GATE_WEIGHT = 10; // weight on [0, 1] gate fractions in step-size caps (Δgate = 0.1 ↔ 1 thermal unit)
const BRANCH_TOL = 1e-10; // corrector Newton tolerance in thermal units (~2.5e-12 V)
const MIN_ARC = 1e-6;
const MAX_ARC = 1e5;
const ALPHA_RE = 0.75; // nearer-end mode rule: |ΔRe λ| ≤ ALPHA_RE · max(min(|a|, |b|), kBar · err, 1e-3 · ω)
const MAX_REFINE = 10; // max halvings to explain a Morse-index flip, undecided end, or fold-adjacent crossing
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

// Dimensionless state unknown indices (uniform RMS across all state unknowns + floating V/VT).
function stateMetricIndices(s) {
  const { n, M, rix: R } = s;
  const idx = [];
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (s.nodeConductor[g] < 0 && s.loc[b * M] >= 0) idx.push(R[b * M]);
    for (let i = 0; i < n; i++) {
      if (s.present[g * n + i] && s.loc[b * M + 1 + i] >= 0) idx.push(R[b * M + 1 + i]);
    }
    const k = s.surfAt[g];
    if (k >= 0) {
      for (let q = 0; q < s.model.surfaces[k].species.length; q++) {
        if (s.loc[b * M + 1 + n + q] >= 0) idx.push(R[b * M + 1 + n + q]);
      }
    }
  }
  for (const gt of s.gateList) if (s.loc[gt.o] >= 0) idx.push(R[gt.o]);
  for (const d of s.diffList) {
    const bf = s.blockOfFace[d.f];
    if (s.loc[d.o] >= 0) idx.push(R[d.o]);
    for (const q of d.drawSlot) if (q >= 0 && s.loc[bf * M + q] >= 0) idx.push(R[bf * M + q]);
  }
  return Int32Array.from(idx);
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

function stateOnlyNorm(potIdx, invN, VT, tx, tV) {
  let sum = 0;
  for (let i = 0; i < potIdx.length; i++) {
    const v = tx[potIdx[i]];
    sum += invN * v * v;
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tV.length; a++) sum += invVT2 * tV[a] * tV[a];
  return Math.sqrt(sum);
}

function tangentInner(potIdx, invN, wsEff, VT, txA, tVA, tsA, txB, tVB, tsB) {
  let sum = wsEff * wsEff * tsA * tsB;
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i];
    sum += invN * txA[j] * txB[j];
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tVA.length; a++) sum += invVT2 * tVA[a] * tVB[a];
  return sum;
}

function tangentNorm(potIdx, invN, wsEff, VT, tx, tV, ts) {
  return Math.sqrt(tangentInner(potIdx, invN, wsEff, VT, tx, tV, ts, tx, tV, ts));
}

function buildWeights(potIdx, invN, wsEff, VT, tx, tV, ts, wx, wV) {
  wx.fill(0);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i];
    wx[j] = invN * tx[j];
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < tV.length; a++) wV[a] = invVT2 * tV[a];
  return wsEff * wsEff * ts;
}

function stateDistW(s, ctx, uA, uLoA, termVA, sA, uB, uLoB, termVB, sB) {
  const { potIdx, invN, wsEff } = ctx;
  const fl = s.floating, VT = s.VT;
  let sum = wsEff * wsEff * (sB - sA) * (sB - sA);
  for (let i = 0; i < potIdx.length; i++) {
    const o = s.fullOf[potIdx[i]];
    const d = (uB[o] - uA[o]) + (uLoB[o] - uLoA[o]);
    sum += invN * d * d;
  }
  const invVT2 = 1 / (VT * VT);
  for (let a = 0; a < fl.length; a++) {
    const k = fl[a];
    const d = termVB[k] - termVA[k];
    sum += invVT2 * d * d;
  }
  return Math.sqrt(sum);
}

// Check whether the Hermite cubic interpolating (s_A, ts_A) and (s_B, ts_B) over arclength d
// reverses direction inside (0, 1) (an S-curve / double fold hidden inside one step) or whether
// the unit chord departs sharply from either endpoint tangent.
function stepHasHiddenS(s, ctx, snapA, sA, tAx, tAV, tAs, uB, uLoB, termVB, sB, tBx, tBV, tBs) {
  const { potIdx, invN, wsEff } = ctx;
  const fl = s.floating, invVT2 = 1 / (s.VT * s.VT);
  const dS = sB - sA;
  let d2 = wsEff * wsEff * dS * dS;
  let dotA = wsEff * wsEff * tAs * dS;
  let dotB = wsEff * wsEff * tBs * dS;
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i], o = s.fullOf[j];
    const d = (uB[o] - snapA.u[o]) + (uLoB[o] - snapA.uLo[o]);
    d2 += invN * d * d;
    dotA += invN * tAx[j] * d;
    dotB += invN * tBx[j] * d;
  }
  for (let a = 0; a < fl.length; a++) {
    const k = fl[a];
    const d = termVB[k] - snapA.termV[k];
    d2 += invVT2 * d * d;
    dotA += invVT2 * tAV[a] * d;
    dotB += invVT2 * tBV[a] * d;
  }
  const dLen = Math.sqrt(d2);
  if (!(dLen > 1e-12)) return false;
  if (dotA < 0.85 * dLen || dotB < 0.85 * dLen) return true;
  if (tangentInner(potIdx, invN, wsEff, s.VT, tAx, tAV, tAs, tBx, tBV, tBs) < 0.90) return true;
  if (tAs * tBs > 0 && wsEff * Math.max(Math.abs(tAs), Math.abs(tBs), Math.abs(dS) / dLen) > 1e-6) {
    const sgn = tAs > 0 ? 1 : -1;
    const a = sgn * tAs, b = sgn * tBs, c = sgn * (dS / dLen);
    if (c <= 0) return true;
    const den = a + b - 2 * c;
    if (den > 1e-12) {
      const num = 2 * a + b - 3 * c;
      const tau = num / (3 * den);
      if (tau > 0 && tau < 1 && a - (num * num) / (3 * den) <= 0) return true;
    }
  }
  return false;
}

// Drive round-off floor estimated from the terminal current's gross species fluxes (matching
// _missesKirchhoff in src/steady.js) and Newton tolerance. Turns below this floor are plateau noise.
function driveRoundoffFloor(s, termIdx, p) {
  const d = s.terms[termIdx].drive;
  const { z } = s;
  let mx = 0, gross = 0;
  for (const I of s.termI) mx = Math.max(mx, Math.abs(I));
  for (const flux of [s.contactFlux.left, s.contactFlux.right, ...s.portFlux]) {
    for (let i = 0; i < flux.length; i++) gross += Math.abs(z[i] * FARADAY * flux[i]);
  }
  const { area } = s.model.grid;
  const iFloor = Math.max(
    1e-11 * Math.max(Math.abs(p), mx),
    1e-11 * gross,
    1e-18 * Math.max(area[0], area[area.length - 1]),
  );
  if (d.kind === 'I') return iFloor;
  const vFloor = 1e-11 * Math.max(Math.abs(p), s.VT);
  return d.R > 0 ? Math.max(vFloor, d.R * iFloor) : vFloor;
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

// Evaluate a linearised interior current cut (from s._currentCuts()) on compact vector v.
function dotCut(s, cut, sgn, v) {
  let t = 0;
  const { loc, rix } = s;
  const idx = cut.idx, g = cut.g;
  for (let q = 0; q < idx.length; q++) {
    const o = idx[q];
    if (loc[o] >= 0) t += g[q] * v[rix[o]];
  }
  return sgn * t;
}

// Raw contact/port flux current entering the device at steady state (never echoing passedCurrent).
function contactCurrent(s, k) {
  const flux = k === 0 ? s.contactFlux.left : k === 1 ? s.contactFlux.right : s.portFlux[k - 2];
  let I = 0;
  for (let i = 0; i < s.n; i++) I += FARADAY * s.z[i] * flux[i];
  return k === 1 ? -I : I;
}

// Look up a cached linearised interior current cut by segment index (computing _currentCuts once per assembly).
function getCut(s, ctx, seg) {
  let allCuts = ctx.cutsCache;
  if (!allCuts) {
    allCuts = s._currentCuts();
    ctx.cutsCache = allCuts;
  }
  for (let ci = 0; ci < allCuts.length; ci++) if (allCuts[ci].s === seg) return allCuts[ci];
  return null;
}

// Exact nonlinear steady conduction current entering from contact fk (sgn = +1 left, -1 right)
// across interior segment `seg`, matching dotCut's linearisation identically on every Newton step.
function cutCurrent(s, seg, sgn) {
  const { model, n, z } = s;
  const { grid, materials, regions } = model;
  let r = 0;
  while (r + 1 < regions.length && seg >= grid.regionEnd[r]) r++;
  const reg = regions[r], mat = materials[reg.material];
  const bL = seg + r, bR = bL + 1, h = grid.segLength[seg] / grid.segArea[seg], vel = reg.velocity;
  let I = 0;
  for (let i = 0; i < n; i++) {
    if (z[i] === 0 || !mat.present[i] || mat.D[i] === 0) continue;
    const N = mat.ideal ? s._sgSegment(seg, bL, i, h, vel)[0] : s._fluxConcentrated(seg, bL, bR, mat, h, vel, i).N;
    I += z[i] * FARADAY * N;
  }
  return sgn * I;
}

// Solve the (P + K + 1)-bordered system reusing the factored s.sys.
// When a floating contact's own flux derivative cancels >3 digits (e.g. heavily doped ohmic
// contacts where termC·col cancels bulk shift), substitutes the least-cancelling valid interior
// segment cut from s._currentCuts() (the same continuity identity impedance() uses), locked once
// per corrector or tangent solve in `cutLock`, and when `useCutRes` is true (corrector Newton)
// substitutes the matching nonlinear cut residual.
function solveBordered(s, ctx, work, rhs, pinRhs, circRhs, gRhs, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws, useCutRes = false, cutLock = null) {
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
  let dIdV = 0;
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
      // Check for catastrophic cancellation in the contact's own flux derivative; if >3 digits are
      // lost, read the steady current across a single locked interior segment cut instead.
      if (fk < 2 && ctx && ctx.cutEligible[fk]) {
        let bestCut = null;
        if (cutLock && cutLock.decided[fk]) {
          if (cutLock.seg[fk] >= 0) bestCut = getCut(s, ctx, cutLock.seg[fk]);
        } else {
          const C = s.termC[fk], colA = cols[a];
          let plainGross = Math.abs(s.termDI[fk]);
          for (let j = 0; j < N; j++) plainGross += Math.abs(C[j] * colA[j]);
          if (plainGross > 1e3 * Math.abs(A[a][a])) {
            const allCuts = ctx.cutsCache ?? (ctx.cutsCache = s._currentCuts());
            const lo = ctx.cutLo, hi = ctx.cutHi;
            const R = s.terms[fk].drive.R;
            const diCut = R > 0 ? 1 / R : 0;
            let bestGross = plainGross;
            for (let ci = 0; ci < allCuts.length; ci++) {
              const cut = allCuts[ci];
              if (fk === 0 ? !(cut.s + 1 <= lo) : !(cut.s >= hi)) continue;
              let cg = diCut;
              const idx = cut.idx, g = cut.g;
              for (let q = 0; q < idx.length; q++) {
                const o = idx[q];
                if (s.loc[o] >= 0) cg += Math.abs(g[q] * colA[s.rix[o]]);
              }
              if (cg < bestGross) { bestGross = cg; bestCut = cut; }
            }
            if (!(bestCut && bestGross < 0.01 * plainGross)) bestCut = null;
          }
          if (cutLock) {
            cutLock.decided[fk] = true;
            cutLock.seg[fk] = bestCut ? bestCut.s : -1;
          }
        }
        if (bestCut) {
          const sgn = fk === 0 ? 1 : -1;
          const R = s.terms[fk].drive.R;
          const diCut = R > 0 ? 1 / R : 0;
          for (let b = 0; b < S; b++) A[a][b] = (b < P ? 1 : -1) * dotCut(s, bestCut, sgn, cols[b]);
          A[a][a] += diCut;
          A[a][P + K] += circDs[a - P];
          let rhsCut = circRhs[a - P];
          if (useCutRes) {
            const d = s.terms[fk].drive;
            const src = s.sourceOverride.has(fk) ? s.sourceOverride.get(fk) : sourceAt(d.src, s.sourceTime, s.sourceBefore);
            const iCut = cutCurrent(s, bestCut.s, sgn);
            rhsCut = d.kind === 'I' ? iCut - src : iCut - (src - s.termV[fk]) / d.R;
          }
          A[a][S] = rhsCut - dotCut(s, bestCut, sgn, delta);
        }
      }
      if (ctx && fk === ctx.termIdx) dIdV = Math.abs(A[a][a]);
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
  return { ds: x[P + K], dIdV };
}

// GMRES(m) refinement of a bordered linear solve using difference-form J·v (_assembleDifference),
// preconditioned by the already-factored solveBordered with the exact same locked cut rows so
// operator, preconditioner and right-hand side belong to one consistent system.
function refineBordered(s, ctx, delta, deltaV, dsInit, wx, wV, ws, isFixed, cutLock = null) {
  const { termIdx, dpSpan, Bs, pinDs, circDs, work } = ctx;
  const N = s.sys.size, fl = s.floating, K = fl.length, pins = s.pins, P = pins.length, R = s.rix;
  const dim = N + K + (isFixed ? 0 : 1);
  const cutRows = new Array(K).fill(null);
  if (cutLock) {
    for (let a = 0; a < K; a++) {
      const fk = fl[a];
      if (fk < 2 && cutLock.decided[fk] && cutLock.seg[fk] >= 0) {
        const cut = getCut(s, ctx, cutLock.seg[fk]);
        if (cut) {
          const Rk = s.terms[fk].drive.R;
          cutRows[a] = { cut, sgn: fk === 0 ? 1 : -1, diCut: Rk > 0 ? 1 / Rk : 0 };
        }
      }
    }
  }
  const save = {
    res: Float64Array.from(s.res),
    termRes: Float64Array.from(s.termRes),
    termB: s.termB.map((b) => Float64Array.from(b)),
    termC: s.termC.map((c) => Float64Array.from(c)),
    termDI: Float64Array.from(s.termDI),
    termDX: s.termDX.map((x) => Float64Array.from(x)),
    transformed: s.transformed,
    pinRes: pins.map((p) => p.res),
  };
  const restore = () => {
    s.res.set(save.res);
    s.termRes.set(save.termRes);
    save.termB.forEach((b, k) => s.termB[k].set(b));
    save.termC.forEach((c, k) => s.termC[k].set(c));
    s.termDI.set(save.termDI);
    save.termDX.forEach((x, k) => s.termDX[k].set(x));
    s.transformed = save.transformed;
    pins.forEach((p, a) => (p.res = save.pinRes[a]));
  };
  const { rest, lin } = s._assembleDifference(Infinity);
  const BkU = fl.map((k) => Float64Array.from(s.termB[k]));
  const CkU = fl.map((k) => Float64Array.from(s.termC[k]));
  const DkU = fl.map((k) => s.termDI[k]);
  const DXU = s.termDX.map((x) => Float64Array.from(x));
  const BsU = new Float64Array(N + 1);
  const pinDsU = new Float64Array(P);
  const circDsU = new Float64Array(K);
  computeParamDs(s, termIdx, dpSpan, BsU, pinDsU, circDsU);
  const rowDiv = s.rowDiv;
  restore();

  const pRhsW = new Float64Array(Math.max(1, P));
  const cRhsW = new Float64Array(K);
  const rhsW = new Float64Array(N + 1);
  const outX = new Float64Array(N + 1);
  const outV = new Float64Array(K);
  const precond = (r) => {
    for (let j = 0; j < N; j++) rhsW[j] = r[j];
    for (let a = 0; a < P; a++) {
      const idx = R[pins[a].row];
      pRhsW[a] = rhsW[idx];
      rhsW[idx] = 0;
    }
    for (let a = 0; a < K; a++) cRhsW[a] = r[N + a];
    outV.fill(0);
    const sol = isFixed
      ? solveBordered(s, ctx, work, rhsW, pRhsW, cRhsW, 0, outX, outV, Bs, pinDs, circDs, work.zeroWx, work.zeroWv, 1, false, cutLock)
      : solveBordered(s, ctx, work, rhsW, pRhsW, cRhsW, r[N + K], outX, outV, Bs, pinDs, circDs, wx, wV, ws, false, cutLock);
    const z = new Float64Array(dim);
    for (let j = 0; j < N; j++) z[j] = outX[j];
    for (let a = 0; a < K; a++) z[N + a] = outV[a];
    if (!isFixed) z[N + K] = sol.ds;
    return z;
  };

  const v = new Float64Array(N + 1), jv = new Float64Array(N);
  const op = (x) => {
    for (let j = 0; j < N; j++) v[j] = x[j];
    v[N] = 0;
    rest.multiply(v, jv);
    if (save.transformed) s._untransform(v);
    lin.apply(v, jv, 1, 0);
    for (const p of pins) {
      if (!p.col) continue;
      for (const o of p.col) if (o !== p.row && s.loc[o] >= 0) jv[R[o]] -= x[R[p.row]];
    }
    const dsVal = isFixed ? 0 : x[N + K];
    const out = new Float64Array(dim);
    for (let j = 0; j < N; j++) {
      let t = jv[j] + BsU[j] * dsVal;
      for (let a = 0; a < K; a++) t += BkU[a][j] * x[N + a];
      out[j] = t / rowDiv[j];
    }
    for (let a = 0; a < K; a++) {
      const cr = cutRows[a];
      if (cr) {
        out[N + a] = cr.diCut * x[N + a] + circDsU[a] * dsVal + dotCut(s, cr.cut, cr.sgn, v);
      } else {
        let t = DkU[a] * x[N + a] + circDsU[a] * dsVal;
        for (let b = 0; b < K; b++) if (b !== a) t += DXU[fl[a]][fl[b]] * x[N + b];
        for (let j = 0; j < N; j++) t += CkU[a][j] * x[j];
        out[N + a] = t;
      }
    }
    for (let a = 0; a < P; a++) {
      const p = pins[a];
      let t = pinDsU[a] * dsVal;
      for (let j = 0; j < p.len; j++) t += p.w[j] * x[p.idx[j]];
      if (p.termSlope) for (let b = 0; b < K; b++) t += (p.termSlope.get(fl[b]) ?? 0) * x[N + b];
      out[R[p.row]] = t;
    }
    if (!isFixed) {
      let t = ws * dsVal;
      for (let j = 0; j < N; j++) if (wx[j] !== 0) t += wx[j] * x[j];
      for (let a = 0; a < K; a++) t += wV[a] * x[N + a];
      out[N + K] = t;
    }
    return precond(out);
  };

  // pb = precond(b) is the incoming plain solveBordered solution (delta, deltaV, dsInit).
  const pb = new Float64Array(dim);
  for (let j = 0; j < N; j++) pb[j] = delta[j];
  for (let a = 0; a < K; a++) pb[N + a] = deltaV[a];
  if (!isFixed) pb[N + K] = dsInit;
  const x = Float64Array.from(pb);

  const dot = (a, c) => { let t = 0; for (let j = 0; j < dim; j++) t += a[j] * c[j]; return t; };
  const bn = Math.sqrt(dot(pb, pb)) || 1;
  const m = 8, tol = 1e-10;
  try {
    const ax = op(x);
    const r0 = Float64Array.from(pb, (val, j) => val - ax[j]);
    const beta = Math.sqrt(dot(r0, r0));
    if (!(beta > tol * bn)) return { ok: true, ds: dsInit, moved: 0 };
    const V = [r0.map((val) => val / beta)];
    const H = Array.from({ length: m + 1 }, () => new Float64Array(m));
    const cs = [], sn = [], g = new Float64Array(m + 1);
    g[0] = beta;
    let k = 0, converged = false;
    for (; k < m; k++) {
      const w = op(V[k]);
      for (let i = 0; i <= k; i++) {
        const h = dot(V[i], w);
        H[i][k] = h;
        for (let j = 0; j < dim; j++) w[j] -= h * V[i][j];
      }
      const hn = Math.sqrt(dot(w, w));
      V.push(w.map((val) => val / (hn || 1)));
      for (let i = 0; i < k; i++) {
        const a = H[i][k], c = H[i + 1][k];
        H[i][k] = cs[i] * a + sn[i] * c;
        H[i + 1][k] = -sn[i] * a + cs[i] * c;
      }
      const t = Math.hypot(H[k][k], hn), c = t > 0 ? H[k][k] / t : 1, sv = t > 0 ? hn / t : 0;
      cs.push(c);
      sn.push(sv);
      H[k][k] = t;
      g[k + 1] = -sv * g[k];
      g[k] *= c;
      if (!(Math.abs(g[k + 1]) > tol * bn) || !(hn > 0)) { k++; converged = true; break; }
    }
    const y = new Float64Array(k);
    for (let i = k - 1; i >= 0; i--) {
      let val = g[i];
      for (let j = i + 1; j < k; j++) val -= H[i][j] * y[j];
      y[i] = val / H[i][i];
    }
    for (let i = 0; i < k; i++) for (let j = 0; j < dim; j++) x[j] += y[i] * V[i][j];
    if (!converged) return { ok: false, ds: dsInit, moved: Infinity };
    let maxDiff = 0;
    for (let j = 0; j < dim; j++) {
      if (!Number.isFinite(x[j])) return { ok: false, ds: dsInit, moved: Infinity };
      if (j < N + K) {
        const d = Math.abs(x[j] - pb[j]);
        if (d > maxDiff) maxDiff = d;
      }
    }
    for (let j = 0; j < N; j++) delta[j] = x[j];
    for (let a = 0; a < K; a++) deltaV[a] = x[N + a];
    return { ok: true, ds: isFixed ? 0 : x[N + K], moved: maxDiff };
  } catch {
    return { ok: false, ds: dsInit, moved: Infinity };
  }
}

// Compute unit tangent (tx, tV, ts) at current converged state + sVal, oriented as a whole by
// the sign of its inner product with (prevTx, prevTV, prevTs) in the trace's metric, or with
// ts > 0 when prevTx is null.
function computeTangent(s, ctx, sVal, prevTx, prevTV, prevTs, tx, tV, reuseFactor = false, refine = false) {
  const { termIdx, p0, dpSpan, potIdx, invN, wsEff, Bs, pinDs, circDs, wx, wV, work } = ctx;
  const N = s.sys.size, K = s.floating.length;
  if (!reuseFactor) {
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    s.islandsOn = false;
    s.assemble(Infinity);
    ctx.cutsCache = null;
    equilibrate(s);
    computeParamDs(s, termIdx, dpSpan, Bs, pinDs, circDs);
    s._factor();
  }
  let ws;
  if (prevTx) {
    ws = buildWeights(potIdx, invN, wsEff, s.VT, prevTx, prevTV, prevTs, wx, wV);
  } else {
    wx.fill(0);
    wV.fill(0);
    ws = 1;
  }
  const cutLock = { decided: [false, false], seg: [-1, -1] };
  const res = solveBordered(s, ctx, work, work.zero, work.zeroPin, work.zeroCirc, 1, tx, tV, Bs, pinDs, circDs, wx, wV, ws, false, cutLock);
  let tsRaw = res.ds;
  if (refine || s.sys.perturbedPivots > 0 || s.sys.lostColumns > 0) {
    const ref = refineBordered(s, ctx, tx, tV, tsRaw, wx, wV, ws, false, cutLock);
    if (ref.ok) tsRaw = ref.ds;
  }
  if (!Number.isFinite(tsRaw)) return { ok: false };
  const nrm = tangentNorm(potIdx, invN, wsEff, s.VT, tx, tV, tsRaw);
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
    const dotPrev = tangentInner(potIdx, invN, wsEff, s.VT, prevTx, prevTV, prevTs, tx, tV, ts);
    if (dotPrev < 0) {
      for (let j = 0; j < N; j++) tx[j] = -tx[j];
      for (let a = 0; a < K; a++) tV[a] = -tV[a];
      ts = -ts;
    }
    cosTheta = Math.abs(dotPrev);
  }
  return { ok: true, ts, cosTheta };
}

function kirchhoffSum(s) {
  let sum = 0;
  for (const I of s.termI) sum += I;
  return Math.abs(sum);
}

// Corrector Newton on F(x, V_fl, s) = 0 and wx·(x − xPred) + wV·(V − vPred) + ws·(s − sPred) = 0.
function correctorNewton(s, ctx, sStart, uPred, uLoPred, vPred, sPred, wx, wV, ws, { maxIter = 40, tol = BRANCH_TOL, maxStep = 10, fixedS = false } = {}) {
  const { termIdx, p0, dpSpan, potIdx, wsEff, Bs, pinDs, circDs, work, budget } = ctx;
  budget?.check();
  const fl = s.floating, K = fl.length;
  const delta = work.delta, deltaV = work.deltaV;
  const pinRhs = work.pinRhs, circRhs = work.circRhs;
  let sVal = fixedS ? sPred : sStart;

  const runLoop = (islandsOn, isFixed, limit) => {
    const { isCurrentDrive } = ctx;
    const cutLock = { decided: [false, false], seg: [-1, -1] };
    let lastStep = Infinity;
    let refining = false;
    let refineFails = 0;
    let usedRefine = false;
    let eqFailed = false;
    for (let it = 1; it <= limit; it++) {
      budget?.check();
      s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
      s.islandsOn = islandsOn;
      try {
        s.assemble(Infinity);
        ctx.cutsCache = null;
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
          solved = solveBordered(s, ctx, work, s.res, pinRhs, circRhs, 0, delta, deltaV, Bs, pinDs, circDs, work.zeroWx, work.zeroWv, 1, true, cutLock);
        } else {
          solved = solveBordered(s, ctx, work, s.res, pinRhs, circRhs, g, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws, true, cutLock);
        }
      } catch (err) {
        if (!(err instanceof SolverError)) throw err;
        return { converged: false, iterations: it, s: sVal, error: err.message };
      }
      let ds = isFixed ? 0 : solved.ds;
      if (s.sys.perturbedPivots > 0 || s.sys.lostColumns > 0) refining = true;
      if (refining && refineFails < 2) {
        const ref = isFixed
          ? refineBordered(s, ctx, delta, deltaV, 0, work.zeroWx, work.zeroWv, 1, true, cutLock)
          : refineBordered(s, ctx, delta, deltaV, ds, wx, wV, ws, false, cutLock);
        if (ref.ok && ref.moved < 10) {
          ds = isFixed ? 0 : ref.ds;
          usedRefine = true;
        } else {
          refineFails++;
        }
      }
      clipLevels(s, delta, maxStep);
      const stateStep = maxPotentialStep(s, delta, deltaV);
      const stateOwn = maxPotentialStep(s, delta, deltaV, false);
      const step = Math.max(stateStep, wsEff * Math.abs(ds));
      const own = Math.max(stateOwn, wsEff * Math.abs(ds));
      if (!Number.isFinite(step) || !Number.isFinite(ds)) return { converged: false, iterations: it, s: sVal, error: 'non-finite update' };
      if (it > 1 && own > 1e5) return { converged: false, iterations: it, s: sVal, error: 'diverging' };
      const alpha = own > maxStep ? maxStep / own : 1;
      addToState(s, delta, -alpha);
      for (let a = 0; a < K; a++) s.termV[fl[a]] -= alpha * deltaV[a];
      if (!isFixed) sVal -= alpha * ds;
      budget?.markAttempted();
      const pAfter = p0 + sVal * dpSpan;
      const pScale = isCurrentDrive
        ? Math.max(Math.abs(pAfter), solved.dIdV * s.VT, 1e-300)
        : Math.max(Math.abs(pAfter), s.VT);
      const paramRel = isFixed ? 0 : Math.abs(ds * dpSpan) / pScale;
      const convMetric = Math.max(stateStep, 0.1 * paramRel);
      if (alpha === 1 && stateStep < tol && paramRel <= 1e-9) {
        const checkEq = fl.includes(termIdx) && (s.islands.length === 0 || islandsOn);
        try {
          s.sourceOverride.set(termIdx, pAfter);
          s.computeConcentrations();
          if (checkEq) s._assembleBookkeeping(Infinity);
        } catch (err) {
          if (!(err instanceof SolverError)) throw err;
          return { converged: false, iterations: it, s: sVal, modelError: err.message };
        }
        let eqOk = true;
        if (checkEq) {
          const d = s.terms[termIdx].drive;
          const iTarget = d.kind === 'I' ? pAfter : (pAfter - s.termV[termIdx]) / d.R;
          const iRo = d.kind === 'I' ? 0 : (4 * Number.EPSILON * Math.max(Math.abs(pAfter), Math.abs(s.termV[termIdx]))) / Math.abs(d.R);
          const iScale = Math.max(Math.abs(iTarget), solved.dIdV * s.VT, 1e-300);
          const hasCut = termIdx < 2 && cutLock.decided[termIdx] && cutLock.seg[termIdx] >= 0;
          if (hasCut) {
            const iCut = cutCurrent(s, cutLock.seg[termIdx], termIdx === 0 ? 1 : -1);
            if (Math.abs(iCut - iTarget) > 1e-8 * Math.max(iScale, Math.abs(iCut)) + 10 * iRo) eqOk = false;
          } else {
            const iRead = contactCurrent(s, termIdx);
            if (Math.abs(iRead - iTarget) > 1e-5 * Math.max(iScale, Math.abs(iRead)) + 10 * iRo) eqOk = false;
          }
          if (eqOk && s.terms.length === 2 && s.islands.length === 0) {
            const iOther = contactCurrent(s, 1 - termIdx);
            const kTol = 1e-5 * Math.max(iScale, Math.abs(iOther)) + 10 * iRo;
            if (Math.abs(iOther + iTarget) > kTol) eqOk = false;
          }
        }
        if (!eqOk) {
          eqFailed = true;
          refining = true;
          lastStep = convMetric;
          continue;
        }
        return {
          converged: true,
          iterations: it,
          s: sVal,
          reusable: !checkEq && !islandsOn && s.sys.perturbedPivots === 0 && s.sys.lostColumns === 0,
          refined: usedRefine,
        };
      }
      if (!refining && it >= 3 && alpha === 1 && convMetric > 0.25 * lastStep && convMetric < 1e-3) refining = true;
      lastStep = convMetric;
    }
    return { converged: false, iterations: limit, s: sVal, error: eqFailed ? 'reported terminal equation not satisfied' : 'Newton did not converge' };
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
      return { converged: true, iterations: Math.max(r1.iterations, 6), s: sVal, reusable: false };
    }
    s._restore(preIslSnap);
    s.sourceOverride.set(termIdx, p0 + sVal * dpSpan);
    const qPol = s._steadySteps(budget ? { maxSteps: 1, budget } : { maxSteps: 1 });
    if (qPol.converged && okKirchhoff()) {
      return { converged: true, iterations: Math.max(r1.iterations, 6), s: sVal, reusable: false };
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

// Unique modes (im >= 0), collapsing duplicate copies returned across multiple Arnoldi shifts.
function uniqueModes(st) {
  if (!st || !st.modes) return [];
  const out = [];
  for (const m of st.modes) {
    if (m.im < 0) continue;
    let dupIdx = -1;
    for (let k = 0; k < out.length; k++) {
      const o = out[k];
      const sc = Math.max(Math.hypot(o.re, o.im), Math.hypot(m.re, m.im), 1e-12);
      if (Math.hypot(o.re - m.re, o.im - m.im) <= 1e-3 * sc) {
        const ov = vecOverlap(o, m);
        if (ov < 0 || ov >= 0.98) { dupIdx = k; break; }
      }
    }
    if (dupIdx >= 0) {
      if (m.error < out[dupIdx].error) out[dupIdx] = m;
    } else {
      out.push(m);
    }
  }
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
    if (Number.isFinite(ov) && ov <= 4) return Math.min(1, ov);
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

// Count verified strictly growing eigenvalues (complex pairs counted as 2, real modes as 1).
function countGrowing(st, kBar, foldBelow = 0) {
  const modes = uniqueModes(st);
  // At a located fold, J is singular: one real mode is passing through zero, and the fold locator
  // leaves it at a small Re λ of either sign. It isn't counted. It is the real mode nearest zero,
  // and it must lie nearer zero than the nearest real mode at either end of the step (`foldBelow`),
  // since it moves through zero across the step.
  let foldIdx = -1;
  if (foldBelow > 0) {
    let best = foldBelow;
    for (let i = 0; i < modes.length; i++) {
      const m = modes[i];
      if (!(m.im > kBar * m.error) && Math.abs(m.re) < best) { best = Math.abs(m.re); foldIdx = i; }
    }
  }
  let count = 0;
  for (let i = 0; i < modes.length; i++) {
    if (i === foldIdx) continue;
    const m = modes[i];
    if (m.re > kBar * m.error) count += m.im > kBar * m.error ? 2 : 1;
  }
  return count;
}

// |Re λ| of the real mode nearest zero (Infinity if none listed).
function nearestRealAbs(st, kBar) {
  let best = Infinity;
  for (const m of uniqueModes(st)) if (!(m.im > kBar * m.error)) best = Math.min(best, Math.abs(m.re));
  return best;
}

// Nearer-end guard across all matched oscillatory pairs (im > kBar * error):
//   |Re λ_B − Re λ_A| ≤ ALPHA_RE · max(min(|Re λ_A|, |Re λ_B|), kBar · err)
// on the same side of Re λ = 0. Real modes deep in the LHP that merely scale exponentially with
// Boltzmann depletion are not braked here; real crossings are caught by folds and index refinement.
function modeStepSmooth(stPrev, stNew, kBar, stats) {
  if (!stPrev || !stNew) return { ok: true };
  const { pairs } = matchModes(stPrev, stNew);
  for (const { a, b } of pairs) {
    if (!(a.im > kBar * a.error || b.im > kBar * b.error)) continue;
    if ((a.re < 0) !== (b.re < 0)) continue;
    const err = Math.max(a.error, b.error);
    const omega = Math.max(a.im, b.im);
    const dNear = Math.max(Math.min(Math.abs(a.re), Math.abs(b.re)), kBar * err, 1e-3 * omega);
    if (Math.abs(b.re - a.re) > ALPHA_RE * dNear) {
      if (stats) stats.modeNearEnd++;
      return { ok: false, reason: 'oscillatory Re λ moved past nearer-end scale' };
    }
  }
  return { ok: true };
}

// Chord-hyperplane solve at fraction theta in (0, 1) between two converged states A and B.
function solveOnChord(s, ctx, A, B, theta) {
  const { potIdx, invN, wsEff, wx, wV, uPred, uLoPred, vPred } = ctx;
  const N = s.sys.size, fl = s.floating, K = fl.length;
  const chX = ctx.work.chordX, chV = ctx.work.chordV;
  chX.fill(0);
  for (let i = 0; i < potIdx.length; i++) {
    const j = potIdx[i], o = s.fullOf[j];
    chX[j] = (B.snap.u[o] - A.snap.u[o]) + (B.snap.uLo[o] - A.snap.uLo[o]);
  }
  for (let a = 0; a < K; a++) chV[a] = B.snap.termV[fl[a]] - A.snap.termV[fl[a]];
  const chS = B.s - A.s;
  const nrm = tangentNorm(potIdx, invN, wsEff, s.VT, chX, chV, chS);
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
  const ws = buildWeights(potIdx, invN, wsEff, s.VT, chX, chV, chSu, wx, wV);
  const cor = correctorNewton(s, ctx, sPred, uPred, uLoPred, vPred, sPred, wx, wV, ws);
  if (!cor.converged) return { ok: false, chordLen: nrm };
  const tan = computeTangent(s, ctx, cor.s, A.tx, A.tV, A.ts, ctx.work.locTx, ctx.work.locTV, !!cor.reusable, !!cor.refined);
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
  const floor = driveRoundoffFloor(s, ctx.termIdx, drive);
  const located = converged && solvedAny;
  const note = located ? undefined : (failNote ?? 'fold location did not reach tolerance (using best bracket point)');
  return { drive, terminals, s: best.s, state: best, floor, located, note };
}

function findComplexPair(st, omegaEst, refMode = null, kBar = 10) {
  let bestVec = null, bestOv = 0.45;
  if (refMode && refMode._vr) {
    for (const m of st.modes) {
      if (!(m.im > kBar * m.error)) continue;
      const ov = vecOverlap(refMode, m);
      if (ov > bestOv) { bestOv = ov; bestVec = m; }
    }
    if (bestVec) return bestVec;
  }
  let best = null, bestDist = Infinity;
  for (const m of st.modes) {
    if (!(m.im > kBar * m.error)) continue;
    const dist = Math.abs(m.im - omegaEst);
    if (dist < bestDist) { bestDist = dist; best = m; }
  }
  return best;
}

function evalHopfMode(s, ctx, M, omegaEst, refMode = null, extraOmega = 0) {
  const drive = ctx.p0 + M.s * ctx.dpSpan;
  s._restore(M.snap);
  const kBar = ctx.stabOpts.kBar ?? 10;
  const fHz = omegaEst / (2 * Math.PI);
  const baseOpts = ctx.stabOpts;
  const freqs = baseOpts.frequencies ? [...baseOpts.frequencies, fHz] : [fHz];
  if (extraOmega > 0 && Math.abs(extraOmega - omegaEst) > 0.1 * Math.max(extraOmega, omegaEst)) {
    freqs.push(extraOmega / (2 * Math.PI));
  }
  const fastOpts = { ...baseOpts, shifts: [Math.max(1, 0.25 * omegaEst)], frequencies: freqs };
  let st = evaluateStability(s, ctx.termIdx, drive, fastOpts, ctx.isCurrentDrive, ctx.budget);
  let pair = findComplexPair(st, omegaEst, refMode, kBar);
  if (!pair) {
    st = evaluateStability(s, ctx.termIdx, drive, baseOpts, ctx.isCurrentDrive, ctx.budget);
    pair = findComplexPair(st, omegaEst, refMode, kBar);
  }
  return pair;
}

// Locate a Hopf bifurcation (Re λ = 0 for complex conjugate pair with Im λ > kBar * error) between
// A0 and B0, using purely relative drive and chord tolerances.
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
  let sawHidden = false;
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
    if (A0.ts * B0.ts > 0 && M.ts * A0.ts <= 0) {
      sawHidden = true;
      failNote = 'fold encountered during Hopf bisection (using best bracket point)';
      break;
    }
    const refM = theta <= 0.5 ? mA : mB;
    const omegaTry = (1 - theta) * mA.im + theta * mB.im;
    const mM = evalHopfMode(s, ctx, M, omegaTry > 0 ? omegaTry : omega, refM, omega);
    if (!mM) {
      sawHidden = true;
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
  return { drive, terminals, omega: best.mode.im, located, sawHidden, note };
}

function detectHopfBracket(stA, stB, kBar, stats = null) {
  const sig = (m) => (m.re > kBar * m.error ? 1 : m.re < -kBar * m.error ? -1 : 0);
  const hasIm = (m) => m.im > kBar * m.error;
  const matched = matchModes(stA, stB);
  for (const { a: ma, b: mb } of matched.pairs) {
    if (!(hasIm(ma) && hasIm(mb))) continue;
    const sa = sig(ma), sb = sig(mb);
    if (sa !== 0 && sb !== 0 && sa !== sb) {
      if (stats) stats.hopfByVecMatch++;
      return { pairA: ma, pairB: mb, omega: 0.5 * (ma.im + mb.im) };
    }
  }
  return null;
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
        rejectAngle: 0,
        modeNearEnd: 0,
        refineFlip: 0,
        refineFoldIndex: 0,
        roundoffFoldSuppressed: 0,
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
    let pendingFold = null;
    // A provisional fold still pending when the trace ends (maxSteps, voltageBound, budgetMs) is
    // reported when its located drive lies beyond BOTH of its step's accepted drives by more than
    // the round-off floor: the accepted points then bracket a turnover, though no later point
    // retreated to confirm it. It carries a note saying so, and goes in event order. Otherwise it
    // is dropped.
    const emitPendingAtEnd = () => {
      const pend = pendingFold;
      pendingFold = null;
      if (!pend || !pend.strict) return;
      const pf = pend.preFold;
      const ev = {
        kind: 'fold',
        drive: pf.drive,
        terminals: pf.terminals,
        between: pend.between,
        note: 'the trace ended before a later point confirmed the turn',
      };
      let at = events.findIndex((e) => e.between[0] >= pend.between[0]);
      if (at < 0) at = events.length;
      events.splice(at, 0, ev);
    };

    const isDecided = (v) => v === 'no growing mode found' || v === 'unstable';

    // Post-trace consistency audit between geometric events and stability marks.
    const runAudit = (allowLocate) => {
      if (!doStab || points.length < 2) return;

      let prevDec = -1;
      for (let i = 0; i < points.length; i++) {
        const v = points[i].stability?.verdict;
        if (!isDecided(v)) continue;
        if (prevDec >= 0 && points[prevDec].stability.verdict !== v) {
          const covered = events.some((ev) => {
            if (ev.between[0] < prevDec || ev.between[1] > i) return false;
            if (ev.kind === 'hopf') return true;
            if (ev.kind === 'fold') {
              const vL = points[ev.between[0]]?.stability?.verdict;
              const vR = points[ev.between[1]]?.stability?.verdict;
              return vL !== 'undecided' && vR !== 'undecided';
            }
            return false;
          });
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
      const potIdx = stateMetricIndices(this);
      const invN = potIdx.length > 0 ? 1 / potIdx.length : 1;
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
      const chargedPort = (p) => p.capacitance || p.reactions.length > 0 || this.model.species.some((_, i) => this.z[i] !== 0 && p.species[i].type !== 'blocked');
      const portWindows = this.model.ports.filter(chargedPort).flatMap((p) => Array.from(p.nodes));
      const cutLo = portWindows.length ? Math.min(...portWindows) : Infinity;
      const cutHi = portWindows.length ? Math.max(...portWindows) : -Infinity;
      const cutEligible = [
        !this.constraints.some((cs) => cs.circuit === 0),
        !this.constraints.some((cs) => cs.circuit === 1),
      ];
      ctx = {
        termIdx,
        p0,
        dpSpan,
        potIdx,
        invN,
        wsEff: 1,
        cutLo,
        cutHi,
        cutEligible,
        cutsCache: null,
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

      // 1. Initial unit tangent at s = 0; scale wsEff = max(1, stateSlope) so active devices start
      //    at a 50/50 state/drive balance while insensitive terminals retain dimensionless s scaling.
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
      if (Math.abs(tan0.ts) > 0) {
        const stateSlope = stateOnlyNorm(potIdx, invN, this.VT, tx, tV) / Math.abs(tan0.ts);
        if (stateSlope > 1 && stateSlope < 1e100 && Number.isFinite(stateSlope)) {
          ctx.wsEff = stateSlope;
          tan0 = computeTangent(this, ctx, 0, null, null, null, tx, tV, true);
        }
      }

      let ts = tan0.ts;
      let sCur = 0;

      // Rebalance `wsEff` whenever the drive component drifts above 0.85 of the unit tangent while
      // the state responds, restoring 50/50 balance so opposite-state branches never look parallel.
      const rebalanceMetric = () => {
        if (!(ctx.wsEff * Math.abs(ts) > 0.8)) return 1;
        const stNorm = stateOnlyNorm(potIdx, invN, this.VT, tx, tV);
        if (!(stNorm > 1e-12 && Math.abs(ts) > 0)) return 1;
        const nextWs = Math.min(1e100, Math.max(1, stNorm / Math.abs(ts)));
        if (nextWs >= 0.99 * ctx.wsEff) return 1;
        ctx.wsEff = nextWs;
        const newNrm = tangentNorm(potIdx, invN, ctx.wsEff, this.VT, tx, tV, ts);
        for (let j = 0; j < N; j++) tx[j] /= newNrm;
        for (let a = 0; a < K; a++) tV[a] /= newNrm;
        ts /= newNrm;
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

      // Initial stride in the balanced metric.
      const maxPotPerStep = 6;
      const maxSPerStep = step !== undefined ? step / absSpan : 0.2;
      let dsArc = 2;
      if (Math.abs(ts) > 0) dsArc = Math.min(dsArc, maxSPerStep / Math.abs(ts));
      const initStateRate = maxStateStep(this, tx, tV);
      if (initStateRate > 1e-12) dsArc = Math.min(dsArc, 5 / initStateRate);
      dsArc = Math.max(MIN_ARC, Math.min(MAX_ARC, dsArc));

      let consecRefine = 0;
      let flippedStart = false;
      let stepCount = 0;
      let lastDriveDir = Math.sign(dpSpan);
      pendingFold = null;

      // Resolve any provisional fold against an accepted point's drive `pAcc` before recording
      // events from the new step. A provisional fold stands only if a subsequent accepted drive
      // retreats beyond the drive round-off floor; advancing forward drops it.
      const settlePendingFold = (pAcc, floorAcc) => {
        if (!pendingFold) return;
        const fTol = Math.max(pendingFold.floor, floorAcc);
        const delta = pendingFold.dir * (pAcc - pendingFold.peakDrive);
        if (delta < -fTol) {
          const pf = pendingFold.preFold;
          events.push({
            kind: 'fold',
            drive: pf.drive,
            terminals: pf.terminals,
            between: pendingFold.between,
            ...(pf.located ? {} : { located: false, note: pf.note }),
          });
          lastDriveDir = -pendingFold.dir;
          pendingFold = null;
        } else if (delta > fTol) {
          pendingFold = null;
        }
      };

      // Evaluate stability on a candidate step and verify both the oscillatory nearer-end rule and
      // Morse-index consistency (rejecting to halve when an unexplained flip, undecided mark,
      // fold-adjacent Hopf, or hidden fold/coalescence inside a Hopf bracket is detected, up to
      // MAX_REFINE halvings).
      const checkStepStability = (stPrev, candDrive, candStateObj, preFold, dsStep) => {
        budget?.check();
        const stCand = evaluateStability(this, termIdx, candDrive, stabOpts, isCurrentDrive, budget);
        const sm = modeStepSmooth(stPrev, stCand, kBar, _stats);
        if (!sm.ok) return { ok: false, reason: sm.reason, isRefine: false, stCand };

        let preHopf = null;
        const canRefine = consecRefine < MAX_REFINE && dsStep > MIN_ARC * 2;
        if (canRefine && stPrev) {
          const prevUndec = stPrev.verdict === 'undecided';
          const candUndec = stCand.verdict === 'undecided';
          if (preFold) {
            if (prevUndec || candUndec) {
              if (_stats) _stats.refineFoldIndex++;
              return { ok: false, reason: 'fold step has undecided stability mark', isRefine: true, stCand };
            }
            if (isDecided(stPrev.verdict) && isDecided(stCand.verdict)) {
              // Simple fold invariant: |nB - nA| === 1 AND N_+(fold) === min(nA, nB).
              // Violated whenever a Hopf sits just before or after the fold!
              const nA = countGrowing(stPrev, kBar);
              const nB = countGrowing(stCand, kBar);
              this._restore(preFold.state.snap);
              const stFold = evaluateStability(this, termIdx, preFold.drive, stabOpts, isCurrentDrive, budget);
              this._restore(candStateObj.snap);
              if (stFold.verdict === 'undecided') {
                if (_stats) _stats.refineFoldIndex++;
                return { ok: false, reason: 'fold state has undecided stability mark', isRefine: true, stCand };
              }
              const nF = countGrowing(stFold, kBar, Math.min(nearestRealAbs(stPrev, kBar), nearestRealAbs(stCand, kBar)));
              if (Math.abs(nB - nA) !== 1 || nF !== Math.min(nA, nB)) {
                if (_stats) _stats.refineFoldIndex++;
                return { ok: false, reason: 'fold step contains additional stability crossing', isRefine: true, stCand };
              }
            }
          } else {
            const nA = countGrowing(stPrev, kBar);
            const nB = countGrowing(stCand, kBar);
            const hopfBracket = detectHopfBracket(stPrev, stCand, kBar);
            if (hopfBracket) {
              if (prevUndec || candUndec) {
                if (_stats) _stats.refineFlip++;
                return { ok: false, reason: 'Hopf step has undecided stability mark', isRefine: true, stCand };
              }
              const prevStateObj = states[points.length - 1];
              const hopf = locateHopf(this, ctx, prevStateObj, candStateObj, hopfBracket.pairA, hopfBracket.pairB, hopfBracket.omega);
              this._restore(candStateObj.snap);
              if (!hopf || hopf.sawHidden) {
                if (_stats) _stats.refineFlip++;
                return { ok: false, reason: 'Hopf bracket contains hidden fold or coalescence', isRefine: true, stCand };
              }
              preHopf = hopf;
            } else if (nA !== nB) {
              if (_stats) _stats.refineFlip++;
              return { ok: false, reason: 'unexplained stability flip', isRefine: true, stCand };
            }
          }
        }
        return { ok: true, stCand, preHopf };
      };

      const acceptPoint = (sNew, newTx, newTV, newTs, preSt = null, preFold = null, preHopf = null, foldConfirm = 'none') => {
        const kPrev = points.length - 1;
        const prevState = states[kPrev];
        const pPrevAcc = points[kPrev].drive;
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
        const floorNew = driveRoundoffFloor(this, termIdx, driveNew);
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

        settlePendingFold(driveNew, floorNew);

        if (preFold) {
          if (foldConfirm === 'confirm') {
            events.push({
              kind: 'fold',
              drive: preFold.drive,
              terminals: preFold.terminals,
              between: [kPrev, kNew],
              ...(preFold.located ? {} : { located: false, note: preFold.note }),
            });
            lastDriveDir = -lastDriveDir;
            pendingFold = null;
          } else if (foldConfirm === 'provisional') {
            const peakDrive = lastDriveDir * (driveNew - pPrevAcc) >= 0 ? driveNew : pPrevAcc;
            pendingFold = {
              preFold,
              dir: lastDriveDir,
              peakDrive,
              between: [kPrev, kNew],
              floor: Math.max(preFold.floor, floorNew),
              strict: !!preFold.located
                && lastDriveDir * (preFold.drive - pPrevAcc) > Math.max(preFold.floor, floorNew)
                && lastDriveDir * (preFold.drive - driveNew) > Math.max(preFold.floor, floorNew),
            };
          }
        }

        if (preHopf) {
          if (_stats) _stats.hopfByVecMatch++;
          events.push({
            kind: 'hopf',
            drive: preHopf.drive,
            terminals: preHopf.terminals,
            omega: preHopf.omega,
            between: [kPrev, kNew],
            ...(preHopf.located ? {} : { located: false, note: preHopf.note }),
          });
        } else if (doStab && stFullList[kPrev] && stNew) {
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
        ctx.wx.fill(0);
        ctx.wV.fill(0);
        const corLand = correctorNewton(this, ctx, sTarget, ctx.uPred, ctx.uLoPred, ctx.vPred, sTarget, ctx.wx, ctx.wV, 1, { fixedS: true });
        if (!corLand.converged) return null;

        for (let a = 0; a < K; a++) {
          if (Math.abs(this.termV[fl[a]] - vStart[a]) > voltageBound) return null;
        }
        const lastPt = states[points.length - 1];
        if (step !== undefined && Math.abs((sTarget - lastPt.s) * dpSpan) > step * (1 + 1e-6)) return null;

        let tanLand = computeTangent(this, ctx, sTarget, refTx, refTV, refTs, work.landTx, work.landTV, !!corLand.reusable, !!corLand.refined);
        if (!tanLand.ok && !corLand.refined && (this.pins.length > 0 || this.islands.length > 0)) {
          try {
            const tanRef = computeTangent(this, ctx, sTarget, refTx, refTV, refTs, work.landTx, work.landTV, true, true);
            if (tanRef.ok) tanLand = tanRef;
          } catch (err) {
            if (!(err instanceof Error)) throw err;
          }
        }
        if (!tanLand.ok) return null;

        return { s: sTarget, tx: Float64Array.from(work.landTx), tV: Float64Array.from(work.landTV), ts: tanLand.ts, tan: tanLand, snap: this._snapshot() };
      };

      // 2. Main pseudo-arclength continuation loop.
      while (stepCount < maxSteps) {
        budget?.check();

        let dsStep = dsArc;
        if (Math.abs(ts) > 0) dsStep = Math.min(dsStep, maxSPerStep / Math.abs(ts));
        const stateRate = maxStateStep(this, tx, tV);
        if (stateRate > 1e-12) dsStep = Math.min(dsStep, maxPotPerStep / stateRate);
        if (dsStep < MIN_ARC) dsStep = MIN_ARC;

        if (sCur < 1 && ts > 0 && sCur + dsStep * ts > 1.02) {
          dsStep = Math.max(MIN_ARC, (1.01 - sCur) / ts);
        }

        const beforeSnap = this._snapshot();
        const prevStateObj = states[points.length - 1];
        const sPred = sCur + dsStep * ts;

        addToState(this, tx, dsStep);
        for (let a = 0; a < K; a++) this.termV[fl[a]] += dsStep * tV[a];
        ctx.uPred.set(this.u);
        ctx.uLoPred.set(this.uLo);
        for (let a = 0; a < K; a++) ctx.vPred[a] = this.termV[fl[a]];
        const ws = buildWeights(potIdx, invN, ctx.wsEff, this.VT, tx, tV, ts, ctx.wx, ctx.wV);

        const cor = correctorNewton(this, ctx, sPred, ctx.uPred, ctx.uLoPred, ctx.vPred, sPred, ctx.wx, ctx.wV, ws, {
          maxIter: 35,
        });

        const rejectStep = (factor, reasonIfMin, isModeOnly = false) => {
          this._restore(beforeSnap);
          if (!isModeOnly) consecRefine = 0;
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

        // Orient the first accepted step by its actual drive, not noisy `ts` at s = 0; and never
        // step out of the span (`sNext < 0`) from `sCur <= 0`.
        if (stepCount === 0 && sNext < 0 && !flippedStart) {
          this._restore(beforeSnap);
          for (let j = 0; j < N; j++) tx[j] = -tx[j];
          for (let a = 0; a < K; a++) tV[a] = -tV[a];
          ts = -ts;
          states[0].tx = Float64Array.from(tx);
          states[0].tV = Float64Array.from(tV);
          states[0].ts = ts;
          flippedStart = true;
          continue;
        }
        if (sCur <= 0 && sNext < 0) {
          if (!rejectStep(0.5, `continuation stopped: could not advance into span from ${formatDrive(p0, kind, hasGeometry)}`)) break;
          continue;
        }

        // Guard 1: Hard `step` cap and pointwise state-step bound (potentials + 10× gates).
        const actualDp = Math.abs((sNext - sCur) * dpSpan);
        const stepExceeded = step !== undefined && actualDp > step * (1 + 1e-6);
        const actualStateMove = maxStateDiff(this, beforeSnap, this.u, this.uLo, this.termV);
        if (stepExceeded || actualStateMove > maxPotPerStep * 1.25) {
          const factor = stepExceeded && actualDp > 0 ? Math.min(0.5, 0.8 * (step / actualDp)) : 0.5;
          if (!rejectStep(factor, `continuation stopped: step control could not bound state step near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`)) break;
          continue;
        }

        // Guard 2: Tangent continuity in the balanced metric.
        prevTx.set(tx);
        prevTV.set(tV);
        const prevTs = ts;
        let tanNext;
        try {
          tanNext = computeTangent(this, ctx, sNext, prevTx, prevTV, prevTs, tx, tV, !!cor.reusable, !!cor.refined);
        } catch (err) {
          if (!(err instanceof Error)) throw err;
          tanNext = { ok: false };
        }
        const badTan = !tanNext.ok || stepHasHiddenS(this, ctx, beforeSnap, sCur, prevTx, prevTV, prevTs, this.u, this.uLo, this.termV, sNext, tx, tV, tanNext.ts);
        if (badTan) {
          if (_stats) _stats.rejectAngle++;
          tx.set(prevTx);
          tV.set(prevTV);
          ts = prevTs;
          if (!rejectStep(0.5, `continuation stopped: tangent turned sharply or became singular near ${formatDrive(p0 + sCur * dpSpan, kind, hasGeometry)}`)) break;
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

        const candState = {
          snap: this._snapshot(),
          s: sNext,
          tx: Float64Array.from(tx),
          tV: Float64Array.from(tV),
          ts: tanNext.ts,
        };
        const pPrev = p0 + sCur * dpSpan;
        const pCand = p0 + sNext * dpSpan;
        const candFloor = driveRoundoffFloor(this, termIdx, pCand);

        // Effective travel direction entering this step (accounting for any pending fold that
        // pCand would confirm on acceptance, without mutating state until accepted).
        let effDir = lastDriveDir;
        if (pendingFold && pendingFold.dir * (pCand - pendingFold.peakDrive) < -Math.max(pendingFold.floor, candFloor)) {
          effDir = -pendingFold.dir;
        }
        const dirSign = effDir * dpSpan > 0 ? 1 : -1;

        // If a provisional fold is waiting for confirmation and `ts` flips back without pCand
        // retreating yet, check whether this step straddled a second narrow fold (an S-curve
        // middle branch narrower than one stride); if so, halve to resolve the middle branch.
        if (pendingFold && effDir === lastDriveDir && prevTs * dirSign < 0 && tanNext.ts * dirSign > 0 && consecRefine < MAX_REFINE && dsStep > MIN_ARC * 2) {
          const fBack = locateFold(this, ctx, prevStateObj, candState);
          this._restore(candState.snap);
          const backTol = Math.max(fBack.floor, pendingFold.floor, candFloor);
          if ((-pendingFold.dir) * (fBack.drive - pPrev) > backTol || (-pendingFold.dir) * (fBack.drive - pCand) > backTol) {
            consecRefine++;
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            if (!rejectStep(0.5, `continuation stopped: could not resolve narrow fold pair near ${formatDrive(pPrev, kind, hasGeometry)}`, true)) break;
            continue;
          }
        }

        // Check for fold FIRST so the One Span Rule handles folds past `to` and immediate return
        // past `p0` uniformly. A fold is located when `ts` turns against the current travel
        // direction and accepted only when accepted drives retreat beyond the drive round-off
        // floor (immediately if this step itself retreats, otherwise provisionally until a
        // subsequent accepted point retreats).
        let stepFold = prevTs * dirSign > 0 && tanNext.ts * dirSign < 0;
        let preFold = null;
        let foldConfirm = 'none';
        if (stepFold) {
          try {
            preFold = locateFold(this, ctx, prevStateObj, candState);
          } catch (err) {
            if (err === BUDGET_SPENT && !doStab) {
              const kPrev = points.length - 1;
              this._restore(candState.snap);
              const ptData = readPointData(this, termIdx, pCand, solutions);
              points.push({ drive: pCand, terminals: ptData.terminals, ...(solutions ? { solution: ptData.solution } : {}) });
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

          const exPrev = effDir * (preFold.drive - pPrev);
          const exCand = effDir * (preFold.drive - pCand);
          const floor = Math.max(preFold.floor, candFloor);
          const clearTurnover = preFold.located
            ? (exPrev >= 0 && exCand >= 0 && Math.max(exPrev, exCand) > floor)
            : (exPrev > floor && exCand > floor);
          const stepRetreated = effDir * (pCand - pPrev) < -floor;
          if (!clearTurnover && !stepRetreated && !(exPrev >= -floor && exCand >= -floor)) {
            if (_stats) _stats.roundoffFoldSuppressed++;
            stepFold = false;
            preFold = null;
          } else {
            foldConfirm = stepRetreated ? 'confirm' : 'provisional';
            if (clearTurnover && effDir * dpSpan > 0 && sCur < 1 && preFold.s >= 1) {
              const landed = landExact(prevStateObj, preFold.state, 1, prevTx, prevTV, prevTs);
              if (landed) {
                let stLand = null, hopfLand = null;
                if (doStab) {
                  const chk = checkStepStability(stFullList[points.length - 1], to, landed, null, dsStep);
                  if (!chk.ok) {
                    if (chk.isRefine) consecRefine++;
                    tx.set(prevTx);
                    tV.set(prevTV);
                    ts = prevTs;
                    if (!rejectStep(0.5, `continuation stopped: mode guard (${chk.reason}) near ${formatDrive(to, kind, hasGeometry)}`, true)) break;
                    continue;
                  }
                  stLand = chk.stCand;
                  hopfLand = chk.preHopf;
                }
                acceptPoint(1, landed.tx, landed.tV, landed.tan.ts, stLand, null, hopfLand);
                end = `reached to (${formatDrive(to, kind, hasGeometry)})`;
                break;
              }
              tx.set(prevTx);
              tV.set(prevTV);
              ts = prevTs;
              if (!rejectStep(0.5, `continuation stopped: could not land on target ${formatDrive(to, kind, hasGeometry)}`)) break;
              continue;
            }
            if (clearTurnover && effDir * dpSpan < 0 && sCur > 0 && preFold.s <= 0) {
              const landed = landExact(prevStateObj, preFold.state, 0, prevTx, prevTV, prevTs);
              if (landed) {
                let stLand = null, hopfLand = null;
                if (doStab) {
                  const chk = checkStepStability(stFullList[points.length - 1], p0, landed, null, dsStep);
                  if (!chk.ok) {
                    if (chk.isRefine) consecRefine++;
                    tx.set(prevTx);
                    tV.set(prevTV);
                    ts = prevTs;
                    if (!rejectStep(0.5, `continuation stopped: mode guard (${chk.reason}) near ${formatDrive(p0, kind, hasGeometry)}`, true)) break;
                    continue;
                  }
                  stLand = chk.stCand;
                  hopfLand = chk.preHopf;
                }
                acceptPoint(0, landed.tx, landed.tV, landed.tan.ts, stLand, null, hopfLand);
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
            if (!clearTurnover && !stepRetreated) preFold.quiet = true;
          }
        }

        // One Span Rule — edge s = 1 (`to`).
        if (sCur < 1 && sNext >= 1) {
          const landed = landExact(prevStateObj, candState, 1, prevTx, prevTV, prevTs);
          if (landed) {
            let stLand = null, hopfLand = null;
            if (doStab) {
              const chk = checkStepStability(stFullList[points.length - 1], to, landed, null, dsStep);
              if (!chk.ok) {
                if (chk.isRefine) consecRefine++;
                tx.set(prevTx);
                tV.set(prevTV);
                ts = prevTs;
                if (!rejectStep(0.5, `continuation stopped: mode guard (${chk.reason}) near ${formatDrive(to, kind, hasGeometry)}`, true)) break;
                continue;
              }
              stLand = chk.stCand;
              hopfLand = chk.preHopf;
            }
            consecRefine = 0;
            acceptPoint(1, landed.tx, landed.tV, landed.tan.ts, stLand, null, hopfLand);
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
            let stLand = null, hopfLand = null;
            if (doStab) {
              const chk = checkStepStability(stFullList[points.length - 1], p0, landed, preFold, dsStep);
              if (!chk.ok) {
                if (chk.isRefine) consecRefine++;
                tx.set(prevTx);
                tV.set(prevTV);
                ts = prevTs;
                if (!rejectStep(0.5, `continuation stopped: mode guard (${chk.reason}) near ${formatDrive(p0, kind, hasGeometry)}`, true)) break;
                continue;
              }
              stLand = chk.stCand;
              hopfLand = chk.preHopf;
            }
            consecRefine = 0;
            acceptPoint(0, landed.tx, landed.tV, landed.tan.ts, stLand, preFold, hopfLand, preFold ? 'confirm' : 'none');
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

        // Guard 4: Oscillatory nearer-end rule + Morse-index refinement.
        let stCand = null, preHopf = null;
        if (doStab) {
          const chk = checkStepStability(stFullList[points.length - 1], pCand, candState, preFold && !preFold.quiet ? preFold : null, dsStep);
          if (!chk.ok) {
            if (chk.isRefine) consecRefine++;
            tx.set(prevTx);
            tV.set(prevTV);
            ts = prevTs;
            if (!rejectStep(0.5, `continuation stopped: mode guard (${chk.reason}) near ${formatDrive(pPrev, kind, hasGeometry)}`, true)) break;
            continue;
          }
          stCand = chk.stCand;
          preHopf = chk.preHopf;
        }

        consecRefine = 0;

        // Accept this point!
        ts = tanNext.ts;
        sCur = sNext;
        stepCount++;
        acceptPoint(sNext, tx, tV, ts, stCand, preFold, preHopf, foldConfirm);

        // Adapt arclength stride from Newton iteration count and tangent smoothness.
        if (cor.iterations <= 5 && tanNext.cosTheta > 0.94) {
          dsArc = Math.min(MAX_ARC, dsStep * 1.6);
        } else if (cor.iterations <= 8 && tanNext.cosTheta > 0.88) {
          dsArc = Math.min(MAX_ARC, dsStep * 1.3);
        } else if (cor.iterations >= 14 || tanNext.cosTheta < 0.82) {
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
      emitPendingAtEnd();
      runAudit(true);

      return { terminal, kind, points, events, end, warnings };
    } catch (err) {
      if (err === BUDGET_SPENT) {
        emitPendingAtEnd();
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
