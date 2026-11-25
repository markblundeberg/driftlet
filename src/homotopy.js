// Unified pseudo-arclength homotopy continuation for steady solves.
//
// CONTRACT
//
// Functions over a Solver `s` (src/solver.js; its layout is in docs/internals.md), called by
// Solver._steadyFromHere (src/steady.js) on a direct steady device (`s._directSteady()`) unless
// `opts.homotopy` is explicitly false.
//
//   homotopyContinue(s, opts, r, start0) → { converged, steps, iterations, history, ramped?, hitModelWall? }
//     Traces the 1D steady solution manifold F(x, V_fl, s) = 0 from a known or equilibrium start
//     (s = 0) to the target drives (s = 1) without mutating `terms[k].drive` or calling redrive(),
//     so every conserved amount (lawRefs), floating charge (chargeRefs) and immobile sum
//     (localKept) captured at the start of solveSteady() holds along the entire path.
//     - Warm start (s.solvedV != null && !s.stepped): traces from solvedV and the starting state's
//       currents to stay on the same branch; falls back to cold if warm fails.
//     - Cold start: level contacts, capacitive ports at zeroCharge, floating terminals blended via
//       `s.circuitHomotopy[k] = { alpha, beta, v0, gamma }` (r_k = α I_k + β(V_k − v_0) − γ),
//       and generation ramped via `s.generationScale = exp(30·ln(10)·(s − 1))`.
//     - Honours `opts.budget` (`check()` at every arclength step and round, `record()` on every
//       completed corrector round; never swallows non-Error budget interrupts).
//     - Always restores `s.sourceOverride`, `s.generationScale`, `s.allTerminals`, `s.homotopyGen`,
//       `s.homotopyGenTerm` and `s.circuitHomotopy` in `finally`.
//     - Sets `s.unreached` when a single current-driven terminal is bounded in both voltage
//       directions (limiting current, peak/fold such as Frumkin's dip, or nearing at ±20 V), and
//       returns `hitModelWall: true` when repeated model refusals (e.g. Debye–Hückel non-convexity)
//       terminate the path.
//
//   computeHomotopyDs(s, plan, sVal, Bs, pinDs, circDs)
//     Exact ∂F/∂s after `s.assemble(Infinity)` and `equilibrate(s)`: scaled block residual
//     derivative into `Bs`, pin derivatives into `pinDs`, and floating circuit derivatives into
//     `circDs`. Tested against central finite differences in test/homotopy.test.js.
//
//   solveBorderedArclength(s, rhs, pinRhs, circRhs, gRhs, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws) → ds
//     Solves the (P + K + 1)-bordered system reusing the factorised block-tridiagonal `s.sys`
//     (Woodbury), writing δx into `delta`, δV into `deltaV` and returning δs. Tested against the
//     full unfactored augmented system in test/homotopy.test.js.

import { sourceAt } from './device.js';
import { SolverError } from './errors.js';
import { chargeOnly } from './laws.js';
import { addToState, clipLevels, equilibrate, maxPotentialStep } from './newton.js';

const GEN_LOG_SPAN = 30 * Math.LN10; // generationScale(s) = exp(GEN_LOG_SPAN * (s − 1)): 1e-30 at s = 0, 1 at s = 1
const WS = 10; // arclength weight on s (Δs = 0.1 ↔ ~1 thermal unit RMS motion)

// Indices of dimensionless potential/level unknowns (matching maxPotentialStep, excluding fluxes
// and displacements) for a grid-independent arclength metric.
function potentialIndices(s) {
  const { n, M, rix: R } = s, idx = [];
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

// Apply the homotopy parameter s to sourceOverride, circuitHomotopy and generationScale.
function applyS(s, plan, sVal) {
  const { p0, p1, ch, rampGen } = plan;
  for (let k = 0; k < s.terms.length; k++) {
    if (p0[k] !== p1[k]) s.sourceOverride.set(k, p0[k] + sVal * (p1[k] - p0[k]));
    else s.sourceOverride.delete(k);
  }
  if (s.circuitHomotopy) {
    for (let k = 0; k < s.terms.length; k++) {
      const c = ch[k];
      s.circuitHomotopy[k] = c ? {
        alpha: c.a0 + sVal * c.da,
        beta: c.b0 + sVal * c.db,
        v0: c.v0 + sVal * c.dv,
        gamma: c.g0 + sVal * c.dg,
      } : null;
    }
  }
  if (rampGen) s.generationScale = Math.exp(GEN_LOG_SPAN * (sVal - 1));
}

function clearOverrides(s, plan) {
  for (let k = 0; k < s.terms.length; k++) s.sourceOverride.delete(k);
  if (s.circuitHomotopy) s.circuitHomotopy.fill(null);
  if (plan.rampGen) s.generationScale = 1;
  s._refreshSources();
}

// Compute exact ∂F/∂s after assemble(Infinity) and equilibrate(s):
//   Bs[j]:     ∂(scaled block residual j)/∂s
//   pinDs[q]:  ∂(pin q residual)/∂s
//   circDs[a]: ∂(floating terminal a's circuit residual)/∂s
export function computeHomotopyDs(s, plan, sVal, Bs, pinDs, circDs) {
  const N = s.sys.size, fl = s.floating, pins = s.pins;
  const { p0, p1, ch, rampGen } = plan;
  Bs.fill(0);
  pinDs.fill(0);
  circDs.fill(0);
  // Held voltage terminals: termB[k] is already scaled by rowDiv after equilibrate(s).
  // termDX[fl[a]][k] already includes the α factor from _circuit(fl[a]).
  for (let k = 0; k < s.terms.length; k++) {
    const dp = p1[k] - p0[k];
    const d = s.terms[k].drive;
    if (d.kind === 'V' && !(d.R > 0)) {
      if (dp === 0) continue;
      const Bk = s.termB[k];
      for (let j = 0; j < N; j++) if (Bk[j] !== 0) Bs[j] += dp * Bk[j];
      for (let q = 0; q < pins.length; q++) if (pins[q].dV) pinDs[q] += dp * pins[q].dV[k];
      for (let a = 0; a < fl.length; a++) circDs[a] += dp * s.termDX[fl[a]][k];
    } else {
      const a = fl.indexOf(k);
      if (a < 0) continue;
      const c = ch[k];
      if (c) {
        const beta = c.b0 + sVal * c.db, v = c.v0 + sVal * c.dv;
        circDs[a] += c.da * s.termI[k] + c.db * (s.termV[k] - v) - beta * c.dv - c.dg;
      } else if (dp !== 0) {
        if (d.kind === 'I') circDs[a] -= dp;
        else circDs[a] -= dp / d.R;
      }
    }
  }
  // Generation reactions: rate is linear in generationScale = exp(GEN_LOG_SPAN * (s − 1)).
  if (rampGen && s.homotopyGen) {
    const hGen = s.homotopyGen, rowDiv = s.rowDiv;
    for (let j = 0; j < N; j++) if (hGen[j] !== 0) Bs[j] += (GEN_LOG_SPAN * hGen[j]) / rowDiv[j];
    for (let a = 0; a < fl.length; a++) {
      const fk = fl[a], c = ch[fk];
      const alpha = c ? c.a0 + sVal * c.da : 1;
      circDs[a] += alpha * GEN_LOG_SPAN * s.homotopyGenTerm[fk];
    }
  }
}

// Solve the (P + K + 1)-bordered linear system reusing the already-factorised block-tridiagonal
// Jacobian s.sys (Woodbury):
//   T δx − Σ e_q μ_q + Σ B_k δV_k + B_s δs = rhs
//   W_q·δx + Σ D_pin δV + D_pin,s δs       = pinRhs[q]
//   C_k·δx + Σ D_fl δV  + D_fl,s δs        = circRhs[a]
//   wx·δx  + wV·δV      + ws δs            = gRhs
// Returns δs; writes δx into `delta` and δV into `deltaV`.
export function solveBorderedArclength(s, rhs, pinRhs, circRhs, gRhs, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws) {
  const N = s.sys.size, fl = s.floating, pins = s.pins, P = pins.length, K = fl.length, S = P + K + 1;
  s._solveLinear(rhs, delta);
  const e = s.dWork ?? (s.dWork = new Float64Array(N + 1));
  const cols = s.borderCols ?? (s.borderCols = []);
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
  for (let col = 0; col < S; col++) {
    let p = col;
    for (let r = col + 1; r < S; r++) if (Math.abs(A[r][col]) > Math.abs(A[p][col])) p = r;
    [A[col], A[p]] = [A[p], A[col]];
    if (A[col][col] === 0) throw new SolverError('the bordered arclength system is singular');
    for (let r = col + 1; r < S; r++) {
      const f = A[r][col] / A[col][col];
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
    const sa = a < P ? x[a] : -x[a], col = cols[a];
    for (let j = 0; j < N; j++) delta[j] += sa * col[j];
    if (a >= P && a < P + K) deltaV[a - P] = x[a];
  }
  return x[P + K];
}

// Weighted arclength norm and inner product weights for tangent (tx, tV, ts).
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

function pickVref(s, p1) {
  const nt = s.terms.length;
  for (let k = 0; k < Math.min(2, nt); k++) {
    const d = s.terms[k].drive, ct = s.model.contacts[k === 0 ? 'left' : 'right'];
    if (d.kind === 'V' && !(d.R > 0) && ct.phi.type !== 'capacitive') return p1[k];
  }
  for (let k = 0; k < Math.min(2, nt); k++) {
    const d = s.terms[k].drive;
    if (d.kind === 'V' && !(d.R > 0)) return p1[k];
  }
  for (let k = 2; k < nt; k++) {
    const d = s.terms[k].drive;
    if (d.kind === 'V' && !(d.R > 0) && !s.model.ports[k - 2].capacitance) return p1[k];
  }
  for (let k = 0; k < nt; k++) {
    if (s.terms[k].drive.kind === 'V') return p1[k];
  }
  return 0;
}

// Build a warm or cold homotopy plan:
//   mode = 'warm'       from solvedV / current terminal state (same branch)
//   mode = 'cold'       standard cold start (R-backed starts at R_eff = 0, V = V_target; I with
//                       generation or I = 0 uses shunt-to-current; dark nonzero I ramps 0 → I)
//   mode = 'cold-level' all contacts level at Vref (R-backed at Vref, I uses shunt-to-current)
function buildPlan(s, mode) {
  const nt = s.terms.length;
  const p1 = new Float64Array(nt), p0 = new Float64Array(nt);
  const ch = new Array(nt).fill(null);
  for (let k = 0; k < nt; k++) p1[k] = sourceAt(s.terms[k].drive.src, s.time);

  if (mode === 'warm') {
    if (s.solvedV == null) return null;
    for (let k = 0; k < nt; k++) {
      const d = s.terms[k].drive;
      if (d.kind === 'V' && !(d.R > 0)) p0[k] = k < 2 ? s.solvedV[k] : s.termV[k];
      else p0[k] = p1[k];
    }
    for (let k = 0; k < nt; k++) {
      const d = s.terms[k].drive;
      if (d.kind === 'V' && !(d.R > 0)) s.sourceOverride.set(k, p0[k]);
    }
    try {
      s.computeConcentrations();
      s._assembleBookkeeping(Infinity);
      for (let k = 0; k < nt; k++) {
        const d = s.terms[k].drive;
        if (d.kind === 'I') p0[k] = s.termI[k];
        else if (d.kind === 'V' && d.R > 0) p0[k] = s.termV[k] + d.R * s.termI[k];
      }
    } catch (err) {
      if (!(err instanceof Error)) throw err;
      return null;
    } finally {
      for (let k = 0; k < nt; k++) s.sourceOverride.delete(k);
      s._refreshSources();
    }
    const moves = p0.some((v, k) => v !== p1[k]);
    if (!moves) return null;
    return { p0, p1, ch, rampGen: false, warm: true };
  }

  const Vref = pickVref(s, p1);
  const rampGen = s.hasGeneration;
  const levelR = mode === 'cold-level';
  const shuntI = mode === 'cold-level';
  const hasAnchor = s.terms.some((t) => t.drive.kind === 'V' && !(t.drive.R > 0));
  // If every other terminal is driven by a current and one terminal is behind R, Kirchhoff gives
  // its steady voltage V_src − I·R (so starting R_eff = 0 there moves nothing when R turns on).
  const rKs = s.floating.filter((k) => s.terms[k].drive.kind === 'V' && s.terms[k].drive.R > 0);
  const kirchhoffR = rKs.length === 1 && s.terms.every((t, j) => j === rKs[0] || t.drive.kind === 'I');

  for (let k = 0; k < nt; k++) {
    const d = s.terms[k].drive;
    if (d.kind === 'I') {
      if (s._bound(k)) {
        p0[k] = p1[k];
      } else if (hasAnchor && (shuntI || rampGen)) {
        // Hold V_k = Vref at s = 0 via virtual shunt (1 − s) G_0 (V_k − Vref) + s (I_k − I_target) = 0.
        const G0 = Math.max(1, Math.abs(p1[k]) / s.VT);
        ch[k] = { a0: 0, da: 1, b0: G0, db: -G0, v0: Vref, dv: 0, g0: 0, dg: p1[k] };
        p0[k] = p1[k];
      } else {
        p0[k] = 0;
      }
    } else if (d.R > 0) {
      if (s._bound(k)) {
        p0[k] = p1[k];
      } else if (chargeOnly(s, k)) {
        p0[k] = Vref;
      } else {
        // R_eff(s) = s·R (with O(1) initial conductance when light also ramps on a high-R device):
        //   r_k = s·I_k + ((1 − s)·b0 + s/R)·(V_k − v(s)) = 0
        let vStart = levelR ? Vref : p1[k];
        if (!levelR && kirchhoffR && k === rKs[0]) {
          const Itot = -s.terms.reduce((t, tm, j) => (j === k ? t : t + p1[j]), 0);
          vStart = p1[k] - Itot * d.R;
        }
        const b1 = 1 / d.R, b0 = rampGen ? Math.max(1, b1) : b1;
        ch[k] = { a0: 0, da: 1, b0, db: b1 - b0, v0: vStart, dv: p1[k] - vStart, g0: 0, dg: 0 };
        p0[k] = p1[k];
      }
    } else if (k >= 2 && s.model.ports[k - 2].capacitance) {
      p0[k] = s.model.ports[k - 2].capacitance.zeroCharge;
    } else {
      p0[k] = Vref;
    }
  }
  const moves = rampGen || p0.some((v, k) => v !== p1[k]) || ch.some((c) => c !== null);
  if (!moves) return null;
  const fastS0 = (!levelR && rKs.some((k) => ch[k] && ch[k].v0 !== Vref)) || (levelR && !rampGen && rKs.length === 0);
  return { p0, p1, ch, rampGen, warm: false, fastS0 };
}

// Lay out a cold start at s = 0, treating any floating terminal currently voltage-held by
// circuitHomotopy (alpha === 0) as held at its v0 so initFromComposition sets electrodes there.
function layOutAtS0(s, plan, budget) {
  const savedFl = s.floating;
  const heldAt0 = savedFl.filter((k) => plan.ch[k]?.a0 === 0);
  if (heldAt0.length === 0) return s._layOut(budget);
  s.floating = savedFl.filter((k) => plan.ch[k]?.a0 !== 0);
  for (const k of heldAt0) s.termV[k] = plan.ch[k].v0;
  try {
    return s._layOut(budget);
  } finally {
    s.floating = savedFl;
    for (const k of heldAt0) s.termV[k] = plan.ch[k].v0;
  }
}

// Compute unit tangent (tx, tV, ts) at the current converged point (u, uLo, termV, sVal),
// oriented to agree with (prevTx, prevTV, prevTs) (or with ts > 0 / tV[0] > 0 if prevTx is null).
// When `reuseFactor` is true, reuses the already-factored s.sys and (Bs, pinDs, circDs) from a
// just-converged correctorNewton step (within 1e-6 of the solution).
function computeTangent(s, plan, sVal, potIdx, Bs, pinDs, circDs, wx, wV, prevTx, prevTV, prevTs, tx, tV, orientV, reuseFactor = false) {
  const N = s.sys.size, K = s.floating.length;
  if (!reuseFactor) {
    applyS(s, plan, sVal);
    s.islandsOn = false;
    s.assemble(Infinity);
    equilibrate(s);
    computeHomotopyDs(s, plan, sVal, Bs, pinDs, circDs);
    s._factor();
  }
  const zero = new Float64Array(N + 1);
  const zeroPin = new Float64Array(s.pins.length);
  const zeroCirc = new Float64Array(K);
  let ws;
  if (prevTx) {
    ws = buildWeights(potIdx, s.VT, prevTx, prevTV, prevTs, wx, wV);
  } else if (orientV) {
    wx.fill(0);
    wV.fill(0);
    wV[0] = 1;
    ws = 0;
  } else {
    wx.fill(0);
    wV.fill(0);
    ws = 1;
  }
  const tsRaw = solveBorderedArclength(s, zero, zeroPin, zeroCirc, 1, tx, tV, Bs, pinDs, circDs, wx, wV, ws);
  if (!Number.isFinite(tsRaw)) return NaN;
  const nrm = tangentNorm(potIdx, s.VT, tx, tV, tsRaw);
  if (!(nrm > 0 && Number.isFinite(nrm))) return NaN;
  for (let j = 0; j < N; j++) tx[j] /= nrm;
  for (let a = 0; a < K; a++) tV[a] /= nrm;
  let ts = tsRaw / nrm;
  if (!prevTx) {
    const flip = orientV ? tV[0] < 0 : ts < 0;
    if (flip) {
      for (let j = 0; j < N; j++) tx[j] = -tx[j];
      for (let a = 0; a < K; a++) tV[a] = -tV[a];
      ts = -ts;
    }
  }
  return ts;
}

// Corrector Newton solve on F(x, V_fl, s) = 0 and wx·(x − x_pred) + wV·(V − V_pred) + ws·(s − s_pred) = 0.
function correctorNewton(s, plan, sStart, potIdx, uPred, uLoPred, vPred, sPred, wx, wV, ws, Bs, pinDs, circDs, { maxIter = 60, tol = 1e-6, maxStep = 10, budget } = {}) {
  budget?.check();
  const done = (res) => {
    budget?.record(res);
    return res;
  };
  const N = s.sys.size, fl = s.floating, K = fl.length;
  const delta = s.delta, deltaV = s.deltaV ?? (s.deltaV = new Float64Array(K));
  const pinRhs = new Float64Array(Math.max(1, s.pins.length));
  const circRhs = new Float64Array(K);
  let sVal = sStart;
  for (let it = 1; it <= maxIter; it++) {
    if (plan.rampGen && (sVal < -0.5 || sVal > 1.5)) return done({ converged: false, iterations: it, s: sVal, error: 's out of range' });
    applyS(s, plan, sVal);
    s.islandsOn = false;
    try {
      s.assemble(Infinity);
    } catch (err) {
      if (!(err instanceof SolverError)) throw err;
      s.modelError = err.message;
      return done({ converged: false, iterations: it, s: sVal, modelError: err.message });
    }
    equilibrate(s);
    computeHomotopyDs(s, plan, sVal, Bs, pinDs, circDs);
    try {
      s._factor();
    } catch (err) {
      if (!(err instanceof Error)) throw err;
      return done({ converged: false, iterations: it, s: sVal, error: err.message });
    }
    let g = ws * (sVal - sPred);
    for (let i = 0; i < potIdx.length; i++) {
      const j = potIdx[i], o = s.fullOf[j];
      g += wx[j] * ((s.u[o] - uPred[o]) + (s.uLo[o] - uLoPred[o]));
    }
    for (let a = 0; a < K; a++) g += wV[a] * (s.termV[fl[a]] - vPred[a]);
    const pins = s.pins;
    for (let q = 0; q < pins.length; q++) pinRhs[q] = pins[q].res;
    for (let a = 0; a < K; a++) circRhs[a] = s.termRes[fl[a]];
    let ds;
    try {
      deltaV.fill(0);
      ds = solveBorderedArclength(s, s.res, pinRhs, circRhs, g, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws);
    } catch (err) {
      if (!(err instanceof SolverError)) throw err;
      return done({ converged: false, iterations: it, s: sVal, error: err.message });
    }
    clipLevels(s, delta, maxStep);
    const step = Math.max(maxPotentialStep(s, delta, deltaV), WS * Math.abs(ds));
    const own = Math.max(maxPotentialStep(s, delta, deltaV, false), WS * Math.abs(ds));
    if (!Number.isFinite(step)) return done({ converged: false, iterations: it, s: sVal, error: 'non-finite update' });
    if (it > 1 && own > 1e5) return done({ converged: false, iterations: it, s: sVal, error: 'diverging' });
    const alpha = own > maxStep ? maxStep / own : 1;
    addToState(s, delta, -alpha);
    for (let a = 0; a < K; a++) s.termV[fl[a]] -= alpha * deltaV[a];
    sVal -= alpha * ds;
    if (alpha === 1 && step < tol) {
      try {
        s.computeConcentrations();
      } catch (err) {
        if (!(err instanceof SolverError)) throw err;
        s.modelError = err.message;
        return done({ converged: false, iterations: it, s: sVal, modelError: err.message });
      }
      return done({ converged: true, iterations: it, s: sVal, reusable: s.sys.perturbedPivots === 0 && s.sys.lostColumns === 0 });
    }
  }
  return done({ converged: false, iterations: maxIter, s: sVal });
}

// Run one homotopy path from s = 0 to s = 1 for `plan`.
function tracePath(s, plan, opts, r, wasAll) {
  const N = s.sys.size, fl = s.floating, K = fl.length;
  const potIdx = potentialIndices(s);
  const Bs = new Float64Array(N + 1);
  const pinDs = new Float64Array(Math.max(1, s.constraints.length + s.affinities.length + s.islands.length + 4));
  const circDs = new Float64Array(K);
  const wx = new Float64Array(N + 1), wV = new Float64Array(K);
  const tx = new Float64Array(N + 1), tV = new Float64Array(K);
  const prevTx = new Float64Array(N + 1), prevTV = new Float64Array(K);
  const uPred = new Float64Array(s.u.length), uLoPred = new Float64Array(s.uLo.length), vPred = new Float64Array(K);
  const sub = { ...opts, maxSteps: 1 };

  let steps = r.steps ?? 0, iterations = r.iterations ?? 0;
  const history = (r.history ?? []).slice();

  // 1. Solve at s = 0. When a speculative cold start holds R-backed terminals at their target
  // voltage (not Vref), try only a direct solve before falling back to cold-level so a wall at
  // the target doesn't spend pseudo-transient steps.
  applyS(s, plan, 0);
  for (const k of fl) if (plan.ch[k]?.a0 === 0) s.termV[k] = plan.ch[k].v0;
  const q0Opts = plan.warm || plan.fastS0 ? sub : { ...opts, maxSteps: Math.min(opts.maxSteps ?? 80, 20) };
  let q0 = s._steadySteps(q0Opts);
  steps += q0.steps;
  iterations += q0.iterations;
  if (!q0.converged && !plan.warm && s.unkept.length === 0) {
    if (layOutAtS0(s, plan, opts.budget)) {
      applyS(s, plan, 0);
      q0 = s._steadySteps(q0Opts);
      steps += q0.steps;
      iterations += q0.iterations;
    }
  }
  if (!q0.converged) return { converged: false, steps, iterations, history, s0Failed: true };

  // Pure single current-driven terminal (no virtual shunt, no generation ramp): march both V
  // directions from s = 0 so both limiting currents and sign slips are detected cleanly.
  const singleI = fl.length === 1 && s.terms[fl[0]].drive.kind === 'I' && !s._bound(fl[0]) && !plan.ch[fl[0]] && plan.p0.every((v, k) => k === fl[0] || v === plan.p1[k]) && !plan.rampGen ? fl[0] : -1;
  const seen = [];
  if (singleI >= 0) seen.push([s.termV[singleI], s.termI[singleI]]);

  // 2. Initial unit tangent at s = 0.
  const orientV = singleI >= 0;
  let ts;
  opts.budget?.check();
  try {
    ts = computeTangent(s, plan, 0, potIdx, Bs, pinDs, circDs, wx, wV, null, null, null, tx, tV, orientV);
  } catch (err) {
    if (!(err instanceof Error)) throw err;
    return { converged: false, steps, iterations, history };
  }
  if (!Number.isFinite(ts)) return { converged: false, steps, iterations, history };

  const start0 = s._levels();
  const vStart = Float64Array.from(fl, (k) => s.termV[k]);
  let limited = 0, nearing = null, hitModelWall = false;

  // Try landing on the true target (s = 1). When `checkMonotone` is true (a single current-driven
  // terminal), verify the tangent at s = 1 still has ts > 0 along the arriving direction so an
  // N-shaped I–V curve can never leap over a peak onto the unstable descending branch. On devices
  // with floating islands, reject a coarse jump that misses Kirchhoff so homotopy steps closer.
  const tryLand = (predScale, checkMonotone = false, strictKirchhoff = false) => {
    const before = s._levels();
    if (predScale !== 0) {
      const maxPt = maxPotentialStep(s, tx, tV) * Math.abs(predScale);
      const a = maxPt > 15 ? (15 / maxPt) * predScale : predScale;
      addToState(s, tx, a);
      for (let a0 = 0; a0 < K; a0++) s.termV[fl[a0]] += a * tV[a0];
    }
    let q;
    if (plan.keepChOnLand) {
      applyS(s, plan, 1);
      q = s._steadySteps(sub);
    } else {
      clearOverrides(s, plan);
      // Run the landing steady round in the device's own terminal mode so its diagnostics
      // (islandsFailed, islandsWorse, islandImbalance, conditioning) never come from all-terminals mode.
      s.allTerminals = wasAll;
      try {
        q = s._steadySteps(sub);
      } finally {
        s.allTerminals = true;
      }
    }
    steps += q.steps;
    iterations += q.iterations;
    history.push({ homotopy: 1, converged: q.converged });
    if (q.converged) {
      if (strictKirchhoff && s.islands.length > 0 && s._missesKirchhoff()) {
        s._putLevels(before);
        return null;
      }
      if (checkMonotone && singleI >= 0) {
        let tsLand = NaN;
        try {
          tsLand = computeTangent(s, plan, 1, potIdx, Bs, pinDs, circDs, wx, wV, tx, tV, ts, prevTx, prevTV, false);
        } catch (err) {
          if (!(err instanceof Error)) throw err;
          tsLand = NaN;
        }
        if (plan.keepChOnLand) applyS(s, plan, 1);
        else clearOverrides(s, plan);
        if (Number.isFinite(tsLand) && tsLand <= 0) {
          s._putLevels(before);
          return null;
        }
      }
      return q;
    }
    s._putLevels(before);
    return null;
  };
  if (singleI < 0 && s.islands.length === 0 && ts > 0.0995) {
    const landed = tryLand(0);
    if (landed) return { converged: true, steps, iterations, history, ramped: true };
  }

  const dirs = singleI >= 0 ? (ts >= 0 ? [1, -1] : [-1, 1]) : [1];
  for (let di = 0; di < dirs.length; di++) {
    const dir = dirs[di];
    if (di > 0) {
      opts.budget?.check();
      s._putLevels(start0);
      try {
        ts = computeTangent(s, plan, 0, potIdx, Bs, pinDs, circDs, wx, wV, null, null, null, tx, tV, orientV);
      } catch (err) {
        if (!(err instanceof Error)) throw err;
        break;
      }
    }
    if (dir === -1) {
      for (let j = 0; j < N; j++) tx[j] = -tx[j];
      for (let a = 0; a < K; a++) tV[a] = -tV[a];
      ts = -ts;
    }

    let sCur = 0;
    let dsArc = Math.abs(ts) > 1e-5 ? 0.15 / Math.abs(ts) : 4;
    if (singleI >= 0 && Math.abs(tV[0]) > 1e-9) dsArc = Math.min(dsArc, 0.1 / Math.abs(tV[0]));
    if (plan.rampGen && Math.abs(ts) > 1e-6) dsArc = Math.min(dsArc, 0.08 / Math.abs(ts));
    const minArc = 1e-4, maxArc = 1e6, maxHomSteps = 200;
    let foldSeen = false, sMax = 0, prevAbsDist = Math.abs(1 - sCur), modelFails = 0, consecFails = 0, recentFail = 0;

    for (let hs = 0; hs < maxHomSteps; hs++) {
      opts.budget?.check();
      // Cap Δs when generation ramps logarithmically (~3.5 decades/step max) to avoid exp overflow.
      if (plan.rampGen && Math.abs(ts) > 1e-6) dsArc = Math.min(dsArc, 0.12 / Math.abs(ts));
      // On a single current-driven terminal, bound ΔV per step so an N-shaped I–V curve resolves
      // any peak rather than leaping across it onto an unstable descending branch.
      if (singleI >= 0 && Math.abs(tV[0]) > 1e-9) {
        const dVsoFar = Math.abs(s.termV[singleI] - vStart[0]);
        const maxDV = dVsoFar >= 2 && Math.abs(sCur) < 0.01 ? 2.0 : 0.35;
        dsArc = Math.min(dsArc, maxDV / Math.abs(tV[0]));
      }

      if (sCur < 1 && ts > 1e-9 && sCur + dsArc * ts >= 1 && !foldSeen) {
        const predScale = (1 - sCur) / ts;
        if (singleI < 0 || Math.abs(predScale * tV[0]) <= 0.35) {
          const landed = tryLand(predScale, singleI >= 0, sCur < 0.99);
          if (landed) return { converged: true, steps, iterations, history, ramped: true };
          dsArc = Math.min(dsArc * 0.5, (1 - sCur) / (2 * ts));
          if (dsArc < minArc) break;
        }
      }

      const before = s._levels();
      const sPred = sCur + dsArc * ts;
      // Target hyperplane passes through z_pred = z_cur + dsArc·t, while the Newton starting guess
      // clamps pointwise potential extrapolation to ≤15 V_T so moving fronts (e.g. a 40 V depletion
      // sweep) don't overshoot into exponential overflow on iteration 1.
      const maxPt = maxPotentialStep(s, tx, tV) * dsArc;
      const predAlpha = maxPt > 10 ? 10 / maxPt : 1;
      addToState(s, tx, dsArc);
      for (let a = 0; a < K; a++) s.termV[fl[a]] += dsArc * tV[a];
      uPred.set(s.u);
      uLoPred.set(s.uLo);
      for (let a = 0; a < K; a++) vPred[a] = s.termV[fl[a]];
      if (predAlpha < 1) {
        s._putLevels(before);
        addToState(s, tx, dsArc * predAlpha);
        for (let a = 0; a < K; a++) s.termV[fl[a]] += dsArc * predAlpha * tV[a];
      }
      const sInit = sCur + dsArc * predAlpha * ts;
      const ws = buildWeights(potIdx, s.VT, tx, tV, ts, wx, wV);

      let cor = correctorNewton(s, plan, sInit, potIdx, uPred, uLoPred, vPred, sPred, wx, wV, ws, Bs, pinDs, circDs, { budget: opts.budget });
      steps++;
      iterations += cor.iterations;
      history.push({ homotopy: cor.s, converged: cor.converged });

      // On devices with weakly held islands or severe cancellation (where plain factorisation loses
      // digits that newton()'s island rounds and GMRES refinement recover), try a fixed-s steady
      // step at sPred (or s = 1) before shrinking.
      const illCond = s.islands.length > 0 || s.sys.perturbedPivots > 0 || s.sys.lostColumns > 0 || (s.sys.minPivotRatio > 0 && s.sys.minPivotRatio < 1e-10);
      if (!cor.converged && !cor.modelError && singleI < 0 && illCond && sPred > sCur && sPred <= 1) {
        s._putLevels(before);
        const sTry = Math.min(1, sPred);
        applyS(s, plan, sTry);
        const qIsl = s._steadySteps(sub);
        steps += qIsl.steps;
        iterations += qIsl.iterations;
        history.push({ homotopy: sTry, converged: qIsl.converged });
        if (qIsl.converged) {
          if (sTry >= 1) {
            const landed = tryLand(0);
            if (landed) return { converged: true, steps, iterations, history, ramped: true };
          }
          cor = { converged: true, iterations: qIsl.iterations, s: sTry };
        }
      }

      if (!cor.converged) {
        s._putLevels(before);
        if (cor.modelError) modelFails++;
        else modelFails = 0;
        consecFails++;
        recentFail = 3;
        dsArc *= consecFails > 1 ? 0.25 : 0.5;
        if (dsArc < minArc || (modelFails >= 5 && dsArc * Math.abs(ts) < 1e-3)) {
          if (cor.modelError && modelFails >= 4) hitModelWall = true;
          if (singleI >= 0 && (Math.abs(1 - sCur) > prevAbsDist || (Math.abs(s.termV[singleI] - vStart[0]) >= 1 && sCur < 1e-3))) limited++;
          break;
        }
        continue;
      }
      modelFails = 0;
      consecFails = 0;

      const sNext = cor.s;
      if (singleI >= 0) seen.push([s.termV[singleI], s.termI[singleI]]);

      if ((sCur < 1 && sNext >= 1 && !(singleI >= 0 && foldSeen)) || (singleI < 0 && sCur > 1 && sNext <= 1)) {
        const frac = (1 - sNext) / (sCur - sNext);
        if (frac > 0 && frac < 1) {
          for (let j = 0; j < s.u.length; j++) {
            s.u[j] = s.u[j] + frac * (before.snap.u[j] - s.u[j]);
            s.uLo[j] = s.uLo[j] + frac * (before.snap.uLo[j] - s.uLo[j]);
          }
          for (let a = 0; a < K; a++) s.termV[fl[a]] += frac * (before.snap.termV[fl[a]] - s.termV[fl[a]]);
        }
        const landed = tryLand(0, singleI >= 0, sCur < 0.99);
        if (landed) return { converged: true, steps, iterations, history, ramped: true };
        s._putLevels(before);
        dsArc *= 0.25;
        if (dsArc < minArc) break;
        continue;
      }

      prevTx.set(tx);
      prevTV.set(tV);
      const prevTs = ts;
      let nextTs = NaN;
      try {
        nextTs = computeTangent(s, plan, sNext, potIdx, Bs, pinDs, circDs, wx, wV, prevTx, prevTV, prevTs, tx, tV, false, !!cor.reusable);
      } catch (err) {
        if (!(err instanceof Error)) throw err;
        nextTs = NaN;
      }
      if (!Number.isFinite(nextTs) || (singleI < 0 && illCond && sNext > sCur && sNext < 1 && s.sys.minPivotRatio > 0 && s.sys.minPivotRatio < 1e-13)) {
        // On a device with weakly held islands or severe cancellation, tangent factorisation
        // without GMRES/islands can lose all digits even though newton() converges: advance s
        // monotonically with fixed-s steady steps from sNext to 1.
        if (singleI < 0 && illCond && sNext > sCur && sNext < 1) {
          let sM = sNext, dM = Math.min(1 - sNext, Math.max(0.05, sNext - sCur));
          while (sM < 1) {
            opts.budget?.check();
            const nxt = Math.min(1, sM + dM), bef = s._levels();
            applyS(s, plan, nxt);
            const qM = s._steadySteps(sub);
            steps += qM.steps;
            iterations += qM.iterations;
            history.push({ homotopy: nxt, converged: qM.converged });
            if (qM.converged) {
              sM = nxt;
              dM = Math.min(1 - sM, dM * 1.5);
              if (sM >= 1) {
                const landed = tryLand(0);
                if (landed) return { converged: true, steps, iterations, history, ramped: true };
                break;
              }
            } else {
              s._putLevels(bef);
              dM = (nxt - sM) * 0.25;
              if (dM < 1e-5) break;
            }
          }
          break;
        }
        s._putLevels(before);
        dsArc *= 0.25;
        if (dsArc < minArc) break;
        continue;
      }
      ts = nextTs;

      if (sNext > sMax) sMax = sNext;
      const distNow = Math.abs(1 - sNext);

      if (prevTs > 0 && ts <= 0 && sMax > 0.05) foldSeen = true;
      if (foldSeen && (sNext < -0.1 || (singleI >= 0 && sNext <= 0.05 * sMax))) {
        if (singleI >= 0) limited++;
        break;
      }

      if (singleI >= 0) {
        const dVtot = Math.abs(s.termV[singleI] - vStart[0]);
        if (sNext <= -1.5) {
          // Passed >1.5× the target current in the opposite direction: clearly limited in this
          // direction, and `seen` now passes −target for the sign-slip hint.
          limited++;
          break;
        }
        if (dVtot >= 5 && Math.abs(sNext - sCur) < 1e-8 * Math.max(1e-6, Math.abs(1 - sNext))) {
          // Completely saturated (dI/dV ≈ 0) after moving ≥5 V away from start.
          limited++;
          break;
        }
        if (dVtot >= 20) {
          const dVstep = Math.abs(s.termV[singleI] - before.snap.termV[singleI]);
          const dsStep = sNext - sCur;
          if (dsStep > 1e-6 * (1 - sMax + 1e-12) && dVstep > 0) {
            const ahead = ((1 - sNext) * dVstep) / dsStep;
            if (ahead <= 10 * dVtot) nearing = { at: s.termV[singleI], guess: s.termV[singleI] + Math.sign(s.termV[singleI] - vStart[0]) * ahead };
            else limited++;
          } else limited++;
          break;
        }
      }

      prevAbsDist = Math.min(prevAbsDist, distNow);
      sCur = sNext;

      if (recentFail > 0) {
        recentFail--;
        if (cor.iterations <= 14) dsArc = Math.min(maxArc, dsArc * 1.25);
      } else if (cor.iterations <= 10) dsArc = Math.min(maxArc, dsArc * 2.0);
      else if (cor.iterations <= 30) dsArc = Math.min(maxArc, dsArc * 1.5);
      else if (cor.iterations >= 45) dsArc = Math.max(minArc, dsArc * 0.7);
    }
    if (singleI >= 0 && Math.abs(s.termV[singleI] - vStart[0]) >= 5 && Math.abs(sCur) < 0.01 && limited === 0 && !nearing) limited++;
  }

  if (singleI >= 0 && seen.length > 1 && limited + (nearing === null ? 0 : 1) === 2) {
    const Vs = seen.map(([V]) => V), Is = seen.map(([, I]) => I);
    s.unreached = { terminal: s.terms[singleI].name, target: plan.p1[singleI], V: [Math.min(...Vs), Math.max(...Vs)], I: [Math.min(...Is), Math.max(...Is)], nearing };
  }
  return { converged: false, steps, iterations, history, hitModelWall };
}

// Entry point for unified steady-state homotopy continuation.
export function homotopyContinue(s, opts, r, start0) {
  if (!s._directSteady()) return { converged: false, steps: r.steps ?? 0, iterations: r.iterations ?? 0, history: r.history ?? [] };
  const wasAll = s.allTerminals;
  s.allTerminals = true;
  if (s.hasGeneration) {
    s.homotopyGen = new Float64Array(s.sys.size + 1);
    s.homotopyGenTerm = new Float64Array(s.terms.length);
  }
  s.circuitHomotopy = new Array(s.terms.length).fill(null);
  try {
    let cur = r;
    const runMode = (mode) => {
      opts.budget?.check();
      s._putLevels(start0);
      const plan = buildPlan(s, mode);
      if (!plan) return null;
      try {
        return tracePath(s, plan, opts, cur, wasAll);
      } finally {
        clearOverrides(s, plan);
      }
    };
    if (s.solvedV != null && !s.stepped) {
      const w = runMode('warm');
      if (w) {
        if (w.converged) return w;
        cur = w;
      }
    }
    const c = runMode('cold');
    if (c) {
      if (c.converged || c.hitModelWall || s.unreached) return c;
      cur = c;
    }
    // If cold had nothing to move or failed at s = 0 (e.g. V_src far from equilibrium behind R,
    // or dark I = 0 singular on a fine grid), try the level-contacts / shunt start.
    if (!c || c.s0Failed) {
      const cl = runMode('cold-level');
      if (cl) {
        if (cl.converged || cl.hitModelWall || s.unreached) return cl;
        cur = cl;
      }
    }
    // Decoupled 2-stage fallback for lit devices with floating terminals (e.g. Schottky + pn solar
    // cell whose coupled shunt+light path folds): first ramp generation to 1 with floating
    // terminals held at Vref via circuitHomotopy (no drive mutation), then trace from that lit
    // state to the target current/R at full light.
    if (s.hasGeneration && s.floating.some((k) => !s._bound(k))) {
      opts.budget?.check();
      s._putLevels(start0);
      const nt = s.terms.length;
      const p1 = new Float64Array(nt);
      for (let k = 0; k < nt; k++) p1[k] = sourceAt(s.terms[k].drive.src, s.time);
      const Vref = pickVref(s, p1);
      const p0 = new Float64Array(nt), ch1 = new Array(nt).fill(null);
      for (let k = 0; k < nt; k++) {
        const d = s.terms[k].drive;
        if (s.floating.includes(k) && !s._bound(k) && !chargeOnly(s, k)) {
          ch1[k] = { a0: 0, da: 0, b0: 1, db: 0, v0: Vref, dv: 0, g0: 0, dg: 0 };
          p0[k] = p1[k];
        } else if (d.kind === 'V' && !(d.R > 0)) {
          p0[k] = k >= 2 && s.model.ports[k - 2].capacitance ? s.model.ports[k - 2].capacitance.zeroCharge : Vref;
        } else {
          p0[k] = p1[k];
        }
      }
      // Phase 1: ramp generation (and held voltages) to s = 1 while floating terminals stay held at Vref.
      const plan1 = { p0, p1: Float64Array.from(p0), ch: ch1, rampGen: true, warm: false, keepChOnLand: true };
      let g1;
      try {
        g1 = tracePath(s, plan1, opts, cur, wasAll);
      } finally {
        clearOverrides(s, plan1);
      }
      cur = g1;
      if (g1.converged) {
        // Phase 2: at full light (generationScale = 1), trace from lit Vref state to target I / R / V.
        const savedSolvedV = s.solvedV;
        s.solvedV = [s.termV[0], s.termV[1]];
        const plan2 = buildPlan(s, 'warm');
        s.solvedV = savedSolvedV;
        if (plan2) {
          try {
            const g2 = tracePath(s, plan2, opts, cur, wasAll);
            if (g2.converged || g2.hitModelWall || s.unreached) return g2;
            cur = g2;
          } finally {
            clearOverrides(s, plan2);
          }
        }
      }
    }
    // Held-V 1D root-finder fallback when a single floating terminal's bordered solve loses
    // precision (e.g. a stiff chain behind a large R where termC·T⁻¹ amplifies cancellation):
    // holds V_k via circuitHomotopy (α = 0, β = 1, zero drive mutation) and finds I(V) = want(V)
    // by bracketed Illinois, then polishes floated once at V*.
    if (s.floating.length === 1 && !s._bound(s.floating[0]) && !chargeOnly(s, s.floating[0]) && !s.hasGeneration) {
      opts.budget?.check();
      const k = s.floating[0], d = s.terms[k].drive, target = sourceAt(d.src, s.time);
      const want = d.kind === 'I' ? () => target : (V) => (target - V) / d.R;
      const sub = { ...opts, maxSteps: 1 };
      let steps = cur.steps ?? 0, iterations = cur.iterations ?? 0;
      const history = (cur.history ?? []).slice();
      const evalAt = (V) => {
        opts.budget?.check();
        const bef = s._levels();
        s.circuitHomotopy[k] = { alpha: 0, beta: 1, v0: V, gamma: 0 };
        s.termV[k] = V;
        const q = s._steadySteps(sub);
        steps += q.steps;
        iterations += q.iterations;
        if (q.converged) return s.termI[k] - want(V);
        s._putLevels(bef);
        return NaN;
      };
      s._putLevels(start0);
      const V0 = d.kind === 'V' ? target : s.termV[k];
      let f0 = evalAt(V0);
      if (Number.isFinite(f0)) {
        const step0 = d.kind === 'V' ? Math.max(-0.5, Math.min(0.5, -f0 * d.R)) : 0.02;
        const dirs = step0 !== 0 ? [Math.sign(step0), -Math.sign(step0)] : [1, -1];
        let Va = V0, fa = f0, Vb = NaN, fb = NaN;
        const snap0 = s._levels();
        for (const sg of dirs) {
          s._putLevels(snap0);
          let V = V0, f = f0, dv = Math.max(0.02, Math.min(0.5, Math.abs(step0) || 0.02)) * sg;
          for (let j = 0; j < 25 && Math.abs(V - V0) < 20; j++) {
            const Vn = V + dv, fn = evalAt(Vn);
            if (!Number.isFinite(fn)) { dv *= 0.25; if (Math.abs(dv) < 1e-4) break; continue; }
            if (Math.sign(fn) !== Math.sign(f) || fn === 0) { Va = V; fa = f; Vb = Vn; fb = fn; break; }
            V = Vn; f = fn; dv *= 1.6;
          }
          if (Number.isFinite(fb)) break;
        }
        if (Number.isFinite(fb)) {
          let side = 0;
          const scale = Math.max(d.kind === 'I' ? Math.abs(target) : 0, Math.abs(fa + want(Va)), Math.abs(fb + want(Vb)), 1e-18);
          for (let it = 0; it < 40; it++) {
            const Vm = fb === fa ? Va : Va - (fa * (Vb - Va)) / (fb - fa), fm = evalAt(Vm);
            if (!Number.isFinite(fm)) break;
            if (Math.abs(fm) <= 1e-10 * scale) {
              s.circuitHomotopy[k] = null;
              s.allTerminals = wasAll;
              let pol;
              try {
                pol = s._steadySteps(sub);
              } finally {
                s.allTerminals = true;
              }
              steps += pol.steps;
              iterations += pol.iterations;
              if (pol.converged) return { converged: true, steps, iterations, history, ramped: true };
              return { converged: true, steps, iterations, history, residual: Math.abs(fm), ramped: true };
            }
            if (Math.sign(fm) === Math.sign(fb)) { Vb = Vm; fb = fm; if (side === -1) fa /= 2; side = -1; }
            else { Va = Vm; fa = fm; if (side === 1) fb /= 2; side = 1; }
          }
        }
      }
      s.circuitHomotopy[k] = null;
      cur = { ...cur, steps, iterations, history };
    }
    return cur;
  } finally {
    s.allTerminals = wasAll;
    s.homotopyGen = null;
    s.homotopyGenTerm = null;
    s.circuitHomotopy = null;
  }
}
