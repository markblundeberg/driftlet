// Newton's iteration: the equations at one time (a step's end, or a steady state at dt = Infinity)
// solved from the state as it stands, its linear solves refined by GMRES where the factorised
// Jacobian loses a mode.
//
// CONTRACT
//
// newton(s, dt, { maxIter = 60, tol = 1e-10, maxStep = 10, islands = false }) iterates on a Solver
// `s` (src/solver.js; its layout is in docs/internals.md): its unknowns u (+ uLo, compensated) and
// the voltages termV of its floating terminals. It's called by step() in src/timeSteps.js, which
// sets what the equations are (s.cOld, the sources' time, s.steady, s.combining) and restores the
// state where it fails. Each iteration calls s.assemble(dt) (the residual s.res, the Jacobian
// s.sys, the floating terminals' terms termRes, termB, termC, termDI, termDX, and in a steady
// solve the pins s.pins), then s._factor(), solveBordered (src/lawRows.js) and, where the
// assembly transformed the unknowns (s.transformed), s._untransform.
//
// Each iteration, in order:
//   1. Assembled; every row of the Jacobian (its three blocks), with its residual and its termB
//      entries, divided by its largest entry (s.rowDiv; a row of zeros is left). A SolverError
//      from the assembly or the solve, or any error from the factorisation, fails the call with
//      its message.
//   2. Solved for the update; refined (7) where refining is on.
//   3. Its size: `step`, the largest |Δ| over φ̂ and each present η at the grid's nodes (a metal's
//      carrier, not its flux slot), each coverage and gate fraction, each diffuse layer's drop and
//      its drawn species' plane levels, and each floating terminal's |ΔV|/V_T; `own`, the same
//      without the terminals and the diffuse layers, each plane level's own update clipped to
//      ±maxStep first.
//   4. Fails where step isn't finite ('non-finite update'), or where own > 1e4 after the first
//      iteration, or own > 10× the first's after the sixth ('diverging').
//   5. Added, damped by α = min(1, maxStep/own) (the terminal voltages by the same α), u and uLo
//      in compensated arithmetic.
//   6. Converged, undamped, where step < tol; or in a time step (it ≥ 2; not on a lost-mode
//      device) where the contraction θ = step/(the last step) < 0.1, step < 1e-6 and
//      θ/(1 − θ)·step < tol; or at round-off (it ≥ 4; refined at least once where refining is on)
//      where step and the last are both under 1e-6 and neither shrank 4× from the one before it:
//      `roundoff: true`. On a lost-mode device's time step (7), only on a refined update, unless
//      its refinement failed three times in the call. Where a pivot was perturbed this iteration
//      (s.sys.perturbedPivots) or dt = Infinity, a convergence stands only if the block system's
//      residual (s.res; not the floating terminals' circuit rows nor the pins), reassembled and
//      scaled at the new state, is under 1e-6 in every row; else the call fails, saying which.
//   7. Refining (refineSolve): on from an iteration with a perturbed pivot, or once Newton stalls
//      (it ≥ 3, undamped, step > ¼ of the last and either under 1e-4 or grown twice running), or
//      throughout a lost-mode device's time steps. GMRES on the exact J·v (the dilute terms in
//      difference form, s._assembleDifference), preconditioned by the factorised solve, from that
//      solve's update; the refined update is taken where GMRES converged and moved it by under one
//      thermal unit, else refining stops: for the rest of the call, or on a lost-mode device for
//      two iterations, three times at most. A lost-mode device (s.lostMode, set for good): in a
//      time step, a refined correction that the next update kept (under a tenth of it), where the
//      correction was of a plain update under 1e-4, by more than its size and 1e-6, or a probe's.
//      A probe: on a device not yet lost-mode, a time step's 1st, 2nd, 4th, 8th… and every 32nd
//      convergence (counted in s.accepted) is followed by one refined update at the accepted
//      state, reassembled; over 1e-8 thermal units, it isn't accepted, and Newton goes on
//      refining.
//   8. Islands (s.steady and dt = Infinity, with islands): converged, the terminal currents off
//      Kirchhoff by more than 1e-9 of the largest and 1e-18 A/m², the islands' summed rows go in
//      (s.islandsOn) and Newton carries on, refining from the second iteration after; at most
//      two rounds, one where the call began with the rows in (`islands`). If it then fails, or
//      goes 12 iterations (or to maxIter) without converging, the solution it converged to before
//      that round is restored, with s.islandsOn as it was, and returned as converged.
//
// Returns { converged, iterations, history (each iteration's step), residual (the largest scaled
// residual before the last update; converged only), roundoff (true where 6 accepted at
// round-off), error (a message, on some failures) }.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. Converged: the last update met 6, and where 6 reassembles, every scaled residual row of the
//      block system at the returned state is under 1e-6.
//   2. Only u, uLo and termV are moved: by step 5, by 8's restore, and in a steady solve by each
//      assembly's snapFlat (src/lawRows.js); on return s.c and the coverages are those of the
//      returned state (computed there, or by the assembly at it), but after a SolverError from
//      computeConcentrations, which leaves them part done.
//   3. Failed: the state is where the last iteration left it (8's restore aside); the caller
//      restores it.
//   4. refineSolve leaves the residual, the terminal terms, the pins' residuals, s.transformed and
//      the factorised Jacobian as it found them, and the update as the plain solve gave it unless
//      it returns true.
//   5. Deterministic: the same state and options give the same iterates, bit for bit.
// Also written: s.islandsOn (as asked, or on by 8); s.lostMode, s.accepted (7); s.delta, s.deltaV,
// s.rowDiv; what an assembly writes, at the returned state where a probe or 6's check ran;
// s.conditioning (a steady solve's failure: the worst cancellation a factorisation showed, for
// its message); s.modelError, the message of a SolverError a statistics model threw at a state
// Newton reached (in the assembly, or computing the concentrations where it stops), which fails
// the call rather than leaving it (cleared by solveSteady and integrate, for their failures).
// Tests: every solve and step goes through here; test/newton.test.js (the islands round, the
// residual check's failure) and the MOS gate step in test/ports.test.js (a lost mode, to 1e-6);
// test/jacobian.test.js checks the system it solves. The devices that needed each branch are in
// its comments.

import { solveBordered } from './lawRows.js';
import { SolverError } from './solver.js';

// GMRES(m), restarted, for a real system op(x) = b, op already preconditioned; x holds the
// starting guess and receives the answer. Returns whether it converged.
function gmresReal(op, b, x, { m = 10, restarts = 3, tol = 1e-12 } = {}) {
  const N = b.length;
  const dot = (a, c) => {
    let t = 0;
    for (let j = 0; j < N; j++) t += a[j] * c[j];
    return t;
  };
  const bn = Math.sqrt(dot(b, b)) || 1;
  for (let cycle = 0; cycle < restarts; cycle++) {
    const ax = op(x), r = Float64Array.from(b, (v, j) => v - ax[j]);
    const beta = Math.sqrt(dot(r, r));
    if (!(beta > tol * bn)) return true;
    const V = [r.map((v) => v / beta)], H = Array.from({ length: m + 1 }, () => new Float64Array(m));
    const cs = [], sn = [], g = new Float64Array(m + 1);
    g[0] = beta;
    let k = 0, converged = false;
    for (; k < m; k++) {
      const w = op(V[k]);
      for (let i = 0; i <= k; i++) {
        const h = dot(V[i], w);
        H[i][k] = h;
        for (let j = 0; j < N; j++) w[j] -= h * V[i][j];
      }
      const hn = Math.sqrt(dot(w, w));
      V.push(w.map((v) => v / (hn || 1)));
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
      if (!(Math.abs(g[k + 1]) > tol * bn) || !(hn > 0)) {
        k++;
        converged = true;
        break;
      }
    }
    const y = new Float64Array(k);
    for (let i = k - 1; i >= 0; i--) {
      let v = g[i];
      for (let j = i + 1; j < k; j++) v -= H[i][j] * y[j];
      y[i] = v / H[i][i];
    }
    for (let i = 0; i < k; i++) for (let j = 0; j < N; j++) x[j] += y[i] * V[i][j];
    if (converged) return true;
  }
  return false;
}

// Scale every row by its largest Jacobian entry (in place, residual too), dividing by it: its
// reciprocal overflows where it's subnormal (a co-ion crowded out of a double layer, ~1e-310).
function equilibrate(s) {
  const { nB, sys, res, termB } = s, live = s._liveTerms();
  const { A, B, C, sizes, offA, offB, offC, offX } = sys;
  for (let b = 0; b < nB; b++) {
    const m = sizes[b], mp = b > 0 ? sizes[b - 1] : 0, mn = b < nB - 1 ? sizes[b + 1] : 0;
    for (let r = 0; r < m; r++) {
      const oa = offA[b] + r * mp, ob = offB[b] + r * m, oc = offC[b] + r * mn;
      let mx = 0;
      for (let k = 0; k < mp; k++) mx = Math.max(mx, Math.abs(A[oa + k]));
      for (let k = 0; k < m; k++) mx = Math.max(mx, Math.abs(B[ob + k]));
      for (let k = 0; k < mn; k++) mx = Math.max(mx, Math.abs(C[oc + k]));
      if (mx === 0) {
        s.rowDiv[offX[b] + r] = 1;
        continue;
      }
      const d = mx;
      s.rowDiv[offX[b] + r] = d;
      for (let k = 0; k < mp; k++) A[oa + k] /= d;
      for (let k = 0; k < m; k++) B[ob + k] /= d;
      for (let k = 0; k < mn; k++) C[oc + k] /= d;
      res[offX[b] + r] /= d;
      for (const k of live) termB[k][offX[b] + r] /= d; // (each terminal's ∂/∂V)
    }
  }
}

// A Newton solve refined by GMRES: the factorised Jacobian as the preconditioner, and J·x
// exact, with the dilute kernels' terms in difference form (see DifferenceTerms), as the
// impedance does. The assembled J holds a flux's dependence on η at its two ends as two
// entries; where they're huge (an inversion layer's 0.1 nm cells, a conductance ~1e11 against
// the layer's own storage and the trickle from the bulk), eliminating them loses the layer's
// overall level to round-off, and Newton rattles there. delta and deltaV come in as the plain
// solve's and leave refined (in the transformed unknowns, as solveBordered gives them); false
// if refining didn't help, leaving them as they came.
function refineSolve(s, dt, delta, deltaV) {
  const { res, termRes, termB, termC, termDI, rowDiv, floating: fl } = s;
  const N = s.sys.size, K = fl.length, pins = s.pins, R = s.rix;
  const save = { res: Float64Array.from(res), termRes: Float64Array.from(termRes), termB: termB.map((b) => Float64Array.from(b)), termC: termC.map((c) => Float64Array.from(c)), termDI: Float64Array.from(termDI), termDX: s.termDX.map((x) => Float64Array.from(x)), transformed: s.transformed, pinRes: s.pins.map((p) => p.res) };
  const restoreArrays = () => {
    res.set(save.res);
    termRes.set(save.termRes);
    save.termB.forEach((b, k) => termB[k].set(b));
    save.termC.forEach((c, k) => termC[k].set(c));
    termDI.set(save.termDI);
    save.termDX.forEach((x, k) => s.termDX[k].set(x));
    s.transformed = save.transformed;
    s.pins.forEach((p, a) => (p.res = save.pinRes[a]));
  };
  // The same system, its terms kept apart: the rest, and each floating terminal's (unscaled).
  const { rest, lin } = s._assembleDifference(dt);
  const Bk = fl.map((k) => Float64Array.from(termB[k])), Ck = fl.map((k) => Float64Array.from(termC[k])), Dk = fl.map((k) => termDI[k]);
  restoreArrays();
  const v = new Float64Array(N + 1), jv = new Float64Array(N);
  const op = (x) => {
    for (let j = 0; j < N; j++) v[j] = x[j];
    v[N] = 0;
    rest.multiply(v, jv);
    if (save.transformed) s._untransform(v); // (the dilute terms in the plain unknowns)
    lin.apply(v, jv, 1, 1 / dt);
    // A flat stretch's rows hold each level to its first node's (see flatRows in src/lawRows.js).
    for (const p of pins) {
      if (!p.col) continue;
      for (const o of p.col) if (o !== p.row && s.loc[o] >= 0) jv[R[o]] -= x[R[p.row]];
    }
    const out = new Float64Array(N + K);
    for (let j = 0; j < N; j++) {
      let t = jv[j];
      for (let a = 0; a < K; a++) t += Bk[a][j] * x[N + a];
      out[j] = t / rowDiv[j];
    }
    for (let a = 0; a < K; a++) {
      let t = Dk[a] * x[N + a];
      for (let b = 0; b < K; b++) if (b !== a) t += s.termDX[fl[a]][fl[b]] * x[N + b];
      for (let j = 0; j < N; j++) t += Ck[a][j] * x[j];
      out[N + a] = t;
    }
    // A pinned row is its constraint's (a conserved amount, an island's summed balance).
    for (const p of pins) {
      let t = 0;
      for (let j = 0; j < p.len; j++) t += p.w[j] * x[p.idx[j]];
      out[R[p.row]] = t;
    }
    return precondition(out);
  };
  const rhs = new Float64Array(N + 1), out = new Float64Array(N + 1), outV = new Float64Array(K);
  const precondition = (r) => {
    for (let j = 0; j < N; j++) rhs[j] = r[j];
    fl.forEach((k, a) => (termRes[k] = r[N + a]));
    for (const p of pins) {
      p.res = rhs[R[p.row]];
      rhs[R[p.row]] = 0;
    }
    solveBordered(s, rhs, out, outV, pins);
    const z = new Float64Array(N + K);
    for (let j = 0; j < N; j++) z[j] = out[j];
    for (let a = 0; a < K; a++) z[N + a] = outV[a];
    return z;
  };
  try {
    const b = new Float64Array(N + K);
    for (let j = 0; j < N; j++) b[j] = save.res[j];
    fl.forEach((k, a) => (b[N + a] = save.termRes[k]));
    pins.forEach((p, a) => (b[R[p.row]] = save.pinRes[a]));
    const pb = precondition(b), x = Float64Array.from(pb);
    // A lost mode or two, GMRES finds in a couple of iterations. The refined update is kept
    // where GMRES converged and its correction is under a thermal unit; otherwise newton()
    // stops refining.
    const x0 = Float64Array.from(x);
    const ok = gmresReal(op, pb, x, { m: 8, restarts: 1, tol: 1e-8 });
    // (Measured on the potentials, as Newton's own step is: a rate the plain solve left
    // unsolved, at a perturbed pivot, can need a correction of 41.7 mol/(m²·s).)
    const dx = new Float64Array(N + 1), dV = new Float64Array(K);
    for (let j = 0; j < N; j++) dx[j] = x[j] - x0[j];
    for (let a = 0; a < K; a++) dV[a] = x[N + a] - x0[N + a];
    if (s.transformed) s._untransform(dx);
    const c = maxPotentialStep(s, dx, dV);
    if (!ok || !(c < 1)) return false;
    // A time step's plain update near the solution corrected by more than its own size: either the
    // factorised solves lose a mode here (an inversion layer's level behind its 0.1 nm cells),
    // which neither their updates nor the residual show, or J·v's own round-off in a row of huge
    // entries made the correction (one ulp's worth, 2⁻¹³ or 2⁻¹⁹ exactly, on healthy devices).
    // newton() tells them apart by the next update: a real correction stays.
    s.lostCandidate = 0;
    if (dt !== Infinity && !s.lostMode && c > 1e-6) {
      const px = new Float64Array(N + 1), pV = new Float64Array(K);
      for (let j = 0; j < N; j++) px[j] = x0[j];
      for (let a = 0; a < K; a++) pV[a] = x0[N + a];
      if (s.transformed) s._untransform(px);
      const p = maxPotentialStep(s, px, pV);
      if (p < 1e-4 && c > p) s.lostCandidate = c;
    }
    for (let j = 0; j < N; j++) delta[j] = x[j];
    for (let a = 0; a < K; a++) deltaV[a] = x[N + a];
    return true;
  } finally {
    restoreArrays();
  }
}

// Each diffuse layer's drawn species' plane level, its update clipped to ±maxStep.
function clipLevels(s, delta, maxStep) {
  const R = s.rix;
  for (const d of s.diffList) {
    const bf = s.blockOfFace[d.f];
    for (const q of d.drawSlot) {
      if (q < 0) continue;
      const k = R[bf * s.M + q];
      if (Math.abs(delta[k]) > maxStep) delta[k] = Math.sign(delta[k]) * maxStep;
    }
  }
}

// Largest update among the potential-like unknowns (φ̂ and η at grid nodes; layers false leaves out
// the diffuse layers' drops and plane levels).
export function maxPotentialStep(s, delta, deltaV = s.deltaV, layers = true) { // (deltaV null: the device's potentials only)
  const { n, M } = s;
  const R = s.rix; // (the sink entry of delta is 0)
  let mx = 0;
  if (deltaV) for (let a = 0; a < s.floating.length; a++) mx = Math.max(mx, Math.abs(deltaV[a]) / s.VT);
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (s.nodeConductor[g] < 0) mx = Math.max(mx, Math.abs(delta[R[b * M]])); // (a metal's slot 0 is a flux)
    for (let i = 0; i < n; i++) {
      if (s.present[g * n + i]) mx = Math.max(mx, Math.abs(delta[R[b * M + 1 + i]]));
    }
    const k = s.surfAt[g]; // an electrode surface's coverages, as η
    if (k >= 0) for (let q = 0; q < s.model.surfaces[k].species.length; q++) mx = Math.max(mx, Math.abs(delta[R[b * M + 1 + n + q]]));
  }
  for (const gt of s.gateList) mx = Math.max(mx, Math.abs(delta[R[gt.o]])); // a gate's fraction open
  if (layers) for (const d of s.diffList) { // a diffuse layer's drop, and its drawn species' levels at its plane, as η
    const bf = s.blockOfFace[d.f];
    mx = Math.max(mx, Math.abs(delta[R[d.o]]));
    for (const q of d.drawSlot) if (q >= 0) mx = Math.max(mx, Math.abs(delta[R[bf * M + q]]));
  }
  return mx;
}

export function newton(s, dt, { maxIter = 60, tol = 1e-10, maxStep = 10, islands = false } = {}) {
  const { delta, res } = s;
  const deltaV = s.deltaV ?? (s.deltaV = new Float64Array(s.floating.length));
  const history = [], ownHistory = [];
  // Refined solves (see _refine), once the plain ones stall.
  let refine = false, refined = 0, noRefine = false;
  // A device whose plain solves lose a mode (see refineSolve): its time steps refined throughout,
  // a convergence taken only on a refined update. A refinement that fails is tried again two
  // iterations on (an early iterate's GMRES(8) can fail where a later one converges), three times.
  const lost = () => dt !== Infinity && s.lostMode === true;
  let lastRefined = false, failedAt = 0, failures = 0, candidate = 0, pending = 0;
  if (lost()) refine = true;
  // Islands (see applyIslands in src/lawRows.js): their summed rows go in only once the plain
  // solves converge.
  // Far from the solution, Newton exact along a weakly held level takes it by 1e5 thermal
  // units (a floating base, cold at bias), where the plain solves, which barely see that level,
  // converge. Converged, an island at the wrong level shows: it takes in a different current
  // than it passes on, so the terminal currents don't add up to zero. Then the summed rows go
  // in and Newton carries on from there. Converged again, a current that tiny against the
  // conduction around it rides on differences (~1e-28 thermal units for 5 fA/m² through a bulk
  // of 1e4 mol/m³ at 1e-4 m²/s) that the last update leaves known only to its own round-off:
  // where the currents still don't add up, one more update polishes them. Where the islands'
  // own system is held by static pivots elsewhere (a GaAs stack's minority carriers, at 1e-17
  // A/m² of noise), their rows can throw Newton off instead: if it doesn't converge again
  // within a dozen iterations, the plain solution stands.
  s.islandsOn = islands;
  let rounds = islands ? 1 : 0, kept = null, perturbed = false; // (rows in from the start: only the polishing round)
  // Converged by a small update, but solved at a perturbed pivot, or in a steady solve (where a
  // factorisation that lost a row's digits left a reaction's rate 3.6e7 off it, the update
  // small, an I = 0 electrode at a voltage 0.4 V from its steady one): the residual itself,
  // assembled at the state, must be small too (each row scaled by its largest Jacobian entry).
  // Returns why it isn't solved, or null.
  const unsolved = () => {
    if (!perturbed && dt !== Infinity) return null;
    try {
      s.assemble(dt);
    } catch (err) {
      if (!(err instanceof SolverError)) throw err;
      return err.message;
    }
    equilibrate(s);
    let r = 0;
    for (let k = 0; k < s.SINK; k++) r = Math.max(r, Math.abs(s.res[k]));
    if (r < 1e-6) return null;
    if (dt === Infinity) {
      // (the conditioning at the state it stopped at, the assembly having cleared the last one)
      try {
        s._factor();
      } catch {}
      noteConditioning(s);
    }
    return perturbed ? 'a row left unsolved at a pivot that cancelled exactly' : `its updates converged, but a row's residual is ${r.toExponential(1)} of its largest term`;
  };
  // One refined update at the state as it stands (see refineSolve), its size; NaN where it can't.
  const probe = () => {
    try {
      s.assemble(dt);
      equilibrate(s);
      try {
        s._factor();
      } catch {
        return NaN;
      }
      deltaV.fill(0);
      solveBordered(s, res, delta, deltaV, s.pins);
      if (!refineSolve(s, dt, delta, deltaV)) return NaN;
      if (s.transformed) s._untransform(delta);
      return maxPotentialStep(s, delta, deltaV);
    } catch (err) {
      if (!(err instanceof SolverError)) throw err;
      return NaN;
    }
  };
  // A lost mode that never stalls Newton: every plain update converges quadratically while the
  // mode keeps the predictor's error (a floating gate at rest in inversion, its electrons
  // reaching the back contact only by minority diffusion: 2.3e-8 thermal units at the 16th
  // step, 8e-5 V of drift ten steps on, as the extrapolation grew it). So a time step's 1st,
  // 2nd, 4th, 8th… and every 32nd acceptance on a device not yet known to lose one is probed
  // with a refined update at the accepted state: over the suite, under 1e-8 but for about a
  // dozen (round-off acceptances on strictly neutral and 1e-4 s steps). One over 1e-8 sends
  // Newton on, refining, and the next update decides (see candidate).
  const probed = () => {
    if (dt === Infinity || lost()) return false;
    const k = (s.accepted = (s.accepted ?? 0) + 1);
    if ((k & (k - 1)) !== 0 && k % 32 !== 0) return false;
    const d = probe();
    if (!(d > 1e-8)) return false;
    refine = true;
    noRefine = false;
    pending = d;
    return true;
  };
  const islandsOff = (it, result) => {
    if (dt !== Infinity || !s.steady || s.islands.length === 0 || rounds >= 2 || !(kirchhoff(s) > 1e-9)) return false;
    kept = { it, result, islandsOn: s.islandsOn, u: Float64Array.from(s.u), uLo: Float64Array.from(s.uLo), termV: Float64Array.from(s.termV) };
    s.islandsOn = true;
    rounds++;
    return true;
  };
  // A steady solve that fails: record how nearly singular the system was, and where; or, failing
  // with the islands' rows in, go back to the solution without them.
  // The concentrations at the state Newton stops at, where a statistics model can refuse it (a
  // Debye–Hückel solution past its range, at the converged state of a transient: the error was
  // thrown from here, uncaught, out of solve() and advance()). Its message, for the failure.
  const settle = () => {
    try {
      s.computeConcentrations();
      return null;
    } catch (err) {
      if (!(err instanceof SolverError)) throw err;
      s.modelError = err.message;
      return err.message;
    }
  };
  const fail = (r) => {
    if (kept) {
      s.u.set(kept.u);
      s.uLo.set(kept.uLo);
      s.termV.set(kept.termV);
      s.islandsOn = kept.islandsOn;
      s.computeConcentrations();
      return { ...kept.result, iterations: r.iterations, history };
    }
    if (dt === Infinity) noteConditioning(s);
    return r;
  };
  for (let it = 1; it <= maxIter; it++) {
    if (kept && it - kept.it > 12) return fail({ converged: false, iterations: it - 1, history });
    if (noRefine && failures < 3 && it - failedAt >= 2 && lost()) {
      noRefine = false;
      refine = true;
    }
    if (kept && it === kept.it + 2) {
      // Refined from the islands' second update (their first moves the level, perhaps by
      // more than a refinement may), since the plain solves have shown a mode they lose.
      refine = true;
      noRefine = false;
    }
    try {
      s.assemble(dt);
    } catch (err) {
      // A statistics model can fail far from the solution (e.g. Debye–Hückel beyond its range).
      if (!(err instanceof SolverError)) throw err;
      s.modelError = err.message;
      return fail({ converged: false, iterations: it, history, error: err.message });
    }
    equilibrate(s);
    let rmax = 0;
    for (let k = 0; k < s.SINK; k++) rmax = Math.max(rmax, Math.abs(res[k]));
    try {
      s._factor();
    } catch (err) {
      return fail({ converged: false, iterations: it, history, error: err.message });
    }
    // A pivot perturbed (one that cancelled to exactly nothing) can leave a row unsolved, its
    // update small while its residual isn't: the same couple written twice on an electrode at
    // open circuit converged with one spelling's rate row 41.7 mol/(m²·s) off, the couple
    // 1.5 RT from equilibrium. Refine the solves (exact J·v) wherever one was, and before
    // taking a small update for convergence there, look at the residual itself (see unsolved).
    perturbed = s.sys.perturbedPivots > 0;
    if (perturbed && !noRefine) refine = true;
    try {
      const pins = s.pins;
      deltaV.fill(0);
      solveBordered(s, res, delta, deltaV, pins);
      lastRefined = refine && refineSolve(s, dt, delta, deltaV);
      if (refine && !lastRefined) {
        refine = false; // (GMRES didn't converge, or strayed: the plain solves it is)
        noRefine = true;
        failedAt = it;
        failures++;
      }
      if (s.transformed) s._untransform(delta);
    } catch (err) {
      if (!(err instanceof SolverError)) throw err;
      return fail({ converged: false, iterations: it, history, error: err.message });
    }
    // Convergence counts every unknown, floating terminal voltages too. Damping and the
    // divergence check count only the device's own potentials: a terminal voltage enters
    // linearly (conductance links, held levels, a capacitive face), so a large swing in one is
    // safe, and damping it would throttle the whole update (a current-driven port switching
    // off, its voltage collapsing by hundreds of volts, crawled 10 thermal units at a time).
    // (So do the diffuse layers' drops, set by the potentials, and their drawn species' plane
    // levels, each level's own update clipped to maxStep instead: far from the solution a
    // reaction's rate can ask more of its layer than it can pass, and the level's update then
    // heads for −∞.)
    clipLevels(s, delta, maxStep);
    const step = maxPotentialStep(s, delta, deltaV);
    const own = maxPotentialStep(s, delta, null, false);
    // A refined correction the next update keeps (under a tenth of it) was a lost mode: this
    // device's time steps are refined from now on, a convergence taken only on a refined update.
    if (candidate > 0 && step < 0.1 * candidate) s.lostMode = true;
    if (lost() && !noRefine) refine = true;
    candidate = lastRefined ? pending || s.lostCandidate : 0;
    pending = 0;
    history.push(step);
    ownHistory.push(own);
    if (!Number.isFinite(step)) return fail({ converged: false, iterations: it, history, error: 'non-finite update' });
    // Give up early on clear divergence; the caller will take a smaller step instead.
    if ((it > 1 && own > 1e4) || (it > 6 && own > 10 * ownHistory[0])) {
      s.computeConcentrations();
      return fail({ converged: false, iterations: it, history, error: 'diverging' });
    }
    const alpha = own > maxStep ? maxStep / own : 1;
    addToState(s, delta, -alpha);
    s.floating.forEach((k, a) => (s.termV[k] -= alpha * deltaV[a]));
    // A time step's Newton contracts fast from its predictor (1e-2, 1e-9, 1e-15 thermal
    // units). The next update would be about θ = step/(the last one) times this one, so where
    // θ/(1 − θ) of it (Hairer and Wanner's estimate) is under tol, this update converged it:
    // the iteration that would only show that is skipped, a third of a typical step's work.
    const last = history[history.length - 2], theta = step / last;
    const trusted = !lost() || lastRefined || failures >= 3;
    const contracted = dt !== Infinity && !lost() && it >= 2 && step < 1e-6 && theta < 0.1 && (theta / (1 - theta)) * step < tol;
    if (alpha === 1 && trusted && (step < tol || contracted)) {
      const why = unsolved();
      if (why) return fail({ converged: false, iterations: it, history, error: why });
      const done = { converged: true, iterations: it, history, residual: rmax };
      if (probed()) continue;
      if (islandsOff(it, done)) continue;
      const refused = settle();
      if (refused) return fail({ converged: false, iterations: it, history, error: refused });
      return done;
    }
    // Converged as far as round-off allows: the updates are already tiny (below 1e-6 thermal
    // units, ~26 nV) and have stopped shrinking. A badly conditioned system's floor can sit
    // above tol: a strictly neutral material on a short step, where φ is fixed only through
    // fluxes that the storage term dwarfs, rattles at ~1e-9 after converging quadratically.
    // (Raising this floor lets a weakly held population drift: a MOS capacitor without a
    // channel port then showed a DC leak through its oxide.)
    const [p1, p2] = [history[history.length - 2], history[history.length - 3]], floor = 1e-6;
    // Stalled near the solution (not shrinking quadratically), or growing undamped twice
    // running: the factorised Jacobian may have lost a slow mode to cancellation, and
    // round-off's floor would accept it wrong, or its noise keep Newton from the solution (a
    // bipolar stack's floating base, its holes held ~1e14 more weakly than they move within
    // it). Refine the solves from here.
    if (refine) refined++;
    else if (!noRefine && alpha === 1 && it >= 3 && step > 0.25 * p1 && (step < 1e-4 || (step > p1 && p1 > p2))) refine = true;
    if (alpha === 1 && trusted && it >= 4 && (!refine || refined >= 1 || noRefine) && step < floor && p1 < floor && step > 0.25 * p1 && p1 > 0.25 * p2) {
      const why = unsolved();
      if (why) return fail({ converged: false, iterations: it, history, error: why });
      const done = { converged: true, iterations: it, history, residual: rmax, roundoff: true };
      if (probed()) continue;
      if (islandsOff(it, done)) continue;
      const refused = settle();
      if (refused) return fail({ converged: false, iterations: it, history, error: refused });
      return done;
    }
  }
  const refused = settle();
  return fail({ converged: false, iterations: maxIter, history, ...(refused ? { error: refused } : {}) });
}

// How far the terminal currents (as last assembled) are from adding up to zero, relative to the
// largest of them: Kirchhoff's law, which a steady state keeps.
function kirchhoff(s) {
  let sum = 0, mx = 0;
  for (const I of s.termI) {
    sum += I;
    mx = Math.max(mx, Math.abs(I));
  }
  // (currents adding to under 1e-18 A/m², a few charges per m² a second, are round-off: at
  // equilibrium every terminal's is, 1e-20 to 1e-31, and their sum as large as they are)
  return mx > 0 && Math.abs(sum) > 1e-18 ? Math.abs(sum) / mx : 0;
}

// How many digits the last factorisation lost to cancellation, and where (the worst seen since
// the steady solve began). Only for diagnosing failures: a saturated species, for instance,
// can lose digits harmlessly.
function noteConditioning(s) {
  let c;
  try {
    c = s.sys.cancellation();
  } catch {
    return;
  }
  if (!(c.digits > (s.conditioning?.digits ?? -1))) return;
  const { grid } = s.model;
  let x = NaN, where = '', nodes = [];
  for (let g = 0; g < s.nNodes; g++) {
    if (s.blockOfNode[g] !== c.block) continue;
    x = grid.x[g];
    where = s.model.regions[grid.nodeRegion[g]].name;
    nodes = [g];
  }
  for (let f = 0; f < s.nFaces; f++) {
    if (s.blockOfFace[f] !== c.block) continue;
    x = grid.x[grid.regionEnd[f]];
    where = `interfaces[${f}]`;
    nodes = [grid.regionEnd[f], grid.regionEnd[f] + 1];
  }
  s.conditioning = { digits: c.digits, x, where, nodes };
}

// u += scale·d in compensated arithmetic (Knuth two-sum, then renormalise hi/lo).
function addToState(s, d, scale) {
  const { u, uLo, fullOf } = s;
  for (let j = 0; j < fullOf.length; j++) {
    const k = fullOf[j], a = u[k], b = scale * d[j];
    const sum = a + b, bb = sum - a;
    const err = a - (sum - bb) + (b - bb);
    const lo = uLo[k] + err;
    const hi = sum + lo;
    u[k] = hi;
    uLo[k] = lo - (hi - sum);
  }
}
