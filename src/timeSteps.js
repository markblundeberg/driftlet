// Time steps: one implicit step of the device's equations (backward Euler, or variable-step
// BDF2), a step halved until it converges (advance), and adaptive stepping to a time with the
// local error controlled (integrate).
//
// CONTRACT
//
// These are functions of a Solver `s` (src/solver.js), which they read and write; its layout is
// in docs/internals.md. The Solver's own step, advance and integrate call them (the API, the
// steady solves' time stepping and the kit call those). What they use of it:
//   s.newton(dt, opts)  solves the equations at the step's end, from the state as it stands,
//                       with storage (c − c*)/dt against s.cOld, s.thOld, s.th0Old, s.gateOld,
//                       s.contactDStart, s.portQStart; leaves the state at its last iterate.
//   s._snapshot() / s._restore(snap)  the whole time state (u, uLo, termV, the start-of-step
//                       values, time, lastDt, history, the boundary intake, sourceTime).
//   s._assembleBookkeeping(dt), s.computeConcentrations(), s._accumulateBoundaryIntake(dt, c),
//   s._terminalCurrent(), s._probes(probes), s.fastestTime(), s.terms (drives, breakpoints).
//
// The state: u + uLo (Newton's unknowns, η and φ̂ in thermal units), termV (terminal voltages),
// time. The history: s.history, the accepted steps before this one, newest first, at most 3,
// each { t (its start), dt, u (its start state), c, th, th0, g (gates), D (contacts' D), Q
// (ports' charges) }: what BDF2 and the predictor need.
//
// step(s, dt, opts) → newton's result (+ bdf)
//   From a cold start (s.coldStart, outside a steady solve), the faces between strictly neutral
//   regions first set to their t = 0⁺ equilibrium (equilibrateFaces, src/start.js).
//   One step from s.time to tEnd = s.time + dt (or s.landing exactly, where it's within 1e-9·dt
//   of that). Backward Euler; BDF2 when opts.method is 'bdf2', a previous step exists and the
//   ratio ω = dt/dt_prev ≤ 2: then c* = b1·c_n − b2·c_{n−1} (b1 = (1+ω)/a0, b2 = ω²/((1+ω)a0),
//   a0 = (1+2ω)/(1+ω)) for every stored quantity, and the step is solved as backward Euler of
//   length dt/a0, which is BDF2's storage term. Sources are read at tEnd from within the step
//   (s.sourceTime, s.sourceBefore). If Newton fails: once more with the islands' rows in
//   throughout (a steady solve's dt = ∞ only), and once more with every circuit law read by its
//   change (see src/lawRows.js), each from the step's start. Accepted: the history gains the
//   step, s.time = tEnd, s.lastDt the effective dt, the boundary intake accumulates. Failed:
//   all as before the step (u, uLo, termV, the start-of-step values, atSteady, the sources'
//   time, the concentrations).
//
// advance(s, dt, opts, depth) → { converged, iterations, history, substeps }
//   step(dt); failing, two halves (each by advance), to a depth of 30.
//
// integrate(s, tEnd, opts) → { converged, done, steps, rejected, iterations, trace, stopped? }
//   Steps to the absolute time tEnd (an error if not a number; NaN, or a time not ahead of
//   s.time, does nothing and says why where it's likely a mistake). opts: tol (1e-3: the local
//   error allowed per step, in thermal units, as errorNorm measures it), dtMax, budgetMs (wall
//   time: stops early, done false, converged true), maxSteps, method ('bdf2' | 'be'), dt0,
//   probes. The step size starts from s.dtNext (the last call's), else dt0, else 1e-4 of the
//   interval, and is carried to the next call.
//   - The first step (no history): backward Euler, its error the difference from two half
//     steps, whose result is kept (and the first half step's point is in the trace).
//   - Later: step with Newton started from an extrapolation through the history; BDF2 once two
//     steps are in it; the error |u − pred|·scale against the explicit predictor (quadratic
//     after BDF2 with three, linear otherwise), scale making it the step's own local error.
//   - Accepted where error ≤ tol; the next h = h·factor, factor = 0.9·(tol/err)^(1/(p+1)) in
//     [0.2, 2] (2 at err = 0). Rejected: the state restored, h·factor. A step that doesn't
//     converge: h/4. Below the floor (1e-14 of tEnd, or 1e-2 of the fastest time scale), it
//     stops, converged false, with a message.
//   - A first step whose error doesn't shrink as it's cut (three tries within a factor 2 of
//     each other): a jump at a boundary, self-similar in x/√t, so taken anyway.
//   - Steps land exactly on the sources' breakpoints (and the history is cleared there: the
//     slope jumps) and on tEnd; within round-off of tEnd (1e-10 of it and 1e-9 of the
//     interval) the time snaps to it.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. The state only moves by accepted steps: a rejected step, a failed one, or one stopped
//      by the budget leaves u, uLo, termV, the start-of-step values, time and the history as
//      they were before it.
//   2. time advances by exactly the steps accepted, lands exactly on each breakpoint and on
//      tEnd, and never passes either.
//   3. BDF2's coefficients are those of variable-step BDF2 with ratio ω (exact for quadratics
//      in t), applied to every stored quantity alike (c, coverages, gates, contacts' D, ports'
//      charges), and only with ω ≤ 2 and a history entry of the previous step.
//   4. The error estimate is Milne's device as if the predictor's past states were exact: for
//      BDF2, |u − pred|·Cc/(Cc + Cp) (Cc, Cp the corrector's and predictor's error constants at
//      these step sizes); for backward Euler, h/(2h + h_prev). Through the integrator's own
//      states (which satisfy the method) the true local error is larger: 1 + 1/(2ω) times
//      backward Euler's estimate, and (2/9)/(2/11) ≈ 1.2 times BDF2's at equal steps, so a step
//      accepted at tol errs up to ~1.3 tol (a cold review, on an RC circuit). Known, and left so
//      far: backward Euler's low reading at short steps is what lets the controller step through
//      a state left off its algebraic equations (see ROADMAP). The exponent is 1/(p + 1).
//   5. The history holds only accepted steps of the same trajectory: newest first, at most 3,
//      cleared at a breakpoint, cut to the half step after a first step.
//   6. errorNorm counts every potential that carries state (φ̂ where it's defined, every
//      present species' η, an electrode surface's coverages as η, gates' fractions), in
//      thermal units, and nothing that doesn't.
// Tests: test/transient.test.js, and every test that calls advance(); tools/lawsweep.js (a
// settled transient as the steady state's reference).

import { nextBreakpoint } from './device.js';
import { SolverError } from './solver.js';
import { equilibrateFaces } from './start.js';

// Lagrange weights at x through the points s0, s1, s2.
function lagrange3(x, s0, s1, s2) {
  return [((x - s1) * (x - s2)) / ((s0 - s1) * (s0 - s2)), ((x - s0) * (x - s2)) / ((s1 - s0) * (s1 - s2)), ((x - s0) * (x - s1)) / ((s2 - s0) * (s2 - s1))];
}

/**
 * Advance by dt, splitting the interval into halves (recursively) wherever Newton fails.
 * Converged, the end time is t + dt; failed, wherever the halves that converged got to (the
 * substeps counted are those of the last attempt's halves).
 */
export function advance(s, dt, opts, depth = 0) {
  const r = s.step(dt, opts);
  if (r.converged || depth >= 30) return { ...r, substeps: 1 };
  const a = s.advance(dt / 2, opts, depth + 1);
  if (!a.converged) return a;
  const b = s.advance(dt / 2, opts, depth + 1);
  return {
    converged: b.converged,
    iterations: r.iterations + a.iterations + b.iterations,
    history: b.history,
    substeps: a.substeps + b.substeps,
  };
}

/**
 * One time step of size dt. Backward Euler by default; with `method: 'bdf2'`, variable-step
 * BDF2 once a previous step exists (and the step ratio is at most 2, for stability). On
 * failure the state is restored.
 *
 * BDF2's storage term (a0 c − (1+ω) c_n + ω²/(1+ω) c_{n−1})/dt is written in backward-Euler
 * form (c − c*)/(dt/a0), so assembly is shared: c* and the displacement histories replace
 * the start-of-step values, and dt/a0 replaces dt.
 */
export function step(s, dt, opts = {}) {
  // A transient's first step from a cold start begins at its t = 0⁺ (see equilibrateFaces).
  if (s.coldStart && !s.steady) equilibrateFaces(s);
  // What the step overwrites before Newton runs, for a failed one to put back: the Solution read
  // after it takes its contacts' fluxes against the start-of-step values (a host emptied by its
  // current read its contact's Li⁺ flux 8× off at the same state).
  const was = { cOld: Float64Array.from(s.cOld), thOld: Float64Array.from(s.thOld), th0Old: Float64Array.from(s.th0Old), gateOld: Float64Array.from(s.gateOld), contactDStart: s.contactDStart, portQStart: s.portQStart, atSteady: s.atSteady };
  s.atSteady = false;
  const prev = s.history[0];
  const w = prev ? dt / prev.dt : 0;
  const bdf = opts.method === 'bdf2' && prev !== undefined && w <= 2;
  s.uPrev.set(s.u);
  s.uPrevLo.set(s.uLo);
  s.termVPrev.set(s.termV);
  s.computeConcentrations();
  const cN = Float64Array.from(s.c), thN = Float64Array.from(s.th), th0N = Float64Array.from(s.th0);
  const gN = Float64Array.from(s.gateList, (gt) => s.u[gt.o] + s.uLo[gt.o]);
  // Contact displacement before the step: as it was at the end of the previous step (under
  // the parameters then), so a gate-voltage change shows up as displacement current.
  if (!s.contactDEnd || !s.portQEnd) {
    s._assembleBookkeeping(dt);
    s.contactDEnd = { ...s.contactD };
    s.portQEnd = Float64Array.from(s.portQ);
  }
  const DN = { ...s.contactDEnd }, QN = Float64Array.from(s.portQEnd);
  if (opts.guess) setGuess(s, opts.guess); // start Newton from a predicted state
  let dtEff = dt;
  if (bdf) {
    const a0 = (1 + 2 * w) / (1 + w), b1 = (1 + w) / a0, b2 = (w * w) / (1 + w) / a0;
    const { cOld } = s;
    for (let k = 0; k < cOld.length; k++) cOld[k] = b1 * cN[k] - b2 * prev.c[k];
    for (let k = 0; k < thN.length; k++) s.thOld[k] = b1 * thN[k] - b2 * prev.th[k];
    for (let k = 0; k < th0N.length; k++) s.th0Old[k] = b1 * th0N[k] - b2 * prev.th0[k];
    for (let k = 0; k < gN.length; k++) s.gateOld[k] = b1 * gN[k] - b2 * prev.g[k];
    s.contactDStart = { left: b1 * DN.left - b2 * prev.D.left, right: b1 * DN.right - b2 * prev.D.right };
    s.portQStart = QN.map((q, k) => b1 * q - b2 * prev.Q[k]);
    dtEff = dt / a0;
  } else {
    s.cOld.set(cN);
    s.thOld.set(thN);
    s.th0Old.set(th0N);
    s.gateOld.set(gN);
    s.contactDStart = DN;
    s.portQStart = QN;
  }
  // Implicit: sources at the step's end, as seen from within the step (before any jump there).
  // A step landing on a breakpoint ends on it exactly.
  const tEnd = s.landing !== undefined && Math.abs(s.time + dt - s.landing) <= 1e-9 * dt ? s.landing : s.time + dt;
  s.sourceTime = s.steady ? s.steadyAt : tEnd;
  s.sourceBefore = !s.steady;
  s.combining = true;
  s.lawForm = new Map(); // (see circuitLaws in src/lawRows.js)
  let result;
  try {
    result = s.newton(dtEff, opts);
    if (!result.converged && dtEff === Infinity && s.steady && s.islands.length > 0) {
      // Plain solves that don't converge may be lost along an island's level (two regions
      // conducting 1e20 times less hold it): once more from the start, with the islands'
      // summed rows in throughout (see newton).
      s.u.set(s.uPrev);
      s.uLo.set(s.uPrevLo);
      s.termV.set(s.termVPrev);
      if (opts.guess) setGuess(s, opts.guess);
      const again = s.newton(dtEff, { ...opts, islands: true });
      result = { ...again, iterations: result.iterations + again.iterations };
    }
    if (!result.converged && Number.isFinite(dtEff) && [...s.lawForm.values()].some((v) => !v)) {
      // A terminal read by its flux whose response is lost to round-off (a floating gate's
      // charging, below its metal's G Δη): once more, read by the change (see circuitLaws in
      // src/lawRows.js).
      s.u.set(s.uPrev);
      s.uLo.set(s.uPrevLo);
      s.termV.set(s.termVPrev);
      if (opts.guess) setGuess(s, opts.guess);
      for (const cs of s.lawForm.keys()) s.lawForm.set(cs, true);
      const again = s.newton(dtEff, opts);
      result = { ...again, iterations: result.iterations + again.iterations };
    }
  } finally {
    s.combining = false;
  }
  result.bdf = bdf;
  if (result.converged) {
    s.history.unshift({ t: s.time, dt, u: Float64Array.from(s.uPrev), c: cN, th: thN, th0: th0N, g: gN, D: DN, Q: QN });
    if (s.history.length > 3) s.history.length = 3;
    s.time = tEnd;
    s.lastDt = dtEff;
    s.contactDOld = s.contactDStart;
    s.portQOld = s.portQStart;
    s._accumulateBoundaryIntake(dtEff, cN);
    // (In time, the state is no longer the steady solve's that a continuation would start from.)
    if (!s.steady) {
      s.solvedV = undefined;
      s.stepped = true;
    }
  } else {
    s.u.set(s.uPrev);
    s.uLo.set(s.uPrevLo);
    s.termV.set(s.termVPrev);
    s.cOld.set(was.cOld);
    s.thOld.set(was.thOld);
    s.th0Old.set(was.th0Old);
    s.gateOld.set(was.gateOld);
    s.contactDStart = was.contactDStart;
    s.portQStart = was.portQStart;
    s.atSteady = was.atSteady;
    s.sourceTime = s.steady ? s.steadyAt : s.time;
    s.computeConcentrations();
  }
  s.sourceBefore = false;
  return result;
}

// Newton's unknowns set to a guess (e.g. extrapolated from the history). Only its unknowns: the
// other slots never change (see _activeSlots), and among them is a gauge, an insertion host's
// φ̂, that nothing pulls back. Extrapolated through steps that grow, it wandered off
// geometrically, 147 to 8000 thermal units in a battery's first 100 s and 1e11 by 1e5 s, where
// its species' levels (η − zφ̂, summed back by the host's statistics) lost their digits to it.
function setGuess(s, guess) {
  const { u, uLo } = s;
  for (const k of s.fullOf) {
    u[k] = guess[k];
    uLo[k] = 0;
  }
}

/**
 * Adaptive time stepping to tEnd: variable-step BDF2 (or backward Euler), with the local error
 * estimated against an explicit predictor through the previous states, and controlled to
 * `tol` thermal units per step in every potential (φ̂, each η, a surface's coverages as η). The first
 * step is checked by step doubling. Stops early when `budgetMs` of wall time is used, so an
 * animation can call it once per frame; the step size carries over between calls.
 */
export function integrate(s, tEnd, opts = {}) {
  // tEnd is an absolute time. Not a number at all is an error; a NaN (a first animation frame's
  // clock, or a time read from a field that doesn't exist) does nothing, and says so.
  const absolute = `tEnd is the time to reach in s, absolute (the device is at t = ${s.time}, which a solution gives as its time)`;
  if (typeof tEnd !== 'number') throw new SolverError(`advance(tEnd): ${absolute}; got ${tEnd}`);
  if (Number.isNaN(tEnd)) return { converged: true, done: false, steps: 0, rejected: 0, iterations: 0, trace: { t: [], current: [], voltage: [] }, stopped: `advance(NaN) did nothing: ${absolute}` };
  const { tol = 1e-3, dtMax = Infinity, budgetMs = Infinity, maxSteps = 100000, method = 'bdf2' } = opts;
  const clock = () => (globalThis.performance ? globalThis.performance.now() : Date.now());
  const start = clock();
  const sample = s._probes(opts.probes);
  const trace = { t: [], current: [], voltage: [], ...(sample ? { probes: sample.out } : {}) };
  let steps = 0, rejected = 0, iterations = 0, failed = false;
  s.modelError = null; // (a statistics model's refusal, for the failure: see newton)
  // Already there (an animation frame with no time to add): nothing to do, and the step size
  // carried to the next call stays as it was. A time already passed, or reached again outside
  // an animation, is most likely advance(dt) read as a duration: say so.
  if (!(tEnd > s.time)) {
    const relative = tEnd < s.time || (s.time > 0 && opts.budgetMs === undefined);
    return { converged: true, done: true, steps: 0, rejected: 0, iterations: 0, trace, ...(relative ? { stopped: `advance(${tEnd}) did nothing: ${absolute}. To run for a further dt, advance(sol.time + dt)` } : {}) };
  }
  const t0 = s.time, fromLayout = s.history.length === 0 && !s.stepped && !s.atSteady; // (c0's layout, not a solution)
  // A cold start whose floating double layers start uncharged far from the equilibrium of an ion
  // link beside them (a capacitive contact or port passing K⁺, its outside level set for the
  // charged steady state): the link closes the gap in the first instant, which can be far below
  // any step's reach. Noted here, for a failure to say what to do instead.
  const violent = s.history.length === 0 && !s.stepped ? violentStart(s) : [];
  let dt = s.dtNext ?? opts.dt0 ?? (tEnd - s.time) * 1e-4;
  if (!(dt > 0)) dt = (tEnd - s.time) * 1e-4;
  // The shortest step worth trying: round-off relative to the time, or the device's fastest
  // time scale, whichever is shorter (a cold start can need steps far below a long run's 1e-14).
  const floor = Math.min(1e-14 * Math.max(tEnd, 1e-300), s.fastestTime() * 1e-2);
  const pred = new Float64Array(s.u.length), guess = new Float64Array(s.u.length);
  const firstErrs = []; // the first step's error at each size tried (see below)
  const factor = (err, p) => (err > 0 ? Math.min(2, Math.max(0.2, 0.9 * Math.exp(Math.log(tol / err) / (p + 1)))) : 2);
  while (s.time < tEnd) {
    if (steps + rejected >= maxSteps || clock() - start > budgetMs) break;
    // Within round-off of the end (a target a hair past a breakpoint just landed on, say):
    // snap to it rather than attempt a step of ~1e-13 s, which can't be resolved and, in a
    // strictly neutral device, fails outright. Round-off of the time and negligible beside the
    // interval asked for: by the time alone, a GITT pulse of 0.01 s at t = 1e9 s was skipped
    // whole (the time moved, nothing else did), and from 1e10 s 17% of a pulse train's charge.
    if (tEnd - s.time <= 1e-10 * Math.abs(tEnd) && tEnd - s.time <= 1e-9 * (tEnd - t0)) {
      s.time = tEnd;
      break;
    }
    const remaining = tEnd - s.time;
    let h = Math.min(dt, dtMax);
    let clamped = h >= remaining * (1 - 1e-9);
    if (clamped) h = remaining;
    // Land on the waveforms' breakpoints, where a source's slope jumps.
    const tb = s.terms.reduce((m, t) => Math.min(m, nextBreakpoint(t.drive.src, s.time)), Infinity);
    const atBreak = tb - s.time <= h * (1 + 1e-9) && tb <= tEnd; // (one just past tEnd is the next call's)
    if (atBreak) {
      h = tb - s.time;
      clamped = true;
    }
    s.landing = atBreak ? tb : clamped ? tEnd : undefined;
    const snap = s._snapshot();
    let err, p, half = null;
    if (s.history.length > 0) firstErrs.length = 0;
    if (s.history.length === 0) {
      // First step: backward Euler, checked against two half steps (whose result is kept).
      const full = s.step(h);
      iterations += full.iterations;
      if (full.converged) {
        const uFull = Float64Array.from(s.u);
        s._restore(snap);
        const a = s.step(h / 2);
        // The first half step is kept too, so it goes in the trace (a jump's charge is in it).
        half = a.converged ? { t: s.time, current: s._terminalCurrent(), voltage: s.termV[1] - s.termV[0], probes: sample?.read() } : null;
        const b = a.converged ? s.step(h / 2) : a;
        iterations += a.iterations + (b === a ? 0 : b.iterations);
        if (b.converged) {
          err = errorNorm(s, uFull);
          p = 1;
        }
      }
      if (err === undefined) {
        s._restore(snap);
        rejected++;
        dt = h / 4;
        if (dt < floor) { failed = true; break; }
        continue;
      }
    } else {
      // (Order 2 once three consistent states can check it.)
      const r = s.step(h, { method: s.history.length >= 2 ? method : 'be', guess: extrapolate(s, guess, h) });
      iterations += r.iterations;
      if (!r.converged) {
        s._restore(snap);
        rejected++;
        dt = h / 4;
        if (dt < floor) { failed = true; break; }
        continue;
      }
      p = predict(s, pred, r.bdf);
      err = errorNorm(s, pred) * p.scale;
      p = p.order;
    }
    // A first step whose error doesn't shrink with it: a jump in a level held at a boundary
    // (a bath's, after a voltage step) starts a profile self-similar in x/√t, so a first step
    // of any length errs alike, down to the grid's own diffusion time, below Newton's reach.
    // Three tries in a row, each within a factor 2 of the last, say so: the error is the
    // jump's, not the step's. Take the step (backward Euler damps it as it should); the next is
    // sized from its error as any other. (Counted within one call: under a frame's budget, a
    // first step cut across calls isn't seen to stall.)
    if (half && err > tol) firstErrs.push(err);
    const jump = half && firstErrs.length >= 3 && firstErrs.slice(-3).every((e, j, a) => j === 0 || (e > 0.5 * a[j - 1] && e < 2 * a[j - 1]));
    if (err > tol && !jump) {
      s._restore(snap);
      rejected++;
      dt = h * factor(err, p);
      continue;
    }
    steps++;
    // The state the first step started from needn't satisfy the algebraic equations (φ in a
    // strictly neutral region, an interface's unknowns, after a start or a jump), and they
    // jump in the first instant; a predictor through it would see that jump as error on every
    // step after. The history starts from the half step.
    if (half) s.history.length = 1;
    dt = clamped ? Math.max(dt, h * factor(err, p)) : h * factor(err, p);
    if (atBreak) s.history = []; // the solution's slope jumps here: restart the order
    if (half) {
      trace.t.push(half.t);
      trace.current.push(half.current);
      trace.voltage.push(half.voltage);
      if (sample) sample.push(half.probes);
    }
    trace.t.push(s.time);
    trace.current.push(s._terminalCurrent());
    trace.voltage.push(s.termV[1] - s.termV[0]);
    if (sample) sample.push(sample.read());
  }
  s.dtNext = dt;
  s.landing = undefined;
  // Stopped short, other than on the caller's wall-time budget: say where, and why.
  let stopped;
  const at = `advance stopped at t = ${s.time.toPrecision(6)} s of ${tEnd.toPrecision(6)}`;
  // (What usually does it, most likely first. Round 11's agents read the old list, which led
  // with a drive's jump, as blaming a modest current pulse: a host emptied at its surface, and
  // an electrode at open circuit with no double layer, held only by a slow reaction.)
  const unbalanced = s.history.length === 0 && s.time === t0 ? (s.unbalanced ?? []).map((r) => s.model.regions[r].name ?? `regions[${r}]`) : [];
  if (failed && unbalanced.length > 0 && s.time === t0)
    stopped = `${at}: the start can't be neutral in ${unbalanced.join(', ')}, strictly neutral (ε = 0): the amounts its c0 gives leave a net charge that no species held at a contact's level can balance, and a step can't move a stored amount at once. Give c0 that balances (with any fixed charge), or leave out of c0 the species a contact feeds, so that they start at its level.`;
  else if (failed && violent.length > 0)
    stopped = `${at}: a step did not converge even when shortened to ${dt.toExponential(1)} s. The start is violent: ${violent.join('; ')}. A double layer starts uncharged, and what passes beside it closes such a gap in the first instant, faster than any step can follow. For a start in equilibrium, solve() first (the steady state at the starting drives), then change what the experiment changes, then advance(). Where solve() fails from there too, hold that terminal at a voltage, solve(), then let it float (set its I: 0): it keeps the state the held one gave it.`;
  else if (failed && s.modelError)
    stopped = `${at}: a step did not converge even when shortened to ${dt.toExponential(1)} s, the state it reached refused by a statistics model: ${s.modelError}. The model doesn't hold there: a smaller concentration or drive, or a model for it (lattice for crowding).`;
  else if (failed && fromLayout && steps === 0)
    // (round 15's fuzzer: the pn template advanced from its c0, a step of 3e-24 s failing, 23 of 62
    // plain diodes; nothing said to solve first)
    stopped = `${at}: the first step did not converge even when shortened to ${dt.toExponential(1)} s, from the starting state, which is c0's layout rather than a solution (a junction's depletion not yet formed, say). For a transient from equilibrium, solve() first at the starting drives, then change what the experiment changes, then advance().`;
  else if (failed) stopped = `${at}: a step did not converge even when shortened to ${dt.toExponential(1)} s. Look for: a host or surface driven to its limit (emptied or filled at its surface by a current it can no longer pass); an electrode at open circuit held only by a slow reaction, with no double-layer capacitance to keep its level (give its face one: phi capacitive); a drive's jump too big for one step (ramp it); a driven current the device cannot keep passing (a diode's reverse recovery past its storage time: stop there); or a device the steady solves also find hard.`;
  else if (s.time < tEnd && steps + rejected >= maxSteps) stopped = `${at}: maxSteps (${maxSteps}) used up (${rejected} of them rejected). Call again to go on, or loosen tol.`;
  return { converged: !failed, done: s.time >= tEnd, steps, rejected, iterations, trace, ...(stopped ? { stopped } : {}) };
}

// Quadratic (or linear) extrapolation of the state a time h ahead, from the current state and
// the two most recent history entries: a starting guess for Newton. (In offsets built from the
// steps' own lengths, never differences of absolute times: a second into a run, after a jump
// that a fine grid resolves in steps of 1e-13 s, those differences keep only 3 digits.)
function extrapolate(s, out, h) {
  const [e0, e1] = s.history;
  const u = s.u;
  if (!e0) return undefined;
  if (!e1) {
    const l = h / e0.dt;
    for (let k = 0; k < out.length; k++) out[k] = u[k] + l * (u[k] - e0.u[k]);
    return out;
  }
  const [l0, l1, l2] = lagrange3(h, 0, -e0.dt, -(e0.dt + e1.dt));
  for (let k = 0; k < out.length; k++) out[k] = l0 * u[k] + l1 * e0.u[k] + l2 * e1.u[k];
  return out;
}

// Explicit predictor at the new time through the previous states (quadratic after BDF2,
// linear after backward Euler), and the factor turning |u − pred| into the local error.
function predict(s, pred, bdf) {
  const [e0, e1, e2] = s.history;
  const h = e0.dt;
  if (!e1) {
    pred.set(e0.u);
    return { scale: 1, order: 1 };
  }
  const hp = e1.dt;
  if (bdf && e2) {
    const hpp = e2.dt;
    // Lagrange through the three states before at the new time (offsets as in extrapolate).
    const [l0, l1, l2] = lagrange3(0, -h, -(h + hp), -(h + hp + hpp));
    for (let k = 0; k < pred.length; k++) pred[k] = l0 * e0.u[k] + l1 * e1.u[k] + l2 * e2.u[k];
    const w = h / hp;
    const Cc = (h * h * h * (1 + w) * (1 + w)) / (w * (1 + 2 * w));
    const Cp = h * (h + hp) * (h + hp + hpp);
    return { scale: Cc / (Cc + Cp), order: 2 };
  }
  const l0 = (h + hp) / hp, l1 = 1 - l0;
  for (let k = 0; k < pred.length; k++) pred[k] = l0 * e0.u[k] + l1 * e1.u[k];
  // Backward Euler's error against a linear predictor; after a BDF2 step without enough
  // history this overestimates, which is safe.
  return { scale: bdf ? 1 : h / (2 * h + hp), order: 1 };
}

// Largest difference from ref among the potentials that carry state, in thermal units.
function errorNorm(s, ref) {
  const { n, M, u } = s;
  let mx = 0;
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (!s.phiUndefined[g]) mx = Math.max(mx, Math.abs(u[b * M] - ref[b * M]));
    for (let i = 0; i < n; i++) {
      if (s.present[g * n + i]) mx = Math.max(mx, Math.abs(u[b * M + 1 + i] - ref[b * M + 1 + i]));
    }
    const k = s.surfPort[g]; // an electrode surface's coverages, as η
    if (k >= 0) for (let q = 0; q < s.model.ports[k].surface.length; q++) mx = Math.max(mx, Math.abs(u[b * M + 1 + n + q] - ref[b * M + 1 + n + q]));
  }
  for (const { o } of s.gateList) mx = Math.max(mx, Math.abs(u[o] - ref[o])); // (a fraction: 1e-3 of m³ is 0.3% of P)
  return mx;
}

// The floating double layers a cold start leaves far (0.1 V or more) from the equilibrium of a
// path beside them: a contact with a capacitive φ, or a port with a capacitance, at I = 0 or
// behind a resistance, with a link holding a charged species at an outside level V + offset,
// or (a port) a reaction that takes its carrier. Each as words for a message.
export function violentStart(s) {
  const out = [], { M, n, z, VT, model } = s;
  const check = (what, k, links, g) => {
    for (let i = 0; i < n; i++) {
      const l = links[i];
      if (z[i] === 0 || !['conductance', 'equilibrium'].includes(l.type) || !s.present[g * n + i]) continue;
      const gap = (VT * s.u[s.blockOfNode[g] * M + 1 + i]) / z[i] - (s.termV[k] + (l.offset ?? 0));
      if (Math.abs(gap) >= 0.1) out.push(`${what}'s ${model.species[i].name} link is ${Math.abs(gap).toFixed(2)} V from equilibrium`);
    }
  };
  for (const k of s.floating) {
    const t = s.terms[k];
    if (t.kind === 'contact') {
      const ct = model.contacts[t.side];
      if (ct.phi.type === 'capacitive') check(`contact ${t.side}`, k, ct.species, t.side === 'left' ? 0 : s.nNodes - 1);
      continue;
    }
    const port = model.ports[t.index];
    if (!port.capacitance) continue;
    check(`port ${port.name}`, k, port.species, port.nodes[0]);
    // Its reactions at its first node: the affinity a = fixedA − Σ ν η, in volts by how it moves
    // with the port's voltage (−Σ ν z/V_T over the carrier's part).
    const b = s.blockOfNode[port.nodes[0]];
    port.reactions.forEach((rx, j) => {
      let a = rx.fixedA, slope = 0;
      for (const p of rx.part) {
        if (p.side === 1) {
          a -= p.nu * s.portEta(port, p.i);
          slope -= (p.nu * z[p.i]) / VT;
        } else a -= p.nu * s.u[b * M + (p.side === 2 ? 1 + n + p.s : 1 + p.i)];
      }
      const gap = slope !== 0 ? Math.abs(a / slope) : 0;
      if (gap >= 0.1) out.push(`port ${port.name}'s reactions[${j}] is ${gap.toFixed(2)} V from equilibrium`);
    });
  }
  return out;
}
