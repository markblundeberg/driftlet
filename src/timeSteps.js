// Time steps: one implicit step of the device's equations (backward Euler, or variable-step
// BDF2), a step halved until it converges (advance), and adaptive stepping to a time with the
// local error controlled (integrate).
//
// The names shift by one from the API's: Device.step(dt) calls advance here (a duration, halved
// on failure), Device.advance(tEnd) calls integrate (an absolute time), and step here is a single
// attempt. integrate's messages say advance(), the name its user called.
//
// CONTRACT
//
// These are functions of a Solver `s` (src/solver.js), which they read and write; its layout is
// in docs/internals.md. The Solver's own step, advance and integrate call them (the API, the
// steady solves' time stepping and the kit call those). What they use of it:
//   s.newton(dt, opts)  solves the equations at the step's end, from the state as it stands,
//                       with storage (c − c*)/dt against s.cOld, s.thOld, s.th0Old, s.gamOld (the
//                       diffuse layers' excesses), s.gateOld,
//                       s.contactDStart, s.portQStart; leaves the state at its last iterate.
//   s.jumpedT           [T before, T after] where set() changed T and no step has been taken
//                       since (cleared by an accepted step or a converged solve), for the
//                       first step's failure to name.
//   s._snapshot() / s._restore(snap)  the whole state between calls (STATE in src/solver.js: u,
//                       uLo, termV, termCharge + termChargeLo, the start-of-step values, time, lastDt, history,
//                       the boundary intake, sourceTime, and stepped, solvedV and jumpedT).
//   s._assembleBookkeeping(dt), s.computeConcentrations(), s._accumulateBoundaryIntake(dt, c),
//   s._terminalCurrent(), s._terminalCurrents(), s._probes(probes), s.fastestTime(), s.terms
//   (drives, breakpoints).
//
// The state: u + uLo (Newton's unknowns, η and φ̂ in thermal units), termV (terminal voltages),
// time, termCharge + termChargeLo (each terminal's charge passed into the device in time steps,
// and its rounding). The history: s.history, the accepted steps before this one, newest first, at
// most 3, each { t (its start), dt, u (its start state), V (its start's terminal voltages), c, th,
// th0, gam, g (gates), D (contacts' D), Q (ports' charges), P + PLo (termCharge), and for the
// charges' error: o and t0 (P, PLo and the time where the history began), Y ({ c, y }: each
// terminal's charge since then, whole and less what its contact or port stores) and G (its gross,
// Σ|ΔY.c|) }: what BDF2 and the predictor need.
//
// step(s, dt, opts) → newton's result (+ bdf)
//   From a cold start (s.coldStart, outside a steady solve), the faces between strictly neutral
//   regions first set to their t = 0⁺ equilibrium (equilibrateFaces, src/start.js).
//   One step from s.time to tEnd = s.time + dt (or s.landing exactly, where it's within 1e-9·dt or 4 ulps
//   of that). Backward Euler; BDF2 when opts.method is 'bdf2', a previous step exists and the
//   ratio ω = dt/dt_prev ≤ 2: then c* = b1·c_n − b2·c_{n−1} (b1 = (1+ω)/a0, b2 = ω²/((1+ω)a0),
//   a0 = (1+2ω)/(1+ω)) for every stored quantity, and the step is solved as backward Euler of
//   length dt/a0, which is BDF2's storage term. Sources are read at tEnd from within the step
//   (s.sourceTime, s.sourceBefore). If Newton fails: once more with the islands' rows in
//   throughout (a steady solve's dt = ∞ only), and once more with every circuit law read by its
//   change (see src/lawRows.js), each from the step's start. Accepted: the history gains the
//   step, s.time = tEnd, s.lastDt the effective dt, the boundary intake accumulates, and outside
//   a steady solve each terminal's charge steps by its current as a stored charge does, P = P* +
//   I·dt/a0 (P* from P as c* from c; backward Euler P_n + I·dt; compensated, its rounding in
//   termChargeLo), so that where the device stores it a trace's charge is the stored charge's
//   change exactly, and Kirchhoff's sum holds in it as in the currents. Failed:
//   all as before the step, by s._restore() of a snapshot taken on entry (the whole of STATE,
//   and the concentrations from it), as every other rollback. A step that can't move the clock
//   (s.time + dt === s.time, outside a steady solve) fails untouched.
//
// advance(s, dt, opts, depth) → { converged, iterations, history, substeps }
//   step(dt); failing, two halves (each by advance), to a depth of 30.
//
// integrate(s, tEnd, opts) → { converged, done, steps, rejected, iterations, trace, stopped? }
//   Steps to the absolute time tEnd (an error if not a number; NaN, or a time not ahead of
//   s.time, does nothing and says why where it's likely a mistake). opts: tol (1e-3: the local
//   error allowed per step, in thermal units, as errorNorm measures it), dtMax, budgetMs (wall
//   time: stops early, done false, converged true), maxSteps, method ('bdf2' | 'be'), dt0,
//   probes. What the arguments alone refuse (tEnd ±∞ or not a number, tol, dtMax, maxSteps, a
//   method not 'bdf2' or 'be', a probe) is refused by integrateArgs(s, tEnd, opts), which the
//   API asks before anything is solved. The step size starts from dt0 where given, else
//   s.dtNext (the last call's), else 1e-4 of the interval, and is carried to the next call by a
//   call that attempts a step (one that attempts none leaves it, and the history, as they were).
//   - The first step (no history): backward Euler, its error the difference from two half
//     steps, whose result is kept (and the first half step's point is in the trace).
//   - Later: step with Newton started from an extrapolation through the history; BDF2 once two
//     steps are in it; the error |u − pred|·scale against the explicit predictor (quadratic
//     after BDF2 with three, linear otherwise), scale making it the step's own local error.
//     errorNorm is over φ̂, each η, a surface's coverages, the gates and the floating terminals'
//     voltages (over V_T), the largest of them.
//   - And each terminal's charge, every terminal's, relative (chargeError, chargeDrift): the
//     pair's two quadratures, or the charge against the predictor's (as u), over the larger of
//     the terminal's gross since the history began and (the first step) what its current would
//     pass in the segment's time left; differences below what the charges are known to (an
//     absolute 1e-12 C/m², 1e-15 A/m² over the segment, 1e-12 of what the terminal's boxes hold)
//     are noise. The step's error is the larger of the two.
//   - Accepted where error ≤ tol; the next h = h·factor, factor = 0.9·(tol/err)^(1/(p+1)) in
//     [0.2, 2] (2 at err = 0). Rejected: the state restored, h·factor. A step that doesn't
//     converge: h/4. Below the floor (1e-14 of tEnd, or 1e-2 of the fastest time scale), it
//     stops, converged false, with a message; past the first step, an error rejection whose
//     step can't move the clock (s.time + dt === s.time) does too.
//   - A first step whose error doesn't shrink as it's cut (three tries within a factor 2 of
//     each other): a jump at a boundary, self-similar in x/√t, so taken anyway (the potentials'
//     error only; a charge's must still be within tol).
//   - Steps land exactly on the sources' breakpoints (and the history is cleared there: the
//     slope jumps) and on tEnd; within round-off of tEnd (1e-12 of it and 1e-9 of the
//     interval) the time snaps to it. A step that ends, or a call that starts, within
//     round-off of a breakpoint (1e-12 of the time, where nextBreakpoint() may count it passed)
//     clears the history too: the next step crosses that breakpoint unlanded.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. The state only moves by accepted steps: a rejected step, a failed one, or one stopped
//      by the budget leaves u, uLo, termV, the start-of-step values, time and the history as
//      they were before it.
//   2. time advances by exactly the steps accepted, lands exactly on each breakpoint and on
//      tEnd, and never passes either.
//   3. BDF2's coefficients are those of variable-step BDF2 with ratio ω (exact for quadratics
//      in t), applied to every stored quantity alike (c, coverages, gates, contacts' D, ports'
//      charges), and only with ω ≤ 2 and a history entry of the previous step; but for the
//      gates' target, held to [0, 1], so that a gate's fraction stays there (docs/device.md).
//   4. The error estimate is Milne's device as if the predictor's past states were exact: for
//      BDF2, |u − pred|·Cc/(Cc + Cp) (Cc, Cp the corrector's and predictor's error constants at
//      these step sizes); for backward Euler, h/(2h + h_prev). Through the integrator's own
//      states (which satisfy the method) the true local error is larger: 1 + 1/(2ω) times
//      backward Euler's estimate, and (2/9)/(2/11) ≈ 1.2 times BDF2's at equal steps, so a step
//      accepted at tol errs up to ~1.3 tol (a cold review, on an RC circuit). Known, and left so
//      far: backward Euler's low reading at short steps is what lets the controller step through
//      a state left off its algebraic equations (see ROADMAP). The exponent is 1/(p + 1).
//   5. The history holds only accepted steps of the same trajectory: newest first, at most 3,
//      cleared at a breakpoint (landed on, or within round-off of a step's end or a call's start), cut to the
//      half step after a first step. So no BDF2 step spans a source's jump.
//   6. errorNorm counts every potential that carries state (φ̂ where it's defined, every
//      present species' η, an electrode surface's coverages as η, gates' fractions), in
//      thermal units, and nothing that doesn't; every terminal's charge is in the step's error
//      beside it, relative, so that a path through which a charge flows without being stored is
//      held to tol too.
// Tests: test/transient.test.js, and every test that calls advance(); tools/lawsweep.js (a
// settled transient as the steady state's reference).

import { nextBreakpoint, nearBreakpoint } from './device.js';
import { SolverError } from './errors.js';
import { FARADAY } from './constants.js';
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
  // At a source's breakpoint, BDF2's history spans the kink in what it stores, and a current
  // switched on there lost I·dt/2 of charge for good (Fable round 4's review of this module):
  // the order restarts there, as integrate's does.
  // At every depth: a half that starts on the breakpoint, after its first half failed and its
  // second converged, crossed it with BDF2 and lost a third of its charge (Sol's third review).
  if (s.history.length > 0 && s.terms.some((t) => nearBreakpoint(t.drive.src, s.time))) s.history = [];
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
  // A step that can't move the clock isn't taken: it fails, as before it, untouched. (A first
  // step of 1e-18 s at t = 1 s posted its charge with the clock unmoved, and a step of one ulp
  // halved on a failure posted one too; Sol's re-read of the time integrator, round 7.)
  if (!s.steady && s.time + dt === s.time) return { converged: false, iterations: 0, history: [], bdf: false };
  // A failed step puts back the whole state from before it (STATE in src/solver.js), by the one
  // way every rollback takes. (By a hand-written list: a failed cold one left the projection, and
  // coldStart cleared, behind, Sol's third review; the list then missed what newton's islands'
  // round writes, and a failed step() from a warned floating base lost its "doesn't balance"
  // warning, and the sources' time, rebuilt from the clock after a snap past an edge, read a
  // blocking plate's 0.01 C/m² that never passed; round 7's S15. And the start-of-step values:
  // the Solution read after it takes its contacts' fluxes against them, a host emptied by its
  // current read its contact's Li⁺ flux 8× off at the same state.)
  const snap = s._snapshot();
  // A transient's first step from a cold start begins at its t = 0⁺ (see equilibrateFaces).
  if (s.coldStart && !s.steady) equilibrateFaces(s);
  s.atSteady = false;
  const prev = s.history[0];
  const w = prev ? dt / prev.dt : 0;
  const bdf = opts.method === 'bdf2' && prev !== undefined && w <= 2;
  s.uPrev.set(s.u);
  s.uPrevLo.set(s.uLo);
  s.termVPrev.set(s.termV);
  s.computeConcentrations();
  const cN = Float64Array.from(s.c), thN = Float64Array.from(s.th), th0N = Float64Array.from(s.th0), gamN = Float64Array.from(s.gam);
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
  const a0 = (1 + 2 * w) / (1 + w), b2 = (w * w) / (1 + w) / a0;
  // (As x_n + b2 (x_n − x_{n−1}), which is b1 x_n − b2 x_{n−1} since b1 − b2 = 1: exact on a
  // state at rest, where that form drifted by ε a step, 3.7e-11 of a closed film's amount over
  // 81,000 steps and ε·|P| of a terminal's charge a step; Sol's third review.)
  if (bdf) {
    const { cOld } = s;
    for (let k = 0; k < cOld.length; k++) cOld[k] = cN[k] + b2 * (cN[k] - prev.c[k]);
    for (let k = 0; k < thN.length; k++) s.thOld[k] = thN[k] + b2 * (thN[k] - prev.th[k]);
    for (let k = 0; k < th0N.length; k++) s.th0Old[k] = th0N[k] + b2 * (th0N[k] - prev.th0[k]);
    for (let k = 0; k < gamN.length; k++) s.gamOld[k] = gamN[k] + b2 * (gamN[k] - prev.gam[k]);
    // A gate's target held to [0, 1], a fraction's range: then the step's x = (x* + h α)/(1 +
    // h(α + β)), with α, β ≥ 0 and h = dt/a0, is in [0, 1] too. (Extrapolated past it, BDF2's
    // x* overshoots a fast gate's swing where h(α + β) > 1/2, its roots complex: m went to
    // −2.7e-2 on a clamp from +50 to −100 mV at 6.3 °C in steps of 0.1 ms, where the exact
    // relaxation is at 7.7e-4; round 7's S6. Inside [0, 1], it's BDF2's own.)
    for (let k = 0; k < gN.length; k++) s.gateOld[k] = Math.min(1, Math.max(0, gN[k] + b2 * (gN[k] - prev.g[k])));
    s.contactDStart = { left: DN.left + b2 * (DN.left - prev.D.left), right: DN.right + b2 * (DN.right - prev.D.right) };
    s.portQStart = QN.map((q, k) => q + b2 * (q - prev.Q[k]));
    dtEff = dt / a0;
  } else {
    s.cOld.set(cN);
    s.thOld.set(thN);
    s.th0Old.set(th0N);
    s.gamOld.set(gamN);
    s.gateOld.set(gN);
    s.contactDStart = DN;
    s.portQStart = QN;
  }
  // Implicit: sources at the step's end, as seen from within the step (before any jump there).
  // A step landing on a breakpoint ends on it exactly.
  // (Within a few ulps too: at 1e-9 of a short step, below an ulp of the time, a startup pair's
  // second half ended an ulp off its target, and the call then snapped the clock past the
  // sources' time; Sol's third review. But only a step that ends nearer the landing than half of
  // itself: the pair's first half of a step of two ulps snapped there too, and the second half
  // stepped no time at all, two rows at one time; Sol's re-read, round 7.)
  const near = s.landing !== undefined && Math.abs(s.time + dt - s.landing) <= Math.min(Math.max(1e-9 * dt, 4 * Number.EPSILON * Math.abs(s.landing)), dt / 2);
  const tEnd = near ? s.landing : s.time + dt;
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
    const changeable = (cs) => Number.isFinite(cs.lambda) && cs.lambda !== 0; // (one that counts a current: λ set)
    if (!result.converged && Number.isFinite(dtEff) && [...s.lawForm].some(([cs, v]) => !v && changeable(cs))) {
      // A terminal read by its flux whose response is lost to round-off (a floating gate's
      // charging, below its metal's G Δη): once more, read by the change (see circuitLaws in
      // src/lawRows.js).
      s.u.set(s.uPrev);
      s.uLo.set(s.uPrevLo);
      s.termV.set(s.termVPrev);
      if (opts.guess) setGuess(s, opts.guess);
      // (not a law that counts no current, λ unset: read by its change, it would be NaN, as in
      // the impedance's byChange)
      for (const cs of s.lawForm.keys()) if (changeable(cs)) s.lawForm.set(cs, true);
      const again = s.newton(dtEff, opts);
      result = { ...again, iterations: result.iterations + again.iterations };
    }
  } finally {
    s.combining = false;
  }
  result.bdf = bdf;
  if (result.converged) {
    s.jumpedT = undefined; // (a step taken from the jump: a later failure isn't its; Astra's review)
    s.readV = undefined; // (a step taken at the new drives: see Device's _redefineDrives)
    const PN = s.termCharge, PNLo = s.termChargeLo;
    // (and, for the error in each terminal's charge, see checkedCharges: where the history began
    // (o, t0), the charge checked since then (Y) and its gross (G), through the step's start)
    const o = prev ? prev.o : { P: PN, Lo: PNLo }, YN = checkedCharges(s, PN, PNLo, DN, QN, o);
    const GN = prev ? prev.G.map((g, k) => g + Math.abs(YN.c[k] - prev.Y.c[k])) : new Float64Array(PN.length);
    s.history.unshift({ t: s.time, dt, u: Float64Array.from(s.uPrev), V: Float64Array.from(s.termVPrev), c: cN, th: thN, th0: th0N, gam: gamN, g: gN, D: DN, Q: QN, P: PN, PLo: PNLo, Y: YN, G: GN, o, t0: prev ? prev.t0 : s.time });
    if (s.history.length > 3) s.history.length = 3;
    s.time = tEnd;
    s.lastDt = dtEff;
    s.contactDOld = s.contactDStart;
    s.portQOld = s.portQStart;
    s._accumulateBoundaryIntake(dtEff, cN, gamN);
    // The charge each terminal passed: its current, which is BDF2's rate of its stored charge at
    // the step's end, integrated by the same formula, P = P* + I dt/a0, so that a stored charge
    // passed is the stored charge's change exactly (and Kirchhoff's sum holds in it as in I).
    // Not Σ I·dt: that's a quadrature of the steps' end values, 22% short of an RC's charge at
    // the default tol, whose steps grow to several RC (the trace of a 5 mV step; Fable 5.2's
    // first read of the docs).
    // (In time, the state is no longer the steady solve's that a continuation would start from;
    // and stepped before the currents are read: a current-driven terminal's current is its
    // drive's from within the step, the before side of a jump at tEnd, as the step solved it.
    // Read after, a first fixed step onto a current's edge posted the far side's I·dt for good;
    // Sol's third review.)
    if (!s.steady) {
      s.solvedV = undefined;
      s.stepped = true;
    }
    // (Compensated, the running total's rounding kept in termChargeLo: after 1e12 C/m² had
    // passed, a current of 1 mA/m² for 1 s added nothing, its increments below an ulp; Sol's
    // re-read of the time integrator, round 7.)
    if (!s.steady && Number.isFinite(dtEff)) {
      const I = s._terminalCurrents(), P = new Float64Array(PN.length), Lo = new Float64Array(PN.length);
      for (let k = 0; k < PN.length; k++) {
        const q = PN[k], inc = (bdf ? b2 * (q - prev.P[k] + (PNLo[k] - prev.PLo[k])) : 0) + dtEff * I[k];
        const hi = q + inc, e = Math.abs(q) >= Math.abs(inc) ? q - hi + inc : inc - hi + q, lo = PNLo[k] + e;
        P[k] = hi + lo;
        Lo[k] = lo - (P[k] - hi);
      }
      s.termCharge = P;
      s.termChargeLo = Lo;
    }
  } else {
    s._restore(snap); // (and the concentrations from it)
  }
  s.sourceBefore = false;
  return result;
}

// Each terminal's charge since the history began (o: its P and PLo there), read two ways for its
// error (see chargeError): whole (c), and less what its contact or port stores (y: D, Q, at the
// step's start or now), which is the state's, under the potentials' control. Each reading is
// exact where the other carries a stored charge's round-off: a blocking cell's plates at rest
// rattle by 1e-13 C/m², which the whole carries; a strictly neutral contact's D, read from its
// box's balance, rattles ∝ c/dt, which only the second does. Either cut the steps without end
// once a terminal's charge was in the error, so the error is the smaller of the two. (A
// terminal driven by a current: its drive's integral, both ways.)
function checkedCharges(s, P, Lo, D, Q, o) {
  const c = Float64Array.from(P, (_, k) => P[k] - o.P[k] + (Lo[k] - o.Lo[k]));
  const y = Float64Array.from(c, (ck, k) => {
    if (s.terms[k].drive.kind === 'I') return ck;
    if (k < 2) return k === 0 ? ck - D.left : ck + D.right; // (the right contact's current is out of the device)
    return s.model.ports[k - 2].capacitance ? ck - Q[k - 2] : ck;
  });
  return { c, y };
}
// What a terminal's charge is read to in one step: 1e-12 of the charge its boxes hold (each
// species' |z|F c over the box beside a contact, or a port's window), which its currents are
// balances of, solved to Newton's precision. (A strictly neutral contact's box of 500 mol/m³
// over 10 µm, its current read from that balance, rattled by 5e-11 C/m² a step, 1e-13 of what it
// holds and half of what 1 mA/m² passed in 0.1 µs; a metal's box holds no concentration, and
// reads its current as a conductance's.)
function boxNoise(s, k) {
  const { n, z, c } = s, vol = s.model.grid.vol, nodes = k < 2 ? [k === 0 ? 0 : s.nNodes - 1] : s.model.ports[k - 2].nodes;
  let q = 0;
  for (const g of nodes) for (let i = 0; i < n; i++) if (Number.isFinite(c[g * n + i])) q += Math.abs(z[i] * c[g * n + i]) * vol[g];
  return 1e-12 * FARADAY * q;
}

const checkedNow = (s) => checkedCharges(s, s.termCharge, s.termChargeLo, s.contactDEnd, s.portQEnd, s.history[0].o);

// The first step's two quadratures of what each terminal passed (checkedCharges: Y0 at the
// start, Ymid after the first half step, YFull after the full one), full step against two halves,
// relative to the halves' (as a tol in thermal units bounds a potential's), or to what the
// current at the step's end would pass over the time T left in the drives' segment, if more.
// Every terminal's: on a flow-through path (a resistor, an electrode fed by a bath) nothing in
// the state feels the charge, and the pair, agreeing on the end state, spanned a slow CV's whole
// wave after a rest, 52% high at its vertex (Sol's re-read of the time integrator, round 7). (And
// the segment's share, since a current from zero, a ramp after a rest, errs by a third of its
// charge at any step: measured against the pair's alone, it was cut without end.) A difference
// below quiet is noise, see integrate.
function chargeError(s, Y0, Ymid, YFull, h, T, quiet) {
  const { c, y } = checkedNow(s);
  let e = 0;
  for (let k = 0; k < c.length; k++) {
    const d = Math.min(Math.abs(YFull.c[k] - c[k]), Math.abs(YFull.y[k] - y[k]));
    if (d > Math.max(quiet, boxNoise(s, k))) e = Math.max(e, d / Math.max(Math.abs(c[k] - Y0.c[k]), (2 * T * Math.abs(c[k] - Ymid.c[k])) / h, Number.MIN_VALUE));
  }
  return e;
}

// A later step's error in what each terminal passed (checkedCharges): against the predictor's
// (pred, through the history's), relative to its gross since the history began, this step's
// included (scaled after, as the potentials' error is). Not relative to the step's own: where a
// current crosses zero, that's as small as the error at any step, and the steps fell toward the
// crossing without end.
function chargeDrift(s, pred, quiet) {
  const e0 = s.history[0], { c, y } = checkedNow(s);
  let e = 0;
  for (let k = 0; k < c.length; k++) {
    const d = Math.min(Math.abs(c[k] - pred.c[k]), Math.abs(y[k] - pred.y[k]));
    if (d > Math.max(quiet, boxNoise(s, k))) e = Math.max(e, d / Math.max(e0.G[k] + Math.abs(c[k] - e0.Y.c[k]), Number.MIN_VALUE));
  }
  return e;
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

const absoluteOf = (s) => `tEnd is the time to reach in s, absolute (the device is at t = ${s.time}, which a solution gives as its time)`;

/**
 * What integrate refuses of its arguments alone, before anything is solved: Device.advance asks
 * first, since settling a pending regrid is a steady solve (a refused call settled it, and left
 * the device not as a twin that never made the call; round 7's S16). A NaN tEnd isn't refused:
 * the call does nothing, and says so (see integrate).
 */
export function integrateArgs(s, tEnd, opts = {}) {
  // tEnd is an absolute time. Not a number at all is an error; a NaN (a first animation frame's
  // clock, or a time read from a field that doesn't exist) does nothing, and says so.
  const absolute = absoluteOf(s);
  if (typeof tEnd !== 'number') throw new SolverError(`advance(tEnd): ${absolute}; got ${tEnd}`);
  // (advance(Infinity) set the time to ∞, converged in no steps, and nothing ran after it)
  if (tEnd === Infinity || tEnd === -Infinity) throw new SolverError(`advance(${tEnd}): ${absolute}. For the long-time limit, solve() the steady state`);
  if (Number.isNaN(tEnd)) return;
  const { tol = 1e-3, dtMax = Infinity, maxSteps = 100000, method } = opts;
  // (a negative dtMax ran the clock backwards and converged; a NaN tol stepped by NaN; Sol's
  // third review)
  if (!(tol > 0 && tol < Infinity)) throw new SolverError(`advance(): tol is the local error allowed per step, in thermal units, a positive number; got ${tol}`);
  if (!(dtMax > 0)) throw new SolverError(`advance(): dtMax is the longest step, in s, a positive number; got ${dtMax}`);
  if (!(Number.isInteger(maxSteps) && maxSteps >= 0)) throw new SolverError(`advance(): maxSteps is a whole number ≥ 0; got ${maxSteps}`);
  checkMethod('advance', method);
  s._probes(opts.probes); // (its readers built again by the call)
}

// A time step's method: 'bdf2' or 'be'. (Any other ran backward Euler without a word, step(dt,
// { method: 'rk4' }) among them; round 7's S16.)
export function checkMethod(where, method) {
  if (method !== undefined && method !== 'bdf2' && method !== 'be') throw new SolverError(`${where}(): method is 'bdf2' (variable-step BDF2) or 'be' (backward Euler); got ${typeof method === 'string' ? `'${method}'` : method}`);
}

/**
 * Adaptive time stepping to tEnd: variable-step BDF2 (or backward Euler), with the local error
 * estimated against an explicit predictor through the previous states, and controlled to
 * `tol` thermal units per step in every potential (φ̂, each η, a surface's coverages as η). The first
 * step is checked by step doubling. Stops early when `budgetMs` of wall time is used, so an
 * animation can call it once per frame; the step size carries over between calls.
 */
export function integrate(s, tEnd, opts = {}) {
  integrateArgs(s, tEnd, opts);
  const absolute = absoluteOf(s);
  if (Number.isNaN(tEnd)) return { converged: true, done: false, steps: 0, rejected: 0, iterations: 0, trace: { t: [], current: [], charge: [], voltage: [], terminals: Object.fromEntries(s.terms.map((t) => [t.name, { V: [], current: [], charge: [] }])) }, stopped: `advance(NaN) did nothing: ${absolute}` };
  const { tol = 1e-3, dtMax = Infinity, budgetMs = Infinity, maxSteps = 100000, method = 'bdf2' } = opts;
  const clock = () => (globalThis.performance ? globalThis.performance.now() : Date.now());
  const start = clock();
  const sample = s._probes(opts.probes);
  // (and every terminal's voltage, current and charge passed into the device, by name, as
  // sol.terminals, with the charge's running total: see step)
  const names = s.terms.map((t) => t.name), terminals = Object.fromEntries(names.map((nm) => [nm, { V: [], current: [], charge: [] }]));
  const trace = { t: [], current: [], charge: [], voltage: [], terminals, ...(sample ? { probes: sample.out } : {}) };
  const readTerminals = () => ({ V: Array.from(s.termV), current: s._terminalCurrents(), charge: Array.from(s.termCharge) });
  const pushTerminals = (r) => names.forEach((nm, k) => (terminals[nm].V.push(r.V[k]), terminals[nm].current.push(r.current[k]), terminals[nm].charge.push(r.charge[k])));
  let steps = 0, rejected = 0, iterations = 0, failed = false, stalled = null;
  s.modelError = null; // (a statistics model's refusal, for the failure: see newton)
  s.newtonError = null; // (and Newton's own reason)
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
  // A call that starts within round-off of a breakpoint it didn't land on (the last call ended a
  // hair before a pulse's edge, its target a sum of durations) steps across it from here, since
  // nextBreakpoint() counts it as passed: restart the order, as landing on it would have. BDF2's
  // history across a current pulse's edge lost 1/12 of the step's charge, for good (Fable round
  // 2's hunter: 8.33% of a 1 ms pulse, silent, the plates' charge outside the error norm).
  // (Restarted at the call's first attempt: a call that makes none, maxSteps 0 or a spent budget,
  // leaves the history and the step size carried as they were; cleared here, and dtNext set from
  // its own interval, a later advance(1e-3) took 51 steps where a twin took 32; round 7's S15.)
  let restart = s.history.length > 0 && s.terms.some((t) => nearBreakpoint(t.drive.src, s.time));
  // An explicit dt0 is this call's first step, over the size carried from the last call (which
  // can't see a set() that has just made the first instants fast; Fable round 2's hunter).
  let dt = opts.dt0 ?? s.dtNext ?? (tEnd - s.time) * 1e-4;
  if (!(dt > 0)) dt = (tEnd - s.time) * 1e-4;
  // (and at least what moves the clock: a dt0 of 1e-18 s at t = 1 s was refused as a step the
  // device needed, below the time's resolution; round 7)
  if (s.time + dt === s.time) dt = 4 * Number.EPSILON * Math.abs(s.time);
  // The shortest step worth trying: round-off relative to the time, or the device's fastest
  // time scale, whichever is shorter (a cold start can need steps far below a long run's 1e-14).
  const floor = Math.min(1e-14 * Math.max(tEnd, 1e-300), s.fastestTime() * 1e-2);
  const pred = new Float64Array(s.u.length), guess = new Float64Array(s.u.length), predV = new Float64Array(s.termV.length), predP = { c: new Float64Array(s.terms.length), y: new Float64Array(s.terms.length) };
  // The first step's error at each size tried (see below), kept on the solver across calls: a
  // first step cut across an animation's frames is seen to stall as in one call (it wasn't, and a
  // heavy device's first step took 134 one-attempt frames where one call took 30, against
  // docs/kit.md's "step for step"; Sol's third review). Cleared by an accepted step, and by
  // set() (a redrive, or a rebuilt solver). Each try's size with its error: a try of the same
  // size again isn't a cut (three calls of maxSteps 1, each from the same dt0, took a step at
  // 7e5 times tol, its trace 46% off; Sol's re-read of the time integrator, round 7).
  const firstErrs = (s.firstErrs ??= []);
  // What a terminal's charge is known to, below which a difference is noise: 1e-12 C/m², and what
  // a current of 1e-15 A/m² passes over the span ahead (below; as in solution.js, which reads a
  // sum of currents against 1e-15 A/m²), and its boxes' own (boxNoise). At rest, the currents a
  // held electrode reads (4e-16 A/m², jittering by a tenth from step to step), a floating gate's
  // contacts' 1e-19 to 1e-17 C/m² a step, and a bath's 1e-17 C/m² of displacement in its first
  // instants, each cut the steps without end once the charge was in the error.
  const area = Math.max(...s.model.grid.area), quiet = (T) => area * (1e-12 + 1e-15 * T);
  const factor = (err, p) => (err > 0 ? Math.min(2, Math.max(0.2, 0.9 * Math.exp(Math.log(tol / err) / (p + 1)))) : 2);
  while (s.time < tEnd) {
    if (steps + rejected >= maxSteps || clock() - start > budgetMs) break;
    // Within round-off of the end (a target a hair past a breakpoint just landed on, say):
    // snap to it rather than attempt a step of ~1e-13 s, which can't be resolved and, in a
    // strictly neutral device, fails outright. Round-off of the time and negligible beside the
    // interval asked for: by the time alone, a GITT pulse of 0.01 s at t = 1e9 s was skipped
    // whole (the time moved, nothing else did), and from 1e10 s 17% of a pulse train's charge.
    // Within 1e-12 of the time, a step's reach in a strictly neutral device: at 1e-10, the tail
    // after an edge (5e-6 s at 1e5 s) went with the current through it, for good (Sol's third
    // review); a hair's charge is what's lost now.
    if (tEnd - s.time <= 1e-12 * Math.abs(tEnd) && tEnd - s.time <= 1e-9 * (tEnd - t0)) {
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
    // (what's left of the drives' segment, or of the call where they have none ahead: the time a
    // terminal's charge is weighed over, see chargeError; later, from where the history began, which
    // doesn't vanish at the segment's end)
    const span = (Number.isFinite(tb) ? tb : tEnd) - s.time;
    // (the snapshot holds what a step sets of the steady solves' choices too: a rejected first
    // step left stepped set and solvedV cleared, a different continuation later; Fable round 4)
    if (restart) {
      s.history = [];
      restart = false;
    }
    const snap = s._snapshot();
    const back = () => s._restore(snap);
    let err, errQ = 0, errG = 0, p, half = null;
    if (s.history.length > 0) firstErrs.length = 0; // (a step taken since, by step() say)
    if (s.history.length === 0) {
      // First step: backward Euler, checked against two half steps (whose result is kept).
      const full = s.step(h);
      iterations += full.iterations;
      if (full.converged) {
        const uFull = Float64Array.from(s.u), vFull = Float64Array.from(s.termV), YFull = checkedNow(s);
        back();
        const a = s.step(h / 2);
        // The first half step is kept too, so it goes in the trace (a jump's charge is in it).
        half = a.converged ? { t: s.time, current: s._terminalCurrent(), voltage: s.termV[1] - s.termV[0], terminals: readTerminals(), probes: sample?.read() } : null;
        const b = a.converged ? s.step(h / 2) : a;
        iterations += a.iterations + (b === a ? 0 : b.iterations);
        if (b.converged) {
          // (and the charge each terminal passes, relative: nothing in the state need feel it on a
          // resistive path, where the pair's two quadratures of a current ramp agree at the step's
          // end, err 0, and a ramp's charge came out 50% high; Sol's third review. See chargeError.)
          err = errorNorm(s, uFull, vFull, false);
          errG = gateNorm(s, uFull);
          errQ = chargeError(s, s.history[1].Y, s.history[0].Y, YFull, h, span, quiet(span));
          p = 1;
        }
      }
      if (err === undefined) {
        back();
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
        back();
        rejected++;
        dt = h / 4;
        if (dt < floor) { failed = true; break; }
        continue;
      }
      p = predict(s, pred, r.bdf, predV, predP);
      err = errorNorm(s, pred, predV) * p.scale;
      errQ = p.charge ? chargeDrift(s, predP, quiet(span + snap.time - s.history[0].t0)) * p.scale : 0;
      p = p.order;
    }
    // A first step whose error doesn't shrink with it: a jump in a level held at a boundary
    // (a bath's, after a voltage step) starts a profile self-similar in x/√t, so a first step
    // of any length errs alike, down to the grid's own diffusion time, below Newton's reach.
    // Three tries in a row, each within a factor 2 of the last, say so: the error is the
    // jump's, not the step's. Take the step (backward Euler damps it as it should); the next is
    // sized from its error as any other. (Counted across calls: see firstErrs above.) Not a
    // gate's error, nor a charge's: a gate's is an ODE's at its own rate, which shrinks once the
    // step is below 1/(α + β), and three tries above that took a voltage clamp's first step after
    // an edge to −120 mV at 88 times tol, m 0.442 for 0.281, the Na⁺ tail 4× high for 40 µs (round
    // 7's state chaser); a charge's is the step's own, and shrinks as it's cut.
    if (half && err > tol) firstErrs.push({ h, err });
    const jump = half && firstErrs.length >= 3 && firstErrs.slice(-3).every((e, j, a) => j === 0 || (e.h < a[j - 1].h && e.err > 0.5 * a[j - 1].err && e.err < 2 * a[j - 1].err));
    if ((err > tol && !jump) || errG > tol || errQ > tol) {
      err = jump ? Math.max(errG, errQ) : Math.max(err, errG, errQ);
      back();
      rejected++;
      dt = h * factor(err, p);
      // Past the first step, a stall stops here, by its own signatures: where the error
      // estimate has a floor above tol (Newton's precision on the device), the step fell without
      // end, below the time's ulp, steps "accepted" that didn't move it, until maxSteps ran out
      // with currents of 1e87 beside converged: true (Fable round 4: three reviewers). So it
      // stops where a step can't move the clock. (Not at the non-convergence floor: at tight tol,
      // or late in a long run, steps far below it are legitimate; nor where the error stops
      // shrinking as the step is cut: right after a jump it doesn't, self-similar in x/√t, and
      // a run must go on through it; both from the re-checks.)
      if (!half && s.time + dt === s.time) {
        stalled = { err, h };
        failed = true;
        break;
      }
      continue;
    }
    steps++;
    err = Math.max(err, errG, errQ);
    firstErrs.length = 0;
    s.jumpedT = undefined; // (a step taken: see the first step's failure below)
    // The state the first step started from needn't satisfy the algebraic equations (φ in a
    // strictly neutral region, an interface's unknowns, after a start or a jump), and they
    // jump in the first instant; a predictor through it would see that jump as error on every
    // step after. The history starts from the half step.
    if (half) s.history.length = 1;
    dt = clamped ? Math.max(dt, h * factor(err, p)) : h * factor(err, p);
    // The solution's slope jumps at a breakpoint: restart the order there, landed on or ended
    // within round-off of (a step of regular size can end a few ulps short, outside atBreak's
    // 1e-9 of the step yet inside nextBreakpoint()'s 1e-12 of the time, which then counts the
    // breakpoint passed: BDF2 stepped across a pulse's edge, 1.4% of its charge lost; Fable
    // round 2's hunter, on its re-check).
    if (atBreak || s.terms.some((t) => nearBreakpoint(t.drive.src, s.time))) s.history = [];
    if (half) {
      trace.t.push(half.t);
      trace.current.push(half.current);
      trace.charge.push(-half.terminals.charge[1]);
      trace.voltage.push(half.voltage);
      pushTerminals(half.terminals);
      if (sample) sample.push(half.probes);
    }
    trace.t.push(s.time);
    trace.current.push(s._terminalCurrent());
    trace.charge.push(-s.termCharge[1]);
    trace.voltage.push(s.termV[1] - s.termV[0]);
    pushTerminals(readTerminals());
    if (sample) sample.push(sample.read());
  }
  if (steps + rejected > 0) s.dtNext = dt;
  s.landing = undefined;
  // Stopped short, other than on the caller's wall-time budget: say where, and why.
  let stopped;
  const at = `advance stopped at t = ${s.time.toPrecision(6)} s of ${tEnd.toPrecision(6)}`;
  // (What usually does it, most likely first. Round 11's agents read the old list, which led
  // with a drive's jump, as blaming a modest current pulse: a host emptied at its surface, and
  // an electrode at open circuit with no double layer, held only by a slow reaction.)
  const unbalanced = s.history.length === 0 && s.time === t0 ? (s.unbalanced ?? []).map((r) => s.model.regions[r].name ?? `regions[${r}]`) : [];
  if (failed && !stalled && dt < 8 * Number.EPSILON * Math.abs(s.time))
    // (advance(1e300) blamed the physics for a step of 9.6e-3 s at t = 3e296 s)
    stopped = `${at}: a step of ${dt.toExponential(1)} s, as short as the device needs, is below the time's own resolution at t = ${s.time.toPrecision(3)} s. Advance by intervals nearer the device's time scales, or solve() for the long-time limit.`;
  else if (failed && unbalanced.length > 0 && s.time === t0)
    stopped = `${at}: the start can't be neutral in ${unbalanced.join(', ')}, strictly neutral (ε = 0): the amounts its c0 gives leave a net charge that no species held at a contact's level can balance, and a step can't move a stored amount at once. Give c0 that balances (with any fixed charge), or leave out of c0 the species a contact feeds, so that they start at its level.`;
  else if (failed && violent.length > 0)
    stopped = `${at}: a step did not converge even when shortened to ${dt.toExponential(1)} s. The start is violent: ${violent.join('; ')}. A double layer starts uncharged, and what passes beside it closes such a gap in the first instant, faster than any step can follow. For a start in equilibrium, solve() first (the steady state at the starting drives), then change what the experiment changes, then advance(). Where solve() fails from there too, hold that terminal at a voltage, solve(), then let it float (set its I: 0): it keeps the state the held one gave it.`;
  else if (failed && s.modelError)
    stopped = `${at}: a step did not converge even when shortened to ${dt.toExponential(1)} s, the state it reached refused by a statistics model: ${s.modelError}. The model doesn't hold there: a smaller concentration or drive, or a model for it (lattice for crowding).`;
  else if (failed && s.jumpedT && steps === 0)
    // (Fable round 3's convergence hunt: a Si diode solved at 150 K, set() to 500 K, failed at its
    // first step at 1e-18 s, told to solve() from c0's layout. Its minority carriers, carried over at
    // e⁻⁶⁸ of the new equilibrium, must rise ~68 thermal units in the first instant.)
    stopped = `${at}: the first step did not converge even when shortened to ${dt.toExponential(1)} s, from a state carried across set()'s jump in temperature (${s.jumpedT[0]} → ${s.jumpedT[1]} K): the concentrations carry over, and what the new temperature makes of them in the first instant (a minority carrier far below its new equilibrium) is more than a step can follow. Change T in smaller jumps, advancing a little after each (a Si diode from 150 K to 500 K went in jumps of 50 K), or solve() at the new temperature and advance from there.`;
  else if (failed && fromLayout && steps === 0)
    // (round 15's fuzzer: the pn template advanced from its c0, a step of 3e-24 s failing, 23 of 62
    // plain diodes; nothing said to solve first. And a closed film of A = B with c0.A = 1e-30, where
    // a step of 1e-20 s makes 1e-20 of A, read only of a junction; round 4's chaser.)
    stopped = `${at}: the first step did not converge even when shortened to ${dt.toExponential(1)} s, from the starting state, which is c0's layout rather than a solution (a junction's depletion not yet formed, say, or a species started far below what a reaction or link makes of it in the first instant: a c0 of 1e-30 where it's made at 1 mol/(m³·s)). For a transient from equilibrium, solve() first at the starting drives, then change what the experiment changes, then advance(); for one from a species' absence, start it nearer what the shortest step makes of it (1e-20 of its scale, say).`;
  else if (stalled)
    stopped = `${at}: the error control cut the step to ${dt.toExponential(1)} s, below the time's own resolution at t = ${s.time.toPrecision(3)} s, without bringing the local error estimate (${stalled.err.toExponential(1)} in thermal units) under tol (${tol}). Early in a run that's a floor under the estimate, the precision the steps are solved to (about 1e-7 on steps far shorter than a device's dielectric relaxation, with ε > 0): loosen tol (${Math.min(1e-3, 10 * stalled.err).toExponential(0)} or more, say). Late in a long one, advance by intervals nearer the device's time scales.`;
  else if (failed) stopped = `${at}: a step did not converge even when shortened to ${dt.toExponential(1)} s${s.newtonError ? ` (Newton's last try: ${s.newtonError})` : ''}. Look for: a host or surface driven to its limit (emptied or filled at its surface by a current it can no longer pass); an electrode at open circuit held only by a slow reaction, with no double-layer capacitance to keep its level (give its face one: phi capacitive); a drive's jump too big for one step (ramp it); a driven current the device cannot keep passing (a diode's reverse recovery past its storage time: stop there); or a device the steady solves also find hard.`;
  else if (s.time < tEnd && steps + rejected >= maxSteps) stopped = `${at}: maxSteps (${maxSteps}) used up (${rejected} of them rejected). Call again to go on, or loosen tol.`;
  return { converged: !failed, done: s.time >= tEnd, steps, rejected, iterations, trace, ...(stopped ? { stopped } : {}) };
}

// Quadratic (or linear) extrapolation of the state a time h ahead, from the current state and
// the two most recent history entries: a starting guess for Newton. (In offsets built from the
// steps' own lengths, never differences of absolute times: a second into a run, after a jump
// that a fine grid resolves in steps of 1e-13 s, those differences keep only 3 digits.) Linear
// where the older step is below 1e-12 of the newer, which the quadratic's offsets can't tell
// apart: after a step of 5e-13 s and one of 1e5 s, −(1e5 + 5e-13) is −1e5, the weights were
// 0/0, and every later call failed "non-finite update" (Sol's re-read of the time integrator,
// round 7).
function extrapolate(s, out, h) {
  const [e0, e1] = s.history;
  const u = s.u;
  if (!e0) return undefined;
  if (!e1 || !(e1.dt >= 1e-12 * e0.dt)) {
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
// The terminals' charges too (predP), as stored quantities are (charge: false where there's no
// predictor for them, only the step's start).
// (Through a step below 1e-12 of the one after it, as the start's state alone, and the quadratic
// only where its three offsets are apart: see extrapolate.)
function predict(s, pred, bdf, predV, predP) {
  const [e0, e1, e2] = s.history;
  const h = e0.dt;
  if (!e1 || !(e1.dt >= 1e-12 * h)) {
    pred.set(e0.u);
    predV.set(e0.V);
    return { scale: 1, order: 1, charge: false };
  }
  const hp = e1.dt;
  if (bdf && e2 && e2.dt >= 1e-12 * (h + hp)) {
    const hpp = e2.dt;
    // Lagrange through the three states before at the new time (offsets as in extrapolate).
    const [l0, l1, l2] = lagrange3(0, -h, -(h + hp), -(h + hp + hpp));
    for (let k = 0; k < pred.length; k++) pred[k] = l0 * e0.u[k] + l1 * e1.u[k] + l2 * e2.u[k];
    for (let k = 0; k < predV.length; k++) predV[k] = l0 * e0.V[k] + l1 * e1.V[k] + l2 * e2.V[k];
    for (const w of ['c', 'y']) for (let k = 0; k < predP[w].length; k++) predP[w][k] = l0 * e0.Y[w][k] + l1 * e1.Y[w][k] + l2 * e2.Y[w][k];
    const w = h / hp;
    const Cc = (h * h * h * (1 + w) * (1 + w)) / (w * (1 + 2 * w));
    const Cp = h * (h + hp) * (h + hp + hpp);
    return { scale: Cc / (Cc + Cp), order: 2, charge: true };
  }
  const l0 = (h + hp) / hp, l1 = 1 - l0;
  for (let k = 0; k < pred.length; k++) pred[k] = l0 * e0.u[k] + l1 * e1.u[k];
  for (let k = 0; k < predV.length; k++) predV[k] = l0 * e0.V[k] + l1 * e1.V[k];
  for (const w of ['c', 'y']) for (let k = 0; k < predP[w].length; k++) predP[w][k] = l0 * e0.Y[w][k] + l1 * e1.Y[w][k];
  // Backward Euler's error against a linear predictor; after a BDF2 step without enough
  // history this overestimates, which is safe.
  return { scale: bdf ? 1 : h / (2 * h + hp), order: 1, charge: true };
}

// Largest difference from ref among the potentials that carry state, in thermal units: and the
// floating terminals' voltages (refV), which carry an RC's charging where φ is held beside it
// (a capacitive contact behind R against a bath: stepped blind, tol 1e-3 and 1e-5 alike took 3
// steps from τ to 3τ and erred 1% in the charge; Fable round 2's hunter).
function errorNorm(s, ref, refV, gates = true) {
  const { n, M, u } = s;
  let mx = 0;
  for (let g = 0; g < s.nNodes; g++) {
    const b = s.blockOfNode[g];
    if (!s.phiUndefined[g]) mx = Math.max(mx, Math.abs(u[b * M] - ref[b * M]));
    for (let i = 0; i < n; i++) {
      if (s.present[g * n + i]) mx = Math.max(mx, Math.abs(u[b * M + 1 + i] - ref[b * M + 1 + i]));
    }
    const k = s.surfAt[g]; // an electrode surface's coverages, as η
    if (k >= 0) for (let q = 0; q < s.model.surfaces[k].species.length; q++) mx = Math.max(mx, Math.abs(u[b * M + 1 + n + q] - ref[b * M + 1 + n + q]));
  }
  if (gates) mx = Math.max(mx, gateNorm(s, ref));
  for (const k of s.floating) mx = Math.max(mx, Math.abs(s.termV[k] - refV[k]) / s.VT);
  return mx;
}

// The gates' part of errorNorm (a fraction: 1e-3 of m³ is 0.3% of P).
function gateNorm(s, ref) {
  let mx = 0;
  for (const { o } of s.gateList) mx = Math.max(mx, Math.abs(s.u[o] - ref[o]));
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
