// Polarization curves at an electrode, for driftlet/kit: each reaction's rate with the solution
// beside the electrode held as it is and the metal's level moved, which is what an Evans diagram
// plots. Where the net current crosses zero is the spot's open-circuit (mixed) potential.

import { DeviceError, checkOptions } from './errors.js';
import { FARADAY } from './constants.js';
import { mhcRate } from './marcus.js';
import { powi, powr } from './pow.js';
import { solveDense } from './dense.js';
import { addTerm } from './sum.js';

const fail = (message) => {
  throw new DeviceError(message);
};

// The largest share of its exchange rate a step held at equilibrium may pass (see polarization).
const QUASI = 0.1;

/**
 * Each electrode reaction's rate at a face (a metal's) or at one spot of an electrode port, with
 * the composition beside it as in `sol` and the metal's carrier at each level in `V` (in the
 * terminals' convention, V = μ̄/(zF); against SHE, subtract φ for the usual tables' ions). The
 * rates are the solver's law, written again here from the normalised reaction (the tests hold the
 * two together, so at the electrode's actual level they are `sol.interfaces[f].rates`, or
 * `sol.ports[k].rates` at that node): a new kind of kinetics needs a branch here too, and one it
 * lacks is refused rather than read as Butler–Volmer.
 * @param {import('./index.js').Device} device
 * @param {object} sol a solution of it
 * @param {{ face: number } | { port: string | number, x?: number }} where a face's index, or a port
 *   (name or index) and the position in its window (m, the device's x; default its first node)
 * @param {number | ArrayLike<number>} V the metal's level(s), V
 * @param {{ surface?: 'held' | 'equilibrium' }} [opts] an electrode's surface coverages: held as
 *   in `sol` (the default: this instant's curves), or in equilibrium at each V (the steady-state
 *   curve, with a passive film's active–passive peak): each reaction with surface species at zero
 *   affinity, as far as they're independent (the fastest first, by exchange rate with the surface
 *   as held: where two would set one species, the faster does), and the totals they leave fixed
 *   kept as in `sol` (a tethered couple's). Those reactions are `equilibrated`, their rates what
 *   the surface's steady balance asks of them, given the others' (the others limit the rate: a
 *   quasi-equilibrium curve). One asked for more than a tenth of its exchange rate isn't near
 *   equilibrium, and the curve is refused at that V
 * @returns {{ V: Float64Array, current: Float64Array, bare?: Float64Array, coverage?: Record<string, Float64Array>, reactions: { equation: string, level: number,
 *   rate: Float64Array, current: Float64Array, equilibrated?: boolean }[] }} per reaction: its equilibrium level here (where
 *   its rate is zero; NaN for a chemical step, which takes no electron from the metal and passes no current), its rate (mol/(m²·s) of electrode, forward) and its current (A/m², anodic
 *   positive: into the solution); and their sum, with the surface's bare fraction and coverages
 *   at each V where there is one
 */
export function polarization(device, sol, where, V, opts = {}) {
  checkOptions('polarization', arguments[4], ['surface']);
  const model = device.model, { species, grid, RT } = model, solver = device.solver;
  // The solution's composition with the device's present kinetics: a solution from before a set()
  // that rebuilt the device is of another definition (k0 set 10× after it, the curve read 10× the
  // rates that composition had, and after a regrid its nodes were another grid's, NaN), so it's
  // refused. (A later solve or step is fine: a recorded run's frames each have their curve.) A copy
  // without its mark is held to the grid alone. (Astra's fifth pass)
  const state = device._stateOf?.(sol);
  const stale = state && !state.definition ? state.moved : sol?.x?.length !== grid.nNodes ? `it has ${sol?.x?.length} nodes, the device ${grid.nNodes}` : null;
  if (stale) fail(`polarization: sol isn't of the device's present definition (${stale}): solve() the device as it now stands and pass that solution`);
  const Vs = Float64Array.from(typeof V === 'number' ? [V] : V);
  let reactions, nodeOf, carrier, names, surface = null; // (an electrode port's coverages at the spot)
  let atPlane = () => null;
  if (where?.face !== undefined) {
    const f = where.face, itf = model.interfaces[f];
    if (!itf) fail(`polarization: no face ${f} (the device has ${model.interfaces.length})`);
    const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
    const metal = [gL, gR].findIndex((g) => solver.nodeConductor[g] >= 0);
    if (metal < 0) fail(`polarization: face ${f} has no metal on either side`);
    carrier = { side: metal, i: solver.nodeConductor[metal ? gR : gL] };
    nodeOf = (p) => (p.side ? gR : gL);
    // (beside a diffuse layer, a participant reacts at its plane: its concentration and level there)
    const layers = sol.interfaces[f].diffuse ?? {};
    atPlane = (p, name) => layers[p.side ? 'right' : 'left'] && { c: layers[p.side ? 'right' : 'left'].plane[name], dmu: layers[p.side ? 'right' : 'left'].draw?.[name] ?? 0 };
    reactions = itf.reactions;
    names = device.def.interfaces?.[f]?.reactions;
    // (its surface, an electrode's, as a port's: each species' coverage, and the bare fraction)
    if (itf.surface.length > 0) surface = { sp: itf.surface, theta: itf.surface.map((sp) => sol.interfaces[f].coverage[sp.name]), bare: sol.interfaces[f].bare };
  } else if (where?.port !== undefined) {
    const k = typeof where.port === 'number' ? where.port : model.ports.findIndex((p) => p.name === where.port);
    const port = model.ports[k];
    if (!port) fail(`polarization: no port ${JSON.stringify(where.port)}`);
    if (port.reactions.length === 0) fail(`polarization: port ${port.name} has no reactions`);
    let g = port.nodes[0];
    if (where.x !== undefined) for (const h of port.nodes) if (Math.abs(grid.x[h] - where.x) < Math.abs(grid.x[g] - where.x)) g = h;
    carrier = { side: 1, i: port.terminal };
    nodeOf = () => g;
    if (port.surface.length > 0) {
      const w = port.nodes.indexOf(g), theta = port.surface.map((sp) => sol.ports[k].coverage[sp.name][w]);
      surface = { sp: port.surface, theta, bare: sol.ports[k].bare[w] };
    }
    reactions = port.reactions;
    names = device.def.ports?.[k]?.reactions;
  } else fail('polarization: say where, { face: f } or { port: name, x }');
  if (reactions.length === 0) fail('polarization: no reactions there');

  // Each reaction as a = base + s·V − Σ ν_s ln a_s (the surface species' activities apart), with
  // its prefactor apart from them too, so a surface can be held or re-equilibrated at each V.
  // (The affinity summed as the solver sums it: compensated, base + baseLo; of the device's present
  // state, each level as the solver holds it, η = u + uLo, not μ̄/RT from the solution's doubles;
  // the carrier's as z·V/V_T; and the Butler–Volmer factor by expm1. Summed plainly from the
  // solution's μ̄, an iron port's fast couple, k0 1e4, read 0.56 off the solved rate and at 1e5
  // read 0, its rate a remainder near equilibrium below the doubles' rounding; now it's the solved
  // rate to the last bit at every k0. Astra's fifth pass.)
  const acc = new Float64Array(2), exact = state?.moved === null;
  const parts = reactions.map((rx, k) => {
    if (rx.srh) fail('polarization: SRH kinetics has no electrode level');
    if (!rx.vmax && !rx.mhc && !(rx.alpha >= 0)) fail('polarization: a kind of kinetics it has no branch for (Butler–Volmer, saturating and Marcus–Hush–Chidsey it has)');
    let s = 0, q = 0, pref = rx.vmax ? rx.vmax : rx.k0, lnP = 0, lnQ = 0; // (Marcus–Hush–Chidsey's P and Q apart)
    acc[0] = rx.fixedA;
    acc[1] = 0;
    const surf = [], carriers = [];
    // (a face's surface species are rx.surf, a port's side-2 participants: the same from here on)
    const all = [...rx.part, ...(rx.surf ?? []).map((p) => ({ ...p, side: 2 }))];
    all.forEach((p, x) => {
      if (p.side === 2) {
        // a surface species: μ = μ° + RT ln a, a = θ/θ₀
        addTerm(acc, -p.nu, surface.sp[p.s].mu0 / RT, 0);
        surf.push({ s: p.s, nu: p.nu, e: p.nu < 0 ? -p.nu * (1 - rx.alpha) : p.nu * rx.alpha });
        return;
      }
      if (p.side === carrier.side && p.i === carrier.i) {
        s -= (p.nu * species[p.i].z * FARADAY) / RT;
        carriers.push({ nu: p.nu, z: species[p.i].z });
        return;
      }
      const g = nodeOf(p), name = species[p.i].name, pl = atPlane(p, name);
      if (exact) {
        const o = solver.blockOfNode[g] * solver.M + 1 + p.i;
        addTerm(acc, -p.nu, solver.u[o], solver.uLo[o]);
      } else addTerm(acc, -p.nu, sol.mu[name][g] / RT, 0);
      if (pl) addTerm(acc, -p.nu, pl.dmu / RT, 0);
      q += p.side === carrier.side ? 0 : p.nu * species[p.i].z; // charge the forward reaction puts into the solution
      // (Another conductor's carrier, the far metal's electrons at a metal–metal face, has activity 1,
      // as in the solver: read from its c, NaN in a metal, every rate was NaN; Astra's fifth pass.)
      if (solver.nodeConductor[g] === p.i && !rx.vmax) return;
      const c = pl ? pl.c : sol.c[name][g];
      if (rx.vmax) {
        if (p.nu < 0) pref *= powi(c / (c + rx.K[x]), -p.nu);
      } else if (rx.mhc) {
        const l = Math.log(c / solver.cRef[g * species.length + p.i]);
        if (p.nu < 0) lnP -= p.nu * l;
        else lnQ += p.nu * l;
      } else pref *= powr(c / solver.cRef[g * species.length + p.i], p.nu < 0 ? -p.nu * (1 - rx.alpha) : p.nu * rx.alpha);
    });
    // (A chemical step, s = 0, 2 H = H₂ beside H⁺ + e⁻ = H: no current and no level, its rate at the
    // composition held, or with the coverages re-equilibrated at each V. Once refused.)
    const base = acc[0] + acc[1], baseLo = acc[1] - (base - acc[0]);
    return { rx, k, base, baseLo, s, carriers, q, pref, surf, lnP, lnQ };
  });
  // The surface's activities (ln a) and bare fraction at a level: held as `sol` has them, or in
  // equilibrium there.
  const held = surface && surface.theta.map((t) => Math.log(t / surface.bare));
  if (opts.surface !== undefined && opts.surface !== 'held' && opts.surface !== 'equilibrium') fail("polarization: surface must be 'held' (the default) or 'equilibrium'");
  // A reaction's rate at level v with the surface's activities lnA and bare fraction, or (exchange)
  // its exchange rate there, ∂rate/∂a at zero affinity: the one-way rate at equilibrium, how fast
  // the step can pass what's asked of it near equilibrium.
  const law = ({ rx, base, baseLo, carriers, pref, surf, lnP, lnQ }, v, lnA, bare, exchange = false) => {
    let p = pref * (rx.sites ? powr(bare, rx.sites) : 1);
    acc[0] = base;
    acc[1] = baseLo;
    for (const cr of carriers) addTerm(acc, -cr.nu, (cr.z * v) / solver.VT, 0); // (as the solver's portEta, its rounding too)
    for (const u of surf) {
      addTerm(acc, -u.nu, lnA[u.s], 0);
      p *= Math.exp(u.e * lnA[u.s]) * powr(bare, u.e); // its coverage θ = θ₀ a
    }
    const a = acc[0] + acc[1];
    if (rx.mhc) {
      // (a surface species at θ = θ₀ a, and the bare sites the transition state holds: θ₀^m forward, θ₀^(m−Δ) back)
      let sP = lnP, sQ = lnQ;
      for (const u of surf) (u.nu < 0 ? (sP -= u.nu * (lnA[u.s] + Math.log(bare))) : (sQ += u.nu * (lnA[u.s] + Math.log(bare))));
      if (rx.m !== undefined) {
        // (a face's reaction without surface species has no m; a port's always has, bare: true included)
        sP += rx.m * Math.log(bare);
        sQ += (rx.m - rx.fills) * Math.log(bare);
      }
      return exchange ? mhcRate(rx.mhc, rx.k0, sP, sQ, 0).ra : mhcRate(rx.mhc, rx.k0, sP, sQ, a).rate;
    }
    if (exchange) return p; // (Butler–Volmer's and the saturating law's slope at a = 0)
    if (rx.vmax) return p * -Math.expm1(-a);
    // (e^{αa} − e^{−(1−α)a} as e^{−(1−α)a}·expm1(a), as the solver's bvFactor; far out, the plain difference)
    const e = Math.exp(-(1 - rx.alpha) * a);
    return p * (Math.abs(a) > 50 ? Math.exp(rx.alpha * a) - e : e * Math.expm1(a));
  };
  const equilibrium = surface && opts.surface === 'equilibrium' ? surfaceEquilibrium(parts, surface, held, (pt) => law(pt, 0, held, surface.bare, true)) : null;
  const surfaceAt = (v) => {
    if (!surface) return { lnA: [], bare: 1 };
    return equilibrium ? equilibrium(v) : { lnA: held, bare: surface.bare };
  };
  const at = Array.from(Vs, (v) => surfaceAt(v));
  const out = parts.map((pt) => {
    const { k, base, s, surf, q } = pt;
    const rate = Float64Array.from(Vs, (v, j) => law(pt, v, at[j].lnA, at[j].bare));
    // (its level with the surface as held)
    const aHeld = surf.reduce((t, u) => t - u.nu * (held ? held[u.s] : 0), base);
    return { equation: names?.[k]?.equation ?? `reaction ${k}`, level: s === 0 ? NaN : -aHeld / s, rate, q };
  });
  // A step held at equilibrium isn't idle: at the steady state its rate is what the surface's
  // balance asks of it, every surface species made as fast as it's used, Σ ν r = 0, with the other
  // steps' rates as their kinetics give them. Read off its own law (zero affinity), it was 0: with
  // Volmer and Tafel, Volmer listed first, the curve was 0 where the solve passed 0.192 A/m², and
  // the equilibrated step's current is the whole current (Astra's fifth pass). The steps the
  // equilibrium rests on are independent rows, so their rates follow from the rows' Gram matrix,
  // as the coverages do; what the others make or use has no part along the free directions (each
  // reaction's ν is orthogonal to them), so it balances exactly.
  // That's a quasi-equilibrium only while the step is fast enough for its share: passing y near
  // equilibrium takes it an affinity of about y/i₀ (i₀ its exchange rate), and the activities it
  // sets are about that fraction off the surface's steady state. Past a tenth it's refused, at the
  // first V (as given) where it is: the step isn't fast enough there for what it passes.
  // (Volmer and Tafel at k0 1e-3 each, round 5's A8 device, read 2.6 times the solved current at
  // −0.1 V, unwarned; with Volmer at 1e2 against Tafel's 1e-3 it holds to −0.4 V, θ₀ 5% off.)
  if (equilibrium && equilibrium.rows.length > 0) {
    const { rows, gram } = equilibrium, X = surface.sp.length, equilibrated = rows.map((r) => r.pt.k);
    for (let j = 0; j < Vs.length; j++) {
      const made = new Float64Array(X); // what the other steps make of each surface species
      parts.forEach((pt, k) => {
        if (!equilibrated.includes(k)) for (const u of pt.surf) made[u.s] += u.nu * out[k].rate[j];
      });
      const y = solveDense(gram, rows.map(({ r }) => -r.reduce((t, x, s) => t + x * made[s], 0)));
      rows.forEach(({ pt }, m) => {
        out[pt.k].rate[j] = y[m];
        const i0 = law(pt, Vs[j], at[j].lnA, at[j].bare, true);
        if (!(Math.abs(y[m]) <= QUASI * i0)) {
          fail(
            `polarization: at V = ${Vs[j]} the surface isn't near equilibrium: ${out[pt.k].equation}, held at equilibrium as the fastest step (by exchange rate, ` +
              `the surface as sol has it), would pass ${y[m].toPrecision(3)} mol/(m²·s), ${(Math.abs(y[m]) / i0).toPrecision(3)} times its exchange rate there ` +
              `(${i0.toPrecision(3)}), and a quasi-equilibrium curve needs under ${QUASI}. The surface's steady state there is a solve's (the electrode set to that level); or hold the surface (surface: 'held')`,
          );
        }
      });
    }
    for (const k of equilibrated) out[k].equilibrated = true;
  }
  for (const r of out) {
    r.current = r.rate.map((x) => r.q * FARADAY * x);
    delete r.q;
  }
  const current = new Float64Array(Vs.length);
  for (const r of out) for (let j = 0; j < Vs.length; j++) current[j] += r.current[j];
  const coverage = surface && Object.fromEntries(surface.sp.map((sp, x) => [sp.name, Float64Array.from(at, ({ lnA, bare }) => bare * Math.exp(lnA[x]))]));
  return { V: Vs, current, ...(surface ? { bare: Float64Array.from(at, (x) => x.bare), coverage } : {}), reactions: out };
}

// The surface in equilibrium at each level: each reaction with surface species at zero affinity,
// ν·ln a = base + s·V (ln a = ln θ/θ₀, θ₀ the bare fraction), as far as they're independent, the
// fastest first (by exchange(pt), its exchange rate with the surface as held: where two would set
// one species, the one that exchanges it faster holds it at equilibrium, and the other limits),
// and what they leave free fixed by the totals they keep: a tethered couple's θ_R + θ_O, a species
// no reaction touches. (Once in listing order, single-species steps first: with Tafel's 2 H = H₂,
// k0 1e-3, listed before Volmer, k0 1e2, Tafel was held at equilibrium, θ_H stuck at ½, and the
// curve read 1.2e5 times the solved current, unwarned; Neon's fourth review, R4 3.2.) Those free
// directions, n, are the ones every reaction's ν is orthogonal to, so n·θ is what the surface
// conserves, and with ln a = P(V) + Σ t_k n_k the totals are kept where t minimizes
// F(t) = ln(1 + Σ e^{ln a}) − t·(n·θ_held): strictly convex while any site is bare, its gradient
// n·(θ − θ_held) and its Hessian the covariance of n over the sites (the bare one at n = 0), so
// Newton's method finds it. Done
// plainly, as n·θ − n·θ_held and Σ n n θ − m mᵀ, both are differences of numbers near 1 apart by
// about the bare fraction, and lose −log₁₀ θ₀ digits (2.8% of the Hessian at θ₀ = 1e-15, by
// Neon's count, all of it below 1e-16). A tethered couple 99.99% covered and swept 0.6 V reached
// its root, but there its Newton step was the gradient's rounding over a Hessian of 5e-5,
// 1.1e-11, above the 1e-12 it stopped at, and F's value, searched along, couldn't tell a move
// from rounding: 100 iterations and a throw (Neon's R4 3.1). So here each is summed without the
// cancellation, the step is capped and searched along by F's slope, and the free directions are
// exact (see below).
function surfaceEquilibrium(parts, surface, held, exchange) {
  const X = surface.sp.length;
  const dot = (a, b) => a.reduce((t, x, j) => t + x * b[j], 0);
  // Independent rows by Gram–Schmidt (q: the orthonormal basis so far).
  const q = [], rows = [];
  const orth = (v) => {
    const w = Float64Array.from(v);
    for (const b of q) {
      const c = dot(w, b);
      for (let j = 0; j < X; j++) w[j] -= c * b[j];
    }
    return w;
  };
  // (fastest first; an exchange rate that isn't a positive number, as 0; ties in listing order)
  const candidates = parts
    .filter((pt) => pt.surf.length > 0)
    .map((pt) => ({ pt, speed: exchange(pt) > 0 ? exchange(pt) : 0 }))
    .sort((a, b) => b.speed - a.speed || a.pt.k - b.pt.k);
  for (const { pt } of candidates) {
    const r = new Float64Array(X);
    for (const u of pt.surf) r[u.s] += u.nu;
    const w = orth(r), n = Math.sqrt(dot(w, w));
    if (n <= 1e-9 * Math.sqrt(dot(r, r))) continue;
    q.push(w.map((x) => x / n));
    rows.push({ r, pt });
  }
  // The free directions, the rows' null space, from their reduced echelon form: for whole (or
  // halved) ν it's exact, R n = 0 to the bit, so n·θ is kept exactly. (As Gram–Schmidt's
  // complement, a tethered couple's n_R and n_O were an ulp apart, and 99.9999999999999% covered,
  // its bare fraction moved by half as V took θ_R from 0.002 to 1.)
  const E = rows.map(({ r }) => Array.from(r)), pivots = [];
  const big = Math.max(0, ...E.flat().map(Math.abs));
  for (let c = 0, i = 0; c < X && i < E.length; c++) {
    let p = i;
    for (let j = i + 1; j < E.length; j++) if (Math.abs(E[j][c]) > Math.abs(E[p][c])) p = j;
    if (Math.abs(E[p][c]) <= 1e-9 * big) continue;
    [E[i], E[p]] = [E[p], E[i]];
    const lead = E[i][c];
    E[i] = E[i].map((x) => x / lead);
    for (let j = 0; j < E.length; j++) if (j !== i && E[j][c] !== 0) E[j] = E[j].map((x, y) => x - E[j][c] * E[i][y]);
    pivots.push(c);
    i++;
  }
  const free = [];
  for (let f = 0; f < X; f++) {
    if (pivots.includes(f)) continue;
    const n = new Float64Array(X);
    n[f] = 1;
    pivots.forEach((c, i) => (n[c] = -E[i][f]));
    free.push(n);
  }
  // Where every reaction keeps the number of sites filled, the filled total, n = (1, …, 1), is free
  // (the sum of those above: each has 1 at its own free column), and it goes first, in the last
  // one's place: a nearly full surface is stiff along it by only its bare fraction, and as a
  // direction of its own, its part of F's gradient and Hessian is the bare sites' alone. (In the
  // echelon basis two tethered couples' totals each held part of it; at a bare fraction of 8e-16,
  // the Newton step along it was the other directions' rounding, 2e-8 in ln a, and never shrank.)
  if (free.length > 1 && rows.every(({ r }) => r.reduce((s, x) => s + x, 0) === 0)) {
    free.pop();
    free.unshift(new Float64Array(X).fill(1));
  }
  const k = rows.length, d = free.length, S = X + 1;
  // (the rows' Gram matrix, for the particular solution P = Rᵀ (R Rᵀ)⁻¹ b, the one with no free part;
  // the rows and it are the equilibrated steps' rates too: see polarization)
  const gram = rows.map((a) => rows.map((b) => dot(a.r, b.r)));
  // Each free direction over the sites, the bare one last at 0; the held occupancies, the bare
  // fraction as the solution has it (not 1 − Σθ); and the site the held surface fills most.
  const nOf = free.map((n) => [...n, 0]);
  const heldTheta = [...surface.theta, surface.bare];
  const top = heldTheta.reduce((m, x, s) => (x > heldTheta[m] ? s : m), 0);
  // (t starts as the held surface's, least squares: N t nearest its ln a)
  const lnHeld = held.map((l) => (Number.isFinite(l) ? l : 0));
  let t = d > 0 ? solveDense(free.map((a) => free.map((b) => dot(a, b))), free.map((n) => dot(n, lnHeld))) : [];
  const lnAOf = (P, t) => P.map((p, x) => free.reduce((s, n, j) => s + t[j] * n[x], p));
  const lse = (lnA) => {
    const m = Math.max(0, ...lnA);
    return m + Math.log(lnA.reduce((s, l) => s + Math.exp(l - m), Math.exp(-m)));
  };
  // (the sites' occupancies at t, the bare one last; F's gradient n·(θ − θ_held) with each site's
  // n taken from the top one's, whose own difference, of two numbers near 1, drops out, the
  // occupancies summing to 1 either way)
  const occupancy = (P, t) => {
    const lnA = lnAOf(P, t), Z = lse(lnA);
    return [...lnA.map((l) => Math.exp(l - Z)), Math.exp(-Z)];
  };
  const gradient = (th) => nOf.map((n) => th.reduce((g, x, s) => (s === top ? g : g + (n[s] - n[top]) * (x - heldTheta[s])), 0));
  const at = (v) => {
    const P = new Float64Array(X);
    if (k > 0) {
      const y = solveDense(gram, rows.map(({ pt }) => pt.base + pt.s * v));
      rows.forEach(({ r }, j) => {
        for (let x = 0; x < X; x++) P[x] += y[j] * r[x];
      });
    }
    // (F's slope along dir at t + τ·dir)
    const slope = (tau, dir) => dot(gradient(occupancy(P, t.map((x, j) => x + tau * dir[j]))), dir);
    let cap = 1, best = Infinity, stalled = 0;
    for (let it = 0; d > 0; it++) {
      if (it === 200) fail(`polarization: the surface's equilibrium at V = ${v} didn't converge`);
      const th = occupancy(P, t), g = gradient(th);
      // Newton's step, H z = −g, from H's square root: H = AᵀA, A's rows √(θ_s θ_u)(n_s − n_u)
      // over the pairs of sites, and A = QR by Householder (its rows largest first), so RᵀR z = −g.
      // (Two couples tethered at a bare fraction of 4e-9: H's entries were 0.24 and its least
      // eigenvalue, along their common total, 4e-9; eliminating H subtracts the one from the other,
      // and swept 1.4 V, through a bare fraction of 1e-21, its pivot came out 0.)
      const A = [];
      for (let s = 0; s < S; s++) {
        for (let u = s + 1; u < S; u++) {
          const w = Math.sqrt(th[s]) * Math.sqrt(th[u]), row = nOf.map((n) => w * (n[s] - n[u]));
          if (row.some((x) => x !== 0)) A.push(row);
        }
      }
      A.sort((a, b) => dot(b, b) - dot(a, a));
      const R = householderR(A, d), step = new Array(d);
      for (let i = 0; i < d; i++) step[i] = (-g[i] - R.slice(0, i).reduce((s, row, j) => s + row[i] * step[j], 0)) / R[i][i];
      for (let i = d - 1; i >= 0; i--) step[i] = (step[i] - R[i].slice(i + 1).reduce((s, x, j) => s + x * step[i + 1 + j], 0)) / R[i][i];
      // (its size in ln a: the largest change it makes to a species' ln a)
      const size = Math.max(...free[0].map((_, x) => Math.abs(free.reduce((s, n, j) => s + step[j] * n[x], 0))));
      if (!Number.isFinite(size)) fail(`polarization: the surface's equilibrium at V = ${v} didn't converge (its Hessian is singular: bare fraction ${th[X]})`);
      if (size < 1e-12) break; // (near quadratic here)
      // (at the doubles' rounding floor the step stops shrinking: what's left is theirs)
      if (size < 0.5 * best) (best = size), (stalled = 0);
      else if (size < 1e-8 && ++stalled >= 3) break;
      // A short step is taken whole: F is near its quadratic there, and along a step that small its
      // slope can be the rounding of the stiff directions, not the soft one the step is in (two
      // couples' totals at a bare fraction of 5e-15: the search stalled at a step of 1e-8).
      if (size <= 0.25) {
        t = t.map((x, j) => x + step[j]);
        continue;
      }
      // Longer, the step capped, and along it F's slope rising from below zero: the whole (capped)
      // step where F still falls at its end, the cap then growing, else a point short of where the
      // slope turns, by regula falsi (Illinois), as src/carry.js carries a surface's amounts.
      const cut = Math.min(1, cap / size), dir = step.map((x) => x * cut);
      const g0 = dot(g, dir);
      if (!(g0 < 0)) fail(`polarization: the surface's equilibrium at V = ${v} didn't converge (Newton's step doesn't descend)`);
      let lo = 1, ghi = slope(1, dir);
      if (ghi > 0) {
        let glo = g0, hi = 1, side = 0;
        lo = 0;
        for (let ls = 0; ls < 60 && hi - lo > 1e-12; ls++) {
          const tau = (lo * ghi - hi * glo) / (ghi - glo), gt = slope(tau, dir);
          if (gt <= 0) {
            (lo = tau), (glo = gt);
            if (side === -1) ghi *= 0.5;
            side = -1;
            if (-gt <= -0.5 * g0) break;
          } else {
            (hi = tau), (ghi = gt);
            if (side === 1) glo *= 0.5;
            side = 1;
          }
        }
        if (cut < 1) cap = Math.max(1, cap / 2);
      } else if (cut < 1) cap *= 4;
      if (!(lo > 0)) fail(`polarization: the surface's equilibrium at V = ${v} didn't converge (no point along Newton's step lowers it)`);
      t = t.map((x, j) => x + lo * dir[j]);
    }
    const lnA = lnAOf(P, t);
    return { lnA, bare: Math.exp(-lse(lnA)) };
  };
  return Object.assign(at, { rows, gram });
}

// The d × d triangle R of A = QR (A's m rows of d, copied), by Householder reflections: AᵀA = RᵀR,
// with R's least diagonal as accurate as A's rows hold it (no AᵀA formed). Fewer rows than d, or
// a column of zeros, leaves a zero on R's diagonal.
function householderR(A, d) {
  const M = A.map((row) => Array.from(row)), m = M.length;
  const R = Array.from({ length: d }, () => new Array(d).fill(0));
  for (let c = 0; c < d; c++) {
    if (c >= m) break;
    const scale = M.slice(c).reduce((s, row) => Math.max(s, Math.abs(row[c])), 0);
    if (scale > 0) {
      let norm = 0;
      for (let i = c; i < m; i++) norm += (M[i][c] / scale) * (M[i][c] / scale);
      const alpha = (M[c][c] > 0 ? -1 : 1) * scale * Math.sqrt(norm);
      const v = M.slice(c).map((row) => row[c]);
      v[0] -= alpha;
      const vv = v.reduce((s, x) => s + x * x, 0);
      if (vv > 0) {
        for (let j = c; j < d; j++) {
          const f = (2 * v.reduce((s, x, i) => s + x * M[c + i][j], 0)) / vv;
          for (let i = c; i < m; i++) M[i][j] -= f * v[i - c];
        }
      }
      M[c][c] = alpha;
    }
    for (let j = c; j < d; j++) R[c][j] = M[c][j];
  }
  return R;
}
