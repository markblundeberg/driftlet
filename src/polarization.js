// Polarization curves at an electrode, for driftlet/kit: each reaction's rate with the solution
// beside the electrode held as it is and the metal's level moved, which is what an Evans diagram
// plots. Where the net current crosses zero is the spot's open-circuit (mixed) potential.

import { DeviceError, checkOptions } from './errors.js';
import { FARADAY } from './constants.js';
import { mhcRate } from './marcus.js';
import { powi, powr } from './pow.js';
import { solveDense } from './dense.js';

const fail = (message) => {
  throw new DeviceError(message);
};

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
 *   affinity, as far as they're independent (those making one from the solution alone first), and
 *   the totals they leave fixed kept as in `sol` (a tethered couple's)
 * @returns {{ V: Float64Array, current: Float64Array, bare?: Float64Array, coverage?: Record<string, Float64Array>, reactions: { equation: string, level: number,
 *   rate: Float64Array, current: Float64Array }[] }} per reaction: its equilibrium level here (where
 *   its rate is zero; NaN for a chemical step, which takes no electron from the metal and passes no current), its rate (mol/(m²·s) of electrode, forward) and its current (A/m², anodic
 *   positive: into the solution); and their sum, with the surface's bare fraction and coverages
 *   at each V where there is one
 */
export function polarization(device, sol, where, V, opts = {}) {
  checkOptions('polarization', arguments[4], ['surface']);
  const model = device.model, { species, grid, RT } = model, solver = device.solver;
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
  const parts = reactions.map((rx, k) => {
    if (rx.srh) fail('polarization: SRH kinetics has no electrode level');
    if (!rx.vmax && !rx.mhc && !(rx.alpha >= 0)) fail('polarization: a kind of kinetics it has no branch for (Butler–Volmer, saturating and Marcus–Hush–Chidsey it has)');
    let base = rx.fixedA, s = 0, q = 0, pref = rx.vmax ? rx.vmax : rx.k0, lnP = 0, lnQ = 0; // (Marcus–Hush–Chidsey's P and Q apart)
    const surf = [];
    // (a face's surface species are rx.surf, a port's side-2 participants: the same from here on)
    const all = [...rx.part, ...(rx.surf ?? []).map((p) => ({ ...p, side: 2 }))];
    all.forEach((p, x) => {
      if (p.side === 2) {
        // a surface species: μ = μ° + RT ln a, a = θ/θ₀
        base -= (p.nu * surface.sp[p.s].mu0) / RT;
        surf.push({ s: p.s, nu: p.nu, e: p.nu < 0 ? -p.nu * (1 - rx.alpha) : p.nu * rx.alpha });
        return;
      }
      if (p.side === carrier.side && p.i === carrier.i) {
        s -= (p.nu * species[p.i].z * FARADAY) / RT;
        return;
      }
      const g = nodeOf(p), name = species[p.i].name, pl = atPlane(p, name);
      base -= (p.nu * (sol.mu[name][g] + (pl ? pl.dmu : 0))) / RT;
      q += p.side === carrier.side ? 0 : p.nu * species[p.i].z; // charge the forward reaction puts into the solution
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
    return { rx, k, base, s, q, pref, surf, lnP, lnQ };
  });
  // The surface's activities (ln a) and bare fraction at a level: held as `sol` has them, or in
  // equilibrium there.
  const held = surface && surface.theta.map((t) => Math.log(t / surface.bare));
  if (opts.surface !== undefined && opts.surface !== 'held' && opts.surface !== 'equilibrium') fail("polarization: surface must be 'held' (the default) or 'equilibrium'");
  const equilibrium = surface && opts.surface === 'equilibrium' ? surfaceEquilibrium(parts, surface, held) : null;
  const surfaceAt = (v) => {
    if (!surface) return { lnA: [], bare: 1 };
    return equilibrium ? equilibrium(v) : { lnA: held, bare: surface.bare };
  };
  const at = Array.from(Vs, (v) => surfaceAt(v));
  const out = parts.map(({ rx, k, base, s, q, pref, surf, lnP, lnQ }) => {
    const rate = Float64Array.from(Vs, (v, j) => {
      const { lnA, bare } = at[j];
      let a = base + s * v, p = pref * (rx.sites ? powr(bare, rx.sites) : 1);
      for (const u of surf) {
        a -= u.nu * lnA[u.s];
        p *= Math.exp(u.e * lnA[u.s]) * powr(bare, u.e); // its coverage θ = θ₀ a
      }
      if (rx.mhc) {
        // (a surface species at θ = θ₀ a, and the bare sites the transition state holds: θ₀^m forward, θ₀^(m−Δ) back)
        let sP = lnP, sQ = lnQ;
        for (const u of surf) (u.nu < 0 ? (sP -= u.nu * (lnA[u.s] + Math.log(bare))) : (sQ += u.nu * (lnA[u.s] + Math.log(bare))));
        if (rx.m !== undefined) {
          // (a face's reaction without surface species has no m; a port's always has, bare: true included)
          sP += rx.m * Math.log(bare);
          sQ += (rx.m - rx.fills) * Math.log(bare);
        }
        return mhcRate(rx.mhc, rx.k0, sP, sQ, a).rate;
      }
      return rx.vmax ? p * -Math.expm1(-a) : p * (Math.exp(rx.alpha * a) - Math.exp(-(1 - rx.alpha) * a));
    });
    // (its level with the surface as held)
    const aHeld = surf.reduce((t, u) => t - u.nu * (held ? held[u.s] : 0), base);
    return { equation: names?.[k]?.equation ?? `reaction ${k}`, level: s === 0 ? NaN : -aHeld / s, rate, current: rate.map((r) => q * FARADAY * r) };
  });
  const current = new Float64Array(Vs.length);
  for (const r of out) for (let j = 0; j < Vs.length; j++) current[j] += r.current[j];
  const coverage = surface && Object.fromEntries(surface.sp.map((sp, x) => [sp.name, Float64Array.from(at, ({ lnA, bare }) => bare * Math.exp(lnA[x]))]));
  return { V: Vs, current, ...(surface ? { bare: Float64Array.from(at, (x) => x.bare), coverage } : {}), reactions: out };
}

// The surface in equilibrium at each level: each reaction with surface species at zero affinity,
// ν·ln a = base + s·V (ln a = ln θ/θ₀, θ₀ the bare fraction), as far as they're independent (those
// making one species from the solution alone first, in order, then the rest), and what they leave
// free fixed by the totals they keep: a tethered couple's θ_R + θ_O, a species no reaction touches.
// Those free directions, n, are the ones every reaction's ν is orthogonal to, so n·θ is what the
// surface conserves, and with ln a = P(V) + Σ t_k n_k the totals are kept where t minimizes
// ln(1 + Σ e^{ln a}) − t·(n·θ_held): strictly convex (its Hessian, the covariance of n over the
// sites, is positive definite while any site is bare), so Newton with a line search finds it.
function surfaceEquilibrium(parts, surface, held) {
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
  const candidates = [...parts.filter((pt) => pt.surf.length === 1), ...parts.filter((pt) => pt.surf.length > 1)];
  for (const pt of candidates) {
    const r = new Float64Array(X);
    for (const u of pt.surf) r[u.s] += u.nu;
    const w = orth(r), n = Math.sqrt(dot(w, w));
    if (n <= 1e-9 * Math.sqrt(dot(r, r))) continue;
    q.push(w.map((x) => x / n));
    rows.push({ r, pt });
  }
  const free = [];
  for (let x = 0; x < X && q.length < X; x++) {
    const e = new Float64Array(X);
    e[x] = 1;
    const w = orth(e), n = Math.sqrt(dot(w, w));
    if (n <= 1e-9) continue;
    const b = w.map((v) => v / n);
    q.push(b);
    free.push(b);
  }
  const k = rows.length, d = free.length;
  // (the rows' Gram matrix, for the particular solution P = Rᵀ (R Rᵀ)⁻¹ b, the one with no free part)
  const gram = rows.map((a) => rows.map((b) => dot(a.r, b.r)));
  const theta = surface.theta, target = free.map((n) => dot(n, theta));
  let t = free.map((n) => dot(n, held.map((l) => (Number.isFinite(l) ? l : 0))));
  const lnAOf = (P, t) => P.map((p, x) => free.reduce((s, n, j) => s + t[j] * n[x], p));
  const lse = (lnA) => {
    const m = Math.max(0, ...lnA);
    return m + Math.log(lnA.reduce((s, l) => s + Math.exp(l - m), Math.exp(-m)));
  };
  return (v) => {
    const P = new Float64Array(X);
    if (k > 0) {
      const y = solveDense(gram, rows.map(({ pt }) => pt.base + pt.s * v));
      rows.forEach(({ r }, j) => {
        for (let x = 0; x < X; x++) P[x] += y[j] * r[x];
      });
    }
    const f = (t) => lse(lnAOf(P, t)) - dot(t, target);
    let fx = f(t);
    for (let it = 0; d > 0; it++) {
      if (it === 100) fail(`polarization: the surface's equilibrium at V = ${v} didn't converge`);
      const lnA = lnAOf(P, t), Z = lse(lnA), th = lnA.map((l) => Math.exp(l - Z));
      const m = free.map((n) => dot(n, th)), g = m.map((x, j) => x - target[j]);
      const H = free.map((a, i) => free.map((b, j) => a.reduce((s, ax, x) => s + ax * b[x] * th[x], 0) - m[i] * m[j]));
      const step = solveDense(H, g.map((x) => -x));
      if (Math.max(...step.map(Math.abs)) < 1e-12) break; // (Newton's step, near quadratic here)
      let lam = 1, tn, fn;
      for (;;) {
        tn = t.map((x, j) => x + lam * step[j]);
        fn = f(tn);
        if (fn <= fx + 1e-4 * lam * dot(g, step) || lam < 1e-12) break;
        lam /= 2;
      }
      [t, fx] = [tn, fn];
    }
    const lnA = lnAOf(P, t);
    return { lnA, bare: Math.exp(-lse(lnA)) };
  };
}
