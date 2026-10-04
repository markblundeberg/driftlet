// Polarization curves at an electrode, for driftlet/kit: each reaction's rate with the solution
// beside the electrode held as it is and the metal's level moved, which is what an Evans diagram
// plots. Where the net current crosses zero is the spot's open-circuit (mixed) potential.

import { DeviceError } from './errors.js';
import { FARADAY } from './constants.js';
import { powi, powr } from './pow.js';

const fail = (message) => {
  throw new DeviceError(message);
};

/**
 * Each electrode reaction's rate at a face (a metal's) or at one spot of an electrode port, with
 * the composition beside it as in `sol` and the metal's carrier at each level in `V` (in the
 * terminals' convention, V = μ̄/(zF); against SHE, subtract φ for the usual tables' ions). The
 * rates are the solver's own law, so at the electrode's actual level they are
 * `sol.interfaces[f].rates` (or `sol.ports[k].rates` at that node).
 * @param {import('./index.js').Device} device
 * @param {object} sol a solution of it
 * @param {{ face: number } | { port: string | number, x?: number }} where a face's index, or a port
 *   (name or index) and the position in its window (m, the device's x; default its first node)
 * @param {number | ArrayLike<number>} V the metal's level(s), V
 * @param {{ surface?: 'held' | 'equilibrium' }} [opts] an electrode port's surface coverages: held as
 *   in `sol` (the default: this instant's curves), or in equilibrium at each V (the steady-state
 *   curve, with a passive film's active–passive peak), each species by the reaction that makes it
 * @returns {{ V: Float64Array, current: Float64Array, bare?: Float64Array, reactions: { equation: string, level: number,
 *   rate: Float64Array, current: Float64Array }[] }} per reaction: its equilibrium level here (where
 *   its rate is zero), its rate (mol/(m²·s) of electrode, forward) and its current (A/m², anodic
 *   positive: into the solution); and their sum
 */
export function polarization(device, sol, where, V, opts = {}) {
  const model = device.model, { species, grid, RT } = model, solver = device.solver;
  const Vs = Float64Array.from(typeof V === 'number' ? [V] : V);
  let reactions, nodeOf, carrier, names, surface = null; // (an electrode port's coverages at the spot)
  if (where?.face !== undefined) {
    const f = where.face, itf = model.interfaces[f];
    if (!itf) fail(`polarization: no face ${f} (the device has ${model.interfaces.length})`);
    const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
    const metal = [gL, gR].findIndex((g) => solver.nodeConductor[g] >= 0);
    if (metal < 0) fail(`polarization: face ${f} has no metal on either side`);
    carrier = { side: metal, i: solver.nodeConductor[metal ? gR : gL] };
    nodeOf = (p) => (p.side ? gR : gL);
    reactions = itf.reactions;
    names = device.def.interfaces?.[f]?.reactions;
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
    let base = rx.fixedA, s = 0, q = 0, pref = rx.vmax ? rx.vmax : rx.k0;
    const surf = [];
    rx.part.forEach((p, x) => {
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
      const g = nodeOf(p), name = species[p.i].name;
      base -= (p.nu * sol.mu[name][g]) / RT;
      q += p.side === carrier.side ? 0 : p.nu * species[p.i].z; // charge the forward reaction puts into the solution
      const c = sol.c[name][g];
      if (rx.vmax) {
        if (p.nu < 0) pref *= powi(c / (c + rx.K[x]), -p.nu);
      } else pref *= powr(c / solver.cRef[g * species.length + p.i], p.nu < 0 ? -p.nu * (1 - rx.alpha) : p.nu * rx.alpha);
    });
    if (s === 0) fail(`polarization: reaction ${k} there takes no ${species[carrier.i].name} from the metal`);
    return { rx, k, base, s, q, pref, surf };
  });
  // The surface's activities (ln a) and bare fraction at a level: held as `sol` has them, or in
  // equilibrium there, each species by the first reaction that makes it from the solution alone.
  const held = surface && surface.theta.map((t) => Math.log(t / surface.bare));
  let former = null;
  if (surface && opts.surface === 'equilibrium') {
    former = surface.sp.map((sp, x) => {
      const f = parts.find((pt) => pt.surf.length === 1 && pt.surf[0].s === x);
      if (!f) fail(`polarization: no reaction makes ${sp.name} from the solution alone, so its equilibrium coverage isn't defined`);
      return f;
    });
  } else if (opts.surface !== undefined && opts.surface !== 'held') fail("polarization: surface must be 'held' (the default) or 'equilibrium'");
  const surfaceAt = (v) => {
    if (!surface) return { lnA: [], bare: 1 };
    // base + s v − ν ln a = 0 for each species' former
    const lnA = former ? former.map((f) => (f.base + f.s * v) / f.surf[0].nu) : held;
    if (!former) return { lnA, bare: surface.bare };
    const m = Math.max(0, ...lnA);
    let sum = Math.exp(-m);
    for (const l of lnA) sum += Math.exp(l - m);
    return { lnA, bare: Math.exp(-m) / sum };
  };
  const at = Array.from(Vs, (v) => surfaceAt(v));
  const out = parts.map(({ rx, k, base, s, q, pref, surf }) => {
    const rate = Float64Array.from(Vs, (v, j) => {
      const { lnA, bare } = at[j];
      let a = base + s * v, p = pref * (rx.bare ? bare : 1);
      for (const u of surf) {
        a -= u.nu * lnA[u.s];
        p *= Math.exp(u.e * lnA[u.s]);
      }
      return rx.vmax ? p * -Math.expm1(-a) : p * (Math.exp(rx.alpha * a) - Math.exp(-(1 - rx.alpha) * a));
    });
    // (its level with the surface as held)
    const aHeld = surf.reduce((t, u) => t - u.nu * (held ? held[u.s] : 0), base);
    return { equation: names?.[k]?.equation ?? `reaction ${k}`, level: -aHeld / s, rate, current: rate.map((r) => q * FARADAY * r) };
  });
  const current = new Float64Array(Vs.length);
  for (const r of out) for (let j = 0; j < Vs.length; j++) current[j] += r.current[j];
  return { V: Vs, current, ...(surface ? { bare: Float64Array.from(at, (x) => x.bare) } : {}), reactions: out };
}
