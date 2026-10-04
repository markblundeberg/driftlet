// Polarization curves at an electrode, for driftlet/kit: each reaction's rate with the solution
// beside the electrode held as it is and the metal's level moved, which is what an Evans diagram
// plots. Where the net current crosses zero is the spot's open-circuit (mixed) potential.

import { DeviceError } from './errors.js';
import { FARADAY } from './constants.js';

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
 * @returns {{ V: Float64Array, current: Float64Array, reactions: { equation: string, level: number,
 *   rate: Float64Array, current: Float64Array }[] }} per reaction: its equilibrium level here (where
 *   its rate is zero), its rate (mol/(m²·s) of electrode, forward) and its current (A/m², anodic
 *   positive: into the solution); and their sum
 */
export function polarization(device, sol, where, V) {
  const model = device.model, { species, grid, RT } = model, solver = device.solver;
  const Vs = Float64Array.from(typeof V === 'number' ? [V] : V);
  let reactions, nodeOf, carrier, names;
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
    reactions = port.reactions;
    names = device.def.ports?.[k]?.reactions;
  } else fail('polarization: say where, { face: f } or { port: name, x }');
  if (reactions.length === 0) fail('polarization: no reactions there');

  const out = reactions.map((rx, k) => {
    if (rx.srh) fail('polarization: SRH kinetics has no electrode level');
    // a = A/RT = a0 + s·V: everything but the carrier at its level, and the carrier's part.
    let a0 = rx.fixedA, s = 0, q = 0, pref = rx.vmax ? rx.vmax : rx.k0;
    rx.part.forEach((p, x) => {
      if (p.side === carrier.side && p.i === carrier.i) {
        s -= (p.nu * species[p.i].z * FARADAY) / RT;
        return;
      }
      const g = nodeOf(p), name = species[p.i].name;
      a0 -= (p.nu * sol.mu[name][g]) / RT;
      q += p.side === carrier.side ? 0 : p.nu * species[p.i].z; // charge the forward reaction puts into the solution
      const c = sol.c[name][g];
      if (rx.vmax) {
        if (p.nu < 0) pref *= (c / (c + rx.K[x])) ** -p.nu;
      } else pref *= (c / solver.cRef[g * species.length + p.i]) ** (p.nu < 0 ? -p.nu * (1 - rx.alpha) : p.nu * rx.alpha);
    });
    if (s === 0) fail(`polarization: reaction ${k} there takes no ${species[carrier.i].name} from the metal`);
    const rate = Float64Array.from(Vs, (v) => {
      const a = a0 + s * v;
      return rx.vmax ? pref * -Math.expm1(-a) : pref * (Math.exp(rx.alpha * a) - Math.exp(-(1 - rx.alpha) * a));
    });
    return { equation: names?.[k]?.equation ?? `reaction ${k}`, level: -a0 / s, rate, current: rate.map((r) => q * FARADAY * r) };
  });
  const current = new Float64Array(Vs.length);
  for (const r of out) for (let j = 0; j < Vs.length; j++) current[j] += r.current[j];
  return { V: Vs, current, reactions: out };
}
