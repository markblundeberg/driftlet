// Strictly neutral solutions (epsr 0) meeting at faces, sealed at both ends, under held gates:
// each gate a port whose capacitance takes a charge the solution under it balances. Two to three
// regions; a gate over a whole region or part of one; one or two salts; a face that passes every
// ion or blocks one. Every amount is closed, so the state isn't unique, and with nothing driven
// the steady state is an equilibrium. (Neon's first review found its cold start losing ions here
// at t = 0⁺, 5.2e-3 of the Na⁺ in a sealed cell, and a part-gated region starting off its level.)
import { FARADAY as F, GAS_CONSTANT } from '../../../src/index.js';

export default {
  name: 'neutral',
  make(r) {
    const T = 298.15, RT = GAS_CONSTANT * T;
    const divalent = r.chance(0.4);
    const ions = [['Na+', 1, 1.33e-9], ['Cl-', -1, 2.03e-9], ...(divalent ? [['Ca2+', 2, 0.79e-9]] : [])];
    const species = ions.map(([name, z]) => ({ name, z }));
    const D = Object.fromEntries(ions.map(([name, , d]) => [name, d * r.u(0.7, 1.4)]));
    const materials = { aq: { epsr: 0, species: Object.fromEntries(ions.map(([name]) => [name, { D: D[name], cRef: 1000, mu0: 0 }])) } };
    const nR = r.pick([2, 2, 3]), features = ['strictly neutral', 'sealed', divalent ? 'divalent' : 'one salt'];
    const regions = [], gated = [], neutralC = [];
    for (let g = 0; g < nR; g++) {
      const length = r.log(-6.5, -5.5), name = `r${g}`;
      const na = r.log(1, 2), ca = divalent ? r.log(0.5, 1.5) : 0;
      const gate = r.chance(0.55) || (g === nR - 1 && gated.length === 0);
      // (under a gate over the whole region, the solution needn't be neutral: the gate's charge
      // takes up what's over; drawn always, so whether the gate is whole doesn't move the draws)
      const over = r.u(-0.3, 0.3) * na;
      const whole = r.chance(0.65), excess = gate && whole ? over : 0;
      const c0 = { 'Na+': na, 'Cl-': na + 2 * ca - excess, ...(divalent ? { 'Ca2+': ca } : {}) };
      regions.push({ name, material: 'aq', length, c0 });
      neutralC.push({ 'Na+': na, 'Cl-': na + 2 * ca, ...(divalent ? { 'Ca2+': ca } : {}) });
      if (gate) gated.push({ g, whole, excess });
    }
    const interfaces = Array.from({ length: nR - 1 }, () => {
      const block = r.chance(0.25) ? r.pick(ions.map(([n]) => n)) : null;
      if (block) features.push(`face blocks ${block}`);
      return { phi: 'neutral', species: Object.fromEntries(ions.map(([n]) => [n, n === block ? 'blocked' : 'equilibrium'])) };
    });
    const sealed = () => ({ phi: 'neutral', species: Object.fromEntries(ions.map(([n]) => [n, 'blocked'])) });
    let part = false;
    const ports = gated.map(({ g, whole, excess }, k) => {
      const L = regions[g].length;
      if (!whole) part = true;
      const from = whole ? undefined : r.u(0, 0.4) * L, to = whole ? undefined : from + r.u(0.3, 0.6) * L;
      const p = { name: `g${k}`, region: regions[g].name, ...(whole ? {} : { from, to }), V: r.u(-0.3, 0.3), area: r.log(5, 6), capacitance: { C: r.u(0.2, 1) } };
      // (what's over, no more than the gate takes 0.5 V from its V: drawn freely, a gate 8 V
      // from its solution under it took Newton past what it holds; see ROADMAP)
      const most = (0.5 * p.area * p.capacitance.C) / F;
      if (Math.abs(excess) > most) regions[g].c0['Cl-'] += excess - Math.sign(excess) * most;
      return p;
    });
    features.push(part ? 'a gate over part of a region' : 'gates over whole regions', `${ports.length} gate${ports.length > 1 ? 's' : ''}`);
    // Drawn last, so that the cases above keep their draws: a Donnan face (Na⁺ lower in μ° in the
    // last region, so the face takes a jump in φ at t = 0⁺), a bath at the left, a floating gate
    // (at no current: its charge stays, its voltage follows). (Neon's R3 1.7: a floating gate over a
    // Donnan face's right side took the jump as a charge nothing balanced.)
    const contacts = { left: sealed(), right: sealed() };
    if (r.chance(0.35)) {
      materials.aqD = structuredClone(materials.aq);
      materials.aqD.species['Na+'].mu0 = r.u(2000, 6000);
      regions[nR - 1].material = 'aqD';
      features.push('Donnan face');
    }
    // (not beside a whole gate's excess: the bath's end box, whose φ the bath sets, holds no gate's charge)
    if (r.chance(0.3) && !(gated[0]?.g === 0 && gated[0].excess !== 0)) {
      contacts.left = { V: 0, bath: { c: { ...neutralC[0] } } };
      features.push('bath at the left');
    }
    if (r.chance(0.4)) {
      const g = r.int(0, nR - 1), L = regions[g].length, whole = r.chance(0.5);
      const from = whole ? undefined : r.u(0, 0.4) * L, to = whole ? undefined : from + r.u(0.3, 0.6) * L;
      ports.push({ name: 'float', region: regions[g].name, ...(whole ? {} : { from, to }), I: 0, area: r.log(5, 6), capacitance: { C: r.u(0.2, 1) } });
      features.push('a floating gate');
    }
    const def = { T, species, materials, regions, interfaces, contacts, ports, grid: { minCells: 8 } };
    // The slowest relaxation: diffusion over the whole cell, or a gate charging through the
    // solution under it (C per volume over σ/L²).
    const Ltot = regions.reduce((t, g) => t + g.length, 0), Dmin = Math.min(...Object.values(D));
    const rc = ports.map((p) => {
      const g = regions.find((x) => x.name === p.region), sigma = (F * F * Math.min(...Object.values(g.c0)) * Dmin) / RT;
      return (p.area * p.capacitance.C * Ltot * Ltot) / sigma;
    });
    const slow = Math.max((Ltot * Ltot) / Dmin, ...rc);
    return {
      def,
      params: { features, regions: regions.map((g) => ({ L: g.length, c0: g.c0 })), gates: ports.map((p) => ({ region: p.region, V: p.V, part: p.from !== undefined })) },
      plan: { terminal: ports[0].name, V: ports[0].V + r.u(-0.05, 0.05), dV: 0.01, equilibrium: true, unique: false, slow },
    };
  },
};
