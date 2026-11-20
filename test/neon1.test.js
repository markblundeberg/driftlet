// Neon's first review (20 Nov 2026), area 1: the cold start (src/start.js).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY } from '../src/index.js';
import { equilibrateFaces } from '../src/start.js';
import { build, layer, metal, ohmic, check } from '../src/kit.js';

// Two strictly neutral solutions meeting at a face that passes both ions, sealed, one or both
// under a held gate. The face's t = 0⁺ equilibrium (equilibrateFaces) left the gate's charge out
// of the nodes' neutrality: a gated left node's counter-ions were moved off it, and the first
// step, holding each node neutral, lost them, 5.2e-3 of the ions in a sealed cell, with no
// warning (check() alone said so).
const gated = (regions) => ({
  species: [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }],
  materials: { aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, cRef: 1000, mu0: 0 }, 'Cl-': { D: 2.03e-9, cRef: 1000, mu0: 0 } } } },
  regions: [
    { name: 'L', material: 'aq', length: 1e-6, c0: regions.includes('L') ? { 'Na+': 20, 'Cl-': 10 } : { 'Na+': 20, 'Cl-': 20 } },
    { name: 'R', material: 'aq', length: 1e-6, c0: regions.includes('R') ? { 'Na+': 100, 'Cl-': 110 } : { 'Na+': 100, 'Cl-': 100 } },
  ],
  interfaces: [{ phi: 'neutral', species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium' } }],
  contacts: { left: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } }, right: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } } },
  ports: regions.map((region, k) => ({ name: `g${k}`, region, V: -0.5 + 0.3 * k, area: 1e6, capacitance: { C: 1 } })),
  grid: { minCells: 8 },
});
for (const regions of [['L'], ['R'], ['L', 'R']]) {
  test(`a gated face's t = 0⁺ keeps every ion and the gates' charge (gates on ${regions.join(' and ')}; Neon 1.1)`, () => {
    const steady = new Device(gated(regions)).solve();
    const d = new Device(gated(regions));
    const t = d.advance(0.1);
    assert.ok(t.converged);
    for (let k = 0; k < regions.length; k++) {
      const Q = t.ports[k].charge, Qs = steady.ports[k].charge;
      assert.ok(Math.abs(Q - Qs) < 1e-9 * Math.abs(Qs), `gate ${k}: ${Q} after 0.1 s, ${Qs} C/m² steady`);
    }
    assert.match(check(d, t).text, /ok {3}conservation/);
    // and the face's t = 0⁺ (equilibrateFaces) is neutral at every node, ρ + aσ = 0, gates and all
    const e = new Device(gated(regions)), def = gated(regions);
    e.solution();
    const sv = e.solver, { grid } = sv.model;
    equilibrateFaces(sv);
    for (let g = 0; g < grid.nNodes; g++) {
      const r = grid.nodeRegion[g], b = sv.blockOfNode[g];
      let rho = FARADAY * (sv.c[2 * g] - sv.c[2 * g + 1]);
      def.ports.forEach((p) => {
        if (p.region === def.regions[r].name) rho += p.area * p.capacitance.C * (p.V - sv.VT * sv.u[b * sv.M]);
      });
      assert.ok(Math.abs(rho) < 1e-9 * FARADAY * (sv.c[2 * g] + sv.c[2 * g + 1]), `node ${g} (${def.regions[r].name}): ρ + aσ = ${rho} C/m³ at t = 0⁺`);
    }
  });
}

// Area 2: a group of stretches held flat at its one contact's levels as an equilibrium, where a
// reaction among surface species alone (S = F, a Tafel step) keeps it from one: the check left
// those reactions out, and the surface's terms in the others, so A = S, S = F with F held 25
// kJ/mol below A solved to no flux and A flat at its bath's 1000 mol/m³, with no warning, while
// each time step turned A into F at 9.8e-3 mol/(m²·s) (Neon 2.1).
test('a catalytic face, A = S then S = F to a held product: the steady solve carries what the transient settles to, and only an equilibrium is held flat (Neon 2.1)', () => {
  const lib = { species: [{ name: 'A', z: 0, cRef: 1000 }], materials: { aq: { epsr: 0, species: { A: { D: 1e-9, mu0: 0 } } } } };
  const cell = (dG) =>
    build({
      T: 298.15, library: [lib, metal('Pt')],
      stack: [
        { phi: 'neutral', species: { A: { type: 'equilibrium', mu: 0 } } },
        layer('aq', 1e-6),
        { surface: { S: { mu0: 0, capacity: 1e-5, theta0: 0.1 } }, reactions: [{ equation: 'A = S', k0: 1e-2 }, { equation: 'S = F', fixed: { F: dG }, k0: 1e-2 }] },
        layer('Pt', 1e-7),
        ohmic(0, ['e-']),
      ],
      grid: { minCells: 10 },
    });
  const d = new Device(cell(-25000)), s = d.solve(), t = new Device(cell(-25000)).advance(10);
  assert.ok(s.converged && t.converged);
  const [r, rt] = [s.interfaces[0].rates[1], t.interfaces[0].rates[1]];
  assert.ok(r > 9e-3 && Math.abs(r / rt - 1) < 1e-9, `steady ${r}, after 10 s ${rt} mol/(m²·s)`);
  assert.match(check(d, s).text, /ok {3}balance/);
  const eq = new Device(cell(0)); // (at F's own level the group is at equilibrium: still held flat)
  assert.deepEqual(eq.solver.flatStretches.map((st) => st.species), [0]);
  assert.equal(eq.solve().interfaces[0].rates[1], 0);
});
