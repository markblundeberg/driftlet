// The one owner of strictly neutral charge (src/neutral.js): its projections, and the producers
// that use it (the carry, the transient's cold start). bench/combos' neutralG property checks the
// same after every producer, on every family.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { gaussRows, projectNeutral, strictlyNeutralNodes } from '../src/neutral.js';
import { equilibrateFaces } from '../src/start.js';

const worstOff = (s) => {
  const nodes = strictlyNeutralNodes(s);
  return Math.max(0, ...gaussRows(s, nodes).map((r) => Math.abs(r.G) / r.scale));
};

// A strictly neutral NaCl film under a held gate over half of it, at a bath.
const film = () => ({
  species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
  materials: { aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
  regions: [{ name: 'el', material: 'aq', length: 1e-6, c0: { 'Na+': 100, 'Cl-': 100 } }],
  contacts: { left: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } }, right: { V: 0, bath: { c: { 'Na+': 100, 'Cl-': 100 } } } },
  ports: [{ name: 'gate', region: 'el', from: 0, to: 0.5e-6, area: 1e6, V: 0.3, capacitance: { C: 1, zeroCharge: 0 } }],
  grid: { hmax: 1.25e-7, minCells: 8 },
});

test('projectNeutral brings a node back onto its Gauss row by the unknowns its mode names, and leaves one on it bit for bit', () => {
  const d = new Device(film());
  d.solve();
  const s = d.solver, { M, n } = s, g = 2, b = s.blockOfNode[g] * M;
  const u0 = Float64Array.from(s.u);
  projectNeutral(s, strictlyNeutralNodes(s));
  assert.ok(s.u.every((x, k) => Object.is(x, u0[k])), 'a solved state, left as it was');
  for (const mode of ['charge', 'level', 'jump']) {
    s.u.set(u0);
    s.u[b + 1] += 1e-3; // (Na⁺ at node 2, off neutral by 0.1%)
    s.computeConcentrations();
    const phi = s.u[b], eta = [s.u[b + 1], s.u[b + 2]], c = [s.c[g * n], s.c[g * n + 1]];
    assert.ok(worstOff(s) > 1e-4);
    projectNeutral(s, strictlyNeutralNodes(s).map((nd) => ({ ...nd, mode })));
    assert.ok(worstOff(s) < 1e-14, `${mode}: ${worstOff(s)} off`);
    if (mode === 'charge') assert.equal(s.u[b], phi, 'charge: φ̂ kept');
    if (mode === 'level') assert.deepEqual([s.u[b + 1], s.u[b + 2]], eta, 'level: the levels kept');
    if (mode === 'jump') for (let i = 0; i < 2; i++) assert.ok(Math.abs(s.c[g * n + i] / c[i] - 1) < 1e-13, 'jump: the concentrations kept');
  }
});

// Neon's R2 3.1 device (test/neon2.test.js): Pt | 1 µm of strictly neutral NaCl | Pt, capacitive
// faces, regridded mid-transient. Its ions are closed, so their net charge is kept, and so are the
// faces' charges that balance it: φ̂ stays where it was carried, at the faces too, and the
// concentrations move in the charge direction. (Held to each node's old net charge by φ̂ at the
// levels, the edge's φ̂ moved 1.08 mV and the face's D 2.16e-4 C/m², the first step charging it.)
test('a regrid of a closed strictly neutral region keeps its faces\' D, every node on its Gauss row', () => {
  const d = new Device({
    species: [{ name: 'e-', z: -1, cRef: 1000 }, { name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } }, aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
    regions: [{ name: 'L', material: 'Pt', length: 1e-7 }, { name: 'el', material: 'aq', length: 1e-6, c0: { 'Na+': 100, 'Cl-': 100 } }, { name: 'R', material: 'Pt', length: 1e-7 }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }, { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }],
    contacts: { left: { V: 0.5, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, right: { V: 0, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } } },
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  d.solve();
  d.set({ contacts: { left: { V: 0.51 } } });
  const t1 = d.advance(1e-4);
  d.set({ grid: { hmax: 6.25e-8, minCells: 16 } });
  const c = d.solution();
  for (const f of [0, 1]) assert.ok(Math.abs(c.interfaces[f].D / t1.interfaces[f].D - 1) < 1e-12, `face ${f}: D ${t1.interfaces[f].D} → ${c.interfaces[f].D}`);
  assert.ok(worstOff(d.solver) < 1e-12, `${worstOff(d.solver)} off`);
});

// A capacitive face holding a sheet of charge between two strictly neutral NaCl regions, a bath
// at the left: at c0 the start has no neutral state (the left box would need D = 0, the right
// one D = −σ), so its edge boxes hold the sheet's counter-charge, by φ̂ at their levels, before
// the first step. (Left off, the first step charged them through the boxes beside them.)
test('a transient\'s cold start puts a sheet-charged face\'s counter-charge in its edge boxes', () => {
  const salt = { 'Na+': 100, 'Cl-': 100 };
  const d = new Device({
    species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'aq', length: 1e-6, c0: salt }, { name: 'b', material: 'aq', length: 1e-6, c0: salt }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.01 }, sheetCharge: 0.01, species: { 'Na+': { type: 'permeability', P: 1e-6 }, 'Cl-': { type: 'permeability', P: 1e-6 } } }],
    contacts: { left: { V: 0, bath: { c: salt } }, right: { V: 0, bath: { c: salt } } },
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  const s = d.solver;
  assert.ok(worstOff(s) > 1e-3, 'off at c0');
  equilibrateFaces(s);
  assert.ok(worstOff(s) < 1e-12, `${worstOff(s)} off`);
  assert.ok(d.step(1e-9).converged);
});

// 1 M of NaCl, closed, under a held gate, with H⁺ at a trace open to a contact: its node took φ̂
// at the levels, the closed amounts and the gate chased each other through the carry's 60 rounds,
// and a regrid left Na⁺ and Cl⁻ each 2.4e-10 off what they held. Where the open ions are a trace
// of the screening, the closed ones move in the charge direction alone, the open ones keeping their
// levels. (The reviewer of Neon's third round, reading our account of the owner.)
test('a regrid keeps closed amounts exactly beside an open trace ion', () => {
  for (const trace of [1e-3, 10]) {
    const d = new Device({
      species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'H+', z: 1, cRef: 1000 }],
      materials: { aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 }, 'H+': { D: 9.3e-9, mu0: 0 } } } },
      regions: [{ name: 'el', material: 'aq', length: 1e-6, c0: { 'Na+': 1000, 'Cl-': 1000 + trace, 'H+': trace } }],
      contacts: { left: { V: 0, terminal: 'H+', species: { 'H+': 'equilibrium', 'Na+': 'blocked', 'Cl-': 'blocked' }, phi: 'bulk' }, right: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked', 'H+': 'blocked' } } },
      ports: [{ name: 'gate', region: 'el', from: 0.3e-6, to: 0.8e-6, area: 1e6, V: 0.1, capacitance: { C: 1, zeroCharge: 0 } }],
      grid: { hmax: 1.25e-7, minCells: 8 },
    });
    d.solve();
    const before = d.step(1e-9).conservation;
    d.set({ grid: { hmax: 6.25e-8, minCells: 16 } });
    const after = d.solution().conservation;
    for (const sp of ['Na+', 'Cl-']) {
      const a = before.find((st) => st.species === sp).amount, b = after.find((st) => st.species === sp).amount;
      assert.ok(Math.abs(b / a - 1) < 1e-14, `trace ${trace} mM, ${sp}: ${a} → ${b}`);
    }
    assert.ok(worstOff(d.solver) < 1e-12);
  }
});
