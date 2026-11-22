// Neon's third review (21 Nov 2026): our fixes of its second, one owner for strictly neutral
// charge (a design, in ROADMAP), and its own pick, the time stepper and the statistics.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT } from '../src/index.js';

// R2 1.2's fix left the end box's storage and charge to be summed after the contact read its
// current only at an ideal node with ε > 0. A Bikerman cell (30 cm of 1 M KCl on a lattice, no
// port at all) read Im Z 8% low at 10 MHz and 2.3% off at 1 GHz, relError 1.4e-9; a strictly
// neutral one with a sense port, a pure resistance, read inductive at 1 GHz (6% off) and 52% off at
// 10 GHz, relError 3e-11; all unwarned (Neon's third review, R3 1.1). The cell's conductance doesn't
// depend on the statistics here (the flux is D c ∇μ̄/RT, the state uniform): Z = R/(1 + iωRC).
test("impedance() at a bath contact's end box of any statistics or ε keeps the cell's RC (Neon R3 1.1)", () => {
  const c = 1000, D = 2e-9;
  const cell = (epsr, L, cells, extra) => ({
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr, species: { 'K+': { D, mu0: 0 }, 'Cl-': { D, mu0: 0 } }, ...extra.material } },
    regions: [{ name: 'cell', material: 'aq', length: L, c0: { 'K+': c, 'Cl-': c } }],
    contacts: { left: { V: 0, bath: { c: { 'K+': c, 'Cl-': c } } }, right: { V: 0, bath: { c: { 'K+': c, 'Cl-': c } } } },
    ...(extra.port ? { ports: [{ name: 'sense', region: 'cell', from: 0, to: L, I: 0, species: { 'K+': { type: 'conductance', G: 1e-6, offset: 0 } } }] } : {}),
    grid: { hmax: L / cells, minCells: cells },
  });
  const sigma = (FARADAY * FARADAY * c * 2 * D) / (GAS_CONSTANT * 298.15);
  const cases = [
    { name: 'lattice, no port', def: cell(78.5, 0.3, 10, { material: { statistics: [{ type: 'lattice', species: ['K+', 'Cl-'], cMax: 55000 }] } }), L: 0.3, epsr: 78.5, fs: [1e7, 1e9] },
    { name: 'lattice, a sense port', def: cell(78.5, 0.5, 10, { port: true, material: { statistics: [{ type: 'lattice', species: ['K+', 'Cl-'], cMax: 55000 }] } }), L: 0.5, epsr: 78.5, fs: [1e4, 1e9] },
    { name: 'strictly neutral, a sense port', def: cell(0, 0.05, 8, { port: true, material: {} }), L: 0.05, epsr: 0, fs: [1e9, 1e10] },
  ];
  for (const { name, def, L, epsr, fs } of cases) {
    const R = L / sigma, C = (epsr * EPS0) / L;
    const z = new Device(def).impedance(fs, { terminal: 'right' });
    fs.forEach((f, q) => {
      const wrc = 2 * Math.PI * f * R * C, re = R / (1 + wrc * wrc), im = (-R * wrc) / (1 + wrc * wrc);
      const off = Math.hypot(z.Z.re[q] - re, z.Z.im[q] - im) / Math.hypot(re, im);
      assert.ok(off < 1e-6, `${name}, ${f} Hz: ${z.Z.re[q]} ${z.Z.im[q]}i against ${re} ${im}i (${off.toExponential(2)} off)`);
    });
  }
});

// A Donnan face between two strictly neutral NaCl regions (Na⁺ 5 kJ/mol lower on the right), a
// bath at the left, the right closed, and a floating gate (I = 0, 1 F/m²) over the right region.
// The first step from the start moved the right face node's φ̂ by the Donnan jump, −25.9 mV, its
// ions balanced alone: the floating gate, its voltage left where it was, took that as a charge,
// 1.6e-3 C/m², no ion there balancing it. The step drove 1.6 MA/m² in at the bath while the gate
// passed nothing, locked the gate at that charge for good, and an advance() from the start failed
// at t = 0 (Neon's third review, R3 1.7). Now its voltage follows the φ̂ moved under it.
test('a floating gate over a strictly neutral Donnan face keeps its charge through the first step (Neon R3 1.7)', () => {
  const def = {
    species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: {
      aqL: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } },
      aqR: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 5000 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } },
    },
    regions: [
      { name: 'L', material: 'aqL', length: 5e-7, c0: { 'Na+': 20, 'Cl-': 20 } },
      { name: 'R', material: 'aqR', length: 5e-7, c0: { 'Na+': 100, 'Cl-': 100 } },
    ],
    interfaces: [{ phi: 'neutral', species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium' } }],
    contacts: { left: { V: 0, bath: { c: { 'Na+': 20, 'Cl-': 20 } } }, right: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } } },
    ports: [{ name: 'fgate', region: 'R', from: 0, to: 5e-7, area: 1e6, I: 0, capacitance: { C: 1, zeroCharge: 0 } }],
    grid: { hmax: 1.25e-7, minCells: 4 },
  };
  const s = new Device(def).step(1e-9);
  assert.ok(s.converged);
  assert.ok(Math.abs(s.ports[0].charge) < 1e-15, `the gate's charge ${s.ports[0].charge}`);
  const sum = s.terminals.left.current + s.terminals.fgate.current;
  assert.ok(Math.abs(sum) < 1e-3, `the terminals' currents sum to ${sum} A/m²`);
  const a = new Device(def).advance(1e-3);
  assert.equal(a.warnings.length, 0, a.warnings.join('\n'));
  assert.ok(Math.abs(a.ports[0].charge) < 1e-15, `after 1 ms, the gate's charge ${a.ports[0].charge}`);
});

// A held gate over half of a strictly neutral NaCl film at a bath, and set() adding a second on the
// same window: the carry moved φ̂ under both to keep each node neutral, the second gate's charge
// read from that state, the first's left at what it held before. The first step read the
// difference as the first gate's current, 1.9e5 A/m², nothing else passing it, Kirchhoff broken,
// unwarned (Neon's third review, R3 1.5). Now a gate's hardware changed keeps the concentrations
// and each gate's charge (the new one's, none), φ̂ jumping under them (src/carry.js, src/neutral.js):
// the first step passes the charge between the two gates and the bath, the currents summing to
// zero, and its start is on neutral.
test('a gate added beside another over a strictly neutral film keeps the charges, and the first step Kirchhoff (Neon R3 1.5)', () => {
  const gate = { name: 'gate', region: 'el', from: 0, to: 0.5e-6, area: 1e6, V: 0.5, capacitance: { C: 1, zeroCharge: 0 } };
  const d = new Device({
    species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
    regions: [{ name: 'el', material: 'aq', length: 1e-6, c0: { 'Na+': 100, 'Cl-': 100 } }],
    contacts: { left: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } }, right: { V: 0, bath: { c: { 'Na+': 100, 'Cl-': 100 } } } },
    ports: [gate],
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  d.solve();
  const before = d.step(1e-9), c0 = Float64Array.from(d.solver.c);
  d.set({ ports: [gate, { ...gate, name: 'gate2', V: -0.3 }] });
  assert.ok(d.solver.c.every((c, k) => Math.abs(c - c0[k]) <= 1e-12 * c0[k]), 'the concentrations kept through the set()');
  assert.equal(d.solver.portQEnd[0], before.ports[0].charge);
  assert.equal(d.solver.portQEnd[1], 0);
  const s = d.step(1e-9);
  assert.equal(s.warnings.length, 0, s.warnings.join('\n'));
  const I = s.terminals.right.current + s.terminals.gate.current + s.terminals.gate2.current;
  const gross = Math.abs(s.terminals.right.current) + Math.abs(s.terminals.gate.current) + Math.abs(s.terminals.gate2.current);
  assert.ok(gross > 1e3 && Math.abs(I) < 1e-9 * gross, `the currents sum to ${I} of ${gross} A/m²`);
});
