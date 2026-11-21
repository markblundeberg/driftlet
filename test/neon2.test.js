// Neon's second review (20 Nov 2026): impedance, check(), and its own pick.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { check } from '../src/kit.js';

// Two ports pumping K⁺ in and out of a KCl film at ±1 A/m², Cl⁻ held at both contacts: the K⁺
// law's change is the two ports' together, so a time step reads it by the flux (circuitLaws), and
// impedance() cleared that flux reading for a law no change rows had replaced. Read at a driven
// port, Re Z = −19.3 Ω·m² at 1 mHz (held at its voltage, 0.0036), relError 5e-13, no warning; read
// at a contact, Z moved 9% with the ports' order (a transient sine at 0.1 Hz, its step refined, goes
// to 0.002461 − 0.000636i).
const pumped = (order) => ({
  species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
  materials: { water: { epsr: 78.5, species: { 'K+': { D: 1.96e-9, mu0: 0 }, 'Cl-': { D: 1e-11, mu0: 0 } } } },
  regions: [{ name: 'aq', material: 'water', length: 1e-5, c0: { 'K+': 100, 'Cl-': 100 } }],
  contacts: {
    left: { V: 0, phi: 'bulk', species: { 'Cl-': { type: 'equilibrium', offset: 0 }, 'K+': 'blocked' } },
    right: { V: 0, phi: 'bulk', species: { 'Cl-': { type: 'equilibrium', offset: 0 }, 'K+': 'blocked' } },
  },
  ports: order.map((k) => ({ name: `p${k}`, region: 'aq', from: (k === 1 ? 0.1 : 0.6) * 1e-5, to: (k === 1 ? 0.4 : 0.9) * 1e-5, area: 1e6, capacitance: { C: 0.2, zeroCharge: 0 }, I: k === 1 ? 1 : -1, species: { 'K+': { type: 'conductance', G: 1e8, offset: 0 } } })),
  grid: { hmax: 5e-7 },
});
test("impedance() at a port driven beside another driving the same ion is the port's held impedance, and a contact's doesn't hang on the ports' order (Neon R2 1.1)", () => {
  const fs = [1e-3, 0.1, 10];
  const d = new Device(pumped([1, 2])), s = d.solve(), driven = d.impedance(fs, { terminal: 'p1' });
  const h = new Device(pumped([1, 2]));
  h.solve();
  h.set({ ports: { p1: { V: s.terminals.p1.V } } });
  const held = h.impedance(fs, { terminal: 'p1' });
  fs.forEach((f, q) => {
    const dz = Math.hypot(driven.Z.re[q] - held.Z.re[q], driven.Z.im[q] - held.Z.im[q]), z = Math.hypot(held.Z.re[q], held.Z.im[q]);
    assert.ok(dz < 1e-6 * z && driven.Z.re[q] > 0, `${f} Hz: driven ${driven.Z.re[q]} ${driven.Z.im[q]}i, held ${held.Z.re[q]} ${held.Z.im[q]}i`);
  });
  const [a, b] = [[1, 2], [2, 1]].map((o) => new Device(pumped(o)).impedance([0.1], { terminal: 'right' }).Z);
  assert.ok(Math.abs(a.re[0] / b.re[0] - 1) < 1e-6 && Math.abs(a.im[0] / b.im[0] - 1) < 1e-6, `right: ${a.re[0]} ${a.im[0]}i, ports reversed ${b.re[0]} ${b.im[0]}i`);
  assert.ok(Math.abs(a.re[0] / 0.002461 - 1) < 2e-3 && Math.abs(a.im[0] / -0.000636 - 1) < 3e-3, `right: ${a.re[0]} ${a.im[0]}i against the transient's 0.002461 − 0.000636i`);
});

// 30 cm of 1 M KCl between baths, a sense port over all of it at I: 0 (so the current at a contact
// is read from its end box, no cut between it and a window). The box's ions' storage, in its
// balances, and their charge, in its Gauss row at 1/dt, cancel in the current; summed at
// dt = 1e-30 they left a residue of one ULP of v·F·c = 1.45e6 C/m² on the η columns, against a
// displacement ε V_T/h of 6e-10: Im Z at 10 kHz came out 9× the RC cell's and inductive, 104% off
// at 1 GHz, relError 1e-13, unwarned (Neon's second review, R2 1.2). Without the port, read across
// a segment, it was right.
test("impedance() read in a bath contact's end box, a port's window over it, keeps the cell's displacement (Neon R2 1.2)", () => {
  const c = 1000, D = 2e-9, epsr = 78.5, L = 0.3;
  const def = (port) => ({
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr, species: { 'K+': { D, mu0: 0 }, 'Cl-': { D, mu0: 0 } } } },
    regions: [{ name: 'cell', material: 'aq', length: L, c0: { 'K+': c, 'Cl-': c } }],
    contacts: { left: { V: 0, bath: { c: { 'K+': c, 'Cl-': c } } }, right: { V: 0, bath: { c: { 'K+': c, 'Cl-': c } } } },
    ...(port ? { ports: [{ name: 'sense', region: 'cell', from: 0, to: L, I: 0, species: { 'K+': { type: 'conductance', G: 1e-6, offset: 0 } } }] } : {}),
    grid: { hmax: L / 10 },
  });
  // The cell: R = L/σ, C = ε/L, Z = R/(1 + iωRC).
  const sigma = (FARADAY * FARADAY * c * 2 * D) / (GAS_CONSTANT * 298.15), R = L / sigma, C = (epsr * EPS0) / L;
  const fs = [1e4, 1e9];
  for (const port of [false, true]) {
    const z = new Device(def(port)).impedance(fs, { terminal: 'right' });
    fs.forEach((f, q) => {
      const wrc = 2 * Math.PI * f * R * C, re = R / (1 + wrc * wrc), im = (-R * wrc) / (1 + wrc * wrc);
      assert.ok(Math.abs(z.Z.re[q] / re - 1) < 1e-4 && Math.abs(z.Z.im[q] / im - 1) < 1e-3, `${port ? 'port' : 'no port'}, ${f} Hz: ${z.Z.re[q]} ${z.Z.im[q]}i against ${re} ${im}i`);
    });
  }
});

// Two strictly neutral KCl regions across a membrane (K⁺ at P = 1e-15 m/s, Cl⁻ blocked, 1 mF/m²),
// a bath at the left, a blocked right, and a held port in the right region: a two-terminal cell,
// one Z from either end. Measured at the port, the uniform-shift split that resolves it at the
// contact never ran (it waited on a contact's cuts): relError 2.2e3 at 1e-5 Hz and "isn't resolved"
// at every frequency, where at the contact it was 2e-6 (Neon's second review, R2 1.3).
test('impedance() at a held port with nothing floating splits off the uniform shift as at a contact (Neon R2 1.3)', () => {
  const def = {
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr: 0, species: { 'K+': { D: 1.96e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
    regions: [
      { name: 'L', material: 'aq', length: 1e-5, c0: { 'K+': 100, 'Cl-': 100 } },
      { name: 'R', material: 'aq', length: 1e-5, c0: { 'K+': 100, 'Cl-': 100 } },
    ],
    interfaces: [{ phi: { type: 'capacitive', C: 1e-3 }, species: { 'K+': { type: 'permeability', P: 1e-15 }, 'Cl-': 'blocked' } }],
    contacts: { left: { V: 0, bath: { c: { 'K+': 100, 'Cl-': 100 } } }, right: { phi: 'neutral', species: { 'K+': 'blocked', 'Cl-': 'blocked' } } },
    ports: [{ name: 'we', region: 'R', from: 0.4e-5, to: 0.6e-5, V: 0, species: { 'K+': { type: 'equilibrium', offset: 0 } } }],
    grid: { hmax: 1e-6 },
  };
  const fs = [1e-8, 1e-5, 1e-2], d = new Device(def);
  const left = d.impedance(fs, { terminal: 'left' }), we = d.impedance(fs, { terminal: 'we' });
  assert.deepEqual(we.warnings ?? [], []);
  fs.forEach((f, q) => {
    const dz = Math.hypot(we.Z.re[q] - left.Z.re[q], we.Z.im[q] - left.Z.im[q]) / Math.hypot(left.Z.re[q], left.Z.im[q]);
    assert.ok(we.relError[q] < 1e-2 && dz < 3 * we.relError[q] + 3 * left.relError[q], `${f} Hz: at we ${we.Z.re[q]} ${we.Z.im[q]}i (relError ${we.relError[q]}), at left ${left.Z.re[q]} ${left.Z.im[q]}i`);
  });
});

// Area 2: check()'s grid item. A double layer's charge was compared only as its region's net
// charge, and only where that passed a millionth of the region's gross ionic content: on 1 mm of
// 1 M KCl the gate was 0.26 C/m², over the electrode's 0.058, and D_left came out 6.6% off the
// fine grid's, grid ok; in a symmetric 1 µm cell the walls' layers cancel in the net charge, and
// D_left was 14% off, grid ok (Neon's second review, R2 2.1). Now a contact's D, a capacitive
// face's and a port's charge are compared as they are.
const gridItem = (d, s) => check(d, s).items.find((it) => it.name === 'grid');
test("check()'s grid item compares a contact's charge, not only its region's (Neon R2 2.1)", () => {
  const kcl = (D) => ({ species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }], materials: { aq: { epsr: 78.5, species: { 'K+': { D, mu0: 0 }, 'Cl-': { D: D === 2e-9 ? 2e-9 : 2.03e-9, mu0: 0 } } } } });
  const asym = {
    ...kcl(1.96e-9),
    regions: [{ name: 'el', material: 'aq', length: 1e-3, c0: { 'K+': 1000, 'Cl-': 1000 } }],
    contacts: {
      left: { V: 0.05, phi: { type: 'capacitive', C: 2.0 }, zeroCharge: 0, species: { 'K+': { type: 'conductance', G: 1e-2, offset: 0 }, 'Cl-': 'blocked' } },
      right: { V: 0, bath: { c: { 'K+': 1000, 'Cl-': 1000 } } },
    },
    grid: { hmin: 2.5e-10, hmax: 5e-5, ratio: 2.5 },
  };
  const off = ((GAS_CONSTANT * 298.15) / FARADAY) * Math.log(10 / 1000);
  const link = { 'K+': { type: 'conductance', G: 1e-2, offset: off }, 'Cl-': { type: 'conductance', G: 1e-2, offset: -off } };
  const sym = {
    ...kcl(2e-9),
    regions: [{ name: 'el', material: 'aq', length: 1e-6, c0: { 'K+': 10, 'Cl-': 10 } }],
    contacts: { left: { V: 0.01, phi: 'pinned', zeroCharge: 0, species: link }, right: { V: -0.01, phi: 'pinned', zeroCharge: 0, species: link } },
    grid: { hmin: 2.5e-9, hmax: 2e-7, ratio: 2.5 },
  };
  for (const [name, def] of [['1 mm of 1 M KCl', asym], ['a symmetric 1 µm cell', sym]]) {
    const d = new Device(def), s = d.solve(), g = gridItem(d, s);
    assert.ok(s.converged);
    assert.equal(g.ok, false, `${name}: ${g.summary}`);
    assert.match(g.summary, /contact's D|charge/, name);
  }
});

// A port feeding a fast bulk A = B in strictly neutral KCl, on 8 cells: its flux of A was 12.5×
// the fine grid's, and check() said ok throughout, grid ok at 5.6e-16 (the contacts' fluxes), while
// balance's own ledger listed the port at 1.2e5 (Neon's second review, R2 2.2). Ports' fluxes and
// reactions' totals are compared now.
test("check()'s grid item compares a port's fluxes and the reactions' totals (Neon R2 2.2)", () => {
  const def = {
    species: [{ name: 'A', z: 0, cRef: 1000 }, { name: 'B', z: 0, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr: 0, species: { A: { D: 1e-9, mu0: 0 }, B: { D: 1e-9, mu0: 0 }, 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } } },
    bulkReactions: [{ equation: 'A = B', kf: { aq: 1e7 } }],
    regions: [{ name: 'el', material: 'aq', length: 1e-4, c0: { A: 10, B: 10, 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0.01, bath: { c: { 'K+': 100, 'Cl-': 100 } } }, right: { V: 0, bath: { c: { 'K+': 100, 'Cl-': 100 } } } },
    ports: [{ name: 'feed', region: 'el', from: 0.495e-4, to: 0.505e-4, species: { A: { type: 'equilibrium', mu: 0 }, B: { type: 'equilibrium', mu: -10000 } } }],
    grid: { hmax: 1.25e-5, minCells: 8 },
  };
  const d = new Device(def), s = d.solve(), g = gridItem(d, s);
  assert.equal(g.ok, false, g.summary);
  assert.match(g.summary, /feed|A = B/);
});

