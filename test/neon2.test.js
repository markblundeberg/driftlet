// Neon's second review (20 Nov 2026): impedance, check(), and its own pick.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT } from '../src/index.js';

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
