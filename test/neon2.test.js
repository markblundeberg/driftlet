// Neon's second review (20 Nov 2026): impedance, check(), and its own pick.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';

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
