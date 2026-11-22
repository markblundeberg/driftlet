import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { check } from '../src/checks.js';

// The grid check's floors, after Neon's fourth review (flag 2.1).
const item = (report, name) => report.items.find((it) => it.name === name);
const diff = (report, re) => item(report, 'grid').details.diffs.find((d) => re.test(d.what));

// Flag 2.1: the floor under a current was 1e-16 of what any charged species could carry, a metal's
// carrier by its conductance included. A 1 µm sphere of Pt facing a 1 cm Pt counter electrode
// across 1 cm of electrolyte: the counter's 100 nm of Pt, 4π(1 cm)² across, set a floor of
// 304 nA under the electrode's 1.96 nA, and the current, 9.24% lower on a grid twice as fine
// (19.1% off the converged one), read 6.0e-4: "ok grid", check() ok. A metal's current is an
// unknown of its own (Ohm's law in mixed form), resolved as finely as what feeds it, so its
// conductance is no floor: with the electrolyte's (2.5e-19 A here) the change reads 9.2e-2.
const ume = {
  geometry: { type: 'spherical', r0: 0.9e-6 },
  species: [
    { name: 'e-', z: -1, cRef: 1000 },
    { name: 'A', z: 0, cRef: 1000 },
    { name: 'Red', z: 0, cRef: 1000 },
    { name: 'Ox', z: 1, cRef: 1000 },
    { name: 'Na+', z: 1, cRef: 1000 },
    { name: 'Cl-', z: -1, cRef: 1000 },
  ],
  materials: {
    Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } },
    aq: {
      epsr: 0,
      species: { A: { D: 1e-9, mu0: 0 }, Red: { D: 1e-9, mu0: 0 }, Ox: { D: 1e-9, mu0: 0 }, 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } },
    },
  },
  regions: [
    { name: 'we', material: 'Pt', length: 0.1e-6 },
    { name: 'el', material: 'aq', length: 1e-2 - 1e-6, c0: { A: 1, Red: 1, Ox: 1, 'Na+': 100, 'Cl-': 101 } },
    { name: 'ce', material: 'Pt', length: 1e-7 },
  ],
  bulkReactions: [{ equation: 'A = Red', kf: { aq: 1e3 } }],
  interfaces: [
    { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, reactions: [{ equation: 'Red = Ox + e-', k0: 1e-1, alpha: 0.5 }] },
    { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, reactions: [{ equation: 'Red = Ox + e-', k0: 1e-1, alpha: 0.5 }] },
  ],
  contacts: {
    left: { V: 0.3, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } },
    right: { V: 0, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } },
  },
  grid: { hmin: 2e-6, hmax: 2e-3, minCells: 10 },
};

test('grid check: a metal’s conductance sets no floor under a current it only passes on', () => {
  const dev = new Device(ume), sol = dev.solve(), report = check(dev, sol);
  assert.ok(sol.converged);
  const I = diff(report, /^current$/);
  assert.ok(Math.abs(I.a - 1.959e-9) < 1e-12, `the electrode's current, ${I.a} A`);
  // The change on refinement, against the current itself, not a 304 nA floor.
  const change = Math.abs(I.b - I.a) / Math.max(Math.abs(I.a), Math.abs(I.b));
  assert.ok(change > 0.05, `the current moved by ${change} on refinement`);
  assert.ok(Math.abs(I.rel - change) <= 1e-12, `current read at ${I.rel}, against a change of ${change}`);
  // And the faces' rates (each a metal's carrier's partner) against themselves too.
  const face = diff(report, /^face 0's reaction/);
  assert.ok(face && face.rel > 0.05, `face 0's rate read at ${face?.rel}`);
  assert.equal(item(report, 'grid').ok, false, item(report, 'grid').summary);
  assert.equal(report.ok, false);
  // The order check's currents are compared now as well (under the old floor they were skipped).
  assert.ok(item(report, 'order').details.diffs.some((d) => /current/.test(d.what)));
  assert.equal(item(report, 'order').ok, true, item(report, 'order').summary);
});

// The same electrodes at equilibrium (both at 0 V, one couple): nothing flows, and what the
// solves leave, 1e-30 A, is under the electrolyte's floor. Not a FAIL of round-off.
test('grid check: the same cell at equilibrium has nothing to compare in its currents', () => {
  const dev = new Device({ ...ume, contacts: { ...ume.contacts, left: { ...ume.contacts.left, V: 0 } } });
  const sol = dev.solve(), report = check(dev, sol);
  assert.ok(sol.converged);
  assert.notEqual(item(report, 'grid').ok, false, item(report, 'grid').summary);
  assert.ok(report.ok, report.text);
});

// A device of metals alone has no other scale: its floor is still its metals'.
test('grid check: a device of metals alone passes on its own floor', () => {
  const dev = new Device({
    species: [{ name: 'e-', z: -1 }],
    materials: { Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } }, Cu: { conductor: { species: 'e-', conductivity: 5.96e7 } } },
    regions: [
      { name: 'a', material: 'Pt', length: 1e-3 },
      { name: 'b', material: 'Cu', length: 1e-3 },
    ],
    contacts: { left: { V: 1e-6, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, right: { V: 0, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } } },
  });
  const sol = dev.solve(), report = check(dev, sol);
  assert.ok(sol.converged);
  // Ohm's law through both: 1 µV across 1 mm of Pt and 1 mm of Cu, per m².
  const R = 1e-3 / 9.4e6 + 1e-3 / 5.96e7;
  assert.ok(Math.abs(Math.abs(sol.current) - 1e-6 / R) <= 1e-9 * (1e-6 / R), `${sol.current} A/m² against ${1e-6 / R}`);
  assert.ok(report.ok, report.text);
});
