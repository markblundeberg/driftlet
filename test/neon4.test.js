// Neon's fourth review (22 Nov 2026): the one owner of strictly neutral charge (src/neutral.js) and
// its producers, the carry and the transient's cold start. Each device is the reviewer's
// (~/neon-review-4/review/r4-1.*.mjs at 0de29ae), each assertion what a physicist would ask of it.
// In a device whose every non-metal region is strictly neutral, the net charge is zero at all times
// (every capacitor's two plates are in it), so the terminals' currents sum to zero: a start whose
// charges disagree with its state breaks that sum at the first step, by the disagreement over dt.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY } from '../src/index.js';
import { check } from '../src/kit.js';
import { equilibrateFaces } from '../src/start.js';
import { gaussRows, strictlyNeutralNodes } from '../src/neutral.js';
import { gaussOff, startCharges } from '../bench/combos/charge.js';

const VT = (8.314462618 * 298.15) / FARADAY;
const nacl = (c = 100) => ({ 'Na+': c, 'Cl-': c });
const ions = [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }];
const aq = { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } };
const blocked = { 'Na+': 'blocked', 'Cl-': 'blocked' };
// Kirchhoff's sum over a solution's terminals, against their gross current.
const kirchhoff = (sol) => {
  const I = Object.values(sol.terminals).map((t) => t.current);
  return { sum: I.reduce((a, b) => a + b, 0), gross: I.reduce((a, b) => a + Math.abs(b), 0) };
};
const kirchhoffHolds = (sol, what, rel = 1e-6) => {
  const { sum, gross } = kirchhoff(sol);
  assert.ok(Math.abs(sum) <= rel * gross, `${what}: the terminals' currents sum to ${sum.toExponential(3)} A/m² of ${gross.toExponential(3)}`);
};
const amounts = (sol) => Object.fromEntries(sol.conservation.map((x) => [x.species, x.amount]));

// A capacitive contact is a gate over its end box: unchanged, a regrid moves its charge with the
// box's φ̂ (the end half box halves), and that is its start; changed (its C doubled), its charge
// can't jump, the box's φ̂ does. (Its start was kept from before a regrid, the first 1 ns step's
// currents summing to −1.8e4 A/m² of 1.2 A/m² of charging, unwarned; its C doubled in a sealed
// cell moved 1e-4 of the sealed ions and settled 50% off the fresh device's charge, unwarned;
// Neon's fourth review, R4 1.1.)
test('a capacitive contact beside a strictly neutral end box: its start charge through a regrid and its C changed (Neon R4 1.1)', () => {
  const open = {
    species: ions, materials: { aq },
    regions: [{ name: 'el', material: 'aq', length: 1e-5, c0: nacl() }],
    contacts: { left: { V: 0.2, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, species: blocked }, right: { V: 0, bath: { c: nacl() } } },
    grid: { hmax: 1.25e-6, minCells: 8 },
  };
  const a = new Device(open);
  a.solve();
  a.set({ contacts: { left: { V: 0.21 } } });
  a.advance(1e-5);
  a.set({ grid: { hmax: 6.25e-7, minCells: 16 } });
  assert.ok(gaussOff(a).worst < 1e-12, gaussOff(a).at);
  assert.ok(startCharges(a).worst < 1e-12, startCharges(a).at);
  const s = a.step(1e-9);
  assert.deepEqual(s.warnings, []);
  kirchhoffHolds(s, 'the first step after a regrid');
  assert.ok(Math.abs(s.terminals.right.current) < 2, `the charging current, ${s.terminals.right.current} A/m²`);

  const sealed = (CL) => ({
    species: ions, materials: { aq },
    regions: [{ name: 'el', material: 'aq', length: 1e-5, c0: nacl() }],
    contacts: {
      left: { V: 0.2, phi: { type: 'capacitive', C: CL }, zeroCharge: 0, species: blocked },
      right: { V: 0, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, species: blocked },
    },
    grid: { hmax: 1.25e-6, minCells: 8 },
  });
  const b = new Device(sealed(0.2));
  b.solve();
  b.step(1e-9);
  const before = amounts(b.solution()), Dbefore = b.solution().contacts.left.D;
  b.set({ contacts: { left: { phi: { type: 'capacitive', C: 0.4 } } } });
  const after = b.solution();
  for (const sp of ['Na+', 'Cl-']) assert.ok(Math.abs(amounts(after)[sp] / before[sp] - 1) < 1e-13, `${sp} moved by ${amounts(after)[sp] / before[sp] - 1}`);
  assert.ok(Math.abs(after.contacts.left.D / Dbefore - 1) < 1e-10, `the plate's charge jumped: ${Dbefore} → ${after.contacts.left.D}`);
  assert.ok(gaussOff(b).worst < 1e-12, gaussOff(b).at);
  const first = b.step(1e-9);
  kirchhoffHolds(first, 'the first step after C doubled');
  const end = b.advance(1e-2), fresh = new Device(sealed(0.4)).solve();
  assert.deepEqual(end.warnings, []);
  for (const side of ['left', 'right']) assert.ok(Math.abs(end.contacts[side].D / fresh.contacts[side].D - 1) < 1e-4, `${side}: ${end.contacts[side].D} where a fresh device holds ${fresh.contacts[side].D}`);
});

// A gate over the box beside a capacitive Pt | solution face (a porous electrode's window from
// x = 0) is a gate like any other: a gate added beside it, or a regrid, keeps the gates' charges
// and the first step Kirchhoff. (Its group was dropped, every box of it not a chargeNode: −4.9e4
// A/m² from nowhere at the first step, unwarned; Neon's fourth review, R4 1.2.)
test('a gate whose window reaches a capacitive face keeps its start through a gate added and a regrid (Neon R4 1.2)', () => {
  const g1 = { name: 'g1', region: 'el', from: 0, to: 0.6e-6, area: 1e6, V: 0.3, capacitance: { C: 1, zeroCharge: 0 } };
  const g2 = { name: 'g2', region: 'el', from: 0.2e-6, to: 0.6e-6, area: 1e6, V: 0.1, capacitance: { C: 1, zeroCharge: 0 } };
  const def = (ports) => ({
    species: [{ name: 'e-', z: -1, cRef: 1000 }, ...ions],
    materials: { Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } }, aq },
    regions: [{ name: 'm', material: 'Pt', length: 1e-7 }, { name: 'el', material: 'aq', length: 1e-6, c0: nacl() }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }],
    contacts: { left: { V: 0, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, right: { V: 0, bath: { c: nacl() } } },
    ports,
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  const a = new Device(def([g1]));
  a.solve();
  a.step(1e-9);
  a.set({ ports: [g1, g2] });
  assert.ok(startCharges(a).worst < 1e-12, startCharges(a).at);
  const sa = a.step(1e-9);
  assert.deepEqual(sa.warnings, []);
  kirchhoffHolds(sa, 'a gate added');

  const b = new Device(def([g1]));
  b.solve();
  b.set({ ports: [{ ...g1, V: 0.35 }] });
  b.advance(1e-6);
  b.step(1e-9);
  b.set({ grid: { hmax: 6.25e-8, minCells: 16 } });
  assert.ok(gaussOff(b).worst < 1e-12, gaussOff(b).at);
  assert.ok(startCharges(b).worst < 1e-12, startCharges(b).at);
  const sb = b.step(1e-9);
  assert.deepEqual(sb.warnings, []);
  kirchhoffHolds(sb, 'a regrid');
});

// A gate's hardware changed: a window narrowed keeps the charge over what's left of it; a change
// with a regrid in one set() keeps the group's total; a large floating gate beside a changed one
// keeps its charge, solved with the jump rather than after it. (−8.6e7 A/m², −4.1e4 A/m², each
// warned; and the floating gate left at −1.1e-3 C/m², 89% of its voltage's jump undone after 60
// rounds; Neon's fourth review, R4 1.3.)
test("a gate's hardware: a window narrowed, a change with a regrid, a floating gate beside it (Neon R4 1.3)", () => {
  const film = (ports, grid = { hmax: 1.25e-7, minCells: 8 }) => ({
    species: ions, materials: { aq },
    regions: [{ name: 'el', material: 'aq', length: 1e-6, c0: nacl() }],
    contacts: { left: { phi: 'neutral', species: blocked }, right: { V: 0, bath: { c: nacl() } } },
    ports, grid,
  });
  {
    const g = (to) => ({ name: 'g1', region: 'el', from: 0, to, area: 1e6, V: 0.3, capacitance: { C: 1, zeroCharge: 0 } });
    const d = new Device(film([g(0.6e-6)]));
    d.solve();
    d.step(1e-9);
    d.set({ ports: [g(0.3e-6)] });
    assert.ok(gaussOff(d).worst < 1e-12, gaussOff(d).at);
    assert.ok(startCharges(d).worst < 1e-12, startCharges(d).at);
    const s = d.step(1e-9);
    assert.deepEqual(s.warnings, []);
    kirchhoffHolds(s, 'a window narrowed');
  }
  {
    const g1 = (C) => ({ name: 'g1', region: 'el', from: 0.2e-6, to: 0.6e-6, area: 1e6, V: 0.3, capacitance: { C, zeroCharge: 0 } });
    const g2 = (V) => ({ name: 'g2', region: 'el', from: 0.4e-6, to: 0.8e-6, area: 1e6, V, capacitance: { C: 1, zeroCharge: 0 } });
    const d = new Device(film([g1(1), g2(0.1)]));
    d.solve();
    d.set({ ports: [g1(1), g2(0.4)] });
    d.advance(1e-6);
    d.set({ ports: [g1(2), g2(0.4)], grid: { hmax: 6.25e-8, minCells: 16 } });
    assert.ok(gaussOff(d).worst < 1e-12, gaussOff(d).at);
    assert.ok(startCharges(d).worst < 1e-12, startCharges(d).at);
    const s = d.step(1e-9);
    assert.deepEqual(s.warnings, []);
    kirchhoffHolds(s, 'a C doubled with a regrid');
  }
  {
    const g1 = (C) => ({ name: 'g1', region: 'el', from: 0.2e-6, to: 0.6e-6, area: 1e6, V: 0.3, capacitance: { C, zeroCharge: 0 } });
    const g2 = { name: 'g2', region: 'el', from: 0.2e-6, to: 0.6e-6, area: 1e6, I: 0, capacitance: { C: 10, zeroCharge: 0 } };
    const d = new Device(film([g1(0.01), g2]));
    d.solve();
    d.step(1e-9);
    const q1 = d.solution().ports[0].charge, q2 = d.solution().ports[1].charge, V2 = d.solution().terminals.g2.V;
    d.set({ ports: [g1(0.02), g2] });
    const after = d.solution();
    // (g1 keeps its charge at twice its C, so φ̂ under both jumps by half its drop, and the
    // floating g2, holding its charge, by as much)
    assert.ok(Math.abs(after.ports[1].charge - q2) < 1e-9 * Math.abs(q1), `the floating gate's charge, ${q2} → ${after.ports[1].charge}`);
    // (g1's drop q/(C a L) at 0.01 F/m², over a window of 0.4 µm at 1e6 m⁻¹, then at 0.02)
    const jump = q1 / (0.01 * 1e6 * 0.4e-6) - q1 / (0.02 * 1e6 * 0.4e-6);
    assert.ok(Math.abs(after.terminals.g2.V - V2 - jump) < 1e-9, `the floating gate's voltage moved ${after.terminals.g2.V - V2} V where ${jump} V is the jump`);
    assert.ok(gaussOff(d).worst < 1e-12, gaussOff(d).at);
  }
});

// A diffuse layer's side is 0 (left of its face) or 1: its edge box's φ̂ doesn't move the face's D,
// whose own row ties it to the layer. (Compared with 'left' and 'right', every edge box beside a
// layer took the face's A C V_T as a slope it doesn't have: a gate removed over it moved φ̂ there
// 5 V in 2,400 blind jumps; Neon's fourth review, R4 1.4.) A sheet charge on such a face: its
// counter-charge at t = 0⁺ is in the boxes beside it (the layers start uncharged, as c0 has them,
// and charge later), as beside a face without layers. (The phantom slope put it there; with no
// slope, nothing did, and the first 1 ns step's currents summed to −1e7 A/m², unwarned.)
test("a diffuse layer's edge box has no slope from its face (Neon R4 1.4)", () => {
  const d = new Device({
    species: [{ name: 'e-', z: -1, cRef: 1000 }, ...ions],
    materials: { Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } }, aq },
    regions: [{ name: 'm', material: 'Pt', length: 1e-7 }, { name: 'el', material: 'aq', length: 1e-6, c0: nacl() }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78 } }],
    contacts: { left: { V: 0.3, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, right: { V: 0, bath: { c: nacl() } } },
    ports: [{ name: 'g1', region: 'el', from: 0, to: 0.5e-6, area: 1e6, V: 0.2, capacitance: { C: 1, zeroCharge: 0 } }],
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  d.solve();
  d.step(1e-9);
  const gE = d.solver.model.grid.regionStart[1], phi0 = d.solution().phi[gE];
  d.set({ ports: [] });
  const nodes = strictlyNeutralNodes(d.solver);
  assert.equal(gaussRows(d.solver, nodes).find((r, a) => nodes[a].g === gE).slope, 0);
  assert.ok(Math.abs(d.solution().phi[gE] - phi0) < 0.1, `φ at the edge box moved ${d.solution().phi[gE] - phi0} V`);
  assert.ok(gaussOff(d).worst < 1e-12, gaussOff(d).at);

  const e = new Device({
    species: ions, materials: { aq },
    regions: [{ name: 'a', material: 'aq', length: 1e-6, c0: nacl() }, { name: 'b', material: 'aq', length: 1e-6, c0: nacl() }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.01 }, sheetCharge: 0.01, diffuse: { epsr: 78 }, species: { 'Na+': { type: 'permeability', P: 1e-6 }, 'Cl-': { type: 'permeability', P: 1e-6 } } }],
    contacts: { left: { V: 0, bath: { c: nacl() } }, right: { V: 0, bath: { c: nacl() } } },
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  equilibrateFaces(e.solver);
  assert.ok(gaussOff(e).worst < 1e-12, gaussOff(e).at);
  const se = e.step(1e-9);
  assert.deepEqual(se.warnings, []);
  kirchhoffHolds(se, 'a sheet charge between two layers, the first step');
});

// A supporting electrolyte (1 M NaCl, closed) with an electroactive trace (H⁺ open to both ends,
// 0.5% of the screening; or Fe³⁺ at 2.5 mM, z² c 1.1% of it, under a gate), regridded 4 → 16
// cells: the closed salt keeps its amounts and every box is neutral. (A switch at 1% of the
// screening between holding φ̂ and moving it left the salt 1.2e-7 off below it, a fixed point the
// trace's own regridded charge contradicted, and 8.3e-8 off just above it, the loop converging at
// 1.1% a round; check() failed conservation; Neon's fourth review, R4 1.5.)
test('a trace open ion in closed supporting salt, either side of 1% of its screening, regridded (Neon R4 1.5)', () => {
  const sp3 = (name, z, D) => ({ species: [...ions, { name, z, cRef: 1000 }], materials: { aq: { epsr: 0, species: { ...aq.species, [name]: { D, mu0: 0 } } } } });
  const a = new Device({
    ...sp3('H+', 1, 9.3e-9),
    regions: [{ name: 'el', material: 'aq', length: 1e-6, c0: { 'Na+': 1000, 'Cl-': 1005, 'H+': 5 } }],
    contacts: {
      left: { V: 0, terminal: 'H+', phi: 'neutral', species: { 'H+': { type: 'equilibrium', offset: VT * Math.log(15 / 1000) }, ...blocked } },
      right: { V: 0, terminal: 'H+', phi: 'neutral', species: { 'H+': { type: 'equilibrium', offset: VT * Math.log(0.01 / 1000) }, ...blocked } },
    },
    grid: { hmax: 2.5e-7, minCells: 4 },
  });
  a.solve();
  const b0 = amounts(a.step(1e-9));
  const gate = (V) => ({ name: 'gate', region: 'el', from: 0.25e-6, to: 0.75e-6, area: 1e6, V, capacitance: { C: 1, zeroCharge: 0 } });
  const b = new Device({
    ...sp3('Fe3+', 3, 6e-10),
    regions: [{ name: 'el', material: 'aq', length: 1e-6, c0: { 'Na+': 1000, 'Cl-': 1007.5, 'Fe3+': 2.5 } }],
    contacts: {
      left: { V: 0, terminal: 'Fe3+', phi: 'bulk', species: { 'Fe3+': { type: 'equilibrium', offset: (VT / 3) * Math.log(2.5 / 1000) }, ...blocked } },
      right: { phi: 'neutral', species: { ...blocked, 'Fe3+': 'blocked' } },
    },
    ports: [gate(0.3)],
    grid: { hmax: 2.5e-7, minCells: 4 },
  });
  b.solve();
  b.set({ ports: [gate(0.35)] });
  const c0 = amounts(b.advance(1e-6));
  for (const [d, was, name] of [[a, b0, 'H⁺ at 0.5%'], [b, c0, 'Fe³⁺ at 1.1%']]) {
    d.set({ grid: { hmax: 6.25e-8, minCells: 16 } });
    const now = amounts(d.solution());
    for (const sp of ['Na+', 'Cl-']) assert.ok(Math.abs(now[sp] / was[sp] - 1) < 1e-13, `${name}: ${sp} moved by ${(now[sp] / was[sp] - 1).toExponential(2)}`);
    assert.ok(gaussOff(d).worst < 1e-12, `${name}: ${gaussOff(d).at}`);
    const s = d.step(1e-9);
    assert.ok(check(d, s).ok, `${name}: ${check(d, s).text}`);
  }
});

// A transient's cold start: c0 is its initial condition, and the capacitors start uncharged, φ̂
// beside each where that puts it. (A sealed cell between plates at 0 and 1 V kept φ = 0 at its
// right end, its plate at −0.2 C/m², and the first 1 ns step's currents summed to −2e8 A/m², the
// right plate's charge settling with the wrong sign, check() passing; a floating gate over a held
// one, moved by a Donnan face's start, left it −6.6e7 A/m² off, warned; a floating port over a
// held one without c0 started 1e-5 off neutral; Neon's fourth review, R4 1.6.)
test('a cold start: capacitors uncharged, c0 kept, every box neutral (Neon R4 1.6)', () => {
  const cell = {
    species: ions, materials: { aq },
    regions: [{ name: 'el', material: 'aq', length: 1e-5, c0: nacl() }],
    contacts: {
      left: { V: 0, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, species: blocked },
      right: { V: 1, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, species: blocked },
    },
    grid: { hmax: 1.25e-6, minCells: 8 },
  };
  const a = new Device(cell);
  const s = a.step(1e-9);
  assert.deepEqual(s.warnings, []);
  kirchhoffHolds(s, 'the sealed cell, its first step');
  const b = new Device(cell), end = b.advance(1e-3);
  // (in series, 0.1 F/m² across 1 V: the right plate takes +0.1 C/m² in, the left gives it)
  assert.ok(Math.abs(b.solver.termCharge[1] + b.solver.termCharge[0]) < 1e-9, `the plates' charges, ${b.solver.termCharge[0]} and ${b.solver.termCharge[1]}`);
  assert.ok(b.solver.termCharge[1] > 0.0999 && b.solver.termCharge[1] < 0.1, `the right plate passed ${b.solver.termCharge[1]} C/m²`);
  assert.ok(Math.abs(end.contacts.left.D - end.contacts.right.D) < 1e-9);

  const donnan = new Device({
    species: ions, materials: { aq },
    regions: [
      { name: 'L', material: 'aq', length: 1e-6, c0: nacl() },
      { name: 'R', material: 'aq', length: 1e-6, fixedCharge: -100 * FARADAY, c0: { 'Na+': 200, 'Cl-': 100 } },
    ],
    interfaces: [{ phi: 'neutral', species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium' } }],
    contacts: { left: { V: 0, bath: { c: nacl() } }, right: { phi: 'neutral', species: blocked } },
    ports: [
      { name: 'flt', region: 'R', from: 0, to: 0.5e-6, area: 1e6, I: 0, capacitance: { C: 1, zeroCharge: 0 } },
      { name: 'held', region: 'R', from: 0.25e-6, to: 0.75e-6, area: 1e6, V: 0.2, capacitance: { C: 1, zeroCharge: 0 } },
    ],
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  const sd = donnan.step(1e-9);
  assert.deepEqual(sd.warnings, []);
  kirchhoffHolds(sd, 'the Donnan face with a floating gate over a held one');

  const plain = new Device({
    species: ions, materials: { aq },
    regions: [{ name: 'el', material: 'aq', length: 1e-6 }],
    contacts: { left: { phi: 'neutral', species: blocked }, right: { V: 0, bath: { c: nacl() } } },
    ports: [
      { name: 'held', region: 'el', from: 0, to: 0.5e-6, area: 1e6, V: 0.3, capacitance: { C: 1, zeroCharge: 0 } },
      { name: 'flt', region: 'el', from: 0.25e-6, to: 0.75e-6, area: 1e6, I: 0, capacitance: { C: 1, zeroCharge: 0 } },
    ],
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  equilibrateFaces(plain.solver);
  assert.ok(gaussOff(plain).worst < 1e-12, gaussOff(plain).at);
  const sp = plain.step(1e-9);
  assert.deepEqual(sp.warnings, []);
  kirchhoffHolds(sp, 'a floating port over a held one, no c0');
});

// A gate's C doubled over the box beside a capacitive Pt | solution face: φ̂ jumps there, and the
// face's D with it. The metal's sheet can't move in an instant any more than the gate's charge
// can, so both start from before the jump and the first step passes what moved between them as
// currents, which their terminals' running charges record: each terminal's charge passed less
// what it stores (the gate's charge, the Pt's sheet A·D) is what it was. (Shared out onto the
// gate's start instead, as an impulse at the set(), Kirchhoff held but 0.0115 C/m² passed between
// the gate and the Pt that neither terminal saw; Neon's fourth review, its closing reply.)
test("a capacitor's change beside a face to a metal: what moves between them passes as currents, recorded (Neon R4, closing)", () => {
  const g1 = (C) => ({ name: 'g1', region: 'el', from: 0, to: 0.6e-6, area: 1e6, V: 0.3, capacitance: { C, zeroCharge: 0 } });
  const def = (ports) => ({
    species: [{ name: 'e-', z: -1, cRef: 1000 }, ...ions],
    materials: { Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } }, aq },
    regions: [{ name: 'm', material: 'Pt', length: 1e-7 }, { name: 'el', material: 'aq', length: 1e-6, c0: nacl() }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }],
    contacts: { left: { V: 0, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, right: { V: 0, bath: { c: nacl() } } },
    ports,
    grid: { hmax: 1.25e-7, minCells: 8 },
  });
  // (each terminal's charge passed less what it stores: the Pt's sheet, σ = D on the left of the face)
  const ledger = (d) => {
    const s = d.solver, sol = d.solution(), P = (nm) => s.termCharge[s.terms.findIndex((t) => t.name === nm)];
    return { gate: P('g1') - sol.ports[0].charge, metal: P('left') - s.model.grid.area[0] * sol.interfaces[0].D, D: sol.interfaces[0].D };
  };
  for (const [what, then] of [['alone', null], ['then a drive changed', { contacts: { right: { V: 0.01 } } }], ['then its C again', { ports: [g1(3)] }]]) {
    const d = new Device(def([g1(1)]));
    d.solve();
    d.step(1e-9);
    const before = ledger(d);
    d.set({ ports: [g1(2)] });
    if (then) d.set(then);
    const moved = Math.abs(ledger(d).D - before.D);
    assert.ok(moved > 1e-3, `${what}: the face's D moved by ${moved}`);
    assert.ok(gaussOff(d).worst < 1e-12, gaussOff(d).at);
    assert.ok(startCharges(d).worst < 1e-12, `${what}: ${startCharges(d).at}`);
    const first = d.step(1e-9);
    assert.deepEqual(first.warnings, []);
    kirchhoffHolds(first, `${what}: the first step`);
    d.advance(1e-6);
    const after = ledger(d);
    for (const k of ['gate', 'metal']) assert.ok(Math.abs(after[k] - before[k]) < 1e-12 * moved * 1e3, `${what}: the ${k}'s charge passed less what it stores moved by ${(after[k] - before[k]).toExponential(3)} C/m² (the face's D moved ${moved.toExponential(3)})`);
  }
});
