import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, units } from '../src/index.js';
import { build, combine, layer, ohmic, bath, parseEquation, stoichiometry, half, SHE, level } from '../src/kit.js';

// The kit's stack: a device written left to right as it's drawn, built into the plain
// definition. It must give exactly the definition one would write by hand.

const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
const Si = {
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } } },
};

test('a pn diode from a stack is the definition written by hand', () => {
  const built = build({
    library: [Si],
    stack: [
      ohmic(0),
      layer('Si', 1e-6, { donors: units.perCm3(1e17) }),
      layer('Si', 1e-6, { acceptors: units.perCm3(1e16) }),
      ohmic(0.5),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e13 } }],
    grid: { hmin: 1e-9, hmax: 20e-9 },
  });
  const contact = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
  const byHand = {
    species: Si.species,
    materials: Si.materials,
    regions: [
      { material: 'Si', length: 1e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'Si', length: 1e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
    ],
    contacts: { left: contact(0), right: contact(0.5) },
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e13 } }],
    grid: { hmin: 1e-9, hmax: 20e-9 },
  };
  assert.deepEqual(built, byHand);
  // The equation is the stoichiometry written out.
  const nu = { ...byHand, bulkReactions: [{ nu: { 'e-': -1, 'h+': -1 }, kf: { Si: 1e13 } }] };
  assert.deepEqual(new Device(built).model.reactions, new Device(nu).model.reactions);
  assert.ok(new Device(built).solve().converged);
});

// Silver electrodes in silver nitrate, electrons and the plating reaction at each face.
const ions = {
  species: [
    { name: 'Ag+', z: 1, cRef: 1000 },
    { name: 'NO3-', z: -1, cRef: 1000 },
    { name: 'e-', z: -1 },
  ],
  materials: {
    water: { epsr: 78.5, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } },
    Ag: { conductor: { species: 'e-', conductivity: 6.3e7 } },
  },
};
const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 });
const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 };
const silverCell = (V) =>
  build({
    library: [ions],
    stack: [
      ohmic(0, ['e-']),
      layer('Ag', 1e-6),
      { ...stern, reactions: [{ ...silver, k0: 1e-3, alpha: 0.3 }] },
      layer('water', 10e-6, { c0: { 'Ag+': 10, 'NO3-': 10 } }),
      { ...stern, reactions: [{ ...silver, k0: 1e-3, alpha: 0.3 }] },
      layer('Ag', 1e-6),
      ohmic(V, ['e-']),
    ],
    grid: { hmin: 0.1e-9, hmax: 100e-9, ratio: 1.15 },
  });

// A face reaction's participants by name and side, from the device model.
const parts = (dev, f, k = 0) =>
  Object.fromEntries(dev.model.interfaces[f].reactions[k].part.map((p) => [`${dev.model.species[p.i].name}@${p.side ? 'right' : 'left'}`, p.nu]));

test('a face reaction written as an equation: each participant goes to the side that holds it', () => {
  const def = silverCell(0.05);
  assert.deepEqual(def.interfaces[0].reactions[0], { equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-3, alpha: 0.3 });
  assert.deepEqual(def.contacts.left, { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const dev = new Device(def);
  assert.deepEqual(parts(dev, 0), { 'Ag+@right': -1, 'e-@left': -1 });
  assert.deepEqual(parts(dev, 1), { 'Ag+@left': -1, 'e-@right': -1 });
  // The same as the reaction written out by sides.
  const bySides = { ...def, interfaces: [{ ...stern, reactions: [{ left: { 'e-': -1, 'Ag(s)': 1 }, right: { 'Ag+': -1 }, fixed: { 'Ag(s)': 0 }, k0: 1e-3, alpha: 0.3 }] }, def.interfaces[1]] };
  assert.deepEqual(new Device(bySides).model.interfaces[0].reactions, dev.model.interfaces[0].reactions);
  const sol = dev.solve();
  assert.ok(sol.converged && sol.current < 0, 'silver dissolves at the right electrode (the higher electron voltage) and plates at the left');
  // Both forms at once is an error, not a choice.
  assert.throws(() => new Device({ ...def, interfaces: [{ ...stern, reactions: [{ equation: 'Ag+ + e- = Ag(s)', left: { 'e-': -1 }, fixed: { 'Ag(s)': 0 }, k0: 1 }] }, def.interfaces[1]] }), /not both/);
});

test("at equilibrium, each electrode's Fermi level sits at the silver couple's redox level", () => {
  const sol = new Device(silverCell(0)).solve();
  assert.ok(sol.converged);
  const V = level(sol, silver);
  const inWater = [...sol.region].map((r, g) => (r === 1 ? g : -1)).filter((g) => g >= 0);
  // Ag⁺ + e⁻ ⇌ Ag with μ(Ag) = 0: the level is V_Ag⁺ wherever there's Ag⁺, NaN in the metal.
  for (const g of inWater) assert.ok(Math.abs(V[g] - sol.V['Ag+'][g]) < 1e-12);
  assert.ok(Number.isNaN(V[0]));
  for (const g of inWater) assert.ok(Math.abs(V[g] - sol.V['e-'][0]) < 1e-9, `node ${g}: ${V[g]} vs the metal's ${sol.V['e-'][0]}`);
});

test("the SHE's standard level is φ + μ°_H⁺/F, and its level at pH 3 is 3·ln10·RT/F below that", () => {
  const def = build({
    species: [
      { name: 'H+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 0, species: { 'H+': { D: 9.3e-9, mu0: 0 }, 'Cl-': { D: 2.0e-9, mu0: -131.2e3 } } } },
    stack: [bath({ 'H+': 1, 'Cl-': 1 }, 'Cl-'), layer('water', 1e-6), bath({ 'H+': 1, 'Cl-': 1 }, 'Cl-')],
  });
  const sol = new Device(def).solve();
  const std = level(sol, SHE, { standard: true }), actual = level(sol, SHE);
  const VT = (8.314462618 * 298.15) / FARADAY;
  for (let g = 0; g < sol.x.length; g++) {
    assert.ok(Math.abs(std[g] - sol.phi[g]) < 1e-12);
    assert.ok(Math.abs(std[g] - actual[g] - 3 * Math.LN10 * VT) < 1e-6, `${std[g] - actual[g]}`);
  }
});

test('a species on both sides of a face is labelled with its side or its material', () => {
  const lib = {
    species: [
      { name: 'Li+', z: 1, cRef: 1000 },
      { name: 'PF6-', z: -1, cRef: 1000 },
    ],
    materials: {
      electrolyte: { epsr: 0, species: { 'Li+': { D: 3e-10, mu0: 0 }, 'PF6-': { D: 3e-10, mu0: 0 } } },
      gel: { epsr: 0, species: { 'Li+': { D: 1e-10, mu0: 0 }, 'PF6-': { D: 1e-10, mu0: 0 } } },
    },
  };
  const end = { terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' };
  const def = (equation) =>
    build({
      library: [lib],
      stack: [{ V: 0, ...end }, layer('electrolyte', 1e-6, { c0: { 'PF6-': 1000 } }), { species: { 'Li+': 'blocked' }, reactions: [{ equation, k0: 1e-3 }] }, layer('gel', 1e-6, { c0: { 'PF6-': 1000 } }), { V: 0.01, ...end }],
    });
  const want = { 'Li+@left': -1, 'Li+@right': 1 };
  assert.deepEqual(parts(new Device(def('Li+(left) = Li+(right)')), 0), want);
  assert.deepEqual(parts(new Device(def('Li+(electrolyte) = Li+(gel)')), 0), want);
  assert.throws(() => new Device(def('Li+ = Li+')), (e) => e instanceof DeviceError && /Li\+\(left\)/.test(e.message));
  assert.throws(() => new Device(def('Li+(electrolyte) = Li+(water)')), (e) => e instanceof DeviceError && /neither/.test(e.message));
  assert.ok(new Device(def('Li+(left) = Li+(right)')).solve().converged);
});

test('equations: coefficients, empty sides, and what they refuse', () => {
  assert.deepEqual(parseEquation('2 H+ + 2 e- = H2'), [
    { name: 'H+', nu: -2 },
    { name: 'e-', nu: -2 },
    { name: 'H2', nu: 1 },
  ]);
  assert.deepEqual(stoichiometry('H2O ⇌ H+ + OH-'), { H2O: -1, 'H+': 1, 'OH-': 1 });
  assert.deepEqual(stoichiometry('∅ = e- + h+'), { 'e-': 1, 'h+': 1 });
  // A coefficient needs its space: names may start with digits.
  assert.deepEqual(parseEquation('3He+ + e- = 3He'), [
    { name: '3He+', nu: -1 },
    { name: 'e-', nu: -1 },
    { name: '3He', nu: 1 },
  ]);
  for (const bad of ['Ag+ + e-=Ag', 'Ag+ + e- = Ag = X', 'Ag+ +e- = Ag', '0 = 0', '0 Ag = Ag+']) {
    assert.throws(() => parseEquation(bad), DeviceError, bad);
  }
  assert.throws(() => stoichiometry('A + B = A + C'), /both sides/);
  // Equations are part of the definition itself, so set() takes them too; a forgotten space
  // after a coefficient is pointed out.
  const dev = new Device(build({ library: [Si], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1 }), ohmic(0)] }));
  dev.set({ bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e5 } }] });
  assert.ok(dev.solve().converged);
  assert.throws(() => dev.set({ bulkReactions: [{ equation: '2e- + 2 h+ = 0', kf: { Si: 1 } }] }), /if 2 is a coefficient, put a space after it/);
  assert.throws(() => half('Fe3+ = Fe2+'), /electrons/);
  assert.throws(() => half('Ag+ + e- = Ag', { Au: 0 }), /isn't in/);
});

test('stack and library mistakes say where', () => {
  const ok = [ohmic(0), layer('Si', 1e-6), ohmic(0)];
  assert.throws(() => build({ library: [Si], stack: [ohmic(0), { dipole: 0 }, layer('Si', 1e-6), ohmic(0)] }), /stack\[1\].*before the first layer/);
  assert.throws(() => build({ library: [Si], stack: [ohmic(0), layer('Si', 1e-6), {}, {}, layer('Si', 1e-6), ohmic(0)] }), /stack\[3\].*two faces/);
  assert.throws(() => build({ library: [Si], stack: [layer('Si', 1e-6), layer('Si', 1e-6), ohmic(0)] }), /starts and ends with a contact/);
  assert.throws(() => build({ library: [Si], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1, fixedCharge: 2 }), ohmic(0)] }), /not both/);
  assert.throws(() => build({ library: [Si], stack: ok, interfaces: [] }), /build\.interfaces: unknown field/);
  const other = { species: [{ name: 'e-', z: 1 }] };
  assert.throws(() => combine(Si, other), /species 'e-' is given twice, differently/);
  assert.deepEqual(combine(Si, Si), { species: Si.species, materials: Si.materials });
  // Nested lists flatten, so a helper can return several stack items.
  assert.deepEqual(build({ library: [Si], stack: [[ohmic(0), layer('Si', 1e-6)], ohmic(0)] }), build({ library: [Si], stack: ok }));
});
