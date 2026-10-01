import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, units } from '../src/index.js';

// A valid baseline: an n-type / p-type heterojunction between two ohmic contacts.
function hetero() {
  return {
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: {
      A: {
        epsr: 12.9,
        species: {
          'e-': { D: 2e-2, mu0: units.eV(-4.07), cRef: units.perCm3(4.7e17) },
          'h+': { D: 1e-3, mu0: units.eV(-(-4.07 - 1.42)), cRef: units.perCm3(9e18) },
        },
      },
      B: {
        epsr: 12.2,
        species: {
          'e-': { D: 1e-2, mu0: units.eV(-3.8), cRef: units.perCm3(5e17) },
          'h+': { D: 5e-4, mu0: units.eV(-(-3.8 - 1.8)), cRef: units.perCm3(1e19) },
        },
      },
    },
    regions: [
      { name: 'n', material: 'A', length: 1e-6, fixedCharge: 1e3 },
      { name: 'p', material: 'B', length: 1e-6, fixedCharge: -1e3 },
    ],
    interfaces: [{ step: { species: 'e-', value: units.eV(0.25) } }],
    contacts: {
      left: { species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk', terminal: 'e-' },
      right: { species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk', terminal: 'e-' },
    },
    grid: { hmin: 1e-10, hmax: 2e-8 },
  };
}

const throwsDevice = (def, pattern) =>
  assert.throws(() => new Device(def), (err) => err instanceof DeviceError && pattern.test(err.message));

test('a valid device normalises and builds its grid', () => {
  const dev = new Device(hetero());
  const m = dev.model;
  assert.equal(m.species.length, 2);
  assert.equal(m.regions.length, 2);
  assert.equal(m.interfaces.length, 1);
  assert.equal(dev.grid.nodeRegion[dev.grid.nNodes - 1], 1);
  assert.equal(m.contacts.left.terminal, 0);
  assert.equal(m.T, 298.15);
});

test('a standard-level step alignment becomes the matching interface dipole', () => {
  const def = hetero();
  const m = new Device(def).model;
  // Electron standard level μ° + zFφ must step by +0.25 eV across the face.
  const dmu0 = def.materials.B.species['e-'].mu0 - def.materials.A.species['e-'].mu0;
  const expected = (units.eV(0.25) - dmu0) / (-1 * FARADAY);
  assert.ok(Math.abs(m.interfaces[0].dipole - expected) < 1e-15);
  const stdStep = dmu0 + -1 * FARADAY * m.interfaces[0].dipole;
  assert.ok(Math.abs(stdStep - units.eV(0.25)) < 1e-9);
});

test('an explicit dipole is taken as given', () => {
  const def = hetero();
  def.interfaces = [{ dipole: 0.3 }];
  assert.equal(new Device(def).model.interfaces[0].dipole, 0.3);
});

test('a heterointerface without an alignment is an error, never a vacuum-level default', () => {
  const def = hetero();
  delete def.interfaces;
  throwsDevice(def, /needs an alignment.*no Anderson or Schottky–Mott/);
  def.interfaces = [{ sheetCharge: 0.01 }];
  throwsDevice(def, /needs an alignment/);
});

test('a same-material face needs no alignment but may take one', () => {
  const def = hetero();
  def.regions[1].material = 'A';
  delete def.interfaces;
  assert.equal(new Device(def).model.interfaces[0].dipole, 0);
  def.interfaces = [{ dipole: -0.05 }];
  assert.equal(new Device(def).model.interfaces[0].dipole, -0.05);
});

test('two alignments on one face is an error', () => {
  const def = hetero();
  def.interfaces = [{ dipole: 0.1, step: { species: 'e-', value: 0 } }];
  throwsDevice(def, /exactly one alignment, got dipole and step/);
});

test('step alignment needs a charged species present on both sides', () => {
  const def = hetero();
  def.species.push({ name: 'O2', z: 0 });
  def.materials.A.species.O2 = { D: 1e-9, mu0: 0, cRef: 1 };
  def.materials.B.species.O2 = { D: 1e-9, mu0: 0, cRef: 1 };
  def.interfaces = [{ step: { species: 'O2', value: 0 } }];
  throwsDevice(def, /neutral; alignment needs a charged species/);

  const def2 = hetero();
  delete def2.materials.B.species['h+'];
  def2.contacts.right.species = { 'e-': 'equilibrium' };
  def2.interfaces = [{ step: { species: 'h+', value: 0 } }];
  throwsDevice(def2, /must be present on both sides/);
});

test('unknown names and bad numbers name the offending path', () => {
  let def = hetero();
  def.materials.A.species.Na = { D: 1, mu0: 0, cRef: 1 };
  throwsDevice(def, /materials\.A\.species\.Na: unknown species 'Na'/);

  def = hetero();
  def.regions[0].material = 'C';
  throwsDevice(def, /regions\[0\]\.material: unknown material "C"/);

  def = hetero();
  def.species[0].z = 0.5;
  throwsDevice(def, /species\[0\]\.z .* must be an integer/);

  def = hetero();
  def.species.push({ name: 'e-', z: -1 });
  throwsDevice(def, /duplicate species 'e-'/);

  def = hetero();
  def.materials.A.species['e-'].D = -1;
  throwsDevice(def, /materials\.A\.species\.e-\.D must be a non-negative number/);

  def = hetero();
  def.regions[1].length = 0;
  throwsDevice(def, /regions\[1\]\.length must be a positive number/);

  def = hetero();
  def.contacts.left.species['e-'] = 'ohmic';
  throwsDevice(def, /contacts\.left\.species\.e-\.type must be one of/);
});

test('reference concentrations are never defaulted silently', () => {
  const def = hetero();
  delete def.materials.A.species['e-'].cRef;
  throwsDevice(def, /materials\.A\.species\.e-\.cRef: no reference concentration/);
  def.species[0].cRef = units.perCm3(4.7e17); // a species-level default is fine
  assert.doesNotThrow(() => new Device(def));
});

test('a species present in no material is an error', () => {
  const def = hetero();
  def.species.push({ name: 'Li+', z: 1, cRef: 1000 });
  throwsDevice(def, /'Li\+' is not present in any material/);
});

test('too many interface entries is an error', () => {
  const def = hetero();
  def.interfaces.push({ dipole: 0 });
  throwsDevice(def, /interfaces has 2 entries but there are only 1 faces/);
});

test('a floating island with no gate has no electrostatic anchor', () => {
  const def = hetero();
  def.contacts = {};
  throwsDevice(def, /no electrostatic anchor/);
  // A gate at either end anchors it.
  def.contacts = { right: { phi: { type: 'capacitive', C: 1e-3, zeroCharge: 0 } } };
  const m = new Device(def).model;
  assert.equal(m.contacts.left.phi.type, 'neutral');
  assert.ok(m.contacts.left.species.every((l) => l.type === 'blocked'));
});

test('capacitive links need a capacitance and a zero-charge alignment', () => {
  const def = hetero();
  def.contacts.left.phi = { type: 'capacitive', C: 1e-3 };
  throwsDevice(def, /contacts\.left\.phi\.zeroCharge .* zero-charge alignment/);
  def.contacts.left.phi = { type: 'capacitive', zeroCharge: 0 };
  throwsDevice(def, /contacts\.left\.phi\.C must be a positive number/);
});

test("a 'bulk' contact needs a connected charged species", () => {
  const def = hetero();
  def.contacts.left.species = {};
  delete def.contacts.left.terminal;
  throwsDevice(def, /'bulk' needs at least one connected charged species/);
});

test('only blocking is allowed for a species absent from the end material', () => {
  const def = hetero();
  delete def.materials.A.species['h+'];
  throwsDevice(def, /contacts\.left\.species\.h\+: 'h\+' is absent from the end material 'A'/);
  def.contacts.left.species['h+'] = 'blocked';
  assert.doesNotThrow(() => new Device(def));
});

test('a terminal must be a connected species', () => {
  const def = hetero();
  def.contacts.left.species['e-'] = 'blocked';
  throwsDevice(def, /contacts\.left\.terminal: 'e-' is blocked/);
});

test('definitions are plain data and survive structured cloning', () => {
  const def = structuredClone(hetero());
  assert.doesNotThrow(() => new Device(def));
});

