import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, units } from '../src/index.js';

const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const throwsDevice = (def, pattern) =>
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && pattern.test(e.message), String(pattern));

test('metal region: ohmic, with a Fermi level and no φ or concentration', () => {
  const sigma = 1e6, L = 1e-6, V = 1e-3;
  const sol = new Device({
    species: [{ name: 'e-', z: -1 }],
    materials: { Cu: { metal: { species: 'e-', conductivity: sigma } } },
    regions: [{ material: 'Cu', length: L }],
    contacts: { left: collector(0), right: collector(V) },
  }).solve();
  assert.ok(sol.converged);
  assert.ok(Math.abs(sol.current / ((-sigma * V) / L) - 1) < 1e-12);
  const mid = sol.x.length >> 1;
  assert.ok(Number.isNaN(sol.phi[mid]) && Number.isNaN(sol.c['e-'][mid]) && Number.isNaN(sol.Vstd['e-'][mid]));
  assert.ok(Math.abs(sol.V['e-'].at(-1) - V) < 1e-15);
});

// n-Si with a gold region on the left, against the same metal as a contact.
const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19), ND = units.perCm3(1e16);
const W = 4.75, chi = 4.05;
const silicon = () => ({
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } },
    Au: { metal: { species: 'e-', conductivity: 4e7 } },
  },
  contacts: { right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' } },
  grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
});

test('metal | semiconductor face: a Schottky barrier from work function and electron affinity, as at a contact', () => {
  const asContact = silicon();
  asContact.regions = [{ material: 'Si', length: 2e-6, fixedCharge: ND * FARADAY }];
  asContact.contacts.left = { ...collector(0), phi: { type: 'capacitive', C: 10, vacuum: { outside: W, inside: { anchor: 'e-', offset: chi } } } };
  const asRegion = silicon();
  asRegion.regions = [{ material: 'Au', length: 50e-9 }, { material: 'Si', length: 2e-6, fixedCharge: ND * FARADAY }];
  asRegion.interfaces = [{ phi: { type: 'capacitive', C: 10 }, vacuum: { left: { anchor: 'fermi', offset: W }, right: { anchor: 'e-', offset: chi } } }];
  asRegion.contacts.left = collector(0);
  const a = new Device(asContact), b = new Device(asRegion);
  const g = b.grid.regionStart[1];
  for (const V of [0, 0.1, -0.2]) {
    a.set({ contacts: { right: { V } } });
    b.set({ contacts: { right: { V } } });
    const sa = a.solve(), sb = b.solve();
    assert.ok(sa.converged && sb.converged, `V=${V}`);
    assert.ok(Math.abs(sb.c['e-'][g] / sa.c['e-'][0] - 1) < 1e-9, 'surface electron density');
    // The metal's surface charge is the face's displacement, the same as the contact's.
    assert.ok(Math.abs(sb.interfaces[0].D / sa.contacts.left.D - 1) < 1e-9, 'metal charge');
    // Same current, plus the gold's negligible resistance.
    assert.ok(Math.abs(sb.current - sa.current) <= 1e-9 * Math.abs(sa.current) + 1e-12, `V=${V}: ${sb.current} vs ${sa.current}`);
  }
});

test('metal gate over an oxide region matches a capacitive gate contact (MOS)', () => {
  // Gate | oxide (no mobile species, ε_r = 3.9, 5 nm) | p-Si. As a contact, the oxide and the
  // gate's face capacitance combine in series.
  const tox = 5e-9, Cox = (3.9 * 8.8541878188e-12) / tox, Cface = 50, zc = -0.9;
  const NA = units.perCm3(1e17);
  const asContact = silicon();
  asContact.regions = [{ material: 'Si', length: 0.5e-6, fixedCharge: -NA * FARADAY }];
  asContact.contacts.left = { V: 0, phi: { type: 'capacitive', C: 1 / (1 / Cox + 1 / Cface), zeroCharge: zc } };
  const asRegion = silicon();
  asRegion.materials.SiO2 = { epsr: 3.9, species: {} };
  asRegion.regions = [
    { material: 'Au', length: 20e-9 },
    { material: 'SiO2', length: tox, grid: { minCells: 4 } },
    { material: 'Si', length: 0.5e-6, fixedCharge: -NA * FARADAY },
  ];
  asRegion.interfaces = [{ phi: { type: 'capacitive', C: Cface }, zeroCharge: zc }, { dipole: 0 }];
  asRegion.contacts.left = collector(0);
  const a = new Device(asContact), b = new Device(asRegion);
  const gSi = b.grid.regionStart[2];
  for (const V of [-1, 0, 0.5, 1.5]) {
    a.set({ contacts: { left: { V } } });
    b.set({ contacts: { left: { V } } });
    const sa = a.solve(), sb = b.solve();
    assert.ok(sa.converged && sb.converged, `V=${V}`);
    assert.ok(Math.abs(sb.interfaces[0].D - sa.contacts.left.D) < 1e-9 * Math.abs(sa.contacts.left.D) + 1e-12, `V=${V}: gate charge`);
    assert.ok(Math.abs(sb.phi[gSi] - sb.phi.at(-1) - (sa.phi[0] - sa.phi.at(-1))) < 1e-9, `V=${V}: surface potential`);
  }
});

// Silver nitrate between silver: electrodes as metal regions with reactions at their faces,
// against the same cell with the reactions at contacts.
const ions = [
  { name: 'Ag+', z: 1, cRef: 1000 },
  { name: 'NO3-', z: -1, cRef: 1000 },
];
const water = (epsr) => ({ epsr, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } });
const Ag = { metal: { species: 'e-', conductivity: 6e7 } };
const plating = { reactants: { 'Ag+': 1 }, electrons: 1, products: { Ag: 1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 };
const salt = { 'NO3-': 10, 'Ag+': 10 };

test('electrode reactions at metal faces: an Ag | AgNO₃ | Ag cell, as with contact electrodes', () => {
  for (const [epsr, face] of [
    [0, { phi: 'neutral' }],
    [78.5, { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 }],
  ]) {
    const contactPhi = face.phi === 'neutral' ? 'neutral' : { ...face.phi, zeroCharge: face.zeroCharge };
    for (const V of [0.05, -0.1]) {
      const reference = new Device({
        species: ions,
        materials: { water: water(epsr) },
        regions: [{ material: 'water', length: 20e-6, c0: salt }],
        contacts: { left: { V: 0, phi: contactPhi, reactions: [plating] }, right: { V, phi: contactPhi, reactions: [plating] } },
        grid: { hmin: 0.05e-9, hmax: 100e-9, ratio: 1.1 },
      }).solve();
      const metal = new Device({
        species: [...ions, { name: 'e-', z: -1 }],
        materials: { water: water(epsr), Ag },
        regions: [{ material: 'Ag', length: 1e-6 }, { material: 'water', length: 20e-6, c0: salt }, { material: 'Ag', length: 1e-6 }],
        interfaces: [{ ...face, reactions: [plating] }, { ...face, reactions: [plating] }],
        contacts: { left: collector(0), right: collector(V) },
        grid: { hmin: 0.05e-9, hmax: 100e-9, ratio: 1.1 },
      }).solve();
      assert.ok(reference.converged && metal.converged);
      assert.ok(Math.abs(metal.current / reference.current - 1) < 1e-9, `ε=${epsr}, V=${V}: ${metal.current} vs ${reference.current}`);
    }
  }
});

test('bipolar electrode: a floating metal passes current by reactions on both faces', () => {
  const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 };
  const end = { phi: { type: 'capacitive', C: 0.2, zeroCharge: 0.1 }, reactions: [plating] };
  const dev = new Device({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(78.5), Ag },
    regions: [
      { material: 'water', length: 10e-6, c0: salt },
      { material: 'Ag', length: 2e-6 },
      { material: 'water', length: 10e-6, c0: salt },
    ],
    interfaces: [{ ...stern, reactions: [plating] }, { ...stern, reactions: [plating] }],
    contacts: { left: { V: 0, ...end }, right: { V: 0, ...end } },
    grid: { hmin: 0.05e-9, hmax: 100e-9, ratio: 1.1 },
  });
  const mid = (dev.grid.regionStart[1] + dev.grid.regionEnd[1]) >> 1;
  for (const V of [0, 0.2, 1]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `V=${V}`);
    // By symmetry the floating metal sits halfway, and the same current passes everywhere.
    assert.ok(Math.abs(sol.V['e-'][mid] - V / 2) < 1e-9, `V=${V}: V_F = ${sol.V['e-'][mid]}`);
    assert.ok(Math.abs(sol.contacts.left.current - sol.current) <= 1e-9 * Math.abs(sol.current) + 1e-15);
    if (V === 0) assert.ok(Math.abs(sol.current) < 1e-15);
    for (const st of sol.conservation) if (st.spectator) assert.ok(Math.abs(st.drift) < 1e-12);
  }
});

test('metal definitions are checked', () => {
  const base = () => {
    const def = silicon();
    def.regions = [{ material: 'Au', length: 50e-9 }, { material: 'Si', length: 1e-6, fixedCharge: ND * FARADAY }];
    def.contacts.left = collector(0);
    return def;
  };
  let def = base();
  throwsDevice(def, /needs an explicit φ law/);
  def.interfaces = [{ phi: 'dipole', zeroCharge: 0.5 }];
  throwsDevice(def, /finite capacitance/);
  def.interfaces = [{ phi: { type: 'capacitive', C: 1 }, dipole: 0.1 }];
  throwsDevice(def, /a metal has no φ of its own/);
  def.interfaces = [{ phi: { type: 'capacitive', C: 1 }, vacuum: { left: { anchor: 'e-', offset: W }, right: { anchor: 'e-', offset: chi } } }];
  throwsDevice(def, /anchored to its Fermi level/);
  def.interfaces = [{ phi: { type: 'capacitive', C: 1 }, vacuum: { left: { anchor: 'fermi', offset: W }, right: { anchor: 'fermi', offset: chi } } }];
  throwsDevice(def, /only a metal has a Fermi-level anchor/);
  def = base();
  def.regions[0].c0 = { 'e-': 1 };
  throwsDevice(def, /no composition to give/);
  def = base();
  def.materials.Au.epsr = 1;
  throwsDevice(def, /a metal takes only/);
});
