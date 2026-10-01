import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, units } from '../src/index.js';
import { vacuumZeroCharge } from '../src/kit.js';

const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const throwsDevice = (def, pattern) =>
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && pattern.test(e.message), String(pattern));

test('metal region: ohmic, with a Fermi level and no φ or concentration, on one cell', () => {
  const sigma = 1e6, L = 1e-6, V = 1e-3;
  const dev = new Device({
    species: [{ name: 'e-', z: -1 }],
    materials: { Cu: { conductor: { species: 'e-', conductivity: sigma } } },
    regions: [{ material: 'Cu', length: L }],
    contacts: { left: collector(0), right: collector(V) },
    grid: { hmin: 1e-10 },
  });
  // Nothing is stored inside a metal, so its Fermi level is linear and one cell is exact.
  assert.equal(dev.grid.nNodes, 2);
  const sol = dev.solve();
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
    Au: { conductor: { species: 'e-', conductivity: 4e7 } },
  },
  contacts: { right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' } },
  grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
});

test('metal | semiconductor face: a Schottky barrier from work function and electron affinity, as at a contact', () => {
  const asContact = silicon();
  asContact.regions = [{ material: 'Si', length: 2e-6, fixedCharge: ND * FARADAY }];
  const zc = vacuumZeroCharge(asContact, W, { material: 'Si', anchor: 'e-', offset: chi }); // Schottky–Mott
  asContact.contacts.left = { ...collector(0), phi: { type: 'capacitive', C: 10 }, zeroCharge: zc };
  const asRegion = silicon();
  asRegion.regions = [{ material: 'Au', length: 50e-9 }, { material: 'Si', length: 2e-6, fixedCharge: ND * FARADAY }];
  asRegion.interfaces = [{ phi: { type: 'capacitive', C: 10 }, zeroCharge: zc }];
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
  asContact.contacts.left = { V: 0, phi: { type: 'capacitive', C: 1 / (1 / Cox + 1 / Cface) }, zeroCharge: zc };
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

// Silver nitrate between silver electrodes, as conductor regions with reactions at their faces.
const ions = [
  { name: 'Ag+', z: 1, cRef: 1000 },
  { name: 'NO3-', z: -1, cRef: 1000 },
];
const water = (epsr) => ({ epsr, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } });
const Ag = { conductor: { species: 'e-', conductivity: 6e7 } };
const k0 = 1e-3, alpha = 0.5;
// Ag⁺ + e⁻ ⇌ Ag at a silver face, with the metal on the given side.
const plating = (metal) => ({ [metal]: { 'e-': -1, Ag: 1 }, [metal === 'left' ? 'right' : 'left']: { 'Ag+': -1 }, fixed: { Ag: 0 }, k0, alpha });
const salt = { 'NO3-': 10, 'Ag+': 10 };
const silver = { material: 'Ag', length: 1e-6 };
const RT = 8.314462618 * 298.15;

test('electrode reactions at conductor faces: the rate law at the face state, and one current throughout', () => {
  for (const [epsr, phi] of [
    [0, { phi: 'neutral' }],
    [78.5, { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 }],
  ]) {
    const dev = new Device({
      species: [...ions, { name: 'e-', z: -1 }],
      materials: { water: water(epsr), Ag },
      regions: [silver, { material: 'water', length: 20e-6, c0: salt }, silver],
      interfaces: [{ ...phi, reactions: [plating('left')] }, { ...phi, reactions: [plating('right')] }],
      contacts: { left: collector(0), right: collector(0) },
      grid: { hmin: 0.05e-9, hmax: 100e-9, ratio: 1.1 },
    });
    const eq = dev.solve();
    assert.ok(eq.converged && Math.abs(eq.current) < 1e-12, `ε=${epsr}: zero current at equilibrium`);
    const gM = dev.grid.regionEnd[0], gW = dev.grid.regionStart[1];
    for (const V of [0.05, -0.1]) {
      dev.set({ contacts: { right: { V } } });
      const sol = dev.solve();
      assert.ok(sol.converged);
      // Reduction at the left face (forward) is current toward −x.
      const r = sol.interfaces[0].rates[0];
      assert.ok(Math.abs(-FARADAY * r / sol.current - 1) < 1e-9, `ε=${epsr}, V=${V}: F r = ${-FARADAY * r} vs ${sol.current}`);
      assert.ok(Math.abs(FARADAY * sol.interfaces[1].rates[0] / sol.current - 1) < 1e-9);
      // r = k0 (c/c_ref)^{1−α} (e^{αa} − e^{−(1−α)a}), a = (μ̄_Ag⁺ + μ̄_e⁻ − μ_Ag)/RT
      const a = (sol.mu['Ag+'][gW] + sol.mu['e-'][gM]) / RT;
      const want = k0 * (sol.c['Ag+'][gW] / 1000) ** (1 - alpha) * (Math.exp(alpha * a) - Math.exp(-(1 - alpha) * a));
      assert.ok(Math.abs(r / want - 1) < 1e-9, `rate law: ${r} vs ${want}`);
    }
  }
});

test('bipolar electrode: a floating metal passes current by reactions on both faces', () => {
  const face = (metal) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating(metal)] });
  const dev = new Device({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(78.5), Ag },
    regions: [
      silver,
      { material: 'water', length: 10e-6, c0: salt },
      { material: 'Ag', length: 2e-6 },
      { material: 'water', length: 10e-6, c0: salt },
      silver,
    ],
    interfaces: [face('left'), face('right'), face('left'), face('right')],
    contacts: { left: collector(0), right: collector(0) },
    grid: { hmin: 0.05e-9, hmax: 100e-9, ratio: 1.1 },
  });
  const mid = dev.grid.regionStart[2];
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
  def.interfaces = [{ phi: 'pinned', zeroCharge: 0.5 }];
  throwsDevice(def, /finite capacitance/);
  def.interfaces = [{ phi: { type: 'capacitive', C: 1 }, dipole: 0.1 }];
  throwsDevice(def, /a conductor has no φ of its own/);
  def.interfaces = [{ phi: { type: 'capacitive', C: 1 } }];
  throwsDevice(def, /zeroCharge/);
  def = base();
  def.regions[0].c0 = { 'e-': 1 };
  throwsDevice(def, /no composition to give/);
  def = base();
  def.materials.Au.epsr = 1;
  throwsDevice(def, /a conductor takes only/);
  def = base();
  def.interfaces = [{ phi: { type: 'capacitive', C: 1 }, zeroCharge: 0.5 }];
  new Device(def);
  def.regions[0].grid = { minCells: 10 };
  throwsDevice(def, /single cell/);
  delete def.regions[0].grid;
  def.regions[0].velocity = 1e-3;
  throwsDevice(def, /no flow or mixing/);
  delete def.regions[0].velocity;
  def.ports = [{ region: 0, from: 0, to: 1e-9, V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' } }];
  throwsDevice(def, /attaches to all of it/);
  def.ports = [{ region: 0, V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'conductance', G: 1, offset: 0 } } }];
  throwsDevice(def, /absent|only its carrier/);
});

test('a port on a metal: a wire to the whole metal, with its conductance per area', () => {
  // Copper from a grounded contact, tied to a port at V_p through G per area; the far end is
  // closed. The port's conductance is spread over the metal (half at each node of its one cell),
  // so the far half is reached through the metal's resistance L/σ: G_eff = G/2 + 1/(2/G + L/σ).
  const sigma = 1e6, L = 1e-6, Vp = 1e-3;
  for (const G of [1e10, 1e13]) {
    const sol = new Device({
      species: [{ name: 'e-', z: -1 }],
      materials: { Cu: { conductor: { species: 'e-', conductivity: sigma } } },
      regions: [{ material: 'Cu', length: L }],
      contacts: { left: collector(0), right: { phi: 'neutral' } },
      ports: [{ region: 0, V: Vp, terminal: 'e-', species: { 'e-': { type: 'conductance', G } } }],
    }).solve();
    assert.ok(sol.converged);
    const Geff = G / 2 + 1 / (2 / G + L / sigma);
    assert.ok(Math.abs(sol.ports[0].current / (Geff * Vp) - 1) < 1e-12, `G=${G}: ${sol.ports[0].current} vs ${Geff * Vp}`);
    assert.ok(Math.abs(sol.contacts.left.current + sol.ports[0].current) < 1e-12 * Math.abs(sol.ports[0].current));
  }
});

test('a port grounds a floating metal: a bipolar plate tied to the middle voltage carries no wire current', () => {
  // The symmetric bipolar cell of the test above, with the plate wired to V/2 (where it floats
  // anyway): nothing flows in the wire. Wired to 0 instead, the wire takes the left cell's current.
  const face = (metal) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating(metal)] });
  const cell = (Vwire, V) =>
    new Device({
      species: [...ions, { name: 'e-', z: -1 }],
      materials: { water: water(78.5), Ag },
      regions: [
        silver,
        { material: 'water', length: 10e-6, c0: salt },
        { name: 'plate', material: 'Ag', length: 2e-6 },
        { material: 'water', length: 10e-6, c0: salt },
        silver,
      ],
      interfaces: [face('left'), face('right'), face('left'), face('right')],
      contacts: { left: collector(0), right: collector(V) },
      ports: [{ region: 'plate', V: Vwire, terminal: 'e-', species: { 'e-': { type: 'conductance', G: 1e6 } } }],
      grid: { hmin: 0.05e-9, hmax: 100e-9, ratio: 1.1 },
    }).solve();
  const mid = cell(0.1, 0.2);
  assert.ok(mid.converged);
  assert.ok(Math.abs(mid.ports[0].current) < 1e-9 * Math.abs(mid.current), `${mid.ports[0].current} vs ${mid.current}`);
  const grounded = cell(0, 0.2);
  assert.ok(grounded.converged);
  // The wire takes the difference between the two cells' currents, and the plate sits near 0.
  assert.ok(Math.abs(grounded.contacts.left.current + grounded.ports[0].current - grounded.contacts.right.current) < 1e-9 * Math.abs(grounded.current));
  assert.ok(Math.abs(grounded.ports[0].current) > 0.1 * Math.abs(grounded.current));
});
