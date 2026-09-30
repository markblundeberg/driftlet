import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// Silicon pn junction with ohmic contacts: e⁻ and h⁺ both fixed at the metal's voltage
// (V_h = V_e), and a locally neutral contact node.
const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19), Eg = 1.12;
const ND = units.perCm3(1e17), NA = units.perCm3(1e16);
const Dn = 36e-4, Dp = 12e-4, L = 2e-6;
const ni2 = Nc * Nv * Math.exp(-Eg / VT);
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
const pn = () => ({
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: { epsr: 11.7, species: { 'e-': { D: Dn, mu0: 0, cRef: Nc }, 'h+': { D: Dp, mu0: units.eV(Eg), cRef: Nv } } },
  },
  regions: [
    { name: 'n', material: 'Si', length: L, fixedCharge: ND * FARADAY },
    { name: 'p', material: 'Si', length: L, fixedCharge: -NA * FARADAY },
  ],
  contacts: { left: ohmic(0), right: ohmic(0) },
  grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
});
const nn = ND / 2 + Math.sqrt((ND * ND) / 4 + ni2);
const pp = NA / 2 + Math.sqrt((NA * NA) / 4 + ni2);
const Vbi = VT * Math.log((nn * pp) / ni2);

test('pn junction: exact built-in potential, flat μ̄, zero current at equilibrium', () => {
  const sol = new Device(pn()).solve();
  assert.ok(sol.converged);
  const vbi = sol.phi[0] - sol.phi[sol.phi.length - 1];
  assert.ok(Math.abs(vbi - Vbi) < 1e-12, `Vbi ${vbi} vs ${Vbi}`);
  for (const name of ['e-', 'h+']) {
    const mu = sol.mu[name];
    assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-12, `${name} μ̄ flat`);
  }
  assert.ok(Math.abs(sol.current) < 1e-12);
  // Ohmic contacts: the contact node is neutral, with n p = n_i².
  assert.ok(Math.abs(sol.c['e-'][0] / nn - 1) < 1e-12);
  assert.ok(Math.abs((sol.c['e-'][0] * sol.c['h+'][0]) / ni2 - 1) < 1e-12);
});

test('pn junction: depletion charge matches the depletion approximation', () => {
  const dev = new Device(pn());
  const sol = dev.solve();
  const g = dev.grid;
  let Qn = 0;
  for (let k = g.regionStart[0]; k <= g.regionEnd[0]; k++) {
    Qn += g.vol[k] * FARADAY * (ND - sol.c['e-'][k] + sol.c['h+'][k]);
  }
  const eps = 11.7 * EPS0;
  const Q = Math.sqrt(2 * FARADAY * eps * (Vbi - 2 * VT) * ((NA * ND) / (NA + ND)));
  assert.ok(Math.abs(Qn / Q - 1) < 0.01, `Qn ${Qn} vs ${Q}`);
});

test('pn junction: short-diode J–V (contacts as the only recombination) and uniform current', () => {
  const dev = new Device(pn());
  dev.solve();
  const eps = 11.7 * EPS0;
  for (const V of [-0.5, -0.2, 0.1, 0.2, 0.3, 0.4, 0.5]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `V=${V}`);
    // Current in steady state is the same through both contacts, even when tiny.
    const IL = sol.contacts.left.current, IR = sol.contacts.right.current;
    assert.ok(Math.abs(IL - IR) <= 1e-8 * Math.abs(IR), `V=${V}: ${IL} vs ${IR}`);
    if (V > 0.45) continue; // high injection: the analytic law no longer applies
    const W = Math.sqrt(((2 * eps * (Vbi - V - 2 * VT)) / FARADAY) * (1 / NA + 1 / ND));
    const xn = (W * NA) / (NA + ND), xp = (W * ND) / (NA + ND);
    const J = FARADAY * ni2 * (Dn / (NA * (L - xp)) + Dp / (ND * (L - xn))) * Math.expm1(V / VT);
    // Forward current flows from p (right) to n (left): negative toward +x.
    assert.ok(Math.abs(-IR / J - 1) < 2e-3, `V=${V}: ${-IR} vs ${J}`);
  }
});

// Salt baths as contacts: each fixes every ion's μ̄ from its composition, through a
// reference species (Cl⁻, as for a Ag/AgCl electrode) at the terminal voltage.
const salt = [
  { name: 'Na+', z: 1, cRef: 1000 },
  { name: 'Cl-', z: -1, cRef: 1000 },
];
const water = { epsr: 78.5, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } };
const bath = (c, V = 0) => ({ V, bath: { c: { 'Na+': c, 'Cl-': c }, reference: 'Cl-' } });

test('bath contact reproduces its composition and puts Cl⁻ at the terminal voltage', () => {
  const dev = new Device({
    species: salt,
    materials: { water },
    regions: [{ material: 'water', length: 100e-9 }],
    contacts: { left: bath(10, 0.2), right: bath(10, 0.2) },
    grid: { hmin: 0.1e-9, hmax: 5e-9 },
  });
  const sol = dev.solve();
  assert.ok(sol.converged);
  for (const g of [0, sol.x.length >> 1, sol.x.length - 1]) {
    assert.ok(Math.abs(sol.c['Na+'][g] / 10 - 1) < 1e-12);
    assert.ok(Math.abs(sol.c['Cl-'][g] / 10 - 1) < 1e-12);
    assert.ok(Math.abs(sol.V['Cl-'][g] - 0.2) < 1e-12);
  }
  assert.ok(Math.abs(sol.current) < 1e-15);
});

test('Donnan equilibrium between two real baths', () => {
  const X = 50; // mol/m³ fixed charge; Cl⁻ is excluded from the membrane
  const M = { epsr: 20, species: { 'Na+': { D: 1e-10, mu0: water.species['Na+'].mu0 + 2 * RT } } };
  const dev = new Device({
    species: salt,
    materials: { water, M },
    regions: [
      { material: 'water', length: 100e-9 },
      { material: 'M', length: 200e-9, fixedCharge: -X * FARADAY },
      { material: 'water', length: 100e-9 },
    ],
    interfaces: [{ dipole: 0 }, { dipole: 0 }],
    contacts: { left: bath(10), right: bath(10) },
    grid: { hmin: 0.05e-9, hmax: 5e-9, ratio: 1.1 },
  });
  const sol = dev.solve();
  assert.ok(sol.converged);
  const g = dev.grid;
  const m = (g.regionStart[1] + g.regionEnd[1]) >> 1;
  // Membrane bulk: neutral, and Na⁺ in equilibrium with a 10 mM bath.
  assert.ok(Math.abs(sol.c['Na+'][m] / X - 1) < 1e-9);
  const dphi = (sol.phi[m] - sol.phi[0]) / VT;
  assert.ok(Math.abs(dphi - (Math.log(10 / X) - 2)) < 1e-9, `Donnan ${dphi}`);
  assert.ok(Math.abs(sol.current) < 1e-15);
});

test('contact definitions: offsets and φ conditions are never silently defaulted', () => {
  const base = () => ({
    species: salt,
    materials: { water },
    regions: [{ material: 'water', length: 1e-7 }],
    contacts: { left: bath(10), right: bath(10) },
  });
  const throwsDevice = (def, pattern) =>
    assert.throws(() => new Device(def), (err) => err instanceof DeviceError && pattern.test(err.message));

  let def = base();
  def.contacts.left = { terminal: 'Cl-', species: { 'Cl-': 'equilibrium', 'Na+': 'equilibrium' }, phi: 'bulk' };
  throwsDevice(def, /contacts\.left\.species\.Na\+\.offset: give V_i − V_terminal/);

  def = base();
  def.contacts.left = { terminal: 'Cl-', species: { 'Cl-': 'equilibrium' } };
  throwsDevice(def, /contacts\.left\.phi: a contact with connected species needs an explicit φ law/);

  def = base();
  def.contacts.left = { bath: { c: { 'Na+': 10, 'Cl-': 9 }, reference: 'Cl-' } };
  throwsDevice(def, /bath: composition is not neutral/);

  def = base();
  def.contacts.left = { ...bath(10), species: { 'Na+': 'blocked' } };
  throwsDevice(def, /either bath or species links, not both/);

  def = base();
  def.contacts.left = { phi: { type: 'capacitive', C: 0.1, zeroCharge: 0, V: 1 } };
  throwsDevice(def, /gate voltage is the contact's terminal voltage/);

  def = base();
  def.contacts.left = { species: { 'Na+': { type: 'equilibrium', mu: 0 } }, terminal: 'Na+', phi: 'bulk' };
  throwsDevice(def, /charged species is held by an offset/);
});

test('pinned (dipole) contact: a Schottky barrier directly on n-Si', () => {
  // Electrons in equilibrium with the metal (V_e = V) and φ pinned at φ_edge = V − zeroCharge.
  // With μ°_e = 0, the surface density is N_c e^{−zeroCharge/V_T}: zeroCharge is the barrier.
  const phiB = 0.7, NDs = units.perCm3(1e16);
  const def = pn();
  def.regions = [{ name: 'n', material: 'Si', length: 2e-6, fixedCharge: NDs * FARADAY }];
  def.contacts.left = { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: { type: 'dipole', zeroCharge: phiB } };
  const dev = new Device(def);
  const eps = 11.7 * EPS0;
  const Vbi = phiB - VT * Math.log(Nc / NDs);
  for (const V of [0, 0.1, -0.2]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `V=${V}`);
    assert.ok(Math.abs(sol.c['e-'][0] / (Nc * Math.exp(-phiB / VT)) - 1) < 1e-12, `V=${V}: surface density`);
    let Q = 0;
    for (let k = 0; k < sol.x.length; k++) Q += dev.grid.vol[k] * FARADAY * (NDs - sol.c['e-'][k] + sol.c['h+'][k]);
    const Qd = Math.sqrt(2 * FARADAY * eps * NDs * (Vbi + V - VT)); // right contact at V: reverse for V > 0
    assert.ok(Math.abs(Q / Qd - 1) < 0.02, `V=${V}: depletion charge ${Q} vs ${Qd}`);
    // Gauss: device charge = D_right − D_left. The metal plate holds −D_left; under bias, the ohmic
    // field at the far (bulk) contact holds the rest.
    assert.ok(Math.abs(Q - (sol.contacts.right.D - sol.contacts.left.D)) < 1e-9 * Q, 'Gauss');
    assert.ok(Math.abs(sol.gates.left.charge - sol.contacts.left.D) < 1e-15, 'plate charge is the contact displacement');
  }
});
