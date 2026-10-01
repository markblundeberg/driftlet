import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// O⁺ + e⁻ ⇌ R at an electrode facing a bath across a 20 µm unstirred layer, in a large
// excess of supporting electrolyte (K⁺Cl⁻), so transport of O⁺ and R is diffusion only.
const L = 20e-6, cO = 1, cR = 1, kp = 1e-5, DO = 1e-9, DR = 0.8e-9, alpha = 0.4;
const mu0 = { 'O+': 10e3, R: -86.485e3, 'K+': -283.3e3, 'Cl-': -131.2e3 };
// The electrode: a platinum region (a conductor for e⁻) behind the reaction face, read through
// a collector contact. O⁺ + e⁻ ⇌ R, with the metal on the given side.
const redoxRx = (metal) => ({ [metal]: { 'e-': -1 }, [metal === 'left' ? 'right' : 'left']: { 'O+': -1, R: 1 }, k0: kp * 1000, alpha });
const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const bath = { bath: { c: { 'O+': cO, R: cR, 'K+': 1000, 'Cl-': 1001 }, reference: 'Cl-' } };
const redox = (side, V, circuit) => ({
  species: [
    { name: 'O+', z: 1, cRef: 1000 },
    { name: 'R', z: 0, cRef: 1000 },
    { name: 'K+', z: 1, cRef: 1000 },
    { name: 'Cl-', z: -1, cRef: 1000 },
    { name: 'e-', z: -1 },
  ],
  materials: {
    water: {
      epsr: 78.5,
      species: {
        'O+': { D: DO, mu0: mu0['O+'] },
        R: { D: DR, mu0: mu0.R },
        'K+': { D: 1.96e-9, mu0: mu0['K+'] },
        'Cl-': { D: 2.03e-9, mu0: mu0['Cl-'] },
      },
    },
    Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } },
  },
  // phi 'neutral': no double layer at the electrode, to compare with the textbook result
  ...(side === 'left'
    ? {
        regions: [{ material: 'Pt', length: 1e-6 }, { material: 'water', length: L }],
        interfaces: [{ phi: 'neutral', reactions: [redoxRx('left')] }],
        contacts: { left: collector(V), right: bath },
      }
    : {
        regions: [{ material: 'water', length: L }, { material: 'Pt', length: 1e-6 }],
        interfaces: [{ phi: 'neutral', reactions: [redoxRx('right')] }],
        contacts: { left: bath, right: collector(V) },
      }),
  circuit,
  grid: { hmin: 0.05e-9, hmax: 200e-9, ratio: 1.15 },
});
// Bath φ from its Cl⁻ reference at 0 V, hence the equilibrium electrode voltage and E°′.
const phiB = (mu0['Cl-'] + RT * Math.log(1001 / 1000)) / FARADAY;
const E0p = (mu0['O+'] - mu0.R) / FARADAY + phiB;
const Veq = E0p + VT * Math.log(cO / cR);

test('Butler–Volmer electrode: zero current at the Nernst potential, mixed-control I–V', () => {
  const dev = new Device(redox('left', Veq));
  const eq = dev.solve();
  assert.ok(eq.converged);
  assert.ok(Math.abs(eq.current) < 1e-12);
  for (const name of ['O+', 'R']) {
    const mu = eq.mu[name].filter(Number.isFinite); // (absent from the platinum)
    assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-12, `${name} flat at equilibrium`);
  }
  // Textbook: r = k′(c_O,s e^{αx} − c_R,s e^{−(1−α)x}) with linear diffusion to the bath,
  // x = (E°′ − V)/V_T, solved in closed form. Reduction at the left electrode is current
  // toward −x.
  for (const eta of [0.02, 0.1, 0.2, -0.05, -0.2]) {
    dev.set({ contacts: { left: { V: Veq + eta } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `η=${eta}`);
    const x = (E0p - (Veq + eta)) / VT;
    const ea = Math.exp(alpha * x), eb = Math.exp(-(1 - alpha) * x);
    const r = (kp * (cO * ea - cR * eb)) / (1 + ((kp * L) / DO) * ea + ((kp * L) / DR) * eb);
    const want = -FARADAY * r;
    assert.ok(Math.abs(sol.current / want - 1) < 5e-4, `η=${eta}: ${sol.current} vs ${want}`);
    assert.ok(Math.abs(sol.contacts.left.current - sol.current) < 1e-9 * Math.abs(want));
  }
});

test('galvanostatic kinetic electrode', () => {
  const dev = new Device(redox('right', Veq, { mode: 'current', I: 2 }));
  const sol = dev.solve();
  assert.ok(sol.converged);
  assert.ok(Math.abs(sol.current / 2 - 1) < 1e-9);
  const V = sol.contacts.right.V;
  dev.set({ circuit: { mode: 'voltage', I: undefined }, contacts: { right: { V } } });
  assert.ok(Math.abs(dev.solve().current / 2 - 1) < 1e-8, 'voltage mode at that V returns the same current');
});

// A single-ion conductor (mobile Li⁺ on a fixed negative background) in two materials, with
// the dipole chosen so both sides hold c = X at equilibrium.
const X = 100, D = 1e-11, Lr = 1e-6;
const ionic = (itf) => ({
  species: [{ name: 'Li+', z: 1, cRef: 1000 }],
  materials: {
    A: { epsr: 30, species: { 'Li+': { D, mu0: 0 } } },
    B: { epsr: 30, species: { 'Li+': { D, mu0: 2 * RT } } },
  },
  regions: [
    { material: 'A', length: Lr, fixedCharge: -X * FARADAY },
    { material: 'B', length: Lr, fixedCharge: -X * FARADAY },
  ],
  interfaces: [{ dipole: -2 * VT, ...itf }],
  contacts: {
    left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' },
    right: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' },
  },
  grid: { hmin: 0.1e-9, hmax: 20e-9, ratio: 1.15 },
});
const Rbulk = (2 * Lr) / ((FARADAY * FARADAY * X * D) / RT);

test('interface conductance adds a series resistance 1/G', () => {
  const G = 50;
  const dev = new Device(ionic({ species: { 'Li+': { type: 'conductance', G } } }));
  dev.solve();
  for (const V of [0.01, -0.05, 0.1]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.current / (-V / (Rbulk + 1 / G)) - 1) < 1e-4, `V=${V}`);
  }
});

test('interface ion transfer (Butler–Volmer): equilibrium, and the rate law at the interface state', () => {
  const k0 = 1e-3, a = 0.3;
  const transfer = { left: { 'Li+': -1 }, right: { 'Li+': 1 }, k0, alpha: a };
  // Li⁺ exists on both sides, so its link must be given: here it crosses only by the reaction.
  assert.throws(() => new Device(ionic({ reactions: [transfer] })), (e) => e instanceof DeviceError && /give its link explicitly/.test(e.message));
  const dev = new Device(ionic({ species: { 'Li+': 'blocked' }, reactions: [transfer] }));
  const eq = dev.solve();
  assert.ok(Math.abs(eq.current) < 1e-15);
  const mu = eq.mu['Li+'];
  assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-12);
  const g = dev.grid, gL = g.regionEnd[0], gR = g.regionStart[1];
  for (const V of [0.01, 0.2, -0.1]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged);
    // The transfer overpotential is carried by thin space-charge layers either side (the φ
    // jump is pinned at the dipole), so check the rate law at the actual interface state.
    const af = (sol.mu['Li+'][gL] - sol.mu['Li+'][gR]) / RT;
    const r =
      k0 * (sol.c['Li+'][gL] / 1000) ** (1 - a) * (sol.c['Li+'][gR] / 1000) ** a * (Math.exp(a * af) - Math.exp(-(1 - a) * af));
    assert.ok(Math.abs((FARADAY * r) / sol.current - 1) < 1e-8, `V=${V}`);
    assert.ok(Math.abs((FARADAY * sol.interfaces[0].rates[0]) / sol.current - 1) < 1e-8);
  }
});

test('kinetic definitions are checked', () => {
  const throwsDevice = (def, pattern) =>
    assert.throws(() => new Device(def), (e) => e instanceof DeviceError && pattern.test(e.message));
  let def = redox('left', 0);
  def.interfaces[0].reactions[0].left['e-'] = -2;
  throwsDevice(def, /charge is not balanced/);
  def = redox('left', 0);
  def.interfaces[0].reactions[0].alpha = 1.2;
  throwsDevice(def, /alpha must be between 0 and 1/);
  def = redox('left', 0);
  def.interfaces[0].reactions[0].left['O+'] = -1;
  throwsDevice(def, /left\.O\+: absent from 'Pt'/);
  def = redox('left', 0);
  def.interfaces[0].reactions[0].right.Pt = 1;
  throwsDevice(def, /not a species, so give its μ/);
  def = ionic({ species: { 'Li+': { type: 'conductance' } } });
  throwsDevice(def, /\.G must be a positive number/);
});
