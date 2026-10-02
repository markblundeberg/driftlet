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
// The collector holds the platinum's electrons at V, or drives a current I into the device.
const collector = (V, I) => ({ ...(I === undefined ? { V } : { I }), terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const bath = { bath: { c: { 'O+': cO, R: cR, 'K+': 1000, 'Cl-': 1001 }, reference: 'Cl-' } };
const redox = (side, V, I) => ({
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
        contacts: { left: bath, right: collector(V, I) },
      }),
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
  const dev = new Device(redox('right', Veq, -2)); // 2 A/m² toward +x, out at the right
  const sol = dev.solve();
  assert.ok(sol.converged);
  assert.ok(Math.abs(sol.current / 2 - 1) < 1e-9);
  const V = sol.contacts.right.V;
  dev.set({ contacts: { right: { V, I: undefined } } });
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

test('a redox couple between inert electrodes conserves its total: equilibrium at its Nernst level, and steady current', () => {
  // Pt | Fe³⁺, Fe²⁺, H⁺, Cl⁻ | Pt, a closed cell. Fe³⁺ + e⁻ ⇌ Fe²⁺ exchanges electrons with
  // the metals but conserves the iron, so its amount stays as given, whatever the electrodes do.
  const c3 = 5, c2 = 20;
  const water = { 'H+': { D: 9.3e-9, mu0: 0 }, 'Fe3+': { D: 0.6e-9, mu0: -4.7e3 }, 'Fe2+': { D: 0.72e-9, mu0: -78.9e3 }, 'Cl-': { D: 2.0e-9, mu0: -131.2e3 } };
  const redoxFace = (metal) => ({ reactions: [{ [metal]: { 'e-': -1 }, [metal === 'left' ? 'right' : 'left']: { 'Fe3+': -1, 'Fe2+': 1 }, k0: 1e-3 }] });
  const cell = (V) =>
    new Device({
      species: [
        { name: 'H+', z: 1, cRef: 1000 },
        { name: 'Fe3+', z: 3, cRef: 1000 },
        { name: 'Fe2+', z: 2, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
        { name: 'e-', z: -1 },
      ],
      materials: { water: { epsr: 0, species: water }, Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } } },
      regions: [
        { material: 'Pt', length: 1e-6 },
        { material: 'water', length: 20e-6, c0: { 'H+': 100, 'Fe3+': c3, 'Fe2+': c2, 'Cl-': 100 + 3 * c3 + 2 * c2 } },
        { material: 'Pt', length: 1e-6 },
      ],
      interfaces: [redoxFace('left'), redoxFace('right')],
      contacts: { left: collector(0), right: collector(V) },
      grid: { minCells: 100 },
    });
  const iron = (sol) => sol.conservation.filter((st) => st.species === 'Fe3+' || st.species === 'Fe2+').reduce((a, st) => a + st.amount, 0);
  const total = (c3 + c2) * 20e-6;

  const eq = cell(0).solve();
  assert.ok(eq.converged);
  assert.ok(Math.abs(iron(eq) / total - 1) < 1e-10, `iron ${iron(eq)} vs ${total}`);
  const g = eq.x.length >> 1;
  // The electrode sits at the couple's level: E° + (RT/F) ln(c₃/c₂) above V°(SHE) = φ.
  const E0 = (-4.7e3 + 78.9e3) / FARADAY;
  assert.ok(Math.abs(eq.V['e-'][0] - eq.phi[g] - (E0 + VT * Math.log(c3 / c2))) < 1e-9);

  const run = cell(0.1).solve();
  assert.ok(run.converged && Math.abs(run.current) > 0);
  assert.ok(Math.abs(iron(run) / total - 1) < 1e-10, `iron ${iron(run)} vs ${total}`);
});

test('with complexation too, the closed cell solves straight to its steady state, from cold or warm', () => {
  // Pt | Fe³⁺, Fe²⁺, FeCl²⁺, K⁺, Cl⁻ | Pt with Fe³⁺ + Cl⁻ ⇌ FeCl²⁺ in the bulk: the total iron and
  // the total chloride are conserved combinations of reacting stretches, solved as constraints.
  const water = {
    'K+': { D: 1.96e-9, mu0: -283.3e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 },
    'Fe3+': { D: 0.6e-9, mu0: -4.7e3 }, 'Fe2+': { D: 0.72e-9, mu0: -78.9e3 }, 'FeCl2+': { D: 0.7e-9, mu0: -155.9e3 },
  };
  const redoxFace = (metal) => ({ reactions: [{ [metal]: { 'e-': -1 }, [metal === 'left' ? 'right' : 'left']: { 'Fe3+': -1, 'Fe2+': 1 }, k0: 1e-3 }] });
  const def = {
    species: [
      { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'Fe3+', z: 3, cRef: 1000 },
      { name: 'Fe2+', z: 2, cRef: 1000 }, { name: 'FeCl2+', z: 2, cRef: 1000 }, { name: 'e-', z: -1 },
    ],
    materials: { water: { epsr: 0, species: water }, Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } } },
    regions: [
      { material: 'Pt', length: 1e-6 },
      { material: 'water', length: 100e-6, c0: { 'K+': 500, 'Cl-': 527, 'Fe3+': 5, 'Fe2+': 5, 'FeCl2+': 1 } },
      { material: 'Pt', length: 1e-6 },
    ],
    interfaces: [redoxFace('left'), redoxFace('right')],
    bulkReactions: [{ equation: 'Fe3+ + Cl- = FeCl2+', kf: { water: 1e3 } }],
    contacts: { left: collector(0), right: collector(0.1) },
    grid: { hmin: 1e-8, hmax: 2e-6 },
  };
  const cold = new Device(def).solve();
  const run = new Device(def);
  run.advance(1e6);
  const warm = run.solve();
  // From its own steady state the solve used to creep through giant steps and report failure.
  assert.ok(cold.converged && warm.converged);
  assert.ok(Math.abs(cold.current / warm.current - 1) < 1e-10, `${cold.current} vs ${warm.current}`);
  const iron = (s) => s.conservation.filter((st) => /^Fe/.test(st.species)).reduce((a, st) => a + st.amount, 0);
  assert.ok(Math.abs(iron(cold) / (11 * 100e-6) - 1) < 1e-10, `iron ${iron(cold)}`);
});
