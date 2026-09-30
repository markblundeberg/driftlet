import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// Long silicon pn diode with band-to-band recombination e⁻ + h⁺ ⇌ ∅, r = k(np − n_i²).
// Minority lifetimes are τ_n = 1/(k N_A) and τ_p = 1/(k N_D); regions are many diffusion
// lengths long.
const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19), Eg = 1.12;
const ND = units.perCm3(1e17), NA = units.perCm3(1e16);
const Dn = 36e-4, Dp = 12e-4, Ln = 5e-6;
const kf = Dn / (Ln * Ln * NA);
const Lp = Math.sqrt(Dp / (kf * ND));
const ni2 = Nc * Nv * Math.exp(-Eg / VT);
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
const diode = () => ({
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: { epsr: 11.7, species: { 'e-': { D: Dn, mu0: 0, cRef: Nc }, 'h+': { D: Dp, mu0: units.eV(Eg), cRef: Nv } } },
  },
  regions: [
    { material: 'Si', length: 40e-6, fixedCharge: ND * FARADAY },
    { material: 'Si', length: 40e-6, fixedCharge: -NA * FARADAY },
  ],
  bulkReactions: [{ reactants: { 'e-': 1, 'h+': 1 }, kf: { Si: kf } }],
  contacts: { left: ohmic(0), right: ohmic(0) },
  grid: { hmin: 0.5e-9, hmax: 100e-9, ratio: 1.1 },
});

test('long pn diode: Shockley J–V with diffusion and depletion-region recombination', () => {
  const dev = new Device(diode());
  const eq = dev.solve();
  assert.ok(eq.converged);
  assert.ok(Math.abs(eq.current) < 1e-12, 'no current at equilibrium');
  for (const name of ['e-', 'h+']) {
    const mu = eq.mu[name];
    assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-12, `${name} flat with recombination on`);
  }
  const nn = ND / 2 + Math.sqrt((ND * ND) / 4 + ni2), pp = NA / 2 + Math.sqrt((NA * NA) / 4 + ni2);
  const Vbi = VT * Math.log((nn * pp) / ni2);
  for (const V of [-0.3, 0.1, 0.2, 0.3, 0.4]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `V=${V}`);
    // With flat quasi-Fermi levels across the depletion width W, it adds k·n_i²·W.
    const W = Math.sqrt(((2 * 11.7 * EPS0 * (Vbi - V - 2 * VT)) / FARADAY) * (1 / NA + 1 / ND));
    const J = FARADAY * ni2 * (Dn / (NA * Ln) + Dp / (ND * Lp) + kf * W) * Math.expm1(V / VT);
    assert.ok(Math.abs(-sol.current / J - 1) < 2e-3, `V=${V}: ${-sol.current} vs ${J}`);
    assert.ok(Math.abs(sol.contacts.left.current / sol.current - 1) < 1e-8);
  }
});

// Water autoionisation H⁺ + OH⁻ ⇌ H₂O in a closed, gated box of dilute NaCl.
const mu0 = { 'H+': 0, 'OH-': -157.24e3, 'Na+': -261.9e3, 'Cl-': -131.2e3 };
const muH2O = -237.13e3;
const Kw = 1e6 * Math.exp((muH2O - mu0['H+'] - mu0['OH-']) / RT); // (mol/m³)², ≈ 1e-8
const waterBox = (c0) => ({
  species: Object.keys(mu0).map((name) => ({ name, z: name.endsWith('+') ? 1 : -1, cRef: 1000 })),
  materials: {
    water: {
      epsr: 78.5,
      species: {
        'H+': { D: 9.31e-9, mu0: mu0['H+'] },
        'OH-': { D: 5.27e-9, mu0: mu0['OH-'] },
        'Na+': { D: 1.33e-9, mu0: mu0['Na+'] },
        'Cl-': { D: 2.03e-9, mu0: mu0['Cl-'] },
      },
    },
  },
  regions: [{ material: 'water', length: 100e-9, c0 }],
  bulkReactions: [{ reactants: { 'H+': 1, 'OH-': 1 }, products: { H2O: 1 }, fixed: { H2O: muH2O }, kf: { water: 1.4e8 } }],
  contacts: {
    left: { phi: { type: 'capacitive', C: 0.2, zeroCharge: 0 } },
    right: { phi: { type: 'capacitive', C: 0.2, zeroCharge: 0 } },
  },
  grid: { hmin: 0.05e-9, hmax: 5e-9, ratio: 1.1 },
});

test('water: mass action c(H⁺)c(OH⁻) = K_w at equilibrium, conserving the H⁺ − OH⁻ moiety', () => {
  assert.ok(Math.abs(Kw / 1e-8 - 1) < 0.02, `K_w from standard potentials: ${Kw}`);
  const c0 = { 'H+': 1e-2, 'OH-': 1e-3, 'Na+': 10, 'Cl-': 10.009 };
  const dev = new Device(waterBox(c0));
  const sol = dev.solve();
  assert.ok(sol.converged);
  const m = sol.x.length >> 1;
  assert.ok(Math.abs((sol.c['H+'][m] * sol.c['OH-'][m]) / Kw - 1) < 1e-9);
  const amount = (name) => sol.conservation.find((st) => st.species === name).amount;
  const L = 100e-9;
  assert.ok(Math.abs((amount('H+') - amount('OH-')) / ((c0['H+'] - c0['OH-']) * L) - 1) < 1e-10);
  // Na⁺ and Cl⁻ don't react: still plain spectators.
  for (const name of ['Na+', 'Cl-']) {
    const st = sol.conservation.find((s) => s.species === name);
    assert.ok(st.spectator && Math.abs(st.drift) < 1e-12);
  }
  assert.ok(!sol.conservation.find((s) => s.species === 'H+').spectator);
});

test('bulk reaction definitions are checked', () => {
  const throwsDevice = (def, pattern) =>
    assert.throws(() => new Device(def), (e) => e instanceof DeviceError && pattern.test(e.message));
  let def = waterBox({ 'H+': 1, 'OH-': 1, 'Na+': 1, 'Cl-': 1 });
  def.bulkReactions = [{ reactants: { 'H+': 1 }, products: { H2O: 1 }, fixed: { H2O: muH2O }, kf: { water: 1 } }];
  throwsDevice(def, /charge is not balanced/);
  def.bulkReactions = [{ reactants: { 'H+': 1, 'OH-': 1 }, products: { H2O: 1 }, kf: { water: 1 } }];
  throwsDevice(def, /H2O: not a species, so give its μ/);
  def.bulkReactions = [{ reactants: { 'H+': 1, 'OH-': 1 }, products: { H2O: 1 }, fixed: { H2O: muH2O }, kf: { ice: 1 } }];
  throwsDevice(def, /kf\.ice: unknown material/);
  def = diode();
  def.bulkReactions = [{ reactants: { 'e-': 1.5, 'h+': 1 }, kf: { Si: 1 } }];
  throwsDevice(def, /must be a positive integer/);
});
