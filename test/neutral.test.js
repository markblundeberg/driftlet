import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

const salt = [
  { name: 'Na+', z: 1, cRef: 1000 },
  { name: 'Cl-', z: -1, cRef: 1000 },
];
const water = (epsr = 78.5) => ({ epsr, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } });
const bath = (c) => ({ bath: { c: { 'Na+': c, 'Cl-': c }, reference: 'Cl-' } });

test('ε = 0: a liquid junction is exactly quasi-neutral (Planck EMF)', () => {
  const Dp = 1.33e-9, Dm = 2.03e-9, tp = Dp / (Dp + Dm);
  const sol = new Device({
    species: salt,
    materials: { water: water(0) },
    regions: [{ material: 'water', length: 10e-6 }],
    contacts: { left: bath(100), right: { ...bath(10), I: 0 } },
    grid: { minCells: 400 },
  }).solve();
  assert.ok(sol.converged);
  const emf = 2 * tp * VT * Math.log(10);
  assert.ok(Math.abs(sol.terminalVoltage / emf - 1) < 1e-5, `${sol.terminalVoltage} vs ${emf}`);
  for (let g = 0; g < sol.x.length; g++) assert.ok(Math.abs(sol.c['Na+'][g] / sol.c['Cl-'][g] - 1) < 1e-12, 'neutral everywhere');
  assert.deepEqual(sol.warnings, []);
});

test('ε = 0: concentration polarization follows tanh(V/4V_T) with no overlimiting', () => {
  const Dp = 1.65e-9, c0 = 10, L = 20e-6;
  const el = (V = 0) => ({ V, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
  const dev = new Device({
    species: [
      { name: 'Ag+', z: 1, cRef: 1000 },
      { name: 'NO3-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 0, species: { 'Ag+': { D: Dp, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
    regions: [{ material: 'water', length: L, c0: { 'NO3-': c0 } }],
    contacts: { left: el(), right: el() },
    grid: { hmin: 1e-9, hmax: 20e-9, ratio: 1.1 }, // ln c varies steeply toward the depleted end
  });
  const iLim = (4 * FARADAY * Dp * c0) / L;
  for (const V of [0.05, 0.2, 0.4]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged);
    assert.ok(Math.abs(-sol.current / iLim - Math.tanh(V / (4 * VT))) < 1e-4, `V=${V}`);
  }
});

// Water | cation-exchange membrane | water between baths.
const X = 50;
const donnan = (epsr, interfaces) => {
  const def = {
    species: salt,
    materials: {
      water: water(epsr),
      M: { epsr: (epsr * 20) / 78.5, species: { 'Na+': { D: 1e-10, mu0: -261.9e3 + 2 * RT } } },
    },
    regions: [
      { material: 'water', length: 100e-9 },
      { material: 'M', length: 200e-9, fixedCharge: -X * FARADAY },
      { material: 'water', length: 100e-9 },
    ],
    contacts: { left: bath(10), right: bath(10) },
    grid: { hmin: 0.05e-9, hmax: 5e-9, ratio: 1.1 },
  };
  if (interfaces) def.interfaces = interfaces;
  return def;
};

test('ε = 0 materials meet at neutral faces: Donnan jump, no alignment needed, neutral right up to the face', () => {
  const dev = new Device(donnan(0));
  assert.equal(dev.model.interfaces[0].phi.type, 'neutral');
  const sol = dev.solve();
  assert.ok(sol.converged);
  const g = dev.grid, gL = g.regionEnd[0], gR = g.regionStart[1];
  assert.ok(Math.abs(sol.c['Na+'][gL] / 10 - 1) < 1e-9 && Math.abs(sol.c['Cl-'][gL] / 10 - 1) < 1e-9);
  assert.ok(Math.abs(sol.c['Na+'][gR] / X - 1) < 1e-9);
  // The jump is the Donnan potential, taken entirely at the face.
  assert.ok(Math.abs((sol.phi[gR] - sol.phi[gL]) / VT - (Math.log(10 / X) - 2)) < 1e-9);
  assert.equal(sol.interfaces[0].D, 0);
});

test('an explicit neutral face at ε > 0 gives the same Donnan jump; the ε → 0 limit is smooth', () => {
  const faces = [{ phi: 'neutral' }, { phi: 'neutral' }];
  for (const epsr of [78.5, 1e-3]) {
    const dev = new Device(donnan(epsr, faces));
    const sol = dev.solve();
    const g = dev.grid, m = (g.regionStart[1] + g.regionEnd[1]) >> 1;
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.c['Na+'][m] / X - 1) < 1e-9);
    assert.ok(Math.abs(sol.interfaces[0].D) < 1e-15);
  }
  assert.throws(
    () => new Device(donnan(78.5, [{ phi: 'neutral', dipole: 0.1 }, { phi: 'neutral' }])),
    (e) => e instanceof DeviceError && /neutral interface has a free φ jump/.test(e.message),
  );
});

test('unresolved double layers at a pinned-dipole face are reported as warnings', () => {
  const fine = new Device(donnan(78.5, [{ dipole: 0 }, { dipole: 0 }])).solve();
  assert.deepEqual(fine.warnings, []);
  const def = donnan(78.5, [{ dipole: 0 }, { dipole: 0 }]);
  def.grid = { minCells: 10 };
  const coarse = new Device(def).solve();
  assert.ok(coarse.warnings.some((w) => /interfaces\[0\].*double layer unresolved/.test(w)), coarse.warnings.join('\n'));
});

test('ε = 0 metal region | n-Si: the pinned dipole is the Schottky barrier alignment', () => {
  // Metal: electrons on a fixed positive background, strictly neutral. Its φ is bookkeeping,
  // slaved to μ̄_e. With μ̄_e continuous and the jump pinned at the dipole d, the electron
  // density at the silicon surface is N_c e^{d/V_T} whatever the bias: barrier φ_B = −d.
  const Nm = units.perCm3(1e22), Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19), ND = units.perCm3(1e16);
  const phiB = 0.7;
  const eps = 11.7 * EPS0;
  const dev = new Device({
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: {
      metal: { epsr: 0, species: { 'e-': { D: 1e-2, mu0: 0, cRef: Nm } } },
      Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } },
    },
    regions: [
      { name: 'metal', material: 'metal', length: 10e-9, fixedCharge: Nm * FARADAY, grid: { minCells: 4 } },
      { name: 'n-Si', material: 'Si', length: 2e-6, fixedCharge: ND * FARADAY },
    ],
    interfaces: [{ dipole: -phiB }],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' },
    },
    grid: { hmin: 0.2e-9, hmax: 20e-9, ratio: 1.1 },
  });
  const Vbi = phiB - VT * Math.log(Nc / ND);
  const g = dev.grid, gR = g.regionStart[1];
  for (const V of [0, -0.2, 0.1]) {
    dev.set({ contacts: { right: { V } } }); // metal relative to Si: forward for V < 0 here
    const sol = dev.solve();
    assert.ok(sol.converged, `V=${V}`);
    // The electrode's surface charge sits in the metal's boundary half-box, nudging its electron
    // density (a finite-volume stand-in for the metal's quantum capacitance), hence 0.5%.
    assert.ok(Math.abs(sol.c['e-'][gR] / (Nc * Math.exp(-phiB / VT)) - 1) < 5e-3, `V=${V}: surface density`);
    assert.ok(sol.phi.slice(0, g.regionEnd[0] + 1).every(Number.isFinite), 'metal φ is defined (bookkeeping)');
    // Depletion charge in the silicon (depletion approximation with the −V_T correction).
    let Q = 0;
    for (let k = gR; k <= g.regionEnd[1]; k++) Q += g.vol[k] * FARADAY * (ND - sol.c['e-'][k] + sol.c['h+'][k]);
    const Qd = Math.sqrt(2 * FARADAY * eps * ND * (Vbi + V - VT));
    assert.ok(Math.abs(Q / Qd - 1) < 0.02, `V=${V}: Q ${Q} vs ${Qd}`);
  }
});

test('a cluster cut off by neutral faces with nothing crossing is reported as floating', () => {
  const def = {
    species: [{ name: 'e-', z: -1 }, ...salt],
    materials: { metal: { epsr: 0, species: { 'e-': { D: 1e-2, mu0: 0, cRef: 1e4 } } }, water: water() },
    regions: [
      { material: 'metal', length: 1e-8, fixedCharge: 1e4 * FARADAY, c0: { 'e-': 1e4 } },
      { material: 'water', length: 1e-7 },
    ],
    interfaces: [{ phi: 'neutral' }],
    contacts: { right: bath(10) },
  };
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /regions\[0\].*electrostatically floating/.test(e.message));
});

test('a coarse grid at a depleted electrode is reported, and a graded one is quiet', () => {
  // Near the limiting current, Ag⁺ is depleted at one electrode and its profile there is
  // steeper than one coarse cell can carry: the current comes out too large (exact: tanh).
  const el = (V = 0) => ({ V, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
  const cell = (grid, V) =>
    new Device({
      species: [
        { name: 'Ag+', z: 1, cRef: 1000 },
        { name: 'NO3-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
      regions: [{ material: 'water', length: 20e-6, c0: { 'NO3-': 10 } }],
      contacts: { left: el(), right: el(V) },
      grid,
    }).solve();
  const coarse = cell({ minCells: 20 }, 0.5);
  assert.ok(coarse.warnings.some((w) => /Ag\+: a steep profile.*too large/.test(w)), coarse.warnings.join('\n'));
  for (const V of [0.2, 0.5, 1]) {
    const fine = cell({ hmin: 2e-9, hmax: 1e-6 }, V);
    assert.deepEqual(fine.warnings, [], `V=${V}`);
  }
});

test('ε = 0: a short step converges cleanly across a wide range of concentrations (3 M KCl beside 10 mM ions)', () => {
  // A zinc electrode, ZnSO₄, then a KCl bridge carrying 10 mM of each other ion, then copper. On a
  // 1 µs step the z-weighted balances nearly duplicate the neutrality row; assembled naively,
  // elimination lost ~10 digits here and Newton rattled at ~1e-6 thermal units, never converging.
  const ions = { 'Zn2+': [2, 0.703e-9, -147.06e3], 'Cu2+': [2, 0.714e-9, 65.49e3], 'K+': [1, 1.957e-9, -283.27e3], 'SO42-': [-2, 1.065e-9, -744.53e3], 'Cl-': [-1, 2.032e-9, -131.228e3] };
  const names = Object.keys(ions), tr = 10;
  const c0A = { 'Zn2+': 100, 'SO42-': 100, 'Cu2+': tr, 'K+': tr, 'Cl-': 3 * tr };
  const c0B = { 'K+': 3000, 'Cl-': 3000 + 2 * tr, 'Zn2+': tr, 'Cu2+': tr, 'SO42-': tr };
  const electrode = (metal) => ({
    phi: 'neutral',
    reactions: [{ equation: `${metal}2+ + 2 e- = ${metal}(s)`, fixed: { [`${metal}(s)`]: 0 }, k0: 1e-3, alpha: 0.5 }],
  });
  const dev = new Device({
    species: [...names.map((name) => ({ name, z: ions[name][0], cRef: 1000 })), { name: 'e-', z: -1 }],
    materials: {
      water: { epsr: 0, species: Object.fromEntries(names.map((name) => [name, { D: ions[name][1], mu0: ions[name][2] }])) },
      Zn: { conductor: { species: 'e-', conductivity: 1.7e7 } },
      Cu: { conductor: { species: 'e-', conductivity: 6e7 } },
    },
    regions: [
      { material: 'Zn', length: 1e-3 },
      { name: 'ZnSO4', material: 'water', length: 0.01, c0: c0A },
      { name: 'KCl', material: 'water', length: 0.01, c0: c0B },
      { material: 'Cu', length: 1e-3 },
    ],
    interfaces: [electrode('Zn'), {}, electrode('Cu')],
    contacts: { left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 1.1015, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
    grid: { hmin: 1e-6, hmax: 2e-4 },
  });
  const r = dev.step(1e-6);
  assert.ok(r.converged, JSON.stringify(r.history));
  assert.ok(r.history.at(-1) < 1e-8, `final update ${r.history.at(-1)}`);
});

test('ε = 0: an initial composition that carries net charge, with nothing to neutralise it, is an error', () => {
  const def = (cl) => ({
    species: salt,
    materials: { water: water(0) },
    regions: [{ material: 'water', length: 1e-6, c0: { 'Na+': 100, 'Cl-': cl } }],
    contacts: {
      left: { V: 0, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 },
      right: { V: 0, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 },
    },
  });
  assert.throws(() => new Device(def(90)).solve(), /net charge of 10\.0 mol\/m³.*strictly neutral/);
  assert.ok(new Device(def(100)).solve().converged);
});
