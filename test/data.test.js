import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, DeviceError, units } from '../src/index.js';
import { IONS, H2O, SEMICONDUCTORS, METALS, aqueous, semiconductor, metal, trap, build, layer, ohmic, bath, describe, check } from '../src/kit.js';

// The data library, checked against independent tables, so a typo in one shows up as a
// disagreement with the other.

const RT = GAS_CONSTANT * 298.15;

test('ions: μ° reproduce the electrochemical series, E° = −Σν μ°/(nF)', () => {
  // Standard potentials vs SHE (CRC Handbook, "Electrochemical series"), each from the ion's
  // μ° alone: the metal and Cl₂ are elements in their standard states (μ = 0).
  const series = [
    ['Li+', -3.0401],
    ['Na+', -2.71],
    ['K+', -2.931],
    ['Ag+', 0.7996],
    ['Cu2+', 0.3419],
    ['Zn2+', -0.7618],
  ];
  for (const [ion, E] of series) {
    const { z, mu0 } = IONS[ion];
    assert.ok(Math.abs(mu0 / (z * FARADAY) - E) < 0.01, `${ion}: ${(mu0 / (z * FARADAY)).toFixed(4)} V vs ${E}`);
  }
  // Cl₂ + 2e⁻ ⇌ 2Cl⁻ and Fe³⁺ + e⁻ ⇌ Fe²⁺.
  assert.ok(Math.abs(-IONS['Cl-'].mu0 / FARADAY - 1.35827) < 0.01);
  assert.ok(Math.abs((IONS['Fe3+'].mu0 - IONS['Fe2+'].mu0) / FARADAY - 0.771) < 0.01);
  // H₂O ⇌ H⁺ + OH⁻: pKw = 14.00 at 25 °C.
  const pKw = (IONS['H+'].mu0 + IONS['OH-'].mu0 - H2O) / (RT * Math.LN10);
  assert.ok(Math.abs(pKw - 14.0) < 0.01, `pKw ${pKw}`);
});

test('ions: D reproduce the limiting conductivities by Nernst–Einstein, λ = z²F²D/RT', () => {
  // Equivalent limiting conductivities λ/|z|, S·cm²/mol (CRC Handbook).
  const lambda = { 'H+': 349.65, 'Li+': 38.66, 'Na+': 50.08, 'K+': 73.48, 'Ag+': 61.9, 'Cu2+': 53.6, 'Zn2+': 52.8, 'Fe2+': 54, 'Fe3+': 68, 'OH-': 198, 'Cl-': 76.31, 'NO3-': 71.42, 'SO42-': 80.0 };
  assert.deepEqual(Object.keys(lambda).sort(), Object.keys(IONS).sort());
  for (const [ion, l] of Object.entries(lambda)) {
    const { z, D } = IONS[ion];
    const fromD = (Math.abs(z) * FARADAY * FARADAY * D) / RT / 1e-4;
    assert.ok(Math.abs(fromD / l - 1) < 0.005, `${ion}: ${fromD.toFixed(2)} vs ${l}`);
  }
});

test('semiconductors: the band data reproduce the quoted intrinsic concentrations (Si excepted, see docs)', () => {
  const kT = (GAS_CONSTANT * 300) / FARADAY;
  const ni = ({ Nc, Nv, Eg }) => Math.sqrt(Nc * Nv) * Math.exp(-Eg / (2 * kT));
  assert.ok(Math.abs(ni(SEMICONDUCTORS.Ge) / 2.4e13 - 1) < 0.1);
  assert.ok(Math.abs(ni(SEMICONDUCTORS.GaAs) / 1.79e6 - 1) < 0.15);
  // Si: Sze's Nc, Nv and Eg imply 6.7e9 cm⁻³, below his quoted 1.45e10 and the measured 9.65e9.
  assert.ok(Math.abs(ni(SEMICONDUCTORS.Si) / 6.7e9 - 1) < 0.02);
});

test('library pieces build devices: a silicon diode, a salt bridge, metal electrodes', () => {
  const si = semiconductor('Si');
  assert.ok(Math.abs(si.materials.Si.species['e-'].D - 1500e-4 * (GAS_CONSTANT * 300) / FARADAY) < 1e-12);
  const diode = build({
    T: 300,
    library: [si],
    stack: [ohmic(0), layer('Si', 1e-6, { donors: units.perCm3(1e17) }), layer('Si', 1e-6, { acceptors: units.perCm3(1e16) }), ohmic(0)],
    grid: { hmin: 1e-9, hmax: 20e-9 },
  });
  assert.ok(new Device(diode).solve().converged);
  assert.doesNotMatch(describe(diode), /warnings/);

  const salt = build({
    library: [aqueous(['Na+', 'Cl-'], { epsr: 0 })],
    stack: [bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-'), layer('water', 1e-5), bath({ 'Na+': 10, 'Cl-': 10 }, 'Cl-', { I: 0 })],
  });
  assert.equal(salt.materials.water.epsr, 0);
  assert.ok(new Device(salt).solve().converged);
  assert.deepEqual(metal('Cu'), { species: [{ name: 'e-', z: -1 }], materials: { Cu: { conductor: { species: 'e-', conductivity: METALS.Cu.conductivity } } } });
  assert.throws(() => aqueous(['Na+', 'Br-']), (e) => e instanceof DeviceError && /no data for 'Br-'/.test(e.message));
  assert.throws(() => semiconductor('InP'), /no data for 'InP'/);
  assert.ok(Object.isFrozen(IONS['Na+']));
});

// Semiconductors at another temperature, semiconductor(name, { T }).

// The table's numbers, as semiconductor() built them before it took a T.
function at300(name) {
  const sc = SEMICONDUCTORS[name], VT = (GAS_CONSTANT * 300) / FARADAY, perCm3 = (n) => units.perCm3(n);
  return {
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: {
      [name]: {
        T: 300,
        epsr: sc.epsr,
        species: {
          'e-': { D: sc.mun * 1e-4 * VT, mu0: 0, cRef: perCm3(sc.Nc) },
          'h+': { D: sc.mup * 1e-4 * VT, mu0: sc.Eg * FARADAY, cRef: perCm3(sc.Nv) },
        },
      },
    },
  };
}
// n_i from a library piece: n_i² = N_c N_v e^(−E_g/kT), in the material's own terms (mol/m³).
const niOf = (lib, T) => {
  const { 'e-': e, 'h+': h } = Object.values(lib.materials)[0].species;
  return Math.sqrt(e.cRef * h.cRef * Math.exp(-(h.mu0 - e.mu0) / (GAS_CONSTANT * T)));
};

test('semiconductors at T: 300 K (or no T) is the table, exactly; outside 150–500 K is refused', () => {
  for (const name of Object.keys(SEMICONDUCTORS)) {
    assert.deepStrictEqual(semiconductor(name), at300(name));
    assert.deepStrictEqual(semiconductor(name, { T: 300 }), at300(name));
  }
  for (const T of [149, 501, NaN]) {
    assert.throws(() => semiconductor('Si', { T }), (e) => e instanceof DeviceError && /T must be the temperature in K, from 150 to 500/.test(e.message));
  }
  // (a string is named as one, not shown as if it were the number)
  assert.throws(() => semiconductor('Si', { T: '350' }), (e) => e instanceof DeviceError && e.message.includes('T must be a number, the temperature in K, got "350" (a string)'));
  assert.ok(semiconductor('Si', { T: 150 }) && semiconductor('GaAs', { T: 500 }));
});

test('semiconductors at T: the gap by Varshni, N_c and N_v as T^1.5, the mobilities as lattice-limited power laws', () => {
  // Silicon's n_i grows 8.3%/K at 298 K: d ln n_i/dT = 3/(2T) + E_g0/(2kT²), with E_g0 the gap
  // extrapolated linearly to 0 K from its slope at T (Varshni, α = 4.73e-4 eV/K, β = 636 K).
  const T = 298, dT = 0.5, VT = (GAS_CONSTANT * T) / FARADAY;
  const ni = (T) => niOf(semiconductor('Si', { T }), T);
  const growth = Math.log(ni(T + dT) / ni(T - dT)) / (2 * dT);
  const slope = -(4.73e-4 * T * (T + 2 * 636)) / (T + 636) ** 2; // dE_g/dT, eV/K
  const Eg = 1.12 - ((4.73e-4 * T * T) / (T + 636) - (4.73e-4 * 300 * 300) / 936);
  const expected = 1.5 / T + (Eg - T * slope) / (2 * VT * T);
  assert.ok(Math.abs(growth / expected - 1) < 1e-4, `${growth} vs ${expected}`);
  assert.ok(Math.abs(growth - 0.083) < 0.002, `n_i grows ${(100 * growth).toFixed(2)}%/K`);
  // At 400 K: the gap narrows by Varshni's shift from 300 K, and each piece moves by its power.
  const t = 400 / 300, kT = (GAS_CONSTANT * 400) / FARADAY;
  for (const [name, [Eg0, a, b], gn, gp] of [['Si', [1.17, 4.73e-4, 636], 2.42, 2.2], ['Ge', [0.7437, 4.774e-4, 235], 1.66, 2.33], ['GaAs', [1.519, 5.405e-4, 204], 1.0, 2.1]]) {
    const sc = SEMICONDUCTORS[name], m = semiconductor(name, { T: 400 }).materials[name].species;
    const varshni = (T) => Eg0 - (a * T * T) / (T + b);
    const close = (x, y, what) => assert.ok(Math.abs(x / y - 1) < 1e-12, `${name} ${what}: ${x} vs ${y}`);
    close(m['h+'].mu0 / FARADAY, sc.Eg + varshni(400) - varshni(300), 'E_g');
    close(m['e-'].cRef, units.perCm3(sc.Nc) * t ** 1.5, 'N_c');
    close(m['h+'].cRef, units.perCm3(sc.Nv) * t ** 1.5, 'N_v');
    close(m['e-'].D, sc.mun * 1e-4 * t ** -gn * kT, 'D_n = μ_n kT/q');
    close(m['h+'].D, sc.mup * 1e-4 * t ** -gp * kT, 'D_p = μ_p kT/q');
  }
});

test('a material carries the T its data are for: a device at another T warns (and after set({ T })), hand-written ones don\'t', () => {
  const T = (w) => w.filter((m) => /its data are for/.test(m));
  const slab = (lib, more = {}) => build({ ...more, library: [lib], stack: [ohmic(0), layer('Si', 1e-6), ohmic(0)] });
  // Matched: silent. Every semiconductor() carries its T, the default 300 K too.
  assert.equal(semiconductor('Si').materials.Si.T, 300);
  assert.equal(semiconductor('Si', { T: 350 }).materials.Si.T, 350);
  assert.deepEqual(T(new Device(slab(semiconductor('Si'), { T: 300 })).solution().warnings), []);
  assert.deepEqual(T(new Device(slab(semiconductor('Si', { T: 350 }), { T: 350 })).solution().warnings), []);
  // A device at 350 K with 300 K data: intrinsic n_i 36% low, and now it says so, with the cure.
  const dev = new Device(slab(semiconductor('Si'), { T: 350 }));
  const sol = dev.solve();
  const [w] = T(sol.warnings);
  assert.ok(w && w.includes('materials.Si') && w.includes('300 K') && w.includes('350 K') && w.includes('T: 300') && w.includes("semiconductor(name, { T: 350 })"), w);
  const said = check(dev, sol, { refine: false }).items.find((c) => c.name === 'warnings');
  assert.ok(said.ok === null && said.summary.includes('materials.Si: its data are for 300 K'), said.summary);
  assert.match(describe(slab(semiconductor('Si'), { T: 350 })), /warnings:\n {2}materials\.Si: its data are for 300 K/);
  // No T on the device: the default 298.15 K, which the warning names as the default.
  assert.match(T(new Device(slab(semiconductor('Si'))).solution().warnings)[0], /298\.15 K \(the default, as it gives no T\)/);
  // set({ T }) alone moves RT, not the data: it warns; giving the material for T too cures it,
  // and the intrinsic density is then the 350 K data's n_i.
  const live = new Device(slab(semiconductor('Si'), { T: 300 }));
  live.solve();
  live.set({ T: 350 });
  assert.equal(T(live.solve().warnings).length, 1);
  live.set({ T: 350, materials: semiconductor('Si', { T: 350 }).materials });
  const hot = live.solve();
  assert.ok(hot.converged);
  assert.deepEqual(T(hot.warnings), []);
  const ni = niOf(semiconductor('Si', { T: 350 }), 350), mid = hot.x.findIndex((x) => x >= 0.5e-6);
  assert.ok(Math.abs(hot.c['e-'][mid] / ni - 1) < 1e-3, `n ${hot.c['e-'][mid]} vs n_i ${ni}`);
  // Two materials made at different temperatures: the one that differs from the device's warns,
  // and only where a region uses it.
  const two = (stack) => build({ T: 300, library: [semiconductor('Si'), semiconductor('Si', { material: 'hot', T: 350 })], stack });
  const both = T(new Device(two([ohmic(0), layer('Si', 1e-6), { step: { species: 'e-', value: 0 } }, layer('hot', 1e-6), ohmic(0)])).solution().warnings);
  assert.equal(both.length, 1);
  assert.match(both[0], /^materials\.hot: its data are for 350 K, but the device is at 300 K/);
  assert.deepEqual(T(new Device(two([ohmic(0), layer('Si', 1e-6), ohmic(0)])).solution().warnings), []);
  // Hand-written materials say nothing of their temperature: no warning (it's unknown).
  const { T: _, ...plain } = semiconductor('Si').materials.Si;
  assert.deepEqual(T(new Device(slab({ species: semiconductor('Si').species, materials: { Si: plain } }, { T: 350 })).solution().warnings), []);
  // A material's T is checked like any number, and a trap takes the material's T.
  assert.throws(() => new Device(slab({ ...semiconductor('Si'), materials: { Si: { ...plain, T: NaN } } })), (e) => e instanceof DeviceError && /materials\.Si\.T .* got NaN/.test(e.message));
  const hotSi = semiconductor('Si', { T: 350 }), n = niOf(hotSi, 350);
  const tr = trap(hotSi, { density: 1, tauN: 1e-6, tauP: 1e-6 });
  assert.equal(tr.library.materials.Si.T, 350);
  // (midgap at 350 K: X⁻'s μ° sits an electron's μ at n_i(350 K) above X⁰'s)
  const e = hotSi.materials.Si.species['e-'];
  assert.ok(Math.abs(tr.library.materials.Si.species['X-'].mu0 - (e.mu0 + GAS_CONSTANT * 350 * Math.log(n / e.cRef))) < 1e-9);
  assert.throws(() => trap(hotSi, { density: 1, tauN: 1e-6, tauP: 1e-6, T: 300 }), (e) => e instanceof DeviceError && /material 'Si' is for 350 K/.test(e.message));
});

test('temperature errors show the value given: NaN and Infinity as themselves, a string as a string', () => {
  const def = build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6), ohmic(0)] });
  for (const [T, shown] of [[NaN, 'NaN'], [Infinity, 'Infinity'], [-1, '-1'], ['350', '"350" (a string)']]) {
    assert.throws(() => new Device({ ...def, T }), (e) => e instanceof DeviceError && e.message === `T (the temperature, K) must be a positive number, got ${shown}`, `Device T ${String(T)}`);
    assert.throws(() => new Device(def).set({ T }), (e) => e instanceof DeviceError && e.message.endsWith(`got ${shown}`), `set T ${String(T)}`);
    assert.throws(() => build({ T, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6), ohmic(0)] }), (e) => e instanceof DeviceError && e.message === `build.T must be the temperature in K, a positive number, got ${shown}`, `build T ${String(T)}`);
  }
});

test('a silicon cell under light: dV_oc/dT ≈ −2.1 mV/K, as Green\'s −(E_g0/q − V_oc + γkT/q)/T', () => {
  // A 200 µm p base (10¹⁶ cm⁻³, τ = 100 µs) under a 0.3 µm n⁺ emitter, contacts on both faces,
  // light making pairs uniformly. V_oc = (kT/q) ln(J_sc/J_0), and J_0 ∝ T^γ e^(−E_g0/kT), with
  // E_g0 the gap extrapolated linearly to 0 K (1.196 eV here), gives Green's
  // dV_oc/dT = −(E_g0/q − V_oc + γkT/q)/T (M. A. Green, Prog. Photovolt. 11, 333, 2003). γ is 3,
  // from n_i² ∝ N_c N_v, plus the power of T in the minority carriers' transport: D ∝ μT ∝
  // T^−1.42 in the base (T^−1.2 in the emitter), so J_0 ∝ n_i² D/W in a short base gives
  // γ = 1.58, and n_i² √(D/τ) in a long one 2.29; 1.5 ≤ γ ≤ 3 bounds both, and the emitter's
  // share. Without N_c, N_v ∝ T^1.5, γ would fall by 3, and the slope by 0.26 mV/K, outside.
  const NA = units.perCm3(1e16), tau = 100e-6;
  const cell = (T, right) =>
    build({
      T,
      library: [semiconductor('Si', { T })],
      stack: [ohmic(0), layer('Si', units.um(0.3), { donors: units.perCm3(1e19) }), layer('Si', units.um(200), { acceptors: NA }), ohmic(right)],
      bulkReactions: [
        { equation: 'e- + h+ = 0', kf: { Si: 1 / (tau * NA) } },
        { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: 20 } }, // G, mol/(m³·s)
      ],
      grid: { hmin: units.nm(1), hmax: units.um(2) },
    });
  const Voc = (T) => new Device(cell(T, { I: 0 })).solve().terminals.right.V;
  const Jsc = (T) => new Device(cell(T, 0)).solve().current;
  const T = 300, dT = 5, VT = (GAS_CONSTANT * T) / FARADAY;
  const V = Voc(T), dVdT = (Voc(T + dT) - Voc(T - dT)) / (2 * dT);
  const slope = -(4.73e-4 * T * (T + 2 * 636)) / (T + 636) ** 2, Eg0 = 1.12 - T * slope;
  const green = (gamma) => -(Eg0 - V + gamma * VT) / T;
  assert.ok(dVdT < green(1.5) && dVdT > green(3), `dV_oc/dT = ${(1e3 * dVdT).toFixed(3)} mV/K, Green's ${(1e3 * green(1.5)).toFixed(3)} to ${(1e3 * green(3)).toFixed(3)}`);
  // Sharper: the γ this cell implies, with J_sc's own small drift taken out, against J_0's two
  // parts. The base's ∝ n_i²(D/W + W/3τ) (coth for W/L = 0.32), γ = 3 − 1.42·0.967 = 1.63; the
  // transparent emitter's ∝ n_i² D_p/(N_D W_E), γ = 1.8, a fifth of J_0: γ ≈ 1.66 in all.
  const dlnJ = Math.log(Jsc(T + dT) / Jsc(T - dT)) / (2 * dT);
  const gamma = (-(dVdT - VT * dlnJ) * T - Eg0 + V) / VT;
  assert.ok(Math.abs(gamma - 1.66) < 0.1, `γ = ${gamma.toFixed(3)}`);
  assert.ok(Math.abs(V - 0.598) < 0.01, `V_oc = ${V}`);
});
