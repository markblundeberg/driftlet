import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, DeviceError, units } from '../src/index.js';
import { IONS, H2O, SEMICONDUCTORS, METALS, aqueous, semiconductor, metal, build, layer, ohmic, bath, describe } from '../src/kit.js';

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
