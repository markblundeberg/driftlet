import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, units } from '../src/index.js';
import { vacuumDipole, vacuumLevel, vacuumZeroCharge } from '../src/kit.js';

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

const throwsDevice = (fn, pattern) => assert.throws(fn, (e) => e instanceof DeviceError && pattern.test(e.message), String(pattern));

test('vacuum-level alignment: anchors and offsets reproduce Anderson, ionisation energies and gates', () => {
  // Anderson: electron affinities on the e⁻ standard level (the conduction band) give
  // ΔE_c = χ_A − χ_B, the same dipole as the equivalent step.
  const def = hetero();
  const chiA = 4.07, chiB = 3.8;
  const dv = vacuumDipole(def, { material: 'A', anchor: 'e-', offset: chiA }, { material: 'B', anchor: 'e-', offset: chiB });
  const ds = new Device({ ...def, interfaces: [{ step: { species: 'e-', value: units.eV(chiA - chiB) } }] }).model.interfaces[0].dipole;
  assert.ok(Math.abs(dv - ds) < 1e-12, `${dv} vs ${ds}`);
  assert.equal(new Device({ ...def, interfaces: [{ dipole: dv }] }).model.interfaces[0].dipole, dv);
  // The same alignment through ionisation energies on the valence band (h⁺ standard level).
  const egA = 1.42, egB = 1.8;
  const di = vacuumDipole(def, { material: 'A', anchor: 'h+', offset: chiA + egA }, { material: 'B', anchor: 'h+', offset: chiB + egB });
  assert.ok(Math.abs(di - ds) < 1e-12);
  // φ as the anchor: the offsets are surface potentials χ (V_vac = φ − χ), so φ_R − φ_L = χ_R − χ_L.
  assert.ok(Math.abs(vacuumDipole(def, { material: 'A', anchor: 'phi', offset: 0.3 }, { material: 'B', anchor: 'phi', offset: 0.1 }) + 0.2) < 1e-12);
  // A gate (work function W) on material A: zeroCharge = W − χ − μ°_e/F, the flat-band form.
  const zc = vacuumZeroCharge(def, 4.5, { material: 'A', anchor: 'e-', offset: chiA });
  assert.ok(Math.abs(zc - (4.5 - chiA - def.materials.A.species['e-'].mu0 / FARADAY)) < 1e-12);
  assert.equal(vacuumLevel(def, { material: 'A', anchor: 'phi', offset: 4 }), -4);
});

test('vacuum helpers check their inputs', () => {
  const def = hetero();
  throwsDevice(() => vacuumDipole(def, { material: 'A', anchor: 'X', offset: 1 }, { material: 'B', anchor: 'e-', offset: 1 }), /anchor must be 'phi' or a charged species/);
  throwsDevice(() => vacuumDipole(def, { material: 'C', anchor: 'e-', offset: 1 }, { material: 'B', anchor: 'e-', offset: 1 }), /unknown material/);
  throwsDevice(() => vacuumZeroCharge(def, NaN, { material: 'A', anchor: 'e-', offset: 1 }), /workFunction/);
  def.materials.Au = { conductor: { species: 'e-', conductivity: 4e7 } };
  throwsDevice(() => vacuumLevel(def, { material: 'Au', anchor: 'e-', offset: 5 }), /conductor, with no φ/);
});

test('the kit and plots refuse an option they do not take, and suggest the one meant', async () => {
  const { bandDiagram } = await import('../src/plot.js');
  const { pulse, photogeneration, semiconductor, build, layer, ohmic } = await import('../src/kit.js');
  const sol = new Device(build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donors: 16.6 }), ohmic(0)], grid: { hmin: 1e-8, hmax: 1e-7 } })).solve();
  assert.throws(() => bandDiagram(sol, { engery: true }), /bandDiagram: no option 'engery' \(did you mean 'energy'\?\)/);
  assert.throws(() => bandDiagram(sol, { Range: [0, 1] }), /did you mean 'range'/);
  assert.throws(() => bandDiagram(sol, { levels: 'x' }), /levels must be a list/);
  assert.throws(() => pulse({ start: 0, width: 1, amplitude: 1, rise_time: 1 }), /pulse: no option 'rise_time'/);
  assert.throws(() => photogeneration({ material: 'Si', flux: 1, alpha: 1, mu: 1, form: 0 }), /did you mean 'from'/);
  assert.throws(() => semiconductor('Si', { T: 350 }), /data are for 300 K/);
  assert.throws(() => build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donor: 16.6 }), ohmic(0)] }), /stack\[1\]: no option 'donor' \(did you mean 'donors'\?\)/);
});
