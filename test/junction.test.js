import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { henderson, planck } from '../src/kit.js';

// Liquid-junction potentials, φ_R − φ_L between two neutral solutions, three ways: Henderson's
// formula (concentrations mixing linearly), Planck's (the steady state of constrained diffusion
// between two held solutions), and a free-diffusion junction grown from a sharp boundary. driftlet
// simulates the last two as devices: a steady solve between two baths, and a transient.
//
// The ion sets and conductivities are LJPcalc's (github.com/saglag/LJPcalc, MIT): equivalent
// conductivities from Lange's Handbook (15th ed., Table 8.32), and relative mobilities from
// JPCalc (Barry and Lynch, J. Membr. Biol. 121, 101 (1991)) for gluconate and HEPES, each D =
// λ·RT/(|z|F²) at the set's temperature. The published values: LJPcalc's and JLJP's stationary
// Nernst–Planck solutions, JPCalc's Henderson ones, Ng and Barry's measurements (J. Neurosci.
// Methods 56, 37 (1995)).

const LAMBDA = { // charge, equivalent conductivity (10⁻⁴ m² S/mol)
  K: [1, 73.5], Na: [1, 50.11], Li: [1, 38.69], Cs: [1, 77.3], Zn: [2, 52.8], Mg: [2, 53.06], Ca: [2, 59.5],
  Cl: [-1, 76.31], F: [-1, 55.4], HCO3: [-1, 44.5], H2PO4: [-1, 33], SO4: [-2, 80], Gluconate: [-1, 0.33 * 73.5], HEPES: [-1, 0.3 * 73.5],
};
const ionsAt = (T) => Object.fromEntries(Object.entries(LAMBDA).map(([n, [z, l]]) => [n, { z, D: (l * 1e-4 * GAS_CONSTANT * T) / (Math.abs(z) * FARADAY * FARADAY) }]));
// [name, °C, left (c0: the pipette, say), right (cL: the bath), published value (mV), what it is]
const SETS = [
  ['ZnCl₂ | KCl (JLJP)', 25, { Zn: 9, Cl: 18 }, { Zn: 0.0284, K: 3, Cl: 3.0568 }, -20.79558643, 'stationary'],
  ['NaCl | KCl, 50 mM (Ng and Barry)', 25, { Na: 50, Cl: 50 }, { K: 50, Cl: 50 }, -4.3, 'measured'],
  ['CaCl₂ + MgCl₂ | LiCl (Ng and Barry)', 25, { Ca: 50, Mg: 50, Cl: 200 }, { Li: 100, Cl: 100 }, -8.2, 'measured'],
  ['CsF | NaCl (JPCalc manual)', 20, { Na: 10, Cl: 10, Cs: 135, F: 135 }, { Na: 145, Cl: 145 }, 8.74, 'henderson'],
  ['K-gluconate pipette | bath (Figl et al. 2003)', 20, { K: 145, Na: 13, Mg: 1, HEPES: 5, Gluconate: 145, Cl: 10 }, { K: 2.81, Na: 144.98, Mg: 2, Ca: 1, HEPES: 5, Cl: 148.79 }, 16.05, 'stationary'],
  ['K-gluconate pipette | bath (Figl et al. 2003), JPCalc', 20, { K: 145, Na: 13, Mg: 1, HEPES: 5, Gluconate: 145, Cl: 10 }, { K: 2.81, Na: 144.98, Mg: 2, Ca: 1, HEPES: 5, Cl: 148.79 }, 15.6, 'henderson'],
  ['K-gluconate pipette | ACSF (Walz)', 25, { K: 125, Na: 10, Mg: 1, Cl: 2, Gluconate: 125, HEPES: 10 }, { K: 3, Na: 152.2, Mg: 1.5, Ca: 2.4, Cl: 133.8, H2PO4: 1.2, HCO3: 25, SO4: 1.5 }, 15.729, 'stationary'],
];

// A device: the two solutions as baths, a junction zone between, either held in steady state
// (100 µm of water, ε > 0, quasi-neutral at this scale) or grown from a boundary 1 µm wide
// (strictly neutral water, 400 µm, with the boundary mid-way).
// Ions absent from one side are there as a 1 nM trace (balanced on the most concentrated ion of
// the other sign), which shifts nothing measurable.
function junction(left, right, T, { free = false } = {}) {
  const names = [...new Set([...Object.keys(left), ...Object.keys(right)])];
  const ions = ionsAt(T);
  const side = (s) => {
    const c = Object.fromEntries(names.map((n) => [n, Math.max(s[n] ?? 0, 1e-6)]));
    const q = names.reduce((a, n) => a + ions[n].z * c[n], 0);
    if (q !== 0) {
      const fix = names.filter((n) => Math.sign(ions[n].z) === -Math.sign(q)).sort((a, b) => c[b] - c[a])[0];
      c[fix] += q / -ions[fix].z;
    }
    return c;
  };
  const cl = side(left), cr = side(right), L = free ? 400e-6 : 100e-6, w = 0.5e-6;
  const ref = names.find((n) => ions[n].z < 0);
  const region = { material: 'water', length: L };
  if (free) region.c0 = Object.fromEntries(names.map((n) => [n, { x: [0, L / 2 - w, L / 2 + w, L], values: [cl[n], cl[n], cr[n], cr[n]] }]));
  return new Device({
    T,
    species: names.map((n) => ({ name: n, z: ions[n].z, cRef: 1000 })),
    materials: { water: { epsr: free ? 0 : 78.5, species: Object.fromEntries(names.map((n) => [n, { D: ions[n].D, mu0: 0 }])) } },
    regions: [region],
    contacts: { left: { bath: { c: cl, reference: ref } }, right: { bath: { c: cr, reference: ref }, I: 0 } },
    grid: free ? { hmin: 1e-6, hmax: 1e-6 } : { hmin: 0.05e-6, hmax: 0.05e-6 },
  });
}
const dphi = (s) => s.phi.at(-1) - s.phi[0];

test('Henderson and Planck agree for a single salt, with the diffusion potential (t₊ − t₋)V_T ln(c₁/c₂), and every junction does', () => {
  const T = 298.15, ions = ionsAt(T), VT = (GAS_CONSTANT * T) / FARADAY;
  const tp = ions.Na.D / (ions.Na.D + ions.Cl.D), expected = (2 * tp - 1) * VT * Math.log(100 / 10);
  assert.ok(Math.abs(henderson({ Na: 100, Cl: 100 }, { Na: 10, Cl: 10 }, ions) / expected - 1) < 1e-12);
  assert.ok(Math.abs(planck({ Na: 100, Cl: 100 }, { Na: 10, Cl: 10 }, ions) / expected - 1) < 1e-9);
  assert.ok(Math.abs(dphi(junction({ Na: 100, Cl: 100 }, { Na: 10, Cl: 10 }, T).solve()) / expected - 1) < 1e-6);
});

test("Henderson's formula is JPCalc's, and Planck's (by shooting) is driftlet's steady junction and the stationary Nernst–Planck codes'", () => {
  for (const [name, C, left, right, published, kind] of SETS) {
    const T = 273.15 + C, ions = ionsAt(T);
    const H = 1000 * henderson(left, right, ions, { T }), P = 1000 * planck(left, right, ions, { T });
    const steady = junction(left, right, T).solve();
    assert.ok(steady.converged);
    const S = 1000 * dphi(steady);
    // Two independent ways to the same steady state: finite volumes and shooting.
    assert.ok(Math.abs(S - P) < 2e-4, `${name}: driftlet ${S} vs Planck ${P} mV`);
    if (kind === 'henderson') assert.ok(Math.abs(H - published) < 0.015, `${name}: Henderson ${H} vs JPCalc's ${published} mV`);
    // LJPcalc and JLJP solve the same equations; the ZnCl₂ set's 300-fold fall in Zn²⁺ leaves
    // theirs 0.02 mV off (driftlet's converges with the grid, and shooting agrees with it).
    if (kind === 'stationary') assert.ok(Math.abs(P - published) < 0.02, `${name}: Planck ${P} vs ${published} mV`);
    if (kind === 'measured') assert.ok(Math.abs(P - published) < 0.25 && Math.abs(H - published) < 0.25, `${name}: ${P}, ${H} vs measured ${published} mV`);
  }
});

test('a free-diffusion junction grown from a sharp boundary holds a constant potential, between Henderson\'s and Planck\'s for these mixtures', () => {
  for (const [name, C, left, right] of SETS.filter((s) => /JLJP|CsF/.test(s[0]))) {
    const T = 273.15 + C, ions = ionsAt(T), dev = junction(left, right, T, { free: true });
    const E = [0.1, 0.5].map((t) => {
      const s = dev.advance(t, { tol: 1e-4 });
      assert.ok(s.converged);
      return 1000 * dphi(s);
    });
    // Self-similar (√(Dt) from 14 to 30 µm, against a boundary 1 µm wide and baths 200 µm away).
    assert.ok(Math.abs(E[1] - E[0]) < 5e-3, `${name}: ${E.join(', ')} mV`);
    const H = 1000 * henderson(left, right, ions, { T }), P = 1000 * planck(left, right, ions, { T });
    const [lo, hi] = [Math.min(H, P), Math.max(H, P)];
    assert.ok(E[1] > lo + 0.01 && E[1] < hi - 0.01, `${name}: free ${E[1]} not between Henderson ${H} and Planck ${P}`);
  }
});

test('junction formulas check their inputs', () => {
  const ions = ionsAt(298.15);
  assert.throws(() => henderson({ Na: 1, Cl: 2 }, { Na: 1, Cl: 1 }, ions), (e) => e instanceof DeviceError && /left solution carries a net charge/.test(e.message));
  assert.throws(() => planck({ Na: 1, Cl: 1 }, { Rb: 1, Cl: 1 }, ions), /no \{ z, D \} for 'Rb'/);
});
