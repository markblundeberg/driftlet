import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor } from '../src/kit.js';

// Auger recombination as a rate law: r = (C_n n + C_p p) n p (1 − e^{−A/RT}), which is
// (C_n n + C_p p)(n p − n_i²) for Boltzmann carriers.

const T = 300, RT = GAS_CONSTANT * T;
const Si = semiconductor('Si');
const { 'e-': el, 'h+': ho } = Si.materials.Si.species;
const ni2 = el.cRef * ho.cRef * Math.exp(-(ho.mu0 + el.mu0) / RT);
const Cn = units.cm6PerS(2.8e-31), Cp = units.cm6PerS(9.9e-32); // silicon's (Dziewior & Schmid 1977)

const root = (f, lo, hi) => {
  for (let k = 0; k < 200; k++) {
    const m = 0.5 * (lo + hi);
    if (f(lo) * f(m) <= 0) hi = m;
    else lo = m;
  }
  return 0.5 * (lo + hi);
};

// An n-type slab lit uniformly, its holes blocked at both contacts: no gradients, so at every node
// G equals the recombination, and neutrality gives n = N_D + p.
const slab = (ND, G, reactions) =>
  new Device(
    build({
      T,
      library: [Si],
      stack: [ohmic(0, ['e-']), layer('Si', 10e-6, { donors: ND, c0: { 'h+': 1e-12 } }), ohmic(0, ['e-'])],
      bulkReactions: [...reactions, { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: G } }],
    }),
  ).solve();

test('Auger under uniform light: G = (C_n n + C_p p)(n p − n_i²) with n = N_D + p, low injection to high', () => {
  const ND = units.perCm3(1e18);
  for (const G of [units.perCm3(1e20), units.perCm3(1e24), units.perCm3(1e27)]) {
    const s = slab(ND, G, [{ equation: 'e- + h+ = 0', auger: { Si: { Cn, Cp } } }]);
    assert.ok(s.converged);
    const R = (p) => {
      const n = ND + p;
      return (Cn * n + Cp * p) * (n * p - ni2);
    };
    const p = root((q) => R(q) - G, 0, 1e4 * ND);
    const mid = s.x.length >> 1;
    assert.ok(Math.abs(s.c['h+'][mid] / p - 1) < 1e-6, `G ${G}: p ${s.c['h+'][mid]} vs ${p}`);
    const [auger, light] = s.bulkReactions;
    assert.ok(Math.abs(auger.rate[mid] / G - 1) < 1e-6);
    assert.ok(Math.abs(auger.total / light.total - 1) < 1e-9);
  }
  // At low injection in n⁺ silicon, the minority lifetime is 1/(C_n N_D²): 3.6 µs at 1e18 cm⁻³.
  const G = units.perCm3(1e18), s = slab(ND, G, [{ equation: 'e- + h+ = 0', auger: { Si: { Cn, Cp } } }]);
  const tau = s.c['h+'][s.x.length >> 1] / G;
  assert.ok(Math.abs(tau * Cn * ND * ND - 1) < 1e-3, `τ ${tau} vs ${1 / (Cn * ND * ND)}`);
});

test('Auger beside SRH and radiative recombination: the three rates add, each its own law', () => {
  const ND = units.perCm3(1e17), G = units.perCm3(1e25), tau = 1e-6, B = units.cm3PerS(1e-14);
  const s = slab(ND, G, [
    { equation: 'e- + h+ = 0', auger: { Si: { Cn, Cp } } },
    { equation: 'e- + h+ = 0', srh: { Si: { tauN: tau, tauP: tau } } },
    { equation: 'e- + h+ = 0', kf: { Si: B } },
  ]);
  assert.ok(s.converged);
  const ni = Math.sqrt(ni2);
  const R = (p) => {
    const n = ND + p, x = n * p - ni2;
    return (Cn * n + Cp * p) * x + x / (tau * (n + ni) + tau * (p + ni)) + B * x;
  };
  const p = root((q) => R(q) - G, 0, 1e4 * ND), mid = s.x.length >> 1;
  assert.ok(Math.abs(s.c['h+'][mid] / p - 1) < 1e-6, `p ${s.c['h+'][mid]} vs ${p}`);
  const n = ND + p, x = n * p - ni2;
  const [a, srh, rad] = s.bulkReactions;
  assert.ok(Math.abs(a.rate[mid] / ((Cn * n + Cp * p) * x) - 1) < 1e-6);
  assert.ok(Math.abs(rad.rate[mid] / (B * x) - 1) < 1e-6);
  assert.ok(Math.abs((a.total + srh.total + rad.total) / s.bulkReactions[3].total - 1) < 1e-9);
});

test('Auger in the dark: no rate in equilibrium, a pn diode under bias recombining what its current carries', () => {
  const diode = (V) =>
    new Device(
      build({
        T,
        library: [Si],
        stack: [ohmic(V), layer('Si', 2e-6, { donors: units.perCm3(1e19) }), { dipole: 0 }, layer('Si', 2e-6, { acceptors: units.perCm3(1e18) }), ohmic(0)],
        bulkReactions: [{ equation: 'e- + h+ = 0', auger: { Si: { Cn, Cp } } }],
      }),
    ).solve();
  const eq = diode(0);
  assert.ok(eq.converged);
  const [r0] = eq.bulkReactions;
  assert.ok(Math.max(...r0.rate.map(Math.abs)) < 1e-9 * Math.max(...r0.forward), `rate at rest ${Math.max(...r0.rate.map(Math.abs))}`);
  // Forward-biased (the n side's level per charge pulled down): the current is what Auger
  // recombines plus the minority carriers the ohmic contacts take, so 0 < R·F ≤ |I|.
  const on = diode(-0.6);
  assert.ok(on.converged);
  const [r] = on.bulkReactions;
  assert.ok(r.total > 0 && r.total * 96485.33212 <= Math.abs(on.current) * (1 + 1e-9), `R·F ${r.total * 96485.33212} vs I ${on.current}`);
});

test('auger: refuses a reaction that isn\'t e⁻ + h⁺ = 0, both coefficients zero, and a second rate law', () => {
  const def = (rx) => build({ T, library: [Si], stack: [ohmic(0), layer('Si', 1e-6, { donors: units.perCm3(1e16) }), ohmic(0)], bulkReactions: [rx] });
  assert.throws(() => new Device(def({ equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, auger: { Si: { Cn } } })), DeviceError);
  assert.throws(() => new Device(def({ equation: 'e- + h+ = 0', auger: { Si: { Cn: 0, Cp: 0 } } })), /Cn or Cp/);
  assert.throws(() => new Device(def({ equation: 'e- + h+ = 0', auger: { Si: { Cn } }, kf: { Si: 1 } })), /one of them/);
  assert.throws(() => new Device(def({ equation: 'e- + h+ = 0', auger: { Si: { Cn, tau: 1 } } })), DeviceError);
});

test('units: cm³/s and cm⁶/s per particle into per-mol rate constants', () => {
  assert.ok(Math.abs(units.cm3PerS(1e-14) / (1e-20 * 6.02214076e23) - 1) < 1e-15);
  assert.ok(Math.abs(units.cm6PerS(1e-30) / (1e-42 * 6.02214076e23 * 6.02214076e23) - 1) < 1e-15);
});
