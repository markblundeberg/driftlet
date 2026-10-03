import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor } from '../src/kit.js';

// SRH (trap-assisted) kinetics as a rate law: r = n p (1 − e^{−A/RT}) / (τ_p (n + n₁) + τ_n (p + p₁)),
// with p₁ = n p e^{−A/RT}/n₁ from the state, so equilibrium (A = 0) stays exact.

const T = 300, RT = GAS_CONSTANT * T;
const Si = semiconductor('Si');
const { 'e-': el, 'h+': ho } = Si.materials.Si.species;
const ni2 = el.cRef * ho.cRef * Math.exp(-(ho.mu0 - el.mu0) / RT);

// The root of f on [lo, hi] by bisection (f changes sign there).
const root = (f, lo, hi) => {
  for (let k = 0; k < 200; k++) {
    const m = 0.5 * (lo + hi);
    if (f(lo) * f(m) <= 0) hi = m;
    else lo = m;
  }
  return 0.5 * (lo + hi);
};

test('bulk SRH under uniform light: G = R_SRH(n, p) with n = N_D + p, for a midgap trap and a shallow one, low injection to high', () => {
  // An n-type slab lit uniformly, its holes blocked at both contacts: no gradients, so at every
  // node the generation G equals the SRH rate, and neutrality gives n = N_D + p.
  const ND = units.perCm3(1e15), tauN = 2e-6, tauP = 1e-6;
  for (const n1 of [undefined, units.perCm3(1e13)]) {
    for (const G of [units.perCm3(1e19), units.perCm3(1e22)]) {
      const s = new Device(
        build({
          T,
          library: [Si],
          stack: [ohmic(0, ['e-']), layer('Si', 10e-6, { donors: ND, c0: { 'h+': 1e-12 } }), ohmic(0, ['e-'])],
          bulkReactions: [
            { equation: 'e- + h+ = 0', srh: { Si: { tauN, tauP, ...(n1 === undefined ? {} : { n1 }) } } },
            { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: G } },
          ],
        }),
      ).solve();
      assert.ok(s.converged);
      const k1 = n1 ?? Math.sqrt(ni2), p1 = ni2 / k1;
      const R = (p) => {
        const n = ND + p;
        return (n * p - ni2) / (tauP * (n + k1) + tauN * (p + p1));
      };
      const p = root((q) => R(q) - G, 0, 1e3 * (ND + G * tauP));
      const mid = s.x.length >> 1;
      assert.ok(Math.abs(s.c['h+'][mid] / p - 1) < 1e-6, `n1 ${n1}, G ${G}: p ${s.c['h+'][mid]} vs ${p}`);
      // The solution reports each bulk reaction's rate: SRH matching G at every node, and in all.
      const [srh, light] = s.bulkReactions;
      assert.ok(Math.abs(srh.rate[mid] / G - 1) < 1e-6);
      assert.ok(Math.abs(srh.total / light.total - 1) < 1e-9, `R ${srh.total} vs G ${light.total} in all`);
      assert.ok(Math.abs(light.total / (G * 10e-6) - 1) < 1e-12);
      assert.ok(Math.abs(s.c['e-'][mid] / (ND + p) - 1) < 1e-6);
    }
  }
});

test('SRH at a face: zero in the dark, the law at the two edge nodes under light, and every pair either collected or recombined there', () => {
  // n⁺ | p silicon whose holes can't cross the face, only recombine there with the n⁺ side's
  // electrons: e⁻(left) + h⁺(right) = 0 with SRH velocities. Light makes pairs in the p side.
  const vn = 1e3, vp = 1e3, G = 10, Lp = 5e-6; // m/s, m/s, mol/(m³·s), m
  // (The p side is its own material, lit, so that no holes are made where they could neither
  // leave nor recombine.)
  const cell = (V, light) =>
    build({
      T,
      species: Si.species,
      materials: { n: Si.materials.Si, p: Si.materials.Si },
      stack: [
        ohmic(0, ['e-']),
        layer('n', 1e-6, { donors: units.perCm3(1e17), c0: { 'h+': 1e-12 } }),
        { dipole: 0, species: { 'e-': 'equilibrium', 'h+': 'blocked' }, reactions: [{ equation: 'e-(left) + h+(right) = 0', srh: { vn, vp } }] },
        layer('p', Lp, { acceptors: units.perCm3(1e16) }),
        ohmic(V, ['h+']),
      ],
      bulkReactions: light ? [{ equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { p: G } }] : [],
    });
  const dark = new Device(cell(0, false)).solve();
  assert.ok(dark.converged);
  assert.ok(Math.abs(dark.interfaces[0].rates[0]) < 1e-20, `dark rate ${dark.interfaces[0].rates[0]}`);
  for (const V of [0, 0.4, 0.5]) {
    const s = new Device(cell(V, true)).solve();
    assert.ok(s.converged);
    const gL = s.x.findIndex((x, k) => s.x[k + 1] === x); // the face's left node (it's doubled)
    const n = s.c['e-'][gL], p = s.c['h+'][gL + 1];
    const a = (s.mu['e-'][gL] + s.mu['h+'][gL + 1]) / RT;
    const n1 = Math.sqrt(n * p * Math.exp(-a)); // midgap, against the state
    const law = (n * p * -Math.expm1(-a)) / ((n + n1) / vp + (p + n1) / vn);
    const r = s.interfaces[0].rates[0];
    assert.ok(Math.abs(r / law - 1) < 1e-9, `V ${V}: rate ${r} vs the law ${law}`);
    // Every pair made in the p side is collected (electrons out by the left) or recombines at
    // the face; this is the current. (At 0.4 V the face takes 70% of them; at 0.5 V, what the
    // bias injects outweighs the light.)
    assert.ok(Math.abs(s.current / (FARADAY * (G * Lp - r)) - 1) < 1e-9, `V ${V}: ${s.current} vs ${FARADAY * (G * Lp - r)}`);
  }
});

test('SRH definitions are checked', () => {
  const def = (rx) => build({ T, library: [Si], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1 }), ohmic(0)], bulkReactions: [rx] });
  for (const [rx, message] of [
    [{ equation: 'e- + h+ = 0', kf: { Si: 1 }, srh: { Si: { tauN: 1, tauP: 1 } } }, /kf \(mass action\) or srh/],
    [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1 } } }, /tauP must be/],
    [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1, tauP: 1, vn: 1 } } }, /vn/],
    [{ equation: 'photon = e- + h+', fixed: { photon: 1e5 }, srh: { Si: { tauN: 1, tauP: 1 } } }, /consuming one negative and one positive species/],
  ]) {
    assert.throws(() => new Device(def(rx)), (e) => e instanceof DeviceError && message.test(e.message), JSON.stringify(rx));
  }
});

test('explicit immobile traps in steady state recombine exactly at the SRH rate, and keep their total node by node', () => {
  // The same lit n-type slab, recombining through trap states X⁰/X⁻ that don't move (D = 0):
  // e⁻ + X⁰ = X⁻ (capture coefficient c_n) and X⁻ + h⁺ = X⁰ (c_p). In steady state a single
  // trap level recombines at exactly R_SRH(n, p) with τ_n = 1/(c_n N_t), τ_p = 1/(c_p N_t), and
  // n₁ the electron density that half-fills it. The trap total X⁰ + X⁻ is conserved at each node
  // on its own, which the steady solve now holds as a row of that node's block (no time steps).
  const ND = units.perCm3(1e15), Nt = units.perCm3(1e12), tauN = 2e-6, tauP = 1e-6, G = units.perCm3(1e20);
  const n1 = units.perCm3(1e12);
  const cn = 1 / (tauN * Nt), cp = 1 / (tauP * Nt);
  const species = [...Si.species, { name: 'X0', z: 0, cRef: Nt }, { name: 'X-', z: -1, cRef: Nt }];
  const mat = { ...Si.materials.Si, species: { ...Si.materials.Si.species, X0: { D: 0, mu0: 0 }, 'X-': { D: 0, mu0: el.mu0 + RT * Math.log(n1 / el.cRef) } } };
  const dev = new Device(
    build({
      T,
      species,
      materials: { Si: mat },
      stack: [ohmic(0, ['e-']), layer('Si', 10e-6, { donors: ND, c0: { 'h+': 1e-12, X0: Nt / 2, 'X-': Nt / 2 } }), ohmic(0, ['e-'])],
      bulkReactions: [
        { equation: 'e- + X0 = X-', kf: { Si: cn } },
        { equation: 'X- + h+ = X0', kf: { Si: cp } },
        { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: G } },
      ],
    }),
  );
  const s = dev.solve();
  assert.ok(s.converged);
  assert.ok(s.history.every((h) => !(h.dt < Infinity)), 'solved directly, with no time steps');
  const mid = s.x.length >> 1, n = s.c['e-'][mid], p = s.c['h+'][mid];
  const R = (n * p - ni2) / (tauP * (n + n1) + tauN * (p + ni2 / n1));
  assert.ok(Math.abs(R / G - 1) < 1e-6, `R_SRH ${R} against G ${G}`);
  s.c.X0.forEach((x, g) => assert.ok(Math.abs((x + s.c['X-'][g]) / Nt - 1) < 1e-12, `trap total at node ${g}`));
});
