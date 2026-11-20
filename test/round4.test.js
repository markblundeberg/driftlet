import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous, bath, build, layer } from '../src/kit.js';
import { addTerm } from '../src/sum.js';

// Faults found in the discretised equations by round 4's outside reading passes (a second Sol over
// the kernels), each confirmed by a probe against a second route before it was fixed.

const T = 298.15, RT = GAS_CONSTANT * T;

// B6: an affinity summed term by term lost its digits where the participants' μ°/RT are large.
test('addTerm: an affinity summed exactly, in any order (round 4, B6)', () => {
  // Exact sums of doubles as integers in units of 2^-200 (every term here is a multiple of that).
  const big = (x) => BigInt(x * 2 ** 200);
  // large terms that nearly cancel, as A's and B's μ°/RT do, the last leaving a ~ 1e-9
  const terms = [[1, 40.123456789], [-1, 39.98765432101], [3, 0.30000000000000004], [-2, 0.1], [0.5, 7.77]];
  terms.push([-1, terms.reduce((s, [nu, x]) => s + nu * x, 0) - 1e-9]);
  const want = terms.reduce((s, [nu, x]) => s + big(nu * 2) * big(x), 0n); // ν·2 keeps ½ an integer
  const scale = 2 ** 200, wantA = Number(want) / (2 * scale * scale);
  let plainSpread = 0;
  for (let k = 0; k < 720; k++) {
    // every order (6! of them)
    const order = [], rest = [...terms];
    for (let m = k, j = rest.length; j > 0; j--) order.push(rest.splice(m % j, 1)[0]), (m = Math.floor(m / j));
    const acc = new Float64Array(2);
    let plain = 0;
    for (const [nu, x] of order) addTerm(acc, nu, x, 0), (plain += nu * x);
    const a = acc[0] + acc[1];
    assert.ok(Math.abs(a - wantA) <= Math.abs(wantA) * Number.EPSILON + 1e-28, `order ${k}: ${a} vs ${wantA}`);
    plainSpread = Math.max(plainSpread, Math.abs(plain - wantA) / Math.abs(wantA));
  }
  assert.ok(plainSpread > 1e-8, `a plain sum loses a's digits here (${plainSpread})`); // the test's teeth
  // an ∞ passes through (a participant at c = 0)
  const acc = new Float64Array([2, 0]);
  addTerm(acc, -1, Infinity, 0);
  addTerm(acc, 3, 0.1, 0);
  assert.equal(acc[0] + acc[1], -Infinity);
});

// A + B + E = C + D in a 10 µm film, μ°/RT of 40, −40, 0.3, 0.3, 0: A, B, E, C held at the left at
// their μ°, D drawn off at the right half a thermal unit below. The same physics in any order of the
// equation; fast, it's local equilibrium (K = 1), solved by its reaction invariants: c_A = c_B = c_E
// = a, c_C = 2 − a, c_D = s − a, s linear from 2, (2 − a)(s − a) = a³, c_D = e^{−1/2} at the right.
test('a fast bulk reaction among large μ° comes out the same in any order of its equation (round 4, B6)', () => {
  const mu0 = { A: 40, B: -40, E: 0.3, C: 0.3, D: 0 };
  const film = (equation, kf) => new Device({
    T, species: Object.keys(mu0).map((name) => ({ name, z: 0 })),
    materials: { f: { epsr: 4, species: Object.fromEntries(Object.entries(mu0).map(([k, m]) => [k, { D: 1e-9, cRef: 1, mu0: m * RT }])) } },
    regions: [{ name: 'a', material: 'f', length: 1e-5, c0: { A: 1, B: 1, E: 1, C: 1, D: 1 } }],
    bulkReactions: [{ equation, kf: { f: kf } }],
    contacts: {
      left: { phi: 'neutral', species: { A: { type: 'equilibrium', mu: mu0.A * RT }, B: { type: 'equilibrium', mu: mu0.B * RT }, E: { type: 'equilibrium', mu: mu0.E * RT }, C: { type: 'equilibrium', mu: mu0.C * RT }, D: 'blocked' } },
      right: { phi: 'neutral', species: { A: 'blocked', B: 'blocked', E: 'blocked', C: 'blocked', D: { type: 'equilibrium', mu: (mu0.D - 0.5) * RT } } },
    },
    grid: { minCells: 32 },
  }).solve();
  const cD = Math.exp(-0.5);
  let lo = 0, hi = 2;
  for (let k = 0; k < 200; k++) {
    const a = (lo + hi) / 2;
    if (a * a * a - (2 - a) * cD > 0) hi = a;
    else lo = a;
  }
  const exact = (1e-9 * (2 - ((lo + hi) / 2 + cD))) / 1e-5;
  const N = ['A + B + E = C + D', 'A + E + B = C + D', 'E + A + B = C + D'].map((eq) => {
    const s = film(eq, 1e14);
    assert.ok(s.converged, eq);
    return s.contacts.right.flux.D;
  });
  // (summed term by term, these spread by 4.9e-4)
  assert.ok(Math.max(...N) / Math.min(...N) - 1 < 1e-12, `the orders' fluxes ${N}`);
  for (const n of N) assert.ok(Math.abs(n / exact - 1) < 1e-9, `${n} vs local equilibrium's ${exact}`);
});

// The same at a face: A + B + E on the left = C + D on the right, blocked across it otherwise.
test('a face reaction among large μ° comes out the same in any order of its equation (round 4, B6)', () => {
  const mu0 = { A: 40, B: -40, E: 0.3, C: 0.3, D: 0 };
  const cell = (equation) => new Device({
    T, species: Object.keys(mu0).map((name) => ({ name, z: 0 })),
    materials: { f: { epsr: 4, species: Object.fromEntries(Object.entries(mu0).map(([k, m]) => [k, { D: 1e-9, cRef: 1, mu0: m * RT }])) } },
    regions: ['a', 'b'].map((name) => ({ name, material: 'f', length: 1e-5, c0: { A: 1, B: 1, E: 1, C: 1, D: 1 } })),
    interfaces: [{ phi: 'neutral', species: { A: 'blocked', B: 'blocked', E: 'blocked', C: 'blocked', D: 'blocked' }, reactions: [{ equation, k0: 1e14 }] }],
    contacts: {
      left: { phi: 'neutral', species: { A: { type: 'equilibrium', mu: mu0.A * RT }, B: { type: 'equilibrium', mu: mu0.B * RT }, E: { type: 'equilibrium', mu: mu0.E * RT }, C: 'blocked', D: 'blocked' } },
      right: { phi: 'neutral', species: { A: 'blocked', B: 'blocked', E: 'blocked', C: { type: 'equilibrium', mu: mu0.C * RT }, D: { type: 'equilibrium', mu: (mu0.D - 0.5) * RT } } },
    },
    grid: { minCells: 32 },
  }).solve();
  const N = ['A(left) + B(left) + E(left) = C(right) + D(right)', 'A(left) + E(left) + B(left) = C(right) + D(right)', 'E(left) + A(left) + B(left) = C(right) + D(right)'].map((eq) => {
    const s = cell(eq);
    assert.ok(s.converged, eq);
    return s.contacts.right.flux.D;
  });
  // (summed term by term, these spread by 6.7e-15)
  assert.ok(Math.max(...N) / Math.min(...N) - 1 < 1e-15, `the orders' fluxes ${N}`);
});

// B13: the dominant adsorbate's storage diagonal, (Γθ_q/dt)(1 − θ_q), was 0 or noise where the bare
// fraction is below an ulp of 1; test/coverage.test.js's cell, its film 100 kJ/mol more stable.
test('a surface filling to a bare fraction of 4e-19 completes its transient, onto the isotherm (round 4, B13)', () => {
  const G = 1e-5, L = 10e-6, area = 1e4, mu0 = -100e3, V = 0.05;
  const lib = aqueous(['Na+', 'OH-', 'Fe2+'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const dev = new Device(build({
    T, library: [lib],
    stack: [{}, layer('water', L, { name: 'film' }), bath({ 'Na+': 100, 'OH-': 102, 'Fe2+': 1 }, 'OH-')],
    ports: [{ name: 'metal', region: 'film', V, terminal: 'e-', area, surface: { OHads: { mu0, capacity: G, theta0: 1e-4 } }, reactions: [{ equation: 'OH- = OHads + e-', k0: 1e-4, alpha: 0.5, bare: true }] }],
    grid: { hmin: L / 20, hmax: L / 20 },
  }));
  const s = dev.advance(1e-4, { tol: 1e-6 }); // (it stalled at 7.4e-9 s, its steps cut to 3e-19 s)
  assert.ok(s.converged, s.warnings.join('; '));
  const p = s.ports[0], g = s.x.indexOf(p.x[5]);
  const want = 1 / (1 + Math.exp((s.mu['OH-'][g] + FARADAY * V - mu0) / RT));
  assert.ok(want < 1e-18 && Math.abs(p.bare[5] / want - 1) < 1e-6, `bare ${p.bare[5]} vs the isotherm's ${want}`);
});

// B5: a velocity contact's flux v c (e^x − 1) took e^x − 1 as a difference.
test("a velocity contact's flux is exact for a small drive (round 4, B5)", () => {
  // A neutral tracer film, the velocity link on the left (outside μ = 0), equilibrium on the right at
  // μ = −δ RT: linear diffusion, so N = (1 − e^{−δ})/(1/v + L/D) exactly.
  const L = 1e-6, D = 1e-9, v = D / L, delta = 1e-8;
  const s = new Device({
    T, species: [{ name: 'X', z: 0 }], materials: { f: { epsr: 4, species: { X: { D, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: L, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: { type: 'velocity', v, mu: 0 } } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: -delta * RT } } } },
    grid: { minCells: 8 },
  }).solve();
  const exact = -Math.expm1(-delta) / (1 / v + L / D);
  assert.ok(s.converged);
  assert.ok(Math.abs(s.contacts.left.flux.X / exact - 1) < 1e-12, `N ${s.contacts.left.flux.X} vs ${exact}`); // (4.3e-9 before)
});
