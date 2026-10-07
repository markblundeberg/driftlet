import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, SolverError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, bath, semiconductor, SEMICONDUCTORS, recombination, describe } from '../src/kit.js';

// An initial profile as plain data: c0 tabulated against x, the starting state of a transient
// (even for a species a contact feeds), against closed forms that start from that profile.

const sample = (L, f, n = 400) => {
  const x = Array.from({ length: n + 1 }, (_, k) => (k * L) / n);
  return { x, values: x.map(f) };
};

// ∫ (c − base) w(x) dx by the trapezoid rule on the solution's nodes.
const integral = (s, name, base, w = () => 1) => {
  let m = 0;
  for (let g = 1; g < s.x.length; g++) {
    const h = s.x[g] - s.x[g - 1];
    for (const k of [g - 1, g]) m += (s.c[name][k] - base[k]) * w(s.x[k]) * (h / 2);
  }
  return m;
};

test('Haynes–Shockley from a profile: a Gaussian hole packet laid down at t = 0 drifts, spreads and decays at the ambipolar rates', () => {
  // The same n-Ge bar under 10 V/cm as the injector's test, but the packet is there from the
  // start: rest + 0.1% of n₀ in a Gaussian of σ₀ = 30 µm. Nothing else is set, so the device
  // starts from its cold start (electrons at the left contact's level, made neutral node by node
  // around the holes) and the field sets up within picoseconds. Then exactly: the centre moves
  // μ*E t, the variance is σ₀² + 2D*t, and the amount decays at k_f(n₀ + p₀).
  const ND = units.perCm3(1e15), L = 3e-3, E = 1000, tau = 20e-6, T = 300, xe = 0.5e-3, s0 = 30e-6, t = 5e-6;
  const reaction = recombination({ material: 'Ge', tau, majority: ND });
  const def = (c0) =>
    build({
      T,
      library: [semiconductor('Ge')],
      stack: [ohmic(E * L), layer('Ge', L, { name: 'bar', donors: ND, c0 }), ohmic(0)],
      bulkReactions: [reaction],
      grid: { hmin: 5e-6, hmax: 5e-6 },
    });
  const rest = new Device(def(undefined)).solve();
  const mid = rest.x.length >> 1, n0 = rest.c['e-'][mid], p0 = rest.c['h+'][mid];
  const dp = 0.001 * ND;
  const dev = new Device(def({ 'h+': sample(L, (x) => p0 + dp * Math.exp(-((x - xe) ** 2) / (2 * s0 * s0)), 600) }));
  const moments = (s) => {
    const base = rest.c['h+'], m0 = integral(s, 'h+', base), m1 = integral(s, 'h+', base, (x) => x), m2 = integral(s, 'h+', base, (x) => x * x);
    return { amount: m0, mean: m1 / m0, variance: m2 / m0 - (m1 / m0) ** 2 };
  };
  const a = moments(dev.solution());
  assert.ok(Math.abs(a.amount / (dp * s0 * Math.sqrt(2 * Math.PI)) - 1) < 1e-12, 'the profile is the starting state');
  // The electrons follow the holes at the start: neutral node by node.
  const s0s = dev.solution(), peak = s0s.x.findIndex((x) => Math.abs(x - xe) < 1e-12);
  assert.ok(Math.abs((s0s.c['e-'][peak] - s0s.c['h+'][peak]) / ND - 1) < 1e-9);
  const s = dev.advance(t);
  assert.ok(s.converged);
  const b = moments(s);
  const VT = (GAS_CONSTANT * T) / FARADAY, mun = SEMICONDUCTORS.Ge.mun * 1e-4, mup = SEMICONDUCTORS.Ge.mup * 1e-4;
  const muStar = ((n0 - p0) * mun * mup) / (n0 * mun + p0 * mup), DStar = ((n0 + p0) * mun * mup * VT) / (n0 * mun + p0 * mup);
  const drift = (b.mean - xe) / (muStar * E * t), spread = (b.variance - s0 * s0) / (2 * DStar * t);
  const decay = Math.log(a.amount / b.amount) / (reaction.kf.Ge * (n0 + p0) * t);
  assert.ok(Math.abs(drift - 1) < 3e-4, `drift: ${drift} of μ*E`);
  assert.ok(Math.abs(spread - 1) < 1e-2, `spread: ${spread} of 2D*`);
  assert.ok(Math.abs(decay - 1) < 1e-3, `decay: ${decay} of k_f(n₀ + p₀)`);
});

test('a salt profile relaxes as its diffusion mode: c̄ + a·sin(πx/2L) decays at π²D/4L², with D ambipolar', () => {
  // KCl in a strictly neutral electrolyte between a bath (c̄) and a blocking wall. A binary salt
  // diffuses exactly linearly, with D = 2D₊D₋/(D₊ + D₋), so the profile's one mode decays as a
  // pure exponential.
  const L = 100e-6, cbar = 10, amp = 5, Dp = 1.96e-9, Dm = 2.03e-9;
  const profile = sample(L, (x) => cbar + amp * Math.sin((Math.PI * x) / (2 * L)));
  const dev = new Device(
    build({
      species: [
        { name: 'K+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'K+': { D: Dp, mu0: 0 }, 'Cl-': { D: Dm, mu0: 0 } } } },
      stack: [bath({ 'K+': cbar, 'Cl-': cbar }, 'Cl-'), layer('water', L, { c0: { 'K+': profile, 'Cl-': profile } }), { species: { 'K+': 'blocked', 'Cl-': 'blocked' }, phi: 'neutral' }],
      grid: { hmin: 2e-6, hmax: 2e-6 },
    }),
  );
  const cb = new Float64Array(dev.solution().x.length).fill(cbar);
  const mode = (s) => (2 / L) * integral(s, 'K+', cb, (x) => Math.sin((Math.PI * x) / (2 * L)));
  const a0 = mode(dev.solution()), rate = (Math.PI ** 2 * ((2 * Dp * Dm) / (Dp + Dm))) / (4 * L * L), t = 1 / rate;
  assert.ok(Math.abs(a0 / amp - 1) < 1e-4);
  const s = dev.advance(t, { tol: 1e-7 });
  assert.ok(s.converged);
  const ratio = mode(s) / (a0 * Math.exp(-rate * t));
  assert.ok(Math.abs(ratio - 1) < 2e-5, `mode amplitude: ${ratio} of theory`);
  s.c['K+'].forEach((c, k) => assert.ok(Math.abs(c - s.c['Cl-'][k]) < 1e-12 * cbar));
});

test("a spectator's profile fixes the amount it conserves: the steady state holds its mean", () => {
  // Ag | AgNO₃ | Ag at rest, the nitrate (blocked at both electrodes) laid down as a ramp from 5
  // to 15 mol/m³: the steady state is uniform at the ramp's mean, 10.
  const electrode = { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
  const L = 20e-6;
  const dev = new Device({
    species: [
      { name: 'Ag+', z: 1, cRef: 1000 },
      { name: 'NO3-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
    regions: [{ material: 'water', length: L, c0: { 'NO3-': { x: [0, L], values: [5, 15] } } }],
    contacts: { left: electrode, right: electrode },
    grid: { hmin: 0.1e-6, hmax: 0.5e-6 },
  });
  const start = dev.solution();
  start.c['NO3-'].forEach((c, k) => assert.ok(Math.abs(c - (5 + (10 * start.x[k]) / L)) < 1e-9, 'the ramp as given'));
  const s = dev.solve();
  s.c['NO3-'].forEach((c) => assert.ok(Math.abs(c / 10 - 1) < 1e-9, `${c}`));
});

test("a starting profile the grid can't resolve is warned of, with the cells that would", () => {
  // A liquid junction laid down as a tanh 10 µm wide on 40 µm cells starts as something else,
  // and its potential then seems to grow as it forms: silently, until this warning.
  const L = 4e-3, xs = Array.from({ length: 401 }, (_, k) => (k / 400) * L);
  const step = { x: xs, values: xs.map((x) => 10 + 45 * (1 + Math.tanh((x - L / 2) / 10e-6))) };
  const tube = (grid) =>
    new Device({
      species: [
        { name: 'H+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'H+': { D: 9.3e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
      regions: [{ material: 'water', length: L, c0: { 'H+': step, 'Cl-': step } }],
      contacts: { left: { V: 0, terminal: 'Cl-', species: { 'Cl-': 'equilibrium' }, phi: 'bulk' }, right: { I: 0, terminal: 'Cl-', species: { 'Cl-': 'equilibrium' }, phi: 'bulk' } },
      grid,
    }).solution().warnings;
  const coarse = tube({ hmin: 2e-6, hmax: 40e-6 });
  assert.equal(coarse.length, 1, coarse.join('\n'));
  const [, x, h, want] = /^regions\[0\]\.c0 \(H\+, Cl-\): the grid doesn't resolve .*\(at x = (\S+) m, nodes (\S+) m apart .* cells of about (\S+) m there/.exec(coarse[0]);
  assert.ok(Math.abs(x - L / 2) < 20e-6 && Math.abs(h - 40e-6) < 1e-6, coarse[0]);
  assert.deepEqual(tube({ hmin: 2e-6, hmax: +want }), [], 'the cells it suggests resolve it');
});

test('c0 profiles: checked as data, neutral where ε = 0, a number as much a start, and against the region in describe()', () => {
  const wall = { species: { 'K+': 'blocked', 'Cl-': 'blocked' }, phi: 'neutral' };
  const def = (c0) =>
    build({
      species: [
        { name: 'K+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } } },
      stack: [bath({ 'K+': 10, 'Cl-': 10 }, 'Cl-'), layer('water', 1e-4, { c0 }), wall],
    });
  for (const [c0, message] of [
    [{ 'K+': { x: [0, 1e-4], values: [10] } }, /two arrays of the same length/],
    [{ 'K+': { x: [1e-4, 0], values: [10, 10] } }, /strictly increasing/],
    [{ 'K+': { x: [0, 1e-4], values: [10, 0] } }, /values\[1\] must be a concentration > 0/],
    [{ 'K+': { x: [0], values: [10], y: 1 } }, /y/],
    [{ 'K+': [10, 10] }, /a concentration > 0 \(mol\/m³\) or a profile/],
  ]) {
    assert.throws(() => new Device(def(c0)), (e) => e instanceof DeviceError && message.test(e.message), JSON.stringify(c0));
  }
  // Both ions given, so nothing responds: each node must be neutral.
  assert.throws(
    () => new Device(def({ 'K+': { x: [0, 1e-4], values: [10, 20] }, 'Cl-': { x: [0, 1e-4], values: [10, 10] } })).solution(),
    (e) => e instanceof SolverError && /net charge .* at x = /.test(e.message),
  );
  // One given, and the other follows it.
  const s = new Device(def({ 'K+': { x: [0, 1e-4], values: [10, 20] } })).solution();
  s.c['K+'].forEach((c, k) => assert.ok(Math.abs(c - s.c['Cl-'][k]) < 1e-9));
  // A plain number is the starting state too, though the bath feeds K⁺ at 10.
  const n = new Device(def({ 'K+': 15 })).solution();
  n.c['K+'].forEach((c, k) => assert.ok(Math.abs(c / 15 - 1) < 1e-12 && Math.abs(n.c['Cl-'][k] / 15 - 1) < 1e-9));
  // A profile in the region's own x, or in µm, misses the region.
  assert.match(describe(def({ 'K+': { x: [100, 200], values: [10, 20] } })), /misses the region/);
  assert.doesNotMatch(describe(def({ 'K+': { x: [0, 1e-4], values: [10, 20] } })), /misses the region/);
});
