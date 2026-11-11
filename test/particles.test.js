import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, bath, aqueous, particles, vacuumZeroCharge } from '../src/kit.js';

// The dots are a stochastic sample of a solution: in a long run each cell holds, on average, the
// dots its concentration gives, and the dots crossing each boundary, net, are its flux. Seeded,
// so each run is the same; the tolerances are four standard deviations of the sampling.

const RT = GAS_CONSTANT * 298.15;
const seeded = (seed) => () => ((seed = (seed * 16807) % 2147483647) - 1) / 2147483646;

// A neutral species diffusing through a slab between two baths: c falls linearly from 1 to 0.1
// mol/m³, J = D Δc/L. With `decay`, it also decays in the slab, X → nothing (its equilibrium
// concentration ~e⁻²⁰), and the far side is closed.
const slab = ({ decay = 0 } = {}) => {
  const mu0 = decay ? 20 * RT : 0;
  return new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { gel: { epsr: 0, species: { X: { D: 1e-9, mu0, cRef: 1 } } } },
    regions: [{ material: 'gel', length: 1e-4 }],
    bulkReactions: decay ? [{ equation: 'X = 0', kf: { gel: decay } }] : [],
    contacts: {
      left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: mu0 } } },
      right: decay ? { phi: 'neutral' } : { phi: 'neutral', species: { X: { type: 'equilibrium', mu: mu0 + RT * Math.log(0.1) } } },
    },
    grid: { hmax: 1e-6 },
  }).solve();
};

// Run a swarm, sampling the dots in each cell every step after a warm-up.
const run = (sw, nm, dt, steps) => {
  for (let i = 0; i < steps / 10; i++) sw.step(dt);
  sw.crossed[nm].up.fill(0);
  sw.crossed[nm].down.fill(0);
  const hist = new Float64Array(sw.cells.length);
  for (let i = 0; i < steps; i++) {
    sw.step(dt);
    for (const d of sw.dots) if (d.species === nm) hist[d.cell]++;
  }
  return { mean: hist.map((v) => v / steps), T: dt * steps };
};

test('particles: through a slab, the dots fill each cell as its concentration does, and cross it net at the flux', () => {
  const sol = slab();
  const J = (1e-9 * 0.9) / 1e-4;
  assert.ok(Math.abs(sol.flux.X[10] / J - 1) < 1e-9);
  const sw = particles(sol, { dots: 500, cells: 20, random: seeded(1) });
  const h = sw.cells[0].x1 - sw.cells[0].x0;
  const { mean, T } = run(sw, 'X', (0.2 * h * h) / 1e-9, 20000);
  // The counts are correlated in time (a dot stays ~h²/2D), so allow 3% where Poisson would say 1%.
  const ex = sw.expected('X');
  mean.forEach((m, K) => assert.ok(Math.abs(m / ex[K] - 1) < 0.03, `cell ${K}: ${m} dots on average, expected ${ex[K]}`));
  // Net crossings at the middle boundary: J T / w, give or take √(gross).
  const k = sw.cells.length >> 1, { up, down } = sw.crossed.X;
  const net = up[k] - down[k], want = (J * T) / sw.weight.X;
  assert.ok(Math.abs(net - want) < 4 * Math.sqrt(up[k] + down[k]), `${net} dots across, net; expected ${want}`);
});

test('particles: where a species decays, the dots are unmade at its rate, and they fill the profile it leaves', () => {
  const k = 1e-9 / (20e-6) ** 2; // a decay length of 20 µm in the 100 µm slab
  const sol = slab({ decay: k });
  const sw = particles(sol, { dots: 400, cells: 25, random: seeded(3) });
  const h = sw.cells[0].x1 - sw.cells[0].x0;
  const { mean, T } = run(sw, 'X', (0.2 * h * h) / 1e-9, 20000);
  const ex = sw.expected('X');
  mean.forEach((m, K) => {
    if (ex[K] > 5) assert.ok(Math.abs(m / ex[K] - 1) < 0.05, `cell ${K}: ${m} dots on average, expected ${ex[K]}`);
  });
  // All that enters is unmade: net entry from the bath is the decay, D c₀ tanh(L/λ)/λ.
  const lam = 20e-6, J = (1e-9 * Math.tanh(1e-4 / lam)) / lam;
  assert.ok(Math.abs(sol.contacts.left.flux.X / J - 1) < 0.01);
  const { up, down } = sw.crossed.X, net = up[0] - down[0], want = (J * T) / sw.weight.X;
  assert.ok(Math.abs(net - want) < 4 * Math.sqrt(up[0] + down[0]), `${net} dots in from the bath, net; expected ${want}`);
});

test('particles: across a membrane, the dots crossing each way are its one-way fluxes, in Ussing’s ratio', () => {
  const OUT = { 'K+': 5, 'Cl-': 5 }, IN = { 'K+': 140, 'Cl-': 140 }, VT = RT / FARADAY;
  for (const V of [-0.05, 0.03]) {
    const sol = new Device(
      build({
        library: [aqueous(['K+', 'Cl-'], { epsr: 0 })],
        stack: [
          bath(OUT, 0),
          layer('water', 2e-6, { name: 'out', c0: OUT }),
          { phi: { type: 'capacitive', C: 0.01 }, species: { 'K+': { type: 'permeability', P: 1e-5 }, 'Cl-': 'blocked' } },
          layer('water', 2e-6, { name: 'in', c0: IN }),
          bath(IN, V),
        ],
        grid: { hmin: 20e-9, hmax: 100e-9 },
      }),
    ).solve();
    // The face's one-way fluxes, in Ussing's ratio J_in/J_out = c_out/(c_in e^{F V/RT}) at its edges.
    const [inward, outward] = sol.interfaces[0].oneWay['K+'];
    const gL = sol.region.lastIndexOf(0), c = sol.c['K+'], Vm = sol.phi[gL + 1] - sol.phi[gL];
    assert.ok(Math.abs(inward / outward / (c[gL] / (c[gL + 1] * Math.exp(Vm / VT))) - 1) < 1e-9);
    assert.ok(Math.abs((inward - outward) / sol.flux['K+'][gL] - 1) < 1e-9);
    // The dots crossing each way: the one-way flux × T / w, give or take its square root.
    const sw = particles(sol, { species: ['K+'], dots: 300, cells: 10, random: seeded(7) });
    const h = sw.cells[0].x1 - sw.cells[0].x0;
    const { T } = run(sw, 'K+', (0.3 * h * h) / 1.96e-9, 100000);
    const k = sw.cells.findIndex((cl) => cl.region === 1), w = sw.weight['K+'];
    for (const [count, flux] of [[sw.crossed['K+'].up[k], inward], [sw.crossed['K+'].down[k], outward]]) {
      const want = (flux * T) / w;
      assert.ok(Math.abs(count - want) < 4 * Math.sqrt(want), `V ${V}: ${count} dots across, expected ${want}`);
    }
  }
});

test('particles: in equilibrium, a reaction still makes and unmakes dots, each at its one-way rate', () => {
  // X = nothing in a closed gel, at equilibrium (c = c_ref = 1 mol/m³): forward and backward
  // rates both k·c, which the solution reports, and which the dots show as made and unmade.
  const k = 50;
  const sol = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { gel: { epsr: 0, species: { X: { D: 1e-9, mu0: 0, cRef: 1 } } } },
    regions: [{ material: 'gel', length: 1e-4, c0: { X: 1 } }],
    bulkReactions: [{ equation: 'X = 0', kf: { gel: k } }],
    contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
  }).solve();
  const rx = sol.bulkReactions[0];
  assert.deepEqual(rx.nu, { X: -1 });
  assert.ok(Math.abs(rx.forward[5] - k) < 1e-9 && Math.abs(rx.rate[5]) < 1e-9);
  const sw = particles(sol, { dots: 300, cells: 10, random: seeded(11) });
  const dt = 1e-3, steps = 4000;
  let made = 0, unmade = 0;
  for (let i = 0; i < steps; i++) {
    sw.step(dt);
    for (const ev of sw.events) {
      if (ev.kind === 'made') made++;
      if (ev.kind === 'unmade') unmade++;
    }
  }
  // each about k × (dots) × T
  const want = k * (1e-4 / sw.weight.X) * dt * steps;
  for (const count of [made, unmade]) assert.ok(Math.abs(count - want) < 4 * Math.sqrt(want), `${count} made or unmade, expected ${want}`);
});

test('particles: a step far longer than the hop time is still a sample, density and net crossings, and costs little', () => {
  const sol = slab();
  const J = (1e-9 * 0.9) / 1e-4;
  const sw = particles(sol, { dots: 500, cells: 20, random: seeded(5) });
  const h = sw.cells[0].x1 - sw.cells[0].x0;
  const t0 = performance.now();
  const { mean, T } = run(sw, 'X', (500 * h * h) / 1e-9, 2000); // ~1000 hops a step
  assert.ok((performance.now() - t0) / 2200 < 2, 'a step costs under 2 ms');
  const ex = sw.expected('X');
  mean.forEach((m, K) => assert.ok(Math.abs(m / ex[K] - 1) < 0.05, `cell ${K}: ${m} dots on average, expected ${ex[K]}`));
  const k = sw.cells.length >> 1, { up, down } = sw.crossed.X;
  const net = up[k] - down[k], want = (J * T) / sw.weight.X;
  assert.ok(Math.abs(net - want) < 4 * Math.sqrt(up[k] + down[k]), `${net} dots across, net; expected ${want}`);
});

test('particles: a solution on another grid redraws the lattice; nonsense is an error', () => {
  const a = slab(), sw = particles(a, { dots: 300, random: seeded(9) });
  const b = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { gel: { epsr: 0, species: { X: { D: 1e-9, mu0: 0, cRef: 1 } } } },
    regions: [{ material: 'gel', length: 2e-4 }],
    contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral' } },
    grid: { hmax: 1e-6 },
  }).solve();
  sw.update(b);
  assert.ok(Math.abs(sw.cells.at(-1).x1 - 2e-4) < 1e-12, 'the lattice spans the new device');
  const n = sw.dots.length, want = sw.expected('X').reduce((s, v) => s + v, 0);
  assert.ok(Math.abs(n - want) < 4 * Math.sqrt(want), `${n} dots, expected ${want}`);
  const [t, dots] = [sw.time, sw.dots.length];
  sw.step(NaN).step(-1); // (a first frame; a clock stepped back) nothing moves
  assert.ok(sw.time === t && sw.dots.length === dots);
  assert.throws(() => sw.step(Infinity), /step\(dt\)/);
  assert.throws(() => sw.update(null), /update takes a solution/);
  assert.throws(() => particles(a, { dots: 0 }), /dots must be/);
  assert.throws(() => particles(a, { weight: -1 }), /weight of X/);
  assert.throws(() => particles(a, { weight: 1e-15 }), /would draw/);
  assert.throws(() => sw.expected('Y'), /no species "Y"/);
});

test('particles: about as many cells as asked, and seas by density, so a wider cell among the rest is no sea', () => {
  // Nodes 1 µm apart: cells 1.25 µm wide are drawn one node each, not two (the nearer width).
  const a = slab();
  for (const cells of [80, 120]) {
    const n = particles(a, { cells, dots: 300 }).cells.length;
    assert.ok(n > cells / 1.3 && n < cells * 1.3, `${n} cells for ${cells} asked`);
  }
  // A uniform concentration is one density everywhere: no seas at a cap just above it, though
  // each cell (one node, 1 µm, where 0.75 µm was asked) holds more than cap dots.
  const flat = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { gel: { epsr: 0, species: { X: { D: 1e-9, mu0: 0, cRef: 1 } } } },
    regions: [{ material: 'gel', length: 1e-4 }],
    contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral' } },
    grid: { hmax: 1e-6 },
  }).solve();
  const weight = 1e-8, perTarget = (1 * (1e-4 / 133)) / weight; // dots in a cell of the asked width
  const sw = particles(flat, { cells: 133, weight, cap: 1.2 * perTarget });
  assert.ok(sw.cells.length < 110, 'cells a node wide, wider than asked');
  assert.equal(Array.from(sw.sea.X).reduce((a, b) => a + b), 0, 'no seas');
  assert.ok(Math.abs(sw.dots.length / (1e-4 / weight) - 1) < 0.05, 'all of it in dots');
});

test('particles: a solution off a plane is refused (its fluxes are per device, the amounts per area)', () => {
  const sphere = new Device({
    species: [{ name: 'X', z: 0, cRef: 1000 }],
    materials: { m: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } },
    regions: [{ material: 'm', length: 1e-4, c0: { X: 1 } }],
    contacts: { left: { species: { X: { type: 'equilibrium', mu: 0 } }, phi: 'neutral' }, right: { species: { X: 'blocked' }, phi: 'neutral' } },
    geometry: { type: 'spherical', r0: 1e-6 },
  }).solve();
  assert.equal(sphere.geometry, 'spherical');
  assert.throws(() => particles(sphere), /planar device's solution, not a spherical one/);
});

test('a metal between two semiconductors that share its carrier: no dot crosses it inside the lattice', () => {
  // (left out of the lattice, the metal let the cells either side become neighbours: the right
  // side's dots drained into the left, its first cell 0.1 dots against 33.5 expected; Fable's review)
  const Nc = units.perCm3(2.8e19), ND = units.perCm3(1e16);
  const ohm = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const def = {
    species: [{ name: 'e-', z: -1 }],
    materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc } } }, Au: { conductor: { species: 'e-', conductivity: 4e7 } } },
    regions: [{ material: 'Si', length: 1e-6, fixedCharge: ND * FARADAY }, { material: 'Au', length: 50e-9 }, { material: 'Si', length: 1e-6, fixedCharge: ND * FARADAY }],
    contacts: { left: ohm(0), right: ohm(0.05) },
    grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
  };
  const zc = vacuumZeroCharge(def, 4.2, { material: 'Si', anchor: 'e-', offset: 4.05 });
  def.interfaces = [{ phi: { type: 'capacitive', C: 10 }, zeroCharge: zc }, { phi: { type: 'capacitive', C: 10 }, zeroCharge: zc }];
  const sol = new Device(def).solve();
  let seed = 7;
  const random = () => ((seed = (seed * 16807) % 2147483647) - 1) / 2147483646;
  const sw = particles(sol, { dots: 400, cells: 30, random });
  const metal = sw.cells.findIndex((c) => c.region === 1), first = metal + 1;
  assert.ok(metal > 0 && sw.cells[first].region === 2);
  const h = sw.cells[0].x1 - sw.cells[0].x0, dt = (0.05 * h * h) / 36e-4, steps = 4000;
  let through = 0, mean = 0;
  for (let k = 0; k < steps; k++) {
    const before = new Map(sw.dots.map((d) => [d, d.cell]));
    sw.step(dt);
    for (const d of sw.dots) {
      const b = before.get(d);
      if (b !== undefined && (b < metal) !== (d.cell < metal)) through++;
      if (d.cell === first) mean += 1 / steps;
    }
  }
  assert.equal(through, 0);
  const want = sw.expected('e-')[first];
  assert.ok(Math.abs(mean / want - 1) < 0.25, `the first cell past the metal: ${mean} dots against ${want}`);
});
