import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units, FARADAY } from '../src/index.js';
import { build, layer, ohmic, semiconductor, metal, bath, aqueous } from '../src/kit.js';
import { newton } from '../src/newton.js';
import { choice, device } from '../tools/lawsweep-devices.js';

// Newton's iteration (src/newton.js) against its contract, on what a cold review of it found.
// The lost mode a time step's plain solves can't see is in test/ports.test.js (MOS without a
// port: the gate's charge after the step, to 1e-6 of the steady one).

// Each steady newton() call on a device, with whether the islands' rows were in when it returned.
const steadyCalls = (dev) => {
  const s = dev.solver, calls = [], orig = s.newton.bind(s);
  s.newton = (dt, opts) => {
    const r = orig(dt, opts);
    if (dt === Infinity) calls.push({ ...r, asked: !!opts?.islands, on: s.islandsOn });
    return r;
  };
  return calls;
};

test('newton: at equilibrium the terminal currents are round-off, and the islands round stays off', () => {
  // Every current 1e-20 to 1e-31 A/m², their sum as large as they are: off Kirchhoff by O(1)
  // of the largest, which put the islands' rows in, twice, for nothing (2 of 6 iterations).
  for (const i of [288, 365]) {
    const dev = new Device(device(choice(i), { charged: true }));
    const calls = steadyCalls(dev), sol = dev.solve();
    assert.ok(sol.converged && dev.solver.islands.length > 0, `device ${i}`);
    assert.ok(calls.every((c) => !c.on), `device ${i}: ${JSON.stringify(calls.map((c) => [c.iterations, c.on]))}`);
  }
});

test('newton: a call with the islands\' rows in from the start keeps them', () => {
  // step()'s second try (islands: true) converges at its 59th iteration; the islands round ran
  // there as if the rows had just gone in, into maxIter, and came back with s.islandsOn false on
  // a state solved with them in.
  const dev = new Device(device(choice(200), { charged: true }));
  const calls = steadyCalls(dev), sol = dev.solve();
  assert.ok(sol.converged);
  const asked = calls.filter((c) => c.asked);
  assert.ok(asked.length > 0 && asked.every((c) => c.converged && c.on), JSON.stringify(asked.map((c) => [c.iterations, c.on])));
});

test('newton: a steady solve whose residual isn\'t solved says so, and records how conditioned it was', () => {
  // The residual check reassembles, which cleared the factorisation its conditioning note reads:
  // the note came out empty, under a message about a perturbed pivot there wasn't.
  const dev = new Device(device(choice(5), { charged: true })), s = dev.solver;
  const orig = s.newton.bind(s);
  let r = null;
  s.newton = (dt, opts) => {
    if (r || dt !== Infinity) return orig(dt, opts);
    s.conditioning = null;
    return (r = newton(s, dt, { ...opts, tol: 1e9 })); // (an update under 1e9 converged: the first)
  };
  dev.solve();
  assert.ok(r && !r.converged && s.sys.perturbedPivots === 0, JSON.stringify(r));
  assert.match(r.error, /residual is .* of its largest term/);
  assert.ok(s.conditioning?.digits > 0, JSON.stringify(s.conditioning));
});

test('newton: a floating base, held only by its holes\' leaks through the n layers, settles where they balance', () => {
  // An n/p⁺/n GaAs stack without recombination: the base's holes reach the contacts only through
  // the n layers, where they're n_i²/N ~ 1e-5 per cm³ and conduct ~1e-24 of what they do at the
  // junctions, so the base's level is held by leaks of ~3e-14 A/m². Summed at the junction faces,
  // the island's balance read them as differences of terms that much larger, its round was
  // dropped, and a cold solve at −0.2 V converged with the base at the emitter's level, its two
  // currents 3.3e-11 and −3.9e-11 A/m² (the stress run's semi 220, with a Schottky emitter).
  // Summed out to each n layer's bottleneck, it finds the level; a base 4 V_T above it walks down
  // a thermal unit an update (at −0.3 V, cut off at 12 iterations, 1.6e-9 against −1.9e-9).
  // The expected level, from the solution's φ alone: in each n layer the hole flux is uniform,
  // J = D (e^{v_a} − e^{v_b}) / ∫ e^{(Fφ + μ°)/RT}/N_V dx (Slotboom; a sum of positive terms),
  // and the base sits where what comes in from one side leaves by the other.
  const mat = 'GaAs', def = build({
    T: 300,
    library: [semiconductor(mat)],
    stack: [ohmic(0), layer(mat, 400e-9, { donors: units.perCm3(1e17) }), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), layer(mat, 900e-9, { donors: units.perCm3(3e17) }), ohmic(0)],
  });
  const hp = def.materials[mat].species['h+'], RT = 8.314462618 * 300, VT = RT / 96485.33212;
  const balanced = (s) => {
    const { x, phi, region } = s, w = (g) => (phi[g] + hp.mu0 / 96485.33212) / VT;
    const seg = (g) => {
      const a = w(g), b = w(g + 1), dx = x[g + 1] - x[g];
      return (Math.abs(a - b) < 1e-12 ? dx * Math.exp(a) : (dx * (Math.exp(b) - Math.exp(a))) / (b - a)) / hp.cRef;
    };
    const mid = (region.indexOf(1) + region.lastIndexOf(1)) >> 1, Vh = s.V['h+'];
    let IL = 0, IR = 0;
    for (let g = 0; g < x.length - 1; g++) if (x[g + 1] > x[g]) g < mid ? (IL += seg(g)) : (IR += seg(g));
    const vl = Vh[0] / VT, vr = Vh.at(-1) / VT;
    return { base: Vh[mid], level: VT * Math.log((Math.exp(vl) / IL + Math.exp(vr) / IR) / (1 / IL + 1 / IR)) };
  };
  for (const V of [-0.5, -0.3, -0.2, 0.3]) {
    const d = new Device(def);
    d.set({ contacts: { right: { V } } });
    const s = d.solve(), { left, right } = s.terminals, { base, level } = balanced(s);
    assert.ok(s.converged && s.warnings.length === 0, `${V} V`);
    assert.ok(Math.abs(left.current + right.current) < 1e-6 * Math.abs(left.current), `${V} V: ${left.current} and ${right.current} A/m²`);
    assert.ok(Math.abs(base - level) < 1e-6, `${V} V: the base at ${base} V, its leaks balanced at ${level} V`);
  }
});

test('a steady solve whose terminal currents don\'t add up says so', () => {
  // A fourth layer beyond the floating base (n/p⁺/n/n⁺): its right neighbour isn't the stretch's
  // end piece, so the base's island is summed at that face, and from cold at −0.2 V the currents
  // were 3.29e-11 and −3.83e-11 A/m² (warm, 1.96e-13 both), converged with nothing said. Not
  // warned: equilibrium's round-off (the law sweep's dropped rounds, 1e-18 to 2e-17 A/m² apart),
  // and the three-layer base, which now solves.
  const mat = 'GaAs', n = (N, L) => layer(mat, L, { donors: units.perCm3(N) });
  const stack = (more) => [ohmic(0), n(1e17, 400e-9), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), n(3e17, 600e-9), ...more, ohmic(0)];
  const at = (more, V) => {
    const d = new Device(build({ T: 300, library: [semiconductor(mat)], stack: stack(more) }));
    d.set({ contacts: { right: { V } } });
    return d.solve();
  };
  const off = at([n(4e17, 300e-9)], -0.2), { left, right } = off.terminals;
  assert.ok(Math.abs(left.current + right.current) > 1e-3 * Math.abs(left.current), 'this device now solves: find another for this test');
  assert.ok(off.warnings.some((w) => w.includes("don't add up")), JSON.stringify(off.warnings));
  assert.deepEqual(at([], -0.3).warnings, []);
  const eq = new Device(device(choice(288), { charged: true })).solve();
  assert.ok(eq.converged && !eq.warnings.some((w) => w.includes("don't add up")));
});

test('…and says so where one terminal reads ~0, which is not round-off unless a terminal floats', () => {
  // The same base behind a contact that passes holes only: the true current is the holes' leak,
  // ±3.594e-15 A/m² at 0.3 V (a transient from the cold state, to 1e15 s, and a solve from
  // there, agree). From cold it converged at −1.05e-20 and 4.25e-15: the whole sum on one
  // terminal, which the warning once took for round-off (read only under half the largest
  // current). Round-off of that shape comes from a floating terminal, whose circuit's current
  // and fluxes differ (the liquid-junction tests); with none, it's warned. At 0.1 V it solves.
  const mat = 'GaAs', d = (V) => {
    const dev = new Device(build({ T: 300, library: [semiconductor(mat)], stack: [ohmic(0, ['h+']), layer(mat, 400e-9, { donors: units.perCm3(1e17) }), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), layer(mat, 600e-9, { donors: units.perCm3(3e17) }), ohmic(0)] }));
    dev.set({ contacts: { right: { V } } });
    return dev.solve();
  };
  const off = d(0.3), { left, right } = off.terminals, sum = Math.abs(left.current + right.current);
  assert.ok(sum > 0.5 * Math.max(Math.abs(left.current), Math.abs(right.current)), `this device now solves (${left.current}, ${right.current}): find another for this test`);
  assert.ok(off.warnings.some((w) => w.includes("don't add up")), JSON.stringify(off.warnings));
  const ok = d(0.1);
  assert.deepEqual(ok.warnings, []);
  assert.ok(Math.abs(ok.terminals.left.current + 3.519e-15) < 1e-3 * 3.519e-15 && Math.abs(ok.terminals.left.current + ok.terminals.right.current) < 1e-6 * 3.519e-15, JSON.stringify(ok.terminals));
});

test('newton: a floating gate at rest in inversion stays put (a mode the plain solves never see, found by a probe)', () => {
  // Its electrons reach the back contact only by minority diffusion, so the factorised solves
  // lose the inversion layer's level; Newton converged quadratically on the rest at every step,
  // and the extrapolated guess's error in the level grew: here 5e-6 V by 1e7 s, through 63,000
  // steps (8.5e-5 V by 1e5 s, reached by way of 1e3 s).
  const def = build({
    T: 300,
    library: [semiconductor('Si'), metal('Au')],
    materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [
      ohmic(0),
      layer('Si', units.um(0.5), { name: 'p-Si', acceptors: units.perCm3(1e17) }),
      { dipole: 0 },
      layer('SiO2', units.nm(5), { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }),
      { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.1 },
      layer('Au', units.nm(50)),
      ohmic(0, ['e-']),
    ],
    grid: { hmin: units.nm(0.1), hmax: units.nm(10), ratio: 1.1 },
  });
  const dev = new Device(def).set({ contacts: { right: { V: 1.2 } } });
  dev.solve();
  dev.set({ contacts: { right: { V: undefined, I: 0 } } });
  for (const t of [1e5, 1e7]) {
    const s = dev.advance(t);
    assert.ok(s.converged && s.steps < 100 && Math.abs(s.contacts.right.V - 1.2) < 1e-9, `t = ${t}: ${s.contacts.right.V - 1.2} V, ${s.steps} steps`);
  }
});

test('…but not where the species\' currents cancel: a membrane at no current between two baths', () => {
  // A cation-exchange membrane between 500 and 10 mM KCl, open circuit: K⁺ and Cl⁻ cross at
  // ±93 A/m² and the held end read 5e-13 A/m² (5e-15 of them), warned of on every such solve
  // once only floating terminals' gaps were allowed for (found by a cold user, 4 Nov). Its EMF is
  // Teorell–Meyer–Sievers' (that user checked it to 1e-8 refined).
  const X = 200, lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), w = lib.materials.water.species;
  lib.materials.cem = { epsr: 0, species: { 'K+': { D: 2e-10, mu0: w['K+'].mu0, cRef: 1000 }, 'Cl-': { D: 1e-10, mu0: w['Cl-'].mu0, cRef: 1000 } } };
  const coion = (c) => (-X + Math.sqrt(X * X + 4 * c * c)) / 2;
  for (const [c1, L] of [[500, 50e-6], [500, 5e-6], [1000, 50e-6]]) {
    const s = new Device(build({ T: 298.15, library: [lib], grid: { hmin: 0.2e-6, hmax: 2e-6 }, stack: [bath({ 'K+': c1, 'Cl-': c1 }, 'Cl-', 0, { material: 'water' }), layer('cem', L, { fixedCharge: -X * FARADAY, c0: { 'K+': X + coion(c1), 'Cl-': coion(c1) } }), bath({ 'K+': 10, 'Cl-': 10 }, 'Cl-', { I: 0 }, { material: 'water' })] })).solve();
    assert.ok(Math.abs(FARADAY * s.flux['K+'][s.x.length >> 1]) > 20, 'the salt crosses');
    assert.deepEqual(s.warnings, [], `${c1} mM, ${L} m`);
  }
});
