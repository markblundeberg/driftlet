import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units, FARADAY, GAS_CONSTANT } from '../src/index.js';
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
  // A floating base whose islands round fails stands at the plain solve's level, its currents
  // apart: said. (Until 9 Nov the example was a fourth layer beyond the base, n/p⁺/n/n⁺, its
  // right neighbour not the stretch's end piece: from cold at −0.2 V, 3.29e-11 and −3.83e-11
  // A/m², where warm gave 1.96e-13 both. That one solves now, its sum reaching through the third
  // layer. This one is the three-layer base at 250 K with its species listed e⁻ last: the plain
  // factorisation loses a mode mid-round, 5.1e-14 and −6.3e-14 A/m², see ROADMAP.) Not warned:
  // equilibrium's round-off (the law sweep's dropped rounds, 1e-18 to 2e-17 A/m² apart), and the
  // bases that solve.
  const mat = 'GaAs', n = (N, L) => layer(mat, L, { donors: units.perCm3(N) });
  const stack = (more) => [ohmic(0), n(1e17, 400e-9), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), n(3e17, more.length ? 600e-9 : 900e-9), ...more, ohmic(0)];
  const at = (more, V, T = 300, reversed = false) => {
    const def = build({ T, library: [semiconductor(mat, { T })], stack: stack(more) });
    if (reversed) def.species = [...def.species].reverse();
    const d = new Device(def);
    d.set({ contacts: { right: { V } } });
    return d.solve();
  };
  const off = at([], -0.3, 250, true), { left, right } = off.terminals;
  assert.ok(Math.abs(left.current + right.current) > 1e-3 * Math.abs(left.current), 'this device now solves: find another for this test');
  assert.ok(off.warnings.some((w) => w.includes("don't add up")), JSON.stringify(off.warnings));
  assert.deepEqual(at([], -0.3).warnings, []);
  assert.deepEqual(at([n(4e17, 300e-9)], -0.2).warnings, []);
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

test('islands: a floating base whose carriers recombine is summed with its rates, in either listing', () => {
  // Fable 5.1's hunt A (4 Nov): the floating p⁺ base at 250 K with band-to-band recombination in
  // all three layers. A piece with a bulk reaction got no island, nor a cut through an end piece
  // with one, and a cold solve at −0.2 V converged with the base at −4e-10 V where it's −0.1477,
  // its currents 4.6e-16 and −4.8e-16 A/m² (true: ±5.03e-19), silent: Kirchhoff's miss, 1.7e-17
  // A/m², was under the warning's floor. Now each node's rate is summed (+vol·ν·r, as in its row),
  // the piece's own and the end pieces' out to their bottlenecks. Listed in reverse, the round
  // diverged where each iterate re-chose the bottleneck (the summed nodes moving under Newton), so
  // a round holds them. Independent routes: a warm sweep from 0 V in 20 steps, −0.147711380 V;
  // advance() at −0.2 V from the 0 V state, −0.14770784 V at 3e18 s (settling ~1e16 s).
  const mat = 'GaAs', def = build({
    T: 250,
    library: [semiconductor(mat, { T: 250 })],
    stack: [ohmic(0), layer(mat, 0.483e-6, { name: 'n1', donors: 0.66 }), layer(mat, 0.22e-6, { name: 'p', acceptors: 8.94 }), layer(mat, 0.902e-6, { name: 'n2', donors: 3.5 }), ohmic(-0.2)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { [mat]: 1e8 } }],
    grid: { hmin: 0.5e-9, hmax: 40e-9, ratio: 1.11 },
  });
  const bases = [];
  for (const species of [def.species, [...def.species].reverse()]) {
    const s = new Device({ ...def, species }).solve(), { left, right } = s.terminals;
    const mid = (s.region.indexOf(1) + s.region.lastIndexOf(1)) >> 1;
    bases.push(s.V['h+'][mid]);
    assert.ok(s.converged && !s.warnings.some((w) => w.includes('floating')), JSON.stringify(s.warnings));
    assert.ok(Math.abs(left.current + right.current) < 1e-6 * Math.abs(left.current), `${left.current} and ${right.current} A/m²`);
    assert.ok(Math.abs(left.current - 5.0273e-19) < 1e-3 * 5.0273e-19, `${left.current} A/m²`);
  }
  assert.ok(Math.abs(bases[0] + 0.1477114) < 1e-6 && Math.abs(bases[0] - bases[1]) < 1e-9, `the base at ${bases} V`);
});

test('islands: the round is asked for by each island\'s own balance, not by what else the device carries', () => {
  // Found by a Fable 5.1 hunter (4 Nov). The round once started only where the terminal currents
  // missed by 1e-9 of the largest and 1e-18 A/m²: a second electrode carrying 1 mA/m² beside the
  // law sweep's F5 electrode (a floating Pt, a fast couple and one 1e-20 times slower) hid its
  // 2.9e-13 miss, and the Pt converged at 0.274922 V, 80 mV off, check() passing it. Its own
  // balance, summed as its pin sums it, missed by 0.5 to 1.5 of what crosses its faces; right
  // states read under 3e-9.
  const def = device({ family: 2, fast: 1, slow: 1e-20, face: 1e-17, drive: 'I0', reversed: false });
  const alone = new Device(def).solve();
  def.ports = [{ name: 'aux', region: 'water2', from: 0.5e-6, to: 1e-6, I: 1e-3, terminal: 'e-', area: 1e3, reactions: [{ equation: 'B- = C2+ + 3 e-', k0: 1, fixed: { 'S(s)': 0 } }] }];
  const beside = new Device(def).solve();
  assert.ok(Math.abs(alone.terminals.left.V - 0.195442) < 1e-5, `alone: ${alone.terminals.left.V}`);
  assert.ok(Math.abs(beside.terminals.left.V - alone.terminals.left.V) < 1e-6 && beside.warnings.length === 0, `beside: ${beside.terminals.left.V} V, ${beside.warnings}`);
  // A floating base at 200 K: its leaks are 5e-27 A/m², under the old floor; it solved at a
  // contact's level, 0.16 V off, cold and swept. Now at its leaks' balance (the hole level's).
  const mat = 'GaAs', base = (V) => {
    const d = new Device(build({ T: 200, library: [semiconductor(mat, { T: 200 })], stack: [ohmic(0), layer(mat, 400e-9, { donors: units.perCm3(1e17) }), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), layer(mat, 900e-9, { donors: units.perCm3(3e17) }), ohmic(0)] }));
    d.set({ contacts: { right: { V } } });
    return d.solve();
  };
  const s = base(-0.3), mid = (s.region.indexOf(1) + s.region.lastIndexOf(1)) >> 1;
  assert.ok(s.converged && Math.abs(s.V['h+'][mid] + 0.258546) < 1e-4, `200 K: the base at ${s.V['h+'][mid]} V`);
});

test('islands: a base behind further layers is summed out through them, to the first bottleneck either way', () => {
  // A floating p⁺ base with a fourth layer beyond it (n/p⁺/n/n⁺, or n⁺/n/p⁺/n, or both): the
  // holes' bottleneck can lie in the third layer or the fourth, and the base's sum, taken at the
  // face of the third, read its leak as a difference of terms 1e24 larger. At 200 K, −0.2 V, a
  // cold solve left the base at 0.0002 V where its leaks balance at −0.1570 (warned; a sweep there, at 0 V, unwarned), and in
  // 61 of 360 such solves (five stacks, 150 to 300 K, six biases, cold, swept, and listed in
  // reverse) nothing was said. The sum reaches through the layers between now; a layer it
  // passes isn't summed itself; and each sum's row is its piece's middle node, not its left
  // neighbour's outside row. The check is the solver's own in no part: the leaks' balance from
  // the solution's φ alone (see the three-layer test above).
  const mat = 'GaAs', n = (N, L) => layer(mat, L, { donors: units.perCm3(N) }), p = layer(mat, 200e-9, { acceptors: units.perCm3(5e18) });
  const stacks = { right: [n(1e17, 400e-9), p, n(3e17, 600e-9), n(4e17, 300e-9)], left: [n(2e17, 200e-9), n(1e17, 400e-9), p, n(3e17, 900e-9)], both: [n(2e17, 200e-9), n(1e17, 400e-9), p, n(3e17, 600e-9), n(4e17, 300e-9)], lighter: [n(1e17, 400e-9), p, n(3e17, 600e-9), n(1e17, 300e-9)] };
  const balanced = (s, def, T) => {
    const hp = def.materials[mat].species['h+'], VT = (8.314462618 * T) / 96485.33212, { x, phi, region } = s, w = (g) => (phi[g] + hp.mu0 / 96485.33212) / VT;
    const seg = (g) => {
      const a = w(g), b = w(g + 1), dx = x[g + 1] - x[g];
      return (Math.abs(a - b) < 1e-12 ? dx * Math.exp(a) : (dx * (Math.exp(b) - Math.exp(a))) / (b - a)) / hp.cRef;
    };
    const pr = def.regions.findIndex((r) => r.name === 'base'), mid = (region.indexOf(pr) + region.lastIndexOf(pr)) >> 1, Vh = s.V['h+'];
    let IL = 0, IR = 0;
    for (let g = 0; g < x.length - 1; g++) if (x[g + 1] > x[g]) g < mid ? (IL += seg(g)) : (IR += seg(g));
    const vl = Vh[0] / VT, vr = Vh.at(-1) / VT, m = Math.max(vl, vr);
    return { base: Vh[mid], level: VT * (m + Math.log((Math.exp(vl - m) / IL + Math.exp(vr - m) / IR) / (1 / IL + 1 / IR))) };
  };
  for (const [name, T, V] of [['right', 200, -0.2], ['right', 300, -0.4], ['left', 200, -0.3], ['both', 250, -0.2], ['both', 300, 0.3], ['lighter', 200, -0.4], ['lighter', 250, -0.1]]) {
    const def = build({ T, library: [semiconductor(mat, { T })], stack: [ohmic(0), ...stacks[name].map((l) => (l === p ? { ...l, name: 'base' } : l)), ohmic(0)] });
    const cold = new Device(def), swept = new Device(def);
    cold.set({ contacts: { right: { V } } });
    let w = swept.solve();
    for (let k = 1; k <= 4; k++) {
      swept.set({ contacts: { right: { V: (V * k) / 4 } } });
      w = swept.solve();
    }
    for (const [how, s] of [['cold', cold.solve()], ['swept', w]]) {
      const { base, level } = balanced(s, def, T), { left, right } = s.terminals, at = `${name}, ${T} K, ${V} V, ${how}`;
      assert.ok(s.converged && !s.warnings.some((x) => /balance|add up/.test(x)), `${at}: ${s.warnings}`);
      assert.ok(Math.abs(base - level) < 1e-6, `${at}: the base at ${base} V, its leaks balanced at ${level} V`);
      assert.ok(Math.abs(left.current + right.current) < 1e-6 * Math.abs(left.current), `${at}: ${left.current} and ${right.current} A/m²`);
    }
  }
});

test('islands: a piece the round can\'t sum out, left unbalanced, is warned of', () => {
  // The floating base with a side port on its collector carrying 1 mA/m² (a collector contact):
  // the base's island is summed at its face (the port is in the end piece), the round fails, and
  // the cold solve stands 0.135 V off. Its terminals' miss is 8e-10 of the port's current, under
  // the Kirchhoff warning's 1e-6; the island's own balance misses by 0.94.
  const mat = 'GaAs', def = build({ T: 300, library: [semiconductor(mat)], stack: [ohmic(0), layer(mat, 0.483e-6, { name: 'n1', donors: 0.66 }), layer(mat, 0.22e-6, { name: 'p', acceptors: 8.94 }), layer(mat, 0.902e-6, { name: 'n2', donors: 3.5 }), ohmic(-0.2)], grid: { hmin: 0.5e-9, hmax: 40e-9, ratio: 1.11 } });
  def.ports = [{ name: 'side', region: 'n2', from: 0.75e-6, to: 0.9e-6, I: 1e-3, terminal: 'e-', species: { 'e-': 'equilibrium' } }];
  const s = new Device(def).solve();
  assert.ok(s.warnings.some((w) => w.includes("a floating region doesn't balance")), JSON.stringify(s.warnings));
});

test('a terminal behind a small resistance reports what its fluxes carry, where the circuit\'s reading is rougher', () => {
  // (V_src − V)/R carries a few ε of V over R: behind 1e-10 Ω·m² at −0.2 V, 4e-7 A/m², where a
  // floating base passes 1.8e-13. It read 0, and the gap between it and the fluxes (the whole
  // current) hid a wrong level 0.14 V off behind a fourth layer (found by a Fable 5.1 reviewer).
  const mat = 'GaAs', n = (N, L) => layer(mat, L, { donors: units.perCm3(N) });
  const at = (more, R) => new Device(build({ T: 300, library: [semiconductor(mat)], stack: [ohmic(0), n(1e17, 400e-9), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), n(3e17, more ? 600e-9 : 900e-9), ...(more ? [n(4e17, 300e-9)] : []), ohmic({ V: -0.2, R })] })).solve();
  const three = at(false, 1e-10), { left, right } = three.terminals;
  assert.ok(Math.abs(left.current - 1.78864e-13) < 1e-17 && Math.abs(left.current + right.current) < 1e-6 * left.current && three.warnings.length === 0, JSON.stringify(three.terminals));
  // (that fourth layer solves since 9 Nov; behind each resistance its two readings agree, or it says so)
  for (const R of [1e-10, 1e-8, 1e-6]) {
    const four = at(true, R), t = four.terminals;
    assert.ok(Math.abs(t.left.current + t.right.current) < 1e-4 * Math.abs(t.left.current) || four.warnings.some((w) => w.includes("don't add up") || w.includes("doesn't balance")), `four layers behind ${R} Ω·m²: ${t.left.current} and ${t.right.current} A/m², ${four.warnings}`);
  }
});

test('islands: a round that leaves its island less balanced, or rows satisfied out of balance, is said so', () => {
  // The four-layer base at 250 K, swept to −0.3 V in 0.05 V steps: summed at its face, the
  // base's row reads noise and Newton satisfied it 30 mV off (−0.2157 V; cold and the leaks'
  // balance give −0.2461), its currents 13% apart at 3e-18 A/m², under the Kirchhoff floor, and
  // nothing was said (found by a Fable 5.1 reviewer re-checking 5b8f56e).
  const mat = 'GaAs', n = (N, L) => layer(mat, L, { donors: units.perCm3(N) });
  const d = new Device(build({ T: 250, library: [semiconductor(mat, { T: 250 })], stack: [ohmic(0), n(1e17, 400e-9), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), n(3e17, 600e-9), n(4e17, 300e-9), ohmic(0)] }));
  let s;
  for (let k = 1; k <= 6; k++) {
    d.set({ contacts: { right: { V: (-0.3 * k) / 6 } } });
    s = d.solve();
  }
  const mid = (s.region.indexOf(1) + s.region.lastIndexOf(1)) >> 1;
  assert.ok(Math.abs(s.V['h+'][mid] + 0.2461) < 1e-3 || s.warnings.some((w) => w.includes("doesn't balance") || w.includes("don't add up")), `the base at ${s.V['h+'][mid]} V: ${s.warnings}`);
  // A round undone is warned of only where the terminals miss too (relatively; that base's
  // 3e-18 A/m² is under the warning's floor): 22 of the law sweep's devices, each agreeing with
  // its transient, undid a round whose balance was noise before and after, and were warned.
  for (const k of [361, 452, 455, 462]) assert.deepEqual(new Device(device(choice(k))).solve().warnings.filter((w) => w.includes("doesn't balance")), [], `law sweep ${k}`);
});

test('a device with no weakly held region reads one terminal\'s round-off as such: an insertion electrode at rest, its electrons\' 1e-13 A/m² beside Li⁺\'s 0', () => {
  // (docs/statistics.md's OCV example; c7ed8e2 warned of it, 1.8.0 didn't)
  const RT = GAS_CONSTANT * 298.15, x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];
  const E = x.map((v) => 0.4 - (RT / FARADAY) * Math.log(v / (1 - v)));
  const dev = new Device({
    species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: { host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }] } },
    regions: [{ material: 'host', length: 1e-6 }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { V: 0.4, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
  });
  let seen = 0;
  for (const V of [0.4, 0.3]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged);
    seen = Math.max(seen, Math.abs(sol.terminals.right.current));
    assert.ok(!sol.warnings.some((w) => /don't add up/.test(w)), sol.warnings.join('\n'));
  }
  assert.ok(seen > 1e-15, `${seen}: the round-off this is about`);
});
