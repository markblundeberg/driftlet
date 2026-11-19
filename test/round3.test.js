import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous, bath, build, layer } from '../src/kit.js';
import { sourceAt, nextBreakpoint } from '../src/device.js';

// Faults found in the time stepper and the sources by two outside reading passes (OpenAI's Sol
// on src/timeSteps.js, Astra on the definition's normalization), each confirmed by a probe
// against a second route before it was fixed.

// The blocking cell: KCl between two capacitive contacts, a terminal driven as given.
const blocking = (right, left = { V: 0 }) => {
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2;
  return new Device({
    species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { ...left, phi: { type: 'capacitive', C }, zeroCharge: 0 }, right: { ...right, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
};
const last = (a) => a[a.length - 1];

test("a repeating waveform reads the after side at every edge a step lands on, and lands on each", () => {
  // ((0.5 % 0.1) is 0.0999…98: a square wave of period 0.1 s read 1183 of its first 2000 landed
  // edges from before them, and solve() there the wrong drive)
  const sq = { t: [0, 0, 0.05, 0.05, 0.1], values: [0, 1, 1, 0, 0], repeat: true };
  let t = 0;
  for (let i = 0; i < 2000; i++) {
    t = nextBreakpoint(sq, t);
    const rising = Math.round(t / 0.05) % 2 === 0;
    assert.equal(sourceAt(sq, t), rising ? 1 : 0, `after the edge at ${t}`);
    assert.equal(sourceAt(sq, t, true), rising ? 0 : 1, `before the edge at ${t}`);
  }
  assert.ok(Math.abs(t - 100) < 1e-9);
  // a 1 µs square wave late in a run: every edge, not none (past 1e-12 of t, all were passed)
  const fast = { t: [0, 0, 5e-7, 5e-7, 1e-6], values: [0, 1, 1, 0, 0], repeat: true };
  assert.ok(Math.abs(nextBreakpoint(fast, 2e6) - 2e6 - 5e-7) < 1e-8);
  // a long period: its first step isn't lost to (1 % T) + T rounding to T
  const long = { t: [0, 1, 1, 1e16], values: [0, 0, 1, 1], repeat: true };
  assert.equal(sourceAt(long, 1), 1);
  assert.equal(nextBreakpoint(long, 0), 1);
  // a period below the time's resolution: refused, not a hang
  assert.throws(() => nextBreakpoint({ t: [0, 1], values: [0, 1], repeat: true }, 2 ** 53), /resolution/);
});

test("a first fixed step onto a current's edge posts the before side's charge", () => {
  // (the currents were read before the step was marked stepped: the far side's I·dt, for good)
  const I = { t: [0, 0.01, 0.01], values: [0, 0, 1] };
  for (const method of ['be', 'bdf2']) {
    const d = blocking({ I });
    d.solve();
    d.step(0.01, { method });
    const a = d.advance(0.02);
    assert.ok(Math.abs(last(a.trace.terminals.right.charge) - a.plates.right.charge) < 1e-6 * a.plates.right.charge, `${method}: ledger ${last(a.trace.terminals.right.charge)}, plate ${a.plates.right.charge}`);
  }
});

test("a call ending a short tail after an edge steps it, rather than snapping the clock past it", () => {
  // (within 1e-10 of the time it snapped: 5e-6 s at 1e5 s, the current through it lost for good)
  const T = 1e5, tail = 5e-6, d = blocking({ I: { t: [0, T, T], values: [0, 0, 1] } });
  d.solve();
  const a = d.advance(T + tail, { dt0: T });
  assert.ok(Math.abs(a.plates.right.charge / tail - 1) < 1e-3, `${a.plates.right.charge} for ${tail}`);
  assert.equal(d.solver.sourceTime, d.solver.time);
});

test("a step's halves restart the order on an edge they start on, at any depth", () => {
  // (only the outermost call did: a half starting on the edge crossed it by BDF2, a third short)
  const d = blocking({ I: { t: [0, 0.01, 0.01], values: [0, 0, 1] } });
  d.solve();
  d.step(0.005, { method: 'bdf2' });
  const s = d.solver, real = s.newton;
  let calls = 0;
  s.newton = function (dt, opts) {
    return ++calls === 1 ? { converged: false, iterations: 0, history: [] } : real.call(this, dt, opts);
  };
  const r = d.step(0.01, { method: 'bdf2' });
  s.newton = real;
  assert.ok(Math.abs(r.plates.right.charge / 0.005 - 1) < 1e-3, `${r.plates.right.charge} for 0.005`);
});

test('advance() and step() refuse what they can\'t mean, before touching the state', () => {
  const d = blocking({ V: 0 });
  d.solve();
  for (const opts of [{ dtMax: -0.1 }, { dtMax: 0 }, { tol: NaN }, { tol: -1 }, { maxSteps: 1.5 }]) assert.throws(() => d.advance(1, opts), /advance\(\)/, JSON.stringify(opts));
  assert.throws(() => d.advance('1'), /tEnd/);
  d.advance(1);
  const before = d.solver.time, sol = d.step(1e-17);
  assert.equal(d.solver.time, before);
  assert.ok(sol.warnings.some((w) => /resolution/.test(w)));
});

test("set() keeps its own copies of a patch's drive and links, and the constructor of the definition", () => {
  // (the caller's objects were put in: an edit to one made the next set() of it a no-op, and
  // leaked in at an unrelated rebuild; the constructor's model kept the caller's arrays)
  const d = blocking({ V: 0 });
  d.solve();
  const w = { t: [0, 1], values: [0, 1] };
  d.set({ contacts: { right: { V: w } } });
  w.values[1] = 2;
  assert.deepEqual(d.def.contacts.right.V.values, [0, 1]);
  d.set({ contacts: { right: { V: w } } });
  assert.deepEqual([...d.model.contacts.right.drive.src.values], [0, 2]);
  const def = { species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } }, statistics: [{ type: 'redlich-kister', species: ['X'], cMax: 3, A: [100] }] } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
    grid: { minCells: 8 } };
  const e = new Device(def), mu = e.solution().mu.X[0];
  def.materials.f.statistics[0].A[0] = 1000;
  assert.equal(new Device(e.def).solution().mu.X[0], mu);
  assert.equal(e.solution().mu.X[0], mu);
});

// Which nodes a port or a reaction acts on, and what counts as an anchor (Astra's third pass, and
// Sol's on the same rows).

// X in two unit regions joined by an equilibrium face, closed at both ends, with ports holding it.
const faced = (ports, cells = 8) => ({
  species: [{ name: 'X', z: 0 }],
  materials: { f: { epsr: 4, species: { X: { D: 1, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'a', material: 'f', length: 1, c0: { X: 1 } }, { name: 'b', material: 'f', length: 1, c0: { X: 1 } }],
  contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
  ports,
  grid: { hmin: 1 / cells, hmax: 1 / cells },
});
const holdX = (name, region, from, to, mu) => ({ name, region, from, to, V: 0, species: { X: { type: 'equilibrium', mu } } });

test('two ports either side of an equilibrium face are refused, and one there alone solves', () => {
  // (a's last node is held at b's first: both held one row, the last listed won it, and the other
  // reported 2400/RT as its flux out of a closed device, converged and unwarned)
  const pa = holdX('pa', 'a', 1, 1, 0), pb = holdX('pb', 'b', 0, 0, 2400);
  for (const ports of [[pa, pb], [pb, pa], [{ ...pa, from: 0.5 }, pb]]) assert.throws(() => new Device(faced(ports)), /either side of the face at x = 1\.000 m/);
  // alone, at a level the start already has: its steady solve was "exactly singular" (b, where the
  // hold sits, was taken for an island)
  for (const [from, mu] of [[1, 0], [0.9, 0], [1, 2400]]) {
    const d = new Device(faced([holdX('p', 'a', from, 1, mu)])), sv = d.solver;
    // (and a, the island now, doesn't sum the row the port took at b's first node as a balance)
    const took = new Set(Array.from(sv.model.ports[0].held[0], (g) => sv.blockOfNode[g] * sv.M + 1));
    assert.ok(sv.islands.every((isl) => isl.outside.every((o) => !took.has(o))), JSON.stringify(sv.islands.map((isl) => isl.outside)));
    const s = d.solve(), c = Math.exp(mu / (GAS_CONSTANT * 298.15));
    assert.ok(s.converged, `[${from}, 1] at ${mu}: ${s.warnings}`);
    assert.ok(s.c.X.every((v) => Math.abs(v / c - 1) < 1e-9) && Math.abs(s.ports[0].flux.X) < 1e-12, `${s.c.X[0]} ${s.ports[0].flux.X}`);
  }
});

test("an immobile species keeps its amount at a port's nodes the port doesn't act on", () => {
  // (each node a port's window listed counted as touched: an equilibrium port's fringe nodes, and
  // an electrode's where its area is 0, had no law, and S⁺ there moved with φ; 7× here)
  const fringe = (V) => new Device({
    species: [{ name: 'S+', z: 1 }, { name: 'N-', z: -1 }],
    materials: { f: { epsr: 1e12, species: { 'S+': { D: 0, cRef: 1, mu0: 0 }, 'N-': { D: 1, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1, c0: { 'S+': 1, 'N-': 1 } }],
    contacts: { left: { V, terminal: 'N-', species: { 'N-': 'equilibrium', 'S+': 'blocked' }, phi: 'bulk' }, right: { phi: 'neutral', species: { 'N-': 'blocked', 'S+': 'blocked' } } },
    ports: [{ name: 'p', region: 'a', from: 0.35, to: 0.65, V: 0, terminal: 'S+', species: { 'S+': 'equilibrium' } }],
    grid: { minCells: 4 },
  });
  const d = fringe(0);
  d.solve();
  d.set({ contacts: { left: { V: 0.1 } } });
  const s = d.solve(), t = fringe(0);
  t.solve();
  t.set({ contacts: { left: { V: 0.1 } } });
  const r = t.advance(10); // (D = 0 and nothing at the fringe: the transient keeps it exactly)
  assert.ok(s.converged && r.converged);
  for (const g of [0, 1, 3, 4]) assert.ok(Math.abs(s.c['S+'][g] - 1) < 1e-9, `node ${g}: ${s.c['S+'][g]}`);
  assert.ok(Math.abs(s.c['S+'][2] / r.c['S+'][2] - 1) < 1e-6, `${s.c['S+'][2]} ${r.c['S+'][2]}`);
  // an electrode consuming S over the right half only: at the left half's nodes S stays at 0.5,
  // where it followed the bath (0.0004 to 49)
  const el = (V) => new Device({
    species: [{ name: 'S', z: 1 }, { name: 'P', z: 1 }, { name: 'N', z: -1 }, { name: 'e-', z: -1 }],
    materials: { w: { epsr: 0, species: { S: { D: 0, cRef: 1, mu0: 0 }, P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'w', length: 1, c0: { S: 0.5, P: 0.5, N: 1 } }],
    contacts: { left: { V, terminal: 'P', phi: 'neutral', species: { P: 'equilibrium', N: { type: 'equilibrium', offset: 0 }, S: 'blocked' } }, right: { phi: 'neutral', species: { P: 'blocked', N: 'blocked', S: 'blocked' } } },
    ports: [{ name: 'el', region: 'a', V: 0, terminal: 'e-', area: { x: [0, 0.5, 0.5000001, 1], values: [0, 0, 1, 1] }, reactions: [{ equation: 'S + e- = Q', fixed: { Q: 0 }, k0: 1e-6, alpha: 0.5 }] }],
    grid: { hmin: 1 / 8, hmax: 1 / 8 },
  });
  for (const V of [0.2, -0.2]) {
    const sol = el(V).solve();
    assert.ok(sol.converged);
    for (let g = 0; g < 4; g++) assert.ok(Math.abs(sol.c.S[g] - 0.5) < 1e-9, `${V} V, node ${g}: ${sol.c.S[g]}`);
  }
});

test("a bulk reaction acts only where its rate profile isn't 0", () => {
  // (by material: X⁺, immobile, alternated ±25% about 1 where the profile was 0; and a sealed X with
  // a profile 0 over the whole device lost its amount's law, its steady system singular)
  const L = 1e-7, sp = { 'X+': { D: 0, cRef: 1, mu0: 0 }, 'A-': { D: 1e-9, cRef: 1, mu0: 0 } };
  const s = new Device({
    species: [{ name: 'X+', z: 1 }, { name: 'A-', z: -1 }],
    contacts: { left: { V: 0.05, terminal: 'A-', species: { 'A-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, terminal: 'A-', species: { 'A-': 'equilibrium' }, phi: 'bulk' } },
    materials: { f: { epsr: 10, species: sp } },
    regions: [{ name: 'a', material: 'f', length: L, c0: { 'X+': 1, 'A-': 1 }, grid: { minCells: 16 } }],
    bulkReactions: [{ equation: 'X+ + A- = S', fixed: { S: 2000 }, kf: { f: { x: [0.5 * L - 1e-12, 0.5 * L], values: [1, 0] } } }],
  }).solve();
  assert.ok(s.converged);
  s.x.forEach((x, g) => x > 0.55 * L && assert.ok(Math.abs(s.c['X+'][g] - 1) < 1e-9, `x = ${x}: ${s.c['X+'][g]}`));
  const sealed = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1, c0: { X: { x: [0, 1], values: [0.5, 1.5] } } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
    bulkReactions: [{ equation: 'X = S', fixed: { S: GAS_CONSTANT * 298.15 * Math.log(5) }, kf: { f: { x: [0, 1, 2], values: [0, 0, 1] } } }],
    grid: { minCells: 8 },
  }).solve();
  assert.ok(sealed.converged, `${sealed.warnings}`);
  assert.ok(sealed.c.X.every((v) => Math.abs(v - 1) < 1e-9), `${sealed.c.X}`);
});

test("a port's capacitance anchors φ on a dielectric", () => {
  // (two gates at 0 and 1 V on a charge-free layer: φ NaN and both charges 0, unwarned, where
  // Σσ = 0 between closed ends puts φ at 0.5 V and the charges at ∓0.5 C/m²)
  const s = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 1, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
    ports: [{ name: 'g0', region: 'a', V: 0, area: 1, capacitance: { C: 1 } }, { name: 'g1', region: 'a', V: 1, area: 1, capacitance: { C: 1 } }],
    grid: { minCells: 4 },
  }).solve();
  assert.ok(s.converged && s.phi.every((p) => Math.abs(p - 0.5) < 1e-9), `${s.phi}`);
  assert.ok(Math.abs(s.ports[0].charge + 0.5) < 1e-9 && Math.abs(s.ports[1].charge - 0.5) < 1e-9, `${s.ports[0].charge} ${s.ports[1].charge}`);
});

test('a bare-metal reaction in a port without a surface runs on its own bare metal, beside another\'s film', () => {
  // (it read the overlapping port's θ₀: half its rate at θ = 0.5, 1e-5 of it near full; the iron
  // fixtures of test/coverage.test.js)
  const G = 1e-5, L = 10e-6, area = 1e4;
  const cell = (ports) => {
    const lib = aqueous(['Na+', 'OH-', 'Fe2+', 'Fe3+'], { epsr: 0 });
    lib.species.push({ name: 'e-', z: -1 });
    return new Device(build({ library: [lib], stack: [{}, layer('water', L, { name: 'film' }), bath({ 'Na+': 100, 'OH-': 105, 'Fe2+': 1, 'Fe3+': 1 }, 'OH-')], ports, grid: { hmin: L / 20, hmax: L / 20 } }));
  };
  const probe = cell([{ name: 'p', region: 'film', V: 0, terminal: 'e-', area, surface: { OHads: { mu0: 0, capacity: G } }, reactions: [{ equation: 'OH- = OHads + e-', k0: 1e-4, alpha: 0.5, bare: true }] }]).solve();
  const muFe = probe.mu['Fe2+'][3] + 2 * FARADAY * 0.15, mu0 = muFe + 2 * probe.mu['OH-'][3] + 2 * FARADAY * 0.05;
  const film = { equation: 'Fe(s) + 2 OH- = Fe(OH)2 + 2 e-', fixed: { 'Fe(s)': muFe }, k0: 1e-4, alpha: 0.5, bare: true };
  const A = { name: 'A', region: 'film', V: 0.2, terminal: 'e-', area, surface: { 'Fe(OH)2': { mu0, capacity: G } }, reactions: [film] };
  const irons = [
    ['Butler–Volmer', { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': muFe }, k0: 1e-6, alpha: 0.5, bare: true }],
    ['Marcus–Hush–Chidsey', { equation: 'Fe3+ + e- = Fe2+', k0: 1e-6, lambda: 50000, bare: true }],
  ];
  for (const [kind, iron] of irons) {
    const B = { name: 'B', region: 'film', V: 0.05, terminal: 'e-', area, reactions: [iron] };
    const alone = cell([B]).solve(), s = cell([A, B]).solve();
    assert.ok(alone.converged && s.converged && s.ports[0].coverage['Fe(OH)2'][3] > 0.99, `${kind}: A covered`);
    const [r0, r1] = [alone.ports[0].rates[0][3], s.ports[1].rates[0][3]];
    assert.ok(r0 !== 0 && Math.abs(r1 / r0 - 1) < 1e-6, `${kind}: ${r1} beside A's film, ${r0} alone`);
  }
});

// Two strictly neutral salt regions behind a blocked neutral face, the left one held by port A.
const cut = (B) => ({
  species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }],
  materials: { w: { epsr: 0, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'r1', material: 'w', length: 1, c0: { P: 1, N: 1 } }, { name: 'r2', material: 'w', length: 1, c0: { P: 2, N: 2 } }],
  interfaces: [{ phi: 'neutral', species: { P: 'blocked', N: 'blocked' } }],
  contacts: { left: { phi: 'neutral', species: { P: 'blocked', N: 'blocked' } }, right: { phi: 'neutral', species: { P: 'blocked', N: 'blocked' } } },
  ports: [{ name: 'A', region: 'r1', V: 0, terminal: 'P', species: { P: { type: 'conductance', G: 1 } } }, { name: 'B', region: 'r2', ...B, terminal: 'P', species: { P: { type: 'conductance', G: 1 } } }],
  grid: { minCells: 4 },
});

test('a port driven by a current anchors nothing, at construction or set to it', () => {
  // (B at I: 0, the only tie of r2: accepted, and a steady solve settled at the cold start's gauge)
  assert.throws(() => new Device(cut({ I: 0 })), /regions\[1\] \(r2\) is anchored electrostatically only by ports\[1\], driven by a current/);
  const d = new Device(cut({ V: 0 }));
  assert.ok(d.solve().converged);
  assert.throws(() => d.set({ ports: { B: { I: 0 } } }), DeviceError);
  assert.equal(d.model.ports[1].drive.kind, 'V', 'left as it was');
  // (behind a resistance it's held: its level is the source's)
  assert.doesNotThrow(() => new Device(cut({ V: 0, R: 10 })));
});

test("a face reaction that moves no charge across doesn't couple the sides' φ", () => {
  // (P + N on each side passed for a coupling: the far side floated, and every solve failed with
  // nothing specific shown; one that moves P⁺ across does couple them)
  const sp = { D: 1e-9, cRef: 1, mu0: 0 };
  const def = (equation) => ({
    species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }],
    materials: { m: { epsr: 0, species: { P: sp, N: sp } } },
    regions: [{ name: 'a', material: 'm', length: 1e-6, c0: { P: 1, N: 1 } }, { name: 'b', material: 'm', length: 1e-6, c0: { P: 2, N: 2 } }],
    interfaces: [{ species: { P: 'blocked', N: 'blocked' }, reactions: [{ equation, k0: 1e-3 }] }],
    contacts: { left: { bath: { c: { P: 1, N: 1 } } }, right: { phi: 'neutral' } },
    grid: { minCells: 8 },
  });
  assert.throws(() => new Device(def('P(left) + N(left) = P(right) + N(right)')), /regions\[1\] \(b\) is electrostatically floating/);
  const s = new Device(def('P(left) = P(right)')).solve();
  assert.ok(s.converged && Math.abs(s.phi.at(-2) + ((GAS_CONSTANT * 298.15) / FARADAY) * Math.log(2)) < 1e-6, `${s.phi.at(-2)}`);
});

test("a port that acts on no node is refused, so it never passes for a current's return", () => {
  // (each passed for the held terminal the left contact's current returned through, and every
  // solve failed "from the starting state")
  const L = 1e-6, sp = { D: 1e-9, cRef: 1, mu0: 0 };
  const def = (epsr, port, right = { phi: 'neutral' }) => ({
    species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }],
    materials: { m: { epsr, species: { P: sp, N: sp } } },
    regions: [{ name: 'a', material: 'm', length: L, c0: { P: 1, N: 1 } }],
    contacts: { left: { I: 1e-3, terminal: 'P', species: { P: 'equilibrium', N: 'blocked' }, phi: 'bulk' }, right },
    ports: [port],
    grid: { minCells: 8 },
  });
  assert.throws(() => new Device(def(0, { name: 'g', region: 'a', V: 0, area: { x: [0, L], values: [0, 0] }, capacitance: { C: 1 } })), /ports\[0\]: it acts on no node of the device \(its electrode's area is 0 throughout its window\)/);
  assert.throws(() => new Device(def(1e3, { name: 'g', region: 'a', from: L, to: L, V: 0, area: 1, capacitance: { C: 1 } }, { I: 0, phi: 'pinned', zeroCharge: 0 })), /its capacitance is only in the box of an end node whose contact sets φ/);
  assert.throws(() => new Device(def(0, { name: 'g', region: 'a', from: 0, to: 0, V: 0, terminal: 'P', species: { P: 'equilibrium' } })), /the level it holds is at a device end node, which its contact links/);
  // (the same port inside the region is the return the current needs)
  const s = new Device(def(0, { name: 'g', region: 'a', from: L / 2, to: L / 2, V: 0, terminal: 'P', species: { P: 'equilibrium' } })).advance(1e-3);
  assert.ok(s.converged && Math.abs(s.ports[0].current + 1e-3) < 1e-9, `${s.ports[0].current}`);
});
