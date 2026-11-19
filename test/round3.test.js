import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, GAS_CONSTANT } from '../src/index.js';
import { aqueous, photogeneration } from '../src/kit.js';
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

// Input hygiene, profiles and grids (Astra's third pass over the definition's normalization).
// A neutral film of one species, X, sealed at both ends.
const film = (over = {}) => ({
  species: [{ name: 'X', z: 0 }],
  materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
  contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
  grid: { minCells: 8 },
  ...over,
});
const refused = (def, re) => assert.throws(() => new Device(def), (e) => e instanceof DeviceError && re.test(e.message), re.source);
const RT = GAS_CONSTANT * 298.15;
const hole = (a, k) => { delete a[k]; return a; }; // (an array with a hole at k, as [0, , 2] writes it)

test('a malformed definition is a DeviceError naming its path: ports not an array, an area not a number, a hole in an array', () => {
  // (ports: {} and area: null were TypeErrors; a hole in a waveform's values was accepted as NaN,
  // and one in species said 'X' was in no material)
  refused(film({ ports: {} }), /^ports must be an array/);
  const gate = (area) => film({ ports: [{ name: 'g', region: 'a', V: 0, area, capacitance: { C: 1e-3 } }] });
  refused(gate(null), /ports\[0\]\.area must be the electrode's area per volume/);
  refused(gate('1'), /ports\[0\]\.area must be .* got "1"/);
  const held = (V) => { const d = film(); d.contacts.left = { V, phi: 'neutral', species: { X: 'blocked' } }; return d; };
  refused(held({ t: [0, 1], values: hole([0, 1], 1) }), /contacts\.left\.V\.values\[1\] is a hole/);
  refused(held({ t: hole([0, 1, 2], 1), values: [0, 1, 2] }), /contacts\.left\.V\.t\[1\] is a hole/);
  const def = film();
  def.species = hole([{ name: 'Y', z: 0 }, { name: 'X', z: 0 }], 0);
  refused(def, /^species\[0\] is a hole/);
  const prof = film();
  prof.regions[0].c0 = { X: { x: hole([0, 5e-7, 1e-6], 1), values: [1, 1, 1] } };
  refused(prof, /regions\[0\]\.c0\.X\.x\[1\] is a hole/);
  // and through set(), a drive patch alone
  const d = new Device(held(0));
  d.solve();
  assert.throws(() => d.set({ contacts: { left: { V: { t: [0, 1], values: hole([0, 1], 1) } } } }), /values\[1\] is a hole/);
});

test('region and port names are strings: a region named 0 is not index 0, and ports 1 and "1" are not two', () => {
  // (a port's region: 0 went to index 0, not the region named 0; set({ ports: { 1: … } }) edited
  // both ports named 1 and '1')
  const two = film();
  two.regions = [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }, { name: 0, material: 'f', length: 1e-6, c0: { X: 1 } }];
  refused(two, /regions\[1\]\.name must be a non-empty string, got 0/);
  const port = (name, from) => ({ name, region: 'a', from, to: from + 1e-7, species: { X: { type: 'exchange', k: 1, mu: 0 } } });
  refused(film({ ports: [port(1, 1e-7), port('1', 6e-7)] }), /ports\[0\]\.name must be a non-empty string, got 1/);
  assert.equal(new Device(film({ ports: [port('1', 1e-7)] })).model.ports[0].name, '1');
});

test('a charge number past what the solver holds is refused, not wrapped to 0 or flipped', () => {
  // (|z| ≥ 2³¹ reached the solver's Int32 z as 0, or −2³¹ for +2³¹: "did not converge")
  for (const z of [2147483648, 4294967296, 9007199254740992]) refused(film({ species: [{ name: 'X', z }] }), /species\[0\]\.z \(charge number\) must be at most 1e6/);
});

test("species named like Object.prototype's properties take part in equations as any other, at a face and in bulk", () => {
  // ('__proto__ = X' at a face lost its left side and ran as 0 = X: the total fell from 2 to 1.67,
  // converged; 'toString = X' was refused for a coefficient "function toString() …-1"; in bulk the
  // definition's own copy dropped '__proto__'. Analytic: with μ°_X = 1000 J/mol, K = e^(−1000/RT),
  // the steady state has 2/(1 + K) of the one and 2K/(1 + K) of X.)
  const K = Math.exp(-1000 / RT);
  const named = (def, name) => JSON.parse(JSON.stringify(def).replaceAll('NAME', name)); // (JSON: '__proto__' an own key)
  for (const name of ['__proto__', 'toString', 'constructor']) {
    const face = named({
      species: [{ name: 'NAME', z: 0 }, { name: 'X', z: 0 }],
      materials: { A: { epsr: 4, species: { NAME: { D: 1e-9, cRef: 1, mu0: 0 } } }, B: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 1000 } } } },
      regions: [{ name: 'a', material: 'A', length: 1e-6, c0: { NAME: 1 } }, { name: 'b', material: 'B', length: 1e-6, c0: { X: 1 } }],
      interfaces: [{ phi: 'neutral', reactions: [{ equation: 'NAME = X', k0: 1e-6, alpha: 0.5 }] }],
      contacts: { left: { phi: 'neutral', species: { NAME: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
      grid: { minCells: 8 },
    }, name);
    const s = new Device(face).solve();
    assert.ok(s.converged);
    assert.ok(Math.abs(s.c[name][0] - 2 / (1 + K)) < 1e-4, `${name} at a face: ${s.c[name][0]}`);
    assert.ok(Math.abs(s.c.X.at(-1) - (2 * K) / (1 + K)) < 1e-4, `X beside ${name}: ${s.c.X.at(-1)}`);
    const bulk = named({
      species: [{ name: 'NAME', z: 0 }, { name: 'X', z: 0 }],
      materials: { f: { epsr: 4, species: { NAME: { D: 1e-9, cRef: 1, mu0: 0 }, X: { D: 1e-9, cRef: 1, mu0: 1000 } } } },
      regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { NAME: 1, X: 1 } }],
      contacts: { left: { phi: 'neutral', species: { NAME: 'blocked', X: 'blocked' } }, right: { phi: 'neutral', species: { NAME: 'blocked', X: 'blocked' } } },
      bulkReactions: [{ equation: 'NAME = X', kf: { f: 1 } }],
      grid: { minCells: 8 },
    }, name);
    const b = new Device(bulk).solve();
    assert.ok(Math.abs(b.c[name][3] - 2 / (1 + K)) < 1e-4, `${name} in bulk: ${b.c[name][3]}`);
  }
});

test('a field where it means nothing is refused by its presence, not passed for having no effect', () => {
  // (each was accepted and normalised exactly as without it)
  const closed = film();
  closed.contacts.left = { zeroCharge: 0.4 };
  refused(closed, /contacts\.left\.zeroCharge: only a capacitive or pinned φ law takes an alignment/);
  const supply = (extra) => film({ ports: [{ name: 'o2', region: 'a', species: { X: { type: 'exchange', k: 1e-3, mu: 0 } }, ...extra }] });
  refused(supply({ terminal: 'X' }), /ports\[0\]\.terminal: this port passes no current/);
  refused(supply({ V: 0.3 }), /ports\[0\]\.V: this port passes no current/);
  refused(supply({ R: 10 }), /ports\[0\]\.R: this port passes no current/);
  const d = new Device(supply({}));
  d.solve();
  assert.throws(() => d.set({ ports: { o2: { V: 0.3 } } }), /ports\[0\]\.V: this port passes no current/);
  const rx = (fixed) => film({ bulkReactions: [{ equation: 'X = S', fixed, kf: { f: 1 } }] });
  refused(rx({ S: 0, Ss: 5000 }), /bulkReactions\[0\]\.fixed\.Ss: not a participant of the reaction/);
  refused(rx({ S: 0, X: 5000 }), /bulkReactions\[0\]\.fixed\.X: 'X' is a species of the device/);
  assert.equal(new Device(rx({ S: 0 })).model.reactions[0].fixedA, 0);
  refused(film({ geometry: { type: 'planar', r0: 1e-6 } }), /geometry\.r0: a planar device has no centre/);
  const cu = (reg) => ({
    species: [{ name: 'e-', z: -1 }], materials: { Cu: { conductor: { species: 'e-', conductivity: 1e6 } } },
    regions: [{ material: 'Cu', length: 1e-6, ...reg }],
    contacts: { left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 1e-3, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
  });
  refused(cu({ fixedCharge: 0 }), /regions\[0\]\.fixedCharge: a conductor region is neutral in bulk/);
  refused(cu({ fixedCharge: { x: [0, 1e-6], values: [0, 0] } }), /regions\[0\]\.fixedCharge: a conductor region/);
  refused(cu({ velocity: 0 }), /regions\[0\]\.velocity: a conductor region has no flow/);
  refused(cu({ mixing: 0 }), /regions\[0\]\.mixing: a conductor region has no flow/);
  assert.ok(new Device(cu({})).solve().converged);
});

test('a narrow port window holds its width, a point window (to === from) its node\'s whole box', () => {
  // (a window 5.6e-16 m wide passed for one of no width and exchanged a whole box, 0.25 m of it
  // on 4 cells of a 1 m film: a steady flux of 0.2235 where 5.6e-16 is right)
  const held = (from, to) => {
    const d = new Device({
      species: [{ name: 'X', z: 0 }], materials: { f: { epsr: 4, species: { X: { D: 1, cRef: 1, mu0: 0 } } } },
      regions: [{ name: 'a', material: 'f', length: 1, c0: { X: 1 } }],
      contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
      ports: [{ name: 'p', region: 'a', from, to, species: { X: { type: 'exchange', k: 1, mu: RT } } }],
      grid: { minCells: 4 },
    });
    const p = d.model.ports[0];
    return p.nodes.reduce((s, g, w) => s + d.model.grid.vol[g] * p.fill[w], 0);
  };
  assert.equal(held(0.5, 0.5), 0.25);
  for (const [from, to] of [[0.5, 0.5 + 5e-16], [0.6, 0.6 + 5e-16], [0.5, 0.5 + 2e-15]]) assert.ok(Math.abs(held(from, to) / (to - from) - 1) < 1e-6, `[${from}, ${to}]: ${held(from, to)}`);
});

test('grid: {} (or options set to undefined) keeps the default grading', () => {
  // (any grid object turned it off: 8 uniform cells for 54 graded ones, unwarned)
  const cells = (over) => new Device(film(over)).model.grid.nNodes - 1;
  const graded = cells({ grid: undefined });
  assert.ok(graded > 40, `${graded}`);
  assert.equal(cells({ grid: {} }), graded);
  assert.equal(cells({ grid: { hmax: undefined } }), graded);
  const reg = film({ grid: undefined });
  reg.regions[0].grid = {};
  assert.equal(new Device(reg).model.grid.nNodes - 1, graded);
  assert.equal(cells({ grid: { minCells: 8 } }), 8);
});

test("a thin shell far from a sphere's or cylinder's centre has its volume and area, not round-off", () => {
  // (as differences of cubes and squares: a nm shell at r = 1e6 m was 37% off, converged and
  // unwarned; at 1e16 m its volumes were 0 and the solve failed)
  for (const [type, r0, L] of [['spherical', 1e6, 1e-9], ['spherical', 1e16, 1], ['cylindrical', 1e6, 1e-9], ['cylindrical', 1e16, 1]]) {
    const d = new Device(film({ geometry: { type, r0 }, regions: [{ name: 'a', material: 'f', length: L, c0: { X: 1 } }], grid: { minCells: 4, hmax: L / 4 } }));
    const g = d.model.grid, A = (x) => (type === 'spherical' ? 4 * Math.PI * (r0 + x) * (r0 + x) : 2 * Math.PI * (r0 + x));
    // (to first order in h/r, ~1e-15 here: A at a box's middle times its width, a segment's A at its middle)
    for (let k = 1; k < g.nNodes - 1; k++) assert.ok(Math.abs(g.vol[k] / (A(g.x[k]) * (L / 4)) - 1) < 1e-9, `${type} r0 ${r0}: vol[${k}] ${g.vol[k]}`);
    for (let k = 0; k < g.nNodes - 1; k++) assert.ok(Math.abs(g.segArea[k] / A((g.x[k] + g.x[k + 1]) / 2) - 1) < 1e-9, `${type} r0 ${r0}: segArea[${k}] ${g.segArea[k]}`);
    assert.ok(d.solve().converged);
  }
});

test("a constant rate profile is its value in every box, wherever the table's point lies", () => {
  // (taken as the difference of two integrals from the table's first point: from x = −1e16 m, a
  // constant 1 averaged to 0 over every box of a µm film)
  const def = film({ bulkReactions: [{ equation: 'photon = X', fixed: { photon: 3e5 }, kf: { f: { x: [-1e16], values: [1] } } }] });
  assert.deepEqual([...new Device(def).solver.rxs[0].kfNode], new Array(9).fill(1));
});

test('profile warnings: a rate profile zero over all of its material is warned of; a grading as written by the device is not', () => {
  // (light written past the device's end generated nothing, unwarned; a region on [0.1, 1.1] graded
  // by the device's x from a table on [0, 1] was told it wasn't graded as written)
  const off = (kf) => new Device(film({ bulkReactions: [{ equation: 'photon = X', fixed: { photon: 3e5 }, kf: { f: kf } }] })).solution().warnings;
  const nowhere = /kf\.f: the profile .* is zero wherever f is, so the reaction runs nowhere/;
  assert.ok(off({ x: [2e-6, 3e-6, 4e-6], values: [0, 1, 0] }).some((w) => nowhere.test(w)));
  assert.ok(off({ x: [-3e-6, -2e-6, -1e-6], values: [0, 1, 0] }).some((w) => nowhere.test(w)));
  assert.ok(off(photogeneration({ material: 'f', flux: 1, alpha: 1e7, mu: 3e5, from: 2e-6, to: 3e-6 }).kf.f).some((w) => nowhere.test(w)));
  assert.ok(!off({ x: [0.5e-6, 1e-6, 1.5e-6], values: [0, 1, 0] }).some((w) => /kf\.f/.test(w)));
  const two = (mu0) => film({
    regions: [{ name: 'a', material: 'f', length: 0.1, c0: { X: 1 } }, { name: 'b', material: 'f', length: 1, c0: { X: 1 }, grade: { X: { mu0 } } }],
    interfaces: [{ phi: 'neutral', species: { X: 'equilibrium' } }], grid: { minCells: 10 },
  });
  const warned = (mu0) => new Device(two(mu0)).solution().warnings.some((w) => /isn't graded as written/.test(w));
  assert.ok(!warned({ x: [0, 1], values: [0, RT] }));
  assert.ok(warned({ x: [0, 0.1], values: [0, RT] })); // (flat across region b: it grades none of it)
});
