import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { aqueous, bath, build, layer, photogeneration, recorder } from '../src/kit.js';
import { sourceAt, nextBreakpoint } from '../src/device.js';
import { electrodeUnbalanced } from '../src/equation.js';

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

test('a solution read right after a set() of drives alone is the state as solved, until it moves', () => {
  // (read at the new voltage, a patch from 0 to 0.05 V read a plate charge of 1e-2 C/m² and 1e4
  // A/m² that never flowed, against Kirchhoff, with no warning)
  const d = blocking({ V: 0 });
  d.solve();
  const a = d.step(1e-6);
  d.set({ contacts: { right: { V: 0.05 } } });
  const b = d.solution();
  assert.equal(b.plates.right.charge, a.plates.right.charge);
  assert.equal(b.terminals.right.current, a.terminals.right.current);
  const c = d.step(1e-12); // (the step response: ohmic, 0.05 V over the electrolyte, ~749 A/m²)
  assert.ok(Math.abs(c.terminals.right.current / 749 - 1) < 0.01, `${c.terminals.right.current}`);
  assert.equal(c.terminals.right.V, 0.05);
  const e = blocking({ V: 0 });
  e.solve();
  e.set({ contacts: { right: { V: 0.05 } } });
  assert.equal(e.solve().terminals.right.V, 0.05); // (a solve is at the new drives)
});

test('two terminals driven by currents into one stretch: the conservation ledger is theirs together', () => {
  // (the ledger took the contact's current alone: a leak of I·dt/F a step, and check() failed it)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 });
  const d = new Device({
    species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { phi: 'neutral', species: { 'K+': { type: 'conductance', G: 1, offset: 0 }, 'Cl-': 'blocked' }, I: 0 }, right: { V: 0, terminal: 'Cl-', phi: 'neutral', species: { 'Cl-': 'equilibrium', 'K+': 'blocked' } } },
    ports: [{ name: 'p', region: 0, from: 20e-6, to: 40e-6, species: { 'K+': { type: 'conductance', G: 1e6, offset: 0 } }, I: 0 }],
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
  d.solve();
  d.set({ contacts: { left: { I: 1 } }, ports: { p: { I: -1 } } });
  const sol = d.step(0.01);
  assert.ok(sol.converged);
  const K = sol.conservation.find((c) => c.species === 'K+');
  assert.ok(Math.abs(K.drift) < 1e-9, JSON.stringify(K));
});

test("a current ramp's passed charge on a resistive path, from a long first step", () => {
  // (the first step's pair of quadratures agree at its end there, err 0: 1.5 for 1, silently)
  const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const d = new Device({
    species: [{ name: 'e-', z: -1 }], materials: { Cu: { conductor: { species: 'e-', conductivity: 1e6 } } },
    regions: [{ material: 'Cu', length: 1e-6 }],
    contacts: { left: collector(0), right: { ...collector(0), V: undefined, I: { t: [0, 1], values: [0, 2] } } },
    grid: { hmin: 1e-10 },
  });
  d.solve();
  const sol = d.advance(1, { dt0: 1 });
  assert.ok(Math.abs(last(sol.trace.terminals.right.charge) - 1) < 3e-3, `${last(sol.trace.terminals.right.charge)} for 1`);
});

test('a first step cut across calls is seen to stall as in one call: frames step for step as one go', () => {
  // (the stall's count was per call: three calls of one attempt each made no progress where one
  // call of three took the step, and a heavy device's animation sat for frames)
  const mk = () => {
    const d = new Device({ species: [{ name: 'X', z: 0 }],
      materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
      regions: [{ name: 's', material: 'f', length: 100e-6, c0: { X: 1 } }],
      contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
      grid: { hmin: 1e-11, hmax: 1e-6 } });
    d.solve();
    d.set({ contacts: { left: { species: { X: { type: 'equilibrium', mu: 2000 } } } } });
    return d;
  };
  const one = mk().advance(1e-3, { dt0: 1e-6, maxSteps: 12 });
  const d = mk();
  let r = d.advance(1e-3, { dt0: 1e-6, maxSteps: 1 });
  for (let k = 1; k < 12; k++) r = d.advance(1e-3, { maxSteps: 1 });
  assert.ok(one.time > 0);
  assert.equal(r.time, one.time);
});

test("the recorder joins each terminal's charge passed, across its pieces", () => {
  const d = blocking({ V: { t: [0, 1e-9, 1], values: [0, 5e-3, 5e-3] } });
  d.solve();
  const rec = recorder(d, { every: 0.1 });
  const sol = rec.advance(1);
  const q = rec.trace.terminals.right.charge;
  assert.equal(q.length, rec.trace.t.length);
  assert.equal(last(q), last(sol.trace.terminals.right.charge));
});

test('a terminal behind a resistance so small that I·R is below round-off: its trace reads the current that passed', () => {
  // (the trace read 0 where 1e9 A/m² passed, and the charges broke Kirchhoff by 1e7 C/m²)
  const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const d = new Device({
    species: [{ name: 'e-', z: -1 }], materials: { Cu: { conductor: { species: 'e-', conductivity: 1e6 } } },
    regions: [{ material: 'Cu', length: 1e-6 }], contacts: { left: collector(1), right: collector(1) }, grid: { hmin: 1e-10 },
  });
  d.solve();
  d.set({ contacts: { right: { V: 1.001, R: 1e-30 } } });
  const sol = d.advance(0.01);
  const L = last(sol.trace.terminals.left.charge), R = last(sol.trace.terminals.right.charge);
  assert.ok(Math.abs(R / 1e7 - 1) < 1e-6 && Math.abs(L + R) < 1e-6 * R, `left ${L}, right ${R}`);
});

test('a step that throws, or a cold step that fails, leaves the device as it found it', () => {
  // (a throw inside a step left the state at the refused trial, 0.7% off later; a failed cold step
  // left the t = 0⁺ projection behind, coldStart cleared)
  const d = blocking({ V: 0.01 });
  const s = d.solver, real = s.newton, u0 = Float64Array.from(s.u);
  s.newton = function () { this.u[0] += 1; throw new Error('a model of the user\'s'); };
  assert.throws(() => d.step(1e-3), /user/);
  assert.throws(() => d.advance(1e-3), /user/);
  assert.deepEqual(s.u, u0);
  s.newton = function () { return { converged: false, iterations: 0, history: [] }; };
  assert.equal(d.step(1e-3).converged, false);
  assert.deepEqual(s.u, u0);
  assert.equal(s.coldStart, true);
  s.newton = real;
  const twin = blocking({ V: 0.01 });
  assert.equal(d.advance(1e-3).plates.right.charge, twin.advance(1e-3).plates.right.charge);
});

// The cold start (src/start.js): each a start that sat where its mirror or order-permuted twin,
// or the analytic value, didn't, the first step or solve then getting there from it.

test('a bath without a reference reads its cold voltage by φ, whatever order def.species lists its ions in', () => {
  // (its terminal was the first listed: a salt ramped 1 → 9 against baths at 1 read the floating
  // end at +56 or −56 mV by the order, the first step at 0 in both)
  const L = 1e-6, run = (order) => {
    const d = new Device({
      species: order.map((nm) => ({ name: nm, z: nm === 'P' ? 1 : -1 })),
      materials: { w: { epsr: 0, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
      regions: [{ name: 'a', material: 'w', length: L, c0: { P: { x: [0, L], values: [1, 9] }, N: { x: [0, L], values: [1, 9] } } }],
      contacts: { left: { V: 0, bath: { c: { P: 1, N: 1 } } }, right: { I: 0, bath: { c: { P: 1, N: 1 } } } },
      grid: { minCells: 16 },
    });
    return [d.solution().contacts.right.V, d.advance(1e-9).contacts.right.V];
  };
  for (const [cold, first] of [run(['P', 'N']), run(['N', 'P'])]) assert.ok(Math.abs(cold) < 1e-9 && Math.abs(first) < 1e-9, `cold ${cold} V, after 1 ns ${first} V`);
});

test('a bath without a reference outside a cation exchanger is accepted in either order of def.species', () => {
  // (its terminal, the first listed, was N⁻, which the membrane doesn't hold: refused as blocked)
  for (const order of [['P', 'N'], ['N', 'P']]) {
    const d = new Device({
      species: order.map((nm) => ({ name: nm, z: nm === 'P' ? 1 : -1 })),
      materials: { water: { epsr: 0, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } }, mem: { epsr: 0, species: { P: { D: 1e-10, cRef: 1, mu0: 0 } } } },
      regions: [{ name: 'm', material: 'mem', length: 1e-6, fixedCharge: -100 * FARADAY, c0: { P: 100 } }],
      contacts: { left: { V: 0, bath: { c: { P: 1, N: 1 }, material: 'water' } }, right: { V: 0, bath: { c: { P: 1, N: 1 }, material: 'water' } } },
      grid: { minCells: 16 },
    });
    const s = d.solve(), donnan = -d.solver.VT * Math.log(100);
    assert.ok(s.converged && Math.abs(s.phi[0] - donnan) < 1e-9, `${order}: φ ${s.phi[0]} V for ${donnan}`);
  }
});

test('a floating capacitance through a window beside a pinned end starts uncharged', () => {
  // (its starting V averaged φ over the end box too, whose φ the contact sets and whose charge
  // the port doesn't hold: 0.0087 C/m² on two cells, falling with the end box)
  const d = new Device({
    species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }],
    materials: { f: { epsr: 80, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1, c0: { N: { x: [0, 1], values: [1, 9] } } }],
    contacts: { left: { V: 0, terminal: 'P', phi: 'pinned', zeroCharge: 0, species: { P: 'equilibrium', N: 'blocked' } }, right: { phi: 'neutral', species: { P: 'blocked', N: 'blocked' } } },
    ports: [{ name: 'cap', region: 'a', I: 0, area: 1, capacitance: { C: 1 } }],
    grid: { minCells: 2, hmin: 0.5, hmax: 0.5 },
  });
  const Q = d.solution().ports[0].charge;
  assert.ok(Math.abs(Q) < 1e-12, `${Q} C/m²`);
});

test('an electrode at the far end of a ramped salt starts balanced with the salt beside it, as its mirror does', () => {
  // (the balance took the salt at the region's start: a metal on the right of Ag⁺ ramped 10 → 90
  // started 2.2 RT (ln 9) out of equilibrium, the mirror at it; a port's window likewise)
  const T = 298.15, RT = GAS_CONSTANT * T;
  const species = [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }];
  const water = { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } };
  const L = 1e-5, salt = (a, b) => ({ 'Ag+': { x: [0, L], values: [a, b] }, 'NO3-': { x: [0, L], values: [a, b] } });
  const closed = { phi: 'neutral', species: { 'Ag+': 'blocked', 'NO3-': 'blocked', 'e-': 'blocked' } };
  const right = new Device({
    T, species, materials: { water, Ag: { conductor: { species: 'e-', conductivity: 6e7 } } },
    regions: [{ material: 'water', length: L, c0: salt(10, 90) }, { material: 'Ag', length: 1e-6 }],
    interfaces: [{ phi: 'neutral', reactions: [{ right: { 'e-': -1, Ag: 1 }, left: { 'Ag+': -1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 }] }],
    contacts: { left: closed, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } }, grid: { minCells: 8 },
  });
  const g = right.model.grid, s = right.solution(), a = (s.mu['Ag+'][g.regionEnd[0]] + s.mu['e-'][g.regionStart[1]]) / RT;
  assert.ok(Math.abs(a) < 1e-9, `affinity ${a} RT at the metal on the right`);
  // a window at either end: the cold rates the mirror of each other's
  const port = (mirror) => {
    const d = new Device({
      T, species, materials: { water }, regions: [{ name: 'w', material: 'water', length: L, c0: mirror ? salt(90, 10) : salt(10, 90) }],
      contacts: { left: closed, right: closed },
      ports: [{ name: 'Ag', region: 'w', from: mirror ? 0 : 0.9 * L, to: mirror ? 0.1 * L : L, V: 0, terminal: 'e-', area: 1e4, reactions: [{ equation: 'Ag+ + e- = Ag', fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 }] }],
      grid: { minCells: 20, hmin: L / 20, hmax: L / 20 },
    });
    d.solver.assemble(Infinity);
    return Array.from(d.solver.portRates[0][0]);
  };
  const p = port(false), q = port(true).reverse(), scale = Math.max(...p.map(Math.abs));
  p.forEach((r, k) => assert.ok(Math.abs(r - q[k]) < 1e-9 * scale, `rates ${p} against the mirror's ${q}`));
});

test("a floating electrode with a conductance beside its reactions starts where the whole terminal passes what it's set", () => {
  // (the start balanced its reactions alone: 20 A/m² through the link at I: 0; and its search,
  // unchecked, ran out to 1457 V for a drive beyond the reactions' plateau)
  const def = (rx, I, link = true) => ({
    species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }, { name: 'e-', z: -1 }],
    materials: { w: { epsr: 0, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'w', length: 1, c0: { P: 1, N: 1 } }],
    contacts: { left: { V: 0, terminal: 'P', phi: 'neutral', species: { P: 'equilibrium', N: 'blocked' } }, right: { phi: 'neutral', species: { P: 'blocked', N: 'blocked' } } },
    ports: [{ name: 'el', region: 'a', I, terminal: 'e-', area: 1, reactions: [rx], ...(link ? { species: { N: { type: 'conductance', G: 100, offset: 0.2 } } } : {}) }],
    grid: { minCells: 4 },
  });
  const bv = { equation: 'P + e- = S', fixed: { S: 0 }, k0: 1e-6, alpha: 0.5 }, mhc = { equation: 'P + e- = S', fixed: { S: 0 }, k0: 1e-6, lambda: 20000 };
  for (const [rx, I] of [[bv, 0], [mhc, 0], [mhc, 5]]) {
    const d = new Device(def(rx, I)), S = d.solver, V0 = S.termV[2];
    S.assemble(Infinity);
    assert.ok(Math.abs(S.termI[2] - I) < 1e-6, `${rx.lambda ? 'MHC' : 'BV'} at ${I}: cold ${S.termI[2]} A/m² at ${V0} V`);
    if (I === 0) assert.ok(Math.abs(new Device(def(rx, I)).advance(1e-3).ports[0].V - V0) < 1e-6, 'and the first ms stays there');
  }
  // beyond what the reactions alone can pass: left where it was, not at the search's bound
  assert.ok(Math.abs(new Device(def(mhc, 5, false)).solver.termV[2]) < 1, 'MHC beyond its plateau');
});

test("a species held by a port behind a blocked face starts at the port's level, not the right contact's", () => {
  // (its stretch doesn't reach the right contact, which supplied its level all the same: 2.63
  // for the port's 1, kept for ~c/k of the exchange)
  const RT = GAS_CONSTANT * 298.15;
  const d = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6 }, { name: 'b', material: 'f', length: 1e-6 }],
    interfaces: [{ species: { X: 'blocked' } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 2400 } } } },
    ports: [{ name: 'feed', region: 'a', V: 0, species: { X: { type: 'exchange', k: 1e-6, mu: 0 } } }],
    grid: { minCells: 8 },
  });
  const s = d.solution(), last = s.x.length - 1;
  assert.ok(Math.abs(s.c.X[0] - 1) < 1e-12, `cold X in a: ${s.c.X[0]}`);
  assert.ok(Math.abs(s.c.X[last] / Math.exp(2400 / RT) - 1) < 1e-12, `and in b, the contact's: ${s.c.X[last]}`);
  assert.ok(Math.abs(d.advance(1e-3).c.X[0] - 1) < 1e-9);
});

test("an insertion host's wire starts at the open circuit whatever order lists its ion and electrons", () => {
  // (the host's levels were aligned on its first charged species shared across the face: e⁻,
  // listed first and blocked there, left its Li⁺ 0.4 V out of equilibrium and the wire at 0 V)
  const VT = (GAS_CONSTANT * 298.15) / FARADAY;
  const x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95], E = x.map((v) => 0.4 - VT * Math.log(v / (1 - v)));
  const sp = { 'Li+': { name: 'Li+', z: 1 }, 'e-': { name: 'e-', z: -1 } };
  for (const order of [['Li+', 'e-'], ['e-', 'Li+']]) {
    const d = new Device({
      species: order.map((k) => sp[k]),
      materials: {
        mc: { epsr: 0, species: { 'Li+': { D: 1e-9, mu0: 0, cRef: 1000 }, 'e-': { D: 1e-6, mu0: 0, cRef: 1000 } } },
        host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }] },
      },
      regions: [{ name: 'm', material: 'mc', length: 1e-6, c0: { 'Li+': 1000, 'e-': 1000 } }, { name: 'h', material: 'host', length: 1e-6, c0: { 'Li+': 15000 } }],
      interfaces: [{ species: { 'Li+': 'equilibrium', 'e-': 'blocked' } }],
      contacts: {
        left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium', 'e-': { type: 'equilibrium', offset: 0.5 } }, phi: 'bulk' },
        right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      },
      grid: { hmin: 5e-8, hmax: 1e-7 },
    });
    const V = d.solution().terminals.right.V;
    assert.ok(Math.abs(V - 0.4) < 1e-6, `[${order}]: cold ${V} V`);
  }
});

test('a strictly neutral region under a held capacitance starts with ρ + aσ = 0, as the solve has it', () => {
  // (the start left aσ out: a region of P⁺ alone under a gate was refused for the charge the gate
  // balances, though the solve at any εr > 0 gives φ = F c/(a C); and a level-fed one started
  // off it, its first 1 ps step jumping the gate's charge from 0.28 to 0.0035 C/m²)
  const lone = new Device({
    species: [{ name: 'P+', z: 1 }],
    materials: { f: { epsr: 0, species: { 'P+': { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { 'P+': 1 } }],
    contacts: { left: { phi: 'neutral', species: { 'P+': 'blocked' } }, right: { phi: 'neutral', species: { 'P+': 'blocked' } } },
    ports: [{ name: 'g', region: 'a', V: 0, area: 1e7, capacitance: { C: 1 } }],
    grid: { minCells: 8 },
  });
  for (const s of [lone.solution(), lone.solve()]) for (const phi of s.phi) assert.ok(Math.abs(phi - FARADAY / 1e7) < 1e-12, `φ ${phi} V for F/(aC)`);
  const bar = new Device({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * FARADAY }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'dl', region: 'bar', V: 0.2, area: 1e6, capacitance: { C: 0.2 } }],
    grid: { minCells: 40 },
  });
  const Q0 = bar.solution().ports[0].charge, Q1 = bar.step(1e-12).ports[0].charge, Q = bar.solve().ports[0].charge;
  assert.ok(Math.abs(Q1 / Q0 - 1) < 1e-9 && Math.abs(Q / Q0 - 1) < 1e-9, `cold ${Q0}, after 1 ps ${Q1}, steady ${Q} C/m²`);
});

// Definitions accepted that should have been refused, and refused that should have been
// accepted: Astra's third pass over the definition's normalization.

test('a contact at an end of no cross-section (a sphere centre, a cylinder axis) holds no level or φ there', () => {
  // (an equilibrium link at a sphere's centre passed 1.08 mol/s on 4 cells, 0.016 on 256, where
  // nothing passes a point; a pinned φ there stored 0.4 of a capacitor of radius h/2)
  const def = (geometry, left) => ({ geometry, species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1, c0: { X: 1 } }],
    contacts: { left, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } } },
    grid: { minCells: 4, hmin: 0.25, hmax: 0.25 } });
  const held = { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 2400 } } };
  for (const type of ['spherical', 'cylindrical']) {
    assert.throws(() => new Device(def({ type, r0: 0 }, held)), /contacts\.left: the cross-section is zero at this end.*holds 'X's level/);
    assert.throws(() => new Device(def({ type, r0: 0 }, { phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } })), /holds φ there \('pinned'\)/);
  }
  // At a radius, it holds; at the centre a rate link passes what the end's zero area makes of it, nothing.
  assert.ok(new Device(def({ type: 'spherical', r0: 0.1 }, held)).solve().converged);
  const ex = new Device(def({ type: 'spherical', r0: 0 }, { phi: 'neutral', species: { X: { type: 'exchange', k: 1, mu: 2400 } } })).solve();
  assert.equal(ex.contacts.left.flux.X, 0);
});

test("a velocity link on a conductor's carrier is refused, as a permeability is", () => {
  // (a conductor's c is its carrier's excess, 0 in the bulk, so v·c·(…) passed nothing at any v,
  // 1e30 m/s included: current 0 where an equilibrium link passes σV/L = 1 A/m²)
  const def = (link) => ({ species: [{ name: 'e-', z: -1 }], materials: { Cu: { conductor: { species: 'e-', conductivity: 1 } } },
    regions: [{ material: 'Cu', length: 1 }],
    contacts: { left: { V: 1, terminal: 'e-', species: { 'e-': link }, phi: 'neutral' }, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' } } });
  assert.throws(() => new Device(def({ type: 'velocity', v: 1e6 })), (e) => e instanceof DeviceError && /contacts\.left\.species\.e-: a velocity link .* 'Cu' is a conductor/.test(e.message));
  assert.ok(Math.abs(Math.abs(new Device(def('equilibrium')).solve().current) - 1) < 1e-12);
});

test("a bath on a conductor is refused at a port, and called a conductor's at a contact", () => {
  // (a port's bath read the conductor's placeholder μ° = 0 and cRef = 1: e⁻ at 100 put the wire
  // VT·ln 100 = 0.118 V off the port's V, 0.118 A/m² passing between two grounded terminals)
  const def = (port, left) => ({ species: [{ name: 'e-', z: -1 }], materials: { Au: { conductor: { species: 'e-', conductivity: 1 } } },
    regions: [{ name: 'w', material: 'Au', length: 1 }],
    contacts: { left: left ?? { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [port] });
  const port = { name: 'p', region: 'w', V: 0, species: { 'e-': 'equilibrium' } };
  assert.throws(() => new Device(def({ ...port, bath: { c: { 'e-': 100 } } })), /ports\[0\]\.bath: w \(Au\) is a conductor/);
  assert.throws(() => new Device(def({ ...port, terminal: 'e-' }, { bath: { c: { 'e-': 100 } } })), (e) => /contacts\.left\.bath: Au at this end is a conductor/.test(e.message) && !/insertion host/.test(e.message));
  assert.equal(new Device(def({ ...port, terminal: 'e-' })).solve().ports[0].current, 0);
});

test("a reacting port's electrode side is its carrier and fixed participants, whatever else is declared", () => {
  // (a Y in fixed that another material declared went on the electrode's side as the carrier:
  // 'P+ = Y' passed 2.2 mA/m² with no e⁻ and its 9000 J/mol dropped; 'P+ + e- = Y' failed on a NaN)
  const def = (declared, equation) => ({
    species: [{ name: 'P+', z: 1 }, { name: 'N-', z: -1 }, { name: 'e-', z: -1 }, ...(declared ? [{ name: 'Y', z: declared.z }] : [])],
    materials: {
      w: { epsr: 0, species: { 'P+': { D: 1e-9, cRef: 1, mu0: 0 }, 'N-': { D: 1e-9, cRef: 1, mu0: 0 } } },
      ...(declared ? { o: { epsr: 0, species: { Y: { D: 1e-9, cRef: 1, mu0: 0 } } } } : {}),
    },
    regions: [{ name: 'a', material: 'w', length: 1e-5, c0: { 'P+': 1, 'N-': 1 } }],
    contacts: { left: { V: 0, terminal: 'P+', species: { 'P+': 'equilibrium', 'N-': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' }, right: { phi: 'neutral' } },
    ports: [{ name: 'el', region: 'a', V: 0, terminal: 'e-', area: 1, reactions: [{ equation, fixed: { Y: 9000 }, k0: 1e-3, alpha: 0.5 }] }],
    grid: { minCells: 8 },
  });
  assert.throws(() => new Device(def({ z: 1 }, 'P+ = Y')), /ports\[0\]\.reactions\[0\]: no 'e-' from the electrode takes part/);
  const twin = new Device(def(null, 'P+ + e- = Y')), d = new Device(def({ z: 0 }, 'P+ + e- = Y'));
  assert.equal(d.model.ports[0].reactions[0].fixedA, twin.model.ports[0].reactions[0].fixedA);
  assert.ok(Math.abs(d.model.ports[0].reactions[0].fixedA + 9000 / d.model.RT) < 1e-12);
  assert.deepEqual(d.model.ports[0].reactions[0].part, twin.model.ports[0].reactions[0].part);
  const s = d.solve();
  assert.ok(s.converged && s.ports[0].current !== 0);
  assert.equal(s.ports[0].current, twin.solve().ports[0].current);
});

test("Marcus–Hush–Chidsey's one electron counts charge: a divalent carrier with ν = 1 is refused", () => {
  // (B²⁻ with ν = −1 passed as one electron, at a face and at a port: two charges through one
  // electron's Fermi-level integral)
  const lambda = units.eV(0.5);
  const def = (zc) => {
    const car = zc === -2 ? 'B2-' : 'e-', ox = zc === -2 ? 'O2+' : 'O+';
    return {
      species: [{ name: ox, z: -zc, cRef: 1000 }, { name: 'R', z: 0, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: car, z: zc }],
      materials: {
        water: { epsr: 0, species: { [ox]: { D: 1e-9, mu0: 10e3 }, R: { D: 0.8e-9, mu0: -86.485e3 }, 'K+': { D: 1.96e-9, mu0: -283.3e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } },
        Pt: { conductor: { species: car, conductivity: 9.4e6 } },
      },
      regions: [{ material: 'Pt', length: 1e-6 }, { material: 'water', length: 20e-6 }],
      interfaces: [{ phi: 'neutral', reactions: [{ left: { [car]: -1 }, right: { [ox]: -1, R: 1 }, k0: 1e-5, lambda }] }],
      contacts: { left: { V: 0, terminal: car, species: { [car]: 'equilibrium' }, phi: 'bulk' }, right: { bath: { c: { [ox]: 1, R: 1, 'K+': 1000, 'Cl-': 1000 - zc }, reference: 'Cl-' } } },
      grid: { minCells: 8 },
    };
  };
  const onPort = (d, zc) => {
    const car = zc === -2 ? 'B2-' : 'e-', ox = zc === -2 ? 'O2+' : 'O+';
    d.interfaces = [{ phi: 'neutral' }];
    d.ports = [{ name: 'wire', region: 1, from: 0, to: 1e-6, V: 0, terminal: car, area: 1e4, reactions: [{ equation: `${ox} + ${car} = R`, k0: 1e-5, lambda }] }];
    return d;
  };
  assert.ok(new Device(def(-1)).model.interfaces[0].reactions[0].mhc);
  assert.ok(new Device(onPort(def(-1), -1)).model.ports[0].reactions[0].mhc);
  assert.throws(() => new Device(def(-2)), /interfaces\[0\].*Marcus–Hush–Chidsey kinetics is one electron to or from a metal, its carrier with \|ν z\| = 1/);
  assert.throws(() => new Device(onPort(def(-2), -2)), /ports\[0\]\.reactions\[0\]: Marcus–Hush–Chidsey kinetics is one electron to or from the electrode, its carrier with \|ν z\| = 1/);
});

test('a port bath may list one of Fermi–Dirac bands declared together: they are each on its own', () => {
  // (refused as coupled, "each one's level depends on the others'", and accepted declared apart)
  const def = (statistics) => ({ species: [{ name: 'A', z: 0 }, { name: 'B', z: 0 }],
    materials: { f: { epsr: 4, species: { A: { D: 1e-9, cRef: 1, mu0: 0 }, B: { D: 1e-9, cRef: 1, mu0: 0 } }, statistics } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { A: 1, B: 1 } }],
    contacts: { left: { phi: 'neutral', species: { A: 'blocked', B: 'blocked' } }, right: { phi: 'neutral', species: { A: 'blocked', B: 'blocked' } } },
    ports: [{ name: 'p', region: 'a', from: 0, to: 0.5e-6, bath: { c: { A: 3 } }, species: { A: { type: 'exchange', k: 1 } } }],
    grid: { minCells: 8 } });
  const together = new Device(def([{ type: 'fermi-dirac', species: ['A', 'B'] }]));
  const apart = new Device(def([{ type: 'fermi-dirac', species: ['A'] }, { type: 'fermi-dirac', species: ['B'] }]));
  assert.equal(together.model.ports[0].species[0].mu, apart.model.ports[0].species[0].mu);
  // (5251.40 J/mol: ζ = 2.11839 for 𝓕½(ζ) = 3, by an independent quadrature, to its 3e-3 J/mol)
  assert.ok(Math.abs(together.model.ports[0].species[0].mu - 5251.403) < 0.01);
  assert.ok(Math.abs(together.solve().c.A[2] - 3) < 1e-6);
  // A lattice still couples them: listing one is refused.
  assert.throws(() => new Device(def([{ type: 'lattice', species: ['A', 'B'], cMax: 10 }])), /share lattice statistics .* list them all/);
});

test('a species belongs to one statistics model, ideal among them, whatever the order', () => {
  // ([ideal, lattice] was accepted as the lattice, and [lattice, ideal] refused)
  const def = (statistics) => ({ species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } }, statistics } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
    grid: { minCells: 8 } });
  const ideal = { type: 'ideal', species: ['X'] }, lattice = { type: 'lattice', species: ['X'], cMax: 2 };
  for (const st of [[ideal, lattice], [lattice, ideal], [ideal, ideal]]) {
    assert.throws(() => new Device(def(st)), /statistics\[1\]\.species\[0\]: 'X' already belongs to another statistics model/, JSON.stringify(st.map((m) => m.type)));
  }
  assert.equal(new Device(def([ideal])).solution().mu.X[0], 0);
  assert.ok(Math.abs(new Device(def([lattice])).solution().mu.X[0] - GAS_CONSTANT * 298.15 * Math.log(2)) < 1e-9);
});

test("an OCV table's insertion species take no μ° grade: the table sets their level", () => {
  // (graded to a constant 10 kJ/mol, the table read 3.86766 V where it says 3.97130, with the
  // material's μ° at 10 kJ/mol as at 0: Δμ°/F = 104 mV off a staircase of such layers)
  const x = [0.1, 0.5, 0.9], E = [4.2, 3.9, 3.6], L = 1e-6;
  const def = (mu0, grade, ocv = { x, E, muRef: 0 }) => ({ species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: { h: { epsr: 0, species: { 'Li+': { D: 1e-14, cRef: 1, mu0 }, 'e-': { D: 1e-8, cRef: 1, mu0: 0 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 1000, ...(ocv ? { ocv } : { A: [] }) }] } },
    regions: [{ name: 'a', material: 'h', length: L, c0: { 'Li+': 400, 'e-': 400 }, ...(grade ? { grade } : {}) }],
    contacts: { left: { V: 0, terminal: 'e-', phi: 'neutral', species: { 'Li+': 'blocked', 'e-': 'equilibrium' } }, right: { phi: 'neutral', species: { 'Li+': 'blocked', 'e-': 'blocked' } } },
    grid: { minCells: 8 } });
  const comb = (d) => { const s = new Device(d).solve(); return s.mu['Li+'][3] + s.mu['e-'][3]; };
  assert.ok(Math.abs(comb(def(0)) - comb(def(10e3))) < 1e-6);
  for (const sp of ['Li+', 'e-']) {
    assert.throws(() => new Device(def(0, { [sp]: { mu0: { x: [0, L], values: [10e3, 10e3] } } })), new RegExp(`regions\\[0\\]\\.grade\\.${sp.replace('+', '\\+')}\\.mu0: .* OCV table`));
  }
  // Its D may be graded, and a Redlich–Kister host's μ° (no table: μ° is its level's).
  assert.ok(new Device(def(0, { 'Li+': { D: { x: [0, L], values: [1e-14, 2e-14] } } })).solve().converged);
  assert.ok(new Device(def(0, { 'Li+': { mu0: { x: [0, L], values: [10e3, 10e3] } } }, null)).solve().converged);
});

test("a port's bath reads a graded μ° in its window, refused only where it varies there", () => {
  // (refused for any grade object at all: one knot, a constant, one beyond the window, and of a
  // species whose link gave its own level)
  const VT = (GAS_CONSTANT * 298.15) / FARADAY, sp = { D: 1e-9, cRef: 1, mu0: 0 };
  const def = (grade, link = 'equilibrium', statistics) => ({ species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }],
    materials: { m: { epsr: 0, species: { P: sp, N: sp }, ...(statistics ? { statistics } : {}) } },
    regions: [{ name: 'a', material: 'm', length: 1, c0: { P: 1, N: 1 }, ...(grade ? { grade } : {}) }],
    contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'p', region: 'a', from: 0.1, to: 0.2, V: 0, species: { P: link }, bath: { c: { P: 2, N: 2 } } }],
    grid: { minCells: 10 } });
  const offset = (d) => new Device(d).model.ports[0].species[0].offset;
  assert.ok(Math.abs(offset(def(null)) - VT * Math.log(2)) < 1e-15);
  assert.equal(offset(def({ P: { mu0: { x: [0.5], values: [0] } } })), offset(def(null)));
  assert.equal(offset(def({ P: { mu0: { x: [0, 0.8, 1], values: [0, 0, 1000] } } })), offset(def(null)));
  // a constant grade is the level's μ° there, not the material's
  assert.ok(Math.abs(offset(def({ P: { mu0: { x: [0, 1], values: [1000, 1000] } } })) - (VT * Math.log(2) + 1000 / FARADAY)) < 1e-15);
  // varying across the window, it is refused for the level the bath would give, not for a link's own
  const varying = { P: { mu0: { x: [0, 1], values: [0, 1000] } } };
  assert.throws(() => new Device(def(varying)), /ports\[0\]\.species\.P: a grades 'P's mu0 or cRef, so a composition outside gives no one level across the window/);
  assert.equal(offset(def(varying, { type: 'equilibrium', offset: 0.01 })), 0.01);
  // under a lattice, listing both (as it must) with P's level by its link
  const lattice = [{ type: 'lattice', species: ['P', 'N'], cMax: 10 }];
  assert.equal(offset(def(varying, { type: 'equilibrium', offset: 0.01 }, lattice)), 0.01);
});

test("a bath's reference electrode counts each participant once, its ν summed over repeated terms", () => {
  // (M⁺ + e⁻ + e⁻ = M(s) passed, read by its first e⁻; M⁺ + M⁺ + 2 e⁻ = M₂(s) was refused, and
  // the kit's bath() refused both, as needing electrons on one side and having 'M+' and 'M+')
  const species = [{ name: 'M+', z: 1 }, { name: 'A-', z: -1 }];
  const flagged = (eq) => electrodeUnbalanced({ electrode: eq, reference: 'M+' }, species);
  assert.equal(flagged('M+ + M+ + 2 e- = M2(s)'), null);
  assert.equal(flagged('2 M+ + 2 e- = M2(s)'), null);
  assert.match(flagged('M+ + e- + e- = M(s)'), /isn't charge-balanced \(1 M\+ of charge 1 to 2 e-\)/);
  const c = { 'M+': 1, 'A-': 1 }, fixed = { 'M2(s)': -1000 };
  const repeated = bath(c, { equation: 'M+ + M+ + 2 e- = M2(s)', fixed }), once = bath(c, { equation: '2 M+ + 2 e- = M2(s)', fixed });
  assert.equal(repeated.bath.offset, once.bath.offset);
  assert.ok(Math.abs(once.bath.offset + 1000 / (2 * FARADAY)) < 1e-15); // μ_M2/(2F)
  assert.deepEqual({ ...repeated.bath, electrode: null }, { ...once.bath, electrode: null });
  assert.ok(bath(c, { equation: 'M+ + e- + e- = M(s)', fixed: { 'M(s)': 0 } }).bath.electrode);
});

test('SRH and Auger take no fixed participants: they consume one carrier of each sign and make nothing', () => {
  // ('e- + h+ = S' with S fixed at 5 kJ/mol was accepted, in the bulk and at a face, its n·p at
  // e^(5000/RT) = 7.52 n_i²; mass action takes it, as before)
  const semi = { s: { epsr: 0, species: { 'e-': { D: 1e-3, cRef: 1, mu0: 0 }, 'h+': { D: 1e-3, cRef: 1, mu0: 0 } } } };
  const region = (name) => ({ name, material: 's', length: 1e-6, fixedCharge: 100 * FARADAY, c0: { 'e-': 100.01, 'h+': 0.01 } });
  const contacts = { left: { V: 0, terminal: 'e-', phi: 'bulk', species: { 'e-': 'equilibrium', 'h+': 'blocked' } }, right: { phi: 'neutral', species: { 'e-': 'blocked', 'h+': 'blocked' } } };
  const bulk = (rx) => ({ species: [{ name: 'e-', z: -1 }, { name: 'h+', z: 1 }], materials: semi, regions: [region('a')], bulkReactions: [rx], contacts, grid: { minCells: 8 } });
  const face = (rx) => ({ ...bulk(rx), bulkReactions: undefined, regions: [region('a'), region('b')], interfaces: [{ phi: 'neutral', species: { 'e-': 'equilibrium', 'h+': 'blocked' }, reactions: [rx] }] });
  const srh = { s: { tauN: 1e-6, tauP: 1e-6 } }, auger = { s: { Cn: 1e-3, Cp: 1e-3 } };
  for (const [rx, message] of [
    [{ equation: 'e- + h+ = S', fixed: { S: 5000 }, srh }, /bulkReactions\[0\]\.srh\.s: SRH kinetics .*no fixed participants/],
    [{ equation: 'e- + h+ + S = 0', fixed: { S: 0 }, srh }, /bulkReactions\[0\]\.srh\.s: SRH kinetics/],
    [{ equation: 'e- + h+ = S', fixed: { S: 5000 }, auger }, /bulkReactions\[0\]\.auger\.s: Auger kinetics .*no fixed participants/],
  ]) assert.throws(() => new Device(bulk(rx)), message);
  assert.throws(() => new Device(face({ equation: 'e-(left) + h+(right) = S', fixed: { S: 5000 }, srh: { vn: 1, vp: 1 } })), /interfaces\[0\].*reactions\[0\]\.srh: SRH kinetics .*no fixed participants/);
  assert.ok(new Device(face({ equation: 'e-(left) + h+(right) = 0', srh: { vn: 1, vp: 1 } })).solve().converged);
  for (const rx of [{ equation: 'e- + h+ = 0', srh }, { equation: 'e- + h+ = 0', auger }]) assert.ok(new Device(bulk(rx)).solve().converged);
  const s = new Device(bulk({ equation: 'e- + h+ = S', fixed: { S: 5000 }, kf: { s: 1 } })).solve();
  assert.ok(Math.abs((s.c['e-'][4] * s.c['h+'][4]) / Math.exp(5000 / (GAS_CONSTANT * 298.15)) - 1) < 1e-6);
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
  assert.ok(new Device(supply({ V: 0 }))); // (the default, as 1.17 took it)
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
  refused(cu({ fixedCharge: 1 }), /regions\[0\]\.fixedCharge: a conductor region is neutral in bulk/);
  for (const zero of [{ fixedCharge: 0 }, { velocity: 0 }, { mixing: 0 }]) assert.ok(new Device(cu(zero)), JSON.stringify(zero)); // (the defaults, as 1.17 took them)
  refused(cu({ fixedCharge: { x: [0, 1e-6], values: [0, 0] } }), /regions\[0\]\.fixedCharge: a conductor region/);
  refused(cu({ velocity: 1e-3 }), /regions\[0\]\.velocity: a conductor region has no flow/);
  refused(cu({ mixing: 1e-9 }), /regions\[0\]\.mixing: a conductor region has no flow/);
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
