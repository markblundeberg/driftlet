import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { aqueous, recorder } from '../src/kit.js';
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
