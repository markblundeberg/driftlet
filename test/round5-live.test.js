import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Worker } from 'node:worker_threads';
import { Device, DeviceError } from '../src/index.js';
import { build, layer, metal, ohmic, live, recorder } from '../src/kit.js';
import { Session, between } from '../src/live.js';
import { WHOLE } from '../src/merge.js';

// Round 5's live() and recorder() faults (Astra's fifth pass), each with the case that found it.

const wire = (right, left = ohmic(0, ['e-'])) => build({ library: [metal('Cu')], stack: [left, layer('Cu', 100e-6), right] });
// A neutral species in a film, held at a level on the left: cheap, and slow enough to step.
const film = (cells = 8) => ({
  species: [{ name: 'X', z: 0 }],
  materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
  contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 2000 } } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
  grid: { minCells: cells },
});
// A Session whose next n warm solves are made to fail, so its ramp runs.
const failNext = (s, n = 1) => {
  const solve = s.device.solve.bind(s.device);
  s.device.solve = (...a) => (n-- > 0 ? { ...solve(...a), converged: false } : solve(...a));
};
const nodeWorker = () => {
  const w = new Worker(new URL('./fixtures/web-worker.js', import.meta.url));
  return {
    set onmessage(fn) {
      w.on('message', (data) => fn({ data }));
    },
    postMessage: (m) => w.postMessage(m),
    terminate: () => w.terminate(),
  };
};
const race = (p, ms = 2000) => Promise.race([p.then(() => 'settled'), new Promise((r) => setTimeout(() => r('pending'), ms))]);

test('live ramp: a kit contact stays whole along the way, and the last step is the target itself (A2, A24)', () => {
  // A contact behind a series R replaced by the kit's plain one: ramped to after a failed warm
  // solve, the ramp's plain copies merged onto the old contact, the R survived, and the device
  // read 1.7e-8 V where 0.01 is right (the bad definition then kept as the good one).
  const target = wire(ohmic(0.01, ['e-']));
  const s = new Session(wire(ohmic({ V: 0.05, R: 1e-6 }, ['e-'])));
  assert.ok(!s.update([]).info.failed);
  failNext(s);
  const r = s.update([target]);
  const direct = new Device(target).solve();
  assert.ok(!r.info.failed && r.info.ramp === 2, JSON.stringify(r.info));
  assert.equal(r.solution.terminals.right.V, 0.01);
  assert.ok(Math.abs(r.solution.current / direct.current - 1) < 1e-9, `${r.solution.current} vs ${direct.current}`);
  assert.equal(s.device.def.contacts.right.R, undefined);
  // between() keeps the mark, and lands on b at s = 1 (0.05 → 1e-18 rounded to 0).
  assert.equal(between(target.contacts.right, target.contacts.right, 0.5)[WHOLE], true);
  assert.equal(between(0.05, 1e-18, 1), 1e-18);
  assert.equal(between(1e-3, 1e-21, 1), 1e-21);
  assert.equal(between(0.05, 1e-18, 0), 0.05);
  // A drive ramped to a value that the old interpolation rounded to 0: the device holds the target.
  const t = new Session(film());
  t.update([{}]);
  failNext(t);
  const u = t.update([{ contacts: { left: { species: { X: { mu: 1e-18 } } } } }]);
  assert.ok(!u.info.failed && u.info.ramp === 2, JSON.stringify(u.info));
  assert.equal(t.device.def.contacts.left.species.X.mu, 1e-18);
});

test('live set(): the patch is copied when given, functions kept (A17)', async () => {
  const dev = live(wire(ohmic(0, ['e-'])));
  await dev.ready;
  const patch = { contacts: { right: { V: 0.01 } } };
  const p = dev.set(patch);
  patch.contacts.right.V = 0.02; // (edited while it waited: was solved at 0.02)
  assert.equal((await p).solution.terminals.right.V, 0.01);
});

test('live: when nothing has converged, ready and onsolution get no solution, not the unconverged try (A33a)', async () => {
  const real = Device.prototype.solve;
  Device.prototype.solve = function (...a) {
    return { ...real.apply(this, a), converged: false };
  };
  const seen = [];
  try {
    const dev = live(film(), { onsolution: (sol, info) => seen.push([sol, info]) });
    const r = await dev.ready;
    await new Promise((res) => setTimeout(res, 10));
    assert.equal(r.solution, null);
    assert.ok(r.info.failed);
    assert.equal(dev.solution, null);
    assert.deepEqual(seen, []);
  } finally {
    Device.prototype.solve = real;
  }
});

test('live in a worker: a change that cannot be posted is refused alone, and the others made with it are solved (A33b)', async () => {
  const dev = live(wire(ohmic(0, ['e-'])), { worker: nodeWorker });
  try {
    await dev.ready;
    const busy = dev.set({ contacts: { right: { V: 0.004 } } }); // (posted; the next two wait and go together)
    const valid = dev.set({ contacts: { right: { V: 0.01 } } });
    const fn = dev.set({ materials: { Cu: { statistics: [{ type: 'custom', species: ['e-'], f: (x) => x }] } } });
    await busy;
    const [a, b] = await Promise.all([valid, fn]);
    assert.ok(!a.info.failed, JSON.stringify(a.info));
    assert.equal(a.solution.terminals.right.V, 0.01); // (was null, the change lost for good)
    assert.ok(b.info.failed && /can't be sent to a worker/.test(b.info.error), JSON.stringify(b.info));
    assert.equal(b.solution, a.solution);
    assert.equal((await dev.set({})).solution.terminals.right.V, 0.01); // (was 0.004)
  } finally {
    await dev.close();
  }
});

test('live in a worker: close() settles every promise owed, and a set() after it (A50)', async () => {
  let terminated = 0;
  const stub = () => ({ postMessage() {}, terminate() { terminated++; }, onmessage: null, onerror: null }); // (busy: never replies)
  const dev = live(film(), { worker: stub });
  const later = dev.set({ T: 300 });
  dev.close();
  assert.equal(terminated, 1);
  assert.equal(await race(dev.ready), 'settled');
  assert.equal(await race(later), 'settled');
  const after = dev.set({ T: 310 });
  assert.equal(await race(after), 'settled');
  const r = await after;
  assert.ok(r.info.failed && /closed/.test(r.info.error), JSON.stringify(r.info));
  assert.equal(r.solution, null);
});

test('recorder: a frame due at the call end within round-off is kept there (A12)', () => {
  // every 0.1 from 0: the third frame is due at 3·0.1 = 0.30000000000000004, which advance(0.3)
  // stopped short of, keeping no frame; frame(0.3) gave the state at 0.2.
  const d = new Device(wire(ohmic({ I: 1 }, ['e-'])));
  d.solve();
  const rec = recorder(d, { every: 0.1 });
  rec.advance(0.3);
  assert.deepEqual(rec.frames.map((f) => f.time), [0, 0.1, 0.2, 0.3]);
  assert.equal(rec.frame(0.3).time, 0.3);
  rec.advance(0.7);
  assert.equal(rec.frames.length, 8);
  assert.equal(rec.frame(0.6), rec.frames[6]); // (at 0.6000000000000001)
  assert.equal(rec.frame(0.7).time, 0.7);
});

test("recorder: probes given to advance() are refused, before a step (A27)", () => {
  const d = new Device(film());
  d.solve();
  const rec = recorder(d, { every: 1e-3, probes: [{ x: 0.5e-6, species: 'X' }] });
  const t = d.solver.time;
  assert.throws(() => rec.advance(2.5e-3, { probes: [{ x: 0.2e-6, species: 'X' }] }), (e) => e instanceof DeviceError && /give.*recorder\(|given to recorder\(/.test(e.message));
  assert.equal(d.solver.time, t);
  assert.equal(rec.trace.t.length, 1);
});

test('recorder: a change in the terminals is refused before the device advances (A34)', () => {
  const bar = (ports) => ({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * 96485.33212 }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports,
    grid: { hmax: 0.5e-6 },
  });
  const a = { name: 'a', region: 'bar', from: 1e-6, to: 3e-6, V: 0.05, area: 1e6, capacitance: { C: 0.2 } };
  const b = { name: 'b', region: 'bar', from: 6e-6, to: 8e-6, V: 0.02, area: 1e6, capacitance: { C: 0.05 } };
  for (const [from, to] of [[[a], [a, b]], [[a, b], [a]]]) {
    const d = new Device(bar(from));
    d.solve();
    d.set({ ports: { a: { V: 0.1 } } });
    const rec = recorder(d);
    rec.advance(1e-3);
    const n = rec.trace.t.length;
    d.set({ ports: to });
    const t = d.solver.time;
    assert.throws(() => rec.advance(2e-3), (e) => e instanceof DeviceError && /terminals changed/.test(e.message) && /new recorder/.test(e.message));
    assert.equal(d.solver.time, t);
    assert.equal(rec.trace.t.length, n);
    for (const nm in rec.trace.terminals) assert.equal(rec.trace.terminals[nm].V.length, n, nm);
  }
});

test("recorder: a call's maxSteps and dt0 are the call's, not each piece's (A37)", () => {
  // maxSteps 1 over a call cut at frames: as one advance(), it stops after its one step.
  const one = new Device(film());
  one.solve();
  const plain = one.advance(2.2e-3, { maxSteps: 1, dt0: 1e-3 });
  const d = new Device(film());
  d.solve();
  const rec = recorder(d, { every: 1e-3 });
  const s = rec.advance(2.2e-3, { maxSteps: 1, dt0: 1e-3 });
  assert.equal(s.time, plain.time); // (was 2.2e-3, done: a step for each of three pieces)
  assert.equal(s.done, false);
  assert.equal(s.steps, 1);
  assert.ok(s.warnings.some((w) => /maxSteps \(1\) used up/.test(w)), JSON.stringify(s.warnings));
  // Run out inside a later piece: the call's limit and end said, once.
  const e = new Device(film());
  e.solve();
  const rec5 = recorder(e, { every: 1e-5 }), r5 = rec5.advance(5e-3, { maxSteps: 5, dt0: 1e-6 });
  assert.ok(!r5.done && r5.steps + r5.rejected === 5 && rec5.frames.length === 2, `${r5.steps} + ${r5.rejected}, ${rec5.frames.length} frames`);
  const told = r5.warnings.filter((w) => /maxSteps/.test(w));
  assert.equal(told.length, 1, JSON.stringify(told));
  assert.match(told[0], /of 0\.00500000: maxSteps \(5\) used up/);
  // dt0 the call's first step only: the run is the one with dt0 given in a call of its own
  // first (was 43 steps where 18 did, dt0 restarting after every frame).
  const steps = (calls) => {
    const dev = new Device(film(40));
    dev.solve();
    const r = recorder(dev, { every: 1e-3 });
    let n = 0;
    for (const [t, o] of calls) n += r.advance(t, o).steps;
    return { n, t: r.trace.t };
  };
  const once = steps([[3e-3, { dt0: 1e-7 }]]), twin = steps([[1e-7, { dt0: 1e-7 }], [3e-3]]);
  assert.equal(once.n, twin.n); // (17 each)
  assert.deepEqual(once.t, twin.t);
});

test('recorder: frames due before the device was advanced on its own are skipped, not kept as copies (A38)', () => {
  const d = new Device(film());
  d.solve();
  const rec = recorder(d, { every: 1 });
  d.advance(2.5);
  rec.advance(3);
  assert.deepEqual(rec.frames.map((f) => f.time), [0, 3]); // (was 0, 2.5, 2.5, 3)
});

test('a concentration probe in a metal reads NaN, as the solution does there (A40)', () => {
  const dev = new Device(wire(ohmic(0.01, ['e-'])));
  dev.solve();
  const rec = recorder(dev, { every: 1e-9, probes: [{ x: 50e-6, species: 'e-' }, { x: 50e-6, species: 'e-', quantity: 'V' }] });
  rec.advance(2e-9);
  assert.ok(rec.trace.probes[0].every(Number.isNaN), String(rec.trace.probes[0].slice(0, 3))); // (was 0)
  assert.ok(Math.abs(rec.trace.probes[1][0] - 0.005) < 1e-12); // (the level still reads)
});
