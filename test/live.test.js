import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Worker } from 'node:worker_threads';
import { Device, units } from '../src/index.js';
import { build, layer, ohmic, live } from '../src/kit.js';
import { Session } from '../src/live.js';

// The live wrapper: changes merged while a solve runs, warm starts, ramps from the last good
// state, the last good solution kept on failure, and the same results from a worker.

const Si = {
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: units.perCm3(2.8e19) }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: units.perCm3(1.04e19) } } } },
};
const pn = (Nd, Na, V) =>
  build({
    library: [Si],
    stack: [ohmic(0), layer('Si', 1e-6, { donors: units.perCm3(Nd) }), layer('Si', 1e-6, { acceptors: units.perCm3(Na) }), ohmic(V)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e13 } }],
    grid: { hmin: 1e-9, hmax: 20e-9 },
  });

test('changes made during a solve are merged and solved once, and match a direct solve', async () => {
  const seen = [];
  const dev = live(pn(1e17, 1e16, 0), { onsolution: (sol) => seen.push(sol.terminals.right.V) });
  const promises = [];
  for (let k = 1; k <= 10; k++) promises.push(dev.set({ contacts: { right: { V: 0.04 * k } } }));
  const results = await Promise.all(promises);
  // The first solve (the definition as given) was under way; the ten changes made one more.
  assert.deepEqual(seen, [0, 0.4]);
  for (const r of results) assert.equal(r.solution, results[9].solution);
  const direct = new Device(pn(1e17, 1e16, 0.4)).solve();
  assert.ok(Math.abs(results[9].solution.current / direct.current - 1) < 1e-9);
  assert.equal(dev.solution, results[9].solution);
});

test('a jump that fails from the warm start is ramped to; an invalid change keeps the last good solution', async () => {
  // From light doping at equilibrium to heavy doping at bias. The solver's own continuations
  // now make that jump, so its first solve is made to fail, as a harder one would.
  const session = new Session(pn(1e15, 1e15, 0));
  session.update([]);
  const solve = session.device.solve.bind(session.device);
  let failures = 1;
  session.device.solve = (opts) => {
    const sol = solve(opts);
    return failures-- > 0 ? { ...sol, converged: false } : sol;
  };
  const jump = session.update([pn(1e19, 1e19, 0.6)]);
  assert.ok(!jump.info.failed && jump.info.ramp >= 2, JSON.stringify(jump.info));
  const direct = new Device(pn(1e19, 1e19, 0.6)).solve();
  assert.ok(Math.abs(jump.solution.current / direct.current - 1) < 1e-8);

  const dev = live(pn(1e15, 1e15, 0));
  await dev.ready;
  const jumped = await dev.set(pn(1e19, 1e19, 0.6));

  const bad = await dev.set({ contacts: { right: { V: 'high' } } });
  assert.ok(bad.info.failed && /contacts\.right\.V/.test(bad.info.error));
  assert.equal(bad.solution, jumped.solution);
  const next = await dev.set({ contacts: { right: { V: 0.5 } } });
  assert.ok(!next.info.failed && next.info.ramp === 0 && next.solution.converged);
});

test('an invalid definition reports its error instead of a solution', async () => {
  const dev = live({ species: [] });
  const r = await dev.ready;
  assert.equal(r.solution, null);
  assert.ok(r.info.failed && /species/.test(r.info.error));
});

test('in a worker, the same solutions', async () => {
  const url = new URL('./fixtures/web-worker.js', import.meta.url);
  const nodeWorker = () => {
    const w = new Worker(url);
    return {
      set onmessage(fn) {
        w.on('message', (data) => fn({ data }));
      },
      postMessage: (m) => w.postMessage(m),
      terminate: () => w.terminate(),
    };
  };
  const there = live(pn(1e17, 1e16, 0), { worker: nodeWorker });
  const here = live(pn(1e17, 1e16, 0));
  for (const V of [0.2, 0.5]) {
    const [a, b] = await Promise.all([there.set({ contacts: { right: { V } } }), here.set({ contacts: { right: { V } } })]);
    assert.equal(a.solution.current, b.solution.current);
    assert.deepEqual(Array.from(a.solution.phi), Array.from(b.solution.phi));
  }
  await there.close();
});

// Each failure resolves, and leaves the device exactly as it was (found in review).
const within = (p, ms = 20000) => Promise.race([p, new Promise((_, reject) => setTimeout(() => reject(new Error('hung')), ms))]);

test('a change that fails only when the solver is built is refused, and the device carries on', async () => {
  const dev = live(pn(1e17, 1e16, 0));
  await dev.ready;
  // Valid as a definition, but h⁺ then reaches no contact and has no c0: the solver refuses it.
  const blocked = { species: { 'e-': 'equilibrium', 'h+': 'blocked' } };
  const bad = await within(dev.set({ contacts: { left: blocked, right: blocked } }));
  assert.ok(bad.info.failed && /c0/.test(bad.info.error), JSON.stringify(bad.info));
  const next = await within(dev.set({ contacts: { right: { V: 0.3 } } }));
  assert.ok(!next.info.failed && Math.abs(next.solution.terminals.right.V - 0.3) < 1e-12);
  const never = live({ ...pn(1e17, 1e16, 0), contacts: { left: blocked, right: blocked } });
  const r = await within(never.ready);
  assert.ok(r.info.failed && r.solution === null);
});

test('a failed change of drive is rolled back exactly: no leftover I or R', async () => {
  const dev = live(pn(1e17, 1e16, 0), { maxRamp: 2 });
  await dev.ready;
  const huge = await within(dev.set({ contacts: { right: { V: 1000, R: 1e-12 } } }));
  assert.ok(huge.info.failed);
  const back = await within(dev.set({ contacts: { right: { V: 0.2 } } }));
  assert.ok(!back.info.failed);
  const direct = new Device(pn(1e17, 1e16, 0.2)).solve();
  assert.ok(Math.abs(back.solution.current / direct.current - 1) < 1e-9, 'no series resistance left behind');
  const driven = await within(dev.set({ contacts: { right: { V: undefined, I: -1e9 } } }));
  const held = await within(dev.set({ contacts: { right: { I: undefined, V: 0.1 } } }));
  assert.ok(driven.info && !held.info.failed && Math.abs(held.solution.terminals.right.V - 0.1) < 1e-12);
});

test("an invalid change made with valid ones is refused alone, and a throwing onsolution doesn't stop the device", async () => {
  let calls = 0;
  const dev = live(pn(1e17, 1e16, 0), {
    onsolution() {
      calls++;
      if (calls === 1) throw new Error('drawing failed');
    },
  });
  const uncaught = [];
  const saved = globalThis.reportError;
  globalThis.reportError = (e) => uncaught.push(e.message);
  try {
    const [, good, bad] = await within(
      Promise.all([dev.ready, dev.set({ contacts: { right: { V: 0.3 } } }), dev.set({ bulkReactions: [{ nu: { 'e-': -1, 'h+': -1 }, kf: { Si: 'fast' } }] })]),
    );
    assert.ok(!good.info.failed && Math.abs(good.solution.terminals.right.V - 0.3) < 1e-12);
    assert.ok(bad.info.failed && /kf/.test(bad.info.error));
    await new Promise((r) => setTimeout(r, 10));
    assert.deepEqual(uncaught, ['drawing failed']);
  } finally {
    globalThis.reportError = saved;
  }
});
