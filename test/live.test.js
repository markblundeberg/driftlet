import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Worker } from 'node:worker_threads';
import { Device, units } from '../src/index.js';
import { build, layer, ohmic, live } from '../src/kit.js';

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
  const dev = live(pn(1e15, 1e15, 0));
  await dev.ready;
  // From light doping at equilibrium to heavy doping at bias: Newton fails from the old state.
  const warm = new Device(pn(1e15, 1e15, 0));
  warm.solve();
  assert.ok(!warm.set(pn(1e19, 1e19, 0.6)).solve().converged, 'needs the ramp (else this test checks nothing)');
  const jump = await dev.set(pn(1e19, 1e19, 0.6));
  assert.ok(!jump.info.failed && jump.info.ramp >= 2, JSON.stringify(jump.info));
  const direct = new Device(pn(1e19, 1e19, 0.6)).solve();
  assert.ok(Math.abs(jump.solution.current / direct.current - 1) < 1e-8);

  const bad = await dev.set({ contacts: { right: { V: 'high' } } });
  assert.ok(bad.info.failed && /contacts\.right\.V/.test(bad.info.error));
  assert.equal(bad.solution, jump.solution);
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
