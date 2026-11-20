// Faults the combinations harness (bench/combos) found, each with the case that shows it.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';

const electrode = () => ({ V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
const cell = () => ({
  T: 298.15,
  species: [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }],
  materials: { water: { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
  regions: [{ name: 'w', material: 'water', length: 20e-6, c0: { 'NO3-': 10 } }],
  contacts: { left: electrode(), right: electrode() },
  grid: { hmin: 1e-8, hmax: 0.5e-6 },
});
// The device's n-th linear solve from now throws (as a user's statistics model might, mid-solve).
const faultAt = (dev, n) => {
  const s = dev.solver, orig = s._solveLinear;
  let k = 0;
  s._solveLinear = function (...a) {
    if (++k < n) return orig.apply(this, a);
    s._solveLinear = orig;
    throw new Error('injected fault');
  };
};
// How many linear solves the steady solve inside an impedance() takes, from this state.
const steadySolves = (make) => {
  const dev = make(), s = dev.solver, orig = s._solveLinear, steady = s.solveSteady;
  let k = 0, n = 0;
  s._solveLinear = function (...a) { k++; return orig.apply(this, a); };
  s.solveSteady = function (...a) { const r = steady.apply(this, a); n = k; return r; };
  dev.impedance([10]);
  return n;
};

test('impedance(): a throw inside its steady solve leaves the device as it was', () => {
  const make = () => {
    const d = new Device(cell());
    d.solve();
    d.set({ contacts: { right: { V: 0.02 } } }); // (a step pending: the next advance is its response)
    return d;
  };
  const n = steadySolves(make);
  assert.ok(n >= 2);
  const dev = make(), before = dev.solution();
  faultAt(dev, n); // (the steady solve's last: the state has moved by then)
  assert.throws(() => dev.impedance([10]), /injected fault/);
  const after = dev.solution();
  assert.equal(after.time, before.time);
  assert.deepEqual(Array.from(after.mu['Ag+']), Array.from(before.mu['Ag+']));
  // …and what it does next is what a twin that never called does
  const twin = make();
  const a = dev.advance(1e-3), b = twin.advance(1e-3);
  assert.equal(a.steps, b.steps);
  assert.equal(a.current, b.current);
});

test('advance() after a regrid: a throw while it settles the steady state there leaves the device as it was', () => {
  const make = () => {
    const d = new Device(cell());
    d.solve();
    d.set({ grid: { hmin: 1e-8, hmax: 0.25e-6 } }); // (carried onto a new grid: settled at the next advance)
    return d;
  };
  // (the settle's linear solves, counted on a twin)
  const probe = make(), s = probe.solver, orig = s._solveLinear, steady = s.solveSteady;
  let k = 0, n = 0;
  s._solveLinear = function (...a) { k++; return orig.apply(this, a); };
  s.solveSteady = function (...a) { const r = steady.apply(this, a); n ||= k; return r; };
  probe.advance(1e-4);
  assert.ok(n >= 2);
  const dev = make(), before = dev.solution();
  faultAt(dev, n);
  assert.throws(() => dev.advance(1e-4), /injected fault/);
  assert.deepEqual(Array.from(dev.solution().mu['Ag+']), Array.from(before.mu['Ag+']));
  const a = dev.advance(1e-4), b = make().advance(1e-4);
  assert.equal(a.steps, b.steps);
  assert.equal(a.current, b.current);
});
