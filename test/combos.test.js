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

// A gate floated after a held solve, over a strictly neutral region whose c0 the held gate's charge
// balanced: Newton from the held state fails here (an 8 V drop across the other gate: see ROADMAP),
// and the fallback's cold start from c0, which no held gate balances now, was refused, the refusal
// of a c0 the call never touched escaping solve() as a throw (the neutral family, case 36). Now
// it's a failure that says so, and the device is left as it was.
test('a fallback cold start the definition refuses under the drives now is a failure that says why, not a throw', () => {
  const ions = ['Na+', 'Cl-'];
  const def = {
    T: 298.15,
    species: [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { aq: { epsr: 0, species: { 'Na+': { D: 9.355768282487988e-10, cRef: 1000, mu0: 0 }, 'Cl-': { D: 1.9574095568053886e-9, cRef: 1000, mu0: 0 } } } },
    regions: [[5.093560956669155e-7, 14.010739005746752, 10.332847252985113], [3.4283965857435864e-7, 17.87106323219284, 17.87106323219284], [9.9437550633902e-7, 97.56118457059594, 76.58787514581395]].map(([length, na, cl], k) => ({ name: `r${k}`, material: 'aq', length, c0: { 'Na+': na, 'Cl-': cl } })),
    interfaces: [{ phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'equilibrium' } }, { phi: 'neutral', species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium' } }],
    contacts: Object.fromEntries(['left', 'right'].map((s) => [s, { phi: 'neutral', species: Object.fromEntries(ions.map((i) => [i, 'blocked'])) }])),
    ports: [['g0', 'r0', 0.039164065523073055, 560347.6791942775, 0.7281814346089959], ['g1', 'r2', 0.1886057707015425, 298670.0141357943, 0.81847935449332]].map(([name, region, V, area, C]) => ({ name, region, V, area, capacitance: { C } })),
    grid: { minCells: 8 },
  };
  const d = new Device(def), held = d.solve();
  assert.ok(held.converged);
  d.set({ ports: { g0: { I: 0 } } });
  const s = d.solve();
  if (s.converged) return; // (should Newton learn to float it, nothing's left to say: at these digits it doesn't; rounded to 4, it does)
  assert.match(s.warnings.join(' '), /cold start its later stages fall back on can't be laid out under the drives now: regions\[0\] \(r0\): its initial composition carries a net charge/);
  assert.ok(Math.abs(d.solution().ports[0].charge - held.ports[0].charge) < 1e-12, 'left as it was');
});
