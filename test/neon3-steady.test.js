// Neon's third review, flag 1.6: a steady solve's fallback cold starts, every one through _layOut,
// so that a definition refused under the drives now is a failure that says so, not a throw out of
// solve(); a refused layout's solve counted once; and the failure's warning the likeliest cause.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY } from '../src/index.js';

// A hole contact pinned with zeroCharge 20 V behind 1e6 Ω·m²: its edge density cRef·e^(F·zeroCharge/RT)
// is e^778 at any drive (docs/device.md), so the definition has no steady state; behind R, the
// constructor's check of a pinned contact passes it by (a floating one's), and the resistance's
// held drive, R taken off, laid it out with initFromComposition and threw its refusal out of solve().
const pinned20 = () => ({
  species: [{ name: 'h+', z: 1, cRef: 10 }, { name: 'e-', z: -1, cRef: 10 }],
  materials: { semi: { epsr: 10, species: { 'h+': { D: 1e-4, mu0: 0 }, 'e-': { D: 1e-4, mu0: 0 } } } },
  regions: [{ name: 'layer', material: 'semi', length: 2e-7, c0: { 'h+': 1, 'e-': 1 } }],
  contacts: {
    left: { V: 20, R: 1e6, phi: 'pinned', zeroCharge: 20, species: { 'h+': { type: 'equilibrium', offset: 0 }, 'e-': 'blocked' } },
    right: { V: 0, phi: 'pinned', zeroCharge: 0, species: { 'h+': 'blocked', 'e-': { type: 'equilibrium', offset: 0 } } },
  },
  grid: { hmax: 2e-8 },
});

test('a cold start the definition refuses at a resistance\'s held drive: a failure that names the refusal, not a throw', () => {
  const d = new Device(pinned20());
  const before = d.solution();
  const s = d.solve();
  assert.equal(s.converged, false);
  assert.match(s.warnings.join(' '), /can't be laid out under the drives now: contacts\.left: holds h\+ at cRef·e\^778/, s.warnings.join('\n'));
  // (left as it was)
  assert.deepEqual(Array.from(d.solution().phi), Array.from(before.phi));
});

test('a cold continuation whose layout is refused counts its failed solve once', () => {
  // (q stayed the failed solve, its steps and iterations added before the layout and again after;
  // pinned to homotopy: false since _continuation still serves non-direct devices and rollbacks)
  const d = new Device(pinned20()), sv = d.solver;
  let steps = 0, iterations = 0, refusedIn = [];
  const steadySteps = sv._steadySteps.bind(sv), layOut = sv._layOut.bind(sv), from = sv._steadyFromHere.bind(sv);
  sv._steadyFromHere = (opts) => from({ ...opts, homotopy: false });
  sv._steadySteps = (opts) => {
    const q = steadySteps(opts);
    steps += q.steps;
    iterations += q.iterations;
    return q;
  };
  sv._layOut = (budget) => {
    const ok = layOut(budget);
    if (!ok) refusedIn.push(/at Solver\.(\w+)/.exec(new Error().stack.split('\n')[2])?.[1]);
    return ok;
  };
  const s = d.solve();
  assert.equal(s.converged, false);
  assert.ok(refusedIn.includes('_continuation'), `refused in ${refusedIn}`);
  assert.deepEqual([s.steps, s.iterations], [steps, iterations], 'reported against taken');

  // And on the default homotopy path, a refused layout at s = 0 also counts steps and iterations once.
  const dh = new Device(pinned20()), sh = dh.solver;
  let hSteps = 0, hIters = 0;
  const hSteadySteps = sh._steadySteps.bind(sh);
  sh._steadySteps = (opts) => {
    const q = hSteadySteps(opts);
    hSteps += q.steps;
    hIters += q.iterations;
    return q;
  };
  const shSol = dh.solve();
  assert.equal(shSol.converged, false);
  assert.deepEqual([shSol.steps, shSol.iterations], [hSteps, hIters], 'homotopy reported against taken');
});

// A reaction slower than the capped step (steady.test.js's settling case): the ramp runs out of
// steps still settling, and the maxSteps the warning names finishes it.
const slow = () => new Device({
  T: 300,
  species: [{ name: 'e-', z: -1 }, { name: 'Y', z: 0 }, { name: 'Z', z: 0 }],
  materials: {
    m1: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 1e-9, mu0: 0, cRef: 1e3 }, Z: { D: 1e-9, mu0: 0, cRef: 1e3 } } },
    m2: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 0, mu0: 0, cRef: 1e3 }, Z: { D: 1e-9, mu0: 0, cRef: 1e3 } } },
  },
  regions: [{ name: 'A', material: 'm1', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 1, Z: 0.01 } }, { name: 'B', material: 'm2', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 2, Z: 0.01 } }],
  interfaces: [{ dipole: 0 }],
  bulkReactions: [{ equation: 'Y = Z', kf: { m1: 1e-10, m2: 1e-10 } }],
  contacts: { left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
  grid: { minCells: 40 },
});

test('a ramp out of steps while still settling is said ahead of a refused layout, a refused trial state or a violent start', () => {
  // Each read beside the settling ramp, as a stage before it would have left it (a stage 1 whose
  // cold start was refused, then a pseudo-transient that ran out settling, was told only of the
  // refusal); the settling ramp is the one measured, and its maxSteps finishes it.
  const marks = {
    layoutRefused: (sv) => (sv.layoutRefused = 'regions[0] (A): its initial composition carries a net charge of 1 mol/m³'),
    modelError: (sv) => (sv.modelError = 'debye-huckel: not convex'),
    violent: (sv) => (sv.violent = ['contacts.left: 1 V from its balance']),
  };
  let n;
  for (const [name, mark] of Object.entries(marks)) {
    const d = slow(), sv = d.solver, from = sv._steadyFromHere.bind(sv);
    sv._steadyFromHere = (opts) => {
      mark(sv);
      return from(opts);
    };
    const s = d.solve();
    assert.equal(s.converged, false);
    const w = s.warnings.find((x) => /ran out of time steps/.test(x));
    assert.ok(w, `${name}: ${s.warnings.join(' ')}`);
    assert.ok(!s.warnings.some((x) => /can't be laid out|refused by a statistics model|from a start far from equilibrium/.test(x)), `${name}: ${s.warnings.join(' ')}`);
    n = +/maxSteps: (\d+)/.exec(w)[1];
  }
  assert.ok(slow().solve({ maxSteps: n }).converged, `maxSteps ${n}`);
});
