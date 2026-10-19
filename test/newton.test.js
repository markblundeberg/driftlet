import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, metal } from '../src/kit.js';
import { newton } from '../src/newton.js';
import { choice, device } from '../tools/lawsweep-devices.js';

// Newton's iteration (src/newton.js) against its contract, on what a cold review of it found.
// The lost mode a time step's plain solves can't see is in test/ports.test.js (MOS without a
// port: the gate's charge after the step, to 1e-6 of the steady one).

// Each steady newton() call on a device, with whether the islands' rows were in when it returned.
const steadyCalls = (dev) => {
  const s = dev.solver, calls = [], orig = s.newton.bind(s);
  s.newton = (dt, opts) => {
    const r = orig(dt, opts);
    if (dt === Infinity) calls.push({ ...r, asked: !!opts?.islands, on: s.islandsOn });
    return r;
  };
  return calls;
};

test('newton: at equilibrium the terminal currents are round-off, and the islands round stays off', () => {
  // Every current 1e-20 to 1e-31 A/m², their sum as large as they are: off Kirchhoff by O(1)
  // of the largest, which put the islands' rows in, twice, for nothing (2 of 6 iterations).
  for (const i of [288, 365]) {
    const dev = new Device(device(choice(i)));
    const calls = steadyCalls(dev), sol = dev.solve();
    assert.ok(sol.converged && dev.solver.islands.length > 0, `device ${i}`);
    assert.ok(calls.every((c) => !c.on), `device ${i}: ${JSON.stringify(calls.map((c) => [c.iterations, c.on]))}`);
  }
});

test('newton: a call with the islands\' rows in from the start keeps them', () => {
  // step()'s second try (islands: true) converges at its 59th iteration; the islands round ran
  // there as if the rows had just gone in, into maxIter, and came back with s.islandsOn false on
  // a state solved with them in.
  const dev = new Device(device(choice(200)));
  const calls = steadyCalls(dev), sol = dev.solve();
  assert.ok(sol.converged);
  const asked = calls.filter((c) => c.asked);
  assert.ok(asked.length > 0 && asked.every((c) => c.converged && c.on), JSON.stringify(asked.map((c) => [c.iterations, c.on])));
});

test('newton: a steady solve whose residual isn\'t solved says so, and records how conditioned it was', () => {
  // The residual check reassembles, which cleared the factorisation its conditioning note reads:
  // the note came out empty, under a message about a perturbed pivot there wasn't.
  const dev = new Device(device(choice(5))), s = dev.solver;
  const orig = s.newton.bind(s);
  let r = null;
  s.newton = (dt, opts) => {
    if (r || dt !== Infinity) return orig(dt, opts);
    s.conditioning = null;
    return (r = newton(s, dt, { ...opts, tol: 1e9 })); // (an update under 1e9 converged: the first)
  };
  dev.solve();
  assert.ok(r && !r.converged && s.sys.perturbedPivots === 0, JSON.stringify(r));
  assert.match(r.error, /residual is .* of its largest term/);
  assert.ok(s.conditioning?.digits > 0, JSON.stringify(s.conditioning));
});

test('newton: a floating gate at rest in inversion stays put (a mode the plain solves never see, found by a probe)', () => {
  // Its electrons reach the back contact only by minority diffusion, so the factorised solves
  // lose the inversion layer's level; Newton converged quadratically on the rest at every step,
  // and the extrapolated guess's error in the level grew: here 5e-6 V by 1e7 s, through 63,000
  // steps (8.5e-5 V by 1e5 s, reached by way of 1e3 s).
  const def = build({
    T: 300,
    library: [semiconductor('Si'), metal('Au')],
    materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [
      ohmic(0),
      layer('Si', units.um(0.5), { name: 'p-Si', acceptors: units.perCm3(1e17) }),
      { dipole: 0 },
      layer('SiO2', units.nm(5), { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }),
      { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.1 },
      layer('Au', units.nm(50)),
      ohmic(0, ['e-']),
    ],
    grid: { hmin: units.nm(0.1), hmax: units.nm(10), ratio: 1.1 },
  });
  const dev = new Device(def).set({ contacts: { right: { V: 1.2 } } });
  dev.solve();
  dev.set({ contacts: { right: { V: undefined, I: 0 } } });
  for (const t of [1e5, 1e7]) {
    const s = dev.advance(t);
    assert.ok(s.converged && s.steps < 100 && Math.abs(s.contacts.right.V - 1.2) < 1e-9, `t = ${t}: ${s.contacts.right.V - 1.2} V, ${s.steps} steps`);
  }
});
