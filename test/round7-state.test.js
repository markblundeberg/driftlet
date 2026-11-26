import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { aqueous, bath, build, check, hodgkinHuxley, layer, ohmic, semiconductor } from '../src/kit.js';
import { gateRate } from '../src/boundaries.js';
import { STATE } from '../src/solver.js';

// Round 7: faults in the state, the readout and the rollback around the time integrator, from an
// outside reading pass (OpenAI's Sol, on src/timeSteps.js), each confirmed by a chaser's probe
// against a twin or a second route before it was fixed.

const mV = 1e-3;
// A wire between two collectors, the right one driven as given.
const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const wire = (right = {}) =>
  new Device({
    species: [{ name: 'e-', z: -1 }],
    materials: { Cu: { conductor: { species: 'e-', conductivity: 1e6 } } },
    regions: [{ material: 'Cu', length: 1e-6 }],
    contacts: { left: collector(0), right: { ...collector(0), ...right } },
    grid: { hmin: 1e-7 },
  });
// KCl between two blocking capacitive plates, the right one driven as given.
const blocking = (right = { V: 0 }, grid = {}) => {
  const { species, materials } = aqueous(['K+', 'Cl-'], { epsr: 0 });
  const plate = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 };
  return new Device({
    species, materials,
    regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, ...plate }, right: { ...right, ...plate } },
    grid: { hmin: 2e-6, hmax: 5e-6, ...grid },
  });
};
// The solver's state between calls, field by field, exactly (STATE in src/solver.js).
const canon = (v) => (v === undefined ? 'undefined' : ArrayBuffer.isView(v) ? `[${Array.from(v).join()}]` : JSON.stringify(v, (k, x) => (ArrayBuffer.isView(x) ? Array.from(x) : x)));
const stateOf = (s) => STATE.map((k) => `${k}=${canon(s[k])}`);
// Newton made to fail: the real one runs (and moves the state), then says it didn't converge.
const breakNewton = (d) => {
  const s = d.solver, real = s.newton.bind(s);
  s.newton = (...a) => ({ ...real(...a), converged: false });
  return () => delete s.newton;
};

test('S6: gate fractions stay in [0, 1] under fixed BDF2, near the exact relaxation (the clamp fixture)', () => {
  // test/gates.test.js's clamp: a tiny n⁴ link between Cl⁻-referenced baths, the kit's HH gates
  // at 6.3 °C, clamped from +50 to −100 mV and stepped by 0.1 ms. BDF2's extrapolated target
  // overshot m's fast swing (h(α + β) = 2.8, its roots complex): m went to −2.7e-2 at the third
  // step, where the exact relaxation is at 7.7e-4.
  const IONS = ['Na+', 'K+', 'Cl-'], S = { 'Na+': 100, 'K+': 100, 'Cl-': 200 }, T = 273.15 + 6.3;
  const G = hodgkinHuxley({ inside: { 'Na+': 50, 'K+': 400, 'Cl-': 52 }, outside: { 'Na+': 440, 'K+': 20, 'Cl-': 450 }, T }).gates;
  const tiny = { phi: { type: 'capacitive', C: 0.01 }, gates: G, species: { 'K+': { type: 'permeability', P: 1e-12, gates: { n: 4 } }, 'Na+': 'blocked', 'Cl-': 'blocked' } };
  const d = new Device(build({ T, library: [aqueous(IONS, { epsr: 0 })], stack: [bath(S, 'Cl-'), layer('water', 1e-6, { c0: S }), tiny, layer('water', 1e-6, { c0: S }), bath(S, 'Cl-', 50 * mV)], grid: { hmin: 50e-9, hmax: 200e-9 } }));
  const x0 = d.solve().interfaces[0].gates;
  d.set({ contacts: { right: { V: -100 * mV } } });
  const dt = 1e-4;
  for (let k = 1; k <= 10; k++) {
    const sol = d.step(dt, { method: 'bdf2' }), V = sol.interfaces[0].V, x = sol.interfaces[0].gates;
    assert.ok(sol.converged);
    for (const g of ['m', 'h', 'n']) assert.ok(x[g] >= 0 && x[g] <= 1, `step ${k}: ${g} = ${x[g]}, a fraction`);
    const a = gateRate(G.m.alpha, V)[0], b = gateRate(G.m.beta, V)[0], inf = a / (a + b);
    const exact = inf + (x0.m - inf) * Math.exp(-(a + b) * k * dt);
    // (from the third step: the first is backward Euler's, 0.26 where 0.062 is exact, a fast gate
    // a step can't follow; by the third BDF2 is within 4.2e-4, and 7e-5 by the fifth)
    if (k >= 3) assert.ok(Math.abs(x.m - exact) < 1e-3, `step ${k}: m ${x.m} against the exact relaxation's ${exact}`);
  }
});

test("S4: a capacitive contact's displacement isn't posted as the ionic intake of the stretch it feeds", () => {
  // A KCl cell, its left contact a capacitor that also passes K⁺, driven by 1 A/m²: the current
  // is the plate's charging and the K⁺ intake together. Posted whole as the stretch's, the ledger
  // showed 1.04e-7 mol/m² taken in where the K⁺ amount moved by 4.9e-9, and check() failed.
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 });
  const d = new Device({
    species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: {
      left: { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, species: { 'K+': { type: 'conductance', G: 1, offset: 0 }, 'Cl-': 'blocked' }, I: 0 },
      right: { V: 0, terminal: 'Cl-', phi: 'neutral', species: { 'Cl-': 'equilibrium', 'K+': 'blocked' } },
    },
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
  const K0 = d.solve().conservation.find((c) => c.species === 'K+');
  d.set({ contacts: { left: { I: 1 } } });
  let sol, flux = 0;
  for (let k = 0; k < 10; k++) {
    sol = d.step(1e-3);
    flux += sol.contacts.left.flux['K+'] * 1e-3;
  }
  const K = sol.conservation.find((c) => c.species === 'K+');
  const moved = K.amount - K0.amount;
  assert.ok(moved > 1e-9 && moved < 1e-8, `K⁺ taken in: ${moved} mol/m² (most of I·dt/F = 1.04e-7 charges the plate)`);
  assert.ok(Math.abs(K.intake - moved) < 1e-6 * moved, `the ledger's intake ${K.intake} against the amount's change ${moved}`);
  assert.ok(Math.abs(K.intake - flux) < 1e-6 * moved, `and against the contact's K⁺ flux, summed: ${flux}`);
  assert.ok(Math.abs(K.drift) < 1e-12, `drift ${K.drift}`);
  assert.ok(check(d, sol, { refine: false, order: false }).items.find((i) => i.name === 'conservation').ok);
});

test('S5: two drive patches with no step between leave a floating plate where it was solved', () => {
  // I: 0 → V: 0.05 → I: 0 on a floating capacitive plate: the second patch floated it from the
  // first's 0.05 V, which never acted, and the plate read 0.01 C/m² it never took.
  const mk = () => {
    const d = blocking({ I: 0 });
    d.solve();
    d.step(1e-3);
    return d;
  };
  const a = mk(), twin = mk();
  a.set({ contacts: { right: { V: 0.05 } } });
  a.set({ contacts: { right: { I: 0 } } });
  const sa = a.solution(), st = twin.solution();
  assert.equal(sa.terminals.right.V, st.terminals.right.V);
  assert.equal(sa.contacts.right.D, st.contacts.right.D);
  assert.ok(Math.abs(sa.contacts.right.D) < 1e-20, `plate charge ${sa.contacts.right.D} C/m²`);
  assert.deepEqual(Array.from(a.solver.termV), Array.from(twin.solver.termV));
  const na = a.step(1e-3), nt = twin.step(1e-3);
  assert.ok(Math.abs(na.terminals.right.V - nt.terminals.right.V) < 1e-15);
});

test("S14: a read leaves a held terminal's voltage as the step left it, and a rebuild that releases it at an edge floats it from there", () => {
  const edge = { t: [0, 1, 1, 2], values: [0, 0, 0.1, 0.1] };
  const mk = () => {
    const d = wire({ V: edge });
    d.solve();
    d.advance(1); // (onto the 0 → 0.1 V edge, solved from within the step: at 0)
    return d;
  };
  // The read re-read the held source on the far side of the edge, and left it there.
  const d = mk(), was = Array.from(d.solver.termV);
  assert.equal(d.solution().terminals.right.V, 0);
  assert.deepEqual(Array.from(d.solver.termV), was);
  assert.equal(was[1], 0);
  // A rebuild (a new grid; a new conductivity) releasing the terminal there floated it from
  // 0.1 V where the state was solved at 0, and read 0 A/m² across a 1e6 S/m wire; the release
  // by the drive alone reads 0.
  const alone = mk();
  alone.set({ contacts: { right: { I: 0 } } });
  for (const patch of [{ grid: { hmin: 0.9e-7 } }, { materials: { Cu: { conductor: { conductivity: 1.0000001e6 } } } }]) {
    const a = mk();
    a.set({ contacts: { right: { I: 0 } }, ...patch });
    assert.equal(a.solution().terminals.right.V, alone.solution().terminals.right.V, JSON.stringify(patch));
    assert.equal(a.solver.termV[1], 0, JSON.stringify(patch));
  }
});

test('S15: a zero-attempt advance leaves the history and the step size carried as they were', () => {
  // At a knot step() landed on (history 1 there): maxSteps 0 cleared the history and set dtNext
  // from its own interval.
  const wave = { t: [0, 1e-3, 1e-3, 1], values: [0, 0, 0.05, 0.05] };
  const mk = () => {
    const d = blocking({ V: wave });
    d.solve();
    d.step(1e-3);
    return d;
  };
  const a = mk(), twin = mk();
  assert.equal(a.solver.history.length, 1);
  const r = a.advance(2e-3, { maxSteps: 0 });
  assert.equal(r.steps, 0);
  assert.deepEqual(stateOf(a.solver), stateOf(twin.solver));
  // And from a fresh drive: dtNext 1e-13, from advance(1e-9, { maxSteps: 0 }), made a later
  // advance(1e-3) take 51 steps where its twin took 32.
  const mk2 = () => {
    const d = blocking({ V: 0 });
    d.solve();
    d.set({ contacts: { right: { V: 0.05 } } });
    return d;
  };
  const b = mk2(), twin2 = mk2();
  b.advance(1e-9, { maxSteps: 0 });
  assert.equal(b.solver.dtNext, undefined);
  assert.equal(b.advance(1e-3).steps, twin2.advance(1e-3).steps);
});

test("S15: a failed step() puts back the whole state, the islands' warning and the sources' time with it", () => {
  // test/round4-newton's warned floating base: a GaAs n-p-n whose side port can't balance. A
  // fixed step whose Newton fails restored a hand-written list that missed the islands' fields,
  // and the restored steady state lost its "floating region doesn't balance" warning.
  const mat = 'GaAs';
  const def = () => {
    const d = build({ T: 250, library: [semiconductor(mat, { T: 250 })], stack: [ohmic(0), layer(mat, 0.483e-6, { name: 'n1', donors: 0.66 }), layer(mat, 0.22e-6, { name: 'p', acceptors: 8.94 }), layer(mat, 0.902e-6, { name: 'n2', donors: 3.5 }), ohmic(-0.2)], grid: { hmin: 0.5e-9, hmax: 40e-9, ratio: 1.11 } });
    d.ports = [{ name: 'side', region: 'n2', from: 0.75e-6, to: 0.9e-6, I: 1e-3, terminal: 'e-', species: { 'e-': 'equilibrium' } }];
    return d;
  };
  const island = (sol) => sol.warnings.filter((w) => w.includes("floating region doesn't balance")).length;
  const a = new Device(def()), twin = new Device(def());
  assert.equal(island(a.solve()), 1);
  twin.solve();
  const undo = breakNewton(a);
  const r = a.step(1e-12);
  undo();
  assert.equal(r.converged, false);
  assert.deepEqual(stateOf(a.solver), stateOf(twin.solver));
  assert.equal(island(a.solution()), 1);
  // After a snap past an edge (the sources' time left at the edge, the clock a hair past it), a
  // failed step rebuilt the sources' time from the clock: a blocking plate then read the edge's
  // far side, 0.05 V and 0.01 C/m² that never passed.
  const T = 1e5;
  const mk = () => {
    const d = blocking({ V: { t: [0, T, T], values: [0, 0, 0.05] } });
    d.solve();
    d.step(T - 1000);
    d.advance(T + 5e-8, { dt0: 1000 });
    return d;
  };
  const b = mk(), twin2 = mk();
  const undo2 = breakNewton(b);
  const rb = b.step(1e-3);
  undo2();
  assert.equal(rb.converged, false);
  assert.deepEqual(stateOf(b.solver), stateOf(twin2.solver));
  assert.equal(b.solution().contacts.right.D, twin2.solution().contacts.right.D);
});

test('S16: a refused advance() leaves a pending regrid unsettled, and an unknown method is refused', () => {
  // A steady state carried onto a new grid is solved there by the next call that steps. Each of
  // these settled it before it was refused (or took no step), and left the device not as a twin
  // that never made the call.
  const mk = () => {
    const d = blocking({ V: 0.05 });
    d.solve();
    d.set({ grid: { hmax: 3e-6 } });
    return d;
  };
  const calls = {
    'advance(Infinity)': (d) => d.advance(Infinity),
    'tol 0': (d) => d.advance(1e-6, { tol: 0 }),
    'dtMax −1': (d) => d.advance(1e-6, { dtMax: -1 }),
    'maxSteps 0.5': (d) => d.advance(1e-6, { maxSteps: 0.5 }),
    'maxSteps 0': (d) => d.advance(1e-6, { maxSteps: 0 }),
    'a bad probe': (d) => d.advance(1e-6, { probes: [{ x: 1 }] }),
    "advance's method 'rk4'": (d) => d.advance(1e-6, { method: 'rk4' }),
    "step's method 'rk4'": (d) => d.step(1e-6, { method: 'rk4' }),
  };
  for (const [name, call] of Object.entries(calls)) {
    const d = mk(), twin = mk();
    assert.equal(d.solver.regriddedSteady, true);
    try {
      call(d);
    } catch {
      // (refused: what's checked is what it left)
    }
    assert.equal(d.solver.regriddedSteady, true, name);
    assert.deepEqual(stateOf(d.solver), stateOf(twin.solver), name);
  }
  // ({ method: 'rk4' } ran backward Euler, silently)
  assert.throws(() => blocking().step(1e-6, { method: 'rk4' }), /method is 'bdf2' .* or 'be' .*; got 'rk4'/);
  assert.throws(() => blocking().advance(1e-6, { method: 'rk4' }), /method is 'bdf2' .* or 'be' .*; got 'rk4'/);
});
