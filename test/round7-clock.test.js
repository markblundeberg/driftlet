import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, SolverError } from '../src/index.js';
import { aqueous, bath, build, half, hodgkinHuxley, IONS, layer, metal, ohmic, square } from '../src/kit.js';
import { sourceAt } from '../src/device.js';
import { gateRate } from '../src/boundaries.js';

// Faults in the time integrator (src/timeSteps.js and its callers) found by round 7's outside
// reading pass (OpenAI's Sol, re-reading it), each confirmed by a chaser's probe before it was
// fixed.

const E = Number.EPSILON;
const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
// 1 Ω·m²: σ = 1e-6 S/m over 1 µm, so I = V and a ramp's charge is its integral.
const resistor = () =>
  new Device({
    species: [{ name: 'e-', z: -1 }],
    materials: { R: { conductor: { species: 'e-', conductivity: 1e-6 } } },
    regions: [{ material: 'R', length: 1e-6 }],
    contacts: { left: collector(0), right: collector(0) },
    grid: { hmin: 1e-7 },
  });
const wire = (right = {}) =>
  new Device({
    species: [{ name: 'e-', z: -1 }],
    materials: { Cu: { conductor: { species: 'e-', conductivity: 1e6 } } },
    regions: [{ material: 'Cu', length: 1e-6 }],
    contacts: { left: collector(0), right: { ...collector(0), ...right } },
    grid: { hmin: 1e-7 },
  });

test('S1: a resistor ramped after a rest passes the ramp\'s charge, though nothing in the state feels it', () => {
  // (The first step's pair agreed on the end state and spanned the ramp: 0.75 for 0.5 C/m².)
  for (const [drive, exact, opts] of [
    [{ V: { t: [0, 10, 11], values: [0, 0, 1] } }, 0.5, {}],
    [{ V: { t: [0, 1], values: [0, 1] } }, 0.5, { dt0: 1 }],
    [{ V: undefined, I: { t: [0, 10, 11], values: [0, 0, 2] } }, 1, {}],
  ]) {
    const d = resistor();
    d.solve();
    d.set({ contacts: { right: drive } });
    const r = d.advance(drive.V ? drive.V.t.at(-1) : drive.I.t.at(-1), opts);
    const q = r.trace.terminals.right.charge.at(-1);
    assert.ok(r.converged && Math.abs(q / exact - 1) < 1e-4, `${JSON.stringify(drive)}: ${q} vs ${exact}`);
  }
});

test('S1: a current ramp under backward Euler passes its charge to first order (26% high before), and a zero crossing doesn\'t stall', () => {
  const d = resistor();
  d.solve();
  d.set({ contacts: { right: { V: undefined, I: { t: [0, 1], values: [0, 2] } } } });
  const r = d.advance(1, { method: 'be' });
  assert.ok(r.converged && r.done && Math.abs(r.trace.terminals.right.charge.at(-1) - 1) < 0.03, `${r.trace.terminals.right.charge.at(-1)}`);
  // A triangle 0 → 1 → −1 V: the current crosses zero inside a segment, where a charge's error
  // relative to the step's own charge never shrinks. Its net charge is 0.5 C/m² (before: −0.27
  // and −0.40, the wrong sign).
  for (const [method, within] of [['be', 0.05], ['bdf2', 0.005]]) {
    const e = resistor();
    e.solve();
    e.set({ contacts: { right: { V: { t: [0, 1, 3], values: [0, 1, -1] } } } });
    const s = e.advance(3, { method });
    assert.ok(s.converged && s.done && s.steps < 1000 && Math.abs(s.trace.terminals.right.charge.at(-1) - 0.5) < within, `${method}: ${s.trace.terminals.right.charge.at(-1)} in ${s.steps} steps`);
  }
});

test('S1: a slow CV after a rest passes the charge the quasi-steady currents do (52% high at its vertex before, the cycle\'s net of the wrong sign)', () => {
  // A redox couple at Pt under 3.16 µm of stagnant electrolyte, fed by a bath: L²/D ~ 0.01 s,
  // so at 1 mV/s the current is the steady one at each potential. 24136 C/m² to the first vertex
  // by two routes (BDF2 at period/2e4, tol 1e-6: 24136.05; solve() at 401 potentials and Simpson:
  // 24135.80), and 5497 C/m² net over the cycle (the first route); round 7's chaser.
  const iron = half('Fe3+ + e- = Fe2+'), T = 298.15, RT = GAS_CONSTANT * T, L = 3.16e-6;
  const b = { 'K+': 500, 'Cl-': 525, 'Fe3+': 5, 'Fe2+': 5 };
  const she = -(IONS['Cl-'].mu0 + RT * Math.log(b['Cl-'] / 1000)) / FARADAY;
  const electrode = (drive) =>
    build({
      library: [aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 }), metal('Pt')],
      stack: [ohmic(drive, ['e-']), layer('Pt', 1e-6), { reactions: [{ ...iron, k0: 1e-3, alpha: 0.5 }] }, layer('water', L, { name: 'stagnant layer' }), bath(b, 'Cl-', 0, { offset: she })],
      grid: { hmin: 10e-9, hmax: L / 40 },
    });
  const Erest = new Device(electrode({ I: 0 })).solve().terminals.left.V;
  const hold = 1000, period = 2000;
  const d = new Device(electrode(Erest));
  d.solve();
  d.set({ contacts: { left: { V: { t: [0, hold, hold + period / 4, hold + (3 * period) / 4, hold + period], values: [Erest, Erest, Erest + 0.5, Erest - 0.5, Erest] } } } });
  const r = d.advance(hold + period);
  const at = (t) => r.trace.terminals.left.charge[r.trace.t.findIndex((x) => x >= t)];
  assert.ok(r.converged && Math.abs(at(hold + period / 4) / 24136 - 1) < 2e-3, `first vertex ${at(hold + period / 4)}`);
  assert.ok(Math.abs(at(hold + period) / 5497 - 1) < 2e-2, `net ${at(hold + period)}`);
});

test('S1: at rest, the terminals\' charges are noise, and cut no steps', () => {
  // (A blocking cell's plates rattle by 1e-13 C/m², a held electrode's currents by a tenth of
  // 4e-16 A/m², a strictly neutral contact's box by 5e-11 C/m²: each cut the steps without end
  // once the charge was in the error, until it was weighed against what it's known to.)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), plate = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 };
  const d = new Device({
    species: lib.species, materials: lib.materials,
    regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, ...plate }, right: { V: 0, ...plate } },
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
  d.solve();
  for (const tol of [1e-3, 1e-8]) {
    const r = d.advance(d.solution().time + 10, { tol });
    assert.ok(r.converged && r.steps < 40 && r.rejected === 0, `tol ${tol}: ${r.steps} steps, ${r.rejected} rejected`);
  }
});

test('S13: a small current after a large charge passed still adds up (compensated)', () => {
  // (1e8 C/m² posted, then 1 mA/m² in 1 ms steps: 1.6e-3 short; after 1e12, nothing at all.)
  for (const [big, tBig, small, rel] of [[1e4, 1e4, 1e-3, 1e-5], [1e12, 1, 1e-3, 0.13]]) {
    const d = wire({ V: undefined, I: big });
    d.solve();
    d.step(tBig);
    const P0 = d.solver.termCharge[1];
    d.set({ contacts: { right: { I: small } } });
    const r = d.advance(tBig + 1, { dt0: 1e-3, dtMax: 1e-3 });
    const dq = r.trace.terminals.right.charge.at(-1) - P0;
    // (to the running total's ulp: 1.2e-4 at 1e12, a tenth of what passed)
    assert.ok(Math.abs(dq / small - 1) < rel, `after ${big * tBig}: ${dq} vs ${small}`);
  }
});

test('S3: a solution read a few ulps past a current\'s edge reports what the state passed', () => {
  // (A square wave advanced by half periods: a call's end snapped a hair past an edge, and
  // sol.terminals read the drive there, 1 A/m², where the state, solved before it, passed none.)
  const d = wire({ V: undefined, I: 0 });
  d.solve();
  d.set({ contacts: { right: { I: square({ high: 1, low: 0, period: 0.1, duty: 0.5 }) } } });
  let t = 0, snapped = 0;
  for (let k = 0; k < 240; k++) {
    t += 0.05;
    const sol = d.advance(t);
    if (d.solver.sourceTime !== d.solver.time) snapped++;
    assert.ok(Math.abs(sol.terminals.right.current - sol.flux['e-'].at(-1) * FARADAY) < 1e-6, `t = ${t}: terminal ${sol.terminals.right.current}, flux ${sol.flux['e-'].at(-1) * FARADAY}`);
  }
  assert.ok(snapped > 0, 'some calls ended a hair past an edge');
});

test('S7: a long step after one below 1e-12 of it extrapolates linearly, and the next call runs', () => {
  // (History of 1e5 s and 5e-13 s: the quadratic's offsets coincided, its weights 0/0, and every
  // later call failed "non-finite update".)
  for (const dt0 of [1e5, 8192]) {
    const d = wire();
    d.solve();
    d.advance(1e-12);
    d.advance(1e5, { dt0 });
    const r = d.advance(2e5);
    assert.ok(r.converged && r.done, `dt0 ${dt0}: ${r.warnings.join(' ')}`);
  }
});

test('S11: a first step tried again at the same size is not a step that can\'t be cut', () => {
  // (Three calls of maxSteps 1, each from dt0, pushed one error three times, and the third took
  // the step at 7e5 times tol, its trace 13% to 46% off.)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), plate = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 };
  const d = new Device({
    species: lib.species, materials: lib.materials,
    regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, ...plate }, right: { V: 0, ...plate } },
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
  d.solve();
  d.set({ contacts: { right: { V: 0.05 } } });
  for (let k = 0; k < 3; k++) d.advance(1e-4, { dt0: 1e-5, tol: 1e-7, maxSteps: 1 });
  assert.equal(d.solver.time, 0);
});

test('S11: a stiff gate\'s first step after a clamp edge is held to tol, not taken as a jump', () => {
  // (HH gates at 6.3 °C, clamped from +100 to −120 mV: the step-doubling error of m doesn't shrink
  // until (α + β)h < 1, three tries looked like a jump, and the step was taken at 88 times tol,
  // m 0.442 for 0.281 at 6.015 ms; round 7's state chaser.)
  const mV = 1e-3, S = { 'Na+': 100, 'K+': 100, 'Cl-': 200 }, T = 279.45;
  const G = hodgkinHuxley({ inside: { 'Na+': 50, 'K+': 400, 'Cl-': 52 }, outside: { 'Na+': 440, 'K+': 20, 'Cl-': 450 }, T }).gates;
  const membrane = { phi: { type: 'capacitive', C: 0.01 }, gates: G, species: { 'K+': { type: 'permeability', P: 1e-12, gates: { n: 4 } }, 'Na+': 'blocked', 'Cl-': 'blocked' } };
  const wave = { t: [0, 1e-3, 1e-3, 6e-3, 6e-3, 1], values: [-65 * mV, -65 * mV, 100 * mV, 100 * mV, -120 * mV, -120 * mV] };
  const d = new Device(build({ T, library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 0 })], stack: [bath(S, 'Cl-'), layer('water', 1e-6, { c0: S }), membrane, layer('water', 1e-6, { c0: S }), bath(S, 'Cl-', wave)], grid: { hmin: 50e-9, hmax: 200e-9 } }));
  d.solve();
  const probes = [{ interface: 0, gate: 'm' }];
  const m0 = d.advance(6e-3, { probes }).trace.probes[0].at(-1);
  const s = d.advance(7e-3, { probes });
  const V = s.interfaces[0].V, a = gateRate(G.m.alpha, V)[0], b = gateRate(G.m.beta, V)[0];
  for (const tt of [6.005e-3, 6.015e-3, 6.03e-3]) {
    const i = s.trace.t.findIndex((t) => t >= tt), t = s.trace.t[i];
    const exact = a / (a + b) + (m0 - a / (a + b)) * Math.exp(-(a + b) * (t - 6e-3));
    assert.ok(Math.abs(s.trace.probes[0][i] - exact) < 0.01, `t ${t}: m ${s.trace.probes[0][i]} vs ${exact}`);
  }
});

test('S8, S9: a step that can\'t move the clock isn\'t taken, nor two rows put at one time', () => {
  // (A first step of 1e-18 s at t = 1 s posted its charge with the clock unmoved; the pair's
  // first half of a step of two ulps snapped to its landing, and the second stepped no time.)
  const d = wire({ V: undefined, I: 0 });
  d.solve();
  d.step(1);
  d.set({ contacts: { right: { I: 1 } } });
  const P0 = d.solver.termCharge[1];
  const r = d.advance(1.00001, { dt0: 1e-18, maxSteps: 1 });
  assert.ok(r.trace.t.every((t) => t > 1), `trace.t ${r.trace.t}`);
  assert.ok(Math.abs(d.solver.termCharge[1] - P0 - (d.solver.time - 1)) < 1e-20, `charge ${d.solver.termCharge[1] - P0} over ${d.solver.time - 1} s`);
  // Fixed: step(ε) at t = 1, its first Newton made to fail: the halves can't move the clock.
  const f = wire({ V: undefined, I: 0 });
  f.solve();
  f.step(1);
  f.set({ contacts: { right: { I: 1 } } });
  const s = f.solver, newton = s.newton.bind(s);
  let n = 0;
  s.newton = (dt, o) => (++n === 1 ? { converged: false, iterations: 1 } : newton(dt, o));
  const Q0 = s.termCharge[1], sol = f.step(E);
  s.newton = newton;
  assert.ok(!sol.converged && s.time === 1 && s.termCharge[1] === Q0, `converged ${sol.converged}, t ${s.time}, charge ${s.termCharge[1] - Q0}`);
  // S9: two ulps in one call.
  const g = wire({ V: undefined, I: 0 });
  g.solve();
  g.step(1);
  g.set({ contacts: { right: { I: 1 } } });
  const t = g.advance(1 + 2 * E, { dt0: 2 * E }).trace.t;
  assert.ok(t.every((x, k) => k === 0 || x > t[k - 1]), `trace.t ${t.map((x) => (x - 1) / E)} ulps past 1`);
});

test('S10, S12, S18: late clocks: a pulse inside a fixed step is warned of, a long period\'s left limit read, a step past the largest time refused', () => {
  const pulse = { t: [0, 1e9 - 6e-4, 1e9 - 6e-4, 1e9 - 2e-4, 1e9 - 2e-4], values: [0, 0, 1, 1, 0] };
  const d = wire({ V: undefined, I: 0 });
  d.solve();
  d.step(1e9 - 1);
  d.set({ contacts: { right: { I: pulse } } });
  assert.ok(d.step(1).warnings.some((w) => /inside this step/.test(w)));
  // (A period of 1e16 s read its left limit 18 s back, in its first cycle.)
  assert.ok(Math.abs(sourceAt({ t: [0, 1, 1, 1e16], values: [0, 0, 1, 0], repeat: true }, 2, true) - 1) < 1e-12);
  const e = wire();
  e.solve();
  e.step(1e308);
  assert.throws(() => e.step(1e308), (err) => err instanceof SolverError && /past the largest time/.test(err.message));
  assert.equal(e.solver.time, 1e308);
});
