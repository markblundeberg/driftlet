import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0 } from '../src/index.js';

// Terminals: the two contacts and every port, each held at a voltage (behind a resistance, if
// given) or driven by a current, with V and I as numbers or waveforms. Currents are into the
// device at each terminal.

const Dp = 1.65e-9, Dm = 1.9e-9, c0 = 10, L = 20e-6;
const links = { terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
const silverNitrate = (extra = {}) => ({
  species: [
    { name: 'Ag+', z: 1, cRef: 1000 },
    { name: 'NO3-', z: -1, cRef: 1000 },
  ],
  materials: { water: { epsr: 0, species: { 'Ag+': { D: Dp, mu0: 77.1e3 }, 'NO3-': { D: Dm, mu0: -111.3e3 } } } },
  regions: [{ name: 'cell', material: 'water', length: L, c0: { 'NO3-': c0 } }],
  contacts: { left: { V: 0, ...links }, right: { V: 0.05, ...links } },
  grid: { minCells: 200 },
  ...extra,
});

test('a reference electrode (a port at zero current) reads the local level and changes nothing', () => {
  const plain = new Device(silverNitrate()).solve();
  const ref = { name: 'ref', region: 'cell', from: L / 2 - 1e-9, to: L / 2 + 1e-9, I: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium' } };
  const dev = new Device(silverNitrate({ ports: [ref] }));
  const sol = dev.solve();
  assert.ok(sol.converged);
  const g = dev.model.ports[0].nodes[0];
  assert.equal(dev.model.ports[0].nodes.length, 1);
  // Silver in equilibrium with the solution's Ag⁺ there: V_ref = V_Ag⁺(x).
  assert.ok(Math.abs(sol.terminals.ref.V - sol.V['Ag+'][g]) < 1e-12);
  assert.ok(Math.abs(sol.terminals.ref.current) < 1e-12 * Math.abs(plain.current));
  assert.ok(Math.abs(sol.current / plain.current - 1) < 1e-10);
  // Kirchhoff: the currents into the device sum to zero.
  const sum = Object.values(sol.terminals).reduce((a, t) => a + t.current, 0);
  assert.ok(Math.abs(sum) < 1e-10 * Math.abs(plain.current));
});

test('a port held at a voltage, or behind a resistance: I = (V_source − V)/R, and Kirchhoff', () => {
  const port = (drive) => ({ name: 'p', region: 'cell', from: 0.4 * L, to: 0.6 * L, ...drive, terminal: 'Ag+', species: { 'Ag+': { type: 'conductance', G: 1e7 } } });
  const held = new Device(silverNitrate({ ports: [port({ V: 0.02 })] })).solve();
  assert.ok(held.converged);
  const sum = (sol) => Object.values(sol.terminals).reduce((a, t) => a + t.current, 0);
  assert.ok(Math.abs(sum(held)) < 1e-10 * Math.abs(held.terminals.p.current));
  for (const R of [1e-3, 1e-2]) {
    const sol = new Device(silverNitrate({ ports: [port({ V: 0.02, R })] })).solve();
    assert.ok(sol.converged);
    const t = sol.terminals.p;
    assert.ok(Math.abs(t.current - (0.02 - t.V) / R) < 1e-10 * Math.abs(t.current), `R=${R}`);
    assert.ok(Math.abs(sum(sol)) < 1e-10 * Math.abs(t.current));
  }
});

test('a triangle wave on a series capacitor: I = C dV/dt exactly, with steps on the breakpoints', () => {
  // A dielectric between two capacitive plates (a neutral, immobile spectator keeps it a
  // device): C_s = 1/(1/C₁ + L/ε + 1/C₂), and the current follows dV/dt piece by piece.
  const C1 = 0.2, C2 = 0.5, Ld = 50e-9, epsr = 10, period = 2e-3;
  const Cs = 1 / (1 / C1 + Ld / (epsr * EPS0) + 1 / C2);
  const dev = new Device({
    species: [{ name: 'X', z: 0, cRef: 1 }],
    materials: { oxide: { epsr, species: { X: { D: 0, mu0: 0 } } } },
    regions: [{ material: 'oxide', length: Ld, c0: { X: 1 } }],
    contacts: {
      left: { V: { t: [0, period / 2, period], values: [0, 1, 0], repeat: true }, phi: { type: 'capacitive', C: C1 }, zeroCharge: 0 },
      right: { V: 0, phi: { type: 'capacitive', C: C2 }, zeroCharge: 0 },
    },
    grid: { minCells: 4 },
  });
  dev.solve();
  const run = dev.advance(1.75 * period);
  assert.ok(run.converged && run.done);
  const { t, current } = run.trace;
  for (const tb of [period / 2, period, 1.5 * period]) assert.ok(t.some((ti) => Math.abs(ti - tb) < 1e-15), `a step lands on ${tb}`);
  let prev = 0;
  t.forEach((ti, k) => {
    // The slope over the step that ended at ti (each step lies within one ramp).
    const mid = (prev + ti) / 2, rising = mid % period < period / 2;
    const want = Cs * (rising ? 1 : -1) / (period / 2);
    assert.ok(Math.abs(current[k] / want - 1) < 1e-9, `t=${ti}: ${current[k]} vs ${want}`);
    prev = ti;
  });
});

test('a sawtooth and a step: advance() runs through jumps, with I = C dV/dt between them and the jump charge after', () => {
  const C1 = 0.2, C2 = 0.5, Ld = 50e-9, epsr = 10, period = 1e-3;
  const Cs = 1 / (1 / C1 + Ld / (epsr * EPS0) + 1 / C2);
  const device = (V) =>
    new Device({
      species: [{ name: 'X', z: 0, cRef: 1 }],
      materials: { oxide: { epsr, species: { X: { D: 0, mu0: 0 } } } },
      regions: [{ material: 'oxide', length: Ld, c0: { X: 1 } }],
      contacts: { left: { V, phi: { type: 'capacitive', C: C1 }, zeroCharge: 0 }, right: { V: 0, phi: { type: 'capacitive', C: C2 }, zeroCharge: 0 } },
      grid: { minCells: 4 },
    });
  // A sawtooth: 0 → 1 V over each period, then straight back to 0 (the repeat wraps).
  const saw = device({ t: [0, period], values: [0, 1], repeat: true });
  saw.solve();
  const run = saw.advance(2.5 * period);
  assert.ok(run.converged && run.done);
  let prev = 0, charge = 0;
  run.trace.t.forEach((t, k) => {
    const wrapped = prev > 0 && Math.abs(prev / period - Math.round(prev / period)) < 1e-9; // this step began at a wrap
    if (!wrapped) assert.ok(Math.abs(run.trace.current[k] / (Cs / period) - 1) < 1e-9, `t=${t}`);
    charge += run.trace.current[k] * (t - prev);
    prev = t;
  });
  // Over whole periods the charge returns to where it started; after 2.5 the ramp is half up.
  assert.ok(Math.abs(charge / (Cs * 0.5) - 1) < 1e-9, `${charge} vs ${Cs * 0.5}`);
  // A step: two points at one time.
  const step = device({ t: [0, 1e-3, 1e-3], values: [0, 0, 0.7] });
  step.solve();
  const s = step.advance(3e-3);
  assert.ok(s.converged && s.done && s.trace.t.some((t) => t === 1e-3), 'a step lands on the jump');
  let q = 0, last = 0;
  s.trace.t.forEach((t, k) => {
    q += s.trace.current[k] * (t - last);
    last = t;
  });
  assert.ok(Math.abs(q / (Cs * 0.7) - 1) < 1e-9, `${q} vs ${Cs * 0.7}`);
  assert.throws(() => device({ t: [0, 1, 1, 1], values: [0, 1, 2, 3] }), /at most two points at one time/);
});

test('impedance at either terminal of a two-terminal device is the same', () => {
  const dev = new Device(silverNitrate());
  const fs = [0.01, 1, 100];
  const a = dev.impedance(fs, { terminal: 'right' }), b = dev.impedance(fs, { terminal: 'left' });
  fs.forEach((f, q) => {
    assert.ok(Math.abs(a.Z.re[q] / b.Z.re[q] - 1) < 1e-8 && Math.abs(a.Z.im[q] / b.Z.im[q] - 1) < 1e-6, `f=${f}`);
  });
});

test('set() with only new sources keeps the solver, and gives what a rebuild would', () => {
  const dev = new Device(silverNitrate());
  dev.solve();
  const solver = dev.solver;
  dev.set({ contacts: { right: { V: 0.08 } } });
  assert.equal(dev.solver, solver, 'drives updated in place');
  const fast = dev.solve();
  const fresh = new Device(silverNitrate({ contacts: { left: { V: 0, ...links }, right: { V: 0.08, ...links } } })).solve();
  assert.ok(Math.abs(fast.current / fresh.current - 1) < 1e-10);
  // From held to driven and back, in place too: a new kind of drive replaces the old.
  dev.set({ contacts: { right: { I: fast.terminals.right.current } } });
  assert.equal(dev.solver, solver);
  const driven = dev.solve();
  assert.ok(Math.abs(driven.terminals.right.V - 0.08) < 1e-9);
  dev.set({ contacts: { right: { V: 0.08 } } });
  assert.ok(Math.abs(dev.solve().current / fast.current - 1) < 1e-10);
  dev.set({ regions: [{ name: 'cell', material: 'water', length: 2 * L, c0: { 'NO3-': c0 } }] });
  assert.notEqual(dev.solver, solver, 'a structural change rebuilds');
});

test('set() of a parameter mid-run goes on from where the run was: its sources and its step size', () => {
  // A page that changes a rate or a diffusivity every frame rebuilds the device each time; the
  // rebuilt one must read its waveforms at the same time, and not restart from a tiny step.
  const wave = { t: [0, 1, 2], values: [0, 0.05, 0], repeat: true };
  const run = (change) => {
    const dev = new Device(silverNitrate({ contacts: { left: { V: 0, ...links }, right: { V: wave, ...links } } }));
    dev.solve();
    dev.advance(0.5, { tol: 1e-3 });
    let steps = 0;
    for (let f = 0; f < 30; f++) {
      const solver = dev.solver, V = solver.termV[1];
      if (change) {
        dev.set({ materials: { water: { species: { 'Ag+': { D: Dp * (1 + 1e-3 * (1 + (f % 2))) } } } } });
        assert.notEqual(dev.solver, solver, 'a material change rebuilds');
        assert.equal(dev.solver.termV[1], V, 'the held source is read where it was');
      }
      steps += dev.advance(0.5 + (f + 1) / 60, { tol: 1e-3 }).steps;
    }
    return steps;
  };
  const still = run(false), changing = run(true);
  assert.ok(changing <= 1.5 * still, `${changing} steps with a change every frame, ${still} without`);
});
