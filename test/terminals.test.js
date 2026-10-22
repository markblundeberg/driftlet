import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Device, EPS0, units } from '../src/index.js';
import { aqueous, build, half, layer, metal, ohmic, recombination, semiconductor } from '../src/kit.js';

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

test('impedance at a held terminal with every other one driven by a current is refused (no return path), and works at the driven one', () => {
  const dev = new Device(silverNitrate({ contacts: { left: { I: 0, ...links }, right: { V: 0, ...links } } }));
  dev.solve();
  assert.throws(() => dev.impedance([1]), /every other terminal \(left\) is driven by a current.*infinite/);
  const { Z } = dev.impedance([1, 1e3], { terminal: 'left' });
  assert.ok(Z.re.every((r) => r > 0), 'measured where the current is driven, it is finite');
});

// From round 14's fuzzer: with no terminal held (a metal electrode driven by a current, the other
// end behind a resistance), the driven one started at NaN, and its metal's electrons asked for a
// c0 that a conductor can't be given. It starts level with the resistance's source now.
test('no terminal held: a current-driven electrode beside one behind a resistance solves, IR above the held case', () => {
  const lib = aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 });
  const c = { 'K+': 100, 'Cl-': 135, 'Fe3+': 5, 'Fe2+': 10 };
  const def = (right) => ({
    species: [...lib.species, { name: 'e-', z: -1 }],
    materials: { ...lib.materials, ...metal('Pt').materials },
    regions: [{ name: 'm', material: 'Pt', length: 1e-6 }, { name: 'w', material: 'water', length: 20e-6, c0: c }],
    interfaces: [{ reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 1e-3 }] }],
    contacts: { left: { I: 1, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { ...right, bath: { c, reference: 'Cl-' } } },
    grid: { hmin: 10e-9, hmax: 1e-6 },
  });
  const held = new Device(def({ V: 0 })).solve(), behind = new Device(def({ V: 0, R: 1e-3 })).solve();
  assert.ok(held.converged && behind.converged);
  assert.ok(Math.abs(behind.terminals.left.V - (held.terminals.left.V + 1e-3)) < 1e-9, `${behind.terminals.left.V} for ${held.terminals.left.V} + 1 mV`);
});

// From round 14's fuzzer (its device sc 20, a Schottky-contacted diode): a terminal behind a
// resistance, solved cold, failed where the same device held at the voltage it passes converged.
// A cold solve now finds its load line's crossing by held solves (held at its source, then the
// resistance put back, was the first fix).
test('a terminal behind a resistance, solved cold, reaches the state the held one is in', () => {
  const def = JSON.parse(readFileSync(new URL('./fixtures/series-resistance-cold.json', import.meta.url), 'utf8'));
  const base = new Device(def).solve();
  assert.ok(base.converged);
  for (const side of ['left', 'right']) {
    const I = base.terminals[side].current, R = 0.01 / Math.max(Math.abs(I), 1e-3);
    const d = new Device(def);
    d.set({ contacts: { [side]: { V: base.terminals[side].V + I * R, R } } });
    const s = d.solve();
    assert.ok(s.converged, side);
    assert.ok(Math.abs(s.terminals[side].V - base.terminals[side].V) < 1e-9, `${side}: ${s.terminals[side].V} for ${base.terminals[side].V}`);
  }
});

// From round 15's fuzzer: a silver cell passing 36 A/m² at 50 mV, one end behind a resistance
// that drops volts. Held at 0 V on one side, 16 V behind R on the other: holding it at its source
// put 16 V across the cell, and the solve failed (now: the load line's crossing, found by held
// solves). Behind R beside a current: held at V_src − I R, Kirchhoff's I. And a terminal behind R
// started at 0 V, so a −30 V source beside a current began a transient with 30 V across it.
test('a terminal behind a resistance dropping volts solves cold, and starts level', () => {
  const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 }), L = 50e-6;
  const make = (dL, dR) => build({
    library: [aqueous(['Ag+', 'NO3-'], { epsr: 0 }), metal('Ag')],
    stack: [ohmic(dL, ['e-']), layer('Ag', 1e-6), { reactions: [{ ...silver, k0: 1e-2 }] }, layer('water', L, { c0: { 'Ag+': 10, 'NO3-': 10 } }), { reactions: [{ ...silver, k0: 1e-2 }] }, layer('Ag', 1e-6), ohmic(dR, ['e-'])],
    grid: { hmin: L / 1000, hmax: L / 20 },
  });
  const held = new Device(make(0, 0.05)).solve(), I = held.terminals.right.current;
  const R = 16 / Math.abs(I);
  const a = new Device(make({ V: 0 }, { V: 0.05 + I * R, R })).solve();
  assert.ok(a.converged && Math.abs(a.terminals.right.V - 0.05) < 1e-7, `held + R: ${a.converged} ${a.terminals.right.V}`);
  const R8 = 8 / Math.abs(I), b = new Device(make({ V: held.terminals.left.current * R8, R: R8 }, { I })).solve();
  assert.ok(b.converged && Math.abs(b.terminals.right.V - b.terminals.left.V - 0.05) < 1e-7, `R + I: ${b.converged}`);
  const c = new Device(make({ V: -30, R: 0.1 }, { I: 36.3 }));
  assert.equal(c.solver.termV[0], c.solver.termV[1]); // (level at the start)
  assert.ok(c.advance(1e-6).converged);
});

// From a scan of round 15's fuzzer diodes: a lit GaAs Schottky diode passing 0.9 A/m², its
// Schottky side behind R. Held at its source first, the state sat I R from the answer and the
// pseudo-transient relaxed it over R C: 1622 steps behind 1e4 Ω·m², 15458 behind 1e5, steps ∝ R.
// Now the load line's crossing, by held solves: the same 36 steps at any R.
test('a terminal behind a large resistance solves cold in as many steps as a small one', () => {
  const make = (left) => build({
    T: 300, library: [semiconductor('GaAs')],
    stack: [
      { V: left.V, R: left.R, terminal: 'e-', species: { 'e-': { type: 'velocity', v: 2e4 }, 'h+': { type: 'velocity', v: 1e4, offset: 0 } }, phi: 'pinned', zeroCharge: 0.86 },
      layer('GaAs', 1.53e-6, { donors: units.perCm3(7.4e16) }), layer('GaAs', 4.23e-6, { acceptors: units.perCm3(2.4e16) }), ohmic(0.224)],
    bulkReactions: [recombination({ material: 'GaAs', tau: 1.8e-9, majority: units.perCm3(7.4e16) }),
      { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { GaAs: 3.2 } }],
    grid: { hmin: units.nm(1), hmax: 1.53e-6 / 20 },
  });
  const held = new Device(make({ V: 0 })).solve(), I = held.terminals.left.current;
  assert.ok(held.converged && I > 0.5, `${I}`);
  for (const R of [1, 1e4]) {
    const s = new Device(make({ V: I * R, R })).solve();
    assert.ok(s.converged && Math.abs(s.terminals.left.V) < 1e-9, `R ${R}: ${s.converged} ${s.terminals.left.V}`);
    assert.ok(s.steps < 100, `R ${R}: ${s.steps} steps`);
  }
});
