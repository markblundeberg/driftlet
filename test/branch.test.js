// Pseudo-arclength branch continuation with linear stability marking (src/branch.js):
//   1. Fe3+/Fe2+ thin-layer cell under current drive: fold at I* = -108.1199 A/m², V* = -1.12216 V,
//      returns through I = -62.9 at V = -2.0357 V, stability flips across fold, every point on I(V).
//   2. Same Fe cell under voltage drive 0 -> -2.5 V: no fold, all points stable, exact I(V).
//   3. Hodgkin–Huxley space-clamped Hopf at 9.7–9.8 µA/cm², ω ≈ 586 rad/s, verified by fixed-step
//      kicked BDF2 runs resolving 1/ω on both sides of the Hopf.
//   4. Series-resistance terminal, named port terminal, and `solutions: true`.
//   5. Contract invariants: bit-for-bit non-destructiveness, zero clock reads at unlimited budgetMs,
//      finite budgetMs stop, stability: false, voltageBound, maxSteps, to === p0, refusals & errors.
//   6. Fuzz: short traces across random combo devices verified by neighbour re-solve to 1e-9 V.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, SolverError, units } from '../src/index.js';
import { check } from '../src/checks.js';
import { AgAgCl, aqueous, bath, build, half, hodgkinHuxley, layer, metal, ohmic, pulse, semiconductor } from '../src/kit.js';
import { caseRng } from '../bench/combos/rng.js';
import basicFamily from '../bench/combos/families/basic.js';
import echemFamily from '../bench/combos/families/echem.js';
import semiFamily from '../bench/combos/families/semi.js';
import neutralFamily from '../bench/combos/families/neutral.js';
import membraneFamily from '../bench/combos/families/membrane.js';
import mixedfilmFamily from '../bench/combos/families/mixedfilm.js';
import redoxgateFamily from '../bench/combos/families/redoxgate.js';
import wetjunctionFamily from '../bench/combos/families/wetjunction.js';

function feSpec(drive = { I: -62.90054296647188 }) {
  return {
    T: 300,
    species: [
      { name: 'Fe3+', z: 3, cRef: 1000 },
      { name: 'Fe2+', z: 2, cRef: 1000 },
      { name: 'Na+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
      { name: 'e-', z: -1 },
    ],
    materials: {
      water: {
        epsr: 78.3,
        species: {
          'Fe3+': { D: 6.04e-10, mu0: -4700 },
          'Fe2+': { D: 7.19e-10, mu0: -78900 },
          'Na+': { D: 1.334e-9, mu0: -261905 },
          'Cl-': { D: 2.032e-9, mu0: -131228 },
        },
        T: 298.15,
        mu0T: 298.15,
      },
      Pt: { conductor: { species: 'e-', conductivity: 9523809.523809524 } },
    },
    regions: [
      { material: 'Pt', length: 1e-6 },
      {
        c0: { 'Fe3+': 49.26617526505966, 'Fe2+': 2.8287313025525, 'Na+': 14.483483670977007, 'Cl-': 167.93947207126098 },
        material: 'water',
        length: 3.029892868410713e-6,
      },
      { material: 'Pt', length: 1e-6 },
    ],
    interfaces: [
      { phi: { type: 'capacitive', C: 0.28032231191173196 }, zeroCharge: 0.18351268186233938, reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 0.044465444935077036, alpha: 0.4318889715708792 }] },
      { phi: { type: 'capacitive', C: 0.28032231191173196 }, zeroCharge: 0.18351268186233938, reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 0.044465444935077036, alpha: 0.4318889715708792 }] },
    ],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk', ...drive },
    },
    grid: { hmin: 1.4284780382255565e-10, hmax: 6.391660743434045e-8, ratio: 1.1937939652428031 },
  };
}

function hhDevice(I_uAcm2, { L = 1e-6, hmin = 20e-9, hmax = 0.2e-6, g, right } = {}) {
  const T = 279.45;
  const VT = (GAS_CONSTANT * T) / FARADAY;
  const inside = { 'Na+': 50, 'K+': 400, 'Cl-': 40, 'A-': 410 };
  const out = {
    'Na+': 50 * Math.exp(0.05 / VT),
    'K+': 400 * Math.exp(-0.077 / VT),
    'Cl-': 40 * Math.exp(0.054387 / VT),
  };
  out['A-'] = out['Na+'] + out['K+'] - out['Cl-'];
  const water = aqueous(['Na+', 'K+', 'Cl-'], { epsr: 0, T });
  water.species.push({ name: 'A-', z: -1, cRef: 1000 });
  water.materials.water.species['A-'] = { D: 0.5e-9, mu0: 0 };
  const hh = hodgkinHuxley({ T, linear: true, ...(g ? { g } : {}) });
  return new Device(
    build({
      T,
      library: [water],
      stack: [
        bath(out, 0),
        layer('water', L, { c0: out }),
        { phi: { type: 'capacitive', C: 0.01 }, gates: hh.gates, species: { ...hh.species, 'A-': 'blocked' } },
        layer('water', L, { c0: inside }),
        bath(inside, right ?? { I: I_uAcm2 * 1e-2 }),
      ],
      grid: { hmin, hmax },
    }),
  );
}

// ============================================================================
// 1. Fe3+/Fe2+ thin-layer cell under current drive (fold + stability flip)
// ============================================================================

test('Fe thin-layer cell under current drive rounds the fold at I* = -108.1199 A/m², V* = -1.12216 V and flips stability', () => {
  const I0 = -62.90054296647188;
  const dev = new Device(feSpec({ V: -0.95 }));
  assert.equal(dev.solve().converged, true);
  dev.set({ contacts: { right: { I: I0 } } });
  assert.equal(dev.solve().converged, true);

  const br = dev.branch({ terminal: 'right', to: -120 });
  assert.equal(br.terminal, 'right');
  assert.equal(br.kind, 'I');
  assert.match(br.end, /returned to starting drive.*after fold/);
  assert.equal(br.events.length, 1, `expected single fold event, got ${JSON.stringify(br.events)}`);

  const fold = br.events[0];
  assert.equal(fold.kind, 'fold');
  assert.ok(Math.abs(fold.drive - -108.11992) / 108.11992 < 1e-4, `fold I*=${fold.drive}`);
  assert.ok(Math.abs(fold.terminals.right.V - -1.12216) < 1e-3, `fold V*=${fold.terminals.right.V}`);

  // Returns to I0 on the unstable branch at V = -2.0357 V.
  const last = br.points[br.points.length - 1];
  assert.ok(Math.abs(last.drive - I0) < 1e-9);
  assert.ok(Math.abs(last.terminals.right.V - -2.0357) < 1e-3, `far-side V=${last.terminals.right.V}`);
  assert.equal(last.stability.verdict, 'unstable');
  assert.equal(last.stability.im, 0);
  assert.ok(Math.abs(last.stability.rate - 150.65) < 1, `far-side rate=${last.stability.rate}`);

  // Marks before fold are stable, after fold are unstable with one real growing mode.
  const [kBefore, kAfter] = fold.between;
  for (let i = 0; i <= kBefore; i++) {
    assert.equal(br.points[i].stability.verdict, 'no growing mode found', `point ${i} before fold`);
  }
  for (let i = kAfter; i < br.points.length; i++) {
    assert.equal(br.points[i].stability.verdict, 'unstable', `point ${i} after fold`);
    assert.equal(br.points[i].stability.im, 0, `point ${i} should have real growing mode`);
    assert.ok(br.points[i].stability.rate > 0);
  }

  // Every point lies on the exact voltage-driven I(V) curve to < 1e-6 relative, and a sample of
  // points matches standalone stability() calls verbatim.
  const vProbe = new Device(feSpec({ V: br.points[0].terminals.right.V }));
  vProbe.solve();
  for (let i = 0; i < br.points.length; i++) {
    const pt = br.points[i];
    vProbe.set({ contacts: { right: { V: pt.terminals.right.V } } });
    const sv = vProbe.solve();
    assert.equal(sv.converged, true);
    const relErr = Math.abs(pt.drive - sv.terminals.right.current) / Math.abs(pt.drive);
    assert.ok(relErr < 1e-6, `point ${i} relErr=${relErr}`);
    if (i === 0 || i === kBefore || i === kAfter || i === br.points.length - 1) {
      const iProbe = new Device(feSpec({ V: pt.terminals.right.V }));
      iProbe.solve();
      iProbe.set({ contacts: { right: { I: pt.drive } } });
      iProbe.solve();
      const st = iProbe.stability();
      assert.equal(pt.stability.verdict, st.verdict);
      assert.ok(Math.abs(pt.stability.rate - st.modes[0].re) <= 1e-3 * (1 + Math.abs(st.modes[0].re)));
    }
  }
});

// ============================================================================
// 2. Same Fe cell under voltage drive 0 -> -2.5 V (no fold, all stable)
// ============================================================================

test('Fe thin-layer cell under voltage drive 0 -> -2.5 V has no fold and is stable everywhere', () => {
  const dev = new Device(feSpec({ V: 0 }));
  assert.equal(dev.solve().converged, true);
  const br = dev.branch({ terminal: 'right', to: -2.5 });
  assert.equal(br.terminal, 'right');
  assert.equal(br.kind, 'V');
  assert.match(br.end, /^reached to/);
  assert.equal(br.events.length, 0);
  assert.ok(br.points.length >= 10);
  assert.equal(br.points[0].drive, 0);
  assert.equal(br.points[br.points.length - 1].drive, -2.5);
  assert.ok(Math.abs(br.points[br.points.length - 1].terminals.right.current - -48.4741) < 1e-3);
  for (let i = 0; i < br.points.length; i++) {
    assert.equal(br.points[i].stability.verdict, 'no growing mode found', `V=${br.points[i].drive}`);
  }
});

// ============================================================================
// 3. Hodgkin–Huxley Hopf bifurcation + fixed-step kicked BDF2 verification
// ============================================================================

test('Hodgkin–Huxley space-clamped locates Hopf at 9.7–9.8 µA/cm², ω ≈ 586 rad/s, confirmed by fixed-step kicked runs', () => {
  const dev = hhDevice(0);
  assert.equal(dev.solve().converged, true);
  const br = dev.branch({ terminal: 'right', to: 0.12 });
  assert.match(br.end, /^reached to/);

  const hopfs = br.events.filter((e) => e.kind === 'hopf');
  assert.equal(hopfs.length, 1, `expected 1 Hopf event, got ${JSON.stringify(br.events)}`);
  const hopf = hopfs[0];
  const I_uA = hopf.drive * 100;
  assert.ok(I_uA >= 9.7 && I_uA <= 9.8, `Hopf I* = ${I_uA} µA/cm² outside [9.7, 9.8]`);
  assert.ok(Math.abs(hopf.omega - 586) / 586 < 0.02, `Hopf ω* = ${hopf.omega} rad/s`);

  const [kL, kR] = hopf.between;
  assert.equal(br.points[kL].stability.verdict, 'no growing mode found');
  assert.equal(br.points[kR].stability.verdict, 'unstable');
  assert.equal(br.events.filter((e) => e.kind === 'stability mismatch').length, 0);

  // Fixed-step kicked BDF2 runs (`method: 'bdf2'`, second-order; numerical damping is O((ω·dt)³)/dt
  // ≈ 0.16 /s at dt = 5e-5 s, far below |Re λ| ≈ 9.4 /s at I* ± 0.5 µA/cm², whereas plain backward
  // Euler `step(dt)` would damp by ω²dt/2 ≈ 8.6 /s and mask the growth).
  const measureEnvelope = (I_uAcm2) => {
    const d = hhDevice(I_uAcm2);
    const s0 = d.solve();
    const V0 = s0.terminals.right.V;
    const dt = 5e-5; // 1/ω ≈ 1.7 ms, so ω·dt ≈ 0.029 rad
    // Small current kick for 2 steps, then release back to I_uAcm2
    d.set({ contacts: { right: { I: I_uAcm2 * 1e-2 + 1e-4 } } });
    for (let j = 0; j < 2; j++) d.step(dt, { method: 'bdf2' });
    d.set({ contacts: { right: { I: I_uAcm2 * 1e-2 } } });
    let earlyPeak = 0, latePeak = 0;
    for (let j = 0; j < 500; j++) {
      const sol = d.step(dt, { method: 'bdf2' });
      const amp = Math.abs(sol.terminals.right.V - V0);
      if (j >= 40 && j < 200 && amp > earlyPeak) earlyPeak = amp;
      if (j >= 320 && j < 500 && amp > latePeak) latePeak = amp;
    }
    return latePeak / earlyPeak;
  };
  const ratioBelow = measureEnvelope(I_uA - 0.5);
  const ratioAbove = measureEnvelope(I_uA + 0.5);
  assert.ok(ratioBelow < 0.92, `below Hopf should decay: ratio=${ratioBelow}`);
  assert.ok(ratioAbove > 1.08, `above Hopf should grow: ratio=${ratioAbove}`);
});

// ============================================================================
// 4. Series-resistance terminal, named port terminal, and solutions: true
// ============================================================================

test('branch() traces behind a series resistance, on a named port, and with solutions: true', () => {
  const diodeR = new Device(
    build({
      T: 300,
      library: [semiconductor('Si')],
      stack: [
        ohmic(0),
        layer('Si', units.um(1), { donors: units.perCm3(1e16) }),
        layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }),
        ohmic({ V: 0, R: 1e-2 }),
      ],
      bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
      grid: { hmin: units.nm(5), hmax: units.nm(100) },
    }),
  );
  assert.equal(diodeR.solve().converged, true);
  const brR = diodeR.branch({ terminal: 'right', to: 0.5, step: 0.1, solutions: true });
  assert.equal(brR.kind, 'V');
  assert.match(brR.end, /^reached to/);
  assert.ok(brR.points.length >= 5);
  for (const pt of brR.points) {
    assert.ok(pt.solution && pt.solution.converged === true);
    assert.ok(Math.abs(pt.solution.terminals.right.V - pt.terminals.right.V) < 1e-12);
    assert.ok(Math.abs(pt.solution.terminals.right.current - pt.terminals.right.current) < 1e-12);
  }
  // Verify last point matches standalone solve at V = 0.5 behind R = 1e-2.
  diodeR.set({ contacts: { right: { V: 0.5, R: 1e-2 } } });
  const sEnd = diodeR.solve();
  assert.ok(Math.abs(sEnd.terminals.right.V - brR.points.at(-1).terminals.right.V) < 1e-10);

  // Named port terminal continuation.
  const portDev = new Device({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * FARADAY }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { phi: 'neutral' } },
    ports: [{ name: 'drain', region: 'bar', V: 0, terminal: 'Li+', species: { 'Li+': { type: 'conductance', G: 1e9 } } }],
    grid: { minCells: 20 },
  });
  assert.equal(portDev.solve().converged, true);
  const brP = portDev.branch({ terminal: 'drain', to: 0.05 });
  assert.equal(brP.terminal, 'drain');
  assert.equal(brP.kind, 'V');
  assert.match(brP.end, /^reached to/);
  assert.ok(Math.abs(brP.points.at(-1).terminals.drain.V - 0.05) < 1e-12);
});

// ============================================================================
// 5. Non-destructiveness, budgetMs, stability: false, stops, refusals & errors
// ============================================================================

test('branch() is bit-for-bit non-destructive, deterministic, honours budgetMs/maxSteps/voltageBound/stability:false, and refuses invalid calls', (t) => {
  const mkDiode = () =>
    new Device(
      build({
        T: 300,
        library: [semiconductor('Si')],
        stack: [
          ohmic(0),
          layer('Si', units.um(1), { donors: units.perCm3(1e16) }),
          layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }),
          ohmic(0),
        ],
        bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
        grid: { hmin: units.nm(10), hmax: units.nm(100) },
      }),
    );

  // Refuses unsolved device & invalid options.
  const raw = mkDiode();
  assert.throws(() => raw.branch({ to: 0.2 }), SolverError);
  raw.solve();
  assert.throws(() => raw.branch(), SolverError);
  assert.throws(() => raw.branch({ to: 0.2, bogus: 1 }), SolverError);
  assert.throws(() => raw.branch({ terminal: 'nope', to: 0.2 }), SolverError);
  assert.throws(() => raw.branch({ to: NaN }), SolverError);
  assert.throws(() => raw.branch({ to: 0.2, step: 0 }), SolverError);
  assert.throws(() => raw.branch({ to: 0.2, maxSteps: 0 }), SolverError);
  assert.throws(() => raw.branch({ to: 0.2, voltageBound: -1 }), SolverError);
  assert.throws(() => raw.branch({ to: 0.2, budgetMs: -1 }), SolverError);
  assert.throws(() => raw.branch({ to: 0.2, stability: [] }), SolverError);
  assert.throws(() => raw.branch({ to: 0.2, solutions: 1 }), SolverError);

  // to === p0 returns 1 point cleanly.
  const br0 = raw.branch({ to: 0 });
  assert.equal(br0.points.length, 1);
  assert.match(br0.end, /^reached to/);

  // Bit-for-bit non-destructiveness + unlimited budget never reads performance.now().
  const dA = mkDiode(), dB = mkDiode();
  dA.solve();
  dB.solve();
  const uBefore = Float64Array.from(dA.solver.u);
  const vBefore = Float64Array.from(dA.solver.termV);
  let clockReads = 0;
  const mock = t.mock.method(performance, 'now', () => (clockReads++, 0));
  const brNoStab = dA.branch({ terminal: 'right', to: 0.4, stability: false });
  assert.equal(clockReads, 0, 'unlimited budgetMs must never read performance.now()');
  mock.mock.restore();
  assert.equal(brNoStab.points[0].stability, undefined);
  assert.deepEqual(dA.solver.u, uBefore);
  assert.deepEqual(dA.solver.termV, vBefore);
  // Subsequent advance() matches untouched twin bit-for-bit.
  dA.set({ contacts: { right: { V: 0.2 } } });
  dB.set({ contacts: { right: { V: 0.2 } } });
  const aAdv = dA.advance(1e-7);
  const bAdv = dB.advance(1e-7);
  assert.deepEqual(aAdv.phi, bAdv.phi);
  assert.equal(aAdv.current, bAdv.current);

  // Finite budgetMs: 0 allows at least one step attempt and stops cleanly without mutating state.
  const dBud = mkDiode();
  dBud.solve();
  const uBud = Float64Array.from(dBud.solver.u);
  const brBud = dBud.branch({ terminal: 'right', to: 0.4, budgetMs: 0, stability: false });
  assert.match(brBud.end, /^budget spent/);
  assert.deepEqual(dBud.solver.u, uBud);

  // maxSteps stops cleanly.
  const brMax = dBud.branch({ terminal: 'right', to: 0.4, step: 0.05, maxSteps: 2, stability: false });
  assert.equal(brMax.points.length, 3);
  assert.match(brMax.end, /^maxSteps reached/);

  // voltageBound stops a current-driven branch that wanders in voltage.
  const dFe = new Device(feSpec({ V: 0 }));
  dFe.solve();
  dFe.set({ contacts: { right: { I: 0 } } });
  dFe.solve();
  const brVB = dFe.branch({ terminal: 'right', to: -120, voltageBound: 0.5, stability: false });
  assert.match(brVB.end, /^left voltage bound/);
  for (const pt of brVB.points) assert.ok(Math.abs(pt.terminals.right.V) <= 0.5 + 1e-9);

  // Charging current-driven terminal is refused cleanly with `end: 'refused: ...'`.
  const cap = new Device({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * FARADAY }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'dl', region: 'bar', I: 0, area: 1e6, capacitance: { C: 0.2 } }],
    grid: { hmax: 0.5e-6 },
  });
  assert.equal(cap.solve().converged, true);
  const brCap = cap.branch({ terminal: 'dl', to: 1 });
  assert.match(brCap.end, /^refused:/);
  assert.equal(brCap.points.length, 1);

  // Non-autonomous waveform drive emits a warning.
  const dWave = mkDiode();
  dWave.set({ contacts: { right: { V: pulse({ start: 1e-6, width: 1e-6, amplitude: 0.3, rise: 1e-7, base: 0 }) } } });
  assert.equal(dWave.solve().converged, true);
  const brW = dWave.branch({ terminal: 'right', to: 0.2, stability: false });
  assert.ok(brW.warnings.some((w) => w.includes('waveform')));
});

// ============================================================================
// 6. Fuzz: short traces on random combo devices re-solved from neighbour to 1e-9 V
// ============================================================================

test('Fuzz: short traces from random combo devices match neighbour solve() to 1e-9 V with no unannounced jumps', () => {
  const fams = [
    ['basic', basicFamily],
    ['echem', echemFamily],
    ['semi', semiFamily],
    ['neutral', neutralFamily],
    ['membrane', membraneFamily],
  ];
  let checkedDevices = 0;
  let checkedPoints = 0;
  for (const [name, fam] of fams) {
    for (let k = 0; k < 4; k++) {
      let def;
      try {
        const c = fam.make(caseRng(name, k));
        def = c.def ?? c;
      } catch {
        continue;
      }
      let dev;
      try {
        dev = new Device(def);
      } catch {
        continue;
      }
      const s0 = dev.solve({ budgetMs: 2000 });
      if (!s0.converged) continue;

      const dRight = dev.solver.terms[1].drive;
      const p0 = dRight.src.value ?? 0;
      const delta = dRight.kind === 'I' ? (Math.abs(p0) > 1e-6 ? 0.2 * p0 : 1e-2) : 0.15;
      const uBefore = Float64Array.from(dev.solver.u);
      const br = dev.branch({ terminal: 'right', to: p0 + delta, maxSteps: 5, budgetMs: 2000, stability: false });
      assert.ok(typeof br.end === 'string' && br.end.length > 0);
      assert.deepEqual(dev.solver.u, uBefore, `${name}_${k}: branch() mutated state`);
      checkedDevices++;

      if (br.points.length < 2) continue;
      // Re-solve sequentially from neighbour's state at each point's drive.
      const dCheck = new Device(def);
      assert.equal(dCheck.solve({ budgetMs: 2000 }).converged, true);
      for (let i = 1; i < br.points.length; i++) {
        const pt = br.points[i];
        const patch =
          dRight.kind === 'I'
            ? { contacts: { right: { I: pt.drive } } }
            : dRight.R > 0
              ? { contacts: { right: { V: pt.drive, R: dRight.R } } }
              : { contacts: { right: { V: pt.drive } } };
        dCheck.set(patch);
        const sRe = dCheck.solve({ budgetMs: 2000 });
        assert.equal(sRe.converged, true, `${name}_${k} pt ${i}: re-solve failed`);
        for (const [tn, rec] of Object.entries(pt.terminals)) {
          const dv = Math.abs(sRe.terminals[tn].V - rec.V);
          assert.ok(dv < 1e-9, `${name}_${k} pt ${i} terminal ${tn}: dV=${dv}`);
        }
        checkedPoints++;
      }
    }
  }
  assert.ok(checkedDevices >= 10, `expected >=10 fuzz devices, got ${checkedDevices}`);
  assert.ok(checkedPoints >= 20, `expected >=20 fuzz points, got ${checkedPoints}`);
});

// ============================================================================
// 7. Round-2 regression gates: far `to`, narrow double folds, span rule, true
//    max `step`, tight `budgetMs` + unlocated fold fallback, floating islands
// ============================================================================

test('HH with far `to` (1 to 50 A/m²) never skips the unstable window and finds both Hopfs', () => {
  for (const to of [1, 5, 15, 50]) {
    const dev = hhDevice(0);
    assert.equal(dev.solve().converged, true);
    const br = dev.branch({ terminal: 'right', to });
    assert.match(br.end, /^reached to/, `to=${to}: ${br.end}`);
    const hopfs = br.events.filter((e) => e.kind === 'hopf');
    const expectedCount = to < 1.54 ? 1 : 2;
    assert.equal(hopfs.length, expectedCount, `to=${to}: expected ${expectedCount} Hopf(s), got ${JSON.stringify(br.events)}`);
    assert.ok(Math.abs(hopfs[0].drive - 0.097288) < 5e-4, `to=${to}: Hopf1=${hopfs[0].drive}`);
    if (expectedCount === 2) {
      assert.ok(Math.abs(hopfs[1].drive - 1.54295) < 5e-3, `to=${to}: Hopf2=${hopfs[1].drive}`);
      assert.ok(Math.abs(hopfs[1].omega - 1062) / 1062 < 0.02, `to=${to}: Hopf2 ω=${hopfs[1].omega}`);
    }
    assert.equal(br.events.filter((e) => e.kind === 'stability mismatch').length, 0);
  }
});

test('Fe behind R just above Rcrit finds both narrow folds (including step: 2) with zero false mismatches', () => {
  for (const [R, opts] of [
    [0.01631074728741804, {}],
    [0.01631074728741804, { step: 2 }],
    [0.01628633000105963, {}],
    [1.005 * 0.01627819090560683, {}],
  ]) {
    const dev = new Device(feSpec({ V: -0.9, R }));
    assert.equal(dev.solve().converged, true);
    const br = dev.branch({ terminal: 'right', to: -3.5, ...opts });
    assert.match(br.end, /^reached to/);
    const folds = br.events.filter((e) => e.kind === 'fold');
    assert.equal(folds.length, 2, `R=${R} ${JSON.stringify(opts)}: expected 2 folds, got ${JSON.stringify(br.events)}`);
    assert.equal(br.events.filter((e) => e.kind === 'stability mismatch').length, 0, `R=${R}: false stability mismatch`);
  }
});

test('One span rule: `to` just short of a fold lands cleanly and starting right next to a fold returns to p0', () => {
  const I0 = -62.90054296647188;
  const dev = new Device(feSpec({ V: -0.95 }));
  dev.solve();
  dev.set({ contacts: { right: { I: I0 } } });
  dev.solve();

  // Reachable `to` just short of fold (-108.11992) must land on `to` without turning past it.
  for (const to of [-108.1, -108.119, -108.1199]) {
    const br = dev.branch({ terminal: 'right', to, stability: false });
    assert.match(br.end, /^reached to/, `to=${to}: ${br.end}`);
    assert.equal(br.events.length, 0, `to=${to}: unexpected events`);
    assert.equal(br.points.at(-1).drive, to);
  }

  // Starting right next to the fold must turn it and land back on I0 inside the span.
  for (const startI of [-108, -108.1, -108.119]) {
    const d = new Device(feSpec({ V: -1.0 }));
    d.solve();
    d.set({ contacts: { right: { I: startI } } });
    d.solve();
    const br = d.branch({ terminal: 'right', to: -120, stability: false, maxSteps: 60 });
    assert.match(br.end, /^returned to starting drive.*after fold/, `I0=${startI}: ${br.end}`);
    assert.equal(br.points.at(-1).drive, startI);
    assert.equal(br.events.filter((e) => e.kind === 'fold').length, 1);
  }
});

test('`step` is a true maximum |Δp| on every accepted step (including step 1 on HH to=50)', () => {
  const cases = [
    ['HH 50', () => { const d = hhDevice(0); d.solve(); return d; }, { to: 50, step: 1 }],
    ['HH 15', () => { const d = hhDevice(0); d.solve(); return d; }, { to: 15, step: 0.5 }],
    ['Fe I', () => { const d = new Device(feSpec({ V: -0.95 })); d.solve(); d.set({ contacts: { right: { I: -62.9 } } }); d.solve(); return d; }, { to: -120, step: 2 }],
  ];
  for (const [label, mk, opts] of cases) {
    const d = mk();
    const br = d.branch({ terminal: 'right', stability: false, ...opts });
    for (let i = 1; i < br.points.length; i++) {
      const dp = Math.abs(br.points[i].drive - br.points[i - 1].drive);
      assert.ok(dp <= opts.step * (1 + 1e-6), `${label} step ${i}: |Δp|=${dp} > step=${opts.step}`);
    }
  }
});

test('budgetMs stays bounded inside stability() and emits an unlocated fold event if cut short during fold location', (t) => {
  const d = new Device(feSpec({ V: -0.95 }));
  d.solve();
  d.set({ contacts: { right: { I: -62.90054296647188 } } });
  d.solve();

  // Tight budget with stability: true must stop well under one full 90 ms stability() call.
  const t0 = performance.now();
  const brFast = d.branch({ terminal: 'right', to: -120, budgetMs: 20, stability: true });
  const elapsed = performance.now() - t0;
  assert.match(brFast.end, /^budget spent/);
  assert.ok(elapsed < 60, `budgetMs=20 with stability overran to ${elapsed.toFixed(1)} ms`);

  // Simulated budget expiry right inside fold location: fold event is still emitted with located: false
  // and valid `between` indices in `points`.
  let ticks = 0;
  let clock = 0;
  // Advance clock sharply once fold bisection begins (after ~35 Newton/check reads).
  const mock = t.mock.method(performance, 'now', () => {
    ticks++;
    if (ticks > 35) clock += 100;
    return clock;
  });
  try {
    // Run with stability: false and budgetMs: 50; adjust threshold to hit mid-fold.
    let sawUnlocated = false;
    for (const limit of [25, 35, 45, 55, 65, 75, 85]) {
      ticks = 0;
      clock = 0;
      mock.mock.mockImplementation(() => {
        ticks++;
        if (ticks > limit) clock = 100;
        return clock;
      });
      const br = d.branch({ terminal: 'right', to: -120, budgetMs: 50, stability: false });
      const unloc = br.events.find((e) => e.kind === 'fold' && e.located === false);
      if (unloc) {
        assert.ok(unloc.between[0] >= 0 && unloc.between[1] < br.points.length, `invalid between ${unloc.between} for ${br.points.length} pts`);
        assert.match(unloc.note ?? '', /budgetMs/);
        sawUnlocated = true;
        break;
      }
    }
    assert.ok(sawUnlocated, 'expected an unlocated fold event when budget expired mid-fold');
  } finally {
    mock.mock.restore();
  }
});

test('Floating-island devices (echem_4, semi_1, membrane_2) trace cleanly to target without losing Kirchhoff balance', () => {
  for (const [name, fam, k, side, delta] of [
    ['echem_4', echemFamily, 4, 'left', 1.5],
    ['semi_1', semiFamily, 1, 'left', 0.3],
    ['membrane_2', membraneFamily, 2, 'left', 1.5],
  ]) {
    const def = fam.make(caseRng(name.split('_')[0], k)).def;
    const dev = new Device(structuredClone(def));
    assert.equal(dev.solve({ budgetMs: 3000 }).converged, true);
    const tIdx = dev.solver.terms.findIndex((tm) => tm.name === side);
    const p0 = dev.solver.terms[tIdx].drive.src.value ?? 0;
    const to = p0 + delta;
    const br = dev.branch({ terminal: side, to, maxSteps: 40, stability: false });
    assert.match(br.end, /^reached to/, `${name}: ${br.end}`);
    assert.ok(br.points.length >= 2);
    for (let i = 0; i < br.points.length; i++) {
      let sum = 0, mx = 0;
      for (const rec of Object.values(br.points[i].terminals)) {
        sum += rec.current;
        mx = Math.max(mx, Math.abs(rec.current));
      }
      if (mx > 1e-18 && name !== 'semi_1') {
        assert.ok(Math.abs(sum) <= 1e-6 * mx, `${name} pt ${i}: Kirchhoff |ΣI|/max=${(Math.abs(sum) / mx).toExponential(2)}`);
      }
    }
    // Verify endpoint matches plain solve().
    const chk = new Device(structuredClone(def));
    chk.solve({ budgetMs: 3000 });
    chk.set({ contacts: { [side]: { V: to } } });
    const sEnd = chk.solve({ budgetMs: 3000 });
    assert.equal(sEnd.converged, true);
    for (const [tn, rec] of Object.entries(br.points.at(-1).terminals)) {
      assert.ok(Math.abs(sEnd.terminals[tn].V - rec.V) < 1e-8, `${name} ${tn} dV=${Math.abs(sEnd.terminals[tn].V - rec.V)}`);
    }
  }
});

// Round-3 regression tests (every one fails on Round 2 cc4f8c1):

test('Downward HH traces from high current (stock and narrowed g.Na=850) follow secondary modes and locate both Hopfs', () => {
  // I0 in A/m² = I_uAcm2 / 100.
  const cases = [
    ['stock 120->-5', null, 12000, -5, [1.5429, 0.09729]],
    ['stock 3->-3 (large ω shift)', null, 300, -3, [1.5429, 0.09729]],
    ['g.Na=850 120->-5', { Na: 850 }, 12000, -5, [0.7683, 0.3799]],
    ['g.Na=850 5->-5', { Na: 850 }, 500, -5, [0.7683, 0.3799]],
  ];
  for (const [label, g, I0uA, to, [hiRef, loRef]] of cases) {
    const dev = hhDevice(I0uA, g ? { g } : {});
    assert.equal(dev.solve().converged, true);
    const br = dev.branch({ terminal: 'right', to });
    assert.match(br.end, /^reached to/, `${label}: ${br.end}`);
    const hopfs = br.events.filter((e) => e.kind === 'hopf');
    const mism = br.events.filter((e) => e.kind === 'stability mismatch');
    assert.equal(mism.length, 0, `${label}: unexpected stability mismatch ${JSON.stringify(mism)}`);
    assert.equal(hopfs.length, 2, `${label}: expected 2 Hopfs, got ${JSON.stringify(br.events)}`);
    assert.ok(Math.abs(hopfs[0].drive - hiRef) < 5e-3, `${label}: upper Hopf=${hopfs[0].drive} vs ${hiRef}`);
    assert.ok(Math.abs(hopfs[1].drive - loRef) < 5e-3, `${label}: lower Hopf=${hopfs[1].drive} vs ${loRef}`);
  }
});

test('redoxgate_0 left 0->1.5 emits zero false folds and stops honestly if Jacobian loses rank', () => {
  const def = redoxgateFamily.make(caseRng('redoxgate', 0)).def;
  const dev = new Device(structuredClone(def));
  assert.equal(dev.solve({ budgetMs: 3000 }).converged, true);
  const br = dev.branch({ terminal: 'left', to: 1.5, maxSteps: 120 });
  assert.match(br.end, /^(reached to|continuation stopped:.*solve\(\) may continue from here)/, `redoxgate_0: ${br.end}`);
  assert.equal(br.events.filter((e) => e.kind === 'fold').length, 0, `unexpected false folds: ${JSON.stringify(br.events)}`);
});

test('Depleting combo devices (mixedfilm_0/1, semi_5) resolve smoothly or stop honestly without false folds', () => {
  for (const [label, fam, k, term, to] of [
    ['mixedfilm_0', mixedfilmFamily, 0, 'right', -1.5],
    ['mixedfilm_1', mixedfilmFamily, 1, 'right', -1.5],
  ]) {
    const def = fam.make(caseRng('mixedfilm', k)).def;
    const dev = new Device(structuredClone(def));
    assert.equal(dev.solve({ budgetMs: 3000 }).converged, true);
    const br = dev.branch({ terminal: term, to, maxSteps: 120, budgetMs: 20000 });
    assert.match(br.end, /^reached to/, `${label}: ${br.end}`);
    assert.equal(br.events.filter((e) => e.kind === 'fold').length, 0, `${label}: unexpected false fold`);
  }
  const semiDef = semiFamily.make(caseRng('semi', 5)).def;
  const semiDev = new Device(structuredClone(semiDef));
  assert.equal(semiDev.solve({ budgetMs: 3000 }).converged, true);
  const sbr = semiDev.branch({ terminal: 'left', to: -1.5, maxSteps: 40, budgetMs: 20000 });
  assert.match(sbr.end, /^reached to/, `semi_5: ${sbr.end}`);
  assert.equal(sbr.events.filter((e) => e.kind === 'fold').length, 0);
});

test('points always holds the start point even when budgetMs expires during start-point stability()', () => {
  const dev = hhDevice(0);
  assert.equal(dev.solve().converged, true);
  const br = dev.branch({ terminal: 'right', to: 5, budgetMs: 0, stability: true });
  assert.match(br.end, /^budget spent/);
  assert.equal(br.points.length, 1, `expected start point kept, got ${br.points.length}`);
  assert.equal(br.points[0].drive, 0);
  assert.equal(br.points[0].stability.verdict, 'undecided');
  assert.equal(br.points[0].stability.checked, false);
  assert.match(br.points[0].stability.reason, /budgetMs/);
});

// Round-4 regression tests (every one fails on Round 3 5e92bc1):

test('Bistable HH membrane traces L→M→H in monotone V order with zero reversals, both folds, and Hopf (plus 1200/10 stationary mode)', () => {
  const g = { Na: 3000, K: 60, leak: 3 };
  // Voltage-clamp reference I∞(V) to verify every traced point lies on the true curve.
  const ref = hhDevice(-300, { g });
  assert.equal(ref.solve().converged, true);
  const Iof = (V) => {
    ref.set({ contacts: { right: { V } } });
    const s = ref.solve();
    assert.equal(s.converged, true);
    return s.terminals.right.current;
  };
  for (const stab of [false, true]) {
    const dev = hhDevice(-300, { g });
    assert.equal(dev.solve().converged, true);
    const br = dev.branch({ terminal: 'right', to: 1, stability: stab, maxSteps: 400 });
    assert.match(br.end, /^reached to/, `3000/60 stab=${stab}: ${br.end}`);
    const P = br.points;
    for (let i = 2; i < P.length; i++) {
      const a = P[i - 1].terminals.right.V - P[i - 2].terminals.right.V;
      const b = P[i].terminals.right.V - P[i - 1].terminals.right.V;
      assert.ok(a * b >= 0, `3000/60 stab=${stab}: V_right reversed at point ${i - 1} (${P[i - 2].terminals.right.V} -> ${P[i - 1].terminals.right.V} -> ${P[i].terminals.right.V})`);
    }
    if (!stab) {
      let worst = 0;
      for (const p of P) worst = Math.max(worst, Math.abs(Iof(p.terminals.right.V) - p.drive));
      assert.ok(worst < 1e-11, `off-manifold residual ${worst}`);
    }
    const folds = br.events.filter((e) => e.kind === 'fold');
    assert.equal(folds.length, 2, `3000/60 stab=${stab}: expected 2 folds, got ${JSON.stringify(br.events)}`);
    assert.ok(Math.abs(folds[0].drive - -0.050402) < 1e-4, `fold 0: ${folds[0].drive}`);
    assert.ok(Math.abs(folds[1].drive - -1.215051) < 1e-4, `fold 1: ${folds[1].drive}`);
    if (stab) {
      const hopfs = br.events.filter((e) => e.kind === 'hopf');
      assert.equal(hopfs.length, 1, `3000/60 stab=true: expected 1 Hopf, got ${JSON.stringify(br.events)}`);
      assert.ok(Math.abs(hopfs[0].drive - -0.4722) < 5e-3, `Hopf: ${hopfs[0].drive}`);
    }
  }
  // 3000/20 from -40 -> 20 A/m² exercises dynamic metric rebalancing across wide slope change.
  const dev30 = hhDevice(-4000, { g: { Na: 3000, K: 20, leak: 3 } });
  assert.equal(dev30.solve().converged, true);
  const br30 = dev30.branch({ terminal: 'right', to: 20, stability: false, maxSteps: 400 });
  assert.match(br30.end, /^reached to/, `3000/20 -40->20: ${br30.end}`);
  for (let i = 2; i < br30.points.length; i++) {
    const a = br30.points[i - 1].terminals.right.V - br30.points[i - 2].terminals.right.V;
    const b = br30.points[i].terminals.right.V - br30.points[i - 1].terminals.right.V;
    assert.ok(a * b >= 0, `3000/20 -40->20: V_right reversed at point ${i - 1}`);
  }
  assert.equal(br30.events.filter((e) => e.kind === 'fold').length, 2);
  // 1200/10 downward from 20 -> -20 A/m² exercises stationary degenerate mode at I = 20.
  const dev12 = hhDevice(2000, { g: { Na: 1200, K: 10, leak: 3 } });
  assert.equal(dev12.solve().converged, true);
  const br12 = dev12.branch({ terminal: 'right', to: -20, maxSteps: 400 });
  assert.match(br12.end, /^reached to/, `1200/10 20->-20: ${br12.end}`);
  assert.equal(br12.events.filter((e) => e.kind === 'fold').length, 2);
  assert.equal(br12.events.filter((e) => e.kind === 'hopf').length, 1);
});

test('Nearer-end mode rule catches held-out narrow HH windows (g.Na=835, 1100) with zero false events at g.Na=820', () => {
  // g.Na = 835 has a narrow unstable window (0.445–0.67 A/m², peak Re λ ≈ +6.9 /s) skipped in Round 3.
  for (const [label, gNa, I0uA, to] of [
    ['g.Na=835 6.5->0.05', 835, 650, 0.05],
    ['g.Na=835 65->-0.5', 835, 6500, -0.5],
    ['g.Na=1100 300->0.05', 1100, 30000, 0.05],
  ]) {
    const dev = hhDevice(I0uA, { g: { Na: gNa } });
    assert.equal(dev.solve().converged, true);
    const br = dev.branch({ terminal: 'right', to });
    assert.match(br.end, /^reached to/, `${label}: ${br.end}`);
    const hopfs = br.events.filter((e) => e.kind === 'hopf');
    assert.equal(hopfs.length, 2, `${label}: expected 2 Hopfs, got ${JSON.stringify(br.events)}`);
    assert.ok(br.points.some((p) => p.stability.verdict === 'unstable'), `${label}: expected unstable points`);
  }
  // Negative control: g.Na = 820 has no unstable window and must report 0 events.
  const dev820 = hhDevice(650, { g: { Na: 820 } });
  assert.equal(dev820.solve().converged, true);
  const br820 = dev820.branch({ terminal: 'right', to: 0.05 });
  assert.match(br820.end, /^reached to/);
  assert.equal(br820.events.length, 0, `g.Na=820 false events: ${JSON.stringify(br820.events)}`);
  assert.ok(br820.points.every((p) => p.stability.verdict === 'no growing mode found'));
});

test('Scale-free drive tolerances: amperes drives (A = 1e-6..1e-10) on Fe, HH, bistable membrane, and far to = -1e8', () => {
  const hhArea = (I, A, g) => {
    const T = 279.45;
    const VT = (GAS_CONSTANT * T) / FARADAY;
    const L = 1e-6;
    const inside = { 'Na+': 50, 'K+': 400, 'Cl-': 40, 'A-': 410 };
    const out = {
      'Na+': 50 * Math.exp(0.05 / VT),
      'K+': 400 * Math.exp(-0.077 / VT),
      'Cl-': 40 * Math.exp(0.054387 / VT),
    };
    out['A-'] = out['Na+'] + out['K+'] - out['Cl-'];
    const water = aqueous(['Na+', 'K+', 'Cl-'], { epsr: 0, T });
    water.species.push({ name: 'A-', z: -1, cRef: 1000 });
    water.materials.water.species['A-'] = { D: 0.5e-9, mu0: 0 };
    const h = hodgkinHuxley({ T, linear: true, ...(g ? { g } : {}) });
    const spec = build({
      T,
      library: [water],
      stack: [
        bath(out, 0),
        layer('water', L, { c0: out }),
        { phi: { type: 'capacitive', C: 0.01 }, gates: h.gates, species: { ...h.species, 'A-': 'blocked' } },
        layer('water', L, { c0: inside }),
        bath(inside, { I: I * A }),
      ],
      grid: { hmin: 20e-9, hmax: 0.2e-6 },
    });
    spec.geometry = { area: { x: [0, 2 * L], values: [A, A] } };
    return new Device(spec);
  };

  // 1. Bistable membrane at A = 1e-9 m² (both stab=false and stab=true): zero branch hop, both folds located.
  const gBi = { Na: 3000, K: 60, leak: 3 };
  for (const stab of [false, true]) {
    const d = hhArea(-3, 1e-9, gBi);
    assert.equal(d.solve().converged, true);
    const br = d.branch({ to: 1e-9, stability: stab, maxSteps: 400 });
    assert.match(br.end, /^reached to/, `bistable 1e-9 stab=${stab}: ${br.end}`);
    const folds = br.events.filter((e) => e.kind === 'fold');
    assert.equal(folds.length, 2, `bistable 1e-9 stab=${stab}: expected 2 folds, got ${JSON.stringify(br.events)}`);
    assert.notEqual(folds[0].located, false);
    assert.notEqual(folds[1].located, false);
    assert.ok(Math.abs(folds[0].drive / 1e-9 - -0.050402) < 1e-4);
    assert.ok(Math.abs(folds[1].drive / 1e-9 - -1.215051) < 1e-4);
  }

  // 2. HH at A = 1e-6 (duplicate multi-shift mode dedup at I=0) and A = 1e-9 (exact Hopf location, not 21% off).
  for (const A of [1e-6, 1e-9]) {
    const d = hhArea(0, A);
    assert.equal(d.solve().converged, true);
    const br = d.branch({ to: 5 * A });
    assert.match(br.end, /^reached to/, `HH A=${A}: ${br.end}`);
    const hopfs = br.events.filter((e) => e.kind === 'hopf');
    assert.equal(hopfs.length, 2, `HH A=${A}: ${JSON.stringify(br.events)}`);
    assert.notEqual(hopfs[0].located, false);
    assert.notEqual(hopfs[1].located, false);
    assert.ok(Math.abs(hopfs[0].drive / A - 0.0972858) < 1e-4, `HH A=${A} lower Hopf: ${hopfs[0].drive / A}`);
    assert.ok(Math.abs(hopfs[1].drive / A - 1.54295) < 5e-4, `HH A=${A} upper Hopf: ${hopfs[1].drive / A}`);
  }

  // 3. Fe cell at A = 1e-8, 1e-10 and planar with far to = -1e8: fold located to 6+ digits.
  for (const A of [1e-8, 1e-10]) {
    const spec = feSpec({ V: -0.95 });
    spec.geometry = { area: { x: [0, 5.029892868410713e-6], values: [A, A] } };
    const d = new Device(spec);
    assert.equal(d.solve().converged, true);
    d.set({ contacts: { right: { I: -62.90054296647188 * A } } });
    assert.equal(d.solve().converged, true);
    const br = d.branch({ terminal: 'right', to: -120 * A });
    const folds = br.events.filter((e) => e.kind === 'fold');
    assert.equal(folds.length, 1, `Fe A=${A}: ${JSON.stringify(br.events)}`);
    assert.notEqual(folds[0].located, false, `Fe A=${A} unlocated: ${JSON.stringify(folds[0])}`);
    assert.ok(Math.abs(folds[0].drive / A - -108.11992) < 1e-3, `Fe A=${A}: ${folds[0].drive / A}`);
  }
  {
    const d = new Device(feSpec({ V: -0.95 }));
    assert.equal(d.solve().converged, true);
    d.set({ contacts: { right: { I: -62.90054296647188 } } });
    assert.equal(d.solve().converged, true);
    const br = d.branch({ terminal: 'right', to: -1e8, stability: false });
    const folds = br.events.filter((e) => e.kind === 'fold');
    assert.equal(folds.length, 1, `Fe far to=-1e8: ${JSON.stringify(br.events)}`);
    assert.notEqual(folds[0].located, false);
    assert.ok(Math.abs(folds[0].drive - -108.11992) < 1e-4, `Fe far to=-1e8: ${folds[0].drive}`);
  }
});

test('Stability flip with no fold/Hopf emits stability mismatch and continues (echem_25); budget end reports live elapsed ms', () => {
  // echem_25 right: real mode flips sign through infinity near -1.2367 V with no fold.
  const def25 = echemFamily.make(caseRng('echem', 25)).def;
  const d25 = new Device(structuredClone(def25));
  assert.equal(d25.solve({ budgetMs: 3000 }).converged, true);
  const p0 = d25.solver.terms[1].drive.src.value;
  const br25 = d25.branch({ terminal: 'right', to: p0 - 1.5, maxSteps: 200 });
  assert.match(br25.end, /^reached to/, `echem_25: ${br25.end}`);
  assert.ok(br25.events.some((e) => e.kind === 'stability mismatch'), `echem_25 expected stability mismatch: ${JSON.stringify(br25.events)}`);
  assert.ok(br25.points.some((p) => p.stability.verdict === 'unstable'), 'echem_25 should expose the unstable side');

  // echem_23 budgetMs: reported elapsed time in `end` reflects actual wall time (not stale cached time).
  const def23 = echemFamily.make(caseRng('echem', 23)).def;
  const d23 = new Device(structuredClone(def23));
  assert.equal(d23.solve({ budgetMs: 3000 }).converged, true);
  const p23 = d23.solver.terms[1].drive.src.value;
  const budgetMs = 400;
  const t0 = performance.now();
  const br23 = d23.branch({ terminal: 'right', to: p23 - 1.5, budgetMs });
  const wall = performance.now() - t0;
  const m = br23.end.match(/^budget spent after ([\d.e+-]+) ms \(budgetMs: 400\)$/);
  assert.ok(m, `echem_23 end: ${br23.end}`);
  const reported = Number(m[1]);
  assert.ok(reported >= budgetMs * 0.9 && reported <= Math.max(wall * 1.25, budgetMs + 200),
    `stale budget elapsed: reported=${reported} ms, wall=${wall.toFixed(0)} ms, budgetMs=${budgetMs}`);
});

// Round-5 regression & oracle tests (every new test fails on Round 4 016386a):

test('Reference event-list diff without branch(): Fe, HH, bistable membrane (I & behind R), and round-off plateau', () => {
  // Independent oracle: builds the fold & Hopf event list from voltage-clamped I∞(V) + stability()
  // alone (zero branch() calls), then diffs branch().events against it across spans & step options.
  const hhVClamp = (g) => {
    const d = hhDevice(0, g ? { g } : {});
    d.set({ contacts: { right: { V: -0.065 } } });
    assert.equal(d.solve().converged, true);
    return d;
  };
  const goldenExtremum = (f, a, b, isMax) => {
    const phi = 0.5 * (Math.sqrt(5) - 1);
    let x1 = b - phi * (b - a), x2 = a + phi * (b - a);
    let f1 = f(x1), f2 = f(x2);
    for (let it = 0; it < 38; it++) {
      if (isMax ? f1 < f2 : f1 > f2) { a = x1; x1 = x2; f1 = f2; x2 = a + phi * (b - a); f2 = f(x2); }
      else { b = x2; x2 = x1; f2 = f1; x1 = b - phi * (b - a); f1 = f(x1); }
    }
    const v = 0.5 * (a + b);
    return { v, val: f(v) };
  };
  const findFoldsFromVClamp = (driveOfV, vMin, vMax, nGrid) => {
    const vs = [], ps = [];
    for (let i = 0; i <= nGrid; i++) {
      const v = vMin + (vMax - vMin) * (i / nGrid);
      vs.push(v);
      ps.push(driveOfV(v));
    }
    const folds = [];
    for (let i = 1; i < nGrid; i++) {
      if ((ps[i] > ps[i - 1] && ps[i] > ps[i + 1]) || (ps[i] < ps[i - 1] && ps[i] < ps[i + 1])) {
        const isMax = ps[i] > ps[i - 1];
        const ex = goldenExtremum(driveOfV, vs[i - 1], vs[i + 1], isMax);
        folds.push({ kind: 'fold', v: ex.v, drive: ex.val });
      }
    }
    return folds;
  };
  const findHopfsFromVClamp = (vDev, driveAtV, vMin, vMax, nGrid) => {
    const hopfRe = (v) => {
      vDev.set({ contacts: { right: { V: v, R: undefined } } });
      assert.equal(vDev.solve().converged, true);
      const I = vDev.solver._currentInto(1);
      const drive = driveAtV(v, I);
      const snap = vDev.solver._snapshot();
      vDev.set({ contacts: { right: drive } });
      vDev.solver._restore(snap);
      vDev.solver.termV[1] = v;
      vDev.solver.atSteady = true;
      const st = vDev.stability();
      vDev.set({ contacts: { right: { V: v, R: undefined } } });
      vDev.solver._restore(snap);
      vDev.solver.atSteady = true;
      const osc = st.modes.find((m) => m.im > 10 * m.error);
      return { osc, driveVal: drive.I ?? drive.V };
    };
    const hopfs = [];
    let prev = null;
    for (let i = 0; i <= nGrid; i++) {
      const v = vMin + (vMax - vMin) * (i / nGrid);
      const cur = { v, ...hopfRe(v) };
      if (prev && prev.osc && cur.osc && prev.osc.re * cur.osc.re < 0) {
        let a = prev.v, b = cur.v, ra = prev.osc.re;
        let lastDrive = cur.driveVal;
        for (let it = 0; it < 20; it++) {
          const mid = 0.5 * (a + b);
          const m = hopfRe(mid);
          lastDrive = m.driveVal;
          if (!m.osc) break;
          if (m.osc.re * ra > 0) { a = mid; ra = m.osc.re; } else b = mid;
        }
        hopfs.push({ kind: 'hopf', v: 0.5 * (a + b), drive: lastDrive });
      }
      prev = cur;
    }
    return hopfs;
  };
  const assertEventDiff = (label, br, refEvents, relTol = 2e-4) => {
    const mism = br.events.filter((e) => e.kind === 'stability mismatch');
    assert.equal(mism.length, 0, `${label}: unexpected stability mismatch ${JSON.stringify(mism)}`);
    assert.equal(br.events.length, refEvents.length,
      `${label}: event count ${br.events.length} vs ref ${refEvents.length}\n  got: ${JSON.stringify(br.events.map((e) => [e.kind, e.drive]))}\n  ref: ${JSON.stringify(refEvents.map((e) => [e.kind, e.drive]))}`);
    for (let i = 0; i < refEvents.length; i++) {
      assert.equal(br.events[i].kind, refEvents[i].kind, `${label} ev[${i}] kind`);
      const refD = refEvents[i].drive;
      let minGap = Infinity;
      for (let j = 0; j < refEvents.length; j++) {
        if (j === i) continue;
        const g = Math.abs(refEvents[j].drive - refD);
        if (g > 0 && g < minGap) minGap = g;
      }
      const cap = Math.min(
        relTol * Math.max(Math.abs(refD), 1e-12),
        minGap < Infinity ? 0.25 * minGap : Infinity,
      );
      assert.ok(Math.abs(br.events[i].drive - refD) <= cap,
        `${label} ev[${i}] (${refEvents[i].kind}): got ${br.events[i].drive} vs ref ${refD} (tol ${cap.toExponential(2)})`);
    }
  };

  // 1. Fe cell (current drive): reference extremum from V-clamp I(V), checked downward & upward (no crawl!).
  {
    const feV = new Device(feSpec({ V: 0 }));
    assert.equal(feV.solve().converged, true);
    const IofV = (V) => { feV.set({ contacts: { right: { V } } }); assert.equal(feV.solve().converged, true); return feV.solver._currentInto(1); };
    const feFolds = findFoldsFromVClamp(IofV, -1.6, 1.6, 64);
    assert.equal(feFolds.length, 2);
    const negFold = feFolds.find((f) => f.drive < 0);
    const posFold = feFolds.find((f) => f.drive > 0);
    // Upward Fe 0 -> +300: must find positive fold cleanly without unstable-branch crawl (was 201 pts in R4).
    const dUp = new Device(feSpec({ I: 0 }));
    assert.equal(dUp.solve().converged, true);
    const brUp = dUp.branch({ terminal: 'right', to: 300, maxSteps: 100 });
    assert.ok(!/^maxSteps/.test(brUp.end), `Fe 0->+300: ${brUp.end}`);
    assert.ok(brUp.points.length < 60, `Fe 0->+300 crawled: ${brUp.points.length} pts`);
    assertEventDiff('Fe 0->+300', brUp, [posFold]);
    assertEventDiff('Fe -62.9->-120', (() => {
      const d = new Device(feSpec({ V: -0.95 }));
      assert.equal(d.solve().converged, true);
      d.set({ contacts: { right: { I: -62.90054296647188 } } });
      assert.equal(d.solve().converged, true);
      return d.branch({ terminal: 'right', to: -120 });
    })(), [negFold]);
  }

  // 2. HH stock (0->5, 5->0, and step: 0.5): diff against V-clamp Hopf reference.
  {
    const hhV = hhVClamp(null);
    const refH = findHopfsFromVClamp(hhV, (_v, I) => ({ I }), -0.065, -0.04, 40).sort((a, b) => a.v - b.v);
    assert.equal(refH.length, 2);
    for (const [label, I0uA, to, extra, rev] of [
      ['HH 0->5', 0, 5, {}, false],
      ['HH 5->0', 500, 0, {}, true],
      ['HH 0->5 step=0.5', 0, 5, { step: 0.5 }, false],
    ]) {
      const d = hhDevice(I0uA);
      assert.equal(d.solve().converged, true);
      const br = d.branch({ terminal: 'right', to, ...extra });
      assert.match(br.end, /^reached to/, `${label}: ${br.end}`);
      assertEventDiff(label, br, rev ? [...refH].reverse() : refH);
    }
  }

  // 3. Bistable membrane behind R (g = 2500/45/3, R = 1.02 Rc and 2 Rc): lower-branch Hopf sits right
  //    before the fold; Morse-index & sawHidden refinement catch all 4 events with 0 mismatches.
  {
    const g = { Na: 2500, K: 45, leak: 3 };
    const vDev = hhVClamp(g);
    const IofV = (V) => { vDev.set({ contacts: { right: { V } } }); assert.equal(vDev.solve().converged, true); return vDev.solver._currentInto(1); };
    // Estimate Rcrit from steepest negative dI/dV on I∞(V).
    let minSlope = 0;
    let vPrev = -0.075, iPrev = IofV(vPrev);
    for (let k = 1; k <= 80; k++) {
      const v = -0.075 + 0.06 * (k / 80);
      const curI = IofV(v);
      minSlope = Math.min(minSlope, (curI - iPrev) / (v - vPrev));
      vPrev = v; iPrev = curI;
    }
    const Rcrit = -1 / minSlope;
    for (const [mult, runs] of [
      [1.02, [[-0.606, 0.606, false], [0.606, -0.606, true], [0.16, -0.24, true]]],
      [2.0, [[-0.36, 0.24, false], [0.24, -0.36, true]]],
    ]) {
      const R = mult * Rcrit;
      const EofV = (V) => V + R * IofV(V);
      const refFolds = findFoldsFromVClamp(EofV, -0.07, -0.02, 120);
      const refHopfs = findHopfsFromVClamp(vDev, (v, I) => ({ V: v + R * I, R }), -0.07, -0.02, 70);
      const refAll = [...refFolds, ...refHopfs].sort((a, b) => a.v - b.v);
      assert.equal(refFolds.length, 2, `R=${mult}Rc ref folds`);
      assert.equal(refHopfs.length, 2, `R=${mult}Rc ref hopfs`);
      for (const [E0, to, rev] of runs) {
        const d = hhDevice(0, { g });
        d.set({ contacts: { right: { V: E0, R } } });
        assert.equal(d.solve().converged, true);
        const br = d.branch({ terminal: 'right', to, maxSteps: 300 });
        assert.match(br.end, /^reached to/, `2500/45/3 R=${mult}Rc ${E0}->${to}: ${br.end}`);
        assertEventDiff(`2500/45/3 R=${mult}Rc ${E0}->${to}`, br, rev ? [...refAll].reverse() : refAll);
      }
    }
  }

  // 4. Plain I-drive lower-branch Hopf 1.25e-4 A/m² before fold (g = 4000/100/5): all 4 events, 0 mismatches.
  {
    const g = { Na: 4000, K: 100, leak: 5 };
    const vDev = hhVClamp(g);
    const IofV = (V) => { vDev.set({ contacts: { right: { V } } }); assert.equal(vDev.solve().converged, true); return vDev.solver._currentInto(1); };
    const refFolds = findFoldsFromVClamp(IofV, -0.07, -0.01, 100);
    const refHopfs = [
      ...findHopfsFromVClamp(vDev, (_v, I) => ({ I }), refFolds[0].v - 0.002, refFolds[0].v, 16),
      ...findHopfsFromVClamp(vDev, (_v, I) => ({ I }), -0.035, -0.02, 24),
    ];
    const refAll = [...refFolds, ...refHopfs].sort((a, b) => a.v - b.v);
    assert.equal(refAll.length, 4, `4000/100/5 ref events: ${JSON.stringify(refAll)}`);
    const d = hhDevice(-400, { g });
    assert.equal(d.solve().converged, true);
    const br = d.branch({ terminal: 'right', to: 2, maxSteps: 300 });
    assert.match(br.end, /^reached to/, `4000/100/5 -4->2: ${br.end}`);
    assertEventDiff('4000/100/5 -4->2', br, refAll);
  }

  // 5. Round-off plateau vs genuine nA fold: echem_32 ref walks past the plateau to voltageBound
  //    with 0 folds (no invented plateau stop); echem_29 right still locates its genuine nA fold.
  {
    const def32 = echemFamily.make(caseRng('echem', 32)).def;
    for (const to of [0.05, -0.05]) {
      const d = new Device(structuredClone(def32));
      assert.equal(d.solve({ budgetMs: 3000 }).converged, true);
      const br = d.branch({ terminal: 'ref', to, maxSteps: 200 });
      assert.equal(br.events.filter((e) => e.kind === 'fold').length, 0, `echem_32 ref ->${to}: false round-off folds ${JSON.stringify(br.events)}`);
      assert.match(br.end, /^left voltage bound /, `echem_32 ref ->${to}: ${br.end}`);
    }
    const def29 = echemFamily.make(caseRng('echem', 29)).def;
    const d29 = new Device(structuredClone(def29));
    assert.equal(d29.solve({ budgetMs: 3000 }).converged, true);
    const br29 = d29.branch({ terminal: 'right', to: -0.05, maxSteps: 100 });
    const f29 = br29.events.filter((e) => e.kind === 'fold');
    assert.equal(f29.length, 1, `echem_29 right: expected 1 genuine nA fold, got ${JSON.stringify(br29.events)}`);
    assert.ok(Math.abs(f29[0].drive - -9.4859e-9) < 1e-10, `echem_29 fold drive=${f29[0].drive}`);
  }
});

test('Metric rebalancing, insensitive terminals, lost-mode GMRES, landing guard, and held-out reach', () => {
  // 1. Rebalance metric (pins rebalanceMetric): 2500/45/3 -4 -> 2 stab=false has 0 V-reversals & 2 folds.
  {
    const d = hhDevice(-400, { g: { Na: 2500, K: 45, leak: 3 } });
    assert.equal(d.solve().converged, true);
    const br = d.branch({ terminal: 'right', to: 2, stability: false, maxSteps: 400 });
    assert.match(br.end, /^reached to/, `2500/45/3 -4->2: ${br.end}`);
    for (let i = 2; i < br.points.length; i++) {
      const a = br.points[i - 1].terminals.right.V - br.points[i - 2].terminals.right.V;
      const b = br.points[i].terminals.right.V - br.points[i - 1].terminals.right.V;
      assert.ok(a * b >= 0, `2500/45/3 V_right reversed at pt ${i - 1}`);
    }
    assert.equal(br.events.filter((e) => e.kind === 'fold').length, 2);
  }

  // 2. Insensitive terminals (echem_30 link, echem_31 cap): must not jump in 1 step or refuse `step: 0.1`.
  for (const [k, term] of [[30, 'link'], [31, 'cap']]) {
    const def = echemFamily.make(caseRng('echem', k)).def;
    const d1 = new Device(structuredClone(def));
    assert.equal(d1.solve({ budgetMs: 3000 }).converged, true);
    const tIdx = d1.solver.terms.findIndex((tm) => tm.name === term);
    const to = (d1.solver.terms[tIdx].drive.src.value ?? 0) + 1.5;
    const br1 = d1.branch({ terminal: term, to, maxSteps: 40, stability: false });
    assert.match(br1.end, /^reached to/, `echem_${k} ${term}: ${br1.end}`);
    assert.ok(br1.points.length >= 5, `echem_${k} ${term} default took only ${br1.points.length} pts`);
    const br2 = d1.branch({ terminal: term, to, step: 0.1, maxSteps: 40, stability: false });
    assert.match(br2.end, /^reached to/, `echem_${k} ${term} step=0.1: ${br2.end}`);
    assert.ok(br2.points.length >= 16, `echem_${k} ${term} step=0.1 took only ${br2.points.length} pts`);
  }

  // 3. Depleted-pool & lost-mode combos (redoxgate_18, semi_37) + wetjunction_16 (no deep-LHP real-mode crawl).
  for (const [label, fam, k, term, to, maxPts] of [
    ['redoxgate_18', redoxgateFamily, 18, 'left', 1.5, 20],
    ['semi_37', semiFamily, 37, 'left', 1.5, 22],
    ['wetjunction_16', wetjunctionFamily, 16, 'right', -1.5, 25],
  ]) {
    const def = fam.make(caseRng(fam.name, k)).def;
    const d = new Device(structuredClone(def));
    assert.equal(d.solve({ budgetMs: 3000 }).converged, true);
    const br = d.branch({ terminal: term, to, maxSteps: 40, budgetMs: 20000 });
    assert.match(br.end, /^reached to/, `${label}: ${br.end}`);
    assert.ok(br.points.length <= maxPts, `${label}: ${br.points.length} pts > ${maxPts}`);
    assert.equal(br.events.filter((e) => e.kind === 'fold').length, 0, `${label}: false folds`);
  }

  // 4. Mode guard on landing step: narrow g.Na=835 window [0.445, 0.67] inside the landing step 0.35 -> 0.72.
  {
    const d = hhDevice(35, { g: { Na: 835 } });
    assert.equal(d.solve().converged, true);
    const br = d.branch({ terminal: 'right', to: 0.72, step: 0.4 });
    assert.match(br.end, /^reached to/, `g.Na=835 landing: ${br.end}`);
    const hopfs = br.events.filter((e) => e.kind === 'hopf');
    assert.equal(hopfs.length, 2, `g.Na=835 landing: expected 2 Hopfs, got ${JSON.stringify(br.events)}`);
  }
});

// Round-6 regression & oracle tests (every test fails on Round 5 84eaabd):

test('Si pn diode in reverse under current drive (0 -> -1e4): 0 folds, monotone V & I, reaches ±20 V bound', () => {
  const def = build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [
      ohmic(0),
      layer('Si', units.um(1), { donors: units.perCm3(1e16) }),
      layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }),
      ohmic({ I: 0 }),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    grid: { hmin: units.nm(5), hmax: units.nm(100) },
  });
  const d = new Device(def);
  assert.equal(d.solve().converged, true);
  for (const stability of [false, true]) {
    const br = d.branch({ terminal: 'right', to: -1e4, stability, maxSteps: 300 });
    assert.equal(br.events.length, 0, `Si pn reverse stab=${stability}: unexpected events ${JSON.stringify(br.events)}`);
    assert.match(br.end, /^left voltage bound /, `Si pn reverse stab=${stability}: ${br.end}`);
    assert.ok(br.points.length <= 200, `Si pn reverse stab=${stability}: crawled (${br.points.length} pts)`);
    for (let i = 1; i < br.points.length; i++) {
      assert.ok(br.points[i].drive <= br.points[i - 1].drive, `Si pn reverse I rose at pt ${i}`);
      assert.ok(br.points[i].terminals.right.V < br.points[i - 1].terminals.right.V, `Si pn reverse V rose at pt ${i}`);
      const kErr = Math.abs(br.points[i].terminals.left.current + br.points[i].terminals.right.current) / Math.abs(br.points[i].drive);
      assert.ok(kErr < 1e-6, `Si pn reverse stab=${stability} pt ${i}: Kirchhoff err ${kErr.toExponential(2)}`);
    }
    assert.ok(br.points.at(-1).terminals.right.V < -19, `Si pn reverse stopped short: V=${br.points.at(-1).terminals.right.V}`);
  }
});

test('Two-wave electrode walks across Fe3+ limiting-current plateau onto second wave (0 -> -100 A/m²)', () => {
  const water = aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 });
  water.species.push({ name: 'Ox2+', z: 2, cRef: 1000 }, { name: 'Red+', z: 1, cRef: 1000 });
  water.materials.water.species['Ox2+'] = { D: 0.7e-9, mu0: 0 };
  water.materials.water.species['Red+'] = { D: 0.7e-9, mu0: 2.0 * FARADAY };
  const def = build({
    library: [water, metal('Pt')],
    stack: [
      ohmic({ I: 0 }, ['e-']),
      layer('Pt', 1e-6),
      { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1 }, { equation: 'Ox2+ + e- = Red+', k0: 1e-12 }] },
      layer('water', 20e-6),
      bath({ 'K+': 500, 'Cl-': 635 + 1e-6, 'Fe3+': 5, 'Fe2+': 10, 'Ox2+': 50, 'Red+': 1e-6 }, AgAgCl),
    ],
    grid: { hmin: 10e-9, hmax: 0.5e-6 },
  });
  const d = new Device(def);
  assert.equal(d.solve().converged, true);
  for (const stability of [false, true]) {
    const br = d.branch({ terminal: 'left', to: -100, stability, maxSteps: 200 });
    assert.match(br.end, /^reached to/, `two-wave stab=${stability}: ${br.end}`);
    assert.equal(br.events.length, 0, `two-wave stab=${stability}: ${JSON.stringify(br.events)}`);
    const last = br.points.at(-1);
    assert.ok(Math.abs(last.drive - -100) < 1e-6, `two-wave drive=${last.drive}`);
    assert.ok(Math.abs(last.terminals.left.V - -3.47907) < 1e-3, `two-wave V=${last.terminals.left.V}`);
  }
});

test('Held-out membrane 5000/40/10 at 1.5 Rc refines undecided fold end and locates all 4 events in both directions', () => {
  const g = { Na: 5000, K: 40, leak: 10 };
  const vDev = hhDevice(0, { g, right: { V: -0.09 } });
  assert.equal(vDev.solve().converged, true);
  const IofV = (V) => { vDev.set({ contacts: { right: { V } } }); assert.equal(vDev.solve().converged, true); return vDev.solver._currentInto(1); };
  let minSlope = 0;
  for (let v = -0.095; v <= -0.005; v += 0.0005) {
    const sl = (IofV(v + 0.0005) - IofV(v)) / 0.0005;
    if (sl < minSlope) minSlope = sl;
  }
  const R = 1.5 * (-1 / minSlope);
  const EofV = (V) => V + R * IofV(V);
  const golden = (f, a, b, isMax) => {
    const phi = 0.5 * (Math.sqrt(5) - 1);
    let x1 = b - phi * (b - a), x2 = a + phi * (b - a), f1 = f(x1), f2 = f(x2);
    for (let it = 0; it < 38; it++) {
      if (isMax ? f1 < f2 : f1 > f2) { a = x1; x1 = x2; f1 = f2; x2 = a + phi * (b - a); f2 = f(x2); }
      else { b = x2; x2 = x1; f2 = f1; x1 = b - phi * (b - a); f1 = f(x1); }
    }
    const v = 0.5 * (a + b);
    return { v, val: f(v) };
  };
  const vs = [], es = [];
  for (let i = 0; i <= 120; i++) { const v = -0.07 + 0.05 * (i / 120); vs.push(v); es.push(EofV(v)); }
  const refFolds = [];
  for (let i = 1; i < 120; i++) {
    if ((es[i] > es[i - 1] && es[i] > es[i + 1]) || (es[i] < es[i - 1] && es[i] < es[i + 1])) {
      const ex = golden(EofV, vs[i - 1], vs[i + 1], es[i] > es[i - 1]);
      refFolds.push({ kind: 'fold', v: ex.v, drive: ex.val });
    }
  }
  const hopfRe = (v) => {
    vDev.set({ contacts: { right: { V: v, R: undefined } } });
    assert.equal(vDev.solve().converged, true);
    const I = vDev.solver._currentInto(1), E = v + R * I;
    const snap = vDev.solver._snapshot();
    vDev.set({ contacts: { right: { V: E, R } } });
    vDev.solver._restore(snap); vDev.solver.termV[1] = v; vDev.solver.atSteady = true;
    const st = vDev.stability();
    vDev.set({ contacts: { right: { V: v, R: undefined } } });
    vDev.solver._restore(snap); vDev.solver.atSteady = true;
    return { osc: st.modes.find((m) => m.im > 10 * m.error), drive: E };
  };
  const refHopfs = [];
  let prev = null;
  for (let i = 0; i <= 70; i++) {
    const v = -0.07 + 0.05 * (i / 70);
    const cur = { v, ...hopfRe(v) };
    if (prev && prev.osc && cur.osc && prev.osc.re * cur.osc.re < 0) {
      let a = prev.v, b = cur.v, ra = prev.osc.re, lastD = cur.drive;
      for (let it = 0; it < 20; it++) {
        const mid = 0.5 * (a + b), m = hopfRe(mid); lastD = m.drive;
        if (!m.osc) break;
        if (m.osc.re * ra > 0) { a = mid; ra = m.osc.re; } else b = mid;
      }
      refHopfs.push({ kind: 'hopf', v: 0.5 * (a + b), drive: lastD });
    }
    prev = cur;
  }
  const refAll = [...refFolds, ...refHopfs].sort((a, b) => a.v - b.v);
  assert.equal(refAll.length, 4, `5000/40/10 ref events: ${JSON.stringify(refAll)}`);
  const Emax = Math.max(...refAll.map((e) => e.drive)), Emin = Math.min(...refAll.map((e) => e.drive));
  for (const [label, E0, to, rev] of [
    ['5000/40/10 1.5Rc up', Emin - 0.375, Emax + 0.375, false],
    ['5000/40/10 1.5Rc down', Emax + 0.375, Emin - 0.375, true],
    ['5000/40/10 1.5Rc down (undecided fold)', 0.31949482, -0.438868587, true],
  ]) {
    const d = hhDevice(0, { g, right: { V: E0, R } });
    assert.equal(d.solve().converged, true);
    const br = d.branch({ terminal: 'right', to, maxSteps: 300 });
    assert.match(br.end, /^reached to/, `${label}: ${br.end}`);
    const exp = rev ? [...refAll].reverse() : refAll;
    assert.equal(br.events.length, 4, `${label}: got ${JSON.stringify(br.events.map((e) => [e.kind, e.drive]))}`);
    for (let i = 0; i < 4; i++) {
      assert.equal(br.events[i].kind, exp[i].kind, `${label} ev[${i}]`);
      let minGap = Infinity;
      for (let j = 0; j < 4; j++) if (j !== i) minGap = Math.min(minGap, Math.abs(exp[j].drive - exp[i].drive));
      const cap = Math.min(2e-4 * Math.abs(exp[i].drive), 0.25 * minGap);
      assert.ok(Math.abs(br.events[i].drive - exp[i].drive) <= cap,
        `${label} ev[${i}]: got ${br.events[i].drive} vs ref ${exp[i].drive} (tol ${cap.toExponential(2)})`);
    }
  }
});

test('GaAs pn forward (0 -> 10) and floating-base Si npn (check() clean, 0 -> 1e6): 0 folds and 0 points outside span', () => {
  const gaas = new Device(build({
    T: 300,
    library: [semiconductor('GaAs')],
    stack: [
      ohmic(0),
      layer('GaAs', units.um(1), { donors: units.perCm3(1e16) }),
      layer('GaAs', units.um(1), { acceptors: units.perCm3(1e16) }),
      ohmic({ I: 0 }),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { GaAs: 1e8 } }],
    grid: { hmin: units.nm(5), hmax: units.nm(100) },
  }));
  assert.equal(gaas.solve().converged, true);
  const brG = gaas.branch({ terminal: 'right', to: 10, stability: false });
  assert.match(brG.end, /^reached to/, `GaAs forward: ${brG.end}`);
  assert.equal(brG.events.length, 0, `GaAs forward events: ${JSON.stringify(brG.events)}`);
  assert.ok(Math.abs(brG.points.at(-1).terminals.right.V - 0.7597) < 1e-3, `GaAs V=${brG.points.at(-1).terminals.right.V}`);

  const npn = new Device(build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [
      ohmic(0),
      layer('Si', units.um(1), { donors: units.perCm3(1e18) }),
      layer('Si', units.um(1), { acceptors: units.perCm3(1e17) }),
      layer('Si', units.um(3), { donors: units.perCm3(1e15) }),
      ohmic({ I: 0 }),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e6 } }],
    grid: { hmin: units.nm(2), hmax: units.nm(50) },
  }));
  const solNpn = npn.solve();
  assert.equal(solNpn.converged, true);
  assert.equal(check(npn, solNpn).ok, true, 'Si npn should pass check()');
  const brN = npn.branch({ terminal: 'right', to: 1e6, stability: false, maxSteps: 300 });
  assert.match(brN.end, /^left voltage bound /, `Si npn: ${brN.end}`);
  assert.equal(brN.events.length, 0, `Si npn events: ${JSON.stringify(brN.events)}`);
  const outside = brN.points.filter((p) => p.drive < 0 || p.drive > 1e6);
  assert.equal(outside.length, 0, `Si npn points outside [0, 1e6]: ${outside.length}`);
});

// Round-7 regression, self-check, to-independence, and guard tests (each fails on Round 6 134ea98):

test('Semiconductor stacks under current drive at 3 targets (-1e4, -1e8, -1e12) from 2 starts: 0 folds, Kirchhoff & voltage-held match, and to-independence', () => {
  const mkStack = (mat, donors, acceptors, kf, grid) => (right) => new Device(build({
    T: 300,
    library: [semiconductor(mat)],
    stack: [
      ohmic(0),
      layer(mat, units.um(1), { donors: units.perCm3(donors) }),
      layer(mat, units.um(1), { acceptors: units.perCm3(acceptors) }),
      ohmic(right),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { [mat]: kf } }],
    grid: { hmin: units.nm(grid[0]), hmax: units.nm(grid[1]) },
  }));
  const stacks = [
    ['Si p+n+ fine', mkStack('Si', 1e19, 1e19, 1e6, [0.5, 25])],
    ['GaAs pn', mkStack('GaAs', 1e16, 1e16, 1e8, [5, 250])],
    ['suite Si pn', mkStack('Si', 1e16, 1e16, 1e8, [5, 100])],
  ];
  for (const [name, mk] of stacks) {
    const vRef = mk(0);
    assert.equal(vRef.solve().converged, true);
    vRef.set({ contacts: { right: { V: -0.5 } } });
    const sRev = vRef.solve();
    assert.equal(sRev.converged, true);
    const Irev = sRev.terminals.right.current;
    for (const [startLabel, startV] of [['open circuit', null], ['reverse -0.5 V', -0.5]]) {
      for (const to of [-1e4, -1e8, -1e12]) {
        let d;
        if (startV === null) {
          d = mk({ I: 0 });
          assert.equal(d.solve().converged, true);
        } else {
          d = mk(startV);
          assert.equal(d.solve().converged, true);
          const I0 = d.solver._currentInto(1);
          d.set({ contacts: { right: { I: I0 } } });
          assert.equal(d.solve().converged, true);
        }
        const br = d.branch({ terminal: 'right', to, stability: false, maxSteps: 300 });
        assert.equal(br.events.filter((e) => e.kind === 'fold').length, 0, `${name} (${startLabel}) to ${to}: invented fold ${JSON.stringify(br.events)}`);
        if (startV === null) {
          assert.match(br.end, /^left voltage bound /, `${name} (${startLabel}) to ${to}: ${br.end}`);
        } else {
          assert.match(br.end, /^(left voltage bound |maxSteps reached )/, `${name} (${startLabel}) to ${to}: ${br.end}`);
          assert.ok(br.points.at(-1).terminals.right.V < -3, `${name} (${startLabel}) to ${to}: stopped short at V=${br.points.at(-1).terminals.right.V}`);
        }
        for (let i = 1; i < br.points.length; i++) {
          const pt = br.points[i];
          const kRel = Math.abs(pt.terminals.left.current + pt.terminals.right.current) / Math.max(Math.abs(pt.drive), 1e-300);
          assert.ok(kRel < 1e-6, `${name} (${startLabel}) to ${to} pt ${i}: Kirchhoff rel ${kRel.toExponential(2)}`);
          // Spot-check voltage-held solve() across the whole curve (start, every 8th point, and last point).
          if (i <= 3 || i % 8 === 0 || i === br.points.length - 1) {
            vRef.set({ contacts: { right: { V: pt.terminals.right.V } } });
            const sol = vRef.solve();
            assert.equal(sol.converged, true);
            const vRel = Math.abs(sol.terminals.right.current - pt.drive) / Math.max(Math.abs(pt.drive), Math.abs(sol.terminals.right.current), 1e-300);
            assert.ok(vRel < 1e-6, `${name} (${startLabel}) to ${to} pt ${i} (V=${pt.terminals.right.V.toFixed(4)}): off voltage-held by ${vRel.toExponential(2)}`);
          }
        }
      }
    }
  }
  // Points must not depend on `to`: same device gives matching points for to: -1 and to: -1e4.
  const mkFine = stacks[0][1];
  const d1 = mkFine({ I: 0 }); d1.solve();
  const d2 = mkFine({ I: 0 }); d2.solve();
  const b1 = d1.branch({ terminal: 'right', to: -1, stability: false, maxSteps: 300 });
  const b2 = d2.branch({ terminal: 'right', to: -1e4, stability: false, maxSteps: 300 });
  assert.equal(b1.points.length, b2.points.length, `Si p+n+ fine point count depended on to: ${b1.points.length} vs ${b2.points.length}`);
  for (let i = 1; i < b1.points.length; i++) {
    assert.ok(Math.abs(b1.points[i].terminals.right.V - b2.points[i].terminals.right.V) < 1e-9, `to-independence V diff at pt ${i}`);
    const relI = Math.abs(b1.points[i].drive - b2.points[i].drive) / Math.abs(b1.points[i].drive);
    assert.ok(relI < 1e-8, `to-independence drive diff at pt ${i}: ${relI.toExponential(2)}`);
  }
});

test('150 K cold Si diode forward (0 -> 10 A/m²): reaches target, 0 folds, Kirchhoff & voltage-held match across 27 decades', () => {
  const T = 150;
  const mk = (right) => new Device(build({
    T,
    library: [semiconductor('Si', { T })],
    stack: [
      ohmic(0),
      layer('Si', units.um(1), { donors: units.perCm3(1e16) }),
      layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }),
      ohmic(right),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    grid: { hmin: units.nm(5), hmax: units.nm(250) },
  }));
  const vRef = mk(0);
  assert.equal(vRef.solve().converged, true);
  for (const stability of [false, true]) {
    const d = mk({ I: 0 });
    assert.equal(d.solve().converged, true);
    const br = d.branch({ terminal: 'right', to: 10, stability, maxSteps: 300 });
    assert.match(br.end, /^reached to/, `150K Si stab=${stability}: ${br.end}`);
    assert.equal(br.events.length, 0, `150K Si stab=${stability}: unexpected events ${JSON.stringify(br.events)}`);
    assert.ok(br.points.length >= 15, `150K Si stab=${stability}: too few points (${br.points.length})`);
    for (let i = 1; i < br.points.length; i++) {
      const pt = br.points[i];
      const kRel = Math.abs(pt.terminals.left.current + pt.terminals.right.current) / Math.abs(pt.drive);
      assert.ok(kRel < 1e-6, `150K Si stab=${stability} pt ${i}: Kirchhoff rel ${kRel.toExponential(2)}`);
      if (i <= 5 || i % 10 === 0 || i === br.points.length - 1) {
        vRef.set({ contacts: { right: { V: pt.terminals.right.V } } });
        const sol = vRef.solve();
        assert.equal(sol.converged, true);
        const vRel = Math.abs(sol.terminals.right.current - pt.drive) / Math.max(Math.abs(pt.drive), Math.abs(sol.terminals.right.current));
        assert.ok(vRel < 1e-6, `150K Si stab=${stability} pt ${i}: off voltage-held by ${vRel.toExponential(2)}`);
      }
    }
  }
});

test('Provisional fold confirmation drops endpoint false fold (Si npn to -1e-6) and tangent-turn bound catches 1800/20/2.5 at 1.5 Rc without marks', () => {
  const npn = new Device(build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [
      ohmic(0),
      layer('Si', units.um(1), { donors: units.perCm3(1e19) }),
      layer('Si', units.um(1), { acceptors: units.perCm3(1e17) }),
      layer('Si', units.um(3), { donors: units.perCm3(1e15) }),
      ohmic({ I: 0 }),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e6 } }],
    grid: { hmin: units.nm(5), hmax: units.nm(250) },
  }));
  assert.equal(npn.solve().converged, true);
  const brNpn = npn.branch({ terminal: 'right', to: -1e-6, stability: false, maxSteps: 300 });
  assert.equal(brNpn.events.filter((e) => e.kind === 'fold').length, 0, `Si npn -1e-6 invented fold: ${JSON.stringify(brNpn.events)}`);

  const mem = hhDevice(0, { g: { Na: 1800, K: 20, leak: 2.5 }, right: { V: 0.317, R: 1.5 * 0.0188418 } });
  assert.equal(mem.solve().converged, true);
  const brMem = mem.branch({ terminal: 'right', to: -0.438803, stability: false, maxSteps: 400 });
  assert.match(brMem.end, /^reached to/, `1800/20/2.5 1.5Rc: ${brMem.end}`);
  const folds = brMem.events.filter((e) => e.kind === 'fold');
  assert.equal(folds.length, 2, `1800/20/2.5 1.5Rc no marks: expected 2 folds, got ${JSON.stringify(brMem.events)}`);
  assert.ok(Math.abs(folds[0].drive - -0.0638027) < 1e-4, `fold 0 drive ${folds[0].drive}`);
  assert.ok(Math.abs(folds[1].drive - -0.0598846) < 1e-4, `fold 1 drive ${folds[1].drive}`);
});

test('Defense-in-depth guards: law-occupied circuit (R05), start flip (R08), span guard (R09), undecided-fold audit (R11), plateau floor (R15), and acceptance equation + provisional-fold guards', () => {
  // 1. R05: cutEligible refuses to substitute an interior current cut on a floating contact whose circuit row holds a conservation law (cs.circuit === 1, natural gross/diag > 2e3).
  const blocked = new Device(build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [
      ohmic(0, ['e-']),
      layer('Si', units.um(1), { donors: units.perCm3(1e19) }),
      layer('Si', units.um(5), { acceptors: units.perCm3(1e19) }),
      ohmic({ I: 0 }, ['h+']),
    ],
    grid: { hmin: units.nm(5), hmax: units.nm(100) },
  }));
  assert.equal(blocked.solve().converged, true);
  assert.ok(blocked.solver.constraints.some((cs) => cs.circuit === 1), 'expected conservation law on circuit 1');
  let cutCalls5 = 0;
  const origCuts5 = blocked.solver._currentCuts.bind(blocked.solver);
  blocked.solver._currentCuts = () => { cutCalls5++; return origCuts5(); };
  const brBlk = blocked.branch({ terminal: 'left', to: 0.5, stability: false });
  blocked.solver._currentCuts = origCuts5;
  assert.equal(cutCalls5, 0, `R05 cutEligible failed: _currentCuts called ${cutCalls5} times on law-occupied circuit 1`);
  assert.match(brBlk.end, /^reached to/, `R05 blocked cell: ${brBlk.end}`);
  const dVRight = brBlk.points.at(-1).terminals.right.V - brBlk.points[0].terminals.right.V;
  assert.ok(Math.abs(dVRight - 0.5) < 1e-5, `R05 floating open-circuit shift expected +0.5 V, got ${dVRight}`);

  // 2. R08: fault-injected inverted initial tangent (sNext < 0 on step 0) flips cleanly and completes all 10 steps.
  const d8 = new Device(build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [
      ohmic(0),
      layer('Si', units.um(1), { donors: units.perCm3(1e16) }),
      layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }),
      ohmic({ I: 0 }),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    grid: { hmin: units.nm(5), hmax: units.nm(100) },
  }));
  assert.equal(d8.solve().converged, true);
  const s8 = d8.solver;
  const origSolve8 = s8._solveLinear.bind(s8);
  const origCuts8 = s8._currentCuts.bind(s8);
  let zc8 = 0, cutsRef8 = null, flipped8 = false;
  s8._solveLinear = function (rhs, out) {
    let isZ = true;
    for (let j = 0; j < this.sys.size; j++) if (rhs[j] !== 0) { isZ = false; break; }
    if (isZ) {
      zc8++;
      if (zc8 === 3 && cutsRef8) {
        for (const c of cutsRef8) for (let j = 0; j < c.g.length; j++) c.g[j] *= -0.001;
        flipped8 = true;
      }
    }
    return origSolve8(rhs, out);
  };
  s8._currentCuts = function () {
    const cuts = origCuts8();
    if (zc8 === 2) {
      cutsRef8 = cuts;
      for (const c of cuts) for (let j = 0; j < c.g.length; j++) c.g[j] *= 10;
    }
    return cuts;
  };
  const br8 = d8.branch({ terminal: 'right', to: -1e4, stability: false, maxSteps: 10 });
  s8._solveLinear = origSolve8;
  s8._currentCuts = origCuts8;
  assert.equal(flipped8, true, 'R08 fault injection did not trigger');
  assert.equal(br8.points.length, 11, `R08 start flip failed to continue: pts=${br8.points.length}, end=${br8.end}`);
  assert.ok(br8.points.at(-1).terminals.right.V < -1, `R08 final V=${br8.points.at(-1).terminals.right.V}`);

  // 3. R09: starting at the exact Fe fold and asking to trace past the fold into the empty region stops at point 0 without stepping out of span.
  const dFeRef = new Device(feSpec({ I: -62.90054296647188 })); dFeRef.solve();
  const foldRef = dFeRef.branch({ terminal: 'right', to: -120, stability: false }).events.find((e) => e.kind === 'fold');
  const dFeV = new Device(feSpec({ V: foldRef.terminals.right.V })); dFeV.solve();
  const IofV = (V) => { dFeV.set({ contacts: { right: { V } } }); dFeV.solve(); return dFeV.solver._currentInto(1); };
  let va = foldRef.terminals.right.V - 1e-4, vb = foldRef.terminals.right.V + 1e-4;
  const phi = 0.5 * (Math.sqrt(5) - 1);
  let vx1 = vb - phi * (vb - va), vx2 = va + phi * (vb - va), f1 = IofV(vx1), f2 = IofV(vx2);
  for (let it = 0; it < 28; it++) {
    if (f1 > f2) { va = vx1; vx1 = vx2; f1 = f2; vx2 = va + phi * (vb - va); f2 = IofV(vx2); }
    else { vb = vx2; vx2 = vx1; f2 = f1; vx1 = vb - phi * (vb - va); f1 = IofV(vx1); }
  }
  const Vexact = 0.5 * (va + vb);
  const Ifold = IofV(Vexact);
  const d9 = new Device(feSpec({ I: Ifold }));
  d9.solver._restore(dFeV.solver._snapshot());
  d9.solver.termV[1] = Vexact;
  d9.solver.computeConcentrations();
  d9.solver._lawsNow();
  d9.solver._assembleBookkeeping(Infinity);
  d9.solver.atSteady = true;
  d9.solver.regriddedSteady = false;
  const s9 = d9.solver;
  const origSolve9 = s9._solveLinear.bind(s9);
  let zc9 = 0;
  s9._solveLinear = function (rhs, out) {
    let isZ = true;
    for (let j = 0; j < this.sys.size; j++) if (rhs[j] !== 0) { isZ = false; break; }
    if (isZ && ++zc9 === 1) {
      for (let j = 0; j < this.sys.size; j++) this.termC[1][j] *= 1e-100;
      this.termDI[1] *= 1e-100;
    }
    return origSolve9(rhs, out);
  };
  const br9 = d9.branch({ terminal: 'right', to: -120, stability: false });
  s9._solveLinear = origSolve9;
  assert.equal(br9.points.length, 1, `R09 span guard allowed out-of-span step: pts=${br9.points.length}, end=${br9.end}`);
  assert.match(br9.end, /could not advance into span/, `R09 unexpected end: ${br9.end}`);

  // 4. R11: when a fold's endpoints remain 'undecided' even after MAX_REFINE halvings, runAudit refuses to credit the fold and emits 'stability mismatch'.
  const d11 = new Device(feSpec({ I: -62.90054296647188 })); d11.solve();
  const origStab11 = d11.solver.stability.bind(d11.solver);
  d11.solver.stability = function (opts) {
    const V = this.termV[1];
    if (V < -1.08 && V > -1.30) return { verdict: 'undecided', modes: [], reason: 'injected undecided around fold' };
    return origStab11(opts);
  };
  const br11 = d11.branch({ terminal: 'right', to: -120 });
  d11.solver.stability = origStab11;
  assert.ok(br11.events.some((e) => e.kind === 'fold'), `R11 expected fold event: ${JSON.stringify(br11.events)}`);
  assert.ok(br11.events.some((e) => e.kind === 'stability mismatch'), `R11 audit failed to emit stability mismatch for undecided-sided fold: ${JSON.stringify(br11.events)}`);

  // 5. R15: drive round-off floor suppresses 16-digit plateau jitter on echem_32 ref (0 folds, reaches voltage bound).
  const e32Built = echemFamily.make(caseRng('echem', 32));
  const d15 = new Device(e32Built.def ?? e32Built);
  assert.equal(d15.solve().converged, true);
  const br15 = d15.branch({ terminal: 'ref', to: 1, stability: false, maxSteps: 200 });
  assert.equal(br15.events.length, 0, `R15 echem_32 ref plateau invented fold: ${JSON.stringify(br15.events)}`);
  assert.match(br15.end, /^left voltage bound /, `R15 echem_32 ref: ${br15.end}`);

  // 6. Acceptance reported-equation / 2-terminal Kirchhoff guard rejects a corrupted terminal flux instead of accepting it silently.
  const dEq = hhDevice(0);
  assert.equal(dEq.solve().converged, true);
  const origBkEq = dEq.solver._assembleBookkeeping.bind(dEq.solver);
  dEq.solver._assembleBookkeeping = function (dt) {
    origBkEq(dt);
    if (this.sourceOverride.has(1)) this.contactFlux.left[0] += 0.01 / FARADAY;
  };
  const brEq = dEq.branch({ terminal: 'right', to: 0.05, stability: false });
  dEq.solver._assembleBookkeeping = origBkEq;
  assert.equal(brEq.points.length, 1, `acceptance equation guard failed to reject corrupted Kirchhoff: pts=${brEq.points.length}, end=${brEq.end}`);
  assert.match(brEq.end, /reported terminal equation not satisfied/, `acceptance equation guard end: ${brEq.end}`);

  // 7. Provisional fold confirmation drops a single-step tangent-noise fold whose next accepted point keeps advancing forward.
  const dPf = new Device(build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [
      ohmic(0),
      layer('Si', units.um(1), { donors: units.perCm3(1e16) }),
      layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }),
      ohmic({ I: 0 }),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    grid: { hmin: units.nm(5), hmax: units.nm(100) },
  }));
  assert.equal(dPf.solve().converged, true);
  const sPf = dPf.solver;
  const origSolvePf = sPf._solveLinear.bind(sPf);
  const origCutsPf = sPf._currentCuts.bind(sPf);
  let zcPf = 0, flipNextCut = false, firedPf = false;
  sPf._solveLinear = function (rhs, out) {
    let isZ = true;
    for (let j = 0; j < this.sys.size; j++) if (rhs[j] !== 0) { isZ = false; break; }
    if (isZ && ++zcPf === 8) flipNextCut = true;
    return origSolvePf(rhs, out);
  };
  sPf._currentCuts = function () {
    const cuts = origCutsPf();
    if (flipNextCut) {
      flipNextCut = false;
      firedPf = true;
      for (const c of cuts) for (let j = 0; j < c.g.length; j++) c.g[j] = -c.g[j];
    }
    return cuts;
  };
  const brPf = dPf.branch({ terminal: 'right', to: -1e4, stability: false, maxSteps: 20 });
  sPf._solveLinear = origSolvePf;
  sPf._currentCuts = origCutsPf;
  assert.equal(firedPf, true, 'provisional-fold fault injection did not trigger');
  assert.equal(brPf.events.length, 0, `provisional-fold confirmation failed to drop endpoint noise fold: ${JSON.stringify(brPf.events)}`);
});

test('A fold on the trace\'s last step is kept when the accepted points bracket the turn, with a note', () => {
  // Fable r7 finding 1: a fold step that advances in drive leaves the fold provisional; cut the
  // trace at that step (maxSteps = b) and round 7 dropped a real fold. The full trace is the truth.
  const mk = () => { const d = hhDevice(0, { g: { Na: 3000, K: 60, leak: 4 }, right: { I: -2.2 } }); assert.equal(d.solve().converged, true); return d; };
  const full = mk().branch({ to: 0.94, stability: false });
  const folds = full.events.filter((e) => e.kind === 'fold');
  assert.equal(folds.length, 2, JSON.stringify(full.events));
  let advancing = 0;
  for (const f of folds) {
    const [a, b] = f.between;
    const cut = mk().branch({ to: 0.94, stability: false, maxSteps: b });
    assert.match(cut.end, /maxSteps/);
    const got = cut.events.filter((e) => e.kind === 'fold' && e.between[0] === a && e.between[1] === b);
    assert.equal(got.length, 1, `fold between [${a},${b}] dropped when the trace ends on its step: ${JSON.stringify(cut.events)}`);
    assert.equal(got[0].drive, f.drive);
    const adv = Math.sign(full.points[b].drive - full.points[a].drive) === Math.sign(full.points[a].drive - (a > 0 ? full.points[a - 1].drive : -2.2));
    if (adv) { advancing++; assert.match(got[0].note ?? '', /ended before a later point confirmed/); }
    // events stay in trace order
    for (let i = 1; i < cut.events.length; i++) assert.ok(cut.events[i].between[0] >= cut.events[i - 1].between[0]);
  }
  assert.ok(advancing >= 1, 'no fold step advanced in drive, so the provisional path was not exercised');

  // Fe cell with marks on, cut on its fold step: the fold is reported and no mismatch claims none was located.
  const I0 = -62.90054296647188;
  const fe = () => { const d = new Device(feSpec({ V: -0.95 })); d.solve(); d.set({ contacts: { right: { I: I0 } } }); assert.equal(d.solve().converged, true); return d; };
  const feFull = fe().branch({ terminal: 'right', to: -120 });
  const ff = feFull.events.find((e) => e.kind === 'fold');
  assert.ok(ff, JSON.stringify(feFull.events));
  const feCut = fe().branch({ terminal: 'right', to: -120, maxSteps: ff.between[1] });
  assert.ok(feCut.events.some((e) => e.kind === 'fold' && e.drive === ff.drive), `Fe fold dropped: ${JSON.stringify(feCut.events)}`);
  assert.equal(feCut.events.filter((e) => e.kind === 'stability mismatch').length, 0, JSON.stringify(feCut.events));
});

test('A fold on an already-unstable branch: the fold\'s own critical mode is not counted (Frumkin electrode behind R)', () => {
  // S2O8²⁻ at Hg in dilute NaF (demos/bistable.html's device) behind R = 0.05 Ω·m²: the busy branch
  // goes unstable at a Hopf point, and its two real growing modes meet the upper fold. At the
  // located fold one of them is ~0 of either sign (+8e-5 /s here, with an error bar of 4e-11), and
  // counting it made every fold step "contain an additional crossing": round 7 halved to a stop
  // at the fold. The folds must be where the load line touches the potentiostatic curve.
  const F = FARADAY, T = 298.15, RT = GAS_CONSTANT * T, NA = 6.02214076e23;
  const ions = { 'Na+': [1, 1.334e-9, -261905], 'F-': [-1, 1.475e-9, -278790], 'S2O82-': [-2, 1.145e-9, -1114900], 'SO42-': [-2, 1.065e-9, -744530] };
  const names = Object.keys(ions), cMax = 1 / (NA * 1e-30), L = 20e-6;
  const bulk = { 'Na+': 3 + 4e-5, 'F-': 1, 'S2O82-': 1, 'SO42-': 2e-5 };
  const she = -(ions['F-'][2] + RT * Math.log(bulk['F-'] / 1000 / (1 - Object.values(bulk).reduce((s, c) => s + c, 0) / cMax))) / F;
  const electrode = (U, R) => build({
    T,
    library: [{
      species: [{ name: 'e-', z: -1 }, ...names.map((name) => ({ name, z: ions[name][0], cRef: 1000 }))],
      materials: {
        Hg: { conductor: { species: 'e-', conductivity: 1.04e6 } },
        water: { epsr: 0, species: Object.fromEntries(names.map((n) => [n, { D: ions[n][1], mu0: ions[n][2] }])), statistics: [{ type: 'lattice', species: names, cMax }] },
      },
    }],
    stack: [
      { ...ohmic(U, ['e-']), ...(R > 0 ? { R } : {}) },
      layer('Hg', 1e-6),
      { phi: { type: 'capacitive', C: 0.17 }, zeroCharge: -0.19, diffuse: { epsr: 78.3 }, reactions: [{ ...half('S2O82- + 2 e- = 2 SO42-'), k0: 2.5e-9, alpha: 0.125 }] },
      layer('water', L, { c0: bulk }),
      bath(bulk, 'F-', 0, { offset: she }),
    ],
    grid: { hmin: 100e-9, hmax: L / 30, ratio: 1.2 },
  });
  const R = 0.05;
  const d = new Device(electrode(0.3, R));
  assert.equal(d.solve().converged, true);
  const br = d.branch({ terminal: 'left', to: -1.6, maxSteps: 800 });
  assert.match(br.end, /^reached to/);
  const folds = br.events.filter((e) => e.kind === 'fold');
  assert.equal(folds.length, 2, JSON.stringify(br.events.map((e) => [e.kind, e.drive])));
  assert.equal(br.events.filter((e) => e.kind === 'stability mismatch').length, 0);
  // Independent: the potentiostatic curve, and U = E + R·I along it (I as reported), turns where
  // the folds are.
  const h = new Device(electrode(0.3, 0)); h.solve();
  const u = [];
  for (let E = -0.15; E >= -0.45; E -= 0.0025) { h.set({ contacts: { left: { V: E } } }); const s = h.solve(); assert.ok(s.converged); u.push(E + R * s.terminals.left.current); }
  const turns = [];
  for (let k = 1; k + 1 < u.length; k++) if ((u[k] - u[k - 1]) * (u[k + 1] - u[k]) < 0) turns.push(u[k]);
  assert.equal(turns.length, 2, `load-line turns ${turns}`);
  for (const [k, f] of folds.entries()) assert.ok(Math.abs(f.drive - turns[k]) < 2e-3, `fold ${f.drive} against load-line turn ${turns[k]}`);
});
