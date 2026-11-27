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
import { aqueous, bath, build, hodgkinHuxley, layer, metal, ohmic, pulse, semiconductor } from '../src/kit.js';
import { caseRng } from '../bench/combos/rng.js';
import basicFamily from '../bench/combos/families/basic.js';
import echemFamily from '../bench/combos/families/echem.js';
import semiFamily from '../bench/combos/families/semi.js';
import neutralFamily from '../bench/combos/families/neutral.js';
import membraneFamily from '../bench/combos/families/membrane.js';

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

function hhDevice(I_uAcm2, { L = 1e-6, hmin = 20e-9, hmax = 0.2e-6 } = {}) {
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
  const hh = hodgkinHuxley({ T, linear: true });
  return new Device(
    build({
      T,
      library: [water],
      stack: [
        bath(out, 0),
        layer('water', L, { c0: out }),
        { phi: { type: 'capacitive', C: 0.01 }, gates: hh.gates, species: { ...hh.species, 'A-': 'blocked' } },
        layer('water', L, { c0: inside }),
        bath(inside, { I: I_uAcm2 * 1e-2 }),
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

  // Fixed-step kicked BDF2 runs resolving 1/ω (dt = 5e-5 s, ~34 steps per radian) on both sides
  // of the Hopf confirm damping below I* and growth above I*.
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
