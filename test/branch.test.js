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
import mixedfilmFamily from '../bench/combos/families/mixedfilm.js';
import redoxgateFamily from '../bench/combos/families/redoxgate.js';

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

function hhDevice(I_uAcm2, { L = 1e-6, hmin = 20e-9, hmax = 0.2e-6, g } = {}) {
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
  for (const [name, fam, k, side, sgn] of [
    ['echem_4', echemFamily, 4, 'left', 1],
    ['semi_1', semiFamily, 1, 'left', 1],
    ['membrane_2', membraneFamily, 2, 'left', 1],
  ]) {
    const def = fam.make(caseRng(name.split('_')[0], k)).def;
    const dev = new Device(structuredClone(def));
    assert.equal(dev.solve({ budgetMs: 3000 }).converged, true);
    const tIdx = dev.solver.terms.findIndex((tm) => tm.name === side);
    const p0 = dev.solver.terms[tIdx].drive.src.value ?? 0;
    const to = p0 + sgn * 0.3;
    const br = dev.branch({ terminal: side, to, maxSteps: 20, stability: false });
    assert.match(br.end, /^reached to/, `${name}: ${br.end}`);
    assert.ok(br.points.length >= 2);
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
  assert.match(sbr.end, /solve\(\) may continue from here/, `semi_5: ${sbr.end}`);
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
