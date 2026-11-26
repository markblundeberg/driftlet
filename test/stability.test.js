// Linear stability analysis (src/stability.js) and transposed block-tridiagonal solves
// (src/blockTridiagonal.js): four physics acceptance gates (Hodgkin–Huxley Hopf vs advance(),
// RC + Sol's C_face=1e16 cancellation cut, blocked diffusion discrete cosine spectrum + immobile
// local constraints, Fe³⁺/Fe²⁺ d26594 falling vs rising branch vs kicked advance()) plus contract
// checks (non-steady refusal, bit-for-bit state preservation, profiles, debug hook, verdicts).

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, SolverError } from '../src/index.js';
import { BlockTridiagonal, ComplexBlockTridiagonal } from '../src/blockTridiagonal.js';
import { build, layer, ohmic, metal, aqueous, bath, hodgkinHuxley, pulse } from '../src/kit.js';

// --- Helper: deterministic pseudo-random generator for linear algebra unit tests ---
function lcg(seed) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 2147483648 - 1;
  };
}

// --- Helper: Hodgkin–Huxley classic space-clamped fixture (E_Na=+50, E_K=-77, E_Cl=-54.387 mV) ---
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

// --- Helper: Fe3+/Fe2+ d26594 N-shaped I-V fixture (from test/homotopy.test.js) ---
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

// ============================================================================
// 0. BlockTridiagonal.solveT and ComplexBlockTridiagonal.solveT / solveH
// ============================================================================

test('BlockTridiagonal.solveT matches Mᵀ x = b for non-uniform blocks (both left and right sweeps)', () => {
  const sizes = [2, 3, 1, 2];
  const rand = lcg(0x123456);
  for (const forceRight of [false, true]) {
    const bt = new BlockTridiagonal(sizes.length, sizes);
    bt.bothWays = forceRight;
    const N = bt.size;
    const M = Array.from({ length: N }, () => new Float64Array(N));
    for (let i = 0; i < sizes.length; i++) {
      const m = sizes[i];
      if (i > 0) {
        const mp = sizes[i - 1];
        for (let r = 0; r < m; r++) for (let c = 0; c < mp; c++) {
          const v = rand() * 0.4;
          bt.A[bt.offA[i] + r * mp + c] = v;
          M[bt.offX[i] + r][bt.offX[i - 1] + c] = v;
        }
      }
      for (let r = 0; r < m; r++) for (let c = 0; c < m; c++) {
        // Make block 0 ill-conditioned when forceRight is true so bothWays picks fromRight.
        const diag = forceRight && i === 0 ? 1e-6 : 3;
        const v = (r === c ? diag : 0) + rand() * 0.4;
        bt.B[bt.offB[i] + r * m + c] = v;
        M[bt.offX[i] + r][bt.offX[i] + c] = v;
      }
      if (i < sizes.length - 1) {
        const mn = sizes[i + 1];
        for (let r = 0; r < m; r++) for (let c = 0; c < mn; c++) {
          const v = rand() * 0.4;
          bt.C[bt.offC[i] + r * mn + c] = v;
          M[bt.offX[i] + r][bt.offX[i + 1] + c] = v;
        }
      }
    }
    bt.factor();
    const b = Float64Array.from({ length: N }, () => rand());
    const x = bt.solveT(b);
    let maxErr = 0;
    for (let c = 0; c < N; c++) {
      let acc = 0;
      for (let r = 0; r < N; r++) acc += M[r][c] * x[r];
      maxErr = Math.max(maxErr, Math.abs(acc - b[c]));
    }
    assert.ok(maxErr < 1e-11, `solveT residual ${maxErr} (fromRight=${bt.fromRight})`);
  }
});

test('ComplexBlockTridiagonal.solveT and solveH match Mᵀ x = b and Mᴴ x = b', () => {
  const sizes = [2, 3, 1, 2];
  const rand = lcg(0x654321);
  const cbt = new ComplexBlockTridiagonal(sizes.length, sizes);
  const N = cbt.size;
  const Mr = Array.from({ length: N }, () => new Float64Array(N));
  const Mi = Array.from({ length: N }, () => new Float64Array(N));
  for (let i = 0; i < sizes.length; i++) {
    const m = sizes[i];
    if (i > 0) {
      const mp = sizes[i - 1];
      for (let r = 0; r < m; r++) for (let c = 0; c < mp; c++) {
        const vr = rand() * 0.3, vi = rand() * 0.3;
        cbt.Ar[cbt.offA[i] + r * mp + c] = vr;
        cbt.Ai[cbt.offA[i] + r * mp + c] = vi;
        Mr[cbt.offX[i] + r][cbt.offX[i - 1] + c] = vr;
        Mi[cbt.offX[i] + r][cbt.offX[i - 1] + c] = vi;
      }
    }
    for (let r = 0; r < m; r++) for (let c = 0; c < m; c++) {
      const vr = (r === c ? 3 : 0) + rand() * 0.3, vi = (r === c ? 1.5 : 0) + rand() * 0.3;
      cbt.Br[cbt.offB[i] + r * m + c] = vr;
      cbt.Bi[cbt.offB[i] + r * m + c] = vi;
      Mr[cbt.offX[i] + r][cbt.offX[i] + c] = vr;
      Mi[cbt.offX[i] + r][cbt.offX[i] + c] = vi;
    }
    if (i < sizes.length - 1) {
      const mn = sizes[i + 1];
      for (let r = 0; r < m; r++) for (let c = 0; c < mn; c++) {
        const vr = rand() * 0.3, vi = rand() * 0.3;
        cbt.Cr[cbt.offC[i] + r * mn + c] = vr;
        cbt.Ci[cbt.offC[i] + r * mn + c] = vi;
        Mr[cbt.offX[i] + r][cbt.offX[i + 1] + c] = vr;
        Mi[cbt.offX[i] + r][cbt.offX[i + 1] + c] = vi;
      }
    }
  }
  cbt.factor();
  const br = Float64Array.from({ length: N }, () => rand());
  const bi = Float64Array.from({ length: N }, () => rand());

  const [xtr, xti] = cbt.solveT(br, bi);
  let errT = 0;
  for (let c = 0; c < N; c++) {
    let ar = 0, ai = 0;
    for (let r = 0; r < N; r++) {
      ar += Mr[r][c] * xtr[r] - Mi[r][c] * xti[r];
      ai += Mr[r][c] * xti[r] + Mi[r][c] * xtr[r];
    }
    errT = Math.max(errT, Math.hypot(ar - br[c], ai - bi[c]));
  }
  assert.ok(errT < 1e-13, `solveT residual ${errT}`);

  const [xhr, xhi] = cbt.solveH(br, bi);
  let errH = 0;
  for (let c = 0; c < N; c++) {
    let ar = 0, ai = 0;
    for (let r = 0; r < N; r++) {
      // Mᴴ_{c,r} = conj(M_{r,c}) = Mr[r][c] - i Mi[r][c]
      ar += Mr[r][c] * xhr[r] + Mi[r][c] * xhi[r];
      ai += Mr[r][c] * xhi[r] - Mi[r][c] * xhr[r];
    }
    errH = Math.max(errH, Math.hypot(ar - br[c], ai - bi[c]));
  }
  assert.ok(errH < 1e-13, `solveH residual ${errH}`);
});

// ============================================================================
// 1. Acceptance Gate 1: Hodgkin–Huxley Hopf bifurcation vs advance()
// ============================================================================

test('Hodgkin–Huxley space-clamped: Hopf pair crosses Re λ = 0 between 9.7 and 9.8 µA/cm² and matches advance()', () => {
  const expected = [
    { I: 9.0, verdict: 'no growing mode found', reApprox: -13.913, imApprox: 578.14 },
    { I: 9.7, verdict: 'no growing mode found', reApprox: -0.5385, imApprox: 585.26 },
    { I: 9.8, verdict: 'unstable', reApprox: 1.3435, imApprox: 586.22 },
  ];
  for (const { I, verdict, reApprox, imApprox } of expected) {
    const dev = hhDevice(I);
    const sol = dev.solve();
    assert.equal(sol.converged, true);
    const st = dev.stability();
    assert.equal(st.verdict, verdict, `I=${I}: ${st.reason}`);
    assert.notEqual(st.verdict, 'stable');
    const [m0, m1] = st.modes;
    assert.ok(m0 && m1, 'expected complex conjugate pair');
    assert.ok(Math.abs(m0.re - reApprox) < 0.05, `I=${I}: Re λ=${m0.re} vs ${reApprox}`);
    assert.ok(Math.abs(Math.abs(m0.im) - imApprox) < 0.5, `I=${I}: Im λ=${m0.im} vs ${imApprox}`);
    assert.equal(m0.re, m1.re);
    assert.ok(Math.abs(m0.im + m1.im) < 1e-9);
    assert.ok(m0.residual < 1e-8, `I=${I}: residual=${m0.residual}`);
    assert.ok(m0.error < 1e-3, `I=${I}: error=${m0.error}`);
    assert.equal(st.warnings.length, 0);

    // Independent time-domain check via small-kick advance() (BDF2)
    const V0 = sol.terminals.right.V;
    dev.set({ contacts: { right: { I: I * 1e-2 + 1e-5 } } });
    dev.advance(1e-4, { method: 'bdf2', dtMax: 2e-5, tol: 1e-9 });
    dev.set({ contacts: { right: { I: I * 1e-2 } } });
    const dtSample = 2e-4;
    const ts = [], vs = [];
    for (let k = 1; k <= 250; k++) {
      const t = 1e-4 + k * dtSample;
      dev.advance(t, { method: 'bdf2', dtMax: 2.5e-5, tol: 1e-9 });
      ts.push(t);
      vs.push(dev.solution().terminals.right.V - V0);
    }
    const peaks = [];
    for (let k = 1; k < vs.length - 1; k++) {
      if (ts[k] > 0.008 && vs[k] > vs[k - 1] && vs[k] > vs[k + 1] && vs[k] > 0) {
        const a = vs[k - 1], b = vs[k], c = vs[k + 1];
        const delta = (0.5 * (a - c)) / (a - 2 * b + c);
        peaks.push({ t: ts[k] + delta * dtSample, v: b - 0.25 * (a - c) * delta });
      }
    }
    assert.ok(peaks.length >= 3, `I=${I}: expected >=3 oscillation peaks, got ${peaks.length}`);
    const p0 = peaks[0], p1 = peaks[peaks.length - 1];
    const rate = Math.log(p1.v / p0.v) / (p1.t - p0.t);
    const omega = (2 * Math.PI * (peaks.length - 1)) / (p1.t - p0.t);
    assert.ok(Math.abs(rate - m0.re) < 0.02, `I=${I}: advance rate ${rate} vs pencil ${m0.re}`);
    assert.ok(Math.abs(omega - Math.abs(m0.im)) < 0.2, `I=${I}: advance ω ${omega} vs pencil ${Math.abs(m0.im)}`);
  }
});

test('Hodgkin–Huxley with complex frequency shift recovers Hopf pair directly (no real shifts)', () => {
  const dev = hhDevice(9.8, { L: 50e-6 });
  dev.solve();
  const freqOnly = dev.stability({ shifts: [], frequencies: [93] });
  assert.equal(freqOnly.verdict, 'unstable');
  assert.equal(freqOnly.window.rHi, 0);
  assert.ok(freqOnly.window.omegaHi > 500);
  assert.ok(freqOnly.modes[0].re > 2.5 && freqOnly.modes[0].re < 4.0);
  assert.ok(Math.abs(Math.abs(freqOnly.modes[0].im) - 585.3) < 2);
  assert.ok(freqOnly.modes[0].residual < 1e-11, `residual=${freqOnly.modes[0].residual}`);
});

// ============================================================================
// 2. Acceptance Gate 2: RC & Sol's C_face=1e16 cancellation cut
// ============================================================================

test('RC: resistor-driven terminal + capacitive face + ohmic electrolyte matches −1/(RC)', () => {
  const aq = aqueous(['K+', 'Cl-'], { epsr: 0 });
  const Rext = 10, Cdl = 0.05;
  const dev = new Device(
    build({
      library: [metal('Au'), aq],
      stack: [
        ohmic({ V: 0, R: Rext }, ['e-']),
        layer('Au', 1e-6),
        { phi: { type: 'capacitive', C: Cdl }, zeroCharge: 0, species: { 'K+': 'blocked', 'Cl-': 'blocked' } },
        layer('water', 10e-6, { c0: { 'K+': 1000, 'Cl-': 1000 } }),
        bath({ 'K+': 1000, 'Cl-': 1000 }, 0),
      ],
      grid: { hmin: 20e-9, hmax: 1e-6 },
    }),
  );
  dev.solve();
  const st = dev.stability();
  assert.equal(st.verdict, 'no growing mode found');
  const T = 300, VT = (GAS_CONSTANT * T) / FARADAY;
  const sigma = (FARADAY / VT) * (aq.materials.water.species['K+'].D + aq.materials.water.species['Cl-'].D) * 1000;
  const Relec = 10e-6 / sigma;
  const exact = -1 / ((Rext + Relec) * Cdl);
  assert.ok(Math.abs(st.modes[0].re - exact) / Math.abs(exact) < 2e-3, `got ${st.modes[0].re}, expected ${exact}`);
  assert.equal(st.modes[0].im, 0);
  assert.ok(st.modes[0].residual < 1e-8);
});

test('RC dielectric + Sol C_face=1e16 cancellation cut + open/held drives', () => {
  const L = 1e-6, cgeo = (2 * EPS0) / L, R = 5;
  const cap = (C, drive) =>
    new Device({
      species: [{ name: 'X', z: 0, cRef: 1 }],
      materials: { d: { epsr: 2, species: { X: { D: 0, mu0: 0 } } } },
      regions: [{ material: 'd', length: L, c0: { X: 1 } }],
      grid: { minCells: 4, hmin: L / 4, hmax: L / 4 },
      contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0 }, right: { ...drive, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
    });

  // Series capacitance C_geom/2
  const d1 = cap(cgeo, { V: 0, R });
  d1.solve();
  const r1 = d1.stability();
  const exact1 = -1 / (R * (cgeo / 2));
  assert.equal(r1.verdict, 'no growing mode found');
  assert.equal(r1.modes.length, 1);
  assert.ok(Math.abs(r1.modes[0].re - exact1) / Math.abs(exact1) < 1e-12);

  // Sol's cancellation fixture: C_face = 1e16 cancels without cut-current reading
  const d2 = cap(1e16, { V: 0, R });
  d2.solve();
  const r2 = d2.stability({ debug: true });
  const exact2 = -1 / (R * cgeo);
  assert.equal(r2.verdict, 'no growing mode found');
  assert.equal(r2.modes.length, 1);
  assert.ok(r2.debug.cutUsed[0] >= 0, 'expected cut-current reading to engage');
  assert.ok(Math.abs(r2.modes[0].re - exact2) / Math.abs(exact2) < 1e-12, `got ${r2.modes[0].re}, exact ${exact2}`);

  // Open-circuit (I = 0): chargeOnly terminal deflated, no finite modes left
  const dOpen = cap(cgeo, { I: 0 });
  dOpen.solve();
  const rOpen = dOpen.stability();
  assert.equal(rOpen.verdict, 'no growing mode found');
  assert.equal(rOpen.modes.length, 0);
  assert.ok(rOpen.deflated.some((d) => d.kind === 'charge'), 'expected chargeOnly terminal deflated');

  // Held voltage (V = 0): no finite modes left
  const dHeld = cap(cgeo, { V: 0 });
  dHeld.solve();
  const rHeld = dHeld.stability();
  assert.equal(rHeld.verdict, 'no growing mode found');
  assert.equal(rHeld.modes.length, 0);
});

// ============================================================================
// 3. Acceptance Gate 3: Blocked neutral diffusion + immobile local invariants
// ============================================================================

test('Blocked neutral diffusion: exact discrete cosine roots −4D sin²(nπ/2N)/h² with n=0 deflated', () => {
  const L = 1e-6, D = 1e-9;
  for (const N of [4, 8, 16]) {
    const dev = new Device({
      species: [{ name: 'X', z: 0, cRef: 1 }],
      materials: { d: { epsr: 0, species: { X: { D, mu0: 0 } } } },
      regions: [{ material: 'd', length: L, c0: { X: 1 } }],
      grid: { minCells: N, hmin: L / N, hmax: L / N },
    });
    dev.solve();
    const st = dev.stability();
    assert.equal(st.verdict, 'no growing mode found');
    assert.equal(st.deflated.length, 1);
    assert.equal(st.deflated[0].kind, 'constraint');
    const h = L / N;
    for (let k = 1; k <= 3; k++) {
      const s = Math.sin((k * Math.PI) / (2 * N));
      const exact = (-4 * D * s * s) / (h * h);
      const got = st.modes[k - 1];
      assert.ok(Math.abs(got.re - exact) / Math.abs(exact) < 1e-12, `N=${N} n=${k}: got ${got.re} vs ${exact}`);
      assert.equal(got.im, 0);
    }
  }
});

test('Immobile species (D=0) deflates all node-local invariants; neutral equal-D salt matches diffusion spectrum', () => {
  const L = 1e-6, N = 8, D = 1e-9;
  const dImmob = new Device({
    species: [{ name: 'X', z: 0, cRef: 1 }],
    materials: { d: { epsr: 0, species: { X: { D: 0, mu0: 0 } } } },
    regions: [{ material: 'd', length: L, c0: { X: 1 } }],
    grid: { minCells: 4, hmin: L / 4, hmax: L / 4 },
  });
  dImmob.solve();
  const rImmob = dImmob.stability();
  assert.equal(rImmob.verdict, 'no growing mode found');
  assert.equal(rImmob.modes.length, 0);
  assert.equal(rImmob.deflated.length, 1);
  assert.equal(rImmob.deflated[0].kind, 'local');
  assert.equal(rImmob.deflated[0].count, 5);

  const el = { terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
  const dSalt = new Device({
    species: [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }],
    materials: { water: { epsr: 0, species: { 'Ag+': { D, mu0: 77.1e3 }, 'NO3-': { D, mu0: -111.3e3 } } } },
    regions: [{ material: 'water', length: L, c0: { 'NO3-': 10 } }],
    grid: { minCells: N, hmin: L / N, hmax: L / N },
    contacts: { left: { ...el, V: 0 }, right: { ...el, I: 0 } },
  });
  dSalt.solve();
  const rSalt = dSalt.stability();
  assert.equal(rSalt.verdict, 'no growing mode found');
  assert.equal(rSalt.deflated.length, 1);
  const h = L / N;
  for (let k = 1; k <= 3; k++) {
    const s = Math.sin((k * Math.PI) / (2 * N));
    const exact = (-4 * D * s * s) / (h * h);
    assert.ok(Math.abs(rSalt.modes[k - 1].re - exact) / Math.abs(exact) < 1e-8);
  }
});

// ============================================================================
// 4. Acceptance Gate 4: Fe³⁺/Fe²⁺ d26594 falling vs rising branch vs advance()
// ============================================================================

test('Fe3+/Fe2+ d26594: rising branch & V-held falling branch have no growing mode; I-held falling branch is unstable and matches kicked advance()', () => {
  // 1. Rising branch under I drive
  const dRise = new Device(feSpec());
  dRise.solve();
  const rRise = dRise.stability();
  assert.equal(rRise.verdict, 'no growing mode found');
  assert.equal(rRise.deflated.length, 3);
  assert.ok(rRise.modes[0].re < -700 && rRise.modes[0].re > -740, `rising slowest=${rRise.modes[0].re}`);

  // 2. Falling branch under held V (-2.03569 V)
  const dFall = new Device(feSpec({ V: -2.03569 }));
  dFall.solve();
  const rFallV = dFall.stability();
  assert.equal(rFallV.verdict, 'no growing mode found');
  assert.ok(rFallV.modes[0].re < -1000 && rFallV.modes[0].re > -1100, `fall-V slowest=${rFallV.modes[0].re}`);

  // 3. Switch to I drive on the falling branch (preserves branch!)
  dFall.set({ contacts: { right: { V: undefined, I: -62.90054296647188 } } });
  const solFallI = dFall.solve();
  assert.equal(solFallI.converged, true);
  assert.ok(solFallI.terminals.right.V < -2.0, `expected falling branch V<-2, got ${solFallI.terminals.right.V}`);

  const rFallI = dFall.stability();
  assert.equal(rFallI.verdict, 'unstable');
  assert.equal(rFallI.deflated.length, 3);
  const lam = rFallI.modes[0];
  assert.ok(lam.re > 149 && lam.re < 152, `expected ~+150.65 /s, got ${lam.re}`);
  assert.equal(lam.im, 0);
  assert.ok(lam.re > 10 * lam.error);
  assert.ok(lam.residual < 1e-8);
  // Second mode must be genuinely negative (no duplicate of +150.65)
  assert.ok(rFallI.modes[1].re < -800, `expected second mode near -853 /s, got ${rFallI.modes[1].re}`);

  // 4. Independent kicked advance() check of exponential growth rate
  const V0 = solFallI.terminals.right.V;
  dFall.set({ contacts: { right: { I: -62.90054296647188 + 1e-3 } } });
  const kick = dFall.step(1e-4);
  assert.equal(kick.converged, true);
  dFall.set({ contacts: { right: { I: -62.90054296647188 } } });
  dFall.advance(1e-2, { method: 'bdf2', dtMax: 1e-4, tol: 1e-9 });
  const dv1 = dFall.solution().terminals.right.V - V0;
  dFall.advance(2e-2, { method: 'bdf2', dtMax: 1e-4, tol: 1e-9 });
  const dv2 = dFall.solution().terminals.right.V - V0;
  assert.ok(dv1 * dv2 > 0 && Math.abs(dv2) > Math.abs(dv1), 'expected monotonic exponential departure');
  const rate = Math.log(dv2 / dv1) / (2e-2 - 1e-2);
  assert.ok(Math.abs(rate - lam.re) / lam.re < 0.005, `advance growth rate ${rate} vs pencil ${lam.re}`);
});

// ============================================================================
// 5. Contract: refusal, non-steady undecided, bit-for-bit state preservation,
//    profiles, debug hook, waveform warning, determinism
// ============================================================================

test('stability() refuses non-steady / unsolved / regridded states and invalid options', () => {
  const L = 1e-6;
  const mk = () =>
    new Device({
      species: [{ name: 'X', z: 0, cRef: 1 }],
      materials: { d: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } },
      regions: [{ material: 'd', length: L, c0: { X: 1 } }],
      grid: { minCells: 8, hmin: L / 8, hmax: L / 8 },
    });

  // Unsolved
  const d0 = mk();
  assert.throws(() => d0.stability(), SolverError);

  // Invalid options
  d0.solve();
  assert.throws(() => d0.stability({ shifts: [0] }), SolverError);
  assert.throws(() => d0.stability({ shifts: [-1] }), SolverError);
  assert.throws(() => d0.stability({ frequencies: [0] }), SolverError);
  assert.throws(() => d0.stability({ modes: 0 }), SolverError);
  assert.throws(() => d0.stability({ bogus: true }), /no option 'bogus'/);

  // After step()
  d0.step(1e-6);
  assert.throws(() => d0.stability(), SolverError);

  // After regrid (regriddedSteady)
  const d1 = mk();
  d1.solve();
  d1.set({ grid: { minCells: 12, hmin: L / 12, hmax: L / 12 } });
  assert.throws(() => d1.stability(), SolverError);
});

test('stability() returns undecided when state drifts under an un-kicked dynamic step', () => {
  const L = 1e-6;
  const dev = new Device({
    species: [{ name: 'X', z: 0, cRef: 1 }],
    materials: { d: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } },
    regions: [{ material: 'd', length: L, c0: { X: 1 } }],
    grid: { minCells: 8, hmin: L / 8, hmax: L / 8 },
  });
  dev.solve();
  // Artificially perturb state while keeping atSteady=true to simulate a non-steady branch
  dev.solver.u[1] += 0.2;
  const res = dev.stability();
  assert.equal(res.verdict, 'undecided');
  assert.ok(res.drift > 1e-5, `expected drift > 1e-5, got ${res.drift}`);
  assert.equal(res.modes.length, 0);
});

test('stability() is completely non-destructive, deterministic, emits profiles/debug, and warns only on waveforms', () => {
  const twin1 = hhDevice(9.8);
  const twin2 = hhDevice(9.8);
  twin1.solve();
  twin2.solve();

  const learnedBefore = {
    lostMode: twin1.solver.lostMode,
    accepted: twin1.solver.accepted,
    firstErrs: twin1.solver.firstErrs?.slice(),
  };

  const rA = twin1.stability({ profiles: true, debug: true });
  const rB = twin1.stability({ profiles: true });

  // Deterministic across repeated calls
  assert.equal(rA.verdict, rB.verdict);
  assert.equal(rA.modes.length, rB.modes.length);
  for (let i = 0; i < rA.modes.length; i++) {
    assert.equal(rA.modes[i].re, rB.modes[i].re);
    assert.equal(rA.modes[i].im, rB.modes[i].im);
    assert.deepEqual(rA.modes[i].profiles.phi.re, rB.modes[i].profiles.phi.re);
    assert.deepEqual(rA.modes[i].profiles.phi.im, rB.modes[i].profiles.phi.im);
  }
  assert.ok(rA.debug && Array.isArray(rA.debug.J) && Array.isArray(rA.debug.M));
  assert.equal(rA.debug.J.length, rA.debug.M.length);
  assert.equal(rA.warnings.length, 0, 'constant drive must not trigger waveform warning');

  // LEARNED preserved
  assert.equal(twin1.solver.lostMode, learnedBefore.lostMode);
  assert.equal(twin1.solver.accepted, learnedBefore.accepted);
  assert.deepEqual(twin1.solver.firstErrs, learnedBefore.firstErrs);

  // Bit-identical subsequent transient vs twin that never called stability()
  twin1.set({ contacts: { right: { I: 9.81e-2 } } });
  twin2.set({ contacts: { right: { I: 9.81e-2 } } });
  const s1 = twin1.step(1e-4);
  const s2 = twin2.step(1e-4);
  assert.deepEqual(s1.phi, s2.phi);
  assert.deepEqual(s1.terminals, s2.terminals);

  // Waveform drive triggers non-autonomous warning
  const L = 1e-6, cgeo = (2 * EPS0) / L;
  const dw = new Device({
    species: [{ name: 'X', z: 0, cRef: 1 }],
    materials: { d: { epsr: 2, species: { X: { D: 0, mu0: 0 } } } },
    regions: [{ material: 'd', length: L, c0: { X: 1 } }],
    grid: { minCells: 4, hmin: L / 4, hmax: L / 4 },
    contacts: {
      left: { V: 0, phi: 'pinned', zeroCharge: 0 },
      right: { V: pulse({ start: 1e-3, width: 1e-3, amplitude: 0.1 }), R: 5, phi: { type: 'capacitive', C: cgeo }, zeroCharge: 0 },
    },
  });
  dw.solve();
  const rw = dw.stability();
  assert.ok(rw.warnings.some((w) => w.includes('non-autonomous')), `expected waveform warning, got ${JSON.stringify(rw.warnings)}`);
});
