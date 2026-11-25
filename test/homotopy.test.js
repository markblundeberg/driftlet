import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Device, FARADAY, GAS_CONSTANT, AVOGADRO, units } from '../src/index.js';
import { build, aqueous, metal, semiconductor, layer, ohmic, bath, half, photogeneration } from '../src/kit.js';
import { computeHomotopyDs, solveBorderedArclength } from '../src/homotopy.js';
import { equilibrate } from '../src/newton.js';
import { LEARNED } from '../src/solver.js';

function withLadder(fn) {
  const prev = globalThis.DRIFTLET_HOMOTOPY;
  globalThis.DRIFTLET_HOMOTOPY = false;
  try {
    return fn();
  } finally {
    if (prev === undefined) delete globalThis.DRIFTLET_HOMOTOPY;
    else globalThis.DRIFTLET_HOMOTOPY = prev;
  }
}

test('computeHomotopyDs matches central finite differences for held V, circuitHomotopy and generation', () => {
  const Si = semiconductor('Si'), NA = units.perCm3(1e16), W = units.um(60.5);
  const dev = new Device(build({
    T: 300,
    library: [Si],
    stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e19) }), layer('Si', units.um(60), { acceptors: NA }), { ...ohmic(0), V: undefined, I: 0 }],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1 / (60e-6 * NA) } }, photogeneration({ material: 'Si', flux: 3e-5, alpha: 1630, mu: units.eV(1.18), to: W })],
    grid: { hmin: units.nm(1), hmax: units.um(1) },
  }));
  const s = dev.solver;
  s.steady = true;
  s.allTerminals = true;
  s.homotopyGen = new Float64Array(s.sys.size + 1);
  s.homotopyGenTerm = new Float64Array(s.terms.length);
  s.circuitHomotopy = new Array(s.terms.length).fill(null);

  const plan = {
    p0: Float64Array.from([0.1, 0]),
    p1: Float64Array.from([0.4, 0]),
    ch: [null, { a0: 0, da: 1, b0: 5, db: -5, v0: 0.1, dv: 0.2, g0: 0, dg: 12 }],
    rampGen: true,
  };
  const GEN_LOG_SPAN = 30 * Math.LN10;
  const setS = (sv) => {
    s.sourceOverride.set(0, plan.p0[0] + sv * (plan.p1[0] - plan.p0[0]));
    const c = plan.ch[1];
    s.circuitHomotopy[1] = { alpha: c.a0 + sv * c.da, beta: c.b0 + sv * c.db, v0: c.v0 + sv * c.dv, gamma: c.g0 + sv * c.dg };
    s.generationScale = Math.exp(GEN_LOG_SPAN * (sv - 1));
  };
  const s0 = 0.85, eps = 1e-7;
  s.termV[1] = 0.25;
  setS(s0);
  s.assemble(Infinity);
  equilibrate(s);
  const rowDiv = Float64Array.from(s.rowDiv);
  const Bs = new Float64Array(s.sys.size + 1);
  const pinDs = new Float64Array(Math.max(1, s.pins.length));
  const circDs = new Float64Array(s.floating.length);
  computeHomotopyDs(s, plan, s0, Bs, pinDs, circDs);

  setS(s0 + eps);
  s.assemble(Infinity);
  const rPlus = Float64Array.from(s.res), cPlus = s.termRes[1];
  setS(s0 - eps);
  s.assemble(Infinity);
  const rMinus = Float64Array.from(s.res), cMinus = s.termRes[1];

  let maxErrB = 0;
  for (let j = 0; j < s.sys.size; j++) {
    const fd = (rPlus[j] - rMinus[j]) / (2 * eps * rowDiv[j]);
    const scale = Math.max(1, Math.abs(Bs[j]), Math.abs(fd));
    maxErrB = Math.max(maxErrB, Math.abs(Bs[j] - fd) / scale);
  }
  const fdC = (cPlus - cMinus) / (2 * eps);
  const errC = Math.abs(circDs[0] - fdC) / Math.max(1, Math.abs(circDs[0]), Math.abs(fdC));
  assert.ok(maxErrB < 1e-7, `Bs FD error ${maxErrB}`);
  assert.ok(errC < 1e-7, `circDs FD error ${errC}`);
});

test('solveBorderedArclength satisfies the full augmented linear system (pins + floating + arclength)', () => {
  // Electrolyte with a blocked spectator (pin) and a current-driven contact (floating terminal).
  const end = { terminal: 'K+', species: { 'K+': 'equilibrium', 'Cl-': 'blocked' }, phi: 'neutral' };
  const dev = new Device({
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { w: { epsr: 0, species: { 'K+': { D: 1e-9, cRef: 1, mu0: 0 }, 'Cl-': { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ material: 'w', length: 1e-5, c0: { 'K+': 10, 'Cl-': 10 } }],
    contacts: { left: { V: 0, ...end }, right: { I: 1, ...end } },
    grid: { minCells: 10 },
  });
  const s = dev.solver;
  s.steady = true;
  s.allTerminals = true;
  s.assemble(Infinity);
  equilibrate(s);
  const N = s.sys.size, P = s.pins.length, K = s.floating.length;
  assert.ok(P >= 1 && K === 1, `expected pins >= 1 and K = 1, got P=${P}, K=${K}`);

  const Bs = new Float64Array(N + 1), wx = new Float64Array(N + 1);
  const rhs = new Float64Array(N + 1), pinRhs = new Float64Array(P), circRhs = new Float64Array(K);
  const pinDs = new Float64Array(P), circDs = new Float64Array(K), wV = new Float64Array(K);
  const pinnedRows = new Set();
  for (const p of s.pins) for (const o of p.col ?? [p.row]) pinnedRows.add(s.rix[o]);
  for (let j = 0; j < N; j++) {
    Bs[j] = pinnedRows.has(j) ? 0 : Math.sin(j + 1);
    wx[j] = Math.cos(j + 1) / N;
    rhs[j] = pinnedRows.has(j) ? 0 : Math.sin(2 * j + 3);
  }
  for (let q = 0; q < P; q++) { pinRhs[q] = 0.3 + q; pinDs[q] = -0.2 * (q + 1); }
  circRhs[0] = 0.7;
  circDs[0] = -1.5;
  wV[0] = 0.4;
  const ws = 2.5, gRhs = 1.1;

  // Save unfactored matrix-vector product via sys.matvec before factoring.
  const matvec = (v) => {
    const out = new Float64Array(N);
    const { A, B, C, sizes, offA, offB, offC, offX } = s.sys, nBlocks = sizes.length;
    for (let b = 0; b < nBlocks; b++) {
      const m = sizes[b], r0 = offX[b];
      for (let r = 0; r < m; r++) {
        let t = 0;
        for (let c = 0; c < m; c++) t += B[offB[b] + r * m + c] * v[r0 + c];
        if (b > 0) for (let c = 0, mp = sizes[b - 1], c0 = offX[b - 1]; c < mp; c++) t += A[offA[b] + r * mp + c] * v[c0 + c];
        if (b < nBlocks - 1) for (let c = 0, mn = sizes[b + 1], c0 = offX[b + 1]; c < mn; c++) t += C[offC[b] + r * mn + c] * v[c0 + c];
        out[r0 + r] = t;
      }
    }
    return out;
  };

  // Compute exact unfactored matvec after saving A, B, C copies.
  const saved = { A: Float64Array.from(s.sys.A), B: Float64Array.from(s.sys.B), C: Float64Array.from(s.sys.C) };
  s._factor();
  const delta = new Float64Array(N + 1), deltaV = new Float64Array(K);
  const ds = solveBorderedArclength(s, rhs, pinRhs, circRhs, gRhs, delta, deltaV, Bs, pinDs, circDs, wx, wV, ws);

  // Restore unfactored blocks and check every row of the augmented system.
  s.sys.A.set(saved.A);
  s.sys.B.set(saved.B);
  s.sys.C.set(saved.C);
  const Tv = matvec(delta);
  const fk = s.floating[0];
  let maxBlockErr = 0;
  for (let j = 0; j < N; j++) {
    if (pinnedRows.has(j)) continue;
    const lhs = Tv[j] + s.termB[fk][j] * deltaV[0] + Bs[j] * ds;
    maxBlockErr = Math.max(maxBlockErr, Math.abs(lhs - rhs[j]));
  }
  assert.ok(maxBlockErr < 1e-10, `block residual ${maxBlockErr}`);
  for (let q = 0; q < P; q++) {
    const p = s.pins[q];
    let lhs = pinDs[q] * ds;
    for (let j = 0; j < p.len; j++) lhs += p.w[j] * delta[p.idx[j]];
    assert.ok(Math.abs(lhs - pinRhs[q]) < 1e-10, `pin ${q} residual ${Math.abs(lhs - pinRhs[q])}`);
  }
  let cLhs = s.termDI[fk] * deltaV[0] + circDs[0] * ds;
  for (let j = 0; j < N; j++) cLhs += s.termC[fk][j] * delta[j];
  assert.ok(Math.abs(cLhs - circRhs[0]) < 1e-10, `circuit residual ${Math.abs(cLhs - circRhs[0])}`);
  let gLhs = wV[0] * deltaV[0] + ws * ds;
  for (let j = 0; j < N; j++) gLhs += wx[j] * delta[j];
  assert.ok(Math.abs(gLhs - gRhs) < 1e-11, `arclength residual ${Math.abs(gLhs - gRhs)}`);
});

test('homotopy handles multi-decade R, immobile Q at open circuit, galvanic OCV behind R, and warm jumps', () => {
  const Si = semiconductor('Si'), NA = units.perCm3(1e16), W = units.um(60.5);
  // 1. Multi-decade series resistance (1e-6 to 1e9 Ω·m²) on a lit Si solar cell.
  let prevV = -Infinity;
  for (const R of [1e-6, 1, 1e3, 1e6, 1e9]) {
    const sol = new Device(build({
      T: 300, library: [Si],
      stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e19) }), layer('Si', units.um(60), { acceptors: NA }), ohmic({ V: 0, R })],
      bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1 / (60e-6 * NA) } }, photogeneration({ material: 'Si', flux: 3e-5, alpha: 1630, mu: units.eV(1.18), to: W })],
      grid: { hmin: units.nm(1), hmax: units.um(1) },
    })).solve();
    assert.ok(sol.converged, `R=${R} converged`);
    assert.ok(sol.terminals.right.V >= prevV - 1e-9, `V monotonic with R=${R}`);
    prevV = sol.terminals.right.V;
  }
  assert.ok(Math.abs(prevV - 0.38442194) < 1e-5, `high-R approaches V_oc: ${prevV}`);

  // 2. Lit Si cell at open circuit with immobile species Q: conserved node-by-node, matches ladder.
  const defQ = build({
    T: 300, library: [Si],
    stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e19) }), layer('Si', units.um(60), { acceptors: NA }), { ...ohmic(0), V: undefined, I: 0 }],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1 / (60e-6 * NA) } }, photogeneration({ material: 'Si', flux: 3e-5, alpha: 1630, mu: units.eV(1.18), to: W })],
    grid: { hmin: units.nm(0.25), hmax: units.um(1) },
  });
  defQ.species.push({ name: 'Q', z: 0 });
  defQ.materials.Si.species.Q = { D: 0, mu0: 0, cRef: 1000 };
  defQ.regions.forEach((r) => (r.c0 = { Q: 1 }));
  const solQ = new Device(defQ).solve();
  assert.ok(solQ.converged);
  assert.ok(Math.abs(solQ.terminals.right.V - 0.38442194) < 1e-5);
  assert.ok(solQ.c.Q.every((v) => Math.abs(v - 1) < 1e-12), 'immobile Q kept exact');

  // 3. Closed unlike cell (Cu | water | Zn) with Stern layers behind R at OCV.
  const RT = GAS_CONSTANT * 298.15;
  const c = { 'Cu2+': 10, 'Zn2+': 1.5, 'Na+': 100, 'Cl-': 123 };
  const lib = aqueous(['Cu2+', 'Zn2+', 'Na+', 'Cl-'], { epsr: 0 });
  const mu = (name) => lib.materials.water.species[name].mu0 + RT * Math.log(c[name] / 1000);
  const ocv = (mu('Zn2+') - mu('Cu2+')) / (2 * FARADAY), phi = -mu('Cu2+') / (2 * FARADAY);
  const solZn = new Device(build({
    library: [lib, metal('Cu'), metal('Zn')],
    stack: [
      ohmic(0, ['e-']), layer('Cu', 1e-6),
      { phi: { type: 'capacitive', C: 0.45 }, zeroCharge: -phi, reactions: [{ equation: 'Cu2+ + 2 e- = Cu(s)', fixed: { 'Cu(s)': 0 }, k0: 1e-3, alpha: 0.4 }] },
      layer('water', 12.6e-6, { c0: c }),
      { phi: { type: 'capacitive', C: 0.22 }, zeroCharge: ocv - phi, reactions: [{ equation: 'Zn2+ + 2 e- = Zn(s)', fixed: { 'Zn(s)': 0 }, k0: 7.6e-5, alpha: 0.49 }] },
      layer('Zn', 1e-6), ohmic({ V: ocv, R: 1e-4 }, ['e-']),
    ],
    grid: { hmin: 1.35e-8, hmax: 1.3e-7, ratio: 1.1 },
  })).solve();
  assert.ok(solZn.converged);
  assert.ok(Math.abs(solZn.terminals.right.V - ocv) < 1e-9 && Math.abs(solZn.current) < 1e-12);
});

test('homotopy refuses folds (Frumkin dip past peak) and model walls (Debye–Hückel) cleanly and fast', () => {
  // 1. Frumkin dip past its peak current: traces over the fold and refuses with unreached.
  const c = { FeIII: 1, FeII: 1, 'Cl-': 100, 'K+': 107 }, cMax = 1 / (AVOGADRO * 0.71e-9 ** 3);
  const unit = (FARADAY * 0.72e-9 * 1) / 30e-6;
  const over = new Device(build({
    library: [{ species: [{ name: 'FeIII', z: -3, cRef: 1000 }, { name: 'FeII', z: -4, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
      materials: { water: { epsr: 0, species: { FeIII: { D: 0.72e-9, mu0: 0 }, FeII: { D: 0.67e-9, mu0: -FARADAY * 0.36 }, 'K+': { D: 1.957e-9, mu0: 0 }, 'Cl-': { D: 2.032e-9, mu0: 0 } }, statistics: [{ type: 'lattice', species: ['FeIII', 'FeII', 'K+', 'Cl-'], cMax }] } } }, metal('Pt')],
    stack: [ohmic({ I: -0.006 * unit }, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.6, diffuse: { epsr: 78.4 }, reactions: [{ ...half('FeIII + e- = FeII'), k0: 5e-2, alpha: 0.5 }] }, layer('water', 30e-6, { c0: c }), bath(c, 'Cl-')],
    grid: { hmin: 30e-9, hmax: 2e-6, ratio: 1.2 },
  })).solve();
  assert.equal(over.converged, false);
  assert.match(over.warnings.join(' '), /no steady state passes the driven current/);

  // 2. Debye–Hückel fixture past non-convexity boundary: refuses cleanly in < 500 steps (vs ~8,357 on the ladder).
  const spec = JSON.parse(readFileSync(new URL('./fixtures/debye-huckel-past-range.json', import.meta.url), 'utf8'));
  const dh = new Device(spec).solve();
  assert.equal(dh.converged, false);
  assert.match(dh.warnings.join(' '), /debye-huckel: the activity model is not convex/);
  assert.ok(dh.steps < 500, `expected < 500 steps, got ${dh.steps}`);
});

test('homotopy honours budgetMs mid-trace on the Debye–Hückel fixture and restores STATE', (t) => {
  const spec = JSON.parse(readFileSync(new URL('./fixtures/debye-huckel-past-range.json', import.meta.url), 'utf8'));
  const d = new Device(spec), s = d.solver;
  const before = s._snapshot(), learned = LEARNED.map((k) => s[k]);
  const drives = s.terms.map((tm) => tm.drive);
  // Expire the clock at factorisation 100 (deep in correctorNewton's arclength loop, where no
  // _steadySteps check runs): only homotopy.js's per-step budget checks can stop it promptly.
  let elapsed = 0, factors = 0, afterExpiry = 0;
  const factor = s._factor.bind(s);
  t.mock.method(performance, 'now', () => elapsed);
  t.mock.method(s, '_factor', () => {
    if (elapsed) afterExpiry++;
    if (s.allTerminals && ++factors >= 100) elapsed = 50;
    return factor();
  });
  const sol = d.solve({ budgetMs: 10 });
  assert.equal(sol.converged, false);
  assert.equal(sol.budgetSpent, true);
  assert.ok(sol.iterations > 0);
  assert.ok(factors >= 100, `expected budget to fire deep in homotopy trace, got ${factors} factorisations`);
  assert.ok(afterExpiry <= 10, `homotopy kept running after budget expired (${afterExpiry} factorisations after expiry)`);
  assert.match(sol.warnings.join('\n'), /budget spent/);
  assert.doesNotMatch(sol.warnings.join('\n'), /debye-huckel: the activity model is not convex/);
  assert.deepEqual(s._snapshot(), before);
  assert.deepEqual(s.terms.map((tm) => tm.drive), drives);
  LEARNED.forEach((k, i) => assert.equal(s[k], learned[i], k));
  assert.ok(!s.allTerminals);
  assert.equal(s.circuitHomotopy, null);
});

test('homotopy leaves no spurious island warning on mirrored membrane21', () => {
  const T = 298.15, RT = GAS_CONSTANT * T;
  const c = (Na, A, NaA) => ({ 'Na+': Na, 'A-': A, 'Cl-': Na - A, NaA });
  const bathC = c(16.893326496828085, 15.47362164868072, 1.940883084130149);
  const def = {
    T,
    species: [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }, { name: 'A-', z: -1 }, { name: 'NaA', z: 0 }].map((sp) => ({ ...sp, cRef: 1000 })),
    materials: { water: { epsr: 0, species: { 'Na+': { D: 1.334e-9, mu0: -261905 }, 'Cl-': { D: 2.032e-9, mu0: -131228 }, 'A-': { D: 0.7e-9, mu0: -50000 }, NaA: { D: 0.8e-9, mu0: -316874.9194977297 } } } },
    regions: [
      { name: 'far', length: 20.47576830701132e-6, material: 'water', c0: c(51.01561905373098, 13.98690479103286, 58.948681717972356) },
      { name: 'near', length: 18.844268738599778e-6, material: 'water', c0: c(10.038076413435569, 4.308332373487203, 25.50203622444877) },
    ],
    contacts: { right: { V: -0.007051, bath: { c: bathC, reference: 'Cl-' } } },
    interfaces: [{ phi: 'neutral', species: { 'Na+': 'blocked' } }],
    bulkReactions: [{ equation: 'Na+ + A- = NaA', kf: { water: 0.005104701794038667 } }],
    grid: { hmin: 0.4e-6, hmax: 3e-6 },
  };
  const Kd = 1000 * Math.exp((-316874.9194977297 + 261905 + 50000) / RT);
  const want = (16.893326496828085 * 15.47362164868072) / Kd;
  const sol = new Device(def).solve();
  assert.equal(sol.converged, true);
  assert.deepEqual(sol.warnings, []);
  assert.ok(Math.abs(sol.c.NaA[0] / want - 1) < 1e-10, `NaA ${sol.c.NaA[0]} for ${want}`);
});

test('homotopy stays on the stable branch of an N-shaped I–V curve (d26594)', () => {
  const spec = {
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
      { c0: { 'Fe3+': 49.26617526505966, 'Fe2+': 2.8287313025525, 'Na+': 14.483483670977007, 'Cl-': 167.93947207126098 }, material: 'water', length: 3.029892868410713e-6 },
      { material: 'Pt', length: 1e-6 },
    ],
    interfaces: [
      { phi: { type: 'capacitive', C: 0.28032231191173196 }, zeroCharge: 0.18351268186233938, reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 0.044465444935077036, alpha: 0.4318889715708792 }] },
      { phi: { type: 'capacitive', C: 0.28032231191173196 }, zeroCharge: 0.18351268186233938, reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 0.044465444935077036, alpha: 0.4318889715708792 }] },
    ],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk', I: -62.90054296647188 },
    },
    grid: { hmin: 1.4284780382255565e-10, hmax: 6.391660743434045e-8, ratio: 1.1937939652428031 },
  };
  const ref = withLadder(() => new Device(structuredClone(spec)).solve());
  const hom = new Device(structuredClone(spec)).solve();
  assert.equal(ref.converged, true);
  assert.equal(hom.converged, true);
  assert.ok(Math.abs(hom.terminals.right.V - ref.terminals.right.V) < 1e-9, `expected V ≈ ${ref.terminals.right.V}, got ${hom.terminals.right.V}`);
  assert.ok(hom.terminals.right.V > -1.0 && hom.terminals.right.V < -0.9, `expected stable branch near -0.9276 V, got ${hom.terminals.right.V}`);
});

test('homotopy preserves Kirchhoff and species-order invariance on floating-island GaAs stack (d3152)', () => {
  const base = {
    T: 300,
    species: [{ name: 'h+', z: 1 }, { name: 'e-', z: -1 }],
    materials: {
      GaAs: {
        T: 300,
        epsr: 13.1,
        species: {
          'e-': { D: 0.021974199818470205, mu0: 0, cRef: 0.7804533615717079 },
          'h+': { D: 0.0010340799914574212, mu0: 137395.11294359344, cRef: 11.623773470216928 },
        },
      },
    },
    regions: [
      { name: 'L0', material: 'GaAs', length: 4.83135684487053e-7, fixedCharge: 63323.73389198883 },
      { name: 'L1', material: 'GaAs', length: 2.199261404093173e-7, fixedCharge: -862558.9713666966 },
      { name: 'L2', material: 'GaAs', length: 9.02424004639285e-7, fixedCharge: 337977.752668189 },
    ],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: { type: 'capacitive', C: 10 }, zeroCharge: 0.36519684155471577 },
      right: { V: -0.20181845263578002, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' },
    },
    bulkReactions: [],
    grid: { hmin: 4.3408686129575477e-10, hmax: 4.2208602152428585e-8, ratio: 1.1093198062386365 },
  };
  const rev = structuredClone(base);
  rev.species.reverse();

  const ref = withLadder(() => new Device(structuredClone(base)).solve());
  const fwd = new Device(structuredClone(base)).solve();
  const bwd = new Device(structuredClone(rev)).solve();

  assert.equal(ref.converged, true);
  assert.equal(fwd.converged, true);
  assert.equal(bwd.converged, true);
  assert.ok(!fwd.warnings.some((w) => w.includes("don't add up")), `fwd Kirchhoff warning: ${fwd.warnings.join('; ')}`);
  assert.ok(!bwd.warnings.some((w) => w.includes("don't add up")), `bwd Kirchhoff warning: ${bwd.warnings.join('; ')}`);
  const Iref = ref.terminals.right.current;
  for (const [label, sol] of [['fwd', fwd], ['bwd', bwd]]) {
    assert.ok(Math.abs(sol.terminals.left.current + sol.terminals.right.current) < 1e-20, `${label} Kirchhoff sum ${sol.terminals.left.current + sol.terminals.right.current}`);
    assert.ok(Math.abs(sol.terminals.right.current - Iref) < 1e-4 * Math.abs(Iref), `${label} I=${sol.terminals.right.current} vs ref ${Iref}`);
  }
});
