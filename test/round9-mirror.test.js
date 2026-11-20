// Devices that solved cold one way round and were refused as their mirror images (round 9's
// mirror property, from Codex's sweep of bench/combos, 20 Nov 2026: membrane cases 0 and 21, semi
// case 20, cut down). The mirrored equations have the same solution; what differed was the
// solver's left-to-right habits. Each test's comment gives the old numbers.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { check } from '../src/kit.js';
import { BlockTridiagonal } from '../src/blockTridiagonal.js';

// A bath feeding two strictly neutral compartments; free Na⁺ blocked at the face between them,
// its complex NaA passing, so the reaction ties the closed side's Na⁺ to the bath. With the bath
// on the right, the block Thomas elimination began in the closed compartment, whose φ nothing
// but the face beyond it holds (the face's unknowns are fluxes): a pivot ratio of 3e-40 where
// the other way round has 7e-4, Newton's first update 1e29 thermal units, and every stage failed
// (converged false, NaA at the closed end 58.9 for mass action's 1.9409).
const T = 298.15, RT = 8.31446261815324 * T;
const water = (species) => ({ epsr: 0, species });
const membrane21 = (mirror) => {
  const c = (Na, A, NaA) => ({ 'Na+': Na, 'A-': A, 'Cl-': Na - A, NaA });
  const bath = c(16.893326496828085, 15.47362164868072, 1.940883084130149);
  const reservoir = { V: -0.007051, bath: { c: bath, reference: 'Cl-' } };
  const regions = [
    { name: 'near', length: 18.844268738599778e-6, material: 'water', c0: c(10.038076413435569, 4.308332373487203, 25.50203622444877) },
    { name: 'far', length: 20.47576830701132e-6, material: 'water', c0: c(51.01561905373098, 13.98690479103286, 58.948681717972356) },
  ];
  return {
    T, species: [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }, { name: 'A-', z: -1 }, { name: 'NaA', z: 0 }].map((sp) => ({ ...sp, cRef: 1000 })),
    materials: { water: water({ 'Na+': { D: 1.334e-9, mu0: -261905 }, 'Cl-': { D: 2.032e-9, mu0: -131228 }, 'A-': { D: 0.7e-9, mu0: -50000 }, NaA: { D: 0.8e-9, mu0: -316874.9194977297 } }) },
    regions: mirror ? regions.toReversed() : regions,
    contacts: mirror ? { right: reservoir } : { left: reservoir },
    interfaces: [{ phi: 'neutral', species: { 'Na+': 'blocked' } }],
    bulkReactions: [{ equation: 'Na+ + A- = NaA', kf: { water: 0.005104701794038667 } }],
    grid: { hmin: 0.4e-6, hmax: 3e-6 },
  };
};

test('a closed compartment behind a face solves cold with its bath on either side (membrane 21)', () => {
  // Mass action at the bath's levels: NaA = Na⁺·A⁻/K_d, K_d = c° exp((μ°_NaA − μ°_Na − μ°_A)/RT).
  const Kd = 1000 * Math.exp((-316874.9194977297 + 261905 + 50000) / RT);
  const want = (16.893326496828085 * 15.47362164868072) / Kd;
  for (const mirror of [false, true]) {
    const d = new Device(membrane21(mirror)), s = d.solve();
    assert.equal(s.converged, true, `mirror ${mirror}: ${s.warnings.join('; ')}`);
    const end = mirror ? 0 : s.x.length - 1;
    assert.ok(Math.abs(s.c.NaA[end] / want - 1) < 1e-10, `mirror ${mirror}: NaA ${s.c.NaA[end]} for ${want}`);
    assert.ok(check(d, s, { refine: false }).ok);
  }
});

// Three compartments, two walls, three complexations (case 0 without its isomer): the same
// elimination from the closed end met two such levels, one behind each wall (pivot ratio 5e-26
// for 1.7e-3; converged false, NaA at the closed end 37.4 for 2.3795).
const membrane0 = (mirror) => {
  const names = ['K+', 'Na+', 'Cl-', 'A-', 'KA', 'NaA', 'B', 'BA-'], z = [1, 1, -1, -1, 0, 0, 0, -1];
  const D = [1.957, 1.334, 2.032, 0.7, 0.9, 0.8, 1, 0.6];
  const mu0 = [-283270, -261905, -131228, -50000, -340327.5511434275, -318929.61630775867, -5000, -61201.53132902209];
  const by = (v) => Object.fromEntries(names.map((n, i) => [n, v[i]]));
  const bath = by([5.732779930714038, 88.70109830030232, 92.77760045508795, 1.5772618518101975, 0.15584703216163354, 2.379536273779321, 4.105343877679436, 0.0790159241182222]);
  const regions = [
    [45.2756361147569, [44.9821503493127, 9.478800771838495, 4, 3.9769860862029183, 47.15523808317377, 14.38455007677451, 4.309757924763024, 46.48396503494828]],
    [38.10304485343702, [12.91349325040753, 45.43067482320767, 21.134298022216853, 28.433169764448387, 46.84334515664665, 5.910396791305683, 17.228841568581966, 8.776700286949955]],
    [24.903992679768364, [72.5285765252501, 11.627699340119122, 4, 46.59552858435215, 13.60976605575123, 37.42557766250286, 29.197236443430697, 33.56074728101707]],
  ].map(([um, c], i) => ({ name: `r${i}`, length: um * 1e-6, material: 'water', c0: by(c) }));
  const wall = (blocked) => ({ phi: 'neutral', species: Object.fromEntries(blocked.map((n) => [n, 'blocked'])) });
  const faces = [wall(['A-', 'NaA', 'B', 'BA-']), wall(['K+', 'B', 'BA-'])];
  const reservoir = { V: 0.01396, bath: { c: bath, reference: 'Cl-' } }, closed = wall(names);
  return {
    T, species: names.map((name, i) => ({ name, z: z[i], cRef: 1000 })),
    materials: { water: water(by(D.map((v, i) => ({ D: v * 1e-9, mu0: mu0[i] })))) },
    regions: mirror ? regions.toReversed() : regions,
    interfaces: mirror ? faces.toReversed() : faces,
    contacts: mirror ? { left: closed, right: reservoir } : { left: reservoir, right: closed },
    bulkReactions: ['K+ + A- = KA', 'Na+ + A- = NaA', 'B + A- = BA-'].map((equation) => ({ equation, kf: { water: 1.08032564139258e-5 } })),
    grid: { hmin: 0.4e-6, hmax: 3e-6 },
  };
};

test('three compartments behind two walls solve cold with the bath on either side (membrane 0)', () => {
  const s0 = new Device(membrane0(false)).solve(), d1 = new Device(membrane0(true)), s1 = d1.solve();
  assert.equal(s0.converged, true);
  assert.equal(s1.converged, true, s1.warnings.join('; '));
  const n = s0.x.length;
  for (const sp of ['NaA', 'BA-', 'K+', 'B']) {
    for (let g = 0; g < n; g++) {
      const a = s0.c[sp][g], b = s1.c[sp][n - 1 - g];
      assert.ok(Math.abs(a - b) <= 1e-9 * Math.abs(a) + 1e-15, `${sp} at node ${g}: ${a} and ${b}`);
    }
  }
  assert.ok(check(d1, s1, { refine: false }).ok);
});

// The block solver itself: a chain whose first block is singular alone (its second unknown held
// only through the next block, as a closed region's level is held only through the face beyond
// it), the whole system regular. Eliminated from the left, the pivot is perturbed and the solve
// is off; with bothWays it's done from the right, to round-off.
test('BlockTridiagonal factorises from the right where the left leaves a pivot of round-off', () => {
  const sizes = [2, 3, 1, 2];
  const build = (bothWays) => {
    const sys = new BlockTridiagonal(sizes.length, sizes);
    sys.staticPivots = true;
    sys.bothWays = bothWays;
    let seed = 7;
    const rand = () => ((seed = (seed * 16807) % 2147483647) / 2147483647) - 0.5;
    for (let i = 0; i < sys.n; i++) {
      const m = sizes[i], mp = i > 0 ? sizes[i - 1] : m, mn = i < sys.n - 1 ? sizes[i + 1] : m;
      for (let k = 0; k < m * mp; k++) if (i > 0) sys.A[sys.offA[i] + k] = rand();
      for (let k = 0; k < m * mn; k++) if (i < sys.n - 1) sys.C[sys.offC[i] + k] = rand();
      for (let k = 0; k < m * m; k++) sys.B[sys.offB[i] + k] = rand();
      for (let r = 0; r < m; r++) sys.B[sys.offB[i] + r * m + r] += 3 * m;
    }
    // Block 0: its second unknown in no row of its own block, nor its second row in any column.
    sys.B.fill(0, sys.offB[0] + 1, sys.offB[0] + 4);
    sys.B[sys.offB[0]] = 5;
    return sys;
  };
  const x = Float64Array.from({ length: 8 }, (_, k) => k + 1);
  const ref = build(false), b = ref.multiply(x);
  const err = (sys) => {
    sys.factor();
    const y = sys.solve(b);
    return Math.max(...y.map((v, k) => Math.abs(v - x[k])));
  };
  const one = build(false), two = build(true);
  const e1 = err(one), e2 = err(two);
  assert.equal(one.fromRight, false);
  assert.equal(two.fromRight, true);
  assert.ok(one.perturbedPivots > 0 && e1 > 1e-3, `from the left alone: error ${e1}`);
  assert.ok(e2 < 1e-12, `both ways: error ${e2}`);
  assert.equal(two.perturbedPivots, 0);
  assert.ok(two.minPivotRatio > 1e-6);
  // and the cancellation it reports is the kept factorisation's
  assert.ok(two.cancellation().digits < 6);
});
