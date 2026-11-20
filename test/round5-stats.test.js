import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, AVOGADRO } from '../src/index.js';
import { normalizeStatistics } from '../src/statistics.js';

// Faults found by round 5's outside reading pass on the statistics (Sol's fifth pass), each
// confirmed by a probe against a second route before it was fixed.

const T = 298.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;
const logit = (x) => Math.log(x / (1 - x));

// A material's models built straight from its statistics (species [{ name, z }], sp: { name: cRef }).
function models(species, sp, statistics, epsr = 0) {
  const n = species.length, speciesIndex = new Map(species.map((s, i) => [s.name, i]));
  const present = new Uint8Array(n), cRef = new Float64Array(n), mu0 = new Float64Array(n);
  for (const [nm, c] of Object.entries(sp)) {
    const i = speciesIndex.get(nm);
    present[i] = 1;
    cRef[i] = c;
  }
  const need = (ok, msg) => {
    if (!ok) throw new Error(msg);
  };
  const h = { need, finite: (v, p) => (need(Number.isFinite(v), p), v), positive: (v, p) => (need(v > 0, p), v), dense: (a) => a };
  return normalizeStatistics({ epsr, present, cRef, mu0, statistics }, 'm', species, speciesIndex, RT, h);
}
const evaluate = (md, zeta) => {
  const k = md.idx.length, c = new Float64Array(k), K = new Float64Array(k * k);
  md.evaluate(Float64Array.from(zeta), c, K, 0);
  return { c, K };
};
// ζ of the first species at which it holds c, the others' ζ at 0
const zetaAt = (md, c) => {
  const k = md.idx.length, z = new Float64Array(k), fixed = new Uint8Array(k), t = new Float64Array(k);
  fixed[0] = 1;
  t[0] = c;
  md.invert(z, fixed, t, 0);
  return z[0];
};
const lithium = [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }];

// The inversion's bracket was set from r(x) sampled on a grid of 4001 x's, which steps over a knot
// at x = 1e-8 whose residual stands 1.1 above the rest: x(ζ) there was clipped at the bracket's
// edge, e^0.1 = 1.105 times the table's filling, and a device held at that knot's E was 2.57 mV off
// its own table, with no warning. The bracket now holds every knot's r.
test('an OCV table gives its own filling at and near a knot the sampling grid misses', () => {
  const x = [1e-12, 1e-8, 1e-6, 0.9], r = [0, 1.1, 0, 0];
  const E = x.map((v, i) => -VT * (logit(v) + r[i]));
  const st = [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }];
  const md = models(lithium, { 'Li+': 30000, 'e-': 30000 }, st).models[0];
  assert.ok(Math.abs(evaluate(md, [logit(1e-8) + 1.1, 0]).c[0] / 30000 / 1e-8 - 1) < 1e-12);
  for (const xx of [3e-9, 1e-8, 3e-8, 1e-7]) {
    const back = evaluate(md, [zetaAt(md, xx * 30000), 0]).c[0] / 30000;
    assert.ok(Math.abs(back / xx - 1) < 1e-12, `x = ${xx}: x(ζ(x)) = ${back}`);
  }
  const sol = new Device({
    species: lithium,
    materials: { host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } }, statistics: st } },
    regions: [{ material: 'host', length: 1e-6 }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { V: E[1], terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
    grid: { minCells: 4 },
  }).solve();
  for (const c of sol.c['Li+']) assert.ok(Math.abs(c / 30000 / 1e-8 - 1) < 1e-12, `x = ${c / 30000}`);
});

// Knots 1.6e-4 apart, E strictly falling at every one, passed the monotonicity check (the grid's
// points are 2.5e-4 apart), while the interpolated ζ(x) turned back across a third of the interval
// between them (dζ/dy = −0.125 at its middle, by hand and by the model): x(ζ(x)) was 7e-5 off.
test('an OCV table whose curve turns back between two close knots is refused', () => {
  const x = [0.3, 0.500045, 0.500205, 0.7], r = [0, 0, -0.00048, -0.00048];
  const E = x.map((v, i) => -VT * (logit(v) + r[i]));
  assert.ok(E.every((e, i) => i === 0 || e < E[i - 1]));
  const st = (rr) => [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 1, ocv: { x, E: x.map((v, i) => -VT * (logit(v) + rr[i])), muRef: 0 } }];
  assert.throws(() => models(lithium, { 'Li+': 1, 'e-': 1 }, st(r)), /ocv: near x = 0\.5001\d* the curve is flat or rising/);
  // a tenth of that step turns nothing back, and passes
  assert.doesNotThrow(() => models(lithium, { 'Li+': 1, 'e-': 1 }, st(r.map((v) => v / 10))));
});

// For an ion of charge z, μ_ion + ν μ_carrier = muRef − zFE rises with x: E rises for an anion.
// An anion host's physical, ideal-solution table (z = −2, a hole carrier) was refused for not
// falling, and the same table reversed for its slope.
test('an OCV table for an anion host rises with x, and is read as the isotherm', () => {
  const sp = [{ name: 'A2-', z: -2 }, { name: 'h+', z: 1 }], x = [0.1, 0.5, 0.9];
  const E = x.map((v) => (VT / 2) * logit(v));
  const st = (EE) => [{ type: 'insertion', species: ['A2-', 'h+'], cMax: 1, ocv: { x, E: EE, muRef: 0 } }];
  const md = models(sp, { 'A2-': 1, 'h+': 1 }, st(E)).models[0];
  for (const xx of [0.1, 0.3, 0.5, 0.9]) assert.ok(Math.abs(zetaAt(md, xx) - logit(xx)) < 1e-12, `x = ${xx}`);
  assert.throws(() => models(sp, { 'A2-': 1, 'h+': 1 }, st(E.map((e) => -e))), /ocv\.E must be strictly increasing with x for an ion of charge -2/);
});

// K = diag(c) − c cᵀ/c_max: at a saturated species, c_A − c_A²/c_max rounded to 0 (and K's
// determinant went negative). By hand, w = e^ζ, S = 1 + Σw: K_AA = c_A(1 + w_B)/S = 4.277e-18.
test('a lattice keeps the vacancy\'s share of a saturated species\' capacitance', () => {
  const md = models([{ name: 'A', z: 1 }, { name: 'B', z: -1 }], { A: 1, B: 1 }, [{ type: 'lattice', species: ['A', 'B'], cMax: 1 }], 1).models[0];
  const zeta = [40, -5], w = zeta.map(Math.exp), S = 1 + w[0] + w[1], c = w.map((v) => v / S);
  const { K } = evaluate(md, zeta);
  const want = [(c[0] * (1 + w[1])) / S, -c[0] * c[1], -c[0] * c[1], (c[1] * (1 + w[0])) / S];
  K.forEach((v, j) => assert.ok(Math.abs(v / want[j] - 1) < 1e-12, `K[${j}] = ${v}, by hand ${want[j]}`));
  // and where nothing saturates, still c_a(1 − c_a/c_max) on the diagonal
  const mid = evaluate(md, [0.3, -0.7]);
  for (let a = 0; a < 2; a++) assert.ok(Math.abs(mid.K[a * 3] / (mid.c[a] * (1 - mid.c[a])) - 1) < 1e-14);
});

// A bar held 40 RT up holds c = c_max exactly (to round-off), and a set() of its μ° then refused
// to carry it ("exceeds the site density", "outside (0, c_max)"), though nothing was wrong.
test('set() carries a lattice or a Redlich–Kister population that rounds to full', () => {
  for (const st of [{ type: 'lattice', species: ['X'], cMax: 1000 }, { type: 'redlich-kister', species: ['X'], cMax: 1000, A: [] }]) {
    for (const lev of [37, 40]) {
      const held = { phi: 'neutral', species: { X: { type: 'equilibrium', mu: lev * RT } } };
      const dev = new Device({
        T,
        species: [{ name: 'X', z: 0 }],
        materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1000, mu0: 0 } }, statistics: [st] } },
        regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 500 } }],
        contacts: { left: held, right: held },
        grid: { minCells: 8 },
      });
      assert.equal(dev.solve().c.X[3], 1000);
      dev.set({ materials: { f: { species: { X: { mu0: 10 } } } } });
      const s = dev.solution();
      for (const c of s.c.X) assert.ok(Math.abs(c / 1000 - 1) < 1e-15, `${st.type} at ${lev} RT: c = ${c}`);
      assert.ok(dev.solve().converged);
    }
  }
  // still refused where a carried amount does exceed the sites
  const md = models([{ name: 'X', z: 0 }], { X: 1 }, [{ type: 'lattice', species: ['X'], cMax: 1 }]).models[0];
  assert.throws(() => zetaAt(md, 1 + 1e-12), /exceed the site density/);
});

// In Debye–Hückel theory κ's sum runs over every ion in the solution, but a model sees only its
// own species. 10 mM NaCl + 90 mM KNO₃ as one group gives NaCl's (μ̄₊ + μ̄₋)/F at the textbook
// value, −249.629 mV at I = 0.1 M; as two groups, or with KNO₃ ideal, −241.890 mV (I = 0.01 M,
// NaCl's alone), 7.74 mV off, without a word. Those now warn.
test('Debye–Hückel warns where its model leaves out an ion of the solution', () => {
  const ions = [['Na+', 1], ['Cl-', -1], ['K+', 1], ['NO3-', -1]];
  const c0 = { 'Na+': 10, 'Cl-': 10, 'K+': 90, 'NO3-': 90 }, a = 4.5e-10, epsr = 78.5;
  const run = (groups) =>
    new Device({
      species: ions.map(([name, z]) => ({ name, z })),
      materials: { w: { epsr: 0, species: Object.fromEntries(ions.map(([nm]) => [nm, { D: 1e-9, mu0: 0, cRef: 1000 }])), statistics: groups.map((g) => ({ type: 'debye-huckel', species: g, epsr, a })) } },
      regions: [{ name: 'a', material: 'w', length: 1e-6, c0 }],
      contacts: { left: { bath: { c: c0, reference: 'Cl-' } }, right: { bath: { c: c0, reference: 'Cl-' } } },
      grid: { minCells: 4 },
    }).solve();
  const eps = epsr * EPS0, lB = (FARADAY * FARADAY) / (4 * Math.PI * eps * AVOGADRO * RT);
  const k = Math.sqrt((FARADAY * FARADAY * 200) / (eps * RT)), hand = 2 * VT * (Math.log(10 / 1000) - (lB * k) / (2 * (1 + k * a)));
  const one = run([['Na+', 'Cl-', 'K+', 'NO3-']]);
  assert.ok(Math.abs((one.mu['Na+'][2] + one.mu['Cl-'][2]) / FARADAY - hand) < 1e-9, `${hand}`);
  assert.deepEqual(one.warnings, []);
  const two = run([['Na+', 'Cl-'], ['K+', 'NO3-']]);
  assert.equal(two.warnings.filter((w) => /ionic strength runs over every ion in the solution/.test(w)).length, 2);
  assert.ok(two.warnings.some((w) => /statistics\[0\].*'K\+', 'NO3-' are charged, in this material, and left out of it/.test(w)), two.warnings.join('\n'));
  const spectator = run([['Na+', 'Cl-']]);
  assert.ok(spectator.warnings.some((w) => /statistics\[0\].*left out of it: list every ion of the solution in one debye-huckel model/.test(w)), spectator.warnings.join('\n'));
});

// The cold start balanced a held reacting port spread mid-film against the material at the
// region's start: in a region graded in the cation's c_ref (2D Fermi–Dirac, 1000 → 2000 mol/m³) it
// started 0.595 RT from equilibrium (ln(e^(1000/1500) − 1) − ln(e − 1), by hand), the port passing
// 9.5e-4 A cold and a made-up current in the first transient. The material is now read where the
// port is, as the uniform film's (and its mirror's) cold state is: at rest.
test('a held reacting port in a graded region starts at its own equilibrium', () => {
  const L = 1e-6;
  for (const values of [[1000, 2000], [2000, 1000]]) {
    const dev = new Device({
      T,
      species: [{ name: 'M+', z: 1 }, { name: 'A-', z: -1 }, { name: 'e-', z: -1 }],
      materials: { f: { epsr: 0, species: { 'M+': { D: 1e-9, cRef: 1500, mu0: 0 }, 'A-': { D: 1e-9, cRef: 1000, mu0: 0 } }, statistics: [{ type: 'fermi-dirac', order: 0, species: ['M+'] }] } },
      regions: [{ name: 'a', material: 'f', length: L, c0: { 'M+': 1000, 'A-': 1000 }, grade: { 'M+': { cRef: { x: [0, L], values } } } }],
      contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
      ports: [{ name: 'el', region: 'a', from: 0.49e-6, to: 0.51e-6, V: 0, terminal: 'e-', area: 1e4, reactions: [{ equation: 'M+ + e- = M(s)', fixed: { 'M(s)': 0 }, k0: 1e-4, alpha: 0.5 }] }],
      grid: { minCells: 40 },
    });
    const cold = dev.solution();
    const g = cold.x.reduce((b, xx, q) => (Math.abs(xx - L / 2) < Math.abs(cold.x[b] - L / 2) ? q : b), 0);
    assert.ok(Math.abs(cold.mu['M+'][g] / RT) < 1e-9, `${values}: μ̄(M+)/RT = ${cold.mu['M+'][g] / RT} cold`);
    assert.ok(Math.abs(cold.ports[0].current) < 1e-12, `${values}: cold current ${cold.ports[0].current}`);
  }
});
