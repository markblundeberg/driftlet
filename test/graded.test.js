import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';

// Graded regions: a species' μ° and cRef varying within a region, as profiles { x, values }
// against the device's x (regions[r].grade). The flux along a segment drifts down the
// whole change in zφ̂ + μ°/RT − ln cRef, so a level μ̄ still carries exactly nothing.
const T = 298.15, RT = GAS_CONSTANT * T;
const B = (x) => (Math.abs(x) < 1e-12 ? 1 - x / 2 : x / Math.expm1(x));

// One neutral species A between two contacts that hold its μ, in a region whose μ° rises
// linearly by dmu across it. In a linear potential the Scharfetter–Gummel flux is exact on any
// grid: N = (D/L)[B(Δ) c_L − B(−Δ) c_R], Δ = dmu/RT.
const D = 1e-9, L = 10e-6, cRef = 1000;
const neutral = ({ mu0, muL, muR, grade, grid = { hmin: 5e-9, hmax: 0.5e-6, ratio: 1.2 } }) => ({
  T,
  species: [{ name: 'A', z: 0, cRef }],
  materials: { m: { epsr: 0, species: { A: { D, mu0: 0 } } } },
  regions: [{ material: 'm', length: L, grade: grade ?? { A: { mu0: { x: [0, L], values: [0, mu0] } } } }],
  contacts: {
    left: { species: { A: { type: 'equilibrium', mu: muL } }, phi: 'neutral' },
    right: { species: { A: { type: 'equilibrium', mu: muR } }, phi: 'neutral' },
  },
  grid,
});
const exactFlux = (dmu, muL, muR) => {
  const d = dmu / RT, cL = cRef * Math.exp(muL / RT), cR = cRef * Math.exp((muR - dmu) / RT);
  return (D / L) * (B(d) * cL - B(-d) * cR);
};

test('graded μ°: a neutral species down a linear μ° carries the exact Scharfetter–Gummel flux, on a graded grid', () => {
  for (const [dmu, muL, muR] of [[5 * RT, 0, -RT], [-3 * RT, 0, RT], [20 * RT, 2 * RT, -RT], [0.01 * RT, 0, -0.5 * RT]]) {
    const sol = new Device(neutral({ mu0: dmu, muL, muR })).solve();
    assert.ok(sol.converged, sol.warnings.join('; '));
    const N = sol.contacts.left.flux.A, want = exactFlux(dmu, muL, muR);
    assert.ok(Math.abs(N / want - 1) < 1e-9, `Δμ° ${dmu / RT} RT: ${N} vs ${want}`);
    // and c at each node is cRef e^{(μ − μ°(x))/RT}, with μ°(x) the profile's
    const g = sol.x.length >> 1, mu0 = (dmu * sol.x[g]) / L;
    assert.ok(Math.abs(sol.muStd.A[g] - mu0) < 1e-9 * RT, `μ° at the middle: ${sol.muStd.A[g]} vs ${mu0}`);
    assert.ok(Math.abs(sol.c.A[g] / (cRef * Math.exp((sol.mu.A[g] - mu0) / RT)) - 1) < 1e-12);
  }
});

test('graded μ°: at a level μ the flux is exactly zero and c follows e^{−μ°(x)/RT}, whatever the profile', () => {
  const grade = { A: { mu0: { x: [0, 0.3 * L, 0.5 * L, L], values: [0, 8 * RT, -4 * RT, 2 * RT] }, cRef: { x: [0, L], values: [1000, 50] } } };
  const sol = new Device(neutral({ muL: RT, muR: RT, grade })).solve();
  assert.ok(sol.converged);
  assert.equal(sol.contacts.left.flux.A, 0);
  for (let g = 0; g < sol.x.length; g++) {
    const x = sol.x[g], mu0 = sol.muStd.A[g];
    const cr = 1000 + ((50 - 1000) * x) / L;
    assert.ok(Math.abs(sol.c.A[g] / (cr * Math.exp((RT - mu0) / RT)) - 1) < 1e-12, `x = ${x}`);
  }
});

test('graded μ°: the same flux as a staircase of uniform layers, in the limit of many', () => {
  // A staircase of n layers, each at its middle's μ°, joined by faces where A passes in
  // equilibrium: its flux converges on the graded region's (the second route, through faces).
  const dmu = 6 * RT, muL = 0, muR = -2 * RT;
  const graded = new Device(neutral({ mu0: dmu, muL, muR })).solve().contacts.left.flux.A;
  const errs = [4, 16, 64].map((k) => {
    const materials = {}, regions = [];
    for (let j = 0; j < k; j++) {
      materials[`m${j}`] = { epsr: 0, species: { A: { D, mu0: (dmu * (j + 0.5)) / k } } };
      regions.push({ material: `m${j}`, length: L / k });
    }
    const def = { ...neutral({ mu0: 0, muL, muR }), materials, regions, grid: { hmin: 5e-9, hmax: L / k / 4, ratio: 1.2 } };
    const sol = new Device(def).solve();
    assert.ok(sol.converged);
    return Math.abs(sol.contacts.left.flux.A / graded - 1);
  });
  assert.ok(errs[2] < 1e-3, `64 layers: ${errs[2]}`);
  assert.ok(errs[1] < errs[0] / 8 && errs[2] < errs[1] / 8, `second order in 1/n: ${errs}`); // (the midpoint rule's)
});

// A graded-gap semiconductor: e⁻ and h⁺ with μ° and cRef graded at different rates, so the band
// edges tilt by different amounts (the quasi-fields). At equilibrium the Fermi level is flat,
// nothing flows, and n p = n_i²(x) at every node.
const si = (grade, extra = {}) => ({
  T,
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: { s: { epsr: 11.7, species: { 'e-': { D: 3.6e-3, mu0: 0, cRef: 0.047 }, 'h+': { D: 1.2e-3, mu0: 108e3, cRef: 0.031 } } } },
  regions: [{ material: 's', length: 2e-6, fixedCharge: -FARADAY * 0.0166, grade, ...extra }],
  contacts: {
    left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' },
    right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' },
  },
  grid: { hmin: 1e-9, hmax: 20e-9 },
});
// (profiles are against the device's x: for a layer starting at x0, gapFrom(x0))
const gapFrom = (x0) => ({
  'e-': { mu0: { x: [x0, x0 + 2e-6], values: [0, -15e3] }, cRef: { x: [x0, x0 + 2e-6], values: [0.047, 0.03] } },
  'h+': { mu0: { x: [x0, x0 + 2e-6], values: [108e3, 120e3] } },
});
const gap = gapFrom(0);

test('graded gap: at equilibrium the Fermi level is flat, nothing flows, and n p = n_i²(x)', () => {
  const sol = new Device(si(gap)).solve();
  assert.ok(sol.converged, sol.warnings.join('; '));
  assert.ok(Math.abs(sol.current) < 1e-9, `current ${sol.current}`);
  for (let g = 0; g < sol.x.length; g++) {
    assert.ok(Math.abs(sol.V['e-'][g]) < 1e-12 && Math.abs(sol.V['h+'][g]) < 1e-12, `x = ${sol.x[g]}: ${sol.V['e-'][g]}, ${sol.V['h+'][g]}`);
    const f = sol.x[g] / 2e-6, ni2 = (0.047 + (0.03 - 0.047) * f) * 0.031 * Math.exp(-(108e3 + 12e3 * f - 15e3 * f) / RT);
    assert.ok(Math.abs((sol.c['e-'][g] * sol.c['h+'][g]) / ni2 - 1) < 1e-9, `n p at x = ${sol.x[g]}`);
  }
  // the two standard levels' tilts add to their gradings' sum, −15 + 12 kJ/mol: φ's share of each
  // (z_e⁻Fφ and z_h⁺Fφ) cancels between them
  const last = sol.x.length - 1;
  const tilt = (name) => sol.muStd[name][last] - sol.muStd[name][0];
  assert.ok(Math.abs(tilt('e-') + tilt('h+') - (-15e3 + 12e3)) < 1e-6, `${tilt('e-')} + ${tilt('h+')}`);
});

test('graded gap: a face beside a graded region aligns by the standard levels at its edge', () => {
  // A uniform layer, then the graded one: step.value is the jump of e⁻'s standard level μ°+zFφ
  // across the face, read at the graded region's first node.
  const def = si(gap);
  def.regions.unshift({ material: 's', length: 1e-6, fixedCharge: -FARADAY * 0.0166 });
  def.regions[1].grade = gapFrom(1e-6);
  def.interfaces = [{ step: { species: 'e-', value: 2000 } }];
  const sol = new Device(def).solve();
  assert.ok(sol.converged, sol.warnings.join('; '));
  const g = sol.x.findIndex((x, k) => k > 0 && x === sol.x[k - 1]); // the face's second node
  assert.ok(Math.abs(sol.muStd['e-'][g] - sol.muStd['e-'][g - 1] - 2000) < 1e-6, `${sol.muStd['e-'][g] - sol.muStd['e-'][g - 1]}`);
});

test('grade: what it refuses, and says why', () => {
  const bad = (grade, msg) => assert.throws(() => new Device(neutral({ grade, muL: 0, muR: 0 })), (e) => e instanceof DeviceError && msg.test(e.message), JSON.stringify(grade));
  bad({ B: { mu0: { x: [0], values: [0] } } }, /grade\.B: unknown species/);
  bad({ A: { D: { x: [0], values: [1e-9] } } }, /grade\.A.*D/);
  bad({ A: { cRef: { x: [0, L], values: [1, 0] } } }, /grade\.A\.cRef\.values\[1\] must be a concentration > 0/);
  bad({ A: { mu0: 3 } }, /grade\.A\.mu0 must be a profile/);
  bad({ A: 'x' }, /grade\.A must be/);
});

// A pn diode whose p side is graded: its band gap widens toward the right contact (E_c rises,
// E_v falls). The other routes through the solver read the graded drift too: the impedance's
// linearisation (its 1 Hz Z is the steady state's dV/dI) and a transient settling where the
// steady solve goes.
const diode = (V) => {
  const def = si(gap);
  def.regions = [
    { material: 's', length: 1e-6, name: 'n', fixedCharge: FARADAY * 0.166 },
    { material: 's', length: 2e-6, name: 'p', fixedCharge: -FARADAY * 0.0166, grade: gapFrom(1e-6) },
  ];
  def.interfaces = [{ dipole: 0 }];
  def.contacts.right.V = V;
  def.bulkReactions = [{ equation: 'e- + h+ = 0', kf: { s: 1e8 } }];
  return def;
};

test('graded gap: the impedance at low frequency is the steady state’s dV/dI', () => {
  const V = 0.45, dV = 1e-4;
  const I = (v) => new Device(diode(v)).solve().terminals.right.current;
  const Rdc = (2 * dV) / (I(V + dV) - I(V - dV));
  const Z = new Device(diode(V)).impedance([1e-3]);
  assert.ok(Math.abs(Z.Z.re[0] / Rdc - 1) < 1e-5, `${Z.Z.re[0]} vs ${Rdc}`);
  // (and grading moved it: the same diode ungraded differs)
  const flat = diode(V);
  delete flat.regions[1].grade;
  assert.ok(Math.abs(new Device(flat).impedance([1e-3]).Z.re[0] / Rdc - 1) > 0.05);
});

test('graded gap: a transient settles where the steady solve goes', () => {
  // (from equilibrium, as the docs advise: a junction started from its doping's layout is a
  // violent start, graded or not)
  const dev = new Device(diode(0));
  dev.solve();
  dev.set({ contacts: { right: { V: 0.4 } } });
  const settled = dev.advance(1e-2), steady = new Device(diode(0.4)).solve();
  assert.ok(settled.done && settled.converged, settled.warnings.join('; '));
  assert.ok(Math.abs(settled.terminals.right.current / steady.terminals.right.current - 1) < 1e-5, `${settled.terminals.right.current} vs ${steady.terminals.right.current}`);
});

test('graded μ°: a closed box keeps the amount its c0 gives, and settles to c ∝ cRef e^{−μ°/RT}', () => {
  // A neutral species shut in (both ends blocked), started uniform at c0 = 5: the start must be
  // those amounts at every node (a transient's initial state), and the steady state holds the
  // same total in the Boltzmann profile of the graded μ° and cRef.
  const grade = { A: { mu0: { x: [0, L], values: [0, 3 * RT] }, cRef: { x: [0, L], values: [1000, 400] } } };
  const def = {
    ...neutral({ grade, muL: 0, muR: 0 }),
    contacts: { left: { species: {} }, right: { species: {} } },
  };
  def.regions[0].c0 = { A: 5 };
  const dev = new Device(def);
  const t0 = dev.advance(1e-12);
  for (let g = 0; g < t0.x.length; g++) assert.ok(Math.abs(t0.c.A[g] / 5 - 1) < 1e-6, `start at x = ${t0.x[g]}: ${t0.c.A[g]}`);
  const sol = new Device(def).solve();
  assert.ok(sol.converged, sol.warnings.join('; '));
  const shape = (x) => (1000 - 600 * (x / L)) * Math.exp(-3 * (x / L));
  // ∫ shape dx exactly, by parts: shape = (a − b x/L) e^{−k x/L}
  const k = 3, a = 1000, b = 600, I = L * ((a / k) * (1 - Math.exp(-k)) - b * (1 / (k * k) - Math.exp(-k) * (1 / k + 1 / (k * k))));
  for (const g of [0, sol.x.length >> 1, sol.x.length - 1]) {
    const want = (5 * L * shape(sol.x[g])) / I;
    assert.ok(Math.abs(sol.c.A[g] / want - 1) < 2e-3, `x = ${sol.x[g]}: ${sol.c.A[g]} vs ${want}`); // (the grid's quadrature of the total)
  }
});

test('graded gap: the diode is the limit of a staircase of uniform layers joined with no dipole', () => {
  // What a graded region means for charged species: μ° and cRef against the region's φ, with
  // D = −ε∇φ (no built-in displacement inside a region). A staircase of thin uniform layers,
  // each at its middle's values, joined by faces with no dipole (φ and D continuous), converges
  // on it.
  // (on a grid fine enough that its own error, ~5e-4 at hmax 20 nm, doesn't hide the limit)
  const V = 0.45, fine = { hmin: 1e-9, hmax: 5e-9 };
  const graded = new Device({ ...diode(V), grid: fine }).solve().terminals.right.current;
  const at = (p, x) => p.values[0] + ((p.values[1] - p.values[0]) * (x - p.x[0])) / (p.x[1] - p.x[0]);
  const errs = [4, 16, 64].map((k) => {
    const def = diode(V), base = def.materials.s, x0 = 1e-6, w = 2e-6 / k;
    def.regions = [def.regions[0]];
    def.interfaces = [];
    for (let j = 0; j < k; j++) {
      const xm = x0 + (j + 0.5) * w, sp = structuredClone(base.species);
      sp['e-'].mu0 = at(gap['e-'].mu0, xm - x0);
      sp['e-'].cRef = at(gap['e-'].cRef, xm - x0);
      sp['h+'].mu0 = at(gap['h+'].mu0, xm - x0);
      def.materials[`s${j}`] = { ...base, species: sp };
      def.regions.push({ material: `s${j}`, length: w, fixedCharge: -FARADAY * 0.0166 });
      def.interfaces.push({ dipole: 0 });
    }
    const sol = new Device({ ...def, grid: fine }).solve();
    assert.ok(sol.converged, sol.warnings.join('; '));
    return Math.abs(sol.terminals.right.current / graded - 1);
  });
  assert.ok(errs[2] < 1e-4, `64 layers: ${errs}`); // (1.3e-2, 3.8e-4, 3.4e-5)
  assert.ok(errs[1] < errs[0] / 8 && errs[2] < errs[1] / 8, `converging: ${errs}`);
});

test('graded gap: φ → φ + f with μ°_i → μ°_i − z_i F f is a symmetry where faces carry the dipoles', () => {
  // The gauge freedom of a graded medium (ESBD's inhomog: a compensating built-in displacement
  // D_inhomog absorbs any change of φ's convention). driftlet takes D_inhomog only at faces (a
  // dipole: φ steps, D doesn't), so on a staircase it's exact: the same diode with each face given
  // a dipole and each layer's μ° shifted to match carries the same current. Within one graded
  // region there is no D_inhomog, which fixes the convention there (docs/device.md).
  const V = 0.45, k = 16, fine = { hmin: 1e-9, hmax: 5e-9 };
  const at = (p, x) => p.values[0] + ((p.values[1] - p.values[0]) * (x - p.x[0])) / (p.x[1] - p.x[0]);
  const stair = (dipoles) => {
    const def = diode(V), base = def.materials.s, x0 = 1e-6, w = 2e-6 / k;
    def.regions = [def.regions[0]];
    def.interfaces = [];
    let f = 0; // φ's shift in this layer: the dipoles crossed so far, V
    for (let j = 0; j < k; j++) {
      f += dipoles[j];
      const xm = x0 + (j + 0.5) * w, sp = structuredClone(base.species);
      sp['e-'].mu0 = at(gap['e-'].mu0, xm - x0) + FARADAY * f; // μ° − zFf, z = −1
      sp['e-'].cRef = at(gap['e-'].cRef, xm - x0);
      sp['h+'].mu0 = at(gap['h+'].mu0, xm - x0) - FARADAY * f;
      def.materials[`s${j}`] = { ...base, species: sp };
      def.regions.push({ material: `s${j}`, length: w, fixedCharge: -FARADAY * 0.0166 });
      def.interfaces.push({ dipole: dipoles[j] });
    }
    const sol = new Device({ ...def, grid: fine }).solve();
    assert.ok(sol.converged, sol.warnings.join('; '));
    return sol;
  };
  const plain = stair(Array(k).fill(0)), shifted = stair(Array.from({ length: k }, (_, j) => 0.05 * Math.sin(j)));
  assert.ok(Math.abs(shifted.terminals.right.current / plain.terminals.right.current - 1) < 1e-9, `${shifted.terminals.right.current} vs ${plain.terminals.right.current}`);
  for (const g of [10, plain.x.length >> 1, plain.x.length - 10]) assert.ok(Math.abs(shifted.V['e-'][g] - plain.V['e-'][g]) < 1e-9 && Math.abs(shifted.Vstd['e-'][g] - plain.Vstd['e-'][g]) < 1e-9, `levels at node ${g}`);
});
