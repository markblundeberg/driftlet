import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { alloy, build, describe, layer, ohmic, semiconductor } from '../src/kit.js';

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

test('graded D: a neutral species through a D that varies is a series of resistances, N = Δc / ∫dx/D', () => {
  // Uniform μ°, so no drift: N = −D(x) dc/dx is the same everywhere, and Δc = N ∫dx/D. D rises
  // linearly 10-fold, ∫dx/D = L ln 10 / (D₁ − D₀); the grid's sum of h/D at the midpoints is
  // second order in h.
  const D0 = 1e-10, D1 = 1e-9, cL = cRef, cR = cRef * Math.exp(-1);
  const want = (cL - cR) / ((L * Math.log(D1 / D0)) / (D1 - D0));
  const grade = { A: { D: { x: [0, L], values: [D0, D1] } } };
  const errs = [L / 20, L / 80, L / 320].map((h) => {
    const sol = new Device(neutral({ grade, muL: 0, muR: -RT, grid: { hmin: h, hmax: h } })).solve();
    assert.ok(sol.converged, sol.warnings.join('; '));
    const N = sol.contacts.left.flux.A;
    assert.ok(Math.abs(sol.flux.A[sol.x.length >> 1] / N - 1) < 1e-9); // (the same through every segment)
    return Math.abs(N / want - 1);
  });
  assert.ok(errs[2] < 3e-5, `${errs}`); // (3.5e-3, 2.3e-4, 1.4e-5)
  assert.ok(errs[1] < errs[0] / 12 && errs[2] < errs[1] / 12, `second order in h: ${errs}`);
  // and sol.D reads the profile at each node
  const sol = new Device(neutral({ grade, muL: 0, muR: -RT })).solve(), g = sol.x.length >> 1;
  assert.ok(Math.abs(sol.D.A[g] / (D0 + ((D1 - D0) * sol.x[g]) / L) - 1) < 1e-12);
});

test('graded ε: a dielectric between two metals holds D = ε₀ ΔV / ∫dx/ε_r', () => {
  // An insulator with no carriers, its ε_r rising linearly 3.9 → 25 across it, between two metal
  // plates the displacement is the same through it, and the plates' charge is the series
  // capacitance's, ε₀ / ∫dx/ε_r = ε₀ (e₁ − e₀) / (d ln(e₁/e₀)), times the voltage.
  const d = 20e-9, e0 = 3.9, e1 = 25, dV = 0.3, EPS0 = 8.8541878128e-12;
  // (the plates meet it through a Helmholtz layer each, C_H in series)
  const CH = 100, face = { phi: { type: 'capacitive', C: CH }, zeroCharge: 0 };
  const want = dV / ((d * Math.log(e1 / e0)) / (EPS0 * (e1 - e0)) + 2 / CH);
  const def = (h, grade) => ({
    T,
    species: [{ name: 'e-', z: -1 }],
    materials: { Al: { conductor: { species: 'e-', conductivity: 3.5e7 } }, ox: { epsr: e0, species: {} } },
    regions: [
      { material: 'Al', length: 50e-9 },
      { material: 'ox', length: d, grade },
      { material: 'Al', length: 50e-9 },
    ],
    interfaces: [face, face],
    contacts: { left: { V: dV, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' }, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' } },
    grid: { hmin: h, hmax: h },
  });
  const grade = { epsr: { x: [50e-9, 50e-9 + d], values: [e0, e1] } };
  const errs = [d / 10, d / 40, d / 160].map((h) => {
    const sol = new Device(def(h, grade)).solve();
    assert.ok(sol.converged, sol.warnings.join('; '));
    const D = sol.interfaces[0].D;
    assert.ok(Math.abs(sol.interfaces[1].D / D - 1) < 1e-9, `${sol.interfaces[1].D} vs ${D}`);
    return Math.abs(Math.abs(D) / want - 1);
  });
  assert.ok(errs[2] < 5e-5, `${errs}`); // (6.1e-3, 4.0e-4, 2.5e-5)
  assert.ok(errs[1] < errs[0] / 12 && errs[2] < errs[1] / 12, `second order in h: ${errs}`);
  // and its impedance at 1 kHz is that capacitance's, 1/(iωC)
  const dev = new Device(def(d / 160, grade)), C = Math.abs(dev.solve().interfaces[0].D) / dV, Z = dev.impedance([1e3]).Z;
  assert.ok(Math.abs(-1 / (2 * Math.PI * 1e3 * Z.im[0]) / C - 1) < 1e-6 && Math.abs(Z.re[0]) < 1e-6 * Math.abs(Z.im[0]), `Z (${Z.re[0]}, ${Z.im[0]}) for C ${C}`);
  // (ungraded, the plain parallel plate's ε₀ε_r ΔV/d)
  const flat = new Device(def(d / 10)).solve();
  assert.ok(Math.abs(Math.abs(flat.interfaces[0].D) / (dV / (d / (EPS0 * e0) + 2 / CH)) - 1) < 1e-9, `${flat.interfaces[0].D}`);
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
  bad({ A: { D: { x: [0, L], values: [1e-9, 0] } } }, /grade\.A\.D\.values\[1\] must be a diffusion coefficient > 0/);
  bad({ epsr: { x: [0, L], values: [2, 3] } }, /grade\.epsr: 'm' is strictly neutral/);
  bad({ A: { D: { x: [0], values: [1e-9] }, eps: 1 } }, /grade\.A.*eps/);
  bad({ A: { cRef: { x: [0, L], values: [1, 0] } } }, /grade\.A\.cRef\.values\[1\] must be a concentration > 0/);
  bad({ A: { mu0: 3 } }, /grade\.A\.mu0 must be a profile/);
  bad({ A: 'x' }, /grade\.A must be/);
});

// A pn diode whose p side is graded: its band gap widens toward the right contact (E_c rises,
// E_v falls). The other routes through the solver read the graded drift too: the impedance's
// linearisation (its 1 Hz Z is the steady state's dV/dI) and a transient settling where the
// steady solve goes.
const diode = (V, grade = gapFrom(1e-6)) => {
  const def = si(gap);
  def.regions = [
    { material: 's', length: 1e-6, name: 'n', fixedCharge: FARADAY * 0.166 },
    { material: 's', length: 2e-6, name: 'p', fixedCharge: -FARADAY * 0.0166, grade },
  ];
  def.interfaces = [{ dipole: 0 }];
  def.contacts.right.V = V;
  def.bulkReactions = [{ equation: 'e- + h+ = 0', kf: { s: 1e8 } }];
  return def;
};

// The same p layer with its permittivity and both carriers' diffusion coefficients graded too, as
// an alloy's are (ε up 11.7 → 13, D_e⁻ down by half, D_h⁺ up by a third).
const allFrom = (x0) => {
  const g = gapFrom(x0), at = [x0, x0 + 2e-6];
  g['e-'].D = { x: at, values: [3.6e-3, 1.8e-3] };
  g['h+'].D = { x: at, values: [1.2e-3, 1.6e-3] };
  g.epsr = { x: at, values: [11.7, 13] };
  return g;
};

test('graded gap: the impedance at low frequency is the steady state’s dV/dI', () => {
  const V = 0.45, dV = 1e-4;
  const Rs = [gapFrom(1e-6), allFrom(1e-6)].map((grade) => {
    const I = (v) => new Device(diode(v, grade)).solve().terminals.right.current;
    const Rdc = (2 * dV) / (I(V + dV) - I(V - dV));
    const Z = new Device(diode(V, grade)).impedance([1e-3]);
    assert.ok(Math.abs(Z.Z.re[0] / Rdc - 1) < 1e-5, `${Z.Z.re[0]} vs ${Rdc}`);
    return Rdc;
  });
  // (and grading moved it: the same diode ungraded differs, and grading ε and D too moves it again)
  const flat = diode(V);
  delete flat.regions[1].grade;
  assert.ok(Math.abs(new Device(flat).impedance([1e-3]).Z.re[0] / Rs[0] - 1) > 0.05);
  assert.ok(Math.abs(Rs[1] / Rs[0] - 1) > 0.05, `${Rs}`);
});

test('graded gap: a transient settles where the steady solve goes', () => {
  // (from equilibrium, as the docs advise: a junction started from its doping's layout is a
  // violent start, graded or not)
  for (const grade of [gapFrom(1e-6), allFrom(1e-6)]) {
    const dev = new Device(diode(0, grade));
    dev.solve();
    dev.set({ contacts: { right: { V: 0.4 } } });
    const settled = dev.advance(1e-2), steady = new Device(diode(0.4, grade)).solve();
    assert.ok(settled.done && settled.converged, settled.warnings.join('; '));
    assert.ok(Math.abs(settled.terminals.right.current / steady.terminals.right.current - 1) < 1e-5, `${settled.terminals.right.current} vs ${steady.terminals.right.current}`);
  }
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

// A profile's value at x, piecewise linear (flat beyond its ends).
const atX = ({ x: xs, values: v }, x) => {
  if (x <= xs[0]) return v[0];
  for (let q = 1; q < xs.length; q++) if (x <= xs[q]) return v[q - 1] + ((v[q] - v[q - 1]) * (x - xs[q - 1])) / (xs[q] - xs[q - 1]);
  return v[v.length - 1];
};
// The diode's p layer graded between flat ends (a quarter micron each): a staircase's end layers
// then hold the values its contact and its face read (with gradings that run to the ends, they
// hold the values half a layer in, a first-order error in 1/k that isn't the graded region's).
const rampFrom = (x0, what) => {
  const at = [x0, x0 + 0.25e-6, x0 + 1.75e-6, x0 + 2e-6], two = (a, b) => ({ x: at, values: [a, a, b, b] });
  const g = { 'e-': {}, 'h+': {} };
  if (what.gap) Object.assign(g['e-'], { mu0: two(0, -15e3), cRef: two(0.047, 0.03) }), (g['h+'].mu0 = two(108e3, 120e3));
  if (what.D) (g['e-'].D = two(3.6e-3, 1.8e-3)), (g['h+'].D = two(1.2e-3, 1.6e-3));
  return g;
};
// The diode with its p layer as k uniform layers, each at its middle's values, joined by faces
// with no dipole (φ and D continuous), recombining as the graded layer does; its current.
const staircase = (V, grade, k, fine) => {
  const def = diode(V), base = def.materials.s, x0 = 1e-6, w = 2e-6 / k;
  def.regions = [def.regions[0]];
  def.interfaces = [];
  for (let j = 0; j < k; j++) {
    const xm = x0 + (j + 0.5) * w, sp = structuredClone(base.species);
    for (const [nm, gr] of Object.entries(grade)) for (const [f, p] of Object.entries(gr)) sp[nm][f] = atX(p, xm);
    def.materials[`s${j}`] = { ...base, species: sp };
    def.regions.push({ material: `s${j}`, length: w, fixedCharge: -FARADAY * 0.0166 });
    def.interfaces.push({ dipole: 0 });
  }
  def.bulkReactions = [{ equation: 'e- + h+ = 0', kf: Object.fromEntries(Object.keys(def.materials).map((m) => [m, 1e8])) }];
  const sol = new Device({ ...def, grid: fine }).solve();
  assert.ok(sol.converged, sol.warnings.join('; '));
  return sol.terminals.right.current;
};

test('graded gap: the diode is the limit of a staircase of uniform layers joined with no dipole', () => {
  // What a graded region means for charged species: μ° and cRef against the region's φ, with
  // D = −ε∇φ (no built-in displacement inside a region).
  // (on a grid fine enough that its own error, ~5e-4 at hmax 20 nm, doesn't hide the limit)
  const V = 0.45, fine = { hmin: 1e-9, hmax: 5e-9 }, grade = rampFrom(1e-6, { gap: true });
  const graded = new Device({ ...diode(V, grade), grid: fine }).solve().terminals.right.current;
  const errs = [4, 16, 64].map((k) => Math.abs(staircase(V, grade, k, fine) / graded - 1));
  assert.ok(errs[2] < 3e-4, `64 layers: ${errs}`); // (1.8e-2, 1.5e-3, 1.3e-4)
  assert.ok(errs[1] < errs[0] / 8 && errs[2] < errs[1] / 8, `second order in 1/k: ${errs}`);
});

test('graded D: the diode is again the limit of a staircase, each layer at its middle’s D', () => {
  const V = 0.45, fine = { hmin: 1e-9, hmax: 5e-9 }, grade = rampFrom(1e-6, { D: true });
  const graded = new Device({ ...diode(V, grade), grid: fine }).solve().terminals.right.current;
  const errs = [4, 16, 64].map((k) => Math.abs(staircase(V, grade, k, fine) / graded - 1));
  assert.ok(errs[2] < 1e-4, `64 layers: ${errs}`); // (7.5e-3, 3.3e-4, 8.5e-6)
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

// Graded doping: fixedCharge as a profile. The diode's p layer with its acceptors falling 100×
// across it, exponentially between flat ends (a drift transistor's base: Kroemer, 1954), against
// a staircase of uniform layers each at its middle's doping.
const dopedFrom = (x0) => {
  const xs = [x0, ...Array.from({ length: 31 }, (_, j) => x0 + 0.25e-6 + (j * 1.5e-6) / 30), x0 + 2e-6];
  return { x: xs, values: xs.map((x) => -FARADAY * 0.166 * Math.exp(Math.log(0.01) * Math.min(1, Math.max(0, (x - x0 - 0.25e-6) / 1.5e-6)))) };
};
const dopedDiode = (V, fixed, grid) => {
  const def = diode(V, {});
  def.regions[1].fixedCharge = fixed;
  return grid ? { ...def, grid } : def;
};
test('graded doping: fixedCharge as a profile is the limit of a staircase of uniformly doped layers', () => {
  const V = 0.45, fine = { hmin: 1e-9, hmax: 5e-9 }, prof = dopedFrom(1e-6);
  const graded = new Device(dopedDiode(V, prof, fine)).solve();
  assert.ok(graded.converged, graded.warnings.join('; '));
  const stair = (k) => {
    const def = dopedDiode(V, 0, fine), w = 2e-6 / k;
    def.regions = [def.regions[0], ...Array.from({ length: k }, (_, j) => ({ material: 's', length: w, fixedCharge: atX(prof, 1e-6 + (j + 0.5) * w) }))];
    def.interfaces = def.regions.slice(1).map(() => ({ dipole: 0 }));
    const sol = new Device(def).solve();
    assert.ok(sol.converged, sol.warnings.join('; '));
    return sol.terminals.right.current;
  };
  const errs = [4, 16, 64].map((k) => Math.abs(stair(k) / graded.terminals.right.current - 1));
  assert.ok(errs[2] < 1e-4, `64 layers: ${errs}`); // (9.0e-2, 3.2e-3, 1.9e-5)
  assert.ok(errs[1] < errs[0] / 8 && errs[2] < errs[1] / 8, `second order in 1/k: ${errs}`);
});

test('graded doping: another gauge’s built-in polarization P0(x) is its bound charge, −dP0/dx, as fixed charge', () => {
  // φ → φ + f with each μ°_i → μ°_i − z_i F f is the same device where a built-in displacement
  // D_inhomog = ε f′ (a static polarization P0) takes up the change (ESBD's inhomog). driftlet
  // has no D_inhomog inside a region, but Poisson's equation sees only ρ − dP0/dx, so the same
  // device in that gauge is its μ° shifted and −(ε f′)′ added to its fixed charge (and a face's
  // sheetCharge where P0 jumps; none here). Every level and current is
  // the plain diode's, φ shifted by f, to the grid's second order (the box's second difference of
  // f against f″): exact in the limit, so a reader of the profile anywhere else fails it.
  const V = 0.45, x0 = 1e-6, W = 2e-6, A = 0.1, eps = 11.7 * EPS0, q0 = -FARADAY * 0.0166;
  // f = A sin⁴(π(x − x0)/W): f, f′ and f″ all 0 at both ends, so the n side, the face and the
  // contact see nothing (a 'bulk' contact holds its end neutral with the fixed charge counted, so
  // a bound charge reaching it moves the contact: true neutrality for a real polarization, and
  // not a gauge's).
  const xs = Array.from({ length: 4001 }, (_, j) => x0 + (j * W) / 4000), u = (x) => (Math.PI * (x - x0)) / W, k2 = (Math.PI / W) ** 2;
  const f = (x) => A * Math.sin(u(x)) ** 4, f2 = (x) => 4 * A * k2 * Math.sin(u(x)) ** 2 * (3 * Math.cos(u(x)) ** 2 - Math.sin(u(x)) ** 2);
  const prof = (fn) => ({ x: xs, values: xs.map(fn) }), muE = prof((x) => FARADAY * f(x));
  const run = (h) => {
    const grid = { hmin: h / 10, hmax: h };
    const plain = new Device({ ...diode(V, {}), grid }).solve();
    const def = { ...diode(V, { 'e-': { mu0: muE }, 'h+': { mu0: prof((x) => 108e3 - FARADAY * f(x)) } }), grid };
    def.regions[1].fixedCharge = prof((x) => q0 - eps * f2(x));
    const shifted = new Device(def).solve();
    assert.ok(plain.converged && shifted.converged, shifted.warnings.join('; '));
    let dphi = 0, dV = 0;
    for (let g = 0; g < plain.x.length; g++) {
      const shift = plain.x[g] > x0 ? atX(muE, plain.x[g]) / FARADAY : 0; // (f as the profile gives it)
      dphi = Math.max(dphi, Math.abs(shifted.phi[g] - plain.phi[g] - shift));
      for (const sp of ['e-', 'h+']) dV = Math.max(dV, Math.abs(shifted.V[sp][g] - plain.V[sp][g]));
    }
    return [dphi, dV, Math.abs(shifted.terminals.right.current / plain.terminals.right.current - 1)];
  };
  const coarse = run(10e-9), fine = run(5e-9);
  assert.ok(fine[0] < 4e-7 && fine[1] < 6e-8 && fine[2] < 3e-7, `φ − f, levels (V), current: ${fine}`); // (2.1e-7 V, 3.3e-8 V, 1.2e-7)
  // (falling 4× a halving from here: 5.1e-8, 7.6e-9, 3.1e-8 at 2.5 nm; 1.4e-8, 1.9e-9, 8.0e-9 at 1.25 nm)
  for (let k = 0; k < 3; k++) assert.ok(fine[k] < coarse[k] / 2.5, `second order: ${coarse} → ${fine}`);
});

test('graded doping: what build() writes from donor and acceptor profiles, and what is refused, warned of and described', () => {
  // Their difference on both profiles' points (exact: each is linear between its points, and
  // constant beyond its ends).
  const lib = semiconductor('Si', { T: 300 }), si1 = (opts, grid) => build({ T: 300, library: [lib], stack: [ohmic(0), layer('Si', 1e-6, opts), ohmic(0)], ...(grid ? { grid } : {}) });
  const def = si1({ donors: { x: [0, 1e-6], values: [2, 1] }, acceptors: { x: [0.5e-6, 1e-6], values: [0.5, 1.5] } });
  assert.deepEqual(def.regions[0].fixedCharge.x, [0, 0.5e-6, 1e-6]);
  def.regions[0].fixedCharge.values.forEach((q, k) => assert.ok(Math.abs(q / FARADAY - [1.5, 1, -0.5][k]) < 1e-12, `${k}: ${q / FARADAY}`));
  assert.match(describe(def), /fixed charge graded, donor-like 9\.0\de\+17 cm⁻³ to acceptor-like 3\.0\de\+17 cm⁻³/);
  assert.throws(() => si1({ donors: { x: [1e-6, 0], values: [1, 2] } }), /donors\.x must be finite and strictly increasing/);
  assert.throws(() => si1({ acceptors: { x: [0, 1e-6], values: [1, -1] } }), /acceptors must be a concentration ≥ 0.*at each point of its profile/);
  // A step the cells can't follow is warned of, as a starting profile's is.
  assert.match(describe(si1({ donors: { x: [0.5e-6, 0.5e-6 + 1e-10], values: [1, 100] } }, { hmax: 1e-7 })), /regions\[0\]\.fixedCharge: the grid doesn't resolve the profile/);
  // A metal holds no bulk charge; an insertion host's fixed charge sets its background, which is uniform.
  const metal = { species: [{ name: 'e-', z: -1 }], materials: { Al: { conductor: { species: 'e-', conductivity: 3.5e7 } } }, regions: [{ material: 'Al', length: 1e-6, fixedCharge: { x: [0, 1e-6], values: [0, 1] } }], contacts: { left: { V: 0 }, right: { V: 0 } } };
  assert.throws(() => new Device(metal), /fixedCharge: a conductor region is neutral in bulk/);
  const host = {
    species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: { host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 2e4 }, 'e-': { D: 1e-4, mu0: 0, cRef: 2e4 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 2e4, A: [0] }] } },
    regions: [{ material: 'host', length: 1e-6, fixedCharge: { x: [0, 1e-6], values: [-1, -2] }, c0: { 'Li+': 1e4, 'e-': 1e4 } }],
    contacts: { left: { phi: 'neutral' }, right: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' } },
  };
  assert.throws(() => new Device(host), /insertion host the fixed charge sets its background carriers, which are uniform/);
});

test('graded μ° and cRef under Fermi–Dirac statistics: the diode is the limit of a staircase of uniform layers', () => {
  // Every model reads ζ = η − μ°/RT − zφ̂ at the node's own μ°, so a graded μ° needs nothing of
  // the statistics; a Fermi–Dirac c = c_ref 𝓕(ζ) scales exactly with the node's own c_ref (other
  // models' cRef is refused). Both carriers Fermi–Dirac: the n side degenerate (n ≈ 10 N_c), the
  // p side's holes at about half their N_v. (Without the scaling, the staircase settles 21% off.)
  const V = 0.45, fine = { hmin: 1e-9, hmax: 5e-9 }, x0 = 1e-6, at = [x0, x0 + 0.25e-6, x0 + 1.75e-6, x0 + 2e-6], two = (a, b) => ({ x: at, values: [a, a, b, b] });
  const grade = { 'e-': { mu0: two(0, -15e3), cRef: two(0.047, 0.03) }, 'h+': { mu0: two(108e3, 120e3), cRef: two(0.031, 0.06) } };
  const fd = () => {
    const def = diode(V, grade);
    def.materials.s.statistics = [{ type: 'fermi-dirac', species: ['e-'] }, { type: 'fermi-dirac', species: ['h+'] }];
    def.regions[0].fixedCharge = FARADAY * 0.5;
    return { ...def, grid: fine };
  };
  const graded = new Device(fd()).solve();
  assert.ok(graded.converged, graded.warnings.join('; '));
  const stair = (k) => {
    const def = fd(), base = def.materials.s, w = 2e-6 / k;
    def.regions = [def.regions[0]];
    def.interfaces = [];
    for (let j = 0; j < k; j++) {
      const xm = x0 + (j + 0.5) * w, m = structuredClone(base);
      m.species['e-'].mu0 = atX(grade['e-'].mu0, xm);
      m.species['h+'].mu0 = atX(grade['h+'].mu0, xm);
      m.species['e-'].cRef = atX(grade['e-'].cRef, xm);
      m.species['h+'].cRef = atX(grade['h+'].cRef, xm);
      def.materials[`s${j}`] = m;
      def.regions.push({ material: `s${j}`, length: w, fixedCharge: -FARADAY * 0.0166 });
      def.interfaces.push({ dipole: 0 });
    }
    def.bulkReactions = [{ equation: 'e- + h+ = 0', kf: Object.fromEntries(Object.keys(def.materials).map((m) => [m, 1e8])) }];
    const sol = new Device(def).solve();
    assert.ok(sol.converged, sol.warnings.join('; '));
    return sol.terminals.right.current;
  };
  const errs = [4, 16, 64].map((k) => Math.abs(stair(k) / graded.terminals.right.current - 1));
  assert.ok(errs[2] < 1e-3 && errs[1] < errs[0] / 8 && errs[2] < errs[1] / 8, `second order in 1/k: ${errs}`); // (5.6e-2, 4.5e-3, 3.4e-4)
  const def = fd();
  def.materials.s.statistics = [{ type: 'lattice', species: ['e-'], cMax: 1 }];
  assert.throws(() => new Device(def), /grade\.e-\.cRef: 'e-' has lattice statistics in 's'; only an ideal or Fermi–Dirac species' cRef can be graded so far/);
});

// driftlet/kit's alloy(): a layer between two library materials, graded by composition.
const algaas = (x, name = 'AlGaAs') => {
  // docs/data.md's estimates: E_g = 1.424 + 1.247x eV, m_e = 0.067 + 0.083x
  const lib = semiconductor('GaAs', { material: name, T: 300 }), m = lib.materials[name], VT = (GAS_CONSTANT * 300) / FARADAY;
  m.species['h+'] = { ...m.species['h+'], mu0: units.eV(1.424 + 1.247 * x) };
  m.species['e-'] = { ...m.species['e-'], D: 1000e-4 * VT, cRef: m.species['e-'].cRef * Math.pow((0.067 + 0.083 * x) / 0.067, 1.5) };
  m.epsr = 13.1 - 3 * x;
  return lib;
};
const dEc = 0.62 * 1.247 * 0.3; // eV, AlGaAs's conduction band above GaAs's

test('alloy: held at one end, it is that material, lined up as a face would line it up', () => {
  // n-AlGaAs | p-GaAs, joined by a face with a conduction-band step, against n-AlGaAs | an alloy
  // from AlGaAs to GaAs held at fraction 1 throughout, given the same step: the alloy folds the
  // face's dipole into GaAs's μ° (φ → φ + f with μ° → μ° − zFf), so it's the same device.
  const lib = [semiconductor('GaAs', { T: 300 }), algaas(0.3)], kf = { GaAs: 1e8, AlGaAs: 1e8 };
  const step = { species: 'e-', value: units.eV(-dEc) };
  const def = (p) => build({ T: 300, library: lib, stack: [ohmic(0), layer('AlGaAs', 1e-6, { donors: units.perCm3(1e17) }), ...p, ohmic(1.0)], bulkReactions: [{ equation: 'e- + h+ = 0', kf }], grid: { hmin: 1e-9, hmax: 20e-9 } });
  const faced = new Device(def([{ step }, layer('GaAs', 1e-6, { acceptors: units.perCm3(1e17) })])).solve();
  const alloyed = new Device(def([alloy('AlGaAs', 'GaAs', 1e-6, { fraction: [1, 1], step, acceptors: units.perCm3(1e17) })])).solve();
  assert.ok(faced.converged && alloyed.converged);
  assert.ok(Math.abs(alloyed.current / faced.current - 1) < 1e-7, `${alloyed.current} vs ${faced.current}`);
  // the same levels, and the same band edges (φ differs by the folded dipole on the right)
  for (const g of [5, faced.x.length >> 1, faced.x.length - 5]) {
    for (const nm of ['e-', 'h+']) assert.ok(Math.abs(alloyed.V[nm][g] - faced.V[nm][g]) < 1e-7 && Math.abs(alloyed.Vstd[nm][g] - faced.Vstd[nm][g]) < 1e-7, `${nm} at node ${g}`);
  }
  assert.ok(Math.abs(alloyed.phi[faced.x.length - 5] - faced.phi[faced.x.length - 5]) > 0.1); // (and φ's convention did move)
});

test('alloy: what build() writes, against the device’s x', () => {
  const lib = [semiconductor('GaAs', { T: 300 }), algaas(0.3)];
  const step = { species: 'e-', value: units.eV(dEc) };
  const def = build({ T: 300, library: lib, stack: [ohmic(0), layer('GaAs', 0.5e-6), alloy('GaAs', 'AlGaAs', 1e-6, { step, fraction: { x: [0, 0.4e-6, 1e-6], values: [0, 0.5, 1] } }), ohmic(0)] });
  const g = def.regions[1].grade, A = lib[0].materials.GaAs, B = lib[1].materials.AlGaAs;
  assert.deepEqual(g['e-'].mu0.x, [0.5e-6, 0.9e-6, 1.5e-6]); // (the layer's knots, moved to where it starts)
  // e⁻'s standard level rises by the step across it, half of it at the composition's midpoint
  assert.ok(Math.abs(g['e-'].mu0.values[2] - g['e-'].mu0.values[0] - units.eV(dEc)) < 1e-9);
  assert.ok(Math.abs(g['e-'].mu0.values[1] - g['e-'].mu0.values[0] - units.eV(dEc) / 2) < 1e-9);
  // h⁺'s by the rest of the gap's widening: E_v = −(μ°_h⁺ + Fφ) falls by ΔE_g − ΔE_c
  const dEg = B.species['h+'].mu0 - A.species['h+'].mu0 + B.species['e-'].mu0 - A.species['e-'].mu0;
  assert.ok(Math.abs(g['h+'].mu0.values[2] - g['h+'].mu0.values[0] - (dEg - units.eV(dEc))) < 1e-6);
  [A.epsr, (A.epsr + B.epsr) / 2, B.epsr].forEach((e, k) => assert.ok(Math.abs(g.epsr.values[k] - e) < 1e-12, `ε_r at knot ${k}`));
  assert.equal(g['h+'].D, undefined); // (what doesn't differ isn't graded)
  assert.equal(def.regions[1].material, 'GaAs');
  assert.ok(new Device(def).solve().converged);
});

test('alloy: a bowed value, less bow·f(1 − f), the ends unmoved', () => {
  // E_g(x) = (1 − x)E_A + x E_B − C x(1 − x): a band gap's bowing, split between its band edges
  // as the alignment is (here all of it in the valence band's, h⁺'s μ°), and the permittivity's.
  const lib = [semiconductor('GaAs', { T: 300 }), algaas(0.3)], step = { species: 'e-', value: units.eV(dEc) }, C = units.eV(0.4);
  const at = (bow) => build({ T: 300, library: lib, stack: [ohmic(0), layer('GaAs', 0.5e-6), alloy('GaAs', 'AlGaAs', 1e-6, { step, bow }), ohmic(0)] }).regions[1].grade;
  const plain = at(undefined), bowed = at({ 'h+': { mu0: C }, epsr: 0.5 }), mu = bowed['h+'].mu0, n = mu.x.length;
  assert.equal(n, 33);
  const lin = (p, x) => p.values[0] + ((p.values.at(-1) - p.values[0]) * (x - p.x[0])) / (p.x.at(-1) - p.x[0]);
  assert.ok(Math.abs(mu.values[0] - plain['h+'].mu0.values[0]) < 1e-9 && Math.abs(mu.values[n - 1] - plain['h+'].mu0.values.at(-1)) < 1e-9);
  assert.ok(Math.abs(mu.values[16] - (lin(plain['h+'].mu0, mu.x[16]) - C / 4)) < 1e-9, `${mu.values[16]}`); // (f = 1/2: C/4 below the line)
  assert.ok(Math.abs(bowed.epsr.values[8] - (lin(plain.epsr, bowed.epsr.x[8]) - 0.5 * 0.25 * 0.75)) < 1e-12);
  assert.ok(Math.abs(bowed['e-'].mu0.values[16] - lin(plain['e-'].mu0, bowed['e-'].mu0.x[16])) < 1e-9); // (e⁻'s unbowed)
  const dev = new Device(build({ T: 300, library: lib, stack: [ohmic(0), layer('GaAs', 0.5e-6), alloy('GaAs', 'AlGaAs', 1e-6, { step, bow: { 'h+': { mu0: C } } }), ohmic(0)] })), s = dev.solve();
  assert.ok(s.converged && !s.warnings.some((w) => w.includes("doesn't resolve the grading")), s.warnings.join('; '));
  const refuse = (bow, re) => assert.throws(() => at(bow), re);
  refuse({ 'h+': { cRef: 2 } }, /bow\.h\+ must be \{ mu0 \}/);
  refuse({ X: { mu0: 1 } }, /bow\.X: no species 'X'/);
  refuse(3, /bow must be/);
});

test('alloy: what it refuses, and says why', () => {
  const lib = [semiconductor('GaAs', { T: 300 }), algaas(0.3)];
  const bad = (item, msg, library = lib) => assert.throws(() => build({ T: 300, library, stack: [ohmic(0), item, ohmic(0)] }), (e) => e instanceof DeviceError && msg.test(e.message), msg.source);
  const step = { species: 'e-', value: 0 };
  bad(alloy('GaAs', 'AlGaAs', 1e-6), /there's no default: give step/);
  bad(alloy('GaAs', 'AlGaAs', 1e-6, { step, dipole: 0 }), /step or dipole, not both/);
  bad(alloy('GaAs', 'AlGaAs', 1e-6, { step, fraction: [0, 1.2] }), /between 0 \(GaAs\) and 1 \(AlGaAs\)/);
  bad(alloy('GaAs', 'InAs', 1e-6, { step }), /no material 'InAs'/);
  bad(alloy('GaAs', 'AlGaAs', 1e-6, { step: { species: 'x', value: 0 } }), /step must be \{ species, value \}/);
  const hot = algaas(0.3);
  hot.materials.AlGaAs.T = 350;
  bad(alloy('GaAs', 'AlGaAs', 1e-6, { step }), /differ in T/, [lib[0], hot]);
  assert.throws(() => alloy('GaAs', 'AlGaAs', 1e-6, { step, grade: {} }), /writes the layer's grade itself/);
});

test('graded cRef: the degeneracy warning reads the graded band’s density of states', () => {
  // n-type at 0.08 mol/m³, above the material's cRef for e⁻ (0.047) but below the graded one (0.5)
  const def = (grade) => si(grade, { fixedCharge: FARADAY * 0.08 });
  const degenerate = (sol) => sol.warnings.some((w) => /e- reaches .* above its cRef/.test(w));
  assert.ok(degenerate(new Device(def(undefined)).solve()));
  assert.ok(!degenerate(new Device(def({ 'e-': { cRef: { x: [0, 2e-6], values: [0.5, 0.5] } } })).solve()));
  assert.ok(degenerate(new Device(def({ 'e-': { cRef: { x: [0, 2e-6], values: [0.5, 0.06] } } })).solve())); // (0.06 at the right end)
});

test('graded D: the time scales that size the steps see it, slowest and fastest', () => {
  // (slowestTime sizes the giant steps a steady solve by time steps takes, 1e6·τ and on: a D
  // graded far below the material's left them short by that ratio; Fable's smell review)
  const dev = (values) => new Device(neutral({ grade: { A: { D: { x: [0, L], values } } }, muL: 0, muR: 0 }));
  assert.equal(dev([D, 1e-19]).solver.slowestTime(), (L * L) / 1e-19);
  const fast = dev([D, 10]).solver;
  assert.ok(fast.fastestTime() < new Device(neutral({ mu0: 0, muL: 0, muR: 0 })).solver.fastestTime() / 1e9);
  // The material's own D, where a grading replaces it everywhere, isn't counted (a material's
  // 1e-22 made the slowest time 1e12 s for 0.1 s, and describe() said so; Fable's re-check).
  const def = neutral({ grade: { A: { D: { x: [0, L], values: [D, 2 * D] } } }, muL: 0, muR: 0 });
  def.materials.m.species.A.D = 1e-22;
  assert.equal(new Device(def).solver.slowestTime(), (L * L) / D);
  assert.match(describe(def), /diffusion across region 0: L²\/D up to 100 ms/);
});

test('graded μ°: a port’s bath of a species graded where the port is, is refused', () => {
  // (it read the material's ungraded μ°, holding c = 0.1 and 0.01 at the window's ends for c 1;
  // Fable's review)
  const def = neutral({ grade: { A: { mu0: { x: [0, L], values: [RT * Math.log(10), RT * Math.log(100)] } } }, muL: 0, muR: 0 });
  def.ports = [{ region: 0, V: 0, bath: { c: { A: 1 } }, species: { A: 'equilibrium' } }];
  assert.throws(() => new Device(def), /grades 'A's mu0 or cRef, so a composition outside gives no one level/);
});

// A grading the grid doesn't follow, or one that misses its region: from Fable's review of graded
// regions (12 Nov 2026). Each was converged and unwarned.
test('graded: a grading the cells miss, or a profile that misses its region, is warned of; a resolved one is not', () => {
  const unresolved = (sol) => sol.warnings.some((w) => /grade: the grid doesn't resolve the grading/.test(w));
  // A μ° spike 20 RT high and 2 nm wide in the middle of 10 µm: no node on it, the flux was
  // 4853× the exact one (D·cRef(1 − e⁻¹)/∫e^{μ°/RT}dx). With a region boundary at each side,
  // which the grid grades toward, it's resolved.
  const spike = { x: [0, 4.999e-6, 5e-6, 5.001e-6, L], values: [0, 0, 20 * RT, 0, 0] };
  const s = new Device(neutral({ muL: 0, muR: -RT, grade: { A: { mu0: spike } } })).solve();
  assert.ok(unresolved(s), s.warnings.join('; '));
  const sharp = neutral({ muL: 0, muR: -RT });
  sharp.regions = [[0, 4.999e-6], [4.999e-6, 5.001e-6], [5.001e-6, L]].map(([a, b]) => ({ material: 'm', length: b - a, grade: { A: { mu0: spike } } }));
  sharp.interfaces = [{}, {}];
  sharp.grid = { hmin: 1e-11, hmax: 0.5e-6, ratio: 1.2 };
  assert.ok(!unresolved(new Device(sharp).solve()));
  // cRef linear across six decades: the flux takes ln cRef as linear in each cell (0.094 of
  // exact on the default grid).
  const fall = new Device(neutral({ muL: 0, muR: -RT, grade: { A: { cRef: { x: [0, L], values: [1000, 1e-3] } } } })).solve();
  assert.ok(unresolved(fall));
  // A gentle grading on the default grid isn't flagged.
  assert.ok(!unresolved(new Device(neutral({ mu0: 5 * RT, muL: 0, muR: -RT })).solve()));
  // A profile written from the layer's own left end, for a layer that starts at 10 µm: constant
  // beyond its ends, it graded nothing (4.9× off).
  const two = neutral({ muL: 0, muR: -RT });
  two.regions = [{ material: 'm', length: L }, { material: 'm', length: L, grade: { A: { mu0: { x: [0, L], values: [0, 5 * RT] } } } }];
  two.interfaces = [{}];
  assert.ok(new Device(two).solve().warnings.some((w) => /regions\[1\]\.grade \(A\.mu0\): the profile runs from x = 0\.00e\+0 to 1\.00e-5 m, but region 1 is from 1\.00e-5/.test(w)));
});

test('graded ε under Debye–Hückel that reads the material’s ε is refused', () => {
  // (Poisson took the graded ε and the activities the material's: ln γ the same at both ends.)
  const def = {
    T,
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { w: { epsr: 78.4, species: { 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } }, statistics: [{ type: 'debye-huckel', species: ['K+', 'Cl-'], a: 4e-10 }] } },
    regions: [{ material: 'w', length: 1e-6, grade: { epsr: { x: [0, 1e-6], values: [78.4, 20] } } }],
    contacts: { left: { bath: { c: { 'K+': 100, 'Cl-': 100 } } }, right: { bath: { c: { 'K+': 100, 'Cl-': 100 } } } },
  };
  assert.throws(() => new Device(def), /grade\.epsr: 'w' has Debye–Hückel statistics, which read the material's epsr; give them their own/);
  def.materials.w.statistics[0].epsr = 78.4;
  assert.ok(new Device(def));
});
