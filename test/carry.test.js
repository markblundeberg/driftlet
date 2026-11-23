// What set() carries across a change of definition (src/carry.js): by what each value means,
// never by where it sat. Each test holds one family of silent wrong answers from Astra's review of
// set() (OpenAI, 17 Nov 2026), confirmed by four helpers' scripts; the old number is in each
// comment, the right one from a second route (an analytic value, a fresh device at the new
// settings, the same change split, or an untouched twin).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous, bath, build, half, hodgkinHuxley, layer, metal, ohmic, pulse, semiconductor } from '../src/kit.js';
import { sgInterpolate } from '../src/carry.js';
import { gateRate } from '../src/boundaries.js';

const rel = (a, b) => Math.abs(a / b - 1);
const worst = (xs, f) => Math.max(...Array.from(xs, f));
const L = 1e-6, RT300 = GAS_CONSTANT * 300, mV = 1e-3;

// --- shared fixtures
// A neutral 1 µm film of X (and Y), closed unless given contacts; statistics on request.
const film = ({ statistics, c0 = { X: 50 }, cells = 8, species = ['X'], contacts = {}, ...extra } = {}) => ({
  species: species.map((name) => ({ name, z: 0 })),
  materials: { film: { epsr: 4, species: Object.fromEntries(species.map((s) => [s, { D: 1e-9, cRef: 1, mu0: 0 }])), ...(statistics ? { statistics } : {}) } },
  regions: [{ material: 'film', length: L, c0 }],
  contacts,
  grid: { minCells: cells },
  ...extra,
});
const lattice = (species, cMax) => [{ type: 'lattice', species, cMax }];
// X held at 1 (left) and 10 (right): steady c = 1 + 9x/L.
const held = (muRight = RT300 * Math.log(10)) => ({
  left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: { type: 'equilibrium', mu: 0 } } },
  right: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: { type: 'equilibrium', mu: muRight } } },
});
// A 10 nm oxide between a capacitive plate (left) and a pinned one.
const oxide = (V = 1, extra = {}) => ({
  species: [{ name: 'X', z: 0 }], materials: { oxide: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'ox', material: 'oxide', length: 10e-9, c0: { X: 1 } }],
  contacts: { left: { V, zeroCharge: 0, phi: { type: 'capacitive', C: 0.01 } }, right: { V: 0, zeroCharge: 0, phi: 'pinned' } },
  grid: { minCells: 8 },
  ...extra,
});
// Hodgkin–Huxley gates on a capacitive membrane face between two baths of water.
const HH = {
  m: { alpha: { type: 'expLinear', rate: 1000, midpoint: -40 * mV, scale: 10 * mV }, beta: { type: 'exp', rate: 4000, midpoint: -65 * mV, scale: -18 * mV } },
  h: { alpha: { type: 'exp', rate: 70, midpoint: -65 * mV, scale: -20 * mV }, beta: { type: 'sigmoid', rate: 1000, midpoint: -35 * mV, scale: 10 * mV } },
  n: { alpha: { type: 'expLinear', rate: 100, midpoint: -55 * mV, scale: 10 * mV }, beta: { type: 'exp', rate: 125, midpoint: -65 * mV, scale: -80 * mV } },
};
const S = { 'Na+': 100, 'K+': 100, 'Cl-': 200 }, water3 = aqueous(['Na+', 'K+', 'Cl-'], { epsr: 0 });
const membrane = (gates, extra = {}) => ({
  phi: { type: 'capacitive', C: 0.01 }, gates,
  species: { 'K+': { type: 'permeability', P: 1e-8, gates: { n: 4 } }, 'Na+': gates.m ? { type: 'permeability', P: 1e-8, gates: { m: 3, h: 1 } } : 'blocked', 'Cl-': 'blocked' },
  ...extra,
});
const cell = (faces, V) => build({
  T: 273.15 + 6.3, library: [water3],
  stack: [bath(S, 'Cl-'), ...faces.flatMap((f) => [layer('water', 1e-6, { c0: S }), f]), layer('water', 1e-6, { c0: S }), bath(S, 'Cl-', V)],
  grid: { hmin: 50e-9, hmax: 200e-9 },
});
const gatesOf = (s, k = 0) => s.interfaces[k].gates;
// A squid axon, sealed, its HH membrane a port over 0..to of 6 mm.
const axon = (to, h = 2e-4) => {
  const T = 273.15 + 18.5, VT = (GAS_CONSTANT * T) / FARADAY, am = 2 / 238e-6;
  const z = { 'Na+': 1, 'K+': 1, 'Cl-': -1, 'A-': -1 }, E = { 'Na+': 50 * mV, 'K+': -77 * mV, 'Cl-': -54.387 * mV };
  const IN = { 'Na+': 50, 'K+': 400, 'Cl-': 52.5, 'A-': 397.5 };
  const OUT = Object.fromEntries(['Na+', 'K+', 'Cl-'].map((s) => [s, IN[s] * Math.exp((z[s] * E[s]) / VT)]));
  const G = { 'Na+': 1200, 'K+': 360, 'Cl-': 3 }, offset = (s) => (VT / z[s]) * Math.log(OUT[s] / 1000);
  const sealed = { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' };
  const sp = (D) => ({ D, mu0: 0, cRef: 1000 });
  return {
    T, species: Object.entries(z).map(([name, zz]) => ({ name, z: zz })),
    materials: { axoplasm: { epsr: 0, species: { 'Na+': sp(2e-9), 'K+': sp(2e-9), 'Cl-': sp(2e-9), 'A-': sp(1e-11) } } },
    regions: [{ name: 'axon', material: 'axoplasm', length: 6e-3, c0: IN }], contacts: { left: sealed, right: sealed },
    ports: [{
      name: 'membrane', region: 'axon', from: 0, to, V: 0, terminal: 'K+', area: am, capacitance: { C: 0.01, zeroCharge: 0 },
      gates: hodgkinHuxley({ inside: IN, outside: { ...OUT, 'A-': 1 }, T }).gates,
      species: { 'Na+': { type: 'conductance', G: G['Na+'] * am, offset: offset('Na+'), gates: { m: 3, h: 1 } }, 'K+': { type: 'conductance', G: G['K+'] * am, offset: offset('K+'), gates: { n: 4 } }, 'Cl-': { type: 'conductance', G: G['Cl-'] * am, offset: offset('Cl-') } },
    }],
    grid: { hmin: h, hmax: h },
  };
};
// A ferrocene couple tethered to gold in perchloric acid.
const acidLib = { species: [{ name: 'H+', z: 1, cRef: 1000 }, { name: 'ClO4-', z: -1, cRef: 1000 }], materials: { acid: { epsr: 0, species: { 'H+': { D: 9.311e-9, mu0: 0 }, 'ClO4-': { D: 1.792e-9, mu0: -8.52e3 } } } } };
const acid = { 'H+': 1000, 'ClO4-': 1000 }, Fc = (mu0 = 0) => ({ mu0, capacity: 1e-6, theta0: 0.3 }), FcClO4 = { mu0: -10e3, capacity: 1e-6, theta0: 0.1 };
const couple = (surface = { Fc: Fc(), FcClO4 }) => ({ surface, reactions: [{ equation: 'Fc + ClO4- = FcClO4 + e-', k0: 1e-6, alpha: 0.5 }] });
const tethered = (V = 0, face = couple()) => build({ T: 298.15, library: [acidLib, metal('Au')], stack: [ohmic(V, ['e-']), layer('Au', 1e-6), face, layer('acid', 10e-6, { c0: acid }), bath(acid, 'ClO4-')], grid: { hmin: 1e-8, hmax: 0.5e-6 } });
const cov = (s) => s.interfaces[0].coverage;
// Electrode ports with an adsorbate S in a 1 µm cell of water.
const surfaced = (ports, T = 298.15, hmax = 50e-9) => ({
  T, species: [{ name: 'A+', z: 1, cRef: 1000 }, { name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }],
  materials: { water: { epsr: 0, species: { 'A+': { D: 1e-9, mu0: 0 }, 'Na+': { D: 1e-9, mu0: 0 }, 'Cl-': { D: 1e-9, mu0: 0 } } } },
  regions: [{ material: 'water', length: 1e-6, c0: { 'A+': 10, 'Na+': 100, 'Cl-': 110 } }],
  contacts: { left: { V: 0, terminal: 'Na+', species: { 'Na+': 'equilibrium', 'A+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' }, right: {} },
  ports, grid: { hmin: hmax, hmax },
});
const adsorbing = (name, from, to, V) => ({ name, region: 0, from, to, V, terminal: 'e-', area: 1e6, surface: { S: { mu0: 2000, capacity: 1e-5, theta0: 1e-3 } }, reactions: [{ equation: 'A+ + e- = S', k0: 1e-4, alpha: 0.5 }] });
const theta = (s, name) => s.ports.find((p) => p.name === name).coverage.S;
// A Li⁺ bar with a capacitive port 'a' (state.test.js's capacitor, windowed).
const bar = () => ({
  species: [{ name: 'Li+', z: 1, cRef: 1000 }],
  materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
  regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * 96485.33212 }],
  contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
  ports: [{ name: 'a', region: 'bar', from: 1e-6, to: 3e-6, V: 0.05, area: 1e6, capacitance: { C: 0.2 } }],
  grid: { hmax: 0.25e-6 },
});

test('carry: a remesh with as many nodes as before carries the profile by place, not by index', () => {
  // (The mirrored mesh had 12 nodes too, and the vector was copied by index: c off by up to 3.62
  // of 10 against 1 + 9x/L; a cell moved between two regions, 4+6 → 5+5, put 0.665 at a node
  // where 4.61 is right.)
  const gA = { minCells: 4, hLeft: 4e-8, hRight: 1.5e-7 }, gB = { minCells: 4, hLeft: 1.5e-7, hRight: 4e-8 };
  const d = new Device(film({ T: 300, contacts: held(), c0: {}, grid: gA }));
  d.solve();
  const n0 = d.solver.u.length;
  const s = d.set({ grid: gB }).solution();
  assert.equal(d.solver.u.length, n0); // (the equal-length case)
  const err = worst(s.x, (x, g) => Math.abs(s.c.X[g] - (1 + (9 * x) / L)));
  assert.ok(err < 1e-6, `max |c − (1 + 9x/L)| ${err}`);
  const two = (n0, n1) => ({
    T: 300, species: [{ name: 'X', z: 0 }],
    materials: { a: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } }, b: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 1000 } } } },
    regions: [{ material: 'a', length: L / 2, grid: { minCells: n0, hmax: L / 2 / n0 } }, { material: 'b', length: L / 2, grid: { minCells: n1, hmax: L / 2 / n1 } }],
    interfaces: [{ species: { X: 'equilibrium' }, dipole: 0 }], contacts: held(),
  });
  const e = new Device(two(4, 6));
  e.solve();
  const t = e.set({ regions: two(5, 5).regions }).solution(), f = new Device(two(5, 5)).solve();
  const off = worst(t.c.X, (c, g) => Math.abs(c - f.c.X[g]));
  assert.ok(off < 1e-6, `4+6 → 5+5: max |c − fresh| ${off}`);
});

// A closed film of A and B (90 and 10), with A ⇌ B or without.
const AB = (rx, cells = 8) => ({
  species: [{ name: 'A', z: 0 }, { name: 'B', z: 0 }],
  materials: { film: { epsr: 4, species: { A: { D: 1e-9, cRef: 1, mu0: 0 }, B: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ material: 'film', length: L, c0: { A: 90, B: 10 } }], contacts: {},
  bulkReactions: rx ? [{ equation: 'A = B', kf: { film: 1 } }] : [], grid: { minCells: cells },
});
const amounts = (s) => Object.fromEntries(s.conservation.filter((c) => c.species).map((c) => [c.species, c.amount]));

test('carry: a reaction removed leaves each species what it holds now', () => {
  // (Only A + B was conserved before; set() kept each one's construction amount, and the next
  // solve() went back to 90/10 from 50/50: A 9e-5, B 1e-5 mol/m², converged.)
  const d = new Device(AB(true));
  d.solve();
  const a = amounts(d.set({ bulkReactions: [] }).solve()), fresh = amounts(new Device({ ...AB(false), regions: [{ material: 'film', length: L, c0: { A: 50, B: 50 } }] }).solve());
  for (const k of ['A', 'B']) assert.ok(rel(a[k], 5e-5) < 1e-9 && rel(a[k], fresh[k]) < 1e-9, `${k}: ${a[k]} vs ${fresh[k]}`);
});

test('carry: a regrid keeps a running reaction where it is, and a dissolved–adsorbed couple’s total', () => {
  // (Refined while A ⇌ B ran, the film read 90/10 again at once and the transient restarted
  // from there: 8.92e-5 of A after 10 ms, where 5e-5 holds.)
  const d = new Device(AB(true));
  d.solve();
  d.step(1e-6);
  d.set({ grid: { minCells: 16 } });
  for (const s of [d.solution(), d.step(1e-9)]) for (const [k, v] of Object.entries(amounts(s))) assert.ok(rel(v, 5e-5) < 1e-9, `${k}: ${v}`);
  // A closed X(aq) ⇌ X(ads) couple, refined after adsorption: the bulk took back its construction
  // amount and the couple held 1.880e-6 mol/m², +86%, for good (c0 gives 1.010e-6).
  const lib = { species: [{ name: 'X', z: 0, cRef: 1 }], materials: { water: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } } };
  const ads = (cells) => build({
    library: [lib, metal('Au')],
    stack: [ohmic(0, ['e-']), layer('Au', 1e-7), { surface: { Xads: { mu0: -10e3, capacity: 1e-6, theta0: 0.01 } }, reactions: [{ equation: 'X = Xads', k0: 1e-6 }] }, layer('water', 1e-6, { c0: { X: 1 }, grid: { minCells: cells } }), { species: { X: 'blocked' } }],
  });
  const total = (s) => s.conservation.find((c) => c.species === 'X').amount + s.interfaces[0].coverage.Xads * 1e-6;
  const e = new Device(ads(8));
  e.solve();
  const s = e.set({ regions: e.def.regions.map((r) => (r.material === 'water' ? { ...r, grid: { minCells: 16 } } : r)) }).solve();
  const f = new Device(ads(16)).solve();
  assert.ok(rel(total(s), 1.01e-6) < 1e-9 && rel(cov(s).Xads, cov(f).Xads) < 1e-6, `total ${total(s)}, θ ${cov(s).Xads} vs fresh ${cov(f).Xads}`);
});

test('carry: statistics switched to ideal, or a species taken off a lattice, keeps its concentration', () => {
  // (The lattice's activity was read as an ideal concentration: X 50 → 100, a closed amount
  // doubled at the next step; taken off a shared lattice, 30 → 75.)
  const d = new Device(film({ statistics: lattice(['X'], 100) }));
  d.solve();
  const s = d.set({ materials: { film: { statistics: [] } } }).step(1e-6);
  assert.ok(s.converged && rel(s.c.X[4], 50) < 1e-12 && rel(s.conservation[0].amount, 5e-5) < 1e-12, `${s.c.X[4]}`);
  const e = new Device(film({ species: ['X', 'Y'], c0: { X: 30, Y: 30 }, statistics: lattice(['X', 'Y'], 100) }));
  e.solve();
  const t = e.set({ materials: { film: { statistics: lattice(['Y'], 100) } } }).step(1e-6);
  assert.ok(rel(t.c.X[4], 30) < 1e-12 && rel(t.c.Y[4], 30) < 1e-12, `X ${t.c.X[4]}, Y ${t.c.Y[4]}`);
});

test('carry: species sharing a lattice keep their amounts on a new grid, in any listing order', () => {
  // (Regrid and cMax 100 → 200 together: X 30, Y 40 read 37.7/46.2, or 33.7/50.3 with the
  // species listed the other way; a nonuniform pair, regrid only: X's amount −8.4e-4; an open
  // lattice profile's inserted node 14.2 where the fine grid's steady state has 20.97.)
  for (const order of [['X', 'Y'], ['Y', 'X']]) {
    const d = new Device(film({ species: order, c0: { X: 30, Y: 40 }, statistics: lattice(order, 100) }));
    d.solve();
    d.step(1e-9);
    const s = d.set({ grid: { minCells: 16 }, materials: { film: { statistics: lattice(order, 200) } } }).step(1e-6);
    assert.ok(rel(s.c.X[3], 30) < 1e-12 && rel(s.c.Y[3], 40) < 1e-12, `${order}: X ${s.c.X[3]}, Y ${s.c.Y[3]}`);
  }
  const d = new Device(film({ species: ['X', 'Y'], c0: { X: { x: [0, L], values: [5, 60] }, Y: { x: [0, L], values: [80, 20] } }, statistics: lattice(['X', 'Y'], 100) }));
  const was = amounts(d.step(1e-12)), now = amounts(d.set({ grid: { minCells: 16 } }).solution());
  for (const k of ['X', 'Y']) assert.ok(rel(now[k], was[k]) < 1e-12, `${k}: ${now[k] / was[k] - 1}`);
  const open = (cells) => film({ cells, statistics: lattice(['X'], 100), contacts: { left: { species: { X: { type: 'equilibrium', mu: GAS_CONSTANT * 298.15 * Math.log(5.263) } }, phi: 'neutral' }, right: { species: { X: { type: 'equilibrium', mu: GAS_CONSTANT * 298.15 * Math.log(1900) } }, phi: 'neutral' } } });
  const e = new Device(open(8));
  e.solve();
  e.step(1e-15);
  const s = e.set({ grid: { minCells: 16 } }).solution(), f = new Device(open(16)).solve();
  const off = worst(s.c.X, (c, g) => rel(c, f.c.X[g]));
  assert.ok(off < 1e-9, `open lattice: worst |c/c_fine − 1| ${off}`);
});

test('carry: a set() that throws leaves the device as it was', () => {
  // (The carried 50 overflowed the new cMax 25 after the new definition was installed: the
  // device kept cMax 25, c0 10 and a clock at 0, and its next step read 20 where its twin has 50.)
  const make = () => {
    const d = new Device(film({ statistics: lattice(['X'], 100) }));
    d.solve();
    d.step(1e-6);
    return d;
  };
  const d = make(), twin = make(), solver = d.solver, def = JSON.stringify(d.def);
  assert.throws(() => d.set({ regions: [{ material: 'film', length: L, c0: { X: 10 } }], materials: { film: { statistics: lattice(['X'], 25) } } }), /site density/);
  assert.ok(d.solver === solver && JSON.stringify(d.def) === def, 'the definition and solver are the old ones');
  const a = d.step(1e-6), b = twin.step(1e-6);
  assert.ok(a.time === b.time && a.c.X[4] === b.c.X[4], `t ${a.time} vs ${b.time}, c ${a.c.X[4]} vs ${b.c.X[4]}`);
});

test('carry: a face’s reaction removed ahead of its gate leaves the gate, and blocks carry by meaning', () => {
  // (The widest face fixed the block width, so the slots kept their length and were copied by
  // index: face 0's gate n read its old reaction rate, 3.9e-9 where 0.0022 holds.)
  const face0 = (rx) => membrane({ n: HH.n }, rx === undefined ? {} : { reactions: [{ equation: 'Na+(left) = Na+(right)', k0: rx }] });
  const V = pulse({ start: 1e-4, width: 1, amplitude: -100 * mV, base: -130 * mV });
  const run = (patch) => {
    const d = new Device(cell([face0(1e-9), membrane(HH)], V));
    d.solve();
    d.advance(3e-4);
    d.set({ interfaces: [patch, membrane(HH)] });
    return gatesOf(d.solution()).n;
  };
  const removed = run(face0()), negligible = run(face0(1e-30)); // (the documented way to keep a reaction's state)
  assert.ok(rel(removed, negligible) < 1e-9 && negligible > 1e-3, `${removed} vs ${negligible}`);
  // The same vector length from other blocks (6+6 cells by a link → 3+4 and a face reaction):
  // the old steady profile, 1 + 9x/L, read 1 1 3.25 1 1 6.25 1 8.5 1.
  const reg = (cells) => ({ material: 'film', length: L / 2, grid: { minCells: cells, hmax: L / 2 / cells } });
  const d = new Device(film({ T: 300, contacts: held(), c0: {}, regions: [reg(6), reg(6)], interfaces: [{ species: { X: 'equilibrium' }, dipole: 0 }], grid: undefined }));
  d.solve();
  const n0 = d.solver.u.length, s = d.set({ regions: [reg(3), reg(4)], interfaces: [{ species: { X: 'blocked' }, dipole: 0, reactions: [{ equation: 'X(left) = X(right)', k0: 1e-3 }] }] }).solution();
  assert.equal(d.solver.u.length, n0);
  const err = worst(s.x, (x, g) => Math.abs(s.c.X[g] - (1 + (9 * x) / L)));
  assert.ok(err < 1e-9, `max |c − (1 + 9x/L)| ${err}`);
  // (and two regions to one at the same length threw a TypeError after committing)
  d.set({ regions: [{ ...reg(10), grid: { minCells: 10, hmax: L / 20 } }], interfaces: [] });
  const t = d.solve(), f = new Device(d.def).solve();
  assert.ok(worst(t.c.X, (c, g) => Math.abs(c - f.c.X[g])) < 1e-9);
});

test('carry: a terminal added keeps every other terminal’s charge', () => {
  // (Port 'b' added with 'a's C doubled: 'a's stored charge was recomputed at the new C, and the
  // next step failed at 9.3e-19 s, Q_a −0.00732 where its −0.00366 is kept at I = 0; a floating
  // plate's charge doubled at the next step, 0.00523 where 0.00262 is.)
  const b = { name: 'b', region: 'bar', from: 6e-6, to: 8e-6, V: 0, area: 1e6, capacitance: { C: 0.05 } };
  const a1 = { ...bar().ports[0], V: undefined, I: 0 }, a2 = { ...a1, capacitance: { C: 0.4 } };
  const charged = () => {
    const d = new Device(bar());
    d.solve();
    d.set({ ports: { a: { I: 0 } } }).advance(1e-4);
    return d;
  };
  const Qa = (s) => s.ports.find((p) => p.name === 'a').charge;
  const one = charged(), split = charged(), q0 = Qa(one.solution());
  one.set({ ports: [a2, b] });
  split.set({ ports: [a2] }).set({ ports: [a2, b] });
  const s = one.advance(1e-4 + 1e-9), t = split.advance(1e-4 + 1e-9);
  // (V to 1e-7: the window's boxes are put back to neutrality once on one route, twice on the other)
  assert.ok(s.converged && rel(Qa(s), q0) < 1e-9 && rel(s.terminals.a.V, t.terminals.a.V) < 1e-7, `Q_a ${Qa(s)} vs ${q0}, V_a ${s.terminals.a.V} vs ${t.terminals.a.V}`);
  const p = { name: 'p', region: 'ox', from: 4e-9, to: 6e-9, V: 0, area: 1e6, capacitance: { C: 1e-6 } };
  const d = new Device(oxide());
  d.solve();
  d.set({ contacts: { left: { I: 0 } } }).step(1e-9);
  const q = d.solution().plates.left.charge, u = d.set({ contacts: { left: { phi: { C: 0.02 } } }, ports: [p] }).step(1e-12);
  assert.ok(u.converged && rel(u.plates.left.charge, q) < 1e-9, `${u.plates.left.charge} vs ${q}`);
});

test('carry: a surface μ° changed moves no adsorbate; T and a regrid keep the couple whole', () => {
  // (The surface level was kept and read under the new μ°: Fc 0.00696 → 0.00466 at zero time,
  // Fc + FcClO₄ 0.4 → 0.3986, kept by the next solve(), which read 0.004659 where a device built
  // so has 0.004676; T 298.15 → 310 with a regrid copied the levels unconverted: the couple 0.364.)
  const d = new Device(tethered());
  const s0 = d.solve();
  const s1 = d.set({ interfaces: [couple({ Fc: Fc(1000), FcClO4 })] }).solution();
  assert.ok(rel(cov(s1).Fc, cov(s0).Fc) < 1e-12 && rel(cov(s1).Fc + cov(s1).FcClO4, 0.4) < 1e-12, `${cov(s1).Fc}`);
  const s2 = d.solve(), f = new Device(tethered(0, couple({ Fc: Fc(1000), FcClO4 }))).solve();
  assert.ok(rel(cov(s2).Fc, cov(f).Fc) < 1e-6, `${cov(s2).Fc} vs fresh ${cov(f).Fc}`);
  const e = new Device(tethered());
  e.solve();
  e.step(1e-12);
  const t = e.set({ T: 310, grid: { hmin: 0.5e-8, hmax: 0.25e-6 } }).step(1e-9);
  assert.ok(t.converged && rel(cov(t).Fc + cov(t).FcClO4, 0.4) < 1e-12, `${cov(t).Fc + cov(t).FcClO4}`);
});

test('carry: a species newly present starts at its c0, on a new grid or the same one', () => {
  // (Its slots were interpolated from an old solver that never had it: X 133 where c0 is 7;
  // charged, on the same grid, a level from the cold start's φ: 0.14 to 7; newly on a lattice
  // shared with a carried species, 1.8 to 9.6.)
  const water = { epsr: 78, species: { 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } };
  const def = (z, statistics, waterStats) => ({
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'X', z, cRef: 1000 }],
    materials: { water: { ...water, ...(waterStats ? { statistics: waterStats } : {}) }, water2: { epsr: 78, species: { ...water.species, X: { D: 1e-9, mu0: 5e3 } }, ...(statistics ? { statistics } : {}) } },
    regions: [{ name: 'w', material: 'water', length: 1e-7, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0 }, right: { V: 0.1, phi: 'pinned', zeroCharge: 0 } },
    grid: { hmin: 0.2e-9, hmax: 2e-9 },
  });
  const withX = { name: 'w', material: 'water2', length: 1e-7, c0: { 'K+': 100, 'Cl-': 100, X: 7 } };
  for (const [tag, d0, grid] of [
    ['neutral, regridded', def(0), { hmin: 0.1e-9, hmax: 1e-9 }],
    ['charged, same grid', def(1)],
    ['on a shared lattice', def(0, lattice(['K+', 'X'], 300), lattice(['K+'], 300))],
  ]) {
    const d = new Device(d0);
    d.solve();
    d.step(1e-9);
    const s = d.set({ regions: [withX], ...(grid ? { grid } : {}) }).solution();
    const off = worst(s.c.X, (c) => Math.abs(c - 7));
    assert.ok(off < 1e-9, `${tag}: max |X − 7| ${off}`);
  }
});

test('carry: a read before the first drive patch doesn’t change the first step', () => {
  // (A bare solution() made the patch take the drive shortcut, which kept the V = 0 cold start:
  // the first 100 ns read −73846 A/m² where the series capacitance passes +26154.)
  const eps0 = 8.8541878128e-12, Cs = 1 / (10e-9 / (4 * eps0) + 1 / 0.01);
  const first = (read) => {
    const d = new Device(oxide(0));
    if (read) d.solution();
    return d.set({ contacts: { left: { V: 1 } } }).step(100e-9).terminals.left.current;
  };
  const a = first(true), b = first(false);
  assert.ok(rel(a, b) < 1e-9 && rel(a, Cs / 100e-9) < 1e-6, `${a} vs ${b}, ${Cs / 100e-9}`);
});

// Two Pt electrodes behind Stern layers and diffuse layers across 10 µm of KCl.
const stern = (epsr) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr } });
const kcl = (V, epsr = 78.4) => build({ library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')], stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), stern(epsr), layer('water', 10e-6, { c0: { 'K+': 100, 'Cl-': 100 } }), stern(epsr), layer('Pt', 1e-6), ohmic(0, ['e-'])], grid: { hmin: 0.5e-6, hmax: 0.5e-6 } });
const inventory = (d, s, name) => {
  let a = 0;
  s.c[name].forEach((c, g) => Number.isFinite(c) && (a += d.model.grid.vol[g] * c));
  for (const itf of s.interfaces) for (const l of Object.values(itf.diffuse ?? {})) a += l.excess[name];
  return a;
};

test('carry: a diffuse layer’s ε or T changed keeps each closed ion’s total, layers and all', () => {
  // (The drop was kept and Γ recomputed under the new law, and taken as the step's start: at 1 V,
  // K⁺ and Cl⁻ each +2.35e-4 of their closed 1e-3 mol/m², kept through 1 ms; T 298.15 → 330,
  // −9.7e-5.)
  for (const patch of [(d) => ({ interfaces: d.def.interfaces.map((i) => (i.diffuse ? { ...i, diffuse: { epsr: 156.8 } } : i)) }), () => ({ T: 330 })]) {
    const d = new Device(kcl(1));
    d.solve();
    d.step(1e-9);
    d.set(patch(d));
    for (const s of [d.solution(), d.step(1e-9)]) for (const k of ['K+', 'Cl-']) assert.ok(rel(inventory(d, s, k), 1e-3) < 1e-9, `${JSON.stringify(patch(d)).slice(0, 20)} ${k}: ${inventory(d, s, k)}`);
  }
});

test('carry: gates by name, a widened gated window, surface species reordered', () => {
  // (Gates went by slot: listed n, h, m, the face read n 0.207 and m 0.334, swapped, and its Na⁺
  // flux 10 µs on was 3.5×. A gated window widened: the new nodes' gates at 0. A face's surface
  // species reordered: coverages reset to θ0, 4.03 A/m² from nothing; and ports reordered with a
  // T change moved their coverages, 0.7576 → 0.7632.)
  const V = pulse({ start: 1e-4, width: 1, amplitude: 60 * mV, base: -65 * mV });
  const run = (gates) => {
    const d = new Device(cell([membrane(HH)], V));
    d.solve();
    d.advance(1.5e-4);
    d.set({ materials: { water: { species: { 'Cl-': { D: water3.materials.water.species['Cl-'].D * (1 + 1e-9) } } } }, interfaces: [membrane(gates)] });
    const now = gatesOf(d.solution());
    return { now, Na: d.advance(1.6e-4).interfaces[0].N['Na+'] };
  };
  const same = run(HH), reversed = run({ n: HH.n, h: HH.h, m: HH.m });
  for (const k of ['m', 'h', 'n']) assert.ok(rel(reversed.now[k], same.now[k]) < 1e-12, `${k}: ${reversed.now[k]} vs ${same.now[k]}`);
  assert.ok(rel(reversed.Na, same.Na) < 1e-9, `Na⁺ ${reversed.Na} vs ${same.Na}`);
  const d = new Device(axon(3e-3));
  const p0 = d.advance(5e-3).ports[0], p = d.set({ ports: { membrane: { to: 6e-3 } } }).solution().ports[0];
  assert.ok(p.x.length > p0.x.length);
  // (the old part keeps what it held, its window uniform but for its edge; the new part, new
  // membrane, starts uncharged, φ under it jumping to the clamp's (a gate's hardware changed:
  // src/carry.js), its gates at rest at the voltage across it, α/(α + β): 0 V)
  for (const k of ['m', 'h', 'n']) {
    const lo = Math.min(...p0.gates[k]), hi = Math.max(...p0.gates[k]);
    p.x.forEach((x, w) => {
      if (x <= p0.x.at(-1)) assert.ok(p.gates[k][w] >= lo * (1 - 1e-12) && p.gates[k][w] <= hi * (1 + 1e-12), `${k} at ${x}: ${p.gates[k][w]} vs ${lo}..${hi}`);
    });
  }
  for (const gt of d.solver.gateList) {
    if (d.solver.model.grid.x[gt.g] <= p0.x.at(-1)) continue;
    const V = d.solver._gateVoltage(gt), a = gateRate(gt.gate.alpha, V)[0], b = gateRate(gt.gate.beta, V)[0];
    assert.ok(rel(d.solver.u[gt.o], a / (a + b)) < 1e-12 && Math.abs(V) < 1e-6, `${gt.gate.name} at ${d.solver.model.grid.x[gt.g]}: ${d.solver.u[gt.o]} vs ${a / (a + b)} at ${V} V`);
  }
  const e = new Device(tethered(0.15));
  e.solve();
  const before = cov(e.step(1e-9));
  const s = e.set({ materials: { acid: { species: { 'H+': { D: 9.311e-9 * (1 + 1e-9) } } } }, interfaces: [couple({ FcClO4, Fc: Fc() })] }).solution();
  for (const k of ['Fc', 'FcClO4']) assert.ok(rel(cov(s)[k], before[k]) < 1e-9, `${k}: ${cov(s)[k]} vs ${before[k]}`);
  const P = adsorbing('p', 0.1e-6, 0.4e-6, -0.05), Q = adsorbing('q', 0.6e-6, 0.9e-6, 0.1);
  const g = new Device(surfaced([P, Q]));
  const t0 = g.solve(), t = g.set({ T: 310, ports: [Q, P] }).solution();
  for (const nm of ['p', 'q']) assert.ok(rel(theta(t, nm)[1], theta(t0, nm)[1]) < 1e-9, `${nm}: ${theta(t, nm)[1]} vs ${theta(t0, nm)[1]}`);
});

test('carry: a regrid reads a surface’s or a gate’s state from within its own window only', () => {
  // (Interpolated across the window's edge from the empty slots beside it: an edge coverage 0.541
  // where 0.758 is uniform, the adsorbed amount −0.91%, taken as the new reference; a gated
  // window's last nodes, n 0.216 and 0.144 where 0.287 is uniform.)
  const d = new Device(surfaced([adsorbing('p', 0.06e-6, 0.94e-6, -0.05)], 298.15, 100e-9));
  d.solve();
  const th0 = theta(d.step(1e-9), 'p')[1], amt0 = d.solver.surfaceAmount(0);
  const s = d.set({ grid: { hmin: 25e-9, hmax: 25e-9 } }).solution();
  assert.ok(worst(theta(s, 'p'), (v) => rel(v, th0)) < 1e-9 && rel(d.solver.surfaceAmount(0), amt0) < 1e-12, `${Array.from(theta(s, 'p')).slice(0, 2)}; ${d.solver.surfaceAmount(0) / amt0 - 1}`);
  const e = new Device(axon(3.09e-3));
  e.advance(5e-3);
  const n0 = e.step(1e-6).ports[0].gates.n.at(-1), n = e.set({ grid: { hmin: 0.05e-3, hmax: 0.05e-3 } }).solution().ports[0].gates.n;
  assert.ok(worst(n, (v) => rel(v, n0)) < 1e-6, `${Array.from(n).slice(-3)} vs ${n0}`); // (uniform within 4e-9 of itself)
});

test('carry: a flowing cell refined follows the convection–diffusion profile', () => {
  // (Interpolated by diffusion alone, the last old cell's midpoint read 5.328 where
  // 1 + 8(e^{Pe x/L} − 1)/(e^{Pe} − 1) gives 3.292; the next outflow 0.0770 where 0.0200 holds.)
  const v = 0.02, Pe = (v * L) / 1e-9, exact = (x) => 1 + (8 * Math.expm1((Pe * x) / L)) / Math.expm1(Pe);
  const flow = (cells) => ({ ...film({ cells, c0: { X: 1 }, contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: RT300 * Math.log(9) } } } } }), T: 300, regions: [{ material: 'film', length: L, velocity: v, c0: { X: 1 } }] });
  const d = new Device(flow(8));
  d.solve();
  d.step(1e-12);
  const s = d.set({ grid: { minCells: 16 } }).solution();
  const err = worst(s.x, (x, g) => rel(s.c.X[g], exact(x)));
  assert.ok(err < 1e-6, `worst |c/exact − 1| ${err}`);
  const J = d.step(1e-9).contacts.right.flux.X;
  assert.ok(rel(J, 0.02) < 1e-6, `${J}`);
});

test('carry: a spherical shell refined follows a + b/r', () => {
  // (Weighted by distance, not ∫dx/A: an inserted node 4.077 where a + b/r gives 5.706.)
  const r0 = 1e-6, r1 = 10e-6, b = (1 - 9) / (1 / r0 - 1 / r1), a = 1 - b / r0, exact = (x) => a + b / (r0 + x);
  const shell = (cells) => ({ ...film({ cells, c0: { X: 1 }, contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: RT300 * Math.log(9) } } } } }), T: 300, geometry: { type: 'spherical', r0 }, regions: [{ material: 'film', length: r1 - r0, c0: { X: 1 } }] });
  const d = new Device(shell(4));
  d.solve();
  d.step(1e-12);
  const s = d.set({ grid: { minCells: 8 } }).solution();
  const err = worst(s.x, (x, g) => rel(s.c.X[g], exact(x)));
  assert.ok(err < 1e-6, `worst |c/(a + b/r) − 1| ${err}`);
});

test('carry: a strictly neutral transient refined keeps its closed ions', () => {
  // (Each species interpolated on its own left charge, and neutrality was reimposed at the next
  // step: the closed Ca²⁺ gained 9.35e-5 of itself, kept.)
  const Lg = 10e-6, sp = (D) => ({ D, mu0: 0, cRef: 1000 }), prof = (p, q) => ({ x: [0, Lg], values: [p, q] });
  const d = new Device({
    species: [{ name: 'Ca++', z: 2 }, { name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { gel: { epsr: 0, species: { 'Ca++': sp(0.8e-9), 'Na+': sp(1.3e-9), 'Cl-': sp(2e-9) } } },
    regions: [{ material: 'gel', length: Lg, c0: { 'Ca++': prof(1, 60), 'Na+': prof(150, 30), 'Cl-': prof(152, 150) } }],
    contacts: { left: { V: 0, terminal: 'Cl-', species: { 'Cl-': 'equilibrium', 'Ca++': 'blocked', 'Na+': 'blocked' }, phi: 'bulk' }, right: {} },
    grid: { minCells: 8 },
  });
  d.advance(2e-3);
  const ca = amounts(d.set({ grid: { minCells: 32 } }).step(1e-7))['Ca++'];
  assert.ok(rel(ca, 30.5 * Lg) < 1e-10, `${ca} vs ${30.5 * Lg}`); // (∫c0 dx: (1 + 60)/2 · L)
});

test('carry: regions swapped and one renamed: the survivor keeps its own, the new one starts from its c0 and says so', () => {
  // (The patch passed reorder detection and was carried by place: 'rich' read 10 where it held 50,
  // and the renamed region took 'rich's 50, unwarned.)
  const sp = (D) => ({ D, mu0: 0, cRef: 1000 }), sealed = { V: 0, terminal: 'Cl-', species: { 'K+': 'blocked', 'Cl-': 'equilibrium', 'A-': 'blocked' }, phi: 'bulk' };
  const poor = { name: 'poor', material: 'gel', length: 5e-6, c0: { 'K+': 110, 'Cl-': 100, 'A-': 10 } }, rich = { name: 'rich', material: 'gel', length: 5e-6, c0: { 'K+': 150, 'Cl-': 100, 'A-': 50 } };
  const d = new Device({
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }, { name: 'A-', z: -1 }],
    materials: { gel: { epsr: 0, species: { 'K+': sp(1e-9), 'Cl-': sp(1e-9), 'A-': sp(1e-10) } } },
    regions: [poor, rich], interfaces: [{ species: { 'A-': 'blocked' }, phi: 'neutral' }], contacts: { left: sealed, right: sealed }, grid: { minCells: 20 },
  });
  d.solve();
  const s = d.set({ regions: [rich, { ...poor, name: 'poor2', c0: { 'K+': 120, 'Cl-': 100, 'A-': 20 } }] }).solution();
  const A = (r) => s.c['A-'][s.region.indexOf(r) + 3];
  assert.ok(Math.abs(A(0) - 50) < 1e-6 && Math.abs(A(1) - 20) < 1e-6, `rich ${A(0)}, poor2 ${A(1)}`);
  assert.ok(s.warnings.some((w) => w.includes("'poor2'") && w.includes('c0')), `${s.warnings}`);
});

test('carry: a new cross-section turns the stored charges into totals', () => {
  // (Its per-m² displacement was read as a total in A: the first step after the patch passed
  // −2615 A where a device built with the area passes 1e-18; floated, its solve() failed.)
  const geometry = { area: { x: [0, 10e-9], values: [1e-6, 1e-6] } };
  const d = new Device(oxide());
  d.solve();
  const I = d.set({ geometry }).step(1e-6).terminals.left.current;
  const t = new Device(oxide(1, { geometry }));
  t.solve();
  const It = t.step(1e-6).terminals.left.current;
  assert.ok(Math.abs(I - It) < 1e-12, `${I} vs ${It} A`);
  const e = new Device(oxide());
  e.solve();
  const s = e.set({ contacts: { left: { I: 0 } } }).set({ geometry }).solve();
  const u = new Device(oxide(1, { geometry }));
  u.solve();
  const q = u.set({ contacts: { left: { I: 0 } } }).solve().plates.left.charge;
  assert.ok(s.converged && rel(s.plates.left.charge, q) < 1e-9, `${s.converged}: ${s.plates.left.charge} vs ${q}`);
});

test('carry: a call that does nothing, or fails, leaves c0 the start', () => {
  // (step(0), advance(0), a refused step, a failed solve, a rolled-back solve each marked the
  // device as run: a c0 patch then read the old 50 where 5 is given.)
  const patch = (d) => d.set({ regions: [{ ...d.def.regions[0], c0: { X: 5 } }] }).solution().c.X[3];
  const twin = new Device(film());
  twin.solution();
  const right = patch(twin);
  assert.equal(right, 5);
  for (const [tag, call] of Object.entries({
    'step(0)': (d) => d.step(0),
    'advance(0)': (d) => d.advance(0),
    'step(-1), refused': (d) => assert.throws(() => d.step(-1)),
    'a failed solve': (d) => ((d.solver.solveSteady = () => ({ converged: false, iterations: 0 })), assert.ok(!d.solve().converged)),
    'a rolled-back solve': (d) => {
      const cp = d._checkpoint();
      d.solve();
      d._rollback(cp);
    },
  })) {
    const d = new Device(film());
    call(d);
    assert.equal(patch(d), right, tag);
  }
});

test('carry: an empty patch, or the same drive again, changes nothing', () => {
  // (Taken as a drive change: the BDF2 history emptied, 70.785 A/m² where the untouched run steps
  // 70.038; and the steady flag dropped, so a regrid after it wasn't settled: outflow 0.0770
  // where 0.0200 holds.)
  const electrode = () => {
    const aq = aqueous(['Fe3+', 'Fe2+', 'K+', 'Cl-'], { epsr: 0 }), au = metal('Au');
    return {
      species: [...aq.species, { name: 'e-', z: -1 }], materials: { ...aq.materials, ...au.materials },
      regions: [{ name: 'Au', material: 'Au', length: 1e-6 }, { name: 'sol', material: 'water', length: 10e-6 }],
      interfaces: [{ reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 1e-3, alpha: 0.5 }] }],
      contacts: { left: { V: 0.02, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { bath: { c: { 'Fe3+': 10, 'Fe2+': 10, 'K+': 1000, 'Cl-': 1050 }, reference: 'Cl-' } } },
      grid: { hmax: 0.5e-6 },
    };
  };
  const run = (patch) => {
    const d = new Device(electrode());
    d.solve();
    d.set({ contacts: { left: { V: 0.04 } } });
    for (let k = 0; k < 3; k++) d.step(1e-6, { method: 'bdf2' });
    if (patch) d.set(patch);
    return d.step(1e-6, { method: 'bdf2' }).terminals.left.current;
  };
  const untouched = run(null);
  for (const patch of [{}, { contacts: { left: { V: 0.04 } } }]) assert.ok(rel(run(patch), untouched) < 1e-12, `${JSON.stringify(patch)}: ${run(patch)} vs ${untouched}`);
  const flow = (cells) => ({ ...film({ cells, c0: { X: 1 }, contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: RT300 * Math.log(9) } } } } }), T: 300, regions: [{ material: 'film', length: L, velocity: 0.02, c0: { X: 1 } }] });
  const d = new Device(flow(8));
  d.solve();
  const J = d.set({}).set({ grid: { minCells: 16 } }).step(1e-9).contacts.right.flux.X;
  assert.ok(rel(J, 0.02) < 1e-6, `${J}`);
});

test('carry: a transient landing on a waveform’s step reads as its last step saw it', () => {
  // (advance() to a jump of V 0 → 1 read the far side: V 1, plate 0.01 C/m² and 72214 A/m² where
  // nothing had flowed; a drive-only set() there read V 1 too.)
  const tj = 1e-6, wave = { t: [0, tj, tj, 1], values: [0, 0, 1, 1] };
  const d = new Device(oxide(wave));
  d.solve();
  const a = d.advance(tj);
  const e = new Device(oxide(wave));
  e.solve();
  const b = e.advance(tj * (1 - 1e-6)); // (just before the step: its twin)
  assert.ok(a.terminals.left.V === b.terminals.left.V && a.plates.left.charge === b.plates.left.charge && Math.abs(a.terminals.left.current) < 1e-6, `V ${a.terminals.left.V}, plate ${a.plates.left.charge}, I ${a.terminals.left.current}`);
  assert.equal(d.set({ contacts: { right: { V: 1e-9 } } }).solution().terminals.left.V, 0);
});

test('carry: a named port patch with a typed link replaces the link whole', () => {
  // (Deep-merged: { type: 'conductance', G } over one with an offset kept the 5 mV offset, the
  // drain current halved, −5.0e-5 A/m² where −1.0e-4 is right; to 'equilibrium', the stale G threw.)
  const def = (link) => ({ ...structuredClone(bar()), ports: [{ name: 'drain', region: 'bar', V: 0, terminal: 'Li+', species: { 'Li+': link } }], contacts: { left: { V: 0.01, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { phi: 'neutral' } }, grid: { minCells: 40 } });
  const drain = (d) => d.solve().terminals.drain.current;
  for (const [from, to] of [[{ type: 'conductance', G: 1e3, offset: 0.005 }, { type: 'conductance', G: 1e3 }], [{ type: 'conductance', G: 1e3 }, { type: 'equilibrium' }]]) {
    const d = new Device(def(from));
    d.solve();
    d.set({ ports: { drain: { species: { 'Li+': to } } } });
    assert.deepEqual(d.def.ports[0].species['Li+'], to);
    const I = drain(d), fresh = drain(new Device(def(to)));
    assert.ok(rel(I, fresh) < 1e-9, `${JSON.stringify(to)}: ${I} vs ${fresh}`);
  }
});

test('carry: a contacts container the constructor refuses, set() refuses too', () => {
  // (The drive shortcut stored contacts: [] or 7, and the next rebuild threw.)
  for (const bad of [[], 7]) {
    const d = new Device(film());
    d.solve();
    assert.throws(() => d.set({ contacts: bad }), /contacts must be an object/);
    assert.deepEqual(d.def.contacts, {});
    assert.ok(d.set({ grid: { minCells: 9 } }).solve().converged);
  }
});

test('carry: a steep cell refined neither rounds nor underflows', () => {
  // (Weights e^{Δθ w} rounded at 8 V over four cells: inserted nodes 1.25e-23 where 1e-6 is held
  // throughout, the next current half of a fresh fine device's; at 30 V, 4e-70, and the step failed.)
  const ox = (cells, VL) => ({
    species: [{ name: 'X+', z: 1, cRef: 1e-6 }], materials: { ox: { epsr: 4, species: { 'X+': { D: 1e-12, mu0: 0 } } } },
    regions: [{ material: 'ox', length: 100e-9, c0: { 'X+': 1e-6 } }],
    contacts: { left: { V: VL, terminal: 'X+', species: { 'X+': 'equilibrium' }, phi: 'pinned', zeroCharge: 0 }, right: { V: 0, terminal: 'X+', species: { 'X+': 'equilibrium' }, phi: 'pinned', zeroCharge: 0 } },
    grid: { minCells: cells, hmax: 100e-9 / cells },
  });
  for (const VL of [8, 30]) {
    const d = new Device(ox(4, VL));
    d.solve();
    d.step(1e-9);
    const s = d.set({ grid: { minCells: 8, hmax: 100e-9 / 8 } }).solution();
    assert.ok(worst(s.c['X+'], (c) => rel(c, 1e-6)) < 1e-5, `${VL} V: ${Array.from(s.c['X+']).slice(0, 3)}`);
    const t = new Device(ox(8, VL));
    t.solve();
    t.step(1e-9);
    const a = d.step(1e-7), b = t.step(1e-7);
    assert.ok(a.converged && rel(a.current, b.current) < 1e-3, `${VL} V: ${a.current} vs ${b.current}`);
  }
});

test('carry: a layer resized upstream carries the next one’s profile with it', () => {
  // (Carried by absolute x, the downstream tracer was read 1 µm upstream of itself: its closed
  // amount 5.5e-6 → 1.0e-5 mol/m², +82%, with no drift reported.)
  const L1 = 1e-6, regions = [{ name: 'up', material: 'film', length: 1e-6, c0: { X: 1 } }, { name: 'down', material: 'film', length: L1, c0: { X: { x: [1e-6, 1e-6 + L1], values: [1, 10] } } }];
  const d = new Device({
    T: 300, species: [{ name: 'X', z: 0 }], materials: { film: { epsr: 4, species: { X: { D: 1e-12, cRef: 1, mu0: 0 } } } },
    regions, interfaces: [{ species: { X: 'blocked' }, dipole: 0 }],
    contacts: { left: held(0).left, right: { V: 0, phi: 'pinned', zeroCharge: 0 } }, grid: { hmin: 0.1e-6, hmax: 0.1e-6 },
  });
  const downstream = (s) => Array.from(s.c.X).filter((_, g) => s.region[g] === 1);
  const was = downstream(d.step(1e-6)); // (D 1e-12: the tracer relaxes over ~1 s)
  const now = downstream(d.set({ regions: [{ ...regions[0], length: 2e-6 }, regions[1]] }).solution());
  const err = worst(now, (c, k) => rel(c, was[k])), down = d.step(1e-6).conservation.find((c) => c.species === 'X' && c.regions[0] === 1).amount;
  assert.ok(err < 1e-12 && rel(down, 5.5e-6) < 1e-12, `max off ${err}, amount ${down}`);
});

test('carry: a diffuse-layer face at a driven current steps after a small change of T', () => {
  // (The step after set({ T: T + 0.1 }) failed, even shortened to 1e-15 s; solve() first worked.)
  const ferri = (T) => {
    const c = { FeIII: 5, FeII: 5, 'Cl-': 1000, 'K+': 1035 };
    return build({
      T, library: [{ species: [{ name: 'FeIII', z: -3, cRef: 1000 }, { name: 'FeII', z: -4, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
        materials: { water: { epsr: 0, species: { FeIII: { D: 0.72e-9, mu0: 0 }, FeII: { D: 0.67e-9, mu0: -FARADAY * 0.36 }, 'K+': { D: 1.957e-9, mu0: 0 }, 'Cl-': { D: 2.032e-9, mu0: 0 } } } } }, metal('Pt')],
      stack: [ohmic({ I: 10 }, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 }, reactions: [{ ...half('FeIII + e- = FeII'), k0: 5e-2, alpha: 0.5 }] }, layer('water', 30e-6, { c0: c }), bath(c, 'Cl-')],
      grid: { hmin: 30e-9, hmax: 2e-6, ratio: 1.2 },
    });
  };
  const run = (dT) => {
    const d = new Device(ferri(298.15));
    d.solve();
    d.step(1e-9);
    if (dT) d.set({ T: 298.15 + dT });
    return d.step(1e-12);
  };
  const a = run(0.1), b = run(0);
  assert.ok(a.converged, a.warnings[0]);
  assert.ok(Math.abs(a.terminals.left.V - b.terminals.left.V) < 1e-3, `${a.terminals.left.V} vs ${b.terminals.left.V}`); // (0.1 K moves it by ~RT/F · 0.1/298)
});

// --- Astra's second pass (17 Nov 2026): ~45 more wrong answers in set(), confirmed by helpers'
// scripts (~/.cache/driftlet-claude/astra2-chase/verify). Each test is one family, its pointer
// number first; the old number is in its comment, the right one from a second route.
// Needs, beyond the imports above:
//   import { sgInterpolate } from '../src/carry.js';
//   import { semiconductor } from '../src/kit.js';   (add to the kit import)

// A peaked closed profile over a region of length Ld: ∫ = 13.375 Ld/1 (trapezoids on eight cells).
const peak = (Ld) => ({ X: { x: [0, Ld / 8, Ld / 4, Ld], values: [1, 100, 1, 1] } });

test('carry: 1 · a strictly neutral region with fixed charge, regridded, keeps its ions finite and its closed totals', () => {
  // (Neutrality was put back with ρ_fixed in C/m³ added to z·c in mol/m³: every inserted Li⁺ node
  // read Infinity and φ −2479 V; a closed K⁺ against −10 C/m³ of fixed charge lost 9.84% of itself
  // at the next step, kept through 10 ms.)
  const d = new Device({ ...bar(), ports: [] });
  d.solve();
  const s = d.set({ grid: { hmax: 0.125e-6 } }).solution(), f = new Device({ ...bar(), ports: [], grid: { hmax: 0.125e-6 } }).solve();
  assert.equal(s.x.length, f.x.length);
  const off = worst(s.c['Li+'], (c, g) => rel(c, f.c['Li+'][g]));
  assert.ok(off < 1e-12, `Li⁺: worst |c/fresh − 1| ${off}`);
  const sp = (D) => ({ D, mu0: 0, cRef: 1000 }), rho = -10;
  const gel = (cells) => ({
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { gel: { epsr: 0, species: { 'K+': sp(1e-9), 'Cl-': sp(2e-9) } } },
    regions: [{ material: 'gel', length: 10e-6, fixedCharge: rho, c0: { 'K+': 100 - rho / FARADAY, 'Cl-': 100 } }],
    contacts: { left: { V: 0, terminal: 'Cl-', species: { 'Cl-': 'equilibrium', 'K+': 'blocked' }, phi: 'bulk' }, right: {} }, grid: { minCells: cells },
  });
  const run = (refine) => {
    const e = new Device(gel(16));
    e.solve();
    e.set({ contacts: { left: { V: 0.01 } } }).step(1e-4);
    if (refine) e.set({ grid: { minCells: 32 } });
    return amounts(e.step(1e-6))['K+'];
  };
  const k = run(true), twin = run(false);
  assert.ok(rel(k, twin) < 1e-12, `closed K⁺ ${k} vs ${twin} unrefined`);
});

test('carry: 2 · species reordered keep their own concentrations and fluxes', () => {
  // (Node slots and face fluxes were copied by index: X read 7 and Y 2; an open X held at 1 and 10
  // read 3 3 3 with no flux, where 1 … 10 carries −0.009.)
  const d = new Device(film({ species: ['X', 'Y'], c0: { X: 2, Y: 7 } }));
  d.solve();
  const s = d.set({ species: [...d.def.species].reverse() }).solution();
  assert.ok(rel(s.c.X[3], 2) < 1e-12 && rel(s.c.Y[3], 7) < 1e-12, `X ${s.c.X[3]}, Y ${s.c.Y[3]}`);
  const h = held();
  h.left.species.Y = h.right.species.Y = 'blocked';
  const e = new Device(film({ T: 300, species: ['X', 'Y'], c0: { X: 1, Y: 3 }, contacts: h }));
  const t0 = e.solve(), t = e.set({ species: [...e.def.species].reverse() }).solution();
  assert.ok(rel(t.flux.X[4], t0.flux.X[4]) < 1e-9 && worst(t.c.Y, (c) => rel(c, 3)) < 1e-12 && worst(t.c.X, (c, g) => rel(c, t0.c.X[g])) < 1e-12, `flux X ${t.flux.X[4]} vs ${t0.flux.X[4]}`);
});

test('carry: 3 · a coarsening onto old nodes keeps the amount, and a nested one is settled', () => {
  // (Dyadic lengths: every new node an old one, so the carry skipped the amounts: 13.375 → 1, −92.5%,
  // its reference kept, until a solve(); and the nested 8 → 4 coarsening of a steady state wasn't
  // solved again: the first step's inflow 1.147484e-6 where the 4-cell steady state gives 1.147358e-6.)
  const one = { material: 'film', length: 1, c0: peak(1) };
  const d = new Device(film({ regions: [one] }));
  const a0 = amounts(d.step(1e-12)).X, s = d.set({ grid: { minCells: 4 } }).solution(), c = s.conservation[0];
  assert.ok(rel(a0, 13.375) < 1e-12 && rel(c.amount, a0) < 1e-12 && rel(c.reference, a0) < 1e-12, `${a0} → ${c.amount} (ref ${c.reference})`);
  const RT = GAS_CONSTANT * 300, end = (muX) => ({ V: 0, phi: 'pinned', zeroCharge: 0, species: { X: { type: 'equilibrium', mu: muX }, Y: { type: 'equilibrium', mu: 0 } } });
  const XY = (cells) => ({
    ...film({ T: 300, cells, species: ['X', 'Y'], c0: { X: 1, Y: 1 }, contacts: { left: end(RT * Math.log(10)), right: end(0) } }),
    regions: [{ material: 'film', length: 1, c0: { X: 1, Y: 1 } }], bulkReactions: [{ equation: 'X = Y', kf: { film: 1e-6 } }],
  });
  const e = new Device(XY(8));
  e.solve();
  const t = new Device(XY(4));
  t.solve();
  const J = e.set({ grid: { minCells: 4 } }).step(1e-3).contacts.left.flux.X, Jt = t.step(1e-3).contacts.left.flux.X;
  assert.ok(rel(J, Jt) < 1e-9, `${J} vs ${Jt} solved on 4`);
});

test('carry: 4, 13 · a face made neutral or blocked keeps no old displacement or flux', () => {
  // (The face's inactive slots were copied: made neutral, its old D −1.77e-5 C/m² stayed in Gauss's
  // law for good, the films' field 5e5 V/m where a neutral face leaves them flat; made blocked, its
  // N.X read −0.0045 after a converged solve, where the segment flux reads 0.)
  const two = (face, contacts) => film({ T: 300, c0: { X: 1 }, regions: ['a', 'b'].map((name) => ({ name, material: 'film', length: L, c0: { X: 1 } })), interfaces: [face], contacts });
  const pins = { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } }, right: { V: 1, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } } };
  const d = new Device(two({ phi: 'pinned', dipole: 0 }, pins));
  d.solve();
  const s = d.set({ interfaces: [{ phi: 'neutral' }] }).solve(), f = new Device(two({ phi: 'neutral' }, pins)).solve();
  const off = worst(s.phi, (p, g) => Math.abs(p - f.phi[g]));
  assert.ok(s.converged && off < 1e-9, `neutral face: max |φ − fresh| ${off} V`);
  const e = new Device(two({ species: { X: 'equilibrium' } }, held()));
  assert.ok(e.solve().interfaces[0].N.X < -1e-3);
  for (const t of [e.set({ interfaces: [{ species: { X: 'blocked' } }] }).solution(), e.solve()]) {
    assert.ok(t.interfaces[0].N.X === 0 && t.flux.X[e.grid.regionEnd[0]] === 0, `blocked face: N.X ${t.interfaces[0].N.X}`);
  }
});

test('carry: 6 · one fresh region in a stretch keeps the matched regions’ amounts', () => {
  // (Any fresh region dropped every target in its stretch: A coarsened 8 → 5 with B appended held
  // 1.592e-5 mol/m² where A's 1.3375e-5 and B's 7e-6 make 2.0375e-5, and took that as its reference.)
  const A = { name: 'A', material: 'film', length: L, c0: peak(L) }, B = { name: 'B', material: 'film', length: L, c0: { X: 7 } };
  const d = new Device(film({ regions: [A] }));
  const a0 = amounts(d.step(1e-12)).X;
  const c = d.set({ grid: { minCells: 5 }, regions: [A, B] }).solution().conservation[0];
  const e = new Device(film({ regions: [A] }));
  e.step(1e-12);
  const split = e.set({ regions: [A, B] }).set({ grid: { minCells: 5 } }).solution().conservation[0].amount;
  assert.ok(rel(c.amount, a0 + 7 * L) < 1e-12 && rel(c.reference, a0 + 7 * L) < 1e-12 && rel(split, a0 + 7 * L) < 1e-12, `${c.amount} (ref ${c.reference}), split ${split}, want ${a0 + 7 * L}`);
});

test('carry: 7 · a Debye–Hückel law reading ε from its material moves no concentration when ε changes', () => {
  // (The node was copied bit for bit, its activity read under the new ε: K⁺ and Cl⁻ 10 → 10.570 of
  // a closed 10, kept through 1 ms, where ε given in the law itself kept 10.)
  const salt = (epsr, inLaw) => ({
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { w: { epsr, species: { 'K+': { D: 1e-9, cRef: 10, mu0: 0 }, 'Cl-': { D: 1e-9, cRef: 10, mu0: 0 } }, statistics: [{ type: 'debye-huckel', species: ['K+', 'Cl-'], a: 4e-10, ...(inLaw ? { epsr } : {}) }] } },
    regions: [{ material: 'w', length: 1e-6, c0: { 'K+': 10, 'Cl-': 10 } }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { 'K+': 'blocked', 'Cl-': 'blocked' } }, right: { phi: 'neutral', species: { 'K+': 'blocked', 'Cl-': 'blocked' } } },
    grid: { minCells: 8 },
  });
  for (const inLaw of [false, true]) {
    const d = new Device(salt(80, inLaw));
    d.solve();
    d.set({ materials: { w: salt(60, inLaw).materials.w } });
    for (const s of [d.solution(), d.advance(1e-3)]) for (const k of ['K+', 'Cl-']) assert.ok(worst(s.c[k], (c) => rel(c, 10)) < 1e-12, `${inLaw ? 'ε in the law' : 'ε from the material'} ${k}: ${s.c[k][4]}`);
  }
});

test('carry: 8 · a strictly neutral region regridded keeps a port’s countercharge', () => {
  // (Neutrality ignored the port's spread charge aσ: 1832 C/m³ unbalanced at the new nodes, the
  // port discharged in one step, Q −6.9e-11 where −3.6626e-3 stays; the axon refined fired a spike
  // to +37.5 mV within 0.2 ms where its untouched twin sits at −65 mV.)
  const def = (hmax) => ({
    ...bar(), species: [...bar().species, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 }, 'Cl-': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, c0: { 'Li+': 100, 'Cl-': 100 } }], grid: { hmax },
  });
  const d = new Device(def(0.25e-6)), f = new Device(def(0.125e-6));
  for (const e of [d, f]) e.solve(), e.step(1e-12);
  const q = d.set({ grid: { hmax: 0.125e-6 } }).step(1e-12).ports[0].charge, qf = f.step(1e-12).ports[0].charge;
  assert.ok(rel(q, qf) < 1e-6, `port Q ${q} vs fresh fine ${qf}`);
  const a = new Device(axon(3e-3)), t = new Device(axon(3e-3));
  a.advance(5e-3), t.advance(5e-3);
  a.set({ grid: { hmin: 1e-4, hmax: 1e-4 } });
  const s = a.advance(5.1e-3), u = t.advance(5.1e-3);
  assert.ok(rel(s.ports[0].charge, u.ports[0].charge) < 1e-3 && Math.max(...s.phi) < -60 * mV, `axon: Q ${s.ports[0].charge} vs ${u.ports[0].charge}, max φ ${Math.max(...s.phi)}`);
});

test('carry: 9 · a sphere’s or cylinder’s centre (r0 = 0) carries its charges without NaN', () => {
  // (Its stored charges were scaled by the area ratio 0/0: the centre's running charge and its
  // next current read NaN after set({ T }).)
  for (const type of ['spherical', 'cylindrical']) {
    const def = (T) => film({ T, c0: { X: 1 }, geometry: { type, r0: 0 } });
    const d = new Device(def(300));
    d.solve();
    const I = d.set({ T: 310 }).step(1e-9).terminals.left.current, f = new Device(def(310));
    f.solve();
    assert.ok(Number.isFinite(I) && Math.abs(I - f.step(1e-9).terminals.left.current) < 1e-12, `${type}: ${I}`);
  }
});

test('carry: 10 · a double layer’s charge is held after the amounts are put back', () => {
  // (Its drop was found for the old charge, then the closed ions were rescaled at the edge, which
  // moved it: q −0.066199 → −0.066204 C/m² under the face's old D, and the first 1 ps passed
  // 4.66e6 A/m² where the next passes 8.30e3.)
  const d = new Device(kcl(0.8));
  d.solve();
  const q0 = d.step(1e-9).interfaces[0].diffuse.right.charge;
  const s = d.set({ interfaces: d.def.interfaces.map((i) => (i.diffuse ? { ...i, diffuse: { epsr: 313.6 } } : i)) }).solution();
  const q = s.interfaces[0].diffuse.right.charge, I1 = d.step(1e-12).current, I2 = d.step(1e-12).current;
  assert.ok(rel(q, q0) < 1e-12 && rel(I1, I2) < 1e-3, `q ${q} vs ${q0}; first current ${I1}, next ${I2}`);
});

test('carry: 11 · a diffuse layer removed puts its excess back in the bulk', () => {
  // (Only new layers were visited: the removed layers' excess went missing, each closed ion 0.99960
  // of its 1e-3 mol/m² until a solve().)
  const d = new Device(kcl(0.8));
  d.solve();
  d.step(1e-9);
  d.set({ interfaces: d.def.interfaces.map(({ diffuse, ...rest }) => rest) });
  for (const s of [d.solution(), d.step(1e-9)]) for (const k of ['K+', 'Cl-']) assert.ok(rel(inventory(d, s, k), 1e-3) < 1e-9, `${k}: ${inventory(d, s, k)}`);
});

test('carry: 12 · a waveform held “forever” by a far last knot still reads its near jump', () => {
  // (The before-side read backed off 8 ulps of the whole span, 1.8 ms for a knot at 1e12 s: at 2 µs,
  // 1 µs past a 0 → 1 V jump, the step and the snapshot read 0 V and no charge.)
  const at2us = (tLast) => {
    const d = new Device(oxide({ t: [0, 1e-6, 1e-6, tLast], values: [0, 0, 1, 1] }));
    d.solve();
    const s = d.step(2e-6);
    return [s.terminals.left.V, s.plates.left.charge, d.solution().terminals.left.V];
  };
  const far = at2us(1e12), near = at2us(1);
  assert.ok(far.every((v, k) => v === near[k]) && near[0] === 1, `${far} vs ${near}`);
});

test('carry: 14, 21d · species sharing sites keep each amount in either listing order', () => {
  // (Corrected one at a time, a later species moved the earlier ones: on a lattice 98.5% full, X
  // −9.15e-7 of itself listed first, Y +9.07e-7 listed first; two adsorbates on one window, S1
  // −2.29% listed first, exact listed second.)
  const prof = (a, b, c) => ({ x: [0, L / 2, L], values: [a, b, c] });
  for (const order of [['X', 'Y'], ['Y', 'X']]) {
    const d = new Device(film({ species: order, cells: 4, c0: { X: prof(99, 98, 99), Y: prof(0.99, 1.99, 0.99) }, statistics: lattice(order, 100) }));
    const was = amounts(d.step(1e-16)), now = amounts(d.set({ grid: { minCells: 5 } }).solution());
    for (const k of order) assert.ok(rel(now[k], was[k]) < 1e-12, `${order} ${k}: ${now[k] / was[k] - 1}`);
  }
  const want = { S1: [0.1, 0.1, 0.1, 0.75, 0.1, 0.1, 0.1], S2: [0.895, 0.895, 0.895, 0.05, 0.895, 0.895, 0.895] };
  for (const order of [['S1', 'S2'], ['S2', 'S1']]) {
    const port = { ...adsorbing('p', 0.1e-6, 0.4e-6, 0), surface: Object.fromEntries(order.map((nm) => [nm, { mu0: 2000, capacity: 1e-5, theta0: 0.1 }])), reactions: order.map((nm) => ({ equation: `A+ + e- = ${nm}`, k0: 1e-4, alpha: 0.5 })) };
    const d = new Device(surfaced([port]));
    d.step(1e-15);
    const s = d.solver, sf = s.model.surfaces[0];
    sf.nodes.forEach((g, w) => {
      const b = s.blockOfNode[g], bare = 1 - want.S1[w] - want.S2[w];
      sf.species.forEach((sp, q) => ((s.u[b * s.M + 1 + s.n + q] = sp.mu0 / s.model.RT + Math.log(want[sp.name][w] / bare)), (s.uLo[b * s.M + 1 + s.n + q] = 0)));
    });
    s.computeConcentrations();
    const amt = (v) => Object.fromEntries(v.surfCols.map(([k, q], j) => [v.model.surfaces[k].species[q].name, v.surfaceAmount(j)]));
    const before = amt(s), after = amt(d.set({ grid: { hmin: 70e-9, hmax: 70e-9 } }).solver);
    for (const k of order) assert.ok(rel(after[k], before[k]) < 1e-12, `${order} ${k}: ${after[k] / before[k] - 1}`);
  }
});

test('carry: 16 · a per-region grid settles and redrives like a top-level one', () => {
  // (regions[].grid wasn't a grid change to set(): 4 → 16 cells of a solved double layer wasn't solved
  // again, the first 0.1 ns passed 9.24e6 A/m² where nothing flows; with a bias step in the same
  // patch, 2.15e6 where refine, solve, redrive gives −7.09e6.)
  const def = (regionGrid, grid, V = 0.1) => ({
    T: 300, species: [{ name: 'A+', z: 1, cRef: 1000 }, { name: 'B-', z: -1, cRef: 1000 }],
    materials: { w: { epsr: 80, species: { 'A+': { D: 1e-9, mu0: 0 }, 'B-': { D: 1e-9, mu0: 0 } } } },
    regions: [{ name: 'e', material: 'w', length: 100e-9, c0: { 'A+': 1, 'B-': 1 }, ...(regionGrid ? { grid: regionGrid } : {}) }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { 'A+': 'blocked', 'B-': 'blocked' } }, right: { V, phi: 'pinned', zeroCharge: 0, species: { 'A+': 'blocked', 'B-': 'blocked' } } },
    grid,
  });
  const run = (perRegion, patch = {}) => {
    const d = new Device(perRegion ? def({ minCells: 4 }, {}) : def(undefined, { minCells: 4 }));
    d.solve();
    d.set({ ...(perRegion ? { regions: [{ ...d.def.regions[0], grid: { minCells: 16 } }] } : { grid: { minCells: 16 } }), ...patch });
    return d.step(1e-10);
  };
  const a = run(true), b = run(false);
  assert.ok(Math.abs(a.contacts.right.current - b.contacts.right.current) < 1e-3 && worst(a.phi, (p, g) => Math.abs(p - b.phi[g])) < 1e-12, `${a.contacts.right.current} vs ${b.contacts.right.current}`);
  const bias = { contacts: { right: { V: 0.2 } } }, c = run(true, bias), e = new Device(def(undefined, { minCells: 4 }));
  e.solve();
  e.set({ grid: { minCells: 16 } }).solve();
  const I = e.set(bias).step(1e-10).contacts.right.current;
  assert.ok(rel(c.contacts.right.current, I) < 1e-9, `with a bias step: ${c.contacts.right.current} vs ${I} refined, solved, redriven`);
});

test('carry: 17, 18 · a floating capacitor keeps its charge through a change of C, with or without a new geometry', () => {
  // (Its old voltage was kept and read under the new C: the plate's snapshot charge doubled, 2.615e-3
  // → 5.231e-3 C/m², until a step; the port's too; and with a new cross-section in the same patch the
  // doubled charge was kept for good, −0.014648 C where geometry and C apart give −0.0073238.)
  const d = new Device(oxide());
  d.solve();
  const q = d.set({ contacts: { left: { I: 0 } } }).solve().plates.left.charge;
  const s = d.set({ contacts: { left: { phi: { C: 0.02 } } } }).solution(), Cs = 1 / (1 / 0.02 + 10e-9 / (4 * 8.8541878128e-12));
  assert.ok(rel(s.plates.left.charge, q) < 1e-9 && rel(s.terminals.left.V, q / Cs) < 1e-9, `plate ${s.plates.left.charge} vs ${q}, V ${s.terminals.left.V} vs Q/C ${q / Cs}`);
  const floated = () => {
    const b = new Device(bar());
    b.solve();
    b.set({ ports: { a: { I: 0 } } }).solve();
    return b;
  };
  const C = { ports: { a: { capacitance: { C: 0.4 } } } }, geometry = { area: { x: [0, 10e-6], values: [2, 2] } };
  const b = floated(), qa = b.solution().ports[0].charge;
  assert.ok(rel(b.set(C).solution().ports[0].charge, qa) < 1e-6, `port ${b.solution().ports[0].charge} vs ${qa}`);
  const one = floated().set({ geometry, ...C }).solve().ports[0].charge, two = floated().set({ geometry });
  two.solve();
  const split = two.set(C).solve().ports[0].charge;
  assert.ok(rel(one, split) < 1e-9 && rel(split, 2 * qa) < 1e-8, `together ${one}, apart ${split}, 2·Q ${2 * qa}`);
});

test('carry: 22 · a region restarted mid-ramp starts at the sources of now', () => {
  // (Built at t = 0's sources: the fresh films started at φ 0 beside a plate held at 0.1 V, and the
  // next 1 ps passed 4.3e-7 A/m² where a device built at 0.1 V passes 2.6e8.)
  const two = { regions: ['a', 'b'].map((name) => ({ name, material: 'oxide', length: 5e-9, c0: { X: 1 } })), interfaces: [{ species: { X: 'equilibrium' }, dipole: 0 }] };
  const d = new Device(oxide({ t: [0, 1e-6], values: [0, 0.1] }));
  d.solve();
  d.advance(1e-6);
  const s = d.set(two).solution(), f = new Device({ ...oxide(0.1), ...two });
  assert.ok(worst(s.phi, (p, g) => Math.abs(p - f.solution().phi[g])) < 1e-12, `start φ ${s.phi[0]} vs ${f.solution().phi[0]}`);
  const I = d.step(1e-12).terminals.left.current, If = f.step(1e-12).terminals.left.current;
  assert.ok(rel(I, If) < 1e-6, `${I} vs ${If}`);
});

test('carry: 23 · held → floating, or rebuilt, at a waveform’s landing reads the side before the jump', () => {
  // (A rebuild reset the last step's length, so the snapshot read the jump's far side: V 1 and a
  // plate of 0.01 C/m² where nothing had flowed; switched to I = 0 there, it started floating at 1 V.)
  const tj = 1e-6, wave = { t: [0, tj, tj, 1], values: [0, 0, 1, 1] };
  const landed = (t = tj) => {
    const d = new Device(oxide(wave));
    d.solve();
    d.advance(t);
    return d;
  };
  const twin = landed(tj * (1 - 1e-9)).set({ contacts: { left: { I: 0 } } }).solution();
  for (const [tag, patch] of [['floated', { contacts: { left: { I: 0 } } }], ['rebuilt', { materials: { oxide: { species: { X: { D: 2e-9 } } } } }]]) {
    const s = landed().set(patch).solution();
    assert.ok(Math.abs(s.terminals.left.V - twin.terminals.left.V) < 1e-9 && Math.abs(s.plates.left.charge - twin.plates.left.charge) < 1e-11, `${tag}: V ${s.terminals.left.V}, plate ${s.plates.left.charge}`);
  }
});

test('carry: 25 · a face’s reactions are matched by everything that takes part', () => {
  // (Keyed by mobile participants and order alone, two pathways Fc + ClO₄⁻ = FcClO₄ or FcB + e⁻
  // swapped rates when listed the other way, 1.57e-7 ↔ 5.73e-6, until a step; the first removed,
  // the second inherited its rate.)
  const R1 = { equation: 'Fc + ClO4- = FcClO4 + e-', k0: 1e-6, alpha: 0.5 }, R2 = { equation: 'Fc + ClO4- = FcB + e-', k0: 1e-4, alpha: 0.5 };
  const surface = { Fc: Fc(), FcClO4, FcB: { mu0: -5e3, capacity: 1e-6, theta0: 0.1 } };
  const run = (patch) => {
    const d = new Device(tethered(0, { surface, reactions: [R1, R2] }));
    d.solve();
    d.set({ contacts: { left: { V: 0.25 } } }).advance(1e-4);
    if (patch) d.set(patch(d));
    const s = d.solution();
    return Object.fromEntries(d.def.interfaces[0].reactions.map((r, k) => [r.equation, s.interfaces[0].rates[k]]));
  };
  const unused = (d) => ({ ...d.def.materials, unused: { epsr: 1, species: {} } }); // (a rebuild, no physics)
  const same = run((d) => ({ materials: unused(d), interfaces: [{ ...d.def.interfaces[0], reactions: [R1, R2] }] }));
  const swapped = run((d) => ({ materials: unused(d), interfaces: [{ ...d.def.interfaces[0], reactions: [R2, R1] }] }));
  for (const r of [R1, R2]) assert.ok(rel(swapped[r.equation], same[r.equation]) < 1e-12, `${r.equation}: ${swapped[r.equation]} vs ${same[r.equation]}`);
  const { FcClO4: _, ...rest } = surface, removed = run(() => ({ interfaces: [{ surface: rest, reactions: [R2] }] }));
  assert.ok(rel(removed[R2.equation], same[R2.equation]) < 1e-12, `R2 alone ${removed[R2.equation]} vs ${same[R2.equation]}`);
});

test('carry: 26 · a metal turned semiconductor starts over, and says so', () => {
  // (The metal's flux slot was read as φ̂: the silicon started at φ 0 with n 46.5 mol/m³ where its
  // own start has 0.166, no warning, and every step failed.)
  const lib = [semiconductor('Si'), metal('Au')], Nd = 1e23 / 6.02214076e23; // (1e17 cm⁻³)
  const si = build({ T: 300, library: lib, stack: [ohmic(0), layer('Si', 1e-6, { donors: Nd }), ohmic(0.01)], grid: { minCells: 16 } });
  const d = new Device(build({ T: 300, library: lib, stack: [ohmic(0, ['e-']), layer('Au', 1e-6), ohmic(1e-9, ['e-'])], grid: { minCells: 16 } }));
  d.solve();
  const s = d.set({ materials: si.materials, regions: si.regions, contacts: si.contacts }).solution(), f = new Device(si).solution();
  assert.ok(worst(s.phi, (p, g) => Math.abs(p - f.phi[g])) < 1e-12 && worst(s.c['e-'], (c, g) => rel(c, f.c['e-'][g])) < 1e-12, `φ ${s.phi[8]} vs ${f.phi[8]}`);
  assert.ok(s.warnings.some((w) => w.includes('metal')), `${s.warnings}`);
  assert.ok(d.step(1e-15).converged);
});

test('carry: 27b · a clamped port window widened without a regrid keeps its charge', () => {
  // (The new part's charge was saved but its old part's taken as the whole: the membrane gave up
  // three quarters of its charge in 0.1 ms, Q 0.0078 where a device built with the wide window holds
  // 0.0330, φ −30.9 mV where it sits at −65.4. Then the new part was charged in an instant, by
  // nothing: 0.0165 C/m² put on it with no current through the clamp. Now the window widened is
  // new membrane, uncharged, which the clamp charges through the axoplasm over its RC time: the
  // membrane's charge is what it held at the set(), and 40 ms on, the device is the wide one's.)
  const d = new Device(axon(3e-3)), w = new Device(axon(6e-3));
  const q0 = d.advance(5e-3).ports[0].charge;
  w.advance(5e-3);
  // (to the round-off of the ions' net charge under it, 1e-16 of the ions', which are 1e6 its)
  assert.ok(rel(d.set({ ports: { membrane: { to: 6e-3 } } }).solution().ports[0].charge, q0) < 1e-8);
  assert.equal(d.solver.portQEnd[0], q0);
  const s = d.advance(40e-3), t = w.advance(40e-3);
  assert.equal(s.warnings.length, 0, s.warnings.join('\n'));
  assert.ok(rel(s.ports[0].charge, t.ports[0].charge) < 1e-4 && Math.abs(s.phi[10] - t.phi[10]) < 0.1 * mV, `Q ${s.ports[0].charge} vs ${t.ports[0].charge}, φ ${s.phi[10]} vs ${t.phi[10]}`);
});

test('carry: 30 · an immobile species or an eddy-mixed one is refined linearly', () => {
  // (Both were interpolated along the molecular drift: an immobile graded Y⁺ in 2.3 V read 1.976e-3
  // where its fixed profile has 1.5e-3, its first half's amount +7.0%, kept by solve(); a tracer
  // mixed at 100 D, 1.204 where the fine steady state has 1.455, first outflow −1.174 for −0.712.)
  const imm = (cells) => ({
    T: 300, species: [{ name: 'Y', z: 1 }], materials: { d: { epsr: 4, species: { Y: { D: 0, cRef: 1, mu0: 0 } } } },
    regions: [{ material: 'd', length: L, c0: { Y: { x: [0, L], values: [1e-3, 9e-3] } } }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { Y: 'blocked' } }, right: { V: 1, phi: 'pinned', zeroCharge: 0, species: { Y: 'blocked' } } },
    grid: { minCells: cells },
  });
  const d = new Device(imm(8));
  d.solve();
  const s = d.set({ grid: { minCells: 16 } }).solution(), off = worst(s.x, (x, g) => rel(s.c.Y[g], 1e-3 + (8e-3 * x) / L));
  assert.ok(off < 1e-12, `immobile: worst |c/c0 − 1| ${off}`);
  const RT = RT300, open = { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: RT * Math.log(9) } } } };
  const flow = (cells) => ({ ...film({ cells, c0: { X: 1 }, contacts: open }), T: 300, regions: [{ material: 'film', length: L, velocity: 0.02, c0: { X: 1 }, mixing: 1e-7 }] });
  const e = new Device(flow(8)), f = new Device(flow(16));
  e.solve();
  e.step(1e-12);
  const t = e.set({ grid: { minCells: 16 } }).solution(), u = f.solve(), mixOff = worst(t.c.X, (c, g) => rel(c, u.c.X[g]));
  const J = e.step(1e-9).contacts.right.flux.X, Jf = f.step(1e-9).contacts.right.flux.X;
  assert.ok(mixOff < 1e-3 && rel(J, Jf) < 1e-3, `mixed: worst |c/fine − 1| ${mixOff}, outflow ${J} vs ${Jf}`); // (8 cells' steady state against 16's: ~2e-4)
});

test('carry: 31 · an underflowed concentration stays 0, a full surface stays full', () => {
  // (A carried 0 took the new start's level: X behind a 2000 RT step read 1 after a 0.01 K change,
  // where 0 holds, and the next step failed; a surface at θ = 1 to the last bit was refused for
  // holding more than its sites.)
  const RT = GAS_CONSTANT * 300;
  const two = (T) => ({
    T, species: [{ name: 'X', z: 0 }],
    materials: { a: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } }, b: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 2000 * RT } } } },
    regions: [{ name: 'a', material: 'a', length: L, c0: { X: 1 } }, { name: 'b', material: 'b', length: L, c0: { X: 1 } }],
    interfaces: [{ species: { X: 'equilibrium' }, dipole: 0 }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: { type: 'equilibrium', mu: 0 } } }, right: { V: 0, phi: 'pinned', zeroCharge: 0 } },
    grid: { minCells: 4 },
  });
  const d = new Device(two(300));
  d.solve();
  const inB = (s) => Array.from(s.c.X).filter((_, g) => s.region[g] === 1);
  const s = d.set({ T: 300.01, materials: two(300.01).materials }).solution(), t = d.step(1e-9);
  assert.ok(inB(s).every((c) => c === 0) && t.converged && inB(t).every((c) => c === 0), `X in b ${inB(s)}`);
  const e = new Device(surfaced([{ ...adsorbing('p', 0.1e-6, 0.4e-6, 0), surface: { S: { mu0: -200e3, capacity: 1e-5, theta0: 0.5 } } }]));
  e.solve();
  assert.equal(theta(e.solution(), 'p')[0], 1);
  assert.equal(theta(e.set({ T: 298.2 }).solution(), 'p')[0], 1);
});

test('carry: 33, 34 · terminals’ charges scale with area one way, and still sum to nothing', () => {
  // (A restart kept the running charge per m² while a carry scaled it: species + 3 m² in one patch,
  // the next charge 1.046e-3 C where geometry then species gives 1.569e-3; and contacts' charges were
  // scaled where a port's wasn't, the terminals' sum −0.0398 C where a closed device passes 0.)
  const geometry = { area: { x: [0, 10e-9], values: [3, 3] } };
  const withY = (def) => ({ species: [...def.species, { name: 'Y', z: 0 }], materials: { oxide: { species: { Y: { D: 1e-9, cRef: 1, mu0: 0 } } } }, regions: [{ ...def.regions[0], c0: { X: 1, Y: 1 } }] });
  const charged = () => {
    const d = new Device(oxide(0));
    d.solve();
    d.set({ contacts: { left: { V: 0.1 } } }).advance(1e-6);
    return d;
  };
  const next = (d) => d.advance(d.solution().time + 1e-9).trace.terminals.left.charge.at(-1);
  const a = charged(), b = charged();
  a.set({ ...withY(a.def), geometry });
  b.set({ geometry }).set(withY(b.def));
  const qa = next(a), qb = next(b);
  assert.ok(rel(qa, qb) < 1e-9, `one patch ${qa} vs geometry then species ${qb}`);
  const d = new Device(bar());
  d.solve();
  d.set({ ports: { a: { V: 0.1 } } }).advance(1e-3);
  const tr = d.set({ geometry: { area: { x: [0, 10e-6], values: [3, 3] } } }).advance(1e-3 + 1e-6).trace.terminals;
  const [l, r, p] = ['left', 'right', 'a'].map((k) => tr[k].charge.at(-1));
  assert.ok(Math.abs(l + r + p) < 1e-9 * Math.abs(l) && Math.abs(l) > 0.05, `left ${l} + right ${r} + port ${p} = ${l + r + p}`);
});

test('carry: 36 · the device keeps a copy of the caller’s definition', () => {
  // (Kept by reference: an edit to it in place, then set() with it, compared equal to itself and did
  // nothing, D 1e-9 where 3e-9 is given; an edit with no set() came in at the next unrelated rebuild.)
  const f = film({ T: 300, contacts: held() }), d = new Device(f);
  d.solve();
  f.materials.film.species.X.D = 3e-9;
  assert.equal(d.set({ materials: f.materials }).solve().D.X[0], 3e-9);
  const g = film({ T: 300, contacts: held() }), e = new Device(g);
  e.solve();
  g.materials.film.species.X.D = 5e-9;
  assert.equal(e.set({ grid: { minCells: 9 } }).solve().D.X[0], 1e-9);
});

test('carry: 37 · a current read at t = 0 before any step is the drive at t = 0', () => {
  // (A jump 0 → 1 A/m² at t = 0 on an untouched device: solution() read 0 while the circuit row drove
  // 1, and a constant 1 read 1.)
  const read = (I) => {
    const d = new Device(oxide());
    return d.set({ contacts: { left: { I } } }).solution().terminals.left.current;
  };
  const jump = read({ t: [-1, 0, 0, 1], values: [0, 0, 1, 1] });
  assert.ok(jump === read(1) && jump === 1, `${jump}`);
});

test('carry: 38 · the SG weight keeps a small positive fraction', () => {
  // (Formed through 1 − w, w = 1e-17 rounded away: sgInterpolate(1, 1e20, 1, 1e-17) read 1 where
  // e^{−dw}(c_a + s(c_b e^d − c_a)), s = expm1(dw)/expm1(d), gives 1582.98; at 1e-16, +11%.)
  const ref = (ca, cb, d, w) => Math.exp(-d * w) * (ca + (Math.expm1(d * w) / Math.expm1(d)) * (cb * Math.exp(d) - ca));
  for (const w of [1e-17, 1e-16, 1e-12]) assert.ok(rel(sgInterpolate(1, 1e20, 1, w), ref(1, 1e20, 1, w)) < 1e-12, `w ${w}: ${sgInterpolate(1, 1e20, 1, w)} vs ${ref(1, 1e20, 1, w)}`);
});

test('carry: 40 · rescaling a closed species leaves an open partner on its lattice where it was', () => {
  // (X rescaled on sites it shares with Y, held by a contact: Y read 30.20 far from the peak where
  // it was 30, its amount −4.3%, in either listing.)
  const RT = GAS_CONSTANT * 298.15, prof = (...v) => ({ x: v.map((_, k) => (k * L) / (v.length - 1)), values: v });
  for (const order of [['X', 'Y'], ['Y', 'X']]) {
    const d = new Device(film({
      species: order, statistics: lattice(order, 100), c0: { X: prof(5, 90, 5, 5, 5, 5, 5, 5, 5), Y: prof(30, 5, 30, 30, 30, 30, 30, 30, 30) },
      contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { Y: { type: 'equilibrium', mu: RT * Math.log(30 / 65) }, X: 'blocked' } }, right: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked', Y: 'blocked' } } },
    }));
    const x0 = amounts(d.step(1e-12)).X, s = d.set({ grid: { minCells: 5 } }).solution();
    const off = worst(s.x, (x, g) => (x >= 0.4 * L ? rel(s.c.Y[g], 30) : 0));
    assert.ok(off < 1e-9 && rel(amounts(s).X, x0) < 1e-12, `${order}: worst |Y/30 − 1| past the peak ${off}`);
  }
});

test('carry: 41 · impedance settles a pending regrid, and a call that takes no step doesn’t', () => {
  // (impedance() solved the steady state but left the regrid pending, so the next drive patch solved
  // it again; advance(now) took no step yet settled it, the device then calling itself steady.)
  const regridded = () => {
    const d = new Device(oxide());
    d.solve();
    return d.set({ grid: { minCells: 13 } });
  };
  const d = regridded();
  d.impedance([10]);
  let solves = 0;
  const f = d.solver.solveSteady.bind(d.solver);
  d.solver.solveSteady = (...a) => (solves++, f(...a));
  d.set({ contacts: { right: { V: 0.01 } } });
  assert.ok(!d.solver.regriddedSteady && solves === 0, `pending ${d.solver.regriddedSteady}, solves ${solves}`);
  const e = regridded(), twin = regridded();
  e.advance(e.solution().time);
  assert.ok(e.solver.regriddedSteady && e.solution().steady === twin.solution().steady, `pending ${e.solver.regriddedSteady}, steady ${e.solution().steady}`);
});

test('carry: 42 · keys that join names can’t collide', () => {
  // (Joined with ':' and NUL: port 'a:b' gate 'c' and port 'a' gate 'b:c' were one gate, its window's
  // edge n 0.257 where 0.054 is right; species A, B and one named 'A:0\0B' looked the same, so the
  // new one took A's 2 where its c0 says 5.)
  const gated = (names) => {
    const def = axon(6e-3), [p] = def.ports, { n } = p.gates;
    const port = ([pn, gn], from, to, V) => ({ ...p, name: pn, from, to, V, gates: { [gn]: n }, species: { 'K+': { ...p.species['K+'], gates: { [gn]: 4 } } } });
    const d = new Device({ ...def, ports: [port(names[0], 0, 2.1e-3, 0), port(names[1], 2.2e-3, 6e-3, -0.09)] });
    d.advance(2e-3);
    d.set({ grid: { hmin: 1.5e-4, hmax: 1.5e-4 } });
    return d.solver.gateList.map((gt) => d.solver.u[gt.o]);
  };
  const a = gated([['a:b', 'c'], ['a', 'b:c']]), b = gated([['p', 'c'], ['q', 'd']]);
  assert.ok(a.length === b.length && worst(a, (v, k) => Math.abs(v - b[k])) === 0, `worst gate difference ${worst(a, (v, k) => Math.abs(v - b[k]))}`);
  const sp = (species, c0) => ({ ...film({ T: 300, species: [], c0 }), species, materials: { film: { epsr: 4, species: Object.fromEntries(species.map((s) => [s.name, { D: 1e-9, cRef: 1, mu0: 0 }])) } } });
  const d = new Device(sp([{ name: 'A', z: 0 }, { name: 'B', z: 0 }], { A: 2, B: 7 }));
  d.solve();
  const nul = 'A:0\u0000B';
  d._rebuild(sp([{ name: nul, z: 0 }], { [nul]: 5 }), d._solver); // (set() can't drop material film's A and B)
  assert.equal(d.solution().c[nul][3], 5);
});

test('carry: 45 · holes or trailing undefined in a patched array are refused as the constructor refuses them', () => {
  // (Arrays were compared by their keys alone, so new Array(2) equalled []: the patch returned as a
  // no-op, the definition kept, where new Device() refuses it.)
  const r = film().regions[0];
  for (const [base, patch, msg] of [[{ interfaces: [] }, { interfaces: new Array(2) }, /interfaces has 2 entries/], [{ interfaces: [] }, { interfaces: [undefined] }, /interfaces has 1 entries/], [{}, { regions: [r, undefined] }, /regions\[1\] must be an object/]]) {
    assert.throws(() => new Device(film({ ...base, ...patch })), msg);
    const d = new Device(film(base));
    d.solve();
    assert.throws(() => d.set(patch), msg);
  }
});

test('a near-full lattice refined keeps its amount and fits its sites (Astra, ninth pass)', () => {
  // A gated film's dopant in two valence states on shared sites, 99.6% full at the peak. Refined, the
  // amount was restored by multiplying Ox⁺ (4.2302 asked of 4.2289 sites) and set() refused it; a fresh
  // solve on the finer grid holds 4.2112. Now its odds are scaled, and the carried state solves to it.
  const L = 1.1438352041877807e-7, n = 0.34338815638381825, total = 0.7365310063519617, fraction = 0.6550324019044638;
  const cMax = 4.228913112594661, RT = 300 * GAS_CONSTANT, ox = total * fraction, red = total - ox, mu0 = RT * Math.log1p(-total / cMax);
  const def = (cells) => ({
    T: 300,
    species: [{ name: 'e-', z: -1 }, { name: 'Ox+', z: 1 }, { name: 'Red', z: 0 }],
    materials: {
      film: {
        epsr: 25.545087549835443,
        species: { 'e-': { D: 1.2272493014945945e-6, mu0: RT * Math.log(50), cRef: n * 50 }, 'Ox+': { D: 2.2801449272726107e-13, mu0, cRef: ox }, Red: { D: 3.878570472044756e-13, mu0, cRef: red } },
        statistics: [{ type: 'lattice', species: ['Ox+', 'Red'], cMax }],
      },
    },
    regions: [{ material: 'film', length: L, fixedCharge: FARADAY * (n - ox), c0: { 'e-': n, 'Ox+': ox, Red: red } }],
    contacts: {
      left: { V: 0, phi: 'pinned', zeroCharge: 0, terminal: 'e-', species: { 'e-': 'equilibrium', Red: { type: 'equilibrium', mu: 0 } } },
      right: { V: -0.1926582992170006, phi: 'pinned', zeroCharge: 0 },
    },
    bulkReactions: [{ equation: 'Ox+ + e- = Red', kf: { film: 199343.51796813446 } }],
    grid: { hmax: L / cells },
  });
  const coarse = new Device(def(12));
  assert.ok(coarse.solve().converged);
  coarse.set({ grid: { hmax: L / 24 } });
  const carried = coarse.solve(), fresh = new Device(def(24)).solve();
  assert.ok(carried.converged && fresh.converged);
  const peak = (s) => Math.max(...s.c.Red.map((v, k) => v + s.c['Ox+'][k]));
  assert.ok(peak(fresh) < cMax && Math.abs(peak(carried) - peak(fresh)) < 1e-12, `${peak(carried)} vs ${peak(fresh)}`);
  assert.ok(carried.c['Ox+'].every((v, k) => Math.abs(v - fresh.c['Ox+'][k]) < 1e-12 * cMax));
});

// A held gate over a region with ε > 0 (no strictly neutral group) whose capacitance changes: its
// plate's charge can't jump in the set() any more than a gate's over neutral boxes can, so it
// starts from what it held, and the first step passes the difference through its terminal,
// recorded. (Read off the carried state at the new C instead, the charge jumped unrecorded: a
// thin-film transistor's gate C doubled, its first 22 ms step's currents summed to −31 A/m² of
// 31, both contacts passing it out with the gate at −0.73; bench/combos' kirchhoff property,
// proposed in Neon's review of the suite's shape.)
test("carry: a held gate over a region with ε > 0, its C doubled: it starts from its charge, and the first step's currents and ledger sum to zero", () => {
  const def = bar();
  def.materials.solid.epsr = 10;
  const d = new Device(def);
  assert.ok(d.solve().converged);
  const k = d.solver.terms.findIndex((t) => t.name === 'a');
  const before = d.solver.portQ[0];
  d.set({ ports: { a: { capacitance: { C: 0.4 } } } });
  const s = d.solver;
  assert.equal(s.portQEnd[0], before, 'the gate starts from the charge it held');
  const q0 = Float64Array.from(s.termCharge);
  const st = d.step(1e-9);
  const I = s.terms.map((t) => st.terminals[t.name].current), gross = I.reduce((a, b) => a + Math.abs(b), 0);
  assert.ok(gross > 0, 'something flows');
  assert.ok(Math.abs(I.reduce((a, b) => a + b, 0)) <= 1e-6 * gross, `Kirchhoff: ${I.map((x) => x.toExponential(3)).join(', ')}`);
  const dq = Array.from(s.termCharge, (q, j) => q - q0[j]), moved = dq.reduce((a, b) => a + Math.abs(b), 0);
  assert.ok(Math.abs(dq.reduce((a, b) => a + b, 0)) <= 1e-6 * moved, `the ledger: ${dq.map((x) => x.toExponential(3)).join(', ')}`);
  // (to 1e-6: a step from rest passes the steady solve's round-off, 1e-6 A/m² here; 422 A/m²
  // pass, summing to 2.5e-6) And the gate's terminal recorded what its plate took: from what it held to what it holds)
  assert.ok(rel(dq[k], s.portQ[0] - before) < 1e-9, `the gate passed ${dq[k].toExponential(6)} C, its plate took ${(s.portQ[0] - before).toExponential(6)}`);
});
