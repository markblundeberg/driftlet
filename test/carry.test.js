// What set() carries across a change of definition (src/carry.js): by what each value means,
// never by where it sat. Each test holds one family of silent wrong answers from Astra's review of
// set() (OpenAI, 17 Nov 2026), confirmed by four helpers' scripts; the old number is in each
// comment, the right one from a second route (an analytic value, a fresh device at the new
// settings, the same change split, or an untouched twin).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous, bath, build, half, hodgkinHuxley, layer, metal, ohmic, pulse } from '../src/kit.js';
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
  assert.ok(s.converged && rel(Qa(s), q0) < 1e-9 && rel(s.terminals.a.V, t.terminals.a.V) < 1e-9, `Q_a ${Qa(s)} vs ${q0}, V_a ${s.terminals.a.V} vs ${t.terminals.a.V}`);
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
  // (the old part keeps what it held, its window uniform but for its edge; the new part starts at
  // rest at the voltage the carried state puts across it, α/(α + β))
  for (const k of ['m', 'h', 'n']) {
    const lo = Math.min(...p0.gates[k]), hi = Math.max(...p0.gates[k]);
    p.x.forEach((x, w) => {
      if (x <= p0.x.at(-1)) assert.ok(p.gates[k][w] >= lo * (1 - 1e-12) && p.gates[k][w] <= hi * (1 + 1e-12), `${k} at ${x}: ${p.gates[k][w]} vs ${lo}..${hi}`);
    });
  }
  for (const gt of d.solver.gateList) {
    if (d.solver.model.grid.x[gt.g] <= p0.x.at(-1)) continue;
    const V = d.solver._gateVoltage(gt), a = gateRate(gt.gate.alpha, V)[0], b = gateRate(gt.gate.beta, V)[0];
    assert.ok(rel(d.solver.u[gt.o], a / (a + b)) < 1e-12 && Math.abs(V + 0.0656) < 0.002, `${gt.gate.name} at ${d.solver.model.grid.x[gt.g]}: ${d.solver.u[gt.o]} vs ${a / (a + b)} at ${V} V`);
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
