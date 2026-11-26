// Round 4, Newton, the giant steps and the continuations (Astra's fourth pass, chased by probes):
// each test failed before its fix.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { aqueous, bath, build, half, layer, metal, ohmic, photogeneration, pulse, semiconductor } from '../src/kit.js';

const VT = (GAS_CONSTANT * 298.15) / FARADAY;

// A21c: a current that is zero at the steady solve passes nothing, a waveform's as much as a
// constant's. The open-circuit electrode feeding a closed Ag⁺ pool (test/steady.test.js), driven
// by a pulse not yet begun, kept the silver couple off its affinity row: Q/K 0.0475, 78 mV off.
test('a waveform current at zero at the solve passes nothing: the closed pool at equilibrium', () => {
  const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 }), iron = half('Fe3+ + e- = Fe2+');
  const bathC = { 'K+': 5000, 'Fe3+': 10, 'Fe2+': 10, 'NO3-': 5050 };
  const cell = (I) => new Device(build({
    library: [aqueous(['K+', 'NO3-', 'Fe3+', 'Fe2+', 'Ag+'], { epsr: 0 }), metal('Pt')],
    stack: [ohmic({ I }, ['e-']), layer('Pt', 1e-6), { reactions: [{ ...iron, k0: 1e-20 }, { ...silver, k0: 1e-3 }] },
      layer('water', 20e-6, { c0: { ...bathC, 'Ag+': 1, 'NO3-': 5051 } }), bath(bathC, 'NO3-')],
    grid: { hmin: 10e-9, hmax: 0.5e-6 },
  }));
  const K = Math.exp(-(-78.9e3 + 77.107e3 + 4.7e3) / (VT * FARADAY)); // Fe³⁺ + Ag(s) = Fe²⁺ + Ag⁺
  const d = cell(pulse({ start: 100, width: 1, amplitude: 1 })), s = d.solve(), V0 = cell(0).solve().terminals.left.V;
  assert.ok(s.converged && s.warnings.length === 0, s.warnings.join(' '));
  const g = s.c['Ag+'].findIndex((x) => x > 0), c = (sp) => s.c[sp][g];
  const q = (c('Fe2+') * c('Ag+')) / (c('Fe3+') * 1000) / K;
  assert.ok(Math.abs(q - 1) < 1e-9, `Q/K ${q}`);
  assert.ok(Math.abs(s.terminals.left.V - V0) < 1e-9, `V ${s.terminals.left.V} for ${V0}`);
  // and inside the pulse, the same drive passes current: the laws found again at the solve's time
  d.solver.time = 100.5;
  d.solver._lawsNow();
  assert.equal(d.solver.drivenKey, 'I,V');
  d.solver.time = 200;
  d.solver._lawsNow();
  assert.equal(d.solver.drivenKey, 'I0,V');
});

// A11: whether the state's floating pieces balance is the state's, and a failed solve's restore
// puts it back with the rest. A floating side port's region, warned of, lost its warning after a
// solve that failed.
test('a failed solve after a warned one: the restored state keeps its island warning', () => {
  const mat = 'GaAs', def = build({ T: 250, library: [semiconductor(mat, { T: 250 })], stack: [ohmic(0), layer(mat, 0.483e-6, { name: 'n1', donors: 0.66 }), layer(mat, 0.22e-6, { name: 'p', acceptors: 8.94 }), layer(mat, 0.902e-6, { name: 'n2', donors: 3.5 }), ohmic(-0.2)], grid: { hmin: 0.5e-9, hmax: 40e-9, ratio: 1.11 } });
  def.ports = [{ name: 'side', region: 'n2', from: 0.75e-6, to: 0.9e-6, I: 1e-3, terminal: 'e-', species: { 'e-': 'equilibrium' } }];
  const d = new Device(def), island = (s) => s.warnings.filter((w) => w.includes("floating region doesn't balance")).length;
  assert.equal(island(d.solve()), 1);
  const f = d.solver._factor;
  d.solver._factor = () => { throw new Error('injected'); };
  try {
    assert.equal(d.solve({ maxSteps: 1, continuation: false }).converged, false);
  } finally {
    d.solver._factor = f;
  }
  assert.equal(island(d.solution()), 1);
});

// A13: the last fallback, the pseudo-transient ramp after continuation fails, is a solve like the
// others: a stepped state solved there is no longer stepped. Forced (stage 1 and continuation fail
// untouched), the solution read the right contact's waveform before the step it was at, 0.1 V for 0.
test('the last fallback ramp marks the state solved: no longer stepped, read at its drives', () => {
  const NA = units.perCm3(1e16), d = new Device(build({ T: 300, library: [semiconductor('Si')],
    stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e18) }), layer('Si', units.um(5), { acceptors: NA }), { ...ohmic(0), V: { t: [0, 1e-3, 1e-3], values: [0.1, 0.1, 0] } }],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1 / (1e-6 * NA) } }, photogeneration({ material: 'Si', flux: 1e-4, alpha: 1e5, mu: units.eV(3), to: units.um(5.5) })] }));
  d.solve();
  d.advance(1e-3); // (onto the waveform's step: stepped, at t = 1e-3)
  const sv = d.solver, orig = sv._steadySteps;
  let outerCalls = 0, homCalls = 0;
  sv._steadySteps = function (o) {
    if (sv.allTerminals) { homCalls++; return { converged: false, steps: 1, iterations: 1, history: [] }; }
    if (++outerCalls <= 2) return { converged: false, steps: 1, iterations: 1, history: [] };
    return orig.call(this, o);
  };
  const s = d.solve();
  sv._steadySteps = orig;
  assert.equal(outerCalls, 3);
  assert.ok(homCalls >= 1, `expected homotopy fallback attempts, got ${homCalls}`);
  assert.ok(s.converged);
  assert.equal(sv.stepped, false);
  assert.equal(s.terminals.right.V, 0);
});

// A device solved by giant steps (Y still in part of its stretch: laws.js's unkept), with a gate
// behind a resistance whose RC is 6e5 s: its steady state has the gate at its source.
const stillGate = (R) => ({
  T: 300,
  species: [{ name: 'e-', z: -1 }, { name: 'Y', z: 0 }, { name: 'Z', z: 0 }],
  materials: {
    m1: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 1e-9, mu0: 0, cRef: 1e3 }, Z: { D: 1e-9, mu0: 0, cRef: 1e3 } } },
    m2: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 0, mu0: 0, cRef: 1e3 } } },
  },
  regions: [
    { name: 'A', material: 'm1', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 1, Z: 0.01 } },
    { name: 'B', material: 'm2', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 2 } },
  ],
  interfaces: [{ dipole: 0 }],
  bulkReactions: [{ equation: 'Y = Z', kf: { m1: 1 } }],
  contacts: {
    left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
    right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
  },
  ports: [{ name: 'gate', region: 'A', from: 0.2e-6, to: 0.8e-6, V: 0, R, area: 1, capacitance: { C: 1e-6, zeroCharge: 0 } }],
  grid: { minCells: 40 },
});

// A12: a giant step is accepted as steady only where the terminals' voltages moved under tol over
// it too, not where Newton's last update to them did: the gate, one step 0.26 VT from where it
// began and 3e-15 in Newton's last update, was left at −6.2 mV for its 1 V.
test('giant steps: a gate behind a huge resistance is followed to its source, not left where one step put it', () => {
  const d = new Device(stillGate(1e18));
  assert.equal(d.solver._directSteady(), false);
  const s0 = d.solve();
  assert.ok(s0.converged && Math.abs(s0.terminals.gate.V) < 1e-9, `at 0 V: ${s0.terminals.gate.V}`);
  const s1 = d.set({ ports: { gate: { V: 1 } } }).solve();
  assert.ok(s1.converged && Math.abs(s1.terminals.gate.V - 1) < 1e-9, `at 1 V: ${s1.terminals.gate.V}`);
});

// A14: the same device without the gate, its Y fed at the left contact (1e-5 mol/m² of Y + Z in
// A), then blocked there: Y + Z is closed, unkept by any row, and the giant steps must keep it.
// Where stage 1 and the continuation's first cold solve failed (forced), the continuation laid the
// device out again from c0 (1.03e-6), and each later call took that as its baseline.
test('giant steps keep an unkept amount from when the solve began, through a continuation', () => {
  const RT = GAS_CONSTANT * 300, def = stillGate(1);
  delete def.ports;
  def.contacts.left.species.Y = { type: 'equilibrium', mu: RT * Math.log(5 / 1000) };
  const d = new Device(def);
  d.solve();
  d.set({ contacts: { left: { species: { Y: 'blocked' } } } });
  d.set({ contacts: { right: { V: 0.05 } } });
  const sv = d.solver, totA = () => { let t = 0; for (let g = 0; g < sv.nNodes; g++) if (sv.model.grid.nodeRegion[g] === 0) t += sv.model.grid.vol[g] * (sv.c[g * sv.n + 1] + sv.c[g * sv.n + 2]); return t; };
  assert.equal(sv._directSteady(), false);
  const before = totA(), orig = sv._steadySteps;
  const solveWith = (fake) => {
    let calls = 0;
    sv._steadySteps = function (o) { return fake(++calls, this) ?? orig.call(this, o); };
    try {
      return d.solve();
    } finally {
      sv._steadySteps = orig;
    }
  };
  // stage 1 and the continuation's cold first solve failing: no layout from c0 after them
  let s = solveWith((k) => (k <= 2 ? { converged: false, steps: 1, iterations: 1, history: [] } : undefined));
  assert.ok(s.converged);
  assert.ok(Math.abs(totA() / before - 1) < 1e-6, `Y + Z ${totA()} for ${before}`);
  // and a stage that moved it anyway (laid out from c0 here): the giant steps after it don't take it
  d.set({ contacts: { right: { V: 0.1 } } });
  s = solveWith((k, solver) => {
    if (k === 1) solver.initFromComposition();
    return undefined;
  });
  assert.ok(!s.converged || Math.abs(totA() / before - 1) < 1e-6, `Y + Z ${totA()} for ${before}`);
});

// A25: the held-voltage march of a current-driven iron couple (its floated solve made to fail
// until the march begins). One held solve failing, once, ended that direction: the march went
// the other way to −23 V, and the solve failed claiming a limiting current where 3.3775 V passes
// 0.89 A/m². Nor is a march its held solves stopped short of ±20 V read as one.
test('the floating continuation steps past a failed held sample, and claims no limit it never reached', () => {
  const iron = () => new Device(build({
    library: [aqueous(['Fe3+', 'Fe2+', 'SO42-'], { epsr: 0 }), metal('Pt')],
    stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1.4e-5, alpha: 0.33 }] },
      layer('water', 1.97e-5, { c0: { 'Fe3+': 27.66, 'Fe2+': 0.435, 'SO42-': 41.925 } }),
      bath({ 'Fe3+': 27.66, 'Fe2+': 0.435, 'SO42-': 41.925 }, 'SO42-', { I: 0.89 })],
    grid: { hmin: 1.3e-8, hmax: 4.2e-7, ratio: 1.17 } }));
  const faulty = (fails) => {
    const d = iron(), s = d.solver, from = s._steadyFromHere, steps = s._steadySteps;
    let held = 0;
    s._steadyFromHere = function (opts) {
      if (this.terms[1].drive.kind === 'V' && fails(++held)) return { converged: false, steps: 0, iterations: 0, history: [] };
      return from.call(this, { ...opts, homotopy: false });
    };
    s._steadySteps = function (opts) { // (floated, failing until the march begins)
      return this.terms[1].drive.kind === 'I' && held === 0 ? { converged: false, steps: 1, iterations: 0, history: [] } : steps.call(this, opts);
    };
    return d;
  };
  const s = faulty((n) => n === 4).solve();
  assert.ok(s.converged && Math.abs(s.terminals.right.V - 3.3774722) < 1e-6, `${s.converged} ${s.terminals.right.V}: ${s.warnings.join(' ')}`);
  const d = faulty((n) => n > 3), f = d.solve();
  assert.equal(f.converged, false);
  assert.equal(d.solver.unreached, null);
  assert.ok(!f.warnings.some((w) => w.includes('no steady state passes')), f.warnings.join(' '));
});

// A7: the floated solve failing, the held one at the target's current by regula falsi, accepted
// within 1e-10 of a scale: of every current the march saw, the far end's 7e9 A/m², it took
// −5.3e-7 for a target of +1e-12. (Held solves stubbed by an exponential I(V).)
test('the floating continuation\'s regula falsi meets a small target, not the march\'s largest current', () => {
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), V0 = 0.033333, I = (V) => Math.expm1((V - V0) / 0.025);
  const end = { terminal: 'Cl-', phi: 'neutral', species: { 'Cl-': 'equilibrium', 'K+': 'blocked' } };
  const d = new Device({ species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, ...end }, right: { I: 1e-12, ...end } }, grid: { minCells: 8 } });
  const s = d.solver;
  s.termV[1] = 0.6;
  s._steadyFromHere = function () { const V = this.terms[1].drive.src.value; this.termV[1] = V; this.termI[1] = I(V); return { converged: true, steps: 1, iterations: 1, history: [] }; };
  s._steadySteps = () => ({ converged: false, steps: 1, iterations: 1, history: [] });
  s._bound = () => false;
  const r = s._floatingContinuation({}, { steps: 0, iterations: 0, history: [] }, 1);
  assert.ok(r.converged);
  assert.ok(Math.abs(I(s.termV[1]) - 1e-12) < 1e-13, `I ${I(s.termV[1])} for 1e-12`);
});

// A26: Newton's own reason for failing is kept for the failure's message (newtonError), where a
// singular border ended in "nothing specific showed why". (A solve's message is made in
// src/solution.js; advance's here, in src/timeSteps.js.)
test('Newton\'s own failure reason is kept for the failure, and advance() says it', () => {
  const d = new Device(build({ library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1 }), ohmic(0.01)] }));
  d.solve();
  const f = d.solver._factor, why = 'the bordered system is singular: a conserved amount or a floating terminal is not determined';
  d.solver._factor = () => { throw new Error(why); };
  try {
    assert.equal(d.solve({ maxSteps: 1, continuation: false }).converged, false);
    assert.equal(d.solver.newtonError, why);
    const a = d.advance(1e-6);
    assert.equal(a.converged, false);
    assert.ok(a.warnings.some((w) => w.includes(`Newton's last try: ${why}`)), a.warnings.join(' '));
  } finally {
    d.solver._factor = f;
  }
  assert.ok(d.solve().converged);
  assert.equal(d.solver.newtonError, null);
});

// A29: solve()'s tol and maxSteps checked as advance()'s are (tol NaN with maxSteps Infinity
// stepped without end); a giant step that overflows refused, naming the species (D = 1e-315 made
// it ∞, and the solve converged with 18% of Y + Z lost); and a contact's ramp over a span whose
// eighth is 0 (5e-324 V) ends.
test('solve() refuses a NaN tol, an endless maxSteps and a giant step that overflows; a vanishing ramp ends', () => {
  const d = new Device(stillGate(1));
  for (const opts of [{ tol: NaN }, { tol: Infinity }, { tol: 0 }, { tol: NaN, maxSteps: Infinity }, { maxSteps: Infinity }, { maxSteps: 0 }, { maxSteps: 2.5 }]) {
    assert.throws(() => d.solve(opts), /solve\(\): (tol|maxSteps) is/, JSON.stringify(opts));
  }
  const def = stillGate(1);
  delete def.ports;
  def.species.push({ name: 'W', z: 0 });
  def.materials.m1.species.W = { D: 1e-315, mu0: 0, cRef: 1e3 };
  def.regions[0].c0.W = 1e-3;
  const w = new Device(def), u = Float64Array.from(w.solver.u);
  assert.throws(() => w.solve(), /materials\.m1\.species\['W'\]\.D: 1e-315 m²\/s makes the steady solve's giant time step/);
  assert.deepEqual(w.solver.u, u);
  const film = new Device({
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { w: { epsr: 0, species: { 'K+': { D: 1e-9, cRef: 1, mu0: 0 }, 'Cl-': { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'w', length: 1e-6, c0: { 'K+': 1, 'Cl-': 1 } }],
    contacts: Object.fromEntries(['left', 'right'].map((side) => [side, { V: 0, terminal: 'K+', species: { 'K+': 'equilibrium', 'Cl-': 'blocked' }, phi: 'neutral' }])),
    grid: { minCells: 8 },
  });
  film.solve();
  const s = film.solver, orig = s._steadySteps;
  let calls = 0;
  s._steadySteps = function (o) { assert.ok(++calls < 50, 'the ramp ends'); return orig.call(this, o); };
  s.steady = true;
  try {
    assert.ok(s._continuation({}, 1, 0, 5e-324, { converged: true, steps: 0, iterations: 0, history: [], ramped: true }, true).converged);
  } finally {
    s.steady = false;
    s._steadySteps = orig;
  }
});

// A1: a neutral tracer W fed only by an exchange link at the left contact (its reservoir at
// 2 mol/m³), on the device solved by giant steps. At k = 1e-30 mol/(m²·s), far too slow for them,
// the solve converged with W at its c0 of 1, unwarned; the steady state has it at 2.
test('a contact link too slow for the giant steps, pushing one way, is warned of', () => {
  const make = (k) => {
    const def = stillGate(1);
    delete def.ports;
    delete def.interfaces;
    def.species.push({ name: 'W', z: 0 });
    def.materials.m1.epsr = def.materials.m2.epsr = 0;
    def.materials.m1.species.W = { D: 1e-9, mu0: 0, cRef: 1e3 };
    def.regions[0].c0 = { Y: 1, Z: 1, W: 1 };
    def.contacts.left.species.W = { type: 'exchange', k, mu: GAS_CONSTANT * 300 * Math.log(2 / 1000) };
    return new Device(def);
  };
  const slow = make(1e-30).solve(), link = "contacts.left.species['W']";
  assert.ok(slow.warnings.some((w) => w.startsWith(`${link}: out of equilibrium by 17.9 mV`)), slow.warnings.join(' '));
  const fast = make(1e-4).solve();
  assert.ok(fast.converged && Math.abs(fast.c.W[0] - 2) < 1e-9 && !fast.warnings.some((w) => w.includes(link)), fast.warnings.join(' '));
});

// A10: two closed spectators sharing a lattice's sites, A's level moved 0.2 RT: renormalised in
// one pass in listing order, B (listed first) was left 4.5% off by A's shift. Either order now.
test('spectators on a shared lattice are each renormalised to their amounts, in either listing order', () => {
  for (const order of [['A', 'B'], ['B', 'A']]) {
    const def = stillGate(1);
    delete def.ports;
    def.species.push(...order.map((name) => ({ name, z: 0 })));
    Object.assign(def.materials.m1.species, { A: { D: 1e-9, mu0: 0, cRef: 1 }, B: { D: 1e-9, mu0: 0, cRef: 1 } });
    def.materials.m1.statistics = [{ type: 'lattice', species: ['A', 'B'], cMax: 10 }];
    Object.assign(def.regions[0].c0, { A: 4, B: 4 });
    const sv = new Device(def).solver, M = sv.M;
    const at = (name) => sv.stretches.findIndex((st) => sv.model.species[st.species].name === name);
    const [kA, kB] = [at('A'), at('B')];
    for (const k of [kA, kB]) {
      assert.ok(sv.stretches[k].spectator);
      sv.referenceAmounts[k] = sv.amount(sv.stretches[k]);
    }
    for (let g = sv.stretches[kA].nodes[0]; g <= sv.stretches[kA].nodes[1]; g++) sv.u[sv.blockOfNode[g] * M + 1 + sv.stretches[kA].species] += 0.2;
    sv.computeConcentrations();
    sv._renormalizeSpectators();
    for (const k of [kA, kB]) assert.ok(Math.abs(sv.amount(sv.stretches[k]) / sv.referenceAmounts[k] - 1) < 1e-10, `${order}: ${sv.model.species[sv.stretches[k].species].name} ${sv.amount(sv.stretches[k]) / sv.referenceAmounts[k]}`);
  }
});
