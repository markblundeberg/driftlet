import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous } from '../src/kit.js';
import { sourceAt, nextBreakpoint } from '../src/device.js';

// Faults found in the time stepper and the sources by two outside reading passes (OpenAI's Sol
// on src/timeSteps.js, Astra on the definition's normalization), each confirmed by a probe
// against a second route before it was fixed.

// The blocking cell: KCl between two capacitive contacts, a terminal driven as given.
const blocking = (right, left = { V: 0 }) => {
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2;
  return new Device({
    species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { ...left, phi: { type: 'capacitive', C }, zeroCharge: 0 }, right: { ...right, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
};
const last = (a) => a[a.length - 1];

test("a repeating waveform reads the after side at every edge a step lands on, and lands on each", () => {
  // ((0.5 % 0.1) is 0.0999…98: a square wave of period 0.1 s read 1183 of its first 2000 landed
  // edges from before them, and solve() there the wrong drive)
  const sq = { t: [0, 0, 0.05, 0.05, 0.1], values: [0, 1, 1, 0, 0], repeat: true };
  let t = 0;
  for (let i = 0; i < 2000; i++) {
    t = nextBreakpoint(sq, t);
    const rising = Math.round(t / 0.05) % 2 === 0;
    assert.equal(sourceAt(sq, t), rising ? 1 : 0, `after the edge at ${t}`);
    assert.equal(sourceAt(sq, t, true), rising ? 0 : 1, `before the edge at ${t}`);
  }
  assert.ok(Math.abs(t - 100) < 1e-9);
  // a 1 µs square wave late in a run: every edge, not none (past 1e-12 of t, all were passed)
  const fast = { t: [0, 0, 5e-7, 5e-7, 1e-6], values: [0, 1, 1, 0, 0], repeat: true };
  assert.ok(Math.abs(nextBreakpoint(fast, 2e6) - 2e6 - 5e-7) < 1e-8);
  // a long period: its first step isn't lost to (1 % T) + T rounding to T
  const long = { t: [0, 1, 1, 1e16], values: [0, 0, 1, 1], repeat: true };
  assert.equal(sourceAt(long, 1), 1);
  assert.equal(nextBreakpoint(long, 0), 1);
  // a period below the time's resolution: refused, not a hang
  assert.throws(() => nextBreakpoint({ t: [0, 1], values: [0, 1], repeat: true }, 2 ** 53), /resolution/);
});

test("a first fixed step onto a current's edge posts the before side's charge", () => {
  // (the currents were read before the step was marked stepped: the far side's I·dt, for good)
  const I = { t: [0, 0.01, 0.01], values: [0, 0, 1] };
  for (const method of ['be', 'bdf2']) {
    const d = blocking({ I });
    d.solve();
    d.step(0.01, { method });
    const a = d.advance(0.02);
    assert.ok(Math.abs(last(a.trace.terminals.right.charge) - a.plates.right.charge) < 1e-6 * a.plates.right.charge, `${method}: ledger ${last(a.trace.terminals.right.charge)}, plate ${a.plates.right.charge}`);
  }
});

test("a call ending a short tail after an edge steps it, rather than snapping the clock past it", () => {
  // (within 1e-10 of the time it snapped: 5e-6 s at 1e5 s, the current through it lost for good)
  const T = 1e5, tail = 5e-6, d = blocking({ I: { t: [0, T, T], values: [0, 0, 1] } });
  d.solve();
  const a = d.advance(T + tail, { dt0: T });
  assert.ok(Math.abs(a.plates.right.charge / tail - 1) < 1e-3, `${a.plates.right.charge} for ${tail}`);
  assert.equal(d.solver.sourceTime, d.solver.time);
});

test("a step's halves restart the order on an edge they start on, at any depth", () => {
  // (only the outermost call did: a half starting on the edge crossed it by BDF2, a third short)
  const d = blocking({ I: { t: [0, 0.01, 0.01], values: [0, 0, 1] } });
  d.solve();
  d.step(0.005, { method: 'bdf2' });
  const s = d.solver, real = s.newton;
  let calls = 0;
  s.newton = function (dt, opts) {
    return ++calls === 1 ? { converged: false, iterations: 0, history: [] } : real.call(this, dt, opts);
  };
  const r = d.step(0.01, { method: 'bdf2' });
  s.newton = real;
  assert.ok(Math.abs(r.plates.right.charge / 0.005 - 1) < 1e-3, `${r.plates.right.charge} for 0.005`);
});

test('advance() and step() refuse what they can\'t mean, before touching the state', () => {
  const d = blocking({ V: 0 });
  d.solve();
  for (const opts of [{ dtMax: -0.1 }, { dtMax: 0 }, { tol: NaN }, { tol: -1 }, { maxSteps: 1.5 }]) assert.throws(() => d.advance(1, opts), /advance\(\)/, JSON.stringify(opts));
  assert.throws(() => d.advance('1'), /tEnd/);
  d.advance(1);
  const before = d.solver.time, sol = d.step(1e-17);
  assert.equal(d.solver.time, before);
  assert.ok(sol.warnings.some((w) => /resolution/.test(w)));
});

test("set() keeps its own copies of a patch's drive and links, and the constructor of the definition", () => {
  // (the caller's objects were put in: an edit to one made the next set() of it a no-op, and
  // leaked in at an unrelated rebuild; the constructor's model kept the caller's arrays)
  const d = blocking({ V: 0 });
  d.solve();
  const w = { t: [0, 1], values: [0, 1] };
  d.set({ contacts: { right: { V: w } } });
  w.values[1] = 2;
  assert.deepEqual(d.def.contacts.right.V.values, [0, 1]);
  d.set({ contacts: { right: { V: w } } });
  assert.deepEqual([...d.model.contacts.right.drive.src.values], [0, 2]);
  const def = { species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } }, statistics: [{ type: 'redlich-kister', species: ['X'], cMax: 3, A: [100] }] } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
    grid: { minCells: 8 } };
  const e = new Device(def), mu = e.solution().mu.X[0];
  def.materials.f.statistics[0].A[0] = 1000;
  assert.equal(new Device(e.def).solution().mu.X[0], mu);
  assert.equal(e.solution().mu.X[0], mu);
});

// The cold start (src/start.js): each a start that sat where its mirror or order-permuted twin,
// or the analytic value, didn't, the first step or solve then getting there from it.

test('a bath without a reference reads its cold voltage by φ, whatever order def.species lists its ions in', () => {
  // (its terminal was the first listed: a salt ramped 1 → 9 against baths at 1 read the floating
  // end at +56 or −56 mV by the order, the first step at 0 in both)
  const L = 1e-6, run = (order) => {
    const d = new Device({
      species: order.map((nm) => ({ name: nm, z: nm === 'P' ? 1 : -1 })),
      materials: { w: { epsr: 0, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
      regions: [{ name: 'a', material: 'w', length: L, c0: { P: { x: [0, L], values: [1, 9] }, N: { x: [0, L], values: [1, 9] } } }],
      contacts: { left: { V: 0, bath: { c: { P: 1, N: 1 } } }, right: { I: 0, bath: { c: { P: 1, N: 1 } } } },
      grid: { minCells: 16 },
    });
    return [d.solution().contacts.right.V, d.advance(1e-9).contacts.right.V];
  };
  for (const [cold, first] of [run(['P', 'N']), run(['N', 'P'])]) assert.ok(Math.abs(cold) < 1e-9 && Math.abs(first) < 1e-9, `cold ${cold} V, after 1 ns ${first} V`);
});

test('a bath without a reference outside a cation exchanger is accepted in either order of def.species', () => {
  // (its terminal, the first listed, was N⁻, which the membrane doesn't hold: refused as blocked)
  for (const order of [['P', 'N'], ['N', 'P']]) {
    const d = new Device({
      species: order.map((nm) => ({ name: nm, z: nm === 'P' ? 1 : -1 })),
      materials: { water: { epsr: 0, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } }, mem: { epsr: 0, species: { P: { D: 1e-10, cRef: 1, mu0: 0 } } } },
      regions: [{ name: 'm', material: 'mem', length: 1e-6, fixedCharge: -100 * FARADAY, c0: { P: 100 } }],
      contacts: { left: { V: 0, bath: { c: { P: 1, N: 1 }, material: 'water' } }, right: { V: 0, bath: { c: { P: 1, N: 1 }, material: 'water' } } },
      grid: { minCells: 16 },
    });
    const s = d.solve(), donnan = -d.solver.VT * Math.log(100);
    assert.ok(s.converged && Math.abs(s.phi[0] - donnan) < 1e-9, `${order}: φ ${s.phi[0]} V for ${donnan}`);
  }
});

test('a floating capacitance through a window beside a pinned end starts uncharged', () => {
  // (its starting V averaged φ over the end box too, whose φ the contact sets and whose charge
  // the port doesn't hold: 0.0087 C/m² on two cells, falling with the end box)
  const d = new Device({
    species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }],
    materials: { f: { epsr: 80, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1, c0: { N: { x: [0, 1], values: [1, 9] } } }],
    contacts: { left: { V: 0, terminal: 'P', phi: 'pinned', zeroCharge: 0, species: { P: 'equilibrium', N: 'blocked' } }, right: { phi: 'neutral', species: { P: 'blocked', N: 'blocked' } } },
    ports: [{ name: 'cap', region: 'a', I: 0, area: 1, capacitance: { C: 1 } }],
    grid: { minCells: 2, hmin: 0.5, hmax: 0.5 },
  });
  const Q = d.solution().ports[0].charge;
  assert.ok(Math.abs(Q) < 1e-12, `${Q} C/m²`);
});

test('an electrode at the far end of a ramped salt starts balanced with the salt beside it, as its mirror does', () => {
  // (the balance took the salt at the region's start: a metal on the right of Ag⁺ ramped 10 → 90
  // started 2.2 RT (ln 9) out of equilibrium, the mirror at it; a port's window likewise)
  const T = 298.15, RT = GAS_CONSTANT * T;
  const species = [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }];
  const water = { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } };
  const L = 1e-5, salt = (a, b) => ({ 'Ag+': { x: [0, L], values: [a, b] }, 'NO3-': { x: [0, L], values: [a, b] } });
  const closed = { phi: 'neutral', species: { 'Ag+': 'blocked', 'NO3-': 'blocked', 'e-': 'blocked' } };
  const right = new Device({
    T, species, materials: { water, Ag: { conductor: { species: 'e-', conductivity: 6e7 } } },
    regions: [{ material: 'water', length: L, c0: salt(10, 90) }, { material: 'Ag', length: 1e-6 }],
    interfaces: [{ phi: 'neutral', reactions: [{ right: { 'e-': -1, Ag: 1 }, left: { 'Ag+': -1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 }] }],
    contacts: { left: closed, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } }, grid: { minCells: 8 },
  });
  const g = right.model.grid, s = right.solution(), a = (s.mu['Ag+'][g.regionEnd[0]] + s.mu['e-'][g.regionStart[1]]) / RT;
  assert.ok(Math.abs(a) < 1e-9, `affinity ${a} RT at the metal on the right`);
  // a window at either end: the cold rates the mirror of each other's
  const port = (mirror) => {
    const d = new Device({
      T, species, materials: { water }, regions: [{ name: 'w', material: 'water', length: L, c0: mirror ? salt(90, 10) : salt(10, 90) }],
      contacts: { left: closed, right: closed },
      ports: [{ name: 'Ag', region: 'w', from: mirror ? 0 : 0.9 * L, to: mirror ? 0.1 * L : L, V: 0, terminal: 'e-', area: 1e4, reactions: [{ equation: 'Ag+ + e- = Ag', fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 }] }],
      grid: { minCells: 20, hmin: L / 20, hmax: L / 20 },
    });
    d.solver.assemble(Infinity);
    return Array.from(d.solver.portRates[0][0]);
  };
  const p = port(false), q = port(true).reverse(), scale = Math.max(...p.map(Math.abs));
  p.forEach((r, k) => assert.ok(Math.abs(r - q[k]) < 1e-9 * scale, `rates ${p} against the mirror's ${q}`));
});

test("a floating electrode with a conductance beside its reactions starts where the whole terminal passes what it's set", () => {
  // (the start balanced its reactions alone: 20 A/m² through the link at I: 0; and its search,
  // unchecked, ran out to 1457 V for a drive beyond the reactions' plateau)
  const def = (rx, I, link = true) => ({
    species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }, { name: 'e-', z: -1 }],
    materials: { w: { epsr: 0, species: { P: { D: 1e-9, cRef: 1, mu0: 0 }, N: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'w', length: 1, c0: { P: 1, N: 1 } }],
    contacts: { left: { V: 0, terminal: 'P', phi: 'neutral', species: { P: 'equilibrium', N: 'blocked' } }, right: { phi: 'neutral', species: { P: 'blocked', N: 'blocked' } } },
    ports: [{ name: 'el', region: 'a', I, terminal: 'e-', area: 1, reactions: [rx], ...(link ? { species: { N: { type: 'conductance', G: 100, offset: 0.2 } } } : {}) }],
    grid: { minCells: 4 },
  });
  const bv = { equation: 'P + e- = S', fixed: { S: 0 }, k0: 1e-6, alpha: 0.5 }, mhc = { equation: 'P + e- = S', fixed: { S: 0 }, k0: 1e-6, lambda: 20000 };
  for (const [rx, I] of [[bv, 0], [mhc, 0], [mhc, 5]]) {
    const d = new Device(def(rx, I)), S = d.solver, V0 = S.termV[2];
    S.assemble(Infinity);
    assert.ok(Math.abs(S.termI[2] - I) < 1e-6, `${rx.lambda ? 'MHC' : 'BV'} at ${I}: cold ${S.termI[2]} A/m² at ${V0} V`);
    if (I === 0) assert.ok(Math.abs(new Device(def(rx, I)).advance(1e-3).ports[0].V - V0) < 1e-6, 'and the first ms stays there');
  }
  // beyond what the reactions alone can pass: left where it was, not at the search's bound
  assert.ok(Math.abs(new Device(def(mhc, 5, false)).solver.termV[2]) < 1, 'MHC beyond its plateau');
});

test("a species held by a port behind a blocked face starts at the port's level, not the right contact's", () => {
  // (its stretch doesn't reach the right contact, which supplied its level all the same: 2.63
  // for the port's 1, kept for ~c/k of the exchange)
  const RT = GAS_CONSTANT * 298.15;
  const d = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6 }, { name: 'b', material: 'f', length: 1e-6 }],
    interfaces: [{ species: { X: 'blocked' } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 2400 } } } },
    ports: [{ name: 'feed', region: 'a', V: 0, species: { X: { type: 'exchange', k: 1e-6, mu: 0 } } }],
    grid: { minCells: 8 },
  });
  const s = d.solution(), last = s.x.length - 1;
  assert.ok(Math.abs(s.c.X[0] - 1) < 1e-12, `cold X in a: ${s.c.X[0]}`);
  assert.ok(Math.abs(s.c.X[last] / Math.exp(2400 / RT) - 1) < 1e-12, `and in b, the contact's: ${s.c.X[last]}`);
  assert.ok(Math.abs(d.advance(1e-3).c.X[0] - 1) < 1e-9);
});

test("an insertion host's wire starts at the open circuit whatever order lists its ion and electrons", () => {
  // (the host's levels were aligned on its first charged species shared across the face: e⁻,
  // listed first and blocked there, left its Li⁺ 0.4 V out of equilibrium and the wire at 0 V)
  const VT = (GAS_CONSTANT * 298.15) / FARADAY;
  const x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95], E = x.map((v) => 0.4 - VT * Math.log(v / (1 - v)));
  const sp = { 'Li+': { name: 'Li+', z: 1 }, 'e-': { name: 'e-', z: -1 } };
  for (const order of [['Li+', 'e-'], ['e-', 'Li+']]) {
    const d = new Device({
      species: order.map((k) => sp[k]),
      materials: {
        mc: { epsr: 0, species: { 'Li+': { D: 1e-9, mu0: 0, cRef: 1000 }, 'e-': { D: 1e-6, mu0: 0, cRef: 1000 } } },
        host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }] },
      },
      regions: [{ name: 'm', material: 'mc', length: 1e-6, c0: { 'Li+': 1000, 'e-': 1000 } }, { name: 'h', material: 'host', length: 1e-6, c0: { 'Li+': 15000 } }],
      interfaces: [{ species: { 'Li+': 'equilibrium', 'e-': 'blocked' } }],
      contacts: {
        left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium', 'e-': { type: 'equilibrium', offset: 0.5 } }, phi: 'bulk' },
        right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      },
      grid: { hmin: 5e-8, hmax: 1e-7 },
    });
    const V = d.solution().terminals.right.V;
    assert.ok(Math.abs(V - 0.4) < 1e-6, `[${order}]: cold ${V} V`);
  }
});

test('a strictly neutral region under a held capacitance starts with ρ + aσ = 0, as the solve has it', () => {
  // (the start left aσ out: a region of P⁺ alone under a gate was refused for the charge the gate
  // balances, though the solve at any εr > 0 gives φ = F c/(a C); and a level-fed one started
  // off it, its first 1 ps step jumping the gate's charge from 0.28 to 0.0035 C/m²)
  const lone = new Device({
    species: [{ name: 'P+', z: 1 }],
    materials: { f: { epsr: 0, species: { 'P+': { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { 'P+': 1 } }],
    contacts: { left: { phi: 'neutral', species: { 'P+': 'blocked' } }, right: { phi: 'neutral', species: { 'P+': 'blocked' } } },
    ports: [{ name: 'g', region: 'a', V: 0, area: 1e7, capacitance: { C: 1 } }],
    grid: { minCells: 8 },
  });
  for (const s of [lone.solution(), lone.solve()]) for (const phi of s.phi) assert.ok(Math.abs(phi - FARADAY / 1e7) < 1e-12, `φ ${phi} V for F/(aC)`);
  const bar = new Device({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * FARADAY }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'dl', region: 'bar', V: 0.2, area: 1e6, capacitance: { C: 0.2 } }],
    grid: { minCells: 40 },
  });
  const Q0 = bar.solution().ports[0].charge, Q1 = bar.step(1e-12).ports[0].charge, Q = bar.solve().ports[0].charge;
  assert.ok(Math.abs(Q1 / Q0 - 1) < 1e-9 && Math.abs(Q / Q0 - 1) < 1e-9, `cold ${Q0}, after 1 ps ${Q1}, steady ${Q} C/m²`);
});
