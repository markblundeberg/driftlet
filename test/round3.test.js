import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { aqueous, bath } from '../src/kit.js';
import { sourceAt, nextBreakpoint } from '../src/device.js';
import { electrodeUnbalanced } from '../src/equation.js';

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

// Definitions accepted that should have been refused, and refused that should have been
// accepted: Astra's third pass over the definition's normalization.

test('a contact at an end of no cross-section (a sphere centre, a cylinder axis) holds no level or φ there', () => {
  // (an equilibrium link at a sphere's centre passed 1.08 mol/s on 4 cells, 0.016 on 256, where
  // nothing passes a point; a pinned φ there stored 0.4 of a capacitor of radius h/2)
  const def = (geometry, left) => ({ geometry, species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1, c0: { X: 1 } }],
    contacts: { left, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } } },
    grid: { minCells: 4, hmin: 0.25, hmax: 0.25 } });
  const held = { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 2400 } } };
  for (const type of ['spherical', 'cylindrical']) {
    assert.throws(() => new Device(def({ type, r0: 0 }, held)), /contacts\.left: the cross-section is zero at this end.*holds 'X's level/);
    assert.throws(() => new Device(def({ type, r0: 0 }, { phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } })), /holds φ there \('pinned'\)/);
  }
  // At a radius, it holds; at the centre a rate link passes what the end's zero area makes of it, nothing.
  assert.ok(new Device(def({ type: 'spherical', r0: 0.1 }, held)).solve().converged);
  const ex = new Device(def({ type: 'spherical', r0: 0 }, { phi: 'neutral', species: { X: { type: 'exchange', k: 1, mu: 2400 } } })).solve();
  assert.equal(ex.contacts.left.flux.X, 0);
});

test("a velocity link on a conductor's carrier is refused, as a permeability is", () => {
  // (a conductor's c is its carrier's excess, 0 in the bulk, so v·c·(…) passed nothing at any v,
  // 1e30 m/s included: current 0 where an equilibrium link passes σV/L = 1 A/m²)
  const def = (link) => ({ species: [{ name: 'e-', z: -1 }], materials: { Cu: { conductor: { species: 'e-', conductivity: 1 } } },
    regions: [{ material: 'Cu', length: 1 }],
    contacts: { left: { V: 1, terminal: 'e-', species: { 'e-': link }, phi: 'neutral' }, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' } } });
  assert.throws(() => new Device(def({ type: 'velocity', v: 1e6 })), (e) => e instanceof DeviceError && /contacts\.left\.species\.e-: a velocity link .* 'Cu' is a conductor/.test(e.message));
  assert.ok(Math.abs(Math.abs(new Device(def('equilibrium')).solve().current) - 1) < 1e-12);
});

test("a bath on a conductor is refused at a port, and called a conductor's at a contact", () => {
  // (a port's bath read the conductor's placeholder μ° = 0 and cRef = 1: e⁻ at 100 put the wire
  // VT·ln 100 = 0.118 V off the port's V, 0.118 A/m² passing between two grounded terminals)
  const def = (port, left) => ({ species: [{ name: 'e-', z: -1 }], materials: { Au: { conductor: { species: 'e-', conductivity: 1 } } },
    regions: [{ name: 'w', material: 'Au', length: 1 }],
    contacts: { left: left ?? { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [port] });
  const port = { name: 'p', region: 'w', V: 0, species: { 'e-': 'equilibrium' } };
  assert.throws(() => new Device(def({ ...port, bath: { c: { 'e-': 100 } } })), /ports\[0\]\.bath: w \(Au\) is a conductor/);
  assert.throws(() => new Device(def({ ...port, terminal: 'e-' }, { bath: { c: { 'e-': 100 } } })), (e) => /contacts\.left\.bath: Au at this end is a conductor/.test(e.message) && !/insertion host/.test(e.message));
  assert.equal(new Device(def({ ...port, terminal: 'e-' })).solve().ports[0].current, 0);
});

test("a reacting port's electrode side is its carrier and fixed participants, whatever else is declared", () => {
  // (a Y in fixed that another material declared went on the electrode's side as the carrier:
  // 'P+ = Y' passed 2.2 mA/m² with no e⁻ and its 9000 J/mol dropped; 'P+ + e- = Y' failed on a NaN)
  const def = (declared, equation) => ({
    species: [{ name: 'P+', z: 1 }, { name: 'N-', z: -1 }, { name: 'e-', z: -1 }, ...(declared ? [{ name: 'Y', z: declared.z }] : [])],
    materials: {
      w: { epsr: 0, species: { 'P+': { D: 1e-9, cRef: 1, mu0: 0 }, 'N-': { D: 1e-9, cRef: 1, mu0: 0 } } },
      ...(declared ? { o: { epsr: 0, species: { Y: { D: 1e-9, cRef: 1, mu0: 0 } } } } : {}),
    },
    regions: [{ name: 'a', material: 'w', length: 1e-5, c0: { 'P+': 1, 'N-': 1 } }],
    contacts: { left: { V: 0, terminal: 'P+', species: { 'P+': 'equilibrium', 'N-': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' }, right: { phi: 'neutral' } },
    ports: [{ name: 'el', region: 'a', V: 0, terminal: 'e-', area: 1, reactions: [{ equation, fixed: { Y: 9000 }, k0: 1e-3, alpha: 0.5 }] }],
    grid: { minCells: 8 },
  });
  assert.throws(() => new Device(def({ z: 1 }, 'P+ = Y')), /ports\[0\]\.reactions\[0\]: no 'e-' from the electrode takes part/);
  const twin = new Device(def(null, 'P+ + e- = Y')), d = new Device(def({ z: 0 }, 'P+ + e- = Y'));
  assert.equal(d.model.ports[0].reactions[0].fixedA, twin.model.ports[0].reactions[0].fixedA);
  assert.ok(Math.abs(d.model.ports[0].reactions[0].fixedA + 9000 / d.model.RT) < 1e-12);
  assert.deepEqual(d.model.ports[0].reactions[0].part, twin.model.ports[0].reactions[0].part);
  const s = d.solve();
  assert.ok(s.converged && s.ports[0].current !== 0);
  assert.equal(s.ports[0].current, twin.solve().ports[0].current);
});

test("Marcus–Hush–Chidsey's one electron counts charge: a divalent carrier with ν = 1 is refused", () => {
  // (B²⁻ with ν = −1 passed as one electron, at a face and at a port: two charges through one
  // electron's Fermi-level integral)
  const lambda = units.eV(0.5);
  const def = (zc) => {
    const car = zc === -2 ? 'B2-' : 'e-', ox = zc === -2 ? 'O2+' : 'O+';
    return {
      species: [{ name: ox, z: -zc, cRef: 1000 }, { name: 'R', z: 0, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: car, z: zc }],
      materials: {
        water: { epsr: 0, species: { [ox]: { D: 1e-9, mu0: 10e3 }, R: { D: 0.8e-9, mu0: -86.485e3 }, 'K+': { D: 1.96e-9, mu0: -283.3e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } },
        Pt: { conductor: { species: car, conductivity: 9.4e6 } },
      },
      regions: [{ material: 'Pt', length: 1e-6 }, { material: 'water', length: 20e-6 }],
      interfaces: [{ phi: 'neutral', reactions: [{ left: { [car]: -1 }, right: { [ox]: -1, R: 1 }, k0: 1e-5, lambda }] }],
      contacts: { left: { V: 0, terminal: car, species: { [car]: 'equilibrium' }, phi: 'bulk' }, right: { bath: { c: { [ox]: 1, R: 1, 'K+': 1000, 'Cl-': 1000 - zc }, reference: 'Cl-' } } },
      grid: { minCells: 8 },
    };
  };
  const onPort = (d, zc) => {
    const car = zc === -2 ? 'B2-' : 'e-', ox = zc === -2 ? 'O2+' : 'O+';
    d.interfaces = [{ phi: 'neutral' }];
    d.ports = [{ name: 'wire', region: 1, from: 0, to: 1e-6, V: 0, terminal: car, area: 1e4, reactions: [{ equation: `${ox} + ${car} = R`, k0: 1e-5, lambda }] }];
    return d;
  };
  assert.ok(new Device(def(-1)).model.interfaces[0].reactions[0].mhc);
  assert.ok(new Device(onPort(def(-1), -1)).model.ports[0].reactions[0].mhc);
  assert.throws(() => new Device(def(-2)), /interfaces\[0\].*Marcus–Hush–Chidsey kinetics is one electron to or from a metal, its carrier with \|ν z\| = 1/);
  assert.throws(() => new Device(onPort(def(-2), -2)), /ports\[0\]\.reactions\[0\]: Marcus–Hush–Chidsey kinetics is one electron to or from the electrode, its carrier with \|ν z\| = 1/);
});

test('a port bath may list one of Fermi–Dirac bands declared together: they are each on its own', () => {
  // (refused as coupled, "each one's level depends on the others'", and accepted declared apart)
  const def = (statistics) => ({ species: [{ name: 'A', z: 0 }, { name: 'B', z: 0 }],
    materials: { f: { epsr: 4, species: { A: { D: 1e-9, cRef: 1, mu0: 0 }, B: { D: 1e-9, cRef: 1, mu0: 0 } }, statistics } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { A: 1, B: 1 } }],
    contacts: { left: { phi: 'neutral', species: { A: 'blocked', B: 'blocked' } }, right: { phi: 'neutral', species: { A: 'blocked', B: 'blocked' } } },
    ports: [{ name: 'p', region: 'a', from: 0, to: 0.5e-6, bath: { c: { A: 3 } }, species: { A: { type: 'exchange', k: 1 } } }],
    grid: { minCells: 8 } });
  const together = new Device(def([{ type: 'fermi-dirac', species: ['A', 'B'] }]));
  const apart = new Device(def([{ type: 'fermi-dirac', species: ['A'] }, { type: 'fermi-dirac', species: ['B'] }]));
  assert.equal(together.model.ports[0].species[0].mu, apart.model.ports[0].species[0].mu);
  // (5251.40 J/mol: ζ = 2.11839 for 𝓕½(ζ) = 3, by an independent quadrature, to its 3e-3 J/mol)
  assert.ok(Math.abs(together.model.ports[0].species[0].mu - 5251.403) < 0.01);
  assert.ok(Math.abs(together.solve().c.A[2] - 3) < 1e-6);
  // A lattice still couples them: listing one is refused.
  assert.throws(() => new Device(def([{ type: 'lattice', species: ['A', 'B'], cMax: 10 }])), /share lattice statistics .* list them all/);
});

test('a species belongs to one statistics model, ideal among them, whatever the order', () => {
  // ([ideal, lattice] was accepted as the lattice, and [lattice, ideal] refused)
  const def = (statistics) => ({ species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } }, statistics } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
    grid: { minCells: 8 } });
  const ideal = { type: 'ideal', species: ['X'] }, lattice = { type: 'lattice', species: ['X'], cMax: 2 };
  for (const st of [[ideal, lattice], [lattice, ideal], [ideal, ideal]]) {
    assert.throws(() => new Device(def(st)), /statistics\[1\]\.species\[0\]: 'X' already belongs to another statistics model/, JSON.stringify(st.map((m) => m.type)));
  }
  assert.equal(new Device(def([ideal])).solution().mu.X[0], 0);
  assert.ok(Math.abs(new Device(def([lattice])).solution().mu.X[0] - GAS_CONSTANT * 298.15 * Math.log(2)) < 1e-9);
});

test("an OCV table's insertion species take no μ° grade: the table sets their level", () => {
  // (graded to a constant 10 kJ/mol, the table read 3.86766 V where it says 3.97130, with the
  // material's μ° at 10 kJ/mol as at 0: Δμ°/F = 104 mV off a staircase of such layers)
  const x = [0.1, 0.5, 0.9], E = [4.2, 3.9, 3.6], L = 1e-6;
  const def = (mu0, grade, ocv = { x, E, muRef: 0 }) => ({ species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: { h: { epsr: 0, species: { 'Li+': { D: 1e-14, cRef: 1, mu0 }, 'e-': { D: 1e-8, cRef: 1, mu0: 0 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 1000, ...(ocv ? { ocv } : { A: [] }) }] } },
    regions: [{ name: 'a', material: 'h', length: L, c0: { 'Li+': 400, 'e-': 400 }, ...(grade ? { grade } : {}) }],
    contacts: { left: { V: 0, terminal: 'e-', phi: 'neutral', species: { 'Li+': 'blocked', 'e-': 'equilibrium' } }, right: { phi: 'neutral', species: { 'Li+': 'blocked', 'e-': 'blocked' } } },
    grid: { minCells: 8 } });
  const comb = (d) => { const s = new Device(d).solve(); return s.mu['Li+'][3] + s.mu['e-'][3]; };
  assert.ok(Math.abs(comb(def(0)) - comb(def(10e3))) < 1e-6);
  for (const sp of ['Li+', 'e-']) {
    assert.throws(() => new Device(def(0, { [sp]: { mu0: { x: [0, L], values: [10e3, 10e3] } } })), new RegExp(`regions\\[0\\]\\.grade\\.${sp.replace('+', '\\+')}\\.mu0: .* OCV table`));
  }
  // Its D may be graded, and a Redlich–Kister host's μ° (no table: μ° is its level's).
  assert.ok(new Device(def(0, { 'Li+': { D: { x: [0, L], values: [1e-14, 2e-14] } } })).solve().converged);
  assert.ok(new Device(def(0, { 'Li+': { mu0: { x: [0, L], values: [10e3, 10e3] } } }, null)).solve().converged);
});

test("a port's bath reads a graded μ° in its window, refused only where it varies there", () => {
  // (refused for any grade object at all: one knot, a constant, one beyond the window, and of a
  // species whose link gave its own level)
  const VT = (GAS_CONSTANT * 298.15) / FARADAY, sp = { D: 1e-9, cRef: 1, mu0: 0 };
  const def = (grade, link = 'equilibrium', statistics) => ({ species: [{ name: 'P', z: 1 }, { name: 'N', z: -1 }],
    materials: { m: { epsr: 0, species: { P: sp, N: sp }, ...(statistics ? { statistics } : {}) } },
    regions: [{ name: 'a', material: 'm', length: 1, c0: { P: 1, N: 1 }, ...(grade ? { grade } : {}) }],
    contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'p', region: 'a', from: 0.1, to: 0.2, V: 0, species: { P: link }, bath: { c: { P: 2, N: 2 } } }],
    grid: { minCells: 10 } });
  const offset = (d) => new Device(d).model.ports[0].species[0].offset;
  assert.ok(Math.abs(offset(def(null)) - VT * Math.log(2)) < 1e-15);
  assert.equal(offset(def({ P: { mu0: { x: [0.5], values: [0] } } })), offset(def(null)));
  assert.equal(offset(def({ P: { mu0: { x: [0, 0.8, 1], values: [0, 0, 1000] } } })), offset(def(null)));
  // a constant grade is the level's μ° there, not the material's
  assert.ok(Math.abs(offset(def({ P: { mu0: { x: [0, 1], values: [1000, 1000] } } })) - (VT * Math.log(2) + 1000 / FARADAY)) < 1e-15);
  // varying across the window, it is refused for the level the bath would give, not for a link's own
  const varying = { P: { mu0: { x: [0, 1], values: [0, 1000] } } };
  assert.throws(() => new Device(def(varying)), /ports\[0\]\.species\.P: a grades 'P's mu0 or cRef, so a composition outside gives no one level across the window/);
  assert.equal(offset(def(varying, { type: 'equilibrium', offset: 0.01 })), 0.01);
  // under a lattice, listing both (as it must) with P's level by its link
  const lattice = [{ type: 'lattice', species: ['P', 'N'], cMax: 10 }];
  assert.equal(offset(def(varying, { type: 'equilibrium', offset: 0.01 }, lattice)), 0.01);
});

test("a bath's reference electrode counts each participant once, its ν summed over repeated terms", () => {
  // (M⁺ + e⁻ + e⁻ = M(s) passed, read by its first e⁻; M⁺ + M⁺ + 2 e⁻ = M₂(s) was refused, and
  // the kit's bath() refused both, as needing electrons on one side and having 'M+' and 'M+')
  const species = [{ name: 'M+', z: 1 }, { name: 'A-', z: -1 }];
  const flagged = (eq) => electrodeUnbalanced({ electrode: eq, reference: 'M+' }, species);
  assert.equal(flagged('M+ + M+ + 2 e- = M2(s)'), null);
  assert.equal(flagged('2 M+ + 2 e- = M2(s)'), null);
  assert.match(flagged('M+ + e- + e- = M(s)'), /isn't charge-balanced \(1 M\+ of charge 1 to 2 e-\)/);
  const c = { 'M+': 1, 'A-': 1 }, fixed = { 'M2(s)': -1000 };
  const repeated = bath(c, { equation: 'M+ + M+ + 2 e- = M2(s)', fixed }), once = bath(c, { equation: '2 M+ + 2 e- = M2(s)', fixed });
  assert.equal(repeated.bath.offset, once.bath.offset);
  assert.ok(Math.abs(once.bath.offset + 1000 / (2 * FARADAY)) < 1e-15); // μ_M2/(2F)
  assert.deepEqual({ ...repeated.bath, electrode: null }, { ...once.bath, electrode: null });
  assert.ok(bath(c, { equation: 'M+ + e- + e- = M(s)', fixed: { 'M(s)': 0 } }).bath.electrode);
});

test('SRH and Auger take no fixed participants: they consume one carrier of each sign and make nothing', () => {
  // ('e- + h+ = S' with S fixed at 5 kJ/mol was accepted, in the bulk and at a face, its n·p at
  // e^(5000/RT) = 7.52 n_i²; mass action takes it, as before)
  const semi = { s: { epsr: 0, species: { 'e-': { D: 1e-3, cRef: 1, mu0: 0 }, 'h+': { D: 1e-3, cRef: 1, mu0: 0 } } } };
  const region = (name) => ({ name, material: 's', length: 1e-6, fixedCharge: 100 * FARADAY, c0: { 'e-': 100.01, 'h+': 0.01 } });
  const contacts = { left: { V: 0, terminal: 'e-', phi: 'bulk', species: { 'e-': 'equilibrium', 'h+': 'blocked' } }, right: { phi: 'neutral', species: { 'e-': 'blocked', 'h+': 'blocked' } } };
  const bulk = (rx) => ({ species: [{ name: 'e-', z: -1 }, { name: 'h+', z: 1 }], materials: semi, regions: [region('a')], bulkReactions: [rx], contacts, grid: { minCells: 8 } });
  const face = (rx) => ({ ...bulk(rx), bulkReactions: undefined, regions: [region('a'), region('b')], interfaces: [{ phi: 'neutral', species: { 'e-': 'equilibrium', 'h+': 'blocked' }, reactions: [rx] }] });
  const srh = { s: { tauN: 1e-6, tauP: 1e-6 } }, auger = { s: { Cn: 1e-3, Cp: 1e-3 } };
  for (const [rx, message] of [
    [{ equation: 'e- + h+ = S', fixed: { S: 5000 }, srh }, /bulkReactions\[0\]\.srh\.s: SRH kinetics .*no fixed participants/],
    [{ equation: 'e- + h+ + S = 0', fixed: { S: 0 }, srh }, /bulkReactions\[0\]\.srh\.s: SRH kinetics/],
    [{ equation: 'e- + h+ = S', fixed: { S: 5000 }, auger }, /bulkReactions\[0\]\.auger\.s: Auger kinetics .*no fixed participants/],
  ]) assert.throws(() => new Device(bulk(rx)), message);
  assert.throws(() => new Device(face({ equation: 'e-(left) + h+(right) = S', fixed: { S: 5000 }, srh: { vn: 1, vp: 1 } })), /interfaces\[0\].*reactions\[0\]\.srh: SRH kinetics .*no fixed participants/);
  assert.ok(new Device(face({ equation: 'e-(left) + h+(right) = 0', srh: { vn: 1, vp: 1 } })).solve().converged);
  for (const rx of [{ equation: 'e- + h+ = 0', srh }, { equation: 'e- + h+ = 0', auger }]) assert.ok(new Device(bulk(rx)).solve().converged);
  const s = new Device(bulk({ equation: 'e- + h+ = S', fixed: { S: 5000 }, kf: { s: 1 } })).solve();
  assert.ok(Math.abs((s.c['e-'][4] * s.c['h+'][4]) / Math.exp(5000 / (GAS_CONSTANT * 298.15)) - 1) < 1e-6);
});
