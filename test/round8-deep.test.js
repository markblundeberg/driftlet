// The steady solve (src/steady.js, src/newton.js) against what Astra's re-read of it pointed at
// and a chaser confirmed (round 8, 20 Nov 2026): rollbacks that kept part of what they undid,
// converged returns that skipped a check, a march that claimed a limit it never showed, and two
// absurd inputs that ran without end. Each test's comment gives the old numbers.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Device, FARADAY, units } from '../src/index.js';
import { aqueous, bath, build, layer, ohmic, semiconductor } from '../src/kit.js';
import { newton } from '../src/newton.js';
import { SolverError } from '../src/errors.js';
import { STATE } from '../src/solver.js';

// A16: a closed Fe²⁺/Fe³⁺ film on a finite-step device (a neutral Y = Z, Y still in a second
// region: an amount no row keeps), a floating capacitive electrode p1 at I = 0 running the
// couple, and a capacitive gate p2 at 0.3 V. Stage 1 forced to run out after 8 accepted steps:
// its restart put back u, uLo and termV, not the plate's charge (portQEnd) the steps had moved,
// and the next stage's steps stored against it: Fe³⁺ + Q₁/F 4.2e-3 of itself off (5.8e-3 after
// 12 steps), p1 1.1e-4 V from a twin whose stage 1 took no step, converged and unwarned.
const film = () => {
  const l = aqueous(['Na+', 'Cl-', 'Fe2+', 'Fe3+'], { epsr: 0 });
  l.species.push({ name: 'e-', z: -1 });
  const d = build({
    T: 298.15, library: [l],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'Cl-': 105, 'Fe2+': 1, 'Fe3+': 1 } }), bath({ 'Na+': 100, 'Cl-': 102, 'Fe2+': 1 }, 'Cl-')],
    ports: [
      { name: 'p1', region: 'film', from: 0, to: 4e-6, I: 0, terminal: 'e-', area: 1e4, capacitance: { C: 0.2, zeroCharge: 0 }, reactions: [{ equation: 'Fe2+ = Fe3+ + e-', k0: 1e-3, alpha: 0.5 }] },
      { name: 'p2', region: 'film', from: 6e-6, to: 10e-6, V: 0.3, terminal: 'e-', area: 1e4, capacitance: { C: 0.2, zeroCharge: 0 } },
    ],
    grid: { hmin: 5e-7, hmax: 5e-7 },
  });
  d.species.push({ name: 'Y', z: 0 }, { name: 'Z', z: 0 });
  const w = d.materials.water;
  w.species.Y = { D: 1e-9, mu0: 0, cRef: 1000 };
  w.species.Z = { D: 1e-9, mu0: 0, cRef: 1000 };
  d.materials.still = { ...w, species: { ...w.species, Y: { D: 0, mu0: 0, cRef: 1000 } } };
  delete d.materials.still.species.Z;
  d.regions[0].c0.Y = 1;
  d.regions[0].c0.Z = 0.01;
  d.regions.unshift({ name: 'still', material: 'still', length: 1e-6, c0: { 'Na+': 100, 'Cl-': 105, 'Fe2+': 1, 'Fe3+': 1, Y: 2 } });
  d.bulkReactions = [{ equation: 'Y = Z', kf: { water: 1 } }];
  return new Device(d);
};

test('a steady stage that fails puts back the capacitors\' charges with the levels', () => {
  const pool = (s) => {
    s._assembleBookkeeping(Infinity);
    const i = s.model.species.findIndex((sp) => sp.name === 'Fe3+');
    let t = 0;
    for (let g = 0; g < s.nNodes; g++) t += s.model.grid.vol[g] * s.c[g * s.n + i];
    return t + s.portQ[0] / FARADAY;
  };
  const run = (N) => {
    const d = film(), s = d.solver, steps = s._steadySteps, S0 = pool(s);
    let calls = 0, moved = NaN;
    s._steadySteps = function (opts) {
      if (++calls > 1) return steps.call(this, opts);
      if (N === 0) return { converged: false, steps: 0, iterations: 0, history: [] };
      const q = steps.call(this, { ...opts, maxSteps: N });
      moved = this.portQEnd[0];
      return { ...q, converged: false }; // (out of steps, the state where they left it)
    };
    const sol = d.solve();
    s._steadySteps = steps;
    return { sol, V: sol.terminals.p1.V, drift: (pool(s) - S0) / S0, moved };
  };
  const twin = run(0), r = run(8);
  assert.ok(twin.sol.converged && r.sol.converged, r.sol.warnings.join('; '));
  assert.ok(Math.abs(r.moved) > 1e-3, `the forced stage moved p1's charge (${r.moved} C/m²)`);
  assert.ok(Math.abs(r.drift) < 1e-12, `Fe³⁺ + Q₁/F drifted by ${r.drift} of itself`);
  assert.ok(Math.abs(r.V - twin.V) < 1e-12, `p1 at ${r.V} V, its twin at ${twin.V} V`);
});

// A14(b): an islands round that fails restored the state it converged to before the round and
// returned it as converged, without the checks a converged return passes after the round
// (rateHeld, settle). Forced (every factorisation with the islands' rows in refused), the
// myelinated axon's cold solve came back at its false root, φ +0.424 V, Na⁺ at 2.07e10 mol/m³,
// which rateHeld refuses, converged.
test('a failed islands round\'s restored state passes the checks every converged return does', () => {
  const def = JSON.parse(readFileSync(new URL('./fixtures/myelin-k-electrode.json', import.meta.url), 'utf8'));
  const run = (force) => {
    const d = new Device(structuredClone(def)), v = d.solver, factor = v._factor, newton = v.newton, calls = [];
    v._factor = function (...a) {
      if (force && this.islandsOn && this.steady) throw new Error('refused with the islands\' rows in');
      return factor.apply(this, a);
    };
    v.newton = function (dt, opts) {
      const r = newton.call(this, dt, opts);
      if (dt === Infinity) calls.push({ converged: r.converged, error: r.error, failed: this.islandsFailed });
      return r;
    };
    const s = d.solve();
    return { s, calls, mid: s.x.length >> 1 };
  };
  const twin = run(false), f = run(true);
  const first = f.calls[0];
  assert.ok(first.failed && !first.converged && /Na\+ .* loses .* through its links/.test(first.error), JSON.stringify(first));
  assert.ok(f.s.converged, f.s.warnings.join('; '));
  assert.ok(Math.abs(f.s.phi[f.mid] - twin.s.phi[twin.mid]) < 1e-9, `φ ${f.s.phi[f.mid]} V, the twin's ${twin.s.phi[twin.mid]} V`);
  assert.ok(Math.abs(f.s.c['Na+'][f.mid] / twin.s.c['Na+'][twin.mid] - 1) < 1e-9, `Na⁺ ${f.s.c['Na+'][f.mid]}`);
});

// A27(b): a refusal (SolverError) in the islands round's own assembly (imbalance) escaped
// accept(), past fail(), out of solve() (which restored the state), and left the islands' rows
// in. It fails the Newton call instead, which the steady solve then gets past. (The GaAs n/p⁺/n
// stack at −0.3 V, warm from −0.2 V.)
const stack = () => {
  const mat = 'GaAs', n = (N, L) => layer(mat, L, { donors: units.perCm3(N) });
  return new Device(build({ T: 300, library: [semiconductor(mat)], stack: [ohmic(0), n(1e17, 400e-9), layer(mat, 200e-9, { acceptors: units.perCm3(5e18) }), n(3e17, 900e-9), ohmic(0)] }));
};

test('a refusal in the islands round\'s assembly fails its Newton call, not the solve', () => {
  const [A, B] = [stack(), stack()];
  for (const d of [A, B]) {
    d.set({ contacts: { right: { V: -0.2 } } });
    d.solve();
    d.set({ contacts: { right: { V: -0.3 } } });
  }
  const s = A.solver, assemble = s.assemble, newton = s.newton;
  let thrown = 0, after = null;
  // (The first assembly with the rows in is the round's measure of the balance: the rows go in
  // only after it.)
  s.assemble = function (dt) {
    if (this.islandsOn && thrown === 0) {
      thrown++;
      throw new SolverError('refused: the islands\' assembly');
    }
    return assemble.call(this, dt);
  };
  s.newton = function (dt, opts) {
    const was = thrown, r = newton.call(this, dt, opts);
    if (was === 0 && thrown === 1) after = { converged: r.converged, error: r.error, on: this.islandsOn, asked: !!opts?.islands };
    return r;
  };
  const sa = A.solve(), sb = B.solve();
  assert.equal(thrown, 1);
  assert.ok(after && !after.converged && after.error === 'refused: the islands\' assembly' && after.on === after.asked, JSON.stringify(after));
  assert.ok(sa.converged && sb.converged);
  let dphi = 0;
  for (let j = 0; j < sa.phi.length; j++) dphi = Math.max(dphi, Math.abs(sa.phi[j] - sb.phi[j]));
  assert.ok(dphi < 1e-9, `${dphi} V from the twin`);
});

// A27(c): impedance() made the steady state's readouts (a failed one's warnings, a solved one's)
// outside the try that restores: a throw there left the device at the steady state it had moved
// to, u 2.46 off, 15 STATE fields from a twin's, its next advance −7.4e-7 A/m² for the twin's
// 4.0e-6.
test('impedance() leaves the device as it was when a readout of its steady state throws', () => {
  const differ = (a, b) => {
    const j = (v) => JSON.stringify(v, (_, w) => (ArrayBuffer.isView(w) ? Array.from(w) : w));
    return STATE.filter((k) => j(a[k]) !== j(b[k]));
  };
  for (const failed of [false, true]) {
    const [A, B] = [stack(), stack()];
    for (const d of [A, B]) {
      d.solve();
      d.set({ contacts: { right: { V: 0.1 } } });
      d.advance(1e-12);
    }
    const s = A.solver, fluxes = s.segmentFluxes, steady = s.solveSteady;
    let armed = false, thrown = 0;
    s.solveSteady = function (opts) {
      const r = steady.call(this, opts);
      armed = true;
      return failed ? { ...r, converged: false } : r;
    };
    s.segmentFluxes = function () {
      if (armed && thrown++ === 0) throw new SolverError('refused: a readout');
      return fluxes.call(this);
    };
    assert.throws(() => A.impedance([1]), /refused: a readout/);
    s.segmentFluxes = fluxes;
    s.solveSteady = steady;
    assert.deepEqual(differ(A.solver, B.solver), [], failed ? 'after a failed steady state' : 'after a solved one');
    assert.equal(A.advance(1e-9).current, B.advance(1e-9).current);
  }
});

// A25: the floating march's claims. Held solves of the round-4 stub (src/steady.js's
// _floatingContinuation, its held solves replaced by an I(V); the floated ones failing).
const stub = (I, target, V0) => {
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), end = { terminal: 'Cl-', phi: 'neutral', species: { 'Cl-': 'equilibrium', 'K+': 'blocked' } };
  const d = new Device({
    species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, ...end }, right: { I: target, ...end } }, grid: { minCells: 8 },
  });
  const s = d.solver, held = [];
  s.termV[1] = V0;
  s._steadyFromHere = function () {
    const V = this.terms[1].drive.src.value;
    held.push(V);
    if (held.length > 500) throw new Error('runaway');
    this.termV[1] = V;
    this.termI[1] = I(V);
    return { converged: Number.isFinite(I(V)), steps: 1, iterations: 1, history: [] };
  };
  s._steadySteps = () => ({ converged: false, steps: 1, iterations: 1, history: [] });
  s._bound = () => false;
  return { s, held, r: s._floatingContinuation({}, { steps: 0, iterations: 0, history: [] }, 1) };
};

test('the floating march: an exact root at its start, both ends matching, and the range of currents it passed', () => {
  // A start that passes the target exactly, its neighbours failing: refused (it noticed fa = 0
  // only after a neighbour converged).
  const a = stub((V) => (V === 0.5 ? 1 : NaN), 1, 0.5);
  assert.ok(a.r.converged && a.s.termV[1] === 0.5, JSON.stringify(a.r));
  // A flat stretch matching the target, the march starting inside it: regula falsi divided by
  // fb − fa = 0 and held NaN V.
  const flat = (V) => (V < 0.4 ? V - 0.4 : V > 0.6 ? V - 0.6 : 0) + 1, b = stub(flat, 1, 0.5);
  assert.ok(b.r.converged && b.s.termV[1] === 0.5 && b.held.every(Number.isFinite), JSON.stringify(b.held));
  // A peak of 10 A/m² at 1 V, for a target of 20: what was passed read "0 to 0", the currents at
  // the two extreme voltages.
  const c = stub((V) => 10 * Math.exp(-(((V - 1) / 0.5) ** 2) / 2), 20, 0);
  assert.ok(!c.r.converged && c.s.unreached);
  const { V, I, nearing } = c.s.unreached;
  assert.ok(V[0] < -20 && V[1] > 20 && I[0] === 0 && Math.abs(I[1] - 10) < 0.01 && nearing === null, JSON.stringify(c.s.unreached));
});

test('the floating march claims no limiting current where its current was still nearing the target at 20 V', () => {
  // A 10 µm Si resistor (donors 1 mol/m³), driven at what 100 V passes, its floated solves failing
  // until the march began: refused as "a limiting current, or kinetics too slow", where held at
  // 100 V it passes it (1.456e11 A/m², L/σ from the electrons' D).
  const def = (right) => build({ library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 10e-6, { donors: 1 }), ohmic(right)] });
  const I100 = new Device(def(100)).solve().terminals.right.current;
  const d = new Device(def({ I: I100 })), s = d.solver, from = s._steadyFromHere, steps = s._steadySteps;
  let held = 0;
  s._steadyFromHere = function (opts) {
    if (this.terms[1].drive.kind === 'V') held++;
    return from.call(this, opts);
  };
  s._steadySteps = function (opts) {
    return this.terms[1].drive.kind === 'I' && held === 0 ? { converged: false, steps: 1, iterations: 0, history: [] } : steps.call(this, opts);
  };
  const r = d.solve();
  assert.equal(r.converged, false);
  assert.ok(s.unreached && s.unreached.nearing.at > 20 && Math.abs(s.unreached.nearing.guess - 100) < 1, JSON.stringify(s.unreached));
  const w = r.warnings.find((x) => x.startsWith('right:')) ?? '';
  assert.ok(/still nearing it at 26\.2 V \(at its last stride's pace, it would pass it near 100 V\)/.test(w) && !/limiting current, or kinetics/.test(w), w);
});

// A13: the residual check after a converged update reads the block system (s.res) only; a pin's
// residual (a conserved amount's) and a floating terminal's circuit row are solved beside it, in
// the border. Disturbed on that check's assembly alone, each was taken as solved.
test('newton: a converged update with a pin or a circuit row left unsolved is refused, and says which', () => {
  const end = { terminal: 'K+', species: { 'K+': 'equilibrium', 'Cl-': 'blocked' }, phi: 'neutral' };
  const cases = [
    // (Cl⁻ blocked at both ends: its amount, a pin)
    ['a pin', new Device({
      species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
      materials: { w: { epsr: 0, species: { 'K+': { D: 1e-9, cRef: 1, mu0: 0 }, 'Cl-': { D: 1e-9, cRef: 1, mu0: 0 } } } },
      regions: [{ name: 'a', material: 'w', length: 1e-6, c0: { 'K+': 1, 'Cl-': 1 } }],
      contacts: { left: { V: 0, ...end }, right: { V: 0.01, ...end } }, grid: { minCells: 8 },
    }), (s) => {
      const p = s.pins.find((q) => !s.pinnedIslands.includes(q));
      let g = 0;
      for (let j = 0; j < p.len; j++) g += Math.abs(p.w[j]);
      p.res += 1e-3 * g;
    }],
    // (the right contact driven by a current: its circuit row)
    ['a circuit row', new Device(build({ library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 10e-6, { donors: 1 }), ohmic({ I: 1e6 })] })), (s) => {
      let g = 0;
      for (const x of s.termC[1]) g += Math.abs(x);
      s.termRes[1] += 1e-3 * (g + Math.abs(s.termDI[1]) * s.VT);
    }],
  ];
  for (const [what, d, disturb] of cases) {
    assert.ok(d.solve().converged, what);
    const s = d.solver, orig = s.newton, assemble = s.assemble;
    let r = null, n = 0;
    s.newton = function (dt, opts) {
      if (r || dt !== Infinity) return orig.call(this, dt, opts);
      s.assemble = function (t) {
        assemble.call(this, t);
        if (++n === 2) disturb(this); // (the residual check's assembly, after the first update)
      };
      try {
        return (r = newton(s, dt, { ...opts, tol: 1e9 }));
      } finally {
        s.assemble = assemble;
      }
    };
    d.solve();
    s.newton = orig;
    assert.ok(r && !r.converged, `${what}: ${JSON.stringify(r)}`);
    assert.match(r.error, /a row the bordered solve holds .* is 1\.0e-3 thermal units from solved/, what);
  }
});

// A29: a contact's ramp over a span whose millionth underflows (1e-320 V), its target always
// failing, quartered its stride to 0 and held there without end; and a giant step within 1e6 of
// overflowing (1e303 s, a species at D = 4e-309) grew ×10 to ∞, a step without the storage term
// that keeps what no row does.
test('a ramp whose stride underflows ends, and a giant step grows no further than the largest finite one', () => {
  const end = { terminal: 'K+', species: { 'K+': 'equilibrium', 'Cl-': 'blocked' }, phi: 'neutral' };
  const d = new Device({
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { w: { epsr: 0, species: { 'K+': { D: 1e-9, cRef: 1, mu0: 0 }, 'Cl-': { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'w', length: 1e-6, c0: { 'K+': 1, 'Cl-': 1 } }],
    contacts: { left: { V: 0, ...end }, right: { V: 0, ...end } }, grid: { minCells: 8 },
  });
  d.solve();
  for (const span of [5e-324, 1e-320]) {
    const s = d.solver, steps = s._steadySteps;
    let calls = 0;
    s._steadySteps = function (opts) {
      if (++calls > 200) throw new Error('runaway');
      return this.sourceOverride.get(1) === span ? { converged: false, steps: 1, iterations: 0, history: [] } : steps.call(this, opts);
    };
    s.steady = true;
    try {
      assert.equal(s._continuation({}, 1, 0, span, { converged: true, steps: 0, iterations: 0, history: [], ramped: true }, true).converged, false);
    } finally {
      s.steady = false;
      s._steadySteps = steps;
    }
  }
  const gate = new Device({
    T: 300,
    species: [{ name: 'e-', z: -1 }, { name: 'Y', z: 0 }, { name: 'Z', z: 0 }, { name: 'W', z: 0 }],
    materials: {
      m1: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 1e-9, mu0: 0, cRef: 1e3 }, Z: { D: 1e-9, mu0: 0, cRef: 1e3 }, W: { D: 4e-309, mu0: 0, cRef: 1e3 } } },
      m2: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 0, mu0: 0, cRef: 1e3 } } },
    },
    regions: [
      { name: 'A', material: 'm1', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 1, Z: 1, W: 1e-3 } },
      { name: 'B', material: 'm2', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 1 } },
    ],
    interfaces: [{ dipole: 0 }],
    bulkReactions: [{ equation: 'Y = Z', kf: { m1: 1 } }],
    contacts: { left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0.1, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
    grid: { minCells: 40 },
  });
  const s = gate.solver, step = s.step, dts = [];
  s.step = function (dt, opts) {
    dts.push(dt);
    return step.call(this, dt, opts);
  };
  // (a tol no step meets: dt grows on to its cap, every step of it finite)
  const r = gate.solve({ tol: 1e-300, maxSteps: 12 });
  s.step = step;
  assert.equal(r.converged, false);
  assert.ok(dts.length > 0 && dts.every(Number.isFinite) && dts.includes(Number.MAX_VALUE), dts.slice(0, 12).join(' '));
});
