import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, metal, square, bath, aqueous } from '../src/kit.js';
import { choice, device } from '../tools/lawsweep-devices.js';
import { equilibrateFaces } from '../src/start.js';

const RT = GAS_CONSTANT * 298.15;

// For a linear single-pole response with time constant τ, backward Euler with step dt gives
// exactly x_n = x_∞ (1 − (1 + dt/τ)^(−n)), so time-stepping error drops out of the comparison.
const beStep = (n, dt, tau) => 1 - (1 + dt / tau) ** -n;

test('floating island charging after a gate step: RC time constant, Gauss, conservation', () => {
  const Dp = 1.33e-9, Dm = 2.03e-9, c = 10, L = 3e-6, Cg = 0.2, epsr = 78.5;
  const gate = (V) => ({ V, phi: { type: 'capacitive', C: Cg }, zeroCharge: 0 });
  const def = {
    species: [
      { name: 'Na+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr, species: { 'Na+': { D: Dp, mu0: -261.9e3 }, 'Cl-': { D: Dm, mu0: -131.2e3 } } } },
    regions: [{ material: 'water', length: L, c0: { 'Na+': c, 'Cl-': c } }],
    contacts: { left: gate(0), right: gate(0) },
    grid: { hmin: 0.05e-9, hmax: 20e-9, ratio: 1.1 },
  };
  const dev = new Device(def);
  dev.solve();
  // Small step on the left gate. Series circuit: gate—C_s—R_bulk—C_s—gate, with
  // C_s = C_g ∥series∥ ε/λ_D and R_bulk = L/σ, so τ = R_bulk·C_s/2.
  const dV = 1e-3;
  dev.set({ contacts: { left: { V: dV } } });
  const qInf = new Device(dev.def).solve().plates.left.charge;
  const lam = Math.sqrt((epsr * EPS0 * RT) / (2 * FARADAY * FARADAY * c));
  const Cs = 1 / (1 / Cg + lam / (epsr * EPS0));
  const sigma = (FARADAY * FARADAY * c * (Dp + Dm)) / RT;
  const tau = ((L / sigma) * Cs) / 2;
  assert.ok(Math.abs(qInf / ((Cs * dV) / 2) - 1) < 1e-3);

  const dt = tau / 20;
  for (let n = 1; n <= 80; n++) {
    const sol = dev.step(dt);
    assert.ok(sol.converged);
    const q = sol.plates.left.charge / qInf;
    // Distributed double-layer charging adds a small correction (~λ_D/L and ~τ_D/τ).
    assert.ok(Math.abs(q - beStep(n, dt, tau)) < 1e-2, `n=${n}: ${q} vs ${beStep(n, dt, tau)}`);
    const gross = FARADAY * 2 * c * L;
    assert.ok(Math.abs(sol.charge + sol.plates.left.charge + sol.plates.right.charge) < 1e-12 * gross);
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `${st.species} drift ${st.drift}`);
    // Gate current is displacement current, equal at both ends. (Newton stops at δφ̂ ~ 1e-10,
    // which bounds the displacement difference quotient to ~1e-7 relative.)
    assert.ok(Math.abs(sol.contacts.left.current - sol.contacts.right.current) < 1e-6 * Math.abs(sol.contacts.left.current));
  }
});

test('water: homogeneous relaxation of excess H⁺/OH⁻ at rate k (c_H + c_OH)', () => {
  const kf = 1.4e8; // m³/(mol·s), i.e. 1.4e11 M⁻¹s⁻¹
  const muH2O = -237.13e3, mu0OH = -157.24e3;
  const Kw = 1e6 * Math.exp((muH2O - mu0OH) / RT);
  const ceq = Math.sqrt(Kw);
  const excess = 1e-3 * ceq; // small, so the response is linear
  const def = {
    species: [
      { name: 'H+', z: 1, cRef: 1000 },
      { name: 'OH-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 78.5, species: { 'H+': { D: 9.31e-9, mu0: 0 }, 'OH-': { D: 5.27e-9, mu0: mu0OH } } } },
    regions: [{ material: 'water', length: 50e-9, c0: { 'H+': ceq + excess, 'OH-': ceq + excess } }],
    bulkReactions: [{ nu: { 'H+': -1, 'OH-': -1, H2O: 1 }, fixed: { H2O: muH2O }, kf: { water: kf } }],
    contacts: {
      left: { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 },
      right: { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 },
    },
    grid: { minCells: 10 },
  };
  const dev = new Device(def);
  const tau = 1 / (kf * 2 * ceq);
  const dt = tau / 10;
  const m = dev.grid.nNodes >> 1;
  for (let n = 1; n <= 30; n++) {
    const sol = dev.step(dt);
    assert.ok(sol.converged);
    const left = (sol.c['H+'][m] - ceq) / excess; // fraction of the excess remaining
    const want = 1 - beStep(n, dt, tau);
    assert.ok(Math.abs(left - want) < 2e-3, `n=${n}: ${left} vs ${want}`);
  }
});

test('open system: amounts track the time-integrated contact fluxes after a bias step', () => {
  const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
  const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
  const dev = new Device({
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } } },
    regions: [
      { material: 'Si', length: 2e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'Si', length: 2e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
    ],
    contacts: { left: ohmic(0), right: ohmic(0) },
    grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
  });
  dev.solve();
  const start = dev.def;
  for (const method of ['be', 'bdf2']) {
    dev.set(start);
    dev.solve();
    dev.set({ contacts: { right: { V: 0.4 } } });
    let dt = 1e-13;
    for (let n = 0; n < 40; n++, dt *= 1.5) {
      const sol = dev.step(dt, { method });
      assert.ok(sol.converged, `step ${n}`);
      for (const st of sol.conservation) {
        assert.ok(st.connected);
        assert.ok(Math.abs(st.drift) < 1e-11, `${method} step ${n} ${st.species}: drift ${st.drift}`);
      }
    }
  }
});

// A gated island charging after a small gate step (as above, smaller for speed).
function chargingIsland() {
  const gate = (V) => ({ V, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 });
  const dev = new Device({
    species: [
      { name: 'Na+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 78.5, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } } },
    regions: [{ material: 'water', length: 1e-6, c0: { 'Na+': 10, 'Cl-': 10 } }],
    contacts: { left: gate(0), right: gate(0) },
    grid: { hmin: 0.1e-9, hmax: 20e-9, ratio: 1.15 },
  });
  dev.solve();
  dev.set({ contacts: { left: { V: 0.02 } } });
  return dev;
}
const T = 1e-6; // a few RC times
let reference; // gate charge at T from fine BDF2 steps
const referenceCharge = () => {
  if (reference === undefined) {
    const dev = chargingIsland();
    let sol;
    for (let k = 0; k < 800; k++) sol = dev.step(T / 800, { method: 'bdf2' });
    reference = sol.plates.left.charge;
  }
  return reference;
};

test('BDF2 is second order in time, backward Euler first order', () => {
  const run = (n, method) => {
    const dev = chargingIsland();
    let sol;
    for (let k = 0; k < n; k++) sol = dev.step(T / n, { method });
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `${method}: ${st.species} drift ${st.drift}`);
    return sol.plates.left.charge;
  };
  const ref = referenceCharge();
  const err = (n, method) => Math.abs(run(n, method) / ref - 1);
  const be = err(20, 'be') / err(40, 'be'), bdf2 = err(20, 'bdf2') / err(40, 'bdf2');
  assert.ok(be > 1.8 && be < 2.2, `BE ratio ${be}`);
  assert.ok(bdf2 > 3.5 && bdf2 < 4.5, `BDF2 ratio ${bdf2}`);
});

test('adaptive advance: error control, exact landing, traces, and frame budgets', () => {
  const ref = referenceCharge();
  const errs = [1e-3, 1e-5].map((tol) => {
    const sol = chargingIsland().advance(T, { tol });
    assert.ok(sol.converged && sol.done && sol.time === T);
    assert.ok(sol.trace.t.every((t, k) => k === 0 || t > sol.trace.t[k - 1]) && sol.trace.t.at(-1) === T);
    assert.ok(Math.abs(sol.trace.current.at(-1) - sol.current) < 1e-12 * Math.abs(sol.current));
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12);
    return Math.abs(sol.plates.left.charge / ref - 1);
  });
  assert.ok(errs[0] < 2e-2 && errs[1] < 1e-3 && errs[1] < errs[0] / 5, `errors ${errs}`);

  // A frame budget stops early; later calls pick up where it left off.
  const dev = chargingIsland();
  let sol = dev.advance(T, { budgetMs: 0 });
  assert.ok(!sol.done && sol.time < T);
  let frames = 1;
  while (!sol.done && frames < 1000) {
    sol = dev.advance(T, { budgetMs: 1 });
    frames++;
  }
  assert.ok(sol.done && sol.time === T && frames > 1);
  assert.ok(Math.abs(sol.plates.left.charge / ref - 1) < 2e-2);
});

test('a jump in a strictly neutral cell: advance() finds its first step (longer, not shorter, steps help there)', async () => {
  // Zn | ZnSO₄ (ε = 0) | Zn, stepped by 0.1 V. Right after a jump, a strictly neutral material's
  // field is fixed only through fluxes, which a short step's storage term dwarfs.
  const { build, layer, ohmic, half, aqueous, metal } = await import('../src/kit.js');
  const zinc = half('Zn2+ + 2 e- = Zn(s)', { 'Zn(s)': 0 });
  const face = { reactions: [{ ...zinc, k0: 1e-3 }] };
  const dev = new Device(
    build({
      library: [aqueous(['Zn2+', 'SO42-'], { epsr: 0 }), metal('Zn')],
      stack: [ohmic(0, ['e-']), layer('Zn', 1e-6), face, layer('water', 100e-6, { c0: { 'Zn2+': 100, 'SO42-': 100 } }), face, layer('Zn', 1e-6), ohmic(0, ['e-'])],
      grid: { hmin: 10e-9, hmax: 2e-6 },
    }),
  );
  dev.solve();
  dev.set({ contacts: { right: { V: 0.1 } } });
  const run = dev.advance(1e-3);
  assert.ok(run.converged && run.done, JSON.stringify({ steps: run.steps, rejected: run.rejected }));
});

test('advance() to where it already is changes nothing, and the next call still moves', () => {
  // As an animation frame with no time to add would ask.
  const dev = new Device({
    species: [{ name: 'X', z: 0, cRef: 1 }],
    materials: { m: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } },
    regions: [{ material: 'm', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { species: { X: { type: 'equilibrium', mu: 1000 } }, phi: 'neutral' } },
    grid: { minCells: 10 },
  });
  const still = dev.advance(0);
  assert.ok(still.converged && still.done && still.steps === 0);
  const moved = dev.advance(1e-4);
  assert.ok(moved.converged && moved.done && moved.steps > 0 && dev.solver.time === 1e-4);
  // A time that isn't one (dev.time, which doesn't exist) is an error; a NaN (a first frame's
  // clock, dev.time + 1) does nothing, and says so.
  assert.throws(() => dev.advance(dev.time), /tEnd is the time to reach in s, absolute/);
  const none = dev.advance(dev.time + 1);
  assert.ok(none.converged && !none.done && none.steps === 0 && dev.solver.time === 1e-4);
  assert.match(none.warnings.join(), /advance\(NaN\) did nothing/);
  assert.deepEqual(still.warnings, []);
  // advance(dt) read as a duration: called again, it would silently do nothing, so it says so.
  assert.match(dev.advance(1e-4).warnings.join(), /advance\(0\.0001\) did nothing.*advance\(sol\.time \+ dt\)/);
  assert.match(dev.advance(5e-5).warnings.join(), /did nothing/);
  assert.deepEqual(dev.advance(1e-4, { budgetMs: 10 }).warnings, []); // an animation frame with no time to add
});

test('cyclic voltammetry in a strictly neutral cell runs through its turns (round-off at breakpoints)', async () => {
  // Pt | Fe³⁺, Fe²⁺ in KCl (ε = 0) through 1 mm to a bath, swept by a triangle wave. Short steps
  // after each turn are badly conditioned (φ fixed only through fluxes), so Newton converges to a
  // round-off floor above its tolerance there; and a frame's target can land a hair past a turn.
  const { build, layer, ohmic, bath, half, aqueous, metal } = await import('../src/kit.js');
  const iron = half('Fe3+ + e- = Fe2+');
  const cell = (drive) =>
    build({
      library: [aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 }), metal('Pt')],
      stack: [ohmic(drive, ['e-']), layer('Pt', 40e-6), { reactions: [{ ...iron, k0: 1e-3 }] }, layer('water', 1e-3), bath({ 'K+': 500, 'Cl-': 525, 'Fe3+': 5, 'Fe2+': 5 }, 'Cl-')],
      grid: { hmin: 10e-9, hmax: 25e-6 },
    });
  const E0 = new Device(cell({ I: 0 })).solve().terminals.left.V, period = 0.2;
  const dev = new Device(cell(E0));
  dev.solve();
  dev.set({ contacts: { left: { V: { t: [0, period / 4, (3 * period) / 4, period], values: [E0, E0 + 0.5, E0 - 0.5, E0], repeat: true } } } });
  assert.ok(dev.advance(period / 4).converged, 'up to the first turn');
  assert.ok(dev.step(7e-6).converged, 'a short step right after it');
  // Frames of period/360, as an animation takes them: one lands 1e-13 past a turn.
  for (let f = 0; f < 720; f++) {
    const run = dev.advance(dev.solver.time + period / 360);
    assert.ok(run.converged, `frame ${f}, t = ${dev.solver.time}`);
  }
});

test('a potential step a second into a run, on a 5 nm grid: Cottrell (step lengths, not differences of times near 1 s)', async () => {
  // Pt | Fe³⁺ in KCl (ε = 0) through 1 mm to a bath, held where nothing happens, then stepped
  // 0.8 V down at t = 1 s: Fe³⁺ is reduced at the diffusion limit. Just after the jump the
  // profile is self-similar in x/√t and steps of 1e-13 s resolve it, a second into the run.
  const { build, layer, ohmic, bath, half, aqueous, metal, IONS } = await import('../src/kit.js');
  const iron = half('Fe3+ + e- = Fe2+');
  const c = 1, bulk = { 'K+': 1000, 'Cl-': 1003.00002, 'Fe3+': c, 'Fe2+': 1e-5 };
  const cell = (drive) =>
    build({
      library: [aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 }), metal('Pt')],
      stack: [ohmic(drive, ['e-']), layer('Pt', 1e-6), { reactions: [{ ...iron, k0: 10 }] }, layer('water', 1e-3), bath(bulk, 'Cl-')],
      grid: { hmin: 5e-9, hmax: 20e-6, ratio: 1.15 },
    });
  const Eeq = new Device(cell({ I: 0 })).solve().terminals.left.V; // Fe³⁺/Fe²⁺ at 1e5 : 1
  const up = Eeq + 0.2, down = Eeq - 0.6;
  const dev = new Device(cell({ t: [0, 1, 1, 2], values: [up, up, down, down] }));
  dev.solve();
  assert.ok(dev.advance(1).done);
  const run = dev.advance(1.3);
  assert.ok(run.converged && run.done && run.warnings.length === 0, JSON.stringify({ steps: run.steps, warnings: run.warnings }));
  const D = IONS['Fe3+'].D, { t, current } = run.trace;
  for (const after of [1e-3, 1e-2, 0.1, 0.3]) {
    const k = t.findIndex((tk) => tk - 1 >= after * (1 - 1e-9));
    const cottrell = FARADAY * c * Math.sqrt(D / (Math.PI * (t[k] - 1)));
    assert.ok(Math.abs(-current[k] / cottrell - 1) < 0.01, `${after} s after the step: ${-current[k]} A/m² against Cottrell's ${cottrell}`);
  }
  // Stopped short: the solution says where, and why.
  const short = dev.advance(2, { maxSteps: 3 });
  assert.ok(!short.done && short.warnings.some((w) => /advance stopped at t = .*maxSteps \(3\)/.test(w)), short.warnings.join('\n'));
});

test('a probe at a face reads the side with the species; a species absent there, or both sides, is an error', async () => {
  const { build, layer, ohmic, bath, half, aqueous, metal } = await import('../src/kit.js');
  const bulk = { 'K+': 500, 'Cl-': 525, 'Fe3+': 5, 'Fe2+': 5 };
  const dev = new Device(
    build({
      library: [aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 }), metal('Pt')],
      stack: [ohmic({ I: -10 }, ['e-']), layer('Pt', 1e-6), { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1, alpha: 0.5 }] }, layer('water', 1e-4), bath(bulk, 'Cl-')],
      grid: { hmin: 1e-8, hmax: 2e-6 },
    }),
  );
  // At the Pt | water face (x = 1 µm), Fe³⁺ is only in the water: read there, not NaN.
  const run = dev.advance(1e-3, { probes: [{ x: 1e-6, species: 'Fe3+' }, { x: 1e-6, species: 'Fe3+', region: 1 }] });
  const [face, told] = run.trace.probes;
  assert.ok(face.every(Number.isFinite) && face.every((v, k) => v === told[k]), 'the water side');
  assert.ok(face.at(-1) < 5, 'Fe³⁺ drawn down at the electrode');
  assert.throws(() => dev.advance(2e-3, { probes: [{ x: 0.5e-6, species: 'Fe3+' }] }), /'Fe3\+' is absent from region/);
  // Two layers of water: both sides have it, so say which.
  const two = new Device(
    build({
      library: [aqueous(['K+', 'Cl-'], { epsr: 0 })],
      stack: [bath({ 'K+': 100, 'Cl-': 100 }, 'Cl-'), layer('water', 1e-5, { name: 'a' }), layer('water', 1e-5, { name: 'b' }), bath({ 'K+': 10, 'Cl-': 10 }, 'Cl-')],
    }),
  );
  assert.throws(() => two.advance(1e-3, { probes: [{ x: 1e-5, species: 'K+' }] }), /on the face between "a" and "b": give region/);
});

// A terminal at no current that alone feeds what it holds (a host's electrons, a floating gate's)
// keeps it exactly however long the steps. A host filled linearly from x = 0.1 to 0.9 relaxes to
// x̄ = 0.5, where its OCV table (symmetric about it) gives 0.4 V; read by the flux through its
// contact, it had drifted to 0.3995 V by 1e10 s and failed by 1e14 s.
test('a closed host at open circuit keeps its charge through any step: its OCV at the mean filling, to 1e14 s', () => {
  const x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95], E = x.map((v) => 0.4 - (RT / FARADAY) * Math.log(v / (1 - v)));
  const d = new Device({
    species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: { host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }] } },
    regions: [{ material: 'host', length: 1e-6, c0: { 'Li+': { x: [0, 1e-6], values: [3000, 27000] } } }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
    grid: { hmin: 2e-8, hmax: 5e-8 },
  });
  d.advance(1e6, { tol: 1e-4 });
  for (const t of [1e10, 1e14]) {
    const s = d.advance(t, { tol: 1e-4 });
    assert.ok(s.done && s.steps < 30, `t ${t}: ${s.steps} steps`);
    assert.ok(Math.abs(s.terminals.right.V - 0.4) < 1e-12, `t ${t}: V ${s.terminals.right.V}`);
    assert.ok(Math.abs(s.conservation.find((c) => c.species === 'e-').drift) < 1e-12);
  }
});

// A cell at rest: two hosts, Li⁺ crossing into each from an electrolyte whose salt gradient
// relaxes, the cathode's collector at no current. Each host keeps its Li (no current), so it ends
// at its OCV at its own filling, and the cell at their difference. The hosts' φ̂ is a gauge that
// Newton never moves; extrapolated into each step's first guess, it wandered off geometrically
// as the steps grew (1e11 thermal units by 1e5 s), and the levels measured from it lost their
// digits: every step past ~1e5 s failed, and past ~400 s at tol 1e-8.
test('a battery at rest settles to the difference of its hosts\' OCVs, out to 1e9 s at any tolerance', () => {
  const VT = RT / FARADAY, x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];
  const ocv = (E0) => (v) => E0 - 1.5 * VT * Math.log(v / (1 - v)) - 0.1 * v;
  const host = (E0, cMax) => ({
    epsr: 0,
    species: { 'Li+': { D: 1e-14, mu0: 0, cRef: cMax }, 'e-': { D: 1e-4, mu0: 0, cRef: cMax } },
    statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax, ocv: { x, E: x.map(ocv(E0)), muRef: 0 } }],
  });
  const transfer = { species: { 'Li+': 'blocked' }, reactions: [{ equation: 'Li+(left) = Li+(right)', k0: 1e-3, alpha: 0.5 }] };
  const salt = { x: [2e-6, 7e-6], values: [500, 1500] };
  const cell = () => new Device({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }, { name: 'e-', z: -1, cRef: 1000 }, { name: 'X-', z: -1, cRef: 1000 }],
    materials: { anode: host(0.2, 30000), cathode: host(3.8, 50000), elyte: { epsr: 0, species: { 'Li+': { D: 1e-10, mu0: 0 }, 'X-': { D: 1.5e-10, mu0: 0 } } } },
    regions: [
      { material: 'anode', length: 2e-6, c0: { 'Li+': 0.8 * 30000 } },
      { material: 'elyte', length: 5e-6, c0: { 'Li+': salt, 'X-': salt } },
      { material: 'cathode', length: 2e-6, c0: { 'Li+': 0.4 * 50000 } },
    ],
    interfaces: [transfer, transfer],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
    },
    grid: { hmin: 2e-8, hmax: 2e-7 },
  });
  const V = ocv(3.8)(0.4) - ocv(0.2)(0.8);
  for (const tol of [1e-4, 1e-8]) {
    const d = cell();
    const s = d.advance(1e9, { tol, maxSteps: 5000 });
    assert.ok(s.done, `tol ${tol}: ${s.warnings.join(' ')}`);
    assert.ok(s.rejected < 50, `tol ${tol}: ${s.rejected} steps rejected`);
    assert.ok(Math.abs(s.terminals.right.V - V) < 1e-12, `tol ${tol}: V ${s.terminals.right.V} for ${V}`);
  }
});

// Short calls late in a run: a GITT pulse of 0.01 s at t = 1e9 s. advance() snapped any remainder
// under 1e-10 of the absolute time, as if it were round-off, so the whole pulse was skipped (the
// time moved, the state didn't), and from 1e10 s a pulse train lost 17% of its charge.
test('advance() integrates short calls late in a run: the charge of every pulse, from t = 1e10 s', () => {
  const cMax = 2e4, L = 5e-6, x = [0.05, 0.3, 0.7, 0.95];
  const d = new Device({
    species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: { host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: cMax }, 'e-': { D: 1e-4, mu0: 0, cRef: cMax } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax, ocv: { x, E: [0.6, 0.2, 0.15, 0.1], muRef: 0 } }] } },
    regions: [{ material: 'host', length: L, c0: { 'Li+': 0.3 * cMax } }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
  });
  d.advance(1e10, { tol: 1e-3 });
  d.set({ contacts: { right: { I: -1 } } });
  const li = (s) => s.conservation.find((c) => c.species === 'Li+').amount, a0 = li(d.solution());
  let total = 0;
  for (const dt of [0.01, 0.1, 0.3, 1, 3, 10]) {
    const s = d.advance(d.solution().time + dt, { tol: 1e-3 });
    assert.ok(s.done, `dt ${dt}`);
    total += dt;
  }
  const gained = li(d.solution()) - a0;
  assert.ok(Math.abs(gained / (total / FARADAY) - 1) < 1e-6, `gained ${gained} of ${total / FARADAY} mol/m²`);
});

// A MOS capacitor's gate, settled at 1 V and floated: nothing moves, and its charge stays put.
// Read by the flux into its metal, the charging current fell below that flux's round-off at
// steps of ~0.03 s, every longer step came out singular, and 100 s took 5000 steps.
test('a floating gate at rest steps as a held one does, its charge kept', () => {
  const lib = { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } };
  const def = build({
    T: 300,
    library: [semiconductor('Si'), metal('Al'), lib],
    stack: [ohmic(0, ['e-']), layer('Al', units.nm(20)), { phi: { type: 'capacitive', C: 10 }, zeroCharge: -0.95 }, layer('SiO2', units.nm(10)), { dipole: 0 }, layer('Si', units.um(1), { name: 'Si', acceptors: units.perCm3(1e17) }), ohmic(0)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    ports: [{ name: 'channel', region: 'Si', from: 0, to: units.nm(1), V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' } }],
    grid: { hmin: units.nm(0.25), hmax: units.nm(50), ratio: 1.15 },
  });
  const run = (float) => {
    const d = new Device(def);
    d.set({ contacts: { left: { V: 1 } } });
    const s0 = d.solve();
    if (float) d.set({ contacts: { left: { I: 0 } } });
    let steps = 0, rejected = 0, s;
    for (const t of [1e-6, 1e-3, 1, 10, 100]) {
      s = d.advance(t);
      steps += s.steps;
      rejected += s.rejected;
    }
    return { steps, rejected, s, s0 };
  };
  const held = run(false), floating = run(true);
  assert.ok(floating.steps <= held.steps && floating.rejected === 0, `floating ${floating.steps} steps (${floating.rejected} rejected), held ${held.steps}`);
  const D = (s) => s.interfaces[0].D;
  assert.ok(Math.abs(D(floating.s) / D(floating.s0) - 1) < 1e-9 && Math.abs(floating.s.terminals.left.V - 1) < 1e-9, `D ${D(floating.s)} vs ${D(floating.s0)}, V ${floating.s.terminals.left.V}`);
});

// A metal gate behind a resistance R, across an oxide from a held plate: an RC circuit,
// C = 1/(2/C_face + t_ox/ε). Backward Euler's response to a 1 V step is exactly
// Q_n = Q_∞ (1 − (1 + dt/τ)^−n), τ = RC. Read by the flux into its metal, the charging current
// was lost below that flux's round-off: fixed steps were 0.9% off, silently, and advance() took
// 3280 steps for 5τ, half rejected.
test('a gate behind a resistance charges as RC: exactly backward Euler on steps of τ/20', () => {
  const tox = units.nm(10), Cf = 10, C = 1 / (2 / Cf + tox / (3.9 * EPS0)), R = 10 / C, tau = R * C;
  const face = { phi: { type: 'capacitive', C: Cf }, zeroCharge: 0 };
  const def = (V) => build({ T: 300, library: [metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }], stack: [ohmic(V, ['e-']), layer('Al', units.nm(20)), face, layer('SiO2', tox), face, layer('Al', units.nm(20)), ohmic(0, ['e-'])] });
  const qInf = new Device(def(1)).solve().interfaces[0].D;
  assert.ok(Math.abs(qInf / C - 1) < 1e-12, `Q∞ ${qInf} vs C ${C}`);
  const start = () => {
    const d = new Device(def(0));
    d.solve();
    return d.set({ contacts: { left: { V: 1, R } } });
  };
  const d = start(), dt = tau / 20;
  let D = 0;
  for (let n = 1; n <= 40; n++) {
    const s = d.step(dt);
    assert.ok(s.converged && Math.abs(s.interfaces[0].D / qInf - beStep(n, dt, tau)) < 1e-12, `n=${n}: ${s.interfaces[0].D / qInf} vs ${beStep(n, dt, tau)}`);
    // Its current, (1 V − V)/R, is the charging of each step exactly, however few iterations
    // Newton took (read as the flux into its metal, 1.5e-4 off where it stopped early).
    const I = s.terminals.left.current, dD = (s.interfaces[0].D - D) / dt;
    D = s.interfaces[0].D;
    assert.ok(Math.abs(I / dD - 1) < 1e-9 && Math.abs(I + s.terminals.right.current) < 1e-12 * I, `n=${n}: I ${I}, dD/dt ${dD}, right ${s.terminals.right.current}`);
  }
  // Then let go: { I: 0 } drops the resistance with the voltage source.
  d.set({ contacts: { left: { I: 0 } } });
  assert.ok(Math.abs(d.advance(d.solution().time + 100 * tau).interfaces[0].D - D) < 1e-12 * D);
  const r = start().advance(5 * tau);
  assert.ok(r.steps < 100 && r.rejected === 0, `${r.steps} steps, ${r.rejected} rejected`);
  assert.ok(Math.abs(r.interfaces[0].D / qInf - (1 - Math.exp(-5))) < 1e-3);
});

test('a mistyped option or contact is an error, never ignored', async () => {
  const { build, layer, ohmic, semiconductor } = await import('../src/kit.js');
  const dev = new Device(build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donors: 16.6 }), layer('Si', 1e-6, { acceptors: 1.66 }), ohmic(0.5)], grid: { hmin: 1e-9, hmax: 2e-8 } }));
  const I = dev.solve().current;
  assert.throws(() => dev.advance(1e-6, { dtmax: 1e-9 }), /advance: no option 'dtmax' \(did you mean 'dtMax'\?\)/);
  assert.throws(() => dev.solve({ tolerance: 1 }), /solve: no option 'tolerance'; the options are maxSteps, tol, continuation/);
  assert.throws(() => dev.impedance([1], { termnal: 'left' }), /impedance: no option 'termnal'/);
  // After a solve a drive changes in place; a misspelt contact must still be refused, and leave
  // the device as it was (the next, valid change works).
  assert.throws(() => dev.set({ contacts: { rigth: { V: 0.6 } } }), /contacts\.rigth/);
  assert.ok(Math.abs(dev.solve().current / I - 1) < 1e-12);
  dev.set({ contacts: { right: { V: 0.6 } } });
  assert.ok(Math.abs(dev.solve().current) > 10 * Math.abs(I));
});

test('two strictly neutral regions out of equilibrium at their face: a transient starts where it reaches it at t = 0⁺, each ion\'s amount kept, and runs', () => {
  // (The law sweep's device 843: an immobile S⁺ at 1 and 0.1 mol/m³ either side, so c0's
  // c_Na·c_Cl differs across a face that holds both in equilibrium.) The first step, made to
  // take that jump, met an exactly singular pivot. The first step from a cold start takes it
  // (equilibrateFaces), a steady solve doesn't: there it was the worse start (stress react 381).
  const def = device(choice(843), { charged: true }), d = new Device(def), v = d.model.grid.vol;
  assert.equal(d.solver.coldStart, true);
  const e = new Device(def);
  equilibrateFaces(e.solver);
  const s0 = e.solution();
  const gL = d.model.grid.regionEnd[0], gR = gL + 1, total = (name, c) => c.reduce((t, x, g) => t + v[g] * x, 0);
  const given = (name) => def.regions.reduce((t, reg, r) => t + (reg.c0[name] ?? 0) * d.model.grid.vol.slice(d.model.grid.regionStart[r], d.model.grid.regionEnd[r] + 1).reduce((a, b) => a + b, 0), 0);
  for (const name of ['Na+', 'K+', 'Cl-', 'S+']) assert.ok(Math.abs(total(name, s0.c[name]) / given(name) - 1) < 1e-12, `${name}: ${total(name, s0.c[name])} for ${given(name)}`);
  const level = (name, g) => s0.mu[name][g];
  for (const name of ['Na+', 'K+', 'Cl-']) assert.ok(Math.abs(level(name, gL) - level(name, gR)) < 1e-6 * 8.314 * 298, `${name}'s μ̄ across the face`);
  assert.equal(s0.c['S+'][gL], 1);
  // The first step starts there (its start state is the history's u): the face's levels equal.
  d.step(1e-9);
  const { M, n, blockOfNode, model } = d.solver, u0 = d.solver.history[0].u;
  for (const name of ['Na+', 'K+', 'Cl-']) {
    const i = model.species.findIndex((sp) => sp.name === name);
    assert.ok(Math.abs(u0[blockOfNode[gL] * M + 1 + i] - u0[blockOfNode[gR] * M + 1 + i]) < 1e-6, `${name}'s level across the face at the first step's start`);
  }
  assert.ok(n > 0 && d.solver.coldStart === false);
  const s = d.advance(100, { tol: 1e-9 });
  assert.ok(s.converged && s.time === 100, s.warnings?.join('; '));
});

test('a step that fails leaves the state as it was: the start-of-step values too, which a solution reads its contacts\' fluxes against', () => {
  // (From a cold review of src/timeSteps.js.) A failed step left its start-of-step values (c*)
  // behind, and the solution read after it gave a host's contact flux 8× off at the same state.
  const def = {
    species: [{ name: 'X', z: 0, cRef: 1 }],
    materials: { m: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } },
    regions: [{ material: 'm', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { species: { X: { type: 'exchange', k: 1e-3, mu: 3000 } }, phi: 'neutral' } },
    grid: { minCells: 20 },
  };
  for (const method of ['be', 'bdf2']) {
    const d = new Device(def);
    for (let k = 0; k < 5; k++) d.step(2e-5, { method });
    const before = d.solution().contacts.left.flux.X;
    assert.equal(d.solver.step(1e3, { method, maxIter: 1 }).converged, false); // (made to fail: one Newton iteration)
    assert.equal(d.solution().contacts.left.flux.X, before, method);
  }
});

test('a breakpoint just past the time asked for: advance stops at that time, not at the breakpoint', () => {
  // (From a cold review of src/timeSteps.js.) A square wave's edge at 0.7000000000000001 took
  // advance(0.7) to it, and the next frame's advance(0.7) said it did nothing.
  // An MOS-like stack, its left contact behind R driven by the wave.
  const tox = units.nm(10), Cf = 10, C = 1 / (2 / Cf + tox / (3.9 * EPS0)), face = { phi: { type: 'capacitive', C: Cf }, zeroCharge: 0 };
  const dev = new Device(build({ T: 300, library: [metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }], stack: [ohmic(0, ['e-']), layer('Al', units.nm(20)), face, layer('SiO2', tox), face, layer('Al', units.nm(20)), ohmic(0, ['e-'])] }));
  dev.solve();
  dev.set({ contacts: { left: { V: square({ high: 1, low: 0, period: 0.1 }), R: 0.05 / C } } });
  for (let k = 1; k <= 10; k++) {
    const t = k / 10, s = dev.advance(t);
    assert.ok(s.converged && s.time === t, `advance(${t}) ended at ${s.time}`);
  }
});

test('a capacitive contact that starts uncharged far from the equilibrium of the ion it passes: the transient fails, saying why and what to do', () => {
  // (The law sweep's device 841: a contact at I = 0 with a double layer and a K⁺ link whose
  // outside level suits the charged steady state, 1.47 V from where the uncharged start is.) The
  // link closes the gap within ~1e-13 s, below any step. Asked for, that's right to fail on; the
  // message says how to start in equilibrium instead, and that works.
  const def = device(choice(841), { charged: true }), s = new Device(def).advance(1e3, { tol: 1e-9 });
  assert.ok(!s.converged);
  assert.match(s.warnings.join(' '), /contact left's K\+ link is 1\.\d\d V from equilibrium.*solve\(\) first/);
  const d = new Device(def);
  assert.ok(d.solve().converged);
  const late = d.advance(1e3, { tol: 1e-9 });
  assert.ok(late.converged && late.time === 1e3);
});

test('a floating capacitive electrode whose starting coverages balance its reactions far from where it starts: solve() finds it, keeping its sites and charge; advance() fails, saying so and what to do, which works', () => {
  // (Gemini 4 Argon's law8.) Uncharged at −1.69 V, its surface balances 1.6–1.7 V away. The steady
  // solve failed too, until the contact's current took the port's ∂/∂V (termDX): its Newton
  // never found the state, which keeps the sites (Aads + Bads) and the charge the surface gave
  // back (F × the Bads and OHads it reduced).
  const lib = aqueous(['Na+', 'OH-', 'Fe2+'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const G = 1e-5, L = 10e-6;
  const def = (drive) =>
    build({
      T: 298.15,
      library: [lib],
      stack: [{}, layer('water', L, { name: 'film' }), bath({ 'Na+': 100, 'OH-': 102, 'Fe2+': 1 }, 'OH-')],
      ports: [{
        name: 'metal', region: 'film', ...drive, terminal: 'e-', area: 1e4, capacitance: { C: 0.2, zeroCharge: 0 },
        surface: { Aads: { mu0: 0, capacity: G, theta0: 0.2 }, Bads: { mu0: 3e3, capacity: G, theta0: 0.1 }, OHads: { mu0: 0, capacity: G, theta0: 0.05 } },
        reactions: [{ equation: 'Aads + OH- = Bads + e-', k0: 1e-4, alpha: 0.5 }, { equation: 'OH- = OHads + e-', k0: 1e-4, alpha: 0.5, bare: true }],
      }],
      grid: { hmin: L / 20, hmax: L / 20 },
    });
  const cold = new Device(def({ I: 0 })), s = cold.solve(), sv = cold.solver;
  assert.ok(s.converged);
  const [A, B, OH] = [0, 1, 2].map((j) => sv.surfaceAmount(j));
  assert.ok(Math.abs((A + B) / (G * L * 1e4 * 0.3) - 1) < 1e-9, `${A + B}`);
  assert.ok(Math.abs(sv.portQ[0] / (FARADAY * (G * L * 1e4 * 0.15 - B - OH)) - 1) < 1e-6, `${sv.portQ[0]}`);
  const a = new Device(def({ I: 0 })).advance(10, { tol: 1e-6 });
  assert.ok(!a.converged);
  assert.match(a.warnings.join(' '), /The start is violent: port metal's reactions/);
  const d = new Device(def({ V: -0.1 }));
  assert.ok(d.solve().converged);
  d.set({ ports: { metal: { I: 0 } } });
  const f = d.solve();
  assert.ok(f.converged && Math.abs(f.terminals.metal.V + 0.1) < 1e-6, `${f.terminals.metal.V}`);
});

// From a cold review of src/start.js: what a start must be.
test('the start: a floating electrode behind a resistance starts where its reactions pass what the circuit says', () => {
  // Started at no reaction current, both solves failed for R ≤ 10 Ω·m².
  const lib = aqueous(['Fe2+', 'H+', 'Cl-'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const c = { 'Fe2+': 10, 'H+': 10, 'Cl-': 30 };
  const def = build({
    T: 298.15, library: [lib], stack: [bath(c, 'Cl-'), layer('water', 10e-6, { name: 'film' }), {}],
    ports: [{ name: 'iron', region: 'film', V: 0, R: 1, terminal: 'e-', area: 1e3, reactions: [
      { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 1e-7, alpha: 0.5 },
      { equation: '2 H+ + 2 e- = H2', fixed: { H2: 0 }, k0: 1e-9, alpha: 0.5 },
    ] }],
    grid: { hmin: 0.5e-6, hmax: 0.5e-6 },
  });
  const s = new Device(def).solve(), a = new Device(def).advance(10);
  assert.ok(s.converged && a.converged);
  for (const x of [s, a]) assert.ok(Math.abs(x.terminals.iron.current + x.terminals.iron.V / 1) < 1e-6 * Math.abs(x.terminals.iron.V), 'I = (0 − V)/R');
});

test('the start: a strictly neutral region whose c0 no held level can balance can\'t start a transient, and says so', () => {
  const salt = [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }];
  const sp = { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } };
  const def = {
    species: salt, materials: { w: { epsr: 0, species: sp }, sol: { epsr: 0, species: sp } },
    regions: [{ name: 'gel', material: 'w', length: 10e-6, fixedCharge: -100 * FARADAY, c0: { 'Na+': 50 } }],
    contacts: { left: { V: 0, bath: { c: { 'Na+': 100, 'Cl-': 100 }, reference: 'Cl-', material: 'sol' } }, right: { phi: 'neutral' } },
    grid: { hmin: 0.5e-6, hmax: 0.5e-6 },
  };
  const a = new Device(def).advance(1);
  assert.ok(!a.converged);
  assert.match(a.warnings.join(' '), /the start can't be neutral in gel/);
  assert.ok(new Device(def).solve().converged); // (the steady state takes Na⁺ from its bath)
});

test('the start: a species fed by a contact\'s conductance or exchange link starts at its outside level, no c0 asked for', () => {
  const d = new Device({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ material: 'solid', length: 1e-6, fixedCharge: -100 * FARADAY }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': { type: 'conductance', G: 1e3, offset: 0 } }, phi: 'bulk' }, right: { phi: 'neutral' } },
  });
  assert.ok(Math.abs(d.solution().c['Li+'][0] / 100 - 1) < 1e-6);
});

test('the start: a metal beside a semiconductor isn\'t shifted as if it were an insertion host', () => {
  // It started 0.5 V from the contact holding it, reading −1.9e6 A/m² through a gate that passes none.
  const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19), ND = units.perCm3(1e16);
  const st = new Device({
    species: [{ name: 'e-', z: -1 }, { name: 'h+', z: 1 }],
    materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } }, Au: { conductor: { species: 'e-', conductivity: 4e7 } } },
    regions: [{ material: 'Au', length: 50e-9 }, { material: 'Si', length: 1e-6, fixedCharge: ND * FARADAY }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.01 }, zeroCharge: 0, species: { 'e-': 'blocked' } }],
    contacts: { left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0.5, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' } },
    grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
  }).solution();
  assert.equal(st.contacts.left.current, 0);
  assert.ok(Math.abs(st.V['e-'][0]) < 1e-12, `${st.V['e-'][0]}`);
});

// From round 15's fuzzer: a transient straight from c0's layout (a 1e18 | 1e18 cm⁻³ Si junction,
// its depletion not yet formed) failed at t = 0 with a list of causes that never said to solve
// first; now it does, and from the solved state the same transient runs.
test('a transient that fails at its first step from c0\'s layout says to solve() first, and runs from there', () => {
  const def = build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donors: units.perCm3(1e18) }), layer('Si', 1e-6, { acceptors: units.perCm3(1e18) }), ohmic(0)], grid: { hmin: units.nm(1), hmax: units.nm(20) } });
  const cold = new Device(def).advance(1e-9);
  if (!cold.converged) assert.ok(cold.warnings.some((w) => /solve\(\) first/.test(w)), cold.warnings.join('\n'));
  const d = new Device(def);
  assert.ok(d.solve().converged);
  assert.ok(d.advance(1e-9).converged);
});

test('a call ending a hair before a pulse\'s edge: the next steps across it from a restarted order, and the pulse\'s charge is kept', () => {
  // (Fable round 2's hunter: BDF2's history across the edge lost 1/12 of the first step's charge,
  // 8.33% of a 1 ms pulse, silently; backward Euler kept it.)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2, w = 1e-3;
  const ulpsBelow = (x, n) => {
    const f = new Float64Array([x]), u = new BigUint64Array(f.buffer);
    u[0] -= BigInt(n);
    return f[0];
  };
  for (const tb of [1, 1e5]) {
    for (const tPrev of [tb, ulpsBelow(tb, 2), tb * (1 - 1e-13)]) {
      const dev = new Device({
        species: lib.species,
        materials: lib.materials,
        regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
        contacts: { left: { V: 0, phi: { type: 'capacitive', C }, zeroCharge: 0 }, right: { I: { t: [0, tb, tb, tb + w, tb + w], values: [0, 0, 1, 1, 0] }, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
        grid: { hmin: 2e-6, hmax: 5e-6 },
      });
      dev.solve();
      dev.advance(tb - 0.01);
      dev.advance(tPrev);
      const sol = dev.advance(tb + 2e-3);
      assert.ok(sol.converged);
      const Q = sol.plates.right.charge;
      assert.ok(Math.abs(Q / w - 1) < 2e-5, `tb ${tb}, ${tb - tPrev} s short: ${Q} for ${w}`);
    }
  }
});

test('a step of regular size ending a few ulps before a pulse\'s edge restarts the order too, within one call', () => {
  // (Fable round 2's hunter, on its re-check: outside atBreak's 1e-9 of the step, inside
  // nextBreakpoint()'s 1e-12 of the time; BDF2 stepped across the edge, 1.4% of the charge lost)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2, w = 1e-3, T0 = 1e5;
  const up = (x, n) => {
    const f = new Float64Array([x]), u = new BigUint64Array(f.buffer);
    u[0] += BigInt(n);
    return f[0];
  };
  let tsum = T0;
  for (let k = 0; k < 5; k++) tsum += 1e-3; // where five steps of dtMax land
  for (const tb of [tsum, up(tsum, 1), up(tsum, 3)]) {
    const dev = new Device({
      species: lib.species,
      materials: lib.materials,
      regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
      contacts: { left: { V: 0, phi: { type: 'capacitive', C }, zeroCharge: 0 }, right: { I: { t: [0, tb, tb, tb + w, tb + w], values: [0, 0, 1, 1, 0] }, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
      grid: { hmin: 2e-6, hmax: 5e-6 },
    });
    dev.solve();
    dev.advance(T0);
    const sol = dev.advance(T0 + 8e-3, { dtMax: 1e-3 });
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.plates.right.charge / w - 1) < 1e-6, `tb − tsum ${tb - tsum}: ${sol.plates.right.charge}`);
  }
});

test('dt0 sets a later call\'s first step too, over the size carried from the last call', () => {
  // (Fable round 2's hunter: after a set(), dt0 was ignored, so the first instants of a new
  // drive were stepped over at the carried size)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2;
  const dev = new Device({
    species: lib.species,
    materials: lib.materials,
    regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, phi: { type: 'capacitive', C }, zeroCharge: 0 }, right: { V: 0, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
  dev.solve();
  dev.advance(10);
  dev.set({ contacts: { right: { V: 0.01 } } });
  const sol = dev.advance(11, { dt0: 1e-9 });
  assert.ok(sol.converged);
  const first = sol.trace.t.find((t) => t > 10);
  assert.ok(first - 10 <= 1.01e-9, `first step ${first - 10} s`);
});

test('a trace\'s charge is Σ current × the step ending at each point, across a pulse\'s edges, as the docs say', () => {
  // (Fable round 2's hunter: a point at a jump has the value just before it, so a trapezoid
  // miscounts; this is the rule the docs give instead)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2;
  for (const method of ['be', 'bdf2']) {
    const dev = new Device({
      species: lib.species,
      materials: lib.materials,
      regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
      contacts: { left: { V: 0, phi: { type: 'capacitive', C }, zeroCharge: 0 }, right: { I: { t: [0, 0.01, 0.01, 0.011, 0.011], values: [0, 0, 1, 1, 0] }, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
      grid: { hmin: 2e-6, hmax: 5e-6 },
    });
    dev.solve();
    const sol = dev.advance(0.02, { method });
    const { t, current } = sol.trace;
    let Q = current[0] * t[0];
    for (let k = 1; k < t.length; k++) Q += current[k] * (t[k] - t[k - 1]);
    assert.ok(Math.abs(-Q / sol.plates.right.charge - 1) < 1e-6, `${method}: ${-Q} against ${sol.plates.right.charge}`);
  }
});

test('a floating terminal\'s voltage is in the error norm: an RC whose φ is held beside it charges as e^{−t/RC}, tighter with tol', () => {
  // (Fable round 2's hunter: a capacitive contact behind R, against a bath that holds φ, has its
  // dynamics in the terminal's voltage alone; stepped blind to it, tol 1e-3 and 1e-5 took the
  // same 3 steps from τ to 3τ and erred 1% in the charge)
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2, R = 1, tau = R * C, salt = { 'K+': 100, 'Cl-': 100 };
  const run = (tol) => {
    const dev = new Device(build({ library: [lib], stack: [bath(salt, 'Cl-'), layer('water', 100e-6, { c0: salt }), { V: 0, R, phi: { type: 'capacitive', C }, zeroCharge: 0 }], grid: { hmin: 2e-6, hmax: 5e-6 } }));
    const s0 = dev.solve();
    dev.set({ contacts: { right: { V: 0.01 } } });
    dev.advance(tau, { tol });
    const s = dev.advance(3 * tau, { tol });
    return { q: (s.plates.right.charge - s0.plates.right.charge) / (0.01 * C), steps: s.steps };
  };
  const loose = run(1e-3), tight = run(1e-5), want = -Math.expm1(-3);
  assert.ok(tight.steps > 2 * loose.steps, `steps ${loose.steps} at 1e-3, ${tight.steps} at 1e-5`);
  assert.ok(Math.abs(tight.q - want) < 1e-3, `at 1e-5: ${tight.q} against ${want}`);
  assert.ok(Math.abs(loose.q - want) < 1e-2, `at 1e-3: ${loose.q} against ${want}`);
});

// Fable round 3's convergence hunt: a Si diode solved at 150 K and set() to 500 K failed its
// first step, told it was advancing from c0's layout and to solve() first, which it had. Its
// minority carriers carry over at e⁻⁸⁰ of the new equilibrium; the failure now names the jump,
// and the two ways round it work.
test("advance() right after a far jump in T: the failure names the jump, and smaller jumps or a solve() get there", () => {
  const Si = (T) => semiconductor('Si', { T });
  const diode = (T) => build({
    T,
    library: [Si(T)],
    stack: [ohmic(0), layer('Si', units.um(1), { donors: units.perCm3(1e17) }), layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }), ohmic(0)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    grid: { hmin: 1e-9, hmax: 20e-9 },
  });
  const d = new Device(diode(150));
  assert.ok(d.solve().converged);
  d.set({ T: 500, materials: Si(500).materials });
  const a = d.advance(1e-9), said = a.warnings.join(' | ');
  assert.ok(!a.converged && /jump in temperature \(150 → 500 K\)/.test(said) && !/c0's layout/.test(said), said);
  // In jumps of 50 K, each advanced a little, it gets there; and so does solve() at 500 K.
  const e = new Device(diode(150));
  e.solve();
  for (let T = 200; T <= 500; T += 50) {
    e.set({ T, materials: Si(T).materials });
    assert.ok(e.advance(1e-9).converged, `at ${T} K`);
  }
  const f = new Device(diode(150));
  f.solve();
  f.set({ T: 500, materials: Si(500).materials });
  assert.ok(f.solve().converged && f.advance(1e-9).converged);
});

// Fable round 4's review of src/timeSteps.js: step() with BDF2 across a current source's edge
// kept a history that spans the kink in the stored charge, and lost I·dt/2 of it for good,
// silently (−10% at 2 ms steps); and advance(Infinity) stamped the state t = ∞ in no steps.
test('step() with BDF2 across a current switched on: the order restarts at the edge, and the charge is exact', () => {
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), C = 0.2;
  const mk = () => new Device({
    species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, phi: { type: 'capacitive', C }, zeroCharge: 0 }, right: { I: { t: [0, 0.01, 0.01], values: [0, 0, 1] }, phi: { type: 'capacitive', C }, zeroCharge: 0 } },
    grid: { hmin: 2e-6, hmax: 5e-6 },
  });
  for (const dt of [2e-3, 1e-3]) {
    const d = mk();
    d.solve();
    let s;
    for (let k = 0; k < Math.round(0.02 / dt); k++) s = d.step(dt, { method: 'bdf2' });
    assert.ok(Math.abs(s.plates.right.charge / 1e-2 - 1) < 1e-9, `dt ${dt}: Q ${s.plates.right.charge} vs 1e-2`);
  }
  const d = mk();
  d.solve();
  assert.throws(() => d.advance(Infinity), /solve\(\) the steady state/);
  assert.equal(d.advance(1e-3).converged, true);
});

// Fable round 4 (three reviewers apart): where the error estimate has a floor above tol, error
// rejections had none, so the step fell below the time's ulp, steps were "accepted" that didn't
// move it, and maxSteps ran out at a standstill with currents of 1e87 beside converged: true.
// And advance(1e300) blamed the physics for a step below the time's resolution.
test('advance() whose error estimate has a floor above tol stops there, converged false, and says to loosen tol', () => {
  const ions = [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }];
  const start = () => {
    const d = new Device({
      species: ions, materials: { water: { epsr: 0, species: { 'Na+': { D: 1.3e-9, mu0: -261.9e3 }, 'Cl-': { D: 2e-9, mu0: -131.2e3 } } } },
      regions: [{ material: 'water', length: 1e-3, c0: { 'Na+': 500, 'Cl-': 500 } }],
      contacts: { left: { V: 0, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }, right: { phi: 'neutral' } },
      ports: [{ name: 'dl', region: 0, from: 0.5e-3, V: 0, area: 1e4, capacitance: { C: 0.2 } }],
      grid: { hmin: 10e-6, hmax: 50e-6 },
    });
    d.solve();
    d.set({ ports: { dl: { V: 0.2 } } });
    return d;
  };
  const tight = start().advance(0.5, { tol: 1e-8 }), said = tight.warnings.join(' | ');
  assert.ok(!tight.converged && !tight.done && tight.steps + tight.rejected < 2000, `${tight.converged} ${tight.steps} + ${tight.rejected}`);
  assert.match(said, /below the time's own resolution.*loosen tol/);
  const ok = start().advance(0.5, { tol: 1e-5 });
  assert.ok(ok.converged && ok.done && Math.abs(ok.ports[0].charge / (0.2 / 6) - 1) < 1e-3, `${ok.ports[0].charge}`);
  // A run to 1e300 s: the step it needs is below the time's resolution, and it says so.
  const far = start().advance(1e300);
  assert.ok(!far.converged && /below the time's own resolution/.test(far.warnings.join(' | ')), far.warnings.join(' | '));
});

// Fable round 4's review of src/timeSteps.js: a first step whose full step converged and whose half
// step failed was restored, but left stepped set and solvedV cleared, steering a later steady
// solve's continuation as if a step had been taken.
test('a rejected first step leaves the steady solves\' choices as they were', () => {
  const tox = units.nm(10), Cf = 10, C = 1 / (2 / Cf + tox / (3.9 * EPS0)), R = 10 / C;
  const face = { phi: { type: 'capacitive', C: Cf }, zeroCharge: 0 };
  const d = new Device(build({ T: 300, library: [metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }], stack: [ohmic(0, ['e-']), layer('Al', units.nm(20)), face, layer('SiO2', tox), face, layer('Al', units.nm(20)), ohmic(0, ['e-'])] }));
  d.solve();
  d.set({ contacts: { left: { V: 1, R } } });
  const s = d.solver, before = { stepped: s.stepped, solvedV: s.solvedV && Array.from(s.solvedV) };
  const newton = s.newton;
  let calls = 0;
  s.newton = function (dt, opts) {
    const r = newton.call(this, dt, opts);
    return ++calls === 2 ? { ...r, converged: false } : r; // (the first half step fails)
  };
  d.advance(R * C, { maxSteps: 1 });
  s.newton = newton;
  assert.deepEqual({ stepped: s.stepped, solvedV: s.solvedV && Array.from(s.solvedV) }, before);
});
