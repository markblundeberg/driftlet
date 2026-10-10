import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, metal } from '../src/kit.js';

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
  const qInf = new Device(dev.def).solve().gates.left.charge;
  const lam = Math.sqrt((epsr * EPS0 * RT) / (2 * FARADAY * FARADAY * c));
  const Cs = 1 / (1 / Cg + lam / (epsr * EPS0));
  const sigma = (FARADAY * FARADAY * c * (Dp + Dm)) / RT;
  const tau = ((L / sigma) * Cs) / 2;
  assert.ok(Math.abs(qInf / ((Cs * dV) / 2) - 1) < 1e-3);

  const dt = tau / 20;
  for (let n = 1; n <= 80; n++) {
    const sol = dev.step(dt);
    assert.ok(sol.converged);
    const q = sol.gates.left.charge / qInf;
    // Distributed double-layer charging adds a small correction (~λ_D/L and ~τ_D/τ).
    assert.ok(Math.abs(q - beStep(n, dt, tau)) < 1e-2, `n=${n}: ${q} vs ${beStep(n, dt, tau)}`);
    const gross = FARADAY * 2 * c * L;
    assert.ok(Math.abs(sol.charge + sol.gates.left.charge + sol.gates.right.charge) < 1e-12 * gross);
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
    reference = sol.gates.left.charge;
  }
  return reference;
};

test('BDF2 is second order in time, backward Euler first order', () => {
  const run = (n, method) => {
    const dev = chargingIsland();
    let sol;
    for (let k = 0; k < n; k++) sol = dev.step(T / n, { method });
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `${method}: ${st.species} drift ${st.drift}`);
    return sol.gates.left.charge;
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
    return Math.abs(sol.gates.left.charge / ref - 1);
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
  assert.ok(Math.abs(sol.gates.left.charge / ref - 1) < 2e-2);
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
  for (let n = 1; n <= 40; n++) {
    const s = d.step(dt);
    assert.ok(s.converged && Math.abs(s.interfaces[0].D / qInf - beStep(n, dt, tau)) < 1e-12, `n=${n}: ${s.interfaces[0].D / qInf} vs ${beStep(n, dt, tau)}`);
  }
  const r = start().advance(5 * tau);
  assert.ok(r.steps < 100 && r.rejected === 0, `${r.steps} steps, ${r.rejected} rejected`);
  assert.ok(Math.abs(r.interfaces[0].D / qInf - (1 - Math.exp(-5))) < 1e-3);
});
