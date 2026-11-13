// The solver's state between calls, made explicit (STATE in src/solver.js): every field a run
// writes is state (snapshotted, restored by a failed call or a rollback), or scratch (written
// before it's read, within each call), or else structure, built with the solver and never written
// after. And the behaviour that rests on it: a call that fails or throws, a read, and a rollback
// each leave a device doing exactly (bit for bit) what a twin that never made the call does.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units } from '../src/index.js';
import { CARRY, LEARNED, STATE } from '../src/solver.js';
import { AgAgCl, aqueous, bath, metal, build, layer, ohmic, semiconductor, check } from '../src/kit.js';

// --- four small devices: a port along a bar, an electrode's face reaction against a bath, a
// capacitive port, a diode
const devices = {
  port: () => ({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * 96485.33212 }],
    contacts: { left: { V: 0.01, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { phi: 'neutral' } },
    ports: [{ name: 'drain', region: 'bar', V: 0, terminal: 'Li+', species: { 'Li+': { type: 'conductance', G: 1e9 } } }],
    grid: { minCells: 40 },
  }),
  electrode: () => {
    const aq = aqueous(['Fe3+', 'Fe2+', 'K+', 'Cl-'], { epsr: 0 }), au = metal('Au');
    return {
      species: [...aq.species, { name: 'e-', z: -1 }],
      materials: { ...aq.materials, ...au.materials },
      regions: [{ name: 'Au', material: 'Au', length: 1e-6 }, { name: 'sol', material: 'water', length: 10e-6 }],
      interfaces: [{ reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 1e-3, alpha: 0.5 }] }],
      contacts: {
        left: { V: 0.02, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
        right: { bath: { c: { 'Fe3+': 10, 'Fe2+': 10, 'K+': 1000, 'Cl-': 1050 }, reference: 'Cl-' } },
      },
      grid: { hmax: 0.5e-6 },
    };
  },
  // (a capacitive port's charge, which a step's bookkeeping aliases: portQOld = portQStart)
  capacitor: () => ({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * 96485.33212 }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'dl', region: 'bar', V: 0, area: 1e6, capacitance: { C: 0.2 } }],
    grid: { hmax: 0.25e-6 },
  }),
  diode: () =>
    build({
      library: [semiconductor('Si')],
      stack: [ohmic(), layer('Si', 1e-6, { donors: units.perCm3(1e16) }), layer('Si', 1e-6, { acceptors: units.perCm3(1e16) }), ohmic({ V: 0.1 })],
      grid: { hmax: 2e-8 },
    }),
};

// A dark Si cell driven by a current past its reverse saturation: no steady state, so its solve
// fails, after the floating continuation's held solves.
const darkCell = () => {
  const NA = units.perCm3(1e16);
  return build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e19) }), layer('Si', units.um(60), { acceptors: NA }), { ...ohmic(0), V: undefined, I: -50 }],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1 / (60e-6 * NA) } }],
    grid: { hmin: units.nm(1), hmax: units.um(1) },
  });
};

// --- canonical dumps of the solver's fields
function canon(v, seen = new Set(), depth = 0) {
  if (v === null || v === undefined || typeof v !== 'object') return typeof v === 'number' && Object.is(v, -0) ? '-0' : typeof v === 'function' ? 'fn' : String(v);
  if (ArrayBuffer.isView(v)) return `[${Array.from(v, (x) => (Object.is(x, -0) ? '-0' : String(x))).join()}]`;
  if (seen.has(v) || depth > 8) return '<…>';
  seen.add(v);
  const out = Array.isArray(v) ? `[${v.map((x) => canon(x, seen, depth + 1)).join()}]`
    : v instanceof Map ? `Map{${[...v].map(([k, x]) => `${canon(k, seen, depth + 1)}:${canon(x, seen, depth + 1)}`).join()}}`
    : v instanceof Set ? `Set{${[...v].map((x) => canon(x, seen, depth + 1)).join()}}`
    : `{${Object.keys(v).map((k) => `${k}:${canon(v[k], seen, depth + 1)}`).join()}}`;
  seen.delete(v);
  return out;
}
const fields = (s) => new Map(Object.keys(s).filter((k) => k !== 'model').map((k) => [k, canon(s[k])]));
const stateOf = (s) => STATE.map((k) => `${k}=${canon(s[k])}`);
// What a solution says, exactly.
const seen = (sol) => canon({ converged: sol.converged, time: sol.time, warnings: sol.warnings, c: sol.c, phi: sol.phi, terminals: sol.terminals });
const attempt = (fn) => {
  try {
    return fn();
  } catch (e) {
    return { threw: e.message };
  }
};

// A run of what users do: solves, a drive changed and stepped, adaptive time, a rebuild (a
// diffusivity), a solve that fails, a check refined on a finer grid, a current drive.
const work = (side) => [
  (d) => d.solve(),
  (d) => d.set({ contacts: { [side]: { V: d.def.contacts[side].V + 0.01 } } }).solution(),
  (d) => d.step(1e-6),
  (d) => d.step(1e-6, { method: 'bdf2' }),
  (d) => d.advance(d.solver.time + 1e-4),
  (d) => d.solve({ maxSteps: 1 }),
  (d) => (check(d, d.solution(), { refine: true }), d.solution()),
  (d) => d.set({ materials: Object.fromEntries(Object.entries(d.def.materials).slice(0, 1).map(([m, mat]) => [m, { epsr: (mat.epsr ?? 0) + 1e-3 }])) }).solve(),
  (d) => d.set({ contacts: { [side]: { I: 0 } } }).solve(),
  (d) => d.advance(d.solver.time + 1e-4),
];
const sideOf = { port: 'left', electrode: 'left', capacitor: 'left', diode: 'right' };

// Scratch: fields a call writes and reads within itself, each written before it's read.
const SCRATCH = new Set([
  // assembly and Newton's work: the residual, the Jacobian and its factors, the update, the fluxes
  // and charges the last assembly read
  'sys', 'res', 'delta', 'deltaV', 'rowDiv', 'sg', 'dA', 'dB', 'jL', 'jR', 'dWr', 'aWr', 'dWork', 'termB', 'termC', 'termDI', 'termDX', 'termRes', 'termI',
  'contactFlux', 'contactD', 'portFlux', 'portQ', 'qg', 'qgSlope', 'c', 'th', 'th0', 'gam', 'gate', 'ex', 'K', 'dK', 'zeta', 'scratch', 'lin', 'dtNow', 'transformed',
  'chargePivots', 'sysRest', 'diffTerms', 'lostCandidate', 'uPrev', 'uPrevLo', 'termVPrev', 'borderCols', 'islandRow', 'linRow',
  // a step's own: where it lands, its source reading, the laws read as changes this step
  'landing', 'sourceBefore', 'combining', 'lawForm', 'steady', 'sourceOverride',
  // the laws a steady solve captures as it begins (captureLaws), and its stages' bookkeeping
  'lawRefs', 'chargeRefs', 'pins', 'pinnedIslands', 'islandsOn', 'steadyAt', 'carried', 'generationScale', 'dimmer',
  // derived from the drives: the laws and constraints, rebuilt by redrive() (_laws)
  'constraints', 'islands', 'localConstraints', 'affinities', 'carriers', 'contactLaw', 'contactQ', 'floating', 'drivenKey',
  'stretches', 'flatStretches', // (each stretch's conserved flag)
  // (and not the solver's: the model's terminals, whose drives set() replaces, put back by a
  // rollback from the checkpoint's)
  'terms',
  // what a call's solution reads of how it went (diagnostics): set by the call, read just after it
  'violent', 'unreached', 'charging', 'drifted', 'settling', 'modelError', 'conditioning', 'islandImbalance', 'islandsFailed', 'islandsWorse',
]);

test('every field a run writes is in STATE, LEARNED or scratch; STATE names real fields', () => {
  const stray = new Set();
  for (const [name, make] of Object.entries(devices)) {
    const d = new Device(make());
    const solvers = [d.solver];
    let before = fields(d.solver);
    for (const op of work(sideOf[name])) {
      const s0 = d.solver;
      attempt(() => op(d));
      if (d.solver !== s0) {
        solvers.push(d.solver);
        before = fields(new Device(d.def).solver); // (a rebuild: against a solver fresh from its definition)
      }
      const after = fields(d.solver);
      for (const [k, v] of after) if (before.get(k) !== v && !STATE.includes(k) && !SCRATCH.has(k) && !LEARNED.includes(k)) stray.add(`${name}: ${k}`);
      before = after;
    }
    for (const s of solvers) for (const k of STATE) assert.ok(k in s, `STATE names ${k}, which the solver doesn't have`);
  }
  // (and a steady solve that fails deep in its continuations)
  const d = new Device(darkCell());
  let before = fields(d.solver);
  for (const op of [(x) => x.advance(1e-9), (x) => x.solve(), (x) => x.set({ contacts: { right: { V: -20 } } }).solve()]) {
    attempt(() => op(d));
    const after = fields(d.solver);
    for (const [k, v] of after) if (before.get(k) !== v && !STATE.includes(k) && !SCRATCH.has(k) && !LEARNED.includes(k)) stray.add(`dark cell: ${k}`);
    before = after;
  }
  assert.deepEqual([...stray], [], 'fields written between calls that are in none of STATE, LEARNED (src/solver.js) or scratch (here): which is each?');
});

// The twin: the same device put through the same calls, less the one under test.
const run = (make, ops) => {
  const d = new Device(make());
  for (const op of ops) attempt(() => op(d));
  return d;
};
const after = (d, side) => [d.step(1e-6), d.step(1e-6, { method: 'bdf2' }), d.solve(), d.set({ contacts: { [side]: { V: 0.03 } } }).advance(d.solver.time + 1e-4)].map(seen);

test('a rollback undoes a change of drive, a rebuild and the steps since: the device then does what its twin does', () => {
  // (a drive changed in place, then a rebuild, left the checkpoint's model driven at the new value)
  for (const [name, make] of Object.entries(devices)) {
    const side = sideOf[name], prep = work(side).slice(0, 4);
    const d = run(make, prep), twin = run(make, prep);
    const cp = d._checkpoint(), was = stateOf(d.solver);
    for (const op of work(side).slice(1)) attempt(() => op(d));
    d._rollback(cp);
    assert.deepEqual(stateOf(d.solver), was, `${name}: the state after the rollback`);
    assert.deepEqual(d.model.terminals.map((t) => t.drive), twin.model.terminals.map((t) => t.drive), `${name}: the drives after the rollback`);
    assert.deepEqual(after(d, side), after(twin, side), `${name}: what the device does after the rollback`);
  }
});

test('a solve that fails, a call that throws and a read each leave the state as it was', () => {
  const calls = {
    'solve({ maxSteps: 1 })': (d) => d.solve({ maxSteps: 1 }),
    'advance(Infinity)': (d) => d.advance(Infinity),
    'step(-1)': (d) => d.step(-1),
    'a set() that fails': (d) => d.set({ regions: [{ material: 'no such material', length: 1e-6 }] }),
    'check(d) refined': (d) => check(d, d.solution(), { refine: true }),
  };
  for (const [name, make] of Object.entries(devices)) {
    const side = sideOf[name];
    for (const prep of [work(side).slice(0, 1), work(side).slice(0, 4)]) {
      for (const [what, call] of Object.entries(calls)) {
        const d = run(make, prep), twin = run(make, prep), was = stateOf(d.solver);
        const r = attempt(() => call(d));
        if (what.startsWith('solve') && r.converged) continue; // (a solve that converged anyway)
        assert.deepEqual(stateOf(d.solver), was, `${name}, ${prep.length} calls in: ${what}`);
        assert.deepEqual(after(d, side), after(twin, side), `${name}, ${prep.length} calls in: what the device does after ${what}`);
      }
    }
  }
});

test('a steady solve that fails deep in its continuations (no steady state) leaves the state as it was', () => {
  for (const prep of [[], [(d) => d.advance(1e-9)]]) {
    const d = run(darkCell, prep), twin = run(darkCell, prep), was = stateOf(d.solver);
    assert.equal(d.solve().converged, false);
    assert.deepEqual(stateOf(d.solver), was, `${prep.length} calls in`);
    assert.deepEqual(after(d, 'right'), after(twin, 'right'), `${prep.length} calls in: what the device does after the failed solve`);
  }
});

// set(), from Fable's code-smell review of the front door (10 Nov 2026): three silent wrong answers
// where set() carries the state by hand.
import { triangle, ramp } from '../src/kit.js';
const pn = () => build({ library: [semiconductor('Si')], bulkReactions: [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1e-7, tauP: 1e-7 } } }],
  stack: [ohmic(0), layer('Si', 2e-6, { donors: units.perCm3(1e16) }), layer('Si', 2e-6, { acceptors: units.perCm3(1e16) }), ohmic(0.5)] });

test('set: a rebuild carries each terminal’s running charge, as it carries the clock', () => {
  // (zeroed, a charge read across a slider move that changed anything but a drive was short by
  // all that had passed before it)
  const dev = new Device(pn());
  const I = dev.solve().terminals.right.current;
  dev.advance(1e-6);
  const D = dev.def.materials.Si.species['e-'].D;
  dev.set({ materials: { Si: { species: { 'e-': { D: D * (1 + 1e-7) } } } } }); // (a rebuild, same unknowns)
  const q = dev.advance(2e-6).trace.terminals.right.charge.at(-1);
  assert.ok(Math.abs(q / (I * 2e-6) - 1) < 1e-4, `${q} vs ${I * 2e-6}`);
  dev.set({ grid: { hmin: 1e-9, hmax: 5e-8 } }); // (and a regrid)
  const after = dev.advance(3e-6), q3 = after.trace.terminals.right.charge.at(-1);
  // (what passed since is the new grid's current times 1 µs, and a little more as the interpolated
  // state settles on it, 0.6%: reset, the total would read a third of this)
  assert.ok(Math.abs((q3 - q) / (after.terminals.right.current * 1e-6) - 1) < 1e-2, `${q3 - q} vs ${after.terminals.right.current * 1e-6}`);
});

test('set: a waveform replaces the old one whole', () => {
  // (merged key by key, a ramp set over a triangle kept the triangle's repeat)
  const dev = new Device(pn());
  dev.set({ contacts: { right: { V: triangle({ from: 0, to: 0.2, period: 1 }) } } });
  dev.set({ contacts: { right: { V: ramp({ from: 0, to: 0.2, duration: 1 }) } } });
  assert.equal(dev.def.contacts.right.V.repeat, undefined);
  assert.deepEqual(dev.def.contacts.right.V, ramp({ from: 0, to: 0.2, duration: 1 }));
});

test('set: the state carries over only where the species are the same, in the same order', () => {
  // (tested by count, a reordered list swapped two species' amounts: A 5, B 1 where c0 had A 1, B 5)
  const def = {
    species: [{ name: 'A', z: 0 }, { name: 'B', z: 0 }],
    materials: { m: { epsr: 0, species: { A: { D: 1e-9, mu0: 0, cRef: 1000 }, B: { D: 1e-9, mu0: 0, cRef: 1000 } } } },
    regions: [{ material: 'm', length: 1e-6, c0: { A: 1, B: 5 } }],
    contacts: { left: { species: {} }, right: { species: {} } },
  };
  const dev = new Device(def);
  dev.solve();
  dev.set({ species: [{ name: 'B', z: 0 }, { name: 'A', z: 0 }] });
  const sol = dev.solve(), g = sol.x.length >> 1;
  assert.ok(Math.abs(sol.c.A[g] - 1) < 1e-9 && Math.abs(sol.c.B[g] - 5) < 1e-9, `A ${sol.c.A[g]}, B ${sol.c.B[g]}`);
});

test('set: every STATE field has a policy for a rebuild, and the kept ones are kept', () => {
  // (set() carries the state by hand: a field added to STATE after it was written, a terminal's
  // running charge, was dropped. CARRY in src/solver.js names what set() does with each.)
  assert.deepEqual(Object.keys(CARRY).sort(), [...STATE].sort());
  for (const [k, v] of Object.entries(CARRY)) assert.ok(['kept', 'rebased', 'restarted', 'set'].includes(v), `${k}: ${v}`);
  // On three states that between them give every kept field a value of its own: a diode
  // stepped; the same solved (solvedV); an electrode with an adsorbate (surfaceRef), floating
  // (a floating terminal's termV). (On the first alone, removing solvedV, termV's floating carry or
  // surfaceRef from set() passed; Fable's re-check.)
  const salt = { 'K+': 100, 'Cl-': 100 };
  const sorbing = () => build({
    library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')],
    stack: [ohmic(0.6, ['e-']), layer('Pt', 1e-6), { surface: { Clads: { mu0: 58e3, capacity: 2e-5, theta0: 0.1 } }, reactions: [{ equation: 'Cl- = Clads + e-', k0: 1e-3, alpha: 0.5 }] }, layer('water', 10e-6, { c0: salt }), bath(salt, 'Cl-')],
    grid: { hmin: 1e-8, hmax: 0.5e-6 },
  });
  const cases = [
    [pn, (d) => d.solve() && d.advance(1e-6), { materials: { Si: { species: { 'e-': { D: 0.0036 * (1 + 1e-7) } } } } }],
    [pn, (d) => d.solve(), { materials: { Si: { species: { 'e-': { D: 0.0036 * (1 + 1e-7) } } } } }],
    [sorbing, (d) => d.solve() && d.set({ contacts: { left: { I: 0 } } }).advance(1e-6), { materials: { water: { species: { 'K+': { D: 1.957e-9 * (1 + 1e-7) } } } } }],
  ];
  for (const [make, run, patch] of cases) {
    const dev = new Device(make());
    run(dev);
    const old = dev.solver, before = Object.fromEntries(STATE.map((k) => [k, canon(old[k])]));
    dev.set(patch); // (a rebuild, same unknowns)
    assert.notEqual(dev.solver, old);
    for (const k of STATE) {
      if (CARRY[k] !== 'kept') continue;
      assert.equal(canon(dev.solver[k]), before[k], `${k} is marked kept, and set() changed it`);
    }
    // (a closed stretch keeps its amount, its intake folded in; an open one takes the state's)
    old.referenceAmounts.forEach((a, j) => assert.equal(dev.solver.referenceAmounts[j], old.stretches[j].connected ? dev.solver.amount(dev.solver.stretches[j]) : a + old.boundaryIntake[j]));
  }
});

// --- set() across a new grid, and what it carries: from Fable's re-check of CARRY (12 Nov 2026).
// A silver cell, Ag | AgNO₃ | Ag, its NO₃⁻ closed (c0·L = 2e-4 mol/m²): i_lim tanh(V/4V_T) is
// its current, 143.743 A/m² at 0.05 V and 298.15 K, 138.989 at 310 K.
const silver = (grid, extra = {}) => {
  const links = { terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
  return {
    species: [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }],
    materials: { water: { epsr: 78.5, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
    regions: [{ material: 'water', length: 20e-6, c0: { 'NO3-': 10 } }],
    contacts: { left: { V: 0, ...links }, right: { V: 0.05, ...links } },
    grid,
    ...extra,
  };
};
const coarse = { hmin: 0.2e-9, hmax: 400e-9, ratio: 1.15 }, fine = { hmin: 0.2e-9, hmax: 100e-9, ratio: 1.15 };
const no3 = (sol) => sol.conservation.find((q) => q.species === 'NO3-').amount;

test('set: a new grid with a new T or μ° in the same set() carries concentrations, and keeps a closed amount', () => {
  // (Interpolated as potentials over RT, set({ T: 310, grid }) lost 82% of the NO₃⁻, and a μ°
  // change 99%: 24.98 and 1.506 A/m², converged and unwarned.)
  for (const [patch, I] of [[{ T: 310 }, 138.989], [{ materials: { water: { species: { 'NO3-': { mu0: -100e3 } } } } }, 143.743]]) {
    const d = new Device(silver(coarse));
    d.solve();
    const s = d.set({ ...patch, grid: fine }).solve();
    assert.ok(s.converged && Math.abs(no3(s) / 2e-4 - 1) < 1e-9 && Math.abs(s.terminals.right.current / I - 1) < 1e-5, `${JSON.stringify(patch)}: ${no3(s)}, ${s.terminals.right.current}`);
  }
});

test('set: a closed amount is the same on every grid, through many regrids', () => {
  // (Read off the interpolated state, it ratcheted down at each, −0.47% over 60.)
  const d = new Device(silver(coarse));
  d.solve();
  let s;
  for (let k = 0; k < 12; k++) s = d.set({ grid: { ...coarse, hmax: (k % 2 ? 400e-9 : 150e-9) * (1 + k / 50) } }).solve();
  const fresh = new Device(silver(d.def.grid)).solve();
  assert.ok(Math.abs(no3(s) / no3(fresh) - 1) < 1e-12 && Math.abs(s.terminals.right.current / fresh.terminals.right.current - 1) < 1e-9, `${no3(s)} vs ${no3(fresh)}`);
});

test('set: a regridded state follows each old cell’s profile, and a transient keeps its closed amount', () => {
  // The electrode stepped 0.02 → 0.04 V, refined 2× at t = 10 µs: interpolated linearly in η, Fe²⁺
  // beside the face came out at the geometric mean (0.041 between 0.0034 and 0.50), and the next
  // current read 22.24 A/m² where the same grid gives 69.06.
  const run = (regrid) => {
    const d = new Device(devices.electrode());
    d.solve();
    d.set({ contacts: { left: { V: 0.04 } } }).advance(1e-5);
    if (regrid) d.set({ grid: { hmax: 0.25e-6 } });
    return d.advance(1.1e-5).terminals.left.current;
  };
  const same = run(false), moved = run(true);
  assert.ok(Math.abs(moved / same - 1) < 1e-3, `${moved} vs ${same}`); // (69.062 vs 69.060)
  const d = new Device(silver(coarse));
  d.solve();
  d.set({ contacts: { right: { V: 0.3 } } }).advance(2e-3);
  d.set({ grid: fine });
  const s = d.advance(4e-3);
  assert.ok(Math.abs(no3(s) / 2e-4 - 1) < 1e-12, `${no3(s)}`); // (−0.019% before: the transient keeps what it's given)
});

test('set: a regridded steady state is solved there before another rebuild, and before step()', () => {
  // (Its flag went with the old solver: a second set() then stepped the interpolated state,
  // 23.7 A/m² where 68.8 is right; step() never settled it, 27.0.)
  const steady = new Device(devices.electrode()).solve().terminals.left.current;
  const twice = new Device(devices.electrode());
  twice.solve();
  twice.set({ grid: { hmax: 0.25e-6 } }).set({ materials: { water: { species: { 'K+': { D: 1.957e-9 * (1 + 1e-9) } } } } });
  const stepped = new Device(devices.electrode());
  stepped.solve();
  stepped.set({ grid: { hmax: 0.25e-6 } });
  for (const I of [twice.advance(1e-5).terminals.left.current, stepped.step(1e-5).terminals.left.current]) assert.ok(Math.abs(I / steady - 1) < 1e-3, `${I} vs ${steady}`);
});

test('set: ports reordered keep their own charges; regions reordered restart, and say so', () => {
  const def = devices.capacitor();
  def.ports = [{ name: 'a', region: 'bar', from: 1e-6, to: 3e-6, V: 0, area: 1e6, capacitance: { C: 0.2 } }, { name: 'b', region: 'bar', from: 6e-6, to: 8e-6, V: 0.05, area: 1e6, capacitance: { C: 0.05 } }];
  const passed = (swap) => {
    const d = new Device(structuredClone(def));
    d.solve();
    d.set({ contacts: { left: { V: 0.02 } } }).advance(1e-5);
    const q0 = Object.fromEntries(d.solver.terms.map((t, k) => [t.name, d.solver.termCharge[k]]));
    d.set({ ports: swap ? [d.def.ports[1], d.def.ports[0]] : [{ ...d.def.ports[0], area: 1e6 * (1 + 1e-12) }, d.def.ports[1]] }).advance(2e-5);
    return Object.fromEntries(d.solver.terms.map((t, k) => [t.name, d.solver.termCharge[k] - q0[t.name]]));
  };
  const kept = passed(false), swapped = passed(true);
  // (by index, 'b' read +0.0233 C/m² and 'a' −0.0238, where −1.0e-6 and −4.5e-4 are right)
  for (const p of ['a', 'b']) assert.ok(Math.abs(swapped[p] / kept[p] - 1) < 1e-6, `${p}: ${swapped[p]} vs ${kept[p]}`);
  // Two closed compartments of A⁻, 'poor' (10) and 'rich' (50), swapped: carried by place,
  // 'rich' read 10 and 'poor' 50, unwarned.
  const sp = (D) => ({ D, mu0: 0, cRef: 1000 }), sealed = { V: 0, terminal: 'Cl-', species: { 'K+': 'blocked', 'Cl-': 'equilibrium', 'A-': 'blocked' }, phi: 'bulk' };
  const poor = { name: 'poor', material: 'gel', length: 5e-6, c0: { 'K+': 110, 'Cl-': 100, 'A-': 10 } }, rich = { name: 'rich', material: 'gel', length: 5e-6, c0: { 'K+': 150, 'Cl-': 100, 'A-': 50 } };
  const gel = (regions) => ({
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }, { name: 'A-', z: -1 }],
    materials: { gel: { epsr: 0, species: { 'K+': sp(1e-9), 'Cl-': sp(1e-9), 'A-': sp(1e-10) } } },
    regions, interfaces: [{ species: { 'A-': 'blocked' }, phi: 'neutral' }], contacts: { left: sealed, right: sealed }, grid: { minCells: 20 },
  });
  const d = new Device(gel([poor, rich]));
  d.solve();
  const s = d.set({ regions: [rich, poor] }).solve(), A = (r) => s.c['A-'][s.region.indexOf(r) + 3];
  assert.ok(Math.abs(A(0) - 50) < 1e-6 && Math.abs(A(1) - 10) < 1e-6 && s.warnings.some((w) => w.includes('the regions were reordered')), `${A(0)}, ${A(1)}: ${s.warnings}`);
});

test('set: a restart is said through a further rebuild; a c0 given before anything ran is the start', () => {
  const d = new Device(devices.port());
  d.solve();
  d.set({ species: [{ name: 'Li+', z: 1, cRef: 1000 }, { name: 'X', z: 0, cRef: 1000 }], materials: { solid: { species: { X: { D: 1e-11, mu0: 0 } } } }, regions: [{ ...d.def.regions[0], c0: { X: 1 } }] });
  d.set({ materials: { solid: { species: { 'Li+': { D: 2e-11 } } } } });
  assert.ok(d.solution().warnings.some((w) => w.includes("the state couldn't carry over (the species changed)")));
  // (after a bare solution(), the old c0 was carried as if it were a state: 143.74 A/m², not 287.49)
  const e = new Device(silver(fine));
  e.solution();
  const s = e.set({ regions: [{ ...e.def.regions[0], c0: { 'NO3-': 20 } }] }).solve();
  assert.ok(Math.abs(s.terminals.right.current / 287.486 - 1) < 1e-5, `${s.terminals.right.current}`);
});

test('set: a contact the kit wrote, a link or a φ law given with its type, replaces the old one whole', () => {
  // (Merged, bath(c, 'Cl-') over an Ag/AgCl bath kept its 1.138 V offset: −6863.6 A/m² where a
  // fresh device reads −340.06; a velocity link over an equilibrium one kept its offset.)
  const c = { 'K+': 100, 'Cl-': 100 }, aq = aqueous(['K+', 'Cl-']);
  const cell = (left) => build({ library: [aq], stack: [left, layer('water', 10e-6), bath(c, 'Cl-', 0.01)], grid: { hmax: 0.5e-6 } });
  const d = new Device(cell(bath(c, AgAgCl)));
  d.solve();
  const I = d.set({ contacts: { left: bath(c, 'Cl-') } }).solve().terminals.left.current, fresh = new Device(cell(bath(c, 'Cl-'))).solve().terminals.left.current;
  assert.ok(Math.abs(I / fresh - 1) < 1e-9 && d.def.contacts.left.bath.offset === undefined, `${I} vs ${fresh}`);
  const e = new Device(devices.diode());
  e.set({ contacts: { right: { species: { 'e-': { type: 'equilibrium', offset: 0.2 } }, phi: { type: 'capacitive', C: 0.1 }, zeroCharge: 0 } } });
  e.set({ contacts: { right: { species: { 'e-': { type: 'velocity', v: 1e5 } }, phi: { type: 'pinned' } } } });
  assert.deepEqual(e.def.contacts.right.species['e-'], { type: 'velocity', v: 1e5 });
  assert.deepEqual(e.def.contacts.right.phi, { type: 'pinned' });
});
