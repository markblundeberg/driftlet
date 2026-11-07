// The solver's state between calls, made explicit (STATE in src/solver.js): every field a run
// writes is state (snapshotted, restored by a failed call or a rollback), or scratch (written
// before it's read, within each call), or else structure, built with the solver and never written
// after. And the behaviour that rests on it: a call that fails or throws, a read, and a rollback
// each leave a device doing exactly (bit for bit) what a twin that never made the call does.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units } from '../src/index.js';
import { STATE } from '../src/solver.js';
import { aqueous, metal, build, layer, ohmic, semiconductor, check } from '../src/kit.js';

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
  'chargePivots', 'accepted', 'sysRest', 'diffTerms', 'lostCandidate', 'lostMode', 'uPrev', 'uPrevLo', 'termVPrev', 'borderCols', 'islandRow', 'linRow',
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

test('every field a run writes is in STATE or is scratch; STATE names real fields', () => {
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
      for (const [k, v] of after) if (before.get(k) !== v && !STATE.includes(k) && !SCRATCH.has(k)) stray.add(`${name}: ${k}`);
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
    for (const [k, v] of after) if (before.get(k) !== v && !STATE.includes(k) && !SCRATCH.has(k)) stray.add(`dark cell: ${k}`);
    before = after;
  }
  assert.deepEqual([...stray], [], 'fields written between calls that are neither in STATE (src/solver.js) nor scratch (here): which is each?');
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
