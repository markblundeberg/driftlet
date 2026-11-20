// Properties that hold for any device, each comparing two routes to one answer. They run on
// every case of every family in families/ (or the ones named):
//
//   node bench/combos/run.js [family ...] [--n 20] [--only prop,prop] [--case k] [--verbose]
//   npm run combos -- --n 8      (as CI runs it: exits 1 on any failure)
//
//   order      key order and the order of lists that mean nothing (species, reactions, ports)
//              change no answer
//   relax      a long transient at a drive ends where solve() does, warm from it (the steady and
//              the time-stepped equations are one discretisation) and, where the state is
//              unique, where a cold solve does
//   history    a device taken through other drives, transients, an impedance and a parameter
//              changed and changed back ends, at plan.V, where a fresh device does (unique only)
//   rollback   a call that throws partway (a fault injected into the linear solves) leaves the
//              device as it was: what it does next is what a twin that never made the call does,
//              to the last bit
//   noop       a set() that changes nothing physical, mid-transient, changes nothing after
//   split      advance(t) in one call and in several ends at the same state, to its tolerance
//   charge     the charge each terminal passes over a slow ramp after a rest, in one default
//              advance(), is what a run with its steps capped passes, to 3%
//   dcZ        the impedance at low frequency is the steady dI/dV at plan.V (unique only)
//   area       a constant cross-section A scales every terminal current by A, and moves no level
//   stale      check() refuses a solution of the device before a set()
//   mirror     reflect positions, face laws and drives: the reflected device's steady state is the reflection
//   gauge      shift the common voltage reference and compare concentrations, currents and level differences
//   scale      dilate the device and its transport coefficients without changing the dimensionless system
//
// Each failure line names the family, the case and the property: rerun one with
//   node bench/combos/run.js FAMILY --case K --only PROP --verbose
import { readdirSync } from 'node:fs';
import { Device } from '../../src/index.js';
import { FARADAY } from '../../src/constants.js';
import { check } from '../../src/kit.js';
import { caseRng, rng, hash } from './rng.js';

const argv = process.argv.slice(2);
const flag = (name, dflt) => {
  const k = argv.indexOf(`--${name}`);
  return k >= 0 ? argv[k + 1] : dflt;
};
const takesValue = new Set(['--n', '--only', '--case']);
const positional = argv.filter((a, k) => !a.startsWith('--') && !(k > 0 && takesValue.has(argv[k - 1])));
const N = +flag('n', 20), only = flag('only')?.split(','), onlyCase = flag('case'), verbose = argv.includes('--verbose');
const dir = new URL('./families/', import.meta.url);
const familyNames = positional.length ? positional : readdirSync(dir).filter((f) => f.endsWith('.js')).map((f) => f.slice(0, -3));

// --- helpers
const clone = (x) => structuredClone(x);
const isTerminalContact = (p) => p.terminal === 'left' || p.terminal === 'right';
const driveAt = (plan, V) => (isTerminalContact(plan) ? { contacts: { [plan.terminal]: { V } } } : { ports: { [plan.terminal]: { V } } });
const RT = (sol) => 8.314462618 * (sol.T ?? 300);

// The largest differences between two solutions of one device: in any species' level (in RT),
// in φ (in V_T), and in any terminal's current (relative to the largest terminal current) and V.
// (scale: a current to measure the differences against, where the solutions' own are a decayed
// transient's tail: the run's peak)
function diff(a, b, scale = 0) {
  const out = { mu: 0, phi: 0, I: 0, V: 0, where: '' };
  const rt = RT(a);
  for (const name of Object.keys(a.mu ?? {})) {
    const x = a.mu[name], y = b.mu?.[name];
    if (!y || y.length !== x.length) { out.mu = Infinity; out.where = `mu.${name} missing or regridded`; continue; }
    for (let i = 0; i < x.length; i++) {
      if (!Number.isFinite(x[i]) && !Number.isFinite(y[i])) continue;
      const d = Math.abs(x[i] - y[i]) / rt;
      if (!(d <= out.mu)) { out.mu = Number.isNaN(d) ? Infinity : d; out.where = `mu.${name}[${i}]`; }
    }
  }
  const gross = Math.max(1e-300, scale, ...Object.values(a.terminals).map((t) => Math.abs(t.current)), ...Object.values(b.terminals ?? {}).map((t) => Math.abs(t.current)));
  for (const [name, t] of Object.entries(a.terminals)) {
    const u = b.terminals?.[name];
    if (!u) { out.I = Infinity; continue; }
    const dI = Math.abs(t.current - u.current) / gross, dV = Math.abs(t.V - u.V);
    if (!(dI <= out.I)) { out.I = Number.isNaN(dI) ? Infinity : dI; out.Iat = `${name}: ${t.current} vs ${u.current}`; }
    if (!(dV <= out.V)) { out.V = Number.isNaN(dV) ? Infinity : dV; out.Vat = `${name}: ${t.V} vs ${u.V}`; }
  }
  out.gross = gross;
  return out;
}
const fmt = (d) => `levels ${d.mu.toExponential(1)} RT (${d.where}), currents ${d.I.toExponential(1)} of ${d.gross.toExponential(2)} A/m² (${d.Iat ?? ''}), V ${d.V.toExponential(1)} (${d.Vat ?? ''})`;
// Within: levels to muTol RT, currents to ITol of the largest (plus an absolute floor), V to VTol.
const within = (d, { mu = 1e-6, I = 1e-6, Iabs = 1e-10, V = 1e-6 } = {}) => d.mu <= mu && (d.I <= I || d.I * d.gross <= Iabs) && d.V <= V;

// Shuffle every plain object's keys, and reverse the lists whose order means nothing.
function reorder(x, r, path = '') {
  if (Array.isArray(x)) {
    const items = x.map((v, k) => reorder(v, r, `${path}[${k}]`));
    return /(^species|bulkReactions|^ports|reactions)$/.test(path) ? items.reverse() : items;
  }
  if (x && typeof x === 'object' && Object.getPrototypeOf(x) === Object.prototype) {
    const keys = Object.keys(x);
    for (let k = keys.length - 1; k > 0; k--) {
      const j = Math.floor(r.next() * (k + 1));
      [keys[k], keys[j]] = [keys[j], keys[k]];
    }
    return Object.fromEntries(keys.map((k) => [k, reorder(x[k], r, path ? `${path}.${k}` : k)]));
  }
  return x;
}
const reverseOrders = (def) => reorder(def, rng(1));

// The largest terminal current any of these transients' traces passed.
const peak = (...runs) => Math.max(0, ...runs.flatMap((r) => Object.values(r?.trace?.terminals ?? {}).flatMap((t) => Array.from(t.current, Math.abs))));
const skip = (why) => `skip: ${why}`;
const converged = (s) => s && s.converged;

// The families' spatial profiles use device x. A port window instead uses distance from
// its own region's left end. OCV tables also have x, but that x is a filling fraction.
function mapSpatialProfiles(def, map) {
  // A generator can reuse one profile under several keys (e.g. a light source).
  const replace = (object, key) => {
    const p = object?.[key];
    if (p?.x) object[key] = map(p);
  };
  for (const region of def.regions) {
    replace(region, 'fixedCharge');
    for (const name of Object.keys(region.c0 ?? {})) replace(region.c0, name);
    replace(region.grade, 'epsr');
    for (const [name, grade] of Object.entries(region.grade ?? {})) {
      if (name === 'epsr') continue;
      for (const key of Object.keys(grade)) replace(grade, key);
    }
  }
  for (const rx of def.bulkReactions ?? []) for (const name of Object.keys(rx.kf ?? {})) replace(rx.kf, name);
  for (const port of def.ports ?? []) replace(port, 'area');
  replace(def.geometry, 'area');
}

const swapEnds = (name) => (name === 'left' ? 'right' : name === 'right' ? 'left' : name);
const regionOf = (def, port) => (Number.isInteger(port.region) ? def.regions[port.region] : def.regions.find((r) => r.name === port.region));

function mirrorDefinition(def) {
  const length = def.regions.reduce((sum, r) => sum + r.length, 0);
  mapSpatialProfiles(def, (p) => ({ x: p.x.map((x) => length - x).reverse(), values: [...p.values].reverse() }));
  for (const port of def.ports ?? []) {
    const { length } = regionOf(def, port);
    if (port.from !== undefined || port.to !== undefined) [port.from, port.to] = [length - (port.to ?? length), length - (port.from ?? 0)];
    if (Number.isInteger(port.region)) port.region = def.regions.length - 1 - port.region;
  }
  for (const grid of [def.grid, ...def.regions.map((r) => r.grid)]) {
    if (grid) [grid.hLeft, grid.hRight] = [grid.hRight, grid.hLeft];
  }
  for (const region of def.regions) if (region.velocity !== undefined) region.velocity *= -1;
  def.regions.reverse();
  def.interfaces?.reverse();
  for (const face of def.interfaces ?? []) {
    if (!face) continue;
    if (face.dipole !== undefined) face.dipole *= -1;
    if (face.step) face.step.value *= -1;
    for (const gate of Object.values(face.gates ?? {})) {
      for (const rate of [gate.alpha, gate.beta]) {
        rate.midpoint *= -1;
        rate.scale *= -1;
      }
    }
    for (const rx of face.reactions ?? []) {
      if (rx.equation) rx.equation = rx.equation.replace(/\((left|right)\)/g, (_, name) => `(${swapEnds(name)})`);
      else [rx.left, rx.right] = [rx.right, rx.left];
    }
  }
  def.contacts = { left: def.contacts?.right, right: def.contacts?.left };
}

function mapDrive(term, key, map, fallback) {
  const drive = term[key] ?? fallback;
  if (drive === undefined) return;
  term[key] = typeof drive === 'number' ? map(drive) : { ...drive, values: drive.values.map(map) };
}

function shiftDefinition(def, shift) {
  def.contacts ??= {};
  for (const name of ['left', 'right']) {
    def.contacts[name] ??= {};
    const contact = def.contacts[name];
    if (contact.I === undefined) mapDrive(contact, 'V', (v) => v + shift, 0);
  }
  const charged = new Set(def.species.filter((s) => s.z !== 0).map((s) => s.name));
  for (const port of def.ports ?? []) {
    if (port.I !== undefined) continue;
    const electrical = port.capacitance || port.reactions?.length || Object.entries(port.species ?? {}).some(([name, link]) => charged.has(name) && (link.type ?? link) !== 'blocked');
    if (electrical) mapDrive(port, 'V', (v) => v + shift, 0);
  }
}

function scaleFields(object, keys, factor) {
  for (const key of keys) if (object?.[key] !== undefined) object[key] *= factor;
}

function scaleDefinition(def, factor) {
  const squared = factor ** 2;
  const radial = def.geometry?.type ?? def.geometry;
  const areaFactor = radial === 'spherical' ? squared : radial === 'cylindrical' ? factor : 1;
  const currentFactor = factor * areaFactor;
  mapSpatialProfiles(def, (p) => ({ x: p.x.map((x) => x * factor), values: p.values }));
  scaleFields(def.geometry, ['r0'], factor);
  for (const material of Object.values(def.materials)) {
    const models = Array.isArray(material.statistics) ? material.statistics : [material.statistics];
    for (const model of models) if (model?.type === 'debye-huckel') model.epsr ??= material.epsr;
    scaleFields(material, ['epsr'], squared);
    scaleFields(material.conductor, ['conductivity'], squared);
    for (const sp of Object.values(material.species ?? {})) scaleFields(sp, ['D'], squared);
  }
  for (const grid of [def.grid, ...def.regions.map((r) => r.grid)]) scaleFields(grid, ['hmin', 'hmax', 'hLeft', 'hRight'], factor);
  for (const region of def.regions) {
    scaleFields(region, ['length', 'velocity'], factor);
    scaleFields(region, ['mixing'], squared);
    if (region.grade?.epsr) region.grade.epsr.values = region.grade.epsr.values.map((v) => v * squared);
    for (const grade of Object.values(region.grade ?? {})) if (grade.D) grade.D.values = grade.D.values.map((v) => v * squared);
  }
  for (const face of def.interfaces ?? []) {
    if (!face) continue;
    scaleFields(face, ['sheetCharge'], factor);
    scaleFields(face.phi, ['C'], factor);
    scaleFields(face.diffuse, ['epsr'], squared);
    for (const sp of Object.values(face.species ?? {})) scaleFields(sp, ['G', 'P'], factor);
    for (const sp of Object.values(face.surface ?? {})) scaleFields(sp, ['capacity'], factor);
    for (const rx of face.reactions ?? []) {
      scaleFields(rx, ['k0', 'vmax'], factor);
      scaleFields(rx.srh, ['vn', 'vp'], factor);
    }
  }
  for (const contact of Object.values(def.contacts ?? {})) {
    scaleFields(contact.phi, ['C'], factor);
    for (const sp of Object.values(contact.species ?? {})) scaleFields(sp, ['G', 'k', 'v'], factor);
  }
  for (const port of def.ports ?? []) {
    scaleFields(port, ['from', 'to'], factor);
    if (def.materials[regionOf(def, port).material].conductor) {
      for (const sp of Object.values(port.species ?? {})) scaleFields(sp, ['G'], factor);
    }
  }
  for (const term of [...Object.values(def.contacts ?? {}), ...(def.ports ?? [])]) {
    scaleFields(term, ['R'], 1 / currentFactor);
    mapDrive(term, 'I', (i) => i * currentFactor);
  }
  return currentFactor;
}

function matchingGrid(a, b, { mirror = false, lengthFactor = 1 }) {
  const length = a.at(-1);
  return a.length === b.length && a.every((x, i) => {
    const mapped = mirror ? length - b[b.length - 1 - i] : b[i] / lengthFactor;
    return Math.abs(x - mapped) <= 1e-11 * length;
  });
}

function compareSteady(a, b, { mirror = false, shift = 0, currentFactor = 1, currentFloor, neutralPorts }) {
  // Use RT for levels, a relative concentration error, and the existing current/voltage
  // conventions. Current errors are measured against the largest current, with a noise floor.
  let worst = 1, failure = '';
  const compare = (label, x, y, tol) => {
    if (Number.isNaN(x) && Number.isNaN(y)) return;
    const error = Math.abs(x - y) / tol;
    if (error <= worst) return;
    worst = Number.isFinite(error) ? error : Infinity;
    failure = `${label}: original ${x} vs mapped ${y} (absolute tolerance ${tol.toExponential(2)})`;
  };
  const nameInB = (name) => (mirror ? swapEnds(name) : name);
  const z = Object.fromEntries(a.species.map((sp) => [sp.name, sp.z]));
  const readouts = (s) => ({ ...s.contacts, ...Object.fromEntries(s.ports.map((p) => [p.name, p])) });
  const [aRead, bRead] = [readouts(a), readouts(b)];
  const chargeFlux = (t) => Object.entries(t.flux).reduce((sum, [name, flux]) => sum + Math.abs(z[name] * FARADAY * flux), 0);
  for (const [name, t] of Object.entries(a.terminals)) {
    const mapped = nameInB(name), u = b.terminals[mapped];
    if (!neutralPorts.has(name)) compare(`${name} voltage (V)`, t.V, u.V - shift, 1e-6);
    const scaled = u.current / currentFactor;
    // Large opposing carrier currents leave a small difference whose round-off scales with
    // the carriers. Use the same 1e-12 of gross charge flux as the steady current-sum warning.
    const gross = Math.max(chargeFlux(aRead[name]), chargeFlux(bRead[mapped]) / currentFactor);
    const tol = Math.max(currentFloor, 1e-6 * Math.max(Math.abs(t.current), Math.abs(scaled)), 1e-12 * gross);
    compare(`${name} current (A or A/m²)`, t.current, scaled, tol);
  }
  for (const { name, z } of a.species) {
    for (let i = 0; i < a.x.length; i++) {
      const j = mirror ? b.x.length - 1 - i : i;
      compare(`mu.${name}[${i}] (J/mol)`, a.mu[name][i], b.mu[name][j] - z * FARADAY * shift, 1e-6 * RT(a));
      const x = a.c[name][i], y = b.c[name][j];
      compare(`c.${name}[${i}] (mol/m³)`, x, y, 1e-5 * Math.max(1e-12, Math.abs(x), Math.abs(y)));
    }
  }
  return failure && `${failure}; warning counts ${a.warnings.length} vs ${b.warnings.length}`;
}

function compareDefinitions(c, def, plan, transform) {
  const original = new Device(clone(c.def)).set(driveAt(c.plan, c.plan.V));
  const a = original.solve();
  if (!converged(a)) return skip('no steady state at plan.V');
  const mapped = new Device(def).set(driveAt(plan, plan.V));
  if (!matchingGrid(a.x, mapped.model.grid.x, transform)) return skip('the generated grids do not map node for node');
  const b = mapped.solve();
  if (!converged(b)) return `steady convergence: ${Number(a.converged)} vs ${Number(b.converged)}; ${b.warnings.at(-1) ?? 'no warning'}`;
  const { grid, ports } = original.model;
  const neutralPorts = new Set(ports.filter((p) => !p.passes).map((p) => p.name));
  const currentFloor = 1e-10 * Math.max(grid.area[0], grid.area.at(-1));
  return compareSteady(a, b, { ...transform, currentFloor, neutralPorts });
}

// --- the properties: each takes a case and returns '' (holds), 'skip: …', or what failed
const props = {
  order(c) {
    const a = new Device(clone(c.def)).set(driveAt(c.plan, c.plan.V)).solve();
    if (!converged(a)) return skip('no steady state at plan.V');
    const def2 = reorder(clone(c.def), rng(hash(JSON.stringify(c.params)) + 7));
    const b = new Device(def2).set(driveAt(c.plan, c.plan.V)).solve();
    if (!converged(b)) return 'reordered definition: not converged, where the original did';
    const d = diff(a, b);
    return within(d, { mu: 1e-7, I: 1e-7, V: 1e-7 }) ? '' : `reordered: ${fmt(d)}`;
  },

  relax(c) {
    const d = new Device(clone(c.def));
    if (!converged(d.solve())) return skip('no steady state at the drives');
    d.set(driveAt(c.plan, c.plan.V));
    const t1 = 40 * c.plan.slow;
    const tr = d.advance(t1, { tol: 1e-5, budgetMs: 20000, maxSteps: 50000 });
    if (!tr.converged) return skip(`transient failed at t = ${tr.time.toPrecision(3)} s`);
    if (!tr.done) return skip(`transient too slow: t = ${tr.time.toPrecision(3)} of ${t1.toPrecision(3)} s in ${tr.steps} steps`);
    // (settled? a family's slow is an estimate: one was 4× short, its transient still 0.12 RT from
    // the end at 40 × slow and exactly there by 640 ×; another 40 × must move nothing)
    const tr2 = d.advance(2 * t1, { tol: 1e-5, budgetMs: 20000, maxSteps: 50000 });
    if (!tr2.converged || !tr2.done) return skip('transient too slow to check it settled');
    if (!within(diff(tr, tr2, peak(tr)), { mu: 1e-6, I: 1e-6, V: 1e-7 })) return skip('not settled at 40 × plan.slow (slow underestimated)');
    const w = d.solve();
    if (!converged(w)) return 'solve() from the relaxed state: not converged';
    const dw = diff(tr2, w, peak(tr));
    // (40 relaxation times: what's left of the transient is e^-40 of the step)
    if (!within(dw, { mu: 1e-5, I: 1e-5, V: 1e-6 })) return `40 × slow later, the transient vs solve() from there: ${fmt(dw)}`;
    if (!c.plan.unique) return '';
    const cold = new Device(clone(c.def)).set(driveAt(c.plan, c.plan.V)).solve();
    if (!converged(cold)) return skip('relaxed; a cold solve at plan.V failed');
    const dc = diff(w, cold);
    return within(dc, { mu: 1e-6, I: 1e-6, V: 1e-6 }) ? '' : `relaxed then solved vs cold at plan.V: ${fmt(dc)}`;
  },

  history(c, r) {
    if (!c.plan.unique) return skip('not unique');
    const d = new Device(clone(c.def)), log = [];
    if (!converged(d.solve())) return skip('no steady state at the drives');
    // A numeric leaf of a material to change and change back: a diffusivity.
    const leaves = [];
    for (const [m, mat] of Object.entries(c.def.materials ?? {})) for (const [s, sp] of Object.entries(mat.species ?? {})) if (typeof sp.D === 'number' && sp.D > 0) leaves.push([m, s, sp.D]);
    const ops = r.int(3, 6);
    for (let j = 0; j < ops; j++) {
      const op = r.pick(['bias', 'bias', 'transient', 'impedance', 'material', 'grid']);
      log.push(op);
      if (op === 'bias') {
        const s = d.set(driveAt(c.plan, c.plan.V + r.u(-0.2, 0.2))).solve();
        if (!converged(s)) return skip(`history (${log.join(', ')}): a bias failed`);
      } else if (op === 'transient') {
        d.set(driveAt(c.plan, c.plan.V + r.u(-0.1, 0.1)));
        d.advance(d.solution().time + r.log(-2, 0) * c.plan.slow, { budgetMs: 5000 });
      } else if (op === 'impedance') {
        try {
          d.impedance([1, 1e3], { terminal: c.plan.terminal });
        } catch (e) {
          log.push(`(refused: ${e.message.split('\n')[0].slice(0, 60)})`);
        }
      } else if (op === 'material' && leaves.length) {
        const [m, s, D] = r.pick(leaves), f = r.pick([0.5, 2, 10]);
        d.set({ materials: { [m]: { species: { [s]: { D: D * f } } } } });
        if (r.chance(0.5) && !converged(d.solve())) return skip(`history (${log.join(', ')}): solve with D × ${f} failed`);
        d.set({ materials: { [m]: { species: { [s]: { D } } } } });
      } else if (op === 'grid' && c.def.grid && typeof c.def.grid.hmax === 'number') {
        d.set({ grid: { hmax: c.def.grid.hmax / 2 } });
        if (r.chance(0.5)) d.solve();
        d.set({ grid: { hmax: c.def.grid.hmax } });
      }
    }
    const a = d.set(driveAt(c.plan, c.plan.V)).solve();
    if (!converged(a)) return `after ${log.join(', ')}: not converged at plan.V, where a fresh device may be`;
    const b = new Device(clone(c.def)).set(driveAt(c.plan, c.plan.V)).solve();
    if (!converged(b)) return skip('a fresh device at plan.V failed');
    const dd = diff(a, b);
    // (1e-9 A/m² of currents read as small differences of large terms: a lit thin-film transistor's
    // 8e-10 read 1e-10 after a history, its levels the same to 2e-15 RT)
    return within(dd, { mu: 1e-6, I: 1e-6, Iabs: 1e-9, V: 1e-6 }) ? '' : `after ${log.join(', ')}, vs fresh: ${fmt(dd)}`;
  },

  rollback(c, r) {
    const make = () => {
      const d = new Device(clone(c.def));
      return converged(d.solve()) ? d : null;
    };
    const D = make(), T = make();
    if (!D || !T) return skip('no steady state at the drives');
    const op = r.pick(['solve', 'advance', 'step', 'impedance']);
    const V1 = c.plan.V, t = () => D.solution().time;
    const call = (dev) => {
      dev.set(driveAt(c.plan, V1));
      if (op === 'solve') return dev.solve();
      if (op === 'advance') return dev.advance(dev.solution().time + c.plan.slow, { budgetMs: 10000 });
      if (op === 'step') return dev.step(c.plan.slow / 10);
      return dev.impedance([10], { terminal: c.plan.terminal });
    };
    // (the drive set on both before the call: set() itself isn't the call under test)
    T.set(driveAt(c.plan, V1));
    // Count the twin's linear solves for the call, to inject a fault somewhere inside it.
    let calls = 0;
    {
      const P = make();
      const s = P.solver, orig = s._solveLinear;
      s._solveLinear = function (...a) { calls++; return orig.apply(this, a); };
      try { call(P); } catch { /* (refused: no fault to inject into) */ }
    }
    if (calls < 2) return skip(`${op}: ${calls} linear solves`);
    const at = r.int(1, calls - 1);
    const s = D.solver, orig = s._solveLinear;
    let n = 0;
    s._solveLinear = function (...a) {
      if (++n === at) throw new Error('injected fault');
      return orig.apply(this, a);
    };
    let threw = false;
    try { call(D); } catch (e) { threw = e.message === 'injected fault'; if (!threw) return skip(`${op} refused: ${e.message.slice(0, 80)}`); }
    s._solveLinear = orig;
    if (!threw) return skip(`${op}: the fault at linear solve ${at} of ${calls} wasn't reached`);
    // Then both do the same: a transient at another drive, and a solve.
    const next = (dev) => {
      dev.set(driveAt(c.plan, V1 + c.plan.dV));
      const a = dev.advance(dev.solution().time + c.plan.slow, { maxSteps: 3000 });
      return [dev.solution(), a, dev.solve()];
    };
    let x, y;
    try { y = next(T); } catch (e) { return skip(`the twin's next calls threw: ${e.message.slice(0, 80)}`); }
    try { x = next(D); } catch (e) { return `after ${op} threw at linear solve ${at} of ${calls}: the next calls threw (${e.message.split('\n')[0].slice(0, 100)}), where the twin's didn't`; }
    if (x[0].time !== y[0].time) return `after ${op} threw at linear solve ${at} of ${calls}: time ${x[0].time} vs the twin's ${y[0].time}`;
    if (x[1].steps !== y[1].steps) return `after ${op} threw at linear solve ${at} of ${calls}: the next advance took ${x[1].steps} steps, the twin's ${y[1].steps}`;
    for (let k = 1; k < 3; k++) {
      const dd = diff(x[k], y[k]);
      // (to round-off, not the bit: what Newton learned on the call that threw is kept, LEARNED in
      // src/solver.js, so a later solve may take another path to the same state)
      if (!within(dd, { mu: 1e-9, I: 1e-9, Iabs: 1e-12, V: 1e-12 })) return `after ${op} threw at linear solve ${at} of ${calls}, then ${k === 1 ? 'advance' : 'solve'}: ${fmt(dd)}`;
    }
    return '';
  },

  noop(c, r) {
    // Mid-transient, a set() that changes nothing physical: the same drive again, a material's D
    // given as it is, a drive moved and moved back before any step. The device goes on as a twin
    // that was never patched does, to the last bit.
    const make = () => {
      const d = new Device(clone(c.def));
      if (!converged(d.solve())) return null;
      d.set(driveAt(c.plan, c.plan.V + c.plan.dV));
      d.advance(c.plan.slow / 3, { maxSteps: 3000 });
      return d;
    };
    const D = make(), T = make();
    if (!D) return skip('no steady state at the drives');
    const kind = r.pick(['same drive', 'same D', 'there and back']);
    if (kind === 'same drive') D.set(driveAt(c.plan, c.plan.V + c.plan.dV));
    else if (kind === 'same D') {
      const leaves = [];
      for (const [m, mat] of Object.entries(c.def.materials ?? {})) for (const [sp, v] of Object.entries(mat.species ?? {})) if (typeof v.D === 'number') leaves.push([m, sp, v.D]);
      if (!leaves.length) return skip('no D to give again');
      const [m, sp, Dv] = r.pick(leaves);
      D.set({ materials: { [m]: { species: { [sp]: { D: Dv } } } } });
    } else {
      D.set(driveAt(c.plan, c.plan.V - c.plan.dV));
      D.set(driveAt(c.plan, c.plan.V + c.plan.dV));
    }
    const [x0, y0] = [D.solution(), T.solution()];
    const d0 = diff(x0, y0);
    if (!within(d0, { mu: 0, I: 0, Iabs: 0, V: 0 })) return `${kind}: solution() right after: ${fmt(d0)}`;
    const t1 = T.solution().time + c.plan.slow;
    // (there and back, both at tol 1e-6: the restart's first-order step is then far below the
    // allowance, which at the default tol a gated film's tail exceeded by itself, 1.07e-3 RT, its
    // twin and itself converging together as tol tightened; Astra's ninth pass)
    const opts = kind === 'there and back' ? { maxSteps: 30000, tol: 1e-6 } : { maxSteps: 3000 };
    const [x, y] = [D.advance(t1, opts), T.advance(t1, opts)];
    if (kind !== 'there and back' && x.steps !== y.steps) return `${kind}: the next advance took ${x.steps} steps, the twin's ${y.steps}`;
    const d1 = diff(x, y, peak(x, y));
    // (a drive moved restarts the time steps' order, as any discontinuity does, so there and
    // back costs a first-order step: equal to the transient's tolerance, 1e-3 thermal units per step
    // by default, not to the bit)
    // (and its currents then are a tolerance-level tail, compared through the levels alone)
    const tol = kind === 'there and back' ? { mu: 1e-3, I: Infinity, V: 1e-5 } : { mu: 0, I: 0, Iabs: 0, V: 0 };
    return within(d1, tol) ? '' : `${kind}: after the next advance: ${fmt(d1)}`;
  },

  split(c, r) {
    const run = (pieces) => {
      const d = new Device(clone(c.def));
      if (!converged(d.solve())) return null;
      d.set(driveAt(c.plan, c.plan.V + c.plan.dV));
      const t1 = 3 * c.plan.slow;
      let s;
      for (let j = 1; j <= pieces; j++) {
        s = d.advance((t1 * j) / pieces, { tol: 1e-6, budgetMs: 10000, maxSteps: 30000 });
        if (!s.converged || !s.done) return { failed: `${pieces} pieces: stopped at t = ${s.time.toPrecision(3)}` };
      }
      return s;
    };
    const one = run(1);
    if (!one) return skip('no steady state at the drives');
    if (one.failed) return skip(one.failed);
    const many = run(r.int(3, 9));
    if (many.failed) return `one call reached the end, pieces didn't: ${many.failed}`;
    const dd = diff(one, many, peak(one, many));
    // (tol 1e-6 per step, local: a few hundred steps' worth)
    return within(dd, { mu: 1e-3, I: 1e-3, Iabs: 1e-9, V: 1e-4 }) ? '' : `one advance vs pieces: ${fmt(dd)}`;
  },

  charge(c) {
    // A rest long enough that the step size carried out of it is large, then a slow ramp of the
    // drive: the charge each terminal passes, in one default advance(), against the same run with
    // its steps capped. (A flow-through path's charge was out of the error control: a CV after a
    // rest read 52% high at its vertex; round 7, Sol's re-read of the time integrator.)
    if (!isTerminalContact(c.plan)) return skip('drive is a port');
    const run = (opts) => {
      const d = new Device(clone(c.def));
      if (!converged(d.solve())) return null;
      d.set(driveAt(c.plan, c.plan.V)); // (resting at plan.V: the ramp starts there)
      d.advance(100 * c.plan.slow, { budgetMs: 5000 });
      const t0 = d.solution().time, ramp = 20 * c.plan.slow;
      if (!(t0 >= 100 * c.plan.slow * (1 - 1e-9))) return { failed: 'rest too slow' };
      d.set(driveAt(c.plan, { t: [t0, t0 + ramp], values: [0, c.plan.dV].map((v, k) => (k ? c.plan.V + v : c.plan.V)) }));
      // (the running totals the trace reads, as the solver holds them before the ramp)
      const names = d.solver.terms.map((t) => t.name), before = Float64Array.from(d.solver.termCharge);
      const a = d.advance(t0 + ramp, { budgetMs: 20000, maxSteps: 40000, ...opts(ramp) });
      if (!a.converged || !a.done) return { failed: `stopped at t = ${a.time.toPrecision(3)}` };
      const q = {}, moved = {};
      names.forEach((n, k) => {
        const tr = a.trace.terminals[n].charge;
        q[n] = tr.at(-1) - before[k];
        moved[n] = tr.reduce((sum, v, j) => sum + Math.abs(v - (j ? tr[j - 1] : before[k])), 0);
      });
      return { q, moved, a };
    };
    const one = run(() => ({}));
    if (!one) return skip('no steady state at the drives');
    if (one.failed) return skip(one.failed);
    const ref = run((ramp) => ({ dtMax: ramp / 400, tol: 1e-6 }));
    if (ref.failed) return skip(`reference ${ref.failed}`);
    let worst = 0, at = '';
    // (against the charge each terminal moved either way, not the net: a terminal passing charge
    // both ways over the ramp leaves a net remainder, and its error read 3.2% and 4.5% of that
    // where it was 0.6% and 0.7% of what passed, each confirmed by a separate ledger; Astra's ninth
    // pass. A flow-through charge outside the error control was 52% of it, as of the net.)
    const gross = Math.max(1e-300, ...Object.values(ref.moved));
    for (const [n, q] of Object.entries(ref.q)) {
      const e = Math.abs(one.q[n] - q) / gross;
      if (e > worst) (worst = e), (at = `${n}: ${one.q[n].toExponential(4)} vs ${q.toExponential(4)} C/m² (${one.a.steps} steps vs ${ref.a.steps})`);
    }
    // (3%: a diode's stored charge at the default tol came 1.8% off, falling under 1% at 1e-4, as
    // the potentials' tolerance allows; the flow-through faults were 5 to 90%. And no finer than
    // 1e-6 of the run's peak current over the ramp: a membrane's 2e-11 C/m² at a current near zero
    // is noise.)
    // (a closed cell peaking at 3e-13 A/m², a full cell at open circuit at 1.6e-20: nothing passes)
    if (peak(ref.a) < 1e-9) return skip(`nothing passes (peak ${peak(ref.a).toExponential(1)} A/m²)`);
    const floor = Math.max(1e-12, 1e-6 * peak(ref.a) * 20 * c.plan.slow);
    return worst <= 3e-2 || worst * gross < floor ? '' : `charge passed over a slow ramp after a rest, default advance() vs steps capped: ${(100 * worst).toFixed(1)}% (${at})`;
  },

  dcZ(c) {
    if (!c.plan.unique) return skip('not unique');
    const d = new Device(clone(c.def)).set(driveAt(c.plan, c.plan.V));
    const s = d.solve();
    if (!converged(s)) return skip('no steady state at plan.V');
    const term = c.plan.terminal, dV = 1e-4;
    const I = (V) => {
      const e = new Device(clone(c.def)).set(driveAt(c.plan, V)).solve();
      return converged(e) ? e.terminals[term].current : NaN;
    };
    const G = (I(c.plan.V + dV) - I(c.plan.V - dV)) / (2 * dV);
    if (!Number.isFinite(G)) return skip('no steady solves either side');
    const f = Math.min(1e-6, 1e-4 / c.plan.slow);
    let z;
    try {
      z = d.impedance([f], { terminal: term });
    } catch (e) {
      return skip(`impedance refused: ${e.message.split('\n')[0].slice(0, 80)}`);
    }
    if (!(z.relError[0] < 1e-3)) return skip(`Z unresolved at ${f.toPrecision(2)} Hz (relError ${z.relError[0]})`);
    const re = z.Z.re[0], im = z.Z.im[0], m = re * re + im * im, Yre = re / m, Yim = -im / m;
    if (Math.abs(Yim) > 1e-2 * Math.abs(Yre)) return skip(`still reactive at ${f.toPrecision(2)} Hz`);
    return Math.abs(Yre - G) <= 1e-3 * Math.abs(G) + 1e-12 ? '' : `Re Y(${f.toPrecision(2)} Hz) = ${Yre.toExponential(4)}, steady dI/dV = ${G.toExponential(4)} S/m²`;
  },

  area(c, r) {
    // (a series R or an I drive is per area in a slab and a total with a cross-section: scaled
    // with it; a port's electrode area is its own, so devices with ports are left out, as are
    // devices with a geometry already)
    if (c.def.geometry) return skip('has a geometry');
    if ((c.def.regions ?? []).some((g) => g.velocity) || Object.values(c.def.materials ?? {}).some((m) => m.velocity)) return skip('has flow (planar only)');
    if ((c.def.ports ?? []).length) return skip('has ports');
    const a = new Device(clone(c.def)).set(driveAt(c.plan, c.plan.V)).solve();
    if (!converged(a)) return skip('no steady state at plan.V');
    const A = r.pick([1e-6, 3.7e-3, 42]), L = c.def.regions.reduce((t, g) => t + g.length, 0);
    const def2 = clone(c.def);
    def2.geometry = { area: { x: [0, L], values: [A, A] } };
    for (const side of ['left', 'right']) {
      const k = def2.contacts?.[side];
      if (!k) continue;
      if (typeof k.R === 'number') k.R /= A;
      if (typeof k.I === 'number') k.I *= A;
      else if (k.I && Array.isArray(k.I.values)) k.I.values = k.I.values.map((v) => v * A);
    }
    const b = new Device(def2).set(driveAt(c.plan, c.plan.V)).solve();
    if (!converged(b)) return `with a cross-section of ${A} m²: not converged, where the slab did`;
    const scaled = { ...b, terminals: Object.fromEntries(Object.entries(b.terminals).map(([n, t]) => [n, { ...t, current: t.current / A }])) };
    const dd = diff(a, scaled);
    return within(dd, { mu: 1e-7, I: 1e-7, V: 1e-7 }) ? '' : `a constant cross-section of ${A} m², currents divided by it: ${fmt(dd)}`;
  },

  mirror(c) {
    const def = clone(c.def);
    const geometry = def.geometry?.type ?? def.geometry;
    if (geometry === 'spherical' || geometry === 'cylindrical') return skip('radial geometry needs a centre on the other side');
    if (def.interfaces?.some((f) => f?.sheetCharge)) return skip('the reflected sheet charge would have to sit on the left plate');
    mirrorDefinition(def);
    return compareDefinitions(c, def, { ...c.plan, terminal: swapEnds(c.plan.terminal) }, { mirror: true });
  },

  gauge(c, r) {
    // A common voltage shift moves φ by the same amount and μ_i by z_i F times that
    // amount. Charge-balanced affinities, capacitance drops and neutral levels stay put.
    const shift = r.pick([-0.125, 0.125]);
    const def = clone(c.def);
    shiftDefinition(def, shift);
    return compareDefinitions(c, def, { ...c.plan, V: c.plan.V + shift }, { shift });
  },

  scale(c, r) {
    // At x' = s x, keep t and the concentrations. D and ε scale by s², so the bulk
    // flux and displacement scale by s and their divergences are unchanged. Face flux
    // and capacity must scale by s too. A distributed port's laws per volume stay put.
    const factor = r.pick([0.5, 2]);
    const def = clone(c.def);
    const currentFactor = scaleDefinition(def, factor);
    return compareDefinitions(c, def, c.plan, { lengthFactor: factor, currentFactor });
  },

  stale(c) {
    const d = new Device(clone(c.def));
    const s = d.solve();
    if (!converged(s)) return skip('no steady state at the drives');
    const leaves = [];
    for (const [m, mat] of Object.entries(c.def.materials ?? {})) for (const [sp, v] of Object.entries(mat.species ?? {})) if (typeof v.D === 'number' && v.D > 0) leaves.push([m, sp, v.D]);
    if (!leaves.length) return skip('no D to change');
    const [m, sp, D] = leaves[0];
    d.set({ materials: { [m]: { species: { [sp]: { D: 3 * D } } } } });
    let res;
    try {
      res = check(d, s, { refine: false });
    } catch {
      return ''; // (refused loudly: fine)
    }
    return res.ok ? `check() passed a solution from before set(D × 3)` : '';
  },
};

// --- run
const which = only ?? Object.keys(props);
for (const p of which) if (!props[p]) throw new Error(`no property '${p}' (known: ${Object.keys(props).join(', ')})`);
const tally = {}, failures = [];
const t0 = performance.now();
for (const name of familyNames) {
  const fam = (await import(new URL(`${name}.js`, dir))).default;
  const ks = onlyCase !== undefined ? [+onlyCase] : Array.from({ length: N }, (_, k) => k);
  for (const k of ks) {
    let c;
    try {
      c = fam.make(caseRng(name, k));
    } catch (e) {
      failures.push(`${name} ${k} · make: ${e.message}`);
      continue;
    }
    if (c.plan.transientOnly) continue;
    if (verbose) console.log(JSON.stringify({ params: c.params, plan: c.plan }));
    for (const p of which) {
      const key = `${name}: ${p}`;
      tally[key] ??= { ok: 0, skip: 0, fail: 0 };
      let res;
      const tp = performance.now();
      try {
        res = props[p](c, caseRng(`${name}/${p}`, k));
      } catch (e) {
        res = `threw: ${String(e.message).split('\n')[0].slice(0, 200)}`;
      }
      if (verbose) console.log(`  ${p}: ${res || 'ok'} (${((performance.now() - tp) / 1000).toFixed(1)} s)`);
      if (res === '') tally[key].ok++;
      else if (res.startsWith('skip: ')) tally[key].skip++;
      else tally[key].fail++, failures.push(`${name} ${k} · ${p}: ${res}`);
    }
  }
}
const width = Math.max(...Object.keys(tally).map((s) => s.length));
for (const [key, v] of Object.entries(tally)) console.log(`${key.padEnd(width)}  ok ${String(v.ok).padStart(3)}  skip ${String(v.skip).padStart(3)}  fail ${String(v.fail).padStart(3)}`);
console.log(`${((performance.now() - t0) / 1000).toFixed(1)} s`);
if (failures.length) console.log(`\nFailures:\n  ${failures.join('\n  ')}`);
process.exitCode = failures.length ? 1 : 0; // (CI runs a few cases of every family: any failure is new)
