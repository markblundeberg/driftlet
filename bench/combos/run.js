// Properties that hold for any device, each comparing two routes to one answer. They run on
// every case of every family in families/ (or the ones named):
//
//   node bench/combos/run.js [family ...] [--n 20] [--only prop,prop] [--case k] [--verbose]
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
//   split      advance(t) in one call and in several ends at the same state, to its tolerance
//   dcZ        the impedance at low frequency is the steady dI/dV at plan.V (unique only)
//   stale      check() refuses a solution of the device before a set()
//
// Each failure line names the family, the case and the property: rerun one with
//   node bench/combos/run.js FAMILY --case K --only PROP --verbose
import { readdirSync } from 'node:fs';
import { Device } from '../../src/index.js';
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
function diff(a, b) {
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
  const gross = Math.max(1e-300, ...Object.values(a.terminals).map((t) => Math.abs(t.current)), ...Object.values(b.terminals ?? {}).map((t) => Math.abs(t.current)));
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

const skip = (why) => `skip: ${why}`;
const converged = (s) => s && s.converged;

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
    const w = d.solve();
    if (!converged(w)) return 'solve() from the relaxed state: not converged';
    const dw = diff(tr, w);
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
    return within(dd, { mu: 1e-6, I: 1e-6, V: 1e-6 }) ? '' : `after ${log.join(', ')}, vs fresh: ${fmt(dd)}`;
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
      if (!within(dd, { mu: 1e-12, I: 1e-12, Iabs: 0, V: 1e-12 })) return `after ${op} threw at linear solve ${at} of ${calls}, then ${k === 1 ? 'advance' : 'solve'}: ${fmt(dd)}`;
    }
    return '';
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
    const dd = diff(one, many);
    // (tol 1e-6 per step, local: a few hundred steps' worth)
    return within(dd, { mu: 1e-3, I: 1e-3, Iabs: 1e-9, V: 1e-4 }) ? '' : `one advance vs pieces: ${fmt(dd)}`;
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
