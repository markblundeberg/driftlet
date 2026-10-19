// A sweep for steady states that keep the wrong amounts: small devices combining the features at
// the edges of the steady system's conservation laws (floating and capacitive electrodes, two of
// them, surfaces, immobile species, face reactions, capacitive contacts), each solved directly
// and run long as a transient from the same start. Where the transient has settled (it moves no
// further over 100× as long), the two must agree.
//   node tools/lawsweep.js [--from k] [--to k] [--only i,j,…] [--show k] [--out file.json]
// Prints one line per device that disagrees (or fails), and a tally. --show k prints device k's
// definition and both answers.
import { writeFileSync } from 'node:fs';
import { Device } from '../src/index.js';
import { families, total, choice, device } from './lawsweep-devices.js';

const args = process.argv.slice(2), opt = (name, d) => (args.includes(name) ? args[args.indexOf(name) + 1] : d);

// What a device keeps: each species' amount, each surface species' sites, and each floating
// terminal's voltage.
function measures(d, s) {
  const v = d.model.grid.vol, m = {};
  for (const sp of d.model.species) {
    const c = s.c[sp.name];
    if (!c) continue;
    m[`n ${sp.name}`] = c.reduce((t, x, g) => t + (Number.isFinite(x) ? v[g] * x : 0), 0);
  }
  s.ports.forEach((p) => {
    for (const [name, th] of Object.entries(p.coverage ?? {})) m[`θ ${name}`] = p.x.reduce((t, x, j) => t + v[s.x.indexOf(x)] * p.area[j] * th[j], 0);
  });
  for (const side of ['left', 'right']) if (s.contacts[side].D) m[`D ${side}`] = s.contacts[side].D;
  for (const k of d.solver.floating) {
    const name = k < 2 ? ['left', 'right'][k] : d.model.ports[k - 2].name;
    m[`V ${name}`] = s.terminals[name].V;
  }
  return m;
}
// (Voltages to a mV; charges to a µC/m²; amounts against themselves, floored at 1e-10 of the
// device's whole inventory, below which they're round-off.)
const differ = (a, b) => {
  let worst = { what: null, rel: 0 };
  const inventory = Object.keys(a).filter((k) => k.startsWith('n ')).reduce((t, k) => t + Math.abs(a[k]), 0);
  for (const k of Object.keys(a)) {
    const scale = k.startsWith('V ') ? 1e-3 : k.startsWith('D ') ? Math.max(Math.abs(a[k]), Math.abs(b[k]), 1e-6) : Math.max(Math.abs(a[k]), Math.abs(b[k]), 1e-10 * inventory);
    const rel = Math.abs(a[k] - b[k]) / scale;
    if (rel > worst.rel) worst = { what: k, rel, a: a[k], b: b[k] };
  }
  return worst;
};

const from = +opt('--from', 0), to = Math.min(+opt('--to', total), total), show = opt('--show');
const only = opt('--only')?.split(',').map(Number) ?? (show !== undefined ? [+show] : Array.from({ length: Math.max(0, to - from) }, (_, j) => from + j));
const tally = { ok: 0, skipped: 0, invalid: 0, transientFailed: 0, unsettled: 0, solveFailed: 0, differ: 0 }, found = [];
for (const idx of only) {
  const c = choice(idx), def = device(c);
  if (!def) {
    tally.skipped++;
    continue;
  }
  const tag = `${idx} ${families[c.family].keys.map((k) => c[k]).join(' ')}`;
  try {
    new Device(def);
  } catch (e) {
    tally.invalid++;
    if (show !== undefined) console.log(tag, 'invalid:', e.message);
    continue;
  }
  try {
    const d = new Device(def), s = d.solve();
    const e = new Device(def), late = e.advance(1e6, { tol: 1e-9 });
    const ms = measures(d, s), ml = measures(e, late);
    if (show !== undefined) console.log(JSON.stringify(def, null, 1), '\nsolve', ms, '\nadvance', ml);
    if (!late.converged || !late.done) {
      tally.transientFailed++; // (nothing to compare with)
      continue;
    }
    if (!s.converged) {
      tally.solveFailed++;
      console.log(`${tag}: solve() did not converge`);
      found.push({ idx, c, kind: 'solveFailed' });
      continue;
    }
    const w = differ(ms, ml);
    if (w.rel < 1e-5) {
      tally.ok++;
      continue;
    }
    // Settled? 100× as long again moves it no further.
    const last = e.advance(1e8, { tol: 1e-9 }), later = measures(e, last), drift = differ(ml, later);
    if (show !== undefined) console.log('transient from 1e6 to 1e8 s:', last.converged, drift);
    if (!last.converged || !last.done || drift.rel > 1e-6) {
      tally.unsettled++;
      continue;
    }
    tally.differ++;
    console.log(`${tag}: ${w.what} ${w.a.toPrecision(6)} solved, ${w.b.toPrecision(6)} transient (${w.rel.toExponential(1)})`);
    found.push({ idx, c, kind: 'differ', worst: w });
  } catch (err) {
    tally.solveFailed++;
    console.log(`${tag}: threw ${err.message.slice(0, 120)}`);
    found.push({ idx, c, kind: 'threw', message: err.message });
  }
}
console.log(`of ${only.length} (${total} in all):`, JSON.stringify(tally));
if (opt('--out')) writeFileSync(opt('--out'), JSON.stringify(found, null, 1));
