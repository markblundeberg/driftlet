// node bench/combos/try.js FAMILY [N] [k]: build each case cold, solve at its drives and at plan.V,
// and report how many converge and pass check(). A family is worth fuzzing with when nearly all do.
import { Device } from '../../src/index.js';
import { check } from '../../src/kit.js';
import { caseRng } from './rng.js';
const [name, N = '20', only] = process.argv.slice(2);
const fam = (await import(`./families/${name}.js`)).default;
const ks = only !== undefined ? [+only] : Array.from({ length: +N }, (_, k) => k);
const tally = { built: 0, cold: 0, checked: 0, biased: 0 }, bad = [];
const t0 = performance.now();
for (const k of ks) {
  let c;
  try { c = fam.make(caseRng(name, k)); } catch (e) { bad.push(`${k} make: ${e.message.split('\n')[0]}`); continue; }
  if (only !== undefined) console.log(JSON.stringify({ params: c.params, plan: c.plan, def: c.def }, null, 1));
  try {
    const d = new Device(c.def); tally.built++;
    const s = d.solve();
    if (!s.converged) { bad.push(`${k} cold: not converged ${s.warnings?.[0] ?? ''}`.slice(0, 200)); continue; }
    tally.cold++;
    const ch = check(d, s, { refine: false });
    if (ch.ok) tally.checked++; else bad.push(`${k} check: ${ch.items.filter((i) => i.ok === false).map((i) => i.name + ': ' + i.summary).join('; ')}`.slice(0, 240));
    const p = c.plan;
    d.set(p.terminal === 'left' || p.terminal === 'right' ? { contacts: { [p.terminal]: { V: p.V } } } : { ports: { [p.terminal]: { V: p.V } } });
    const b = d.solve();
    if (b.converged) tally.biased++; else bad.push(`${k} biased: not converged`);
  } catch (e) { bad.push(`${k} threw: ${e.message.split('\n')[0]}`.slice(0, 200)); }
}
console.log(name, JSON.stringify(tally), `of ${ks.length}`, `${((performance.now() - t0) / 1000).toFixed(1)} s`);
for (const b of bad) console.log('  ' + b);
