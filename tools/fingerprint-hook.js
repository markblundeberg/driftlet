// The preload of tools/fingerprint.js: a fixed clock, and a hash of the solver's state after each
// public solve, appended to FINGERPRINT_OUT/<test file>.txt.
import { appendFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { basename } from 'node:path';
import { pathToFileURL } from 'node:url';

let now = 0;
Object.defineProperty(globalThis, 'performance', { value: { now: () => (now += 0.1) }, configurable: true, writable: true });

const { Solver } = await import(pathToFileURL(`${process.env.FINGERPRINT_SRC}/solver.js`).href);
const out = `${process.env.FINGERPRINT_OUT}/${basename(process.argv[1] ?? 'repl')}.txt`;
const hash = (s, r) => {
  const h = createHash('sha1');
  for (const a of [s.u, s.uLo, s.termV]) if (a) h.update(new Uint8Array(a.buffer, a.byteOffset, a.byteLength));
  if (r && typeof r === 'object') h.update(JSON.stringify(r, (k, v) => (ArrayBuffer.isView(v) ? Array.from(v) : v)) ?? '');
  return h.digest('hex').slice(0, 16);
};
for (const m of ['solveSteady', 'advance', 'integrate', 'impedance']) {
  const orig = Solver.prototype[m];
  Solver.prototype[m] = function (...a) {
    let r;
    try {
      r = orig.apply(this, a);
    } catch (e) {
      appendFileSync(out, `${m} throws ${String(e.message).slice(0, 60)} ${hash(this)}\n`);
      throw e;
    }
    appendFileSync(out, `${m} ${hash(this, r)}\n`);
    return r;
  };
}
