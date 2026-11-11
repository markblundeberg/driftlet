import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readdirSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join } from 'node:path';

// Every module in src/ loads as the first one imported, as a visitor's test of one module does.
// (SolverError and profileAt lived in solver.js, which the modules they served import back: entered
// anywhere but solver.js or index.js, five modules threw "Cannot access 'X' before initialization";
// Fable's review of the solver's assembly, 10 Nov 2026.)
const src = join(dirname(fileURLToPath(import.meta.url)), '..', 'src');
test('every module in src/ loads when imported first', () => {
  const failed = [];
  for (const f of readdirSync(src).filter((n) => n.endsWith('.js'))) {
    if (f === 'worker.js') continue; // (a Worker's entry: it listens on self)
    const r = spawnSync(process.execPath, ['--input-type=module', '-e', `await import(${JSON.stringify(pathToFileURL(join(src, f)).href)})`], { encoding: 'utf8' });
    if (r.status !== 0) failed.push(`${f}: ${(r.stderr.match(/^\w*Error.*$/m) ?? [r.stderr.trim()])[0]}`);
  }
  assert.deepEqual(failed, []);
});
