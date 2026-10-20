// Fingerprints of every solve in the test suite, for a change meant to change nothing (a module
// moved, a refactor): each test file runs with a preload that hashes the solver's whole state
// (u, uLo, termV, and what the call returned) after every solveSteady, advance, integrate and
// impedance, one line per call. Two runs, before and after, should match line for line. Wall
// time is replaced by a clock that ticks 0.1 ms per reading (about a step's), so a run that
// stops on its budget (budgetMs) stops at the same step every time.
//
//   node tools/fingerprint.js record DIR [--root CHECKOUT] [files…]   (default: this checkout)
//   node tools/fingerprint.js compare DIR_A DIR_B
//
// To compare against an older commit: `git worktree add /tmp/old <commit>`, then record with
// --root /tmp/old, using this copy of the tool. A test that fails under the fixed clock is
// reported; its fingerprints still compare.

import { spawnSync } from 'node:child_process';
import { mkdirSync, readdirSync, readFileSync, existsSync } from 'node:fs';
import { join, resolve, basename } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = fileURLToPath(new URL('.', import.meta.url));
const [cmd, ...args] = process.argv.slice(2);

if (cmd === 'record') {
  const dir = resolve(args[0]), at = args.indexOf('--root'), root = resolve(at >= 0 ? args[at + 1] : join(here, '..'));
  const files = args.slice(1).filter((a, i, all) => a !== '--root' && all[i - 1] !== '--root');
  const tests = files.length > 0 ? files : readdirSync(join(root, 'test')).filter((f) => f.endsWith('.test.js')).sort().map((f) => `test/${f}`);
  // (Each call appends to DIR/<test>.txt: into a folder an earlier run used, the comparison
  // would carry that run's lines too, and can pass on them.)
  if (existsSync(dir) && readdirSync(dir).length > 0) {
    console.log(`record: ${dir} isn't empty; give a new folder`);
    process.exit(2);
  }
  mkdirSync(dir, { recursive: true });
  let failed = 0;
  for (const t of tests) {
    const r = spawnSync(process.execPath, ['--import', join(here, 'fingerprint-hook.js'), '--test', '--test-concurrency=1', t], {
      cwd: root,
      env: { ...process.env, FINGERPRINT_OUT: dir, FINGERPRINT_SRC: join(root, 'src') },
      encoding: 'utf8',
    });
    if (r.status !== 0) {
      failed++;
      console.log(`fails under the fixed clock: ${t}`);
    }
  }
  const lines = readdirSync(dir).reduce((n, f) => n + readFileSync(join(dir, f), 'utf8').split('\n').length - 1, 0);
  console.log(`${tests.length} files, ${lines} fingerprints in ${dir}${failed ? `, ${failed} files failing` : ''}`);
} else if (cmd === 'compare') {
  const [a, b] = args.map((d) => resolve(d)), names = new Set([...readdirSync(a), ...readdirSync(b)]);
  let same = 0, differ = 0;
  for (const f of [...names].sort()) {
    const A = existsSync(join(a, f)) ? readFileSync(join(a, f), 'utf8').split('\n') : [], B = existsSync(join(b, f)) ? readFileSync(join(b, f), 'utf8').split('\n') : [];
    const k = A.findIndex((line, i) => line !== B[i]);
    if (k < 0 && A.length === B.length) same++;
    else {
      differ++;
      console.log(`${basename(f, '.txt')}: first differs at solve ${(k < 0 ? Math.min(A.length, B.length) : k) + 1} (${A.length - 1} and ${B.length - 1} solves)`);
    }
  }
  console.log(`${same} files the same, ${differ} differ`);
  process.exitCode = differ > 0 ? 1 : 0;
} else {
  console.log('usage: node tools/fingerprint.js record DIR [--root CHECKOUT] [files…] | compare DIR_A DIR_B');
  process.exitCode = 2;
}
