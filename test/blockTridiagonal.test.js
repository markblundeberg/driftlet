import { test } from 'node:test';
import assert from 'node:assert/strict';
import { BlockTridiagonal } from '../src/blockTridiagonal.js';

// Small seeded PRNG (mulberry32) so failures are reproducible.
function rng(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296 - 0.5;
  };
}

// Fill a random, block-diagonally dominant system; return the RHS.
function randomSystem(sys, rand) {
  const { n, m } = sys;
  const mm = m * m;
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < mm; k++) {
      if (i > 0) sys.A[i * mm + k] = rand();
      if (i < n - 1) sys.C[i * mm + k] = rand();
      sys.B[i * mm + k] = rand();
    }
    for (let r = 0; r < m; r++) sys.B[i * mm + r * m + r] += 2 * m + 1;
  }
  return Float64Array.from({ length: n * m }, rand);
}

// Dense reference: expand to a full matrix and solve by Gaussian elimination with pivoting.
function denseSolve(sys, d) {
  const { n, m, A, B, C } = sys;
  const N = n * m, mm = m * m;
  const M = Array.from({ length: N }, () => new Float64Array(N + 1));
  for (let i = 0; i < n; i++) {
    for (let r = 0; r < m; r++) {
      const row = M[i * m + r];
      for (let c = 0; c < m; c++) {
        row[i * m + c] = B[i * mm + r * m + c];
        if (i > 0) row[(i - 1) * m + c] = A[i * mm + r * m + c];
        if (i < n - 1) row[(i + 1) * m + c] = C[i * mm + r * m + c];
      }
      row[N] = d[i * m + r];
    }
  }
  for (let k = 0; k < N; k++) {
    let p = k;
    for (let r = k + 1; r < N; r++) if (Math.abs(M[r][k]) > Math.abs(M[p][k])) p = r;
    [M[k], M[p]] = [M[p], M[k]];
    for (let r = k + 1; r < N; r++) {
      const f = M[r][k] / M[k][k];
      for (let c = k; c <= N; c++) M[r][c] -= f * M[k][c];
    }
  }
  const x = new Float64Array(N);
  for (let r = N - 1; r >= 0; r--) {
    let s = M[r][N];
    for (let c = r + 1; c < N; c++) s -= M[r][c] * x[c];
    x[r] = s / M[r][r];
  }
  return x;
}

const maxAbsDiff = (a, b) => a.reduce((mx, v, k) => Math.max(mx, Math.abs(v - b[k])), 0);

test('matches a dense solve for many sizes', () => {
  const rand = rng(1);
  for (const n of [1, 2, 3, 7, 40]) {
    for (const m of [1, 2, 3, 5, 8]) {
      const sys = new BlockTridiagonal(n, m);
      const d = randomSystem(sys, rand);
      sys.factor();
      const x = sys.solve(d);
      const ref = denseSolve(sys, d);
      assert.ok(maxAbsDiff(x, ref) < 1e-12, `n=${n} m=${m}: ${maxAbsDiff(x, ref)}`);
      const resid = maxAbsDiff(sys.multiply(x), d);
      assert.ok(resid < 1e-12, `n=${n} m=${m} residual ${resid}`);
    }
  }
});

test('handles zero diagonal entries by pivoting within blocks', () => {
  // Permute the rows within every block row of a well-posed system. The solution is
  // unchanged, but the diagonal blocks now have zeros (or tiny values) on their diagonals.
  const rand = rng(2);
  const n = 20, m = 5, mm = m * m;
  const sys = new BlockTridiagonal(n, m);
  const d = randomSystem(sys, rand);
  const ref = denseSolve(sys, d);
  const perm = [3, 0, 4, 1, 2];
  for (const M of [sys.A, sys.B, sys.C]) {
    for (let i = 0; i < n; i++) {
      const blk = M.slice(i * mm, (i + 1) * mm);
      for (let r = 0; r < m; r++) for (let c = 0; c < m; c++) M[i * mm + r * m + c] = blk[perm[r] * m + c];
    }
  }
  const dp = new Float64Array(d.length);
  for (let i = 0; i < n; i++) for (let r = 0; r < m; r++) dp[i * m + r] = d[i * m + perm[r]];
  // Make some diagonal entries exactly zero.
  for (let i = 0; i < n; i++) sys.B[i * mm] = 0;
  const ref2 = denseSolve(sys, dp);
  sys.factor();
  const x = sys.solve(dp);
  assert.ok(maxAbsDiff(x, ref2) < 1e-11);
  assert.ok(maxAbsDiff(sys.multiply(x), dp) < 1e-11);
  assert.ok(ref.length === x.length);
});

test('solve may overwrite its right-hand side in place', () => {
  const rand = rng(3);
  const sys = new BlockTridiagonal(10, 3);
  const d = randomSystem(sys, rand);
  sys.factor();
  const ref = sys.solve(d);
  const buf = Float64Array.from(d);
  sys.solve(buf, buf);
  assert.ok(maxAbsDiff(buf, ref) === 0);
});

test('refactoring after clear() and reassembly gives fresh results', () => {
  const rand = rng(4);
  const sys = new BlockTridiagonal(8, 4);
  randomSystem(sys, rand);
  sys.factor();
  sys.clear();
  assert.equal(sys.factored, false);
  const d = randomSystem(sys, rand);
  sys.factor();
  assert.ok(maxAbsDiff(sys.solve(d), denseSolve(sys, d)) < 1e-12);
});

test('reports an exactly singular diagonal block', () => {
  const sys = new BlockTridiagonal(3, 2);
  sys.B.set([1, 0, 0, 1], 0);
  sys.B.set([1, 2, 2, 4], 4); // singular
  sys.B.set([1, 0, 0, 1], 8);
  assert.throws(() => sys.factor(), /block 1 is singular/);
  assert.throws(() => sys.solve(new Float64Array(6)), /factor\(\) before solve/);
});

test('tracks the smallest pivot ratio as a conditioning hint', () => {
  const sys = new BlockTridiagonal(2, 2);
  sys.B.set([1, 0, 0, 1e-9], 0);
  sys.B.set([1, 0, 0, 1], 4);
  sys.factor();
  assert.ok(Math.abs(sys.minPivotRatio - 1e-9) < 1e-24);
});

test('rejects bad sizes', () => {
  assert.throws(() => new BlockTridiagonal(0, 2), RangeError);
  assert.throws(() => new BlockTridiagonal(3, 1.5), RangeError);
});
