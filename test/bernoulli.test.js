import { test } from 'node:test';
import assert from 'node:assert/strict';
import { bernoulli, bernoulliDerivative } from '../src/bernoulli.js';

const relErr = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

// Reference Taylor series with more terms than the implementation uses (valid |x| ≤ 0.1).
function seriesB(x) {
  const c = [1, -1 / 2, 1 / 12, 0, -1 / 720, 0, 1 / 30240, 0, -1 / 1209600, 0, 1 / 47900160, 0, -691 / 1307674368000];
  let s = 0;
  for (let n = c.length - 1; n >= 0; n--) s = s * x + c[n];
  return s;
}

const sweep = [];
for (let e = -300; e <= 3; e += 0.25) sweep.push(10 ** e, -(10 ** e));
sweep.push(0, 0.1, -0.1, 0.0999999999, 0.1000000001, 700, 710, 800, -700, -800, 1e4, -1e4);

test('B(0) = 1 and B′(0) = −1/2', () => {
  assert.equal(bernoulli(0), 1);
  assert.equal(bernoulliDerivative(0), -0.5);
});

test('B matches a longer Taylor series near 0', () => {
  for (const x of [1e-12, 1e-6, 1e-3, 0.01, 0.05, 0.0999, -0.0999, -0.05, -1e-6]) {
    assert.ok(relErr(bernoulli(x), seriesB(x)) < 4e-16, `x=${x}`);
  }
});

test('identity B(−x) = B(x) + x holds across the range', () => {
  for (const x of sweep) {
    if (Math.abs(x) > 700) continue; // B(x) underflows; the identity is then trivially B(−x) ≈ x
    // Absolute tolerance scales with the largest term, since B(x) + x can cancel.
    const lhs = bernoulli(-x), rhs = bernoulli(x) + x;
    assert.ok(Math.abs(lhs - rhs) <= 4e-16 * Math.max(1, Math.abs(x)), `x=${x}: ${lhs} vs ${rhs}`);
  }
});

test('B and B′ are continuous across the series/closed-form switch at |x| = 0.1', () => {
  // Adjacent doubles either side of the switch: series on one, closed form on the other.
  for (const s of [1, -1]) {
    const a = s * 0.09999999999999999, b = s * 0.1;
    assert.ok(relErr(bernoulli(a), bernoulli(b)) < 1e-15);
    // The closed-form B′ cancels (1 − B − x ≈ −x/2), costing a few ulps at the switch.
    assert.ok(relErr(bernoulliDerivative(a), bernoulliDerivative(b)) < 1e-14);
  }
});

test('B′ matches a central finite difference', () => {
  for (const x of [-30, -5, -1, -0.3, -0.1, -0.05, 1e-3, 0.05, 0.1, 0.3, 1, 5, 30]) {
    const h = 1e-5 * Math.max(1, Math.abs(x));
    const fd = (bernoulli(x + h) - bernoulli(x - h)) / (2 * h);
    assert.ok(Math.abs(bernoulliDerivative(x) - fd) < 1e-8 * Math.max(1, Math.abs(fd)), `x=${x}`);
  }
});

test('derivative identity B′(−x) = −B′(x) − 1', () => {
  for (const x of sweep) {
    if (Math.abs(x) > 700) continue;
    const lhs = bernoulliDerivative(-x), rhs = -bernoulliDerivative(x) - 1;
    assert.ok(Math.abs(lhs - rhs) < 1e-14, `x=${x}: ${lhs} vs ${rhs}`);
  }
});

test('no NaN or infinity for any argument; correct large-|x| limits', () => {
  for (const x of sweep) {
    assert.ok(Number.isFinite(bernoulli(x)), `B(${x})`);
    assert.ok(Number.isFinite(bernoulliDerivative(x)), `B′(${x})`);
    assert.ok(bernoulli(x) >= 0, `B(${x}) ≥ 0`);
  }
  assert.equal(bernoulli(-800), 800);
  assert.ok(bernoulli(800) < 1e-300);
  assert.ok(Math.abs(bernoulliDerivative(-800) + 1) < 1e-15);
  assert.ok(Math.abs(bernoulliDerivative(800)) < 1e-300);
});
