// The Bernoulli function B(x) = x / (eˣ − 1) and its derivative, evaluated stably.
//
// These are the weights of Scharfetter–Gummel fluxes. Near x = 0 both are computed from
// their Taylor series; elsewhere B uses expm1 (accurate for small |x|, and it overflows
// gracefully to B → 0 for large positive x). Identity used: B(−x) = B(x) + x.

const SERIES_LIMIT = 0.1;

/** B(x) = x / (eˣ − 1), with B(0) = 1. */
export function bernoulli(x) {
  if (Math.abs(x) < SERIES_LIMIT) {
    const x2 = x * x;
    // 1 − x/2 + x²/12 − x⁴/720 + x⁶/30240 − x⁸/1209600
    return 1 - x / 2 + x2 * (1 / 12 + x2 * (-1 / 720 + x2 * (1 / 30240 - x2 / 1209600)));
  }
  return x / Math.expm1(x);
}

/** dB/dx, with B'(0) = −1/2. */
export function bernoulliDerivative(x) {
  if (Math.abs(x) < SERIES_LIMIT) {
    const x2 = x * x;
    // −1/2 + x/6 − x³/180 + x⁵/5040 − x⁷/151200
    return -0.5 + x * (1 / 6 + x2 * (-1 / 180 + x2 * (1 / 5040 - x2 / 151200)));
  }
  const b = x / Math.expm1(x);
  // B' = (B/x)(1 − B(−x)) = (B/x)(1 − B − x)
  return (b / x) * (1 - b - x);
}
