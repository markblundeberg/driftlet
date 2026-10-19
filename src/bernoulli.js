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

// Past a drift of SG_FAR thermal units across a segment, a Scharfetter–Gummel flux's B(d) heads
// for underflow and its expm1(Δη) for overflow (both by d ≈ 710, where their product is 0·∞):
// the product is then taken whole.
export const SG_FAR = 300;

/**
 * A Scharfetter–Gummel flux N = −g·B(d)·c_L·expm1(Δη) and its derivatives, into `out`: [N, ∂N/∂d,
 * g·B(d)·c_L (= ∂N/∂η_L), g·B(d)·c_L·e^{Δη} (= −∂N/∂η_R), −g·expm1(Δη)·(B(d) + B′(d)) (the
 * non-ideal kernel's ∂N/∂ζ_L per unit K_L)]. d is the drift, Δη the difference of levels, so
 * that c_R = c_L e^{Δη−d}. From d = SG_FAR on, the same flux from the downstream end, N =
 * g·B(−d)·c_R·expm1(−Δη), whose factors stay finite; ∂N/∂d = N·B′(d)/B(d).
 */
export function sgFlux(g, d, deta, cL, cR, out) {
  if (d < SG_FAR) {
    const E = Math.expm1(deta), B = bernoulli(d), Bp = bernoulliDerivative(d), gBc = g * B * cL;
    out[0] = -gBc * E;
    out[1] = -g * Bp * cL * E;
    out[2] = gBc;
    out[3] = gBc * (E + 1);
    out[4] = -g * E * (B + Bp);
    return out;
  }
  const B = bernoulli(d), up = g * bernoulli(-d) * cR, L = (1 - B - d) / d; // L = B′/B
  out[0] = up * Math.expm1(-deta);
  out[1] = out[0] * L;
  out[2] = g * B * cL;
  out[3] = up;
  // (B·expm1(Δη), taken whole: d e^{Δη−d}(1 − e^{−Δη})/(1 − e^{−d}))
  out[4] = -g * ((d * Math.exp(deta - d) * -Math.expm1(-deta)) / -Math.expm1(-d)) * (1 + L);
  return out;
}
