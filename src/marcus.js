// Marcus–Hush–Chidsey electron transfer at a metal: the Marcus rate for each electron energy in
// the metal, weighted by the Fermi function and summed over the band (Chidsey, Science 251, 919
// (1991)), which bends the Tafel line over and saturates it instead of running on exponentially.
//
// CONTRACT
//
// mhcIntegral(lam, a) is a pure function of its arguments: lam = λ/RT > 0, the reorganization
// energy in thermal units, and a the reaction's affinity in thermal units with every reactant and
// product at its standard state (a = 0 with the metal's Fermi level at the couple's standard level;
// a > 0 drives the forward reaction). It returns { I, dI }:
//   I(λ, a) = ∫ exp(−(x − λ + a)²/(4λ)) / (1 + eˣ) dx over all x (x the electron's energy from the
//   Fermi level, in RT), and dI = ∂I/∂a. The forward rate constant is proportional to I(λ, a), the
//   backward to I(λ, −a), and I(λ, −a) = e^{−a} I(λ, a) exactly (detailed balance; it's how a < 0
//   is computed, so the identity holds to round-off). I → √(4πλ) as a → ∞ (every occupied level
//   below the Gaussian: the rate saturates); for λ ≫ |a|, I(λ, a)/I(λ, 0) → exp(a/2 − a²/(4λ)),
//   Marcus's parabola, whose slope at a = 0 is Butler–Volmer's α = ½.
// The integrand is analytic in the strip |Im x| < π (the Fermi function's poles), so the trapezoid
// rule on a uniform grid converges geometrically: step h ≤ ½ (error ~e^{−2π²/h} < 1e−17) and h ≤
// √λ/2 for the Gaussian, over the range where the integrand is within e^{−46} of its peak.
// Relative accuracy ~1e−12 or better (test/marcus.test.js, against mpmath at 30 digits).

const K = 46; // the integrand's range: down to e^{−K} of its peak

// ln(1 + eˣ), without overflow.
const softplus = (x) => (x > 0 ? x + Math.log1p(Math.exp(-x)) : Math.log1p(Math.exp(x)));

// I(λ, b) and ∂I/∂b for b ≥ 0, where the Gaussian sits at or below λ.
function forward(lam, b) {
  const x0 = lam - b, w = 2 * Math.sqrt(lam * K);
  // The peak is at x0 when that's below the Fermi level, else between 0 and x0 (x0 ≤ λ keeps it
  // there): the Gaussian's tail to the left, and above 0 a decay at least as fast as e^{−x/2}, or
  // the Gaussian's own where that's narrower.
  const lo = x0 > 0 ? x0 - Math.sqrt(x0 * x0 + 4 * lam * K) : x0 - w;
  const hi = x0 > 0 ? Math.min(2 * K, x0 + w) : x0 + w;
  const h0 = Math.min(0.5, Math.sqrt(lam) / 2), m = Math.ceil((hi - lo) / h0), h = (hi - lo) / m;
  let I = 0, J = 0;
  for (let k = 0; k <= m; k++) {
    const x = lo + k * h, d = x - x0, g = Math.exp(-(d * d) / (4 * lam) - softplus(x));
    I += g;
    J += d * g;
  }
  return { I: I * h, dI: (-J * h) / (2 * lam) };
}

/** I(λ, a) and ∂I/∂a (see the contract above). */
export function mhcIntegral(lam, a) {
  a = Math.max(-1e6, Math.min(1e6, a)); // (an infinite affinity, from a concentration underflowed to 0, as a huge one)
  if (a >= 0) return forward(lam, a);
  // I(λ, a) = e^{a} I(λ, −a), so ∂I/∂a = e^{a} (I(λ, −a) − I′(λ, −a))
  const { I, dI } = forward(lam, -a), e = Math.exp(a);
  return { I: e * I, dI: e * (I - dI) };
}
