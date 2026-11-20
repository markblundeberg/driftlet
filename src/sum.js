// An affinity, a = fixedA + Σ ν η over a reaction's participants, summed without losing its
// digits. Each η carries its species' μ°/RT, which can be tens, while the affinity a fast reaction
// settles at is its rate over kf: summed term by term in doubles, the rounding of the largest
// partial sum (7e-15 at 40) is the floor under a, and kf times that floor is a rate. In a film of
// A + B + E = C + D with μ°/RT of 40, −40 and 0.3, kf 1e14, the flux out moved by 5e-4 with the
// order the equation was written in (E's digits dropped into A's partial sum before A cancelled
// against B; round 4, B6). So each term's product and each partial sum keep their rounding in a
// low word, and the sum is good to about an ulp of a itself, in any order.
//
// The accumulator is two doubles, acc[0] the high word and acc[1] the low, so that a + its
// rounding is acc[0] + acc[1]: set acc[0] = the start (fixedA), acc[1] = 0, add each term with
// addTerm, and read acc[0] + acc[1]. A non-finite term or sum keeps no rounding (a participant at
// c = 0 has η = −∞, and the sum is that ∞).

const SPLIT = 134217729; // 2^27 + 1, Dekker's split of a double into two 26-bit halves

/** acc[0] + acc[1] += nu·(x + xLo), keeping the roundings of nu·x and of the sum in acc[1]. */
export function addTerm(acc, nu, x, xLo) {
  const hi = acc[0], p = nu * x, s = hi + p;
  acc[0] = s;
  if (s - s !== 0) return; // ±∞ or NaN
  // TwoSum (Neumaier's ordering): hi + p = s + e exactly
  let e = Math.abs(hi) >= Math.abs(p) ? hi - s + p : p - s + hi;
  if (nu !== 1 && nu !== -1) {
    // TwoProduct (Dekker): nu·x = p + its rounding, exactly
    const cx = SPLIT * x, xh = cx - (cx - x), xl = x - xh;
    const cn = SPLIT * nu, nh = cn - (cn - nu), nl = nu - nh;
    e += nh * xh - p + nh * xl + nl * xh + nl * xl;
  }
  acc[1] += e + nu * xLo;
}
