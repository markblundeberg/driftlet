// Powers that come out the same on every JavaScript engine. Math.pow (and **) doesn't: V8's
// changed between Node 22 and 24 by an ulp here and there, even for integer exponents, which was
// enough to flip a few accept/reject decisions in long adaptive transients. Integer exponents
// go by repeated squaring, half-integers add a correctly rounded square root, and anything else
// goes through exp and log (which agree across engine versions).

/** x^n for an integer n. */
export function powi(x, n) {
  if (n < 0) return 1 / powi(x, -n);
  let r = 1, b = x;
  while (n > 0) {
    if (n & 1) r *= b;
    b *= b;
    n >>= 1;
  }
  return r;
}

/** x^y for x ≥ 0 and any real y. */
export function powr(x, y) {
  if (Number.isInteger(y)) return powi(x, y);
  if (Number.isInteger(2 * y)) return powi(x, (2 * y - 1) / 2) * Math.sqrt(x);
  return x === 0 ? (y > 0 ? 0 : Infinity) : Math.exp(y * Math.log(x));
}
