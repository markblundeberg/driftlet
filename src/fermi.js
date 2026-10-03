// Complete Fermi–Dirac integrals of order 1/2 and −1/2, normalised so that 𝓕_j(x) → e^x as
// x → −∞:
//   𝓕_j(x) = (1/Γ(j+1)) ∫₀^∞ t^j / (1 + e^{t−x}) dt,   d𝓕_j/dx = 𝓕_{j−1}.
// Relative accuracy ~1e-14 everywhere:
//   x < −2       alternating series Σ (−1)^{k+1} e^{kx} / k^{j+1}
//   −2 ≤ x ≤ 50  piecewise Chebyshev fits, built on first use from an exact quadrature
//   x > 50       Sommerfeld's asymptotic expansion (its error is ~e^{−x})

import { powi, powr } from './pow.js';

const INTERVALS = [-2, 1, 5, 12, 25, 50];
const NCHEB = 24;
const SQRT_PI = Math.sqrt(Math.PI);

// Reference values by quadrature. With t = u² the integrand is even and smooth in u, so the
// trapezoid rule converges geometrically; h = 0.02 already reaches round-off.
function quadrature(j, x) {
  const h = 0.02;
  const umax = Math.sqrt(Math.max(x, 0) + 45);
  const f = (u) => {
    const p = j === 0.5 ? 2 * u * u : 2;
    const e = u * u - x;
    return e > 0 ? (p * Math.exp(-e)) / (1 + Math.exp(-e)) : p / (1 + Math.exp(e));
  };
  let s = f(0) / 2;
  for (let k = 1; k * h <= umax; k++) s += f(k * h);
  return (s * h) / (j === 0.5 ? SQRT_PI / 2 : SQRT_PI);
}

function chebFit(j, a, b) {
  const v = new Float64Array(NCHEB);
  for (let k = 0; k < NCHEB; k++) {
    const t = Math.cos((Math.PI * (k + 0.5)) / NCHEB);
    v[k] = quadrature(j, 0.5 * (a + b) + 0.5 * (b - a) * t);
  }
  const c = new Float64Array(NCHEB);
  for (let m = 0; m < NCHEB; m++) {
    let s = 0;
    for (let k = 0; k < NCHEB; k++) s += v[k] * Math.cos((Math.PI * m * (k + 0.5)) / NCHEB);
    c[m] = (2 / NCHEB) * s;
  }
  c[0] /= 2;
  return c;
}

let fits = null; // [j = 1/2, j = −1/2] × intervals
function tables() {
  if (!fits) {
    fits = [0.5, -0.5].map((j) => INTERVALS.slice(0, -1).map((a, k) => chebFit(j, a, INTERVALS[k + 1])));
  }
  return fits;
}

function clenshaw(c, a, b, x) {
  const t = (2 * x - a - b) / (b - a);
  let b1 = 0, b2 = 0;
  for (let m = c.length - 1; m >= 1; m--) {
    const tmp = 2 * t * b1 - b2 + c[m];
    b2 = b1;
    b1 = tmp;
  }
  return t * b1 - b2 + c[0];
}

function series(j, x) {
  const e = Math.exp(x);
  let term = e, s = 0;
  for (let k = 1; k <= 40; k++) {
    const t = term / powr(k, j + 1);
    s += k % 2 ? t : -t;
    if (t < 1e-17 * s) break;
    term *= e;
  }
  return s;
}

// η(2k) = (1 − 2^{1−2k}) ζ(2k), the Dirichlet eta function at even arguments.
const ETA = [powi(Math.PI, 2) / 12, (7 * powi(Math.PI, 4)) / 720, (31 * powi(Math.PI, 6)) / 30240, 0.99623300185264789922, 0.99903950759827156564, 0.99975768514385819085];

function asymptotic(j, x) {
  // 𝓕_j(x) ~ x^{j+1}/Γ(j+2) · [1 + Σ_k 2η(2k) (j+1)(j)…(j+2−2k) x^{−2k}]
  let s = 1, fall = 1;
  const x2 = 1 / (x * x);
  let p = 1;
  for (let k = 1; k <= ETA.length; k++) {
    fall *= (j + 1 - (2 * k - 2)) * (j + 1 - (2 * k - 1));
    p *= x2;
    s += 2 * ETA[k - 1] * fall * p;
  }
  const gamma = j === 0.5 ? (3 * SQRT_PI) / 4 : SQRT_PI / 2; // Γ(j + 2)
  return (powr(x, j + 1) / gamma) * s;
}

function fermi(j, x) {
  if (x < INTERVALS[0]) return series(j, x);
  if (x > INTERVALS[INTERVALS.length - 1]) return asymptotic(j, x);
  const t = tables()[j === 0.5 ? 0 : 1];
  let k = 0;
  while (x > INTERVALS[k + 1]) k++;
  return clenshaw(t[k], INTERVALS[k], INTERVALS[k + 1], x);
}

/** 𝓕_{1/2}(x), normalised so that it tends to e^x for x → −∞. */
export const fermiHalf = (x) => fermi(0.5, x);
/** 𝓕_{−1/2}(x) = d𝓕_{1/2}/dx. */
export const fermiMinusHalf = (x) => fermi(-0.5, x);
/** Reference quadrature, for tests. */
export const fermiQuadrature = quadrature;
