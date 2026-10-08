// Statistics: how a material's concentrations follow from its species' potentials.
//
// Every species i in a material has a reduced chemical potential
//   ζ_i = (μ̄_i − z_i F φ − μ°_i) / RT,
// and a statistics model maps the ζ of the species it covers to their concentrations c(ζ) and
// the Jacobian K = ∂c/∂ζ (the chemical capacitance matrix, per RT). K is the Hessian of a
// convex potential, so it's symmetric and positive definite. Ideal (dilute) statistics are
// c = c_ref e^ζ, K = diag(c), and every model here reduces to them in the dilute limit, so μ°
// and c_ref keep one meaning (the dilute, Henry's-law reference) whatever the model.
//
// A model object covers `idx` (device species indices, k of them) and provides
//   evaluate(zeta, c, K, bg)       fill c (k) and K (k×k, row-major) at ζ (k)
//   invert(zeta, fixed, target, bg) for each a with fixed[a], set ζ_a so that c_a = target[a],
//                                   the other ζ held
// `bg` is the per-node background carrier density of an insertion host (0 elsewhere).

import { FARADAY, EPS0, AVOGADRO } from './constants.js';
import { fermiHalf, fermiMinusHalf } from './fermi.js';
import { powi } from './pow.js';
import { SolverError } from './solver.js';

const MODEL_TYPES = ['ideal', 'fermi-dirac', 'lattice', 'redlich-kister', 'debye-huckel', 'insertion', 'custom'];

const logistic = (y) => (y >= 0 ? 1 / (1 + Math.exp(-y)) : Math.exp(y) / (1 + Math.exp(y)));

// Safeguarded Newton for an increasing function: f(y) = target on [lo, hi]. A Newton step that
// leaves the bracket, or that wouldn't halve the last step (Newton zigzagging across an
// inflection, its bracket shrinking by a little each time: a steep tabulated OCV's ran out of
// iterations 0.2 thermal units from its root), bisects instead.
function solveIncreasing(f, target, lo, hi, guess) {
  let y = Math.min(hi, Math.max(lo, guess)), last = Infinity;
  for (let it = 0; it < 200; it++) {
    const { v, d } = f(y);
    const r = v - target;
    if (r > 0) hi = y;
    else lo = y;
    let next = y - r / d;
    if (!(next > lo && next < hi) || Math.abs(next - y) > 0.5 * last) next = 0.5 * (lo + hi);
    last = Math.abs(next - y);
    if (Math.abs(next - y) <= 1e-15 * Math.max(1, Math.abs(y)) || hi - lo <= 1e-15 * Math.max(1, Math.abs(y))) return next;
    y = next;
  }
  return y;
}

// --- occupancy of a single set of sites: x ∈ (0, 1) against a reduced potential ζ
//
// ζ(x) is strictly increasing. RK: ζ = ln(x/(1−x)) + ln(c_max/c_ref) + (g(x) − g(0))/RT, with
// g the Redlich–Kister excess chemical potential. Table: ζ(x) interpolated from data.

function redlichKisterOccupancy(A, cMax, cRef, RT, fail) {
  // Excess Gibbs energy per site x(1−x) Σ A_k u^k with u = 1 − 2x; g = its x-derivative.
  const P = (u, d) => {
    // d-th derivative of Σ A_k u^k
    let s = 0;
    for (let k = A.length - 1; k >= d; k--) {
      let coef = A[k];
      for (let q = 0; q < d; q++) coef *= k - q;
      s = s * u + coef;
    }
    // Horner over the shifted polynomial above builds Σ coef_k u^{k−d}
    return s;
  };
  const g = (x) => {
    const u = 1 - 2 * x;
    return u * P(u, 0) - ((1 - u * u) / 2) * P(u, 1);
  };
  const gp = (x) => {
    const u = 1 - 2 * x;
    return -2 * P(u, 0) - 4 * u * P(u, 1) + (1 - u * u) * P(u, 2);
  };
  const g0 = g(0), c0 = Math.log(cMax / cRef);
  // Convexity: dζ/dy = 1 + x(1−x) g′(x)/RT > 0 on (0, 1), where y = ln(x/(1−x)).
  let gmin = 0, gmax = 0;
  for (let k = 0; k <= 4000; k++) {
    const x = k / 4000;
    const slope = 1 + (x * (1 - x) * gp(x)) / RT;
    if (!(slope > 1e-9)) {
      fail(
        `the free energy is not convex near x = ${x.toFixed(3)} (a miscibility gap). Phase separation needs a ` +
          'Cahn–Hilliard model, which is on the roadmap',
      );
    }
    const e = (g(x) - g0) / RT;
    gmin = Math.min(gmin, e);
    gmax = Math.max(gmax, e);
  }
  return {
    zeta(x) {
      return Math.log(x / (1 - x)) + c0 + (g(x) - g0) / RT;
    },
    // x(ζ), 1 − x, and dx/dζ
    x(zeta) {
      const f = (y) => {
        const x = logistic(y);
        const omx = logistic(-y);
        return { v: y + (g(x) - g0) / RT, d: 1 + (x * omx * gp(x)) / RT };
      };
      const t = zeta - c0;
      const y = A.length === 0 ? t : solveIncreasing(f, t, t - gmax - 1, t - gmin + 1, t);
      const x = logistic(y), omx = logistic(-y);
      return { x, omx, dxdz: (x * omx) / f(y).d };
    },
  };
}

// Tabulated ζ(x), as ζ = ln(x/(1−x)) + r(x): the residual r (the non-ideality) is a
// shape-preserving cubic (PCHIP) through the data, continued linearly beyond the ends. Ideal
// and regular-solution tables are reproduced exactly, and x stays within (0, 1) for any ζ.
function tableOccupancy(xs, zs, fail) {
  const N = xs.length;
  const logit = (x) => Math.log(x / (1 - x));
  const rs = xs.map((x, k) => zs[k] - logit(x));
  const h = [], del = [];
  for (let k = 0; k < N - 1; k++) {
    h.push(xs[k + 1] - xs[k]);
    del.push((rs[k + 1] - rs[k]) / h[k]);
  }
  const d = new Float64Array(N);
  for (let k = 1; k < N - 1; k++) {
    if (del[k - 1] * del[k] <= 0) continue; // local extremum: flat
    const w1 = 2 * h[k] + h[k - 1], w2 = h[k] + 2 * h[k - 1];
    d[k] = (w1 + w2) / (w1 / del[k - 1] + w2 / del[k]);
  }
  const end = (dl, dn, hl, hn) => {
    if (N === 2) return dl;
    const s = ((2 * hl + hn) * dl - hl * dn) / (hl + hn);
    if (s * dl <= 0) return 0;
    if (dl * dn <= 0 && Math.abs(s) > 3 * Math.abs(dl)) return 3 * dl;
    return s;
  };
  d[0] = end(del[0], del[1], h[0], h[1]);
  d[N - 1] = end(del[N - 2], del[N - 3], h[N - 2], h[N - 3]);
  // r(x) and r′(x)
  const r = (x) => {
    if (x <= xs[0]) return { v: rs[0] + d[0] * (x - xs[0]), d: d[0] };
    if (x >= xs[N - 1]) return { v: rs[N - 1] + d[N - 1] * (x - xs[N - 1]), d: d[N - 1] };
    let lo = 0, hi = N - 2;
    while (lo < hi) {
      const m = (lo + hi + 1) >> 1;
      if (xs[m] <= x) lo = m;
      else hi = m - 1;
    }
    const k = lo, t = (x - xs[k]) / h[k], t2 = t * t, t3 = t2 * t;
    return {
      v: (2 * t3 - 3 * t2 + 1) * rs[k] + (t3 - 2 * t2 + t) * h[k] * d[k] + (-2 * t3 + 3 * t2) * rs[k + 1] + (t3 - t2) * h[k] * d[k + 1],
      d: ((6 * t2 - 6 * t) * rs[k] + (3 * t2 - 4 * t + 1) * h[k] * d[k] + (-6 * t2 + 6 * t) * rs[k + 1] + (3 * t2 - 2 * t) * h[k] * d[k + 1]) / h[k],
    };
  };
  // In y = ln(x/(1−x)): ζ = y + r(x), dζ/dy = 1 + x(1−x) r′(x), which must stay positive.
  const f = (y) => {
    const x = logistic(y), q = r(x);
    return { v: y + q.v, d: 1 + x * logistic(-y) * q.d };
  };
  let rmin = Infinity, rmax = -Infinity;
  for (let k = 0; k <= 4000; k++) {
    const x = k / 4000;
    const q = r(x);
    if (x > 0 && x < 1 && !(1 + x * (1 - x) * q.d > 1e-9)) {
      fail(`ocv: the interpolated curve isn't monotone near x = ${x.toFixed(3)} (E must fall steadily with x; add points or smooth the data)`);
    }
    rmin = Math.min(rmin, q.v);
    rmax = Math.max(rmax, q.v);
  }
  return {
    zeta: (x) => logit(x) + r(x).v,
    x(zeta) {
      const y = solveIncreasing(f, zeta, zeta - rmax - 1, zeta - rmin + 1, zeta - r(0.5).v);
      const x = logistic(y), omx = logistic(-y);
      return { x, omx, dxdz: (x * omx) / f(y).d };
    },
  };
}

// --- models

function fermiDiracModel(idx, cRef, order) {
  const F = order === 0.5 ? fermiHalf : (x) => (x > 0 ? x + Math.log1p(Math.exp(-x)) : Math.log1p(Math.exp(x)));
  const dF = order === 0.5 ? fermiMinusHalf : logistic;
  const k = idx.length;
  return {
    type: 'fermi-dirac',
    idx,
    evaluate(zeta, c, K) {
      K.fill(0);
      for (let a = 0; a < k; a++) {
        c[a] = cRef[a] * F(zeta[a]);
        K[a * k + a] = cRef[a] * dF(zeta[a]);
      }
    },
    invert(zeta, fixed, target) {
      for (let a = 0; a < k; a++) {
        if (!fixed[a]) continue;
        const r = target[a] / cRef[a];
        // ln 𝓕(ζ) is increasing with slope 𝓕′/𝓕 ∈ (0, 1]
        const f = (y) => ({ v: Math.log(F(y)), d: dF(y) / F(y) });
        const guess = r < 1 ? Math.log(r) : order === 0.5 ? Math.cbrt(powi(0.75 * Math.sqrt(Math.PI) * r, 2)) : r;
        zeta[a] = solveIncreasing(f, Math.log(r), Math.log(r) - 1, Math.max(guess, Math.log(r)) + 10, guess);
      }
    },
  };
}

function latticeModel(idx, cRef, cMax) {
  const k = idx.length;
  const lref = cRef.map((c) => Math.log(c / cMax));
  return {
    type: 'lattice',
    idx,
    // c_a = c_max w_a / (1 + Σ w), w_a = (c_ref,a / c_max) e^{ζ_a}; K = diag(c) − c cᵀ/c_max.
    evaluate(zeta, c, K) {
      let m = 0;
      for (let a = 0; a < k; a++) m = Math.max(m, zeta[a] + lref[a]);
      let S = Math.exp(-m);
      for (let a = 0; a < k; a++) {
        c[a] = Math.exp(zeta[a] + lref[a] - m);
        S += c[a];
      }
      for (let a = 0; a < k; a++) c[a] *= cMax / S;
      for (let a = 0; a < k; a++) {
        for (let b = 0; b < k; b++) K[a * k + b] = (a === b ? c[a] : 0) - (c[a] * c[b]) / cMax;
      }
    },
    invert(zeta, fixed, target) {
      let T = 0, m = 0;
      for (let a = 0; a < k; a++) {
        if (fixed[a]) T += target[a] / cMax;
        else m = Math.max(m, zeta[a] + lref[a]);
      }
      if (!(T < 1)) throw new SolverError(`lattice: concentrations exceed the site density ${cMax}`);
      let s = Math.exp(-m);
      for (let a = 0; a < k; a++) if (!fixed[a]) s += Math.exp(zeta[a] + lref[a] - m);
      const ln1W = m + Math.log(s); // ln(1 + Σ_unfixed w)
      for (let a = 0; a < k; a++) if (fixed[a]) zeta[a] = Math.log(target[a] / cRef[a]) + ln1W - Math.log1p(-T);
    },
  };
}

// Single species on its own sites, with a Redlich–Kister excess (Langmuir when A is empty).
function redlichKisterModel(idx, cRef, cMax, occ) {
  return {
    type: 'redlich-kister',
    idx,
    evaluate(zeta, c, K) {
      const o = occ.x(zeta[0]);
      c[0] = cMax * o.x;
      K[0] = cMax * o.dxdz;
    },
    invert(zeta, fixed, target) {
      if (!fixed[0]) return;
      const x = target[0] / cMax;
      if (!(x > 0 && x < 1)) throw new SolverError(`redlich-kister: concentration ${target[0]} is outside (0, c_max)`);
      zeta[0] = occ.zeta(x);
    },
  };
}

// Extended Debye–Hückel: ln γ_i = −z_i² ℓ_B κ / (2(1 + κa)), with κ² = F² Σ z² c / (ε RT) and
// c_i = c_ref,i e^{ζ_i} / γ_i solved by Newton in ln c. The excess Hessian is rank one, so
// K = diag(c) + A (s∘c)(s∘c)ᵀ / (1 − A Σ s² c), with s = z² and
// A = ℓ_B F² / (4 ε RT κ (1 + κa)²). Convex while 1 − A Σ s² c > 0.
function debyeHuckelModel(idx, cRef, z, epsr, a, RT) {
  const k = idx.length;
  const eps = epsr * EPS0;
  const lB = (FARADAY * FARADAY) / (4 * Math.PI * eps * AVOGADRO * RT);
  const s = z.map((zi) => zi * zi);
  const y = new Float64Array(k), cc = new Float64Array(k), lng = new Float64Array(k);
  const state = (c) => {
    let q = 0;
    for (let b = 0; b < k; b++) q += s[b] * c[b];
    const kappa = Math.sqrt((FARADAY * FARADAY * q) / (eps * RT));
    const A = kappa > 0 ? (lB * FARADAY * FARADAY) / (4 * eps * RT * kappa * (1 + kappa * a) * (1 + kappa * a)) : 0;
    let s2c = 0;
    for (let b = 0; b < k; b++) s2c += s[b] * s[b] * c[b];
    for (let b = 0; b < k; b++) lng[b] = (-s[b] * lB * kappa) / (2 * (1 + kappa * a));
    const denom = 1 - A * s2c;
    if (!(denom > 0)) {
      throw new SolverError(
        `debye-huckel: the activity model is not convex at κ = ${kappa.toExponential(3)} /m (ionic strength too high for Debye–Hückel)`,
      );
    }
    return { A, denom };
  };
  return {
    type: 'debye-huckel',
    idx,
    evaluate(zeta, c, K) {
      for (let b = 0; b < k; b++) y[b] = zeta[b] + Math.log(cRef[b]);
      let st;
      for (let it = 0; it < 100; it++) {
        for (let b = 0; b < k; b++) cc[b] = Math.exp(y[b]);
        st = state(cc);
        // R = y + ln γ − ζ − ln c_ref; J = I − A s (s∘c)ᵀ; J⁻¹R = R + A s ((s∘c)·R) / denom
        let dot = 0;
        const R = c; // scratch
        for (let b = 0; b < k; b++) {
          R[b] = y[b] + lng[b] - zeta[b] - Math.log(cRef[b]);
          dot += s[b] * cc[b] * R[b];
        }
        let mx = 0;
        for (let b = 0; b < k; b++) {
          let dy = R[b] + (st.A * s[b] * dot) / st.denom;
          if (dy > 1) dy = 1;
          else if (dy < -1) dy = -1;
          y[b] -= dy;
          mx = Math.max(mx, Math.abs(dy));
        }
        if (mx < 1e-14) break;
      }
      for (let b = 0; b < k; b++) c[b] = Math.exp(y[b]);
      st = state(c);
      const f = st.A / st.denom;
      for (let p = 0; p < k; p++) {
        for (let q = 0; q < k; q++) K[p * k + q] = (p === q ? c[p] : 0) + f * s[p] * c[p] * s[q] * c[q];
      }
    },
    invert(zeta, fixed, target) {
      for (let b = 0; b < k; b++) cc[b] = fixed[b] ? target[b] : cRef[b] * Math.exp(zeta[b]);
      for (let it = 0; it < 500; it++) {
        state(cc);
        let change = 0;
        for (let b = 0; b < k; b++) {
          if (fixed[b]) continue;
          const next = cRef[b] * Math.exp(zeta[b] - lng[b]);
          if (cc[b] > 0) change = Math.max(change, Math.abs(next / cc[b] - 1));
          cc[b] = next;
        }
        if (!(change > 1e-15)) break;
      }
      state(cc);
      for (let b = 0; b < k; b++) if (fixed[b]) zeta[b] = Math.log(target[b] / cRef[b]) + lng[b];
    },
  };
}

// A neutral combination ion + ν carrier in a strictly neutral host (e.g. Li⁺ + e⁻). The host's
// content s depends only on ζ_ion + ν ζ_carrier (the combination's chemical potential, i.e. the
// OCV curve); the carrier is s ν plus the host's own background carriers.
//   c_ion = s,  c_carrier = ν s + bg,  K = s′ [1 ν; ν ν²]
function insertionModel(idx, nu, cMax, occ) {
  return {
    type: 'insertion',
    idx,
    evaluate(zeta, c, K, bg) {
      const o = occ.x(zeta[0] + nu * zeta[1]);
      const s = cMax * o.x, ds = cMax * o.dxdz;
      c[0] = s;
      c[1] = nu * s + bg;
      K[0] = ds;
      K[1] = K[2] = nu * ds;
      K[3] = nu * nu * ds;
    },
    invert(zeta, fixed, target, bg) {
      if (fixed[0]) {
        const x = target[0] / cMax;
        if (!(x > 0 && x < 1)) throw new SolverError(`insertion: concentration ${target[0]} is outside (0, c_max)`);
        const zc = occ.zeta(x);
        if (fixed[1] && Math.abs(target[1] - (nu * target[0] + bg)) > 1e-9 * target[1]) {
          throw new SolverError(`insertion: carrier concentration ${target[1]} is inconsistent with ${nu}·${target[0]} + background ${bg}`);
        }
        zeta[0] = zc - nu * zeta[1];
      } else if (fixed[1]) {
        const x = (target[1] - bg) / (nu * cMax);
        if (!(x > 0 && x < 1)) throw new SolverError(`insertion: carrier concentration ${target[1]} is outside the host's range`);
        zeta[1] = (occ.zeta(x) - zeta[0]) / nu;
      }
    },
  };
}

// User function: evaluate(zeta) → { c: number[], dcdzeta: number[][] }.
function customModel(idx, fn, fail) {
  const k = idx.length;
  const call = (zeta, c, K) => {
    const out = fn(Array.from(zeta));
    for (let a = 0; a < k; a++) {
      c[a] = out.c[a];
      for (let b = 0; b < k; b++) K[a * k + b] = out.dcdzeta[a][b];
    }
  };
  // Construction check at ζ = 0: finite, positive, symmetric, positive definite.
  const c = new Float64Array(k), K = new Float64Array(k * k);
  let ok = true;
  try {
    call(new Float64Array(k), c, K);
  } catch (err) {
    fail(`evaluate failed at ζ = 0: ${err.message}`);
  }
  for (let a = 0; a < k; a++) {
    if (!(c[a] > 0 && Number.isFinite(c[a]))) ok = false;
    for (let b = 0; b < k; b++) {
      if (!Number.isFinite(K[a * k + b])) ok = false;
      if (Math.abs(K[a * k + b] - K[b * k + a]) > 1e-8 * (Math.abs(K[a * k + a]) + Math.abs(K[b * k + b]))) {
        fail("evaluate: ∂c/∂ζ must be symmetric (it's the Hessian of a potential)");
      }
    }
  }
  if (!ok) fail('evaluate must return positive, finite c and a finite dcdzeta');
  if (!choleskyOk(K, k)) fail('evaluate: ∂c/∂ζ must be positive definite (the potential must be convex)');
  const r = new Float64Array(k), J = new Float64Array(k * k), dz = new Float64Array(k);
  return {
    type: 'custom',
    idx,
    evaluate: call,
    invert(zeta, fixed, target) {
      const cc = new Float64Array(k), KK = new Float64Array(k * k);
      const f = [];
      for (let a = 0; a < k; a++) if (fixed[a]) f.push(a);
      if (f.length === 0) return;
      for (let it = 0; it < 200; it++) {
        call(zeta, cc, KK);
        const m = f.length;
        for (let p = 0; p < m; p++) {
          r[p] = cc[f[p]] - target[f[p]];
          for (let q = 0; q < m; q++) J[p * m + q] = KK[f[p] * k + f[q]];
        }
        denseSolve(J, r, dz, m);
        let mx = 0;
        for (let p = 0; p < m; p++) mx = Math.max(mx, Math.abs(dz[p]));
        const scale = mx > 2 ? 2 / mx : 1;
        for (let p = 0; p < m; p++) zeta[f[p]] -= scale * dz[p];
        if (mx < 1e-14) return;
      }
      throw new SolverError('custom statistics: could not find ζ for the given concentrations');
    },
  };
}

function choleskyOk(K, k) {
  const L = Float64Array.from(K);
  for (let j = 0; j < k; j++) {
    let d = L[j * k + j];
    for (let p = 0; p < j; p++) d -= L[j * k + p] * L[j * k + p];
    if (!(d > 0)) return false;
    d = Math.sqrt(d);
    L[j * k + j] = d;
    for (let i = j + 1; i < k; i++) {
      let v = L[i * k + j];
      for (let p = 0; p < j; p++) v -= L[i * k + p] * L[j * k + p];
      L[i * k + j] = v / d;
    }
  }
  return true;
}

// Gaussian elimination with partial pivoting (small dense systems).
function denseSolve(A, b, x, m) {
  const M = Float64Array.from(A.subarray(0, m * m)), v = Float64Array.from(b.subarray(0, m));
  for (let c = 0; c < m; c++) {
    let p = c;
    for (let r = c + 1; r < m; r++) if (Math.abs(M[r * m + c]) > Math.abs(M[p * m + c])) p = r;
    if (p !== c) {
      for (let q = 0; q < m; q++) [M[c * m + q], M[p * m + q]] = [M[p * m + q], M[c * m + q]];
      [v[c], v[p]] = [v[p], v[c]];
    }
    for (let r = c + 1; r < m; r++) {
      const f = M[r * m + c] / M[c * m + c];
      for (let q = c; q < m; q++) M[r * m + q] -= f * M[c * m + q];
      v[r] -= f * v[c];
    }
  }
  for (let r = m - 1; r >= 0; r--) {
    let s = v[r];
    for (let q = r + 1; q < m; q++) s -= M[r * m + q] * x[q];
    x[r] = s / M[r * m + r];
  }
}

/**
 * Validate a material's `statistics` list and build its models.
 * @returns {{ models: object[], modelOf: Int32Array, warnings: string[], phiFree: boolean }}
 */
export function normalizeStatistics(mat, path, species, speciesIndex, RT, h) {
  const { need, finite, positive } = h;
  const n = species.length;
  const modelOf = new Int32Array(n).fill(-1);
  const models = [];
  const warnings = [];
  const list = mat.statistics ?? [];
  need(Array.isArray(list), `${path}.statistics must be an array of statistics models`);
  list.forEach((sdef, m) => {
    const sp = `${path}.statistics[${m}]`;
    const fail = (msg) => need(false, `${sp}: ${msg}`);
    need(sdef !== null && typeof sdef === 'object' && !Array.isArray(sdef), `${sp} must be an object`);
    need(MODEL_TYPES.includes(sdef.type), `${sp}.type must be one of ${MODEL_TYPES.join(', ')}`);
    const allowed = {
      'fermi-dirac': ['order'],
      lattice: ['cMax'],
      'redlich-kister': ['cMax', 'A'],
      'debye-huckel': ['epsr', 'a'],
      insertion: ['cMax', 'A', 'ocv'],
      custom: ['evaluate'],
    }[sdef.type];
    for (const k of Object.keys(sdef)) need(k === 'type' || k === 'species' || allowed.includes(k), `${sp}.${k}: not a field of a ${sdef.type} model (type, species, ${allowed.join(', ')})`);
    need(Array.isArray(sdef.species) && sdef.species.length > 0, `${sp}.species must be a non-empty array of species names`);
    const idx = sdef.species.map((name, a) => {
      need(speciesIndex.has(name), `${sp}.species[${a}]: unknown species ${JSON.stringify(name)}`);
      const i = speciesIndex.get(name);
      need(mat.present[i], `${sp}.species[${a}]: '${name}' is absent from this material`);
      need(modelOf[i] === -1, `${sp}.species[${a}]: '${name}' already belongs to another statistics model`);
      modelOf[i] = models.length;
      return i;
    });
    const cRef = idx.map((i) => mat.cRef[i]);
    const one = () => need(idx.length === 1, `${sp}.species: a ${sdef.type} model covers exactly one species`);
    let model;
    switch (sdef.type) {
      case 'ideal':
        idx.forEach((i) => (modelOf[i] = -1));
        return;
      case 'fermi-dirac': {
        const order = sdef.order ?? 0.5;
        need(order === 0.5 || order === 0, `${sp}.order must be 1/2 (3D parabolic band) or 0 (2D)`);
        model = fermiDiracModel(idx, cRef, order);
        break;
      }
      case 'lattice':
        model = latticeModel(idx, cRef, positive(sdef.cMax, `${sp}.cMax`));
        break;
      case 'redlich-kister': {
        one();
        const cMax = positive(sdef.cMax, `${sp}.cMax`);
        const A = sdef.A ?? [];
        need(Array.isArray(A), `${sp}.A must be an array of Redlich–Kister coefficients (J/mol)`);
        A.forEach((v, q) => finite(v, `${sp}.A[${q}]`));
        model = redlichKisterModel(idx, cRef, cMax, redlichKisterOccupancy(A, cMax, cRef[0], RT, fail));
        break;
      }
      case 'debye-huckel': {
        idx.forEach((i) => need(species[i].z !== 0, `${sp}.species: '${species[i].name}' is neutral, so Debye–Hückel doesn't apply`));
        const epsr = positive(sdef.epsr ?? (mat.epsr > 0 ? mat.epsr : undefined), `${sp}.epsr (the solvent permittivity for the theory)`);
        const a = sdef.a === undefined ? 0 : finite(sdef.a, `${sp}.a`);
        need(a >= 0, `${sp}.a (distance of closest approach, m) must be non-negative`);
        if (mat.epsr > 0) {
          warnings.push(
            `${sp}: Debye–Hückel activity coefficients split single-ion activities by convention, but a material ` +
              'with ε > 0 resolves single-ion potentials; this mixes the two. Prefer it in ε = 0 materials.',
          );
        }
        model = debyeHuckelModel(idx, cRef, idx.map((i) => species[i].z), epsr, a, RT);
        break;
      }
      case 'insertion': {
        need(idx.length === 2, `${sp}.species must be [ion, carrier]: the inserted ion and its compensating electronic carrier`);
        const [ion, car] = idx;
        const zi = species[ion].z, zc = species[car].z;
        const nu = -zi / zc;
        need(zi !== 0 && zc !== 0 && Number.isInteger(nu) && nu > 0, `${sp}.species: the ion's charge must be a positive multiple of the carrier's, with opposite sign`);
        need(mat.epsr === 0, `${sp}: insertion statistics describe a neutral combination, so the material needs epsr: 0`);
        need(
          species.every((s, i) => !mat.present[i] || s.z === 0 || i === ion || i === car),
          `${sp}: the host's only charged species must be the inserted ion and its carrier`,
        );
        const cMax = positive(sdef.cMax, `${sp}.cMax`);
        let occ;
        if (sdef.ocv !== undefined) {
          need(sdef.A === undefined, `${sp}: give either A (Redlich–Kister) or ocv (a table), not both`);
          const t = sdef.ocv;
          need(t !== null && typeof t === 'object' && Array.isArray(t.x) && Array.isArray(t.E), `${sp}.ocv must be { x: [...], E: [...], muRef }`);
          for (const k of Object.keys(t)) need(['x', 'E', 'muRef'].includes(k), `${sp}.ocv.${k}: not a field (x, E, muRef)`);
          need(t.x.length >= 2 && t.x.length === t.E.length, `${sp}.ocv: x and E need the same length, at least 2`);
          t.x.forEach((v, q) => {
            finite(v, `${sp}.ocv.x[${q}]`);
            need(v > 0 && v < 1, `${sp}.ocv.x[${q}] must be a site fraction in (0, 1)`);
            need(q === 0 || v > t.x[q - 1], `${sp}.ocv.x must be strictly increasing`);
          });
          t.E.forEach((v, q) => {
            finite(v, `${sp}.ocv.E[${q}]`);
            need(q === 0 || v < t.E[q - 1], `${sp}.ocv.E must be strictly decreasing with x (plateaus are phase separation, not yet supported)`);
          });
          const muRef = finite(t.muRef, `${sp}.ocv.muRef (the combination's chemical potential in the reference electrode, J/mol)`);
          // μ_ion + ν μ_carrier = muRef − z_ion F E, and ζ_comb = that minus the standard potentials.
          const mu0 = mat.mu0[ion] + nu * mat.mu0[car];
          occ = tableOccupancy(t.x, t.E.map((E) => (muRef - zi * FARADAY * E - mu0) / RT), fail);
        } else {
          const A = sdef.A ?? [];
          need(Array.isArray(A), `${sp}.A must be an array of Redlich–Kister coefficients (J/mol)`);
          A.forEach((v, q) => finite(v, `${sp}.A[${q}]`));
          occ = redlichKisterOccupancy(A, cMax, cRef[0], RT, fail);
        }
        model = insertionModel(idx, nu, cMax, occ);
        model.carrier = car;
        break;
      }
      case 'custom':
        need(typeof sdef.evaluate === 'function', `${sp}.evaluate must be a function (zeta) => ({ c, dcdzeta })`);
        model = customModel(idx, sdef.evaluate, fail);
        break;
    }
    models.push(model);
  });
  // Which ζ positions the models' species occupy, and remap modelOf to the kept models.
  const charged = species.map((s, i) => mat.present[i] && s.z !== 0);
  const inCombination = (i) => modelOf[i] >= 0 && models[modelOf[i]].type === 'insertion';
  // Strictly neutral, and every charged species in a neutral combination: nothing depends on φ.
  const phiFree = mat.epsr === 0 && charged.some(Boolean) && charged.every((ch, i) => !ch || inCombination(i));
  return { models, modelOf, warnings, phiFree };
}
