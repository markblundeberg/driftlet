// Block-tridiagonal linear solver (block Thomas algorithm).
//
// Solves, for node i = 0 … n−1 with blocks of size m,
//     A_i x_{i−1} + B_i x_i + C_i x_{i+1} = d_i
// (A_0 and C_{n−1} are ignored). Blocks are stored row-major in flat Float64Arrays: block i
// occupies [i·m², (i+1)·m²). All storage is allocated once; factor() and solve() allocate
// nothing.
//
// Within each diagonal block, LU uses partial pivoting, so zero diagonal entries (such as a
// slot that carries an interface flux instead of a value) are fine. There is no pivoting
// *between* blocks: that is the usual block Thomas trade-off, safe for the diagonally
// dominant-ish systems that finite-volume assembly produces.

export class BlockTridiagonal {
  /**
   * @param {number} n number of block rows (nodes)
   * @param {number} m block size (unknowns per node)
   */
  constructor(n, m) {
    if (!(Number.isInteger(n) && n >= 1)) throw new RangeError(`n must be a positive integer, got ${n}`);
    if (!(Number.isInteger(m) && m >= 1)) throw new RangeError(`m must be a positive integer, got ${m}`);
    this.n = n;
    this.m = m;
    const mm = m * m;
    /** Sub-diagonal blocks A_i (coupling to node i−1). */
    this.A = new Float64Array(n * mm);
    /** Diagonal blocks B_i. */
    this.B = new Float64Array(n * mm);
    /** Super-diagonal blocks C_i (coupling to node i+1). */
    this.C = new Float64Array(n * mm);

    this._lu = new Float64Array(n * mm); // LU factors of W_i = B_i − A_i C'_{i−1}
    this._piv = new Int32Array(n * m);
    this._cp = new Float64Array(n * mm); // C'_i = W_i⁻¹ C_i
    this._y = new Float64Array(n * m);
    this._tmp = new Float64Array(m);
    this.factored = false;
    /** Smallest |pivot| / (largest |entry| of its block) seen in the last factor(). */
    this.minPivotRatio = Infinity;
  }

  /** Zero all of A, B, C (ready for fresh assembly). */
  clear() {
    this.A.fill(0);
    this.B.fill(0);
    this.C.fill(0);
    this.factored = false;
  }

  /**
   * Factorise the current A, B, C. A, B, C are left untouched.
   * Throws if a diagonal block is exactly singular.
   */
  factor() {
    const { n, m, A, B, C } = this;
    const mm = m * m;
    const lu = this._lu, piv = this._piv, cp = this._cp;
    let minRatio = Infinity;

    for (let i = 0; i < n; i++) {
      const o = i * mm;
      // W_i = B_i − A_i C'_{i−1}
      if (i === 0) {
        for (let k = 0; k < mm; k++) lu[o + k] = B[o + k];
      } else {
        const op = o - mm;
        for (let r = 0; r < m; r++) {
          for (let c = 0; c < m; c++) {
            let s = B[o + r * m + c];
            for (let k = 0; k < m; k++) s -= A[o + r * m + k] * cp[op + k * m + c];
            lu[o + r * m + c] = s;
          }
        }
      }
      const ratio = luFactor(lu, o, piv, i * m, m);
      if (!(ratio > 0)) {
        this.factored = false;
        throw new Error(`BlockTridiagonal: diagonal block ${i} is singular`);
      }
      if (ratio < minRatio) minRatio = ratio;
      // C'_i = W_i⁻¹ C_i, column by column
      if (i < n - 1) {
        const tmp = this._tmp;
        for (let c = 0; c < m; c++) {
          for (let r = 0; r < m; r++) tmp[r] = C[o + r * m + c];
          luSolve(lu, o, piv, i * m, m, tmp);
          for (let r = 0; r < m; r++) cp[o + r * m + c] = tmp[r];
        }
      }
    }
    this.minPivotRatio = minRatio;
    this.factored = true;
  }

  /**
   * Solve with the current factorisation. `rhs` and `out` are Float64Arrays of length n·m
   * (node-major); they may be the same array.
   */
  solve(rhs, out = new Float64Array(this.n * this.m)) {
    if (!this.factored) throw new Error('BlockTridiagonal: call factor() before solve()');
    const { n, m, A } = this;
    const mm = m * m;
    const lu = this._lu, piv = this._piv, cp = this._cp, y = this._y;

    // Forward: y_i = W_i⁻¹ (d_i − A_i y_{i−1})
    for (let i = 0; i < n; i++) {
      const o = i * mm, v = i * m;
      if (i === 0) {
        for (let r = 0; r < m; r++) y[r] = rhs[r];
      } else {
        const vp = v - m;
        for (let r = 0; r < m; r++) {
          let s = rhs[v + r];
          for (let k = 0; k < m; k++) s -= A[o + r * m + k] * y[vp + k];
          y[v + r] = s;
        }
      }
      luSolveAt(lu, o, piv, v, m, y, v);
    }
    // Backward: x_i = y_i − C'_i x_{i+1}
    const last = (n - 1) * m;
    for (let r = 0; r < m; r++) out[last + r] = y[last + r];
    for (let i = n - 2; i >= 0; i--) {
      const o = i * mm, v = i * m, vn = v + m;
      for (let r = 0; r < m; r++) {
        let s = y[v + r];
        for (let k = 0; k < m; k++) s -= cp[o + r * m + k] * out[vn + k];
        out[v + r] = s;
      }
    }
    return out;
  }

  /** out = M·x using the stored (unfactored) blocks. For residual checks and tests. */
  multiply(x, out = new Float64Array(this.n * this.m)) {
    const { n, m, A, B, C } = this;
    const mm = m * m;
    for (let i = 0; i < n; i++) {
      const o = i * mm, v = i * m;
      for (let r = 0; r < m; r++) {
        let s = 0;
        for (let k = 0; k < m; k++) s += B[o + r * m + k] * x[v + k];
        if (i > 0) for (let k = 0; k < m; k++) s += A[o + r * m + k] * x[v - m + k];
        if (i < n - 1) for (let k = 0; k < m; k++) s += C[o + r * m + k] * x[v + m + k];
        out[v + r] = s;
      }
    }
    return out;
  }
}

// In-place LU with partial pivoting of the m×m block at a[o…]. Row swaps are recorded in
// piv[p…]. Returns min |pivot| / max |entry| (0 if singular).
function luFactor(a, o, piv, p, m) {
  let maxEntry = 0;
  for (let k = 0; k < m * m; k++) {
    const v = Math.abs(a[o + k]);
    if (v > maxEntry) maxEntry = v;
  }
  if (maxEntry === 0) return 0;
  let minPivot = Infinity;
  for (let k = 0; k < m; k++) {
    let best = k, bestVal = Math.abs(a[o + k * m + k]);
    for (let r = k + 1; r < m; r++) {
      const v = Math.abs(a[o + r * m + k]);
      if (v > bestVal) { best = r; bestVal = v; }
    }
    piv[p + k] = best;
    if (bestVal === 0) return 0;
    if (bestVal < minPivot) minPivot = bestVal;
    if (best !== k) {
      for (let c = 0; c < m; c++) {
        const t = a[o + k * m + c];
        a[o + k * m + c] = a[o + best * m + c];
        a[o + best * m + c] = t;
      }
    }
    const inv = 1 / a[o + k * m + k];
    for (let r = k + 1; r < m; r++) {
      const f = (a[o + r * m + k] *= inv);
      if (f !== 0) for (let c = k + 1; c < m; c++) a[o + r * m + c] -= f * a[o + k * m + c];
    }
  }
  return minPivot / maxEntry;
}

// Solve (LU) x = b in place, b in vec[0…m).
function luSolve(a, o, piv, p, m, vec) {
  luSolveAt(a, o, piv, p, m, vec, 0);
}

// Solve (LU) x = b in place, b in vec[s…s+m).
function luSolveAt(a, o, piv, p, m, vec, s) {
  for (let k = 0; k < m; k++) {
    const j = piv[p + k];
    if (j !== k) {
      const t = vec[s + k];
      vec[s + k] = vec[s + j];
      vec[s + j] = t;
    }
  }
  for (let r = 1; r < m; r++) {
    let v = vec[s + r];
    for (let c = 0; c < r; c++) v -= a[o + r * m + c] * vec[s + c];
    vec[s + r] = v;
  }
  for (let r = m - 1; r >= 0; r--) {
    let v = vec[s + r];
    for (let c = r + 1; c < m; c++) v -= a[o + r * m + c] * vec[s + c];
    vec[s + r] = v / a[o + r * m + r];
  }
}

/**
 * Complex block-tridiagonal solver, for small-signal (frequency-domain) problems. Same layout
 * and algorithm as BlockTridiagonal, with real and imaginary parts in separate arrays
 * (Ar/Ai, Br/Bi, Cr/Ci; right-hand sides and solutions as re/im pairs).
 */
export class ComplexBlockTridiagonal {
  constructor(n, m) {
    this.n = n;
    this.m = m;
    const N = n * m * m;
    for (const k of ['Ar', 'Ai', 'Br', 'Bi', 'Cr', 'Ci', '_lr', '_li', '_cr', '_ci']) this[k] = new Float64Array(N);
    this._piv = new Int32Array(n * m);
    this._yr = new Float64Array(n * m);
    this._yi = new Float64Array(n * m);
    this._tr = new Float64Array(m);
    this._ti = new Float64Array(m);
  }

  clear() {
    for (const k of ['Ar', 'Ai', 'Br', 'Bi', 'Cr', 'Ci']) this[k].fill(0);
  }

  factor() {
    const { n, m, Ar, Ai, Br, Bi, Cr, Ci, _lr: lr, _li: li, _cr: cr, _ci: ci, _piv: piv, _tr: tr, _ti: ti } = this;
    const mm = m * m;
    for (let i = 0; i < n; i++) {
      const o = i * mm;
      for (let r = 0; r < m; r++) {
        for (let c = 0; c < m; c++) {
          let sr = Br[o + r * m + c], si = Bi[o + r * m + c];
          if (i > 0) {
            for (let k = 0; k < m; k++) {
              const ar = Ar[o + r * m + k], ai = Ai[o + r * m + k];
              const pr = cr[o - mm + k * m + c], pi = ci[o - mm + k * m + c];
              sr -= ar * pr - ai * pi;
              si -= ar * pi + ai * pr;
            }
          }
          lr[o + r * m + c] = sr;
          li[o + r * m + c] = si;
        }
      }
      if (!cluFactor(lr, li, o, piv, i * m, m)) throw new Error(`ComplexBlockTridiagonal: diagonal block ${i} is singular`);
      if (i < n - 1) {
        for (let c = 0; c < m; c++) {
          for (let r = 0; r < m; r++) {
            tr[r] = Cr[o + r * m + c];
            ti[r] = Ci[o + r * m + c];
          }
          cluSolve(lr, li, o, piv, i * m, m, tr, ti, 0);
          for (let r = 0; r < m; r++) {
            cr[o + r * m + c] = tr[r];
            ci[o + r * m + c] = ti[r];
          }
        }
      }
    }
  }

  /** Solve for (re, im) right-hand sides; results written to outRe/outIm. */
  solve(bRe, bIm, outRe, outIm) {
    const { n, m, Ar, Ai, _lr: lr, _li: li, _cr: cr, _ci: ci, _piv: piv, _yr: yr, _yi: yi } = this;
    const mm = m * m;
    for (let i = 0; i < n; i++) {
      const o = i * mm, v = i * m;
      for (let r = 0; r < m; r++) {
        let sr = bRe[v + r], si = bIm[v + r];
        if (i > 0) {
          for (let k = 0; k < m; k++) {
            const ar = Ar[o + r * m + k], ai = Ai[o + r * m + k];
            sr -= ar * yr[v - m + k] - ai * yi[v - m + k];
            si -= ar * yi[v - m + k] + ai * yr[v - m + k];
          }
        }
        yr[v + r] = sr;
        yi[v + r] = si;
      }
      cluSolve(lr, li, o, piv, v, m, yr, yi, v);
    }
    const last = (n - 1) * m;
    for (let r = 0; r < m; r++) {
      outRe[last + r] = yr[last + r];
      outIm[last + r] = yi[last + r];
    }
    for (let i = n - 2; i >= 0; i--) {
      const o = i * mm, v = i * m;
      for (let r = 0; r < m; r++) {
        let sr = yr[v + r], si = yi[v + r];
        for (let k = 0; k < m; k++) {
          const pr = cr[o + r * m + k], pi = ci[o + r * m + k];
          const xr = outRe[v + m + k], xi = outIm[v + m + k];
          sr -= pr * xr - pi * xi;
          si -= pr * xi + pi * xr;
        }
        outRe[v + r] = sr;
        outIm[v + r] = si;
      }
    }
  }
}

// Complex LU with partial pivoting (by modulus) of the block at o; false if singular.
function cluFactor(ar, ai, o, piv, p, m) {
  for (let k = 0; k < m; k++) {
    let best = k, bestVal = Math.hypot(ar[o + k * m + k], ai[o + k * m + k]);
    for (let r = k + 1; r < m; r++) {
      const v = Math.hypot(ar[o + r * m + k], ai[o + r * m + k]);
      if (v > bestVal) {
        best = r;
        bestVal = v;
      }
    }
    piv[p + k] = best;
    if (bestVal === 0) return false;
    if (best !== k) {
      for (let c = 0; c < m; c++) {
        let t = ar[o + k * m + c];
        ar[o + k * m + c] = ar[o + best * m + c];
        ar[o + best * m + c] = t;
        t = ai[o + k * m + c];
        ai[o + k * m + c] = ai[o + best * m + c];
        ai[o + best * m + c] = t;
      }
    }
    const dr = ar[o + k * m + k], di = ai[o + k * m + k], d2 = dr * dr + di * di;
    const invr = dr / d2, invi = -di / d2;
    for (let r = k + 1; r < m; r++) {
      const xr = ar[o + r * m + k], xi = ai[o + r * m + k];
      const fr = xr * invr - xi * invi, fi = xr * invi + xi * invr;
      ar[o + r * m + k] = fr;
      ai[o + r * m + k] = fi;
      if (fr === 0 && fi === 0) continue;
      for (let c = k + 1; c < m; c++) {
        const ur = ar[o + k * m + c], ui = ai[o + k * m + c];
        ar[o + r * m + c] -= fr * ur - fi * ui;
        ai[o + r * m + c] -= fr * ui + fi * ur;
      }
    }
  }
  return true;
}

function cluSolve(ar, ai, o, piv, p, m, vr, vi, s) {
  for (let k = 0; k < m; k++) {
    const j = piv[p + k];
    if (j !== k) {
      let t = vr[s + k];
      vr[s + k] = vr[s + j];
      vr[s + j] = t;
      t = vi[s + k];
      vi[s + k] = vi[s + j];
      vi[s + j] = t;
    }
  }
  for (let r = 1; r < m; r++) {
    let xr = vr[s + r], xi = vi[s + r];
    for (let c = 0; c < r; c++) {
      const lr = ar[o + r * m + c], li = ai[o + r * m + c];
      xr -= lr * vr[s + c] - li * vi[s + c];
      xi -= lr * vi[s + c] + li * vr[s + c];
    }
    vr[s + r] = xr;
    vi[s + r] = xi;
  }
  for (let r = m - 1; r >= 0; r--) {
    let xr = vr[s + r], xi = vi[s + r];
    for (let c = r + 1; c < m; c++) {
      const ur = ar[o + r * m + c], ui = ai[o + r * m + c];
      xr -= ur * vr[s + c] - ui * vi[s + c];
      xi -= ur * vi[s + c] + ui * vr[s + c];
    }
    const dr = ar[o + r * m + r], di = ai[o + r * m + r], d2 = dr * dr + di * di;
    vr[s + r] = (xr * dr + xi * di) / d2;
    vi[s + r] = (xi * dr - xr * di) / d2;
  }
}
