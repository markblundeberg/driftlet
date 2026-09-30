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
