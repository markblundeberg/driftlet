// Block-tridiagonal linear solver (block Thomas algorithm).
//
// Solves, for node i = 0 … n−1 with blocks of size m_i,
//     A_i x_{i−1} + B_i x_i + C_i x_{i+1} = d_i
// where B_i is m_i × m_i, A_i is m_i × m_{i−1} and C_i is m_i × m_{i+1}. Blocks are stored
// row-major, one after another, in flat Float64Arrays; `offA[i]`, `offB[i]`, `offC[i]` give
// where each starts and `offX[i]` where node i's entries start in a vector. A_0 and C_{n−1} are
// ignored, but stored as if the end blocks repeated (m_{−1} = m_0, m_n = m_{n−1}), so with a
// single size m every block i sits at i·m², as for a uniform array of blocks. A size may be 0
// (the chain then splits there). All storage is allocated once; factor() and solve() allocate
// nothing.
//
// Within each diagonal block, LU uses partial pivoting, so zero diagonal entries (such as a
// slot that carries an interface flux instead of a value) are fine. There is no pivoting
// *between* blocks: that is the usual block Thomas trade-off, safe for the diagonally
// dominant-ish systems that finite-volume assembly produces. What it can't survive is a mode
// held only by the blocks after it: eliminated from the left, a closed strictly neutral region's
// φ is free until the face beyond it (the face's unknowns are fluxes, so the blocks before it
// alone leave the level undetermined), and a whole pivot column cancels to round-off there. With
// `bothWays`, a factorisation from the left that loses a column that way is done from the right
// too (a copy of the system with its blocks in reverse order), and that one kept where it loses
// nothing: a device and its mirror image are then factorised the same way round (see factor).

// Offsets of every block and vector segment for block sizes `sizes`.
function layout(n, m, name) {
  if (!(Number.isInteger(n) && n >= 1)) throw new RangeError(`n must be a positive integer, got ${n}`);
  let sizes;
  if (typeof m === 'number') {
    if (!(Number.isInteger(m) && m >= 1)) throw new RangeError(`m must be a positive integer, got ${m}`);
    sizes = new Int32Array(n).fill(m);
  } else {
    if (!(m && m.length === n)) throw new RangeError(`${name}: give one block size per node (${n})`);
    sizes = Int32Array.from(m);
    for (const v of sizes) if (!(v >= 0)) throw new RangeError(`block sizes must be non-negative integers, got ${v}`);
  }
  const offA = new Int32Array(n + 1), offB = new Int32Array(n + 1), offC = new Int32Array(n + 1), offX = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) {
    const mi = sizes[i], mp = i > 0 ? sizes[i - 1] : mi, mn = i < n - 1 ? sizes[i + 1] : mi;
    offA[i + 1] = offA[i] + mi * mp;
    offB[i + 1] = offB[i] + mi * mi;
    offC[i + 1] = offC[i] + mi * mn;
    offX[i + 1] = offX[i] + mi;
  }
  const maxM = sizes.reduce((a, v) => Math.max(a, v), 0);
  return { sizes, offA, offB, offC, offX, maxM, uniform: typeof m === 'number' ? m : null };
}

export class BlockTridiagonal {
  /**
   * @param {number} n number of block rows (nodes)
   * @param {number | ArrayLike<number>} m block size (unknowns per node), or one size per node
   */
  constructor(n, m) {
    const L = layout(n, m, 'BlockTridiagonal');
    this.n = n;
    /** The block size if every block has the same size, else null. */
    this.m = L.uniform;
    /** Block size of each node. */
    this.sizes = L.sizes;
    this.offA = L.offA;
    this.offB = L.offB;
    this.offC = L.offC;
    this.offX = L.offX;
    /** Total number of unknowns. */
    this.size = L.offX[n];
    /** Sub-diagonal blocks A_i (coupling to node i−1). */
    this.A = new Float64Array(L.offA[n]);
    /** Diagonal blocks B_i. */
    this.B = new Float64Array(L.offB[n]);
    /** Super-diagonal blocks C_i (coupling to node i+1). */
    this.C = new Float64Array(L.offC[n]);

    this._lu = new Float64Array(L.offB[n]); // LU factors of W_i = B_i − A_i C'_{i−1}
    this._piv = new Int32Array(this.size);
    this._cp = new Float64Array(L.offC[n]); // C'_i = W_i⁻¹ C_i
    this._y = new Float64Array(this.size);
    this._tmp = new Float64Array(L.maxM);
    this.factored = false;
    /** Smallest |pivot| / (largest |entry| of its block) seen in the last factor(). */
    this.minPivotRatio = Infinity;
    /** The block found exactly singular by the last factor(), or −1. */
    this.singularBlock = -1;
    /**
     * Static pivoting: an exactly zero pivot (all its candidates cancelled to nothing) is
     * replaced by 1e-15 of its block's largest entry rather than thrown on, for a caller that
     * refines its solves (and counted in `perturbedPivots`).
     */
    this.staticPivots = false;
    /** How many pivots the last factor() replaced. */
    this.perturbedPivots = 0;
    /**
     * Where the factorisation from the left loses a column to round-off (or a block is
     * singular), factorise from the right as well, and keep that one where it loses nothing (see
     * factor; with staticPivots, which is where columns are looked at).
     */
    this.bothWays = false;
    /** How many of the last factor()'s pivot columns cancelled to round-off (see factor). */
    this.lostColumns = 0;
    /** Whether the last factor() kept the factorisation from the right (solve() then uses it). */
    this.fromRight = false;
    this._rev = null; // (the system with its blocks in reverse order, made when first needed)
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
   * Throws if a diagonal block is exactly singular (unless `staticPivots`).
   * With `bothWays`, where the factorisation from the left lost a column (lostColumns: a
   * pivot's every candidate cancelled to 1e-15 of the magnitudes it was formed from, so that the
   * blocks up to it are singular to working precision) or met a singular block, the system is
   * factorised from the right too, and that one kept (fromRight) where it loses nothing: no
   * column, no perturbed pivot, every pivot ratio at least 1e-12 (or the other failed).
   * minPivotRatio, perturbedPivots, lostColumns and singularBlock are then the kept one's (its
   * blocks numbered as here). A pivot that's small but formed from small entries (a population
   * held weakly, which no direction holds better) leaves the left one standing, its perturbed
   * pivots turning on the refined solves (src/newton.js) such a population needs: kept by pivot
   * ratio alone, the right one had a floating MOS gate's inversion layer drift 9e-8 V from rest
   * over 1e7 s.
   */
  factor() {
    this.fromRight = false;
    let failed = null;
    try {
      this._factorFromLeft();
    } catch (err) {
      if (!this.bothWays) throw err;
      failed = err;
    }
    if (!this.bothWays || (!failed && !(this.lostColumns > 0))) return;
    // Eliminated from the left, a mode that only the blocks after it hold (a closed region's
    // level behind a face, whose unknowns are fluxes) leaves a pivot of round-off, which then
    // scales all that follows: a device solved, its mirror image
    // refused (round 9's mirror property: a reacting solution with its bath on the right, pivot
    // ratio 3e-40 where the other way round has 7e-4, and Newton's first update 1e29 thermal
    // units).
    const rev = this._reversed();
    try {
      rev.factor();
    } catch (err) {
      if (failed) throw failed;
      return;
    }
    if (!failed && !(rev.lostColumns === 0 && rev.perturbedPivots === 0 && rev.minPivotRatio >= 1e-12)) return;
    this.fromRight = true;
    this.factored = true;
    this.minPivotRatio = rev.minPivotRatio;
    this.perturbedPivots = rev.perturbedPivots;
    this.lostColumns = rev.lostColumns;
    this.singularBlock = -1;
  }

  // The same system with its blocks in reverse order (block i as n − 1 − i: B as it is, A and C
  // swapped), copied in; its storage made on first use.
  _reversed() {
    const { n, A, B, C, sizes, offA, offB, offC } = this;
    const rev = (this._rev ??= new BlockTridiagonal(n, this.m ?? Array.from(sizes).reverse()));
    rev.staticPivots = this.staticPivots;
    for (let i = 0; i < n; i++) {
      const j = n - 1 - i;
      rev.B.set(B.subarray(offB[i], offB[i + 1]), rev.offB[j]);
      rev.A.set(C.subarray(offC[i], offC[i + 1]), rev.offA[j]); // (m_i × m_{i+1}: rev's A at j)
      rev.C.set(A.subarray(offA[i], offA[i + 1]), rev.offC[j]); // (m_i × m_{i−1}: rev's C at j)
    }
    return rev;
  }

  // From the left: the block Thomas factorisation itself.
  _factorFromLeft() {
    const { n, A, B, C, sizes, offA, offB, offC, offX } = this;
    const lu = this._lu, piv = this._piv, cp = this._cp, tmp = this._tmp;
    let minRatio = Infinity;
    this.singularBlock = -1;
    this.perturbedPivots = 0;
    this.lostColumns = 0;

    for (let i = 0; i < n; i++) {
      const m = sizes[i], o = offB[i];
      if (m === 0) continue;
      // W_i = B_i − A_i C'_{i−1}
      if (i === 0 || sizes[i - 1] === 0) {
        for (let k = 0; k < m * m; k++) lu[o + k] = B[o + k];
      } else {
        // Row by row, skipping A's zeros (a gate or a local unknown has none to its neighbours):
        // the same sums, subtracted in the same order.
        const mp = sizes[i - 1], oa = offA[i], op = offC[i - 1]; // C'_{i−1}: mp × m
        for (let r = 0; r < m; r++) {
          const or = o + r * m;
          for (let c = 0; c < m; c++) lu[or + c] = B[or + c];
          for (let k = 0; k < mp; k++) {
            const a = A[oa + r * mp + k];
            if (a === 0) continue;
            const ok = op + k * m;
            for (let c = 0; c < m; c++) lu[or + c] -= a * cp[ok + c];
          }
        }
      }
      const perturbedBefore = this.perturbedPivots;
      let ratio = luFactor(lu, o, piv, offX[i], m, this);
      // A pivot far below the entries it's formed from may be round-off alone (the level of a
      // region that only the blocks after it hold): its value, 7e-32 or 3e-16 or 0 as the
      // unknowns happen to be listed, then scales all that follows (by 1e31, and the next block
      // lost its own pivot outright). Redone with each entry's magnitudes carried (see
      // cancellation), a pivot 1e-18 of them, all its digits lost, is perturbed like one that
      // cancelled to exactly nothing. (At 1e-15 it can still hold a digit Newton needs: a port's
      // steady solve wandered instead of converging.)
      // (the block redone from scratch: its pivots counted there, not twice; Neon's review, 3.1)
      if (this.staticPivots && ratio < 1e-12) {
        this.perturbedPivots = perturbedBefore;
        ratio = this._noisyPivots(i);
      }
      if (!(ratio > 0)) {
        this.factored = false;
        this.singularBlock = i;
        throw new Error(`BlockTridiagonal: diagonal block ${i} is singular`);
      }
      if (ratio < minRatio) minRatio = ratio;
      // C'_i = W_i⁻¹ C_i, column by column
      if (i < n - 1) {
        const mn = sizes[i + 1], oc = offC[i];
        for (let c = 0; c < mn; c++) {
          let any = false;
          for (let r = 0; r < m; r++) if ((tmp[r] = C[oc + r * mn + c]) !== 0) any = true;
          if (any) luSolve(lu, o, piv, offX[i], m, tmp); // (a zero column stays zero)
          for (let r = 0; r < m; r++) cp[oc + r * mn + c] = tmp[r];
        }
      }
    }
    this.minPivotRatio = minRatio;
    this.factored = true;
  }

  // Block i again, its W and LU with the magnitudes each entry is formed from (as cancellation()
  // does), perturbing pivots 1e-18 of theirs. Returns the block's pivot ratio.
  _noisyPivots(i) {
    const { A, B, sizes, offA, offB, offC, offX } = this, cp = this._cp, lu = this._lu, piv = this._piv;
    const m = sizes[i], mp = i > 0 ? sizes[i - 1] : 0, oa = offA[i], o = offB[i], oc = i > 0 ? offC[i - 1] : 0;
    const S = (this._mag ??= new Float64Array(this._tmp.length * this._tmp.length));
    for (let r = 0; r < m; r++) {
      for (let c = 0; c < m; c++) {
        let v = B[o + r * m + c], sum = Math.abs(v);
        for (let k = 0; k < mp; k++) {
          const a = A[oa + r * mp + k];
          if (a === 0) continue;
          const t = a * cp[oc + k * m + c];
          v -= t;
          sum += Math.abs(t);
        }
        lu[o + r * m + c] = v;
        S[r * m + c] = sum;
      }
    }
    let maxEntry = 0;
    for (let k = 0; k < m * m; k++) maxEntry = Math.max(maxEntry, Math.abs(lu[o + k]));
    if (maxEntry === 0) return 0;
    let minPivot = Infinity;
    for (let k = 0; k < m; k++) {
      let best = k, bestVal = Math.abs(lu[o + k * m + k]), col = S[k * m + k];
      for (let r = k + 1; r < m; r++) {
        const v = Math.abs(lu[o + r * m + k]);
        if (v > bestVal) { best = r; bestVal = v; }
        if (S[r * m + k] > col) col = S[r * m + k];
      }
      // (every candidate round-off against what the column was formed from: the blocks so far
      // singular, whichever row the pivot is taken from; see factor. Counted from 1e-15, not the
      // 1e-18 below which a pivot is perturbed: a closed compartment's level lost to round-off
      // left pivots 2e-17 to 6e-15 of their terms, solves of backward error up to 0.24, and the
      // other way round was never tried, a mirrored membrane cell taking 112 factorisations for
      // the 104 written the other way (Neon's review, 3.1). Not higher: from 1e-14 a floating MOS
      // gate's inversion layer, a population held weakly, counts as lost, and the right one kept
      // lets it drift; see factor.)
      if (bestVal <= 1e-15 * col) this.lostColumns++;
      piv[offX[i] + k] = best;
      if (best !== k) {
        for (let c = 0; c < m; c++) {
          let t = lu[o + k * m + c];
          lu[o + k * m + c] = lu[o + best * m + c];
          lu[o + best * m + c] = t;
          t = S[k * m + c];
          S[k * m + c] = S[best * m + c];
          S[best * m + c] = t;
        }
      }
      const sum = S[k * m + k];
      if (bestVal <= 1e-18 * sum || bestVal === 0) {
        const size = sum > 0 ? 1e-15 * sum : 1e-15 * maxEntry;
        bestVal = size;
        lu[o + k * m + k] = lu[o + k * m + k] < 0 ? -size : size;
        this.perturbedPivots++;
      }
      if (bestVal < minPivot) minPivot = bestVal;
      const inv = 1 / lu[o + k * m + k];
      for (let r = k + 1; r < m; r++) {
        const f = (lu[o + r * m + k] *= inv);
        if (f !== 0) {
          const af = Math.abs(f);
          for (let c = k + 1; c < m; c++) {
            lu[o + r * m + c] -= f * lu[o + k * m + c];
            S[r * m + c] += af * S[k * m + c];
          }
        }
      }
    }
    return minPivot / maxEntry;
  }

  /**
   * Solve with the current factorisation. `rhs` and `out` are Float64Arrays of length `size`
   * (node-major); they may be the same array.
   */
  solve(rhs, out = new Float64Array(this.size)) {
    if (!this.factored) throw new Error('BlockTridiagonal: call factor() before solve()');
    if (this.fromRight) {
      // (through the reversed copy: rhs into its order, its solution back into ours)
      const rev = this._rev, { n, offX } = this, v = (this._revVec ??= new Float64Array(this.size));
      for (let i = 0; i < n; i++) v.set(rhs.subarray(offX[i], offX[i + 1]), rev.offX[n - 1 - i]);
      rev.solve(v, v);
      for (let i = 0; i < n; i++) out.set(v.subarray(rev.offX[n - 1 - i], rev.offX[n - i]), offX[i]);
      return out;
    }
    const { n, A, sizes, offA, offB, offC, offX } = this;
    const lu = this._lu, piv = this._piv, cp = this._cp, y = this._y;

    // Forward: y_i = W_i⁻¹ (d_i − A_i y_{i−1})
    for (let i = 0; i < n; i++) {
      const m = sizes[i], v = offX[i];
      if (m === 0) continue;
      const mp = i > 0 ? sizes[i - 1] : 0, oa = offA[i], vp = v - mp;
      for (let r = 0; r < m; r++) {
        let s = rhs[v + r];
        for (let k = 0; k < mp; k++) s -= A[oa + r * mp + k] * y[vp + k];
        y[v + r] = s;
      }
      luSolveAt(lu, offB[i], piv, v, m, y, v);
    }
    // Backward: x_i = y_i − C'_i x_{i+1}
    for (let i = n - 1; i >= 0; i--) {
      const m = sizes[i], v = offX[i];
      const mn = i < n - 1 ? sizes[i + 1] : 0, oc = offC[i], vn = v + m;
      for (let r = 0; r < m; r++) {
        let s = y[v + r];
        for (let k = 0; k < mn; k++) s -= cp[oc + r * mn + k] * out[vn + k];
        out[v + r] = s;
      }
    }
    return out;
  }

  /**
   * After factor(): the digits lost to cancellation in the factorisation, at worst, and the block
   * where that happens. A running error bound: alongside each entry of W_i = B_i − A_i C′_{i−1}
   * and of its LU (with the same partial pivoting), the sum of the magnitudes it was formed from
   * is carried, and each pivot is compared with its sum: log₁₀(Σ|terms| / |pivot|). Many digits
   * lost means a nearly singular system, such as a stiff chain held only weakly at its ends,
   * whose level then rests on the lost digits. Costs about as much as a factorisation. After a
   * factorisation that failed on a singular block, it covers the blocks up to that one.
   * @returns {{ digits: number, block: number }}
   */
  cancellation() {
    if (this.fromRight && this.factored) {
      const { digits, block } = this._rev.cancellation();
      return { digits, block: block < 0 ? -1 : this.n - 1 - block };
    }
    if (!this.factored && !(this.singularBlock >= 0)) throw new Error('BlockTridiagonal: call factor() before cancellation()');
    const { A, B, sizes, offA, offB, offC } = this, cp = this._cp;
    const n = this.factored ? this.n : this.singularBlock + 1;
    const mx = sizes.reduce((a, v) => Math.max(a, v), 0);
    const w = new Float64Array(mx * mx), S = new Float64Array(mx * mx);
    let worst = 1, block = -1;
    for (let i = 0; i < n; i++) {
      const m = sizes[i], mp = i > 0 ? sizes[i - 1] : 0, oa = offA[i], ob = offB[i], oc = i > 0 ? offC[i - 1] : 0;
      for (let r = 0; r < m; r++) {
        for (let c = 0; c < m; c++) {
          let v = B[ob + r * m + c], sum = Math.abs(v);
          for (let k = 0; k < mp; k++) {
            const t = A[oa + r * mp + k] * cp[oc + k * m + c];
            v -= t;
            sum += Math.abs(t);
          }
          w[r * m + c] = v;
          S[r * m + c] = sum;
        }
      }
      for (let k = 0; k < m; k++) {
        let p = k;
        for (let r = k + 1; r < m; r++) if (Math.abs(w[r * m + k]) > Math.abs(w[p * m + k])) p = r;
        if (p !== k) {
          for (let c = 0; c < m; c++) {
            [w[k * m + c], w[p * m + c]] = [w[p * m + c], w[k * m + c]];
            [S[k * m + c], S[p * m + c]] = [S[p * m + c], S[k * m + c]];
          }
        }
        const piv = Math.abs(w[k * m + k]), sum = S[k * m + k];
        if (sum > worst * piv) {
          worst = piv === 0 ? Infinity : sum / piv;
          block = i;
        }
        if (piv === 0) continue;
        for (let r = k + 1; r < m; r++) {
          const f = w[r * m + k] / w[k * m + k];
          if (f === 0) continue;
          for (let c = k + 1; c < m; c++) {
            w[r * m + c] -= f * w[k * m + c];
            S[r * m + c] += Math.abs(f) * S[k * m + c];
          }
        }
      }
    }
    return { digits: Math.log10(worst), block };
  }

  /** out = M·x using the stored (unfactored) blocks. For residual checks and tests. */
  multiply(x, out = new Float64Array(this.size)) {
    const { n, A, B, C, sizes, offA, offB, offC, offX } = this;
    for (let i = 0; i < n; i++) {
      const m = sizes[i], v = offX[i];
      const mp = i > 0 ? sizes[i - 1] : 0, mn = i < n - 1 ? sizes[i + 1] : 0;
      for (let r = 0; r < m; r++) {
        let s = 0;
        for (let k = 0; k < m; k++) s += B[offB[i] + r * m + k] * x[v + k];
        for (let k = 0; k < mp; k++) s += A[offA[i] + r * mp + k] * x[v - mp + k];
        for (let k = 0; k < mn; k++) s += C[offC[i] + r * mn + k] * x[v + m + k];
        out[v + r] = s;
      }
    }
    return out;
  }
}

// In-place LU with partial pivoting of the m×m block at a[o…]. Row swaps are recorded in
// piv[p…]. Returns min |pivot| / max |entry| (0 if singular).
function luFactor(a, o, piv, p, m, perturb = null) {
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
    if (bestVal === 0) {
      if (!perturb?.staticPivots) return 0;
      bestVal = a[o + best * m + k] = 1e-15 * maxEntry;
      perturb.perturbedPivots++;
    }
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
 * and algorithm as BlockTridiagonal (including per-node block sizes), with real and imaginary
 * parts in separate arrays (Ar/Ai, Br/Bi, Cr/Ci; right-hand sides and solutions as re/im
 * pairs).
 */
export class ComplexBlockTridiagonal {
  /**
   * @param {number} n number of block rows (nodes)
   * @param {number | ArrayLike<number>} m block size, or one size per node
   */
  constructor(n, m) {
    const L = layout(n, m, 'ComplexBlockTridiagonal');
    this.n = n;
    this.m = L.uniform;
    this.sizes = L.sizes;
    this.offA = L.offA;
    this.offB = L.offB;
    this.offC = L.offC;
    this.offX = L.offX;
    this.size = L.offX[n];
    for (const k of ['Ar', 'Ai']) this[k] = new Float64Array(L.offA[n]);
    for (const k of ['Br', 'Bi', '_lr', '_li']) this[k] = new Float64Array(L.offB[n]);
    for (const k of ['Cr', 'Ci', '_cr', '_ci']) this[k] = new Float64Array(L.offC[n]);
    this._piv = new Int32Array(this.size);
    this._yr = new Float64Array(this.size);
    this._yi = new Float64Array(this.size);
    this._tr = new Float64Array(L.maxM);
    this._ti = new Float64Array(L.maxM);
  }

  clear() {
    for (const k of ['Ar', 'Ai', 'Br', 'Bi', 'Cr', 'Ci']) this[k].fill(0);
  }

  factor() {
    const { n, Ar, Ai, Br, Bi, Cr, Ci, _lr: lr, _li: li, _cr: cr, _ci: ci, _piv: piv, _tr: tr, _ti: ti } = this;
    const { sizes, offA, offB, offC, offX } = this;
    for (let i = 0; i < n; i++) {
      const m = sizes[i], o = offB[i];
      if (m === 0) continue;
      const mp = i > 0 ? sizes[i - 1] : 0, oa = offA[i], op = i > 0 ? offC[i - 1] : 0;
      for (let r = 0; r < m; r++) {
        for (let c = 0; c < m; c++) {
          let sr = Br[o + r * m + c], si = Bi[o + r * m + c];
          for (let k = 0; k < mp; k++) {
            const ar = Ar[oa + r * mp + k], ai = Ai[oa + r * mp + k];
            const pr = cr[op + k * m + c], pi = ci[op + k * m + c];
            sr -= ar * pr - ai * pi;
            si -= ar * pi + ai * pr;
          }
          lr[o + r * m + c] = sr;
          li[o + r * m + c] = si;
        }
      }
      if (!cluFactor(lr, li, o, piv, offX[i], m)) throw new Error(`ComplexBlockTridiagonal: diagonal block ${i} is singular`);
      if (i < n - 1) {
        const mn = sizes[i + 1], oc = offC[i];
        for (let c = 0; c < mn; c++) {
          for (let r = 0; r < m; r++) {
            tr[r] = Cr[oc + r * mn + c];
            ti[r] = Ci[oc + r * mn + c];
          }
          cluSolve(lr, li, o, piv, offX[i], m, tr, ti, 0);
          for (let r = 0; r < m; r++) {
            cr[oc + r * mn + c] = tr[r];
            ci[oc + r * mn + c] = ti[r];
          }
        }
      }
    }
  }

  /** Solve for (re, im) right-hand sides; results written to outRe/outIm. */
  solve(bRe, bIm, outRe, outIm) {
    const { n, Ar, Ai, _lr: lr, _li: li, _cr: cr, _ci: ci, _piv: piv, _yr: yr, _yi: yi } = this;
    const { sizes, offA, offB, offC, offX } = this;
    for (let i = 0; i < n; i++) {
      const m = sizes[i], v = offX[i];
      if (m === 0) continue;
      const mp = i > 0 ? sizes[i - 1] : 0, oa = offA[i], vp = v - mp;
      for (let r = 0; r < m; r++) {
        let sr = bRe[v + r], si = bIm[v + r];
        for (let k = 0; k < mp; k++) {
          const ar = Ar[oa + r * mp + k], ai = Ai[oa + r * mp + k];
          sr -= ar * yr[vp + k] - ai * yi[vp + k];
          si -= ar * yi[vp + k] + ai * yr[vp + k];
        }
        yr[v + r] = sr;
        yi[v + r] = si;
      }
      cluSolve(lr, li, offB[i], piv, v, m, yr, yi, v);
    }
    for (let i = n - 1; i >= 0; i--) {
      const m = sizes[i], v = offX[i];
      const mn = i < n - 1 ? sizes[i + 1] : 0, oc = offC[i], vn = v + m;
      for (let r = 0; r < m; r++) {
        let sr = yr[v + r], si = yi[v + r];
        for (let k = 0; k < mn; k++) {
          const pr = cr[oc + r * mn + k], pi = ci[oc + r * mn + k];
          const xr = outRe[vn + k], xi = outIm[vn + k];
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
