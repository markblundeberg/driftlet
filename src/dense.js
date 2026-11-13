// Small dense linear systems: Gaussian elimination with partial pivoting, the one routine the
// statistics models' Newton, the junction helper and the polarization fit share. (Three copies,
// with three conventions for A, before: rows, columns and flat; a visitor copying the nearest one
// got the transpose wrong. lawRows.js keeps its own, which pivots on zero; impedance.js its
// complex one.)
//
// solveFlat(A, b, x, m): A row-major in a flat array (m × m, its first m·m entries read), b's
// first m entries; x (length ≥ m) gets the solution. A and b are left as they were. A zero
// pivot gives a non-finite x, which the callers' convergence checks see.
export function solveFlat(A, b, x, m) {
  const M = Float64Array.from(A.subarray ? A.subarray(0, m * m) : A.slice(0, m * m)), v = Float64Array.from(b.subarray ? b.subarray(0, m) : b.slice(0, m));
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
  return x;
}

// A given as rows (A[i][j]), or with byColumns as columns (A[j][i]): the solution as a new array.
export function solveDense(A, b, { byColumns = false } = {}) {
  const m = b.length, flat = new Float64Array(m * m);
  for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) flat[i * m + j] = byColumns ? A[j][i] : A[i][j];
  return Array.from(solveFlat(flat, Float64Array.from(b), new Float64Array(m), m));
}
