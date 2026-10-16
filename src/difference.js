// The dilute kernels' Jacobian terms kept in difference form, for an exact J·v.
//
// Assembled, a Scharfetter–Gummel flux's dependence on η at its two ends is two entries, G and
// −G(1 + E), that cancel for a nearly uniform v; where G is huge (an inversion layer's 0.1 nm
// cells), J·v loses the flux to round-off. Kept as terms, each is evaluated on differences first:
//   δN = a·(δη_L − z δφ̂_L) + b·(δη_R − δη_L) + d·(δφ̂_R − δφ̂_L),
// with a = −G E (through c_L), b = −G(1 + E) and d = z ∂N/∂Δ (through the drift). A node's
// storage and space charge depend on δη_i − z_i δφ̂ alone; the displacement on δφ̂_R − δφ̂_L. A
// uniform shift of a region (z_i s on each η_i, s on φ̂) then gives exactly zero, as it must.
//
// Non-ideal statistics couple the species: c_i depends on every ζ_j, so a node's storage of i
// has a term in each δη_j − z_j δφ̂ (K_ij), and a flux has side terms in those at either end (its
// c_L and the excess potentials in its drift). Each is a term in the same differences, so the
// shift still gives exactly zero.
//
// Indices are compact (the solver's rix), with the sink (index N) for an unknown that isn't one:
// vectors passed in are read there as zero.
export class DifferenceTerms {
  constructor(N, nodes, segments, sides = 0) {
    this.N = N;
    this.dead = new Uint8Array(N + 1); // rows replaced since (a contact's level, a held stretch)
    this.nR = new Int32Array(nodes); // node term: its balance row, η column, φ̂ column (= Gauss row)
    this.nI = new Int32Array(nodes);
    this.nP = new Int32Array(nodes);
    this.nZ = new Float64Array(nodes);
    this.nMass = new Float64Array(nodes); // storage per unit (δη − z δφ̂), times 1/dt (or iω)
    this.nCharge = new Float64Array(nodes); // into the Gauss row, per unit (δη − z δφ̂)
    this.sI = new Int32Array(4 * segments); // segment term: η_L, η_R, φ̂_L, φ̂_R
    this.sZ = new Float64Array(segments);
    this.sC = new Float64Array(3 * segments); // a, b, d
    this.dI = new Int32Array(2 * segments); // displacement: φ̂_L, φ̂_R
    this.dK = new Float64Array(segments);
    this.xI = new Int32Array(4 * sides); // a flux's side term: its two rows, an η column, its φ̂ column
    this.xZ = new Float64Array(sides);
    this.xC = new Float64Array(sides);
    this.combos = []; // [rk, ri, w]: row rk += w × row ri (strictly neutral nodes' charge rows)
    this.t = new Float64Array(N + 1);
    this.x = new Float64Array(N + 1);
    this.reset();
  }

  // Afresh for an assembly; `rate` multiplies storage in captured rows (1/dt).
  reset(rate = 0) {
    this.rate = rate;
    this.nn = this.ns = this.nd = this.nx = 0;
    this.dead.fill(0);
    this.combos.length = 0;
  }

  // Storage in row r (mass) and charge in Gauss row p, per unit δη_i − z δφ̂ (row r = i, unless
  // non-ideal statistics couple the species).
  node(i, p, z, mass, charge, r = i) {
    const k = this.nn++;
    this.nR[k] = r;
    this.nI[k] = i;
    this.nP[k] = p;
    this.nZ[k] = z;
    this.nMass[k] = mass;
    this.nCharge[k] = charge;
  }

  segment(iL, iR, pL, pR, z, a, b, d) {
    const k = this.ns++;
    this.sI[4 * k] = iL;
    this.sI[4 * k + 1] = iR;
    this.sI[4 * k + 2] = pL;
    this.sI[4 * k + 3] = pR;
    this.sZ[k] = z;
    this.sC[3 * k] = a;
    this.sC[3 * k + 1] = b;
    this.sC[3 * k + 2] = d;
  }

  // A flux from row rL to row rR, per unit δη_i − z δφ̂ at one end.
  side(rL, rR, i, p, z, w) {
    const k = this.nx++;
    this.xI[4 * k] = rL;
    this.xI[4 * k + 1] = rR;
    this.xI[4 * k + 2] = i;
    this.xI[4 * k + 3] = p;
    this.xZ[k] = z;
    this.xC[k] = w;
  }

  displacement(pL, pR, k) {
    const q = this.nd++;
    this.dI[2 * q] = pL;
    this.dI[2 * q + 1] = pR;
    this.dK[q] = k;
  }

  kill(row) {
    this.dead[row] = 1;
  }

  combine(rk, ri, w) {
    this.combos.push([rk, ri, w]);
  }

  // into[col] += w × row's coefficients, the steady ones plus `mass` times the storage's (a row
  // captured as plain entries: a contact's current, read before its row is replaced).
  capture(row, w, mass, into) {
    const N = this.N, add = (col, v) => {
      if (col < N && v !== 0) into[col] += w * v;
    };
    for (let k = 0; k < this.nn; k++) {
      const i = this.nI[k], p = this.nP[k], z = this.nZ[k];
      let f = 0;
      if (row === this.nR[k]) f += mass * this.nMass[k];
      if (row === p) f += this.nCharge[k];
      if (f === 0) continue;
      add(i, f);
      add(p, -z * f);
    }
    for (let k = 0; k < this.ns; k++) {
      const iL = this.sI[4 * k], iR = this.sI[4 * k + 1];
      const s = row === iL ? 1 : row === iR ? -1 : 0;
      if (s === 0) continue;
      const z = this.sZ[k], a = this.sC[3 * k], b = this.sC[3 * k + 1], d = this.sC[3 * k + 2];
      add(iL, s * (a - b));
      add(iR, s * b);
      add(this.sI[4 * k + 2], s * (-z * a - d));
      add(this.sI[4 * k + 3], s * d);
    }
    for (let k = 0; k < this.nx; k++) {
      const s = row === this.xI[4 * k] ? 1 : row === this.xI[4 * k + 1] ? -1 : 0;
      if (s === 0) continue;
      add(this.xI[4 * k + 2], s * this.xC[k]);
      add(this.xI[4 * k + 3], -s * this.xZ[k] * this.xC[k]);
    }
    for (let q = 0; q < this.nd; q++) {
      const pL = this.dI[2 * q], pR = this.dI[2 * q + 1], s = row === pL ? 1 : row === pR ? -1 : 0;
      if (s === 0) continue;
      add(pL, s * this.dK[q]);
      add(pR, -s * this.dK[q]);
    }
  }

  // out += (steady terms) × v + mass × (storage terms) × v, over the live rows. v is read over
  // its first N entries; transformed (the solver's charge rows), it's mapped back by the caller.
  apply(v, out, steady = 1, mass = 0) {
    const { N, t, x } = this;
    for (let j = 0; j < N; j++) x[j] = v[j];
    x[N] = 0;
    t.fill(0);
    const { nR, nI, nP, nZ, nMass, nCharge } = this;
    for (let k = 0; k < this.nn; k++) {
      const i = nI[k], p = nP[k], w = x[i] - nZ[k] * x[p];
      t[nR[k]] += mass * nMass[k] * w;
      t[p] += steady * nCharge[k] * w;
    }
    if (steady !== 0) {
      const { sI, sZ, sC, dI, dK } = this;
      for (let k = 0; k < this.ns; k++) {
        const iL = sI[4 * k], iR = sI[4 * k + 1], pL = sI[4 * k + 2], pR = sI[4 * k + 3];
        const dN = steady * (sC[3 * k] * (x[iL] - sZ[k] * x[pL]) + sC[3 * k + 1] * (x[iR] - x[iL]) + sC[3 * k + 2] * (x[pR] - x[pL]));
        t[iL] += dN;
        t[iR] -= dN;
      }
      const { xI, xZ, xC } = this;
      for (let k = 0; k < this.nx; k++) {
        const dN = steady * xC[k] * (x[xI[4 * k + 2]] - xZ[k] * x[xI[4 * k + 3]]);
        t[xI[4 * k]] += dN;
        t[xI[4 * k + 1]] -= dN;
      }
      for (let q = 0; q < this.nd; q++) {
        const pL = dI[2 * q], pR = dI[2 * q + 1], dD = -steady * dK[q] * (x[pR] - x[pL]);
        t[pL] += dD;
        t[pR] -= dD;
      }
    }
    for (let j = 0; j < N; j++) if (this.dead[j]) t[j] = 0;
    for (const [rk, ri, w] of this.combos) t[rk] += w * t[ri];
    for (let j = 0; j < N; j++) out[j] += t[j];
    return out;
  }
}
