// Liquid-junction potentials from formulas, for driftlet/kit: the diffusion potential φ_R − φ_L
// between two neutral solutions, as the classical theories give it. A device simulates the actual
// junction (steady between two baths, it is Planck's; stepping from a sharp boundary, a free
// diffusion junction); these are for comparison, and for quick estimates.

import { DeviceError } from './errors.js';
import { FARADAY, GAS_CONSTANT } from './constants.js';

const fail = (message) => {
  throw new DeviceError(message);
};

// The ions both solutions name, with their charge and diffusivity, and both concentrations.
function prepare(left, right, ions, what) {
  const names = [...new Set([...Object.keys(left), ...Object.keys(right)])];
  return names.map((name) => {
    const ion = ions[name];
    if (!ion || !Number.isFinite(ion.z) || !(ion.D > 0)) fail(`${what}: no { z, D } for '${name}' in the ions given`);
    const cL = left[name] ?? 0, cR = right[name] ?? 0;
    if (!(cL >= 0 && cR >= 0)) fail(`${what}: concentrations must be ≥ 0 (mol/m³), got ${name} ${cL} and ${cR}`);
    return { name, z: ion.z, D: ion.D, cL, cR };
  });
}
function neutral(list, side, what) {
  let q = 0, s = 0;
  for (const p of list) {
    q += p.z * p[side];
    s += Math.abs(p.z * p[side]);
  }
  if (Math.abs(q) > 1e-9 * s) fail(`${what}: the ${side === 'cL' ? 'left' : 'right'} solution carries a net charge (Σ z·c = ${q.toPrecision(3)} mol/m³)`);
}

/**
 * Henderson's liquid-junction potential, φ_R − φ_L (V): the concentrations taken to mix linearly
 * through the junction. Exact for a single binary salt, a close estimate otherwise (it's what
 * JPCalc and most electrophysiology corrections use).
 * @param {Record<string, number>} left mol/m³ by ion
 * @param {Record<string, number>} right mol/m³ by ion
 * @param {Record<string, { z: number, D: number }>} ions charge and diffusivity (m²/s), e.g. IONS
 * @param {{ T?: number }} [opts] K (default 298.15)
 */
export function henderson(left, right, ions, { T = 298.15 } = {}) {
  const list = prepare(left, right, ions, 'henderson');
  neutral(list, 'cL', 'henderson');
  neutral(list, 'cR', 'henderson');
  let a = 0, b = 0, gL = 0, gR = 0;
  for (const { z, D, cL, cR } of list) {
    a += z * D * (cR - cL);
    b += z * z * D * (cR - cL);
    gL += z * z * D * cL;
    gR += z * z * D * cR;
  }
  const VT = (GAS_CONSTANT * T) / FARADAY;
  if (b === 0) return 0;
  return (-VT * a * Math.log(gR / gL)) / b;
}

/**
 * Planck's liquid-junction potential, φ_R − φ_L (V): the steady state of constrained diffusion,
 * a junction zone held between the two solutions, electroneutral, carrying no current. Solved by
 * shooting on the Nernst–Planck equations across it (its length drops out). Exact for that
 * junction; equal to Henderson's for a single binary salt.
 * @param {Record<string, number>} left mol/m³ by ion
 * @param {Record<string, number>} right mol/m³ by ion
 * @param {Record<string, { z: number, D: number }>} ions charge and diffusivity (m²/s), e.g. IONS
 * @param {{ T?: number, steps?: number }} [opts] K (default 298.15); RK4 steps across (default 4000)
 */
export function planck(left, right, ions, { T = 298.15, steps = 4000 } = {}) {
  const list = prepare(left, right, ions, 'planck');
  neutral(list, 'cL', 'planck');
  neutral(list, 'cR', 'planck');
  const n = list.length, z = list.map((p) => p.z), D = list.map((p) => p.D), Dmax = Math.max(...D);
  const scale = Math.max(...list.map((p) => Math.max(p.cL, p.cR)));
  const cL = list.map((p) => p.cL / scale), cR = list.map((p) => p.cR / scale);
  // Across the zone (x from 0 to 1, D's scale dropping out with its length), each ion's flux J_i is
  // constant: c_i' = −J_i/D_i − z_i c_i ψ', with ψ = Fφ/RT, and neutrality fixes the field,
  // ψ' = −Σ(z_i J_i/D_i)/Σ(z_i² c_i). Integrated from the left; the fluxes are chosen so the
  // right comes out as given and no current flows.
  const deriv = (c, J, out) => {
    let num = 0, den = 0;
    for (let i = 0; i < n; i++) {
      num += (z[i] * J[i]) / D[i];
      den += z[i] * z[i] * Math.max(c[i], 0);
    }
    const dpsi = -num / den;
    for (let i = 0; i < n; i++) out[i] = -J[i] / D[i] - z[i] * c[i] * dpsi;
    return dpsi;
  };
  const shoot = (J) => {
    const c = cL.slice(), k1 = new Float64Array(n), k2 = new Float64Array(n), k3 = new Float64Array(n), k4 = new Float64Array(n), t = new Float64Array(n);
    const h = 1 / steps;
    let psi = 0;
    for (let s = 0; s < steps; s++) {
      const p1 = deriv(c, J, k1);
      for (let i = 0; i < n; i++) t[i] = c[i] + (h / 2) * k1[i];
      const p2 = deriv(t, J, k2);
      for (let i = 0; i < n; i++) t[i] = c[i] + (h / 2) * k2[i];
      const p3 = deriv(t, J, k3);
      for (let i = 0; i < n; i++) t[i] = c[i] + h * k3[i];
      const p4 = deriv(t, J, k4);
      for (let i = 0; i < n; i++) c[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
      psi += (h / 6) * (p1 + 2 * p2 + 2 * p3 + p4);
    }
    return { c, psi };
  };
  // Residuals: n − 1 ions reaching the right (the last follows by neutrality), and no current.
  const residual = (J) => {
    const { c, psi } = shoot(J), r = new Float64Array(n);
    for (let i = 0; i < n - 1; i++) r[i] = c[i] - cR[i];
    let I = 0;
    for (let i = 0; i < n; i++) I += z[i] * J[i];
    r[n - 1] = I / Dmax; // scaled like the concentrations
    return { r, psi };
  };
  let J = list.map((_, i) => -D[i] * (cR[i] - cL[i])); // plain diffusion to start
  let { r, psi } = residual(J);
  for (let it = 0; it < 50; it++) {
    const norm = Math.max(...r.map(Math.abs));
    if (norm < 1e-12) break;
    // Newton with a finite-difference Jacobian.
    const Jac = [];
    for (let k = 0; k < n; k++) {
      const dJ = 1e-7 * Math.max(Math.abs(J[k]), 1e-3 * D[k]), Jp = J.slice();
      Jp[k] += dJ;
      const rp = residual(Jp).r;
      Jac.push(Array.from(rp, (v, i) => (v - r[i]) / dJ));
    }
    const step = solveDense(Jac, r, n);
    let lambda = 1, next;
    for (;;) {
      const Jn = J.map((v, k) => v - lambda * step[k]);
      next = residual(Jn);
      if (Math.max(...next.r.map(Math.abs)) < norm || lambda < 1e-3) {
        J = Jn;
        break;
      }
      lambda /= 2;
    }
    ({ r, psi } = next);
  }
  if (!(Math.max(...r.map(Math.abs)) < 1e-8)) fail('planck: the shooting did not converge');
  return ((GAS_CONSTANT * T) / FARADAY) * psi;
}

// Solves Σ_k A[k][i] x_k = b_i (A given column by column), by Gaussian elimination with pivoting.
function solveDense(cols, b, n) {
  const M = Array.from({ length: n }, (_, i) => [...cols.map((col) => col[i]), b[i]]);
  for (let p = 0; p < n; p++) {
    let best = p;
    for (let i = p + 1; i < n; i++) if (Math.abs(M[i][p]) > Math.abs(M[best][p])) best = i;
    [M[p], M[best]] = [M[best], M[p]];
    for (let i = p + 1; i < n; i++) {
      const f = M[i][p] / M[p][p];
      for (let j = p; j <= n; j++) M[i][j] -= f * M[p][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}
