// Liquid-junction potentials from formulas, for driftlet/kit: the diffusion potential φ_R − φ_L
// between two neutral solutions, as the classical theories give it. A device simulates the actual
// junction (steady between two baths, it is Planck's; stepping from a sharp boundary, a free
// diffusion junction); these are for comparison, and for quick estimates.

import { DeviceError, checkOptions } from './errors.js';
import { FARADAY, GAS_CONSTANT } from './constants.js';
import { solveDense } from './dense.js';

const fail = (message) => {
  throw new DeviceError(message);
};

// The ions both solutions name, with their charge and diffusivity, and both concentrations. A
// neutral solute carries no current and feels no field, so it has no part in the potential: it's
// dropped. (Kept, listed last in planck's, its concentration was the one neutrality fixes, and the
// shooting "did not converge"; Astra's fifth pass.)
function prepare(left, right, ions, what, T) {
  if (!(T > 0 && T < Infinity)) fail(`${what}: T is the temperature in K, a positive number; got ${T}`);
  const names = [...new Set([...Object.keys(left), ...Object.keys(right)])];
  return names.flatMap((name) => {
    const ion = ions[name];
    if (!ion || !Number.isFinite(ion.z) || !(ion.D > 0)) fail(`${what}: no { z, D } for '${name}' in the ions given`);
    const cL = left[name] ?? 0, cR = right[name] ?? 0;
    if (!(cL >= 0 && cR >= 0)) fail(`${what}: concentrations must be ≥ 0 (mol/m³), got ${name} ${cL} and ${cR}`);
    return ion.z === 0 ? [] : [{ name, z: ion.z, D: ion.D, cL, cR }];
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
function holds(list, what) {
  for (const side of ['cL', 'cR']) if (!list.some((p) => p[side] > 0)) fail(`${what}: the ${side === 'cL' ? 'left' : 'right'} solution holds no ions`);
}

/**
 * Henderson's liquid-junction potential, φ_R − φ_L (V): the junction taken to be the two solutions
 * mixed in one proportion at every point (every ion's mixing fraction the same curve). Exact for a single binary salt, a close estimate otherwise (it's what
 * JPCalc and most electrophysiology corrections use).
 * @param {Record<string, number>} left mol/m³ by ion
 * @param {Record<string, number>} right mol/m³ by ion
 * @param {Record<string, { z: number, D: number }>} ions charge and diffusivity (m²/s), e.g. IONS
 * @param {{ T?: number }} [opts] K (default 298.15)
 */
export function henderson(left, right, ions, { T = 298.15 } = {}) {
  checkOptions('henderson', arguments[3], ['T']);
  const list = prepare(left, right, ions, 'henderson', T);
  neutral(list, 'cL', 'henderson');
  neutral(list, 'cR', 'henderson');
  holds(list, 'henderson');
  let a = 0, b = 0, gL = 0, gR = 0;
  for (const { z, D, cL, cR } of list) {
    a += z * D * (cR - cL);
    b += z * z * D * (cR - cL);
    gL += z * z * D * cL;
    gR += z * z * D * cR;
  }
  // ln(g_R/g_L)/(g_R − g_L) = ln(1 + r)/(r g_L), r = b/g_L: 1/g_L at equal conductivities, where the
  // plain form is 0/0; near them, log1p; and beyond, the ratio's own log (1 + r rounded, 100 mM
  // against 1e-15 read −∞, and at 1e-12 64 µV off; Astra's fifth pass).
  const VT = (GAS_CONSTANT * T) / FARADAY, r = b / gL;
  return (-VT * a * (Math.abs(r) < 1e-8 ? 1 - r / 2 : Math.abs(r) < 0.5 ? Math.log1p(r) / r : Math.log(gR / gL) / r)) / gL;
}

/**
 * Planck's liquid-junction potential, φ_R − φ_L (V): the steady state of constrained diffusion,
 * a junction zone held between the two solutions, electroneutral, carrying no current. Solved by
 * shooting on the Nernst–Planck equations across it (its length drops out). Exact for that
 * junction; equal to Henderson's for a single binary salt.
 * @param {Record<string, number>} left mol/m³ by ion
 * @param {Record<string, number>} right mol/m³ by ion
 * @param {Record<string, { z: number, D: number }>} ions charge and diffusivity (m²/s), e.g. IONS
 * @param {{ T?: number, steps?: number }} [opts] K (default 298.15); RK4 steps across (default 4000,
 *   graded toward the dilute side; with half as many the answer must agree to 1 µV, or it throws)
 */
export function planck(left, right, ions, { T = 298.15, steps = 4000 } = {}) {
  checkOptions('planck', arguments[3], ['T', 'steps']);
  if (!(Number.isInteger(steps) && steps >= 1)) fail(`planck: steps is the number of integration steps across the junction, a whole number ≥ 1; got ${steps}`);
  const list = prepare(left, right, ions, 'planck', T);
  neutral(list, 'cL', 'planck');
  neutral(list, 'cR', 'planck');
  holds(list, 'planck');
  const n = list.length, z = list.map((p) => p.z), D = list.map((p) => p.D), Dmax = Math.max(...D);
  const scale = Math.max(...list.map((p) => Math.max(p.cL, p.cR)));
  // Shot from the dilute side, the sign turned back at the end: aimed at a dilute side, the
  // shooting can't hit concentrations near its tolerance (100 mM against 1 nM didn't converge),
  // and from it, it holds to 1e12:1 (Astra's fifth pass).
  const g = (side) => list.reduce((t, p) => t + p.z * p.z * p[side], 0), sign = g('cL') > g('cR') ? -1 : 1;
  const [from, to] = sign > 0 ? ['cL', 'cR'] : ['cR', 'cL'];
  const cL = list.map((p) => p[from] / scale), cR = list.map((p) => p[to] / scale);
  // Across the zone (x from 0 to 1, D's scale dropping out with its length), each ion's flux J_i is
  // constant: c_i' = −J_i/D_i − z_i c_i ψ', with ψ = Fφ/RT, and neutrality fixes the field,
  // ψ' = −Σ(z_i J_i/D_i)/Σ(z_i² c_i). Integrated from the dilute side (cL here); the fluxes are
  // chosen so the other side (cR) comes out as given and no current flows.
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
  // The field goes as 1/Σz²c, so toward a dilute side it steepens without bound: steps even in x
  // missed it, 100 mM against 0.1 mM −274 mV where −73.6 is right, and against 10 µM already
  // 0.13 mV off (Astra's fifth pass). So the steps are even in ξ, with Σz²c taken to go
  // geometrically across, g(x) = g_L (g_R/g_L)^ξ: x = expm1(Lξ)/expm1(L), L = ln(g_R/g_L), each
  // step in x as long as g is large there, which follows a binary salt's profile (c linear in x)
  // to within the steps' resolution of its logarithm.
  const gOf = (cs) => cs.reduce((t, c, i) => t + z[i] * z[i] * c, 0), L = Math.log(gOf(cR) / gOf(cL));
  const dxdxi = (xi) => (L === 0 ? 1 : (L * Math.exp(L * xi)) / Math.expm1(L));
  const shoot = (J, N = steps) => {
    const c = cL.slice(), k1 = new Float64Array(n), k2 = new Float64Array(n), k3 = new Float64Array(n), k4 = new Float64Array(n), t = new Float64Array(n);
    const h = 1 / N;
    // (d/dξ = dx/dξ · d/dx at each stage's own ξ)
    const at = (xi, cs, out) => {
      const w = dxdxi(xi), p = deriv(cs, J, out);
      for (let i = 0; i < n; i++) out[i] *= w;
      return p * w;
    };
    let psi = 0;
    for (let s = 0; s < N; s++) {
      const xi = s * h;
      const p1 = at(xi, c, k1);
      for (let i = 0; i < n; i++) t[i] = c[i] + (h / 2) * k1[i];
      const p2 = at(xi + h / 2, t, k2);
      for (let i = 0; i < n; i++) t[i] = c[i] + (h / 2) * k2[i];
      const p3 = at(xi + h / 2, t, k3);
      for (let i = 0; i < n; i++) t[i] = c[i] + h * k3[i];
      const p4 = at(xi + h, t, k4);
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
    const step = solveDense(Jac, r, { byColumns: true });
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
  // Resolved: with half the steps, the same potential to 1 µV (RK4's error falls 16× a halving, so
  // the difference is about 15× this answer's own error). Unresolved, it says so.
  const VT = (GAS_CONSTANT * T) / FARADAY;
  if (steps >= 2) {
    const coarse = shoot(J, Math.ceil(steps / 2)).psi;
    if (!(VT * Math.abs(coarse - psi) <= 1e-6)) fail(`planck: the integration across the junction isn't resolved: with ${steps} steps it gives ${(1e3 * sign * VT * psi).toPrecision(6)} mV, with ${Math.ceil(steps / 2)} ${(1e3 * sign * VT * coarse).toPrecision(6)} mV; pass more steps`);
  }
  return sign * VT * psi;
}
