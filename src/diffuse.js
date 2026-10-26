// A diffuse (Gouy–Chapman) layer treated analytically: the ions' Boltzmann layer beside a face,
// between the face's plane of closest approach and a strictly neutral bulk, held as one number
// per face side instead of on the grid.
//
// CONTRACT
//
// diffuseLayer(c, z, n, s, kappaUnits, out) is a pure function of its arguments.
//   c (mol/m³): the bulk concentrations beside the layer, n of them (species absent there at 0);
//   z: their charges; s = F ψ/RT, the diffuse drop ψ = φ(plane) − φ(bulk) in thermal units;
//   kappaUnits = { A: √(ε RT/2)/F (mol/m² per √(mol/m³)), B: √(2 ε RT) (C/m² per √(mol/m³)) },
//   ε the layer's permittivity (F/m).
// It writes into out (a scratch object it fills: { q, dqds, dqdc[n], gam[n], dgds[n], dgdc[n×n] }):
//   q    the layer's charge per area (C/m²), −B s √R(s), R(t) = Σ_j c_j z_j² ρ(z_j t),
//        ρ(x) = (e^{−x} − 1 + x)/x², so R(t) t² = Σ c_j (e^{−z_j t} − 1 + z_j t), the first
//        integral of Poisson–Boltzmann (exact where the bulk is neutral, F Σ z c + ρ_fixed = 0,
//        with any fixed charge uniform through the layer as in the bulk; written with its linear
//        term so it stays ≥ 0 off neutrality, as Newton's iterates are);
//   gam  each species' surface excess ∫ (c_i(x) − c_i) dx (mol/m²), c_i A ∫₀ˢ E_i(t)/√R(t) dt,
//        E_i(t) = (e^{−z_i t} − 1)/t: smooth through s = 0, where everything vanishes;
//   and their derivatives: dqds = ∂q/∂s, dqdc[j] = ∂q/∂c_j, dgds[i] = ∂Γ_i/∂s,
//   dgdc[i·n + j] = ∂Γ_i/∂c_j.
// Invariants (each says what a reviewer should be able to check here):
//   1. F Σ z_i Γ_i = q at a neutral bulk (Gauss's law across the layer), to the quadrature's
//      accuracy (≤ 1e-13 relative).
//   2. q, Γ_i are odd in (s, z) together: flipping every charge and s gives −q and the same Γ.
//   3. The derivatives are those of the functions (the quadrature's error aside): ∂/∂s exact by
//      the fundamental theorem, ∂/∂c_j by the same rule as Γ.
//   4. Against the Poisson–Boltzmann problem solved in x (scipy solve_bvp), and a z:z salt's
//      closed forms q = −√(8εRTc) sinh(zs/2), Γ_± = (2c/κ)(e^{∓zs/2} − 1): within 1e-9 and 1e-14
//      (test/diffuse.test.js, the table from the answer sheet).
// Ideal ions only: the excess and the charge are those of a dilute Boltzmann layer, which
// overstates both once the plane's concentration passes a few molar (no crowding); a Stern layer
// in series is what keeps s modest in practice. Arguments past |z s| = 600 are held there (e^600
// is the end of double precision), so absurd Newton iterates stay finite.

// Gauss–Legendre, 8 points on [0, 1].
const GL_X = [0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751, 0.5917173212478249, 0.7627662049581645, 0.8983332387068134, 0.9801449282487681];
const GL_W = [0.05061426814518813, 0.11119051722668724, 0.15685332293894363, 0.18134189168918100, 0.18134189168918100, 0.15685332293894363, 0.11119051722668724, 0.05061426814518813];
const CAP = 600;

// ρ(x) = (e^{−x} − 1 + x)/x² (1/2 at 0), by its series where the difference would cancel.
function rho(x) {
  if (x > CAP) x = CAP;
  else if (x < -CAP) x = -CAP;
  if (Math.abs(x) < 0.5) {
    // Σ_{k≥0} (−x)^k/(k+2)!
    let term = 0.5, sum = 0.5;
    for (let k = 1; k < 18; k++) {
      term *= -x / (k + 2);
      sum += term;
    }
    return sum;
  }
  return (Math.expm1(-x) + x) / (x * x);
}

// e(x) = (1 − e^{−x})/x (1 at 0), so E_i(t) = −z_i e(z_i t).
function eps1(x) {
  if (x > CAP) x = CAP;
  else if (x < -CAP) x = -CAP;
  return x === 0 ? 1 : -Math.expm1(-x) / x;
}

/** The diffuse layer's charge, excesses and their derivatives (see the contract above). */
export function diffuseLayer(c, z, n, s, { A, B }, out) {
  const { dqdc, gam, dgds, dgdc } = out;
  let zmax = 0;
  for (let i = 0; i < n; i++) if (c[i] > 0) zmax = Math.max(zmax, Math.abs(z[i]));
  gam.fill(0);
  dgds.fill(0);
  dgdc.fill(0);
  dqdc.fill(0);
  if (zmax === 0) {
    out.q = 0;
    out.dqds = 0;
    return out;
  }
  // R and E at s itself: the charge and the ∂/∂s.
  let Rs = 0;
  for (let j = 0; j < n; j++) if (c[j] > 0) Rs += c[j] * z[j] * z[j] * rho(z[j] * s);
  const sqR = Math.sqrt(Rs);
  out.q = -B * s * sqR;
  let dq = 0;
  for (let i = 0; i < n; i++) {
    if (!(c[i] > 0) || z[i] === 0) continue;
    const E = -z[i] * eps1(z[i] * s);
    dgds[i] = (c[i] * A * E) / sqR;
    dq += c[i] * z[i] * E;
  }
  out.dqds = (B * dq) / (2 * sqR);
  for (let j = 0; j < n; j++) if (c[j] > 0 || z[j] !== 0) dqdc[j] = (-B * s * z[j] * z[j] * rho(z[j] * s)) / (2 * sqR);
  if (s === 0) return out;
  // ∫₀ˢ by panels a unit of z_max·t wide (each spans at most an e-fold of the integrand), 8 points each.
  const panels = Math.max(1, Math.ceil(Math.min(Math.abs(s), CAP) * zmax)), h = s / panels;
  const Rj = out.Rj ?? (out.Rj = new Float64Array(n)), I = out.I ?? (out.I = new Float64Array(n));
  I.fill(0);
  for (let p = 0; p < panels; p++) {
    for (let k = 0; k < 8; k++) {
      const t = (p + GL_X[k]) * h, w = GL_W[k] * h;
      let R = 0;
      for (let j = 0; j < n; j++) {
        Rj[j] = c[j] > 0 || z[j] !== 0 ? z[j] * z[j] * rho(z[j] * t) : 0;
        if (c[j] > 0) R += c[j] * Rj[j];
      }
      const r12 = 1 / Math.sqrt(R), r32 = r12 / R;
      for (let i = 0; i < n; i++) {
        if (z[i] === 0) continue;
        const E = -z[i] * eps1(z[i] * t);
        I[i] += w * E * r12;
        // ∂/∂c_j of c_i ∫ E_i R^{-1/2}: −½ c_i ∫ E_i R_j R^{-3/2}
        const f = -0.5 * w * E * r32 * (c[i] > 0 ? c[i] : 0);
        if (f !== 0) for (let j = 0; j < n; j++) dgdc[i * n + j] += f * Rj[j];
      }
    }
  }
  for (let i = 0; i < n; i++) {
    if (z[i] === 0) continue;
    gam[i] = (c[i] > 0 ? c[i] : 0) * A * I[i];
    for (let j = 0; j < n; j++) dgdc[i * n + j] *= A;
    dgdc[i * n + i] += A * I[i];
  }
  return out;
}

/** A scratch object for diffuseLayer over n species. */
export function diffuseScratch(n) {
  return { q: 0, dqds: 0, dqdc: new Float64Array(n), gam: new Float64Array(n), dgds: new Float64Array(n), dgdc: new Float64Array(n * n) };
}

/** The constants A = √(εRT/2)/F and B = √(2εRT) for a layer of permittivity eps (F/m) at RT (J/mol). */
export function diffuseUnits(eps, RT, F) {
  return { A: Math.sqrt((eps * RT) / 2) / F, B: Math.sqrt(2 * eps * RT) };
}
