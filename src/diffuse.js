// A diffuse (Gouy–Chapman) layer treated analytically: the ions' Boltzmann layer beside a face,
// between the face's plane of closest approach and a strictly neutral bulk, held as one number
// per face side instead of on the grid.
//
// CONTRACT
//
// diffuseLayer(c, z, n, s, kappaUnits, out, lattice, draw) is a pure function of its arguments.
//   c (mol/m³): the bulk concentrations beside the layer, n of them (species absent there at 0);
//   z: their charges; s = F ψ/RT, the diffuse drop ψ = φ(plane) − φ(bulk) in thermal units;
//   kappaUnits = { A: √(ε RT/2)/F (mol/m² per √(mol/m³)), B: √(2 ε RT) (C/m² per √(mol/m³)) },
//   ε the layer's permittivity (F/m); lattice null (dilute ions) or { cMax, on }, the bulk's lattice
//   statistics (src/statistics.js), on[i] = 1 for the species on its sites (every charged one);
//   draw null, or draw[i] = 1 for the species whose transport across the layer is wanted (K below).
// Dilute, c_i(s) = c_i e^{−z_i s}; crowded (Bikerman: the same sites as the bulk's lattice), divided
// by D(s) = 1 + Σ_on θ_j (e^{−z_j s} − 1), θ_j = c_j/cMax, and the first integral is then
// G(s) = cMax ln D + s Σ c_j z_j, its neutral species displaced from the plane (an excess < 0).
// It writes into out (a scratch object it fills: { q, dqds, dqdc[n], gam[n], dgds[n], dgdc[n×n],
// K[n], dKds[n], dKdc[n×n], D, dDds, dDdc[n] }: D(s) and its derivatives, 1 and 0 for a dilute layer):
//   q    the layer's charge per area (C/m²), −B s √R(s), R = G/t²; dilute, R(t) = Σ_j c_j z_j² ρ(z_j t),
//        ρ(x) = (e^{−x} − 1 + x)/x², so R(t) t² = Σ c_j (e^{−z_j t} − 1 + z_j t), the first
//        integral of Poisson–Boltzmann (exact where the bulk is neutral, F Σ z c + ρ_fixed = 0,
//        with any fixed charge uniform through the layer as in the bulk; written with its linear
//        term so it stays ≥ 0 off neutrality, as Newton's iterates are);
//   gam  each species' surface excess ∫ (c_i(x) − c_i) dx (mol/m²), c_i A ∫₀ˢ E_i(t)/√R(t) dt,
//        E_i(t) = (c_i(t)/c_i − 1)/t: smooth through s = 0, where everything vanishes (crowded,
//        written so that no two large terms cancel, past |X| = 1 by G whole);
//   K    for each drawn species, the layer's resistance to it in excess of the bulk's it displaces,
//        as a length (m): K_i = ∫ (c_i/c_i(x) − 1) dx = A ∫₀ˢ H_i(t)/√R(t) dt, H_i = (c_i/c_i(t) − 1)/t.
//        A flux N_i toward the plane (mol/(m²·s)), through the layer's profile, lowers its level
//        there: a_plane = a_edge (1 − N_i K_i/(D_i c_i)), exact for a species too dilute to shape the
//        layer (the reactant a charged electrode repels, the case that matters: K_i grows as
//        c_i/c_i(plane)), first order in N otherwise; < 0 for one the layer gathers;
//   and their derivatives: dqds = ∂q/∂s, dqdc[j] = ∂q/∂c_j, dgds[i] = ∂Γ_i/∂s,
//   dgdc[i·n + j] = ∂Γ_i/∂c_j, dKds[i] = ∂K_i/∂s, dKdc[i·n + j] = ∂K_i/∂c_j (0 where not drawn).
// Invariants (each says what a reviewer should be able to check here):
//   1. F Σ z_i Γ_i = q at a neutral bulk (Gauss's law across the layer), to the quadrature's
//      accuracy (≤ 1e-13 relative).
//   2. q, Γ_i are odd in (s, z) together: flipping every charge and s gives −q and the same Γ.
//   3. The derivatives are those of the functions (the quadrature's error aside): ∂/∂s exact by
//      the fundamental theorem, ∂/∂c_j by the same rule as Γ (and K).
//   4. Against the Poisson–Boltzmann problem solved in x (scipy solve_bvp), and a z:z salt's
//      closed forms q = −√(8εRTc) sinh(zs/2), Γ_± = (2c/κ)(e^{∓zs/2} − 1): within 1e-9 and 1e-14
//      (test/diffuse.test.js, the tables from the answer sheets); crowded, against the same
//      problem with the lattice's concentrations, and Kilic–Bazant–Ajdari's q for a z:z salt;
//      K, against ∫ (e^{z φ̃(x)} − 1) dx on that problem's profile, and a z:z salt's closed form
//      (2/κ)(e^{±zs/2} − 1) (its own Γ with the charge flipped: c/c(x) = e^{+z φ̃}).
// A dilute layer overstates the charge and the excesses once the plane's concentration passes a
// few molar; a lattice bounds it at cMax. Arguments past |z s| = 600 are held there (e^600 is
// the end of double precision), so absurd Newton iterates stay finite.

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

// λ(X) = (log1p(X) − X)/X² (−1/2 at 0), by its series where the difference would cancel.
function lam(X) {
  if (Math.abs(X) < 0.1) {
    // Σ_{k≥0} (−1)^{k+1} X^k/(k + 2)
    let p = 1, sum = 0;
    for (let k = 0; k < 24; k++, p *= -X) sum -= p / (k + 2);
    return sum;
  }
  return (Math.log1p(X) - X) / (X * X);
}

// The layer's functions at t (a scratch `at` to fill): R(t) = G(t)/t², each species' E_i(t) =
// (c_i(t)/c_i − 1)/t, and their ∂/∂c_j, R_j and (crowded) Ec_ij = ∂E_i/∂c_j; and D(t), its ∂/∂t and
// ∂/∂c_j. Dilute: c_i(t) = c_i e^{−z_i t}; on a lattice of cMax sites, divided by D(t) = 1 + X(t),
// X = Σ_lattice θ_j (e^{−z_j t} − 1), θ_j = c_j/cMax, and G = cMax ln D + t Σ c_j z_j (Bikerman).
// Drawn (draw[i]), H_i(t) = (c_i/c_i(t) − 1)/t too, the resistance's integrand, and (crowded)
// Hc_ij = ∂H_i/∂c_j.
function functionsAt(c, z, n, t, lat, at, draw) {
  const { Rj, E, Ec, dDdc, H, Hc } = at;
  let R = 0;
  for (let j = 0; j < n; j++) {
    Rj[j] = z[j] !== 0 ? z[j] * z[j] * rho(z[j] * t) : 0;
    if (c[j] > 0) R += c[j] * Rj[j];
    E[j] = z[j] !== 0 ? -z[j] * eps1(z[j] * t) : 0;
    if (draw) H[j] = draw[j] && z[j] !== 0 ? z[j] * eps1(-z[j] * t) : 0; // (e^{z t} − 1)/t
  }
  at.D = 1;
  at.dDdt = 0;
  if (!lat) {
    at.R = R;
    return at;
  }
  // Y = X/t, smooth through 0; D' = −Σ θ z e^{−z t}
  const { cMax, on } = lat;
  let Y = 0, dD = 0;
  for (let j = 0; j < n; j++) {
    dDdc[j] = 0;
    if (!on[j] || z[j] === 0) continue;
    const x = Math.max(-CAP, Math.min(CAP, -z[j] * t));
    Y -= ((c[j] > 0 ? c[j] : 0) / cMax) * z[j] * eps1(z[j] * t);
    dD -= ((c[j] > 0 ? c[j] : 0) / cMax) * z[j] * Math.exp(x);
    dDdc[j] = Math.expm1(x) / cMax;
  }
  const X = Y * t, D = 1 + X;
  at.D = D;
  at.dDdt = dD;
  if (Math.abs(X) <= 1) {
    // the dilute layer's R and a correction, cMax (log1p X − X)/t² = cMax Y² λ(X); its ∂/∂c_j, Y z_j e(z_j t)/D
    at.R = R + cMax * Y * Y * lam(X);
    for (let j = 0; j < n; j++) if (on[j] && z[j] !== 0) Rj[j] += (Y * z[j] * eps1(z[j] * t)) / D;
  } else {
    // (far out, where the two would cancel: G whole, and ∂G/∂c_j = expm1(−z_j t)/D + z_j t)
    let Gl = cMax * Math.log1p(X);
    for (let j = 0; j < n; j++) if (on[j] && z[j] !== 0 && c[j] > 0) Gl += c[j] * z[j] * t;
    at.R = Gl / (t * t);
    for (let j = 0; j < n; j++) if (on[j] && z[j] !== 0) Rj[j] = (Math.expm1(Math.max(-CAP, Math.min(CAP, -z[j] * t))) / D + z[j] * t) / (t * t);
  }
  // E_i = (e^{−z_i t}/D − 1)/t = (−z_i e(z_i t) − Y)/D, a lattice's neutral species too; ∂E_i/∂c_j = e^{−z_i t} z_j e(z_j t)/(cMax D²)
  for (let i = 0; i < n; i++) {
    if (!on[i]) continue;
    E[i] = (E[i] - Y) / D;
    const w = Math.exp(Math.max(-CAP, Math.min(CAP, -z[i] * t))) / (cMax * D * D);
    for (let j = 0; j < n; j++) Ec[i * n + j] = on[j] && z[j] !== 0 ? w * z[j] * eps1(z[j] * t) : 0;
    // H_i = (e^{z_i t} D − 1)/t = z_i e(−z_i t) + e^{z_i t} Y; ∂H_i/∂c_j = −e^{z_i t} z_j e(z_j t)/cMax
    if (!draw?.[i]) continue;
    const up = Math.exp(Math.max(-CAP, Math.min(CAP, z[i] * t)));
    H[i] += up * Y;
    for (let j = 0; j < n; j++) Hc[i * n + j] = on[j] && z[j] !== 0 ? (-up * z[j] * eps1(z[j] * t)) / cMax : 0;
  }
  return at;
}

/** The diffuse layer's charge, excesses and their derivatives (see the contract above). */
export function diffuseLayer(c, z, n, s, { A, B }, out, lat = null, draw = null) {
  const { dqdc, gam, dgds, dgdc, K, dKds, dKdc } = out;
  let zmax = 0;
  for (let i = 0; i < n; i++) if (c[i] > 0) zmax = Math.max(zmax, Math.abs(z[i]));
  gam.fill(0);
  dgds.fill(0);
  dgdc.fill(0);
  dqdc.fill(0);
  K.fill(0);
  dKds.fill(0);
  dKdc.fill(0);
  out.D = 1;
  out.dDds = 0;
  out.dDdc.fill(0);
  if (zmax === 0) {
    out.q = 0;
    out.dqds = 0;
    return out;
  }
  const holds = (i) => z[i] !== 0 || (lat && lat.on[i]); // (a lattice's neutral species is displaced)
  const at = out.at ?? (out.at = { Rj: new Float64Array(n), E: new Float64Array(n), Ec: new Float64Array(n * n), dDdc: new Float64Array(n), H: new Float64Array(n), Hc: new Float64Array(n * n) });
  if (draw) at.Hc.fill(0);
  // At s itself: the charge and the ∂/∂s.
  functionsAt(c, z, n, s, lat, at, draw);
  const sqR = Math.sqrt(at.R);
  out.q = -B * s * sqR;
  out.D = at.D;
  out.dDds = at.dDdt;
  out.dDdc.set(at.dDdc);
  let dq = 0;
  for (let i = 0; i < n; i++) {
    if (!(c[i] > 0) || !holds(i)) continue;
    dgds[i] = (c[i] * A * at.E[i]) / sqR;
    dq += c[i] * z[i] * at.E[i];
  }
  if (draw) for (let i = 0; i < n; i++) if (draw[i]) dKds[i] = (A * at.H[i]) / sqR;
  out.dqds = (B * dq) / (2 * sqR);
  for (let j = 0; j < n; j++) dqdc[j] = (-B * s * at.Rj[j]) / (2 * sqR);
  if (s === 0) return out;
  // ∫₀ˢ by panels a unit of z_max·t wide (each spans at most an e-fold of the integrand), 8 points each.
  const panels = Math.max(1, Math.ceil(Math.min(Math.abs(s), CAP) * zmax)), h = s / panels;
  const I = out.I ?? (out.I = new Float64Array(n));
  I.fill(0);
  for (let p = 0; p < panels; p++) {
    for (let k = 0; k < 8; k++) {
      const t = (p + GL_X[k]) * h, w = GL_W[k] * h;
      functionsAt(c, z, n, t, lat, at, draw);
      const r12 = 1 / Math.sqrt(at.R), r32 = r12 / at.R;
      if (draw) {
        // K_i = A ∫ H_i R^{-1/2}; its ∂/∂c_j, A ∫ (Hc_ij R^{-1/2} − ½ H_i R_j R^{-3/2})
        for (let i = 0; i < n; i++) {
          if (!draw[i]) continue;
          K[i] += w * at.H[i] * r12;
          const f = -0.5 * w * at.H[i] * r32;
          for (let j = 0; j < n; j++) dKdc[i * n + j] += f * at.Rj[j] + w * at.Hc[i * n + j] * r12;
        }
      }
      for (let i = 0; i < n; i++) {
        if (!holds(i)) continue;
        I[i] += w * at.E[i] * r12;
        // ∂/∂c_j of c_i ∫ E_i R^{-1/2}: c_i ∫ (∂E_i/∂c_j R^{-1/2} − ½ E_i R_j R^{-3/2})
        if (!(c[i] > 0)) continue;
        const f = -0.5 * w * at.E[i] * r32 * c[i];
        for (let j = 0; j < n; j++) dgdc[i * n + j] += f * at.Rj[j] + (lat && lat.on[i] ? w * c[i] * at.Ec[i * n + j] * r12 : 0);
      }
    }
  }
  for (let i = 0; i < n; i++) {
    if (!holds(i)) continue;
    gam[i] = (c[i] > 0 ? c[i] : 0) * A * I[i];
    for (let j = 0; j < n; j++) dgdc[i * n + j] *= A;
    dgdc[i * n + i] += A * I[i];
  }
  if (draw) {
    for (let i = 0; i < n; i++) K[i] *= A;
    for (let k = 0; k < n * n; k++) dKdc[k] *= A;
  }
  return out;
}

/** A scratch object for diffuseLayer over n species. */
export function diffuseScratch(n) {
  return { q: 0, dqds: 0, dqdc: new Float64Array(n), gam: new Float64Array(n), dgds: new Float64Array(n), dgdc: new Float64Array(n * n), K: new Float64Array(n), dKds: new Float64Array(n), dKdc: new Float64Array(n * n), D: 1, dDds: 0, dDdc: new Float64Array(n) };
}

/** The constants A = √(εRT/2)/F and B = √(2εRT) for a layer of permittivity eps (F/m) at RT (J/mol). */
export function diffuseUnits(eps, RT, F) {
  return { A: Math.sqrt((eps * RT) / 2) / F, B: Math.sqrt(2 * eps * RT) };
}
