// The one owner of strictly neutral charge: what a strictly neutral node's Gauss row says, and
// how a state is brought onto it.
//
// A node of a strictly neutral (ε = 0) region is an index-1 DAE: its balances, with an algebraic
//   G_g(u) = −v_g (ρ_fixed + F Σ z_i c_i) − v_g Σ_k a_k σ_k + (what box g holds of a face's or a
//            contact's displacement) = 0,
// σ_k = C_k (V_k − zeroCharge_k − φ) a port's charge per area of electrode where it holds one
// (Solver._portHolds). That row is assembled in src/regions.js and src/boundaries.js and nowhere
// else: here it's read from an assembly at the state (gaussRows), so nothing that brings a state
// onto it can disagree with what the step then solves. (Each producer summing its own subset of
// it was the seam of Neon's R1 1.1, R2 3.1, R3 1.5 and R3 1.7.)
//
// Contract:
//   1. strictlyNeutralNodes(s): every node of a strictly neutral non-metal region where φ is
//      defined, but an end whose contact pins φ (its row is the pin); each brought onto G = 0 by
//      default by 'level' (φ̂ moved, the levels η kept: an open solution's levels are what it has,
//      φ̂ what neutrality makes of them). Or by 'charge' (φ̂ kept, each charged species' η_i moved
//      by −z_i δ: the concentrations in the charge direction only, no face's D or port's σ moving),
//      which src/carry.js uses where the ions whose amounts are kept carry most of a node's
//      screening (Σ z²c), and so hold its net charge; a node's `only` (a mask by species) then
//      names the ions that move, the open ones keeping their levels.
//   2. gaussRows(s, nodes, { syncFaces }): one assemble(Infinity) with s.steady off (so no law
//      takes a row's place); with syncFaces, each capacitive face's D first put on its own row (its
//      slot u, which the Gauss rows read: the face's row is linear in it with slope 1), but beside
//      a diffuse layer, whose row ties D to the layer's charge. (The carry syncs, its faces carried
//      on their rows; a cold start doesn't: its D is set to balance its edge boxes, the face's own
//      row left to the first step, and synced there the boxes were 0.59 of their charge off.) For each node: G (row 0's residual; at a
//      bulk contact's end, whose row 0 is the contact's, the ions alone: the bulk law's G), ions
//      v(ρ + F Σ z c), rest = G + ions (everything but the ions), slope ∂rest/∂φ̂ at the node (the
//      ports' v Σ aC V_T, a capacitive face's or contact's A C V_T; the other side held), and
//      scale (the gross charges). Writes what assemble writes (res, sys, the readouts).
//   3. projectNeutral(s, nodes, { tol }): each node, by its mode, brought onto G = 0 to 1e-15 of
//      its scale, ideal nodes in closed form, others through their statistics, rows read once a
//      round, rounds until every node is within tol of its scale (or 40 rounds). A node already
//      within 1e-12 of its scale is left bit for bit. 'jump' (φ̂ moved with each η_i by z_i δ: the
//      concentrations kept) is for a set() of a gate's hardware: see src/carry.js. Returns the
//      worst |G|/scale left and marks each node it moved in `moved`.
// Tests: test/neutralOwner.test.js; bench/combos' neutralG property (every family).

import { FARADAY } from './constants.js';

export function strictlyNeutralNodes(s) {
  const out = [];
  const end = (g) => (g === 0 ? s.model.contacts.left : g === s.nNodes - 1 ? s.model.contacts.right : null);
  for (let g = 0; g < s.nNodes; g++) {
    const mat = s.model.materials[s.nodeMaterial[g]];
    if (mat.conductor || mat.epsr !== 0 || s.phiUndefined[g] || s.nodeConductor[g] >= 0) continue;
    const ct = end(g), law = ct?.phi.type;
    if (law === 'pinned') continue;
    out.push({ g, mode: 'level', bulk: law === 'bulk' });
  }
  return out;
}

// A capacitive face's D, or a capacitive contact's, beside node g: the slope of what box g holds of
// it by φ̂ there (A C V_T, the other side held).
function boundarySlope(s, g) {
  const { model, VT } = s, grid = model.grid;
  let k = 0;
  const r = grid.nodeRegion[g];
  if (g === grid.regionEnd[r] && r < model.interfaces.length) {
    const itf = model.interfaces[r];
    if (itf.phi.type === 'capacitive' && !s.diffList.some((d) => d.f === r && d.side === 'left')) k += grid.area[g] * itf.phi.C * VT;
  }
  if (g === grid.regionStart[r] && r > 0) {
    const itf = model.interfaces[r - 1];
    if (itf.phi.type === 'capacitive' && !s.diffList.some((d) => d.f === r - 1 && d.side === 'right')) k += grid.area[grid.regionEnd[r - 1]] * itf.phi.C * VT;
  }
  const side = g === 0 ? 'left' : g === s.nNodes - 1 ? 'right' : null;
  if (side && model.contacts[side].phi.type === 'capacitive') k += grid.area[g] * model.contacts[side].phi.C * VT;
  return k;
}

export function gaussRows(s, nodes, { syncFaces = false } = {}) {
  const { n, M, z, c, present, VT } = s, R = s.rix, grid = s.model.grid, F = FARADAY;
  const steady = s.steady;
  s.steady = false;
  try {
    s.assemble(Infinity);
    // A capacitive face's D on its own row (linear in its slot, slope 1), and the Gauss rows beside
    // it moved by what that moved. (Left as it was where it's within round-off of it: a state that
    // holds it, carried, stays bit for bit.)
    if (syncFaces) s.model.interfaces.forEach((itf, f) => {
      // (not beside a diffuse layer, whose own row ties D to the layer's charge: there the face's
      // row is φ̂'s to meet; synced, a layer whose ε changed lost its charge to the box beside it)
      if (itf.phi.type !== 'capacitive' || s.diffList.some((d) => d.f === f)) return;
      const bf = s.blockOfFace[f], l = s.loc[bf * M];
      if (l < 0) return;
      const r = s.res[R[bf * M]], D = s.u[bf * M];
      if (!(Math.abs(r) > 1e-13 * (Math.abs(D) + itf.phi.C * VT))) return;
      s.u[bf * M] = D - r;
      const gL = grid.regionEnd[f], gR = gL + 1, Af = grid.area[gL];
      if (s.nodeConductor[gL] < 0) s.res[R[s.blockOfNode[gL] * M]] -= Af * r;
      if (s.nodeConductor[gR] < 0) s.res[R[s.blockOfNode[gR] * M]] += Af * r;
    });
  } finally {
    s.steady = steady;
  }
  return nodes.map(({ g, bulk }) => {
    const v = grid.vol[g], b = s.blockOfNode[g];
    let q = s.rhoFixed[g], gross = Math.abs(q);
    for (let i = 0; i < n; i++) if (present[g * n + i]) (q += F * z[i] * c[g * n + i]), (gross += F * Math.abs(z[i] * c[g * n + i]));
    const ions = v * q, G = bulk ? -ions : s.res[R[b * M]];
    const rest = G + ions;
    const slope = bulk ? 0 : v * (s.qgSlope[g] ?? 0) * VT + boundarySlope(s, g);
    return { G, ions, rest, slope, scale: v * (gross + Math.abs(s.qg[g] ?? 0)) + Math.abs(rest + v * (s.qg[g] ?? 0)) };
  });
}

// Node g's ions v(ρ + F Σ z c) and their slope by d, with each moved species' ζ = η − zφ̂ moved
// by −z d (ideal: c e^{−z d}; otherwise through its statistics, ∂c_i/∂d = −Σ_j K_ij z_j over the
// moved j, left at d). `only`, where given, is which species move (a mixed node's closed ions);
// otherwise every charged one does.
function ionsAt(s, g, b, d, c0, only) {
  const { n, M, z, c, present, K, u } = s, F = FARADAY, v = s.model.grid.vol[g];
  const moves = (i) => z[i] !== 0 && (!only || only[i] === 1);
  let q = s.rhoFixed[g], dq = 0;
  if (s.nodeIdeal[g]) {
    for (let i = 0; i < n; i++) {
      if (!present[g * n + i] || z[i] === 0) continue;
      const ci = moves(i) ? c0[i] * Math.exp(-z[i] * d) : c0[i];
      q += F * z[i] * ci;
      if (moves(i)) dq -= F * z[i] * z[i] * ci;
    }
  } else {
    if (only) {
      for (let i = 0; i < n; i++) if (present[g * n + i] && moves(i)) u[b * M + 1 + i] -= z[i] * d;
      s._nodeStatistics(g, b, u[b * M]);
      for (let i = 0; i < n; i++) if (present[g * n + i] && moves(i)) u[b * M + 1 + i] += z[i] * d;
    } else s._nodeStatistics(g, b, u[b * M] + d);
    const Kg = g * n * n;
    for (let i = 0; i < n; i++) {
      if (!present[g * n + i]) continue;
      q += F * z[i] * c[g * n + i];
      let Kz = 0;
      for (let j = 0; j < n; j++) if (moves(j)) Kz += K[Kg + i * n + j] * z[j];
      dq -= F * z[i] * Kz;
    }
  }
  return [v * q, v * dq];
}

export function projectNeutral(s, nodes, { tol = 1e-13, moved = null, syncFaces = false } = {}) {
  if (nodes.length === 0) return 0;
  const { n, M, z, present, c } = s;
  let worst = 0;
  for (let round = 0; round < 40; round++) {
    const rows = gaussRows(s, nodes, { syncFaces });
    worst = 0;
    for (let a = 0; a < nodes.length; a++) {
      const { g, mode } = nodes[a], row = rows[a], b = s.blockOfNode[g];
      const off = Math.abs(row.G) / Math.max(row.scale, 1e-300);
      worst = Math.max(worst, off);
      // (a jump to round-off: what it moves is a gate's charge, which can be 1e-4 of the ions'
      // that the scale holds; left at 1e-12 of that, 1e-9 of a gate's charge read as a current)
      if (!(off > (mode === 'jump' ? 1e-15 : 1e-12))) continue;
      // The unknown d: 'level' moves φ̂ by d, its levels kept (each ζ = η − zφ̂ by −z d, and the
      // ports' and faces' terms by slope·d); 'charge' moves each charged η by −z d at the φ̂ it
      // has (the ions alike, the rest not at all); 'jump' moves φ̂ by d and each η by z d, the
      // concentrations kept. G(d) = rest + slope·d − ions(d) (no slope for 'charge'), rising in d.
      const c0 = Float64Array.from({ length: n }, (_, i) => c[g * n + i]);
      const kRest = mode === 'charge' ? 0 : row.slope;
      let d = 0;
      if (mode === 'jump') {
        if (!(row.slope > 0)) continue;
        d = -row.G / row.slope;
      } else {
        for (let it = 0; it < 60; it++) {
          const [ions, dIons] = ionsAt(s, g, b, d, c0, mode === 'charge' ? nodes[a].only : undefined);
          const f = row.rest + kRest * d - ions, df = kRest - dIons;
          if (Math.abs(f) <= 1e-15 * row.scale || !(df > 0)) break;
          d -= f / df;
          if (!Number.isFinite(d)) break;
        }
        if (!Number.isFinite(d)) continue;
      }
      if (d === 0) continue;
      if (mode === 'charge') {
        const only = nodes[a].only;
        for (let i = 0; i < n; i++) if (present[g * n + i] && z[i] !== 0 && (!only || only[i] === 1)) s.u[b * M + 1 + i] -= z[i] * d;
      } else if (mode === 'level') s.u[b * M] += d;
      else {
        s.u[b * M] += d;
        for (let i = 0; i < n; i++) if (present[g * n + i] && z[i] !== 0) s.u[b * M + 1 + i] += z[i] * d;
      }
      if (moved) moved[g] = 1;
    }
    s.computeConcentrations();
    if (worst <= tol) break;
  }
  return worst;
}
