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
//      φ̂ what neutrality makes of them), or by 'jump' (φ̂ moved with each η_i by z_i δ: the
//      concentrations kept), where a capacitor's charge can't jump and its box's φ̂ does (a set()
//      of its hardware, a cold start's uncharged capacitors). projectNeutral also has 'charge' (φ̂
//      kept, each charged η_i moved by −z_i δ, or only the species a node's `only` names).
//   2. gaussRows(s, nodes, { syncFaces }): one assemble(Infinity) with s.steady off (so no law
//      takes a row's place); with syncFaces, each capacitive face's D first put on its own row (its
//      slot u, which the Gauss rows read: the face's row is linear in it with slope 1), but beside
//      a diffuse layer, whose row ties D to the layer's charge. (The carry syncs, its faces carried
//      on their rows; a cold start doesn't: its D is set to balance its edge boxes, the face's own
//      row left to the first step, and synced there the boxes were 0.59 of their charge off.) For
//      each node: G (row 0's residual; at a bulk contact's end, whose row 0 is the contact's, the
//      ions alone: the bulk law's G), ions v(ρ + F Σ z c), rest = G + ions (everything but the
//      ions), slope ∂rest/∂φ̂ at the node (the ports' v Σ aC V_T, a capacitive contact's A C V_T, a
//      capacitive face's where synced and without a diffuse layer; the other side held), and scale
//      (the gross charges). Writes what assemble writes (res, sys, the readouts).
//   3. projectNeutral(s, nodes, { tol }): each node alone, by its mode, onto G = 0 to 1e-15 of its
//      scale, rows read once a round, rounds until every node is within tol of its scale (or 40).
//      A node already within 1e-12 of its scale is left bit for bit. Returns the worst |G|/scale
//      left and marks each node it moved in `moved`. (No producer uses it: they hold more than
//      neutrality at once, which by turns didn't converge; see 4.)
//   4. solveNeutral(s, nodes, { closed, floating, syncFaces }): Newton on every node's d (by its
//      mode), each closed stretch's ℓ (its species' η shifted over it, but where a node jumps;
//      with closedAtNodes, only over the listed nodes, keeping the cold start's c0 elsewhere)
//      and each floating capacitor's V together, onto G = 0, each closed amount its `want` and
//      each floating charge its `want`; the nodes' block tridiagonal (a synced capacitive face
//      couples its two boxes, −A C V_T), the border by its Schur complement. A 'jump' node nothing
//      holds by its φ̂ (no capacitor in its box) moves at its levels instead; a node nothing holds
//      at all is left. A state that holds already is left bit for bit. Returns { worst }, the worst
//      of the three residuals, each relative to its scale; marks moved nodes. The carry (each
//      regrid, hardware change) and the transient's cold start (src/start.js, neutralStart) use it.
// Tests: test/neutralOwner.test.js, test/neon4.test.js; bench/combos' neutralG property (every family).

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
function boundarySlope(s, g, faces = true) {
  const { model, VT } = s, grid = model.grid;
  let k = 0;
  const r = grid.nodeRegion[g];
  // (a diffuse layer's side is 0, left of its face, or 1: compared with 'left' and 'right', every
  // box beside a layer took this slope, and a jump there moved φ̂ 5 V; Neon's fourth review, R4 1.4)
  if (faces && g === grid.regionEnd[r] && r < model.interfaces.length) {
    const itf = model.interfaces[r];
    if (itf.phi.type === 'capacitive' && !s.diffList.some((d) => d.f === r && d.side === 0)) k += grid.area[g] * itf.phi.C * VT;
  }
  if (faces && g === grid.regionStart[r] && r > 0) {
    const itf = model.interfaces[r - 1];
    if (itf.phi.type === 'capacitive' && !s.diffList.some((d) => d.f === r - 1 && d.side === 1)) k += grid.area[grid.regionEnd[r - 1]] * itf.phi.C * VT;
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
    const slope = bulk ? 0 : v * (s.qgSlope[g] ?? 0) * VT + boundarySlope(s, g, syncFaces);
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

// The capacitive faces whose D follows φ̂ on both sides between two neighbouring nodes of the list
// (synced: see gaussRows): their coupling, −A C V_T, in each other's rows.
function faceCoupling(s, nodes, syncFaces) {
  const out = new Float64Array(nodes.length), grid = s.model.grid;
  if (!syncFaces) return out;
  for (let a = 0; a + 1 < nodes.length; a++) {
    const g = nodes[a].g;
    if (nodes[a + 1].g !== g + 1 || nodes[a].bulk || nodes[a + 1].bulk) continue;
    const r = grid.nodeRegion[g];
    if (g !== grid.regionEnd[r] || r >= s.model.interfaces.length) continue;
    const itf = s.model.interfaces[r];
    if (itf.phi.type === 'capacitive' && !s.diffList.some((d) => d.f === r)) out[a] = -grid.area[g] * itf.phi.C * s.VT;
  }
  return out;
}

// ∂c_i/∂η_j at node g (K for a non-ideal node, c_i δ_ij for an ideal one).
const Kat = (s, g, i, j) => (s.nodeIdeal[g] ? (i === j ? s.c[g * s.n + i] : 0) : s.K[g * s.n * s.n + i * s.n + j]);

// Solve a·x = b in place (a is m×m, row-major; partial pivoting); a pivot below `tiny` of the
// largest leaves its unknown at 0 (a direction nothing in the system holds).
function denseSolve(a, b, m) {
  const x = new Float64Array(m), perm = Int32Array.from({ length: m }, (_, k) => k);
  let big = 0;
  for (let k = 0; k < m * m; k++) big = Math.max(big, Math.abs(a[k]));
  const dead = new Uint8Array(m);
  for (let k = 0; k < m; k++) {
    let p = k;
    for (let r = k + 1; r < m; r++) if (Math.abs(a[r * m + k]) > Math.abs(a[p * m + k])) p = r;
    if (p !== k) {
      for (let c = 0; c < m; c++) [a[k * m + c], a[p * m + c]] = [a[p * m + c], a[k * m + c]];
      [b[k], b[p]] = [b[p], b[k]];
      [perm[k], perm[p]] = [perm[p], perm[k]];
    }
    if (!(Math.abs(a[k * m + k]) > 1e-13 * big)) {
      dead[k] = 1;
      continue;
    }
    for (let r = k + 1; r < m; r++) {
      const f = a[r * m + k] / a[k * m + k];
      if (f === 0) continue;
      for (let c = k; c < m; c++) a[r * m + c] -= f * a[k * m + c];
      b[r] -= f * b[k];
    }
  }
  for (let k = m - 1; k >= 0; k--) {
    if (dead[k]) {
      x[k] = 0;
      continue;
    }
    let t = b[k];
    for (let c = k + 1; c < m; c++) t -= a[k * m + c] * x[c];
    x[k] = t / a[k * m + k];
  }
  return x;
}

// A tridiagonal solve (diag, the symmetric off-diagonal `off` between a and a+1) of rhs, in place.
function thomas(diag, off, rhs) {
  const N = diag.length, cp = new Float64Array(N), x = rhs;
  let den = diag[0];
  for (let a = 0; a < N; a++) {
    if (a > 0) {
      den = diag[a] - off[a - 1] * cp[a - 1];
      x[a] -= off[a - 1] * x[a - 1];
    }
    cp[a] = a + 1 < N ? off[a] / den : 0;
    x[a] /= den;
  }
  for (let a = N - 2; a >= 0; a--) x[a] -= cp[a] * x[a + 1];
  return x;
}

// Every strictly neutral node onto its Gauss row together with what holds each closed stretch's
// amount and each floating capacitor's charge, by Newton on all of it at once:
//   unknowns: d at each node (by its mode: 'level', φ̂ moved at its levels; 'jump', φ̂ moved with
//   each η by z d, its concentrations kept), ℓ for each closed stretch (its species' η shifted by
//   ℓ over the stretch, but at a node that jumps), V for each floating terminal;
//   equations: G = 0 at each node, each closed stretch's amount its `want`, each floating
//   terminal's charge its `want`.
// The nodes' block is tridiagonal (a capacitive face synced on both sides couples its two boxes),
// the rest a border solved by its Schur complement. Each projection by turns met its own target and
// undid another's: the amounts and the small capacitors chased each other at their ratio, 1e-5 a
// round; a floating gate 500 times the held one beside it at 0.998 a round; and where a trace open
// ion's regridded charge contradicted the closed amounts there was no fixed point at all (Neon's
// fourth review, R4 1.3 and 1.5). A state that holds already is left bit for bit.
//   closed: [{ st, want }]; floating: [{ k, want }] (k a terminal: 0, 1 the contacts, 2 + a port).
// Returns { worst } (the worst of the three, each relative to its scale), marking moved nodes.
export function solveNeutral(s, nodes, { closed = [], floating = [], moved = null, syncFaces = false, closedAtNodes = false, tol = 1e-14, rounds = 30 } = {}) {
  const { n, M, z, present, VT } = s, F = FARADAY, grid = s.model.grid, N = nodes.length;
  if (N === 0 && closed.length === 0) return { worst: 0 };
  const posOf = new Int32Array(s.nNodes).fill(-1);
  nodes.forEach((nd, a) => (posOf[nd.g] = a));
  const off = faceCoupling(s, nodes, syncFaces);
  const K = closed.length, m = K + floating.length;
  const charge = (k) => (k >= 2 ? s.portQ[k - 2] : k === 0 ? s.contactD.left : -s.contactD.right);
  // (a closed stretch's ℓ acts where its species is, but not at a conductor node or one that jumps)
  const acts = (g, i) => s.nodeConductor[g] < 0 && present[g * n + i] && (posOf[g] < 0 ? !closedAtNodes : nodes[posOf[g]].mode !== 'jump');
  let worst = Infinity, best = Infinity;
  for (let it = 0; it <= rounds; it++) {
    const rows = gaussRows(s, nodes, { syncFaces });
    // (a box whose charge nothing holds by its φ̂, a gate taken away over it, can't jump: its ions
    // move, at its levels; a jump's G is a gate's charge, which can be 1e-4 of the ions' that the
    // scale holds, so it's met ten times closer)
    if (it === 0) nodes.forEach((nd, a) => nd.mode === 'jump' && !(rows[a].slope > 0) && (nd.mode = 'level'));
    const rG = rows.map((r) => r.G);
    let wG = 0;
    rows.forEach((r, a) => (wG = Math.max(wG, (Math.abs(r.G) / Math.max(r.scale, 1e-300)) * (nodes[a].mode === 'jump' ? 10 : 1))));
    const rm = new Float64Array(m), sm = new Float64Array(m);
    closed.forEach(({ st, want }, k) => {
      rm[k] = s.amount(st) - want;
      sm[k] = Math.abs(want) || 1e-300;
    });
    const q0 = floating.map(({ k }) => charge(k));
    floating.forEach(({ want }, f) => (rm[K + f] = q0[f] - want));
    // ∂G/∂V and ∂q/∂V for each floating terminal: by moving it (σ is linear in V)
    const dGdV = [], dqdV = new Float64Array(floating.length);
    floating.forEach(({ k }, f) => {
      const V0 = s.termV[k], dV = 1e-3 * Math.max(1, Math.abs(V0));
      s.termV[k] = V0 + dV;
      const r1 = gaussRows(s, nodes, { syncFaces });
      dqdV[f] = (charge(k) - q0[f]) / dV;
      dGdV.push(Float64Array.from(r1, (r, a) => (r.G - rG[a]) / dV));
      s.termV[k] = V0;
      sm[K + f] = Math.abs(floating[f].want) + Math.abs(dqdV[f]) * 1e-6 || 1e-300;
    });
    if (floating.length) gaussRows(s, nodes, { syncFaces }); // (the readouts at the state again)
    let wm = 0;
    for (let j = 0; j < m; j++) wm = Math.max(wm, Math.abs(rm[j]) / sm[j]);
    worst = Math.max(wG, wm);
    if (it === 0 && wG <= 1e-12 && wm <= tol) break; // (holds already: left as it is)
    if (worst <= tol || it === rounds || (it > 3 && !(worst < 0.5 * best) && worst <= 1e-12)) break;
    best = Math.min(best, worst);
    // The node block: ∂G_a/∂d_a, and the faces between neighbours.
    const diag = new Float64Array(N), Bm = Array.from({ length: m }, () => new Float64Array(N)), Cm = Array.from({ length: m }, () => new Float64Array(N));
    const Dm = new Float64Array(m * m);
    for (let a = 0; a < N; a++) {
      const { g, mode, bulk } = nodes[a], v = grid.vol[g];
      let dIons = 0; // ∂(ions)/∂d at levels kept
      if (mode !== 'jump') for (let i = 0; i < n; i++) {
        if (!present[g * n + i] || z[i] === 0) continue;
        let Kz = 0;
        for (let j = 0; j < n; j++) if (present[g * n + j] && z[j] !== 0) Kz += Kat(s, g, i, j) * z[j];
        dIons -= v * F * z[i] * Kz;
      }
      diag[a] = (bulk ? 0 : rows[a].slope) - dIons;
    }
    closed.forEach(({ st }, k) => {
      const i = st.species;
      for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
        if (!acts(g, i)) continue;
        const v = grid.vol[g], a = posOf[g];
        // the amount row: by ℓ of each closed stretch over the same node, and by d there
        closed.forEach(({ st: s2 }, k2) => {
          if (g >= s2.nodes[0] && g <= s2.nodes[1] && acts(g, s2.species)) Dm[k * m + k2] += v * Kat(s, g, i, s2.species);
        });
        if (a < 0) continue;
        let q = 0, Kz = 0;
        for (let j = 0; j < n; j++) if (present[g * n + j] && z[j] !== 0) {
          q += z[j] * Kat(s, g, j, i);
          Kz += Kat(s, g, i, j) * z[j];
        }
        Bm[k][a] = -v * F * q;
        if (nodes[a].mode !== 'jump') Cm[k][a] = -v * Kz;
      }
    });
    floating.forEach((_, f) => {
      const j = K + f;
      for (let a = 0; a < N; a++) {
        Bm[j][a] = dGdV[f][a];
        Cm[j][a] = VT * dGdV[f][a];
      }
      Dm[j * m + j] = dqdV[f];
    });
    // (a node nothing holds — no ion moves it, no capacitor — stays where it is)
    const dead = Uint8Array.from(diag, (x) => (x > 0 ? 0 : 1));
    for (let a = 0; a < N; a++) if (dead[a]) {
      diag[a] = 1;
      rG[a] = 0;
      if (a > 0) off[a - 1] = 0;
      if (a < N - 1) off[a] = 0;
      for (let j = 0; j < m; j++) Bm[j][a] = 0;
    }
    // Schur: d = T⁻¹(−rG − B y); (D − C T⁻¹ B) y = −rm − C T⁻¹(−rG)
    const x0 = thomas(diag, off, rG.map((x) => -x));
    const X = Bm.map((col) => thomas(diag, off, Float64Array.from(col)));
    const S = Float64Array.from(Dm), rhs = Float64Array.from(rm, (x) => -x);
    for (let j = 0; j < m; j++) {
      for (let a = 0; a < N; a++) rhs[j] -= Cm[j][a] * x0[a];
      for (let k = 0; k < m; k++) {
        let t = 0;
        for (let a = 0; a < N; a++) t += Cm[j][a] * X[k][a];
        S[j * m + k] -= t;
      }
    }
    const y = m > 0 ? denseSolve(S, rhs, m) : new Float64Array(0);
    const d = Float64Array.from(x0, (x, a) => {
      let t = x;
      for (let k = 0; k < m; k++) t -= X[k][a] * y[k];
      return dead[a] ? 0 : t;
    });
    // (a Newton step of at most 2 V_T in any potential or level: the ions are exponential in them)
    let big = 0;
    for (let a = 0; a < N; a++) big = Math.max(big, Math.abs(d[a]));
    for (let k = 0; k < K; k++) big = Math.max(big, Math.abs(y[k]));
    const al = big > 2 ? 2 / big : 1;
    for (let a = 0; a < N; a++) {
      if (d[a] === 0) continue;
      const { g, mode } = nodes[a], b = s.blockOfNode[g] * M;
      s.u[b] += al * d[a];
      if (mode === 'jump') for (let i = 0; i < n; i++) if (present[g * n + i] && z[i] !== 0) s.u[b + 1 + i] += z[i] * al * d[a];
      if (moved) moved[g] = 1;
    }
    closed.forEach(({ st }, k) => {
      if (y[k] === 0) return;
      for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
        if (!acts(g, st.species)) continue;
        s.u[s.blockOfNode[g] * M + 1 + st.species] += al * y[k];
        if (moved) moved[g] = 1;
      }
    });
    floating.forEach(({ k }, f) => {
      const dV = al * y[K + f];
      if (dV === 0) return;
      s.termV[k] += dV;
      if (moved && k >= 2) for (const g of s.model.ports[k - 2].nodes) moved[g] = 1;
    });
    s.computeConcentrations();
  }
  return { worst };
}
