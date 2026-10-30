// The equations at one time, assembled: every row's residual and its Jacobian, for Newton's
// iteration (src/newton.js), the time steps' bookkeeping and the small signal.
//
// CONTRACT
//
// These are methods of the Solver (src/solver.js installs them on its prototype), as are the
// kernels they call: the regions' in src/regions.js and the faces' and terminals' in
// src/boundaries.js, each with its own contract under this one. The layout they write into
// (blocks, slots, row o = block·M + slot, its compact index R = rix[o], loc[o] < 0 for an
// inactive slot, u + uLo, terminals 0, 1, 2 + k) is in docs/internals.md.
//
// The discretisation. Every region boundary is a doubled grid node (one per side) with a
// zero-volume flux block between them (the face's), so the Jacobian stays block-tridiagonal. A
// node's slots: φ̂ = Fφ/RT (in a metal, its carrier's flux J along the segment to its right),
// η_i = μ̄_i/RT, and at a port's nodes its coverages and gates. A face's: its displacement D, each
// linked species' flux N_i (toward +x, per area), each of its reactions' rate, each gate's
// fraction. Geometry is in the grid (src/grid.js): a node's box vol[g] = ∫A dx, the area A at a
// node, a segment's h = length/segArea. The rows are totals over a box (mol/s and C/s, or C):
//   Gauss (φ̂ row):   Σ D out of the box − vol·ρ = 0, ρ = F Σ z_i c_i + ρ_fixed (+ a port's
//                     capacitive charge, a face's sheet charge)
//   balance (η_i):    vol·(c_i − c_i,old)/dt + Σ N_i out − Σ N_i in + vol·Σ ν r (consumed: a
//                     reactant's ν > 0) − (made) − (a port's or contact's supply) = 0
//   a face's slot:    unknown − its law (D, N_i, a rate), or the law itself (a pinned jump, μ̄
//                     continuous), or a gate's (x − x_old)/dt − α(1 − x) + βx
// A flux along a segment is computed once and added with opposite signs to its two nodes, so
// sums of balance rows over boxes telescope exactly. At dt = Infinity there is no storage term
// (a steady system), and a gate's or surface's storage drops out the same way.
//
// assemble(dt), in order:
//   1. sys cleared and res zeroed; dtNow = dt; in a steady solve at dt = ∞, snapFlat (src/lawRows.js).
//   2. computeConcentrations (c, K, ex, coverages, from u + uLo).
//   3. Each region by its kind: a conductor's nodes and segments; else its nodes (dilute or
//      concentrated), its bulk reactions, then per segment the displacement, eddy mixing, a
//      concentrated region's fluxes; a dilute region's fluxes last.
//   4. Each face (_face).
//   5. The terminals (_terminals): the ports, the contacts, the floating ones' circuit rows.
//   6. applyLaws(dt) (src/lawRows.js): rows the laws take in place of the rows above.
//   7. On a time step with strictly neutral nodes (combining), _chargeRows: those nodes' rows
//      and unknowns transformed (transformed = true; Newton maps the update back).
//
// What it writes:
//   res (compact): each active row's residual as above; sys: its Jacobian, ∂res/∂(unknowns),
//     by _j into the three blocks; inactive rows and columns go nowhere (rix = N, the sink).
//   termB[k]: ∂res/∂V_k (per volt) for terminal k; termC[k]: ∂I_k/∂(unknowns); termI[k]: the
//     current into the device through terminal k (A, planar: A/m²); termDI[k], termDX[k][j]:
//     ∂ of k's circuit residual by V_k and V_j; termRes[k]: a floating k's circuit residual. Kept
//     exactly only for _liveTerms() (the floating terminals, or all where allTerminals).
//   The terminals' readouts (src/boundaries.js): contactFlux, contactD, portFlux, portQ, qg,
//     qgSlope, portRates, contactQ; and s.pins (applyLaws).
// Where this.lin is set (only by _assembleDifference), the dilute kernels' Jacobian terms go into
// it, in difference form (src/difference.js), and the rest into sys, which is then not Newton's.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. Every Jacobian entry is the derivative of the residual as assembled: in (φ̂, η, the face
//      and port slots), by V_k in termB, and in the transformed unknowns after 7.
//      (test/jacobian.test.js checks it against central differences, on the devices it lists.)
//   2. res, sys and the live terminals' terms depend only on the state (u, uLo, termV), dt, the
//      sources' time and what a step or solve set (cOld, the gates' and coverages' old values,
//      the charges at the step's start, steady, combining, the laws' references), never on an
//      earlier assembly's leftovers: assemble(dt) twice at one state gives the same, bit for bit.
//      (A held terminal's termB and termC aren't cleared, so they aren't kept.)
//   3. A kernel adds into res and sys (+=) and sets a row only where the row is its own (a face's
//      slots, a metal's flux slot, a port's gates) or it replaces the row whole (a contact's or
//      a port's held level, a contact's φ law; the laws in 6), and no kernel adds to a replaced
//      row after it is replaced. A row read as a flux before it is replaced (a contact's, a held
//      level's) has every other kernel's terms in it by then.
//   4. _assembleDifference(dt): res and the terminal terms as assemble(dt) gives them; its
//      matrix plus its difference terms is assemble(dt)'s Jacobian (J·v equal to round-off);
//      this.sys and this.lin as they were.
//   5. _assembleBookkeeping(dt): contactFlux, contactD, portFlux, portQ and termI bit for bit
//      as assemble(dt) would set them at the same state (one a steady assembly has already
//      snapped flat: it doesn't snap), with sys, res, termB, termC and termDI elsewhere left
//      partial, not to be read; it doesn't recompute concentrations (they must be current).
// Tests: every test (every solve assembles); test/jacobian.test.js (invariant 1), test/layout.test.js
// (difference form, the bookkeeping's readouts).

import { BlockTridiagonal } from './blockTridiagonal.js';
import { DifferenceTerms } from './difference.js';
import { applyLaws, snapFlat } from './lawRows.js';

export class Assembly {
  // Add ∂(row rb, slot rs)/∂(block cb, slot cs) to the Jacobian (nothing if either slot isn't
  // an unknown).
  _j(rb, rs, cb, cs, v) {
    const M = this.M, r = this.loc[rb * M + rs], c = this.loc[cb * M + cs];
    if (r < 0 || c < 0) return;
    const sys = this.sys, sz = sys.sizes;
    if (cb === rb) sys.B[sys.offB[rb] + r * sz[rb] + c] += v;
    else if (cb === rb - 1) sys.A[sys.offA[rb] + r * sz[cb] + c] += v;
    else if (cb === rb + 1) sys.C[sys.offC[rb] + r * sz[cb] + c] += v;
    else throw new Error(`internal: non-tridiagonal coupling ${rb}→${cb}`);
  }

  /** Assemble residual and Jacobian for a backward-Euler step of size dt. */
  assemble(dt) {
    const { model, res, sys } = this;
    const { grid, materials, regions } = model;
    sys.clear();
    res.fill(0);
    this.dtNow = dt;
    if (this.lin) this.lin.reset(1 / dt);
    if (this.steady && dt === Infinity) snapFlat(this);
    this.computeConcentrations();

    // Regions, each assembled by the kernel for its kind: a conductor (its carrier only), a
    // dilute region (ideal statistics, the fast path) or a concentrated one (any statistics).
    for (let r = 0; r < regions.length; r++) {
      const reg = regions[r], mat = materials[reg.material];
      const g0 = grid.regionStart[r], g1 = grid.regionEnd[r];
      if (mat.conductor) {
        for (let g = g0; g <= g1; g++) this._nodeConductor(g, dt);
        for (let s = g0; s < g1; s++) this._segmentConductor(s, s + r, s + r + 1, mat, grid.segLength[s] / grid.segArea[s]);
        continue;
      }
      if (mat.ideal) this._nodesDilute(g0, g1, dt);
      else for (let g = g0; g <= g1; g++) this._nodeConcentrated(g, dt);
      const rxs = this.rxsIn[reg.material];
      if (rxs.length > 0) for (let g = g0; g <= g1; g++) this._bulkReactions(g, rxs);
      for (let s = g0; s < g1; s++) {
        const bL = s + r, bR = bL + 1, h = grid.segLength[s] / grid.segArea[s]; // (a length over the cross-section: fluxes are totals)
        this._segmentDisplacement(s, bL, bR, mat, h);
        if (reg.mixing > 0) this._segmentMixing(s, bL, bR, reg.mixing, h, mat);
        if (!mat.ideal) this._segmentConcentrated(s, bL, bR, mat, h, reg.velocity);
      }
      if (mat.ideal) this._segmentsDilute(g0, g1, r, mat, reg.velocity);
    }

    // Faces: the flux node carries D, the linked species' fluxes and the reaction rates; and a
    // face's surface its coverages' storage (a port's comes with its terminal's terms).
    for (let f = 0; f < this.nFaces; f++) this._face(f);
    if (Number.isFinite(dt)) for (const sf of model.surfaces) if (sf.face >= 0) this._surfaceStorage(sf, dt);

    // Terminals: ports (after every other term at their nodes, so a held level can read its
    // flux), then contacts, then the circuit rows of the floating ones.
    this._terminals(dt);
    // The laws' rows (src/lawRows.js): a time step's circuit laws; a steady system's levels held
    // flat, affinities, conserved amounts, islands' balances, and the pins held aside of them.
    applyLaws(this, dt);
    this.transformed = this.combining && dt !== Infinity;
    if (this.transformed) this._chargeRows(dt);
  }

  // Zero one row of the Jacobian (all three blocks), ready to be replaced.
  _replaceRow(b, rs) {
    const l = this.loc[b * this.M + rs];
    if (l < 0) return;
    if (this.lin) this.lin.kill(this.rix[b * this.M + rs]);
    const { sys, nB } = this, sz = sys.sizes, m = sz[b];
    sys.B.fill(0, sys.offB[b] + l * m, sys.offB[b] + (l + 1) * m);
    if (b > 0) sys.A.fill(0, sys.offA[b] + l * sz[b - 1], sys.offA[b] + (l + 1) * sz[b - 1]);
    if (b < nB - 1) sys.C.fill(0, sys.offC[b] + l * sz[b + 1], sys.offC[b] + (l + 1) * sz[b + 1]);
  }

  // The Jacobian at dt with the dilute kernels' terms kept aside in difference form (see
  // DifferenceTerms), the rest assembled into a matrix of its own. The residual and terminal terms
  // come out as from assemble(dt); the Newton matrix is left as it was.
  _assembleDifference(dt) {
    const keep = this.sys;
    const rest = this.sysRest ?? (this.sysRest = new BlockTridiagonal(keep.n, keep.sizes));
    const n = this.n, nN = this.nNodes;
    const lin = this.diffTerms ?? (this.diffTerms = new DifferenceTerms(keep.size, nN * n * n, nN * n, 2 * nN * n * n));
    this.sys = rest;
    this.lin = lin;
    try {
      this.assemble(dt);
    } finally {
      this.sys = keep;
      this.lin = null;
    }
    return { rest, lin };
  }

  // The contact and port readouts (fluxes, displacements, the last segment's current) at the
  // current state, from only the boxes they're read from. Each of those boxes gets the same terms
  // in the same order as in assemble(), so the readouts are identical; other residuals and the
  // Jacobian are left partial (the next assemble() starts afresh). Concentrations must be current.
  _assembleBookkeeping(dt) {
    const { model } = this;
    this.dtNow = dt;
    const { grid, materials, regions } = model;
    this.res.fill(0);
    for (let r = 0; r < regions.length; r++) {
      const { nodes, segs } = this.bookkeeping[r];
      if (nodes.length === 0 && segs.length === 0) continue;
      const reg = regions[r], mat = materials[reg.material];
      if (mat.conductor) {
        for (const g of nodes) this._nodeConductor(g, dt);
        for (const s of segs) this._segmentConductor(s, s + r, s + r + 1, mat, grid.segLength[s] / grid.segArea[s]);
        continue;
      }
      for (const g of nodes) {
        if (mat.ideal) this._nodesDilute(g, g, dt);
        else this._nodeConcentrated(g, dt);
      }
      const rxs = this.rxsIn[reg.material];
      if (rxs.length > 0) for (const g of nodes) this._bulkReactions(g, rxs);
      for (const s of segs) {
        const bL = s + r, bR = bL + 1, h = grid.segLength[s] / grid.segArea[s]; // (a length over the cross-section: fluxes are totals)
        this._segmentDisplacement(s, bL, bR, mat, h);
        if (reg.mixing > 0) this._segmentMixing(s, bL, bR, reg.mixing, h, mat);
        if (!mat.ideal) this._segmentConcentrated(s, bL, bR, mat, h, reg.velocity);
      }
      if (mat.ideal) for (const s of segs) this._segmentsDilute(s, s + 1, r, mat, reg.velocity);
    }
    for (const f of this.bookkeepingFaces) this._face(f);
    this._terminals(dt);
  }
}
