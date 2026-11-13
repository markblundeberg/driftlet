// The laws in the Newton system: each law src/laws.js finds, evaluated on the state and written
// into the equations in place of the row it takes, and the linear solve that holds the rows too
// dense for the block-tridiagonal matrix aside.
//
// CONTRACT
//
// These are functions of a Solver `s` (src/solver.js), which they read and write; its layout
// (blocks, slots, row o = block·M + slot, its compact index R = s.rix[o], loc[o] < 0 for an
// inactive slot, u + uLo, terminals 0, 1, 2 + k) is in docs/internals.md. The laws are as
// Solver._laws() set them from findLaws: s.constraints, affinities, localConstraints, islands,
// timeLaws, contactLaw, and s.flatStretches, each meaning what src/laws.js's contract says. The
// solver's own helpers they call: s._j(rb, r, cb, c, v) adds v at row (rb, r), column (cb, c) of
// the block matrix s.sys; s._replaceRow(b, r) clears that row (all three blocks, its terms in
// difference form, s.lin, its terminal terms, termB, and its residual); s._dc(g, i, d) puts ∂c_i/∂(slot) at node g into d; s._solveLinear(rhs,
// out) solves the factorised s.sys; s._assembleBookkeeping(dt) sets the terminals' charges and
// currents (portQ, contactD, termI); s.amount(st), s.surfaceAmount(col) a stretch's or surface
// column's amount.
//
// The Newton system (see Solver.newton): J δ = res, the update u −= δ (and termV[k] −= δV for a
// floating k). It is s.sys (the blocks, residual s.res), and beside it, for each floating
// terminal k, its column termB[k] (∂ of every row by V_k), its circuit row termC[k], termRes[k]
// and its ∂/∂V's, termDI[k] (own) and termDX[k][j] (terminal j's); and s.pins.
//
// A pin: a law too dense for the block matrix, { row, idx, w, len, res } (and `col`, see below).
// Its row in s.sys is the identity with residual 0, and the law, Σ_j w_j δ[idx_j] = res, is
// solved beside it by solveBordered: the pin's multiplier is free, so whatever else that row
// holds (a terminal's termB entry, say) is absorbed by it. Not scaled by equilibrate() (src/newton.js),
// which scales s.sys's rows, s.res and termB, and runs between assembly and the solve.
//
// The functions, and when each runs:
//   captureLaws(s)   when a steady solve begins (Solver.solveSteady): s.lawRefs, the amount each
//                    constraint fed by a terminal driven by a current holds now (by its key);
//                    s.chargeRefs, the charge (portQ, or the contact's ±D) of each terminal
//                    driven by a current that passes current only by charging. (A steady solve
//                    has such a terminal only at I = 0: Solver._solveSteadyAll stops at any
//                    other.) Calls _assembleBookkeeping where it needs the charges.
//   captureLocal(s)  when a steady solve's steady system begins (Solver._steadyFromHere): each
//                    local constraint's reference[j], its sum at node nodes[j] now.
//   snapFlat(s)      in a steady assembly (dt = ∞), before the concentrations: each flat
//                    stretch's level (u and uLo) set at every node to its first node's.
//   applyLaws(s, dt) at the end of Solver.assemble(dt), after every kernel and the terminals:
//                    dt finite: each circuit row a law of the constraints or timeLaws may read by
//                    the change of its amount (circuitLaws). Steady (s.steady and dt = ∞): the
//                    rows of the flat stretches, affinities, constraints, local constraints and,
//                    once Newton turns them on (s.islandsOn), the islands, in that order; and
//                    s.pins, those of them held aside (constraints in a balance row, affinities,
//                    islands). Otherwise s.pins is empty.
//   solveBordered(s, rhs, delta, deltaV, pins)  each Newton solve: δ (delta, compact) and the
//                    floating terminals' δV (deltaV, by s.floating's order) with s.sys
//                    factorised.
//
// What each row says, as written (the steady rows replace a row whole, with _replaceRow first):
//   flat stretch     at each of its nodes, u_o − level = 0, Jacobian 1 on itself (the level:
//                    its contact's, the terminal's ∂/∂V in termB; reaching none, its first
//                    node's own u, so that row is 0 = 0, and the stretch's amount, a constraint
//                    with `col` = its nodes' rows, moves it whole: a unit pin there shifts every
//                    row of `col`). A face inside it whose link is an equilibrium: its flux
//                    unknown = 0.
//   affinity         pin: fixedA − Σ ν η at the face's two nodes = 0, ∂ = −ν.
//   constraint       in a balance row, a pin: amount − ref = 0 (ref: its reference without
//                    terminals, else lawRefs; with no reference captured there, none, the amount
//                    held where it is), its ∂ from lawAmount. A stretch's amount is Σ vol·c over
//                    its nodes and A_f Γ_i of each face's diffuse layer at them (s.gam; its ∂
//                    by the layer's drop and its edge node's slots). In terminal k's circuit row
//                    (in place of what Solver._terminals wrote): the same, its contacts' charges
//                    λ·D_in/F added from s.contactQ[t] (as Solver._contact captured them in this
//                    assembly), termDI and termDX from the electrodes' charges' ∂/∂V.
//   local constraint at each node of `nodes`, Σ w c_i (+ Σ w Γ a θ_q where port k's surface is)
//                    − reference = 0, its ∂ in the node's block (a coverage's,
//                    θ_q (δ_qx − θ_x) over its surface's species).
//   island           pin: Σ weight·u over the `flux` rows (faces' fluxes, and face reactions'
//                    rates) + Σ vol·ν·r over its `bulk` entries (each node's rate and its ∂ from
//                    the solver's _bulkRateInto, as that node's row has them) + Σ the `outside`
//                    rows' residuals − Σ I/(zF) over `currents` (a
//                    metal's floating contacts: I_set, or (V_src − V)/R, its ∂/∂V in `termSlope` for
//                    solveBordered) = 0, its ∂ the same sums of their
//                    Jacobian rows (rowInto), gathered by column. A cut whose reach has a segment
//                    conducting the species under 1e-6 of its first (∂N/∂η_L at the state) takes
//                    the least of those before the conductance has risen 1e6 times again (the
//                    first bottleneck outward: past it another pool begins): that face's flux and
//                    outside row are left out, and sign·N of the
//                    segment (its Scharfetter–Gummel flux, exactly as _segmentsDilute's) is in,
//                    with its ∂ by both ends' η and φ̂, and vol·ν·r of each `span` piece's `rxs` at
//                    each node between the face and that segment. The segments cuts take are chosen
//                    once a round (chooseCuts: cutAt, cleared by Newton at a call's start and where
//                    it measures the imbalance) and held through it, and with them the islands the
//                    round leaves out (`absorbed`, their rows their nodes' balances): one summed
//                    past its piece whose own node conducts under 1e-6 of the best in its span (its
//                    row isn't in the pool it sums), and, most conductive first, one whose node a
//                    standing island of its stretch already sums. Islands that share a rate are
//                    then recombined (separateRates) so that each rate is in one pin only, the
//                    fastest (k0) first: the same constraints, each still in its own row.
//   circuit law      (time steps) terminal k's circuit row, where it reads by the change: F (S −
//                    S*)/(λ dt) − (its drive's current: the source, or (V_src − V)/R), with S*
//                    the amount at the step's start (lawStart), chosen once per step per law
//                    (s.lawForm) by comparing the round-off of the two readings; only where k
//                    is the law's one current-driven terminal passing a current (with more, the
//                    change is theirs together, and k reads its flux).
//
// Invariants (each says what a reviewer should be able to check here):
//   1. A row is replaced only where its slot is active (loc ≥ 0), and every entry a row held
//      before is gone: _replaceRow clears the three blocks, s.lin's, its termB entries (where a
//      terminal's ∂/∂V was in it: a port reaction's, in the slots it touches) and its residual,
//      pins included.
//   2. Each law's residual and its derivatives are of the same function of the state: every
//      term in amount (lawAmount, and the rows above) has its ∂ in idx/w, dV or the block row,
//      and none is counted twice (a metal's carriers by its charged face's D, not its η).
//   3. lawAmount at a stretch, surface or electrode sums the same quantity lawStart sums at the
//      step's start, and captureLaws at a solve's start: what a step's balance rows change is
//      what a circuit law reads.
//   4. solveBordered solves the full system: the block rows, each pin's law, and each floating
//      terminal's circuit row, with δV, by the Woodbury identity written at it, exact given the
//      block solves (T⁻¹ applied by the factorised s.sys): its accuracy is theirs.
//   5. Nothing here writes the laws' structure (their parts, rows, nodes): only the workspace
//      fields src/laws.js leaves to the solver (idx, w, len, res, dV, reference; an island's
//      deep and cutAt, the segment each cut took or −1, absorbed, fluxScale, Σ|flux| over its
//      faces not cut and its cuts' segments, cutSlope, Σ|∂N/∂η| over those segments, bulkScale,
//      Σ|vol·ν·r| over its rates,
//      and termSlope; idx and w grown where its row needs more room), the solver's maps lawRefs,
//      chargeRefs, lawForm, s.pins and s.pinnedIslands, scratch arrays kept on s (islandRow,
//      islandDr, linRow, dWork, borderCols), and, in snapFlat alone, the state (u, uLo).
// Tests: test/laws.test.js, test/steady.test.js and every test that solves a steady state;
// tools/lawsweep.js (solve() against a settled advance()).

import { FARADAY } from './constants.js';
import { sourceAt } from './device.js';
import { carrierZ, chargeOnly } from './laws.js';
import { SolverError } from './errors.js';

// The laws' rows at the end of an assembly (see the contract above).
export function applyLaws(s, dt) {
  if (Number.isFinite(dt) && (s.constraints.length > 0 || s.timeLaws.length > 0)) circuitLaws(s, dt);
  const steady = s.steady && dt === Infinity;
  if (steady) flatRows(s);
  if (steady && s.affinities.length > 0) applyAffinities(s);
  if (steady && s.constraints.length > 0) applyConstraints(s);
  if (steady && s.localConstraints.length > 0) applyLocalConstraints(s);
  s.pinnedIslands = steady && s.islandsOn ? applyIslands(s) : [];
  s.pins = steady ? [...s.constraints.filter((cs) => cs.circuit < 0), ...s.affinities, ...s.pinnedIslands] : [];
}

// What each combination fed by a terminal driven at no current holds now, which a steady solve
// keeps: a closed battery at open circuit stays at the charge it had.
export function captureLaws(s) {
  // A terminal at no current that passes current only by charging a capacitance keeps its charge.
  s.chargeRefs.clear();
  const charged = s.floating.filter((k) => s.terms[k].drive.kind === 'I' && chargeOnly(s, k));
  if (charged.length > 0) {
    s._assembleBookkeeping(Infinity);
    for (const k of charged) s.chargeRefs.set(k, k >= 2 ? s.portQ[k - 2] : k === 0 ? s.contactD.left : -s.contactD.right);
  }
  s.lawRefs.clear();
  if (s.constraints.some((cs) => cs.contacts.length > 0 || cs.carriers.some(({ k }) => s.model.ports[k - 2].capacitance)) && charged.length === 0) s._assembleBookkeeping(Infinity);
  for (const cs of s.constraints) {
    if (cs.terminals.length === 0) continue;
    let a = 0;
    for (const { stretch, w } of cs.parts) a += w * s.amount(s.stretches[stretch]);
    for (const { col, w } of cs.surface) a += w * s.surfaceAmount(col);
    for (const { k, w } of cs.carriers) if (s.model.ports[k - 2].capacitance) a += (w * s.portQ[k - 2]) / (carrierZ(s, k) * FARADAY);
    for (const { t, lambda } of cs.contacts) a += (lambda * (t === 0 ? s.contactD.left : -s.contactD.right)) / FARADAY;
    s.lawRefs.set(cs.key, a);
  }
}

// Each local constraint's weighted sum at each node, as the state holds it now: what a steady
// solve keeps.
export function captureLocal(s) {
  const { n, c, th, nSurf } = s;
  for (const lc of s.localConstraints) {
    lc.nodes.forEach((g, j) => {
      let t = 0;
      for (const { i, q, k, w } of lc.parts) {
        if (i !== undefined) t += w * c[g * n + i];
        else if (s.surfAt[g] === k) t += w * sites(s, g) * th[g * nSurf + q];
      }
      lc.reference[j] = t;
    });
  }
}

// The islands' summed balances (steady solves only): Σ A_f u over the faces into each, its
// row (its middle node's) pinned in the factorised matrix and the sum bordered (solveBordered) like a
// conserved amount's. The pieces of a flat stretch are left as they are (its rows are flatRows').
function applyIslands(s) {
  const { u, res } = s, R = s.rix, out = [];
  if (s.islands.some((isl) => !isl.cutAt)) chooseCuts(s);
  for (const isl of s.islands) {
    if (s.flatStretches.includes(isl.stretch) || isl.absorbed) continue;
    const b0 = Math.floor(isl.row / s.M), r0 = isl.row % s.M;
    if (s.loc[isl.row] < 0) continue;
    // Gathered by column (the faces' u, then each outside row's entries, which cancel those u
    // exactly where they meet), then the nonzeros kept.
    const acc = s.islandRow ?? (s.islandRow = new Float64Array(s.sys.size + 1));
    const cols = [];
    const add = (col, v) => {
      if (acc[col] === 0) cols.push(col);
      acc[col] += v;
    };
    let t = 0;
    // (a face cut deeper, at its end piece's bottleneck: that face's flux and outside row left
    // out, the segment's flux in. Chosen once a round, at the state the round starts from, and
    // held through it: re-chosen at each iterate, the summed nodes moved under Newton as it
    // walked the level (38, 41, 47, 48 in hunt A's base, listed in reverse), and an update of
    // 4e6 thermal units on the end piece's minority holes, held by nothing, ended the round.)
    const deep = isl.cutAt, skipFace = new Set(), skipRow = new Set(), { grid } = s.model;
    isl.deep = deep;
    deep.forEach((sg, k) => sg >= 0 && (skipFace.add(isl.cuts[k].face), skipRow.add(isl.cuts[k].outside)));
    let fluxScale = 0;
    isl.flux.forEach((o, j) => {
      if (s.loc[o] < 0 || skipFace.has(j)) return;
      t += isl.weight[j] * u[o];
      fluxScale += Math.abs(isl.weight[j] * u[o]);
      add(R[o], isl.weight[j]);
    });
    for (const o of isl.outside) {
      if (s.loc[o] < 0 || skipRow.has(o)) continue;
      t += res[R[o]];
      rowInto(s, o, add);
    }
    isl.cutSlope = 0;
    deep.forEach((sg, k) => {
      if (sg < 0) return;
      const N = segmentInto(s, isl.stretch.species, sg, grid.nodeRegion[sg], isl.cuts[k].sign, add);
      t += N;
      fluxScale += Math.abs(N);
      isl.cutSlope += Math.max(Math.abs(s.sg[2]), Math.abs(s.sg[3]));
    });
    isl.fluxScale = fluxScale;
    // (its bulk reactions, node by node: the piece's own, and an end piece's out to the
    // bottleneck, the nodes between it and the face; each +vol·ν·r, as in the node's row)
    let bulkScale = 0;
    const bulkInto = (g, rx, nu) => {
      const b = s.blockOfNode[g], w = grid.vol[g] * nu, dr = s.islandDr ?? (s.islandDr = new Float64Array(s.M));
      const r = s._bulkRateInto(rx, g, dr);
      t += w * r;
      bulkScale += Math.abs(w * r);
      for (let k = 0; k < s.M; k++) if (dr[k] !== 0 && s.loc[b * s.M + k] >= 0) add(R[b * s.M + k], w * dr[k]);
    };
    for (const [g, rx, nu] of isl.bulk ?? []) bulkInto(g, rx, nu);
    deep.forEach((sg, k) => {
      const cut = isl.cuts[k];
      if (sg < 0) return;
      for (const { region, rxs } of cut.span) {
        if (rxs.length === 0) continue;
        const g0 = cut.sign < 0 ? Math.max(sg + 1, grid.regionStart[region]) : grid.regionStart[region], g1 = cut.sign < 0 ? grid.regionEnd[region] : Math.min(sg, grid.regionEnd[region]);
        for (let g = g0; g <= g1; g++) for (const [rx, nu] of rxs) bulkInto(g, rx, nu);
      }
    });
    isl.bulkScale = bulkScale;
    // (a metal's: less what its floating contacts pass, I/(zF) each: I_set, or (V_src − V)/R,
    // whose ∂/∂V the bordered system takes as the pin's own, termSlope)
    isl.termSlope = null;
    for (const k of isl.currents ?? []) {
      const d = s.terms[k].drive, src = sourceAt(d.src, s.sourceTime, s.sourceBefore);
      if (d.kind === 'I') t -= src / (isl.z * FARADAY);
      else {
        t -= (src - s.termV[k]) / (d.R * isl.z * FARADAY);
        (isl.termSlope ??= new Map()).set(k, 1 / (d.R * isl.z * FARADAY));
      }
    }
    // (bulk reactions add a column per slot of each node they sum over: room grown here)
    if (cols.length > isl.idx.length) {
      isl.idx = new Int32Array(cols.length);
      isl.w = new Float64Array(cols.length);
    }
    let len = 0;
    for (const col of cols) {
      if (col < s.sys.size && acc[col] !== 0) {
        isl.idx[len] = col;
        isl.w[len++] = acc[col];
      }
      acc[col] = 0;
    }
    isl.len = len;
    isl.res = t;
    s._replaceRow(b0, r0);
    s._j(b0, r0, b0, r0, 1);
    res[R[isl.row]] = 0;
    out.push(isl);
  }
  separateRates(s, out);
  return out;
}

// A segment's flux of species i toward +x (Scharfetter–Gummel, as the assembly has it: N, then
// ∂N/∂Δ, ∂/∂η_L, −∂/∂η_R, in s.sg), segment sg of region r.
function segmentFlux(s, i, sg, r) {
  const { grid, regions } = s.model;
  return s._sgSegment(sg, sg + r, i, grid.segLength[sg] / grid.segArea[sg], regions[r].velocity);
}

// Where the pieces a cut reaches through conduct the species least (∂N/∂η_L, the conductance at
// the state), if that's under 1e-6 of the first segment's: the segment to sum the island out to,
// or −1. The first such bottleneck outward: past one, where the conductance has risen 1e6 times
// again, another pool begins (a second base, n/p⁺/n/p⁺/n), and its far side isn't this one's.
function bottleneck(s, isl, cut) {
  const { nodeRegion } = s.model.grid;
  let first = NaN, min = Infinity, at = -1;
  for (const sg of cut.reach) {
    const G = segmentFlux(s, isl.stretch.species, sg, nodeRegion[sg])[2];
    if (first !== first) first = G;
    if (G < min) [min, at] = [G, sg];
    else if (min < 1e-6 * first && G > 1e6 * min) break;
  }
  return min < 1e-6 * first ? at : -1;
}

// Each island's cuts for a round, at the state it starts from (cutAt), and which islands the
// round leaves out (absorbed): a piece whose first node another island's sum takes in, out to
// its bottleneck, is part of that island's pool, and its own sum would be the same one (the n
// layer between a floating base and a fourth layer: its holes' sum, reaching back through the
// base, has the base's two bottlenecks). The island whose piece conducts most stands, the pool
// itself; the rows of those it takes in stay their nodes' balances, as an end piece's do.
function chooseCuts(s) {
  const { nodeRegion } = s.model.grid, cutting = [];
  for (const isl of s.islands) {
    isl.cutAt = (isl.cuts ?? []).map((cut) => bottleneck(s, isl, cut));
    isl.absorbed = false;
    if (!isl.piece || s.flatStretches.includes(isl.stretch) || s.loc[isl.row] < 0) continue;
    let [g0, g1] = isl.piece, peak = 0;
    const G = (sg) => (nodeRegion[sg] === nodeRegion[sg + 1] ? segmentFlux(s, isl.stretch.species, sg, nodeRegion[sg])[2] : 0);
    for (let sg = g0; sg < g1; sg++) peak = Math.max(peak, G(sg));
    isl.cutAt.forEach((sg, k) => sg >= 0 && (isl.cuts[k].sign < 0 ? (g0 = sg + 1) : (g1 = sg)));
    // (summed out past its piece, the row it replaces must be in the pool it sums: its node
    // conducting within 1e-6 of the span's best. An n layer's holes beside a floating base,
    // summed through the base to the far bottleneck, replaced a row 1e29 times less conductive
    // than the base's; that node's own balance, left to the difference of the sum and the rest,
    // was lost, and the flux through it jumped 29-fold, the base 12 mV off with nothing said.)
    if (isl.cutAt.some((sg) => sg >= 0)) {
      let best = 0;
      for (let sg = g0; sg < g1; sg++) best = Math.max(best, G(sg));
      if (Math.max(G(isl.node - 1), G(isl.node)) < 1e-6 * best) {
        isl.absorbed = true;
        continue;
      }
    }
    cutting.push({ isl, g0, g1, peak });
  }
  cutting.sort((a, b) => b.peak - a.peak);
  const kept = [];
  for (const c of cutting) {
    const g = c.isl.node;
    if (kept.some((k) => k.isl.stretch === c.isl.stretch && g >= k.g0 && g <= k.g1)) c.isl.absorbed = true;
    else kept.push(c);
  }
}

// That segment's flux into an island's sum, sign·N (−1: the island to its right), its ∂ into
// `add` by column; returns the term.
function segmentInto(s, i, sg, r, sign, add) {
  const { M, z, loc } = s, R = s.rix, bL = sg + r, bR = bL + 1;
  const out = segmentFlux(s, i, sg, r), N = out[0], dNdd = out[1];
  const put = (o, v) => loc[o] >= 0 && v !== 0 && add(R[o], sign * v);
  put(bL * M + 1 + i, out[2]);
  put(bR * M + 1 + i, -out[3]);
  if (z[i] !== 0) {
    put(bL * M, -z[i] * dNdd - z[i] * N);
    put(bR * M, z[i] * dNdd);
  }
  return sign * N;
}

// Islands that share a face reaction's rate (a fast couple between a floating metal and the pool
// beside it, round 11's F5) are, in the bordered system, each dominated by its response to that
// rate, k0 c per thermal unit, and nearly dependent: the slow rates that tell them apart were lost
// in the comparison ("the bordered system is singular"). The same constraints in another basis:
// each rate kept in one sum only, the fastest first (Gauss–Jordan over the rate columns, rate by
// rate), so a fast couple drops out exactly from all but one (e⁻ + 3 B⁻: what it doesn't change).
// Each sum still replaces its own row; they span what the replaced rows did.
function separateRates(s, isls) {
  if (isls.length < 2) return;
  const R = s.rix, rows = isls.map((isl) => {
    const m = new Map();
    for (let j = 0; j < isl.len; j++) m.set(isl.idx[j], (m.get(isl.idx[j]) ?? 0) + isl.w[j]);
    return m;
  });
  const speed = new Map();
  for (const isl of isls) for (const [o, k] of isl.speed ?? []) if (s.loc[o] >= 0) speed.set(R[o], Math.max(speed.get(R[o]) ?? 0, k));
  const order = [...speed.keys()].sort((a, b) => speed.get(b) - speed.get(a));
  const pivoted = new Set();
  let changed = false;
  for (const col of order) {
    const holders = isls.map((_, q) => q).filter((q) => (rows[q].get(col) ?? 0) !== 0);
    const p = holders.find((q) => !pivoted.has(q));
    if (p === undefined || holders.length < 2) {
      if (p !== undefined) pivoted.add(p);
      continue;
    }
    pivoted.add(p);
    const wp = rows[p].get(col);
    for (const q of holders) {
      if (q === p) continue;
      const f = rows[q].get(col) / wp;
      for (const [c, w] of rows[p]) rows[q].set(c, (rows[q].get(c) ?? 0) - f * w);
      rows[q].delete(col); // (exactly)
      isls[q].res -= f * isls[p].res;
      for (const [k, v] of isls[p].termSlope ?? []) (isls[q].termSlope ??= new Map()).set(k, (isls[q].termSlope.get(k) ?? 0) - f * v);
      changed = true;
    }
  }
  if (!changed) return;
  isls.forEach((isl, q) => {
    const entries = [...rows[q]].filter(([, w]) => w !== 0);
    if (entries.length > isl.idx.length) [isl.idx, isl.w] = [new Int32Array(entries.length), new Float64Array(entries.length)];
    entries.forEach(([c, w], j) => {
      isl.idx[j] = c;
      isl.w[j] = w;
    });
    isl.len = entries.length;
  });
}

// Row o's Jacobian entries, as add(compact column, value): the assembled ones, and the dilute
// kernels' kept aside in difference form (see Solver._assembleDifference) when they are.
function rowInto(s, o, add) {
  const { sys, M, nB } = s, b = Math.floor(o / M), l = s.loc[o];
  const { A, B, C, sizes, offA, offB, offC, offX } = sys, m = sizes[b];
  for (let k = 0; k < m; k++) if (B[offB[b] + l * m + k] !== 0) add(offX[b] + k, B[offB[b] + l * m + k]);
  if (b > 0) for (let k = 0, mp = sizes[b - 1]; k < mp; k++) if (A[offA[b] + l * mp + k] !== 0) add(offX[b - 1] + k, A[offA[b] + l * mp + k]);
  if (b < nB - 1) for (let k = 0, mn = sizes[b + 1]; k < mn; k++) if (C[offC[b] + l * mn + k] !== 0) add(offX[b + 1] + k, C[offC[b] + l * mn + k]);
  if (s.lin) {
    const into = s.linRow ?? (s.linRow = new Float64Array(sys.size + 1));
    s.lin.capture(s.rix[o], 1, 0, into);
    for (let k = 0; k < sys.size; k++) {
      if (into[k] !== 0) add(k, into[k]);
      into[k] = 0;
    }
  }
}

// An electrode surface's sites per volume at node g, Γ a (mol/m³): what a coverage of 1 holds,
// beside the concentrations of a local constraint. (Summed as θ alone, an immobile S that adsorbs
// as A(ads) kept S + θ_A, mol/m³ and a fraction.)
function sites(s, g) {
  const k = s.surfAt[g];
  return s.model.surfaces[k].species[0].capacity * s.surfArea[k][s.surfW[g]];
}

// The local constraints' rows (steady solves only): Σ w c_i − reference at each node.
function applyLocalConstraints(s) {
  const { n, M, c, th, nSurf, res, dA: d } = s;
  for (const lc of s.localConstraints) {
    lc.nodes.forEach((g, j) => {
      const b = s.blockOfNode[g];
      if (s.loc[b * M + lc.row] < 0) return;
      s._replaceRow(b, lc.row); // (a port reaction's ∂/∂V too, where one touched this slot)
      let t = 0;
      for (const { i, q, k, w } of lc.parts) {
        if (i !== undefined) {
          t += w * c[g * n + i];
          s._dc(g, i, d);
          for (let r = 0; r < M; r++) if (d[r] !== 0) s._j(b, lc.row, b, r, w * d[r]);
          continue;
        }
        // A coverage, with ∂θ_q/∂η_x = θ_q (δ_qx − θ_x) over its surface's species (where it is).
        if (s.surfAt[g] !== k) continue;
        const tq = th[g * nSurf + q], ns = s.model.surfaces[k].species.length, ws = w * sites(s, g);
        t += ws * tq;
        for (let x = 0; x < ns; x++) s._j(b, lc.row, b, 1 + n + x, ws * tq * ((q === x ? 1 : 0) - th[g * nSurf + x]));
      }
      res[s.rix[b * M + lc.row]] = t - lc.reference[j];
    });
  }
}

// Every constraint's row (steady solves only): a spectator's conserved amount, Σ v c_i = amount
// over its stretch, and the other constraints and affinities in the rows findLaws gave them,
// balance rows and circuit rows alike. A stretch's balance row (its first node's) is in steady
// state the negative sum of the others. That row is dense, so it's kept aside: the factorised matrix gets a pin
// (identity row) there instead, and solveBordered restores the constraint.
function applyConstraints(s) {
  const { M, res } = s;
  for (const cs of s.constraints) {
    const { amount, reference, dV } = lawAmount(s, cs);
    cs.res = amount - (cs.terminals.length > 0 ? (s.lawRefs.get(cs.key) ?? amount) : reference);
    if (cs.circuit >= 0) {
      // In place of the terminal's circuit row (its I = 0).
      const C = s.termC[cs.circuit];
      C.fill(0);
      for (let j = 0; j < cs.len; j++) C[cs.idx[j]] += cs.w[j];
      for (const { t, lambda } of cs.contacts) {
        // A contact's charge, D_in, as Solver._contact had it, with its derivatives.
        const q = s.contactQ[t], f = lambda / FARADAY;
        cs.res += f * q.D;
        for (let j = 0; j < C.length; j++) C[j] += f * q.row[j];
        dV[t] += f * q.dV;
      }
      s.termRes[cs.circuit] = cs.res;
      s.termDI[cs.circuit] = dV[cs.circuit];
      dV.forEach((x, k) => k !== cs.circuit && (s.termDX[cs.circuit][k] = x));
      continue;
    }
    const b0 = Math.floor(cs.row / M), r0 = cs.row % M;
    s._replaceRow(b0, r0);
    s._j(b0, r0, b0, r0, 1);
    res[s.rix[cs.row]] = 0;
  }
}

// A conserved combination's amount as the state holds it (and its reference, closed), with its
// derivatives in cs.idx/w (and dV, by the voltage of a terminal whose electrode's charge it
// counts).
function lawAmount(s, cs) {
  const { n, M, c, dA: d } = s;
  const vol = s.model.grid.vol;
  let amount = 0, reference = 0, len = 0, gross = 0;
  for (const { stretch, w } of cs.parts) {
    const st = s.stretches[stretch], i = st.species;
    if (cs.terminals.length === 0) reference += w * s.referenceAmounts[stretch];
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
      amount += w * vol[g] * c[g * n + i];
      gross += Math.abs(w * vol[g] * c[g * n + i]);
      if (s.nodeConductor[g] >= 0) {
        // A metal's carriers: the sheet at a charged face, vol·c = ±D_f·A/(zF), held by the
        // face's displacement.
        const f = s.sheetFace[g];
        if (f >= 0 && s.loc[s.blockOfFace[f] * M] >= 0) {
          cs.idx[len] = s.rix[s.blockOfFace[f] * M];
          cs.w[len++] = (w * s.sheetSign[g] * s.model.grid.area[g]) / (s.z[i] * FARADAY);
        }
        continue;
      }
      s._dc(g, i, d);
      const b = s.blockOfNode[g];
      for (let r = 0; r < M; r++) {
        if (d[r] === 0 || s.loc[b * M + r] < 0) continue;
        cs.idx[len] = s.rix[b * M + r];
        cs.w[len++] = w * vol[g] * d[r];
      }
    }
  }
  // What the faces' diffuse layers hold of it (src/diffuse.js): A_f Γ_i at their edge nodes, Γ a
  // function of the layer's drop and the edge node's c (through its ∂c_j/∂(φ̂, η): _diffuseChain).
  for (const { stretch, w } of cs.parts) {
    const st = s.stretches[stretch], i = st.species;
    s.diffList.forEach((dl, k) => {
      const g = dl.g;
      if (g < st.nodes[0] || g > st.nodes[1] || !dl.holds[i]) return;
      const f = w * s.model.grid.area[g], b = dl.b, D = dl.dl, dc = s._diffuseChain(k);
      amount += f * s.gam[k * n + i];
      gross += Math.abs(f * s.gam[k * n + i]);
      cs.idx[len] = s.rix[dl.o];
      cs.w[len++] = f * D.dgds[i];
      for (let r = 0; r <= n; r++) {
        if (s.loc[b * M + r] < 0) continue;
        let v = 0;
        for (let j = 0; j < n; j++) v += D.dgdc[i * n + j] * dc[j * M + r];
        if (v === 0) continue;
        cs.idx[len] = s.rix[b * M + r];
        cs.w[len++] = f * v;
      }
    });
  }
  // What the electrodes' surfaces hold of it: Γ Σ v·a·θ, with ∂θ_s/∂η_x = θ_s (δ_sx − θ_x).
  for (const { col, w } of cs.surface) {
    const [k, q] = s.surfCols[col], sf = s.model.surfaces[k], G = sf.species[0].capacity, ns = sf.species.length;
    if (cs.terminals.length === 0) reference += w * s.surfaceRef[col];
    sf.nodes.forEach((g, j) => {
      const f = w * G * vol[g] * s.surfArea[k][j], b = s.blockOfNode[g], tq = s.th[g * s.nSurf + q];
      amount += f * tq;
      gross += Math.abs(f * tq);
      for (let x = 0; x < ns; x++) {
        cs.idx[len] = s.rix[b * M + 1 + n + x];
        cs.w[len++] = f * tq * ((q === x ? 1 : 0) - s.th[g * s.nSurf + x]);
      }
    });
  }
  // What an electrode holds, its charge Q/(zF) on the port's capacitance (none without one),
  // Q = Σ vol·a·C (V − zeroCharge − φ) over its nodes (as Solver._portCapacitance has it).
  const dV = cs.dV ?? (cs.dV = new Float64Array(s.terms.length));
  dV.fill(0);
  for (const { k, w } of cs.carriers) {
    const port = s.model.ports[k - 2];
    if (!port.capacitance) continue;
    const f = w / (carrierZ(s, k) * FARADAY), { C: Cs, zeroCharge } = port.capacitance, V = s.termV[k];
    const sets = (side) => ['pinned', 'bulk'].includes(s.model.contacts[side].phi.type);
    port.nodes.forEach((g, j) => {
      if ((g === 0 && sets('left')) || (g === s.nNodes - 1 && sets('right')) || s.phiUndefined[g]) return;
      const b = s.blockOfNode[g], cap = vol[g] * s.portArea[k - 2][j] * Cs;
      amount += f * cap * (V - zeroCharge - s.VT * (s.u[b * M] + s.uLo[b * M]));
      gross += Math.abs(f * cap * (V - zeroCharge - s.VT * (s.u[b * M] + s.uLo[b * M])));
      dV[k] += f * cap;
      cs.idx[len] = s.rix[b * M];
      cs.w[len++] = -f * cap * s.VT;
    });
  }
  // (A typed array drops writes past its end without a word: the law's derivative would lose them.)
  if (len > cs.idx.length) throw new Error(`internal: a law's derivative needs ${len} entries, its workspace holds ${cs.idx.length}`);
  cs.len = len;
  return { amount, reference, dV, gross }; // (dV by terminal)
}

// The same amount at the step's start, as the storage terms have it (BDF2's c* included).
function lawStart(s, cs) {
  const { n, nSurf } = s, vol = s.model.grid.vol;
  let a = 0;
  for (const { stretch, w } of cs.parts) {
    const st = s.stretches[stretch], i = st.species;
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) a += w * vol[g] * s.cOld[g * n + i];
    s.diffList.forEach((dl, k) => dl.g >= st.nodes[0] && dl.g <= st.nodes[1] && (a += w * s.model.grid.area[dl.g] * s.gamOld[k * n + i]));
  }
  for (const { col, w } of cs.surface) {
    const [k, q] = s.surfCols[col], sf = s.model.surfaces[k], G = sf.species[0].capacity;
    sf.nodes.forEach((g, j) => (a += w * G * vol[g] * s.surfArea[k][j] * s.thOld[g * nSurf + q]));
  }
  for (const { k, w } of cs.carriers) if (s.model.ports[k - 2].capacitance) a += (w * s.portQStart[k - 2]) / (carrierZ(s, k) * FARADAY);
  return a;
}

// Time steps: a terminal driven by a current that alone feeds a conserved combination (the
// only one of its `terminals` passing a current; with more, each step reads the flux) can read
// its current as the combination's change, F (S − S*)/(λ dt), which the step's balance rows sum
// to exactly, rather than as the flux through it. Both are exact; they differ in round-off and
// in what Newton can see. Through a metal, the flux is G Δη, its two η all but equal: a floating
// gate's charging current fell below that at steps of ~0.03 s, and the step came out singular.
// Over long steps the change is the better kept (a closed host at open circuit, read by its
// flux, drifted from its charge by 3 mV over 5e11 s, and its steps failed by 4e12 s). Over
// short ones it's the amount's round-off over dt (a battery's cathode holds 0.04 mol/m² of
// electrons: 1e-6 A/m² at steps of 1e-6 s, enough to reject steps at tol 1e-8), where the flux,
// its levels carried in two words (u + uLo), is far cleaner than ε G. So a step reads the change
// where that's ε·F·S/(λ dt) below 1e-6 of the flux's ε G, as its start has them, keeps to it
// throughout, and where a step read by the flux fails, it's tried once more read by the change
// (see step).
function circuitLaws(s, dt) {
  for (const cs of [...s.constraints, ...s.timeLaws]) {
    const k = cs.circuit;
    if (k < 0 || !s.floating.includes(k) || cs.contacts?.length > 0) continue; // (a contact's charge: its own rows)
    // (Fed by another current-driven terminal passing a current, the change is theirs together,
    // not k's alone: read by the flux instead. Read by the change, K⁺ driven in at a contact and
    // out at a port at equal currents gained I·dt/F a step, the contact passing twice its
    // current, silently; Fable round 3's review of src/lawRows.js. Others at I = 0 add nothing.)
    if (cs.terminals?.some((t) => t !== k && sourceAt(s.terms[t].drive.src, s.sourceTime, s.sourceBefore) !== 0)) continue;
    const { amount, dV, gross } = lawAmount(s, cs), f = FARADAY / (cs.lambda * dt), C = s.termC[k];
    if (!s.lawForm.has(cs)) {
      let G = 0;
      for (let j = 0; j < C.length; j++) G = Math.max(G, Math.abs(C[j]));
      s.lawForm.set(cs, Math.abs(f) * gross < 1e-6 * G);
    }
    if (!s.lawForm.get(cs)) continue;
    C.fill(0);
    for (let j = 0; j < cs.len; j++) C[cs.idx[j]] += f * cs.w[j];
    const d = s.terms[k].drive, src = sourceAt(d.src, s.sourceTime, s.sourceBefore);
    s.termRes[k] = f * (amount - lawStart(s, cs)) - (d.kind === 'I' ? src : (src - s.termV[k]) / d.R);
    s.termDI[k] = f * dV[k] + (d.kind === 'I' ? 0 : 1 / d.R);
    dV.forEach((x, j) => j !== k && (s.termDX[k][j] = f * x));
  }
}

// Solve the full Newton system, J δ = rhs, with the extra unknowns and rows that don't fit the
// block-tridiagonal matrix T: the floating terminals' voltages (with their circuit rows), and
// the rows pinned in T: the constraints and affinities (see applyConstraints) and the
// islands' summed balances (see applyIslands). By low-rank
// updates of T (Woodbury): δ = y + Σ_q Q_q μ_q − Σ_k X_k δV_k with y = T⁻¹ rhs, Q_q = T⁻¹ e_q
// (the response to a unit pin) and X_k = T⁻¹ B_k (to a unit change of V_k), then a small dense
// system for the μ (pins) and δV (terminals):
//   pin q:      Σ (W_q·Q_q′) μ_q′ − Σ (W_q·X_k) δV_k = res_q − W_q·y
//   terminal k: Σ (C_k·Q_q) μ_q + Σ (δ_kk′ ∂I_k/∂V_k − C_k·X_k′) δV_k′ = res_k − C_k·y
export function solveBordered(s, rhs, delta, deltaV, pins) {
  const N = s.sys.size, fl = s.floating, P = pins.length, K = fl.length, S = P + K;
  s._solveLinear(rhs, delta);
  if (S === 0) return;
  const e = s.dWork ?? (s.dWork = new Float64Array(N + 1));
  const cols = s.borderCols ?? (s.borderCols = []);
  for (let a = 0; a < S; a++) {
    if (!cols[a]) cols[a] = new Float64Array(N + 1);
    if (a < P) {
      // The response to a unit pin, or where the pin moves a flat stretch, to its shift.
      const col = pins[a].col ?? [pins[a].row];
      for (const o of col) e[s.rix[o]] = 1;
      s._solveLinear(e, cols[a]);
      for (const o of col) e[s.rix[o]] = 0;
    } else s._solveLinear(s.termB[fl[a - P]], cols[a]);
  }
  // Row a of the small system, as a dot product with a compact vector (sparse for pins).
  const dot = (a, v) => {
    let t = 0;
    if (a < P) {
      const p = pins[a];
      for (let j = 0; j < p.len; j++) t += p.w[j] * v[p.idx[j]];
    } else {
      const C = s.termC[fl[a - P]];
      for (let j = 0; j < N; j++) t += C[j] * v[j];
    }
    return t;
  };
  const A = Array.from({ length: S }, () => new Float64Array(S + 1));
  for (let a = 0; a < S; a++) {
    for (let b = 0; b < S; b++) A[a][b] = (b < P ? 1 : -1) * dot(a, cols[b]);
    if (a < P && pins[a].termSlope) for (let b = P; b < S; b++) A[a][b] += pins[a].termSlope.get(fl[b - P]) ?? 0; // (a metal island behind a resistance)
    if (a >= P) {
      A[a][a] += s.termDI[fl[a - P]];
      for (let b = P; b < S; b++) if (b !== a) A[a][b] += s.termDX[fl[a - P]][fl[b - P]];
    }
    A[a][S] = (a < P ? pins[a].res : s.termRes[fl[a - P]]) - dot(a, delta);
    let mx = 0;
    for (let b = 0; b < S; b++) mx = Math.max(mx, Math.abs(A[a][b]));
    if (mx > 0) for (let b = 0; b <= S; b++) A[a][b] /= mx;
  }
  // Small dense solve with partial pivoting.
  for (let col = 0; col < S; col++) {
    let p = col;
    for (let r = col + 1; r < S; r++) if (Math.abs(A[r][col]) > Math.abs(A[p][col])) p = r;
    [A[col], A[p]] = [A[p], A[col]];
    if (A[col][col] === 0) throw new SolverError('the bordered system is singular: a conserved amount or a floating terminal is not determined');
    for (let r = col + 1; r < S; r++) {
      const f = A[r][col] / A[col][col];
      for (let k = col; k <= S; k++) A[r][k] -= f * A[col][k];
    }
  }
  const x = new Float64Array(S);
  for (let r = S - 1; r >= 0; r--) {
    let v = A[r][S];
    for (let k = r + 1; k < S; k++) v -= A[r][k] * x[k];
    x[r] = v / A[r][r];
  }
  for (let a = 0; a < S; a++) {
    const sa = a < P ? x[a] : -x[a], col = cols[a];
    for (let j = 0; j < N; j++) delta[j] += sa * col[j];
    if (a >= P) deltaV[a - P] = x[a];
  }
}

// A flat stretch's level: its contact's, through the terminal k that holds it, or (reaching
// none) its first node's, which its amount fixes (that row is the amount's, a pin whose
// response moves the whole stretch: see applyConstraints).
function flatLevel(s, st) {
  const { M, u, uLo, z, VT, model } = s, i = st.species;
  if (st.leftOpen || st.rightOpen) {
    const side = st.leftOpen ? 'left' : 'right', link = model.contacts[side].species[i], k = side === 'left' ? 0 : 1;
    return { k, hi: z[i] === 0 ? link.mu / model.RT : (z[i] * (s.termV[k] + link.offset)) / VT, lo: 0 };
  }
  const o = s.blockOfNode[st.nodes[0]] * M + 1 + i;
  return { k: -1, hi: u[o], lo: uLo[o] };
}

// Each flat stretch's nodes set to its first node's level exactly (where Newton's damped
// updates have taken it, on the way to its contact's), so that it carries exactly no flux (the
// expm1 form of a zero difference). Flat only to round-off, a level moved by its contact had
// read as a current: a MOS capacitor's holes, their contact jumped to −6 V, as
// 3e-4 A/m² through a bulk that conducts ~1e10 S/m² (the true current was 3.5e-7).
export function snapFlat(s) {
  const { M, u, uLo } = s;
  for (const st of s.flatStretches) {
    const f = s.blockOfNode[st.nodes[0]] * M + 1 + st.species, hi = u[f], lo = uLo[f];
    for (let g = st.nodes[0] + 1; g <= st.nodes[1]; g++) {
      const o = s.blockOfNode[g] * M + 1 + st.species;
      if (s.loc[o] < 0) continue;
      u[o] = hi;
      uLo[o] = lo;
    }
  }
}

// Steady solves: each flat stretch's level held at its contact's throughout (see
// flatStretches), in place of its balances, as an equilibrium link holds it at the end node.
function flatRows(s) {
  const { M, u, uLo, res, z, VT, model } = s, R = s.rix;
  for (const st of s.flatStretches) {
    const i = st.species, open = st.leftOpen || st.rightOpen;
    const { k, hi, lo } = flatLevel(s, st);
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
      const b = s.blockOfNode[g], o = b * M + 1 + i;
      if (s.loc[o] < 0) continue;
      s._replaceRow(b, 1 + i); // (the row's old terminal terms too: the end node's own link)
      s._j(b, 1 + i, b, 1 + i, 1);
      res[R[o]] = u[o] - hi + (uLo[o] - lo); // (the big terms first: hi + lo rounds to hi)
      if (open && z[i] !== 0) s.termB[k][R[o]] += -z[i] / VT;
    }
    // A face inside the stretch whose link holds the level continuous carries the species' flux
    // as an unknown that only the edge balances, now replaced, determined: it's zero.
    for (let f = st.regions[0]; f < st.regions[1]; f++) {
      if (model.interfaces[f].links[i].type !== 'equilibrium') continue;
      const b = s.blockOfFace[f], o = b * M + 1 + i;
      if (s.loc[o] < 0) continue;
      s._replaceRow(b, 1 + i);
      s._j(b, 1 + i, b, 1 + i, 1);
      res[R[o]] = u[o] + uLo[o];
    }
  }
}

// Steady solves: each reaction that alone changes a combination (see src/laws.js) held at its
// equilibrium, a = fixedA − Σ ν η = 0, in place of one of that combination's balance rows (a pin,
// like an amount's).
function applyAffinities(s) {
  const { M, u, uLo, res } = s, grid = s.model.grid;
  for (const af of s.affinities) {
    const b0 = Math.floor(af.row / M), r0 = af.row % M;
    s._replaceRow(b0, r0);
    s._j(b0, r0, b0, r0, 1);
    res[s.rix[af.row]] = 0;
    const gL = grid.regionEnd[af.f], gR = grid.regionStart[af.f + 1];
    let hi = af.rx.fixedA, lo = 0, len = 0;
    const sSide = s.model.interfaces[af.f].surfaceSide;
    const terms = [...af.rx.part.map((p) => [s.blockOfNode[p.side ? gR : gL] * M + 1 + p.i, p.nu]), ...af.rx.surf.map((p) => [s.blockOfNode[sSide ? gR : gL] * M + 1 + s.n + p.s, p.nu])];
    for (const [o, nu] of terms) {
      hi -= nu * u[o];
      lo -= nu * uLo[o];
      if (s.loc[o] < 0) continue;
      af.idx[len] = s.rix[o];
      af.w[len++] = -nu;
    }
    af.len = len;
    af.res = hi + lo;
  }
}
