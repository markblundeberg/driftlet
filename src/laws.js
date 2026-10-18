// The steady system's conservation laws: which amounts a steady solve keeps, and which of its
// equations each takes the place of.
//
// A steady state's balance rows are dependent wherever an amount is conserved: the weighted sum
// of a conserved combination's balances is zero for every state, so one of them says nothing,
// and the amount itself (which the state has, and a steady solve must not lose) is what goes in
// its place. Likewise a reaction whose rate a steady state sets to zero is better said by its
// affinity, a species nothing carries by its amount node by node, and a piece of a stretch fed
// only through its faces by its summed balance. This module finds them all, from the device's
// structure alone; src/lawRows.js evaluates them on the state and writes them into the Newton
// system.
//
// CONTRACT
//
// Input: a structure `s` (Solver._structure()), read and never written:
//   model                the normalised device (regions, materials, interfaces, ports, contacts,
//                        reactions, species, grid; ports[k] is terminal k + 2)
//   n, M, z, nNodes      species, slots per block (1 + n + the most extra slots), charges, nodes
//   stretches            the connected runs of regions where a species is present, each
//                        { species, regions: [r0, r1], nodes: [g0, g1], leftOpen, rightOpen,
//                        linked (ports holding or exchanging it), ports, reactive, spectator,
//                        contactFed, mobile, connected }
//   stretchOf            stretch index by region·n + species, −1 where absent
//   flatStretches        stretches a steady solve holds at a flat level (Solver._findStretches)
//   rough                each stretch's starting amount (Infinity where a level sets it)
//   surfCols, surfColOf  electrode surface species as columns: surfCols[j] = [port, slot];
//                        surfColOf[port] its first column (−1 for none)
//   terms, floating      terminals (0 left, 1 right, 2 + k port k) with their drives as they are
//                        now; the terminals whose voltage is an unknown
//   blockOfNode, blockOfFace, rxsIn, nodeMaterial, phiUndefined
//                        the block layout (row o = block·M + slot), bulk reactions by material,
//                        each node's material, nodes with no φ
// (The small predicates below, moves … chargeOnly, read only the fields they destructure; the
// solver also calls them on itself.)
//
// Column spaces. A law's weights w are over C columns: stretches [0, S), electrode surfaces'
// species [S, X) (column k is surface column k − S), electrodes' charges [X, C) (column k is
// terminal carriers[k − X], counting its charge as Q/(zF), z = carrierZ). The rows a law can
// take, v, are over S + T: stretch k's balance row (at its first node) [0, S), then terminal t's
// circuit row, S + t.
//
// Output of findLaws(s): { constraints, affinities, localConstraints, islands, timeLaws,
// contactLaw, carriers, conserved, unkept }:
//   constraints[]  a kept amount, Σ w·n_stretch + Σ w·(surface's Γ a θ) + Σ w·Q_k/(z_k F) (ports
//                  with a capacitance) + Σ λ·D_t/F (contacts' charges), from
//                  { parts: [{stretch, w}], surface: [{col, w}], carriers: [{k, w}],
//                  contacts: [{t, lambda}] }. It replaces either one balance row, `row` (block·M
//                  + slot, circuit −1), or terminal `circuit`'s row, I = 0, where the amount
//                  changes by λ I/F (`lambda`). `terminals` are those driven by a current that
//                  feed it: with none it keeps its reference (the amounts at construction), with
//                  some what the state holds when a solve begins (captured under `key`). `col`,
//                  for a stretch held flat: the rows its pin moves together. idx, w, len, res, dV:
//                  workspace for the solver.
//   affinities[]   { rx, f, row }: face reaction rx at interface f, its affinity Σ ν η = 0 in place
//                  of balance row `row`; idx, w, len, res: workspace.
//   localConstraints[]  { parts: [{i, w}] or [{q, k, w}], row, nodes, reference }: at each node
//                  g of `nodes`, Σ w c_i + Σ w Γ a θ_q (a surface's only where port k is), as it
//                  was, in place of slot `row` of g's block. Only of still stretches (D = 0 in
//                  one region: nothing carries them between nodes) that no terminal feeds, and
//                  surfaces. One law can have several entries, a row per set of nodes.
//   islands[]      { stretch, row, flux, weight, outside }: a piece of a stretch between faces fed
//                  by nothing else, its balance summed over its nodes (Σ weight·u over the `flux`
//                  rows, plus the `outside` rows' balances) in place of balance row `row`.
//   timeLaws[]     time steps only: as constraints (no row, key or contacts), each for a terminal
//                  behind a resistance that alone feeds it.
//   contactLaw[t]  whether some constraint keeps contact t's charge.
//   carriers       the terminals with an electrode's-charge column, in column order.
//   conserved      the stretches some constraint or local constraint keeps (indices).
//   unkept         the conserved combinations (basis vectors) that no row here can keep. With
//                  any, or any stretch neither reached by a contact nor conserved, the solver
//                  doesn't solve the steady system directly but by time steps, which keep every
//                  amount themselves (Solver._directSteady).
//
// Invariants (each says what a reviewer should be able to check here):
//   1. conservedMoieties' basis spans the null space of its rows: every reaction's stoichiometry,
//      and what each terminal feeds (freely, held at a voltage; in proportion to charge, driven
//      by a current, with `drives`). Every basis vector w has w·row = 0 for every row.
//   2. No equation is replaced twice and none is lost: the rows the constraints and affinities
//      take are chosen by elimination over v (one map from w, rowsOf), each law's v reduced by
//      those already chosen, so their v are independent. The local constraints' rows are chosen
//      the same way node by node, over the slots there. Islands take rows no other law has.
//   3. A balance row taken is a moving stretch's (a floating metal's own amount takes its
//      carriers'). An affinity's is never a metal's carriers', nor a stretch's that a terminal
//      feeds, nor a circuit row.
//   4. A charge is kept only in a circuit row (its pin there carries the voltage it depends on):
//      an electrode's, at weight w/z on its carrier column, and a contact's (D_in, no column), at
//      λ = w/z on what it passes. The first law through each takes a circuit row (the first
//      terminal's in order whose v is left), and the laws after it have its multiple taken off.
//      One with a charge left and no circuit row is unkept.
//   5. Every basis vector of conservedMoieties(s, { drives: true }) is accounted for: a
//      spectator's own amount, a constraint, local constraints node by node, unkept, or dropped
//      for one of the reasons marked DROP in findLaws. Each DROP states why nothing is lost.
//   6. The pivots and the affinities' rows don't depend on the order species, reactions or ports
//      are listed in: they're chosen by the reactions' own terms and by name. A constraint's row
//      can (the total iron in Fe²⁺'s row or Fe³⁺'s): the same law in another of its rows, with
//      the same solution.
//   7. Nothing in `s` is written.
// Tests: test/laws.test.js (one device per law once dropped); tools/lawsweep.js (810 devices,
// solve() against a settled advance()).

// Whether a stretch's species moves through it (mobile in every region, or a metal's carrier),
// or stays still: D = 0 in its one region, so nothing carries it between nodes (trap states,
// conserved node by node). A stretch over two regions is joined by a face that passes it, which
// carries it between the face's nodes: still in neither sense. (Kept node by node, S⁺ linked
// across a face stayed at 1 and 0.1 on its two nodes, where the transient evens them to 0.4.)
export function moves({ model }, st) {
  const { regions, materials } = model;
  for (let q = st.regions[0]; q <= st.regions[1]; q++) {
    const mat = materials[regions[q].material];
    if (!(mat.conductor || mat.D[st.species] > 0)) return false;
  }
  return true;
}

export function still({ model }, st) {
  const { regions, materials } = model;
  if (st.regions[0] !== st.regions[1]) return false;
  for (let q = st.regions[0]; q <= st.regions[1]; q++) {
    const mat = materials[regions[q].material];
    if (mat.conductor || mat.D[st.species] > 0) return false;
  }
  return true;
}

// The charge number an electrode's carrier column counts its charge in, Q/(zF): its terminal
// species', or 1 for a port with none (a membrane's capacitance: Q/F).
export function carrierZ({ model, z }, k) {
  const t = model.ports[k - 2].terminal;
  return t === null || t === undefined ? 1 : z[t];
}

// Whether contact k (0, 1) holds a charge its current counts, D_in on its φ law (capacitive,
// pinned or bulk, where φ is defined at its end node).
export function contactCharges({ model, nNodes, phiUndefined }, k) {
  const g = k === 0 ? 0 : nNodes - 1;
  return ['capacitive', 'pinned', 'bulk'].includes(model.contacts[k === 0 ? 'left' : 'right'].phi.type) && !phiUndefined[g];
}

// Whether terminal k (a contact, or port k − 2) feeds stretch st: links its species there.
export function feeds(k, st) {
  return k === 0 ? st.leftOpen : k === 1 ? st.rightOpen : st.linked.includes(k - 2);
}

// Whether a terminal is driven by a current held at zero throughout (no resistance).
export function noCurrent(t) {
  return t.drive.kind === 'I' && !(t.drive.R > 0) && t.drive.src.value === 0;
}

// Whether terminal k passes current only by charging a capacitance: no charged species through
// it (nor, a port, reactions), and a capacitance to charge (a gate's φ law, a port's own).
export function chargeOnly({ model, terms }, k) {
  const t = terms[k], { species } = model;
  const ions = (links) => links.some((l, i) => l.type !== 'blocked' && species[i].z !== 0);
  if (t.kind === 'port') {
    const p = model.ports[t.index];
    return p.reactions.length === 0 && !ions(p.species) && !!p.capacitance;
  }
  const c = model.contacts[t.side];
  return !ions(c.species) && (c.phi.type === 'capacitive' || c.phi.type === 'pinned');
}

// A basis of the conserved combinations of stretch amounts: the null space of the
// stoichiometry (a row per reaction, over the stretches it touches) together with what the
// terminals feed. A terminal held at a voltage (or behind a resistance, which holds it at
// steady state) feeds each stretch it links freely: a unit row. With `drives`, one driven by a
// current feeds its charged species only as that current allows, in proportion to their
// charges, so a combination weighted by z there is conserved while it passes none: a battery's
// cathode, its electrons behind an open circuit, keeps its charge. Exact integer data, so plain
// elimination will do. Returns the basis (each vector over the columns, w.own its free column),
// with basis.carriers (the terminals with a charge column) and basis.skipped.
export function conservedMoieties(s, { drives = false, alsoCurrent = [], skip = null } = {}) {
  const { model, n, z, stretches, stretchOf, terms, surfCols, surfColOf, rough } = s, S = stretches.length;
  // With `drives`, an electrode port driven by a current whose reactions take its carrier from
  // the electrode has a column of its own, last: what the electrode holds, its charge Q/(zF),
  // which the current feeds and the reactions draw on (a port at no current with M⁺ + e⁻ = M(s)
  // keeps M⁺ + Q/F, or with no capacitance, Q = 0, the M⁺ alone).
  const byCurrent = (t, k) => t.drive.kind === 'I' || alsoCurrent.includes(k); // (with `alsoCurrent`, terminals counted as such)
  // (Or a port with a capacitance that passes ions, a membrane's: what it holds is Q/F. One that
  // only charges keeps its charge alone, see lawRows.js, captureLaws.)
  const carriers = drives ? terms.flatMap((t, k) => (byCurrent(t, k) && k >= 2 && (model.ports[k - 2].reactions.length > 0 || (model.ports[k - 2].capacitance && !chargeOnly(s, k))) ? [k] : [])) : [];
  const X = S + surfCols.length, C = X + carriers.length;
  // (columns: stretches [0, S), surfaces' species [S, X), electrodes' charges [X, C), as in findLaws)
  const isStretch = (k) => k < S, isSurface = (k) => k >= S && k < X;
  // (With `skip`, a reaction left out: the combinations that it alone changes come out too, and
  // its own rows, over the columns, are handed back as `basis.skipped`.)
  const rows = [], skipped = [];
  model.interfaces.forEach((itf, f) => {
    for (const rx of itf.reactions) {
      const row = new Float64Array(C);
      for (const p of rx.part) row[stretchOf[(f + p.side) * n + p.i]] += p.nu;
      (rx === skip ? skipped : rows).push(row);
    }
  });
  // An electrode surface's species are columns too, after the stretches (S + j for surface
  // column j), so a combination that passes through them (A⁺ + e⁻ = S, S + e⁻ = B⁻) is conserved
  // with what the surface holds.
  model.ports.forEach((port, k) => {
    for (const rx of port.reactions) {
      const row = new Float64Array(C); // (the electrode's carrier is outside the device)
      if (rx === skip) skipped.push(row);
      for (const p of rx.part) {
        if (p.side === 0) row[stretchOf[port.region * n + p.i]] += p.nu;
        else if (p.side === 2) row[S + surfColOf[k] + p.s] += p.nu;
        else if (carriers.includes(k + 2)) row[X + carriers.indexOf(k + 2)] += p.nu;
      }
      if (rx !== skip) rows.push(row);
    }
  });
  model.regions.forEach((reg, q) => {
    for (const rx of model.reactions) {
      if (!(rx.kf[reg.material] > 0)) continue;
      const row = new Float64Array(C);
      for (const p of rx.reactants) row[stretchOf[q * n + p.i]] -= p.nu;
      for (const p of rx.products) row[stretchOf[q * n + p.i]] += p.nu;
      (rx === skip ? skipped : rows).push(row);
    }
  });
  const unit = (k) => {
    const row = new Float64Array(C);
    row[k] = 1;
    rows.push(row);
  };
  terms.forEach((t, k) => {
    const fed = stretches.flatMap((st, j) => (feeds(k, st) ? [j] : []));
    if (!(drives && byCurrent(t, k))) return fed.forEach(unit);
    const charged = fed.filter((j) => z[stretches[j].species] !== 0).map((j) => [j, z[stretches[j].species]]);
    if (carriers.includes(k)) charged.push([X + carriers.indexOf(k), carrierZ(s, k)]);
    fed.filter((j) => z[stretches[j].species] === 0).forEach(unit);
    for (let a = 1; a < charged.length; a++) {
      const row = new Float64Array(C);
      row[charged[0][0]] = 1 / charged[0][1];
      row[charged[a][0]] = -1 / charged[a][1];
      rows.push(row);
    }
  });
  // Reduced row echelon form; the free columns give the null space. Immobile stretches' columns
  // go first, and electrode surfaces' (as local), so that a combination of them alone (trap
  // states, X⁰ + X⁻; adsorbates, A(ads) + B(ads)) comes out on its own, to be kept node by node,
  // apart from any that moves (e⁻ − X⁰, the electrons the traps hold). (Pivoted after a closed
  // Cl⁻ in A(ads) + Cl⁻ = B(ads) + e⁻, the surface's own came out mixed with it, and was lost.)
  // Then the stretches by their starting amounts, least first, and by name: a trace species
  // pivoted after an abundant one comes out added to it (A + B ⇌ C with B at 1 M and A at 1e-14
  // gave B + C and A − B, where A + C was lost to round-off), and the basis followed the listing.
  const isStill = (k) => ((isStretch(k) && still(s, stretches[k])) || isSurface(k) ? 1 : 0);
  const name = (k) => model.species[stretches[k].species].name;
  const order = [...Array(C).keys()].sort((a, b) => {
    if (isStill(a) !== isStill(b)) return isStill(b) - isStill(a);
    if (!isStretch(a) || !isStretch(b)) return a - b;
    if (rough[a] !== rough[b]) return rough[a] - rough[b];
    const na = name(a), nb = name(b);
    return na < nb ? -1 : na > nb ? 1 : a - b;
  });
  const pivots = [];
  let r = 0;
  for (const col of order) {
    if (r >= rows.length) break;
    let best = r;
    for (let i = r + 1; i < rows.length; i++) if (Math.abs(rows[i][col]) > Math.abs(rows[best][col])) best = i;
    if (Math.abs(rows[best][col]) < 1e-9) continue;
    [rows[r], rows[best]] = [rows[best], rows[r]];
    const pr = rows[r], pv = pr[col];
    for (let j = 0; j < C; j++) pr[j] /= pv;
    for (let i = 0; i < rows.length; i++) {
      if (i === r || rows[i][col] === 0) continue;
      const fct = rows[i][col];
      for (let j = 0; j < C; j++) rows[i][j] -= fct * pr[j];
    }
    pivots.push(col);
    r++;
  }
  const isPivot = new Uint8Array(C);
  for (const col of pivots) isPivot[col] = 1;
  const basis = [];
  for (let free = 0; free < C; free++) {
    if (isPivot[free]) continue;
    const w = new Float64Array(C);
    w[free] = 1;
    pivots.forEach((col, i) => {
      const v = -rows[i][free];
      w[col] = Math.abs(v) < 1e-9 ? 0 : v;
    });
    w.own = free; // (its weight is 1 there, and no other vector's is: the row it replaces)
    basis.push(w);
  }
  basis.skipped = skipped;
  basis.carriers = carriers;
  return basis;
}

// What a steady solve keeps (under the terminals' present drives, so again after a change of
// drive): each conserved amount, solved directly in place of one (redundant) balance row, and
// the islands' summed balances. See the contract above.
export function findLaws(s) {
  const { model, n, M, z, nNodes, stretches, stretchOf, flatStretches, surfCols, terms, floating, blockOfNode, blockOfFace, rxsIn, nodeMaterial } = s;
  const { regions, materials, grid } = model;
  const laws = conservedMoieties(s, { drives: true }), carriersOf = laws.carriers;
  const conserved = new Set();
  // Conserved amounts solved directly, each in place of one (redundant) balance row: a
  // spectator's own amount (a floating metal's: the charge on its faces), each conserved
  // combination of reacting stretches (the total iron of Fe³⁺, Fe²⁺ and FeCl²⁺, say), whose
  // weighted balances sum to zero at steady state, and what a terminal driven at no current
  // keeps behind it. A combination may take in still stretches (the electrons trap states hold),
  // their amounts summed like any other; one of still stretches alone conserves node by node
  // instead (below), and one of a stretch still in part only is left to giant time steps.
  const constraints = [];
  const conductor = (st) => materials[regions[st.regions[0]].material].conductor;
  // The terminals driven by a current that feed a stretch: a combination through one keeps what
  // the state holds when a solve begins (see lawRows.js, captureLaws), and has no steady state at
  // any other current.
  const driven = (st) => terms.flatMap((t, k) => (t.drive.kind === 'I' && feeds(k, st) ? [k] : []));
  const add = (parts, rowStretch, surface = [], circuit = -1, carriers = []) => {
    const st = stretches[rowStretch];
    // (the terminal's weight in it, w/z on what that feeds: the amount changes by λ I/F)
    const own = carriers.find((p) => p.k === circuit);
    const lambda = circuit < 0 ? 0 : ((own && own.w / carrierZ(s, circuit)) || parts.map((p) => (driven(stretches[p.stretch]).includes(circuit) ? p.w / z[stretches[p.stretch].species] : 0)).find((x) => x !== 0));
    // In a circuit row, what passes through that terminal's capacitance is kept too, λ Q/F, as
    // its current counts dQ/dt: a port or contact at no current that passes K⁺ and charges keeps
    // n_K + Q/F. (Its K⁺ alone held, a membrane port's voltage settled 7 mV from the transient's.)
    // A port's charge, as an electrode's carrier; a contact's, its displacement D_in, below.
    if (circuit >= 2 && model.ports[circuit - 2].capacitance && !carriers.some((p) => p.k === circuit)) carriers = [...carriers, { k: circuit, w: lambda * carrierZ(s, circuit) }];
    // Every contact driven by a current that feeds it and holds a charge, D_in, kept beside the
    // ions it passes at its own λ (two contacts at no current passing K⁺, one a capacitor: its
    // charge was left out where the law took the other's row).
    const lam = (t) => parts.map((p) => (driven(stretches[p.stretch]).includes(t) && z[stretches[p.stretch].species] !== 0 ? p.w / z[stretches[p.stretch].species] : 0)).find((x) => x !== 0);
    const contacts = circuit < 0 ? [] : [0, 1].filter((t) => lam(t) && contactCharges(s, t)).map((t) => ({ t, lambda: lam(t) }));
    let nn = parts.reduce((m, p) => m + stretches[p.stretch].nodes[1] - stretches[p.stretch].nodes[0] + 1, 0);
    for (const p of surface) nn += model.ports[surfCols[p.col][0]].nodes.length * model.ports[surfCols[p.col][0]].surface.length;
    for (const p of carriers) nn += model.ports[p.k - 2].nodes.length;
    for (const p of parts) conserved.add(p.stretch);
    constraints.push({
      parts,
      surface,
      carriers, // (what a port's electrode holds, Q/(zF) at weight w: see conservedMoieties)
      contacts, // (contacts whose charges it keeps, λ D_in/F each)
      terminals: [...new Set([...parts.flatMap((p) => driven(stretches[p.stretch])), ...carriers.map((p) => p.k)])],
      key: parts.map((p) => `${p.stretch}:${p.w}`).join() + surface.map((p) => `|${p.col}:${p.w}`).join() + carriers.map((p) => `|T${p.k}:${p.w}`).join(),
      circuit, // (the terminal whose circuit row it stands in for, or −1 for a balance row's, `row`)
      lambda,
      row: circuit >= 0 ? -1 : blockOfNode[st.nodes[0]] * M + 1 + st.species,
      idx: new Int32Array(nn * M),
      w: new Float64Array(nn * M),
      len: 0,
      res: 0,
    });
  };
  stretches.forEach((st, k) => {
    if (!st.spectator || !moves(s, st)) return;
    add([{ stretch: k, w: 1 }], k);
    // Held flat (see lawRows.js, flatRows), it moves as a whole: the pin's response is the whole
    // stretch's shift, its amount saying how far.
    if (flatStretches.includes(st)) constraints.at(-1).col = Array.from({ length: st.nodes[1] - st.nodes[0] + 1 }, (_, j) => blockOfNode[st.nodes[0] + j] * M + 1 + st.species);
  });
  const S = stretches.length, T = terms.length, X = S + surfCols.length;
  // Two column spaces (see the contract). A law's weights w: stretches [0, S), electrode
  // surfaces' species [S, X), electrodes' charges [X, …). The rows a law can take, v: stretches'
  // balance rows [0, S), then terminals' circuit rows (terminal j − S).
  const isStretch = (k) => k < S, isSurface = (k) => k >= S && k < X, isCarrier = (k) => k >= X;
  const circuitOf = (j) => (j >= S ? j - S : -1);
  const zc = (k) => carrierZ(s, k);
  // The rows replaced: a moving stretch's balance (one of still stretches alone is kept node by
  // node, below), or the circuit row of a terminal driven at no current, which says nothing the
  // combination's balance rows don't (their weighted sum is that terminal's current, weighted
  // by w/z). Together they must be independent, each combination's weights over the rows
  // replaced a nonsingular matrix, or one balance is lost and another said twice (A + B = AB and
  // A + C = AC: −A + B + C in place of B's row, B + AB of AB's and A − B + AC of A's said B and
  // A twice and C not at all, and the solve converged to a state that wasn't steady). So they're
  // chosen by elimination, each combination's weights reduced by the rows already chosen: a
  // circuit row where one is left, else the first stretch with weight 1, else the basis
  // vector's own (its weight is 1 there, and no other vector's is), else the largest.
  const rowOf = (k) => blockOfNode[stretches[k].nodes[0]] * M + 1 + stretches[k].species;
  const chosen = []; // (a column and the reduced weights, per row replaced: S + t for terminal t's)
  const spectators = new Set(constraints.map((cs) => cs.row));
  // A law's weights over the rows it can take: a moving stretch's balance row (not a spectator's,
  // taken above), and the circuit row of each terminal driven by a current that feeds it, at w/z
  // (w/z of an electrode's charge too). The same map for every law, affinities and constraints
  // alike, so that the elimination compares like with like.
  const rowsOf = (w) => {
    const v = new Float64Array(S + T);
    for (let k = 0; k < S; k++) {
      if (w[k] === 0) continue;
      const st = stretches[k];
      if (moves(s, st) && !spectators.has(rowOf(k))) v[k] = w[k];
      if (z[st.species] !== 0) for (const t of driven(st)) v[S + t] = w[k] / z[st.species];
    }
    carriersOf.forEach((k, j) => w[X + j] !== 0 && (v[S + k] = w[X + j] / zc(k)));
    return v;
  };
  // What an electrode holds goes only in its terminal's circuit row (its charge depends on the
  // terminal's voltage, which a balance row's pin can't carry): the combinations through it come
  // first, so that the first takes that row, and the rest have its multiple taken off.
  const through = (w) => carriersOf.some((k, j) => w[X + j] !== 0);
  const held = new Map(); // (terminal → the combination in place of its circuit row)
  // Reactions' equilibria first. A face reaction that alone changes some combination (without
  // it, w would be conserved: Ag⁺ closed in a cell where Ag⁺ + e⁻ = Ag(s) is its one source, or
  // Ag⁺ + AgNO₃ where it also complexes) has that combination's summed balance rows equal to its
  // rate, w·ν times it, which a steady state sets to zero. The rate is a poor row: k₀ times a
  // function of the levels, beside diffusion's terms in the same balances. An Ag⁺/Ag couple
  // 1e-20 as fast as Fe³⁺/Fe²⁺ on the same platinum was lost to their round-off, the
  // open-circuit potential settling where silver wasn't (−0.4661 V, the joint equilibrium
  // −0.5004 V), and 1e-15 to 1e-18 failed. Its affinity says the same at any k₀ (each kinetic
  // form vanishes only there), so it stands in for one of those rows (see lawRows.js,
  // applyAffinities), chosen by the same elimination, and a stretch's: never a circuit row.
  // The slowest first, whose rate rows are the poorest, and each row picked by the reaction's
  // own terms, never by where its species are listed: the same rows in any listing. (Picked by
  // listing, an electrode's two couples took different rows with the species reversed, and
  // Newton converged in one listing and crawled in the other.)
  const affinities = [];
  const faceRx = model.interfaces.flatMap((itf, f) => itf.reactions.map((rx) => ({ f, rx, key: rx.part.map((p) => `${model.species[p.i].name}:${p.side}:${p.nu}`).sort().join() })));
  faceRx.sort((a, b) => a.rx.k0 - b.rx.k0 || (a.key < b.key ? -1 : a.key > b.key ? 1 : a.f - b.f));
  for (const { f, rx } of faceRx) {
    const nu = new Float64Array(S);
    let inFlat = false;
    for (const p of rx.part) {
      const k = stretchOf[(f + p.side) * n + p.i];
      if (k < 0 || flatStretches.includes(stretches[k])) inFlat = true;
      else nu[k] += p.nu;
    }
    if (inFlat) continue; // (held flat at a contact's levels, in equilibrium already)
    // (Not one fed by a terminal driven by a current: its sum is the rate plus what that
    // passes. Silver's, e⁻ + Fe²⁺ with the electrons fed, or Ag⁺ + AgNO₃, which nothing is.
    // A current held at zero passes nothing: an electrode at open circuit feeding a closed
    // pool of Ag⁺ beside a slow Fe³⁺/Fe²⁺ couple settled 78 mV off, its iron not stationary.)
    const passes = (k) => !noCurrent(terms[k]);
    const fed = (w) => w.some((x, k) => x !== 0 && ((isCarrier(k) && passes(carriersOf[k - X])) || (isStretch(k) && driven(stretches[k]).some(passes))));
    // Any of the combinations it alone changes will do, with any conserved amount added (whose
    // summed balance is zero too): one of them needs a row of its own left to take. (Silver's,
    // e⁻ − Fe³⁺ beside a closed iron couple that had Fe³⁺'s, is e⁻ + Fe²⁺ with the iron added.)
    const news = conservedMoieties(s, { drives: true, skip: rx }).filter((w) => Math.abs(nu.reduce((t, x, k) => t + x * w[k], 0)) > 1e-9 && !fed(w));
    // (One it changes only through immobile species has its rate alone in their balance at the
    // face's node already, a row with no diffusion beside it to lose it to: left as it is. Stood
    // in for by a mobile partner, Ag⁺ + Cl⁻ = S(immobile) beside a floating Ag electrode said its
    // rate twice and lost the electrode's I = 0.)
    if (news.some((w) => [...w.keys()].every((k) => !isStretch(k) || w[k] === 0 || still(s, stretches[k])))) continue;
    const kept = laws.filter((l) => !fed(l));
    const tries = [...news, ...news.flatMap((w) => kept.flatMap((l) => [1, -1].map((sg) => w.map((x, m) => x + sg * l[m]))))];
    for (const w of tries) {
      const v = rowsOf(w);
      for (const [col, row] of chosen) {
        const fc = v[col];
        if (fc !== 0) for (let j = 0; j < S + T; j++) v[j] -= fc * row[j];
      }
      const scale = v.reduce((m, x) => Math.max(m, Math.abs(x)), 0);
      if (!(scale > 1e-9)) continue;
      // (Never a metal's carriers, nor a stretch a terminal feeds, even at no current: their rows
      // carry the terminal's current, and Pt's electrons standing for an iron couple's balance
      // put the electrode at 602 V.)
      const ok = (j) => j !== undefined && j >= 0 && j < S && Math.abs(v[j]) > 1e-9 * scale && !conductor(stretches[j]) && terms.every((_, t) => !feeds(t, stretches[j]));
      // (One of its own species first, then weight 1, the largest, and by name.)
      const own = (k) => rx.part.some((p) => p.i === stretches[k].species && stretchOf[(f + p.side) * n + p.i] === k);
      const key = (k) => [own(k) ? 0 : 1, Math.abs(v[k]) === 1 ? 0 : 1, -Math.abs(v[k]), model.species[stretches[k].species].name, stretches[k].regions[0]];
      const less = (a, b) => { const ka = key(a), kb = key(b); for (let j = 0; j < ka.length; j++) if (ka[j] !== kb[j]) return ka[j] < kb[j]; return false; };
      const col = [...Array(S).keys()].filter(ok).reduce((b, k) => (b < 0 || less(k, b) ? k : b), -1);
      if (col < 0) continue;
      chosen.push([col, v.map((x) => x / v[col])]);
      affinities.push({ rx, f, row: rowOf(col), idx: new Int32Array(rx.part.length), w: new Float64Array(rx.part.length), len: 0, res: 0 });
      break;
    }
  }
  // A contact's charge has no column: a law through what a contact driven by a current passes,
  // where it holds a charge, counts λ D_in/F of it at λ = w/z on what it passes (see add). Kept
  // once, in a circuit row, like an electrode's: the first law through it takes one, and the rest
  // have its multiple taken off. (Two laws through K⁺ at a capacitor at I = 0, K⁺ − L⁻ and K⁺ + KL
  // under K⁺ + L⁻ = KL: the second took KL's balance row without the charge, keeping n_K + n_KL,
  // which the transient moves 15×, and the solve failed.)
  const lamOf = (w, t) => (contactCharges(s, t) ? (stretches.map((st, k) => (w[k] !== 0 && z[st.species] !== 0 && driven(st).includes(t) ? w[k] / z[st.species] : 0)).find((x) => x !== 0) ?? 0) : 0);
  const heldContact = new Map(); // (contact → the law that keeps its charge)
  const unkept = [];
  const reduce = (w, c, f) => {
    const r = w.map((x, m) => (Math.abs(x - f * c[m]) < 1e-9 ? 0 : x - f * c[m]));
    r.own = w.own;
    return r;
  };
  for (let w of [...laws.filter(through), ...laws.filter((w) => !through(w))]) {
    for (const [j, k] of carriersOf.entries()) {
      const c = held.get(k);
      if (w[X + j] === 0 || !c) continue;
      w = reduce(w, c, w[X + j] / c[X + j]);
    }
    for (const [t, c] of heldContact) if (lamOf(w, t) !== 0) w = reduce(w, c, lamOf(w, t) / lamOf(c, t));
    const parts = [...w.keys()].filter((k) => isStretch(k) && w[k] !== 0).map((k) => ({ stretch: k, w: w[k] }));
    const surface = [...w.keys()].filter((k) => isSurface(k) && w[k] !== 0).map((k) => ({ col: k - S, w: w[k] }));
    const carriers = carriersOf.flatMap((k, j) => (w[X + j] !== 0 ? [{ k, w: w[X + j] }] : []));
    // DROP: a surface's own is kept node by node, below. (One with an electrode's charge in it
    // too is the circuit's, here: below skips it, and a floating electrode whose adsorbates gave
    // it electrons settled 15 mV from where the transient ends.)
    if (parts.length === 0 && carriers.length === 0) continue;
    // UNKEPT: through a stretch still in part only, which no row here can carry (a known gap).
    if (parts.some((p) => !moves(s, stretches[p.stretch]) && !still(s, stretches[p.stretch]))) {
      unkept.push(w);
      continue;
    }
    // DROP: a spectator's own amount: added above where it moves, node by node below where it's
    // still.
    if (parts.length === 1 && surface.length === 0 && stretches[parts[0].stretch].spectator) continue;
    const v = rowsOf(w);
    for (const [col, row] of chosen) {
      const f = v[col];
      if (f !== 0) for (let j = 0; j < S + T; j++) v[j] -= f * row[j];
    }
    const scale = v.reduce((m, x) => Math.max(m, Math.abs(x)), 0);
    // DROP: no row of its own left once those already chosen are taken off: of still stretches
    // alone (kept node by node, below), or said by the laws already kept.
    if (!(scale > 1e-9)) continue;
    const ok = (j) => j !== undefined && Math.abs(v[j]) > 1e-9 * scale;
    const col = [...Array(T).keys()].map((t) => S + t).find(ok) ?? [parts.find((p) => w[p.stretch] === 1)?.stretch, w.own].find(ok) ?? v.findIndex((x) => Math.abs(x) === scale);
    // UNKEPT: a charge (an electrode's on its capacitance, or a contact's) in a balance row,
    // whose pin can't carry the terminal's voltage it depends on, nor (a contact's) hold it at
    // all: none left to take off, no circuit row left to take.
    if (circuitOf(col) < 0 && (carriers.some(({ k }) => model.ports[k - 2].capacitance) || [0, 1].some((t) => lamOf(w, t) !== 0))) {
      unkept.push(w);
      continue;
    }
    // (Other electrodes' charges in a circuit row not their own: their voltages enter it through
    // termDX. Two floating electrodes sharing a closed Fe³⁺, their law dropped, settled at 1.3 V
    // for −0.65.)
    chosen.push([col, v.map((x) => x / v[col])]);
    if (circuitOf(col) >= 0) held.set(circuitOf(col), w);
    for (const t of [0, 1]) if (circuitOf(col) >= 0 && lamOf(w, t) !== 0 && !heldContact.has(t)) heldContact.set(t, w);
    add(parts, circuitOf(col) < 0 ? col : parts[0]?.stretch, surface, circuitOf(col), carriers);
  }
  const contactLaw = [0, 1].map((k) => constraints.some((cs) => cs.contacts.some((c) => c.t === k)));
  // A combination that nothing carries anywhere is conserved node by node: of still stretches
  // (trap states: X⁰ + X⁻ under e⁻ + X⁰ = X⁻), or of an electrode's surface species (adsorbates
  // that only turn into each other). At each node its weighted sum stays what it was, a row local
  // to that node's block, in place of one of their balances.
  const localConstraints = [];
  const span = ([g0, g1]) => Array.from({ length: g1 - g0 + 1 }, (_, j) => g0 + j);
  const picked = new Map(); // (node → the rows chosen there: [slot, the reduced weights by slot])
  for (const w of laws) {
    const parts = [...w.keys()].filter((k) => w[k] !== 0);
    if (parts.some(isCarrier)) continue; // (through an electrode: above)
    // (Fed by a terminal driven by a current: in its circuit row, above. A port spread over
    // several nodes carries it between them too. Kept here as well, S⁺ that a port at I = 0
    // links was held twice, and the port's I = 0 lost: the solve settled at 0 V for −1.42 V.)
    if (parts.some((k) => isStretch(k) && terms.some((_, t) => feeds(t, stretches[k])))) continue;
    const sets = parts.map((k) => (isSurface(k) ? model.ports[surfCols[k - S][0]].nodes : still(s, stretches[k]) ? span(stretches[k].nodes) : null));
    if (sets.includes(null)) continue; // (through a moving stretch: above)
    // Its nodes: the immobile stretches', side by side, with any surface's window among them (a
    // surface's term only where it is: an immobile S that adsorbs on an electrode spread over
    // part of its region keeps S + Γaθ there and S alone elsewhere; required to match, that law
    // was lost, and a floating electrode's coverage drifted). Or a surface's alone.
    const own = parts.filter(isStretch), base = own.length > 0 ? sets[parts.indexOf(own[0])] : sets[0];
    const same = (set) => set.length === base.length && set.every((g, j) => g === base[j]);
    // UNKEPT: immobile stretches over different nodes, or a surface outside them: a face
    // reaction joining the two nodes of a face, whose sum it keeps and neither alone.
    if (!parts.every((k, j) => (isStretch(k) ? same(sets[j]) : sets[j].every((g) => base.includes(g))))) {
      unkept.push(w);
      continue;
    }
    for (const k of own) conserved.add(k);
    const slot = (k) => (isStretch(k) ? 1 + stretches[k].species : 1 + n + surfCols[k - S][1]);
    // Its row at each node, chosen by elimination as the constraints' are: its weights over the
    // slots there (a surface's only where it is) reduced by the rows already chosen there, then
    // the slot it prefers (a stretch's where it has one, present at every node, its own column's
    // first; else the surface's) where that's left, else the first that is. One said already
    // there, by those before it, is left out there. (Both given the same row, immobile S⁺ with
    // S⁺ + A(ads) + e⁻ = B(ads): A − S⁺ and B + S⁺ took S⁺'s, one was lost, and θ_A settled at 1
    // for the transient's 0.3.)
    const rowPart = own.length > 0 ? (own.includes(w.own) ? w.own : own[0]) : parts.includes(w.own) ? w.own : parts[0];
    const rowsAt = new Map(); // (slot → its nodes)
    for (const g of base) {
      const v = new Map();
      parts.forEach((k, j) => sets[j].includes(g) && v.set(slot(k), w[k]));
      for (const [sl, row] of picked.get(g) ?? []) {
        const f = v.get(sl);
        if (f) for (const [sx, x] of row) v.set(sx, (v.get(sx) ?? 0) - f * x);
      }
      const scale = [...v.values()].reduce((m, x) => Math.max(m, Math.abs(x)), 0);
      if (!(scale > 1e-9)) continue;
      const ok = (sl) => Math.abs(v.get(sl) ?? 0) > 1e-9 * scale;
      const sl = ok(slot(rowPart)) ? slot(rowPart) : [...v.keys()].sort((a, b) => a - b).find(ok);
      picked.set(g, [...(picked.get(g) ?? []), [sl, new Map([...v].map(([sx, x]) => [sx, x / v.get(sl)]))]]);
      rowsAt.set(sl, [...(rowsAt.get(sl) ?? []), g]);
    }
    for (const [row, nodes] of rowsAt) {
      localConstraints.push({
        parts: parts.map((k) => (isStretch(k) ? { i: stretches[k].species, w: w[k] } : { q: surfCols[k - S][1], k: surfCols[k - S][0], w: w[k] })),
        row,
        nodes,
        reference: new Float64Array(nodes.length),
      });
    }
  }
  // An immobile stretch at a node where nothing touches it (no reaction there, no port or contact
  // linking it) keeps its concentration: its balance there is 0 = 0 in a steady state, which
  // holds nothing, and a charged one's level stayed where it was while φ moved (S⁺ beside a face
  // that makes it, the bath moved 0.2 V: S⁺ fell to e^{−0.2/V_T} of itself away from the face).
  for (const st of stretches) {
    if (!still(s, st)) continue;
    const i = st.species, touched = new Set();
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) if (rxsIn[nodeMaterial[g]].some((rx) => rx.sp.includes(i))) touched.add(g);
    model.interfaces.forEach((itf, f) => {
      for (const rx of itf.reactions) for (const p of rx.part) if (p.i === i) touched.add(p.side === 0 ? grid.regionEnd[f] : grid.regionStart[f + 1]);
    });
    for (const port of model.ports) {
      if (port.species[i].type !== 'blocked' || port.reactions.some((rx) => rx.part.some((p) => p.side === 0 && p.i === i))) for (const g of port.nodes) touched.add(g);
    }
    if (st.leftOpen) touched.add(0);
    if (st.rightOpen) touched.add(nNodes - 1);
    for (const lc of localConstraints) if (lc.row === 1 + i) for (const g of lc.nodes) touched.add(g);
    const free = span(st.nodes).filter((g) => !touched.has(g));
    if (free.length > 0) localConstraints.push({ parts: [{ i, w: 1 }], row: 1 + i, nodes: free, reference: new Float64Array(free.length) });
  }
  // Islands: the pieces of a stretch between its faces. One that nothing else feeds (no
  // contact, port or reaction) holds its level only through its faces' fluxes, which can be
  // ~1e-24 of its own conduction (a tiny conductance at the face, or a neighbour that barely
  // conducts): eliminated, the level is lost to round-off, and Newton settles at a wrong one with
  // nothing to show for it. In steady solves the piece's balance summed over its nodes, where
  // the internal fluxes cancel exactly and the faces' alone are left, replaces its first node's
  // (see lawRows.js, applyIslands), once the plain solves converge (see Solver.newton).
  const islands = [];
  const used = new Set([...constraints.map((cs) => cs.row).filter((o) => o >= 0), ...affinities.map((af) => af.row)]);
  for (const st of stretches) {
    const i = st.species, [s0, s1] = st.regions;
    if (!st.mobile || conductor(st)) continue;
    for (let q0 = s0, q1 = s0; q1 <= s1; q1++) {
      const p0 = q0;
      q0 = q1 + 1;
      if (p0 === s0 && q1 === s1) continue; // (uncut)
      if ((p0 === s0 && st.leftOpen) || (q1 === s1 && st.rightOpen)) continue;
      const inside = (q) => q >= p0 && q <= q1;
      if (model.ports.some((port) => inside(port.region) && (port.species[i].type !== 'blocked' || port.reactions.some((rx) => rx.part.some((p) => p.side === 0 && p.i === i))))) continue;
      let reactive = false;
      for (let q = p0; q <= q1; q++) {
        for (const rx of model.reactions) if (rx.kf[regions[q].material] > 0 && [...rx.reactants, ...rx.products].some((p) => p.i === i)) reactive = true;
      }
      for (let f = Math.max(0, p0 - 1); f <= Math.min(model.interfaces.length - 1, q1); f++) {
        for (const rx of model.interfaces[f].reactions) if (rx.part.some((p) => p.i === i && inside(f + p.side))) reactive = true;
      }
      const row = blockOfNode[grid.regionStart[p0]] * M + 1 + i;
      if (reactive || used.has(row)) continue;
      used.add(row);
      // The faces' fluxes into the piece: −A_f u at its left face, +A_f u at its right. Across
      // a face that holds μ̄ level, u is set only by the edge balances, and eliminated, it can
      // come out of the piece's own (G Δη, its two η all but equal: a neighbour 1e20 times
      // less conductive, and the piece's level was seen through one face of two, and Newton
      // swung it back and forth across both). There the outside edge node's balance is added:
      // u cancels, and the outside's first segment carries the flux instead.
      const faces = [], outside = [];
      if (p0 > s0) faces.push([p0 - 1, -1, grid.regionEnd[p0 - 1]]);
      if (q1 < s1) faces.push([q1, 1, grid.regionStart[q1 + 1]]);
      for (const [f, , g] of faces) {
        if (model.interfaces[f].links[i].type !== 'equilibrium' || g === 0 || g === nNodes - 1) continue;
        if (model.ports.some((port) => port.nodes.includes(g))) continue;
        outside.push(blockOfNode[g] * M + 1 + i);
      }
      islands.push({
        stretch: st,
        row,
        flux: faces.map(([f]) => blockOfFace[f] * M + 1 + i),
        weight: faces.map(([f, sg]) => sg * grid.area[grid.regionEnd[f]]),
        outside,
        idx: new Int32Array(faces.length + 3 * M * outside.length),
        w: new Float64Array(faces.length + 3 * M * outside.length),
        len: 0,
        res: 0,
      });
    }
  }
  for (const isl of islands) isl.outside = isl.outside.filter((o) => !used.has(o));

  // Time steps only: a terminal behind a resistance that alone feeds a combination reads its
  // current by that combination's change too (see lawRows.js, circuitLaws), its combinations found
  // with it counted as driven by a current. (Not in the steady system, where it holds its level.)
  const timeLaws = [];
  const behind = floating.filter((k) => terms[k].drive.kind !== 'I');
  if (behind.length > 0) {
    const more = conservedMoieties(s, { drives: true, alsoCurrent: behind }), X2 = S + surfCols.length;
    for (const w of more) {
      const parts = [...w.keys()].filter((k) => k < S && w[k] !== 0).map((k) => ({ stretch: k, w: w[k] }));
      if (parts.length === 0 || parts.some((p) => !moves(s, stretches[p.stretch]) && !still(s, stretches[p.stretch]))) continue;
      const surface = [...w.keys()].filter((k) => k >= S && k < X2 && w[k] !== 0).map((k) => ({ col: k - S, w: w[k] }));
      const carriers = more.carriers.flatMap((k, j) => (w[X2 + j] !== 0 ? [{ k, w: w[X2 + j] }] : []));
      // Fed by one terminal behind a resistance and nothing else that passes current.
      const fedBy = new Set([...carriers.map((p) => p.k), ...parts.flatMap((p) => floating.filter((k) => feeds(k, stretches[p.stretch])))]);
      if (fedBy.size !== 1) continue;
      const [k] = fedBy;
      if (!behind.includes(k) || timeLaws.some((cs) => cs.circuit === k)) continue;
      // (Nor a contact that holds a charge: its current is what it passes and dD_in/dt, which the
      // change leaves out. Read by it, a capacitive contact behind a resistance passing a trace
      // of K⁺ failed every step.)
      if (k < 2 && contactCharges(s, k)) continue;
      const zk = k >= 2 ? carrierZ(s, k) : 0;
      const lambda = carriers.find((p) => p.k === k)?.w / zk || parts.map((p) => (feeds(k, stretches[p.stretch]) && z[stretches[p.stretch].species] !== 0 ? p.w / z[stretches[p.stretch].species] : 0)).find((x) => x !== 0);
      if (!lambda) continue;
      let nn = parts.reduce((m, p) => m + stretches[p.stretch].nodes[1] - stretches[p.stretch].nodes[0] + 1, 0);
      for (const p of surface) nn += model.ports[surfCols[p.col][0]].nodes.length * model.ports[surfCols[p.col][0]].surface.length;
      for (const p of carriers) nn += model.ports[p.k - 2].nodes.length;
      timeLaws.push({ parts, surface, carriers, terminals: [k], circuit: k, lambda, idx: new Int32Array(nn * M), w: new Float64Array(nn * M), len: 0 });
    }
  }
  return { constraints, affinities, localConstraints, islands, timeLaws, contactLaw, carriers: carriersOf, conserved, unkept };
}
