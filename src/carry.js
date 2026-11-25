// Carrying a solved state across Device.set() (docs/device.md, "set(patch)"): from the old solver
// onto a new one built at the new definition, by what each value means, never by where it sat.
//
// The contract:
//   1. The new solver starts from its own start (its c0, gates at rest, coverages at θ₀, at the
//      sources then), and every value that can be matched by meaning overwrites it: a region by its
//      name (given once; by place where it has none, or was renamed: the same count, the same
//      index), a node by its fraction along its region, a species by its name, a face by the
//      regions on its two sides (in order), a face's reaction by its participants, a gate or a
//      surface species by its name, a surface or a port by its owner's name, a terminal by its
//      name; a slot only where both solvers solve it. What can't be matched keeps its own start
//      (a gate at rest at the carried voltage); a region that had a state and gets none (or turns
//      from a metal to a non-metal, or back) says so.
//   2. What's carried is physical: concentrations, φ, coverages, gates' fractions, a face's fluxes,
//      rates and displacement, a double layer's charge (its drop found), a floating capacitor's
//      charge (its voltage found), a metal's flux and its carrier's level. Each is written back
//      through the new definition's own relations (its statistics, inverted; its μ°, c_ref, T), so
//      a change of any of them moves no concentration.
//   3. Between old nodes, a species follows the profile Scharfetter–Gummel takes within the old
//      cell (c e^θ̃ linear in ∫e^θ̃, θ̃ the drift potential the flux uses: zφ̂, the grading, the
//      excess potential, less a flow's Péclet shift), in the cell's own coordinate (∫dx/A), and
//      in logarithms; an immobile or eddy-mixed one linearly. Coverages and gates are linear,
//      within their own window (one moved at its width carries them along; resized, in place).
//   4. A node that nothing changed (the same place, to a few ulps; material data, statistics, T)
//      is copied bit for bit, low parts and all.
//   5. Amounts: each stretch's, region by region, what the old state held over a region matched to
//      one that held it (where resized or reshaped, concentrations kept: its volume's part scaled
//      node by node, a diffuse layer's excess by its plane's area), its own start's where nothing carries; a stretch closed on its own over the same space keeps its reference
//      exactly. The state is put at those amounts (a species scaled, the others held) unless a held
//      terminal sets its level; a surface's carried nodes likewise, all its species at once (in any
//      listing order, the same state), where its window and electrode are unchanged, refused where
//      the sites can't hold it (by more than its amounts' round-off). A strictly neutral node is brought onto its Gauss row with the closed amounts and
//      the floating capacitors' charges in one solve (src/neutral.js, solveNeutral): at its levels,
//      including the shared-site derivatives on a regrid of a neutral lattice where every partner
//      has an amount target; or, under a capacitor whose hardware changed, at its concentrations
//      (φ̂ jumping). Each held
//      capacitor's start charge (a gate's, a capacitive contact's, a metal's sheet at a capacitive
//      face, sheetDEnd) is what it held before the jump where its hardware's group changed (and a
//      held gate's whose own hardware changed over a region with ε > 0, where nothing jumps), else
//      what the final state holds at its old drive; with a regrid, the group's starts are shared
//      onto what its ions balance, by capacitance. Terminals' running charges
//      and stored charges scale with their area (a contact's, a port's electrode where its
//      capacitance holds charge); a port's stored charge, node by node (its σ per area kept), and a
//      floating one's voltage holds it on the final state.
//
// carryState(old, s, { oldDef, newDef, pure }) writes s's state (and its clock, terminals and
// references); it changes nothing else, so a set() that throws from here leaves the device as it
// was. `pure`: the definitions differ only in the grid (a steady state then stands for the old
// definition's, solved again on the new grid before the next step: see Device.advance).

import { crossSection } from './grid.js';
import { DeviceError } from './device.js';
import { gateRate } from './boundaries.js';
import { FARADAY } from './constants.js';
import { deepEqual } from './merge.js';
import { strictlyNeutralNodes, solveNeutral } from './neutral.js';

// [ln(1 − s), ln s] for s = expm1(d w)/expm1(d), the weight of the right end at fraction w of a cell
// whose drift potential rises by d: exact in logarithms at any d.
export function sgLogWeights(d, w) {
  if (w <= 0) return [0, -Infinity];
  if (w >= 1) return [-Infinity, 0];
  if (Math.abs(d) < 1e-12) return [Math.log1p(-w), Math.log(w)];
  if (d > 0) {
    // 1 − s = (1 − e^{−d(1−w)})/(1 − e^{−d}), s = e^{d(w−1)}(1 − e^{−dw})/(1 − e^{−d}): directly (through
    // 1 − w, a small w's weight rounded away; Astra's second pass)
    const den = Math.log(-Math.expm1(-d));
    return [Math.log(-Math.expm1(-d * (1 - w))) - den, d * (w - 1) + Math.log(-Math.expm1(-d * w)) - den];
  }
  // d < 0: 1 − s = e^{dw}(1 − e^{d(1−w)})/(1 − e^d), s = (1 − e^{dw})/(1 − e^d)
  const den = Math.log(-Math.expm1(d));
  return [d * w + Math.log(-Math.expm1(d * (1 - w))) - den, Math.log(-Math.expm1(d * w)) - den];
}

// c at fraction w of a cell from c_a to c_b with drift rise d (see the contract, 3).
export function sgInterpolate(ca, cb, d, w) {
  const [l0, l1] = sgLogWeights(d, w);
  const p = ca > 0 ? l0 + Math.log(ca) : -Infinity, q = cb > 0 ? l1 + Math.log(cb) + d : -Infinity;
  const m = Math.max(p, q);
  if (m === -Infinity) return 0;
  return Math.exp(m + Math.log(Math.exp(p - m) + Math.exp(q - m)) - d * w);
}

// (a name only where it's given, and given once: duplicates are matched by place)
const regionNamed = (def, r) => {
  const nm = def?.regions?.[r]?.name;
  return nm !== undefined && def.regions.filter((rg) => rg.name === nm).length === 1;
};

// Each new region's old region (−1: none). By an explicit name first; then, where the counts are
// equal, by place (a renamed or unnamed region), unless that old region is claimed by name.
export function matchRegions(oldModel, newModel, oldDef, newDef) {
  const map = new Array(newModel.regions.length).fill(-1), claimed = new Set();
  newModel.regions.forEach((rg, r) => {
    if (!regionNamed(newDef, r)) return;
    const q = oldModel.regions.findIndex((o, k) => regionNamed(oldDef, k) && o.name === rg.name);
    if (q >= 0 && !claimed.has(q)) (map[r] = q), claimed.add(q);
  });
  if (oldModel.regions.length === newModel.regions.length) {
    map.forEach((q, r) => {
      if (q >= 0 || claimed.has(r)) return;
      // (a name that's in the old list elsewhere isn't a rename: it moved, and was matched above)
      if (regionNamed(newDef, r) && oldModel.regions.some((o, k) => regionNamed(oldDef, k) && o.name === newModel.regions[r].name)) return;
      map[r] = r;
      claimed.add(r);
    });
  }
  return map;
}

// A face's reaction, by what it does: its participants (species, side, ν), and which of
// several alike it is.
const reactionKeys = (itf, species) => {
  const seen = new Map();
  return itf.reactions.map((rx) => {
    // (its surface participants and fixed share too: two pathways alike in their mobile species,
    // Fc → FcClO₄ and Fc → FcB, matched by order, swapped rates; Astra's second pass)
    const k = JSON.stringify([rx.part.map((p) => [species[p.i]?.name ?? p.i, p.side, p.nu]).sort(), (rx.surf ?? []).map((p) => [itf.surface?.[p.s]?.name ?? p.s, p.nu]).sort(), rx.fixedA ?? 0]);
    const j = seen.get(k) ?? 0;
    seen.set(k, j + 1);
    return JSON.stringify([k, j]);
  });
};

// Where a new node sits in the old region: its old position, and the bracketing old nodes a ≤ b
// with the weight w of b in the cell's own coordinate (∫dx/A), and wx in x.
function locate(old, q, xo) {
  const og = old.model.grid, g0 = og.regionStart[q], g1 = og.regionEnd[q];
  if (xo <= og.x[g0]) return { a: g0, b: g0, w: 0, wx: 0 };
  if (xo >= og.x[g1]) return { a: g1, b: g1, w: 0, wx: 0 };
  let a = g0;
  while (a + 1 < g1 && og.x[a + 1] < xo) a++;
  const b = a + 1, h = og.x[b] - og.x[a];
  if (!(h > 0) || xo === og.x[a]) return { a, b: a, w: 0, wx: 0 };
  if (xo === og.x[b]) return { a: b, b, w: 0, wx: 0 };
  const wx = (xo - og.x[a]) / h;
  let w = wx;
  const geo = old.model.geometry;
  if (geo && geo.type !== 'planar' && og.area[a] > 0 && og.area[b] > 0) {
    const A = crossSection(geo), all = A.resistance(og.x[a], og.x[b]);
    if (all > 0) w = A.resistance(og.x[a], xo) / all;
  }
  return { a, b, w, wx };
}

// Where a node of a new window sits in the old window of the same owner, among the old window's
// nodes only (beyond the last of them, still within the window, at it): { a, b, wx }, or null
// where the old window doesn't reach (its own start). The windows, [lo, hi] in x, as declared. A window
// of the same width, moved, carries its coverages and gates with it (an electrode moved); one
// widened or narrowed keeps them where they were (`toOld`: the node's place in the old device),
// its new part at its own start.
function windowLocate(og, oldNodes, [lo, hi], ng, [nlo, nhi], g, toOld) {
  let xo;
  if (lo === nlo && hi === nhi) xo = ng.x[g];
  else if (Math.abs(nhi - nlo - (hi - lo)) <= 1e-12 * Math.max(hi - lo, nhi - nlo)) xo = lo + (ng.x[g] - nlo);
  else {
    xo = toOld(g);
    if (!(xo >= lo - 1e-12 * (hi - lo)) || !(xo <= hi + 1e-12 * (hi - lo))) return null;
  }
  const order = Array.from(oldNodes).sort((p, q) => og.x[p] - og.x[q]);
  if (order.length === 0) return null;
  if (xo <= og.x[order[0]]) return { a: order[0], b: order[0], wx: 0 };
  for (let k = 0; k + 1 < order.length; k++) {
    const a = order[k], b = order[k + 1];
    if (xo === og.x[a]) return { a, b: a, wx: 0 };
    if (xo === og.x[b]) return { a: b, b, wx: 0 };
    if (xo < og.x[b]) return og.x[b] > og.x[a] ? { a, b, wx: (xo - og.x[a]) / (og.x[b] - og.x[a]) } : { a, b: a, wx: 0 };
  }
  const last = order[order.length - 1];
  return { a: last, b: last, wx: 0 };
}

// The model's statistics for a material, as the definition gave them (custom functions by
// identity): what decides whether a node's concentrations read the same off its levels.
const statsOf = (def, model, m) => def?.materials?.[model.materials[m].name]?.statistics;

/**
 * Carry old's state onto s (see the contract above). Returns what it matched and copied: `map`,
 * each new region's old one (−1: none, started from its c0); `copied`, the nodes copied bit for bit.
 */
export function carryState(old, s, { oldDef, newDef, pure = false }) {
  const { n, M } = s, oldM = old.M, species = s.model.species;
  // each species' index in the old solver, by name (the same set, in any order: see Device._rebuild)
  const sp = species.map((x) => old.model.species.findIndex((y) => y.name === x.name));
  const sameOrder = sp.every((o, i) => o === i);
  const VTo = old.model.RT, VTn = s.model.RT; // (RT: φ̂ = Fφ/RT, so φ = φ̂·RT/F; the F cancels)
  const sameT = VTo === VTn;
  old.computeConcentrations();
  old._assembleBookkeeping(Infinity); // (its ports' charges, portQ, as the old state holds them: see portChargeAcross)
  s.computeConcentrations();
  const fresh = Float64Array.from(s.u), freshC = Float64Array.from(s.c), freshTh = Float64Array.from(s.th);
  const map = matchRegions(old.model, s.model, oldDef, newDef);
  const ng = s.model.grid, og = old.model.grid;
  const sameGeometry = deepEqual(old.model.geometry, s.model.geometry);
  const copied = new Uint8Array(s.nNodes); // (nodes carried bit for bit)
  const moved = new Uint8Array(s.nNodes); // (nodes whose concentrations were interpolated, not carried as they were)
  const target = new Float64Array(s.nNodes * n).fill(NaN); // (concentrations to write, where carried)
  const phiAt = new Float64Array(s.nNodes).fill(NaN);
  const levelAt = new Float64Array(s.nNodes * n).fill(NaN); // (each species' level μ̄/RT·RT_old, linear: for a level that's its own)
  const sameStats = new Map(), kindChanged = new Set();
  const statsSame = (mo, mn) => {
    const k = `${mo}:${mn}`;
    if (!sameStats.has(k)) sameStats.set(k, old.model.materials[mo].name === s.model.materials[mn].name && deepEqual(statsOf(oldDef, old.model, mo), statsOf(newDef, s.model, mn)) && deepEqual(old.model.materials[mo].conductor, s.model.materials[mn].conductor) && (old.model.materials[mo].ideal || old.model.materials[mo].epsr === s.model.materials[mn].epsr));
    return sameStats.get(k);
  };

  // --- nodes
  for (let r = 0; r < s.model.regions.length; r++) {
    const q = map[r];
    if (q < 0) continue;
    const L = s.model.regions[r].length, Lo = old.model.regions[q].length;
    const x0 = ng.x[ng.regionStart[r]], xo0 = og.x[og.regionStart[q]];
    for (let g = ng.regionStart[r]; g <= ng.regionEnd[r]; g++) {
      const t = L > 0 ? (ng.x[g] - x0) / L : 0, xo = xo0 + t * Lo;
      let { a, b, w, wx } = locate(old, q, xo);
      // 4. an old node at the very same place (by its offset in the region), where nothing changed: copy
      let at = -1;
      // (equal to a few ulps of the region: an upstream region resized moves the subtraction's rounding)
      const tolX = 4 * Number.EPSILON * Math.max(Math.abs(ng.x[g]), Math.abs(og.x[a]), Math.abs(og.x[b]), L);
      if (L === Lo) for (const k of [a, b]) if (Math.abs(og.x[k] - xo0 - (ng.x[g] - x0)) <= tolX) at = k;
      if (at >= 0) (a = b = at), (w = wx = 0);
      const bg = s.blockOfNode[g], ba = old.blockOfNode[a], bb = old.blockOfNode[b];
      // (a metal made a semiconductor, or back: nothing physical carries; its own start, and said,
      // its old net charge too: kept, a Pt made water had the metal's sheet charge put on its new
      // ions, K⁺ 100.10 / Cl⁻ 99.90 at its edge, while the solution said it started from its c0;
      // the inventory of duplicated quantities)
      if ((s.nodeConductor[g] >= 0) !== (old.nodeConductor[a] >= 0) || (s.nodeConductor[g] >= 0) !== (old.nodeConductor[b] >= 0)) {
        kindChanged.add(r);
        continue;
      }
      let exact = at >= 0 && sameT && statsSame(old.nodeMaterial[a], s.nodeMaterial[g]) && old.nodeConductor[a] === s.nodeConductor[g] && old.nodeIdeal[a] === s.nodeIdeal[g] && old.background[a] === s.background[g];
      if (exact) for (let i = 0; i < n && exact; i++) {
        const ko = a * n + sp[i], kn = g * n + i;
        if (old.present[ko] !== s.present[kn] || (s.present[kn] && (!Object.is(old.cRef[ko], s.cRef[kn]) || !Object.is(old.mu0hat[ko], s.mu0hat[kn])))) exact = false;
      }
      if (exact) {
        s.u[bg * M] = old.u[ba * oldM];
        s.uLo[bg * M] = old.uLo[ba * oldM];
        for (let i = 0; i < n; i++) {
          s.u[bg * M + 1 + i] = old.u[ba * oldM + 1 + sp[i]];
          s.uLo[bg * M + 1 + i] = old.uLo[ba * oldM + 1 + sp[i]];
        }
        copied[g] = 1;
        continue;
      }
      // a metal: its flux J (linear), its carrier's level μ̄ (linear in volts)
      const im = s.nodeConductor[g];
      if (im >= 0) {
        if (old.nodeConductor[a] === sp[im] && old.nodeConductor[b] === sp[im]) {
          const v = (bk, k) => old.u[bk * oldM + k] + old.uLo[bk * oldM + k];
          s.u[bg * M] = (1 - wx) * v(ba, 0) + wx * v(bb, 0);
          s.u[bg * M + 1 + im] = (((1 - wx) * v(ba, 1 + sp[im]) + wx * v(bb, 1 + sp[im])) * VTo) / VTn;
          s.uLo[bg * M] = s.uLo[bg * M + 1 + im] = 0;
        }
        continue;
      }
      const ph = (bk) => old.u[bk * oldM] + old.uLo[bk * oldM];
      const phi = ((1 - w) * ph(ba) + w * ph(bb)) * VTo; // (in RT/F units of the old: physical φ·F)
      phiAt[g] = phi;
      const seg = a; // (the old segment from a to b = a + 1)
      const rg = old.model.regions[q], vel = rg.velocity ?? 0, h = b > a ? og.segLength[seg] / og.segArea[seg] : 0;
      for (let i = 0; i < n; i++) {
        if (!s.present[g * n + i]) continue;
        const oi = sp[i];
        if (!old.present[a * n + oi] || !old.present[b * n + oi]) continue; // (newly present: its own start)
        const ca = old.c[a * n + oi], cb = old.c[b * n + oi];
        levelAt[g * n + i] = ((1 - w) * (old.u[ba * oldM + 1 + oi] + old.uLo[ba * oldM + 1 + oi]) + w * (old.u[bb * oldM + 1 + oi] + old.uLo[bb * oldM + 1 + oi])) * VTo;
        let c = ca;
        if (b > a) {
          moved[g] = 1;
          // (eddy mixing adds to diffusion what it adds: the flow's Péclet number over D + D_mix, the
          // drift's share D/(D + D_mix), exact for a neutral species; an immobile one has no profile
          // to follow, linear, which keeps local sums linear too: along SG, immobile A ⇌ B moved its
          // local totals; Astra's second pass)
          const k = seg * n + oi, D = old.segD[k], mix = rg.mixing > 0 ? rg.mixing : 0, Dt = D + mix;
          const d = (D > 0 ? (D / Dt) * (old.z[oi] * (ph(bb) - ph(ba)) + old.segShift[k] + old.ex[k + n] - old.ex[k]) : 0) - (Dt > 0 ? (vel * h) / Dt : 0);
          c = D > 0 ? sgInterpolate(ca, cb, d, w) : (1 - wx) * ca + wx * cb;
        }
        target[g * n + i] = c;
      }
    }
  }
  // Write the carried concentrations through the new definition: φ̂, then each species' level
  // (ideal: ln(c/c_ref) + μ°/RT + zφ̂; otherwise each model inverted at them, the new species'
  // own start among them).
  for (let g = 0; g < s.nNodes; g++) {
    if (copied[g] || Number.isNaN(phiAt[g])) continue;
    const bg = s.blockOfNode[g];
    const phiHat = phiAt[g] / VTn;
    s.u[bg * M] = phiHat;
    s.uLo[bg * M] = 0;
    for (let i = 0; i < n; i++) {
      const k = g * n + i;
      if (!s.present[k]) continue;
      if (Number.isNaN(target[k])) target[k] = freshC[k]; // (its own start, as a concentration, at the carried φ)
      s.uLo[bg * M + 1 + i] = 0;
      if (target[k] > 0) s.u[bg * M + 1 + i] = Math.log(target[k] / s.cRef[k]) + s.mu0hat[k] + s.z[i] * phiHat;
      // (nothing there, an underflow: its carried level, not its start, which refilled it from c0)
      else s.u[bg * M + 1 + i] = Number.isFinite(levelAt[k]) ? levelAt[k] / VTn : fresh[bg * M + 1 + i];
    }
    if (!s.nodeIdeal[g]) {
      // (an insertion host's carrier: its level, linear, the host's content its ion's)
      for (const md of s.model.materials[s.nodeMaterial[g]].models) {
        const k = g * n + md.idx[1];
        if (md.type === 'insertion' && Number.isFinite(levelAt[k])) s.u[bg * M + 1 + md.idx[1]] = levelAt[k] / VTn;
      }
      s._invertAt(g, target);
    }
  }

  // --- surfaces: by owner, each species by name, linear within the window; moles kept
  // (a face by its two regions, in the old solver's numbering)
  const faceKey = (sv, f) => (sv === old ? f : map[f] >= 0 && map[f + 1] === map[f] + 1 ? map[f] : `new${f}`);
  const surfKey = (sv, sf) => JSON.stringify(sf.port >= 0 ? ['port', sf.name] : ['face', faceKey(sv, sf.face)]);
  const oldSurf = new Map(old.model.surfaces.map((sf, j) => [surfKey(old, sf), j]));
  // (a port's window as declared, in x; a face's surface, its node)
  const windowOf = (sv, sf) => {
    if (sf.port < 0) return [sv.model.grid.x[sf.nodes[0]], sv.model.grid.x[sf.nodes[0]]];
    const p = sv.model.ports[sf.port], x0 = sv.model.grid.x[sv.model.grid.regionStart[p.region]];
    return [x0 + p.from, x0 + p.to];
  };
  // (a new node's place in the old device, through its region)
  const toOld = (g) => {
    const r = ng.nodeRegion[g], qr = map[r];
    if (qr < 0) return NaN;
    const L = s.model.regions[r].length, Lo = old.model.regions[qr].length;
    return og.x[og.regionStart[qr]] + (L > 0 ? ((ng.x[g] - ng.x[ng.regionStart[r]]) / L) * Lo : 0);
  };
  const surfCarried = new Set(), surfNodes = new Set(); // (surfaces with anything carried, and the nodes carried)
  let surfExact = true;
  s.model.surfaces.forEach((sf, j) => {
    const jo = oldSurf.get(surfKey(s, sf));
    if (jo === undefined) return;
    const of = old.model.surfaces[jo], capRatio = of.species[0].capacity / sf.species[0].capacity;
    const qOf = sf.species.map((sp) => of.species.findIndex((o) => o.name === sp.name));
    sf.nodes.forEach((g) => {
      const loc = windowLocate(og, of.nodes, windowOf(old, of), ng, windowOf(s, sf), g, toOld);
      if (!loc) return;
      const { a, b, wx } = loc;
      const bg = s.blockOfNode[g], th = new Float64Array(sf.species.length);
      if (copied[g] && a === b && sameT && sf.species.length === of.species.length && sf.species.every((sp, qq) => qOf[qq] === qq && sp.mu0 === of.species[qq].mu0 && sp.capacity === of.species[qq].capacity)) {
        const ba = old.blockOfNode[a];
        for (let qq = 0; qq < sf.species.length; qq++) (s.u[bg * M + 1 + n + qq] = old.u[ba * oldM + 1 + n + qq]), (s.uLo[bg * M + 1 + n + qq] = old.uLo[ba * oldM + 1 + n + qq]);
        surfCarried.add(j);
        surfNodes.add(g);
        return;
      }
      surfExact = false;
      let sum = 0;
      sf.species.forEach((sp, qq) => {
        const qo = qOf[qq];
        th[qq] = qo < 0 ? freshTh[g * s.nSurf + qq] : ((1 - wx) * old.th[a * old.nSurf + qo] + wx * old.th[b * old.nSurf + qo]) * capRatio;
        sum += th[qq];
      });
      // (more than the new sites hold: refused, the device left as it was; kept as θ instead, the
      // adsorbed amount would change with nothing to carry it)
      // (the bare fraction carried directly where nothing scales: a coverage within round-off of 1
      // made 1 − Σθ zero, and refused a patch that changed nothing there; Astra's second pass)
      const allMatched = qOf.every((qo) => qo >= 0) && sf.species.length === of.species.length;
      let bare = capRatio === 1 && allMatched ? (1 - wx) * old.th0[a] + wx * old.th0[b] : 1 - sum;
      if (!(bare > 0) && !(sum < 1)) throw new DeviceError(`set(): surface '${sf.name}' now has ${sf.species[0].capacity} mol/m² of sites, fewer than it holds (its coverage there would be ${sum.toPrecision(4)}); keep the capacity above what's adsorbed, or build a new Device to start over`);
      if (!(bare > 0)) bare = 1 - sum;
      sf.species.forEach((sp, qq) => {
        s.u[bg * M + 1 + n + qq] = th[qq] > 0 ? sp.mu0 / s.model.RT + Math.log(th[qq] / bare) : fresh[bg * M + 1 + n + qq];
        s.uLo[bg * M + 1 + n + qq] = 0;
      });
      surfCarried.add(j);
      surfNodes.add(g);
    });
  });

  // --- gates: by owner and name
  const gateKey = (sv, gt) => JSON.stringify(gt.port === undefined ? ['face', faceKey(sv, gt.f), gt.gate.name] : ['port', sv.model.ports[gt.port].name, gt.gate.name]);
  const oldGates = new Map();
  old.gateList.forEach((gt) => {
    const k = gateKey(old, gt);
    if (!oldGates.has(k)) oldGates.set(k, []);
    oldGates.get(k).push(gt);
  });
  const gatesAtRest = []; // (gates with nothing to carry: at rest at the carried voltage, below)
  for (const gt of s.gateList) {
    const list = oldGates.get(gateKey(s, gt));
    if (!list) {
      gatesAtRest.push(gt);
      continue;
    }
    if (gt.port === undefined) {
      s.u[gt.o] = old.u[list[0].o];
      s.uLo[gt.o] = old.uLo[list[0].o];
      continue;
    }
    const at = new Map(list.map((o) => [o.g, o.o]));
    const po = old.model.ports.findIndex((q) => q.name === s.model.ports[gt.port].name);
    const loc = windowLocate(og, [...at.keys()], windowOf(old, { port: po }), ng, windowOf(s, { port: gt.port }), gt.g, toOld);
    if (!loc) {
      gatesAtRest.push(gt);
      continue;
    }
    const { a, b, wx } = loc;
    if (a === b) (s.u[gt.o] = old.u[at.get(a)]), (s.uLo[gt.o] = old.uLo[at.get(a)]);
    else (s.u[gt.o] = (1 - wx) * (old.u[at.get(a)] + old.uLo[at.get(a)]) + wx * (old.u[at.get(b)] + old.uLo[at.get(b)])), (s.uLo[gt.o] = 0);
  }

  // --- faces: by their two regions; each slot by meaning
  const r2 = VTo / VTn;
  s.model.interfaces.forEach((itf, f) => {
    const qf = map[f];
    if (qf < 0 || map[f + 1] !== qf + 1) return;
    const oitf = old.model.interfaces[qf], bf = s.blockOfFace[f], bo = old.blockOfFace[qf];
    // (a slot copied whole, or scaled: its low part then folded in)
    const put = (k, ko, scale = 1) => {
      // (an unknown in both: a slot the new face doesn't solve keeps its own start, as a face made
      // neutral kept its old D for good, a sheet of charge; Astra's second pass)
      if (!s.active[bf * M + k] || !old.active[bo * oldM + ko]) return;
      if (scale === 1) (s.u[bf * M + k] = old.u[bo * oldM + ko]), (s.uLo[bf * M + k] = old.uLo[bo * oldM + ko]);
      else (s.u[bf * M + k] = (old.u[bo * oldM + ko] + old.uLo[bo * oldM + ko]) * scale), (s.uLo[bf * M + k] = 0);
    };
    put(0, 0); // (its displacement, C/m²)
    for (let i = 0; i < n; i++) put(1 + i, 1 + sp[i]); // (each species' flux through it, by name)
    const ok = reactionKeys(oitf, old.model.species), nk = reactionKeys(itf, species);
    nk.forEach((key, k) => {
      const ko = ok.indexOf(key);
      if (ko >= 0) put(1 + n + k, 1 + n + ko);
    });
    // (its gates were carried above, by name)
    for (const d of s.diffList) {
      if (d.f !== f) continue;
      const od = old.diffList.find((e) => e.f === qf && e.side === d.side);
      if (!od) continue;
      put(d.slot, od.slot, r2); // (the drop, in volts)
      d.drawSlot.forEach((slot, i) => {
        if (slot >= 0 && od.drawSlot[sp[i]] >= 0) put(slot, od.drawSlot[sp[i]], r2); // (a level at the plane, μ̄)
      });
    }
  });
  s.computeConcentrations();
  const info = { map, kindChanged, copied, moved, surfCarried, surfNodes, surfExact, sameGeometry, sameOrder, sp, faceKey, surfKey, gatesAtRest, levelAt };
  carryRest(old, s, info, { pure });
  return info;
}

// The clock, the terminals, neutrality, and what each stretch and surface conserves from here.
function carryRest(old, s, info, { pure }) {
  const { map, copied, sameGeometry } = info, n = s.n, ng = s.model.grid, og = old.model.grid;
  // (every new node copied, and no old one gone: a coarsening's surviving nodes are all old ones, and
  // its amounts still change; Astra's second pass)
  const allCopied = copied.every((v) => v === 1) && info.surfExact && og.nNodes === ng.nNodes && info.sameOrder;
  // The clock carries on (an advance(t) after a change at t still reaches t), held sources read
  // then; an adaptive transient goes on at its step size; carried as cold, as stepped, as solved.
  s.time = old.time;
  s.sourceTime = old.sourceTime;
  s.dtNext = old.dtNext;
  s.coldStart = old.coldStart;
  s.stepped = old.stepped;
  s.solvedV = old.solvedV;
  // (a jump in T, for a transient's first step to name if it fails: see integrate in src/timeSteps.js)
  const T0 = old.jumpedT?.[0] ?? old.model.T;
  s.jumpedT = T0 !== s.model.T ? [T0, s.model.T] : undefined;
  s.regriddedSteady = !!old.atSteady && pure && !allCopied;
  // Terminals by name: a floating one's voltage, each one's running charge; a contact's or port's
  // charge at the end of the last solve or step (a contact's per area: scaled by the contact's
  // area where the geometry changed; a port's per area of electrode, node by node).
  // (a contact of no area, a sphere's centre, holds no charge to scale: NaN there poisoned the state)
  const ratio = (k) => {
    const a = k === 0 ? og.area[0] : og.area[og.nNodes - 1], b = k === 0 ? ng.area[0] : ng.area[ng.nNodes - 1];
    return a > 0 ? b / a : 1;
  };
  s.terms.forEach((t, k) => {
    const o = old.terms.findIndex((u) => u.name === t.name);
    if (o >= 0 && s.floating.includes(k)) s.termV[k] = old.termV[o];
  });
  carryTermCharges(old, s);
  s._refreshSources();
  // What each stretch holds from here (stretchTargets), and the strictly neutral nodes whose every
  // ion's amount the carry keeps (rescale): there the ions' net charge is kept, and with it what
  // the capacitances balancing it hold (see the neutrality below).
  const targets = stretchTargets(old, s, info);
  // The strictly neutral boxes (src/neutral.js), and the groups capacitors join over them
  // (src/regions.js, _gateGroupsOf): a port's window, a capacitive contact's end box, a face.
  const nodes = strictlyNeutralNodes(s), neutralAt = new Uint8Array(s.nNodes);
  for (const nd of nodes) if (!nd.bulk) neutralAt[nd.g] = 1;
  const groupOfPort = (k) => s.gateGroups.find((gr) => gr.ports.includes(k));
  const groupOfContact = (side) => s.gateGroups.find((gr) => gr.contacts.includes(side));
  const sides = ['left', 'right'], endOf = (side) => (side === 'left' ? 0 : s.nNodes - 1), kOf = (side) => (side === 'left' ? 0 : 1);
  info.hardware = new Set();
  info.contactHw = new Set();
  info.heldGates = [];
  info.heldContacts = [];
  if (old.contactDEnd || old.portQEnd) {
    s._assembleBookkeeping(Infinity); // (each from the carried state, for any not carried)
    // (a contact beside a region that restarted holds its own start's charge: carried by name, the
    // plate's old charge beside a fresh region passed as a current of nothing; Astra's second pass)
    // Each terminal's charge carries what it held, unless the terminal itself changed (its φ law,
    // capacitance, window, area) or its region restarted: then it's what its drive puts on the
    // carried state, read once that state is final (a window widened under a clamp, charged by it).
    // (Read off an interpolated state, a held plate's charge was off by its φ's interpolation, and
    // the next step read that as a current.)
    const lastR = s.model.regions.length - 1, freshAt = (r) => map[r] < 0 || info.kindChanged.has(r);
    const contactSame = (side) => deepEqual(old.model.contacts[side].phi, s.model.contacts[side].phi) && old.model.contacts[side].zeroCharge === s.model.contacts[side].zeroCharge;
    // (one driven by a current, or behind a resistance, keeps its charge through any change of its
    // own: it can't jump, its voltage does)
    const driven = (k) => s.floating.includes(k);
    const capacitive = (m, side) => m.contacts[side].phi.type === 'capacitive';
    // A capacitive contact beside a strictly neutral end box is a gate over that box. Changed (its
    // C or zeroCharge, or newly capacitive), its plate's charge can't jump (a new plate's is none)
    // and the box's φ̂ does, at its concentrations (see the neutrality below). Unchanged, its
    // charge is what the final state holds at its old drive, as an unchanged gate's. (Kept from
    // before a regrid, the first 1 ns step's currents summed to −1.8e4 A/m² of 1.2 A/m², unwarned;
    // its C doubled in a sealed cell, φ̂ held, moved 1e-4 of the sealed ions and settled 50% off a
    // fresh device's charge, unwarned; Neon's fourth review, R4 1.1.)
    const own = {};
    for (const side of sides) {
      const k = kOf(side), r = k === 0 ? 0 : lastR, beside = neutralAt[endOf(side)] && capacitive(s.model, side);
      if (!old.contactDEnd || freshAt(r)) own[side] = true;
      else if (driven(k)) own[side] = false;
      else if (beside && !contactSame(side)) {
        own[side] = false;
        info.contactHw.add(side);
      } else if (beside) {
        own[side] = false;
        info.heldContacts.push(side);
      } else own[side] = !contactSame(side);
    }
    const ratioOf = (side) => ratio(kOf(side));
    s.contactDEnd = Object.fromEntries(sides.map((side) => [side, own[side] ? s.contactD[side] : info.contactHw.has(side) && !capacitive(old.model, side) ? 0 : old.contactDEnd[side] * ratioOf(side)]));
    info.ownCharge = { left: own.left, right: own.right, ports: s.model.ports.map(() => true) };
    s.portQEnd = Float64Array.from(s.portQ);
    // (each port's by name, a total: kept only in one geometry, a charge read at the new C and the
    // old V stayed for good; Astra's second pass)
    info.perNode = [];
    // A gate over strictly neutral boxes (src/regions.js, _gateGroupsOf) whose hardware changed
    // (its capacitance, window or area; or new): the concentrations under it can't move in an
    // instant, nor can its charge, so it keeps what it held over what's left of its window (a new
    // one, or a window's new part, nothing) as its start, and φ̂ under it jumps until each box is
    // neutral again at those concentrations (see the neutrality below). The first step then reads
    // the charge passing between its gates as currents that sum to zero. (Read off the state
    // instead, a gate added beside another took a charge its ions didn't move for: Kirchhoff broken
    // by 1.9e5 A/m², unwarned; Neon's third review, R3 1.5. Its whole old charge kept, a window
    // narrowed by half broke it by 8.6e7 A/m²; Neon's fourth review, R4 1.3.) Its unchanged
    // neighbours in the group start from what the state holds at their old drives before the jump.
    // Over a region with ε > 0 (no group: φ there is Poisson's, carried, and nothing jumps) a held
    // gate's charge can't jump either: it starts from what it held, and the first step passes the
    // difference through its terminal, recorded. (Read off the carried state at the new C, a
    // thin-film transistor's gate C doubled: the first 22 ms step's currents summed to −31 A/m² of
    // 31, both contacts passing it out, unwarned; bench/combos' kirchhoff, its first run.)
    const allCopied = (k) => s.model.ports[k].nodes.every((g) => info.copied[g] === 1);
    if (old.portQEnd) s.model.ports.forEach((p, k) => {
      const o = old.model.ports.findIndex((q) => q.name === p.name);
      // (in the same geometry, the charge it held: a window widened over the same nodes adds
      // uncharged membrane, its voltage moving if it floats; a new geometry keeps it per area of
      // electrode, node by node: see portChargeAcross)
      const op = old.model.ports[o];
      const same = o >= 0 && !freshAt(p.region) && (driven(2 + k) || (deepEqual(op.capacitance, p.capacitance) && op.from === p.from && op.to === p.to && deepEqual(op.area, p.area)));
      if (!same && p.capacitance && !freshAt(p.region) && !driven(2 + k)) {
        info.ownCharge.ports[k] = false;
        s.portQEnd[k] = o >= 0 && op.capacitance ? portChargeKept(old, s, o, k, info) : 0;
        info.hardware.add(k);
        return;
      }
      if (!same) return;
      // (an unchanged held gate over strictly neutral boxes: what it holds at its old drives, read
      // off the final state, whose φ̂ neutrality sets there; kept as a total, the new grid's sum of
      // its σ was 5e-5 of it off, read as a current)
      info.ownCharge.ports[k] = false;
      if (!driven(2 + k) && p.capacitance && groupOfPort(k)) info.heldGates.push([o, k, allCopied(k)]);
      // (a capacitance taken away or added under a current drive holds none to spread: its old total, 0)
      if (sameGeometry || !op.capacitance || !p.capacitance) s.portQEnd[k] = old.portQEnd[o] * (sameGeometry ? 1 : portRatio(old, s, o, k));
      else (s.portQEnd[k] = portChargeAcross(old, s, o, k)), info.perNode.push([o, k]);
    });
    // (beside a changed gate or plate, an unchanged one starts from what the state holds at its old
    // drives before φ̂ jumps under them both: the jump moves it)
    const hwGroup = (gr) => !!gr && (gr.ports.some((q) => info.hardware.has(q)) || gr.contacts.some((sd) => info.contactHw.has(sd)));
    info.heldGates = info.heldGates.filter(([o, k]) => {
      if (!hwGroup(groupOfPort(k))) return true;
      s.portQEnd[k] = portChargeAcross(old, s, o, k);
      return false;
    });
    info.heldContacts = info.heldContacts.filter((side) => {
      if (!hwGroup(groupOfContact(side))) return true;
      s.contactDEnd[side] = contactChargeAcross(old, s, side, ratioOf(side));
      return false;
    });
    info.hwGroups = s.gateGroups.filter(hwGroup);
    const listed = (k) => info.heldGates.some(([, q]) => q === k);
    info.perNode = info.perNode.filter(([, k]) => !info.hwGroups.some((gr) => gr.ports.includes(k)) && !listed(k));
  }

  // Neutrality (src/neutral.js, solveNeutral): each strictly neutral node onto its Gauss row,
  // G = 0, as assembled from the carried state with the new definition's ports, faces and
  // contacts, together with each closed stretch's amount (exactly, through the statistics) and
  // each floating capacitor's charge (by its voltage), in one Newton solve. Each node by its φ̂ at
  // the levels carried (the levels are what an open solution has; φ̂ is what neutrality makes of
  // them, a face's D and a held gate's σ with it, their start charges read off the final state);
  // under a gate or plate whose hardware changed (and where a removed gate's window was), at its
  // concentrations instead. Then each double layer's charge, again, as each moves the others a
  // little. (Held to each node's old net charge, not its Gauss row, a capacitive face's
  // counter-charge was lost, Neon's R2 3.1; with φ̂ held instead, a regridded gate's charge kept the
  // coarse grid's error, 8e-6 of it, where by its levels the state is the fine grid's to 1e-8. And
  // by turns, each projection undoing another's, a trace open ion in closed salt had no fixed point
  // and a floating gate beside a changed one stalled 89% short; Neon's fourth review, R4 1.3, 1.5.)
  // The held terminals at their old voltages meanwhile: a drive changed in the same set() is the
  // next step's to pass, as a drive changed alone is.
  const moved = info.moved;
  const hw = new Uint8Array(s.nNodes);
  for (const k of info.hardware) for (const g of s.model.ports[k].nodes) hw[g] = 1;
  for (const side of info.contactHw) hw[endOf(side)] = 1;
  // (every box of a group a changed capacitor is in: its unchanged neighbours' start charges were
  // read before the jump, so nothing else may move φ̂ under them)
  for (const gr of info.hwGroups ?? []) for (const g of gr.nodes) hw[g] = 1;
  if (info.hardware.size > 0) old.model.ports.forEach((op, o) => {
    if (!op.capacitance || s.model.ports.some((p) => p.name === op.name)) return;
    const r = map.indexOf(op.region);
    if (r < 0) return;
    const x0 = ng.x[ng.regionStart[r]], L = s.model.regions[r].length, Lo = old.model.regions[op.region].length;
    for (let g = ng.regionStart[r]; g <= ng.regionEnd[r]; g++) {
      const xo = L > 0 ? ((ng.x[g] - x0) / L) * Lo : 0;
      if (xo >= op.from && xo <= op.to) hw[g] = 1;
    }
  });
  for (const nd of nodes) nd.mode = hw[nd.g] && !nd.bulk ? 'jump' : 'level';
  const V0 = Float64Array.from(s.termV), heldV = [];
  s.terms.forEach((t, k) => {
    const o = old.terms.findIndex((u) => u.name === t.name);
    if (o >= 0 && t.drive.kind === 'V' && !s.floating.includes(k)) (s.termV[k] = old.termV[o]), heldV.push(k);
  });
  const closedAll = targets.filter((t) => t.rescale && t.want > 0);
  // On a neutral lattice with every partner's amount carried, their shared-site derivatives
  // belong in the same solve as Gauss. In turns the amounts undid its correction: 2.13e-9 of
  // the ions after 30 rounds, still 2.13e-9 after 150. With an open partner, rescale in
  // concentrations to keep that partner; on an insertion host the carrier follows its ion
  // (the echem combos' case 29, neutralG; Astra, 25 Nov).
  const joint = [], apart = [];
  for (const target of closedAll) {
    const st = target.st;
    let separate = false;
    for (let r = st.regions[0]; r <= st.regions[1]; r++) {
      const mat = s.model.materials[s.model.regions[r].material], md = mat.models[mat.modelOf[st.species]];
      if (md?.type === 'insertion') separate = true;
      if (md?.type === 'lattice' && (!pure || mat.epsr !== 0 || md.idx.some((i) => {
        const t = targets[s.stretchOf[r * s.n + i]];
        return !t.rescale || !(t.want > 0);
      }))) separate = true;
    }
    (separate ? apart : joint).push(target);
  }
  const floating = floatingTargets(s, neutralAt);
  const layers = diffuseTargets(old, s, info);
  for (let round = 0; round < 30; round++) {
    let { worst: off } = solveNeutral(s, nodes, { closed: joint, floating, moved, syncFaces: true });
    off = Math.max(off, holdDiffuse(s, layers, moved));
    off = Math.max(off, rescaleAll(s, apart, moved));
    // (a floating port's voltage again, as neutrality moved the φ under it: held once, before, a
    // port over a bar made 1 → 3 kept 9e-9 C/m² off what it held, and the first step drove that
    // through the bar, its voltage jumping 0.58 mV in 1 ns, 58 mV in 10 ps; the inventory's item 1)
    off = Math.max(off, holdFloatingCharges(s, moved));
    if (round > 0 && off < 1e-13) break;
  }
  const surfaces = surfaceTargets(old, s, info);
  rescaleSurfaces(s, surfaces);
  s.computeConcentrations();
  // A gate with nothing to carry (new, or on a window's new part) starts at rest at the voltage
  // the carried state puts across it, α/(α + β), as a start does at its own (from the cold
  // start's, a widened window's new membrane had its gates open at 0 mV beside −65; Astra's review).
  for (const k of heldV) s.termV[k] = V0[k];
  for (const gt of info.gatesAtRest) {
    const V = s._gateVoltage(gt), a = gateRate(gt.gate.alpha, V)[0], b = gateRate(gt.gate.beta, V)[0];
    s.u[gt.o] = a + b > 0 ? a / (a + b) : 0;
    s.uLo[gt.o] = 0;
  }
  // The charges held terminals take from the state, read now it's final (read before the amounts
  // and layers were put right, a held plate's step read the difference as a current: 6.3e-3 where
  // 3.0e-4 is right).
  if (info.ownCharge) {
    s._assembleBookkeeping(Infinity);
    if (info.ownCharge.left) s.contactDEnd.left = s.contactD.left;
    if (info.ownCharge.right) s.contactDEnd.right = s.contactD.right;
    info.ownCharge.ports.forEach((own, k) => own && (s.portQEnd[k] = s.portQ[k]));
    // (a held port's charge across a new geometry, from the final state too: its φ moved since; a
    // floating one's voltage was set to hold what it held, above, and stays. Where nothing under
    // it moved, the charge it held, bit for bit.)
    const still = (ns) => ns.every((g) => !moved[g]);
    for (const [o, k] of info.perNode) if (!s.floating.includes(2 + k)) s.portQEnd[k] = portChargeAcross(old, s, o, k);
    for (const [o, k, copied] of info.heldGates) if (!(copied && sameGeometry && still(s.model.ports[k].nodes))) s.portQEnd[k] = portChargeAcross(old, s, o, k);
    for (const side of info.heldContacts) if (!(sameGeometry && info.copied[endOf(side)] === 1 && still([endOf(side)]))) s.contactDEnd[side] = contactChargeAcross(old, s, side, ratio(kOf(side)));
    // (a group with a changed capacitor in a regrid: its start charges together are what the new
    // grid's ions under it balance, the interpolation's difference shared out by capacitance, as
    // φ̂ would share it. Kept as they were, a C doubled with a regrid broke Kirchhoff by 4.1e4 A/m²;
    // Neon's fourth review, R4 1.3.)
    // A group with a changed capacitor at a capacitive face to a metal, on the same nodes: the
    // jump moved the face's D, and the metal's sheet with it, which can't move in an instant any
    // more than a gate's charge can; so the sheet starts from the D before the jump (sheetDEnd),
    // as the changed gate starts from its old charge, and the first step passes what moved
    // between them as currents through their terminals, which record it. (Shared out onto the
    // held sources' starts instead, as an impulse at the set(), Kirchhoff held but no terminal's
    // running charge saw 0.0115 C/m² pass between a gate and its Pt; Neon's fourth review, its
    // closing reply.) A start still to pass from an earlier set() carries by its difference from
    // the state, as a port's (portChargeAcross). Elsewhere (a regrid, a face to a non-metal), the
    // share.
    const nF = s.model.interfaces.length, sameFaces = old.model.interfaces.length === nF;
    const faceD = (sv, f) => sv.u[sv.blockOfFace[f] * sv.M];
    const sheet = (f) => sameFaces && s.sheetFace.includes(f) && old.sheetFace.includes(f);
    const pending = (f) => (sameFaces && Number.isFinite(old.sheetDEnd?.[f]) ? old.sheetDEnd[f] - faceD(old, f) : 0);
    let sheetD = null;
    const put = (f, D) => ((sheetD ??= new Float64Array(nF).fill(NaN))[f] = D);
    for (let f = 0; f < nF; f++) if (sheet(f) && pending(f) !== 0) put(f, faceD(s, f) + pending(f));
    for (const gr of info.hwGroups ?? []) {
      const copied = gr.nodes.every((g) => info.copied[g] === 1);
      if (copied && gr.faces.every(({ f, side }) => side === 2 || sheet(f))) {
        for (const { f, side } of gr.faces) if (side !== 2) put(f, faceD(old, f) + pending(f));
      } else if (gr.faces.length > 0 || !copied) {
        if (sheetD) for (const { f } of gr.faces) sheetD[f] = NaN;
        shareGroupCharge(s, gr);
      }
    }
    s.sheetDEnd = sheetD;
  }
  // (a carried state that isn't finite is refused, the device left as it was: a crowded double layer
  // whose ε fell twentyfold came back NaN and converged, unwarned; Astra's second pass)
  for (let k = 0; k < s.u.length; k++) {
    if (!Number.isFinite(s.u[k]) || !Number.isFinite(s.uLo[k])) {
      const b = Math.floor(k / s.M), g = s.blockOfNode.indexOf(b);
      throw new DeviceError(`set(): the state couldn't be carried onto the new definition${g >= 0 ? ` (at x = ${s.model.grid.x[g]})` : ' (at a face)'}: it came out undefined. Make the change in smaller steps, or solve() at the new definition from a new Device`);
    }
  }
  // What each stretch conserves from here (see the contract, 5).
  s.referenceAmounts = s.stretches.map((st, k) => {
    const t = targets[k];
    if (t.keep !== undefined) return t.keep;
    return !st.connected && Number.isFinite(t.want) ? t.want : s.amount(st);
  });
  s.boundaryIntake.fill(0);
  s.surfaceRef = Float64Array.from(s.surfCols, (_, j) => (Number.isFinite(surfaces[j].want) ? surfaces[j].want : s.surfaceAmount(j)));
  // A region that had a state and got none: said in the next solutions.
  const fresh = s.model.regions.flatMap((rg, r) => (map[r] < 0 || info.kindChanged.has(r) ? [rg.name] : []));
  if (fresh.length > 0) {
    const all = fresh.length === s.model.regions.length, many = fresh.length > 1 || all;
    s.restarted = `set(): ${all ? 'no region' : fresh.map((nm) => `'${nm}'`).join(', ')} had ${all ? '' : many ? 'no region' : 'no region'} before ${all ? 'it' : many ? 'them' : 'it'} to carry a state from (matched by name; or, as many regions as before, by place; a metal and a non-metal don't carry), so ${all ? 'the device' : many ? 'they' : 'it'} started from ${many ? 'their' : 'its'} c0 at t = ${old.time} (the clock carries on). (Said in every solution up to the next steady solve.)`;
  } else if (old.restarted) s.restarted = old.restarted;
}

// A port's electrode area where its capacitance holds charge (Σ box · area per volume over the
// nodes Solver._portHolds keeps), new over old: what its running charge scales by. (1 where either
// has none.)
function portArea(sv, k) {
  const port = sv.model.ports[k], a = sv.portArea?.[k];
  if (!port || !a) return 0;
  let t = 0;
  port.nodes.forEach((g, w) => sv._portHolds(g) && (t += sv.model.grid.vol[g] * a[w]));
  return t;
}
export function portRatio(old, s, o, k) {
  const a = portArea(old, o), b = portArea(s, k);
  return a > 0 && b > 0 ? b / a : 1;
}

// A port's stored charge across a new geometry: what it held per area of electrode, node by node,
// over the new electrode, read through Solver._portCapacitance's own sum (s.portQ, at the carried
// state and s's V, C and zeroCharge, as last assembled). σ is linear in V, so at the old port's
// C_o, zeroCharge z_o and V_o it's (C_o/C) s.portQ + C_o (V_o − z_o − V + z) A, A the new electrode
// (portArea), plus what the old one held beyond its state's own charge (a drive changed since its
// last step), spread evenly per area. (One ratio of areas for it all was exact only for a uniform σ:
// a bar's area made 1 → 3 carried 1.9% off, the first 1 ns step reading +5.59e5 A/m² where −1.08e-2
// is right; the inventory of duplicated quantities.)
function portChargeAcross(old, s, o, k) {
  const po = old.model.ports[o].capacitance, pn = s.model.ports[k].capacitance;
  const Ao = portArea(old, o), An = portArea(s, k);
  const held = Ao > 0 ? ((old.portQEnd[o] - old.portQ[o]) * An) / Ao : 0;
  return (po.C / pn.C) * s.portQ[k] + po.C * (old.termV[2 + o] - po.zeroCharge - s.termV[2 + k] + pn.zeroCharge) * An + held;
}

// A terminal's running charge, by name, scaled as its stored charge is (a contact's by its area, a
// port's by its electrode's; one rule for all, so the terminals' charges still add up to what the
// device holds; Astra's second pass). Used by the carry and by a restart.
export function carryTermCharges(old, s) {
  const og = old.model.grid, ng = s.model.grid, sameGeometry = deepEqual(old.model.geometry, s.model.geometry);
  const contact = (k) => {
    const a = k === 0 ? og.area[0] : og.area[og.nNodes - 1], b = k === 0 ? ng.area[0] : ng.area[ng.nNodes - 1];
    return a > 0 ? b / a : 1;
  };
  s.terms.forEach((t, k) => {
    const o = old.terms.findIndex((u) => u.name === t.name);
    if (o < 0) return;
    const r = k < 2 ? (o < 2 ? contact(k) : 1) : o >= 2 && !sameGeometry ? portRatio(old, s, o - 2, k - 2) : 1;
    s.termCharge[k] = old.termCharge[o] * r;
    s.termChargeLo[k] = old.termChargeLo[o] * r;
  });
}

// Each floating terminal at no current that stores charge (a capacitive contact's plate, a port's
// spread capacitance): its voltage set so the state holds the charge it ended with (linear in V).
// Returns how far off the worst was, relative to its charge or a µV's worth of it.
function holdFloatingCharges(s, moved) {
  if (!s.contactDEnd && !s.portQEnd) return 0;
  // (driven by a current and holding a capacitance's charge, which can't jump: a capacitive
  // contact's plate, a port's spread capacitance; ionic leaks beside it don't change that)
  const stores = (k) => (k >= 2 ? !!s.model.ports[k - 2].capacitance : s.model.contacts[k === 0 ? 'left' : 'right'].phi.type === 'capacitive');
  const held = s.floating.filter(stores);
  if (held.length === 0) return 0;
  const charge = (k) => (k >= 2 ? s.portQ[k - 2] : k === 0 ? s.contactD.left : -s.contactD.right);
  const want = (k) => (k >= 2 ? s.portQEnd?.[k - 2] : k === 0 ? s.contactDEnd?.left : -s.contactDEnd?.right);
  let worst = 0;
  for (const k of held) {
    const target = want(k);
    if (!Number.isFinite(target)) continue;
    s._assembleBookkeeping(Infinity);
    const q0 = charge(k), V0 = s.termV[k], dV = 1e-6 * Math.max(1, Math.abs(V0));
    s.termV[k] = V0 + dV;
    s._assembleBookkeeping(Infinity);
    const slope = (charge(k) - q0) / dV;
    // (against its charge, or a µV's worth of it where it holds none: a thermal volt's worth, 1e5
    // times a gate's charge, let the carry stop with 1.5e-9 of it unheld, a current at the first step)
    worst = Math.max(worst, Math.abs(target - q0) / (Math.abs(target) + Math.abs(slope) * 1e-6));
    s.termV[k] = Math.abs(slope) > 0 && Math.abs(target - q0) > 1e-14 * Math.abs(target) ? V0 + (target - q0) / slope : V0;
    // (a port's voltage moved: its spread charge moved in every box of its window, which a strictly
    // neutral region then balances)
    if (s.termV[k] !== V0 && k >= 2) for (const g of s.model.ports[k - 2].nodes) moved[g] = 1;
  }
  s._assembleBookkeeping(Infinity);
  return worst;
}

// A gate whose window changed (narrowed, moved, widened; its C or area too): what it held over the
// part of its new window its old one covered, per area of electrode, node by node, at the state as
// carried (before the neutrality moves φ̂: the old φ̂ there), and the old window's charge beyond its
// state's own (a drive changed since its last step) spread over it by electrode area. Where the
// window is the same, its whole charge (scaled with its electrode across a new geometry).
function portChargeKept(old, s, o, k, info) {
  const op = old.model.ports[o], p = s.model.ports[k], sameGeometry = deepEqual(old.model.geometry, s.model.geometry);
  if (op.from === p.from && op.to === p.to) return old.portQEnd[o] * (sameGeometry ? 1 : portRatio(old, s, o, k));
  if (op.region !== info.map[p.region]) return 0;
  const ng = s.model.grid, og = old.model.grid, po = op.capacitance, V = old.termV[2 + o], Ao = portArea(old, o);
  const sigmaOld = (sv, g) => po.C * (V - po.zeroCharge - sv.VT * (sv.u[sv.blockOfNode[g] * sv.M] + sv.uLo[sv.blockOfNode[g] * sv.M]));
  let q = 0, kept = 0;
  // (on the same nodes: each box both windows hold keeps its whole old charge, the box's φ̂ jumping
  // at its concentrations, which balance just that; a box the new one doesn't hold has none, its
  // ions moving)
  if (sameGeometry && p.nodes.every((g) => info.copied[g] === 1)) {
    const at = new Map();
    op.nodes.forEach((g, w) => old._portHolds(g) && at.set(og.x[g], og.vol[g] * old.portArea[o][w]));
    for (const g of p.nodes) {
      const e = at.get(ng.x[g]);
      if (e === undefined || !s._portHolds(g)) continue;
      const go = og.x.indexOf(ng.x[g]);
      q += e * sigmaOld(old, go);
      kept += e;
    }
  } else {
    // (across a regrid: per area of electrode over the part of each new box the old window covered,
    // at the φ̂ carried there)
    const x0 = ng.x[ng.regionStart[p.region]], L = s.model.regions[p.region].length, Lo = old.model.regions[op.region].length;
    const toNew = (xo) => x0 + (Lo > 0 ? (xo / Lo) * L : 0), a = toNew(op.from), b = toNew(op.to);
    p.nodes.forEach((g, w) => {
      if (!s._portHolds(g) || !(p.hi[w] > p.lo[w])) return;
      const part = Math.max(0, Math.min(p.hi[w], b) - Math.max(p.lo[w], a)) / (p.hi[w] - p.lo[w]);
      if (!(part > 0)) return;
      const e = ng.vol[g] * s.portArea[k][w] * part;
      q += e * sigmaOld(s, g);
      kept += e;
    });
  }
  return q + (Ao > 0 ? ((old.portQEnd[o] - old.portQ[o]) * kept) / Ao : 0);
}

// A held capacitive contact's charge on the carried state at its old drive (its plate's D_in =
// A C (V − zeroCharge − V_T φ̂), linear in V), plus what it held beyond its old state's own (a drive
// changed since its last step), as contactD (the displacement toward +x).
function contactChargeAcross(old, s, side, ratio) {
  const sg = side === 'left' ? 1 : -1, k = side === 'left' ? 0 : 1;
  const po = old.model.contacts[side], pn = s.model.contacts[side];
  const A = s.model.grid.area[side === 'left' ? 0 : s.nNodes - 1];
  const Din = sg * s.contactD[side];
  const at = (po.phi.C / pn.phi.C) * Din + A * po.phi.C * (old.termV[k] - po.phi.zeroCharge - s.termV[k] + pn.phi.zeroCharge);
  return sg * at + (old.contactDEnd[side] - old.contactD[side]) * ratio;
}

// The floating capacitors whose charge sits on strictly neutral boxes, for solveNeutral to hold by
// their voltages: a floating port over any of them, a floating capacitive contact beside one.
function floatingTargets(s, neutralAt) {
  if (!s.contactDEnd && !s.portQEnd) return [];
  const out = [];
  for (const k of s.floating) {
    if (k >= 2) {
      const p = s.model.ports[k - 2];
      if (p.capacitance && Number.isFinite(s.portQEnd?.[k - 2]) && p.nodes.some((g) => neutralAt[g] && s._portHolds(g))) out.push({ k, want: s.portQEnd[k - 2] });
    } else {
      const side = k === 0 ? 'left' : 'right', g = k === 0 ? 0 : s.nNodes - 1;
      if (s.model.contacts[side].phi.type === 'capacitive' && neutralAt[g] && Number.isFinite(s.contactDEnd?.[side])) out.push({ k, want: k === 0 ? s.contactDEnd.left : -s.contactDEnd.right });
    }
  }
  return out;
}

// A group's start charges brought onto what its state holds in all (Σ of its ports' and plates'),
// the difference shared by each one's capacitance over the group's boxes (a floating one's charge
// is its own: it keeps it). See the start charges in carryRest.
function shareGroupCharge(s, gr) {
  s._assembleBookkeeping(Infinity);
  const vol = s.model.grid.vol, grid = s.model.grid;
  let want = 0, have = 0, gross = 0, wsum = 0;
  const share = [];
  for (const k of gr.ports) {
    want += s.portQEnd[k];
    have += s.portQ[k];
    gross += Math.abs(s.portQ[k]);
    if (s.floating.includes(2 + k)) continue;
    let w = 0;
    s.model.ports[k].nodes.forEach((g, j) => s._portHolds(g) && (w += vol[g] * s.portArea[k][j] * s.model.ports[k].capacitance.C));
    share.push(['port', k, w]);
    wsum += w;
  }
  for (const side of gr.contacts) {
    const sg = side === 'left' ? 1 : -1;
    want += sg * s.contactDEnd[side];
    have += sg * s.contactD[side];
    gross += Math.abs(s.contactD[side]);
    if (s.floating.includes(side === 'left' ? 0 : 1)) continue;
    const w = grid.area[side === 'left' ? 0 : s.nNodes - 1] * s.model.contacts[side].phi.C;
    share.push(['contact', side, w]);
    wsum += w;
  }
  const diff = have - want;
  if (!(wsum > 0) || !(Math.abs(diff) > 1e-15 * gross)) return;
  for (const [kind, k, w] of share) {
    const d = (diff * w) / wsum;
    if (kind === 'port') s.portQEnd[k] += d;
    else s.contactDEnd[k] += (k === 'left' ? 1 : -1) * d;
  }
}

// What the old state held of species i over solver sv's region q (its faces' diffuse layers'
// share included), on its own grid.
function regionAmount(sv, q, i) {
  const gr = sv.model.grid, n = sv.n;
  let a = 0;
  for (let g = gr.regionStart[q]; g <= gr.regionEnd[q]; g++) a += gr.vol[g] * sv.c[g * n + i];
  sv.diffList.forEach((dl, k) => {
    if (dl.g >= gr.regionStart[q] && dl.g <= gr.regionEnd[q]) a += gr.area[dl.g] * sv.gam[k * n + i];
  });
  return a;
}

// A region's volume (Σ box), per area where planar.
function regionVolume(sv, q) {
  const gr = sv.model.grid;
  let v = 0;
  for (let g = gr.regionStart[q]; g <= gr.regionEnd[q]; g++) v += gr.vol[g];
  return v;
}

// What old region q held of species oi, carried onto new region r resized or in a new geometry,
// concentrations kept: its volume's part, Σ box·c, scaled node by node, by Σ box·c over Σ box̂·c at
// the carried c, box̂ a new box measured where it sits in the old region (its fraction along it) in
// the old cross-section; and each diffuse layer's, Γ per area, by its plane's area. (One ratio of
// volumes for both was exact only for a uniform c and no layer: a 10% longer KCl film beside two
// diffuse layers moved its bulk by 6.0e-3 relative at 1 mM, the layers' excess spread through it;
// a bar's area made 1 → 3 moved its concentrations 3e-5, the first step reading the φ that moved
// with them as a port's charging; the inventory of duplicated quantities.)
function resizedAmount(old, s, info, r, q, i, oi) {
  const og = old.model.grid, ng = s.model.grid, n = s.n, geo = old.model.geometry;
  const A = geo && geo.type !== 'planar' ? crossSection(geo) : null;
  const L = s.model.regions[r].length, Lo = old.model.regions[q].length;
  const x0 = ng.x[ng.regionStart[r]], xo0 = og.x[og.regionStart[q]], toOld = (x) => xo0 + (L > 0 ? ((x - x0) / L) * Lo : 0);
  let bulk = 0, held = 0, num = 0, den = 0;
  for (let g = og.regionStart[q]; g <= og.regionEnd[q]; g++) bulk += og.vol[g] * old.c[g * old.n + oi];
  for (let g = ng.regionStart[r]; g <= ng.regionEnd[r]; g++) {
    // (its box as applyGeometry draws it, src/grid.js, mapped into the old region)
    const lo = g > ng.regionStart[r] ? ng.x[g] - ng.segLength[g - 1] / 2 : ng.x[g], hi = g < ng.regionEnd[r] ? ng.x[g] + ng.segLength[g] / 2 : ng.x[g];
    const a = toOld(lo), b = toOld(hi), c = s.present[g * n + i] ? s.c[g * n + i] : 0;
    num += ng.vol[g] * c;
    den += (b > a ? (A ? A.integral(a, b) : b - a) : 0) * c;
  }
  const scale = num > 0 && den > 0 ? num / den : regionVolume(s, r) / regionVolume(old, q);
  old.diffList.forEach((od, k) => {
    if (od.g < og.regionStart[q] || od.g > og.regionEnd[q]) return;
    const d = s.diffList.find((e) => info.map[e.f] === od.f && info.map[e.f + 1] === od.f + 1 && e.side === od.side);
    const amt = og.area[od.g] * old.gam[k * old.n + oi];
    // (a layer gone, its excess spread through the region as before)
    held += d && og.area[od.g] > 0 ? amt * (ng.area[d.g] / og.area[od.g]) : amt * scale;
  });
  return bulk * scale + held;
}

// Each new stretch: its amount from here, want, region by region: what the old state held over a
// region matched to one that held the species (scaled by its volume where resized or reshaped), its own start's
// elsewhere (a fresh region, one newly holding it); whether that amount is to be written into the
// state (rescale: not where a held terminal sets its level, which a shift would move off); and its
// reference kept exactly (keep: closed on its own before and after, over the same space). (All or
// nothing, one fresh region dropped every matched one's amount; a resize adopted the new one;
// Astra's second pass.)
function stretchTargets(old, s, info) {
  const { map, sameGeometry, sp } = info;
  const held = (k) => s.terms[k] && s.terms[k].drive.kind === 'V' && !(s.terms[k].drive.R > 0);
  return s.stretches.map((st) => {
    const oi = sp[st.species];
    let want = 0, same = sameGeometry;
    const regs = [];
    for (let r = st.regions[0]; r <= st.regions[1]; r++) {
      const q = map[r];
      regs.push(q);
      const carried = oi >= 0 && q >= 0 && !info.kindChanged.has(r) && old.model.materials[old.model.regions[q].material].present[oi];
      // (a region resized or in a new geometry keeps its concentrations: its amount scales with
      // its volume; a planar device's amounts are per m², so moles kept would have meant nothing)
      const resized = carried && (!sameGeometry || old.model.regions[q].length !== s.model.regions[r].length);
      if (!carried) want += regionAmount(s, r, st.species);
      else if (!resized) want += regionAmount(old, q, oi) * (regionVolume(s, r) / regionVolume(old, q));
      else want += resizedAmount(old, s, info, r, q, st.species, oi);
      if (!carried || old.model.regions[q].length !== s.model.regions[r].length) same = false;
    }
    const fedHeld = (st.leftOpen && held(0)) || (st.rightOpen && held(1)) || st.linked.some((k) => held(2 + k));
    let keep;
    if (same && st.spectator && regs.every((q, j) => q === regs[0] + j)) {
      const k = old.stretches.findIndex((o) => o.species === oi && o.regions[0] === regs[0] && o.regions[1] === regs[regs.length - 1]);
      if (k >= 0 && old.stretches[k].spectator) keep = old.referenceAmounts[k] + old.boundaryIntake[k];
    }
    return { st, want, rescale: !fedHeld, keep };
  });
}

// Scale each stretch's species until it holds its `want`, all together: an ideal node's levels
// shifted (c ∝ e^η), a non-ideal node's concentrations put back with each species scaled by its
// stretch's factor and the others held (a level shifted on shared sites moved an open partner too;
// and one species at a time, an overshoot past a lattice's capacity refused a state that fits;
// Astra's second pass). Returns how far off the worst was (relative); round-off isn't worth moving.
export function rescaleAll(s, list, moved) {
  const { n, M } = s;
  let first = 0;
  for (let it = 0; it < 40; it++) {
    const f = list.map(({ st, want }) => {
      const have = s.amount(st);
      return want > 0 && have > 0 ? want / have : 1;
    });
    const off = f.reduce((m, x) => Math.max(m, Math.abs(x - 1)), 0);
    if (it === 0) first = off;
    if (off < (it === 0 ? 1e-13 : 1e-14)) break;
    const tgt = Float64Array.from(s.c), inv = new Set();
    list.forEach(({ st }, j) => {
      if (f[j] === 1) return;
      const i = st.species;
      for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
        if (s.nodeConductor[g] >= 0 || !s.present[g * n + i]) continue;
        moved[g] = 1;
        if (s.nodeIdeal[g]) s.u[s.blockOfNode[g] * M + 1 + i] += Math.log(f[j]);
        else {
          // (an insertion host's carrier follows its ion, c = ν s + bg: its amount is moved through
          // the ion's; held, a regridded host's closed electrons drifted 2.6e-4)
          const md = s.model.materials[s.nodeMaterial[g]].models.find((m) => m.type === 'insertion' && m.idx[1] === i);
          if (md) {
            const ion = md.idx[0], ce = s.c[g * n + i], ci = s.c[g * n + ion], bg = s.background[g], nu = ci > 0 ? (ce - bg) / ci : 0;
            if (nu > 0) tgt[g * n + ion] = (f[j] * ce - bg) / nu;
          } else {
            // On a lattice, a species scaled up past (nearly) the room the others leave it has its
            // odds scaled instead, c/(room − c) by the factor: never past the sites, and the rest of
            // the stretch takes up what that node can't over the next rounds. (Multiplied, a
            // near-full lattice refined was refused as over its capacity though the amount fitted:
            // 4.2302 asked of 4.2289 sites where a fresh solve on the finer grid held 4.2112; Astra's
            // ninth pass. Odds everywhere moved a node near full too slowly to converge in 40 rounds.)
            const lat = s.model.materials[s.nodeMaterial[g]].models.find((m) => m.type === 'lattice' && m.idx.includes(i));
            const ci = s.c[g * n + i];
            let room = lat?.cMax;
            if (lat) for (const b of lat.idx) if (b !== i) room -= s.c[g * n + b];
            if (lat && f[j] > 1 && ci < room && ci * f[j] > room * (1 - 1e-3)) {
              const odds = (f[j] * ci) / (room - ci);
              tgt[g * n + i] = (room * odds) / (1 + odds);
            } else tgt[g * n + i] *= f[j];
          }
          inv.add(g);
        }
      }
    });
    try {
      for (const g of inv) s._invertAt(g, tgt);
    } catch (err) {
      throw new DeviceError(`set(): the carried amounts don't fit the new definition (${err.message}); keep the capacity above what's held, or build a new Device to start over`);
    }
    s.computeConcentrations();
  }
  return first;
}

// Each diffuse layer carried (its face matched, its side there before): the charge it held, which
// it keeps (charges don't move in an instant; its potential does).
function diffuseTargets(old, s, info) {
  const { map } = info;
  return s.diffList.flatMap((d) => {
    const qf = map[d.f];
    if (qf < 0 || map[d.f + 1] !== qf + 1) return [];
    const od = old.diffList.find((e) => e.f === qf && e.side === d.side);
    return od ? [{ d, want: od.dl.q }] : [];
  });
}

// Each carried layer's drop set so it holds its charge: bracketed, then secant within the bracket
// (the charge rises with the drop). Returns how far off the worst was (relative). (Held once, the
// amounts' shift undid it: the first step read 4.7e7 A/m²; and a secant alone stopped short or
// went to NaN on a crowded layer; Astra's second pass.)
function holdDiffuse(s, layers, moved) {
  let worst = 0;
  for (const { d, want } of layers) {
    s.computeConcentrations();
    const q = () => d.dl.q;
    let x = s.u[d.o], f = q() - want;
    const tol = 1e-14 * Math.max(Math.abs(want), 1e-300);
    if (!(Math.abs(f) > tol)) continue;
    worst = Math.max(worst, Math.abs(f) / Math.max(Math.abs(want), 1e-300));
    moved[d.g] = 1; // (its excesses moved: the amounts it shares are put back)
    const at = (v) => {
      s.u[d.o] = v;
      s.computeConcentrations();
      return q() - want;
    };
    // bracket: step away from x on both sides, doubling, until the charge passes what it held
    // (which way it rises depends on the layer's sign)
    let lo = x, flo = f, hi = x, fhi = f, found = false;
    for (let it = 0, step = Math.max(1e-6, Math.abs(x) * 1e-6); it < 80 && !found; it++, step *= 2) {
      for (const dir of [1, -1]) {
        const v = x + dir * step, fv = at(v);
        if (Number.isFinite(fv) && Math.sign(fv) !== Math.sign(f)) {
          (hi = v), (fhi = fv), (found = true);
          // (the nearer side of the bracket: the last point on this side with f's sign)
          const back = x + dir * (step / 2), fb = it > 0 ? at(back) : f;
          if (it > 0 && Number.isFinite(fb) && Math.sign(fb) === Math.sign(f)) (lo = back), (flo = fb);
          break;
        }
      }
    }
    if (!(Number.isFinite(fhi) && Math.sign(fhi) !== Math.sign(flo))) {
      at(x);
      continue; // (no bracket: left as carried)
    }
    // Illinois regula falsi within [lo, hi]
    let side = 0;
    for (let it = 0; it < 100; it++) {
      const m = (lo * fhi - hi * flo) / (fhi - flo), fm = at(m);
      if (!Number.isFinite(fm) || Math.abs(fm) <= tol) break;
      if (Math.sign(fm) === Math.sign(fhi)) {
        (hi = m), (fhi = fm);
        if (side === -1) flo /= 2;
        side = -1;
      } else {
        (lo = m), (flo = fm);
        if (side === 1) fhi /= 2;
        side = 1;
      }
    }
  }
  s.computeConcentrations();
  return worst;
}

// Each surface column's amount from here, want: what the old surface held over the nodes carried,
// plus the new nodes' own start (a widened window's new part, at its start, adds what it holds);
// written into the carried nodes' coverages below, or refused where they can't hold it.
function surfaceTargets(old, s, info) {
  return s.surfCols.map(([k, q], j) => {
    const sf = s.model.surfaces[k];
    const jo = old.model.surfaces.findIndex((o) => info.surfKey(old, o) === info.surfKey(s, sf));
    if (jo < 0) return { want: NaN };
    const qo = old.model.surfaces[jo].species.findIndex((x) => x.name === sf.species[q].name);
    const co = old.surfCols.findIndex(([kk, qq]) => kk === jo && qq === qo);
    if (co < 0) return { want: NaN };
    // (held to what it held only where its window and electrode area are as they were, over the
    // same geometry: a window resized keeps its coverages where they were, its amount following)
    const of = old.model.surfaces[jo], same = info.sameGeometry && sameWindow(old, of, s, sf);
    if (!same) return { want: NaN };
    let fresh = 0;
    const vol = s.model.grid.vol;
    sf.nodes.forEach((g, w) => {
      if (!info.surfNodes.has(g)) fresh += vol[g] * s.surfArea[k][w] * s.th[g * s.nSurf + q];
    });
    return { want: old.surfaceAmount(co) + sf.species[0].capacity * fresh, k, q, nodes: Array.from(sf.nodes).filter((g) => info.surfNodes.has(g)) };
  });
}

// Whether a surface's window and electrode are as they were: a face's, its face (a point); a port's,
// its declared window and area per volume.
function sameWindow(old, of, s, sf) {
  if (of.port < 0 || sf.port < 0) return of.port < 0 && sf.port < 0;
  const a = old.model.ports[of.port], b = s.model.ports[sf.port];
  return a.from === b.from && a.to === b.to && deepEqual(a.area, b.area) && old.model.regions[a.region].length === s.model.regions[b.region].length;
}

// Each surface's carried coverages put at the amounts surfaceTargets asks, every species of a surface
// at once, so the answer can't hang on the order they're listed in (one at a time, a growing species
// listed first couldn't fit in the bare sites the shrinking one still held, while the other listing
// order passed; Neon's third review, 1.4). A column is held where it has an amount to keep and its
// carried nodes hold some of it; the others keep their coverages (multiplied) or levels (scaled).
//   - Multiplied, each coverage by its factor, where that leaves every carried node 1e-3 of its sites
//     bare: exact in one pass.
//   - Otherwise (near full) each held species' odds over the bare sites scaled, θ_q' ∝ λ_q θ_q with
//     the bare sites and the species not held at weight 1: a shift of ln λ_q in q's level, the same
//     at every carried node. Of the profiles that keep the amounts, it's the one nearest the old in
//     relative entropy (Σ_g s_g Σ θ' ln(θ'/θ), the bare sites counted as a species), and its
//     ℓ = ln λ minimise the convex F(ℓ) = Σ_g s_g ln(r_g + Σ_q e^{ℓ_q} θ_qg) − Σ_q ℓ_q W_q (s_g a
//     node's sites, r_g its bare and unheld share, W_q what the carried nodes must hold), whose
//     gradient is the amounts' misfit and whose Hessian, Σ_g s_g (diag θ' − θ'θ'ᵀ), is positive
//     definite. So the minimum exists, and is unique, exactly when each W_q > 0 and Σ W_q is below
//     the carried sites. Newton finds it, each step a descent, cut to where F stops falling along
//     it. The bare fraction comes from the odds, never as 1 − Σθ (which cancelled to 0 at bare 1e-17,
//     and refused a pure regrid).
//   - Where every node is that full, the bare sites' total is the sites less the amounts, to the
//     amounts' round-off: at bare 1e-17 everywhere, two species' wants came to 1.4e-16 more than the
//     sites, and below that, the bare share a solve lands on is round-off's. So F is given a penalty,
//     ½·1e-14 Σ_q W_q L_q² on each level's whole shift L_q over the rounds: an amount kept to 1e-14
//     of itself per unit of L, and a want that's over the sites by its round-off fits (that one's
//     bare sites moved 4% at most; over by 1e-14, an e-fold). Over by more: refused.
function rescaleSurfaces(s, targets) {
  const { n, M } = s, vol = s.model.grid.vol, nS = s.nSurf, RT = s.model.RT;
  const sites = (k, g) => s.model.surfaces[k].species[0].capacity * vol[g] * s.surfArea[k][s.surfW[g]];
  // a held column's amount, A on its carried nodes, and W, what they must hold (its want, less its
  // window's other nodes, held at their own start)
  const columns = () =>
    targets.map((t, j) => {
      if (!(t.want > 0) || !t.nodes || t.nodes.length === 0) return null;
      let fresh = 0;
      for (const g of s.model.surfaces[t.k].nodes) if (!t.nodes.includes(g)) fresh += sites(t.k, g) * s.th[g * nS + t.q];
      const A = s.surfaceAmount(j) - fresh;
      return A > 0 ? { ...t, A, W: t.want - fresh } : null;
    });
  const misfit = (cols) => cols.reduce((m, c) => (c ? Math.max(m, Math.abs(c.W / c.A - 1)) : m), 0);

  const multiplied = (k, held, nodes, bareAfter) => {
    const sf = s.model.surfaces[k], th = new Float64Array(sf.species.length);
    for (const g of nodes) {
      const b = s.blockOfNode[g], bare = bareAfter(g);
      for (let qq = 0; qq < th.length; qq++) th[qq] = s.th[g * nS + qq];
      for (const c of held) th[c.q] *= c.W / c.A;
      // (each species' level is its odds over the bare sites)
      for (let qq = 0; qq < th.length; qq++) if (th[qq] > 0) (s.u[b * M + 1 + n + qq] = sf.species[qq].mu0 / RT + Math.log(th[qq] / bare)), (s.uLo[b * M + 1 + n + qq] = 0);
    }
  };

  const KAPPA = 1e-14, shifts = new Map(); // (each surface's levels' whole shifts L, over the rounds)
  const scaled = (k, held, nodes) => {
    const sf = s.model.surfaces[k], ns = sf.species.length, T = held.length, m1 = T + 1;
    if (!shifts.has(k)) shifts.set(k, new Float64Array(ns));
    const L = shifts.get(k);
    const bad = held.find((c) => !(c.W > 0));
    if (bad) throw new DeviceError(`set(): surface '${sf.name}' can't hold what it held on the new grid (its ${sf.species[bad.q].name} would need ${bad.W.toPrecision(4)} where it was carried): keep its window as it was, or build a new Device to start over`);
    // (each species' place among the held, −1 if not; each node's sites and log odds ζ = η − μ°/RT)
    const pos = new Int32Array(ns).fill(-1);
    held.forEach((c, t) => (pos[c.q] = t));
    const W = Float64Array.from(held, (c) => c.W), sg = Float64Array.from(nodes, (g) => sites(k, g));
    const zeta = new Float64Array(nodes.length * ns);
    nodes.forEach((g, w) => {
      const b = s.blockOfNode[g];
      for (let qq = 0; qq < ns; qq++) zeta[w * ns + qq] = s.u[b * M + 1 + n + qq] + s.uLo[b * M + 1 + n + qq] - sf.species[qq].mu0 / RT;
    });
    const S = sg.reduce((a, x) => a + x, 0), want = W.reduce((a, x) => a + x, 0);
    // (more than the sites hold, by more than its round-off: refused, the device left as it was)
    if (!(want < S * (1 + KAPPA))) throw new DeviceError(`set(): surface '${sf.name}' can't hold what it held (its species would need ${((want - S) / S).toPrecision(2)} more than all its sites): its sites shrank below what's adsorbed; keep them larger, or build a new Device to start over`);
    // The carried nodes' amounts at ℓ, and with wts, F's Hessian as a graph's weights (vertex 0 the
    // ground: wts[(t+1)·m1] = Σ s r θ'_t, to the bare and unheld share; wts[(t+1)·m1 + p+1] = Σ s θ'_t θ'_p).
    const p = new Float64Array(T);
    const amounts = (ell, wts) => {
      const A = new Float64Array(T);
      if (wts) wts.fill(0);
      for (let w = 0; w < nodes.length; w++) {
        let mx = 0;
        for (let qq = 0; qq < ns; qq++) mx = Math.max(mx, zeta[w * ns + qq] + (pos[qq] >= 0 ? ell[pos[qq]] : 0));
        let rest = Math.exp(-mx), den = rest;
        for (let qq = 0; qq < ns; qq++) {
          const e = Math.exp(zeta[w * ns + qq] + (pos[qq] >= 0 ? ell[pos[qq]] : 0) - mx);
          if (pos[qq] >= 0) p[pos[qq]] = e;
          else rest += e;
          den += e;
        }
        for (let t = 0; t < T; t++) {
          const th = p[t] / den;
          A[t] += sg[w] * th;
          if (!wts) continue;
          wts[(t + 1) * m1] += sg[w] * th * (rest / den);
          for (let q2 = 0; q2 < T; q2++) if (q2 !== t) wts[(t + 1) * m1 + q2 + 1] += sg[w] * th * (p[q2] / den);
        }
      }
      return A;
    };
    // (the penalised F's gradient at ℓ, given its amounts there, and its slope along d)
    const grad = (A, ell) => A.map((x, t) => x - W[t] + KAPPA * W[t] * (L[held[t].q] + ell[t]));
    const slope = (G, d) => G.reduce((a, x, t) => a + x * d[t], 0);
    const ell = new Float64Array(T), wts = new Float64Array(m1 * m1), trial = new Float64Array(T);
    const at = (tau, d) => {
      for (let t = 0; t < T; t++) trial[t] = ell[t] + tau * d[t];
      return slope(grad(amounts(trial), trial), d);
    };
    let best = Infinity, stalled = 0;
    for (let it = 0; it < 60 && stalled < 3; it++) {
      const G = grad(amounts(ell, wts), ell);
      const res = G.reduce((m, x, t) => Math.max(m, Math.abs(x) / W[t]), 0);
      if (res < 1e-15) break;
      if (res < 0.5 * best) (best = res), (stalled = 0);
      else stalled++;
      for (let t = 0; t < T; t++) wts[t + 1] = wts[(t + 1) * m1] += KAPPA * W[t];
      const d = solveGrounded(wts, G.map((x) => -x), T);
      // (along d, F's slope rises from below zero: the whole step where F still falls at its end,
      // else a point short of where the slope turns, found by regula falsi, Illinois)
      const g0 = slope(G, d);
      if (!(g0 < 0)) break;
      let lo = 1;
      let ghi = at(1, d);
      if (ghi > 0) {
        let glo = g0, hi = 1, side = 0;
        lo = 0;
        for (let ls = 0; ls < 60 && hi - lo > 1e-12; ls++) {
          const tau = (lo * ghi - hi * glo) / (ghi - glo), g = at(tau, d);
          if (g <= 0) {
            (lo = tau), (glo = g);
            if (side === -1) ghi *= 0.5;
            side = -1;
            if (-g <= -0.5 * g0) break;
          } else {
            (hi = tau), (ghi = g);
            if (side === 1) glo *= 0.5;
            side = 1;
          }
        }
      }
      if (!(lo > 0)) break;
      for (let t = 0; t < T; t++) ell[t] += lo * d[t];
    }
    // each held species' level shifted by its ℓ at every carried node (its low word kept: TwoSum)
    held.forEach((c, t) => (L[c.q] += ell[t]));
    for (const g of nodes) {
      const b = s.blockOfNode[g];
      held.forEach((c, t) => {
        const o = b * M + 1 + n + c.q, x = s.u[o], l = ell[t];
        if (l === 0 || !Number.isFinite(x)) return;
        const sum = x + l;
        s.uLo[o] += Math.abs(x) >= Math.abs(l) ? x - sum + l : l - sum + x;
        s.u[o] = sum;
      });
    }
  };

  let worst = Infinity;
  for (let round = 0; round < 8; round++) {
    s.computeConcentrations();
    const cols = columns();
    worst = misfit(cols);
    if (worst < 1e-14) break;
    // (a surface at a time, every column of it from the same state)
    for (const k of new Set(cols.filter((c) => c && c.W !== c.A).map((c) => c.k))) {
      const held = cols.filter((c) => c && c.k === k), nodes = held[0].nodes;
      const bareAfter = (g) => held.reduce((b, c) => b - s.th[g * nS + c.q] * (c.W / c.A - 1), s.th0[g]);
      if (held.every((c) => c.W > 0) && nodes.every((g) => bareAfter(g) >= 1e-3)) multiplied(k, held, nodes, bareAfter);
      else scaled(k, held, nodes);
    }
  }
  s.computeConcentrations();
  // (what the surfaces held, to the last few digits, or refused: the device left as it was)
  if (!(worst < 1e-14)) worst = misfit(columns());
  if (!(worst < 1e-9)) throw new DeviceError(`set(): a surface can't hold what it held on the new grid (its amount ${worst.toPrecision(2)} off after 8 rounds, its sites near full): keep them larger, or build a new Device to start over`);
}

// L x = b for L the Laplacian of the weights w (an m × m array, m = T + 1, symmetric, vertex 0 the
// ground, x₀ = 0), eliminated from the last vertex down with each pivot summed from its weights,
// never differenced (Grassmann–Taksar–Heyman): a near-full surface's smallest pivot, its bare share,
// is a sum of small positive terms, where one found by subtraction cancels to round-off. Returns
// x₁ … x_T; overwrites w and b.
function solveGrounded(w, b, T) {
  const m = T + 1, d = new Float64Array(m), x = new Float64Array(m), r = new Float64Array(m);
  for (let k = 1; k <= T; k++) r[k] = b[k - 1];
  for (let k = T; k >= 1; k--) {
    let dk = 0;
    for (let j = 0; j < k; j++) dk += w[k * m + j];
    d[k] = dk;
    for (let i = 0; i < k; i++) {
      const wik = w[i * m + k];
      if (wik === 0) continue;
      if (i > 0) r[i] += (wik * r[k]) / dk;
      for (let j = 0; j < k; j++) if (j !== i) w[i * m + j] += (wik * w[k * m + j]) / dk;
    }
  }
  for (let k = 1; k <= T; k++) {
    let v = r[k];
    for (let j = 1; j < k; j++) v += w[k * m + j] * x[j];
    x[k] = v / d[k];
  }
  return x.subarray(1);
}
