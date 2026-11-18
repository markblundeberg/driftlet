// Carrying a solved state across Device.set() (docs/device.md, "set(patch)"): from the old solver
// onto a new one built at the new definition, by what each value means, never by where it sat.
//
// The contract:
//   1. The new solver starts from its own start (its c0, its gates and coverages at rest), and
//      every value that can be matched by meaning overwrites it: a region by its name (by place
//      where it has none, or was renamed: the same count, the same index), a node by where it
//      sits within its region (as a fraction of the region's length), a species by its name, a
//      face by the regions on its two sides, a face's reaction by its stoichiometry, a gate or a
//      surface species by its name, a surface by its owner (a port's name, or its face), a
//      terminal by its name. What can't be matched keeps its own start; a region that had a
//      state and gets none says so (restarted).
//   2. What's carried is physical: concentrations, φ, coverages (as moles: a capacity changed
//      leaves the adsorbed amount, the new sites empty), gates' fractions, a face's fluxes,
//      rates and displacement, a diffuse layer's drop in volts, a metal's flux and its carrier's
//      level. Each is written back through the new definition's own relations (its statistics,
//      inverted; its μ°, c_ref and T), so a change of any of them moves no concentration.
//   3. Between old nodes, a species follows the profile Scharfetter–Gummel takes within the old
//      cell (c e^θ̃ linear in ∫e^θ̃, θ̃ the drift potential the flux uses: zφ̂, the grading, the
//      excess potential, less a flow's Péclet shift), in the cell's own coordinate (∫dx/A), and
//      in logarithms, so a steep cell (8 V across one) neither rounds nor underflows. Coverages
//      and gates are linear, within their own window only.
//   4. A node that nothing changed (the same place, material data, statistics and T) is copied
//      bit for bit, low parts and all.
//   5. Amounts: a stretch that was closed on its own and still is, over the same space, keeps
//      its reference exactly; any other stretch with a carried amount (in a conserved combination,
//      newly closed, split, regridded, resized) takes what the old state held over its regions, on
//      the old grid; and where the state was interpolated or converted, the carried state is put
//      at exactly that amount (a uniform shift of the species' level along the stretch, through
//      its statistics), as are the surfaces'. A strictly neutral node interpolated is put back to
//      neutrality first (its φ̂).
//
// carryState(old, s, { oldDef, newDef, pure }) writes s's state (and its clock, terminals and
// references); it changes nothing else, so a set() that throws from here leaves the device as it
// was. `pure`: the definitions differ only in the grid (a steady state then stands for the old
// definition's, solved again on the new grid before the next step: see Device.advance).

import { crossSection } from './grid.js';
import { DeviceError } from './device.js';
import { gateRate } from './boundaries.js';
import { deepEqual } from './merge.js';

// [ln(1 − s), ln s] for s = expm1(d w)/expm1(d), the weight of the right end at fraction w of a cell
// whose drift potential rises by d: exact in logarithms at any d.
export function sgLogWeights(d, w) {
  if (w <= 0) return [0, -Infinity];
  if (w >= 1) return [-Infinity, 0];
  if (Math.abs(d) < 1e-12) return [Math.log1p(-w), Math.log(w)];
  if (d > 0) {
    const [a, b] = sgLogWeights(-d, 1 - w);
    return [b, a];
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
    const k = rx.part.map((p) => `${species[p.i]?.name ?? `@${p.i}`}:${p.side}:${p.nu}`).sort().join();
    const j = seen.get(k) ?? 0;
    seen.set(k, j + 1);
    return `${k}#${j}`;
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
 * Carry old's state onto s (see the contract above). Returns { unmatched }: the new regions that
 * had no old region to carry from (started from their own c0).
 */
export function carryState(old, s, { oldDef, newDef, pure = false }) {
  const { n, M } = s, oldM = old.M, species = s.model.species;
  const VTo = old.model.RT, VTn = s.model.RT; // (RT: φ̂ = Fφ/RT, so φ = φ̂·RT/F; the F cancels)
  const sameT = VTo === VTn;
  old.computeConcentrations();
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
  const sameStats = new Map();
  const statsSame = (mo, mn) => {
    const k = `${mo}:${mn}`;
    if (!sameStats.has(k)) sameStats.set(k, old.model.materials[mo].name === s.model.materials[mn].name && deepEqual(statsOf(oldDef, old.model, mo), statsOf(newDef, s.model, mn)) && deepEqual(old.model.materials[mo].conductor, s.model.materials[mn].conductor));
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
      if (L === Lo) for (const k of [a, b]) if (og.x[k] - xo0 === ng.x[g] - x0) at = k;
      if (at >= 0) (a = b = at), (w = wx = 0);
      const bg = s.blockOfNode[g], ba = old.blockOfNode[a], bb = old.blockOfNode[b];
      let exact = at >= 0 && sameT && statsSame(old.nodeMaterial[a], s.nodeMaterial[g]) && old.nodeConductor[a] === s.nodeConductor[g] && old.nodeIdeal[a] === s.nodeIdeal[g] && old.background[a] === s.background[g];
      if (exact) for (let i = 0; i < n && exact; i++) {
        const ko = a * n + i, kn = g * n + i;
        if (old.present[ko] !== s.present[kn] || (s.present[kn] && (!Object.is(old.cRef[ko], s.cRef[kn]) || !Object.is(old.mu0hat[ko], s.mu0hat[kn])))) exact = false;
      }
      if (exact) {
        s.u[bg * M] = old.u[ba * oldM];
        s.uLo[bg * M] = old.uLo[ba * oldM];
        for (let i = 0; i < n; i++) {
          s.u[bg * M + 1 + i] = old.u[ba * oldM + 1 + i];
          s.uLo[bg * M + 1 + i] = old.uLo[ba * oldM + 1 + i];
        }
        copied[g] = 1;
        continue;
      }
      // a metal: its flux J (linear), its carrier's level μ̄ (linear in volts)
      const im = s.nodeConductor[g];
      if (im >= 0) {
        if (old.nodeConductor[a] === im && old.nodeConductor[b] === im) {
          const v = (bk, k) => old.u[bk * oldM + k] + old.uLo[bk * oldM + k];
          s.u[bg * M] = (1 - wx) * v(ba, 0) + wx * v(bb, 0);
          s.u[bg * M + 1 + im] = (((1 - wx) * v(ba, 1 + im) + wx * v(bb, 1 + im)) * VTo) / VTn;
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
        if (!old.present[a * n + i] || !old.present[b * n + i]) continue; // (newly present: its own start)
        const ca = old.c[a * n + i], cb = old.c[b * n + i];
        levelAt[g * n + i] = ((1 - w) * (old.u[ba * oldM + 1 + i] + old.uLo[ba * oldM + 1 + i]) + w * (old.u[bb * oldM + 1 + i] + old.uLo[bb * oldM + 1 + i])) * VTo;
        let c = ca;
        if (b > a) {
          moved[g] = 1;
          const k = seg * n + i, D = old.segD[k], pe = D > 0 ? (vel * h) / D : 0;
          const d = old.z[i] * (ph(bb) - ph(ba)) + old.segShift[k] + old.ex[k + n] - old.ex[k] - pe;
          c = sgInterpolate(ca, cb, d, w);
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
      else s.u[bg * M + 1 + i] = fresh[bg * M + 1 + i];
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
  const faceKey = (sv, f) => (sv === old ? `${f}` : map[f] >= 0 && map[f + 1] === map[f] + 1 ? `${map[f]}` : `new${f}`);
  const surfKey = (sv, sf) => (sf.port >= 0 ? `port:${sf.name}` : `face:${faceKey(sv, sf.face)}`);
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
  const surfCarried = new Set();
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
      if (!(sum < 1)) throw new DeviceError(`set(): surface '${sf.name}' now has ${sf.species[0].capacity} mol/m² of sites, fewer than it holds (its coverage there would be ${sum.toPrecision(4)}); keep the capacity above what's adsorbed, or build a new Device to start over`);
      const bare = 1 - sum;
      sf.species.forEach((sp, qq) => {
        s.u[bg * M + 1 + n + qq] = th[qq] > 0 ? sp.mu0 / s.model.RT + Math.log(th[qq] / bare) : fresh[bg * M + 1 + n + qq];
        s.uLo[bg * M + 1 + n + qq] = 0;
      });
      surfCarried.add(j);
    });
  });

  // --- gates: by owner and name
  const gateKey = (sv, gt) => (gt.port === undefined ? `face:${faceKey(sv, gt.f)}:${gt.gate.name}` : `port:${sv.model.ports[gt.port].name}:${gt.gate.name}`);
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
      if (scale === 1) (s.u[bf * M + k] = old.u[bo * oldM + ko]), (s.uLo[bf * M + k] = old.uLo[bo * oldM + ko]);
      else (s.u[bf * M + k] = (old.u[bo * oldM + ko] + old.uLo[bo * oldM + ko]) * scale), (s.uLo[bf * M + k] = 0);
    };
    put(0, 0); // (its displacement, C/m²)
    for (let i = 0; i < n; i++) put(1 + i, 1 + i); // (each species' flux through it)
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
        if (slot >= 0 && od.drawSlot[i] >= 0) put(slot, od.drawSlot[i], r2); // (a level at the plane, μ̄)
      });
    }
  });
  s.computeConcentrations();
  // A diffuse layer carried where its law changed (T, its ε, the edge's composition): its charge
  // held, its drop found (charges don't move in an instant; its potential does). (Held in volts,
  // a T raised by 0.1 K changed the layer's charge against the face's, and a cell at a driven
  // current couldn't take a step of 1e-15 s; Astra's review, chased.)
  for (const d of s.diffList) {
    const qf = map[d.f];
    if (qf < 0 || map[d.f + 1] !== qf + 1) continue;
    const od = old.diffList.find((e) => e.f === qf && e.side === d.side);
    if (!od) continue;
    const want = od.dl.q;
    let x0 = s.u[d.o], q0 = d.dl.q;
    if (Math.abs(q0 - want) <= 1e-14 * Math.abs(want)) continue;
    copied[d.g] = 0; // (its excesses moved: the amounts it shares are put back below)
    let x1 = x0 * (q0 !== 0 ? want / q0 : 1) + (q0 === 0 ? 1e-6 : 0);
    for (let it = 0; it < 60; it++) {
      s.u[d.o] = x1;
      s.computeConcentrations();
      const q1 = d.dl.q;
      if (Math.abs(q1 - want) <= 1e-14 * Math.abs(want) || q1 === q0) break;
      const next = x1 + ((want - q1) * (x1 - x0)) / (q1 - q0);
      (x0 = x1), (q0 = q1), (x1 = next);
    }
  }
  const info = { map, copied, moved, surfCarried, surfExact, sameGeometry, faceKey, surfKey, gatesAtRest };
  carryRest(old, s, info, { pure });
  return info;
}

// The clock, the terminals, neutrality, and what each stretch and surface conserves from here.
function carryRest(old, s, info, { pure }) {
  const { map, copied, sameGeometry } = info, n = s.n, ng = s.model.grid, og = old.model.grid;
  const allCopied = copied.every((v) => v === 1) && info.surfExact;
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
  // area where the geometry changed; a port's, the same geometry only).
  const ratio = (k) => (k === 0 ? ng.area[0] / og.area[0] : ng.area[ng.nNodes - 1] / og.area[og.nNodes - 1]);
  s.terms.forEach((t, k) => {
    const o = old.terms.findIndex((u) => u.name === t.name);
    if (o < 0) return;
    if (s.floating.includes(k)) s.termV[k] = old.termV[o];
    s.termCharge[k] = old.termCharge[o] * (k < 2 && o < 2 ? ratio(k) : 1);
  });
  s._refreshSources();
  // A gate with nothing to carry (new, or on a window's new part) starts at rest at the voltage
  // the carried state puts across it, α/(α + β), as a start does at its own (from the cold
  // start's, a widened window's new membrane had its gates open at 0 mV beside −65; Astra's review).
  for (const gt of info.gatesAtRest) {
    const V = s._gateVoltage(gt), a = gateRate(gt.gate.alpha, V)[0], b = gateRate(gt.gate.beta, V)[0];
    s.u[gt.o] = a + b > 0 ? a / (a + b) : 0;
    s.uLo[gt.o] = 0;
  }
  if (old.contactDEnd || old.portQEnd) {
    s._assembleBookkeeping(Infinity); // (each from the carried state, for any not carried)
    s.contactDEnd = old.contactDEnd ? { left: old.contactDEnd.left * ratio(0), right: old.contactDEnd.right * ratio(1) } : { ...s.contactD };
    s.portQEnd = Float64Array.from(s.portQ);
    if (old.portQEnd && sameGeometry) s.model.ports.forEach((p, k) => {
      const o = old.model.ports.findIndex((q) => q.name === p.name);
      if (o >= 0) s.portQEnd[k] = old.portQEnd[o];
    });
  }
  if (allCopied && !anyFresh(map)) return referencesKept(old, s, info);

  // Neutrality where a strictly neutral node was interpolated or converted (its φ̂, the levels
  // held), then each carried amount exactly; again, as each moves the other a little.
  // (every node of a strictly neutral material, ε = 0, whose concentrations moved: interpolated,
  // or shifted with a carried amount; one carried as it was keeps its own, neutral as solved)
  const moved = info.moved;
  const neutralNodes = () => {
    const out = [];
    for (let g = 0; g < s.nNodes; g++) if (moved[g] && s.nodeIdeal[g] && s.nodeConductor[g] < 0 && !s.phiUndefined[g] && s.model.materials[s.nodeMaterial[g]].epsr === 0) out.push(g);
    return out;
  };
  const targets = stretchTargets(old, s, info);
  for (let round = 0; round < 60; round++) {
    neutralize(s, neutralNodes());
    s.computeConcentrations();
    let off = 0;
    targets.forEach(({ st, want, rescale }) => {
      if (!rescale) return;
      const o = rescaleStretch(s, st, want);
      if (o > 0) for (let g = st.nodes[0]; g <= st.nodes[1]; g++) moved[g] = 1;
      off = Math.max(off, o);
    });
    if (off < 1e-14) break;
  }
  neutralize(s, neutralNodes());
  rescaleSurfaces(old, s, info);
  s.computeConcentrations();
  // What each stretch conserves from here (see the contract, 5).
  s.referenceAmounts = s.stretches.map((st, k) => {
    const t = targets[k];
    if (t.keep !== undefined) return t.keep;
    return !st.connected && Number.isFinite(t.want) ? t.want : s.amount(st);
  });
  s.boundaryIntake.fill(0);
  s.surfaceRef = Float64Array.from(s.surfCols, (_, j) => surfaceTarget(old, s, info, j).keep ?? s.surfaceAmount(j));
  // A region that had a state and got none: said in the next solutions.
  const fresh = s.model.regions.flatMap((rg, r) => (map[r] < 0 ? [rg.name] : []));
  if (fresh.length > 0) s.restarted = `set(): ${fresh.length === s.model.regions.length ? 'no region' : `${fresh.map((nm) => `'${nm}'`).join(', ')}`} matched one before it (by name; or, as many regions as before, by place), so ${fresh.length === s.model.regions.length ? 'the device' : fresh.length > 1 ? 'they' : 'it'} started from ${fresh.length > 1 || fresh.length === s.model.regions.length ? 'their' : 'its'} c0 at t = ${old.time} (the clock carries on). (Said in every solution up to the next steady solve.)`;
  else if (old.restarted) s.restarted = old.restarted;
}

const anyFresh = (map) => map.some((q) => q < 0);

// Every node and surface copied whole: each reference as it was, by stretch and column.
function referencesKept(old, s, info) {
  const targets = stretchTargets(old, s, info);
  s.referenceAmounts = s.stretches.map((st, k) => targets[k].keep ?? (!st.connected && Number.isFinite(targets[k].want) ? targets[k].want : s.amount(st)));
  s.boundaryIntake.fill(0);
  s.surfaceRef = Float64Array.from(s.surfCols, (_, j) => surfaceTarget(old, s, info, j).keep ?? s.surfaceAmount(j));
  if (old.restarted) s.restarted = old.restarted;
}

// What the old state held of species i over old region q (its faces' diffuse layers' share
// included), on the old grid.
function regionAmount(old, q, i) {
  const og = old.model.grid, n = old.n;
  let a = 0;
  for (let g = og.regionStart[q]; g <= og.regionEnd[q]; g++) a += og.vol[g] * old.c[g * n + i];
  old.diffList.forEach((dl, k) => {
    if (dl.g >= og.regionStart[q] && dl.g <= og.regionEnd[q]) a += og.area[dl.g] * old.gam[k * n + i];
  });
  return a;
}

// Each new stretch: the amount the old state held over its regions (want: where every region
// matched one of the same length, in the same geometry, holding it), whether that amount is to be
// written into the state (rescale), and its reference kept exactly (keep: closed on its own before
// and after, over the same space).
function stretchTargets(old, s, info) {
  const { map, copied, sameGeometry } = info;
  const held = (k) => s.terms[k] && s.terms[k].drive.kind === 'V' && !(s.terms[k].drive.R > 0);
  return s.stretches.map((st) => {
    const name = s.model.species[st.species].name, oi = old.model.species.findIndex((sp) => sp.name === name);
    let want = NaN, touched = false;
    const regs = [];
    for (let r = st.regions[0]; r <= st.regions[1]; r++) regs.push(map[r]);
    const same = oi >= 0 && sameGeometry && regs.every((q, j) => q >= 0 && old.model.regions[q].length === s.model.regions[st.regions[0] + j].length && old.model.materials[old.model.regions[q].material].present[oi]);
    if (same) want = regs.reduce((a, q) => a + regionAmount(old, q, oi), 0);
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) if (!copied[g]) touched = true;
    // (open to a held terminal: its level is the terminal's, and a shift would move it off)
    const fedHeld = (st.leftOpen && held(0)) || (st.rightOpen && held(1)) || st.linked.some((k) => held(2 + k));
    const rescale = same && touched && !fedHeld;
    let keep;
    if (same && st.spectator) {
      const k = old.stretches.findIndex((o) => o.species === oi && o.regions[0] === regs[0] && o.regions[1] === regs[regs.length - 1]);
      const contiguous = regs.every((q, j) => q === regs[0] + j);
      if (k >= 0 && contiguous && old.stretches[k].spectator) keep = old.referenceAmounts[k] + old.boundaryIntake[k];
    }
    // (kept, the state is put at that amount where it was interpolated)
    if (keep !== undefined) want = keep;
    return { st, want, rescale, keep };
  });
}

// Shift a stretch's species level uniformly (η + δ at its nodes) until it holds `want`, through
// its statistics (exact in one step where ideal). Returns how far off it was (relative).
function rescaleStretch(s, st, want) {
  if (!(want > 0)) return 0;
  const { n, M } = s, i = st.species;
  let have = s.amount(st);
  const off = Math.abs(have / want - 1);
  if (!(have > 0) || off < 1e-15) return 0;
  const shift = (d) => {
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) if (s.nodeConductor[g] < 0 && s.present[g * n + i]) s.u[s.blockOfNode[g] * M + 1 + i] += d;
    s.computeConcentrations();
    return s.amount(st);
  };
  // secant on ln(amount) against δ, from the ideal step
  let d0 = 0, l0 = Math.log(have), d1 = Math.log(want / have), at = 0;
  for (let it = 0; it < 60; it++) {
    have = shift(d1 - at);
    at = d1;
    const l1 = Math.log(have);
    if (Math.abs(have / want - 1) < 1e-14) break;
    const slope = (l1 - l0) / (d1 - d0);
    const next = d1 + (Math.log(want) - l1) / (slope > 1e-6 ? slope : 1e-6);
    (d0 = d1), (l0 = l1), (d1 = next);
  }
  return off;
}

// A strictly neutral ideal node's φ̂ moved (its levels held) until its charge balances.
function neutralize(s, nodes) {
  const { n, M, c, z } = s;
  for (const g of nodes) {
    const b = s.blockOfNode[g];
    let d = 0;
    for (let it = 0; it < 60; it++) {
      let f = s.rhoFixed[g], df = 0, scale = Math.abs(s.rhoFixed[g]);
      for (let i = 0; i < n; i++) {
        if (!s.present[g * n + i] || z[i] === 0) continue;
        const ci = c[g * n + i] * Math.exp(-z[i] * d);
        f += z[i] * ci;
        df -= z[i] * z[i] * ci;
        scale += Math.abs(z[i]) * ci;
      }
      if (!(df < 0) || Math.abs(f) <= 1e-15 * scale) break;
      d -= f / df;
    }
    if (d === 0 || !Number.isFinite(d)) continue;
    s.u[b * M] += d;
    for (let i = 0; i < n; i++) if (s.present[g * n + i]) c[g * n + i] *= Math.exp(-z[i] * d);
  }
}

// Surface column j's carried amount: kept exactly where its surface was copied whole over the same
// space, else what the old surface held (moles), written into the coverages.
function surfaceTarget(old, s, info, j) {
  const [k, q] = s.surfCols[j], sf = s.model.surfaces[k];
  const jo = old.model.surfaces.findIndex((o) => info.surfKey(old, o) === info.surfKey(s, sf));
  if (jo < 0) return {};
  const qo = old.model.surfaces[jo].species.findIndex((sp) => sp.name === sf.species[q].name);
  const co = old.surfCols.findIndex(([kk, qq]) => kk === jo && qq === qo);
  if (co < 0 || !info.sameGeometry) return {};
  // (what it held: its own reference is stale wherever it's in a conserved combination, as a
  // surface species bound to a dissolved one is, which counts only the sum; Astra's review)
  const nodesSame = sf.nodes.length === old.model.surfaces[jo].nodes.length && sf.nodes.every((g) => info.copied[g]);
  return { keep: old.surfaceAmount(co), carried: !(nodesSame && info.surfExact) };
}

function rescaleSurfaces(old, s, info) {
  const { n, M } = s;
  s.surfCols.forEach((col, j) => {
    const t = surfaceTarget(old, s, info, j);
    if (!t.carried) return;
    const have = s.surfaceAmount(j);
    if (!(have > 0) || !(t.keep > 0)) return;
    const f = t.keep / have;
    if (Math.abs(f - 1) < 1e-15) return;
    const [k, q] = col, sf = s.model.surfaces[k];
    for (const g of sf.nodes) {
      const b = s.blockOfNode[g];
      let sum = 0;
      for (let qq = 0; qq < sf.species.length; qq++) sum += s.th[g * s.nSurf + qq] * (qq === q ? f : 1);
      if (!(sum < 1)) continue;
      const bare = 1 - sum;
      for (let qq = 0; qq < sf.species.length; qq++) {
        const th = s.th[g * s.nSurf + qq] * (qq === q ? f : 1);
        if (th > 0) (s.u[b * M + 1 + n + qq] = sf.species[qq].mu0 / s.model.RT + Math.log(th / bare)), (s.uLo[b * M + 1 + n + qq] = 0);
      }
    }
    s.computeConcentrations();
  });
}
