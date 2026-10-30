
// The numerical core: the Solver's state (laid out in docs/internals.md), its setup, and the
// drivers around the modules that state their contracts: the equations (src/assemble.js, with
// src/regions.js and src/boundaries.js, whose methods are the Solver's own), Newton
// (src/newton.js), time steps (src/timeSteps.js), steady solves (src/steady.js) and the small
// signal (src/impedance.js), the Solver's own too, the cold start (src/start.js) and the laws
// (src/laws.js, src/lawRows.js).

import { Assembly } from './assemble.js';
import { BlockTridiagonal } from './blockTridiagonal.js';
import { bernoulli, SG_FAR, sgFlux } from './bernoulli.js';
import { BoundaryTerms } from './boundaries.js';
import { EPS0, FARADAY } from './constants.js';
import { sourceAt } from './device.js';
import { diffuseScratch, diffuseUnits } from './diffuse.js';
import { DeviceError } from './errors.js';
import { conservedMoieties, findLaws, noCurrent } from './laws.js';
import { newton as newtonFn } from './newton.js';
import { RegionTerms } from './regions.js';
import { SteadySolves } from './steady.js';
import { SmallSignal } from './impedance.js';
import { advance as advanceFn, integrate as integrateFn, step as stepFn } from './timeSteps.js';
import { initFromComposition as initFromCompositionFn } from './start.js';

// A profile's value at x (c0, kf): piecewise linear, constant beyond its ends.
export function profileAt({ x, values: c }, at) {
  if (at <= x[0]) return c[0];
  const last = x.length - 1;
  if (at >= x[last]) return c[last];
  let lo = 0, hi = last;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (x[mid] <= at) lo = mid;
    else hi = mid;
  }
  return c[lo] + ((c[hi] - c[lo]) * (at - x[lo])) / (x[hi] - x[lo]);
}

// A profile's mean over [a, b], exactly (its integral is piecewise quadratic).
function profileMean(p, a, b) {
  if (!(b > a)) return profileAt(p, a);
  const { x, values: v } = p, last = x.length - 1;
  // ∫ from x[0] to t.
  const integral = (t) => {
    if (t <= x[0]) return v[0] * (t - x[0]);
    let sum = 0;
    for (let k = 0; k < last; k++) {
      if (t <= x[k + 1]) return sum + ((v[k] + profileAt(p, t)) / 2) * (t - x[k]);
      sum += ((v[k] + v[k + 1]) / 2) * (x[k + 1] - x[k]);
    }
    return sum + v[last] * (t - x[last]);
  };
  return (integral(b) - integral(a)) / (b - a);
}

export class SolverError extends Error {
  constructor(message, details) {
    super(message);
    this.name = 'SolverError';
    this.details = details;
  }
}

// The species a face's diffuse layer on that side (0 left, 1 right) passes to its reactions:
// each participant there that the layer holds (charged, or on its lattice), but SRH's (refused
// beside a layer), in species order. A flux across the layer lowers each one's level at the plane.
function drawSpecies(itf, side, species) {
  const lat = itf.diffuse.lattice[side], out = new Set();
  for (const rx of itf.reactions) {
    if (rx.srh) continue;
    for (const p of rx.part) if (p.side === side && (species[p.i].z !== 0 || lat?.on[p.i])) out.add(p.i);
  }
  return [...out].sort((a, b) => a - b);
}

export class Solver {
  constructor(model) {
    this.model = model;
    const { grid, species, materials, regions, contacts } = model;
    const n = species.length;
    // Slots per block: φ̂ (or D, or a conductor's segment flux J), one per species, and at a face
    // one per face reaction (its rate), one per gate (its fraction open) and one per diffuse
    // layer (its drop) and one per species that layer passes to a face reaction (its level at the
    // plane: drawSpecies), as many as the busiest face has.
    const nDiff = (itf) => (itf.diffuse ? [0, 1].reduce((m, side) => m + ((side ? itf.diffuse.right : itf.diffuse.left) ? 1 + drawSpecies(itf, side, species).length : 0), 0) : 0);
    const nRx = model.interfaces.reduce((m, itf) => Math.max(m, itf.reactions.length + itf.gates.length + nDiff(itf)), 0);
    // (and at an electrode surface's nodes, a port's or a face's, one per surface species: its coverage's η)
    const nSurf = model.surfaces.reduce((m, sf) => Math.max(m, sf.species.length), 0);
    // (and a membrane port's gates, after them: each one's fraction open)
    const nPortSlots = Math.max(nSurf, model.ports.reduce((m, port) => Math.max(m, port.surface.length + port.gates.length), 0));
    const M = n + 1 + Math.max(nRx, nPortSlots);
    const nNodes = grid.nNodes;
    const nFaces = regions.length - 1;
    const nB = nNodes + nFaces;
    this.n = n;
    this.M = M;
    this.nRx = nRx;
    this.nSurf = nSurf;
    this.nNodes = nNodes;
    this.nFaces = nFaces;
    this.nB = nB;
    this.VT = model.RT / FARADAY; // thermal voltage, V


    // Solver block of each grid node and of each interface flux node.
    this.blockOfNode = new Int32Array(nNodes);
    for (let g = 0; g < nNodes; g++) this.blockOfNode[g] = g + grid.nodeRegion[g];
    this.blockOfFace = new Int32Array(nFaces);
    for (let f = 0; f < nFaces; f++) this.blockOfFace[f] = grid.regionEnd[f] + f + 1;

    // Electrode surfaces (model.surfaces: a port's, through its window, or a face's, at one node):
    // each node's surface (−1: none) and its index among the surface's nodes; the coverages there,
    // θ [g·nSurf + s], now and at the step's start.
    this.surfAt = new Int32Array(nNodes).fill(-1);
    this.surfW = new Int32Array(nNodes);
    model.surfaces.forEach((sf, j) => sf.nodes.forEach((g, w) => ((this.surfAt[g] = j), (this.surfW[g] = w))));
    // Each species of each surface, as a column of what's conserved: [surface, species].
    this.surfCols = model.surfaces.flatMap((sf, j) => sf.species.map((_, q) => [j, q]));
    this.surfColOf = model.surfaces.map((sf, j) => this.surfCols.findIndex(([jj]) => jj === j));
    this.th = new Float64Array(nNodes * nSurf);
    this.thOld = new Float64Array(nNodes * nSurf);
    // The bare fraction θ₀, directly (1 − Σθ cancels as θ → 1); 1 where there's no surface.
    this.th0 = new Float64Array(nNodes).fill(1);
    this.th0Old = new Float64Array(nNodes).fill(1);
    // Gates: each one's slot (at a face, after its reaction rates; at a membrane port's nodes, after
    // the surface's coverages) and value at the step's start. Faces' first, then ports', node by
    // node through each window.
    this.gatePort = new Int32Array(nNodes).fill(-1);
    model.ports.forEach((port, k) => port.gates.length > 0 && port.nodes.forEach((g) => (this.gatePort[g] = k)));
    this.gateList = [
      ...model.interfaces.flatMap((itf, f) => itf.gates.map((gate, q) => ({ f, q, gate, o: this.blockOfFace[f] * M + 1 + n + itf.reactions.length + q }))),
      ...model.ports.flatMap((port, k) => Array.from(port.nodes).flatMap((g) => port.gates.map((gate, q) => ({ port: k, g, q, gate, o: this.blockOfNode[g] * M + 1 + n + port.surface.length + q })))),
    ];
    this.gateOld = new Float64Array(this.gateList.length);
    this.gateIndex = new Int32Array(nFaces); // face f's first gate in gateList
    let j = 0;
    for (let f = 0; f < nFaces; j += model.interfaces[f].gates.length, f++) this.gateIndex[f] = j;
    this.portGateIndex = model.ports.map((port) => ((j += port.nodes.length * port.gates.length), j - port.nodes.length * port.gates.length));
    // Diffuse layers (src/diffuse.js): one per strictly neutral side of a face that has them, in
    // face order, left side first. Each: its face f, side (0 left, 1 right), edge node g and its
    // block b, its drop s = Fψ/RT in slot o of the face's block (after the gates), the layer's
    // constants, its lattice (crowded: { cMax, on }, the species on it; or null, dilute), which
    // species it holds an excess of (the charged, and a lattice's neutral ones, displaced), and
    // its scratch (q, Γ and their derivatives, D, from computeConcentrations). Γ
    // [d·n + i], the excesses, now and at the step's start; diffAt[f] = [left d, right d] (−1: none).
    // And the species it passes to the face's reactions (the drawn: each one's level at the plane
    // in slot drawSlot[i], after both sides' drops, −1 where not drawn; draw[i] = 1 for them), and
    // their diffusivities there, Dm.
    this.diffList = [];
    this.diffAt = model.interfaces.map(() => [-1, -1]);
    model.interfaces.forEach((itf, f) => {
      if (!itf.diffuse) return;
      let slot = 1 + n + itf.reactions.length + itf.gates.length;
      const sides = [0, 1].filter((side) => (side ? itf.diffuse.right : itf.diffuse.left));
      let level = slot + sides.length;
      for (const side of sides) {
        const g = side ? grid.regionStart[f + 1] : grid.regionEnd[f];
        this.diffAt[f][side] = this.diffList.length;
        const lat = itf.diffuse.lattice[side], drawSlot = new Int32Array(n).fill(-1);
        for (const i of drawSpecies(itf, side, species)) drawSlot[i] = level++;
        const draw = drawSlot.some((q) => q >= 0) ? Uint8Array.from(drawSlot, (q) => (q >= 0 ? 1 : 0)) : null;
        const Dm = model.materials[regions[grid.nodeRegion[g]].material].D;
        this.diffList.push({ f, side, g, b: this.blockOfNode[g], slot, o: this.blockOfFace[f] * M + slot, units: diffuseUnits(itf.diffuse.epsr * EPS0, model.RT, FARADAY), dl: diffuseScratch(n), lat, holds: Uint8Array.from(species, (sp, i) => (sp.z !== 0 || lat?.on[i] ? 1 : 0)), drawSlot, draw, Dm });
        slot++;
      }
    });
    this.gam = new Float64Array(this.diffList.length * n);
    this.gamOld = new Float64Array(this.diffList.length * n);

    // Per-node material data, flattened [g·n + i].
    this.present = new Uint8Array(nNodes * n);
    this.cRef = new Float64Array(nNodes * n);
    this.mu0hat = new Float64Array(nNodes * n); // μ°/RT
    this.rhoFixed = new Float64Array(nNodes);
    this.nodeMaterial = new Int32Array(nNodes);
    this.z = Int32Array.from(species, (s) => s.z);
    for (let g = 0; g < nNodes; g++) {
      const reg = regions[grid.nodeRegion[g]];
      const mat = materials[reg.material];
      this.rhoFixed[g] = reg.fixedCharge;
      this.nodeMaterial[g] = reg.material;
      for (let i = 0; i < n; i++) {
        this.present[g * n + i] = mat.present[i];
        this.cRef[g * n + i] = mat.cRef[i];
        this.mu0hat[g * n + i] = mat.mu0[i] / model.RT;
      }
    }

    // Regions where φ is undefined: no charged species, and either ε = 0 or nothing couples the
    // region electrostatically (neutral faces, no gate). φ there gets an identity row.
    // φ is undefined in a conductor or an insertion host, and in a charge-free region (no
    // charged species) unless its electrostatic cluster, the regions joined to it by faces that
    // aren't neutral, reaches something that fixes φ: a region with charged species, a gate or
    // pinned contact, or a face that isn't neutral to a conductor or to a strictly neutral region
    // with ions. (Else φ there is a free constant.)
    this.phiUndefined = new Uint8Array(nNodes);
    const charged = (mat) => species.some((sp, i) => mat.present[i] && sp.z !== 0);
    const pins = (ct) => ct.phi.type === 'capacitive' || ct.phi.type === 'pinned';
    const free = (r) => materials[regions[r].material].phiFree || materials[regions[r].material].epsr === 0;
    const anchored = new Uint8Array(regions.length);
    for (let r0 = 0; r0 < regions.length; r0++) {
      if (free(r0) || anchored[r0]) continue;
      let r1 = r0;
      while (r1 + 1 < regions.length && model.interfaces[r1].phi.type !== 'neutral' && !free(r1 + 1)) r1++;
      let anchor = (r0 === 0 && pins(contacts.left)) || (r1 === regions.length - 1 && pins(contacts.right));
      for (let r = r0; r <= r1; r++) if (charged(materials[regions[r].material])) anchor = true;
      // A neighbour across a face that isn't neutral: a conductor (a gate), or a strictly neutral
      // region whose ions fix its φ.
      const fixes = (r) => materials[regions[r].material].conductor || (!materials[regions[r].material].phiFree && charged(materials[regions[r].material]));
      if (r0 > 0 && model.interfaces[r0 - 1].phi.type !== 'neutral' && fixes(r0 - 1)) anchor = true;
      if (r1 < regions.length - 1 && model.interfaces[r1].phi.type !== 'neutral' && fixes(r1 + 1)) anchor = true;
      for (let r = r0; r <= r1; r++) anchored[r] = anchor ? 1 : 2;
      r0 = r1;
    }
    regions.forEach((reg, r) => {
      const mat = materials[reg.material];
      if (mat.phiFree || (!charged(mat) && (mat.epsr === 0 || anchored[r] !== 1))) {
        for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) this.phiUndefined[g] = 1;
      }
    });
    // A face with φ undefined on both sides holds no charge, whatever its law says.
    model.interfaces.forEach((itf, f) => {
      if (itf.phi.type !== 'neutral' && this.phiUndefined[grid.regionEnd[f]] && this.phiUndefined[grid.regionEnd[f] + 1]) itf.phi = { type: 'neutral' };
      // (a diffuse side has ions, so φ is defined there: the face above never has one)
    });

    // Unknowns as compensated double-doubles, u + uLo. Only differences of η need the extra
    // precision: a majority carrier carrying a small current has a quasi-Fermi step between
    // nodes far below the ulp of η itself (e.g. 1e-19 vs 3.5e-15), and would otherwise carry
    // exactly zero current there. The Jacobian and linear solve stay in plain doubles.
    this.u = new Float64Array(nB * M);
    this.uLo = new Float64Array(nB * M);
    this.uPrev = new Float64Array(nB * M);
    this.uPrevLo = new Float64Array(nB * M);
    this.c = new Float64Array(nNodes * n);
    this.cOld = new Float64Array(nNodes * n);
    // Non-ideal statistics: per node, K = ∂c/∂ζ (n×n) and the excess ζ − ln(c/c_ref), which
    // enters the fluxes as an extra potential. Ideal nodes skip both.
    this.nodeIdeal = new Uint8Array(nNodes);
    this.background = new Float64Array(nNodes);
    for (let g = 0; g < nNodes; g++) {
      const reg = regions[grid.nodeRegion[g]];
      this.nodeIdeal[g] = materials[reg.material].ideal ? 1 : 0;
      this.background[g] = reg.background;
    }
    this.anyNonIdeal = !this.nodeIdeal.every((v) => v === 1);
    // Metal nodes: the carrier's index, and for a node at a charged face of its metal, the face
    // whose displacement is the metal's surface charge (σ = ±D_f), held as a sheet of excess
    // carriers in that node's half-box.
    this.nodeConductor = new Int32Array(nNodes).fill(-1);
    this.sheetFace = new Int32Array(nNodes).fill(-1);
    this.sheetSign = new Int8Array(nNodes);
    // A metal node's φ slot (φ is undefined there) carries instead the carrier flux J through the
    // segment to its right: Ohm's law and continuity in mixed form stay well conditioned however
    // large σ is (eliminating a stiff ohmic chain directly would cancel catastrophically).
    this.conductorLast = new Uint8Array(nNodes);
    for (let g = 0; g < nNodes; g++) {
      const mat = materials[regions[grid.nodeRegion[g]].material];
      if (mat.conductor) {
        this.nodeConductor[g] = mat.conductor.i;
        this.nodeIdeal[g] = 1;
        this.conductorLast[g] = g === grid.regionEnd[grid.nodeRegion[g]] ? 1 : 0;
      }
    }
    model.interfaces.forEach((itf, f) => {
      if (!itf.conductor || itf.conductor.side === 'both' || itf.phi.type === 'neutral') return;
      const g = itf.conductor.side === 'left' ? grid.regionEnd[f] : grid.regionStart[f + 1];
      this.sheetFace[g] = f;
      this.sheetSign[g] = itf.conductor.side === 'left' ? 1 : -1; // σ_metal = D_f on the left, −D_f on the right
    });
    this._activeSlots();
    this.K = this.anyNonIdeal ? new Float64Array(nNodes * n * n) : null;
    // Crowded nodes under eddy mixing: K's slope too (_statSlopes), for its mobility's Jacobian.
    this.nodeDK = new Uint8Array(nNodes);
    regions.forEach((reg, r) => {
      if (reg.mixing > 0) for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) if (!this.nodeIdeal[g]) this.nodeDK[g] = 1;
    });
    this.dK = this.nodeDK.some((v) => v) ? new Float64Array(nNodes * n * n * n) : null;
    if (!this.dK) this.nodeDK = null;
    this.ex = new Float64Array(nNodes * n);
    this.zeta = new Float64Array(n);
    this.scratch = new Map(); // per statistics model: ζ, c, K work arrays
    this.dA = new Float64Array(M); // derivative work vectors over one block's slots
    this.dB = new Float64Array(M);
    this.jL = new Float64Array(M);
    this.jR = new Float64Array(M);
    this.sg = new Float64Array(5); // a Scharfetter–Gummel flux and its derivatives (sgFlux)
    this.time = 0;
    this.lastDt = Infinity;
    this.atSteady = false; // whether the state is a converged steady solve
    this.lin = null; // (set while assembling for DifferenceTerms)
    // Contact bookkeeping, filled by assemble(): particle flux toward +x through each contact,
    // and the displacement there (the metal's surface charge for a neutral link).
    this.contactFlux = { left: new Float64Array(n), right: new Float64Array(n) };
    this.portFlux = model.ports.map(() => new Float64Array(n)); // into the device, mol/(m²·s)
    this.contactD = { left: 0, right: 0 };
    // A capacitive port's charge (its side's, Σ vol·a·σ over its window), as the contacts' D: now,
    // at the step's start and at the previous step's end.
    this.portQ = new Float64Array(model.ports.length);
    this.portQStart = new Float64Array(model.ports.length);
    this.portQOld = new Float64Array(model.ports.length);
    this.portQEnd = null;
    // And per node, the charge per volume the capacitances hold against the window (Σ aσ, over
    // every port with one there: two gates on one channel, say; the ions hold −Σ aσ), its slope
    // Σ a·C per volt, and the ports [terminal, index in its window].
    this.qg = new Float64Array(grid.nNodes);
    this.qgSlope = new Float64Array(grid.nNodes);
    this.qgTerms = Array.from({ length: grid.nNodes }, () => []);
    model.ports.forEach((port, k) => port.capacitance && port.nodes.forEach((g, w) => this.qgTerms[g].push([2 + k, w])));
    this.contactDOld = { left: 0, right: 0 };
    this.contactDStart = { left: 0, right: 0 };
    this.lawForm = new Map(); // (each circuit law: read as its change this step, see src/lawRows.js)
    // Terminals (the two contacts, then the ports): each one's voltage (V), held by its source or
    // floating (driven by a current, or behind a resistance), and its current into the device.
    // Sources are read at sourceTime: a step's end, or now for a steady solve.
    this.terms = model.terminals;
    this.termV = new Float64Array(this.terms.length);
    this.termVPrev = new Float64Array(this.terms.length);
    this.termI = new Float64Array(this.terms.length);
    this.sourceTime = 0;
    this.sourceOverride = new Map(); // terminal → voltage, while a steady solve ramps it
    this.steady = false; // inside a steady solve: sources at the present time, not a step's end
    // The floating ones are extra unknowns, solved with the grid's by bordering (solveBordered in
    // src/lawRows.js).
    this.floating = this.terms.flatMap((t, k) => (t.drive.kind === 'I' || t.drive.R > 0 ? [k] : []));
    this.drivenKey = this._drivenKey();
    this.lawRefs = new Map(); // what a steady solve keeps behind a terminal at no current (src/lawRows.js)
    this.chargeRefs = new Map(); // …and the charge a terminal at no current that only charges keeps
    // Per terminal: B = ∂res/∂V (compact column), C = ∂I/∂x (compact row), ∂I/∂V, and for a
    // floating one its circuit residual (I − I_set, or I − (V_src − V)/R).
    this.termB = this.terms.map(() => new Float64Array(this.sys.size + 1));
    this.termC = this.terms.map(() => new Float64Array(this.sys.size + 1));
    this.termDI = new Float64Array(this.terms.length);
    // (and a circuit row's derivatives by other terminals' voltages: a law in one's row that counts
    // another electrode's charge)
    this.termDX = this.terms.map(() => new Float64Array(this.terms.length));
    this.termRes = new Float64Array(this.terms.length);
    this._refreshSources();
    // A floating terminal starts level with the first held one (no current); where none is held,
    // with the source behind the first resistance (it was NaN there, and a metal contact's
    // carrier, its level from that V, asked for a c0 a conductor can't be given), else 0.
    const heldK = this.terms.findIndex((_, k) => !this.floating.includes(k)), behindR = this.terms.find((t) => t.drive.R > 0);
    const level = heldK >= 0 ? this.termV[heldK] : behindR ? sourceAt(behindR.drive.src, this.time) : 0;
    // (One behind a resistance too: at its own source, a silver cell held at 0 V and 16 V behind R
    // began with 16 V across it; left at 0 V, one behind a −30 V source beside a current began
    // with 30 V across it. Level, the circuit rows move them from a device at rest.)
    for (const k of this.floating) this.termV[k] = level;
    // Accepted steps, most recent first: start time, size, start state (for BDF2 and the
    // error estimate of adaptive stepping).
    this.history = [];
    this.dtNext = undefined;

    // Bulk reactions as flat participant lists: reactants with +ν, products with −ν.
    // Each has its rate constant at every node: its material's, or its profile's mean over the
    // node's box, so that what a sharply varying profile (strongly absorbed light) generates
    // in all is exact on any grid.
    const { x: gx, segLength } = model.grid;
    this.rxs = model.reactions.map((rx) => ({
      kf: rx.kf,
      srh: rx.srh, // per material: SRH kinetics in place of mass action, or null
      kfNode: Float64Array.from(gx, (x, g) => {
        const m = this.nodeMaterial[g], p = rx.kfProfile[m];
        return p ? profileMean(p, x - (g > 0 ? segLength[g - 1] : 0) / 2, x + (g < segLength.length ? segLength[g] : 0) / 2) : rx.kf[m];
      }),
      generation: rx.generation,
      fixedA: rx.fixedA,
      sp: Int32Array.from([...rx.reactants, ...rx.products], (p) => p.i),
      nu: Float64Array.from([...rx.reactants.map((p) => p.nu), ...rx.products.map((p) => -p.nu)]),
    }));

    // Each reacting port's electrode area per volume at its window's nodes: a profile's mean over
    // each node's box, as for a rate constant, with no electrode outside the window (the profile's
    // mean over the part inside, times that part's share of the box, port.fill). And its
    // reactions' rates there, mol/(m²·s).
    this.portArea = model.ports.map((port) =>
      port.area === null ? null : Float64Array.from(port.nodes, (g, w) => (port.area.value ?? profileMean(port.area, port.lo[w], port.hi[w])) * port.fill[w]),
    );
    this.portRates = model.ports.map((port) => port.reactions.map(() => new Float64Array(port.nodes.length)));
    // Each surface's electrode area per volume at its nodes: a port's, or a face's area over its node's box.
    this.surfArea = model.surfaces.map((sf) => (sf.port >= 0 ? this.portArea[sf.port] : Float64Array.from(sf.nodes, (g) => model.grid.area[g] / model.grid.vol[g])));

    // Scales the generation reactions' rates during a steady solve's continuation (else 1).
    this.generationScale = 1;
    this.hasGeneration = this.rxs.some((rx) => rx.generation);

    // The bulk reactions running in each material.
    this.rxsIn = materials.map((_, m) => this.rxs.filter((rx) => rx.kf[m] > 0));

    // The boxes whose residuals the bookkeeping reads (contact and port fluxes: the end nodes and
    // every port's window), with the segments and faces that touch them, per region. After each
    // step only these are evaluated (_assembleBookkeeping), not the whole device.
    const need = new Uint8Array(nNodes);
    need[0] = need[nNodes - 1] = 1;
    for (const port of model.ports) for (const g of port.nodes) need[g] = 1;
    // Nodes of strictly neutral (ε = 0) regions that stay neutral, where φ is defined: interior
    // ones, and edges at neutral faces (a capacitive or pinned face's edge node holds the face's
    // charge), away from the contacts and from ports that hold a level (which replaces a balance
    // row). On a transient step these are solved in better-conditioned terms (see _chargeRows).
    // A port that only exchanges (an O₂ supply, a leak) is solved there too: left out, its window
    // fell back on the ill-conditioned terms, and a closed neutral electrolyte whose potential is
    // set only by a face reaction couldn't take a step.
    // An end node is solved there too if its contact is closed to every ion and leaves φ alone.
    const held = new Uint8Array(nNodes);
    const closed = (side) => model.contacts[side].phi.type === 'neutral' && model.contacts[side].species.every((l, i) => species[i].z === 0 ? l.type !== 'equilibrium' : l.type === 'blocked');
    this.closedEnd = { left: closed('left'), right: closed('right') };
    held[0] = this.closedEnd.left ? 0 : 1;
    held[nNodes - 1] = this.closedEnd.right ? 0 : 1;
    for (const port of model.ports) if (port.species.some((l) => l.type === 'equilibrium')) port.nodes.forEach((g, w) => port.inside[w] && (held[g] = 1));
    this.chargeNode = new Uint8Array(nNodes);
    const neutralFace = (f) => (f === -1 ? this.closedEnd.left : f === model.interfaces.length ? this.closedEnd.right : model.interfaces[f].phi.type === 'neutral');
    regions.forEach((reg, r) => {
      const mat = materials[reg.material];
      if (mat.conductor || mat.epsr !== 0) return;
      const g0 = grid.regionStart[r], g1 = grid.regionEnd[r];
      for (let g = g0; g <= g1; g++) {
        const edgeOk = (g > g0 || neutralFace(r - 1)) && (g < g1 || neutralFace(r));
        if (edgeOk && !held[g] && !this.phiUndefined[g]) this.chargeNode[g] = 1;
      }
    });
    this.chargeNodes = Int32Array.from([...this.chargeNode.keys()].filter((g) => this.chargeNode[g]));
    // An end node of a strictly neutral region that its contact holds (leaving it neutral: φ law
    // bulk or neutral): the contact's current is the ions' charge passing, and their storage
    // carries none under neutrality, nor does the box's net charge change (a bulk law's "D").
    // Read from the box's balance, both come with round-off ∝ c/dt that cancels only roughly,
    // which on short steps swamps the flux (a floating bath on a dilute side rattled). So there,
    // on a transient step, the ions' storage goes in only after the contact has read its current,
    // and the box's charge isn't counted as displacement (see _contact).
    this.lateStorage = new Uint8Array(nNodes);
    [0, nNodes - 1].forEach((g, k) => {
      const mat = materials[regions[grid.nodeRegion[g]].material], law = contacts[k === 0 ? 'left' : 'right'].phi.type;
      if (!mat.conductor && mat.epsr === 0 && (law === 'bulk' || law === 'neutral') && !this.chargeNode[g] && !this.phiUndefined[g] && this.nodeIdeal[g]) this.lateStorage[g] = 1;
    });
    // An end node's species held by a port (its contact blocks it, so the level is the port's):
    // its row is the port's, storage and all, and the contact passes none of it.
    this.portHeldEnd = new Uint8Array(2 * n);
    model.ports.forEach((port) => {
      [0, nNodes - 1].forEach((g, k) => {
        const w = port.nodes.indexOf(g);
        if (w < 0 || !port.inside[w]) return;
        const ct = contacts[k === 0 ? 'left' : 'right'];
        for (let i = 0; i < n; i++) if (port.species[i].type === 'equilibrium' && ct.species[i].type === 'blocked') this.portHeldEnd[k * n + i] = 1;
      });
    });
    this.combining = false;
    this.bookkeeping = regions.map((_, r) => {
      const g0 = grid.regionStart[r], g1 = grid.regionEnd[r], nodes = [], segs = [];
      for (let g = g0; g <= g1; g++) if (need[g]) nodes.push(g);
      for (let q = g0; q < g1; q++) if (need[q] || need[q + 1]) segs.push(q);
      return { nodes, segs };
    });
    this.bookkeepingFaces = [];
    for (let f = 0; f < nFaces; f++) if (need[grid.regionEnd[f]] || need[grid.regionStart[f + 1]]) this.bookkeepingFaces.push(f);

    this._findStretches();
    this.initFromComposition();
    this.referenceAmounts = this.stretches.map((st) => this.amount(st));
    // And what each surface species held then (a conserved combination through one counts it).
    this.surfaceRef = Float64Array.from(this.surfCols, (_, j) => this.surfaceAmount(j));
    // ∫ (flux in − flux out) dt through the contacts, per stretch, since the reference.
    this.boundaryIntake = new Float64Array(this.stretches.length);
  }

  // The slots that are unknowns, block by block; the rest only ever have identity rows and a
  // zero residual, so they're left out of the linear system (and never change).
  _activeSlots() {
    const { n, M, nB, model } = this;
    const active = (this.active = new Uint8Array(nB * M));
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g], im = this.nodeConductor[g];
      if (im >= 0) {
        active[b * M] = this.conductorLast[g] ? 0 : 1; // the segment flux J
        active[b * M + 1 + im] = 1;
        continue;
      }
      active[b * M] = this.phiUndefined[g] ? 0 : 1;
      for (let i = 0; i < n; i++) active[b * M + 1 + i] = this.present[g * n + i];
      const k = this.surfAt[g];
      if (k >= 0) for (let q = 0; q < model.surfaces[k].species.length; q++) active[b * M + 1 + n + q] = 1;
      const kg = this.gatePort[g], pg = model.ports[kg];
      if (kg >= 0) for (let q = 0; q < pg.gates.length; q++) active[b * M + 1 + n + pg.surface.length + q] = 1;
    }
    model.interfaces.forEach((itf, f) => {
      const bf = this.blockOfFace[f];
      active[bf * M] = itf.phi.type === 'neutral' ? 0 : 1;
      for (let i = 0; i < n; i++) active[bf * M + 1 + i] = itf.links[i].type === 'blocked' ? 0 : 1;
      for (let k = 0; k < itf.reactions.length + itf.gates.length; k++) active[bf * M + 1 + n + k] = 1;
    });
    for (const d of this.diffList) {
      active[d.o] = 1;
      for (const q of d.drawSlot) if (q >= 0) active[this.blockOfFace[d.f] * M + q] = 1;
    }
    const sizes = new Int32Array(nB), loc = (this.loc = new Int32Array(nB * M).fill(-1));
    for (let b = 0; b < nB; b++) for (let r = 0; r < M; r++) if (active[b * M + r]) loc[b * M + r] = sizes[b]++;
    const sys = (this.sys = new BlockTridiagonal(nB, sizes));
    // A pivot that cancels to exactly nothing (a population held only through conductances
    // 1e-16 of its neighbours', as GaAs's minority carriers ~1 per m³ beside a face) is
    // perturbed, not thrown on: Newton then converges, its small updates meaning small residuals
    // (the factorisation only amplifies), refining its solves where they stall.
    sys.staticPivots = true;
    const N = sys.size;
    this.SINK = N;
    const rix = (this.rix = new Int32Array(nB * M).fill(N)), fullOf = (this.fullOf = new Int32Array(N));
    for (let k = 0; k < nB * M; k++) {
      if (loc[k] < 0) continue;
      rix[k] = sys.offX[Math.floor(k / M)] + loc[k];
      fullOf[rix[k]] = k;
    }
    this.res = new Float64Array(N + 1);
    this.delta = new Float64Array(N + 1);
    this.rowDiv = new Float64Array(N + 1).fill(1); // what each row is divided by in the factorised system
  }

  _factor() {
    this.sys.factor();
  }

  // delta = J⁻¹ rhs (compact vectors).
  _solveLinear(rhs, delta) {
    this.sys.solve(rhs, delta);
  }

  // Connected stretches of regions where a species is present. A stretch not connected to
  // a contact (and, later, with no reaction) is a conserved spectator.
  _findStretches() {
    const { model, n } = this;
    const { regions, materials, contacts, grid } = model;
    const last = regions.length - 1;
    this.stretches = [];
    this.stretchOf = new Int32Array(regions.length * n).fill(-1);
    for (let i = 0; i < n; i++) {
      let r = 0;
      while (r <= last) {
        if (!materials[regions[r].material].present[i]) { r++; continue; }
        const r0 = r;
        while (r + 1 <= last && materials[regions[r + 1].material].present[i] && model.interfaces[r].links[i].type !== 'blocked') r++;
        const touches = (ct) => ct.species[i].type !== 'blocked';
        const leftOpen = r0 === 0 && touches(contacts.left);
        const rightOpen = r === last && touches(contacts.right);
        let reactive = false;
        // made or consumed by a reaction at any face it touches
        for (let f = Math.max(0, r0 - 1); f <= Math.min(model.interfaces.length - 1, r); f++) {
          for (const rx of model.interfaces[f].reactions) {
            if (rx.part.some((p) => p.i === i && f + p.side >= r0 && f + p.side <= r)) reactive = true;
          }
        }
        for (let q = r0; q <= r; q++) {
          this.stretchOf[q * n + i] = this.stretches.length;
          const m = regions[q].material;
          for (const rx of model.reactions) {
            if (rx.kf[m] > 0 && [...rx.reactants, ...rx.products].some((p) => p.i === i)) reactive = true;
          }
        }
        const inside = (port) => port.region >= r0 && port.region <= r;
        const linked = model.ports.flatMap((port, k) => (inside(port) && port.species[i].type !== 'blocked' ? [k] : []));
        // made or consumed by a port's reactions (an electrode spread through its window)
        const reacting = model.ports.flatMap((port, k) => (inside(port) && port.reactions.some((rx) => rx.part.some((p) => p.side === 0 && p.i === i)) ? [k] : []));
        if (reacting.length > 0) reactive = true;
        const connected = leftOpen || rightOpen || linked.length > 0;
        this.stretches.push({
          leftOpen,
          rightOpen,
          linked, // ports that hold or exchange it
          ports: [...new Set([...linked, ...reacting])], // ports that bring any in
          species: i,
          regions: [r0, r],
          nodes: [grid.regionStart[r0], grid.regionEnd[r]],
          connected,
          reactive,
          // Conserved on its own: not fed by a contact, not made or consumed by a reaction.
          spectator: !connected && !reactive,
          contactFed: connected, // reached directly by a contact
        });
        r++;
      }
    }
    // Which stretches are fed: those in no conserved combination. A conserved combination is a
    // weighting w of the stretches' amounts that no reaction changes (w·ν = 0 for every face and
    // bulk reaction) and that nothing outside feeds (w = 0 on stretches reached by a contact or
    // a port). So Ag⁺ between silver electrodes is fed (Ag⁺ + e⁻ ⇌ Ag(s), the electrons fed by the
    // contacts), while Fe³⁺ and Fe²⁺ between platinum electrodes aren't (Fe³⁺ + e⁻ ⇌ Fe²⁺
    // conserves the iron, whatever the electrons do).
    this.rough = this.stretches.map((st) => this._roughAmount(st));
    this.moieties = conservedMoieties(this._structure());
    this.stretches.forEach((st, k) => {
      st.connected = st.contactFed || this.moieties.every((w) => w[k] === 0);
    });
    for (const st of this.stretches) {
      st.spectator = !st.connected && !st.reactive;
      // Mobile throughout: its steady state is then a single level fixed by its amount.
      st.mobile = true;
      for (let q = st.regions[0]; q <= st.regions[1]; q++) if (!(materials[regions[q].material].D[st.species] > 0)) st.mobile = false;
    }
    // Stretches that carry nothing in a steady state, so whose level is flat there: a species that
    // no reaction or port touches, reached by one contact at most, moving through ideal regions
    // without flow or mixing (where carrying nothing is a level flat). Found through their own
    // conduction, which can be all but nothing, such a level can't be moved (a MOS capacitor's
    // inversion electrons, which reach the back contact only through a bulk with ~1e3 of them per
    // cm³; blocked SO₄²⁻ driven from a zinc cathode's extended space charge, down to 1e-25 mol/m³
    // there, whose chain of conductances lost the system 14 digits). So steady solves hold it
    // flat (flatRows in src/lawRows.js): at its contact's level, or reaching none, at the level
    // its amount fixes.
    // (Not through flow or mixing, nor a concentrated material's cross-diffusion, where carrying
    // nothing isn't a level flat.)
    const plain = (st) => {
      if (st.ports.length > 0 || !st.mobile || (st.leftOpen && st.rightOpen)) return false;
      for (let q = st.regions[0]; q <= st.regions[1]; q++) {
        const reg = regions[q], mat = materials[reg.material];
        if (mat.conductor || !mat.ideal || reg.velocity !== 0 || reg.mixing > 0) return false;
      }
      return true;
    };
    this.flatStretches = this.stretches.filter((st) => !st.reactive && plain(st));
    // Likewise a group of stretches that react only among themselves, each reaching the same one
    // contact and nothing else, their reactions balanced at that contact's levels (no light, no
    // other reservoir): in a steady state, one reservoir drives nothing, so the group is at
    // equilibrium, each level flat at the contact's. A MOS capacitor's electrons and the holes
    // they recombine with reach only its back contact: moved, slow generation through a bulk of
    // GaAs with ~1e-4 electrons per cm³ couldn't carry the inversion layer's level with it.
    const group = this.stretches.map((_, k) => k), root = (k) => (group[k] === k ? k : (group[k] = root(group[k])));
    const join = (ks) => ks.forEach((k) => (group[root(k)] = root(ks[0])));
    const rxns = []; // [stretches, A/RT at the contact's levels, as a function of the side]
    regions.forEach((reg, q) => {
      for (const rx of this.rxsIn[reg.material]) {
        const ks = rx.sp.map((i) => this.stretchOf[q * n + i]);
        join(ks);
        rxns.push([ks, (side) => rx.sp.reduce((a, i, p) => a + rx.nu[p] * this.contactEta(side, i), rx.fixedA)]);
      }
    });
    model.interfaces.forEach((itf, f) => {
      for (const rx of itf.reactions) {
        const ks = rx.part.map((p) => this.stretchOf[(f + p.side) * n + p.i]);
        if (ks.some((k) => k < 0)) continue;
        join(ks);
        rxns.push([ks, (side) => rx.part.reduce((a, p) => a - p.nu * this.contactEta(side, p.i), rx.fixedA)]); // (a face's ν: products +, so − Σ ν η)
      }
    });
    const groups = new Map();
    this.stretches.forEach((st, k) => {
      if (!st.reactive) return;
      if (!groups.has(root(k))) groups.set(root(k), []);
      groups.get(root(k)).push(k);
    });
    for (const ks of groups.values()) {
      const sts = ks.map((k) => this.stretches[k]), side = sts[0].leftOpen ? 'left' : 'right';
      if (!sts.every((st) => plain(st) && st.leftOpen !== st.rightOpen && (side === 'left' ? st.leftOpen : st.rightOpen))) continue;
      if (!rxns.every(([rk, A]) => !rk.some((k) => ks.includes(k)) || Math.abs(A(side)) < 1e-9)) continue;
      this.flatStretches.push(...sts);
    }
    this._laws();
  }

  // A stretch's amount from its regions' starting concentrations (c0), per unit area, before
  // any state exists: Infinity where one isn't given (set by a level, or doping), to order it last.
  _roughAmount(st) {
    const { regions, grid } = this.model;
    let t = 0;
    for (let q = st.regions[0]; q <= st.regions[1]; q++) {
      const reg = regions[q], c0 = reg.c0Profile[st.species] ? profileAt(reg.c0Profile[st.species], grid.x[grid.regionStart[q]]) : reg.c0[st.species];
      if (!(c0 >= 0)) return Infinity;
      t += c0 * (grid.x[grid.regionEnd[q]] - grid.x[grid.regionStart[q]]);
    }
    return t;
  }

  // The laws a steady solve keeps (src/laws.js finds them, from the structure alone), set on the
  // solver with the state that goes with them; again after a change of drive.
  _laws() {
    const found = findLaws(this._structure());
    for (const [k, st] of this.stretches.entries()) st.conserved = found.conserved.has(k);
    this.constraints = found.constraints;
    this.affinities = found.affinities;
    this.localConstraints = found.localConstraints;
    this.islands = found.islands;
    this.timeLaws = found.timeLaws;
    this.contactLaw = found.contactLaw;
    this.carriers = found.carriers;
    this.unkept = found.unkept; // (laws no row here can keep: see _directSteady)
    // (room in each law's derivative for the diffuse layers' excesses, which src/laws.js, reading
    // the structure alone, doesn't count: each layer's drop and its edge node's slots)
    if (this.diffList.length > 0) {
      const more = this.diffList.length * (this.M + 1);
      for (const cs of [...this.constraints, ...this.timeLaws]) {
        const idx = new Int32Array(cs.idx.length + more), w = new Float64Array(cs.w.length + more);
        idx.set(cs.idx);
        w.set(cs.w);
        cs.idx = idx;
        cs.w = w;
      }
    }
    this.contactQ = [null, null]; // (each contact's charge as _contact captures it, for contactLaw)
    this.pinnedIslands = [];
    this.islandsOn = false; // (newton() turns them on, see there)
    this.pins = [];
  }

  // What src/laws.js reads (see its contract): the device's structure, never the state.
  _structure() {
    const { model, n, M, z, nNodes, stretches, stretchOf, flatStretches, rough, surfCols, surfColOf, terms, floating } = this;
    const { blockOfNode, blockOfFace, rxsIn, nodeMaterial, phiUndefined } = this;
    return { model, n, M, z, nNodes, stretches, stretchOf, flatStretches, rough, surfCols, surfColOf, terms, floating, blockOfNode, blockOfFace, rxsIn, nodeMaterial, phiUndefined };
  }

  // An electrode surface at its starting coverages at node g: η = μ°/RT + ln(θ/θ₀).
  _surfaceStart(g) {
    const { n, M, u, model } = this, surf = model.surfaces[this.surfAt[g]].species, b = this.blockOfNode[g];
    const bare = 1 - surf.reduce((t, sp) => t + sp.theta0, 0);
    surf.forEach((sp, q) => (u[b * M + 1 + n + q] = sp.mu0 / model.RT + Math.log(sp.theta0 / bare)));
  }

  /** Amount (mol per unit area) of surface column j's species, Γ Σ v·a·θ, in the current state. */
  surfaceAmount(j) {
    const [k, q] = this.surfCols[j], sf = this.model.surfaces[k], vol = this.model.grid.vol;
    let s = 0;
    sf.nodes.forEach((g, w) => (s += vol[g] * this.surfArea[k][w] * this.th[g * this.nSurf + q]));
    return sf.species[0].capacity * s;
  }

  /** Total amount (mol per unit area) of a stretch's species in the current state, its faces' diffuse layers' included. */
  amount(stretch) {
    const { n, c, model } = this;
    const vol = model.grid.vol;
    let s = 0;
    for (let g = stretch.nodes[0]; g <= stretch.nodes[1]; g++) s += vol[g] * c[g * n + stretch.species];
    return s + this._diffuseHeld(stretch, this.gam);
  }

  // What the faces' diffuse layers hold of a stretch's species (mol, per area if planar), by the
  // excesses gam (now, or gamOld at the step's start).
  _diffuseHeld(stretch, gam) {
    let s = 0;
    this.diffList.forEach((dl, k) => {
      if (dl.g >= stretch.nodes[0] && dl.g <= stretch.nodes[1]) s += this.model.grid.area[dl.g] * gam[k * this.n + stretch.species];
    });
    return s;
  }

  // After the terminals' drives change (Device.set): which are floating, and the held values. A
  // change of source is a discontinuity, so time stepping restarts its order.
  redrive() {
    this.atSteady = false;
    const before = this.floating, driven = this.drivenKey;
    this.floating = this.terms.flatMap((t, k) => (t.drive.kind === 'I' || t.drive.R > 0 ? [k] : []));
    if (this.floating.length !== before.length) this.deltaV = null;
    // A terminal newly driven by a current (or no longer) can close a combination behind it.
    // (and one put behind a resistance, or taken from it, its time steps' reading: see src/laws.js)
    this.drivenKey = this._drivenKey();
    if (this.drivenKey !== driven) this._laws();
    this._refreshSources();
    this.history = [];
  }

  // How each terminal is driven, as far as the laws care (see src/laws.js): by a voltage or a current,
  // behind a resistance or not, and a current that's held at zero throughout (which feeds a
  // combination nothing, so a reaction's affinity may stand in for its summed balance).
  _drivenKey() {
    return this.terms.map((t) => `${t.drive.kind}${t.drive.R > 0 ? 'R' : ''}${noCurrent(t) ? '0' : ''}`).join();
  }

  // Held terminals take their sources' values at sourceTime (a floating one keeps its voltage,
  // an unknown; one behind a resistance reads its source voltage in _circuit).
  _refreshSources() {
    this.terms.forEach((t, k) => {
      if (t.drive.kind !== 'V' || t.drive.R > 0) return;
      this.termV[k] = this.sourceOverride.has(k) ? this.sourceOverride.get(k) : sourceAt(t.drive.src, this.sourceTime, this.sourceBefore);
    });
  }

  /** η_i (μ̄/RT) of a port's outside level for species i (its electrode's carrier: at V itself). */
  portEta(port, i) {
    const link = port.species[i];
    const V = this.termV[2 + this.model.ports.indexOf(port)];
    return this.z[i] === 0 ? link.mu / this.model.RT : (this.z[i] * (V + (link.offset ?? 0))) / this.VT;
  }

  /** η_i/RT that a fixed contact link imposes (or NaN if the link isn't fixed). */
  contactEta(side, i) {
    const link = this.model.contacts[side].species[i];
    if (link.type !== 'equilibrium' && link.type !== 'velocity') return NaN;
    const z = this.z[i];
    return z === 0 ? link.mu / this.model.RT : (z * (this.termV[side === 'left' ? 0 : 1] + link.offset)) / this.VT;
  }

  // The cold start (src/start.js).
  initFromComposition() {
    return initFromCompositionFn(this);
  }

  // Newton's iteration on the equations at dt (src/newton.js).
  newton(dt, opts) {
    return newtonFn(this, dt, opts);
  }

  /**
   * Advance by dt, splitting the interval into halves (recursively) wherever Newton fails.
   * The end time is always t + dt; the result reports the substeps taken.
   */
  // A time step, halved until it converges (src/timeSteps.js).
  advance(dt, opts, depth = 0) {
    return advanceFn(this, dt, opts, depth);
  }

  // One time step (src/timeSteps.js).
  step(dt, opts = {}) {
    return stepFn(this, dt, opts);
  }

  // Adaptive time stepping to tEnd (src/timeSteps.js).
  integrate(tEnd, opts = {}) {
    return integrateFn(this, tEnd, opts);
  }

  // Probes for a transient's trace: { x, species, quantity: 'c' | 'V', region } read after every
  // accepted step, linearly between the two nodes of x's region around it. At an interface x
  // belongs to both regions, so `region` (a name or index) picks the side; it defaults to the
  // first region that holds x.
  _probes(probes) {
    if (probes === undefined) return null;
    const { grid, regions, species } = this.model;
    const fail = (m) => {
      throw new DeviceError(m);
    };
    if (!Array.isArray(probes)) fail('advance: probes must be an array of { x, species, quantity, region }, or { interface, species } / { interface, gate }');
    const plan = probes.map((p, k) => {
      const path = `advance: probes[${k}]`;
      if (p?.interface !== undefined) {
        // At a face: a species' flux through it (mol/(m²·s), toward +x), or a gate's fraction open.
        const f = p.interface, itf = this.model.interfaces[f];
        if (!(Number.isInteger(f) && itf)) fail(`${path}.interface: no interface ${JSON.stringify(f)} (they're numbered from 0, left to right)`);
        if (p.gate !== undefined) {
          const q = itf.gates.findIndex((gt) => gt.name === p.gate);
          if (q < 0) fail(`${path}.gate: interfaces[${f}] has no gate ${JSON.stringify(p.gate)}${itf.gates.length ? ` (it has ${itf.gates.map((gt) => gt.name).join(', ')})` : ''}`);
          return { o: this._gateSlot(f, q) };
        }
        const i = species.findIndex((sp) => sp.name === p.species);
        if (i < 0) fail(`${path}.species: no species ${JSON.stringify(p.species)}`);
        return { o: this.blockOfFace[f] * this.M + 1 + i };
      }
      const quantity = p?.quantity ?? 'c';
      if (quantity !== 'c' && quantity !== 'V' && quantity !== 'phi') fail(`${path}.quantity must be 'c' (mol/m³), 'V' (the species voltage) or 'phi' (φ, V)`);
      const i = quantity === 'phi' ? -1 : species.findIndex((sp) => sp.name === p?.species);
      if (quantity !== 'phi' && i < 0) fail(`${path}.species: no species ${JSON.stringify(p?.species)}`);
      if (quantity === 'V' && species[i].z === 0) fail(`${path}: '${p.species}' is neutral, so it has no voltage; read 'c'`);
      const x = grid.x, inside = (r) => p.x >= x[grid.regionStart[r]] && p.x <= x[grid.regionEnd[r]];
      let r;
      if (p.region !== undefined) {
        r = typeof p.region === 'number' ? p.region : regions.findIndex((reg) => reg.name === p.region);
        if (!(r >= 0 && r < regions.length)) fail(`${path}.region: no region ${JSON.stringify(p.region)}`);
        if (!inside(r)) fail(`${path}.x (${p.x} m) is outside region ${JSON.stringify(p.region)}`);
      } else {
        // At a face, the side where the probe has something to read; both, and it must be told.
        const has = (rr) => {
          const mat = this.model.materials[regions[rr].material];
          return i < 0 ? !mat.conductor && !mat.phiFree : mat.present[i];
        };
        const sides = regions.flatMap((_, rr) => (inside(rr) ? [rr] : []));
        if (!(Number.isFinite(p.x) && sides.length)) fail(`${path}.x must be a position in the device, in m (got ${p.x})`);
        const able = sides.filter(has);
        if (able.length > 1) fail(`${path}.x (${p.x} m) is on the face between ${able.map((rr) => JSON.stringify(regions[rr].name)).join(' and ')}: give region, the side to read`);
        r = able[0] ?? sides[0];
      }
      if (i >= 0 && !this.model.materials[regions[r].material].present[i]) fail(`${path}: '${p.species}' is absent from region ${JSON.stringify(regions[r].name)}, so there's nothing to read at x = ${p.x} m`);
      let g = grid.regionStart[r];
      while (g + 1 < grid.regionEnd[r] && x[g + 1] < p.x) g++;
      const w = x[g + 1] > x[g] ? (p.x - x[g]) / (x[g + 1] - x[g]) : 0;
      return { i, g, w, quantity };
    });
    const out = plan.map(() => []);
    const { n, M } = this;
    const at = (q, g) => {
      if (q.quantity === 'phi') return this.phiUndefined[g] ? NaN : this.VT * (this.u[this.blockOfNode[g] * M] + this.uLo[this.blockOfNode[g] * M]);
      if (!this.present[g * n + q.i]) return NaN;
      if (q.quantity === 'c') return this.c[g * n + q.i];
      const o = this.blockOfNode[g] * M + 1 + q.i;
      return (this.VT * (this.u[o] + this.uLo[o])) / this.z[q.i];
    };
    return {
      out,
      read: () => plan.map((q) => (q.o !== undefined ? (this.loc[q.o] < 0 ? 0 : this.u[q.o]) : (1 - q.w) * at(q, q.g) + q.w * at(q, q.g + 1))),
      push: (values) => values.forEach((v, k) => out[k].push(v)),
    };
  }

  // Start from another solver's state for the same device on another grid (a refined one): each
  // node's unknowns interpolated linearly in x within its region, each face's copied, and the
  // terminals' voltages. A warm start for a steady solve; nothing else is carried over.
  _warmFrom(src) {
    const { M, u, uLo } = this, grid = this.model.grid, sg = src.model.grid;
    if (src.M !== M || sg.regionStart.length !== grid.regionStart.length) throw new Error('_warmFrom: a different device');
    const value = (b, k) => src.u[b * M + k] + src.uLo[b * M + k];
    for (let r = 0; r < grid.regionStart.length; r++) {
      let a = sg.regionStart[r];
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const x = grid.x[g];
        while (a + 1 < sg.regionEnd[r] && sg.x[a + 1] < x) a++;
        const b = Math.min(a + 1, sg.regionEnd[r]), h = sg.x[b] - sg.x[a];
        const w = h > 0 ? Math.min(1, Math.max(0, (x - sg.x[a]) / h)) : 0;
        const ba = src.blockOfNode[a], bb = src.blockOfNode[b], bg = this.blockOfNode[g];
        for (let k = 0; k < M; k++) {
          u[bg * M + k] = (1 - w) * value(ba, k) + w * value(bb, k);
          uLo[bg * M + k] = 0;
        }
      }
    }
    for (let f = 0; f < this.nFaces; f++) {
      for (let k = 0; k < M; k++) {
        u[this.blockOfFace[f] * M + k] = value(src.blockOfFace[f], k);
        uLo[this.blockOfFace[f] * M + k] = 0;
      }
    }
    this.termV.set(src.termV);
    this.solvedV = src.solvedV;
    this.coldStart = src.coldStart;
    this.time = src.time;
    this.computeConcentrations();
  }

  _snapshot() {
    return {
      u: Float64Array.from(this.u),
      uLo: Float64Array.from(this.uLo),
      cOld: Float64Array.from(this.cOld),
      thOld: Float64Array.from(this.thOld),
      th0Old: Float64Array.from(this.th0Old),
      gamOld: Float64Array.from(this.gamOld),
      gateOld: Float64Array.from(this.gateOld),
      time: this.time,
      lastDt: this.lastDt,
      atSteady: this.atSteady,
      contactDStart: this.contactDStart,
      contactDOld: this.contactDOld,
      contactDEnd: this.contactDEnd,
      portQStart: this.portQStart,
      portQOld: this.portQOld,
      portQEnd: this.portQEnd,
      boundaryIntake: Float64Array.from(this.boundaryIntake),
      history: this.history.slice(),
      termV: Float64Array.from(this.termV),
      sourceTime: this.sourceTime,
      coldStart: this.coldStart,
    };
  }

  _restore(s) {
    this.u.set(s.u);
    this.uLo.set(s.uLo);
    this.cOld.set(s.cOld);
    this.thOld.set(s.thOld);
    this.th0Old.set(s.th0Old);
    this.gamOld.set(s.gamOld);
    this.gateOld.set(s.gateOld);
    this.time = s.time;
    this.lastDt = s.lastDt;
    this.atSteady = s.atSteady;
    this.contactDStart = s.contactDStart;
    this.contactDOld = s.contactDOld;
    this.contactDEnd = s.contactDEnd;
    this.portQStart = s.portQStart;
    this.portQOld = s.portQOld;
    this.portQEnd = s.portQEnd;
    this.boundaryIntake.set(s.boundaryIntake);
    this.history = s.history.slice();
    this.termV.set(s.termV);
    this.sourceTime = s.sourceTime;
    this.coldStart = s.coldStart;
    this.computeConcentrations();
  }

  // What a floating terminal k passes into the device, which its circuit row holds exactly once
  // converged: the current it's driven by, or (V_source − V)/R behind a resistance (the source
  // as the last step saw it, before any jump at its end). Read as a flux through a metal instead,
  // it carried that flux's round-off, G Δη with G huge: a gate behind R read 1.5e-4 off wherever
  // Newton stopped an update early, and below ~1e-16 A/m², nothing but noise.
  passedCurrent(k) {
    const d = this.terms[k].drive, src = sourceAt(d.src, this.time, !this.atSteady);
    return d.kind === 'I' ? src : (src - this.termV[k]) / d.R;
  }

  // Terminal current toward +x at the right contact, from the last assembly at this state.
  _terminalCurrent() {
    if (this.floating.includes(1)) return -this.passedCurrent(1);
    let I = 0;
    for (let i = 0; i < this.n; i++) I += FARADAY * this.z[i] * this.contactFlux.right[i];
    const dt = this.lastDt;
    return I + (Number.isFinite(dt) ? (this.contactD.right - this.contactDOld.right) / dt : 0);
  }

  // Each species' flux toward +x across every segment, at the current state (concentrations
  // computed): totals through the cross-section, mol/s (mol/(m²·s) for a planar device). A face's
  // segment, between its pair of nodes, carries the face's flux; a conductor's, its carrier's.
  // These are the fluxes the balances add up (displacement aside), from the same formulas.
  segmentFluxes() {
    const { model, n, M, u, uLo, c, z, ex } = this;
    const { grid, materials, regions, interfaces } = model;
    const out = Array.from({ length: n }, () => new Float64Array(grid.nNodes - 1));
    for (let r = 0; r < regions.length; r++) {
      const reg = regions[r], mat = materials[reg.material], vel = reg.velocity;
      for (let s = grid.regionStart[r]; s < grid.regionEnd[r]; s++) {
        const bL = s + r, bR = bL + 1, h = grid.segLength[s] / grid.segArea[s];
        if (mat.conductor) {
          out[mat.conductor.i][s] = u[bL * M];
          continue;
        }
        for (let i = 0; i < n; i++) {
          if (!mat.present[i] || mat.D[i] === 0) continue;
          const pe = (vel * h) / mat.D[i];
          const d = z[i] * (u[bR * M] - u[bL * M]) + ex[(s + 1) * n + i] - ex[s * n + i] - pe;
          const deta = u[bR * M + 1 + i] - u[bL * M + 1 + i] + (uLo[bR * M + 1 + i] - uLo[bL * M + 1 + i]) - pe;
          out[i][s] = sgFlux(mat.D[i] / h, d, deta, c[s * n + i], c[(s + 1) * n + i], this.sg)[0];
        }
        if (reg.mixing > 0) {
          // Eddy mixing (see _segmentMixing): −(D_mix/h) c̄_i (Δη_i − z_i q), c̄ the logarithmic
          // mean; between crowded ends, _mixCrowded's.
          const on = [...Array(n).keys()].filter((i) => mat.present[i] && mat.D[i] > 0);
          if (!this.nodeIdeal[s] || !this.nodeIdeal[s + 1]) {
            const { N } = this._mixBase(s, s + 1, bL, bR, on);
            on.forEach((i, p) => (out[i][s] += (reg.mixing / h) * N[p]));
            continue;
          }
          const cb = [], de = [];
          let S = 0, Q = 0;
          for (const k of on) {
            const a = c[s * n + k], b = c[(s + 1) * n + k];
            cb[k] = a === b ? a : (b - a) / Math.log(b / a);
            de[k] = u[bR * M + 1 + k] - u[bL * M + 1 + k] + (uLo[bR * M + 1 + k] - uLo[bL * M + 1 + k]);
            S += z[k] * z[k] * cb[k];
            Q += z[k] * cb[k] * de[k];
          }
          const q = S > 0 ? Q / S : 0;
          for (const i of on) out[i][s] -= (reg.mixing / h) * cb[i] * (de[i] - z[i] * q);
        }
      }
    }
    interfaces.forEach((itf, f) => {
      const gL = grid.regionEnd[f], b = this.blockOfFace[f];
      for (let i = 0; i < n; i++) if (itf.links[i].type !== 'blocked') out[i][gL] = grid.area[gL] * u[b * M + 1 + i];
    });
    return out;
  }

  // A face's one-way fluxes for each species crossing it by permeability, per area: toward +x
  // and toward −x, whose difference is its net flux (Ussing's unidirectional fluxes; GHK's two
  // terms between like solutions). null for a species crossing any other way.
  faceOneWay(f) {
    const { model, n, M, u, uLo, z } = this, itf = model.interfaces[f], grid = model.grid;
    const bL = this.blockOfFace[f] - 1, bR = bL + 2, gl = grid.regionEnd[f], gr = gl + 1;
    return itf.links.map((link, i) => {
      if (link.type !== 'permeability') return null;
      const kl = gl * n + i, kr = gr * n + i;
      const shift = this.mu0hat[kr] - this.mu0hat[kl] - Math.log(this.cRef[kr] / this.cRef[kl]);
      // (between the diffuse layers' planes, where the face has them: see _face)
      const [dL, dR] = this.diffAt[f], sL = dL < 0 ? 0 : u[this.diffList[dL].o], sR = dR < 0 ? 0 : u[this.diffList[dR].o];
      const cl = this.c[kl] * Math.exp(-z[i] * sL), cr = this.c[kr] * Math.exp(-z[i] * sR);
      const d = z[i] * (u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) + sR - sL) + shift;
      const deta = u[bR * M + 1 + i] - u[bL * M + 1 + i] + (uLo[bR * M + 1 + i] - uLo[bL * M + 1 + i]);
      const gated = link.gates ?? [];
      const P = link.P * (gated.length ? this._open(gated, (q) => u[this._gateSlot(f, q)]) : 1);
      const right = P * bernoulli(d) * cl;
      // (far up a drift, right × e^Δη is 0·∞: the backward one is then sgFlux's, whole)
      return [right, d < SG_FAR ? right * Math.exp(deta) : sgFlux(P, d, deta, cl, cr, this.sg)[3]];
    });
  }

  // Add this step's contact fluxes (at the converged state) to each stretch's intake.
  // (A BDF2 step also moves each amount by Σ v (c* − c_n), its history term.)
  _accumulateBoundaryIntake(dt, cN, gamN = this.gamOld) {
    this._assembleBookkeeping(dt);
    this.contactDEnd = { ...this.contactD };
    this.portQEnd = Float64Array.from(this.portQ);
    if (!Number.isFinite(dt)) return;
    const last = this.model.regions.length - 1;
    const { n, cOld } = this, vol = this.model.grid.vol;
    // A stretch alone in what a terminal driven by a current feeds takes in what that current
    // says, which the step kept exactly (see circuitLaws in src/lawRows.js); the flux through a
    // metal is round-off.
    const driven = new Map();
    for (const cs of this.constraints) {
      if (cs.circuit < 0 || this.terms[cs.circuit].drive.kind !== 'I' || cs.parts.length !== 1 || cs.surface.length || cs.carriers.length) continue;
      driven.set(cs.parts[0].stretch, (cs.lambda * sourceAt(this.terms[cs.circuit].drive.src, this.sourceTime, this.sourceBefore)) / (cs.parts[0].w * FARADAY));
    }
    this.stretches.forEach((st, k) => {
      if (!st.connected) return;
      let q = 0;
      if (st.regions[0] === 0) q += this.contactFlux.left[st.species];
      if (st.regions[1] === last) q -= this.contactFlux.right[st.species];
      for (const p of st.ports) q += this.portFlux[p][st.species];
      if (driven.has(k)) q = driven.get(k);
      let hist = 0;
      for (let g = st.nodes[0]; g <= st.nodes[1]; g++) hist += vol[g] * (cOld[g * n + st.species] - cN[g * n + st.species]);
      hist += this._diffuseHeld(st, this.gamOld) - this._diffuseHeld(st, gamN);
      this.boundaryIntake[k] += q * dt + hist;
    });
  }
}

// The assembly's methods, the steady solves' and the small signal's, defined with their contracts
// in src/assemble.js, src/regions.js, src/boundaries.js, src/steady.js and src/impedance.js, are
// the Solver's own.
for (const part of [Assembly, RegionTerms, BoundaryTerms, SteadySolves, SmallSignal]) {
  for (const [k, d] of Object.entries(Object.getOwnPropertyDescriptors(part.prototype))) if (k !== 'constructor') Object.defineProperty(Solver.prototype, k, d);
}
