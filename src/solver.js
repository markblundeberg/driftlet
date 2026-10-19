
// Solve the K×K complex system (Ar + i Ai) x = (last column), partial pivoting by modulus.
function complexSolve(Ar, Ai, K) {
  for (let col = 0; col < K; col++) {
    let p = col;
    for (let r = col + 1; r < K; r++) if (Math.hypot(Ar[r][col], Ai[r][col]) > Math.hypot(Ar[p][col], Ai[p][col])) p = r;
    [Ar[col], Ar[p]] = [Ar[p], Ar[col]];
    [Ai[col], Ai[p]] = [Ai[p], Ai[col]];
    const dr = Ar[col][col], di = Ai[col][col], d2 = dr * dr + di * di;
    for (let r = col + 1; r < K; r++) {
      const fr = (Ar[r][col] * dr + Ai[r][col] * di) / d2, fi = (Ai[r][col] * dr - Ar[r][col] * di) / d2;
      for (let k = col; k <= K; k++) {
        const ur = Ar[col][k], ui = Ai[col][k];
        Ar[r][k] -= fr * ur - fi * ui;
        Ai[r][k] -= fr * ui + fi * ur;
      }
    }
  }
  const re = new Float64Array(K), im = new Float64Array(K);
  for (let r = K - 1; r >= 0; r--) {
    let vr = Ar[r][K], vi = Ai[r][K];
    for (let k = r + 1; k < K; k++) {
      vr -= Ar[r][k] * re[k] - Ai[r][k] * im[k];
      vi -= Ar[r][k] * im[k] + Ai[r][k] * re[k];
    }
    const dr = Ar[r][r], di = Ai[r][r], d2 = dr * dr + di * di;
    re[r] = (vr * dr + vi * di) / d2;
    im[r] = (vi * dr - vr * di) / d2;
  }
  return { re, im };
}

// GMRES(m), restarted, for a complex system op(x) = b, with vectors as [re, im] pairs and op
// already preconditioned; x holds the starting guess and receives the answer. A cycle's own
// residual estimate can drift far from the true one (a preconditioner nearly singular at low
// frequency, or an op whose round-off the preconditioner amplifies), so each cycle is judged by
// the true residual, at the next one's start: it returns once that's below tol, and gives up when
// two cycles in a row fail to lower it, leaving x at the best iterate. Returns the true residuals (relative
// to b) at the start and at that iterate.
function gmres(op, b, x, { m = 12, restarts = 8, tol = 1e-8 } = {}) {
  const N = b[0].length;
  const cdot = (a, c) => {
    let re = 0, im = 0;
    for (let j = 0; j < N; j++) {
      re += a[0][j] * c[0][j] + a[1][j] * c[1][j];
      im += a[0][j] * c[1][j] - a[1][j] * c[0][j];
    }
    return [re, im]; // conj(a)·c
  };
  const norm = (a) => Math.sqrt(cdot(a, a)[0]);
  const bn = norm(b) || 1;
  const best = [Float64Array.from(x[0]), Float64Array.from(x[1])];
  let initial = Infinity, least = Infinity, idle = 0;
  for (let cycle = 0; ; cycle++) {
    const ax = op(x), r = [Float64Array.from(b[0], (v, j) => v - ax[0][j]), Float64Array.from(b[1], (v, j) => v - ax[1][j])];
    const beta = norm(r) / bn;
    if (cycle === 0) initial = beta;
    idle = beta < least ? 0 : idle + 1;
    if (idle === 0) {
      least = beta;
      best[0].set(x[0]);
      best[1].set(x[1]);
    }
    if (!(beta > tol) || cycle === restarts || idle === 2) break;
    const V = [[r[0].map((v) => v / (beta * bn)), r[1].map((v) => v / (beta * bn))]];
    const H = Array.from({ length: m + 1 }, () => Array.from({ length: m }, () => [0, 0]));
    const cs = [], sn = [], g = Array.from({ length: m + 1 }, () => [0, 0]);
    g[0] = [beta * bn, 0];
    let k = 0;
    for (; k < m; k++) {
      const w = op(V[k]);
      for (let i = 0; i <= k; i++) {
        const [hr, hi] = cdot(V[i], w);
        H[i][k] = [hr, hi];
        for (let j = 0; j < N; j++) {
          w[0][j] -= hr * V[i][0][j] - hi * V[i][1][j];
          w[1][j] -= hr * V[i][1][j] + hi * V[i][0][j];
        }
      }
      const hn = norm(w);
      H[k + 1][k] = [hn, 0];
      V.push([w[0].map((v) => v / (hn || 1)), w[1].map((v) => v / (hn || 1))]);
      // The earlier rotations, then a new one zeroing H[k+1][k]: x' = c x + s y, y' = −s̄ x + c y.
      for (let i = 0; i < k; i++) {
        const [xr, xi] = H[i][k], [yr, yi] = H[i + 1][k], c = cs[i], [sr, si] = sn[i];
        H[i][k] = [c * xr + sr * yr - si * yi, c * xi + sr * yi + si * yr];
        H[i + 1][k] = [-(sr * xr + si * xi) + c * yr, -(sr * xi - si * xr) + c * yi];
      }
      const [ar, ai] = H[k][k], am = Math.hypot(ar, ai), t = Math.hypot(am, hn);
      const c = t > 0 ? am / t : 1, s = am > 0 ? [(ar / am) * (hn / t), (ai / am) * (hn / t)] : [1, 0];
      cs.push(c);
      sn.push(s);
      H[k][k] = am > 0 ? [(ar / am) * t, (ai / am) * t] : [hn, 0];
      H[k + 1][k] = [0, 0];
      const [gr, gi] = g[k];
      g[k] = [c * gr, c * gi];
      g[k + 1] = [-(s[0] * gr + s[1] * gi), -(s[0] * gi - s[1] * gr)];
      if (!(Math.hypot(g[k + 1][0], g[k + 1][1]) > tol * bn)) {
        k++;
        break;
      }
    }
    // Back-substitute H y = g over the first k columns, then x += V y.
    const y = Array.from({ length: k }, () => [0, 0]);
    for (let i = k - 1; i >= 0; i--) {
      let [vr, vi] = g[i];
      for (let j = i + 1; j < k; j++) {
        const [hr, hi] = H[i][j], [yr, yi] = y[j];
        vr -= hr * yr - hi * yi;
        vi -= hr * yi + hi * yr;
      }
      const [dr, di] = H[i][i], d2 = dr * dr + di * di;
      y[i] = [(vr * dr + vi * di) / d2, (vi * dr - vr * di) / d2];
    }
    for (let i = 0; i < k; i++) {
      const [yr, yi] = y[i];
      for (let j = 0; j < N; j++) {
        x[0][j] += yr * V[i][0][j] - yi * V[i][1][j];
        x[1][j] += yr * V[i][1][j] + yi * V[i][0][j];
      }
    }
  }
  x[0].set(best[0]);
  x[1].set(best[1]);
  return { initial, residual: least };
}

// The numerical core: the Solver's state (laid out in docs/internals.md), its setup, and the
// drivers around the modules that state their contracts: the equations (src/assemble.js, with
// src/regions.js and src/boundaries.js, whose methods are the Solver's own), Newton
// (src/newton.js), time steps (src/timeSteps.js), the cold start (src/start.js) and the laws
// (src/laws.js, src/lawRows.js).

import { Assembly } from './assemble.js';
import { BlockTridiagonal, ComplexBlockTridiagonal } from './blockTridiagonal.js';
import { bernoulli, bernoulliDerivative } from './bernoulli.js';
import { BoundaryTerms } from './boundaries.js';
import { EPS0, FARADAY } from './constants.js';
import { sourceAt } from './device.js';
import { DeviceError } from './errors.js';
import { chargeOnly, conservedMoieties, feeds, findLaws, noCurrent } from './laws.js';
import { captureLaws, captureLocal } from './lawRows.js';
import { maxPotentialStep, newton as newtonFn } from './newton.js';
import { RegionTerms } from './regions.js';
import { advance as advanceFn, integrate as integrateFn, step as stepFn, violentStart } from './timeSteps.js';
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

export class Solver {
  constructor(model) {
    this.model = model;
    const { grid, species, materials, regions, contacts } = model;
    const n = species.length;
    // Slots per block: φ̂ (or D, or a conductor's segment flux J), one per species, and at a face
    // one per face reaction (its rate) and one per gate (its fraction open), as many as the
    // busiest face has.
    const nRx = model.interfaces.reduce((m, itf) => Math.max(m, itf.reactions.length + itf.gates.length), 0);
    // (and at an electrode port's nodes, one per surface species: its coverage's η)
    const nSurf = model.ports.reduce((m, port) => Math.max(m, port.surface.length), 0);
    // (and a membrane port's gates, after them: each one's fraction open)
    const nPortSlots = model.ports.reduce((m, port) => Math.max(m, port.surface.length + port.gates.length), 0);
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

    // Electrode ports' surfaces: each node's port (−1: none) and its index in the port's window;
    // the coverages there, θ [g·nSurf + s], now and at the step's start.
    this.surfPort = new Int32Array(nNodes).fill(-1);
    this.surfW = new Int32Array(nNodes);
    model.ports.forEach((port, k) => {
      if (port.surface.length > 0) port.nodes.forEach((g, w) => ((this.surfPort[g] = k), (this.surfW[g] = w)));
    });
    // Each surface species of each port, as a column of what's conserved: [port, species].
    this.surfCols = model.ports.flatMap((port, k) => port.surface.map((_, q) => [k, q]));
    this.surfColOf = model.ports.map((port, k) => this.surfCols.findIndex(([kk]) => kk === k));
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
    this.ex = new Float64Array(nNodes * n);
    this.zeta = new Float64Array(n);
    this.scratch = new Map(); // per statistics model: ζ, c, K work arrays
    this.dA = new Float64Array(M); // derivative work vectors over one block's slots
    this.dB = new Float64Array(M);
    this.jL = new Float64Array(M);
    this.jR = new Float64Array(M);
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
    // A floating terminal starts level with the first held one (no current).
    const level = this.termV[this.terms.findIndex((_, k) => !this.floating.includes(k))];
    for (const k of this.floating) this.termV[k] = this.terms[k].drive.kind === 'V' ? this.termV[k] : level;
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

    // Each reacting port's electrode area per volume at its window's nodes (a profile's mean over
    // each node's box, as for a rate constant), and its reactions' rates there, mol/(m²·s).
    this.portArea = model.ports.map((port) =>
      port.area === null ? null : Float64Array.from(port.nodes, (g) => port.area.value ?? profileMean(port.area, gx[g] - (g > 0 ? segLength[g - 1] : 0) / 2, gx[g] + (g < segLength.length ? segLength[g] : 0) / 2)),
    );
    this.portRates = model.ports.map((port) => port.reactions.map(() => new Float64Array(port.nodes.length)));

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
    for (const port of model.ports) if (port.species.some((l) => l.type === 'equilibrium')) for (const g of port.nodes) held[g] = 1;
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
        if (!port.nodes.includes(g)) return;
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
      const k = this.surfPort[g];
      if (k >= 0) for (let q = 0; q < model.ports[k].surface.length; q++) active[b * M + 1 + n + q] = 1;
      const kg = this.gatePort[g], pg = model.ports[kg];
      if (kg >= 0) for (let q = 0; q < pg.gates.length; q++) active[b * M + 1 + n + pg.surface.length + q] = 1;
    }
    model.interfaces.forEach((itf, f) => {
      const bf = this.blockOfFace[f];
      active[bf * M] = itf.phi.type === 'neutral' ? 0 : 1;
      for (let i = 0; i < n; i++) active[bf * M + 1 + i] = itf.links[i].type === 'blocked' ? 0 : 1;
      for (let k = 0; k < itf.reactions.length + itf.gates.length; k++) active[bf * M + 1 + n + k] = 1;
    });
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

  // Every stretch either reaches a contact or is kept by a law, and no law is left unkept: the
  // steady equations can be solved directly (dt = ∞), with the laws in place of their rows. Else
  // the steady solve goes by time steps, which keep every amount themselves.
  _directSteady() {
    return this.unkept.length === 0 && this.stretches.every((st) => st.connected || st.conserved);
  }

  // An electrode surface at its starting coverages at node g: η = μ°/RT + ln(θ/θ₀).
  _surfaceStart(g) {
    const { n, M, u, model } = this, surf = model.ports[this.surfPort[g]].surface, b = this.blockOfNode[g];
    const bare = 1 - surf.reduce((t, sp) => t + sp.theta0, 0);
    surf.forEach((sp, q) => (u[b * M + 1 + n + q] = sp.mu0 / model.RT + Math.log(sp.theta0 / bare)));
  }

  /** Amount (mol per unit area) of surface column j's species, Γ Σ v·a·θ, in the current state. */
  surfaceAmount(j) {
    const [k, q] = this.surfCols[j], port = this.model.ports[k], vol = this.model.grid.vol;
    let s = 0;
    port.nodes.forEach((g, w) => (s += vol[g] * this.portArea[k][w] * this.th[g * this.nSurf + q]));
    return port.surface[0].capacity * s;
  }

  /** Total amount (mol per unit area) of a stretch's species in the current state. */
  amount(stretch) {
    const { n, c, model } = this;
    const vol = model.grid.vol;
    let s = 0;
    for (let g = stretch.nodes[0]; g <= stretch.nodes[1]; g++) s += vol[g] * c[g * n + stretch.species];
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

  // After a steady solve: each reaction that alone changes some combination nothing else feeds
  // (see src/laws.js), whose total rate a steady state therefore sets to zero, but whose affinity has one
  // sign wherever it acts, so that its total can't be zero. Too slow to show, its rate lost to
  // round-off, it was left where the solve began: a slow complexation 43 mV off, an electrode
  // port's slow couple at its starting composition, each converged without a word. (A face
  // reaction stands on its affinity, so it's only caught here if that couldn't be arranged.)
  stuckReactions() {
    const { n, M, u, uLo, model } = this, grid = model.grid, S = this.stretches.length, X = S + this.surfCols.length, out = [];
    const passes = (k) => !noCurrent(this.terms[k]);
    const fed = (w) => w.some((x, k) => x !== 0 && ((k >= X && passes(this.carriers[k - X])) || (k < S && this.terms.some((t, j) => t.drive.kind === 'I' && passes(j) && feeds(j, this.stretches[k])))));
    const needed = (rx) => {
      const basis = conservedMoieties(this._structure(), { drives: true, skip: rx });
      return basis.some((w) => !fed(w) && basis.skipped.some((row) => Math.abs(row.reduce((t, x, k) => t + x * w[k], 0)) > 1e-9));
    };
    const eta = (g, slot) => u[this.blockOfNode[g] * M + slot] + uLo[this.blockOfNode[g] * M + slot];
    const report = (where, as) => {
      if (as.length === 0 || !(as.every((a) => a > 1e-6) || as.every((a) => a < -1e-6))) return;
      out.push({ where, affinity: as.reduce((m, a) => (Math.abs(a) < Math.abs(m) ? a : m)) });
    };
    model.interfaces.forEach((itf, f) => {
      const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
      itf.reactions.forEach((rx, q) => {
        if (this.affinities.some((af) => af.rx === rx) || !needed(rx)) return;
        report(`interfaces[${f}].reactions[${q}]`, [rx.part.reduce((a, p) => a - p.nu * eta(p.side ? gR : gL, 1 + p.i), rx.fixedA)]);
      });
    });
    model.reactions.forEach((mrx, x) => {
      const rx = this.rxs[x];
      if (rx.generation || !needed(mrx)) return;
      const as = [];
      for (let g = 0; g < this.nNodes; g++) if (rx.kfNode[g] > 0 || rx.srh?.[this.nodeMaterial[g]]) as.push(rx.sp.reduce((a, i, p) => a + rx.nu[p] * eta(g, 1 + i), rx.fixedA));
      report(`bulkReactions[${x}]`, as);
    });
    model.ports.forEach((port, k) => {
      port.reactions.forEach((rx, q) => {
        if (!needed(rx)) return;
        let aV = rx.fixedA;
        for (const p of rx.part) if (p.side === 1) aV -= p.nu * this.portEta(port, p.i);
        report(`ports[${k}].reactions[${q}]`, port.nodes.map((g) => rx.part.reduce((a, p) => (p.side === 1 ? a : a - p.nu * eta(g, p.side === 2 ? 1 + n + p.s : 1 + p.i)), aV)));
      });
    });
    return out;
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
    };
  }

  _restore(s) {
    this.u.set(s.u);
    this.uLo.set(s.uLo);
    this.cOld.set(s.cOld);
    this.thOld.set(s.thOld);
    this.th0Old.set(s.th0Old);
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

  /**
   * Small-signal impedance about the current (steady) state, seen at one terminal. Linearising
   * gives (J + iωM)·δx + Σ B_k δV_k = rhs for the grid's rows, and (C_k + iωC′_k)·δx +
   * (∂I_k/∂V_k)·δV_k = δI_set,k for each floating terminal's circuit row, with J the steady
   * Jacobian and M the storage part (from J(dt) = J + M/dt), B_k = ∂res/∂V_k and C_k = ∂I_k/∂x.
   * At the measured terminal, a held voltage is perturbed (rhs = −B δV) and its current read
   * from C, or a driven current is perturbed and its voltage solved for. The other terminals
   * keep their drives: held ones at AC ground, driven ones at AC open circuit.
   * Z = δV/δI with I into the device (positive real part for a passive device).
   * @param {ArrayLike<number>} frequencies Hz
   */
  impedance(frequencies, opts) {
    const was = this.allTerminals;
    this.allTerminals = true; // (it reads its terminal's ∂res/∂V and ∂I/∂x)
    try {
      return this._impedance(frequencies, opts);
    } finally {
      this.allTerminals = was;
    }
  }

  _impedance(frequencies, { profiles = false, terminal = 'right' } = {}) {
    const { M, nB, sys, VT } = this;
    const kT = this.terms.findIndex((t) => t.name === terminal);
    if (kT < 0) throw new SolverError(`impedance: no terminal named '${terminal}' (${this.terms.map((t) => t.name).join(', ')})`);
    const dT = this.terms[kT].drive;
    if (dT.R > 0) throw new SolverError(`impedance: terminal '${terminal}' has a series resistance; that belongs to the external circuit`);
    // A small signal returns through the terminals held at a voltage: with every other one driven
    // by a current (held open), none comes back, and the impedance is infinite (else round-off).
    const others = this.terms.filter((t, k) => k !== kT);
    if (others.length && others.every((t) => t.drive.kind === 'I')) {
      throw new SolverError(
        `impedance: every other terminal (${others.map((t) => t.name).join(', ')}) is driven by a current, so it is open to a small signal and no current returns through '${terminal}': ` +
          `its impedance is infinite. Measure at one of those instead (impedance(f, { terminal: '${others[0].name}' }), its voltage against '${terminal}'), or hold one of them at a voltage`,
      );
    }
    const N = nB * M, fl = this.floating, K = fl.length;
    this.sourceTime = this.time;
    this.computeConcentrations();
    this.cOld.set(this.c);
    this.thOld.set(this.th);
    this.th0Old.set(this.th0);
    this.gateList.forEach((gt, j) => (this.gateOld[j] = this.u[gt.o]));
    this.assemble(Infinity);
    this.contactDStart = { ...this.contactD };
    this.portQStart = Float64Array.from(this.portQ);
    this.assemble(Infinity);
    // The steady parts, then the storage parts from an assembly at a tiny dt.
    const J = { A: sys.A.slice(), B: sys.B.slice(), C: sys.C.slice() };
    const Bt = this.termB.map((v) => v.slice()), Ct = this.termC.map((v) => v.slice()), DI = Float64Array.from(this.termDI);
    const dts = 1e-30;
    this.assemble(dts);
    const S = {};
    for (const X of ['A', 'B', 'C']) S[X] = sys[X].map((v, k) => (v - J[X][k]) * dts);
    const Cs = this.termC.map((v, k) => v.map((x, j) => (x - Ct[k][j]) * dts));
    const DIs = this.termDI.map((x, k) => (x - DI[k]) * dts);
    // Where the measured contact's current may be read across a segment: anywhere without
    // ports; with them, only between it and the nearest window (a port's exchange takes current
    // in or out along it). Never for a port's own terminal.
    let cuts = kT < 2 ? this._currentCuts() : null;
    if (cuts && this.model.ports.length > 0) {
      const windows = this.model.ports.flatMap((p) => Array.from(p.nodes));
      const lo = Math.min(...windows), hi = Math.max(...windows);
      cuts = cuts.filter((cut) => (kT === 0 ? cut.s + 1 <= lo : cut.s >= hi));
    }
    // J·v for GMRES, exact: the assembled J holds a flux's dependence on η_L and η_R as two
    // entries, and where both are huge (an inversion layer) and v nearly uniform, J·v loses the
    // flux to round-off. So the dilute kernels' terms are kept in difference form (see
    // DifferenceTerms), with the rest from a matrix of their own, steady (Jr) and storage (Sr).
    this._assembleDifference(Infinity);
    const Jr = { A: this.sysRest.A.slice(), B: this.sysRest.B.slice(), C: this.sysRest.C.slice() };
    const { lin } = this._assembleDifference(dts);
    const Sr = {};
    for (const X of ['A', 'B', 'C']) Sr[X] = this.sysRest[X].map((v, k) => (v - Jr[X][k]) * dts);
    this.assemble(Infinity);

    const { sizes, offA, offB, offC, offX } = sys;
    const csys = new ComplexBlockTridiagonal(nB, sizes);
    const Nc = csys.size;
    const rr = new Float64Array(Nc), ri = new Float64Array(Nc);
    const yr = new Float64Array(Nc), yi = new Float64Array(Nc);
    const Xr = fl.map(() => new Float64Array(Nc)), Xi = fl.map(() => new Float64Array(Nc));
    const xr = new Float64Array(N), xi = new Float64Array(N);
    const rRes = new Float64Array(Nc), iRes = new Float64Array(Nc);
    const jr = new Float64Array(Nc), ji = new Float64Array(Nc), sr = new Float64Array(Nc), si = new Float64Array(Nc);
    const out = { f: Float64Array.from(frequencies), Z: { re: new Float64Array(frequencies.length), im: new Float64Array(frequencies.length) } };
    if (profiles) out.profiles = [];
    // (C + iωC′)·(ar + i ai) for terminal k
    const dot = (k, w, ar, ai) => {
      let re = 0, im = 0;
      const c = Ct[k], cs = Cs[k];
      for (let j = 0; j < Nc; j++) {
        re += c[j] * ar[j] - w * cs[j] * ai[j];
        im += c[j] * ai[j] + w * cs[j] * ar[j];
      }
      return [re, im];
    };
    Array.from(frequencies).forEach((f, q) => {
      const w = 2 * Math.PI * f;
      if (!(w > 0)) throw new SolverError('impedance: frequencies must be positive');
      // Rows scaled by their largest entry (and the terminal columns with them).
      const scale = new Float64Array(Nc);
      for (let blk = 0; blk < nB; blk++) {
        const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
        for (let r = 0; r < m; r++) {
          const rows = [['A', offA[blk] + r * mp, mp], ['B', offB[blk] + r * m, m], ['C', offC[blk] + r * mn, mn]];
          let mx = 0;
          for (const [X, o, len] of rows) for (let c = 0; c < len; c++) mx = Math.max(mx, Math.abs(J[X][o + c]), w * Math.abs(S[X][o + c]));
          const sc = mx > 0 ? 1 / mx : 1;
          scale[offX[blk] + r] = sc;
          for (const [X, o, len] of rows) {
            const JX = J[X], SX = S[X], Xre = csys[X + 'r'], Xim = csys[X + 'i'];
            for (let c = 0; c < len; c++) {
              Xre[o + c] = JX[o + c] * sc;
              Xim[o + c] = w * SX[o + c] * sc;
            }
          }
        }
      }
      csys.factor();
      const sx = (X, vr, out) => {
        for (let blk = 0; blk < nB; blk++) {
          const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
          for (let r = 0; r < m; r++) {
            let acc = 0;
            for (let c = 0; c < mp; c++) acc += X.A[offA[blk] + r * mp + c] * vr[offX[blk - 1] + c];
            for (let c = 0; c < m; c++) acc += X.B[offB[blk] + r * m + c] * vr[offX[blk] + c];
            for (let c = 0; c < mn; c++) acc += X.C[offC[blk] + r * mn + c] * vr[offX[blk + 1] + c];
            out[offX[blk] + r] = acc;
          }
        }
      };
      // A x = −σ ∂res/∂V_k, x = outr + i outi (σ = 1 for the response to the measured terminal,
      // −1 for a floating one's unit δV): GMRES on the system preconditioned by the factorised one
      // (whose error is concentrated in a few slow modes, which it finds in a few iterations).
      const zero = new Float64Array(Nc);
      const precondition = (ar, ai) => {
        const outR = new Float64Array(Nc), outI = new Float64Array(Nc);
        for (let j = 0; j < Nc; j++) {
          rRes[j] = ar[j] * scale[j];
          iRes[j] = ai[j] * scale[j];
        }
        csys.solve(rRes, iRes, outR, outI);
        return [outR, outI];
      };
      // (J + iωS)·v, exactly (see DifferenceTerms), into (ar, ai).
      const jv = (vr, vi, ar, ai) => {
        sx(Jr, vr, jr);
        sx(Jr, vi, ji);
        sx(Sr, vi, si);
        sx(Sr, vr, sr);
        for (let j = 0; j < Nc; j++) {
          ar[j] = jr[j] - w * si[j];
          ai[j] = ji[j] + w * sr[j];
        }
        lin.apply(vr, ar, 1, 0);
        lin.apply(vi, ar, 0, -w);
        lin.apply(vi, ai, 1, 0);
        lin.apply(vr, ai, 0, w);
      };
      const op = ([vr, vi]) => {
        const ar = new Float64Array(Nc), ai = new Float64Array(Nc);
        jv(vr, vi, ar, ai);
        return precondition(ar, ai);
      };
      const solve = (br, bi, outr, outi, k, sigma) => {
        const rhs = k >= 0 ? precondition(Bt[k].map((v) => -sigma * v), zero) : [new Float64Array(Nc), new Float64Array(Nc)];
        outr.set(rhs[0]);
        outi.set(rhs[1]);
        // (GMRES keeps its best iterate, which is never worse than the factorised solve.)
        const { residual } = gmres(op, rhs, [outr, outi]);
        accuracy = Math.max(accuracy, residual);
      };
      let accuracy = 1e-14; // (relative, of the solves: the readouts' error is weighed with it)
      // y: the response with the floating terminals held; X_k: to a unit δV_k.
      const heldT = dT.kind === 'V';
      for (let j = 0; j < Nc; j++) {
        rr[j] = heldT ? -Bt[kT][j] * scale[j] : 0;
        ri[j] = 0;
      }
      solve(rr, ri, yr, yi, heldT ? kT : -1, heldT ? 1 : 0);
      fl.forEach((k, a) => {
        for (let j = 0; j < Nc; j++) rr[j] = Bt[k][j] * scale[j];
        ri.fill(0);
        solve(rr, ri, Xr[a], Xi[a], k, -1);
      });
      // The floating terminals' circuit rows, a small complex system for their δV.
      const Ar = Array.from({ length: K }, () => new Float64Array(K + 1)), Ai = Array.from({ length: K }, () => new Float64Array(K + 1));
      fl.forEach((k, a) => {
        for (let b = 0; b < K; b++) {
          const [re, im] = dot(k, w, Xr[b], Xi[b]);
          Ar[a][b] = -re;
          Ai[a][b] = -im;
        }
        Ar[a][a] += DI[k];
        Ai[a][a] += w * DIs[k];
        const [yre, yim] = dot(k, w, yr, yi);
        Ar[a][K] = (k === kT ? 1 : 0) - yre;
        Ai[a][K] = -yim;
      });
      const dV = complexSolve(Ar, Ai, K);
      // δx = y − Σ X_k δV_k
      for (let j = 0; j < Nc; j++) {
        let re = yr[j], im = yi[j];
        for (let a = 0; a < K; a++) {
          re -= Xr[a][j] * dV.re[a] - Xi[a][j] * dV.im[a];
          im -= Xr[a][j] * dV.im[a] + Xi[a][j] * dV.re[a];
        }
        yr[j] = re;
        yi[j] = im;
        xr[this.fullOf[j]] = re;
        xi[this.fullOf[j]] = im;
      }
      // The current into terminal kT per volt, (C + iωC′)·δx + (∂I/∂V)(1 + iω′), with its error.
      // At the contact it's a sum whose terms can be far larger than itself (in a junction's
      // neutral ends, a conductance times a δη that cancels its neighbours' to 1e-17), and its
      // round-off reads as a constant conductance. Between the two contacts of a device without
      // ports, the total current (conduction and displacement) is the same through every cut, so
      // it's also read across each segment, and the reading with the least error taken. Each
      // reading's error is its terms' round-off and the solves' error in δx (to `accuracy` of its
      // largest entry) through its coefficients: a contact's or a cut's flux is a difference of
      // levels that the solved rows don't pin.
      // `parts`: the response as a sum, each part compact (r, i) and full-index (fr, fi), with the
      // absolute error of its entries.
      const read = (parts) => {
        const eps = 1.1e-16, c = Ct[kT], cs = Cs[kT];
        let Ir = DI[kT], Ii = w * DIs[kT], err = eps * (Math.abs(DI[kT]) + w * Math.abs(DIs[kT]));
        for (const p of parts) {
          const [re, im] = dot(kT, w, p.r, p.i);
          Ir += re;
          Ii += im;
          for (let j = 0; j < Nc; j++) err += (Math.abs(c[j]) + w * Math.abs(cs[j])) * (eps * Math.hypot(p.r[j], p.i[j]) + p.floor);
        }
        if (!cuts) return { Ir, Ii, err };
        const sgn = kT === 0 ? 1 : -1;
        for (const cut of cuts) {
          let re = 0, im = 0, e = 0;
          for (const p of parts) {
            for (let a = 0; a < cut.idx.length; a++) {
              const x = cut.idx[a], gk = cut.g[a], dk = w * cut.d[a];
              // (g + iωd)(xr + i xi)
              re += gk * p.fr[x] - dk * p.fi[x];
              im += gk * p.fi[x] + dk * p.fr[x];
              e += (Math.abs(gk) + Math.abs(dk)) * (eps * Math.hypot(p.fr[x], p.fi[x]) + p.floor);
            }
          }
          if (e < err) [err, Ir, Ii] = [e, sgn * re, sgn * im];
        }
        return { Ir, Ii, err };
      };
      const largest = (r, i) => {
        let m = 0;
        for (let j = 0; j < r.length; j++) m = Math.max(m, Math.hypot(r[j], i[j]));
        return m;
      };
      const full = (r, i) => {
        const fr = new Float64Array(N), fi = new Float64Array(N);
        for (let j = 0; j < Nc; j++) {
          fr[this.fullOf[j]] = r[j];
          fi[this.fullOf[j]] = i[j];
        }
        return [fr, fi];
      };
      let Zr, Zi;
      if (heldT) {
        // δI per volt: (C + iωC′)·δx + (∂I/∂V)(1 + iω′) at the terminal…
        let { Ir, Ii, err } = read([{ r: yr, i: yi, fr: xr, fi: xi, floor: accuracy * largest(yr, yi) }]);
        // …and where the solves leave that unresolved (the estimate past 1% of the current), the
        // response is split: a uniform shift of each region (z_i s on every η_i, s on φ̂), which
        // changes no flux and no concentration, so that J·e vanishes inside each region and is
        // taken at its edges; and what's left, small, solved for to GMRES's floor, its error now
        // relative to itself. An electrolyte at its open circuit needs it: the response to its
        // terminal is nearly such a shift, and its current the slope of levels uniform to 1e-14,
        // more than a double can hold. (A port's exchange sees the shift: its window's rows are
        // kept with the edges'.) So does a MOS capacitor's back contact, its channel held by a port.
        if (K === 0 && cuts && err > 1e-2 * Math.hypot(Ir, Ii)) {
          const [er, ei] = this._regionShift(yr, yi);
          const tr = new Float64Array(Nc), ti = new Float64Array(Nc);
          jv(er, ei, tr, ti);
          const edge = this._regionEdgeRows();
          const ar = new Float64Array(Nc), ai = new Float64Array(Nc);
          for (let j = 0; j < Nc; j++) {
            ar[j] = -Bt[kT][j] - (edge[j] ? tr[j] : 0);
            ai[j] = edge[j] ? -ti[j] : 0;
          }
          const rhs = precondition(ar, ai);
          const dr = yr.map((v, j) => v - er[j]), di = yi.map((v, j) => v - ei[j]);
          const { residual } = gmres(op, rhs, [dr, di], { tol: 1e-11 });
          const [efr, efi] = full(er, ei), [dfr, dfi] = full(dr, di);
          const r2 = read([
            { r: er, i: ei, fr: efr, fi: efi, floor: 0 },
            { r: dr, i: di, fr: dfr, fi: dfi, floor: Math.max(1e-14, residual) * largest(dr, di) },
          ]);
          if (r2.err < err) ({ Ir, Ii, err } = r2);
          if (profiles) {
            for (let j = 0; j < N; j++) {
              xr[j] = efr[j] + dfr[j];
              xi[j] = efi[j] + dfi[j];
            }
          }
        }
        const d2 = Ir * Ir + Ii * Ii; // Z = 1/δI
        Zr = Ir / d2;
        Zi = -Ii / d2;
      } else {
        // δV per unit current
        const a = fl.indexOf(kT);
        Zr = dV.re[a];
        Zi = dV.im[a];
      }
      out.Z.re[q] = Zr;
      out.Z.im[q] = Zi;
      if (profiles) out.profiles.push(this._smallSignalProfiles(xr, xi));
    });
    return out;
  }

  // A uniform shift of each region, fitted to a small-signal response (compact, complex): φ̂ by
  // s_r and each species' η by z_i s_r, with s_r the region's mean δφ̂ (or, where φ̂ isn't an
  // unknown, its carriers' mean δη/z). Coverages and flux unknowns aren't shifted.
  _regionShift(yr, yi) {
    const { n, M, z, loc, model } = this, R = this.rix, { grid } = model;
    const er = new Float64Array(yr.length), ei = new Float64Array(yr.length);
    for (let r = 0; r < model.regions.length; r++) {
      let sr = 0, si = 0, k = 0;
      const g0 = grid.regionStart[r], g1 = grid.regionEnd[r];
      for (let g = g0; g <= g1; g++) {
        const o = this.blockOfNode[g] * M;
        if (loc[o] >= 0) {
          sr += yr[R[o]];
          si += yi[R[o]];
          k++;
        }
      }
      if (k === 0) {
        for (let g = g0; g <= g1; g++) {
          const o = this.blockOfNode[g] * M;
          for (let i = 0; i < n; i++) {
            if (z[i] === 0 || loc[o + 1 + i] < 0) continue;
            sr += yr[R[o + 1 + i]] / z[i];
            si += yi[R[o + 1 + i]] / z[i];
            k++;
          }
        }
      }
      if (k === 0) continue;
      sr /= k;
      si /= k;
      for (let g = g0; g <= g1; g++) {
        const o = this.blockOfNode[g] * M;
        if (loc[o] >= 0) {
          er[R[o]] = sr;
          ei[R[o]] = si;
        }
        for (let i = 0; i < n; i++) {
          if (loc[o + 1 + i] < 0) continue;
          er[R[o + 1 + i]] = z[i] * sr;
          ei[R[o + 1 + i]] = z[i] * si;
        }
      }
    }
    return [er, ei];
  }

  // Compact rows a region's uniform shift can touch: faces, the nodes beside them, the two ends
  // and the ports' windows (whatever holds a level from outside). Elsewhere inside a region,
  // every row is unchanged by it.
  _regionEdgeRows() {
    const { M, loc, nB } = this, R = this.rix;
    const edge = new Uint8Array(this.sys.size);
    const mark = (b) => {
      if (b < 0 || b >= nB) return;
      for (let r = 0; r < M; r++) if (loc[b * M + r] >= 0) edge[R[b * M + r]] = 1;
    };
    mark(0);
    mark(nB - 1);
    for (let f = 0; f < this.nFaces; f++) {
      const b = this.blockOfFace[f];
      mark(b - 1);
      mark(b);
      mark(b + 1);
    }
    for (const port of this.model.ports) for (const g of port.nodes) mark(this.blockOfNode[g]);
    return edge;
  }

  // The cuts where the impedance may read a small-signal current. Each is the linearised total
  // current toward +x through one segment s of a region, Σ (g + iω d)·δu over the unknowns idx:
  // the Scharfetter–Gummel fluxes' derivatives (as in _segmentsDilute, or _segmentConcentrated
  // for non-ideal statistics) and the displacement's. Without ports it's the same through every cut; with them, only between a
  // contact and the nearest port's window (see _impedance).
  _currentCuts() {
    const { model, n, M, u, uLo, c, z } = this, F = FARADAY;
    const { grid, materials, regions } = model;
    const cuts = [];
    for (let r = 0; r < regions.length; r++) {
      const reg = regions[r], mat = materials[reg.material];
      if (mat.conductor || reg.mixing > 0) continue;
      const vel = reg.velocity;
      for (let s = grid.regionStart[r]; s < grid.regionEnd[r]; s++) {
        const bL = s + r, bR = bL + 1, h = grid.segLength[s] / grid.segArea[s];
        const terms = new Map();
        const add = (x, g, d) => {
          const t = terms.get(x) || [0, 0];
          t[0] += g;
          t[1] += d;
          terms.set(x, t);
        };
        const phiL = u[bL * M], phiR = u[bR * M], hasPhi = !this.phiUndefined[s];
        for (let i = 0; i < n; i++) {
          if (z[i] === 0 || !mat.present[i] || mat.D[i] === 0) continue;
          if (!mat.ideal) {
            this._fluxConcentrated(s, bL, bR, mat, h, vel, i);
            for (let k = 0; k < M; k++) {
              if (k === 0 && !hasPhi) continue;
              if (this.jL[k] !== 0) add(bL * M + k, z[i] * F * this.jL[k], 0);
              if (this.jR[k] !== 0) add(bR * M + k, z[i] * F * this.jR[k], 0);
            }
            continue;
          }
          const zi = z[i], r1 = 1 + i, cL = c[s * n + i], g = mat.D[i] / h;
          const pe = (vel * h) / mat.D[i];
          const d = zi * (phiR - phiL) - pe;
          const E = Math.expm1(u[bR * M + r1] - u[bL * M + r1] + (uLo[bR * M + r1] - uLo[bL * M + r1]) - pe);
          const gBc = g * bernoulli(d) * cL;
          const dNdd = -g * bernoulliDerivative(d) * cL * E;
          add(bL * M + r1, zi * F * gBc, 0);
          add(bR * M + r1, -zi * F * gBc * (E + 1), 0);
          if (hasPhi) {
            add(bL * M, zi * F * (-zi * dNdd + zi * gBc * E), 0);
            add(bR * M, zi * F * zi * dNdd, 0);
          }
        }
        // The displacement toward +x, ε(φ_L − φ_R)/h (in φ̂ = φ/V_T).
        const k = hasPhi ? (mat.epsr * EPS0 * this.VT) / h : 0;
        if (k > 0) {
          add(bL * M, 0, k);
          add(bR * M, 0, -k);
        }
        if (terms.size === 0) continue;
        const idx = [...terms.keys()];
        cuts.push({ s, idx, g: idx.map((x) => terms.get(x)[0]), d: idx.map((x) => terms.get(x)[1]) });
      }
    }
    return cuts;
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
          out[i][s] = -(mat.D[i] / h) * bernoulli(d) * c[s * n + i] * Math.expm1(deta);
        }
        if (reg.mixing > 0) {
          // Eddy mixing (see _segmentMixing): −(D_mix/h) c̄_i (Δη_i − z_i q), c̄ the logarithmic mean.
          const on = [...Array(n).keys()].filter((i) => mat.present[i] && mat.D[i] > 0);
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
      const d = z[i] * (u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M])) + shift;
      const deta = u[bR * M + 1 + i] - u[bL * M + 1 + i] + (uLo[bR * M + 1 + i] - uLo[bL * M + 1 + i]);
      const gated = link.gates ?? [];
      const P = link.P * (gated.length ? this._open(gated, (q) => u[this._gateSlot(f, q)]) : 1);
      const right = P * bernoulli(d) * this.c[kl];
      return [right, right * Math.exp(deta)];
    });
  }

  // Complex profiles of δφ (V), δμ̄ (J/mol) and δc (mol/m³) from a small-signal solution.
  _smallSignalProfiles(xr, xi) {
    const { n, M, z, VT, model } = this;
    const nN = this.nNodes, RT = model.RT;
    const pair = () => ({ re: new Float64Array(nN), im: new Float64Array(nN) });
    const phi = pair(), mu = {}, c = {};
    for (const sp of model.species) {
      mu[sp.name] = pair();
      c[sp.name] = pair();
    }
    for (let g = 0; g < nN; g++) {
      const b = this.blockOfNode[g];
      const undef = this.phiUndefined[g];
      phi.re[g] = undef ? NaN : VT * xr[b * M];
      phi.im[g] = undef ? NaN : VT * xi[b * M];
      for (let i = 0; i < n; i++) {
        const name = model.species[i].name;
        if (!this.present[g * n + i]) {
          mu[name].re[g] = mu[name].im[g] = c[name].re[g] = c[name].im[g] = NaN;
          continue;
        }
        mu[name].re[g] = RT * xr[b * M + 1 + i];
        mu[name].im[g] = RT * xi[b * M + 1 + i];
        if (this.nodeConductor[g] >= 0) {
          // a metal has a Fermi level but no concentration (as a solution reports it)
          c[name].re[g] = c[name].im[g] = NaN;
          continue;
        }
        // δc_i = Σ_j K_ij (δη_j − z_j δφ̂)
        let sr = 0, si = 0;
        for (let j = 0; j < n; j++) {
          const Kij = this.nodeIdeal[g] ? (i === j ? this.c[g * n + i] : 0) : this.K[g * n * n + i * n + j];
          if (Kij === 0) continue;
          sr += Kij * (xr[b * M + 1 + j] - z[j] * xr[b * M]);
          si += Kij * (xi[b * M + 1 + j] - z[j] * xi[b * M]);
        }
        c[name].re[g] = sr;
        c[name].im[g] = si;
      }
    }
    return { phi, mu, c };
  }

  // Add this step's contact fluxes (at the converged state) to each stretch's intake.
  // (A BDF2 step also moves each amount by Σ v (c* − c_n), its history term.)
  _accumulateBoundaryIntake(dt, cN) {
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
      this.boundaryIntake[k] += q * dt + hist;
    });
  }

  // Slowest diffusion time across the device, used to size "giant" steps.
  slowestTime() {
    const { model } = this;
    let dMin = Infinity;
    for (const mat of model.materials) {
      for (let i = 0; i < this.n; i++) if (mat.present[i] && mat.D[i] > 0) dMin = Math.min(dMin, mat.D[i]);
    }
    const L = model.grid.length;
    return Number.isFinite(dMin) ? (L * L) / dMin : 1;
  }

  // The fastest diffusion time: across the smallest cell, at the largest D.
  fastestTime() {
    let dMax = 0;
    for (const mat of this.model.materials) {
      for (let i = 0; i < this.n; i++) if (mat.present[i]) dMax = Math.max(dMax, mat.D[i]);
    }
    let h = Infinity;
    for (const s of this.model.grid.segLength) if (s > 0) h = Math.min(h, s);
    return dMax > 0 && Number.isFinite(h) ? (h * h) / dMax : this.slowestTime();
  }

  /**
   * Steady state: where a transient from the present state ends, found as the solution of the
   * steady equations (dt = ∞, no storage term), with no slow modes to wait out, however slow the
   * physics (e.g. exponentially scarce minority carriers). What the storage term held is held by
   * the steady system's own rows instead (see src/laws.js): each conserved amount in place of a balance
   * row that says nothing new, levels that carry nothing held flat, islands' summed balances.
   * Where Newton from here fails: continuation along the drives, then the pseudo-transient ramp.
   * The clock is not advanced, and open-system conservation bookkeeping restarts here.
   */
  // Steady state from the present one; `atSteady` says whether the state is one (until a step or
  // a change of drive).
  solveSteady(opts = {}) {
    this.violent = this.history.length === 0 && !this.stepped ? violentStart(this) : []; // (for a failure to say: see src/timeSteps.js)
    this.unreached = null; // (a driven current no held voltage reaches: see _floatingContinuation)
    this.charging = null; // (a driven current that only charges: see _solveSteadyAll)
    captureLaws(this);
    this.steady = true; // (the steady system's rows; sources at the present time)
    this.steadyAt = this.time; // (which a pseudo-transient's steps don't move on)
    try {
      return this._solveSteadyAll(opts);
    } finally {
      this.steady = false;
    }
  }

  _solveSteadyAll(opts) {
    // A terminal driven by a current that can only store what comes in (a battery's cathode
    // behind it, closed but for that terminal; a gate, or a capacitance with no species through
    // it) has no steady state but at no current, where it keeps what it holds.
    const { species, regions } = this.model;
    for (const k of this.floating) {
      const I = sourceAt(this.terms[k].drive.src, this.time);
      if (this.terms[k].drive.kind !== 'I' || I === 0) continue;
      const laws = this.constraints.filter((cs) => cs.terminals.includes(k));
      if (laws.length === 0 && !chargeOnly(this, k)) continue;
      const what = laws.flatMap((cs) => cs.parts.map(({ stretch }) => {
        const st = this.stretches[stretch], [r0, r1] = st.regions;
        return `${species[st.species].name} in ${regions[r0].name}${r1 > r0 ? `…${regions[r1].name}` : ''}`;
      }));
      this.charging = { terminal: this.terms[k].name, current: I, what };
      return { converged: false, steps: 0, iterations: 0, history: [] };
    }
    const r = this._steadyFromHere(opts);
    this.atSteady = r.converged;
    return r;
  }

  // Steady state from the present one, in three stages: Newton from here; failing that, from a
  // solved state along the drives (the light ramped up, a current's terminal held at the voltages
  // that pass it, the voltages ramped from where they were last solved); failing that, the
  // pseudo-transient ramp, through time.
  _steadyFromHere(opts) {
    this.conditioning = null;
    this.computeConcentrations();
    captureLocal(this); // immobile combinations keep, node by node, what they hold now
    this.sourceTime = this.time;
    // Continuation applies where the right terminal is held at a voltage, against a held left one.
    const [dl, dr] = [this.terms[0].drive, this.terms[1].drive];
    const held = (d) => d.kind === 'V' && !(d.R > 0);
    const target = sourceAt(dr.src, this.time), level = sourceAt(dl.src, this.time);
    const allowed = opts.continuation !== false;
    const canContinue = allowed && held(dl) && held(dr) && target !== level;
    const direct = this._directSteady();
    const u0 = Float64Array.from(this.u), u0Lo = Float64Array.from(this.uLo), v0 = Float64Array.from(this.termV);
    const restart = () => {
      this.u.set(u0);
      this.uLo.set(u0Lo);
      this.termV.set(v0);
      this.computeConcentrations();
    };
    const solved = (q) => {
      if (q.converged) {
        this.solvedV = [level, target];
        this.stepped = false;
      }
      return q;
    };
    const after = (q, r) => ({ ...q, steps: q.steps + r.steps, iterations: q.iterations + r.iterations });

    // 1. Newton from here. Where a continuation could follow, one direct attempt only, rather
    // than the pseudo-transient ramp straight away.
    const quick = (canContinue || this.hasGeneration) && direct && allowed;
    let r = this._steadySteps(quick ? { ...opts, maxSteps: 1 } : opts);
    if (r.converged) return solved(r);
    // A transient's state can be a worse start than a cold one: 1 ns after a lit MOS capacitor's
    // gate jumps, its surface is in deep depletion, minority carriers ten orders below where they
    // end, and Newton fails from there where it converges from the layout for the drives (and
    // the stages below took 9000 iterations from there, 2000 from the layout). What a steady
    // solve keeps it took from the state already, so from here on it goes as a cold device would.
    if (this.stepped && direct && allowed) {
      this._refreshSources();
      this.initFromComposition();
      u0.set(this.u);
      u0Lo.set(this.uLo);
      const q = this._steadySteps(quick ? { ...opts, maxSteps: 1 } : opts);
      if (q.converged) return solved(after(q, r));
      restart();
      r = after(q, r);
    }

    // 2. Continuation. Generation (e.g. light) holding the device far from equilibrium: ramp it up
    // from nearly nothing, each solve warm from the last.
    let dimmer = false;
    if (this.hasGeneration && allowed) {
      restart();
      this.dimmer = false;
      const g = this._generationContinuation(opts, r);
      dimmer = this.dimmer;
      if (g.converged) return solved(g);
      r = g;
    }
    // A terminal driven by a current (open circuit, say): hold it at a voltage instead, march the
    // voltage until the current crosses its target, and float it from there.
    if (allowed && this.floating.some((k) => this.terms[k].drive.kind === 'I' && !this._bound(k))) {
      restart();
      const q = this._floatingContinuation(opts, r);
      if (q.converged) return q;
      r = q;
    }
    if (canContinue) {
      // The right terminal's voltage ramped to its target: from the contacts' voltages at the last
      // converged solve when the state is that solution (a warm start), the left one first where
      // it moved, the right held where it was; otherwise from level terminals, consistent with a
      // cold start.
      restart();
      const was = this.solvedV;
      let c;
      if (was) {
        c = was[0] === level ? { ...r, converged: true, ramped: false } : this._continuation(opts, 0, was[0], level, r, true, was[1]);
        // (Ramping the right from level terminals, solve there first: a sweep's start.)
        if (c.converged) c = this._continuation(opts, 1, was[1], target, c, c.ramped || was[1] !== level);
      } else c = this._continuation(opts, 1, level, target, r, false);
      if (c.converged) return solved(c);
      // 3. The pseudo-transient ramp, from the start.
      restart();
      return solved(after(this._steadySteps(opts), c));
    }
    // 3. The pseudo-transient ramp (where the first attempt wasn't it already): from the dimmer
    // light's solution where the light's continuation got anywhere (a lit floating base charges
    // from there as it would in time), else from the start.
    if (!quick) return r;
    if (!dimmer) restart();
    else this.computeConcentrations();
    return after(this._steadySteps(opts), r);
  }

  // Whether terminal k, driven by a current, feeds a combination that it alone can change, or
  // only charges a capacitance (so that its steady states are one per amount or charge, not one
  // per voltage).
  _bound(k) {
    return this.constraints.some((cs) => cs.terminals.includes(k)) || chargeOnly(this, k);
  }

  // A cold steady solve with a terminal driven by a current, which can fail where the same state
  // is easy warm (a solar cell's open circuit, from cold on a fine grid): the terminal held at a
  // voltage instead, whose steady solves have their own continuation, the voltage marched from
  // the start's until the current crosses its target, the crossing narrowed by bisection, and
  // the terminal floated again from beside it.
  _floatingContinuation(opts, r) {
    const k = this.floating.find((j) => this.terms[j].drive.kind === 'I' && !this._bound(j));
    const drive = this.terms[k].drive, target = sourceAt(drive.src, this.time);
    let steps = r.steps ?? 0, iterations = r.iterations ?? 0;
    const history = (r.history ?? []).slice();
    const hold = (V) => {
      this.terms[k].drive = { kind: 'V', src: { value: V }, R: 0 };
      this.redrive();
      const q = this._steadyFromHere(opts);
      steps += q.steps ?? 0;
      iterations += q.iterations ?? 0;
      history.push({ held: V, converged: q.converged });
      if (q.converged) seen.push([V, this.termI[k]]);
      return q.converged ? this.termI[k] - target : NaN;
    };
    const fail = () => ({ converged: false, steps, iterations, history });
    const seen = []; // [V, I] of each held solve that converged
    let Va = this.termV[k], fa, Vb, fb;
    try {
      fa = hold(Va);
      if (!Number.isFinite(fa)) {
        Va = 0;
        fa = hold(Va);
      }
      if (!Number.isFinite(fa)) return fail();
      // March, the step growing: up where the current is below its target, since raising a
      // terminal's voltage raises the current into a passive device (the other way, from the
      // start again, if 20 V brings no crossing).
      const start = Va, f0 = fa;
      [Vb, fb] = [Va, fa];
      for (const dir of [fa < 0 ? 1 : -1, fa < 0 ? -1 : 1]) {
        [Va, fa] = [start, f0];
        if (dir !== (f0 < 0 ? 1 : -1)) hold(start); // (back to the start, for a warm march)
        let dV = 0.02;
        Vb = Va + dir * dV;
        fb = hold(Vb);
        while (Number.isFinite(fb) && Math.sign(fb) === Math.sign(fa) && fa !== 0 && Math.abs(Vb - start) < 20) {
          [Va, fa] = [Vb, fb];
          dV *= 1.5;
          Vb = Va + dir * dV;
          fb = hold(Vb);
        }
        if (Number.isFinite(fb) && (Math.sign(fb) !== Math.sign(fa) || fa === 0)) break;
      }
      if (!Number.isFinite(fb) || (Math.sign(fb) === Math.sign(fa) && fa !== 0)) {
        // What the held voltages passed, at the two extremes reached.
        if (seen.length > 1) {
          seen.sort((a, b) => a[0] - b[0]);
          this.unreached = { terminal: this.terms[k].name, target, low: seen[0], high: seen.at(-1) };
        }
        return fail();
      }
      // Bisect (each held solve warm from the last) to a tenth of a millivolt.
      for (let it = 0; it < 40 && Math.abs(Vb - Va) > 1e-4; it++) {
        const Vm = (Va + Vb) / 2, fm = hold(Vm);
        if (!Number.isFinite(fm)) return fail();
        if (Math.sign(fm) === Math.sign(fa)) [Va, fa] = [Vm, fm];
        else [Vb, fb] = [Vm, fm];
      }
      hold(Math.abs(fa) < Math.abs(fb) ? Va : Vb);
    } finally {
      this.terms[k].drive = drive;
      this.redrive();
    }
    const q = this._steadySteps(opts);
    steps += q.steps ?? 0;
    iterations += q.iterations ?? 0;
    history.push(...(q.history ?? []));
    if (q.converged) return { ...q, steps, iterations, history };
    // Floated, the system can lose what held it doesn't (a closed redox cell driven near its
    // limit: its bordered row 33 digits short), but a current-driven steady state is the held
    // one at the voltage that passes the target: found by held solves, by regula falsi
    // (Illinois) within the bracket, to 1e-10 of the target, or at open circuit of the largest
    // current the march saw (a target of 0 is never met exactly).
    const scale = Math.max(Math.abs(target), ...seen.map(([, I]) => Math.abs(I)));
    try {
      let side = 0;
      for (let it = 0; it < 60; it++) {
        const Vm = Va - (fa * (Vb - Va)) / (fb - fa), fm = hold(Vm);
        if (!Number.isFinite(fm)) return fail();
        if (Math.abs(fm) <= 1e-10 * scale) return { converged: true, steps, iterations, history, residual: Math.abs(fm) };
        if (Math.sign(fm) === Math.sign(fb)) {
          [Vb, fb] = [Vm, fm];
          if (side === -1) fa /= 2;
          side = -1;
        } else {
          [Va, fa] = [Vm, fm];
          if (side === 1) fb /= 2;
          side = 1;
        }
      }
      return fail();
    } finally {
      this.terms[k].drive = drive;
      this.redrive();
    }
  }

  // Steady solves with the generation reactions' rates scaled from 1e-12 up to 1, ×100 a step
  // while they converge, smaller steps where they don't.
  _generationContinuation(opts, r) {
    const sub = this._directSteady() ? { ...opts, maxSteps: 1 } : opts;
    let s = 1e-30, factor = 100, steps = r.steps, iterations = r.iterations;
    const history = r.history.slice();
    try {
      this.generationScale = s;
      // (nearly dark: a cold start can need the full ramp, as a dark device's does)
      let q = this._steadySteps(sub);
      steps += q.steps;
      iterations += q.iterations;
      if (!q.converged) return { converged: false, steps, iterations, history };
      this.dimmer = true; // (the state is now a solution under dimmer light: a better start than cold)
      while (s < 1) {
        const next = Math.min(1, s * factor);
        const u1 = Float64Array.from(this.u), u1Lo = Float64Array.from(this.uLo), v1 = Float64Array.from(this.termV);
        this.generationScale = next;
        q = this._steadySteps(sub);
        steps += q.steps;
        iterations += q.iterations;
        history.push({ generation: next, converged: q.converged });
        if (q.converged) {
          s = next;
          factor = Math.min(100, factor * factor);
        } else {
          this.u.set(u1);
          this.uLo.set(u1Lo);
          this.termV.set(v1);
          this.computeConcentrations();
          factor = Math.sqrt(factor);
          if (factor < 1.5) return { converged: false, steps, iterations, history }; // (crawling: time steps do better)
        }
      }
      return { converged: true, steps, iterations, history };
    } finally {
      this.generationScale = 1;
    }
  }

  // Terminal k (a contact) ramped from `level` to `target` through steady solves, each warm from
  // the last (`hold`: the other contact's voltage meanwhile, if not its own). A warm ramp that has
  // nowhere to go converges only if what it started from did (`ramped`): the state it carries
  // must be some ramp's solution, never just the one a failed solve left.
  _continuation(opts, k, level, target, r, warm, hold) {
    const sub = this._directSteady() ? { ...opts, maxSteps: 1 } : opts;
    let V = level, dV = (target - level) / 8, steps = r.steps, iterations = r.iterations;
    const history = r.history.slice();
    if (warm && V === target) return { ...r, converged: r.ramped === true, steps, iterations, history, ramped: r.ramped };
    try {
      this.sourceOverride.set(k, V);
      if (hold !== undefined) this.sourceOverride.set(1 - k, hold);
      let q;
      if (!warm) {
        q = this._steadySteps(opts);
        if (!q.converged) {
          // A cold device's start was laid out for the target (a set() before the first solve
          // builds it there): lay it out again, level.
          steps += q.steps;
          iterations += q.iterations;
          this._refreshSources();
          this.initFromComposition();
          q = this._steadySteps(opts);
        }
        // (lit, the level start may need the light ramped up too)
        if (!q.converged && this.hasGeneration) q = this._generationContinuation(opts, q);
        steps += q.steps;
        iterations += q.iterations;
        if (!q.converged) return { ...r, steps, iterations };
      }
      while (V !== target) {
        const next = Math.abs(target - V) <= Math.abs(dV) ? target : V + dV;
        const u1 = Float64Array.from(this.u), u1Lo = Float64Array.from(this.uLo), v1 = Float64Array.from(this.termV);
        this.sourceOverride.set(k, next);
        q = this._steadySteps(sub);
        steps += q.steps;
        iterations += q.iterations;
        history.push({ continuation: next, converged: q.converged });
        if (q.converged) {
          V = next;
          dV *= 1.5;
        } else {
          this.u.set(u1);
          this.uLo.set(u1Lo);
          this.termV.set(v1);
          this.computeConcentrations();
          dV /= 4;
          if (Math.abs(dV) < 1e-6 * Math.abs(target - level)) return { converged: false, steps, iterations, history };
        }
      }
      return { converged: true, steps, iterations, history, ramped: true };
    } finally {
      this.sourceOverride.delete(k);
      if (hold !== undefined) this.sourceOverride.delete(1 - k);
      this._refreshSources();
    }
  }

  _steadySteps({ maxSteps = 80, tol = 1e-11 } = {}) {
    const time = this.time;
    const tau = this.slowestTime();
    const direct = this._directSteady();
    const giant = direct ? Infinity : 1e6 * tau;
    let dt = giant, grow = 10; // (after a failure, dt grows by less, then back up to ×10)
    let triedDirect = 0;
    let totalIter = 0, steps = 0, converged = false;
    const history = [];
    while (steps < maxSteps) {
      steps++;
      // Restore conserved amounts exactly before each huge step; the step then re-solves, so
      // the final state satisfies every equation.
      if (!direct && dt >= giant) this._renormalizeSpectators();
      let r = this.step(dt);
      totalIter += r.iterations;
      history.push({ dt, converged: r.converged, iterations: r.iterations });
      if (!r.converged && dt === Infinity && steps === 1) {
        // A large change can overshoot Newton's usual damping: retry with a tighter limit.
        r = this.step(dt, { maxStep: 3, maxIter: 80 });
        totalIter += r.iterations;
        history.push({ dt, converged: r.converged, iterations: r.iterations, maxStep: 3 });
      }
      if (!r.converged) {
        // Pseudo-transient continuation: down from a slow time scale until a step converges,
        // as far as the fastest (slow ions mustn't stop it short of what a cold start needs).
        dt = dt >= giant ? tau * 1e-6 : dt / 4;
        grow = 2;
        if (dt < Math.min(tau * 1e-15, this.fastestTime() * 1e-2)) break;
        continue;
      }
      if (dt === Infinity) {
        converged = true; // the steady equations themselves were solved
        break;
      }
      // Past the slowest diffusion time, the direct solve may already reach from here: try it
      // once a decade (the state kept if it fails), rather than step on to the giant dt through
      // steps that can fail for slow ions' sake.
      if (direct && dt < giant && dt >= tau && dt >= 10 * triedDirect) {
        triedDirect = dt;
        const snap = this._snapshot();
        const d = this.step(Infinity);
        totalIter += d.iterations;
        history.push({ dt: Infinity, converged: d.converged, iterations: d.iterations });
        if (d.converged) {
          converged = true;
          break;
        }
        this._restore(snap);
      }
      if (dt >= giant && maxPotentialStep(this, this._diff()) < tol) {
        converged = true;
        if (dt > giant) {
          // Finish at the base giant step, where the storage term pins amounts tightly.
          this._renormalizeSpectators();
          const f = this.step(giant);
          totalIter += f.iterations;
          converged = f.converged;
        }
        break;
      }
      if (dt < giant) {
        dt *= grow; // ramping up after a failure
        grow = Math.min(10, grow * 1.5);
        if (dt > 1e6 * tau) dt = giant; // then the direct steady solve (or the giant step)
      } else if (dt < 1e6 * giant) {
        dt *= 10; // conserved amounts present and still moving: let dt keep growing (capped)
      }
    }
    this.time = time;
    this.history = []; // a steady solve isn't a trajectory: transients restart from here
    this.dtNext = undefined;
    if (converged) {
      this.stretches.forEach((st, k) => {
        if (!st.connected) return;
        this.referenceAmounts[k] = this.amount(st);
        this.boundaryIntake[k] = 0;
      });
    }
    return { converged, steps, iterations: totalIter, history };
  }

  // Huge steps pin each conserved amount only through a tiny storage term, so round-off can let
  // it creep. Shift each spectator's level uniformly to restore its amount exactly (in one step
  // for ideal statistics, where c ∝ e^η; by Newton on the shift otherwise).
  _renormalizeSpectators() {
    const { n, M, u } = this;
    let changed = false;
    this.stretches.forEach((st, k) => {
      // A floating conductor holds its charge on its faces, not as an amount: the huge steps keep
      // it, and there's nothing here to shift.
      if (!st.spectator || this.model.materials[this.model.regions[st.regions[0]].material].conductor) return;
      const want = this.referenceAmounts[k];
      const ideal = this.nodeIdeal.subarray(st.nodes[0], st.nodes[1] + 1).every((v) => v === 1);
      for (let it = 0; it < (ideal ? 1 : 30); it++) {
        const now = this.amount(st);
        if (!(now > 0) || now === want) return;
        let shift;
        if (ideal) shift = Math.log(want / now);
        else {
          let dA = 0;
          for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
            const dc = this.nodeIdeal[g] ? this.c[g * n + st.species] : this.K[g * n * n + st.species * (n + 1)];
            dA += this.model.grid.vol[g] * dc;
          }
          shift = (want - now) / dA;
          if (!(Math.abs(shift) > 1e-16)) return;
        }
        for (let g = st.nodes[0]; g <= st.nodes[1]; g++) u[this.blockOfNode[g] * M + 1 + st.species] += shift;
        changed = true;
        if (!ideal) this.computeConcentrations();
      }
    });
    if (changed) this.computeConcentrations();
  }

  _diff() {
    const d = this.delta; // reuse as scratch
    const { fullOf } = this;
    for (let j = 0; j < fullOf.length; j++) d[j] = this.u[fullOf[j]] - this.uPrev[fullOf[j]];
    return d;
  }
}

// The assembly's methods, defined with their contracts in src/assemble.js, src/regions.js and
// src/boundaries.js, are the Solver's own.
for (const part of [Assembly, RegionTerms, BoundaryTerms]) {
  for (const [k, d] of Object.entries(Object.getOwnPropertyDescriptors(part.prototype))) if (k !== 'constructor') Object.defineProperty(Solver.prototype, k, d);
}
