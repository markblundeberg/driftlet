// Butler–Volmer factor g(a) = e^{αa} − e^{−(1−α)a} = e^{−(1−α)a}·expm1(a) (precise near a = 0),
// and its derivative g′(a) = α e^{αa} + (1−α) e^{−(1−α)a}.
function bvFactor(a, alpha) {
  const e = Math.exp(-(1 - alpha) * a);
  // (far from equilibrium, the plain difference: the product would give 0·∞ there)
  const g = Math.abs(a) > 50 ? Math.exp(alpha * a) - e : e * Math.expm1(a);
  return { g, gp: alpha * Math.exp(alpha * a) + (1 - alpha) * e };
}

// GMRES(m), restarted, for a real system op(x) = b, op already preconditioned; x holds the
// starting guess and receives the answer. Returns whether it converged.
function gmresReal(op, b, x, { m = 10, restarts = 3, tol = 1e-12 } = {}) {
  const N = b.length;
  const dot = (a, c) => {
    let t = 0;
    for (let j = 0; j < N; j++) t += a[j] * c[j];
    return t;
  };
  const bn = Math.sqrt(dot(b, b)) || 1;
  for (let cycle = 0; cycle < restarts; cycle++) {
    const ax = op(x), r = Float64Array.from(b, (v, j) => v - ax[j]);
    const beta = Math.sqrt(dot(r, r));
    if (!(beta > tol * bn)) return true;
    const V = [r.map((v) => v / beta)], H = Array.from({ length: m + 1 }, () => new Float64Array(m));
    const cs = [], sn = [], g = new Float64Array(m + 1);
    g[0] = beta;
    let k = 0, converged = false;
    for (; k < m; k++) {
      const w = op(V[k]);
      for (let i = 0; i <= k; i++) {
        const h = dot(V[i], w);
        H[i][k] = h;
        for (let j = 0; j < N; j++) w[j] -= h * V[i][j];
      }
      const hn = Math.sqrt(dot(w, w));
      V.push(w.map((v) => v / (hn || 1)));
      for (let i = 0; i < k; i++) {
        const a = H[i][k], c = H[i + 1][k];
        H[i][k] = cs[i] * a + sn[i] * c;
        H[i + 1][k] = -sn[i] * a + cs[i] * c;
      }
      const t = Math.hypot(H[k][k], hn), c = t > 0 ? H[k][k] / t : 1, sv = t > 0 ? hn / t : 0;
      cs.push(c);
      sn.push(sv);
      H[k][k] = t;
      g[k + 1] = -sv * g[k];
      g[k] *= c;
      if (!(Math.abs(g[k + 1]) > tol * bn) || !(hn > 0)) {
        k++;
        converged = true;
        break;
      }
    }
    const y = new Float64Array(k);
    for (let i = k - 1; i >= 0; i--) {
      let v = g[i];
      for (let j = i + 1; j < k; j++) v -= H[i][j] * y[j];
      y[i] = v / H[i][i];
    }
    for (let i = 0; i < k; i++) for (let j = 0; j < N; j++) x[j] += y[i] * V[i][j];
    if (converged) return true;
  }
  return false;
}

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
// already preconditioned; x holds the starting guess and receives the answer.
function gmres(op, b, x, { m = 12, restarts = 4, tol = 1e-8 } = {}) {
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
  for (let cycle = 0; cycle < restarts; cycle++) {
    const ax = op(x), r = [Float64Array.from(b[0], (v, j) => v - ax[0][j]), Float64Array.from(b[1], (v, j) => v - ax[1][j])];
    const beta = norm(r);
    if (!(beta > tol * bn)) return;
    const V = [[r[0].map((v) => v / beta), r[1].map((v) => v / beta)]];
    const H = Array.from({ length: m + 1 }, () => Array.from({ length: m }, () => [0, 0]));
    const cs = [], sn = [], g = Array.from({ length: m + 1 }, () => [0, 0]);
    g[0] = [beta, 0];
    let k = 0, converged = false;
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
        converged = true;
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
    if (converged) return; // (by its own residual estimate)
  }
}

// Discretisation and nonlinear solver.
//
// Unknowns, per solver block (M = 1 + nSpecies slots):
//   grid node: [φ̂, η_1 … η_n], with φ̂ = Fφ/RT and η_i = μ̄_i/RT
//   flux node: [D, N_1 … N_n, r_1 … r_K], the displacement and particle fluxes through an
//              interface, and the rate of each reaction at it
// Every region boundary is a doubled grid node (one per side) with a zero-volume flux node
// between them in the linear system, so the Jacobian stays block-tridiagonal.
//
// Not every slot is an unknown everywhere: an absent species, φ where it's undefined, a blocked
// interface flux. The state u keeps every slot (one that isn't an unknown just keeps its value),
// but the residual, the update and the Jacobian hold only the unknowns: block b has as many as
// it has active slots, and the linear system's block sizes vary from node to node. `loc` maps a
// slot to its row within its block (or −1), and `rix` to its index in the compact vectors (or
// the sink, one spare entry at the end that absorbs writes for slots that aren't unknowns).
//
// Balance rows (node g, box volume v per unit area):
//   φ:   D_out − D_in − v·(F Σ z_i c_i + ρ_fixed) = 0
//   i:   v·(c_i − c_i,old)/dt + N_out − N_in = 0
// Fluxes along segments are Scharfetter–Gummel. Each flux is computed once per segment and
// added with opposite signs to both neighbours, so sums over boxes telescope exactly.
//
// Concentrations come from each material's statistics, c(ζ) with ζ_i = η_i − μ°_i/RT − z_i φ̂,
// and K = ∂c/∂ζ. Nodes of ideal materials (c = c_ref e^ζ) take a fast path throughout.

import { BlockTridiagonal, ComplexBlockTridiagonal } from './blockTridiagonal.js';
import { bernoulli, bernoulliDerivative } from './bernoulli.js';
import { EPS0, FARADAY } from './constants.js';
import { nextBreakpoint, sourceAt } from './device.js';
import { DeviceError } from './errors.js';
import { powi, powr } from './pow.js';

// A profile's value at x (c0, kf): piecewise linear, constant beyond its ends.
function profileAt({ x, values: c }, at) {
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
    // one per face reaction (its rate), as many as the busiest face has.
    const nRx = model.interfaces.reduce((m, itf) => Math.max(m, itf.reactions.length), 0);
    // (and at an electrode port's nodes, one per surface species: its coverage's η)
    const nSurf = model.ports.reduce((m, port) => Math.max(m, port.surface.length), 0);
    const M = n + 1 + Math.max(nRx, nSurf);
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
    // pinned contact, or a capacitive face to a conductor. (Else φ there is a free constant.)
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
      if (r0 > 0 && model.interfaces[r0 - 1].phi.type !== 'neutral' && materials[regions[r0 - 1].material].conductor) anchor = true;
      if (r1 < regions.length - 1 && model.interfaces[r1].phi.type !== 'neutral' && materials[regions[r1 + 1].material].conductor) anchor = true;
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
    // The floating ones are extra unknowns, solved with the grid's by bordering (_solveBordered).
    this.floating = this.terms.flatMap((t, k) => (t.drive.kind === 'I' || t.drive.R > 0 ? [k] : []));
    // Per terminal: B = ∂res/∂V (compact column), C = ∂I/∂x (compact row), ∂I/∂V, and for a
    // floating one its circuit residual (I − I_set, or I − (V_src − V)/R).
    this.termB = this.terms.map(() => new Float64Array(this.sys.size + 1));
    this.termC = this.terms.map(() => new Float64Array(this.sys.size + 1));
    this.termDI = new Float64Array(this.terms.length);
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
    }
    model.interfaces.forEach((itf, f) => {
      const bf = this.blockOfFace[f];
      active[bf * M] = itf.phi.type === 'neutral' ? 0 : 1;
      for (let i = 0; i < n; i++) active[bf * M + 1 + i] = itf.links[i].type === 'blocked' ? 0 : 1;
      for (let k = 0; k < itf.reactions.length; k++) active[bf * M + 1 + n + k] = 1;
    });
    const sizes = new Int32Array(nB), loc = (this.loc = new Int32Array(nB * M).fill(-1));
    for (let b = 0; b < nB; b++) for (let r = 0; r < M; r++) if (active[b * M + r]) loc[b * M + r] = sizes[b]++;
    const sys = (this.sys = new BlockTridiagonal(nB, sizes));
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
    this.rowScale = new Float64Array(N + 1).fill(1); // each row's scaling in the factorised system
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
    this.moieties = this._conservedMoieties();
    this.stretches.forEach((st, k) => {
      st.connected = st.contactFed || this.moieties.every((w) => w[k] === 0);
    });
    for (const st of this.stretches) {
      st.spectator = !st.connected && !st.reactive;
      // Mobile throughout: its steady state is then a single level fixed by its amount.
      st.mobile = true;
      for (let q = st.regions[0]; q <= st.regions[1]; q++) if (!(materials[regions[q].material].D[st.species] > 0)) st.mobile = false;
    }
    // Stretches whose steady state is known outright: a species that no reaction or port touches,
    // reached by one contact only (the other end blocked), carries no flux at steady state, so its
    // level is flat at that contact's. Steady solves pin it there (_flatRows) rather than find it
    // through its own conduction, which can be all but nothing: a MOS capacitor's inversion
    // electrons reach the back contact only through a bulk with ~1e3 of them per cm³. (Not
    // through flow or mixing, nor a concentrated material's cross-diffusion, where zero flux
    // isn't a flat level.)
    this.flatStretches = this.stretches.filter((st) => {
      if (st.reactive || st.ports.length > 0 || !st.mobile || st.leftOpen === st.rightOpen) return false;
      for (let q = st.regions[0]; q <= st.regions[1]; q++) {
        const reg = regions[q], mat = materials[reg.material];
        if (mat.conductor || !mat.ideal || reg.velocity !== 0 || reg.mixing > 0) return false;
      }
      return true;
    });
    this.flattening = false;
    // Conserved amounts solved directly, each in place of one (redundant) balance row: a
    // spectator's own amount, and each conserved combination of reacting stretches (the total
    // iron of Fe³⁺, Fe²⁺ and FeCl²⁺, say), whose weighted balances sum to zero at steady state.
    // Immobile stretches conserve node by node instead, and a floating conductor holds charge on
    // its faces; both are left to giant time steps.
    this.constraints = [];
    const conductor = (st) => materials[regions[st.regions[0]].material].conductor;
    const add = (parts, rowStretch, surface = []) => {
      const st = this.stretches[rowStretch];
      let nNodes = parts.reduce((m, p) => m + this.stretches[p.stretch].nodes[1] - this.stretches[p.stretch].nodes[0] + 1, 0);
      for (const p of surface) nNodes += model.ports[this.surfCols[p.col][0]].nodes.length * model.ports[this.surfCols[p.col][0]].surface.length;
      for (const p of parts) this.stretches[p.stretch].conserved = true;
      this.constraints.push({
        parts,
        surface,
        row: this.blockOfNode[st.nodes[0]] * this.M + 1 + st.species,
        idx: new Int32Array(nNodes * this.M),
        w: new Float64Array(nNodes * this.M),
        len: 0,
        res: 0,
        q: new Float64Array(this.sys.size + 1),
      });
    };
    this.stretches.forEach((st, k) => {
      if (st.spectator && st.mobile) add([{ stretch: k, w: 1 }], k);
    });
    const S = this.stretches.length;
    for (const w of this.moieties) {
      const parts = [...w.keys()].filter((k) => k < S && w[k] !== 0).map((k) => ({ stretch: k, w: w[k] }));
      const surface = [...w.keys()].filter((k) => k >= S && w[k] !== 0).map((k) => ({ col: k - S, w: w[k] }));
      if (parts.length === 0) continue; // (a surface's own: left to giant time steps)
      if (parts.length === 1 && surface.length === 0 && this.stretches[parts[0].stretch].spectator) continue; // (above)
      if (parts.some((p) => !this.stretches[p.stretch].mobile || conductor(this.stretches[p.stretch]))) continue;
      // The row replaced: the basis vector's own stretch (its weight is 1, and no other vector has one).
      add(parts, parts.find((p) => w[p.stretch] === 1)?.stretch ?? parts[0].stretch, surface);
    }
    // A combination of immobile stretches (trap states: X⁰ + X⁻ under e⁻ + X⁰ = X⁻) is conserved
    // node by node, since nothing carries it anywhere: at each node its weighted sum stays what
    // it was, a row local to that node's block, in place of the balance row of one of them.
    this.localConstraints = [];
    const immobile = (st) => {
      for (let q = st.regions[0]; q <= st.regions[1]; q++) if (materials[regions[q].material].D[st.species] > 0) return false;
      return true;
    };
    for (const w of this.moieties) {
      const parts = [...w.keys()].filter((k) => w[k] !== 0).map((k) => ({ stretch: k, w: w[k] }));
      const sts = parts.map((p) => this.stretches[p.stretch]);
      if (sts.some((st) => st.conserved || !immobile(st) || conductor(st))) continue;
      if (sts.some((st) => st.nodes[0] !== sts[0].nodes[0] || st.nodes[1] !== sts[0].nodes[1])) continue; // (side by side only)
      for (const st of sts) st.conserved = true;
      const row = parts.find((p) => p.w === 1) ?? parts[0];
      this.localConstraints.push({
        parts: parts.map((p) => ({ i: this.stretches[p.stretch].species, w: p.w })),
        row: 1 + this.stretches[row.stretch].species,
        nodes: sts[0].nodes,
        reference: new Float64Array(sts[0].nodes[1] - sts[0].nodes[0] + 1),
      });
    }
    this.constrained = false;
  }

  // Each local constraint's weighted sum at each node, as the state holds it now: what a steady
  // solve keeps.
  _captureLocal() {
    const { n, c } = this;
    for (const lc of this.localConstraints) {
      for (let g = lc.nodes[0]; g <= lc.nodes[1]; g++) {
        let t = 0;
        for (const { i, w } of lc.parts) t += w * c[g * n + i];
        lc.reference[g - lc.nodes[0]] = t;
      }
    }
  }

  // The local constraints' rows (steady solves only): Σ w c_i − reference at each node.
  _applyLocalConstraints() {
    const { n, M, c, res, dA: d } = this;
    for (const lc of this.localConstraints) {
      for (let g = lc.nodes[0]; g <= lc.nodes[1]; g++) {
        const b = this.blockOfNode[g];
        if (this.loc[b * M + lc.row] < 0) continue;
        this._replaceRow(b, lc.row);
        let t = 0;
        for (const { i, w } of lc.parts) {
          t += w * c[g * n + i];
          this._dc(g, i, d);
          for (let s = 0; s < M; s++) if (d[s] !== 0) this._j(b, lc.row, b, s, w * d[s]);
        }
        res[this.rix[b * M + lc.row]] = t - lc.reference[g - lc.nodes[0]];
      }
    }
  }

  // A basis of the conserved combinations of stretch amounts: the null space of the
  // stoichiometry (a row per reaction, over the stretches it touches) together with a unit row
  // for each stretch fed from outside. Exact integer data, so plain elimination will do.
  _conservedMoieties() {
    const { model, n } = this, S = this.stretches.length, C = S + this.surfCols.length;
    const rows = [];
    model.interfaces.forEach((itf, f) => {
      for (const rx of itf.reactions) {
        const row = new Float64Array(C);
        for (const p of rx.part) row[this.stretchOf[(f + p.side) * n + p.i]] += p.nu;
        rows.push(row);
      }
    });
    // An electrode surface's species are columns too, after the stretches (S + j for surface
    // column j), so a combination that passes through them (A⁺ + e⁻ = S, S + e⁻ = B⁻) is conserved
    // with what the surface holds.
    model.ports.forEach((port, k) => {
      for (const rx of port.reactions) {
        const row = new Float64Array(C); // (the electrode's carrier is outside the device)
        for (const p of rx.part) {
          if (p.side === 0) row[this.stretchOf[port.region * n + p.i]] += p.nu;
          else if (p.side === 2) row[S + this.surfColOf[k] + p.s] += p.nu;
        }
        rows.push(row);
      }
    });
    model.regions.forEach((reg, q) => {
      for (const rx of model.reactions) {
        if (!(rx.kf[reg.material] > 0)) continue;
        const row = new Float64Array(C);
        for (const p of rx.reactants) row[this.stretchOf[q * n + p.i]] -= p.nu;
        for (const p of rx.products) row[this.stretchOf[q * n + p.i]] += p.nu;
        rows.push(row);
      }
    });
    this.stretches.forEach((st, k) => {
      if (!st.contactFed) return;
      const row = new Float64Array(C);
      row[k] = 1;
      rows.push(row);
    });
    // Reduced row echelon form; the free columns give the null space.
    const pivots = [];
    let r = 0;
    for (let col = 0; col < C && r < rows.length; col++) {
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
      basis.push(w);
    }
    return basis;
  }

  // Every stretch either reaches a contact or is a mobile spectator: the steady equations can
  // be solved directly (dt = ∞), with the spectators' amounts as constraints.
  _directSteady() {
    return this.stretches.every((st) => st.connected || st.conserved);
  }

  // Conservation rows for the spectators (steady solves only): Σ v c_i = amount over the
  // stretch replaces the balance row of its first node, which in steady state is the negative
  // sum of the others. That row is dense, so it's kept aside: the factorised matrix gets a pin
  // (identity row) there instead, and _solveBordered restores the constraint.
  _applyConstraints() {
    const { n, M, res, c, dA: d } = this;
    const vol = this.model.grid.vol;
    for (const cs of this.constraints) {
      let amount = 0, reference = 0, len = 0;
      for (const { stretch, w } of cs.parts) {
        const st = this.stretches[stretch], i = st.species;
        reference += w * this.referenceAmounts[stretch];
        for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
          amount += w * vol[g] * c[g * n + i];
          this._dc(g, i, d);
          const b = this.blockOfNode[g];
          for (let r = 0; r < M; r++) {
            if (d[r] === 0 || this.loc[b * M + r] < 0) continue;
            cs.idx[len] = this.rix[b * M + r];
            cs.w[len++] = w * vol[g] * d[r];
          }
        }
      }
      // What the electrodes' surfaces hold of it: Γ Σ v·a·θ, with ∂θ_s/∂η_x = θ_s (δ_sx − θ_x).
      for (const { col, w } of cs.surface) {
        const [k, q] = this.surfCols[col], port = this.model.ports[k], G = port.surface[0].capacity, ns = port.surface.length;
        reference += w * this.surfaceRef[col];
        port.nodes.forEach((g, j) => {
          const f = w * G * vol[g] * this.portArea[k][j], b = this.blockOfNode[g], tq = this.th[g * this.nSurf + q];
          amount += f * tq;
          for (let x = 0; x < ns; x++) {
            cs.idx[len] = this.rix[b * M + 1 + n + x];
            cs.w[len++] = f * tq * ((q === x ? 1 : 0) - this.th[g * this.nSurf + x]);
          }
        });
      }
      cs.len = len;
      cs.res = amount - reference;
      const b0 = Math.floor(cs.row / M), r0 = cs.row % M;
      this._replaceRow(b0, r0);
      this._j(b0, r0, b0, r0, 1);
      res[this.rix[cs.row]] = 0;
    }
  }

  // Solve the full Newton system, J δ = rhs, with the extra unknowns and rows that don't fit the
  // block-tridiagonal matrix T: the floating terminals' voltages (with their circuit rows), and
  // the spectators' conservation rows (pinned in T, see _applyConstraints). By low-rank
  // updates of T (Woodbury): δ = y + Σ_q Q_q μ_q − Σ_k X_k δV_k with y = T⁻¹ rhs, Q_q = T⁻¹ e_q
  // (the response to a unit pin) and X_k = T⁻¹ B_k (to a unit change of V_k), then a small dense
  // system for the μ (pins) and δV (terminals):
  //   pin q:      Σ (W_q·Q_q′) μ_q′ − Σ (W_q·X_k) δV_k = res_q − W_q·y
  //   terminal k: Σ (C_k·Q_q) μ_q + Σ (δ_kk′ ∂I_k/∂V_k − C_k·X_k′) δV_k′ = res_k − C_k·y
  _solveBordered(rhs, delta, deltaV, pins) {
    const N = this.sys.size, fl = this.floating, P = pins.length, K = fl.length, S = P + K;
    this._solveLinear(rhs, delta);
    if (S === 0) return;
    const e = this.dWork ?? (this.dWork = new Float64Array(N + 1));
    const cols = this.borderCols ?? (this.borderCols = []);
    for (let a = 0; a < S; a++) {
      if (!cols[a]) cols[a] = new Float64Array(N + 1);
      if (a < P) {
        const row = this.rix[pins[a].row];
        e[row] = 1;
        this._solveLinear(e, cols[a]);
        e[row] = 0;
      } else this._solveLinear(this.termB[fl[a - P]], cols[a]);
    }
    // Row a of the small system, as a dot product with a compact vector (sparse for pins).
    const dot = (a, v) => {
      let t = 0;
      if (a < P) {
        const p = pins[a];
        for (let j = 0; j < p.len; j++) t += p.w[j] * v[p.idx[j]];
      } else {
        const C = this.termC[fl[a - P]];
        for (let j = 0; j < N; j++) t += C[j] * v[j];
      }
      return t;
    };
    const A = Array.from({ length: S }, () => new Float64Array(S + 1));
    for (let a = 0; a < S; a++) {
      for (let b = 0; b < S; b++) A[a][b] = (b < P ? 1 : -1) * dot(a, cols[b]);
      if (a >= P) A[a][a] += this.termDI[fl[a - P]];
      A[a][S] = (a < P ? pins[a].res : this.termRes[fl[a - P]]) - dot(a, delta);
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
    const before = this.floating;
    this.floating = this.terms.flatMap((t, k) => (t.drive.kind === 'I' || t.drive.R > 0 ? [k] : []));
    if (this.floating.length !== before.length) this.deltaV = null;
    this._refreshSources();
    this.history = [];
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
    if (link.type !== 'equilibrium') return NaN;
    const z = this.z[i];
    return z === 0 ? link.mu / this.model.RT : (z * (this.termV[side === 'left' ? 0 : 1] + link.offset)) / this.VT;
  }

  /**
   * Cold start. Each species takes its region's c0 where one is given (for a spectator, it fixes
   * the conserved amount), else its contact's level. φ in each region is then chosen for local neutrality (node by node where a
   * profile varies), or continued across the interface dipole if nothing there responds.
   */
  initFromComposition() {
    const { model, n, M, u, z } = this;
    const { regions, species, interfaces, materials } = model;
    const grid = model.grid;
    u.fill(0);
    this.uLo.fill(0);
    const eta = new Float64Array(n), cFix = new Float64Array(n), mode = new Int8Array(n); // 1 level, 2 amount
    // φ̂ runs from the left: from a pinned (or gated) left contact's φ at its zero charge, else 0.
    const leftPhi = model.contacts.left.phi;
    let phiHat = leftPhi.type === 'pinned' || leftPhi.type === 'capacitive' ? (this.termV[0] - leftPhi.zeroCharge) / this.VT : 0;
    for (let r = 0; r < regions.length; r++) {
      if (r > 0) phiHat += interfaces[r - 1].dipole / this.VT;
      const reg = regions[r], mat = materials[reg.material];
      mode.fill(0);
      const profiled = reg.c0Profile.some((p) => p !== null);
      for (let i = 0; i < n; i++) {
        if (!mat.present[i]) continue;
        const st = this.stretches[this.stretchOf[r * n + i]];
        if (reg.c0Profile[i] || reg.c0[i] > 0) {
          // A given c0 is the starting state, even where a contact feeds the species.
          mode[i] = 2;
          cFix[i] = reg.c0Profile[i] ? profileAt(reg.c0Profile[i], grid.x[grid.regionStart[r]]) : reg.c0[i];
        } else if (!st.contactFed && mat.conductor) {
          // A conductor away from the contacts: start uncharged, with its carrier's level in
          // equilibrium with the first reaction on its left face that takes it, else at the
          // running φ.
          eta[i] = z[i] * phiHat;
          const rx = r > 0 ? interfaces[r - 1].reactions.find((x) => x.part.some((p) => p.side === 1 && p.i === i)) : undefined;
          if (rx) {
            const bn = this.blockOfNode[grid.regionEnd[r - 1]];
            let a = rx.fixedA, nu = 0; // a = fixedA − Σ ν η = 0, solved for the carrier's η
            for (const p of rx.part) {
              if (p.side === 1) nu += p.nu;
              else a -= p.nu * u[bn * M + 1 + p.i];
            }
            eta[i] = a / nu;
          }
          mode[i] = 1;
        } else if (!st.contactFed) {
          if (!(reg.c0[i] > 0)) {
            throw new SolverError(
              `regions[${r}].c0.${species[i].name}: a species that doesn't reach a contact needs its initial concentration` +
                (st.connected ? ' (reactions make and consume it, so this is only where the solve starts, not an amount it keeps: a tiny value will do)' : ', which fixes the amount it conserves'),
            );
          }
          mode[i] = 2;
          cFix[i] = reg.c0[i];
        } else {
          const left = st.regions[0] === 0 ? this.contactEta('left', i) : NaN;
          eta[i] = Number.isFinite(left) ? left : this.contactEta('right', i);
          if (!Number.isFinite(eta[i]) && st.linked.length > 0) eta[i] = this.portEta(model.ports[st.linked[0]], i);
          mode[i] = 1;
          if (!Number.isFinite(eta[i])) {
            // Fed only through reactions (no level held at a contact): start from c0.
            if (!(reg.c0[i] > 0)) {
              throw new SolverError(
                `regions[${r}].c0.${species[i].name}: a species fed only through reactions needs its initial concentration`,
              );
            }
            mode[i] = 2;
            cFix[i] = reg.c0[i];
          }
        }
      }
      // Net charge (mol/m³) at trial φ̂; decreasing in φ̂ wherever a level-fixed ion responds.
      const zeta = new Float64Array(n), cc = new Float64Array(n);
      const charge = (ph) => {
        if (!mat.ideal) {
          this._materialAt(mat, reg.background, eta, mode, cFix, ph, zeta, cc);
          let q = reg.fixedCharge / FARADAY;
          for (let i = 0; i < n; i++) if (mat.present[i]) q += z[i] * cc[i];
          return q;
        }
        let q = reg.fixedCharge / FARADAY;
        for (let i = 0; i < n; i++) {
          if (mode[i] === 2) q += z[i] * cFix[i];
          else if (mode[i] === 1 && z[i] !== 0) {
            const ex = Math.min(700, Math.max(-700, eta[i] - mat.mu0[i] / model.RT - z[i] * ph));
            q += z[i] * mat.cRef[i] * Math.exp(ex);
          }
        }
        return q;
      };
      const responds = mode.some((m, i) => m === 1 && z[i] !== 0);
      const mustBalance = !responds && !mat.conductor && !mat.phiFree && mat.epsr === 0;
      const checkNeutral = (where) => {
        // Strictly neutral, with every charged species' amount given: its c0 must be neutral.
        let q = reg.fixedCharge / FARADAY, scale = Math.abs(q);
        for (let i = 0; i < n; i++) {
          if (mode[i] !== 2 || z[i] === 0) continue;
          q += z[i] * cFix[i];
          scale += Math.abs(z[i] * cFix[i]);
        }
        if (Math.abs(q) > 1e-9 * scale) {
          throw new SolverError(
            `regions[${r}] (${reg.name}): its initial composition carries a net charge of ${q.toPrecision(3)} mol/m³${where} (Σ z·c0, with any fixed charge), ` +
              'but ε = 0 makes it strictly neutral; adjust c0 so that it balances',
          );
        }
      };
      if (mustBalance && !profiled) checkNeutral('');
      if (!mat.phiFree && !responds && mat.ideal && r > 0 && materials[regions[r - 1].material].conductor) {
        // Nothing here fixes φ, but an electrode on the left does: start with its first reaction
        // that takes a species from this side at equilibrium (the electrode at its open-circuit
        // level), rather than φ carried over from the metal, which can be volts away.
        const bn = this.blockOfNode[grid.regionEnd[r - 1]];
        for (const rx of interfaces[r - 1].reactions) {
          let a = rx.fixedA, s = 0;
          for (const p of rx.part) {
            if (p.side === 0) a -= p.nu * u[bn * M + 1 + p.i];
            else if (mode[p.i] === 2) {
              a -= p.nu * (Math.log(cFix[p.i] / mat.cRef[p.i]) + mat.mu0[p.i] / model.RT);
              s += p.nu * z[p.i];
            } else s = NaN;
          }
          if (s !== 0 && Number.isFinite(s)) {
            phiHat = a / s; // a − s·φ̂ = 0
            break;
          }
        }
      }
      const electrode = model.ports.find((port) => port.region === r && port.reactions.length > 0);
      if (!mat.phiFree && !responds && mat.ideal && electrode) {
        // Nor here, but an electrode spread through it does: start with the port's first reaction
        // that can be balanced at equilibrium with the port's level (the electrode at open circuit).
        for (const rx of electrode.reactions) {
          let a = rx.fixedA, s = 0;
          for (const p of rx.part) {
            if (p.side === 1) a -= p.nu * this.portEta(electrode, p.i);
            else if (mode[p.i] === 2) {
              a -= p.nu * (Math.log(cFix[p.i] / mat.cRef[p.i]) + mat.mu0[p.i] / model.RT);
              s += p.nu * z[p.i];
            } else s = NaN;
          }
          if (s !== 0 && Number.isFinite(s)) {
            phiHat = a / s;
            break;
          }
        }
      }
      const neutralPhi = (ph) => {
        let lo = ph - 1, hi = ph + 1;
        while (charge(lo) < 0 && lo > -1e4) lo -= 2 * (hi - lo);
        while (charge(hi) > 0 && hi < 1e4) hi += 2 * (hi - lo);
        // Out of reach (every responding carrier has one sign and nothing balances it, as in an
        // undoped layer that holds only holes): keep the running φ̂ rather than run it to ±∞.
        if (charge(lo) < 0 || charge(hi) > 0) return ph;
        for (let it = 0; it < 200 && hi - lo > 1e-12; it++) {
          const m = 0.5 * (lo + hi);
          if (charge(m) > 0) lo = m;
          else hi = m;
        }
        return 0.5 * (lo + hi);
      };
      if (!mat.phiFree && responds) phiHat = neutralPhi(phiHat);
      if (!mat.ideal) this._materialAt(mat, reg.background, eta, mode, cFix, phiHat, zeta, cc);
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        if (profiled) {
          for (let i = 0; i < n; i++) if (reg.c0Profile[i]) cFix[i] = profileAt(reg.c0Profile[i], grid.x[g]);
          if (mustBalance) checkNeutral(` at x = ${grid.x[g].toPrecision(4)} m`);
          if (!mat.phiFree && responds) phiHat = neutralPhi(phiHat);
          if (!mat.ideal) this._materialAt(mat, reg.background, eta, mode, cFix, phiHat, zeta, cc);
        }
        const b = this.blockOfNode[g];
        u[b * M] = mat.conductor ? 0 : phiHat; // a metal's slot 0 is its segment flux
        for (let i = 0; i < n; i++) {
          if (mode[i] === 1) u[b * M + 1 + i] = eta[i];
          else if (mode[i] === 2) {
            const zt = mat.ideal ? Math.log(cFix[i] / mat.cRef[i]) : zeta[i];
            u[b * M + 1 + i] = zt + mat.mu0[i] / model.RT + z[i] * phiHat;
          }
        }
      }
    }
    // A floating contact facing a held one starts where the start's own composition puts it: the
    // held voltage plus the difference between the two terminal species' levels beside each
    // contact (c0 sets them, and φ was chosen for neutrality, not to match either contact), so
    // that a transient's first reading means something.
    const level = (side) => {
      const ct = model.contacts[side], i = ct.terminal;
      if (i === null || ct.species[i].type !== 'equilibrium' || z[i] === 0) return NaN;
      const b = this.blockOfNode[side === 'left' ? 0 : grid.nNodes - 1];
      return (this.VT * u[b * M + 1 + i]) / z[i] - ct.species[i].offset;
    };
    for (const k of this.floating) {
      const t = this.terms[k], held = k === 0 ? 1 : 0;
      if (t.kind === 'port' || held >= this.terms.length || this.floating.includes(held)) continue;
      const shift = level(t.side) - level(this.terms[held].side);
      if (Number.isFinite(shift)) this.termV[k] = this.termV[held] + shift;
    }
    for (let g = 0; g < this.nNodes; g++) if (this.surfPort[g] >= 0) this._surfaceStart(g);
    this.computeConcentrations();
    // A floating capacitance starts uncharged (σ = 0 at the window's mean φ), as the starting
    // composition, neutral without it, assumes; with reactions too, its double layer then charges
    // toward their mixed potential.
    for (const k of this.floating) {
      const t = this.terms[k], port = t.kind === 'port' ? model.ports[t.index] : null;
      if (!port?.capacitance) continue;
      let w = 0, sum = 0;
      port.nodes.forEach((g, j) => {
        if (this.phiUndefined[g]) return;
        const a = this.model.grid.vol[g] * this.portArea[t.index][j];
        w += a;
        sum += a * this.VT * u[this.blockOfNode[g] * M];
      });
      if (w > 0) this.termV[k] = port.capacitance.zeroCharge + sum / w;
    }
    // A floating electrode spread through a port (without a double layer) starts where its
    // reactions pass the current it's set (none, behind a resistance): at its mixed potential in
    // the start's composition, not level with a held terminal, which can be volts away and pass
    // an absurd current.
    for (const k of this.floating) {
      const t = this.terms[k], port = t.kind === 'port' ? model.ports[t.index] : null;
      if (!port || port.reactions.length === 0 || port.capacitance) continue;
      const want = t.drive.kind === 'I' ? sourceAt(t.drive.src, this.time) : 0;
      const current = (V) => {
        this.termV[k] = V;
        let I = 0;
        port.reactions.forEach((rx) => {
          let q = 0;
          for (const p of rx.part) if (p.side === 0) q += p.nu * z[p.i];
          port.nodes.forEach((g, w) => (I += q * FARADAY * this.model.grid.vol[g] * this.portArea[t.index][w] * this._portRate(port, rx, g)));
        });
        return I;
      };
      // Current into the device rises with V (oxidation); bracket, then bisect.
      let lo = this.termV[k] - 1, hi = this.termV[k] + 1;
      while (current(lo) > want && lo > -1e3) lo -= 2 * (hi - lo);
      while (current(hi) < want && hi < 1e3) hi += 2 * (hi - lo);
      for (let it = 0; it < 200 && hi - lo > 1e-12; it++) {
        const m = 0.5 * (lo + hi);
        if (current(m) < want) lo = m;
        else hi = m;
      }
      this.termV[k] = 0.5 * (lo + hi);
    }
  }

  // A port reaction's rate per area at node g (mol/(m²·s), forward), at the present state.
  _portRate(port, rx, g) {
    const { n, M, u, uLo, c, nSurf, th } = this, b = this.blockOfNode[g], RT = this.model.RT;
    let pref = rx.k0, a = rx.fixedA;
    for (const p of rx.part) {
      if (p.side === 1) {
        a -= p.nu * this.portEta(port, p.i);
        continue;
      }
      const o = b * M + (p.side === 2 ? 1 + n + p.s : 1 + p.i), e = p.nu < 0 ? -p.nu * (1 - rx.alpha) : p.nu * rx.alpha;
      a -= p.nu * (u[o] + uLo[o]);
      pref *= p.side === 2 ? Math.exp(e * (u[o] + uLo[o] - port.surface[p.s].mu0 / RT)) : powr(c[g * n + p.i] / this.cRef[g * n + p.i], e);
    }
    if (rx.bare) pref *= this.th0[g];
    return pref * bvFactor(a, rx.alpha).g;
  }

  // A non-ideal material's composition at φ̂ = ph: level-fixed species (mode 1) at their η,
  // amount-fixed ones (mode 2) at cFix. Fills ζ and c for every present species.
  _materialAt(mat, background, eta, mode, cFix, ph, zeta, cc) {
    const { n, z } = this;
    const RT = this.model.RT;
    for (let i = 0; i < n; i++) {
      if (!mat.present[i]) continue;
      zeta[i] = mode[i] === 2 ? Math.log(cFix[i] / mat.cRef[i]) : eta[i] - mat.mu0[i] / RT - z[i] * ph;
      if (mat.modelOf[i] < 0) cc[i] = mode[i] === 2 ? cFix[i] : mat.cRef[i] * Math.exp(Math.min(700, zeta[i]));
    }
    for (const md of mat.models) {
      const w = this._work(md), idx = md.idx, k = idx.length;
      for (let a = 0; a < k; a++) {
        w.z[a] = zeta[idx[a]];
        w.fixed[a] = mode[idx[a]] === 2 ? 1 : 0;
        w.t[a] = cFix[idx[a]];
      }
      md.invert(w.z, w.fixed, w.t, background);
      md.evaluate(w.z, w.c, w.K, background);
      for (let a = 0; a < k; a++) {
        zeta[idx[a]] = w.z[a];
        cc[idx[a]] = w.c[a];
      }
    }
  }

  // Each electrode surface's coverages from its η: Langmuir on shared sites, θ_s = e^{ζ_s}/(1 + Σ e^{ζ}),
  // ζ = η − μ°/RT (so θ_s/θ₀ = e^{ζ_s}), at every node of a surfaced window.
  _coverages() {
    const { n, M, u, uLo, nSurf, th } = this, RT = this.model.RT;
    for (let g = 0; g < this.nNodes; g++) {
      const k = this.surfPort[g];
      if (k < 0) continue;
      const surf = this.model.ports[k].surface, b = this.blockOfNode[g];
      let m = 0;
      for (let q = 0; q < surf.length; q++) m = Math.max(m, u[b * M + 1 + n + q] + uLo[b * M + 1 + n + q] - surf[q].mu0 / RT);
      let sum = Math.exp(-m);
      for (let q = 0; q < surf.length; q++) sum += (th[g * nSurf + q] = Math.exp(u[b * M + 1 + n + q] + uLo[b * M + 1 + n + q] - surf[q].mu0 / RT - m));
      for (let q = 0; q < surf.length; q++) th[g * nSurf + q] /= sum;
      this.th0[g] = Math.exp(-m) / sum;
    }
  }

  computeConcentrations() {
    const { n, M, u, uLo, c, z } = this;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      const phiHat = u[b * M];
      if (!this.nodeIdeal[g]) {
        this._nodeStatistics(g, b, phiHat);
        continue;
      }
      const im = this.nodeConductor[g];
      if (im >= 0) {
        // Excess carriers: zero in the bulk, the surface sheet at a charged face.
        for (let i = 0; i < n; i++) c[g * n + i] = 0;
        const f = this.sheetFace[g];
        if (f >= 0) c[g * n + im] = (this.sheetSign[g] * u[this.blockOfFace[f] * M] * this.model.grid.area[g]) / (z[im] * FARADAY * this.model.grid.vol[g]);
        continue;
      }
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        c[k] = this.present[k] ? this.cRef[k] * Math.exp(u[b * M + 1 + i] + uLo[b * M + 1 + i] - this.mu0hat[k] - z[i] * phiHat) : 0;
      }
    }
    this._coverages();
  }

  _work(md) {
    let w = this.scratch.get(md);
    if (!w) {
      const k = md.idx.length;
      w = { z: new Float64Array(k), c: new Float64Array(k), K: new Float64Array(k * k), fixed: new Uint8Array(k), t: new Float64Array(k) };
      this.scratch.set(md, w);
    }
    return w;
  }

  // c, K and the excess at a node of a non-ideal material.
  _nodeStatistics(g, b, phiHat) {
    const { n, M, u, uLo, c, z, K, ex, zeta } = this;
    const mat = this.model.materials[this.nodeMaterial[g]];
    const Kg = g * n * n;
    K.fill(0, Kg, Kg + n * n);
    for (let i = 0; i < n; i++) {
      const k = g * n + i;
      ex[k] = 0;
      if (!this.present[k]) {
        c[k] = 0;
        continue;
      }
      zeta[i] = u[b * M + 1 + i] + uLo[b * M + 1 + i] - this.mu0hat[k] - z[i] * phiHat;
      if (mat.modelOf[i] < 0) {
        c[k] = this.cRef[k] * Math.exp(zeta[i]);
        K[Kg + i * n + i] = c[k];
      }
    }
    for (const md of mat.models) {
      const w = this._work(md), idx = md.idx, k = idx.length;
      for (let a = 0; a < k; a++) w.z[a] = zeta[idx[a]];
      md.evaluate(w.z, w.c, w.K, this.background[g]);
      for (let a = 0; a < k; a++) {
        const i = idx[a], ci = w.c[a];
        c[g * n + i] = ci;
        ex[g * n + i] = ci > 0 ? zeta[i] - Math.log(ci / this.cRef[g * n + i]) : 0;
        for (let q = 0; q < k; q++) K[Kg + i * n + idx[q]] = w.K[a * k + q];
      }
    }
  }

  // Σ_j K_ij z_j at node g (= z_i c_i for ideal statistics).
  _Kz(g, i) {
    const { n, z } = this;
    if (this.nodeIdeal[g]) return z[i] * this.c[g * n + i];
    let s = 0;
    const o = g * n * n + i * n;
    for (let j = 0; j < n; j++) s += this.K[o + j] * z[j];
    return s;
  }

  /** zᵀKz at node g (per RT): the charge response that sets the screening length. */
  screening(g) {
    let s = 0;
    for (let i = 0; i < this.n; i++) s += this.z[i] * this._Kz(g, i);
    return s;
  }

  // Add w·∂(ln c_i)/∂(slot) at node g into d (slot 0 = φ̂, 1 + j = η_j).
  _dlnc(g, i, w, d) {
    const { n, z } = this;
    if (this.nodeIdeal[g]) {
      d[1 + i] += w;
      d[0] -= w * z[i];
      return;
    }
    const ci = this.c[g * n + i];
    if (!(ci > 0)) return;
    const o = g * n * n + i * n, f = w / ci;
    for (let j = 0; j < n; j++) d[1 + j] += f * this.K[o + j];
    d[0] -= f * this._Kz(g, i);
  }

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

    // Faces: the flux node carries D, the linked species' fluxes and the reaction rates.
    for (let f = 0; f < this.nFaces; f++) this._face(f);

    // Terminals: ports (after every other term at their nodes, so a held level can read its
    // flux), then contacts, then the circuit rows of the floating ones.
    this._terminals(dt);
    if (this.flattening && dt === Infinity) this._flatRows();
    if (this.constrained && dt === Infinity && this.constraints.length > 0) this._applyConstraints();
    if (this.constrained && dt === Infinity && this.localConstraints.length > 0) this._applyLocalConstraints();
    this.transformed = this.combining && dt !== Infinity;
    if (this.transformed) this._chargeRows(dt);
  }

  // A conductor node: only its carrier's balance (the segment flux J has its own row, in
  // _segmentConductor). Its bulk stores nothing; a charged face's sheet sits in the edge node.
  _nodeConductor(g, dt) {
    const { n, M, res, c, cOld, z } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g], im = this.nodeConductor[g];
    const r = 1 + im, k = g * n + im;
    res[R[b * M + r]] += (v * (c[k] - cOld[k])) / dt;
    const f = this.sheetFace[g];
    if (f >= 0) this._j(b, r, this.blockOfFace[f], 0, (this.sheetSign[g] * this.model.grid.area[g]) / (z[im] * FARADAY * dt));
  }

  // Dilute nodes g0…g1 (ideal statistics, c = c_ref e^ζ): storage of each species and the space
  // charge in the Gauss row, written straight into the diagonal blocks by local index. A
  // dielectric is the case with no species.
  _nodesDilute(g0, g1, dt) {
    const { n, M, res, c, cOld, z, sys, loc, present, rhoFixed, blockOfNode } = this, R = this.rix, F = FARADAY;
    const JB = sys.B, sizes = sys.sizes, offB = sys.offB, vol = this.model.grid.vol;
    for (let g = g0; g <= g1; g++) {
      const b = blockOfNode[g], v = vol[g];
      const m = sizes[b], oB = offB[b], lb = b * M, p = loc[lb]; // p: φ's row, −1 where undefined
      if (this.combining && this.chargeNode[g] === 1 && dt !== Infinity) {
        // Storage and neutrality go in later, in _chargeRows; a neutral species' storage here.
        for (let i = 0; i < n; i++) if (present[g * n + i] && z[i] === 0) this._storage(g, i, dt);
        continue;
      }
      let q = rhoFixed[g], dq = 0;
      const late = this.combining && this.lateStorage[g] === 1 && dt !== Infinity; // (the ions' storage: in _contact)
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        if (!present[k]) continue;
        const ck = c[k], l = loc[lb + 1 + i];
        if (!(late && z[i] !== 0)) {
          res[R[lb + 1 + i]] += (v * (ck - cOld[k])) / dt;
          JB[oB + l * m + l] += (v * ck) / dt;
          if (p >= 0 && z[i] !== 0) JB[oB + l * m + p] += (-v * z[i] * ck) / dt;
        }
        q += F * z[i] * ck;
        dq += F * z[i] * z[i] * ck;
        if (p >= 0 && z[i] !== 0) JB[oB + p * m + l] += -v * F * z[i] * ck;
      }
      if (p >= 0) {
        res[R[lb]] -= v * q;
        JB[oB + p * m + p] += v * dq;
      }
    }
  }

  // A concentrated node (any statistics): ∂c_i/∂η_j = K_ij and ∂c_i/∂φ̂ = −(Kz)_i.
  _nodeConcentrated(g, dt) {
    const { n, M, res, c, z } = this, R = this.rix, F = FARADAY;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g];
    if (this.combining && this.chargeNode[g] === 1 && dt !== Infinity) {
      // Storage and neutrality go in later, in _chargeRows; a neutral species' storage here.
      for (let i = 0; i < n; i++) if (this.present[g * n + i] && z[i] === 0) this._storage(g, i, dt);
      return;
    }
    let q = this.rhoFixed[g], dq = 0;
    for (let i = 0; i < n; i++) {
      const k = g * n + i, r = 1 + i;
      if (!this.present[k]) continue;
      const Kz = this._Kz(g, i);
      this._storage(g, i, dt, Kz);
      q += F * z[i] * c[k];
      dq += F * z[i] * Kz;
      if (Kz !== 0) this._j(b, 0, b, r, -v * F * Kz); // K symmetric: ∂(Σ z c)/∂η_i = (Kz)_i
    }
    // (where no charge responds and no field reaches, φ isn't defined, nor an unknown)
    if (!this.phiUndefined[g]) {
      res[R[b * M]] -= v * q;
      this._j(b, 0, b, 0, v * dq);
    }
  }

  // Species i's storage at node g, v(c − c_old)/dt, in its balance row and the Jacobian
  // (∂c_i/∂η_j = K_ij, ∂c_i/∂φ̂ = −(Kz)_i; for ideal statistics K is diagonal, c_i).
  _storage(g, i, dt, Kz = this._Kz(g, i)) {
    const { n, M, res, c, cOld } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g], k = g * n + i, r = 1 + i;
    res[R[b * M + r]] += (v * (c[k] - cOld[k])) / dt;
    if (this.nodeIdeal[g]) this._j(b, r, b, r, (v * c[k]) / dt);
    else {
      const Kg = g * n * n;
      for (let j = 0; j < n; j++) {
        const Kij = this.K[Kg + i * n + j];
        if (Kij !== 0) this._j(b, r, b, 1 + j, (v * Kij) / dt);
      }
    }
    if (Kz !== 0 && !this.phiUndefined[g]) this._j(b, r, b, 0, (-v * Kz) / dt);
  }

  // Strictly neutral nodes on a transient step, solved in better-conditioned terms. There, a
  // change of φ̂ with every η_i shifted by z_i times it leaves every concentration as it was, so
  // storage and neutrality don't see it; but in (φ̂, η) they see it as pairs of huge entries that
  // cancel, and round-off in that cancellation (storage/flux ~ h²/(D·dt), times the range of
  // concentrations: a trace ion beside 3 M KCl) swamps the fluxes that do fix it. So, at each node
  // that stays neutral (see chargeNode), with storage and neutrality left out of assembly:
  // - rows: the balance of the most abundant charged species (by z²c) becomes Σ (z_i/z_k) × each
  //   balance, which without storage is current continuity;
  // - columns: the unknowns become φ̂' and η'_i = η_i − z_i φ̂ (newton() maps the update back),
  //   in which storage and neutrality have no φ̂' term at all;
  // - then storage (on the other balances) and neutrality go in, exactly, in those terms.
  // A start-of-step net charge (round-off in a solved state) isn't carried over: the neutrality
  // row holds the new state neutral.
  _chargeRows(dt) {
    const { n, M, res, c, cOld, z, sys, loc, present, termB, termC } = this, R = this.rix, F = FARADAY;
    const { A, B, C, sizes, offA, offB, offC } = sys, vol = this.model.grid.vol, last = this.nB - 1;
    for (const g of this.chargeNodes) {
      const b = this.blockOfNode[g], lb = b * M, m = sizes[b], mp = b > 0 ? sizes[b - 1] : 0, mn = b < last ? sizes[b + 1] : 0;
      const p = loc[lb], gn = g * n;
      let k = -1, best = -1;
      for (let i = 0; i < n; i++) {
        if (!present[gn + i] || z[i] === 0) continue;
        const w = z[i] * z[i] * c[gn + i];
        if (w > best) {
          k = i;
          best = w;
        }
      }
      if (k < 0 || p < 0) continue;
      // Rows: the pivot's balance becomes the z-weighted sum (all still without storage).
      const lk = loc[lb + 1 + k], rk = R[lb + 1 + k];
      for (let i = 0; i < n; i++) {
        if (i === k || !present[gn + i] || z[i] === 0) continue;
        const li = loc[lb + 1 + i], w = z[i] / z[k], ri = R[lb + 1 + i];
        for (let q = 0, o = offA[b] + lk * mp, s = offA[b] + li * mp; q < mp; q++) A[o + q] += w * A[s + q];
        for (let q = 0, o = offB[b] + lk * m, s = offB[b] + li * m; q < m; q++) B[o + q] += w * B[s + q];
        for (let q = 0, o = offC[b] + lk * mn, s = offC[b] + li * mn; q < mn; q++) C[o + q] += w * C[s + q];
        res[rk] += w * res[ri];
        for (let t = 0; t < termB.length; t++) termB[t][rk] += w * termB[t][ri];
      }
      // Columns: φ̂' = φ̂ with η fixed becomes φ̂ with η_i shifted by z_i, wherever node g's
      // unknowns appear (its own rows, its neighbours', a terminal's current).
      for (let i = 0; i < n; i++) {
        if (!present[gn + i] || z[i] === 0) continue;
        const li = loc[lb + 1 + i], zi = z[i];
        for (let r = 0, o = offB[b]; r < m; r++) B[o + r * m + p] += zi * B[o + r * m + li];
        if (b > 0) for (let r = 0, mr = mp, o = offC[b - 1]; r < mr; r++) C[o + r * m + p] += zi * C[o + r * m + li];
        if (b < last) for (let r = 0, mr = mn, o = offA[b + 1]; r < mr; r++) A[o + r * m + p] += zi * A[o + r * m + li];
        for (let t = 0; t < termC.length; t++) termC[t][R[lb]] += zi * termC[t][R[lb + 1 + i]];
      }
      // Storage on the other charged balances, and neutrality: no φ̂' terms, exactly.
      const v = vol[g], oB = offB[b];
      let q = this.rhoFixed[g];
      if (this.nodeIdeal[g]) {
        for (let i = 0; i < n; i++) {
          if (!present[gn + i] || z[i] === 0) continue;
          const ci = c[gn + i], li = loc[lb + 1 + i];
          if (i !== k) {
            res[R[lb + 1 + i]] += (v * (ci - cOld[gn + i])) / dt;
            B[oB + li * m + li] += (v * ci) / dt;
          }
          q += F * z[i] * ci;
          B[oB + p * m + li] += -v * F * z[i] * ci;
        }
      } else {
        for (let i = 0; i < n; i++) {
          if (!present[gn + i]) continue;
          if (i !== k && z[i] !== 0) this._storage(g, i, dt, 0);
          q += F * z[i] * c[gn + i];
          const Kz = this._Kz(g, i); // (nonzero for a neutral species only on a shared lattice)
          if (Kz !== 0) this._j(b, 0, b, 1 + i, -v * F * Kz);
        }
      }
      res[R[lb]] -= v * q;
      // The pivot's row (charge continuity) left storage out, the ions' net charge being constant
      // under neutrality. Against a capacitance it isn't: neutrality makes it −(aσ + ρ_fixed)/F,
      // so their storage is that less what the step started with, Σ z c_old, over dt (in units
      // of the pivot's balance).
      if (this.qgTerms[g].length > 0) {
        const f = -v / (FARADAY * z[k] * dt);
        let qOld = this.rhoFixed[g];
        for (let i = 0; i < n; i++) if (present[gn + i]) qOld += F * z[i] * cOld[gn + i];
        res[rk] += f * (this.qg[g] + qOld);
        B[offB[b] + lk * m + p] += -f * this.qgSlope[g] * this.VT;
        for (const [kq, w] of this.qgTerms[g]) termB[kq][rk] += f * this.portArea[kq - 2][w] * this.model.ports[kq - 2].capacitance.C;
      }
    }
  }

  // Newton's update in the transformed unknowns (see _chargeRows) back in η: η_i = η'_i + z_i φ̂'.
  _untransform(delta) {
    const { n, M, z, loc, present } = this, R = this.rix;
    for (const g of this.chargeNodes) {
      const lb = this.blockOfNode[g] * M;
      if (loc[lb] < 0) continue;
      const dphi = delta[R[lb]];
      for (let i = 0; i < n; i++) if (present[g * n + i] && z[i] !== 0) delta[R[lb + 1 + i]] += z[i] * dphi;
    }
  }

  // Bulk reactions at node g: r = k_f Π c_R^ν · (−expm1(−a)), with a = A/RT from the
  // (compensated) η. Reactants are consumed (+v·ν·r in their balance), products made (−v·ν·r).
  _bulkReactions(g, rxs) {
    const { n, M, u, uLo, res, c } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g];
    const dr = this.dA;
    for (let x = 0; x < rxs.length; x++) {
      const rx = rxs[x];
      const law = rx.srh?.[this.nodeMaterial[g]];
      if (law) {
        // SRH: dr is ∂rate/∂slot directly, from its sensitivities to ln n, ln p and a.
        let aHi = rx.fixedA, aLo = 0;
        for (let p = 0; p < rx.sp.length; p++) {
          aHi += rx.nu[p] * u[b * M + 1 + rx.sp[p]];
          aLo += rx.nu[p] * uLo[b * M + 1 + rx.sp[p]];
        }
        const { rate, Cn, Cp, Ca } = this._srhRate(law, c[g * n + law.n.i], c[g * n + law.p.i], aHi + aLo);
        dr.fill(0);
        this._dlnc(g, law.n.i, Cn, dr);
        this._dlnc(g, law.p.i, Cp, dr);
        for (let p = 0; p < rx.sp.length; p++) dr[1 + rx.sp[p]] += Ca * rx.nu[p];
        for (let p = 0; p < rx.sp.length; p++) {
          const row = 1 + rx.sp[p], w = v * rx.nu[p];
          res[R[b * M + row]] += w * rate;
          for (let k = 0; k < M; k++) if (dr[k] !== 0) this._j(b, row, b, k, w * dr[k]);
        }
        continue;
      }
      let P = rx.generation ? rx.kfNode[g] * this.generationScale : rx.kfNode[g], aHi = rx.fixedA, aLo = 0;
      dr.fill(0); // ∂ ln P / ∂slot
      // Participants: reactants with ν > 0, then products with ν < 0.
      for (let p = 0; p < rx.sp.length; p++) {
        const i = rx.sp[p], nu = rx.nu[p];
        if (nu > 0) {
          const ci = c[g * n + i];
          P *= nu === 1 ? ci : powi(ci, nu);
          this._dlnc(g, i, nu, dr);
        }
        aHi += nu * u[b * M + 1 + i];
        aLo += nu * uLo[b * M + 1 + i];
      }
      const a = aHi + aLo;
      const f = -Math.expm1(-a); // 1 − e^{−a}
      const rate = P * f;
      const Pd = P * (1 - f); // P·df/da, df/da = e^{−a}
      // ∂rate/∂slot = rate·∂lnP + P·f′·∂a
      for (let k = 0; k < M; k++) dr[k] *= rate;
      for (let p = 0; p < rx.sp.length; p++) dr[1 + rx.sp[p]] += Pd * rx.nu[p];
      for (let p = 0; p < rx.sp.length; p++) {
        const row = 1 + rx.sp[p], w = v * rx.nu[p];
        res[R[b * M + row]] += w * rate;
        for (let k = 0; k < M; k++) if (dr[k] !== 0) this._j(b, row, b, k, w * dr[k]);
      }
    }
  }

  // Displacement along a segment where φ is defined: D = −ε(φ_R − φ_L)/h into the two Gauss rows.
  _segmentDisplacement(s, bL, bR, mat, h) {
    const { M, u, res, sys, loc } = this, R = this.rix;
    const k = this.phiUndefined[s] ? 0 : (mat.epsr * EPS0 * this.VT) / h; // ε = 0: no displacement
    if (k === 0) return;
    const D = -k * (u[bR * M] - u[bL * M]);
    const mL = sys.sizes[bL], mR = sys.sizes[bR], pL = loc[bL * M], pR = loc[bR * M];
    res[R[bL * M]] += D;
    res[R[bR * M]] -= D;
    sys.B[sys.offB[bL] + pL * mL + pL] += k;
    sys.C[sys.offC[bL] + pL * mR + pR] -= k;
    sys.A[sys.offA[bR] + pR * mL + pL] -= k;
    sys.B[sys.offB[bR] + pR * mR + pR] += k;
  }

  // Scharfetter–Gummel fluxes along a dilute region's segments (nodes g0…g1), written
  // by local index. The flux
  // N = g[B(Δ)c_L − B(−Δ)c_R] is rewritten with B(−Δ) = B(Δ)e^Δ and c_R e^Δ = c_L e^{Δη} as
  // N = −g·B(Δ)·c_L·expm1(Δη). This is precise relative to the quasi-Fermi difference Δη, so tiny
  // fluxes (e.g. majority carriers carrying a small current) don't vanish in the cancellation of
  // two huge drift and diffusion terms.
  _segmentsDilute(g0, g1, region, mat, vel) {
    const { n, M, u, uLo, res, c, z, sys, loc } = this, R = this.rix, F = FARADAY;
    const { A: JA, B: JB, C: JC, sizes, offA, offB, offC } = sys;
    const { segLength, segArea } = this.model.grid;
    for (let s = g0; s < g1; s++) {
      const bL = s + region, bR = bL + 1, h = segLength[s] / segArea[s]; // (fluxes are totals through the cross-section)
      const mL = sizes[bL], mR = sizes[bR], pL = loc[bL * M], pR = loc[bR * M];
      const phiL = u[bL * M], phiR = u[bR * M];
      for (let i = 0; i < n; i++) {
        if (!mat.present[i] || mat.D[i] === 0) continue;
        const r = 1 + i, zi = z[i];
        const cL = c[s * n + i];
        const g = mat.D[i] / h;
        const pe = (vel * h) / mat.D[i]; // advection: a Péclet shift of the drift potential
        const d = zi * (phiR - phiL) - pe;
        const deta = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]) - pe;
        const E = Math.expm1(deta);
        const gBc = g * bernoulli(d) * cL;
        const N = -gBc * E;
        const dNdd = -g * bernoulliDerivative(d) * cL * E;
        const dNdEtaL = gBc;
        const dNdEtaR = -gBc * (E + 1);
        const dNdPhiL = -zi * dNdd + zi * gBc * E; // via Δ, and via c_L ∝ e^{−zφ̂_L}
        const dNdPhiR = zi * dNdd;
        res[R[bL * M + r]] += N;
        res[R[bR * M + r]] -= N;
        // Row r of block bL (couplings to itself in B, to bR in C) and of bR (to bL in A).
        const lL = loc[bL * M + r], lR = loc[bR * M + r];
        const oBL = offB[bL] + lL * mL, oCL = offC[bL] + lL * mR, oAR = offA[bR] + lR * mL, oBR = offB[bR] + lR * mR;
        JB[oBL + lL] += dNdEtaL;
        JC[oCL + lR] += dNdEtaR;
        JA[oAR + lL] -= dNdEtaL;
        JB[oBR + lR] -= dNdEtaR;
        if (zi !== 0 && pL >= 0) {
          JB[oBL + pL] += dNdPhiL;
          JC[oCL + pR] += dNdPhiR;
          JA[oAR + pL] -= dNdPhiL;
          JB[oBR + pR] -= dNdPhiR;
        }
      }
    }
  }

  // A face: its flux node's rows are the interface laws (φ law, species links, reactions), and
  // its fluxes enter the two edge nodes' balances.
  _face(f) {
    const { model, n, M, u, uLo, res, z, VT } = this, R = this.rix, F = FARADAY;
    const grid = model.grid, interfaces = model.interfaces;
    const bf = this.blockOfFace[f], bL = bf - 1, bR = bf + 1;
    const itf = interfaces[f];
    // φ law: pinned jump (dipole), Helmholtz capacitor, or no charge at all (neutral).
    const law = itf.phi.type;
    const condSide = itf.conductor && law !== 'neutral' ? itf.conductor.side : null;
    if (condSide) {
      // Against a metal: the other side's φ is tied to the metal's Fermi level V_F = V_T η/z.
      const im = itf.conductor.i, bm = condSide === 'left' ? bL : bR, bo = condSide === 'left' ? bR : bL;
      const sg = condSide === 'left' ? 1 : -1;
      const vf = (u[bm * M + 1 + im] + uLo[bm * M + 1 + im]) / z[im]; // V_F / V_T
      const gap = vf - itf.zeroCharge / VT - (u[bo * M] + uLo[bo * M]); // (V_F − zeroCharge − φ_edge)/V_T
      if (law === 'pinned') {
        res[R[bf * M]] = gap;
        this._j(bf, 0, bm, 1 + im, 1 / z[im]);
        this._j(bf, 0, bo, 0, -1);
      } else {
        // D toward the other side = C (V_F − zeroCharge − φ_edge); toward +x that's sg times it.
        const kC = itf.phi.C * VT;
        res[R[bf * M]] = u[bf * M] - sg * kC * gap;
        this._j(bf, 0, bf, 0, 1);
        this._j(bf, 0, bm, 1 + im, (-sg * kC) / z[im]);
        this._j(bf, 0, bo, 0, sg * kC);
      }
    } else if (law === 'neutral') {
      // D = 0 (not an unknown); the jump is whatever each side's neutrality needs
    } else {
      const jump = u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) - itf.dipole / VT;
      if (law === 'pinned') {
        res[R[bf * M]] = jump;
        this._j(bf, 0, bR, 0, 1);
        this._j(bf, 0, bL, 0, -1);
      } else {
        // D = −C (φ_R − φ_L − dipole): displacement toward +x drops across the layer.
        const kC = itf.phi.C * VT;
        res[R[bf * M]] = u[bf * M] + kC * jump;
        this._j(bf, 0, bf, 0, 1);
        this._j(bf, 0, bR, 0, kC);
        this._j(bf, 0, bL, 0, -kC);
      }
    }
    // D_f enters each side's Gauss row; a metal side holds it as surface carriers instead. The
    // flux node's unknowns are per area, so what they carry is that times the face's area.
    const gL = grid.regionEnd[f], gR = gL + 1, Af = grid.area[gL];
    if (this.nodeConductor[gL] < 0) {
      res[R[bL * M]] += Af * u[bf * M];
      this._j(bL, 0, bf, 0, Af);
    }
    if (this.nodeConductor[gR] < 0) {
      res[R[bR * M]] -= Af * (u[bf * M] + itf.sheetCharge);
      this._j(bR, 0, bf, 0, -Af);
    }
    for (let i = 0; i < n; i++) {
      const r = 1 + i, o = bf * M + r;
      const type = itf.links[i].type;
      const deta = u[bL * M + r] - u[bR * M + r] + (uLo[bL * M + r] - uLo[bR * M + r]); // η_L − η_R
      if (type === 'blocked') continue;
      if (type === 'equilibrium') {
        res[R[o]] = -deta; // μ̄ continuous
        this._j(bf, r, bR, r, 1);
        this._j(bf, r, bL, r, -1);
      } else if (type === 'permeability') {
        // Electrodiffusion through a thin membrane in a constant field: Scharfetter–Gummel across
        // the face with D/h → P, its drift taken over the face's whole jump in the standard level
        // zφ̂ + μ°/RT − ln c_ref (the φ jump of a capacitive face, and any step between the two
        // materials). Between like solutions it's Goldman–Hodgkin–Katz's flux,
        //   N = P·zu·(c_L − c_R e^{zu})/(e^{zu} − 1),  u = (φ_R − φ_L)/V_T,
        // and it's exactly zero where μ̄ is level.
        const gl = grid.regionEnd[f], gr = gl + 1, kl = gl * n + i, kr = gr * n + i, zi = z[i], P = itf.links[i].P;
        const shift = this.mu0hat[kr] - this.mu0hat[kl] - Math.log(this.cRef[kr] / this.cRef[kl]);
        const d = zi * (u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M])) + shift;
        const E = Math.expm1(-deta); // η_R − η_L
        const gBc = P * bernoulli(d) * this.c[kl];
        const dNdd = -P * bernoulliDerivative(d) * this.c[kl] * E;
        res[R[o]] = u[o] + gBc * E;
        this._j(bf, r, bf, r, 1);
        this._j(bf, r, bL, r, -gBc);
        this._j(bf, r, bR, r, gBc * (E + 1));
        if (zi !== 0) {
          this._j(bf, r, bL, 0, zi * dNdd - zi * gBc * E);
          this._j(bf, r, bR, 0, -zi * dNdd);
        }
      } else {
        // conductance: J = G (V_L − V_R), V = V_T η / z  ⇒  N = G V_T (η_L − η_R) / (z² F)
        const gG = (itf.links[i].G * VT) / (z[i] * z[i] * F);
        res[R[o]] = u[o] - gG * deta;
        this._j(bf, r, bf, r, 1);
        this._j(bf, r, bL, r, -gG);
        this._j(bf, r, bR, r, gG);
      }
      res[R[bL * M + r]] += Af * u[o];
      this._j(bL, r, bf, r, Af);
      res[R[bR * M + r]] -= Af * u[o];
      this._j(bR, r, bf, r, -Af);
    }
    itf.reactions.forEach((rx, k) => this._faceReaction(rx, k, f, bf, bL, bR, Af));
  }

  // Scharfetter–Gummel with non-ideal statistics. The excess ex = ζ − ln(c/c_ref) acts as an
  // extra potential, linear along the segment like φ, so Δ = zΔφ̂ + Δex and
  //   N = −(D/h)·B(Δ)·c_L·expm1(η_R − η_L),
  // still exactly zero at equilibrium. c_L and ex depend on every ζ at their node through K.
  _segmentConcentrated(s, bL, bR, mat, h, vel) {
    const R = this.rix;
    const { n, M, u, uLo, res, c, z, K, ex, jL, jR } = this;
    const gL = s, gR = s + 1, KL = gL * n * n, KR = gR * n * n;
    for (let i = 0; i < n; i++) {
      if (!mat.present[i] || mat.D[i] === 0) continue;
      const r = 1 + i, zi = z[i];
      const cL = c[gL * n + i], cR = c[gR * n + i];
      const g = mat.D[i] / h;
      const pe = (vel * h) / mat.D[i];
      const d = zi * (u[bR * M] - u[bL * M]) + ex[gR * n + i] - ex[gL * n + i] - pe;
      const deta = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]) - pe;
      const E = Math.expm1(deta);
      const B = bernoulli(d), Bp = bernoulliDerivative(d);
      const gBc = g * B * cL;
      const N = -gBc * E;
      const dNdd = -g * Bp * cL * E;
      // Left node: ∂N/∂ζ_Lj = −gE(B + B′)K_L,ij − dNdd·δ_ij, plus the direct η_L,i term.
      const fL = -g * E * (B + Bp);
      jL[0] = -fL * this._Kz(gL, i);
      for (let j = 0; j < n; j++) jL[1 + j] = fL * K[KL + i * n + j];
      jL[r] += -dNdd + gBc * (E + 1);
      // Right node: ∂N/∂ζ_Rj = dNdd·(δ_ij − K_R,ij/c_R), plus the direct η_R,i term.
      jR.fill(0);
      if (cR > 0) {
        const fR = -dNdd / cR;
        jR[0] = -fR * this._Kz(gR, i);
        for (let j = 0; j < n; j++) jR[1 + j] = fR * K[KR + i * n + j];
      } else {
        jR[0] = dNdd * zi; // K/c → δ as c → 0
        jR[r] -= dNdd;
      }
      jR[r] += dNdd - gBc * (E + 1);
      res[R[bL * M + r]] += N;
      res[R[bR * M + r]] -= N;
      for (let k = 0; k < M; k++) {
        if (jL[k] !== 0) {
          this._j(bL, r, bL, k, jL[k]);
          this._j(bR, r, bL, k, -jL[k]);
        }
        if (jR[k] !== 0) {
          this._j(bL, r, bR, k, jR[k]);
          this._j(bR, r, bR, k, -jR[k]);
        }
      }
    }
  }

  // Ohmic conduction in a metal, in mixed form. The flux J (slot 0 of the left node) obeys
  //   η_R − η_L + J/g = 0,  g = σ RT/(z²F² h)   (Ohm's law, J = −g Δη),
  // and enters the two nodes' carrier balances as outflow and inflow.
  _segmentConductor(s, bL, bR, mat, h) {
    const R = this.rix;
    const { M, u, uLo, res, z } = this;
    const i = mat.conductor.i, r = 1 + i;
    const g = (mat.conductor.sigma * this.model.RT) / (z[i] * z[i] * FARADAY * FARADAY * h);
    const J = u[bL * M];
    res[R[bL * M]] = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]) + J / g;
    this._j(bL, 0, bR, r, 1);
    this._j(bL, 0, bL, r, -1);
    this._j(bL, 0, bL, 0, 1 / g);
    res[R[bL * M + r]] += J;
    res[R[bR * M + r]] -= J;
    this._j(bL, r, bL, 0, 1);
    this._j(bR, r, bL, 0, -1);
  }

  // Eddy mixing: N = −(D_mix/RT) P ∇μ̄ with P = C − (Cz)(Cz)ᵀ/(zᵀCz), C = diag(c). It mixes
  // composition without carrying current (zᵀP = 0) and vanishes exactly at equilibrium. On a
  // segment, N_i = −(D_mix/h) Σ_j P̄_ij Δη_j, with P̄ from the logarithmic mean of each c, so a
  // neutral species gets exactly −D_mix Δc/h. Only mobile species (D > 0) take part.
  _segmentMixing(s, bL, bR, Dm, h, mat) {
    const R = this.rix;
    const { n, M, u, uLo, res, c, z } = this;
    const gL = s, gR = s + 1, k0 = Dm / h;
    const on = mat.mobile ?? (mat.mobile = [...Array(n).keys()].filter((i) => mat.present[i] && mat.D[i] > 0));
    if (on.length === 0) return;
    const w = this.mixWork ?? (this.mixWork = [0, 1, 2, 3].map(() => new Float64Array(n)));
    const [cb, dLa, dLb, de] = w;
    let S = 0, Q = 0;
    for (const k of on) {
      const a = c[gL * n + k], b = c[gR * n + k];
      // l = ln(c_R/c_L), from the compensated potentials where statistics are ideal
      const l = this.nodeIdeal[gL]
        ? u[bR * M + 1 + k] - u[bL * M + 1 + k] + (uLo[bR * M + 1 + k] - uLo[bL * M + 1 + k]) - z[k] * (u[bR * M] - u[bL * M])
        : Math.log(b / a);
      // L = a·f(l), f = expm1(l)/l; ∂L/∂a = f − f′, ∂L/∂b = a f′/b
      let f, fp;
      if (Math.abs(l) < 1e-3) {
        f = 1 + l / 2 + (l * l) / 6 + (l * l * l) / 24;
        fp = 0.5 + l / 3 + (l * l) / 8 + (l * l * l) / 30;
      } else {
        const em = Math.expm1(l);
        f = em / l;
        fp = (l * (em + 1) - em) / (l * l);
      }
      cb[k] = a * f;
      dLa[k] = f - fp;
      dLb[k] = b > 0 ? (a * fp) / b : 0.5;
      de[k] = u[bR * M + 1 + k] - u[bL * M + 1 + k] + (uLo[bR * M + 1 + k] - uLo[bL * M + 1 + k]);
      S += z[k] * z[k] * cb[k];
      Q += z[k] * cb[k] * de[k];
    }
    const q = S > 0 ? Q / S : 0;
    const { jL, jR } = this;
    const dcL = this.dA, dcR = this.dB; // ∂c_k/∂slot at each end
    for (const i of on) {
      const r = 1 + i;
      const N = -k0 * cb[i] * (de[i] - z[i] * q);
      jL.fill(0);
      jR.fill(0);
      for (const j of on) {
        // direct dependence on Δη_j
        const P = cb[i] * ((i === j ? 1 : 0) - (S > 0 ? (z[i] * z[j] * cb[j]) / S : 0));
        jL[1 + j] += k0 * P;
        jR[1 + j] -= k0 * P;
      }
      // through c̄_k: ∂a_i/∂c̄_k = δ_ik (Δη_i − z_i q) − c̄_i z_i (z_k Δη_k − q z_k²)/S
      for (const k of on) {
        let da = (i === k ? de[i] - z[i] * q : 0) - (S > 0 ? (cb[i] * z[i] * (z[k] * de[k] - q * z[k] * z[k])) / S : 0);
        if (da === 0) continue;
        da *= -k0;
        this._dc(gL, k, dcL);
        this._dc(gR, k, dcR);
        for (let t = 0; t < M; t++) {
          jL[t] += da * dLa[k] * dcL[t];
          jR[t] += da * dLb[k] * dcR[t];
        }
      }
      res[R[bL * M + r]] += N;
      res[R[bR * M + r]] -= N;
      for (let t = 0; t < M; t++) {
        if (jL[t] !== 0) {
          this._j(bL, r, bL, t, jL[t]);
          this._j(bR, r, bL, t, -jL[t]);
        }
        if (jR[t] !== 0) {
          this._j(bL, r, bR, t, jR[t]);
          this._j(bR, r, bR, t, -jR[t]);
        }
      }
    }
  }

  // ∂c_k/∂slot at node g into d (slot 0 = φ̂, 1 + j = η_j).
  _dc(g, k, d) {
    const { n } = this;
    d.fill(0);
    if (this.nodeIdeal[g]) {
      const ck = this.c[g * n + k];
      d[1 + k] = ck;
      d[0] = -this.z[k] * ck;
      return;
    }
    const o = g * n * n + k * n;
    for (let j = 0; j < n; j++) d[1 + j] = this.K[o + j];
    d[0] = -this._Kz(g, k);
  }

  // A bulk reaction's rate at node g (mol/(m³·s), forward), as assembled: for the solution.
  bulkRate(rx, g) {
    const { n, M, u, uLo, c } = this, b = this.blockOfNode[g];
    let aHi = rx.fixedA, aLo = 0;
    for (let p = 0; p < rx.sp.length; p++) {
      aHi += rx.nu[p] * u[b * M + 1 + rx.sp[p]];
      aLo += rx.nu[p] * uLo[b * M + 1 + rx.sp[p]];
    }
    const law = rx.srh?.[this.nodeMaterial[g]];
    if (law) return this._srhRate(law, c[g * n + law.n.i], c[g * n + law.p.i], aHi + aLo).rate;
    let P = rx.generation ? rx.kfNode[g] * this.generationScale : rx.kfNode[g];
    for (let p = 0; p < rx.sp.length; p++) if (rx.nu[p] > 0) P *= rx.nu[p] === 1 ? c[g * n + rx.sp[p]] : powi(c[g * n + rx.sp[p]], rx.nu[p]);
    return P * -Math.expm1(-(aHi + aLo));
  }

  // The larger of a bulk reaction's two one-way rates at node g (mol/(m³·s)): the scale its net
  // rate is read against (in equilibrium the net rate is their round-off).
  bulkOneWay(rx, g) {
    const { n, M, u, uLo, c } = this, b = this.blockOfNode[g];
    let a = rx.fixedA;
    for (let p = 0; p < rx.sp.length; p++) a += rx.nu[p] * (u[b * M + 1 + rx.sp[p]] + uLo[b * M + 1 + rx.sp[p]]);
    const law = rx.srh?.[this.nodeMaterial[g]];
    if (law) return this._srhRate(law, c[g * n + law.n.i], c[g * n + law.p.i], a).oneWay;
    let P = rx.generation ? rx.kfNode[g] * this.generationScale : rx.kfNode[g];
    for (let p = 0; p < rx.sp.length; p++) if (rx.nu[p] > 0) P *= rx.nu[p] === 1 ? c[g * n + rx.sp[p]] : powi(c[g * n + rx.sp[p]], rx.nu[p]);
    return P * Math.max(1, Math.exp(-a));
  }

  // SRH kinetics (see srhLaw in device.js): the rate, and its sensitivities to ln n, ln p and
  // a = A/RT. With no n₁ given, the trap is midgap: n₁ = p₁ = n_i = √(n p e^{−a}).
  _srhRate(law, cn, cp, a) {
    const { tn, tp } = law, lnnp = Math.log(cn) + Math.log(cp);
    let n1, kN = 0, kP = 0, kA = 0; // ln n₁ = kN ln n + kP ln p + kA a + const
    if (Number.isNaN(law.n1)) {
      n1 = Math.exp((lnnp - a) / 2);
      kN = kP = 0.5;
      kA = -0.5;
    } else n1 = law.n1;
    const p1 = Math.exp(lnnp - a - Math.log(n1)), np = cn * cp;
    const den = tp * (cn + n1) + tn * (cp + p1), f = -Math.expm1(-a);
    const rate = (np * f) / den;
    const dDn = tp * cn + tp * n1 * kN + tn * p1 * (1 - kN);
    const dDp = tp * n1 * kP + tn * cp + tn * p1 * (1 - kP);
    const dDa = tp * n1 * kA + tn * p1 * (-1 - kA);
    return { rate, oneWay: (np / den) * Math.max(1, Math.exp(-a)), Cn: rate * (1 - dDn / den), Cp: rate * (1 - dDp / den), Ca: (np / den) * (1 - f) - (rate * dDa) / den };
  }

  // A face reaction (see normalizeFaceReactions). Its rate r_k is an unknown of the face block,
  // with the row r_k − rate(u_L, u_R) = 0, and each participant's edge node takes ν·r_k (made
  // there when ν > 0). The rate couples the two edge nodes only through the face block between
  // them, which keeps the system block-tridiagonal.
  _faceReaction(rx, k, f, bf, bL, bR, Af) {
    const { n, M, u, uLo, c, res } = this;
    const R = this.rix;
    const grid = this.model.grid;
    const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
    const al = rx.alpha;
    const dL = this.dA.fill(0), dR = this.dB.fill(0); // ∂ ln(prefactor)/∂slot, per side
    let pref = rx.k0, aHi = rx.fixedA, aLo = 0, rate;
    for (const p of rx.part) {
      const g = p.side ? gR : gL, b = p.side ? bR : bL, o = b * M + 1 + p.i;
      aHi -= p.nu * u[o];
      aLo -= p.nu * uLo[o];
      if (rx.srh || rx.vmax || this.nodeConductor[g] === p.i) continue; // a conductor's carrier has activity 1
      const e = p.nu < 0 ? -p.nu * (1 - al) : p.nu * al;
      pref *= powr(c[g * n + p.i] / this.cRef[g * n + p.i], e);
      this._dlnc(g, p.i, e, p.side ? dR : dL);
    }
    if (rx.vmax) {
      // Saturating: r = vmax Π (c/(c + K))^{|ν|} (1 − e^{−a}); ∂ln S/∂ln c = |ν| K/(c + K).
      let S = 1;
      rx.part.forEach((p, q) => {
        if (p.nu >= 0) return;
        const g = p.side ? gR : gL, cc = c[g * n + p.i], K = rx.K[q];
        S *= powi(cc / (cc + K), -p.nu);
        this._dlnc(g, p.i, (-p.nu * K) / (cc + K), p.side ? dR : dL);
      });
      const a = aHi + aLo, back = Math.exp(-a);
      rate = rx.vmax * S * -Math.expm1(-a);
      for (let s = 0; s < M; s++) {
        dL[s] *= rate;
        dR[s] *= rate;
      }
      for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= rx.vmax * S * back * p.nu;
    } else if (rx.srh) {
      // SRH, with n and p each at its own side's edge node.
      const { n: pn, p: pp } = rx.srh, gn = pn.side ? gR : gL, gq = pp.side ? gR : gL;
      const s = this._srhRate(rx.srh, c[gn * n + pn.i], c[gq * n + pp.i], aHi + aLo);
      rate = s.rate;
      this._dlnc(gn, pn.i, s.Cn, pn.side ? dR : dL);
      this._dlnc(gq, pp.i, s.Cp, pp.side ? dR : dL);
      for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= s.Ca * p.nu;
    } else {
      const { g: bv, gp } = bvFactor(aHi + aLo, al);
      rate = pref * bv;
      // ∂rate/∂slot = rate·∂ln(prefactor) + pref·g′·∂a, with ∂a/∂η_i = −ν
      for (let s = 0; s < M; s++) {
        dL[s] *= rate;
        dR[s] *= rate;
      }
      for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= pref * gp * p.nu;
    }
    const slot = 1 + n + k, o = bf * M + slot;
    res[R[o]] = u[o] - rate;
    this._j(bf, slot, bf, slot, 1);
    for (let s = 0; s < M; s++) {
      if (dL[s] !== 0) this._j(bf, slot, bL, s, -dL[s]);
      if (dR[s] !== 0) this._j(bf, slot, bR, s, -dR[s]);
    }
    for (const p of rx.part) {
      const b = p.side ? bR : bL;
      res[R[b * M + 1 + p.i]] -= Af * p.nu * u[o];
      this._j(b, 1 + p.i, bf, slot, -Af * p.nu);
    }
  }

  // Every terminal's terms: its voltage held or floating, how the residual depends on it
  // (termB), its current into the device (termI) and how that depends on the state (termC).
  _terminals(dt) {
    this._refreshSources();
    for (let k = 0; k < this.terms.length; k++) {
      this.termB[k].fill(0);
      this.termC[k].fill(0);
      this.termDI[k] = 0;
      this.termI[k] = 0;
    }
    this.qg.fill(0);
    this.qgSlope.fill(0);
    this.model.ports.forEach((port, k) => this._port(port, this.portFlux[k], 2 + k, dt, false));
    this.model.ports.forEach((port, k) => this._port(port, this.portFlux[k], 2 + k, dt, true));
    this._contact('left', dt);
    this._contact('right', dt);
    for (const k of this.floating) this._circuit(k);
  }

  // A floating terminal's circuit law: I − I_set = 0 (driven by a current), or I − (V_src − V)/R = 0
  // (a source behind a resistance).
  _circuit(k) {
    const d = this.terms[k].drive, src = sourceAt(d.src, this.sourceTime, this.sourceBefore);
    if (d.kind === 'I') this.termRes[k] = this.termI[k] - src;
    else {
      this.termRes[k] = this.termI[k] - (src - this.termV[k]) / d.R;
      this.termDI[k] += 1 / d.R;
    }
  }

  // into += w × (the Jacobian row of block b, slot s), over compact columns.
  _captureRow(b, s, w, into) {
    const l = this.loc[b * this.M + s];
    if (l < 0 || w === 0) return;
    const { sys, nB } = this, { sizes, offA, offB, offC, offX } = sys, m = sizes[b];
    for (const [X, nb, off] of [['A', b - 1, offA], ['B', b, offB], ['C', b + 1, offC]]) {
      if (nb < 0 || nb >= nB) continue;
      const mc = sizes[nb], o = off[b] + l * mc, arr = sys[X], c0 = offX[nb];
      for (let c = 0; c < mc; c++) if (arr[o + c] !== 0) into[c0 + c] += w * arr[o + c];
    }
    return m;
  }

  // An internal port (terminal k): its exchange with each node of its window, as a source per
  // volume. A held ('equilibrium') level replaces the node's balance row; the source is then that
  // row's residual, read just before. Its current into the device is Σ z F × the sources. In two
  // passes over the ports: everything else, then the held levels (`holds`), so that a level reads
  // every other port's terms at its nodes, and none is added into its row after.
  _port(port, flux, k, dt, holds) {
    const R = this.rix;
    const { n, M, u, uLo, res, z, VT } = this;
    const F = FARADAY, vol = this.model.grid.vol;
    const B = this.termB[k], C = this.termC[k];
    if (!holds) flux.fill(0);
    for (let i = 0; i < n; i++) {
      const link = port.species[i];
      if (link.type === 'blocked' || (link.type === 'equilibrium') !== holds) continue;
      const r = 1 + i, target = this.portEta(port, i), zF = z[i] * F;
      for (const g of port.nodes) {
        const b = this.blockOfNode[g], o = b * M + r, v = vol[g];
        const deta = target - (u[o] + uLo[o]); // (μ̄_out − μ̄)/RT
        if (link.type === 'equilibrium') {
          if (g === 0 || g === this.nNodes - 1) continue; // a device end node's level is its contact's business
          flux[i] += res[R[o]];
          this._captureRow(b, r, zF, C);
          this._replaceRow(b, r);
          this._j(b, r, b, r, 1);
          res[R[o]] = -deta;
          B[R[o]] += -z[i] / VT;
          continue;
        }
        // conductance: s = G V_T (η_out − η)/(z² F); exchange: s = k (η_out − η)
        const kk = link.type === 'conductance' ? (link.G * VT) / (z[i] * z[i] * F) : link.k;
        res[R[o]] -= v * kk * deta;
        this._j(b, r, b, r, v * kk);
        flux[i] += v * kk * deta;
        if (z[i] !== 0) {
          B[R[o]] += (-v * kk * z[i]) / VT;
          C[R[o]] += -zF * v * kk;
          this.termDI[k] += (zF * v * kk * z[i]) / VT;
        }
      }
    }
    if (holds) {
      for (let i = 0; i < n; i++) this.termI[k] += z[i] * F * flux[i];
      return;
    }
    port.reactions.forEach((rx, x) => this._portReaction(port, rx, this.portArea[k - 2], this.portRates[k - 2][x], flux, k));
    if (port.capacitance) this._portCapacitance(port, k, dt);
    if (port.surface.length > 0 && Number.isFinite(dt)) this._surfaceStorage(port, dt);
  }

  // A port's reaction at each node of its window: Butler–Volmer per area, as at a face, against
  // the electrode's carrier at the port's level (activity 1), times the electrode's area per
  // volume. Each node's rate is explicit (no unknown of its own): it enters the balances of the
  // species it makes and consumes, and the port's current.
  _portReaction(port, rx, area, rates, flux, k) {
    const { n, M, u, uLo, c, res, z, VT, nSurf, th } = this, R = this.rix, F = FARADAY, vol = this.model.grid.vol, RT = this.model.RT;
    const B = this.termB[k], C = this.termC[k], al = rx.alpha, d = this.dA, surf = port.surface;
    // The carrier's part of the affinity, and how it moves with the port's voltage.
    let aV = rx.fixedA, daV = 0, q = 0; // q: charge the forward reaction brings into the device, per event
    for (const p of rx.part) {
      if (p.side === 1) {
        aV -= p.nu * this.portEta(port, p.i);
        daV -= (p.nu * z[p.i]) / VT;
      } else if (p.side === 0) q += p.nu * z[p.i];
    }
    const slot = (p) => (p.side === 2 ? 1 + n + p.s : 1 + p.i);
    port.nodes.forEach((g, w) => {
      const b = this.blockOfNode[g], s = vol[g] * area[w];
      d.fill(0);
      let pref = rx.k0, aHi = aV, aLo = 0;
      for (const p of rx.part) {
        if (p.side === 1) continue;
        const o = b * M + slot(p);
        aHi -= p.nu * u[o];
        aLo -= p.nu * uLo[o];
        const e = p.nu < 0 ? -p.nu * (1 - al) : p.nu * al;
        if (p.side === 2) {
          // a surface species at activity θ/θ₀ = e^ζ
          pref *= Math.exp(e * (u[o] + uLo[o] - surf[p.s].mu0 / RT));
          d[1 + n + p.s] += e;
        } else {
          pref *= powr(c[g * n + p.i] / this.cRef[g * n + p.i], e);
          this._dlnc(g, p.i, e, d);
        }
      }
      if (rx.bare) {
        // on bare metal only: the free fraction θ₀ = 1 − Σθ, with ∂ln θ₀/∂η_s = −θ_s
        for (let x = 0; x < surf.length; x++) d[1 + n + x] -= th[g * nSurf + x];
        pref *= this.th0[g];
      }
      const { g: bv, gp } = bvFactor(aHi + aLo, al);
      const rate = pref * bv;
      rates[w] = rate;
      // ∂rate/∂slot = rate·∂ln(prefactor) + pref·g′·∂a, with ∂a/∂η = −ν
      for (let t = 0; t < M; t++) d[t] *= rate;
      for (const p of rx.part) if (p.side !== 1) d[slot(p)] -= pref * gp * p.nu;
      const dV = pref * gp * daV;
      for (const p of rx.part) {
        if (p.side === 1) continue;
        // the region's species per volume of the window; a surface's per area of electrode
        const row = slot(p), o = b * M + row, f = p.side === 2 ? 1 : s;
        res[R[o]] -= f * p.nu * rate;
        if (p.side === 0) flux[p.i] += s * p.nu * rate;
        for (let t = 0; t < M; t++) if (d[t] !== 0) this._j(b, row, b, t, -f * p.nu * d[t]);
        B[R[o]] -= f * p.nu * dV;
      }
      for (let t = 0; t < M; t++) if (d[t] !== 0) C[R[b * M + t]] += q * F * s * d[t];
      this.termDI[k] += q * F * s * dV;
    });
  }

  // A capacitance spread through a port's window: per area of electrode, σ = C (V − zeroCharge − φ)
  // on the port's side, so each node's charge balance (its φ row) gains vol·a·σ, and the port
  // passes the charging current d(Σ vol·a·σ)/dt (none in a steady state). An end node whose φ its
  // contact sets is left to the contact.
  _portCapacitance(port, k, dt) {
    const { M, u, uLo, res, VT } = this, R = this.rix, vol = this.model.grid.vol, area = this.portArea[k - 2];
    const { C: Cs, zeroCharge } = port.capacitance, V = this.termV[k], B = this.termB[k], Ct = this.termC[k];
    const sets = (side) => ['pinned', 'bulk'].includes(this.model.contacts[side].phi.type);
    let Q = 0, dQdV = 0;
    port.nodes.forEach((g, w) => {
      if ((g === 0 && sets('left')) || (g === this.nNodes - 1 && sets('right')) || this.phiUndefined[g]) return;
      const b = this.blockOfNode[g], s = vol[g] * area[w] * Cs;
      const q = s * (V - zeroCharge - VT * (u[b * M] + uLo[b * M]));
      this.qg[g] += area[w] * Cs * (V - zeroCharge - VT * (u[b * M] + uLo[b * M]));
      this.qgSlope[g] += area[w] * Cs;
      Q += q;
      dQdV += s;
      res[R[b * M]] -= q;
      this._j(b, 0, b, 0, s * VT);
      B[R[b * M]] -= s;
      if (Number.isFinite(dt)) Ct[R[b * M]] += (-s * VT) / dt;
    });
    this.portQ[k - 2] = Q;
    if (Number.isFinite(dt)) {
      this.termI[k] += (Q - this.portQStart[k - 2]) / dt;
      this.termDI[k] += dQdV / dt;
    }
  }

  // An electrode surface's storage, per area of electrode: Γ (θ_s − θ_s,old)/dt in each surface
  // species' row, with ∂θ_s/∂η_x = θ_s (δ_sx − θ_x).
  _surfaceStorage(port, dt) {
    const { n, M, res, nSurf, th, thOld } = this, R = this.rix, surf = port.surface, G = surf[0].capacity;
    for (const g of port.nodes) {
      const b = this.blockOfNode[g];
      // Near full coverage θ − θ_old cancels (both ~1, the change ~θ₀): the dominant species'
      // change comes from the small ones instead, Δθ = −Δθ₀ − Σ_others Δθ.
      let top = 0;
      for (let q = 1; q < surf.length; q++) if (th[g * nSurf + q] > th[g * nSurf + top]) top = q;
      for (let q = 0; q < surf.length; q++) {
        const tq = th[g * nSurf + q];
        let change = tq - thOld[g * nSurf + q];
        if (q === top && tq > 0.5) {
          change = -(this.th0[g] - this.th0Old[g]);
          for (let x = 0; x < surf.length; x++) if (x !== q) change -= th[g * nSurf + x] - thOld[g * nSurf + x];
        }
        res[R[b * M + 1 + n + q]] += (G * change) / dt;
        for (let x = 0; x < surf.length; x++) this._j(b, 1 + n + q, b, 1 + n + x, ((G * tq) / dt) * ((q === x ? 1 : 0) - th[g * nSurf + x]));
      }
    }
  }

  // One contact (terminal 0 or 1): record the flux through its outer face and the current into
  // the device (read from the end box before anything here is added: whatever the box needs
  // comes through the contact), then add its exchange terms, its equilibrium links and its φ law.
  _contact(side, dt) {
    const R = this.rix;
    const { model, n, M, u, uLo, res, c, z, VT } = this;
    const F = FARADAY;
    const ct = model.contacts[side];
    const k = side === 'left' ? 0 : 1, Vt = this.termV[k];
    const B = this.termB[k], C = this.termC[k];
    const g = side === 'left' ? 0 : this.nNodes - 1;
    const b = this.blockOfNode[g];
    const sgn = side === 'left' ? 1 : -1;
    const flux = this.contactFlux[side];
    const Ac = model.grid.area[g]; // the links' laws are per area: what they pass is that times A here
    const dyn = Number.isFinite(dt);

    // Before anything is added, each balance residual is the flux into the device here. (A
    // closed end solved in charge rows has its ions' storage left out until later: none of them
    // passes, so their flux is zero.)
    const closedIons = this.combining && dyn && this.chargeNode[g] === 1;
    let I = 0;
    for (let i = 0; i < n; i++) {
      const active = this.loc[b * M + 1 + i] >= 0 && !(closedIons && z[i] !== 0);
      const nIn = active ? res[R[b * M + 1 + i]] : 0;
      flux[i] = sgn * nIn; // (toward +x)
      if (active && z[i] !== 0) {
        I += z[i] * F * nIn;
        this._captureRow(b, 1 + i, z[i] * F, C);
      }
    }
    // The ions' storage, left out until now (see lateStorage): in each balance, and in what
    // each species brings in, but not in the current.
    if (this.combining && dyn && this.lateStorage[g] === 1) {
      const v = model.grid.vol[g], { cOld } = this;
      for (let i = 0; i < n; i++) {
        if (z[i] === 0 || this.loc[b * M + 1 + i] < 0) continue;
        const ci = c[g * n + i], st = (v * (ci - cOld[g * n + i])) / dt;
        res[R[b * M + 1 + i]] += st;
        flux[i] += sgn * st;
        this._j(b, 1 + i, b, 1 + i, (v * ci) / dt);
        this._j(b, 1 + i, b, 0, (-v * z[i] * ci) / dt);
      }
    }

    // Conductance links: J (toward the device) = G (V_out − V_i), V_out = V_t + offset; exchange
    // links (neutral species): N_in = k (μ_out − μ)/RT.
    for (let i = 0; i < n; i++) {
      const link = ct.species[i];
      const o = b * M + 1 + i;
      if (link.type === 'exchange') {
        res[R[o]] -= Ac * link.k * (link.mu / model.RT - (u[o] + uLo[o]));
        this._j(b, 1 + i, b, 1 + i, Ac * link.k);
      } else if (link.type === 'conductance') {
        const Vi = (VT * (u[o] + uLo[o])) / z[i];
        res[R[o]] -= (Ac * link.G * (Vt + link.offset - Vi)) / (z[i] * F);
        this._j(b, 1 + i, b, 1 + i, (Ac * link.G * VT) / (z[i] * z[i] * F));
        B[R[o]] += (-Ac * link.G) / (z[i] * F);
      }
    }

    // Equilibrium links: Dirichlet on the known outside level.
    for (let i = 0; i < n; i++) {
      const link = ct.species[i];
      if (link.type !== 'equilibrium') continue;
      const o = b * M + 1 + i;
      this._replaceRow(b, 1 + i);
      this._j(b, 1 + i, b, 1 + i, 1);
      res[R[o]] = u[o] + uLo[o] - this.contactEta(side, i);
      if (z[i] !== 0) B[R[o]] += -z[i] / VT;
    }

    // φ law (none where φ is undefined: nothing at the end node responds to it). The displacement
    // into the device, D_in = sgn·D, adds (D_in − D_in at the step's start)/dt to the current.
    const link = ct.phi, Dstart = sgn * this.contactDStart[side];
    if (this.phiUndefined[g]) {
      this.contactD[side] = 0;
    } else if (link.type === 'capacitive') {
      // Gate, or metal across a Stern layer, at φ_g = V_t − zeroCharge. D toward +x.
      const Din = Ac * link.C * (Vt - link.zeroCharge - VT * u[b * M]);
      res[R[b * M]] -= Din;
      this._j(b, 0, b, 0, Ac * link.C * VT);
      B[R[b * M]] += -Ac * link.C;
      this.contactD[side] = sgn * Din;
      if (dyn) {
        I += (Din - Dstart) / dt;
        C[R[b * M]] += (-Ac * link.C * VT) / dt;
        this.termDI[k] += (Ac * link.C) / dt;
      }
    } else if (link.type === 'pinned' || link.type === 'bulk') {
      // The Poisson residual is the outside's charge, D_in.
      const Din = res[R[b * M]];
      this.contactD[side] = sgn * Din;
      if (dyn && !(this.combining && this.lateStorage[g] === 1)) {
        I += (Din - Dstart) / dt;
        this._captureRow(b, 0, 1 / dt, C);
      }
      this._replaceRow(b, 0);
      if (link.type === 'pinned') {
        // φ_edge = V_t − zeroCharge (the C → ∞ limit).
        this._j(b, 0, b, 0, 1);
        res[R[b * M]] = u[b * M] + uLo[b * M] - (Vt - link.zeroCharge) / VT;
        B[R[b * M]] += -1 / VT;
      } else {
        // Plain bulk: local neutrality.
        let q = this.rhoFixed[g], dq = 0;
        for (let i = 0; i < n; i++) {
          const kk = g * n + i;
          if (!this.present[kk]) continue;
          const Kz = this._Kz(g, i);
          q += F * z[i] * c[kk];
          dq -= F * z[i] * Kz;
          if (Kz !== 0) this._j(b, 0, b, 1 + i, F * Kz);
        }
        this._j(b, 0, b, 0, dq);
        res[R[b * M]] = q;
      }
    } else {
      this.contactD[side] = 0;
    }
    this.termI[k] = I;
  }




  // Steady solves: each flat stretch's level held at its contact's throughout (see
  // flatStretches), in place of its balances, as an equilibrium link holds it at the end node.
  _flatRows() {
    const { M, u, uLo, res, z, VT, model } = this, R = this.rix;
    for (const st of this.flatStretches) {
      const i = st.species, side = st.leftOpen ? 'left' : 'right', k = side === 'left' ? 0 : 1;
      const link = model.contacts[side].species[i];
      const level = z[i] === 0 ? link.mu / model.RT : (z[i] * (this.termV[k] + link.offset)) / VT;
      for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
        const b = this.blockOfNode[g], o = b * M + 1 + i;
        if (this.loc[o] < 0) continue;
        this._replaceRow(b, 1 + i);
        for (const B of this.termB) B[R[o]] = 0; // (the row's old terminal terms: the end node's own link)
        this._j(b, 1 + i, b, 1 + i, 1);
        res[R[o]] = u[o] + uLo[o] - level;
        if (z[i] !== 0) this.termB[k][R[o]] += -z[i] / VT;
      }
    }
  }

  // Zero one row of the Jacobian (all three blocks), ready to be replaced.
  _replaceRow(b, rs) {
    const l = this.loc[b * this.M + rs];
    if (l < 0) return;
    const { sys, nB } = this, sz = sys.sizes, m = sz[b];
    sys.B.fill(0, sys.offB[b] + l * m, sys.offB[b] + (l + 1) * m);
    if (b > 0) sys.A.fill(0, sys.offA[b] + l * sz[b - 1], sys.offA[b] + (l + 1) * sz[b - 1]);
    if (b < nB - 1) sys.C.fill(0, sys.offC[b] + l * sz[b + 1], sys.offC[b] + (l + 1) * sz[b + 1]);
  }

  // Scale every row by its largest Jacobian entry (in place, residual too).
  _equilibrate() {
    const { nB, sys, res, termB } = this;
    const { A, B, C, sizes, offA, offB, offC, offX } = sys;
    for (let b = 0; b < nB; b++) {
      const m = sizes[b], mp = b > 0 ? sizes[b - 1] : 0, mn = b < nB - 1 ? sizes[b + 1] : 0;
      for (let r = 0; r < m; r++) {
        const oa = offA[b] + r * mp, ob = offB[b] + r * m, oc = offC[b] + r * mn;
        let mx = 0;
        for (let k = 0; k < mp; k++) mx = Math.max(mx, Math.abs(A[oa + k]));
        for (let k = 0; k < m; k++) mx = Math.max(mx, Math.abs(B[ob + k]));
        for (let k = 0; k < mn; k++) mx = Math.max(mx, Math.abs(C[oc + k]));
        if (mx === 0) {
          this.rowScale[offX[b] + r] = 1;
          continue;
        }
        const s = 1 / mx;
        this.rowScale[offX[b] + r] = s;
        for (let k = 0; k < mp; k++) A[oa + k] *= s;
        for (let k = 0; k < m; k++) B[ob + k] *= s;
        for (let k = 0; k < mn; k++) C[oc + k] *= s;
        res[offX[b] + r] *= s;
        for (let k = 0; k < termB.length; k++) termB[k][offX[b] + r] *= s; // (each terminal's ∂/∂V)
      }
    }
  }

  // A Newton solve refined by GMRES: the factorised Jacobian as the preconditioner, and J·x
  // taken from the residual itself, by a central difference (η moved in the state's low word),
  // as the impedance does. The assembled J holds a flux's dependence on η at its two ends as two
  // entries; where they're huge (an inversion layer's 0.1 nm cells, a conductance ~1e11 against
  // the layer's own storage and the trickle from the bulk), eliminating them loses the layer's
  // overall level to round-off, and Newton rattles there. The residual takes η differences
  // first, in double-double, so it keeps that level. delta and deltaV come in as the plain
  // solve's and leave refined (in the transformed unknowns, as _solveBordered gives them); false
  // if refining didn't help, leaving them as they came.
  _refine(dt, delta, deltaV) {
    const { res, termRes, termB, termC, termDI, rowScale, floating: fl } = this;
    const N = this.sys.size, K = fl.length;
    const save = { res: Float64Array.from(res), termRes: Float64Array.from(termRes), termB: termB.map((b) => Float64Array.from(b)), termC: termC.map((c) => Float64Array.from(c)), termDI: Float64Array.from(termDI), transformed: this.transformed };
    const u0 = Float64Array.from(this.u), u0Lo = Float64Array.from(this.uLo), V0 = Float64Array.from(this.termV);
    const restoreArrays = () => {
      res.set(save.res);
      termRes.set(save.termRes);
      save.termB.forEach((b, k) => termB[k].set(b));
      save.termC.forEach((c, k) => termC[k].set(c));
      termDI.set(save.termDI);
      this.transformed = save.transformed;
    };
    const restoreState = () => {
      this.u.set(u0);
      this.uLo.set(u0Lo);
      this.termV.set(V0);
    };
    const fd = this.sysFD ?? (this.sysFD = new BlockTridiagonal(this.sys.n, this.sys.sizes));
    const at = (x, sgn, hh, out) => {
      restoreState();
      const v = this.dWork2 ?? (this.dWork2 = new Float64Array(N + 1));
      for (let j = 0; j < N; j++) v[j] = x[j];
      v[N] = 0;
      if (save.transformed) this._untransform(v);
      this._addToState(v, sgn * hh);
      fl.forEach((k, a) => (this.termV[k] += sgn * hh * x[N + a]));
      const keep = this.sys;
      this.sys = fd;
      try {
        this.assemble(dt);
      } finally {
        this.sys = keep;
      }
      for (let j = 0; j < N; j++) out[j] = res[j] * rowScale[j];
      fl.forEach((k, a) => (out[N + a] = termRes[k]));
    };
    const plus = new Float64Array(N + K), minus = new Float64Array(N + K);
    const rhs = new Float64Array(N + 1), out = new Float64Array(N + 1), outV = new Float64Array(K);
    const precondition = (r) => {
      restoreArrays();
      for (let j = 0; j < N; j++) rhs[j] = r[j];
      fl.forEach((k, a) => (termRes[k] = r[N + a]));
      this._solveBordered(rhs, out, outV, []);
      const z = new Float64Array(N + K);
      for (let j = 0; j < N; j++) z[j] = out[j];
      for (let a = 0; a < K; a++) z[N + a] = outV[a];
      return z;
    };
    const op = (x) => {
      let mx = 0;
      for (let j = 0; j < N + K; j++) mx = Math.max(mx, Math.abs(x[j]));
      if (mx === 0) return new Float64Array(N + K);
      const hh = 1e-6 / mx;
      at(x, 1, hh, plus);
      at(x, -1, hh, minus);
      const jx = new Float64Array(N + K);
      for (let j = 0; j < N + K; j++) jx[j] = (plus[j] - minus[j]) / (2 * hh);
      return precondition(jx);
    };
    try {
      const b = new Float64Array(N + K);
      for (let j = 0; j < N; j++) b[j] = save.res[j];
      fl.forEach((k, a) => (b[N + a] = save.termRes[k]));
      const pb = precondition(b), x = Float64Array.from(pb);
      // A lost mode or two, GMRES finds in a couple of iterations. The refined update is kept if
      // its correction dwarfs the plain update (which, along a lost mode, comes out as noise)
      // and is under a thermal unit. A correction no bigger than the plain update means the
      // system is near-singular more broadly (a slow ion over a long step), its residual there
      // round-off: solving it exactly would only chase that, so newton() stops refining.
      const x0 = Float64Array.from(x);
      const ok = gmresReal(op, pb, x, { m: 8, restarts: 1, tol: 1e-8 });
      let c = 0, size = 0;
      for (let j = 0; j < N + K; j++) {
        c = Math.max(c, Math.abs(x[j] - x0[j]));
        size = Math.max(size, Math.abs(x0[j]));
      }
      if (!ok || !(c > 10 * size) || !(c < 1)) {
        x.set(x0);
        return false;
      }
      for (let j = 0; j < N; j++) delta[j] = x[j];
      for (let a = 0; a < K; a++) deltaV[a] = x[N + a];
      return true;
    } finally {
      restoreState();
      restoreArrays();
      this.computeConcentrations();
    }
  }

  // Largest update among the potential-like unknowns (φ̂ and η at grid nodes).
  _maxPotentialStep(delta, deltaV = this.deltaV) { // (deltaV null: the device's potentials only)
    const { n, M } = this;
    const R = this.rix; // (the sink entry of delta is 0)
    let mx = 0;
    if (deltaV) for (let a = 0; a < this.floating.length; a++) mx = Math.max(mx, Math.abs(deltaV[a]) / this.VT);
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      if (this.nodeConductor[g] < 0) mx = Math.max(mx, Math.abs(delta[R[b * M]])); // (a metal's slot 0 is a flux)
      for (let i = 0; i < n; i++) {
        if (this.present[g * n + i]) mx = Math.max(mx, Math.abs(delta[R[b * M + 1 + i]]));
      }
      const k = this.surfPort[g]; // an electrode surface's coverages, as η
      if (k >= 0) for (let q = 0; q < this.model.ports[k].surface.length; q++) mx = Math.max(mx, Math.abs(delta[R[b * M + 1 + n + q]]));
    }
    return mx;
  }

  /**
   * Newton iteration for one backward-Euler step from the current cOld.
   * @returns {{converged: boolean, iterations: number, history: number[]}}
   */
  newton(dt, { maxIter = 60, tol = 1e-10, maxStep = 10 } = {}) {
    const { delta, res } = this;
    const deltaV = this.deltaV ?? (this.deltaV = new Float64Array(this.floating.length));
    const history = [], ownHistory = [];
    // Refined solves (see _refine), once the plain ones stall.
    let refine = false, refined = 0, noRefine = false;
    // A steady solve that fails: record how nearly singular the system was, and where.
    const fail = (r) => {
      if (dt === Infinity) this._noteConditioning();
      return r;
    };
    for (let it = 1; it <= maxIter; it++) {
      try {
        this.assemble(dt);
      } catch (err) {
        // A statistics model can fail far from the solution (e.g. Debye–Hückel beyond its range).
        if (!(err instanceof SolverError)) throw err;
        return fail({ converged: false, iterations: it, history, error: err.message });
      }
      this._equilibrate();
      let rmax = 0;
      for (let k = 0; k < this.SINK; k++) rmax = Math.max(rmax, Math.abs(res[k]));
      try {
        this._factor();
      } catch (err) {
        return fail({ converged: false, iterations: it, history, error: err.message });
      }
      try {
        const pins = this.constrained && dt === Infinity ? this.constraints : [];
        deltaV.fill(0);
        this._solveBordered(res, delta, deltaV, pins);
        if (refine && pins.length === 0 && !this._refine(dt, delta, deltaV)) {
          refine = false; // (nothing lost: the plain solves are as good as any here)
          noRefine = true;
        }
        if (this.transformed) this._untransform(delta);
      } catch (err) {
        if (!(err instanceof SolverError)) throw err;
        return fail({ converged: false, iterations: it, history, error: err.message });
      }
      // Convergence counts every unknown, floating terminal voltages too. Damping and the
      // divergence check count only the device's own potentials: a terminal voltage enters
      // linearly (conductance links, held levels, a capacitive face), so a large swing in one is
      // safe, and damping it would throttle the whole update (a current-driven port switching
      // off, its voltage collapsing by hundreds of volts, crawled 10 thermal units at a time).
      const step = this._maxPotentialStep(delta, deltaV);
      const own = this._maxPotentialStep(delta, null);
      history.push(step);
      ownHistory.push(own);
      if (!Number.isFinite(step)) return fail({ converged: false, iterations: it, history, error: 'non-finite update' });
      // Give up early on clear divergence; the caller will take a smaller step instead.
      if (own > 1e4 || (it > 6 && own > 10 * ownHistory[0])) {
        this.computeConcentrations();
        return fail({ converged: false, iterations: it, history, error: 'diverging' });
      }
      const alpha = own > maxStep ? maxStep / own : 1;
      this._addToState(delta, -alpha);
      this.floating.forEach((k, a) => (this.termV[k] -= alpha * deltaV[a]));
      if (alpha === 1 && step < tol) {
        this.computeConcentrations();
        return { converged: true, iterations: it, history, residual: rmax };
      }
      // Converged as far as round-off allows: the updates are already tiny (below 1e-6 thermal
      // units, ~26 nV) and have stopped shrinking. A badly conditioned system's floor can sit
      // above tol: a strictly neutral material on a short step, where φ is fixed only through
      // fluxes that the storage term dwarfs, rattles at ~1e-9 after converging quadratically.
      // (Raising this floor lets a weakly held population drift: a MOS capacitor without a
      // channel port then showed a DC leak through its oxide.)
      const [p1, p2] = [history[history.length - 2], history[history.length - 3]], floor = 1e-6;
      // Stalled near the solution (not shrinking quadratically): the factorised Jacobian may
      // have lost a slow mode to cancellation, and round-off's floor would accept it wrong.
      // Refine the solves from here.
      if (refine) refined++;
      else if (!noRefine && alpha === 1 && it >= 3 && step < 1e-4 && step > 0.25 * p1) refine = true;
      if (alpha === 1 && it >= 4 && (!refine || refined >= 1 || noRefine) && step < floor && p1 < floor && step > 0.25 * p1 && p1 > 0.25 * p2) {
        this.computeConcentrations();
        return { converged: true, iterations: it, history, residual: rmax, roundoff: true };
      }
    }
    this.computeConcentrations();
    return fail({ converged: false, iterations: maxIter, history });
  }

  // How many digits the last factorisation lost to cancellation, and where (the worst seen since
  // the steady solve began). Only for diagnosing failures: a saturated species, for instance,
  // can lose digits harmlessly.
  _noteConditioning() {
    let c;
    try {
      c = this.sys.cancellation();
    } catch {
      return;
    }
    if (!(c.digits > (this.conditioning?.digits ?? -1))) return;
    const { grid } = this.model;
    let x = NaN, where = '';
    for (let g = 0; g < this.nNodes; g++) {
      if (this.blockOfNode[g] !== c.block) continue;
      x = grid.x[g];
      where = this.model.regions[grid.nodeRegion[g]].name;
    }
    for (let f = 0; f < this.nFaces; f++) {
      if (this.blockOfFace[f] !== c.block) continue;
      x = grid.x[grid.regionEnd[f]];
      where = `interfaces[${f}]`;
    }
    this.conditioning = { digits: c.digits, x, where };
  }

  /**
   * Advance by dt, splitting the interval into halves (recursively) wherever Newton fails.
   * The end time is always t + dt; the result reports the substeps taken.
   */
  advance(dt, opts, depth = 0) {
    const r = this.step(dt, opts);
    if (r.converged || depth >= 30) return { ...r, substeps: 1 };
    const a = this.advance(dt / 2, opts, depth + 1);
    if (!a.converged) return a;
    const b = this.advance(dt / 2, opts, depth + 1);
    return {
      converged: b.converged,
      iterations: r.iterations + a.iterations + b.iterations,
      history: b.history,
      substeps: a.substeps + b.substeps,
    };
  }

  // u += scale·d in compensated arithmetic (Knuth two-sum, then renormalise hi/lo).
  _addToState(d, scale) {
    const { u, uLo, fullOf } = this;
    for (let j = 0; j < fullOf.length; j++) {
      const k = fullOf[j], a = u[k], b = scale * d[j];
      const s = a + b, bb = s - a;
      const err = a - (s - bb) + (b - bb);
      const lo = uLo[k] + err;
      const hi = s + lo;
      u[k] = hi;
      uLo[k] = lo - (hi - s);
    }
  }

  /**
   * One time step of size dt. Backward Euler by default; with `method: 'bdf2'`, variable-step
   * BDF2 once a previous step exists (and the step ratio is at most 2, for stability). On
   * failure the state is restored.
   *
   * BDF2's storage term (a0 c − (1+ω) c_n + ω²/(1+ω) c_{n−1})/dt is written in backward-Euler
   * form (c − c*)/(dt/a0), so assembly is shared: c* and the displacement histories replace
   * the start-of-step values, and dt/a0 replaces dt.
   */
  step(dt, opts = {}) {
    this.atSteady = false;
    const prev = this.history[0];
    const w = prev ? dt / prev.dt : 0;
    const bdf = opts.method === 'bdf2' && prev !== undefined && w <= 2;
    this.uPrev.set(this.u);
    this.uPrevLo.set(this.uLo);
    this.termVPrev.set(this.termV);
    this.computeConcentrations();
    const cN = Float64Array.from(this.c), thN = Float64Array.from(this.th), th0N = Float64Array.from(this.th0);
    // Contact displacement before the step: as it was at the end of the previous step (under
    // the parameters then), so a gate-voltage change shows up as displacement current.
    if (!this.contactDEnd || !this.portQEnd) {
      this._assembleBookkeeping(dt);
      this.contactDEnd = { ...this.contactD };
      this.portQEnd = Float64Array.from(this.portQ);
    }
    const DN = { ...this.contactDEnd }, QN = Float64Array.from(this.portQEnd);
    if (opts.guess) {
      // Start Newton from a predicted state (e.g. extrapolated from the history).
      this.u.set(opts.guess);
      this.uLo.fill(0);
    }
    let dtEff = dt;
    if (bdf) {
      const a0 = (1 + 2 * w) / (1 + w), b1 = (1 + w) / a0, b2 = (w * w) / (1 + w) / a0;
      const { cOld } = this;
      for (let k = 0; k < cOld.length; k++) cOld[k] = b1 * cN[k] - b2 * prev.c[k];
      for (let k = 0; k < thN.length; k++) this.thOld[k] = b1 * thN[k] - b2 * prev.th[k];
      for (let k = 0; k < th0N.length; k++) this.th0Old[k] = b1 * th0N[k] - b2 * prev.th0[k];
      this.contactDStart = { left: b1 * DN.left - b2 * prev.D.left, right: b1 * DN.right - b2 * prev.D.right };
      this.portQStart = QN.map((q, k) => b1 * q - b2 * prev.Q[k]);
      dtEff = dt / a0;
    } else {
      this.cOld.set(cN);
      this.thOld.set(thN);
      this.th0Old.set(th0N);
      this.contactDStart = DN;
      this.portQStart = QN;
    }
    // Implicit: sources at the step's end, as seen from within the step (before any jump there).
    // A step landing on a breakpoint ends on it exactly.
    const tEnd = this.landing !== undefined && Math.abs(this.time + dt - this.landing) <= 1e-9 * dt ? this.landing : this.time + dt;
    this.sourceTime = this.steady ? this.time : tEnd;
    this.sourceBefore = !this.steady;
    this.combining = true;
    let result;
    try {
      result = this.newton(dtEff, opts);
    } finally {
      this.combining = false;
    }
    result.bdf = bdf;
    if (result.converged) {
      this.history.unshift({ t: this.time, dt, u: Float64Array.from(this.uPrev), c: cN, th: thN, th0: th0N, D: DN, Q: QN });
      if (this.history.length > 3) this.history.length = 3;
      this.time = tEnd;
      this.lastDt = dtEff;
      this.contactDOld = this.contactDStart;
      this.portQOld = this.portQStart;
      this._accumulateBoundaryIntake(dtEff, cN);
    } else {
      this.u.set(this.uPrev);
      this.uLo.set(this.uPrevLo);
      this.termV.set(this.termVPrev);
      this.sourceTime = this.time;
      this.computeConcentrations();
    }
    this.sourceBefore = false;
    return result;
  }

  /**
   * Adaptive time stepping to tEnd: variable-step BDF2 (or backward Euler), with the local error
   * estimated against an explicit predictor through the previous states, and controlled to
   * `tol` thermal units per step in every potential (φ̂, each η, a surface's coverages as η). The first
   * step is checked by step doubling. Stops early when `budgetMs` of wall time is used, so an
   * animation can call it once per frame; the step size carries over between calls.
   */
  integrate(tEnd, opts = {}) {
    const { tol = 1e-3, dtMax = Infinity, budgetMs = Infinity, maxSteps = 100000, method = 'bdf2' } = opts;
    const clock = () => (globalThis.performance ? globalThis.performance.now() : Date.now());
    const start = clock();
    const sample = this._probes(opts.probes);
    const trace = { t: [], current: [], voltage: [], ...(sample ? { probes: sample.out } : {}) };
    let steps = 0, rejected = 0, iterations = 0, failed = false;
    // Already there (an animation frame with no time to add): nothing to do, and the step size
    // carried to the next call stays as it was.
    if (!(tEnd > this.time)) return { converged: true, done: true, steps: 0, rejected: 0, iterations: 0, trace };
    let dt = this.dtNext ?? opts.dt0 ?? (tEnd - this.time) * 1e-4;
    if (!(dt > 0)) dt = (tEnd - this.time) * 1e-4;
    // The shortest step worth trying: round-off relative to the time, or the device's fastest
    // time scale, whichever is shorter (a cold start can need steps far below a long run's 1e-14).
    const floor = Math.min(1e-14 * Math.max(tEnd, 1e-300), this.fastestTime() * 1e-2);
    const pred = new Float64Array(this.u.length), guess = new Float64Array(this.u.length);
    const firstErrs = []; // the first step's error at each size tried (see below)
    const factor = (err, p) => (err > 0 ? Math.min(2, Math.max(0.2, 0.9 * Math.exp(Math.log(tol / err) / (p + 1)))) : 2);
    while (this.time < tEnd) {
      if (steps + rejected >= maxSteps || clock() - start > budgetMs) break;
      // Within round-off of the end (a target a hair past a breakpoint just landed on, say):
      // snap to it rather than attempt a step of ~1e-13 s, which can't be resolved and, in a
      // strictly neutral device, fails outright.
      if (tEnd - this.time <= 1e-10 * Math.abs(tEnd)) {
        this.time = tEnd;
        break;
      }
      const remaining = tEnd - this.time;
      let h = Math.min(dt, dtMax);
      let clamped = h >= remaining * (1 - 1e-9);
      if (clamped) h = remaining;
      // Land on the waveforms' breakpoints, where a source's slope jumps.
      const tb = this.terms.reduce((m, t) => Math.min(m, nextBreakpoint(t.drive.src, this.time)), Infinity);
      const atBreak = tb - this.time <= h * (1 + 1e-9);
      if (atBreak) {
        h = tb - this.time;
        clamped = true;
      }
      this.landing = atBreak ? tb : clamped ? tEnd : undefined;
      const snap = this._snapshot();
      let err, p, half = null;
      if (this.history.length > 0) firstErrs.length = 0;
      if (this.history.length === 0) {
        // First step: backward Euler, checked against two half steps (whose result is kept).
        const full = this.step(h);
        iterations += full.iterations;
        if (full.converged) {
          const uFull = Float64Array.from(this.u);
          this._restore(snap);
          const a = this.step(h / 2);
          // The first half step is kept too, so it goes in the trace (a jump's charge is in it).
          half = a.converged ? { t: this.time, current: this._terminalCurrent(), voltage: this.termV[1] - this.termV[0], probes: sample?.read() } : null;
          const b = a.converged ? this.step(h / 2) : a;
          iterations += a.iterations + (b === a ? 0 : b.iterations);
          if (b.converged) {
            err = this._errorNorm(uFull);
            p = 1;
          }
        }
        if (err === undefined) {
          this._restore(snap);
          rejected++;
          dt = h / 4;
          if (dt < floor) { failed = true; break; }
          continue;
        }
      } else {
        // (Order 2 once three consistent states can check it.)
        const r = this.step(h, { method: this.history.length >= 2 ? method : 'be', guess: this._extrapolate(guess, this.time + h) });
        iterations += r.iterations;
        if (!r.converged) {
          rejected++;
          dt = h / 4;
          if (dt < floor) { failed = true; break; }
          continue;
        }
        p = this._predict(pred, r.bdf);
        err = this._errorNorm(pred) * p.scale;
        p = p.order;
      }
      // A first step whose error doesn't shrink with it: a jump in a level held at a boundary
      // (a bath's, after a voltage step) starts a profile self-similar in x/√t, so a first step
      // of any length errs alike, down to the grid's own diffusion time, below Newton's reach.
      // Three tries in a row, each within a factor 2 of the last, say so: the error is the
      // jump's, not the step's. Take the step (backward Euler damps it as it should) and grow.
      if (half && err > tol) firstErrs.push(err);
      const jump = half && firstErrs.length >= 3 && firstErrs.slice(-3).every((e, j, a) => j === 0 || (e > 0.5 * a[j - 1] && e < 2 * a[j - 1]));
      if (err > tol && !jump) {
        this._restore(snap);
        rejected++;
        dt = h * factor(err, p);
        continue;
      }
      steps++;
      // The state the first step started from needn't satisfy the algebraic equations (φ in a
      // strictly neutral region, an interface's unknowns, after a start or a jump), and they
      // jump in the first instant; a predictor through it would see that jump as error on every
      // step after. The history starts from the half step.
      if (half) this.history.length = 1;
      dt = clamped ? Math.max(dt, h * factor(err, p)) : h * factor(err, p);
      if (atBreak) this.history = []; // the solution's slope jumps here: restart the order
      if (half) {
        trace.t.push(half.t);
        trace.current.push(half.current);
        trace.voltage.push(half.voltage);
        if (sample) sample.push(half.probes);
      }
      trace.t.push(this.time);
      trace.current.push(this._terminalCurrent());
      trace.voltage.push(this.termV[1] - this.termV[0]);
      if (sample) sample.push(sample.read());
    }
    this.dtNext = dt;
    this.landing = undefined;
    return { converged: !failed, done: this.time >= tEnd, steps, rejected, iterations, trace };
  }

  // Quadratic (or linear) extrapolation of the state to time t, from the current state and the
  // two most recent history entries: a starting guess for Newton.
  _extrapolate(out, t) {
    const [e0, e1] = this.history;
    const t0 = this.time, u = this.u;
    if (!e0) return undefined;
    if (!e1) {
      const l = (t - t0) / (t0 - e0.t);
      for (let k = 0; k < out.length; k++) out[k] = u[k] + l * (u[k] - e0.u[k]);
      return out;
    }
    const t1 = e0.t, t2 = e1.t;
    const l0 = ((t - t1) * (t - t2)) / ((t0 - t1) * (t0 - t2));
    const l1 = ((t - t0) * (t - t2)) / ((t1 - t0) * (t1 - t2));
    const l2 = ((t - t0) * (t - t1)) / ((t2 - t0) * (t2 - t1));
    for (let k = 0; k < out.length; k++) out[k] = l0 * u[k] + l1 * e0.u[k] + l2 * e1.u[k];
    return out;
  }

  // Explicit predictor at the new time through the previous states (quadratic after BDF2,
  // linear after backward Euler), and the factor turning |u − pred| into the local error.
  _predict(pred, bdf) {
    const [e0, e1, e2] = this.history;
    const t = this.time, h = e0.dt;
    if (!e1) {
      pred.set(e0.u);
      return { scale: 1, order: 1 };
    }
    const hp = e1.dt;
    if (bdf && e2) {
      const hpp = e2.dt;
      // Lagrange through (t0, u0), (t1, u1), (t2, u2) at t.
      const t0 = e0.t, t1 = e1.t, t2 = e2.t;
      const l0 = ((t - t1) * (t - t2)) / ((t0 - t1) * (t0 - t2));
      const l1 = ((t - t0) * (t - t2)) / ((t1 - t0) * (t1 - t2));
      const l2 = ((t - t0) * (t - t1)) / ((t2 - t0) * (t2 - t1));
      for (let k = 0; k < pred.length; k++) pred[k] = l0 * e0.u[k] + l1 * e1.u[k] + l2 * e2.u[k];
      const w = h / hp;
      const Cc = (h * h * h * (1 + w) * (1 + w)) / (w * (1 + 2 * w));
      const Cp = h * (h + hp) * (h + hp + hpp);
      return { scale: Cc / (Cc + Cp), order: 2 };
    }
    const l0 = (t - e1.t) / (e0.t - e1.t), l1 = 1 - l0;
    for (let k = 0; k < pred.length; k++) pred[k] = l0 * e0.u[k] + l1 * e1.u[k];
    // Backward Euler's error against a linear predictor; after a BDF2 step without enough
    // history this overestimates, which is safe.
    return { scale: bdf ? 1 : h / (2 * h + hp), order: 1 };
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
    if (!Array.isArray(probes)) fail('advance: probes must be an array of { x, species, quantity, region }');
    const plan = probes.map((p, k) => {
      const path = `advance: probes[${k}]`;
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
        r = regions.findIndex((_, rr) => inside(rr));
        if (!(Number.isFinite(p.x) && r >= 0)) fail(`${path}.x must be a position in the device, in m (got ${p.x})`);
      }
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
      read: () => plan.map((q) => (1 - q.w) * at(q, q.g) + q.w * at(q, q.g + 1)),
      push: (values) => values.forEach((v, k) => out[k].push(v)),
    };
  }

  // Largest difference from ref among the potentials that carry state, in thermal units.
  _errorNorm(ref) {
    const { n, M, u } = this;
    let mx = 0;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      if (!this.phiUndefined[g]) mx = Math.max(mx, Math.abs(u[b * M] - ref[b * M]));
      for (let i = 0; i < n; i++) {
        if (this.present[g * n + i]) mx = Math.max(mx, Math.abs(u[b * M + 1 + i] - ref[b * M + 1 + i]));
      }
      const k = this.surfPort[g]; // an electrode surface's coverages, as η
      if (k >= 0) for (let q = 0; q < this.model.ports[k].surface.length; q++) mx = Math.max(mx, Math.abs(u[b * M + 1 + n + q] - ref[b * M + 1 + n + q]));
    }
    return mx;
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

  // Terminal current toward +x at the right contact, from the last assembly at this state.
  _terminalCurrent() {
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
  impedance(frequencies, { profiles = false, terminal = 'right' } = {}) {
    const { M, nB, sys, VT } = this;
    const kT = this.terms.findIndex((t) => t.name === terminal);
    if (kT < 0) throw new SolverError(`impedance: no terminal named '${terminal}' (${this.terms.map((t) => t.name).join(', ')})`);
    const dT = this.terms[kT].drive;
    if (dT.R > 0) throw new SolverError(`impedance: terminal '${terminal}' has a series resistance; that belongs to the external circuit`);
    const N = nB * M, fl = this.floating, K = fl.length;
    this.sourceTime = this.time;
    this.computeConcentrations();
    this.cOld.set(this.c);
    this.thOld.set(this.th);
    this.th0Old.set(this.th0);
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
    // J·v + (∂res/∂V_k)σ from the residual itself, by a central difference: the assembled J holds
    // a flux's dependence on η_L and η_R as two entries, and where both are huge (an inversion
    // layer) and v nearly uniform, J·v loses the flux to round-off; the residual takes the η
    // difference first, in double-double. So η moves in the low word, φ̂ and the flux unknowns
    // in the high one (φ̂ only scales expm1(Δη), near zero there).
    const u0 = Float64Array.from(this.u), uLo0 = Float64Array.from(this.uLo), V0 = Float64Array.from(this.termV);
    // The unknowns perturbed in the low word: a node's η (not φ̂, a metal's segment flux, or a
    // face's unknowns, which are read from the high word).
    const low = new Uint8Array(this.fullOf.length);
    const nodeBlock = new Uint8Array(nB);
    for (let g = 0; g < this.nNodes; g++) if (this.nodeConductor[g] < 0) nodeBlock[this.blockOfNode[g]] = 1;
    for (let j = 0; j < low.length; j++) {
      const full = this.fullOf[j];
      low[j] = nodeBlock[Math.floor(full / M)] && full % M !== 0 ? 1 : 0;
    }
    const resP = new Float64Array(this.res.length), resM = new Float64Array(this.res.length);
    const derivative = (v, k, sigma, out) => {
      let mx = Math.abs(sigma);
      for (let j = 0; j < v.length; j++) mx = Math.max(mx, Math.abs(v[j]));
      const h = 1e-4 / (mx || 1);
      for (const [sgn, into] of [[1, resP], [-1, resM]]) {
        this.u.set(u0);
        this.uLo.set(uLo0);
        this.termV.set(V0);
        for (let j = 0; j < v.length; j++) {
          const full = this.fullOf[j];
          if (low[j]) this.uLo[full] += sgn * h * v[j];
          else this.u[full] += sgn * h * v[j];
        }
        // A held terminal's voltage comes from its source at each assembly; a floating one's is
        // an unknown.
        if (k >= 0 && fl.includes(k)) this.termV[k] += sgn * h * sigma;
        else if (k >= 0) this.sourceOverride.set(k, V0[k] + sgn * h * sigma);
        this.assemble(Infinity);
        if (k >= 0) this.sourceOverride.delete(k);
        into.set(this.res);
      }
      for (let j = 0; j < out.length; j++) out[j] = (resP[j] - resM[j]) / (2 * h);
    };
    const restore = () => {
      this.u.set(u0);
      this.uLo.set(uLo0);
      this.termV.set(V0);
      this._refreshSources();
      this.computeConcentrations();
    };

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
      // Each solve uses that residual-based J·v (S·v is safe from the matrix): in an inversion
      // layer fed by minority carriers that number ~1e3 cm⁻³ in the bulk, the assembled J alone
      // gives the layer an exchange path many orders too fast.
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
      const op = ([vr, vi]) => {
        derivative(vr, -1, 0, jr);
        derivative(vi, -1, 0, ji);
        sx(S, vi, si);
        sx(S, vr, sr);
        const ar = new Float64Array(Nc), ai = new Float64Array(Nc);
        for (let j = 0; j < Nc; j++) {
          ar[j] = jr[j] - w * si[j];
          ai[j] = ji[j] + w * sr[j];
        }
        return precondition(ar, ai);
      };
      const solve = (br, bi, outr, outi, k, sigma) => {
        const rhs = k >= 0 ? precondition(Bt[k].map((v) => -sigma * v), zero) : [new Float64Array(Nc), new Float64Array(Nc)];
        outr.set(rhs[0]);
        outi.set(rhs[1]);
        gmres(op, rhs, [outr, outi]);
        restore();
      };
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
      let Zr, Zi;
      if (heldT) {
        // δI per volt: (C + iωC′)·δx + (∂I/∂V)(1 + iω′) at the terminal
        let [Ir, Ii] = dot(kT, w, yr, yi);
        Ir += DI[kT];
        Ii += w * DIs[kT];
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

  // The contact and port readouts (fluxes, displacements, the last segment's current) at the
  // current state, from only the boxes they're read from. Each of those boxes gets the same terms
  // in the same order as in assemble(), so the readouts are identical; other residuals and the
  // Jacobian are left partial (the next assemble() starts afresh). Concentrations must be current.
  _assembleBookkeeping(dt) {
    const { model } = this;
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

  // Add this step's contact fluxes (at the converged state) to each stretch's intake.
  // (A BDF2 step also moves each amount by Σ v (c* − c_n), its history term.)
  _accumulateBoundaryIntake(dt, cN) {
    this._assembleBookkeeping(dt);
    this.contactDEnd = { ...this.contactD };
    this.portQEnd = Float64Array.from(this.portQ);
    if (!Number.isFinite(dt)) return;
    const last = this.model.regions.length - 1;
    const { n, cOld } = this, vol = this.model.grid.vol;
    this.stretches.forEach((st, k) => {
      if (!st.connected) return;
      let q = 0;
      if (st.regions[0] === 0) q += this.contactFlux.left[st.species];
      if (st.regions[1] === last) q -= this.contactFlux.right[st.species];
      for (const p of st.ports) q += this.portFlux[p][st.species];
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
   * Steady state.
   * - If every species stretch is fed by a contact, nothing is conserved on its own, so the
   *   true steady equations (dt = ∞, no storage term) are solved directly: no slow modes to
   *   wait out, however slow the physics (e.g. exponentially scarce minority carriers).
   * - Otherwise backward-Euler steps at a huge dt, whose storage term pins each conserved
   *   amount exactly; dt keeps growing ×10 while the state still moves.
   * If Newton fails, dt ramps up from a small value (pseudo-transient continuation) instead.
   * The clock is not advanced, and open-system conservation bookkeeping restarts here.
   */
  // Steady state from the present one; `atSteady` says whether the state is one (until a step or
  // a change of drive).
  solveSteady(opts = {}) {
    this.flattening = true;
    try {
      return this._solveSteadyAll(opts);
    } finally {
      this.flattening = false;
    }
  }

  _solveSteadyAll(opts) {
    // A terminal that passes current only by charging (a gate, a capacitance with no species
    // through it) has no steady state under a current drive: it charges for ever, or at I = 0
    // keeps whatever charge it started with, which a steady solve doesn't know.
    const { species } = this.model, ions = (links) => links.some((l, i) => l.type !== 'blocked' && species[i].z !== 0);
    for (const k of this.floating) {
      const t = this.terms[k];
      if (t.drive.kind !== 'I') continue;
      const only = t.kind === 'port' ? (p) => p.reactions.length === 0 && !ions(p.species) : (c) => !ions(c.species);
      if (only(t.kind === 'port' ? this.model.ports[t.index] : this.model.contacts[t.side])) {
        throw new DeviceError(`${t.kind === 'port' ? `ports[${t.index}]` : `contacts.${t.side}`} passes current only by charging, so driven by a current it has no steady state; hold it at V (or a source behind R), or advance() in time`);
      }
    }
    const r = this._steadyFromHere(opts);
    this.atSteady = r.converged;
    return r;
  }

  _steadyFromHere(opts) {
    this.conditioning = null;
    this.computeConcentrations();
    this._captureLocal(); // immobile combinations keep, node by node, what they hold now
    this.sourceTime = this.time;
    // Continuation applies where the right terminal is held at a voltage, against a held left one.
    const [dl, dr] = [this.terms[0].drive, this.terms[1].drive];
    const held = (d) => d.kind === 'V' && !(d.R > 0);
    const target = sourceAt(dr.src, this.time), level = sourceAt(dl.src, this.time);
    const canContinue = opts.continuation !== false && held(dl) && held(dr) && target !== level;
    const direct = this._directSteady();
    const u0 = Float64Array.from(this.u), u0Lo = Float64Array.from(this.uLo), v0 = Float64Array.from(this.termV);
    // Where a direct solve applies and continuation is possible, don't spend long on the
    // pseudo-transient ramp: one direct attempt first.
    const quick = (canContinue || this.hasGeneration) && direct && opts.continuation !== false;
    let r = this._solveSteady(quick ? { ...opts, maxSteps: 1 } : opts);
    if (r.converged) this.solvedV = target;
    if (r.converged) return r;
    // Generation (e.g. light) holding the device far from equilibrium: ramp it up from nearly
    // nothing, each solve warm from the last.
    if (this.hasGeneration && opts.continuation !== false) {
      this.u.set(u0);
      this.uLo.set(u0Lo);
      this.termV.set(v0);
      this.computeConcentrations();
      const g = this._generationContinuation(opts, r);
      if (g.converged) {
        this.solvedV = target;
        return g;
      }
      r = g;
    }
    if (!canContinue) {
      if (quick) {
        this.u.set(u0);
        this.uLo.set(u0Lo);
        this.termV.set(v0);
        this.computeConcentrations();
        const full = this._solveSteady(opts); // the full pseudo-transient ramp
        return { ...full, steps: full.steps + r.steps, iterations: full.iterations + r.iterations };
      }
      return r;
    }
    // Source continuation: solve with both terminals level (consistent with a cold start), then
    // ramp the right terminal's voltage to its target in adaptive steps.
    const restart = () => {
      this.u.set(u0);
      this.uLo.set(u0Lo);
      this.termV.set(v0);
      this.computeConcentrations();
    };
    restart();
    // Ramp from the last converged voltage when the state is that solution (a warm start),
    // otherwise from level terminals.
    const from = this.solvedV !== undefined ? this.solvedV : level;
    const c = this._continuation(opts, from, target, r, from !== level);
    if (c.converged) {
      this.solvedV = target;
      return c;
    }
    restart();
    r = this._solveSteady(opts); // the full pseudo-transient ramp
    if (r.converged) this.solvedV = target;
    return { ...r, steps: r.steps + c.steps, iterations: r.iterations + c.iterations };
  }

  // Steady solves with the generation reactions' rates scaled from 1e-12 up to 1, ×100 a step
  // while they converge, smaller steps where they don't.
  _generationContinuation(opts, r) {
    const sub = this._directSteady() ? { ...opts, maxSteps: 1 } : opts;
    let s = 1e-12, factor = 100, steps = r.steps, iterations = r.iterations;
    const history = r.history.slice();
    try {
      this.generationScale = s;
      let q = this._solveSteady(sub);
      steps += q.steps;
      iterations += q.iterations;
      if (!q.converged) return { converged: false, steps, iterations, history };
      while (s < 1) {
        const next = Math.min(1, s * factor);
        const u1 = Float64Array.from(this.u), u1Lo = Float64Array.from(this.uLo), v1 = Float64Array.from(this.termV);
        this.generationScale = next;
        q = this._solveSteady(sub);
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
          if (factor < 1.01) return { converged: false, steps, iterations, history };
        }
      }
      return { converged: true, steps, iterations, history };
    } finally {
      this.generationScale = 1;
    }
  }

  _continuation(opts, level, target, r, warm) {
    const sub = this._directSteady() ? { ...opts, maxSteps: 1 } : opts;
    let V = level, dV = (target - level) / 8, steps = r.steps, iterations = r.iterations;
    const history = r.history.slice();
    try {
      this.sourceOverride.set(1, V);
      let q;
      if (!warm) {
        q = this._solveSteady(opts);
        steps += q.steps;
        iterations += q.iterations;
        if (!q.converged) return { ...r, steps, iterations };
      }
      while (V !== target) {
        const next = Math.abs(target - V) <= Math.abs(dV) ? target : V + dV;
        const u1 = Float64Array.from(this.u), u1Lo = Float64Array.from(this.uLo), v1 = Float64Array.from(this.termV);
        this.sourceOverride.set(1, next);
        q = this._solveSteady(sub);
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
      return { converged: true, steps, iterations, history };
    } finally {
      this.sourceOverride.delete(1);
      this._refreshSources();
    }
  }

  _solveSteady(opts = {}) {
    this.constrained = true; // spectators' amounts as constraints in the dt = ∞ solves
    this.steady = true;
    try {
      return this._steadySteps(opts);
    } finally {
      this.constrained = false;
      this.steady = false;
    }
  }

  _steadySteps({ maxSteps = 80, tol = 1e-11 } = {}) {
    const time = this.time;
    const tau = this.slowestTime();
    const direct = this._directSteady();
    const giant = direct ? Infinity : 1e6 * tau;
    if (direct) this._renormalizeSpectators(); // a starting point with the right amounts
    let dt = giant;
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
        if (dt < Math.min(tau * 1e-15, this.fastestTime() * 1e-2)) break;
        continue;
      }
      if (dt === Infinity) {
        converged = true; // the steady equations themselves were solved
        break;
      }
      if (dt >= giant && this._maxPotentialStep(this._diff()) < tol) {
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
        dt *= 10; // ramping up after a failure
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
