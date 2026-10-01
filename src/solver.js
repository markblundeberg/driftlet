// Butler–Volmer factor g(a) = e^{αa} − e^{−(1−α)a} = e^{−(1−α)a}·expm1(a) (precise near a = 0),
// and its derivative g′(a) = α e^{αa} + (1−α) e^{−(1−α)a}.
function bvFactor(a, alpha) {
  const e = Math.exp(-(1 - alpha) * a);
  return { g: e * Math.expm1(a), gp: alpha * Math.exp(alpha * a) + (1 - alpha) * e };
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
    const M = n + 1 + nRx;
    const nNodes = grid.nNodes;
    const nFaces = regions.length - 1;
    // A floating terminal with no fixed species to read its voltage from gets its own block.
    const terminalUnknown = model.circuit.mode !== 'voltage' && model.circuit.terminalUnknown;
    const nB = nNodes + nFaces + (terminalUnknown ? 1 : 0);
    this.terminalBlock = terminalUnknown ? nB - 1 : -1;
    this.n = n;
    this.M = M;
    this.nRx = nRx;
    this.nNodes = nNodes;
    this.nFaces = nFaces;
    this.nB = nB;
    this.VT = model.RT / FARADAY; // thermal voltage, V


    // Solver block of each grid node and of each interface flux node.
    this.blockOfNode = new Int32Array(nNodes);
    for (let g = 0; g < nNodes; g++) this.blockOfNode[g] = g + grid.nodeRegion[g];
    this.blockOfFace = new Int32Array(nFaces);
    for (let f = 0; f < nFaces; f++) this.blockOfFace[f] = grid.regionEnd[f] + f + 1;

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
    this.phiUndefined = new Uint8Array(nNodes);
    regions.forEach((reg, r) => {
      const mat = materials[reg.material];
      if (mat.phiFree) {
        for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) this.phiUndefined[g] = 1;
        return;
      }
      if (species.some((sp, i) => mat.present[i] && sp.z !== 0)) return;
      const pins = (ct) => ct.phi.type === 'capacitive' || ct.phi.type === 'dipole';
      const leftOpen = r === 0 ? pins(contacts.left) : model.interfaces[r - 1].phi.type !== 'neutral';
      const rightOpen =
        r === regions.length - 1 ? pins(contacts.right) : model.interfaces[r].phi.type !== 'neutral';
      if (mat.epsr === 0 || !(leftOpen || rightOpen)) {
        for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) this.phiUndefined[g] = 1;
      }
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
    // Contact bookkeeping, filled by assemble(): particle flux toward +x through each contact,
    // and the displacement there (the metal's surface charge for a neutral link).
    this.contactFlux = { left: new Float64Array(n), right: new Float64Array(n) };
    this.portFlux = model.ports.map(() => new Float64Array(n)); // into the device, mol/(m²·s)
    this.contactD = { left: 0, right: 0 };
    this.contactDOld = { left: 0, right: 0 };
    // Total current through the last segment and its derivatives (for current/load circuits).
    this.segI = 0;
    this.segIJac = new Float64Array(2 * M); // [∂/∂(block last−1 slots), ∂/∂(block last slots)]
    this.segDOld = 0;
    this._termI = 0;
    this._termIJac = new Float64Array(M + 1);
    this.contactDStart = { left: 0, right: 0 };
    // Accepted steps, most recent first: start time, size, start state (for BDF2 and the
    // error estimate of adaptive stepping).
    this.history = [];
    this.dtNext = undefined;

    // Bulk reactions as flat participant lists: reactants with +ν, products with −ν.
    this.rxs = model.reactions.map((rx) => ({
      kf: rx.kf,
      fixedA: rx.fixedA,
      sp: Int32Array.from([...rx.reactants, ...rx.products], (p) => p.i),
      nu: Float64Array.from([...rx.reactants.map((p) => p.nu), ...rx.products.map((p) => -p.nu)]),
    }));

    // The bulk reactions running in each material.
    this.rxsIn = materials.map((_, m) => this.rxs.filter((rx) => rx.kf[m] > 0));

    this._findStretches();
    this.initFromComposition();
    this.referenceAmounts = this.stretches.map((st) => this.amount(st));
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
    }
    model.interfaces.forEach((itf, f) => {
      const bf = this.blockOfFace[f];
      active[bf * M] = itf.phi.type === 'neutral' ? 0 : 1;
      for (let i = 0; i < n; i++) active[bf * M + 1 + i] = itf.links[i].type === 'blocked' ? 0 : 1;
      for (let k = 0; k < itf.reactions.length; k++) active[bf * M + 1 + n + k] = 1;
    });
    if (this.terminalBlock >= 0) active[this.terminalBlock * M] = 1;
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
        const ports = model.ports.flatMap((port, k) => (port.region >= r0 && port.region <= r && port.species[i].type !== 'blocked' ? [k] : []));
        const connected = leftOpen || rightOpen || ports.length > 0;
        this.stretches.push({
          ports,
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
    // A stretch exchanging with a fed stretch through a face reaction is fed too (e.g. ions
    // between two electrodes that each reach a contact): nothing is conserved there.
    for (let changed = true; changed; ) {
      changed = false;
      model.interfaces.forEach((itf, f) => {
        for (const rx of itf.reactions) {
          const touching = rx.part.map((p) => this.stretches[this.stretchOf[(f + p.side) * n + p.i]]);
          if (touching.some((st) => st.connected) && touching.some((st) => !st.connected)) {
            for (const st of touching) st.connected = true;
            changed = true;
          }
        }
      });
    }
    for (const st of this.stretches) {
      st.spectator = !st.connected && !st.reactive;
      // Mobile throughout: its steady state is then a single level fixed by its amount.
      st.mobile = true;
      for (let q = st.regions[0]; q <= st.regions[1]; q++) if (!(materials[regions[q].material].D[st.species] > 0)) st.mobile = false;
    }
    // Spectators whose steady state is solved directly, with the conservation of their amount
    // in place of one (redundant) balance row. Immobile ones conserve node by node instead, and
    // are left to giant time steps.
    this.constraints = [];
    this.stretches.forEach((st, k) => {
      if (!st.spectator || !st.mobile) return;
      const nNodes = st.nodes[1] - st.nodes[0] + 1;
      this.constraints.push({
        stretch: k,
        row: this.blockOfNode[st.nodes[0]] * this.M + 1 + st.species,
        idx: new Int32Array(nNodes * this.M),
        w: new Float64Array(nNodes * this.M),
        len: 0,
        res: 0,
        q: new Float64Array(this.sys.size + 1),
      });
    });
    this.constrained = false;
  }

  // Every stretch either reaches a contact or is a mobile spectator: the steady equations can
  // be solved directly (dt = ∞), with the spectators' amounts as constraints.
  _directSteady() {
    return this.stretches.every((st) => st.connected || (st.spectator && st.mobile));
  }

  // Conservation rows for the spectators (steady solves only): Σ v c_i = amount over the
  // stretch replaces the balance row of its first node, which in steady state is the negative
  // sum of the others. That row is dense, so it's kept aside: the factorised matrix gets a pin
  // (identity row) there instead, and _solveConstrained restores the constraint.
  _applyConstraints() {
    const { n, M, res, c, dA: d } = this;
    const vol = this.model.grid.vol;
    for (const cs of this.constraints) {
      const st = this.stretches[cs.stretch], i = st.species;
      let amount = 0, len = 0;
      for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
        amount += vol[g] * c[g * n + i];
        this._dc(g, i, d);
        const b = this.blockOfNode[g];
        for (let r = 0; r < M; r++) {
          if (d[r] === 0 || this.loc[b * M + r] < 0) continue;
          cs.idx[len] = this.rix[b * M + r];
          cs.w[len++] = vol[g] * d[r];
        }
      }
      cs.len = len;
      cs.res = amount - this.referenceAmounts[cs.stretch];
      const b0 = Math.floor(cs.row / M), r0 = cs.row % M;
      this._replaceRow(b0, r0);
      this._j(b0, r0, b0, r0, 1);
      res[this.rix[cs.row]] = 0;
    }
  }

  // Solve J δ = rhs where J has the constraint rows (bordered system): δ = p + Σ_k q_k s_k, with
  // p and q_k from the pinned matrix (q_k the response to a unit pin at row k, the stretch's
  // free level) and s from the k×k system W·δ = residuals.
  _solveConstrained(rhs, delta) {
    const cons = this.constraints, K = cons.length;
    this._solveLinear(rhs, delta);
    const e = this.dWork ?? (this.dWork = new Float64Array(this.sys.size + 1));
    for (const cs of cons) {
      e[this.rix[cs.row]] = 1;
      this._solveLinear(e, cs.q);
      e[this.rix[cs.row]] = 0;
    }
    const S = Array.from({ length: K }, () => new Float64Array(K + 1));
    cons.forEach((cj, j) => {
      let wp = 0, mx = 0;
      for (let a = 0; a < cj.len; a++) wp += cj.w[a] * delta[cj.idx[a]];
      for (let k = 0; k < K; k++) {
        let s = 0;
        for (let a = 0; a < cj.len; a++) s += cj.w[a] * cons[k].q[cj.idx[a]];
        S[j][k] = s;
        mx = Math.max(mx, Math.abs(s));
      }
      S[j][K] = cj.res - wp;
      if (mx > 0) for (let k = 0; k <= K; k++) S[j][k] /= mx;
    });
    // Small dense solve with partial pivoting.
    for (let col = 0; col < K; col++) {
      let p = col;
      for (let r = col + 1; r < K; r++) if (Math.abs(S[r][col]) > Math.abs(S[p][col])) p = r;
      [S[col], S[p]] = [S[p], S[col]];
      if (S[col][col] === 0) throw new SolverError('steady state: a conserved amount is not determined by its level');
      for (let r = col + 1; r < K; r++) {
        const f = S[r][col] / S[col][col];
        for (let k = col; k <= K; k++) S[r][k] -= f * S[col][k];
      }
    }
    const sol = new Float64Array(K);
    for (let r = K - 1; r >= 0; r--) {
      let v = S[r][K];
      for (let k = r + 1; k < K; k++) v -= S[r][k] * sol[k];
      sol[r] = v / S[r][r];
    }
    for (let k = 0; k < K; k++) {
      const q = cons[k].q, sk = sol[k];
      for (let a = 0; a < delta.length; a++) delta[a] += sk * q[a];
    }
  }

  /** Total amount (mol per unit area) of a stretch's species in the current state. */
  amount(stretch) {
    const { n, c, model } = this;
    const vol = model.grid.vol;
    let s = 0;
    for (let g = stretch.nodes[0]; g <= stretch.nodes[1]; g++) s += vol[g] * c[g * n + stretch.species];
    return s;
  }

  /** η_i (μ̄/RT) of a port's outside level for species i. */
  portEta(port, i) {
    const link = port.species[i];
    return this.z[i] === 0 ? link.mu / this.model.RT : (this.z[i] * (port.V + link.offset)) / this.VT;
  }

  /** η_i/RT that a fixed contact link imposes (or NaN if the link isn't fixed). */
  contactEta(side, i) {
    const ct = this.model.contacts[side], link = ct.species[i];
    if (link.type !== 'equilibrium') return NaN;
    const z = this.z[i];
    return z === 0 ? link.mu / this.model.RT : (z * (ct.V + link.offset)) / this.VT;
  }

  /**
   * Cold start. Species connected to a contact take that contact's level; spectators take
   * their region's c0 (which fixes their conserved amount). φ in each region is then chosen
   * for local neutrality, or continued across the interface dipole if nothing there responds.
   */
  initFromComposition() {
    const { model, n, M, u, z } = this;
    const { regions, species, interfaces, materials } = model;
    const grid = model.grid;
    u.fill(0);
    this.uLo.fill(0);
    const eta = new Float64Array(n), cFix = new Float64Array(n), mode = new Int8Array(n); // 1 level, 2 amount
    let phiHat = 0;
    for (let r = 0; r < regions.length; r++) {
      if (r > 0) phiHat += interfaces[r - 1].dipole / this.VT;
      const reg = regions[r], mat = materials[reg.material];
      mode.fill(0);
      for (let i = 0; i < n; i++) {
        if (!mat.present[i]) continue;
        const st = this.stretches[this.stretchOf[r * n + i]];
        if (!st.contactFed && mat.conductor) {
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
              `regions[${r}].c0.${species[i].name}: a species that doesn't reach a contact needs its initial concentration`,
            );
          }
          mode[i] = 2;
          cFix[i] = reg.c0[i];
        } else {
          const left = st.regions[0] === 0 ? this.contactEta('left', i) : NaN;
          eta[i] = Number.isFinite(left) ? left : this.contactEta('right', i);
          if (!Number.isFinite(eta[i]) && st.ports.length > 0) eta[i] = this.portEta(model.ports[st.ports[0]], i);
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
      if (!mat.phiFree && mode.some((m, i) => m === 1 && z[i] !== 0)) {
        let lo = phiHat - 1, hi = phiHat + 1;
        while (charge(lo) < 0 && lo > -1e4) lo -= 2 * (hi - lo);
        while (charge(hi) > 0 && hi < 1e4) hi += 2 * (hi - lo);
        for (let it = 0; it < 200 && hi - lo > 1e-12; it++) {
          const m = 0.5 * (lo + hi);
          if (charge(m) > 0) lo = m;
          else hi = m;
        }
        phiHat = 0.5 * (lo + hi);
      }
      if (!mat.ideal) this._materialAt(mat, reg.background, eta, mode, cFix, phiHat, zeta, cc);
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
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
    if (this.terminalBlock >= 0) u[this.terminalBlock * M] = model.contacts.right.V / this.VT;
    this.computeConcentrations();
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
        if (f >= 0) c[g * n + im] = (this.sheetSign[g] * u[this.blockOfFace[f] * M]) / (z[im] * FARADAY * this.model.grid.vol[g]);
        continue;
      }
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        c[k] = this.present[k] ? this.cRef[k] * Math.exp(u[b * M + 1 + i] + uLo[b * M + 1 + i] - this.mu0hat[k] - z[i] * phiHat) : 0;
      }
    }
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
    const R = this.rix;
    const { model, n, M, u, uLo, res, c, cOld, z, VT, sys } = this;
    const { grid, materials, regions, interfaces, contacts } = model;
    const F = FARADAY;
    sys.clear();
    res.fill(0);
    this.computeConcentrations();
    // Current through the last segment (for the circuit): displacement, then conduction.
    this.segI = 0;
    this.segD = 0;
    this.segIJac.fill(0);

    // Regions, each assembled by the kernel for its kind: a conductor (its carrier only), a
    // dilute region (ideal statistics, the fast path) or a concentrated one (any statistics).
    for (let r = 0; r < regions.length; r++) {
      const reg = regions[r], mat = materials[reg.material];
      const g0 = grid.regionStart[r], g1 = grid.regionEnd[r];
      if (mat.conductor) {
        for (let g = g0; g <= g1; g++) this._nodeConductor(g, dt);
        for (let s = g0; s < g1; s++) this._segmentConductor(s, s + r, s + r + 1, mat, grid.segLength[s], s === this.nNodes - 2);
        continue;
      }
      if (mat.ideal) this._nodesDilute(g0, g1, dt);
      else for (let g = g0; g <= g1; g++) this._nodeConcentrated(g, dt);
      const rxs = this.rxsIn[reg.material];
      if (rxs.length > 0) for (let g = g0; g <= g1; g++) this._bulkReactions(g, rxs);
      for (let s = g0; s < g1; s++) {
        const bL = s + r, bR = bL + 1, h = grid.segLength[s], lastSeg = s === this.nNodes - 2;
        this._segmentDisplacement(s, bL, bR, mat, h, lastSeg, dt);
        if (reg.mixing > 0) this._segmentMixing(s, bL, bR, reg.mixing, h, mat, lastSeg);
        if (!mat.ideal) this._segmentConcentrated(s, bL, bR, mat, h, lastSeg, reg.velocity);
      }
      if (mat.ideal) this._segmentsDilute(g0, g1, r, mat, reg.velocity);
    }

    // Interfaces: the flux node carries D and N_i; its rows are the interface laws.
    for (let f = 0; f < this.nFaces; f++) {
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
        if (law === 'dipole') {
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
        res[R[bf * M]] = u[bf * M]; // D = 0; the jump is whatever each side's neutrality needs
        this._j(bf, 0, bf, 0, 1);
      } else {
        const jump = u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) - itf.dipole / VT;
        if (law === 'dipole') {
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
      // D_f enters each side's Gauss row; a metal side holds it as surface carriers instead.
      const gL = grid.regionEnd[f], gR = gL + 1;
      if (this.nodeConductor[gL] < 0) {
        res[R[bL * M]] += u[bf * M];
        this._j(bL, 0, bf, 0, 1);
      }
      if (this.nodeConductor[gR] < 0) {
        res[R[bR * M]] -= u[bf * M] + itf.sheetCharge;
        this._j(bR, 0, bf, 0, -1);
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
        } else {
          // conductance: J = G (V_L − V_R), V = V_T η / z  ⇒  N = G V_T (η_L − η_R) / (z² F)
          const gG = (itf.links[i].G * VT) / (z[i] * z[i] * F);
          res[R[o]] = u[o] - gG * deta;
          this._j(bf, r, bf, r, 1);
          this._j(bf, r, bL, r, -gG);
          this._j(bf, r, bR, r, gG);
        }
        res[R[bL * M + r]] += u[o];
        this._j(bL, r, bf, r, 1);
        res[R[bR * M + r]] -= u[o];
        this._j(bR, r, bf, r, -1);
      }
      itf.reactions.forEach((rx, k) => this._faceReaction(rx, k, f, bf, bL, bR));
    }

    // Internal ports (after every other term at their nodes, so a held level can read its flux).
    model.ports.forEach((port, k) => this._port(port, this.portFlux[k]));

    // Contacts.
    const circuit = model.circuit;
    const tb = this.terminalBlock; // extra block carrying a floating terminal voltage, or −1
    this._termI = 0;
    this._termIJac.fill(0); // [∂/∂(last node slots), ∂/∂V_t]
    for (const side of ['left', 'right']) this._contact(side, dt);
    if (tb >= 0) {
      // Circuit law for the floating terminal: I_contact(V_t, last node) − I_circuit(V_t) = 0.
      const Vt = this.terminalV;
      let I = circuit.I, dIdV = 0;
      if (circuit.mode === 'load') {
        I = (Vt - contacts.left.V - circuit.V) / circuit.R;
        dIdV = 1 / circuit.R;
      }
      res[R[tb * M]] = this._termI - I;
      for (let r = 0; r < M; r++) this._j(tb, 0, tb - 1, r, this._termIJac[r]);
      this._j(tb, 0, tb, 0, (this._termIJac[M] - dIdV) * VT);
      for (let r = 1; r < M; r++) {
        res[R[tb * M + r]] = u[tb * M + r];
        this._j(tb, r, tb, r, 1);
      }
    }
    if (this.constrained && dt === Infinity && this.constraints.length > 0) this._applyConstraints();
  }

  // A conductor node: only its carrier's balance (the segment flux J has its own row, in
  // _segmentConductor). Its bulk stores nothing; a charged face's sheet sits in the edge node.
  _nodeConductor(g, dt) {
    const { n, M, res, c, cOld, z } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g], im = this.nodeConductor[g];
    const r = 1 + im, k = g * n + im;
    res[R[b * M + r]] += (v * (c[k] - cOld[k])) / dt;
    const f = this.sheetFace[g];
    if (f >= 0) this._j(b, r, this.blockOfFace[f], 0, this.sheetSign[g] / (z[im] * FARADAY * dt));
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
      let q = rhoFixed[g], dq = 0;
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        if (!present[k]) continue;
        const ck = c[k], l = loc[lb + 1 + i];
        res[R[lb + 1 + i]] += (v * (ck - cOld[k])) / dt;
        JB[oB + l * m + l] += (v * ck) / dt;
        q += F * z[i] * ck;
        dq += F * z[i] * z[i] * ck;
        if (p >= 0 && z[i] !== 0) {
          JB[oB + l * m + p] += (-v * z[i] * ck) / dt;
          JB[oB + p * m + l] += -v * F * z[i] * ck;
        }
      }
      if (p >= 0) {
        res[R[lb]] -= v * q;
        JB[oB + p * m + p] += v * dq;
      }
    }
  }

  // A concentrated node (any statistics): ∂c_i/∂η_j = K_ij and ∂c_i/∂φ̂ = −(Kz)_i.
  _nodeConcentrated(g, dt) {
    const { n, M, res, c, cOld, z } = this, R = this.rix, F = FARADAY;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g], Kg = g * n * n;
    let q = this.rhoFixed[g], dq = 0;
    for (let i = 0; i < n; i++) {
      const k = g * n + i, r = 1 + i;
      if (!this.present[k]) continue;
      res[R[b * M + r]] += (v * (c[k] - cOld[k])) / dt;
      for (let j = 0; j < n; j++) {
        const Kij = this.K[Kg + i * n + j];
        if (Kij !== 0) this._j(b, r, b, 1 + j, (v * Kij) / dt);
      }
      const Kz = this._Kz(g, i);
      this._j(b, r, b, 0, (-v * Kz) / dt);
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

  // Bulk reactions at node g: r = k_f Π c_R^ν · (−expm1(−a)), with a = A/RT from the
  // (compensated) η. Reactants are consumed (+v·ν·r in their balance), products made (−v·ν·r).
  _bulkReactions(g, rxs) {
    const { n, M, u, uLo, res, c } = this, R = this.rix;
    const b = this.blockOfNode[g], v = this.model.grid.vol[g], m = this.nodeMaterial[g];
    const dr = this.dA;
    for (let x = 0; x < rxs.length; x++) {
      const rx = rxs[x];
      let P = rx.kf[m], aHi = rx.fixedA, aLo = 0;
      dr.fill(0); // ∂ ln P / ∂slot
      // Participants: reactants with ν > 0, then products with ν < 0.
      for (let p = 0; p < rx.sp.length; p++) {
        const i = rx.sp[p], nu = rx.nu[p];
        if (nu > 0) {
          const ci = c[g * n + i];
          P *= nu === 1 ? ci : ci ** nu;
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
  _segmentDisplacement(s, bL, bR, mat, h, lastSeg, dt) {
    const { M, u, res, sys, loc } = this, R = this.rix;
    const k = this.phiUndefined[s] ? 0 : (mat.epsr * EPS0 * this.VT) / h; // ε = 0: no displacement
    const D = -k * (u[bR * M] - u[bL * M]);
    if (lastSeg) {
      this.segD = D;
      if (Number.isFinite(dt)) {
        this.segI += (D - this.segDOld) / dt;
        this.segIJac[0] = k / dt;
        this.segIJac[M] = -k / dt;
      }
    }
    if (k === 0) return;
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
    const segLength = this.model.grid.segLength, last = this.nNodes - 2;
    for (let s = g0; s < g1; s++) {
      const bL = s + region, bR = bL + 1, h = segLength[s], lastSeg = s === last;
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
        if (lastSeg) {
          const q = F * zi;
          this.segI += q * N;
          this.segIJac[r] += q * dNdEtaL;
          this.segIJac[M + r] += q * dNdEtaR;
          this.segIJac[0] += q * dNdPhiL;
          this.segIJac[M] += q * dNdPhiR;
        }
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

  // Scharfetter–Gummel with non-ideal statistics. The excess ex = ζ − ln(c/c_ref) acts as an
  // extra potential, linear along the segment like φ, so Δ = zΔφ̂ + Δex and
  //   N = −(D/h)·B(Δ)·c_L·expm1(η_R − η_L),
  // still exactly zero at equilibrium. c_L and ex depend on every ζ at their node through K.
  _segmentConcentrated(s, bL, bR, mat, h, lastSeg, vel) {
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
      if (lastSeg) {
        const q = FARADAY * zi;
        this.segI += q * N;
        for (let k = 0; k < M; k++) {
          this.segIJac[k] += q * jL[k];
          this.segIJac[M + k] += q * jR[k];
        }
      }
    }
  }

  // Ohmic conduction in a metal, in mixed form. The flux J (slot 0 of the left node) obeys
  //   η_R − η_L + J/g = 0,  g = σ RT/(z²F² h)   (Ohm's law, J = −g Δη),
  // and enters the two nodes' carrier balances as outflow and inflow.
  _segmentConductor(s, bL, bR, mat, h, lastSeg) {
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
    if (lastSeg) {
      const q = FARADAY * z[i];
      this.segI += q * J;
      this.segIJac[0] += q;
    }
  }

  // Eddy mixing: N = −(D_mix/RT) P ∇μ̄ with P = C − (Cz)(Cz)ᵀ/(zᵀCz), C = diag(c). It mixes
  // composition without carrying current (zᵀP = 0) and vanishes exactly at equilibrium. On a
  // segment, N_i = −(D_mix/h) Σ_j P̄_ij Δη_j, with P̄ from the logarithmic mean of each c, so a
  // neutral species gets exactly −D_mix Δc/h. Only mobile species (D > 0) take part.
  _segmentMixing(s, bL, bR, Dm, h, mat, lastSeg) {
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
      if (lastSeg) {
        const qz = FARADAY * z[i];
        this.segI += qz * N;
        for (let t = 0; t < M; t++) {
          this.segIJac[t] += qz * jL[t];
          this.segIJac[M + t] += qz * jR[t];
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

  // A face reaction (see normalizeFaceReactions). Its rate r_k is an unknown of the face block,
  // with the row r_k − rate(u_L, u_R) = 0, and each participant's edge node takes ν·r_k (made
  // there when ν > 0). The rate couples the two edge nodes only through the face block between
  // them, which keeps the system block-tridiagonal.
  _faceReaction(rx, k, f, bf, bL, bR) {
    const { n, M, u, uLo, c, res } = this;
    const R = this.rix;
    const grid = this.model.grid;
    const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
    const al = rx.alpha;
    const dL = this.dA.fill(0), dR = this.dB.fill(0); // ∂ ln(prefactor)/∂slot, per side
    let pref = rx.k0, aHi = rx.fixedA, aLo = 0;
    for (const p of rx.part) {
      const g = p.side ? gR : gL, b = p.side ? bR : bL, o = b * M + 1 + p.i;
      aHi -= p.nu * u[o];
      aLo -= p.nu * uLo[o];
      if (this.nodeConductor[g] === p.i) continue; // a conductor's carrier has activity 1
      const e = p.nu < 0 ? -p.nu * (1 - al) : p.nu * al;
      pref *= (c[g * n + p.i] / this.cRef[g * n + p.i]) ** e;
      this._dlnc(g, p.i, e, p.side ? dR : dL);
    }
    const { g: bv, gp } = bvFactor(aHi + aLo, al);
    const rate = pref * bv;
    // ∂rate/∂slot = rate·∂ln(prefactor) + pref·g′·∂a, with ∂a/∂η_i = −ν
    for (let s = 0; s < M; s++) {
      dL[s] *= rate;
      dR[s] *= rate;
    }
    for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= pref * gp * p.nu;
    const slot = 1 + n + k, o = bf * M + slot;
    res[R[o]] = u[o] - rate;
    this._j(bf, slot, bf, slot, 1);
    for (let s = 0; s < M; s++) {
      if (dL[s] !== 0) this._j(bf, slot, bL, s, -dL[s]);
      if (dR[s] !== 0) this._j(bf, slot, bR, s, -dR[s]);
    }
    for (const p of rx.part) {
      const b = p.side ? bR : bL;
      res[R[b * M + 1 + p.i]] -= p.nu * u[o];
      this._j(b, 1 + p.i, bf, slot, -p.nu);
    }
  }

  // An internal port's exchange with each node of its window, as a source per volume. A held
  // ('equilibrium') level replaces the node's balance row; the source is then that row's residual.
  // Ports come before the contacts, so a contact's flux readout includes a port's source there.
  _port(port, flux) {
    const R = this.rix;
    const { n, M, u, uLo, res, z, VT } = this;
    const F = FARADAY, vol = this.model.grid.vol;
    flux.fill(0);
    for (let i = 0; i < n; i++) {
      const link = port.species[i];
      if (link.type === 'blocked') continue;
      const r = 1 + i, target = this.portEta(port, i);
      for (const g of port.nodes) {
        const b = this.blockOfNode[g], o = b * M + r, v = vol[g];
        const deta = target - (u[o] + uLo[o]); // (μ̄_out − μ̄)/RT
        if (link.type === 'equilibrium') {
          if (g === 0 || g === this.nNodes - 1) continue; // a device end node's level is its contact's business
          flux[i] += res[R[o]];
          this._replaceRow(b, r);
          this._j(b, r, b, r, 1);
          res[R[o]] = -deta;
          continue;
        }
        // conductance: s = G V_T (η_out − η)/(z² F); exchange: s = k (η_out − η)
        const kk = link.type === 'conductance' ? (link.G * VT) / (z[i] * z[i] * F) : link.k;
        res[R[o]] -= v * kk * deta;
        this._j(b, r, b, r, v * kk);
        flux[i] += v * kk * deta;
      }
    }
  }

  // One contact: record the flux through its outer face, then add its exchange terms
  // (conductance links), then apply equilibrium links and the φ law.
  _contact(side, dt) {
    const R = this.rix;
    const { model, n, M, u, uLo, res, c, z, VT } = this;
    const F = FARADAY;
    const contacts = model.contacts, circuit = model.circuit;
    const ct = contacts[side];
    const g = side === 'left' ? 0 : this.nNodes - 1;
    const b = this.blockOfNode[g];
    const sgn = side === 'left' ? 1 : -1;
    const flux = this.contactFlux[side];
    const floating = side === 'right' && circuit.mode !== 'voltage';
    const readout = floating && !circuit.terminalUnknown;
    const tb = this.terminalBlock;
    const t = ct.terminal;

    // Terminal voltage, and how it depends on the unknowns.
    let Vt = ct.V;
    if (readout) {
      const ot = b * M + 1 + t;
      Vt = (VT * (u[ot] + uLo[ot])) / z[t] - ct.species[t].offset;
    } else if (floating) {
      Vt = VT * u[tb * M];
    }
    if (floating) this.terminalV = Vt;
    const addVt = (rs, val) => {
      if (readout) this._j(b, rs, b, 1 + t, (val * VT) / z[t]);
      else if (floating) this._j(b, rs, tb, 0, val * VT);
    };
    const toTerminal = floating && !readout; // accumulate the terminal-current row
    const termJ = this._termIJac;

    // Before anything is added, each balance residual is the flux through this face.
    for (let i = 0; i < n; i++) flux[i] = this.loc[b * M + 1 + i] < 0 ? 0 : sgn * res[R[b * M + 1 + i]];

    // Conductance links: J (toward the device) = G (V_out − V_i), V_out = V_t + offset; exchange
    // links (neutral species): N_in = k (μ_out − μ)/RT.
    for (let i = 0; i < n; i++) {
      const link = ct.species[i];
      if (link.type === 'exchange') {
        const o = b * M + 1 + i;
        res[R[o]] -= link.k * (link.mu / this.model.RT - (u[o] + uLo[o]));
        this._j(b, 1 + i, b, 1 + i, link.k);
        continue;
      }
      if (link.type !== 'conductance') continue;
      const o = b * M + 1 + i;
      const Vi = (VT * (u[o] + uLo[o])) / z[i];
      const Nin = (link.G * (Vt + link.offset - Vi)) / (z[i] * F); // particles entering
      const k = (link.G * VT) / (z[i] * z[i] * F);
      res[R[o]] -= Nin;
      this._j(b, 1 + i, b, 1 + i, k);
      addVt(1 + i, -link.G / (z[i] * F));
      if (toTerminal) {
        // Current toward +x leaving through the right face: −z F N_in.
        this._termI += -z[i] * F * Nin;
        termJ[1 + i] += link.G * VT / z[i];
        termJ[M] += -link.G;
      }
    }

    // Equilibrium links: Dirichlet on the known outside level.
    for (let i = 0; i < n; i++) {
      if (ct.species[i].type !== 'equilibrium') continue;
      const o = b * M + 1 + i;
      this._replaceRow(b, 1 + i);
      if (readout && i === t) {
        // I_segment − I_circuit(V_t) = 0. The last segment's total current equals the
        // terminal current exactly (box balance plus Poisson, differenced in time).
        const I = circuit.mode === 'current' ? circuit.I : (Vt - contacts.left.V - circuit.V) / circuit.R;
        res[R[o]] = this.segI - I;
        for (let r = 0; r < M; r++) {
          this._j(b, 1 + i, b - 1, r, this.segIJac[r]);
          this._j(b, 1 + i, b, r, this.segIJac[M + r]);
        }
        if (circuit.mode === 'load') this._j(b, 1 + i, b, 1 + t, -VT / (z[t] * circuit.R));
      } else if (floating && ct.species[i].mu === undefined) {
        // Charged species tied to the floating terminal: η_i = z_i (V_t + offset_i)/V_T.
        res[R[o]] = u[o] + uLo[o] - (z[i] * (Vt + ct.species[i].offset)) / VT;
        this._j(b, 1 + i, b, 1 + i, 1);
        addVt(1 + i, -z[i] / VT);
      } else {
        this._j(b, 1 + i, b, 1 + i, 1);
        res[R[o]] = u[o] - this.contactEta(side, i) + uLo[o];
      }
    }

    // φ link (none where φ is undefined: nothing at the end node responds to it).
    const link = ct.phi;
    if (this.phiUndefined[g]) {
      this.contactD[side] = 0;
    } else if (link.type === 'capacitive') {
      // Gate, or metal across a Stern layer, at φ_g = V_t − zeroCharge. D toward +x.
      const phiG = Vt - link.zeroCharge;
      const D = link.C * (phiG - VT * u[b * M]) * sgn;
      res[R[b * M]] -= sgn * D;
      this._j(b, 0, b, 0, link.C * VT);
      addVt(0, -sgn * sgn * link.C); // ∂(−sgn·D)/∂V_t = −C
      this.contactD[side] = D;
      if (toTerminal && Number.isFinite(dt)) {
        this._termI += (D - this.contactDStart.right) / dt;
        termJ[0] += (link.C * VT) / dt;
        termJ[M] += -link.C / dt;
      }
    } else if (link.type === 'dipole') {
      // Pinned: φ_edge = V_t − zeroCharge (the C → ∞ limit). The residual is the metal's charge.
      this.contactD[side] = sgn * res[R[b * M]];
      this._replaceRow(b, 0);
      this._j(b, 0, b, 0, 1);
      res[R[b * M]] = u[b * M] + uLo[b * M] - (Vt - link.zeroCharge) / VT;
      addVt(0, -1 / VT);
    } else if (link.type === 'bulk') {
      this.contactD[side] = sgn * res[R[b * M]];
      this._replaceRow(b, 0);
      let q = this.rhoFixed[g];
      let dq = 0;
      for (let i = 0; i < n; i++) {
        const k = g * n + i;
        if (!this.present[k]) continue;
        const Kz = this._Kz(g, i);
        q += F * z[i] * c[k];
        dq -= F * z[i] * Kz;
        if (Kz !== 0) this._j(b, 0, b, 1 + i, F * Kz);
      }
      this._j(b, 0, b, 0, dq);
      res[R[b * M]] = q;
    } else {
      this.contactD[side] = 0;
    }
  }

  // Displacement toward +x through the last segment, from the current state.
  _lastSegmentD() {
    const { model, M, u, VT } = this;
    const grid = model.grid, s = this.nNodes - 2;
    const mat = model.materials[model.regions[grid.segRegion[s]].material];
    const b = this.blockOfNode[s];
    return (-(mat.epsr * EPS0 * VT) / grid.segLength[s]) * (u[(b + 1) * M] - u[b * M]);
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
    const { nB, sys, res } = this;
    const { A, B, C, sizes, offA, offB, offC, offX } = sys;
    for (let b = 0; b < nB; b++) {
      const m = sizes[b], mp = b > 0 ? sizes[b - 1] : 0, mn = b < nB - 1 ? sizes[b + 1] : 0;
      for (let r = 0; r < m; r++) {
        const oa = offA[b] + r * mp, ob = offB[b] + r * m, oc = offC[b] + r * mn;
        let mx = 0;
        for (let k = 0; k < mp; k++) mx = Math.max(mx, Math.abs(A[oa + k]));
        for (let k = 0; k < m; k++) mx = Math.max(mx, Math.abs(B[ob + k]));
        for (let k = 0; k < mn; k++) mx = Math.max(mx, Math.abs(C[oc + k]));
        if (mx === 0) continue;
        const s = 1 / mx;
        for (let k = 0; k < mp; k++) A[oa + k] *= s;
        for (let k = 0; k < m; k++) B[ob + k] *= s;
        for (let k = 0; k < mn; k++) C[oc + k] *= s;
        res[offX[b] + r] *= s;
      }
    }
  }

  // Largest update among the potential-like unknowns (φ̂ and η at grid nodes).
  _maxPotentialStep(delta) {
    const { n, M } = this;
    const R = this.rix; // (the sink entry of delta is 0)
    let mx = this.terminalBlock >= 0 ? Math.abs(delta[R[this.terminalBlock * M]]) : 0;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      if (this.nodeConductor[g] < 0) mx = Math.max(mx, Math.abs(delta[R[b * M]])); // (a metal's slot 0 is a flux)
      for (let i = 0; i < n; i++) {
        if (this.present[g * n + i]) mx = Math.max(mx, Math.abs(delta[R[b * M + 1 + i]]));
      }
    }
    return mx;
  }

  /**
   * Newton iteration for one backward-Euler step from the current cOld.
   * @returns {{converged: boolean, iterations: number, history: number[]}}
   */
  newton(dt, { maxIter = 60, tol = 1e-10, maxStep = 10 } = {}) {
    const { u, delta, res } = this;
    const history = [];
    for (let it = 1; it <= maxIter; it++) {
      try {
        this.assemble(dt);
      } catch (err) {
        // A statistics model can fail far from the solution (e.g. Debye–Hückel beyond its range).
        if (!(err instanceof SolverError)) throw err;
        return { converged: false, iterations: it, history, error: err.message };
      }
      this._equilibrate();
      let rmax = 0;
      for (let k = 0; k < this.SINK; k++) rmax = Math.max(rmax, Math.abs(res[k]));
      try {
        this._factor();
      } catch (err) {
        return { converged: false, iterations: it, history, error: err.message };
      }
      try {
        if (this.constrained && dt === Infinity && this.constraints.length > 0) this._solveConstrained(res, delta);
        else this._solveLinear(res, delta);
      } catch (err) {
        if (!(err instanceof SolverError)) throw err;
        return { converged: false, iterations: it, history, error: err.message };
      }
      const step = this._maxPotentialStep(delta);
      history.push(step);
      if (!Number.isFinite(step)) return { converged: false, iterations: it, history, error: 'non-finite update' };
      // Give up early on clear divergence; the caller will take a smaller step instead.
      if (step > 1e4 || (it > 6 && step > 10 * history[0])) {
        this.computeConcentrations();
        return { converged: false, iterations: it, history, error: 'diverging' };
      }
      const alpha = step > maxStep ? maxStep / step : 1;
      this._addToState(delta, -alpha);
      if (alpha === 1 && step < tol) {
        this.computeConcentrations();
        return { converged: true, iterations: it, history, residual: rmax };
      }
    }
    this.computeConcentrations();
    return { converged: false, iterations: maxIter, history };
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
    const prev = this.history[0];
    const w = prev ? dt / prev.dt : 0;
    const bdf = opts.method === 'bdf2' && prev !== undefined && w <= 2;
    this.uPrev.set(this.u);
    this.uPrevLo.set(this.uLo);
    this.computeConcentrations();
    const cN = Float64Array.from(this.c);
    const segDN = this._lastSegmentD();
    // Contact displacement before the step: as it was at the end of the previous step (under
    // the parameters then), so a gate-voltage change shows up as displacement current.
    if (!this.contactDEnd) {
      this.assemble(dt);
      this.contactDEnd = { ...this.contactD };
    }
    const DN = { ...this.contactDEnd };
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
      this.segDOld = b1 * segDN - b2 * prev.segD;
      this.contactDStart = { left: b1 * DN.left - b2 * prev.D.left, right: b1 * DN.right - b2 * prev.D.right };
      dtEff = dt / a0;
    } else {
      this.cOld.set(cN);
      this.segDOld = segDN;
      this.contactDStart = DN;
    }
    const result = this.newton(dtEff, opts);
    result.bdf = bdf;
    if (result.converged) {
      this.history.unshift({ t: this.time, dt, u: Float64Array.from(this.uPrev), c: cN, segD: segDN, D: DN });
      if (this.history.length > 3) this.history.length = 3;
      this.time += dt;
      this.lastDt = dtEff;
      this.contactDOld = this.contactDStart;
      this._accumulateBoundaryIntake(dtEff, cN);
    } else {
      this.u.set(this.uPrev);
      this.uLo.set(this.uPrevLo);
      this.computeConcentrations();
    }
    return result;
  }

  /**
   * Adaptive time stepping to tEnd: variable-step BDF2 (or backward Euler), with the local error
   * estimated against an explicit predictor through the previous states, and controlled to
   * `tol` thermal units per step in every potential (φ̂, each η, a floating terminal). The first
   * step is checked by step doubling. Stops early when `budgetMs` of wall time is used, so an
   * animation can call it once per frame; the step size carries over between calls.
   */
  integrate(tEnd, opts = {}) {
    const { tol = 1e-3, dtMax = Infinity, budgetMs = Infinity, maxSteps = 100000, method = 'bdf2' } = opts;
    const clock = () => (globalThis.performance ? globalThis.performance.now() : Date.now());
    const start = clock();
    const trace = { t: [], current: [], voltage: [] };
    let steps = 0, rejected = 0, iterations = 0, failed = false;
    let dt = this.dtNext ?? opts.dt0 ?? (tEnd - this.time) * 1e-4;
    const pred = new Float64Array(this.u.length), guess = new Float64Array(this.u.length);
    const factor = (err, p) => (err > 0 ? Math.min(2, Math.max(0.2, 0.9 * (tol / err) ** (1 / (p + 1)))) : 2);
    while (this.time < tEnd) {
      if (steps + rejected >= maxSteps || clock() - start > budgetMs) break;
      const remaining = tEnd - this.time;
      let h = Math.min(dt, dtMax);
      const clamped = h >= remaining * (1 - 1e-9);
      if (clamped) h = remaining;
      const snap = this._snapshot();
      let err, p;
      if (this.history.length === 0) {
        // First step: backward Euler, checked against two half steps (whose result is kept).
        const full = this.step(h);
        iterations += full.iterations;
        if (full.converged) {
          const uFull = Float64Array.from(this.u);
          this._restore(snap);
          const a = this.step(h / 2), b = a.converged ? this.step(h / 2) : a;
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
          if (dt < 1e-14 * Math.max(tEnd, 1e-300)) { failed = true; break; }
          continue;
        }
      } else {
        const r = this.step(h, { method, guess: this._extrapolate(guess, this.time + h) });
        iterations += r.iterations;
        if (!r.converged) {
          rejected++;
          dt = h / 4;
          if (dt < 1e-14 * Math.max(tEnd, 1e-300)) { failed = true; break; }
          continue;
        }
        p = this._predict(pred, r.bdf);
        err = this._errorNorm(pred) * p.scale;
        p = p.order;
      }
      if (err > tol) {
        this._restore(snap);
        rejected++;
        dt = h * factor(err, p);
        continue;
      }
      steps++;
      dt = clamped ? Math.max(dt, h * factor(err, p)) : h * factor(err, p);
      trace.t.push(this.time);
      trace.current.push(this._terminalCurrent());
      trace.voltage.push(this.model.circuit.mode !== 'voltage' ? this.terminalV - this.model.contacts.left.V : this.model.contacts.right.V - this.model.contacts.left.V);
    }
    this.dtNext = dt;
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
      const Cc = (h ** 3 * (1 + w) ** 2) / (w * (1 + 2 * w));
      const Cp = h * (h + hp) * (h + hp + hpp);
      return { scale: Cc / (Cc + Cp), order: 2 };
    }
    const l0 = (t - e1.t) / (e0.t - e1.t), l1 = 1 - l0;
    for (let k = 0; k < pred.length; k++) pred[k] = l0 * e0.u[k] + l1 * e1.u[k];
    // Backward Euler's error against a linear predictor; after a BDF2 step without enough
    // history this overestimates, which is safe.
    return { scale: bdf ? 1 : h / (2 * h + hp), order: 1 };
  }

  // Largest difference from ref among the potentials that carry state, in thermal units.
  _errorNorm(ref) {
    const { n, M, u } = this;
    let mx = this.terminalBlock >= 0 ? Math.abs(u[this.terminalBlock * M] - ref[this.terminalBlock * M]) : 0;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      if (!this.phiUndefined[g]) mx = Math.max(mx, Math.abs(u[b * M] - ref[b * M]));
      for (let i = 0; i < n; i++) {
        if (this.present[g * n + i]) mx = Math.max(mx, Math.abs(u[b * M + 1 + i] - ref[b * M + 1 + i]));
      }
    }
    return mx;
  }

  _snapshot() {
    return {
      u: Float64Array.from(this.u),
      uLo: Float64Array.from(this.uLo),
      cOld: Float64Array.from(this.cOld),
      time: this.time,
      lastDt: this.lastDt,
      segDOld: this.segDOld,
      contactDStart: this.contactDStart,
      contactDOld: this.contactDOld,
      contactDEnd: this.contactDEnd,
      boundaryIntake: Float64Array.from(this.boundaryIntake),
      history: this.history.slice(),
    };
  }

  _restore(s) {
    this.u.set(s.u);
    this.uLo.set(s.uLo);
    this.cOld.set(s.cOld);
    this.time = s.time;
    this.lastDt = s.lastDt;
    this.segDOld = s.segDOld;
    this.contactDStart = s.contactDStart;
    this.contactDOld = s.contactDOld;
    this.contactDEnd = s.contactDEnd;
    this.boundaryIntake.set(s.boundaryIntake);
    this.history = s.history.slice();
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
   * Small-signal impedance about the current (steady) state. Linearising the balance equations
   * gives (J + iωM)·δx = −b·δs, with J the steady Jacobian, M the storage Jacobian (from
   * J(dt) = J + M/dt) and b the residual's derivative in the source δs: the right terminal's
   * voltage (voltage mode) or the circuit current (current mode). The terminal current response
   * is read from the last segment, conduction plus iω times its displacement.
   * Z = −δV/δI: the impedance seen at the terminals (positive real part for a passive device,
   * since current toward +x leaves through the right terminal).
   * @param {ArrayLike<number>} frequencies Hz
   */
  impedance(frequencies, { profiles = false } = {}) {
    const { M, nB, sys, model, n, VT } = this;
    const circuit = model.circuit;
    if (circuit.mode === 'load') {
      throw new SolverError("impedance: use circuit mode 'voltage' or 'current' (a load resistor belongs to the external circuit)");
    }
    const N = nB * M;
    this.computeConcentrations();
    this.cOld.set(this.c);
    this.segDOld = this._lastSegmentD();
    this.assemble(Infinity);
    this.contactDStart = { ...this.contactD };
    this.assemble(Infinity);
    const J = { A: sys.A.slice(), B: sys.B.slice(), C: sys.C.slice() };
    const segJ = Float64Array.from(this.segIJac);
    const dts = 1e-30;
    this.assemble(dts);
    const S = {};
    for (const X of ['A', 'B', 'C']) S[X] = sys[X].map((v, k) => (v - J[X][k]) * dts);
    // Source derivative by central differences (exact where the residual is linear in it).
    const ct = model.contacts.right;
    const src = circuit.mode === 'voltage' ? { obj: ct, key: 'V', d: 1e-6 } : { obj: circuit, key: 'I', d: 1 };
    const s0 = src.obj[src.key], res = [];
    for (const sign of [1, -1]) {
      src.obj[src.key] = s0 + sign * src.d;
      this.assemble(Infinity);
      res.push(Float64Array.from(this.res));
    }
    src.obj[src.key] = s0;
    this.assemble(Infinity);
    const b = res[0].map((v, k) => (v - res[1][k]) / (2 * src.d));

    // The complex system has the Jacobian's layout (only the unknowns); slots that aren't
    // unknowns stay at δx = 0 in the full-width profiles.
    const { sizes, offA, offB, offC, offX } = sys;
    const csys = new ComplexBlockTridiagonal(nB, sizes);
    const Nc = csys.size;
    const rr = new Float64Array(Nc), ri = new Float64Array(Nc), cr = new Float64Array(Nc), ci = new Float64Array(Nc);
    const xr = new Float64Array(N), xi = new Float64Array(N);
    const gL = this.nNodes - 2, bL = this.blockOfNode[gL], bR = bL + 1;
    const mat = model.materials[model.regions[model.grid.segRegion[gL]].material];
    const kSeg = this.phiUndefined[gL] ? 0 : (mat.epsr * EPS0 * VT) / model.grid.segLength[gL];
    const out = { f: Float64Array.from(frequencies), Z: { re: new Float64Array(frequencies.length), im: new Float64Array(frequencies.length) } };
    if (profiles) out.profiles = [];
    Array.from(frequencies).forEach((f, q) => {
      const w = 2 * Math.PI * f;
      if (!(w > 0)) throw new SolverError('impedance: frequencies must be positive');
      // Rows scaled by their largest entry.
      for (let blk = 0; blk < nB; blk++) {
        const m = sizes[blk], mp = blk > 0 ? sizes[blk - 1] : 0, mn = blk < nB - 1 ? sizes[blk + 1] : 0;
        for (let r = 0; r < m; r++) {
          const rows = [['A', offA[blk] + r * mp, mp], ['B', offB[blk] + r * m, m], ['C', offC[blk] + r * mn, mn]];
          let mx = 0;
          for (const [X, o, len] of rows) for (let c = 0; c < len; c++) mx = Math.max(mx, Math.abs(J[X][o + c]), w * Math.abs(S[X][o + c]));
          const sc = mx > 0 ? 1 / mx : 1;
          for (const [X, o, len] of rows) {
            const JX = J[X], SX = S[X], Xr = csys[X + 'r'], Xi = csys[X + 'i'];
            for (let c = 0; c < len; c++) {
              Xr[o + c] = JX[o + c] * sc;
              Xi[o + c] = w * SX[o + c] * sc;
            }
          }
          rr[offX[blk] + r] = -b[offX[blk] + r] * sc;
          ri[offX[blk] + r] = 0;
        }
      }
      csys.factor();
      csys.solve(rr, ri, cr, ci);
      for (let k = 0; k < Nc; k++) {
        xr[this.fullOf[k]] = cr[k];
        xi[this.fullOf[k]] = ci[k];
      }
      let Zr, Zi;
      if (circuit.mode === 'voltage') {
        // δI = (∂I_cond/∂x)·δx + iω δD_segment, per volt
        let Ir = 0, Ii = 0;
        for (let r = 0; r < M; r++) {
          Ir += segJ[r] * xr[bL * M + r] + segJ[M + r] * xr[bR * M + r];
          Ii += segJ[r] * xi[bL * M + r] + segJ[M + r] * xi[bR * M + r];
        }
        const dDr = -kSeg * (xr[bR * M] - xr[bL * M]), dDi = -kSeg * (xi[bR * M] - xi[bL * M]);
        Ir -= w * dDi;
        Ii += w * dDr;
        const d2 = Ir * Ir + Ii * Ii; // Z = −1/δI
        Zr = -Ir / d2;
        Zi = Ii / d2;
      } else {
        // δV per unit current
        const tb = this.terminalBlock;
        const t = model.contacts.right.terminal;
        const [vr, vi] = tb >= 0 ? [VT * xr[tb * M], VT * xi[tb * M]] : [(VT * xr[bR * M + 1 + t]) / this.z[t], (VT * xi[bR * M + 1 + t]) / this.z[t]];
        Zr = -vr;
        Zi = -vi;
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

  // Add this step's contact fluxes (at the converged state) to each stretch's intake.
  // (A BDF2 step also moves each amount by Σ v (c* − c_n), its history term.)
  _accumulateBoundaryIntake(dt, cN) {
    this.assemble(dt);
    this.contactDEnd = { ...this.contactD };
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
  solveSteady(opts = {}) {
    const ct = this.model.contacts.right, target = ct.V, level = this.model.contacts.left.V;
    const canContinue = opts.continuation !== false && this.model.circuit.mode === 'voltage' && target !== level;
    const direct = this._directSteady();
    const u0 = Float64Array.from(this.u), u0Lo = Float64Array.from(this.uLo);
    // Where a direct solve applies and continuation is possible, don't spend long on the
    // pseudo-transient ramp: one direct attempt first.
    let r = this._solveSteady(canContinue && direct ? { ...opts, maxSteps: 1 } : opts);
    if (r.converged) this.solvedV = target;
    if (r.converged || !canContinue) return r;
    // Source continuation: solve with both terminals level (consistent with a cold start), then
    // ramp the right terminal's voltage to its target in adaptive steps.
    const restart = () => {
      this.u.set(u0);
      this.uLo.set(u0Lo);
      this.computeConcentrations();
    };
    restart();
    // Ramp from the last converged voltage when the state is that solution (a warm start),
    // otherwise from level terminals.
    const from = this.solvedV !== undefined ? this.solvedV : level;
    const c = this._continuation(opts, ct, from, target, r, from !== level);
    if (c.converged) {
      this.solvedV = target;
      return c;
    }
    restart();
    r = this._solveSteady(opts); // the full pseudo-transient ramp
    if (r.converged) this.solvedV = target;
    return { ...r, steps: r.steps + c.steps, iterations: r.iterations + c.iterations };
  }

  _continuation(opts, ct, level, target, r, warm) {
    const sub = this._directSteady() ? { ...opts, maxSteps: 1 } : opts;
    let V = level, dV = (target - level) / 8, steps = r.steps, iterations = r.iterations;
    const history = r.history.slice();
    try {
      ct.V = V;
      let q;
      if (!warm) {
        q = this._solveSteady(opts);
        steps += q.steps;
        iterations += q.iterations;
        if (!q.converged) return { ...r, steps, iterations };
      }
      while (V !== target) {
        const next = Math.abs(target - V) <= Math.abs(dV) ? target : V + dV;
        const u1 = Float64Array.from(this.u), u1Lo = Float64Array.from(this.uLo);
        ct.V = next;
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
          this.computeConcentrations();
          dV /= 4;
          if (Math.abs(dV) < 1e-6 * Math.abs(target - level)) return { converged: false, steps, iterations, history };
        }
      }
      return { converged: true, steps, iterations, history };
    } finally {
      ct.V = target;
    }
  }

  _solveSteady(opts = {}) {
    this.constrained = true; // spectators' amounts as constraints in the dt = ∞ solves
    try {
      return this._steadySteps(opts);
    } finally {
      this.constrained = false;
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
        dt = dt >= giant ? tau * 1e-6 : dt / 4;
        if (dt < tau * 1e-15) break;
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
      if (!st.spectator) return;
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
