// Butler–Volmer factor g(a) = e^{αa} − e^{−(1−α)a} = e^{−(1−α)a}·expm1(a) (precise near a = 0),
// and its derivative g′(a) = α e^{αa} + (1−α) e^{−(1−α)a}.
function bvFactor(a, alpha) {
  const e = Math.exp(-(1 - alpha) * a);
  return { g: e * Math.expm1(a), gp: alpha * Math.exp(alpha * a) + (1 - alpha) * e };
}

// Discretisation and nonlinear solver.
//
// Unknowns, per solver block (block size M = 1 + nSpecies):
//   grid node: [φ̂, η_1 … η_n], with φ̂ = Fφ/RT and η_i = μ̄_i/RT
//   flux node: [D, N_1 … N_n], the displacement and particle fluxes through an interface
// Every region boundary is a doubled grid node (one per side) with a zero-volume flux node
// between them in the linear system, so the Jacobian stays block-tridiagonal.
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
    const M = n + 1;
    const nNodes = grid.nNodes;
    const nFaces = regions.length - 1;
    // A floating terminal with no fixed species to read its voltage from gets its own block.
    const terminalUnknown = model.circuit.mode !== 'voltage' && model.circuit.terminalUnknown;
    const nB = nNodes + nFaces + (terminalUnknown ? 1 : 0);
    this.terminalBlock = terminalUnknown ? nB - 1 : -1;
    this.n = n;
    this.M = M;
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

    this.sys = new BlockTridiagonal(nB, M);
    // Unknowns as compensated double-doubles, u + uLo. Only differences of η need the extra
    // precision: a majority carrier carrying a small current has a quasi-Fermi step between
    // nodes far below the ulp of η itself (e.g. 1e-19 vs 3.5e-15), and would otherwise carry
    // exactly zero current there. The Jacobian and linear solve stay in plain doubles.
    this.u = new Float64Array(nB * M);
    this.uLo = new Float64Array(nB * M);
    this.res = new Float64Array(nB * M);
    this.delta = new Float64Array(nB * M);
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
    this.nodeMetal = new Int32Array(nNodes).fill(-1);
    this.sheetFace = new Int32Array(nNodes).fill(-1);
    this.sheetSign = new Int8Array(nNodes);
    // A metal node's φ slot (φ is undefined there) carries instead the carrier flux J through the
    // segment to its right: Ohm's law and continuity in mixed form stay well conditioned however
    // large σ is (eliminating a stiff ohmic chain directly would cancel catastrophically).
    this.metalLast = new Uint8Array(nNodes);
    for (let g = 0; g < nNodes; g++) {
      const mat = materials[regions[grid.nodeRegion[g]].material];
      if (mat.metal) {
        this.nodeMetal[g] = mat.metal.i;
        this.nodeIdeal[g] = 1;
        this.metalLast[g] = g === grid.regionEnd[grid.nodeRegion[g]] ? 1 : 0;
      }
    }
    model.interfaces.forEach((itf, f) => {
      if (!itf.metal || itf.metal.side === 'both' || itf.phi.type === 'neutral') return;
      const g = itf.metal.side === 'left' ? grid.regionEnd[f] : grid.regionStart[f + 1];
      this.sheetFace[g] = f;
      this.sheetSign[g] = itf.metal.side === 'left' ? 1 : -1; // σ_metal = D_f on the left, −D_f on the right
    });
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

    this._findStretches();
    this.initFromComposition();
    this.referenceAmounts = this.stretches.map((st) => this.amount(st));
    // ∫ (flux in − flux out) dt through the contacts, per stretch, since the reference.
    this.boundaryIntake = new Float64Array(this.stretches.length);
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
        const touches = (ct) => ct.species[i].type !== 'blocked' || ct.reactions.some((rx) => [...rx.reactants, ...rx.products].some((p) => p.i === i));
        const leftOpen = r0 === 0 && touches(contacts.left);
        const rightOpen = r === last && touches(contacts.right);
        let reactive = false;
        for (const f of [r0 - 1, r]) {
          const itf = model.interfaces[f];
          if (!itf || itf.electrode.length === 0) continue;
          if (itf.metal.i === i || itf.electrode.some((rx) => [...rx.reactants, ...rx.products].some((p) => p.i === i))) reactive = true;
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
    // A stretch exchanging with a fed stretch through an electrode reaction at a metal face is
    // fed too (e.g. ions between two metals that each reach a contact): nothing is conserved there.
    for (let changed = true; changed; ) {
      changed = false;
      model.interfaces.forEach((itf, f) => {
        if (itf.electrode.length === 0) return;
        const involved = new Set([itf.metal.i]);
        for (const rx of itf.electrode) for (const p of [...rx.reactants, ...rx.products]) involved.add(p.i);
        const touching = this.stretches.filter((st) => involved.has(st.species) && (st.regions[1] === f || st.regions[0] === f + 1));
        if (touching.some((st) => st.connected) && touching.some((st) => !st.connected)) {
          for (const st of touching) st.connected = true;
          changed = true;
        }
      });
    }
    for (const st of this.stretches) st.spectator = !st.connected && !st.reactive;
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
        if (!st.contactFed && mat.metal) {
          // A metal away from the contacts: start uncharged, with its Fermi level in equilibrium
          // with the electrode reaction on its left face if it has one, else at the running φ.
          eta[i] = z[i] * phiHat;
          const rx = r > 0 ? interfaces[r - 1].electrode[0] : undefined;
          if (rx) {
            const bn = this.blockOfNode[grid.regionEnd[r - 1]];
            let a = rx.fixedA;
            for (const p of rx.reactants) a += p.nu * u[bn * M + 1 + p.i];
            for (const p of rx.products) a -= p.nu * u[bn * M + 1 + p.i];
            eta[i] = -a / rx.electrons; // A = 0
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
        u[b * M] = mat.metal ? 0 : phiHat; // a metal's slot 0 is its segment flux
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
      const im = this.nodeMetal[g];
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

  // Add ∂(row rb, slot rs)/∂(block cb, slot cs) to the Jacobian.
  _j(rb, rs, cb, cs, v) {
    const M = this.M, o = rb * M * M + rs * M + cs;
    if (cb === rb) this.sys.B[o] += v;
    else if (cb === rb - 1) this.sys.A[o] += v;
    else if (cb === rb + 1) this.sys.C[o] += v;
    else throw new Error(`internal: non-tridiagonal coupling ${rb}→${cb}`);
  }

  /** Assemble residual and Jacobian for a backward-Euler step of size dt. */
  assemble(dt) {
    const { model, n, M, u, uLo, res, c, cOld, z, VT, sys } = this;
    const { grid, materials, regions, interfaces, contacts } = model;
    const F = FARADAY;
    sys.clear();
    res.fill(0);
    this.computeConcentrations();

    // Node terms: storage, space charge, identity rows for absent species.
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g], v = grid.vol[g];
      let q = this.rhoFixed[g], dq = 0;
      const im = this.nodeMetal[g];
      if (im >= 0) {
        for (let i = 0; i < n; i++) {
          const r = 1 + i;
          if (i !== im) {
            this._j(b, r, b, r, 1);
            continue;
          }
          const k = g * n + i;
          res[b * M + r] += (v * (c[k] - cOld[k])) / dt;
          const f = this.sheetFace[g];
          if (f >= 0) this._j(b, r, this.blockOfFace[f], 0, this.sheetSign[g] / (z[im] * F * dt));
        }
        if (this.metalLast[g]) this._j(b, 0, b, 0, 1); // no segment to its right: the J slot is idle
        continue;
      }
      if (this.nodeIdeal[g]) {
        for (let i = 0; i < n; i++) {
          const k = g * n + i, r = 1 + i;
          if (!this.present[k]) {
            this._j(b, r, b, r, 1);
            continue;
          }
          const ck = c[k];
          res[b * M + r] += (v * (ck - cOld[k])) / dt;
          this._j(b, r, b, r, (v * ck) / dt);
          this._j(b, r, b, 0, (-v * z[i] * ck) / dt);
          q += F * z[i] * ck;
          dq += F * z[i] * z[i] * ck;
          if (z[i] !== 0) this._j(b, 0, b, r, -v * F * z[i] * ck);
        }
      } else {
        // ∂c_i/∂η_j = K_ij and ∂c_i/∂φ̂ = −(Kz)_i.
        const Kg = g * n * n;
        for (let i = 0; i < n; i++) {
          const k = g * n + i, r = 1 + i;
          if (!this.present[k]) {
            this._j(b, r, b, r, 1);
            continue;
          }
          res[b * M + r] += (v * (c[k] - cOld[k])) / dt;
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
      }
      if (this.phiUndefined[g]) {
        this._j(b, 0, b, 0, 1); // no charge responds and no field reaches: φ is not defined here
      } else {
        res[b * M] -= v * q;
        this._j(b, 0, b, 0, v * dq);
      }

      // Bulk reactions: r = k_f Π c_R^ν · (−expm1(−a)), with a = A/RT from the (compensated) η.
      const m = this.nodeMaterial[g];
      const dr = this.dA, Bm = sys.B, bo = b * M * M;
      for (const rx of this.rxs) {
        const kf = rx.kf[m];
        if (!(kf > 0)) continue;
        let P = kf, aHi = rx.fixedA, aLo = 0;
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
        // Reactants are consumed (+v·ν·r in their balance), products made (−v·ν·r).
        for (let p = 0; p < rx.sp.length; p++) {
          const row = 1 + rx.sp[p], w = v * rx.nu[p], o = bo + row * M;
          res[b * M + row] += w * rate;
          for (let k = 0; k < M; k++) if (dr[k] !== 0) Bm[o + k] += w * dr[k];
        }
      }
    }

    // Ordinary segments: displacement and Scharfetter–Gummel fluxes.
    for (let s = 0; s < this.nNodes - 1; s++) {
      const reg = grid.segRegion[s];
      if (reg < 0) continue;
      const mat = materials[regions[reg].material];
      const h = grid.segLength[s];
      const bL = this.blockOfNode[s], bR = bL + 1;
      const phiL = u[bL * M], phiR = u[bR * M];

      const k = this.phiUndefined[s] ? 0 : (mat.epsr * EPS0 * VT) / h; // ε = 0: no displacement
      const D = -k * (phiR - phiL);
      const lastSeg = s === this.nNodes - 2;
      if (lastSeg) {
        // Displacement current through this segment; conduction is added per species below.
        this.segD = D;
        this.segI = Number.isFinite(dt) ? (D - this.segDOld) / dt : 0;
        this.segIJac.fill(0);
        if (Number.isFinite(dt)) {
          this.segIJac[0] = k / dt;
          this.segIJac[M] = -k / dt;
        }
      }
      res[bL * M] += D;
      res[bR * M] -= D;
      this._j(bL, 0, bL, 0, k);
      this._j(bL, 0, bR, 0, -k);
      this._j(bR, 0, bL, 0, -k);
      this._j(bR, 0, bR, 0, k);

      if (this.nodeMetal[s] >= 0) {
        this._segmentMetal(s, bL, bR, mat, h, lastSeg);
        continue;
      }
      const vel = regions[reg].velocity;
      if (regions[reg].mixing > 0) this._segmentMixing(s, bL, bR, regions[reg].mixing, h, mat, lastSeg);
      if (!this.nodeIdeal[s]) {
        this._segmentNonIdeal(s, bL, bR, mat, h, lastSeg, vel);
        continue;
      }
      for (let i = 0; i < n; i++) {
        if (!mat.present[i] || mat.D[i] === 0) continue;
        const r = 1 + i, zi = z[i];
        // SG flux N = g[B(Δ)c_L − B(−Δ)c_R], rewritten with B(−Δ) = B(Δ)e^Δ and
        // c_R e^Δ = c_L e^{Δη} as N = −g·B(Δ)·c_L·expm1(Δη). This is precise relative to the
        // quasi-Fermi difference Δη, so tiny fluxes (e.g. majority carriers carrying a small
        // current) don't vanish in the cancellation of two huge drift and diffusion terms.
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
        res[bL * M + r] += N;
        res[bR * M + r] -= N;
        if (lastSeg) {
          const q = F * zi;
          this.segI += q * N;
          this.segIJac[r] += q * dNdEtaL;
          this.segIJac[M + r] += q * dNdEtaR;
          this.segIJac[0] += q * dNdPhiL;
          this.segIJac[M] += q * dNdPhiR;
        }
        // Row r of block bL (couplings to itself in B, to bR in C) and of bR (to bL in A).
        const oL = bL * M * M + r * M, oR = bR * M * M + r * M;
        sys.B[oL + r] += dNdEtaL;
        sys.C[oL + r] += dNdEtaR;
        sys.B[oL] += dNdPhiL;
        sys.C[oL] += dNdPhiR;
        sys.A[oR + r] -= dNdEtaL;
        sys.B[oR + r] -= dNdEtaR;
        sys.A[oR] -= dNdPhiL;
        sys.B[oR] -= dNdPhiR;
      }
    }

    // Interfaces: the flux node carries D and N_i; its rows are the interface laws.
    for (let f = 0; f < this.nFaces; f++) {
      const bf = this.blockOfFace[f], bL = bf - 1, bR = bf + 1;
      const itf = interfaces[f];
      // φ law: pinned jump (dipole), Helmholtz capacitor, or no charge at all (neutral).
      const law = itf.phi.type;
      const metalSide = itf.metal && law !== 'neutral' ? itf.metal.side : null;
      if (metalSide) {
        // Against a metal: the other side's φ is tied to the metal's Fermi level V_F = V_T η/z.
        const im = itf.metal.i, bm = metalSide === 'left' ? bL : bR, bo = metalSide === 'left' ? bR : bL;
        const sg = metalSide === 'left' ? 1 : -1;
        const vf = (u[bm * M + 1 + im] + uLo[bm * M + 1 + im]) / z[im]; // V_F / V_T
        const gap = vf - itf.zeroCharge / VT - (u[bo * M] + uLo[bo * M]); // (V_F − zeroCharge − φ_edge)/V_T
        if (law === 'dipole') {
          res[bf * M] = gap;
          this._j(bf, 0, bm, 1 + im, 1 / z[im]);
          this._j(bf, 0, bo, 0, -1);
        } else {
          // D toward the other side = C (V_F − zeroCharge − φ_edge); toward +x that's sg times it.
          const kC = itf.phi.C * VT;
          res[bf * M] = u[bf * M] - sg * kC * gap;
          this._j(bf, 0, bf, 0, 1);
          this._j(bf, 0, bm, 1 + im, (-sg * kC) / z[im]);
          this._j(bf, 0, bo, 0, sg * kC);
        }
      } else if (law === 'neutral') {
        res[bf * M] = u[bf * M]; // D = 0; the jump is whatever each side's neutrality needs
        this._j(bf, 0, bf, 0, 1);
      } else {
        const jump = u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) - itf.dipole / VT;
        if (law === 'dipole') {
          res[bf * M] = jump;
          this._j(bf, 0, bR, 0, 1);
          this._j(bf, 0, bL, 0, -1);
        } else {
          // D = −C (φ_R − φ_L − dipole): displacement toward +x drops across the layer.
          const kC = itf.phi.C * VT;
          res[bf * M] = u[bf * M] + kC * jump;
          this._j(bf, 0, bf, 0, 1);
          this._j(bf, 0, bR, 0, kC);
          this._j(bf, 0, bL, 0, -kC);
        }
      }
      // D_f enters each side's Gauss row; a metal side holds it as surface carriers instead.
      const gL = grid.regionEnd[f], gR = gL + 1;
      if (this.nodeMetal[gL] < 0) {
        res[bL * M] += u[bf * M];
        this._j(bL, 0, bf, 0, 1);
      }
      if (this.nodeMetal[gR] < 0) {
        res[bR * M] -= u[bf * M] + itf.sheetCharge;
        this._j(bR, 0, bf, 0, -1);
      }
      for (let i = 0; i < n; i++) {
        const r = 1 + i, o = bf * M + r;
        const type = itf.links[i].type;
        const deta = u[bL * M + r] - u[bR * M + r] + (uLo[bL * M + r] - uLo[bR * M + r]); // η_L − η_R
        if (type === 'blocked') {
          res[o] = u[o];
          this._j(bf, r, bf, r, 1);
          continue;
        }
        if (type === 'equilibrium') {
          res[o] = -deta; // μ̄ continuous
          this._j(bf, r, bR, r, 1);
          this._j(bf, r, bL, r, -1);
        } else if (type === 'conductance') {
          // J = G (V_L − V_R), V = V_T η / z  ⇒  N = G V_T (η_L − η_R) / (z² F)
          const gG = (itf.links[i].G * VT) / (z[i] * z[i] * F);
          res[o] = u[o] - gG * deta;
          this._j(bf, r, bf, r, 1);
          this._j(bf, r, bL, r, -gG);
          this._j(bf, r, bR, r, gG);
        } else {
          res[o] = u[o]; // kinetic: N − Σ ν r, the rates are subtracted below
          this._j(bf, r, bf, r, 1);
        }
        // (a species on one side only, e.g. at an electrode reaction, flows on that side only)
        if (this.present[gL * n + i]) {
          res[bL * M + r] += u[o];
          this._j(bL, r, bf, r, 1);
        }
        if (this.present[gR * n + i]) {
          res[bR * M + r] -= u[o];
          this._j(bR, r, bf, r, -1);
        }
      }
      for (const tr of itf.transfers) this._transfer(tr, f, bf, bL, bR);
      for (const rx of itf.electrode) this._electrodeAtFace(rx, itf.metal, f, bf, bL, bR);
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
      res[tb * M] = this._termI - I;
      for (let r = 0; r < M; r++) this._j(tb, 0, tb - 1, r, this._termIJac[r]);
      this._j(tb, 0, tb, 0, (this._termIJac[M] - dIdV) * VT);
      for (let r = 1; r < M; r++) {
        res[tb * M + r] = u[tb * M + r];
        this._j(tb, r, tb, r, 1);
      }
    }
  }

  // Scharfetter–Gummel with non-ideal statistics. The excess ex = ζ − ln(c/c_ref) acts as an
  // extra potential, linear along the segment like φ, so Δ = zΔφ̂ + Δex and
  //   N = −(D/h)·B(Δ)·c_L·expm1(η_R − η_L),
  // still exactly zero at equilibrium. c_L and ex depend on every ζ at their node through K.
  _segmentNonIdeal(s, bL, bR, mat, h, lastSeg, vel) {
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
      res[bL * M + r] += N;
      res[bR * M + r] -= N;
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
  _segmentMetal(s, bL, bR, mat, h, lastSeg) {
    const { M, u, uLo, res, z } = this;
    const i = mat.metal.i, r = 1 + i;
    const g = (mat.metal.sigma * this.model.RT) / (z[i] * z[i] * FARADAY * FARADAY * h);
    const J = u[bL * M];
    res[bL * M] = u[bR * M + r] - u[bL * M + r] + (uLo[bR * M + r] - uLo[bL * M + r]) + J / g;
    this._j(bL, 0, bR, r, 1);
    this._j(bL, 0, bL, r, -1);
    this._j(bL, 0, bL, 0, 1 / g);
    res[bL * M + r] += J;
    res[bR * M + r] -= J;
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
      res[bL * M + r] += N;
      res[bR * M + r] -= N;
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

  // Electrode reaction at a face between a metal and another material, written as reduction:
  //   Σ ν_R R + n e⁻(metal, at its Fermi level) ⇌ Σ ν_P P,
  // with the same standard-rate-constant form as at contacts, participants at the other side's
  // edge node (behind any double layer, so Frumkin effects arise). Each species flows to or from
  // the face on its own side: with s = +1 for a metal on the left, the face fluxes toward +x are
  // −sνr (reactants), +sνr (products) and +s·n·r (the metal's electrons).
  _electrodeAtFace(rx, metal, f, bf, bL, bR) {
    const { n, M, u, uLo, c, res } = this;
    const grid = this.model.grid;
    const sg = metal.side === 'left' ? 1 : -1;
    const bm = sg > 0 ? bL : bR, bo = sg > 0 ? bR : bL;
    const go = sg > 0 ? grid.regionStart[f + 1] : grid.regionEnd[f];
    const al = rx.alpha, ie = metal.i;
    const dr = this.dA.fill(0);
    let pref = rx.k0, aHi = rx.fixedA + rx.electrons * u[bm * M + 1 + ie], aLo = rx.electrons * uLo[bm * M + 1 + ie];
    for (const { i, nu } of rx.reactants) {
      pref *= (c[go * n + i] / this.cRef[go * n + i]) ** (nu * (1 - al));
      this._dlnc(go, i, nu * (1 - al), dr);
      aHi += nu * u[bo * M + 1 + i];
      aLo += nu * uLo[bo * M + 1 + i];
    }
    for (const { i, nu } of rx.products) {
      pref *= (c[go * n + i] / this.cRef[go * n + i]) ** (nu * al);
      this._dlnc(go, i, nu * al, dr);
      aHi -= nu * u[bo * M + 1 + i];
      aLo -= nu * uLo[bo * M + 1 + i];
    }
    const { g, gp } = bvFactor(aHi + aLo, al);
    const rate = pref * g;
    for (let k = 0; k < M; k++) dr[k] *= rate;
    for (const { i, nu } of rx.reactants) dr[1 + i] += pref * gp * nu;
    for (const { i, nu } of rx.products) dr[1 + i] -= pref * gp * nu;
    const dE = pref * gp * rx.electrons; // ∂rate/∂η_e(metal)
    const rows = [...rx.reactants.map(({ i, nu }) => [i, -sg * nu]), ...rx.products.map(({ i, nu }) => [i, sg * nu]), [ie, sg * rx.electrons]];
    for (const [i, w] of rows) {
      const row = 1 + i;
      res[bf * M + row] -= w * rate;
      for (let k = 0; k < M; k++) if (dr[k] !== 0) this._j(bf, row, bo, k, -w * dr[k]);
      this._j(bf, row, bm, 1 + ie, -w * dE);
    }
  }

  // Kinetic transfer across interface f (Butler–Volmer form, forward = left to right):
  //   r = k0 Π [(c_L/c_ref,L)^{ν(1−α)} (c_R/c_ref,R)^{να}] (e^{αa} − e^{−(1−α)a}),
  //   a = Σ ν (η_L − η_R).  Each species' flux-node row gets −ν·r.
  _transfer(tr, f, bf, bL, bR) {
    const { n, M, u, uLo, c, res } = this;
    const gL = this.model.grid.regionEnd[f], gR = this.model.grid.regionStart[f + 1];
    const al = tr.alpha;
    const dL = this.dA.fill(0), dR = this.dB.fill(0); // ∂ ln(prefactor) per slot, each side
    let pref = tr.k0, aHi = 0, aLo = 0;
    for (const { i, nu } of tr.species) {
      pref *= (c[gL * n + i] / this.cRef[gL * n + i]) ** (nu * (1 - al)) * (c[gR * n + i] / this.cRef[gR * n + i]) ** (nu * al);
      this._dlnc(gL, i, nu * (1 - al), dL);
      this._dlnc(gR, i, nu * al, dR);
      aHi += nu * (u[bL * M + 1 + i] - u[bR * M + 1 + i]);
      aLo += nu * (uLo[bL * M + 1 + i] - uLo[bR * M + 1 + i]);
    }
    const { g, gp } = bvFactor(aHi + aLo, al);
    const rate = pref * g;
    for (let k = 0; k < M; k++) {
      dL[k] *= rate;
      dR[k] *= rate;
    }
    for (const { i, nu } of tr.species) {
      dL[1 + i] += pref * gp * nu;
      dR[1 + i] -= pref * gp * nu;
    }
    for (const { i, nu } of tr.species) {
      const row = 1 + i;
      res[bf * M + row] -= nu * rate;
      for (let k = 0; k < M; k++) {
        if (dL[k] !== 0) this._j(bf, row, bL, k, -nu * dL[k]);
        if (dR[k] !== 0) this._j(bf, row, bR, k, -nu * dR[k]);
      }
    }
  }

  // An internal port's exchange with each node of its window, as a source per volume. A held
  // ('equilibrium') level replaces the node's balance row; the source is then that row's residual.
  // Ports come before the contacts, so a contact's flux readout includes a port's source there.
  _port(port, flux) {
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
          flux[i] += res[o];
          this._replaceRow(b, r);
          this._j(b, r, b, r, 1);
          res[o] = -deta;
          continue;
        }
        // conductance: s = G V_T (η_out − η)/(z² F); exchange: s = k (η_out − η)
        const kk = link.type === 'conductance' ? (link.G * VT) / (z[i] * z[i] * F) : link.k;
        res[o] -= v * kk * deta;
        this._j(b, r, b, r, v * kk);
        flux[i] += v * kk * deta;
      }
    }
  }

  // One contact: record the flux through its outer face, then add its exchange terms (electrode
  // reactions, conductance links), then apply equilibrium links and the φ law.
  _contact(side, dt) {
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
    for (let i = 0; i < n; i++) flux[i] = sgn * res[b * M + 1 + i];

    // Electrode reactions: Σ ν_R R + n e⁻(metal, μ̄ = −F V_t) ⇌ Σ ν_P P.
    const dr = this.dA;
    for (const rx of ct.reactions) {
      const al = rx.alpha;
      let pref = rx.k0, aHi = rx.fixedA - (rx.electrons * Vt) / VT, aLo = 0;
      dr.fill(0);
      for (const { i, nu } of rx.reactants) {
        pref *= (c[g * n + i] / this.cRef[g * n + i]) ** (nu * (1 - al));
        this._dlnc(g, i, nu * (1 - al), dr);
        aHi += nu * u[b * M + 1 + i];
        aLo += nu * uLo[b * M + 1 + i];
      }
      for (const { i, nu } of rx.products) {
        pref *= (c[g * n + i] / this.cRef[g * n + i]) ** (nu * al);
        this._dlnc(g, i, nu * al, dr);
        aHi -= nu * u[b * M + 1 + i];
        aLo -= nu * uLo[b * M + 1 + i];
      }
      const { g: gf, gp } = bvFactor(aHi + aLo, al);
      const rate = pref * gf;
      const dVt = (pref * gp * -rx.electrons) / VT;
      // ∂rate/∂slot = rate·∂ln(prefactor) + pref·g′·∂a
      for (let k = 0; k < M; k++) dr[k] *= rate;
      for (const { i, nu } of rx.reactants) dr[1 + i] += pref * gp * nu;
      for (const { i, nu } of rx.products) dr[1 + i] -= pref * gp * nu;
      const rows = [...rx.reactants.map(({ i, nu }) => [i, nu]), ...rx.products.map(({ i, nu }) => [i, -nu])];
      for (const [i, w] of rows) {
        const rs = 1 + i; // consumed (w > 0) or produced (w < 0) at this node
        res[b * M + rs] += w * rate;
        for (let k = 0; k < M; k++) if (dr[k] !== 0) this._j(b, rs, b, k, w * dr[k]);
        addVt(rs, w * dVt);
      }
      if (toTerminal) {
        // Electrons taken from the metal at the right: current toward +x of n F r.
        const q = rx.electrons * F;
        this._termI += q * rate;
        for (let k = 0; k < M; k++) termJ[k] += q * dr[k];
        termJ[M] += q * dVt;
      }
    }

    // Conductance links: J (toward the device) = G (V_out − V_i), V_out = V_t + offset.
    for (let i = 0; i < n; i++) {
      const link = ct.species[i];
      if (link.type !== 'conductance') continue;
      const o = b * M + 1 + i;
      const Vi = (VT * (u[o] + uLo[o])) / z[i];
      const Nin = (link.G * (Vt + link.offset - Vi)) / (z[i] * F); // particles entering
      const k = (link.G * VT) / (z[i] * z[i] * F);
      res[o] -= Nin;
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
        res[o] = this.segI - I;
        for (let r = 0; r < M; r++) {
          this._j(b, 1 + i, b - 1, r, this.segIJac[r]);
          this._j(b, 1 + i, b, r, this.segIJac[M + r]);
        }
        if (circuit.mode === 'load') this._j(b, 1 + i, b, 1 + t, -VT / (z[t] * circuit.R));
      } else if (floating && ct.species[i].mu === undefined) {
        // Charged species tied to the floating terminal: η_i = z_i (V_t + offset_i)/V_T.
        res[o] = u[o] + uLo[o] - (z[i] * (Vt + ct.species[i].offset)) / VT;
        this._j(b, 1 + i, b, 1 + i, 1);
        addVt(1 + i, -z[i] / VT);
      } else {
        this._j(b, 1 + i, b, 1 + i, 1);
        res[o] = u[o] - this.contactEta(side, i) + uLo[o];
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
      res[b * M] -= sgn * D;
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
      this.contactD[side] = sgn * res[b * M];
      this._replaceRow(b, 0);
      this._j(b, 0, b, 0, 1);
      res[b * M] = u[b * M] + uLo[b * M] - (Vt - link.zeroCharge) / VT;
      addVt(0, -1 / VT);
    } else if (link.type === 'bulk') {
      this.contactD[side] = sgn * res[b * M];
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
      res[b * M] = q;
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
  _replaceRow(b, r) {
    const M = this.M, o = b * M * M + r * M;
    for (let k = 0; k < M; k++) {
      this.sys.A[o + k] = 0;
      this.sys.B[o + k] = 0;
      this.sys.C[o + k] = 0;
    }
  }

  // Scale every row by its largest Jacobian entry (in place, residual too).
  _equilibrate() {
    const { M, nB, sys, res } = this;
    const { A, B, C } = sys;
    const mm = M * M;
    for (let b = 0; b < nB; b++) {
      for (let r = 0; r < M; r++) {
        const o = b * mm + r * M;
        let mx = 0;
        for (let k = 0; k < M; k++) {
          mx = Math.max(mx, Math.abs(A[o + k]), Math.abs(B[o + k]), Math.abs(C[o + k]));
        }
        if (mx === 0) continue;
        const s = 1 / mx;
        for (let k = 0; k < M; k++) {
          A[o + k] *= s;
          B[o + k] *= s;
          C[o + k] *= s;
        }
        res[b * M + r] *= s;
      }
    }
  }

  // Largest update among the potential-like unknowns (φ̂ and η at grid nodes).
  _maxPotentialStep(delta) {
    const { n, M } = this;
    let mx = this.terminalBlock >= 0 ? Math.abs(delta[this.terminalBlock * M]) : 0;
    for (let g = 0; g < this.nNodes; g++) {
      const b = this.blockOfNode[g];
      if (this.nodeMetal[g] < 0) mx = Math.max(mx, Math.abs(delta[b * M])); // (a metal's slot 0 is a flux)
      for (let i = 0; i < n; i++) {
        if (this.present[g * n + i]) mx = Math.max(mx, Math.abs(delta[b * M + 1 + i]));
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
      for (let k = 0; k < res.length; k++) rmax = Math.max(rmax, Math.abs(res[k]));
      try {
        this.sys.factor();
      } catch (err) {
        return { converged: false, iterations: it, history, error: err.message };
      }
      this.sys.solve(res, delta);
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
    const { u, uLo } = this;
    for (let k = 0; k < u.length; k++) {
      const a = u[k], b = scale * d[k];
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
    const N = nB * M, mm = M * M;
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

    const csys = new ComplexBlockTridiagonal(nB, M);
    const rr = new Float64Array(N), ri = new Float64Array(N), xr = new Float64Array(N), xi = new Float64Array(N);
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
        for (let r = 0; r < M; r++) {
          const o = blk * mm + r * M;
          let mx = 0;
          for (const X of ['A', 'B', 'C']) for (let c = 0; c < M; c++) mx = Math.max(mx, Math.abs(J[X][o + c]), w * Math.abs(S[X][o + c]));
          const sc = mx > 0 ? 1 / mx : 1;
          for (const X of ['A', 'B', 'C']) {
            for (let c = 0; c < M; c++) {
              csys[X + 'r'][o + c] = J[X][o + c] * sc;
              csys[X + 'i'][o + c] = w * S[X][o + c] * sc;
            }
          }
          rr[blk * M + r] = -b[blk * M + r] * sc;
          ri[blk * M + r] = 0;
        }
      }
      csys.factor();
      csys.solve(rr, ri, xr, xi);
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
    const direct = this.stretches.every((st) => st.connected);
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
    const sub = this.stretches.every((st) => st.connected) ? { ...opts, maxSteps: 1 } : opts;
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

  _solveSteady({ maxSteps = 80, tol = 1e-11 } = {}) {
    const time = this.time;
    const tau = this.slowestTime();
    const direct = this.stretches.every((st) => st.connected);
    const giant = direct ? Infinity : 1e6 * tau;
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
    for (let k = 0; k < d.length; k++) d[k] = this.u[k] - this.uPrev[k];
    return d;
  }
}
