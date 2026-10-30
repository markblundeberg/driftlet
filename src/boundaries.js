// The equations at the boundaries: each face between regions (its φ law, the species' links
// across it, its reactions and gates), and the terminals (the ports inside the device, the two
// contacts, and the circuit rows of the floating ones), with what each passes into the device.
//
// CONTRACT
//
// Methods of the Solver, called by assemble() (src/assemble.js, whose contract sets the rows'
// conventions) and its bookkeeping, after the regions' terms (src/regions.js), whose
// concentrations and derivatives (_dlnc, _Kz) they use. Units: a face's slots are per area (its
// flux N_i toward +x, its D, its reactions' rates), and enter the edge nodes' rows times the
// face's area Af; a port's terms are per volume of its window (times vol·fill, the part of each
// node's box inside the window: port.fill, src/device.js), its reactions' per area of electrode
// (times vol and the electrode's area per volume of each node's box, portArea, fill in it); a contact's
// links are per area, times the area at the end node. V_T = RT/F; a terminal's voltage V (termV)
// is the outside's level per charge.
//
// A face f (_face; its flux block bf between edge nodes bL, bR):
//   φ law: pinned, φ̂_R − φ̂_L − dipole/V_T = 0; capacitive, D + C V_T (φ̂_R − φ̂_L − dipole/V_T)
//     = 0; neutral, no D. Against a metal (itf.conductor, not neutral), the other side's edge φ
//     is tied to the metal's Fermi level V_F = V_T η/z instead: capacitive, D toward the other
//     side = C (V_F − zeroCharge − φ) (pinned, V_F − zeroCharge − φ = 0, is written too, but
//     normalizeDevice refuses it: a metal holds its surface charge in a finite capacitance). D enters the left Gauss row
//     as +Af·D, the right one's as −Af·(D + sheetCharge); a metal side's charge is its carriers'
//     (src/regions.js).
//   links, per species: blocked, no flux; equilibrium, η_L = η_R (the flux whatever the
//     balances need); permeability, N = −P B(d) c_L expm1(η_R − η_L), d = z Δφ̂ + the step in
//     μ°/RT − ln c_ref (Goldman–Hodgkin–Katz between like solutions); conductance, N = G V_T
//     (η_L − η_R)/(z² F). Gated, P or G times Π x^p over its gates, each x taken in [0, 1]. N
//     enters the left node's balance as +Af·N, the right's as −Af·N.
//   reactions (_faceReaction): its rate r an unknown, row r − rate = 0, with a = fixedA − Σ ν η
//     over its participants (ν > 0 made, at either side's edge node); rate Butler–Volmer, k0
//     Π (c/c_ref)^e g(a) (e = |ν|(1 − α) for what it consumes, να for what it makes; a metal's
//     carrier at activity 1), or Marcus–Hush–Chidsey (lambda: k0 Π_{ν<0} (c/c_ref)^{|ν|} F(a°)
//     (1 − e^{−a}), src/marcus.js), or saturating (vmax), or SRH; −Af·ν·r into each participant's row.
//   gates (_faceGates): (x − x_old)/dt − α(V)(1 − x) + β(V)x = 0, V = φ_R − φ_L (no storage at
//     dt = ∞).
//   diffuse layers (_diffuseTerms; itf.diffuse, on each strictly neutral side: s.diffList, its
//     q, Γ and their derivatives from computeConcentrations, src/diffuse.js): each a slot after
//     the gates, its drop s = Fψ/RT, which puts that side's plane at φ̂_edge + s: the φ law (each
//     form above) runs to the planes, φ̂ → φ̂ + s. Its row: q(s, c_edge) − Q, Q the charge the face
//     gives that side (D on the left, −(D + sheetCharge) on the right); −Af·q in the edge node's
//     Gauss row, and Af (Γ_i − Γ_i,old)/dt in its balances (none at dt = ∞). A reaction's
//     participant on that side reacts at the plane: its c times e^{−z s} in the rate (Butler–Volmer
//     and saturating; SRH beside a layer is refused by normalizeDevice), and its level there w
//     above the edge node's (a slot after the drops, _diffuseDraws: its η + w, its c times e^w;
//     row 1 − e^w − N K/(D c), N what the reactions take from that side, K the layer's
//     resistance to it; where K < 0, −w − X/√(1 + X²), X = N K/(D c)); a
//     permeability link runs between the planes (c e^{−z s} at each, its drift over φ̂ + s), and
//     gates follow the voltage between them (_faceVoltage).
//
// The terminals (_terminals), in order: the sources refreshed (_refreshSources); the live
// terminals' termB, termC zeroed, and termDI, termDX, termI, qg, qgSlope; each port's terms but
// its held levels, then each port's held levels; the left contact, the right; each floating
// terminal's circuit row.
//   port k (terminal 2 + k, _port), at each node of its window: an equilibrium link holds η at
//     the port's level (portEta) at the nodes inside the window (port.inside), replacing the
//     balance row, whose residual just before is the flux it passes (not at a device end node
//     whose contact links that species); conductance
//     and exchange links supply vol·fill·k(η_out − η), k = G V_T/(z² F) or the link's k (gated as at
//     a face); reactions (_portReaction) as at a face against the electrode's carrier at the
//     port's level, explicit (no unknown), −vol·area·ν·rate into the region's species' rows,
//     −ν·rate into a surface species' (its rows are per area of electrode); a capacitance (_portCapacitance), σ =
//     C (V − zeroCharge − φ) per electrode area, −vol·area·σ into each Gauss row (not at an end
//     node whose contact sets φ, nor where φ is undefined), its charge portQ and per-node qg;
//     gates as at a face, in V_m = φ − (V − zeroCharge); a surface's storage, Γ(θ − θ_old)/dt.
//     termI: Σ z F × what it passes, plus its capacitance's charging, (Q − Q at the step's
//     start)/dt (steady at no current: Q − its captured charge).
//   contact (_contact): first, at the end node, each active balance residual (all of the box's
//     other terms) is the flux the contact passes in: contactFlux (toward +x), and the current
//     Σ z F × that; then, at a lateStorage node on a time step, the ions' storage (in the rows
//     and contactFlux, not the current); then its links: exchange, A k (μ_out/RT − η); conductance,
//     A G (V + offset − V_i)/(zF), V_i = V_T η/z; velocity, A v c (e^{η_out − η} − 1) (dilute,
//     A v (c_eq − c)); equilibrium, the row replaced by η = η_out (contactEta), any other
//     terminal's ∂/∂V in it cleared. Its φ law: capacitive, D_in = A C (V − zeroCharge − V_T φ̂)
//     (φ̂ as u alone, as the displacement reads it) out of the Gauss row; pinned or bulk, D_in =
//     the Gauss row's residual so far, the row replaced by φ̂ = (V − zeroCharge)/V_T or by
//     neutrality (F Σ z c + ρ_fixed = 0, its ∂/∂η_i = F (Kz)_i, K being symmetric); the current gains (D_in − D_in at the step's start)/dt (not at a lateStorage
//     node on a time step; steady at no current: D_in − its captured charge); contactD = ±D_in
//     toward +x. termI: that current.
//   A terminal whose current is read from rows (a contact's from its end box, a port's held
//     level from the rows it replaces) takes their ∂/∂V's with it (_captureRow): another live
//     terminal j's termB entries there into termDX[k][j] (a port's window in a contact's end
//     box), its own into termDI[k].
//   circuit (_circuit), floating k: termRes = I − I_set, or I − (V_src − V)/R.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. Every Jacobian entry is its residual's derivative; termB[k] is ∂res/∂V_k, termC[k]
//      ∂termI_k/∂(unknowns), and termDI[k], termDX[k][j] ∂termI_k/∂V_k and /∂V_j (with a circuit
//      row's own), for the live terminals (_liveTerms).
//   2. termI[k] is the current into the device through terminal k: what its links and reactions
//      pass, times z F, and its displacement current; a contact's equals what the end box needs
//      from outside, read from its residual (so the current is conserved through the device).
//   3. A face's flux enters its two edge nodes with opposite signs, times Af.
//   4. Where μ̄ is level across a link (or a = 0 for a reaction), its flux (rate) is exactly 0.
//   4b. A diffuse layer's side holds the face's charge in the layer: at its row's solution the
//      edge box's own charge (vol·ρ) is what its other terms put there, as without the face.
//   5. The readouts (contactFlux, contactD, portFlux, portQ, portRates) are those of the state
//      assembled, written at every assembly.
// Tests: test/contacts.test.js, test/ports.test.js, test/terminals.test.js, test/circuit.test.js,
// test/gates.test.js, test/membrane.test.js, test/kinetics.test.js, test/coverage.test.js,
// test/capacitance.test.js, test/metal.test.js, test/diffuse.test.js, test/jacobian.test.js.

import { sgFlux } from './bernoulli.js';
import { FARADAY } from './constants.js';
import { sourceAt } from './device.js';
import { diffuseLayer } from './diffuse.js';
import { mhcIntegral } from './marcus.js';
import { powi, powr } from './pow.js';

// Butler–Volmer factor g(a) = e^{αa} − e^{−(1−α)a} = e^{−(1−α)a}·expm1(a) (precise near a = 0),
// and its derivative g′(a) = α e^{αa} + (1−α) e^{−(1−α)a}.
function bvFactor(a, alpha) {
  const e = Math.exp(-(1 - alpha) * a);
  // (far from equilibrium, the plain difference: the product would give 0·∞ there)
  const g = Math.abs(a) > 50 ? Math.exp(alpha * a) - e : e * Math.expm1(a);
  return { g, gp: alpha * Math.exp(alpha * a) + (1 - alpha) * e };
}

// A gate's rate law (see normalizeGates in device.js) at voltage V: [rate, ∂rate/∂V], 1/s and 1/(s·V).
export function gateRate(law, V) {
  const { rate, midpoint, scale } = law;
  if (law.type === 'exp') {
    // (held finite where Newton's iterates stray to absurd voltages: β_m at −13 V would overflow)
    const y = (V - midpoint) / scale, v = rate * Math.exp(Math.min(y, 700));
    return [v, y > 700 ? 0 : v / scale];
  }
  if (law.type === 'sigmoid') {
    const e = Math.exp((midpoint - V) / scale), d = 1 + e;
    return Number.isFinite(e) ? [rate / d, (rate * e) / (scale * d * d)] : [0, 0];
  }
  // expLinear: rate·y/(1 − e^−y), smooth through y = 0 (rate there) and → rate·y for large y
  const y = (V - midpoint) / scale;
  if (Math.abs(y) < 1e-6) return [rate * (1 + y / 2), rate / (2 * scale)];
  if (y < -700) return [0, 0];
  const em = -Math.expm1(-y), ey = Math.exp(-y); // 1 − e^−y, e^−y
  return [(rate * y) / em, (rate * (em - y * ey)) / (em * em * scale)];
}

export class BoundaryTerms {
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
      pref *= p.side === 2 ? powr(th[g * nSurf + p.s], e) : powr(c[g * n + p.i] / this.cRef[g * n + p.i], e);
    }
    if (rx.sites) pref *= powr(this.th0[g], rx.sites);
    return pref * bvFactor(a, rx.alpha).g;
  }

  // A face: its flux node's rows are the interface laws (φ law, species links, reactions), and
  // its fluxes enter the two edge nodes' balances.
  _face(f) {
    const { model, n, M, u, uLo, res, z, VT } = this, R = this.rix, F = FARADAY;
    const grid = model.grid, interfaces = model.interfaces;
    const bf = this.blockOfFace[f], bL = bf - 1, bR = bf + 1;
    const itf = interfaces[f];
    // A diffuse layer's drop s on a side (src/diffuse.js) puts that side's plane, where the Stern
    // layer ends, at φ̂_edge + s: [left, right] layer indices (−1: none) and slots.
    const [dlL, dlR] = this.diffAt[f], dl = this.diffList;
    const sOf = (d) => (d < 0 ? 0 : u[dl[d].o]);
    const jS = (rs, d, v) => d >= 0 && this._j(bf, rs, bf, dl[d].slot, v);
    // φ law: pinned jump (dipole), Helmholtz capacitor, or no charge at all (neutral).
    const law = itf.phi.type;
    const condSide = itf.conductor && law !== 'neutral' ? itf.conductor.side : null;
    if (condSide) {
      // Against a metal: the other side's φ is tied to the metal's Fermi level V_F = V_T η/z.
      const im = itf.conductor.i, bm = condSide === 'left' ? bL : bR, bo = condSide === 'left' ? bR : bL;
      const sg = condSide === 'left' ? 1 : -1, dlo = condSide === 'left' ? dlR : dlL;
      const vf = (u[bm * M + 1 + im] + uLo[bm * M + 1 + im]) / z[im]; // V_F / V_T
      const gap = vf - itf.zeroCharge / VT - (u[bo * M] + uLo[bo * M]) - sOf(dlo); // (V_F − zeroCharge − φ_plane)/V_T
      if (law === 'pinned') {
        res[R[bf * M]] = gap;
        this._j(bf, 0, bm, 1 + im, 1 / z[im]);
        this._j(bf, 0, bo, 0, -1);
        jS(0, dlo, -1);
      } else {
        // D toward the other side = C (V_F − zeroCharge − φ_plane); toward +x that's sg times it.
        const kC = itf.phi.C * VT;
        res[R[bf * M]] = u[bf * M] - sg * kC * gap;
        this._j(bf, 0, bf, 0, 1);
        this._j(bf, 0, bm, 1 + im, (-sg * kC) / z[im]);
        this._j(bf, 0, bo, 0, sg * kC);
        jS(0, dlo, sg * kC);
      }
    } else if (law === 'neutral') {
      // D = 0 (not an unknown); the jump is whatever each side's neutrality needs
    } else {
      const jump = u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) + sOf(dlR) - sOf(dlL) - itf.dipole / VT;
      if (law === 'pinned') {
        res[R[bf * M]] = jump;
        this._j(bf, 0, bR, 0, 1);
        this._j(bf, 0, bL, 0, -1);
        jS(0, dlR, 1);
        jS(0, dlL, -1);
      } else {
        // D = −C (φ_R − φ_L − dipole): displacement toward +x drops across the layer.
        const kC = itf.phi.C * VT;
        res[R[bf * M]] = u[bf * M] + kC * jump;
        this._j(bf, 0, bf, 0, 1);
        this._j(bf, 0, bR, 0, kC);
        this._j(bf, 0, bL, 0, -kC);
        jS(0, dlR, kC);
        jS(0, dlL, -kC);
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
    if (dlL >= 0) this._diffuseTerms(dlL, bf, Af);
    if (dlR >= 0) this._diffuseTerms(dlR, bf, Af);
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
        // Beside a diffuse layer, between its plane's: c e^{−z s} at each, the drift over the
        // drop between the planes (s acting as φ̂ does on its side).
        const gl = grid.regionEnd[f], gr = gl + 1, kl = gl * n + i, kr = gr * n + i, zi = z[i], link = itf.links[i];
        const shift = this.mu0hat[kr] - this.mu0hat[kl] - Math.log(this.cRef[kr] / this.cRef[kl]);
        const sL = sOf(dlL), sR = sOf(dlR), cl = this.c[kl] * (sL === 0 ? 1 : Math.exp(-zi * sL)), cr = this.c[kr] * (sR === 0 ? 1 : Math.exp(-zi * sR));
        const d = zi * (u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) + sR - sL) + shift;
        // Gated: P × Π x^p, each x taken in [0, 1] (Newton's iterates can stray past its ends).
        const gated = link.gates ?? [], gx = (q) => u[this._gateSlot(f, q)];
        const N1 = gated.length ? sgFlux(link.P, d, -deta, cl, cr, this.sg)[0] : 0; // (fully open)
        const P = link.P * (gated.length ? this._open(gated, gx) : 1);
        const sg = sgFlux(P, d, -deta, cl, cr, this.sg); // (−deta = η_R − η_L)
        const N = sg[0], dNdd = sg[1];
        res[R[o]] = u[o] - N;
        this._j(bf, r, bf, r, 1);
        for (const [q, p] of gated) {
          const dO = this._dOpen(gated, q, p, gx);
          if (dO !== 0) this._j(bf, r, bf, 1 + n + itf.reactions.length + q, -N1 * dO);
        }
        this._j(bf, r, bL, r, -sg[2]);
        this._j(bf, r, bR, r, sg[3]);
        if (zi !== 0) {
          this._j(bf, r, bL, 0, zi * dNdd + zi * N);
          this._j(bf, r, bR, 0, -zi * dNdd);
          jS(r, dlL, zi * dNdd + zi * N);
          jS(r, dlR, -zi * dNdd);
        }
      } else {
        // conductance: J = G (V_L − V_R), V = V_T η / z  ⇒  N = G V_T (η_L − η_R) / (z² F); gated,
        // G times Π x^p
        const gated = itf.links[i].gates ?? [], gx = (q) => u[this._gateSlot(f, q)];
        const gG0 = (itf.links[i].G * VT) / (z[i] * z[i] * F), gG = gG0 * (gated.length ? this._open(gated, gx) : 1);
        res[R[o]] = u[o] - gG * deta;
        this._j(bf, r, bf, r, 1);
        this._j(bf, r, bL, r, -gG);
        this._j(bf, r, bR, r, gG);
        for (const [q, p] of gated) {
          const dO = this._dOpen(gated, q, p, gx);
          if (dO !== 0) this._j(bf, r, bf, 1 + n + itf.reactions.length + q, -gG0 * deta * dO);
        }
      }
      res[R[bL * M + r]] += Af * u[o];
      this._j(bL, r, bf, r, Af);
      res[R[bR * M + r]] -= Af * u[o];
      this._j(bR, r, bf, r, -Af);
    }
    itf.reactions.forEach((rx, k) => this._faceReaction(rx, k, f, bf, bL, bR, Af));
    if (itf.gates.length > 0) this._faceGates(f, bf, bL, bR);
    if (dlL >= 0) this._diffuseDraws(dlL, bf);
    if (dlR >= 0) this._diffuseDraws(dlR, bf);
  }

  // Layer d's drawn species (those its face's reactions take or make at its plane): each one's
  // level there, w = ln(a_plane/a_edge), below the edge node's by the flux N its reactions draw
  // across the layer's resistance K (src/diffuse.js), D c (1 − e^w) = N K: its row, 1 − e^w − N K/(D c).
  // Where it repels the species, K grows as c/c(plane), and the plane can run out of it (at
  // N = D c/K). Where it gathers one, K < 0, bounded by about the layer's width (the bulk the grid
  // counts in its place conducts it less well): a correction of order λ_D/L, and a negative
  // resistance in series, which far from the solution can draw a product's plane level down the
  // more of it is made; its row is bounded there (below).
  _diffuseDraws(d, bf) {
    const { n, M, u, c, res } = this, R = this.rix;
    const { b, g, side, f, dl, drawSlot, Dm } = this.diffList[d], itf = this.model.interfaces[f];
    if (!drawSlot.some((q) => q >= 0)) return;
    const dc = this._diffuseChain(d), K0 = 1 + n;
    for (let i = 0; i < n; i++) {
      const q = drawSlot[i];
      if (q < 0) continue;
      // N, what the reactions take from this side (mol/(m²·s)), from their rates' slots
      let N = 0;
      itf.reactions.forEach((rx, k) => {
        for (const p of rx.part) if (p.side === side && p.i === i) N -= p.nu * u[bf * M + K0 + k];
      });
      // X = N K/(D c); where the layer gathers the species (K < 0) the row is w = −X/√(1 + X²)
      // instead, the same to O(X³) where it's solved (|X| about λ_D/L there) and bounded far off
      const K = dl.K[i], w = u[bf * M + q], cc = Math.max(c[g * n + i], 1e-300), Dc = Dm[i] * cc, X = (N * K) / Dc;
      const gathered = K < 0, root = Math.sqrt(1 + X * X), gX = gathered ? 1 / (root * root * root) : 1;
      res[R[bf * M + q]] = gathered ? -w - X / root : -Math.expm1(w) - X;
      this._j(bf, q, bf, q, gathered ? -1 : -Math.exp(Math.min(700, w)));
      itf.reactions.forEach((rx, k) => {
        let nu = 0;
        for (const p of rx.part) if (p.side === side && p.i === i) nu += p.nu;
        if (nu !== 0) this._j(bf, q, bf, K0 + k, (gX * nu * K) / Dc);
      });
      this._j(bf, q, bf, this.diffList[d].slot, (-gX * N * dl.dKds[i]) / Dc);
      // through the edge node's c: K's ∂/∂c_j, and 1/c_i's
      for (let r = 0; r <= n; r++) {
        let v = (K * dc[i * M + r]) / cc;
        for (let j = 0; j < n; j++) v -= dl.dKdc[i * n + j] * dc[j * M + r];
        if (v !== 0) this._j(bf, q, b, r, (gX * N * v) / Dc);
      }
    }
  }

  // Diffuse layer d's terms (src/diffuse.js; its q and Γ from computeConcentrations): its row,
  // q(s, c) equal to the charge its side holds (D on the left, −(D + sheetCharge) on the right);
  // its charge in the edge box's Gauss row and its excesses in the box's balances, stored like the
  // box's own ions (so the box's bulk stays neutral, and what the face's charge draws in is
  // counted). The edge node's c enters through its ∂c_j/∂(φ̂, η) (_diffuseChain).
  _diffuseTerms(d, bf, Af) {
    const { n, M, u, res, gam, gamOld } = this, R = this.rix, dt = this.dtNow;
    const { b, slot, o, side, f, dl, holds } = this.diffList[d], dc = this._diffuseChain(d);
    const sheet = this.model.interfaces[f].sheetCharge, w = this.dW ?? (this.dW = new Float64Array(M));
    // ∂(Σ_j a_j c_j)/∂(the edge node's slot r), into w
    const chain = (a) => {
      w.fill(0);
      for (let j = 0; j < n; j++) if (a[j] !== 0) for (let r = 0; r <= n; r++) w[r] += a[j] * dc[j * M + r];
      return w;
    };
    // its row (per area): q − Q
    res[R[o]] = dl.q - (side ? -(u[bf * M] + sheet) : u[bf * M]);
    this._j(bf, slot, bf, slot, dl.dqds);
    this._j(bf, slot, bf, 0, side ? 1 : -1);
    chain(dl.dqdc);
    for (let r = 0; r <= n; r++) {
      if (w[r] === 0) continue;
      this._j(bf, slot, b, r, w[r]);
      this._j(b, 0, b, r, -Af * w[r]); // (and the box's charge: −Af·q in its Gauss row)
    }
    res[R[b * M]] -= Af * dl.q;
    this._j(b, 0, bf, slot, -Af * dl.dqds);
    // the box's store of each species it holds: Af (Γ_i − Γ_i,old)/dt
    if (!Number.isFinite(dt)) return;
    const k = Af / dt;
    for (let i = 0; i < n; i++) {
      if (!this.present[this.diffList[d].g * n + i] || !holds[i]) continue;
      res[R[b * M + 1 + i]] += k * (gam[d * n + i] - gamOld[d * n + i]);
      this._j(b, 1 + i, bf, slot, k * dl.dgds[i]);
      chain(dl.dgdc.subarray(i * n, i * n + n));
      for (let r = 0; r <= n; r++) if (w[r] !== 0) this._j(b, 1 + i, b, r, k * w[r]);
    }
  }

  // Layer d's edge node's ∂c_j/∂(slot r) for r = φ̂, η_1…η_n, [j·M + r] (src/regions.js, _dc).
  _diffuseChain(d) {
    const { n, M } = this, g = this.diffList[d].g, out = this.diffList[d].chain ?? (this.diffList[d].chain = new Float64Array(n * M));
    const t = this.dV2 ?? (this.dV2 = new Float64Array(M));
    for (let j = 0; j < n; j++) {
      if (!this.present[g * n + j]) {
        out.fill(0, j * M, (j + 1) * M);
        continue;
      }
      this._dc(g, j, t);
      out.set(t, j * M);
    }
    return out;
  }

  // The diffuse layers' charges and excesses at the present state (after the concentrations).
  _diffuseLayers() {
    const { n, u, c, z, gam } = this;
    this.diffList.forEach((d, k) => {
      diffuseLayer(c.subarray(d.g * n, d.g * n + n), z, n, u[d.o], d.units, d.dl, d.lat, d.draw);
      gam.set(d.dl.gam, k * n);
    });
  }

  // Full index of face f's gate q (its fraction open).
  _gateSlot(f, q) {
    return this.blockOfFace[f] * this.M + 1 + this.n + this.model.interfaces[f].reactions.length + q;
  }

  // The voltage a gate follows: across its face, or across its port's capacitance at its node.
  _gateVoltage(gt) {
    return gt.port === undefined ? this._faceVoltage(gt.f) : this._portVoltage(gt.port, gt.g);
  }

  // The voltage across port k's capacitance at node g, φ − (V − zeroCharge) (V): a membrane's,
  // inside (the region) minus outside.
  _portVoltage(k, g) {
    const { M, u, uLo, VT } = this, b = this.blockOfNode[g], port = this.model.ports[k];
    return VT * (u[b * M] + uLo[b * M]) - (this.termV[2 + k] - port.capacitance.zeroCharge);
  }

  // The voltage across face f, φ_right − φ_left (V): between its diffuse layers' planes, where it
  // has them (what the face's own layer, a membrane, has across it).
  _faceVoltage(f) {
    const { M, u, uLo, VT } = this, bf = this.blockOfFace[f], bL = bf - 1, bR = bf + 1, [dL, dR] = this.diffAt[f];
    const planes = (dR < 0 ? 0 : u[this.diffList[dR].o]) - (dL < 0 ? 0 : u[this.diffList[dL].o]);
    return VT * (u[bR * M] - u[bL * M] + (uLo[bR * M] - uLo[bL * M]) + planes);
  }

  // Face f's gates, each a row (x − x_old)/dt = α(1 − x) − βx in the voltage across the face.
  _faceGates(f, bf, bL, bR) {
    const { res, VT } = this, R = this.rix, itf = this.model.interfaces[f], dt = this.dtNow;
    const V = this._faceVoltage(f), k0 = 1 + this.n + itf.reactions.length;
    itf.gates.forEach((gate, q) => {
      const o = this._gateSlot(f, q), x = this.u[o];
      const [a, da] = gateRate(gate.alpha, V), [b, db] = gateRate(gate.beta, V);
      const store = Number.isFinite(dt) ? 1 / dt : 0, j = this.gateIndex[f] + q;
      res[R[o]] = store * (x - this.gateOld[j]) - (a * (1 - x) - b * x);
      this._j(bf, k0 + q, bf, k0 + q, store + a + b);
      const dV = -(da * (1 - x) - db * x) * VT; // per unit φ̂_R, and minus that per unit φ̂_L
      this._j(bf, k0 + q, bR, 0, dV);
      this._j(bf, k0 + q, bL, 0, -dV);
      // (and per unit of each diffuse layer's drop, as of its side's φ̂)
      const [dlL, dlR] = this.diffAt[f];
      if (dlR >= 0) this._j(bf, k0 + q, bf, this.diffList[dlR].slot, dV);
      if (dlL >= 0) this._j(bf, k0 + q, bf, this.diffList[dlL].slot, -dV);
    });
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
    // Beside a diffuse layer (src/diffuse.js), a participant reacts at the plane, where its
    // concentration is c e^{−z s}/D(s) (its μ̄ level through the layer; D = 1 unless the layer is
    // crowded): dS[side] = ∂ ln(prefactor)/∂s, and a crowded D's ∂/∂c_j into the side's slots.
    const dla = this.diffAt[f], dS = [0, 0];
    const layer = (p) => (dla[p.side] < 0 ? null : this.diffList[dla[p.side]]), atPlane = (p) => layer(p) && (this.z[p.i] !== 0 || layer(p).lat?.on[p.i]);
    // A drawn participant's level at the plane is w below its edge node's (_diffuseDraws): its c
    // there e^w times as much, its η w more; the face-block slot of that w, or −1.
    const wSlot = (p) => (atPlane(p) ? layer(p).drawSlot[p.i] : -1);
    const wOf = (p) => (wSlot(p) < 0 ? 0 : u[bf * M + wSlot(p)]);
    const plane = (p) => (atPlane(p) ? Math.exp(Math.max(-700, Math.min(700, -this.z[p.i] * u[layer(p).o] + wOf(p)))) / layer(p).dl.D : 1);
    const dW = this.dWr ?? (this.dWr = new Float64Array(M)), aW = this.aWr ?? (this.aWr = new Float64Array(M)); // ∂ ln(prefactor), ∂a by each face slot w
    dW.fill(0);
    aW.fill(0);
    // ∂ ln(e^{−z s}/D) by s, and w·(−∂ ln D/∂c_j) through the edge node's slots into d
    const planeSlope = (p) => -this.z[p.i] - layer(p).dl.dDds / layer(p).dl.D;
    const crowdInto = (p, w, d) => {
      const L = layer(p);
      if (!L.lat || w === 0) return;
      const dc = this._diffuseChain(dla[p.side]);
      for (let j = 0; j < n; j++) {
        const a = (-w * L.dl.dDdc[j]) / L.dl.D;
        if (a !== 0) for (let r = 0; r <= n; r++) d[r] += a * dc[j * M + r];
      }
    };
    let pref = rx.k0, aHi = rx.fixedA, aLo = 0, rate;
    for (const p of rx.part) {
      const g = p.side ? gR : gL, b = p.side ? bR : bL, o = b * M + 1 + p.i;
      aHi -= p.nu * (u[o] + wOf(p));
      aLo -= p.nu * uLo[o];
      if (wSlot(p) >= 0) aW[wSlot(p)] -= p.nu;
      if (rx.srh || rx.vmax || rx.mhc || this.nodeConductor[g] === p.i) continue; // a conductor's carrier has activity 1
      const e = p.nu < 0 ? -p.nu * (1 - al) : p.nu * al;
      pref *= powr((c[g * n + p.i] * plane(p)) / this.cRef[g * n + p.i], e);
      this._dlnc(g, p.i, e, p.side ? dR : dL);
      if (atPlane(p)) {
        dS[p.side] += e * planeSlope(p);
        crowdInto(p, e, p.side ? dR : dL);
        if (wSlot(p) >= 0) dW[wSlot(p)] += e;
      }
    }
    if (rx.vmax) {
      // Saturating: r = vmax Π (c/(c + K))^{|ν|} (1 − e^{−a}); ∂ln S/∂ln c = |ν| K/(c + K).
      let S = 1;
      rx.part.forEach((p, q) => {
        if (p.nu >= 0) return;
        const g = p.side ? gR : gL, cc = c[g * n + p.i] * plane(p), K = rx.K[q];
        S *= powi(cc / (cc + K), -p.nu);
        this._dlnc(g, p.i, (-p.nu * K) / (cc + K), p.side ? dR : dL);
        if (atPlane(p)) {
          dS[p.side] += ((-p.nu * K) / (cc + K)) * planeSlope(p);
          crowdInto(p, (-p.nu * K) / (cc + K), p.side ? dR : dL);
          if (wSlot(p) >= 0) dW[wSlot(p)] += (-p.nu * K) / (cc + K);
        }
      });
      const a = aHi + aLo, back = Math.exp(-a);
      rate = rx.vmax * S * -Math.expm1(-a);
      for (let s = 0; s < M; s++) {
        dL[s] *= rate;
        dR[s] *= rate;
        dW[s] = rate * dW[s] + rx.vmax * S * back * aW[s];
      }
      for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= rx.vmax * S * back * p.nu;
    } else if (rx.srh) {
      // SRH, with n and p each at its own side's edge node.
      const { n: pn, p: pp } = rx.srh, gn = pn.side ? gR : gL, gq = pp.side ? gR : gL;
      const s = this._srhRate(rx.srh, c[gn * n + pn.i], c[gq * n + pp.i], aHi + aLo, ...this._srhDegen(gn, gq, { n: pn, p: pp }));
      rate = s.rate;
      this._srhInto(gn, pn.i, s.Cn, s.Cxn, pn.side ? dR : dL);
      this._srhInto(gq, pp.i, s.Cp, s.Cxp, pp.side ? dR : dL);
      for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= s.Ca * p.nu;
    } else if (rx.mhc) {
      // Marcus–Hush–Chidsey: r = k0 P F(a°) (1 − e^{−a}), P = Π_{ν<0} (c/c_ref)^{|ν|}, Q = Π_{ν>0}
      // (c/c_ref)^ν (at the plane beside a diffuse layer), a° = a − ln P + ln Q, F = I(λ, a°)/I(λ, 0).
      let lnP = 0, lnQ = 0;
      for (const p of rx.part) {
        const g = p.side ? gR : gL;
        if (this.nodeConductor[g] === p.i) continue;
        const l = Math.log((c[g * n + p.i] * plane(p)) / this.cRef[g * n + p.i]);
        if (p.nu < 0) lnP -= p.nu * l;
        else lnQ += p.nu * l;
      }
      const a = aHi + aLo, { I, dI } = mhcIntegral(rx.mhc.lam, a - lnP + lnQ);
      const kP = rx.k0 * Math.exp(lnP), F = I / rx.mhc.I0, Fp = dI / rx.mhc.I0;
      const back = Math.exp(Math.min(-a, 700)), fwd = 1 - back; // (held finite where Newton's iterates stray)
      rate = kP * F * fwd;
      // ∂rate/∂a, ∂rate/∂ln P and ∂rate/∂ln Q, each participant's ln c weighted by its |ν| in one of them
      const ra = kP * (Fp * fwd + F * back), rq = kP * Fp * fwd, rp = rate - rq;
      for (const p of rx.part) {
        const g = p.side ? gR : gL;
        if (this.nodeConductor[g] === p.i) continue;
        const e = p.nu < 0 ? -p.nu * rp : p.nu * rq;
        this._dlnc(g, p.i, e, p.side ? dR : dL);
        if (atPlane(p)) {
          dS[p.side] += e * planeSlope(p);
          crowdInto(p, e, p.side ? dR : dL);
          if (wSlot(p) >= 0) dW[wSlot(p)] += e;
        }
      }
      for (let s = 0; s < M; s++) dW[s] += ra * aW[s];
      for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= ra * p.nu;
    } else {
      const { g: bv, gp } = bvFactor(aHi + aLo, al);
      rate = pref * bv;
      // ∂rate/∂slot = rate·∂ln(prefactor) + pref·g′·∂a, with ∂a/∂η_i = −ν
      for (let s = 0; s < M; s++) {
        dL[s] *= rate;
        dR[s] *= rate;
        dW[s] = rate * dW[s] + pref * gp * aW[s];
      }
      for (const p of rx.part) (p.side ? dR : dL)[1 + p.i] -= pref * gp * p.nu;
    }
    const slot = 1 + n + k, o = bf * M + slot;
    res[R[o]] = u[o] - rate;
    this._j(bf, slot, bf, slot, 1);
    // (dS is ∂ ln(prefactor)/∂s, except for Marcus–Hush–Chidsey's, which is ∂rate/∂s already)
    for (const sd of [0, 1]) if (dS[sd] !== 0) this._j(bf, slot, bf, this.diffList[dla[sd]].slot, -(rx.mhc ? 1 : rate) * dS[sd]);
    if (!rx.srh) for (let q = 0; q < M; q++) if (dW[q] !== 0) this._j(bf, slot, bf, q, -dW[q]);
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

  // The terminals whose ∂res/∂V (termB) and ∂I/∂x (termC) are kept, exactly: the floating ones,
  // which the Newton solve borders on, or all of them (allTerminals: for the impedance, which reads
  // its terminal's). A held terminal's are otherwise left as they fall (written but not kept
  // consistent), which spares a device with dozens of ports (an axon's nodes of Ranvier) the cost
  // of carrying each one's two dense vectors through every assembly.
  _liveTerms() {
    if (this.allTerminals) return this._allTerms ?? (this._allTerms = this.terms.map((_, k) => k));
    return this.floating;
  }

  // Every terminal's terms: its voltage held or floating, how the residual depends on it
  // (termB), its current into the device (termI) and how that depends on the state (termC).
  _terminals(dt) {
    this._refreshSources();
    for (const k of this._liveTerms()) {
      this.termB[k].fill(0);
      this.termC[k].fill(0);
    }
    this.termDI.fill(0);
    for (const x of this.termDX) x.fill(0);
    this.termI.fill(0);
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

  // into += w × (the Jacobian row of block b, slot s), over compact columns; and where a terminal
  // k reads that row as its current, its ∂/∂V's too: w × each live terminal's termB entry in it,
  // into termDI[k] (its own) or termDX[k][j].
  _captureRow(b, s, w, into, k = -1) {
    const l = this.loc[b * this.M + s];
    if (l < 0 || w === 0) return;
    if (k >= 0) {
      const o = this.rix[b * this.M + s];
      for (const j of this._liveTerms()) {
        const v = this.termB[j][o];
        if (v === 0) continue;
        if (j === k) this.termDI[k] += w * v;
        else this.termDX[k][j] += w * v;
      }
    }
    if (this.lin) this.lin.capture(this.rix[b * this.M + s], w, this.lin.rate, into);
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
      for (let w = 0; w < port.nodes.length; w++) {
        const g = port.nodes[w], b = this.blockOfNode[g], o = b * M + r, v = vol[g] * port.fill[w];
        const deta = target - u[o] - uLo[o]; // (μ̄_out − μ̄)/RT
        if (link.type === 'equilibrium') {
          if (!port.inside[w]) continue; // (a held level holds the nodes inside the window)
          // a device end node's level is its contact's business, where the contact links it
          if ((g === 0 && this.model.contacts.left.species[i].type !== 'blocked') || (g === this.nNodes - 1 && this.model.contacts.right.species[i].type !== 'blocked')) continue;
          flux[i] += res[R[o]];
          this._captureRow(b, r, zF, C, k);
          this._replaceRow(b, r);
          for (const Bk of this.termB) Bk[R[o]] = 0; // (other ports' terms in the row replaced)
          this._j(b, r, b, r, 1);
          res[R[o]] = -deta;
          B[R[o]] += -z[i] / VT;
          continue;
        }
        // conductance: s = G V_T (η_out − η)/(z² F); exchange: s = k (η_out − η). Gated, G times
        // Π x^p, the gates at this node (each taken in [0, 1]).
        const kk0 = link.type === 'conductance' ? (link.G * VT) / (z[i] * z[i] * F) : link.k;
        const gated = link.gates ?? [], q0 = 1 + n + port.surface.length;
        const open = gated.length ? this._open(gated, (q) => u[b * M + q0 + q]) : 1, kk = kk0 * open;
        res[R[o]] -= v * kk * deta;
        this._j(b, r, b, r, v * kk);
        flux[i] += v * kk * deta;
        if (z[i] !== 0) {
          B[R[o]] += (-v * kk * z[i]) / VT;
          C[R[o]] += -zF * v * kk;
          this.termDI[k] += (zF * v * kk * z[i]) / VT;
        }
        for (const [q, p] of gated) {
          const d = this._dOpen(gated, q, p, (qq) => u[b * M + q0 + qq]);
          if (d === 0) continue;
          this._j(b, r, b, q0 + q, -v * kk0 * deta * d);
          if (z[i] !== 0) C[R[b * M + q0 + q]] += zF * v * kk0 * deta * d;
        }
      }
    }
    if (holds) {
      for (let i = 0; i < n; i++) this.termI[k] += z[i] * F * flux[i];
      return;
    }
    port.reactions.forEach((rx, x) => this._portReaction(port, rx, this.portArea[k - 2], this.portRates[k - 2][x], flux, k));
    if (port.capacitance) this._portCapacitance(port, k, dt);
    if (port.gates.length > 0) this._portGates(port, k, dt);
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
          // a surface species at its coverage θ = θ₀e^ζ, with ∂ln θ_s/∂η_x = δ_sx − θ_x
          pref *= powr(th[g * nSurf + p.s], e);
          d[1 + n + p.s] += e;
          for (let x = 0; x < surf.length; x++) d[1 + n + x] -= e * th[g * nSurf + x];
        } else {
          pref *= powr(c[g * n + p.i] / this.cRef[g * n + p.i], e);
          this._dlnc(g, p.i, e, d);
        }
      }
      if (rx.sites) {
        // the bare sites its transition state holds: θ₀ = 1 − Σθ, with ∂ln θ₀/∂η_s = −θ_s
        for (let x = 0; x < surf.length; x++) d[1 + n + x] -= rx.sites * th[g * nSurf + x];
        pref *= powr(this.th0[g], rx.sites);
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

  // A membrane port's gates, at each node of its window: (x − x_old)/dt = α(1 − x) − βx in the
  // voltage across its capacitance there, φ − (V − zeroCharge).
  _portGates(port, k, dt) {
    const { n, res, VT } = this, R = this.rix, B = this.termB[k], store = Number.isFinite(dt) ? 1 / dt : 0, q0 = 1 + n + port.surface.length;
    port.nodes.forEach((g, w) => {
      const b = this.blockOfNode[g], V = this._portVoltage(k - 2, g);
      port.gates.forEach((gate, q) => {
        const o = b * this.M + q0 + q, x = this.u[o], j = this.portGateIndex[k - 2] + w * port.gates.length + q;
        const [a, da] = gateRate(gate.alpha, V), [bb, db] = gateRate(gate.beta, V);
        res[R[o]] = store * (x - this.gateOld[j]) - (a * (1 - x) - bb * x);
        this._j(b, q0 + q, b, q0 + q, store + a + bb);
        const s = -(da * (1 - x) - db * x); // ∂/∂V_m, with V_m = V_T φ̂ − V + zeroCharge
        this._j(b, q0 + q, b, 0, s * VT);
        B[R[o]] += -s;
      });
    });
  }

  // Π x^p over a link's gates (each x taken in [0, 1]), and its derivative in gate q.
  _open(gated, x) {
    let open = 1;
    for (const [q, p] of gated) open *= powi(Math.min(1, Math.max(0, x(q))), p);
    return open;
  }

  _dOpen(gated, q, p, x) {
    const xq = x(q);
    if (!(xq > 0 && xq < 1)) return 0;
    let others = 1;
    for (const [q2, p2] of gated) if (q2 !== q) others *= powi(Math.min(1, Math.max(0, x(q2))), p2);
    return p * powi(xq, p - 1) * others;
  }

  // A capacitance spread through a port's window: per area of electrode, σ = C (V − zeroCharge − φ)
  // on the port's side, so each node's charge balance (its φ row) gains vol·a·σ, and the port
  // passes the charging current d(Σ vol·a·σ)/dt (none in a steady state). An end node whose φ its
  // contact sets is left to the contact.
  _portCapacitance(port, k, dt) {
    const { M, u, uLo, res, VT } = this, R = this.rix, vol = this.model.grid.vol, area = this.portArea[k - 2];
    const { C: Cs, zeroCharge } = port.capacitance, V = this.termV[k], B = this.termB[k], Ct = this.termC[k];
    const sets = (side) => ['pinned', 'bulk'].includes(this.model.contacts[side].phi.type);
    // Its current, the charging (Q − Q at the step's start)/dt; or, steady at no current, the
    // charge it keeps in place of that (see captureLaws in src/lawRows.js).
    const rate = Number.isFinite(dt) ? 1 / dt : this.steady && this.chargeRefs.has(k) ? 1 : 0;
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
      if (rate > 0) Ct[R[b * M]] += -s * VT * rate;
    });
    this.portQ[k - 2] = Q;
    if (rate > 0) {
      this.termI[k] += (Q - (dt < Infinity ? this.portQStart[k - 2] : this.chargeRefs.get(k))) * rate;
      this.termDI[k] += dQdV * rate;
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
    const held = this.portHeldEnd.subarray(k * n, (k + 1) * n);
    for (let i = 0; i < n; i++) {
      const active = this.loc[b * M + 1 + i] >= 0 && !(closedIons && z[i] !== 0) && !held[i];
      const nIn = active ? res[R[b * M + 1 + i]] : 0;
      flux[i] = sgn * nIn; // (toward +x)
      if (active && z[i] !== 0) {
        I += z[i] * F * nIn;
        this._captureRow(b, 1 + i, z[i] * F, C, k);
      }
    }
    // The ions' storage, left out until now (see lateStorage): in each balance, and in what
    // each species brings in, but not in the current.
    if (this.combining && dyn && this.lateStorage[g] === 1) {
      const v = model.grid.vol[g], { cOld } = this;
      for (let i = 0; i < n; i++) {
        if (z[i] === 0 || this.loc[b * M + 1 + i] < 0 || held[i]) continue;
        const ci = c[g * n + i], st = (v * (ci - cOld[g * n + i])) / dt;
        res[R[b * M + 1 + i]] += st;
        flux[i] += sgn * st;
        this._j(b, 1 + i, b, 1 + i, (v * ci) / dt);
        this._j(b, 1 + i, b, 0, (-v * z[i] * ci) / dt);
      }
    }

    // Conductance links: J (toward the device) = G (V_out − V_i), V_out = V_t + offset; exchange
    // links (neutral species): N_in = k (μ_out − μ)/RT; velocity links, v (c_eq − c).
    for (let i = 0; i < n; i++) {
      const link = ct.species[i];
      const o = b * M + 1 + i;
      if (link.type === 'exchange') {
        res[R[o]] -= Ac * link.k * (link.mu / model.RT - u[o] - uLo[o]);
        this._j(b, 1 + i, b, 1 + i, Ac * link.k);
      } else if (link.type === 'conductance') {
        const Vi = (VT * (u[o] + uLo[o])) / z[i];
        res[R[o]] -= (Ac * link.G * (Vt + link.offset - Vi)) / (z[i] * F);
        this._j(b, 1 + i, b, 1 + i, (Ac * link.G * VT) / (z[i] * z[i] * F));
        B[R[o]] += (-Ac * link.G) / (z[i] * F);
      } else if (link.type === 'velocity') {
        // N_in = v c (e^(η_out − η) − 1): v (c_eq − c), c_eq in equilibrium with the outside.
        const ci = c[g * n + i], E = Math.exp(Math.min(this.contactEta(side, i) - u[o] - uLo[o], 700));
        const N = Ac * link.v * ci * (E - 1), d = this.dA;
        res[R[o]] -= N;
        d.fill(0);
        this._dlnc(g, i, 1, d);
        for (let s = 0; s <= n; s++) if (d[s] !== 0) this._j(b, 1 + i, b, s, -N * d[s]);
        this._j(b, 1 + i, b, 1 + i, Ac * link.v * ci * E);
        if (z[i] !== 0) B[R[o]] += (-Ac * link.v * ci * E * z[i]) / VT;
      }
    }

    // Equilibrium links: Dirichlet on the known outside level.
    for (let i = 0; i < n; i++) {
      const link = ct.species[i];
      if (link.type !== 'equilibrium') continue;
      const o = b * M + 1 + i;
      this._replaceRow(b, 1 + i);
      for (const Bk of this.termB) Bk[R[o]] = 0; // (a port's terms in the row replaced)
      this._j(b, 1 + i, b, 1 + i, 1);
      res[R[o]] = u[o] - this.contactEta(side, i) + uLo[o]; // (the big terms first: u + uLo rounds to u)
      if (z[i] !== 0) B[R[o]] += -z[i] / VT;
    }

    // φ law (none where φ is undefined: nothing at the end node responds to it). The displacement
    // into the device, D_in = sgn·D, adds (D_in − D_in at the step's start)/dt to the current.
    const link = ct.phi, Dstart = sgn * this.contactDStart[side];
    // (In a steady solve, what a law in this circuit row counts of D_in: see applyConstraints in
    // src/lawRows.js.)
    const q = !dyn && this.steady && this.contactLaw[k] ? (this.contactQ[k] ??= { row: new Float64Array(this.sys.size + 1), D: 0, dV: 0 }) : null;
    if (q) {
      q.row.fill(0);
      q.D = 0;
      q.dV = 0;
    }
    if (this.phiUndefined[g]) {
      this.contactD[side] = 0;
    } else if (link.type === 'capacitive') {
      // Gate, or metal across a Stern layer, at φ_g = V_t − zeroCharge. D toward +x.
      const Din = Ac * link.C * (Vt - link.zeroCharge - VT * u[b * M]);
      res[R[b * M]] -= Din;
      this._j(b, 0, b, 0, Ac * link.C * VT);
      B[R[b * M]] += -Ac * link.C;
      this.contactD[side] = sgn * Din;
      if (q && this.loc[b * M] >= 0) [q.D, q.row[R[b * M]], q.dV] = [Din, -Ac * link.C * VT, Ac * link.C];
      // Its current, the charging (D_in − D_in at the step's start)/dt; or, steady at no current,
      // the charge it keeps in place of that (see captureLaws in src/lawRows.js).
      const rate = dyn ? 1 / dt : this.steady && this.chargeRefs.has(k) ? 1 : 0;
      if (rate > 0) {
        I += (Din - (dyn ? Dstart : this.chargeRefs.get(k))) * rate;
        C[R[b * M]] += -Ac * link.C * VT * rate;
        this.termDI[k] += Ac * link.C * rate;
      }
    } else if (link.type === 'pinned' || link.type === 'bulk') {
      // The Poisson residual is the outside's charge, D_in.
      const Din = res[R[b * M]];
      this.contactD[side] = sgn * Din;
      if (q) {
        q.D = Din;
        this._captureRow(b, 0, 1, q.row);
      }
      if (dyn && !(this.combining && this.lateStorage[g] === 1)) {
        I += (Din - Dstart) / dt;
        this._captureRow(b, 0, 1 / dt, C, k);
      } else if (!dyn && this.steady && this.chargeRefs.has(k)) {
        I += Din - this.chargeRefs.get(k);
        this._captureRow(b, 0, 1, C, k);
      }
      this._replaceRow(b, 0);
      for (const Bk of this.termB) Bk[R[b * M]] = 0;
      if (link.type === 'pinned') {
        // φ_edge = V_t − zeroCharge (the C → ∞ limit).
        this._j(b, 0, b, 0, 1);
        res[R[b * M]] = u[b * M] - (Vt - link.zeroCharge) / VT + uLo[b * M];
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
}
