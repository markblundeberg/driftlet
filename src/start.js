// The cold start: the state a device starts from, built from its definition (each region's c0,
// the contacts' and ports' levels) before any solve or step.
//
// CONTRACT
//
// initFromComposition(s) writes a Solver `s` (src/solver.js; its layout is in docs/internals.md)
// once constructed: its unknowns u (every slot: φ̂, each present species' η, an electrode
// surface's coverages as η, gates' fractions) with uLo = 0, the voltages termV of its floating
// terminals, and, through s.computeConcentrations(), its concentrations and coverages. It reads
// the model (regions' c0 and c0Profile, materials, interfaces, contacts, ports), the stretches
// (s.stretches, s.stretchOf), the terminals (s.terms, s.floating) and the held terminals'
// voltages, and calls s.contactEta, s.portEta, s._surfaceStart, s._portRate, s._gateVoltage.
//
// In order:
//   1. Each region, left to right, each species present in it:
//      - its c0 (or c0Profile) where given: an amount, mode 2, the starting state even where a
//        contact feeds the species;
//      - else a conductor's carrier away from the contacts: at equilibrium with the first
//        reaction on its left face that takes it, else at the running φ̂ (level, mode 1);
//      - else not fed by a contact: c0 is required (an error without it);
//      - else at the level of what feeds it: the left contact's where its stretch reaches the
//        left end, else the right's, else its first linked port's; fed by none of them
//        (reactions only), from c0, required.
//      φ̂ runs from the left: a pinned or capacitive left contact's (V − zeroCharge)/V_T, else 0,
//      plus each face's dipole. In a region where a level-fixed charged species responds to φ̂,
//      φ̂ is chosen for neutrality (a bisection; node by node where a profile varies);
//      otherwise it's the running φ̂, or where a conductor on the left or an electrode spread
//      through the region has a reaction that can balance, the φ̂ that balances it. A strictly
//      neutral region (ε = 0) with every charged species' amount given must have a neutral c0:
//      an error otherwise. A non-ideal material's composition comes from its statistics
//      (materialAt).
//   2. A face between two strictly neutral ideal regions, its φ type neutral, whose links hold
//      some charged species in equilibrium: its two nodes are set to that equilibrium (each
//      passing species' μ̄ equal across it, the left node neutral, each passing species' amount
//      over the two nodes kept; the jump in φ̂ between them the one unknown), everything else at
//      them unchanged in concentration.
//   3. An insertion host (a φ-free material): its φ̂ and every η shifted by z·Δ so that its first
//      charged species shared with a neighbour across a face starts at that neighbour's level.
//   4. A floating contact facing a held one: V = the held V + the difference between the two
//      terminal species' levels beside each (where both are equilibrium links of charged
//      species).
//   5. Electrode surfaces at their θ₀ (s._surfaceStart); concentrations computed.
//   6. Floating double layers start uncharged: a capacitive contact at zeroCharge + V_T·φ̂
//      beside it; a port with a capacitance at zeroCharge + its window's mean V_T·φ̂ (by
//      volume·area, where φ is defined).
//   7. A floating port linked to its carrier, with no reactions or capacitance: at the carrier's
//      level across its window, relative to a held contact's (where there is one).
//   8. Gates at α/(α + β) at the starting voltage across them.
//   9. A floating port with reactions and no capacitance: at the voltage where its reactions
//      pass the current it's driven at (0 behind a resistance), by bisection.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. A species given by its amount starts at c0 (or its profile) at every node, but for the
//      two nodes of a face equilibrated in step 2, where its amount over the two is kept.
//   2. A species given by its level starts at the level of what feeds it, as step 1 orders them.
//   3. Where a level-fixed charged species responds, each node starts neutral (Σ z c + the fixed
//      charge = 0) to the bisection's tolerance; a strictly neutral region given by amounts is
//      neutral by its c0, checked.
//   4. Steps 2 and 3 change no concentration but as they say: step 3 is a gauge (φ̂ and each η by
//      z·Δ), step 2 moves the passing species between a face's two nodes only.
//   5. Each floating terminal's V is set by the last of steps 4, 6, 7, 9 that applies to it;
//      held terminals' V are untouched.
//   6. Only u, uLo, termV and what computeConcentrations derives are written; the result
//      depends on the model and the held voltages only, deterministically.
// Known, and not an invariant: the gauge (φ̂ and each η by z·Δ) runs from the left, so behind a
// closed or floating left end the start sits away from a right contact's level (55 RT for a
// bath referenced to Cl⁻); see ROADMAP.
// Tests: every test (each device starts here); test/transient.test.js on what a start must be.

import { FARADAY } from './constants.js';
import { sourceAt } from './device.js';
import { gateRate, profileAt, SolverError } from './solver.js';

/**
 * Cold start. Each species takes its region's c0 where one is given (for a spectator, it fixes
 * the conserved amount), else its contact's level. φ in each region is then chosen for local neutrality (node by node where a
 * profile varies), or continued across the interface dipole if nothing there responds.
 */
export function initFromComposition(s) {
  const { model, n, M, u, z } = s;
  const { regions, species, interfaces, materials } = model;
  const grid = model.grid;
  u.fill(0);
  s.uLo.fill(0);
  const eta = new Float64Array(n), cFix = new Float64Array(n), mode = new Int8Array(n); // 1 level, 2 amount
  // φ̂ runs from the left: from a pinned (or gated) left contact's φ at its zero charge, else 0.
  const leftPhi = model.contacts.left.phi;
  let phiHat = leftPhi.type === 'pinned' || leftPhi.type === 'capacitive' ? (s.termV[0] - leftPhi.zeroCharge) / s.VT : 0;
  for (let r = 0; r < regions.length; r++) {
    if (r > 0) phiHat += interfaces[r - 1].dipole / s.VT;
    const reg = regions[r], mat = materials[reg.material];
    mode.fill(0);
    const profiled = reg.c0Profile.some((p) => p !== null);
    for (let i = 0; i < n; i++) {
      if (!mat.present[i]) continue;
      const st = s.stretches[s.stretchOf[r * n + i]];
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
          const bn = s.blockOfNode[grid.regionEnd[r - 1]];
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
        const left = st.regions[0] === 0 ? s.contactEta('left', i) : NaN;
        eta[i] = Number.isFinite(left) ? left : s.contactEta('right', i);
        if (!Number.isFinite(eta[i]) && st.linked.length > 0) eta[i] = s.portEta(model.ports[st.linked[0]], i);
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
        materialAt(s, mat, reg.background, eta, mode, cFix, ph, zeta, cc);
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
      const bn = s.blockOfNode[grid.regionEnd[r - 1]];
      for (const rx of interfaces[r - 1].reactions) {
        let a = rx.fixedA, slope = 0;
        for (const p of rx.part) {
          if (p.side === 0) a -= p.nu * u[bn * M + 1 + p.i];
          else if (mode[p.i] === 2) {
            a -= p.nu * (Math.log(cFix[p.i] / mat.cRef[p.i]) + mat.mu0[p.i] / model.RT);
            slope += p.nu * z[p.i];
          } else slope = NaN;
        }
        if (slope !== 0 && Number.isFinite(slope)) {
          phiHat = a / slope; // a − slope·φ̂ = 0
          break;
        }
      }
    }
    const electrode = model.ports.find((port) => port.region === r && port.reactions.length > 0);
    if (!mat.phiFree && !responds && mat.ideal && electrode) {
      // Nor here, but an electrode spread through it does: start with the port's first reaction
      // that can be balanced at equilibrium with the port's level (the electrode at open circuit).
      const bare = 1 - electrode.surface.reduce((t, sp) => t + sp.theta0, 0);
      for (const rx of electrode.reactions) {
        let a = rx.fixedA, slope = 0;
        for (const p of rx.part) {
          if (p.side === 1) a -= p.nu * s.portEta(electrode, p.i);
          else if (p.side === 2) {
            const sp = electrode.surface[p.s]; // (at its starting coverage)
            a -= p.nu * (sp.mu0 / model.RT + Math.log(sp.theta0 / bare));
          } else if (mode[p.i] === 2) {
            a -= p.nu * (Math.log(cFix[p.i] / mat.cRef[p.i]) + mat.mu0[p.i] / model.RT);
            slope += p.nu * z[p.i];
          } else slope = NaN;
        }
        if (slope !== 0 && Number.isFinite(slope)) {
          phiHat = a / slope;
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
    if (!mat.ideal) materialAt(s, mat, reg.background, eta, mode, cFix, phiHat, zeta, cc);
    for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
      if (profiled) {
        for (let i = 0; i < n; i++) if (reg.c0Profile[i]) cFix[i] = profileAt(reg.c0Profile[i], grid.x[g]);
        if (mustBalance) checkNeutral(` at x = ${grid.x[g].toPrecision(4)} m`);
        if (!mat.phiFree && responds) phiHat = neutralPhi(phiHat);
        if (!mat.ideal) materialAt(s, mat, reg.background, eta, mode, cFix, phiHat, zeta, cc);
      }
      const b = s.blockOfNode[g];
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
  // Two strictly neutral regions meeting at a face whose links hold ions in equilibrium: where
  // c0 leaves the two sides out of it (a different fixed charge or immobile ion on each, so
  // c₊c₋ differs across it), the face reaches it at t = 0⁺, before anything diffuses, and a
  // first step made to take that jump met an exactly singular pivot. Start there: the ions the
  // face passes shared between its two nodes in equilibrium (μ̄ equal across it), each one's
  // amount over the two nodes kept, everything else as c0 has it. The one unknown is the jump
  // in φ̂ across it (a Donnan potential), set by neutrality on one side; the other's follows,
  // the two nodes' charge summing to what it was, zero.
  for (let f = 0; f < interfaces.length; f++) {
    const itf = interfaces[f], mL = materials[regions[f].material], mR = materials[regions[f + 1].material];
    const strict = (m) => m.epsr === 0 && !m.conductor && !m.phiFree && m.ideal;
    if (itf.phi.type !== 'neutral' || !strict(mL) || !strict(mR)) continue;
    const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1], bL = s.blockOfNode[gL] * M, bR = s.blockOfNode[gR] * M;
    const pass = [];
    for (let i = 0; i < n; i++) if (itf.links[i].type === 'equilibrium' && mL.present[i] && mR.present[i]) pass.push(i);
    if (!pass.some((i) => z[i] !== 0)) continue;
    const RT = model.RT, vL = grid.vol[gL], vR = grid.vol[gR];
    const conc = (m, b, i) => m.cRef[i] * Math.exp(u[b + 1 + i] - m.mu0[i] / RT - z[i] * u[b]);
    const amount = pass.map((i) => vL * conc(mL, bL, i) + vR * conc(mR, bR, i));
    let qOther = regions[f].fixedCharge / FARADAY; // (the left node's charge from all it doesn't share)
    for (let i = 0; i < n; i++) if (mL.present[i] && !pass.includes(i) && z[i] !== 0) qOther += z[i] * conc(mL, bL, i);
    // At a jump D = φ̂_R − φ̂_L: c_L/c_R = (cRef_L/cRef_R)·exp((μ°_R − μ°_L)/RT + z D), the amount fixing both.
    const split = (D) =>
      pass.map((i, j) => {
        const K = (mL.cRef[i] / mR.cRef[i]) * Math.exp((mR.mu0[i] - mL.mu0[i]) / RT + z[i] * D), cR = amount[j] / (vL * K + vR);
        return [K * cR, cR];
      });
    const qL = (D) => split(D).reduce((q, [cL], j) => q + z[pass[j]] * cL, qOther); // (increasing in D)
    const D0 = u[bR] - u[bL];
    let lo = D0 - 1, hi = D0 + 1;
    while (qL(lo) > 0 && lo > D0 - 1e3) lo -= 2 * (hi - lo);
    while (qL(hi) < 0 && hi < D0 + 1e3) hi += 2 * (hi - lo);
    if (!(qL(lo) <= 0 && qL(hi) >= 0)) continue; // (nothing the face shares can balance it: leave c0)
    for (let it = 0; it < 200 && hi - lo > 1e-13 * Math.max(1, Math.abs(D0)); it++) {
      const m = 0.5 * (lo + hi);
      if (qL(m) < 0) lo = m;
      else hi = m;
    }
    const D = 0.5 * (lo + hi), cs = split(D);
    const keep = Float64Array.from({ length: n }, (_, i) => (mR.present[i] ? conc(mR, bR, i) : 0)); // (the right node's, before φ̂ moves)
    u[bR] = u[bL] + D;
    for (let i = 0; i < n; i++) if (mR.present[i] && !pass.includes(i)) u[bR + 1 + i] = Math.log(keep[i] / mR.cRef[i]) + mR.mu0[i] / RT + z[i] * u[bR];
    pass.forEach((i, j) => {
      u[bL + 1 + i] = Math.log(cs[j][0] / mL.cRef[i]) + mL.mu0[i] / RT + z[i] * u[bL];
      u[bR + 1 + i] = Math.log(cs[j][1] / mR.cRef[i]) + mR.mu0[i] / RT + z[i] * u[bR];
    });
  }
  // An insertion host's φ is a gauge: c0 sets only its neutral combination's level, not how that
  // splits between its ion and its electrons. Start its ion level with the same ion's beside it
  // across a face (the equilibrium a link or a transfer reaction heads for), so its electrons,
  // and a floating contact on them, start where the cell's open-circuit voltage puts them, not
  // volts away. (A shift of φ̂ and each η by z·Δ: the concentrations don't change.)
  regions.forEach((reg, r) => {
    const mat = materials[reg.material];
    if (!mat.phiFree) return;
    const sides = [[r - 1, grid.regionStart[r], r > 0 ? grid.regionEnd[r - 1] : -1], [r + 1, grid.regionEnd[r], r + 1 < regions.length ? grid.regionStart[r + 1] : -1]];
    for (const [nr, gIn, gOut] of sides) {
      if (gOut < 0) continue;
      const other = materials[regions[nr].material];
      if (other.phiFree || other.conductor) continue;
      const i = z.findIndex((zi, j) => zi !== 0 && mat.present[j] && other.present[j]);
      if (i < 0) continue;
      const d = (u[s.blockOfNode[gOut] * M + 1 + i] - u[s.blockOfNode[gIn] * M + 1 + i]) / z[i];
      if (!Number.isFinite(d)) continue;
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const b = s.blockOfNode[g];
        u[b * M] += d;
        for (let j = 0; j < n; j++) if (mat.present[j]) u[b * M + 1 + j] += z[j] * d;
      }
      return;
    }
  });
  // A floating contact facing a held one starts where the start's own composition puts it: the
  // held voltage plus the difference between the two terminal species' levels beside each
  // contact (c0 sets them, and φ was chosen for neutrality, not to match either contact), so
  // that a transient's first reading means something.
  const level = (side) => {
    const ct = model.contacts[side], i = ct.terminal;
    if (i === null || ct.species[i].type !== 'equilibrium' || z[i] === 0) return NaN;
    const b = s.blockOfNode[side === 'left' ? 0 : grid.nNodes - 1];
    return (s.VT * u[b * M + 1 + i]) / z[i] - ct.species[i].offset;
  };
  for (const k of s.floating) {
    const t = s.terms[k], held = k === 0 ? 1 : 0;
    if (t.kind === 'port' || held >= s.terms.length || s.floating.includes(held)) continue;
    const shift = level(t.side) - level(s.terms[held].side);
    if (Number.isFinite(shift)) s.termV[k] = s.termV[held] + shift;
  }
  for (let g = 0; g < s.nNodes; g++) if (s.surfPort[g] >= 0) s._surfaceStart(g);
  s.computeConcentrations();
  // A floating capacitance starts uncharged (σ = 0 at the window's mean φ), as the starting
  // composition, neutral without it, assumes; with reactions too, its double layer then charges
  // toward their mixed potential.
  for (const k of s.floating) {
    const t = s.terms[k], port = t.kind === 'port' ? model.ports[t.index] : null;
    const gate = t.kind === 'contact' ? model.contacts[t.side].phi : null;
    if (gate?.type === 'capacitive') {
      // A gate contact likewise: at its zero-charge level from the φ beside it (a pinned
      // one's charge is the field's, already zero against a neutral start).
      const g = t.side === 'left' ? 0 : s.nNodes - 1;
      if (!s.phiUndefined[g]) s.termV[k] = gate.zeroCharge + s.VT * u[s.blockOfNode[g] * M];
      continue;
    }
    if (!port?.capacitance) continue;
    let w = 0, sum = 0;
    port.nodes.forEach((g, j) => {
      if (s.phiUndefined[g]) return;
      const a = s.model.grid.vol[g] * s.portArea[t.index][j];
      w += a;
      sum += a * s.VT * u[s.blockOfNode[g] * M];
    });
    if (w > 0) s.termV[k] = port.capacitance.zeroCharge + sum / w;
  }
  // A floating port linked to its carrier (no reactions, no double layer: a wire into an
  // insertion host) starts at that carrier's level across its window, not at 0 V.
  for (const k of s.floating) {
    const t = s.terms[k], port = t.kind === 'port' ? model.ports[t.index] : null;
    const i = port?.terminal, link = i === null || i === undefined ? null : port.species[i];
    if (!port || port.reactions.length > 0 || port.capacitance || !link || link.type === 'blocked' || z[i] === 0) continue;
    let w = 0, sum = 0;
    port.nodes.forEach((g, j) => {
      const a = s.model.grid.vol[g] * (s.portArea[t.index]?.[j] ?? 1); // (a window with no area: by volume)
      w += a;
      sum += (a * s.VT * u[s.blockOfNode[g] * M + 1 + i]) / z[i];
    });
    // (Levels at the start mean something only against a held terminal's, as for a contact.)
    const held = [0, 1].find((h) => !s.floating.includes(h) && Number.isFinite(level(s.terms[h].side)));
    if (w > 0 && held !== undefined) s.termV[k] = s.termV[held] + sum / w - (link.offset ?? 0) - level(s.terms[held].side);
  }
  // Gates start open as far as the starting voltage across their face holds them, α/(α + β).
  for (const gt of s.gateList) {
    const { gate, o } = gt, V = s._gateVoltage(gt), a = gateRate(gate.alpha, V)[0], b = gateRate(gate.beta, V)[0];
    u[o] = a + b > 0 ? a / (a + b) : 0;
  }
  // A floating electrode spread through a port (without a double layer) starts where its
  // reactions pass the current it's set (none, behind a resistance): at its mixed potential in
  // the start's composition, not level with a held terminal, which can be volts away and pass
  // an absurd current.
  for (const k of s.floating) {
    const t = s.terms[k], port = t.kind === 'port' ? model.ports[t.index] : null;
    if (!port || port.reactions.length === 0 || port.capacitance) continue;
    const want = t.drive.kind === 'I' ? sourceAt(t.drive.src, s.time) : 0;
    const current = (V) => {
      s.termV[k] = V;
      let I = 0;
      port.reactions.forEach((rx) => {
        let q = 0;
        for (const p of rx.part) if (p.side === 0) q += p.nu * z[p.i];
        port.nodes.forEach((g, w) => (I += q * FARADAY * s.model.grid.vol[g] * s.portArea[t.index][w] * s._portRate(port, rx, g)));
      });
      return I;
    };
    // Current into the device rises with V (oxidation); bracket, then bisect.
    let lo = s.termV[k] - 1, hi = s.termV[k] + 1;
    while (current(lo) > want && lo > -1e3) lo -= 2 * (hi - lo);
    while (current(hi) < want && hi < 1e3) hi += 2 * (hi - lo);
    for (let it = 0; it < 200 && hi - lo > 1e-12; it++) {
      const m = 0.5 * (lo + hi);
      if (current(m) < want) lo = m;
      else hi = m;
    }
    s.termV[k] = 0.5 * (lo + hi);
  }
}

// A non-ideal material's composition at φ̂ = ph: level-fixed species (mode 1) at their η,
// amount-fixed ones (mode 2) at cFix. Fills ζ and c for every present species.
function materialAt(s, mat, background, eta, mode, cFix, ph, zeta, cc) {
  const { n, z } = s;
  const RT = s.model.RT;
  for (let i = 0; i < n; i++) {
    if (!mat.present[i]) continue;
    zeta[i] = mode[i] === 2 ? Math.log(cFix[i] / mat.cRef[i]) : eta[i] - mat.mu0[i] / RT - z[i] * ph;
    if (mat.modelOf[i] < 0) cc[i] = mode[i] === 2 ? cFix[i] : mat.cRef[i] * Math.exp(Math.min(700, zeta[i]));
  }
  for (const md of mat.models) {
    const w = s._work(md), idx = md.idx, k = idx.length;
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
