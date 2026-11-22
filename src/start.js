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
// (s.stretches, s.stretchOf), the terminals (s.terms, s.floating), every terminal's voltage as
// it stands (held or floating: a floating one's sets the levels it feeds) and s.time (the
// drives' sources), and calls s.contactEta, s.portEta, s._surfaceStart, s._portRate,
// s._gateVoltage, s._open.
//
// In order, but for step 2, which is last (the state is marked cold after steps 3–9, and
// equilibrateFaces runs from the first time step, after all of them; a code-smell review):
//   1. Each region, left to right, each species present in it:
//      - its c0 (or c0Profile) where given: an amount, mode 2, the starting state even where a
//        contact feeds the species;
//      - else a conductor's carrier away from the contacts: at equilibrium with the first
//        reaction on its left face that takes it, else at the running φ̂ (level, mode 1);
//      - else not fed by a contact: c0 is required (an error without it);
//      - else at the level of what feeds it: a contact's where its stretch reaches that end, the
//        left one's where it reaches both, but the right one's there where the right contact
//        alone has φ 'bulk'; else a linked port's, a held one before a floating one; a contact's link held or exchanged at
//        a rate (conductance, exchange) alike; fed by none of them (reactions only), from c0,
//        required.
//      φ̂ runs from the left: a held (not floating) pinned or capacitive left contact's
//      (V − zeroCharge)/V_T, else 0,
//      plus each face's dipole. In a region where a level-fixed charged species responds to φ̂,
//      φ̂ is chosen for neutrality (a bisection; node by node where a profile varies);
//      otherwise it's the running φ̂, or where a conductor on the left or an electrode spread
//      through the region has a reaction that can balance, the φ̂ that balances it, each species
//      in it at its activity where the electrode is (a c0 profile's at the region's end for a
//      metal on its right, at the window's middle for a port; a non-ideal material's from its
//      statistics) and each surface species at its starting coverage (a face's or the port's).
//      A conductor on the right counts where the conductor on the left balances nothing: in the
//      first region, or further in where nothing has set φ̂ yet (a pinned or capacitive left
//      contact, a region's neutrality or capacitance, an electrode's balance), and then the
//      regions before it that took the running φ̂ are moved with it (φ̂ and each η by z·Δ); so,
//      likewise, where a held electrode spread through a region is the first to set it. Of a
//      region's electrodes, a held one counts before a floating one, whatever order lists them;
//      a floating port's electrode counts only where nothing has set φ̂ (its V, found in step 9
//      or 6, is not known yet). A pinned, held contact whose own density, cRef e^(η − μ°/RT −
//      z φ̂_edge), is beyond e^700 is refused (SolverError).
//      In a strictly neutral region, a held port's capacitance through its window counts in
//      each node's neutrality there, ρ + aσ = 0 with σ = C (V − zeroCharge − V_T φ̂) as the
//      solve has it (not at an end box whose φ its contact sets), so φ̂ there responds even
//      where no level-fixed species does, and is chosen node by node (the neutrality is by
//      bisection; a trial a non-ideal model refuses counts its ideal charge, which has the sign
//      there); where nothing else in the region responds, the nodes outside the window take
//      the window's φ̂, as does what follows. Where
//      neutrality is out of reach (c0's amounts leave a charge no held level balances): in a
//      strictly neutral region the running φ̂, the region noted in s.unbalanced (for a transient
//      that fails from it to say why); otherwise the φ̂ nearest the running one, toward less
//      charge, where the charge is at most ε V_T/(F L²) (L the region's length: what it holds with
//      a thermal voltage across it), or where it stops falling (a fixed charge). A strictly
//      neutral region (ε = 0) with every charged species' amount given must have a neutral c0
//      outside such a window: an error otherwise. A non-ideal material's composition comes from its statistics
//      (materialAt).
//   2. Marks the state cold (s.coldStart). equilibrateFaces(s), which the first time step from it
//      calls (src/timeSteps.js) and a steady solve doesn't: a face between two strictly neutral
//      ideal regions, its φ type neutral, whose links hold some charged species in equilibrium,
//      its two nodes set to that equilibrium, its t = 0⁺ (each passing species' μ̄ equal across
//      it, each node neutral with a held gate's aσ counted in it, each passing species' amount
//      over the two nodes kept, and so the gates' charge over them; the jump in φ̂ between them
//      the one unknown, the left node's φ̂ moving with it only where the right node is gated),
//      everything else at them unchanged in concentration; a floating gate over either node
//      keeping its σ, its voltage and the rest of its window's φ̂ (concentrations kept) moving
//      with the φ̂ moved under it (floatSigma, followFloating);
//      s.coldStart false. A solve that converges, a warm start and a carried state set it as
//      their state is.
//   3. An insertion host (a φ-free material, not a conductor): its φ̂ and every η shifted by z·Δ
//      so that its first charged species shared with a neighbour that the face between them
//      passes (by a link or a reaction taking it across; else its first shared) starts at that
//      neighbour's level.
//   4. A floating contact facing a held one: V = the held V + the difference between the two
//      terminal species' levels beside each (where both are equilibrium links of charged
//      species; a bath without a reference by the φ beside it, its V being its φ).
//   5. Electrode surfaces at their θ₀ (s._surfaceStart); concentrations computed.
//   6. Floating double layers start uncharged: a capacitive contact at zeroCharge + V_T·φ̂
//      beside it; a port with a capacitance at zeroCharge + its window's mean V_T·φ̂ (by
//      volume·area, over the boxes whose charge it holds: where φ is defined, and not an end box
//      whose φ its contact sets).
//   7. A floating port linked to its carrier, with no reactions or capacitance: at the carrier's
//      level across its window, relative to a held contact's (where there is one).
//   8. Gates at α/(α + β) at the starting voltage across them.
//   9. A floating port with reactions and no capacitance: at the voltage where its reactions and
//      its conductance and exchange links pass the current it's driven at, or behind a
//      resistance (V_src − V)/R, by bisection; where none within ±1 kV does, at the V it had.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. A species given by its amount starts at c0 (or its profile) at every node; at a face
//      equilibrated by step 2 (a transient's first step), its amount over the two nodes is kept.
//   2. A species given by its level starts at the level of what feeds it, as step 1 orders them.
//   3. Where a level-fixed charged species responds, each node starts neutral (Σ z c + the fixed
//      charge = 0, plus aσ/F under a held capacitance in a strictly neutral region) to the
//      bisection's tolerance, unless neutrality is out of reach (step 1: then its charge at most
//      ε V_T/(F L²), or the least φ̂ gives, to the bisection's tolerance); a strictly neutral
//      region given by amounts is neutral by its c0, checked (with its aσ under such a window,
//      which φ̂ then sets).
//   4. Steps 2 and 3 change no concentration but as they say: step 3 is a gauge (φ̂ and each η by
//      z·Δ), step 2 moves the passing species between a face's two nodes only.
//   5. Each floating terminal's V is set by the last of steps 4, 6, 7, 9 that applies to it;
//      held terminals' V are untouched.
//   6. Only u, uLo, termV, s.unbalanced, s.coldStart and what computeConcentrations derives are written; the
//      result depends on the model, every terminal's voltage and s.time, deterministically.
// Known, and not an invariant (a cold review; see ROADMAP): the gauge (φ̂ and each η by z·Δ)
// runs from the left, so behind a closed or floating left end the start sits away from a right
// contact's level (55 RT for a bath referenced to Cl⁻), as does a c0 given for a contact-fed
// carrier (φ̂ then from the left's own start); an insertion host's shift moves contact-held
// levels in it too; a terminal behind a resistance starts level with its held neighbour (Solver's
// constructor), which in an unlike cell is its open-circuit voltage away (a closed Cu | Zn cell's
// zinc 1.13 V from rest: the steady solve's held stage lays the start out again at the source,
// src/steady.js); and the φ̂ beside an electrode passes
// over a reaction with a neutral level-held participant. Each the first step absorbs. Where
// ε > 0, a held port's capacitance through its window is left out of the neutrality, so the
// first step moves φ̂ to discharge it and the ions charge it again (a gated Li⁺ bar with εr = 1:
// 0.28 C/m² cold, 0.0035 after 1 ps, back to 0.28 by diffusion; Astra's third pass).
// Tests: every test (each device starts here); test/transient.test.js on what a start must be;
// test/round7-solve.test.js (φ̂ behind a bath or a blocking metal, a floating port's electrode).

import { gateRate } from './boundaries.js';
import { EPS0, FARADAY } from './constants.js';
import { fixedAt, gradedAt, isGraded, profileAt, sourceAt } from './device.js';
import { SolverError } from './errors.js';
import { gaussRows, projectNeutral, strictlyNeutralNodes } from './neutral.js';
import { rescaleAll } from './carry.js';

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
  // The outside level a contact's link offers species i: held there, or exchanged with at a rate
  // (conductance, exchange), whose level a start can take as well (a link of either kind had
  // asked for c0 as if it fed the species through reactions alone).
  const held = (side, i) => {
    const l = model.contacts[side].species[i];
    if (l.type !== 'conductance' && l.type !== 'exchange') return s.contactEta(side, i);
    return z[i] === 0 ? l.mu / model.RT : (z[i] * (s.termV[side === 'left' ? 0 : 1] + l.offset)) / s.VT;
  };
  // Which contact's level a species fed at both ends starts at (see below).
  const anchorSide = model.contacts.right.phi.type === 'bulk' && model.contacts.left.phi.type !== 'bulk' ? 'right' : 'left';
  s.unbalanced = []; // (strictly neutral regions whose start can't be neutral: see neutralPhi)
  // A pinned contact's densities, cRef e^(η − μ°/RT − z φ̂_edge), φ̂_edge = (V − zeroCharge)/V_T:
  // beyond e^700, past what a double holds, the definition can't be solved, and says so here
  // (zeroCharge = V had put a hole contact at 10·e^778 mol/m³, and the solve failed on a
  // singular block without a word why; Neon's review, 1.2).
  for (const [side, k, g] of [['left', 0, 0], ['right', 1, grid.nNodes - 1]]) {
    const ph = model.contacts[side].phi;
    if (ph.type !== 'pinned' || s.floating.includes(k)) continue;
    const mat = gradedAt(materials[regions[grid.nodeRegion[g]].material], regions[grid.nodeRegion[g]], grid.x[g]);
    for (let i = 0; i < n; i++) {
      if (!mat.present[i] || z[i] === 0) continue;
      const ex = s.contactEta(side, i) - mat.mu0[i] / model.RT - (z[i] * (s.termV[k] - ph.zeroCharge)) / s.VT;
      if (ex > 700) {
        throw new SolverError(
          `contacts.${side}: holds ${species[i].name} at cRef·e^${ex.toFixed(0)} at its end, beyond any density (a double holds e^709): ` +
            `its level is ${((s.VT * ex) / z[i]).toPrecision(4)} V from its band at φ_edge = V − zeroCharge = ${(s.termV[k] - ph.zeroCharge).toPrecision(4)} V. ` +
            'zeroCharge is V − φ_edge at zero charge (a barrier, a pzc): check it against the band offsets',
        );
      }
    }
  }
  s.uLo.fill(0);
  const eta = new Float64Array(n), cFix = new Float64Array(n), mode = new Int8Array(n); // 1 level, 2 amount
  // φ̂ runs from the left: from a pinned (or gated) left contact's φ at its zero charge, else 0.
  const leftPhi = model.contacts.left.phi;
  // (Not a floating one, whose V is the start's to find, below: its 0 V as it stands set φ̂, and a
  // cell between a floating gated contact and a platinum held at 2 V started 0.79 V off the
  // platinum's couple, its mirror image at it; Neon's review, 1.3.)
  const leftSets = (leftPhi.type === 'pinned' || leftPhi.type === 'capacitive') && !s.floating.includes(0);
  let phiHat = leftSets ? (s.termV[0] - leftPhi.zeroCharge) / s.VT : 0;
  // Whether anything has set φ̂ yet (a pinned or capacitive left contact, a region's neutrality,
  // an electrode's balance), and the regions before that, which took the running φ̂: 0 behind a
  // bath or an ohmic metal, whatever their voltages, a gauge of its own. Where an electrode on
  // the right is the first thing to set it, those regions move with it (see below).
  let phiSet = leftSets;
  const unset = [];
  const shiftUnset = (d) => {
    for (const q of unset) {
      for (let g = grid.regionStart[q]; g <= grid.regionEnd[q]; g++) {
        const b = s.blockOfNode[g];
        u[b * M] += d;
        for (let i = 0; i < n; i++) if (materials[regions[q].material].present[i]) u[b * M + 1 + i] += z[i] * d;
      }
    }
  };
  // A held port's capacitance through its window, per node: aC (F/m³) and aC (V − zeroCharge),
  // so that aσ = capB − capA V_T φ̂. In a strictly neutral region it's part of each node's
  // neutrality, ρ + aσ = 0: without it, a region of P⁺ alone under a gate was refused for the
  // charge the gate balances, and a level-fed one started off it, its first 1 ps step jumping
  // the gate's charge from 0.28 to 0.0035 C/m² (Astra's third pass). Not at an end box whose φ
  // its contact sets, as _portCapacitance has it.
  const { capA, capB } = heldCapacitance(s);
  const setsPhi = (side) => ['pinned', 'bulk'].includes(model.contacts[side].phi.type);
  for (let r = 0; r < regions.length; r++) {
    if (r > 0) phiHat += interfaces[r - 1].dipole / s.VT;
    // (a graded region's material as it is at each node: here at its start, below node by node)
    const reg = regions[r], base = materials[reg.material], graded = isGraded(reg);
    let mat = graded ? gradedAt(base, reg, grid.x[grid.regionStart[r]]) : base;
    mode.fill(0);
    const profiled = reg.c0Profile.some((p) => p !== null) || reg.fixedProfile !== null;
    let capped = false; // (a held capacitance in a strictly neutral region: then node by node)
    if (base.epsr === 0 && !base.conductor && !base.phiFree) for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) if (capA[g] > 0) capped = true;
    let aC = 0, aCV = 0; // (capA, capB at the node, where capped)
    let qFixed = fixedAt(reg, grid.x[grid.regionStart[r]]) / FARADAY; // (mol/m³; per node where it's a profile)
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
          // (the face's surface species at their starting coverages, as below: left out, a cold
          // solve at a floating Pt with H(ads) beside Fe³⁺/Fe²⁺ started 93 RT off, and converged or
          // failed by which reaction was listed first; Argon's review of 1.6.1)
          const surf = interfaces[r - 1].surface ?? [], bare = 1 - surf.reduce((t, sp) => t + sp.theta0, 0);
          for (const q of rx.surf ?? []) a -= q.nu * (surf[q.s].mu0 / model.RT + Math.log(surf[q.s].theta0 / bare));
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
        // (a contact's only where the stretch reaches it: a species held by a port behind a
        // blocked face took the right contact's level, 2.63 for the port's 1, and kept it for
        // ~c/k of the exchange, an early transient made up; Astra's third pass)
        // Fed at both ends, from the end that's plain neutral bulk at its contact's levels (an ohmic
        // contact, a bath), where only one is: a pinned or gated end is a barrier or a double layer,
        // the region beyond it at the other end's levels. Taken from the left whatever the ends,
        // an illuminated GaAs | InGaAs junction solved cold with its Schottky contact on the right,
        // and with it on the left started the GaAs's electrons at the Schottky's level, 59 thermal
        // units from their own ohmic contact's, and was refused (round 9's mirror property).
        const first = anchorSide === 'right' ? ['right', st.rightOpen] : ['left', st.leftOpen];
        const second = anchorSide === 'right' ? ['left', st.leftOpen] : ['right', st.rightOpen];
        const a = first[1] ? held(first[0], i) : NaN;
        eta[i] = Number.isFinite(a) ? a : second[1] ? held(second[0], i) : NaN;
        // (a held port's level before a floating one's, whatever order lists them)
        const port = [...st.linked].sort((a, b) => s.floating.includes(2 + a) - s.floating.includes(2 + b))[0];
        if (!Number.isFinite(eta[i]) && port !== undefined) eta[i] = s.portEta(model.ports[port], i);
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
    let clamped = false; // (whether charge() just clamped an exponent: its value then a plateau, not a floor)
    const charge = (ph) => {
      if (!mat.ideal) {
        // A model can refuse a trial φ̂ far out, where one ion is many orders denser than the rest
        // (Debye–Hückel's limiting law threw out of the start, for KCl at 10 mM, its neutral φ̂
        // 57 units from where the search began): there the ideal charge has the sign, and the
        // search goes on by it.
        let refused = false;
        try {
          materialAt(s, mat, reg.background, eta, mode, cFix, ph, zeta, cc);
        } catch (err) {
          if (!(err instanceof SolverError)) throw err;
          refused = true;
        }
        if (!refused) {
          let q = qFixed + (aCV - aC * s.VT * ph) / FARADAY;
          for (let i = 0; i < n; i++) if (mat.present[i]) q += z[i] * cc[i];
          return q;
        }
      }
      let q = qFixed + (aCV - aC * s.VT * ph) / FARADAY;
      clamped = false;
      for (let i = 0; i < n; i++) {
        if (mode[i] === 2) q += z[i] * cFix[i];
        else if (mode[i] === 1 && z[i] !== 0) {
          const raw = eta[i] - mat.mu0[i] / model.RT - z[i] * ph, ex = Math.min(700, Math.max(-700, raw));
          if (ex !== raw) clamped = true;
          q += z[i] * mat.cRef[i] * Math.exp(ex);
        }
      }
      return q;
    };
    const responds = mode.some((m, i) => m === 1 && z[i] !== 0);
    const mustBalance = !responds && !mat.conductor && !mat.phiFree && mat.epsr === 0;
    const checkNeutral = (where) => {
      // Strictly neutral, with every charged species' amount given: its c0 must be neutral.
      let q = qFixed, scale = Math.abs(q);
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
    if (mustBalance && !profiled && !capped) checkNeutral('');
    // Each amount-given species' ζ (its activity's ln), for a balance below: no level-fixed charged
    // species responding, they don't depend on φ̂ (a non-ideal model's from its statistics).
    const activity = (i) => (mat.ideal ? Math.log(cFix[i] / mat.cRef[i]) : zeta[i]);
    // φ̂ where the first of an electrode's reactions that can be balanced is at equilibrium: its
    // participants on this region's side (`inside`) at their amounts, the metal's at levelOf(p)
    // (NaN where it has none), its surface species at their starting coverages (a face's in
    // rx.surf, a port's in part with side 2); or ph where none can be (`hit` says which). (Written
    // out for each kind of electrode before, a face's coverage fix made in one copy at a time; a
    // code-smell review.)
    let hit = false;
    const balanced = (reactions, surface, inside, levelOf, ph) => {
      const bare = 1 - surface.reduce((t, sp) => t + sp.theta0, 0);
      const lnSurf = (q) => surface[q].mu0 / model.RT + Math.log(surface[q].theta0 / bare);
      for (const rx of reactions) {
        let a = rx.fixedA, slope = 0;
        for (const q of rx.surf ?? []) a -= q.nu * lnSurf(q.s);
        for (const p of rx.part) {
          if (p.side === 2) a -= p.nu * lnSurf(p.s);
          else if (p.side !== inside) a -= p.nu * levelOf(p);
          else if (mode[p.i] === 2) {
            a -= p.nu * (activity(p.i) + mat.mu0[p.i] / model.RT);
            slope += p.nu * z[p.i];
          } else slope = NaN;
        }
        if (slope !== 0 && Number.isFinite(slope) && Number.isFinite(a)) {
          hit = true;
          return a / slope; // a − slope·φ̂ = 0
        }
      }
      return ph;
    };
    // The same, with the amounts and the material where the electrode is (a c0 profile's there, a
    // graded region's μ° and c_ref there), not at the region's start (a metal on the right of a salt
    // ramped 10 to 90 mM started 2.2 RT from its equilibrium, the mirror device at it, and a window
    // at the far end passing 8e-4 mol/(m²·s) cold; Astra's third pass. A port mid-film in a region
    // graded in c_ref started 0.6 RT off, its first transient's current made up; Sol's fifth pass).
    const balancedAt = (x, reactions, surface, inside, levelOf, ph) => {
      if (!graded && !reg.c0Profile.some((p) => p)) return balanced(reactions, surface, inside, levelOf, ph);
      const keep = [cFix.slice(), zeta.slice(), cc.slice()], keepMat = mat;
      if (graded) mat = gradedAt(base, reg, x);
      for (let i = 0; i < n; i++) if (reg.c0Profile[i]) cFix[i] = profileAt(reg.c0Profile[i], x);
      if (!mat.ideal) materialAt(s, mat, reg.background, eta, mode, cFix, ph, zeta, cc);
      const out = balanced(reactions, surface, inside, levelOf, ph);
      [cFix, zeta, cc].forEach((a, k) => a.set(keep[k]));
      mat = keepMat;
      return out;
    };
    if (!mat.phiFree && !responds && !mat.ideal) materialAt(s, mat, reg.background, eta, mode, cFix, phiHat, zeta, cc);
    hit = false;
    if (!mat.phiFree && !responds && r > 0 && materials[regions[r - 1].material].conductor) {
      // Nothing here fixes φ, but an electrode on the left does: start with its first reaction
      // that takes a species from this side at equilibrium (the electrode at its open-circuit
      // level), rather than φ carried over from the metal, which can be volts away.
      const bn = s.blockOfNode[grid.regionEnd[r - 1]];
      phiHat = balanced(interfaces[r - 1].reactions, interfaces[r - 1].surface ?? [], 1, (p) => u[bn * M + 1 + p.i], phiHat);
    }
    if (!hit && !mat.phiFree && !responds && (r === 0 || !phiSet) && r + 1 < regions.length && materials[regions[r + 1].material].conductor) {
      if (graded) mat = gradedAt(base, reg, grid.x[grid.regionEnd[r]]); // (the electrode's face is this region's end)
      // Likewise an electrode on the right, its level its contact's, where nothing on the left
      // sets φ̂ (the stack written the other way round started at φ̂ = 0, volts from the
      // electrode's level, and failed; Argon's review). Further in, φ̂ carried over from the left
      // is kept where something there set it: a Daniell cell's KCl bridge balanced against its
      // copper instead started with a jump at the junction, and a 1 µs step failed to settle.
      // Where nothing did, the regions before take this one's φ̂ too (a gauge: φ̂ and each η by
      // z·Δ, no concentration changes): behind a bath, two solutions in a row started at φ̂ = 0,
      // 2 V from their copper's level when every voltage was moved by 2 V, and failed where one
      // solution, or the mirror image, converged at any shift; and behind a blocking platinum,
      // which balances nothing, a solution started 0.83 V from its zinc's level (the echem
      // combos' cases 212 and 341, cut down).
      // (its carrier at its contact's level, the metal fed by the right contact)
      const fed = (p) => {
        const st = s.stretches[s.stretchOf[(r + 1) * n + p.i]];
        return st?.rightOpen && st.contactFed ? held('right', p.i) : NaN; // (held: a conductance or exchange link's level too)
      };
      const was = phiHat;
      phiHat = balancedAt(grid.x[grid.regionEnd[r]], interfaces[r].reactions, interfaces[r].surface ?? [], 0, fed, phiHat);
      if (hit && !phiSet) shiftUnset(phiHat - was);
    }
    // (A held one before a floating one, whatever order lists them: the first listed, floating at
    // its 0 V, started a solution 2 V from the held electrode's couple, and at 5 V was refused,
    // where listed the other way it solved in one iteration; Neon's review, 1.3.)
    const electrode = model.ports.filter((port) => port.region === r && port.reactions.length > 0).sort((a, b) => s.floating.includes(2 + model.ports.indexOf(a)) - s.floating.includes(2 + model.ports.indexOf(b)))[0];
    // (A floating one's voltage is the start's to find, below, from this φ̂: where something has
    // set φ̂ already, its V as it stands, 0 V, doesn't override it. It did, and a porous Fe³⁺/Fe²⁺
    // electrode at I = 0 beside a zinc moved the solution's φ̂ from the zinc's 0.827 V to
    // −0.725, the couple at 0 V, its mixed potential then found there: the solve, started 1.55 V
    // from the zinc's level, failed, where the same port held at its couple's level converged;
    // so with Ag/AgCl, a double layer and a Pt beside it, the echem combos' case 65.)
    const floats = electrode && s.floating.includes(2 + model.ports.indexOf(electrode));
    if (!mat.phiFree && !responds && electrode && !(floats && (hit || phiSet))) {
      // Nor here, but an electrode spread through it does: start with the port's first reaction
      // that can be balanced at equilibrium with the port's level (the electrode at open circuit),
      // at its window's middle.
      const mid = 0.5 * (grid.x[electrode.nodes[0]] + grid.x[electrode.nodes[electrode.nodes.length - 1]]);
      const was = phiHat;
      phiHat = balancedAt(mid, electrode.reactions, electrode.surface, 0, (p) => s.portEta(electrode, p.i), phiHat);
      // (and the regions before that nothing set move with it, as for an electrode on the right:
      // two solutions with the electrode in the second started the first at φ̂ = 0, 1.2 V from
      // the couple's level, and at 10 V were refused, where the mirror image solved in 7; Neon's
      // review, 1.3)
      if (hit && !phiSet) shiftUnset(phiHat - was);
    }
    if (hit || (!mat.phiFree && (responds || capped))) phiSet = true;
    else if (!phiSet && !mat.phiFree) unset.push(r);
    const neutralPhi = (ph) => {
      let lo = ph - 1, hi = ph + 1;
      while (charge(lo) < 0 && lo > -1e4) lo -= 2 * (hi - lo);
      while (charge(hi) > 0 && hi < 1e4) hi += 2 * (hi - lo);
      // Out of reach (every responding carrier has one sign and nothing balances it, as in an
      // undoped layer that holds only holes): not run to ±∞, but where the charge is no more than
      // the region holds with a thermal voltage across it, ε V_T/(F L²), or as little as it gets
      // (a fixed charge left). Kept at the running φ̂ (a contact's level), a hole-only layer 100
      // nm thick started full of the anode's 1e27 m⁻³ holes, 7 kV of space charge: with traps in
      // it, every step failed from there.
      if (charge(lo) < 0 || charge(hi) > 0) {
        // In a strictly neutral region that's a start that can't be: the amounts c0 gives leave
        // a charge no level-held species can balance. Noted, for a transient that fails from it
        // to say so (a steady solve, its levels from the contacts, can still find the state).
        if (mat.epsr === 0) {
          s.unbalanced.push(r);
          return ph;
        }
        const L = grid.x[grid.regionEnd[r]] - grid.x[grid.regionStart[r]];
        const most = (mat.epsr * EPS0 * model.RT) / (FARADAY * FARADAY * L * L);
        const q0 = charge(ph), dir = q0 > 0 ? 1 : -1; // (the charge falls as φ̂ rises)
        if (!(Math.abs(q0) > most)) return ph;
        let a = ph, qa = Math.abs(q0);
        for (let step = 1; step < 1e4; step *= 2) {
          let b = a + dir * step;
          const qb = Math.abs(charge(b));
          if (qb <= most) {
            for (let it = 0; it < 100 && Math.abs(b - a) > 1e-9; it++) {
              const m = 0.5 * (a + b);
              if (Math.abs(charge(m)) > most) a = m;
              else b = m;
            }
            return b;
          }
          // (a fixed charge: as little as it gets; but not where the charge is clamped flat, beyond
          // e^700: a hole-only layer with 20 V across it, fed at its band's level, stopped there,
          // on the first step, and started at c = ∞; Neon's review, 1.2)
          if (!(qb < qa * (1 - 1e-9)) && !clamped) return a;
          [a, qa] = [b, qb];
        }
        return a;
      }
      for (let it = 0; it < 200 && hi - lo > 1e-12; it++) {
        const m = 0.5 * (lo + hi);
        if (charge(m) > 0) lo = m;
        else hi = m;
      }
      return 0.5 * (lo + hi);
    };
    if (!mat.phiFree && responds) phiHat = neutralPhi(phiHat);
    // A held gate over part of a region where nothing responds: its window's neutrality sets φ̂,
    // and the nodes either side of it take that (nothing else sets theirs), as does what follows.
    // (Taken node by node, the nodes before the window kept the running φ̂ and those after it the
    // gate's, 2 V apart, which side by which way the stack was written; Neon's review, 1.5.)
    let phiGate = NaN;
    if (capped && !responds && !mat.phiFree) {
      const gw = grid.regionStart[r] + capA.subarray(grid.regionStart[r], grid.regionEnd[r] + 1).findIndex((a) => a > 0);
      if (graded) mat = gradedAt(base, reg, grid.x[gw]);
      qFixed = fixedAt(reg, grid.x[gw]) / FARADAY;
      [aC, aCV] = [capA[gw], capB[gw]];
      for (let i = 0; i < n; i++) if (reg.c0Profile[i]) cFix[i] = profileAt(reg.c0Profile[i], grid.x[gw]);
      phiGate = neutralPhi(phiHat); // (the window's own nodes still find theirs from the running φ̂, below)
      [aC, aCV] = [0, 0];
    }
    if (!mat.ideal) materialAt(s, mat, reg.background, eta, mode, cFix, phiHat, zeta, cc);
    for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
      if (profiled || graded || capped) {
        if (graded) mat = gradedAt(base, reg, grid.x[g]);
        qFixed = fixedAt(reg, grid.x[g]) / FARADAY;
        if (capped) [aC, aCV] = [capA[g], capB[g]];
        for (let i = 0; i < n; i++) if (reg.c0Profile[i]) cFix[i] = profileAt(reg.c0Profile[i], grid.x[g]);
        if (mustBalance && !(aC > 0)) checkNeutral(` at x = ${grid.x[g].toPrecision(4)} m`);
        if (!mat.phiFree && (responds || aC > 0)) phiHat = neutralPhi(phiHat);
      }
      // (outside a gate's window, its φ̂; but not at an end box whose φ its contact sets or gates)
      const own = (side) => setsPhi(side) || model.contacts[side].phi.type === 'capacitive';
      const ends = (g === 0 && own('left')) || (g === grid.nNodes - 1 && own('right'));
      const ph = Number.isFinite(phiGate) && !(aC > 0) && !ends && Math.abs(phiGate - phiHat) > 1e-6 ? phiGate : phiHat; // (within round-off, the running one)
      if ((profiled || graded || capped) && !mat.ideal) materialAt(s, mat, reg.background, eta, mode, cFix, ph, zeta, cc);
      const b = s.blockOfNode[g];
      u[b * M] = mat.conductor ? 0 : ph; // a metal's slot 0 is its segment flux
      for (let i = 0; i < n; i++) {
        if (mode[i] === 1) u[b * M + 1 + i] = eta[i];
        else if (mode[i] === 2) {
          const zt = mat.ideal ? Math.log(cFix[i] / mat.cRef[i]) : zeta[i];
          u[b * M + 1 + i] = zt + mat.mu0[i] / model.RT + z[i] * ph;
        }
      }
    }
  }
  // An insertion host's φ is a gauge: c0 sets only its neutral combination's level, not how that
  // splits between its ion and its electrons. Start its ion level with the same ion's beside it
  // across a face (the equilibrium a link or a transfer reaction heads for), so its electrons,
  // and a floating contact on them, start where the cell's open-circuit voltage puts them, not
  // volts away. (A shift of φ̂ and each η by z·Δ: the concentrations don't change.) Not through
  // a metal collector between the host and the contact: the contact then starts at 0 V, and the
  // first step or the steady solve brings it to the open circuit (a cold review's p9: 4.09 V
  // either way, in the same 14 steps).
  regions.forEach((reg, r) => {
    const mat = materials[reg.material];
    if (!mat.phiFree || mat.conductor) return; // (a metal is φ-free too, but no host: its level is its contact's or its faces')
    const sides = [[r - 1, grid.regionStart[r], r > 0 ? grid.regionEnd[r - 1] : -1], [r + 1, grid.regionEnd[r], r + 1 < regions.length ? grid.regionStart[r + 1] : -1]];
    for (const [nr, gIn, gOut] of sides) {
      if (gOut < 0) continue;
      const other = materials[regions[nr].material];
      if (other.phiFree || other.conductor) continue;
      // (the first the face passes, by a link or a reaction that takes it from one side to the
      // other, else the first shared: with e⁻ listed first and blocked at the face, a host beside
      // a mixed conductor started its Li⁺ 0.4 V out of equilibrium, and its wire at 0 V for the
      // 0.4 V it starts at in the other order; Astra's third pass)
      const itf = interfaces[Math.min(r, nr)], shared = (j) => z[j] !== 0 && mat.present[j] && other.present[j];
      const passes = (j) => itf.links[j].type !== 'blocked' || itf.reactions.some((rx) => rx.part.some((p) => p.i === j && p.side === 0) && rx.part.some((p) => p.i === j && p.side === 1));
      let i = z.findIndex((_, j) => shared(j) && passes(j));
      if (i < 0) i = z.findIndex((_, j) => shared(j));
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
    const ct = model.contacts[side], i = ct.terminal, g = side === 'left' ? 0 : grid.nNodes - 1, b = s.blockOfNode[g];
    // A bath without a reference: its V is its φ, so the end's φ, not its terminal's level (the
    // first listed of its species: a salt whose c0 differs from the bath's read the floating end
    // at +56 or −56 mV by the order of def.species, the first step at 0; Astra's third pass).
    if (ct.she && !ct.she.reference && !s.phiUndefined[g]) return s.VT * u[b * M];
    if (i === null || ct.species[i].type !== 'equilibrium' || z[i] === 0) return NaN;
    return (s.VT * u[b * M + 1 + i]) / z[i] - ct.species[i].offset;
  };
  for (const k of s.floating) {
    const t = s.terms[k], held = k === 0 ? 1 : 0;
    if (t.kind === 'port' || held >= s.terms.length || s.floating.includes(held)) continue;
    const shift = level(t.side) - level(s.terms[held].side);
    if (Number.isFinite(shift)) s.termV[k] = s.termV[held] + shift;
  }
  for (let g = 0; g < s.nNodes; g++) if (s.surfAt[g] >= 0) s._surfaceStart(g);
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
    // (over the boxes whose charge it holds: an end box whose φ its contact sets is the
    // contact's, as _portCapacitance has it; counted, a port beside a pinned end started with an
    // O(h) charge, 0.0087 C/m² on two cells; Astra's third pass)
    let w = 0, sum = 0;
    port.nodes.forEach((g, j) => {
      if (!s._portHolds(g)) return;
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
      const a = s.model.grid.vol[g] * (s.portArea[t.index]?.[j] ?? port.fill[j]); // (a window with no area: by volume)
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
    // What it must pass: its driven current, or behind a resistance what the circuit says,
    // (V_src − V)/R (taken as 0 there, both solves failed for R ≤ 10 Ω·m²: a cold review).
    const src = sourceAt(t.drive.src, s.time), want = (V) => (t.drive.kind === 'I' ? src : (src - V) / t.drive.R);
    const current = (V) => {
      s.termV[k] = V;
      let I = 0;
      port.reactions.forEach((rx) => {
        let q = 0;
        for (const p of rx.part) if (p.side === 0) q += p.nu * z[p.i];
        port.nodes.forEach((g, w) => (I += q * FARADAY * s.model.grid.vol[g] * s.portArea[t.index][w] * s._portRate(port, rx, g)));
      });
      // and what its links to the region's species pass, as the port's terms have it (left out, a
      // floating electrode with a conductance beside its reaction started passing 20 A/m² where
      // it was set to none; Astra's third pass)
      for (let i = 0; i < n; i++) {
        const link = port.species[i];
        if (z[i] === 0 || (link.type !== 'conductance' && link.type !== 'exchange')) continue;
        const kk = link.type === 'conductance' ? (link.G * s.VT) / (z[i] * z[i] * FARADAY) : link.k, gated = link.gates ?? [], q0 = 1 + n + port.surface.length;
        port.nodes.forEach((g, w) => {
          const b = s.blockOfNode[g], open = gated.length ? s._open(gated, (q) => u[b * M + q0 + q]) : 1;
          I += z[i] * FARADAY * s.model.grid.vol[g] * port.fill[w] * kk * open * (s.portEta(port, i) - u[b * M + 1 + i]);
        });
      }
      return I;
    };
    // Current into the device rises with V (oxidation); bracket, then bisect. Where no V passes
    // it (a drive beyond what the reactions and links can carry), it stays where it was: a
    // search run out to its bound had started an MHC electrode set beyond its plateau at 1457 V.
    const V0 = s.termV[k];
    let lo = V0 - 1, hi = V0 + 1;
    while (current(lo) > want(lo) && lo > -1e3) lo -= 2 * (hi - lo);
    while (current(hi) < want(hi) && hi < 1e3) hi += 2 * (hi - lo);
    if (!(current(lo) <= want(lo) && current(hi) >= want(hi))) {
      s.termV[k] = V0;
      continue;
    }
    for (let it = 0; it < 200 && hi - lo > 1e-12; it++) {
      const m = 0.5 * (lo + hi);
      if (current(m) < want(m)) lo = m;
      else hi = m;
    }
    s.termV[k] = 0.5 * (lo + hi);
  }
  s.coldStart = true; // (until a solve or a step moves it: see equilibrateFaces)
}

/**
 * The t = 0⁺ state at the faces between strictly neutral regions, taken by the first time step
 * from a cold start (src/timeSteps.js) and not by a steady solve, which needs no t = 0⁺: there it
 * was the worse start, its face nodes moved while the rest of a region stayed (a stack whose
 * middle region started 2.2 V from its neighbour's level failed every steady solve, where the
 * cold start converged in one).
 */
export function equilibrateFaces(s) {
  const { model, n, M, u, z } = s;
  const { regions, interfaces, materials } = model;
  const grid = model.grid;
  s.coldStart = false;
  const { capA, capB } = heldCapacitance(s);
  // Two strictly neutral regions meeting at a face whose links hold ions in equilibrium: where
  // c0 leaves the two sides out of it (a different fixed charge or immobile ion on each, so
  // c₊c₋ differs across it), the face reaches it at t = 0⁺, before anything diffuses, and a
  // first step made to take that jump met an exactly singular pivot. Start there: the ions the
  // face passes shared between its two nodes in equilibrium (μ̄ equal across it), each one's
  // amount over the two nodes kept, everything else as c0 has it. The one unknown is the jump
  // in φ̂ across it (a Donnan potential), set by neutrality on one side; the other's follows,
  // the two nodes' charge summing to what it was, zero.
  for (let f = 0; f < interfaces.length; f++) {
    const itf = interfaces[f], gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
    const mL = gradedAt(materials[regions[f].material], regions[f], grid.x[gL]), mR = gradedAt(materials[regions[f + 1].material], regions[f + 1], grid.x[gR]);
    const strict = (m) => m.epsr === 0 && !m.conductor && !m.phiFree && m.ideal;
    if (itf.phi.type !== 'neutral' || !strict(mL) || !strict(mR)) continue;
    const bL = s.blockOfNode[gL] * M, bR = s.blockOfNode[gR] * M;
    const pass = [];
    for (let i = 0; i < n; i++) if (itf.links[i].type === 'equilibrium' && mL.present[i] && mR.present[i]) pass.push(i);
    if (!pass.some((i) => z[i] !== 0)) continue;
    const RT = model.RT, vL = grid.vol[gL], vR = grid.vol[gR];
    const conc = (m, b, i) => m.cRef[i] * Math.exp(u[b + 1 + i] - m.mu0[i] / RT - z[i] * u[b]);
    const amount = pass.map((i) => vL * conc(mL, bL, i) + vR * conc(mR, bR, i));
    let qOther = (fixedAt(regions[f], grid.x[gL]) + floatSigma(s, gL)) / FARADAY; // (the left node's charge from all it doesn't share)
    for (let i = 0; i < n; i++) if (mL.present[i] && !pass.includes(i) && z[i] !== 0) qOther += z[i] * conc(mL, bL, i);
    // A held gate on either node: its charge is part of that node's neutrality, ρ + aσ = 0. The
    // ions' total over the two nodes is kept, so the gates' is too: vL·aC_L·φ̂_L + vR·aC_R·φ̂_R
    // stays as it was, φ̂_L moving with the jump where the right node is gated, and with the left
    // node neutral the right one is. (The ions alone balanced, a gated left node gave up its
    // gate's counter-ions to the ungated right one, and the first step lost them; Neon's review, 1.1.)
    const aL = capA[gL], aR = capA[gR], S = vL * aL + vR * aR, phiL0 = u[bL], D0 = u[bR] - u[bL];
    const phiLAt = (D) => (S > 0 ? phiL0 - ((vR * aR) / S) * (D - D0) : phiL0); // (falls as D rises)
    const gateL = (D) => (capB[gL] - aL * s.VT * phiLAt(D)) / FARADAY; // (aσ/F at the left node: rises with D)
    // At a jump D = φ̂_R − φ̂_L: c_L/c_R = K = (cRef_L/cRef_R)·exp((μ°_R − μ°_L)/RT + z D), the
    // amount fixing both; in logs, ln c = ln amount − ln(v_L + v_R/K) and so on, so that a K
    // beyond a double's range (a divalent ion as the bracket widens) gives c → 0, not ∞·0 = NaN,
    // which had ended the search and left the face as c0 had it (Neon's review, 1.4).
    const lae = (a, b) => Math.max(a, b) + Math.log1p(Math.exp(-Math.abs(a - b))); // ln(e^a + e^b)
    const lnSplit = (D) =>
      pass.map((i, j) => {
        const lnK = Math.log(mL.cRef[i] / mR.cRef[i]) + (mR.mu0[i] - mL.mu0[i]) / RT + z[i] * D, lnA = Math.log(amount[j]);
        return [lnA - lae(Math.log(vL), Math.log(vR) - lnK), lnA - lae(Math.log(vL) + lnK, Math.log(vR))];
      });
    const split = (D) => lnSplit(D).map(([a, b]) => [Math.exp(a), Math.exp(b)]);
    const qL = (D) => split(D).reduce((q, [cL], j) => q + z[pass[j]] * cL, qOther + gateL(D)); // (increasing in D)
    let lo = D0 - 1, hi = D0 + 1;
    while (qL(lo) > 0 && lo > D0 - 1e3) lo -= 2 * (hi - lo);
    while (qL(hi) < 0 && hi < D0 + 1e3) hi += 2 * (hi - lo);
    if (!(qL(lo) <= 0 && qL(hi) >= 0)) continue; // (nothing the face shares can balance it: leave c0)
    for (let it = 0; it < 200 && hi - lo > 1e-13 * Math.max(1, Math.abs(D0)); it++) {
      const m = 0.5 * (lo + hi);
      if (qL(m) < 0) lo = m;
      else hi = m;
    }
    const D = 0.5 * (lo + hi), cs = lnSplit(D); // (ln c_L, ln c_R)
    const keep = Float64Array.from({ length: n }, (_, i) => (mR.present[i] ? conc(mR, bR, i) : 0)); // (the right node's, before φ̂ moves)
    const phiR0 = u[bR];
    const dL = phiLAt(D) - phiL0; // (the left node's φ̂, moved where a gate on the right asks it: what it doesn't share kept)
    u[bL] += dL;
    for (let i = 0; i < n; i++) if (mL.present[i] && !pass.includes(i)) u[bL + 1 + i] += z[i] * dL;
    u[bR] = u[bL] + D;
    for (let i = 0; i < n; i++) if (mR.present[i] && !pass.includes(i)) u[bR + 1 + i] = Math.log(keep[i] / mR.cRef[i]) + mR.mu0[i] / RT + z[i] * u[bR];
    pass.forEach((i, j) => {
      if (!(amount[j] > 0)) return; // (nothing of it at either node: its levels as they were)
      u[bL + 1 + i] = cs[j][0] - Math.log(mL.cRef[i]) + mL.mu0[i] / RT + z[i] * u[bL];
      u[bR + 1 + i] = cs[j][1] - Math.log(mR.cRef[i]) + mR.mu0[i] / RT + z[i] * u[bR];
    });
    // A floating gate over either node keeps the charge it started with (none: it passes no
    // current, at t = 0⁺ as after): its voltage moves with the φ̂ the face moved under it, and the
    // rest of its window's φ̂ with it, concentrations kept, so its σ is as it was at every node and
    // each stays neutral. (Left where it was, a floating gate over a Donnan face's right side took
    // the jump as a charge its node's ions didn't balance: the first step drove 1.6 MA/m² through
    // the bath, and an advance() from the start failed at t = 0; Neon's third review, R3 1.7.)
    followFloating(s, gL, dL);
    followFloating(s, gR, u[bR] - phiR0);
  }
  s.computeConcentrations();
  neutralStart(s);
}

// Every strictly neutral node onto its Gauss row (src/neutral.js) where the start left one off:
// by its φ̂ at its levels, each closed stretch's amount held at what the start gave it (rescaleAll,
// src/carry.js). (A capacitive face's sheet charge between two strictly neutral regions has no
// neutral start at c0: its edge boxes hold its counter-charge. Left off, the first step charged
// them, 1.6e-2 C/m², through the box beside them.)
// Whether node g's box holds a capacitive face's or contact's displacement.
function boundaryCapacitive(s, g) {
  const grid = s.model.grid, r = grid.nodeRegion[g], itfs = s.model.interfaces;
  if (g === grid.regionEnd[r] && r < itfs.length && itfs[r].phi.type === 'capacitive') return true;
  if (g === grid.regionStart[r] && r > 0 && itfs[r - 1].phi.type === 'capacitive') return true;
  const side = g === 0 ? 'left' : g === s.nNodes - 1 ? 'right' : null;
  return !!side && s.model.contacts[side].phi.type === 'capacitive';
}

function neutralStart(s) {
  const nodes = strictlyNeutralNodes(s);
  // (only where a box beside a capacitive face or contact is off, past the start's own round-off:
  // a box whose own c0 can't be neutral is the start's to refuse, and says so; and its solves leave
  // 1e-12 of the gross, which isn't worth moving a closed amount for)
  if (nodes.length === 0 || !gaussRows(s, nodes).some((r, a) => r.slope > 0 && Math.abs(r.G) > 1e-10 * r.scale && boundaryCapacitive(s, nodes[a].g))) return;
  const u0 = Float64Array.from(s.u);
  const list = s.stretches.flatMap((st) => (st.connected || st.reactive ? [] : [{ st, want: s.amount(st) }]));
  // (where every ion is closed, their net charge is too: there the potential stays and the
  // concentrations move in the charge direction, as in the carry; by φ̂ at the levels, the
  // amounts and neutrality chased each other, a closed one 7e-4 off after 60 rounds)
  const closed = new Set(list.map(({ st }) => st));
  for (const nd of nodes) {
    const r = s.model.grid.nodeRegion[nd.g];
    let all = !nd.bulk;
    for (let i = 0; i < s.n && all; i++) if (s.present[nd.g * s.n + i] && s.z[i] !== 0 && !closed.has(s.stretches[s.stretchOf[r * s.n + i]])) all = false;
    if (all) nd.mode = 'charge';
  }
  const moved = new Uint8Array(s.nNodes);
  for (let round = 0; round < 60; round++) {
    projectNeutral(s, nodes, { moved });
    if (!(rescaleAll(s, list, moved) >= 1e-13)) break;
  }
  projectNeutral(s, nodes, { moved });
  // (a closed amount it can't keep exactly with every node neutral: the start as it was, the
  // first step left to take the difference, as before)
  if (list.some(({ st, want }) => !(Math.abs(s.amount(st) - want) <= 1e-12 * Math.abs(want)))) {
    s.u.set(u0);
    s.computeConcentrations();
  }
}

// The charge per volume a floating capacitive port puts at node g, Σ aσ (C/m³), where it holds
// one (see _portHolds).
function floatSigma(s, g) {
  const { model, M } = s;
  let q = 0;
  for (const [k, w] of s.qgTerms[g]) {
    const port = model.ports[k - 2];
    if (!s.floating.includes(k) || !s._portHolds(g)) continue;
    q += s.portArea[k - 2][w] * port.capacitance.C * (s.termV[k] - port.capacitance.zeroCharge - s.VT * s.u[s.blockOfNode[g] * M]);
  }
  return q;
}

// Each floating capacitive port over node g, whose φ̂ moved by d: its voltage moved by V_T d, and
// φ̂ at the rest of its window too, each η_i by z_i d (its concentrations kept), so that its σ is
// as it was everywhere.
function followFloating(s, g, d) {
  if (d === 0) return;
  const { model, M, n, z, u } = s;
  for (const [k] of s.qgTerms[g]) {
    if (!s.floating.includes(k) || !s._portHolds(g)) continue;
    s.termV[k] += s.VT * d;
    for (const h of model.ports[k - 2].nodes) {
      if (h === g || !s._portHolds(h)) continue;
      const b = s.blockOfNode[h] * M;
      u[b] += d;
      for (let i = 0; i < n; i++) u[b + 1 + i] += z[i] * d;
    }
  }
}

// A held port's capacitance through its window, per node: capA = aC (F/m³) and capB = aC (V −
// zeroCharge), so that aσ = capB − capA V_T φ̂. Not at an end box whose φ its contact sets, as
// _portCapacitance has it, nor a floating port's (its V is the start's to find).
function heldCapacitance(s) {
  const { model } = s, grid = model.grid;
  const capA = new Float64Array(grid.nNodes), capB = new Float64Array(grid.nNodes);
  model.ports.forEach((port, k) => {
    if (!port.capacitance || s.floating.includes(2 + k)) return;
    port.nodes.forEach((g, w) => {
      if (!s._portHolds(g)) return; // (the boxes _portCapacitance puts its charge in: one predicate)
      const aC = s.portArea[k][w] * port.capacitance.C;
      capA[g] += aC;
      capB[g] += aC * (s.termV[2 + k] - port.capacitance.zeroCharge);
    });
  });
  return { capA, capB };
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
