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
//   1. Anchor-outward layout, mirror-symmetric and gauge-covariant:
//      a) Species modes and levels, for each region and species present:
//         - its c0 (or c0Profile) where given: an amount, mode 2, the starting state even where a
//           contact feeds the species;
//         - else a conductor's carrier away from the contacts: interior conductor, resolved in (d)
//           from whichever adjacent face(s) have a reaction taking it, else the adjacent φ̂;
//         - else not fed by a contact: c0 is required (an error without it);
//         - else fed by contacts/ports (held or exchanged at a rate: conductance, exchange alike):
//           held sources before floating ones, contacts before ports, a 'bulk' contact before a
//           non-'bulk' one; where two equal-authority contacts feed a pure-conductor stretch (or a
//           single uniformly-doped region between two 'bulk' contacts) at different levels, the
//           level interpolates ohmically in ∫dx/σ, and on any other non-conductor stretch takes
//           their symmetric midpoint (avoiding an artificial bulk current gradient across blocking
//           junctions); multiple ports interpolate piecewise in ∫dx/σ (charged) or ∫dx/D (neutral);
//           fed by none of them (reactions only), from c0, required.
//         A pinned, held contact whose own density, cRef e^(η − μ°/RT − z φ̂_edge), is beyond e^700
//         is refused (SolverError). A strictly neutral region (ε = 0) with every charged species'
//         amount given must have a neutral c0 outside any held port capacitance's window: an error
//         otherwise.
//      b) Held insertion hosts (φ-free, non-conductor, with a charged species held by a held
//         contact or port): evaluated at φ̂ = 0 (its ion level is fixed by its held carrier and c0).
//      c) Intrinsic φ̂ anchors (in dipole-reduced coordinate φ̃ = φ̂ − Σ_{f<r} dipole_f / V_T):
//         - Tier 1a: every node of a region where a level-fixed charged species responds to φ̂, by
//           neutrality (bisection; node by node where profiled, graded, capped, or level-ramped;
//           a trial a non-ideal model refuses counts its ideal charge; where neutrality is out of
//           reach, in a strictly neutral region left unanchored and noted in s.unbalanced, else
//           toward less charge to at most ε V_T/(F L²) or where it stops falling, from the region's
//           boundary reference).
//         - Tier 1b: every node of a strictly neutral region under a held port's capacitance
//           (capA > 0, excluding end boxes whose φ a contact sets), by ρ + aσ = 0.
//         - Tier 1c: a held contact's end node where the end region doesn't respond — (V −
//           zeroCharge)/V_T for a species-linked 'pinned' contact, or neutrality at the contact's
//           outside levels for 'bulk' or a charged equilibrium link (so a c0 given for a
//           contact-fed carrier keeps the contact's gauge). Capacitive contacts and blocking
//           'pinned' contacts wait for the blocking-face tier below.
//         - Tier 2 (on regions neither responding, capped, nor directly clamped on a charged
//           stretch by a held 'bulk'/equilibrium contact): a held electrode on either face (metal
//           at its carrier's level, or a held port electrode at its window), balancing its first
//           reaction that depends on φ̂ (including neutral level-held participants, non-ideal
//           activities, profiles/grading at the electrode, and surface θ₀), or a held insertion
//           host across a face passing a charged species (by link or transfer reaction).
//         - Blocking capacitive/pinned metal face or contact (V − zeroCharge)/V_T: anchors a
//           charge-free dielectric directly (e.g. a MOS oxide), and otherwise serves only as a
//           fallback for chains with no Tier 1–2 anchor.
//      d) Outward propagation and unheld/interior conductors: within each contiguous non-φ-free
//         chain, φ̃ propagates constantly outward from the outer anchors; an all-c0 conducting chain
//         between two boundary anchors takes their flat symmetric midpoint (avoiding a fake drift
//         field across amount-fixed spectators), and otherwise φ̃ interpolates in ∫dx/σ (or ∫dx/ε
//         across insulators) between consecutive anchors, restoring face dipoles.
//         Unheld and interior conductors take their carrier's equilibrium level from whichever
//         adjacent anchored face(s) carry a reaction for it (the mean if both do), and may in turn
//         anchor unanchored neighbours; if a component still has no held anchor, blocking metal
//         faces and then floating electrodes serve as fallbacks, else φ̃ = 0.
//      A non-ideal material's composition comes from its statistics (materialAt).
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
//      where a capacitor's box needs counter-charge, neutralStart puts it on its Gauss row,
//      each closed amount kept and concentrations outside the projected boxes unchanged;
//      s.coldStart false. A solve that converges, a warm start and a carried state set it as
//      their state is.
//   3. An unheld insertion host (a φ-free material, not a conductor, with no charged species held
//      by a held contact or port): its φ̂ and every η shifted by z·Δ so that its first charged
//      species shared with an adjacent non-φ-free neighbour that the face between them passes (by a
//      link or a reaction taking it across; else its first shared; mean if both sides pass) starts
//      at that neighbour's level; any adjacent unheld conductor sharing its carrier across the face
//      takes the host's shifted carrier level. A held insertion host is never shifted off its held
//      contact.
//   4. A floating contact: facing a held contact, V = the held V + the difference between the two
//      terminal species' levels beside each (a bath without a reference by the φ beside it, its V
//      being its φ); else (anchored by a held port) directly its own level.
//   5. Electrode surfaces at their θ₀ (s._surfaceStart); concentrations computed.
//   6. Floating double layers start uncharged: a capacitive contact at zeroCharge + V_T·φ̂
//      beside it; a port with a capacitance at zeroCharge + its window's mean V_T·φ̂ (by
//      volume·area, over the boxes whose charge it holds: where φ is defined, and not an end box
//      whose φ its contact sets).
//   7. A floating port linked to its carrier, with no reactions or capacitance: at the carrier's
//      level across its window, relative to the held contact(s)' (mean where both are held).
//   8. Gates at α/(α + β) at the starting voltage across them.
//   9. A floating port with reactions and no capacitance: at the voltage where its reactions and
//      its conductance and exchange links pass the current it's driven at, or behind a
//      resistance (V_src − V)/R, by bisection; where none within ±1 kV does, at the V it had.
//
// Invariants (each says what a reviewer should be able to check here):
//   1. A species given by its amount starts at c0 (or its profile) at every node; at a face
//      equilibrated by step 2 (a transient's first step), its amount over the two nodes is kept.
//   2. A species given by its level starts at the level of what feeds it, as step 1(a) orders them.
//   3. Where a level-fixed charged species responds, each node starts neutral (Σ z c + the fixed
//      charge = 0, plus aσ/F under a held capacitance in a strictly neutral region) to the
//      bisection's tolerance, unless neutrality is out of reach (step 1(c): then its charge at most
//      ε V_T/(F L²), or the least φ̂ gives, to the bisection's tolerance); a strictly neutral
//      region given by amounts is neutral by its c0, checked (with its aσ under such a window,
//      which φ̂ then sets).
//   4. Steps 2 and 3 change no concentration but as they say: step 3 is a gauge (φ̂ and each η by
//      z·Δ) on unheld insertion hosts alone. Step 2 shares passing species across an equilibrium face
//      and corrects the capacitor boxes' charge at fixed closed amounts.
//   5. Each floating terminal's V is set by the last of steps 4, 6, 7, 9 that applies to it;
//      held terminals' V are untouched.
//   6. Only u, uLo, termV, s.unbalanced, s.coldStart and what computeConcentrations derives are written; the
//      result depends on the model, every terminal's voltage and s.time, deterministically, and is
//      mirror-symmetric on mirror-image definitions.
// Known, and not an invariant (see ROADMAP): a floating metal contact behind a resistance (or
// current-driven) starts level with its held neighbour (Solver's constructor), which in an unlike
// cell is its open-circuit voltage away (a closed Cu | Zn cell's zinc 1.13 V from rest: the steady
// solve's held stage lays the start out again at the source, src/steady.js; balancing a floating
// contact's face reactions needs a face-rate evaluator covering Stern, diffuse, MHC and surfaces).
// Where ε > 0, a held port's capacitance through its window is left out of the neutrality, so the
// first step moves φ̂ to discharge it and the ions charge it again (a gated Li⁺ bar with εr = 1:
// 0.28 C/m² cold, 0.0035 after 1 ps, back to 0.28 by diffusion; Astra's third pass).
// Tests: every test (each device starts here); test/start.test.js (anchor-outward contract,
// mirror symmetry of the cold state, and the resolved Known items); test/transient.test.js on what
// a start must be; test/round7-solve.test.js (φ̂ behind a bath or a blocking metal, a floating
// port's electrode); test/neon1.test.js; test/round9-mirror.test.js.

import { gateRate } from './boundaries.js';
import { EPS0, FARADAY } from './constants.js';
import { fixedAt, gradedAt, isGraded, profileAt, sourceAt } from './device.js';
import { diffuseLayer } from './diffuse.js';
import { SolverError } from './errors.js';
import { gaussRows, solveNeutral, strictlyNeutralNodes } from './neutral.js';

/**
 * Cold start. Each species takes its region's c0 where one is given (for a spectator, it fixes
 * the conserved amount), else its contact's or port's level. Intrinsic φ̂ anchors (neutrality in
 * responding regions, held port capacitances in strictly neutral regions, held contacts, held
 * electrodes and held insertion hosts) are found symmetrically and propagated outward (constant
 * outside, ohmic ∫dx/σ or dielectric ∫dx/ε between consecutive anchors, restoring face dipoles).
 */
export function initFromComposition(s) {
  const { model, n, M, u, z } = s;
  const { regions, species, interfaces, materials, ports } = model;
  const grid = model.grid;
  const R = regions.length, N = grid.nNodes;
  u.fill(0);
  s.uLo.fill(0);
  s.unbalanced = [];

  const isHeldTerm = (k) => !s.floating.includes(k);
  // The outside level a contact's link offers species i: held there, or exchanged with at a rate
  // (conductance, exchange), whose level a start can take as well.
  const contactLevel = (side, i) => {
    const l = model.contacts[side].species[i];
    if (l.type === 'blocked') return NaN;
    if (l.type !== 'conductance' && l.type !== 'exchange') return s.contactEta(side, i);
    return z[i] === 0 ? l.mu / model.RT : (z[i] * (s.termV[side === 'left' ? 0 : 1] + l.offset)) / s.VT;
  };
  const contactHeld = (side, i) => z[i] === 0 || isHeldTerm(side === 'left' ? 0 : 1);
  const portHeld = (p, i) => z[i] === 0 || isHeldTerm(2 + p);

  // A pinned contact's densities, cRef e^(η − μ°/RT − z φ̂_edge), φ̂_edge = (V − zeroCharge)/V_T:
  // beyond e^700, past what a double holds, the definition can't be solved, and says so here
  // (Neon's review, 1.2).
  for (const [side, k, g] of [['left', 0, 0], ['right', 1, N - 1]]) {
    const ph = model.contacts[side].phi;
    if (ph.type !== 'pinned' || !isHeldTerm(k)) continue;
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

  // Cumulative interface dipoles (in V_T): φ̃ = φ̂ − cumDip[r] is continuous across every face
  // under constant propagation and reverses cleanly under mirror reflection.
  const cumDip = new Float64Array(R);
  for (let r = 1; r < R; r++) cumDip[r] = cumDip[r - 1] + interfaces[r - 1].dipole / s.VT;

  const { capA, capB } = heldCapacitance(s);

  // --- Pass 1: species modes, cFix, and feeding sources (order-independent, mirror-symmetric) ---
  const regMode = Array.from({ length: R }, () => new Int8Array(n));
  const regEta = Array.from({ length: R }, () => new Float64Array(n));
  const regEtaNode = Array.from({ length: R }, () => new Array(n).fill(null));
  const regCFix = Array.from({ length: R }, () => new Float64Array(n));
  const regHeldLevel = Array.from({ length: R }, () => new Uint8Array(n));
  const interiorMetal = new Uint8Array(R);
  const heldMetal = new Uint8Array(R);
  const hostHeld = new Uint8Array(R);
  const responds = new Uint8Array(R);
  const capped = new Uint8Array(R);
  const bulkFed = new Uint8Array(R);
  const regRamped = new Uint8Array(R);

  // Classify each stretch's feeding sources once so multi-region stretches share one rule.
  const stretchSpec = s.stretches.map((st) => {
    const i = st.species;
    const sources = [];
    if (st.leftOpen) {
      const e = contactLevel('left', i);
      if (Number.isFinite(e)) sources.push({ kind: 'contact', side: 'left', g: 0, held: contactHeld('left', i), bulk: model.contacts.left.phi.type === 'bulk', eta: e });
    }
    if (st.rightOpen) {
      const e = contactLevel('right', i);
      if (Number.isFinite(e)) sources.push({ kind: 'contact', side: 'right', g: N - 1, held: contactHeld('right', i), bulk: model.contacts.right.phi.type === 'bulk', eta: e });
    }
    for (const p of st.linked) {
      const port = ports[p], e = s.portEta(port, i);
      if (Number.isFinite(e)) {
        const gMid = port.nodes[port.nodes.length >> 1];
        sources.push({ kind: 'port', p, g: gMid, x: 0.5 * (grid.x[port.nodes[0]] + grid.x[port.nodes[port.nodes.length - 1]]), held: portHeld(p, i), bulk: false, eta: e });
      }
    }
    if (sources.length === 0) return { kind: 'none', held: false };
    const anyHeld = sources.some((src) => src.held);
    const active = anyHeld ? sources.filter((src) => src.held) : sources;
    const cts = active.filter((src) => src.kind === 'contact');
    if (cts.length === 1) return { kind: 'const', eta: cts[0].eta, held: anyHeld };
    if (cts.length === 2) {
      if (cts[0].bulk !== cts[1].bulk) {
        const winner = cts[0].bulk ? cts[0] : cts[1];
        return { kind: 'const', eta: winner.eta, held: anyHeld };
      }
      if (cts[0].eta === cts[1].eta) return { kind: 'const', eta: cts[0].eta, held: anyHeld };
      // Across a pure conductor stretch, or a single uniformly-doped region between two bulk
      // contacts (e.g. a resistor), ohmic ramp is exact; across a multi-region or junction stretch
      // (which may contain blocking pn junctions or floating bases), the symmetric midpoint avoids
      // imposing a huge artificial bulk current gradient while staying mirror-symmetric.
      let allMetal = true;
      for (let r = st.regions[0]; r <= st.regions[1]; r++) if (!materials[regions[r].material].conductor) { allMetal = false; break; }
      const singleBulk = st.regions[0] === st.regions[1] && cts[0].bulk && cts[1].bulk && !regions[st.regions[0]].fixedProfile;
      if (allMetal || singleBulk) return { kind: 'ramp2', gL: st.nodes[0], gR: st.nodes[1], etaL: cts[0].eta, etaR: cts[1].eta, held: anyHeld };
      return { kind: 'const', eta: 0.5 * (cts[0].eta + cts[1].eta), held: anyHeld };
    }
    // Ports only: merge ports at the same x (commutative), then constant or piecewise ramp.
    const byX = active.slice().sort((a, b) => a.x - b.x);
    const pts = [];
    for (const src of byX) {
      if (pts.length > 0 && Math.abs(pts[pts.length - 1].x - src.x) <= 1e-12 * Math.max(1, Math.abs(src.x))) {
        const prev = pts[pts.length - 1];
        prev.eta = (prev.eta * prev.cnt + src.eta) / (prev.cnt + 1);
        prev.cnt++;
      } else {
        pts.push({ g: src.g, x: src.x, eta: src.eta, cnt: 1 });
      }
    }
    if (pts.length === 1 || pts.every((p) => p.eta === pts[0].eta)) return { kind: 'const', eta: pts[0].eta, held: anyHeld };
    return { kind: 'rampPts', pts, held: anyHeld };
  });

  // First determine mode (amount vs level) and validate missing c0 / strict neutrality.
  for (let r = 0; r < R; r++) {
    const reg = regions[r], base = materials[reg.material];
    if (base.epsr === 0 && !base.conductor && !base.phiFree) {
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) if (capA[g] > 0) capped[r] = 1;
    }
    for (let i = 0; i < n; i++) {
      if (!base.present[i]) continue;
      const stIdx = s.stretchOf[r * n + i], st = s.stretches[stIdx], spec = stretchSpec[stIdx];
      if (reg.c0Profile[i] || reg.c0[i] > 0) {
        regMode[r][i] = 2;
        regCFix[r][i] = reg.c0Profile[i] ? profileAt(reg.c0Profile[i], grid.x[grid.regionStart[r]]) : reg.c0[i];
      } else if (!st.contactFed && base.conductor) {
        regMode[r][i] = 1;
        interiorMetal[r] = 1;
      } else if (!st.contactFed) {
        if (!(reg.c0[i] > 0)) {
          throw new SolverError(
            `regions[${r}].c0.${species[i].name}: a species that doesn't reach a contact needs its initial concentration` +
              (st.connected ? ' (reactions make and consume it, so this is only where the solve starts, not an amount it keeps: a tiny value will do)' : ', which fixes the amount it conserves'),
          );
        }
        regMode[r][i] = 2;
        regCFix[r][i] = reg.c0[i];
      } else if (spec.kind === 'none') {
        if (!(reg.c0[i] > 0)) {
          throw new SolverError(
            `regions[${r}].c0.${species[i].name}: a species fed only through reactions needs its initial concentration`,
          );
        }
        regMode[r][i] = 2;
        regCFix[r][i] = reg.c0[i];
      } else {
        regMode[r][i] = 1;
        if (spec.kind === 'const') regEta[r][i] = spec.eta;
        if (spec.held) regHeldLevel[r][i] = 1;
      }
    }
    if (base.conductor) {
      const ci = base.conductor.i;
      if (!interiorMetal[r] && regHeldLevel[r][ci]) heldMetal[r] = 1;
    } else if (base.phiFree) {
      for (let i = 0; i < n; i++) if (regMode[r][i] === 1 && z[i] !== 0 && regHeldLevel[r][i]) hostHeld[r] = 1;
    } else {
      for (let i = 0; i < n; i++) if (regMode[r][i] === 1 && z[i] !== 0) responds[r] = 1;
    }
  }

  // Strictly neutral regions with every charged species' amount given: check neutrality of c0
  // outside any held port capacitance's window (exact same check and error as before).
  for (let r = 0; r < R; r++) {
    const reg = regions[r], base = materials[reg.material];
    if (responds[r] || base.conductor || base.phiFree || base.epsr !== 0) continue;
    const profiled = reg.c0Profile.some((p) => p !== null) || reg.fixedProfile !== null;
    const checkAt = (g, where) => {
      let q = fixedAt(reg, grid.x[g]) / FARADAY, scale = Math.abs(q);
      for (let i = 0; i < n; i++) {
        if (regMode[r][i] !== 2 || z[i] === 0) continue;
        const ci = reg.c0Profile[i] ? profileAt(reg.c0Profile[i], grid.x[g]) : regCFix[r][i];
        q += z[i] * ci;
        scale += Math.abs(z[i] * ci);
      }
      if (Math.abs(q) > 1e-9 * scale) {
        throw new SolverError(
          `regions[${r}] (${reg.name}): its initial composition carries a net charge of ${q.toPrecision(3)} mol/m³${where} (Σ z·c0, with any fixed charge), ` +
            'but ε = 0 makes it strictly neutral; adjust c0 so that it balances',
        );
      }
    };
    if (!profiled && !capped[r]) checkAt(grid.regionStart[r], '');
    else {
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        if (!(capA[g] > 0)) checkAt(g, ` at x = ${grid.x[g].toPrecision(4)} m`);
      }
    }
  }

  // --- Neutrality helper (ideal + non-ideal, with clamped-exponent and out-of-reach handling) ---
  const zetaScratch = new Float64Array(n), ccScratch = new Float64Array(n);
  const solveNeutralAt = (r, mat, bg, qFixed, etaArr, modeArr, cFixArr, aC, aCV, phRef, noteUnbalanced) => {
    // Pre-sum the φ-independent fixed + amount-fixed (mode 2) charge before adding level-fixed
    // carriers so exact neutrality between fixedCharge and c0 (e.g. perovskite qF + z·c0 = 0)
    // doesn't lose ~12 digits of float64 precision when added to tiny carrier densities.
    let qFix = qFixed;
    for (let i = 0; i < n; i++) if (mat.present[i] && modeArr[i] === 2) qFix += z[i] * cFixArr[i];
    let clamped = false;
    const charge = (ph) => {
      if (!mat.ideal) {
        let refused = false;
        try {
          materialAt(s, mat, bg, etaArr, modeArr, cFixArr, ph, zetaScratch, ccScratch);
        } catch (err) {
          if (!(err instanceof SolverError)) throw err;
          refused = true;
        }
        if (!refused) {
          let q = qFix + (aCV - aC * s.VT * ph) / FARADAY;
          for (let i = 0; i < n; i++) if (mat.present[i] && modeArr[i] !== 2) q += z[i] * ccScratch[i];
          return q;
        }
      }
      let q = qFix + (aCV - aC * s.VT * ph) / FARADAY;
      clamped = false;
      for (let i = 0; i < n; i++) {
        if (modeArr[i] === 1 && z[i] !== 0) {
          const raw = etaArr[i] - mat.mu0[i] / model.RT - z[i] * ph, ex = Math.min(700, Math.max(-700, raw));
          if (ex !== raw) clamped = true;
          q += z[i] * mat.cRef[i] * Math.exp(ex);
        }
      }
      return q;
    };
    let lo = -1, hi = 1;
    while (charge(lo) < 0 && lo > -1e4) lo -= 2 * (hi - lo);
    while (charge(hi) > 0 && hi < 1e4) hi += 2 * (hi - lo);
    if (charge(lo) < 0 || charge(hi) > 0) {
      if (mat.epsr === 0) {
        if (noteUnbalanced && !s.unbalanced.includes(r)) s.unbalanced.push(r);
        return { ph: phRef, ok: false };
      }
      const L = grid.x[grid.regionEnd[r]] - grid.x[grid.regionStart[r]];
      const most = (mat.epsr * EPS0 * model.RT) / (FARADAY * FARADAY * L * L);
      const q0 = charge(phRef), dir = q0 > 0 ? 1 : -1;
      if (!(Math.abs(q0) > most)) return { ph: phRef, ok: true };
      let a = phRef, qa = Math.abs(q0);
      for (let step = 1; step < 1e4; step *= 2) {
        let b = a + dir * step;
        const qb = Math.abs(charge(b));
        if (qb <= most) {
          for (let it = 0; it < 100 && Math.abs(b - a) > 1e-9; it++) {
            const m = 0.5 * (a + b);
            if (Math.abs(charge(m)) > most) a = m;
            else b = m;
          }
          return { ph: b, ok: true };
        }
        if (!(qb < qa * (1 - 1e-9)) && !clamped) return { ph: a, ok: true };
        [a, qa] = [b, qb];
      }
      return { ph: a, ok: true };
    }
    const tol = 4 * Number.EPSILON * Math.max(1, Math.abs(lo), Math.abs(hi));
    for (let it = 0; it < 200 && hi - lo > tol; it++) {
      const m = 0.5 * (lo + hi), qm = charge(m);
      if (qm === 0) return { ph: m, ok: true };
      if (qm > 0) lo = m;
      else hi = m;
    }
    return { ph: 0.5 * (lo + hi), ok: true };
  };

  // --- Cell conductivity σ (S/m) and permittivity ε (F/m) for ohmic / dielectric interpolation ---
  // Shared across charged species so neutral combinations (e⁻ + h⁺, Na⁺ + Cl⁻) stay bounded.
  const segSigma = new Float64Array(Math.max(0, N - 1));
  const segEps = new Float64Array(Math.max(0, N - 1));
  const etaMidBuf = new Float64Array(n), cFixMidBuf = new Float64Array(n);
  const f2OverRT = (FARADAY * FARADAY) / model.RT;
  for (let r = 0; r < R; r++) {
    const reg = regions[r], base = materials[reg.material], graded = isGraded(reg);
    const mode = regMode[r];
    // For an undoped, un-c0 responding region (e.g. intrinsic i-Si in a p-i-n diode), estimate its
    // intrinsic carrier scale at the midpoint of its feeding levels once so σ reflects n_i, not cRef.
    let intrinsicC = null;
    if (!base.conductor && !base.phiFree && responds[r]) {
      let hasZeroFixed = !reg.fixedProfile && reg.fixedCharge === 0;
      if (hasZeroFixed) {
        for (let i = 0; i < n; i++) {
          if (mode[i] === 1) {
            const sp = stretchSpec[s.stretchOf[r * n + i]];
            etaMidBuf[i] = sp.kind === 'const' ? sp.eta : sp.kind === 'ramp2' ? 0.5 * (sp.etaL + sp.etaR) : sp.kind === 'rampPts' ? 0.5 * (sp.pts[0].eta + sp.pts[sp.pts.length - 1].eta) : 0;
          } else if (mode[i] === 2) {
            cFixMidBuf[i] = reg.c0Profile[i] ? profileAt(reg.c0Profile[i], 0.5 * (grid.x[grid.regionStart[r]] + grid.x[grid.regionEnd[r]])) : regCFix[r][i];
          }
        }
        const ph0 = solveNeutralAt(r, base, reg.background, 0, etaMidBuf, mode, cFixMidBuf, 0, 0, 0, false).ph;
        intrinsicC = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          if (mode[i] === 1 && z[i] !== 0) {
            intrinsicC[i] = base.cRef[i] * Math.exp(Math.min(700, Math.max(-700, etaMidBuf[i] - base.mu0[i] / model.RT - z[i] * ph0)));
          }
        }
      }
    }
    for (let g = grid.regionStart[r]; g < grid.regionEnd[r]; g++) {
      const xMid = 0.5 * (grid.x[g] + grid.x[g + 1]);
      const mat = graded ? gradedAt(base, reg, xMid) : base;
      segEps[g] = mat.epsr * EPS0;
      if (mat.conductor) {
        segSigma[g] = mat.conductor.conductivity;
        continue;
      }
      const qFix = fixedAt(reg, xMid) / FARADAY;
      let nCounter = 0;
      if (qFix !== 0) {
        for (let i = 0; i < n; i++) if (mode[i] === 1 && z[i] !== 0 && mat.D[i] > 0 && z[i] * qFix < 0) nCounter++;
      }
      let sig = 0;
      for (let i = 0; i < n; i++) {
        if (!mat.present[i] || z[i] === 0 || !(mat.D[i] > 0)) continue;
        let ci = 0;
        if (mode[i] === 2) ci = reg.c0Profile[i] ? profileAt(reg.c0Profile[i], xMid) : regCFix[r][i];
        else if (mode[i] === 1) {
          if (qFix !== 0 && z[i] * qFix < 0 && nCounter > 0) ci = Math.abs(qFix) / (Math.abs(z[i]) * nCounter);
          else if (intrinsicC) ci = intrinsicC[i];
          else ci = mat.cRef[i];
        }
        sig += z[i] * z[i] * mat.D[i] * ci;
      }
      segSigma[g] = f2OverRT * sig;
    }
  }

  // Cumulative interpolation fraction [0, 1] from node g0 to g1 along ∫dx/σ (or ∫dx/ε across
  // insulators, or ∫dx/D for a neutral species). Mirror-symmetric by construction.
  const interpFractions = (g0, g1, neutralSpecies = -1) => {
    const len = g1 - g0 + 1, frac = new Float64Array(len);
    if (len <= 1) return frac;
    const dw = new Float64Array(len - 1);
    if (neutralSpecies >= 0) {
      for (let k = 0; k < len - 1; k++) {
        const g = g0 + k, r = grid.nodeRegion[g];
        if (r !== grid.nodeRegion[g + 1]) continue; // interface: dx = 0
        const dx = grid.x[g + 1] - grid.x[g];
        const reg = regions[r], base = materials[reg.material];
        const mat = isGraded(reg) ? gradedAt(base, reg, 0.5 * (grid.x[g] + grid.x[g + 1])) : base;
        const D = mat.D[neutralSpecies];
        dw[k] = D > 0 ? dx / D : dx;
      }
    } else {
      let hasInsulator = false;
      for (let k = 0; k < len - 1; k++) {
        const g = g0 + k;
        if (grid.nodeRegion[g] !== grid.nodeRegion[g + 1]) continue;
        if (!(segSigma[g] > 0) && segEps[g] > 0) { hasInsulator = true; break; }
      }
      for (let k = 0; k < len - 1; k++) {
        const g = g0 + k;
        if (grid.nodeRegion[g] !== grid.nodeRegion[g + 1]) continue;
        const dx = grid.x[g + 1] - grid.x[g];
        if (hasInsulator) dw[k] = !(segSigma[g] > 0) && segEps[g] > 0 ? dx / segEps[g] : 0;
        else if (segSigma[g] > 0) dw[k] = dx / segSigma[g];
      }
    }
    let total = 0;
    for (let k = 0; k < len - 1; k++) total += dw[k];
    if (!(total > 0)) {
      // Fallback: geometric distance (or uniform steps if zero length)
      for (let k = 0; k < len - 1; k++) {
        const g = g0 + k;
        dw[k] = grid.nodeRegion[g] === grid.nodeRegion[g + 1] ? grid.x[g + 1] - grid.x[g] : 0;
        total += dw[k];
      }
    }
    if (!(total > 0)) {
      for (let k = 0; k < len; k++) frac[k] = k / (len - 1);
      return frac;
    }
    let acc = 0;
    for (let k = 0; k < len - 1; k++) {
      acc += dw[k];
      frac[k + 1] = acc / total;
    }
    frac[len - 1] = 1;
    return frac;
  };

  // Fill ramped levels for stretches fed by two equal-authority contacts (or multiple ports).
  for (let stIdx = 0; stIdx < s.stretches.length; stIdx++) {
    const spec = stretchSpec[stIdx];
    if (spec.kind !== 'ramp2' && spec.kind !== 'rampPts') continue;
    const st = s.stretches[stIdx], i = st.species, neutralSp = z[i] === 0 ? i : -1;
    const gStart = st.nodes[0], gEnd = st.nodes[1];
    const full = new Float64Array(gEnd - gStart + 1);
    if (spec.kind === 'ramp2') {
      const f = interpFractions(gStart, gEnd, neutralSp);
      for (let k = 0; k < full.length; k++) full[k] = (1 - f[k]) * spec.etaL + f[k] * spec.etaR;
    } else {
      const pts = spec.pts;
      for (let g = gStart; g <= pts[0].g; g++) full[g - gStart] = pts[0].eta;
      for (let g = pts[pts.length - 1].g; g <= gEnd; g++) full[g - gStart] = pts[pts.length - 1].eta;
      for (let p = 0; p + 1 < pts.length; p++) {
        const a = pts[p], b = pts[p + 1], f = interpFractions(a.g, b.g, neutralSp);
        for (let g = a.g; g <= b.g; g++) full[g - gStart] = (1 - f[g - a.g]) * a.eta + f[g - a.g] * b.eta;
      }
    }
    for (let r = st.regions[0]; r <= st.regions[1]; r++) {
      if (regMode[r][i] !== 1) continue;
      const a = grid.regionStart[r], b = grid.regionEnd[r];
      regEtaNode[r][i] = full.slice(a - gStart, b - gStart + 1);
      regEta[r][i] = regEtaNode[r][i][0];
      regRamped[r] = 1;
    }
  }

  const etaAt = (r, g, i) => (regEtaNode[r][i] ? regEtaNode[r][i][g - grid.regionStart[r]] : regEta[r][i]);
  const fillNodeInputs = (r, g, etaOut, cFixOut) => {
    const reg = regions[r], mode = regMode[r];
    for (let i = 0; i < n; i++) {
      if (mode[i] === 1) etaOut[i] = etaAt(r, g, i);
      else if (mode[i] === 2) cFixOut[i] = reg.c0Profile[i] ? profileAt(reg.c0Profile[i], grid.x[g]) : regCFix[r][i];
    }
  };

  // --- Pass 1b: populate contact-fed conductors and held insertion hosts ---
  const etaBuf = new Float64Array(n), cFixBuf = new Float64Array(n);
  const writeNode = (r, g, ph, mat, mode, etaArr, cFixArr) => {
    if (!mat.ideal) materialAt(s, mat, regions[r].background, etaArr, mode, cFixArr, ph, zetaScratch, ccScratch);
    const b = s.blockOfNode[g] * M;
    u[b] = mat.conductor ? 0 : ph;
    for (let i = 0; i < n; i++) {
      if (mode[i] === 1) u[b + 1 + i] = etaArr[i];
      else if (mode[i] === 2) {
        const zt = mat.ideal ? Math.log(cFixArr[i] / mat.cRef[i]) : zetaScratch[i];
        u[b + 1 + i] = zt + mat.mu0[i] / model.RT + z[i] * ph;
      }
    }
  };
  for (let r = 0; r < R; r++) {
    const reg = regions[r], base = materials[reg.material];
    if (base.conductor && !interiorMetal[r]) {
      const ci = base.conductor.i;
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const b = s.blockOfNode[g] * M;
        u[b] = 0;
        u[b + 1 + ci] = etaAt(r, g, ci);
      }
    } else if (base.phiFree && !base.conductor && hostHeld[r]) {
      const graded = isGraded(reg);
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const mat = graded ? gradedAt(base, reg, grid.x[g]) : base;
        fillNodeInputs(r, g, etaBuf, cFixBuf);
        writeNode(r, g, 0, mat, regMode[r], etaBuf, cFixBuf);
      }
    }
  }

  // --- Electrode reaction balance helper (now handles neutral mode-1 participants too) ---
  const balancedAt = (r, x, gForEta, reactions, surface, inside, levelOf) => {
    const reg = regions[r], base = materials[reg.material], mode = regMode[r];
    const mat = isGraded(reg) ? gradedAt(base, reg, x) : base;
    for (let i = 0; i < n; i++) {
      if (mode[i] === 1) etaBuf[i] = etaAt(r, gForEta, i);
      else if (mode[i] === 2) cFixBuf[i] = reg.c0Profile[i] ? profileAt(reg.c0Profile[i], x) : regCFix[r][i];
    }
    if (!mat.ideal) materialAt(s, mat, reg.background, etaBuf, mode, cFixBuf, 0, zetaScratch, ccScratch);
    const activity = (i) => (mat.ideal ? Math.log(cFixBuf[i] / mat.cRef[i]) : zetaScratch[i]);
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
        } else if (mode[p.i] === 1 && z[p.i] === 0) {
          a -= p.nu * etaBuf[p.i];
        } else slope = NaN;
      }
      if (slope !== 0 && Number.isFinite(slope) && Number.isFinite(a)) return a / slope;
    }
    return NaN;
  };

  // --- Pass 2: intrinsic φ̃ anchors (completely symmetric) ---
  // Tier 1: responds (neutrality); Tier 2: capped (ρ + aσ = 0); Tier 3: held contact;
  // Tier 4: held electrode / held insertion host; Tier 5: fallback floating electrode.
  const ancTier = new Int8Array(N);
  const ancPhiTilde = new Float64Array(N);
  const ancCount = new Int32Array(N);
  const setAnchor = (g, tier, phiHatVal) => {
    if (!Number.isFinite(phiHatVal)) return;
    const pt = phiHatVal - cumDip[grid.nodeRegion[g]];
    if (ancTier[g] === 0 || tier < ancTier[g]) {
      ancTier[g] = tier;
      ancPhiTilde[g] = pt;
      ancCount[g] = 1;
    } else if (tier === ancTier[g]) {
      ancPhiTilde[g] = (ancPhiTilde[g] * ancCount[g] + pt) / (ancCount[g] + 1);
      ancCount[g]++;
    }
  };

  // Direct Dirichlet bulk/equilibrium clamp check on each contact.
  const directClamp = (side, k) => {
    if (!isHeldTerm(k)) return false;
    const ct = model.contacts[side];
    return ct.phi.type === 'bulk' || ct.species.some((l, i) => z[i] !== 0 && l.type === 'equilibrium');
  };
  const clampL = directClamp('left', 0), clampR = directClamp('right', 1);

  // Compute each held contact's own boundary φ̂ (pinned/capacitive, or bulk/equilibrium neutrality).
  const modeTrial = new Int8Array(n), etaTrial = new Float64Array(n);
  const contactAnchorPhi = (side, k, gEnd, rEnd) => {
    const reg = regions[rEnd], base = materials[reg.material];
    if (!isHeldTerm(k) || base.phiFree) return NaN;
    const ct = model.contacts[side];
    if (ct.phi.type === 'pinned' || ct.phi.type === 'capacitive') return (s.termV[k] - ct.phi.zeroCharge) / s.VT;
    if (!directClamp(side, k)) return NaN;
    const mat = isGraded(reg) ? gradedAt(base, reg, grid.x[gEnd]) : base;
    fillNodeInputs(rEnd, gEnd, etaBuf, cFixBuf);
    let anyLinked = false;
    for (let i = 0; i < n; i++) {
      const cl = contactLevel(side, i);
      if (mat.present[i] && z[i] !== 0 && Number.isFinite(cl)) {
        modeTrial[i] = 1;
        etaTrial[i] = cl;
        anyLinked = true;
      } else {
        const sp = mat.present[i] && z[i] !== 0 ? stretchSpec[s.stretchOf[rEnd * n + i]] : null;
        if (sp && sp.held && sp.kind === 'const') {
          modeTrial[i] = 1;
          etaTrial[i] = sp.eta;
        } else {
          modeTrial[i] = regMode[rEnd][i];
          etaTrial[i] = etaBuf[i];
        }
      }
    }
    if (!anyLinked) return NaN;
    const qF = fixedAt(reg, grid.x[gEnd]) / FARADAY;
    // A bath without a reference (V = φ_bath) on a region with no fixed charge has exact φ̂ = V/V_T
    // whenever the bath's own levels are already neutral there (avoiding ~1e-15 bisection noise).
    if (ct.she && !ct.she.reference && qF === 0) {
      const ph0 = s.termV[k] / s.VT;
      let q0 = 0, sc = 0, ok0 = true;
      if (!mat.ideal) {
        try {
          materialAt(s, mat, reg.background, etaTrial, modeTrial, cFixBuf, ph0, zetaScratch, ccScratch);
          for (let i = 0; i < n; i++) if (mat.present[i]) { q0 += z[i] * ccScratch[i]; sc += Math.abs(z[i] * ccScratch[i]); }
        } catch (err) {
          if (!(err instanceof SolverError)) throw err;
          ok0 = false;
        }
      } else {
        for (let i = 0; i < n; i++) {
          if (modeTrial[i] === 2) { q0 += z[i] * cFixBuf[i]; sc += Math.abs(z[i] * cFixBuf[i]); }
          else if (modeTrial[i] === 1 && z[i] !== 0) {
            const ex = etaTrial[i] - mat.mu0[i] / model.RT - z[i] * ph0;
            if (Math.abs(ex) < 700) {
              const ci = mat.cRef[i] * Math.exp(ex);
              q0 += z[i] * ci;
              sc += Math.abs(z[i] * ci);
            } else ok0 = false;
          }
        }
      }
      if (ok0 && Number.isFinite(sc) && Math.abs(q0) <= 1e-11 * Math.max(1, sc)) return ph0;
    }
    const res = solveNeutralAt(rEnd, mat, reg.background, qF, etaTrial, modeTrial, cFixBuf, 0, 0, 0, false);
    if (res.ok) return res.ph;
    if (!mat.ideal) materialAt(s, mat, reg.background, etaBuf, regMode[rEnd], cFixBuf, 0, zetaScratch, ccScratch);
    const ti = ct.terminal !== null && modeTrial[ct.terminal] === 1 && regMode[rEnd][ct.terminal] === 2
      ? ct.terminal
      : z.findIndex((zi, i) => zi !== 0 && modeTrial[i] === 1 && regMode[rEnd][i] === 2);
    if (ti < 0) return NaN;
    const zt = mat.ideal ? Math.log(cFixBuf[ti] / mat.cRef[ti]) : zetaScratch[ti];
    return (etaTrial[ti] - mat.mu0[ti] / model.RT - zt) / z[ti];
  };
  const endPhiL = contactAnchorPhi('left', 0, 0, 0);
  const endPhiR = contactAnchorPhi('right', 1, N - 1, R - 1);

  // Symmetric boundary reference for out-of-reach neutrality (hole-only layers, or a resistor with
  // c0 given on its majority carrier alone): uses held contact φ̂ when available.
  const refPhiAt = (g) => {
    const d = cumDip[grid.nodeRegion[g]];
    if (Number.isFinite(endPhiL) && Number.isFinite(endPhiR)) {
      const f = grid.x[N - 1] > grid.x[0] ? (grid.x[g] - grid.x[0]) / (grid.x[N - 1] - grid.x[0]) : 0.5;
      return (1 - f) * endPhiL + f * (endPhiR - cumDip[R - 1]) + d;
    }
    if (Number.isFinite(endPhiL)) return endPhiL + d;
    if (Number.isFinite(endPhiR)) return endPhiR - cumDip[R - 1] + d;
    return d;
  };

  for (let r = 0; r < R; r++) {
    const base = materials[regions[r].material];
    if (base.phiFree) continue;
    for (let i = 0; i < n; i++) {
      if (!base.present[i] || z[i] === 0) continue;
      const st = s.stretches[s.stretchOf[r * n + i]];
      if ((st.leftOpen && clampL) || (st.rightOpen && clampR)) { bulkFed[r] = 1; break; }
    }
  }

  // Tier 1a (responds) and Tier 1b (capped in strictly neutral regions).
  for (let r = 0; r < R; r++) {
    const reg = regions[r], base = materials[reg.material];
    if (base.phiFree) continue;
    if (responds[r]) {
      const profiled = reg.c0Profile.some((p) => p !== null) || reg.fixedProfile !== null;
      const graded = isGraded(reg);
      const vary = profiled || graded || capped[r] || regRamped[r];
      if (!vary) {
        const g0 = grid.regionStart[r], g1 = grid.regionEnd[r];
        fillNodeInputs(r, g0, etaBuf, cFixBuf);
        const qF = fixedAt(reg, grid.x[g0]) / FARADAY;
        const res = solveNeutralAt(r, base, reg.background, qF, etaBuf, regMode[r], cFixBuf, 0, 0, refPhiAt((g0 + g1) >> 1), true);
        if (res.ok) for (let g = g0; g <= g1; g++) setAnchor(g, 1, res.ph);
      } else {
        let phGuess = refPhiAt(grid.regionStart[r]);
        for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
          const mat = graded ? gradedAt(base, reg, grid.x[g]) : base;
          fillNodeInputs(r, g, etaBuf, cFixBuf);
          const qF = fixedAt(reg, grid.x[g]) / FARADAY;
          const aC = capped[r] ? capA[g] : 0, aCV = capped[r] ? capB[g] : 0;
          const res = solveNeutralAt(r, mat, reg.background, qF, etaBuf, regMode[r], cFixBuf, aC, aCV, phGuess, true);
          if (res.ok) {
            phGuess = res.ph;
            setAnchor(g, 1, res.ph);
          }
        }
      }
    } else if (capped[r]) {
      const graded = isGraded(reg);
      let phGuess = refPhiAt(grid.regionStart[r]);
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        if (!(capA[g] > 0)) continue;
        const mat = graded ? gradedAt(base, reg, grid.x[g]) : base;
        fillNodeInputs(r, g, etaBuf, cFixBuf);
        const qF = fixedAt(reg, grid.x[g]) / FARADAY;
        const res = solveNeutralAt(r, mat, reg.background, qF, etaBuf, regMode[r], cFixBuf, capA[g], capB[g], phGuess, false);
        if (res.ok) {
          phGuess = res.ph;
          setAnchor(g, 2, res.ph);
        }
      }
    }
  }

  // Tier 1c (held contacts on non-responding, non-φ-free end regions; capacitive contacts and
  // blocking pinned contacts wait for Tier 5 so they never impose an artificial bulk ramp across a
  // conducting electrolyte anchored by a bath or electrode).
  const isBlockingCapPin = (side) => {
    const ct = model.contacts[side];
    if (ct.phi.type === 'capacitive') return true;
    if (ct.phi.type === 'pinned') return !ct.species.some((l, i) => z[i] !== 0 && l.type !== 'blocked');
    return false;
  };
  if (!responds[0] && Number.isFinite(endPhiL) && !isBlockingCapPin('left')) setAnchor(0, 3, endPhiL);
  if (!responds[R - 1] && Number.isFinite(endPhiR) && !isBlockingCapPin('right')) setAnchor(N - 1, 3, endPhiR);

  // Helper: anchor region r from a held insertion host at nr (r - 1 or r + 1).
  const anchorFromHost = (r, nr, gIn, gOut) => {
    const reg = regions[r], base = materials[reg.material];
    const other = materials[regions[nr].material];
    const itf = interfaces[Math.min(r, nr)];
    const shared = (j) => z[j] !== 0 && base.present[j] && other.present[j] && regMode[r][j] === 2;
    const passesLink = (j) => itf.links[j].type !== 'blocked';
    const passesRx = (j) => itf.reactions.some((rx) => rx.part.some((p) => p.i === j && p.side === 0) && rx.part.some((p) => p.i === j && p.side === 1));
    // Prefer a direct link, else a face reaction (which balancedAt handles with any fixedA), else first shared.
    let i = z.findIndex((_, j) => shared(j) && passesLink(j));
    if (i >= 0) {
      const mat = isGraded(reg) ? gradedAt(base, reg, grid.x[gIn]) : base;
      fillNodeInputs(r, gIn, etaBuf, cFixBuf);
      if (!mat.ideal) materialAt(s, mat, reg.background, etaBuf, regMode[r], cFixBuf, 0, zetaScratch, ccScratch);
      const zt = mat.ideal ? Math.log(cFixBuf[i] / mat.cRef[i]) : zetaScratch[i];
      return (u[s.blockOfNode[gOut] * M + 1 + i] - mat.mu0[i] / model.RT - zt) / z[i];
    }
    const inside = r > nr ? 1 : 0;
    const bOut = s.blockOfNode[gOut] * M;
    const phRx = balancedAt(r, grid.x[gIn], gIn, itf.reactions, itf.surface ?? [], inside, (p) => u[bOut + 1 + p.i]);
    if (Number.isFinite(phRx)) return phRx;
    i = z.findIndex((_, j) => shared(j) && passesRx(j));
    if (i < 0) i = z.findIndex((_, j) => shared(j));
    if (i < 0) return NaN;
    const mat = isGraded(reg) ? gradedAt(base, reg, grid.x[gIn]) : base;
    fillNodeInputs(r, gIn, etaBuf, cFixBuf);
    if (!mat.ideal) materialAt(s, mat, reg.background, etaBuf, regMode[r], cFixBuf, 0, zetaScratch, ccScratch);
    const zt = mat.ideal ? Math.log(cFixBuf[i] / mat.cRef[i]) : zetaScratch[i];
    return (u[bOut + 1 + i] - mat.mu0[i] / model.RT - zt) / z[i];
  };

  // Tier 2 (held electrodes on either face, held port electrodes, held adjacent insertion hosts).
  // A blocking capacitive/pinned metal face (Tier 5) anchors a dielectric directly (e.g. Al | SiO2
  // in a MOS), and otherwise serves only as a fallback for chains with no Tier 1–4 anchor, so it
  // never introduces a fake ohmic drop across a conducting electrolyte anchored by a Faradaic
  // electrode (e.g. echem 36). Tier 6: fallback floating electrodes on unanchored chains.
  const metalCapPhi = (itf, gOut, ci) => {
    if (itf.phi.type !== 'capacitive' && itf.phi.type !== 'pinned') return NaN;
    const Vm = (s.VT * u[s.blockOfNode[gOut] * M + 1 + ci]) / z[ci];
    return (Vm - itf.zeroCharge) / s.VT;
  };
  const isDielectric = (r) => {
    const g0 = grid.regionStart[r], g1 = grid.regionEnd[r];
    if (g1 <= g0) return false;
    for (let g = g0; g < g1; g++) if (segSigma[g] > 0) return false;
    return true;
  };
  const applyBlockingMetalFaces = (r, allowConducting) => {
    const reg = regions[r], base = materials[reg.material];
    if (base.phiFree || responds[r] || capped[r] || bulkFed[r]) return false;
    if (!allowConducting && !isDielectric(r)) return false;
    let added = false;
    if (r === 0 && isBlockingCapPin('left') && Number.isFinite(endPhiL)) {
      setAnchor(0, 5, endPhiL);
      added = true;
    }
    if (r === R - 1 && isBlockingCapPin('right') && Number.isFinite(endPhiR)) {
      setAnchor(N - 1, 5, endPhiR);
      added = true;
    }
    if (r > 0) {
      const lm = materials[regions[r - 1].material];
      if (lm.conductor && metalResolved[r - 1]) {
        const cp = metalCapPhi(interfaces[r - 1], grid.regionEnd[r - 1], lm.conductor.i);
        if (Number.isFinite(cp)) { setAnchor(grid.regionStart[r], 5, cp); added = true; }
      }
    }
    if (r + 1 < R) {
      const rm = materials[regions[r + 1].material];
      if (rm.conductor && metalResolved[r + 1]) {
        const cp = metalCapPhi(interfaces[r], grid.regionStart[r + 1], rm.conductor.i);
        if (Number.isFinite(cp)) { setAnchor(grid.regionEnd[r], 5, cp); added = true; }
      }
    }
    return added;
  };
  const applyElectrodes = (r, allowFloating, tier) => {
    const reg = regions[r], base = materials[reg.material];
    if (base.phiFree || responds[r] || capped[r] || bulkFed[r]) return false;
    let added = false;
    // Left metal or held insertion host
    if (r > 0) {
      const leftMat = materials[regions[r - 1].material];
      const gIn = grid.regionStart[r], gOut = grid.regionEnd[r - 1];
      if (leftMat.conductor && !interiorMetal[r - 1] && (heldMetal[r - 1] || allowFloating)) {
        const bOut = s.blockOfNode[gOut] * M;
        const ph = balancedAt(r, grid.x[gIn], gIn, interfaces[r - 1].reactions, interfaces[r - 1].surface ?? [], 1, (p) => u[bOut + 1 + p.i]);
        if (Number.isFinite(ph)) { setAnchor(gIn, tier, ph); added = true; }
      } else if (leftMat.phiFree && !leftMat.conductor && hostHeld[r - 1] && !allowFloating) {
        const ph = anchorFromHost(r, r - 1, gIn, gOut);
        if (Number.isFinite(ph)) { setAnchor(gIn, tier, ph); added = true; }
      }
    }
    // Right metal or held insertion host
    if (r + 1 < R) {
      const rightMat = materials[regions[r + 1].material];
      const gIn = grid.regionEnd[r], gOut = grid.regionStart[r + 1];
      if (rightMat.conductor && !interiorMetal[r + 1] && (heldMetal[r + 1] || allowFloating)) {
        const bOut = s.blockOfNode[gOut] * M;
        const ph = balancedAt(r, grid.x[gIn], gIn, interfaces[r].reactions, interfaces[r].surface ?? [], 0, (p) => u[bOut + 1 + p.i]);
        if (Number.isFinite(ph)) { setAnchor(gIn, tier, ph); added = true; }
      } else if (rightMat.phiFree && !rightMat.conductor && hostHeld[r + 1] && !allowFloating) {
        const ph = anchorFromHost(r, r + 1, gIn, gOut);
        if (Number.isFinite(ph)) { setAnchor(gIn, tier, ph); added = true; }
      }
    }
    // Port electrodes in this region (held in Tier 4, floating only in fallback Tier 6)
    for (let p = 0; p < ports.length; p++) {
      const port = ports[p];
      if (port.region !== r || port.reactions.length === 0) continue;
      if (isHeldTerm(2 + p) === allowFloating) continue;
      const gMid = port.nodes[port.nodes.length >> 1];
      const mid = 0.5 * (grid.x[port.nodes[0]] + grid.x[port.nodes[port.nodes.length - 1]]);
      const ph = balancedAt(r, mid, gMid, port.reactions, port.surface, 0, (q) => s.portEta(port, q.i));
      if (Number.isFinite(ph)) {
        for (const g of port.nodes) setAnchor(g, tier, ph);
        added = true;
      }
    }
    return added;
  };
  // --- Pass 3: outward propagation, ohmic interpolation, and unheld/interior conductors ---
  const regionResolved = new Uint8Array(R);
  const metalResolved = new Uint8Array(R);
  for (let r = 0; r < R; r++) if (materials[regions[r].material].conductor && heldMetal[r]) metalResolved[r] = 1;
  for (let r = 0; r < R; r++) {
    applyElectrodes(r, false, 4);
    applyBlockingMetalFaces(r, false);
  }

  // Contiguous chains of non-φ-free regions.
  const chains = [];
  for (let r = 0; r < R; ) {
    if (materials[regions[r].material].phiFree) { r++; continue; }
    const r0 = r;
    while (r + 1 < R && !materials[regions[r + 1].material].phiFree) r++;
    chains.push([r0, r, grid.regionStart[r0], grid.regionEnd[r]]);
    r++;
  }

  const phiTilde = new Float64Array(N);
  const writeChain = (r0, r1, g0, g1, useZero) => {
    if (useZero) {
      for (let g = g0; g <= g1; g++) phiTilde[g] = 0;
    } else {
      const anchors = [];
      for (let g = g0; g <= g1; g++) if (ancTier[g] > 0) anchors.push(g);
      if (anchors.length === 0) return false;
      let hasInsulator = false;
      for (let g = g0; g < g1; g++) {
        if (grid.nodeRegion[g] === grid.nodeRegion[g + 1] && !(segSigma[g] > 0) && segEps[g] > 0) { hasInsulator = true; break; }
      }
      if (!hasInsulator && anchors.length === 2 && anchors[0] === g0 && anchors[1] === g1) {
        // An all-c0 conducting chain between two boundary anchors (e.g. two electrodes or two
        // baths): keep φ̃ flat at their symmetric midpoint so amount-fixed spectators don't start
        // with an artificial bulk drift field into blocked boundaries.
        const pMid = 0.5 * (ancPhiTilde[g0] + ancPhiTilde[g1]);
        for (let g = g0; g <= g1; g++) phiTilde[g] = pMid;
      } else {
        const first = anchors[0], last = anchors[anchors.length - 1];
        for (let g = g0; g < first; g++) phiTilde[g] = ancPhiTilde[first];
        for (let g = last + 1; g <= g1; g++) phiTilde[g] = ancPhiTilde[last];
        for (let k = 0; k < anchors.length; k++) phiTilde[anchors[k]] = ancPhiTilde[anchors[k]];
        for (let k = 0; k + 1 < anchors.length; k++) {
          const a = anchors[k], b = anchors[k + 1];
          if (b - a <= 1) continue;
          const f = interpFractions(a, b);
          const pa = ancPhiTilde[a], pb = ancPhiTilde[b];
          for (let g = a + 1; g < b; g++) phiTilde[g] = (1 - f[g - a]) * pa + f[g - a] * pb;
        }
      }
    }
    for (let r = r0; r <= r1; r++) {
      const reg = regions[r], base = materials[reg.material], graded = isGraded(reg);
      const profiled = reg.c0Profile.some((p) => p !== null);
      const ga = grid.regionStart[r], gb = grid.regionEnd[r], dip = cumDip[r];
      let uniform = !graded && !profiled && !regRamped[r];
      if (uniform) {
        const p0 = phiTilde[ga];
        for (let g = ga + 1; g <= gb; g++) if (phiTilde[g] !== p0) { uniform = false; break; }
      }
      if (uniform) {
        const ph = phiTilde[ga] + dip;
        fillNodeInputs(r, ga, etaBuf, cFixBuf);
        if (!base.ideal) materialAt(s, base, reg.background, etaBuf, regMode[r], cFixBuf, ph, zetaScratch, ccScratch);
        for (let g = ga; g <= gb; g++) {
          const b = s.blockOfNode[g] * M;
          u[b] = ph;
          for (let i = 0; i < n; i++) {
            if (regMode[r][i] === 1) u[b + 1 + i] = etaBuf[i];
            else if (regMode[r][i] === 2) {
              const zt = base.ideal ? Math.log(cFixBuf[i] / base.cRef[i]) : zetaScratch[i];
              u[b + 1 + i] = zt + base.mu0[i] / model.RT + z[i] * ph;
            }
          }
        }
      } else {
        for (let g = ga; g <= gb; g++) {
          const mat = graded ? gradedAt(base, reg, grid.x[g]) : base;
          fillNodeInputs(r, g, etaBuf, cFixBuf);
          writeNode(r, g, phiTilde[g] + dip, mat, regMode[r], etaBuf, cFixBuf);
        }
      }
      regionResolved[r] = 1;
    }
    return true;
  };

  // Reaction equilibrium level for an interior conductor carrier on one of its faces.
  const metalFaceEta = (itf, metalSide, gAdj, carrier) => {
    const rx = itf.reactions.find((x) => x.part.some((p) => p.side === metalSide && p.i === carrier));
    if (!rx) return NaN;
    const bn = s.blockOfNode[gAdj] * M;
    let a = rx.fixedA, nu = 0;
    const surf = itf.surface ?? [], bare = 1 - surf.reduce((t, sp) => t + sp.theta0, 0);
    for (const q of rx.surf ?? []) a -= q.nu * (surf[q.s].mu0 / model.RT + Math.log(surf[q.s].theta0 / bare));
    for (const p of rx.part) {
      if (p.side === metalSide) nu += p.nu;
      else a -= p.nu * u[bn + 1 + p.i];
    }
    return nu !== 0 && Number.isFinite(a) ? a / nu : NaN;
  };

  for (let pass = 0; pass <= R + 1; pass++) {
    let progress = false;
    for (const [r0, r1, g0, g1] of chains) {
      if (writeChain(r0, r1, g0, g1, false)) {
        // Check if any unresolved interior metal beside this chain can now be set.
      }
    }
    // Resolve unheld conductors (interior metals AND metals attached only to a floating contact)
    // from whichever adjacent face(s) are resolved and carry a reaction for their carrier.
    for (let r = 0; r < R; r++) {
      const base = materials[regions[r].material];
      if (!base.conductor || metalResolved[r]) continue;
      const ci = base.conductor.i;
      const etaL = r > 0 && regionResolved[r - 1] ? metalFaceEta(interfaces[r - 1], 1, grid.regionEnd[r - 1], ci) : NaN;
      const etaR = r + 1 < R && regionResolved[r + 1] ? metalFaceEta(interfaces[r], 0, grid.regionStart[r + 1], ci) : NaN;
      let etaVal = NaN;
      if (Number.isFinite(etaL) && Number.isFinite(etaR)) etaVal = 0.5 * (etaL + etaR);
      else if (Number.isFinite(etaL)) etaVal = etaL;
      else if (Number.isFinite(etaR)) etaVal = etaR;
      if (!Number.isFinite(etaVal)) continue;
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const b = s.blockOfNode[g] * M;
        u[b] = 0;
        u[b + 1 + ci] = etaVal;
      }
      metalResolved[r] = 1;
      progress = true;
      // Can this newly resolved metal anchor an unanchored non-φ-free neighbour (by reaction, or
      // across a blocking capacitive/pinned face)?
      for (const [nr, itf, inside, gIn, gOut] of [
        [r - 1, r > 0 ? interfaces[r - 1] : null, 0, r > 0 ? grid.regionEnd[r - 1] : -1, grid.regionStart[r]],
        [r + 1, r + 1 < R ? interfaces[r] : null, 1, r + 1 < R ? grid.regionStart[r + 1] : -1, grid.regionEnd[r]],
      ]) {
        if (nr < 0 || nr >= R || regionResolved[nr] || materials[regions[nr].material].phiFree || responds[nr] || capped[nr] || bulkFed[nr]) continue;
        const bOut = s.blockOfNode[gOut] * M;
        const ph = balancedAt(nr, grid.x[gIn], gIn, itf.reactions, itf.surface ?? [], inside, (p) => u[bOut + 1 + p.i]);
        if (Number.isFinite(ph)) { setAnchor(gIn, 4, ph); progress = true; }
        else {
          const cp = metalCapPhi(itf, gOut, ci);
          if (Number.isFinite(cp)) { setAnchor(gIn, 5, cp); progress = true; }
        }
      }
    }
    if (progress) continue;
    // Fallback Tier 5: blocking capacitive/pinned metal faces on still-unresolved chains.
    let fellBack = false;
    for (const [r0, r1] of chains) {
      if (regionResolved[r0]) continue;
      for (let r = r0; r <= r1; r++) if (applyBlockingMetalFaces(r, true)) fellBack = true;
    }
    if (fellBack) continue;
    // Fallback Tier 6: floating electrodes on still-unresolved chains.
    for (const [r0, r1] of chains) {
      if (regionResolved[r0]) continue;
      for (let r = r0; r <= r1; r++) if (applyElectrodes(r, true, 6)) fellBack = true;
    }
    if (fellBack) continue;
    // Any chain still without an anchor starts at φ̃ = 0 (face dipoles still honoured).
    for (const [r0, r1, g0, g1] of chains) {
      if (!regionResolved[r0]) writeChain(r0, r1, g0, g1, true);
    }
    // Any interior conductor still without a balancing reaction takes the mean adjacent φ̂ (else 0).
    for (let r = 0; r < R; r++) {
      if (!interiorMetal[r] || metalResolved[r]) continue;
      const ci = materials[regions[r].material].conductor.i;
      let sum = 0, cnt = 0;
      if (r > 0 && !materials[regions[r - 1].material].phiFree) { sum += u[s.blockOfNode[grid.regionEnd[r - 1]] * M]; cnt++; }
      if (r + 1 < R && !materials[regions[r + 1].material].phiFree) { sum += u[s.blockOfNode[grid.regionStart[r + 1]] * M]; cnt++; }
      const etaVal = z[ci] * (cnt > 0 ? sum / cnt : 0);
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const b = s.blockOfNode[g] * M;
        u[b] = 0;
        u[b + 1 + ci] = etaVal;
      }
      metalResolved[r] = 1;
    }
    break;
  }

  // --- Pass 4: unheld insertion hosts (gauge-shifted to neighbour's passing ion) + unheld collectors ---
  for (let r = 0; r < R; r++) {
    const reg = regions[r], base = materials[reg.material];
    if (!base.phiFree || base.conductor || hostHeld[r]) continue;
    const graded = isGraded(reg);
    for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
      const mat = graded ? gradedAt(base, reg, grid.x[g]) : base;
      fillNodeInputs(r, g, etaBuf, cFixBuf);
      writeNode(r, g, 0, mat, regMode[r], etaBuf, cFixBuf);
    }
  }
  for (let r = 0; r < R; r++) {
    const reg = regions[r], mat = materials[reg.material];
    if (!mat.phiFree || mat.conductor || hostHeld[r]) continue;
    const sides = [
      [r - 1, grid.regionStart[r], r > 0 ? grid.regionEnd[r - 1] : -1],
      [r + 1, grid.regionEnd[r], r + 1 < R ? grid.regionStart[r + 1] : -1],
    ];
    let dSum = 0, dCnt = 0;
    for (const [nr, gIn, gOut] of sides) {
      if (gOut < 0) continue;
      const other = materials[regions[nr].material];
      if (other.phiFree || other.conductor) continue;
      const itf = interfaces[Math.min(r, nr)], shared = (j) => z[j] !== 0 && mat.present[j] && other.present[j];
      const passes = (j) => itf.links[j].type !== 'blocked' || itf.reactions.some((rx) => rx.part.some((p) => p.i === j && p.side === 0) && rx.part.some((p) => p.i === j && p.side === 1));
      let i = z.findIndex((_, j) => shared(j) && passes(j));
      if (i < 0) i = z.findIndex((_, j) => shared(j));
      if (i < 0) continue;
      const d = (u[s.blockOfNode[gOut] * M + 1 + i] - u[s.blockOfNode[gIn] * M + 1 + i]) / z[i];
      if (Number.isFinite(d)) { dSum += d; dCnt++; }
    }
    if (dCnt > 0) {
      const d = dSum / dCnt;
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const b = s.blockOfNode[g] * M;
        u[b] += d;
        for (let j = 0; j < n; j++) if (mat.present[j]) u[b + 1 + j] += z[j] * d;
      }
      // Propagate the host's shifted carrier level into an adjacent unheld metal collector sharing
      // that carrier (e.g. cathode | Cu | floating I = 0 starts at OCV, not 0 V).
      for (const [nr, gHost] of [[r - 1, grid.regionStart[r]], [r + 1, grid.regionEnd[r]]]) {
        if (nr < 0 || nr >= R) continue;
        const other = materials[regions[nr].material];
        if (!other.conductor || heldMetal[nr]) continue;
        const ci = other.conductor.i;
        if (!mat.present[ci] || interfaces[Math.min(r, nr)].links[ci].type === 'blocked') continue;
        const etaCar = u[s.blockOfNode[gHost] * M + 1 + ci];
        for (let g = grid.regionStart[nr]; g <= grid.regionEnd[nr]; g++) u[s.blockOfNode[g] * M + 1 + ci] = etaCar;
      }
    }
  }

  // --- Step 4: floating contacts ---
  const level = (side) => {
    const ct = model.contacts[side], i = ct.terminal, g = side === 'left' ? 0 : N - 1, b = s.blockOfNode[g];
    if (ct.she && !ct.she.reference && !s.phiUndefined[g]) return s.VT * u[b * M];
    if (i === null || ct.species[i].type !== 'equilibrium' || z[i] === 0) return NaN;
    return (s.VT * u[b * M + 1 + i]) / z[i] - ct.species[i].offset;
  };
  for (const k of s.floating) {
    const t = s.terms[k], other = k === 0 ? 1 : 0;
    if (t.kind === 'port') continue;
    const lv = level(t.side);
    if (!Number.isFinite(lv)) continue;
    if (other < s.terms.length && isHeldTerm(other)) {
      const lvHeld = level(s.terms[other].side);
      if (Number.isFinite(lvHeld)) { s.termV[k] = s.termV[other] + (lv - lvHeld); continue; }
    }
    // Anchored by a held port (no held opposite contact): lv is already in the held gauge.
    if (s.terms.some((_, h) => h >= 2 && isHeldTerm(h))) s.termV[k] = lv;
  }
  for (let g = 0; g < s.nNodes; g++) if (s.surfAt[g] >= 0) s._surfaceStart(g);
  s.computeConcentrations();

  // --- Step 6: floating double layers start uncharged ---
  for (const k of s.floating) {
    const t = s.terms[k], port = t.kind === 'port' ? ports[t.index] : null;
    const gate = t.kind === 'contact' ? model.contacts[t.side].phi : null;
    if (gate?.type === 'capacitive') {
      const g = t.side === 'left' ? 0 : s.nNodes - 1;
      if (!s.phiUndefined[g]) s.termV[k] = gate.zeroCharge + s.VT * u[s.blockOfNode[g] * M];
      continue;
    }
    if (!port?.capacitance) continue;
    let w = 0, sum = 0;
    port.nodes.forEach((g, j) => {
      if (!s._portHolds(g)) return;
      const a = grid.vol[g] * s.portArea[t.index][j];
      w += a;
      sum += a * s.VT * u[s.blockOfNode[g] * M];
    });
    if (w > 0) s.termV[k] = port.capacitance.zeroCharge + sum / w;
  }

  // --- Step 7: floating port linked to its carrier (no reactions, no capacitance) ---
  for (const k of s.floating) {
    const t = s.terms[k], port = t.kind === 'port' ? ports[t.index] : null;
    const i = port?.terminal, link = i === null || i === undefined ? null : port.species[i];
    if (!port || port.reactions.length > 0 || port.capacitance || !link || link.type === 'blocked' || z[i] === 0) continue;
    let w = 0, sum = 0;
    port.nodes.forEach((g, j) => {
      const a = grid.vol[g] * (s.portArea[t.index]?.[j] ?? port.fill[j]);
      w += a;
      sum += (a * s.VT * u[s.blockOfNode[g] * M + 1 + i]) / z[i];
    });
    if (!(w > 0)) continue;
    let refSum = 0, refCnt = 0;
    for (const h of [0, 1]) {
      if (h < s.terms.length && isHeldTerm(h)) {
        const lv = level(s.terms[h].side);
        if (Number.isFinite(lv)) { refSum += s.termV[h] - lv; refCnt++; }
      }
    }
    if (refCnt > 0) s.termV[k] = sum / w - (link.offset ?? 0) + refSum / refCnt;
    else if (s.terms.some((_, h) => h >= 2 && isHeldTerm(h))) s.termV[k] = sum / w - (link.offset ?? 0);
  }

  // --- Step 8: gates at α/(α + β) at the starting voltage ---
  for (const gt of s.gateList) {
    const { gate, o } = gt, V = s._gateVoltage(gt), a = gateRate(gate.alpha, V)[0], b = gateRate(gate.beta, V)[0];
    u[o] = a + b > 0 ? a / (a + b) : 0;
  }

  // --- Step 9: floating reactive port without capacitance, by bisection ---
  for (const k of s.floating) {
    const t = s.terms[k], port = t.kind === 'port' ? ports[t.index] : null;
    if (!port || port.reactions.length === 0 || port.capacitance) continue;
    const src = sourceAt(t.drive.src, s.time), want = (V) => (t.drive.kind === 'I' ? src : (src - V) / t.drive.R);
    const current = (V) => {
      s.termV[k] = V;
      let I = 0;
      port.reactions.forEach((rx) => {
        let q = 0;
        for (const p of rx.part) if (p.side === 0) q += p.nu * z[p.i];
        port.nodes.forEach((g, w) => (I += q * FARADAY * grid.vol[g] * s.portArea[t.index][w] * s._portRate(port, rx, g)));
      });
      for (let i = 0; i < n; i++) {
        const link = port.species[i];
        if (z[i] === 0 || (link.type !== 'conductance' && link.type !== 'exchange')) continue;
        const kk = link.type === 'conductance' ? (link.G * s.VT) / (z[i] * z[i] * FARADAY) : link.k, gated = link.gates ?? [], q0 = 1 + n + port.surface.length;
        port.nodes.forEach((g, w) => {
          const b = s.blockOfNode[g], open = gated.length ? s._open(gated, (q) => u[b * M + q0 + q]) : 1;
          I += z[i] * FARADAY * grid.vol[g] * port.fill[w] * kk * open * (s.portEta(port, i) - u[b * M + 1 + i]);
        });
      }
      return I;
    };
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
 * Steady-solve diffuse-layer warm-up: balances each one-sided diffuse layer's 1D Stern–Gouy–Chapman
 * relation q_d(s) = Q(s) at the cold start's bulk concentrations and electrode potentials, setting
 * its plane drop s and face displacement D before the first direct steady Newton step. Called from
 * src/steady.js AFTER captureLaws / referenceAmounts so conserved amounts and transient cold starts
 * keep uncharged double layers (s = 0, D = 0, Γ = 0).
 */
export function initSteadyDiffuse(s) {
  if (s.diffList.length === 0) return;
  const { model, n, M, u, z, VT } = s;
  const { interfaces } = model;
  let moved = false;
  for (const dl of s.diffList) {
    const itf = interfaces[dl.f], law = itf.phi.type;
    if (law !== 'capacitive' && law !== 'pinned') continue;
    const [dlL, dlR] = s.diffAt[dl.f];
    if (dlL >= 0 && dlR >= 0) continue; // two-sided diffuse face: keep s = 0
    const bf = s.blockOfFace[dl.f], bL = bf - 1, bR = bf + 1;
    let gap0, sg;
    if (itf.conductor) {
      const condSide = itf.conductor.side, im = itf.conductor.i;
      const bm = condSide === 'left' ? bL : bR, bo = condSide === 'left' ? bR : bL;
      sg = condSide === 'left' ? 1 : -1;
      gap0 = u[bm * M + 1 + im] / z[im] - itf.zeroCharge / VT - u[bo * M];
    } else {
      sg = dl.side ? -1 : 1;
      const dip = itf.dipole / VT;
      gap0 = dl.side ? u[bL * M] - u[bR * M] + dip : u[bR * M] - u[bL * M] - dip;
    }
    if (!Number.isFinite(gap0)) continue;
    const cSlice = s.c.subarray(dl.g * n, dl.g * n + n);
    if (law === 'pinned') {
      u[dl.o] = gap0;
      diffuseLayer(cSlice, z, n, gap0, dl.units, dl.dl, dl.lat, null);
      u[bf * M] = dl.side ? -(dl.dl.q + itf.sheetCharge) : dl.dl.q;
      moved = true;
      continue;
    }
    const kC = itf.phi.C * VT, sheet = dl.side ? itf.sheetCharge : 0;
    const fVal = (sh) => {
      diffuseLayer(cSlice, z, n, sh, dl.units, dl.dl, dl.lat, null);
      return dl.dl.q + kC * (gap0 - sh) + sheet;
    };
    let lo = Math.min(0, gap0) - 1, hi = Math.max(0, gap0) + 1;
    while (fVal(lo) < 0 && lo > -1e4) lo -= 2 * (hi - lo);
    while (fVal(hi) > 0 && hi < 1e4) hi += 2 * (hi - lo);
    if (!(fVal(lo) >= 0 && fVal(hi) <= 0)) continue;
    for (let it = 0; it < 60 && hi - lo > 1e-12 * Math.max(1, Math.abs(gap0)); it++) {
      const m = 0.5 * (lo + hi);
      if (fVal(m) > 0) lo = m;
      else hi = m;
    }
    const sh = 0.5 * (lo + hi);
    u[dl.o] = sh;
    u[bf * M] = sg * kC * (gap0 - sh);
    moved = true;
  }
  if (moved) s.computeConcentrations();
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

// Every strictly neutral node onto its Gauss row (src/neutral.js) where the start left one off
// beside a capacitor, at t = 0⁺ (solveNeutral). c0 is the initial condition and the capacitors
// start uncharged: a box a port's or a contact's capacitance holds charge in keeps its c0 and its
// φ̂ jumps until it's neutral, a floating capacitor keeping the charge it started with by its
// voltage. A box beside a capacitive face holds what a sheet charge there puts on it at t = 0⁺ (a
// double layer charges later): its φ̂ at its levels, each closed stretch's amount held at what the
// start gave it. (Left off, the first step charged a sheet's edge boxes, 1.6e-2 C/m², through the
// box beside them. φ̂ held where every ion was closed, a sealed cell between plates at 0 and 1 V
// couldn't be neutral, the start was put back as it was, and its first 1 ns step's currents summed
// to −2e8 A/m², the right plate settling with the wrong sign, unwarned; a floating port over a held
// one without c0 wasn't seen at all; Neon's fourth review, R4 1.6.) An interior box whose own c0
// can't be neutral is the start's to refuse, and says so at the first step.
// Whether node g's box holds a capacitive face's or contact's displacement.
function boundaryCapacitive(s, g) {
  const grid = s.model.grid, r = grid.nodeRegion[g], itfs = s.model.interfaces;
  if (g === grid.regionEnd[r] && r < itfs.length && itfs[r].phi.type === 'capacitive') return true;
  if (g === grid.regionStart[r] && r > 0 && itfs[r - 1].phi.type === 'capacitive') return true;
  const side = g === 0 ? 'left' : g === s.nNodes - 1 ? 'right' : null;
  return !!side && s.model.contacts[side].phi.type === 'capacitive';
}

function neutralStart(s) {
  const all = strictlyNeutralNodes(s);
  if (all.length === 0) return;
  const rows = gaussRows(s, all);
  // (a port's or a contact's capacitance in the box: it jumps; beside a face: by its levels)
  const end = (g) => (g === 0 ? 'left' : g === s.nNodes - 1 ? 'right' : null);
  const held = (g) => (s.qgSlope[g] ?? 0) > 0 || (!!end(g) && s.model.contacts[end(g)].phi.type === 'capacitive');
  const nodes = all.filter((nd) => !nd.bulk && (held(nd.g) || boundaryCapacitive(s, nd.g)));
  // (only where one is off past the start's own round-off: its solves leave 1e-12 of the gross)
  if (!nodes.some((nd) => Math.abs(rows[all.indexOf(nd)].G) > 1e-10 * rows[all.indexOf(nd)].scale)) return;
  for (const nd of nodes) nd.mode = held(nd.g) ? 'jump' : 'level';
  const closed = s.stretches.flatMap((st) => (st.connected || st.reactive ? [] : [{ st, want: s.amount(st) }])).filter((t) => t.want > 0);
  const neutralAt = new Uint8Array(s.nNodes);
  for (const nd of nodes) neutralAt[nd.g] = 1;
  s._assembleBookkeeping(Infinity);
  const floating = [];
  for (const k of s.floating) {
    if (k >= 2) {
      const p = s.model.ports[k - 2];
      if (p.capacitance && p.nodes.some((g) => neutralAt[g] && s._portHolds(g))) floating.push({ k, want: s.portQ[k - 2] });
    } else if (neutralAt[k === 0 ? 0 : s.nNodes - 1] && s.model.contacts[k === 0 ? 'left' : 'right'].phi.type === 'capacitive') floating.push({ k, want: k === 0 ? s.contactD.left : -s.contactD.right });
  }
  // The sheet's counter-charge belongs at the face. Shifting a closed level outside these boxes
  // changed A⁻ alone throughout the neutral interior by 0.017 mol/m³ (the membrane combos'
  // case 34, neutralG; Astra, 25 Nov).
  solveNeutral(s, nodes, { closed, floating, closedAtNodes: true });
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
