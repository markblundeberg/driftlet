# Roadmap

Where driftlet is, and where it's going. Scope stays sharp: 1D, local physics only (everything
must keep the Jacobian block-tridiagonal), pure JavaScript with no dependencies, fast enough
for live demos.

## Done

- **1.0** (2026-10-13), agent-first: the API audited and frozen (every option name checked, with
  the one meant suggested), llms.txt complete and within one read, warnings that say what to do,
  `check()` covering what an agent can't judge, and two rounds of cold agents building demos that
  agree with theory. From here, work follows users' reports.
- Solver core: graded grid with doubled interface nodes, block Thomas, Scharfetter–Gummel in
  quasi-Fermi (expm1) form, compensated double-double $`\bar\mu`$, damped Newton.
- Equilibrium and steady transport: direct steady solves when nothing is conserved, and
  exactly conserving giant steps when something is.
- Materials, regions and explicit per-interface alignment; interface $`\phi`$ laws (dipole, neutral,
  Helmholtz); strictly neutral ($`\varepsilon = 0`$) materials.
- Contacts: fixed levels with offsets, ohmic contacts, baths, gates and Stern layers,
  conductance links.
- Interfaces: blocking, conductance, and Butler–Volmer reactions with signed participants on
  either side (each rate an unknown of the face block); electrodes as conductor regions.
- Circuits: voltage, galvanostatic (incl. open circuit), load resistor.
- Bulk mass-action reactions.
- Strictly neutral regions solved in better-conditioned unknowns on transient steps (current
  continuity for one balance, $`\hat\phi'`$ and $`\eta - z\hat\phi`$ as unknowns): a 3 M KCl salt bridge
  with 1 µM traces of every other ion, between eddy-mixed solutions, runs from its unmixed start
  through days.
- Conserved combinations of reacting stretches (moieties) as constraint rows in the direct
  steady solve, as spectators already were: a closed cell with complexation solves straight to
  its steady state, cold or warm.
- Backward-Euler transients with step halving, and conservation bookkeeping.
- Statistics: one general per-material interface, $`c(\zeta)`$ with its capacitance matrix, and the
  excess-potential Scharfetter–Gummel flux. Fermi–Dirac, lattice gas (Bikerman crowding, site
  filling), Redlich–Kister, Debye–Hückel, custom functions.
- Neutral-combination (insertion) statistics for $`\varepsilon = 0`$ hosts, from an isotherm or a tabulated
  OCV curve, with $`\phi`$ undefined inside.
- Time and frequency: variable-step BDF2; adaptive `advance(tEnd, { tol, budgetMs })` for
  animation; small-signal impedance with complex block Thomas and complex profiles.
- Source continuation for cold starts at bias.
- A first performance pass: a benchmark suite (`npm run bench`, with deterministic counters checked
  against a committed baseline in CI), assembly without per-node
  allocation (about 2× faster), Newton started from the extrapolated state in adaptive
  transients (a third fewer iterations), and a tighter-damping retry for large jumps.
- Conductor regions (metals, fast ion conductors): the carrier's level and a conductivity only (mixed-form Ohm's law, surface charge as
  sheets), capacitive faces aligned by work function, electrode reactions at internal faces,
  bipolar electrodes.
- Internal ports: reservoirs attached to a window of interior nodes (held levels, volumetric
  conductance or exchange), e.g. grounding a 1D MOS channel, or a wire to a whole metal. And
  electrodes spread through a window (a port with Butler–Volmer reactions against its carrier, per
  area of electrode): the metal under a thin film, a crevice's walls, a porous electrode's matrix,
  held or floating at its mixed potential.
- A second-law validation: in steady state the free energy brought in through every terminal,
  $`\sum N\bar\mu`$ over every species, equals the dissipation inside (fluxes down their own $`\bar\mu`$
  drops, reactions by their affinities, interface and conductor resistance), each term $`\ge 0`$,
  exactly for the discrete equations.
- Terminals: the contacts and every port, each held at a voltage, behind a resistance or driven
  by a current (reference electrodes, three-electrode cells, grounded floating electrodes), with
  piecewise-linear waveforms (cyclic voltammetry), impedance at any terminal, and a fast `set()`
  for sources alone. Floating terminal voltages are bordered unknowns.
- Variable block sizes: each block holds only the unknowns that exist at its node, assembled
  straight into compact blocks by region kernels (conductor, dilute, concentrated); a face's
  block holds its displacement, its linked fluxes and its reaction rates. Bookkeeping after each
  step is read from the end boxes and port windows alone. Conductor regions are a single cell.
- Spectators solved directly in steady state, with their amounts as constraints (a bordered
  solve) instead of giant time steps: 5–8× fewer iterations where they occur.
- Transport extras: imposed advection (exact exponential fitting) and current-free eddy mixing
  (an Onsager term projected to carry no current).
- Validation suite (see the README), and docs.
- Shipped as 0.1.0: on npm and jsDelivr, with JSDoc types and a `.d.ts`, live demos on GitHub
  Pages (screenshots in the README), and CI on Node 22, 24 and 26.
- Releases by trusted publishing: a version tag stages the package on npm from GitHub Actions,
  with signed provenance and no stored token, and the owner approves it with 2FA (from 0.1.1).
- A real-browser smoke test of the demos in CI (`npm run smoke`): headless Chrome over the
  DevTools protocol, every control moved to its extremes, failing on any exception, console
  error, failed request or NaN drawn.
- Position-dependent rate constants (`kf` as a profile against x, averaged over each node's box
  so the total is exact on any grid), and Beer–Lambert `photogeneration()` in the kit, validated
  against collection theory from blue to red.
- Membranes as faces: a capacitor with Goldman–Hodgkin–Katz permeabilities, and pumps as
  saturating reactions on it (validated against GHK, Mullins–Noda, the static head and Donnan),
  from the membrane agent's debrief.
- Checks as results: `check(device, sol)` reports each species' ledger (J = F(G − R) under
  light), conservation, and the same device on a grid twice as fine, asked for by the perovskite
  agent.
- Curated devices with reference results from an independent code (`perovskiteCell()` and
  IonMonger's scans of it), and a transient recorder (`recorder()`: the trace across calls,
  frames at set times for scrubbing), both asked for by the perovskite agent.
- Screenshots retaken by script (`npm run shots`): each demo held at a set state (`?t=…`,
  `?phase=…`) on virtual time, so an unchanged page comes out byte-identical.
- Particle cartoons consistent with the model (`particles()` in the kit): dots drawn as a
  sample of a solution's concentrations, hopping at its one-way rates, so their density and net
  crossings agree with the fields at any step.

## Next

0. **One owner for strictly neutral charge (the next design piece).** A strictly neutral node is
   an index-1 DAE: its balances with an algebraic G(u) = 0, ions plus fixed charge plus a port's
   spread charge plus what its box holds of a face's or contact's D. That G is assembled in
   src/regions.js and src/boundaries.js. But src/start.js (equilibrateFaces' qOther, neutralPhi)
   and src/carry.js (netOf, neutralize) each sum their own subset of it. And _chargeRows, assuming
   every step starts on G = 0, deletes any disagreement from the pivot ion. Neon's 1.1 and R2 3.1
   (and tonight's _layOut) were that seam, and so were the 37 tests an exact _chargeRows broke.
   Neon's reading (its round-2 closing note, 20 Nov 2026): they are the map of every path relying
   on the silent deletion. The plan:
   - one gaussResidual(g)/projectNeutral(nodes), from the same terms the Gauss row assembles;
   - every state producer (initFromComposition, equilibrateFaces, carryState, the set()s that move
     charge) projecting through it;
   - then _chargeRows exact (the old charge at every node), with a check that no node enters a step
     off G = 0, the 37 tests read one by one as what each was relying on.

   Done first (20 Nov), the last step:
   - The exact rows were the 37's 16 by then, the other doors closed one by one in the meantime.
     All 16 were round-off: a solved state's net charge of ~1e-16 of its ions', over a step of
     1e-18 s, a current swamping the fluxes. So the rows carry a start's net charge exactly,
     except below 1e-12 of the ions' (round-off) where no gate is.
   - Probed through the whole suite, every real path enters its steps neutral to Newton's
     tolerance (3e-10 at most). The only starts further off were a test that perturbs u on purpose
     and a start that is refused.
   - Past 1e-8, a solution says so.
   - Neon's R2 3.1 device, with its carry fix undone, now keeps both ions to 1e-16 and fails its
     step, saying where. Before, it lost Q_dl/2F unwarned.

   Still open:
   - The producers' single owner, so the starts are right by construction.
   - A start far off neutral can't be carried in one step after a solve: it needs a current of
     δq/dt, which shortening only raises. (Below: Neon finds no physical path that needs it but a
     step in fixed charge.)
   - So does a start at round-off, in proportion: the first step after a solve reads the start's
     Gauss residual, ~1 to 10 ε of the stored ionic charge, as δq/dt through the terminals and the
     ledger alike: 3.0e-8 A/m² in an 86 ns step (the membrane combos' case 18), and 1e-3 A/m² it
     would be at 1 ps. bench/combos' kirchhoff allows exactly that offset (bathStartRoundoff in
     bench/combos/charge.js; Sol's chase, 23 Nov). Removing it would mean a start exact in Gauss.
   - A solve past a statistics model's range takes as long to refuse as its fallback ladder takes
     to search: the Debye–Hückel fixture past its range (test/fixtures) spends 35 s, and 159 of its
     222 held solves converge, so the time is real search, not wasted retries. Two orders of the
     same nested search run (an R-backed contact outside a floating one, then inside it); the
     second has rescued nothing in the suite, combos --n 40 or stress, but nothing there needs it
     either, so it isn't proved idle. A wall detector on the ramp's failures was ruled out: a
     passable basin 1 µV wide gave 18 failures before the ramp reached its target (Astra's chase,
     23 Nov; Neon's analysis before it). What would help: a solvable fixture that needs the second
     order, then sharing what the first order learned with it.

   Neon's third review (21 Nov 2026) audited every sum of this charge (13 places) and designed the
   owner; its plan replaces the one above where they differ:
   - The owner is built from shared term functions (a contact's and a face's capacitive D, a port's
     σ at a node) that both assembly and a gaussResidual(g) call, not read from the assembled Gauss
     row: the producers evaluate it at trial φ̂ many times a node, and a face's D is its own unknown.
     A standing check holds the two equal on every combos family.
   - Each producer's projection is its own:
     - the cold start, φ̂ at fixed η where a bath sets the levels, at fixed c where c0 does;
     - the face start, each passing ion's amount over the pair kept and a floating gate following
       the φ̂ moved under it (done: R3 1.7);
     - a regrid, φ̂ held (so a face's D and a gate's σ are kept) and each η_i shifted by −z_i δ, the
       concentrations moving in the charge direction only;
     - a set() of a gate's C, or a gate added: the concentrations kept and φ̂ jumping at the gates'
       total charge, each port's charge before it kept as its start (a new one's, zero), so the
       first step reads the transfer between gates as currents that sum to zero (R3 1.5: until
       then, the step says its start was off);
     - a step in a gate's V: nothing; the start is on G = 0 at the old drives, and the step solves
       the RC.
   - The t = 0⁺ impulse by transference numbers is dropped. A gate's step has no impulse (φ jumps,
     c and the gate's charge are continuous, the current finite over the RC time); and a local
     t_i split makes ions in a sealed cell, drops i·∇t_i, and can turn a co-ion negative. Only a
     step in a strictly neutral region's fixed charge needs one: refuse it, or solve the
     migration-only projection, −∇·(κ∇Φ) = Δρ_fixed with Φ the potential's impulse.
   - A bench/combos property, neutralG: after each producer (cold start, regrid, a drive step, a
     port's change), every strictly neutral node on G = 0 to 1e-12 of its scale; then one short
     step, no offNeutral warning.

   Done in round 3: the face start's floating gates (R3 1.7), and a step's check of each gated
   strictly neutral stretch, its ions against its gates' starting charges (R3 1.5 now says so).

   Done after it (22 Nov), Neon's plan with one change:
   - The check first: bench/combos' neutralG, which fails on R3 1.7, R2 3.1 and R3 1.5 with their
     fixes undone, and found two more (a regrid's gate start charges 1e-5 off; a sheet-charged
     face's cold start).
   - The owner is the assembled Gauss row itself, read once a round (src/neutral.js), not shared
     term functions: a projection that holds φ̂ or moves it linearly needs the row's non-ion part
     and its slope, not the row at every trial φ̂, so one assembly does, and nothing is written
     twice. A face's D enters as its slot, synced to its own row by the carry.
   - The carry's projections: a gate's hardware, φ̂ jumping at fixed c, the gates' charges kept (R3
     1.5 fixed, its warning kept as a net); a regrid, by stretch: where every ion is closed, φ̂ held
     and c in the charge direction, a held gate's total kept (Neon's rule: the net charge is kept,
     so the capacitors' is); where one is open, φ̂ at the carried levels and a held gate's charge
     read off the state. The change from Neon's plan is that second half: φ̂ held where the ions
     are open kept the coarse grid's error in the gate's charge, 8e-6, where by the levels the
     regridded state is the fine grid's to 1e-8. And φ̂ at the levels in a closed stretch never
     converges: the amounts and the small capacitors chase each other at 1e-5 a round.
   - The transient's cold start projects a box beside a capacitive face or contact that it left off.

   Done after Neon's fourth review (22 Nov), which read the owner and found six holes in it:
   - The joint solve (src/neutral.js, solveNeutral): each node's φ̂ (at its levels, or jumping at
     its concentrations), each closed stretch's level shift and each floating capacitor's voltage
     in one Newton solve. No switch by screening, no 'charge' mode in the producers, no held
     gate's total kept by hand: a regrid's face D lands on the fine grid's own history to 4e-9.
   - Capacitive contacts are capacitors like gates (start charge read off the state at the old
     drive; a hardware change keeps the plate's charge, the box jumping). The groups capacitors
     join are built over every strictly neutral box, with their plates and faces, and each step
     checks each group's start.
   - The cold start: capacitors uncharged, c0 kept (φ̂ jumping), a sheet's counter-charge in its
     edge boxes, closed amounts and floating charges held, never put back.
   - A capacitor's change beside a face to a metal: the metal's sheet starts from its charge before
     the jump (Neon's design, its closing reply), so the first step passes the charge between them
     through their terminals, recorded.

   Still open:
   - A group beside a diffuse layer isn't checked at the step (its rows read the layer's charge).
   - A drive changed in the same set() as a regrid: the solve holds the drives at their old
     values and the step passes the change, as for a drive alone; not yet in the neutralG property.
   - A step in a strictly neutral region's fixed charge: refuse it, or solve the migration-only
     projection (Neon, R3 DESIGN §3.3).

   Seen while fixing round 3, not reproduced yet: check()'s conservation floor and a steady
   solution's Kirchhoff check are sized by the largest area in the grid (src/checks.js
   conservation, src/solution.js), as the stepper's charge floor was (R3 3.1); and the start's
   pinned-contact density check skips a floating contact, though for an equilibrium link the
   density doesn't depend on V (a contact behind R at zeroCharge 20 V is refused only at its
   first solve).


   A second seam, after: assemble emitting steady and storage parts into separate accumulators,
   as DifferenceTerms does, so impedance() no longer takes (X(1e-30) − X(∞))·1e-30. Neon's R2 1.1
   and 1.2 both sat there.

1. **Harden the core** (ongoing), so the porcelain doesn't force refactors.
   - Transients with a capacitive contact beside a held electrode (found by `npm run lawsweep`):
     with a membrane capacitor too, the conserved K⁺ total drifted in proportion to the tolerance
     (fixed in 1.2.7: a contact's current lacked the membrane port's ∂/∂V, found by a cold review
     of `src/boundaries.js`; 28 of the failing transients pass with it too, 87 of 870 fail now),
     and some such devices fail their first step at one tolerance and not another. The sweep's second family shows it plainer: 34 of its 60
     transients fail, nearly all with a capacitive contact passing K⁺ (at I = 0 or behind a
     resistance) on a thin polar skin, most within 1e-7 s of the start, before and after the
     fixes of the cold review alike. Two mechanisms found (18 Oct 2026), neither fixed yet: the start's φ̂ runs
     from the left, so behind a closed or floating left end it sits 55 RT from the gauge a bath
     on the right holds (device 846's first step, carrying the whole film there by damped
     updates across a face, fails); and where c0 leaves a face out of equilibrium (S⁺ 1 and 0.1
     either side, so c_Na·c_Cl differs), a first step from the right gauge meets an exactly
     singular pivot (device 843). Starting at the bath's gauge fixed the first and broke the
     second, so it was reverted. The face start (the two nodes at such a face start in its
     equilibrium, each ion's amount kept) fixes both, the gauge offset being a bystander: 32 of
     the 60 fail now, every one with the capacitive contact on its polar skin, where starting at
     the bath's gauge made it worse (36). Its mechanism: the contact starts uncharged, 57 RT
     from the K⁺ level its link holds, and drives it across the skin in the first instant;
     started with the link at equilibrium instead (and the bath's gauge), its capacitor starts
     charged against a skin whose double layer doesn't exist yet, and the failures moved (817
     fixed, 869 broken, those behind a resistance now failing at t = 0). Either start breaks
     Poisson's equation or the link. A projection at t = 0 (every amount held, the potentials
     free) showed it isn't an inconsistent start but a real transient: the link closes 1.47 V
     within ~1e-13 s, below any step. Settled with the owner: a start asked for that violently is
     right to fail on, if it says why and what to do. It does (f1fcc37), and `solve()` first
     then `advance()` runs all 16 tried. A floating electrode's reactions likewise (Argon's
     law8, below). The sweep's devices no longer ask for it (2 Nov): their capacitors' zero
     charge is measured from the bath's φ (the K⁺ offset was right for a capacitor already charged
     1.47 V, not for one starting uncharged, and port p1's capacitance likewise 0.77 V), so each
     starts near its steady state, and all 90 transients that failed settle and agree with
     `solve()` (the sweep: all 735 compared, none differing); the old devices stay as the tests' fixtures,
     `device(c, { charged: true })`. Starting at the bath's gauge was tried again and again left
     alone: 64 of those 90 settle without it, 62 with it, and it moved the round-off-set F5 below
     by 1e-4 V between listings.
   - A floating double layer on a thin polar skin beside a neutral face (the sweep's second
     family, `capK`/`RcapK`): the skin's absolute level is held only through fluxes the storage
     term dwarfs, so at short steps its φ̂ and every η rattle together by ~1e-12/dt thermal units,
     no concentration moving. The error control counts that: at tol 1e-9 it chases it down to the
     shortest step at ~5e-5 s (28 devices; at 1e-8 and looser they run, and the sweep retries at
     1e-7). And `advance(1e-9)` from the cold start fails its first step even at the default tol
     (devices 812, 823), where `advance(1e6)` runs. An error norm on ln c and fields, rather than
     on φ̂ and η, would ignore it, a change to every transient's steps.
   - A port's window counted each of its nodes' whole box (a cold review of `src/boundaries.js`),
     a half cell past each end: up to 27% off on a coarse grid. Fixed in 1.2.8: each node whose box
     reaches into the window takes the part inside, so a window holds `to − from` on any grid
     (whole-region windows were right already). Still open: at a device end whose contact sets φ,
     the end node's share of a capacitance goes to the contact (0.0175 of 0.0200 F/m² for a
     0.2 µm window on 0.05 µm cells), O(h) there.
   - The time stepper's error control (a cold review of `src/timeSteps.js`): its estimate is
     Milne's device as if past states were exact, and reads low through the integrator's own
     (backward Euler by 1 + 1/(2ω), BDF2 by ~1.2), so accepted steps err up to ~1.3 tol. The
     right scales (1/2; (h + hp)/(3h + 4hp + 2hpp)) were tried: backward Euler's at 1/2 stalls
     runs whose state a start leaves off its algebraic equations (it's the low reading at short
     steps that lets them through), and BDF2's alone moved two MOS tests' trajectories. Also from
     the review: each call boundary costs a rejected step (a clamped last step carries a size
     whose ratio refuses BDF2), the backward-Euler estimate just after a jump runs 170–430× its
     true error, and the end-of-call snap (1e-10 of tEnd) is wider than round-off. The consistent
     start would take the first of these with it.
   - From Fable round 2's reviewers, still open: a trace has one point at a drive's jump (the
     value just before it), so a trapezoid of it miscounts; the docs give the right sum (each
     point's current times the step ending there), and a second point at each jump would make
     a trapezoid right too. And a cold `solve()` of an n-Si photoanode in Fe²⁺/Fe³⁺ at 0 V vs SHE
     fails from its own layout: its pseudo-transient ramp stops at the step floor after 23 steps,
     where the same ramp at 0 V from the layout built for 0.5 V converges in 52 (and 0.5 V → 0 V
     warm in 19 iterations). The start's layout, not the solve: one more of the starts below.
     (Fixed in 1.10: the error norm's floating terminals, `dt0` after a call, port traces, each
     face's production, `photogeneration()` toward −x, the return to level contacts,
     `impedance()` re-solving.)
   - The state audit (7 Nov 2026): the solver's state between calls is one list (`STATE` in
     src/solver.js), and test/state.test.js holds each field a run writes to it or to scratch,
     and a failed, refused or read-only call and a rollback to leave a device doing what a twin
     does, bit for bit. Still open: `step()` in src/timeSteps.js keeps its own list of what a
     failed step puts back (the start-of-step arrays), a second hand list beside `STATE`; and a
     `set()` of a drive to the value it already has still restarts the time steps' order (harmless,
     a first-order step).
   - From Fable round 5 (6 Nov 2026), still open:
     - Any grid option (`ratio` alone, `minCells`) drops the default grading from hmin = L/1000:
       a uniform 8-cell grid, unwarned. Documented, but the docs invite tuning `ratio` alone;
       filling the defaults per field would move every result of such a device.
     - The grid check can't see an error that doesn't shrink as hmin and hmax halve (a port's
       window falling between nodes: 20% at 9 nodes, 0.06% change on refinement).
     - Small acceptances: unused keys in `fixed` (a check broke face and port reactions, whose
       `fixed` is read through their sides), a bath material's temperature unwarned; a strictly neutral region's unbalanced `c0` refused
       at the first solve rather than at construction.
     - A regrid loses 7 ppm of a closed stretch's amount to interpolation (documented). (The
       2e-12 by which an `advance()` after a failed `solve()` differed from its twin's was the
       twin's: it lacked the two `set()` calls of drive around the solve, each of which restarts
       the time steps' order. With them, bit for bit the same; the state audit, 7 Nov.)
   - From Fable round 4 (6 Nov 2026), still open:
     - `advance()` at a `tol` its error estimate can't reach now stops where a step can't move
       the clock, `converged: false` (it stalled at a standstill, `converged: true`). The floor is
       the φ̂/η split's round-off at ε > 0 on steps far below the dielectric relaxation (~1e-7,
       measured, docs/numerics.md); the charge rows remove that mode at ε = 0, and the same
       transform at ε > 0 nodes on short steps might lift it. A run with a small `maxSteps` can
       still end crawling above the clock's resolution (a sphere at tol 1e-7, 4000 steps:
       `converged: true`, current 3× off); a stop on an error that won't shrink was tried and
       broke runs through a jump, whose error is self-similar.
     - The Jacobian test now allows for the differences' own round-off (ε of the row over h), so
       a lit device checks clean, at a light dim enough to check its hole rows (1e19 cm⁻³s⁻¹; at
       1e22, 34 rows would pass a 10% error). A brighter device's rows need a larger step or an
       analytic check.
     - A contact's fluxes are still read from the end box's residuals, a small quantity as the
       difference of large ones on a short step; only a closed contact (every link blocked) now
       passes exactly nothing. Where some links are blocked, a blocked species' slot can hold a
       charge row or a law's row with the linked ions' fluxes in it, so the simple fix (skip it)
       broke a diffuse layer charged through a contact; the link's own law at the solution would
       be exact. With the closed contact fixed, blocking electrodes resolve Z down to 1 µHz, and
       `impedance()`'s "Z isn't resolved" warnings have no device in the tests that triggers them.
     - The lattice model's K at a full lattice cancels to zero (c − c²/c_max); carrying the vacancy
       fraction would keep it. No device effect found (a crowded double layer at 60 kT agrees with
       Kilic–Bazant–Ajdari to 1e-4).
   - From Fable round 3 (5 Nov 2026), still open:
     - `impedance()`: Z measured at a port with no floating terminal gets neither cuts nor the
       region-shift split, so the MOS capacitor's channel port is unresolved (warned) where its
       gate contact resolves to 1e-6. (A floating terminal's δV now carries its readings' error,
       and the floating responses are split too where a reading is unresolved: the Pt island
       behind a bath's 1 Ω·m², once 0.3% off at relError 5e-12, reads Z − R within 4e-10 of the
       held device, and behind 1 kΩ·m² within 2e-7 at relError 4e-6; a MOS gate behind 1e3–1e6
       Ω·m², read at the back contact, Re Z − R within 2e-7 of the held gate's 2869 at 10 nHz–1 mHz;
       Sol's fourth pass.)
     - `stuckReactions` now checks a bulk reaction node by node where it alone changes a
       participant that can't move; a port's reaction out of equilibrium both ways across its
       window, with an immobile surface species, is still checked only for one sign.
     - A strictly neutral (ε = 0) bipolar membrane in reverse bias past ~1.2 V fails cold, warm and
       from a converged transient's state, "lost 12 of its digits" at the CEM | water face, and the
       message's remedy (anchor the region) doesn't help; a 37 s transient and the resolved ε > 0
       model agree on the state (−6.495 A/m²). The neutral junction's co-ions must fall to e^{−45}.
     - The start fills a trap level half full: `trap().c0For()` gives each layer's equilibrium
       filling, which a pn junction's cold solve needs; the start could fill immobile states itself,
       node by node, from the carriers it lays out.
     - `advance()` fails at its first step right after a `set({ T, materials })` jump far in T
       (150 → 500 K at −5 V; 150 → 300 K at 0 V): the minority carriers carry over at e⁻⁶⁸ of the
       new equilibrium (e⁻⁶⁹ at −5 V), so their levels must rise ~68 thermal units within the step, and Newton's divergence
       guard stops the climb (a guard that let shrinking updates go on rescued only one case of
       five). The failure says so now, and jumps of 50 K, or a `solve()` first, get there.
     - Demo pages at a phone's width: the membrane page's three dot-panel captions overprint, the
       myelin level diagram's 40 node lines crowd it; the rde page's Koutecký–Levich axis goes
       negative at E° + 0.1 V; the redox page could show its no-iR separation beside the readout's.
   - Starts needlessly far from where the physics goes (a cold review of `src/start.js`; each
     the first step absorbs, listed under Known in its contract): the gauge from the left; an
     insertion host's shift moving contact-held levels (a two-host battery 137 mV from both held
     contacts; behind a Cu collector at I = 0 it starts at 0 V for 4.09); a contact-fed carrier given a c0 (a resistor at rest
     0.2 V off both contacts, 17× the steps); a terminal behind R at 0 V; a reaction with a
     neutral level-held participant passed over beside an electrode (0.68 V off). Their common
     shape: the start's gauge and levels should come from what holds them, not from the left.
     (One more was found and fixed in 2026-10: a non-ideal solution beside an electrode skipped the
     balance an ideal one gets, so a resolved crowded double layer at −0.6 V failed every cold
     solve; it now starts balanced and converges cold.)
   - A capacitive electrode with adsorbates, at no current from an uncharged start 1.65 V from
     where its surface balances (Argon's law8): steady solve and transient both fail, and now
     say why and what to do (hold it at a voltage, solve, then let it float), which works.
   - Floating conductors in $`\eta`$ form. A floating metal is fixed by the mixed form, but the same
     cancellation awaits a floating mixed conductor or semiconductor region held only by face
     reactions or recombination, in steady state (in transients, storage anchors it): a fast
     species across many segments, with a mode only weakly held from outside. Candidates are
     a Grassmann–Taksar–Heyman-style elimination that carries leakage separately, or the mixed
     form for the fast species only. (A failed steady solve already says where it lost its digits.
     Newton's refined solves, which recover a MOS inversion layer's lost level, may cover it.)
     `npm run stress` found the same thing in a bipolar stack's floating base (a p⁺ base between
     light n layers, its holes held ~1e14 weaker than they move within it: 25 digits lost even
     at 10 mV), which refined solves now reach, and in a closed Fe³⁺/Fe²⁺ cell driven by a
     current near its limit, where Fe³⁺ falls 20 orders below Fe²⁺ (33 digits floated), now found
     by held solves instead. Past ~1e-22, a
     solve converged to the wrong level outright (a region behind face conductances 1e-24 of
     its own: the plain solves gave no update along its level at all, so nothing stalled). Such
     regions are now deflated: each one's balance summed over it, in which its own fluxes cancel
     exactly, is bordered like a conserved amount, once the plain solves converge with terminal
     currents that don't add up (and since 2 Nov, one fed only by its faces' reactions, and a
     floating metal: F5 below; since 7 Nov, one with bulk reactions, each node's rate summed).
     Not covered yet: a population held that weakly within a region, or one fed by a port. Measured (2 Nov): an islands round tried and dropped
     isn't the sign of a wrong steady state. It happened 19 times in the law sweep, every one on a
     device that agrees with its transient, none in the suite, and in the stress run only on
     semi 220 (below). Their currents were 1e-18 to 2e-17 A/m² from adding up (equilibrium's
     noise), against 1e-15 to 1e-9 on the floating bases that were wrong. So since 3 Nov a steady
     solution whose terminal currents miss Kirchhoff says so (since 4 Nov, by more than 1e-15
     A/m², 1e-6 of the largest, the floating terminals' gaps and 1e-12 of what the species carry:
     numerics.md); and since 5 Nov the round is also asked for by each island's own balance, and
     a piece the round couldn't sum out, left unbalanced, is warned of.
     Open: warned, not solved, a base whose end piece has a port (a collector contact carrying
     1 mA/m²: 0.135 V off). And a floating base with the species listed holes first, below
     300 K (the three-layer base at 250 K, −0.3 V: 0.246 V off, its currents 5.1e-14 and
     −6.3e-14 A/m², warned): its round walks the level down a thermal unit an update, as it
     should, and midway (the sixth update there) the factorised matrix's plain solve returns
     2e24 at a node of the emitter's minority holes between the cut and the junction, GMRES on
     it doesn't converge, and the round is dropped. The same device listed electrons first
     solves. It belongs with the listings item below. At 150 K (leaks 1e-42 mol/(m²·s)) a base
     at −0.4 V is left 0.16 V short (its round stops partway down; why isn't traced), and some
     end microvolts off;
     both warned. (Fixed 9 Nov: a base behind further layers, n/p⁺/n/n⁺, silent at 200 K and
     0.16 to 0.26 V off, and the three-layer base at 150 K, 0.1 V off under the imbalance's
     floor: the sum reaches through the layers between, a layer it passes isn't summed itself,
     each sum's row is its middle node's, and the imbalance is read against the cuts' fluxes;
     numerics.md. A survey of 360 such solves, ~/.cache/driftlet-claude/base4/survey.mjs, had
     61 silently off and has none.) (Fixed 7 Nov: a floating base whose
     piece has a bulk reaction, hunt A, GaAs at 250 K with band-to-band recombination, which
     converged silently at −4e-10 V where it's −0.1477: each node's rate is in the sum now, and
     a cut passes an end piece's reactions; numerics.md. Not covered: separateRates takes face
     rates only, so a fast bulk reaction shared by two islands could leave the bordered system
     nearly singular, as F5's fast couple did; hunt A's base has two such islands, e⁻ and h⁺,
     sharing a slow one, and solves.)
   - Kinetic steady states finer than a fast couple's round-off. An electrode at open circuit
     with a fast couple (k0 = 1) and a slow one, fed through a slow transfer at a face (k0 1e-17),
     settles at a balance of ~1e-18 mol/(m²·s), which needs the fast couple's net rate to 1e-20
     of its gross. Steady solves converged silently to states 20–80 mV off (and different in
     each listing of the species), their terminal currents not adding up by ~3e-13 A/m², which
     is about the fast couple's round-off. Found by round 11's fuzzer (its F5). Since 1.1.1 (the
     contact and flat-stretch rows keep a level's low part) its two listings agree to round-off
     and its ledger closes to 4e-16, where they were 4 mV apart. But it isn't the steady state:
     advance() from the solution (0.274922 V) moves at once (the fast couple's rate falls from
     1e-18 to 1.6e-21 mol/(m²·s) within 1 ms), drifts ~20 mV a decade from 1e9 s, and settles by
     1e15 s at 0.195442 V, unmoved to 1e18 s; a cold advance() ends there too. Fixed (2 Nov):
     the Pt's electron balance and the B⁻ pool's, summed over their nodes to their rates alone,
     are pinned like islands', each rate in one sum only (numerics.md), and the solve lands on
     0.195442 V to 1e-9 in either listing, and driven at ±2e-16 or 1e-14 A/m² too. Not covered:
     a time step so long that storage falls below the same round-off (here dt ~1e20 s, from
     1e18) lands on the old state, islands being steady-only; and such a pool or metal fed by a
     bulk reaction or an electrode port. The law sweep's third family (108 such devices: fast and
     slow couples, the pool's transfer, I = 0, a set current or a resistance, either listing) had
     20 settling off before; all agree now (a terminal behind a resistance taking the sum's ∂/∂V).
     Seen on the way, a different case: the stress run's semi 220 (three GaAs layers, no
     recombination, its electrons degenerate past Boltzmann, warned) gave two currents by path at
     −0.2 V, −1.73e-13 A/m² warm and −3.17e-11 cold, unwarned. Fixed (3 Nov): it's a floating p⁺
     base whose holes leak through the n layers at ~4e-15 A/m², ~1e-24 of what they conduct at
     the junctions. An independent leak integral over each state's φ put the warm state's leaks
     in balance to 3e-7, while the cold one's base, at the emitter's level, drained 8.6e-13 A/m²,
     which was the gap between its terminal currents. Its islands round, summed at the faces,
     diverged. Summed out to the n layers' bottlenecks (numerics.md), it now lands on the warm
     answer from cold. The same round was cut off at its dozen updates while still converging,
     walking the base down a thermal unit at a time; those updates no longer count (test:
     test/newton.test.js, a floating base checked against its leaks' balance at four biases).
     Not covered: a base whose neighbour isn't the stretch's end piece. With a fourth layer
     (n/p⁺/n/n⁺) the island is still summed at that face, and cold solves at −0.4, −0.2, −0.1
     and 0.1 V miss Kirchhoff (−0.2 V: 3.29e-11 and −3.83e-11 A/m², where warm gives 1.96e-13);
     they were warned; since 9 Nov the sum reaches through the fourth layer, and they solve.
     Nor a neighbour that has a port. Nor a contact passing holes
     only (n/p⁺/n behind ohmic(0, ['h+'])): from cold at 0.3 V it put −1e-20 A/m² on that terminal
     and 4.25e-15 on the other (true: ±3.59e-15), and with the fourth layer cold solves at −0.2,
     −0.1, 0.1 and 0.3 V and warm ones at 0.1 and 0.3 V miss too, each confirmed by a transient to
     1e18 s. All warned (since 4 Nov: the warning had passed a sum on one terminal as round-off).
   - A decay's tail in a trace. `tol` is absolute, in thermal units, so a step of ΔV is followed
     only to about tol·V_T/ΔV of itself: a double layer stepped by 5 to 250 mV erred in its
     current by 3 to 8 × tol·V_T/ΔV of the peak, and at the default tol the 5 mV step's current
     changed sign past 5.4 RC (the Fable 5.2 docs review, 9 Nov 2026). The docs give a tol for a
     tail (e^−n·ΔV/(80 V_T), within 8% at all three). Candidate: an error control that also
     weighs each terminal's stored charge against what it has passed, so that a small signal
     isn't left to the user's tol (it would change every run's steps: measure the cost first).
     (Its charge, Σ I·dt 22% short there, is `trace.charge` since 9 Nov, exact.)
   - Steady solves that depend on the order species and reactions are listed in. Of random
     reacting devices (round 11's fuzzer: ions, complexations, couples on metals, I = 0 and
     resistive terminals, k0 log-uniform 1e-25 to 1e3), 1.3% of the cells and ~10% of those with
     electrode ports converge in one listing and not another, or to different states, and ~1%
     have a steady state (a long advance() settles to it) that solve() doesn't find. Four causes
     are gone: a pivot left at round-off without being exactly zero, affinity rows picked by
     listing (cells 1.8% before), conserved combinations eliminated in listing order (a trace
     species added to an abundant one, its amount lost to round-off; 1.1.1, from Argon's
     review), and a level's low part dropped at contacts and flat stretches (F5 above, 1.1.1). What's left: Newton paths that differ with the order of a
     block's unknowns (one path crawls or lands on a dead state, which the stuck-reaction warning
     flags); slow port reactions (warned, below); weakly held regions that lose all their digits;
     and devices with more than one steady state. A pair of port couples settled either near the
     fast one's equilibrium or far from it, its reactant at 1e-18 mol/m³, both passing every
     check; which state a solve finds then depends on its path, and advance() from the start
     finds one of them. In that pair, the second state is unstable: kicked by 10 mV and advanced
     with dtMax capped, it left for the first within 2e5 s. check() solves the device again with its
     species reversed, and fails where the two differ (on 2600 of the fuzzer's devices: every
     real disagreement, no false alarm).
   - A stability check after a steady solve. Tried: the sign of the time step's Jacobian
     determinant, J + σM, along σ = 1/dt from 1e12 to 1e-12 per s (it changes at each real growth
     rate), from the block factorisation's pivots and the bordered rows' small system, with each
     conserved amount in place of one of its rows to remove its zero root. It found the unstable
     port state (between 10 and 1e7 s; ~3000 s by a kick), but 3 of its 4 flags on 600 fuzzer
     devices were stable by a kick, and 9–40% of devices came out undecided: the bordered rows
     (a floating terminal's current response) lose their sign to cancellation at long steps.
     Also tried: power iteration on kicked backward-Euler steps at fixed h (a growing mode
     multiplies by 1/(1 − gh)), right for that state (1.032 per step at h = 100 s, 1.031
     predicted) but with spurious growth of 10–1e4 per step on many port devices, unexplained.
     Either needs its noise understood before it can fail a check; meanwhile numerics.md says
     how to kick a state by hand.
   - Bulk reactions too slow to see. A face reaction that alone changes a combination now stands
     on its affinity in steady solves, at any k0; a bulk reaction can't, since a zero integral
     of its rate doesn't make it zero at every node. K⁺ + A⁻ = KA in a cell whose only reservoir
     is a K⁺ contact solves at kf = 1e-12 but fails at 1e-20 and at 1e-30 silently keeps the
     starting composition (K = 0.1 for 7.5), now with a warning that it isn't a steady state.
     Such a device is at equilibrium throughout (one reservoir drives nothing), which a steady
     solve could recognise and hold, every affinity zero node by node; but a field set up by
     other species' currents could drive a closed charged species round a reaction cycle, so
     the rule needs care. Electrode ports' reactions likewise (spread over a window).
   - Strictly neutral regions on very short steps, where storage dwarfs fluxes. (Round 11's GITT
     agent met it from the other side: two insertion hosts with a background charge, a current
     switched on, at tol ≤ 1e-7 the first step's error control shortens it below ~1e-7 s, where
     Newton no longer converges, and the run stops at t = 0; tol 3e-7 runs.) Interior nodes and
     the edges of neutral faces are solved in better-conditioned unknowns ($`\hat\phi'`$,
     $`\eta - z\hat\phi`$), but at a neutral face the two edge nodes still pass every species' flux
     through its own unknown, which carries storage-sized values, so round-off returns at ~1e-5
     thermal units on steps of ~1e-11 s with 3 M against 1 µM (a sharp junction between such
     solutions now starts, even at tol 1e-6; `npm run hard` keeps it so). The fix may be a
     supernode: each region still supplies its own edge block, and the interface code merges the two
     edge nodes and the face's unknowns into one larger block, so every quantity at the face couples
     to every other. Each species' two edge balances can then be summed, with the face fluxes
     cancelling exactly, into a row whose outside couplings reach only the nodes just beyond. The
     block at each face grows (about threefold) and so do its off-diagonal neighbours; the rest of
     the matrix is unchanged, and the layout stays block-tridiagonal. It may simplify the face code
     generally. Tried as merging alone (each face and its edge nodes factorised as one dense
     block, LU pivoting across all three): it made static pivoting unnecessary (GaAs stacks
     without recombination solved without it), but partial pivoting across the merged block
     picked pivots from other nodes' rows for a nearly empty carrier's column, and Newton diverged
     in most MOS capacitors (462 of 500 stress scenarios); it also hid where a weakly held region
     loses its digits from the conditioning warning. So the summed rows would need a structured
     elimination at the face, not generic pivoting.
   - Fewer Newton iterations where the benchmarks show many: the pn transient (about 1000
     factorisations for 100 ns) and large warm jumps. The pn transient's 295 steps are set by
     BDF2's error control (about 50 a decade over six decades), and almost every one takes three
     updates, 3.5e-3, 7e-6, 3e-11 thermal units, the third predicted by the quadratic estimate
     to within 10%. Stopping at the second saves a third of the factorisations, but leaves each
     step's equations solved only to ~3e-11, and conservation then drifts ~1e-11 a step (the
     tests hold it to 1e-12). A chord update with the second's factors (a residual and a
     back-substitution, no factorisation) would keep it exact, but saves little: assembly costs
     about twice a factorisation here (~100 µs against ~55 µs for the pn diode's 265 blocks of 3,
     and for n-Si | KCl's 315 of up to 5), so ~10% of a transient. Assembly is the larger cost.
   - A port's window reaching an end node whose φ its contact sets leaves that node's half box
     out (its capacitance and exchange), so a porous electrode against a bath has h/2 less
     electrode than drawn: an O(h/L) shortfall in its capacitance and impedance (1.25e-3 at 400
     cells, found by a porous-electrode agent). Including it means the contact supplying that
     box's counter-charge, and its current reading counting it.
2. **Porcelain, for one-shot demos.** An optional layer, the `driftlet/kit` subpath export
   (started, with the vacuum-level alignment helpers), still dependency-free, that writes plain
   specs, so users and LLM agents start from something
   correct. It never hides a physical choice: everything it produces is ordinary, inspectable
   spec data.
   - Done: `build({ library, stack })`, the device as drawn (contact, layers and faces, contact)
     into the plain definition; `ohmic()`, `bath()`, `layer()` with doping as donors and
     acceptors; reactions written as equations (`'Ag+ + e- = Ag(s)'`), with a label for a species
     on both sides of a face; half-reactions as plain data, and `level(sol, half)`, the redox
     level they imply (with `SHE`); `traces(sol)`, the level diagram as plot-ready data, and
     `driftlet/plot`'s `bandDiagram(sol)`, an SVG string; `live(def)` for sliders (merging
     changes made during a solve, warm starts, ramps from the last good state, Worker-backed);
     `describe(def)`, a readable summary with unit-slip and resolution warnings; a small data
     library with a source on every entry, checked against independent tables (13 aqueous
     ions, Si, Ge and GaAs, seven metals); an agent guide, [`llms.txt`](llms.txt), with templates
     and common mistakes, whose code (and a slider page, against a stub document) runs in the
     docs test. See [the kit](docs/kit.md) and [data](docs/data.md).
   - Built-in and open-circuit potentials in `describe()` (they need a solve).
   - Whether the stack form should replace `regions` and `interfaces` in the plain definition
     itself, once it has been used for a while.
   - Half-reactions as vacuum-alignment anchors.
   - Done: voltage-gated channels, Hodgkin–Huxley gates as unknowns of a membrane face (and of
     each node of a membrane port's window) scaling its permeabilities or conductances
     (`hodgkinHuxley()` in the kit), probes at a face, a port's `bath`; the
     [action-potential](demos/axon.html) (with Frankenhaeuser and Hodgkin's periaxonal space),
     [propagation](demos/propagation.html) and [myelin](demos/myelin.html) demos; validated
     against the space-clamped HH equations, the linearised ones (the excitable membrane's
     impedance), the cable equation (18.73 m/s) and a compartmental cable (saltatory
     conduction). Still to come: Markov-scheme channels (states with transitions, for channels
     that HH's independent gates can't describe), calcium (its channels, buffering, and the
     Ca²⁺-activated K⁺ channels behind bursting), and an extracellular space with its own
     resistance along the axon (two lines, ephaptic coupling between neighbouring axons).
     Devices with dozens of ports are slower per step than their grid suggests.
   - Electron transport and the protonmotive force: a mitochondrial or thylakoid membrane as a
     face, with V_H⁺ the quantity that couples the respiratory chain to ATP synthase.
   - From the corrosion agent's review:
     - Done: a cross-section that varies along x (`geometry`), spherical and cylindrical shells
       included. With electrode ports, that's the radial Evans drop, A = 2πr·h(r).
     - Done: rates at a hypothetical level, `polarization()`, for Evans diagrams.
     - Done: an electrode's surface (Langmuir coverages, reactions on bare metal), for passive
       films, and a capacitance through a port's window (a double layer, a gate along a channel),
       which carries a transient through a film's passivation, and a thin-film transistor in the
       saturation demo.
     - Done: a bath's SHE level against its terminal in `describe()`. An agent held the iron at
       V = 0 beside a bath referenced on Cl⁻ (NBS data), 1.37 V anodic, and lost time finding out
       why. Still to come: the same in failure messages, and each terminal against SHE in a
       solution.
     - Done: cheap mid-run parameter changes. `set()` of a reaction's k0 rebuilds the device
       (about 0.4 ms) but now goes on at the step size it had reached, where it used to restart
       from a tiny step (7× the steps in a cyclic voltammogram changing k0 every frame).
   - Done: sources as data (`pulse`, `square`, `triangle`, `ramp`, an `injector` port,
     `recombination` from a lifetime), probes in a transient's trace, `units` back out of SI.
   - Done: an initial profile as plain data, `c0` tabulated against x, so a transient starts
     from an already-injected or relaxing state without simulating how it got there (and
     without a back door around the conservation bookkeeping: it's the definition's start).
   - From the perovskite agent's review: the pieces every demo rebuilds.
     - **More checks.** `check()` covers each species' ledger, conservation and grid
       convergence; still to come are limits that need a device family's theory (V_oc against
       detailed balance, an ionic RC time), a transient's time step against the device's fast
       and slow time scales, and the second law's ledger (free energy in = dissipation ≥ 0, which
       the test suite checks) as a reported item.
     - Done: energy-up band diagrams in `driftlet/plot` (`energy: true`) for electrons and holes,
       with species voltages kept for what that picture can't show: ions, and equilibrium between
       species.
     - **Inputs in the source's vocabulary**: a trap level against a band edge (SRH takes n₁ as a
       concentration today), and `describe()` reporting derived quantities to check against intent
       (built-in voltage, n_i, band offsets, the equilibrium levels against vacuum).
     - **Frozen species** as a named switch (today D → 1e-40) for separating time scales.
     - **More curated devices** with a reference result, as `perovskiteCell()` has IonMonger's
       scans: a silicon solar cell (with silicon's absorption spectrum as data), an organic
       bilayer, an electrochemical cell, each where an independent code or measurement gives
       numbers to check against.
     - **Failure messages that prescribe the next step** (continuation, a smaller first step), as
       the lost-digits warning already does.
   - Done: worked examples: the Daniell cell's voltage ledger, and the electrode pages (Frumkin,
     potential steps, the rotating disk, Randles, corrosion, lithium plating, Chidsey's
     monolayer, hydrogen evolution), each checked against theory solved apart from the solver in
     the page's own script (test/demos.test.js). A ledger's division is specific to its
     device family, so it belongs in a worked example, not the library.

## 2.0

What a major version would clean up. Anything 1.x keeps only for compatibility goes here when it
is kept, so it isn't forgotten; when keeping something would take a contortion, the answer is a
2.0, not the contortion. With no known users yet (no issues, no dependents), a 2.0 is cheap: cut
one when this list is worth it.

- Refuse the 1.17 defaults that 1.18 still accepts on fields that mean nothing there: `V: 0` on a
  port that passes no current, and `fixedCharge`, `velocity` or `mixing` of 0 on a conductor.
- Frozen species as a named switch, replacing the D → 1e-40 idiom (additive in 1.x; 2.0 can refuse
  the idiom).
- With any rename: an optional `version` field in definitions, and errors that name the renames
  ("`foo` became `bar` in 2.0"), since agents copy old definitions long after the API moves on.

### The shape: simpler and harder to break

The API list below is the easy half. The other half is internal: 2.0 is the time to reshape the code
so whole classes of fault can't happen, rather than fixing them one review at a time. Notes for
that, not a plan yet.

**Where the faults come from.** The outside reviews of November 2026 (rounds 3–6, about a hundred
confirmed faults fixed, many of them silent wrong answers) sort into four sources, roughly by size:

1. **Features meeting features.** The biggest pile: a port turning off the immobile check of a
   species it doesn't touch, a floating gate behind a resistance, a Hodgkin–Huxley leak on the gated
   ion, kit builders merging their arguments the wrong way round. Every physics feature has its
   analytic test, but tests rarely put two features together, and the pairs grow with the square of
   the feature count.
2. **State that goes stale.** `Device` is a single mutable object holding many caches, and every
   method that changes it has to remember to invalidate, snapshot or roll each one back. A new
   feature adds a path that forgets one: check() passing a solution from before a `set()`, a refused
   `impedance()` using up the pending transient, a throw inside `solve()` leaving the clock moved,
   island flags left out of the snapshot. The 1.19 solution stamp turned one "remember to" into a
   mechanism.
3. **One law written twice.** The steady solver and the time stepper each assemble some boundary and
   interface terms their own way, and the copies drift: the steady pin read the edge levels where the
   other path used the plane (a current 3× off), the collinear-reaction skip lived in two files, and
   the start ignored grades that the regions applied.
4. **Numerical traps.** Integer arrays, sums losing digits, fixed-step quadrature, round-off cut-offs.
   These are tricky in any code; the fixes (compensated sums, graded grids, self-checks) are now
   shared helpers.

**Directions.**

- **One state, one transaction.** Everything that changes the device (`set`, `solve`, `step`,
  `advance`, `impedance`, …) runs as begin, then commit or roll back, over one plain state object.
  Caches are derived from it and stamped, never invalidated by hand, and a solution knows which
  state it came from.
- **One definition per law.** The steady and transient paths, the start and the readout all get
  their face, port and pin terms from a single assembly (the lawRows contract, extended), so a law
  can't mean two things.
- **Normalize once.** Everything after normalization reads only the canonical model, never the raw
  definition; the kit returns plain, deep-copied data.
- **set() as three verbs** (from a question Mark put to Neon, 21 Nov 2026). set() does three jobs
  with different answers to "what's held fixed": a drive (V, I, R, light; the state as it was,
  exact), a regrid (nothing physical changed: c, φ, amounts and charges kept), and a hardware edit
  (a thickness, a doping, ε, a capacitance, a port added). Most carry faults were hardware edits
  where two conventions met. A rule for those: each part keeps its natural variables, those its
  boundary conditions fix, the ensemble solve() already uses (a region open to a bath its levels,
  a sealed compartment its amounts, a floating terminal its charge, a held one its V), then the
  state is projected onto neutrality by the single owner (## Next, item 0). It's the carry, not only
  a guess, for solve() too: with sealed amounts, floating charges or two steady states, the carried
  state picks which one. A hardware edit mid-transient ends the transient (clock kept, history
  cleared, marked restarted) rather than inventing a t = 0⁺.
- **Invariants as the gate.** Before any of this, a randomized suite of properties that hold for any
  device: a steady solve agrees with a long time-stepped run; Z at low frequency agrees with dV/dI
  from two steady solves; one `advance(t)` agrees with many small ones; a definition's key and
  listing order don't change the answer; check() refuses any stale solution. This suite covers the
  combinations (source 1) that reviews only sample, and it is what proves a refactor kept the
  physics.
- **One owner per quantity, from an inventory of the copies** (20 Nov 2026; the inventory and its
  probes are in the maintainer's notes). 17 quantities are computed in more than one place. Two of
  them disagreed when they were found: a capacitive port's charge carried through a new cross-section,
  and a region changing kind. Seven readers decide which side of a drive's jump the state is on, by
  four rules. A terminal's current has six copies, which agree except at a strictly neutral end on
  short steps and right after a drive-only `set()`. Nine predicates decide what a contact's φ law
  sets. In the shape this suggests:
  - Facts are decided once and stored as data: a role record per contact and port (sets φ, holds a
    charge, closed), a transport descriptor per region, and each law's circuit reading.
  - Each quantity is computed once per state, by the kernel that owns its law, with its derivative,
    into stamped arrays: charges, rates, a segment's flux and its ∂, face affinities, the terminal
    currents as solved. Readers never evaluate a law again.
  - The state carries the sources it was solved at, `{ t, before }`, one record in place of
    `sourceTime`, `sourceBefore`, `readV` and `readI`.
  - The order: the invariants first (bench/combos), then the role records (mechanical, checked by
    tools/fingerprint.js), the sources record together with the one transaction, the per-state
    arrays, and carry last, re-expressed as asking the new solver what it must match.
- **What each feature costs.** An honest audit of what the demos and docs use. Something rarely used
  that touches every path is a candidate to move out of the core.

From a fresh-eyes read of the docs and types (19 Nov 2026; its report and probes in the maintainer's
notes), ranked by how much each trips an agent writing its first device:

- **One sign rule for terminal readings.** `terminals`, `ports[k]`, `trace.terminals` and the I
  drives are into the device; `sol.current`, `contacts.*.current/flux`, `trace.current/charge` are
  toward +x; `terminalVoltage` and `trace.voltage` are two names for V_right − V_left. The pn template
  reads one current twice with opposite signs. What belongs to a terminal is into the device, what
  belongs to a place in x is toward +x; drop the duplicates that `terminals` already hold.
- **Grid options keep their own defaults.** Given any grid option, the graded defaults are all off,
  so `{ hmax }` alone makes the end cells 50× coarser, unwarned (a 1.4% limiting-current overshoot at
  an electrode, by device.md). Each option left out should keep its default.
- **`bath()` takes an options object** (`bath(c, { reference, V | I | R, offset, material })`): the
  second argument is today a reference, a half-reaction or a drive, and one ordering once built a
  bath without its reference, silently.
- **One x frame.** A port's `from`/`to` (and `alloy`'s fraction) are metres into the region, every
  profile (a port's own `area` among them) is device x. Device x throughout, or a key that names its
  frame (`window`).
- **Names that don't collide:** `D` is a diffusivity in `sol.D` and a displacement at faces, contacts
  and plates; a face's flux is `N`; `charge` is stored in `sol.charge` but passed in `trace.charge`.
  `displacement`, `flux`, `passed`.
- **`advance(tEnd)` reads as a duration.** `advanceTo(t)` can come in 1.x; 2.0 drops the old name.
- **A changed `c0` in a `set()` patch is refused** (today accepted and ignored on a carried state; a
  warning is the 1.x half), pointing to a new Device.
- **`set()`'s merge-or-replace rule is visible in the data**, not a hidden mark on kit contacts.
- **A port's `zeroCharge` is required**, as at faces and contacts (the docs' first page promises
  alignments are never defaulted; a port's 0 is worth 1.4 V of charge), and sits where theirs does.
- Smaller: `check(dev)` rather than `check(dev, sol)` with two things that must agree; ports
  addressed one way (by name) in the definition, `set()` and the solution; one shape for reaction
  results, `oneWay` (not `forward`) for a one-way rate; T required, or taken from the materials when
  they agree (298.15 K beside 300 K semiconductor data); `'neutral'` as a φ law renamed for what it
  says (`'uncharged'`); the device's internals private, with `dev.time` added in 1.x for the demos
  that read `dev.solver.time`; units in the data tables' keys (`SEMICONDUCTORS` is in eV and cm⁻³);
  `grade: { species, epsr }`; `ocv: { x, values }`; a conductor port's `G` under its own key (S/m²,
  not S/m³); `solve()`'s ramp options named apart from `advance()`'s (`tol`, `maxSteps` mean
  different things in each).

## Later

- **Open from Neon's first review (20 Nov 2026).** The cold start still sweeps left to right: a
  metal away from the contacts takes its level from its left face's reaction only, so one whose
  only reaction is on its right face starts at the running φ̂ (1.6, item 2). A start that finds its
  anchors first (held contacts, held electrodes, gates, responding regions) and runs outward from
  them would make every such case symmetric by construction, where each has been patched as found
  (round 9's mirror, Neon 1.3 and 1.5). And a time step's charge rows assume the step starts
  neutral at an ungated node: the pivot species' storage is left out there, so a start that isn't
  neutral there loses the difference from that species, unseen (how Neon 1.1 lost its ions).
  Taking the old charge into that storage at every node, as gated nodes have it, broke 37 tests:
  understand what relies on the present form before changing it. And advance()'s diagnosis of a
  stop at the time's resolution (the error estimate's floor, "loosen tol", or a step as short as
  the device needs) turns on round-off: a start moved by 1e-11 V turned the one into the other in
  test/transient.test.js's floor case. The block elimination is unstable on a well-conditioned
  system that every counter passes (Neon's 4×4, [[ε,1,0,0],[1,1,½,0],[0,½,1,1],[0,0,1,ε]], κ₂ = 4:
  at ε = 1e-15 an error of 0.125 against GEPP's 3e-16, lost 0, perturbed 0, ratio 1): no device in
  a 30-run sweep reached it, but the conditioning note could report growth (max|C′| against the
  blocks' entries), and a fix needs pivoting across blocks or a twisted factorisation.
- **Ideas from Neon's architecture read (25 Nov 2026).** These are notes, not plans. Neon read the
  source for smells; its counts are its own, not rechecked. Two of its five points are already
  here: storage apart from transport (## Next, item 0's second seam) and set() as three verbs
  (## 2.0). The rest:
  - *Model, State, Workspace.* The Solver carries about 300 fields on one object: physics state,
    history, Jacobian scratch, mode flags and diagnostics together. The split would be an immutable
    model, a plain-value state that a failed solve simply drops, and a scratch workspace per call.
    That is ## 2.0's "one state, one transaction" carried through to the data structures. Mark is
    happy with a rewrite for tidiness alone.
  - *Continuation through folds.* The steady ladder walks one parameter at a time (a voltage, a
    current, a resistance). Where the I–V curve folds back, as Debye–Hückel's spinodal near
    V* ≈ 1.1694 V seems to, such a walk can only halve its stride against the fold: that is the
    Zeno walk. Pseudo-arclength continuation follows the curve's length instead, so it can turn
    the fold, or report that the branch turns. Neon's bolder form replaces the ladder with one
    homotopy, from dark equilibrium at s = 0 to every drive on at s = 1. The ladder exists because
    simpler paths failed on particular devices (0.11, 0.16). One s ties the whole path to its
    stiffest part, so the moonshot is to try it and see where it breaks.
  - *Fewer rescue layers.* Compensated levels (u + uLo), charge rows, the difference form with GMRES,
    islands and the steady polish each answer a loss of digits. Neon suggests a pivoted or twisted
    block factorisation, or the conductors' mixed flux unknowns used more widely, could replace
    some of them. Caution: as I understand them, most of these losses happen when the residual is
    evaluated, before any factorisation, which pivoting can't recover. Its point that fluxes and
    reactions are evaluated in several places, in more than one layout, is ## 2.0's "one owner per
    quantity" again.

- **Open from Neon's second review (20 Nov 2026).** In a transient, check()'s conservation item
  reads only unreacting stretches, so where every species reacts it says "nothing to conserve",
  though the conserved moieties (the laws' combinations) are there to check (R2 2.3). And
  _chargeRows' assumption (above) has now been reached through two doors, Neon 1.1 at the cold
  start and R2 3.1 through set({ grid }). Each is closed where it entered, but the charge rows
  taking the old charge at every node would close the room.

- **Open from bench/combos' neutral family (20 Nov 2026; repros in the maintainer's notes).** A
  strictly neutral region whose c0 carries more than its held gate takes within a volt puts φ̂ many
  volts from the gate's V (8–11 V, 300–430 thermal units, at a gate of 7e4–2e5 F/m³). Two things
  fail there. A divalent salt in an ungated region left of the gated one fails cold, saying nothing
  more specific: the left-first sweep starts it at φ̂ = 0 where its steady state follows the gate,
  and monovalent salts recover from that start, where Ca²⁺ doesn't. And a gate floated after a
  held solve fails from the held state, an exact solution of the floated problem (residual 4e-14),
  and converges with a third of the excess. The family keeps its gates within 0.5 V; the start
  from anchors outward (above) is the fix for the first.

- **Open from bench/combos and round 7 (20 Nov 2026; repros in the maintainer's notes).** A floating
  metal contact (behind R or driven) whose metal meets a face reaction should start at its open-circuit
  level, where its face reactions pass the drive's current, as a port's start does: that needs a
  face-rate evaluator for the start (Stern and diffuse layers, MHC, surfaces). A bulk bath on the left
  could set the running φ̂ itself, making the start gauge-covariant without an electrode on the right.
  Newton's islands round accepts a state after two rounds with an island still 7.7e-5 out of balance,
  unwarned. A reaction balance among mobile species loses its conserved total on long steps as traps
  did (they had no per-node row; a mobile combination has none to fall back on, since its row would
  need the fluxes assembled apart from the rates): measured by Astra's ninth pass, a sealed film whose
  dopant light ionizes (Red = Ox⁺ + e⁻, neither crossing a terminal) gained 1.8e-8 of its total in a
  day and lost 12 ppm in a year on a warning-free grid, converged and passing check(); steps capped
  at 60 s cut it 20-fold, and in the dark it's zero. check() has no ledger for such a closed moiety,
  nor for a trap's per-node totals (it needs a reference carried through set()).

- A sub-grid diffuse layer's draw from everything that crosses it: an ion with a permeability link
  and a face reaction at once is drawn by the reactions alone, its channel at the undrawn plane
  (negligible at a membrane's fluxes, two-fold at 1.4 mol/(m²·s); Argon's review of 1.6.1). And a
  draw that respects a crowded layer's sites (now warned of where the plane overfills).

- **Steady solves the fuzzers still beat (round 15, each with a repro in
  `~/.cache/driftlet-claude/round15/fuzzB/work/`).** A reverse-biased diode behind a large R
  converges cold now (the load line), but its steady solve can stall and be taken as converged at
  round-off (Newton's updates near 1e-6 thermal units, not shrinking), a few µV off with Kirchhoff
  off 1e-6 to 1e-4 (behind 1e7 Ω·m² and 3.2e4 Ω·m²): 4000× and 200× below check()'s grid error
  there, and check()'s balance now fails the worse one. The round-off acceptance could ask the
  residual's true floor (a held solve gets Kirchhoff to 1e-20 A/m²; the scaled rows stop at 3e-14 of
  terms 1e10 times the current). Two terminals behind resistances with nothing held, past ~16 V
  across a silver cell; a resistance beside a current at 30 V (16 and 100 converge); a slow
  reaction's failure unnamed, and whether it converges, varying with the listing order (Y = Z, Y
  mobile in part); and the open-circuit march's sign hint, wrong beside a resistance.

- **From round 9's mirror property (2026-11-20; repros in
  `~/.cache/driftlet-claude/codex/runs/sol9/driftlet/out/repro/`).** Its three refusals are fixed (a
  factorisation from the right where the left loses a column; a carrier fed at both ends started
  from the 'bulk' contact). Left: the same GaAs | InGaAs junction in the dark, its Schottky on the
  left, is still refused cold, since the voltage continuation ramps the right contact, the ohmic one
  there, and a warm step of an ohmic contact at that bias stalls at Newton's round-off acceptance
  (updates rattling at 1e-6 to 8e-6 thermal units) in either orientation, where stepping the
  Schottky contact solves to 1.26 V (ramping the non-'bulk' contact fixes it, untried on the suites);
  the start's level for a carrier fed at two 'bulk' (or two non-'bulk') contacts is still the left's;
  and a device held only in its middle, closed at both ends, would need a twisted factorisation.
  Also seen: membrane 21 seeded at its exact steady state, its bath on the left, warns "a floating region doesn't balance"
  (imbalance ~1, its balance and the fluxes crossing it both round-off, ~1e-23 A/m²), as 1.19.0 did.
- **From round 17's fuzzer (2026-10-30; repros in `~/.cache/driftlet-claude/round17/fuzz/repro/`).** No
  silent wrong answer in ~1,840 devices; the cold start's missing face surface species is fixed. Left:
  a failed solve's warning names the reaction out of equilibrium at the failed iterate, not at the
  start, and its advice (solve without it) can be the wrong one; a metal on the right with c0 given
  starts ~55 RT from the bath (start.js's Known list: 20 of 273 mirrored face devices fail cold,
  converging with c0 dropped); a closed end reports ~1e-8 A/m² in transients, its end box's balance
  residual (steady states read 0); a port refuses chemical steps (no electron from the electrode).
- **The impedance read by Kirchhoff through held terminals too.** Where a terminal's current is
  lost to round-off but another held terminal's isn't, the measured one's is minus the sum of
  the others': a MOS whose back contact sees the silicon only through a channel port, its gate
  floating, reads Z at the back contact as noise (relError ~100, now warned of) while the
  channel's own reads to 5e-11. Today only terminals that all float are summed this way, from
  their voltages. (From the impedance module's cold review, 2026-10-20.)

- **Water's μ° at other temperatures.** `aqueous(list, { T })` moves D (Stokes–Einstein on
  water's viscosity) and ε_r (Malmberg–Maryott), 0–50 °C, since 5 Nov 2026. Still to come: μ°ᵢ(T)
  from the standard entropies (and heat capacities where tabulated; NBS/CODATA, conventional for
  single ions, which cancel in any neutral combination), checked against K_w(T) and an E°(T) from
  the tables. Also H⁺'s and OH⁻'s D(T), which proton hopping keeps off Stokes–Einstein (22% at
  0 °C), from measured limiting conductivities, and a check of everyone's D above 25 °C. Until
  then μ° are 25 °C's, warned of where a reaction reads them.

- **More interface kinetics.** (Marcus–Hush–Chidsey rates at faces and ports and surface species at
  faces done 2026-10-29.) Surface species carrying charge (in the face's
  φ law) and lateral interactions between them (Frumkin isotherms); and a custom forward rate
  $`r_f(\mathrm{state})`$ as an escape hatch, which the solver multiplies by $`(1 - e^{-A/RT})`$ so
  equilibrium stays exact. Not a bare current–overpotential curve, which loses the concentration
  dependence and can break detailed balance.


- **Cross-species transport coefficients together with cross chemical capacitances.** They're
  the two halves of one Onsager / Jamnik–Maier network, so one shouldn't come without the other.
  This is also the general mixed conductor with $`n > 2`$ mobile species, beyond the metal ($`n = 1`$)
  and insertion-host ($`n = 2`$) cases.
- **Concentration unknowns for dilute regions, as an option.** A hybrid basis: conductor regions
  keep their Fermi level, but a dilute region may ask to be solved internally in $`\{\phi, c_i\}`$
  instead of $`\{\hat\phi, \eta_i\}`$, where storage and neutrality are linear and a strongly diluted
  species (1 µM against 3 M) doesn't cost digits in its level. Transparent to the API: the solution
  still reports $`\bar\mu`$, $`V_i`$ and $`c`$, and faces still align levels.
- **More performance**, measured against the committed benchmark baseline: Newton-iteration
  counts (above), and specialised small-block elimination.
- **Concentration-dependent diffusivities** D(c), beyond what non-ideal statistics already give
  the flux through the chemical potential.
- **More statistics.** Gaussian and exponential densities of states (disordered and organic
  semiconductors), species occupying several lattice sites, several models on one species
  (e.g. crowding plus activity coefficients), and cross-model shared sites.
- **Graded materials.** Done (10 Nov 2026): a region's $`\mu^\circ(x)`$, $`c_{\mathrm{ref}}(x)`$,
  $`D(x)`$ and $`\varepsilon(x)`$ (`regions[r].grade`), with the full $`z\hat\phi + \mu^\circ/RT - \ln c_{\mathrm{ref}}`$ change
  across each segment, checked against the exact flux down a linear $`\mu^\circ`$ and against
  staircases of uniform layers, and $`D`$ and $`\varepsilon`$ against a series resistance and
  capacitance (test/graded.test.js). $`\phi`$'s gauge freedom is a continuous
  function there, with a built-in displacement $`D_{\mathrm{inhomog}}`$ moving to match (ESBD's
  inhomog page); driftlet takes $`D_{\mathrm{inhomog}}`$ only at faces, so within a region it's
  0 and the split of a grading between the charged species is an input. (26ff538 called this
  entry's "continuous gauge freedom" wrong; it wasn't: the model fixes the gauge, the physics
  doesn't.) A grading's $`\mu^\circ`$ given in another gauge, one with a built-in polarization
  $`P_0(x)`$ ($`D_{\mathrm{inhomog}}`$) inside the region, needs no knob of its own: Poisson sees
  only $`\rho_{\mathrm{fixed}} - \nabla\cdot P_0`$, so those $`\mu^\circ`$ as given, with
  $`-dP_0/dx`$ as fixed charge and each jump of $`P_0`$ as a face's `sheetCharge`, is the same
  device (every $`\bar\mu`$, $`\phi`$ and current; only $`D`$'s split of free and bound differs).
  Done 12 Nov (planned with the owner on the 11th): `fixedCharge` as a profile (`{ x, values }`), which also
  covers graded doping; held to that equivalence (the plain diode against its μ° moved by
  $`z_iFf`$ with $`-\varepsilon f''`$ as fixed charge, second order in the cell) and to a
  staircase of uniformly doped layers. Still to come: a graded $`c_{\mathrm{ref}}`$ under the crowded and activity
  models (lattice, Debye–Hückel, insertion); a non-ideal $`\mu^\circ`$ and a Fermi–Dirac
  $`c_{\mathrm{ref}}`$ grade since 13 Nov. Done with it: a bowed alloy (`alloy(…, { bow })`, 13 Nov), the kit's
  `alloy()` (between two library materials, its alignment given as at a face) and a graded-base
  transistor demo (demos/hbt.html, its transit time held to Kroemer's formula).
- **More of the diffuse layer below the grid** (`interfaces[f].diffuse`, src/diffuse.js, done
  2026-10-25 for electrodes, charged faces, membranes' permeabilities and gates, crowding on the
  bulk's lattice, and a reactant's transport across the layer). Not yet: a contact's own (a
  capacitive contact holds its charge in its end cell, warned of); other statistics in the layer
  (Debye–Hückel, Fermi–Dirac); SRH carriers read at the plane; a permeability link's ions'
  transport across it; and the extended space charge of a salt run out at the electrode (past an
  unsupported electrolyte's limiting current), which no thin layer holds.
- **Driven species fluxes.** A terminal driving one species' molar flux instead of the charge
  current (a gas feed at a fixed rate, a neutral species injected), the flux-side twin of a held
  $`\mu`$ for neutral species.
- **Interface states** (charge that depends on the local $`\bar\mu_{\mathrm{e}^-}`$, i.e. Fermi-level
  pinning) and thermionic-emission links at heterojunctions.
- Non-isothermal transport (Soret, Seebeck/Peltier, a heat equation), which stays local.

- **Not planned, though not ruled out: machinery heavier than a library for live demos usually
  needs.** Worth taking on if a demo, or the physics, really calls for it.
  - **Solid phases with complementarity**: a phase at unit activity while any is present, its
    amount never negative (n ≥ 0, saturation index ≤ 0, one of them zero), so it appears,
    dissolves away and reappears: a pit's salt film, a real precipitate (rust rather than an ideal
    solute). It needs semismooth Newton and a stepper that finds when a phase appears or
    vanishes, which is where reactive-transport codes spend their effort. Surface coverage
    (Langmuir, smooth, bounded by its own statistics) stands in for a passive film.
  - **Phase separation (Cahn–Hilliard)**: a non-convex free energy (LFP's miscibility gap) with a
    gradient-energy term, which is still local and block-tridiagonal. A regular-solution OCV
    through the gap is the stand-in.

## Out of scope, for good

2D/3D; non-local couplings such as recirculation or one well-mixed reservoir feeding several
faces (a region with large mixing is the local stand-in); anything that breaks the
block-tridiagonal structure.
