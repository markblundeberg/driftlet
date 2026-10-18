# Changelog

driftlet follows semantic versioning: from 1.0, a change that breaks working code comes only with a
major version.

## 1.2.2 (2026-10-18)

The steady system's conservation laws now live in `src/laws.js`, found from the device's structure
alone, with a contract at its head (what it reads, what it returns, what must hold;
[internals](docs/internals.md)). A cold review of that module against its contract found five
more laws kept wrongly, each now fixed and in `test/laws.test.js`:
- An immobile species linked across a face was kept node by node, though the face carries it
  between its two nodes (solve() left them at 1 and 0.1 for the transient's 0.4 and 0.4).
- Two node-by-node laws could take the same row, and one was lost (an immobile S⁺ with
  S⁺ + A(ads) + e⁻ = B(ads): θ_A settled at 1 for 0.3). Their rows are chosen by elimination,
  node by node.
- An immobile ion that a port at no current links was kept both node by node and in the port's
  circuit row, and the port's I = 0 was lost (0 V for −1.42 V).
- Two laws through what a capacitive contact at no current passes: the second left out the
  contact's charge, and the solve failed. The charge is kept once, in a circuit row, as an
  electrode's is.
- A capacitive contact behind a resistance: time steps reading its current by the change of
  what it passes left out its charging current, and failed. It's read by the flux.

A law that no row of the steady system can keep is set aside, and the steady solve goes by time
steps instead, where before it was dropped.

How the laws enter the Newton system (each evaluated on the state and written in place of the row
it takes, and the bordered solve that holds the dense ones aside) moves from the solver to
`src/lawRows.js`, with a contract of its own; every solve in the test suite comes out bit for bit
as before. A cold review of it against that contract found one slip: a node-by-node law's row, in
the window of an electrode at no current, kept the port reaction's ∂/∂V from the row it replaced.
The answer was right (the residual was), but Newton took an extra iteration to it.

## 1.2.1 (2026-10-18)

Three more conservation laws a steady solve kept wrongly, found the same day (two by Gemini 4
Argon's third round, one by a new sweep), each converging where a long `advance()` doesn't go:
- Two contacts at no current passing an ion, one a capacitor: the law kept only the charge of
  the contact whose circuit row it took. It keeps every feeding contact's now.
- A charged immobile species: at the nodes where nothing touches it, its steady balance is
  0 = 0, and its level stayed put while φ moved, so its concentration didn't (S⁺ fell to
  e^{−0.2/V_T} of itself with the bath moved 0.2 V, and an electrode settled 19 mV off). Those
  nodes now hold their concentrations.
- An immobile species adsorbing on an electrode spread over part of its region: its local law
  (S + Γaθ where the electrode is, S alone elsewhere) needed the two to cover the same nodes,
  and was dropped.

`npm run lawsweep` (tools/lawsweep.js) builds 810 small devices from these edges and compares
each steady solve with a settled transient; see docs/reliability.md. The law code names its
column spaces rather than doing arithmetic on raw indices.

## 1.2.0 (2026-10-18)

Fixes from a third outside review by Gemini 4 Argon: a breadth-first pass that flagged 35
things, each checked here, then a pointed hunt for steady states that drop a conservation law.

Wrong answers, fixed:
- **Steady states that kept the wrong amounts.** A steady solve replaces one redundant row per
  conserved combination with the combination's amount; seven combinations at the edges of that
  bookkeeping were dropped or weighed wrongly, each converging where a long `advance()` doesn't
  go: a floating electrode whose adsorbates gave it electrons (−53 mV for −37.8); two floating
  electrodes sharing a closed ion (1.3 V for −0.65); adsorbates reacting with a closed ion (the
  sites' total 0.125 for 0.3); a capacitive port or contact at no current passing an ion (the ion
  kept without the capacitor's charge, Q/F); an immobile species that adsorbs (a coverage added to
  a concentration without the sites per volume, Γa); and a face reaction making an immobile
  species beside a floating electrode (its rate said twice, the electrode's I = 0 lost). Each
  device is now a test against the transient.
- A face reaction's affinity at a contact took +Σνη where every other one takes −Σνη, so a
  reaction with a fixed participant could be held at equilibrium when it isn't: a face ran 3×
  too fast.
- `henderson()` returned 0 when the two solutions conduct equally; it's the limit, −V_T Σ z D Δc
  / g (8.56 mV in the example, beside Planck's 8.66). An empty solution is refused by both
  formulas.
- `trap()`'s midgap n₁ used √(N_c N_v e^{−(μ°_h − μ°_e)/RT}), not the sum: invisible with the
  library's materials (μ°_e = 0), 1.1% in a diode's current with φ's reference moved 1 V.
- `live()`: a change whose ramp passes through an invalid state (a bath neutral at both ends but
  not between) failed without trying the fresh solve that works.
- `ghkCurrent()` was NaN past u = −710; it's finite everywhere.
- The nerve demos dropped each call's first step (the trace is that call's own steps, all of
  them): myelin.html kept one in eight. And propagation.html paired the wrong spikes when two were in flight (17.7 m/s for 3.90).

Additions: grids take `hLeft` and `hRight`, as the types said; a solution says its `geometry`;
`themeColor()` takes φ's role. The types (`types.js`) now cover gates, permeability, SRH and
saturating face kinetics, baths on contacts and ports, and the solution's flux, D, T and rates.

Sharper errors: `combine()` given trap()'s whole result, a `bath()` option typo, a species on
both sides of a face equation, an unbalanced half-reaction's level, `hodgkinHuxley()`'s `g` and
`names` typos, `perovskiteCell({ light: -1 })`, `hysteresis()` with a turn outside the trace,
`themeColor()` with no such role, `speciesRole()` of a neutral or absent species, and
`particles()` off a plane (its fluxes per device, its amounts per area). A Debye–Hückel
composition that doesn't converge says so, and a partly given one inverts by Newton.

Smaller: `check()` gave two grid items with `order: false`; `describe()` printed α undefined;
a negative-time breakpoint stalled; the demos' readouts follow a slider while it's dragged, a
`?t=` link plays live after a restart, the particle loops rest offscreen, and cell.html and
axon.html take water's 25 °C data as their own temperature's, having checked that only D counts
there and moves nothing shown.

## 1.1.1 (2026-10-16)

Fixes from a second outside review by Gemini 4 Argon (in Antigravity, running its own scripts),
each finding rerun and checked before it was taken.

- The impedance with non-ideal statistics (Fermi–Dirac, lattice, Debye–Hückel, regular
  solutions) lost its small currents to round-off: those regions' storage and fluxes were
  assembled directly, where ideal ones are kept in difference form for an exact J·v, and their
  segments weren't cuts where the current could be read. A Fermi–Dirac MOS capacitor in
  inversion, its surface only 0.03 N_C, had Re Z = −10.8 Ω·m² against its Boltzmann twin's
  +2870, and its low-frequency step two decades off; it matches to 0.2% now. A blocking
  electrode in a crowded electrolyte had a negative real part at 1 mHz; it's positive, and a
  lattice gas with a vast site density matches ideal ions to six digits in Im Z.
- `aqueous()`'s water says its data are 25 °C's (`T: 298.15`), so a device at another
  temperature warns, as one with 300 K silicon does: its μ° (so K_w, ~6× off at 6.3 °C, and
  electrode potentials), permittivity and diffusivities don't hold there. Scaling them is on the
  roadmap.
- SRH with Fermi–Dirac carriers split their degeneracy evenly between n₁ and p₁ (or, with a
  given n1, put it all on p₁): a degenerate n⁺ region's γ_n reduced its holes' emission too.
  Detailed balance scales each carrier's by its own γ = 𝓕₁/₂(ζ)/e^ζ, and now so does the rate;
  the product n₁p₁ = n p e^{−a} was always right. It hardly changes a device (where a carrier is
  degenerate, it outnumbers n₁ by far). The rate also stays finite for a quasi-Fermi split past
  ~709 kT (e^{−a} overflowed), though no device we tried gets near one: a diode reverse-biased
  to 25 V keeps its split under 30 kT, generation refilling it.
- Conserved combinations came from eliminating the stretches in listing order, so a trace
  species could come out added to an abundant one: A + B ⇌ C with B at 1 M and A and C at
  1e-14 mol/m³, listed B first, gave B + C and A − B, where A + C was lost to round-off, and the
  solve didn't converge or converged with A + C 65% off. They're eliminated by starting amount,
  least first, then by name, and give the same combinations in any listing.
- The extra precision kept for each level (its low part, `uLo`) never reached a contact's row
  or a flat stretch's: they computed u + uLo − target, and u + uLo rounds to u. They subtract
  first now. An electrode at open circuit with a fast couple and a slow one (round 11's F5),
  whose steady states were 4 mV apart by the species' listing, now gives the same in both,
  its ledger closing to 4e-16. check()'s order test is now tested on a solve made to differ.
- A double layer across more than ~360 thermal voltages (9.25 V at 25 °C, in a crowded
  lattice-gas electrolyte) didn't solve, and its warning didn't say why: its co-ion, crowded out
  to ~1e-310 mol/m³, made a row whose largest entry was subnormal, and scaling the row by that
  entry's reciprocal overflowed to Infinity. Rows are divided by it now. The same layer solves to
  12 V and beyond, past where the co-ion underflows to zero.

## 1.1.0 (2026-10-16)

- `semiconductor(name, { T })` takes a temperature, 150 to 500 K, where it was refused: the gap
  moves by Varshni's shift from 300 K, N_c and N_v as T^1.5, and the mobilities by
  lattice-limited power laws (light doping), with D = μkT/q at T. At 300 K, or with no T, it
  returns the table exactly, as before. Checked against a silicon cell's open-circuit voltage,
  which falls 2.14 mV/K at 300 K, within Green's dV_oc/dT = −(E_g0/q − V_oc + γkT/q)/T for its
  transport's γ (it implies γ = 1.67, against 1.66 from its two parts of J_0).
- The kit's `trap()` writes a trap level that holds charge: its two states as immobile species,
  filled and emptied by capture, from the `srh` law's `{ tauN, tauP, n1 }` and a density. In
  steady state it recombines at the same rate, and its charge is in Poisson's equation, which
  `srh` leaves out (checked against n = N_D − N_t n/(n + n₁) for acceptor-like traps
  compensating donors, and the donor-like mirror). Explicit traps always worked; nothing in
  llms.txt said how.
- A material may carry the temperature its data are for, `T` (K), and `semiconductor()` always
  sets it (300 K by default). A device at another temperature now warns, in its solutions'
  warnings (so `check()`) and in `describe()`, naming both temperatures and the cure. Before,
  a device at 350 K on 300 K silicon had an intrinsic density 36% low with no word, and one left
  at the default 298.15 K, 13% low; `set({ T })` changed kT but not the data. To change T on a
  live device, `set({ T, materials: semiconductor('Si', { T }).materials })`. `trap()` takes its
  material's T. Hand-written materials without `T` are taken as written, as before.
- `describe()` lists the definition's own warnings too (the temperature one, an unresolved
  starting profile, a statistics model's conventions), as every solution already did.
- Error messages show the value given: a T of NaN or Infinity said "got null", and a string
  where a number belongs is named as one (`got "350" (a string)`). `build()` checks its `T`.
  A conductor material now refuses fields it doesn't take, as other materials do.
- A failed steady solve's warnings say why it failed, and leave out those judging a solution's
  accuracy (the grid, double layers, space charge at a contact): the state it stopped in isn't a
  solution, and a diode driven at a current of the wrong sign, stopped at 26 V, came with two
  spurious "space charge reaches the contact" warnings.
- `check()`: two bulk reactions with the same equation (radiative and SRH recombination, both
  `e- + h+ = 0`) were merged into one ledger term; each now has its own, named by its place.
- Docs: where N_c and N_v live in a material, χ in eV throughout, the device's 298.15 K default
  against the table's 300 K, a diode's J₀ ∝ T^(4−γ) e^(−E_g/kT) and its Arrhenius slope, junctions
  that fade when n_i passes the doping, and a recipe for an alloy (AlGaAs from GaAs).
- llms.txt: the junction template is a short diode (its Shockley check takes the quasi-neutral
  widths, not L), and the plot module draws profiles in x, not I–V curves.

## 1.0.1 (2026-10-13)

Fixes from an outside review of the physics and numerics by another model family (Gemini 4
Argon), each finding checked by a script before it was taken. Nothing in the API changes; the
surface-kinetics fix changes the rates of reactions with surface species that don't run only on
bare metal, which were wrong.

- A charge-free dielectric between strictly neutral (ε = 0) electrolytes, a film between two
  baths, had φ undefined: no solve, and φ NaN in the film. An ε = 0 region with ions now fixes φ
  across a face that isn't neutral, as a conductor does, and the film charges as a capacitor.
- A contact exchanging only a neutral species counted as passing current, so it could be the
  voltage another terminal's current drive was read against, and the solve then failed without
  saying why. It's refused at construction now.
- A gate contact at open circuit (`I: 0`, capacitive) started at 0 V rather than at its
  zero-charge level, so `solve()` kept the charge of that difference: C·1 V on a gate beside p-Si,
  the silicon inverted at open circuit. It starts uncharged now, as a floating port's capacitance
  did.
- Surface species' kinetics are Langmuir throughout. The Butler–Volmer prefactor took each at its
  activity θ/θ₀, so a reaction that needs no free site sped up as the surface filled: an adsorbed
  couple `'Oads + e- = Rads + Cl-'` converted 76× faster for ninefold coverages, where Langmuir
  gives 9×. Each now counts at its coverage θ, and a reaction's transition state holds the bare
  sites it fills (θ₀^m forward, θ₀^(m−Δ) back). `bare: true` keeps its meaning, a reaction that
  needs bare metal, and the passive film behaves as before; on an adsorption reaction it's no
  longer needed. A reaction with surface species and no `bare` runs at a different rate than in
  1.0.0. Equilibria are unchanged.
- A port holding an ion over a window that ends at a contact blocking it: the contact put the
  end node's deferred ion storage into the port's held-level row, pulling the level off (16 mV
  for 50 on a 1 µs step) and reporting a flux through the blocked contact. That ion's row is the
  port's now, storage and all.
- `impedance(f, { profiles: true })` gave numbers for concentrations inside a metal, where a
  solution has NaN.
- docs: a velocity link's flux is v c (e^(Δμ̄/RT) − 1), which is v (c_eq − c) only for dilute
  species.

## 1.0.0 (2026-10-13)

The API is stable from here. 1.0 was aimed at its likely users, coding agents asked to simulate
something or build a demo: getting a correct result on the first try, and noticing when they
haven't. Two rounds of agents built demos cold from llms.txt alone (a pn diode, a battery, an
impedance page, a solar cell, a cyclic voltammogram), each checked against theory (to 1% or
better, the impedance to 0.4–5%, where circuit theory leaves out migration), and what they tripped
on is fixed below.

- Every option name is checked. `solve()`, `step()`, `advance()`, `impedance()` and every kit and
  plot helper refuse an option they don't take and suggest the one meant (`advance(t, { dtmax })`
  ran at the default step size; `bandDiagram(sol, { engery: true })` drew the other picture). A
  misspelt contact in `set()` after a solve was accepted, ignored and kept; it's refused now.
- `check()`'s balance called a diode passing 0.1 A/m² "equilibrium": flows below 1e-8 of D·c/L of
  the majority carriers were taken for round-off. Flows that balance, as real ones do and
  round-off doesn't, now count down to 1e-20 of it.
- `hodgkinHuxley()` requires `T`: its rates scale with it, so it must be the device's (their
  6.3 °C under a body-temperature device doubled a spike's peak, silently).
  `semiconductor(name, { T })` is refused (the data are for 300 K).
- `impedance()` returns the steady state's `warnings` too.
- A face reaction beside a strictly neutral edge cell that holds the face's charge is warned of
  past a 1% shift: a Frumkin effect set by the mesh, which refining makes worse (an electrode's
  R_dc fell 16% from h = 300 nm to 1 nm, silently until ~3 nm).
- A floating port wired to an insertion host's carrier starts at the host's OCV, so a reading
  before the first step (or after a failed solve) means something; it read 0 V.
- An OCV curve too flat for its interpolation (a plateau) says what passes instead of "add
  points".
- llms.txt: a battery-particle template (a sphere discharged at two rates, against the quasi-steady
  profile), and notes found missing by agents building demos cold; repeats cut. Its six
  specialised templates (Schottky, thin-film transistor, corrosion, voltage clamp, action
  potential, axon) move to docs/templates.md, so llms.txt (~20k tokens) fits an agent's
  single read.
- Toward 1.0, the exports are only what's documented. `driftlet` no longer exports
  `BlockTridiagonal`, `bernoulli`, `bernoulliDerivative`, `buildGrid`, `gradedCells` or
  `normalizeDevice`, and `driftlet/kit` no longer exports `parseEquation` or `stoichiometry`.
- A solution's `gates.left/right` (the charge on a gate or Stern plate at a capacitive contact)
  is now `plates.left/right`: `gates` are a membrane's channel gates everywhere else.
- `advance(t)` read as a duration says so: called with a time already passed (or reached, outside
  an animation's `budgetMs`), it warns that t is absolute, rather than silently doing nothing.
- A grid of more than a million cells is refused at once, with the region and its length. A
  length in µm given as metres (`layer('Si', 1)`) used to build a grid of 1e8 cells, silently,
  until memory ran out.
- `build()` says what to write when a contact holds a carrier its layer lacks: `ohmic()` on a
  metal, whose default holds h⁺ too, now asks for `ohmic(V, ['e-'])`.

## 0.16.4 (2026-10-12)

- Steady solves depend less on the order species and reactions are listed in. A pivot that was
  round-off alone without being exactly zero (strictly neutral water whose level only the fast
  couple at its far face held: 7e-32 in one listing, 3e-16 in the other) scaled all that followed
  by 1e31, and the solve failed in that listing only; such pivots are now perturbed and refined
  like exact zeros. Which balance row a reaction's equilibrium stands in for was picked by where
  its species were listed, and an electrode's two couples solved in one listing and failed in
  another; the rows are now the same in any listing, the slowest reaction first. Of 600 random
  reacting cells, those that converged in only some listings or to different states fell from
  11 to 8. Devices with electrode ports remain at ~10%, much of it a second steady state (see the
  roadmap).
- `check()` solves the device again with its species listed in reverse, and fails where the two
  differ: one of them isn't steady, though its solve converged.
- `check()`'s balance fails where the solve left a reaction stuck (one sign of affinity wherever
  it acts) rather than call it equilibrium.

## 0.16.3 (2026-10-12)

Mostly found by round 11 of agents trying what's new cold (batteries at rest, slow couples and
mixed potentials, floating and resistive gates, an invariant fuzzer), whose results otherwise
matched theory: relaxed GITT voltages to 1e-9 V, Wagner–Traud mixed potentials to 1e-5 mV.

- A battery at rest steps out to any time. Each step's first guess extrapolated an insertion
  host's φ̂ too, a gauge that nothing pulls back, and through steps that double it ran off
  geometrically (1e11 thermal units by 1e5 s), until the levels measured from it lost their digits
  and every step failed (from 1e5 s at tol 1e-6, 400 s at 1e-8). It now reaches 1e9 s at any
  tolerance, the voltage at its hosts' OCV difference to 1e-12 V.
- A face reaction that alone changes what nothing else feeds (Ag⁺ in a closed solution, beside a
  faster couple on the same electrode) is held at equilibrium in steady solves by its affinity,
  not its rate, so it's right at any k0. At k0 1e-20 the rate was lost to diffusion's round-off:
  the solve settled where silver wasn't (−0.4661 V for −0.5004), silently, and 1e-15 to 1e-18
  failed. That covers an electrode at open circuit feeding a closed pool too (78 mV off before),
  and a couple whose species also complexes.
- A bulk or electrode-port reaction a steady state needs at equilibrium, too slow to see, is
  warned of rather than passed off as steady (a slow complexation had kept its starting
  composition, K = 0.1 for 7.5).
- A steady solve converges only where its residual is small, not just its update. A pivot that
  cancelled exactly, or a factorisation that lost a reaction's row, left a rate row unsolved
  under a small update: one couple written twice on an electrode at open circuit sat at 0.7293 V
  for 0.6906, silently, and in one species order only.
- `advance()` integrates short calls late in a run. It skipped any remainder under 1e-10 of the
  absolute time as round-off: a 0.01 s GITT pulse at 1e9 s was skipped whole (the time moved,
  nothing else did), and from 1e10 s a pulse train lost 17% of its charge.
- A floating terminal reports what its circuit passes, exactly: a gate behind R read its current
  as the flux into its metal, 1.5e-4 off wherever Newton stopped an update early.
- `set({ contacts: { left: { I: 0 } } })` on a contact `{ V, R }` drops the resistance too (it
  threw).
- A time step that fails says what usually does it, most likely first: a host or surface at its
  limit; an electrode at open circuit with no double layer, held only by a slow reaction (give
  its face a capacitance).

## 0.16.2 (2026-10-09)

- A terminal driven by a current that alone feeds what it keeps (a floating gate's metal, a
  host's electrons behind an open circuit) reads its current in a time step as that amount's
  change, where that's the quieter reading. A floating gate at rest took 5000 steps per 100 s,
  half rejected, its charging current lost below its metal's round-off (it now takes the 44 a
  held gate does); a closed host at open circuit drifted from its charge (to 0.3995 V from
  0.4 V by 1e10 s), and its steps failed by 1e14 s. A terminal behind a resistance reads its
  current the same way: a gate charged through one was 0.9% off its RC response, silently, and
  took 3280 steps for 5τ, half rejected; it now follows backward Euler's exactly, in 60.
- A steady solve from a transient's state that Newton can't take goes on as a cold device would,
  from the layout for its drives: a lit MOS capacitor a nanosecond after its gate jumped failed
  where a cold solve takes 27 iterations; it now takes 32. A continuation no longer starts from
  the last steady solve's voltages once a transient has moved the state.
- `set({ grid })` carries the state over, interpolated onto the new nodes, where the regions
  and terminals are the same: a floating gate kept nothing of its charge through a regrid.
- A solution warns of a species denser than any material (1e6 mol/m³): a contact holding one
  ion at a level far from its neighbours' held it at 1e37 mol/m³, its currents then reading
  1e28 A/m², with no warning.
- A species that can't move on either side of a face (D = 0: trap states) is blocked there by
  default. Crossing freely, its two edge boxes were held at one level, whatever flux that took:
  two regions' traps mixed, and steady solves failed.
- `check()` leaves a contact that blocks a species out of that species' ledger: the flux it read
  there was the end box's own imbalance, and counted as a source it had balanced the books of a
  state that wasn't steady.
- The stress test solves each case again with its species listed in reverse, and has a family of
  reacting solutions in several regions, with an electrode at no current: before 0.16.1's fixes
  it failed 237 of 1496 checks.

## 0.16.1 (2026-10-08)

Found by a round of agents trying the new steady solves cold (GITT on an insertion host,
floating gates, closed thin-layer cells), whose results otherwise matched theory, a closed
battery's rested voltage to 1e-11 V:

- A steady solve of a device driven by a waveform read it at its pseudo-transient steps' clock,
  where that ramp ran: a triangle wave's steady state came out at its value seconds later (1 V
  instead of 0), silently. It's read at the present time now.
- A level held flat (a species reaching one contact alone) is flat exactly, every node at its
  first's. Flat only to round-off after a jump, a MOS capacitor's bulk holes, behind a floated
  gate with the back contact jumped to −6 V, had read as a current 1000 times the true leakage,
  silently: the terminal currents didn't add up.
- A film whose ions are all given by `c0` starts its electrostatic potential balanced against
  its electrode's surface species too. An adsorbed couple had started 1.4 V off: cold solves
  failed or took the pseudo-transient ramp, and `advance()` from the start stalled.
- An insertion host's occupancy is inverted exactly on a steep tabulated OCV (a 60 mV step over
  Δx ≈ 0.1): the root solver zigzagged across the step's inflections, ran out of iterations at
  points that weren't roots, and time steps failed down to 1e-14 s.

A second agent compared steady solves with long transients across random devices, and found two
more, both silent, both now matching to round-off in every species order:

- Conserved amounts overlapping (two complexations sharing a species, A + B = AB and
  A + C = AC) could stand in for rows that weren't independent, one balance said twice and
  another not at all, depending on the order the species were listed in: c_A 7% off, a phantom
  flux through a contact that blocks it, an I = 0 terminal passing 20 A/m². It hit 8–13% of
  random reacting multi-region devices. The rows are chosen by elimination now.
- An electrode port at no current whose reactions use up a closed ion (M⁺ + e⁻ = M(s)) keeps
  that ion together with the charge on its capacitance (or alone, with none). Unseen, the steady
  solve had no law for it and landed anywhere: c(M⁺) 2.2e8 mol/m³ where it's 10.

## 0.16.0 (2026-10-08)

The steady solve, reworked around one idea: every row it changes from a time step's is a sum of
balance rows, written out exactly (a conserved amount, a level held flat, an island's balance),
with its conservation laws now worked out under the terminals' drives. [Numerics](docs/numerics.md#steady-state)
describes it that way now.

Silent wrong answers, fixed:

- A closed battery at open circuit keeps its state of charge: a cold `solve()` had converged to
  another of its steady states (4.145 V where its `c0` said 3.958 V). What a terminal driven at no
  current alone feeds is conserved, at what the state holds when the solve begins, so after a
  transient the steady state is where that transient relaxes. An insertion host starts with its
  ion level at the electrolyte's (its φ is a gauge), so the floating terminal of a closed cell
  starts at its OCV, and transients from `c0` run at any current (the agent guide's workaround is
  gone).
- Adsorbates on an electrode that only turn into each other keep their total on every site; such
  a device had failed to build.
- A MOS capacitor's back contact moved under its biased gate: where its carriers recombine (so
  weren't held flat), its inversion layer couldn't follow through a scarce-minority bulk, and of
  500 random MOS capacitors 38 failed and 5 converged to a gate charge up to 2× off. Species that
  react only among themselves, each reaching the same one contact and nothing else, are now held
  flat at its levels: in a steady state they're at equilibrium. The stress test moves the back
  contact now too.

New:

- A gate or a capacitance left floating (`I: 0`) keeps its charge in a steady solve, where it
  threw. Driven at a current, a terminal that only stores what comes in (a gate, a host filling)
  makes `solve()` fail at once with a warning saying so, rather than throwing or marching a held
  voltage 20 V each way.
- Steady solves go direct (rather than through huge time steps) for a floating metal (its amount
  the charge on its faces), and for trap states holding a closed population's electrons.
- A blocked ion's level is held flat across faces too: a closed zinc cell past its limiting
  current, its electrolyte in two regions, failed at 0.9 V and now solves as one region does.

## 0.15.1 (2026-10-08)

From a further round of cold agents (an ion-selective electrode, a diode's ideality, Hebb–Wagner
polarisation, a Gärtner photoanode), all matching their theory:

- The new warning of a space charge reaching a contact no longer fires on a quasi-neutral layer
  split off a junction's (it read round-off as a cut-off layer).
- `check()` compares a current too small for its ledgers to call moving (a diode at low bias,
  1e-4 A/m², changed 3.6% on a finer grid while `check()` said ok).
- Electrons or holes above their band's density of states under Boltzmann statistics are warned
  of (degenerate; a photoanode driven past its couple's diffusion limit converged that way).
- `describe()` names a bath without a reference species as read through its φ, and flags an SRH
  n₁ given per cm³.
- Docs: an exchange or partition constant as a μ° difference; `dt0` for a transient's first
  instants; a contact holding several species fixes its end's composition; a bath's composition
  `set()` mid-run.

## 0.15.0 (2026-10-08)

Mostly from four more rounds of cold agents (fresh tasks, from the agent guide alone): sixteen
of them, every device matching its theory once set up, and the places where it silently didn't.

Silent wrong answers, fixed:

- `photogeneration()` generates nothing outside `[from, to]`: light given to one layer of a
  material that's also elsewhere (a p⁺-i-n⁺ diode of one silicon) generated in the layers beside
  it too, 10–25% too much photocurrent. And an αL of exactly 50 no longer fails.
- A probe at a face reads the side where its species is (it read an electrode's metal side,
  NaN), and asks for `region` where both sides have it.
- `impedance()` at a held terminal with every other one driven by a current (no way back for the
  signal) is an error saying to measure at a driven one; it returned round-off.
- A port's window at a device end node holds its species there unless that contact links them:
  a wire on a particle's surface node alone now charges it (it passed nothing).
- `set()` of a parameter mid-run reads held waveforms at the run's time at once (they read t = 0
  until the next step).

New:

- A bath can name its own material (`bath: { c, material }`, or the kit's `bath(c, ref, drive,
  { material })`): a solution straight against an ion exchanger or a gel, whose levels the end
  node takes through a Donnan step (tested against Donnan exclusion and the Teorell–Meyer–Sievers
  salt flux, to 1e-8). Without it, a bath is the end material at that composition, fixed charge
  and all, and the error now says so and what to do.
- `set()` of a parameter mid-run (a rate constant, a diffusivity) goes on at the step size an
  adaptive transient had reached, rather than restarting from a tiny one: 7× fewer steps for a
  page changing k0 every frame.
- `particles()` stays a sample at any step: a cell a dot would hop out of many times in one step
  is a sea for that step and then drawn afresh, so a long step no longer over-fills cells (a
  checkerboard) or costs in proportion. It draws as many cells as asked, as nearly as the nodes
  allow (it drew 52 of 80 on a 1 µm grid); `cap` is a density; dots sit within a cell by its
  concentration; a solution on another grid redraws the lattice; bad options are errors, and a
  NaN or negative step (a first frame) moves nothing.
- The membrane page draws its ions as dots at its foot: chloride thinning across the membrane,
  Donnan-excluded, against seas of Na⁺.

Warnings and errors that say what happened:

- A space charge that reaches a contact holding its end neutral (a double layer or depletion
  longer than its region, cut off there: σ 58% low, silently) is warned of.
- Beside a strictly neutral material, a capacitive face's charge sits in the cell next to it,
  a diffuse layer that wide in series with C, so a fine grid lowered the double-layer
  capacitance silently (to 0.036 F/m² from 0.2 with 0.1 nm cells): warned of past 2%.
- `describe()` says where a bath's SHE level sits against its terminal when the bath is read
  through a reference ion (1.37 V below a 0.5 M Cl⁻ terminal on table μ°); warns of a bath out
  of equilibrium with a bulk reaction at its contact (H⁺ and OH⁻ a little off K_w, a current at
  zero bias); gives an ion exchanger's fixed charge in mol/m³; and flags a large k0 as an
  exchange current only at an electrode.
- A steady solve that fails where an ion is too scarce for double precision (a minority swept
  out of a junction) says so, instead of suggesting a floating region.
- `advance(NaN)` does nothing and says so in its warnings; `advance()` of a non-number is an
  error. A capacitance-only port given a terminal species that no material holds is told it
  needs none.

Faster:

- Block elimination skips the coupling blocks' zeros (bit-identical results): the nerve demos'
  factorisations are ~15% faster. The benchmarks gain a squid-axon spike.

Docs:

- The agent guide gains ion-exchange membranes, excitons in an organic cell, and how to start a
  closed battery (hold it at its OCV, then drive it).
- How fine a grid must be where nothing warns: a transient's diffusion layer, a spike's front, an
  impedance's reach into a porous electrode, a layer depleted past a limiting current.
- A MOS capacitor without a port solves, and its impedance shows the high-frequency C–V above the
  minority supply's knee; the library's K_w; an OCV curve far flatter than ideal at an end; the
  double layer's voltage limit (~350 V_T); conservation's drift is relative.

## 0.14.0

- `particles()` in the kit (experimental): a solution drawn as dots of a fixed amount each,
  hopping between display cells at one-way rates taken from the solution, so their density is
  its concentration and their net crossing rate its flux (tested against Fick's slab, a decaying
  species and Ussing's flux ratio). Dense species are seas that hold no dots but exchange them,
  so a junction shows its majority carriers as seas and its minority carriers one by one. The
  pn demo draws its carriers this way, at the bottom, on the band diagram, and the axon demo
  the ions crossing its membrane each way, at their unidirectional fluxes.
- Solutions carry each species' flux across every segment (`flux`), its diffusivity (`D`), the
  temperature (`T`), how each species crosses each contact and face (`links`), a permeable
  face's one-way fluxes (`oneWay`), and each bulk reaction's one-way forward rate (`forward`)
  and stoichiometry (`nu`).
- A visual identity: the logo's wordmark redrawn in Recursive's casual sans (as outlines, so it
  looks the same everywhere), and the demo site in paper and ink with the logo's vermilion, a
  gallery whose sections show each field's cast, charges drifting along a potential at its foot
  (cations down it, anions up), and a stamp on each card naming what its page is checked against.
- A starting `c0` profile steeper than the grid's cells under it is warned of in every
  solution, with where and the cells that would resolve it: the start was silently the profile
  sampled at the nodes (a liquid junction's potential then seemed to grow as it formed).
- `check()` in a transient says its grid isn't checked (`?`), rather than leaving the grid item
  out: a transient's answer depends on the whole run, so it's rerun on a finer grid by hand.
- A time step's Newton iteration stops once its contraction shows the next update would be under
  tolerance, rather than taking that update to show it: about a third fewer factorisations where
  steps are capped by `dtMax` (the nerve demos: a bare axon's step 4.4 → 3.3 ms), 4–8% fewer in
  the benchmarks' adaptive transients.
- The agent guide is shorter (a map of the templates up top, tighter front matter, the
  level-diagram and mistakes sections halved), with notes from a third round of cold agents: a
  BJT's base contact, the semiconductor data, `vacuumDipole` with `build()`, a cyclic
  voltammogram from the stirred-electrode template.

## 0.13.0 (2026-10-07)

- A port can take its outside's composition, `bath: { c }`: its V is then the outside's φ, and
  each linked species' level there follows from the composition (an offset no longer has to be
  worked out). `hodgkinHuxley({ T, area })` in the kit writes a membrane port's links, Hodgkin
  and Huxley's linear conductances gated as theirs.
- The action-potential demo can put the membrane behind Frankenhaeuser and Hodgkin's 30 nm
  periaxonal space: K⁺ piles up there with each spike (about 3.5 mM) and clears over 50 ms, and
  E_K, the rest and the after-hyperpolarisation follow it, which Hodgkin and Huxley's fixed E_K
  can't. A test checks that such a space clears as one compartment, e^(−tP/θ), to 2e-3.
- A contact `bath` without a reference species, `bath(c, drive)` in the kit: its terminal voltage
  is its φ (an ideal salt bridge), so between two such baths V_right − V_left is the membrane
  potential a voltage clamp sets. `hodgkinHuxley({ T, linear: true })` writes Hodgkin and
  Huxley's own linear channels for a face (gated conductance links); clamped, the K⁺ current
  is g_K n∞⁴ (V − E_K) to 2e-3.
- A contact link for a surface velocity, `{ type: 'velocity', v }`: N_in = v (c_eq − c), between
  blocked and held. Thermionic emission over a Schottky barrier (v = A*T²/(F N_c)), within a
  few percent of J = A*T² e^(−φ_B/V_T)(e^(V/V_T) − 1), or a contact's surface recombination
  velocity: minority electrons through a base reach it as n₀(e^(V/V_T) − 1)/(W/D + 1/S) to 1e-5.
- Fixed: after a change of the left contact's voltage that the direct solve couldn't reach, the
  continuation (which ramped only the right contact, from where it was last solved) had nothing
  to ramp and returned the old state as converged. It ramps whichever contact moved now, and a
  ramp that runs no solve doesn't claim one.
- A solve that fails leaves the device as it was before it (it used to leave the failed iterate,
  from which nothing could go on: a battery host filling at a current, which has no steady
  state, couldn't then be advanced). An unreachable driven current's warning says when the
  opposite sign would have passed (a terminal's current is into the device), and when nothing
  passes at any voltage (a device that only stores charge: advance it in time instead).
- Fixed: a transient a second or more into a run, after a jump that a fine grid resolves in
  steps of 1e-13 s, stalled and stopped (a potential step at t = 1 s on a 5 nm grid): the error
  estimate was built from differences of absolute times, which keep only a few digits there.
  It uses the steps' own lengths now. An `advance()` that stops short, other than on its
  `budgetMs`, says where and why in its warnings.
- Fixed: a current-driven solve at I = 0 could find its voltage and still report no
  convergence (its last stage asked for a current of exactly zero). A solve that doesn't
  converge always carries a warning now.
- Fixed: with any port, `impedance()` read a contact's current at the contact, where behind a
  metal region it's a large conductance times round-off: a MOS capacitor with a metal gate and
  its channel held by a port read Re Z = −21.7 Ω·m² at its gate, at every bias. It's read across
  the device where no port's window intervenes now, and the back contact's nearly uniform
  response is split off as an open electrolyte's is.
- Fixed: a port holding a level inside another port's window left that port's ∂res/∂V in the
  row it replaced (a stale entry, ~1e-6 of the row: it touched only that port's impedance, or
  its circuit if floating). A gate's exponential rate law is held finite at absurd voltages,
  and `describe()` shows gated conductance links' gates.

## 0.12.0 (2026-10-07)

- Voltage-gated channels: a face's `gates`, each a fraction open with Hodgkin–Huxley kinetics in
  the voltage across the face (α and β in NeuroML's three standard forms), and permeabilities
  scaled by them (`gates: { m: 3, h: 1 }`). They're unknowns of the face's block; solutions report
  `interfaces[f].gates` and `.V`. `hodgkinHuxley()` in the kit writes the squid axon's channels,
  with GHK permeabilities matched to HH's conductances at rest. Under a voltage clamp each gate
  sits at α/(α + β) and relaxes at α + β; a squid axon's action potential follows the
  space-clamped HH equations with GHK currents to within 1 mV, and its impedance at rest (the
  inductance Cole measured, and a resonance near 60 Hz) the linearised equations' to 5e-4.
- Probes at a face: `{ interface: f, species }` (the flux through it) and `{ interface: f, gate }`.
- Gates on a port too: a membrane all along a region (an axon's), each node's gates following
  the voltage across the port's capacitance there, scaling its conductance links (Hodgkin and
  Huxley's linear channels, in species voltages). A squid axon's action potential propagates
  at 18.73 m/s, the cable equation's speed to 2e-3, out of the ions' drift along the axon
  (Hodgkin and Huxley computed 18.8). Face conductance links can be gated as well. A myelinated
  axon, nodes of Ranvier and internodes as regions with ports of their own, conducts
  saltatorily, each node's arrival within 1 µs of a compartmental cable's.
- A device with many held terminals (an axon's nodes of Ranvier, each a port) steps faster: only
  the floating terminals' ∂res/∂V and ∂I/∂x are carried through each assembly (all of them while
  `impedance()` runs).
- A demo, [myelin](demos/myelin.html): a bare and a myelinated axon of the same radius racing, the
  myelinated spike jumping node to node over ten times faster on a tenth of the current; myelin
  thickness, node spacing, nodal channel density.
- A demo, [propagation](demos/propagation.html): a spike travelling along a squid axon at Hodgkin
  and Huxley's speed, its local circuit, the ions that carry the current along the axon (mostly
  K⁺), and each ion's driving force along it; radius (speed as √a), temperature (heat block at
  35 °C), Na⁺ channel block.
- A demo, [action potential](demos/axon.html): a squid axon's spikes, their currents and gates,
  and each ion's driving force through them; the stimulus, temperature, and Na⁺ and K⁺ channel
  block.

## 0.11.4 (2026-10-06)

- Islands behind neighbours that conduct 1e20 times less than they do (across faces that hold μ̄
  level) solve too: there the summed row takes in the outside edge node's balance, so the flux
  is carried by the outside's own first segment rather than by the face flux unknown, which
  elimination could leave to the island's cancellation. Seen through one face of two, the
  island's level had been overshot twofold and swung between the two contacts', and the solves
  kept the plain, wrong level.

## 0.11.3 (2026-10-06)

- A region held only through its faces (a conductance there 1e-24 of its own, or neighbours that
  conduct 1e16 times less) solves to the right level in steady state: its balance summed over
  it, where its own fluxes cancel exactly, goes in as a bordered row once the plain solves
  converge with terminal currents that don't add up. Such solves had converged silently to a
  wrong level, their currents 2× or 16× off; now right to ~1e-15, and directly.

## 0.11.2 (2026-10-06)

- J·v is exact where a solve needs it (the impedance's GMRES and Newton's refined solves): the
  dilute kernels' terms are kept in difference form, each flux on its η difference, instead of a
  central difference of the residual, whose truncation swamped a nearly uniform response. The
  impedance is 15–40% faster, and keeps GMRES's answer always (its fallback to the factorised
  solve, for the old operator's noise, is gone). Newton keeps every converged refinement: a GaAs
  junction without recombination, its ~1.5e-7 A/m² carried by minority carriers from the
  contacts, had cold and warm solves 1e-3 apart, now 1e-14; a dim solar cell's open circuit on a
  0.25 nm grid takes 258 factorisations (790).
- Newton also refines its solves when its undamped updates grow twice running: a bipolar
  stack's floating base (its holes held ~1e14 more weakly than they move within it) now solves
  in steady state, cold and swept, and a MOS capacitor's gate step over 1000 s takes 146
  factorisations (408).
- A current drive whose floated solve fails is found by held solves instead: the voltage that
  passes the target, by regula falsi within the bracket the continuation found, to 1e-10 of the
  current. A closed Fe³⁺/Fe²⁺ cell driven near its limit (Fe³⁺ 20 orders below Fe²⁺, the floated
  system 33 digits short) now solves. Stress: 11,380 of 11,381.

## 0.11.1 (2026-10-05)

- The impedance splits its response into a uniform shift of each region and the rest, where a
  reading is uncertain: a redox electrode at its open circuit against a bath had its DC
  conductance 2–40% low (the current is the slope of levels uniform to 1e-14 of the response);
  now it matches the steady dI/dV to six digits. GMRES's answer is kept where it cuts the
  residual tenfold (a thousandfold had thrown away a real improvement at 0.01 Hz). The stress
  test now passes every impedance at 500 devices per family (11,375 of 11,381 scenarios).

## 0.11.0 (2026-10-05)

Robustness, found by a new stress test of random devices and fixed where it pointed.

**Testing, and what it says**
- `npm run stress`: random but plausible devices in five families (semiconductor stacks, MOS
  capacitors, electrolyte cells, electrodes, liquid junctions), each solved cold, swept warm,
  driven at open circuit or a current, stepped in time and probed by its impedance, every result
  judged by `check()` and by invariants (flat levels at equilibrium, warm and cold agreeing, the
  impedance passive and its DC limit the steady dI/dV). At 500 per family, 11,371 of 11,381
  pass. [docs/reliability.md](docs/reliability.md) has the families, the results and the known
  limits, and the README says plainly what it's tested on, linking there.

**Steady solves**
- A MOS capacitor without a channel port is robust: steady solves pin a species' level flat at
  its contact's where one contact alone reaches it and nothing else touches it, rather than find
  it through the bulk's ~1e3 minority electrons per cm³, and Newton refines its solves by GMRES
  (J·v from the residual) when the factorisation loses a mode, as it does an inversion layer's
  level. A 10 mV gate step reaches the low-frequency charge in ~70 steps instead of 80,000.
- A blocked species mobile through one region is held flat too, at the level its amount fixes:
  a closed Zn | ZnSO₄ | Zn cell past its limiting current, its sulfate excluded to 1e-25 mol/m³,
  solves at any voltage in 6 iterations, where it had failed from 0.7 V.
- The block factorisation perturbs a pivot that cancels to exactly zero (static pivoting): GaAs
  stacks without recombination (minority carriers ~1 per m³) solve, and so does a region held
  only by face conductances 1e17 weaker than its insides. `BlockTridiagonal` takes it as
  `staticPivots` (off by default).
- A terminal driven by a current (open circuit) is held at a voltage and marched to the
  crossing: a solar cell's V_oc in dim light on a 0.25 nm grid, a lit Schottky diode's, a redox
  electrode driven at a current needing 3 V. A current no voltage reaches fails with a warning
  saying so, and what the voltages passed.
- Cold starts: the light's continuation starts at 1e-30 of it and hands a stall to the
  pseudo-transient from its dimmer solution; Newton's divergence check spares the first update
  (a flooded scarce population predicts 1e10 thermal units, which damping covers); the
  pseudo-transient's dt grows gently after a failure and tries the direct solve once past the
  slowest diffusion time; a level start that fails from a layout made for the target voltage is
  laid out again. A lit GaAs phototransistor solves in ~1000 iterations (it failed after
  20,000), a cold start with immobile traps in 191 (280).
- Fixed: a flat level through a face between two regions left the face's flux undetermined.

**Impedance**
- The current is read where its error is least: between the contacts of a device without
  ports, the total current is the same through every cut, so it's read across each segment
  too. A p⁺n⁺ junction without recombination had its capacitance 9× off at 100 Hz (the contact's
  sum of huge terms read a constant conductance of round-off); now its DC limit is the steady
  dI/dV and its capacitance dQ/dV. Where a reading's estimated error passes 1%, the solve is
  continued to GMRES's floor: a redox cell's DC conductance had been 50% off.
- GMRES judges convergence by the true residual (its running estimate had passed unconverged
  solves: n-Si | KCl had a negative resistance at 0.01 Hz), keeps its best iterate, and falls
  back to the factorised solve where it can't converge (a strictly neutral electrolyte above
  1 MHz had gone non-passive).

**Transients**
- A first step whose error doesn't shrink with it (a jump starting a self-similar profile) is
  taken rather than shrunk until Newton fails: transients from a jump take 5–30% fewer
  factorisations. A bath beside a strictly neutral region resolves short steps (a sharp
  3 M | 1 µM junction starts at tol 1e-6). Haynes–Shockley takes 7% fewer factorisations.

**API and checks**
- `set({ ports: { gate: { V: 0.2 } } })` patches ports by name; a new kind of drive drops the
  old, as for contacts.
- `check()` counts a species at rest below 1e-8 of what it or its reactions' (non-metal)
  partners could carry: a supporting ion's round-off flux, or a nearly absent species' reaction
  noise, had read as unsteady.
- Fixed: a `set()` that widened a surfaced window started the new nodes at θ = 0.5, not `theta0`.

**Demos and docs**
- The MOS demo: C–V at any frequency from 10 µHz to 1 MHz against the ideal low- and
  high-frequency curves, the lifetime, the capacitance against frequency, and where the AC
  charge is answered. The demos' charts take `logx`.
- From an outside agent building an OECT from llms.txt: a volumetric capacitance C* as
  `area: 1/t, C: C*·t`; Bernards–Malliaras as the charge-sheet formula without its V_T term;
  `R` and the impedance in Ω with a geometry; which current `trace` records; tightening `tol`
  before fitting a decay's time constant.

## 0.10.0 (2026-10-04)

- The saturation demo has a thin-film transistor below the silver cell: the gate a capacitance
  along the channel, its output characteristic against the charge-sheet model, and the bands
  pinching off at the drain under the gate's flat Fermi level. Each part's controls stay with it
  as you scroll.
- The live demos have a pause button, and pause themselves while scrolled out of view or in a
  hidden tab, so a page doesn't keep a CPU busy unseen.
- Transients in strictly neutral materials start better. The state a first step starts from
  needn't satisfy the algebraic equations (φ, an interface's unknowns), which jump in the first
  instant, and the error estimate had read that jump as error on every later step; the history
  now starts after it. A closed device end is solved in the better-conditioned charge rows like
  an interior node, so an electrode spread through neutral water up to one takes its first step
  at tol 1e-8. `npm run hard` runs a set of hard cases, with a baseline of what still fails.
- Fixed: an electrode surface's coverages weren't under the time step's error control (a Langmuir
  filling ran 4% off its exponential at any tolerance).
- Fixed: species fed only through an electrode surface (A⁺ + e⁻ = S, S + e⁻ = B⁻) looked fed from
  outside, so a steady solve didn't conserve their total; it now does, with what the surface
  holds.
- Fixed: a reaction marked `bare` on a port without surface species ran at zero.
- Fixed: a held port level beside another port's electrode reaction in the same window (the hold
  listed first) left the currents unbalanced.
- Fixed: capacitances on overlapping windows overwrote each other; they add (two gates on a
  channel).
- Fixed: `set()` to a new grid or geometry carried the old conserved amounts (a slab's mol/m² as
  a sphere's mol).
- Fixed: `check()` on a non-planar device misread a face's per-area rates against the totals, and
  says mol/s, A and C there; `describe()` likewise for current drives.
- A capacitance alone driven by a current starts uncharged and charges in a transient; `solve()`
  says it has no steady state instead of running its voltage off. `ports[k].sigma` is zero where a
  contact holds that charge, as the solver counts it. Unit warnings cover ports' k0 and C.
- A floating electrode with a double layer starts uncharged (as the neutral starting composition
  has it) and charges toward its mixed potential; it had started at the mixed potential, charged,
  and the first step failed. A terminal driven by a current where every held terminal passes
  nothing (closed ends) is now an error that says to hold one that does, where Newton diverged.
- device.md's ports section is arranged by what a port does, its outputs in one table; llms.txt
  has templates for corrosion under a drop (an electrode port, a double layer, a cross-section, an
  Evans diagram) and a thin-film transistor.

## 0.9.1 (2026-10-04)

- Fixed: with a capacitance through a port in strictly neutral water, the most abundant ion leaked
  on time steps (about 1e-6 of it, so `check()` failed conservation): the row that stands for
  charge continuity there assumed the ions' net charge constant, but against a double layer it's
  −aσ. It now counts that change.
- `polarization(..., { surface: 'equilibrium' })`: an electrode port's coverages re-equilibrated at
  each potential, the steady-state curve with a passive film's active–passive peak.
- A solution warns when an electrode's surface is covered to a bare fraction below 1e-10, where the
  solver struggles; the docs' film example forms in minutes, and they say a film's μ° is fitted.

## 0.9.0 (2026-10-04)

- A capacitance spread through a port's window, `capacitance: { C, zeroCharge }` with `area`:
  σ = C (V − zeroCharge − φ) per area of electrode, aσ in each node's charge balance (neutrality
  at ε = 0), and the charging current through the port. A gate along a channel (the
  gradual-channel model of a thin-film transistor), an organic electrochemical transistor's
  volumetric capacitance, a porous electrode's double layer. With it an electrode port's potential
  swings continuously when a passive film covers the last active patch. Validated: a TFT against
  the charge-sheet model from below threshold through saturation, the de Levie impedance, the
  gate's low-frequency impedance against dQ/dV (test/capacitance.test.js).
- An electrode port's surface: species on its sites with a coverage at each node (Langmuir,
  μ = μ° + RT ln(θ/θ₀), Γ mol of sites per m² of electrode), named in its reactions like species,
  and `bare: true` for a reaction that runs only on bare metal (its rate times θ₀). Adsorbed
  intermediates, and passive films as a coverage: a film forms where potential and pH favour it,
  and blocks dissolution, the active–passive curve. Solutions report `ports[k].coverage`, and
  `polarization()` holds the coverages with the solution. Validated: the Langmuir isotherm, the
  charge to fill the surface, the blocking (test/coverage.test.js). The bare fraction is carried to
  full precision, and a nearly full surface's storage computed from its small complements. When a
  film passivates an electrode's last active patch the mixed potential must jump; give the port a
  capacitance (its double layer) and the swing is continuous.
- An electrode spread through a port now requires its region strictly neutral (ε = 0): its double
  layers are below the grid, and with ε > 0 its reactions left the solution charged (steady solves
  converged, transients didn't).
- `build()` passes `geometry` through. The docs give the shapes of a port's outputs.

## 0.8.0 (2026-10-03)

- Devices with a cross-section A(x): `geometry: { type: 'spherical' | 'cylindrical', r0 }` (shells
  about a centre r0 to the left of x = 0) or `{ area: { x, values } }`. The finite volumes carry
  totals through A: boxes hold ∫A dx, faces and contacts pass their flux times A there, and each
  segment is weighted by its length over ∫dx/A, so steady diffusion between nodes is exact for
  any A. Currents, fluxes, amounts and charges are then totals (planar devices are unchanged, bit
  for bit, with A = 1 m²). Validated in test/geometry.test.js: steady diffusion to a sphere and a
  cylinder (exact on any grid), Cottrell's transient with the spherical term, uptake by a sphere
  from its surface (centre at x = 0), Debye–Hückel screening around a charged sphere.
- Electrodes spread through a window: a port with `reactions` (Butler–Volmer per area, as at a
  face) against its terminal species, a metal's carrier at the port's level, and an `area` of
  electrode per volume (a number or a profile). Used for the metal under a thin film of electrolyte,
  a crevice's walls, or a porous electrode with a well-conducting matrix. Held, it's a
  potentiostat; floating (or with nothing else to carry current) it settles at the mixed potential,
  each spot a net anode or cathode. Solutions report each reaction's rate along the window
  (`ports[k].x`, `.rates`, `.area`). Validated against the transmission line (linear kinetics along
  a bar) and the Wagner–Traud mixed potential of two couples. A floating electrode port starts at
  its mixed potential.
- `polarization(device, sol, where, V)` in the kit: an electrode's reactions (at a face, or at a
  spot of an electrode port) with the solution beside it held and the metal's level moved. It
  returns each reaction's rate, partial current (anodic positive) and redox level, and the net
  current: the curves of an Evans diagram. It uses the solver's own law, so at the metal's actual
  level it gives the solved rates (test/polarization.test.js, with the Tafel slopes).
- Fixed: a port that only exchanges (an O₂ supply from the air, say) in a strictly neutral region
  could stop a transient at its first step. With a face reaction setting the electrolyte's
  potential, as in a drop of salt water on iron, no step converged. The port's window is now solved
  in the same well-conditioned terms as the rest of the region. A port that exchanges only neutral
  species carries no current, and driving one by a current is now an error.
- The liquid-junction demo says what the minority ions show: Mg²⁺'s plateau, and K⁺ riding the
  front.

## 0.7.0 (2026-10-03)

- Membranes as faces: a species law `{ type: 'permeability', P }` (m/s), electrodiffusion through a
  thin membrane in a constant field (Goldman–Hodgkin–Katz), exactly zero where μ̄ is level. With a
  capacitive φ law, a cell membrane is one face. Validated: the GHK resting potential with the
  solutions strictly neutral or their diffuse layers resolved, the GHK current–voltage curve, the
  agreement with a resolved 5 nm lipid layer, and charging at the GHK current over C.
- Saturating kinetics for face reactions, `vmax` and `K` in place of `k0` and `alpha`:
  r = vmax Π (c/(c + K))^|ν| (1 − e^(−A/RT)), Michaelis–Menten in each substrate and exactly zero
  at A = 0, for pumps and transporters. The Na⁺/K⁺-ATPase on a membrane face, in a closed cell,
  gives the Mullins–Noda potential, stalls at its static head (3Δμ̄_Na − 2Δμ̄_K = ΔG_ATP) with
  nothing leaking back, and leaves Donnan equilibrium when off.
- `check()` keeps a ledger per compartment: a species a face lets through only by a law or a
  reaction gets one on each side, with what crosses as a term (in a cell: Na⁺ leaks in as fast as
  it's pumped out). `describe()` lists what crosses a face by a law, and saturating kinetics.
- Liquid-junction potentials in the kit: `henderson()` (the closed form JPCalc uses) and `planck()`
  (the steady state of constrained diffusion, by shooting on the Nernst–Planck equations), for any
  mixture and valences. Validated in test/junction.test.js against JPCalc's published values,
  against driftlet's own steady junction (to 2e-4 mV) and the stationary codes LJPcalc and JLJP
  (to 0.02 mV: a patch pipette's K-gluconate, ZnCl₂ | KCl), and a free-diffusion junction grown
  from a sharp boundary, constant in time and between the two.
- A liquid-junction demo: a patch pipette's junction growing from contact against Henderson,
  Planck and the published values (and why patch-clampers correct for it), ZnCl₂ | KCl, and a
  3 M KCl salt bridge whose charged frit turns permselective against a dilute sample.
- A charged-nanochannel demo: a channel 10 nm across with charged walls, resolved along its axis;
  Donnan layers at each mouth that overlap in dilute salt (where Poisson's equation, not
  electroneutrality, sets what it does), K⁺ selectivity, the conductance plateau at low salt
  against Teorell–Meyer–Sievers, the zero-current voltage across a salt gradient, and
  rectification by a one-sided or bipolar charge, explained as the pn junction it is. Validated in
  test/pore.test.js (TMS to 1% in a long pore, its end correction ∝ 1/L, rectification).
- A resting-potential demo: a closed cell's leaks and the Na⁺/K⁺ pump on a membrane face, V_m
  against Mullins–Noda, each ion's species voltage stepping at the membrane by its driving force
  V_m − E_i, and the run-down to Donnan equilibrium over hours when the pump stops.
- Probes read φ (`{ x, quantity: 'phi' }`), so a transient's trace can carry a membrane potential.
- A floating contact starts where the start's own composition puts it, so a cold-start
  transient's first trace point is right (it read the held terminal's voltage, off by a Nernst
  term or more).
- Docs: a cell-membrane template in the agent guide; a bath terminal reads its reference
  species' level, not φ; which contact a species without `c0` starts from; `set()` from
  `onsolution`; the worker from a CDN; the data's 25 °C. All from the membrane agent's debrief.

## 0.6.0 (2026-10-03)

- `check(device, sol)` in the kit: checks on a solved device, each saying what it compared, as
  `{ ok, items, text }`. Convergence; warnings and likely unit slips (marked `?`, a prompt to
  look); in a steady state, each species' ledger (what every terminal and reaction brings in or
  uses, summing to zero: J = F(G − R) under light, and where each carrier went), and the device
  solved again on a grid twice as fine (current, floating voltages and region charges to 1 %, the
  change estimating the grid's own error); in a transient, conservation. Validated on a lit
  silicon cell (the ledger term by term, and the grid change against a 16× finer solve) and run
  over every steady solve in the test suite without false alarms. The solar demo shows it, and
  its grid is finer for it (hmin 1 nm had left the emitter's depletion charge 2 % off).
- Solutions say whether they're a steady state, `sol.steady`.
- Curated devices in the kit: `perovskiteCell({ V, light, bulkSRH })`, IonMonger's default
  perovskite solar cell as a plain definition, and `PEROVSKITE_SCANS`, IonMonger's own J–V scans
  of it (current at 0 V, each sweep's P_max and V_oc, hysteresis index, at seven rates, with and
  without bulk SRH), so a demo can start from benchmarked numbers and show its agreement;
  `hysteresis(trace, turn)` reads the same figures from a run. The IonMonger test and the
  perovskite demo use them. Asked for by the perovskite agent.
- `recorder(device, { every, times, probes })` in the kit: a transient run in pieces (or in
  animation-frame budgets, step for step the same run), with the trace joined up from the start
  and whole solutions kept as frames at set times, `frame(t)` to scrub. Validated on a
  diffusion mode decaying frame by frame at its analytic rate.
- The IonMonger benchmark measures the loop as a current difference where it's gentle and as a
  voltage offset where it's steep: within 0.010 mA/cm² and 0.31 mV at all seven rates (the
  0.06–0.23 mA/cm² reported before was a sub-millivolt shift read on a near-vertical curve). The
  perovskite agent spotted this in its own comparison.
- Each solution reports its bulk reactions' rates, `bulkReactions[k]`: per node, per region and in
  all, so charge balance J = F(G − R) is a few lines (validated to 1e-6 in the IonMonger cell).
- `describe()`: an SRH reaction is listed where it runs; a region's Debye length counts only the
  mobile charge `c0` leaves unbalanced on top of `c0` itself (0.5.0 counted a perovskite's
  vacancies twice, once as themselves and once as balancing their background).

## 0.5.0 (2026-10-03)

- SRH (trap-assisted) recombination as a rate law, without trap species: in the bulk,
  `srh: { material: { tauN, tauP, n1 } }` in place of `kf`; at a face, `srh: { vn, vp, n1 }` in place
  of `k0` and `alpha`, with n and p each from its own side. Equilibrium stays exact (p₁ comes from
  the state), and at a face it saturates at one carrier's capture, as interface recombination in
  perovskite cells does. Validated against the bulk closed form, and at a face against the law
  and the device's bookkeeping.
- A benchmark against an independent code: perovskite J–V hysteresis (mobile iodide vacancies,
  interface SRH, with and without bulk SRH) against IonMonger at seven scan rates from 1 mV/s to
  1 kV/s. Hysteresis index within 5e-4, maximum power within 0.02 mW/cm², V_oc within 1 mV, whole
  loops within 0.06 mA/cm² (0.23 at the fastest-changing sweep). Three species and two face SRH
  laws; each scan takes 0.3 s against IonMonger's 50.
- Immobile reacting species (trap states, fixed charge states) in steady solves: a combination
  of them that reactions conserve is now held node by node, as a row of that node's block, so
  the steady solve goes direct instead of through giant time steps (which stalled when the
  traps relaxed slowly). Validated: explicit immobile traps recombine at exactly the SRH rate
  in steady state, their totals kept at every node.
- A cold start with slow ions alongside fast physics (a perovskite's vacancies with carriers and
  charged traps) solves: the steady solve's fallback time steps, and `advance()`'s first step,
  shrink as far as the device's fastest time scale, not only to a fraction of the slowest (or
  of the run's length). Both used to give up a few decades short. Reported, with a minimal
  repro and the diagnosis, by an agent building a perovskite demo.
- `traces()` and `bandDiagram()` take `energy: true`: the familiar band diagram, energy up (E_c,
  E_v and the quasi-Fermi levels in eV). The semiconductor demos (pn, solar, MOS, organic) now use
  it; species voltages stay where ions share the diagram (electrochemistry, the perovskite's
  vacancies). The docs say which picture to use when.
- Demo: a perovskite solar cell's J–V hysteresis from mobile iodide vacancies (IonMonger's default
  cell), with a scan-rate control, the vacancies' layers at each face, the steady and frozen
  limits, and the hysteresis index against scan rate with IonMonger's own points.
- `set()`: a contact given `V` or `I` drops the other, so `{ I: 0 }` switches a held contact to
  open circuit (it used to need `V: undefined` alongside).
- A face reaction's α may be 0 or 1 (plain mass action one way), not only strictly between.
- Each solution's `interfaces[f]` names the regions it joins, `left` and `right`.
- `describe()`: a region's Debye length counts the carriers balancing its doping, not only what
  `c0` lists (a trace in `c0` made a doped layer's read in µm instead of nm).

## 0.4.0 (2026-10-03)

- A bulk reaction's `kf` can vary with position: a profile `{ x, values }` against the device's
  x, averaged over each node's box, so the total rate is exact on any grid.
- `photogeneration({ material, flux, alpha, mu, from, to })` in the kit: Beer–Lambert absorption
  as a generation reaction. Validation: J_sc in an n⁺p cell against collection theory, from blue
  light lost in the emitter to red collected deep in the base.
- `photogeneration({ makes })`: light can make something other than e⁻ + h⁺, such as an
  exciton. Validation: an organic bilayer, with excitons as a third, neutral species diffusing to
  the donor/acceptor interface and splitting there (J_sc against exciton-diffusion theory), and
  V_oc from charge-transfer recombination by detailed balance.
- Solar demo: a wavelength control (silicon's α from Green 2008), the spectral response against
  collection theory, and where the light makes pairs and which are collected.
- Demo: an organic bilayer solar cell, where light makes neutral excitons that diffuse to the
  donor/acceptor interface and split; J_sc against donor thickness and V_oc by detailed balance,
  each against theory.
- The cold start of an undoped layer holding one carrier (which can't be neutral) keeps the
  running φ instead of driving it to an extreme, and φ starts from a pinned or gated left
  contact. Such devices used to fail on the first Newton step.
- A species fed only by reactions still needs a `c0`, and the error now says that it's only a
  starting point, not an amount it keeps.
- Breaking: a `c0` profile is `{ x, values }` (it was `{ x, c }`), the same shape as a `kf`
  profile and, against t, a waveform.

## 0.3.0 (2026-10-02)

- `c0` profiles: a region's initial concentration can be tabulated against x (`{ x, c }`,
  piecewise linear), so a transient starts from a state it didn't simulate reaching (a packet
  already injected, a gradient laid down). Species without a `c0` start at their contact's
  level, with $`\phi`$ chosen node by node for local neutrality, so a counter-ion or majority
  carrier follows the profile. `describe()` warns when a profile misses its region.
- A `c0` given for a species a contact feeds is now its starting state (it was ignored, the
  species starting at the contact's level); without one, nothing changes.
- Validation: a Gaussian hole packet laid down as a profile drifting, spreading and decaying at
  the ambipolar rates, and a salt's diffusion mode decaying at $`\pi^2 D/4L^2`$.
- Results identical across JavaScript engines: powers go through `powi`/`powr` instead of `**`,
  whose last bit differs between Node 22 and 24 (enough to flip a step in a long adaptive run).
- Demos: a real-browser smoke test in CI (`npm run smoke`: every page in headless Chrome, every
  control moved to its extremes), and screenshots retaken by script (`npm run shots`), each
  page held at a set state (`?t=…`, and `?phase=…` on the impedance page).

## 0.2.0 (2026-10-02)

- Kit sources: waveforms (`pulse`, `square`, `triangle`, `ramp`), `injector()` (a port that
  spreads a driven current evenly) and `recombination()` (from a minority lifetime).
- `advance(t, { probes })`: a species' concentration or voltage at points inside the device, in
  the trace after every accepted step.
- `units` back out of SI: `toPerCm3`, `toMolar`, `toCm2PerS`, `toUm`, `toNm`.
- Validation: a Haynes–Shockley pulse against the ambipolar drift, spread and decay.
- Newton damps only the device's potentials, not floating terminal voltages (which enter
  linearly): a current-driven port switching off, its voltage collapsing by hundreds of volts,
  no longer crawls. A strong Haynes–Shockley pulse takes 55% fewer factorisations.
- Device reference: what a driven port's conductance actually does (it shares the current; a
  small one spreads it evenly), correcting 0.1.2's note.
- `llms.txt`: good habits (draw the closed form alongside, stay honestly in 1D).

## 0.1.2 (2026-10-02)

- The npm package ships `docs/`, so the guide's references to them work from the tarball.
- `llms.txt`: a transient template (a pulse injected through a current-driven port, sampled
  downstream), and what makes transients expensive.
- Device reference: what a driven port's conductance does, and that a device's time starts at 0
  when it's made (`solve()` doesn't advance it).
- `SPEC.md`, the original design brief, retired: it lives on in the docs and in git history.

## 0.1.1 (2026-10-02)

The first release published by trusted publishing from GitHub Actions, with provenance.

- README: what driftlet is for, links to the live demos, npm and docs, and every demo in the
  table of fields.
- Guide and demos: the current a species carries flows down its species voltage $`V_i`$, so
  electrons and anions themselves move uphill (they said every species drifts down); $`V_i`$ is
  described as a display of $`\bar\mu_i`$, which also covers neutral species.

## 0.1.0 (2026-10-02)

First release: the solver (`driftlet`), the porcelain (`driftlet/kit`) and level-diagram plots
(`driftlet/plot`), with TypeScript declarations, the validation suite, and live demos.
