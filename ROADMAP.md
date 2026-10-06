# Roadmap

Where driftlet is, and where it's going. Scope stays sharp: 1D, local physics only (everything
must keep the Jacobian block-tridiagonal), pure JavaScript with no dependencies, fast enough
for live demos.

## Done

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

## Next

1. **Harden the core** (ongoing), so the porcelain doesn't force refactors.
   - Floating conductors in $`\eta`$ form. A floating metal is fixed by the mixed form, but the same
     cancellation awaits a floating mixed conductor or semiconductor region held only by face
     reactions or recombination, in steady state (in transients, storage anchors it): a fast
     species across many segments, with a mode only weakly held from outside. Candidates are
     a Grassmann–Taksar–Heyman-style elimination that carries leakage separately, or the mixed
     form for the fast species only. (A failed steady solve already says where it lost its digits.
     Newton's refined solves, which recover a MOS inversion layer's lost level, may cover it.)
     `npm run stress` finds the same thing in a bipolar stack's floating base (a p⁺ base between
     light n layers, its holes held ~1e14 weaker than they move within it: 25 digits lost even
     at 10 mV), and in a closed Fe³⁺/Fe²⁺ cell driven near its limit, where Fe³⁺ falls 20 orders
     below Fe²⁺ (33 digits). The level of such a population could be deflated: a bordered unknown
     for its uniform shift, its row the summed balance, read from the residual as the impedance
     reads its current.
   - Strictly neutral regions on very short steps, where storage dwarfs fluxes. Interior nodes and
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
     generally.
   - The impedance's δx in double-double. A redox electrode at its open circuit against a bath
     responds to its terminal with a nearly uniform shift of every level, and its current is the
     slope, which some cells need to better than 1e-14 of the response: plain doubles can't hold
     it, and the finite-difference J·v rounds at ~1e-11 (φ̂ in the high word, the exponentials'
     arguments summed before exp). Their DC conductance can come out 2–40% low, a few in a
     thousand of `npm run stress`'s electrodes. Solving for the deviation from a per-region
     uniform shift, whose J·v is exactly zero inside each region, may be the way.
   - Fewer Newton iterations where the benchmarks show many: the pn transient (about 1000
     factorisations for 100 ns) and large warm jumps.
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
   - Voltage-gated channels (an action potential): a permeability that depends on gating
     variables with their own kinetics (Hodgkin–Huxley), which would be extra unknowns at a face.
     Today a page can change `P` with `set()` between steps.
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
     - **Each terminal against SHE** in `describe()` and in failure messages, in the solution beside
       it. An agent held the iron at V = 0 beside a bath referenced on Cl⁻ (NBS data), 1.37 V
       anodic, and lost time finding out why.
     - **Cheap mid-run parameter changes**: `set()` of a reaction's k0 between steps restarts the
       stepper (a page updating it every frame ran 30× slower). Drive changes already have a fast path.
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
     - **Energy-up band diagrams** in `driftlet/plot` for electrons and holes (E_c, E_v, quasi-Fermi
       levels in eV), with species voltages kept for what that picture can't show: ions, and
       equilibrium between species.
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
   - Before 1.0, definitions as a contract: an optional `version` field, and errors that name
     renames ("`foo` became `bar` in 0.3"), since agents will copy old definitions long after the
     API moves on.
   - A written worked example of an electrochemical cell with its voltage ledger (the Daniell
     demo has one: the couples' levels, the electrolyte's share, each electrode's
     overpotential, summing exactly to the terminal voltage). A ledger's division is specific
     to its device family, so it belongs in a worked example, not the library; the spatial
     diagram stays the primary picture.

## Later

- **More interface kinetics.** Marcus–Hush–Chidsey rates (curved Tafel plots, saturation at large
  overpotential); surface species (adsorbed intermediates with a coverage), so multi-step mechanisms
  such as hydrogen evolution emerge from elementary steps; and a custom forward rate
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
- **Demo extras.** Concentration/flux sampling to drive particle animations consistently with
  the model.
- **Graded materials.** $`\mu^\circ(x)`$, $`c_{\mathrm{ref}}(x)`$, $`\varepsilon(x)`$, $`D(x)`$ varying
  within a region, with the full $`z\hat\phi + \mu^\circ/RT - \ln c_{\mathrm{ref}}`$ change across
  each segment. $`\phi`$'s gauge freedom then becomes a continuous function.
- **A sub-grid Gouy–Chapman interface law.** Diffuse layers treated analytically where
  $`\lambda_D \ll h`$, for macroscopic devices with real double-layer charge on any grid. It tends to
  `neutral` as $`\varepsilon \to 0`$.
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
