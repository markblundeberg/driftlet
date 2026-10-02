# Roadmap

Where driftlet is, and where it's going. Scope stays sharp: 1D, local physics only (everything
must keep the Jacobian block-tridiagonal), pure JavaScript with no dependencies, fast enough
for live demos.

## Done

- Solver core: graded grid with doubled interface nodes, block Thomas, Scharfetter–Gummel in
  quasi-Fermi (expm1) form, compensated double-double μ̄, damped Newton.
- Equilibrium and steady transport: direct steady solves when nothing is conserved, and
  exactly conserving giant steps when something is.
- Materials, regions and explicit per-interface alignment; interface φ laws (dipole, neutral,
  Helmholtz); strictly neutral (ε = 0) materials.
- Contacts: fixed levels with offsets, ohmic contacts, baths, gates and Stern layers,
  conductance links.
- Interfaces: blocking, conductance, and Butler–Volmer reactions with signed participants on
  either side (each rate an unknown of the face block); electrodes as conductor regions.
- Circuits: voltage, galvanostatic (incl. open circuit), load resistor.
- Bulk mass-action reactions.
- Strictly neutral regions solved in better-conditioned unknowns on transient steps (current
  continuity for one balance, φ̂' and η − zφ̂ as unknowns): a 3 M KCl salt bridge with 1 µM traces
  of every other ion, between eddy-mixed solutions, runs from its unmixed start through days.
- Conserved combinations of reacting stretches (moieties) as constraint rows in the direct
  steady solve, as spectators already were: a closed cell with complexation solves straight to
  its steady state, cold or warm.
- Backward-Euler transients with step halving, and conservation bookkeeping.
- Statistics: one general per-material interface, c(ζ) with its capacitance matrix, and the
  excess-potential Scharfetter–Gummel flux. Fermi–Dirac, lattice gas (Bikerman crowding, site
  filling), Redlich–Kister, Debye–Hückel, custom functions.
- Neutral-combination (insertion) statistics for ε = 0 hosts, from an isotherm or a tabulated
  OCV curve, with φ undefined inside.
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
  conductance or exchange), e.g. grounding a 1D MOS channel, or a wire to a whole metal.
- A second-law validation: in steady state the free energy brought in through every terminal,
  Σ N μ̄ over every species, equals the dissipation inside (fluxes down their own μ̄ drops,
  reactions by their affinities, interface and conductor resistance), each term ≥ 0, exactly
  for the discrete equations.
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

## Next

1. **Harden the core** (ongoing), so the porcelain doesn't force refactors.
   - Floating conductors in η form. A floating metal is fixed by the mixed form, but the same
     cancellation awaits a floating mixed conductor or semiconductor region held only by face
     reactions or recombination, in steady state (in transients, storage anchors it): a fast
     species across many segments, with a mode only weakly held from outside. Candidates are
     a Grassmann–Taksar–Heyman-style elimination that carries leakage separately, or the mixed
     form for the fast species only. (A failed steady solve already says where it lost its digits.)
   - A minority population held only weakly: the MOS capacitor without its channel port,
     whose inversion electrons reach the back contact only through a p-type bulk of ~1e3 cm⁻³
     (minutes to fill, so a high-frequency C–V at 1 Hz). The impedance now gets this right
     (GMRES with the residual's J·v), but the steady solve is fragile: stepping the gate up warm
     from depletion fails at the onset of inversion ("exactly singular"), and cold starts
     converge at some gate voltages and grids, not others (a grid shifted by round-off can flip
     it). Its transients struggle too: after a 10 mV gate step, `advance()` takes ~66,000 steps
     to reach 1e5 s, and the gate charge overshoots the low-frequency value (1.01 against 0.885
     × C_ox·δV) with a residual current of ~1e-9 A/m² through the oxide. Once robust, the MOS
     demo can get a frequency control sweeping from the low-frequency C–V to the high-frequency
     one.
   - Strictly neutral regions on very short steps, where storage dwarfs fluxes. Interior nodes
     and the edges of neutral faces are solved in better-conditioned unknowns (φ̂', η − zφ̂),
     but at a neutral face the two edge nodes still pass every species' flux through its own
     unknown, which carries storage-sized values, so round-off returns at ~1e-5 thermal units
     on steps of ~1e-11 s with 3 M against 1 µM. A perfectly sharp junction between such
     solutions then can't start (a pre-mixed zone a few layers wide can). The fix is a
     supernode: each region still supplies its own edge block, and the interface code merges
     the two edge nodes and the face's unknowns into one larger block, so every quantity at
     the face couples to every other. Each species' two edge balances can then be summed, with
     the face fluxes cancelling exactly, into a row whose outside couplings reach only the
     nodes just beyond. The block at each face grows (about threefold) and so do its
     off-diagonal neighbours; the rest of the matrix is unchanged, and the layout stays
     block-tridiagonal. It may simplify the face code generally.
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
   - A worked membrane-potential example with a validation against Goldman–Hodgkin–Katz in its
     own regime (a thin membrane with ε > 0 and few ions, so the field is constant), not the
     electroneutral one a guinea pig reached for.
   - A written worked example of an electrochemical cell with its voltage ledger (the Daniell
     demo has one: the couples' levels, the electrolyte's share, each electrode's
     overpotential, summing exactly to the terminal voltage). A ledger's division is specific
     to its device family, so it belongs in a worked example, not the library; the spatial
     diagram stays the primary picture.

3. **Ship 0.1.0.** npm package and jsdelivr (the agent guide's pinned CDN version is checked
   against `package.json` once it's set), live demos (code sandboxes, screenshots in the
   README), JSDoc types and a `.d.ts`, CI. (Tests already check what browsers and CDNs need:
   relative `.js` imports only, no Node or DOM globals, and a definition solving identically on
   a worker thread. A real-browser smoke test remains.)

## Later

- **Generation profiles.** Position-dependent generation (Beer–Lambert absorption) for solar
  cells and photodetectors: today a reaction's rate constant is uniform within each material.
- **More interface kinetics.** Marcus–Hush–Chidsey rates (curved Tafel plots, saturation at large
  overpotential); surface species (adsorbed intermediates with a coverage), so multi-step
  mechanisms such as hydrogen evolution emerge from elementary steps; and a custom forward rate
  r_f(state) as an escape hatch, which the solver multiplies by (1 − e^{−A/RT}) so equilibrium
  stays exact. Not a bare current–overpotential curve, which loses the concentration dependence
  and can break detailed balance.

- **Cross-species transport coefficients together with cross chemical capacitances.** They're
  the two halves of one Onsager / Jamnik–Maier network, so one shouldn't come without the other.
  This is also the general mixed conductor with n > 2 mobile species, beyond the metal (n = 1)
  and insertion-host (n = 2) cases.
- **Concentration unknowns for dilute regions, as an option.** A hybrid basis: conductor
  regions keep their Fermi level, but a dilute region may ask to be solved internally in
  {φ, c_i} instead of {φ̂, η_i}, where storage and neutrality are linear and a strongly
  diluted species (1 µM against 3 M) doesn't cost digits in its level. Transparent to the
  API: the solution still reports μ̄, V_i and c, and faces still align levels.
- **More performance**, measured against the committed benchmark baseline: Newton-iteration
  counts (above), and specialised small-block elimination.
- **More statistics.** Gaussian and exponential densities of states (disordered and organic
  semiconductors), species occupying several lattice sites, several models on one species
  (e.g. crowding plus activity coefficients), and cross-model shared sites.
- **Demo extras.** A tiny optional plotting entry point (`driftlet/plot`), a deterministic mode
  for screenshots, and concentration/flux sampling to drive particle animations consistently
  with the model.
- **Phase separation (Cahn–Hilliard).** A non-convex free energy (e.g. LFP's miscibility gap)
  makes c(ζ) multivalued. It needs c as an extra unknown plus a gradient-energy term, which is
  still local and block-tridiagonal.
- **Graded materials.** μ°(x), c_ref(x), ε(x), D(x) varying within a region, with the full
  `zφ̂ + μ°/RT − ln c_ref` change across each segment. φ's gauge freedom then becomes a
  continuous function, as in ESBD's inhomogeneities topic.
- **A sub-grid Gouy–Chapman interface law.** Diffuse layers treated analytically where λ_D ≪ h,
  for macroscopic devices with real double-layer charge on any grid. It tends to `neutral` as
  ε → 0.
- **Driven species fluxes.** A terminal driving one species' molar flux instead of the charge
  current (a gas feed at a fixed rate, a neutral species injected), the flux-side twin of a held
  μ for neutral species.
- **Interface states** (charge that depends on the local μ̄_e, i.e. Fermi-level pinning) and
  thermionic-emission links at heterojunctions.
- Non-isothermal transport (Soret, Seebeck/Peltier, a heat equation), which stays local.

## Out of scope, for good

2D/3D; non-local couplings such as recirculation or one well-mixed reservoir feeding several
faces (a region with large mixing is the local stand-in); anything that breaks the
block-tridiagonal structure.
