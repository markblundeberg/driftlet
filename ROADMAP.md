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
  conductance links, Butler–Volmer electrode reactions.
- Interfaces: blocking, conductance, Butler–Volmer ion/electron transfer.
- Circuits: voltage, galvanostatic (incl. open circuit), load resistor.
- Bulk mass-action reactions.
- Backward-Euler transients with step halving, and conservation bookkeeping.
- Statistics: one general per-material interface, c(ζ) with its capacitance matrix, and the
  excess-potential Scharfetter–Gummel flux. Fermi–Dirac, lattice gas (Bikerman crowding, site
  filling), Redlich–Kister, Debye–Hückel, custom functions.
- Neutral-combination (insertion) statistics for ε = 0 hosts, from an isotherm or a tabulated
  OCV curve, with φ undefined inside.
- Time and frequency: variable-step BDF2; adaptive `advance(tEnd, { tol, budgetMs })` for
  animation; small-signal impedance with complex block Thomas and complex profiles.
- Source continuation for cold starts at bias.
- A first performance pass: a benchmark script (`npm run bench`), assembly without per-node
  allocation (about 2× faster), Newton started from the extrapolated state in adaptive
  transients (a third fewer iterations), and a tighter-damping retry for large jumps.
- Metal regions: a Fermi level and a conductivity only (mixed-form Ohm's law, surface charge as
  sheets), capacitive faces aligned by work function, electrode reactions at internal faces,
  bipolar electrodes.
- Internal ports: reservoirs attached to a window of interior nodes (held levels, volumetric
  conductance or exchange), e.g. grounding a 1D MOS channel.
- Transport extras: imposed advection (exact exponential fitting) and current-free eddy mixing
  (an Onsager term projected to carry no current).
- Validation suite (see the README), and docs.

## Next

1. **Porcelain, for one-shot demos.** An optional layer (e.g. a `driftlet/kit` subpath export,
   still dependency-free) that writes plain specs, so users and LLM agents start from something
   correct. It never hides a physical choice: everything it produces is ordinary, inspectable
   spec data.
   - Contact and charge shorthands: `ohmic(V)`, `bath(…)`, `electrode(reaction, …)`; region
     charge as `{ donors }`, `{ acceptors }`, `{ fixed: { c, z } }` instead of `… * FARADAY`
     with a sign to remember.
   - A live wrapper for sliders: frame throttling, warm starts, dropping stale requests,
     ramping across big jumps, keeping the last good solution; optionally Worker-backed.
   - `describe()` and unit-slip warnings: Debye lengths against the grid, time constants,
     conductivities, built-in or open-circuit potentials, "this D looks like cm²/s".
   - Half-reactions as objects (participants, electrons, fixed participants' μ given
     explicitly). They build electrode reactions, serve as vacuum-alignment anchors, and give
     the spatial profile of their implied electronic level (`sol.level(halfReaction)`), with SHE
     as one instance: the ESBD redox-level view in one line.
   - A voltage ledger: the terminal voltage as named contributions that sum exactly. The honest
     form sums species-voltage steps along the current path, region by region and face by face;
     φ-based textbook terms (ohmic drop vs junction potential) are flagged as conventional.
   - Plot-ready traces: segments broken where a species is absent, region bands, label-ready
     interface steps, suggested ranges.
   - A small data library with a source on every entry (aqueous ions: z, D, μ° at 25 °C on
     table conventions; common semiconductors; metals), kept small and vetted, since a wrong
     library number is worse than a wrong demo number. Then recipes built from it (pn junction,
     galvanic and concentration cells, double layer, membrane), returning editable specs.
   - An agent guide (`AGENT_GUIDE.md` or `llms.txt`): conventions, a pinned CDN import,
     copy-paste demo templates, common mistakes. Its code blocks run in the docs test, and error
     messages say how to fix the problem.
2. **Ship 0.1.0.** npm package and jsdelivr, live demos (code sandboxes, screenshots in the
   README), JSDoc types and a `.d.ts`, CI.

## Later

- **Cross-species transport coefficients together with cross chemical capacitances.** They're
  the two halves of one Onsager / Jamnik–Maier network, so one shouldn't come without the other.
  This is also the general mixed conductor with n > 2 mobile species, beyond the metal (n = 1)
  and insertion-host (n = 2) cases.
- **More performance.** A variable block size where species are absent (identity rows cost a
  full block today; a metal node needs only its Fermi level and current), and benchmarks
  tracked across versions in CI.
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
- **Ports as circuit terminals.** Internal ports exist (fixed-voltage reservoirs); letting a
  port's voltage follow a terminal, or float in the circuit, would come next.
- **Interface states** (charge that depends on the local μ̄_e, i.e. Fermi-level pinning) and
  thermionic-emission links at heterojunctions.
- Non-isothermal transport (Soret, Seebeck/Peltier, a heat equation), which stays local.

## Out of scope, for good

2D/3D; non-local couplings such as recirculation or one well-mixed reservoir feeding several
faces (a region with large mixing is the local stand-in); anything that breaks the
block-tridiagonal structure.
