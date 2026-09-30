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
- Validation suite (see the README), and docs.

## Next

1. **Time and frequency.**
   - Adaptive time steps from a BE/BDF2 error estimate, and BDF2 itself.
   - `advance(tEnd, { budgetMs })` for animation: as many steps as fit a frame budget.
   - Small-signal impedance, (J + iωM)·δx = b, with complex block Thomas, reporting Z(ω) and
     complex profiles (Warburg, double-layer and Maxwell–Wagner demos).
2. **Transport extras.** Imposed advection v(x), and current-free eddy mixing: an Onsager
   term −(D_mix/RT) P ∇μ̄, with P = C − C z zᵀ C/(zᵀ C z), that mixes composition without
   conducting and vanishes exactly at equilibrium.
3. **Performance pass.** Benchmarks tracked across versions. Cheaper warm re-solves (2–5
   Newton iterations for small parameter changes), allocation-free assembly, and a variable
   block size for absent species.
4. **Ship 0.1.0.** npm package and jsdelivr, live demos (code sandboxes, screenshots in the
   README), JSDoc types and a `.d.ts`, CI.

## Later

- **Cross-species transport coefficients together with cross chemical capacitances.** They're
  the two halves of one Onsager / Jamnik–Maier network, so one shouldn't come without the other.
- **More statistics.** Gaussian and exponential densities of states (disordered and organic
  semiconductors), species occupying several lattice sites, several models on one species
  (e.g. crowding plus activity coefficients), and cross-model shared sites.
- **Phase separation (Cahn–Hilliard).** A non-convex free energy (e.g. LFP's miscibility gap)
  makes c(ζ) multivalued. It needs c as an extra unknown plus a gradient-energy term, which is
  still local and block-tridiagonal.
- **Graded materials.** μ°(x), c_ref(x), ε(x), D(x) varying within a region, with the full
  `zφ̂ + μ°/RT − ln c_ref` change across each segment. φ's gauge freedom then becomes a
  continuous function, as in ESBD's inhomogeneities topic.
- **A sub-grid Gouy–Chapman interface law.** Diffuse layers treated analytically where λ_D ≪ h,
  for macroscopic devices with real double-layer charge on any grid. It tends to `neutral` as
  ε → 0.
- **Internal ports.** An outside phase attached at an interior node, exchanging species through
  the same laws as a contact (e.g. injecting salt mid-solution). It's local, because it only
  adds terms to that node's balance. Start with reservoir exchange; extra circuit terminals
  would come later.
- **Interface states** (charge that depends on the local μ̄_e, i.e. Fermi-level pinning) and
  thermionic-emission links at heterojunctions.
- Non-isothermal transport (Soret, Seebeck/Peltier, a heat equation), which stays local.

## Out of scope, for good

2D/3D; non-local couplings such as recirculation or one well-mixed reservoir feeding several
faces (a region with large mixing is the local stand-in); anything that breaks the
block-tridiagonal structure.
