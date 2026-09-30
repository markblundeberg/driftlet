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
- Validation suite (see the README), and docs.

## Next

1. **Statistics.** One general per-material interface: concentrations c(ζ) from reduced
   potentials, with their Jacobian, the chemical capacitance matrix. It must be the Hessian of
   a convex potential, so symmetric and positive definite, which is checked at construction.
   - Built-ins: ideal, lattice gas (shared-site groups: Bikerman-type crowding in double
     layers, site filling in intercalation hosts), Fermi–Dirac (degenerate carriers),
     regular solution / Redlich–Kister, extended Debye–Hückel, and tabulated monotone c(ζ)
     (e.g. measured OCV curves).
   - An escape hatch for user JavaScript functions. Such a device isn't serialisable, so it's
     built inside a Worker.
   - Fluxes via the excess-chemical-potential generalisation of Scharfetter–Gummel, which
     keeps equilibrium exact.
   - Validation: Fermi-function filling, the Kilic–Bazant–Ajdari crowded double layer, an OCV
     table round trip.
2. **Neutral-combination statistics** on ε = 0 regions. For e⁻ + Li⁺ in a host, the
   composition depends only on μ_Li = μ̄_Li⁺ + μ̄_e⁻ (the OCV curve), and the mutual chemical
   capacitance is the OCV slope. That's the natural description of battery electrodes.
3. **Time and frequency.**
   - Adaptive time steps from a BE/BDF2 error estimate, and BDF2 itself.
   - `advance(tEnd, { budgetMs })` for animation: as many steps as fit a frame budget.
   - Small-signal impedance, (J + iωM)·δx = b, with complex block Thomas, reporting Z(ω) and
     complex profiles (Warburg, double-layer and Maxwell–Wagner demos).
4. **Transport extras.** Imposed advection v(x), and current-free eddy mixing: an Onsager
   term −(D_mix/RT) P ∇μ̄, with P = C − C z zᵀ C/(zᵀ C z), that mixes composition without
   conducting and vanishes exactly at equilibrium.
5. **Performance pass.** Benchmarks tracked across versions. Cheaper warm re-solves (2–5
   Newton iterations for small parameter changes), allocation-free assembly, and a variable
   block size for absent species.
6. **Ship 0.1.0.** npm package and jsdelivr, live demos (code sandboxes, screenshots in the
   README), JSDoc types and a `.d.ts`, CI.

## Later

- **Cross-species transport coefficients together with cross chemical capacitances.** They're
  the two halves of one Onsager / Jamnik–Maier network, so one shouldn't come without the other.
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
