# driftlet: design spec

*A starting point for the agent building this library. Treat it as a strong draft: where it's
wrong or awkward, improve it and write down why.*

## 1. What it is

**driftlet** is a small, dependency-free JavaScript library that solves **1D drift–diffusion–
reaction problems**: Poisson–Nernst–Planck transport of any mix of charged species (ions,
electrons, holes), with bulk reactions and interfacial kinetics. It is fast enough to sit behind
a slider in a web page. Drag a bias or a concentration, and the solution updates at interactive
rates.

It covers electrochemical cells, membranes, electrolyte junctions, double layers, mixed
ionic–electronic conductors and simple semiconductor devices, all in one formulation.

**It is deliberately not** a TCAD or battery-modelling package. Its niche is teaching, live
demos and quick exploration. The README should point anyone needing more (2D/3D, thermal,
concentrated-solution transport, fitting) to serious tools: ChargeTransport.jl, Driftfusion,
SIMsalabim, PyBaMM, COMSOL.

**Why it should exist:** as of September 2026 there is no browser/npm drift-diffusion or PNP
solver. The existing codes are MATLAB, Julia, Pascal, Python or C++, and none of them can be
dropped into a web page.

## 2. Hard constraints

- **Pure JavaScript, ES module, zero runtime dependencies.** It must load from a CDN
  (`https://cdn.jsdelivr.net/npm/driftlet@x.y.z/+esm`) inside sandboxed pages that allow scripts
  only from a few CDNs and block `fetch`. So no WASM side files and no runtime fetches.
- **1D only, local physics only.** The Jacobian must stay **block-tridiagonal**. That is what
  keeps it fast and simple, and it's the whole point of an interactive toy. Every physics model
  must be local (a node couples only to its neighbours), and so is everything else: circuits
  are boundary conditions at one end, and conserved inventories come from the storage term
  (or, as a fallback, a running-integral unknown). There is no bordering and no dense row.
  Features that would break this structure are out of scope, not deferred.
- **Speed:** with ~300 nodes and up to 6 species, one Newton iteration should take well under
  1 ms in V8, and a warm-started re-solve after a small parameter change should take under
  ~5 ms. Inner loops should be allocation-free (`Float64Array`, reused buffers).
- **It must work inside a Web Worker** (no DOM access anywhere in the core).
- **Size budget:** aim for ≲ 30 kB minified.

## 3. Physics and conventions

Document all of this in the README; conventions are where users go wrong.

**Principles: honest thermodynamics.** These are the owner's firm requirements. The rest of
this section implements them.
- **Every control and every measurable is an electrochemical potential** (or a difference of
  them). A bias voltage sets a difference in μ̄_e⁻ between terminals (`V = −Δμ̄_e⁻/F`), and a
  gate sets its metal's μ̄_e⁻. Nothing a user sets is an electrostatic potential.
- **φ is bookkeeping: an arbitrarily anchored point on each material's standard-level ladder.**
  Any consistent φ gives identical physics (ESBD's "any consistent φ works equally well"). Each
  material's μ° values fix where its φ sits relative to its own V°_i ladder, so φ **jumps at
  every interface between distinct materials**, by that interface's dipole. Outputs and docs
  label φ as such.
- **Bulk standard chemical potentials are only meaningful in neutral combinations**
  (`μ°_e⁻ + μ°_h⁺`, `μ°_Na⁺ + μ°_Cl⁻`, `μ°_Li⁺ + μ°_e⁻`). How charged levels line up across a
  boundary between different materials is a property of *that interface* and must be given
  explicitly. No Anderson rule, no Schottky–Mott rule, no implied common vacuum level in the
  API. (The docs give those rules as explicit recipes users may apply by hand, with warnings;
  see §7.) Those rules imply that alignments add up transitively (so methanol | water |
  methanol would have cancelling dipoles). driftlet assumes no such thing: each interface's alignment is
  independent, even when the same pair of materials meets twice, since for example different
  adsorbates can sit on each interface. (A uniform alignment across each interface is an
  assumption only because the model is 1D; patch potentials and Schottky-barrier
  inhomogeneity are out of scope.)
- **A work function is a property of a surface, not a material.** A vacuum level
  (`φ_vac = V_e⁻ − W/e` just outside) is only reported at a free surface for which the user
  gave W, and it never feeds back into alignment.

**Species.** Each species *i* has a charge number `z` (any integer, **including 0**: neutral
mobile species such as dissolved gases or water are supported from v0), and, per material, a
diffusivity `D` and a **standard chemical potential** `μ°_i`. The primary variable is the
**electrochemical potential** `μ̄_i`. Ideal (Boltzmann/dilute) statistics by default:

- `μ̄_i = μ°_i(material) + z_i F φ + RT ln(c_i / c_ref,i)`, so
  `c_i = c_ref,i · exp( (μ̄_i − μ°_i − z_i F φ) / RT )`.
- The API speaks μ̄_i (J/mol, with eV helpers). The core works in `η_i = μ̄_i / RT`.
- Output conveniences for charged species: the species voltage `V_i = μ̄_i / (z_i F)` and the
  standard level `V°_i = φ + μ°_i / (z_i F)` (the band-edge analogue), as in the owner's
  [ESBD](https://marklundeberg.com/esbd/) diagrams. These are optional views, not the
  interface.
- Electrons and holes are ordinary species (`z = ∓1`), with `c_ref` = the effective density of
  states and `μ°` = the band edges (as molar energies). Plotted as voltages (V_i upward), the
  conduction band sits below the valence band; that's expected.
- **Lattice saturation (early optional add-on, per material).** Species sharing a site
  lattice of density c_max:
  `μ̄_i = μ°_i + z_i F φ + RT ln(c_i / (c_max − Σ_{j∈lattice} c_j))`. It inverts in closed form
  (`c_i = c_max·a_i / (1 + Σ_j a_j)`, with `a_i = exp((μ̄_i − μ°_i − z_i F φ)/RT)`), so it stays
  local and cheap in η. It matters wherever Boltzmann gives absurd numbers: double layers
  beyond a few RT/F (Bikerman-type crowding, instead of unbounded counter-ion pile-up) and
  intercalation hosts (Li in graphite or LFP, where site filling sets the OCV curve). A single
  species on its own lattice is Fermi–Dirac for one level.
  - Transport keeps the mobility form `N = −(D c/RT)∇μ̄`. Document that D is then a mobility
    coefficient, not the Fickian D, and offer the `(1 − θ)` hopping factor as an option.
  - Plain SG assumes Boltzmann. Use the excess-chemical-potential generalisation of SG (as in
    ChargeTransport.jl / VoronoiFVM), which still preserves equilibrium exactly.
- Degenerate band statistics (Fermi–Dirac integrals, e.g. Blakemore or Joyce–Dixon) for
  electrons and holes are a later add-on; they use the same generalised-SG machinery.

**Transport.** Particle flux per species:
`N_i = −(D_i c_i / RT) ∇μ̄_i + c_i v(x) + N_i^mix`, with continuity
`∂c_i/∂t + ∂N_i/∂x = R_i`.
- The first term is drift + diffusion in one (Einstein relation automatic). For charged species
  the charge current is `J_i = z_i F N_i = −σ_i ∇V_i`, with `σ_i = z_i² F² D_i c_i / RT`.
- `v(x)`: an optional **imposed advection velocity** (e.g. the axial flow at a rotating disk
  electrode, or flow through a porous layer). It's local, so it keeps the Jacobian
  block-tridiagonal, and SG handles drift + advection as one effective velocity.
- `N_i^mix`: optional **eddy mixing** with diffusivity `D_mix(x)`, large in a stirred bulk and
  falling to zero at walls. It models stirring locally and produces Nernst diffusion layers
  naturally. Stirring moves composition without conducting, so it must not be done by raising
  the D_i (that would also raise σ). It is a current-free Onsager term:
  `N_i^mix = −(D_mix/RT) Σ_j P_ij ∇μ̄_j`, with `P = C − C z zᵀ C / (zᵀ C z)` and `C = diag(c)`.
  - `P z = 0`, so it carries **no current, anywhere**, and ∇φ drops out (no migration).
  - It is driven by ∇μ̄, so it vanishes **exactly at equilibrium**, even inside double layers.
    The naive `−D_mix ∇c_i` breaks equilibrium there, and it also carries a spurious current
    wherever `Σ z_i c_i` varies (e.g. graded fixed charge).
  - P is symmetric positive semidefinite, so entropy production is ≥ 0.
  - In a neutral region with uniform fixed charge it reduces to `−D_mix ∇c_i` for every
    species: species-blind, like real advection. An inhomogeneous blob homogenises all
    species at the same rate, and mixing generates no diffusion potential of its own
    (molecular diffusion still does).
  - The zero-current projection isn't unique. The c-weighting is the species-blind
    (advective) choice; a D·c weighting would instead reproduce ambipolar molecular
    diffusion, which is wrong for stirring.
  - Neutral species are simply mixed: `−D_mix ∇c` within a region.
  Discretely, use the logarithmic mean of c at faces. Then `c_face Δln c = Δc` exactly, and
  equilibrium is preserved exactly for any face values. It couples species within a face
  (dense within-block entries), which the block structure already allows. Default 0.
- Out of scope, permanently (non-local): recirculation loops, or one well-mixed reservoir
  feeding several faces. The local stand-in for a well-mixed volume is a region with large
  D_mix.

**Electrostatics.** `−∂/∂x (ε ∂φ/∂x) = F Σ z_i c_i + ρ_fixed(x)`, where `ε` belongs to the
material and the fixed charge `ρ_fixed` (doping, ionomer charge, …) to the region.

**Materials, regions and interfaces.** A *material* holds bulk properties: ε, and per species
D, c_ref and μ°_i (or "absent"). Its μ°_i are defined only up to a charge gauge
`μ°_i → μ°_i + z_i F s`, since only neutral combinations are bulk-measurable, so a material's
absolute charged levels mean nothing on their own. A *region* is a material plus a length and
a fixed charge ρ_fixed (doping, ionomer). Every face between regions of *different* materials
takes exactly one alignment number of its own. Given both materials' μ° values, this is
equivalent to the φ jump (dipole) at that face. It is one of:
- the step in the standard level of a charged species present on both sides (e.g. ΔE_c for a
  heterojunction; ΔE_v then follows from the two bulk gaps);
- the standard free energy of a charge-transfer reaction across the face (for an electrode,
  the standard potential E° *is* the alignment);
- the dipole (φ jump) itself, in the materials' own φ anchoring.
Faces between regions of the same material default to no dipole (a homojunction, e.g. pn),
and can take one if wanted (e.g. a grain boundary). Alignments are per face and independent:
nothing requires them to be consistent around the device, and the same material pair can
meet twice with different alignments. The overall offset of all μ̄ is left free and fixed by
a terminal. A missing or doubled alignment at a heterointerface is a construction error, never
a silent default. Faces may also carry per-species links (as for contacts, below: e.g. thermionic
emission or an interface resistance instead of local equilibrium) and a fixed sheet charge.
Species that can't enter a region are blocked at its face (zero flux).

**Absent species.** Where a species is absent, its η_i is undefined and has no equation. Two
implementation options: keep a fixed block size with decoupled identity rows, or let the block
size vary per region (block Thomas copes with non-square off-diagonal blocks). Start with
whichever is simpler and benchmark the other; the choice must not leak into the API. Outputs
for an absent species are `NaN` at those nodes, so plots show a gap rather than a fake value.

**Conserved inventories (blocked species).** A species blocked on every side of some connected
stretch of the device (no contact link, no reaction) is a *conserved spectator*: its amount
there never changes, and it just drifts and screens. Transients conserve it automatically, and
so does the giant-time-step steady solve (§4). The pure dt = ∞ steady equations can't determine
it: its μ̄ level floats, and the balance equations are rank-deficient by one. Fallback, if the
giant step proves numerically fragile: make `∫c dx = amount` local with a running integral.
Over that stretch the species gets a second per-node unknown `Q_k = Σ_{j≤k} c_j·vol_j`, with
rows `Q_k − Q_{k−1} − c_k·vol_k = 0` (starting from zero). At the stretch's last node, the
redundant balance row is replaced by `Q = amount`, a condition at the far end. The amount is
whatever the current state holds, or is given explicitly. With bulk reactions, the conserved quantities
are moieties (left null vectors of the stoichiometry restricted to that stretch) rather than
species. A device can be *entirely* blocked (a floating island) provided φ is anchored
electrostatically somewhere (a capacitive gate link, below); then the island simply holds a fixed charge.
With no electrostatic anchor at all, the problem is singular and must be rejected up front.

**Bulk reactions (mass action, thermodynamically consistent).** Each reaction has a
stoichiometry over mobile species plus optional fixed-activity neutral participants (e.g. H₂O).
Rate: `r = k_f Π c_reactants − k_b Π c_products`. Derive `k_b` from `k_f` and the standard
chemical potentials, so that equilibrium is *exactly* the point where the affinity
`A = −Σ ν_i μ̄_i` (over all participants, charged or neutral) vanishes. Canonical cases to support in v0: e⁻ + h⁺ ⇌ ∅
(generation/recombination) and H⁺ + OH⁻ ⇌ H₂O.

**Interfacial reactions (electrode kinetics).** A reaction located at a region boundary, e.g.
Mⁿ⁺(electrolyte) + n e⁻(metal) ⇌ M(s). Use Butler–Volmer written in the reaction's
overpotential, `η = affinity / (nF)`:
`i = i0 · [ exp(α n F η/RT) − exp(−(1−α) n F η/RT) ]`,
with optional concentration dependence of `i0`. The rate enters as a boundary flux for each
participant.

**Contacts (domain ends) as links.** Each end couples every species and φ to the outside
world through a *link*. This is the model from the owner's earlier linear prototype, and it is
more general than a zoo of contact types:
- **species links**, one per species: *fixed* (μ̄_i pinned to an outside value: the normal
  case), *conductance* G_i to an outside μ̄_i (an ohmic interface resistance), *kinetic*
  (Butler–Volmer or thermionic emission, i.e. a nonlinear conductance), or *blocked* (G_i = 0);
- **electrostatic link**: *neutral* (the boundary node is locally neutral given its fixed μ̄_i:
  an ohmic contact or a well-stirred bath, where no alignment is needed because the outside
  composition is given), *capacitive* (a capacitance per area C to an outside conductor whose
  μ̄_e⁻ is set, e.g. a gate across an insulator or a metal across its Stern layer), or *free*
  (zero field). A capacitive link needs its own zero-charge alignment, an interface property:
  the flat-band voltage for a gate, the potential of zero charge for a metal/electrolyte
  interface. It is never computed from a work function and an electron affinity. C → ∞ is the
  ideal limit: a fixed μ̄_e⁻ plus the alignment, which is all a "fixed φ" boundary ever
  honestly meant. The gate "sets" φ only through this alignment; users never set φ directly.
Familiar contacts are presets: bath = all species fixed + neutral; semiconductor ohmic = e⁻/h⁺
fixed + neutral; Schottky = e⁻/h⁺ thermionic with a given barrier height + capacitive; metal
electrode = e⁻ terminal, ions blocked or reacting (alignment E°), capacitive Stern layer
(alignment pzc); ideal insulator = all blocked + capacitive gate. Infinite G or C must be
imposed exactly (eliminate the unknown), not with the large-penalty trick the prototype used.

**Contact values and the circuit (as implemented).** Each contact has a terminal voltage V
(default 0). A fixed charged species there sits at `V_i = V + offset_i`, i.e.
`μ̄_i = z_i F (V + offset_i)`. The offset is an interface property: 0 for the contact's named
terminal species (e⁻ at a metal; h⁺ at the same metal is also 0, since V_h = V_e when
e⁻ + h⁺ ⇌ ∅ is in equilibrium there), E° for an ion at a reversible electrode. It is required
for every fixed species except the terminal one. A fixed neutral species takes an absolute μ̄.
Ohmic contacts are simply fixed V_e⁻ and/or V_h⁺: both gives an infinite-recombination contact
(np = n_i² there), one alone a selective contact. A **bath** is a helper, not a link type: from
a neutral composition and a charged reference species (Cl⁻ for Ag/AgCl) it computes every
species' offset (neutral species: μ̄), with a neutral φ link. It passes current through its
fixed species; open-circuit use needs the current-mode circuit.
- Numerically, a fixed link replaces the boundary node's balance row with a Dirichlet row, and a
  neutral φ link replaces its Poisson row with local neutrality (total charge, mobile plus
  fixed, zero). Before replacement, each row's residual *is* the flux (or displacement)
  through the contact, so contact fluxes and terminal currents are recovered exactly, with no
  extra unknowns. (A left-end "contact flux node" would have a singular first block for block
  Thomas, so that route was rejected.)
- A contact with any connected species must state its φ condition (neutral, capacitive,
  free). A capacitive link's gate voltage is the contact's V.
- Circuit modes (implemented): `voltage` (each contact at its V), `current` (fixed I through
  the right terminal; I = 0 is open circuit), `load` (I = (V_right − V_left − V_src)/R). In
  the floating modes, the right terminal voltage is read off the terminal species' own η at
  the contact node. Its row becomes `I_segment − I_circuit(V_t) = 0`, where I_segment is the
  total current (conduction plus displacement) through the last grid segment. That equals
  the terminal current exactly, since the box balance plus Poisson, differenced in time,
  telescope. Other fixed species there follow V_t. No extra unknowns are needed; this is
  the owner's "Robin condition at the rightmost node".

**Temperature:** a single uniform T in v0 (a parameter, default 298.15 K). Non-isothermal
transport is a possible later extension and does *not* require non-ideal solutions: it needs
∇(1/T) as an extra Onsager force with heats of transport (Soret, Seebeck/Peltier) plus a heat
equation (conduction, Joule and Peltier heating), with T as one more per-node unknown. Those
couplings are all local, so the structure stays block-tridiagonal, but it's many more
user-supplied coefficients.

**Units:** SI throughout (m, mol/m³, J/mol, V, A/m², F/m, s). Provide small helpers for
mol/L, µm/nm, and eV ↔ J/mol. **Current sign:** positive current = positive
charge flowing toward +x.

## 4. Numerics

- **Grid:** 1D non-uniform node grid. Every material interface is a *doubled node* (below), so
  each grid segment lies inside one material. Include a grid
  builder with geometric/tanh grading toward interfaces and contacts, since Debye lengths
  (nm) are much smaller than devices (µm–mm). Given a target resolution near interfaces, it
  should pick spacing automatically.
- **Discretization:** a finite-volume (box) method with **Scharfetter–Gummel** fluxes (the same
  exponential fitting as Il'in / Allen–Southwell for convection–diffusion). For species *i* the
  "drift" across a face combines `z_i F Δφ + Δμ°_i` (over RT) with the advective Péclet number;
  for neutral species with no flow it reduces to plain diffusion. Because segments never
  straddle an interface, D, ε and μ° are constant along every SG segment; steps only happen at
  doubled nodes. Evaluate the
  Bernoulli function `B(x) = x/(eˣ−1)` stably (series near 0, asymptotics for large |x|).
- **Sharp interfaces: doubled nodes plus a flux node.** In the grid, each interface is two
  nodes at the same x, one per side. In the linear system, a zero-volume *flux node* sits
  between them. Its unknowns are the interface fluxes (one per species, plus the displacement
  for φ), and its rows are the interface laws. For each variable:
  - *fixed offset* (the default: local equilibrium for species, i.e. offset 0; for φ, the
    interface dipole): the row is `u_R − u_L = offset`;
  - *linked*: the row is `flux − f(u_L, u_R) = 0` (Butler–Volmer, thermionic emission,
    interface resistance; for φ, an interface capacitance `D − C·Δφ`);
  - *blocked* or absent on one side: `flux = 0`.
  The left node's balance row subtracts the interface flux and the right node's adds it. Each
  row then couples only neighbours, with a fixed block size and no penalty terms.
  (Correction to an earlier draft, which put the flux in the right copy's slot and read its
  value from the left node. That makes the right node's other segment reach two nodes back,
  which breaks the tridiagonal structure.) A continuity row has a zero diagonal block, but
  block Thomas still sees a non-singular Schur complement there, since within-block pivoting
  and the preceding elimination supply it. Cost: two extra solver nodes per interface. The
  interface fluxes come out as unknowns, which is exactly what `sol.interfaces` and electrode
  currents need. Fixed sheet charge and, later, interface states enter the adjacent balance
  rows.
- **Precision of small currents (found in the pn junction).** A majority carrier carrying a
  small current has a quasi-Fermi step between nodes far below the ulp of η (≈3e-19 against
  3.5e-15), and in the concentration form the same current is a cancellation of two huge
  terms. Two fixes, both in place:
  - SG is evaluated as `N = −(D/h)·B(Δ)·c_L·expm1(η_R − η_L)`, algebraically identical but
    precise relative to Δη, and exactly zero at equilibrium;
  - η (with φ̂) is stored as a compensated double-double (hi + lo). Only differences use the
    lo part, updates use an error-free two-sum, and the Jacobian and linear solve stay in
    plain doubles.
  Result: left and right contact currents agree to ~1e-10 relative down to µA/m² leakage,
  where they previously disagreed by 3%. This is the same problem commercial TCAD tools
  address with extended precision.
- **Unknowns (per node):** the dimensionless potential `φ̂ = Fφ/RT` and the dimensionless
  electrochemical potentials `η_i = μ̄_i / RT`. Scale everything internally (thermal energy,
  reference concentration, Debye length). This is a trade-off, not a free win:
  - *For η:* positivity of c is automatic, so there's no clamping, which would itself break
    conservation. Flat-μ̄ equilibrium is exact and trivial. The physics is posed in μ̄ anyway:
    fixed-μ̄ links, interface offsets, affinities and Butler–Volmer are all linear or natural in
    η.
  - *Against η:* the storage term is nonlinear, so conservation holds only to the Newton
    residual (see Conservation). A species' rows scale with its c, so where it is scarce
    (minority carriers deep in a depletion zone) its η is weakly determined and Newton can
    swing it wildly. That is harmless physically but can stall damping. So the damping limit
    must not be dominated by δη at nodes where that species is negligible: limit δφ̂ strictly,
    and weight δη limits by relevance. Also row-scale the Jacobian.
  - *(φ, c) variables,* the TCAD default, give linear storage (exact conservation every
    iteration) and polynomial reactions. But Newton overshoots to negative c, and the fixes
    (damping, clamping) are fiddly or non-conservative. Slotboom variables `u = exp(η)`
    linearise SG but overflow a double at modest voltages. Stay with η unless benchmarks say
    otherwise.
- **Solver:** damped Newton, limiting each update to a few thermal voltages per iteration (or
  Bank–Rose damping). The Jacobian is block-tridiagonal with block size `1 + nSpecies`; solve it
  with block Thomas (a small dense LU per block). Nothing is bordered (§2): circuits are
  boundary rows and inventories ride on the storage term, so the solver is plain block Thomas.
- **Transient first; steady state is a giant time step.** This is the owner's linear
  prototype's trick, generalised. A backward-Euler step with dt far beyond the slowest time
  constant *is* the steady-state problem, except that the storage term `(c − c_old)·vol/dt`
  stays in. That keeps the Jacobian nonsingular, and it conserves every spectator inventory
  exactly: summing a spectator's rows telescopes the fluxes away and leaves
  `Σ (c − c_old)·vol = 0` for any dt. `solve()` is therefore one or two backward-Euler steps at
  a huge dt (repeat until nothing changes; at that fixed point the steady equations hold
  exactly). Nonlinearity is handled by Newton on each step. When Newton struggles (first
  solves, big jumps), fall back to pseudo-transient continuation: start at a modest dt and
  grow it as the residual falls (switched evolution relaxation, Kelley–Keyes Ψtc). That
  follows a physical path, so it's robust, and it becomes plain Newton as dt → ∞.
  - Initial guess for a cold start: local neutrality region by region. The equilibrium problem
    (flat μ̄, nonlinear Poisson–Boltzmann) is convex and converges reliably with step limiting.
  - `solve({ warm: previous })` reuses the previous solution; that is what makes slider
    interaction fast (typically 2–5 Newton iterations in all).
  - Risk to test early: the spectator level is pinned only by the tiny storage term, so with
    dt/τ_fastest ≳ 1/ε_machine, round-off could let the inventory drift (stiff devices have
    τ_slow/τ_fast ~ 1e12). Pick dt relative to τ_slowest, not arbitrarily huge. Check the
    floating-island inventory after a giant step. If it drifts, use the running-integral
    fallback (§3).
- **Transient and time-step control.** Backward Euler with a BDF2 option, finite-volume mass
  matrix, and step halving when Newton fails. Real devices span nanoseconds (dielectric
  relaxation) to kiloseconds (diffusion across a millimetre), so dt must adapt:
  - estimate local error from the BE–BDF2 difference (or step doubling) and grow or shrink dt
    to hit a tolerance;
  - log-spaced output times for "watch it relax" views;
  - an animation-friendly call, `advance(tEnd, { budgetMs })`, that does as many internal
    steps as fit in the budget and returns the state reached, so a UI can call it once per
    frame and stay responsive.
- **Conservation: exact equations, iterate only as good as Newton.**
  - *Discretely exact.* Summed over boxes, each face flux cancels exactly in floating point,
    provided it is **computed once per face and added with ± to both neighbours**, never
    recomputed per node. The same goes for the interface flux unknowns at doubled nodes.
    Reaction sources are written as `ν_i · r` with one r per reaction per node, so every
    moiety (mᵀν = 0) cancels exactly too. BE and BDF2 are both conservative multistep forms.
    So the *converged* discrete solution conserves every inventory, and total charge, to
    round-off.
  - *The iterate conserves only to the Newton residual.* With c = exp(η…) nonlinear, a Newton
    update conserves the linearised inventory exactly, but the true one only to second order
    in the update. Per step, the drift equals the summed residual of that species' rows
    (times dt). This is the one real difference from the linear prototype. It's harmless if
    the stopping test is stated in conservation units: stop when the summed per-species
    residual implies an inventory drift below ~1e-13 of that inventory per step, not just a
    generic norm. Quadratic convergence makes this cost about one extra iteration.
  - *Optional exact fix for spectators:* after each accepted step, shift a spectator's η level
    uniformly by `ln(amount_target / amount_now)` (a ~1e-13 nudge), which makes its inventory
    exact to round-off whatever Newton did. Measure first; add it only if drift shows up in
    long runs.
  - *Diagnostics:* every solution reports conservation bookkeeping: per species (or moiety),
    the mismatch between the change in inventory and the time-integrated boundary flux, and
    likewise for charge. "Never silently wrong" applies to conservation too.
- **Small-signal impedance.** Linearise about a steady state and solve
  `(J + iωM)·δx = b` at each frequency: the same block-tridiagonal Jacobian J and mass matrix
  M, in complex arithmetic (or real blocks of doubled size). It's local, so plain block Thomas
  applies; one factorisation per frequency. The output is the terminal impedance Z(ω) plus
  the complex profiles δμ̄_i(x), δφ(x), for Warburg, double-layer and Maxwell–Wagner demos.
- **Grid under parameter changes.** Warm starts need the same grid, so the grid stays fixed
  across `set()` unless the geometry changes (then regrid and interpolate the warm state).
  After every solve, check resolution against the local Debye length (and diffusion-layer
  width). If a parameter change has left a layer under-resolved, say so in the solution's
  diagnostics rather than returning a quietly wrong profile.
- **Quasi-neutral mode (after v0):** replace Poisson by local neutrality
  `F Σ z_i c_i + ρ_fixed = 0` at each node. This is much cheaper and has no Debye layers; φ then
  jumps at interfaces (Donnan, junction potentials): at doubled nodes the two φ copies are
  simply independent, each fixed by neutrality on its own side.
- **Interface alignment and unresolved double layers:** when the Debye length is far below the
  grid spacing (metals, concentrated electrolytes), the double layer collapses into the
  doubled-node pair. The alignment then only sets that sub-grid charge, and it drops out of
  the μ̄ profiles automatically, as it should. Check this in tests (metal | metal, and a
  quasi-neutral limit).
- **φ in outputs** is exactly the φ solved for: anchored per material by the user's μ° values,
  with a jump at each doubled node given by its alignment (plus any interface-capacitance
  charging). The size of a dipole is as subjective as the anchoring, correctly so; a user who
  has data for true mean inner potentials can anchor to those. Output arrays list each doubled
  node twice (x repeated), so c, μ° and φ plot as true vertical steps.
- **Failure behaviour:** never return a silently wrong answer. Report non-convergence with the
  residual history, and expose `converged`, `iterations` and `residual`.

## 5. API sketch (to be refined)

```js
import { Device } from 'driftlet';

const dev = new Device({
  T: 298.15,
  species: [
    { name: 'Zn2+', z: 2, cRef: 1000 },       // cRef in mol/m³ (1 M)
    { name: 'SO4 2-', z: -2, cRef: 1000 },
  ],
  materials: {
    water: { eps: 78.5 * 8.854e-12,            // μ°s on the usual SHE-based table convention
      species: { 'Zn2+': { D: 0.70e-9, mu0: -147.1e3 }, 'SO4 2-': { D: 1.07e-9, mu0: -744.5e3 } } },
  },
  regions: [ { name: 'electrolyte', material: 'water', length: 1e-3, fixedCharge: 0 } ],
  interfaces: [ /* per face between different materials: exactly one alignment, optional links */ ],
  bulkReactions: [ /* stoichiometry, kf, neutral activities */ ],
  contacts: {                                  // links (§3); presets expand to these
    left: { preset: 'metalElectrode', terminal: 'e-',
            reaction: { eq: 'Zn2+ + 2e- = Zn(s)', E0: -0.762, i0: 10, alpha: 0.5 }, // E0 = alignment
            stern: { C: 0.2, pzc: -0.63 } },   // F/m²; pzc = zero-charge alignment, V (same convention)
    right: { /* same */ },
  },
  circuit: { mode: 'load', R: 50 },            // or { mode: 'voltage', V } / { mode: 'current', I }
  grid: { refineNear: 'interfaces', minSpacing: 1e-9, maxSpacing: 5e-6 },
});

let sol = dev.solve();                          // steady state
dev.set({ circuit: { mode: 'load', R: 20 } });
sol = dev.solve({ warm: sol });                 // fast re-solve

// Plotting-ready output:
sol.x;                                           // Float64Arrays
sol.mu['Zn2+']; sol.c['Zn2+']; sol.N['Zn2+'];     // μ̄ (J/mol), c, particle flux
sol.muStd['Zn2+'];                                // standard level μ° + zFφ (band-edge analogue)
sol.phi;                                          // bookkeeping φ (per-material anchor; jumps at interfaces)
sol.V['Zn2+']; sol.Vstd['Zn2+']; sol.J['Zn2+'];    // optional voltage-scaled views
sol.current; sol.terminalVoltage;
sol.interfaces;   // per face: steps in each μ̄_i and V_i (overpotentials, junction/Donnan potentials)
sol.converged; sol.iterations; sol.residual;
sol.conservation;  // per species/moiety and charge: Δinventory vs ∫boundary flux dt

dev.step(dt, { from: sol });                     // one transient step
dev.advance(tEnd, { from: sol, budgetMs: 8 });   // adaptive steps within a frame budget
dev.impedance(freqs, { about: sol });            // small-signal Z(ω) and profiles
```

Design goals for the API: plain objects in, typed arrays out, no classes users must subclass.
Everything should be serializable, so a device definition can be posted to a Worker.

## 6. Validation suite (the real deliverable)

Each item is an automated test with a stated tolerance. The README carries a table of them.

1. **Equilibrium invariance:** no bias and no imposed advection v ⇒ every μ̄_i flat to ≤1e-9·RT and
   all fluxes zero, for every example device. **Gauge invariance:** shifting one material's
   μ°_i by z_i F s (alignments given as relative steps, so unchanged) leaves every μ̄_i, c_i,
   flux and current identical to round-off; only the bookkeeping φ moves.
2. **Nonlinear Poisson–Boltzmann / Gouy–Chapman:** diffuse layer at a charged blocking wall in a
   1:1 electrolyte; the potential profile and differential capacitance match the analytic
   formulas. Its linear limit gives Debye screening.
3. **Donnan potential** across a fixed-charge membrane between two baths vs the analytic value.
4. **Floating island:** every species blocked, φ anchored only by capacitive gate links at both ends.
   Each inventory is conserved to round-off through a gate sweep; Gauss's law holds exactly
   (island charge = −(total charge on the gates)); in the linear limit the island's response matches
   the series network C_g + C_diffuse (Debye) at each face. Also the singular case (no
   electrostatic anchor) is rejected with a clear error.
5. **Planck liquid junction:** steady diffusion junction of a binary salt between two
   reservoirs. The diffusion potential is `(RT/F)(t₊ − t₋)·ln(c₁/c₂)` for a 1:1 salt.
6. **Concentration polarization / limiting current:** a symmetric metal | binary salt | metal
   cell with deposition/dissolution. Match the analytic quasi-neutral profile, the I–V
   relation and the limiting current `i_lim = z F D_salt c₀ · 2 / ((1 − t₊) L)` (derive and
   check the exact form in the test).
7. **Butler–Volmer interface:** a single electrode facing a reservoir; the current–overpotential
   curve is exact.
8. **pn junction:** built-in potential `(kT/q) ln(N_A N_D / n_i²)` exactly; depletion width
   within the depletion approximation's accuracy; long-diode Shockley J–V within a few % at low
   injection. The J–V part needs recombination (e⁻ + h⁺ ⇌ ∅, M2). Mass action gives the
   radiative-like `k(np − n_i²)`, so the minority lifetime is `1/(k·N_A)` at low injection; it is
   not SRH, and the README should say so.
9. **Heterojunction:** an abrupt pn heterojunction with a user-given ΔE_c. The built-in
   potential and the split of band bending between the sides match the depletion approximation
   for *that* ΔE_c, and omitting the alignment is a construction error (no Anderson default).
   Also an A | B | A stack (methanol | water | methanol, say) with *unequal* alignments at its
   two faces: equilibrium holds, and each face's double-layer charge and potential split match
   its own analytic value.
10. **Bulk reaction:** the mass-action law holds at equilibrium (e.g. `c_H+ c_OH− = K_w`), and
   homogeneous relaxation matches the analytic time constant.
11. **Advection–diffusion:** Levich limiting current at a rotating disk electrode (imposed
   axial velocity profile), and a Nernst layer from an eddy-diffusivity profile. Mixing carries
   zero current to round-off everywhere (including with graded fixed charge), leaves bulk
   conductivity unchanged (the same IR drop with and without mixing at uniform composition),
   and leaves equilibrium double layers untouched (D_mix nonzero right up to a charged wall,
   test 1 still passes). A non-uniform neutral blob homogenises at rate D_mix for all species.
12. **Transient:** a blocking electrolyte cell charging with the analytic time constant
   (≈ λ_D·L/D in the small-signal limit). Conservation: per-step drift ≤ 1e-13 relative, and
   ≤ 1e-9 relative accumulated over a 1e4-step run, for every inventory and for total charge;
   the reported diagnostics agree with an independent recount.
13. **Grid convergence:** observed order ≈ 2 on smooth problems.
14. **Cross-check (optional, documented):** reproduce one published example from
    ChargeTransport.jl or Driftfusion within stated tolerance.
15. **Performance benchmark:** time per Newton iteration and per warm re-solve on reference
    devices, tracked across versions.
16. **Lattice saturation:** a single-species lattice gas reproduces the Fermi-function
    filling curve exactly. A Bikerman-type double layer matches the analytic
    Kilic–Bazant–Ajdari profile and capacitance (camel shape, then saturation at large drop).
    Equilibrium invariance (test 1) still holds with generalised SG.
17. **Impedance:** a blocking electrolyte cell matches its analytic R–C response. A reversible
    electrode with finite-length diffusion matches the finite Warburg element, including the
    45° region and both low-frequency limits. The ω → 0 impedance matches the slope of the
    steady I–V curve.

## 7. Deliverables and packaging

- **Repo layout:** `src/` (core: grid, SG fluxes, assembly, block-tridiagonal solver, Newton,
  transient); `test/` (node's built-in test runner, no deps); `bench/`; `examples/` (standalone
  HTML demos importing from jsdelivr, built as they become relevant; no fixed list).
- **Docs:** README (what it is and isn't, conventions, API, validation table, pointers to
  heavier tools), a ROADMAP, and JSDoc types (optionally a `.d.ts`). This design brief is
  temporary: it dissolves into README, ROADMAP, docs and code. Code, tests and docs never
  cite it.
- **Alignment guide (docs page):** explains per-interface alignment and bookkeeping φ, then
  gives precise, worked recipes for turning vacuum- or reference-based data into driftlet
  alignments, for users who have nothing better (it's a legitimate best guess) and as a
  reference example. Recipes: Anderson (ΔE_c from electron affinities); Schottky–Mott
  (barrier height from W_metal − χ); a gate's flat-band voltage from W_gate − W_semiconductor;
  pzc via the "absolute" SHE (4.44 V) and a metal's work function; single-ion transfer
  energies between solvents via an extrathermodynamic assumption (e.g. TATB). Each recipe
  states its hidden assumptions (vacuum levels line up, alignments add up transitively and are
  laterally uniform) and cites where it is known to fail (covalent semiconductors, Fermi-level
  pinning). These stay **recipes in the docs, never API**: no helper computes an alignment from
  vacuum quantities.
- **npm:** ESM with an `exports` map, `sideEffects: false`, semver starting at 0.1.0, and a
  CHANGELOG.
- **CI:** GitHub Actions running tests and the benchmark on push.
- **Licence:** 0BSD (public-domain-equivalent; SPDX `0BSD`).

## 8. Milestones

Each milestone ends green: its tests pass, and the benchmark has no regressions.

- **M0: skeleton and solver core.** package.json (ESM, exports map), `node --test` setup,
  stable Bernoulli function, block Thomas (tested against a dense solve, including
  near-singular blocks), grid builder with grading and doubled interface nodes, and device
  construction with validation (missing or doubled alignment, no electrostatic anchor,
  unknown species: all clear errors). No physics yet; everything unit-tested.
- **M1: equilibrium through the time-step machinery.** Assembly (Poisson, SG fluxes, storage;
  each face flux computed once), end links (fixed, blocked, neutral, capacitive, free),
  doubled-node fixed-offset interfaces, backward Euler + damped Newton, giant-step `solve()`
  with Ψtc fallback, warm starts, conservation diagnostics, and the inventory round-off check.
  Tests 1, 2, 3, 4, 8 (built-in potential), 9.
- **M2: transport, circuits, first transients.** Bias, galvanostatic, load-resistor and
  open-circuit boundary conditions; conductance links; steady transport; plain fixed-dt
  transients. Tests 5, 6, 12, 13, and a first benchmark baseline (15), so performance is
  tracked from here on.
- **M3: reactions and statistics.** Bulk mass action (moiety-exact sources), Butler–Volmer and
  thermionic links, lattice saturation with generalised SG. Tests 7, 8 (J–V), 10, 16.
- **M4: time and frequency.** Adaptive dt, BDF2, `advance()` with a frame budget, small-signal
  impedance, advection and current-free mixing, grid-resolution diagnostics. Tests 11, 17.
- **M5: ship 0.1.0.** README, ROADMAP, alignment guide, JSDoc/.d.ts, CI, whatever examples
  exist by then, npm publish.
- **Later (ROADMAP):** quasi-neutral mode, degenerate band statistics, interface states
  (charge that depends on μ̄_e⁻, i.e. Fermi-level pinning), concentration-dependent D,
  non-isothermal transport, optional cross-check (14).

## 9. Questions to settle with the owner early

Settled (2026-09-29): primary branch is `master`; honest-thermodynamics principles, with the
API in μ̄ and V_i as optional views (§3); blocked species are conserved spectators and
all-blocked islands are valid with gate anchoring (§3); local physics only (§2); demos are
built as they become relevant.

- npm: which account owns the package (personal, or an npm org for co-maintainers), whether
  to reserve the name `driftlet` now (free as of 2026-09-29) with a placeholder publish, and
  whether releases are published by hand or from CI with a token.
