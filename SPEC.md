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
  (or, as a fallback, a running-integral unknown). There is no bordering and no dense row. Features that would break this structure
  are out of scope, not deferred.
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
  API. (The docs give those rules as explicit recipes users may apply by hand, with warnings; see
  §7.) Those
  rules imply that alignments add up transitively (so methanol | water | methanol would have
  cancelling dipoles). driftlet assumes no such thing: each interface's alignment is
  independent, even when the same pair of materials meets twice, since for example different
  adsorbates can sit on each interface. (A uniform alignment across each interface is an
  assumption only because the model is 1D; patch potentials and Schottky-barrier
  inhomogeneity are out of scope.)
- **A work function is a property of a surface, not a material.** A vacuum level
  (`φ_vac = V_e⁻ − W/e` just outside) is only reported at a free surface for which the user
  gave W, and it never feeds back into alignment.

**Species.** Each species *i* has a charge number `z` (any integer, **including 0**: neutral
mobile species such as dissolved gases or water are supported from v0), a diffusivity `D` (per
region), and, per region, a **standard chemical potential** `μ°_i`. The primary variable is the
**electrochemical potential** `μ̄_i`. Ideal (Boltzmann/dilute) statistics:

- `μ̄_i = μ°_i(region) + z_i F φ + RT ln(c_i / c_ref,i)`, so
  `c_i = c_ref,i · exp( (μ̄_i − μ°_i − z_i F φ) / RT )`.
- The API speaks μ̄_i (J/mol, with eV helpers). The core works in `η_i = μ̄_i / RT`.
- Output conveniences for charged species: the species voltage `V_i = μ̄_i / (z_i F)` and the
  standard level `V°_i = φ + μ°_i / (z_i F)` (the band-edge analogue), as in the owner's
  [ESBD](https://marklundeberg.com/esbd/) diagrams. These are optional views, not the
  interface.
- Electrons and holes are ordinary species (`z = ∓1`), with `c_ref` = the effective density of
  states and `μ°` = the band edges (as molar energies). Plotted as voltages (V_i upward), the
  conduction band sits below the valence band; that's expected.
- Degenerate (Fermi–Dirac) statistics are a later, optional add-on (e.g. the Blakemore
  approximation).

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

**Electrostatics.** `−∂/∂x (ε ∂φ/∂x) = F Σ z_i c_i + ρ_fixed(x)`, where `ε` and the fixed
charge `ρ_fixed` (doping, ionomer charge, …) are set per region.

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

**Circuit.** A terminal is the μ̄ of a named species at a contact: μ̄_e⁻ of the metal in the
normal case. For an ion-only cell, an ion (e.g. Cl⁻) stands in for an ideal reversible
electrode for that ion (Ag/AgCl). Terminal voltage is `−Δμ̄_e⁻/F`, or `Δμ̄_i/(z_i F)` for an ion
terminal. The left terminal is ground (fixed μ̄), and the circuit is a boundary condition on the
rightmost node, so it stays local:
- applied voltage: Dirichlet (fixed μ̄);
- applied current (galvanostatic): a Neumann flux;
- a load resistor R: a Robin condition relating the terminal μ̄ to the terminal current.
Open circuit is the current mode with I = 0 (liquid junctions, Donnan). This is exact in
transients too, because in 1D the total current (conduction plus displacement) is uniform in
x: what leaves one terminal enters the other. When several species share a terminal (e⁻ and h⁺
at an ohmic contact), the terminal node's rows use their summed current. Their μ̄ are tied
together, so only one terminal level is unknown.

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
- **Sharp interfaces are doubled nodes:** two nodes at the same x, one per side, joined by a
  zero-width face. The zero-width face's flux law is the interface link, not SG. For each
  variable (each η_i, and φ) the pair is one of:
  - *linked*: both copies are unknowns and the face flux is a function of them: a kinetic or
    conductance link for species (Butler–Volmer, thermionic emission, interface resistance), or
    an interface capacitance for φ (displacement = C·Δφ). Blocked means zero flux.
  - *fixed offset* (the default: local equilibrium for species, i.e. offset 0; for φ, the
    interface dipole): the right copy's slot holds the **interface flux** (for φ, the
    displacement) as its unknown instead of a duplicate value. The left row is left-box balance
    minus that flux. The right row is that flux minus right-box balance, reading its own value
    as the left node's value plus the fixed offset.
    This couples only neighbours, so it stays block-tridiagonal with a fixed block size and
    no penalty terms. Summing the two rows gives exactly the classic shared-node box method,
    each half-box using its own material's μ°, ε and D.
  The interface fluxes then come out as unknowns, which is exactly what `sol.interfaces`
  (and electrode currents) need. Fixed sheet charge and, later, interface states enter the
  φ balance at the pair. Where a species is present on one side only, its copy on the other
  side is an absent-species row. Cost is one extra node per interface.
- **Unknowns (per node):** the dimensionless potential `φ̂ = Fφ/RT` and the dimensionless
  electrochemical potentials `η_i = μ̄_i / RT`. Concentrations can span 40 orders of magnitude
  (minority carriers), and these log-like variables keep them well conditioned. At equilibrium
  every `μ̄_i` is exactly flat. Scale everything internally (thermal energy, reference
  concentration, Debye length).
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
- **Transient:** backward Euler with a BDF2 option, finite-volume mass matrix, and step halving
  when Newton fails.
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

dev.step(dt, { from: sol });                     // transient
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

## 7. Deliverables and packaging

- **Repo layout:** `src/` (core: grid, SG fluxes, assembly, block-tridiagonal solver, Newton,
  transient); `test/` (node's built-in test runner, no deps); `bench/`; `examples/` (standalone
  HTML demos importing from jsdelivr: pn junction with a bias slider; electrolyte cell with a
  load slider; double layer at a blocking electrode).
- **Docs:** README (what it is and isn't, conventions, API, validation table, pointers to
  heavier tools) and JSDoc types (optionally a `.d.ts`).
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

- **M0: equilibrium.** Grid, regions, links, conserved inventories, nonlinear
  Poisson–Boltzmann with flat μ̄_i fixed by links or inventories; materials and interface
  alignment. Tests 1, 2, 3, 4, 8 (built-in potential), 9.
- **M1: steady transport.** SG fluxes, Newton on (φ̂, η_i), voltage/current/load circuit, warm
  starts. Tests 5, 6, 13.
- **M2: reactions.** Bulk mass-action and interfacial Butler–Volmer. Tests 7, 8 (J–V), 10.
- **M3: transient.** Test 12 and the conservation checks.
- **M4: ship 0.1.0.** Examples, docs, benchmark, npm publish.
- **Later:** quasi-neutral mode, degenerate statistics, concentration-dependent D, non-ideal
  activities, non-isothermal transport, impedance (small-signal AC) via the same Jacobian.
  (Advection v(x) and eddy mixing D_mix(x) are cheap and could land in M1 or M2.)

## 9. Questions to settle with the owner early

Settled (2026-09-29): primary branch is `master`; honest-thermodynamics principles, with the
API in μ̄ and V_i as optional views (§3); blocked
species are conserved spectators and all-blocked islands are valid with gate anchoring (§3);
milestones are a brainstorm, free to reshuffle.

- Repo location and npm ownership (who publishes; the name `driftlet` was free on npm as of
  2026-09-29; the fallback is a scoped name).
- Which example devices matter most for the first demos (this sets which of advection/mixing,
  reactions, transients get polished first).
