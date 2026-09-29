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
- **1D only.** The Jacobian must stay **block-tridiagonal** (plus at most a few bordering scalar
  unknowns for circuit constraints). That is what keeps it fast and simple. Features that would
  break this structure are out of scope.
- **Speed:** with ~300 nodes and up to 6 species, one Newton iteration should take well under
  1 ms in V8, and a warm-started re-solve after a small parameter change should take under
  ~5 ms. Inner loops should be allocation-free (`Float64Array`, reused buffers).
- **It must work inside a Web Worker** (no DOM access anywhere in the core).
- **Size budget:** aim for ≲ 30 kB minified.

## 3. Physics and conventions

Document all of this in the README; conventions are where users go wrong.

**Species.** Each species *i* has a charge number `z` (any integer, **including 0**: neutral
mobile species such as dissolved gases or water are supported from v0), a diffusivity `D` (per
region), and, per region, a **standard chemical potential** `μ°_i`. The primary variable is the
**electrochemical potential** `μ̄_i`. Ideal (Boltzmann/dilute) statistics:

- `μ̄_i = μ°_i(region) + z_i F φ + RT ln(c_i / c_ref,i)`, so
  `c_i = c_ref,i · exp( (μ̄_i − μ°_i − z_i F φ) / RT )`.
- **Species voltage is the user-facing language.** For charged species, `V_i = μ̄_i / (z_i F)`
  *is* a voltage (for electrons it's exactly what a voltmeter reads; for ions it's the
  analogous quantity a reversible electrode reads), with the standard level
  `V°_i = φ + μ°_i / (z_i F)` playing the part of a band edge. This is the owner's
  [ESBD](https://marklundeberg.com/esbd/) framing, and the API should speak it: contacts,
  circuits and outputs are stated in V_i for charged species (μ̄_i only for neutral ones).
  The core works in `η_i = μ̄_i / RT` (well defined for z = 0) and converts at the edges.
- Electrons and holes are ordinary species (`z = ∓1`), with `c_ref` = the effective density of
  states and `μ°` = the band edges (as molar energies). Plotted as voltages (V_i upward), the
  conduction band sits below the valence band; that's expected.
- Degenerate (Fermi–Dirac) statistics are a later, optional add-on (e.g. the Blakemore
  approximation).

**Transport.** Particle flux per species:
`N_i = −(D_i c_i / RT) ∇μ̄_i + c_i v(x) − D_mix(x) ∇c_i`, with continuity
`∂c_i/∂t + ∂N_i/∂x = R_i`.
- The first term is drift + diffusion in one (Einstein relation automatic). For charged species
  the charge current is `J_i = z_i F N_i = −σ_i ∇V_i`, with `σ_i = z_i² F² D_i c_i / RT`.
- `v(x)`: an optional **imposed advection velocity** (e.g. the axial flow at a rotating disk
  electrode, or flow through a porous layer). It's local, so it keeps the Jacobian
  block-tridiagonal, and SG handles drift + advection as one effective velocity.
- `D_mix(x)`: an optional **eddy (mixing) diffusivity**, large in a stirred bulk and falling to
  zero at walls. It models stirring locally and produces Nernst diffusion layers naturally. It
  acts identically on every species and has **no associated mobility** (it mixes, it doesn't
  conduct: in neutral regions its net current vanishes). Default 0 everywhere. Caveat: it acts
  on ∇c, not ∇μ̄, so it would drive a spurious flux through an *equilibrium* double layer
  (flat μ̄, varying c). Profiles must fall to zero before the Debye layers at walls; the
  grid/region setup should warn when D_mix overlaps a strongly non-neutral zone.
- Out of scope (non-local): recirculation loops or one well-mixed reservoir feeding several
  faces. A few such couplings could be added later via bordering, like the circuit unknowns.

**Electrostatics.** `−∂/∂x (ε ∂φ/∂x) = F Σ z_i c_i + ρ_fixed(x)`, where `ε` and the fixed
charge `ρ_fixed` (doping, ionomer charge, …) are set per region.

**Regions.** A device is a sequence of regions, each with its own ε, fixed charge, and
per-species parameters (D, `μ°_i`, or "absent/blocked"). A step in `μ°_i` between regions is a
band offset or solvation step. Species that can't enter a region are blocked at its face
(zero flux).

**Absent species.** Where a species is absent, its η_i is undefined and has no equation. Two
implementation options: keep a fixed block size with decoupled identity rows, or let the block
size vary per region (block Thomas copes with non-square off-diagonal blocks). Start with
whichever is simpler and benchmark the other; the choice must not leak into the API. Outputs
for an absent species are `NaN` at those nodes, so plots show a gap rather than a fake value.

**Conserved inventories (blocked species).** A species blocked on every side of some connected
stretch of the device (no contact link, no reaction) is a *conserved spectator*: its amount
there never changes, and it just drifts and screens. Transients conserve it automatically. A
steady-state or equilibrium solve cannot determine it (its μ̄ level floats), so each such
inventory is one extra scalar unknown (its level) with one constraint (`∫c dx` = the given
amount), handled by bordering like the circuit unknowns. The amount defaults to what the
initial/warm state holds, or is given explicitly. With bulk reactions, the conserved quantities
are moieties (left null vectors of the stoichiometry restricted to that stretch) rather than
species. A device can be *entirely* blocked (a floating island) provided φ is anchored
electrostatically somewhere (a gate link, below); then the island simply holds a fixed charge.
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
- **species links**, one per species: *fixed* (V_i, or μ̄_i if neutral, pinned at the boundary:
  the normal case), *conductance* G_i to an outside V_i (an ohmic interface resistance),
  *kinetic* (Butler–Volmer, below, i.e. a nonlinear conductance), or *blocked* (G_i = 0);
- **electrostatic link**: *neutral* (φ at the boundary node follows from local neutrality,
  given the fixed V_i: an ohmic contact or well-stirred bath), *gate* (a capacitance per area
  C_g to a gate at V_g, with a gate offset for its work function: `ε ∂φ/∂n = C_g(φ_g − φ)`;
  this also models a Stern/Helmholtz layer), *fixed* φ (the C_g → ∞ limit), or *free*
  (C_g = 0, zero field).
Familiar contacts are presets: bath = all species fixed + neutral; semiconductor ohmic = e⁻/h⁺
fixed + neutral; metal electrode = ions blocked or kinetic, gate link through the Stern layer to
the metal's V_e⁻; ideal insulator = all blocked + gate. Infinite G or C must be imposed exactly
(eliminate the unknown), not with the large-penalty trick the prototype used.

**Circuit.** The terminals are the V of a named species at each end (default e⁻ at a metal; for
an ion-only cell, the ion a reversible reference electrode would sense, e.g. Cl⁻ for Ag/AgCl).
This settles what "terminal voltage" means with no electrons in the device. Choose one:
applied voltage; applied current (galvanostatic); or an external load resistor R
(V_term = I·R). The latter two add one scalar unknown, handled by bordering the block system.
Open circuit (I = 0) is the current mode with I = 0 (liquid junctions, Donnan).

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

- **Grid:** 1D non-uniform node grid; region boundaries fall on cell faces. Include a grid
  builder with geometric/tanh grading toward interfaces and contacts, since Debye lengths
  (nm) are much smaller than devices (µm–mm). Given a target resolution near interfaces, it
  should pick spacing automatically.
- **Discretization:** a finite-volume (box) method with **Scharfetter–Gummel** fluxes (the same
  exponential fitting as Il'in / Allen–Southwell for convection–diffusion). For species *i* the
  "drift" across a face combines `z_i F Δφ + Δμ°_i` (over RT) with the advective Péclet number;
  for neutral species with no flow it reduces to plain diffusion. Offset jumps between regions must be handled
  correctly: SG across a face with a step is exact in the exponential-fitting sense. Evaluate the
  Bernoulli function `B(x) = x/(eˣ−1)` stably (series near 0, asymptotics for large |x|).
- **Unknowns (per node):** the dimensionless potential `φ̂ = Fφ/RT` and the dimensionless
  electrochemical potentials `η_i = μ̄_i / RT`. Concentrations can span 40 orders of magnitude
  (minority carriers), and these log-like variables keep them well conditioned. At equilibrium
  every `μ̄_i` is exactly flat. Scale everything internally (thermal energy, reference
  concentration, Debye length).
- **Solver:** damped Newton, limiting each update to a few thermal voltages per iteration (or
  Bank–Rose damping). The Jacobian is block-tridiagonal with block size `1 + nSpecies`; solve it
  with block Thomas (a small dense LU per block). Circuit constraints and conserved
  inventories are extra scalar unknowns, handled by bordering (Schur complement: one extra
  block-Thomas back-substitution per scalar, reusing the factorisation).
- **Steady-state strategy:** first solve for equilibrium (flat V_i fixed by the contact links
  or, for spectators, by their inventories: a nonlinear Poisson–Boltzmann problem), then use
  continuation in bias or current.
  `solve({ warm: previous })` must reuse the previous solution; that is what makes slider
  interaction fast (typically 2–5 iterations).
- **Transient:** backward Euler with a BDF2 option, finite-volume mass matrix, and step halving
  when Newton fails. Total charge and mass must be conserved to round-off (test it).
- **Quasi-neutral mode (after v0):** replace Poisson by local neutrality
  `F Σ z_i c_i + ρ_fixed = 0` at each node. This is much cheaper and has no Debye layers; φ then
  jumps at interfaces (Donnan, junction potentials), which SG fluxes handle.
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
  regions: [
    { name: 'electrolyte', length: 1e-3, eps: 78.5 * 8.854e-12, fixedCharge: 0,
      species: { 'Zn2+': { D: 0.70e-9, mu0: -147.1e3 }, 'SO4 2-': { D: 1.07e-9, mu0: -744.5e3 } } },
  ],
  interfaceReactions: [ /* at a region face or contact: stoichiometry, i0, alpha */ ],
  bulkReactions: [ /* stoichiometry, kf, neutral activities */ ],
  contacts: {                                  // links (§3); presets expand to these
    left:  { species: { 'Zn2+': { fixed: true }, 'SO4 2-': 'blocked' },
             phi: { gate: { C: 0.2 } },        // Stern layer to the metal, F/m²
             terminal: 'Zn2+' },               // terminal = V of this species here
    right: { species: { 'Zn2+': { fixed: true }, 'SO4 2-': 'blocked' },
             phi: { gate: { C: 0.2 } }, terminal: 'Zn2+' },
  },
  circuit: { mode: 'load', R: 50 },            // or { mode: 'voltage', V } / { mode: 'current', I }
  grid: { refineNear: 'interfaces', minSpacing: 1e-9, maxSpacing: 5e-6 },
});

let sol = dev.solve();                          // steady state
dev.set({ circuit: { mode: 'load', R: 20 } });
sol = dev.solve({ warm: sol });                 // fast re-solve

// Plotting-ready output:
sol.x; sol.phi;                                  // Float64Arrays
sol.mu['Zn2+']; sol.c['Zn2+']; sol.N['Zn2+'];     // μ̄ (J/mol), c, particle flux
sol.V['Zn2+']; sol.Vstd['Zn2+']; sol.J['Zn2+'];    // voltage-scaled views (charged species)
sol.current; sol.terminalVoltage;
sol.interfaces;   // per face: steps in each μ̄_i and V_i (overpotentials, junction/Donnan potentials)
sol.converged; sol.iterations; sol.residual;

dev.step(dt, { from: sol });                     // transient
```

Design goals for the API: plain objects in, typed arrays out, no classes users must subclass.
Everything should be serializable, so a device definition can be posted to a Worker.

## 6. Validation suite (the real deliverable)

Each item is an automated test with a stated tolerance. The README carries a table of them.

1. **Equilibrium invariance:** no bias and no imposed flow ⇒ every μ̄_i flat to ≤1e-9·RT and
   all fluxes zero, for every example device.
2. **Nonlinear Poisson–Boltzmann / Gouy–Chapman:** diffuse layer at a charged blocking wall in a
   1:1 electrolyte; the potential profile and differential capacitance match the analytic
   formulas. Its linear limit gives Debye screening.
3. **Donnan potential** across a fixed-charge membrane between two baths vs the analytic value.
4. **Floating island:** every species blocked, φ anchored only by gate links at both ends.
   Each inventory is conserved to round-off through a gate sweep; Gauss's law holds exactly
   (island charge = −Σ C_g(φ_g − φ_edge)); in the linear limit the island's response matches
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
9. **Bulk reaction:** the mass-action law holds at equilibrium (e.g. `c_H+ c_OH− = K_w`), and
   homogeneous relaxation matches the analytic time constant.
10. **Advection–diffusion:** Levich limiting current at a rotating disk electrode (imposed
   axial velocity profile), and a Nernst layer from an eddy-diffusivity profile. The eddy
   diffusivity must add no current in a neutral bulk.
11. **Transient:** a blocking electrolyte cell charging with the analytic time constant
   (≈ λ_D·L/D in the small-signal limit), with exact conservation of total charge and mass.
12. **Grid convergence:** observed order ≈ 2 on smooth problems.
13. **Cross-check (optional, documented):** reproduce one published example from
    ChargeTransport.jl or Driftfusion within stated tolerance.
14. **Performance benchmark:** time per Newton iteration and per warm re-solve on reference
    devices, tracked across versions.

## 7. Deliverables and packaging

- **Repo layout:** `src/` (core: grid, SG fluxes, assembly, block-tridiagonal solver, Newton,
  transient); `test/` (node's built-in test runner, no deps); `bench/`; `examples/` (standalone
  HTML demos importing from jsdelivr: pn junction with a bias slider; electrolyte cell with a
  load slider; double layer at a blocking electrode).
- **Docs:** README (what it is and isn't, conventions, API, validation table, pointers to
  heavier tools) and JSDoc types (optionally a `.d.ts`).
- **npm:** ESM with an `exports` map, `sideEffects: false`, semver starting at 0.1.0, and a
  CHANGELOG.
- **CI:** GitHub Actions running tests and the benchmark on push.
- **Licence:** 0BSD (public-domain-equivalent; SPDX `0BSD`).

## 8. Milestones

- **M0: equilibrium.** Grid, regions, links, conserved inventories, nonlinear
  Poisson–Boltzmann with flat V_i fixed by links or inventories. Tests 1, 2, 3, 4, 8 (built-in
  potential).
- **M1: steady transport.** SG fluxes, Newton on (φ̂, η_i), voltage/current/load circuit, warm
  starts. Tests 5, 6, 12.
- **M2: reactions.** Bulk mass-action and interfacial Butler–Volmer. Tests 7, 8 (J–V), 9.
- **M3: transient.** Test 11 and the conservation checks.
- **M4: ship 0.1.0.** Examples, docs, benchmark, npm publish.
- **Later:** quasi-neutral mode, degenerate statistics, concentration-dependent D, non-ideal
  activities, non-isothermal transport, impedance (small-signal AC) via the same Jacobian.
  (Advection v(x) and eddy mixing D_mix(x) are cheap and could land in M1 or M2.)

## 9. Questions to settle with the owner early

Settled (2026-09-29): primary branch is `master`; V_i is the API's language (§3); blocked
species are conserved spectators and all-blocked islands are valid with gate anchoring (§3);
milestones are a brainstorm, free to reshuffle.

- Repo location and npm ownership (who publishes; the name `driftlet` was free on npm as of
  2026-09-29; the fallback is a scoped name).
- Which example devices matter most for the first demos (this sets which of advection/mixing,
  reactions, transients get polished first).
