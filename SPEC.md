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

- `μ̄_i = μ°_i(region) + z_i F ψ + RT ln(c_i / c_ref,i)`, so
  `c_i = c_ref,i · exp( (μ̄_i − μ°_i − z_i F ψ) / RT )`.
- Output conveniences for charged species: the species voltage `V_i = μ̄_i / (z_i F)` and the
  standard level `V°_i = ψ + μ°_i / (z_i F)`. Users who think in voltages scale by z_iF; the
  core never needs to.
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
  conduct: in neutral regions its net current vanishes). Default 0 everywhere.
- Out of scope (non-local): recirculation loops or one well-mixed reservoir feeding several
  faces. A few such couplings could be added later via bordering, like the circuit unknowns.

**Electrostatics.** `−∂/∂x (ε ∂ψ/∂x) = F Σ z_i c_i + ρ_fixed(x)`, where `ε` and the fixed
charge `ρ_fixed` (doping, ionomer charge, …) are set per region.

**Regions.** A device is a sequence of regions, each with its own ε, fixed charge, and
per-species parameters (D, `μ°_i`, or "absent/blocked"). A step in `μ°_i` between regions is a
band offset or solvation step. Species that can't enter a region are blocked at its face
(zero flux).

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

**Contacts (domain ends).**
- **reservoir/bath:** every species' μ̄_i and c_i fixed (a well-stirred solution);
- **ohmic (semiconductor):** local equilibrium + neutrality fix ψ and all μ̄_i, given the applied
  voltage;
- **electrode (metal):** the electron μ̄ (i.e. terminal voltage) is set by the circuit; ions are blocked or take
  part in an interfacial reaction;
- **blocking:** zero flux for the listed species, with a surface charge or potential condition
  for ψ.

**Circuit.** Choose one: applied voltage; applied current (galvanostatic); or an external load
resistor R (V_term = I·R). The latter two add one scalar unknown, handled by bordering the
block system.

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
  "drift" across a face combines `z_i F Δψ + Δμ°_i` (over RT) with the advective Péclet number;
  for neutral species with no flow it reduces to plain diffusion. Offset jumps between regions must be handled
  correctly: SG across a face with a step is exact in the exponential-fitting sense. Evaluate the
  Bernoulli function `B(x) = x/(eˣ−1)` stably (series near 0, asymptotics for large |x|).
- **Unknowns (per node):** the dimensionless potential `ψ̂ = Fψ/RT` and the dimensionless
  electrochemical potentials `η_i = μ̄_i / RT`. Concentrations can span 40 orders of magnitude
  (minority carriers), and these log-like variables keep them well conditioned. At equilibrium
  every `μ̄_i` is exactly flat. Scale everything internally (thermal energy, reference
  concentration, Debye length).
- **Solver:** damped Newton, limiting each update to a few thermal voltages per iteration (or
  Bank–Rose damping). The Jacobian is block-tridiagonal with block size `1 + nSpecies`; solve it
  with block Thomas (a small dense LU per block). Circuit constraints are extra scalar
  unknowns, handled by bordering (Schur complement).
- **Steady-state strategy:** first solve for equilibrium (flat V_i fixed by the contacts, a
  nonlinear Poisson–Boltzmann problem), then use continuation in bias or current.
  `solve({ warm: previous })` must reuse the previous solution; that is what makes slider
  interaction fast (typically 2–5 iterations).
- **Transient:** backward Euler with a BDF2 option, finite-volume mass matrix, and step halving
  when Newton fails. Total charge and mass must be conserved to round-off (test it).
- **Quasi-neutral mode (after v0):** replace Poisson by local neutrality
  `F Σ z_i c_i + ρ_fixed = 0` at each node. This is much cheaper and has no Debye layers; ψ then
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
  contacts: { left: { type: 'electrode', /* … */ }, right: { type: 'electrode' } },
  circuit: { mode: 'load', R: 50 },            // or { mode: 'voltage', V } / { mode: 'current', I }
  grid: { refineNear: 'interfaces', minSpacing: 1e-9, maxSpacing: 5e-6 },
});

let sol = dev.solve();                          // steady state
dev.set({ circuit: { mode: 'load', R: 20 } });
sol = dev.solve({ warm: sol });                 // fast re-solve

// Plotting-ready output:
sol.x; sol.psi;                                  // Float64Arrays
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

1. **Equilibrium invariance:** no bias ⇒ every μ̄_i flat to ≤1e-9·RT and all fluxes zero, for
   every example device.
2. **Nonlinear Poisson–Boltzmann / Gouy–Chapman:** diffuse layer at a charged blocking wall in a
   1:1 electrolyte; the potential profile and differential capacitance match the analytic
   formulas. Its linear limit gives Debye screening.
3. **Donnan potential** across a fixed-charge membrane between two baths vs the analytic value.
4. **Planck liquid junction:** steady diffusion junction of a binary salt between two
   reservoirs. The diffusion potential is `(RT/F)(t₊ − t₋)·ln(c₁/c₂)` for a 1:1 salt.
5. **Concentration polarization / limiting current:** a symmetric metal | binary salt | metal
   cell with deposition/dissolution. Match the analytic quasi-neutral profile, the I–V
   relation and the limiting current `i_lim = z F D_salt c₀ · 2 / ((1 − t₊) L)` (derive and
   check the exact form in the test).
6. **Butler–Volmer interface:** a single electrode facing a reservoir; the current–overpotential
   curve is exact.
7. **pn junction:** built-in potential `(kT/q) ln(N_A N_D / n_i²)` exactly; depletion width
   within the depletion approximation's accuracy; long-diode Shockley J–V within a few % at low
   injection.
8. **Bulk reaction:** the mass-action law holds at equilibrium (e.g. `c_H+ c_OH− = K_w`), and
   homogeneous relaxation matches the analytic time constant.
9. **Advection–diffusion:** Levich limiting current at a rotating disk electrode (imposed
   axial velocity profile), and a Nernst layer from an eddy-diffusivity profile. The eddy
   diffusivity must add no current in a neutral bulk.
10. **Transient:** a blocking electrolyte cell charging with the analytic time constant
   (≈ λ_D·L/D in the small-signal limit), with exact conservation of total charge and mass.
11. **Grid convergence:** observed order ≈ 2 on smooth problems.
12. **Cross-check (optional, documented):** reproduce one published example from
    ChargeTransport.jl or Driftfusion within stated tolerance.
13. **Performance benchmark:** time per Newton iteration and per warm re-solve on reference
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

- **M0: equilibrium.** Grid, regions, nonlinear Poisson–Boltzmann with flat V_i fixed by
  contacts. Tests 1, 2, 3, 7 (built-in potential).
- **M1: steady transport.** SG fluxes, Newton on (ψ̂, η_i), contacts, voltage/current/load
  circuit, warm starts. Tests 4, 5, 7 (J–V), 11.
- **M2: reactions.** Bulk mass-action and interfacial Butler–Volmer. Tests 6, 8.
- **M3: transient.** Test 10 and the conservation checks.
- **M4: ship 0.1.0.** Examples, docs, benchmark, npm publish.
- **Later:** quasi-neutral mode, degenerate statistics, concentration-dependent D, non-ideal
  activities, non-isothermal transport, impedance (small-signal AC) via the same Jacobian.
  (Advection v(x) and eddy mixing D_mix(x) are cheap and could land in M1 or M2.)

## 9. Questions to settle with the owner early

- Repo location and npm ownership (who publishes; the name `driftlet` was free on npm as of
  2026-09-29; the fallback is a scoped name).
- Which example devices matter most for the first demos (this sets which of advection/mixing,
  reactions, transients get polished first).
