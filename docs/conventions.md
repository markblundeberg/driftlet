# Conventions

Conventions are where users of any drift–diffusion code go wrong, so here they all are in one
place. driftlet's are chosen to be thermodynamically honest, and some of them are deliberately
stricter than you may be used to.

## Principles

1. **Every control and every measurable is an electrochemical potential** μ̄, or a
   difference of them. A bias voltage sets a difference in μ̄ of electrons between terminals,
   `V = −Δμ̄_e⁻/F`. A gate sets its metal's μ̄_e⁻. Nothing you set is an electrostatic potential.
2. **φ is bookkeeping.** Any consistent φ gives identical physics. Each material's standard
   chemical potentials fix where its φ sits relative to its own levels, so φ jumps at every
   interface between different materials, by that interface's dipole. How big a dipole is
   depends on how each material was anchored, which is why it's arbitrary, correctly so. If you
   do have true mean inner potentials, you may anchor to them, but nothing requires it.
   Extrathermodynamic conventions are one way people settle a "real" φ per material. TATB, for
   instance, fixes single-ion transfer energies between solvents by assuming
   ΔG_tr(Ph₄As⁺) = ΔG_tr(Ph₄B⁻). That's a statement about bulk single-ion energies, i.e. about
   how each solvent's μ° are anchored. It says nothing about the dipole at an actual interface,
   which stays a separate input.
3. **Bulk standard potentials only mean something in neutral combinations**
   (μ°_e⁻ + μ°_h⁺, μ°_Na⁺ + μ°_Cl⁻, μ°_Li⁺ + μ°_e⁻). How charged levels line up across a
   boundary between different materials is a property of that interface, and you must give it.
   There's no Anderson rule, no Schottky–Mott rule, and no implied common vacuum level. Each
   interface's alignment is independent, even when the same pair of materials meets twice
   (e.g. with different adsorbates). When vacuum-level estimates are all you have, helpers in
   `driftlet/kit` turn anchors and offsets into an alignment (see the
   [alignment guide](alignment.md)).
4. **A work function belongs to a surface, not a material.**

## Electrochemical potential and its views

For species *i* with charge number z_i, in a material where its standard chemical potential
is μ°_i, the electrochemical potential splits as

```
μ̄_i = μ°_i + z_i F φ + RT ζ_i
```

The reduced chemical potential ζ_i is what a material's **statistics** turn into a
concentration, c_i = c_i(ζ). The default is ideal (dilute) statistics,

```
c_i = c_ref,i · e^{ζ_i},   i.e.   μ̄_i = μ°_i + z_i F φ + RT ln(c_i / c_ref,i)
```

with c_ref,i the reference concentration that μ°_i refers to. Non-ideal statistics
(degenerate carriers, crowding, intercalation hosts, activity coefficients) are chosen per
material; see [statistics](statistics.md). They all reduce to the ideal form in the dilute
limit, so μ° and c_ref mean the same thing in every model: the dilute (Henry's-law) reference.

The API speaks μ̄ in J/mol (`units.eV` converts per-particle energies). Two voltage views are
provided for charged species, as in [ESBD](https://marklundeberg.com/esbd/) diagrams:

- the **species voltage** `V_i = μ̄_i / (z_i F)`, which is what a voltmeter reads for
  electrons, and the analogous quantity a reversible electrode reads for an ion;
- the **standard level** `V°_i = φ + μ°_i / (z_i F)`, which plays the part of a band edge.

Plotted as voltages, a semiconductor's conduction-band level V°_e⁻ sits *below* its
valence-band level V°_h⁺. That's expected, since V = −E/e for electrons.

Electrons and holes are ordinary species: z = ∓1, c_ref = effective density of states, and
μ° = the band edges as molar energies (conduction band for e⁻, minus the valence band for h⁺).
With Fermi–Dirac statistics the same parameters give c = N_c 𝓕_{1/2}(ζ).

## Transport

Each species moves down its own electrochemical potential:

```
N_i = −(D_i c_i / RT) ∇μ̄_i        (particle flux, mol/(m²·s))
J_i = z_i F N_i = −σ_i ∇V_i,       σ_i = z_i² F² D_i c_i / RT
```

Drift and diffusion are one law, and the Einstein relation is automatic. D is the coefficient
of c∇μ̄/RT, i.e. RT times the mobility per mole, and it's the same parameter whatever the
statistics. For ideal statistics it's the ordinary diffusivity. With non-ideal statistics the
measured chemical diffusivity is D times the thermodynamic factor, ∂ln a/∂ln c: the
generalised Einstein relation, as for degenerate electrons.

Two optional transport terms describe a moving fluid:

- **Advection** by an imposed velocity v adds v·c_i to every mobile species' flux.
- **Eddy mixing** adds `−(D_mix/RT) P ∇μ̄`, with `P = C − C z zᵀ C / (zᵀ C z)` and
  C = diag(c).
  - It's an Onsager term: symmetric, positive semi-definite, and zero at equilibrium (flat
    μ̄), so double layers at rest are untouched.
  - The projection removes any current (zᵀP = 0). Eddies move neutral fluid, so they mix salt
    and stir concentration gradients but conduct nothing.
  - For a neutral dilute species it's exactly −D_mix ∇c. For a binary salt it adds D_mix to
    the salt's ambipolar diffusivity and leaves diffusion potentials unchanged.

## Materials, regions, interfaces

A **material** holds bulk properties (ε, and D, μ°, c_ref per species). A **region** is a
material plus a length, fixed charge and initial composition. Every **interface** between
different materials takes exactly one alignment (see principle 3), unless its electrostatic law
is `neutral`, in which case the alignment has no effect and is not allowed. The three laws:

- `dipole`: φ jumps by exactly the alignment, and the double layers on either side are
  resolved by the grid.
- `neutral`: the macroscopic idealisation. No charge sits at the face, and φ jumps by whatever
  local neutrality on each side requires (a Donnan or Galvani step).
- `capacitive`: a Helmholtz layer of capacitance C, charged by any departure of the jump
  from the alignment.

A **conductor** (a metal, or a fast ion conductor) is a material of its own kind: its one level
is its carrier's (a metal's Fermi level), it has no φ, and a face beside it ties the other side's
φ to that level through a capacitance (see the
[device reference](device.md#faces-next-to-a-conductor)).

A face between a strictly neutral (ε = 0) material and one with ε > 0, such as an electrode and
an electrolyte, is naturally `dipole`. The neutral side holds the surface charge, and the
alignment plays the role of a potential of zero charge or work function.

## Contacts and terminals

A contact is an interface with an **outside phase** whose levels are known, like one more region
with a fully known node. The outside phase's charged levels form a rigid ladder,
`V_i = V + offset_i`, which the circuit slides by the **terminal voltage** V (the voltage of a
named terminal species). The laws joining it to the device (equilibrium, blocked, conductance;
and for φ, bulk, neutral, capacitive, dipole) are the same as at internal faces. Reactions sit
at faces between regions: an electrode is a conductor region with reactions at its face.

- **Species in equilibrium** with the outside phase sit at `V_i = V + offset_i`. The offset belongs to the outside phase: it's
  the chemical potential, per charge, of whatever neutral combination the species forms with
  the terminal species there. It's 0 for the terminal species itself. For an ion at a
  reversible electrode, Mⁿ⁺ + n e⁻ ⇌ M(s) gives `μ_M/(nF)`, which is 0 on the usual table
  convention: E° is carried by the ion's μ° in the solution, not by the offset. It's never
  defaulted for any other species. Neutral species are fixed by an absolute μ̄.
- **Ohmic contacts** to a semiconductor are just V_e⁻ and/or V_h⁺ in equilibrium with the metal:
  - both, at the same metal, gives V_h⁺ = V_e⁻, an infinite-recombination contact
    (n·p = n_i² there);
  - one alone gives a selective contact.
- **Baths** are a shorthand. From a neutral composition and a reference species (the ion a
  reversible reference electrode would sense, e.g. Cl⁻ for Ag/AgCl), driftlet works out every
  species' equilibrium level.
- **Gates and Stern layers** are capacitive φ laws to a conductor at the terminal voltage.
  Their `zeroCharge` is the value of V − φ_edge at which the interface is uncharged: the
  flat-band voltage, or the potential of zero charge. Like every alignment, it's a property
  of that interface. The pinned `dipole` law is the C → ∞ limit (e.g. a Schottky barrier),
  which is what a "fixed φ" boundary honestly means.
- **`bulk`** says the end node is plain bulk (locally neutral, no double layer), while
  **`neutral`** says no charge sits at the face (D = 0), as for internal faces.
- **Electrode reactions** take their electrons from the conductor region at its Fermi level,
  which the contact behind it holds at μ̄_e = −F·V.

The electronic level of a standard hydrogen electrode in a solution is
V°_e⁻(SHE) = φ + (μ°_H⁺ − ½μ°_H₂)/F. The usual tables set μ°_H⁺ = 0 (so that E°(SHE) = 0) and
μ°_H₂ = 0 for the element in its standard state, and then V°_e⁻(SHE) = φ. So with ion μ° from
such tables, V − φ in a solution is a potential against SHE. A metal's V − φ_edge at zero
charge (a `zeroCharge`, i.e. a pzc) reads directly on the SHE scale.

## Circuits and signs

- The left terminal is the reference.
- **Current is positive toward +x.** A forward-biased pn junction with its p side on the right
  carries negative current.
- The terminal voltage is `V_right − V_left`.
- Galvanostatic and load modes float the right terminal.

## Units

SI throughout: m, mol/m³, J/mol, V, A/m², F/m² (per-area capacitance), C/m³ (fixed charge),
C/m² (sheet charge), S/m² (interface conductance), mol/(m²·s) (surface rate constants),
s. Temperature is uniform (default 298.15 K).
