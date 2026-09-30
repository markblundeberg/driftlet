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
3. **Bulk standard potentials only mean something in neutral combinations**
   (μ°_e⁻ + μ°_h⁺, μ°_Na⁺ + μ°_Cl⁻, μ°_Li⁺ + μ°_e⁻). How charged levels line up across a
   boundary between different materials is a property of that interface, and you must give it.
   There's no Anderson rule, no Schottky–Mott rule, and no implied common vacuum level. Each
   interface's alignment is independent, even when the same pair of materials meets twice
   (e.g. with different adsorbates). The [alignment guide](alignment.md) gives vacuum-level
   recipes for when that's all you have.
4. **A work function belongs to a surface, not a material.**

## Electrochemical potential and its views

For species *i* with charge number z_i, in a material with standard potential μ°_i and
reference concentration c_ref,i, with ideal (dilute) statistics:

```
μ̄_i = μ°_i + z_i F φ + RT ln(c_i / c_ref,i)
```

The API speaks μ̄ in J/mol (`units.eV` converts per-particle energies). Two voltage views are
provided for charged species, as in [ESBD](https://marklundeberg.com/esbd/) diagrams:

- the **species voltage** `V_i = μ̄_i / (z_i F)`, which is what a voltmeter reads for
  electrons, and the analogous quantity a reversible electrode reads for an ion;
- the **standard level** `V°_i = φ + μ°_i / (z_i F)`, which plays the part of a band edge.

Plotted as voltages, a semiconductor's conduction-band level V°_e⁻ sits *below* its
valence-band level V°_h⁺. That's expected, since V = −E/e for electrons.

Electrons and holes are ordinary species: z = ∓1, c_ref = effective density of states, and
μ° = the band edges as molar energies (conduction band for e⁻, minus the valence band for h⁺).

## Transport

Each species moves down its own electrochemical potential:

```
N_i = −(D_i c_i / RT) ∇μ̄_i        (particle flux, mol/(m²·s))
J_i = z_i F N_i = −σ_i ∇V_i,       σ_i = z_i² F² D_i c_i / RT
```

Drift and diffusion are one law, and the Einstein relation is automatic. D here multiplies
c∇μ̄/RT. For ideal statistics it's the ordinary diffusivity. With non-ideal statistics
(coming), it becomes a mobility coefficient, and the measured chemical diffusivity is D times
the thermodynamic factor.

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

A face between a strictly neutral (ε = 0) material and one with ε > 0, such as an electrode and
an electrolyte, is naturally `dipole`. The neutral side holds the surface charge, and the
alignment plays the role of a potential of zero charge or work function.

## Contacts and terminals

A contact has a **terminal voltage** V, the μ̄ of a named terminal species, and a link for every
species plus one for φ.

- **Fixed species** sit at `V_i = V + offset_i`. The offset belongs to the outside phase: it's
  the chemical potential, per charge, of whatever neutral combination the species forms with
  the terminal species there. It's 0 for the terminal species itself. For an ion at a
  reversible electrode, Mⁿ⁺ + n e⁻ ⇌ M(s) gives `μ_M/(nF)`, which is 0 on the usual table
  convention: E° is carried by the ion's μ° in the solution, not by the offset. It's never
  defaulted for any other species. Neutral species are fixed by an absolute μ̄.
- **Ohmic contacts** to a semiconductor are just fixed V_e⁻ and/or fixed V_h⁺:
  - both fixed at the same metal gives V_h⁺ = V_e⁻, an infinite-recombination contact
    (n·p = n_i² there);
  - one alone gives a selective contact.
- **Baths** are a shorthand. From a neutral composition and a reference species (the ion a
  reversible reference electrode would sense, e.g. Cl⁻ for Ag/AgCl), driftlet works out every
  species' fixed level.
- **Gates and Stern layers** are capacitive φ links to a conductor at the terminal voltage.
  Their `zeroCharge` is the value of V − φ_edge at which the interface is uncharged: the
  flat-band voltage, or the potential of zero charge. Like every alignment, it's a property
  of that interface.
- **Electrode reactions** take their electrons from the metal at μ̄_e = −F·V.

With ion μ° from the usual SHE-based tables, the solution's bookkeeping φ is on the SHE scale:
a standard hydrogen electrode in that solution would sit at V = φ. So electrode voltages and
`zeroCharge` values can be quoted against SHE directly.

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
