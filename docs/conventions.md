# Conventions

Conventions are where users of any drift–diffusion code go wrong, so here they all are in one
place. driftlet's are chosen to be thermodynamically honest, and some of them are deliberately
stricter than you may be used to.

## Principles

1. **Every control and every measurable is an electrochemical potential** $`\bar\mu`$, or a
   difference of them. A bias voltage sets a difference in $`\bar\mu`$ of electrons between
   terminals, $`V = -\Delta\bar\mu_{\mathrm{e}^-}/F`$. A gate sets its metal's $`\bar\mu_{\mathrm{e}^-}`$.
   Nothing you set is an electrostatic potential.
2. **$`\phi`$ is bookkeeping.** Any consistent $`\phi`$ gives identical physics. Each material's standard
   chemical potentials fix where its $`\phi`$ sits relative to its own levels, so $`\phi`$ jumps at
   every interface between different materials, by that interface's dipole. How big a dipole is
   depends on how each material was anchored, which is why it's arbitrary, correctly so. If you
   do have true mean inner potentials, you may anchor to them, but nothing requires it.
   Extrathermodynamic conventions are one way people settle a "real" $`\phi`$ per material. TATB,
   for instance, fixes single-ion transfer energies between solvents by assuming
   $`\Delta G_{\mathrm{tr}}(\mathrm{Ph_4As}^+) = \Delta G_{\mathrm{tr}}(\mathrm{Ph_4B}^-)`$. That's a
   statement about bulk single-ion energies, i.e. about how each solvent's $`\mu^\circ`$ are
   anchored. It says nothing about the dipole at an actual interface, which stays a separate
   input.
3. **Bulk standard potentials only mean something in neutral combinations**
   ($`\mu^\circ_{\mathrm{e}^-} + \mu^\circ_{\mathrm{h}^+}`$,
   $`\mu^\circ_{\mathrm{Na}^+} + \mu^\circ_{\mathrm{Cl}^-}`$,
   $`\mu^\circ_{\mathrm{Li}^+} + \mu^\circ_{\mathrm{e}^-}`$). How charged levels line up across a
   boundary between different materials is a property of that interface, and you must give it.
   There's no Anderson rule, no Schottky–Mott rule, and no implied common vacuum level. Each
   interface's alignment is independent, even when the same pair of materials meets twice
   (e.g. with different adsorbates). When vacuum-level estimates are all you have, helpers in
   `driftlet/kit` turn anchors and offsets into an alignment (see the
   [alignment guide](alignment.md)).
4. **A work function belongs to a surface, not a material.**

## Electrochemical potential and its views

For species *i* with charge number $`z_i`$, in a material where its standard chemical potential
is $`\mu^\circ_i`$, the electrochemical potential splits as

```math
\bar\mu_i = \mu^\circ_i + z_i F \phi + RT \zeta_i
```

The reduced chemical potential $`\zeta_i`$ is what a material's **statistics** turn into a
concentration, $`c_i = c_i(\zeta)`$. The default is ideal (dilute) statistics,

```math
c_i = c_{\mathrm{ref},i} \cdot e^{\zeta_i}, \quad\text{i.e.}\quad \bar\mu_i = \mu^\circ_i + z_i F \phi + RT \ln(c_i / c_{\mathrm{ref},i})
```

with $`c_{\mathrm{ref},i}`$ the reference concentration that $`\mu^\circ_i`$ refers to. Non-ideal
statistics (degenerate carriers, crowding, intercalation hosts, activity coefficients) are chosen
per material; see [statistics](statistics.md). They all reduce to the ideal form in the dilute
limit, so $`\mu^\circ`$ and $`c_{\mathrm{ref}}`$ mean the same thing in every model: the dilute
(Henry's-law) reference.

The API speaks $`\bar\mu`$ in J/mol (`units.eV` converts per-particle energies). Two voltage views
are provided for charged species (what they show, and why they're drawn, is in
[reading level diagrams](visualization.md)):

- the **species voltage** $`V_i = \bar\mu_i / (z_i F)`$, which is what a voltmeter reads for
  electrons, and the analogous quantity a reversible electrode reads for an ion;
- the **standard level** $`V^\circ_i = \phi + \mu^\circ_i / (z_i F)`$, which plays the part of a
  band edge.

Plotted as voltages, a semiconductor's conduction-band level $`V^\circ_{\mathrm{e}^-}`$ sits
*below* its valence-band level $`V^\circ_{\mathrm{h}^+}`$. That's expected, since $`V = -E/e`$ for
electrons.

Electrons and holes are ordinary species: $`z = \mp 1`$, $`c_{\mathrm{ref}}`$ = effective density of
states, and $`\mu^\circ`$ = the band edges as molar energies (conduction band for e⁻, minus the
valence band for h⁺). With Fermi–Dirac statistics the same parameters give
$`c = N_c \mathcal{F}_{1/2}(\zeta)`$.

## Transport

Each species moves down its own electrochemical potential:

```math
\begin{aligned}
N_i &= -(D_i c_i / RT) \nabla\bar\mu_i && \quad\text{(particle flux, mol/(m}^2\cdot\text{s))} \\
J_i &= z_i F N_i = -\sigma_i \nabla V_i, \qquad \sigma_i = z_i^2 F^2 D_i c_i / RT
\end{aligned}
```

Drift and diffusion are one law, and the Einstein relation is automatic. $`D`$ is the coefficient
of $`c\nabla\bar\mu/RT`$, i.e. $`RT`$ times the mobility per mole, and it's the same parameter
whatever the statistics. For ideal statistics it's the ordinary diffusivity. With non-ideal
statistics the measured chemical diffusivity is $`D`$ times the thermodynamic factor,
$`\partial \ln a/\partial \ln c`$: the generalised Einstein relation, as for degenerate electrons.

Two optional transport terms describe a moving fluid:

- **Advection** by an imposed velocity $`v`$ adds $`v \cdot c_i`$ to every mobile species' flux.
- **Eddy mixing** adds $`-(D_{\mathrm{mix}}/RT) P \nabla\bar\mu`$, with
  $`P = C - C z z^\mathsf{T} C / (z^\mathsf{T} C z)`$ and $`C = \mathrm{diag}(c)`$.
  - It's an Onsager term: symmetric, positive semi-definite, and zero at equilibrium (flat
    $`\bar\mu`$), so double layers at rest are untouched.
  - The projection removes any current ($`z^\mathsf{T} P = 0`$). Eddies move neutral fluid, so
    they mix salt and stir concentration gradients but conduct nothing.
  - For a neutral dilute species it's exactly $`-D_{\mathrm{mix}} \nabla c`$. For a binary salt it
    adds $`D_{\mathrm{mix}}`$ to the salt's ambipolar diffusivity and leaves diffusion potentials
    unchanged.

## Materials, regions, interfaces

A **material** holds bulk properties ($`\varepsilon`$, and $`D`$, $`\mu^\circ`$, $`c_{\mathrm{ref}}`$ per
species). A **region** is a material plus a length, fixed charge and initial composition. Every
**interface** between different materials takes exactly one alignment (see principle 3), unless
its electrostatic law is `neutral`, in which case the alignment has no effect and is not allowed.
The three laws:

- `pinned`: $`\phi`$ jumps by exactly the alignment, and the double layers on either side are
  resolved by the grid.
- `neutral`: the macroscopic idealisation. No charge sits at the face, and $`\phi`$ jumps by
  whatever local neutrality on each side requires (a Donnan or Galvani step).
- `capacitive`: a Helmholtz layer of capacitance $`C`$, charged by any departure of the jump
  from the alignment.

A **conductor** (a metal, or a fast ion conductor) is a material of its own kind: its one level
is its carrier's (a metal's Fermi level), it has no $`\phi`$, and a face beside it ties the other
side's $`\phi`$ to that level through a capacitance (see the
[device reference](device.md#faces-next-to-a-conductor)).

A face between a strictly neutral ($`\varepsilon = 0`$) material and one with $`\varepsilon > 0`$,
such as an electrode and an electrolyte, is naturally `pinned`. The neutral side holds the
surface charge, and the alignment plays the role of a potential of zero charge or work function.

## Contacts and terminals

A contact is an interface with an **outside phase** whose levels are known, like one more region
with a fully known node. The outside phase's charged levels form a rigid ladder,
$`V_i = V + \mathtt{offset}_i`$, which the external circuit slides by the **terminal voltage** $`V`$
(the voltage of a named terminal species). The laws joining it to the device (equilibrium,
blocked, conductance; and for $`\phi`$, bulk, neutral, capacitive, dipole) are the same as at
internal faces. Reactions sit at faces between regions: an electrode is a conductor region with
reactions at its face.

- **Species in equilibrium** with the outside phase sit at $`V_i = V + \mathtt{offset}_i`$. The
  offset belongs to the outside phase: it's the chemical potential, per charge, of whatever
  neutral combination the species forms with the terminal species there. It's 0 for the
  terminal species itself. For an ion at a reversible electrode,
  $`\mathrm{M}^{n+} + n\,\mathrm{e}^- \rightleftharpoons \mathrm{M(s)}`$ gives
  $`\mu_{\mathrm{M}}/(nF)`$, which is 0 on the usual table convention: $`E^\circ`$ is carried by the
  ion's $`\mu^\circ`$ in the solution, not by the offset. It's never defaulted for any other
  species. Neutral species are fixed by an absolute $`\bar\mu`$.
- **Ohmic contacts** to a semiconductor are just $`V_{\mathrm{e}^-}`$ and/or $`V_{\mathrm{h}^+}`$ in
  equilibrium with the metal:
  - both, at the same metal, gives $`V_{\mathrm{h}^+} = V_{\mathrm{e}^-}`$, an
    infinite-recombination contact ($`n \cdot p = n_i^2`$ there);
  - one alone gives a selective contact.
- **Baths** are a shorthand. From a neutral composition and a reference species (the ion a
  reversible reference electrode would sense, e.g. Cl⁻ for Ag/AgCl), driftlet works out every
  species' equilibrium level.
- **Gates and Stern layers** are capacitive $`\phi`$ laws to a conductor at the terminal voltage.
  Their `zeroCharge` is the value of $`V - \phi_{\mathrm{edge}}`$ at which the interface is
  uncharged. For an electrode in a solution with table $`\mu^\circ`$ (where $`V - \phi`$ reads
  against SHE) that's the potential of zero charge; for a gate it sets the flat-band voltage,
  $`V_{\mathrm{FB}} = \mathtt{zeroCharge} + \phi_{\mathrm{bulk}}`$, with $`\phi_{\mathrm{bulk}}`$ in
  the semiconductor's own anchoring (−1.00 V in p-type silicon at 1e17 cm⁻³ with
  $`\mu^\circ_{\mathrm{e}^-} = 0`$). Like every alignment, it's a property of that interface. The
  `pinned` law is the $`C \to \infty`$ limit (e.g. a Schottky barrier), which is what a
  "fixed $`\phi`$" boundary honestly means.
- **`bulk`** says the end node is plain bulk (locally neutral, no double layer), while
  **`neutral`** says no charge sits at the face ($`D = 0`$), as for internal faces.
- **Electrode reactions** take their electrons from the conductor region at its Fermi level,
  which the contact behind it holds at $`\bar\mu_{\mathrm{e}^-} = -F \cdot V`$.

The electronic level of a standard hydrogen electrode in a solution is
$`V^\circ_{\mathrm{e}^-}(\mathrm{SHE}) = \phi + (\mu^\circ_{\mathrm{H}^+} - \tfrac{1}{2}\mu^\circ_{\mathrm{H_2}})/F`$.
The usual tables set $`\mu^\circ_{\mathrm{H}^+} = 0`$ (so that $`E^\circ(\mathrm{SHE}) = 0`$) and
$`\mu^\circ_{\mathrm{H_2}} = 0`$ for the element in its standard state, and then
$`V^\circ_{\mathrm{e}^-}(\mathrm{SHE}) = \phi`$. So with ion $`\mu^\circ`$ from such tables, $`V - \phi`$
in a solution is a potential against SHE. A metal's $`V - \phi_{\mathrm{edge}}`$ at zero charge (a
`zeroCharge`, i.e. a pzc) reads directly on the SHE scale.

## Terminals and signs

- A terminal's voltage $`V`$ is the shift of its outside phase's ladder, measured as the voltage of
  its terminal species. Every voltage is absolute (no terminal is special), and at least one
  terminal must be held at one.
- **A terminal's current is into the device** (conduction plus displacement), so in steady
  state they sum to zero. $`V \cdot I`$ is the electrical work done through it; an outside phase
  that exchanges several species also trades chemical free energy through its offsets.
- **`current` is positive toward +x** through the device (into it at the left, out at the
  right). A forward-biased pn junction with its p side on the right carries negative current.
- `terminalVoltage` is $`V_{\mathrm{right}} - V_{\mathrm{left}}`$.

## Units

SI throughout: m, mol/m³, J/mol, V, A/m², F/m² (per-area capacitance), C/m³ (fixed charge),
C/m² (sheet charge), S/m² (interface conductance), mol/(m²·s) (surface rate constants),
s. Temperature is uniform (default 298.15 K).
