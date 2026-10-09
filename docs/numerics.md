# Numerics

How driftlet turns the physics into equations and solves them. The guiding constraint: every
equation couples a node only to its neighbours, so the Jacobian is **block-tridiagonal** and
each Newton iteration is a single block-Thomas sweep. That constraint is what makes
interactive speeds possible, so non-local physics is out of scope.

## The equations

In one dimension, inside each region, every species $`i`$ present there obeys a balance, a flux
law and its material's statistics, and the electrostatic potential obeys Poisson's equation:

```math
\begin{aligned}
&\frac{\partial c_i}{\partial t} + \frac{\partial N_i}{\partial x} = \sum_r \nu_{ir}\, r_r
  && \quad\text{(conservation, with bulk reactions } r\text{)} \\
&N_i = -\frac{D_i c_i}{RT}\frac{\partial\bar\mu_i}{\partial x} + v\, c_i
  - \frac{D_{\mathrm{mix}}}{RT}\Big(P\,\frac{\partial\bar\mu}{\partial x}\Big)_i
  && \quad\text{(drift and diffusion, advection, eddy mixing)} \\
&\bar\mu_i = \mu^\circ_i + z_i F\phi + RT\zeta_i, \qquad c_i = c_i(\zeta)
  && \quad\text{(the material's statistics)} \\
&-\frac{\partial}{\partial x}\Big(\varepsilon\frac{\partial\phi}{\partial x}\Big) = F\sum_i z_i c_i + \rho_{\mathrm{fixed}}
  && \quad\text{(Poisson)}
\end{aligned}
```

Each reaction runs at $`r = k_f \prod c^\nu\,(1 - e^{-A/RT})`$, its affinity $`A`$ from the reactants'
and products' $`\bar\mu`$. In a strictly neutral material ($`\varepsilon = 0`$) Poisson's equation becomes
the constraint $`F\sum_i z_i c_i + \rho_{\mathrm{fixed}} = 0`$. Some materials need less:

- A **conductor region** (a metal, a fast ion conductor) has only its carrier, with no storage in
  its bulk: $`\partial J/\partial x = 0`$ with Ohm's law $`J = -\sigma\,\partial V/\partial x`$, and no $`\phi`$.
- An **insertion host** holds its ion and carrier as one neutral combination, so $`\phi`$ cancels
  from everything there and is left undefined.
- A species absent from a material has no equations in it.

At each interface the laws are local to the face: per species, continuity of $`\bar\mu`$, a
conductance, a blocked face, or face reactions with their own rates; for $`\phi`$, an alignment
(dipole), neutrality or a Helmholtz capacitance, with the displacement continuous apart from any
sheet charge. Contacts and ports join the device to outside phases with known levels, and the
terminals to their circuits ([the device reference](device.md) has every option). The rest of
this page is how these become a block-tridiagonal Newton system.

## Unknowns

At a node of an ordinary region (a dielectric, a semiconductor, a resolved electrolyte), the
unknowns are:

- $`\hat\phi = F\phi/RT`$, the dimensionless (bookkeeping) electrostatic potential;
- $`\eta_i = \bar\mu_i/RT`$, the dimensionless electrochemical potential of each species present.

Each block holds only the unknowns that exist at its node, up to $`1 + (\text{number of species})`$.
A conductor region's nodes hold its carrier's $`\eta`$ and the flux through its single segment
([metal regions](#metal-regions)); an insertion host has no $`\phi`$ ([statistics](#statistics)), nor
does a charge-free region that nothing fixes it in (an oxide between silicon and a gate keeps its
$`\phi`$); a species absent from a material has no $`\eta`$ there; and each
interface adds a block of its own, for what crosses it ([interfaces](#interfaces-doubled-nodes-plus-a-flux-node)).

Concentrations follow from these through each material's statistics, $`c = c(\zeta)`$ with
$`\zeta_i = \eta_i - \mu^\circ_i/RT - z_i \hat\phi`$ (ideal: $`c_i = c_{\mathrm{ref},i} \cdot e^{\zeta_i}`$;
see [Statistics](#statistics) below). Working in $`\bar\mu`$ rather than $`c`$ is a deliberate
trade-off:

- **For $`\bar\mu`$:** concentrations are positive automatically, with no clamping (which would itself
  break conservation). Equilibrium is exactly flat $`\bar\mu`$. Contacts, alignments, affinities and
  Butler–Volmer are all natural in $`\bar\mu`$.
- **Against $`\bar\mu`$:**
  - Storage is nonlinear in the unknowns, so conservation holds to Newton's tolerance rather
    than trivially (see [Conservation](#conservation)).
  - A species' equations scale with its own concentration, so where it's scarce (minority
    carriers), its $`\bar\mu`$ is weakly determined and Newton can swing it wildly. Damping and time-step
    control deal with that.

Concentration variables (the usual TCAD choice) would conserve every iteration exactly, but
their positivity fixes are fiddly or non-conservative. Slotboom variables $`\exp(\eta)`$ overflow at
modest voltages.

### Extra precision where it matters: compensated μ̄

A majority carrier carrying a small current needs a quasi-Fermi step between neighbouring nodes
far below the resolution of $`\eta`$ itself. For example, holes in 1e16 cm⁻³ p-silicon carrying a
µA/m² reverse current need $`\Delta\eta \approx`$ 3e-19 across a cell, while $`\eta \approx 20`$ has an
ulp of 3.5e-15. In plain doubles that current would come out as exactly zero. (Concentration variables
hit the same wall, as cancellation between huge drift and diffusion terms.)

So every unknown is stored as a compensated double-double, `u + uLo`:

- Only *differences* of unknowns use the low part:
  $`\Delta\eta = (\mathrm{hi}_R - \mathrm{hi}_L) + (\mathrm{lo}_R - \mathrm{lo}_L)`$. Neighbouring
  values are close, so $`\mathrm{hi}_R - \mathrm{hi}_L`$ is exact (Sterbenz), and $`\Delta\eta`$ is
  recovered to ~1e-30 absolute.
- Updates use an error-free two-sum.
- The Jacobian and linear solve stay in plain doubles. Newton's fixed point is set by the
  residual, not the Jacobian, so an inexact Jacobian only slows convergence. This is the
  standard mixed-precision idea: residual in high precision, correction in low.

The cost is a few flops per unknown per iteration, lost in the noise. With it, left and right
contact currents agree to ~1e-10 down to µA/m² leakage; without it they disagreed by 3%.

## Grid and geometry

Nodes sit at positions $`x_k`$, and segments join consecutive nodes. Each node owns a box reaching
halfway to each neighbour (a vertex-centred finite volume, the "box method"):
$`\mathrm{vol}_k = (h_{k-1} + h_k)/2`$.

Every region boundary is a **doubled node**: two nodes at the same $`x`$, one per region, each
owning only a half-box on its own side. So no segment ever straddles a material step, and every
box is entirely one material.

```
region A                                   region B
 o-------o------o----o--o-o-o|o-o-o--o----o------o-------o
                             ^ doubled node: L and R at the same x
```

Grids grade geometrically from `hmin` at each region end up to `hmax`. They're never refined
automatically, to avoid performance surprises. Instead, solutions warn when a double layer the
model is meant to resolve is coarser than the local Debye length.

### A cross-section

With a [`geometry`](device.md#geometry), the device has a cross-section $`A(x)`$ (spherical
shells $`4\pi r^2`$, cylindrical $`2\pi r`$ with $`r = r_0 + x`$, or a profile), and every balance
is for totals through it:
$`\partial_t (A c_i) + \partial_x (A N_i) = A \sum_r \nu_{ir} r_r`$, and Gauss's law
$`-\partial_x (A\varepsilon\,\partial_x\phi) = A\rho`$. The box method carries this without writing out any
$`1/r`$ terms: a box's volume becomes $`\int A\,dx`$ over it (exactly, for each shape), and what a
segment, a face or a contact passes is multiplied by the area it passes through. A segment's
$`1/h`$ becomes $`1/\int dx/A`$ over it, its resistance per unit conductivity, which is the exact
weighting for steady diffusion (a constant total flux $`AN`$, so $`\Delta c = AN\int dx/A\,/D`$).
Scharfetter–Gummel keeps its form with that weighting, and so does the displacement. Where $`A`$
vanishes at an end (a sphere's centre), $`\int dx/A`$ diverges, so that segment takes $`A`$ at its
middle. A flux node's unknowns stay per area, entering each side's balance times the face's
area. A planar device has $`A = 1`$ and every factor is exactly 1, so it's unchanged to the last bit.

### What lives where

| Place | Holds |
|---|---|
| **Node** (box centre) | the unknowns $`\hat\phi`$, $`\eta_i`$, and hence $`c_i`$. The node's region supplies $`\mu^\circ`$, $`c_{\mathrm{ref}}`$, presence and fixed charge. Everything volumetric, times the box volume: storage $`\mathrm{vol}\cdot(c - c_{\mathrm{old}})/dt`$, space charge $`\mathrm{vol}\cdot(F\sum z c + \rho_{\mathrm{fixed}})`$, bulk reaction rates $`\mathrm{vol}\cdot\nu\cdot r`$. One Poisson row and one balance row per species. |
| **Segment** (always within one material) | particle fluxes $`N_i`$ (Scharfetter–Gummel, with the segment's $`D`$) and displacement $`D = -\varepsilon\Delta\phi/h`$ (the segment's $`\varepsilon`$). Each is computed **once** and added with + to the left box and − to the right box, which is what makes conservation telescope exactly. |
| **Interface** (zero width) | a zero-volume **flux node** between the doubled nodes (below). Its unknowns are what crosses ($`D_f`$ and each $`N_{f,i}`$); its rows are the surface's laws. |
| **Contact** | conditions on the outer face of the end node's half-box (below). |

Output arrays are indexed by node, so each interface $`x`$ appears twice and steps plot as true
vertical lines. Nothing is averaged across a material step: no harmonic means of $`D`$ or $`\varepsilon`$
across an interface, because no segment crosses one.

## Fluxes

Within a segment (constant $`D`$, $`\mu^\circ`$, $`\varepsilon`$), the Scharfetter–Gummel flux is the exact
constant-flux solution for a linear potential:

```math
N = \frac{D}{h} \left[ B(\Delta)\, c_L - B(-\Delta)\, c_R \right], \quad \Delta = z (\hat\phi_R - \hat\phi_L), \quad B(x) = \frac{x}{e^x - 1}
```

It's evaluated in an algebraically identical, better-conditioned form. Using
$`B(-\Delta) = B(\Delta)e^\Delta`$ and $`c_R e^\Delta = c_L e^{\Delta\eta}`$:

```math
N = -\frac{D}{h} \cdot B(\Delta) \cdot c_L \cdot \mathrm{expm1}(\eta_R - \eta_L)
```

This is exactly zero at equilibrium (flat $`\bar\mu`$), and it's precise relative to the quasi-Fermi
difference instead of being a cancellation of two large terms. $`B`$ and $`B'`$ are evaluated from a
Taylor series near 0 and with `expm1` elsewhere, finite for all arguments.

**Advection** at velocity $`v`$ adds to the drift: in the segment's Scharfetter–Gummel flux, the
potential and the quasi-Fermi difference are both shifted by the cell Péclet number
$`\mathrm{Pe} = v h/D`$, $`\Delta \to \Delta - \mathrm{Pe}`$, $`\Delta\eta \to \Delta\eta - \mathrm{Pe}`$.
That's the exact constant-flux solution of $`N = v c - D(c' + c \psi')`$ across the cell, so pure
convection–diffusion is exact on any grid, upwinding itself automatically at large $`\mathrm{Pe}`$.

**Eddy mixing** is discretised on each segment as
$`N_i = -(D_{\mathrm{mix}}/h) \sum_j \bar P_{ij} \Delta\eta_j`$, where $`\bar P`$ is built from the
logarithmic mean of each concentration across the segment. The logarithmic mean makes a neutral
species' mixing flux exactly $`-D_{\mathrm{mix}} \Delta c/h`$. The flux is exactly zero at
equilibrium, and exactly current-free ($`z^\mathsf{T}\bar P = 0`$ for any $`\bar P`$ built this way). Its
full Jacobian, including the dependence of $`\bar P`$ on both ends, keeps Newton quadratic.

With non-ideal statistics, $`\ln c = \zeta - \mathrm{ex}`$, where the excess
$`\mathrm{ex} = \zeta - \ln(c/c_{\mathrm{ref}})`$ is zero for ideal statistics. The excess enters
exactly like an extra potential, taken linear along the segment as $`\phi`$ is, so the same formula
holds with $`\Delta = z(\hat\phi_R - \hat\phi_L) + \mathrm{ex}_R - \mathrm{ex}_L`$. That's the
excess-chemical-potential generalisation of Scharfetter–Gummel: still exactly zero at
equilibrium, since the expm1 factor is untouched. It's exact for a single diffusing species on
a lattice (whose excess is proportional to its grand potential, linear in $`x`$ at steady state)
and second order otherwise. Its Jacobian couples each flux to every species at both nodes through
$`K`$, which fills the blocks but keeps them tridiagonal.

## Interfaces: doubled nodes plus a flux node

In the linear system, a zero-volume flux node sits between the two sides' nodes:

```
 … [node L] — [flux node f] — [node R] …
    φ̂_L, η_L     D_f, N_f,i      φ̂_R, η_R
```

Node L's balance rows treat $`D_f`$ and $`N_{f,i}`$ as outflow, and node R's treat them as inflow. The
flux node's own rows are the interface laws:

| Law | Row |
|---|---|
| $`\phi`$, `pinned` | $`\hat\phi_R - \hat\phi_L = \mathtt{dipole}/V_T`$ |
| $`\phi`$, `neutral` | $`D_f = 0`$ |
| $`\phi`$, `capacitive` | $`D_f + C \cdot V_T \cdot (\hat\phi_R - \hat\phi_L - \mathtt{dipole}/V_T) = 0`$ |
| species, equilibrium | $`\eta_R - \eta_L = 0`$ |
| species, blocked or absent on one side | $`N_f = 0`$ |
| species, conductance | $`N_f - G V_T (\eta_L - \eta_R)/(z^2 F) = 0`$ |
| face reaction $`k`$ | $`r_k - \mathrm{rate}_k(\mathrm{state}_L, \mathrm{state}_R) = 0`$, and each participant's edge node takes $`\nu \cdot r_k`$ |
| gate $`q`$ | $`(x_q - x_q^{\mathrm{old}})/dt - \alpha(V)(1 - x_q) + \beta(V)\,x_q = 0`$, $`V = V_T(\hat\phi_R - \hat\phi_L)`$; a gated permeability's $`P`$ takes $`\prod x^p`$ |

Every row couples only neighbours, with no penalty terms. A face reaction's rate is an unknown
of the face block, so the two edge nodes it couples meet only through the block between them;
so is a gate's fraction open, its row coupling it to the two edge nodes' $`\hat\phi`$, and the gated
fluxes to it, all within the face's block and its neighbours'. Newton's convergence and the time
step's error count a gate's change as they count a potential's in thermal units (1e-3 of $`m`$ is
0.3% of $`m^3`$), and a gate is evaluated in [0, 1] where Newton's iterates stray past its ends. A
membrane port's gates are unknowns of each window node's block, after its surface's coverages,
with the same row in the voltage across the port's capacitance there.
A species can take part in a reaction at a face and also cross it by a link (its flux slot)
at the same time; where it exists on both sides, that link must be given explicitly. A continuity
row has a zero diagonal block, but block Thomas still sees a non-singular block there once the
preceding elimination has run. Interface fluxes come out as unknowns, which is what the
solution reports per interface.

(An earlier design kept only the doubled nodes and put the flux in the right node's slot, with
the value read from the left node. That breaks tridiagonality, because the right node's other
segment then reaches two nodes back.)

A **sheet charge** at a face is booked in node R's Poisson row. With a `pinned` law that's
immaterial: the jump is fixed either way, and only the reported $`D_f`$ shifts by $`\sigma`$. With a
capacitive law it would matter, since the charge then sits on one plate or the other.

## Contacts

The end node's balance rows are complete except for the flux through the outer face. So
**before** a row is replaced by a contact condition, its residual *is* that flux (or
displacement). Contact fluxes and terminal currents are recovered exactly this way, with no
extra unknowns. Then:

- **Species in equilibrium with the outside phase:** the balance row becomes the Dirichlet row
  $`\eta = \mathtt{target}`$ (the known outside level).
- **$`\phi`$ `bulk`:** the Poisson row becomes local neutrality (total charge, mobile plus fixed, is
  zero). The recorded residual is the outside's surface charge.
- **$`\phi`$ `pinned`:** the Poisson row becomes the Dirichlet row
  $`\phi_{\mathrm{edge}} = V - \mathtt{zeroCharge}`$, again with the residual recorded as the outside's
  charge.
- **$`\phi`$ `capacitive`:** a term $`C \cdot ((V - \mathtt{zeroCharge}) - \phi_{\mathrm{edge}})`$ joins the
  Poisson row.
- **$`\phi`$ `neutral`:** nothing is added ($`D = 0`$ at the face).
- **Conductance links** add their fluxes to the balance rows. (Contacts carry no reactions:
  an electrode is a conductor region with reactions at its face.)

A left-end "contact flux node" would have a singular first block for block Thomas, which is why
contact fluxes are read from residuals instead. After each step these readouts (contact and port
fluxes, contact displacements, the last segment's current) are evaluated at the converged state
from only the boxes they come from, the two end nodes and the port windows with the segments and
faces that touch them, in the same order as a full assembly, so they're identical to one.

### Terminals

Every terminal (a contact, or a port) has a voltage $`V`$, held by its source or floating, and a
current into the device, $`I`$. Assembly records, for each:
- **$`B = \partial\,\mathrm{res}/\partial V`$**, a column: where its outside levels enter (Dirichlet
  rows, conductance links, a gate's or pinned $`\phi`$ law);
- **$`I`$ and $`C = \partial I/\partial x`$**, a row. A contact's current is read from its end box
  before the contact's own terms go in: the box's balance residuals are what must come through
  the contact (as for the flux readouts), plus the displacement through its $`\phi`$ law,
  $`(D_{\mathrm{in}} - D_{\mathrm{in},0})/dt`$. A port's is its sources over its window, with a held
  level's read from the row it replaces.

A floating terminal's voltage is an extra unknown, with its circuit law as an extra row:
$`I - I_{\mathrm{set}} = 0`$ (driven by a current) or $`I - (V_{\mathrm{src}} - V)/R = 0`$ (behind a
resistance). These don't fit the block-tridiagonal matrix $`T`$ (a port couples to its whole
window), so they're solved by bordering, together with a steady solve's dense rows, the
conserved amounts' and islands' (see [the steady system](#the-steady-system)): with
$`y = T^{-1}\,\mathrm{rhs}`$, $`X_k = T^{-1} B_k`$ and $`Q_q = T^{-1} e_q`$ (the response to a unit
pin; for a flat stretch's amount, to its shift), the update is
$`\delta = y + \sum Q_q \mu_q - \sum X_k\, \delta V_k`$, and the $`\mu`$ (pins) and $`\delta V`$ (terminals)
come from a small dense system of the extra rows. Each costs one more back-substitution per
Newton iteration.

A terminal driven by a current that alone feeds a conserved combination (a host's electrons
behind an open circuit, a floating gate's metal) can read its current another way in a time
step: as the combination's change, $`F\,(S - S^*)/(\lambda\, dt)`$, which the step's balance rows
sum to exactly (λ is its weight over the terminal's charge, $`w/z`$). Both readings are exact;
they differ in round-off, and in what Newton can see. Through a metal, the flux is $`G\,\Delta\eta`$
with its two $`\eta`$ all but equal: a floating gate's charging current fell below that at steps
of ~0.03 s, and the bordered system came out exactly singular (a gate at rest took 5000 steps
per 100 s, half rejected, where a held one takes 44). Over long steps the flux lets a closed host
drift from its charge (3 mV over 5e11 s, and steps failed by 4e12 s), where the change keeps
it. Over short ones the change is the amount's round-off over $`dt`$ (a cathode's 0.04 mol/m² of
electrons, 1e-6 A/m² at 1 µs), and the flux, its levels carried in two words, is far cleaner. So
a step reads the change where its round-off, $`\varepsilon F S/(\lambda\, dt)`$, is below 1e-6 of
the flux's estimate, $`\varepsilon\,\max|C|`$, and a step read by the flux that fails is tried once
more read by the change.

Sources are read at the end of a step (implicit), or at the present time in a steady solve.
`advance()` lands on every waveform breakpoint and restarts its order there.

## Statistics

Each node evaluates its material's statistics once per Newton iteration: $`c(\zeta)`$, the Jacobian
$`K = \partial c/\partial\zeta`$ ($`n \times n`$, symmetric positive definite), and the excess. Every
$`c`$-dependent term then differentiates through $`K`$:

- storage and space charge: $`\partial c_i/\partial\eta_j = K_{ij}`$,
  $`\partial c_i/\partial\hat\phi = -(Kz)_i`$;
- reaction prefactors: $`\partial \ln c_i/\partial\eta_j = K_{ij}/c_i`$;
- the $`\phi`$ row's response, and the Debye length in resolution warnings, use $`z^\mathsf{T} K z`$, the
  charge capacitance, in place of $`\sum z^2 c`$.

Nodes of ideal materials keep a separate fast path with $`K = \mathrm{diag}(c)`$ implied.

Model evaluation:

- **Explicit models** (Fermi–Dirac, lattice gas, insertion) evaluate $`c(\zeta)`$ in closed form.
- **Implicit models** solve a small problem per node: Redlich–Kister and tabulated OCVs by
  safeguarded Newton in $`\mathrm{logit}(x)`$; Debye–Hückel by Newton in $`\ln c`$, whose excess
  Hessian is rank one, so each step is a Sherman–Morrison update.
- **Fermi–Dirac integrals** $`\mathcal{F}_{1/2}`$ and $`\mathcal{F}_{-1/2}`$ use:
  - an alternating series below $`x = -2`$;
  - Sommerfeld's expansion above 50;
  - piecewise Chebyshev fits (degree 23, five intervals) in between, built on first use
    (~4 ms) from a trapezoid quadrature in $`u = \sqrt{t}`$.

  Relative error is ~6e-15.

Where a *concentration* is prescribed rather than a potential (bath compositions, spectators'
initial amounts), each model inverts itself for the $`\zeta`$ of those species, holding the others.
The lattice gas does this in closed form, and the others with the same scalar or small Newton
iterations.

An insertion host's species depend on $`\phi`$ only through their neutral combination, where it
cancels, so its $`\phi`$ rows are identity rows ($`\phi`$ undefined, like a region with no charged
species). The current-continuity constraint between ion and carrier comes from their balance rows.

## Metal regions

A metal node's only unknown is the carrier's $`\eta`$, its Fermi level. Its $`\phi`$ slot is free, and
it carries the carrier flux $`J`$ through the segment to the node's right instead. Each metal
segment then has two rows:

```math
\begin{aligned}
&\eta_R - \eta_L + J/g = 0 && \quad\text{(Ohm's law, } g = \sigma RT/(z^2F^2h)\text{)} \\
&\dots + J\ (\text{out of } L),\ -J\ (\text{into } R) && \quad\text{in the carrier balances}
\end{aligned}
```

This mixed form matters. A floating metal, held only by weak conductances $`G`$ such as its
electrode reactions, has its overall level set by $`G`$ alone, while every entry of an $`\eta`$-only
Jacobian is of order $`g`$. Eliminating it then computes the last pivot as $`g - g^2/(g + G) \approx G`$
with an absolute error of order $`\varepsilon \cdot g`$: digits are lost in proportion to $`g/G`$, and
beyond $`g/G \approx`$ 1e15 the pivot is noise or zero. ($`G`$ is already lost when the diagonal
$`g + G`$ is stored.) In mixed form the matrix holds $`1/g`$ instead, every coefficient is O(1), and
the limit $`\sigma \to \infty`$ is simply the constraint $`\eta_R = \eta_L`$. A
Grassmann–Taksar–Heyman-style elimination, which carries each row's leakage separately so that
pivots are built as sums, would also keep $`\eta`$ alone exact, but it needs to know which rows form
a conductance network.

A metal region is a single cell. Nothing is stored in its bulk, so $`\eta`$ is linear across it in
transients too, and one segment is exact. That also keeps $`g = \sigma RT/(z^2F^2L)`$ as small as the
region allows.

A metal's bulk holds no charge. Its surface charge at a capacitive face is the face's
displacement, booked as a sheet of excess carriers in the edge node's half-box. That gives the
edge node a storage term in $`D_f`$ (the neighbouring block), so transients, conservation and
impedance all see the charging current. A `pinned` law isn't offered at a metal face.
With nothing on the metal side depending on $`D_f`$, the steady system's flux block would be
singular to left-to-right elimination, and the charging current would need the charge read two
blocks away.

Electrode reactions at a conductor's face take its carrier at its edge node, like any face
reaction. Whether a stretch is fed is decided by its conserved combinations: weightings $`w`$ of
the stretches' amounts that no reaction changes ($`w \cdot \nu = 0`$ for every face and bulk
reaction) and that nothing outside feeds ($`w = 0`$ on stretches reached by a contact or a port),
the null space of that stoichiometry. A stretch in none of them is fed. So Ag⁺ between silver
electrodes is fed (Ag⁺ + e⁻ ⇌ Ag(s), the electrons fed by the contacts), and the steady equations
are solved directly; Fe³⁺ and Fe²⁺ between platinum electrodes aren't, since Fe³⁺ + e⁻ ⇌ Fe²⁺
conserves the iron whatever the electrons do, so their total is kept by huge steps (below). An
electrode surface's species are columns too, holding $`\Gamma \sum v\,a\,\theta`$, so species fed only
through one (A⁺ + e⁻ = S, S + e⁻ = B⁻) conserve their total with what the surface holds. A
floating conductor starts uncharged, with its carrier's level in equilibrium with the first
reaction on its left face that takes it.

## Internal ports

A port adds a source per volume to the balance rows of the nodes in its window, after every
other term at those nodes and before the contacts. A held (`'equilibrium'`) level replaces the
balance row with a Dirichlet row, after every port's other terms (another port's electrode can
react in the same window). The port's flux is then that row's residual, read just
before replacement, exactly as at a contact. Conductance and exchange links are linear in
$`(\bar\mu_{\mathrm{out}} - \bar\mu)`$. Everything stays on the node's own block. A contact's flux
readout at an end node already includes any port source there, so contact and port fluxes always
balance. A stretch reached by a port counts as fed, and its conservation intake includes the
port's flux. On a metal, a port spans both nodes of its single cell, with $`G`$ per area spread as
$`G/L`$ per volume, so each node gets half.

## Bulk reactions

$`r = k_f \prod c_R^\nu \cdot (-\mathrm{expm1}(-a))`$, with $`a = A/RT`$ computed from the
compensated $`\eta`$. That's mass action with the reverse rate implied by the standard potentials:
exactly zero at $`A = 0`$, and free of cancellation near equilibrium. Sources enter each balance as
$`\nu \cdot r`$, with one $`r`$ per reaction per node, so every moiety (a combination that no reaction
changes) telescopes exactly.

## Assembly

Each region is assembled by the kernel for its kind, over all its nodes and segments at once:
- **conductor:** its carrier's balance, and Ohm's law in mixed form (below);
- **dilute** (ideal statistics, the fast path): storage, space charge and Scharfetter–Gummel
  fluxes, written straight into the blocks by local index; a dielectric is the case with no
  species;
- **concentrated** (any statistics): the same through $`K = \partial c/\partial\zeta`$ and the
  excess potential.

Faces, ports and contacts follow. Every kernel writes only the unknowns its region has.

## Linear algebra

Block Thomas with partial pivoting *within* each diagonal block (zero diagonal entries are
fine) but not *between* blocks. Rows are equilibrated (each divided by its largest Jacobian
entry) before factorisation. Storage is flat `Float64Array`s, and factor and solve allocate
nothing. For 300 nodes × 7 unknowns, one factor and solve takes ~0.4 ms in Node 22.

Block sizes differ from node to node: each block holds only the unknowns that exist there. A
node has $`M = 1 + n`$ slots ($`\hat\phi`$ and each species' $`\eta`$), but these aren't unknowns:
- a species absent from its material;
- $`\phi`$ where it's undefined;
- a conductor node's other species, and the flux slot of its last node;
- a blocked interface flux, or the displacement of a `'neutral'` face;
- a floating terminal's spare slots.

The state vector keeps every slot (one that isn't an unknown keeps its value), while the
residual, the update and the Jacobian are compact: assembly writes each entry straight into
blocks of $`m_b \times m_b`$, $`m_b \times m_{b-1}`$ and $`m_b \times m_{b+1}`$, through each slot's row
within its block, and drops terms in slots that aren't unknowns. Elimination costs $`\sum m_b^3`$
instead of $`n \cdot M^3`$, which pays wherever species are confined to some regions: for n-Si
against KCl, three unknowns of five at every node, it's a fifth of the work. A finite-difference
test checks the compact Jacobian column by column on devices that cover every assembly path.

### Exact J·v

The assembled matrix loses what a flux's two entries share. A Scharfetter–Gummel flux depends on
$`\eta`$ at its two ends through $`G`$ and $`-G(1 + E)`$, with $`E = \operatorname{expm1}(\Delta\eta)`$;
for a nearly uniform $`v`$ they cancel, and $`J \cdot v`$ keeps only round-off of $`G|v|`$, which in
an inversion layer's 0.1 nm cells ($`G`$ ~ 1e11) dwarfs the flux. Where a solve needs
$`J \cdot v`$ itself (GMRES, in the impedance and in Newton's refined solves), the dilute kernels'
terms are kept aside instead, each evaluated on differences first:

```math
\delta N = a\,(\delta\eta_L - z\,\delta\hat\phi_L) + b\,(\delta\eta_R - \delta\eta_L) + d\,(\delta\hat\phi_R - \delta\hat\phi_L),
```

with $`a = -GE`$ (through $`c_L`$), $`b = -G(1 + E)`$ and $`d`$ the drift's. A node's storage and space
charge depend on $`\delta\eta_i - z_i\,\delta\hat\phi`$ alone, the displacement on
$`\delta\hat\phi_R - \delta\hat\phi_L`$. Everything else (faces, contacts, ports, reactions,
concentrated and metal regions) is assembled into a matrix of its own. Rows that a contact or a
steady solve replaces drop the kept terms; the transformed rows of strictly neutral nodes combine
them as they combine the matrix's. A region's uniform shift ($`z_i s`$ on each $`\eta_i`$, $`s`$ on
$`\hat\phi`$) then gives exactly zero inside it, as it must. It replaced a central difference of
the residual, whose truncation (~1e-9 of $`|J||v|`$) swamped the answer for such a shift, and
which cost four assemblies per product.

## Newton

- The update of the device's potentials ($`\hat\phi`$, $`\eta`$) is limited to 10 thermal units per
  iteration by uniform scaling, which scales a floating terminal voltage's update too. That
  voltage doesn't count toward the limit: it enters only linearly (conductance links, held
  levels, a capacitive face), so a large swing in it is safe. When it did count, a port driven by
  a current pulse, its voltage collapsing by hundreds of volts as the current switched off,
  dragged every unknown along 10 thermal units at a time (and a swing past 1e4 tripped the
  divergence check): the strong-injection Haynes–Shockley benchmark took 55% more
  factorisations.
- Converged when a full, undamped update is below 1e-10 (thermal units). Quadratic
  convergence makes the remaining residual negligible.
- Or converged as far as round-off allows: updates already below 1e-6 thermal units (~26 nV)
  that have stopped shrinking for two iterations. A badly conditioned system's round-off floor
  can sit above 1e-10: a strictly neutral material on a short step, where $`\phi`$ is fixed only
  through fluxes that the storage term dwarfs (a condition number of about $`h^2/(D \cdot dt)`$,
  much reduced by the change of variables below), converges quadratically to ~1e-9 and then
  rattles there. Short steps come right after every waveform breakpoint, so without this a
  cyclic voltammogram in a neutral electrolyte stalled at its turns.
- **Refined solves, where the factorisation loses a mode.** The assembled $`J`$ holds a flux's
  dependence on $`\eta`$ at its two ends as two entries; where they're huge, eliminating them loses
  what they share to round-off. In a MOS capacitor's inversion layer the electrons' cells, 0.1 nm
  wide, are coupled about $`h^2/(D\,dt)`$ ~ 1e-18 times more tightly than they store charge over a
  second, so the layer's overall level, set by its storage and the trickle of minority electrons
  from the bulk, is lost, and Newton's update there is noise. When Newton stalls short of
  convergence (updates below 1e-4 that stopped shrinking quadratically), or its undamped updates
  grow twice running (a bipolar stack's floating base, whose holes are held ~1e14 more weakly
  than they move within it, kept Newton wandering a thermal unit at a time), the solve is refined:
  GMRES on $`J\delta = r`$, the factorised system as the preconditioner and $`J \cdot v`$ exact, each
  flux on its $`\eta`$ difference (see [exact J·v](#exact-jv)), as the impedance does. A lost mode
  takes GMRES a couple of iterations. The refined update is kept where GMRES converged and it
  moves no level by a thermal unit or more; otherwise refining stops for the rest of the solve.
  A MOS capacitor's gate step now reaches the low-frequency charge in ~70 steps (it took 80,000,
  and ended 1.5% off). Once refined, a solve stays refined: a GaAs junction without
  recombination, whose minority carriers reach it only from the contacts, passes ~1.5e-7 A/m²,
  which plain solves past the first refinement left to round-off (cold and warm solves 1e-3
  apart), and refined ones give to 1e-14. (When $`J \cdot v`$ was a central difference of the
  residual, its noise made refining past the first lost mode chase round-off, so it was kept
  only where its correction dwarfed the plain update.)
- Clear divergence (device updates beyond 1e4 after the first iteration, or ten times the
  first after six) bails out early, so the caller can take a smaller step. The first update is
  exempt: light flooding a population that starts at ~1e3 per m³ predicts, linearised in
  $`\eta`$, a change of 1e10 thermal units, which damping covers in a few iterations; bailing
  out there had a lit phototransistor's pseudo-transient fail at every step above 1e-21 s.
- A pivot that cancels to exactly zero in the block factorisation is perturbed to 1e-15 of its
  block's largest entry (static pivoting) rather than failing the solve. GaAs minority carriers,
  ~1 per m³ beside a face, are held through conductances 1e-16 of the face flux's coefficient in
  their row, and the elimination cancels them exactly: such devices now solve, even at
  equilibrium. Newton's updates stay meaningful (the factorisation only amplifies the residual,
  so a small update still means a small residual), and its refined solves correct the rest.
- When a steady solve fails, the solution's warnings say how nearly singular the system was,
  and where. A running error bound through the factorisation compares each pivot with the
  magnitudes it was formed from, $`\log_{10}(\sum|\mathrm{terms}| / |\mathrm{pivot}|)`$: the digits
  lost. A part of the device held only weakly, such as a floating region coupled through tiny
  conductances, loses about $`\log_{10}(g/G)`$ and fails near 15. It's not reported for solves that
  converge, since digits can be lost harmlessly too (a saturated species whose $`\eta`$ nothing
  depends on).
- The last iteration usually only confirms convergence (updates go like 4e-3, 7e-6, 3e-11).
  Stopping one iteration early on a convergence-rate estimate would save a third of the work,
  but leave ~1e-11 in $`\eta`$ per step, which shows up as conservation drift, so it isn't done.
- Tried and rejected for large jumps (the retry at 3 thermal units per iteration takes 30–50
  iterations for a 1.4 V jump on a pn diode): per-component logarithmic damping, as in
  SPICE junction limiting, and an adaptive limit that grows while updates shrink. Both
  diverged more often and cost more overall.

## Time stepping

- **Backward Euler** is the default for `step(dt)`. It recursively halves the interval wherever
  Newton fails, always landing exactly on $`t + dt`$. Minority carriers that must rise by many
  orders of magnitude after a bias step are the typical cause.
- **Contact displacement** at the start of a step is carried over from the end of the previous
  one, including across `set()`, so a step change in a gate voltage shows up as displacement
  current.
- **BDF2** (`method: 'bdf2'`) is the variable-step second-order backward difference formula.
  With $`\omega = dt/dt_{\mathrm{prev}}`$, its storage term is
  $`[a_0 c - (1+\omega) c_n + \omega^2/(1+\omega)\, c_{n-1}]/dt`$, where $`a_0 = (1+2\omega)/(1+\omega)`$.
  - It's written in backward-Euler form, $`(c - c^*)/(dt/a_0)`$ with the history reference
    $`c^* = [(1+\omega) c_n - \omega^2/(1+\omega)\, c_{n-1}]/a_0`$. So assembly is shared, and the
    displacement histories (for terminal currents) take the same combination.
  - Conservation still telescopes. Closed amounts stay exact, since $`c^*`$ has the same amount
    as $`c_n`$ when $`c_n`$ and $`c_{n-1}`$ did. Open stretches add the history term
    $`\sum v (c^* - c_n)`$ to their intake.
  - A step ratio above 2 falls back to backward Euler for that step, to stay safely within
    variable-step BDF2's zero-stability limit ($`1 + \sqrt{2}`$).
- **Adaptive stepping** (`advance`):
  - The local error of each step is estimated against an explicit predictor through the
    previous states: quadratic after a BDF2 step, scaled by $`C_c/(C_c + C_p)`$, with
    $`C_c = h^3(1+\omega)^2/(\omega(1+2\omega))`$ and $`C_p = h(h+h_1)(h+h_1+h_2)`$; linear after
    backward Euler, scaled by $`h/(2h + h_1)`$. The first step, with no history, is checked by step
    doubling.
  - A first step whose error doesn't shrink with it is taken anyway. A jump in a level held at a
    boundary (a bath's, after a voltage step) starts a profile self-similar in $`x/\sqrt{t}`$, so a
    first step of any length errs alike (0.11 thermal units over ten decades of step, in a
    0.6 M electrolyte with resolved double layers), down to the grid's own diffusion time,
    where Newton runs out of digits. Three tries in a row, each within a factor 2 of the last,
    mark the error as the jump's; backward Euler damps it as it should, and the steps grow
    from there. Transients that start with a jump take 5–30% fewer factorisations for it.
  - The state a first step starts from needn't satisfy the algebraic equations ($`\hat\phi`$ in a
    strictly neutral material, an interface's unknowns), after a start or a jump; they jump in
    the first instant. So the history starts from the first step's half step, not from there,
    and the second step is backward Euler: a predictor through the starting state saw the jump
    as error that shrank only as fast as the step, and a junction from a sharp boundary at a
    tight tolerance stepped down into round-off.
  - The error is measured in thermal units over every state potential ($`\hat\phi`$ where defined,
    each present $`\eta`$, an electrode surface's coverages as their $`\eta`$; a floating terminal's
    voltage follows from these). That includes species far below their
    largest concentration, such as minority carriers ahead of a diffusion front. Weighting
    them down (a mixed relative/absolute tolerance,
    $`\delta\eta \cdot c/(c + 10^{-10} c_{\mathrm{max}})`$) took a third fewer steps on a pn turn-on,
    but the tail's errors are swept into the front, and at matched accuracy of the terminal
    current it cost 50% more.
  - A step is rejected above `tol`. The next step size is
    $`h \cdot \min(2, \max(0.2, 0.9\, (\mathtt{tol}/\mathrm{err})^{1/(p+1)}))`$.
  - Newton failure quarters the step. (Strictly neutral materials had tried longer first steps
    instead, shorter ones being worse conditioned there; with their ends and starts fixed, no
    case needs it.)
  - A target within round-off of the present time (1e-10 relative), such as an animation
    frame's that lands a hair past a breakpoint just reached, is snapped to rather than
    stepped to: a step of 1e-13 s can't be resolved.
  - Newton starts each step from the state extrapolated through the last three, which saves
    about a third of the iterations.
  - The wall-clock budget is checked between steps.

## Small-signal impedance

About a steady state, a small sinusoidal excitation at one terminal gives, to first order,

```math
\begin{aligned}
&(J + i\omega M)\, \delta x + \sum_k B_k\, \delta V_k = -B_T\, \delta V_T && \quad\text{(the grid's rows)} \\
&(C_k + i\omega C'_k)\, \delta x + (\partial I_k/\partial V_k + i\omega\dots)\, \delta V_k = \delta I_T\, \delta_{kT} && \quad\text{(each floating terminal's circuit row)}
\end{aligned}
```

- **$`J`$** is the steady Jacobian, **$`M`$** the Jacobian of the time-derivative terms (storage
  $`v \cdot K`$, displacement), and likewise $`C`$ and $`C'`$ for the terminal currents. The
  time-derivative parts are read off two assemblies, $`J(dt) = J + M/dt`$, with a tiny $`dt`$ so that
  the subtraction loses nothing.
- **At the measured terminal $`T`$,** a held voltage is perturbed ($`\delta V_T = 1`$, and $`\delta I_T`$
  read from its $`C`$ row), or a driven current ($`\delta I_T = 1`$, and $`\delta V_T`$ solved for).
  $`Z = \delta V_T/\delta I_T`$, with $`I`$ into the device. The other terminals keep their drives:
  held ones at AC ground, driven ones open.

Each frequency costs one complex block-Thomas factorisation (rows scaled by their largest
entry), plus one back-substitution per floating terminal for the bordering. Conserved (blocked)
species make $`J`$ singular, but $`J + i\omega M`$ isn't for $`\omega > 0`$: at low frequency a blocking
device looks like a capacitor, as it should.

The factorised system is used as a preconditioner, not trusted alone: each solve runs GMRES with
an exact $`J \cdot v`$ (see [exact J·v](#exact-jv)). The assembled $`J`$ holds a flux's dependence
on $`\eta_L`$ and $`\eta_R`$ as two entries; in an inversion layer both are huge (~1e11), and for a
nearly uniform $`\delta\eta`$ they cancel to round-off of order 1e-5. That round-off had given a
MOS capacitor's inversion layer, fed by minority diffusion from a bulk with ~1e3 electrons per cm³
(a time constant of minutes), an exchange path that followed the gate at 1 Hz. Where the plain
solve is already accurate, one application of the operator confirms it. GMRES judges convergence only by the true residual,
recomputed at each restart: its own running estimate can drift far from it where the
preconditioner is nearly singular (an inversion layer's exchange, lost in the assembled
$`J`$, at a frequency of µHz), and had declared solves converged that weren't. It keeps its
best iterate, never worse than the factorised solve.

The current is read where its error is smallest. At a contact it is a sum of terms that can
be far larger than itself: across a junction without recombination, between heavily doped
sides, what passes is minority carriers of ~1e8 per m³, while the contact's terms are its
majority conductances times a $`\delta\eta`$ cancelling its neighbours' to 1e-17, and their
round-off reads as a constant conductance of order 1 S/m², a silicon p⁺n⁺ junction's capacitance
9× off at 100 Hz. Between the two contacts of a device without ports, the total current
(conduction plus displacement, $`i\omega\,\delta D`$) is the same through every cut, so it is also
read across each segment of an ideal region, linearised as the Scharfetter–Gummel flux is, and
the reading with the least estimated error is taken: in that junction, the depletion region's;
in a MOS capacitor, the oxide's displacement. A reading's error is its terms' round-off, and the
solves' error in $`\delta x`$ (their residual, relative to its largest entry) through its
coefficients: a contact's reading is a flux too (a held contact's balance row is replaced by its
level, so nothing solved pins it), and in a neutral bulk the current crosses a huge conductance
on a tiny gradient.

Where even the best reading's estimated error exceeds 1% of the current, the response is split
in two: a uniform shift of each region (every $`\eta_i`$ by $`z_i s_r`$ and $`\hat\phi`$ by $`s_r`$, with
$`s_r`$ fitted to the solve so far), which changes no flux and no concentration, so its $`J \cdot e`$ is
exactly zero inside each region and is taken at the faces and ends;
and the rest, small, solved for to GMRES's floor with its error now relative to itself, and read
separately. An electrolyte at its open circuit needs it: the response to its terminal is nearly
such a shift, and its current is the slope of levels uniform to 1e-14 of the response, more than
a double holds. A redox electrode's DC conductance had come out 2–40% low; now it matches the
steady dI/dV to six digits. (Not with ports, whose exchange sees the shift.) The real part of a nearly ideal capacitor, 1e-7 of $`|Z|`$, is still
only resolved to ~1e-4 of $`|Z|`$.

## Steady state

A steady state is where a transient from the present state ends. `solve()` finds it without
stepping there, as the solution of the steady equations ($`dt = \infty`$, no storage term): slow
physics can take seconds or hours, far beyond any $`L^2/D`$ estimate (exponentially scarce
minority carriers slowly filling an inversion layer behind a Schottky contact), and stepping
would never finish. Three things make that more than one Newton solve: without the storage term,
the steady equations say nothing of what the transient conserved; some levels are held so weakly
that a factorisation loses them; and Newton needs a start near enough. `solve()` does not advance
the clock. Nor does it ask whether the steady state is stable: an axon held at a current that
makes it fire over and over has a resting point all the same, which `solve()` returns, and only
a transient leaves it.

### The steady system

Every row of a steady solve that differs from a time step's is a sum of balance rows. It's
written out in its telescoped form, never added up from the rows themselves, so what the sum
leaves is exact, however small beside the terms that cancel:

- **A conserved amount.** A combination of stretches (weights $`w`$) that no reaction changes and
  nothing outside feeds has balance rows that sum to zero identically in a steady state: the
  fluxes cancel across its stretches, and the reactions by $`w \cdot \nu = 0`$. One of those rows
  is redundant, and the combination's amount, $`\sum_k w_k \sum v c_k`$, stands in for it. That's
  what a time step's own sum, $`w \cdot (S - S_0)/dt = 0`$, says at any $`dt`$.

  The combinations are the null space of the stoichiometry (a row per reaction, over the
  stretches it touches) together with what the terminals feed, worked out once and again after a
  change of drive. A terminal held at a voltage (or behind a resistance, which holds it at steady
  state) feeds each stretch it links freely. One driven by a current feeds its charged species
  only as that current allows, in proportion to their charges. So:
  - a **spectator** (a species blocked all round) keeps its amount;
  - so does each **moiety**: the total iron of Fe³⁺, Fe²⁺ and FeCl²⁺ in a closed cell with a
    complexation reaction;
  - so does what a **terminal at no current** alone feeds. A closed battery at open circuit
    has one steady state per state of charge, and a solve that imposed only I = 0 had converged
    to another (4.145 V where its start said 3.958 V). A terminal that passes current only by
    charging a capacitance (a gate, a port's capacitance) is the same, its charge in place of an
    amount; and an electrode port's reactions draw on what the electrode holds, its charge
    $`Q/(zF)`$ a column of its own under a current, so a port at no current reducing a closed
    M⁺ (M⁺ + e⁻ = M(s)) keeps M⁺ + Q/F, or with no capacitance, M⁺ alone. At any other current such a terminal only stores what comes in, which has no steady
    state, and the solve says so at once.
  - an **immobile** combination (trap states X⁰ and X⁻ under e⁻ + X⁰ = X⁻, with D = 0) keeps it
    node by node, since nothing carries it anywhere, as do an electrode's surface species
    (adsorbates that only turn into each other keep their total on every site), its row local,
    one per node in its own block. (The null space is taken with the immobile columns first, so
    that these come out apart from the combinations that move.) A combination of immobile
    stretches with a moving one (the electrons the traps hold, $`e^- - X^0`$, in a closed film) is
    an amount like any other, its immobile parts summed over their nodes;
  - a **floating metal** keeps the charge on its faces, its sheets' displacement.

  A closed stretch keeps the amount it was given (its `c0`, or what it held when `set()`
  closed it); a combination behind a terminal at no current keeps what the state holds when the
  solve begins. The amount's row is dense, so it's bordered (see [terminals](#terminals)): the
  factorised matrix gets a pin, an identity row, in place of the balance row, which makes the
  stretch's level a well-conditioned unknown, and the response to a unit pin (one more
  back-substitution) goes in at whatever weight satisfies the amount. Behind a terminal at no
  current, the amount stands in for the terminal's circuit row instead: at no current that row is
  the sum of the combination's balance rows, so it says nothing new. (An electrode's charge, which
  moves with the terminal's voltage, goes only there.) The rows replaced must be independent:
  each combination's weights over them a nonsingular matrix, or one balance is said twice and
  another not at all. Picked one combination at a time, two complexations sharing A
  (A + B = AB, A + C = AC) had −A + B + C in place of B's row, B + AB of AB's and A − B + AC of
  A's, and the solve converged to a state that wasn't steady, depending on the order the species
  were listed in. So they're chosen by elimination.

  Holding an amount through a storage term at a huge $`dt`$ instead is badly conditioned: the
  stretch's level is then held only by $`v \cdot c/dt`$, against internal conductances
  $`D \cdot c/h`$ larger by $`D \cdot dt/(hL)`$, around 1e12 for a micron-scale cell (see
  [metal regions](#metal-regions) for the same cancellation). Newton then needs dozens of
  iterations, or fails. It's left only where a species is still in one region and moves in the
  next: there, backward-Euler steps at a huge $`dt`$ (10⁶ × the slowest diffusion time) keep
  the storage term, growing ×10 (capped) while the state still moves, each spectator's level
  first shifted uniformly to restore its amount exactly (round-off creeps through the vanishing
  storage term), and one last step at the base huge $`dt`$, where pinning is tight.
- **A flat level.** A species that no reaction or port touches, reached by one contact at most,
  carries no flux in a steady state: every sum of its balance rows from its closed end is the
  flux through one segment, zero. In an ideal region without flow or mixing that's a level flat
  across the segment, so its level is flat throughout: its contact's, or reaching none (a
  spectator), the level its amount fixes. Steady solves (not the impedance, where it carries a
  current) hold it there in place of its balances, rather than find it through its own
  conduction, which can be all but nothing. A level found that way can't be moved: a MOS
  capacitor's inversion electrons reach the back contact only through a bulk with ~1e3 of them
  per cm³, and with that contact moved, each solve failed through ~90 pseudo-transient steps;
  past a closed zinc cell's limiting current, extended space charge at the cathode excludes the
  blocked sulfate to ~1e-25 mol/m³, and the chain of conductances through it had lost the steady
  system 14 digits (steady solves failed from 0.7 V, and took 37,000 iterations at 0.6 V). (A
  level that starts flat stays exactly so either way: carrying nothing, the flux's expm1 form is
  exactly zero, whatever $`\phi`$ does.) A spectator's first node keeps its amount's row, a pin
  whose response is the whole stretch's shift, so the stretch moves as a whole. At a face inside
  the stretch whose link holds the level continuous, the face's flux, which only the replaced edge
  balances set, is pinned to zero too. (Not through flow or mixing, nor a concentrated material's
  cross-diffusion, where zero flux isn't a flat level.) Species that react only among themselves
  are held flat the same way where each reaches the same one contact and nothing else, their
  reactions balanced at its levels: one reservoir drives nothing, so in a steady state they're at
  equilibrium. A MOS capacitor's electrons and the holes they recombine with reach only its back
  contact; moved, slow generation through a bulk of GaAs with ~1e-4 electrons per cm³ couldn't
  carry the inversion layer's level with it (of 500 random MOS capacitors, 38 failed and 5
  converged to a gate charge up to 2× off).
- **An island's balance.** A piece of a stretch between its faces that nothing else feeds (no
  contact, port or reaction) holds its level only through its faces' fluxes, which can be far
  weaker than its own conduction: a face conductance of 1e-15 S/m² against a region conducting
  ~4e2 S/m, or neighbours that conduct 1e20 times less than it. Eliminated, the level was lost to
  round-off below anything Newton could see (the plain solves gave no update along it, so nothing
  stalled), and solves converged to a wrong one silently: currents 2× or 16× off. Summed over the
  piece's nodes, its balance rows are its faces' flux unknowns alone, $`A_R u_R - A_L u_L`$: the
  internal fluxes cancel exactly, by construction rather than in round-off. Across a face that
  holds $`\bar\mu`$ level, $`u`$ is set only by the edge balances, and eliminated it can come out
  of the piece's own, as $`G\,\Delta\eta`$ with its two $`\eta`$ all but equal; so there the
  outside edge node's balance is added too, $`u`$ cancels, and the outside's first segment
  carries the flux. (Seen through one face of two, a level behind neighbours 1e20 times less
  conductive was overshot twofold, and Newton swung it from one contact's level to the other's.)
  That sum replaces the first node's row, bordered like an amount. The response to the unit pin
  is the island's uniform shift, with whatever $`\hat\phi`$ and the other species do along with
  it, so no mode needs guessing; refined solves carry the pinned rows too.

  Newton exact along such a level can be worse than Newton blind to it, far from the solution
  (cold at bias, a bipolar stack's floating base took 1e5 thermal units in its first update), so
  the summed rows go in only once the plain solves converge. An island at the wrong level then
  shows: it takes in a different current than it passes on, so the terminal currents don't add
  up to zero (to 1e-9 of the largest). Newton then carries on with the summed rows, refined from
  their second update. A current that tiny against the conduction around it rides on $`\eta`$
  differences (~1e-28 for 5 fA/m² through 1e4 mol/m³ at D = 1e-4 m²/s) that the last update
  leaves known only to its own round-off, so where the currents still don't add up, one more
  update polishes them. If Newton doesn't converge again within a dozen updates (a GaAs stack's
  minority carriers, whose own rows rest on static pivots, at 1e-17 A/m² of noise), the plain
  solution stands; and a steady solve that fails outright is tried once more with the summed rows
  in from the start. Such devices now solve directly, their currents right to ~1e-15.

All three move a weakly held level as a whole, where a factorisation can't: an amount's or an
island's pin by a shift it solves for, a flat level by holding every node to one. What the
factorisation still loses elsewhere, Newton's refined solves recover ([Newton](#newton)).

### Reaching it

A solve that fails falls back stage by stage:

1. **Newton from here.** If it diverges, it's retried once with tighter damping (3 thermal units
   per iteration), enough for most large jumps, such as a cold start at forward bias. Where no
   continuation applies, the first stage goes straight on into the pseudo-transient ramp (3).
   Where the state is a transient's, it can be a worse start than a cold one (a nanosecond after
   a lit MOS capacitor's gate jumps, its surface is in deep depletion, minority carriers ten
   orders below where they end): if Newton fails from there, the solve goes on as a cold device
   would, from the layout for the drives, keeping what it took from the state (its amounts and
   charges). Nor does a continuation start from the last steady solve's voltages once a
   transient has moved the state.
2. **Continuation from a solved state, along the drives,** each step a steady solve warm from
   the last:
   - **The light.** A device with generation reactions (species made only from, or turned only
     into, fixed reservoirs, such as photogeneration from a photon reservoir) can be held far
     from equilibrium even with its terminals level, where a voltage ramp can't help; nor can a
     time step, since light flooding minority carriers 20 orders below their lit level makes even
     a step of 1e-16 s a jump Newton can't take. The generation rates are scaled down to 10⁻³⁰ and
     ramped back up, ×100 a step while each solve converges and by the square root of the factor
     when one doesn't, giving up when that factor falls below 1.5 (crawling, which time steps do
     better); the pseudo-transient ramp then starts from the dimmer light's solution rather than
     cold. (At 10⁻¹², a GaAs layer's generation still exceeded its equilibrium recombination,
     ~3e-24 mol/m³ of minority electrons, a billionfold: the first step was the whole jump.) An
     illuminated 80 µm silicon diode solves cold this way in about 80 iterations. (A lit n⁺pn
     stack with a floating base crawled for 20,000 iterations: at a trillionth of the light, the
     base's level is held by couplings of ~1e-17.)
   - **A current.** A cold steady solve at open circuit can fail where the same state is easy
     warm (a solar cell in dim light on a fine grid). The terminal is held at a voltage instead,
     its steady solves having their own continuation, and the voltage marched from the start's,
     the step growing, until the current crosses its target: up where the current is below it,
     since raising a terminal's voltage raises the current into a passive device (the other way
     if 20 V brings no crossing). Bisection narrows the crossing to 0.1 mV, and the terminal is
     floated from beside it. Floated, the system can lose what held it doesn't (a closed
     Fe³⁺/Fe²⁺ cell driven near its limit, Fe³⁺ 20 orders below Fe²⁺: 33 digits short), but the
     current-driven steady state is the held one at the voltage that passes the target, so where
     the floated solve fails, held solves find that voltage instead, by regula falsi within the
     bracket, to 1e-10 of the current. (Not a terminal at no current that keeps an amount or a
     charge: every voltage passes nothing there.)
   - **The voltages.** The right terminal's voltage ramped to its target: from the contacts'
     voltages at the last converged solve when the state is that solution (the left one first,
     where it moved), and otherwise from level terminals, where a cold start is consistent. (A
     device first solved after a `set()` was laid out for its new voltages; if the level start
     fails from there, it's laid out again, level, and retried.) The step starts at 1/8 of the
     way, grows ×1.5 on success and shrinks ×4 on failure. The ramp's own start, with the
     terminals level, ramps the light up when there is some.
3. **The pseudo-transient ramp.** $`dt`$ ramps up from a small value, through time, ending in
   the direct solve. After a failed step it grows ×2, then ×1.5 more each success, back to ×10.
   Once steps reach the slowest diffusion time, the direct solve is tried from there, once a
   decade of $`dt`$ (the state kept if it fails): with immobile traps and ions of D = 1e-17 m²/s,
   steps from 1e8 to 1e9 s kept failing on the way to the end.

Each of these carries cases that nothing else does (switched off one at a time, each fails
tests or stress cases the others pass).

## Conservation

- **Discretely exact.** Face fluxes are computed once and added ± to both neighbours, the
  interface flux unknowns do the same, and reactions are $`\nu \cdot r`$. So summed over boxes,
  everything telescopes, and the *converged* discrete solution conserves every amount and total
  charge to round-off.
- **The iterate conserves to Newton's tolerance.** Storage is nonlinear in $`\eta`$, so each step's
  drift is its summed residual. Converging to 1e-10 in $`\eta`$ makes that ~1e-15 per step in
  practice; but over steps of 1e10 s, the residual times $`dt`$ adds up, which is why a terminal
  at no current reads its current there by what it keeps ([terminals](#terminals)).
- **Reported.** Every solution's `conservation` lists each species stretch: its amount, its
  reference, the time-integrated flux through its contacts, and the resulting drift.

## Strictly neutral materials and unresolved double layers

A material with $`\varepsilon = 0`$ has local neutrality in place of Poisson, and no displacement
along its segments. $`\phi`$ there is a bookkeeping multiplier, fixed by neutrality given the
$`\bar\mu`$'s. It's reported as NaN if no charged species is present and nothing couples the region
electrostatically. The bulk $`\varepsilon \to 0`$ limit is smooth: the Poisson row simply
degenerates into neutrality. The same Newton count holds from $`\varepsilon_r = 78.5`$ down to 1e-6,
as in a linear prototype that motivated the design.

Interfaces are where it gets subtle. With a `pinned` law and double layers the grid can't
resolve (because $`\varepsilon`$ is tiny, or the cells are coarse), the charge the interface needs is
crammed into the two half-boxes beside it. The bulk stays right, but the interface compositions
and $`D_f`$ become mesh-dependent. The true $`\varepsilon \to 0`$ limit has no interface charge and a
free Donnan-type jump, which is the `neutral` law. The defaults follow from this: `neutral`
between two $`\varepsilon = 0`$ materials, `pinned` otherwise (and `pinned` between two
$`\varepsilon = 0`$ materials is an error: with no field on either side, nothing would determine its
charge). Solutions warn when a resolved-model double layer is under-resolved.

On a transient step, strictly neutral nodes are solved in better-conditioned terms. There, a
change of $`\hat\phi`$ with every $`\eta_i`$ shifted by $`z_i`$ times it leaves every concentration as it
was, so storage and neutrality don't see it. Only the fluxes do. But in $`(\hat\phi, \eta)`$ storage
and neutrality see it as pairs of huge entries that cancel, and round-off in that cancellation, of
order storage/flux ~ $`h^2/(D \cdot dt)`$ times the range of concentrations, swamps the fluxes that
fix it: a trace ion beside 3 M KCl lost about ten digits, and Newton stalled. So at each node that
stays neutral (interior nodes, edges at `neutral` faces, and a device's end where its contact is
closed to every ion and leaves $`\phi`$ alone; not the edge of a capacitive or pinned face, which
holds the face's charge, nor an end whose contact holds a level), assembled without storage and
neutrality:
- rows: the balance of the most abundant charged species (by $`z^2 c`$) becomes
  $`\sum (z_i/z_k) \times`$ each balance, which without storage is current continuity;
- columns: the unknowns become $`\hat\phi'`$ and $`\eta'_i = \eta_i - z_i \hat\phi`$ (the update is
  mapped back), in which storage and neutrality have no $`\hat\phi'`$ term at all;
- then storage (on the other balances) and neutrality go in, exactly, in those terms.

A start-of-step net charge (round-off in a solved state) isn't carried over: the neutrality row
holds the new state neutral. An end node that its contact holds (a bath, an electrode's level)
keeps its plain rows, but the contact's current there is read before the ions' storage goes in,
and a bulk law doesn't count the box's net charge as displacement: under neutrality neither
carries charge, and read from the box's balance their round-off (∝ c/dt) cancelled only roughly.
A bath floating on 1 µM beside 3 M KCl had rattled at 2e-5 thermal units on steps of 1e-14 s. Steady solves and the impedance keep the plain rows (there's no
storage at $`dt = \infty`$, and the impedance reads the storage matrix from them).

A sub-grid Gouy–Chapman law, treating the diffuse layers analytically when $`\lambda_D \ll h`$, is on
the [roadmap](../ROADMAP.md). It would handle macroscopic devices with real double-layer charge at
any grid, and tends to `neutral` as $`\varepsilon \to 0`$.

At an $`\varepsilon = 0`$ | $`\varepsilon > 0`$ face, the neutral side's boundary half-box holds the
face's charge. That slightly perturbs the carrier density there, by an amount that depends on the
mesh (0.3% on a Schottky test with the metal as an $`\varepsilon = 0`$ region). For a metal, a metal
region avoids this: its surface charge is a sheet that leaves its Fermi level untouched.
