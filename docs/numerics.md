# Numerics

How driftlet turns the physics into equations and solves them. The guiding constraint: every
equation couples a node only to its neighbours, so the Jacobian is **block-tridiagonal** and
each Newton iteration is a single block-Thomas sweep. That constraint is what makes
interactive speeds possible, so non-local physics is out of scope.

## Unknowns

At every grid node, with block size M = 1 + (number of species):

- `φ̂ = Fφ/RT`, the dimensionless (bookkeeping) electrostatic potential;
- `η_i = μ̄_i/RT`, the dimensionless electrochemical potential of each species.

Concentrations follow from these: `c_i = c_ref,i · exp(η_i − μ°_i/RT − z_i φ̂)`. Working in μ̄
rather than c is a deliberate trade-off:

- **For μ̄:** concentrations are positive automatically, with no clamping (which would itself
  break conservation). Equilibrium is exactly flat μ̄. Contacts, alignments, affinities and
  Butler–Volmer are all natural in μ̄.
- **Against μ̄:**
  - Storage is nonlinear in the unknowns, so conservation holds to Newton's tolerance rather
    than trivially (see [Conservation](#conservation)).
  - A species' equations scale with its own concentration, so where it's scarce (minority
    carriers), its μ̄ is weakly determined and Newton can swing it wildly. Damping and time-step
    control deal with that.

Concentration variables (the usual TCAD choice) would conserve every iteration exactly, but
their positivity fixes are fiddly or non-conservative. Slotboom variables `exp(η)` overflow at
modest voltages.

### Extra precision where it matters: compensated μ̄

A majority carrier carrying a small current needs a quasi-Fermi step between neighbouring nodes
far below the resolution of η itself. For example, holes in 1e16 cm⁻³ p-silicon carrying a
µA/m² reverse current need Δη ≈ 3e-19 across a cell, while η ≈ 20 has an ulp of 3.5e-15. In
plain doubles that current would come out as exactly zero. (Concentration variables hit the
same wall, as cancellation between huge drift and diffusion terms.)

So every unknown is stored as a compensated double-double, `u + uLo`:

- Only *differences* of unknowns use the low part: `Δη = (hi_R − hi_L) + (lo_R − lo_L)`.
  Neighbouring values are close, so `hi_R − hi_L` is exact (Sterbenz), and Δη is recovered to
  ~1e-30 absolute.
- Updates use an error-free two-sum.
- The Jacobian and linear solve stay in plain doubles. Newton's fixed point is set by the
  residual, not the Jacobian, so an inexact Jacobian only slows convergence. This is the
  standard mixed-precision idea: residual in high precision, correction in low.

The cost is a few flops per unknown per iteration, lost in the noise. With it, left and right
contact currents agree to ~1e-10 down to µA/m² leakage; without it they disagreed by 3%.

## Grid and geometry

Nodes sit at positions x_k, and segments join consecutive nodes. Each node owns a box reaching
halfway to each neighbour (a vertex-centred finite volume, the "box method"):
`vol_k = (h_{k−1} + h_k)/2`.

Every region boundary is a **doubled node**: two nodes at the same x, one per region, each
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

### What lives where

| Place | Holds |
|---|---|
| **Node** (box centre) | the unknowns φ̂, η_i, and hence c_i. The node's region supplies μ°, c_ref, presence and fixed charge. Everything volumetric, times the box volume: storage `vol·(c − c_old)/dt`, space charge `vol·(FΣz c + ρ_fixed)`, bulk reaction rates `vol·ν·r`. One Poisson row and one balance row per species. |
| **Segment** (always within one material) | particle fluxes N_i (Scharfetter–Gummel, with the segment's D) and displacement D = −εΔφ/h (the segment's ε). Each is computed **once** and added with + to the left box and − to the right box, which is what makes conservation telescope exactly. |
| **Interface** (zero width) | a zero-volume **flux node** between the doubled nodes (below). Its unknowns are what crosses (D_f and each N_f,i); its rows are the surface's laws. |
| **Contact** | conditions on the outer face of the end node's half-box (below). |

Output arrays are indexed by node, so each interface x appears twice and steps plot as true
vertical lines. Nothing is averaged across a material step: no harmonic means of D or ε across
an interface, because no segment crosses one.

## Fluxes

Within a segment (constant D, μ°, ε), the Scharfetter–Gummel flux is the exact constant-flux
solution for a linear potential:

```
N = (D/h) [ B(Δ) c_L − B(−Δ) c_R ],   Δ = z (φ̂_R − φ̂_L),   B(x) = x / (eˣ − 1)
```

It's evaluated in an algebraically identical, better-conditioned form. Using `B(−Δ) = B(Δ)e^Δ`
and `c_R e^Δ = c_L e^{Δη}`:

```
N = −(D/h) · B(Δ) · c_L · expm1(η_R − η_L)
```

This is exactly zero at equilibrium (flat μ̄), and it's precise relative to the quasi-Fermi
difference instead of being a cancellation of two large terms. B and B′ are evaluated from a
Taylor series near 0 and with `expm1` elsewhere, finite for all arguments.

## Interfaces: doubled nodes plus a flux node

In the linear system, a zero-volume flux node sits between the two sides' nodes:

```
 … [node L] — [flux node f] — [node R] …
    φ̂_L, η_L     D_f, N_f,i      φ̂_R, η_R
```

Node L's balance rows treat D_f and N_f,i as outflow, and node R's treat them as inflow. The
flux node's own rows are the interface laws:

| Law | Row |
|---|---|
| φ, `dipole` | `φ̂_R − φ̂_L = dipole/V_T` |
| φ, `neutral` | `D_f = 0` |
| φ, `capacitive` | `D_f + C·V_T·(φ̂_R − φ̂_L − dipole/V_T) = 0` |
| species, equilibrium | `η_R − η_L = 0` |
| species, blocked or absent on one side | `N_f = 0` |
| species, conductance | `N_f − G V_T (η_L − η_R)/(z² F) = 0` |
| species, transfer kinetics | `N_f − Σ ν r(state_L, state_R) = 0` |

Every row couples only neighbours, with a fixed block size and no penalty terms. A continuity
row has a zero diagonal block, but block Thomas still sees a non-singular block there once the
preceding elimination has run. Interface fluxes come out as unknowns, which is what the
solution reports per interface.

(An earlier design kept only the doubled nodes and put the flux in the right node's slot, with
the value read from the left node. That breaks tridiagonality, because the right node's other
segment then reaches two nodes back.)

A **sheet charge** at a face is booked in node R's Poisson row. With a pinned dipole that's
immaterial: the jump is fixed either way, and only the reported D_f shifts by σ. With a
capacitive law it would matter, since the charge then sits on one plate or the other.

## Contacts

The end node's balance rows are complete except for the flux through the outer face. So
**before** a row is replaced by a contact condition, its residual *is* that flux (or
displacement). Contact fluxes and terminal currents are recovered exactly this way, with no
extra unknowns. Then:

- **Species in equilibrium with the outside phase:** the balance row becomes the Dirichlet row
  `η = target` (the known outside level).
- **φ `bulk`:** the Poisson row becomes local neutrality (total charge, mobile plus fixed, is
  zero). The recorded residual is the outside's surface charge.
- **φ `dipole`:** the Poisson row becomes the Dirichlet row φ_edge = V − zeroCharge, again with
  the residual recorded as the outside's charge.
- **φ `capacitive`:** a term `C·((V − zeroCharge) − φ_edge)` joins the Poisson row.
- **φ `neutral`:** nothing is added (D = 0 at the face).
- **Conductance links and electrode reactions** add their fluxes to the balance rows.
  Reactions use the node's concentrations (behind any Stern layer, hence Frumkin effects)
  and the metal's electrons at −F·V.

A left-end "contact flux node" would have a singular first block for block Thomas, which is why
contact fluxes are read from residuals instead.

### Floating terminals

In galvanostatic and load modes the right terminal's voltage V_t is unknown.

- **If a charged terminal species is in equilibrium there,** V_t is read off that species' own η at
  the contact node. Its row becomes the circuit law `I_segment − I_circuit(V_t) = 0`.
  I_segment is the total current (conduction plus displacement) through the last grid segment,
  which equals the terminal current exactly. The last box's species balance plus its Poisson
  row, differenced in time, telescope to that identity. So the circuit is a local condition at
  the rightmost node.
- **Otherwise** (a kinetic or conductance electrode), V_t becomes the unknown of one extra
  block after the last node. Its row is the circuit law with the contact current: reaction
  electrons, conductance currents, and the Stern displacement current. All of those depend
  only on V_t and the last node, so the structure stays tridiagonal.

## Bulk reactions

`r = k_f Π c_R^ν · (−expm1(−a))`, with `a = A/RT` computed from the compensated η. That's mass
action with the reverse rate implied by the standard potentials: exactly zero at A = 0, and
free of cancellation near equilibrium. Sources enter each balance as `ν·r`, with one r per
reaction per node, so every moiety (a combination that no reaction changes) telescopes exactly.

## Linear algebra

Block Thomas with partial pivoting *within* each diagonal block (zero diagonal entries are
fine) but not *between* blocks. Rows are equilibrated (each divided by its largest Jacobian
entry) before factorisation. Storage is flat `Float64Array`s, and factor and solve allocate
nothing. For 300 nodes × 7 unknowns, one factor and solve takes ~0.4 ms in Node 22.

## Newton

- The update of potential-like unknowns (φ̂, η, a floating V_t) is limited to 10 thermal
  units per iteration by uniform scaling.
- Converged when a full, undamped update is below 1e-10 (thermal units). Quadratic
  convergence makes the remaining residual negligible.
- Clear divergence (updates beyond 1e4, or ten times the first update after six iterations)
  bails out early, so the caller can take a smaller step.

## Time stepping

Backward Euler. `step(dt)` recursively halves the interval wherever Newton fails, always
landing exactly on t + dt. Minority carriers that must rise by many orders of magnitude after
a bias step are the typical cause. Contact displacement at the start of a step is carried over
from the end of the previous one, including across `set()`, so a step change in a gate voltage
shows up as displacement current.

## Steady state

- **If every species stretch is fed by a contact,** nothing is conserved on its own, and the
  true steady equations are solved directly (dt = ∞, no storage term). This matters because
  slow physics can take seconds, far beyond any L²/D estimate, and stepping would never
  finish. An example is exponentially scarce minority carriers slowly filling an inversion
  layer behind a Schottky contact.
- **Otherwise** (blocked species, or moieties conserved by reactions), backward-Euler steps at
  a huge dt (10⁶ × the slowest diffusion time) keep the storage term, which pins each
  conserved amount exactly: sum a species' rows and the fluxes cancel. dt grows ×10 (capped)
  while the state still moves. Before each huge step, each spectator's level is shifted
  uniformly to restore its amount exactly. That's exact for ideal statistics, and it guards
  against round-off creeping through the vanishing storage term. The solve finishes with one
  step at the base giant dt, where pinning is tight.
- **If Newton fails,** dt ramps up from a small value (pseudo-transient continuation), ending
  in the direct solve where applicable.

`solve()` does not advance the clock.

## Conservation

- **Discretely exact.** Face fluxes are computed once and added ± to both neighbours, the
  interface flux unknowns do the same, and reactions are `ν·r`. So summed over boxes,
  everything telescopes, and the *converged* discrete solution conserves every amount and total
  charge to round-off.
- **The iterate conserves to Newton's tolerance.** Storage is nonlinear in η, so each step's
  drift is its summed residual. Converging to 1e-10 in η makes that ~1e-15 per step in
  practice.
- **Reported.** Every solution's `conservation` lists each species stretch: its amount, its
  reference, the time-integrated flux through its contacts, and the resulting drift.

## Strictly neutral materials and unresolved double layers

A material with ε = 0 has local neutrality in place of Poisson, and no displacement along its
segments. φ there is a bookkeeping multiplier, fixed by neutrality given the μ̄'s. It's reported
as NaN if no charged species is present and nothing couples the region electrostatically. The
bulk ε → 0 limit is smooth: the Poisson row simply degenerates into neutrality. The same Newton
count holds from ε_r = 78.5 down to 1e-6, as in a linear prototype that motivated the design.

Interfaces are where it gets subtle. With a pinned dipole and double layers the grid can't
resolve (because ε is tiny, or the cells are coarse), the charge the interface needs is crammed
into the two half-boxes beside it. The bulk stays right, but the interface compositions and D_f
become mesh-dependent. The true ε → 0 limit has no interface charge and a free Donnan-type jump,
which is the `neutral` law. The defaults follow from this: `neutral` between two ε = 0
materials, `dipole` otherwise. Solutions warn when a resolved-model double layer is under-resolved.

A sub-grid Gouy–Chapman law, treating the diffuse layers analytically when λ_D ≪ h, is on the
[roadmap](../ROADMAP.md). It would handle macroscopic devices with real double-layer charge at
any grid, and tends to `neutral` as ε → 0.

At an ε = 0 | ε > 0 face (an electrode), the neutral side's boundary half-box holds the
electrode's surface charge. That slightly perturbs the carrier density there: a finite-volume
stand-in for the metal's quantum capacitance, mesh-dependent in size (0.3% on the Schottky test).
