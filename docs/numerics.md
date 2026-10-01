# Numerics

How driftlet turns the physics into equations and solves them. The guiding constraint: every
equation couples a node only to its neighbours, so the Jacobian is **block-tridiagonal** and
each Newton iteration is a single block-Thomas sweep. That constraint is what makes
interactive speeds possible, so non-local physics is out of scope.

## Unknowns

At every grid node, with block size M = 1 + (number of species):

- `φ̂ = Fφ/RT`, the dimensionless (bookkeeping) electrostatic potential;
- `η_i = μ̄_i/RT`, the dimensionless electrochemical potential of each species.

Concentrations follow from these through each material's statistics, `c = c(ζ)` with
`ζ_i = η_i − μ°_i/RT − z_i φ̂` (ideal: `c_i = c_ref,i · e^{ζ_i}`; see
[Statistics](#statistics) below). Working in μ̄ rather than c is a deliberate trade-off:

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

**Advection** at velocity v adds to the drift: in the segment's Scharfetter–Gummel flux, the
potential and the quasi-Fermi difference are both shifted by the cell Péclet number
Pe = v h/D, `Δ → Δ − Pe`, `Δη → Δη − Pe`. That's the exact constant-flux solution of
N = v c − D(c′ + c ψ′) across the cell, so pure convection–diffusion is exact on any grid,
upwinding itself automatically at large Pe.

**Eddy mixing** is discretised on each segment as `N_i = −(D_mix/h) Σ_j P̄_ij Δη_j`, where P̄
is built from the logarithmic mean of each concentration across the segment. The logarithmic
mean makes a neutral species' mixing flux exactly −D_mix Δc/h. The flux is exactly zero at
equilibrium, and exactly current-free (zᵀP̄ = 0 for any P̄ built this way). Its full Jacobian,
including the dependence of P̄ on both ends, keeps Newton quadratic.

With non-ideal statistics, ln c = ζ − ex, where the excess `ex = ζ − ln(c/c_ref)` is zero for
ideal statistics. The excess enters exactly like an extra potential, taken linear along the
segment as φ is, so the same formula holds with `Δ = z(φ̂_R − φ̂_L) + ex_R − ex_L`. That's the
excess-chemical-potential generalisation of Scharfetter–Gummel: still exactly zero at
equilibrium, since the expm1 factor is untouched. It's exact for a single diffusing species on
a lattice (whose excess is proportional to its grand potential, linear in x at steady state)
and second order otherwise. Its Jacobian couples each flux to every species at both nodes through
K, which fills the blocks but keeps them tridiagonal.

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

## Statistics

Each node evaluates its material's statistics once per Newton iteration: c(ζ), the Jacobian
K = ∂c/∂ζ (n × n, symmetric positive definite), and the excess. Every c-dependent term then
differentiates through K:

- storage and space charge: `∂c_i/∂η_j = K_ij`, `∂c_i/∂φ̂ = −(Kz)_i`;
- reaction prefactors: `∂ln c_i/∂η_j = K_ij/c_i`;
- the φ row's response, and the Debye length in resolution warnings, use zᵀKz, the charge
  capacitance, in place of Σz²c.

Nodes of ideal materials keep a separate fast path with K = diag(c) implied.

Model evaluation:

- **Explicit models** (Fermi–Dirac, lattice gas, insertion) evaluate c(ζ) in closed form.
- **Implicit models** solve a small problem per node: Redlich–Kister and tabulated OCVs by
  safeguarded Newton in logit(x); Debye–Hückel by Newton in ln c, whose excess Hessian is rank
  one, so each step is a Sherman–Morrison update.
- **Fermi–Dirac integrals** 𝓕_{1/2} and 𝓕_{−1/2} use:
  - an alternating series below x = −2;
  - Sommerfeld's expansion above 50;
  - piecewise Chebyshev fits (degree 23, five intervals) in between, built on first use
    (~4 ms) from a trapezoid quadrature in u = √t.

  Relative error is ~6e-15.

Where a *concentration* is prescribed rather than a potential (bath compositions, spectators'
initial amounts), each model inverts itself for the ζ of those species, holding the others.
The lattice gas does this in closed form, and the others with the same scalar or small Newton
iterations.

An insertion host's species depend on φ only through their neutral combination, where it
cancels, so its φ rows are identity rows (φ undefined, like a region with no charged species).
The current-continuity constraint between ion and carrier comes from their balance rows.

## Metal regions

A metal node's only unknown is the carrier's η, its Fermi level. Its φ slot is free, and it
carries the carrier flux J through the segment to the node's right instead. Each metal segment
then has two rows:

```
η_R − η_L + J/g = 0        (Ohm's law, g = σRT/(z²F²h))
… + J (out of L), − J (into R) in the carrier balances
```

This mixed form matters. Eliminating a stiff ohmic chain in η alone computes
g − g²/(g + G) wherever the chain hangs on a weak conductance G, such as a floating metal held
only by its electrode reactions. For a real metal that cancels catastrophically, with pivots
near 1e-15. In mixed form every coefficient is O(1), whatever σ.

A metal's bulk holds no charge. Its surface charge at a capacitive face is the face's
displacement, booked as a sheet of excess carriers in the edge node's half-box. That gives the
edge node a storage term in D_f (the neighbouring block), so transients, conservation and
impedance all see the charging current. A pinned (`dipole`) law isn't offered at a metal face.
With nothing on the metal side depending on D_f, the steady system's flux block would be
singular to left-to-right elimination, and the charging current would need the charge read two
blocks away.

Electrode reactions at a metal face put each participant's flux on its own side of the face.
The face's flux node carries one flux per species, and a species present on one side only flows
on that side only. A stretch exchanging through such a reaction with one that reaches a contact
counts as fed, so nothing is conserved there and the steady equations are solved directly. A
floating metal starts uncharged, with its Fermi level in equilibrium with the reaction on its
left face.

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

- **Backward Euler** is the default for `step(dt)`. It recursively halves the interval wherever
  Newton fails, always landing exactly on t + dt. Minority carriers that must rise by many
  orders of magnitude after a bias step are the typical cause.
- **Contact displacement** at the start of a step is carried over from the end of the previous
  one, including across `set()`, so a step change in a gate voltage shows up as displacement
  current.
- **BDF2** (`method: 'bdf2'`) is the variable-step second-order backward difference formula.
  With ω = dt/dt_prev, its storage term is
  `[a₀ c − (1+ω) c_n + ω²/(1+ω) c_{n−1}]/dt`, where `a₀ = (1+2ω)/(1+ω)`.
  - It's written in backward-Euler form, `(c − c*)/(dt/a₀)` with the history reference
    `c* = [(1+ω) c_n − ω²/(1+ω) c_{n−1}]/a₀`. So assembly is shared, and the displacement
    histories (for terminal currents) take the same combination.
  - Conservation still telescopes. Closed amounts stay exact, since c* has the same amount
    as c_n when c_n and c_{n−1} did. Open stretches add the history term Σ v (c* − c_n) to
    their intake.
  - A step ratio above 2 falls back to backward Euler for that step, to stay safely within
    variable-step BDF2's zero-stability limit (1 + √2).
- **Adaptive stepping** (`advance`):
  - The local error of each step is estimated against an explicit predictor through the
    previous states: quadratic after a BDF2 step, scaled by C_c/(C_c + C_p), with
    `C_c = h³(1+ω)²/(ω(1+2ω))` and `C_p = h(h+h₁)(h+h₁+h₂)`; linear after backward Euler,
    scaled by h/(2h + h₁). The first step, with no history, is checked by step doubling.
  - The error is measured in thermal units over every state potential (φ̂ where defined,
    each present η, a floating terminal voltage).
  - A step is rejected above `tol`. The next step size is h·min(2, max(0.2,
    0.9 (tol/err)^{1/(p+1)})).
  - Newton failure quarters the step.
  - Newton starts each step from the state extrapolated through the last three, which saves
    about a third of the iterations.
  - The wall-clock budget is checked between steps.

## Small-signal impedance

About a steady state x₀, a small sinusoidal source δs·e^{iωt} gives, to first order,

```
(J + iωM) δx = −b δs
```

- **J** is the steady Jacobian.
- **M** is the Jacobian of the time-derivative terms: storage v·K, displacement in the
  circuit rows. It's read off two assemblies, J(dt) = J + M/dt, with a tiny dt so that the
  subtraction loses nothing.
- **b = ∂(residual)/∂s** comes from central differences in the source: the right terminal's
  voltage in voltage mode, the circuit current in current mode. It's exact wherever the
  residual is linear in the source.

The system keeps the block-tridiagonal structure, so each frequency costs one complex block-Thomas
factorisation (rows scaled by their largest entry). The terminal current comes from the last
grid segment, conduction plus iω times its displacement, by the same identity as the floating
terminal. In current mode, the voltage response is read from the terminal unknown instead.
Conserved (blocked) species make J singular, but J + iωM isn't for ω > 0: at low frequency a
blocking device looks like a capacitor, as it should.

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
  uniformly to restore its amount exactly (in one step for ideal statistics, by Newton on the
  shift otherwise). It guards
  against round-off creeping through the vanishing storage term. The solve finishes with one
  step at the base giant dt, where pinning is tight.
- **If a direct solve diverges,** it's retried once with tighter damping (3 thermal units per
  iteration). That's enough for most large jumps, such as a cold start at forward bias.
- **If that fails too,** source continuation ramps the right terminal's voltage to its target.
  It ramps from the voltage of the last converged solve when the state is that solution, and
  otherwise from level terminals, where a cold start is consistent. The ramp step starts at
  1/8 of the way, grows ×1.5 on success and shrinks ×4 on failure.
- **If Newton still fails,** dt ramps up from a small value (pseudo-transient continuation),
  ending in the direct solve where applicable.

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

At an ε = 0 | ε > 0 face, the neutral side's boundary half-box holds the face's charge. That
slightly perturbs the carrier density there, by an amount that depends on the mesh (0.3% on a
Schottky test with the metal as an ε = 0 region). For a metal, a metal region avoids this: its
surface charge is a sheet that leaves its Fermi level untouched.
