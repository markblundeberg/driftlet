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
| φ, `pinned` | `φ̂_R − φ̂_L = dipole/V_T` |
| φ, `neutral` | `D_f = 0` |
| φ, `capacitive` | `D_f + C·V_T·(φ̂_R − φ̂_L − dipole/V_T) = 0` |
| species, equilibrium | `η_R − η_L = 0` |
| species, blocked or absent on one side | `N_f = 0` |
| species, conductance | `N_f − G V_T (η_L − η_R)/(z² F) = 0` |
| face reaction k | `r_k − rate_k(state_L, state_R) = 0`, and each participant's edge node takes `ν·r_k` |

Every row couples only neighbours, with no penalty terms. A face reaction's rate is an unknown
of the face block, so the two edge nodes it couples meet only through the block between them.
A species can take part in a reaction at a face and also cross it by a link (its flux slot)
at the same time; where it exists on both sides, that link must be given explicitly. A continuity
row has a zero diagonal block, but block Thomas still sees a non-singular block there once the
preceding elimination has run. Interface fluxes come out as unknowns, which is what the
solution reports per interface.

(An earlier design kept only the doubled nodes and put the flux in the right node's slot, with
the value read from the left node. That breaks tridiagonality, because the right node's other
segment then reaches two nodes back.)

A **sheet charge** at a face is booked in node R's Poisson row. With a `pinned` law that's
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
- **φ `pinned`:** the Poisson row becomes the Dirichlet row φ_edge = V − zeroCharge, again with
  the residual recorded as the outside's charge.
- **φ `capacitive`:** a term `C·((V − zeroCharge) − φ_edge)` joins the Poisson row.
- **φ `neutral`:** nothing is added (D = 0 at the face).
- **Conductance links** add their fluxes to the balance rows. (Contacts carry no reactions:
  an electrode is a conductor region with reactions at its face.)

A left-end "contact flux node" would have a singular first block for block Thomas, which is why
contact fluxes are read from residuals instead. After each step these readouts (contact and port
fluxes, contact displacements, the last segment's current) are evaluated at the converged state
from only the boxes they come from, the two end nodes and the port windows with the segments and
faces that touch them, in the same order as a full assembly, so they're identical to one.

### Terminals

Every terminal (a contact, or a port) has a voltage V, held by its source or floating, and a
current into the device, I. Assembly records, for each:
- **B = ∂res/∂V**, a column: where its outside levels enter (Dirichlet rows, conductance links,
  a gate's or pinned φ law);
- **I and C = ∂I/∂x**, a row. A contact's current is read from its end box before the contact's
  own terms go in: the box's balance residuals are what must come through the contact (as for
  the flux readouts), plus the displacement through its φ law, (D_in − D_in,0)/dt. A port's is
  its sources over its window, with a held level's read from the row it replaces.

A floating terminal's voltage is an extra unknown, with its circuit law as an extra row:
I − I_set = 0 (driven by a current) or I − (V_src − V)/R = 0 (behind a resistance). These don't
fit the block-tridiagonal matrix T (a port couples to its whole window), so they're solved by
bordering, together with the spectators' conservation rows (see [steady state](#steady-state)):
with y = T⁻¹ rhs, X_k = T⁻¹ B_k and Q_q = T⁻¹ e_q (the response to a unit pin), the update is
δ = y + Σ Q_q μ_q − Σ X_k δV_k, and the μ (pins) and δV (terminals) come from a small dense
system of the extra rows. Each costs one more back-substitution per Newton iteration.

Sources are read at the end of a step (implicit), or at the present time in a steady solve.
`advance()` lands on every waveform breakpoint and restarts its order there.

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

This mixed form matters. A floating metal, held only by weak conductances G such as its
electrode reactions, has its overall level set by G alone, while every entry of an η-only
Jacobian is of order g. Eliminating it then computes the last pivot as g − g²/(g + G) ≈ G with
an absolute error of order ε·g: digits are lost in proportion to g/G, and beyond g/G ≈ 1e15 the
pivot is noise or zero. (G is already lost when the diagonal g + G is stored.) In mixed form
the matrix holds 1/g instead, every coefficient is O(1), and the limit σ → ∞ is simply the
constraint η_R = η_L. A Grassmann–Taksar–Heyman-style elimination, which carries each row's
leakage separately so that pivots are built as sums, would also keep η alone exact, but it needs
to know which rows form a conductance network.

A metal region is a single cell. Nothing is stored in its bulk, so η is linear across it in
transients too, and one segment is exact. That also keeps g = σRT/(z²F²L) as small as the
region allows.

A metal's bulk holds no charge. Its surface charge at a capacitive face is the face's
displacement, booked as a sheet of excess carriers in the edge node's half-box. That gives the
edge node a storage term in D_f (the neighbouring block), so transients, conservation and
impedance all see the charging current. A `pinned` law isn't offered at a metal face.
With nothing on the metal side depending on D_f, the steady system's flux block would be
singular to left-to-right elimination, and the charging current would need the charge read two
blocks away.

Electrode reactions at a conductor's face take its carrier at its edge node, like any face
reaction. Whether a stretch is fed is decided by its conserved combinations: weightings w of the
stretches' amounts that no reaction changes (w·ν = 0 for every face and bulk reaction) and that
nothing outside feeds (w = 0 on stretches reached by a contact or a port), the null space of that
stoichiometry. A stretch in none of them is fed. So Ag⁺ between silver electrodes is fed (Ag⁺ + e⁻
⇌ Ag(s), the electrons fed by the contacts), and the steady equations are solved directly; Fe³⁺
and Fe²⁺ between platinum electrodes aren't, since Fe³⁺ + e⁻ ⇌ Fe²⁺ conserves the iron whatever
the electrons do, so their total is kept by huge steps (below). A floating
conductor starts uncharged, with its carrier's level in equilibrium with the first reaction on its
left face that takes it.

## Internal ports

A port adds a source per volume to the balance rows of the nodes in its window, after every
other term at those nodes and before the contacts. A held (`'equilibrium'`) level replaces the
balance row with a Dirichlet row. The port's flux is then that row's residual, read just
before replacement, exactly as at a contact. Conductance and exchange links are linear in
(μ̄_out − μ̄). Everything stays on the node's own block. A contact's flux readout at an end node
already includes any port source there, so contact and port fluxes always balance. A stretch
reached by a port counts as fed, and its conservation intake includes the port's flux. On a
metal, a port spans both nodes of its single cell, with G per area spread as G/L per volume, so
each node gets half.

## Bulk reactions

`r = k_f Π c_R^ν · (−expm1(−a))`, with `a = A/RT` computed from the compensated η. That's mass
action with the reverse rate implied by the standard potentials: exactly zero at A = 0, and
free of cancellation near equilibrium. Sources enter each balance as `ν·r`, with one r per
reaction per node, so every moiety (a combination that no reaction changes) telescopes exactly.

## Assembly

Each region is assembled by the kernel for its kind, over all its nodes and segments at once:
- **conductor:** its carrier's balance, and Ohm's law in mixed form (below);
- **dilute** (ideal statistics, the fast path): storage, space charge and Scharfetter–Gummel
  fluxes, written straight into the blocks by local index; a dielectric is the case with no
  species;
- **concentrated** (any statistics): the same through K = ∂c/∂ζ and the excess potential.

Faces, ports and contacts follow. Every kernel writes only the unknowns its region has.

## Linear algebra

Block Thomas with partial pivoting *within* each diagonal block (zero diagonal entries are
fine) but not *between* blocks. Rows are equilibrated (each divided by its largest Jacobian
entry) before factorisation. Storage is flat `Float64Array`s, and factor and solve allocate
nothing. For 300 nodes × 7 unknowns, one factor and solve takes ~0.4 ms in Node 22.

Block sizes differ from node to node: each block holds only the unknowns that exist there. A
node has M = 1 + n slots (φ̂ and each species' η), but these aren't unknowns:
- a species absent from its material;
- φ where it's undefined;
- a conductor node's other species, and the flux slot of its last node;
- a blocked interface flux, or the displacement of a `'neutral'` face;
- a floating terminal's spare slots.

The state vector keeps every slot (one that isn't an unknown keeps its value), while the
residual, the update and the Jacobian are compact: assembly writes each entry straight into
blocks of m_b × m_b, m_b × m_{b−1} and m_b × m_{b+1}, through each slot's row within its block,
and drops terms in slots that aren't unknowns. Elimination costs Σ m_b³ instead of n·M³, which
pays wherever species are confined to some regions: for n-Si against KCl, three unknowns of five
at every node, it's a fifth of the work. A finite-difference test checks the compact Jacobian
column by column on devices that cover every assembly path.

## Newton

- The update of potential-like unknowns (φ̂, η, a floating V_t) is limited to 10 thermal
  units per iteration by uniform scaling.
- Converged when a full, undamped update is below 1e-10 (thermal units). Quadratic
  convergence makes the remaining residual negligible.
- Or converged as far as round-off allows: updates already below 1e-6 thermal units (~26 nV)
  in a steady solve, or 1e-5 on a time step (whose error control works at ~1e-3), that have
  stopped shrinking for two iterations. A badly conditioned system's round-off floor can sit
  above 1e-10: a strictly neutral material on a short step, where φ is fixed only through
  fluxes that the storage term dwarfs (a condition number of about h²/(D·dt), much reduced by
  the change of variables below), converges quadratically to ~1e-9 and then rattles there.
  Short steps come right after every waveform breakpoint, so without this a cyclic voltammogram
  in a neutral electrolyte stalled at its turns.
- Clear divergence (updates beyond 1e4, or ten times the first update after six iterations)
  bails out early, so the caller can take a smaller step.
- When a steady solve fails, the solution's warnings say how nearly singular the system was,
  and where. A running error bound through the factorisation compares each pivot with the
  magnitudes it was formed from, log₁₀(Σ|terms| / |pivot|): the digits lost. A part of the
  device held only weakly, such as a floating region coupled through tiny conductances, loses
  about log₁₀(g/G) and fails near 15. It's not reported for solves that converge, since digits
  can be lost harmlessly too (a saturated species whose η nothing depends on).
- The last iteration usually only confirms convergence (updates go like 4e-3, 7e-6, 3e-11).
  Stopping one iteration early on a convergence-rate estimate would save a third of the work,
  but leave ~1e-11 in η per step, which shows up as conservation drift, so it isn't done.
- Tried and rejected for large jumps (the retry at 3 thermal units per iteration takes 30–50
  iterations for a 1.4 V jump on a pn diode): per-component logarithmic damping, as in
  SPICE junction limiting, and an adaptive limit that grows while updates shrink. Both
  diverged more often and cost more overall.

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
    each present η, a floating terminal voltage). That includes species far below their
    largest concentration, such as minority carriers ahead of a diffusion front. Weighting
    them down (a mixed relative/absolute tolerance, δη·c/(c + 1e-10 c_max)) took a third fewer
    steps on a pn turn-on, but the tail's errors are swept into the front, and at matched
    accuracy of the terminal current it cost 50% more.
  - A step is rejected above `tol`. The next step size is h·min(2, max(0.2,
    0.9 (tol/err)^{1/(p+1)})).
  - Newton failure quarters the step; on the first step after a start or a jump in a device
    with a strictly neutral material, longer steps are tried first, since there shorter ones are
    worse conditioned.
  - A target within round-off of the present time (1e-10 relative), such as an animation
    frame's that lands a hair past a breakpoint just reached, is snapped to rather than
    stepped to: a step of 1e-13 s can't be resolved.
  - Newton starts each step from the state extrapolated through the last three, which saves
    about a third of the iterations.
  - The wall-clock budget is checked between steps.

## Small-signal impedance

About a steady state, a small sinusoidal excitation at one terminal gives, to first order,

```
(J + iωM) δx + Σ_k B_k δV_k = −B_T δV_T                 (the grid's rows)
(C_k + iωC′_k) δx + (∂I_k/∂V_k + iω…) δV_k = δI_T δ_kT  (each floating terminal's circuit row)
```

- **J** is the steady Jacobian, **M** the Jacobian of the time-derivative terms (storage v·K,
  displacement), and likewise C and C′ for the terminal currents. The time-derivative parts are
  read off two assemblies, J(dt) = J + M/dt, with a tiny dt so that the subtraction loses
  nothing.
- **At the measured terminal T,** a held voltage is perturbed (δV_T = 1, and δI_T read from its
  C row), or a driven current (δI_T = 1, and δV_T solved for). Z = δV_T/δI_T, with I into the
  device. The other terminals keep their drives: held ones at AC ground, driven ones open.

Each frequency costs one complex block-Thomas factorisation (rows scaled by their largest
entry), plus one back-substitution per floating terminal for the bordering. Conserved (blocked)
species make J singular, but J + iωM isn't for ω > 0: at low frequency a blocking device looks
like a capacitor, as it should.

## Steady state

- **If every species stretch is fed by a contact,** nothing is conserved on its own, and the
  true steady equations are solved directly (dt = ∞, no storage term). This matters because
  slow physics can take seconds, far beyond any L²/D estimate, and stepping would never
  finish. An example is exponentially scarce minority carriers slowly filling an inversion
  layer behind a Schottky contact.
- **Spectators** (a species blocked all round, mobile throughout its stretch) are solved
  directly too. In steady state the sum of a spectator's balance rows over its stretch is zero
  identically, so one of them (the first node's) is redundant, and it's replaced by the
  conservation of the amount, Σ v c = amount. That row is dense, so the solve is bordered:
  the matrix is factorised with a pin (an identity row) in its place, which makes the stretch's
  level a well-conditioned unknown, and the response to a unit pin (one extra
  back-substitution per spectator) is added in the amount that satisfies the constraint, from
  a k×k system for k spectators.
- **Conserved combinations** (moieties) of reacting stretches are solved the same way: the
  total iron of Fe³⁺, Fe²⁺ and FeCl²⁺ in a closed cell with a complexation reaction, say. The
  basis of combinations w comes from the null space of the stoichiometry. In steady state the
  w-weighted sum of their balance rows is zero identically (the reactions cancel by w·ν = 0,
  and nothing crosses the stretches' ends), so one of them is replaced by Σ_k w_k Σ v c_k =
  the combination's amount, bordered like a spectator's. (Spectators are the one-stretch
  case.)

  Holding the amount through a storage term at huge dt instead is badly conditioned: the
  stretch's level is then held only by v·c/dt, against internal conductances D·c/h larger by
  D·dt/(hL), around 1e12 for a micron-scale cell (see [metal regions](#metal-regions) for the
  same cancellation). Newton then needs dozens of iterations, or fails.
- **Otherwise** (immobile spectators, which conserve node by node, or a combination that
  includes a floating conductor), backward-Euler steps at a huge dt (10⁶ × the slowest diffusion time) keep the storage
  term, which pins each conserved amount exactly: sum a species' rows and the fluxes cancel. dt
  grows ×10 (capped) while the state still moves. Before each huge step, each spectator's level
  is shifted uniformly to restore its amount exactly (in one step for ideal statistics, by
  Newton on the shift otherwise). It guards against round-off creeping through the vanishing
  storage term. The solve finishes with one step at the base giant dt, where pinning is
  tight.
- **If a direct solve diverges,** it's retried once with tighter damping (3 thermal units per
  iteration). That's enough for most large jumps, such as a cold start at forward bias.
- **Generation continuation.** A device with generation reactions (species made only from, or
  turned only into, fixed reservoirs, such as photogeneration from a photon reservoir) can be
  held far from equilibrium even with its terminals level, where bias continuation can't help.
  If a direct solve fails there, the generation rates are scaled down to 10⁻¹² and ramped back
  up, ×100 a step while each solve converges (warm from the last) and by the square root of the
  factor when one doesn't. An illuminated 80 µm silicon diode solves cold this way in about 80
  iterations.
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

Interfaces are where it gets subtle. With a `pinned` law and double layers the grid can't
resolve (because ε is tiny, or the cells are coarse), the charge the interface needs is crammed
into the two half-boxes beside it. The bulk stays right, but the interface compositions and D_f
become mesh-dependent. The true ε → 0 limit has no interface charge and a free Donnan-type jump,
which is the `neutral` law. The defaults follow from this: `neutral` between two ε = 0
materials, `pinned` otherwise (and `pinned` between two ε = 0 materials is an error: with no field
on either side, nothing would determine its charge). Solutions warn when a resolved-model double
layer is under-resolved.

On a transient step, strictly neutral nodes are solved in better-conditioned terms. There, a
change of φ̂ with every η_i shifted by z_i times it leaves every concentration as it was, so
storage and neutrality don't see it. Only the fluxes do. But in (φ̂, η) storage and neutrality
see it as pairs of huge entries that cancel, and round-off in that cancellation, of order
storage/flux ~ h²/(D·dt) times the range of concentrations, swamps the fluxes that fix it: a
trace ion beside 3 M KCl lost about ten digits, and Newton stalled. So at each node that stays
neutral (interior nodes, and edges at `neutral` faces; not the edge of a capacitive or pinned
face, which holds the face's charge), assembled without storage and neutrality:
- rows: the balance of the most abundant charged species (by z²c) becomes Σ (z_i/z_k) × each
  balance, which without storage is current continuity;
- columns: the unknowns become φ̂' and η'_i = η_i − z_i φ̂ (the update is mapped back), in which
  storage and neutrality have no φ̂' term at all;
- then storage (on the other balances) and neutrality go in, exactly, in those terms.

A start-of-step net charge (round-off in a solved state) isn't carried over: the neutrality row
holds the new state neutral. Steady solves and the impedance keep the plain rows (there's no
storage at dt = ∞, and the impedance reads the storage matrix from them).

A sub-grid Gouy–Chapman law, treating the diffuse layers analytically when λ_D ≪ h, is on the
[roadmap](../ROADMAP.md). It would handle macroscopic devices with real double-layer charge at
any grid, and tends to `neutral` as ε → 0.

At an ε = 0 | ε > 0 face, the neutral side's boundary half-box holds the face's charge. That
slightly perturbs the carrier density there, by an amount that depends on the mesh (0.3% on a
Schottky test with the metal as an ε = 0 region). For a metal, a metal region avoids this: its
surface charge is a sheet that leaves its Fermi level untouched.
