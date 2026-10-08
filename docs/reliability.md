# How far to trust it

A drift–diffusion solver is quick to start and slow to finish. The first version solves a pn
junction; the trouble comes after, with what real devices ask of it:

- concentrations spanning 30 orders of magnitude in one device (majority carriers beside
  minority ones, a salt beside the ions its double layer excludes);
- double layers a millionth of the device's length, resolved or not;
- species whose amount is fixed rather than their level (a blocked ion, a floating base, a
  closed cell's salt), which makes the steady equations singular;
- populations held only by couplings of 1e-17 (an inversion layer fed by ~1e3 minority
  carriers per cm³), which a factorisation loses to round-off;
- terminals driven by a current, which may be beyond what the device can pass;
- small-signal currents read from sums whose terms are 1e17 times bigger.

Each of these makes a plain Newton solver fail, or worse, converge to the wrong answer. driftlet
has met them, mostly through the two kinds of testing below, and its numerics are built around
them ([numerics](numerics.md)). This page says what that testing covers and what it found, and
what is known not to work.

## Correct: the validation suite

Every physics feature is checked against an analytic result, or an independent code where
there isn't one: `npm test`, over 300 tests, each a worked device with its closed form. They include
Shockley's diode equation, depletion and MOS charge, collection theory under Beer–Lambert light,
Gouy–Chapman and Kilic–Bazant–Ajdari double layers, Donnan and Teorell–Meyer–Sievers
partitions, Goldman–Hodgkin–Katz, Hodgkin and Huxley's action potential and its conduction
velocity, the limiting current and Butler–Volmer kinetics, the
Macdonald and finite-length Warburg impedances, Cottrell and spherical diffusion, and the
second law (the free energy in through the terminals equals what's dissipated, every term
non-negative). Against other codes: IonMonger's perovskite J–V hysteresis (finite elements,
Octave), and liquid-junction potentials against JPCalc, LJPcalc and JLJP. The
[table below](#the-validation-table) lists each with its tolerance.

Each solution can check itself too: `check()` from `driftlet/kit` balances every species'
ledger (what comes in through each terminal, is made or destroyed), the conservation
bookkeeping of a transient, and the grid (the device solved again on a grid twice as fine).

### The validation table

Each check, with the tolerance it's held to. The tests are worked setups, so they double as
recipes.

| What | Checked against | Tolerance |
|---|---|---|
| Equilibrium of any device | μ̄ of every species flat; zero current | 1e-9 RT or better |
| Gouy–Chapman double layer | analytic charge and full nonlinear profile | 2e-3 |
| Debye screening (linear limit) | gate in series with `ε/λ_D` | 1e-3 |
| Donnan potential | analytic partition, between floating layers and real baths | 1e-8 |
| Heterointerfaces with unequal dipoles (A \| B \| A) | per-face analytic double-layer split | 3e-3 |
| Floating island | Gauss's law; conserved amounts through gate sweeps | 1e-12 |
| pn junction | exact built-in potential; depletion charge; short-diode J–V | 1e-12; 1%; 2e-3 |
| Long pn diode with recombination | Shockley J–V incl. depletion recombination | 2e-3 |
| Haynes–Shockley pulse | an injected hole packet's drift, spread and decay at the ambipolar μ*, D*, 1/τ*; nothing lost while it goes in | 2e-4; 1e-2; 2e-3 |
| A packet laid down as a `c0` profile | the same rates from a Gaussian at t = 0: μ*E t, σ₀² + 2D*t, e^(−t/τ*) | 3e-4; 1e-2; 1e-3 |
| Beer–Lambert photogeneration in an n⁺p cell | J_sc = qΦ∫αe^(−αx)η(x)dx with the emitter's, depletion layer's and base's collection, α from 1e3 to 1e7 /m on one coarse grid | 5e-4 (3e-3 in the emitter) |
| Organic bilayer: excitons | J_sc from excitons diffusing to the donor/acceptor interface (Beer–Lambert, a blocking anode, the splitting reaction's finite velocity); V_oc from charge-transfer detailed balance | 5e-4; 1e-6 V |
| SRH recombination | bulk: G = R_SRH(n, p) with n = N_D + p, midgap and shallow traps, low to high injection; at a face: zero in the dark, the law between the two edge nodes, every pair collected or recombined | 1e-6; 1e-9 |
| Perovskite J–V hysteresis, against IonMonger | an independent code (finite elements, Octave): mobile iodide vacancies, interface SRH, scans from 1 mV/s to 1 kV/s — hysteresis index, P_max, V_oc, the whole loop (as current where it's gentle, as a voltage offset near V_oc, where it falls at up to 1100 mA/cm² per volt) | 1.5e-3; 0.05 mW/cm²; 2 mV; 0.03 mA/cm² and 1 mV (0.010 and 0.31 mV found) |
| Illuminated long pn diode | `J_sc = qG(L_n + L_p + W)` from a cold start; superposition at low injection | 5e-3; 1e-2 |
| Schottky barrier (metal region \| n-Si) | surface density from the alignment; depletion charge | 5e-3; 2% |
| Surface velocity at a contact | minority electrons through a base to a contact of velocity S, n₀(e^(V/V_T) − 1)/(W/D + 1/S), S from 10 to 1e5 m/s; a Schottky diode's thermionic emission, short of A*T²e^(−φ_B/V_T)(e^(V/V_T) − 1) by its depletion region in series, the same swept as built cold | 1e-5; 3–7%; 1e-6 |
| Charged pore or nanochannel (resolved Donnan layers at each mouth) | conductance against salt, down to its plateau, as Teorell–Meyer–Sievers in a long pore (short ones short by end layers ∝ 1/L); the zero-current voltage across a salt gradient; a symmetric pore's I–V odd, rectification with charge on one side or bipolar | 1%; 5e-3 to 1e-4; 1e-9 |
| Cell membrane (a capacitive face with permeabilities) | GHK resting potential (neutral or resolved solutions) and current–voltage curve; the same as a resolved 5 nm lipid layer; charging at I/C; with the Na⁺/K⁺ pump in a closed cell, Mullins–Noda, the static head 3Δμ̄_Na − 2Δμ̄_K = ΔG_ATP, and Donnan with the pump off | 1e-5 V; 3e-3; 2e-5 V; 1%; 1e-9 |
| Voltage-gated channels (Hodgkin–Huxley gates on a membrane face) | each gate at α/(α + β) under a voltage clamp, and relaxing at α + β after a step; HH's voltage clamp with their linear channels between baths, the K⁺ current g_K n∞⁴ (V − E_K); a squid axon's action potential against the space-clamped HH equations with GHK currents; its impedance at rest (inductive, resonant) against the linearised equations | 1e-12; 2e-4; 2e-3; 1 mV; 5e-4 |
| A propagating action potential (a squid axon along x: gated membrane port, ions drifting in the axoplasm) | rest at HH's −65 mV and gate values; conduction velocity against the cable equation solved directly, and HH's 18.8 m/s; a myelinated axon's saltatory conduction, node by node, against a compartmental cable | 0.01 mV, 1e-4; 2e-3; 1%; 1 µs of ~25 per node |
| Liquid junctions in mixtures | Henderson's formula against JPCalc's published values; driftlet's steady junction against Planck by shooting and the stationary Nernst–Planck codes LJPcalc and JLJP, a patch pipette's K-gluconate among them; a free-diffusion junction constant in time, between the two for the mixtures tested | 0.015 mV; 2e-4 and 0.02 mV; 5e-3 mV |
| Liquid junction, open circuit | cell EMF 2t₊(RT/F) ln(c₁/c₂); Planck diffusion potential | 1e-4 |
| Salt diffusion mode from a `c0` profile | c̄ + a·sin(πx/2L) decaying at π²D/4L², D = 2D₊D₋/(D₊ + D₋) | 2e-5 |
| Concentration polarization | `i = i_lim tanh(V/4V_T)`, incl. galvanostatic and load modes | 2e-4 |
| Butler–Volmer electrode | Nernst equilibrium; mixed kinetic/diffusion closed form | 5e-4 |
| Redox couple between inert electrodes | Nernst level of Fe³⁺/Fe²⁺; the couple's total conserved at equilibrium and under current | 1e-9; 1e-10 |
| Interface conductance and ion transfer | series 1/G; BV rate law at the interface state | 1e-4; 1e-8 |
| Mass action | `c(H⁺)c(OH⁻) = K_w` from standard potentials; moiety conservation | 1e-9 |
| Second law | free energy in through the terminals (Σ N μ̄ over every species) equals the dissipation, every term ≥ 0: a pn diode, an open-circuit junction running on chemical input alone, electrodes with a bipolar plate, a port and a contact behind conductances | 1e-10 |
| Transients | RC charging of a gated island; water relaxation rate | 1%; 2e-3 |
| Time integration | BE first order, BDF2 second order; adaptive error control | ratios 2, 4 |
| Advection and eddy mixing | exact convection–diffusion profile; `D + D_mix`; junction EMF unchanged by mixing | 1e-12; 1e-12; 1e-5 |
| Metal regions | ohmic; Schottky face and MOS gate equal their contact forms; Ag \| AgNO₃ \| Ag with metal electrodes; bipolar electrode at V/2 | 1e-12; 1e-9; 1e-9; 1e-9 |
| Internal ports | transmission line σV tanh(L/λ)/λ (O(h²)); held level; MOS low-frequency C–V with a grounded channel, and without one the same equilibrium, and a gate step filling the inversion layer to it | 1e-4; 1e-9; 2e-3; 1e-6; 1e-3 |
| Cross-sections (spherical, cylindrical, any A(x)) | steady diffusion to a sphere and a cylinder (exact on any grid); Cottrell with the spherical term; uptake by a sphere filling from its surface; Debye–Hückel around a charged sphere, potential and charge | 1e-10; 3e-4; 2e-4; 1e-3, 3e-4 |
| A capacitance through a port | a thin-film transistor against the charge-sheet model (below threshold to saturation, second order); de Levie impedance; gate impedance against dQ/dV | 5e-4; 1e-3; 1e-6 |
| Electrode surfaces (coverage) | Langmuir isotherm against potential; the charge to fill a surface; the active–passive curve, blocking exactly (1 − θ) | 1e-9; 1e-3; 1e-9 |
| Electrodes spread through a port | transmission line with the reaction's linear kinetics as the conductance (O(h²)), rates along it; Wagner–Traud mixed potential of two Butler–Volmer couples, floating | 1e-4, 2e-3; 1e-6 V |
| Impedance | Macdonald blocking-electrode spectrum, 100 Hz–1 GHz; finite-length Warburg (Ag \| AgNO₃ \| Ag); DC limit = differential resistance; a p⁺n⁺ junction without recombination, Re Y → dI/dV and C = dQ/dV; a redox electrode against a bath at its open circuit, Re Y → dI/dV; a MOS capacitor with a metal gate and a channel port, C = dQ/dV at the gate, passive at the back | 3e-4; 3e-5; 1e-6; 1e-4, 1e-6; 1e-4; 1e-3 |
| Conservation | per step, and against time-integrated contact fluxes | 1e-11 relative |
| Strictly neutral limit (ε = 0) | Planck EMF; polarization with no overlimiting; Donnan at neutral faces; an ion exchanger between baths (Donnan at the contacts, Teorell–Meyer–Sievers salt flux); a capacitive face's edge cell as a diffuse layer in series | 1e-5; 1e-4; 1e-9; 1e-9, 1e-6; 2e-3 |
| Fermi–Dirac statistics | `𝓕_{±1/2}` vs quadrature; degenerate bulk; accumulation charge via `𝓕_{3/2}` | 1e-13; 1e-12; 1e-3 |
| Crowded double layer (lattice gas) | Kilic–Bazant–Ajdari charge, up to `ψ = 40 V_T`; custom function reproduces it | 5e-4; 1e-10 |
| Non-ideal transport | steady flux −(D/L)ΔP of the grand potential (lattice exact, Redlich–Kister O(h²)) | 1e-12; 1e-4 |
| Debye–Hückel | junction EMF 2t₊(RT/F) ln(a₁/a₂) with activities | 1e-6 |
| Data library | ion μ° against the electrochemical series and `K_w`; D against limiting conductivities; band data against `n_i` | 10 mV; 0.5%; 15% |
| Intercalation host (OCV) | composition vs table and isotherm; chemical diffusion flux and relaxation rate | 1e-12; 1e-5; 1e-4 |

## Robust: random devices

The validation suite shows driftlet is right where it converges. Whether it converges, on
devices nobody tuned it for, is what `npm run stress` tests. It builds random but plausible
devices in five families:

| Family | What's drawn |
|---|---|
| Semiconductor stacks | Si, Ge or GaAs; 1–3 layers n, p or intrinsic, doping 1e14–1e19 cm⁻³, 50 nm–10 µm each; no recombination, band-to-band or SRH (τ from 1 ns to 10 µs); sometimes lit; sometimes a Schottky contact |
| MOS capacitors | Si, Ge or GaAs, doped 1e15–1e18 cm⁻³; oxide 1–20 nm; Au, Al or Pt gate at any flat-band offset; with or without generation |
| Electrolyte cells | one to three salts from the data library's ions, 0.1 mM–1 M, 1 µm–1 mm; double layers resolved or strictly neutral; bath, reversible or blocking electrodes |
| Electrodes | Ag, Cu or Zn deposition, or Fe³⁺/Fe²⁺ on platinum, with Butler–Volmer kinetics over four decades of k₀; with or without a supporting salt; against a second electrode or a bath |
| Liquid junctions | two mixed solutions, 1 µM–1 M, sharp or graded, at open circuit |

Each device goes through what a user would do with it: a cold solve at no bias, and one at a
bias; a warm sweep to that bias; open circuit, where it's lit; a current drive (a fraction of
what passes at 1 V); a transient after a voltage step; and the impedance about equilibrium,
from far below its slowest relaxation to 1 GHz. Every result is judged by `check()`, and by what must hold whatever the device:
levels flat at equilibrium with no current, warm and cold solves agreeing, the current driven
being the one passed, the impedance passive (to 1e-3 of |Z|) and, at low frequency, equal to
the steady dI/dV,
and a MOS capacitance never above its oxide's. Each case is seeded from its family and index, so
any failure reruns alone (`npm run stress -- semi 79`).

At 500 devices per family (about three minutes on a desktop), 11,380 of 11,381 scenarios pass:

| Family | Cold at 0 V | Cold at bias | Warm sweep | Open circuit or current | Transient | Impedance |
|---|---|---|---|---|---|---|
| Semiconductor stacks | 500/500 | 500/500 | 500/500 | 105/105 | 500/500 | 395/395 |
| MOS capacitors | 500/500 | 500/500 | 500/500 | | 500/500 | 500/500 |
| Electrolyte cells | 500/500 | 500/500 | 500/500 | | 500/500 | 381/381 |
| Electrodes | 500/500 | 500/500 | 500/500 | 500/500 | 500/500 | 500/500 |
| Liquid junctions | | | | 500/500 | 499/500 | |

(Impedance only where the device has an equilibrium to linearise about; open circuit only for
lit cells.) The one that fails is among the known limits below: a junction on a grid far too
coarse for its double layers.
`bench/stress.json` keeps the summary.

Most of what driftlet's numerics do differently began as a failure here: the impedance's
current read across the quietest cut of the device rather than at the contact, GMRES judged by
its true residual, blocked species held flat in steady solves, continuations for a current drive
and for light, and Newton refining its solves where a factorisation loses a mode, with an exact
$`J \cdot v`$.

## Hard cases

`npm run hard` keeps a corpus of named cases that were hard once (sharp 3 M | 1 µM junctions on
short steps, an inversion layer filling over minutes, a solar cell's open circuit on a 0.25 nm
grid, corroding drops), with their step and factorisation counts against a committed baseline,
so a change that makes one worse shows.

## Known limits

- **A population held only through couplings far below double precision** is solved in steady
  state by refining Newton's solves, as far as its own equations resolve it: a bipolar stack's
  floating base at bias (a p⁺ base whose holes are held ~1e14 more weakly than they move within
  it) is, and a closed Fe³⁺/Fe²⁺ cell driven by a current near its limit, where Fe³⁺ falls 20
  orders below Fe²⁺, is found as the held voltage that passes it. Where a solve fails, its warning
  says where it lost its digits (a transient gets there). A region held only through its faces (by
  conductances 1e-24 of its own, or by neighbours conducting 1e20 times less) is solved directly
  through its balance summed over it, where its own fluxes cancel exactly
  ([islands](numerics.md#steady-state)): its current comes out right to ~1e-15, where it was
  up to 16× off. A population held that weakly within a region, or fed by a reaction, isn't
  covered, and its solve can still converge to the wrong level; the terminal currents then
  don't add up, which is the thing to check.
- **The impedance is resolved to about 1e-4 of |Z|** where its solves are hardest: the real
  part of a nearly ideal capacitor (a closed redox cell with almost none of one partner), 1e-7 of
  |Z|, can come out slightly negative. Its DC limit matches the steady dI/dV in every stress case.
- **Strictly neutral regions on very short steps** lose digits at a neutral face between very
  different solutions (3 M against 1 µM at steps of 1e-11 s).
- **An unresolved double layer** (a cell coarser than the Debye length) gives a charge that
  depends on the grid, and in a sharp junction between very different solutions can stop a
  transient; the solution warns. So does a steep profile on too coarse a cell.
- **A surface covered to a bare fraction below ~1e-10** barely affects anything any more, and
  the solver struggles with it (the hard corpus's passivating drop without a double layer stops
  partway through its hour); the solution warns, and suggests a less stable film.
- **Currents below ~0.1 nA/m²** (a GaAs junction's leakage without recombination) aren't resolved
  between solves; such a current is zero for any practical purpose.
- **A failed solve says so:** `converged: false`, with a warning naming where the system lost
  its digits, or which current no voltage reaches. In these tests, every steady state and
  transient reported converged has passed `check()` and the invariants above.
- What driftlet doesn't model at all (more than one dimension, cross-diffusion, heat, optics,
  field-dependent mobility) is in the [README](../README.md#what-it-doesnt-do), and open
  numerical work in the [roadmap](../ROADMAP.md).
