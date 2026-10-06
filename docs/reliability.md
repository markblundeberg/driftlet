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
partitions, Goldman–Hodgkin–Katz, the limiting current and Butler–Volmer kinetics, the
Macdonald and finite-length Warburg impedances, Cottrell and spherical diffusion, and the
second law (the free energy in through the terminals equals what's dissipated, every term
non-negative). Against other codes: IonMonger's perovskite J–V hysteresis (finite elements,
Octave), and liquid-junction potentials against JPCalc, LJPcalc and JLJP. The README's
[validation table](../README.md#validation) lists each with its tolerance.

Each solution can check itself too: `check()` from `driftlet/kit` balances every species'
ledger (what comes in through each terminal, is made or destroyed), the conservation
bookkeeping of a transient, and the grid (the device solved again on a grid twice as fine).

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
  conductances 1e-24 of its own, or by neighbours conducting 1e16 times less) is solved directly
  through its balance summed over it, where its own fluxes cancel exactly
  ([islands](numerics.md#steady-state)): its current comes out right to ~1e-15, where it was
  up to 16× off. Past that, a solve can still converge to the wrong level: behind neighbours
  conducting 1e20 times less (the face's own elimination loses the digits), or for a population
  held that weakly within a region, or fed by a reaction. The terminal currents then don't add
  up, which is the thing to check.
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
