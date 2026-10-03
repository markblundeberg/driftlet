# IonMonger reference scans

J–V hysteresis scans of a planar perovskite solar cell computed with
[IonMonger](https://github.com/PerovskiteSCModelling/IonMonger) 2.0 (Courtier, Cave, Walker,
Richardson & Foster, J. Comput. Electron. 18, 1435 (2019); commit `1ae4f9e`), for
[`test/ionmonger.test.js`](../../ionmonger.test.js), which runs the same cell in driftlet.

Each `scan_<rate>.csv` (bulk SRH off) and `full_<rate>.csv` (IonMonger's full defaults, bulk SRH
on) has columns t (s), V (V), J (mA/cm², photocurrent positive); the first row's J is NaN. Rates
are in V/s. The `full_` set is also what the [perovskite demo](../../../demos/perovskite.html)
draws as IonMonger's.

**The cell.** IonMonger's default parameters (`parameters_template.m`): TiO₂ (100 nm) | MAPbI₃
(400 nm) with mobile iodide vacancies over an equal immobile background | spiro-OMeTAD (200 nm),
298 K, Boltzmann statistics, light (1.4e21 m⁻² s⁻¹, α = 1.3e7 m⁻¹) entering through the TiO₂,
interface SRH at both faces; bulk SRH (midgap, τ_n = 3 ns, τ_p = 300 ns) is on in the `full_`
set and turned off (`tn = tp = 0`) in the `scan_` set.

**The protocol.** Steady state at 1.2 V under light, then 1.2 → 0 → 1.2 V at the given rate
(IonMonger's own option of starting from a steady state, in place of its default tanh
preconditioning ramp from V_bi).

**How they were made.** GNU Octave 8.4, with compatibility changes to IonMonger that leave its
equations alone: an `optimoptions` shim onto `optimset`; the `bvpinit` guess used directly
(Octave has no `bvp4c`; `fsolve` refines it); `Vectorized` dropped; a consistent initial slope
for the DAE (solving M y′ = F with the algebraic rows pinned) in `FE_solve.m` and
`precondition.m`; a transpose of the state passed between splits. Those changes were worked out
by a Claude (Opus 5.5) agent building a perovskite demo with driftlet, which first ran this
comparison; these files were regenerated from its changes, reproducing its numbers exactly.
IonMonger is AGPL-3.0, so none of its code is here: only its output.
