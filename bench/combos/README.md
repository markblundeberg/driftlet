# bench/combos: features meeting features

`bench/stress.js` puts a few fixed device families through what a user does with them. The outside
reviews of November 2026 found that most faults come from features combined in ways no test had
tried (ROADMAP, "## 2.0", "The shape"). This harness goes after those faults. Its generators each
build a device from independently chosen features, and its properties compare two routes to one
answer.

## A family

`families/NAME.js` default-exports `{ name, make(r) }`, where `r` is the seeded generator in
`rng.js` (`u`, `log`, `pick`, `chance`, `int`). `make` returns `{ def, params, plan }`:

- `def`: a device definition. Plain data that survives `structuredClone`: no functions, no
  shared objects between two calls.
- `params`: what was chosen, for a failure's report, including which features are on
  (`features: ['traps', 'graded', …]`).
- `plan`:
  - `terminal`: the name of a terminal held at a voltage, which the properties drive (`'right'`,
    or a port's name).
  - `V`: a bias (V) for steady solves at that terminal. It is absolute, so where the device's
    rest isn't at 0 V, add the rest voltage.
  - `dV`: a small step (V) for transients.
  - `equilibrium`: true when the definition, at its own drives, is at equilibrium (flat levels,
    no current).
  - `unique`: true when the steady state at any drive doesn't depend on history. That means no
    closed amounts, no terminal driven at `I: 0` or floating, and no blocked stretch. Some
    properties need it.
  - `slow`: an estimate in seconds of the slowest relaxation (L²/D, RC, 1/k). Transient properties
    run to many times it.

Check a family with `node bench/combos/try.js NAME 40` before fuzzing with it. Nearly every case
should build, solve cold, pass check() and solve at `plan.V`. A family whose devices fail cold
only finds the cold start's known limits.

## Symmetry properties

`gauge` shifts every applied voltage by the same amount. It compares concentrations,
terminal currents, and the electrochemical levels after subtracting the charge times that
voltage shift. `scale` dilates the length at fixed time, concentration and voltage. It changes
transport and electrostatic coefficients together, then compares the mapped steady state and
terminal currents.

`mirror` compares steady states on reflected grids. It swaps contact roles and maps the signs
of directional inputs. A terminal current defined into the device retains its sign at the
swapped contact. Radial geometry, a face with charge fixed to its right plate, and meshes
that do not reflect node for node are skipped. The property is opt-in while reflected cold
solves can refuse cases that solve in the original orientation.

```sh
node bench/combos/run.js --n 30 --only mirror,gauge,scale
```
