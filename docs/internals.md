# Internals

A map of the source for anyone reading or reviewing it: what each module is for, the layout of
the solver's state that the rest assumes, and which modules state a contract. The user-facing
reference is the other docs; this one is about the code.

## The rule: a module states its contract

A module whose correctness can't be seen from its call sites says at its head what it reads,
what it returns or writes, and what must hold of the result, so that it can be checked on its
own, without the rest of the library in mind. Where the contract and the code disagree, one of
them is wrong: fix both, with a test, in the same commit.

| Module | Contract |
| --- | --- |
| `src/laws.js` | yes: the steady system's conservation laws |
| `src/lawRows.js` | yes: the laws written into the Newton system, and the bordered solve |
| `src/timeSteps.js` | yes: one step (backward Euler, BDF2), advance's halving, integrate's error control |
| `src/start.js` | yes: the cold start, from each region's c0 and the contacts' and ports' levels |
| `src/newton.js` | yes: Newton's iteration, its scaled rows, GMRES-refined solves and the islands' round |
| `src/assemble.js` | yes: the equations assembled, the conventions of every row, and the order of the terms |
| `src/regions.js` | yes: statistics, storage, fluxes along segments, bulk reactions, the neutral nodes' charge rows |
| `src/boundaries.js` | yes: faces (φ laws, links, reactions, gates), ports, contacts, circuit rows |
| `src/steady.js` | yes: steady solves, their stages and continuations, the pseudo-transient ramp |
| `src/difference.js` | header: the dilute kernels in difference form |
| `src/blockTridiagonal.js` | header: the block Thomas solver |
| `src/statistics.js` | header: each model's `evaluate` / `invert` |
| `src/solver.js` | not yet (below) |

## From a definition to a solution

- `src/device.js`: `normalizeDevice(def)` validates a plain definition and returns the
  normalised `model` (species, materials, regions, interfaces, contacts, ports, terminals, the
  grid from `src/grid.js`). Every later stage reads the model and never the definition.
- `src/index.js`: `Device`, the public object. It holds the model, builds a `Solver` on first
  use, and turns solver states into `Solution`s (`src/solution.js`).
- `src/solver.js`: `Solver`, the numerical core (below).
- The kit (`src/kit.js` and what it re-exports: stack, data, sources, channels, checks,
  describe, traces, live, recorder, particles, polarization, junction, devices) works on
  definitions and solutions only, never on a solver's internals. `src/plot.js` draws
  solutions.

## The solver's state

- **Blocks.** One per grid node and one per face between regions: node g is block
  `blockOfNode[g]`, face f's flux block is `blockOfFace[f]`, in order along x.
- **Slots.** Each block has `M = 1 + n + extra` slots: slot 0 is φ̂ (at a face its
  displacement D, at a metal's node its carriers' segment flux J), slots `1 … n` the species'
  η = μ̄/RT, then a face's reaction rates and gates' fractions, or at an electrode port's nodes
  its surface's coverages (as η) and a membrane port's gates. Row and unknown `o = b·M + slot`.
- **Compact indices.** Only active slots enter the linear system: `loc[o]` is a slot's place in
  its block (−1: inactive), `rix[o]` its index in compact vectors, with `rix = N` (the sink) for
  inactive ones. `res`, `delta` and the Jacobian (`sys`, a `BlockTridiagonal`) are compact.
- **Unknowns** are `u + uLo`, a double-double, so that differences of η keep their digits; the
  Jacobian and solves are plain doubles. (η is read as `u + uLo` wherever a difference of levels
  is taken; φ̂ in the statistics and the displacement, a face's own slots and the gates as `u`
  alone.) Concentrations `c[g·n + i]` and coverages
  `th[g·nSurf + q]` follow from `u` (`computeConcentrations`).
- **Terminals:** 0 the left contact, 1 the right, `2 + k` port k. A floating terminal (driven by
  a current, or behind a resistance) has its voltage as an extra unknown and its circuit row
  beside the block system, solved by bordering (`solveBordered` in `src/lawRows.js`: `termB`
  its column, `termC` its row, `termDI` and `termDX` its own and other terminals' ∂/∂V).
- **Stretches:** each connected run of regions where a species is present
  (`_findStretches`), the unit that conservation counts in.

## The solver's sections

By method, in file order:

- **Setup:** the constructor (layout, per-node data, φ-undefined regions, terminals),
  `_activeSlots`, `_findStretches` (stretches, flat stretches), `_laws` and `_structure` (the
  laws, from `src/laws.js`).
- **Applying the laws** is `src/lawRows.js` (`captureLaws`, `captureLocal`, `snapFlat`,
  `applyLaws`, `solveBordered`), called from `assemble`, Newton and the steady solves.
- **Starting state** is `src/start.js` (`initFromComposition`, which the Solver's own calls);
  `_surfaceStart`, which the API uses too, stays in the solver.
- **Assembly** is `src/assemble.js` (`assemble`, `_assembleDifference` with `src/difference.js`,
  `_assembleBookkeeping`, and `_j` and `_replaceRow`, which write the Jacobian), with its kernels
  in `src/regions.js` (`computeConcentrations` and the statistics, nodes, segments, bulk
  reactions and their rates, `_chargeRows` and `_untransform`) and `src/boundaries.js` (faces,
  `_terminals`: ports, contacts, circuit rows). Their methods are the Solver's own (installed on
  its prototype), since every term reads and writes its arrays.
- **Newton** is `src/newton.js` (the iteration, its row scaling, the GMRES-refined solves, the
  islands' second round), which the Solver's own `newton` calls; `_factor` stays in the solver.
- **Time steps** are `src/timeSteps.js` (`step`, `advance`, `integrate`, the predictor and the
  error norm), which the Solver's own `step`, `advance` and `integrate` call; `_snapshot` /
  `_restore`, `_warmFrom` and `_probes`, which the API uses too, stay in the solver.
- **Small signal:** `impedance` and the helpers after it.
- **Bookkeeping:** `segmentFluxes`, `faceOneWay`, `_assembleBookkeeping`, the contacts' intake.
- **Steady solves** are `src/steady.js` (`solveSteady`, `_steadyFromHere`, the continuations:
  floating, generation, a contact's voltage; `_steadySteps`, the pseudo-transient ramp, and
  `_renormalizeSpectators`; `_directSteady`, `slowestTime`, `fastestTime`, `stuckReactions`),
  methods of the Solver like the assembly's.

The rest of `src/solver.js` is to be split along these lines, each piece with its contract.

## Reviewing one module

A review of one module reads its contract, its code and its tests, and nothing else unless the
contract points there. Each claim cites `file:line` with the line quoted, so it can be checked
mechanically, and a finding comes with the device or script that shows it. A claim about code
not in view is a question, not a finding. Two instruments check a change to the core beyond the tests'
answers: `test/jacobian.test.js` (every row's derivatives against finite differences of its
residual) and, for a move or refactor meant to change nothing, `tools/fingerprint.js`: a hash of
the solver's state after every solve in the suite (on a fixed clock, so runs stopped by a time
budget repeat), recorded before and after and compared. `tools/check-hook.js`, loaded into any run (`node --import
./tools/check-hook.js …`: the tests, the law sweep, a fuzzer), checks the laws' cheap invariants
as the solver goes: no row taken twice, none inactive, no stray ∂/∂V in a law's row, findLaws'
structure unwritten.
