# AGENTS.md: driftlet

Guide for coding agents working in this repo: one preparing a pull request, or one changing a
fork of driftlet for its own use.

- **Start with the docs:** [README](README.md), [docs/conventions.md](docs/conventions.md),
  [docs/device.md](docs/device.md), [docs/statistics.md](docs/statistics.md),
  [docs/numerics.md](docs/numerics.md), [docs/kit.md](docs/kit.md),
  [docs/visualization.md](docs/visualization.md), [docs/data.md](docs/data.md),
  [docs/reliability.md](docs/reliability.md), [docs/templates.md](docs/templates.md) and [ROADMAP.md](ROADMAP.md). Keep them true as the code changes, since they're the reference.
  For the source itself, [docs/internals.md](docs/internals.md): a module that states a contract at its
  head (src/laws.js) is changed with its contract and a test in the same commit. Math
  in the docs is LaTeX for GitHub: inline as $`…`$ (the backticks keep markdown's emphasis rules
  out of it; a bare $…$ breaks wherever a `_` or `*` pairs up), displays in ```math blocks; the
  README stays Unicode, for npm.
- **Constraints that don't bend:** pure-JS ES module, zero runtime dependencies, 1D,
  block-tridiagonal Jacobian, runs in a Web Worker, loadable from jsdelivr.
- **Tests are the deliverable.** Every physics feature lands with an analytic validation test
  (the validation suite). Use node's built-in test runner (`node --test`), no test dependencies.
  The demos have a browser smoke test, `npm run smoke` (headless Chrome, also in CI): run it after
  touching a demo or anything they import. `npm run shots` retakes the screenshots (each page at
  a set state, on virtual time, so unchanged pages come out byte-identical). A figure in the docs
  is drawn by the example above it (a block fenced ```` ```js figure=NAME ````, its picture at
  docs/figures/NAME.svg): after changing one, `npm run figures` and look at what changed.
- **Contributing:** work on a branch and open a pull request against `master`; the owner reviews
  it and cuts releases. Keep commits small and logically scoped. Before you open it, run
  `npm test`, and `npm run smoke` if you touched a demo or anything one imports. CI also runs
  `npm run bench -- --check` (its deterministic counters, matrix factorisations and the like,
  against bench/baseline.json), `npm run types` with a strict check of the declarations, and
  `npm run stress -- --check`. A CDN pin in the README and llms.txt follows `version` in
  package.json, and a test checks it: leave both to the release.
- **Background:** driftlet spins off the owner's ESBD web book (https://marklundeberg.com/esbd/,
  "electrochemical species band diagrams"), whose central idea is the species voltage `V_i =
  μ̄_i/(z_i F)`. driftlet doesn't require V_i (it's a display of μ̄_i, which also covers neutral
  species; the physics and the API are in μ̄), but it does insist on honest thermodynamics
  ([conventions](docs/conventions.md)): controls are μ̄ differences, φ is bookkeeping, band
  alignment is per interface. V_i is unorthodox, so don't push it: a device of electrons and holes
  alone gets the familiar energy-up band diagram (`energy: true`), and V_i comes in where it shows
  something that picture can't (ions, equilibrium between species), as an extra tool rather than a
  replacement.
