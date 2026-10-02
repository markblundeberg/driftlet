# AGENTS.md: driftlet

Guide for coding agents working in this repo.

- **Start with the docs:** [README](README.md), [docs/conventions.md](docs/conventions.md),
  [docs/device.md](docs/device.md), [docs/statistics.md](docs/statistics.md),
  [docs/numerics.md](docs/numerics.md), [docs/kit.md](docs/kit.md),
  [docs/visualization.md](docs/visualization.md), [docs/data.md](docs/data.md) and
  [ROADMAP.md](ROADMAP.md). Keep them true as the code changes, since they're the reference. Math
  in the docs is LaTeX for GitHub: inline as $`…`$ (the backticks keep markdown's emphasis rules
  out of it; a bare $…$ breaks wherever a `_` or `*` pairs up), displays in ```math blocks; the
  README stays Unicode, for npm.
- **Constraints that don't bend:** pure-JS ES module, zero runtime dependencies, 1D,
  block-tridiagonal Jacobian, runs in a Web Worker, loadable from jsdelivr.
- **Tests are the deliverable.** Every physics feature lands with an analytic validation test
  (the validation suite). Use node's built-in test runner (`node --test`), no test dependencies.
- **Commits:** small and logically scoped. End messages with the Co-Authored-By trailer only;
  never a `Claude-Session:` line (a local commit-msg hook strips it as a backstop). Commit and
  push freely. Releases: add a CHANGELOG.md entry, bump `version` in package.json (and the CDN
  pins in the README and llms.txt, which a test checks), commit it as the bare version, tag
  `vX.Y.Z` and push with `--follow-tags`. The release workflow stages it on npm through trusted
  publishing, and it goes public only when the owner approves it with 2FA, so the owner decides
  when. Primary branch is `master`.
- **Hook setup after a fresh clone** (`.claude/` is gitignored, local infra):
  `git config core.hooksPath .claude/hooks && chmod +x .claude/hooks/commit-msg`
- **Background:** driftlet spins off the owner's ESBD web book (https://marklundeberg.com/esbd/,
  "electrochemical species band diagrams"), whose central idea is the species voltage `V_i =
  μ̄_i/(z_i F)`. driftlet doesn't require V_i (it's a display of μ̄_i, which also covers neutral
  species; the physics and the API are in μ̄), but it does insist on honest thermodynamics
  ([conventions](docs/conventions.md)): controls are μ̄ differences, φ is bookkeeping, band
  alignment is per interface.
