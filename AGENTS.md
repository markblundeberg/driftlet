# AGENTS.md: driftlet

Guide for coding agents working in this repo.

- **Start with the docs:** [README](README.md), [docs/conventions.md](docs/conventions.md),
  [docs/device.md](docs/device.md), [docs/statistics.md](docs/statistics.md),
  [docs/numerics.md](docs/numerics.md), [docs/kit.md](docs/kit.md), [docs/visualization.md](docs/visualization.md), [docs/data.md](docs/data.md) and
  [ROADMAP.md](ROADMAP.md). Keep them true as the code changes, since they're the reference.
  SPEC.md is the original design brief, being retired: consult it for history, but don't cite
  it from code, tests or docs.
- **Constraints that don't bend:** pure-JS ES module, zero runtime dependencies, 1D,
  block-tridiagonal Jacobian, runs in a Web Worker, loadable from jsdelivr.
- **Tests are the deliverable.** Every physics feature lands with an analytic validation test
  (the validation suite). Use node's built-in test runner (`node --test`), no test dependencies.
- **Commits:** small and logically scoped. End messages with the Co-Authored-By trailer only;
  never a `Claude-Session:` line (a local commit-msg hook strips it as a backstop).
  Commit and push freely. Releases: bump `version` in package.json (and the CDN pins in the
  README and llms.txt, which a test checks), commit it as the bare version, tag `vX.Y.Z` and
  push with `--follow-tags`. The release workflow stages it on npm through trusted publishing,
  and it goes public only when the owner approves it with 2FA, so the owner decides when.
  Primary branch is `master`.
- **Hook setup after a fresh clone** (`.claude/` is gitignored, local infra):
  `git config core.hooksPath .claude/hooks && chmod +x .claude/hooks/commit-msg`
- **Background:** driftlet spins off the owner's ESBD web book (https://marklundeberg.com/esbd/,
  "electrochemical species band diagrams"), whose central idea is the species voltage
  `V_i = μ̄_i / (z_i F)`. driftlet does not require V_i, but it does insist on honest
  thermodynamics (SPEC §3 principles): controls are μ̄ differences, φ is bookkeeping, band
  alignment is per interface.
