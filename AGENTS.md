# AGENTS.md: driftlet

Guide for coding agents working in this repo.

- **Start with [SPEC.md](SPEC.md).** It's the design brief: scope, physics conventions, numerics,
  API sketch, validation suite and milestones (M0 → M4). Improve it when it's wrong, and note
  why in the commit.
- **Constraints that don't bend:** pure-JS ES module, zero runtime dependencies, 1D,
  block-tridiagonal Jacobian, runs in a Web Worker, loadable from jsdelivr.
- **Tests are the deliverable.** Every physics feature lands with an analytic validation test
  (SPEC §6). Use node's built-in test runner (`node --test`), no test dependencies.
- **Commits:** small and logically scoped. End messages with the Co-Authored-By trailer only;
  never a `Claude-Session:` line (a local commit-msg hook strips it as a backstop).
  Commit freely; don't push or publish to npm unless the owner asks.
  Primary branch is `master`.
- **Hook setup after a fresh clone** (`.claude/` is gitignored, local infra):
  `git config core.hooksPath .claude/hooks && chmod +x .claude/hooks/commit-msg`
- **Background:** driftlet spins off the owner's ESBD web book (https://marklundeberg.com/esbd/,
  "electrochemical species band diagrams"), whose central idea is the species voltage
  `V_i = μ̄_i / (z_i F)`. driftlet does not require V_i, but it does insist on honest
  thermodynamics (SPEC §3 principles): controls are μ̄ differences, φ is bookkeeping, band
  alignment is per interface.
