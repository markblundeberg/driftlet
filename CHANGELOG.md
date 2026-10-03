# Changelog

driftlet follows semantic versioning; while it's 0.x, a minor version may change the API.

## 0.3.0 (2026-10-02)

- `c0` profiles: a region's initial concentration can be tabulated against x (`{ x, c }`,
  piecewise linear), so a transient starts from a state it didn't simulate reaching (a packet
  already injected, a gradient laid down). Species without a `c0` start at their contact's
  level, with $`\phi`$ chosen node by node for local neutrality, so a counter-ion or majority
  carrier follows the profile. `describe()` warns when a profile misses its region.
- A `c0` given for a species a contact feeds is now its starting state (it was ignored, the
  species starting at the contact's level); without one, nothing changes.
- Validation: a Gaussian hole packet laid down as a profile drifting, spreading and decaying at
  the ambipolar rates, and a salt's diffusion mode decaying at $`\pi^2 D/4L^2`$.
- Results identical across JavaScript engines: powers go through `powi`/`powr` instead of `**`,
  whose last bit differs between Node 22 and 24 (enough to flip a step in a long adaptive run).
- Demos: a real-browser smoke test in CI (`npm run smoke`: every page in headless Chrome, every
  control moved to its extremes), and screenshots retaken by script (`npm run shots`), each
  page held at a set state (`?t=…`, and `?phase=…` on the impedance page).

## 0.2.0 (2026-10-02)

- Kit sources: waveforms (`pulse`, `square`, `triangle`, `ramp`), `injector()` (a port that
  spreads a driven current evenly) and `recombination()` (from a minority lifetime).
- `advance(t, { probes })`: a species' concentration or voltage at points inside the device, in
  the trace after every accepted step.
- `units` back out of SI: `toPerCm3`, `toMolar`, `toCm2PerS`, `toUm`, `toNm`.
- Validation: a Haynes–Shockley pulse against the ambipolar drift, spread and decay.
- Newton damps only the device's potentials, not floating terminal voltages (which enter
  linearly): a current-driven port switching off, its voltage collapsing by hundreds of volts,
  no longer crawls. A strong Haynes–Shockley pulse takes 55% fewer factorisations.
- Device reference: what a driven port's conductance actually does (it shares the current; a
  small one spreads it evenly), correcting 0.1.2's note.
- `llms.txt`: good habits (draw the closed form alongside, stay honestly in 1D).

## 0.1.2 (2026-10-02)

- The npm package ships `docs/`, so the guide's references to them work from the tarball.
- `llms.txt`: a transient template (a pulse injected through a current-driven port, sampled
  downstream), and what makes transients expensive.
- Device reference: what a driven port's conductance does, and that a device's time starts at 0
  when it's made (`solve()` doesn't advance it).
- `SPEC.md`, the original design brief, retired: it lives on in the docs and in git history.

## 0.1.1 (2026-10-02)

The first release published by trusted publishing from GitHub Actions, with provenance.

- README: what driftlet is for, links to the live demos, npm and docs, and every demo in the
  table of fields.
- Guide and demos: the current a species carries flows down its species voltage $`V_i`$, so
  electrons and anions themselves move uphill (they said every species drifts down); $`V_i`$ is
  described as a display of $`\bar\mu_i`$, which also covers neutral species.

## 0.1.0 (2026-10-02)

First release: the solver (`driftlet`), the porcelain (`driftlet/kit`) and level-diagram plots
(`driftlet/plot`), with TypeScript declarations, the validation suite, and live demos.
