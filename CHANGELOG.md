# Changelog

driftlet follows semantic versioning; while it's 0.x, a minor version may change the API.

## Unreleased

- `set()`: a contact given `V` or `I` drops the other, so `{ I: 0 }` switches a held contact to
  open circuit (it used to need `V: undefined` alongside).
- A face reaction's α may be 0 or 1 (plain mass action one way), not only strictly between.
- Each solution's `interfaces[f]` names the regions it joins, `left` and `right`.
- `describe()`: a region's Debye length counts the carriers balancing its doping, not only what
  `c0` lists (a trace in `c0` made a doped layer's read in µm instead of nm).

## 0.4.0 (2026-10-03)

- A bulk reaction's `kf` can vary with position: a profile `{ x, values }` against the device's
  x, averaged over each node's box, so the total rate is exact on any grid.
- `photogeneration({ material, flux, alpha, mu, from, to })` in the kit: Beer–Lambert absorption
  as a generation reaction. Validation: J_sc in an n⁺p cell against collection theory, from blue
  light lost in the emitter to red collected deep in the base.
- `photogeneration({ makes })`: light can make something other than e⁻ + h⁺, such as an
  exciton. Validation: an organic bilayer, with excitons as a third, neutral species diffusing to
  the donor/acceptor interface and splitting there (J_sc against exciton-diffusion theory), and
  V_oc from charge-transfer recombination by detailed balance.
- Solar demo: a wavelength control (silicon's α from Green 2008), the spectral response against
  collection theory, and where the light makes pairs and which are collected.
- Demo: an organic bilayer solar cell, where light makes neutral excitons that diffuse to the
  donor/acceptor interface and split; J_sc against donor thickness and V_oc by detailed balance,
  each against theory.
- The cold start of an undoped layer holding one carrier (which can't be neutral) keeps the
  running φ instead of driving it to an extreme, and φ starts from a pinned or gated left
  contact. Such devices used to fail on the first Newton step.
- A species fed only by reactions still needs a `c0`, and the error now says that it's only a
  starting point, not an amount it keeps.
- Breaking: a `c0` profile is `{ x, values }` (it was `{ x, c }`), the same shape as a `kf`
  profile and, against t, a waveform.

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
