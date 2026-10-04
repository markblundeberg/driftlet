# Changelog

driftlet follows semantic versioning; while it's 0.x, a minor version may change the API.

## Unreleased

- Membranes as faces: a species law `{ type: 'permeability', P }` (m/s), electrodiffusion through a
  thin membrane in a constant field (Goldman–Hodgkin–Katz), exactly zero where μ̄ is level. With a
  capacitive φ law, a cell membrane is one face. Validated: the GHK resting potential with the
  solutions strictly neutral or their diffuse layers resolved, the GHK current–voltage curve, the
  agreement with a resolved 5 nm lipid layer, and charging at the GHK current over C.
- Saturating kinetics for face reactions, `vmax` and `K` in place of `k0` and `alpha`:
  r = vmax Π (c/(c + K))^|ν| (1 − e^(−A/RT)), Michaelis–Menten in each substrate and exactly zero
  at A = 0, for pumps and transporters. The Na⁺/K⁺-ATPase on a membrane face, in a closed cell,
  gives the Mullins–Noda potential, stalls at its static head (3Δμ̄_Na − 2Δμ̄_K = ΔG_ATP) with
  nothing leaking back, and leaves Donnan equilibrium when off.
- `check()` keeps a ledger per compartment: a species a face lets through only by a law or a
  reaction gets one on each side, with what crosses as a term (in a cell: Na⁺ leaks in as fast as
  it's pumped out). `describe()` lists what crosses a face by a law, and saturating kinetics.
- A charged-pore demo (an ion channel resolved along its axis): Donnan layers at each mouth, K⁺
  selectivity, the conductance plateau at low salt against Teorell–Meyer–Sievers, the zero-current
  voltage across a salt gradient, and rectification by a one-sided or bipolar charge. Validated in
  test/pore.test.js (TMS to 1% in a long pore, its end correction ∝ 1/L, rectification).
- A resting-potential demo: a closed cell's leaks and the Na⁺/K⁺ pump on a membrane face, V_m
  against Mullins–Noda, each ion's species voltage stepping at the membrane by its driving force
  V_m − E_i, and the run-down to Donnan equilibrium over hours when the pump stops.
- Probes read φ (`{ x, quantity: 'phi' }`), so a transient's trace can carry a membrane potential.
- A floating contact starts where the start's own composition puts it, so a cold-start
  transient's first trace point is right (it read the held terminal's voltage, off by a Nernst
  term or more).
- Docs: a cell-membrane template in the agent guide; a bath terminal reads its reference
  species' level, not φ; which contact a species without `c0` starts from; `set()` from
  `onsolution`; the worker from a CDN; the data's 25 °C. All from the membrane agent's debrief.

## 0.6.0 (2026-10-03)

- `check(device, sol)` in the kit: checks on a solved device, each saying what it compared, as
  `{ ok, items, text }`. Convergence; warnings and likely unit slips (marked `?`, a prompt to
  look); in a steady state, each species' ledger (what every terminal and reaction brings in or
  uses, summing to zero: J = F(G − R) under light, and where each carrier went), and the device
  solved again on a grid twice as fine (current, floating voltages and region charges to 1 %, the
  change estimating the grid's own error); in a transient, conservation. Validated on a lit
  silicon cell (the ledger term by term, and the grid change against a 16× finer solve) and run
  over every steady solve in the test suite without false alarms. The solar demo shows it, and
  its grid is finer for it (hmin 1 nm had left the emitter's depletion charge 2 % off).
- Solutions say whether they're a steady state, `sol.steady`.
- Curated devices in the kit: `perovskiteCell({ V, light, bulkSRH })`, IonMonger's default
  perovskite solar cell as a plain definition, and `PEROVSKITE_SCANS`, IonMonger's own J–V scans
  of it (current at 0 V, each sweep's P_max and V_oc, hysteresis index, at seven rates, with and
  without bulk SRH), so a demo can start from benchmarked numbers and show its agreement;
  `hysteresis(trace, turn)` reads the same figures from a run. The IonMonger test and the
  perovskite demo use them. Asked for by the perovskite agent.
- `recorder(device, { every, times, probes })` in the kit: a transient run in pieces (or in
  animation-frame budgets, step for step the same run), with the trace joined up from the start
  and whole solutions kept as frames at set times, `frame(t)` to scrub. Validated on a
  diffusion mode decaying frame by frame at its analytic rate.
- The IonMonger benchmark measures the loop as a current difference where it's gentle and as a
  voltage offset where it's steep: within 0.010 mA/cm² and 0.31 mV at all seven rates (the
  0.06–0.23 mA/cm² reported before was a sub-millivolt shift read on a near-vertical curve). The
  perovskite agent spotted this in its own comparison.
- Each solution reports its bulk reactions' rates, `bulkReactions[k]`: per node, per region and in
  all, so charge balance J = F(G − R) is a few lines (validated to 1e-6 in the IonMonger cell).
- `describe()`: an SRH reaction is listed where it runs; a region's Debye length counts only the
  mobile charge `c0` leaves unbalanced on top of `c0` itself (0.5.0 counted a perovskite's
  vacancies twice, once as themselves and once as balancing their background).

## 0.5.0 (2026-10-03)

- SRH (trap-assisted) recombination as a rate law, without trap species: in the bulk,
  `srh: { material: { tauN, tauP, n1 } }` in place of `kf`; at a face, `srh: { vn, vp, n1 }` in place
  of `k0` and `alpha`, with n and p each from its own side. Equilibrium stays exact (p₁ comes from
  the state), and at a face it saturates at one carrier's capture, as interface recombination in
  perovskite cells does. Validated against the bulk closed form, and at a face against the law
  and the device's bookkeeping.
- A benchmark against an independent code: perovskite J–V hysteresis (mobile iodide vacancies,
  interface SRH, with and without bulk SRH) against IonMonger at seven scan rates from 1 mV/s to
  1 kV/s. Hysteresis index within 5e-4, maximum power within 0.02 mW/cm², V_oc within 1 mV, whole
  loops within 0.06 mA/cm² (0.23 at the fastest-changing sweep). Three species and two face SRH
  laws; each scan takes 0.3 s against IonMonger's 50.
- Immobile reacting species (trap states, fixed charge states) in steady solves: a combination
  of them that reactions conserve is now held node by node, as a row of that node's block, so
  the steady solve goes direct instead of through giant time steps (which stalled when the
  traps relaxed slowly). Validated: explicit immobile traps recombine at exactly the SRH rate
  in steady state, their totals kept at every node.
- A cold start with slow ions alongside fast physics (a perovskite's vacancies with carriers and
  charged traps) solves: the steady solve's fallback time steps, and `advance()`'s first step,
  shrink as far as the device's fastest time scale, not only to a fraction of the slowest (or
  of the run's length). Both used to give up a few decades short. Reported, with a minimal
  repro and the diagnosis, by an agent building a perovskite demo.
- `traces()` and `bandDiagram()` take `energy: true`: the familiar band diagram, energy up (E_c,
  E_v and the quasi-Fermi levels in eV). The semiconductor demos (pn, solar, MOS, organic) now use
  it; species voltages stay where ions share the diagram (electrochemistry, the perovskite's
  vacancies). The docs say which picture to use when.
- Demo: a perovskite solar cell's J–V hysteresis from mobile iodide vacancies (IonMonger's default
  cell), with a scan-rate control, the vacancies' layers at each face, the steady and frozen
  limits, and the hysteresis index against scan rate with IonMonger's own points.
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
