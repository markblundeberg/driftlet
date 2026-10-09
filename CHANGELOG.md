# Changelog

driftlet follows semantic versioning; while it's 0.x, a minor version may change the API.

## Unreleased

- A terminal driven by a current that alone feeds what it keeps (a floating gate's metal, a
  host's electrons behind an open circuit) reads its current in a time step as that amount's
  change, where that's the quieter reading. A floating gate at rest took 5000 steps per 100 s,
  half rejected, its charging current lost below its metal's round-off (it now takes the 44 a
  held gate does); a closed host at open circuit drifted from its charge (to 0.3995 V from
  0.4 V by 1e10 s), and its steps failed by 1e14 s. A terminal behind a resistance reads its
  current the same way: a gate charged through one was 0.9% off its RC response, silently, and
  took 3280 steps for 5τ, half rejected; it now follows backward Euler's exactly, in 60.
- A steady solve from a transient's state that Newton can't take goes on as a cold device would,
  from the layout for its drives: a lit MOS capacitor a nanosecond after its gate jumped failed
  where a cold solve takes 27 iterations; it now takes 32. A continuation no longer starts from
  the last steady solve's voltages once a transient has moved the state.
- `set({ grid })` carries the state over, interpolated onto the new nodes, where the regions
  and terminals are the same: a floating gate kept nothing of its charge through a regrid.
- `check()` leaves a contact that blocks a species out of that species' ledger: the flux it read
  there was the end box's own imbalance, and counted as a source it had balanced the books of a
  state that wasn't steady.
- The stress test solves each case again with its species listed in reverse, and has a family of
  reacting solutions in several regions, with an electrode at no current: before 0.16.1's fixes
  it failed 237 of 1496 checks.

## 0.16.1 (2026-10-08)

Found by a round of agents trying the new steady solves cold (GITT on an insertion host,
floating gates, closed thin-layer cells), whose results otherwise matched theory, a closed
battery's rested voltage to 1e-11 V:

- A steady solve of a device driven by a waveform read it at its pseudo-transient steps' clock,
  where that ramp ran: a triangle wave's steady state came out at its value seconds later (1 V
  instead of 0), silently. It's read at the present time now.
- A level held flat (a species reaching one contact alone) is flat exactly, every node at its
  first's. Flat only to round-off after a jump, a MOS capacitor's bulk holes, behind a floated
  gate with the back contact jumped to −6 V, had read as a current 1000 times the true leakage,
  silently: the terminal currents didn't add up.
- A film whose ions are all given by `c0` starts its electrostatic potential balanced against
  its electrode's surface species too. An adsorbed couple had started 1.4 V off: cold solves
  failed or took the pseudo-transient ramp, and `advance()` from the start stalled.
- An insertion host's occupancy is inverted exactly on a steep tabulated OCV (a 60 mV step over
  Δx ≈ 0.1): the root solver zigzagged across the step's inflections, ran out of iterations at
  points that weren't roots, and time steps failed down to 1e-14 s.

A second agent compared steady solves with long transients across random devices, and found two
more, both silent, both now matching to round-off in every species order:

- Conserved amounts overlapping (two complexations sharing a species, A + B = AB and
  A + C = AC) could stand in for rows that weren't independent, one balance said twice and
  another not at all, depending on the order the species were listed in: c_A 7% off, a phantom
  flux through a contact that blocks it, an I = 0 terminal passing 20 A/m². It hit 8–13% of
  random reacting multi-region devices. The rows are chosen by elimination now.
- An electrode port at no current whose reactions use up a closed ion (M⁺ + e⁻ = M(s)) keeps
  that ion together with the charge on its capacitance (or alone, with none). Unseen, the steady
  solve had no law for it and landed anywhere: c(M⁺) 2.2e8 mol/m³ where it's 10.

## 0.16.0 (2026-10-08)

The steady solve, reworked around one idea: every row it changes from a time step's is a sum of
balance rows, written out exactly (a conserved amount, a level held flat, an island's balance),
with its conservation laws now worked out under the terminals' drives. [Numerics](docs/numerics.md#steady-state)
describes it that way now.

Silent wrong answers, fixed:

- A closed battery at open circuit keeps its state of charge: a cold `solve()` had converged to
  another of its steady states (4.145 V where its `c0` said 3.958 V). What a terminal driven at no
  current alone feeds is conserved, at what the state holds when the solve begins, so after a
  transient the steady state is where that transient relaxes. An insertion host starts with its
  ion level at the electrolyte's (its φ is a gauge), so the floating terminal of a closed cell
  starts at its OCV, and transients from `c0` run at any current (the agent guide's workaround is
  gone).
- Adsorbates on an electrode that only turn into each other keep their total on every site; such
  a device had failed to build.
- A MOS capacitor's back contact moved under its biased gate: where its carriers recombine (so
  weren't held flat), its inversion layer couldn't follow through a scarce-minority bulk, and of
  500 random MOS capacitors 38 failed and 5 converged to a gate charge up to 2× off. Species that
  react only among themselves, each reaching the same one contact and nothing else, are now held
  flat at its levels: in a steady state they're at equilibrium. The stress test moves the back
  contact now too.

New:

- A gate or a capacitance left floating (`I: 0`) keeps its charge in a steady solve, where it
  threw. Driven at a current, a terminal that only stores what comes in (a gate, a host filling)
  makes `solve()` fail at once with a warning saying so, rather than throwing or marching a held
  voltage 20 V each way.
- Steady solves go direct (rather than through huge time steps) for a floating metal (its amount
  the charge on its faces), and for trap states holding a closed population's electrons.
- A blocked ion's level is held flat across faces too: a closed zinc cell past its limiting
  current, its electrolyte in two regions, failed at 0.9 V and now solves as one region does.

## 0.15.1 (2026-10-08)

From a further round of cold agents (an ion-selective electrode, a diode's ideality, Hebb–Wagner
polarisation, a Gärtner photoanode), all matching their theory:

- The new warning of a space charge reaching a contact no longer fires on a quasi-neutral layer
  split off a junction's (it read round-off as a cut-off layer).
- `check()` compares a current too small for its ledgers to call moving (a diode at low bias,
  1e-4 A/m², changed 3.6% on a finer grid while `check()` said ok).
- Electrons or holes above their band's density of states under Boltzmann statistics are warned
  of (degenerate; a photoanode driven past its couple's diffusion limit converged that way).
- `describe()` names a bath without a reference species as read through its φ, and flags an SRH
  n₁ given per cm³.
- Docs: an exchange or partition constant as a μ° difference; `dt0` for a transient's first
  instants; a contact holding several species fixes its end's composition; a bath's composition
  `set()` mid-run.

## 0.15.0 (2026-10-08)

Mostly from four more rounds of cold agents (fresh tasks, from the agent guide alone): sixteen
of them, every device matching its theory once set up, and the places where it silently didn't.

Silent wrong answers, fixed:

- `photogeneration()` generates nothing outside `[from, to]`: light given to one layer of a
  material that's also elsewhere (a p⁺-i-n⁺ diode of one silicon) generated in the layers beside
  it too, 10–25% too much photocurrent. And an αL of exactly 50 no longer fails.
- A probe at a face reads the side where its species is (it read an electrode's metal side,
  NaN), and asks for `region` where both sides have it.
- `impedance()` at a held terminal with every other one driven by a current (no way back for the
  signal) is an error saying to measure at a driven one; it returned round-off.
- A port's window at a device end node holds its species there unless that contact links them:
  a wire on a particle's surface node alone now charges it (it passed nothing).
- `set()` of a parameter mid-run reads held waveforms at the run's time at once (they read t = 0
  until the next step).

New:

- A bath can name its own material (`bath: { c, material }`, or the kit's `bath(c, ref, drive,
  { material })`): a solution straight against an ion exchanger or a gel, whose levels the end
  node takes through a Donnan step (tested against Donnan exclusion and the Teorell–Meyer–Sievers
  salt flux, to 1e-8). Without it, a bath is the end material at that composition, fixed charge
  and all, and the error now says so and what to do.
- `set()` of a parameter mid-run (a rate constant, a diffusivity) goes on at the step size an
  adaptive transient had reached, rather than restarting from a tiny one: 7× fewer steps for a
  page changing k0 every frame.
- `particles()` stays a sample at any step: a cell a dot would hop out of many times in one step
  is a sea for that step and then drawn afresh, so a long step no longer over-fills cells (a
  checkerboard) or costs in proportion. It draws as many cells as asked, as nearly as the nodes
  allow (it drew 52 of 80 on a 1 µm grid); `cap` is a density; dots sit within a cell by its
  concentration; a solution on another grid redraws the lattice; bad options are errors, and a
  NaN or negative step (a first frame) moves nothing.
- The membrane page draws its ions as dots at its foot: chloride thinning across the membrane,
  Donnan-excluded, against seas of Na⁺.

Warnings and errors that say what happened:

- A space charge that reaches a contact holding its end neutral (a double layer or depletion
  longer than its region, cut off there: σ 58% low, silently) is warned of.
- Beside a strictly neutral material, a capacitive face's charge sits in the cell next to it,
  a diffuse layer that wide in series with C, so a fine grid lowered the double-layer
  capacitance silently (to 0.036 F/m² from 0.2 with 0.1 nm cells): warned of past 2%.
- `describe()` says where a bath's SHE level sits against its terminal when the bath is read
  through a reference ion (1.37 V below a 0.5 M Cl⁻ terminal on table μ°); warns of a bath out
  of equilibrium with a bulk reaction at its contact (H⁺ and OH⁻ a little off K_w, a current at
  zero bias); gives an ion exchanger's fixed charge in mol/m³; and flags a large k0 as an
  exchange current only at an electrode.
- A steady solve that fails where an ion is too scarce for double precision (a minority swept
  out of a junction) says so, instead of suggesting a floating region.
- `advance(NaN)` does nothing and says so in its warnings; `advance()` of a non-number is an
  error. A capacitance-only port given a terminal species that no material holds is told it
  needs none.

Faster:

- Block elimination skips the coupling blocks' zeros (bit-identical results): the nerve demos'
  factorisations are ~15% faster. The benchmarks gain a squid-axon spike.

Docs:

- The agent guide gains ion-exchange membranes, excitons in an organic cell, and how to start a
  closed battery (hold it at its OCV, then drive it).
- How fine a grid must be where nothing warns: a transient's diffusion layer, a spike's front, an
  impedance's reach into a porous electrode, a layer depleted past a limiting current.
- A MOS capacitor without a port solves, and its impedance shows the high-frequency C–V above the
  minority supply's knee; the library's K_w; an OCV curve far flatter than ideal at an end; the
  double layer's voltage limit (~350 V_T); conservation's drift is relative.

## 0.14.0

- `particles()` in the kit (experimental): a solution drawn as dots of a fixed amount each,
  hopping between display cells at one-way rates taken from the solution, so their density is
  its concentration and their net crossing rate its flux (tested against Fick's slab, a decaying
  species and Ussing's flux ratio). Dense species are seas that hold no dots but exchange them,
  so a junction shows its majority carriers as seas and its minority carriers one by one. The
  pn demo draws its carriers this way, at the bottom, on the band diagram, and the axon demo
  the ions crossing its membrane each way, at their unidirectional fluxes.
- Solutions carry each species' flux across every segment (`flux`), its diffusivity (`D`), the
  temperature (`T`), how each species crosses each contact and face (`links`), a permeable
  face's one-way fluxes (`oneWay`), and each bulk reaction's one-way forward rate (`forward`)
  and stoichiometry (`nu`).
- A visual identity: the logo's wordmark redrawn in Recursive's casual sans (as outlines, so it
  looks the same everywhere), and the demo site in paper and ink with the logo's vermilion, a
  gallery whose sections show each field's cast, charges drifting along a potential at its foot
  (cations down it, anions up), and a stamp on each card naming what its page is checked against.
- A starting `c0` profile steeper than the grid's cells under it is warned of in every
  solution, with where and the cells that would resolve it: the start was silently the profile
  sampled at the nodes (a liquid junction's potential then seemed to grow as it formed).
- `check()` in a transient says its grid isn't checked (`?`), rather than leaving the grid item
  out: a transient's answer depends on the whole run, so it's rerun on a finer grid by hand.
- A time step's Newton iteration stops once its contraction shows the next update would be under
  tolerance, rather than taking that update to show it: about a third fewer factorisations where
  steps are capped by `dtMax` (the nerve demos: a bare axon's step 4.4 → 3.3 ms), 4–8% fewer in
  the benchmarks' adaptive transients.
- The agent guide is shorter (a map of the templates up top, tighter front matter, the
  level-diagram and mistakes sections halved), with notes from a third round of cold agents: a
  BJT's base contact, the semiconductor data, `vacuumDipole` with `build()`, a cyclic
  voltammogram from the stirred-electrode template.

## 0.13.0 (2026-10-07)

- A port can take its outside's composition, `bath: { c }`: its V is then the outside's φ, and
  each linked species' level there follows from the composition (an offset no longer has to be
  worked out). `hodgkinHuxley({ T, area })` in the kit writes a membrane port's links, Hodgkin
  and Huxley's linear conductances gated as theirs.
- The action-potential demo can put the membrane behind Frankenhaeuser and Hodgkin's 30 nm
  periaxonal space: K⁺ piles up there with each spike (about 3.5 mM) and clears over 50 ms, and
  E_K, the rest and the after-hyperpolarisation follow it, which Hodgkin and Huxley's fixed E_K
  can't. A test checks that such a space clears as one compartment, e^(−tP/θ), to 2e-3.
- A contact `bath` without a reference species, `bath(c, drive)` in the kit: its terminal voltage
  is its φ (an ideal salt bridge), so between two such baths V_right − V_left is the membrane
  potential a voltage clamp sets. `hodgkinHuxley({ T, linear: true })` writes Hodgkin and
  Huxley's own linear channels for a face (gated conductance links); clamped, the K⁺ current
  is g_K n∞⁴ (V − E_K) to 2e-3.
- A contact link for a surface velocity, `{ type: 'velocity', v }`: N_in = v (c_eq − c), between
  blocked and held. Thermionic emission over a Schottky barrier (v = A*T²/(F N_c)), within a
  few percent of J = A*T² e^(−φ_B/V_T)(e^(V/V_T) − 1), or a contact's surface recombination
  velocity: minority electrons through a base reach it as n₀(e^(V/V_T) − 1)/(W/D + 1/S) to 1e-5.
- Fixed: after a change of the left contact's voltage that the direct solve couldn't reach, the
  continuation (which ramped only the right contact, from where it was last solved) had nothing
  to ramp and returned the old state as converged. It ramps whichever contact moved now, and a
  ramp that runs no solve doesn't claim one.
- A solve that fails leaves the device as it was before it (it used to leave the failed iterate,
  from which nothing could go on: a battery host filling at a current, which has no steady
  state, couldn't then be advanced). An unreachable driven current's warning says when the
  opposite sign would have passed (a terminal's current is into the device), and when nothing
  passes at any voltage (a device that only stores charge: advance it in time instead).
- Fixed: a transient a second or more into a run, after a jump that a fine grid resolves in
  steps of 1e-13 s, stalled and stopped (a potential step at t = 1 s on a 5 nm grid): the error
  estimate was built from differences of absolute times, which keep only a few digits there.
  It uses the steps' own lengths now. An `advance()` that stops short, other than on its
  `budgetMs`, says where and why in its warnings.
- Fixed: a current-driven solve at I = 0 could find its voltage and still report no
  convergence (its last stage asked for a current of exactly zero). A solve that doesn't
  converge always carries a warning now.
- Fixed: with any port, `impedance()` read a contact's current at the contact, where behind a
  metal region it's a large conductance times round-off: a MOS capacitor with a metal gate and
  its channel held by a port read Re Z = −21.7 Ω·m² at its gate, at every bias. It's read across
  the device where no port's window intervenes now, and the back contact's nearly uniform
  response is split off as an open electrolyte's is.
- Fixed: a port holding a level inside another port's window left that port's ∂res/∂V in the
  row it replaced (a stale entry, ~1e-6 of the row: it touched only that port's impedance, or
  its circuit if floating). A gate's exponential rate law is held finite at absurd voltages,
  and `describe()` shows gated conductance links' gates.

## 0.12.0 (2026-10-07)

- Voltage-gated channels: a face's `gates`, each a fraction open with Hodgkin–Huxley kinetics in
  the voltage across the face (α and β in NeuroML's three standard forms), and permeabilities
  scaled by them (`gates: { m: 3, h: 1 }`). They're unknowns of the face's block; solutions report
  `interfaces[f].gates` and `.V`. `hodgkinHuxley()` in the kit writes the squid axon's channels,
  with GHK permeabilities matched to HH's conductances at rest. Under a voltage clamp each gate
  sits at α/(α + β) and relaxes at α + β; a squid axon's action potential follows the
  space-clamped HH equations with GHK currents to within 1 mV, and its impedance at rest (the
  inductance Cole measured, and a resonance near 60 Hz) the linearised equations' to 5e-4.
- Probes at a face: `{ interface: f, species }` (the flux through it) and `{ interface: f, gate }`.
- Gates on a port too: a membrane all along a region (an axon's), each node's gates following
  the voltage across the port's capacitance there, scaling its conductance links (Hodgkin and
  Huxley's linear channels, in species voltages). A squid axon's action potential propagates
  at 18.73 m/s, the cable equation's speed to 2e-3, out of the ions' drift along the axon
  (Hodgkin and Huxley computed 18.8). Face conductance links can be gated as well. A myelinated
  axon, nodes of Ranvier and internodes as regions with ports of their own, conducts
  saltatorily, each node's arrival within 1 µs of a compartmental cable's.
- A device with many held terminals (an axon's nodes of Ranvier, each a port) steps faster: only
  the floating terminals' ∂res/∂V and ∂I/∂x are carried through each assembly (all of them while
  `impedance()` runs).
- A demo, [myelin](demos/myelin.html): a bare and a myelinated axon of the same radius racing, the
  myelinated spike jumping node to node over ten times faster on a tenth of the current; myelin
  thickness, node spacing, nodal channel density.
- A demo, [propagation](demos/propagation.html): a spike travelling along a squid axon at Hodgkin
  and Huxley's speed, its local circuit, the ions that carry the current along the axon (mostly
  K⁺), and each ion's driving force along it; radius (speed as √a), temperature (heat block at
  35 °C), Na⁺ channel block.
- A demo, [action potential](demos/axon.html): a squid axon's spikes, their currents and gates,
  and each ion's driving force through them; the stimulus, temperature, and Na⁺ and K⁺ channel
  block.

## 0.11.4 (2026-10-06)

- Islands behind neighbours that conduct 1e20 times less than they do (across faces that hold μ̄
  level) solve too: there the summed row takes in the outside edge node's balance, so the flux
  is carried by the outside's own first segment rather than by the face flux unknown, which
  elimination could leave to the island's cancellation. Seen through one face of two, the
  island's level had been overshot twofold and swung between the two contacts', and the solves
  kept the plain, wrong level.

## 0.11.3 (2026-10-06)

- A region held only through its faces (a conductance there 1e-24 of its own, or neighbours that
  conduct 1e16 times less) solves to the right level in steady state: its balance summed over
  it, where its own fluxes cancel exactly, goes in as a bordered row once the plain solves
  converge with terminal currents that don't add up. Such solves had converged silently to a
  wrong level, their currents 2× or 16× off; now right to ~1e-15, and directly.

## 0.11.2 (2026-10-06)

- J·v is exact where a solve needs it (the impedance's GMRES and Newton's refined solves): the
  dilute kernels' terms are kept in difference form, each flux on its η difference, instead of a
  central difference of the residual, whose truncation swamped a nearly uniform response. The
  impedance is 15–40% faster, and keeps GMRES's answer always (its fallback to the factorised
  solve, for the old operator's noise, is gone). Newton keeps every converged refinement: a GaAs
  junction without recombination, its ~1.5e-7 A/m² carried by minority carriers from the
  contacts, had cold and warm solves 1e-3 apart, now 1e-14; a dim solar cell's open circuit on a
  0.25 nm grid takes 258 factorisations (790).
- Newton also refines its solves when its undamped updates grow twice running: a bipolar
  stack's floating base (its holes held ~1e14 more weakly than they move within it) now solves
  in steady state, cold and swept, and a MOS capacitor's gate step over 1000 s takes 146
  factorisations (408).
- A current drive whose floated solve fails is found by held solves instead: the voltage that
  passes the target, by regula falsi within the bracket the continuation found, to 1e-10 of the
  current. A closed Fe³⁺/Fe²⁺ cell driven near its limit (Fe³⁺ 20 orders below Fe²⁺, the floated
  system 33 digits short) now solves. Stress: 11,380 of 11,381.

## 0.11.1 (2026-10-05)

- The impedance splits its response into a uniform shift of each region and the rest, where a
  reading is uncertain: a redox electrode at its open circuit against a bath had its DC
  conductance 2–40% low (the current is the slope of levels uniform to 1e-14 of the response);
  now it matches the steady dI/dV to six digits. GMRES's answer is kept where it cuts the
  residual tenfold (a thousandfold had thrown away a real improvement at 0.01 Hz). The stress
  test now passes every impedance at 500 devices per family (11,375 of 11,381 scenarios).

## 0.11.0 (2026-10-05)

Robustness, found by a new stress test of random devices and fixed where it pointed.

**Testing, and what it says**
- `npm run stress`: random but plausible devices in five families (semiconductor stacks, MOS
  capacitors, electrolyte cells, electrodes, liquid junctions), each solved cold, swept warm,
  driven at open circuit or a current, stepped in time and probed by its impedance, every result
  judged by `check()` and by invariants (flat levels at equilibrium, warm and cold agreeing, the
  impedance passive and its DC limit the steady dI/dV). At 500 per family, 11,371 of 11,381
  pass. [docs/reliability.md](docs/reliability.md) has the families, the results and the known
  limits, and the README says plainly what it's tested on, linking there.

**Steady solves**
- A MOS capacitor without a channel port is robust: steady solves pin a species' level flat at
  its contact's where one contact alone reaches it and nothing else touches it, rather than find
  it through the bulk's ~1e3 minority electrons per cm³, and Newton refines its solves by GMRES
  (J·v from the residual) when the factorisation loses a mode, as it does an inversion layer's
  level. A 10 mV gate step reaches the low-frequency charge in ~70 steps instead of 80,000.
- A blocked species mobile through one region is held flat too, at the level its amount fixes:
  a closed Zn | ZnSO₄ | Zn cell past its limiting current, its sulfate excluded to 1e-25 mol/m³,
  solves at any voltage in 6 iterations, where it had failed from 0.7 V.
- The block factorisation perturbs a pivot that cancels to exactly zero (static pivoting): GaAs
  stacks without recombination (minority carriers ~1 per m³) solve, and so does a region held
  only by face conductances 1e17 weaker than its insides. `BlockTridiagonal` takes it as
  `staticPivots` (off by default).
- A terminal driven by a current (open circuit) is held at a voltage and marched to the
  crossing: a solar cell's V_oc in dim light on a 0.25 nm grid, a lit Schottky diode's, a redox
  electrode driven at a current needing 3 V. A current no voltage reaches fails with a warning
  saying so, and what the voltages passed.
- Cold starts: the light's continuation starts at 1e-30 of it and hands a stall to the
  pseudo-transient from its dimmer solution; Newton's divergence check spares the first update
  (a flooded scarce population predicts 1e10 thermal units, which damping covers); the
  pseudo-transient's dt grows gently after a failure and tries the direct solve once past the
  slowest diffusion time; a level start that fails from a layout made for the target voltage is
  laid out again. A lit GaAs phototransistor solves in ~1000 iterations (it failed after
  20,000), a cold start with immobile traps in 191 (280).
- Fixed: a flat level through a face between two regions left the face's flux undetermined.

**Impedance**
- The current is read where its error is least: between the contacts of a device without
  ports, the total current is the same through every cut, so it's read across each segment
  too. A p⁺n⁺ junction without recombination had its capacitance 9× off at 100 Hz (the contact's
  sum of huge terms read a constant conductance of round-off); now its DC limit is the steady
  dI/dV and its capacitance dQ/dV. Where a reading's estimated error passes 1%, the solve is
  continued to GMRES's floor: a redox cell's DC conductance had been 50% off.
- GMRES judges convergence by the true residual (its running estimate had passed unconverged
  solves: n-Si | KCl had a negative resistance at 0.01 Hz), keeps its best iterate, and falls
  back to the factorised solve where it can't converge (a strictly neutral electrolyte above
  1 MHz had gone non-passive).

**Transients**
- A first step whose error doesn't shrink with it (a jump starting a self-similar profile) is
  taken rather than shrunk until Newton fails: transients from a jump take 5–30% fewer
  factorisations. A bath beside a strictly neutral region resolves short steps (a sharp
  3 M | 1 µM junction starts at tol 1e-6). Haynes–Shockley takes 7% fewer factorisations.

**API and checks**
- `set({ ports: { gate: { V: 0.2 } } })` patches ports by name; a new kind of drive drops the
  old, as for contacts.
- `check()` counts a species at rest below 1e-8 of what it or its reactions' (non-metal)
  partners could carry: a supporting ion's round-off flux, or a nearly absent species' reaction
  noise, had read as unsteady.
- Fixed: a `set()` that widened a surfaced window started the new nodes at θ = 0.5, not `theta0`.

**Demos and docs**
- The MOS demo: C–V at any frequency from 10 µHz to 1 MHz against the ideal low- and
  high-frequency curves, the lifetime, the capacitance against frequency, and where the AC
  charge is answered. The demos' charts take `logx`.
- From an outside agent building an OECT from llms.txt: a volumetric capacitance C* as
  `area: 1/t, C: C*·t`; Bernards–Malliaras as the charge-sheet formula without its V_T term;
  `R` and the impedance in Ω with a geometry; which current `trace` records; tightening `tol`
  before fitting a decay's time constant.

## 0.10.0 (2026-10-04)

- The saturation demo has a thin-film transistor below the silver cell: the gate a capacitance
  along the channel, its output characteristic against the charge-sheet model, and the bands
  pinching off at the drain under the gate's flat Fermi level. Each part's controls stay with it
  as you scroll.
- The live demos have a pause button, and pause themselves while scrolled out of view or in a
  hidden tab, so a page doesn't keep a CPU busy unseen.
- Transients in strictly neutral materials start better. The state a first step starts from
  needn't satisfy the algebraic equations (φ, an interface's unknowns), which jump in the first
  instant, and the error estimate had read that jump as error on every later step; the history
  now starts after it. A closed device end is solved in the better-conditioned charge rows like
  an interior node, so an electrode spread through neutral water up to one takes its first step
  at tol 1e-8. `npm run hard` runs a set of hard cases, with a baseline of what still fails.
- Fixed: an electrode surface's coverages weren't under the time step's error control (a Langmuir
  filling ran 4% off its exponential at any tolerance).
- Fixed: species fed only through an electrode surface (A⁺ + e⁻ = S, S + e⁻ = B⁻) looked fed from
  outside, so a steady solve didn't conserve their total; it now does, with what the surface
  holds.
- Fixed: a reaction marked `bare` on a port without surface species ran at zero.
- Fixed: a held port level beside another port's electrode reaction in the same window (the hold
  listed first) left the currents unbalanced.
- Fixed: capacitances on overlapping windows overwrote each other; they add (two gates on a
  channel).
- Fixed: `set()` to a new grid or geometry carried the old conserved amounts (a slab's mol/m² as
  a sphere's mol).
- Fixed: `check()` on a non-planar device misread a face's per-area rates against the totals, and
  says mol/s, A and C there; `describe()` likewise for current drives.
- A capacitance alone driven by a current starts uncharged and charges in a transient; `solve()`
  says it has no steady state instead of running its voltage off. `ports[k].sigma` is zero where a
  contact holds that charge, as the solver counts it. Unit warnings cover ports' k0 and C.
- A floating electrode with a double layer starts uncharged (as the neutral starting composition
  has it) and charges toward its mixed potential; it had started at the mixed potential, charged,
  and the first step failed. A terminal driven by a current where every held terminal passes
  nothing (closed ends) is now an error that says to hold one that does, where Newton diverged.
- device.md's ports section is arranged by what a port does, its outputs in one table; llms.txt
  has templates for corrosion under a drop (an electrode port, a double layer, a cross-section, an
  Evans diagram) and a thin-film transistor.

## 0.9.1 (2026-10-04)

- Fixed: with a capacitance through a port in strictly neutral water, the most abundant ion leaked
  on time steps (about 1e-6 of it, so `check()` failed conservation): the row that stands for
  charge continuity there assumed the ions' net charge constant, but against a double layer it's
  −aσ. It now counts that change.
- `polarization(..., { surface: 'equilibrium' })`: an electrode port's coverages re-equilibrated at
  each potential, the steady-state curve with a passive film's active–passive peak.
- A solution warns when an electrode's surface is covered to a bare fraction below 1e-10, where the
  solver struggles; the docs' film example forms in minutes, and they say a film's μ° is fitted.

## 0.9.0 (2026-10-04)

- A capacitance spread through a port's window, `capacitance: { C, zeroCharge }` with `area`:
  σ = C (V − zeroCharge − φ) per area of electrode, aσ in each node's charge balance (neutrality
  at ε = 0), and the charging current through the port. A gate along a channel (the
  gradual-channel model of a thin-film transistor), an organic electrochemical transistor's
  volumetric capacitance, a porous electrode's double layer. With it an electrode port's potential
  swings continuously when a passive film covers the last active patch. Validated: a TFT against
  the charge-sheet model from below threshold through saturation, the de Levie impedance, the
  gate's low-frequency impedance against dQ/dV (test/capacitance.test.js).
- An electrode port's surface: species on its sites with a coverage at each node (Langmuir,
  μ = μ° + RT ln(θ/θ₀), Γ mol of sites per m² of electrode), named in its reactions like species,
  and `bare: true` for a reaction that runs only on bare metal (its rate times θ₀). Adsorbed
  intermediates, and passive films as a coverage: a film forms where potential and pH favour it,
  and blocks dissolution, the active–passive curve. Solutions report `ports[k].coverage`, and
  `polarization()` holds the coverages with the solution. Validated: the Langmuir isotherm, the
  charge to fill the surface, the blocking (test/coverage.test.js). The bare fraction is carried to
  full precision, and a nearly full surface's storage computed from its small complements. When a
  film passivates an electrode's last active patch the mixed potential must jump; give the port a
  capacitance (its double layer) and the swing is continuous.
- An electrode spread through a port now requires its region strictly neutral (ε = 0): its double
  layers are below the grid, and with ε > 0 its reactions left the solution charged (steady solves
  converged, transients didn't).
- `build()` passes `geometry` through. The docs give the shapes of a port's outputs.

## 0.8.0 (2026-10-03)

- Devices with a cross-section A(x): `geometry: { type: 'spherical' | 'cylindrical', r0 }` (shells
  about a centre r0 to the left of x = 0) or `{ area: { x, values } }`. The finite volumes carry
  totals through A: boxes hold ∫A dx, faces and contacts pass their flux times A there, and each
  segment is weighted by its length over ∫dx/A, so steady diffusion between nodes is exact for
  any A. Currents, fluxes, amounts and charges are then totals (planar devices are unchanged, bit
  for bit, with A = 1 m²). Validated in test/geometry.test.js: steady diffusion to a sphere and a
  cylinder (exact on any grid), Cottrell's transient with the spherical term, uptake by a sphere
  from its surface (centre at x = 0), Debye–Hückel screening around a charged sphere.
- Electrodes spread through a window: a port with `reactions` (Butler–Volmer per area, as at a
  face) against its terminal species, a metal's carrier at the port's level, and an `area` of
  electrode per volume (a number or a profile). Used for the metal under a thin film of electrolyte,
  a crevice's walls, or a porous electrode with a well-conducting matrix. Held, it's a
  potentiostat; floating (or with nothing else to carry current) it settles at the mixed potential,
  each spot a net anode or cathode. Solutions report each reaction's rate along the window
  (`ports[k].x`, `.rates`, `.area`). Validated against the transmission line (linear kinetics along
  a bar) and the Wagner–Traud mixed potential of two couples. A floating electrode port starts at
  its mixed potential.
- `polarization(device, sol, where, V)` in the kit: an electrode's reactions (at a face, or at a
  spot of an electrode port) with the solution beside it held and the metal's level moved. It
  returns each reaction's rate, partial current (anodic positive) and redox level, and the net
  current: the curves of an Evans diagram. It uses the solver's own law, so at the metal's actual
  level it gives the solved rates (test/polarization.test.js, with the Tafel slopes).
- Fixed: a port that only exchanges (an O₂ supply from the air, say) in a strictly neutral region
  could stop a transient at its first step. With a face reaction setting the electrolyte's
  potential, as in a drop of salt water on iron, no step converged. The port's window is now solved
  in the same well-conditioned terms as the rest of the region. A port that exchanges only neutral
  species carries no current, and driving one by a current is now an error.
- The liquid-junction demo says what the minority ions show: Mg²⁺'s plateau, and K⁺ riding the
  front.

## 0.7.0 (2026-10-03)

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
- Liquid-junction potentials in the kit: `henderson()` (the closed form JPCalc uses) and `planck()`
  (the steady state of constrained diffusion, by shooting on the Nernst–Planck equations), for any
  mixture and valences. Validated in test/junction.test.js against JPCalc's published values,
  against driftlet's own steady junction (to 2e-4 mV) and the stationary codes LJPcalc and JLJP
  (to 0.02 mV: a patch pipette's K-gluconate, ZnCl₂ | KCl), and a free-diffusion junction grown
  from a sharp boundary, constant in time and between the two.
- A liquid-junction demo: a patch pipette's junction growing from contact against Henderson,
  Planck and the published values (and why patch-clampers correct for it), ZnCl₂ | KCl, and a
  3 M KCl salt bridge whose charged frit turns permselective against a dilute sample.
- A charged-nanochannel demo: a channel 10 nm across with charged walls, resolved along its axis;
  Donnan layers at each mouth that overlap in dilute salt (where Poisson's equation, not
  electroneutrality, sets what it does), K⁺ selectivity, the conductance plateau at low salt
  against Teorell–Meyer–Sievers, the zero-current voltage across a salt gradient, and
  rectification by a one-sided or bipolar charge, explained as the pn junction it is. Validated in
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
