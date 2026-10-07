# The kit

`driftlet/kit` is the porcelain: helpers that write the plain device definition for you, so a
device can be set up the way it's drawn and from familiar quantities. It never makes a physical
choice. Alignments, offsets, rate constants and every other physical input stay yours to give,
and what the kit produces is an ordinary definition you can print, edit and pass to
`new Device` (see the [device reference](device.md)).

Its data library, with sources, is described in [data](data.md).

## A device as a stack

`build()` takes the device left to right: a contact, then layers and the faces between them,
then a contact.

```js
import { Device, units } from 'driftlet';
import { build, layer, ohmic } from 'driftlet/kit';

const Si = {
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: {
      epsr: 11.7,
      species: {
        'e-': { D: 36e-4, mu0: 0, cRef: units.perCm3(2.8e19) },
        'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: units.perCm3(1.04e19) },
      },
    },
  },
};

const def = build({
  library: [Si],
  stack: [
    ohmic(0),
    layer('Si', units.um(1), { donors: units.perCm3(1e17) }),
    layer('Si', units.um(1), { acceptors: units.perCm3(1e16) }),
    ohmic(0.5),
  ],
  bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
  grid: { hmin: 1e-9, hmax: 20e-9 },
});
console.log(JSON.stringify(def.regions)); // plain data: fixedCharge, in C/m³
const sol = new Device(def).solve();
console.log(`I = ${sol.current.toFixed(1)} A/m²`);
```

- **Layers** are regions: anything with a `material`. `layer(material, length, more)` writes
  one. Doping can be given as `donors` and `acceptors` (mol/m³) in place of `fixedCharge`.
- **Faces** are anything else between layers: an interface's fields, exactly as in
  `interfaces[f]`. Two adjacent layers with no face between them get the default face, which
  between different materials is an error asking for the alignment, as always.
- **Contacts** are the first and last items. `ohmic(drive, carriers)` holds each carrier in
  equilibrium with a metal at the terminal voltage (default `['e-', 'h+']`, an
  infinite-recombination contact; `['e-']` alone is selective, and is the current collector on
  a metal region). `bath(c, reference, drive, { offset })` is a bath whose terminal voltage is
  the level of its `reference` species, as an electrode reversible to it would read (`offset`
  places that species relative to the terminal voltage, e.g. to read a solution's potentials
  against SHE). `bath(c, drive)`, without a reference, is a bath whose terminal voltage is its
  φ: an ideal salt bridge, the convention of membrane potentials, so that between two such baths
  V_right − V_left is the φ difference a voltage clamp sets. A drive is a voltage, a waveform, or
  `{ V }`, `{ I }`, `{ V, R }`.
- **`library`** takes pieces with `species` and `materials`, merged with any given directly
  (`combine()` does the merging on its own). A species or material given twice must be given
  identically.
- Nested lists in the stack are flattened, so a helper can return several items.
- Reactions pass through as written (see below), and so do `T`, `grid`, `ports` and `geometry`.

## Reactions and half-reactions

Reactions are written as equations, `'Ag+ + e- = Ag(s)'`, in the definition itself (see the
[device reference](device.md#interfaces)), so `build()`, `new Device` and `set()` all take them.
A coefficient is an integer and a space (`2 e-`); names may start with digits (`3He`).

`half(equation, fixed)` writes a half-reaction, with electrons `e-` on one side, as plain
`{ equation, fixed }`: spread it into a face reaction with its kinetics. `level(sol, half)` is the
electronic level it implies at each node, as an electron voltage: where an electrode exchanging
electrons by that couple would sit in equilibrium with the local composition: the couple's redox
level (`NaN` where a participant is absent; see [reading level diagrams](visualization.md)). With `{ standard: true }` every species
is at its reference concentration instead. `SHE` is 2 H⁺ + 2 e⁻ ⇌ H₂ with
$`\mu_{\mathrm{H_2}} = 0`$, so its standard level is $`\phi + \mu^\circ_{\mathrm{H}^+}/F`$.

```js
import { Device } from 'driftlet';
import { build, layer, ohmic, half, level, polarization } from 'driftlet/kit';

const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 });
const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 }; // an assumed pzc
const def = build({
  species: [
    { name: 'Ag+', z: 1, cRef: 1000 },
    { name: 'NO3-', z: -1, cRef: 1000 },
    { name: 'e-', z: -1 },
  ],
  materials: {
    water: { epsr: 78.5, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } },
    Ag: { conductor: { species: 'e-', conductivity: 6.3e7 } },
  },
  stack: [
    ohmic(0, ['e-']),
    layer('Ag', 1e-6),
    { ...stern, reactions: [{ ...silver, k0: 1e-3, alpha: 0.5 }] },
    layer('water', 10e-6, { c0: { 'Ag+': 10, 'NO3-': 10 } }),
    { ...stern, reactions: [{ ...silver, k0: 1e-3, alpha: 0.5 }] },
    layer('Ag', 1e-6),
    ohmic(0.05, ['e-']),
  ],
  grid: { hmin: 0.1e-9, hmax: 100e-9, ratio: 1.15 },
});
console.log(JSON.stringify(def.interfaces[0].reactions[0])); // { left: { e-: -1, Ag(s): 1 }, right: { Ag+: -1 }, … }
const dev = new Device(def), sol = dev.solve();
const redox = level(sol, silver); // V_Ag⁺(x) in the solution
console.log(`redox level at mid-cell: ${redox[sol.x.length >> 1].toFixed(4)} V`);
// The left electrode's polarization curve, 0.2 V either side of its own level (see below).
const metalSide = sol.V['e-'][0];
const curve = polarization(dev, sol, { face: 0 }, Array.from({ length: 81 }, (_, j) => metalSide - 0.2 + j * 0.005));
console.log(curve.reactions[0].level, curve.current[40]); // the couple's level here, and the current at the metal's own
```

### Polarization curves: an Evans diagram

`polarization(device, sol, where, V)` evaluates every electrode reaction at one spot with the
solution beside it held as `sol` has it and the metal's level set to each `V` (an electron
voltage, the terminals' convention). It works at a face (`{ face: f }`, the face's metal side) or
at a spot of an [electrode port](device.md#electrodes-spread-through-a-window)
(`{ port: name, x }`, with its surface's coverages there held too). Each reaction comes back with its rate (mol/(m²·s), forward), its partial
current (A/m², anodic positive: into the solution) and its `level`, where the rate vanishes (the
couple's redox level there, as `level()` gives it). The net current is their sum. Plotting
log|current| against V gives the Evans diagram. The net crosses zero at the spot's open-circuit
potential, between the couples' levels, and the metal's actual level is where the solve put it.
The law is the solver's own, so at that level the rates are the solution's. Against SHE, subtract
the SHE level (`sol.phi` with the usual tables' ions).

The example above ends with the left electrode's curve.

An electrode port's surface is held too by default. With `{ surface: 'equilibrium' }` (a fifth
argument) its coverages are re-equilibrated at each V instead, each species by the reaction that
makes it from the solution: the steady-state curve, where a passive film's active–passive peak
shows (`curve.bare` gives the bare fraction along it).

Concentrations are held, so these are the curves of that instant (for a mixed potential, the
kinetic picture). Away from the actual level the composition beside the electrode would change,
which is the transport the full solve accounts for.

## Channels

**`hodgkinHuxley({ inside, outside, T, g, leak, leakZ, at, names, linear, area })`**: Hodgkin and Huxley's
squid-axon channels for a [membrane face](device.md#voltage-gated-channels) with the outside on its
left: `{ gates, species }`, the face's `gates` (m, h and n, their 1952 rates at 6.3 °C scaled to
`T` by Q₁₀ = 3) and links for Na⁺ (P·m³h), K⁺ (P·n⁴) and the leak (one ion, `leak`, of charge
`leakZ`: Cl⁻ by default). The currents are GHK's, which rectify, so each permeability is the one
whose chord conductance at `at` (−65 mV) is HH's `g` (120, 36 and 0.3 mS/cm²) at the
concentrations given. Spread it into the face beside its $`\phi`$ law, adding the other species'
links: `{ phi: { type: 'capacitive', C: 0.01 }, gates: hh.gates, species: { ...hh.species, 'A-': 'blocked' } }`.
With `linear: true` (and no concentrations) the face gets HH's own linear channels instead:
conductance links, `G` = g, gated the same, whose currents are g(V − E) with E each ion's Nernst
level at the concentrations either side; between two `bath`s without a reference species, that's
their voltage clamp. For a [membrane port](device.md#a-membrane-through-a-window-gated-channels) along a region (an axon
along x), give `area`, the membrane per volume (2/a for radius a), in place of `inside`, `outside`
and `at`: the links are then their linear conductances, `G` = g × area, gated the same, each ion's
reversal potential following from the port's `bath` and the concentrations inside.
`ghkCurrent(z, inside, outside, V, T)` is the GHK current (outward, A/m²) at unit permeability.

## Sources

Helpers that write sources as plain data:

- **Waveforms** for a terminal's `V` or `I`: `pulse({ start, width, amplitude, rise, base })`,
  `square({ high, low, period, duty })`, `triangle({ from, to, period })` (a voltammogram's sweep)
  and `ramp({ from, to, duration, start })`. Each returns `{ t, values }` (with `repeat` for the
  periodic ones), times from the device's time 0. A pulse keeps its area,
  `amplitude · width`, whatever its `rise`; with no rise its edges are jumps, which the time
  stepping handles exactly, so it must start after t = 0 (where a device starts, and `solve()`
  reads it).
- **`injector({ name, region, from, to, species, I, headroom })`**: a port that injects a species
  into a window at a driven current (a number or a waveform), spread evenly. Its conductance is
  set small enough that at the peak current the port sits `headroom` volts (default 2.5) above
  the window, so every node takes the same share, like a current source.
- **`recombination({ material, tau, majority })`**: band-to-band recombination e⁻ + h⁺ ⇌ 0 as
  a bulk reaction, from the low-injection minority lifetime it gives at that majority density,
  $`k_f = 1/(\tau \cdot n_{\mathrm{maj}})`$.
- **`photogeneration({ material, flux, alpha, mu, from, to, makes })`**: light absorbed as it goes in
  (Beer–Lambert). A photon flux (mol/(m²·s); a power P at photon energy E is
  $`P/(E N_A)`$) enters the material at x = `from` heading toward +x and generates pairs at
  $`G(x) = \Phi\alpha e^{-\alpha(x - \mathrm{from})}`$ until `to`, the material's far end. It's
  photon ⇌ e⁻ + h⁺ from a photon reservoir at $`\mu`$ = `mu` (J/mol, well above the gap, e.g.
  `units.eV(3)`), with `kf` a profile against x. `makes` replaces the right-hand side: in an
  organic semiconductor, light makes neutral excitons, `makes: 'X'`. The table generates exactly what the light
  loses, $`\Phi(1 - e^{-\alpha(\mathrm{to} - \mathrm{from})})`$, and the solver averages it
  over each node's box, so the total stays right on a grid coarser than $`1/\alpha`$ (where
  carriers start within a cell still needs the grid). Light through several materials is one
  call per material, each with the flux that reaches it.

Read a transient at points inside the device with `advance(t, { probes })` (see the
[device reference](device.md#using-a-device)). And `units` converts back out of SI for display:
`toPerCm3`, `toMolar`, `toCm2PerS`, `toUm`, `toNm`, `toEV`.

```js
import { Device, FARADAY, units } from 'driftlet';
import { build, layer, ohmic, semiconductor, pulse, injector, recombination } from 'driftlet/kit';

// Holes injected for 0.5 µs into n-Ge under 10 V/cm (Haynes–Shockley), read 1 mm downstream.
const ND = units.perCm3(1e15), I = (0.01 * ND * FARADAY * 50e-6) / 0.5e-6; // 1% of n₀ over 50 µm
const dev = new Device(
  build({
    T: 300,
    library: [semiconductor('Ge')],
    stack: [ohmic(3), layer('Ge', 3e-3, { name: 'bar', donors: ND }), ohmic(0)],
    bulkReactions: [recombination({ material: 'Ge', tau: 20e-6, majority: ND })],
    ports: [injector({ region: 'bar', from: 0.475e-3, to: 0.525e-3, species: 'h+', I: pulse({ width: 0.5e-6, amplitude: I, rise: 1e-9 }) })],
    grid: { hmin: 10e-6, hmax: 10e-6 },
  }),
);
dev.solve();
const run = dev.advance(8e-6, { dtMax: 20e-9, probes: [{ x: 1.5e-3, species: 'h+' }] });
const seen = run.trace.probes[0], k = seen.indexOf(Math.max(...seen));
console.log(`the pulse passes at ${(run.trace.t[k] * 1e6).toFixed(2)} µs, Δp ≈ ${units.toPerCm3(seen[k]).toExponential(1)} cm⁻³`);
```

## Liquid junctions

`henderson(left, right, ions, { T })` and `planck(left, right, ions, { T })` give the junction
potential $`\phi_R - \phi_L`$ (V) between two neutral solutions (mol/m³ by ion), with `ions` their
`{ z, D }` (the data library's `IONS` has that shape):

- **Henderson's** formula takes the junction to be the two solutions mixed in one proportion at
  every point (so every ion's mixing fraction follows the same curve). It's
  a closed form, exact for one salt, and what JPCalc and most electrophysiology corrections use.
- **Planck's** is the steady state of diffusion through a zone held between the two solutions,
  electroneutral and carrying no current, solved by shooting on the Nernst–Planck equations. It's
  exact for that junction, and equal to driftlet's own steady state between two baths.

A junction growing freely from first contact is neither: simulate it, as a transient from a
`c0` step. For one salt all three agree; for mixtures they differ by up to a millivolt or so.
Validated against JPCalc (Henderson), and LJPcalc and JLJP (stationary Nernst–Planck), in
test/junction.test.js.

```js
import { henderson, planck, IONS } from 'driftlet/kit';

// 50 mM NaCl against 50 mM KCl, and a K-gluconate pipette (gluconate's D from its relative mobility, 0.33 of K⁺'s).
const ions = { ...IONS, 'gluconate-': { z: -1, D: 0.33 * IONS['K+'].D } };
console.log(1000 * henderson({ 'Na+': 50, 'Cl-': 50 }, { 'K+': 50, 'Cl-': 50 }, ions), 'mV');
const pipette = { 'K+': 140, 'gluconate-': 130, 'Cl-': 10 }, bath = { 'Na+': 145, 'K+': 5, 'Cl-': 150 };
console.log(`Henderson ${(1000 * henderson(pipette, bath, ions)).toFixed(2)} mV, Planck ${(1000 * planck(pipette, bath, ions)).toFixed(2)} mV`);
```

## Recording a transient

`recorder(device, { every, times, probes })` runs a transient in pieces and keeps what a demo
needs to show it: `trace`, the terminal current and voltage (and the probes' readings) at every
accepted step from the start, joined up across calls; and `frames`, whole solutions to scrub back
through or draw side by side. Frames are kept at the start, then at each of `times` or every
`every` seconds, which the stepping lands on exactly; with neither, at the end of each call.
Its `advance(tEnd, opts)` takes `Device.advance()`'s options, so an animation calls it once per
frame with a `budgetMs` (`done` says whether it got there), and a run in budgets is step for
step the run made in one go. `frame(t)` is the frame at or just before t.

```js
import { Device } from 'driftlet';
import { build, layer, bath, recorder } from 'driftlet/kit';

// A salt step relaxing: 100 mM KCl in from a bath, into 10 mM held behind a blocking wall.
const dev = new Device(
  build({
    species: [
      { name: 'K+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 0, species: { 'K+': { D: 1.96e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
    stack: [bath({ 'K+': 100, 'Cl-': 100 }, 'Cl-'), layer('water', 100e-6, { c0: { 'K+': 10, 'Cl-': 10 } }), { species: { 'K+': 'blocked', 'Cl-': 'blocked' }, phi: 'neutral' }],
    grid: { hmin: 2e-6, hmax: 2e-6 },
  }),
);
const rec = recorder(dev, { every: 1, probes: [{ x: 100e-6, species: 'K+' }] });
rec.advance(10, { tol: 1e-4 });
console.log(`${rec.frames.length} frames, ${rec.trace.t.length} steps; at the wall after 5 s: ${rec.frame(5).c['K+'].at(-1).toFixed(2)} mol/m³`);
```

## Curated devices

Whole devices from a published parameter set, with reference results from an independent code
to check a run against, so a demo starts from numbers that someone has measured or
benchmarked rather than invented ones. Each returns an ordinary definition (print it to see
every parameter, or edit it).

- **`perovskiteCell({ V, light, bulkSRH })`**: IonMonger's default planar perovskite solar
  cell (Courtier et al., J. Comput. Electron. 18, 1435 (2019)): TiO₂ (100 nm) | MAPbI₃ (400 nm)
  with mobile iodide vacancies over an equal immobile background | spiro-OMeTAD (200 nm), at
  298 K, lit through the TiO₂ (1.4e21 photons/(m²·s) at α = 1.3e7 /m, times `light`). SRH at
  both faces and, unless `bulkSRH: false`, in the perovskite. `V` drives the spiro side's contact,
  so photocurrent is positive. `PEROVSKITE_SCANS` holds IonMonger's own J–V scans of it (steady
  at 1.2 V under light, then 1.2 → 0 → 1.2 V, at seven rates from 1 mV/s to 1 kV/s): per scan
  the current at 0 V, each sweep's maximum power and V_oc, and the hysteresis index, with bulk
  SRH (`full`) and without (`noBulkSRH`). driftlet's test suite reproduces them (and the whole
  loops, to 0.03 mA/cm² where they're gentle and 1 mV where steep). `hysteresis(trace, turn)`
  reads the same figures from a run's trace.

```js
import { Device } from 'driftlet';
import { perovskiteCell, PEROVSKITE_SCANS, hysteresis, recorder } from 'driftlet/kit';

const rate = 10, half = 1.2 / rate; // V/s
const dev = new Device(perovskiteCell({ V: 1.2 }));
dev.solve();
dev.set({ contacts: { right: { V: { t: [0, half, 2 * half], values: [1.2, 0, 1.2] } } } });
const rec = recorder(dev);
rec.advance(2 * half, { tol: 1e-4, dtMax: (2 * half) / 480 });
const ours = hysteresis(rec.trace, half), theirs = PEROVSKITE_SCANS.full.find((s) => s.rate === rate);
console.log(`P_max ${ours.rev.Pmax.toFixed(1)} W/m² going down (IonMonger ${theirs.rev.Pmax}), hysteresis index ${ours.hi.toFixed(4)} (${theirs.hi})`);
```

## Plotting

Solutions are plain arrays, ready for any plotting toolkit. Two helpers save the bookkeeping of
a level diagram:

- `traces(sol, { species, standard, phi, levels, labels, shifts })` (in `driftlet/kit`) returns
  its lines as data: per charged species its voltage $`V_i`$ (`kind: 'level'`) and standard level
  $`V^\circ_i`$ (`'standard'`, the band edges for e⁻ and h⁺), optionally $`\phi`$ (`phi: true`, or a
  list of regions to draw it in only, by name or material, such as an insulator's:
  `phi: ['SiO2']`; see [what's left out](visualization.md#whats-left-out-and-why)) and
  half-reaction levels (`levels: [{ half, label, standard }]`, kinds `'redox'` and
  `'redox-standard'`). Labels typeset the species as a subscript (`V_{SO₄²⁻}`, the species as
  `typeset(name, z)` writes it); `labels` renames lines by id (`{ 'V:e-': 'Fermi level' }`).
  `shifts` gives species display offsets in volts (`{ 'K+': 3.5 }`), moving all of a species'
  lines together so that widely separated species share one readable plot; the label says so
  (⌇). Each line has a colour `role` (`'electron'` for e⁻, `'cation'` (holes too), `'anion'`,
  `'redox'`, `'phi'`) and a `slot` within it, so a colour follows its species whichever lines are
  shown (`speciesRole(sol, name)` gives a species' own). Undefined values are `NaN` (break the line
  there), and doubled interface nodes share an x, so steps draw as vertical lines. It also gives
  the regions (`{ name, material, x0, x1 }`), the faces' positions and a suggested range.
- `energy: true` (in `traces` and `bandDiagram`) gives the familiar band diagram instead, energy
  up: electrons and holes as electron energies in eV (E_c and E_v from the standard levels, E_Fn
  and E_Fp, or E_F where there are no holes, from the levels), with redox levels and $`\phi`$ (as
  $`-q\phi`$) turned over too and shifts in eV. Ions have no band, so asking for one is an error.
  It's the right picture for semiconductor devices; species voltages earn their place where ions
  share the diagram (see [which to use](visualization.md)).
- `bandDiagram(sol, opts)` in `driftlet/plot` draws them as an SVG string, with no DOM needed, in
  the style [reading level diagrams](visualization.md) explains: species voltages thick and solid, standard levels thin and solid, redox
  levels thick and dashed, standard redox levels thin and dashed, $`\phi`$ thin, dotted and grey;
  electrons steel blue, cations warm colours, anions cool ones, redox levels blue-violets (a
  palette checked for colour-vision deficiency, light and dark); a ⌇ across a shifted species,
  with a note. The colours are CSS custom properties (`--driftlet-electron`,
  `--driftlet-cation-1` …, `--driftlet-anion-1` …, `--driftlet-redox-1` …, `--driftlet-ink`, …)
  with defaults of no specificity, so a rule such as `.figure svg { --driftlet-electron: … }`
  restyles them; `THEME` and `themeColor(role, slot, { dark })` give them to other charts.
  `xlim` (m) zooms into a window, with ticks and range fitted to it; `ytick(v)` writes the y
  axis's tick labels. `levelChart(traces, opts)` draws traces you've edited, or any profiles
  against x ($`\log_{10}`$ concentrations, say, with `ytick` writing powers of ten). Labels, titles and
  the y label mark subscripts TeX-style, `V_{e⁻}` or `C_ox` (a run of letters, digits and
  charge signs), drawn lowered and smaller; tooltips and the accessible name read `V_e⁻`.
  `labelParts(label)` splits a label into its plain and subscript runs, for other renderers.

```js
import { Device } from 'driftlet';
import { build, layer, bath } from 'driftlet/kit';
import { bandDiagram } from 'driftlet/plot';

const def = build({
  species: [
    { name: 'Na+', z: 1, cRef: 1000 },
    { name: 'Cl-', z: -1, cRef: 1000 },
  ],
  materials: { water: { epsr: 78.5, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } } },
  stack: [bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-'), layer('water', 10e-6), bath({ 'Na+': 10, 'Cl-': 10 }, 'Cl-', { I: 0 })],
  grid: { hmin: 0.2e-9, hmax: 50e-9, ratio: 1.15 },
});
const svg = bandDiagram(new Device(def).solve(), { title: 'liquid junction', phi: true, xlim: [0, 1e-6] });
console.log(svg.length); // in a page: element.innerHTML = svg
```

## Checking a definition

`describe(def)` validates a definition and returns a readable summary: each region with its
doping, Debye length and the grid at its ends, each face's law, alignment, blocked species and
reactions, the contacts, bulk reactions, and characteristic times (dielectric relaxation,
diffusion across each region). It ends with warnings:

- numbers that look like unit slips (`unitWarnings(def)` gives these alone): a D beyond 1 m²/s
  (cm²/s?), a nonzero $`\mu^\circ`$ under 100 J/mol (eV or volts?), a capacitance over 10 F/m²
  (µF/cm²?), a concentration over 1000 M, lengths over a metre or under an atom, a T under
  200 K (°C?), and so on;
- double layers the grid won't resolve: end cells coarser than the Debye length (the solution's
  `warnings` check this again against the solved concentrations).

They're heuristics, prompts to check, never errors.

```js
import { units } from 'driftlet';
import { build, layer, ohmic, describe } from 'driftlet/kit';

const def = build({
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: { epsr: 11.7, species: { 'e-': { D: 36, mu0: 0, cRef: units.perCm3(2.8e19) }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: units.perCm3(1.04e19) } } },
  },
  stack: [ohmic(0), layer('Si', 1e-6, { name: 'n', donors: units.perCm3(1e17) }), layer('Si', 1e-6, { name: 'p', acceptors: units.perCm3(1e16) }), ohmic(0)],
  grid: { hmin: 1e-9, hmax: 20e-9 },
});
console.log(describe(def)); // … warnings: materials.Si.species.e-.D: 36 m²/s is beyond any real diffusivity …
```

## Checking a result

`check(device, sol, { refine, tol })` checks a solved device, and says what it compared, so a
page can show that its numbers hold up. Pass the solution `solve()` or `advance()` returned.

- **converged**: the solve converged.
- **warnings**: the solution's warnings (unresolved double layers, steep profiles, a weakly held
  part) and the definition's likely unit slips. They're heuristics, so they're marked `?`, a
  prompt to look, and don't fail the check.
- **balance**, in a steady state: each species' ledger, what comes in through each terminal and
  what each reaction (bulk, or at a face) makes or uses, in mol/(m²·s) and for a charged species
  as a current. The terms must sum to zero. A species that a face lets through only by a law (a
  membrane's permeability, a conductance) or by a reaction has a ledger on each side, with what
  crosses as a term: in a cell, Na⁺ leaks into the cytoplasm as fast as the pump puts it out. Under light it's J = F(G − R), and it says where
  every carrier went: collected, or recombined in the bulk or at which face. A species at rest
  (in equilibrium, flows only of round-off against what it could carry) isn't listed.
- **conservation**, in a transient: every closed, unreacting stretch keeps what it held plus
  what came in.
- **grid**, in a steady state (unless `refine: false`): the device solved again on a grid twice
  as fine everywhere (half of each `hmin` and `hmax`, the square root of each `ratio`), starting
  from this solution interpolated onto it. The current, each floating terminal's voltage (an open-circuit voltage) and each
  region's charge must agree to `tol` (relative; default 1e-2, about what a plot shows; a
  benchmark wants 1e-3 or less). The discretisation is second order, so the change estimates
  this grid's own error (about ¾ of it). It costs one more solve, so run it once, not on every
  slider move. In a transient it's `?`, not checked: the answer there depends on the whole run
  (and on drives set along the way), so run it again yourself on a grid twice as fine.

It returns `{ ok, items, text }`: `ok` when nothing failed, an item per check
(`{ name, ok, summary, details }`, with `ok: null` for a prompt to look), and `text`, a line
per check marked `ok`, `FAIL` or `?`, ready to log or put on the page.

```js
import { Device, units } from 'driftlet';
import { build, layer, ohmic, semiconductor, photogeneration, check } from 'driftlet/kit';

// An n⁺p silicon cell under blue-green light (α = 1e5 /m), at 0.5 V.
const W = units.um(60.5);
const dev = new Device(
  build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e19) }), layer('Si', units.um(60), { acceptors: units.perCm3(1e16) }), ohmic(0.5)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }, photogeneration({ material: 'Si', flux: 3e-3, alpha: 1e5, mu: units.eV(3), to: W })],
    grid: { hmin: units.nm(1), hmax: units.um(1) },
  }),
);
const report = check(dev, dev.solve());
console.log(report.text);
// ok   balance: 2 ledgers sum to zero …
//        e-: photon = e- + h+: +2.99e-3 (289 A/m²); left contact: -2.23e-3 (-215 A/m²); …
// FAIL grid: on a grid twice as fine (172 → 339 nodes), the largest change is region 0 (Si)
//      charge, … (0.021, over 0.01: refine the grid …): the current is fine (1e-4), but the
//      emitter's depletion charge is 2 % off on this grid; hmin: units.nm(0.25) fixes it.
```

## Live demos

`live(def, opts)` wraps a device for sliders. Call `set(patch)` as often as a control moves
(the patch merges into the definition, as `Device.set()` does), and draw what comes back:

```js nocheck
const dev = live(def, { worker: true, onsolution: (sol, info) => (figure.innerHTML = bandDiagram(sol)) });
slider.oninput = () => dev.set({ contacts: { right: { V: +slider.value } } });
```

- Solves run one at a time, each warm-started from the last. Changes made while one runs are
  merged and solved together, so a fast slider never queues up stale work. On the page's own
  thread, each solve is its own task, so the page paints and takes input in between.
- A change that fails from the warm start is ramped to from the last good state, in more and
  more steps (every number that differs is interpolated), then tried from scratch. If all of
  that fails, the last good solution stays, with `info.failed` and the solver's `warnings`.
  An invalid change does the same, with `info.error`. `info.ramp` says how many ramp steps were
  needed, `info.ms` how long it took.
- `set()` returns a promise of `{ solution, info }`, and `ready` is the first solve's.
  `solution` is always the last good one.
- `worker: true` solves in a Web Worker (a module worker), keeping the page responsive however
  long a solve takes. It works with the library loaded from a CDN too.

```js
import { units } from 'driftlet';
import { build, layer, ohmic, live } from 'driftlet/kit';

const dev = live(
  build({
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: {
      Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: units.perCm3(2.8e19) }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: units.perCm3(1.04e19) } } },
    },
    stack: [ohmic(0), layer('Si', 1e-6, { donors: units.perCm3(1e17) }), layer('Si', 1e-6, { acceptors: units.perCm3(1e16) }), ohmic(0)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    grid: { hmin: 1e-9, hmax: 20e-9 },
  }),
);
for (const V of [0.1, 0.2, 0.3, 0.4]) dev.set({ contacts: { right: { V } } }); // as a slider would
const { solution, info } = await dev.set({ contacts: { right: { V: 0.5 } } });
console.log(`I(0.5 V) = ${solution.current.toFixed(1)} A/m², in ${info.ms.toFixed(1)} ms`);
```

## Alignment from vacuum levels

`vacuumLevel`, `vacuumDipole` and `vacuumZeroCharge` turn vacuum-level estimates (electron
affinities, work functions, Trasatti's absolute SHE) into an interface's alignment. See the
[alignment guide](alignment.md).
