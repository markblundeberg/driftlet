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
  a metal region). `bath(c, reference, drive, { offset })` is a bath (`offset` places the reference species
  relative to the terminal voltage, e.g. to read a solution's potentials against SHE). A drive is a voltage, a waveform, or
  `{ V }`, `{ I }`, `{ V, R }`.
- **`library`** takes pieces with `species` and `materials`, merged with any given directly
  (`combine()` does the merging on its own). A species or material given twice must be given
  identically.
- Nested lists in the stack are flattened, so a helper can return several items.
- Reactions pass through as written (see below).

## Reactions and half-reactions

Reactions are written as equations, `'Ag+ + e- = Ag(s)'`, in the definition itself (see the
[device reference](device.md#interfaces)), so `build()`, `new Device` and `set()` all take them.
A coefficient is an integer and a space (`2 e-`); names may start with digits (`3He`).

`half(equation, fixed)` writes a half-reaction, with electrons `e-` on one side, as plain
`{ equation, fixed }`: spread it into a face reaction with its kinetics. `level(sol, half)` is the
electronic level it implies at each node, as an electron voltage: where an electrode exchanging
electrons by that couple would sit in equilibrium with the local composition, the redox level
of ESBD diagrams (`NaN` where a participant is absent). With `{ standard: true }` every species
is at its reference concentration instead. `SHE` is 2 H⁺ + 2 e⁻ ⇌ H₂ with μ(H₂) = 0, so its
standard level is φ + μ°_H⁺/F.

```js
import { Device } from 'driftlet';
import { build, layer, ohmic, half, level } from 'driftlet/kit';

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
const sol = new Device(def).solve();
const redox = level(sol, silver); // V_Ag⁺(x) in the solution
console.log(`redox level at mid-cell: ${redox[sol.x.length >> 1].toFixed(4)} V`);
```

## Plotting

Solutions are plain arrays, ready for any plotting toolkit. Two helpers save the bookkeeping of
a level diagram:

- `traces(sol, { species, standard, phi, levels, labels, shifts })` (in `driftlet/kit`) returns
  its lines as data: per charged species its voltage V_i (`kind: 'level'`) and standard level
  V°_i (`'standard'`, the band edges for e⁻ and h⁺), optionally φ and half-reaction levels
  (`levels: [{ half, label, standard }]`, kinds `'redox'` and `'redox-standard'`). `labels`
  renames lines by id (`{ 'V:e-': 'Fermi level' }`). `shifts` gives species display offsets in
  volts (`{ 'K+': 3.5 }`), moving all of a species' lines together so that widely separated
  species share one readable plot, as the ESBD book does; the label says so (⌇). Each line has a
  colour `role` (`'electron'` for e⁻, `'cation'` (holes too), `'anion'`, `'redox'`, `'phi'`) and a
  `slot` within it, so a colour follows its species whichever lines are shown
  (`speciesRole(sol, name)` gives a species' own). Undefined values are `NaN` (break the line
  there), and doubled interface nodes share an x, so steps draw as vertical lines. It also gives
  the regions (`{ name, material, x0, x1 }`), the faces' positions and a suggested range.
- `bandDiagram(sol, opts)` in `driftlet/plot` draws them as an SVG string, with no DOM needed, in
  the ESBD book's style: species voltages thick and solid, standard levels thin and solid, redox
  levels thick and dashed, standard redox levels thin and dashed, φ thin, dotted and grey;
  electrons steel blue, cations warm colours, anions cool ones, redox levels blue-violets (a
  palette checked for colour-vision deficiency, light and dark); a ⌇ across a shifted species,
  with a note. The colours are CSS custom properties (`--driftlet-electron`,
  `--driftlet-cation-1` …, `--driftlet-anion-1` …, `--driftlet-redox-1` …, `--driftlet-ink`, …)
  with defaults of no specificity, so a rule such as `.figure svg { --driftlet-electron: … }`
  restyles them; `THEME` and `themeColor(role, slot, { dark })` give them to other charts.
  `xlim` (m) zooms into a window, with ticks and range fitted to it. `levelChart(traces, opts)`
  draws traces you've edited.

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
  (cm²/s?), a nonzero μ° under 100 J/mol (eV or volts?), a capacitance over 10 F/m²
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
