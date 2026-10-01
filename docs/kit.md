# The kit

`driftlet/kit` is the porcelain: helpers that write the plain device definition for you, so a
device can be set up the way it's drawn and from familiar quantities. It never makes a physical
choice. Alignments, offsets, rate constants and every other physical input stay yours to give,
and what the kit produces is an ordinary definition you can print, edit and pass to
`new Device` (see the [device reference](device.md)).

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
  bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e13 } }],
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
  a metal region). `bath(c, reference, drive)` is a bath. A drive is a voltage, a waveform, or
  `{ V }`, `{ I }`, `{ V, R }`.
- **`library`** takes pieces with `species` and `materials`, merged with any given directly
  (`combine()` does the merging on its own). A species or material given twice must be given
  identically.
- Nested lists in the stack are flattened, so a helper can return several items.

## Reactions as equations

Anywhere a reaction goes, `equation` can stand in for its stoichiometry: `'Ag+ + e- = Ag(s)'`,
`'2 H+ + 2 e- = H2'`, `'e- + h+ = 0'`. Terms are separated by ` + ` with spaces, which keeps it
apart from charges, a leading integer is a coefficient, and an empty side is `0`.

- In `bulkReactions`, it becomes `nu`. A participant on both sides is an error (give the net
  reaction).
- At a face, each participant goes to the side whose material holds it: on a metal's side, only
  its carrier. Participants that aren't species are fixed-activity neutrals, with their μ in
  `fixed`, as in the plain form. A species held on **both** sides is ambiguous, so it's labelled
  with its side or its material: `'Li+(left) = Li+(right)'`, `'Li+(electrolyte) = Li+(graphite)'`.

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

## Alignment from vacuum levels

`vacuumLevel`, `vacuumDipole` and `vacuumZeroCharge` turn vacuum-level estimates (electron
affinities, work functions, Trasatti's absolute SHE) into an interface's alignment. See the
[alignment guide](alignment.md).
