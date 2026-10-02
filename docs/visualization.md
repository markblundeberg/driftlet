# Reading level diagrams

driftlet draws its results as level diagrams: voltages against position, on one axis, for the
electrons in the metals, for the reactions in the solutions, and for every charged species. This
page explains what each line means and why it's drawn that way. The first half covers the
electronic levels, Fermi levels and redox levels, which connect directly to electrode potentials
and other familiar quantities. The second half covers the species voltages of every charged
species, ions included. For how to draw them, see [the kit's plotting helpers](kit.md#plotting).

## Electron levels

### The Fermi level

An electron's electrochemical potential μ̄_e⁻, written as a voltage, is

```
V_e⁻ = −μ̄_e⁻ / F
```

It's the voltage of ordinary circuits. It's what a voltmeter reads between two terminals, it
equalises between metals wired together, it drops along a resistor by IR, and a cell's voltage
is the difference between its terminals' V_e⁻. In a metal it's the Fermi level, in volts. A
terminal's `V` in driftlet is this voltage, for an electron contact. In a semiconductor out of
equilibrium, electrons and holes each have their own (the quasi-Fermi levels, below).

### Redox levels

A half-reaction Ox + z e⁻ ⇌ Red is in equilibrium when μ̄_Ox + z μ̄_e⁻ = μ̄_Red. That fixes an
electron level, even where no free electrons exist:

```
V_e⁻(Ox/Red) = (μ̄_Ox − μ̄_Red) / (zF)
```

This is the couple's **redox level** (sometimes called its redox Fermi level). At each point
where Ox and Red are present, it's where an electrode exchanging electrons by that couple would
sit if it were in equilibrium with the local composition. `level(sol, half)` evaluates it at
every node.

It belongs on the same axis as the metals' Fermi levels, because then the things that matter
are gaps you can point at:

- **At an electrode**, the step from the metal's V_e⁻ to the couple's level at the face is the
  surface overpotential. It's zero at equilibrium, and its size and sign drive the reaction's
  kinetics. Where the metal sits higher (in volts) than the couple's level, electrons leave the
  solution for the metal: oxidation, an anodic current.
- **Through the solution**, the couple's level slopes wherever current flows. Part of the slope
  is ohmic drop and part is concentration polarisation: the couple's composition changing near
  the electrode.
- **Several couples** in one solution each have their own level. The levels agree only if the
  couples trade electrons fast enough to reach equilibrium with each other, which often they
  don't. An electrode that exchanges with several couples at once settles between their levels:
  a mixed potential.

### Standard redox levels

With every participant at its reference concentration, the same formula gives the couple's
**standard redox level** V°_e⁻(Ox/Red) (`level(sol, half, { standard: true })`). The redox level
sits above or below it by the activity term, which makes the Nernst equation local:

```
V_e⁻(Ox/Red) = V°_e⁻(Ox/Red) + (RT / zF) ln(a_Ox / a_Red)
```

The standard level moves with the solution's φ, so it slopes with the ohmic drop. The gap
between the two lines is purely composition.

### Electrode potentials

The familiar potentials are gaps between these levels, measured from one more standard level:
the hydrogen couple's, V°_e⁻(SHE), for 2 H⁺ + 2 e⁻ ⇌ H₂ at unit activities.

```
E    = V_e⁻(electrode) − V°_e⁻(SHE)        the electrode potential, "vs SHE"
E_eq = V_e⁻(Ox/Red)    − V°_e⁻(SHE)        the couple's equilibrium potential
E°   = V°_e⁻(Ox/Red)   − V°_e⁻(SHE)        its standard potential
η    = E − E_eq                            the overpotential
```

The Nernst equation, E_eq = E° + (RT/zF) ln(a_Ox/a_Red), is the floating form above with the
SHE's level subtracted from both sides.

With the data library's standard potentials, which follow the usual tables (μ°_H⁺ = 0 and
μ(H₂) = 0), the SHE's standard level in an aqueous solution is φ itself. So drawing φ
(`traces(sol, { phi: true })`) draws the SHE's level, and the demos label it that way.

One caution: V°_e⁻(SHE) is a level in the solution, and it slopes wherever current flows. So
"vs SHE" needs a place. A reference electrode samples the level where its tip sits, which is why
practitioners correct for iR drop, use a Luggin capillary, or add supporting electrolyte to
flatten the bulk. In a cell polarised from wall to wall there's no flat bulk and no single E, but
the levels are all still there to read.

For example, platinum at open circuit in a ferric/ferrous solution:

```js
import { Device, GAS_CONSTANT, FARADAY } from 'driftlet';
import { build, layer, ohmic, bath, half, level, aqueous, metal } from 'driftlet/kit';

// Platinum at open circuit in 2 mM Fe³⁺ and 8 mM Fe²⁺, in 0.5 M KCl.
const iron = half('Fe3+ + e- = Fe2+');
const def = build({
  library: [aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 }), metal('Pt')],
  stack: [
    ohmic({ I: 0 }, ['e-']),
    layer('Pt', 1e-6),
    { reactions: [{ ...iron, k0: 1e-3, alpha: 0.5 }] },
    layer('water', 10e-6),
    bath({ 'K+': 500, 'Cl-': 522, 'Fe3+': 2, 'Fe2+': 8 }, 'Cl-', 0),
  ],
  grid: { hmin: 10e-9, hmax: 0.5e-6 },
});
const sol = new Device(def).solve();
const g = sol.x.length - 1; // in the bulk
const fermi = sol.terminals.left.V; // the platinum's V_e⁻
const she = sol.phi[g]; // V°_e⁻(SHE), with the library's data
const E = fermi - she, Eeq = level(sol, iron)[g] - she, E0 = level(sol, iron, { standard: true })[g] - she;
const VT = (GAS_CONSTANT * 298.15) / FARADAY;
console.log(`E = ${E.toFixed(4)} V vs SHE, E_eq = ${Eeq.toFixed(4)} V, Nernst: ${(E0 + VT * Math.log(2 / 8)).toFixed(4)} V`);
```

A cell's voltage splits the same way along the electrons' path. It's the difference between
the two terminals' V_e⁻: at open circuit, the gap between the two couples' levels; under load,
less each electrode's overpotential and the solution's share. The
[Daniell demo](../demos/daniell.html) draws this and adds it up.

### How they're drawn

| line | what it is |
|---|---|
| thick, solid | a species voltage V_i: for electrons, the Fermi level |
| thick, dashed | a redox level V_e⁻(Ox/Red) |
| thin, dashed | a standard redox level V°_e⁻(Ox/Red) |
| thin, solid | a standard level V°_i (below) |
| thin, dotted, grey | φ, on request: across an insulator, or as the SHE's level (below) |

Electrons are steel blue, cations (and holes) warm colours, anions cool colours, and redox levels
blue-violet.

## Species voltages for every charged species

The Fermi level's definition works for any charged species i:

```
V_i  = μ̄_i / (z_i F)            the species voltage
V°_i = φ + μ°_i / (z_i F)       its standard level
```

For electrons, V_i is the Fermi level. For an ion, it's the voltage an electrode reversible to
that ion would read: a silver wire in equilibrium with the solution's Ag⁺ reads V_Ag⁺. The
standard level is where V_i would sit at the reference concentration, so for ideal species

```
V_i − V°_i = (RT / z_i F) ln(c_i / c_ref,i)
```

and for non-ideal ones the activity replaces c_i/c_ref,i.

Drawing every species this way, ions alongside electrons, shows things that are otherwise hard
to see:

- **Each species' driving force is the slope of its own line.** A species moves down slopes in
  its V_i, and in equilibrium its V_i is flat. That one slope is drift and diffusion together, so
  there's no need to split it into an electric field and a concentration gradient (a split that
  depends on how φ was chosen, and so means little on its own).
- **Interfaces read at a glance.** Where a species crosses a face freely, its V_i is continuous.
  A step in V_i means something resists its transfer and dissipates energy there. The standard
  levels step at faces along with φ (Donnan potentials, liquid junctions, dipoles), but those
  steps are bookkeeping: they cancel from anything measured.
- **Reactions are levels lining up.** A redox level is a weighted combination of species
  voltages. For Fe³⁺ + e⁻ ⇌ Fe²⁺, V_e⁻(Fe³⁺/Fe²⁺) = 3V_Fe³⁺ − 2V_Fe²⁺. For a plating couple with
  μ(M) = 0, the couple's level coincides with the ion's own line: V_e⁻(Ag⁺/Ag) = V_Ag⁺.
- **Semiconductor band diagrams are already this picture.** V_e⁻ and V_h⁺ are the quasi-Fermi
  levels, and their standard levels are the band edges. Since V = −E/e for electrons, the
  conduction band sits *below* the valence band in volts.
- **Small-signal response has the same rails.** Jamnik and Maier's exact equivalent circuit for
  transport has a rail for each species' μ̄_i/(z_iF), with resistors along it for transport and
  chemical capacitances from it to a φ rail for storage (J. Jamnik and J. Maier,
  [Phys. Chem. Chem. Phys. 3, 1668 (2001)](https://doi.org/10.1039/b100180i)). The
  [impedance demo](../demos/impedance.html) draws the voltages on those rails.

Different species' levels sit volts apart, because their μ° differ by hundreds of kJ/mol. The gap
between two species' lines carries a constant fixed by convention (below), so moving one species'
lines up or down loses nothing. `traces()` does this with `shifts`, and marks the shifted lines ⌇.

## What's left out, and why

### φ

The electrostatic potential φ isn't drawn by default, though it's in the same volts. Inside a material, φ
is a bookkeeping convention: each material's standard potentials fix where its φ sits, and any
consistent choice gives the same physics (see [conventions](conventions.md)). Its steps between
materials can't be measured, and inside a metal it has no measurable meaning at all. And while
the electron diagram is in volts, every line on it should be an electronically meaningful level;
−eφ isn't the energy of any electron.

What φ would show is drawn anyway, more usefully. The standard levels V°_i = φ + μ°_i/(z_iF) of a
material move together with its φ, so their bending ladder carries the same variation in space
within each material. Unlike φ, the ladder means something across materials too: its offsets at a
face are band offsets and the like, which belong to the interface and can be measured.

There are two places where φ is the right line to draw:

- **In an insulator** with no mobile charged species, there are no species levels to draw. Band
  edges can stand in, unpopulated and with no Fermi level, but when the point is that nothing
  mobile lives there, the plain sloped φ line says it best: its slope is the field. It's also
  harmless to the reader, since there's no interfacial step in φ in view to mistake for something
  measurable, just a lone line connected to nothing else. The [MOS demo](../demos/mos.html) draws φ
  across its oxide this way. It draws φ in its silicon too, deliberately anchored at an
  arbitrary offset from the band edges and stepping at the oxide's face, to show the band edges
  riding on it while its position is a convention.
- **In vacuum**, −eφ is exactly what it claims to be: the energy of an electron at rest there.

`traces(sol, { phi: true })` draws φ everywhere it's defined, and `phi: ['SiO2']` only in the
regions named. Besides the MOS demo, the demos draw it only as the SHE's standard level in water,
which is what it equals with the library's data.

### A reference

Nor is anything drawn against a reference: not E(x) vs SHE, not energies from the vacuum. The
diagrams are reference-free. The simulation grounds one terminal at 0 V because it has to put the
zero somewhere, but nothing depends on where: shift the whole diagram up or down and the physics
is unchanged. So pick your own reference. Every vertical gap on the electron diagram is
physical, so choose a level at a point (the SHE's standard level in the stirred bulk, a reference
electrode's Fermi level, a band edge deep in a contact) and read everything else against it.

That's also why "vs SHE" and "from the vacuum" need a place. The SHE's standard level and the
vacuum level both vary in space wherever there's current or charge: the SHE's slopes with the
ohmic drop in a solution, and the vacuum level follows the electrostatics. Plotting E(x) against
the local SHE level would subtract a sloping line from every other line, mixing changes in the
reference into changes in the thing measured.

On the species diagram there's one technicality. A gap V_i − V_j between two species at a point is
free of φ (it cancels), so its changes in space and time are physical. Its absolute value carries
one constant for that pair, though, fixed once by the chemical convention for the elements' zeros
(μ = 0 for each element in its standard state, as in every table). For instance, V_Ag⁺ − V_e⁻ =
μ_Ag/F at equilibrium, which is zero only because the convention says so. Electron levels carry no
such constant, which is why the redox levels are drawn as electron levels: the conventions cancel
out of V_e⁻(Ag⁺/Ag) = V_Ag⁺ − μ_Ag/F.

These diagrams are developed at length, from what a voltmeter measures through semiconductors,
electrolytes and electrochemical cells, in
[Electrochemical species band diagrams](https://marklundeberg.com/esbd/).
