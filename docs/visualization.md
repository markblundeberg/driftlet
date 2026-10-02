# Reading level diagrams

driftlet draws its results as level diagrams: voltages against position, on one axis, for the
electrons in the metals, for the reactions in the solutions, and for every charged species. This
page explains what each line means and why it's drawn that way. The first half covers the
electronic levels, Fermi levels and redox levels, which connect directly to electrode potentials
and other familiar quantities. The second half covers the species voltages of every charged
species, ions included. For how to draw them, see [the kit's plotting helpers](kit.md#plotting).

## Electron levels

### The Fermi level

An electron's electrochemical potential $`\bar\mu_{\mathrm{e}^-}`$, written as a voltage, is

```math
V_{\mathrm{e}^-} = -\bar\mu_{\mathrm{e}^-} / F
```

It's the voltage of ordinary circuits. It's what a voltmeter reads between two terminals, it
equalises between metals wired together, it drops along a resistor by $`IR`$, and a cell's voltage
is the difference between its terminals' $`V_{\mathrm{e}^-}`$. In a metal it's the Fermi level, in
volts. A terminal's `V` in driftlet is this voltage, for an electron contact. In a semiconductor
out of equilibrium, electrons and holes each have their own (the quasi-Fermi levels, below).

### Redox levels

A half-reaction $`\mathrm{Ox} + z\,\mathrm{e}^- \rightleftharpoons \mathrm{Red}`$ is in equilibrium when
$`\bar\mu_{\mathrm{Ox}} + z \bar\mu_{\mathrm{e}^-} = \bar\mu_{\mathrm{Red}}`$. That fixes an electron
level, even where no free electrons exist:

```math
V_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red}) = (\bar\mu_{\mathrm{Ox}} - \bar\mu_{\mathrm{Red}}) / (zF)
```

This is the couple's **redox level** (sometimes called its redox Fermi level). At each point
where Ox and Red are present, it's where an electrode exchanging electrons by that couple would
sit if it were in equilibrium with the local composition. `level(sol, half)` evaluates it at
every node.

It belongs on the same axis as the metals' Fermi levels, because then the things that matter
are gaps you can point at:

- **At an electrode**, the step from the metal's $`V_{\mathrm{e}^-}`$ to the couple's level at the
  face is the surface overpotential. It's zero at equilibrium, and its size and sign drive the
  reaction's kinetics. Where the metal sits higher (in volts) than the couple's level, electrons
  leave the solution for the metal: oxidation, an anodic current.
- **Through the solution**, the couple's level slopes wherever current flows. Part of the slope
  is ohmic drop and part is concentration polarisation: the couple's composition changing near
  the electrode.
- **Several couples** in one solution each have their own level. The levels agree only if the
  couples trade electrons fast enough to reach equilibrium with each other, which often they
  don't. An electrode that exchanges with several couples at once settles between their levels:
  a mixed potential.

### Standard redox levels

With every participant at its reference concentration, the same formula gives the couple's
**standard redox level** $`V^\circ_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red})`$
(`level(sol, half, { standard: true })`). The redox level sits above or below it by the activity
term, which makes the Nernst equation local:

```math
V_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red}) = V^\circ_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red}) + \frac{RT}{zF} \ln\frac{a_{\mathrm{Ox}}}{a_{\mathrm{Red}}}
```

The standard level moves with the solution's $`\phi`$, so it slopes with the ohmic drop. The gap
between the two lines is purely composition.

### Electrode potentials

The familiar potentials are gaps between these levels, measured from one more standard level:
the hydrogen couple's, $`V^\circ_{\mathrm{e}^-}(\mathrm{SHE})`$, for 2 H⁺ + 2 e⁻ ⇌ H₂ at unit
activities.

```math
\begin{aligned}
E &= V_{\mathrm{e}^-}(\text{electrode}) - V^\circ_{\mathrm{e}^-}(\mathrm{SHE}) && \quad\text{the electrode potential, "vs SHE"} \\
E_{\mathrm{eq}} &= V_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red}) - V^\circ_{\mathrm{e}^-}(\mathrm{SHE}) && \quad\text{the couple's equilibrium potential} \\
E^\circ &= V^\circ_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red}) - V^\circ_{\mathrm{e}^-}(\mathrm{SHE}) && \quad\text{its standard potential} \\
\eta &= E - E_{\mathrm{eq}} && \quad\text{the overpotential}
\end{aligned}
```

The Nernst equation, $`E_{\mathrm{eq}} = E^\circ + (RT/zF) \ln(a_{\mathrm{Ox}}/a_{\mathrm{Red}})`$, is the
floating form above with the SHE's level subtracted from both sides.

With the data library's standard potentials, which follow the usual tables
($`\mu^\circ_{\mathrm{H}^+} = 0`$ and $`\mu_{\mathrm{H_2}} = 0`$), the SHE's standard level in an
aqueous solution is $`\phi`$ itself. So drawing $`\phi`$ (`traces(sol, { phi: true })`) draws the
SHE's level, and the demos label it that way.

One caution: $`V^\circ_{\mathrm{e}^-}(\mathrm{SHE})`$ is a level in the solution, and it slopes
wherever current flows. So "vs SHE" needs a place. A reference electrode samples the level where
its tip sits, which is why practitioners correct for iR drop, use a Luggin capillary, or add
supporting electrolyte to flatten the bulk. In a cell polarised from wall to wall there's no flat
bulk and no single $`E`$, but the levels are all still there to read.

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
the two terminals' $`V_{\mathrm{e}^-}`$: at open circuit, the gap between the two couples' levels;
under load, less each electrode's overpotential and the solution's share. The
[Daniell demo](https://markblundeberg.github.io/driftlet/demos/daniell.html) draws this and adds it up.

### How they're drawn

| line | what it is |
|---|---|
| thick, solid | a species voltage $`V_i`$: for electrons, the Fermi level |
| thick, dashed | a redox level $`V_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red})`$ |
| thin, dashed | a standard redox level $`V^\circ_{\mathrm{e}^-}(\mathrm{Ox}/\mathrm{Red})`$ |
| thin, solid | a standard level $`V^\circ_i`$ (below) |
| thin, dotted, grey | $`\phi`$, on request: across an insulator, or as the SHE's level (below) |

Electrons are steel blue, cations (and holes) warm colours, anions cool colours, and redox levels
blue-violet.

## Species voltages for every charged species

The Fermi level's definition works for any charged species i:

```math
\begin{aligned}
V_i &= \bar\mu_i / (z_i F) && \quad\text{the species voltage} \\
V^\circ_i &= \phi + \mu^\circ_i / (z_i F) && \quad\text{its standard level}
\end{aligned}
```

$`V_i`$ is a way of drawing $`\bar\mu_i`$, not new physics: the equations, the docs and the API work
in $`\bar\mu_i`$, which also covers neutral species (where $`V_i`$ is undefined). For electrons, $`V_i`$
is the Fermi level. For an ion, it's the voltage an electrode reversible to that ion would read: a
silver wire in equilibrium with the solution's Ag⁺ reads $`V_{\mathrm{Ag}^+}`$. The standard level
is where $`V_i`$ would sit at the reference concentration, so for ideal species

```math
V_i - V^\circ_i = \frac{RT}{z_i F} \ln\frac{c_i}{c_{\mathrm{ref},i}}
```

and for non-ideal ones the activity replaces $`c_i/c_{\mathrm{ref},i}`$.

Drawing every species this way, ions alongside electrons, shows things that are otherwise hard
to see:

- **Each species' driving force is the slope of its own line.** A species moves down its own
  $`\bar\mu_i`$, drift and diffusion together, so on the voltage plot the current it carries flows
  down slopes in its $`V_i`$ (cations and holes move downhill, electrons and anions uphill), and in
  equilibrium its $`V_i`$ is flat. There's no need to split that one slope into an electric field
  and a concentration gradient (a split that depends on how $`\phi`$ was chosen, and so means little
  on its own).
- **Interfaces read at a glance.** Where a species crosses a face freely, its $`V_i`$ is
  continuous. A step in $`V_i`$ means something resists its transfer and dissipates energy there.
  The standard levels step at faces along with $`\phi`$ (Donnan potentials, liquid junctions,
  dipoles), but those steps are bookkeeping: they cancel from anything measured.
- **Reactions are levels lining up.** A redox level is a weighted combination of species
  voltages. For Fe³⁺ + e⁻ ⇌ Fe²⁺,
  $`V_{\mathrm{e}^-}(\mathrm{Fe}^{3+}/\mathrm{Fe}^{2+}) = 3V_{\mathrm{Fe}^{3+}} - 2V_{\mathrm{Fe}^{2+}}`$.
  For a plating couple with $`\mu_{\mathrm{M}} = 0`$, the couple's level coincides with the ion's
  own line: $`V_{\mathrm{e}^-}(\mathrm{Ag}^+/\mathrm{Ag}) = V_{\mathrm{Ag}^+}`$.
- **Semiconductor band diagrams are already this picture.** $`V_{\mathrm{e}^-}`$ and
  $`V_{\mathrm{h}^+}`$ are the quasi-Fermi levels, and their standard levels are the band edges.
  Since $`V = -E/e`$ for electrons, the conduction band sits *below* the valence band in volts.
- **Small-signal response has the same rails.** Jamnik and Maier's exact equivalent circuit for
  transport has a rail for each species' $`\bar\mu_i/(z_i F)`$, with resistors along it for
  transport and chemical capacitances from it to a $`\phi`$ rail for storage (J. Jamnik and J. Maier,
  [Phys. Chem. Chem. Phys. 3, 1668 (2001)](https://doi.org/10.1039/b100180i)). The
  [impedance demo](https://markblundeberg.github.io/driftlet/demos/impedance.html) draws the voltages on those rails.

Different species' levels sit volts apart, because their $`\mu^\circ`$ differ by hundreds of
kJ/mol. The gap between two species' lines carries a constant fixed by convention (below), so
moving one species' lines up or down loses nothing. `traces()` does this with `shifts`, and marks
the shifted lines ⌇.

## What's left out, and why

### φ

The electrostatic potential $`\phi`$ isn't drawn by default, though it's in the same volts. Inside a
material, $`\phi`$ is a bookkeeping convention: each material's standard potentials fix where its
$`\phi`$ sits, and any consistent choice gives the same physics (see [conventions](conventions.md)).
Its steps between materials can't be measured, and inside a metal it has no measurable meaning at
all. And while the electron diagram is in volts, every line on it should be an electronically
meaningful level; $`-e\phi`$ isn't the energy of any electron.

What $`\phi`$ would show is drawn anyway, more usefully. The standard levels
$`V^\circ_i = \phi + \mu^\circ_i/(z_i F)`$ of a material move together with its $`\phi`$, so their
bending ladder carries the same variation in space within each material. Unlike $`\phi`$, the ladder
means something across materials too: its offsets at a face are band offsets and the like, which
belong to the interface and can be measured.

There are two places where $`\phi`$ is the right line to draw:

- **In an insulator** with no mobile charged species, there are no species levels to draw. Band
  edges can stand in, unpopulated and with no Fermi level, but when the point is that nothing
  mobile lives there, the plain sloped $`\phi`$ line says it best: its slope is the field. It's also
  harmless to the reader, since there's no interfacial step in $`\phi`$ in view to mistake for
  something measurable, just a lone line connected to nothing else. The
  [MOS demo](https://markblundeberg.github.io/driftlet/demos/mos.html) draws $`\phi`$ across its oxide this way. It draws $`\phi`$ in its
  silicon too, deliberately anchored at an arbitrary offset from the band edges and stepping at
  the oxide's face, to show the band edges riding on it while its position is a convention.
- **In vacuum**, $`-e\phi`$ is exactly what it claims to be: the energy of an electron at rest there.

`traces(sol, { phi: true })` draws $`\phi`$ everywhere it's defined, and `phi: ['SiO2']` only in the
regions named. Besides the MOS demo, the demos draw it only as the SHE's standard level in water,
which is what it equals with the library's data.

### A reference

Nor is anything drawn against a reference: not $`E(x)`$ vs SHE, not energies from the vacuum. The
diagrams are reference-free. The simulation grounds one terminal at 0 V because it has to put the
zero somewhere, but nothing depends on where: shift the whole diagram up or down and the physics
is unchanged. So pick your own reference. Every vertical gap on the electron diagram is
physical, so choose a level at a point (the SHE's standard level in the stirred bulk, a reference
electrode's Fermi level, a band edge deep in a contact) and read everything else against it.

That's also why "vs SHE" and "from the vacuum" need a place. The SHE's standard level and the
vacuum level both vary in space wherever there's current or charge: the SHE's slopes with the
ohmic drop in a solution, and the vacuum level follows the electrostatics. Plotting $`E(x)`$
against the local SHE level would subtract a sloping line from every other line, mixing changes
in the reference into changes in the thing measured.

On the species diagram there's one technicality. A gap $`V_i - V_j`$ between two species at a point
is free of $`\phi`$ (it cancels), so its changes in space and time are physical. Its absolute value
carries one constant for that pair, though, fixed once by the chemical convention for the
elements' zeros ($`\mu = 0`$ for each element in its standard state, as in every table). For
instance, $`V_{\mathrm{Ag}^+} - V_{\mathrm{e}^-} = \mu_{\mathrm{Ag}}/F`$ at equilibrium, which is
zero only because the convention says so. Electron levels carry no such constant, which is why
the redox levels are drawn as electron levels: the conventions cancel out of
$`V_{\mathrm{e}^-}(\mathrm{Ag}^+/\mathrm{Ag}) = V_{\mathrm{Ag}^+} - \mu_{\mathrm{Ag}}/F`$.

These diagrams are developed at length, from what a voltmeter measures through semiconductors,
electrolytes and electrochemical cells, in
[Electrochemical species band diagrams](https://marklundeberg.com/esbd/).
