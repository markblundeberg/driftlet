# Device reference

A device is defined by one plain object (serialisable, so it can be posted to a Web Worker),
passed to `new Device(def)`. Construction validates everything and throws a `DeviceError`
naming the offending path (`contacts.left.species.Na+.offset: …`). A field that doesn't belong
where it's written (a typo, or a value in the wrong object) is an error too, never ignored. Physically meaningful
choices are never defaulted silently: reference concentrations, standard potentials and
interface alignments must all be given.

Units are SI throughout: m, mol/m³, J/mol, V, A/m², F/m², s. Helpers in `units` convert
common ones (`units.molar`, `units.perCm3`, `units.eV`, `units.nm`, …). Constants are exported
too (`FARADAY`, `GAS_CONSTANT`, `EPS0`, …).

## Top level

| Key | Meaning |
|---|---|
| `T` | temperature, K (default 298.15) |
| `species` | list of mobile species |
| `materials` | named materials: bulk properties |
| `regions` | the device, left to right: each a material plus length, fixed charge, composition |
| `interfaces` | one entry per face between consecutive regions (optional where defaults apply) |
| `bulkReactions` | homogeneous reactions |
| `contacts` | `{ left, right }` |
| `ports` | internal ports: outside phases exchanging with a window of nodes |
| `grid` | default spacing for every region |

## Species

```js nocheck
species: [
  { name: 'e-', z: -1 },
  { name: 'Na+', z: 1, cRef: 1000 }, // optional default reference concentration, mol/m³
]
```

`z` is any integer, including 0 for neutral species. A species can be absent from some
materials; its outputs are `NaN` there.

## Materials

```js nocheck
materials: {
  water: {
    epsr: 78.5, // relative permittivity; 0 makes the material strictly neutral
    species: {
      'Na+': { D: 1.33e-9, mu0: -261.9e3, cRef: 1000 },
    },
  },
}
```

Per species present: `D` (m²/s, the coefficient in the flux $`N = -(D c/RT)\nabla\bar\mu`$), `mu0`
(standard chemical potential, J/mol) and `cRef` (the concentration that `mu0` refers to; falls back
to the species' `cRef`, and must exist one way or the other). Species not listed are absent.

In a material, only *neutral combinations* of standard potentials mean anything physically.
Shifting every `mu0` by $`z_i F s`$ just moves the material's $`\phi`$, its bookkeeping anchor. See
[conventions](conventions.md).

`epsr: 0` makes a material **strictly neutral**: Poisson is replaced by local neutrality, and $`\phi`$
there is only a bookkeeping multiplier. Use it for macroscopic systems whose double layers you
don't want to resolve, such as electrolytes and mixed conductors. (Metals have their own kind of
material, below.)

**Conductors** (metals, and fast ion conductors) are their own kind of material:

```js nocheck
materials: { Au: { conductor: { species: 'e-', conductivity: 4.1e7 } } } // S/m
```

A conductor has one mobile carrier, and its only unknown is that carrier's $`\bar\mu`$: for a metal,
the Fermi level. There's no $`\varepsilon`$, no `mu0` or `cRef`, and $`\phi`$ is undefined inside. Its
bulk is neutral and incompressible, conduction is ohmic ($`J = -\sigma\nabla V`$), and any charge it
holds sits as a sheet at a charged face. Since nothing is stored inside, the carrier's level is
exactly linear across a conductor region, so the region is a single grid cell, whatever its length.
Conductor regions take no `fixedCharge`, `c0`, `grid`, `velocity` or `mixing`. A port on a conductor
attaches to all of it (see [ports](#internal-ports)). See
[Faces next to a conductor](#faces-next-to-a-conductor).

`statistics` (optional) lists non-ideal statistics models, each covering named species:
Fermi–Dirac, lattice gas (crowding), Redlich–Kister, Debye–Hückel, insertion hosts (OCV
curves) and custom functions. Species not listed are ideal. See [statistics](statistics.md).

```js nocheck
statistics: [{ type: 'lattice', species: ['Na+', 'Cl-'], cMax: 5000 }]
```

## Regions

```js nocheck
regions: [
  { name: 'n', material: 'Si', length: 2e-6, fixedCharge: 16.6 * FARADAY }, // C/m³
  { material: 'water', length: 1e-6, c0: { 'Na+': 10, 'Cl-': 10 } },
]
```

- `fixedCharge`: immobile charge density (doping, ionomer), C/m³. Default 0. In an insertion
  host it's balanced by background electronic carriers.
- `c0`: initial concentrations, mol/m³. Species connected to a contact start from that contact's
  level and don't need it. Any other species (blocked everywhere, or only made and consumed by
  reactions) does: its `c0` fixes the amount it conserves. Concentrations are positive ($`\bar\mu`$ is
  logarithmic in them); to have none of a species in a region, leave it out of that region's
  material. In a strictly neutral region, $`\sum z \cdot c_0`$ must be zero. Where nothing else sets
  $`\phi`$ at the start, a region beside an electrode starts with the electrode's first reaction at
  equilibrium (the electrode at its open-circuit level); otherwise $`\phi`$ carries over from the
  left.
- `velocity`: imposed flow toward +x, m/s (default 0). It carries every mobile species
  ($`D > 0`$) along with the fluid. It's uniform within the region: in strict 1D, incompressible
  flow is the same everywhere. A step in velocity at a face means solvent enters or leaves
  there sideways *without* its solutes, which still cross the face by flux continuity: an ideal
  ultrafiltration membrane, where solutes pile up (concentration polarisation). It isn't a
  sideways inflow or outflow of solution.
- `mixing`: eddy (turbulent) mixing diffusivity $`D_{\mathrm{mix}}`$, m²/s (default 0). It mixes
  composition without carrying current, and does nothing at equilibrium (see
  [conventions](conventions.md#transport)). A strongly mixed region is the local stand-in for a
  well-stirred bath.
- `grid`: per-region override of the grid options (below).

## Interfaces

One entry per face, `interfaces[f]` sitting between `regions[f]` and `regions[f + 1]`. All
fields are optional where the defaults apply.

```js nocheck
interfaces: [
  {
    step: { species: 'e-', value: units.eV(0.25) }, // alignment: or { dipole: volts }
    phi: 'pinned', // 'pinned' | 'neutral' | { type: 'capacitive', C }
    sheetCharge: 0, // C/m²
    species: { 'Cl-': 'blocked', 'Li+': { type: 'conductance', G: 50 }, 'Na+': 'blocked' },
    reactions: [{ left: { 'Na+': -1 }, right: { 'Na+': 1 }, k0: 1e-3, alpha: 0.5 }], // Na⁺ crosses only by this
  },
]
```

**Alignment** (exactly one, for a face between *different* materials under a `pinned` or
`capacitive` law):
- `step: { species, value }`: the step in that charged species' standard level across the face,
  right minus left, $`(\mu^\circ_R + zF\phi_R) - (\mu^\circ_L + zF\phi_L)`$, in J/mol. For electrons
  that's the conduction-band offset.
- `dipole`: the $`\phi`$ jump, right minus left, in each material's own anchoring, in volts.
There is no default: omitting it is an error. When vacuum-level estimates are all you have,
`vacuumDipole` from `driftlet/kit` turns them into a `dipole` (see the
[alignment guide](alignment.md)). A face between regions of the same material
defaults to no dipole.

**Electrostatic law** `phi`:
- `'pinned'`: $`\phi`$ jumps by the alignment. The default, and exact when the grid resolves the
  double layers on both sides. Not between two $`\varepsilon = 0`$ materials: with no field on either side,
  nothing would determine the face's charge.
- `'neutral'`: no charge at the face and a free jump, set by neutrality on each side (Donnan).
  The alignment drops out and must not be given. This is the default between two $`\varepsilon = 0`$
  materials, and next to an insertion host (where $`\phi`$ is undefined).
- `{ type: 'capacitive', C }`: a Helmholtz layer, $`D = -C \cdot (\Delta\phi - \mathtt{dipole})`$, C
  in F/m².

**Species laws** `species` (default: local equilibrium, $`\bar\mu`$ continuous, where the species is
present on both sides; blocked otherwise): `'equilibrium'`, `'blocked'`, or
`{ type: 'conductance', G }` ($`J = G \cdot (V_L - V_R)`$, G in S/m², charged species). A species that
takes part in a reaction at the face and exists on both sides has no default: give its link
(`'blocked'` if it crosses only through the reaction), since free crossing alongside would
short-circuit the kinetics.

**Reactions** `reactions` at the face, each written as an equation or as its participants on
each side, with signed stoichiometric coefficients ($`\nu < 0`$ consumed, $`\nu > 0`$ produced by the forward
reaction):

```js nocheck
{ equation: 'Fe3+ + e- = Fe2+', k0: 1e-4, alpha: 0.5 }                          // a metal on one side
{ equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-3 }               // Ag⁺ + e⁻ ⇌ Ag(s)
{ equation: 'Li+(left) = Li+(right)', k0: 1e-3 }                                // Li⁺ transfer
{ left: { 'e-': -1 }, right: { 'Fe3+': -1, 'Fe2+': 1 }, k0: 1e-4, alpha: 0.5 }  // by sides
```

In an equation, terms are separated by ` + ` (with spaces), sides by `=`, a coefficient is an
integer and a space (`2 e-`), and an empty side is `0`. Each participant goes to the side whose
material holds it (on a conductor's side, only its carrier); a species held on both sides is
labelled with its side, `'Li+(left)'`, or its material, `'Li+(graphite)'`.

For example, an intercalation electrode taking Li⁺ from the electrolyte through a
desolvation (or SEI) step has Li⁺ on both sides, with `species: { 'Li+': 'blocked' }` and the
transfer reaction above.

A participant is a species present on its side (on a conductor's side, only its carrier) or a
fixed-activity neutral, given by its $`\mu`$ in `fixed` (J/mol; the side doesn't matter). Charge must
balance. The rate per area, with $`a = A/RT = -\sum \nu\bar\mu/RT`$ over all participants, is
Butler–Volmer,

```math
r = k_0 \prod_{\nu<0} (c/c_{\mathrm{ref}})^{|\nu|(1-\alpha)} \prod_{\nu>0} (c/c_{\mathrm{ref}})^{\nu\alpha} \left(e^{\alpha a} - e^{-(1-\alpha)a}\right)
```

with $`k_0`$ in mol/(m²·s) and $`\alpha`$ (default 0.5) between 0 and 1. A conductor's carrier has
activity 1 and no factor. That's mass action with rate constants that depend on the electrical part
of the affinity: exactly zero at $`A = 0`$, whatever $`k_0`$ and $`\alpha`$. Each participant is made or
consumed at its side's edge node, just behind any Stern layer, so Frumkin effects arise by
themselves. A face can carry several reactions; a conductor coupled to two couples settles at their
mixed potential. Each solution reports the rates as `interfaces[f].rates`.

### Faces next to a conductor

A conductor has no $`\phi`$, so a face beside it can't take a `dipole` or `step` alignment. Its $`\phi`$
law is:

- `'neutral'` (no charge at the face): the default next to an $`\varepsilon = 0`$ material, and the only
  choice between two conductors; or
- `{ type: 'capacitive', C }` with `zeroCharge`: the other side's $`\phi`$ is tied to the conductor's
  level $`V_F`$ (the carrier's $`\bar\mu`$ as a voltage), with displacement
  $`C \cdot (V_F - \mathtt{zeroCharge} - \phi_{\mathrm{edge}})`$ toward the other side, which is the
  conductor's surface charge. `zeroCharge` is $`V_F - \phi_{\mathrm{edge}}`$ at zero charge
  (`vacuumZeroCharge` in `driftlet/kit` estimates it from a work function).

A `pinned` law isn't offered here: an internal conductor holds its surface charge in a
finite capacitance (a large C approaches the pinned limit). A pinned barrier is still available
with the metal as a contact.

The conductor's carrier can continue across as a species (e⁻ into a semiconductor,
`'equilibrium'` by default). An electrode is a conductor region with `reactions` at its face,
the metal's electrons taking part on its side at its Fermi level, behind a contact that holds
the electrons (`{ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }`). A
floating conductor region with reactions on both faces is a bipolar electrode.

## Bulk reactions

```js nocheck
bulkReactions: [
  { equation: 'e- + h+ = 0', kf: { Si: 2e8 } },
  { equation: 'H+ + OH- = H2O', fixed: { H2O: -237.13e3 }, kf: { water: 1.4e8 } },
  { nu: { 'H+': -1, 'OH-': -1, H2O: 1 }, fixed: { H2O: -237.13e3 }, kf: { water: 1.4e8 } }, // the same, as nu
]
```

`equation`, or `nu` with signed stoichiometric coefficients ($`\nu < 0`$ consumed, $`\nu > 0`$ produced by
the forward reaction); a participant on both sides of an equation is an error (give the net
reaction). Participants that aren't species are fixed-activity neutrals, given by their $`\mu`$ in
`fixed` (J/mol). Charge must balance. The rate is mass action,
$`r = k_f \prod_{\nu<0} c^{|\nu|} (1 - e^{-A/RT})`$ with $`A = -\sum \nu\bar\mu`$ the affinity. That's
$`k_f \prod c_R - k_b \prod c_P`$, with $`k_b`$ fixed by the standard potentials, so equilibrium is exactly
$`A = 0`$. `kf` maps material
names to forward rate constants (units making r mol/(m³·s)), and the reaction runs only in
those materials.

**Generation** is a reaction from a reservoir: photogeneration is
$`\text{photon} \to \mathrm{e}^- + \mathrm{h}^+`$, with the photons a fixed participant whose $`\mu`$
sits well above the gap, so that $`e^{-A/RT}`$ is negligible and the rate is `kf` itself, the
generation rate (mol/(m³·s)). It's uniform within each material listed in `kf` (there's no optical
absorption profile). The photons' $`\mu`$ is the honest part: light is a reservoir far from the
device's temperature, and the rate still vanishes at $`A = 0`$.

```js nocheck
bulkReactions: [{ equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: 0.1 } }]
```

## Contacts

A contact is an interface whose far side is an **outside phase with known levels**: think of it as
one more region whose node is fully known. The outside phase's levels form a rigid ladder,
$`V_i = V + \mathtt{offset}_i`$ for charged species ($`\bar\mu`$ given directly for neutral ones). The
offsets are the outside phase's own chemistry, and the external circuit slides the whole ladder by
the terminal voltage V: held, or floating under a current (see [terminals](#terminals)). The laws
joining it to the device are the same as at internal interfaces.

```js nocheck
contacts: {
  left: {
    V: 0, // terminal voltage, V
    terminal: 'e-', // the species whose voltage V is
    species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 }, 'Cl-': 'blocked' },
    phi: 'bulk', // 'bulk' | 'neutral' | 'pinned' | { type: 'capacitive', C }, with zeroCharge beside the last two
  },
  right: { bath: { c: { 'Na+': 10, 'Cl-': 10 }, reference: 'Cl-' } },
}
```

**Species laws** (default `'blocked'`):
- `'equilibrium'`: the species is in equilibrium with the outside phase, so it's held at
  $`V_i = V + \mathtt{offset}`$ (charged) or $`\bar\mu = \mathtt{mu}`$ (neutral, J/mol).
  - The offset defaults to 0 only for the terminal species. Any other species needs one,
    because it's a property of the outside phase: the chemical potential, per charge, of the
    neutral combination the species forms with the terminal species there.
  - For an ion at a reversible electrode,
    $`\mathrm{M}^{n+} + n\,\mathrm{e}^- \rightleftharpoons \mathrm{M(s)}`$ gives
    $`\mathtt{offset} = \mu_{\mathrm{M}}/(nF)`$, the metal's own chemical potential per charge, which
    is 0 on the usual table convention. ($`E^\circ`$ isn't the offset: it's already carried by the
    ion's $`\mu^\circ`$ in the solution.)
  - At a metal, holes and electrons share the metal's voltage: `offset: 0`.
- `{ type: 'conductance', G, offset }` (charged species): ohmic exchange with the outside level
  at $`V + \mathtt{offset}`$, $`J = G \cdot (V_{\mathrm{out}} - V_i)`$, G in S/m².
- `{ type: 'exchange', k, mu }` (neutral species):
  $`N_{\mathrm{in}} = k \cdot (\mu_{\mathrm{out}} - \mu)/RT`$, k in mol/(m²·s).
- `'blocked'`: no flux.

**`bath`** (instead of `species` and `phi`): the outside phase is a neutral composition `c`
held in place by a charged `reference` species (as for a reversible reference electrode, e.g.
Cl⁻ for Ag/AgCl), which is the terminal. Every bath species is in equilibrium at the level its
composition implies, and the $`\phi`$ law is `'bulk'`. An optional `offset` places the reference
species relative to V.

**`phi`**, required whenever any species or reaction connects at the contact (default
`'neutral'`):
- `'bulk'`: the end node is plain bulk, locally neutral, with no double layer at the contact
  (ohmic contacts, baths). The outside takes whatever surface charge that needs.
- `'neutral'`: no charge at the face ($`D = 0`$), as for an internal `'neutral'` face.
- `{ type: 'capacitive', C }` with the contact's `zeroCharge`: a gate or Stern layer to a conductor
  at the terminal voltage V. The displacement into the device is
  $`C \cdot ((V - \mathtt{zeroCharge}) - \phi_{\mathrm{edge}})`$, so `zeroCharge` is the value of
  $`V - \phi_{\mathrm{edge}}`$ at which the interface carries no charge: the potential of zero charge
  (pzc) of an electrode, or for a gate the flat-band voltage less the semiconductor's bulk $`\phi`$
  ($`V_{\mathrm{FB}} = \mathtt{zeroCharge} + \phi_{\mathrm{bulk}}`$).
- `'pinned'` with `zeroCharge`: the $`C \to \infty`$ limit,
  $`\phi_{\mathrm{edge}} = V - \mathtt{zeroCharge}`$. This is what a "fixed $`\phi`$" boundary honestly
  means. For example, a Schottky barrier $`\phi_B`$ on n-type material with
  $`\mu^\circ_{\mathrm{e}^-} = 0`$ is `zeroCharge: φ_B`.

Like every alignment, `zeroCharge` is a property of that interface. To estimate it from vacuum
levels, `vacuumZeroCharge(def, W, inside)` from `driftlet/kit` takes the work function W of the
conductor at the terminal voltage and the inside material's anchor and offset (see the
[alignment guide](alignment.md)).

Contacts carry no reactions. An electrode reaction belongs at the face of a conductor region,
with the contact behind it holding the conductor's electrons (see
[faces next to a conductor](#faces-next-to-a-conductor)).

## Internal ports

A port is an outside phase with known levels, like a contact's, attached to a window of nodes
inside one region instead of at an end. It's the 1D stand-in for whatever feeds or drains a
species sideways: source and drain grounding a MOS channel, salt injected mid-solution, a
reference electrode's reservoir.

```js nocheck
ports: [{
  name: 'channel',
  region: 'Si',              // name or index
  from: 495e-9, to: 500e-9,  // window, m from the region's left end (default: the whole region)
  V: 0,                      // the port's terminal voltage (or I, or V and R: see terminals)
  terminal: 'e-',
  species: {
    'e-': 'equilibrium',                         // μ̄ held at V_i = V + offset throughout the window
    'Na+': { type: 'conductance', G: 1e9, offset: 0.2 }, // source G (V_out − V_i)/(zF) per volume, G in S/m³
    O2: { type: 'exchange', k: 1e-3, mu: -2e3 },  // neutral: source k (μ_out − μ)/RT per volume, k in mol/(m³·s)
  },
}]
```

Links and offsets are those of contacts, per volume instead of per area: 0 by default only for
the terminal species, and an absolute `mu` for neutral species. A port is a terminal like a
contact: held at a voltage, behind a resistance, or driven by a current (a reference electrode
is a port at `I: 0`). Each solution reports `ports[k]`, with `{ name, V, flux, current }`: what
the port brings into the device.

A held (`'equilibrium'`) level leaves the device's two end nodes to their contacts.

A port driven by a current (`I`, or a waveform) floats to whatever voltage delivers it, so a
conductance link's `G` doesn't set how much enters, only how it's shared across the window: in
proportion to each node's level below the port's. Any large value spreads it evenly. At `I: 0`
the port still ties that species' level together across its window (a large `G` shorts the
window for it), so keep such windows narrow.

On a [conductor region](#materials), a port is a wire to the whole conductor: it takes no window,
only the conductor's carrier, and a conductance link's `G` is per area (S/m², a resistance
$`R \cdot A`$ to the port's voltage as $`G = 1/(R \cdot A)`$). That's how a floating electrode is tied
to ground through a resistor.

For example, in a 1D MOS capacitor without generation, inversion electrons can only arrive by
minority-carrier diffusion from the back contact, which can take weeks. The inversion layer is
then effectively floating, and its Fermi level is undetermined to within round-off, so no
steady state can be computed there. A port holding the electrons beside the oxide at the
channel's potential anchors it, as source and drain would, and the device then shows the
low-frequency C–V.

## Terminals

The two contacts and every port are the device's terminals, named `left`, `right` and by each
port's `name`. Each is driven on its own:

```js nocheck
contacts: {
  left: { V: 0, ... },               // held at a voltage
  right: { I: -2, ... },             // driven by a current into the device, A/m² (I: 0 is open circuit)
},
ports: [{ name: 'ref', I: 0, ... },  // a reference electrode: no current, its voltage read off
        { name: 'wire', V: 0, R: 1e-3, ... }], // a source V behind a series resistance R (Ω·m²)
```

- **`V`**: the terminal's voltage, the shift of its outside phase's ladder, as the voltage of
  its `terminal` species (a gate's is its own metal's electrons). Neither `V` nor `I` given
  means `V: 0`.
- **`I`**: the current into the device through this terminal (conduction plus displacement);
  the voltage floats and is solved for.
- **`V` and `R`**: a source behind a resistance, $`I = (V - V_{\mathrm{terminal}})/R`$.

At least one terminal must be held at a voltage, or the device's overall level floats. In
steady state the terminal currents sum to zero. A contact can be driven by a current only if
something passes it (a linked species, or a gate's displacement).

**Waveforms.** `V` and `I` can be piecewise linear in time:
`{ t: [0, 1, 2], values: [0, 0.5, 0], repeat: true }` (s and V, or s and A/m²), constant beyond the
points, or periodic with period $`t_{\mathrm{last}} - t_0`$ when `repeat` is set: a triangle wave is a
cyclic voltammogram. A **step** is two points at one time,
`{ t: [0, 1e-3, 1e-3], values: [0, 0, 0.5] }`, and a repeating waveform whose last value differs
from its first jumps back each period (a sawtooth); a square wave is steps both ways. Write a jump
as a step rather than as a very fast ramp, which the time-step control has to resolve. A time step
ending at t sees the sources as they are just before t (so it ends on a jump's near side, and the
next step carries the jump); `solve()` and `impedance()` use the values at the present time, after
any jump there. `advance()` lands on every breakpoint and restarts its time stepping there. A
device's time starts at 0 when it's made, and waveforms are read against it: `solve()` doesn't
advance it, `advance()` and `step()` do.

## Grid

`{ hmin, hmax, ratio = 1.2, minCells = 8 }`, per device, overridable per region. With no grid
options at all, each region is graded from `hmin` = 1/1000 to `hmax` = 1/20 of its length. Cells grow
geometrically from `hmin` at both ends of each region up to `hmax`, scaled to fit exactly. Each
region boundary becomes a pair of nodes at the same x, one per side. The grid is never refined
automatically. If a double layer that the model resolves is coarser than the local Debye
length, the solution's `warnings` say so.

Grade the grid toward wherever a profile is steep: double layers, and electrodes where a species
is depleted. Near a limiting current the depleted species' profile is steep in a thin layer at
the electrode, and a uniform grid there overshoots: a 100 µm silver nitrate cell on 100 uniform
cells exceeds the limiting current by 1.4% at 0.5 V, while `{ hmin: 10e-9, hmax: 2e-6 }` stays
within 0.1% with as many nodes. Where a species carrying the current is steep across a region's
end cells in this way, the solution's `warnings` say so, with a rough estimate of the excess.

## Using a device

```js nocheck
const dev = new Device(def);
const sol = dev.solve();                      // steady state (or equilibrium)
dev.set({ contacts: { right: { V: 0.3 } } }); // deep-merged change; the state is kept as a warm start
const sol2 = dev.solve();
const tr = dev.step(1e-6);                    // one backward-Euler step (auto-subdivided if needed)
const tr2 = dev.step(1e-6, { method: 'bdf2' }); // second-order BDF2 once a previous step exists
const run = dev.advance(1e-3, { tol: 1e-3 });  // adaptive steps to t = 1 ms
const frame = dev.advance(t1, { budgetMs: 8 }); // …or as far as 8 ms of compute allows
const Z = dev.impedance([1, 10, 100, 1e3]);    // small-signal impedance about the steady state
const now = dev.solution();                   // snapshot of the current state
```

- `solve()` finds the steady state from the current state without advancing time. If every
  species is fed by a contact, it solves the steady equations directly. Otherwise conserved
  amounts (blocked species, reactive moieties) are kept exactly. If a direct solve fails at a
  bias (a cold start far from equilibrium), it solves with both terminals level and ramps the
  right terminal's voltage to its target; with generation reactions (below), it ramps their
  rates up from nearly nothing.
- `step(dt, { method })` advances the transient by dt seconds, halving internally where Newton
  needs it. `method` is `'be'` (backward Euler, the default) or `'bdf2'`.
- `advance(tEnd, opts)` integrates adaptively to `tEnd` with variable-step BDF2, controlling
  the local error per step to `tol` (default 1e-3) in thermal units of every potential ($`\phi`$ and
  each $`\bar\mu/RT`$): roughly 0.1% in concentrations. It lands exactly on `tEnd`. Options: `tol`,
  `dt0` (first step), `dtMax`, `budgetMs` (return after this much wall time, with
  `done: false`), `maxSteps`, `method`. The step size carries over between calls, so an
  animation can call `advance(tNext, { budgetMs })` once per frame. The solution adds `done`,
  `rejected`, and a `trace` of terminal current and voltage after every accepted step.
- `impedance(frequencies, { terminal, profiles })` solves the steady state, then linearises about
  it: $`Z(f) = \delta V/\delta I`$ in Ω·m² at one terminal (`'right'` by default), with I into the
  device. A held terminal's voltage is perturbed, or a driven one's current; the other terminals
  keep their drives (held ones at AC ground, driven ones open). A terminal behind a resistance can't
  be the one measured: the resistance belongs to the external circuit. With `profiles: true`, each
  frequency also returns complex profiles of $`\delta\phi`$, $`\delta\bar\mu`$ and $`\delta c`$ per unit
  excitation.
- `set(patch)` merges plain objects deeply (arrays are replaced). A patch that changes only the
  terminals' drives (`V`, `I`, `R`) updates them in place, cheaply, keeping everything else.
  Otherwise the device is rebuilt: the current state carries over while the grid and species
  are unchanged, or restarts from the regions' `c0`. A carried state keeps what it holds: a
  stretch that was closed stays at its amount, and one that a change closes (or splits) keeps
  what it holds at that moment. So a new `c0` doesn't apply to a carried state; to start over
  from `c0`, make a new `Device`. Either way the time stepping restarts its
  order, as after any discontinuity.

For example, a silver nitrate cell between silver electrodes: its impedance spectrum, then the
current transient after a voltage step.

```js
import { Device } from 'driftlet';

const electrode = (V) => ({ V, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
const dev = new Device({
  species: [
    { name: 'Ag+', z: 1, cRef: 1000 },
    { name: 'NO3-', z: -1, cRef: 1000 },
  ],
  materials: { water: { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
  regions: [{ material: 'water', length: 20e-6, c0: { 'NO3-': 10 } }],
  contacts: { left: electrode(0), right: electrode(0) },
  grid: { hmin: 10e-9, hmax: 0.5e-6 }, // graded: fine at the electrodes, where Ag⁺ depletes
});

const { f, Z } = dev.impedance([0.01, 1, 100, 1e4]);
f.forEach((fk, k) => console.log(`${fk} Hz: Z = ${Z.re[k].toExponential(3)} ${Z.im[k].toExponential(3)}i Ω·m²`));

dev.set({ contacts: { right: { V: 0.05 } } });
const run = dev.advance(1, { tol: 1e-3 }); // one second, adaptively
console.log(`${run.steps} steps; I(0.01 s) ≈ ${run.trace.current[run.trace.t.findIndex((t) => t >= 0.01)].toFixed(1)} A/m², I(1 s) = ${run.current.toFixed(2)} A/m²`);
```

## Solutions

| Field | Contents |
|---|---|
| `x`, `region` | node positions and region index; each interface position appears twice (one node per side) |
| `species`, `regions` | `{ name, z }` for each species; `{ name, material, x0, x1 }` for each region |
| `phi` | bookkeeping $`\phi`$, V (`NaN` where undefined) |
| `c[name]`, `mu[name]`, `muStd[name]` | concentration, $`\bar\mu`$, standard level $`\mu^\circ + zF\phi`$ (`NaN` where absent) |
| `V[name]`, `Vstd[name]` | species voltage $`\bar\mu/(zF)`$ and standard level as a voltage (charged species) |
| `current`, `terminalVoltage` | current toward +x through the device (A/m²) and $`V_{\mathrm{right}} - V_{\mathrm{left}}`$ |
| `terminals[name]` | `{ V, current }` for each terminal (contacts and ports), current into the device |
| `contacts.left/right` | `{ V, flux: {name}, D, current }` at each contact |
| `gates.left/right` | charge on a gate or Stern plate, where the contact is capacitive |
| `ports[k]` | `{ name, V, flux: {name}, current }`: what each internal port brings into the device |
| `interfaces[f]` | `{ dipole, sheetCharge, D, N: {name}, rates }`: what crosses each face by its links, and each face reaction's rate (mol/(m²·s)) |
| `charge` | total charge in the device, C/m² |
| `conservation` | per species stretch: amount, reference, intake through contacts, drift |
| `warnings` | e.g. unresolved double layers, conventions a statistics model relies on, and for a failed solve, where the system is nearly singular |
| `converged`, `iterations`, `steps`, `substeps`, `history`, `time` | solver bookkeeping |
| `done`, `rejected`, `trace` | from `advance()`: whether `tEnd` was reached; rejected steps; `{ t, current, voltage }` per accepted step |

A stretch is a run of regions in which a species is present and connected. Its `drift`
compares its amount with the reference amount plus everything that came in through the
contacts. It's `NaN` where a reaction also makes or consumes the species.
