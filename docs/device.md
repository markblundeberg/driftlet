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
| `geometry` | the cross-section $`A(x)`$: planar (default), spherical or cylindrical shells, or a profile (see [geometry](#geometry)) |

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
- `c0`: initial concentrations, mol/m³: the starting state of each species given. Species
  connected to a contact don't need it (without it, they start from that contact's level: the
  left contact's if the species reaches it, else the right's, else a port's). Any
  other species (blocked everywhere, or only made and consumed by
  reactions) does: its `c0` fixes the amount it conserves. Concentrations are positive ($`\bar\mu`$ is
  logarithmic in them); to have none of a species in a region, leave it out of that region's
  material. In a strictly neutral region, $`\sum z \cdot c_0`$ must be zero.
  An entry can instead be a profile, `{ x: [...], values: [...] }`, against the device's x (m, not
  the region's own), piecewise linear between its points and constant beyond its ends, so a
  transient can start from a state it didn't simulate reaching: a packet already injected, a
  gradient laid down. Skip `solve()`, which goes straight to the steady state, and `advance()`
  from it. Species without a `c0` start from their contact's level, and $`\phi`$ is chosen node
  by node for local neutrality, so a mobile counter-ion follows the profile (electrons follow an
  injected hole packet); what the cold start leaves unbalanced settles within the first steps.
  That fill-in is right in a quasi-neutral bulk (it's what dielectric relaxation would do), not
  where the rest state carries charge: a packet laid into a junction starts with the junction
  itself unformed. The grid isn't refined to fit a profile: where it's steeper than the cells
  under it, the start is the profile sampled at the nodes, something else, and every solution
  warns of it with the cells that would resolve it (grade toward it, or put a region boundary
  there).
  For a spectator, the profile's integral is the amount it conserves. Where nothing else sets
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
  in F/m². Beside a strictly neutral ($`\varepsilon = 0`$) material, the face's charge is held in
  the cell next to it, which acts as a diffuse layer half that cell wide, $`F^2 \sum_i z_i^2 c_i \, h/2RT`$
  in series with C. So keep that cell coarse rather than fine (a warning says when it costs C more
  than 2%, and how coarse is enough).

**Species laws** `species` (default: local equilibrium, $`\bar\mu`$ continuous, where the species is
present on both sides; blocked otherwise): `'equilibrium'`, `'blocked'`, or
`{ type: 'conductance', G }` ($`J = G \cdot (V_L - V_R)`$, G in S/m², charged species), or
`{ type: 'permeability', P }` (P in m/s; see [membranes](#membranes)). A species that
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

with $`k_0`$ in mol/(m²·s) and $`\alpha`$ (default 0.5) from 0 to 1 inclusive ($`\alpha = 0`$ makes the
forward rate plain mass action, independent of the energetics). A conductor's carrier has
activity 1 and no factor. That's mass action with rate constants that depend on the electrical part
of the affinity: exactly zero at $`A = 0`$, whatever $`k_0`$ and $`\alpha`$. Each participant is made or
consumed at its side's edge node, just behind any Stern layer, so Frumkin effects arise by
themselves. A face can carry several reactions; a conductor coupled to two couples settles at their
mixed potential. Each solution reports the rates as `interfaces[f].rates`.

**Interface recombination through traps** (SRH at a face) replaces `k0` and `alpha` with
`srh: { vn, vp, n1 }`, on a reaction consuming one negative and one positive species and making
none, such as `'e-(left) + h+(right) = 0'`:

```math
r = \frac{n\,p\,(1 - e^{-A/RT})}{(n + n_1)/v_p + (p + p_1)/v_n}, \qquad p_1 = \frac{n\,p\,e^{-A/RT}}{n_1}
```

with n and p each at its own side's edge node, the capture velocities $`v_n`$ and $`v_p`$ in m/s,
and $`n_1`$ (mol/m³) the negative species' concentration with its level at the trap (default:
$`n_1 = p_1 = n_i`$, a midgap trap). Taking $`p_1`$ from the state keeps equilibrium exact whatever
lies between the two sides (a band offset, a $`\phi`$ jump). Where one carrier is plentiful the
rate saturates at the other's capture, $`r \to v_p\,p`$: holes reaching an electron-rich layer
recombine at $`v_p`$ however many electrons wait there, which plain mass action can't do.

**Saturating kinetics** (a transporter, a pump, an enzyme: something that turns over at most so
fast) replaces `k0` and `alpha` with `vmax` (mol/(m²·s)) and `K`, a half-saturation
concentration (mol/m³) for every species the forward reaction consumes:

```math
r = v_{\max} \prod_{\nu<0} \left(\frac{c}{c + K}\right)^{|\nu|} \left(1 - e^{-A/RT}\right)
```

Michaelis–Menten in each substrate, it runs at $`v_{\max}`$ when they're plentiful and the driving
force is large, and it's still exactly zero at $`A = 0`$, so a pump stalls where the free energy it
spends is used up: its static head. The backward rate is the one detailed balance implies, and
isn't itself saturating, so this is for reactions that run forward.

### Membranes

A membrane thin against everything else (a lipid bilayer, 4–5 nm) is a face. Give it a
capacitive $`\phi`$ law, its capacitance (about 0.01 F/m², 1 µF/cm²), and a species law
`{ type: 'permeability', P }` for each ion that crosses it, through channels or carriers. The
flux is electrodiffusion across the membrane in a constant field (Goldman–Hodgkin–Katz),

```math
N = P\,\frac{z u\,(c_L - c_R\,e^{z u})}{e^{z u} - 1}, \qquad u = \frac{\phi_R - \phi_L}{V_T}
```

toward +x, with $`c_L`$ and $`c_R`$ at the face's two edge nodes and $`\phi_R - \phi_L`$ the capacitor's
voltage. (It's Scharfetter–Gummel across the face, with $`D/h \to P`$; where the two sides are
different materials, the step in standard level between them counts with $`z\Delta\phi`$, so the flux
is still exactly zero where $`\bar\mu`$ is level.) The solutions on either side can be strictly
neutral ($`\varepsilon = 0`$), the membrane's charge then held in each edge node's box, or resolve
their diffuse layers ($`\varepsilon > 0`$); either way the membrane potential comes out the same, to
a few microvolts. Pumps and other
transporters go on the same face as reactions with saturating kinetics, for example the
Na⁺/K⁺-ATPase between the outside (left) and the cytoplasm (right):

```js nocheck
{
  phi: { type: 'capacitive', C: 0.01 },
  species: { 'K+': { type: 'permeability', P: 1e-8 }, 'Na+': { type: 'permeability', P: 4e-10 }, 'Cl-': { type: 'permeability', P: 4.5e-9 }, 'A-': 'blocked' },
  reactions: [{ equation: '3 Na+(right) + 2 K+(left) + ATP = 3 Na+(left) + 2 K+(right) + ADP', fixed: { ATP: 50e3, ADP: 0 }, vmax: 5e-7, K: { 'Na+': 10, 'K+': 1.5 } }],
}
```

`fixed` gives ATP's hydrolysis free energy (here 50 kJ/mol, with ADP and phosphate taken
together). Membrane potentials are $`\phi`$ differences between the two solutions, $`\phi_{\mathrm{in}} - \phi_{\mathrm{out}}`$,
read from `sol.phi` (or a `'phi'` probe in a transient). A bath contact's terminal voltage is
its reference species' level instead, what an electrode reversible to that species would read,
so between two different solutions it differs from $`\Delta\phi`$ by that species' Nernst term;
a bath given no reference has its $`\phi`$ as its terminal voltage, and between two of those the
terminals read $`\Delta\phi`$ itself.

### Voltage-gated channels

A permeability can be gated, as Hodgkin and Huxley's Na⁺ and K⁺ channels are: the face lists its
`gates`, each a fraction $`x`$ in [0, 1] with first-order kinetics in the voltage across the face,
$`V = \phi_{\mathrm{right}} - \phi_{\mathrm{left}}`$ at its two edge nodes,

```math
\frac{dx}{dt} = \alpha(V)\,(1 - x) - \beta(V)\,x ,
```

and a link names the gates that scale its $`P`$ (or a conductance link's $`G`$), with their exponents:
$`P\,m^3 h`$ for `{ m: 3, h: 1 }`. (A membrane along a region instead, an axon's, is a
[port](#a-membrane-through-a-window-gated-channels).)
Each of $`\alpha`$ and $`\beta`$ is one of NeuroML's three standard forms, with `rate` in 1/s and `midpoint`
and `scale` in volts:

| `type` | rate |
| --- | --- |
| `'exp'` | $`\mathtt{rate}\; e^{(V - \mathtt{midpoint})/\mathtt{scale}}`$ |
| `'sigmoid'` | $`\mathtt{rate} / (1 + e^{(\mathtt{midpoint} - V)/\mathtt{scale}})`$ |
| `'expLinear'` | $`\mathtt{rate}\; y/(1 - e^{-y})`$, $`y = (V - \mathtt{midpoint})/\mathtt{scale}`$ |

With the outside on the left, $`V`$ is the membrane potential $`\phi_{\mathrm{in}} - \phi_{\mathrm{out}}`$. Hodgkin and
Huxley's squid axon (at their 6.3 °C):

```js nocheck
{
  phi: { type: 'capacitive', C: 0.01 },
  gates: {
    m: { alpha: { type: 'expLinear', rate: 1000, midpoint: -0.040, scale: 0.010 }, beta: { type: 'exp', rate: 4000, midpoint: -0.065, scale: -0.018 } },
    h: { alpha: { type: 'exp', rate: 70, midpoint: -0.065, scale: -0.020 }, beta: { type: 'sigmoid', rate: 1000, midpoint: -0.035, scale: 0.010 } },
    n: { alpha: { type: 'expLinear', rate: 100, midpoint: -0.055, scale: 0.010 }, beta: { type: 'exp', rate: 125, midpoint: -0.065, scale: -0.080 } },
  },
  species: {
    'Na+': { type: 'permeability', P: 1.16e-6, gates: { m: 3, h: 1 } },
    'K+': { type: 'permeability', P: 1.33e-6, gates: { n: 4 } },
    'Cl-': { type: 'permeability', P: 6.4e-9 }, // the leak
    'A-': 'blocked',
  },
}
```

The gates scale permeabilities and nothing else: they hold no charge (the small gating current
of the channels' voltage sensors is left out, as Hodgkin and Huxley left it out), and each ion
still crosses by electrodiffusion down its own $`\bar\mu`$, so the dissipation stays positive whatever
they do. In steady state each gate is at $`\alpha/(\alpha + \beta)`$; each solution reports them as
`interfaces[f].gates` (`{ m, h, n }`) with the voltage they follow, `interfaces[f].V`.

Two things differ from the textbook equations. The currents are GHK's, which rectify: a channel's
conductance depends on the voltage, so a permeability matches a conductance $`g`$ at one voltage
only (above, the chord conductances at rest, −65 mV, for Hodgkin and Huxley's 120, 36 and 0.3
mS/cm²; matched at the reversal potentials instead, the rest depolarises itself and the axon
fires on its own). Their own linear channels are conductance links: $`G\,(V_{i,L} - V_{i,R})`$
toward +x is $`g\,(E_i - V)`$, with $`E_i`$ the Nernst potential at the edge concentrations, so
`{ type: 'conductance', G: 360, gates: { n: 4 } }` is their K⁺ channel exactly
(`hodgkinHuxley({ T, linear: true })` writes all three). Between two baths whose terminal
voltage is their $`\phi`$, a held step is their voltage clamp, and the K⁺ current at its end is
$`g_K n_\infty^4 (V - E_K)`$ to 2e-3. And a cell without a pump has no resting steady state: `solve()` finds its
Donnan equilibrium, every permeant ion level across the membrane, which the real cell reaches
only over hours. Start a transient from the concentrations given (it settles to rest within
milliseconds), or keep the gradients with the Na⁺/K⁺ pump on the same face. Started from the
same state, a squid axon's action potential follows the space-clamped Hodgkin–Huxley equations
with GHK currents to within 1 mV; and in small signals, with the membrane held at rest between
two baths, `impedance()` gives the gates' inductance that Cole measured (the impedance inductive
below ~50 Hz and resonant near 60 Hz), the linearised equations' to 5e-4
([`gates`](../test/gates.test.js) test).

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

**SRH (trap-assisted) recombination** replaces `kf` with `srh`, a map from material names to
`{ tauN, tauP, n1 }`, on a reaction consuming one negative and one positive species and making
none (`'e- + h+ = 0'`):
$`r = n p (1 - e^{-A/RT}) / (\tau_p (n + n_1) + \tau_n (p + p_1))`$, with lifetimes in s and $`n_1`$ as
for [interface SRH](#interfaces) (default midgap). It's what explicit trap species reduce to when
the traps are few and fast, without them.

```js nocheck
bulkReactions: [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1e-6, tauP: 1e-6 } } }]
```

**Generation** is a reaction from a reservoir: photogeneration is
$`\text{photon} \to \mathrm{e}^- + \mathrm{h}^+`$, with the photons a fixed participant whose $`\mu`$
sits well above the gap, so that $`e^{-A/RT}`$ is negligible and the rate is `kf` itself, the
generation rate (mol/(m³·s)). The photons' $`\mu`$ is the honest part: light is a reservoir far from the
device's temperature, and the rate still vanishes at $`A = 0`$.

```js nocheck
bulkReactions: [{ equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: 0.1 } }]
```

A rate constant can also vary with position: give it as a profile `{ x: [...], values: [...] }`
against the device's x (the same tables as `c0`), piecewise linear and constant beyond its ends,
used only at that material's nodes. Each node takes the profile's mean over its box (not its
value at the node), so the total rate is exact on any grid: light absorbed within a fraction of
a cell still generates what it should. The kit's `photogeneration()` writes Beer–Lambert
absorption this way.

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
- `{ type: 'velocity', v, offset }` (or `mu`, for a neutral species): a surface velocity,
  $`N_{\mathrm{in}} = v\,(c_{\mathrm{eq}} - c)`$ with v in m/s and $`c_{\mathrm{eq}}`$ the end node's
  concentration in equilibrium with the outside level, i.e.
  $`v\,c\,(e^{(\bar\mu_{\mathrm{out}} - \bar\mu)/RT} - 1)`$, which keeps the flux's sign that of the
  level difference whatever the statistics. It lies between `'blocked'` (v → 0) and
  `'equilibrium'` (v → ∞): a contact's surface recombination velocity, or thermionic emission
  over a Schottky barrier (below).
- `'blocked'`: no flux.

**`bath`** (instead of `species` and `phi`): the outside phase is a neutral composition `c`
held in place by a charged `reference` species (as for a reversible reference electrode, e.g.
Cl⁻ for Ag/AgCl), which is the terminal. Every bath species is in equilibrium at the level its
composition implies, and the $`\phi`$ law is `'bulk'`. An optional `offset` places the reference
species relative to V. Without a `reference`, V is the bath's own $`\phi`$: an ideal salt bridge
(no junction potential), the convention by which membrane potentials are $`\phi`$ differences.
Between two such baths, a held $`V_{\mathrm{right}} - V_{\mathrm{left}}`$ is a voltage clamp. The bath is the end
material at that composition (which then balances any fixed charge there, as an insertion host's
background carriers do). A solution against a charged end region (an ion exchanger, a gel) names
its own material, `bath: { c, material: 'water' }`: its composition is then neutral by itself, its
levels are that material's, and the end node takes them through a Donnan step.

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

A Schottky contact is that pinned barrier with its carriers passing at their emission velocities:
`{ V, terminal: 'e-', species: { 'e-': { type: 'velocity', v }, 'h+': { type: 'velocity', v: vp, offset: 0 } }, phi: 'pinned', zeroCharge: φ_B }`,
with Bethe's thermionic-emission velocity $`v = A^* T^2/(F N_c)`$ ($`N_c`$ in mol/m³, the data
library's `cRef` for e⁻; about $`2 \times 10^4`$ m/s in n-Si with $`A^*`$ = 112 A/(cm²K²)). The
forward current is then thermionic emission's, $`J = A^* T^2 e^{-\phi_B/V_T}(e^{V/V_T} - 1)`$,
lowered where the drift–diffusion through the depletion region can't keep up (Crowell and Sze's
combination): 4–6% short of it for a 0.85 eV barrier on 1e16 cm⁻³ n-Si from 0.1 to 0.4 V
([`contacts`](../test/contacts.test.js) test). With `'equilibrium'` links instead (an infinite
velocity) the barrier passes whatever the semiconductor brings to it, which is the diffusion
theory of the Schottky diode: there, 17 to 27 times more current. Neither has image-force
lowering or tunnelling.

Like every alignment, `zeroCharge` is a property of that interface. To estimate it from vacuum
levels, `vacuumZeroCharge(def, W, inside)` from `driftlet/kit` takes the work function W of the
conductor at the terminal voltage and the inside material's anchor and offset (see the
[alignment guide](alignment.md)).

Contacts carry no reactions. An electrode reaction belongs at the face of a conductor region,
with the contact behind it holding the conductor's electrons (see
[faces next to a conductor](#faces-next-to-a-conductor)).

## Internal ports

A port is an outside phase with known levels, like a contact's, attached to a window of nodes
inside one region instead of at an end: the 1D stand-in for whatever feeds, drains or charges the
region sideways. It's a [terminal](#terminals) like a contact (held at `V`, behind `R`, or driven
by `I`; a reference electrode is a port at `I: 0`), and in its window it does any of these, in
any combination:

| Field | What it adds | For |
|---|---|---|
| `species` | [levels held or exchanged](#levels-held-or-exchanged), per volume | source and drain grounding a channel, salt injected mid-solution, O₂ from the air |
| `reactions`, `area` | an [electrode spread through the window](#electrodes-spread-through-a-window) | a thin film on a metal, a crevice, a porous electrode |
| `surface` | [species on that electrode's sites](#an-electrodes-surface-adsorbates-and-passive-films) | adsorbates, a passive film |
| `capacitance`, `area` | [charge held across a capacitance](#a-capacitance-through-a-window) | a gate along a channel, a double layer |
| `gates` | [voltage-gated channels](#a-membrane-through-a-window-gated-channels) on its conductance links | an axon's membrane along it |

```js nocheck
ports: [{
  name: 'channel',
  region: 'Si',              // name or index
  from: 495e-9, to: 500e-9,  // window, m from the region's left end (default: the whole region)
  V: 0,                      // or I, or V and R (see terminals)
  terminal: 'e-',            // the species whose level V is
  species: { 'e-': 'equilibrium' },
}]
```

On a [conductor region](#materials), a port is a wire to the whole conductor: it takes no window,
only the conductor's carrier, and a conductance link's `G` is per area (S/m², a resistance
$`R \cdot A`$ to the port's voltage as $`G = 1/(R \cdot A)`$). That's how a floating electrode is tied
to ground through a resistor.

### Levels held or exchanged

```js nocheck
species: {
  'e-': 'equilibrium',                                // μ̄ held at V_i = V + offset throughout the window
  'Na+': { type: 'conductance', G: 1e9, offset: 0.2 }, // source G (V_out − V_i)/(zF) per volume, G in S/m³
  O2: { type: 'exchange', k: 1e-3, mu: -2e3 },         // neutral: source k (μ_out − μ)/RT per volume, k in mol/(m³·s)
}
```

Links and offsets are those of contacts, per volume instead of per area: 0 by default only for
the terminal species, and an absolute `mu` for neutral species. A held (`'equilibrium'`) level
leaves the device's two end nodes to their contacts. A port that exchanges only neutral species
carries no current, so its voltage means nothing: give it no drive.

A port driven by a current floats to whatever voltage delivers it: with a conductance `G` per
volume over a window of width $`w`$, about $`I/(G w)`$ above the window's level for that species. So
`G` decides how the current is shared across the window, not how much enters. A small `G` spreads
it evenly, like a current source; a large one pins the window near the port's level, and the
current enters wherever the species is drawn away fastest. At `I: 0` the port still ties that
species' level together across its window (a large `G` shorts the window for it), so keep such
windows narrow.

An exchange link's source is linear in $`\mu_{out} - \mu`$, so it grows without bound as the
species runs out ($`\ln c`$). A supply limited by diffusion through a film above the window, such
as O₂ reaching a thin layer of water from the air, is linear in $`c`$ instead: $`k(c_{sat} - c)`$. For
that, use a [bulk reaction](#bulk-reactions) with the outside phase as a fixed participant,
`{ equation: 'Air = O2', fixed: { Air: muSat }, kf: { water: profile } }`, whose rate is
$`k_f (1 - c/c_{sat})`$; a profile gives $`k_f`$ its spatial variation ($`D c_{sat}/h^2`$ under a film
of thickness $`h(x)`$).

A held level can also anchor what nothing else holds. In a 1D MOS capacitor without generation,
inversion electrons can only arrive by minority-carrier diffusion from the back contact (weeks),
so the inversion layer's Fermi level is undetermined to within round-off and no steady state can
be computed. A port holding the electrons beside the oxide, as source and drain would, anchors it,
and the device shows the low-frequency C–V.

### Electrodes spread through a window

A port can hold an electrode: a metal whose carrier sits at the port's level, reacting all
through the window rather than at one face. That's the floor under a thin film of electrolyte, the
walls of a crevice or a pit, or the matrix of a porous electrode, wherever the metal conducts well
enough to be one level. Give the port its `reactions` (Butler–Volmer, as at a face), the carrier as
its `terminal` (the species need be in no material), and the electrode's `area` per volume of the
window:

```js nocheck
ports: [{
  name: 'iron', region: 'film', terminal: 'e-', V: 0,
  area: 1 / 100e-6,           // m²/m³: a film 100 µm thick on the metal; or a profile against x
  reactions: [
    { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 5e-6, alpha: 0.5 },
    { equation: 'O2 + 2 H2O + 4 e- = 4 OH-', fixed: { H2O: -237.13e3 }, k0: 3e-9, alpha: 0.125 },
  ],
}]
```

Each node's rate per area is the face's law, with the region's species at that node and the
carrier at activity 1 and level $`V`$. Times the area per volume, it makes and consumes the
region's species there, and the electrons it takes or gives are the port's current. Held at a
voltage, the port is a potentiostat. At `I: 0` (or with nothing else to carry current), the
reactions settle at the mixed potential, the corrosion potential where anodic and cathodic
currents cancel over the window, while each node can be a net anode or a net cathode; a floating
one starts there, for the start's composition (with a double layer, it starts uncharged and
charges toward it). An electrode that's the device's only path for current is at `I: 0` whatever
its drive: hold it at `V: 0`, which also fixes the level the water's φ is read against. Per volume of the window a rate is
`rates[q][j] * area[j]`, and per area of a film's floor (thickness h, `area` 1/h) the rate itself.

In 1D, a film on a metal is a slice along the metal: the film's thickness enters through `area`,
and through a bulk reaction for anything that reaches the film from above (O₂ from the air, see
above). Where the thickness varies, give the device its [cross-section](#geometry) too, so that ions
travelling along the film squeeze through where it thins: for a drop on a metal, seen from its
centre, $`A = 2\pi r\,h(r)`$ with `area` $`1/h(r)`$. The region must be strictly neutral
($`\varepsilon = 0`$), as in porous-electrode theory: each spot's double layer is below the grid, so
the reactions can't leave the solution charged (give the port a
[capacitance](#a-capacitance-through-a-window) for the double layer's charge).

Validated against the transmission line (linear kinetics along a bar, test/ports.test.js) and the
Wagner–Traud mixed potential of two Butler–Volmer couples.

### An electrode's surface: adsorbates and passive films

An electrode port can carry species on its surface, sharing its sites: an adsorbed intermediate,
or a passive film taken as a coverage. Each has a coverage $`\theta`$ at every node of the window,
Langmuir statistics, $`\mu = \mu^\circ + RT\ln(\theta/\theta_0)`$ with $`\theta_0 = 1 - \sum\theta`$ the
bare fraction, and `capacity` $`\Gamma`$, mol of sites per m² of electrode. Reactions name them like
species, and a reaction with `bare: true` runs only on bare metal: its rate is multiplied by
$`\theta_0`$, both ways, so equilibrium stays exact. Written for adsorption itself, that's Langmuir
kinetics (on at a rate ∝ $`\theta_0`$, off ∝ $`\theta`$).

```js nocheck
ports: [{
  name: 'iron', region: 'film', terminal: 'e-', V: 0, area: 1 / 100e-6,
  surface: { 'Fe(OH)2': { mu0: -490e3, capacity: 2e-5, theta0: 1e-6 } }, // theta0: the start
  reactions: [
    { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 5e-6, alpha: 0.5, bare: true },  // dissolves where bare
    { equation: 'Fe(s) + 2 OH- = Fe(OH)2 + 2 e-', fixed: { 'Fe(s)': 0 }, k0: 1e-4, alpha: 0.5, bare: true }, // the film forms (in minutes)
    { equation: 'O2 + 2 H2O + 4 e- = 4 OH-', fixed: { H2O: -237.13e3 }, k0: 3e-9, alpha: 0.125 }, // on film or not
  ],
}]
```

A film forms where the potential and the pH favour it, through its own thermodynamics, and blocks
what needs bare metal: the dissolution current climbs with potential, peaks, and falls as the film
covers the metal, the active–passive curve. Each node's surface stores $`\Gamma\,d\theta/dt`$ per area
of electrode, so filling it passes a current, and what passes through a surface is conserved with
it. It's a monolayer picture: a real passive oxide is nanometres thick, grows, and breaks down
under chloride, none of which a coverage has.

Cautions:
- **Keep a film's stability within reason.** Blocking the bare fraction to 1e-3–1e-6 is already
  a passive metal; far below that ($`\theta_0 \lesssim 10^{-10}`$) the coverage's η barely affects
  anything and Newton struggles with it (a solution warns).
- **A film's μ° is a fitted number.** A monolayer given a bulk hydroxide's ΔG°f passivates iron far
  into neutral water, so choose it to put the passivation where the metal's chemistry does.
- **Without a double layer, passivation can jump.** With no capacitance the potential follows the
  kinetics instantly, and when the last active patch covers over, the potential where no current
  flows can vanish from the active branch: a discontinuity, where a transient stops. Give the port
  a [capacitance](#a-capacitance-through-a-window) and it's a quick, continuous swing, as on a real
  electrode.

Validated in test/coverage.test.js: the Langmuir isotherm against the electrode's potential (to
1e-9), filling against its exponential and its charge $`F\,\Gamma\,\Delta\theta`$ per area of electrode,
the active–passive curve with blocking exactly $`(1-\theta)`$ times Butler–Volmer, and a steady state
through a surface intermediate where the transient ends.

### A capacitance through a window

A port can hold charge across a capacitance spread through its window: per area of electrode,
$`\sigma = C\,(V - \mathtt{zeroCharge} - \phi)`$ on the port's side (C in F/m²), with `area` per volume
as for reactions, and the window holds the opposite. Each node's charge balance gains $`a\sigma`$:
Gauss's law averaged across the section, $`-\partial_x(\varepsilon\,\partial_x\phi) = \rho + a\sigma`$,
and at $`\varepsilon = 0`$ neutrality $`\rho + a\sigma = 0`$, the mobile species rearranging to supply
it (counter-ions in, co-ions out). The port passes the charging current, $`d(\int a\sigma\,dV)/dt`$,
none in a steady state. Capacitances on overlapping windows add (a channel's top and bottom
gates). A capacitance alone, driven by a current, only charges: it starts uncharged, `advance()`
follows it, and `solve()` refuses (there's no steady state, as for a gate contact driven by one).

```js nocheck
// a thin-film transistor: the gate along a 10 nm channel through 1 mF/m² (35 nm of SiO₂)
ports: [{ name: 'gate', region: 'channel', V: 1.5, area: 1 / 10e-9, capacitance: { C: 1e-3, zeroCharge: 0.3 } }]
// a porous electrode's double layer, alongside its reactions
ports: [{ name: 'iron', region: 'film', terminal: 'e-', V: 0, area: 1e4, capacitance: { C: 0.2, zeroCharge: -0.4 }, reactions: [/* … */] }]
```

A volumetric capacitance $`C^*`$ (F/m³, an organic electrochemical transistor's) goes in per
area of a nominal electrode: `area: 1 / t, capacitance: { C: Cstar * t }` for a film t thick
(`area: 1, C: Cstar` is the same physics, but reads to `describe()` like a unit slip). To read the
gate's voltage as the film's effective one, put the undisturbed film at φ = 0 (its carrier's
`cRef` at its doping, `mu0: 0`) and leave `zeroCharge` at 0. The gate charges the carriers
through their own chemical capacitance too, in series, so $`dQ/dV_G`$ comes out below $`C^*`$
(about 12% for holes at 1e20 cm⁻³ against 40 F/cm³).

That's the gradual-channel approximation of a field-effect transistor, exact for a thin film
with no body (a TFT, an organic or oxide transistor) while the channel is long against the gate
dielectric, and failing where it pinches off (where the 2D field takes over). It's also the
Bernards model of an organic electrochemical transistor (a volumetric capacitance; the
Bernards–Malliaras current is the charge-sheet formula below without its $`V_T`$ term, so they
part near pinch-off, where driftlet keeps the diffusion that formula drops), cable theory
for an axon's membrane, and a porous electrode's double layer. The gate's flat level pins the
channel's standard level through C, as a blocked spectator's flat level does in an electrolyte,
and the current saturates as the carriers run out at the drain, in both.

Validated in test/capacitance.test.js: a thin-film transistor's current against the
charge-sheet model, $`I/W = (\mu/L)[(\sigma_s^2 - \sigma_d^2)/(2C) + V_T(\sigma_s - \sigma_d)]`$, from
below threshold through saturation (second order in the grid); a porous electrode's de Levie
impedance $`\sqrt{r/y}\coth(L\sqrt{ry})`$, the double layer in series with the ions' chemical
capacitance; and the gate's low-frequency impedance against $`dQ/dV_G`$ from steady states.

### A membrane through a window: gated channels

A port with a capacitance can carry [gates](#voltage-gated-channels), as a face can: at each node of
its window, each gate is a fraction open with Hodgkin–Huxley kinetics in the voltage across the
port's capacitance there, $`V_m = \phi - (V - \mathtt{zeroCharge})`$, and its conductance links name
the gates that scale them (`G × m³h`). That's a membrane all along a region, between it and an
outside held by the port: an axon's, with x running along the axon. Per volume of axoplasm, a
cylinder of radius a has $`2/a`$ of membrane (`area`). Give the port the outside's composition as a
`bath`: its V is then the outside's $`\phi`$ (held at 0, with `zeroCharge: 0`), and each linked ion's
level there follows from the composition, $`(\mu^\circ + RT\ln(c_{out}/c_{\mathrm{ref}}))/(zF)`$ against
it, as an `offset` would give it. (Only the linked ions' levels are taken, so the bath needn't list
or balance the rest.) A conductance link's current is then $`G\,(V_{i,out} - V_i)`$, which is
$`g\,a_m (E_i - V_m)`$: Hodgkin and Huxley's linear channel, with $`E_i`$ the Nernst potential at the
inside's concentration there. The kit's `hodgkinHuxley({ T, area })` writes their links so:

```js nocheck
const hh = hodgkinHuxley({ T, area: 2 / 238e-6 }); // a squid giant axon, radius 238 µm
ports: [{
  name: 'membrane', region: 'axon', V: 0, area: 2 / 238e-6,
  capacitance: { C: 0.01, zeroCharge: 0 },
  bath: { c: { 'Na+': 367, 'K+': 18.6, 'Cl-': 459 } }, // sea water, as far as the channels go
  gates: hh.gates,
  species: hh.species, // Na⁺ (G × m³h), K⁺ (G × n⁴), and the leak (Cl⁻), G = g × area
}]
```

With the axoplasm strictly neutral and the outside held (no resistance outside), the ions'
drift along x and the membrane's charging are cable theory, without its being put in: a squid
axon's action potential propagates at Hodgkin and Huxley's speed, 18.73 m/s for their 18.5 °C
axon (they computed 18.8 by hand), the cable equation's own to 2e-3
([`cable`](../test/cable.test.js) test). A myelinated axon is regions of one axoplasm, each node
of Ranvier and each internode its own, with a port each: the nodes' gated, the internodes' a small
capacitance (myelin) and nothing else. The grid refines at every node's faces, and the spike jumps
node to node, arriving at each within 1 µs of a compartmental cable's prediction (the same test).
Each solution's `ports[k]` reports the gates at each node of the window, `gates[name][j]`,
and the voltage they follow, `Vm[j]`. A gated window can't overlap another port's gated or
surfaced one.

### What a solution reports

Each solution's `ports[k]`, as the port has them:

| Field | Contents |
|---|---|
| `name`, `V` | the port's name and voltage |
| `flux[name]`, `current` | what it brings into the device: each species (mol/(m²·s)) and the current (A/m², the charging current included), totals through a [cross-section](#geometry) |
| `x[j]`, `area[j]` | with reactions or a capacitance: the window's node j (m) and the electrode's area per volume there (m²/m³) |
| `rates[q][j]` | reaction q's forward rate at node j, mol/(m²·s) of electrode (reaction first, then node) |
| `coverage[name][j]`, `bare[j]` | with a surface: each species' coverage θ, and the bare fraction θ₀ to full precision |
| `sigma[j]`, `charge` | with a capacitance: σ per area of electrode (C/m²), and the total the port's side holds |
| `gates[name][j]`, `Vm[j]` | with gates: each gate's fraction open, and the voltage across the capacitance that they follow (V) |

The arrays are typed (`Float64Array`), so `Array.from` them before mapping to anything but
numbers. The kit's [`polarization()`](kit.md#polarization-curves-an-evans-diagram) gives each reaction's rate at a spot of the
electrode with its level moved, an Evans diagram.

## Terminals

The two contacts and every port are the device's terminals, named `left`, `right` and by each
port's `name`. Each is driven on its own:

```js nocheck
contacts: {
  left: { V: 0, ... },               // held at a voltage
  right: { I: -2, ... },             // driven by a current into the device, A/m² (I: 0 is open circuit)
},
ports: [{ name: 'ref', I: 0, ... },  // a reference electrode: no current, its voltage read off
        { name: 'wire', V: 0, R: 1e-3, ... }], // a source V behind a series resistance R (Ω·m², or Ω with a geometry)
```

- **`V`**: the terminal's voltage, the shift of its outside phase's ladder, as the voltage of
  its `terminal` species (a gate's is its own metal's electrons). Neither `V` nor `I` given
  means `V: 0`.
- **`I`**: the current into the device through this terminal (conduction plus displacement);
  the voltage floats and is solved for.
- **`V` and `R`**: a source behind a resistance, $`I = (V - V_{\mathrm{terminal}})/R`$.

At least one terminal must be held at a voltage, or the device's overall level floats, and if
one is driven by a current, a held one must pass current (a closed end held at a voltage fixes
nothing). In
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

In a transient, a diffusion layer grows from the electrode as $`\sqrt{Dt}`$, through cells up to
`hmax` wide, and nothing warns while they're coarse: give `hmax` a small fraction of
$`\sqrt{Dt}`$ at the earliest time that matters. In a chronopotentiometry run, cells half
$`\sqrt{D\tau}`$ wide put the surface concentration 0.5% off, a fifteenth 0.06% (the error goes
as $`h^2`$).

## Geometry

A device is planar unless it says otherwise: a slab whose cross-section $`A`$ is 1 m², so every
current, flux, amount and charge reads per m². With `geometry`, the cross-section varies along x,
and the conservation laws hold for totals through it: a box holds $`\int A\,dx`$, and what crosses
a face is the flux there times $`A`$ there.

```js nocheck
geometry: { type: 'spherical', r0: 10e-6 },   // shells, r = r0 + x: a microelectrode of radius r0
geometry: { type: 'cylindrical', r0: 0 },     // a wire or fibre from its axis (per metre of length)
geometry: { area: { x: [0, 1e-3], values: [3e-6, 1e-7] } }, // A (m²) against x: a film thinning, a pit
```

This is how radial problems become 1D, with the right $`1/r^2`$ (or $`1/r`$) in every divergence and in
Gauss's law, without writing them out. The finite volumes integrate over shells, so conservation
is exact as in a slab. Each segment's transport is weighted by its length over $`\int dx/A`$, which
makes steady diffusion between two nodes exact for any $`A(x)`$ (the $`1/r`$ and $`\ln r`$ profiles
on any grid). Where $`A`$ vanishes at an end (a sphere's centre, `r0: 0`), the segment takes $`A`$ at
its middle, and nothing passes the end. Only an end may have $`A = 0`$. Flow (`velocity`) is for
planar devices.

With a cross-section, the solution's currents, contact and port fluxes, terminal currents and
`charge` are totals (A, mol/s; per metre for cylinders), and so are `I` drives, a series `R`
and the impedance (Ω, not Ω·m²). Face laws,
rates, `interfaces[f].N` and `D`, and concentrations stay per area or per volume. Validated against
steady diffusion to a sphere and a cylinder (exact), Cottrell's transient with the spherical term,
uptake by a sphere filling from its surface, and Debye–Hückel screening around a charged sphere
(test/geometry.test.js).

## Using a device

```js nocheck
const dev = new Device(def);
const sol = dev.solve();                      // steady state (or equilibrium)
dev.set({ contacts: { right: { V: 0.3 } } }); // deep-merged change; the state is kept as a warm start
dev.set({ ports: { gate: { V: 0.2 } } });     // ports patched by name (an array replaces them all)
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
  bias, it ramps the contacts' voltages from where they were last solved (the left one first,
  where it moved), or from level terminals on a cold start; with generation reactions (below),
  it ramps their rates up from nearly nothing. A terminal driven by a current (open circuit,
  say) is held at a voltage instead, the voltage marched until the current crosses its target,
  and floated from there. If 20 V either way brings no crossing (a current beyond the limiting
  one, or a device that only stores what comes in, which has no steady state at a current), the
  solve fails with a warning giving the range of currents the held voltages passed. A solve that
  fails reports the state it failed in, then leaves the device as it was before it.
- `step(dt, { method })` advances the transient by dt seconds, halving internally where Newton
  needs it. `method` is `'be'` (backward Euler, the default) or `'bdf2'`.
- `advance(tEnd, opts)` integrates adaptively to `tEnd` with variable-step BDF2, controlling
  the local error per step to `tol` (default 1e-3) in thermal units of every potential ($`\phi`$ and
  each $`\bar\mu/RT`$): roughly 0.1% in concentrations per step. That's local: over a long decay
  the errors add up, and a small current at the end of one can be tens of percent off, so tighten
  `tol` (1e-5, say) before fitting a time constant to a tail. It lands exactly on `tEnd`. Options: `tol`,
  `dt0` (first step), `dtMax`, `budgetMs` (return after this much wall time, with
  `done: false`), `maxSteps`, `method`. The step size carries over between calls, so an
  animation can call `advance(tNext, { budgetMs })` once per frame. The solution adds `done`,
  `rejected`, and a `trace` of terminal current and voltage after every accepted step (the
  current is the right contact's, toward +x, as `sol.current`; a port's is in `sol.terminals`
  after each call; the voltage is the right contact's minus the left's: what a meter between the terminals reads,
  which at a bath is its reference species' level, not $`\phi`$). With
  `probes: [{ x, species, quantity, region }]` the trace also reads inside the device,
  as `trace.probes[k]` beside `trace.t`: a species' concentration (`quantity: 'c'`, the default,
  mol/m³) or species voltage (`'V'`), or $`\phi`$ itself (`'phi'`, no species), linearly between the
  nodes around `x` (at a face, the side where the species is; where both sides have it, `region`
  says which). That's what a detector at `x` sees.
  At a face, `{ interface: f, species }` reads a species' flux through it (mol/(m²·s), toward +x,
  as `interfaces[f].N`), and `{ interface: f, gate }` a gate's fraction open. Probes read at accepted steps, which grow
  as a transient slows, so give `dtMax` to resolve a signal in time.
- `impedance(frequencies, { terminal, profiles })` solves the steady state, then linearises about
  it: $`Z(f) = \delta V/\delta I`$ in Ω·m² at one terminal (`'right'` by default), with I into the
  device. A held terminal's voltage is perturbed, or a driven one's current; the other terminals
  keep their drives (held ones at AC ground, driven ones open). A terminal behind a resistance can't
  be the one measured: the resistance belongs to the external circuit. With `profiles: true`, each
  frequency also returns complex profiles of $`\delta\phi`$, $`\delta\bar\mu`$ and $`\delta c`$ per unit
  excitation.
- `set(patch)` merges plain objects deeply (arrays are replaced), except that a contact given `V`
  or `I` drops the other: `{ contacts: { right: { I: 0 } } }` switches it to open circuit. A patch that changes only the
  terminals' drives (`V`, `I`, `R`) updates them in place, cheaply, keeping everything else.
  Otherwise the device is rebuilt: the current state carries over while the grid and species
  are unchanged, or restarts from the regions' `c0`. A carried state keeps what it holds: a
  stretch that was closed stays at its amount, and one that a change closes (or splits) keeps
  what it holds at that moment. So a new `c0` doesn't apply to a carried state; to start over
  from `c0`, make a new `Device`. Either way the time stepping restarts its
  order, as after any discontinuity, but a carried state goes on at the step size an adaptive
  transient had reached, so a page may change a rate constant or a diffusivity every frame.

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
| `flux[name]` | flux toward +x across each segment, between nodes `g` and `g + 1` (a face's across its pair of nodes), mol/(m²·s) (mol/s through a cross-section): the fluxes the balances add up |
| `D[name]` | diffusivity at each node, m²/s (`NaN` where absent, or in a metal) |
| `current`, `terminalVoltage` | current toward +x through the device (A/m², or A through a [cross-section](#geometry)) and $`V_{\mathrm{right}} - V_{\mathrm{left}}`$ |
| `terminals[name]` | `{ V, current }` for each terminal (contacts and ports), current into the device |
| `contacts.left/right` | `{ V, flux: {name}, D, current, links: {name} }` at each contact; `links` names how each species that crosses it does (`'equilibrium'`, `'velocity'`, …) |
| `gates.left/right` | charge on a gate or Stern plate, where the contact is capacitive |
| `ports[k]` | `{ name, V, flux: {name}, current }`: what each internal port brings into the device; an electrode's rates, coverages and charge [as well](#what-a-solution-reports) |
| `interfaces[f]` | `{ left, right, dipole, sheetCharge, D, N: {name}, rates, links: {name} }`: the names of the regions the face joins, what crosses it by its links, each face reaction's rate (mol/(m²·s)), and how each species that crosses it does; for those crossing by permeability, `oneWay: { name: [toward +x, toward −x] }`, the one-way fluxes whose difference is `N` (Ussing's unidirectional fluxes) |
| `bulkReactions[k]` | `{ rate, forward, nu, regions, total }`: each bulk reaction's net forward rate at every node (mol/(m³·s), `NaN` where it doesn't run), its one-way forward rate (the backward is `forward − rate`), what one forward reaction makes of each species (`nu`, negative for what it consumes), and integrated over each region and the device (mol/(m²·s)), as the balances count it: a charge-balance check is J = F(generation − recombination), which the kit's `check()` does for every species |
| `charge` | total charge in the device, C/m² (C through a cross-section) |
| `conservation` | per species stretch: amount, reference, intake through contacts, drift |
| `warnings` | e.g. unresolved double layers, conventions a statistics model relies on, and for a failed solve, where the system is nearly singular |
| `converged`, `iterations`, `steps`, `substeps`, `history`, `time`, `T` | solver bookkeeping, and the temperature |
| `steady` | whether this is a converged steady state (from `solve()`), not a transient's: what `check()` in the kit decides its checks by |
| `done`, `rejected`, `trace` | from `advance()`: whether `tEnd` was reached; rejected steps; `{ t, current, voltage }` per accepted step |

A stretch is a run of regions in which a species is present and connected. Its `drift`
compares its amount with the reference amount plus everything that came in through the
contacts. It's `NaN` where a reaction also makes or consumes the species.
