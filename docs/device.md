# Device reference

A device is defined by one plain object (serialisable, so it can be posted to a Web Worker),
passed to `new Device(def)`. Construction validates everything and throws a `DeviceError`
naming the offending path (`contacts.left.species.Na+.offset: …`). Physically meaningful
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
| `circuit` | how the terminals are driven (default: each contact at its own voltage) |
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

Per species present: `D` (m²/s, the coefficient in the flux N = −(D c/RT)∇μ̄), `mu0` (standard
chemical potential, J/mol) and `cRef` (the concentration that `mu0` refers to; falls back to
the species' `cRef`, and must exist one way or the other). Species not listed are absent.

In a material, only *neutral combinations* of standard potentials mean anything physically.
Shifting every `mu0` by `z_i F s` just moves the material's φ, its bookkeeping anchor. See
[conventions](conventions.md).

`epsr: 0` makes a material **strictly neutral**: Poisson is replaced by local neutrality, and φ
there is only a bookkeeping multiplier. Use it for metals and for macroscopic systems whose
double layers you don't want to resolve.

## Regions

```js nocheck
regions: [
  { name: 'n', material: 'Si', length: 2e-6, fixedCharge: 16.6 * FARADAY }, // C/m³
  { material: 'water', length: 1e-6, c0: { 'Na+': 10, 'Cl-': 10 } },
]
```

- `fixedCharge`: immobile charge density (doping, ionomer), C/m³. Default 0.
- `c0`: initial concentrations, mol/m³. Species connected to a contact start from that
  contact's level and don't need it. Any other species (blocked everywhere, or only made and
  consumed by reactions) does: its `c0` fixes the amount it conserves.
- `grid`: per-region override of the grid options (below).

## Interfaces

One entry per face, `interfaces[f]` sitting between `regions[f]` and `regions[f + 1]`. All
fields are optional where the defaults apply.

```js nocheck
interfaces: [
  {
    step: { species: 'e-', value: units.eV(0.25) }, // alignment: or { dipole: volts }
    phi: 'dipole', // 'dipole' | 'neutral' | { type: 'capacitive', C }
    sheetCharge: 0, // C/m²
    species: { 'Cl-': 'blocked', 'Li+': { type: 'conductance', G: 50 } },
    reactions: [{ transfer: { 'Na+': 1 }, k0: 1e-3, alpha: 0.5 }],
  },
]
```

**Alignment** (exactly one, for a face between *different* materials under a `dipole` or
`capacitive` law):
- `step: { species, value }`: the step in that charged species' standard level across the
  face, right minus left, `(μ°_R + zFφ_R) − (μ°_L + zFφ_L)`, in J/mol. For electrons that's
  the conduction-band offset.
- `dipole`: the φ jump, right minus left, in each material's own anchoring, in volts.

There is no default: omitting it is an error. A face between regions of the same material
defaults to no dipole.

**Electrostatic law** `phi`:
- `'dipole'`: φ jumps by the alignment. The default, and exact when the grid resolves the
  double layers on both sides.
- `'neutral'`: no charge at the face and a free jump, set by neutrality on each side (Donnan).
  The alignment drops out and must not be given. This is the default between two ε = 0
  materials.
- `{ type: 'capacitive', C }`: a Helmholtz layer, D = −C·(Δφ − dipole), C in F/m².

**Species laws** `species` (default: local equilibrium, μ̄ continuous, where the species is
present on both sides; blocked otherwise): `'equilibrium'`, `'blocked'`, or
`{ type: 'conductance', G }` (J = G·(V_L − V_R), G in S/m², charged species).

**Transfer kinetics** `reactions`: listed species cross left to right (forward) with a
Butler–Volmer rate per area,
`r = k0 Π[(c_L/c_ref,L)^{ν(1−α)} (c_R/c_ref,R)^{να}] (e^{αa} − e^{−(1−α)a})`,
with a = Σν(μ̄_L − μ̄_R)/RT, k0 in mol/(m²·s), and α (default 0.5) between 0 and 1.

## Bulk reactions

```js nocheck
bulkReactions: [
  { reactants: { 'e-': 1, 'h+': 1 }, kf: { Si: 2e8 } }, // e⁻ + h⁺ ⇌ ∅
  { reactants: { 'H+': 1, 'OH-': 1 }, products: { H2O: 1 }, fixed: { H2O: -237.13e3 }, kf: { water: 1.4e8 } },
]
```

Participants that aren't species are fixed-activity neutrals, given by their μ in `fixed`
(J/mol). Charge must balance. The rate is mass action,
`r = k_f Π c_R^ν (1 − e^{−A/RT})` with A the affinity. That's `k_f Π c_R − k_b Π c_P`, with
`k_b` fixed by the standard potentials, so equilibrium is exactly A = 0. `kf` maps material
names to forward rate constants (units making r mol/(m³·s)), and the reaction runs only in
those materials.

## Contacts

```js nocheck
contacts: {
  left: {
    V: 0, // terminal voltage, V
    terminal: 'e-', // the species whose voltage V is
    species: { 'e-': 'fixed', 'h+': { type: 'fixed', offset: 0 }, 'Cl-': 'blocked' },
    phi: 'neutral', // 'neutral' | 'free' | { type: 'capacitive', C, zeroCharge }
    reactions: [ /* electrode reactions, below */ ],
  },
  right: { bath: { c: { 'Na+': 10, 'Cl-': 10 }, reference: 'Cl-' } },
}
```

**Species links** (default `'blocked'`):
- `'fixed'`: the species is held at `V_i = V + offset` (charged) or `μ̄ = mu` (neutral, J/mol).
  The offset defaults to 0 only for the terminal species. Any other fixed species needs one,
  because it's an interface property (e.g. E° for an ion at a reversible electrode). At a metal,
  holes and electrons share the metal's voltage: `offset: 0`.
- `{ type: 'conductance', G, offset }`: ohmic exchange with an outside level at
  `V + offset`, J = G·(V_out − V_i).
- `'blocked'`: no flux.

**`bath`** (instead of `species` and `phi`): a neutral composition `c` and a charged
`reference` species (as for a reversible reference electrode, e.g. Cl⁻ for Ag/AgCl), which is
the terminal. Every bath species is fixed at the level its composition implies, and φ is
neutral. An optional `offset` places the reference species relative to V.

**`phi`**, required whenever any species or reaction connects at the contact:
- `'neutral'`: the contact node is locally neutral (ohmic contacts, baths). The metal takes
  whatever surface charge that needs.
- `{ type: 'capacitive', C, zeroCharge }`: a gate or Stern layer to a conductor at the
  terminal voltage V. The displacement into the device is `C·((V − zeroCharge) − φ_edge)`, so
  `zeroCharge` is the value of V − φ_edge at which the interface carries no charge: the
  flat-band voltage of a gate, or the potential of zero charge (pzc) of an electrode. It's an
  interface property, never computed from work functions.
- `'free'`: zero field.

**Electrode `reactions`**, written as reduction when `electrons > 0`:
`Σν_R R + n e⁻(metal) ⇌ Σν_P P`. The metal's electrons sit at μ̄_e = −F·V, and
non-species participants are fixed-activity (`fixed`, μ in J/mol). The rate is the
standard-rate-constant form
`r = k0 Π_R (c/c_ref)^{ν(1−α)} Π_P (c/c_ref)^{να} (e^{αa} − e^{−(1−α)a})`,
k0 in mol/(m²·s). The reaction sits at the contact node, just behind any Stern layer, so
Frumkin effects arise by themselves.

## Circuit

```js nocheck
circuit: { mode: 'voltage' }        // default: each contact at its own V
circuit: { mode: 'current', I: 2 }   // A/m² toward +x through the right terminal; I = 0 is open circuit
circuit: { mode: 'load', R: 1e-3, V: 0.1 } // I = (V_right − V_left − V)/R, R in Ω·m²
```

The left terminal is the reference. In current and load modes the right terminal's voltage is
solved for, which needs a current path there: a fixed charged terminal species, an electrode
reaction or a conductance link. A floating kinetic electrode also needs a capacitive or free φ
link.

## Grid

`{ hmin, hmax, ratio = 1.2, minCells = 8 }`, per device, overridable per region. Cells grow
geometrically from `hmin` at both ends of each region up to `hmax`, scaled to fit exactly. Each
region boundary becomes a pair of nodes at the same x, one per side. The grid is never refined
automatically. If a double layer that the model resolves is coarser than the local Debye
length, the solution's `warnings` say so.

## Using a device

```js nocheck
const dev = new Device(def);
const sol = dev.solve();                      // steady state (or equilibrium)
dev.set({ contacts: { right: { V: 0.3 } } }); // deep-merged change; the state is kept as a warm start
const sol2 = dev.solve();
const tr = dev.step(1e-6);                    // backward-Euler transient step (auto-subdivided if needed)
const now = dev.solution();                   // snapshot of the current state
```

- `solve()` finds the steady state from the current state without advancing time. If every
  species is fed by a contact, it solves the steady equations directly. Otherwise conserved
  amounts (blocked species, reactive moieties) are kept exactly.
- `step(dt)` advances the transient by dt seconds, halving internally where Newton needs it.
- `set(patch)` merges plain objects deeply (arrays are replaced). The current state carries
  over while the grid and species are unchanged; otherwise it restarts from the regions' `c0`.

## Solutions

| Field | Contents |
|---|---|
| `x`, `region` | node positions and region index; each interface position appears twice (one node per side) |
| `phi` | bookkeeping φ, V (`NaN` where undefined) |
| `c[name]`, `mu[name]`, `muStd[name]` | concentration, μ̄, standard level μ° + zFφ (`NaN` where absent) |
| `V[name]`, `Vstd[name]` | species voltage μ̄/(zF) and standard level as a voltage (charged species) |
| `current`, `terminalVoltage` | terminal current toward +x (A/m²) and V_right − V_left |
| `contacts.left/right` | `{ V, flux: {name}, D, current }` at each contact |
| `gates.left/right` | charge on a gate or Stern plate, where the contact is capacitive |
| `interfaces[f]` | `{ dipole, sheetCharge, D, N: {name} }`: what crosses each face |
| `charge` | total charge in the device, C/m² |
| `conservation` | per species stretch: amount, reference, intake through contacts, drift |
| `warnings` | e.g. unresolved double layers |
| `converged`, `iterations`, `steps`, `substeps`, `history`, `time` | solver bookkeeping |

A stretch is a run of regions in which a species is present and connected. Its `drift`
compares its amount with the reference amount plus everything that came in through the
contacts. It's `NaN` where a reaction also makes or consumes the species.
