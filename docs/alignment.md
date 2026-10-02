# Band alignment from vacuum levels

driftlet never decides how levels line up across an interface. Each interface's alignment is
its own physical property, and you give it explicitly (see [conventions](conventions.md)). But
sometimes vacuum-level reasoning is the only estimate available, and it's a useful worked
reference. This page sets it out as one recipe with two ingredients, and lists what it
assumes.

**Before reaching for it,** ask whether you need an alignment at all. If only macroscopic,
bulk-neutral behaviour matters, a `phi: 'neutral'` interface makes the alignment drop out, and
none is needed.

## The recipe: an anchor and an offset per material

Every vacuum-level number in the literature describes **one material against vacuum**. It
picks a level inside the material (the **anchor**) and says how far above that level, in
electron energy, the vacuum level just outside the surface lies (the **offset**). In
driftlet's voltage units (V = −E/e for electrons):

```
V_vac = V_anchor − offset
```

| Quantity | Anchor | Offset |
|---|---|---|
| electron affinity χ | conduction band = the e⁻ standard level | χ |
| ionisation energy | valence band = the h⁺ standard level | χ + E_g |
| work function W | Fermi level | W |
| Trasatti's absolute SHE, ≈ 4.44 V | the SHE level of a solution | 4.44 V |
| surface potential χ_s | the inner potential φ | χ_s |

Two materials are then aligned by **assuming their vacuum levels coincide** across the
interface. That assumption is the whole heuristic: Anderson's rule, the Schottky–Mott rule, a
gate's flat-band voltage and the pzc estimate are all this one recipe applied to different
pairs. It can combine any pair, for example a semiconductor (electron affinity) against an
aqueous electrolyte (absolute SHE).

### The anchors in driftlet

An anchor has to be a level that driftlet knows in terms of a material's φ:

- **A charged species' standard level**, `V°_i = φ + μ°_i/(z_i F)`. For e⁻ and h⁺ these are
  the band edges, which are electronic levels. For an ion, V°_i is simply that ion's
  standard level: a fixed offset from φ, not by itself the level of any electrons.
- **The SHE level**, the electronic level of the standard hydrogen electrode in a solution:

  ```
  V°_e⁻(SHE) = φ + (μ°_H⁺ − ½ μ°_H₂) / F
  ```

  As an anchor this is `'H+'` (its standard level φ + μ°_H⁺/F), with the H₂ term folded into
  the offset: Trasatti's value becomes `{ anchor: 'H+', offset: 4.44 + μ°_H₂/(2F) }`. On the
  usual conventions (μ°_H₂ = 0 for the element in its standard state, and μ°_H⁺ = 0 so that
  E°(SHE) = 0), that's just offset 4.44 V. Trasatti's 4.44 V then works exactly like a
  semiconductor's electron affinity, with the SHE level in place of the conduction band. It
  needs H⁺ among the species. Other half-reactions could be used as anchors in principle, but
  each brings its own conventions (the neutral element's μ, E°(SHE) = 0), so they're not
  offered as shortcuts here.
- **The Fermi level** of a metal, the proper reference for a work function. It's taken by
  `vacuumZeroCharge`, for a conductor region's face or a contact (below). A semiconductor's Fermi
  level moves with doping and bias, so its work function isn't a material constant: use its
  electron affinity instead.
- **φ itself** (`'phi'`), with a surface potential as the offset. That's discouraged, since φ
  is bookkeeping: its value inside a material depends on how that material's μ° were
  anchored, and surface potentials are not measurable on their own. If you have settled that
  bookkeeping to your own satisfaction, it's available.

## In the API

The recipe lives in `driftlet/kit`, as helpers that read a device definition and return the
spec value to put in it. The spec itself only takes the result (a `dipole` or a `zeroCharge`),
so what a device assumes stays visible in its definition.

```js nocheck
import { vacuumDipole, vacuumZeroCharge } from 'driftlet/kit';

// GaAs | AlGaAs: electron affinities on the e⁻ standard level (the conduction band)
def.interfaces = [{
  dipole: vacuumDipole(def, { material: 'GaAs', anchor: 'e-', offset: 4.07 }, { material: 'AlGaAs', anchor: 'e-', offset: 3.8 }),
}];
```

`vacuumDipole` gives φ_R − φ_L = (V_vac − φ)_L − (V_vac − φ)_R, each side's V_vac − φ coming
from its anchor and offset (`vacuumLevel` gives one side's).

Where a conductor meets a material, at a conductor region's face or at a contact's capacitive
or pinned law, the conductor's side is its work function W from its Fermi level, and the
helper gives `zeroCharge` (V_F − φ_edge at zero charge):

```js nocheck
const zc = vacuumZeroCharge(def, 4.75, { material: 'Si', anchor: 'e-', offset: 4.05 }); // Au | n-Si
def.interfaces = [{ phi: { type: 'capacitive', C: 10 }, zeroCharge: zc }];        // Au as a region
Object.assign(def.contacts.left, { phi: { type: 'capacitive', C: 10 }, zeroCharge: zc }); // or as a contact
// zeroCharge = W + (V_vac − φ)_inside = W − χ − μ°_e/F
```

## Pairwise examples

**Semiconductor heterojunction (Anderson's rule).** Electron affinities on both sides,
`{ material, anchor: 'e-', offset: χ }`, give ΔE_c = χ_L − χ_R. The valence-band offset then
follows from the two gaps, already in your holes' μ°. Ionisation energies on `'h+'` give the same thing from
the other band.

**Metal | semiconductor (Schottky–Mott).** Work function against electron affinity gives the
n-type barrier φ_B = W − χ. With the metal as a contact, a pinned law
`phi: 'pinned', zeroCharge: vacuumZeroCharge(def, W, { material, anchor: 'e-', offset: χ })`
puts the semiconductor's surface electron density at N_c e^{−φ_B/kT}. With the metal as a
region, the face is capacitive with the same `zeroCharge`; a large C approaches the pinned
barrier.

**Gate | semiconductor (flat band).** The same with a capacitive law. Its `zeroCharge` is the
anchor-level form of the textbook V_FB = (W_g − W_s)/e, which refers to the semiconductor's
bulk Fermi level instead. Fixed charge in a real insulator shifts flat band further.

**Metal | electrolyte (potential of zero charge).** Work function against the absolute SHE
gives pzc ≈ W/e − 4.44 V on the SHE scale (on the usual conventions above):
`zeroCharge: vacuumZeroCharge(def, W, { material, anchor: 'H+', offset: 4.44 })`.

**Semiconductor | electrolyte.** Electron affinity against the absolute SHE places the band
edges on the SHE scale, E_c ≈ −(χ − 4.44) eV relative to SHE: the usual photoelectrochemistry
estimate: `vacuumDipole(def, { material: 'TiO2', anchor: 'e-', offset: χ }, { material: 'water', anchor: 'H+', offset: 4.44 })`.

**Solvent | solvent.** Each solvent's own absolute-SHE-type value (its Trasatti-type offset)
aligns them, exactly as Anderson's rule aligns two semiconductors.

## What the heuristic assumes, and where it fails

1. **The vacuum level runs continuously through the interface.** In reality, the interface
   forms its own dipole (charge transfer, bonding, solvent orientation, adsorbates) beyond
   what the two free surfaces imply.
2. **Surface quantities are properties of the material.** They aren't: W, χ and 4.44 V are
   barrier heights of particular surfaces, and change with termination, adsorbates and
   preparation. The 4.44 V is itself the work function of a hydrogen electrode, not a
   universal zero.
3. **Alignments add up transitively** around any set of materials, and are uniform along each
   interface.

At covalent interfaces (common semiconductors, most metal–semiconductor contacts), measured
alignments often depart substantially from these rules, and Fermi-level pinning is common.
Ionic materials tend to follow them better. At metal–solution interfaces, the neglected
metal–solvent terms (orientation, chemisorption) are often a few tenths of a volt.
Whenever a measured offset, barrier, flat-band voltage or pzc exists, use it instead:
`step`, `dipole` and `zeroCharge` take measured values directly.
