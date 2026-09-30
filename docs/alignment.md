# Band alignment: recipes from vacuum levels, and their caveats

driftlet never computes how levels line up across an interface. Each interface's alignment is
its own physical property, and you give it explicitly (see [conventions](conventions.md)). But
sometimes vacuum-level reasoning is the only estimate available. It's also a useful worked
reference. This page gives the standard recipes, translated into driftlet's inputs, each with
its hidden assumptions and known failures. Use them knowingly.

**What every vacuum-level rule assumes:**
1. Each material's levels can be referred to a vacuum level that runs continuously through the
   interface. There's no interface dipole beyond what the materials' separate surfaces imply.
2. Surface quantities (work function W, electron affinity χ) are bulk properties. They aren't:
   they're barrier heights of particular surfaces, and change with termination, adsorbates and
   preparation.
3. As a consequence, alignments add up transitively around any set of materials, and are
   uniform along each interface.

At covalent interfaces (common semiconductors, most metal–semiconductor contacts), measured
alignments often depart substantially from these rules, and Fermi-level pinning is common.
Ionic materials tend to follow them better. See the
[ESBD vacuum-level discussion](https://marklundeberg.com/esbd/vacuum/) for the reasoning.

**Before reaching for a rule,** ask whether you need the alignment at all. If only macroscopic,
bulk-neutral behaviour matters, a `phi: 'neutral'` interface makes the alignment irrelevant:
it drops out, and none is needed.

## Translating energies into driftlet

For electrons (z = −1), the standard level in J/mol, `μ°_e + z F φ = μ°_e − F φ`, *is* the
conduction-band energy (per mole). So a face's `step` alignment for `'e-'` is the
conduction-band offset, right minus left. In eV per particle, use `units.eV(ΔE_c)`.

A metal modelled as a strictly neutral region (`epsr: 0`), whose electron `cRef` equals its
fixed-background electron density, has c = c_ref everywhere. Its electron standard level then
coincides with its Fermi level.

## Semiconductor heterojunctions: Anderson's rule

With electron affinities χ_L, χ_R, the rule puts each conduction band at −χ below a common
vacuum level:

```
ΔE_c = E_c,R − E_c,L = χ_L − χ_R
```

```js nocheck
interfaces: [{ step: { species: 'e-', value: units.eV(chiL - chiR) } }]
```

The valence-band offset then follows from the two bulk gaps, already present in your holes'
μ°. Measured offsets frequently differ from this by tenths of an eV. When data exist (e.g.
photoemission offsets for your pair), use them instead: `step` takes the measured ΔE_c directly.

## Metal–semiconductor barriers: the Schottky–Mott rule

The rule gives the n-type barrier as φ_B,n = W_m − χ_s. With the metal as a neutral region as
above (metal on the left):

```js nocheck
interfaces: [{ step: { species: 'e-', value: units.eV(Wm - chiS) } }] // = the barrier φ_B,n
```

With μ̄_e continuous across the face (the default), the semiconductor's surface electron
density is then N_c·e^{−φ_B/kT}, whatever the bias. Measured barriers on Si and GaAs are
largely insensitive to W_m (strong pinning), so a measured φ_B is far preferable. It goes in the
same `step`.

To model the metal as a contact instead of a region, put `'e-'` in equilibrium at the contact
(offset 0) and pin φ with a `dipole` law:

```js nocheck
contacts: { left: { terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: { type: 'dipole', zeroCharge: phiB - mu0e / FARADAY } } }
```

Here `mu0e` is the semiconductor's own electron standard potential (J/mol). With `mu0e = 0`,
`zeroCharge` is simply the barrier. A capacitive law with finite C adds an interfacial layer.

## Gates: flat-band voltage

A gate at terminal voltage V couples through a capacitive link, with no charge when
`V − φ_edge = zeroCharge`. Vacuum alignment puts the semiconductor's conduction band at the
surface W_g − χ above the gate's Fermi level at flat band, which gives

```
zeroCharge = (W_g − χ)/e − μ°_e/F
```

(μ°_e is the semiconductor's electron standard potential, J/mol). This is the anchor-level form
of the textbook V_FB = (W_g − W_s)/e, which refers to the semiconductor's bulk Fermi level
instead. The same caveats apply, plus any fixed charge in a real insulator, which shifts flat
band further.

## Electrodes: potential of zero charge from the work function

With ion μ° taken from the usual SHE-based tables, the solution's bookkeeping φ is on the SHE
scale (see [conventions](conventions.md)). So a Stern link's `zeroCharge` is simply the
potential of zero charge against SHE. Trasatti's reading of the "absolute" SHE (≈4.44 V) gives
the estimate

```
zeroCharge ≈ W_m/e − 4.44 V
```

This neglects the metal–solvent interaction terms (solvent orientation and chemisorption at the
metal surface), which are often a few tenths of a volt. The 4.44 V is itself a surface
property, the work function of a hydrogen electrode, not a universal zero. Measured pzc values
are preferable.

## Ions between solvents: extrathermodynamic assumptions

Transfer energies of whole salts between solvents are measurable. Single-ion values are not, and
need an extrathermodynamic assumption. The TATB assumption, for example, takes
ΔG_tr(Ph₄As⁺) = ΔG_tr(Ph₄B⁻). To use such values, give each ion in solvent B the standard
potential `μ°_B = μ°_A + ΔG_tr(A→B)` and pin the face with `dipole: 0`. That combination *is*
the assumption: it states that the bookkeeping φ is continuous when single-ion energies are
split that way.

Different assumptions disagree by amounts that can matter at the interface. For bulk
partitioning of salts, a `neutral` face needs no single-ion values at all, since only the
measurable salt combinations enter.
