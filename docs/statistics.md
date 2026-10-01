# Statistics

A material's **statistics** say how its concentrations follow from its species' potentials.
Every species *i* has a reduced chemical potential

```
ζ_i = (μ̄_i − z_i F φ − μ°_i) / RT
```

and a statistics model gives the concentrations c(ζ) of the species it covers. Species that
aren't in any model are ideal, c = c_ref·e^ζ, which is the default and the right choice for
dilute species. Non-ideal models describe degenerate carriers, crowding, filling of a finite
set of sites, activity coefficients, and intercalation hosts.

Three rules hold for every model:

- **The dilute limit is ideal.** Every model tends to c_ref·e^ζ as c → 0, so `mu0` and `cRef`
  keep one meaning (the dilute, Henry's-law reference state) whatever the statistics. Changing
  a material's statistics changes only its non-dilute behaviour.
- **The Jacobian is a chemical capacitance.** K = ∂c/∂ζ is the Hessian of a convex potential
  (the grand potential density P(ζ), with c = ∂P/∂ζ), so it's symmetric and positive definite.
  RT·K is the chemical capacitance matrix per volume. Models whose free energy isn't convex
  (a miscibility gap) are refused, because they need phase separation (Cahn–Hilliard, on the
  [roadmap](../ROADMAP.md)).
- **D keeps its meaning.** Fluxes are always N = −(D c/RT)∇μ̄, with D the dilute-limit
  coefficient (constant mobility). The chemical diffusivity you'd measure is D times the
  thermodynamic factor ∂ln a/∂ln c, the generalised Einstein relation. See
  [conventions](conventions.md#transport).

## Declaring statistics

```js nocheck
materials: {
  Si: {
    epsr: 11.7,
    species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } },
    statistics: [{ type: 'fermi-dirac', species: ['e-', 'h+'] }],
  },
}
```

`statistics` is a list of models, each naming the species it covers. A species belongs to at
most one model.

| Model | Species | c(ζ) | Parameters |
|---|---|---|---|
| `fermi-dirac` | any, each independently | c_ref 𝓕_{1/2}(ζ) (3D band); c_ref ln(1 + e^ζ) with `order: 0` (2D) | `order`: 1/2 (default) or 0 |
| `lattice` | any, sharing sites | c_a = c_max w_a / (1 + Σ w), w_a = (c_ref,a/c_max) e^{ζ_a} | `cMax`, site density, mol/m³ |
| `redlich-kister` | one | on its own sites: ζ = ln(x/(1−x)) + ln(c_max/c_ref) + (g(x) − g(0))/RT, x = c/c_max | `cMax`; `A`, J/mol |
| `debye-huckel` | charged | c_ref e^ζ / γ, ln γ_i = −z_i² ℓ_B κ / (2(1 + κa)) | `epsr` (solvent, for the theory); `a`, m |
| `insertion` | [ion, carrier], ε = 0 host | occupancy from the neutral combination (below) | `cMax`; `A` or `ocv` |
| `custom` | any | your function | `evaluate(zeta) → { c, dcdzeta }` |

### Fermi–Dirac

For electrons and holes, c_ref is the effective density of states and μ° the band edge, as for
ideal statistics, and c = N_c 𝓕_{1/2}(ζ) with the normalised Fermi–Dirac integral
(𝓕_j → e^ζ for ζ → −∞). driftlet evaluates 𝓕_{1/2} and 𝓕_{−1/2} (its derivative) to ~1e-14
relative, so equilibrium carrier densities are exact. `order: 0` is a two-dimensional density
of states, and it's exact in closed form.

### Lattice gas

The listed species share c_max sites per volume and each occupies one, so no concentration can
exceed c_max. With both ions of an electrolyte on one lattice, this is the Bikerman /
Kilic–Bazant–Ajdari model of a crowded double layer. At large potentials the counter-ion
saturates at c_max, and the diffuse charge grows as √ψ instead of exponentially. With one
species it's Langmuir filling of sites. The dilute reference c_ref is kept. To use the
half-filling convention of site models (e^ζ = x/(1−x)), set `cRef` equal to `cMax`.

### Redlich–Kister

One species on its own c_max sites, with the excess Gibbs energy of mixing per site
x(1−x) Σ A_k (1−2x)^k. g(x) in the table is its x-derivative, the excess chemical potential.
A single `A: [Ω]` is the regular solution, μ = μ° + RT ln(x/(1−x)) + Ω(1−2x) + const. It's
convex only for Ω < 2RT, which is checked at construction.

The table form is referred to the dilute limit. The half-filling convention of the battery
literature, μ = μ°_½ + RT ln(x/(1−x)) + Ω(1−2x), corresponds to `cRef = cMax` and
`mu0 = μ°_½ + A₀` (the excess at x = 0, Σ A_k in general).

### Debye–Hückel

Extended Debye–Hückel activity coefficients for dilute electrolytes: ln γ_i = −z_i² ℓ_B κ /
(2(1 + κa)), with the Bjerrum length ℓ_B = e²/(4πε kT) and κ² = F² Σ z² c / (ε RT) from the
solvent permittivity `epsr` (which defaults to the material's, if that's non-zero) and a
distance of closest approach `a` (default 0, the limiting law). It stops being convex at high
ionic strength, and the solver then reports a failure.

Single-ion activity coefficients are a convention: only neutral combinations (the mean
activity γ±) are measurable. In a strictly neutral (ε = 0) material only those combinations
matter, so that's where this model belongs. With ε > 0, the resolved double layers would
depend on the convention, so driftlet warns.

### Insertion hosts: neutral-combination statistics

In an intercalation host, an inserted ion enters together with its compensating electronic
carrier, e.g. Li⁺ + e⁻. Only their neutral combination has a well-defined chemical potential,
μ_Li = μ̄_Li⁺ + μ̄_e⁻, and the host's composition depends on nothing else. That's what an OCV
curve measures, and the mutual chemical capacitance of the pair is the slope of that curve.

The `insertion` model describes exactly this, in a strictly neutral host (`epsr: 0`) whose
only charged species are the pair `[ion, carrier]`. For an ion of charge z and a carrier of
charge −z/ν:

```
ζ_comb = ζ_ion + ν ζ_carrier             (φ cancels)
c_ion = s(ζ_comb),   c_carrier = ν s + background
```

The occupancy s = c_max x comes either from a Redlich–Kister isotherm (`A`, as above, referred
to the ion's dilute reference) or from a table:

```js nocheck
statistics: [{
  type: 'insertion',
  species: ['Li+', 'e-'],
  cMax: 30000, // sites, mol/m³
  ocv: { x: [...], E: [...], muRef: 0 }, // E against a reference where μ_Li = muRef
}]
```

The table gives the open-circuit voltage E(x) of the host against a reference electrode in
which the combination has chemical potential `muRef` (J/mol), so that
μ_ion + ν μ_carrier = muRef − z F E(x). Against Li metal, on the usual table convention,
`muRef` is 0. x must increase and E must strictly decrease. A flat plateau means phase
separation, which isn't supported yet. The table is read as ζ = ln(x/(1−x)) + r(x): the
residual r, the host's non-ideality, is a shape-preserving cubic between points and continues
linearly beyond them. Ideal and regular-solution curves are therefore reproduced exactly, and
x stays within (0, 1) at any potential. The interpolated curve is checked for monotonicity at
construction.

Because every charged species sits in a neutral combination, **φ is undefined** in such a host
and is reported as `NaN`. Its faces must be `neutral` (the default next to such a host), and
its contacts `bulk` or `neutral`. The host's `fixedCharge` is balanced by background carriers,
`background = −fixedCharge/(z_carrier F)`, which must be non-negative. This is a jellium
description of the host's own electrons, which keep it conducting when it's empty. With no
fixed charge, the carriers are only those that came in with the ions. The carrier's `cRef`
plays no part (only the combination's μ° and the ion's `cRef` enter), but it must still be
given.

An insertion host belongs with [metal regions](device.md#materials) more than with the other
statistics here. Both are conductors with no φ of their own:
- **A metal** is the one-carrier case: its Fermi level is its only unknown.
- **An insertion host** is the two-carrier case: only the neutral combination's potential sets
  the composition, and the split between ion and carrier is free.

A general mixed conductor with more mobile species needs fully populated mutual capacitance and
conductance matrices, a different model again (see the [roadmap](../ROADMAP.md)).

### Custom

```js nocheck
statistics: [{
  type: 'custom',
  species: ['A', 'B'],
  evaluate: (zeta) => ({ c: [cA, cB], dcdzeta: [[dA_dA, dA_dB], [dB_dA, dB_dB]] }),
}]
```

`evaluate` gets the species' ζ (in the listed order) and returns their concentrations and the
Jacobian ∂c/∂ζ. The Jacobian must be symmetric and positive definite, and that's checked at
ζ = 0 when the device is built. A device with a function in it isn't plain data, so it can't be
posted to a Web Worker. Build it inside the worker instead.

## Transport with non-ideal statistics

The flux law is N = −(D c/RT)∇μ̄ for every model. Writing ln c = ζ − ex with the excess
ex = ζ − ln(c/c_ref) (zero for ideal statistics), the excess acts like an extra potential. The
Scharfetter–Gummel flux keeps its form with Δ = zΔφ̂ + Δex (see [numerics](numerics.md#fluxes)),
so it's still exactly zero at equilibrium. It's second-order accurate in the cell size (and
exact for a single diffusing species on a lattice).

Reaction rates keep their concentration-based prefactors (mass action in c, and c/c_ref in the
Butler–Volmer forms), while the driving force is always the affinity from μ̄. So equilibrium
is exact for any statistics, and only the kinetics away from it depend on that choice of
prefactor.

Cross-species transport coefficients (an Onsager matrix) aren't included: each species moves
down its own μ̄. They should arrive together with cross chemical capacitances, as two halves of
the same network ([roadmap](../ROADMAP.md)). Note that K itself can already be non-diagonal,
e.g. on a shared lattice.

## Example: a lithium host against its OCV

```js
import { Device, FARADAY, GAS_CONSTANT } from 'driftlet';

const RT = GAS_CONSTANT * 298.15;
const x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];
const E = x.map((v) => 0.4 - (RT / FARADAY) * Math.log(v / (1 - v))); // an ideal-solution OCV, V vs Li

const dev = new Device({
  species: [
    { name: 'Li+', z: 1 },
    { name: 'e-', z: -1 },
  ],
  materials: {
    host: {
      epsr: 0,
      species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } },
      statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }],
    },
  },
  regions: [{ material: 'host', length: 1e-6 }],
  contacts: {
    // Li⁺ from a lithium reference (via an electrolyte) on the left; electrons from a current
    // collector on the right, whose voltage is the electrode potential vs Li.
    left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' },
    right: { V: 0.4, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
  },
});
for (const V of [0.5, 0.4, 0.3]) {
  dev.set({ contacts: { right: { V } } });
  const sol = dev.solve();
  console.log(`E = ${V} V: x = ${(sol.c['Li+'][5] / 30000).toFixed(4)}`); // 0.0200, 0.5000, 0.9800
}
```
