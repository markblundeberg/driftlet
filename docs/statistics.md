# Statistics

A material's **statistics** say how its concentrations follow from its species' potentials.
Every species *i* has a reduced chemical potential

```math
\zeta_i = (\bar\mu_i - z_i F \phi - \mu^\circ_i) / RT
```

and a statistics model gives the concentrations $c(\zeta)$ of the species it covers. Species
that aren't in any model are ideal, $c = c_{\mathrm{ref}} \cdot e^\zeta$, which is the default and
the right choice for dilute species. Non-ideal models describe degenerate carriers, crowding,
filling of a finite set of sites, activity coefficients, and intercalation hosts.

Three rules hold for every model:

- **The dilute limit is ideal.** Every model tends to $c_{\mathrm{ref}} \cdot e^\zeta$ as $c \to 0$,
  so `mu0` and `cRef` keep one meaning (the dilute, Henry's-law reference state) whatever the
  statistics. Changing a material's statistics changes only its non-dilute behaviour.
- **The Jacobian is a chemical capacitance.** $K = \partial c/\partial\zeta$ is the Hessian of a
  convex potential (the grand potential density $P(\zeta)$, with $c = \partial P/\partial\zeta$),
  so it's symmetric and positive definite. $RT \cdot K$ is the chemical capacitance matrix per
  volume. Models whose free energy isn't convex (a miscibility gap) are refused, because they
  need phase separation (Cahn–Hilliard, on the [roadmap](../ROADMAP.md)).
- **D keeps its meaning.** Fluxes are always $N = -(D c/RT)\nabla\bar\mu$, with $D$ the
  dilute-limit coefficient (constant mobility). The chemical diffusivity you'd measure is $D$
  times the thermodynamic factor $\partial\ln a/\partial\ln c$, the generalised Einstein
  relation. See [conventions](conventions.md#transport).

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

| Model | Species | $c(\zeta)$ | Parameters |
|---|---|---|---|
| `fermi-dirac` | any, each independently | $c_{\mathrm{ref}} \mathcal{F}_{1/2}(\zeta)$ (3D band); $c_{\mathrm{ref}} \ln(1 + e^\zeta)$ with `order: 0` (2D) | `order`: 1/2 (default) or 0 |
| `lattice` | any, sharing sites | $c_a = c_{\mathrm{max}} w_a / (1 + \sum w)$, $w_a = (c_{\mathrm{ref},a}/c_{\mathrm{max}}) e^{\zeta_a}$ | `cMax`, site density, mol/m³ |
| `redlich-kister` | one | on its own sites: $\zeta = \ln(x/(1-x)) + \ln(c_{\mathrm{max}}/c_{\mathrm{ref}}) + (g(x) - g(0))/RT$, $x = c/c_{\mathrm{max}}$ | `cMax`; `A`, J/mol |
| `debye-huckel` | charged | $c_{\mathrm{ref}} e^\zeta / \gamma$, $\ln\gamma_i = -z_i^2 \ell_B \kappa / (2(1 + \kappa a))$ | `epsr` (solvent, for the theory); `a`, m |
| `insertion` | [ion, carrier], $\varepsilon = 0$ host | occupancy from the neutral combination (below) | `cMax`; `A` or `ocv` |
| `custom` | any | your function | `evaluate(zeta) → { c, dcdzeta }` |

### Fermi–Dirac

For electrons and holes, $c_{\mathrm{ref}}$ is the effective density of states and $\mu^\circ$ the
band edge, as for ideal statistics, and $c = N_c \mathcal{F}_{1/2}(\zeta)$ with the normalised
Fermi–Dirac integral ($\mathcal{F}_j \to e^\zeta$ for $\zeta \to -\infty$). driftlet evaluates
$\mathcal{F}_{1/2}$ and $\mathcal{F}_{-1/2}$ (its derivative) to ~1e-14 relative, so equilibrium
carrier densities are exact. `order: 0` is a two-dimensional density of states, and it's exact
in closed form.

### Lattice gas

The listed species share $c_{\mathrm{max}}$ sites per volume and each occupies one, so no
concentration can exceed $c_{\mathrm{max}}$. With both ions of an electrolyte on one lattice, this
is the Bikerman / Kilic–Bazant–Ajdari model of a crowded double layer. At large potentials the
counter-ion saturates at $c_{\mathrm{max}}$, and the diffuse charge grows as $\sqrt{\psi}$ instead
of exponentially. With one species it's Langmuir filling of sites. The dilute reference
$c_{\mathrm{ref}}$ is kept. To use the half-filling convention of site models
($e^\zeta = x/(1-x)$), set `cRef` equal to `cMax`.

### Redlich–Kister

One species on its own $c_{\mathrm{max}}$ sites, with the excess Gibbs energy of mixing per site
$x(1-x) \sum A_k (1-2x)^k$. $g(x)$ in the table is its $x$-derivative, the excess chemical
potential. A single `A: [Ω]` is the regular solution,
$\mu = \mu^\circ + RT \ln(x/(1-x)) + \Omega(1-2x) + \mathrm{const}$. It's convex only for
$\Omega < 2RT$, which is checked at construction.

The table form is referred to the dilute limit. The half-filling convention of the battery
literature, $\mu = \mu^\circ_{1/2} + RT \ln(x/(1-x)) + \Omega(1-2x)$, corresponds to
`cRef = cMax` and `mu0` $= \mu^\circ_{1/2} + A_0$ (the excess at $x = 0$, $\sum A_k$ in general).

### Debye–Hückel

Extended Debye–Hückel activity coefficients for dilute electrolytes:
$\ln\gamma_i = -z_i^2 \ell_B \kappa / (2(1 + \kappa a))$, with the Bjerrum length
$\ell_B = e^2/(4\pi\varepsilon kT)$ and $\kappa^2 = F^2 \sum z^2 c / (\varepsilon RT)$ from the
solvent permittivity `epsr` (which defaults to the material's, if that's non-zero) and a
distance of closest approach `a` (default 0, the limiting law). It stops being convex at high
ionic strength, and the solver then reports a failure.

Single-ion activity coefficients are a convention: only neutral combinations (the mean
activity $\gamma_\pm$) are measurable. In a strictly neutral ($\varepsilon = 0$) material only
those combinations matter, so that's where this model belongs. With $\varepsilon > 0$, the
resolved double layers would depend on the convention, so driftlet warns.

### Insertion hosts: neutral-combination statistics

In an intercalation host, an inserted ion enters together with its compensating electronic
carrier, e.g. $\mathrm{Li}^+ + \mathrm{e}^-$. Only their neutral combination has a well-defined
chemical potential, $\mu_{\mathrm{Li}} = \bar\mu_{\mathrm{Li}^+} + \bar\mu_{\mathrm{e}^-}$, and the
host's composition depends on nothing else. That's what an OCV curve measures, and the mutual
chemical capacitance of the pair is the slope of that curve.

The `insertion` model describes exactly this, in a strictly neutral host (`epsr: 0`) whose
only charged species are the pair `[ion, carrier]`. For an ion of charge $z$ and a carrier of
charge $-z/\nu$:

```math
\begin{aligned}
\zeta_{\mathrm{comb}} &= \zeta_{\mathrm{ion}} + \nu \zeta_{\mathrm{carrier}} \quad (\phi\text{ cancels}) \\
c_{\mathrm{ion}} &= s(\zeta_{\mathrm{comb}}), \qquad c_{\mathrm{carrier}} = \nu s + \text{background}
\end{aligned}
```

The occupancy $s = c_{\mathrm{max}} x$ comes either from a Redlich–Kister isotherm (`A`, as above,
referred to the ion's dilute reference) or from a table:

```js nocheck
statistics: [{
  type: 'insertion',
  species: ['Li+', 'e-'],
  cMax: 30000, // sites, mol/m³
  ocv: { x: [...], E: [...], muRef: 0 }, // E against a reference where μ_Li = muRef
}]
```

The table gives the open-circuit voltage $E(x)$ of the host against a reference electrode in
which the combination has chemical potential `muRef` (J/mol), so that
$\mu_{\mathrm{ion}} + \nu \mu_{\mathrm{carrier}} = \mathtt{muRef} - z F E(x)$. Against Li metal,
on the usual table convention, `muRef` is 0. $x$ must increase and $E$ must strictly decrease.
A flat plateau means phase separation, which isn't supported yet. The table is read as
$\zeta = \ln(x/(1-x)) + r(x)$: the residual $r$, the host's non-ideality, is a shape-preserving
cubic between points and continues linearly beyond them. Ideal and regular-solution curves are
therefore reproduced exactly, and $x$ stays within $(0, 1)$ at any potential. The interpolated
curve is checked for monotonicity at construction.

Because every charged species sits in a neutral combination, **$\phi$ is undefined** in such a
host and is reported as `NaN`. Its faces must be `neutral` (the default next to such a host),
and its contacts `bulk` or `neutral`. The host's `fixedCharge` is balanced by background
carriers, $\mathtt{background} = -\mathtt{fixedCharge}/(z_{\mathrm{carrier}} F)$, which must be
non-negative. This is a jellium description of the host's own electrons, which keep it
conducting when it's empty. With no fixed charge, the carriers are only those that came in with
the ions. The carrier's `cRef` plays no part (only the combination's $\mu^\circ$ and the ion's
`cRef` enter), but it must still be given.

An insertion host belongs with [conductor regions](device.md#materials) more than with the other
statistics here. Both are conductors with no $\phi$ of their own:
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

`evaluate` gets the species' $\zeta$ (in the listed order) and returns their concentrations and
the Jacobian $\partial c/\partial\zeta$. The Jacobian must be symmetric and positive definite,
and that's checked at $\zeta = 0$ when the device is built. A device with a function in it isn't
plain data, so it can't be posted to a Web Worker. Build it inside the worker instead.

## Transport with non-ideal statistics

The flux law is $N = -(D c/RT)\nabla\bar\mu$ for every model. Writing $\ln c = \zeta - \mathrm{ex}$
with the excess $\mathrm{ex} = \zeta - \ln(c/c_{\mathrm{ref}})$ (zero for ideal statistics), the
excess acts like an extra potential. The Scharfetter–Gummel flux keeps its form with
$\Delta = z\Delta\hat\phi + \Delta\mathrm{ex}$ (see [numerics](numerics.md#fluxes)), so it's still
exactly zero at equilibrium. It's second-order accurate in the cell size (and exact for a single
diffusing species on a lattice).

Reaction rates keep their concentration-based prefactors (mass action in $c$, and
$c/c_{\mathrm{ref}}$ in the Butler–Volmer forms), while the driving force is always the affinity
from $\bar\mu$. So equilibrium is exact for any statistics, and only the kinetics away from it
depend on that choice of prefactor.

Cross-species transport coefficients (an Onsager matrix) aren't included: each species moves
down its own $\bar\mu$. They should arrive together with cross chemical capacitances, as two
halves of the same network ([roadmap](../ROADMAP.md)). Note that $K$ itself can already be
non-diagonal, e.g. on a shared lattice.

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
