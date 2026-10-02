# Data library

`driftlet/kit` carries a small library of bulk properties, each with its source, so a demo can
start from vetted numbers. It's kept small on purpose, because a wrong library number is worse
than a wrong demo number. The ions and semiconductors are checked in the tests against
independent tables: $\mu^\circ$ against the electrochemical series, D against limiting conductivities,
band data against the quoted intrinsic concentrations.

The library holds **bulk properties only**. How levels line up at an interface, a surface's
work function, a pzc and rate constants all belong to a particular interface or surface, so
you give them yourself (see [conventions](conventions.md)).

## Aqueous ions, 25 °C

`aqueous(names, { epsr, material })` returns water holding the named ions, as a library piece
for `build()`: each ion a species with $c_{\mathrm{ref}}$ = 1000 mol/m³ (1 M), and the material
`water` with its D and $\mu^\circ$. `epsr` defaults to water's 78.30 at 25 °C (Malmberg and
Maryott 1956); `epsr: 0` makes the solution strictly neutral, for macroscopic models without
double layers. The data are in `IONS`, and $\mu^\circ$ of liquid water is `H2O`
(−237.129 kJ/mol), for reactions such as H₂O ⇌ H⁺ + OH⁻.

| Ion | z | D (10⁻⁹ m²/s) | $\mu^\circ$ (kJ/mol) |
|---|---|---|---|
| H+ | +1 | 9.311 | 0 |
| Li+ | +1 | 1.029 | −293.31 |
| Na+ | +1 | 1.334 | −261.905 |
| K+ | +1 | 1.957 | −283.27 |
| Ag+ | +1 | 1.648 | 77.107 |
| Cu2+ | +2 | 0.714 | 65.49 |
| Zn2+ | +2 | 0.703 | −147.06 |
| Fe2+ | +2 | 0.719 | −78.9 |
| Fe3+ | +3 | 0.604 | −4.7 |
| OH- | −1 | 5.273 | −157.244 |
| Cl- | −1 | 2.032 | −131.228 |
| NO3- | −1 | 1.902 | −111.25 |
| SO42- | −2 | 1.065 | −744.53 |

- **D**: diffusivities at infinite dilution, from the CRC Handbook of Chemistry and Physics,
  "Ionic conductivity and diffusion at infinite dilution". They're limiting values, so at
  real concentrations ions move somewhat slower.
- **$\mu^\circ$**: standard Gibbs energies of formation, from Wagman et al., "The NBS tables of
  chemical thermodynamic properties" (J. Phys. Chem. Ref. Data 11, Suppl. 2, 1982), on the usual
  table convention: $\Delta_{\mathrm{f}} G^\circ(\mathrm{H}^+) = 0$ and elements in their standard
  states at 0. With these $\mu^\circ$, the electronic level of a standard hydrogen electrode in
  the solution is $\phi$ itself (see [conventions](conventions.md)), so $V - \phi$ reads against
  SHE.
- **Standard state**: the tables use 1 mol/kg, and the library uses 1 M. In dilute water at 25 °C
  they differ by 0.3%, about 7 J/mol in $\mu^\circ$.
- **Iron**: Fe²⁺'s $\Delta_{\mathrm{f}} G^\circ$ gives $E^\circ(\mathrm{Fe}^{2+}/\mathrm{Fe})$ =
  −0.41 V, where the CRC series lists −0.447 V; sources disagree on iron. The difference between
  Fe³⁺ and Fe²⁺ gives $E^\circ(\mathrm{Fe}^{3+}/\mathrm{Fe}^{2+})$ = 0.769 V, against the
  tabulated 0.771 V, which is what a Fe³⁺/Fe²⁺ redox demo needs.

## Semiconductors, 300 K

`semiconductor(name, { material })` returns electrons `e-` and holes `h+` and the material:
$c_{\mathrm{ref}} = N_c$ and $N_v$, D from the mobilities by the Einstein relation at 300 K,
$\mu^\circ_{\mathrm{e}^-} = 0$ (so $\phi$ sits at the conduction band) and
$\mu^\circ_{\mathrm{h}^+} = E_g \cdot F$. The data are for 300 K, so give the device `T: 300`.
Statistics are ideal (Boltzmann); add Fermi–Dirac statistics to the material for degenerate
doping. The data are in `SEMICONDUCTORS`, in the source's units (cm⁻³, cm²/(V·s), eV).

| | $\varepsilon_r$ | $E_g$ (eV) | $\chi$ (V) | $N_c$ (cm⁻³) | $N_v$ (cm⁻³) | $\mu_n$ (cm²/(V·s)) | $\mu_p$ (cm²/(V·s)) |
|---|---|---|---|---|---|---|---|
| Si | 11.9 | 1.12 | 4.05 | 2.8 × 10¹⁹ | 1.04 × 10¹⁹ | 1500 | 450 |
| Ge | 16.0 | 0.66 | 4.0 | 1.04 × 10¹⁹ | 6.0 × 10¹⁸ | 3900 | 1900 |
| GaAs | 13.1 | 1.424 | 4.07 | 4.7 × 10¹⁷ | 7.0 × 10¹⁸ | 8500 | 400 |

- Source: Sze, *Physics of Semiconductor Devices*, 2nd ed. (1981), appendix "Properties of Ge,
  Si and GaAs at 300 K". The mobilities are for lightly doped material.
- The electron affinity $\chi$ is for vacuum-level estimates with `vacuumDipole` (anchor `'e-'`,
  see the [alignment guide](alignment.md)). It isn't used otherwise, since an alignment is never
  implied.
- **Silicon's intrinsic concentration**: these $N_c$, $N_v$ and $E_g$ give $n_i$ = 6.7 × 10⁹ cm⁻³.
  Sze quotes 1.45 × 10¹⁰ (an older measurement), and the accepted value is 9.65 × 10⁹. The
  library keeps the band data self-consistent rather than matching a quoted $n_i$. For Ge and
  GaAs, the data reproduce the quoted $n_i$ (2.4 × 10¹³ and 1.79 × 10⁶) to within 10–15%.

## Metals, 20 °C

`metal(name, { material })` returns a conductor material carrying `e-`, with its conductivity
$1/\rho$ from the CRC Handbook's resistivities of pure metals at 293 K ("Electrical resistivity of
pure metals"). The data are in `METALS`, as `rho` (10⁻⁸ Ω·m) and `conductivity` (S/m). They
aren't cross-checked against a second table in the tests; a metal's conductivity only sets the
ohmic drop inside it, which is usually negligible.

| | Ag | Cu | Au | Al | Zn | Li | Pt |
|---|---|---|---|---|---|---|---|
| $\rho$ (10⁻⁸ Ω·m) | 1.587 | 1.678 | 2.214 | 2.65 | 5.90 | 9.28 | 10.5 |
| $\sigma$ (10⁷ S/m) | 6.30 | 5.96 | 4.52 | 3.77 | 1.69 | 1.08 | 0.952 |

There are no work functions: a work function belongs to a surface, not a metal, and it varies
by facet and adsorbate. Give a face's `zeroCharge` (or use `vacuumZeroCharge` with a work
function you've chosen for that surface).

## Example

```js
import { Device } from 'driftlet';
import { build, layer, ohmic, half, level, SHE, aqueous, metal } from 'driftlet/kit';

// A platinum electrode in Fe³⁺/Fe²⁺ chloride at equilibrium: its Fermi level sits at the
// couple's redox level, E° against SHE when the two are equal: 0.769 V from the library's NBS values (tables list 0.771 V).
const iron = half('Fe3+ + e- = Fe2+');
const def = build({
  library: [aqueous(['H+', 'Fe3+', 'Fe2+', 'Cl-'], { epsr: 0 }), metal('Pt')],
  stack: [
    ohmic(0, ['e-']),
    layer('Pt', 1e-6),
    { reactions: [{ ...iron, k0: 1e-3 }] },
    layer('water', 1e-6, { c0: { 'H+': 100, 'Fe3+': 10, 'Fe2+': 10, 'Cl-': 150 } }),
    { reactions: [{ ...iron, k0: 1e-3 }] },
    layer('Pt', 1e-6),
    ohmic(0, ['e-']),
  ],
});
const sol = new Device(def).solve();
const g = sol.x.length >> 1;
console.log(`V_e(Pt) − V°(SHE) = ${(sol.V['e-'][0] - level(sol, SHE, { standard: true })[g]).toFixed(4)} V`);
```
