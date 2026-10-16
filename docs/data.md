# Data library

`driftlet/kit` carries a small library of bulk properties, each with its source, so a demo can
start from vetted numbers. It's kept small on purpose, because a wrong library number is worse
than a wrong demo number. The ions and semiconductors are checked in the tests against
independent tables: $`\mu^\circ`$ against the electrochemical series, D against limiting conductivities,
band data against the quoted intrinsic concentrations.

The library holds **bulk properties only**. How levels line up at an interface, a surface's
work function, a pzc and rate constants all belong to a particular interface or surface, so
you give them yourself (see [conventions](conventions.md)).

## Aqueous ions, 25 °C

`aqueous(names, { epsr, material })` returns water holding the named ions, as a library piece
for `build()`: each ion a species with $`c_{\mathrm{ref}}`$ = 1000 mol/m³ (1 M), and the material
`water` with its D and $`\mu^\circ`$. `epsr` defaults to water's 78.30 at 25 °C (Malmberg and
Maryott 1956); `epsr: 0` makes the solution strictly neutral, for macroscopic models without
double layers. The data are in `IONS`, and $`\mu^\circ`$ of liquid water is `H2O`
(−237.129 kJ/mol), for reactions such as H₂O ⇌ H⁺ + OH⁻. Together they give $`K_w`$ =
1.011e-14 at 25 °C (pK_w 13.995): neutral water holds 1.0055e-4 mol/m³ of each, which a bath
with that reaction should hold too.

| Ion | z | D (10⁻⁹ m²/s) | $`\mu^\circ`$ (kJ/mol) |
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
- **$`\mu^\circ`$**: standard Gibbs energies of formation, from Wagman et al., "The NBS tables of
  chemical thermodynamic properties" (J. Phys. Chem. Ref. Data 11, Suppl. 2, 1982), on the usual
  table convention: $`\Delta_{\mathrm{f}} G^\circ(\mathrm{H}^+) = 0`$ and elements in their standard
  states at 0. With these $`\mu^\circ`$, the electronic level of a standard hydrogen electrode in
  the solution is $`\phi`$ itself (see [conventions](conventions.md)), so $`V - \phi`$ reads against
  SHE.
- **Standard state**: the tables use 1 mol/kg, and the library uses 1 M. In dilute water at 25 °C
  they differ by 0.3%, about 7 J/mol in $`\mu^\circ`$.
- **Iron**: Fe²⁺'s $`\Delta_{\mathrm{f}} G^\circ`$ gives $`E^\circ(\mathrm{Fe}^{2+}/\mathrm{Fe})`$ =
  −0.41 V, where the CRC series lists −0.447 V; sources disagree on iron. The difference between
  Fe³⁺ and Fe²⁺ gives $`E^\circ(\mathrm{Fe}^{3+}/\mathrm{Fe}^{2+})`$ = 0.769 V, against the
  tabulated 0.771 V, which is what a Fe³⁺/Fe²⁺ redox demo needs.

## Semiconductors, 300 K and other temperatures

`semiconductor(name, { material, T })` returns electrons `e-` and holes `h+` and the material:
$`c_{\mathrm{ref}} = N_c`$ and $`N_v`$, D from the mobilities by the Einstein relation at T,
$`\mu^\circ_{\mathrm{e}^-} = 0`$ (so $`\phi`$ sits at the conduction band) and
$`\mu^\circ_{\mathrm{h}^+} = E_g \cdot F`$. The table is for 300 K, the default; `T` (150 to 500 K)
moves it as [below](#at-other-temperatures). Statistics are ideal (Boltzmann); add Fermi–Dirac
statistics to the material for degenerate doping. The data are in `SEMICONDUCTORS`, in the
source's units (cm⁻³, cm²/(V·s), eV). In the material they're in SI: $`N_c`$ and $`N_v`$ are the
`cRef` of e⁻ and h⁺ in mol/m³, at `materials.Si.species['e-'].cRef` and `['h+'].cRef` (the
`species` list carries only names and charges).

**Give the device the same `T`.** The material's gap, $`N_c`$, $`N_v`$ and D are for the
temperature it was made at, and it carries that temperature as `materials.Si.T`. A device at
another temperature warns, in `sol.warnings` (so `check()`) and `describe()`, naming both
temperatures. The device's own default is 298.15 K (25 °C, the ions' temperature), not the
table's 300 K, so a semiconductor device says `T: 300`: at 298.15 K the 300 K data's
intrinsic density is 13% low. A material written by hand has no `T` and is never warned about.

To change the temperature of a live device, give it the material for the new T in the same
`set()`. The state carries over as a warm start:

```js nocheck
dev.set({ T: 350, materials: semiconductor('Si', { T: 350 }).materials });
```

A `set({ T })` alone changes kT but not the material's data, and warns. (A material with a
`trap()` in it: remake it with `trap(semiconductor('Si', { T }), …)` the same way. The trap
takes its material's T.)

| | $`\varepsilon_r`$ | $`E_g`$ (eV) | $`\chi`$ (eV) | $`N_c`$ (cm⁻³) | $`N_v`$ (cm⁻³) | $`\mu_n`$ (cm²/(V·s)) | $`\mu_p`$ (cm²/(V·s)) |
|---|---|---|---|---|---|---|---|
| Si | 11.9 | 1.12 | 4.05 | 2.8 × 10¹⁹ | 1.04 × 10¹⁹ | 1500 | 450 |
| Ge | 16.0 | 0.66 | 4.0 | 1.04 × 10¹⁹ | 6.0 × 10¹⁸ | 3900 | 1900 |
| GaAs | 13.1 | 1.424 | 4.07 | 4.7 × 10¹⁷ | 7.0 × 10¹⁸ | 8500 | 400 |

- Source: Sze, *Physics of Semiconductor Devices*, 2nd ed. (1981), appendix "Properties of Ge,
  Si and GaAs at 300 K". The mobilities are for lightly doped material.
- The electron affinity $`\chi`$ (an energy, in eV) is for vacuum-level estimates with
  `vacuumDipole` (anchor `'e-'`, `offset: chi`, the same number in volts; see the
  [alignment guide](alignment.md)). It isn't used otherwise, since an alignment is never implied.
- **Silicon's intrinsic concentration**: these $`N_c`$, $`N_v`$ and $`E_g`$ give $`n_i`$ = 6.7 × 10⁹ cm⁻³.
  Sze quotes 1.45 × 10¹⁰ (an older measurement), and the accepted value is 9.65 × 10⁹. The
  library keeps the band data self-consistent rather than matching a quoted $`n_i`$. For Ge and
  GaAs, the data reproduce the quoted $`n_i`$ (2.4 × 10¹³ and 1.79 × 10⁶) to within 10–15%.

### At other temperatures

With `T`, each 300 K number moves by a standard temperature law:

- **Band gap**, by Varshni's form $`E_g(T) = E_{g0} - \alpha T^2/(T + \beta)`$ (Varshni, Physica 34,
  149, 1967), with Thurmond's parameters as tabulated by Sze (Thurmond, J. Electrochem. Soc. 122,
  1133, 1975). Only the shift from 300 K is used, on the table's gap,
  $`E_g(T) = E_g^{\text{table}} + [E_g^{\text{V}}(T) - E_g^{\text{V}}(300\ \mathrm{K})]`$, so that 300 K is the
  table exactly. Varshni's own value at 300 K differs from the table by a few meV.
- **$`N_c`$ and $`N_v`$** $`\propto T^{3/2}`$, from parabolic bands with fixed effective masses.
- **Mobilities**, by the power laws for scattering by lattice vibrations,
  $`\mu \propto T^{-\gamma}`$ (Sze, ch. 1). They are approximate and hold for lightly doped material.
  Where scattering by ionized impurities matters (heavy doping, or the cold end of the range), the real
  mobility is lower and rises less on cooling. D is then $`\mu kT/q`$ at T, so it falls more slowly
  than μ.
- **Electron affinity** $`\chi`$ and **$`\varepsilon_r`$** stay at their 300 K values.

| | $`E_{g0}`$ (eV) | $`\alpha`$ (10⁻⁴ eV/K) | $`\beta`$ (K) | Varshni at 300 K (eV) | $`\gamma_n`$ | $`\gamma_p`$ |
|---|---|---|---|---|---|---|
| Si | 1.170 | 4.73 | 636 | 1.1245 | 2.42 | 2.20 |
| Ge | 0.7437 | 4.774 | 235 | 0.6634 | 1.66 | 2.33 |
| GaAs | 1.519 | 5.405 | 204 | 1.4225 | 1.0 | 2.1 |

Outside 150–500 K the call is refused: Varshni's fits and the lattice-limited mobilities are
not meant for the cold, where impurities and freeze-out take over. With these laws, silicon's
$`n_i`$ grows 8.3% per kelvin near room temperature. In the tests, a silicon cell's
open-circuit voltage falls 2.14 mV/K at 300 K. That lies within Green's
$`\mathrm{d}V_{oc}/\mathrm{d}T = -(E_{g0}/q - V_{oc} + \gamma kT/q)/T`$ for the 1.5 ≤ γ ≤ 3 that
its transport allows, and the γ = 1.67 it implies matches its two parts of $`J_0`$.

**A diode's saturation current** follows from these laws. With $`n_i^2 \propto T^3 e^{-E_g/kT}`$
and $`D = \mu kT/q \propto T^{1-\gamma}`$ ($`\gamma`$ of the minority carrier that dominates),
a short base, $`J_0 \propto n_i^2 D/W`$, gives

```math
J_0 \propto T^{3 + 1 - \gamma}\, e^{-E_g(T)/kT}, \qquad E_a = -k\,\frac{\mathrm{d}\ln J_0}{\mathrm{d}(1/T)} = E_g - T\frac{\mathrm{d}E_g}{\mathrm{d}T} + (4 - \gamma)\,kT .
```

So an Arrhenius plot of $`J_0`$ doesn't read the gap. $`E_g - T\,\mathrm{d}E_g/\mathrm{d}T`$ is
the gap extrapolated linearly to 0 K, 1.196 eV for silicon at 300 K, and the prefactor adds
$`(4-\gamma)kT`$, 0.041 eV for its electrons, about 1.24 eV in all. A long base,
$`J_0 \propto n_i^2\sqrt{D/\tau}`$, has $`(1-\gamma)/2`$ in place of $`1-\gamma`$, with τ as you
give it (the library doesn't move it with T).

**Hot, lightly doped junctions disappear.** $`n_i`$ grows fast: silicon's is 3 × 10¹³ cm⁻³ at
450 K and 1.6 × 10¹⁴ at 500 K, germanium's 6 × 10¹⁵ at 450 K and 2 × 10¹⁶ at 500 K. Where $`n_i`$
approaches the doping, both sides go intrinsic and the junction (its built-in voltage, its
rectification) fades away. Compare $`n_i(T)`$ with the doping before reading a trend in V(T)
across that range.

### Alloys and other materials from these

A material the library lacks can start from a near relative's entry, with its own numbers put
in. For example, Al₀.₃Ga₀.₇As from GaAs: give the holes' `mu0` the alloy's gap ($`E_g F`$), and set
whatever else differs, such as the mobilities and $`N_c`$, $`N_v`$. A field you don't set keeps
GaAs's value, its electron mobility of 8500 cm²/(V·s) included. Then give each face between the
two its band offset as a `step` on e⁻, $`\Delta E_c`$ right minus left (the valence-band offset
follows from the two gaps). `semiconductor(name, { T })` moves only the named material's own data,
by its laws: numbers you put in stay as given, so at another temperature move them yourself. The
material still carries the `T` it was made at, for the warning.

```js
import { Device, units, GAS_CONSTANT, FARADAY } from 'driftlet';
import { build, layer, ohmic, semiconductor } from 'driftlet/kit';

// The usual estimates for x = 0.3 (E_g = 1.424 + 1.247x eV, ΔE_c ≈ 0.62 ΔE_g,
// m_e = 0.067 + 0.083x); check them against your source.
const T = 300, VT = (GAS_CONSTANT * T) / FARADAY, x = 0.3;
const gap = 1.424 + 1.247 * x, dEc = 0.62 * (gap - 1.424);
const alloy = semiconductor('GaAs', { material: 'AlGaAs', T }), m = alloy.materials.AlGaAs;
m.species['h+'] = { ...m.species['h+'], mu0: units.eV(gap) };
m.species['e-'] = {
  ...m.species['e-'],
  D: 1000e-4 * VT, // μ_n = 1000 cm²/(V·s), a rough figure
  cRef: m.species['e-'].cRef * ((0.067 + 0.083 * x) / 0.067) ** 1.5, // N_c ∝ m_e^1.5
};
// A heterojunction diode, n-AlGaAs | p-GaAs, forward biased.
const dev = new Device(build({
  T,
  library: [semiconductor('GaAs', { T }), alloy],
  stack: [
    ohmic(0), layer('AlGaAs', 1e-6, { donors: units.perCm3(1e17) }),
    { step: { species: 'e-', value: units.eV(-dEc) } }, // GaAs's conduction band is lower
    layer('GaAs', 1e-6, { acceptors: units.perCm3(1e17) }), ohmic(1.0),
  ],
  bulkReactions: [{ equation: 'e- + h+ = 0', kf: { GaAs: 1e8, AlGaAs: 1e8 } }],
}));
const sol = dev.solve();
// The wider gap keeps holes out of the n side: nearly all the current is electrons injected
// into the GaAs.
const { 'e-': e, 'h+': h } = sol.contacts.left.flux;
console.log(sol.warnings.length, (e / -h).toExponential(1)); // 0, 2.5e+7
```

## Metals, 20 °C

`metal(name, { material })` returns a conductor material carrying `e-`, with its conductivity
$`1/\rho`$ from the CRC Handbook's resistivities of pure metals at 293 K ("Electrical resistivity of
pure metals"). The data are in `METALS`, as `rho` (10⁻⁸ Ω·m) and `conductivity` (S/m). They
aren't cross-checked against a second table in the tests; a metal's conductivity only sets the
ohmic drop inside it, which is usually negligible.

| | Ag | Cu | Au | Al | Zn | Li | Pt |
|---|---|---|---|---|---|---|---|
| $`\rho`$ (10⁻⁸ Ω·m) | 1.587 | 1.678 | 2.214 | 2.65 | 5.90 | 9.28 | 10.5 |
| $`\sigma`$ (10⁷ S/m) | 6.30 | 5.96 | 4.52 | 3.77 | 1.69 | 1.08 | 0.952 |

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
