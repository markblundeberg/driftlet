# driftlet

A small, dependency-free JavaScript library for **1D drift–diffusion–reaction** problems:
Poisson–Nernst–Planck transport of ions, electrons and holes, with bulk reactions, electrode
kinetics, interfaces between materials, and external circuits. It's built to be fast enough to
run behind a slider in a web page.

**Status: pre-alpha.** Not yet published to npm, and the API may still change.

What it's for: teaching, live interactive demos, and quick exploration of electrochemical
cells, membranes, junctions, double layers, mixed conductors and simple semiconductor devices,
all in one formulation.

What it's not: a TCAD or battery-modelling package. For 2D/3D, thermal, concentrated-solution
transport or parameter fitting, see
[ChargeTransport.jl](https://github.com/WIAS-PDELib/ChargeTransport.jl),
[Driftfusion](https://github.com/barnesgroupICL/Driftfusion),
[SIMsalabim](https://github.com/kostergroup/SIMsalabim),
[PyBaMM](https://github.com/pybamm-team/PyBaMM) or COMSOL.

## A first example

A silicon pn junction with ohmic contacts, at equilibrium and then under forward bias:

```js
import { Device, units, FARADAY } from 'driftlet';

const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19); // effective densities of states
const ohmic = (V) => ({
  V, // terminal voltage
  terminal: 'e-',
  species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, // V_h+ = V_e- at a metal
  phi: 'bulk',
});

const dev = new Device({
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: {
      epsr: 11.7,
      species: {
        'e-': { D: 36e-4, mu0: 0, cRef: Nc }, // conduction band edge as the standard level
        'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv }, // gap 1.12 eV
      },
    },
  },
  regions: [
    { name: 'n', material: 'Si', length: 2e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
    { name: 'p', material: 'Si', length: 2e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
  ],
  contacts: { left: ohmic(0), right: ohmic(0) },
  grid: { hmin: 0.5e-9, hmax: 20e-9 },
});

const eq = dev.solve();
console.log('built-in potential (V):', eq.phi[0] - eq.phi[eq.phi.length - 1]); // ≈ 0.797

dev.set({ contacts: { right: { V: 0.4 } } }); // p side at +0.4 V: forward bias
const on = dev.solve(); // warm-started from the previous solution
console.log('current (A/m²):', on.current); // ≈ −6.5 (flowing toward −x, from p to n)
```

Solutions are plain objects of `Float64Array`s, ready to plot: `x`, `phi`, and per species `c`,
`mu` (electrochemical potential μ̄), `muStd` (standard level) and the voltage views `V`,
`Vstd`. They also carry currents, contact fluxes, conservation bookkeeping and warnings; see
[the device reference](docs/device.md).

## How to think about it

driftlet insists on thermodynamically honest concepts ([conventions](docs/conventions.md)):

- **Everything you set or measure is an electrochemical potential μ̄** (or a difference of
  them). A bias sets Δμ̄ of electrons between terminals; a gate sets its metal's μ̄_e⁻. Nothing
  you set is an electrostatic potential.
- **φ is bookkeeping.** Each material's standard chemical potentials anchor its own φ, and φ
  jumps at every interface between different materials.
- **Band alignment is a property of each interface** and must be given explicitly: there is no
  Anderson or Schottky–Mott rule built in. If vacuum-level reasoning is all you have, the
  [alignment guide](docs/alignment.md) gives the recipes, and their caveats.
- The species voltage `V_i = μ̄_i / (z_i F)` and standard level `V°_i` are available as views,
  as in [ESBD](https://marklundeberg.com/esbd/) diagrams.

A device is a line of **regions** (each a **material** plus a length, fixed charge and
initial composition), joined at **interfaces**, with a **contact** at each end and a
**circuit** across them. Supported physics:

- any mix of charged and neutral species; per-material diffusivities and standard potentials;
  species absent from some materials;
- Poisson electrostatics, or strictly neutral (ε = 0) materials such as metals and macroscopic
  electrolytes;
- contacts as outside phases with known levels, joined by the same laws as internal faces:
  ohmic contacts, reversible electrodes, baths, pinned barriers, conductance links,
  Butler–Volmer electrode reactions, gates and Stern layers;
- interfaces with explicit alignment, blocking, interface resistance, or Butler–Volmer
  ion/electron transfer, and a choice of electrostatic law (pinned dipole, neutral, Helmholtz);
- bulk reactions with thermodynamically consistent mass action (recombination, water
  autoionisation, …);
- voltage, galvanostatic (including open circuit) and load-resistor circuits;
- steady states and backward-Euler transients, with exact conservation bookkeeping.

How it works numerically is in [numerics](docs/numerics.md).

## Validation

Every physics feature is tested against analytic results (`npm test`, node's built-in runner).

| What | Checked against | Tolerance |
|---|---|---|
| Equilibrium of any device | μ̄ of every species flat; zero current | 1e-9 RT or better |
| Gouy–Chapman double layer | analytic charge and full nonlinear profile | 2e-3 |
| Debye screening (linear limit) | gate in series with ε/λ_D | 1e-3 |
| Donnan potential | analytic partition, between floating layers and real baths | 1e-8 |
| Heterointerfaces with unequal dipoles (A \| B \| A) | per-face analytic double-layer split | 3e-3 |
| Floating island | Gauss's law; conserved amounts through gate sweeps | 1e-12 |
| pn junction | exact built-in potential; depletion charge; short-diode J–V | 1e-12; 1%; 2e-3 |
| Long pn diode with recombination | Shockley J–V incl. depletion recombination | 2e-3 |
| Schottky barrier (metal region \| n-Si) | surface density from the alignment; depletion charge | 5e-3; 2% |
| Liquid junction, open circuit | cell EMF 2t₊(RT/F) ln(c₁/c₂); Planck diffusion potential | 1e-4 |
| Concentration polarization | i = i_lim tanh(V/4V_T), incl. galvanostatic and load modes | 2e-4 |
| Butler–Volmer electrode | Nernst equilibrium; mixed kinetic/diffusion closed form | 5e-4 |
| Interface conductance and ion transfer | series 1/G; BV rate law at the interface state | 1e-4; 1e-8 |
| Mass action | c(H⁺)c(OH⁻) = K_w from standard potentials; moiety conservation | 1e-9 |
| Transients | RC charging of a gated island; water relaxation rate | 1%; 2e-3 |
| Conservation | per step, and against time-integrated contact fluxes | 1e-11 relative |
| Strictly neutral limit (ε = 0) | Planck EMF; polarization with no overlimiting; Donnan at neutral faces | 1e-5; 1e-4; 1e-9 |

## Performance

On a desktop Ryzen 7600X in Node 22, one linear solve for 300 nodes × 7 unknowns takes about 0.4 ms, and a
warm-started steady re-solve after a bias step takes a few milliseconds. The library has no
dependencies and does no DOM access, so it runs in a Web Worker. (A device definition is plain
data and can be posted to a worker. Devices using the future custom-function statistics can't
be, and must be built inside the worker.)

## Development

No dependencies. Tests use node's built-in runner (Node ≥ 20):

```sh
npm test
```

Design notes and plans: [numerics](docs/numerics.md), [conventions](docs/conventions.md),
[roadmap](ROADMAP.md).

## Licence

[0BSD](LICENSE): do anything you like with it.
