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

## Demos

Live pages in [`demos/`](demos/), each solving its device in the browser as you move a slider.
Serve the repository root (e.g. `python3 -m http.server`) and open `/demos/`.

| | |
|---|---|
| [![pn junction](demos/screenshots/pn.png)](demos/pn.html) | [![double layer](demos/screenshots/double-layer.png)](demos/double-layer.html) |
| **pn junction**: species voltages and band edges under bias | **Double layer**: dilute vs crowded (lattice-gas) ions |
| [![cell](demos/screenshots/cell.png)](demos/cell.html) | [![insertion](demos/screenshots/insertion.png)](demos/insertion.html) |
| **Ag \| AgNO₃ \| Ag**: polarization transient and impedance | **Intercalation host**: constant-current lithiation vs the OCV |

## How to think about it

driftlet insists on thermodynamically honest concepts ([conventions](docs/conventions.md)):

- **Everything you set or measure is an electrochemical potential μ̄** (or a difference of
  them). A bias sets Δμ̄ of electrons between terminals; a gate sets its metal's μ̄_e⁻. Nothing
  you set is an electrostatic potential.
- **φ is bookkeeping.** Each material's standard chemical potentials anchor its own φ, and φ
  jumps at every interface between different materials.
- **Band alignment is a property of each interface** and must be given explicitly: nothing is
  assumed by default. If vacuum-level estimates are all you have, give each side an anchor and
  an offset (electron affinity, work function, absolute SHE potential, …) and driftlet lines
  up the vacuum levels. The [alignment guide](docs/alignment.md) covers the recipe and its
  caveats.
- The species voltage `V_i = μ̄_i / (z_i F)` and standard level `V°_i` are available as views,
  as in [ESBD](https://marklundeberg.com/esbd/) diagrams.

A device is a line of **regions** (each a **material** plus a length, fixed charge and
initial composition), joined at **interfaces**, with a **contact** at each end and a
**circuit** across them. Supported physics:

- any mix of charged and neutral species; per-material diffusivities and standard potentials;
  species absent from some materials;
- ideal statistics by default, or per material: Fermi–Dirac (degenerate carriers), lattice
  gas (crowding, site filling), Redlich–Kister, Debye–Hückel, intercalation hosts described by
  their OCV curve, or your own function ([statistics](docs/statistics.md));
- Poisson electrostatics, or strictly neutral (ε = 0) materials such as macroscopic
  electrolytes;
- metal regions with only a Fermi level and a conductivity: Schottky and MOS gates as regions,
  electrode reactions at internal metal faces, bipolar electrodes;
- contacts as outside phases with known levels, joined by the same laws as internal faces:
  ohmic contacts, reversible electrodes, baths, pinned barriers, conductance links,
  Butler–Volmer electrode reactions, gates and Stern layers;
- interfaces with explicit alignment, blocking, interface resistance, or Butler–Volmer
  ion/electron transfer, and a choice of electrostatic law (pinned dipole, neutral, Helmholtz);
- bulk reactions with thermodynamically consistent mass action (recombination, water
  autoionisation, …);
- imposed flow (advection) and current-free eddy mixing;
- internal ports: reservoirs feeding a window of nodes, e.g. grounding a MOS channel;
- voltage, galvanostatic (including open circuit) and load-resistor circuits;
- steady states; transients by backward Euler, BDF2, or adaptive BDF2 with error control and
  frame budgets for animation; exact conservation bookkeeping;
- small-signal impedance spectra Z(f), with complex profiles.

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
| Time integration | BE first order, BDF2 second order; adaptive error control | ratios 2, 4 |
| Advection and eddy mixing | exact convection–diffusion profile; D + D_mix; junction EMF unchanged by mixing | 1e-12; 1e-12; 1e-5 |
| Metal regions | ohmic; Schottky face and MOS gate equal their contact forms; Ag \| AgNO₃ \| Ag with metal electrodes; bipolar electrode at V/2 | 1e-12; 1e-9; 1e-9; 1e-9 |
| Internal ports | transmission line σV tanh(L/λ)/λ (O(h²)); held level; MOS low-frequency C–V with a grounded channel | 1e-4; 1e-9; 2e-3 |
| Impedance | Macdonald blocking-electrode spectrum, 100 Hz–1 GHz; finite-length Warburg (Ag \| AgNO₃ \| Ag); DC limit = differential resistance | 3e-4; 3e-5; 1e-6 |
| Conservation | per step, and against time-integrated contact fluxes | 1e-11 relative |
| Strictly neutral limit (ε = 0) | Planck EMF; polarization with no overlimiting; Donnan at neutral faces | 1e-5; 1e-4; 1e-9 |
| Fermi–Dirac statistics | 𝓕_{±1/2} vs quadrature; degenerate bulk; accumulation charge via 𝓕_{3/2} | 1e-13; 1e-12; 1e-3 |
| Crowded double layer (lattice gas) | Kilic–Bazant–Ajdari charge, up to ψ = 40 V_T; custom function reproduces it | 5e-4; 1e-10 |
| Non-ideal transport | steady flux −(D/L)ΔP of the grand potential (lattice exact, Redlich–Kister O(h²)) | 1e-12; 1e-4 |
| Debye–Hückel | junction EMF 2t₊(RT/F) ln(a₁/a₂) with activities | 1e-6 |
| Intercalation host (OCV) | composition vs table and isotherm; chemical diffusion flux and relaxation rate | 1e-12; 1e-5; 1e-4 |

## Performance

`npm run bench` runs typical interactive workloads and compares them with a committed
baseline. Matrix factorisations (≈ Newton iterations) and block work are deterministic and
checked in CI; times are from a desktop Ryzen 7600X in Node 22:

| Task | Time | Factorisations |
|---|---|---|
| Linear solve, 300 nodes × 7 unknowns (factor + solve) | 0.4 ms | 1 |
| pn diode (264 nodes): cold equilibrium | 5 ms | 9 |
| pn diode: I–V sweep 0 → 0.6 V, 31 points | 34 ms | 137 |
| pn diode: warm jump from +0.4 V to −1 V | 9 ms | 51 |
| pn diode: adaptive transient, 0 → 0.5 V, 100 ns | 200 ms | 1026 |
| pn diode: impedance at 20 frequencies | 13 ms | 21 |
| MOS with a metal gate: C–V sweep, 26 points | 11 ms | 126 |
| Ag \| AgNO₃ \| Ag with double layers: sweep 0 → 0.1 V, 21 points | 27 ms | 170 |
| Ag \| AgNO₃ \| Ag, neutral: adaptive transient over 1 s | 36 ms | 287 |
| Bipolar Ag electrode: sweep 0 → 1 V, 11 points | 300 ms | 555 |

The library has no dependencies and does no DOM access, so it runs in a Web Worker. (A device
definition is plain data and can be posted to a worker. Devices using custom-function
statistics can't be, and must be built inside the worker.)

## Development

No dependencies. Tests use node's built-in runner (Node ≥ 20):

```sh
npm test         # the validation suite
npm run bench    # typical workloads against bench/baseline.json (--save to update it)
npm run types    # TypeScript declarations from the JSDoc, into types/ (fetches TypeScript via npx)
```

Design notes and plans: [numerics](docs/numerics.md), [conventions](docs/conventions.md),
[statistics](docs/statistics.md), [roadmap](ROADMAP.md).

## Licence

[0BSD](LICENSE): do anything you like with it.
