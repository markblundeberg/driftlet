# driftlet

A small, dependency-free JavaScript library for **1D drift–diffusion–reaction** problems:
Poisson–Nernst–Planck transport of ions, electrons and holes, with bulk reactions and electrode
kinetics. It's built to be fast enough to run behind a slider in a web page.

**Status: pre-alpha, nothing to use yet.**

What it's for: teaching, live interactive demos and quick exploration of electrochemical cells,
membranes, junctions, double layers and simple semiconductor devices.

What it's not: a TCAD or battery-modelling package. For 2D/3D, thermal, concentrated-solution
transport or parameter fitting, see [ChargeTransport.jl](https://github.com/WIAS-PDELib/ChargeTransport.jl),
[Driftfusion](https://github.com/barnesgroupICL/Driftfusion),
[SIMsalabim](https://github.com/kostergroup/SIMsalabim) or [PyBaMM](https://github.com/pybamm-team/PyBaMM).

## Development

No dependencies. Tests use node's built-in runner (Node ≥ 20):

```sh
npm test
```

## Licence

[0BSD](LICENSE): do anything you like with it.
