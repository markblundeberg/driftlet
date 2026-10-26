# More templates

Templates for devices outside [llms.txt](../llms.txt)'s core set, written the same way: each
block runs in the test suite and prints its closed-form check. Read llms.txt first; these build
on it.

- semiconductors: [a Schottky diode](#a-schottky-diode), [a thin-film transistor](#a-thin-film-transistor-a-capacitance-along-a-channel),
  [traps in a layer of holes alone](#traps-in-a-layer-of-holes-alone-a-space-charge-limited-diode);
- electrochemistry: [corrosion under a drop](#corrosion-under-a-drop-an-electrode-spread-through-a-film)
  (an electrode spread through a film, and a cross-section that varies),
  [an electrode's double-layer capacitance](#an-electrodes-double-layer-capacitance-gouychapmanstern),
  [a semiconductor electrode's Mott–Schottky plot](#a-semiconductor-electrodes-mottschottky-plot),
  [a concentration cell with transference](#a-concentration-cell-with-transference);
- physiology: [a voltage clamp](#a-voltage-clamp), [an action potential](#an-action-potential),
  [one travelling along an axon](#an-action-potential-travelling-along-an-axon).

## a Schottky diode

A barrier pinned at the metal (`phi: 'pinned'`, `zeroCharge` the barrier height for Si, whose
$`\mu^\circ_{\mathrm{e}^-}`$ is 0), crossed at the emission velocity: a `velocity` link, $`N_{\mathrm{in}} = v\,(c_{\mathrm{eq}} - c)`$. With
`'equilibrium'` links (infinite velocity) it's diffusion theory instead, 20 times the current here.

```js
import { Device, units, FARADAY, GAS_CONSTANT } from 'driftlet';
import { build, layer, ohmic, semiconductor } from 'driftlet/kit';

const T = 300, VT = (GAS_CONSTANT * T) / FARADAY, Si = semiconductor('Si');
const phiB = 0.85, A = 1.12e6; // the barrier (V); Richardson's constant for n-Si, A/(m²K²)
const v = (A * T * T) / (FARADAY * Si.materials.Si.species['e-'].cRef); // Bethe's emission velocity, 2.2e4 m/s (cRef is N_c)
const schottky = (V) => ({
  V,
  terminal: 'e-',
  species: { 'e-': { type: 'velocity', v }, 'h+': { type: 'velocity', v, offset: 0 } },
  phi: 'pinned',
  zeroCharge: phiB, // E_c − E_F at the metal
});
const dev = new Device(build({ T, library: [Si], stack: [schottky(0), layer('Si', units.um(2), { donors: units.perCm3(1e16) }), ohmic(0)], grid: { hmin: units.nm(1), hmax: units.nm(20) } }));
for (const V of [0.1, 0.2, 0.3]) {
  dev.set({ contacts: { left: { V } } }); // the metal positive: forward
  const J = dev.solve().current, te = A * T * T * Math.exp(-phiB / VT) * Math.expm1(V / VT);
  console.log(`${V} V: J ${J.toExponential(3)} A/m², thermionic emission ${te.toExponential(3)}`); // about 4% below: the depletion region in series
}
```

The same `velocity` link at an ohmic contact is its surface recombination velocity.

## a thin-film transistor (a capacitance along a channel)

A **port with a `capacitance`** holds charge $`\sigma = C\,(V - \mathtt{zeroCharge} - \phi)`$ per area of electrode
along its window, balanced by the region's: a gate along a channel, a porous electrode's double
layer. Here x runs from source to drain with the gate spread along it: the gradual-channel
approximation, exact for a thin film until it pinches off.

```js
import { Device, GAS_CONSTANT, FARADAY } from 'driftlet';

const T = 300, VT = (GAS_CONSTANT * T) / FARADAY;
const t = 20e-9, L = 10e-6, W = 100e-6, C = 1e-3, mu = 1e-3; // film, length, width (m); F/m² of gate; m²/(V·s)
const tft = (Vg, Vd) =>
  new Device({
    T,
    species: [{ name: 'e-', z: -1, cRef: 46 }], // an undoped film's electrons (N_c 2.8e19 cm⁻³)
    materials: { film: { epsr: 0, species: { 'e-': { D: mu * VT, mu0: 0 } } } }, // ε = 0: F·n·t balances σ
    regions: [{ name: 'channel', material: 'film', length: L }],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' }, // source
      right: { V: Vd, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' }, // drain
    },
    ports: [{ name: 'gate', region: 'channel', V: Vg, area: 1 / t, capacitance: { C, zeroCharge: 0.3 } }],
    grid: { hmin: L / 2000, hmax: L / 100, ratio: 1.05 }, // fine at the drain, where it pinches off
  });
// The output characteristic at V_G = 2 V: linear, then saturating past V_D ≈ V_G − 0.3 V.
for (const Vd of [0.1, 0.5, 1, 1.5, 2, 3]) {
  const sol = tft(2, Vd).solve();
  const I = -sol.current * t * W; // sol.current is per m² of the film's section
  console.log(`V_D ${Vd.toFixed(1)} V: I_D ${(I * 1e6).toFixed(2)} µA, σ at the drain ${(sol.ports[0].sigma.at(-1) * 1e3).toFixed(3)} mC/m²`);
}
```

To sweep the gate on one device instead, patch the port: `dev.set({ ports: { gate: { V } } })`.
An organic electrochemical transistor is the same with a volumetric capacitance C*:
`area: 1 / t, capacitance: { C: Cstar * t }`.

## traps in a layer of holes alone (a space-charge-limited diode)

An organic layer with holes and no dopant can't be neutral anywhere, so `ohmic()`'s `'bulk'` law
can't hold at its contacts: **pinned contacts** set the hole density at each edge instead
(`zeroCharge` = −barrier for holes of $`\mu^\circ = 0`$, density $`N_v e^{-\mathrm{barrier}/kT}`$). A **trap**
is two immobile species, empty and filled, and a reaction between them; its depth is the filled
one's $`\mu^\circ`$. In a steady state immobile traps sit at their reaction's equilibrium at every
node, under current too, so their filling is exactly $`p/(p + p_1)`$; past the trap-filled limit
the current climbs toward the trap-free space-charge-limited one.

```js
import { Device, units, FARADAY, GAS_CONSTANT } from 'driftlet';
import { build, layer } from 'driftlet/kit';

const T = 300, VT = (GAS_CONSTANT * T) / FARADAY, Nv = units.perCm3(1e21), Nt = units.perCm3(1e17), Et = 0.45;
// Holes alone (an organic semiconductor's), and a trap level Et above the valence band: empty T0,
// filled T+, h⁺ + T0 = T+. The filled one's μ° below the empty's sets the depth.
const lib = {
  species: [{ name: 'h+', z: 1, cRef: Nv }, { name: 'T0', z: 0, cRef: 1 }, { name: 'T+', z: 1, cRef: 1 }],
  materials: { org: { epsr: 3, T, species: { 'h+': { D: 1e-8 * VT, mu0: 0 }, T0: { D: 0, mu0: 0 }, 'T+': { D: 0, mu0: -units.eV(Et) } } } },
};
// Pinned contacts set the hole density at each edge: Nv at the anode, Nv·e^(−0.5 eV/kT) at the cathode.
const contact = (V, barrier) => ({ V, terminal: 'h+', species: { 'h+': 'equilibrium' }, phi: 'pinned', zeroCharge: -barrier });
const diode = (V, traps) => build({
  T, library: [lib],
  stack: [contact(V, 0), layer('org', 100e-9, { c0: { T0: Nt / 2, 'T+': Nt / 2 } }), contact(0, 0.5)],
  bulkReactions: traps ? [{ equation: 'h+ + T0 = T+', kf: { org: 1e6 } }] : [],
  grid: { hmin: 0.01e-9, hmax: 1e-9, ratio: 1.15 },
});
const p1 = Nv * Math.exp(-Et / VT); // the hole density at which a trap is half full
for (const V of [1, 3, 10]) {
  const s = new Device(diode(V, true)).solve(), free = new Device(diode(V, false)).solve();
  // Immobile traps in a steady state sit at their reaction's equilibrium at every node: T+/Nt = p/(p + p1).
  const off = Math.max(...Array.from(s.c['h+'], (p, g) => Math.abs(s.c['T+'][g] / Nt - p / (p + p1))));
  const mid = s.x.length >> 1;
  console.log(`${V} V: J ${s.current.toExponential(3)} A/m², ${(s.current / free.current).toFixed(3)} of the trap-free; traps ${(100 * s.c['T+'][mid] / Nt).toFixed(1)}% full mid-layer; occupancy off p/(p + p1) by ${off.toExponential(1)}`);
}
```

## corrosion under a drop (an electrode spread through a film)

A **port with reactions** is an electrode spread through a strictly neutral region: a metal under
a film of water, a crevice's walls, a porous electrode. Its `terminal` is the metal's carrier,
`area` its area per volume (1/h under a film h thick), each reaction Butler–Volmer per area, and
a `capacitance` the double layer. A **`geometry`** gives the device a cross-section A(x): here x
is the radius of a drop, $`A = 2\pi r\,h(r)`$, and currents and amounts become totals (A, mol).

```js
import { Device, GAS_CONSTANT, FARADAY } from 'driftlet';
import { aqueous, polarization, H2O } from 'driftlet/kit';

// Salt water on iron: a drop 3 mm across. The iron dissolves and reduces O₂, which comes in from
// the air through the drop's top, fastest at the thin rim: the centre corrodes, the rim doesn't.
const T = 298.15, RT = GAS_CONSTANT * T, R = 3e-3, h0 = 1e-3, cO2 = 0.26, D = 2e-9;
const lib = aqueous(['Na+', 'Cl-', 'Fe2+', 'OH-', 'H+'], { epsr: 0 }); // strictly neutral water
lib.species.push({ name: 'O2', z: 0, cRef: 1000 }, { name: 'e-', z: -1 }); // e⁻: the iron's, in no material
lib.materials.water.species.O2 = { D, mu0: 16.4e3 };
const r = Array.from({ length: 61 }, (_, j) => (j / 60) * R);
const h = r.map((x) => Math.max(h0 * (1 - (x / R) * (x / R)), 0.02 * h0)); // the drop's height
const dev = new Device({
  T,
  species: lib.species,
  materials: lib.materials,
  geometry: { area: { x: r, values: r.map((x, j) => 2 * Math.PI * x * h[j]) } }, // ions squeeze through the thin rim
  regions: [{ name: 'drop', material: 'water', length: R, c0: { 'Na+': 500, 'Cl-': 500 + 2e-6, 'Fe2+': 1e-6, 'OH-': 1e-4, 'H+': 1e-4, O2: cO2 } }],
  contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } }, // closed: the centre, and the rim
  bulkReactions: [
    { equation: 'H+ + OH- = H2O', fixed: { H2O }, kf: { water: 1.4e8 } }, // H2O: water's μ, −237.129 kJ/mol
    // O₂ from the air through the drop: a rate k_f (1 − c/c_sat), k_f = D c_sat/h² (not a port's
    // exchange link, whose rate is linear in μ rather than c)
    { equation: 'Air = O2', fixed: { Air: 16.4e3 + RT * Math.log(cO2 / 1000) }, kf: { water: { x: r, values: h.map((y) => (D * cO2) / (y * y)) } } },
  ],
  ports: [{
    name: 'iron', region: 'drop', terminal: 'e-',
    V: 0, // the only metal and the only path for current, so it passes none: held, it fixes the level φ is read against
    area: { x: r, values: h.map((y) => 1 / y) }, // m² of iron per m³ of water
    capacitance: { C: 0.2, zeroCharge: -0.4 }, // the double layer: the potential moves continuously
    reactions: [
      { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 5e-6, alpha: 0.5 },
      { equation: 'O2 + 2 H2O + 4 e- = 4 OH-', fixed: { H2O }, k0: 3e-9, alpha: 0.125 },
    ],
  }],
  grid: { hmin: 5e-6, hmax: 60e-6, ratio: 1.2 },
});
const sol = dev.advance(600); // ten minutes
const iron = sol.ports[0], last = iron.x.length - 1;
// rates[q][j]: reaction q at node j of the window, forward, mol/(m²·s) of iron; anodic current +
const anodic = (j) => 2 * FARADAY * -iron.rates[0][j] - 4 * FARADAY * iron.rates[1][j];
console.log(`iron at ${(iron.V - sol.phi[0]).toFixed(3)} V against the water at the centre (SHE, with table μ°)`);
console.log(`net current on the iron: ${anodic(0).toFixed(2)} A/m² at the centre, ${anodic(last).toFixed(2)} A/m² at the rim`);
// The Evans diagram at the centre: each reaction's current with the iron's level moved, the water held.
const levels = [-0.55, -0.5, -0.45, -0.4].map((v) => v + sol.phi[0]);
const evans = polarization(dev, sol, { port: 'iron', x: 0 }, levels);
for (const rx of evans.reactions) console.log(`${rx.equation}: ${Array.from(rx.current, (i) => i.toFixed(2)).join(', ')} A/m²`);
```

Species can take an electrode's surface sites (`surface`, Langmuir coverages), and a reaction
with `bare: true` runs only on bare metal: a passive film (docs/device.md).

## an electrode's double-layer capacitance (Gouy–Chapman–Stern)

A blocking metal in dilute KCl, the solution strictly neutral (`epsr: 0`) and the face's double
layer below the grid: a Stern layer (`phi: { type: 'capacitive', C }`) and its diffuse layer
(`diffuse: { epsr }`). The capacitance from `impedance()` at a low frequency has the classic
minimum at the potential of zero charge, where the diffuse layer is thickest, and follows Stern and
Gouy–Chapman in series, $`1/C = 1/C_H + 1/(\varepsilon\kappa\cosh(F\psi/2RT))`$, with the drop ψ
the solution reports. (The salt the layer takes in or gives up comes from the bath by diffusion:
at 0.1 Hz across this 100 µm it shows, 2e-5 at ±0.2 V.)

```js
import { Device, EPS0, FARADAY, GAS_CONSTANT } from 'driftlet';
import { build, layer, ohmic, bath, aqueous, metal } from 'driftlet/kit';

const T = 298.15, VT = (GAS_CONSTANT * T) / FARADAY, CH = 0.2, epsr = 78.4; // a Stern layer of 20 µF/cm²
const salt = { 'K+': 10, 'Cl-': 10 }; // 10 mM KCl
const cell = (V) => build({
  T, library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')],
  stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: CH }, zeroCharge: 0, diffuse: { epsr } }, layer('water', 100e-6, { c0: salt }), bath(salt, 'Cl-')],
  grid: { hmin: 1e-6, hmax: 10e-6 },
});
// The potential of zero charge: the solution's φ at the face (zeroCharge 0: V_F − φ_plane there).
const first = new Device(cell(0)).solve(), pzc = first.phi[first.region.indexOf(1)];
const kappa = Math.sqrt((2 * FARADAY ** 2 * 10) / (epsr * EPS0 * GAS_CONSTANT * T)); // 1/Debye length
for (const E of [-0.2, -0.1, 0, 0.1, 0.2]) {
  const dev = new Device(cell(pzc + E));
  const psi = dev.solve().interfaces[0].diffuse.right.psi; // the diffuse layer's drop
  const f = 1e-3, { Z } = dev.impedance([f]), C = -1 / (2 * Math.PI * f * Z.im[0]);
  const gcs = 1 / (1 / CH + 1 / (epsr * EPS0 * kappa * Math.cosh(psi / (2 * VT))));
  console.log(`${E.toFixed(1)} V from the pzc: C ${(100 * C).toFixed(3)} µF/cm², Gouy–Chapman–Stern ${(100 * gcs).toFixed(3)}`);
}
```

## a semiconductor electrode's Mott–Schottky plot

A semiconductor dipped in an electrolyte, the face between them **blocking** (no reactions), with
its **alignment** given (a `dipole`: what puts silicon's levels against the solution's is a
choice, here Si's electron affinity against an absolute SHE, not a datum). Its capacitance from
`impedance()` at a frequency where the space charge is what's measured, $`1/C^2`$ against the
voltage across the depletion range, gives back the doping from the slope and the flat band from
the intercept less kT/q. Past about $`2\phi_F`$ of depletion, holes gather at the face (inversion):
keep the fit below it.

```js
import { Device, units, EPS0, FARADAY, GAS_CONSTANT } from 'driftlet';
import { build, layer, ohmic, bath, semiconductor, aqueous } from 'driftlet/kit';

const T = 298.15, VT = (GAS_CONSTANT * T) / FARADAY, ND = units.perCm3(1e16);
const water = { 'K+': 100, 'Cl-': 100 }; // 0.1 M KCl
const lib = [semiconductor('Si', { T }), aqueous(['K+', 'Cl-'])];
// n-Si | a blocking face | 40 nm of the solution, resolved | the bath, read against an electrode reversible to Cl⁻.
// The face's dipole aligns the two (here Si's electron affinity against the solution: an assumption, not data).
const cell = (V) => build({
  T, library: lib,
  stack: [ohmic(V), layer('Si', units.um(2), { donors: ND }), { dipole: 0.39 }, layer('water', 40e-9, { c0: water }), bath(water, 'Cl-')],
  grid: { hmin: 0.1e-9, hmax: 5e-9, ratio: 1.1 },
});
const pts = [];
for (const V of [-1.3, -1.2, -1.1, -1.0]) { // the silicon positive: its electrons pushed back, a depletion layer
  const dev = new Device(cell(V));
  dev.solve();
  const { f, Z } = dev.impedance([1e5]); // 100 kHz: the space charge, the double layer in series negligible
  const w = 2 * Math.PI * f[0], C = -Z.im[0] / (w * (Z.re[0] ** 2 + Z.im[0] ** 2));
  pts.push([V, 1 / C ** 2]);
}
// 1/C² = 2/(F ε N_D) (V − V_fb − kT/q) in depletion (N_D in mol/m³, as driftlet's are)
const slope = (pts.at(-1)[1] - pts[0][1]) / (pts.at(-1)[0] - pts[0][0]), V0 = pts[0][0] - pts[0][1] / slope;
console.log(`N_D from the slope: ${units.toPerCm3(2 / (FARADAY * 11.9 * EPS0 * slope)).toPrecision(5)} cm⁻³ (1e16 put in)`);
// The model's own flat band: where the silicon holds no charge at the face (interfaces[0].D = 0).
let lo = -2, hi = -1.3;
for (let k = 0; k < 40; k++) {
  const m = (lo + hi) / 2, D = new Device(cell(m)).solve().interfaces[0].D;
  if (D < 0) lo = m;
  else hi = m;
}
console.log(`flat band: ${(V0 - VT).toFixed(4)} V from the intercept less kT/q, ${hi.toFixed(4)} V where the face holds no charge`);
```

## a concentration cell with transference

Two electrodes reversible to Cl⁻ (`bath(c, 'Cl-')`), one held and one at open circuit, with
two strengths of HCl between them meeting at a **liquid junction**. For one salt the voltage
doesn't depend on the junction's profile, so baths at both ends holding it stationary give
the free-diffusion junction's value: $`E = -2t_+ (RT/F)\ln(a_{\pm,2}/a_{\pm,1})`$, the junction's
share $`(2t_+ - 1)/(2t_+)`$. Activity coefficients go on the water's `statistics`.

```js
import { Device, FARADAY, GAS_CONSTANT, EPS0, AVOGADRO } from 'driftlet';
import { build, layer, bath, aqueous, IONS } from 'driftlet/kit';

const T = 298.15, VT = (GAS_CONSTANT * T) / FARADAY;
const dilute = { 'H+': 1, 'Cl-': 1 }, strong = { 'H+': 10, 'Cl-': 10 }; // mol/m³: 1 and 10 mM HCl
const tPlus = IONS['H+'].D / (IONS['H+'].D + IONS['Cl-'].D);
// Ag|AgCl at each end: a bath read against an electrode reversible to Cl⁻. The left one held, the
// right at open circuit; between them a junction held stationary by the two baths (strictly
// neutral water, ε = 0: the junction is 2 mm wide, its double layers nothing to it).
for (const activity of [null, { type: 'debye-huckel', species: ['H+', 'Cl-'], epsr: 78.4, a: 0 }]) {
  const lib = aqueous(['H+', 'Cl-'], { epsr: 0 });
  if (activity) lib.materials.water.statistics = [activity];
  const s = new Device(build({
    T, library: [lib],
    stack: [bath(dilute, 'Cl-'), layer('water', 2e-3, { c0: dilute }), { phi: 'neutral' }, layer('water', 2e-3, { c0: strong }), bath(strong, 'Cl-', { I: 0 })],
    grid: { hmin: 1e-7, hmax: 2e-5 },
  })).solve();
  // E = −2 t₊ (RT/F) ln(a±₂/a±₁), the limiting law's ln γ± = −ℓ_B κ/2 (Bjerrum length, Debye κ)
  const eps = EPS0 * (activity?.epsr ?? 1), lB = FARADAY ** 2 / (4 * Math.PI * eps * GAS_CONSTANT * T * AVOGADRO);
  const lnGamma = (c) => (activity ? (-lB * Math.sqrt((2 * FARADAY ** 2 * c) / (eps * GAS_CONSTANT * T))) / 2 : 0);
  const E = s.terminals.right.V - s.terminals.left.V, theory = -2 * tPlus * VT * (Math.log(10) + lnGamma(10) - lnGamma(1));
  console.log(`${activity ? 'Debye–Hückel' : 'ideal'}: ${(E * 1e3).toFixed(3)} mV, theory ${(theory * 1e3).toFixed(3)}; the junction's share ${(2 * tPlus - 1) / (2 * tPlus)}`);
}
```

## a voltage clamp

A face's `gates` (each a fraction open, $`dx/dt = \alpha(V)(1 - x) - \beta(V)x`$ with $`V = \phi_{\mathrm{right}} - \phi_{\mathrm{left}}`$)
scale its links (`gates: { m: 3, h: 1 }`). `hodgkinHuxley({ T, linear: true })` (kit) writes the
squid axon's as Hodgkin and Huxley had them, conductances $`g(V - E)`$. Between two baths without a
reference species, $`V_{\mathrm{right}} - V_{\mathrm{left}}`$ is $`V_m`$, so a held step is their voltage clamp.

```js
import { Device, FARADAY, GAS_CONSTANT } from 'driftlet';
import { build, layer, bath, aqueous, hodgkinHuxley } from 'driftlet/kit';

const T = 279.45, VT = (GAS_CONSTANT * T) / FARADAY; // HH's 6.3 °C
const out = { 'Na+': 440, 'K+': 20, 'Cl-': 450, 'A-': 10 }, inside = { 'Na+': 50, 'K+': 400, 'Cl-': 52, 'A-': 398 }; // squid, mM
const water = aqueous(['Na+', 'K+', 'Cl-'], { epsr: 0 });
water.species.push({ name: 'A-', z: -1, cRef: 1000 });
water.materials.water.species['A-'] = { D: 0.5e-9, mu0: 0 };
const hh = hodgkinHuxley({ T, linear: true });
const step = -0.01; // from −65 mV to −10 mV at 1 ms
const dev = new Device(
  build({
    T,
    library: [water],
    stack: [
      bath(out), // the outside, on the left
      layer('water', 1e-6, { c0: out }),
      { phi: { type: 'capacitive', C: 0.01 }, gates: hh.gates, species: { ...hh.species, 'A-': 'blocked' } },
      layer('water', 1e-6, { c0: inside }),
      bath(inside, { t: [0, 1e-3, 1e-3, 1], values: [-0.065, -0.065, step, step] }),
    ],
    grid: { hmin: 20e-9, hmax: 0.2e-6 },
  }),
);
dev.solve(); // the rest at −65 mV, the gates at their values there
const run = dev.advance(11e-3, { dtMax: 0.02e-3, probes: [{ interface: 0, species: 'Na+' }, { interface: 0, species: 'K+' }] });
const [na, k] = run.trace.probes.map((N) => N.map((x) => -FARADAY * x)); // outward currents, A/m²
const peak = na.indexOf(Math.min(...na)), EK = VT * Math.log(out['K+'] / inside['K+']);
console.log(`I_Na peaks ${((run.trace.t[peak] - 1e-3) * 1e3).toFixed(2)} ms into the step, ${(na[peak] * 0.1).toFixed(2)} mA/cm²`); // 0.73 ms, −1.48
console.log(`g_K after 10 ms: ${((k.at(-1) / (step - EK)) * 0.1).toFixed(2)} mS/cm²`); // 21.12: HH's 36 n∞⁴ is 21.5, which n nearly reaches by 10 ms
```

## an action potential

Without `linear: true` (given `inside` and `outside`), `hodgkinHuxley` writes GHK permeabilities
matched to HH's conductances at −65 mV, which rectify as real channels do. A pump-less cell has no
steady rest (`solve()` finds Donnan), so start the transient from the concentrations given.

```js
import { Device } from 'driftlet';
import { build, layer, bath, aqueous, injector, pulse, hodgkinHuxley } from 'driftlet/kit';

const T = 279.45; // HH's 6.3 °C
const out = { 'Na+': 440, 'K+': 20, 'Cl-': 450, 'A-': 10 }, inside = { 'Na+': 50, 'K+': 400, 'Cl-': 52, 'A-': 398 }; // squid, mM
const water = aqueous(['Na+', 'K+', 'Cl-'], { epsr: 0 });
water.species.push({ name: 'A-', z: -1, cRef: 1000 });
water.materials.water.species['A-'] = { D: 0.5e-9, mu0: 0 };
const hh = hodgkinHuxley({ inside, outside: out, T }); // GHK channels; outside on the left
const axon = new Device(
  build({
    T,
    library: [water],
    stack: [
      bath(out, 'Cl-'),
      layer('water', 5e-6, { name: 'outside', c0: out }),
      { phi: { type: 'capacitive', C: 0.01 }, gates: hh.gates, species: { ...hh.species, 'A-': 'blocked' } },
      layer('water', 20e-6, { name: 'cytoplasm', c0: inside }),
      { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' },
    ],
    // the stimulus: 10 µA/cm² of K⁺ for 1 ms at 20 ms, from a pipette in the cytoplasm
    ports: [injector({ region: 'cytoplasm', from: 15e-6, to: 20e-6, species: 'K+', I: pulse({ start: 20e-3, width: 1e-3, amplitude: 0.1 }) })],
    grid: { hmin: 50e-9, hmax: 1e-6 },
  }),
);
let peak = -1;
for (let k = 1; k <= 300; k++) {
  const sol = axon.advance(k * 1e-4); // 0.1 ms at a time
  peak = Math.max(peak, sol.interfaces[0].V); // V_m, and sol.interfaces[0].gates = { m, h, n }
}
console.log(`peak ${(peak * 1e3).toFixed(1)} mV`); // +12.6: GHK's currents rectify; matched at rest, at positive V Na⁺'s conductance is lower and K⁺'s higher (linear: true gives +37.9, HH's own)
```

## an action potential travelling along an axon

Here x runs along the axon: the axoplasm is one region, and the membrane a port spanning it, a
capacitance and gated channels to an outside given as a `bath` (the port's V is then the
outside's $`\phi`$). `hodgkinHuxley({ T, area })` writes HH's channels for it, `area` being the membrane
per volume, 2/a. Cable theory isn't put in; it comes out of the ions' drift. `ports[k].Vm[j]` and
`ports[k].gates[name][j]` give the membrane at each node. A myelinated axon is regions, each with
its own port (an internode's just a small capacitance). Nothing warns of a grid too coarse for
the spike's front: the speed comes out 3% low on 2 mm cells, 0.2% on 0.5 mm, for a = 238 µm
(scale $`h`$ with $`\sqrt a`$). The stimulus is a current per volume, so a thinner axon needs it as $`1/a`$.

```js
import { Device, FARADAY, GAS_CONSTANT } from 'driftlet';
import { hodgkinHuxley, injector, pulse } from 'driftlet/kit';

const T = 291.65, VT = (GAS_CONSTANT * T) / FARADAY, a = 238e-6; // 18.5 °C, a squid giant axon
const inside = { 'Na+': 50, 'K+': 400, 'Cl-': 52.5, 'A-': 397.5 }; // mM
// Outside levels that put the Nernst potentials at HH's +50, −77 and −54.4 mV (leak as Cl⁻).
const outside = { 'Na+': 50 * Math.exp(0.05 / VT), 'K+': 400 * Math.exp(-0.077 / VT), 'Cl-': 52.5 * Math.exp(0.054387 / VT) };
// One D for the small ions, giving HH's 35.4 Ω·cm.
const D = 1 / 0.354 / ((FARADAY * FARADAY) / (GAS_CONSTANT * T)) / (50 + 400 + 52.5);
const sp = (d) => ({ D: d, mu0: 0, cRef: 1000 });
const sealed = { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' };
const hh = hodgkinHuxley({ T, area: 2 / a });
const axon = new Device({
  T,
  species: [{ name: 'Na+', z: 1 }, { name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }, { name: 'A-', z: -1 }],
  materials: { axoplasm: { epsr: 0, species: { 'Na+': sp(D), 'K+': sp(D), 'Cl-': sp(D), 'A-': sp(1e-12) } } },
  regions: [{ name: 'axon', material: 'axoplasm', length: 0.06, c0: inside }],
  contacts: { left: sealed, right: sealed },
  ports: [
    { name: 'membrane', region: 'axon', V: 0, area: 2 / a, capacitance: { C: 0.01, zeroCharge: 0 }, bath: { c: outside }, gates: hh.gates, species: hh.species },
    injector({ name: 'stim', region: 'axon', from: 0, to: 1e-3, species: 'K+', I: pulse({ start: 30e-3, width: 0.2e-3, amplitude: 30 }) }),
  ],
  grid: { hmin: 2e-4, hmax: 2e-4 },
});
axon.advance(30e-3); // settle to rest (no pump: start from the concentrations, not solve())
const sol = axon.advance(34e-3, { dtMax: 5e-6, probes: [{ x: 0.02, quantity: 'phi' }, { x: 0.04, quantity: 'phi' }] });
const { t, probes } = sol.trace, at = (v) => t[v.findIndex((x) => x > 0)]; // (φ inside, the outside's held at 0)
console.log(`${(0.02 / (at(probes[1]) - at(probes[0]))).toFixed(1)} m/s`); // 18.8 here (finely resolved, 18.73; HH computed 18.8)
```
