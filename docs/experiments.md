# Against experiment

The validation suite holds driftlet to analytic limits: Cottrell, Debye, Shockley, cable theory. This
page holds it to measurements. Each entry takes a classic measurement with its published numbers,
models it the way a user would, with the library's own data where it has them, and says how close
driftlet lands and why it misses where it does. A miss here is information, not a failure: it marks
where the physics driftlet solves stops being the physics of the experiment.

Every example below runs as part of the test suite, and the numbers in its `// =>` comments are
what it prints.

| Measurement | Measured | driftlet | Where the gap comes from |
|---|---|---|---|
| [KCl's conductivity](#the-conductivity-of-kcl-solutions), 0.1 mol/kg, 25 °C (NIST, 1991) | 1.2825 S/m | 1.483 S/m, 16% high | ions' correlated motion (Debye–Hückel–Onsager), not in the model |
| [Squid axon at rest](#the-resting-potential-of-the-squid-axon), the change on replacing sea water, nine solutions (Hodgkin and Katz, 1949) | −4 to +4 mV | within 1 mV in seven, 2 mV in two | Goldman's constant-field membrane, their model: driftlet matches their calculation |
| [Squid axon conduction velocity](#the-conduction-velocity-of-the-squid-axon) (Hodgkin and Huxley, 1952) | 21.2 m/s | 18.7 m/s, 12% low | Hodgkin and Huxley's own model gives 18.8 |

## The conductivity of KCl solutions

The measurement: Y. C. Wu, W. F. Koch and K. W. Pratt, "Proposed new electrolytic conductivity
primary standards for KCl solutions", *J. Res. Natl. Inst. Stand. Technol.* 96, 191 (1991), table 5:
the electrolytic conductivity κ of 0.01 and 0.1 mol/kg KCl from 0 to 50 °C, to 0.03%. These are
the primary standards conductivity cells are calibrated against.

The model: a 1 mm cell of the solution between two baths of it, with a ±1 mV bias, and κ = |dI/dV|·L.
The ions' D come from the library (`aqueous`, the CRC's limiting values at 25 °C, moved to T with
water's viscosity). Molality becomes concentration with pure water's density, 997.05 kg/m³. That
neglects the salt's own volume, a few tenths of a percent at 0.1 mol/kg.

```js
import { Device } from 'driftlet';
import { aqueous } from 'driftlet/kit';

// κ (µS/cm) from Wu, Koch and Pratt's table 5, for 0.01 and 0.1 mol/kg
const nist = { 0: [772.921, 7116.85], 25: [1408.23, 12824.6], 50: [2123.43, 19180.9] };
const L = 1e-3, M = 0.0745513; // the cell (m), KCl's molar mass (kg/mol)
const kappa = (tC, m) => {
  const T = 273.15 + tC, c = (m * 997.05) / (1 + m * M), bath = { c: { 'K+': c, 'Cl-': c }, reference: 'Cl-' };
  const I = (V) => new Device({
    T, ...aqueous(['K+', 'Cl-'], { T }),
    regions: [{ material: 'water', length: L, c0: { 'K+': c, 'Cl-': c } }],
    contacts: { left: { V: 0, bath }, right: { V, bath } },
  }).solve().terminals.right.current;
  return (Math.abs(I(1e-3) - I(-1e-3)) / 2e-3) * L; // S/m
};
const ratio = (tC, k) => kappa(tC, [0.01, 0.1][k]) / (nist[tC][k] * 1e-4);
console.log(`25 °C: ${ratio(25, 0).toFixed(3)} at 0.01 mol/kg, ${ratio(25, 1).toFixed(3)} at 0.1`); // => 1.060, 1.156
console.log(`0 °C: ${ratio(0, 0).toFixed(3)}, ${ratio(0, 1).toFixed(3)}`); // => 0.960, 1.035
console.log(`50 °C: ${ratio(50, 0).toFixed(3)}, ${ratio(50, 1).toFixed(3)}`); // => 1.144, 1.258
```

**At 25 °C, driftlet is 6% high at 0.01 mol/kg and 16% high at 0.1.** It gives the limiting molar
conductivity at every concentration, 149.8 S·cm²/mol at 25 °C, the Nernst–Einstein sum of the two ions'
D. A real solution conducts less as it gets more concentrated (Kohlrausch's √c law). Each ion
drags an atmosphere of the opposite charge, which lags behind it and pulls it back (the relaxation
effect), and the atmosphere, driven the other way, drags the solvent around it against the ion
(electrophoresis). That is Debye–Hückel–Onsager
theory, a correlation between ions' motions. driftlet's fluxes are each ion's own, so it doesn't
have it. Its `debye-huckel` statistics wouldn't add it either: they change activities, and a
uniform cell has no activity gradient to drive.

**Across temperature, the ratio moves too: from 0.960 at 0 °C to 1.144 at 50 °C, at 0.01 mol/kg.**
This one is the library's data, not the model. `aqueous` moves each D with Stokes–Einstein, Dη/T
held, and [data.md](data.md) says that runs low when cold and is unchecked above 25 °C. Take the
relaxation's share as the same at each temperature. Then each ratio over its own 25 °C value
estimates how far the D scaling is off. That is 9% high at 50 °C at both concentrations, and 10%
low at 0 °C (0.906 and 0.895). That the two concentrations agree is some support for the
assumption, but it isn't proof: the relaxation's share does move with ε and η. For a device away
from 25 °C, give `aqueous` your own `D` at T.

## The resting potential of the squid axon

The measurement: A. L. Hodgkin and B. Katz, "The effect of sodium ions on the electrical activity of
the giant axon of the squid", *J. Physiol.* 108, 37 (1949), table 7. Squid axons in sea water, then
in nine test solutions: potassium-free or potassium-rich artificial sea water (A–C), sea water
diluted with isotonic dextrose (D–G), or with NaCl added (H, I). For each, the average change in
membrane potential, and what Goldman's constant-field equation gives for it with the permeabilities
P_K : P_Na : P_Cl = 1 : 0.04 : 0.45, and inside 345 mM K⁺, 72 mM Na⁺ and 61 mM Cl⁻. These are the
ratios physiology texts still quote, and this is where they come from: Hodgkin and Katz found them
"by trial and error" to fit this table.

The model: the membrane as a face with a permeability for each ion (driftlet's `permeability` link,
which is Goldman's flux), between 1 µm of each solution and a bath of it. The inside is held at no
current, as a cell is, and the impermeant ions that make each solution neutral are blocked. The
potential is read outside minus inside, the paper's sign. The table and its text don't give the
calculation's temperature (the experiments ran at 19–24 °C, table 3); 18.5 °C reproduces their 59 mV
at rest.

```js
import { Device } from 'driftlet';
import { build, layer, bath, aqueous } from 'driftlet/kit';

const T = 273.15 + 18.5, IONS = ['K+', 'Na+', 'Cl-'];
const inside = { 'K+': 345, 'Na+': 72, 'Cl-': 61 }; // their internal concentrations (mM)
const seaWater = { 'K+': 10, 'Na+': 455, 'Cl-': 540 }, artificial = { 'K+': 10, 'Na+': 455, 'Cl-': 587 };
// Table 7's test solutions (K, Na, Cl in mM); A–C replaced artificial sea water, D–I sea water
const tests = { A: [0, 465, 587], B: [15, 450, 587], C: [20, 445, 587], D: [7, 324, 384], E: [5, 227, 270], F: [3, 152, 180], G: [2, 91, 108], H: [10, 573, 658], I: [10, 711, 796] };
const water = aqueous(IONS, { epsr: 0, T });
for (const name of ['A-', 'M+']) { // what makes each solution neutral, impermeant
  water.species.push({ name, z: name === 'A-' ? -1 : 1, cRef: 1000 });
  water.materials.water.species[name] = { D: 5e-10, mu0: 0 };
}
const solution = (c) => {
  const s = Object.fromEntries(IONS.map((k) => [k, Math.max(c[k], 1e-6)])); // (A has no K⁺: a trace)
  const q = s['K+'] + s['Na+'] - s['Cl-'];
  return { ...s, 'A-': Math.max(q, 0) + 1, 'M+': Math.max(-q, 0) + 1 };
};
// The membrane potential, outside − inside (mV), for permeabilities P_K : P_Na : P_Cl
const membrane = (outside, [pK, pNa, pCl]) => {
  const o = solution(outside), i = solution(inside), link = (P) => (P > 0 ? { type: 'permeability', P: P * 1e-8 } : 'blocked');
  const s = new Device(build({
    T, library: [water],
    stack: [
      bath(o, 'Cl-'), layer('water', 1e-6, { c0: o }),
      { phi: { type: 'capacitive', C: 0.01 }, species: { 'K+': link(pK), 'Na+': link(pNa), 'Cl-': link(pCl), 'A-': 'blocked', 'M+': 'blocked' } },
      layer('water', 1e-6, { c0: i }), bath(i, 'Cl-', { I: 0 }),
    ],
    grid: { hmin: 50e-9, hmax: 200e-9 },
  })).solve();
  return (s.phi[0] - s.phi.at(-1)) * 1000;
};
const rest = [1, 0.04, 0.45], active = [1, 20, 0.45], refractory = [1.8, 0, 0.45];
// In sea water: at rest, at the peak of the spike, and in the refractory positive phase
console.log([rest, active, refractory].map((P) => membrane(seaWater, P).toFixed(1)).join(', ')); // => 59.4, -37.8, 74.0
const change = (P) => Object.entries(tests).map(([k, [K, Na, Cl]]) => membrane({ 'K+': K, 'Na+': Na, 'Cl-': Cl }, P) - membrane('ABC'.includes(k) ? artificial : seaWater, P));
console.log(change(rest).map((x) => x.toFixed(1)).join(', ')); // => 4.8, -2.1, -4.0, 0.8, 1.6, 2.5, 3.1, 0.1, 0.2
console.log(change(active).map((x) => x.toFixed(1)).join(', ')); // => -0.5, 0.3, 0.5, 7.6, 15.8, 25.3, 37.6, -5.1, -9.8
```

**driftlet reproduces their calculation.** In sea water they give +59, −38 and +74 mV, and
driftlet 59.4, −37.8 and 74.0: Goldman's equation by hand gives the same to 0.01 mV, and the
unstirred micrometre beside the membrane moves it by microvolts. Of their eighteen calculated
changes at rest and at the peak, driftlet's round to theirs in seventeen. The other, F at rest, is
2.53 against their +2, from 18.5 to 21 °C alike: we can't account for it.

**Against the measurement, at rest:**

| | A | B | C | D | E | F | G | H | I |
|---|---|---|---|---|---|---|---|---|---|
| observed (mV) | +3 | −2 | −4 | 0 | +2 | +2 | +4 | +1 | −2 |
| driftlet | +4.8 | −2.1 | −4.0 | +0.8 | +1.6 | +2.5 | +3.1 | +0.1 | +0.2 |

Seven of the nine are within a millivolt; the two potassium-free and salt-added cases (A, I) are 2 mV
off. That is the constant-field membrane's accuracy, Hodgkin and Katz's "reasonable agreement", and
not driftlet's: driftlet solves their equation. The absolute resting potential, 48 mV observed
against 59, they put down to a junction potential between axoplasm and the sea water in their
microelectrode, which they couldn't measure (a tentative estimate of theirs gives 14 mV).

**At the peak of the spike**, the same equation with P_Na raised to 20 undershoots the dilute
solutions: 25 and 38 mV against the observed 44 and 59 (F, G). Hodgkin and Katz read this as the
active membrane not reaching its full sodium permeability in those solutions. A peak isn't a
steady state, so this row is a test of their picture of excitation, not of a resting membrane.

## The conduction velocity of the squid axon

The measurement: A. L. Hodgkin and A. F. Huxley, "A quantitative description of membrane current and
its application to conduction and excitation in nerve", *J. Physiol.* 117, 500 (1952), p. 528. On one
fibre, radius 238 µm, axoplasm 35.4 Ω·cm, membrane 1.0 µF/cm², at 18.5 °C, the propagated action
potential travelled at 21.2 m/s. Their equations, solved for a travelling wave, gave 18.8 m/s.

The model: the axon as a region of strictly neutral axoplasm whose ions drift along it, and the
membrane as a port along the whole region, with the kit's `hodgkinHuxley` channels and the sea water
as the port's bath. Nothing of cable theory is put in. The small ions share one D, chosen so the
axoplasm has their 35.4 Ω·cm.

```js
import { Device, FARADAY, GAS_CONSTANT } from 'driftlet';
import { hodgkinHuxley, injector, pulse } from 'driftlet/kit';

// Hodgkin and Huxley's fibre: radius 238 µm, axoplasm 35.4 Ω·cm, 1 µF/cm², 18.5 °C
const T = 273.15 + 18.5, RT = GAS_CONSTANT * T, a = 238e-6, am = 2 / a, L = 0.06, h = 2e-4;
const z = { 'Na+': 1, 'K+': 1, 'Cl-': -1, 'A-': -1 };
const IN = { 'Na+': 50, 'K+': 400, 'Cl-': 52.5, 'A-': 397.5 };
// outside, set by their reversal potentials (rest at −65 mV; the leak carried by Cl⁻)
const E = { 'Na+': 0.05, 'K+': -0.077, 'Cl-': -0.054387 };
const OUT = Object.fromEntries(['Na+', 'K+', 'Cl-'].map((s) => [s, IN[s] * Math.exp((z[s] * E[s] * FARADAY) / RT)]));
// the small ions share one D, set so the axoplasm's resistivity is 35.4 Ω·cm
const DA = 1e-11, D = (1 / 0.354 / ((FARADAY * FARADAY) / RT) - DA * IN['A-']) / (IN['Na+'] + IN['K+'] + IN['Cl-']);
const hh = hodgkinHuxley({ T, area: am });
const ion = (Dv) => ({ D: Dv, mu0: 0, cRef: 1000 });
const sealed = { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' };
const axon = new Device({
  T,
  species: Object.entries(z).map(([name, zz]) => ({ name, z: zz })),
  materials: { axoplasm: { epsr: 0, species: { 'Na+': ion(D), 'K+': ion(D), 'Cl-': ion(D), 'A-': ion(DA) } } },
  regions: [{ name: 'axon', material: 'axoplasm', length: L, c0: IN }],
  contacts: { left: sealed, right: sealed },
  ports: [
    { name: 'membrane', region: 'axon', V: 0, area: am, capacitance: { C: 0.01, zeroCharge: 0 }, gates: hh.gates, species: hh.species, bath: { c: OUT } },
    injector({ name: 'stim', region: 'axon', from: 0, to: 1e-3, species: 'K+', I: pulse({ start: 30e-3, width: 0.2e-3, amplitude: 30 }) }),
  ],
  grid: { hmin: h, hmax: h },
});
axon.advance(30e-3, { tol: 1e-4 }); // to rest
const run = axon.advance(34e-3, { tol: 1e-4, dtMax: 5e-6, probes: [{ x: 0.02, quantity: 'phi' }, { x: 0.04, quantity: 'phi' }] });
const { t, probes } = run.trace;
const cross = (v) => { for (let i = 1; i < t.length; i++) if (v[i - 1] < 0 && v[i] >= 0) return t[i - 1] + ((t[i] - t[i - 1]) * -v[i - 1]) / (v[i] - v[i - 1]); };
const v = 0.02 / (cross(probes[1]) - cross(probes[0])); // between 2 and 4 cm
console.log(`${v.toFixed(1)} m/s: ${(v / 18.8).toFixed(3)} of their computed 18.8, ${(v / 21.2).toFixed(3)} of the measured 21.2`); // => 18.7, 0.996, 0.883
```

**driftlet lands on their model, 0.4% below its 18.8 m/s, and 12% below the measurement.** The
first is a check of driftlet: the cable equation comes out of ions drifting along the axon, and
test/cable.test.js holds it to an independent integration of the cable equation, to 0.2%. The
second is Hodgkin and Huxley's model against their own fibre, a gap they reported beside the
computed value.
