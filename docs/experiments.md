# Against experiment

The validation suite holds driftlet to analytic limits: Cottrell, Debye, Shockley, cable theory. This
page holds it to measurements. Each entry takes a classic measurement with its published numbers,
models it the way a user would, with the library's own data where it has them, and says how close
driftlet lands and why it misses where it does. A miss here is information, not a failure: it marks
where the physics driftlet solves stops being the physics of the experiment.

Every example below runs as part of the test suite, and the numbers in its `// =>` comments are
what it prints.

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
