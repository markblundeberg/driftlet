import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer, bath, aqueous, check } from '../src/kit.js';

// A charged pore between two KCl reservoirs (an ion channel, a nanopore), resolved: its walls'
// charge as a fixed charge density X, ε > 0 throughout, so the Donnan layers at each mouth and
// the depletion beside them come out of Poisson's equation. In a long pore the bulk inside is
// neutral, and the Teorell–Meyer–Sievers theory holds: Donnan equilibrium at each mouth, so
// K⁺ = (X + √(X² + 4c²))/2 and Cl⁻ = K⁺ − X inside, and diffusion between. The reservoirs widen
// away from the mouth in 3D, and add almost no resistance; in this 1D slice they're layers as wide
// as the pore, so their D is raised to stand in for the widening (no equilibrium depends on D).
// The dilute, uncharged end (a constant field, Goldman–Hodgkin–Katz) is test/membrane's resolved
// lipid layer.

const T = 298.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;
const D = { 'K+': 1.96e-10, 'Cl-': 2.03e-10 }; // in the pore: a tenth of water's
const library = aqueous(['K+', 'Cl-'], { epsr: 78.5 });
for (const sp of Object.values(library.materials.water.species)) sp.D *= 1e4; // the widening reservoirs
library.materials.pore = { epsr: 78.5, species: Object.fromEntries(['K+', 'Cl-'].map((s) => [s, { D: D[s], mu0: library.materials.water.species[s].mu0 }])) };

const salt = (c) => ({ 'K+': c, 'Cl-': c });
const access = (c) => Math.max(10e-9, (8 * 0.304e-9) / Math.sqrt(c / 1000)); // eight Debye lengths
// The pore in two halves, each with its own charge (mol/m³ of fixed negative charge).
const pore = ({ cl, cr = cl, X1, X2 = X1, L, V = { I: 0 } }) =>
  build({
    library: [library],
    stack: [
      bath(salt(cl), 'Cl-'),
      layer('water', access(cl), { name: 'left' }),
      { dipole: 0 },
      layer('pore', L / 2, { fixedCharge: -X1 * FARADAY }),
      layer('pore', L / 2, { fixedCharge: -X2 * FARADAY }),
      { dipole: 0 },
      layer('water', access(cr), { name: 'right' }),
      bath(salt(cr), 'Cl-', V),
    ],
    grid: { hmin: 0.005e-9, hmax: 0.5e-9, ratio: 1.1 },
  });
const inside = (c, X) => {
  const k = (X + Math.sqrt(X * X + 4 * c * c)) / 2;
  return { 'K+': k, 'Cl-': k - X };
};
// TMS conductance per area at small bias (the Donnan layers stay in equilibrium).
const tmsG = (c, X, L) => {
  const n = inside(c, X);
  return ((FARADAY * FARADAY) / (RT * L)) * (D['K+'] * n['K+'] + D['Cl-'] * n['Cl-']);
};
const G = (c, X, L) => {
  const s = new Device(pore({ cl: c, X1: X, L, V: 1e-4 })).solve();
  assert.ok(s.converged);
  return -s.current / 1e-4; // the right reservoir 0.1 mV up drives current toward −x
};

test('a charged pore conducts as Teorell–Meyer–Sievers says, down to a plateau set by its charge at low salt, short pores less by their end layers (as 1/L)', () => {
  const X = 1000; // 1 M of fixed negative charge
  for (const c of [1, 10, 100, 1000]) {
    const r = G(c, X, 320e-9) / tmsG(c, X, 320e-9);
    assert.ok(Math.abs(r - 1) < 0.01, `${c} mM: ${r} of TMS`);
  }
  // The plateau: a thousandfold less salt, the conductance hardly falls (bulk KCl's would fall
  // a thousandfold): the counter-ions are the wall's, not the salt's.
  assert.ok(G(1, X, 320e-9) / G(1000, X, 320e-9) > 0.4);
  // Shorter pores fall short of TMS by the layers at each end, over a fixed length: the shortfall
  // times L is the same.
  const short = [20e-9, 80e-9, 320e-9].map((L) => (1 - G(30, X, L) / tmsG(30, X, L)) * L);
  assert.ok(short.every((d) => Math.abs(d / short[2] - 1) < 0.1), short.join(', '));
});

test('across a salt gradient, the zero-current potential is TMS: towards a K⁺ electrode as the charge grows', () => {
  // 100 mM | pore | 10 mM: the voltage between Ag/AgCl electrodes in the two reservoirs (each
  // reads Cl⁻'s level), 2t₊V_T ln 10 for a perfectly K⁺-selective pore.
  const Dp = D['K+'], Dm = D['Cl-'];
  const tms = (c1, c2, X) => {
    const s = (c) => (Math.sqrt(X * X + 4 * c * c) - X) / 2, s1 = s(c1), s2 = s(c2);
    return VT * (((Dm - Dp) / (Dp + Dm)) * Math.log(((Dp + Dm) * s2 + Dp * X) / ((Dp + Dm) * s1 + Dp * X)) - Math.log(s2 / s1));
  };
  for (const [X, L, tol] of [[100, 100e-9, 5e-3], [300, 400e-9, 5e-4], [1000, 100e-9, 1e-4]]) {
    const dev = new Device(pore({ cl: 100, cr: 10, X1: X, L })), s = dev.solve();
    const E = s.terminalVoltage;
    assert.ok(Math.abs(E / tms(100, 10, X) - 1) < tol, `X ${X}: ${E} vs TMS ${tms(100, 10, X)}`);
    assert.ok(check(dev, s, { refine: false }).ok);
  }
  // With 1 M of charge K⁺ is in equilibrium across the pore (its level flat), Cl⁻ all but shut out.
  const s = new Device(pore({ cl: 100, cr: 10, X1: 1000, L: 100e-9 })).solve(), g = s.x.length - 1;
  assert.ok(Math.abs(s.V['K+'][g] - s.V['K+'][0]) < 2e-3);
  assert.ok(Math.abs(s.terminalVoltage / (2 * VT * Math.log(10)) - 1) < 0.01);
});

test('a symmetric pore conducts the same both ways; charge on one side rectifies, and a bipolar pore is a diode', () => {
  const I = (X1, X2, V) => new Device(pore({ cl: 100, X1, X2, L: 20e-9, V })).solve().current;
  const sym = [I(1000, 1000, 0.1), I(1000, 1000, -0.1)];
  assert.ok(Math.abs(sym[0] + sym[1]) < 1e-9 * Math.abs(sym[0]));
  // Current toward −x (positive bias on the right) through charge on the left half only, or
  // − then + (cation-selective on the left, anion-selective on the right): one way the pore
  // fills with salt, the other it empties.
  const half = Math.abs(I(1000, 0, -0.1) / I(1000, 0, 0.1)), bipolar = Math.abs(I(1000, -1000, -0.1) / I(1000, -1000, 0.1));
  assert.ok(half > 2 && bipolar > 10 && bipolar > half, `${half}, ${bipolar}`);
});
