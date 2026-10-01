import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, SolverError } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// Complex helpers on [re, im] pairs.
const add = (a, b) => [a[0] + b[0], a[1] + b[1]];
const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const div = (a, b) => {
  const q = b[0] ** 2 + b[1] ** 2;
  return [(a[0] * b[0] + a[1] * b[1]) / q, (a[1] * b[0] - a[0] * b[1]) / q];
};
const sqrt = ([x, y]) => {
  const r = Math.hypot(x, y);
  return [Math.sqrt((r + x) / 2), Math.sign(y || 1) * Math.sqrt((r - x) / 2)];
};
const tanh = ([x, y]) => {
  if (x > 20) return [1, 0];
  const d = Math.cosh(2 * x) + Math.cos(2 * y);
  return [Math.sinh(2 * x) / d, Math.sin(2 * y) / d];
};
const relErr = (Z, q, exact) => Math.hypot(Z.Z.re[q] - exact[0], Z.Z.im[q] - exact[1]) / Math.hypot(...exact);

test('blocking electrodes: the Macdonald impedance of a symmetric electrolyte, 100 Hz to 1 GHz', () => {
  // Equal diffusivities, walls pinned (no Stern layer). Linear PNP gives
  //   Z = 2d/(εDk²) + 2 tanh(kd)/(iωελ²k³),   k² = (1 + iωλ²/D)/λ²,
  // bulk resistance in series with the diffuse layers, turning into the geometric capacitance.
  const D = 1.5e-9, c0 = 1, eps = 78.5 * EPS0, L = 1e-6, d = L / 2;
  const lam = Math.sqrt((eps * RT) / (2 * FARADAY ** 2 * c0));
  const wall = { V: 0, phi: 'pinned', zeroCharge: 0 };
  const dev = new Device({
    species: [
      { name: 'Na+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 78.5, species: { 'Na+': { D, mu0: -261.9e3 }, 'Cl-': { D, mu0: -131.2e3 } } } },
    regions: [{ material: 'water', length: L, c0: { 'Na+': c0, 'Cl-': c0 } }],
    contacts: { left: wall, right: wall },
    grid: { hmin: lam / 200, hmax: lam / 5, ratio: 1.05 },
  });
  const fs = [1e2, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9];
  const Z = dev.impedance(fs);
  fs.forEach((f, q) => {
    const w = 2 * Math.PI * f;
    const k = div(sqrt([1, (w * lam * lam) / D]), [lam, 0]);
    const k2 = mul(k, k), k3 = mul(k2, k);
    const exact = add(div([2 * d, 0], mul([eps * D, 0], k2)), mul(div([2, 0], mul([0, w * eps * lam * lam], k3)), tanh(mul(k, [d, 0]))));
    assert.ok(relErr(Z, q, exact) < 3e-4, `f = ${f}: ${relErr(Z, q, exact)}`);
  });
});

test('reversible electrodes: finite-length Warburg of Ag | AgNO₃ | Ag, in voltage and current modes', () => {
  // ε = 0 electrolyte, Ag⁺ reversible at both ends, NO₃⁻ blocked. Linearising,
  //   Z = (RT/F²c)[L/(D₊+D₋) + 2D₋/(D₊(D₊+D₋)) · tanh(qL/2)/q],  q = √(iω/D_s),
  // with D_s = 2D₊D₋/(D₊+D₋): the bulk resistance plus a salt-diffusion Warburg element.
  const Dp = 1.65e-9, Dm = 1.9e-9, c0 = 10, L = 20e-6;
  const el = { terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
  const cell = (right) =>
    new Device({
      species: [
        { name: 'Ag+', z: 1, cRef: 1000 },
        { name: 'NO3-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'Ag+': { D: Dp, mu0: 77.1e3 }, 'NO3-': { D: Dm, mu0: -111.3e3 } } } },
      regions: [{ material: 'water', length: L, c0: { 'NO3-': c0 } }],
      contacts: { left: { ...el, V: 0 }, right },
      grid: { minCells: 400 },
    });
  const Ds = (2 * Dp * Dm) / (Dp + Dm);
  const exact = (f) => {
    const q = sqrt([0, (2 * Math.PI * f) / Ds]);
    const tq = div(tanh(mul(q, [L / 2, 0])), q);
    const pre = VT / (FARADAY * c0), RD = (2 * Dm) / (Dp * (Dp + Dm));
    return [pre * (L / (Dp + Dm) + RD * tq[0]), pre * RD * tq[1]];
  };
  const fs = [1e-3, 0.1, 1, 10, 100];
  for (const right of [{ ...el, V: 0 }, { ...el, I: 0 }]) {
    const Z = cell(right).impedance(fs, { profiles: true });
    fs.forEach((f, q) => assert.ok(relErr(Z, q, exact(f)) < 3e-5, `${right.I === undefined ? 'held' : 'driven'}, f = ${f}: ${relErr(Z, q, exact(f))}`));
    // Profiles: the electrolyte stays neutral, and at low frequency the salt profile is linear.
    const p = Z.profiles[0];
    for (let g = 0; g < p.c['Ag+'].re.length; g++) assert.ok(Math.abs(p.c['Ag+'].re[g] - p.c['NO3-'].re[g]) < 1e-9 * Math.abs(p.c['Ag+'].re[0]));
    const mid = p.c['Ag+'].re.length >> 1;
    assert.ok(Math.abs(p.c['Ag+'].re[mid]) < 1e-6 * Math.abs(p.c['Ag+'].re[0]));
  }
  assert.throws(() => cell({ ...el, V: 0, R: 1 }).impedance([1]), (e) => e instanceof SolverError && /series resistance/.test(e.message));
});

test('low-frequency limit equals the steady differential resistance (Fermi–Dirac pn diode with recombination)', () => {
  const Nc = 2.8e25 / 6.02214076e23, Nv = 1.04e25 / 6.02214076e23; // mol/m³
  const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
  const dev = new Device({
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: {
      Si: {
        epsr: 11.7,
        species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: 1.12 * FARADAY, cRef: Nv } },
        statistics: [{ type: 'fermi-dirac', species: ['e-', 'h+'] }],
      },
    },
    regions: [
      { material: 'Si', length: 2e-6, fixedCharge: 1e25 / 6.02214076e23 * FARADAY },
      { material: 'Si', length: 2e-6, fixedCharge: -1e22 / 6.02214076e23 * FARADAY },
    ],
    contacts: { left: ohmic(0), right: ohmic(0.45) },
    bulkReactions: [{ nu: { 'e-': -1, 'h+': -1 }, kf: { Si: 1e-6 } }],
    grid: { hmin: 0.5e-9, hmax: 20e-9 },
  });
  const Z = dev.impedance([1e-3, 1e6]);
  const I = (V) => {
    dev.set({ contacts: { right: { V } } });
    return dev.solve().current;
  };
  const dV = 1e-5;
  const R = -(2 * dV) / (I(0.45 + dV) - I(0.45 - dV));
  assert.ok(Math.abs(Z.Z.re[0] / R - 1) < 1e-6, `${Z.Z.re[0]} vs ${R}`);
  assert.ok(Math.abs(Z.Z.im[0]) < 1e-6 * R);
  assert.ok(Z.Z.im[1] < 0, 'capacitive (depletion and diffusion charge) at high frequency');
});
