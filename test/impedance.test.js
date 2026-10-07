import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, SolverError, units } from '../src/index.js';
import { build, layer, ohmic, bath, aqueous, metal, half, semiconductor } from '../src/kit.js';

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

test('an abrupt junction without recombination: Re Y is the steady dI/dV, and C the junction charge dQ/dV', () => {
  // Heavy doping on both sides: what passes is ~1e8 minority carriers per m³ reaching the far
  // contacts, a conductance some 1e17 times the terms of the contact's current, which the
  // impedance reads instead across the depletion region (the same total current).
  const NA = 6.02214076e23;
  const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
  const def = {
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: 2.8e25 / NA }, 'h+': { D: 12e-4, mu0: 1.12 * FARADAY, cRef: 1.04e25 / NA } } } },
    regions: [
      { material: 'Si', length: 1e-6, fixedCharge: (1e24 / NA) * FARADAY },
      { material: 'Si', length: 1e-6, fixedCharge: (-1e24 / NA) * FARADAY },
    ],
    contacts: { left: ohmic(0), right: ohmic(0) },
  };
  const fs = [1e-6, 1e2];
  const Z = new Device(def).impedance(fs).Z;
  // Steady states either side: the current, and the charge left of the junction (its D).
  const at = (V) => {
    const d = new Device(def);
    d.set({ contacts: { right: { V } } });
    const s = d.solve();
    return [s.current, s.interfaces[0].D];
  };
  const dV = 1e-4, [Ip, Dp] = at(dV), [Im, Dm] = at(-dV);
  const G = -(Ip - Im) / (2 * dV), C = -(Dp - Dm) / (2 * dV);
  fs.forEach((f, q) => {
    const m = Z.re[q] ** 2 + Z.im[q] ** 2, Yr = Z.re[q] / m, Cq = -Z.im[q] / m / (2 * Math.PI * f);
    if (q === 0) assert.ok(Math.abs(Yr / G - 1) < 1e-4, `${f} Hz: Re Y ${Yr} vs ${G} S/m²`);
    assert.ok(Math.abs(Cq / C - 1) < 1e-6, `${f} Hz: C ${Cq} vs ${C} F/m²`);
  });
});

test('a redox electrode at its open circuit against a bath: the DC limit is the steady dI/dV', () => {
  // Fe³⁺/Fe²⁺ on platinum, double layer resolved, a bath 10 µm away. The response to the
  // terminal is nearly a uniform shift of every level in the solution, and the current is their
  // slope, uniform to ~1e-14 of the response: solved as a shift of each region plus the rest.
  const c = { 'Fe3+': 9.6, 'Fe2+': 0.536, 'Cl-': 29.872 };
  const def = build({
    T: 300,
    library: [aqueous(['Fe3+', 'Fe2+', 'Cl-'], { epsr: 78.3 }), metal('Pt')],
    stack: [
      ohmic(0, ['e-']),
      layer('Pt', 1e-6),
      { phi: { type: 'capacitive', C: 0.395 }, zeroCharge: -0.119, reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 0.0411, alpha: 0.377 }] },
      layer('water', 9.83e-6, { c0: { ...c } }),
      bath({ ...c }, 'Cl-'),
    ],
    grid: { hmin: 2.6e-10, hmax: 3.7e-7, ratio: 1.136 },
  });
  const oc = new Device(def);
  oc.set({ contacts: { right: { I: 0 } } });
  const V0 = oc.solve().terminals.right.V;
  const at = (V) => {
    const d = new Device(def);
    d.set({ contacts: { right: { V } } });
    return d;
  };
  const Z = at(V0).impedance([1e-8, 1e-4]).Z;
  const I = (V) => at(V).solve().terminals.right.current, dV = 1e-5, G = (I(V0 + dV) - I(V0 - dV)) / (2 * dV);
  for (let q = 0; q < 2; q++) {
    const Y = Z.re[q] / (Z.re[q] ** 2 + Z.im[q] ** 2);
    assert.ok(Math.abs(Y / G - 1) < 1e-4, `Re Y ${Y} vs dI/dV ${G} S/m²`);
  }
});

test('a MOS capacitor with a metal gate and its channel held by a port: the gate is dQ/dV, the back contact passive', () => {
  // The gate's current, read at its contact, is a metal's σ/h (~1e15 S/m²) times a difference
  // of levels equal to round-off: it's read across the oxide instead, which the port's window
  // doesn't reach. The back contact moves the whole silicon nearly uniformly: its current is
  // resolved by splitting off that shift, the port's window kept with the edges.
  const zc = -0.95;
  const def = (V) => ({
    ...build({
      T: 300,
      library: [semiconductor('Si'), metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }],
      stack: [ohmic(V, ['e-']), layer('Al', 100e-9), { phi: { type: 'capacitive', C: 100 }, zeroCharge: zc }, layer('SiO2', 10e-9), { dipole: 0 }, layer('Si', 2e-6, { name: 'Si', acceptors: units.perCm3(1e17) }), ohmic(0)],
      grid: { hmin: 0.25e-9, hmax: 100e-9, ratio: 1.1 },
    }),
    ports: [{ name: 'channel', region: 'Si', from: 0, to: 1e-9, V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' } }],
  });
  const f = [1e-2, 1, 100, 1e4];
  for (const V of [-1.5, 0, 1.5]) {
    const dev = new Device(def(V));
    const gate = dev.impedance(f, { terminal: 'left' }).Z, back = dev.impedance(f, { terminal: 'right' }).Z;
    const dV = 1e-4, Q = (v) => new Device(def(v)).solve().interfaces[0].D;
    const C = (Q(V + dV) - Q(V - dV)) / (2 * dV);
    f.forEach((fq, j) => {
      const mag = Math.hypot(gate.re[j], gate.im[j]);
      assert.ok(Math.abs(gate.re[j]) < 1e-6 * mag, `V = ${V}, ${fq} Hz: gate Re Z ${gate.re[j]} of |Z| ${mag}`);
      if (fq <= 1) assert.ok(Math.abs(-1 / (2 * Math.PI * fq * gate.im[j]) / C - 1) < 1e-3, `V = ${V}, ${fq} Hz: C ${-1 / (2 * Math.PI * fq * gate.im[j])} against dQ/dV ${C}`);
      assert.ok(back.re[j] > -1e-3 * Math.hypot(back.re[j], back.im[j]), `V = ${V}, ${fq} Hz: back contact Re Z ${back.re[j]}`);
    });
  }
});
