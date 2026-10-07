import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { hodgkinHuxley, injector, pulse } from '../src/kit.js';

// A propagating action potential. Along a squid giant axon (x runs along it), the axoplasm is a
// strictly neutral solution whose ions drift in the field along x, and the membrane is a port
// through the whole region: a capacitance (1 µF/cm², 2/a of it per volume of a radius-a axon) and
// Hodgkin and Huxley's channels as conductance links to the sea water outside, gated by m³h and n⁴.
// Nothing of cable theory is put in; with the outside held (no external resistance), it should
// come out. Hodgkin and Huxley's own case (1952): a = 238 µm, R_i = 35.4 Ω·cm, 18.5 °C, where they
// computed 18.8 m/s.

const T = 273.15 + 18.5, RT = GAS_CONSTANT * T, VT = RT / FARADAY, mV = 1e-3;
const a = 238e-6, am = 2 / a, Cm = 0.01, L = 0.06, h = 2e-4;
const z = { 'Na+': 1, 'K+': 1, 'Cl-': -1, 'A-': -1 };
// Their reversal potentials (modern convention: rest −65 mV), the leak's carried by Cl⁻.
const E = { 'Na+': 50 * mV, 'K+': -77 * mV, 'Cl-': -54.387 * mV };
const IN = { 'Na+': 50, 'K+': 400, 'Cl-': 52.5, 'A-': 397.5 };
const OUT = Object.fromEntries(['Na+', 'K+', 'Cl-'].map((s) => [s, IN[s] * Math.exp((z[s] * E[s]) / VT)]));
const G = { 'Na+': 1200, 'K+': 360, 'Cl-': 3 }; // S/m² of membrane (120, 36, 0.3 mS/cm²)
// The small ions share one D, set so the axoplasm's resistivity is theirs, 35.4 Ω·cm (the
// impermeant anions barely move).
const DA = 1e-11;
const D = (1 / 0.354 / ((FARADAY * FARADAY) / RT) - DA * IN['A-']) / (IN['Na+'] + IN['K+'] + IN['Cl-']);
const gates = hodgkinHuxley({ inside: IN, outside: { ...OUT, 'A-': 1 }, T }).gates;

function axon() {
  const species = (Dv) => ({ D: Dv, mu0: 0, cRef: 1000 });
  // Each ion's level outside, against φ there (the port's V, held at 0): its chemical part.
  const offset = (s) => (VT / z[s]) * Math.log(OUT[s] / 1000);
  const sealed = { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' };
  return new Device({
    T,
    species: Object.entries(z).map(([name, zz]) => ({ name, z: zz })),
    materials: { axoplasm: { epsr: 0, species: { 'Na+': species(D), 'K+': species(D), 'Cl-': species(D), 'A-': species(DA) } } },
    regions: [{ name: 'axon', material: 'axoplasm', length: L, c0: IN }],
    contacts: { left: sealed, right: sealed },
    ports: [
      {
        name: 'membrane', region: 'axon', V: 0, terminal: 'K+', area: am, capacitance: { C: Cm, zeroCharge: 0 }, gates,
        species: {
          'Na+': { type: 'conductance', G: G['Na+'] * am, offset: offset('Na+'), gates: { m: 3, h: 1 } },
          'K+': { type: 'conductance', G: G['K+'] * am, offset: offset('K+'), gates: { n: 4 } },
          'Cl-': { type: 'conductance', G: G['Cl-'] * am, offset: offset('Cl-') },
        },
      },
      // the stimulus: K⁺ into the first millimetre, 0.2 ms at 30 A/m² of cross-section
      injector({ name: 'stim', region: 'axon', from: 0, to: 1e-3, species: 'K+', I: pulse({ start: 30e-3, width: 0.2e-3, amplitude: 30 }) }),
    ],
    grid: { hmin: h, hmax: h },
  });
}

const rate = (r, V) => {
  const y = (V - r.midpoint) / r.scale;
  if (r.type === 'exp') return r.rate * Math.exp(y);
  if (r.type === 'sigmoid') return r.rate / (1 + Math.exp(-y));
  return (r.rate * y) / -Math.expm1(-y);
};
const inf = (g, V) => rate(g.alpha, V) / (rate(g.alpha, V) + rate(g.beta, V));
const Iion = (V, m, hh, n) => G['Na+'] * m ** 3 * hh * (V - E['Na+']) + G['K+'] * n ** 4 * (V - E['K+']) + G['Cl-'] * (V - E['Cl-']);

// The cable equation itself, C V_t = (σ/a_m) V_xx − I_ion + injected, by the method of lines
// (explicit, second order) on the same grid, sealed ends: the time its spike crosses 0 mV at x.
function cable(rest) {
  const sigma = ((FARADAY * FARADAY) / RT) * (D * (IN['Na+'] + IN['K+'] + IN['Cl-']) + DA * IN['A-']);
  const Dv = sigma / (am * Cm), Nx = Math.round(L / h) + 1, dt = (0.4 * h * h) / Dv;
  let y = [new Float64Array(Nx).fill(rest), ...['m', 'h', 'n'].map((k) => new Float64Array(Nx).fill(inf(gates[k], rest)))];
  const f = (t, [V, m, hh, n]) => {
    const out = [0, 0, 0, 0].map(() => new Float64Array(Nx));
    for (let i = 0; i < Nx; i++) {
      const l = V[i > 0 ? i - 1 : 1], r = V[i < Nx - 1 ? i + 1 : Nx - 2];
      const inj = i * h <= 1e-3 && t >= 0.1e-3 && t < 0.3e-3 ? 30 / 1e-3 / am : 0;
      out[0][i] = (Dv * (l - 2 * V[i] + r)) / (h * h) + (inj - Iion(V[i], m[i], hh[i], n[i])) / Cm;
      ['m', 'h', 'n'].forEach((k, j) => (out[1 + j][i] = rate(gates[k].alpha, V[i]) * (1 - [m, hh, n][j][i]) - rate(gates[k].beta, V[i]) * [m, hh, n][j][i]));
    }
    return out;
  };
  const at = { 0.02: NaN, 0.04: NaN };
  for (let t = 0; t < 4e-3; t += dt) {
    const k1 = f(t, y), y1 = y.map((v, j) => v.map((x, i) => x + dt * k1[j][i])), k2 = f(t + dt, y1);
    const next = y.map((v, j) => v.map((x, i) => x + (dt / 2) * (k1[j][i] + k2[j][i])));
    for (const x of [0.02, 0.04]) {
      const i = Math.round(x / h), v0 = y[0][i], v1 = next[0][i];
      if (Number.isNaN(at[x]) && v0 < 0 && v1 >= 0) at[x] = t + (dt * -v0) / (v1 - v0);
    }
    y = next;
  }
  return 0.02 / (at[0.04] - at[0.02]);
}

test('a propagating action potential: the cable equation, out of ions drifting along the axon, at Hodgkin and Huxley\'s 18.8 m/s', () => {
  const d = axon();
  // From the concentrations given, the membrane settles to rest, in a damped oscillation over
  // ~20 ms (no pump: it has no steady state).
  const rest = d.advance(30e-3, { tol: 1e-4 }).ports[0], j = rest.Vm.length >> 1;
  const Vr = rest.Vm[j];
  // Rest is where their ionic current vanishes, the gates at their limits: −65.0 mV.
  assert.ok(Math.abs(Iion(Vr, inf(gates.m, Vr), inf(gates.h, Vr), inf(gates.n, Vr))) < 1e-3 && Math.abs(Vr + 65 * mV) < 0.01 * mV, `rest ${Vr}`);
  for (const k of ["m", "h", "n"]) assert.ok(Math.abs(rest.gates[k][j] - inf(gates[k], Vr)) < 1e-4, `${k} at rest: ${rest.gates[k][j]} vs ${inf(gates[k], Vr)}`);
  // The spike crosses 0 mV at 2 cm and then at 4 cm.
  const sol = d.advance(34e-3, { tol: 1e-4, dtMax: 5e-6, probes: [{ x: 0.02, quantity: 'phi' }, { x: 0.04, quantity: 'phi' }] });
  assert.ok(sol.converged);
  const { t, probes } = sol.trace;
  const cross = (v) => {
    for (let i = 1; i < t.length; i++) if (v[i - 1] < 0 && v[i] >= 0) return t[i - 1] + ((t[i] - t[i - 1]) * -v[i - 1]) / (v[i] - v[i - 1]);
    return NaN;
  };
  const v = 0.02 / (cross(probes[1]) - cross(probes[0]));
  const ref = cable(Vr);
  assert.ok(Math.abs(v / ref - 1) < 2e-3, `${v} m/s against the cable equation's ${ref}`);
  assert.ok(Math.abs(v / 18.8 - 1) < 0.01, `${v} m/s against Hodgkin and Huxley's 18.8`);
  // It overshoots, and the ions that carry it barely change the axoplasm (well under 0.1%).
  assert.ok(Math.max(...probes[0]) > 20 * mV);
  const c = d.solution().c;
  for (const s of ['Na+', 'K+']) assert.ok(c[s].every((x) => Math.abs(x / IN[s] - 1) < 1e-3), s);
});
