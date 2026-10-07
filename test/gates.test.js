import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer, bath, aqueous, injector, pulse, hodgkinHuxley } from '../src/kit.js';

// Voltage-gated channels: a membrane face whose permeabilities are scaled by gates with
// Hodgkin–Huxley kinetics, dx/dt = α(V)(1 − x) − β(V)x in the voltage across the face. The squid
// axon of Hodgkin and Huxley (1952), at their 6.3 °C, with V = φ_in − φ_out (inside on the right).

const T = 273.15 + 6.3, RT = GAS_CONSTANT * T, mV = 1e-3;
const IONS = ['Na+', 'K+', 'Cl-'];
const OUT = { 'Na+': 440, 'K+': 20, 'Cl-': 450, 'A-': 10 }, IN = { 'Na+': 50, 'K+': 400, 'Cl-': 52, 'A-': 398 };
const HH = {
  m: { alpha: { type: 'expLinear', rate: 1000, midpoint: -40 * mV, scale: 10 * mV }, beta: { type: 'exp', rate: 4000, midpoint: -65 * mV, scale: -18 * mV } },
  h: { alpha: { type: 'exp', rate: 70, midpoint: -65 * mV, scale: -20 * mV }, beta: { type: 'sigmoid', rate: 1000, midpoint: -35 * mV, scale: 10 * mV } },
  n: { alpha: { type: 'expLinear', rate: 100, midpoint: -55 * mV, scale: 10 * mV }, beta: { type: 'exp', rate: 125, midpoint: -65 * mV, scale: -80 * mV } },
};
// The rate laws written out, as HH did (V in volts, rates in 1/s).
const rate = (r, V) => {
  const y = (V - r.midpoint) / r.scale;
  if (r.type === 'exp') return r.rate * Math.exp(y);
  if (r.type === 'sigmoid') return r.rate / (1 + Math.exp(-y));
  return (r.rate * y) / -Math.expm1(-y);
};
// GHK current (outward, A/m²) per unit permeability.
const ghk = (z, ci, co, V) => {
  const u = (z * V * FARADAY) / RT;
  return (z * FARADAY * u * (ci - co * Math.exp(-u))) / -Math.expm1(-u);
};
// Permeabilities giving HH's chord conductances (Na 120, K 36, leak 0.3 mS/cm²; the leak as Cl⁻) at
// rest, −65 mV.
const chord = (g, z, ci, co) => (g * (-65 * mV - (RT / (z * FARADAY)) * Math.log(co / ci))) / ghk(z, ci, co, -65 * mV);
const P = { Na: chord(1200, 1, 50, 440), K: chord(360, 1, 400, 20), Cl: chord(3, -1, 52, 450) };

const water = () => {
  const w = aqueous(IONS, { epsr: 0 });
  w.species.push({ name: 'A-', z: -1, cRef: 1000 });
  w.materials.water.species['A-'] = { D: 0.5e-9, mu0: 0 };
  return w;
};
const membrane = (gates = HH) => ({
  phi: { type: 'capacitive', C: 0.01 },
  gates,
  species: {
    'Na+': { type: 'permeability', P: P.Na, gates: { m: 3, h: 1 } },
    'K+': { type: 'permeability', P: P.K, gates: { n: 4 } },
    'Cl-': { type: 'permeability', P: P.Cl },
    'A-': 'blocked',
  },
});

test('gates under voltage clamp: α/(α + β) in steady state, and relaxing at α + β after a step', () => {
  // The same solution both sides, held by Cl⁻-referenced baths: the membrane sees the clamp, but
  // for the tens of µV that its charge, held in the strictly neutral solutions' edge cells, shifts
  // their φ by. The gates follow the voltage across the face itself.
  const S = { 'Na+': 100, 'K+': 100, 'Cl-': 200, 'A-': 0 };
  delete S['A-'];
  const w = aqueous(IONS, { epsr: 0 });
  const tiny = { phi: { type: 'capacitive', C: 0.01 }, gates: HH, species: { 'K+': { type: 'permeability', P: 1e-12, gates: { n: 4 } }, 'Na+': 'blocked', 'Cl-': 'blocked' } };
  const clamp = (V) =>
    build({ T, library: [w], stack: [bath(S, 'Cl-'), layer('water', 1e-6, { c0: S }), tiny, layer('water', 1e-6, { c0: S }), bath(S, 'Cl-', V)], grid: { hmin: 50e-9, hmax: 200e-9 } });
  for (const V of [-90 * mV, -65 * mV, -40 * mV, 0, 30 * mV]) {
    const sol = new Device(clamp(V)).solve(), itf = sol.interfaces[0];
    assert.ok(sol.converged && Math.abs(itf.V - V) < 0.1 * mV, `V = ${V}: membrane at ${itf.V}`);
    for (const [k, g] of Object.entries(HH)) {
      const inf = rate(g.alpha, itf.V) / (rate(g.alpha, itf.V) + rate(g.beta, itf.V));
      assert.ok(Math.abs(itf.gates[k] - inf) < 1e-12, `${k} at ${V}: ${itf.gates[k]} vs ${inf}`);
    }
  }
  // A step from −65 to −10 mV at 1 ms: x = x∞ + (x₀ − x∞) e^{−(α+β)t} after it.
  const d = new Device(clamp(pulse({ start: 1e-3, width: 1, amplitude: 55 * mV, base: -65 * mV })));
  const V0 = d.solve().interfaces[0].V;
  for (const t of [1.2e-3, 1.5e-3, 2e-3, 4e-3, 8e-3]) {
    const sol = d.advance(t, { tol: 1e-6 }), itf = sol.interfaces[0], V1 = itf.V; // (charged in ~0.1 µs)
    for (const [k, g] of Object.entries(HH)) {
      const a0 = rate(g.alpha, V0), b0 = rate(g.beta, V0), a1 = rate(g.alpha, V1), b1 = rate(g.beta, V1);
      const x0 = a0 / (a0 + b0), x1 = a1 / (a1 + b1);
      const x = x1 + (x0 - x1) * Math.exp(-(a1 + b1) * (t - 1e-3));
      assert.ok(Math.abs(itf.gates[k] - x) < 2e-4, `${k} at ${t} s: ${itf.gates[k]} vs ${x}`);
    }
  }
});

test('an action potential: the electrodiffusion model against the space-clamped HH equations with GHK currents', () => {
  // Bath | 5 µm outside | membrane | 20 µm of cytoplasm | closed, started from the concentrations
  // given (no pump, so no steady rest: the cell would run down to Donnan over hours) and stimulated
  // at 20 ms with 10 µA/cm² of K⁺ for 1 ms, injected 15–20 µm in.
  const I0 = 0.1;
  const d = new Device(
    build({
      T,
      library: [water()],
      stack: [
        bath(OUT, 'Cl-'),
        layer('water', 5e-6, { name: 'outside', c0: OUT }),
        membrane(),
        layer('water', 20e-6, { name: 'cytoplasm', c0: IN }),
        { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' },
      ],
      ports: [injector({ region: 'cytoplasm', from: 15e-6, to: 20e-6, species: 'K+', I: pulse({ start: 20e-3, width: 1e-3, amplitude: I0 }) })],
      grid: { hmin: 50e-9, hmax: 1e-6 },
    }),
  );
  // The ODE: C dV/dt = I_inj − Σ I_GHK(V) at the bulk concentrations, with the same gates, by RK4.
  const dt = 1e-7, steps = Math.round(30e-3 / dt), ode = new Float64Array(steps + 1);
  const inf = (g, V) => rate(g.alpha, V) / (rate(g.alpha, V) + rate(g.beta, V));
  const f = (t, [V, m, h, n]) => {
    const I = P.Na * m ** 3 * h * ghk(1, IN['Na+'], OUT['Na+'], V) + P.K * n ** 4 * ghk(1, IN['K+'], OUT['K+'], V) + P.Cl * ghk(-1, IN['Cl-'], OUT['Cl-'], V);
    const x = (g, v) => rate(g.alpha, V) * (1 - v) - rate(g.beta, V) * v;
    return [((t >= 20e-3 && t < 21e-3 ? I0 : 0) - I) / 0.01, x(HH.m, m), x(HH.h, h), x(HH.n, n)];
  };
  let y = [1e-12, inf(HH.m, 0), inf(HH.h, 0), inf(HH.n, 0)]; // (V₀ = 0: GHK's 0/0 is avoided)
  const add = (a, b, s) => a.map((v, i) => v + s * b[i]);
  ode[0] = y[0];
  for (let k = 0; k < steps; k++) {
    const t = k * dt, k1 = f(t, y), k2 = f(t + dt / 2, add(y, k1, dt / 2)), k3 = f(t + dt / 2, add(y, k2, dt / 2)), k4 = f(t + dt, add(y, k3, dt));
    y = y.map((v, i) => v + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
    ode[k + 1] = y[0];
  }
  let worst = 0, peak = -1, peakAt = 0, odePeak = -1, odePeakAt = 0;
  const probes = [{ interface: 0, gate: 'm' }, { interface: 0, species: 'Na+' }];
  for (let k = 1; k <= 300; k++) {
    const t = k * 1e-4, sol = d.advance(t, { tol: 1e-5, probes });
    // (the interface probes read what the solution reports, at each step)
    const [m, N] = sol.trace.probes.map((v) => v.at(-1));
    assert.ok(m === sol.interfaces[0].gates.m && N === sol.interfaces[0].N['Na+']);
    assert.ok(sol.converged, `at ${t} s`);
    const V = sol.interfaces[0].V, Vo = ode[Math.round(t / dt)];
    worst = Math.max(worst, Math.abs(V - Vo));
    if (V > peak) [peak, peakAt] = [V, t];
    if (Vo > odePeak) [odePeak, odePeakAt] = [Vo, t];
  }
  // It fires (the peak 70 mV above rest), when and as high as the ODE says: they differ by what
  // the spike's ions do to the concentrations at the membrane (~0.4% of Na⁺ inside) and by the
  // time stepping.
  assert.ok(peak > 0 && peakAt > 21e-3 && peakAt < 24e-3, `peak ${peak} at ${peakAt}`);
  assert.ok(Math.abs(peakAt - odePeakAt) <= 1e-4 && Math.abs(peak - odePeak) < 1 * mV, `peak ${peak} at ${peakAt} vs ${odePeak} at ${odePeakAt}`);
  assert.ok(worst < 1 * mV, `largest difference ${worst * 1e3} mV`);
});

test('gates: what the definition must say', () => {
  const def = (itf) => ({
    T,
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { water: { epsr: 0, species: { 'K+': { D: 2e-9, mu0: 0, cRef: 1000 }, 'Cl-': { D: 2e-9, mu0: 0, cRef: 1000 } } } },
    regions: [{ material: 'water', length: 1e-6 }, { material: 'water', length: 1e-6 }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.01 }, ...itf }],
    contacts: { left: { V: 0, terminal: 'Cl-', species: { 'K+': { type: 'equilibrium', offset: 0 }, 'Cl-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, terminal: 'Cl-', species: { 'K+': { type: 'equilibrium', offset: 0 }, 'Cl-': 'equilibrium' }, phi: 'bulk' } },
  });
  const n = HH.n;
  assert.throws(() => new Device(def({ gates: { n }, species: { 'K+': { type: 'permeability', P: 1e-8, gates: { q: 4 } } } })), /gates\.q: no such gate on this face \(it has n\)/);
  assert.throws(() => new Device(def({ gates: { n }, species: { 'K+': { type: 'permeability', P: 1e-8, gates: { n: 2.5 } } } })), /exponent must be a whole number/);
  assert.throws(() => new Device(def({ gates: { n: { alpha: { ...n.alpha, type: 'linear' }, beta: n.beta } } })), /alpha must be \{ type, rate, midpoint, scale \}/);
  assert.throws(() => new Device(def({ gates: { n: { alpha: { ...n.alpha, scale: 0 }, beta: n.beta } } })), /scale \(V\) must not be zero/);
  assert.throws(() => new Device(def({ gates: { n: { alpha: n.alpha } } })), DeviceError);
  // A gate that scales nothing still runs (and is reported).
  const sol = new Device(def({ gates: { n } })).solve();
  assert.ok(sol.converged && Math.abs(sol.interfaces[0].gates.n - rate(n.alpha, 0) / (rate(n.alpha, 0) + rate(n.beta, 0))) < 1e-12);
});

test("the kit's hodgkinHuxley(): their gates, Q₁₀ = 3, and GHK permeabilities with their chord conductances at rest", () => {
  const hh = hodgkinHuxley({ inside: IN, outside: OUT });
  assert.deepEqual(hh.gates, HH);
  for (const [sp, p] of [['Na+', P.Na], ['K+', P.K], ['Cl-', P.Cl]]) assert.ok(Math.abs(hh.species[sp].P / p - 1) < 1e-12, `${sp}: ${hh.species[sp].P} vs ${p}`);
  assert.deepEqual([hh.species['Na+'].gates, hh.species['K+'].gates, hh.species['Cl-'].gates], [{ m: 3, h: 1 }, { n: 4 }, undefined]);
  // Each one's GHK current at −65 mV is g (V − E) for its g.
  const V = -65 * mV;
  for (const [sp, z, g] of [['Na+', 1, 1200], ['K+', 1, 360], ['Cl-', -1, 3]]) {
    const E = (RT / (z * FARADAY)) * Math.log(OUT[sp] / IN[sp]);
    assert.ok(Math.abs((hh.species[sp].P * ghk(z, IN[sp], OUT[sp], V)) / (g * (V - E)) - 1) < 1e-12, sp);
  }
  // Ten degrees warmer, three times as fast.
  const warm = hodgkinHuxley({ inside: IN, outside: OUT, T: T + 10 });
  for (const k of ['m', 'h', 'n']) for (const ab of ['alpha', 'beta']) assert.ok(Math.abs(warm.gates[k][ab].rate / HH[k][ab].rate - 3) < 1e-12);
  assert.throws(() => hodgkinHuxley({ inside: { 'Na+': 50 }, outside: OUT }), /give K\+'s concentration inside and outside/);
});
