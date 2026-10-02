import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, SEMICONDUCTORS, pulse, square, triangle, ramp, injector, recombination } from '../src/kit.js';

// Sources as plain data (waveforms, an injector port, recombination from a lifetime), and the
// Haynes–Shockley experiment built from them, read with probes, against its analytic solution.

const area = ({ t, values }) => t.reduce((a, tk, k) => (k ? a + ((values[k] + values[k - 1]) / 2) * (tk - t[k - 1]) : 0), 0);

test('waveforms: a pulse keeps its area whatever its rise; square, triangle and ramp', () => {
  assert.deepEqual(pulse({ start: 1e-6, width: 2e-6, amplitude: 5 }), { t: [1e-6, 1e-6, 3e-6, 3e-6], values: [0, 5, 5, 0] });
  for (const rise of [0, 1e-9, 0.5e-6]) assert.ok(Math.abs(area(pulse({ start: 1e-6, width: 2e-6, amplitude: 5, rise })) - 1e-5) < 1e-18);
  assert.throws(() => pulse({ width: 1e-6, amplitude: 1 }), /start after t = 0/);
  assert.throws(() => pulse({ start: 1, width: 1e-6, amplitude: 1, rise: 2e-6 }), DeviceError);
  assert.deepEqual(square({ high: 1, low: -1, period: 2 }), { t: [0, 1, 1, 2], values: [1, 1, -1, -1], repeat: true });
  assert.deepEqual(triangle({ from: 0.2, to: 0.8, period: 4 }), { t: [0, 2, 4], values: [0.2, 0.8, 0.2], repeat: true });
  assert.deepEqual(ramp({ from: 0, to: 1, duration: 2, start: 1 }), { t: [1, 3], values: [0, 1] });
});

test('recombination from a lifetime, and units back out of SI', () => {
  const n0 = units.perCm3(1e16);
  assert.deepEqual(recombination({ material: 'Si', tau: 1e-6, majority: n0 }), { equation: 'e- + h+ = 0', kf: { Si: 1 / (1e-6 * n0) } });
  assert.throws(() => recombination({ material: 'Si', tau: 1e-6 }), /majority/);
  for (const [to, from, v] of [['toPerCm3', 'perCm3', 3e15], ['toMolar', 'molar', 0.1], ['toNm', 'nm', 7], ['toUm', 'um', 7], ['toCm2PerS', 'cm2PerS', 12], ['toEV', 'eV', 1.1]]) {
    assert.ok(Math.abs(units[to](units[from](v)) / v - 1) < 1e-14, to);
  }
});

test('Haynes–Shockley: an injected hole pulse drifts, spreads and recombines at the ambipolar rates, and a probe sees it pass', () => {
  // n-Ge bar under 10 V/cm. An injector puts in 0.1% of n₀ over 50 µm for 0.5 µs; the holes
  // drift, spread and recombine (τ = 20 µs). At low injection the excess moves as one
  // quasi-neutral packet, whatever its shape: between two times its centre moves at
  // μ* = (n₀ − p₀)μₙμₚ/(n₀μₙ + p₀μₚ) times E, its variance grows by 2D*Δt with
  // D* = (n₀ + p₀)DₙDₚ/(n₀Dₙ + p₀Dₚ), and its amount decays at k_f(n₀ + p₀). (In germanium p₀/n₀ is
  // 6e-4, so μ* is 0.9991 μₚ, which the simulation resolves.)
  const ND = units.perCm3(1e15), L = 3e-3, E = 1000, tp = 0.5e-6, W = 50e-6, xe = 0.5e-3, tau = 20e-6, T = 300;
  const I = (0.001 * ND * FARADAY * W) / tp;
  const reaction = recombination({ material: 'Ge', tau, majority: ND });
  const dev = new Device(
    build({
      T,
      library: [semiconductor('Ge')],
      stack: [ohmic(E * L), layer('Ge', L, { name: 'bar', donors: ND }), ohmic(0)],
      bulkReactions: [reaction],
      ports: [injector({ region: 'bar', from: xe - W / 2, to: xe + W / 2, species: 'h+', I: pulse({ width: tp, amplitude: I, rise: 1e-9 }) })],
      grid: { hmin: 5e-6, hmax: 5e-6 },
    }),
  );
  const rest = dev.solve();
  const moments = (s) => {
    const m = [0, 0, 0];
    for (let g = 1; g < s.x.length; g++) {
      for (const k of [g - 1, g]) {
        const dp = (s.c['h+'][k] - rest.c['h+'][k]) * ((s.x[g] - s.x[g - 1]) / 2);
        m[0] += dp;
        m[1] += dp * s.x[k];
        m[2] += dp * s.x[k] ** 2;
      }
    }
    return { amount: m[0], mean: m[1] / m[0], variance: m[2] / m[0] - (m[1] / m[0]) ** 2 };
  };
  const d = 1e-3, t1 = 2e-6, t2 = 7e-6; // the probe 1 mm downstream; t2 before the pulse's tail reaches the right contact
  // dtMax: a probe reads at accepted steps, which grow to ~200 ns here; 10 ns resolves the passing pulse.
  const s1 = dev.advance(t1, { dtMax: 10e-9, probes: [{ x: xe + d, species: 'h+' }] });
  const a = moments(s1);
  const s2 = dev.advance(t2, { dtMax: 10e-9, probes: [{ x: xe + d, species: 'h+' }] });
  const b = moments(s2);
  assert.ok(s1.converged && s2.converged);

  const n0 = rest.c['e-'][rest.x.length >> 1], p0 = rest.c['h+'][rest.x.length >> 1];
  const VT = (GAS_CONSTANT * T) / FARADAY, mun = SEMICONDUCTORS.Ge.mun * 1e-4, mup = SEMICONDUCTORS.Ge.mup * 1e-4;
  const muStar = ((n0 - p0) * mun * mup) / (n0 * mun + p0 * mup), DStar = ((n0 + p0) * mun * mup * VT) / (n0 * mun + p0 * mup);
  const drift = (b.mean - a.mean) / (muStar * E * (t2 - t1)), spread = (b.variance - a.variance) / (2 * DStar * (t2 - t1));
  const decay = Math.log(a.amount / b.amount) / (reaction.kf.Ge * (n0 + p0) * (t2 - t1));
  assert.ok(Math.abs(drift - 1) < 2e-4, `drift: ${drift} of μ*E`);
  assert.ok(Math.abs(spread - 1) < 1e-2, `spread: ${spread} of 2D*`);
  assert.ok(Math.abs(decay - 1) < 2e-3, `decay: ${decay} of k_f(n₀ + p₀)`);
  // Nothing was lost while it went in: the amount at t1 is what was injected, each hole decaying
  // from when it went in.
  const survived = (tau / tp) * (Math.exp(-(t1 - tp) / tau) - Math.exp(-t1 / tau));
  assert.ok(Math.abs(a.amount / (((I * tp) / FARADAY) * survived) - 1) < 2e-3, `amount: ${a.amount / (((I * tp) / FARADAY) * survived)} of what went in`);

  // The probe sees it pass when a drifting, spreading, decaying Gaussian would (to within that
  // picture's own approximation: the injection is a box, not a point).
  const theory = (tt) => {
    const tc = tt - tp / 2, varx = 2 * DStar * tc + W * W / 12 + (muStar * E * tp) ** 2 / 12;
    return Math.exp(-((d - muStar * E * tc) ** 2) / (2 * varx) - tc / tau) / Math.sqrt(varx);
  };
  let tTheory = 0, best = 0;
  for (let tt = tp; tt < t2; tt += 1e-10) if (theory(tt) > best) [best, tTheory] = [theory(tt), tt];
  const probe = s2.trace.probes[0], k = probe.indexOf(Math.max(...probe));
  assert.ok(Math.abs(s2.trace.t[k] / tTheory - 1) < 1e-2, `peak at ${s2.trace.t[k]}, theory ${tTheory}`);
  // And it reads the solution: its last sample is the concentration at its node now.
  const node = s2.x.findIndex((x) => Math.abs(x - (xe + d)) < 1e-12);
  assert.ok(Math.abs(probe.at(-1) / s2.c['h+'][node] - 1) < 1e-12);
});
