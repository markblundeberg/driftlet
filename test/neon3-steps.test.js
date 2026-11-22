// Neon's third review (21 Nov 2026): the time steps' charge floor, and an OCV table's check.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';

const VT = (GAS_CONSTANT * 298.15) / FARADAY;
const logit = (x) => Math.log(x / (1 - x));

// A 1 µm spherical electrode in a 1 cm cell of 100 mM KCl, K⁺ through a conductance G = 10 S/m²:
// 10 ms at rest, then 0 → 10 mV over 100 µs. The charge floor was the grid's largest area (the
// outer sphere's, 1e8 times the electrode's) for every terminal, which kept the electrode's charge
// out of the error: the ramp went in one step, Q 50% high, no warning. What it passes is G A ∫V dt
// to the cell's polarization, 0.02% here (a run at dtMax 1e-6 gives 6.2833e-17 C).
test("a small spherical electrode's charge is held to tol (each terminal's floor by its own area)", () => {
  const off = VT * Math.log(100 / 1000);
  const dev = new Device({
    geometry: { type: 'spherical', r0: 1e-6 },
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr: 0, species: { 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } } },
    regions: [{ name: 'el', material: 'aq', length: 1e-2, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: {
      left: { V: { t: [0, 1e-2, 1.01e-2], values: [0, 0, 0.01] }, phi: 'neutral', species: { 'K+': { type: 'conductance', G: 10, offset: off }, 'Cl-': 'blocked' } },
      right: { V: 0, bath: { c: { 'K+': 100, 'Cl-': 100 } } },
    },
    grid: { hmin: 1e-7, hmax: 1e-3, ratio: 1.5 },
  });
  dev.solve();
  const run = dev.advance(1.01e-2, { tol: 1e-3 });
  assert.ok(run.converged && run.done, run.warnings.join('; '));
  const A = 4 * Math.PI * 1e-12, ideal = A * 10 * 0.005 * 1e-4;
  const q = run.trace.terminals.left.charge.at(-1);
  assert.ok(Math.abs(q / ideal - 1) < 1e-2, `Q = ${q.toExponential(5)} C, G A ∫V dt = ${ideal.toExponential(5)} C`);
  assert.ok(run.trace.t.filter((t) => t > 1e-2).length > 3, 'the ramp in more than one step');
});

// An insertion host's OCV table, E strictly falling, its first knot at x = 1e-4: the end slope
// PCHIP gives there turns ζ(x) back over (0.9e-4, 1.02e-4]. The check sampled a grid 2.5e-4 apart
// and each interval's interior, never a knot nor the ends, and passed it: x(ζ) at the knot's own
// ζ read 8.09e-5 (19% low), a solve held at E[0] 6.7% off the knot. Refused now, at the knot; and
// mirrored (x → 1 − x, r → −r), at the last one.
test('an OCV table whose curve turns back at an end knot, or beyond it, is refused there', () => {
  const host = (xs, rs) => {
    const E = xs.map((x, i) => 3.5 - VT * (logit(x) + rs[i]));
    assert.ok(E.every((e, i) => i === 0 || e < E[i - 1]), 'E strictly falls');
    return () => new Device({
      species: [{ name: 'Li+', z: 1, cRef: 1000 }, { name: 'e-', z: -1, cRef: 1000 }],
      materials: { host: { epsr: 0, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 50000, ocv: { x: xs, E, muRef: 3.5 * FARADAY } }], species: { 'Li+': { D: 1e-12, mu0: 0 }, 'e-': { D: 1e-6, mu0: 0 } } } },
      regions: [{ name: 'cathode', material: 'host', length: 1e-5, c0: { 'Li+': 25000, 'e-': 25000 } }],
      contacts: { left: { V: 0, phi: 'neutral', species: { 'Li+': 'blocked', 'e-': { type: 'equilibrium', offset: 0 } } }, right: { V: -3.6, phi: 'neutral', species: { 'Li+': { type: 'equilibrium', offset: 3.5 }, 'e-': 'blocked' } } },
    });
  };
  const xs = [1e-4, 1.5e-4, 1.6e-4, 0.5, 0.9], rs = [0, -3700 * 5e-5, -3700 * 5e-5 + 20000 * 1e-5, 0, 0];
  assert.throws(host(xs, rs), /ocv: near x = 0\.000100000 the curve is flat or rising/);
  const mx = xs.map((x) => 1 - x).reverse(), mr = rs.map((r) => -r).reverse();
  assert.throws(host(mx, mr), /ocv: near x = 0\.999900 the curve is flat or rising/);
  // the same knots with a tenth of that dip turn nothing back, and pass
  assert.doesNotThrow(host(xs, rs.map((r) => r / 10)));
  assert.doesNotThrow(host(mx, mr.map((r) => r / 10)));
});
