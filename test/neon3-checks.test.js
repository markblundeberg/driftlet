import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY } from '../src/index.js';
import { check } from '../src/checks.js';

// The grid check's comparisons, after Neon's third review (flags 1.2 and 1.3).
const grid = (report) => report.items.find((it) => it.name === 'grid');
const diff = (report, re) => grid(report).details.diffs.find((d) => re.test(d.what));

// Flag 1.2: a port's reaction, summed over its window, by the port's own nodes. A 1 µm catalyst
// layer, finely gridded, beside 1 mm of electrolyte on 5 cells: the face's two nodes share an x,
// and keyed by x the port's last node read the bulk's 100 µm half-box for its own 25 nm one, a
// total 100× too large that halved on refinement (FAIL grid, 0.5) where everything moved by 4e-13.
test('grid check: a port reaching a face sums its reaction over its own boxes, not its neighbour’s', () => {
  const c0 = { O: 100, R: 100, 'K+': 100 };
  const dev = new Device({
    species: [
      { name: 'e-', z: -1, cRef: 1000 },
      { name: 'O', z: 0, cRef: 1000 },
      { name: 'R', z: -1, cRef: 1000 },
      { name: 'K+', z: 1, cRef: 1000 },
    ],
    materials: { aq: { epsr: 0, species: { O: { D: 1e-9, mu0: 0 }, R: { D: 1e-9, mu0: 0 }, 'K+': { D: 1e-9, mu0: 0 } } } },
    regions: [
      { name: 'cat', material: 'aq', length: 1e-6, c0, grid: { hmax: 5e-8, minCells: 20 } },
      { name: 'bulk', material: 'aq', length: 1e-3, c0, grid: { hmax: 2e-4, minCells: 5 } },
    ],
    interfaces: [{ phi: 'neutral', species: { O: 'equilibrium', R: 'equilibrium', 'K+': 'equilibrium' } }],
    contacts: {
      left: { V: 0, phi: 'neutral', species: { O: 'blocked', R: 'blocked', 'K+': 'blocked' } },
      right: { V: 0, bath: { c: c0 } },
    },
    ports: [{ name: 'we', terminal: 'e-', region: 'cat', from: 0, to: 1e-6, area: 1e6, V: -0.02, reactions: [{ equation: 'O + e- = R', k0: 1e-9, alpha: 0.5 }] }],
  });
  const sol = dev.solve(), report = check(dev, sol, { order: false });
  assert.ok(sol.converged);
  // What the reaction turns over is what the port makes of R (O + e- = R, one each), to round-off.
  const rx = diff(report, /port we's reaction/), made = sol.ports[0].flux.R;
  assert.ok(Math.abs(rx.a - made) <= 1e-9 * Math.abs(made), `turnover ${rx.a} against the port's R ${made}`);
  assert.ok(rx.rel < 1e-6, `the turnover moved by ${rx.rel} on refinement`);
  assert.equal(grid(report).ok, true, grid(report).summary);
  assert.ok(Math.abs(made * FARADAY + sol.terminals.we.current) <= 1e-9 * Math.abs(sol.terminals.we.current)); // (the port's current is its reaction's)
});

// Flag 1.3: in a sphere, a face's D (and its rates) are per area, the contacts' D and the regions'
// charges totals. Per area beside them, a 10 µm particle's face D (9.6e-5 C/m², its 1.5e-13 C)
// set the floor, 6.6e8× too high, and the outer contact's double layer, 4.9% off on hmin =
// 2.5 nm < λ_D = 3 nm, passed: "ok grid … face 0's D, 9.6301e-5 → 9.5981e-5 C".
test('grid check: in a sphere a face’s D and rates are compared as totals, and an unresolved double layer fails', () => {
  const VT = (8.314462618 * 298.15) / FARADAY, off = VT * Math.log(10 / 1000);
  const def = {
    geometry: { type: 'spherical', r0: 1e-5 },
    species: [
      { name: 'e-', z: -1, cRef: 1000 },
      { name: 'K+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: {
      Pt: { conductor: { species: 'e-', conductivity: 1e7 } },
      aq: { epsr: 78.5, species: { 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } },
    },
    regions: [
      { name: 'core', material: 'Pt', length: 1e-6 },
      { name: 'shell', material: 'aq', length: 1e-6, c0: { 'K+': 10, 'Cl-': 10 } },
    ],
    interfaces: [{ phi: { type: 'capacitive', C: 0.01 }, zeroCharge: 0 }],
    contacts: {
      left: { V: 0.01, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } },
      right: { V: 0, phi: 'pinned', zeroCharge: 0.01, species: { 'K+': { type: 'conductance', G: 1e-2, offset: off }, 'Cl-': { type: 'conductance', G: 1e-2, offset: -off } } },
    },
    grid: { hmin: 2.5e-9, hmax: 2e-7, ratio: 2.5 },
  };
  const dev = new Device(def), sol = dev.solve(), report = check(dev, sol, { order: false });
  assert.ok(sol.converged);
  const A = dev.model.grid.area[dev.model.grid.regionEnd[0]];
  const D = diff(report, /face 0's D/);
  assert.equal(D.unit, 'C');
  assert.ok(Math.abs(D.a - A * sol.interfaces[0].D) <= 1e-12 * Math.abs(D.a), `face D ${D.a} C against A·D ${A * sol.interfaces[0].D}`);
  // The outer contact's D, on a grid twice as fine by hand, moves by ~5%: the check must see it.
  const fine = new Device({ ...def, grid: { hmin: 1.25e-9, hmax: 1e-7, ratio: 2.5 } }).solve();
  const truth = Math.abs(fine.contacts.right.D - sol.contacts.right.D) / Math.abs(fine.contacts.right.D);
  assert.ok(truth > 0.03, `the double layer's own change, ${truth}`);
  assert.equal(grid(report).ok, false, grid(report).summary);
  assert.ok(diff(report, /right contact's D/).rel > 0.03);

  // A face's reaction in a sphere: its total, A × rate, in mol/s.
  const c0 = { O: 100, R: 100, 'K+': 100 };
  const redox = new Device({
    geometry: { type: 'spherical', r0: 1e-5 },
    species: [
      { name: 'e-', z: -1, cRef: 1000 },
      { name: 'O', z: 0, cRef: 1000 },
      { name: 'R', z: -1, cRef: 1000 },
      { name: 'K+', z: 1, cRef: 1000 },
    ],
    materials: {
      Pt: { conductor: { species: 'e-', conductivity: 1e7 } },
      aq: { epsr: 0, species: { O: { D: 1e-9, mu0: 0 }, R: { D: 1e-9, mu0: 0 }, 'K+': { D: 1e-9, mu0: 0 } } },
    },
    regions: [
      { name: 'core', material: 'Pt', length: 1e-6 },
      { name: 'shell', material: 'aq', length: 1e-5, c0 },
    ],
    interfaces: [{ reactions: [{ equation: 'O + e- = R', k0: 1e-6, alpha: 0.5 }] }],
    contacts: {
      left: { V: -0.02, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } },
      right: { V: 0, bath: { c: c0 } },
    },
  });
  const rs = redox.solve(), rr = check(redox, rs, { order: false });
  assert.ok(rs.converged);
  const Ar = redox.model.grid.area[redox.model.grid.regionEnd[0]], rate = diff(rr, /face 0's reaction/);
  assert.equal(rate.unit, 'mol/s');
  assert.ok(Math.abs(rate.a - Ar * rs.interfaces[0].rates[0]) <= 1e-12 * Math.abs(rate.a), `face rate ${rate.a} mol/s against A·rate ${Ar * rs.interfaces[0].rates[0]}`);
  // (and it is what the face makes of R, which leaves through the right contact, toward +x)
  assert.ok(Math.abs(rate.a - rs.contacts.right.flux.R) <= 1e-6 * Math.abs(rate.a), `${rate.a} against the contact's R ${rs.contacts.right.flux.R}`);
});
