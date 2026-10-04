import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { check } from '../src/kit.js';

// Devices that aren't planar: a cross-section A(x) that the conservation laws carry, flux × A
// through each face of a box and ∫A dx in it. Spherical and cylindrical shells about a centre
// r0 to the left of x = 0, or any A(x). Each segment's transport is weighted by its length over
// ∫dx/A, so steady diffusion between nodes is exact; currents and fluxes are totals through A.

const T = 298.15, RT = GAS_CONSTANT * T;
const Dx = 1e-9;
// A neutral species X between an electrode (left, at r0) that holds it at cIn and the bulk
// (right) at cOut.
const shell = (geometry, L, { cIn, cOut, grid = { hmin: L / 200, hmax: L / 20, ratio: 1.3 }, c0 } = {}) =>
  new Device({
    T,
    geometry,
    species: [{ name: 'X', z: 0, cRef: 1000 }],
    materials: { m: { epsr: 10, species: { X: { D: Dx, mu0: 0 } } } },
    regions: [{ name: 'm', material: 'm', length: L, ...(c0 ? { c0 } : {}) }],
    contacts: {
      left: cIn === undefined ? { phi: 'neutral' } : { species: { X: { type: 'equilibrium', mu: RT * Math.log(cIn / 1000) } }, phi: 'neutral' },
      right: { species: { X: { type: 'equilibrium', mu: RT * Math.log(cOut / 1000) } }, phi: 'neutral' },
    },
    grid,
  });

test('steady diffusion to a sphere, 4πDΔc/(1/r₀ − 1/r₁), and to a cylinder, 2πDΔc/ln(r₁/r₀): exact on any grid', () => {
  const r0 = 5e-6, L = 200e-6, r1 = r0 + L, cOut = 2, cIn = 1e-12;
  const sphere = shell({ type: 'spherical', r0 }, L, { cIn, cOut }), sph = sphere.solve();
  const report = check(sphere, sph);
  assert.ok(report.ok, report.text);
  const cyl = shell({ type: 'cylindrical', r0 }, L, { cIn, cOut }).solve();
  const wantS = (4 * Math.PI * Dx * (cOut - cIn)) / (1 / r0 - 1 / r1);
  const wantC = (2 * Math.PI * Dx * (cOut - cIn)) / Math.log(r1 / r0);
  // (fluxes toward +x: the electrode draws X in, toward −x)
  assert.ok(Math.abs(-sph.contacts.left.flux.X / wantS - 1) < 1e-10, `${sph.contacts.left.flux.X} vs ${wantS} mol/s`);
  assert.ok(Math.abs(-cyl.contacts.left.flux.X / wantC - 1) < 1e-10, `${cyl.contacts.left.flux.X} vs ${wantC} mol/(m·s)`);
  // The profile is the 1/r one, node by node.
  const g = sph.x.length >> 1, r = r0 + sph.x[g];
  const cS = cIn + ((cOut - cIn) * (1 / r0 - 1 / r)) / (1 / r0 - 1 / r1);
  assert.ok(Math.abs(sph.c.X[g] / cS - 1) < 1e-10);
  // Planar is the default, and A = 1 m² reads as before: D Δc/L.
  const flat = shell(undefined, L, { cIn, cOut }).solve();
  assert.ok(Math.abs(-flat.contacts.left.flux.X / ((Dx * (cOut - cIn)) / L) - 1) < 1e-12);
});

test('diffusion into a spherical electrode after a step: Cottrell plus the spherical term, D c (1/√(πDt) + 1/r₀) per area', () => {
  // Semi-infinite: the bulk 400 µm out, far beyond √(Dt) ≈ 3 µm at 10 ms. Halving the cells
  // cuts the error about fourfold.
  const r0 = 10e-6, L = 400e-6, c = 1, A0 = 4 * Math.PI * r0 * r0;
  const err = [1, 2].map((k) => {
    const dev = shell({ type: 'spherical', r0 }, L, { cIn: 1e-12, cOut: c, c0: { X: c }, grid: { hmin: 0.02e-6 / k, hmax: 20e-6 / k, ratio: 1 + 0.1 / k } });
    return [1e-3, 1e-2].map((t) => {
      const s = dev.advance(t, { tol: 1e-6 });
      assert.ok(s.converged);
      const want = Dx * c * (1 / Math.sqrt(Math.PI * Dx * t) + 1 / r0) * A0;
      return Math.abs(-s.contacts.left.flux.X / want - 1);
    });
  });
  assert.ok(err[1].every((e, j) => e < 3e-4 && err[0][j] / e > 3), `errors ${err}`);
});

test('a sphere filling from its surface, centre at x = 0: the uptake 1 − (6/π²) Σ e^(−n²π²Dt/R²)/n²', () => {
  const R = 10e-6, cs = 1;
  const dev = shell({ type: 'spherical', r0: 0 }, R, { cOut: cs, c0: { X: 1e-12 }, grid: { hmin: R / 400, hmax: R / 100, ratio: 1.1 } });
  const full = (4 / 3) * Math.PI * R ** 3 * cs;
  for (const t of [0.005, 0.02, 0.1]) {
    const s = dev.advance(t, { tol: 1e-6 });
    assert.ok(s.converged && check(dev, s).ok);
    let held = 0;
    const x = s.x;
    for (let g = 0; g < x.length; g++) {
      const lo = g > 0 ? (x[g - 1] + x[g]) / 2 : 0, hi = g < x.length - 1 ? (x[g] + x[g + 1]) / 2 : x[g];
      held += ((4 / 3) * Math.PI * (hi ** 3 - lo ** 3)) * s.c.X[g];
    }
    let sum = 0;
    for (let n = 1; n < 200; n++) sum += Math.exp((-n * n * Math.PI * Math.PI * Dx * t) / (R * R)) / (n * n);
    const want = 1 - (6 / (Math.PI * Math.PI)) * sum;
    assert.ok(Math.abs(held / full / want - 1) < 2e-4, `t = ${t}: ${held / full} vs ${want}`);
  }
});

test('Debye–Hückel around a charged sphere: ψ = ψ₀ (a/r) e^(−κ(r−a)), its charge 4πε a (1 + κa) ψ₀', () => {
  // A sphere of radius a held at ψ₀ (small, so linear) in 1 mM salt, the bulk far away.
  const a = 2e-9, c = 1, epsr = 78.5, eps = epsr * EPS0, psi0 = 1e-4;
  const kappa = Math.sqrt((2 * FARADAY * FARADAY * c) / (eps * RT)), L = 20 / kappa;
  const sol = new Device({
    T,
    geometry: { type: 'spherical', r0: a },
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { water: { epsr, species: { 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } } },
    regions: [{ material: 'water', length: L }],
    contacts: {
      left: { V: psi0, phi: 'pinned', zeroCharge: 0 },
      right: { V: 0, terminal: 'K+', species: { 'K+': { type: 'equilibrium', offset: RT / FARADAY * Math.log(c / 1000) }, 'Cl-': { type: 'equilibrium', offset: -RT / FARADAY * Math.log(c / 1000) } }, phi: 'bulk' },
    },
    grid: { hmin: 0.005e-9, hmax: 0.5e-9, ratio: 1.05 },
  }).solve();
  assert.ok(sol.converged);
  const psi = (g) => sol.phi[g] - sol.phi.at(-1);
  for (const r of [a + 0.5 / kappa, a + 2 / kappa]) {
    const g = sol.x.findIndex((x) => a + x >= r), rr = a + sol.x[g];
    const want = (psi0 * a * Math.exp(-kappa * (rr - a))) / rr;
    assert.ok(Math.abs(psi(g) / want - 1) < 1e-3, `r − a = ${(rr - a) * kappa} Debye lengths: ${psi(g)} vs ${want} V`);
  }
  const Q = 4 * Math.PI * eps * a * (1 + kappa * a) * psi0;
  assert.ok(Math.abs(sol.gates.left.charge / Q - 1) < 3e-4, `${sol.gates.left.charge} vs ${Q} C`);
});

test('geometry is checked', () => {
  const def = (geometry) => ({
    geometry,
    species: [{ name: 'X', z: 0, cRef: 1000 }],
    materials: { m: { epsr: 1, species: { X: { D: 1e-9, mu0: 0 } } } },
    regions: [{ material: 'm', length: 1e-6, velocity: 1e-3, c0: { X: 1 } }],
  });
  assert.throws(() => new Device(def('round')), (e) => e instanceof DeviceError && /geometry must be/.test(e.message));
  assert.throws(() => new Device(def({ type: 'spherical', r0: -1 })), /geometry\.r0/);
  assert.throws(() => new Device(def('spherical')), /flow is for a planar device/);
  const still = def({ area: { x: [0, 0.5e-6, 1e-6], values: [1, 0, 1] } });
  delete still.regions[0].velocity;
  assert.throws(() => new Device(still), /vanishes inside the device/);
});
