import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// A neutral species between two reservoirs at c_L and c_R.
const D = 1e-9, L = 1e-5, cL = 2, cR = 1;
const bar = (region) =>
  new Device({
    species: [{ name: 'X', z: 0, cRef: 1000 }],
    materials: { m: { epsr: 10, species: { X: { D, mu0: 0 } } } },
    regions: [{ material: 'm', length: L, ...region }],
    contacts: {
      left: { species: { X: { type: 'equilibrium', mu: RT * Math.log(cL / 1000) } }, phi: 'neutral' },
      right: { species: { X: { type: 'equilibrium', mu: RT * Math.log(cR / 1000) } }, phi: 'neutral' },
    },
    grid: { hmin: L / 10, hmax: L / 10 },
  });

test('advection: exact convection–diffusion profile and flux on a coarse grid, either direction', () => {
  // N = v c − D c′ is constant: c = c_L + (c_R − c_L)(e^{vx/D} − 1)/(e^{Pe} − 1), Pe = vL/D.
  for (const v of [1e-4, -3e-4, 1e-6]) {
    const sol = bar({ velocity: v }).solve();
    assert.ok(sol.converged);
    const Pe = (v * L) / D;
    const want = (v * (cL * Math.exp(Pe) - cR)) / Math.expm1(Pe);
    assert.ok(Math.abs(sol.contacts.left.flux.X / want - 1) < 1e-12, `v=${v}`);
    assert.ok(Math.abs(sol.contacts.right.flux.X / want - 1) < 1e-12);
    for (let g = 0; g < sol.x.length; g++) {
      const exact = cL + ((cR - cL) * Math.expm1((v * sol.x[g]) / D)) / Math.expm1(Pe);
      assert.ok(Math.abs(sol.c.X[g] / exact - 1) < 1e-12);
    }
  }
});

test('eddy mixing: a neutral species diffuses with D + D_mix, exactly', () => {
  for (const mixing of [1e-9, 1e-7]) {
    const sol = bar({ mixing }).solve();
    assert.ok(sol.converged && sol.iterations <= 8);
    assert.ok(Math.abs(sol.contacts.left.flux.X / (((D + mixing) * (cL - cR)) / L) - 1) < 1e-12);
  }
});

test('eddy mixing carries no current: junction EMF unchanged, salt flux raised by D_mix', () => {
  const Dp = 1.33e-9, Dm = 2.03e-9, tp = Dp / (Dp + Dm), c1 = 100, c2 = 10, Lj = 10e-6;
  const bath = (c) => ({ bath: { c: { 'Na+': c, 'Cl-': c }, reference: 'Cl-' } });
  const junction = (mixing) =>
    new Device({
      species: [
        { name: 'Na+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'Na+': { D: Dp, mu0: -261.9e3 }, 'Cl-': { D: Dm, mu0: -131.2e3 } } } },
      regions: [{ material: 'water', length: Lj, mixing }],
      contacts: { left: bath(c1), right: { ...bath(c2), I: 0 } },
      grid: { minCells: 400 },
    }).solve();
  const Ds = (2 * Dp * Dm) / (Dp + Dm), emf = 2 * tp * VT * Math.log(c1 / c2);
  for (const mixing of [0, 1e-8]) {
    const sol = junction(mixing);
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.terminalVoltage / emf - 1) < 1e-5, `EMF with D_mix = ${mixing}`);
    const flux = ((Ds + mixing) * (c1 - c2)) / Lj;
    assert.ok(Math.abs(sol.contacts.left.flux['Na+'] / flux - 1) < 1e-6, `salt flux with D_mix = ${mixing}`);
  }
});

test('eddy mixing leaves equilibrium double layers alone', () => {
  // A gated island with strong mixing: still flat μ̄ and the same gate charge as without it.
  const gate = (V) => ({ V, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 });
  const island = (mixing) =>
    new Device({
      species: [
        { name: 'Na+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 78.5, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } } },
      regions: [{ material: 'water', length: 200e-9, c0: { 'Na+': 10, 'Cl-': 10 }, mixing }],
      contacts: { left: gate(0.3), right: gate(-0.3) },
      grid: { hmin: 0.05e-9, hmax: 5e-9, ratio: 1.1 },
    }).solve();
  const a = island(0), b = island(1e-6);
  assert.ok(b.converged);
  assert.ok(Math.abs(b.gates.left.charge / a.gates.left.charge - 1) < 1e-10);
  for (const mu of Object.values(b.mu)) assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-9);
});

test('flow and mixing are checked', () => {
  const def = bar({}).def;
  def.regions[0].mixing = -1;
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /mixing must be a non-negative number/.test(e.message));
  def.regions[0].mixing = 0;
  def.regions[0].velocity = 'fast';
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /velocity must be a finite number/.test(e.message));
});

test('exchange link at a contact: N = k (μ_out − μ)/RT in series with diffusion', () => {
  // The bar's left end exchanges X with an outside phase at c_out through a rate constant k; the
  // right end holds c_R. Steady: N = D (c_0 − c_R)/L = k ln(c_out/c_0), solved for c_0.
  const k = 1e-4, cOut = 2;
  const dev = bar();
  dev.set({ contacts: { left: { species: { X: { type: 'exchange', k, mu: RT * Math.log(cOut / 1000) } } } } });
  const sol = dev.solve();
  assert.ok(sol.converged);
  let lo = cR, hi = cOut;
  for (let it = 0; it < 200; it++) {
    const c0 = 0.5 * (lo + hi);
    if ((D * (c0 - cR)) / L > k * Math.log(cOut / c0)) hi = c0;
    else lo = c0;
  }
  const want = (D * (lo - cR)) / L;
  assert.ok(Math.abs(sol.contacts.left.flux.X / want - 1) < 1e-10, `${sol.contacts.left.flux.X} vs ${want}`);
  assert.ok(Math.abs(sol.c.X[0] / lo - 1) < 1e-10);
});
