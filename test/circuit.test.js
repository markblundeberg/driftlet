import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';

const VT = (GAS_CONSTANT * 298.15) / FARADAY;

test('liquid junction at open circuit: cell EMF and Planck diffusion potential', () => {
  // Two NaCl baths, each read through Cl⁻ (as with Ag/AgCl electrodes), joined by a
  // diffusion zone. At zero current, E = 2 t₊ (RT/F) ln(c₁/c₂), and the diffusion potential
  // is (t₊ − t₋)(RT/F) ln(c₁/c₂).
  const Dp = 1.33e-9, Dm = 2.03e-9, c1 = 100, c2 = 10;
  const tp = Dp / (Dp + Dm);
  const bath = (c) => ({ bath: { c: { 'Na+': c, 'Cl-': c }, reference: 'Cl-' } });
  const dev = new Device({
    species: [
      { name: 'Na+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 78.5, species: { 'Na+': { D: Dp, mu0: -261.9e3 }, 'Cl-': { D: Dm, mu0: -131.2e3 } } } },
    regions: [{ material: 'water', length: 10e-6 }],
    contacts: { left: bath(c1), right: { ...bath(c2), I: 0 } },
    grid: { hmin: 0.2e-9, hmax: 50e-9, ratio: 1.15 },
  });
  const sol = dev.solve();
  assert.ok(sol.converged);
  const emf = 2 * tp * VT * Math.log(c1 / c2);
  assert.ok(Math.abs(sol.terminalVoltage / emf - 1) < 1e-4, `EMF ${sol.terminalVoltage} vs ${emf}`);
  const diff = (tp - (1 - tp)) * VT * Math.log(c1 / c2);
  const dphi = sol.phi[sol.phi.length - 1] - sol.phi[0];
  assert.ok(Math.abs(dphi / diff - 1) < 1e-4, `diffusion potential ${dphi} vs ${diff}`);
  // Individual ion fluxes are large; the net current cancels.
  const I = sol.contacts.right.current;
  assert.ok(Math.abs(I) < 1e-9 * FARADAY * Math.abs(sol.contacts.right.flux['Na+']));
});

// Ag | AgNO₃ | Ag: Ag⁺ reversible at both electrodes, NO₃⁻ blocked (a conserved spectator).
const Dp = 1.65e-9, Dm = 1.9e-9, c0 = 10, L = 20e-6;
const links = { terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
const electrode = (V = 0) => ({ V, ...links });
const cell = (right = electrode()) => ({
  species: [
    { name: 'Ag+', z: 1, cRef: 1000 },
    { name: 'NO3-', z: -1, cRef: 1000 },
  ],
  materials: { water: { epsr: 78.5, species: { 'Ag+': { D: Dp, mu0: 77.1e3 }, 'NO3-': { D: Dm, mu0: -111.3e3 } } } },
  regions: [{ material: 'water', length: L, c0: { 'NO3-': c0 } }],
  contacts: { left: electrode(), right },
  grid: { hmin: 0.2e-9, hmax: 100e-9, ratio: 1.15 },
});
// Blocked anion ⇒ flat μ̄₋ and a linear salt profile; the cation's current is −2 F D₊ ∂c/∂x.
// Salt is conserved, so the depletion saturates at i_lim = 4 F D₊ c₀ / L, and
// i = i_lim · tanh(V / 4V_T) (current flows toward −x for V > 0 on the right electrode).
const iLim = (4 * FARADAY * Dp * c0) / L;

test('concentration polarization: I–V and limiting current (quasi-neutral analytic)', () => {
  const dev = new Device(cell());
  for (const V of [-0.1, 0.01, 0.05, 0.1, 0.2]) {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `V=${V}`);
    const want = -iLim * Math.tanh(V / (4 * VT));
    assert.ok(Math.abs(sol.current - want) < 2e-4 * iLim, `V=${V}: ${sol.current} vs ${want}`);
    assert.ok(Math.abs(sol.contacts.left.current - sol.current) < 1e-9 * iLim);
    const no3 = sol.conservation.find((st) => st.species === 'NO3-');
    assert.ok(no3.spectator && Math.abs(no3.drift) < 1e-12);
  }
});

test('a current-driven terminal inverts the I–V; open circuit gives zero voltage', () => {
  // (I into the device at the right terminal: −0.5 i_lim is 0.5 i_lim toward +x)
  const dev = new Device(cell({ ...links, I: -0.5 * iLim }));
  const sol = dev.solve();
  assert.ok(sol.converged);
  assert.ok(Math.abs(sol.terminalVoltage - -4 * VT * Math.atanh(0.5)) < 1e-4 * VT);
  assert.ok(Math.abs(sol.current / (0.5 * iLim) - 1) < 1e-9);
  dev.set({ contacts: { right: { I: 0 } } });
  assert.ok(Math.abs(dev.solve().terminalVoltage) < 1e-12);
});

test('a source behind a resistance: I = (V − V_source)/R toward +x', () => {
  for (const [R, Vs] of [
    [1e-3, 0.1],
    [1e-4, -0.05],
  ]) {
    const sol = new Device(cell({ ...links, V: Vs, R })).solve();
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.current - (sol.terminalVoltage - Vs) / R) < 1e-9 * Math.abs(sol.current));
    // And it sits on the cell's own I–V curve.
    const want = -iLim * Math.tanh(sol.terminalVoltage / (4 * VT));
    assert.ok(Math.abs(sol.current - want) < 2e-4 * iLim);
  }
});

test('terminal drives are checked', () => {
  const throwsDevice = (def, pattern) => assert.throws(() => new Device(def), (e) => e instanceof DeviceError && pattern.test(e.message), String(pattern));
  throwsDevice(cell({ I: 1 }), /contacts\.right\.I: this contact passes no current/);
  throwsDevice(cell({ ...links, V: 0, I: 1 }), /give V or I, not both/);
  throwsDevice(cell({ ...links, I: 1, R: 1 }), /a series resistance goes with a voltage source/);
  const def = cell({ ...links, I: 0 });
  def.contacts.left = { ...links, I: 0 };
  throwsDevice(def, /hold at least one at a voltage/);
});
