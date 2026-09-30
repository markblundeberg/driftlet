import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;

// For a linear single-pole response with time constant τ, backward Euler with step dt gives
// exactly x_n = x_∞ (1 − (1 + dt/τ)^(−n)), so time-stepping error drops out of the comparison.
const beStep = (n, dt, tau) => 1 - (1 + dt / tau) ** -n;

test('floating island charging after a gate step: RC time constant, Gauss, conservation', () => {
  const Dp = 1.33e-9, Dm = 2.03e-9, c = 10, L = 3e-6, Cg = 0.2, epsr = 78.5;
  const gate = (V) => ({ V, phi: { type: 'capacitive', C: Cg, zeroCharge: 0 } });
  const def = {
    species: [
      { name: 'Na+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr, species: { 'Na+': { D: Dp, mu0: -261.9e3 }, 'Cl-': { D: Dm, mu0: -131.2e3 } } } },
    regions: [{ material: 'water', length: L, c0: { 'Na+': c, 'Cl-': c } }],
    contacts: { left: gate(0), right: gate(0) },
    grid: { hmin: 0.05e-9, hmax: 20e-9, ratio: 1.1 },
  };
  const dev = new Device(def);
  dev.solve();
  // Small step on the left gate. Series circuit: gate—C_s—R_bulk—C_s—gate, with
  // C_s = C_g ∥series∥ ε/λ_D and R_bulk = L/σ, so τ = R_bulk·C_s/2.
  const dV = 1e-3;
  dev.set({ contacts: { left: { V: dV } } });
  const qInf = new Device(dev.def).solve().gates.left.charge;
  const lam = Math.sqrt((epsr * EPS0 * RT) / (2 * FARADAY * FARADAY * c));
  const Cs = 1 / (1 / Cg + lam / (epsr * EPS0));
  const sigma = (FARADAY * FARADAY * c * (Dp + Dm)) / RT;
  const tau = ((L / sigma) * Cs) / 2;
  assert.ok(Math.abs(qInf / ((Cs * dV) / 2) - 1) < 1e-3);

  const dt = tau / 20;
  for (let n = 1; n <= 80; n++) {
    const sol = dev.step(dt);
    assert.ok(sol.converged);
    const q = sol.gates.left.charge / qInf;
    // Distributed double-layer charging adds a small correction (~λ_D/L and ~τ_D/τ).
    assert.ok(Math.abs(q - beStep(n, dt, tau)) < 1e-2, `n=${n}: ${q} vs ${beStep(n, dt, tau)}`);
    const gross = FARADAY * 2 * c * L;
    assert.ok(Math.abs(sol.charge + sol.gates.left.charge + sol.gates.right.charge) < 1e-12 * gross);
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `${st.species} drift ${st.drift}`);
    // Gate current is displacement current, equal at both ends. (Newton stops at δφ̂ ~ 1e-10,
    // which bounds the displacement difference quotient to ~1e-7 relative.)
    assert.ok(Math.abs(sol.contacts.left.current - sol.contacts.right.current) < 1e-6 * Math.abs(sol.contacts.left.current));
  }
});

test('water: homogeneous relaxation of excess H⁺/OH⁻ at rate k (c_H + c_OH)', () => {
  const kf = 1.4e8; // m³/(mol·s), i.e. 1.4e11 M⁻¹s⁻¹
  const muH2O = -237.13e3, mu0OH = -157.24e3;
  const Kw = 1e6 * Math.exp((muH2O - mu0OH) / RT);
  const ceq = Math.sqrt(Kw);
  const excess = 1e-3 * ceq; // small, so the response is linear
  const def = {
    species: [
      { name: 'H+', z: 1, cRef: 1000 },
      { name: 'OH-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 78.5, species: { 'H+': { D: 9.31e-9, mu0: 0 }, 'OH-': { D: 5.27e-9, mu0: mu0OH } } } },
    regions: [{ material: 'water', length: 50e-9, c0: { 'H+': ceq + excess, 'OH-': ceq + excess } }],
    bulkReactions: [{ reactants: { 'H+': 1, 'OH-': 1 }, products: { H2O: 1 }, fixed: { H2O: muH2O }, kf: { water: kf } }],
    contacts: {
      left: { phi: { type: 'capacitive', C: 0.2, zeroCharge: 0 } },
      right: { phi: { type: 'capacitive', C: 0.2, zeroCharge: 0 } },
    },
    grid: { minCells: 10 },
  };
  const dev = new Device(def);
  const tau = 1 / (kf * 2 * ceq);
  const dt = tau / 10;
  const m = dev.grid.nNodes >> 1;
  for (let n = 1; n <= 30; n++) {
    const sol = dev.step(dt);
    assert.ok(sol.converged);
    const left = (sol.c['H+'][m] - ceq) / excess; // fraction of the excess remaining
    const want = 1 - beStep(n, dt, tau);
    assert.ok(Math.abs(left - want) < 2e-3, `n=${n}: ${left} vs ${want}`);
  }
});

test('open system: amounts track the time-integrated contact fluxes after a bias step', () => {
  const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
  const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
  const dev = new Device({
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } } },
    regions: [
      { material: 'Si', length: 2e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'Si', length: 2e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
    ],
    contacts: { left: ohmic(0), right: ohmic(0) },
    grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
  });
  dev.solve();
  dev.set({ contacts: { right: { V: 0.4 } } });
  let dt = 1e-13;
  for (let n = 0; n < 40; n++, dt *= 1.5) {
    const sol = dev.step(dt);
    assert.ok(sol.converged, `step ${n}`);
    for (const st of sol.conservation) {
      assert.ok(st.connected);
      assert.ok(Math.abs(st.drift) < 1e-11, `step ${n} ${st.species}: drift ${st.drift}`);
    }
  }
});
