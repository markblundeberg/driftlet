import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';

const VT = (GAS_CONSTANT * 298.15) / FARADAY;

// Ag | AgNO₃ | Ag, strictly neutral, with the nitrate blocked at both electrodes. In steady state
// the nitrate carries no flux, so c ∝ e^{φ̂} and the silver flux N = −2D∇c makes c linear about
// its mean c₀, which only the conserved nitrate amount fixes. The electrode voltage is then
// V = 2V_T ln(c_L/c_R), so I = (4 F D c₀ / L) tanh(V / 4V_T).
test('steady state with a conserved spectator, solved directly: I = (4FDc₀/L) tanh(V/4V_T)', () => {
  const D = 1.65e-9, c0 = 10, L = 20e-6;
  const electrode = (V) => ({ V, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
  const cell = (cells) =>
    new Device({
      species: [
        { name: 'Ag+', z: 1, cRef: 1000 },
        { name: 'NO3-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'Ag+': { D, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
      regions: [{ material: 'water', length: L, c0: { 'NO3-': c0 } }],
      contacts: { left: electrode(0), right: electrode(0) },
      grid: { minCells: cells },
    });
  const error = (dev, V) => {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged);
    assert.equal(sol.steps, 1, 'one steady (dt = ∞) solve, with the nitrate amount as a constraint');
    for (const st of sol.conservation) if (st.spectator) assert.ok(Math.abs(st.drift) < 1e-13, `drift ${st.drift}`);
    return Math.abs(sol.current) / (((4 * FARADAY * D * c0) / L) * Math.tanh(V / (4 * VT))) - 1;
  };
  const coarse = cell(400), fine = cell(800);
  for (const V of [0.01, 0.1]) assert.ok(Math.abs(error(coarse, V)) < 1e-5, `V=${V}`);
  // Near depletion (c_R ≈ 0.006 c₀ at 0.3 V) the error is larger, and second order.
  const e1 = error(coarse, 0.3), e2 = error(fine, 0.3);
  assert.ok(Math.abs(e1) < 5e-4 && e1 / e2 > 3.6 && e1 / e2 < 4.4, `${e1}, ${e2}`);
});
