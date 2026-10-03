import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, AVOGADRO, units } from '../src/index.js';
import { build, layer, ohmic } from '../src/kit.js';

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

test('set() that closes a stretch conserves what it holds then, not what it held at the start', () => {
  // O₂ diffuses into a film for a while; then its contact is blocked. The film keeps that O₂.
  const RT = 8.314462618 * 298.15, L = 10e-6;
  const open = (c) => ({ species: { O2: { type: 'equilibrium', mu: RT * Math.log(c / 1000) } }, phi: 'neutral' });
  const def = (left) => ({
    species: [{ name: 'O2', z: 0, cRef: 1000 }],
    materials: { polymer: { epsr: 3, species: { O2: { D: 1e-11, mu0: 0 } } } },
    regions: [{ material: 'polymer', length: L, c0: { O2: 1 } }],
    contacts: { left, right: { phi: 'neutral' } },
    grid: { minCells: 50 },
  });
  const dev = new Device(def(open(1)));
  dev.solve();
  dev.set({ contacts: { left: open(2) } });
  const held = dev.advance(1).conservation[0].amount; // partly filled toward c = 2
  assert.ok(held > 1.3e-5 && held < 1.4e-5);
  dev.set({ contacts: { left: { species: { O2: 'blocked' } } } });
  const closed = dev.solve();
  assert.ok(closed.converged);
  assert.ok(Math.abs(closed.conservation[0].amount / held - 1) < 1e-12, `${closed.conservation[0].amount} vs ${held}`);
  assert.ok(Math.abs(closed.conservation[0].drift) < 1e-12);
});

// Slow ions sharing a layer with charged, reacting, immobile traps (SRH through explicit trap
// states, X⁰ + e⁻ = X⁻, X⁻ + h⁺ = X⁰): a cold start needs steps far below anything the ions'
// time scale suggests. (From a perovskite cell: iodide vacancies, D = 1e-17 m²/s, in MAPbI₃.)
test('a cold start with slow ions and immobile traps: solve() and advance() both reach equilibrium', () => {
  const F = FARADAY, T = 298, RT = GAS_CONSTANT * T, m3 = (n) => n / AVOGADRO;
  const gc = 8.1e24, gv = 5.8e24, Eg = 1.7, N0 = 1.6e25, Nt = 1e20;
  const ni = Math.sqrt(gc * gv) * Math.exp(-(Eg * F) / (2 * RT));
  const def = build({
    T,
    species: [{ name: 'e-', z: -1 }, { name: 'h+', z: 1 }, { name: 'X0', z: 0 }, { name: 'X-', z: -1 }, { name: 'V+', z: 1 }],
    materials: {
      MAPI: {
        epsr: 24.1,
        species: {
          'e-': { D: 1.7e-4, mu0: 0, cRef: m3(gc) },
          'h+': { D: 1.7e-4, mu0: units.eV(Eg), cRef: m3(gv) },
          X0: { D: 0, mu0: 0, cRef: m3(Nt) },
          'X-': { D: 0, mu0: RT * Math.log(ni / gc), cRef: m3(Nt) },
          'V+': { D: 1e-17, mu0: 0, cRef: m3(N0) },
        },
      },
    },
    stack: [ohmic(0), layer('MAPI', 400e-9, { fixedCharge: -F * m3(N0), c0: { X0: m3(Nt) / 2, 'X-': m3(Nt) / 2, 'V+': m3(N0) } }), ohmic(0)],
    bulkReactions: [
      { equation: 'e- + X0 = X-', kf: { MAPI: 1 / (3e-9 * m3(Nt)) } },
      { equation: 'X- + h+ = X0', kf: { MAPI: 1 / (3e-7 * m3(Nt)) } },
    ],
    grid: { hmin: 0.05e-9, hmax: 10e-9, ratio: 1.15 },
  });
  const flat = (s) => {
    for (const name of ['e-', 'h+', 'X-']) {
      const mu = s.mu[name];
      if ((Math.max(...mu) - Math.min(...mu)) / RT > 1e-6) return `${name} not flat`;
    }
    return '';
  };
  const s = new Device(def).solve();
  assert.ok(s.converged, 'cold steady solve');
  assert.equal(flat(s), '');
  const dev = new Device(def);
  const r = dev.advance(1e5, { tol: 1e-4 }); // from the cold start, with the default first step
  assert.ok(r.done && r.converged, `advance: done ${r.done}, ${r.steps} steps`);
  assert.equal(flat(dev.solution()), '');
});
