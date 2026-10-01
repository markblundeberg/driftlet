import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, units } from '../src/index.js';

// Slots that aren't unknowns (absent species, undefined φ, blocked or neutral interface fluxes,
// a metal's idle slots) are left out of the linear system. Leaving them out must change nothing.

const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
const silicon = { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } };
const carriers = [
  { name: 'e-', z: -1 },
  { name: 'h+', z: 1 },
];
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });

// n-Si | KCl across a Helmholtz layer: four species, three unknowns at each node.
const semiconductorElectrolyte = () =>
  new Device({
    species: [...carriers, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { Si: silicon, water: { epsr: 78.5, species: { 'K+': { D: 1.96e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
    regions: [
      { material: 'Si', length: 0.5e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'water', length: 0.5e-6, c0: { 'K+': 100, 'Cl-': 100 } },
    ],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, dipole: 0 }],
    contacts: { left: ohmic(0), right: { V: 0.3, bath: { c: { 'K+': 100, 'Cl-': 100 }, reference: 'Cl-' } } },
    grid: { hmin: 0.1e-9, hmax: 20e-9, ratio: 1.2 },
  });

// A floating silver plate between plating electrodes: metal nodes have two unknowns, the
// electrolyte none for electrons.
const ions = [
  { name: 'Ag+', z: 1, cRef: 1000 },
  { name: 'NO3-', z: -1, cRef: 1000 },
];
const plating = { reactants: { 'Ag+': 1 }, electrons: 1, products: { Ag: 1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 };
const bipolar = () => {
  const water = { epsr: 78.5, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } };
  const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating] };
  const end = { V: 0, phi: { type: 'capacitive', C: 0.2, zeroCharge: 0.1 }, reactions: [plating] };
  const salt = { 'NO3-': 10, 'Ag+': 10 };
  return new Device({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water, Ag: { conductor: { species: 'e-', conductivity: 6e7 } } },
    regions: [
      { material: 'water', length: 5e-6, c0: salt },
      { material: 'Ag', length: 1e-6 },
      { material: 'water', length: 5e-6, c0: salt },
    ],
    interfaces: [stern, stern],
    contacts: { left: end, right: { ...end, V: 0.4 } },
    grid: { hmin: 0.1e-9, hmax: 100e-9, ratio: 1.2 },
  });
};

test('left-out slots have identity rows and zero residuals, and the reduced solve equals the full one', () => {
  for (const make of [semiconductorElectrolyte, bipolar]) {
    const dev = make();
    const sol = dev.solve();
    assert.ok(sol.converged);
    const s = dev.solver;
    assert.ok(s.compact, `${make.name}: expected a reduced system`);
    const { M, nB, sys, active } = s;
    // Away from the solution, so the update is non-trivial.
    s.u[s.blockOfNode[5] * M + 1] += 0.1;
    s.assemble(1e-6);
    s._equilibrate();
    for (let b = 0; b < nB; b++) {
      for (let r = 0; r < M; r++) {
        if (active[b * M + r]) continue;
        const o = b * M * M + r * M;
        for (let c = 0; c < M; c++) {
          assert.equal(sys.B[o + c], c === r ? 1 : 0, `${make.name}: block ${b} slot ${r}`);
          assert.equal(sys.A[o + c], 0);
          assert.equal(sys.C[o + c], 0);
        }
        assert.equal(s.res[b * M + r], 0);
      }
    }
    const reduced = new Float64Array(nB * M), full = new Float64Array(nB * M);
    s._factor();
    s._solveLinear(s.res, reduced);
    sys.factor();
    sys.solve(s.res, full);
    const scale = full.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
    const diff = full.reduce((m, v, k) => Math.max(m, Math.abs(v - reduced[k])), 0);
    assert.ok(scale > 1e-3 && diff < 1e-12 * scale, `${make.name}: ${diff} of ${scale}`);
    const sizes = s.lin.sizes;
    assert.ok(Math.max(...sizes) < M && sizes.reduce((a, m) => a + m ** 3, 0) < 0.5 * nB * M ** 3);
  }
});
