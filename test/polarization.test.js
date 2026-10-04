import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer, ohmic, bath, half, level, aqueous, metal, polarization } from '../src/kit.js';

// Polarization curves at an electrode with the solution beside it held: the curves of an Evans
// diagram. They're the solver's own rate law, so at the metal's actual level they must be the
// rates the solve found, and each couple's level is where its rate vanishes.

const VT = (GAS_CONSTANT * 298.15) / FARADAY;
const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 });

test('at a face: the solved rates at the metal\'s own level, the couple\'s level where the rate vanishes, Tafel slopes (1 − α) and α', () => {
  const alpha = 0.3, face = { reactions: [{ ...silver, k0: 1e-3, alpha }] };
  const dev = new Device(build({
    library: [aqueous(['Ag+', 'NO3-'], { epsr: 0 }), metal('Ag')],
    stack: [ohmic(0, ['e-']), layer('Ag', 1e-6), face, layer('water', 20e-6, { c0: { 'Ag+': 10, 'NO3-': 10 } }), face, layer('Ag', 1e-6), ohmic(0.05, ['e-'])],
    grid: { hmin: 0.2e-6, hmax: 1e-6 },
  }));
  const sol = dev.solve();
  assert.ok(sol.converged);
  const { regionStart, regionEnd } = dev.model.grid;
  for (const [f, metalNode, edge] of [[0, regionEnd[0], regionStart[1]], [1, regionStart[2], regionEnd[1]]]) {
    const Ve = sol.V['e-'][metalNode];
    const p = polarization(dev, sol, { face: f }, [Ve, Ve + 0.3, Ve + 0.4, Ve - 0.3, Ve - 0.4]);
    const r = p.reactions[0];
    assert.ok(Math.abs(r.rate[0] / sol.interfaces[f].rates[0] - 1) < 1e-9, `${r.rate[0]} vs ${sol.interfaces[f].rates[0]}`);
    assert.equal(r.equation, silver.equation);
    // The couple's level at the face is silver's redox level in the solution beside it.
    assert.ok(Math.abs(r.level - level(sol, silver)[edge]) < 1e-9);
    // Anodic (into the solution, positive) above the level; slopes by decade of current.
    assert.ok(p.current[1] > 0 && p.current[3] < 0 && p.current[0] === r.current[0]);
    const slope = (j, k) => Math.log(Math.abs(p.current[k] / p.current[j])) / Math.abs(p.V[k] - p.V[j]);
    assert.ok(Math.abs(slope(1, 2) / ((1 - alpha) / VT) - 1) < 1e-4, `anodic ${slope(1, 2)}`);
    assert.ok(Math.abs(slope(3, 4) / (alpha / VT) - 1) < 1e-4, `cathodic ${slope(3, 4)}`);
  }
  // The cell's current (toward +x) is what the left face passes into the solution: cathodic
  // there, the right electrode's level being 50 mV up (anodic).
  const left = polarization(dev, sol, { face: 0 }, sol.V['e-'][0]).current[0];
  assert.ok(left < 0 && Math.abs(left / sol.current - 1) < 1e-9, `${left} vs ${sol.current}`);
});

test('at an electrode port: two couples\' curves, the net crossing zero at the floating metal\'s mixed potential', () => {
  const lib = aqueous(['Fe2+', 'H+', 'Cl-'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const dev = new Device(build({
    library: [lib],
    stack: [bath({ 'Fe2+': 10, 'H+': 10, 'Cl-': 30 }, 'Cl-'), layer('water', 10e-6, { name: 'film' }), {}],
    ports: [{ name: 'iron', region: 'film', I: 0, terminal: 'e-', area: 1e3, reactions: [
      { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 1e-7, alpha: 0.5 },
      { equation: '2 H+ + 2 e- = H2', fixed: { H2: 0 }, k0: 1e-9, alpha: 0.5 },
    ] }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  }));
  const sol = dev.solve();
  assert.ok(sol.converged);
  const port = sol.ports[0], w = 4, x = port.x[w], E = port.V;
  const p = polarization(dev, sol, { port: 'iron', x }, [E - 0.05, E, E + 0.05]);
  p.reactions.forEach((r, k) => assert.ok(Math.abs(r.rate[1] / port.rates[k][w] - 1) < 1e-9));
  // Iron's level below, hydrogen's above, and the metal between them, its net current ~0 there
  // (not exactly: the whole window balances, each spot only nearly).
  const [fe, h2] = p.reactions;
  assert.ok(fe.level < E && E < h2.level);
  const g = sol.x.indexOf(x);
  assert.ok(Math.abs(fe.level - sol.V['Fe2+'][g]) < 1e-9, `${fe.level} vs ${sol.V['Fe2+'][g]}`);
  assert.ok(p.current[0] < 0 && p.current[2] > 0 && Math.abs(p.current[1]) < 1e-3 * fe.current[1]);
  assert.ok(fe.current[1] > 0 && h2.current[1] < 0, 'iron dissolving, hydrogen evolving');
});

test('polarization checks where it is asked', () => {
  const dev = new Device(build({ library: [aqueous(['Na+', 'Cl-'], { epsr: 0 })], stack: [bath({ 'Na+': 1, 'Cl-': 1 }, 'Cl-'), layer('water', 1e-6), {}] }));
  const sol = dev.solve();
  assert.throws(() => polarization(dev, sol, {}, 0), (e) => e instanceof DeviceError && /say where/.test(e.message));
  assert.throws(() => polarization(dev, sol, { port: 'x' }, 0), /no port/);
});
