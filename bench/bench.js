// Benchmarks for typical interactive workloads. Run: npm run bench
// Prints wall time per task (best of several runs) and Newton iteration counts.

import { BlockTridiagonal, Device, FARADAY, units } from '../src/index.js';

const now = () => performance.now();
function best(fn, runs = 5) {
  let t = Infinity, out;
  for (let k = 0; k < runs; k++) {
    const t0 = now();
    out = fn();
    t = Math.min(t, now() - t0);
  }
  return { ms: t, out };
}

const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
const diode = () =>
  new Device({
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
    bulkReactions: [{ reactants: { 'e-': 1, 'h+': 1 }, kf: { Si: 1e-6 } }],
    grid: { hmin: 0.5e-9, hmax: 20e-9 },
  });
const electrode = (V) => ({ V, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
const cell = (epsr) =>
  new Device({
    species: [
      { name: 'Ag+', z: 1, cRef: 1000 },
      { name: 'NO3-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
    regions: [{ material: 'water', length: 20e-6, c0: { 'NO3-': 10 } }],
    contacts: { left: electrode(0), right: electrode(0) },
    grid: { hmin: 0.1e-9, hmax: 200e-9, ratio: 1.15 },
  });

const rows = [];
const report = (name, { ms, out }, extra = '') => rows.push([name, ms.toFixed(2), extra || (out?.iterations !== undefined ? `${out.iterations} its` : '')]);

{
  const n = 300, m = 7, sys = new BlockTridiagonal(n, m);
  const fill = () => {
    for (let k = 0; k < sys.B.length; k++) {
      sys.A[k] = Math.sin(k);
      sys.C[k] = Math.cos(k);
      sys.B[k] = (k % (m * m)) % (m + 1) === 0 ? 20 : Math.sin(3 * k);
    }
  };
  const rhs = new Float64Array(n * m).fill(1);
  fill();
  report('linear solve 300 × 7 (factor + solve)', best(() => (sys.factor(), sys.solve(rhs)), 50));
}
{
  const dev = diode();
  report('pn diode: cold equilibrium', best(() => diode().solve()));
  dev.solve();
  let V = 0;
  report('pn diode: warm re-solve, +10 mV', best(() => (dev.set({ contacts: { right: { V: (V += 0.01) } } }), dev.solve()), 10));
  report('pn diode: warm re-solve, 0.4 V → −1 V', best(() => {
    dev.set({ contacts: { right: { V: 0.4 } } });
    dev.solve();
    dev.set({ contacts: { right: { V: -1 } } });
    return dev.solve();
  }, 3));
  dev.set({ contacts: { right: { V: 0 } } });
  dev.solve();
  report('pn diode: advance 0 → 0.5 V, 100 ns', best(() => {
    const d = diode();
    d.solve();
    d.set({ contacts: { right: { V: 0.5 } } });
    return d.advance(1e-7);
  }, 3), '');
  report('pn diode: impedance, 20 frequencies', best(() => dev.impedance(Array.from({ length: 20 }, (_, k) => 10 ** (k / 2)))));
}
{
  const dev = cell(78.5);
  report('Ag|AgNO₃|Ag (ε > 0): cold steady', best(() => cell(78.5).solve()));
  dev.solve();
  let V = 0;
  report('Ag|AgNO₃|Ag (ε > 0): warm re-solve, +5 mV', best(() => (dev.set({ contacts: { right: { V: (V += 0.005) } } }), dev.solve()), 10));
  report('Ag|AgNO₃|Ag (ε = 0): advance 1 s after 50 mV', best(() => {
    const d = cell(0);
    d.solve();
    d.set({ contacts: { right: { V: 0.05 } } });
    return d.advance(1);
  }, 3), '');
}

const w = Math.max(...rows.map((r) => r[0].length));
console.log(`driftlet benchmarks (Node ${process.versions.node})\n`);
for (const [name, ms, extra] of rows) console.log(`${name.padEnd(w)}  ${ms.padStart(9)} ms  ${extra}`);
