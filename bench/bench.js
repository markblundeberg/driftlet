// Benchmarks for typical interactive workloads.
//
//   npm run bench               run, and compare with bench/baseline.json
//   npm run bench -- --save     run, and write the results as the new baseline
//   npm run bench -- --check    run, and exit non-zero if any counter got worse (for CI)
//
// Each case reports wall time (best of several runs) and deterministic counters: the number of
// matrix factorisations (≈ Newton iterations, over every step and frequency), and the block
// work Σ m³ summed over them, the leading cost of block elimination. The counters are the same
// on every machine, so they're what CI checks; times are from whichever machine saved the
// baseline, and only reported. A change that improves a counter re-saves the baseline in the
// same commit, so the diff shows the effect.

import { readFileSync, writeFileSync } from 'node:fs';
import { cpus } from 'node:os';
import { BlockTridiagonal, Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { ComplexBlockTridiagonal } from '../src/blockTridiagonal.js';

const args = new Set(process.argv.slice(2));
const baselineUrl = new URL('./baseline.json', import.meta.url);

// --- counters, by wrapping the solvers' factor()
const count = { factor: 0, work: 0 };
const blockWork = (sys) => (sys.sizes ? sys.sizes.reduce((s, m) => s + m ** 3, 0) : sys.n * sys.m ** 3);
for (const Cls of [BlockTridiagonal, ComplexBlockTridiagonal]) {
  const factor = Cls.prototype.factor;
  Cls.prototype.factor = function () {
    count.factor++;
    count.work += blockWork(this);
    return factor.call(this);
  };
}

// --- devices
const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
const silicon = { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } };
const carriers = [
  { name: 'e-', z: -1 },
  { name: 'h+', z: 1 },
];
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const diode = () =>
  new Device({
    species: carriers,
    materials: { Si: silicon },
    regions: [
      { material: 'Si', length: 2e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'Si', length: 2e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
    ],
    contacts: { left: ohmic(0), right: ohmic(0) },
    bulkReactions: [{ reactants: { 'e-': 1, 'h+': 1 }, kf: { Si: 1e-6 } }],
    grid: { hmin: 0.5e-9, hmax: 20e-9 },
  });

// Au gate | 5 nm SiO₂ | p-Si
const mos = () =>
  new Device({
    species: carriers,
    materials: { Si: silicon, SiO2: { epsr: 3.9, species: {} }, Au: { metal: { species: 'e-', conductivity: 4e7 } } },
    regions: [
      { material: 'Au', length: 20e-9 },
      { material: 'SiO2', length: 5e-9, grid: { minCells: 4 } },
      { material: 'Si', length: 0.5e-6, fixedCharge: -units.perCm3(1e17) * FARADAY },
    ],
    interfaces: [{ phi: { type: 'capacitive', C: 50 }, zeroCharge: -0.9 }, { dipole: 0 }],
    contacts: { left: collector(0), right: ohmic(0) },
    grid: { hmin: 0.5e-9, hmax: 20e-9, ratio: 1.1 },
  });

const ions = [
  { name: 'Ag+', z: 1, cRef: 1000 },
  { name: 'NO3-', z: -1, cRef: 1000 },
];
const water = (epsr) => ({ epsr, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } });
const electrode = (V) => ({ V, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
const cell = (epsr) =>
  new Device({
    species: ions,
    materials: { water: water(epsr) },
    regions: [{ material: 'water', length: 20e-6, c0: { 'NO3-': 10 } }],
    contacts: { left: electrode(0), right: electrode(0) },
    grid: { hmin: 0.1e-9, hmax: 200e-9, ratio: 1.15 },
  });

// A floating silver plate between two plating electrodes in AgNO₃
const plating = { reactants: { 'Ag+': 1 }, electrons: 1, products: { Ag: 1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 };
const salt = { 'NO3-': 10, 'Ag+': 10 };
const bipolar = () => {
  const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 };
  const end = { V: 0, phi: { type: 'capacitive', C: 0.2, zeroCharge: 0.1 }, reactions: [plating] };
  return new Device({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(78.5), Ag: { metal: { species: 'e-', conductivity: 6e7 } } },
    regions: [
      { material: 'water', length: 10e-6, c0: salt },
      { material: 'Ag', length: 2e-6 },
      { material: 'water', length: 10e-6, c0: salt },
    ],
    interfaces: [{ ...stern, reactions: [plating] }, { ...stern, reactions: [plating] }],
    contacts: { left: end, right: end },
    grid: { hmin: 0.05e-9, hmax: 100e-9, ratio: 1.1 },
  });
};

// A lithium insertion host against an ideal-solution OCV
const insertion = () => {
  const RT = GAS_CONSTANT * 298.15, x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];
  const E = x.map((v) => 0.4 - (RT / FARADAY) * Math.log(v / (1 - v)));
  return new Device({
    species: [
      { name: 'Li+', z: 1 },
      { name: 'e-', z: -1 },
    ],
    materials: {
      host: {
        epsr: 0,
        species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } },
        statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }],
      },
    },
    regions: [{ material: 'host', length: 1e-6 }],
    contacts: {
      left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' },
      right: { V: 0.4, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
    },
  });
};

const sweep = (dev, side, Vs) => {
  for (const V of Vs) {
    dev.set({ contacts: { [side]: { V } } });
    if (!dev.solve().converged) throw new Error(`no convergence at ${side} V = ${V}`);
  }
};
const range = (a, b, n) => Array.from({ length: n }, (_, k) => a + ((b - a) * k) / (n - 1));
const decades = (n, from = 0) => Array.from({ length: n }, (_, k) => 10 ** (from + k / 2));

// Each case: optional setup (not timed), then run(state).
const cases = [
  {
    name: 'linear solve 300 × 7 (factor + solve)',
    runs: 50,
    setup() {
      const n = 300, m = 7, sys = new BlockTridiagonal(n, m);
      for (let k = 0; k < sys.B.length; k++) {
        sys.A[k] = Math.sin(k);
        sys.C[k] = Math.cos(k);
        sys.B[k] = (k % (m * m)) % (m + 1) === 0 ? 20 : Math.sin(3 * k);
      }
      return { sys, rhs: new Float64Array(n * m).fill(1) };
    },
    run: ({ sys, rhs }) => (sys.factor(), sys.solve(rhs)),
  },
  { name: 'pn diode: cold equilibrium', run: () => diode().solve() },
  {
    name: 'pn diode: I–V sweep 0 → 0.6 V, 31 points',
    setup: () => (d => (d.solve(), d))(diode()),
    run: (d) => sweep(d, 'right', range(0, 0.6, 31)),
  },
  {
    name: 'pn diode: warm jump 0.4 V → −1 V',
    setup: () => (d => (d.solve(), d))(diode()),
    run: (d) => sweep(d, 'right', [0.4, -1]),
  },
  {
    name: 'pn diode: advance 0 → 0.5 V, 100 ns',
    runs: 3,
    setup: () => (d => (d.solve(), d.set({ contacts: { right: { V: 0.5 } } }), d))(diode()),
    run: (d) => {
      d.set({ contacts: { right: { V: 0 } } });
      d.solve();
      d.set({ contacts: { right: { V: 0.5 } } });
      return d.advance(1e-7);
    },
  },
  {
    name: 'pn diode: impedance, 20 frequencies',
    setup: () => (d => (d.solve(), d))(diode()),
    run: (d) => d.impedance(decades(20)),
  },
  {
    name: 'MOS (Au gate): C–V sweep −1 → 1.5 V, 26 points',
    setup: () => (d => (d.solve(), d))(mos()),
    run: (d) => sweep(d, 'left', range(-1, 1.5, 26)),
  },
  { name: 'Ag|AgNO₃|Ag (ε > 0): cold steady', run: () => cell(78.5).solve() },
  {
    name: 'Ag|AgNO₃|Ag (ε > 0): sweep 0 → 0.1 V, 21 points',
    setup: () => (d => (d.solve(), d))(cell(78.5)),
    run: (d) => sweep(d, 'right', range(0, 0.1, 21)),
  },
  {
    name: 'Ag|AgNO₃|Ag (ε > 0): impedance, 20 frequencies',
    setup: () => (d => (d.solve(), d))(cell(78.5)),
    run: (d) => d.impedance(decades(20, -2)),
  },
  {
    name: 'Ag|AgNO₃|Ag (ε = 0): advance 1 s after 50 mV',
    runs: 3,
    run: () => {
      const d = cell(0);
      d.solve();
      d.set({ contacts: { right: { V: 0.05 } } });
      return d.advance(1);
    },
  },
  {
    name: 'bipolar Ag electrode: sweep 0 → 1 V, 11 points',
    runs: 3,
    run: () => sweep(bipolar(), 'right', range(0, 1, 11)),
  },
  {
    name: 'insertion host: OCV sweep 0.55 → 0.25 V, 31 points',
    setup: () => (d => (d.solve(), d))(insertion()),
    run: (d) => sweep(d, 'right', range(0.55, 0.25, 31)),
  },
];

// --- run
const results = {};
for (const c of cases) {
  // The first run is counted; the timed runs follow, each from a fresh setup.
  let state = c.setup?.();
  count.factor = count.work = 0;
  c.run(state);
  const counters = { factorizations: count.factor, work: count.work };
  let ms = Infinity;
  for (let k = 0; k < (c.runs ?? 5); k++) {
    state = c.setup?.();
    const t0 = performance.now();
    c.run(state);
    ms = Math.min(ms, performance.now() - t0);
  }
  results[c.name] = { ms: +ms.toFixed(3), ...counters };
}

// --- report
let baseline = null;
try {
  baseline = JSON.parse(readFileSync(baselineUrl, 'utf8'));
} catch {}
const base = baseline?.results ?? {};
const fmtWork = (w) => (w >= 1e6 ? `${(w / 1e6).toFixed(2)}M` : w >= 1e3 ? `${(w / 1e3).toFixed(1)}k` : String(w));
const delta = (now, then) => (then === undefined ? '' : now === then ? '' : ` (${now > then ? '+' : ''}${(((now - then) / then) * 100).toFixed(0)}%)`);
const rows = [['case', 'ms', 'factorizations', 'block work']];
const worse = [], better = [];
for (const [name, r] of Object.entries(results)) {
  const b = base[name];
  rows.push([name, r.ms.toFixed(2) + delta(r.ms, b?.ms), r.factorizations + delta(r.factorizations, b?.factorizations), fmtWork(r.work) + delta(r.work, b?.work)]);
  if (b) {
    for (const k of ['factorizations', 'work']) {
      if (r[k] > b[k]) worse.push(`${name}: ${k} ${b[k]} → ${r[k]}`);
      if (r[k] < b[k]) better.push(`${name}: ${k} ${b[k]} → ${r[k]}`);
    }
  }
}
const widths = rows[0].map((_, j) => Math.max(...rows.map((row) => row[j].length)));
console.log(`driftlet benchmarks (Node ${process.versions.node}, ${cpus()[0]?.model ?? 'unknown CPU'})`);
if (baseline) console.log(`compared with the baseline from Node ${baseline.node}, ${baseline.cpu} (times are only comparable on the same machine)`);
console.log();
for (const row of rows) console.log(row.map((s, j) => (j === 0 ? s.padEnd(widths[j]) : s.padStart(widths[j]))).join('  '));
console.log();

for (const name of Object.keys(base)) if (!(name in results)) console.log(`(no longer run: ${name})`);
for (const name of Object.keys(results)) if (baseline && !(name in base)) console.log(`(new, not in the baseline: ${name})`);
if (better.length) console.log(`Improved (re-save the baseline with --save):\n  ${better.join('\n  ')}`);
if (worse.length) console.log(`Worse than the baseline:\n  ${worse.join('\n  ')}`);

if (args.has('--save')) {
  writeFileSync(baselineUrl, JSON.stringify({ node: process.versions.node, cpu: cpus()[0]?.model ?? 'unknown', results }, null, 2) + '\n');
  console.log(`Saved ${baselineUrl.pathname}`);
}
if (args.has('--check') && (worse.length || !baseline)) {
  console.log(baseline ? 'Counters regressed.' : 'No baseline to check against.');
  process.exitCode = 1;
}
