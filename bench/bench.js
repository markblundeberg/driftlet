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
import { hodgkinHuxley, injector } from '../src/kit.js';

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
    bulkReactions: [{ nu: { 'e-': -1, 'h+': -1 }, kf: { Si: 1e-6 } }],
    grid: { hmin: 0.5e-9, hmax: 20e-9 },
  });

// Au gate | 5 nm SiO₂ | p-Si
const mos = () =>
  new Device({
    species: carriers,
    materials: { Si: silicon, SiO2: { epsr: 3.9, species: {} }, Au: { conductor: { species: 'e-', conductivity: 4e7 } } },
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

// A floating silver plate between two silver electrodes in AgNO₃, plating at every face
const plating = (metal) => ({ [metal]: { 'e-': -1, Ag: 1 }, [metal === 'left' ? 'right' : 'left']: { 'Ag+': -1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 });
const salt = { 'NO3-': 10, 'Ag+': 10 };
const bipolar = () => {
  const face = (metal) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating(metal)] });
  const Ag = { material: 'Ag', length: 1e-6 };
  return new Device({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(78.5), Ag: { conductor: { species: 'e-', conductivity: 6e7 } } },
    regions: [Ag, { material: 'water', length: 10e-6, c0: salt }, { material: 'Ag', length: 2e-6 }, { material: 'water', length: 10e-6, c0: salt }, Ag],
    interfaces: [face('left'), face('right'), face('left'), face('right')],
    contacts: { left: collector(0), right: collector(0) },
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

// n-Si against KCl across a Helmholtz layer: a Mott–Schottky junction, four species
const msj = () =>
  new Device({
    species: [...carriers, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: {
      Si: silicon,
      water: { epsr: 78.5, species: { 'K+': { D: 1.96e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } },
    },
    regions: [
      { material: 'Si', length: 1e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'water', length: 1e-6, c0: { 'K+': 100, 'Cl-': 100 } },
    ],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, dipole: 0 }],
    contacts: { left: ohmic(0), right: { bath: { c: { 'K+': 100, 'Cl-': 100 }, reference: 'Cl-' } } },
    grid: { hmin: 0.05e-9, hmax: 20e-9, ratio: 1.1 },
  });

// Haynes–Shockley: holes injected for 0.5 µs through a current-driven port into an n-Si bar
// under 30 V/cm, then drifting, spreading and recombining (τ = 20 µs). `injection` is the peak
// Δp at the emitter relative to n₀, roughly.
const pulse = (injection) => {
  const ND = units.perCm3(1e15), L = 3e-3, E = 3000, tp = 0.5e-6, I = (injection * ND * FARADAY * 50e-6) / tp;
  return new Device({
    T: 300,
    species: carriers,
    materials: { Si: silicon },
    regions: [{ material: 'Si', length: L, fixedCharge: FARADAY * ND }],
    contacts: { left: ohmic(E * L), right: ohmic(0) },
    ports: [{ name: 'emitter', region: 0, from: 0.5e-3, to: 0.55e-3, terminal: 'h+', species: { 'h+': { type: 'conductance', G: 1e6 } },
      I: { t: [0, 1e-9, tp, tp + 1e-9], values: [0, I, I, 0] } }],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1 / (20e-6 * ND) } }],
    grid: { hmin: 10e-6, hmax: 10e-6 },
  });
};
// A squid axon, 2 cm of it, in the four-ion cable of the propagation demo: Hodgkin and Huxley's
// channels on a membrane port (three gates per node, so blocks of 8), stimulated at one end.
const axon = (I = 0) => {
  const IN = { 'Na+': 50, 'K+': 400, 'Cl-': 52.5, 'A-': 397.5 }, OUT = { 'Na+': 437, 'K+': 20, 'Cl-': 556 };
  const z = { 'Na+': 1, 'K+': 1, 'Cl-': -1, 'A-': -1 }, a = 238e-6, hh = hodgkinHuxley({ T: 291.65, area: 2 / a });
  const sealed = { species: Object.fromEntries(Object.keys(z).map((s) => [s, 'blocked'])), phi: 'neutral' };
  return {
    T: 291.65,
    species: Object.entries(z).map(([name, zz]) => ({ name, z: zz })),
    materials: { axoplasm: { epsr: 0, species: Object.fromEntries(Object.keys(z).map((s) => [s, { D: s === 'A-' ? 1e-11 : 1.5e-9, mu0: 0, cRef: 1000 }])) } },
    regions: [{ name: 'axon', material: 'axoplasm', length: 0.02, c0: IN }],
    contacts: { left: sealed, right: sealed },
    ports: [
      { name: 'membrane', region: 'axon', V: 0, area: 2 / a, capacitance: { C: 0.01, zeroCharge: 0 }, bath: { c: OUT }, gates: hh.gates, species: hh.species },
      injector({ name: 'stim', region: 'axon', from: 0, to: 1e-3, species: 'K+', I: { t: [0.02, 0.02 + 1e-6, 0.0202, 0.0202 + 1e-6], values: [0, I, I, 0] } }),
    ],
    grid: { hmin: 2e-4, hmax: 2e-4 },
  };
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
    name: 'Haynes–Shockley pulse in n-Si (weak, 5% of n₀): advance 12 µs',
    runs: 3,
    setup: () => (d => (d.solve(), d))(pulse(0.05)),
    run: (d) => d.advance(12e-6),
  },
  {
    name: 'Haynes–Shockley pulse in n-Si (strong, 100% of n₀): advance 12 µs',
    runs: 3,
    setup: () => (d => (d.solve(), d))(pulse(1)),
    run: (d) => d.advance(12e-6),
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
  {
    name: 'n-Si | KCl: Mott–Schottky sweep −0.5 → 1 V, 31 points',
    setup: () => (d => (d.solve(), d))(msj()),
    run: (d) => sweep(d, 'right', range(-0.5, 1, 31)),
  },
  {
    name: 'n-Si | KCl: impedance, 20 frequencies',
    setup: () => (d => (d.set({ contacts: { right: { V: 0.5 } } }), d.solve(), d))(msj()),
    run: (d) => d.impedance(decades(20)),
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
    name: 'squid axon (HH, 100 nodes): a spike, advance 3 ms',
    runs: 3,
    setup: () => {
      const d = new Device(axon());
      d.advance(0.02, { tol: 1e-4 });
      d.set(axon(30));
      return d;
    },
    run: (d) => d.advance(0.023, { tol: 1e-4, dtMax: 1e-5 }),
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
