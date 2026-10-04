// Hard cases: solves and transients that have failed, or nearly, at some point. Each reports
// whether it got there, the steps taken and rejected, the matrix factorisations, and a physical
// result to keep an eye on (a change to the numerics shouldn't move it). Failures are recorded,
// not thrown: a case that fails today is a target, and the baseline shows when it's fixed.
//
//   npm run hard                  run, and compare with bench/hard.json
//   node bench/hard.js --save     run, and write the results as the new baseline
//   node bench/hard.js name…      run only the cases whose names include one of these

import { readFileSync, writeFileSync } from 'node:fs';
import { BlockTridiagonal, Device, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, bath, aqueous, metal, half, semiconductor, photogeneration } from '../src/kit.js';

const argv = process.argv.slice(2);
const save = argv.includes('--save'), only = argv.filter((a) => !a.startsWith('--'));
const baselineUrl = new URL('./hard.json', import.meta.url);

let factorizations = 0;
const factor = BlockTridiagonal.prototype.factor;
BlockTridiagonal.prototype.factor = function () {
  factorizations++;
  return factor.call(this);
};

const T = 298.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;
const MU_H2O = -237.129e3;

// A run's tally: every advance and solve through these, so a case reports its whole effort.
class Tally {
  steps = 0;
  rejected = 0;
  ok = true;
  note = '';
  advance(dev, t, opts) {
    const s = dev.advance(t, opts);
    this.steps += s.steps ?? 0;
    this.rejected += s.rejected ?? 0;
    if (!s.converged || !s.done) this.fail(`stopped at t = ${s.time.toPrecision(3)} s`);
    return s;
  }
  solve(dev, opts) {
    const s = dev.solve(opts);
    if (!s.converged) this.fail('steady solve failed');
    return s;
  }
  fail(why) {
    if (this.ok) this.note = why;
    this.ok = false;
  }
}

// --- liquid junctions grown from a boundary (strictly neutral water): CsF | NaCl
const junctionIons = { Na: [1, 50.11], Cs: [1, 77.3], F: [-1, 55.4], Cl: [-1, 76.31] };
const junction = (w) => {
  const T = 293.15, names = Object.keys(junctionIons), L = 400e-6;
  const D = (n) => (junctionIons[n][1] * 1e-4 * GAS_CONSTANT * T) / (FARADAY * FARADAY);
  const cl = { Na: 10, Cl: 10, Cs: 135, F: 135 }, cr = { Na: 145, Cl: 145, Cs: 1e-6, F: 1e-6 };
  const profile = (n) => (w > 0 ? { x: [0, L / 2 - w / 2, L / 2 + w / 2, L], values: [cl[n], cl[n], cr[n], cr[n]] } : undefined);
  const water = { epsr: 0, species: Object.fromEntries(names.map((n) => [n, { D: D(n), mu0: 0 }])) };
  const regions = w > 0
    ? [{ material: 'water', length: L, c0: Object.fromEntries(names.map((n) => [n, profile(n)])) }]
    : [{ material: 'water', length: L / 2, c0: cl }, { material: 'water', length: L / 2, c0: cr }];
  return new Device({
    T,
    species: names.map((n) => ({ name: n, z: junctionIons[n][0], cRef: 1000 })),
    materials: { water },
    regions,
    contacts: { left: { bath: { c: cl, reference: 'Cl' } }, right: { bath: { c: cr, reference: 'Cl' }, I: 0 } },
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
};
const junctionCase = (w, tol) => (tally) => {
  const dev = junction(w);
  const s = tally.advance(dev, 0.5, { tol });
  return { 'Δφ (mV)': 1000 * (s.phi.at(-1) - s.phi[0]) };
};

// --- an electrode spread through strictly neutral NaOH, its surface filling with OH⁻
const coverageCell = () => {
  const lib = aqueous(['Na+', 'OH-', 'Fe2+'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const L = 10e-6;
  return new Device(build({
    T,
    library: [lib],
    stack: [{}, layer('water', L, { name: 'film' }), bath({ 'Na+': 100, 'OH-': 102, 'Fe2+': 1 }, 'OH-')],
    ports: [{ name: 'metal', region: 'film', V: 0.05, terminal: 'e-', area: 1e4,
      surface: { OHads: { mu0: 0, capacity: 1e-5, theta0: 1e-4 } },
      reactions: [{ equation: 'OH- = OHads + e-', k0: 1e-4, alpha: 0.5, bare: true }] }],
    grid: { hmin: L / 20, hmax: L / 20 },
  }));
};
const coverageCase = (t, tol) => (tally) => {
  const s = tally.advance(coverageCell(), t, { tol });
  return { θ: s.ports[0].coverage.OHads[5] };
};

// --- a silver face in strictly neutral AgNO₃, 50 mV from equilibrium at the start
const silverCell = () =>
  new Device(build({
    library: [aqueous(['Ag+', 'NO3-'], { epsr: 0 }), metal('Ag')],
    stack: [ohmic(0.05, ['e-']), layer('Ag', 1e-6), { reactions: [{ ...half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 }), k0: 1e-3, alpha: 0.5 }] },
      layer('water', 10e-6), bath({ 'Ag+': 10, 'NO3-': 10 }, 'NO3-')],
    grid: { hmin: 0.5e-6, hmax: 0.5e-6 },
  }));
const silverCase = (t, tol) => (tally) => ({ 'I (A/m²)': tally.advance(silverCell(), t, { tol }).current });

// --- a drop of salt water on iron: O₂ from the air through an exchange port, reduced at the
// two faces while the iron dissolves; the water's potential is set only by the faces.
const dropOnFaces = () => {
  const water = aqueous(['Na+', 'Cl-', 'Fe2+', 'OH-'], { epsr: 0 });
  water.species.push({ name: 'O2', z: 0, cRef: 1000 });
  water.materials.water.species.O2 = { D: 2e-9, mu0: 16.4e3 };
  const iron = { species: [{ name: 'e-', z: -1 }], materials: { Fe: { conductor: { species: 'e-', conductivity: 1e7 } } } };
  const face = { reactions: [
    { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 5e-6, alpha: 0.5 },
    { equation: 'O2 + 2 H2O + 4 e- = 4 OH-', fixed: { H2O: MU_H2O }, k0: 3e-9, alpha: 0.125 },
  ] };
  const Lw = 1e-3, cO2 = 0.26;
  return new Device(build({
    library: [water, iron],
    stack: [ohmic(0, ['e-']), layer('Fe', 1e-4), face, layer('water', Lw, { name: 'drop', c0: { 'Na+': 500, 'Cl-': 500, 'Fe2+': 1e-6, 'OH-': 2e-6, O2: cO2 } }), face, layer('Fe', 1e-4), ohmic(0, ['e-'])],
    ports: [{ name: 'air', region: 'drop', from: Lw / 2, species: { O2: { type: 'exchange', k: 0.05, mu: 16.4e3 + RT * Math.log(cO2 / 1000) } } }],
    grid: { hmin: 5e-6, hmax: 50e-6 },
  }));
};
const dropOnFacesCase = (tally) => {
  const dev = dropOnFaces();
  tally.advance(dev, 1);
  const s = tally.advance(dev, 1000);
  return { 'I_Fe (A/m²)': -2 * FARADAY * (s.interfaces[0].rates[0] + s.interfaces[1].rates[0]) };
};

// --- a drop on iron, radially (its height and the iron's area per volume varying with r): the
// iron an electrode spread through the water, dissolving on bare sites and passivating where the
// pH rises, O₂ reduced on it. Without a double layer the mixed potential has to jump when the last
// active patch passivates (no continuous solution); with one it's continuous.
const radialDrop = (C) => {
  const R = 3e-3, h0 = R / 3, cO2 = 0.26, DO2 = 2e-9;
  const lib = aqueous(['Na+', 'Cl-', 'Fe2+', 'OH-', 'H+'], { epsr: 0 });
  lib.species.push({ name: 'O2', z: 0, cRef: 1000 }, { name: 'e-', z: -1 });
  lib.materials.water.species.O2 = { D: DO2, mu0: 16.4e3 };
  const xs = [], A = [], a = [], k = [];
  for (let j = 0; j <= 60; j++) {
    const r = (j / 60) * R, h = Math.max(h0 * (1 - (r / R) ** 2), 0.02 * h0);
    xs.push(r), A.push(2 * Math.PI * r * h), a.push(1 / h), k.push((DO2 * cO2) / (h * h));
  }
  const zc = -0.4;
  return new Device({
    T, species: lib.species, materials: lib.materials, geometry: { area: { x: xs, values: A } },
    regions: [{ name: 'drop', material: 'water', length: R, c0: { 'Na+': 500, 'Cl-': 500 + 2e-6, 'Fe2+': 1e-6, 'OH-': 1e-4, 'H+': 1e-4, O2: cO2 } }],
    contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
    bulkReactions: [
      { equation: 'H+ + OH- = H2O', fixed: { H2O: MU_H2O }, kf: { water: 1.4e8 } },
      { equation: 'Air = O2', fixed: { Air: 16.4e3 + RT * Math.log(cO2 / 1000) }, kf: { water: { x: xs, values: k } } },
    ],
    ports: [{ name: 'iron', region: 'drop', terminal: 'e-', V: 0, area: { x: xs, values: a },
      ...(C ? { capacitance: { C, zeroCharge: zc } } : {}),
      surface: { 'Fe(OH)2': { mu0: 2 * -157244 + 2 * RT * Math.log(1e-3) + 2 * FARADAY * zc, capacity: 2e-5, theta0: 1e-6 } },
      reactions: [
        { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 5e-6, alpha: 0.5, bare: true },
        { equation: 'Fe(s) + 2 OH- = Fe(OH)2 + 2 e-', fixed: { 'Fe(s)': 0 }, k0: 1e-7, alpha: 0.5, bare: true },
        { equation: 'O2 + 2 H2O + 4 e- = 4 OH-', fixed: { H2O: MU_H2O }, k0: 3e-9, alpha: 0.125 },
      ] }],
    grid: { hmin: 5e-6, hmax: 60e-6, ratio: 1.2 },
  });
};
const radialDropCase = (C) => (tally) => {
  const dev = radialDrop(C);
  let s;
  for (const t of [60, 600, 3600]) if (tally.ok) s = tally.advance(dev, t);
  return { 'E − φ₀ (V)': s.ports[0].V - s.phi[0] };
};

// --- a MOS capacitor whose inversion electrons can only come from the back contact
const mosDef = () => {
  const Lsi = units.um(0.5), tox = units.nm(5);
  return build({
    T: 300,
    library: [semiconductor('Si'), metal('Au')],
    materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [
      ohmic(0),
      layer('Si', Lsi, { name: 'p-Si', acceptors: units.perCm3(1e17) }),
      { dipole: 0 },
      layer('SiO2', tox, { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }),
      { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.1 },
      layer('Au', units.nm(50)),
      ohmic(0, ['e-']),
    ],
    grid: { hmin: units.nm(0.1), hmax: units.nm(10), ratio: 1.1 },
  });
};
const Cox = (3.9 * EPS0) / units.nm(5);
const gateCharge = (s) => s.trace.t.reduce((q, t, k) => q + s.trace.current[k] * (t - (s.trace.t[k - 1] ?? 0)), 0);
const mosCases = {
  'MOS without a channel port: cold steady solves, 0 → 1 V in 0.1 V': (tally) => {
    let failed = [];
    for (let V = 0; V <= 1.0001; V += 0.1) {
      const d = new Device(mosDef());
      d.set({ contacts: { right: { V } } });
      if (!d.solve().converged) failed.push(V.toFixed(1));
    }
    if (failed.length) tally.fail(`failed at ${failed.join(', ')} V`);
    return { failures: failed.length };
  },
  'MOS without a channel port: warm steady solves, 0 → 1 V in 0.05 V': (tally) => {
    const d = new Device(mosDef());
    tally.solve(d);
    for (let V = 0.05; V <= 1.0001 && tally.ok; V += 0.05) {
      d.set({ contacts: { right: { V } } });
      if (!d.solve().converged) tally.fail(`failed at ${V.toFixed(2)} V`);
    }
    return {};
  },
  'MOS without a channel port: 10 mV gate step at 0.6 V, advance 1 s': (tally) => {
    const d = new Device(mosDef());
    d.set({ contacts: { right: { V: 0.6 } } });
    tally.solve(d);
    d.set({ contacts: { right: { V: 0.61 } } });
    return { 'ΔQ/(C_ox·δV)': Math.abs(gateCharge(tally.advance(d, 1))) / (Cox * 0.01) };
  },
  'MOS without a channel port: 10 mV gate step at 0.6 V, advance 1000 s': (tally) => {
    const d = new Device(mosDef());
    d.set({ contacts: { right: { V: 0.6 } } });
    tally.solve(d);
    d.set({ contacts: { right: { V: 0.61 } } });
    const s = tally.advance(d, 1000, { maxSteps: 20000 });
    return { 'ΔQ/(C_ox·δV)': Math.abs(gateCharge(s)) / (Cox * 0.01) };
  },
};

// --- the solar demo's n⁺p cell at open circuit in dim infrared light, solved cold
const solarCase = (hmin) => (tally) => {
  const Si = semiconductor('Si'), We = units.um(0.5), Wb = units.um(60), NA = units.perCm3(1e16);
  const alpha = 1630; // 1/m at 1050 nm (Green 2008)
  const dev = new Device(build({
    T: 300,
    library: [Si],
    stack: [ohmic(0), layer('Si', We, { name: 'n⁺', donors: units.perCm3(1e19) }), layer('Si', Wb, { name: 'p', acceptors: NA }), (({ V, ...c }) => ({ ...c, I: 0 }))(ohmic(0))],
    bulkReactions: [
      { equation: 'e- + h+ = 0', kf: { Si: 1 / (60e-6 * NA) } },
      photogeneration({ material: 'Si', flux: 3e-3 * 0.01, alpha, mu: (6.62607015e-34 * 299792458 * 6.02214076e23) / 1050e-9, to: We + Wb }),
    ],
    grid: { hmin, hmax: units.um(1) },
  }));
  return { 'V_oc (V)': tally.solve(dev).terminals.right.V };
};

const cases = {
  'liquid junction from a sharp boundary (ε = 0), tol 1e-4': junctionCase(0, 1e-4),
  'liquid junction from a sharp boundary (ε = 0), tol 1e-3': junctionCase(0, 1e-3),
  'liquid junction from a 1 µm boundary (ε = 0), tol 1e-4': junctionCase(1e-6, 1e-4),
  'liquid junction from a 1 µm boundary (ε = 0), tol 1e-6': junctionCase(1e-6, 1e-6),
  'electrode port filling its surface (ε = 0): to 1e-4 s, tol 1e-8': coverageCase(1e-4, 1e-8),
  'electrode port filling its surface (ε = 0): to 1e-4 s, tol 1e-6': coverageCase(1e-4, 1e-6),
  'electrode port filling its surface (ε = 0): to 10 s, tol 1e-3': coverageCase(10, 1e-3),
  'silver face in AgNO₃ (ε = 0): to 1e-7 s, tol 1e-8': silverCase(1e-7, 1e-8),
  'silver face in AgNO₃ (ε = 0): to 1e-3 s, tol 1e-8': silverCase(1e-3, 1e-8),
  'silver face in AgNO₃ (ε = 0): to 1e-3 s, tol 1e-3': silverCase(1e-3, 1e-3),
  'drop on iron faces, O₂ through an exchange port (ε = 0): to 1000 s': dropOnFacesCase,
  'radial drop on a passivating iron electrode, double layer 0.2 F/m²: to 3600 s': radialDropCase(0.2),
  'radial drop on a passivating iron electrode, no double layer: to 3600 s': radialDropCase(0),
  ...mosCases,
  'solar n⁺p, dim 1050 nm, open circuit from cold, hmin 1 nm': solarCase(units.nm(1)),
  'solar n⁺p, dim 1050 nm, open circuit from cold, hmin 0.25 nm': solarCase(units.nm(0.25)),
};

// --- run
const results = {};
for (const [name, run] of Object.entries(cases)) {
  if (only.length && !only.some((o) => name.includes(o))) continue;
  const tally = new Tally();
  factorizations = 0;
  const t0 = performance.now();
  let value = {};
  try {
    value = run(tally) ?? {};
  } catch (e) {
    tally.fail(`threw: ${e.message.split('\n')[0].slice(0, 80)}`);
  }
  const ms = performance.now() - t0;
  results[name] = { ok: tally.ok, note: tally.note, steps: tally.steps, rejected: tally.rejected, factorizations, ms: Math.round(ms), value };
  const fmt = (v) => (typeof v === 'number' ? (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e4 ? v.toExponential(4) : v.toPrecision(6)) : String(v));
  const vals = Object.entries(value).map(([k, v]) => `${k} ${fmt(v)}`).join(', ');
  console.log(`${tally.ok ? 'ok  ' : 'FAIL'} ${name}\n     ${tally.steps} steps, ${tally.rejected} rejected, ${factorizations} factorizations, ${Math.round(ms)} ms${vals ? `; ${vals}` : ''}${tally.note ? `; ${tally.note}` : ''}`);
}

let baseline = null;
try {
  baseline = JSON.parse(readFileSync(baselineUrl, 'utf8'));
} catch {}
if (baseline) {
  const changes = [];
  for (const [name, r] of Object.entries(results)) {
    const b = baseline[name];
    if (!b) changes.push(`new: ${name}`);
    else if (b.ok !== r.ok) changes.push(`${r.ok ? 'now passes' : 'now FAILS'}: ${name}`);
    else if (r.factorizations !== b.factorizations) changes.push(`${name}: factorizations ${b.factorizations} → ${r.factorizations}`);
  }
  console.log(changes.length ? `\nAgainst the baseline:\n  ${changes.join('\n  ')}` : '\nSame as the baseline.');
}
if (save) {
  writeFileSync(baselineUrl, JSON.stringify(only.length && baseline ? { ...baseline, ...results } : results, null, 2) + '\n');
  console.log(`Saved ${baselineUrl.pathname}`);
}
