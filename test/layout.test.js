import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, bath, ohmic as ohmicKit, aqueous, semiconductor, metal, hodgkinHuxley } from '../src/kit.js';

// The Jacobian, assembled straight into blocks that hold only each node's unknowns, against
// central differences of the residual, column by column. Devices are chosen to cover every
// assembly path: semiconductor carriers with recombination, electrolytes with species confined
// to some regions, conductor regions with electrode reactions at their faces, non-ideal
// statistics, advection and mixing, ports, and a floating terminal.

const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
const silicon = { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } };
const carriers = [
  { name: 'e-', z: -1 },
  { name: 'h+', z: 1 },
];
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
const coarse = { hmin: 2e-9, hmax: 100e-9, ratio: 1.5, minCells: 4 };

const ions = [
  { name: 'Ag+', z: 1, cRef: 1000 },
  { name: 'NO3-', z: -1, cRef: 1000 },
];
const water = (epsr) => ({ epsr, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } });
// Ag⁺ + e⁻ ⇌ Ag at a silver face, with the metal on the given side.
const plating = (metal) => ({ [metal]: { 'e-': -1, Ag: 1 }, [metal === 'left' ? 'right' : 'left']: { 'Ag+': -1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.4 });
const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const salt = { 'NO3-': 10, 'Ag+': 10 };

const devices = {
  'pn diode with recombination, and a held port': () => ({
    species: carriers,
    materials: { Si: silicon },
    regions: [
      { material: 'Si', length: 0.5e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { name: 'p', material: 'Si', length: 0.5e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
    ],
    ports: [{ region: 'p', from: 0.2e-6, to: 0.3e-6, V: 0.2, terminal: 'e-', species: { 'e-': 'equilibrium' } }],
    contacts: { left: ohmic(0), right: ohmic(0.3) },
    bulkReactions: [{ nu: { 'e-': -1, 'h+': -1 }, kf: { Si: 1e-6 } }],
    grid: coarse,
  }),
  'n-Si | KCl, species confined to their regions': () => ({
    species: [...carriers, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { Si: silicon, water: { epsr: 78.5, species: { 'K+': { D: 1.96e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } } },
    regions: [
      { material: 'Si', length: 0.3e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'water', length: 0.3e-6, c0: { 'K+': 100, 'Cl-': 100 } },
    ],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, dipole: 0 }],
    contacts: { left: ohmic(0), right: { V: 0.3, bath: { c: { 'K+': 100, 'Cl-': 100 }, reference: 'Cl-' } } },
    grid: coarse,
  }),
  'silver electrodes and a bipolar plate, with reactions at every face': () => {
    const face = (metal) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating(metal)] });
    const Ag = { material: 'Ag', length: 0.5e-6 };
    return {
      species: [...ions, { name: 'e-', z: -1 }],
      materials: { water: water(78.5), Ag: { conductor: { species: 'e-', conductivity: 6e7 } } },
      regions: [Ag, { material: 'water', length: 1e-6, c0: salt }, Ag, { material: 'water', length: 1e-6, c0: salt }, Ag],
      interfaces: [face('left'), face('right'), face('left'), face('right')],
      contacts: { left: collector(0), right: collector(0.3) },
      grid: coarse,
    };
  },
  'lattice-gas electrolyte with flow and mixing, and a port': () => ({
    species: ions,
    materials: { water: { ...water(78.5), statistics: [{ type: 'lattice', species: ['Ag+', 'NO3-'], cMax: 3000 }] } },
    regions: [{ name: 'w', material: 'water', length: 1e-6, c0: salt, velocity: 1e-4, mixing: 1e-9 }],
    ports: [{ region: 'w', from: 0.4e-6, to: 0.6e-6, V: 0.02, terminal: 'Ag+', species: { 'Ag+': { type: 'conductance', G: 1e6 } } }],
    contacts: {
      left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' },
      right: { V: 0.05, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' },
    },
    grid: coarse,
  }),
  'insertion host (neutral-combination statistics)': () => {
    const RT = GAS_CONSTANT * 298.15, x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];
    return {
      species: [
        { name: 'Li+', z: 1 },
        { name: 'e-', z: -1 },
      ],
      materials: {
        host: {
          epsr: 0,
          species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-8, mu0: 0, cRef: 30000 } },
          statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E: x.map((v) => 0.4 - (RT / FARADAY) * Math.log(v / (1 - v))), muRef: 0 } }],
        },
      },
      regions: [{ material: 'host', length: 1e-6 }],
      contacts: {
        left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' },
        right: { V: 0.42, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      },
      grid: coarse,
    };
  },
  'an electrode spread through a port, floating in current mode, its area varying': () => ({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(0) },
    regions: [{ material: 'water', length: 1e-6, c0: salt }],
    contacts: { left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' }, right: {} },
    ports: [{ region: 0, from: 0.3e-6, I: 2, terminal: 'e-', area: { x: [0.3e-6, 1e-6], values: [1e5, 3e6] }, reactions: [{ equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-3, alpha: 0.4 }] }],
    grid: coarse,
  }),
  'spherical shells: a face with a reaction, a gate, a conductance link in current mode': () => ({
    geometry: { type: 'spherical', r0: 0.3e-6 },
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(78.5), Ag: { conductor: { species: 'e-', conductivity: 6.3e7 } } },
    regions: [{ material: 'Ag', length: 0.2e-6 }, { material: 'water', length: 1e-6, c0: salt }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating('left')] }],
    contacts: {
      left: collector(0),
      right: { I: -1e-12, terminal: 'Ag+', species: { 'Ag+': { type: 'conductance', G: 50 } }, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 },
    },
    grid: coarse,
  }),
  'an electrode port in a profiled cross-section': () => ({
    geometry: { area: { x: [0, 1e-6], values: [2, 0.5] } },
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(0) },
    regions: [{ material: 'water', length: 1e-6, c0: salt }],
    contacts: { left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' }, right: {} },
    ports: [{ region: 0, from: 0.3e-6, I: 1, terminal: 'e-', area: 2e6, reactions: [{ equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-3, alpha: 0.4 }] }],
    grid: coarse,
  }),
  'an electrode port with a surface: two adsorbates, dissolution on bare sites': () => ({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(0) },
    regions: [{ material: 'water', length: 1e-6, c0: salt }],
    contacts: { left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' }, right: {} },
    ports: [{
      // (μ° chosen so both coverages sit near 0.3, where every term is well above round-off)
      region: 0, from: 0.3e-6, V: 0, terminal: 'e-', area: 2e6,
      surface: { 'Ag(ads)': { mu0: 1e3, capacity: 1e-5, theta0: 0.2 }, 'O(ads)': { mu0: 58e3, capacity: 1e-5, theta0: 0.1 } },
      reactions: [
        { equation: 'Ag+ + e- = Ag(ads)', k0: 1e-3, alpha: 0.4, bare: true },
        { equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-3, alpha: 0.6, bare: true },
        { equation: 'O(ads) + NO3- = Ag+ + 2 e- + Ag(ads)', k0: 1e-6, alpha: 0.5 },
      ],
    }],
    grid: coarse,
  }),
  'a capacitance spread along a port (ε > 0), floating behind a resistance': () => ({
    species: ions,
    materials: { water: water(78.5) },
    regions: [{ material: 'water', length: 1e-6, c0: salt }],
    contacts: { left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' }, right: {} },
    ports: [{ name: 'gate', region: 0, from: 0.2e-6, V: 0.02, R: 10, area: { x: [0, 1e-6], values: [1e6, 4e6] }, capacitance: { C: 0.1, zeroCharge: 0.05 } }],
    grid: coarse,
  }),
  // Gates: a membrane face's (GHK permeabilities, m³h and n⁴), and two membrane ports' along the
  // inside, one held, one floating behind a resistance.
  'gated channels on a face and on two ports, one floating': () => {
    const T = 279.45, OUT = { 'Na+': 440, 'K+': 20, 'Cl-': 460 }, IN = { 'Na+': 50, 'K+': 400, 'Cl-': 450 };
    const hh = hodgkinHuxley({ inside: IN, outside: OUT, T }), hp = hodgkinHuxley({ area: 1e4, T });
    const def = build({
      T,
      library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })],
      stack: [bath(OUT, 'Cl-'), layer('water', 1e-6, { name: 'out', c0: OUT }), { phi: { type: 'capacitive', C: 0.01 }, ...hh }, layer('water', 2e-6, { name: 'in', c0: IN }), bath(IN, 'Cl-', -0.02)],
      grid: { hmin: 2e-7, hmax: 2e-7 },
    });
    const port = (name, from, to, drive) => ({ name, region: 'in', from, to, area: 1e4, capacitance: { C: 0.01, zeroCharge: 0.01 }, bath: { c: OUT }, gates: hp.gates, species: hp.species, ...drive });
    def.ports = [port('p1', 0.3e-6, 0.9e-6, { V: -0.01 }), port('p2', 1.1e-6, 1.7e-6, { V: 0.005, R: 30 })];
    return def;
  },
  // A port holding a level inside another's window: its row replaces the node's balance, the
  // other's ∂res/∂V there with it.
  'a held level inside a conductance port\'s window': () => {
    const S = { 'Na+': 100, 'K+': 100, 'Cl-': 200 };
    const def = build({ library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })], stack: [bath(S, 'Cl-'), layer('water', 2e-6, { name: 'in', c0: S }), bath(S, 'Cl-')], grid: { hmin: 2e-7, hmax: 2e-7 } });
    def.ports = [
      { name: 'p1', region: 'in', from: 0.3e-6, to: 1.7e-6, V: 0.01, terminal: 'K+', species: { 'K+': { type: 'conductance', G: 1e6 } } },
      { name: 'eq', region: 'in', from: 0.5e-6, to: 0.9e-6, V: 0.003, terminal: 'K+', R: 5, species: { 'K+': 'equilibrium' } },
    ];
    return def;
  },
  // A MOS capacitor's gate as a metal region (its face to the oxide capacitive), the channel
  // grounded by a port holding the electrons' level.
  'a metal-region gate over an oxide, with a channel port': () => ({
    ...build({
      T: 300,
      library: [semiconductor('Si'), metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }],
      stack: [ohmicKit(0.3, ['e-']), layer('Al', 20e-9), { phi: { type: 'capacitive', C: 10 }, zeroCharge: -0.9 }, layer('SiO2', 10e-9), { dipole: 0 }, layer('Si', 200e-9, { name: 'Si', acceptors: units.perCm3(1e17) }), ohmicKit(0)],
      grid: { hmin: 1e-9, hmax: 20e-9, ratio: 1.5 },
    }),
    ports: [{ name: 'inv', region: 'Si', from: 0, to: 3e-9, V: 0.01, terminal: 'e-', species: { 'e-': 'equilibrium' } }],
  }),
  'floating terminal (a conductance link) in current mode': () => ({
    species: ions,
    materials: { water: water(78.5) },
    regions: [{ material: 'water', length: 1e-6, c0: salt }],
    contacts: {
      left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' },
      right: { I: -5, terminal: 'Ag+', species: { 'Ag+': { type: 'conductance', G: 50 } }, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 },
    },
    grid: coarse,
  }),
};

// The Jacobian entry for compact row i and column k (0 outside the three block diagonals).
function entry(sys, blockOf, i, k) {
  const bi = blockOf[i], bk = blockOf[k], r = i - sys.offX[bi], c = k - sys.offX[bk], sz = sys.sizes;
  if (bk === bi) return sys.B[sys.offB[bi] + r * sz[bi] + c];
  if (bk === bi - 1) return sys.A[sys.offA[bi] + r * sz[bk] + c];
  if (bk === bi + 1) return sys.C[sys.offC[bi] + r * sz[bk] + c];
  return 0;
}

for (const [name, make] of Object.entries(devices)) {
  test(`Jacobian (and each terminal's B, C and ∂I/∂V) equals central differences: ${name}`, () => {
    const dev = new Device(make());
    assert.ok(dev.solve().converged);
    const s = dev.solver, N = s.sys.size;
    const blockOf = new Int32Array(N);
    for (let b = 0; b < s.nB; b++) blockOf.fill(b, s.sys.offX[b], s.sys.offX[b + 1]);
    s.allTerminals = true; // (every terminal's ∂res/∂V and ∂I/∂x kept, to check them all)
    // Away from the solution, so that every term is exercised.
    for (let k = 0; k < N; k++) s.u[s.fullOf[k]] += 0.03 * Math.sin(1 + 7 * k);
    for (const dt of [Infinity, 1e-6]) {
      s.computeConcentrations();
      s.cOld.set(s.c);
      for (let k = 0; k < N; k++) s.u[s.fullOf[k]] += 1e-3 * Math.cos(3 * k); // so storage terms aren't zero
      s.assemble(dt);
      const J = { A: s.sys.A.slice(), B: s.sys.B.slice(), C: s.sys.C.slice(), offA: s.sys.offA, offB: s.sys.offB, offC: s.sys.offC, offX: s.sys.offX, sizes: s.sys.sizes };
      const rowMax = new Float64Array(N);
      for (let i = 0; i < N; i++) for (let k = Math.max(0, i - 3 * s.M); k < Math.min(N, i + 3 * s.M); k++) rowMax[i] = Math.max(rowMax[i], Math.abs(entry(J, blockOf, i, k)));
      // Each terminal: ∂res/∂V (B), ∂I/∂x (C) and ∂I/∂V.
      const T = s.terms.length, Bt = s.termB.map((v) => v.slice(0, N)), Ct = s.termC.map((v) => v.slice(0, N)), DI = Float64Array.from(s.termDI);
      const cMax = Ct.map((row) => row.reduce((m, v) => Math.max(m, Math.abs(v)), 0));
      let worst = 0, at = ''; // |FD − J| / (1e-5 |J| + 1e-8 × the row's largest entry)
      const note = (err, where) => {
        if (err > worst) {
          worst = err;
          at = where;
        }
      };
      for (let k = 0; k < N; k++) {
        const f = s.fullOf[k], u0 = s.u[f], h = 1e-6 * Math.max(1, Math.abs(u0));
        s.u[f] = u0 + h;
        s.assemble(dt);
        const rp = s.res.slice(0, N), Ip = Float64Array.from(s.termI);
        s.u[f] = u0 - h;
        s.assemble(dt);
        const rm = s.res.slice(0, N), Im = Float64Array.from(s.termI);
        s.u[f] = u0;
        for (let t = 0; t < T; t++) {
          const fd = (Ip[t] - Im[t]) / (2 * h), an = Ct[t][k];
          note(Math.abs(fd - an) / (1e-5 * Math.abs(an) + 1e-8 * cMax[t] + 1e-300), `terminal ${s.terms[t].name} current, column ${k}: ${an} vs ${fd}`);
        }
        for (let i = 0; i < N; i++) {
          if (Math.abs(blockOf[i] - blockOf[k]) > 1) continue;
          // Relative to the entry, with a floor at the row's scale (finite-difference noise).
          const fd = (rp[i] - rm[i]) / (2 * h), an = entry(J, blockOf, i, k);
          const err = Math.abs(fd - an) / (1e-5 * Math.abs(an) + 1e-8 * rowMax[i] + 1e-300);
          note(err, `row ${i} (block ${blockOf[i]}), column ${k} (block ${blockOf[k]}, slot ${f % s.M}): ${an} vs ${fd}`);
        }
      }
      for (let t = 0; t < T; t++) {
        const fl = s.floating.includes(t), V0 = s.termV[t], h = 1e-6;
        const setV = (V) => (fl ? (s.termV[t] = V) : s.sourceOverride.set(t, V));
        setV(V0 + h);
        s.assemble(dt);
        const rp = s.res.slice(0, N), Ip = s.termI[t];
        setV(V0 - h);
        s.assemble(dt);
        const rm = s.res.slice(0, N), Im = s.termI[t];
        if (fl) s.termV[t] = V0;
        else s.sourceOverride.delete(t);
        for (let i = 0; i < N; i++) {
          const fd = (rp[i] - rm[i]) / (2 * h), an = Bt[t][i];
          note(Math.abs(fd - an) / (1e-5 * Math.abs(an) + 1e-8 * Math.max(rowMax[i], Math.abs(an)) + 1e-300), `terminal ${s.terms[t].name}: ∂res/∂V, row ${i}: ${an} vs ${fd}`);
        }
        const fd = (Ip - Im) / (2 * h), an = DI[t] - (s.terms[t].drive.R > 0 ? 1 / s.terms[t].drive.R : 0);
        note(Math.abs(fd - an) / (1e-5 * Math.abs(an) + 1e-8 * (cMax[t] / s.VT + Math.abs(an)) + 1e-300), `terminal ${s.terms[t].name}: ∂I/∂V: ${an} vs ${fd}`);
      }
      assert.ok(worst < 1, `dt=${dt}: ${worst.toExponential(1)} at ${at}`);
    }
  });
}

// J·v with the dilute kernels' terms in difference form (for GMRES), the rest from a matrix of
// its own: the same product as the assembled Jacobian's, including replaced rows, captured
// terminal rows and the transformed charge rows of strictly neutral nodes on a step.
test('J·v in difference form equals the assembled Jacobian, and a uniform shift gives exactly zero inside a region', () => {
  for (const [name, make] of Object.entries(devices)) {
    const dev = new Device(make());
    assert.ok(dev.solve().converged);
    const s = dev.solver, N = s.sys.size;
    s.allTerminals = true;
    for (let k = 0; k < N; k++) s.u[s.fullOf[k]] += 0.03 * Math.sin(1 + 7 * k);
    for (const [dt, combining] of [[Infinity, false], [1e-6, false], [1e-6, true]]) {
      s.computeConcentrations();
      s.cOld.set(s.c);
      s.combining = combining;
      s.assemble(dt);
      const J = s.sys, C = s.termC.map((v) => v.slice(0, N)), transformed = s.transformed;
      const v = new Float64Array(N + 1).map((_, k) => (k < N ? Math.cos(5 * k) : 0));
      const ref = J.multiply(v);
      const { rest, lin } = s._assembleDifference(dt);
      const out = rest.multiply(v);
      const w = Float64Array.from(v);
      if (transformed) s._untransform(w);
      lin.apply(w, out, 1, 1 / dt);
      s.combining = false;
      let worst = 0;
      for (let i = 0; i < N; i++) {
        worst = Math.max(worst, Math.abs(out[i] - ref[i]) / (Math.abs(ref[i]) + 1e-300));
      }
      assert.ok(worst < 1e-9, `${name}, dt=${dt}${combining ? ', charge rows' : ''}: ${worst}`);
      // (to round-off of the terminal's largest entry: a full assembly's sum of large terms
      // leaves ~1e-15 of it where the difference form's is exactly 0)
      s.termC.forEach((c, t) => {
        const big = C[t].reduce((m, x) => Math.max(m, Math.abs(x)), 0);
        c.slice(0, N).forEach((x, k) => assert.ok(Math.abs(x - C[t][k]) <= 1e-12 * Math.abs(C[t][k]) + 1e-14 * big, `${name}: terminal ${t}, column ${k} (block ${s.fullOf[k] / s.M | 0}, slot ${s.fullOf[k] % s.M}): ${x} vs ${C[t][k]} (largest ${big}), dt=${dt}`));
      });
    }
  }
  // A pn diode: η_e by −s, η_h by s and φ̂ by s at every node changes no flux, charge or storage.
  const dev = new Device(devices['pn diode with recombination, and a held port']());
  dev.solve();
  const s = dev.solver, N = s.sys.size, M = s.M;
  const { lin } = s._assembleDifference(1e-6);
  const e = new Float64Array(N + 1), out = new Float64Array(N);
  for (let k = 0; k < N; k++) {
    const slot = s.fullOf[k] % M;
    e[k] = slot === 0 ? 1 : slot <= s.n ? s.z[slot - 1] : 0;
  }
  lin.apply(e, out, 1, 1e6);
  let inside = 0;
  for (let k = 0; k < N; k++) if (!lin.dead[k]) inside = Math.max(inside, Math.abs(out[k]));
  assert.equal(inside, 0);
});

// After each step, the contact and port readouts come from only the boxes they're read from.
// They must equal a full assembly's exactly.
test('bookkeeping readouts from the end boxes and port windows equal a full assembly, bit for bit', () => {
  for (const [name, make] of Object.entries(devices)) {
    const dev = new Device(make());
    assert.ok(dev.solve().converged);
    const s = dev.solver;
    for (let k = 0; k < s.sys.size; k++) s.u[s.fullOf[k]] += 0.02 * Math.sin(2 + 5 * k);
    for (const dt of [Infinity, 1e-6]) {
      s.computeConcentrations();
      const read = () => ({
        flux: [...s.contactFlux.left, ...s.contactFlux.right],
        D: [s.contactD.left, s.contactD.right],
        ports: s.portFlux.flatMap((p) => [...p]),
        seg: [s.segI, s.segD],
      });
      s.assemble(dt);
      const full = read();
      s._assembleBookkeeping(dt);
      assert.deepEqual(read(), full, `${name}, dt=${dt}`);
    }
  }
});

test('a region held only weakly: its level found through its summed balance, exactly', () => {
  // A strictly neutral region whose one carrier is held only weakly: by tiny conductances at its
  // faces, ~1e17–1e26 times less than its own conduction, or through equilibrium faces by bulk
  // neighbours that conduct ~1e20 times less. Eliminated, its level is lost to round-off, and
  // Newton used to settle at a wrong one silently (a current 2× or 16× off). Summed over the
  // region, its balance is the faces' fluxes alone, exactly; with that row in place of one node's,
  // the direct solve finds the level, and both ends pass the same current.
  const ohm = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const mat = (D) => ({ epsr: 0, species: { 'e-': { D, mu0: 0, cRef: 1e5 } } });
  const dev = ({ G, Dout = 1e-4, Df = 1e-8, V }) => {
    const face = { phi: 'neutral', species: { 'e-': G ? { type: 'conductance', G } : 'equilibrium' } };
    return new Device({
      species: [{ name: 'e-', z: -1 }],
      materials: { out: mat(Dout), fast: mat(Df) },
      regions: [
        { material: 'out', length: 1e-6, fixedCharge: 1e4 * FARADAY },
        { name: 'fast', material: 'fast', length: 1e-6, fixedCharge: 1e4 * FARADAY, grid: { minCells: 200 } },
        { material: 'out', length: 1e-6, fixedCharge: 1e4 * FARADAY },
      ],
      interfaces: [face, face],
      contacts: { left: ohm(0), right: ohm(V) },
    });
  };
  // The two face conductances in series, G/2 per volt; or the two outer regions' resistances,
  // 1 µm each of σ = F²Dc/RT (the fast region's, 1e20 times less, doesn't show).
  const sigma = (FARADAY * FARADAY * 1e-24 * 1e4) / (GAS_CONSTANT * 298.15);
  for (const [p, I] of [
    [{ G: 1e-6, V: 0.3 }, (-0.3 * 1e-6) / 2],
    [{ G: 1e-15, V: 0.01 }, (-0.01 * 1e-15) / 2],
    [{ Dout: 1e-24, Df: 1e-4, V: 0.3 }, (-0.3 * sigma) / 2e-6],
  ]) {
    const sol = dev(p).solve({ continuation: false, maxSteps: 1 });
    assert.ok(sol.converged, JSON.stringify(p));
    for (const [side, s] of [['left', 1], ['right', -1]]) {
      assert.ok(Math.abs((s * sol.terminals[side].current) / I - 1) < 1e-9, `${JSON.stringify(p)} ${side}: ${sol.terminals[side].current}`);
    }
  }
  // Without that row the direct solve fails, and the warning names the region.
  const d = dev({ G: 1e-6, V: 0.3 });
  d.solver.islands = [];
  const direct = d.solve({ continuation: false, maxSteps: 1 });
  assert.equal(direct.converged, false);
  assert.ok(direct.warnings.some((w) => /\(fast\): part of the device is held only weakly/.test(w)), direct.warnings.join('\n'));
});
