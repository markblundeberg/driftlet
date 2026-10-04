import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';

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
    materials: { water: water(78.5) },
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
    materials: { water: water(78.5) },
    regions: [{ material: 'water', length: 1e-6, c0: salt }],
    contacts: { left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' }, right: {} },
    ports: [{ region: 0, from: 0.3e-6, I: 1, terminal: 'e-', area: 2e6, reactions: [{ equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-3, alpha: 0.4 }] }],
    grid: coarse,
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

test('a failed steady solve says where the system is nearly singular', () => {
  // A strictly neutral region whose one fast carrier is held only by tiny conductances at its
  // faces: in η form its level is held by G against internal conductances ~1e17 times larger,
  // past what doubles can resolve. (A conductor region uses a mixed form that avoids this.)
  const ohm = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const mat = (D) => ({ epsr: 0, species: { 'e-': { D, mu0: 0, cRef: 1e5 } } });
  const weak = { phi: 'neutral', species: { 'e-': { type: 'conductance', G: 1e-6 } } };
  const sol = new Device({
    species: [{ name: 'e-', z: -1 }],
    materials: { out: mat(1e-4), fast: mat(1e-8) },
    regions: [
      { material: 'out', length: 1e-6, fixedCharge: 1e4 * FARADAY },
      { name: 'fast', material: 'fast', length: 1e-6, fixedCharge: 1e4 * FARADAY, grid: { minCells: 200 } },
      { material: 'out', length: 1e-6, fixedCharge: 1e4 * FARADAY },
    ],
    interfaces: [weak, weak],
    contacts: { left: ohm(0), right: ohm(0.01) },
  }).solve();
  assert.equal(sol.converged, false);
  assert.ok(sol.warnings.some((w) => /\(fast\): part of the device is held only weakly/.test(w)), sol.warnings.join('\n'));
});
