// Round 4, the steady solve (Astra's fourth pass, chased by probes): each test failed before its fix.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous } from '../src/kit.js';

// A9 (= Sol's B2): an affinity pin at a face with a diffuse layer held the edge's levels, where the
// kinetics react at the plane: a closed couple beside a reaction drawing O⁺ across the layer.
const drawn = (k1) => ({
  T: 298.15,
  species: [{ name: 'e-', z: -1 }, { name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }, { name: 'O+', z: 1 }, { name: 'R1', z: 0 }, { name: 'R2', z: 0 }],
  materials: {
    M: { conductor: { species: 'e-', conductivity: 1e7 } },
    w: { epsr: 0, species: Object.fromEntries(['K+', 'Cl-', 'O+', 'R1', 'R2'].map((s) => [s, { D: 1e-9, cRef: 1, mu0: 0 }])) },
  },
  regions: [{ name: 'metal', material: 'M', length: 1e-7 }, { name: 'sol', material: 'w', length: 100e-9, c0: { 'K+': 100, 'Cl-': 101, 'O+': 1, R1: 0.1, R2: 0.1 } }],
  interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: -0.45, diffuse: { epsr: 78 }, reactions: [{ equation: 'O+ + e- = R1', k0: k1 }, { equation: 'O+ + e- = R2', k0: 1 }] }],
  contacts: { left: { V: -0.05, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, bath: { c: { 'K+': 100, 'Cl-': 101, 'O+': 1, R2: 0.1 } } } },
  grid: { minCells: 10 },
});
test('a pinned face reaction beside a diffuse layer is at equilibrium at the plane, as the transient settles', () => {
  for (const k1 of [1e-6, 1]) {
    const s = new Device(drawn(k1)).solve(), g = s.region.indexOf(1);
    const t = new Device(drawn(k1));
    let a;
    for (const dt of [1, 1e3, 1e6, 1e8]) a = t.step(dt);
    assert.ok(s.converged && s.warnings.length === 0, s.warnings.join());
    // (the edge's pin gave 1.3052 at k1 = 1e-6, 3.165 at 1, against 0.99766)
    assert.ok(Math.abs(s.c.R1[g] / a.c.R1[g] - 1) < 1e-9, `k1 ${k1}: ${s.c.R1[g]} vs ${a.c.R1[g]}`);
  }
});

// A2: a port reaction's affinities were stored in an Int32Array, 0.3 RT read as 0.
test('a frozen port reaction a fraction of RT off is warned of', () => {
  const cb = 60.42185;
  const def = (port, k0) => ({
    species: ['K+', 'X-', 'M+', 'e-'].map((name) => ({ name, z: name === 'K+' || name === 'M+' ? 1 : -1, cRef: 1000 })),
    materials: { film: { epsr: 0, species: { 'K+': { D: 2e-9, mu0: -283e3 }, 'X-': { D: 2e-9, mu0: -131e3 }, 'M+': { D: 1e-9, mu0: 70e3 } } } },
    regions: [{ name: 'film', material: 'film', length: 36e-6, c0: { 'K+': 60.4, 'X-': 70.4, 'M+': 10 } }],
    contacts: { left: { bath: { c: { 'K+': cb, 'X-': cb }, reference: 'X-' }, V: 0 }, right: { phi: 'neutral' } },
    ports: [{ name: 'el', region: 'film', terminal: 'e-', ...port, area: 1e5, reactions: [{ equation: 'M+ + e- = M(s)', fixed: { 'M(s)': 0 }, k0 }] }],
    grid: { hmin: 1e-6, hmax: 6e-6 },
  });
  const V0 = new Device(def({ I: 0 }, 4e-7)).solve().terminals.el.V, VT = (GAS_CONSTANT * 298.15) / FARADAY;
  const d = new Device(def({ V: V0 }, 1e-30));
  d.solve();
  d.set({ ports: { el: { V: V0 - 0.3 * VT } } });
  const s = d.solve();
  assert.ok(s.warnings.some((w) => /ports\[0\]\.reactions\[0\]: out of equilibrium/.test(w)), JSON.stringify(s.warnings));
});

// A3: a reaction split into two identical channels: each was what the other's check left in.
test('a frozen reaction split into identical channels is warned of, as one is', () => {
  const cell = (kfs) => new Device({
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'A-', z: -1, cRef: 1000 }, { name: 'KA', z: 0, cRef: 1000 }],
    materials: { w: { epsr: 0, species: { 'K+': { D: 1e-9, mu0: 0 }, 'A-': { D: 1e-9, mu0: 0 }, KA: { D: 1e-9, mu0: -5e3 } } } },
    regions: [{ material: 'w', length: 1e-4, c0: { 'K+': 100, 'A-': 100, KA: 1 } }],
    bulkReactions: kfs.map((kf) => ({ equation: 'K+ + A- = KA', kf: { w: kf } })),
    contacts: { left: { V: 0, terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, species: {}, phi: 'neutral' } },
  });
  for (const kfs of [[1e-30], [5e-31, 5e-31]]) {
    const s = cell(kfs).solve();
    assert.ok(s.warnings.some((w) => /bulkReactions\[0\]: out of equilibrium by 111 mV/.test(w)), `${kfs}: ${s.warnings}`);
  }
  // (and solved where it can be, split or not)
  const Q = (s) => (s.c.KA[5] * 1000) / (s.c['K+'][5] * s.c['A-'][5]);
  assert.ok(Math.abs(Q(cell([0.5, 0.5]).solve()) / Q(cell([1]).solve()) - 1) < 1e-9);
});

// A4: +1 A/m² of K⁺ in at a contact and −1 out at a port, one K⁺ pool between them: refused as
// charging, where the transient settles.
test('currents that balance into one conserved pool have a steady state', () => {
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 });
  const make = (IP) => {
    const d = new Device({
      species: lib.species, materials: lib.materials, regions: [{ material: 'water', length: 100e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
      contacts: { left: { phi: 'neutral', species: { 'K+': { type: 'conductance', G: 1, offset: 0 }, 'Cl-': 'blocked' }, I: 0 }, right: { V: 0, terminal: 'Cl-', phi: 'neutral', species: { 'Cl-': 'equilibrium', 'K+': 'blocked' } } },
      ports: [{ name: 'p', region: 0, from: 20e-6, to: 40e-6, species: { 'K+': { type: 'conductance', G: 1e6, offset: 0 } }, I: 0 }],
      grid: { hmin: 2e-6, hmax: 5e-6 },
    });
    d.solve();
    d.set({ contacts: { left: { I: 1 } }, ports: { p: { I: IP } } });
    return d;
  };
  const s = make(-1).solve();
  const t = make(-1);
  let a;
  for (const dt of [0.01, 1, 100, 1e4, 1e4, 1e5]) a = t.step(dt);
  assert.ok(s.converged, s.warnings.join());
  assert.ok(Math.abs(s.terminals.left.V - a.terminals.left.V) < 1e-9, `${s.terminals.left.V} vs ${a.terminals.left.V}`);
  const u = make(-0.5), su = u.solve(); // (unbalanced: still refused)
  assert.ok(!su.converged && u.solver.charging?.terminal === 'left');
});

// A5: a throw inside solve() (a custom statistics' table ending) left the state moved.
test('a solve that throws leaves the device as it was', () => {
  const d = new Device({ species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } }, statistics: [{ type: 'custom', species: ['X'], evaluate: ([z]) => {
      if (z > 0.5) throw new Error('table ends');
      return { c: [Math.exp(z)], dcdzeta: [[Math.exp(z)]] };
    } }] } },
    regions: [{ name: 's', material: 'f', length: 100e-6, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: 0 } } }, right: { phi: 'neutral', species: { X: 'blocked' } } },
    grid: { minCells: 20 } });
  d.solve();
  d.set({ contacts: { left: { species: { X: { type: 'equilibrium', mu: 2000 } } } } });
  const u0 = Float64Array.from(d.solver.u), t0 = d.solver.time;
  assert.throws(() => d.solve(), /table ends/);
  assert.deepEqual(Array.from(d.solver.u), Array.from(u0));
  assert.equal(d.solver.time, t0);
  assert.doesNotThrow(() => d.solution());
});

// A6: an electron-only port over a region switched off the immobile check there for every species.
test('a port acting on other species leaves the frozen-reaction check of an immobile one', () => {
  const dev = (port) => new Device({
    T: 300,
    species: [{ name: 'e-', z: -1 }, { name: 'Y', z: 0 }, { name: 'Z', z: 0 }],
    materials: {
      m1: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 1e-9, mu0: 0, cRef: 1e3 }, Z: { D: 1e-9, mu0: 0, cRef: 1e3 } } },
      m2: { epsr: 10, species: { 'e-': { D: 1e-4, mu0: 0, cRef: 1e3 }, Y: { D: 0, mu0: 0, cRef: 1e3 }, Z: { D: 1e-9, mu0: 0, cRef: 1e3 } } },
    },
    regions: [
      { name: 'A', material: 'm1', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 1, Z: 0.01 } },
      { name: 'B', material: 'm2', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 0.01, Z: 1 } },
    ],
    interfaces: [{ dipole: 0 }],
    bulkReactions: [{ equation: 'Y = Z', kf: { m1: 1e-30, m2: 1e-30 } }],
    contacts: { left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
    ...(port ? { ports: [{ name: 'g', region: 'B', terminal: 'e-', V: 0, species: { 'e-': { type: 'conductance', G: 1, offset: 0 } } }] } : {}),
    grid: { minCells: 40 },
  });
  for (const port of [false, true]) {
    const s = dev(port).solve();
    assert.ok(s.warnings.some((w) => /bulkReactions\[0\]: out of equilibrium/.test(w)), `port ${port}: ${s.warnings}`);
  }
});
