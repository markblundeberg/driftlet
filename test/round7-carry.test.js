// What set() carries where two copies of one quantity disagreed (src/carry.js against the solver's
// own sums): items 1, 5 and 8 of the inventory of duplicated quantities (20 Nov 2026), each settled
// by a probe. The old number is in each test's comment; the right one from a second route (the
// state's own charge, a fresh device at the new definition, the concentrations as they were).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { aqueous, bath, build, layer, metal, ohmic } from '../src/kit.js';

const rel = (a, b) => Math.abs(a / b - 1);

// A strictly neutral LiCl bar with a capacitive port over [from, 3 µm], its left contact holding φ
// (a 'bulk' contact: the port's end node there is left to it).
const bar = (from, port = {}) => ({
  species: [{ name: 'Li+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
  materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 }, 'Cl-': { D: 1e-11, mu0: 0 } } } },
  regions: [{ name: 'bar', material: 'solid', length: 10e-6, c0: { 'Li+': 100, 'Cl-': 100 } }],
  contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium', 'Cl-': 'blocked' }, phi: 'bulk' }, right: { phi: 'neutral' } },
  ports: [{ name: 'a', region: 'bar', from, to: 3e-6, V: 0.05, area: 1e6, capacitance: { C: 0.2 }, ...port }],
  grid: { hmax: 0.5e-6 },
});
const flare = { area: { x: [0, 10e-6], values: [1, 3] } };

test("a capacitive port's charge across a new cross-section is what the carried state holds, node by node (no phantom step)", () => {
  // (carry scaled it by one ratio of Σ box·area over every window node, where the capacitance sums
  // σ node by node and leaves an end node to its contact: off by 3.4e-7 C over [1, 3] µm, the first
  // 1 ns step reading −336 A/m²; by 1.9% from 0, +5.59e5 A/m², where −1.08e-2 is right, the gap left
  // in the running charge for good)
  for (const from of [1e-6, 0]) {
    const d = new Device(bar(from));
    d.solve();
    d.set({ ports: { a: { V: 0.1 } } }).advance(1e-3);
    d.step(1e-9);
    d.set({ geometry: flare });
    const s = d.solver;
    s._assembleBookkeeping(Infinity);
    assert.ok(rel(s.portQEnd[0], s.portQ[0]) < 1e-12, `from ${from}: carried ${s.portQEnd[0]}, the state holds ${s.portQ[0]}`);
    const q = s.termCharge[2];
    d.step(1e-9);
    // (what's left, ~4e-11 C/m², is φ settling at the new cross-section's current balance in a
    // strictly neutral bar: the same at 1e-10 and 1e-11 s)
    assert.ok(Math.abs(s.termCharge[2] - q) < 1e-9, `from ${from}: the first step passed ${s.termCharge[2] - q} C/m²`);
  }
  // A charged port left floating (I = 0) keeps its charge per area of electrode, the state holding
  // it, so its voltage stays where it was through the first step. (One ratio, a total no node held,
  // and the first 1 ns step drove the gap through the bar: the port jumped 43 mV from 1 µm and
  // 23 mV from 0; held only before neutrality moved φ under it, 0.58 mV from 0.)
  for (const from of [1e-6, 0]) {
    const f = new Device(bar(from));
    f.solve();
    f.set({ ports: { a: { V: 0.1 } } }).advance(1e-3);
    f.set({ ports: { a: { I: 0 } } }).advance(2e-3);
    const V0 = f.solver.termV[2];
    f.set({ geometry: flare });
    const s = f.solver;
    s._assembleBookkeeping(Infinity);
    assert.ok(rel(s.portQEnd[0], s.portQ[0]) < 1e-12, `floating from ${from}: carried ${s.portQEnd[0]}, the state holds ${s.portQ[0]}`);
    const V1 = f.step(1e-9).ports[0].V;
    assert.ok(Math.abs(V1 - V0) < 1e-4, `floating from ${from}: ${V0} V, after the first step ${V1}`);
  }
});

test('a resized region beside diffuse layers keeps its bulk concentrations (each layer keeps its excess per area)', () => {
  // (carry scaled the whole amount, the layers' Γ·A with it, by the region's volume: a 10% longer
  // KCl film moved its bulk K⁺ by +1.67e-5 at 100 mM, 0.5 V, and +6.0e-3 at 1 mM, 1 V)
  const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 } };
  for (const [V, c] of [[0.5, 100], [1, 1]]) {
    const def = build({ library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')], stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), stern, layer('water', 10e-6, { c0: { 'K+': c, 'Cl-': c } }), stern, layer('Pt', 1e-6), ohmic(0, ['e-'])], grid: { hmin: 0.5e-6, hmax: 0.5e-6 } });
    const d = new Device(def);
    d.solve();
    const a = d.step(1e-6);
    d.set({ regions: d.def.regions.map((r) => (r.material === 'water' ? { ...r, length: r.length * 1.1 } : r)) });
    const b = d.solution(), mid = (x) => x.c['K+'][Math.floor(x.c['K+'].length / 2)];
    assert.ok(rel(mid(b), mid(a)) < 1e-12, `${c} mM at ${V} V: bulk K+ ${mid(a)} → ${mid(b)}`);
  }
});

test('a region in a new cross-section keeps its concentrations where they vary (the amount scaled node by node)', () => {
  // (one ratio of volumes for the region: a polarised Cl⁻, closed, in a bar made 1 → 3 moved by a
  // uniform factor until it held a total no node held)
  const def = bar(1e-6);
  delete def.ports;
  def.contacts = { left: { I: 5, terminal: 'Li+', species: { 'Li+': 'equilibrium', 'Cl-': 'blocked' }, phi: 'bulk' }, right: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium', 'Cl-': 'blocked' }, phi: 'neutral' } };
  const d = new Device(def);
  d.solve();
  d.advance(1e-2);
  const a = d.solution();
  const spread = Math.max(...a.c['Cl-']) / Math.min(...a.c['Cl-']) - 1;
  assert.ok(spread > 1e-3, `Cl- polarised by ${spread}`);
  d.set({ geometry: flare });
  const b = d.solution();
  let worst = 0;
  for (let g = 0; g < a.x.length; g++) for (const sp of ['Li+', 'Cl-']) worst = Math.max(worst, rel(b.c[sp][g], a.c[sp][g]));
  assert.ok(worst < 1e-10, `concentrations moved by up to ${worst}`);
});

test('a metal made an electrolyte starts from its c0, as its solution says, not with the metal\'s sheet charge on its ions', () => {
  // (the old net charge was read before the kind test, and a strictly neutral region kept it: a Pt
  // made water held K⁺ 100.10 / Cl⁻ 99.90 at its edge)
  const lib = [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')], C = { 'K+': 100, 'Cl-': 100 };
  const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }, grid = { hmin: 0.5e-6, hmax: 0.5e-6 };
  const A = build({ library: lib, stack: [bath(C, 'Cl-'), layer('water', 5e-6, { c0: C }), stern, layer('Pt', 1e-6), stern, layer('water', 5e-6, { c0: C }), bath(C, 'Cl-', 0.1)], grid });
  const B = build({ library: lib, stack: [bath(C, 'Cl-'), layer('water', 5e-6, { c0: C }), {}, layer('water', 1e-6, { c0: C }), {}, layer('water', 5e-6, { c0: C }), bath(C, 'Cl-', 0.1)], grid });
  const d = new Device(A);
  d.solve();
  d.step(1e-6);
  d.set({ materials: B.materials, regions: B.regions, interfaces: B.interfaces, contacts: B.contacts });
  const s = d.solution(), gr = d.solver.model.grid;
  assert.ok(s.warnings.some((w) => /region 1.*c0/.test(w)), s.warnings.join('\n'));
  for (let g = gr.regionStart[1]; g <= gr.regionEnd[1]; g++) for (const sp of ['K+', 'Cl-']) assert.ok(rel(s.c[sp][g], 100) < 1e-12, `${sp} ${s.c[sp][g]} at x = ${s.x[g]}`);
});
