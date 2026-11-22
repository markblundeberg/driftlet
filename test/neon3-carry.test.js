// Neon's third review (21 Nov 2026), flag 1.4: a surface's carried coverages put at the amounts a
// set() keeps (src/carry.js, rescaleSurfaces), on a surface near full.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { check } from '../src/kit.js';

const rel = (a, b) => Math.abs(a / b - 1);
// What the port's surface holds of a species, mol per m² of cross-section (Γ Σ θ·a·v over its window).
const held = (d, sol, name) => sol.ports[0].coverage[name].reduce((t, th, w) => t + 2e-5 * th * sol.ports[0].area[w] * d.model.grid.vol[d.model.grid.x.indexOf(sol.ports[0].x[w])], 0);

// Two adsorbates competing for one porous electrode's sites in a counter-gradient, bare 2.5e-8 to
// 4.2e-7 at −0.25 V, no warnings. Refined, A's amount on the new grid asks 1.3% more of the carried
// nodes and B's less: updated one at a time in listing order, A first couldn't fit in the bare sites
// B still held, its odds' bracket ran to e^1024, and set({ grid }) and check()'s grid item threw
// "coverage would be NaN"; listed B first, both passed.
const competing = (order, V = -0.25) => {
  const sp = { Aads: { mu0: -25000, capacity: 2e-5, theta0: 0.4 }, Bads: { mu0: -18000, capacity: 2e-5, theta0: 0.4 } };
  return {
    species: [{ name: 'e-', z: -1, cRef: 1000 }, { name: 'A+', z: 1, cRef: 1000 }, { name: 'B+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr: 0, species: { 'A+': { D: 1e-9, mu0: 0 }, 'B+': { D: 1e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } } },
    regions: [{ name: 'el', material: 'aq', length: 1e-4, c0: { 'A+': 50, 'B+': 50, 'Cl-': 100 } }],
    contacts: { left: { V: 0, bath: { c: { 'A+': 100, 'B+': 1e-2, 'Cl-': 100.01 } } }, right: { V: 0, bath: { c: { 'A+': 1e-2, 'B+': 100, 'Cl-': 100.01 } } } },
    ports: [{ name: 'cat', terminal: 'e-', region: 'el', from: 0, to: 1e-4, area: 1e5, V, surface: Object.fromEntries(order.map((k) => [k, sp[k]])), reactions: [{ equation: 'B+ + e- = Bads', k0: 1e-4 }, { equation: 'A+ + e- = Aads', k0: 1e-4 }] }],
    grid: { hmin: 5e-7, hmax: 5e-5, ratio: 4 },
  };
};
const refined = { hmin: 1e-6, hmax: 2e-5, ratio: 1.5 };

test("set({ grid }) carries two adsorbates near full in either listing order, to the same state, their amounts kept (Neon R3 1.4)", () => {
  const carried = {};
  for (const order of [['Aads', 'Bads'], ['Bads', 'Aads']]) {
    const d = new Device(competing(order)), s = d.solve();
    assert.ok(s.converged && s.warnings.length === 0 && Math.min(...s.ports[0].bare) < 1e-7, s.warnings.join('; '));
    assert.notEqual(check(d, s).items.find((it) => it.name === 'grid').ok, null, `${order}: grid item solved the refined device`);
    const was = { Aads: held(d, s, 'Aads'), Bads: held(d, s, 'Bads') };
    d.set({ grid: refined });
    const c = d.solution();
    for (const k of ['Aads', 'Bads']) assert.ok(rel(held(d, c, k), was[k]) < 1e-12, `${order}: ${k} ${was[k]} → ${held(d, c, k)}`);
    carried[order.join('')] = c.ports[0];
  }
  // The same carried state: each coverage to 1e-12. The bare fraction is the sites less what's
  // adsorbed, so it's fixed only to the amounts' round-off over its share (1e-16 / 1e-7): 1.4e-9 here.
  const a = carried.AadsBads, b = carried.BadsAads;
  a.x.forEach((x, w) => {
    for (const k of ['Aads', 'Bads']) assert.ok(rel(a.coverage[k][w], b.coverage[k][w]) < 1e-12, `${k} at ${x}: ${a.coverage[k][w]} vs ${b.coverage[k][w]}`);
    assert.ok(rel(a.bare[w], b.bare[w]) < 1e-7, `bare at ${x}: ${a.bare[w]} vs ${b.bare[w]}`);
  });
});

// A porous Pt electrode's adsorbed hydrogen from pH 1 to pH 12 at −0.90 V, bare 1.4e-17 to 1.4e-6:
// its new coverage at x = 0 rounds to 1, and the bare fraction found as 1 − Σθ cancelled to 0, so a
// pure regrid (and check()'s grid item) was refused, "its coverage would be 1.000".
test('set({ grid }) carries a surface whose bare fraction is 1e-17, found from its odds (Neon R3 1.4)', () => {
  const def = {
    species: [{ name: 'e-', z: -1, cRef: 1000 }, { name: 'H+', z: 1, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { aq: { epsr: 0, species: { 'H+': { D: 9.3e-9, mu0: 0 }, 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } } },
    regions: [{ name: 'el', material: 'aq', length: 1e-4, c0: { 'H+': 50, 'K+': 50, 'Cl-': 100 } }],
    contacts: { left: { V: 0, bath: { c: { 'H+': 100, 'K+': 10, 'Cl-': 110 } } }, right: { V: 0, bath: { c: { 'H+': 1e-9, 'K+': 110, 'Cl-': 110 } } } },
    ports: [{ name: 'cat', terminal: 'e-', region: 'el', from: 0, to: 1e-4, area: 1e5, V: -0.9, surface: { Hads: { mu0: -15000, capacity: 2e-5, theta0: 0.5 } }, reactions: [{ equation: 'H+ + e- = Hads', k0: 1e-4, alpha: 0.5 }] }],
  };
  const d = new Device(def), s = d.solve(), bare0 = s.ports[0].bare[0];
  assert.ok(s.converged && bare0 < 1e-16 && bare0 > 0, `bare ${bare0}`);
  assert.notEqual(check(d, s).items.find((it) => it.name === 'grid').ok, null, 'grid item solved the refined device');
  const was = held(d, s, 'Hads');
  d.set({ grid: { hmin: 5e-8, hmax: 2.5e-6 } });
  const c = d.solution(), bare = c.ports[0].bare[0];
  assert.ok(rel(held(d, c, 'Hads'), was) < 1e-12, `Hads ${was} → ${held(d, c, 'Hads')}`);
  // (its level shifted by the 4e-11 more its amount asks of the sites, keeping the moles: 6% less bare)
  assert.ok(bare > 0.8 * bare0 && bare < 1.2 * bare0, `bare at x = 0: ${bare0} → ${bare}`);
});

// The two adsorbates at −0.84 V, bare 2.6e-18 to 4.4e-17 at every node: their wants on the new grid
// came to 1.4e-16 more than all the carried sites, round-off (each amount is good to ~1e-16 of
// itself), and the bare sites' total, the sites less the amounts, is that round-off's. Refused, or
// put wherever round-off pointed, a pure regrid now keeps the amounts and moves the bare sites by
// a few percent at most.
test('set({ grid }) carries two adsorbates full at every node, whose wants exceed the sites by round-off, in either order (Neon R3 1.4)', () => {
  for (const order of [['Aads', 'Bads'], ['Bads', 'Aads']]) {
    const d = new Device(competing(order, -0.84)), s = d.solve();
    assert.ok(s.converged && Math.max(...s.ports[0].bare) < 1e-16);
    const was = { Aads: held(d, s, 'Aads'), Bads: held(d, s, 'Bads') }, b0 = s.ports[0].bare;
    d.set({ grid: refined });
    const c = d.solution(), b1 = c.ports[0].bare;
    for (const k of ['Aads', 'Bads']) assert.ok(rel(held(d, c, k), was[k]) < 1e-12, `${order}: ${k} ${was[k]} → ${held(d, c, k)}`);
    assert.ok(rel(b1[0], b0[0]) < 0.1 && rel(b1[b1.length - 1], b0[b0.length - 1]) < 0.1, `${order}: bare ${b0[0]}…${b0[b0.length - 1]} → ${b1[0]}…${b1[b1.length - 1]}`);
  }
});
