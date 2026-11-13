import { test } from 'node:test';
import assert from 'node:assert/strict';
import { buildGrid, taperedCells } from '../src/grid.js';
import { Device } from '../src/index.js';

const sum = (a) => a.reduce((s, v) => s + v, 0);

test('graded cells fill the region exactly and grade from both ends', () => {
  const L = 1e-3, hmin = 1e-9, hmax = 1e-5, ratio = 1.2;
  const cells = taperedCells(L, { hmin, hmax, ratio });
  assert.ok(Math.abs(sum(cells) - L) < 1e-15 * cells.length);
  // Scaling to fit shrinks cells a little, never grows them.
  assert.ok(cells[0] <= hmin * (1 + 1e-12) && cells[0] > 0.5 * hmin);
  assert.ok(cells.at(-1) <= hmin * (1 + 1e-12) && cells.at(-1) > 0.5 * hmin);
  for (const h of cells) assert.ok(h <= hmax * (1 + 1e-12));
  for (let k = 1; k < cells.length; k++) {
    const q = cells[k] / cells[k - 1];
    assert.ok(q <= ratio * (1 + 1e-9) && q >= 1 / ratio / (1 + 1e-9), `cell ${k}: ratio ${q}`);
  }
  // Roughly logarithmic cost: from 1 nm to 10 µm at ratio 1.2 is ~50 cells per end.
  assert.ok(cells.length < 300, `${cells.length} cells`);
});

test('uniform when no hmin is given, with at least minCells cells', () => {
  const cells = taperedCells(2, { minCells: 10 });
  assert.equal(cells.length, 10);
  for (const h of cells) assert.ok(Math.abs(h - 0.2) < 1e-15);
});

test('one-sided refinement', () => {
  const cells = taperedCells(1, { hLeft: 1e-6, hRight: 0.1, hmax: 0.1, ratio: 1.3 });
  assert.ok(cells[0] < 2e-6);
  assert.ok(cells.at(-1) > 0.05);
  assert.ok(Math.abs(sum(cells) - 1) < 1e-14);
});

test('rejects bad inputs', () => {
  assert.throws(() => taperedCells(0, {}), RangeError);
  assert.throws(() => taperedCells(1, { ratio: 0.9 }), RangeError);
  assert.throws(() => taperedCells(1, { hmin: -1 }), RangeError);
});

test('device grid doubles nodes at every region boundary', () => {
  const regions = [{ length: 1e-6 }, { length: 2e-6 }, { length: 5e-7 }];
  const g = buildGrid(regions, { hmin: 1e-9, hmax: 1e-7 });
  assert.equal(g.x.length, g.nNodes);
  assert.equal(g.x[0], 0);
  assert.equal(g.x[g.nNodes - 1], g.length);
  assert.ok(Math.abs(g.length - 3.5e-6) < 1e-21);

  const links = [];
  for (let s = 0; s < g.nNodes - 1; s++) {
    if (g.segRegion[s] < 0) {
      links.push(s);
      assert.equal(g.segLength[s], 0);
      assert.equal(g.x[s], g.x[s + 1], 'doubled nodes share x exactly');
      assert.equal(g.nodeRegion[s], g.segFace[s]);
      assert.equal(g.nodeRegion[s + 1], g.segFace[s] + 1);
    } else {
      assert.ok(g.segLength[s] > 0);
      assert.equal(g.segFace[s], -1);
      assert.equal(g.nodeRegion[s], g.segRegion[s]);
      assert.equal(g.nodeRegion[s + 1], g.segRegion[s]);
      assert.ok(Math.abs(g.x[s + 1] - g.x[s] - g.segLength[s]) < 1e-18);
    }
  }
  assert.deepEqual(links, [g.regionEnd[0], g.regionEnd[1]]);
  assert.deepEqual([...g.regionStart], [0, g.regionEnd[0] + 1, g.regionEnd[1] + 1]);

  // Box volumes add up to the device length, and each region's boxes to its length.
  assert.ok(Math.abs(sum(g.vol) - g.length) < 1e-20);
  for (let r = 0; r < regions.length; r++) {
    let v = 0;
    for (let k = g.regionStart[r]; k <= g.regionEnd[r]; k++) v += g.vol[k];
    assert.ok(Math.abs(v - regions[r].length) < 1e-20);
  }
});

test('per-region grid options override device defaults', () => {
  const g = buildGrid([{ length: 1, grid: { minCells: 4 } }, { length: 1 }], { minCells: 20 });
  assert.equal(g.regionEnd[0] - g.regionStart[0], 4);
  assert.equal(g.regionEnd[1] - g.regionStart[1], 20);
});

// (a field given as undefined took its default's place: the cell count came out NaN, and the
// device refused with a TypeError; set() of { minCells: undefined } is how a patch removes one)
test('a grid field given as undefined is as if not given', () => {
  assert.deepEqual(buildGrid([{ length: 1, grid: { minCells: undefined } }, { length: 1 }], { minCells: undefined, hmax: undefined }).x, buildGrid([{ length: 1 }, { length: 1 }]).x);
  const def = { species: [{ name: 'X', z: 0, cRef: 1 }], materials: { m: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } }, regions: [{ material: 'm', length: 1e-6 }], contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } }, grid: { minCells: 20 } };
  const d = new Device(def);
  assert.equal(d.grid.nNodes, 21);
  d.set({ grid: { minCells: undefined } });
  assert.equal(d.grid.nNodes, 9);
});

test('a grid far past any demo is refused before it is built: a length in µm given as m', () => {
  // 1 m at hmax 20 nm would be 5e7 cells, and exhaust memory long before the solver said a word.
  const t0 = Date.now();
  assert.throws(() => buildGrid([{ name: 'n', length: 1 }, { name: 'p', length: 1 }], { hmin: 1e-9, hmax: 2e-8 }), /at least 1\.0e\+8 cells.*n is 1 m long at hmax 2e-8 m\. Lengths are in metres/);
  assert.ok(Date.now() - t0 < 100);
  assert.equal(buildGrid([{ length: 1e-3 }], { hmax: 2e-8 }).nNodes, 50001); // 5e4 cells is fine
});

// Fable round 5's review of src/device.js: a ratio alone switched the default grading off and
// gave uniform cells, 125 nm where the default was 1 nm at the ends, unwarned.
test('a grid ratio with neither hmin nor hmax is warned of: the cells are uniform', () => {
  const def = (grid) => ({
    species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: { water: { epsr: 78, species: { 'Na+': { D: 1e-9, mu0: 0 }, 'Cl-': { D: 1e-9, mu0: 0 } } } },
    regions: [{ material: 'water', length: 1e-6, c0: { 'Na+': 100, 'Cl-': 100 } }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0 }, right: { V: 0, phi: 'pinned', zeroCharge: 0 } },
    ...(grid ? { grid } : {}),
  });
  const said = (grid) => new Device(def(grid)).model.warnings.filter((w) => /ratio grades between hmin and hmax/.test(w));
  assert.equal(said({ ratio: 1.05 }).length, 1);
  assert.equal(said().length, 0);
  assert.equal(said({ hmin: 1e-9, hmax: 5e-8, ratio: 1.05 }).length, 0);
  assert.equal(said({ minCells: 40 }).length, 0);
});

test('a grid past the cell limit is refused before it is built, from hmin too; minCells is ≥ 1', () => {
  // (the guard counted length/hmax: { ratio: 1, hmin: L/3e6 } built 3,000,001 nodes, and L/1e9 would
  // have tried 1e9; minCells 0 gave one cell of NaN volume, accepted: Fable's review)
  const t0 = Date.now();
  assert.throws(() => taperedCells(1e-6, { ratio: 1, hmin: 1e-6 / 1e9 }), /more than 1e\+6 cells/);
  assert.ok(Date.now() - t0 < 2000);
  assert.throws(() => taperedCells(1e-6, { minCells: 0 }), /minCells must be ≥ 1/);
  // A device's options, each a positive number: NaN reached the cell guard, whose message then
  // failed with a TypeError's text (a smell review's helper).
  const dev = (grid, regionGrid) => new Device({ species: [{ name: 'X', z: 0, cRef: 1 }], materials: { m: { epsr: 0, species: { X: { D: 1e-9, mu0: 0 } } } }, regions: [{ material: 'm', length: 1e-6, grid: regionGrid }], contacts: { left: { species: { X: { type: 'equilibrium', mu: 0 } }, phi: 'neutral' }, right: { phi: 'neutral' } }, grid });
  assert.throws(() => dev({ minCells: NaN }), /grid: minCells must be a positive number of cells, got NaN/);
  assert.throws(() => dev({}, { hmax: -1 }), /regions\[0\]\.grid\.hmax must be a positive number of metres, got -1/);
  assert.throws(() => dev({ ratio: '2' }), /ratio must be a positive number \(≥ 1\), got "2"/);
});
