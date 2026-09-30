import { test } from 'node:test';
import assert from 'node:assert/strict';
import { buildGrid, gradedCells } from '../src/grid.js';

const sum = (a) => a.reduce((s, v) => s + v, 0);

test('graded cells fill the region exactly and grade from both ends', () => {
  const L = 1e-3, hmin = 1e-9, hmax = 1e-5, ratio = 1.2;
  const cells = gradedCells(L, { hmin, hmax, ratio });
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
  const cells = gradedCells(2, { minCells: 10 });
  assert.equal(cells.length, 10);
  for (const h of cells) assert.ok(Math.abs(h - 0.2) < 1e-15);
});

test('one-sided refinement', () => {
  const cells = gradedCells(1, { hLeft: 1e-6, hRight: 0.1, hmax: 0.1, ratio: 1.3 });
  assert.ok(cells[0] < 2e-6);
  assert.ok(cells.at(-1) > 0.05);
  assert.ok(Math.abs(sum(cells) - 1) < 1e-14);
});

test('rejects bad inputs', () => {
  assert.throws(() => gradedCells(0, {}), RangeError);
  assert.throws(() => gradedCells(1, { ratio: 0.9 }), RangeError);
  assert.throws(() => gradedCells(1, { hmin: -1 }), RangeError);
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
