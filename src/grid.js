// 1D grid with graded spacing and doubled nodes at region boundaries.
//
// Each region gets its own run of nodes from its left end to its right end, so every region
// boundary appears twice: the last node of region r and the first node of region r+1 sit at
// the same x, joined by a zero-length "link" segment. Ordinary segments therefore always lie
// inside a single region, and all material steps happen across links.
//
// Spacing grows geometrically (ratio `ratio`) from `hmin` at each region end, capped at
// `hmax`, then is scaled so the cells fill the region exactly. A region may instead give its
// cell widths outright (`cells`), as a metal region does: one cell, exact for Ohm's law.

const DEFAULTS = Object.freeze({ ratio: 1.2, minCells: 8 });

/**
 * Cell widths for one region: fine at both ends, growing toward the middle.
 * @param {number} length m
 * @param {import('./types.js').GridOptions} [opts]
 * @returns {number[]} widths summing to `length`
 */
export function gradedCells(length, { hmin, hmax, ratio = DEFAULTS.ratio, minCells = DEFAULTS.minCells, hLeft, hRight } = {}) {
  if (!(length > 0)) throw new RangeError(`region length must be > 0, got ${length}`);
  if (!(ratio >= 1)) throw new RangeError(`grid ratio must be ≥ 1, got ${ratio}`);
  const cap = Math.min(hmax ?? length / minCells, length / minCells);
  let hl = Math.min(hLeft ?? hmin ?? cap, cap);
  let hr = Math.min(hRight ?? hmin ?? cap, cap);
  if (!(hl > 0 && hr > 0)) throw new RangeError('grid spacing must be > 0');

  const left = [], right = [];
  let sum = 0;
  // The tolerance stops round-off in the running sum from adding a sliver cell.
  while (sum < length * (1 - 1e-12)) {
    if (hl <= hr) {
      left.push(hl);
      sum += hl;
      hl = Math.min(hl * ratio, cap);
    } else {
      right.push(hr);
      sum += hr;
      hr = Math.min(hr * ratio, cap);
    }
  }
  const cells = left.concat(right.reverse());
  const scale = length / sum;
  for (let k = 0; k < cells.length; k++) cells[k] *= scale;
  return cells;
}

/**
 * Build the device grid.
 * @param {{length: number, grid?: object, cells?: number[]}[]} regions in order, left to right
 * @param {object} [opts] defaults for every region: { hmin, hmax, ratio, minCells }
 */
export function buildGrid(regions, opts = {}) {
  if (!Array.isArray(regions) || regions.length === 0) throw new Error('grid needs at least one region');

  const perRegion = regions.map((reg) => reg.cells ?? gradedCells(reg.length, { ...opts, ...(reg.grid ?? {}) }));
  const nNodes = perRegion.reduce((s, cells) => s + cells.length + 1, 0);
  const nSeg = nNodes - 1;

  const x = new Float64Array(nNodes);
  const nodeRegion = new Int32Array(nNodes);
  const segLength = new Float64Array(nSeg);
  const segRegion = new Int32Array(nSeg); // −1 for a link across a region boundary
  const segFace = new Int32Array(nSeg); // face index for links, −1 otherwise
  const vol = new Float64Array(nNodes); // box volume per unit area
  const regionStart = new Int32Array(regions.length);
  const regionEnd = new Int32Array(regions.length); // inclusive

  let node = 0, x0 = 0;
  for (let r = 0; r < regions.length; r++) {
    const cells = perRegion[r];
    regionStart[r] = node;
    let xi = x0;
    for (let k = 0; k <= cells.length; k++) {
      x[node] = xi;
      nodeRegion[node] = r;
      if (k < cells.length) {
        segLength[node] = cells[k];
        segRegion[node] = r;
        segFace[node] = -1;
        vol[node] += cells[k] / 2;
        vol[node + 1] += cells[k] / 2;
        xi += cells[k];
        node++;
      }
    }
    // Pin the region's right end exactly, so the doubled node's x matches bit-for-bit.
    x0 += regions[r].length;
    x[node] = x0;
    regionEnd[r] = node;
    if (r < regions.length - 1) {
      segLength[node] = 0;
      segRegion[node] = -1;
      segFace[node] = r;
      node++;
    }
  }

  return {
    nNodes,
    x,
    nodeRegion,
    segLength,
    segRegion,
    segFace,
    vol,
    regionStart,
    regionEnd,
    length: x0,
  };
}
