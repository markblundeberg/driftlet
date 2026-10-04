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

// The cross-section A(x) (m²), for a device that isn't planar: the conservation laws hold for
// the totals through it, flux × A, and a box holds its volume ∫A dx. A planar device has A = 1
// everywhere, and its currents read as densities.
//   spherical: A = 4π r², cylindrical: A = 2π r (per metre of length), r = r0 + x;
//   a profile: A(x) piecewise linear through the points given, constant beyond them.
function crossSection(geometry) {
  const { type } = geometry;
  if (type === 'spherical' || type === 'cylindrical') {
    const r0 = geometry.r0, sph = type === 'spherical';
    const sq = (x) => (r0 + x) * (r0 + x), cube = (x) => (r0 + x) * (r0 + x) * (r0 + x);
    return {
      at: (x) => (sph ? 4 * Math.PI * sq(x) : 2 * Math.PI * (r0 + x)),
      // ∫A dx and ∫dx/A from a to b
      integral: (a, b) => (sph ? (4 * Math.PI * (cube(b) - cube(a))) / 3 : Math.PI * (sq(b) - sq(a))),
      resistance: (a, b) => (sph ? (1 / (r0 + a) - 1 / (r0 + b)) / (4 * Math.PI) : Math.log((r0 + b) / (r0 + a)) / (2 * Math.PI)),
    };
  }
  const xs = geometry.x, vs = geometry.values, last = xs.length - 1;
  const at = (x) => {
    if (x <= xs[0]) return vs[0];
    if (x >= xs[last]) return vs[last];
    let k = 0;
    while (xs[k + 1] < x) k++;
    return vs[k] + ((vs[k + 1] - vs[k]) * (x - xs[k])) / (xs[k + 1] - xs[k]);
  };
  // Over [a, b], split at the profile's points, where A is linear on each piece.
  const pieces = (a, b, f) => {
    let s = 0, lo = a;
    for (const xk of [...xs.filter((v) => v > a && v < b), b]) {
      s += f(lo, xk, at(lo), at(xk));
      lo = xk;
    }
    return s;
  };
  return {
    at,
    integral: (a, b) => pieces(a, b, (p, q, A1, A2) => ((A1 + A2) / 2) * (q - p)),
    resistance: (a, b) =>
      pieces(a, b, (p, q, A1, A2) => {
        const d = A2 - A1, m = (A1 + A2) / 2;
        return Math.abs(d) < 1e-6 * m ? (q - p) / m : ((q - p) * Math.log(A2 / A1)) / d;
      }),
  };
}

/**
 * Give a grid its cross-section: `area` (A at each node, m²), `segArea` (each segment's
 * effective area, its length over ∫dx/A, so that a flux D·segArea·Δc/h is exact for steady
 * diffusion between the nodes; A at its middle where A vanishes at an end, as at a sphere's
 * centre) and `vol` (each node's box, ∫A dx). Planar (no geometry): A = 1, and nothing changes.
 * @param {ReturnType<typeof buildGrid>} grid
 * @param {{ type: 'planar' } | { type: 'spherical' | 'cylindrical', r0: number } | { type: 'profile', x: ArrayLike<number>, values: ArrayLike<number> }} geometry
 */
export function applyGeometry(grid, geometry) {
  const { nNodes, x, segLength } = grid;
  grid.area = new Float64Array(nNodes).fill(1);
  grid.segArea = new Float64Array(nNodes - 1).fill(1);
  if (!geometry || geometry.type === 'planar') return grid;
  const A = crossSection(geometry);
  for (let g = 0; g < nNodes; g++) grid.area[g] = A.at(x[g]);
  for (let s = 0; s < nNodes - 1; s++) {
    const h = segLength[s];
    if (h === 0) grid.segArea[s] = grid.area[s];
    else if (grid.area[s] > 0 && grid.area[s + 1] > 0) grid.segArea[s] = h / A.resistance(x[s], x[s + 1]);
    else grid.segArea[s] = A.at(x[s] + h / 2);
  }
  for (let g = 0; g < nNodes; g++) {
    const lo = g > 0 ? x[g] - segLength[g - 1] / 2 : x[g], hi = g < nNodes - 1 ? x[g] + segLength[g] / 2 : x[g];
    grid.vol[g] = hi > lo ? A.integral(lo, hi) : 0;
  }
  return grid;
}
