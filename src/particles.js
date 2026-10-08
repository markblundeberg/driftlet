// Particles: a solution drawn as dots, each a fixed amount of one species, hopping between the
// cells of a display lattice. The hops' one-way rates are set from the solution so that, on
// average, the dots' density is its concentration and their net crossing rate is its flux: the
// dots are a stochastic sample of the drift–diffusion solution, not a separate model. Dots are
// made and unmade where the flux diverges (reactions, generation, ports), and they enter and
// leave through the contacts and any metal.
//
// Across each cell boundary, the one-way fluxes F⁺ (toward +x) and F⁻ are set by two numbers:
// their difference, the solution's net flux J there, and their product, F⁺F⁻ = F₀², the
// boundary's exchange. In equilibrium F⁺ = F⁻ = F₀, detailed balance, so the dots settle into
// the concentration profile whatever F₀ is; F₀ sets only how fast they jitter. Between cells it's
// Scharfetter–Gummel's (whose one-way fluxes over a uniform field have F⁺F⁻ fixed by the drift
// across the boundary): F₀ = (D/h)√(c_K c_K+1)·(d/2)/sinh(d/2), d the step in ψ = μ̄/RT − ln c,
// which gives a dot in a flat potential the diffusivity D. Across a face that a species crosses
// by permeability, the one-way fluxes are the face's own (Ussing's unidirectional fluxes), so the
// dots crossing each way are counted as a tracer experiment would count them.
//
// Where a species is dense, dots would crowd: a cell that would hold more than `cap` of them is
// drawn as a sea instead (its concentration, shaded), which holds no dots but exchanges them with
// its neighbours as a contact's reservoir does: a dot hopping in is absorbed, and dots come out at
// the one-way rate. Crossings between two seas through a face (an ion through a membrane between
// two baths) are events, at the same rates. So a junction shows its majority carriers as seas
// and its minority carriers as dots, each one injected, drifting, diffusing and recombining.

import { GAS_CONSTANT } from './constants.js';

const isObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);

// (d/2)/sinh(d/2) = √(B(d)B(−d)).
const shape = (d) => {
  const a = Math.abs(d) / 2;
  if (a < 1e-4) return 1 - (a * a) / 6;
  if (a > 700) return 0;
  return a / Math.sinh(a);
};

// What lies to one side of a cell boundary, for a species.
const NONE = 0, DOT = 1, OUTSIDE = 2; // nothing (closed); a cell of dots; a sea, a contact's reservoir or a metal

/**
 * A swarm of dots sampling a solution's concentrations and fluxes. Amounts are per area (a planar
 * device's). `step(dt)` moves them on in the device's time, at the rates of the solution last
 * given; `update(sol)` gives a new one.
 * @param {object} sol a solution (from solve(), advance(), live(), …)
 * @param {object} [opts]
 * @param {string[]} [opts.species] which species (default: every species with a concentration)
 * @param {number} [opts.dots] about how many dots per species (default 200), for the default weight
 * @param {number|object} [opts.weight] amount per dot, mol/m²: one number for every species (so
 *   their dots compare), or { name: weight }; default, each species' total over `dots`
 * @param {number} [opts.cap] the most dots a cell holds on average before it's a sea (default ∞)
 * @param {number} [opts.cells] about how many cells across the device (default 120)
 * @param {() => number} [opts.random] a uniform random source in [0, 1) (default Math.random)
 */
export function particles(sol, { species, dots = 200, weight, cap = Infinity, cells = 120, random = Math.random } = {}) {
  const names = species ?? sol.species.map((s) => s.name).filter((nm) => sol.c[nm].some((v) => v > 0));
  for (const nm of names) if (!sol.c[nm]) throw new RangeError(`particles: no species ${JSON.stringify(nm)} in the solution`);
  if (!(dots > 0 && Number.isFinite(dots))) throw new RangeError(`particles: dots must be a positive number, got ${dots}`);
  if (!(cells >= 1 && Number.isFinite(cells))) throw new RangeError(`particles: cells must be a number from 1 up, got ${cells}`);
  if (!(cap >= 0)) throw new RangeError(`particles: cap must be a number from 0 up, got ${cap}`);
  const MAX = 1e6; // the most dots a swarm will draw

  // The lattice: each region's nodes grouped into runs of node boxes about L/cells wide, so a
  // cell boundary is a segment, where the solution's flux is. (A metal's nodes, with no
  // concentration, are left out.) Drawn again for a solution on another grid.
  let x, nn, box, lattice, nc, at;
  function draw(sol) {
    x = Float64Array.from(sol.x);
    nn = x.length;
    box = Float64Array.from(x, (_, g) => ((g > 0 && sol.region[g - 1] === sol.region[g] ? x[g] - x[g - 1] : 0) + (g < nn - 1 && sol.region[g + 1] === sol.region[g] ? x[g + 1] - x[g] : 0)) / 2);
    const target = (x[nn - 1] - x[0]) / cells;
    lattice = [];
    sol.regions.forEach((_, r) => {
      const g0 = sol.region.indexOf(r), g1 = sol.region.lastIndexOf(r);
      if (names.every((nm) => !(sol.c[nm][g0] >= 0))) return;
      const start = lattice.length;
      let a = g0, w = 0;
      for (let g = g0; g <= g1; g++) {
        w += box[g];
        if (w >= target || g === g1) {
          lattice.push({ a, b: g, region: r, h: w });
          [a, w] = [g + 1, 0];
        }
      }
      const last = lattice.length - 1;
      if (last > start && lattice[last].h < target / 2) {
        // a thin last cell joins the one before it
        lattice[last - 1].b = lattice[last].b;
        lattice[last - 1].h += lattice[last].h;
        lattice.pop();
      }
    });
    for (const cl of lattice) {
      cl.x0 = x[cl.a] - (cl.a > 0 && sol.region[cl.a - 1] === cl.region ? (x[cl.a] - x[cl.a - 1]) / 2 : 0);
      cl.x1 = cl.x0 + cl.h;
      let mid = cl.a;
      for (let g = cl.a; g <= cl.b; g++) if (Math.abs(x[g] - (cl.x0 + cl.x1) / 2) < Math.abs(x[mid] - (cl.x0 + cl.x1) / 2)) mid = g;
      cl.mid = mid;
    }
    nc = lattice.length;
    // where each boundary is: k between cells k − 1 and k
    at = Float64Array.from({ length: nc + 1 }, (_, k) => (k < nc ? lattice[k].x0 : lattice[nc - 1].x1));
    swarm.cells = lattice.map(({ x0, x1, region }) => ({ x0, x1, region }));
    swarm.dots = [];
    for (const nm of names) swarm.crossed[nm] = { up: new Float64Array(nc + 1), down: new Float64Array(nc + 1) };
    for (const nm of names) delete S[nm];
  }
  const sameGrid = (s) => s.x.length === nn && s.x.every((v, g) => v === x[g]);

  const swarm = {
    cells: [],
    species: names,
    weight: {}, // mol/m² per dot
    sea: {}, // per species, per cell: 1 where it's a sea
    dots: [], // { species, cell, x }
    // what happened in the last step, beyond hops between cells of dots: { species, kind, x, dir }
    // with kind 'in' (a dot from a sea, contact or metal), 'out' (into one), 'made', 'unmade',
    // or 'cross' (from sea to sea through a face); dir ±1 for a crossing's direction
    events: [],
    crossed: {}, // per species, at each boundary (0 … cells): dots across it, { up, down } (toward +x, −x), counted since you last zeroed them
    time: 0,
    update,
    step,
    expected,
  };
  const S = {}; // per species: amounts, boundaries and rates, from the solution

  function amounts(sol, nm) {
    const c = sol.c[nm];
    return Float64Array.from(lattice, (cl) => {
      let n = 0;
      for (let g = cl.a; g <= cl.b; g++) if (c[g] > 0) n += c[g] * box[g];
      return n;
    });
  }

  function fields(sol, nm, prev) {
    const c = sol.c[nm], mu = sol.mu[nm], D = sol.D[nm], J = sol.flux[nm], RT = GAS_CONSTANT * sol.T;
    const N = amounts(sol, nm), w = swarm.weight[nm];
    const on = Uint8Array.from(lattice, (cl) => (c[cl.mid] >= 0 ? 1 : 0));
    const sea = Uint8Array.from(N, (n, K) => (on[K] && n / w > cap ? 1 : 0));
    const psi = Float64Array.from(lattice, (cl, K) => mu[cl.mid] / RT - Math.log(N[K] / cl.h));
    // Boundary k: what lies to each side, its one-way fluxes (mol/(m²·s)), and whether it's a face
    // with one-way fluxes of its own (a membrane), whose crossings between seas are events.
    const left = new Uint8Array(nc + 1), right = new Uint8Array(nc + 1), face = new Uint8Array(nc + 1);
    const up = new Float64Array(nc + 1), down = new Float64Array(nc + 1);
    const pair = (k, net, F0) => {
      const r = Math.hypot(F0, net / 2);
      up[k] = r + net / 2;
      down[k] = r - net / 2;
    };
    const side = (K) => (sea[K] ? OUTSIDE : DOT);
    for (let k = 0; k <= nc; k++) {
      const A = k > 0 && on[k - 1] ? lattice[k - 1] : null, B = k < nc && on[k] ? lattice[k] : null;
      if (A && B && A.b + 1 === B.a) {
        // two cells: within a region, or across a face
        if (A.region !== B.region) {
          const itf = sol.interfaces[A.region];
          if (!itf.links[nm]) continue;
          [left[k], right[k]] = [side(k - 1), side(k)];
          if (itf.oneWay?.[nm]) {
            face[k] = 1;
            [up[k], down[k]] = itf.oneWay[nm];
            continue;
          }
        }
        [left[k], right[k]] = [side(k - 1), side(k)];
        const hD = A.h / 2 / D[A.mid] + B.h / 2 / D[B.mid];
        pair(k, J[A.b], (Math.sqrt((N[k - 1] / A.h) * (N[k] / B.h)) * shape(psi[k] - psi[k - 1])) / hD);
      } else if (A || B) {
        // out of the cells: through a contact, or a face to a metal
        const cl = A ?? B, g = A ? A.b : B.a;
        let net, open;
        if (g === 0 || g === nn - 1) {
          const ct = sol.contacts[g === 0 ? 'left' : 'right'];
          open = ct.links[nm] !== undefined;
          net = ct.flux[nm]; // toward +x at either end
        } else {
          open = sol.interfaces[A ? A.region : B.region - 1].links[nm] !== undefined;
          net = J[A ? g : g - 1];
        }
        if (!open) continue;
        [left[k], right[k]] = A ? [side(k - 1), OUTSIDE] : [OUTSIDE, side(k)];
        pair(k, net, (2 * D[g] * Math.max(c[g], 0)) / cl.h);
      }
    }
    // Made and unmade in each cell (mol/(m²·s)): by the bulk reactions' one-way rates (so in
    // equilibrium, generation and recombination both go on, in balance), and whatever else the
    // fluxes don't account for, ∂N/∂t − (J_in − J_out) less the reactions' net (ports, face
    // reactions), with ∂N/∂t from the last solution if time has passed since.
    const made = new Float64Array(nc), unmade = new Float64Array(nc);
    for (const rx of sol.bulkReactions ?? []) {
      const nu = rx.nu?.[nm];
      if (!nu) continue;
      lattice.forEach((cl, K) => {
        for (let g = cl.a; g <= cl.b; g++) {
          if (!(rx.forward[g] >= 0)) continue;
          const fwd = rx.forward[g] * box[g], bwd = (rx.forward[g] - rx.rate[g]) * box[g];
          made[K] += nu > 0 ? nu * fwd : -nu * bwd;
          unmade[K] += nu > 0 ? nu * bwd : -nu * fwd;
        }
      });
    }
    const dt = prev && sol.time > prev.time ? sol.time - prev.time : Infinity;
    for (let K = 0; K < nc; K++) {
      if (!on[K]) continue;
      const rest = (Number.isFinite(dt) ? (N[K] - prev.N[K]) / dt : 0) + up[K + 1] - down[K + 1] - (up[K] - down[K]) - (made[K] - unmade[K]);
      if (rest > 0) made[K] += rest;
      else unmade[K] -= rest;
    }
    // Where in a cell a dot sits: in a node's box, chosen by its share of the cell's amount.
    const share = Float64Array.from(x, (_, g) => (c[g] > 0 ? c[g] * box[g] : 0));
    return { N, on, sea, left, right, face, up, down, made, unmade, share, time: sol.time };
  }

  /**
   * A new solution: the dots keep their places and move on at its rates. A cell that has become a
   * sea absorbs its dots, and one that's no longer a sea is filled. With resample (a jump between
   * steady states, which has no path in time to follow), each cell's dots are thinned or added to
   * by the change in what it's expected to hold, m → m′: each kept with probability m′/m if it
   * falls, Poisson(m′ − m) more if it rises, which takes a Poisson sample of the one into a
   * Poisson sample of the other, moving as few dots as can be. A solution on another grid (a
   * resized device) redraws the lattice, and the dots are placed afresh.
   */
  function update(next, { resample = false } = {}) {
    if (!next?.x || !next.c) throw new TypeError('particles: update takes a solution');
    for (const nm of names) if (!next.c[nm]) throw new RangeError(`particles: the new solution has no species ${JSON.stringify(nm)}`);
    sol = next;
    if (!sameGrid(sol)) draw(sol); // (a resized device: the dots placed afresh)
    let total = 0;
    for (const nm of names) total += amounts(sol, nm).reduce((a, n, K) => a + (n / swarm.weight[nm] > cap ? 0 : n), 0) / swarm.weight[nm];
    if (total > MAX) throw new RangeError(`particles: this would draw ${total.toPrecision(2)} dots; give a larger weight, or a cap (dense cells as seas)`);
    for (const nm of names) {
      const was = S[nm];
      S[nm] = fields(sol, nm, was);
      const F = S[nm], w = swarm.weight[nm];
      swarm.sea[nm] = F.sea;
      const before = (K) => (!was || was.sea[K] || !was.on[K] ? 0 : was.N[K] / w);
      swarm.dots = swarm.dots.filter((d) => {
        if (d.species !== nm) return true;
        if (!F.on[d.cell] || F.sea[d.cell]) return false;
        const m = before(d.cell);
        return !(resample && m > 0) || random() * m < F.N[d.cell] / w;
      });
      F.N.forEach((n, K) => {
        if (!F.on[K] || F.sea[K]) return;
        const m = before(K);
        if (m === 0 || resample) for (let j = poisson(Math.max(0, n / w - m)); j > 0; j--) swarm.dots.push(dot(nm, K));
      });
    }
    return swarm;
  }

  /** The dots expected in each cell for a species, at the solution's concentrations. */
  function expected(nm) {
    if (!S[nm]) throw new RangeError(`particles: no species ${JSON.stringify(nm)} in this swarm (${names.join(', ')})`);
    return Float64Array.from(S[nm].N, (n) => n / swarm.weight[nm]);
  }

  // A place in cell K for a dot of species nm: a node's box, by its share of the amount, then
  // anywhere in that box. (A dot follows the profile within its cell, as at a depletion edge.)
  function place(nm, K) {
    const cl = lattice[K], share = S[nm]?.share;
    let g = cl.a;
    if (share && S[nm].N[K] > 0) {
      let u = random() * S[nm].N[K];
      while (g < cl.b && (u -= share[g]) > 0) g++;
    } else return cl.x0 + random() * cl.h;
    const lo = g > cl.a ? (x[g - 1] + x[g]) / 2 : cl.x0, hi = g < cl.b ? (x[g] + x[g + 1]) / 2 : cl.x1;
    return lo + random() * (hi - lo);
  }
  const dot = (nm, K) => ({ species: nm, cell: K, x: place(nm, K) });

  // A Poisson count of mean m (a normal one for large m).
  const gauss = () => Math.sqrt(-2 * Math.log(1 - random())) * Math.cos(2 * Math.PI * random());
  function poisson(m) {
    if (m > 50) return Math.max(0, Math.round(m + Math.sqrt(m) * gauss()));
    const e = Math.exp(-m);
    let k = 0;
    for (let p = random(); p > e; p *= random()) k++;
    return k;
  }

  /**
   * Move the dots on by dt (s, the device's time): each hops, and dots are made and unmade, at
   * the solution's rates. A cell a dot would hop out of more than maxHops/4 times in the step is
   * a sea for the step (whatever crosses its sides counted at the one-way fluxes, dots walking in
   * absorbed), then drawn afresh, a new sample of what it holds: so a step much longer than the
   * hop time costs little, and stays a sample of the solution. A dt that's NaN or negative (a
   * first frame, a clock that stepped back) moves nothing.
   */
  function step(dt, { maxHops = 200 } = {}) {
    // (A first frame or a clock stepping back, NaN or negative, moves nothing.)
    swarm.events = [];
    if (dt === Infinity) throw new RangeError('particles: step(dt) takes a finite time (s), got Infinity');
    if (!(dt > 0)) return swarm;
    const fast = {}, saved = {};
    for (const nm of names) {
      const F = S[nm];
      fast[nm] = Uint8Array.from(F.N, (N, K) => {
        if (!F.on[K] || F.sea[K] || !(N > 0)) return 0;
        const rate = ((F.right[K + 1] ? F.up[K + 1] : 0) + (F.left[K] ? F.down[K] : 0) + F.unmade[K]) / N;
        return rate * dt > maxHops / 4 ? 1 : 0;
      });
      if (!fast[nm].some(Boolean)) continue;
      // For this step, those cells are seas.
      saved[nm] = [F.sea, F.left, F.right];
      F.sea = Uint8Array.from(F.sea, (v, K) => v | fast[nm][K]);
      F.left = Uint8Array.from(F.left, (v, k) => (v && k > 0 && fast[nm][k - 1] ? OUTSIDE : v));
      F.right = Uint8Array.from(F.right, (v, k) => (v && k < nc && fast[nm][k] ? OUTSIDE : v));
    }
    // Births first, each at a random moment of the step: made in a cell, or arriving from beyond
    // (a sea, a contact, a metal). Each then moves for the rest of the step. And crossings
    // between seas: through a face, as events; between cells that are seas only for this step,
    // counted.
    const born = [];
    for (const nm of names) {
      const F = S[nm], w = swarm.weight[nm], cr = swarm.crossed[nm];
      for (let K = 0; K < nc; K++) {
        if (!F.on[K] || F.sea[K]) continue;
        const made = F.made[K] / w;
        const fromLeft = F.left[K] === OUTSIDE ? F.up[K] / w : 0, fromRight = F.right[K + 1] === OUTSIDE ? F.down[K + 1] / w : 0;
        const rate = made + fromLeft + fromRight;
        if (!(rate > 0)) continue;
        for (let j = poisson(rate * dt); j > 0; j--) {
          const u = random() * rate, d = dot(nm, K);
          // (no 'in' from a cell that's a sea only for this step)
          if (u < fromLeft) {
            cr.up[K]++;
            if (!fast[nm][K - 1]) swarm.events.push({ species: nm, kind: 'in', x: at[K], dir: 1 });
          } else if (u < fromLeft + fromRight) {
            cr.down[K + 1]++;
            if (!fast[nm][K + 1]) swarm.events.push({ species: nm, kind: 'in', x: at[K + 1], dir: -1 });
          } else swarm.events.push({ species: nm, kind: 'made', x: d.x });
          born.push([d, random() * dt]);
        }
      }
      for (let k = 0; k <= nc; k++) {
        if (F.left[k] !== OUTSIDE || F.right[k] !== OUTSIDE) continue;
        const events = F.face[k] && !(fast[nm][k - 1] || fast[nm][k]);
        if (!events && !(fast[nm][k - 1] || fast[nm][k])) continue;
        for (const [dir, flux, count] of [[1, F.up[k], cr.up], [-1, F.down[k], cr.down]]) {
          const n = poisson((flux / w) * dt);
          count[k] += n;
          if (events) for (let j = 0; j < n; j++) swarm.events.push({ species: nm, kind: 'cross', x: at[k], dir });
        }
      }
    }
    const keep = [];
    for (const d of swarm.dots) if (!fast[d.species][d.cell] && walk(d, dt, maxHops, fast[d.species])) keep.push(d);
    for (const [d, t] of born) if (walk(d, dt - t, maxHops, fast[d.species])) keep.push(d);
    swarm.dots = keep;
    // The step's seas drawn afresh, and back to cells of dots.
    for (const nm of Object.keys(saved)) {
      const F = S[nm], w = swarm.weight[nm];
      [F.sea, F.left, F.right] = saved[nm];
      fast[nm].forEach((f, K) => {
        if (f) for (let j = poisson(F.N[K] / w); j > 0; j--) swarm.dots.push(dot(nm, K));
      });
    }
    swarm.time += dt;
    return swarm;
  }

  // A dot's hops for a time dt (Gillespie's: exact for rates held fixed): whether it's still in.
  function walk(d, dt, maxHops, fast) {
    const F = S[d.species], cr = swarm.crossed[d.species];
    let t = 0;
    for (let hop = 0; hop < maxHops; hop++) {
      const K = d.cell, N = F.N[K];
      if (!(N > 0)) return true;
      const r = F.right[K + 1] ? F.up[K + 1] / N : 0, l = F.left[K] ? F.down[K] / N : 0, die = F.unmade[K] / N;
      const total = r + l + die;
      if (!(total > 0)) return true;
      t -= Math.log(1 - random()) / total;
      if (t > dt) return true;
      const u = random() * total;
      if (u >= r + l) {
        swarm.events.push({ species: d.species, kind: 'unmade', x: d.x });
        return false;
      }
      const k = u < r ? K + 1 : K, dir = u < r ? 1 : -1;
      (dir > 0 ? cr.up : cr.down)[k]++;
      if ((dir > 0 ? F.right[k] : F.left[k]) === OUTSIDE) {
        if (!fast?.[K + dir]) swarm.events.push({ species: d.species, kind: 'out', x: at[k], dir });
        return false;
      }
      d.cell = K + dir;
      d.x = place(d.species, d.cell);
    }
    return true;
  }

  // The lattice and weights, then the first dots, placed by the concentration.
  draw(sol);
  for (const nm of names) {
    const total = amounts(sol, nm).reduce((a, b) => a + b, 0);
    const w = typeof weight === 'number' ? weight : isObject(weight) && weight[nm] !== undefined ? weight[nm] : total / dots;
    if (!(w > 0 && Number.isFinite(w))) throw new RangeError(`particles: the weight of ${nm} (mol/m² per dot) must be a positive number, got ${w}` + (total > 0 ? '' : ' (it has nothing to draw)'));
    swarm.weight[nm] = w;
  }
  return update(sol);
}
