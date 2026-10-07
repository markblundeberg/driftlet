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
  const nn = sol.x.length, x = sol.x;

  // The lattice: each region's nodes grouped into runs of node boxes about L/cells wide, so a
  // cell boundary is a segment, where the solution's flux is. (A metal's nodes, with no
  // concentration, are left out.)
  const box = Float64Array.from(x, (_, g) => ((g > 0 && sol.region[g - 1] === sol.region[g] ? x[g] - x[g - 1] : 0) + (g < nn - 1 && sol.region[g + 1] === sol.region[g] ? x[g + 1] - x[g] : 0)) / 2);
  const target = (x[nn - 1] - x[0]) / cells;
  const lattice = [];
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
  const nc = lattice.length;
  // where each boundary is: k between cells k − 1 and k
  const at = Float64Array.from({ length: nc + 1 }, (_, k) => (k < nc ? lattice[k].x0 : lattice[nc - 1].x1));

  const swarm = {
    cells: lattice.map(({ x0, x1, region }) => ({ x0, x1, region })),
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
    // Made (or, negative, unmade) in each cell: what the fluxes don't account for,
    // ∂N/∂t − (J_in − J_out), with ∂N/∂t from the last solution if time has passed since.
    const dt = prev && sol.time > prev.time ? sol.time - prev.time : Infinity;
    const make = Float64Array.from(N, (n, K) => (on[K] ? (Number.isFinite(dt) ? (n - prev.N[K]) / dt : 0) + up[K + 1] - down[K + 1] - (up[K] - down[K]) : 0));
    return { N, on, sea, left, right, face, up, down, make, time: sol.time };
  }

  /**
   * A new solution: the dots keep their places and move on at its rates. A cell that has become a
   * sea absorbs its dots, and one that's no longer a sea is filled. With resample (a jump between
   * steady states, which has no path in time to follow), each cell's dots are thinned or added to
   * by the change in what it's expected to hold, m → m′: each kept with probability m′/m if it
   * falls, Poisson(m′ − m) more if it rises, which takes a Poisson sample of the one into a
   * Poisson sample of the other, moving as few dots as can be.
   */
  function update(next, { resample = false } = {}) {
    sol = next;
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
    return Float64Array.from(S[nm].N, (n) => n / swarm.weight[nm]);
  }

  const dot = (nm, K) => ({ species: nm, cell: K, x: lattice[K].x0 + random() * lattice[K].h });

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
   * the solution's rates. A dot hops at most maxHops times in a step.
   */
  function step(dt, { maxHops = 200 } = {}) {
    swarm.events = [];
    // Births first, each at a random moment of the step: made in a cell, or arriving from beyond
    // (a sea, a contact, a metal). Each then moves for the rest of the step. And crossings
    // between seas through a face.
    const born = [];
    for (const nm of names) {
      const F = S[nm], w = swarm.weight[nm], cr = swarm.crossed[nm];
      for (let K = 0; K < nc; K++) {
        if (!F.on[K] || F.sea[K]) continue;
        const made = F.make[K] > 0 ? F.make[K] / w : 0;
        const fromLeft = F.left[K] === OUTSIDE ? F.up[K] / w : 0, fromRight = F.right[K + 1] === OUTSIDE ? F.down[K + 1] / w : 0;
        const rate = made + fromLeft + fromRight;
        if (!(rate > 0)) continue;
        for (let j = poisson(rate * dt); j > 0; j--) {
          const u = random() * rate, d = dot(nm, K);
          if (u < fromLeft) {
            cr.up[K]++;
            swarm.events.push({ species: nm, kind: 'in', x: at[K], dir: 1 });
          } else if (u < fromLeft + fromRight) {
            cr.down[K + 1]++;
            swarm.events.push({ species: nm, kind: 'in', x: at[K + 1], dir: -1 });
          } else swarm.events.push({ species: nm, kind: 'made', x: d.x });
          born.push([d, random() * dt]);
        }
      }
      for (let k = 0; k <= nc; k++) {
        if (!F.face[k] || F.left[k] !== OUTSIDE || F.right[k] !== OUTSIDE) continue;
        for (const [dir, flux, count] of [[1, F.up[k], cr.up], [-1, F.down[k], cr.down]]) {
          for (let j = poisson((flux / w) * dt); j > 0; j--) {
            count[k]++;
            swarm.events.push({ species: nm, kind: 'cross', x: at[k], dir });
          }
        }
      }
    }
    const keep = [];
    for (const d of swarm.dots) if (walk(d, dt, maxHops)) keep.push(d);
    for (const [d, t] of born) if (walk(d, dt - t, maxHops)) keep.push(d);
    swarm.dots = keep;
    swarm.time += dt;
    return swarm;
  }

  // A dot's hops for a time dt (Gillespie's: exact for rates held fixed): whether it's still in.
  function walk(d, dt, maxHops) {
    const F = S[d.species], cr = swarm.crossed[d.species];
    let t = 0;
    for (let hop = 0; hop < maxHops; hop++) {
      const K = d.cell, N = F.N[K];
      if (!(N > 0)) return true;
      const r = F.right[K + 1] ? F.up[K + 1] / N : 0, l = F.left[K] ? F.down[K] / N : 0, die = F.make[K] < 0 ? -F.make[K] / N : 0;
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
        swarm.events.push({ species: d.species, kind: 'out', x: at[k], dir });
        return false;
      }
      d.cell = K + dir;
      d.x = lattice[d.cell].x0 + random() * lattice[d.cell].h;
    }
    return true;
  }

  // Weights, then the first dots, placed by the concentration.
  for (const nm of names) {
    const total = amounts(sol, nm).reduce((a, b) => a + b, 0);
    swarm.weight[nm] = typeof weight === 'number' ? weight : isObject(weight) && weight[nm] > 0 ? weight[nm] : total / dots;
    swarm.crossed[nm] = { up: new Float64Array(nc + 1), down: new Float64Array(nc + 1) };
  }
  return update(sol);
}
