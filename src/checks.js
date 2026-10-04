// Checks on a solved device, for driftlet/kit: what any result should satisfy, worked out from the
// definition and the solution, so a page can show that its numbers hold up without its author
// writing each check by hand. Each check says what it compared and whether it passed.

import { Device } from './index.js';
import { FARADAY } from './constants.js';
import { unitWarnings } from './describe.js';

const sig = (v, p = 3) => (v === 0 || !Number.isFinite(v) ? String(v) : Math.abs(v) >= 1e4 || Math.abs(v) < 1e-2 ? v.toExponential(p - 1) : String(+v.toPrecision(p)));

/**
 * @typedef {object} CheckItem
 * @property {string} name what was checked: 'converged', 'warnings', 'balance', 'conservation' or 'grid'
 * @property {boolean | null} ok whether it passed (null: look at it, a warning or a check that
 *   couldn't be made)
 * @property {string} summary one line
 * @property {object} [details] the numbers behind it
 */

/**
 * Check a device's present state (pass the solution solve() or advance() returned, for its
 * convergence and warnings):
 * - **converged**: the last solve converged.
 * - **warnings**: the solution's warnings (unresolved double layers, steep profiles, weakly held
 *   parts) and the definition's likely unit slips.
 * - **balance** (a steady state): each species' ledger, what comes in through each terminal and
 *   what each reaction (bulk or at a face) makes or uses, in mol/(m²·s) and, for a charged
 *   species, as a current F·|z|·rate. The terms must sum to zero; under light it's J = F(G − R),
 *   and it says where every carrier went.
 * - **conservation** (a transient): every closed, unreacting stretch keeps what it held plus
 *   what came in.
 * - **grid** (a steady state, with `refine`, the default): the device solved again on a grid
 *   twice as fine everywhere, from this solution interpolated onto it. The current, each floating terminal's voltage and each
 *   region's charge must agree to `tol` (relative, default 1e-2, a plot's accuracy; a benchmark
 *   wants less). The discretisation is second order, so the change estimates this grid's error
 *   (about 3/4 of it).
 * @param {import('./index.js').Device} device
 * @param {object} [sol] the solution to check (default: the device's present state)
 * @param {{ refine?: boolean, tol?: number }} [opts]
 * Every check passed when `ok`; an item with `ok: null` (a warning, or a check that couldn't be
 * made) asks you to look, and doesn't fail it. `text` is a line per check, marked ok, FAIL or ?.
 * @returns {{ ok: boolean, items: CheckItem[], text: string }}
 */
export function check(device, sol, { refine = true, tol = 1e-2 } = {}) {
  sol ??= device.solution();
  const { model, solver, def } = device;
  const steady = sol.steady ?? solver.atSteady;
  const items = [];

  items.push(
    sol.converged
      ? { name: 'converged', ok: true, summary: `${steady ? 'steady state' : `transient at t = ${sig(sol.time)} s`}, converged` }
      : { name: 'converged', ok: false, summary: "the solve didn't converge: these numbers aren't a solution" },
  );

  const warnings = [...sol.warnings, ...unitWarnings(def)];
  // Heuristics, so a prompt to look (null) rather than a failure.
  items.push({
    name: 'warnings',
    ok: warnings.length === 0 ? true : null,
    summary: warnings.length === 0 ? 'none' : `${warnings.length}: ${warnings.join(' · ')}`,
    details: { warnings },
  });

  // The rest checks a solution, which an unconverged state isn't.
  if (sol.converged) {
    const ledger = steady ? balance(device, sol) : null;
    items.push(ledger ?? conservation(sol, Math.max(...device.model.grid.area)));
    if (steady && refine) items.push(gridCheck(device, sol, tol, ledger.details.gross));
  }

  const ok = sol.converged && items.every((it) => it.ok !== false);
  const mark = (it) => (it.ok === true ? 'ok  ' : it.ok === false ? 'FAIL' : '?   ');
  const text = items.map((it) => `${mark(it)} ${it.name}: ${it.summary}${it.details?.lines?.length ? `\n${it.details.lines.map((l) => `       ${l}`).join('\n')}` : ''}`).join('\n');
  return { ok, items, text };
}

// Each species' steady ledger: in through the terminals, made by the reactions. Rates in
// mol/(m²·s), positive into the device or made.
function balance(device, sol) {
  const { model, def } = device, { species, regions, interfaces } = model;
  const label = (rdef) => rdef?.equation ?? JSON.stringify(rdef?.nu ?? rdef);
  // Each species' compartments: runs of regions it's present in, joined by faces it crosses
  // freely (in equilibrium). A face it crosses only by a law (a permeability, a conductance) or
  // by a reaction divides them, and what crosses there is a term in each side's ledger.
  const ledgers = [], home = species.map(() => new Int32Array(regions.length).fill(-1));
  species.forEach((sp, i) => {
    let open = null;
    regions.forEach((reg, r) => {
      if (!model.materials[reg.material].present[i]) return void (open = null);
      if (open && interfaces[r - 1].links[i].type === 'equilibrium') open.regions.push(r);
      else ledgers.push((open = { species: sp.name, z: sp.z, i, regions: [r], terms: [] }));
      home[i][r] = ledgers.length - 1;
    });
  });
  const add = (i, r, what, rate) => {
    const l = home[i][r];
    if (!(l >= 0 && rate !== 0 && Number.isFinite(rate))) return;
    const same = ledgers[l].terms.find((t) => t.what === what); // a bulk reaction, region by region
    if (same) same.rate += rate;
    else ledgers[l].terms.push({ what, rate });
  };
  const last = regions.length - 1;
  species.forEach((sp, i) => {
    add(i, 0, 'left contact', sol.contacts.left.flux[sp.name]);
    add(i, last, 'right contact', -sol.contacts.right.flux[sp.name]);
    model.ports.forEach((port, k) => add(i, port.region, `port ${sol.ports[k].name}`, sol.ports[k].flux[sp.name]));
  });
  model.reactions.forEach((rx, k) => {
    const what = label(def.bulkReactions?.[k]);
    sol.bulkReactions[k].regions.forEach((total, r) => {
      for (const { i, nu } of rx.reactants) add(i, r, what, -nu * total);
      for (const { i, nu } of rx.products) add(i, r, what, nu * total);
    });
  });
  interfaces.forEach((itf, f) => {
    // Across the face, where it divides a species' compartments.
    species.forEach((sp, i) => {
      const type = itf.links[i].type;
      if (type === 'equilibrium' || type === 'blocked') return;
      const N = sol.interfaces[f].N[sp.name];
      add(i, f, `across face ${f} (${type})`, -N);
      add(i, f + 1, `across face ${f} (${type})`, N);
    });
    itf.reactions.forEach((rx, k) => {
      const rate = sol.interfaces[f].rates[k], what = `${label(def.interfaces?.[f]?.reactions?.[k])} at face ${f}`;
      for (const { i, nu, side } of rx.part) add(i, f + side, what, nu * rate);
    });
  });
  // A species in one compartment is named alone; one in several, by where.
  const several = species.map((_, i) => ledgers.filter((l) => l.i === i).length > 1);
  for (const l of ledgers) {
    const names = l.regions.map((r) => sol.regions[r].name);
    l.compartment = names.length === 1 ? names[0] : `${names[0]} to ${names.at(-1)}`;
    l.title = several[l.i] ? `${l.species} in ${l.compartment}` : l.species;
  }

  // A species at rest (nothing beyond round-off moving, against what it could carry) isn't reported.
  const size = (l) => Math.max(0, ...l.terms.map((t) => Math.abs(t.rate)));
  const natural = naturalScales(device, sol);
  const moving = ledgers.filter((l) => size(l) > 1e-9 * natural[l.i]);
  let worst = 0;
  const lines = [];
  for (const l of moving) {
    const sum = l.terms.reduce((a, t) => a + t.rate, 0);
    l.residual = sum / size(l);
    worst = Math.max(worst, Math.abs(l.residual));
    // Terms that matter, largest first; a charged species' as a current too.
    const shown = l.terms.filter((t) => Math.abs(t.rate) > 1e-6 * size(l)).sort((a, b) => Math.abs(b.rate) - Math.abs(a.rate));
    const amp = (r) => (l.z !== 0 ? ` (${sig(FARADAY * Math.abs(l.z) * r)} A/m²)` : '');
    lines.push(`${l.title}: ${shown.map((t) => `${t.what}: ${t.rate > 0 ? '+' : ''}${sig(t.rate)}${amp(t.rate)}`).join('; ')}`);
  }
  const ok = worst < 1e-6;
  const summary =
    moving.length === 0
      ? 'nothing flows or reacts (equilibrium)'
      : `${moving.length === 1 ? 'one ledger sums' : `${moving.length} ledgers sum`} to zero ${ok ? `(to ${sig(worst, 1)})` : `only to ${sig(worst, 2)}: not a steady state`}: each species' sources and sinks, mol/(m²·s), + in or made`;
  // The largest current any charged species carries anywhere: the scale a net current is read against.
  const gross = Math.max(0, ...moving.filter((l) => l.z !== 0).map((l) => FARADAY * Math.abs(l.z) * size(l)));
  return { name: 'balance', ok, summary, details: { ledgers, worst, lines, gross } };
}

// What each species could carry, mol/(m²·s): by transport, D·c/L in each region (a metal's
// carrier, σRT/(z²F²L)), and by its bulk reactions, their one-way rates. In equilibrium the net
// flows are round-off of these.
function naturalScales(device, sol) {
  const { model, solver } = device, { species, grid } = model;
  // (totals through the cross-section, as the ledgers are: its largest in each region)
  const area = model.regions.map((_, r) => Math.max(...grid.area.subarray(grid.regionStart[r], grid.regionEnd[r] + 1)));
  return species.map((sp, i) => {
    let S = 0;
    model.regions.forEach((reg, r) => {
      const mat = model.materials[reg.material], L = reg.length;
      if (mat.conductor) {
        if (mat.conductor.i === i) S = Math.max(S, (area[r] * mat.conductor.sigma * model.RT) / (sp.z * sp.z * FARADAY * FARADAY * L));
        return;
      }
      if (!mat.present[i] || !(mat.D[i] > 0)) return;
      let cmax = 0;
      for (let g = 0; g < sol.x.length; g++) if (sol.region[g] === r && Number.isFinite(sol.c[sp.name][g])) cmax = Math.max(cmax, sol.c[sp.name][g]);
      S = Math.max(S, (area[r] * mat.D[i] * cmax) / L);
    });
    solver.rxs.forEach((rx) => {
      const p = rx.sp.indexOf(i);
      if (p < 0) return;
      let gross = 0;
      for (let g = 0; g < sol.x.length; g++) if (rx.kf[solver.nodeMaterial[g]] > 0) gross += grid.vol[g] * solver.bulkOneWay(rx, g);
      S = Math.max(S, Math.abs(rx.nu[p]) * gross);
    });
    return S;
  });
}

// Closed, unreacting stretches keep their amounts in a transient.
function conservation(sol, area = 1) {
  const closed = sol.conservation.filter((st) => !st.connected && !st.reactive);
  if (closed.length === 0) return { name: 'conservation', ok: true, summary: 'nothing closed to conserve (every species reaches a terminal or reacts)' };
  // Relative to the amount, but a drift below 1e-20 mol/m² (1e-15 C/m²) is nothing: a floating
  // metal's excess electrons are themselves ~0. (Amounts are totals through the cross-section.)
  const worst = Math.max(...closed.map((st) => (Math.abs(st.amount - st.reference - st.intake) < 1e-20 * area ? 0 : Math.abs(st.drift))));
  return {
    name: 'conservation',
    ok: worst < 1e-8,
    summary: `${closed.length} closed stretch${closed.length > 1 ? 'es' : ''} (a species held between blocking faces) kept ${worst < 1e-8 ? `to ${sig(worst, 1)}` : `only to ${sig(worst, 2)}`} of what ${closed.length > 1 ? 'they' : 'it'} held plus what came in`,
    details: { worst },
  };
}

// The same device on a grid twice as fine (half the cells' sizes, the square root of the
// grading ratio), solved from this solution interpolated onto it, against this solution.
function gridCheck(device, sol, tol, gross) {
  const { def, model } = device;
  const fine = {
    ...def,
    regions: def.regions.map((reg, r) => {
      if (model.materials[model.regions[r].material].conductor) return reg;
      const g = { ...(def.grid ?? {}), ...(reg.grid ?? {}) };
      if (def.grid === undefined && reg.grid === undefined) Object.assign(g, { hmin: reg.length / 1000, hmax: reg.length / 20 });
      const out = { ratio: Math.sqrt(g.ratio ?? 1.2), minCells: 2 * (g.minCells ?? 8) };
      if (g.hmin !== undefined) out.hmin = g.hmin / 2;
      if (g.hmax !== undefined) out.hmax = g.hmax / 2;
      return { ...reg, grid: out };
    }),
  };
  // Started from this solution, interpolated onto the finer grid; failing that, from cold.
  let other;
  try {
    const dev = new Device(fine);
    dev.solver._warmFrom(device.solver);
    other = dev.solve();
    if (!other.converged) {
      const cold = new Device(fine);
      if (sol.time > 0) cold.solver.time = sol.time; // waveforms at the same time
      other = cold.solve();
    }
  } catch (e) {
    return { name: 'grid', ok: null, summary: `couldn't solve the refined device (${e.message})` };
  }
  if (!other.converged) return { name: 'grid', ok: null, summary: "the device on a twice-finer grid didn't solve; check resolution by hand" };

  const diffs = [];
  // The current, against the largest current any species carries (so at open circuit, or in
  // equilibrium, round-off isn't read as a change).
  const scaleI = Math.max(Math.abs(sol.current), Math.abs(other.current), gross);
  if (gross > 0) diffs.push({ what: 'current', a: sol.current, b: other.current, rel: Math.abs(other.current - sol.current) / scaleI, unit: 'A/m²' });
  // Floating terminals (driven by a current): their voltages, against the larger and the thermal voltage.
  const VT = model.RT / FARADAY;
  model.terminals.forEach((t, k) => {
    if (device.solver.floating.includes(k)) {
      const a = sol.terminals[t.name].V, b = other.terminals[t.name].V;
      diffs.push({ what: `${t.name} voltage`, a, b, rel: Math.abs(b - a) / Math.max(Math.abs(a), Math.abs(b), VT), unit: 'V' });
    }
  });
  // Each region's charge (C/m²), against the largest: summed over the nodes' boxes, as the
  // solver's Gauss law counts it.
  const charges = (s) => s.regions.map((_, r) => {
    let q = 0;
    for (let g = 0; g < s.x.length; g++) {
      if (s.region[g] !== r) continue;
      const left = g > 0 && s.region[g - 1] === r ? s.x[g] - s.x[g - 1] : 0, right = g + 1 < s.x.length && s.region[g + 1] === r ? s.x[g + 1] - s.x[g] : 0;
      q += rho(s, g, model, r) * ((left + right) / 2);
    }
    return q;
  });
  // Charges at round-off (strictly neutral regions) aren't compared: against the gross charge
  // the regions hold, |fixed| + Σ|z|c, they must be more than a millionth.
  const content = Math.max(...sol.regions.map((reg, r) => {
    let q = 0;
    for (let g = 0; g < sol.x.length; g++) {
      if (sol.region[g] !== r) continue;
      let a = Math.abs(model.regions[r].fixedCharge ?? 0);
      for (const sp of model.species) if (Number.isFinite(sol.c[sp.name][g])) a += FARADAY * Math.abs(sp.z * sol.c[sp.name][g]);
      q = Math.max(q, a);
    }
    return q * (reg.x1 - reg.x0);
  }));
  const qa = charges(sol), qb = charges(other), scaleQ = Math.max(0, ...qa.map(Math.abs), ...qb.map(Math.abs));
  if (scaleQ > 1e-6 * content) {
    const named = (reg) => (/^region \d+$/.test(reg.name) ? `${reg.name} (${reg.material})` : reg.name);
    qa.forEach((a, r) => diffs.push({ what: `${named(sol.regions[r])} charge`, a, b: qb[r], rel: Math.abs(qb[r] - a) / scaleQ, unit: 'C/m²' }));
  }
  const worst = diffs.reduce((w, d) => (d.rel > (w?.rel ?? -1) ? d : w), null);
  if (!worst) return { name: 'grid', ok: true, summary: 'nothing to compare (no current, charge or floating voltage)' };
  const ok = worst.rel <= tol;
  return {
    name: 'grid',
    ok,
    summary:
      `on a grid twice as fine (${sol.x.length} → ${other.x.length} nodes), the largest change is ${worst.what}, ${sig(worst.a, 5)} → ${sig(worst.b, 5)} ${worst.unit} ` +
      `(${sig(worst.rel, 2)}${ok ? '' : `, over ${sig(tol, 1)}: refine the grid (hmin, hmax), most likely where profiles are steepest`})`,
    details: { diffs, nodes: [sol.x.length, other.x.length] },
  };
}

// Charge density at a node of region r, C/m³: the species plus the region's fixed charge.
function rho(s, g, model, r) {
  let q = model.regions[r].fixedCharge ?? 0;
  for (const sp of model.species) {
    const c = s.c[sp.name][g];
    if (Number.isFinite(c)) q += FARADAY * sp.z * c;
  }
  return q;
}
