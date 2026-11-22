// Checks on a solved device, for driftlet/kit: what any result should satisfy, worked out from the
// definition and the solution, so a page can show that its numbers hold up without its author
// writing each check by hand. Each check says what it compared and whether it passed.

import { Device } from './index.js';
import { FARADAY } from './constants.js';
import { fixedAt, sourceAt } from './device.js';
import { DEFAULTS, defaultSpacing, gridGiven } from './grid.js';
import { unitWarnings } from './describe.js';
import { carryState } from './carry.js';
import { checkOptions } from './errors.js';

const sig = (v, p = 3) => (v === 0 || !Number.isFinite(v) ? String(v) : Math.abs(v) >= 1e4 || Math.abs(v) < 1e-2 ? v.toExponential(p - 1) : String(+v.toPrecision(p)));

/**
 * @typedef {object} CheckItem
 * @property {string} name what was checked: 'solution' (that it's the device's state), 'converged',
 *   'warnings', 'balance', 'conservation', 'grid' or 'order' (the species listed in reverse)
 * @property {boolean | null} ok whether it passed (null: look at it, a warning or a check that
 *   couldn't be made)
 * @property {string} summary one line
 * @property {object} [details] the numbers behind it
 */

/**
 * Check a device's present state (pass the solution solve() or advance() returned, for its
 * convergence and warnings; one the device has since moved on from, by a set(), a solve or a step,
 * fails as 'solution', by the mark each solution carries, or for a copy by its numbers):
 * - **converged**: the last solve converged.
 * - **warnings**: the solution's warnings (unresolved double layers, steep profiles, weakly held
 *   parts) and the definition's likely unit slips.
 * - **balance** (a steady state): each species' ledger, what comes in through each terminal and
 *   what each reaction (bulk or at a face) makes or uses, in mol/(m²·s) (mol/s if the device
 *   isn't planar) and, for a charged species, as a current F·|z|·rate. The terms must sum to
 *   zero; under light it's J = F(G − R), and it says where every carrier went. An electrode's
 *   surface species has a ledger too, and a term that isn't a number fails it.
 * - **conservation** (a transient): every closed, unreacting stretch keeps what it held plus
 *   what came in.
 * - **order** (a steady state, with `refine` and `order`, the defaults): the device solved from
 *   its definition twice, its species listed as given and in reverse, which must give the same
 *   steady state (each terminal's voltage and current, each species' amount). A solve that
 *   depends on the listing has converged somewhere that isn't steady without noticing (round 11's
 *   fuzzer found it in ~1% of random reacting devices, by exactly this test).
 * - **grid** (with `refine`, the default; in a transient, `ok: null` and a reminder, since the
 *   answer depends on the whole run): in a steady state, the device solved again on a grid
 *   twice as fine everywhere, from this solution interpolated onto it. Every terminal's current,
 *   each neutral species' flux at the contacts, each floating terminal's voltage (against the span
 *   of the terminals') and each region's charge (with its faces' diffuse layers, and each layer's)
 *   must agree to `tol` (with nothing to compare, `ok: null`) (relative, default 1e-2, a plot's accuracy; a benchmark
 *   wants less). The discretisation is second order, so the change estimates this grid's error
 *   (about 3/4 of it).
 * @param {import('./index.js').Device} device
 * @param {object} [sol] the solution to check (default: the device's present state)
 * @param {{ refine?: boolean, tol?: number, order?: boolean }} [opts]
 * Every check passed when `ok`; an item with `ok: null` (a warning, or a check that couldn't be
 * made) asks you to look, and doesn't fail it. `text` is a line per check, marked ok, FAIL or ?.
 * @returns {{ ok: boolean, items: CheckItem[], text: string }}
 */
export function check(device, sol, { refine = true, tol = 1e-2, order = true } = {}) {
  checkOptions('check', arguments[2], ['refine', 'tol', 'order']);
  sol ??= device.solution();
  const { model, solver, def } = device;
  const steady = sol.steady ?? solver.atSteady;
  const items = [];
  // Whether sol is the device's present state, and that state a solution: by the mark each
  // solution carries (see Device._stateOf), or, for a copy without one (or a twin's), by its
  // numbers. (By its numbers alone, a solution from before a set() of a material passed, as did
  // the device's own state after a set() of its drives, and its cold layout given as a solution;
  // Astra's fifth pass.) Unsolved, the state is its c0's layout, or one set() has changed the
  // equations of, which a check passed as "a transient at t = 0, converged" (Fable round 5's
  // review). (Where the numbers show how it moved on, they say so: "its terminal 'right' is at 0.5 V".)
  const state = device._stateOf(sol) ?? { moved: '', unsolved: device._stateOf(device.solution()).unsolved };
  const moved = state.moved === null ? null : movedOn(sol, atDrives(device)) ?? state.moved;
  const unsolved = !moved && state.unsolved;

  items.push(
    unsolved
      ? { name: 'converged', ok: false, summary: `not solved: ${unsolved}; solve() or advance() it, then check that solution` }
      : sol.converged
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

  // A solution from before the device moved on (a set(), another solve, a further advance). The
  // checks below solve the device as it stands, so they'd compare two different problems: a
  // 0.5 V solution checked after set() to 0.4 V read "FAIL grid … refine the grid".
  // (against the device's drives: after a set() of drives alone, its solution reads the state at
  // the voltages it was solved at, until a step or solve moves it, but it stands at the new ones)
  if (moved) {
    items.push({
      name: 'solution',
      ok: false,
      summary: `not the device's present state (${moved}), and the checks compare with the device as it stands: check a solution before the device moves on, or check(device) for its present one`,
    });
  }

  // The rest checks a solution, which an unconverged or unsolved state isn't.
  if (sol.converged && !unsolved && !moved) {
    const ledger = steady ? balance(device, sol) : null;
    // A reaction the steady state needs at equilibrium, left far from it (too slow to see: the
    // solution warns of it), isn't steady, however quietly its ledgers sum: "equilibrium" was
    // printed beside that warning.
    const stuck = steady ? solver.stuckReactions() : [];
    if (ledger && stuck.length > 0) {
      ledger.ok = false;
      const how = (s) => (s.charging ? ` ${sig(s.affinity * solver.VT * 1e3)} mV from its source, which a terminal that only charges reaches in a steady state` : s.face ? ` ${sig(s.affinity * solver.VT * 1e3)} mV apart across the face, where nothing passes in a steady state` : s.links !== undefined ? ` ${sig(s.affinity * solver.VT * 1e3)} mV from where its rate links balance` : ` ${s.immobile === undefined && s.closed === undefined ? '' : 'up to '}${sig(s.affinity * solver.VT * 1e3)} mV of affinity from equilibrium${s.immobile !== undefined ? ` where a participant can't move` : s.closed !== undefined ? ' where no reservoir reaches' : ''}`);
      ledger.summary = `not a steady state: ${stuck.map((s) => `${s.where} is${how(s)}`).join(', ')}, too slow for the solve to see (see the warnings); its ledgers say only that nothing moves fast enough to count`;
    }
    items.push(ledger ?? conservation(sol, Math.max(...device.model.grid.area)));
    if (steady && refine) items.push(gridCheck(device, sol, tol, ledger.details.gross));
    if (steady && refine && order && model.species.length > 1) items.push(orderCheck(device, sol, ledger.details.gross));
    // A transient's answer depends on its whole history (and on drives set along the way), so
    // only a run of its own on a finer grid can check it.
    if (!steady && refine) items.push({ name: 'grid', ok: null, summary: 'not checked in a transient: run it again on a grid twice as fine (hmin and hmax halved), and with a tighter tol (its steps shorter), and compare' });
  }

  const ok = sol.converged && items.every((it) => it.ok !== false);
  const mark = (it) => (it.ok === true ? 'ok  ' : it.ok === false ? 'FAIL' : '?   ');
  const text = items.map((it) => `${mark(it)} ${it.name}: ${it.summary}${it.details?.lines?.length ? `\n${it.details.lines.map((l) => `       ${l}`).join('\n')}` : ''}`).join('\n');
  return { ok, items, text };
}

// The units of totals through the cross-section: per m² of a planar device, or whole.
const unitsOf = (model) => (model.geometry.type === 'planar' ? { rate: 'mol/(m²·s)', current: 'A/m²', charge: 'C/m²', amount: 'mol/m²' } : { rate: 'mol/s', current: 'A', charge: 'C', amount: 'mol' });

// The device's present solution, its held terminals at their drives (see makeSolution's readV).
function atDrives(device) {
  const now = device.solution(), s = device.solver;
  if (!s.readV) return now;
  const terminals = { ...now.terminals };
  s.terms.forEach((t, k) => {
    if (t.drive.kind === 'V' && !(t.drive.R > 0)) terminals[t.name] = { ...terminals[t.name], V: sourceAt(t.drive.src, s.time, s.stepped) };
  });
  return { ...now, terminals };
}

// How a solution differs from the device's present state, in words, or null if it is that state.
function movedOn(sol, now) {
  if (sol.x.length !== now.x.length) return `it has ${sol.x.length} nodes, the device ${now.x.length}`;
  if (sol.time !== now.time) return `it is at t = ${sig(sol.time)} s, the device at ${sig(now.time)} s`;
  for (const [name, t] of Object.entries(now.terminals)) {
    const V = sol.terminals?.[name]?.V;
    if (!(Math.abs(V - t.V) <= 1e-9 * Math.max(1, Math.abs(t.V)))) return `its terminal '${name}' is at ${sig(V)} V, the device's at ${sig(t.V)} V`;
  }
  // (φ is NaN in a metal, and Math.max(d, NaN) masked every node after one: a copy of a solution
  // from before a set() passed; Astra's fifth pass)
  let d = 0;
  for (let g = 0; g < now.phi.length; g++) {
    const a = sol.phi[g], b = now.phi[g];
    if (Number.isNaN(a) && Number.isNaN(b)) continue;
    d = Math.max(d, Math.abs(a - b) || (a === b ? 0 : Infinity));
  }
  if (d > 1e-9) return `its φ differs from the device's by up to ${sig(d)} V`;
  // And each species' concentration: where φ is undefined throughout (neutral species alone), φ
  // compared nothing, and a copy from before a set() of a reaction's kf passed (Neon's review, 4.4).
  for (const [name, c] of Object.entries(now.c ?? {})) {
    const was = sol.c?.[name];
    if (!was) continue;
    for (let g = 0; g < c.length; g++) {
      const a = was[g], b = c[g];
      if (Number.isNaN(a) && Number.isNaN(b)) continue;
      if (!(Math.abs(a - b) <= 1e-9 * Math.max(Math.abs(a), Math.abs(b)) + 1e-300)) return `its ${name} differs from the device's (${a.toPrecision(6)} against ${b.toPrecision(6)} at node ${g})`;
    }
  }
  return null;
}

// Each species' steady ledger: in through the terminals, made by the reactions. Rates in
// mol/(m²·s) (mol/s if not planar), positive into the device or made.
function balance(device, sol) {
  const { model, def } = device, { species, regions, interfaces } = model, units = unitsOf(model);
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
  const partners = new Map(); // a reaction term's species: its net rate is round-off against any of theirs
  // (A term that isn't a number fails the ledger: dropped, a NaN flux at both contacts of a driven
  // wire read "equilibrium"; Astra's fifth pass.)
  const broken = [];
  const add = (i, r, what, rate, species) => {
    if (species) partners.set(what, species);
    const l = home[i][r];
    if (l >= 0 && !Number.isFinite(rate)) broken.push({ l, what, rate });
    if (!(l >= 0 && rate !== 0 && Number.isFinite(rate))) return;
    const same = ledgers[l].terms.find((t) => t.what === what); // a bulk reaction, region by region
    if (same) same.rate += rate;
    else ledgers[l].terms.push({ what, rate });
  };
  const last = regions.length - 1;
  // (A contact that blocks a species is no source of it: what its flux reads there is the end
  // box's own imbalance, which the ledger must show, not balance.)
  const through = (side, i) => model.contacts[side].species[i].type !== 'blocked';
  species.forEach((sp, i) => {
    if (through('left', i)) add(i, 0, 'left contact', sol.contacts.left.flux[sp.name]);
    if (through('right', i)) add(i, last, 'right contact', -sol.contacts.right.flux[sp.name]);
    model.ports.forEach((port, k) => add(i, port.region, `port ${sol.ports[k].name}`, sol.ports[k].flux[sp.name]));
  });
  // (two reactions with one equation, radiative and SRH recombination say, keep a line each,
  // named by their place in the list)
  const unique = (labels, where) => labels.map((w, k) => (labels.indexOf(w) === labels.lastIndexOf(w) ? w : `${w} (${where}[${k}])`));
  const bulkWhat = unique(model.reactions.map((_, k) => label(def.bulkReactions?.[k])), 'bulkReactions');
  model.reactions.forEach((rx, k) => {
    const what = bulkWhat[k];
    sol.bulkReactions[k].regions.forEach((total, r) => {
      const all = [...rx.reactants, ...rx.products].map((p) => p.i);
      for (const { i, nu } of rx.reactants) add(i, r, what, -nu * total, all);
      for (const { i, nu } of rx.products) add(i, r, what, nu * total, all);
    });
  });
  interfaces.forEach((itf, f) => {
    // A face's flux and rates are per area; everything else here is a total, through A there.
    const Af = model.grid.area[model.grid.regionEnd[f]];
    // Across the face, where it divides a species' compartments.
    species.forEach((sp, i) => {
      const type = itf.links[i].type;
      if (type === 'equilibrium' || type === 'blocked') return;
      const N = Af * sol.interfaces[f].N[sp.name];
      add(i, f, `across face ${f} (${type})`, -N);
      add(i, f + 1, `across face ${f} (${type})`, N);
    });
    const faceWhat = unique(itf.reactions.map((_, k) => label(def.interfaces?.[f]?.reactions?.[k])), 'reactions');
    itf.reactions.forEach((rx, k) => {
      const rate = Af * sol.interfaces[f].rates[k], what = `${faceWhat[k]} at face ${f}`;
      const all = rx.part.map((p) => p.i);
      for (const { i, nu, side } of rx.part) add(i, f + side, what, nu * rate, all);
    });
  });
  // A species in one compartment is named alone; one in several, by where.
  const several = species.map((_, i) => ledgers.filter((l) => l.i === i).length > 1);
  for (const l of ledgers) {
    const names = l.regions.map((r) => sol.regions[r].name);
    l.compartment = names.length === 1 ? names[0] : `${names[0]} to ${names.at(-1)}`;
    l.title = several[l.i] ? `${l.species} in ${l.compartment}` : l.species;
  }
  // Each surface species at a face or a port: made as fast as it's used, Σ ν·rate over the
  // reactions there (per area at a face, through the cross-section there; over a port's window).
  // Without them a snapshot whose surface was out of balance by 30% of its rates passed, every
  // volume ledger unchanged (Astra's fifth pass). Read against what the reactions' volume species
  // could carry, as their terms are.
  const surfaceLedger = (sp, where, partOf) => {
    const l = { species: sp.name, z: 0, i: -1, regions: [], terms: [], compartment: where, title: `${sp.name} at ${where}`, partners: partOf };
    const at = ledgers.push(l) - 1;
    return (what, rate) => {
      if (!Number.isFinite(rate)) broken.push({ l: at, what, rate });
      else if (rate !== 0) l.terms.push({ what, rate });
    };
  };
  interfaces.forEach((itf, f) => {
    if (!(itf.surface?.length > 0)) return;
    const Af = model.grid.area[model.grid.regionEnd[f]], faceWhat = unique(itf.reactions.map((_, k) => label(def.interfaces?.[f]?.reactions?.[k])), 'reactions');
    const partOf = [...new Set(itf.reactions.flatMap((rx) => rx.part.map((p) => p.i)))];
    itf.surface.forEach((sp, s) => {
      const put = surfaceLedger(sp, `face ${f}`, partOf);
      itf.reactions.forEach((rx, k) => {
        for (const u of rx.surf ?? []) if (u.s === s) put(`${faceWhat[k]} at face ${f}`, u.nu * Af * sol.interfaces[f].rates[k]);
      });
    });
  });
  model.ports.forEach((port, k) => {
    if (!(port.surface?.length > 0)) return;
    const sp = sol.ports[k], vol = model.grid.vol, portWhat = unique(port.reactions.map((_, j) => label(def.ports?.[k]?.reactions?.[j])), 'reactions');
    const partOf = [...new Set(port.reactions.flatMap((rx) => rx.part.filter((p) => p.side !== 2).map((p) => p.i)))];
    port.surface.forEach((surf, s) => {
      const put = surfaceLedger(surf, `port ${sp.name}`, partOf);
      port.reactions.forEach((rx, j) => {
        for (const u of rx.part) {
          if (u.side !== 2 || u.s !== s) continue;
          let total = 0;
          port.nodes.forEach((g, w) => (total += sp.rates[j][w] * sp.area[w] * vol[g]));
          put(`${portWhat[j]} at port ${sp.name}`, u.nu * total);
        }
      });
    });
  });

  // A species at rest (nothing beyond round-off moving, against what it or its reactions' partners
  // could carry) isn't reported.
  const size = (l) => Math.max(0, ...l.terms.map((t) => Math.abs(t.rate)));
  const natural = naturalScales(device, sol), inSolution = naturalScales(device, sol, false); // (partners' scale: not a metal's conductance)
  const scale = (l) => Math.max(l.i >= 0 ? natural[l.i] : 0, ...(l.partners ?? []).map((j) => inSolution[j]), ...l.terms.flatMap((t) => (partners.get(t.what) ?? []).map((j) => inSolution[j])));
  // Small flows count too where they balance, as real ones do to round-off (a diode at 1 µV,
  // 1e-19 of its scale, to 1e-15), and round-off doesn't. (Above 1e-8 of the scale, any.)
  const balanced = (l) => Math.abs(l.terms.reduce((a, t) => a + t.rate, 0)) <= 1e-9 * size(l);
  const moving = ledgers.filter((l) => size(l) > 1e-8 * scale(l) || (size(l) > 1e-20 * scale(l) && balanced(l)));
  // And a ledger with a reaction's term that a moving one carries as one of its own principal
  // terms, at a like rate (within the stoichiometry's factor): that rate is shown real there, so it
  // is here. (A catalytic face's surface species balanced its two steps at 9.9e-11 mol/(m²·s), 1e-10
  // of its scale, while the solution species the first step took, fed by nothing in a solve that
  // had dropped its flux, sat as "round-off" beside an ok; Neon's reply on 4.1.)
  for (let more = true; more; ) {
    more = false;
    for (const l of ledgers) {
      if (moving.includes(l)) continue;
      const shown = l.terms.some((t) => moving.some((m) => m.terms.some((u) => u.what === t.what && Math.abs(u.rate) >= 1e-3 * size(m) && Math.abs(t.rate) <= 10 * Math.abs(u.rate) && Math.abs(u.rate) <= 10 * Math.abs(t.rate))));
      if (shown) {
        moving.push(l);
        more = true;
      }
    }
  }
  // What's left unbalanced can be round-off or a small real flow the solve didn't settle (a
  // reverse diode behind a resistance, 1.3 µA/m² at 4e-15 of its scale, balanced to 8e-5 by a
  // solve that stalled), and nothing here tells them apart: round-off reaches 8e-12 of the scale
  // in the stress devices, and one in a thousand of those flows balances to 1%. So above ε of the
  // scale it isn't called equilibrium: the largest is named, as what balance can't check.
  const unsure = ledgers.filter((l) => !moving.includes(l) && size(l) > 1e-16 * scale(l)).sort((a, b) => size(b) / scale(b) - size(a) / scale(a))[0];
  let worst = 0;
  const lines = [];
  for (const l of moving) {
    const sum = l.terms.reduce((a, t) => a + t.rate, 0);
    l.residual = sum / size(l);
    worst = Math.max(worst, Math.abs(l.residual));
    // Terms that matter, largest first; a charged species' as a current too.
    const shown = l.terms.filter((t) => Math.abs(t.rate) > 1e-6 * size(l)).sort((a, b) => Math.abs(b.rate) - Math.abs(a.rate));
    const amp = (r) => (l.z !== 0 ? ` (${sig(FARADAY * Math.abs(l.z) * r)} ${units.current})` : '');
    lines.push(`${l.title}: ${shown.map((t) => `${t.what}: ${t.rate > 0 ? '+' : ''}${sig(t.rate)}${amp(t.rate)}`).join('; ')}`);
  }
  const ok = worst < 1e-6 && broken.length === 0;
  // (what it can't tell from round-off, said beside the ledgers that balance too: dropped there,
  // a neutral species' unbalanced 1e-15 went unmentioned; Astra's fifth pass)
  const who = unsure && `${unsure.title} at ${sig(size(unsure))} ${units.rate}${unsure.z !== 0 ? ` (${sig(FARADAY * Math.abs(unsure.z) * size(unsure))} ${units.current})` : ''}`;
  const caveat = unsure && `is ${sig(size(unsure) / scale(unsure), 1)} of what it could carry and doesn't balance (to ${sig(Math.abs(unsure.terms.reduce((a, t) => a + t.rate, 0)) / size(unsure), 1)}), so it's round-off, or a flow too small for this check`;
  const summary = broken.length > 0
    ? `not a number: ${broken.map((b) => `${ledgers[b.l].title}'s ${b.what} reads ${b.rate}`).join(', ')}, so its ledger can't be summed: something upstream of the solution is NaN or infinite`
    : moving.length === 0
      ? unsure
        ? `nothing flows that it can tell from round-off: the largest, ${who}, ${caveat}`
        : 'nothing flows or reacts (equilibrium)'
      : `${moving.length === 1 ? 'one ledger sums' : `${moving.length} ledgers sum`} to zero ${worst < 1e-6 ? `(to ${sig(worst, 1)})` : `only to ${sig(worst, 2)}: not a steady state`}: each species' sources and sinks, ${units.rate}, + in or made${unsure ? `; besides them, ${who} ${caveat}` : ''}`;
  // The largest current any charged species carries anywhere: the scale a net current is read against.
  const gross = Math.max(0, ...moving.filter((l) => l.z !== 0).map((l) => FARADAY * Math.abs(l.z) * size(l)));
  return { name: 'balance', ok, summary, details: { ledgers, worst, lines, gross } };
}

// What each species could carry, mol/(m²·s): by transport, D·c/L in each region (a metal's
// carrier, σRT/(z²F²L)), and by its bulk reactions, their one-way rates. In equilibrium the net
// flows are round-off of these. Without `conductors`, transport through conductor regions isn't
// counted.
function naturalScales(device, sol, conductors = true) {
  const { model, solver } = device, { species, grid } = model;
  // (totals through the cross-section, as the ledgers are: its largest in each region)
  const area = model.regions.map((_, r) => Math.max(...grid.area.subarray(grid.regionStart[r], grid.regionEnd[r] + 1)));
  return species.map((sp, i) => {
    let S = 0;
    model.regions.forEach((reg, r) => {
      const mat = model.materials[reg.material], L = reg.length;
      if (mat.conductor) {
        if (conductors && mat.conductor.i === i) S = Math.max(S, (area[r] * mat.conductor.sigma * model.RT) / (sp.z * sp.z * FARADAY * FARADAY * L));
        return;
      }
      if (!mat.present[i] || !(mat.D[i] > 0)) return;
      let cmax = 0;
      for (let g = 0; g < sol.x.length; g++) if (sol.region[g] === r && Number.isFinite(sol.c[sp.name][g])) cmax = Math.max(cmax, sol.c[sp.name][g]);
      const Dmax = reg.grade?.[i]?.D ? Math.max(...reg.grade[i].D.values) : mat.D[i];
      S = Math.max(S, (area[r] * Dmax * cmax) / L);
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

// The floor under a current the grid and order checks compare, A (or A/m²): 1e-16 of what the
// charged species could carry, above what solves resolve (~1e-19 of it in a GaAs junction). Not
// through a metal unless the device has nothing else to carry a current: a metal's current is an
// unknown of its own (Ohm's law in mixed form), as precise as what feeds it, and its conductance
// is no measure of its round-off. (Counted, a 1 cm Pt counter electrode's set a 304 nA floor
// under a 1 µm electrode's 1.96 nA, and a current 9.2% off on refinement passed at 6.0e-4;
// Neon's fourth review, 2.1.) Its `scales`, per species, are those it was taken from.
function currentFloor(device, sol, natural = naturalScales(device, sol)) {
  const { species } = device.model;
  const charged = (S) => Math.max(0, ...species.map((sp, i) => Math.abs(sp.z) * S[i]).filter(Number.isFinite));
  const inSolution = naturalScales(device, sol, false);
  const scales = charged(inSolution) > 0 ? inSolution : natural;
  return { floorI: 1e-16 * FARADAY * charged(scales), scales };
}

// Closed, unreacting stretches keep their amounts in a transient.
function conservation(sol, area = 1) {
  // Closed stretches (a species held between blocking faces) against what they held plus what
  // came in, and those reaching a terminal against that plus what came in through it (read only
  // for closed ones, a broken intake would have passed; Fable round 5's review).
  const closed = sol.conservation.filter((st) => !st.connected && !st.reactive);
  const open = sol.conservation.filter((st) => st.connected && !st.reactive);
  if (closed.length + open.length === 0) return { name: 'conservation', ok: true, summary: 'nothing to conserve (every species reacts)' };
  // (An unreacting stretch's drift is always a number: one that isn't fails. Left out, a NaN drift
  // in the only open stretch read "nothing to conserve (every species reacts)"; Astra's fifth pass.)
  const broken = [...closed, ...open].filter((st) => ![st.amount, st.reference, st.intake, st.drift].every(Number.isFinite));
  if (broken.length > 0) return { name: 'conservation', ok: false, summary: `not a number: ${broken.map((st) => `${st.species}'s stretch (${['amount', 'reference', 'intake', 'drift'].filter((k) => !Number.isFinite(st[k])).map((k) => `${k} ${st[k]}`).join(', ')})`).join('; ')}, so what it kept can't be read: something upstream of the solution is NaN or infinite`, details: { worst: NaN } };
  // Relative to the amount, but a drift below 1e-20 mol/m² (1e-15 C/m²) is nothing: a floating
  // metal's excess electrons are themselves ~0. (Amounts are totals through the cross-section.)
  const off = (st) => (Math.abs(st.amount - st.reference - st.intake) < 1e-20 * area ? 0 : Math.abs(st.drift));
  const worst = Math.max(...[...closed, ...open].map(off));
  const n = (k, one, many) => `${k} ${k === 1 ? one : many}`;
  const what = [closed.length ? `${n(closed.length, 'closed stretch', 'closed stretches')} (a species held between blocking faces)` : '', open.length ? n(open.length, 'stretch reaching a terminal', 'stretches reaching a terminal') : ''].filter(Boolean).join(' and ');
  return {
    name: 'conservation',
    ok: worst < 1e-8,
    summary: `${what} kept ${worst < 1e-8 ? `to ${sig(worst, 1)}` : `only to ${sig(worst, 2)}`} of what they held plus what came in`,
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
      const defined = (o) => Object.fromEntries(Object.entries(o ?? {}).filter(([, v]) => v !== undefined)); // (as buildGrid reads them)
      const g = { ...defined(def.grid), ...defined(reg.grid) };
      if (!gridGiven(def.grid) && !gridGiven(reg.grid)) Object.assign(g, defaultSpacing(reg.length));
      const out = { ratio: Math.sqrt(g.ratio ?? DEFAULTS.ratio), minCells: 2 * (g.minCells ?? DEFAULTS.minCells) };
      if (g.hmin !== undefined) out.hmin = g.hmin / 2;
      if (g.hmax !== undefined) out.hmax = g.hmax / 2;
      if (g.hLeft !== undefined) out.hLeft = g.hLeft / 2;
      if (g.hRight !== undefined) out.hRight = g.hRight / 2;
      return { ...reg, grid: out };
    }),
  };
  // Started from this solution, interpolated onto the finer grid; failing that, from cold.
  let other, otherModel;
  try {
    const dev = new Device(fine);
    carryState(device.solver, dev.solver, { oldDef: device.def, newDef: fine });
    other = dev.solve();
    otherModel = dev.model;
    // (From cold, the device holds its definition's amounts: where this state's own came from its
    // history (a transient, a floating terminal's charge, a closed stretch's amount), the cold
    // solve is of another state, and a port charged by a run read "FAIL grid … refine the grid",
    // 3.7e-3 → 0 C/m²; Astra's fifth pass.)
    const s = device.solver, history = s.time > 0 || s.stepped || s.history.length > 0 || s.floating.length > 0 || sol.conservation.some((st) => !st.connected);
    if (!other.converged && history) return { name: 'grid', ok: null, summary: "couldn't solve from this state on a twice-finer grid, and from cold it would hold other amounts than this state's history left it: check resolution by hand" };
    if (!other.converged) {
      const cold = new Device(fine);
      if (sol.time > 0) cold.solver.time = sol.time; // waveforms at the same time
      other = cold.solve();
      otherModel = cold.model;
    }
  } catch (e) {
    return { name: 'grid', ok: null, summary: `couldn't solve the refined device (${e.message})` };
  }
  if (!other.converged) return { name: 'grid', ok: null, summary: "the device on a twice-finer grid didn't solve; check resolution by hand" };

  const diffs = [];
  // The current, against the largest current any species carries (so at open circuit, or in
  // equilibrium, round-off isn't read as a change).
  // A current too small for the ledgers to call moving (a diode at low bias: 1e-4 A/m² against
  // majority carriers that could carry 6e8) is compared still, above a floor at 1e-16 of what the
  // species could carry outside metals (currentFloor).
  const natural = naturalScales(device, sol), { floorI, scales } = currentFloor(device, sol, natural);
  // (gross only as a floor, a millionth of it: as the scale, a port's 1 µA/m² doubled beside a
  // contact's 3.8 A/m² changed by 2.9e-7 and passed; Astra's fifth pass)
  const scaleI = Math.max(Math.abs(sol.current), Math.abs(other.current), 1e-6 * gross, floorI);
  if (gross > 0 || Math.max(Math.abs(sol.current), Math.abs(other.current)) > floorI) diffs.push({ what: 'current', a: sol.current, b: other.current, rel: Math.abs(other.current - sol.current) / scaleI, unit: unitsOf(model).current });
  // Every other terminal's current too, and what each species brings in at the contacts (a
  // neutral one's carries no current): the right contact's alone passed a port's current 20% off
  // and a neutral reactant's flux 20× off, "nothing to compare" (Fable round 5's review).
  for (const t of model.terminals) {
    if (t.name === 'right') continue;
    const a = sol.terminals[t.name]?.current, b = other.terminals[t.name]?.current;
    if (!Number.isFinite(a) || !Number.isFinite(b) || !(Math.max(Math.abs(a), Math.abs(b), gross) > floorI)) continue;
    diffs.push({ what: `${t.name} current`, a, b, rel: Math.abs(b - a) / Math.max(Math.abs(a), Math.abs(b), 1e-6 * gross, floorI), unit: unitsOf(model).current });
  }
  const floorN = floorI / FARADAY;
  // Which species a reaction touches (bulk, at a face, at a port): their fluxes are what the
  // reactions turn over, compared as they are; a spectator's against what it could carry too.
  const reacts = new Set();
  for (const rx of device.solver.rxs ?? []) for (const i of rx.sp) reacts.add(i);
  for (const host of [...model.interfaces, ...model.ports]) for (const rx of host.reactions ?? []) for (const p of rx.part ?? []) reacts.add(p.i);
  for (const side of ['left', 'right']) {
    const fa = sol.contacts?.[side]?.flux ?? {}, fb = other.contacts?.[side]?.flux ?? {};
    for (const sp of model.species) {
      // (The charged ones too, not only through the currents: where they cancel (ambipolar flow,
      // open circuit, a symmetric cell), nothing compared moved, and a lit n⁺-Si slab's
      // collected holes 5% off, a salt reacting at its contacts 38× off, passed; Neon's review, 4.2.)
      const a = fa[sp.name], b = fb[sp.name], i = model.species.indexOf(sp);
      // (a neutral one's floor its own, 1e-16 of what it could carry, not the charged carriers'
      // floor in moles: a metal beside it raised that a thousandfold; Neon's reply on 4.3)
      const floor = sp.z === 0 ? (natural[i] > 0 ? 1e-16 * natural[i] : floorN) : Math.max(floorN, (1e-6 * gross) / (FARADAY * Math.abs(sp.z)));
      if (!Number.isFinite(a) || !Number.isFinite(b) || !(Math.max(Math.abs(a), Math.abs(b)) > floor)) continue;
      // (a spectator's against 1e-3 of what it could carry, too: its 1e-5 of it, unresolved in a
      // thin reaction layer, failed a device whose current was right to 1e-3; the re-check. Not
      // a species a reaction turns over: a reactant's flux 12.7% off read as 4.7e-3; Neon 4.3.)
      const scale = Math.max(Math.abs(a), Math.abs(b), floor, reacts.has(i) ? 0 : 1e-3 * (natural[i] || 0));
      diffs.push({ what: `${sp.name} in at the ${side} contact`, a, b, rel: Math.abs(b - a) / scale, unit: model.geometry.type === 'planar' ? 'mol/(m²·s)' : 'mol/s' });
    }
  }
  // What each port brings in, per species, as a contact's is; and what each reaction turns over:
  // a bulk one's total, a face's rate, a port's summed over its window. (A port feeding a fast
  // A = B, on 8 cells, passed 12.5× the fine grid's A with the contacts' fluxes, all compared, at
  // 5.6e-16; Neon's second review, R2 2.2.)
  const fluxUnit = model.geometry.type === 'planar' ? 'mol/(m²·s)' : 'mol/s';
  const speciesFloor = (sp, i) => (sp.z === 0 ? (natural[i] > 0 ? 1e-16 * natural[i] : floorN) : Math.max(floorN, (1e-6 * gross) / (FARADAY * Math.abs(sp.z))));
  model.ports.forEach((port, k) => {
    const fa = sol.ports?.[k]?.flux ?? {}, fb = other.ports?.[k]?.flux ?? {};
    model.species.forEach((sp, i) => {
      const a = fa[sp.name], b = fb[sp.name], floor = speciesFloor(sp, i);
      if (!Number.isFinite(a) || !Number.isFinite(b) || !(Math.max(Math.abs(a), Math.abs(b)) > floor)) return;
      const scale = Math.max(Math.abs(a), Math.abs(b), floor, reacts.has(i) ? 0 : 1e-3 * (natural[i] || 0));
      diffs.push({ what: `${sp.name} in at port ${port.name}`, a, b, rel: Math.abs(b - a) / scale, unit: fluxUnit });
    });
  });
  // (a reaction's floor, too, from what its species carry outside metals: at a metal's face, its
  // carrier's conductance set the floor of the face's rate; Neon's fourth review, 2.1)
  const rateFloor = (idx) => Math.max(floorN, ...idx.map((i) => (scales[i] > 0 ? 1e-16 * scales[i] : 0)));
  const turnover = (what, a, b, idx) => {
    const floor = rateFloor(idx);
    if (!Number.isFinite(a) || !Number.isFinite(b) || !(Math.max(Math.abs(a), Math.abs(b)) > floor)) return;
    diffs.push({ what, a, b, rel: Math.abs(b - a) / Math.max(Math.abs(a), Math.abs(b), floor), unit: fluxUnit });
  };
  (device.solver.rxs ?? []).forEach((rx, k) => turnover(`bulk reaction ${rx.equation ?? k}'s total`, sol.bulkReactions?.[k]?.total, other.bulkReactions?.[k]?.total, Array.from(rx.sp)));
  // Every quantity here is through the whole cross-section, as the contacts', ports' and
  // regions' are: a face's rates and D are per area, so times the face's area. (Per area beside
  // totals, a 10 µm sphere's face D, 9.6e-5 C/m² against its 1.5e-13 C, set a floor that passed
  // the outer contact's D 4.9% off; Neon's third review, 1.3.)
  const faceArea = (m, f) => m.grid.area[m.grid.regionEnd[f]];
  model.interfaces.forEach((itf, f) => (itf.reactions ?? []).forEach((rx, k) => turnover(`face ${f}'s reaction ${rx.equation ?? k}`, faceArea(model, f) * sol.interfaces?.[f]?.rates?.[k], faceArea(otherModel, f) * other.interfaces?.[f]?.rates?.[k], rx.part.map((p) => p.i))));
  // (by the port's own nodes: by coordinate, a face's two nodes at one x collided, and a port
  // reaching the face read its neighbour's box, a 1 mm region's beside a 1 µm one's, 100× the
  // total; Neon's third review, 1.2)
  const windowTotal = (s, m, k, q) => {
    const p = s.ports?.[k], nodes = m.ports[k]?.nodes;
    if (!p?.rates?.[q] || !p.area || !nodes || nodes.length !== p.area.length) return NaN;
    let t = 0;
    nodes.forEach((g, w) => (t += p.rates[q][w] * p.area[w] * m.grid.vol[g]));
    return t;
  };
  model.ports.forEach((port, k) => (port.reactions ?? []).forEach((rx, q) => turnover(`port ${port.name}'s reaction ${rx.equation ?? q}`, windowTotal(sol, model, k, q), windowTotal(other, otherModel, k, q), rx.part.map((p) => p.i))));
  // Floating terminals (driven by a current): their voltages, against the span of the terminals'
  // and the thermal voltage. (Against the larger voltage, the verdict hung on where 0 V was: the
  // same device shifted by 1 V read a 40× smaller change; Astra's fifth pass.)
  const VT = model.RT / FARADAY, spanV = Math.max(VT, span(sol), span(other));
  model.terminals.forEach((t, k) => {
    if (device.solver.floating.includes(k)) {
      const a = sol.terminals[t.name].V, b = other.terminals[t.name].V;
      diffs.push({ what: `${t.name} voltage`, a, b, rel: Math.abs(b - a) / spanV, unit: 'V' });
    }
  });
  // Each region's charge, against the largest: summed over the nodes' boxes, as the solver's
  // Gauss law counts it, with what its faces' diffuse layers hold, below the grid (as the order
  // check counts amounts; left out, an electrode's 0.04 C/m² went uncompared; Astra's fifth pass).
  const layers = (s, grid) => s.interfaces.flatMap((itf, f) => Object.entries(itf.diffuse ?? {}).map(([side, d]) => ({ f, side, r: side === 'left' ? f : f + 1, q: grid.area[side === 'left' ? grid.regionEnd[f] : grid.regionStart[f + 1]] * d.charge })));
  const charges = (s, grid) => {
    const q = s.regions.map((_, r) => {
      let t = 0;
      for (let g = 0; g < s.x.length; g++) if (s.region[g] === r) t += rho(s, g, model, r) * grid.vol[g];
      return t;
    });
    for (const l of layers(s, grid)) q[l.r] += l.q;
    return q;
  };
  // The charges held where the grid is coarsest: each region's (summed over its boxes), each
  // face's diffuse layers, each contact's D, a capacitive face's, each port's. Against
  // themselves, with a floor of a thousandth of the largest and of 1e-9 of the gross charge the
  // regions hold, ∫(|fixed| + F Σ|z|c) dV, below which a strictly neutral region's charge is
  // round-off. (Gated on a millionth of that gross, max |q| × V, 1 mm of 1 M KCl skipped its
  // electrode's 0.058 C/m², 6.6% off the fine grid's; and its region's net charge alone, a
  // symmetric cell's two walls cancel, D 14% off; both grid ok, Neon's second review, R2 2.1.)
  let content = 0;
  for (let g = 0; g < sol.x.length; g++) {
    let a = Math.abs(fixedAt(model.regions[sol.region[g]], sol.x[g]));
    for (const sp of model.species) if (Number.isFinite(sol.c[sp.name][g])) a += FARADAY * Math.abs(sp.z * sol.c[sp.name][g]);
    content += a * model.grid.vol[g];
  }
  // (a contact's D and a port's charge are totals already; a face's D is per area: times its area)
  const held = (s, m) => [
    ...['left', 'right'].map((side) => ({ what: `the ${side} contact's D`, q: model.contacts[side].phi.type === 'neutral' ? NaN : s.contacts?.[side]?.D })),
    ...model.interfaces.map((itf, f) => ({ what: `face ${f}'s D`, q: itf.phi.type === 'capacitive' ? faceArea(m, f) * s.interfaces?.[f]?.D : NaN })),
    ...model.ports.map((port, k) => ({ what: `port ${port.name}'s charge`, q: port.capacitance ? s.ports?.[k]?.charge : NaN })),
  ];
  const la = layers(sol, model.grid), lb = layers(other, otherModel.grid), ha = held(sol, model), hb = held(other, otherModel);
  const qa = charges(sol, model.grid), qb = charges(other, otherModel.grid);
  const finite = (x) => (Number.isFinite(x) ? Math.abs(x) : 0);
  const scaleQ = Math.max(0, ...qa.map(finite), ...qb.map(finite), ...la.map((l) => finite(l.q)), ...lb.map((l) => finite(l.q)), ...ha.map((h) => finite(h.q)), ...hb.map((h) => finite(h.q)));
  const floorQ = Math.max(1e-3 * scaleQ, 1e-9 * content);
  if (scaleQ > floorQ) {
    const named = (reg) => (/^region \d+$/.test(reg.name) ? `${reg.name} (${reg.material})` : reg.name);
    // (against itself, the largest a thousandth as a floor: against the largest alone, a small
    // region's charge could double and pass, as a small terminal's current did; Astra's fifth pass)
    const rel = (a, b) => Math.abs(b - a) / Math.max(Math.abs(a), Math.abs(b), floorQ);
    qa.forEach((a, r) => diffs.push({ what: `${named(sol.regions[r])} charge`, a, b: qb[r], rel: rel(a, qb[r]), unit: unitsOf(model).charge }));
    // (each layer too: an electrode cell's two layers cancel in the electrolyte's charge)
    la.forEach((l, j) => diffs.push({ what: `face ${l.f}'s diffuse layer charge (${l.side})`, a: l.q, b: lb[j].q, rel: rel(l.q, lb[j].q), unit: unitsOf(model).charge }));
    ha.forEach((h, j) => {
      if (Number.isFinite(h.q) && Number.isFinite(hb[j].q)) diffs.push({ what: h.what, a: h.q, b: hb[j].q, rel: rel(h.q, hb[j].q), unit: unitsOf(model).charge });
    });
  }
  const worst = diffs.reduce((w, d) => (d.rel > (w?.rel ?? -1) ? d : w), null);
  if (!worst) return { name: 'grid', ok: null, summary: 'nothing to compare (no current, flux, charge or floating voltage): the grid is unchecked' };
  const ok = worst.rel <= tol;
  return {
    name: 'grid',
    ok,
    summary:
      `on a grid twice as fine (${sol.x.length} → ${other.x.length} nodes), the largest change is ${worst.what}, ${sig(worst.a, 5)} → ${sig(worst.b, 5)} ${worst.unit} ` +
      `(${sig(worst.rel, 2)}${ok ? '' : `, over ${sig(tol, 1)}: refine the grid (hmin, hmax), most likely where profiles are steepest`})`,
    details: { diffs, nodes: [sol.x.length, other.x.length], grids: fine.regions.map((reg) => reg.grid) }, // (grids: the finer one's, by region)
  };
}

// The device solved from its definition with its species listed as given and in reverse: the
// same steady state either way, or one of them isn't steady. (Both from the definition, not this
// state, whose closed amounts a transient or set() may have moved from c0.)
function orderCheck(device, sol, gross) {
  const { def, model } = device;
  const solveAt = (d) => {
    const dev = new Device(d);
    if (sol.time > 0) dev.solver.time = sol.time; // waveforms at the same time
    return dev.solve();
  };
  let a, b;
  try {
    a = solveAt(def);
    b = solveAt({ ...def, species: [...def.species].reverse() });
  } catch (e) {
    return { name: 'order', ok: null, summary: `could not solve the device with its species reversed (${e.message})` };
  }
  if (!a.converged || !b.converged) {
    const which = a.converged ? 'listed in reverse' : b.converged ? 'listed as given' : 'listed either way';
    return { name: 'order', ok: null, summary: `the device solved from its definition did not converge with its species ${which}: compare with a long advance() by hand` };
  }
  const VT = model.RT / FARADAY, diffs = [];
  // (Currents only where the ledgers see something move: in equilibrium they're round-off of
  // whatever exchanges fastest, an electrode's couple, and a wrong state shows in the voltages
  // and amounts instead.)
  // (and above 1e-16 of what the species could carry outside metals, the grid check's floor: a
  // drift at 1e-16 A/m², still moving at 1e12 s, isn't worth reporting)
  const { floorI } = currentFloor(device, sol);
  // (Voltages against the span of the terminals', as the grid check's; currents against their own
  // size, the gross current a floor: against 1e-3 of the larger of the two, a change could reach
  // only ~1e3, and a current alone, 0.38 → 1 A/m², could never fail; Astra's fifth pass.)
  const spanV = Math.max(VT, span(a), span(b));
  for (const [name, t] of Object.entries(a.terminals)) {
    const u = b.terminals[name];
    diffs.push({ what: `${name} voltage`, a: t.V, b: u.V, rel: Math.abs(u.V - t.V) / (1e-6 * spanV), unit: 'V' });
    if (gross > 0 && Math.max(Math.abs(t.current), Math.abs(u.current), gross) > floorI) diffs.push({ what: `${name} current`, a: t.current, b: u.current, rel: Math.abs(u.current - t.current) / Math.max(1e-6 * Math.max(Math.abs(t.current), Math.abs(u.current)), 1e-8 * gross, floorI), unit: unitsOf(model).current });
  }
  // (with what the faces' diffuse layers hold, below the grid)
  const vol = model.grid.vol, area = model.grid.area;
  const held = (s, name) => s.interfaces.reduce((m, itf, f) => m + Object.entries(itf.diffuse ?? {}).reduce((t, [side, d]) => t + area[side === 'left' ? model.grid.regionEnd[f] : model.grid.regionStart[f + 1]] * (d.excess[name] ?? 0), 0), 0);
  const amount = (s, name) => s.x.reduce((m, _, g) => m + (Number.isFinite(s.c[name][g]) ? vol[g] * s.c[name][g] : 0), 0) + held(s, name);
  const amounts = model.species.map((sp) => [sp.name, amount(a, sp.name), amount(b, sp.name)]);
  const most = Math.max(0, ...amounts.map(([, x, y]) => Math.max(Math.abs(x), Math.abs(y))));
  for (const [name, x, y] of amounts) diffs.push({ what: `${name} amount`, a: x, b: y, rel: Math.abs(y - x) / (1e-6 * Math.max(Math.abs(x), Math.abs(y)) + 1e-12 * most), unit: unitsOf(model).amount });
  // Agreeing to near round-off (a µV, a millionth of an amount): ok. Within a thousand times that
  // (a mV, a thousandth), enough for any plot, but a weakly held part of one solve not quite
  // steady: look. Beyond: one of them isn't a steady state.
  const worst = diffs.reduce((w, d) => (d.rel > (w?.rel ?? -1) ? d : w), null);
  const ok = !(worst.rel > 1) ? true : worst.rel > 1e3 ? false : null;
  const what = `${worst.what} ${sig(worst.a, 6)} → ${sig(worst.b, 6)} ${worst.unit}`;
  return {
    name: 'order',
    ok,
    summary:
      ok === true
        ? 'the same steady state with the species listed in reverse'
        : ok === null
          ? `listed in reverse, the species give nearly the same steady state (${what}): not to round-off, so some weakly held part is not quite steady in one of them`
          : `listed in reverse, the species give a different steady state: ${what}. One of them is not steady, though its solve converged, or the device has more than one steady state (one of them perhaps unstable): compare with a long advance() from the start, kick each and watch it with dtMax (see numerics.md, Steady state), and please report the device.`,
    details: { diffs },
  };
}

// The span of a solution's terminal voltages, V: the scale a voltage is read against, whatever is 0 V.
function span(s) {
  const V = Object.values(s.terminals ?? {}).map((t) => t.V).filter(Number.isFinite);
  return V.length > 0 ? Math.max(...V) - Math.min(...V) : 0;
}

// Charge density at a node of region r, C/m³: the species plus the region's fixed charge.
function rho(s, g, model, r) {
  let q = fixedAt(model.regions[r], s.x[g]);
  for (const sp of model.species) {
    const c = s.c[sp.name][g];
    if (Number.isFinite(c)) q += FARADAY * sp.z * c;
  }
  return q;
}
