// Plot-ready traces of a solution, for driftlet/kit: the level diagram's lines with labels,
// line kinds and a colour slot that follows the species, the regions as bands, and a suggested
// range. Plain data, for any plotting toolkit (driftlet/plot draws it as SVG).

import { level } from './equation.js';

const finite = (a) => {
  for (let k = 0; k < a.length; k++) if (Number.isFinite(a[k])) return true;
  return false;
};

/**
 * One line of a level diagram.
 * @typedef {object} Trace
 * @property {string} id e.g. 'V:e-', 'Vstd:e-', 'phi', 'level:0'
 * @property {string} label e.g. 'V e-', 'V° e-'
 * @property {'level' | 'standard' | 'phi' | 'redox' | 'redox-standard'} kind a species voltage, its
 *   standard level, φ, or a half-reaction's level or standard level
 * @property {string} [species]
 * @property {'electron' | 'cation' | 'anion' | 'redox' | 'phi'} role what it is, for its colour:
 *   electrons, cations (holes too), anions, redox levels, φ
 * @property {number} slot its colour within the role: a species' place among the solution's
 *   species of that role, so a colour follows its species whichever lines are shown; a redox
 *   level's place among the distinct half-reactions asked for (a couple's level and standard
 *   level share one)
 * @property {number} [shift] a display offset added to all of a species' lines, V
 * @property {Float64Array} y V, NaN where undefined (a break in the line)
 */

/**
 * A charged species' colour role and slot in level diagrams: electrons (`e-`), cations (holes
 * too), or anions, and its place among the solution's species of that role.
 * @param {import('./types.js').Solution} sol
 * @param {string} name
 * @returns {{ role: 'electron' | 'cation' | 'anion', slot: number }}
 */
export function speciesRole(sol, name) {
  const z = Object.fromEntries(sol.species.map((sp) => [sp.name, sp.z]));
  const role = (n) => (n === 'e-' ? 'electron' : z[n] > 0 ? 'cation' : 'anion');
  const charged = sol.species.filter((sp) => sp.z !== 0).map((sp) => sp.name);
  return { role: role(name), slot: charged.filter((n) => role(n) === role(name)).indexOf(name) };
}

/**
 * The level diagram of a solution as traces: per charged species, its voltage V_i (`level`) and
 * its standard level V°_i (`standard`, the band edge for e⁻ and h⁺), optionally φ and
 * half-reactions' levels. Doubled interface nodes share an x, so steps draw as vertical lines.
 * @param {import('./types.js').Solution} sol
 * @param {{ species?: string[], standard?: boolean, phi?: boolean,
 *   levels?: { half: { equation: string, fixed?: Record<string, number> }, label?: string, standard?: boolean }[],
 *   labels?: Record<string, string>, shifts?: Record<string, number> }} [opts]
 *   species to show (default: every charged one), whether to show standard levels (default true)
 *   and φ (default false), half-reaction levels, labels to use by trace id
 *   ({ 'V:e-': 'Fermi level' }), and display offsets per species (V), which move all of a
 *   species' lines together so that widely separated species can share one readable plot
 * @returns {{ x: Float64Array, series: Trace[], regions: { name: string, material: string, x0: number, x1: number }[],
 *   faces: number[], range: [number, number] }}
 */
export function traces(sol, { species, standard = true, phi = false, levels = [], labels = {}, shifts = {} } = {}) {
  const names = Object.keys(sol.V);
  const role = (name) => speciesRole(sol, name).role;
  const slotOf = (name) => speciesRole(sol, name).slot;
  const shown = species ?? names.filter((name) => finite(sol.V[name]));
  const series = [];
  for (const name of shown) {
    if (!names.includes(name)) throw new Error(`traces: no charged species '${name}' in the solution`);
    const shift = shifts[name] ?? 0;
    const at = (y) => (shift ? y.map((v) => v + shift) : y);
    const common = { species: name, role: role(name), slot: slotOf(name), ...(shift ? { shift } : {}) };
    if (finite(sol.V[name])) series.push({ id: `V:${name}`, label: `V ${name}`, kind: 'level', ...common, y: at(sol.V[name]) });
    if (standard && finite(sol.Vstd[name])) series.push({ id: `Vstd:${name}`, label: `V° ${name}`, kind: 'standard', ...common, y: at(sol.Vstd[name]) });
  }
  if (phi) series.push({ id: 'phi', label: 'φ', kind: 'phi', role: 'phi', slot: 0, y: sol.phi });
  const couples = [...new Set(levels.map((lv) => lv.half.equation))];
  levels.forEach((lv, k) => {
    const y = level(sol, lv.half, { standard: lv.standard });
    series.push({
      id: `level:${k}`,
      label: lv.label ?? `${lv.standard ? 'V° ' : ''}${lv.half.equation}`,
      kind: lv.standard ? 'redox-standard' : 'redox',
      role: 'redox',
      slot: couples.indexOf(lv.half.equation),
      y,
    });
  });
  for (const s of series) {
    if (labels[s.id] !== undefined) s.label = labels[s.id];
    if (s.shift) s.label += ` ⌇${s.shift > 0 ? '+' : '−'}${Math.abs(s.shift)} V`;
  }
  let lo = Infinity, hi = -Infinity;
  for (const s of series) {
    for (const v of s.y) {
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
  }
  if (!(hi > lo)) [lo, hi] = Number.isFinite(lo) ? [lo - 0.5, lo + 0.5] : [-1, 1];
  const pad = 0.06 * (hi - lo);
  const regions = sol.regions.map((r) => ({ ...r }));
  return { x: sol.x, series, regions, faces: regions.slice(1).map((r) => r.x0), range: [lo - pad, hi + pad] };
}
