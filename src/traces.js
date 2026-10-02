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
 * @property {'level' | 'standard' | 'phi' | 'redox'} kind a species voltage (solid), its standard
 *   level (dashed), φ, or a half-reaction's level
 * @property {string} [species]
 * @property {number} slot colour slot: the species' index in the solution, so a colour follows
 *   its species whichever lines are shown; then φ, then the half-reaction levels in order
 * @property {Float64Array} y V, NaN where undefined (a break in the line)
 */

/**
 * The level diagram of a solution as traces: per charged species, its voltage V_i (`level`) and
 * its standard level V°_i (`standard`, the band edge for e⁻ and h⁺), optionally φ and
 * half-reactions' levels. Doubled interface nodes share an x, so steps draw as vertical lines.
 * @param {import('./types.js').Solution} sol
 * @param {{ species?: string[], standard?: boolean, phi?: boolean,
 *   levels?: { half: { equation: string, fixed?: Record<string, number> }, label?: string, standard?: boolean }[],
 *   labels?: Record<string, string> }} [opts]
 *   species to show (default: every charged one), whether to show standard levels (default true)
 *   and φ (default false), half-reaction levels, and labels to use by trace id ({ 'V:e-': 'Fermi level' })
 * @returns {{ x: Float64Array, series: Trace[], regions: { name: string, material: string, x0: number, x1: number }[],
 *   faces: number[], range: [number, number] }}
 */
export function traces(sol, { species, standard = true, phi = false, levels = [], labels = {} } = {}) {
  const names = Object.keys(sol.V);
  const shown = species ?? names.filter((name) => finite(sol.V[name]));
  const series = [];
  for (const name of shown) {
    const slot = names.indexOf(name);
    if (slot < 0) throw new Error(`traces: no charged species '${name}' in the solution`);
    if (finite(sol.V[name])) series.push({ id: `V:${name}`, label: `V ${name}`, kind: 'level', species: name, slot, y: sol.V[name] });
    if (standard && finite(sol.Vstd[name])) series.push({ id: `Vstd:${name}`, label: `V° ${name}`, kind: 'standard', species: name, slot, y: sol.Vstd[name] });
  }
  // φ and the half-reaction levels have slots after the species, the same whether φ is shown.
  if (phi) series.push({ id: 'phi', label: 'φ', kind: 'phi', slot: names.length, y: sol.phi });
  levels.forEach((lv, k) => {
    const y = level(sol, lv.half, { standard: lv.standard });
    series.push({ id: `level:${k}`, label: lv.label ?? `${lv.standard ? 'V° ' : ''}${lv.half.equation}`, kind: 'redox', slot: names.length + 1 + k, y });
  });
  for (const s of series) if (labels[s.id] !== undefined) s.label = labels[s.id];
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
