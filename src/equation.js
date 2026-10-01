// Reaction equations as text, for driftlet/kit: 'Ag+ + e- = Ag(s)', '2 H+ + 2 e- = H2',
// 'e- + h+ = 0'. They compile to the signed stoichiometry of the plain spec (ν < 0 consumed by
// the forward reaction, left to right as written).

import { DeviceError } from './device.js';
import { FARADAY } from './constants.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const isObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);

/**
 * One participant of an equation as written.
 * @typedef {object} Term
 * @property {string} name the participant as written (with any trailing parenthesis)
 * @property {number} nu signed coefficient (< 0 on the left of the equation)
 */

/**
 * Parse an equation into its terms, in order. The two sides are separated by `=` (or `⇌`, `<=>`,
 * `->`, `→`), terms by ` + ` (spaces around the plus, which keeps it apart from charges), and a
 * side with nothing on it is `0` or `∅`. A leading integer is a coefficient (`2 e-` or `2e-`).
 * @param {string} text
 * @param {string} [path] for error messages
 * @returns {Term[]}
 */
export function parseEquation(text, path = 'equation') {
  if (typeof text !== 'string') fail(`${path} must be a string like 'Ag+ + e- = Ag(s)'`);
  const sides = text.trim().split(/\s+(?:=|⇌|<=>|->|→)\s+/);
  if (sides.length !== 2) fail(`${path}: ${JSON.stringify(text)} needs exactly one '=' with a space on either side`);
  const terms = [];
  sides.forEach((side, s) => {
    if (side === '0' || side === '∅') return;
    for (const raw of side.split(/\s+\+\s+/)) {
      const m = /^(?:(\d+)\s*)?(\S+)$/.exec(raw.trim());
      if (!m) fail(`${path}: can't read ${JSON.stringify(raw.trim())} in ${JSON.stringify(text)} (terms are separated by ' + ', with spaces)`);
      const nu = m[1] === undefined ? 1 : Number(m[1]);
      if (nu === 0) fail(`${path}: a zero coefficient in ${JSON.stringify(text)}`);
      terms.push({ name: m[2], nu: s === 0 ? -nu : nu });
    }
  });
  if (terms.length === 0) fail(`${path}: ${JSON.stringify(text)} has no participants`);
  return terms;
}

/**
 * The signed stoichiometry `{ name: ν }` of a homogeneous reaction, as a bulk reaction's `nu`.
 * A participant on both sides is an error: give the net reaction.
 * @param {string} text
 * @param {string} [path]
 * @returns {Record<string, number>}
 */
export function stoichiometry(text, path = 'equation') {
  const nu = {};
  const side = {};
  for (const t of parseEquation(text, path)) {
    if (side[t.name] !== undefined && side[t.name] !== Math.sign(t.nu)) {
      fail(`${path}: '${t.name}' appears on both sides of ${JSON.stringify(text)}; give the net reaction`);
    }
    side[t.name] = Math.sign(t.nu);
    nu[t.name] = (nu[t.name] ?? 0) + t.nu;
  }
  return nu;
}

/**
 * A half-reaction: an equation with electrons `e-` on one side, and the μ (J/mol) of any
 * participant that isn't a species (a fixed-activity neutral, such as a solid metal). It's plain
 * data, `{ equation, fixed }`: spread it into a face reaction with its kinetics,
 * `{ ...silver, k0: 1e-3 }`, or pass it to `level()`.
 * @param {string} equation e.g. 'Ag+ + e- = Ag(s)'
 * @param {Record<string, number>} [fixed] μ of the fixed participants, J/mol
 * @returns {{ equation: string, fixed?: Record<string, number> }}
 */
export function half(equation, fixed) {
  const terms = parseEquation(equation, 'half');
  const electrons = terms.filter((t) => t.name === 'e-');
  if (electrons.length !== 1) fail(`half: ${JSON.stringify(equation)} needs electrons, 'e-', on exactly one side`);
  if (fixed !== undefined && !isObject(fixed)) fail('half: fixed must map participants to their μ (J/mol)');
  for (const [name, mu] of Object.entries(fixed ?? {})) {
    if (!terms.some((t) => t.name === name)) fail(`half: fixed.${name} isn't in ${JSON.stringify(equation)}`);
    if (!Number.isFinite(mu)) fail(`half: fixed.${name} must be a finite μ (J/mol)`);
  }
  return fixed === undefined ? { equation } : { equation, fixed: { ...fixed } };
}

/**
 * The standard hydrogen electrode, 2 H⁺ + 2 e⁻ ⇌ H₂, with μ(H₂) = 0 (the element in its standard
 * state, on the usual table convention). Needs `H+` among the species.
 */
export const SHE = Object.freeze(half('2 H+ + 2 e- = H2', { H2: 0 }));

/**
 * The electronic level a half-reaction implies at each node, as an electron voltage
 * V = −μ̄_e/F: where an electrode exchanging electrons by that reaction would sit in equilibrium
 * with the local composition (the redox level of ESBD diagrams). With `standard: true`, every
 * species is taken at its reference concentration (its standard level μ° + zFφ), so it's the
 * standard level: for `SHE`, V°_e⁻(SHE) = φ + μ°_H⁺/F. `NaN` where a participant is absent.
 * @param {{ mu: Record<string, ArrayLike<number>>, muStd: Record<string, ArrayLike<number>> }} sol a solution
 * @param {{ equation: string, fixed?: Record<string, number> }} halfReaction from `half()`
 * @param {{ standard?: boolean }} [opts]
 * @returns {Float64Array}
 */
export function level(sol, halfReaction, { standard = false } = {}) {
  if (!isObject(halfReaction)) fail('level: give a half-reaction from half()');
  const { equation, fixed = {} } = halfReaction;
  const terms = parseEquation(equation, 'level');
  const source = standard ? sol.muStd : sol.mu;
  let nuE = 0, fixedSum = 0;
  const parts = [];
  for (const t of terms) {
    if (t.name === 'e-') nuE += t.nu;
    else if (fixed[t.name] !== undefined) fixedSum += t.nu * fixed[t.name];
    else if (source[t.name]) parts.push({ nu: t.nu, mu: source[t.name] });
    else fail(`level: '${t.name}' is neither a species of the solution nor in the half-reaction's fixed`);
  }
  if (nuE === 0) fail(`level: ${JSON.stringify(equation)} has no net electrons`);
  const n = sol.mu[Object.keys(sol.mu)[0]].length;
  const V = new Float64Array(n);
  for (let g = 0; g < n; g++) {
    // Σ ν μ̄ = 0 at equilibrium, so μ̄_e = −(Σ_others ν μ̄)/ν_e.
    let s = fixedSum;
    for (const p of parts) s += p.nu * p.mu[g];
    V[g] = s / (nuE * FARADAY);
  }
  return V;
}
