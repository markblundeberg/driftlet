// driftlet/kit: porcelain. Helpers that read and write plain device-definition data, so a
// device can be set up from familiar quantities. Nothing here changes what a device means: each
// helper returns an ordinary spec value (a dipole, a zeroCharge, …) that you can inspect.

import { DeviceError } from './errors.js';
import { FARADAY } from './constants.js';

export { build, combine, layer, ohmic, bath } from './stack.js';
export { parseEquation, stoichiometry, half, SHE, level } from './equation.js';
export { traces, speciesRole, typeset } from './traces.js';
export { live } from './live.js';
export { pulse, square, triangle, ramp, injector, recombination, photogeneration } from './sources.js';
export { describe, unitWarnings } from './describe.js';
export { recorder } from './recorder.js';
export { henderson, planck } from './junction.js';
export { check } from './checks.js';
export { polarization } from './polarization.js';
export { perovskiteCell, PEROVSKITE_SCANS, hysteresis } from './devices.js';
export { IONS, H2O, WATER_EPSR, SEMICONDUCTORS, METALS, aqueous, semiconductor, metal } from './data.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const isObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);

/**
 * One material seen from vacuum: a level inside it (the anchor) and how far above that level,
 * in electron energy, the vacuum level just outside lies (the offset, V). So
 * V_vac = V_anchor − offset.
 * @typedef {object} VacuumSide
 * @property {string} material material name in the definition
 * @property {string} anchor a charged species present there (its standard level
 *   V°_i = φ + μ°_i/(z_i F): a band edge for e⁻ or h⁺), or 'phi' (the inner potential itself)
 * @property {number} offset V: an electron affinity, ionisation energy, Trasatti's 4.44 V, …
 */

/**
 * A material's vacuum level relative to its own φ, V_vac − φ, in volts.
 * @param {object} def device definition (its `species` and `materials` are read)
 * @param {VacuumSide} side
 * @param {string} [path] for error messages
 */
export function vacuumLevel(def, side, path = 'vacuum') {
  if (!isObject(side)) fail(`${path} must be { material, anchor, offset }`);
  const mat = def?.materials?.[side.material];
  if (!isObject(mat)) fail(`${path}.material: unknown material ${JSON.stringify(side.material)}`);
  if (mat.conductor !== undefined) fail(`${path}.material: '${side.material}' is a conductor, with no φ; use vacuumZeroCharge with its work function`);
  if (!Number.isFinite(side.offset)) fail(`${path}.offset must be a finite number (V)`);
  if (side.anchor === 'phi') return -side.offset;
  const sp = (def.species ?? []).find((s) => s.name === side.anchor);
  if (!sp) fail(`${path}.anchor must be 'phi' or a charged species, got ${JSON.stringify(side.anchor)}`);
  if (!(sp.z !== 0)) fail(`${path}.anchor: '${side.anchor}' is neutral, so it has no level in volts`);
  const entry = mat.species?.[side.anchor];
  if (!isObject(entry)) fail(`${path}.anchor: '${side.anchor}' is absent from '${side.material}'`);
  if (!Number.isFinite(entry.mu0)) fail(`${path}.anchor: '${side.anchor}' in '${side.material}' has no mu0`);
  return entry.mu0 / (sp.z * FARADAY) - side.offset;
}

/**
 * The vacuum-level (Anderson-type) estimate of a face's alignment: the two sides' vacuum levels
 * are taken to coincide. Returns the face's `dipole`, φ_R − φ_L in volts.
 * @param {object} def device definition
 * @param {VacuumSide} left
 * @param {VacuumSide} right
 */
export function vacuumDipole(def, left, right) {
  return vacuumLevel(def, left, 'vacuumDipole: left') - vacuumLevel(def, right, 'vacuumDipole: right');
}

/**
 * The vacuum-level (Schottky–Mott-type) estimate of a `zeroCharge` where a conductor meets a
 * material: at a conductor region's face, or a contact's capacitive or dipole law. The conductor's
 * vacuum level sits its work function W beyond its Fermi level, so
 * zeroCharge = V_F − φ_edge at zero charge = W + (V_vac − φ)_inside.
 * @param {object} def device definition
 * @param {number} workFunction V
 * @param {VacuumSide} inside the material across the face
 */
export function vacuumZeroCharge(def, workFunction, inside) {
  if (!Number.isFinite(workFunction)) fail('vacuumZeroCharge: workFunction must be a finite number (V)');
  return workFunction + vacuumLevel(def, inside, 'vacuumZeroCharge: inside');
}
