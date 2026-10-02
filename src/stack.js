// A device as a stack, for driftlet/kit: the layers left to right with the faces between them
// and a contact at each end, as a device is drawn. `build()` turns it into the plain definition
// (species, materials, regions, interfaces, contacts) that `new Device` takes. Nothing is added
// that you didn't write: it only moves things to where the definition keeps them, and converts
// the kit's doping shorthand (donors and acceptors).

import { DeviceError } from './errors.js';
import { FARADAY } from './constants.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const isObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);

function fields(obj, path, allowed) {
  for (const k of Object.keys(obj)) if (obj[k] !== undefined && !allowed.includes(k)) fail(`${path}.${k}: unknown field (expected one of ${allowed.join(', ')})`);
}

const same = (a, b) => {
  if (a === b) return true;
  if (typeof a !== 'object' || typeof b !== 'object' || a === null || b === null || Array.isArray(a) !== Array.isArray(b)) return false;
  const ka = Object.keys(a), kb = Object.keys(b);
  return ka.length === kb.length && ka.every((k) => same(a[k], b[k]));
};

/**
 * Pieces of a definition (species and materials, as in a library entry) merged into one. A
 * species or material given twice must be given identically.
 * @param {...{ species?: object[], materials?: Record<string, object> }} parts
 * @returns {{ species: object[], materials: Record<string, object> }}
 */
export function combine(...parts) {
  const species = [], materials = {};
  parts.forEach((part, k) => {
    if (!isObject(part)) fail(`combine: part ${k} must be { species, materials }`);
    for (const sp of part.species ?? []) {
      const prev = species.find((s) => s.name === sp.name);
      if (!prev) species.push({ ...sp });
      else if (!same(prev, sp)) fail(`combine: species '${sp.name}' is given twice, differently (${JSON.stringify(prev)} and ${JSON.stringify(sp)})`);
    }
    for (const [name, mat] of Object.entries(part.materials ?? {})) {
      if (materials[name] === undefined) materials[name] = mat;
      else if (!same(materials[name], mat)) fail(`combine: material '${name}' is given twice, differently`);
    }
  });
  return { species, materials };
}

/**
 * A layer of the stack: a region of `material`, `length` m long, with any other region fields
 * (name, c0, grid, …) and the kit's doping shorthand `donors`, `acceptors` (mol/m³).
 * @param {string} material
 * @param {number} length m
 * @param {object} [more]
 */
export function layer(material, length, more = {}) {
  return { material, length, ...more };
}

// A drive: a number is a held voltage, a waveform { t, values } a voltage waveform, and an
// object with V, I or R is the drive itself.
function drive(d, path) {
  if (d === undefined) return { V: 0 };
  if (typeof d === 'number') return { V: d };
  if (isObject(d) && Array.isArray(d.t)) return { V: d };
  if (isObject(d)) {
    fields(d, path, ['V', 'I', 'R']);
    return { ...d };
  }
  return fail(`${path} must be a voltage, a waveform, or { V }, { I }, { V, R }`);
}

/**
 * An ohmic contact: each carrier in equilibrium with a metal at the terminal voltage, the
 * first being the terminal species, and every other sharing its voltage (offset 0, as holes and
 * electrons do at a metal). The default, both e⁻ and h⁺, gives n·p = n_i² at the contact
 * (infinite recombination); one carrier alone is a selective contact, and `['e-']` is a
 * current collector on a metal region.
 * @param {number | object} [drive] a voltage, a waveform, or { V }, { I }, { V, R }
 * @param {string[]} [carriers] default ['e-', 'h+']
 */
export function ohmic(d = 0, carriers = ['e-', 'h+']) {
  if (!Array.isArray(carriers) || carriers.length === 0) fail('ohmic: carriers must be a non-empty list of species names');
  const species = {};
  carriers.forEach((name, k) => (species[name] = k === 0 ? 'equilibrium' : { type: 'equilibrium', offset: 0 }));
  return { ...drive(d, 'ohmic'), terminal: carriers[0], species, phi: 'bulk' };
}

/**
 * A bath: a neutral composition `c` (mol/m³) held at the contact, its levels referred to the
 * charged `reference` species (the ion a reference electrode would sense, e.g. Cl⁻ for
 * Ag/AgCl), which is the terminal.
 * @param {Record<string, number>} c
 * @param {string} reference
 * @param {number | object} [drive]
 * @param {{ offset?: number }} [opts] where the reference species sits relative to the terminal
 *   voltage, V_ref = V + offset (V)
 */
export function bath(c, reference, d = 0, { offset } = {}) {
  return { ...drive(d, 'bath'), bath: { c: { ...c }, reference, ...(offset === undefined ? {} : { offset }) } };
}

const STACK_FIELDS = ['T', 'library', 'species', 'materials', 'stack', 'bulkReactions', 'ports', 'grid'];

/**
 * Build a plain device definition from a stack:
 *
 *     build({ library: [Si], stack: [ohmic(0), layer('Si', 1e-6, { donors }), layer(...), ohmic(0.5)] })
 *
 * - `stack`: a contact, then layers (regions: objects with a `material`) and faces (anything
 *   else: an interface's fields) left to right, then a contact. Two adjacent layers with no face
 *   between them get the default face. Nested lists are flattened.
 * - Layers may give doping as `donors` and `acceptors` (mol/m³) in place of `fixedCharge`.
 * - `library`: pieces with species and materials, merged with `species` and `materials`.
 * - `T`, `grid`, `ports` pass through.
 * @param {object} def
 * @returns {object} a device definition for `new Device`
 */
export function build(def) {
  if (!isObject(def)) fail('build: give { library, stack, … }');
  fields(def, 'build', STACK_FIELDS);
  if (def.library !== undefined && !Array.isArray(def.library)) fail('build.library must be a list of { species, materials }');
  const { species, materials } = combine(...(def.library ?? []), { species: def.species ?? [], materials: def.materials ?? {} });

  if (!Array.isArray(def.stack)) fail('build.stack must be a list: a contact, layers and faces, a contact');
  const stack = def.stack.flat(Infinity);
  if (stack.length < 3) fail('build.stack needs a contact, at least one layer, and a contact');
  for (const k of [0, stack.length - 1]) {
    if (!isObject(stack[k]) || stack[k].material !== undefined) fail(`stack[${k}]: the stack starts and ends with a contact`);
  }
  const at = (k) => `stack[${k}]`;
  const regions = [], interfaces = [];
  let pendingFace = null; // the face written since the last layer
  for (let k = 1; k < stack.length - 1; k++) {
    const item = stack[k];
    if (!isObject(item)) fail(`${at(k)} must be a layer (with a material) or a face (an interface's fields)`);
    if (item.material === undefined) {
      if (regions.length === 0) fail(`${at(k)}: a face before the first layer (the contact is the boundary there)`);
      if (pendingFace) fail(`${at(k)}: two faces in a row, with no layer between them`);
      pendingFace = { item, k };
      continue;
    }
    if (regions.length > 0) interfaces.push(pendingFace ? pendingFace.item : {});
    pendingFace = null;
    regions.push(region(item, at(k)));
  }
  if (pendingFace) fail(`${at(pendingFace.k)}: a face after the last layer (the contact is the boundary there)`);
  if (regions.length === 0) fail('build.stack has no layers');

  function region(item, path) {
    const { donors, acceptors, ...r } = item;
    if (donors !== undefined || acceptors !== undefined) {
      if (r.fixedCharge !== undefined) fail(`${path}: give doping as donors and acceptors, or fixedCharge, not both`);
      for (const [k, v] of [['donors', donors], ['acceptors', acceptors]]) {
        if (v !== undefined && !(Number.isFinite(v) && v >= 0)) fail(`${path}.${k} must be a concentration ≥ 0 (mol/m³)`);
      }
      r.fixedCharge = ((donors ?? 0) - (acceptors ?? 0)) * FARADAY;
    }
    return r;
  }
  const out = {};
  if (def.T !== undefined) out.T = def.T;
  out.species = species;
  out.materials = materials;
  out.regions = regions;
  if (interfaces.some((f) => Object.keys(f).length > 0)) out.interfaces = interfaces;
  out.contacts = { left: stack[0], right: stack[stack.length - 1] };
  if (def.bulkReactions !== undefined) out.bulkReactions = def.bulkReactions;
  if (def.ports !== undefined) out.ports = def.ports;
  if (def.grid !== undefined) out.grid = def.grid;
  return out;
}
