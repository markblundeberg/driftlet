// A device as a stack, for driftlet/kit: the layers left to right with the faces between them
// and a contact at each end, as a device is drawn. `build()` turns it into the plain definition
// (species, materials, regions, interfaces, contacts) that `new Device` takes. Nothing is added
// that you didn't write: it only moves things to where the definition keeps them, and converts
// the kit's shorthands (doping as donors and acceptors, reactions as equations).

import { DeviceError } from './device.js';
import { FARADAY } from './constants.js';
import { parseEquation, stoichiometry } from './equation.js';

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
 */
export function bath(c, reference, d = 0) {
  return { ...drive(d, 'bath'), bath: { c: { ...c }, reference } };
}

// Which species a material holds: its listed species, or a conductor's carrier.
function holds(mat, name) {
  if (!isObject(mat)) return false;
  if (mat.conductor) return mat.conductor.species === name;
  return isObject(mat.species) && mat.species[name] !== undefined;
}

// A face reaction written as an equation: each participant goes to the side that holds it. A
// species held on both sides is labelled with its side, 'Li+(left)', or its material's name,
// 'Li+(graphite)'; any other participant must be in `fixed`.
function faceReaction(rx, path, left, right) {
  const { equation, fixed, ...rest } = rx;
  const terms = parseEquation(equation, `${path}.equation`);
  const sides = { left: {}, right: {} };
  const put = (side, name, nu) => (sides[side][name] = (sides[side][name] ?? 0) + nu);
  // A fixed participant's side doesn't matter to the rate; it's written beside a conductor (the
  // solid metal of a plating reaction), else on the left.
  const fixedSide = right.mat?.conductor && !left.mat?.conductor ? 'right' : 'left';
  for (const t of terms) {
    if (fixed?.[t.name] !== undefined) {
      put(fixedSide, t.name, t.nu);
      continue;
    }
    const inL = holds(left.mat, t.name), inR = holds(right.mat, t.name);
    if (inL && inR) {
      fail(
        `${path}.equation: '${t.name}' is in both '${left.name}' and '${right.name}', so say which side: ` +
          `'${t.name}(left)' or '${t.name}(right)'` + (left.name !== right.name ? `, or '${t.name}(${left.name})', '${t.name}(${right.name})'` : ''),
      );
    }
    if (inL || inR) {
      put(inL ? 'left' : 'right', t.name, t.nu);
      continue;
    }
    const m = /^(.+)\(([^()]+)\)$/.exec(t.name);
    if (m) {
      const [, base, label] = m;
      const side = label === 'left' || (label === left.name && label !== right.name) ? 'left' : label === 'right' || (label === right.name && label !== left.name) ? 'right' : null;
      if (side) {
        const at = side === 'left' ? left : right;
        if (!holds(at.mat, base)) fail(`${path}.equation: '${base}' is absent from '${at.name}' (${side} of the face)`);
        put(side, base, t.nu);
        continue;
      }
    }
    fail(`${path}.equation: '${t.name}' is in neither '${left.name}' nor '${right.name}'; if it's a fixed-activity neutral, give its μ in fixed`);
  }
  const out = { ...rest };
  if (Object.keys(sides.left).length) out.left = sides.left;
  if (Object.keys(sides.right).length) out.right = sides.right;
  if (fixed !== undefined) out.fixed = fixed;
  return out;
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
 * - Reactions, in faces and in `bulkReactions`, may be written as an `equation` in place of
 *   their stoichiometry (see `parseEquation`).
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
    if (regions.length > 0) interfaces.push(pendingFace ? face(pendingFace.item, at(pendingFace.k), regions[regions.length - 1], item) : {});
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
  function face(item, path, l, r) {
    if (item.reactions === undefined) return item;
    if (!Array.isArray(item.reactions)) fail(`${path}.reactions must be a list`);
    const left = { name: l.material, mat: materials[l.material] };
    const right = { name: r.material, mat: materials[r.material] };
    return {
      ...item,
      reactions: item.reactions.map((rx, j) => (isObject(rx) && rx.equation !== undefined ? faceReaction(rx, `${path}.reactions[${j}]`, left, right) : rx)),
    };
  }

  const out = {};
  if (def.T !== undefined) out.T = def.T;
  out.species = species;
  out.materials = materials;
  out.regions = regions;
  if (interfaces.some((f) => Object.keys(f).length > 0)) out.interfaces = interfaces;
  out.contacts = { left: stack[0], right: stack[stack.length - 1] };
  if (def.bulkReactions !== undefined) {
    if (!Array.isArray(def.bulkReactions)) fail('build.bulkReactions must be a list');
    out.bulkReactions = def.bulkReactions.map((rx, j) => {
      if (!isObject(rx) || rx.equation === undefined) return rx;
      const { equation, ...rest } = rx;
      return { nu: stoichiometry(equation, `bulkReactions[${j}].equation`), ...rest };
    });
  }
  if (def.ports !== undefined) out.ports = def.ports;
  if (def.grid !== undefined) out.grid = def.grid;
  return out;
}
