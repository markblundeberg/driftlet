// A device as a stack, for driftlet/kit: the layers left to right with the faces between them
// and a contact at each end, as a device is drawn. `build()` turns it into the plain definition
// (species, materials, regions, interfaces, contacts) that `new Device` takes. Nothing is added
// that you didn't write: it only moves things to where the definition keeps them, and converts
// the kit's doping shorthand (donors and acceptors).

import { DeviceError, checkOptions, shown } from './errors.js';
import { FARADAY } from './constants.js';
import { parseEquation } from './equation.js';

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
    // (trap()'s whole result, given for its library, dropped the traps without a word)
    if (part.library !== undefined) fail(`combine: part ${k} is { library, … }, as trap() returns: give its library here (trap(…).library), and its bulkReactions and c0 to the device`);
    fields(part, `combine: part ${k}`, ['species', 'materials']);
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
export function drive(d, path) {
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
 * A bath: a neutral composition `c` (mol/m³) held at the contact, read in one of three ways.
 * - By a reference electrode, given as its half-reaction (`AgAgCl`, or any `half()` exchanging
 *   one of the bath's ions, the rest of it `fixed`): the terminal voltage is that electrode's
 *   own, its electrons' level in equilibrium with the bath, as its wire would read.
 * - By a reference species, a charged one: the terminal voltage is that species' level. That is
 *   an electrode reversible to it whose other participants are elements in their standard states
 *   (`'Cl-'` is a chlorine electrode, 1.36 V above SHE at 1 M; an Ag/AgCl wire sits 1.138 V
 *   below it), which is all the same between two alike.
 * - By nothing, `bath(c, drive, opts)`: the terminal voltage is the bath's φ, an ideal salt
 *   bridge, the convention of membrane potentials.
 * @param {Record<string, number>} c
 * @param {string | { equation: string, fixed?: Record<string, number> } | number | object} [reference]
 *   the reference electrode's half-reaction, the reference species, or the drive
 * @param {number | object} [drive]
 * @param {{ offset?: number, material?: string }} [opts] `offset`: where the reference species
 *   sits relative to the terminal voltage, V_ref = V + offset (V), added to a reference
 *   electrode's own. `material`: the bath's own material (a solution against an ion exchanger),
 *   if not the end layer's
 */
export function bath(c, reference, d = 0, opts = {}) {
  if (!isObject(c)) fail(`bath: the composition comes first, as { species: mol/m³ } (got ${JSON.stringify(c)})`);
  // A reference electrode: Σν·μ̄ = 0 over its half-reaction puts the one ion it exchanges at
  // V − Σν·μ(fixed)/(ν_e·F) when its electrons are at V (the equation balanced as written:
  // ν_e = ν_ion·z_ion). AgCl(s) + e⁻ ⇌ Ag(s) + Cl⁻: Cl⁻ sits −μ°(AgCl)/F = 1.138 V above the wire.
  let electrode = 0;
  if (isObject(reference) && typeof reference.equation === 'string') {
    const terms = parseEquation(reference.equation, 'bath: the reference electrode');
    const fixed = isObject(reference.fixed) ? reference.fixed : {};
    const electrons = terms.filter((t) => t.name === 'e-');
    const ions = terms.filter((t) => t.name !== 'e-' && fixed[t.name] === undefined);
    if (electrons.length !== 1) fail(`bath: the reference electrode ${JSON.stringify(reference.equation)} needs electrons, 'e-', on exactly one side`);
    if (ions.length !== 1 || c[ions[0].name] === undefined) {
      fail(
        `bath: the reference electrode ${JSON.stringify(reference.equation)} must exchange exactly one of the bath's species (${Object.keys(c).join(', ')}); ` +
          `it has ${ions.length === 0 ? 'none' : ions.map((t) => `'${t.name}'`).join(' and ')} besides its electrons and its fixed participants (give the others' μ in half()'s fixed)`,
      );
    }
    let sum = 0;
    for (const t of terms) if (fixed[t.name] !== undefined) sum += t.nu * fixed[t.name];
    electrode = -sum / (electrons[0].nu * FARADAY);
    reference = ions[0].name;
  }
  if (typeof reference !== 'string') {
    // bath(c, drive, opts): a third argument that isn't the options was put in the wrong place
    // (bath(c, { I: 0 }, 'Cl-') once built a bath without its reference, silently).
    if (d !== 0 && d !== undefined && !isObject(d)) fail(`bath: ${JSON.stringify(d)} after the drive: the reference species comes second, bath(c, ${JSON.stringify(d)}, drive), and the options last`);
    [d, opts] = [reference ?? 0, isObject(d) ? d : {}];
  }
  checkOptions('bath', opts, ['offset', 'material']);
  const { material } = opts;
  if (electrode !== 0 && opts.offset !== undefined && !Number.isFinite(opts.offset)) fail(`bath: offset must be a number of volts, got ${shown(opts.offset)}`);
  const offset = electrode !== 0 ? electrode + (opts.offset ?? 0) : opts.offset;
  return {
    ...drive(d, 'bath'),
    bath: {
      c: { ...c },
      ...(typeof reference === 'string' ? { reference } : {}),
      ...(offset === undefined ? {} : { offset }),
      ...(material === undefined ? {} : { material }),
    },
  };
}

const STACK_FIELDS = ['T', 'library', 'species', 'materials', 'stack', 'bulkReactions', 'ports', 'grid', 'geometry'];

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
 * - `T`, `grid`, `ports`, `geometry` pass through.
 * @param {object} def
 * @returns {object} a device definition for `new Device`
 */
export function build(def) {
  if (!isObject(def)) fail('build: give { library, stack, … }');
  fields(def, 'build', STACK_FIELDS);
  if (def.T !== undefined && !(typeof def.T === 'number' && Number.isFinite(def.T) && def.T > 0)) fail(`build.T must be the temperature in K, a positive number, got ${shown(def.T)}`);
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
  // A contact holding a carrier its layer lacks: most often ohmic()'s default on a metal.
  for (const [k, end] of [[0, regions[0]], [stack.length - 1, regions[regions.length - 1]]]) {
    const held = isObject(stack[k].species) ? Object.keys(stack[k].species).filter((s) => stack[k].species[s] !== 'blocked' && stack[k].species[s]?.type !== 'blocked') : [];
    const mat = materials[end.material];
    if (!isObject(mat)) continue;
    const has = (s) => (isObject(mat.conductor) ? mat.conductor.species === s : isObject(mat.species) && mat.species[s] !== undefined);
    const lacks = held.filter((s) => !has(s));
    if (lacks.length === 0) continue;
    const kept = held.filter(has);
    const hint = stack[k].phi === 'bulk' && kept.length > 0 ? `: write ohmic(V, [${kept.map((s) => `'${s}'`).join(', ')}]) (ohmic() holds e- and h+ unless told${isObject(mat.conductor) ? '; a metal has only its carrier' : ''})` : '';
    fail(`${at(k)}: the contact holds ${lacks.map((s) => `'${s}'`).join(' and ')}, which its layer's material '${end.material}' doesn't have${hint}`);
  }

  function region(item, path) {
    checkOptions(path, item, ['name', 'material', 'length', 'fixedCharge', 'c0', 'velocity', 'mixing', 'grid', 'grade', 'donors', 'acceptors']);
    const { donors, acceptors, ...r } = item;
    if (donors !== undefined || acceptors !== undefined) {
      if (r.fixedCharge !== undefined) fail(`${path}: give doping as donors and acceptors, or fixedCharge, not both`);
      for (const [k, v] of [['donors', donors], ['acceptors', acceptors]]) {
        if (v !== undefined && !(Number.isFinite(v) && v >= 0)) fail(`${path}.${k} must be a concentration ≥ 0 (mol/m³)`);
        // Over 1e5 mol/m³ (6e22 cm⁻³, more than a crystal has atoms) it's a count per cm³ given as
        // mol/m³; describe() caught it, but named fixedCharge, which the user hadn't written.
        if (v > 1e5) fail(`${path}.${k}: ${v.toExponential()} mol/m³ is ${(v * 6.02214076e17).toPrecision(2)} cm⁻³, more than a crystal has atoms; ${k} is mol/m³: for ${v.toExponential()} cm⁻³, units.perCm3(${v.toExponential()})`);
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
  if (def.geometry !== undefined) out.geometry = def.geometry;
  return out;
}
