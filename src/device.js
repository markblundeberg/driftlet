// Device definition: validation and normalisation.
//
// A definition is a plain, serialisable object (so it can be posted to a Worker). This module
// checks it thoroughly and turns it into an internal model with index-based lookups. Every
// problem is a DeviceError whose message names the offending path. There are no silent
// defaults for physically meaningful choices: interface alignments, reference concentrations
// and standard potentials must all be given.

import { buildGrid } from './grid.js';
import { FARADAY, GAS_CONSTANT } from './constants.js';

export class DeviceError extends Error {
  constructor(message) {
    super(message);
    this.name = 'DeviceError';
  }
}

const SPECIES_LINK_TYPES = new Set(['blocked', 'fixed', 'conductance', 'kinetic']);
const PHI_LINK_TYPES = new Set(['free', 'neutral', 'capacitive']);

const isObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);
const isFiniteNumber = (v) => typeof v === 'number' && Number.isFinite(v);

function need(cond, message) {
  if (!cond) throw new DeviceError(message);
}

function finite(v, path) {
  need(isFiniteNumber(v), `${path} must be a finite number, got ${JSON.stringify(v)}`);
  return v;
}

function positive(v, path) {
  need(isFiniteNumber(v) && v > 0, `${path} must be a positive number, got ${JSON.stringify(v)}`);
  return v;
}

function nonNegative(v, path) {
  need(isFiniteNumber(v) && v >= 0, `${path} must be a non-negative number, got ${JSON.stringify(v)}`);
  return v;
}

/**
 * Validate and normalise a device definition.
 * @param {object} def
 * @returns {object} internal model
 */
export function normalizeDevice(def) {
  need(isObject(def), 'device definition must be an object');

  const T = def.T === undefined ? 298.15 : positive(def.T, 'T');
  const RT = GAS_CONSTANT * T;

  // --- species
  need(Array.isArray(def.species) && def.species.length > 0, 'species must be a non-empty array');
  const species = [];
  const speciesIndex = new Map();
  def.species.forEach((sp, i) => {
    const path = `species[${i}]`;
    need(isObject(sp), `${path} must be an object`);
    need(typeof sp.name === 'string' && sp.name.length > 0, `${path}.name must be a non-empty string`);
    need(!speciesIndex.has(sp.name), `${path}.name: duplicate species '${sp.name}'`);
    need(Number.isInteger(sp.z), `${path}.z (charge number) must be an integer, got ${JSON.stringify(sp.z)}`);
    const cRef = sp.cRef === undefined ? undefined : positive(sp.cRef, `${path}.cRef`);
    speciesIndex.set(sp.name, i);
    species.push({ name: sp.name, z: sp.z, cRef });
  });
  const nSpecies = species.length;

  // --- materials
  need(isObject(def.materials) && Object.keys(def.materials).length > 0, 'materials must be a non-empty object');
  const materials = [];
  const materialIndex = new Map();
  for (const [mname, mat] of Object.entries(def.materials)) {
    const path = `materials.${mname}`;
    need(isObject(mat), `${path} must be an object`);
    const epsr = positive(mat.epsr, `${path}.epsr`);
    need(isObject(mat.species), `${path}.species must be an object mapping species names to parameters`);
    const present = new Uint8Array(nSpecies);
    const D = new Float64Array(nSpecies);
    const mu0 = new Float64Array(nSpecies);
    const cRef = new Float64Array(nSpecies);
    for (const [sname, p] of Object.entries(mat.species)) {
      const spath = `${path}.species.${sname}`;
      need(speciesIndex.has(sname), `${spath}: unknown species '${sname}'`);
      need(isObject(p), `${spath} must be an object`);
      const i = speciesIndex.get(sname);
      present[i] = 1;
      D[i] = nonNegative(p.D, `${spath}.D`);
      mu0[i] = finite(p.mu0, `${spath}.mu0`);
      const c = p.cRef ?? species[i].cRef;
      need(c !== undefined, `${spath}.cRef: no reference concentration (give it here or on the species)`);
      cRef[i] = positive(c, `${spath}.cRef`);
    }
    materialIndex.set(mname, materials.length);
    materials.push({ name: mname, epsr, present, D, mu0, cRef });
  }
  species.forEach((sp, i) => {
    need(materials.some((m) => m.present[i]), `species '${sp.name}' is not present in any material`);
  });

  // --- regions
  need(Array.isArray(def.regions) && def.regions.length > 0, 'regions must be a non-empty array');
  const regions = def.regions.map((reg, r) => {
    const path = `regions[${r}]`;
    need(isObject(reg), `${path} must be an object`);
    need(materialIndex.has(reg.material), `${path}.material: unknown material ${JSON.stringify(reg.material)}`);
    const length = positive(reg.length, `${path}.length`);
    const fixedCharge = reg.fixedCharge === undefined ? 0 : finite(reg.fixedCharge, `${path}.fixedCharge`);
    if (reg.grid !== undefined) need(isObject(reg.grid), `${path}.grid must be an object`);
    const mat = materials[materialIndex.get(reg.material)];
    // Initial composition: the starting state, and the conserved amount of any spectator.
    const c0 = new Float64Array(nSpecies).fill(NaN);
    if (reg.c0 !== undefined) {
      need(isObject(reg.c0), `${path}.c0 must be an object mapping species names to concentrations`);
      for (const [sname, v] of Object.entries(reg.c0)) {
        need(speciesIndex.has(sname), `${path}.c0.${sname}: unknown species '${sname}'`);
        const i = speciesIndex.get(sname);
        need(mat.present[i], `${path}.c0.${sname}: '${sname}' is absent from material '${mat.name}'`);
        c0[i] = positive(v, `${path}.c0.${sname}`);
      }
    }
    return {
      name: reg.name ?? `region ${r}`,
      material: materialIndex.get(reg.material),
      length,
      fixedCharge,
      c0,
      grid: reg.grid,
    };
  });

  // --- interfaces (one per face between consecutive regions)
  const nFaces = regions.length - 1;
  const idefs = def.interfaces ?? [];
  need(Array.isArray(idefs), 'interfaces must be an array (one entry per face between regions)');
  need(idefs.length <= nFaces, `interfaces has ${idefs.length} entries but there are only ${nFaces} faces between regions`);
  const interfaces = [];
  for (let f = 0; f < nFaces; f++) {
    interfaces.push(normalizeInterface(idefs[f], f, regions, materials, species, speciesIndex));
  }

  // --- contacts
  const cdefs = def.contacts ?? {};
  need(isObject(cdefs), 'contacts must be an object with optional left and right');
  for (const k of Object.keys(cdefs)) need(k === 'left' || k === 'right', `contacts.${k}: expected only 'left' and 'right'`);
  const contacts = {
    left: normalizeContact(cdefs.left, 'left', regions[0], materials, species, speciesIndex, RT),
    right: normalizeContact(cdefs.right, 'right', regions[regions.length - 1], materials, species, speciesIndex, RT),
  };

  // An electrostatic anchor is needed, or φ (and every level with it) floats.
  const anchored = ['left', 'right'].some(
    (side) => contacts[side].phi.type === 'capacitive' || contacts[side].species.some((l) => l.type !== 'blocked'),
  );
  need(
    anchored,
    'device has no electrostatic anchor: every species is blocked at both ends and no contact has a ' +
      'capacitive (gate) link, so the potential is undetermined. Add a gate or connect a species.',
  );

  // --- grid
  if (def.grid !== undefined) need(isObject(def.grid), 'grid must be an object');
  let grid;
  try {
    grid = buildGrid(regions, def.grid ?? {});
  } catch (err) {
    throw new DeviceError(`grid: ${err.message}`);
  }

  return { T, RT, F: FARADAY, species, speciesIndex, materials, materialIndex, regions, interfaces, contacts, grid };
}

function normalizeInterface(idef, f, regions, materials, species, speciesIndex) {
  const left = regions[f], right = regions[f + 1];
  const matL = materials[left.material], matR = materials[right.material];
  const where = `interfaces[${f}] (between ${left.name} [${matL.name}] and ${right.name} [${matR.name}])`;
  const same = left.material === right.material;

  if (idef === undefined || idef === null) {
    need(
      same,
      `${where}: an interface between different materials needs an alignment ` +
        `({ dipole } or { step: { species, value } }). There is no default (no Anderson or Schottky–Mott rule).`,
    );
    return { dipole: 0, sheetCharge: 0 };
  }
  need(isObject(idef), `${where} must be an object`);
  const given = ['dipole', 'step', 'reaction'].filter((k) => idef[k] !== undefined);
  need(given.length <= 1, `${where}: give exactly one alignment, got ${given.join(' and ')}`);
  need(
    given.length === 1 || same,
    `${where}: an interface between different materials needs an alignment ` +
      `({ dipole } or { step: { species, value } }). There is no default (no Anderson or Schottky–Mott rule).`,
  );

  let dipole = 0;
  if (idef.dipole !== undefined) {
    dipole = finite(idef.dipole, `${where}.dipole`);
  } else if (idef.step !== undefined) {
    const st = idef.step;
    need(isObject(st), `${where}.step must be { species, value }`);
    need(speciesIndex.has(st.species), `${where}.step.species: unknown species ${JSON.stringify(st.species)}`);
    const i = speciesIndex.get(st.species);
    const z = species[i].z;
    need(z !== 0, `${where}.step.species: '${st.species}' is neutral; alignment needs a charged species`);
    need(matL.present[i] && matR.present[i], `${where}.step.species: '${st.species}' must be present on both sides`);
    const value = finite(st.value, `${where}.step.value`);
    // value = (μ°_R + zFφ_R) − (μ°_L + zFφ_L), so φ_R − φ_L = (value − Δμ°) / (zF)
    dipole = (value - (matR.mu0[i] - matL.mu0[i])) / (z * FARADAY);
  } else if (idef.reaction !== undefined) {
    throw new DeviceError(`${where}.reaction: reaction-based alignment is not supported yet`);
  }
  const sheetCharge = idef.sheetCharge === undefined ? 0 : finite(idef.sheetCharge, `${where}.sheetCharge`);
  return { dipole, sheetCharge };
}

// A contact has a terminal voltage V (set by the circuit; 0 by default) and a link for every
// species and for φ. A fixed charged species sits at V_i = V + offset_i, i.e.
// μ̄_i = z_i F (V + offset_i); the offset is an interface property (0 for the terminal species,
// e.g. e⁻ at a metal; E° for an ion at a reversible electrode) and is never defaulted for
// the others. A fixed neutral species takes an absolute μ̄. A bath computes the offsets from
// a composition, anchored through its reference species.
function normalizeContact(cdef, side, region, materials, species, speciesIndex, RT) {
  const path = `contacts.${side}`;
  const mat = materials[region.material];
  const links = species.map(() => ({ type: 'blocked' }));
  let phi = { type: 'free' };
  let terminal = null;
  if (cdef === undefined || cdef === null) return { V: 0, species: links, phi, terminal };
  need(isObject(cdef), `${path} must be an object`);
  const V = cdef.V === undefined ? 0 : finite(cdef.V, `${path}.V`);

  if (cdef.terminal !== undefined) {
    need(speciesIndex.has(cdef.terminal), `${path}.terminal: unknown species ${JSON.stringify(cdef.terminal)}`);
    terminal = speciesIndex.get(cdef.terminal);
  }

  if (cdef.species !== undefined) {
    need(isObject(cdef.species), `${path}.species must be an object mapping species names to links`);
    for (const [sname, raw] of Object.entries(cdef.species)) {
      const lpath = `${path}.species.${sname}`;
      need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'`);
      const i = speciesIndex.get(sname);
      const link = typeof raw === 'string' ? { type: raw } : raw;
      need(isObject(link), `${lpath} must be a link type string or an object with a type`);
      need(SPECIES_LINK_TYPES.has(link.type), `${lpath}.type must be one of ${[...SPECIES_LINK_TYPES].join(', ')}`);
      if (link.type !== 'blocked') {
        need(mat.present[i], `${lpath}: '${sname}' is absent from the end material '${mat.name}', so it can only be blocked`);
      }
      if (link.type === 'conductance') positive(link.G, `${lpath}.G`);
      if (link.type === 'fixed') {
        if (species[i].z === 0) {
          links[i] = { type: 'fixed', mu: finite(link.mu, `${lpath}.mu (a neutral species is fixed by its μ̄, J/mol)`) };
          continue;
        }
        need(link.mu === undefined, `${lpath}: a charged species is fixed by an offset from the terminal voltage, not by mu`);
        const offset = link.offset ?? (terminal === i ? 0 : undefined);
        need(
          offset !== undefined,
          `${lpath}.offset: give V_i − V_terminal (V) for this species; only the terminal species defaults to 0`,
        );
        links[i] = { type: 'fixed', offset: finite(offset, `${lpath}.offset`) };
        continue;
      }
      links[i] = { ...link };
    }
  }

  if (cdef.bath !== undefined) {
    need(cdef.species === undefined, `${path}: give either bath or species links, not both`);
    const bath = cdef.bath;
    need(isObject(bath) && isObject(bath.c), `${path}.bath must be { c: { species: concentration }, reference }`);
    need(speciesIndex.has(bath.reference), `${path}.bath.reference: unknown species ${JSON.stringify(bath.reference)}`);
    const r = speciesIndex.get(bath.reference);
    need(species[r].z !== 0, `${path}.bath.reference: the reference species must be charged`);
    need(bath.c[bath.reference] !== undefined, `${path}.bath.reference: '${bath.reference}' must be in the bath`);
    need(cdef.terminal === undefined || terminal === r, `${path}.terminal must be the bath's reference species`);
    terminal = r;
    let charge = region.fixedCharge / FARADAY, scale = Math.abs(charge);
    const cb = new Float64Array(species.length);
    for (const [sname, v] of Object.entries(bath.c)) {
      need(speciesIndex.has(sname), `${path}.bath.c.${sname}: unknown species '${sname}'`);
      const i = speciesIndex.get(sname);
      need(mat.present[i], `${path}.bath.c.${sname}: '${sname}' is absent from the end material '${mat.name}'`);
      cb[i] = positive(v, `${path}.bath.c.${sname}`);
      charge += species[i].z * cb[i];
      scale += Math.abs(species[i].z) * cb[i];
    }
    need(Math.abs(charge) <= 1e-9 * scale, `${path}.bath: composition is not neutral (net ${charge} mol/m³ of charge)`);
    // The reference species pins the bath's φ; every other species follows from composition.
    const level = (i) => mat.mu0[i] + RT * Math.log(cb[i] / mat.cRef[i]); // μ̄ − zFφ_bath
    const refOffset = bath.offset === undefined ? 0 : finite(bath.offset, `${path}.bath.offset`);
    const beta = refOffset - level(r) / (species[r].z * FARADAY); // φ_bath − V
    for (let i = 0; i < species.length; i++) {
      if (!(cb[i] > 0)) continue;
      links[i] =
        species[i].z === 0
          ? { type: 'fixed', mu: level(i) }
          : { type: 'fixed', offset: beta + level(i) / (species[i].z * FARADAY) };
    }
    phi = { type: 'neutral' };
  }

  if (cdef.phi !== undefined) {
    const raw = typeof cdef.phi === 'string' ? { type: cdef.phi } : cdef.phi;
    need(isObject(raw), `${path}.phi must be a link type string or an object with a type`);
    need(PHI_LINK_TYPES.has(raw.type), `${path}.phi.type must be one of ${[...PHI_LINK_TYPES].join(', ')}`);
    if (raw.type === 'capacitive') {
      positive(raw.C, `${path}.phi.C`);
      finite(raw.zeroCharge, `${path}.phi.zeroCharge (the zero-charge alignment: flat-band voltage or pzc)`);
      need(raw.V === undefined, `${path}.phi.V: the gate voltage is the contact's terminal voltage, ${path}.V`);
    }
    phi = { ...raw };
  } else if (cdef.bath === undefined && links.some((l) => l.type !== 'blocked')) {
    throw new DeviceError(
      `${path}.phi: a contact with connected species needs an explicit φ condition ('neutral', 'free', or capacitive)`,
    );
  }
  if (phi.type === 'neutral') {
    need(
      links.some((l, i) => l.type !== 'blocked' && species[i].z !== 0),
      `${path}.phi: a neutral link needs at least one connected charged species at this contact`,
    );
  }

  if (terminal !== null) {
    need(links[terminal].type !== 'blocked', `${path}.terminal: '${species[terminal].name}' is blocked at this contact`);
  }
  return { V, species: links, phi, terminal };
}
