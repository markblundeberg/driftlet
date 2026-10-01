// Device definition: validation and normalisation.
//
// A definition is a plain, serialisable object (so it can be posted to a Worker). This module
// checks it thoroughly and turns it into an internal model with index-based lookups. Every
// problem is a DeviceError whose message names the offending path. There are no silent
// defaults for physically meaningful choices: interface alignments, reference concentrations
// and standard potentials must all be given.

import { buildGrid } from './grid.js';
import { FARADAY, GAS_CONSTANT } from './constants.js';
import { normalizeStatistics } from './statistics.js';

export class DeviceError extends Error {
  constructor(message) {
    super(message);
    this.name = 'DeviceError';
  }
}

// Contacts use the same laws as internal faces: the outside is a phase with known levels.
const SPECIES_LINK_TYPES = new Set(['blocked', 'equilibrium', 'conductance']);
const INTERFACE_LINK_TYPES = new Set(['equilibrium', 'blocked', 'conductance']);
const PHI_LINK_TYPES = new Set(['bulk', 'neutral', 'capacitive', 'dipole']);

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
  const warnings = [];
  need(isObject(def.materials) && Object.keys(def.materials).length > 0, 'materials must be a non-empty object');
  const materials = [];
  const materialIndex = new Map();
  for (const [mname, mat] of Object.entries(def.materials)) {
    const path = `materials.${mname}`;
    need(isObject(mat), `${path} must be an object`);
    if (mat.metal !== undefined) {
      materialIndex.set(mname, materials.length);
      materials.push(normalizeMetal(mat, mname, path, species, speciesIndex));
      continue;
    }
    // ε = 0 makes the material strictly neutral: Poisson becomes local neutrality there.
    const epsr = nonNegative(mat.epsr, `${path}.epsr`);
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
    const m = { name: mname, epsr, present, D, mu0, cRef };
    const st = normalizeStatistics({ ...m, statistics: mat.statistics }, path, species, speciesIndex, RT, { need, finite, positive });
    m.models = st.models;
    m.modelOf = st.modelOf;
    m.ideal = st.models.length === 0;
    m.phiFree = st.phiFree; // every charged species is in a neutral combination: φ is undefined
    warnings.push(...st.warnings);
    materialIndex.set(mname, materials.length);
    materials.push(m);
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
    // Imposed flow (m/s, toward +x) carrying every mobile species, and eddy mixing (m²/s).
    const velocity = reg.velocity === undefined ? 0 : finite(reg.velocity, `${path}.velocity`);
    const mixing = reg.mixing === undefined ? 0 : nonNegative(reg.mixing, `${path}.mixing`);
    if (reg.grid !== undefined) need(isObject(reg.grid), `${path}.grid must be an object`);
    const mat = materials[materialIndex.get(reg.material)];
    if (mat.metal) {
      need(fixedCharge === 0, `${path}.fixedCharge: a metal region is neutral in bulk (its carriers are the conduction electrons)`);
      need(reg.c0 === undefined, `${path}.c0: a metal region has no composition to give, only its Fermi level`);
    }
    if (mat.epsr === 0) {
      need(
        species.some((sp, i) => mat.present[i] && sp.z !== 0) || fixedCharge === 0,
        `${path}: material '${mat.name}' has ε = 0 (strictly neutral) but no mobile charged species to neutralise the fixed charge`,
      );
    }
    // An insertion host's fixed charge is balanced by its own background carriers.
    const host = mat.models.find((md) => md.type === 'insertion');
    const background = host ? -fixedCharge / (species[host.carrier].z * FARADAY) : 0;
    if (host) {
      need(
        background >= 0,
        `${path}.fixedCharge: in an insertion host the fixed charge is balanced by background ${species[host.carrier].name}, ` +
          'so it must have the opposite sign to that carrier',
      );
    }
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
      background,
      velocity,
      mixing,
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
    interfaces.push(normalizeInterface(idefs[f], f, regions, materials, species, speciesIndex, RT));
  }

  // --- bulk reactions
  need(def.bulkReactions === undefined || Array.isArray(def.bulkReactions), 'bulkReactions must be an array');
  const reactions = (def.bulkReactions ?? []).map((rdef, k) =>
    normalizeReaction(rdef, `bulkReactions[${k}]`, species, speciesIndex, materials, materialIndex, regions, RT),
  );

  // --- contacts
  const cdefs = def.contacts ?? {};
  need(isObject(cdefs), 'contacts must be an object with optional left and right');
  for (const k of Object.keys(cdefs)) need(k === 'left' || k === 'right', `contacts.${k}: expected only 'left' and 'right'`);
  const contacts = {
    left: normalizeContact(cdefs.left, 'left', regions[0], materials, species, speciesIndex, RT),
    right: normalizeContact(cdefs.right, 'right', regions[regions.length - 1], materials, species, speciesIndex, RT),
  };

  // Every electrostatically coupled cluster of regions needs an anchor, or its φ (and every
  // charged level with it) floats: shifting φ by s and each η_i by z_i s changes nothing.
  // Regions are coupled across a face by a non-neutral φ law or by any charged species that
  // crosses; a contact anchors its cluster through a gate, a reaction, or a connected ion.
  checkAnchors(regions, materials, interfaces, contacts, species);

  // --- circuit (acts at the right terminal; the left terminal is the reference)
  const circuit = normalizeCircuit(def.circuit, contacts.right, species);

  // --- grid
  if (def.grid !== undefined) need(isObject(def.grid), 'grid must be an object');
  let grid;
  try {
    grid = buildGrid(regions, def.grid ?? {});
  } catch (err) {
    throw new DeviceError(`grid: ${err.message}`);
  }

  return {
    T, RT, F: FARADAY, species, speciesIndex, materials, materialIndex, regions, interfaces, reactions, contacts, circuit, grid, warnings,
  };
}

// Bulk reaction  Σ ν_R R ⇌ Σ ν_P P  with mass-action rate r = k_f Π c_R^ν (1 − e^{−A/RT}), where
// the affinity A = Σ_R ν μ̄ − Σ_P ν μ̄ includes any fixed-activity neutral participants (e.g.
// H₂O, given by its μ). This equals k_f Π c_R − k_b Π c_P with k_b fixed by the standard
// potentials, so equilibrium is exactly A = 0. kf maps material names to rate constants (the
// reaction runs only in those materials).
function normalizeReaction(rdef, path, species, speciesIndex, materials, materialIndex, regions, RT) {
  need(isObject(rdef), `${path} must be an object`);
  const fixed = rdef.fixed ?? {};
  need(isObject(fixed), `${path}.fixed must map fixed-activity participants to their μ (J/mol)`);
  const side = (key) => {
    const m = rdef[key] ?? {};
    need(isObject(m), `${path}.${key} must map participant names to stoichiometric coefficients`);
    const mobile = [];
    let fixedMu = 0, charge = 0;
    for (const [name, nu] of Object.entries(m)) {
      need(Number.isInteger(nu) && nu > 0, `${path}.${key}.${name} must be a positive integer, got ${JSON.stringify(nu)}`);
      if (speciesIndex.has(name)) {
        const i = speciesIndex.get(name);
        mobile.push({ i, nu });
        charge += nu * species[i].z;
      } else {
        need(fixed[name] !== undefined, `${path}.${key}.${name}: not a species, so give its μ in ${path}.fixed (fixed-activity participants are neutral)`);
        fixedMu += nu * finite(fixed[name], `${path}.fixed.${name}`);
      }
    }
    return { mobile, fixedMu, charge };
  };
  const R = side('reactants'), P = side('products');
  need(R.mobile.length + P.mobile.length > 0, `${path}: no mobile participants`);
  need(R.charge === P.charge, `${path}: charge is not balanced (${R.charge} → ${P.charge})`);
  for (const a of R.mobile) {
    need(!P.mobile.some((b) => b.i === a.i), `${path}: '${species[a.i].name}' appears on both sides`);
  }
  need(isObject(rdef.kf), `${path}.kf must map material names to forward rate constants`);
  const kf = new Float64Array(materials.length);
  for (const [mname, v] of Object.entries(rdef.kf)) {
    need(materialIndex.has(mname), `${path}.kf.${mname}: unknown material`);
    const m = materialIndex.get(mname);
    kf[m] = nonNegative(v, `${path}.kf.${mname}`);
    for (const { i } of [...R.mobile, ...P.mobile]) {
      need(materials[m].present[i], `${path}.kf.${mname}: '${species[i].name}' is absent from material '${mname}'`);
    }
  }
  return {
    reactants: R.mobile,
    products: P.mobile,
    fixedA: (R.fixedMu - P.fixedMu) / RT, // fixed participants' share of A/RT
    kf,
  };
}

// Circuit modes. 'voltage' (default): each contact sits at its own V. 'current': a fixed current
// I (A/m², toward +x) leaves through the right terminal, whose voltage floats; I = 0 is open
// circuit. 'load': the right terminal returns to the left one through a resistor R (Ω·m²) in
// series with a source V: I = (V_right − V_left − V) / R.
function normalizeCircuit(cdef, right, species) {
  if (cdef === undefined) return { mode: 'voltage' };
  need(isObject(cdef), 'circuit must be an object with a mode');
  const modes = ['voltage', 'current', 'load'];
  need(modes.includes(cdef.mode), `circuit.mode must be one of ${modes.join(', ')}`);
  if (cdef.mode === 'voltage') return { mode: 'voltage' };
  // The floating terminal voltage is read off a fixed charged terminal species if there is one;
  // otherwise (kinetic or conductance electrode) it becomes an unknown of its own.
  const readout =
    right.terminal !== null && right.species[right.terminal].type === 'equilibrium' && species[right.terminal].z !== 0;
  const exchanges = right.reactions.length > 0 || right.species.some((l) => l.type === 'conductance');
  need(
    readout || exchanges,
    `circuit.mode '${cdef.mode}' needs the right contact to pass current: a charged terminal species in equilibrium, ` +
      'an electrode reaction or a conductance link',
  );
  need(
    readout || right.phi.type === 'capacitive' || right.phi.type === 'neutral',
    `circuit.mode '${cdef.mode}': a floating kinetic electrode needs a capacitive (Stern) or neutral φ law`,
  );
  const base = { terminalUnknown: !readout };
  if (cdef.mode === 'current') return { ...base, mode: 'current', I: finite(cdef.I, 'circuit.I') };
  return {
    ...base,
    mode: 'load',
    R: positive(cdef.R, 'circuit.R'),
    V: cdef.V === undefined ? 0 : finite(cdef.V, 'circuit.V'),
  };
}

function checkAnchors(regions, materials, interfaces, contacts, species) {
  const nR = regions.length;
  const parent = Array.from({ length: nR }, (_, r) => r);
  const find = (r) => (parent[r] === r ? r : (parent[r] = find(parent[r])));
  interfaces.forEach((itf, f) => {
    const coupled = itf.phi.type !== 'neutral' || itf.links.some((l, i) => l.type !== 'blocked' && species[i].z !== 0);
    if (coupled) parent[find(f)] = find(f + 1);
  });
  const anchors = (ct) =>
    ct.phi.type === 'capacitive' ||
    ct.phi.type === 'dipole' ||
    ct.reactions.length > 0 ||
    ct.species.some((l, i) => l.type !== 'blocked' && species[i].z !== 0);
  const anchored = new Set();
  if (anchors(contacts.left)) anchored.add(find(0));
  if (anchors(contacts.right)) anchored.add(find(nR - 1));
  for (let r = 0; r < nR; r++) {
    const mat = materials[regions[r].material];
    const charged = species.some((sp, i) => mat.present[i] && sp.z !== 0);
    if (!charged) continue; // nothing responds to φ there; it is simply not reported
    need(
      anchored.has(find(r)),
      nR === 1 || anchored.size === 0
        ? 'device has no electrostatic anchor: every species is blocked at both ends and no contact has a ' +
            'capacitive (gate) link, so the potential is undetermined. Add a gate or connect a species.'
        : `regions[${r}] (${regions[r].name}) is electrostatically floating: it is cut off by neutral interfaces ` +
            'with no charged species crossing, and nothing anchors its potential.',
    );
  }
}

function normalizeInterface(idef, f, regions, materials, species, speciesIndex, RT) {
  const left = regions[f], right = regions[f + 1];
  const matL = materials[left.material], matR = materials[right.material];
  const where = `interfaces[${f}] (between ${left.name} [${matL.name}] and ${right.name} [${matR.name}])`;
  const same = left.material === right.material;
  if (idef === undefined || idef === null) idef = {};
  need(isObject(idef), `${where} must be an object`);

  // Electrostatic law across the face. 'dipole': φ jumps by the alignment (exact when the grid
  // resolves the double layers). 'neutral': no charge at the face (D = 0) and a free jump, the
  // macroscopic limit; the alignment then drops out. 'capacitive': a Helmholtz layer,
  // D = C (Δφ − dipole). Between two ε = 0 (strictly neutral) materials the default is neutral.
  if (matL.metal || matR.metal) return normalizeMetalInterface(idef, where, matL, matR, species, speciesIndex, RT);
  const bothNeutral = (matL.epsr === 0 && matR.epsr === 0) || matL.phiFree || matR.phiFree;
  const rawPhi = idef.phi ?? (bothNeutral ? 'neutral' : 'dipole');
  const phi = typeof rawPhi === 'string' ? { type: rawPhi } : { ...rawPhi };
  need(['dipole', 'neutral', 'capacitive'].includes(phi.type), `${where}.phi must be 'dipole', 'neutral' or { type: 'capacitive', C }`);
  if (phi.type === 'capacitive') positive(phi.C, `${where}.phi.C`);
  const given = ['dipole', 'step', 'vacuum', 'reaction'].filter((k) => idef[k] !== undefined);
  need(given.length <= 1, `${where}: give exactly one alignment, got ${given.join(' and ')}`);
  if (phi.type === 'neutral') {
    need(given.length === 0, `${where}: a neutral interface has a free φ jump, so an alignment (${given[0]}) would have no effect`);
  } else {
    need(
      given.length === 1 || same,
      `${where}: an interface between different materials needs an alignment ` +
        `({ dipole }, { step: { species, value } } or { vacuum: { left, right } }), or phi: 'neutral' for a macroscopic model. ` +
        'There is no default (no Anderson or Schottky–Mott rule).',
    );
    for (const [mat, reg] of [[matL, left], [matR, right]]) {
      need(
        (mat.epsr > 0 || species.some((sp, i) => mat.present[i] && sp.z !== 0)) && !mat.phiFree,
        `${where}: ${reg.name} has ε = 0 and nothing there responds to φ, so its φ is undefined; use phi: 'neutral'`,
      );
    }
  }

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
  } else if (idef.vacuum !== undefined) {
    // Vacuum-level heuristic: each side's vacuum level sits `offset` volts beyond its anchor
    // (V_vac = V_anchor − offset), and the two vacuum levels are taken to coincide.
    const vac = idef.vacuum;
    need(isObject(vac) && isObject(vac.left) && isObject(vac.right), `${where}.vacuum must be { left: { anchor, offset }, right: { anchor, offset } }`);
    const L = vacuumLevel(vac.left, `${where}.vacuum.left`, matL, species, speciesIndex);
    const R = vacuumLevel(vac.right, `${where}.vacuum.right`, matR, species, speciesIndex);
    dipole = L - R; // V_vac − φ on each side: φ_R − φ_L = (V_vac − φ)_L − (V_vac − φ)_R
  } else if (idef.reaction !== undefined) {
    throw new DeviceError(`${where}.reaction: reaction-based alignment is not supported yet`);
  }
  const sheetCharge = idef.sheetCharge === undefined ? 0 : finite(idef.sheetCharge, `${where}.sheetCharge`);

  // Per-species laws across the face: local equilibrium (default where present on both sides),
  // blocked, or an ohmic interface conductance G (S/m²).
  const links = defaultInterfaceLinks(matL, matR, species);
  if (idef.species !== undefined) {
    need(isObject(idef.species), `${where}.species must map species names to interface links`);
    for (const [sname, raw] of Object.entries(idef.species)) {
      const lpath = `${where}.species.${sname}`;
      need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'`);
      const i = speciesIndex.get(sname);
      const link = typeof raw === 'string' ? { type: raw } : raw;
      need(isObject(link) && INTERFACE_LINK_TYPES.has(link.type), `${lpath}.type must be one of ${[...INTERFACE_LINK_TYPES].join(', ')}`);
      if (link.type !== 'blocked') need(matL.present[i] && matR.present[i], `${lpath}: '${sname}' must be present on both sides`);
      if (link.type === 'conductance') {
        need(species[i].z !== 0, `${lpath}: a conductance link needs a charged species`);
        positive(link.G, `${lpath}.G`);
      }
      links[i] = { ...link };
    }
  }

  // Kinetic transfer across the face (Butler–Volmer form), e.g. ion transfer between two
  // solvents or electron transfer at a heterojunction. Each listed species crosses from the
  // left side to the right side (forward direction) with its coefficient.
  need(idef.reactions === undefined || Array.isArray(idef.reactions), `${where}.reactions must be an array`);
  const transfers = (idef.reactions ?? []).map((rdef, k) => {
    const rpath = `${where}.reactions[${k}]`;
    need(isObject(rdef) && isObject(rdef.transfer), `${rpath} must be { transfer: { species: ν }, k0, alpha }`);
    const list = [];
    for (const [sname, nu] of Object.entries(rdef.transfer)) {
      need(speciesIndex.has(sname), `${rpath}.transfer.${sname}: unknown species`);
      need(Number.isInteger(nu) && nu > 0, `${rpath}.transfer.${sname} must be a positive integer`);
      const i = speciesIndex.get(sname);
      need(matL.present[i] && matR.present[i], `${rpath}.transfer.${sname}: '${sname}' must be present on both sides`);
      need(idef.species?.[sname] === undefined, `${rpath}.transfer.${sname}: '${sname}' also has an interface link`);
      links[i] = { type: 'kinetic' };
      list.push({ i, nu });
    }
    need(list.length > 0, `${rpath}.transfer: no species`);
    return { species: list, k0: positive(rdef.k0, `${rpath}.k0`), alpha: transferCoefficient(rdef.alpha, `${rpath}.alpha`) };
  });
  return { phi, dipole, sheetCharge, links, transfers, electrode: [], metal: null };
}

// A metal: only its conduction carrier, whose single unknown is the Fermi level (μ̄). Its bulk is
// neutral and incompressible, so φ is undefined inside, and transport is ohmic, J = −σ∇V. Any
// charge it holds sits at its surfaces, as a sheet facing a charged interface.
function normalizeMetal(mat, mname, path, species, speciesIndex) {
  const mdef = mat.metal;
  need(isObject(mdef), `${path}.metal must be { species, conductivity }`);
  for (const k of ['epsr', 'species', 'statistics']) {
    need(mat[k] === undefined, `${path}.${k}: a metal takes only { metal: { species, conductivity } }`);
  }
  need(speciesIndex.has(mdef.species), `${path}.metal.species: unknown species ${JSON.stringify(mdef.species)}`);
  const i = speciesIndex.get(mdef.species);
  need(species[i].z !== 0, `${path}.metal.species: the metal's carrier must be charged`);
  const n = species.length;
  const present = new Uint8Array(n);
  present[i] = 1;
  return {
    name: mname,
    epsr: 0,
    present,
    D: new Float64Array(n),
    mu0: new Float64Array(n),
    cRef: new Float64Array(n).fill(1),
    models: [],
    modelOf: new Int32Array(n).fill(-1),
    ideal: true,
    phiFree: true,
    metal: { i, sigma: positive(mdef.conductivity, `${path}.metal.conductivity (S/m)`) },
  };
}

// A face with a metal on one or both sides. The metal has no φ: a capacitive law ties the other
// side's φ to the metal's Fermi level V_F, as at a contact, with zeroCharge = V_F − φ_edge at zero
// charge: D toward the other side = C (V_F − zeroCharge − φ_edge), the metal's surface charge. (No
// pinned 'dipole' law here: an internal metal must hold that charge in a finite capacitance.)
// Electrode reactions take the metal's carriers at its Fermi level.
function normalizeMetalInterface(idef, where, matL, matR, species, speciesIndex, RT) {
  if (idef === undefined || idef === null) idef = {};
  need(isObject(idef), `${where} must be an object`);
  const side = matL.metal ? (matR.metal ? 'both' : 'left') : 'right';
  const other = side === 'left' ? matR : matL;
  const metal = side === 'right' ? matR.metal : matL.metal;
  for (const k of ['dipole', 'step', 'sheetCharge']) {
    need(idef[k] === undefined, `${where}.${k}: a metal has no φ of its own; align with zeroCharge or vacuum (with a 'fermi' anchor on the metal side)`);
  }
  const fieldOutside = side !== 'both' && !other.phiFree && other.epsr > 0;
  const rawPhi = idef.phi ?? (fieldOutside ? undefined : 'neutral');
  need(
    rawPhi !== undefined,
    `${where}.phi: a face between a metal and a material with ε > 0 needs an explicit φ law ` +
      "('neutral', or { type: 'capacitive', C } with zeroCharge or vacuum)",
  );
  const phi = typeof rawPhi === 'string' ? { type: rawPhi } : { ...rawPhi };
  need(
    phi.type !== 'dipole',
    `${where}.phi: a metal region holds its surface charge in a finite capacitance, so use { type: 'capacitive', C } ` +
      '(a large C approaches a pinned barrier), or model the metal as a contact with a dipole law',
  );
  need(['neutral', 'capacitive'].includes(phi.type), `${where}.phi must be 'neutral' or { type: 'capacitive', C }`);
  if (phi.type === 'capacitive') positive(phi.C, `${where}.phi.C`);
  let zeroCharge = 0;
  if (phi.type === 'neutral') {
    need(idef.zeroCharge === undefined && idef.vacuum === undefined, `${where}: a neutral face has no charge, so an alignment would have no effect`);
  } else {
    need(side !== 'both', `${where}.phi: between two metals only 'neutral' applies (neither has a φ)`);
    need(!other.phiFree && (other.epsr > 0 || species.some((sp, i) => other.present[i] && sp.z !== 0)), `${where}.phi: φ is undefined on the non-metal side; use 'neutral'`);
    need((idef.zeroCharge === undefined) !== (idef.vacuum === undefined), `${where}: give exactly one of zeroCharge or vacuum for a ${phi.type} face at a metal`);
    if (idef.zeroCharge !== undefined) zeroCharge = finite(idef.zeroCharge, `${where}.zeroCharge`);
    else {
      const vac = idef.vacuum;
      need(isObject(vac) && isObject(vac.left) && isObject(vac.right), `${where}.vacuum must be { left: { anchor, offset }, right: { anchor, offset } }`);
      const ms = side === 'left' ? vac.left : vac.right, os = side === 'left' ? vac.right : vac.left;
      const mpath = `${where}.vacuum.${side}`, opath = `${where}.vacuum.${side === 'left' ? 'right' : 'left'}`;
      need(ms.anchor === 'fermi', `${mpath}.anchor: a metal's vacuum level is anchored to its Fermi level ('fermi')`);
      const W = finite(ms.offset, `${mpath}.offset (the work function, V)`);
      need(os.anchor !== 'fermi', `${opath}.anchor: only a metal has a Fermi-level anchor`);
      zeroCharge = W + vacuumLevel(os, opath, other, species, speciesIndex); // V_F − φ = W + (V_vac − φ)_other
    }
  }
  // Species laws: the metal's carrier may continue across (e.g. into a semiconductor).
  const links = defaultInterfaceLinks(matL, matR, species);
  if (idef.species !== undefined) {
    need(isObject(idef.species), `${where}.species must map species names to interface links`);
    for (const [sname, raw] of Object.entries(idef.species)) {
      const lpath = `${where}.species.${sname}`;
      need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'`);
      const i = speciesIndex.get(sname);
      const link = typeof raw === 'string' ? { type: raw } : raw;
      need(isObject(link) && INTERFACE_LINK_TYPES.has(link.type), `${lpath}.type must be one of ${[...INTERFACE_LINK_TYPES].join(', ')}`);
      if (link.type !== 'blocked') need(matL.present[i] && matR.present[i], `${lpath}: '${sname}' must be present on both sides`);
      if (link.type === 'conductance') {
        need(species[i].z !== 0, `${lpath}: a conductance link needs a charged species`);
        positive(link.G, `${lpath}.G`);
      }
      links[i] = { ...link };
    }
  }
  need(idef.reactions === undefined || Array.isArray(idef.reactions), `${where}.reactions must be an array`);
  const electrode = (idef.reactions ?? []).map((rdef, k) => {
    const rpath = `${where}.reactions[${k}]`;
    need(side !== 'both', `${rpath}: electrode reactions need a metal on one side only`);
    need(isObject(rdef) && rdef.transfer === undefined, `${rpath}: at a metal face, write an electrode reaction (reactants, electrons, products)`);
    need(species[metal.i].z === -1, `${rpath}: electrode reactions take electrons, so the metal's carrier must have z = −1`);
    const rx = normalizeElectrodeReaction(rdef, rpath, other, species, speciesIndex, RT);
    for (const { i } of [...rx.reactants, ...rx.products]) {
      need(idef.species?.[species[i].name] === undefined, `${rpath}: '${species[i].name}' also has an interface link`);
      links[i] = { type: 'kinetic' };
    }
    links[metal.i] = { type: 'kinetic' };
    return rx;
  });
  return { phi, dipole: 0, zeroCharge, sheetCharge: 0, links, transfers: [], electrode, metal: { side, i: metal.i } };
}

// A material's vacuum level relative to its own φ, V_vac − φ, from an anchor and an offset:
//   anchor = a charged species: its standard level V°_i = φ + μ°_i/(z_i F) (the conduction band
//            for e⁻, the valence band for h⁺, a reversible electrode's level for an ion);
//   anchor = 'phi': the inner potential itself (the offset is then a surface potential).
// offset (V) is the vacuum level's height above the anchor in electron energy, e.g. an electron
// affinity, a work function, an ionisation energy or Trasatti's 4.44 V, so V_vac = V_anchor − offset.
function vacuumLevel(side, path, mat, species, speciesIndex) {
  const offset = finite(side.offset, `${path}.offset (V)`);
  if (side.anchor === 'phi') return -offset;
  need(speciesIndex.has(side.anchor), `${path}.anchor must be 'phi' or a charged species, got ${JSON.stringify(side.anchor)}`);
  const i = speciesIndex.get(side.anchor);
  need(species[i].z !== 0, `${path}.anchor: '${side.anchor}' is neutral, so it has no level in volts`);
  need(mat.present[i], `${path}.anchor: '${side.anchor}' is absent from '${mat.name}'`);
  return mat.mu0[i] / (species[i].z * FARADAY) - offset;
}

function defaultInterfaceLinks(matL, matR, species) {
  return species.map((_, i) => ({ type: matL.present[i] && matR.present[i] ? 'equilibrium' : 'blocked' }));
}

function transferCoefficient(v, path) {
  if (v === undefined) return 0.5;
  need(isFiniteNumber(v) && v > 0 && v < 1, `${path} must be between 0 and 1, got ${JSON.stringify(v)}`);
  return v;
}

// Electrode reaction at a contact, written as reduction when electrons > 0:
//   Σ ν_R R + n e⁻(metal) ⇌ Σ ν_P P
// Species participants live at the contact node; the metal's electrons sit at μ̄_e = −F·V
// (V the contact's terminal voltage); anything else is a fixed-activity participant given by
// its μ. Rate per area, with a = A/RT and standard rate constant k0 (mol/m²/s):
//   r = k0 · Π_R (c/c_ref)^{ν(1−α)} · Π_P (c/c_ref)^{να} · (e^{αa} − e^{−(1−α)a})
// which is mass action with potential-dependent rate constants: exact at A = 0.
function normalizeElectrodeReaction(rdef, path, mat, species, speciesIndex, RT) {
  need(isObject(rdef), `${path} must be an object`);
  const fixed = rdef.fixed ?? {};
  need(isObject(fixed), `${path}.fixed must map fixed-activity participants to their μ (J/mol)`);
  const n = rdef.electrons ?? 0;
  need(Number.isInteger(n), `${path}.electrons must be an integer (electrons taken from the metal)`);
  const side = (key) => {
    const m = rdef[key] ?? {};
    need(isObject(m), `${path}.${key} must map participant names to stoichiometric coefficients`);
    const list = [];
    let fixedMu = 0, charge = 0;
    for (const [name, nu] of Object.entries(m)) {
      need(Number.isInteger(nu) && nu > 0, `${path}.${key}.${name} must be a positive integer`);
      if (speciesIndex.has(name)) {
        const i = speciesIndex.get(name);
        need(mat.present[i], `${path}.${key}.${name}: '${name}' is absent from the end material '${mat.name}'`);
        list.push({ i, nu });
        charge += nu * species[i].z;
      } else {
        need(fixed[name] !== undefined, `${path}.${key}.${name}: not a species, so give its μ in ${path}.fixed`);
        fixedMu += nu * finite(fixed[name], `${path}.fixed.${name}`);
      }
    }
    return { list, fixedMu, charge };
  };
  const R = side('reactants'), P = side('products');
  need(R.charge - n === P.charge, `${path}: charge is not balanced (${R.charge} − ${n} e⁻ → ${P.charge})`);
  need(R.list.length + P.list.length > 0, `${path}: no species participate`);
  return {
    reactants: R.list,
    products: P.list,
    electrons: n,
    fixedA: (R.fixedMu - P.fixedMu) / RT,
    k0: positive(rdef.k0, `${path}.k0`),
    alpha: transferCoefficient(rdef.alpha, `${path}.alpha`),
  };
}

// Reduced potentials ζ of a bath composition in the end material, through its statistics.
function bathZeta(mat, cb, background) {
  const zeta = cb.map((c, i) => (c > 0 ? Math.log(c / mat.cRef[i]) : -800));
  for (const md of mat.models) {
    const z = Float64Array.from(md.idx, (i) => zeta[i]);
    const fixed = md.idx.map((i) => cb[i] > 0);
    const target = Float64Array.from(md.idx, (i) => cb[i]);
    try {
      md.invert(z, fixed, target, background);
    } catch (err) {
      throw new DeviceError(`bath: ${err.message}`);
    }
    md.idx.forEach((i, a) => (zeta[i] = z[a]));
  }
  return zeta;
}

// A contact has a terminal voltage V (set by the circuit; 0 by default) and a link for every
// species and for φ. A fixed charged species sits at V_i = V + offset_i, i.e.
// μ̄_i = z_i F (V + offset_i). The offset belongs to the outside phase: 0 for the terminal
// species (e.g. e⁻ at a metal); μ_M/(nF) for Mⁿ⁺ at a reversible M electrode (0 on table
// conventions, E° being carried by the ion's μ°). It is never defaulted for the others. A fixed neutral species takes an absolute μ̄. A bath computes the offsets from
// a composition, anchored through its reference species.
function normalizeContact(cdef, side, region, materials, species, speciesIndex, RT) {
  const path = `contacts.${side}`;
  const mat = materials[region.material];
  const links = species.map(() => ({ type: 'blocked' }));
  let phi = { type: 'neutral' };
  let terminal = null;
  if (cdef === undefined || cdef === null) return { V: 0, species: links, phi, terminal, reactions: [] };
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
      if (link.type === 'conductance') {
        // Ohmic exchange with an outside reservoir at V_i = V + offset: J = G·(V_out − V_i).
        need(species[i].z !== 0, `${lpath}: a conductance link needs a charged species`);
        const offset = link.offset ?? (terminal === i ? 0 : undefined);
        need(offset !== undefined, `${lpath}.offset: give V_i − V_terminal (V) of the outside reservoir`);
        links[i] = { type: 'conductance', G: positive(link.G, `${lpath}.G`), offset: finite(offset, `${lpath}.offset`) };
        continue;
      }
      if (link.type === 'equilibrium') {
        // In equilibrium with the outside phase, whose levels are known: a Dirichlet condition.
        if (species[i].z === 0) {
          links[i] = { type: 'equilibrium', mu: finite(link.mu, `${lpath}.mu (a neutral species is held at its μ̄, J/mol)`) };
          continue;
        }
        need(link.mu === undefined, `${lpath}: a charged species is held by an offset from the terminal voltage, not by mu`);
        const offset = link.offset ?? (terminal === i ? 0 : undefined);
        need(
          offset !== undefined,
          `${lpath}.offset: give V_i − V_terminal (V) for this species; only the terminal species defaults to 0`,
        );
        links[i] = { type: 'equilibrium', offset: finite(offset, `${lpath}.offset`) };
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
    const zb = bathZeta(mat, cb, region.background);
    const level = (i) => mat.mu0[i] + RT * zb[i]; // μ̄ − zFφ_bath
    const refOffset = bath.offset === undefined ? 0 : finite(bath.offset, `${path}.bath.offset`);
    const beta = refOffset - level(r) / (species[r].z * FARADAY); // φ_bath − V
    for (let i = 0; i < species.length; i++) {
      if (!(cb[i] > 0)) continue;
      links[i] =
        species[i].z === 0
          ? { type: 'equilibrium', mu: level(i) }
          : { type: 'equilibrium', offset: beta + level(i) / (species[i].z * FARADAY) };
    }
    phi = { type: 'bulk' };
  }

  if (cdef.phi !== undefined) {
    const raw = typeof cdef.phi === 'string' ? { type: cdef.phi } : cdef.phi;
    need(isObject(raw), `${path}.phi must be a link type string or an object with a type`);
    need(PHI_LINK_TYPES.has(raw.type), `${path}.phi.type must be one of ${[...PHI_LINK_TYPES].join(', ')}`);
    if (raw.type === 'capacitive' || raw.type === 'dipole') {
      need(!mat.phiFree, `${path}.phi: φ is undefined in '${mat.name}' (only neutral combinations are charged there), so use 'bulk' or 'neutral'`);
      if (raw.type === 'capacitive') positive(raw.C, `${path}.phi.C`);
      need(raw.V === undefined, `${path}.phi.V: the gate voltage is the contact's terminal voltage, ${path}.V`);
      phi = { ...raw };
      if (raw.vacuum !== undefined) {
        // Vacuum-level heuristic: the outside conductor's vacuum level is `outside` volts beyond
        // its terminal level (a work function), so V − W = V_vac = φ_edge + (V_vac − φ)_inside.
        need(raw.zeroCharge === undefined, `${path}.phi: give either zeroCharge or vacuum, not both`);
        const vac = raw.vacuum;
        need(isObject(vac) && isObject(vac.inside), `${path}.phi.vacuum must be { outside, inside: { anchor, offset } }`);
        const W = finite(vac.outside, `${path}.phi.vacuum.outside (the terminal conductor's work function, V)`);
        phi.zeroCharge = W + vacuumLevel(vac.inside, `${path}.phi.vacuum.inside`, mat, species, speciesIndex);
      }
      finite(phi.zeroCharge, `${path}.phi.zeroCharge (the zero-charge alignment: flat-band voltage, pzc or barrier)`);
    } else {
      phi = { ...raw };
    }
  } else if (cdef.bath === undefined && (links.some((l) => l.type !== 'blocked') || cdef.reactions?.length)) {
    throw new DeviceError(
      `${path}.phi: a contact with connected species needs an explicit φ law ('bulk', 'neutral', capacitive or dipole)`,
    );
  }
  need(cdef.reactions === undefined || Array.isArray(cdef.reactions), `${path}.reactions must be an array`);
  const reactions = (cdef.reactions ?? []).map((r, k) =>
    normalizeElectrodeReaction(r, `${path}.reactions[${k}]`, mat, species, speciesIndex, RT),
  );
  if (phi.type === 'bulk') {
    const reacting = (i) => reactions.some((rx) => [...rx.reactants, ...rx.products].some((p) => p.i === i));
    need(
      links.some((l, i) => (l.type !== 'blocked' || reacting(i)) && species[i].z !== 0),
      `${path}.phi: 'bulk' needs at least one connected charged species at this contact`,
    );
  }

  if (terminal !== null) {
    need(links[terminal].type !== 'blocked', `${path}.terminal: '${species[terminal].name}' is blocked at this contact`);
  }
  return { V, species: links, phi, terminal, reactions };
}
