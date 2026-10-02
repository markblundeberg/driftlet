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
const SPECIES_LINK_TYPES = new Set(['blocked', 'equilibrium', 'conductance', 'exchange']);
const INTERFACE_LINK_TYPES = new Set(['equilibrium', 'blocked', 'conductance']);
const PHI_LINK_TYPES = new Set(['bulk', 'neutral', 'capacitive', 'pinned']);
const GRID_FIELDS = ['hmin', 'hmax', 'ratio', 'minCells'];
// The fields of each kind of species link (to an outside phase, or across a face).
const LINK_FIELDS = {
  blocked: ['type'],
  equilibrium: ['type', 'offset', 'mu'],
  conductance: ['type', 'G', 'offset'],
  exchange: ['type', 'k', 'mu'],
};

const isObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);
const isFiniteNumber = (v) => typeof v === 'number' && Number.isFinite(v);

// For an unknown-name error: the names there are.
const known = (index) => ` (the species are ${[...index.keys()].map((k) => `'${k}'`).join(', ')})`;

function need(cond, message) {
  if (!cond) throw new DeviceError(message);
}

// Every key of a spec object must be one it takes: a misspelt or misplaced field would otherwise
// be ignored silently.
function fields(obj, path, allowed) {
  for (const [k, v] of Object.entries(obj)) {
    need(v === undefined || allowed.includes(k), `${path}.${k}: not a field here (${allowed.join(', ')})`);
  }
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
  fields(def, 'device', ['T', 'species', 'materials', 'regions', 'interfaces', 'bulkReactions', 'contacts', 'ports', 'grid']);

  const T = def.T === undefined ? 298.15 : positive(def.T, 'T');
  const RT = GAS_CONSTANT * T;

  // --- species
  need(Array.isArray(def.species) && def.species.length > 0, 'species must be a non-empty array');
  const species = [];
  const speciesIndex = new Map();
  def.species.forEach((sp, i) => {
    const path = `species[${i}]`;
    need(isObject(sp), `${path} must be an object`);
    fields(sp, path, ['name', 'z', 'cRef']);
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
    if (mat.conductor !== undefined) {
      materialIndex.set(mname, materials.length);
      materials.push(normalizeConductor(mat, mname, path, species, speciesIndex));
      continue;
    }
    // ε = 0 makes the material strictly neutral: Poisson becomes local neutrality there.
    fields(mat, path, ['epsr', 'species', 'statistics']);
    const epsr = nonNegative(mat.epsr, `${path}.epsr`);
    need(isObject(mat.species), `${path}.species must be an object mapping species names to parameters`);
    const present = new Uint8Array(nSpecies);
    const D = new Float64Array(nSpecies);
    const mu0 = new Float64Array(nSpecies);
    const cRef = new Float64Array(nSpecies);
    for (const [sname, p] of Object.entries(mat.species)) {
      const spath = `${path}.species.${sname}`;
      need(speciesIndex.has(sname), `${spath}: unknown species '${sname}'${known(speciesIndex)}`);
      need(isObject(p), `${spath} must be an object`);
      fields(p, spath, ['D', 'mu0', 'cRef']);
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
    fields(reg, path, ['name', 'material', 'length', 'fixedCharge', 'c0', 'velocity', 'mixing', 'grid']);
    need(materialIndex.has(reg.material), `${path}.material: unknown material ${JSON.stringify(reg.material)}`);
    const length = positive(reg.length, `${path}.length`);
    const fixedCharge = reg.fixedCharge === undefined ? 0 : finite(reg.fixedCharge, `${path}.fixedCharge`);
    // Imposed flow (m/s, toward +x) carrying every mobile species, and eddy mixing (m²/s).
    const velocity = reg.velocity === undefined ? 0 : finite(reg.velocity, `${path}.velocity`);
    const mixing = reg.mixing === undefined ? 0 : nonNegative(reg.mixing, `${path}.mixing`);
    if (reg.grid !== undefined) {
      need(isObject(reg.grid), `${path}.grid must be an object`);
      fields(reg.grid, `${path}.grid`, GRID_FIELDS);
    }
    const mat = materials[materialIndex.get(reg.material)];
    if (mat.conductor) {
      need(fixedCharge === 0, `${path}.fixedCharge: a conductor region is neutral in bulk`);
      need(reg.c0 === undefined, `${path}.c0: a conductor region has no composition to give, only its carrier's level`);
      need(velocity === 0 && mixing === 0, `${path}: a conductor region has no flow or mixing, only conduction`);
      // Nothing is stored inside a metal (its charge is on its faces), so its Fermi level is
      // exactly linear across it, and one cell is exact.
      need(reg.grid === undefined, `${path}.grid: a conductor region is a single cell (exact for Ohm's law), so it takes no grid options`);
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
        need(speciesIndex.has(sname), `${path}.c0.${sname}: unknown species '${sname}'${known(speciesIndex)}`);
        const i = speciesIndex.get(sname);
        need(mat.present[i], `${path}.c0.${sname}: '${sname}' is absent from material '${mat.name}'`);
        need(
          isFiniteNumber(v) && v > 0,
          `${path}.c0.${sname} must be a concentration > 0 (mol/m³), got ${JSON.stringify(v)}; to have none of '${sname}' here, ` +
            `leave it out of material '${mat.name}' (a region of another material without it, if it's present elsewhere)`,
        );
        c0[i] = v;
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
      cells: mat.conductor ? [length] : undefined,
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

  // --- internal ports (outside phases attached over windows of interior nodes)
  need(def.ports === undefined || Array.isArray(def.ports), 'ports must be an array');
  const ports = (def.ports ?? []).map((pdef, k) => normalizePort(pdef, `ports[${k}]`, regions, materials, species, speciesIndex));

  // Every electrostatically coupled cluster of regions needs an anchor, or its φ (and every
  // charged level with it) floats: shifting φ by s and each η_i by z_i s changes nothing.
  // Regions are coupled across a face by a non-neutral φ law or by any charged species that
  // crosses; a contact anchors its cluster through a gate, a reaction, or a connected ion.
  checkAnchors(regions, materials, interfaces, contacts, species, ports);

  // --- terminals: the two contacts and every port, each held at a voltage or driven by a current
  const terminals = [
    { name: 'left', kind: 'contact', side: 'left', drive: contacts.left.drive },
    { name: 'right', kind: 'contact', side: 'right', drive: contacts.right.drive },
    ...ports.map((port, k) => ({ name: port.name, kind: 'port', index: k, drive: port.drive })),
  ];
  const names = new Set();
  for (const t of terminals) {
    need(!names.has(t.name), `ports: the name '${t.name}' is taken (by a contact or another port)`);
    names.add(t.name);
  }
  // A contact that passes nothing can't be driven by a current.
  for (const side of ['left', 'right']) {
    const ct = contacts[side];
    const passes = ct.species.some((l) => l.type !== 'blocked') || ct.phi.type === 'capacitive' || ct.phi.type === 'pinned';
    need(passes || ct.drive.kind === 'V', `contacts.${side}.I: this contact passes no current (no linked species, no gate)`);
    ct.passes = passes;
  }
  need(
    terminals.some((t) => t.drive.kind === 'V'),
    'every terminal is driven by a current, so the device\'s overall level floats: hold at least one at a voltage V',
  );

  // --- grid
  if (def.grid !== undefined) {
    need(isObject(def.grid), 'grid must be an object');
    fields(def.grid, 'grid', GRID_FIELDS);
  }
  let grid;
  try {
    grid = buildGrid(regions, def.grid ?? {});
  } catch (err) {
    throw new DeviceError(`grid: ${err.message}`);
  }

  // Port windows: the nodes of the port's region within [from, to] of its left end.
  ports.forEach((port, k) => {
    const nodes = [];
    const x0 = grid.x[grid.regionStart[port.region]];
    for (let g = grid.regionStart[port.region]; g <= grid.regionEnd[port.region]; g++) {
      const d = grid.x[g] - x0;
      if (d >= port.from - 1e-15 * port.span && d <= port.to + 1e-15 * port.span) nodes.push(g);
    }
    need(nodes.length > 0, `ports[${k}]: the window [${port.from}, ${port.to}] m holds no grid node; widen it or refine the grid`);
    port.nodes = Int32Array.from(nodes);
  });

  return {
    T, RT, F: FARADAY, species, speciesIndex, materials, materialIndex, regions, interfaces, reactions, contacts, terminals, grid, warnings, ports,
  };
}

// Bulk reaction  Σ ν_R R ⇌ Σ ν_P P  with mass-action rate r = k_f Π c_R^ν (1 − e^{−A/RT}), where
// the affinity A = Σ_R ν μ̄ − Σ_P ν μ̄ includes any fixed-activity neutral participants (e.g.
// H₂O, given by its μ). This equals k_f Π c_R − k_b Π c_P with k_b fixed by the standard
// potentials, so equilibrium is exactly A = 0. kf maps material names to rate constants (the
// reaction runs only in those materials).
function normalizeReaction(rdef, path, species, speciesIndex, materials, materialIndex, regions, RT) {
  need(isObject(rdef), `${path} must be an object`);
  fields(rdef, path, ['nu', 'fixed', 'kf']);
  const st = stoichiometry(rdef.nu, `${path}.nu`, rdef.fixed, `${path}.fixed`, species, speciesIndex, RT);
  need(st.list.length > 0, `${path}: no mobile participants`);
  need(st.charge === 0, `${path}: charge is not balanced (Σ ν z = ${st.charge})`);
  need(isObject(rdef.kf), `${path}.kf must map material names to forward rate constants`);
  const kf = new Float64Array(materials.length);
  for (const [mname, v] of Object.entries(rdef.kf)) {
    need(materialIndex.has(mname), `${path}.kf.${mname}: unknown material`);
    const m = materialIndex.get(mname);
    kf[m] = nonNegative(v, `${path}.kf.${mname}`);
    for (const { i } of st.list) {
      need(materials[m].present[i], `${path}.kf.${mname}: '${species[i].name}' is absent from material '${mname}'`);
    }
  }
  const reactants = st.list.filter((p) => p.nu < 0).map(({ i, nu }) => ({ i, nu: -nu }));
  const products = st.list.filter((p) => p.nu > 0);
  return {
    reactants,
    products,
    fixedA: st.fixedA, // fixed participants' share of A/RT
    kf,
    // Species made only from (or turned only into) fixed reservoirs, such as photogeneration
    // from a photon reservoir: a source that can hold the device far from equilibrium.
    generation: st.nFixed > 0 && (reactants.length === 0 || products.length === 0),
  };
}

// A signed stoichiometry { name: ν } (ν < 0 consumed, ν > 0 produced by the forward reaction).
// Species become { i, nu }; any other name is a fixed-activity neutral whose μ (J/mol) is given
// in `fixed`, and contributes −ν μ/RT to the affinity a = A/RT = −Σ ν μ̄/RT.
function stoichiometry(map, path, fixed, fixedPath, species, speciesIndex, RT) {
  map ??= {};
  fixed ??= {};
  need(isObject(map), `${path} must map participants to signed stoichiometric coefficients (ν < 0 consumed)`);
  need(isObject(fixed), `${fixedPath} must map fixed-activity participants to their μ (J/mol)`);
  const list = [];
  let fixedA = 0, charge = 0, nFixed = 0;
  for (const [name, nu] of Object.entries(map)) {
    need(Number.isInteger(nu) && nu !== 0, `${path}.${name} must be a non-zero integer, got ${JSON.stringify(nu)}`);
    if (speciesIndex.has(name)) {
      const i = speciesIndex.get(name);
      list.push({ i, nu });
      charge += nu * species[i].z;
    } else {
      need(fixed[name] !== undefined, `${path}.${name}: not a species, so give its μ in ${fixedPath} (fixed-activity participants are neutral)`);
      fixedA -= (nu * finite(fixed[name], `${fixedPath}.${name}`)) / RT;
      nFixed++;
    }
  }
  return { list, fixedA, charge, nFixed };
}

// Reactions at a face, Butler–Volmer, with participants on either side (signed stoichiometry
// per side, ν < 0 consumed by the forward reaction). Each side's participants are species
// present there (on a conductor's side, only its carrier) or fixed-activity neutrals. Rate per
// area, with a = A/RT = −Σ ν μ̄/RT and k0 in mol/(m²·s):
//   r = k0 Π_{ν<0} (c/c_ref)^{|ν|(1−α)} Π_{ν>0} (c/c_ref)^{να} (e^{αa} − e^{−(1−α)a})
// (a conductor's carrier has activity 1, so no factor). That's mass action with rate constants
// that depend on the electrical part of the affinity, and exactly zero at A = 0.
function normalizeFaceReactions(idef, where, matL, matR, species, speciesIndex, RT) {
  need(idef.reactions === undefined || Array.isArray(idef.reactions), `${where}.reactions must be an array`);
  return (idef.reactions ?? []).map((rdef, k) => {
    const rpath = `${where}.reactions[${k}]`;
    need(isObject(rdef), `${rpath} must be { left, right, fixed, k0, alpha }`);
    fields(rdef, rpath, ['left', 'right', 'fixed', 'k0', 'alpha']);
    const part = [];
    let fixedA = 0, charge = 0;
    for (const [key, mat] of [['left', matL], ['right', matR]]) {
      const st = stoichiometry(rdef[key], `${rpath}.${key}`, rdef.fixed, `${rpath}.fixed`, species, speciesIndex, RT);
      for (const { i, nu } of st.list) {
        need(mat.present[i], `${rpath}.${key}.${species[i].name}: absent from '${mat.name}'`);
        part.push({ i, nu, side: key === 'left' ? 0 : 1 });
      }
      fixedA += st.fixedA;
      charge += st.charge;
    }
    need(part.length > 0, `${rpath}: no species participate`);
    need(charge === 0, `${rpath}: charge is not balanced (Σ ν z = ${charge})`);
    return { part, fixedA, k0: positive(rdef.k0, `${rpath}.k0`), alpha: transferCoefficient(rdef.alpha, `${rpath}.alpha`) };
  });
}

// A species that takes part in a reaction at a face and exists on both sides needs its link
// there given explicitly: the default ('equilibrium') would cross it freely alongside the
// reaction and short-circuit the kinetics. ('blocked' if it crosses only through the reaction.)
function checkReactingLinks(reactions, idef, where, matL, matR, species) {
  for (const rx of reactions) {
    for (const { i } of rx.part) {
      const name = species[i].name;
      need(
        idef.species?.[name] !== undefined || !(matL.present[i] && matR.present[i]),
        `${where}.species.${name}: '${name}' takes part in a reaction here and exists on both sides, so give its link ` +
          "explicitly ('blocked' if it crosses only through the reaction)",
      );
    }
  }
}

// A terminal's drive: held at a voltage V (behind a series resistance R, Ω·m², if given), or
// driven by a current I (A/m², into the device). V and I are numbers or piecewise-linear
// waveforms { t: [...], values: [...], repeat }. Neither given: held at V = 0.
function normalizeDrive(d, path) {
  need(d.V === undefined || d.I === undefined, `${path}: give V or I, not both`);
  if (d.I !== undefined) {
    need(d.R === undefined, `${path}.R: a series resistance goes with a voltage source V`);
    return { kind: 'I', src: normalizeSource(d.I, `${path}.I`), R: 0 };
  }
  return {
    kind: 'V',
    src: normalizeSource(d.V ?? 0, `${path}.V`),
    R: d.R === undefined ? 0 : positive(d.R, `${path}.R (Ω·m²)`),
  };
}

// A source value: a constant, or a piecewise-linear waveform through the points (t, value),
// constant beyond them, or periodic with period t_last − t_0 when repeat is true.
function normalizeSource(v, path) {
  if (typeof v === 'number') return { value: finite(v, path) };
  need(isObject(v), `${path} must be a number or a waveform { t: [...], values: [...], repeat }`);
  fields(v, path, ['t', 'values', 'repeat']);
  need(Array.isArray(v.t) && Array.isArray(v.values) && v.t.length >= 1 && v.t.length === v.values.length, `${path}: t and values need the same length, at least 1`);
  v.t.forEach((t, k) => {
    finite(t, `${path}.t[${k}]`);
    need(k === 0 || t > v.t[k - 1], `${path}.t must increase strictly`);
  });
  v.values.forEach((x, k) => finite(x, `${path}.values[${k}]`));
  need(v.repeat === undefined || typeof v.repeat === 'boolean', `${path}.repeat must be true or false`);
  need(!v.repeat || v.t.length >= 2, `${path}: a repeating waveform needs at least two points`);
  return { t: Float64Array.from(v.t), values: Float64Array.from(v.values), repeat: v.repeat === true };
}

/**
 * The terminals' drives alone (in terminal order: left, right, then the ports), from a definition
 * whose structure is unchanged: for a fast update of sources.
 */
export function normalizeDrives(def, model) {
  const drives = [
    normalizeDrive(def.contacts?.left ?? {}, 'contacts.left'),
    normalizeDrive(def.contacts?.right ?? {}, 'contacts.right'),
    ...(def.ports ?? []).map((p, k) => normalizeDrive(p, `ports[${k}]`)),
  ];
  for (const [k, side] of [[0, 'left'], [1, 'right']]) {
    need(model.contacts[side].passes || drives[k].kind === 'V', `contacts.${side}.I: this contact passes no current (no linked species, no gate)`);
  }
  need(drives.some((d) => d.kind === 'V'), "every terminal is driven by a current, so the device's overall level floats: hold at least one at a voltage V");
  return drives;
}

/** A source's value at time t (s). */
export function sourceAt(src, t) {
  if (src.value !== undefined) return src.value;
  const { t: ts, values: vs, repeat } = src, n = ts.length;
  if (repeat) {
    const T = ts[n - 1] - ts[0];
    t = ts[0] + ((((t - ts[0]) % T) + T) % T);
  }
  if (t <= ts[0]) return vs[0];
  if (t >= ts[n - 1]) return vs[n - 1];
  let k = 1;
  while (ts[k] < t) k++;
  const w = (t - ts[k - 1]) / (ts[k] - ts[k - 1]);
  return vs[k - 1] + w * (vs[k] - vs[k - 1]);
}

/** The first breakpoint of a source strictly after time t (Infinity if none). */
export function nextBreakpoint(src, t) {
  if (src.value !== undefined) return Infinity;
  const { t: ts, repeat } = src, n = ts.length;
  if (!repeat) {
    for (let k = 0; k < n; k++) if (ts[k] > t * (1 + 1e-12) + 1e-300) return ts[k];
    return Infinity;
  }
  const T = ts[n - 1] - ts[0], cycle = Math.floor((t - ts[0]) / T);
  for (let c = cycle; c <= cycle + 1; c++) {
    for (let k = 0; k < n; k++) {
      const tb = ts[k] + c * T;
      if (tb > t + 1e-12 * Math.max(Math.abs(t), T)) return tb;
    }
  }
  return Infinity;
}

// An internal port: an outside phase with known levels (V_i = V + offset_i, or μ for neutral
// species, as at a contact), exchanging with every node in a window of one region.
// Links per species: 'equilibrium' (μ̄ held at the outside level throughout the window),
// { type: 'conductance', G } for charged species (G in S/m³: a source G (V_out − V_i)/(zF) per
// volume), { type: 'exchange', k, mu } for neutral ones (k in mol/(m³·s): a source
// k (μ_out − μ)/RT per volume), or 'blocked' (the default).
function normalizePort(pdef, path, regions, materials, species, speciesIndex) {
  need(isObject(pdef), `${path} must be an object`);
  fields(pdef, path, ['name', 'region', 'from', 'to', 'V', 'I', 'R', 'terminal', 'species']);
  let r;
  if (Number.isInteger(pdef.region)) r = pdef.region;
  else r = regions.findIndex((reg) => reg.name === pdef.region);
  need(r >= 0 && r < regions.length, `${path}.region: give a region's name or index, got ${JSON.stringify(pdef.region)}`);
  const reg = regions[r], mat = materials[reg.material];
  // A metal is one cell with no interior, so a port attaches to the whole of it, as a wire would.
  if (mat.conductor) need(pdef.from === undefined && pdef.to === undefined, `${path}: a port on a conductor attaches to all of it, so give no window`);
  const from = pdef.from === undefined ? 0 : nonNegative(pdef.from, `${path}.from`);
  const to = pdef.to === undefined ? reg.length : finite(pdef.to, `${path}.to`);
  need(to >= from && to <= reg.length * (1 + 1e-12), `${path}: the window [from, to] must lie within the region (0 to ${reg.length} m)`);
  const drive = normalizeDrive(pdef, path);
  let terminal = null;
  if (pdef.terminal !== undefined) {
    need(speciesIndex.has(pdef.terminal), `${path}.terminal: unknown species ${JSON.stringify(pdef.terminal)}${known(speciesIndex)}`);
    terminal = speciesIndex.get(pdef.terminal);
  }
  need(isObject(pdef.species), `${path}.species must map species names to port links`);
  const links = species.map(() => ({ type: 'blocked' }));
  for (const [sname, raw] of Object.entries(pdef.species)) {
    const lpath = `${path}.species.${sname}`;
    need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'${known(speciesIndex)}`);
    const i = speciesIndex.get(sname);
    const link = typeof raw === 'string' ? { type: raw } : raw;
    need(isObject(link) && ['blocked', 'equilibrium', 'conductance', 'exchange'].includes(link.type), `${lpath}.type must be one of blocked, equilibrium, conductance, exchange`);
    fields(link, lpath, LINK_FIELDS[link.type]);
    if (link.type === 'blocked') continue;
    need(mat.present[i], `${lpath}: '${sname}' is absent from ${reg.name} (material '${mat.name}')`);
    if (mat.conductor) need(i === mat.conductor.i, `${lpath}: a conductor exchanges only its carrier, ${species[mat.conductor.i].name}`);
    const z = species[i].z;
    // The outside level: an offset from V for charged species, an absolute μ for neutral ones.
    let level;
    if (z === 0) {
      need(link.type !== 'conductance', `${lpath}: a neutral species exchanges by { type: 'exchange', k, mu }`);
      level = { mu: finite(link.mu, `${lpath}.mu (the outside μ, J/mol)`) };
    } else {
      need(link.type !== 'exchange', `${lpath}: a charged species exchanges by { type: 'conductance', G }`);
      need(link.mu === undefined, `${lpath}: a charged species is held by an offset from the port voltage, not by mu`);
      const offset = link.offset ?? (terminal === i ? 0 : undefined);
      need(offset !== undefined, `${lpath}.offset: give V_i − V (V); only the terminal species defaults to 0`);
      level = { offset: finite(offset, `${lpath}.offset`) };
    }
    // On a metal, G is a lumped conductance per area (S/m²), spread over its thickness.
    if (link.type === 'conductance') level.G = mat.conductor ? positive(link.G, `${lpath}.G (S/m², for a conductor)`) / reg.length : positive(link.G, `${lpath}.G (S/m³)`);
    if (link.type === 'exchange') level.k = positive(link.k, `${lpath}.k (mol/(m³·s))`);
    links[i] = { type: link.type, ...level };
  }
  need(links.some((l) => l.type !== 'blocked'), `${path}.species: the port exchanges no species`);
  return { name: pdef.name ?? path, region: r, from, to, span: reg.length, drive, terminal, species: links };
}

function checkAnchors(regions, materials, interfaces, contacts, species, ports = []) {
  const nR = regions.length;
  const parent = Array.from({ length: nR }, (_, r) => r);
  const find = (r) => (parent[r] === r ? r : (parent[r] = find(parent[r])));
  interfaces.forEach((itf, f) => {
    // (a reaction moving charge across the face couples the two sides too)
    const charged = (rx, side) => rx.part.some((p) => p.side === side && species[p.i].z !== 0);
    const coupled =
      itf.phi.type !== 'neutral' ||
      itf.links.some((l, i) => l.type !== 'blocked' && species[i].z !== 0) ||
      itf.reactions.some((rx) => charged(rx, 0) && charged(rx, 1));
    if (coupled) parent[find(f)] = find(f + 1);
  });
  const anchors = (ct) =>
    ct.phi.type === 'capacitive' ||
    ct.phi.type === 'pinned' ||
    ct.species.some((l, i) => l.type !== 'blocked' && species[i].z !== 0);
  const anchored = new Set();
  if (anchors(contacts.left)) anchored.add(find(0));
  if (anchors(contacts.right)) anchored.add(find(nR - 1));
  for (const port of ports) if (port.species.some((l, i) => l.type !== 'blocked' && species[i].z !== 0)) anchored.add(find(port.region));
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

  // Electrostatic law across the face. 'pinned': φ jumps by the alignment (exact when the grid
  // resolves the double layers). 'neutral': no charge at the face (D = 0) and a free jump, the
  // macroscopic limit; the alignment then drops out. 'capacitive': a Helmholtz layer,
  // D = C (Δφ − dipole). Between two ε = 0 (strictly neutral) materials the default is neutral.
  if (matL.conductor || matR.conductor) return normalizeConductorInterface(idef, where, matL, matR, species, speciesIndex, RT);
  fields(idef, where, ['phi', 'dipole', 'step', 'sheetCharge', 'species', 'reactions']);
  const bothNeutral = (matL.epsr === 0 && matR.epsr === 0) || matL.phiFree || matR.phiFree;
  const rawPhi = idef.phi ?? (bothNeutral ? 'neutral' : 'pinned');
  const phi = typeof rawPhi === 'string' ? { type: rawPhi } : { ...rawPhi };
  need(isObject(phi), `${where}.phi must be a law name or { type, C }`);
  fields(phi, `${where}.phi`, phi.type === 'capacitive' ? ['type', 'C'] : ['type']);
  need(['pinned', 'neutral', 'capacitive'].includes(phi.type), `${where}.phi must be 'pinned', 'neutral' or { type: 'capacitive', C }`);
  if (phi.type === 'capacitive') positive(phi.C, `${where}.phi.C`);
  const given = ['dipole', 'step'].filter((k) => idef[k] !== undefined);
  need(given.length <= 1, `${where}: give exactly one alignment, got ${given.join(' and ')}`);
  if (phi.type === 'neutral') {
    need(given.length === 0, `${where}: a neutral interface has a free φ jump, so an alignment (${given[0]}) would have no effect`);
  } else {
    need(
      given.length === 1 || same,
      `${where}: an interface between different materials needs an alignment ` +
        "({ dipole } or { step: { species, value } }; driftlet/kit's vacuumDipole estimates one), or phi: 'neutral' for a macroscopic model. " +
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
    fields(st, `${where}.step`, ['species', 'value']);
    need(speciesIndex.has(st.species), `${where}.step.species: unknown species ${JSON.stringify(st.species)}${known(speciesIndex)}`);
    const i = speciesIndex.get(st.species);
    const z = species[i].z;
    need(z !== 0, `${where}.step.species: '${st.species}' is neutral; alignment needs a charged species`);
    need(matL.present[i] && matR.present[i], `${where}.step.species: '${st.species}' must be present on both sides`);
    const value = finite(st.value, `${where}.step.value`);
    // value = (μ°_R + zFφ_R) − (μ°_L + zFφ_L), so φ_R − φ_L = (value − Δμ°) / (zF)
    dipole = (value - (matR.mu0[i] - matL.mu0[i])) / (z * FARADAY);
  }
  const sheetCharge = idef.sheetCharge === undefined ? 0 : finite(idef.sheetCharge, `${where}.sheetCharge`);

  // Per-species laws across the face: local equilibrium (default where present on both sides),
  // blocked, or an ohmic interface conductance G (S/m²).
  const links = defaultInterfaceLinks(matL, matR, species);
  if (idef.species !== undefined) {
    need(isObject(idef.species), `${where}.species must map species names to interface links`);
    for (const [sname, raw] of Object.entries(idef.species)) {
      const lpath = `${where}.species.${sname}`;
      need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'${known(speciesIndex)}`);
      const i = speciesIndex.get(sname);
      const link = typeof raw === 'string' ? { type: raw } : raw;
      need(isObject(link) && INTERFACE_LINK_TYPES.has(link.type), `${lpath}.type must be one of ${[...INTERFACE_LINK_TYPES].join(', ')}`);
      fields(link, lpath, link.type === 'conductance' ? ['type', 'G'] : ['type']);
      if (link.type !== 'blocked') need(matL.present[i] && matR.present[i], `${lpath}: '${sname}' must be present on both sides`);
      if (link.type === 'conductance') {
        need(species[i].z !== 0, `${lpath}: a conductance link needs a charged species`);
        positive(link.G, `${lpath}.G`);
      }
      links[i] = { ...link };
    }
  }

  const reactions = normalizeFaceReactions(idef, where, matL, matR, species, speciesIndex, RT);
  checkReactingLinks(reactions, idef, where, matL, matR, species);
  return { phi, dipole, sheetCharge, links, reactions, conductor: null };
}

// A conductor (a metal, or a fast ion conductor): only its one mobile carrier, whose single
// unknown is its μ̄ (for a metal, the Fermi level). Its bulk is neutral and incompressible, so φ
// is undefined inside, and transport is ohmic, J = −σ∇V. Any charge it holds sits at its
// surfaces, as a sheet facing a charged interface. (The solver still calls it a metal.)
function normalizeConductor(mat, mname, path, species, speciesIndex) {
  const mdef = mat.conductor;
  need(isObject(mdef), `${path}.conductor must be { species, conductivity }`);
  for (const k of ['epsr', 'species', 'statistics']) {
    need(mat[k] === undefined, `${path}.${k}: a conductor takes only { conductor: { species, conductivity } }`);
  }
  need(speciesIndex.has(mdef.species), `${path}.conductor.species: unknown species ${JSON.stringify(mdef.species)}${known(speciesIndex)}`);
  const i = speciesIndex.get(mdef.species);
  need(species[i].z !== 0, `${path}.conductor.species: the conductor's carrier must be charged`);
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
    conductor: { i, sigma: positive(mdef.conductivity, `${path}.conductor.conductivity (S/m)`) },
  };
}

// A face with a metal on one or both sides. The metal has no φ: a capacitive law ties the other
// side's φ to the metal's Fermi level V_F, as at a contact, with zeroCharge = V_F − φ_edge at zero
// charge: D toward the other side = C (V_F − zeroCharge − φ_edge), the metal's surface charge. (No
// 'pinned' law here: an internal conductor must hold that charge in a finite capacitance.)
// Electrode reactions take the metal's carriers at its Fermi level.
function normalizeConductorInterface(idef, where, matL, matR, species, speciesIndex, RT) {
  if (idef === undefined || idef === null) idef = {};
  need(isObject(idef), `${where} must be an object`);
  fields(idef, where, ['phi', 'zeroCharge', 'dipole', 'step', 'sheetCharge', 'species', 'reactions']); // (the alignment fields only to say why not)
  const side = matL.conductor ? (matR.conductor ? 'both' : 'left') : 'right';
  const other = side === 'left' ? matR : matL;
  const metal = side === 'right' ? matR.conductor : matL.conductor;
  for (const k of ['dipole', 'step', 'sheetCharge']) {
    need(idef[k] === undefined, `${where}.${k}: a conductor has no φ of its own; align with zeroCharge`);
  }
  const fieldOutside = side !== 'both' && !other.phiFree && other.epsr > 0;
  const rawPhi = idef.phi ?? (fieldOutside ? undefined : 'neutral');
  need(
    rawPhi !== undefined,
    `${where}.phi: a face between a conductor and a material with ε > 0 needs an explicit φ law ` +
      "('neutral', or { type: 'capacitive', C } with zeroCharge)",
  );
  const phi = typeof rawPhi === 'string' ? { type: rawPhi } : { ...rawPhi };
  need(isObject(phi), `${where}.phi must be a law name or { type, C }`);
  fields(phi, `${where}.phi`, phi.type === 'capacitive' ? ['type', 'C'] : ['type']);
  need(
    phi.type !== 'pinned',
    `${where}.phi: a conductor region holds its surface charge in a finite capacitance, so use { type: 'capacitive', C } ` +
      '(a large C approaches a pinned barrier), or model it as a contact with a pinned law',
  );
  need(['neutral', 'capacitive'].includes(phi.type), `${where}.phi must be 'neutral' or { type: 'capacitive', C }`);
  if (phi.type === 'capacitive') positive(phi.C, `${where}.phi.C`);
  let zeroCharge = 0;
  if (phi.type === 'neutral') {
    need(idef.zeroCharge === undefined, `${where}: a neutral face has no charge, so an alignment would have no effect`);
  } else {
    need(side !== 'both', `${where}.phi: between two conductors only 'neutral' applies (neither has a φ)`);
    need(!other.phiFree && (other.epsr > 0 || species.some((sp, i) => other.present[i] && sp.z !== 0)), `${where}.phi: φ is undefined on the other side; use 'neutral'`);
    zeroCharge = finite(idef.zeroCharge, `${where}.zeroCharge (V_F − φ_edge at zero charge, V)`);
  }
  // Species laws: the metal's carrier may continue across (e.g. into a semiconductor).
  const links = defaultInterfaceLinks(matL, matR, species);
  if (idef.species !== undefined) {
    need(isObject(idef.species), `${where}.species must map species names to interface links`);
    for (const [sname, raw] of Object.entries(idef.species)) {
      const lpath = `${where}.species.${sname}`;
      need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'${known(speciesIndex)}`);
      const i = speciesIndex.get(sname);
      const link = typeof raw === 'string' ? { type: raw } : raw;
      need(isObject(link) && INTERFACE_LINK_TYPES.has(link.type), `${lpath}.type must be one of ${[...INTERFACE_LINK_TYPES].join(', ')}`);
      fields(link, lpath, link.type === 'conductance' ? ['type', 'G'] : ['type']);
      if (link.type !== 'blocked') need(matL.present[i] && matR.present[i], `${lpath}: '${sname}' must be present on both sides`);
      if (link.type === 'conductance') {
        need(species[i].z !== 0, `${lpath}: a conductance link needs a charged species`);
        positive(link.G, `${lpath}.G`);
      }
      links[i] = { ...link };
    }
  }
  const reactions = normalizeFaceReactions(idef, where, matL, matR, species, speciesIndex, RT);
  checkReactingLinks(reactions, idef, where, matL, matR, species);
  return { phi, dipole: 0, zeroCharge, sheetCharge: 0, links, reactions, conductor: { side, i: metal.i } };
}

function defaultInterfaceLinks(matL, matR, species) {
  return species.map((_, i) => ({ type: matL.present[i] && matR.present[i] ? 'equilibrium' : 'blocked' }));
}

function transferCoefficient(v, path) {
  if (v === undefined) return 0.5;
  need(isFiniteNumber(v) && v > 0 && v < 1, `${path} must be between 0 and 1, got ${JSON.stringify(v)}`);
  return v;
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

// A contact has a terminal voltage V (held, or floating under a current drive) and a link for every
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
  if (cdef === undefined || cdef === null) return { drive: normalizeDrive({}, path), species: links, phi, terminal };
  need(isObject(cdef), `${path} must be an object`);
  fields(cdef, path, ['V', 'I', 'R', 'terminal', 'species', 'bath', 'phi', 'zeroCharge']);
  const drive = normalizeDrive(cdef, path);

  if (cdef.terminal !== undefined) {
    need(speciesIndex.has(cdef.terminal), `${path}.terminal: unknown species ${JSON.stringify(cdef.terminal)}${known(speciesIndex)}`);
    terminal = speciesIndex.get(cdef.terminal);
  }

  if (cdef.species !== undefined) {
    need(isObject(cdef.species), `${path}.species must be an object mapping species names to links`);
    for (const [sname, raw] of Object.entries(cdef.species)) {
      const lpath = `${path}.species.${sname}`;
      need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'${known(speciesIndex)}`);
      const i = speciesIndex.get(sname);
      const link = typeof raw === 'string' ? { type: raw } : raw;
      need(isObject(link), `${lpath} must be a link type string or an object with a type`);
      need(SPECIES_LINK_TYPES.has(link.type), `${lpath}.type must be one of ${[...SPECIES_LINK_TYPES].join(', ')}`);
      fields(link, lpath, LINK_FIELDS[link.type]);
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
      if (link.type === 'exchange') {
        // A neutral species exchanging with the outside at its μ: N_in = k (μ_out − μ)/RT.
        need(species[i].z === 0, `${lpath}: a charged species exchanges by { type: 'conductance', G }`);
        links[i] = { type: 'exchange', k: positive(link.k, `${lpath}.k (mol/(m²·s))`), mu: finite(link.mu, `${lpath}.mu (the outside μ, J/mol)`) };
        continue;
      }
      if (link.type === 'conductance') need(species[i].z !== 0, `${lpath}: a neutral species exchanges by { type: 'exchange', k, mu }`);
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
    fields(bath, `${path}.bath`, ['c', 'reference', 'offset']);
    need(speciesIndex.has(bath.reference), `${path}.bath.reference: unknown species ${JSON.stringify(bath.reference)}${known(speciesIndex)}`);
    const r = speciesIndex.get(bath.reference);
    need(species[r].z !== 0, `${path}.bath.reference: the reference species must be charged`);
    need(bath.c[bath.reference] !== undefined, `${path}.bath.reference: '${bath.reference}' must be in the bath`);
    need(cdef.terminal === undefined || terminal === r, `${path}.terminal must be the bath's reference species`);
    terminal = r;
    let charge = region.fixedCharge / FARADAY, scale = Math.abs(charge);
    const cb = new Float64Array(species.length);
    for (const [sname, v] of Object.entries(bath.c)) {
      need(speciesIndex.has(sname), `${path}.bath.c.${sname}: unknown species '${sname}'${known(speciesIndex)}`);
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
    fields(raw, `${path}.phi`, raw.type === 'capacitive' ? ['type', 'C'] : ['type']);
    if (raw.type === 'capacitive' || raw.type === 'pinned') {
      need(!mat.phiFree, `${path}.phi: φ is undefined in '${mat.name}' (only neutral combinations are charged there), so use 'bulk' or 'neutral'`);
      if (raw.type === 'capacitive') positive(raw.C, `${path}.phi.C`);
      // The alignment sits beside the law, as at a conductor's face: V − φ_edge at zero charge.
      phi = { ...raw, zeroCharge: finite(cdef.zeroCharge, `${path}.zeroCharge (V − φ_edge at zero charge: flat-band voltage, pzc or barrier)`) };
    } else {
      need(cdef.zeroCharge === undefined, `${path}.zeroCharge: only a capacitive or pinned φ law takes an alignment`);
      phi = { ...raw };
    }
  } else if (cdef.bath === undefined && links.some((l) => l.type !== 'blocked')) {
    throw new DeviceError(
      `${path}.phi: a contact with connected species needs an explicit φ law ('bulk', 'neutral', capacitive or pinned)`,
    );
  }
  if (phi.type === 'bulk') {
    need(
      links.some((l, i) => l.type !== 'blocked' && species[i].z !== 0),
      `${path}.phi: 'bulk' needs at least one connected charged species at this contact`,
    );
  }

  if (terminal !== null) {
    need(links[terminal].type !== 'blocked', `${path}.terminal: '${species[terminal].name}' is blocked at this contact`);
  }
  return { drive, species: links, phi, terminal };
}
