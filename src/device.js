// Device definition: validation and normalisation.
//
// A definition is a plain, serialisable object (so it can be posted to a Worker). This module
// checks it thoroughly and turns it into an internal model with index-based lookups. Every
// problem is a DeviceError whose message names the offending path. There are no silent
// defaults for physically meaningful choices: interface alignments, reference concentrations
// and standard potentials must all be given.

import { buildGrid, applyGeometry } from './grid.js';
import { FARADAY, GAS_CONSTANT } from './constants.js';
import { normalizeStatistics } from './statistics.js';
import { DeviceError } from './errors.js';
import { stoichiometry as equationStoichiometry, faceSides, parseEquation } from './equation.js';

export { DeviceError };

// Contacts use the same laws as internal faces: the outside is a phase with known levels.
const SPECIES_LINK_TYPES = new Set(['blocked', 'equilibrium', 'conductance', 'exchange']);
const INTERFACE_LINK_TYPES = new Set(['equilibrium', 'blocked', 'conductance', 'permeability']);
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

// A profile against the device's x: { x, values }, piecewise linear between its points and
// constant beyond its ends (read with profileAt).
function profile(v, path, what, ok, okText) {
  fields(v, path, ['x', 'values']);
  need(
    Array.isArray(v.x) && Array.isArray(v.values) && v.x.length > 0 && v.x.length === v.values.length,
    `${path} must be a profile { x, values }: two arrays of the same length, positions (m, the device's x) and ${what}`,
  );
  v.x.forEach((xk, k) => {
    need(isFiniteNumber(xk) && (k === 0 || xk > v.x[k - 1]), `${path}.x must be finite and strictly increasing, but x[${k}] is ${JSON.stringify(xk)}`);
  });
  v.values.forEach((y, k) => need(isFiniteNumber(y) && ok(y), `${path}.values[${k}] must be ${okText}, got ${JSON.stringify(y)}`));
  return { x: Float64Array.from(v.x), values: Float64Array.from(v.values) };
}

// SRH (trap-assisted) kinetics in place of mass action, for a reaction that consumes exactly one
// negative and one positive species and makes none, e⁻ + h⁺ = 0:
//   r = n p (1 − e^{−A/RT}) / (t_p (n + n₁) + t_n (p + p₁)),  p₁ = n p e^{−A/RT} / n₁,
// with t = τ (s) in the bulk, 1/v (m/s) at a face, and n₁ the negative species' concentration
// with its level at the trap (default n_i, a midgap trap). p₁ from the state keeps detailed
// balance exact, whatever lies between the two (a band offset, a φ jump).
function srhLaw(sdef, path, list, species, timeKeys) {
  const [kn, kp] = timeKeys;
  need(isObject(sdef), `${path} must be { ${kn}, ${kp}, n1 }`);
  fields(sdef, path, [kn, kp, 'n1']);
  const reactants = list.filter((p) => p.nu < 0);
  need(
    list.length === 2 && reactants.length === 2 && reactants.every((p) => p.nu === -1) && species[list[0].i].z * species[list[1].i].z < 0,
    `${path}: SRH kinetics needs a reaction consuming one negative and one positive species and making none, such as 'e- + h+ = 0'`,
  );
  const neg = species[list[0].i].z < 0 ? list[0] : list[1], pos = neg === list[0] ? list[1] : list[0];
  const time = (key) => {
    const v = positive(sdef[key], `${path}.${key}`);
    return key.startsWith('v') ? 1 / v : v; // a velocity's reciprocal is the time it stands for
  };
  return { n: neg, p: pos, tn: time(kn), tp: time(kp), n1: sdef.n1 === undefined ? NaN : positive(sdef.n1, `${path}.n1`) };
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
  fields(def, 'device', ['T', 'species', 'materials', 'regions', 'interfaces', 'bulkReactions', 'contacts', 'ports', 'grid', 'geometry']);

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
  // (an electrode's carrier can live outside the device: a reacting port's terminal species)
  const outside = (sp) => (def.ports ?? []).some((p) => p?.terminal === sp.name && Array.isArray(p.reactions) && p.reactions.length > 0);
  species.forEach((sp, i) => {
    const reacting = (def.ports ?? []).some((p) => Array.isArray(p?.reactions) && p.reactions.length > 0);
    need(
      materials.some((m) => m.present[i]) || outside(sp),
      `species '${sp.name}' is not present in any material` + (reacting ? " (an electrode's carrier can instead be a reacting port's terminal)" : ''),
    );
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
    // Initial composition: the starting state, and the conserved amount of any spectator. Each
    // entry is a number or a profile { x, values } against the device's x.
    const c0 = new Float64Array(nSpecies).fill(NaN);
    const c0Profile = new Array(nSpecies).fill(null);
    if (reg.c0 !== undefined) {
      need(isObject(reg.c0), `${path}.c0 must be an object mapping species names to concentrations`);
      for (const [sname, v] of Object.entries(reg.c0)) {
        const at = `${path}.c0.${sname}`;
        need(speciesIndex.has(sname), `${at}: unknown species '${sname}'${known(speciesIndex)}`);
        const i = speciesIndex.get(sname);
        need(mat.present[i], `${at}: '${sname}' is absent from material '${mat.name}'`);
        const none = `; to have none of '${sname}' here, leave it out of material '${mat.name}' (a region of another material without it, if it's present elsewhere)`;
        if (isObject(v)) {
          c0Profile[i] = profile(v, at, 'concentrations (mol/m³)', (ck) => ck > 0, `a concentration > 0 (mol/m³)${none}`);
          continue;
        }
        need(isFiniteNumber(v) && v > 0, `${at} must be a concentration > 0 (mol/m³) or a profile { x, values }, got ${JSON.stringify(v)}${none}`);
        c0[i] = v;
      }
    }
    return {
      name: reg.name ?? `region ${r}`,
      material: materialIndex.get(reg.material),
      length,
      fixedCharge,
      c0,
      c0Profile,
      background,
      velocity,
      mixing,
      // With no grid options at all, each region is graded toward its ends: fine where profiles
      // are steep (double layers, depleted electrodes), coarse inside.
      grid: reg.grid ?? (def.grid === undefined ? { hmin: length / 1000, hmax: length / 20 } : undefined),
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
  const ports = (def.ports ?? []).map((pdef, k) => normalizePort(pdef, `ports[${k}]`, regions, materials, species, speciesIndex, RT));

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
  checkLevel(terminals.map((t) => t.drive), terminals.map((t) => (t.kind === 'port' ? ports[t.index] : contacts[t.side])));

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
  const geometry = normalizeGeometry(def.geometry);
  if (geometry.type !== 'planar') {
    regions.forEach((reg, r) => need(reg.velocity === 0, `regions[${r}].velocity: flow is for a planar device (a uniform velocity through a varying cross-section wouldn't conserve the liquid)`));
  }
  applyGeometry(grid, geometry);
  // (zero only at an end, as at a sphere's centre: nothing passes there)
  const inside = (x) => x > 0 && x < grid.length;
  need(
    grid.area.every((a, g) => a > 0 || g === 0 || g === grid.nNodes - 1) && (geometry.type !== 'profile' || geometry.values.every((a, k) => a > 0 || !inside(geometry.x[k]))),
    'geometry: the cross-section vanishes inside the device; only an end may have A = 0',
  );

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
  // A node's surface belongs to one electrode.
  const surfaced = new Int32Array(grid.nNodes).fill(-1);
  ports.forEach((port, k) => {
    if (port.surface.length === 0) return;
    for (const g of port.nodes) {
      need(surfaced[g] < 0, `ports[${k}]: its window overlaps ports[${surfaced[g]}]'s, and both have a surface; a spot of metal has one surface`);
      surfaced[g] = k;
    }
  });

  return {
    T, RT, F: FARADAY, species, speciesIndex, materials, materialIndex, regions, interfaces, reactions, contacts, terminals, grid, warnings, ports, geometry,
  };
}

// Bulk reaction  Σ ν_R R ⇌ Σ ν_P P  with mass-action rate r = k_f Π c_R^ν (1 − e^{−A/RT}), where
// the affinity A = Σ_R ν μ̄ − Σ_P ν μ̄ includes any fixed-activity neutral participants (e.g.
// H₂O, given by its μ). This equals k_f Π c_R − k_b Π c_P with k_b fixed by the standard
// potentials, so equilibrium is exactly A = 0. kf maps material names to rate constants (the
// reaction runs only in those materials).
function normalizeReaction(rdef, path, species, speciesIndex, materials, materialIndex, regions, RT) {
  need(isObject(rdef), `${path} must be an object`);
  fields(rdef, path, ['nu', 'equation', 'fixed', 'kf', 'srh']);
  need((rdef.nu === undefined) !== (rdef.equation === undefined), `${path}: give the reaction as nu or as an equation, one of them`);
  const nu = rdef.nu ?? equationStoichiometry(rdef.equation, `${path}.equation`);
  const st = stoichiometry(nu, rdef.equation === undefined ? `${path}.nu` : `${path}.equation`, rdef.fixed, `${path}.fixed`, species, speciesIndex, RT);
  need(st.list.length > 0, `${path}: no mobile participants`);
  need(st.charge === 0, `${path}: charge is not balanced (Σ ν z = ${st.charge})`);
  // SRH kinetics instead: srh maps material names to { tauN, tauP, n1 }.
  need((rdef.kf === undefined) !== (rdef.srh === undefined), `${path}: give kf (mass action) or srh (trap-assisted), one of them`);
  if (rdef.srh !== undefined) {
    need(isObject(rdef.srh), `${path}.srh must map material names to { tauN, tauP, n1 }`);
    const srh = new Array(materials.length).fill(null), kf = new Float64Array(materials.length), kfProfile = new Array(materials.length).fill(null);
    for (const [mname, v] of Object.entries(rdef.srh)) {
      need(materialIndex.has(mname), `${path}.srh.${mname}: unknown material`);
      const m = materialIndex.get(mname);
      srh[m] = srhLaw(v, `${path}.srh.${mname}`, st.list, species, ['tauN', 'tauP']);
      for (const { i } of st.list) need(materials[m].present[i], `${path}.srh.${mname}: '${species[i].name}' is absent from material '${mname}'`);
      kf[m] = 1; // (runs here)
    }
    const reactants = st.list.map(({ i, nu }) => ({ i, nu: -nu }));
    return { reactants, products: [], fixedA: st.fixedA, kf, kfProfile, srh, generation: false };
  }
  need(isObject(rdef.kf), `${path}.kf must map material names to forward rate constants`);
  // Each a number, or a profile { x, values } against the device's x (e.g. absorption); kf[m]
  // is then its largest value, which says whether the reaction runs in that material at all.
  const kf = new Float64Array(materials.length), kfProfile = new Array(materials.length).fill(null);
  for (const [mname, v] of Object.entries(rdef.kf)) {
    need(materialIndex.has(mname), `${path}.kf.${mname}: unknown material`);
    const m = materialIndex.get(mname);
    if (isObject(v)) {
      kfProfile[m] = profile(v, `${path}.kf.${mname}`, 'rate constants', (k) => k >= 0, 'a rate constant ≥ 0');
      kf[m] = Math.max(...kfProfile[m].values);
    } else kf[m] = nonNegative(v, `${path}.kf.${mname}`);
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
    kfProfile,
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
      need(
        fixed[name] !== undefined,
        `${path}.${name}: not a species, so give its μ in ${fixedPath} (fixed-activity participants are neutral)` +
          (/^\d+\D/.test(name) ? `; if ${name.match(/^\d+/)[0]} is a coefficient, put a space after it` : ''),
      );
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
// Or saturating (a transporter or an enzyme that turns over at most vmax per area): with K a
// half-saturation concentration for each species the forward reaction consumes,
//   r = vmax Π_{ν<0} (c/(c + K))^{|ν|} (1 − e^{−a}),
// Michaelis–Menten in each substrate, and still exactly zero at A = 0 (the backward rate is the
// one detailed balance implies).
function normalizeFaceReactions(idef, where, matL, matR, species, speciesIndex, RT) {
  need(idef.reactions === undefined || Array.isArray(idef.reactions), `${where}.reactions must be an array`);
  return (idef.reactions ?? []).map((rdef, k) => {
    const rpath = `${where}.reactions[${k}]`;
    need(isObject(rdef), `${rpath} must be { equation, fixed, k0, alpha } or { left, right, fixed, k0, alpha }`);
    fields(rdef, rpath, ['equation', 'left', 'right', 'fixed', 'k0', 'alpha', 'srh', 'vmax', 'K']);
    let sides = rdef;
    if (rdef.equation !== undefined) {
      need(rdef.left === undefined && rdef.right === undefined, `${rpath}: give the reaction as an equation or as left and right, not both`);
      const holds = (mat) => (name) => speciesIndex.has(name) && Boolean(mat.present[speciesIndex.get(name)]);
      sides = faceSides(parseEquation(rdef.equation, `${rpath}.equation`), { name: matL.name, holds: holds(matL), conductor: !!matL.conductor }, { name: matR.name, holds: holds(matR), conductor: !!matR.conductor }, rdef.fixed, `${rpath}.equation`);
    }
    const part = [];
    let fixedA = 0, charge = 0;
    for (const [key, mat] of [['left', matL], ['right', matR]]) {
      const st = stoichiometry(sides[key], rdef.equation === undefined ? `${rpath}.${key}` : `${rpath}.equation`, rdef.fixed, `${rpath}.fixed`, species, speciesIndex, RT);
      for (const { i, nu } of st.list) {
        need(mat.present[i], `${rpath}.${key}.${species[i].name}: absent from '${mat.name}'`);
        part.push({ i, nu, side: key === 'left' ? 0 : 1 });
      }
      fixedA += st.fixedA;
      charge += st.charge;
    }
    need(part.length > 0, `${rpath}: no species participate`);
    need(charge === 0, `${rpath}: charge is not balanced (Σ ν z = ${charge})`);
    if (rdef.vmax !== undefined || rdef.K !== undefined) {
      need(rdef.k0 === undefined && rdef.alpha === undefined && rdef.srh === undefined, `${rpath}: give vmax and K (saturating), srh, or k0 and alpha (Butler–Volmer), one of them`);
      const vmax = positive(rdef.vmax, `${rpath}.vmax (mol/(m²·s))`);
      need(isObject(rdef.K), `${rpath}.K must map each species the forward reaction consumes to its half-saturation concentration (mol/m³)`);
      for (const name of Object.keys(rdef.K)) {
        need(part.some((p) => p.nu < 0 && species[p.i].name === name), `${rpath}.K.${name}: not a species the forward reaction consumes`);
        positive(rdef.K[name], `${rpath}.K.${name} (mol/m³)`);
      }
      for (const p of part) {
        if (p.nu >= 0) continue;
        const mat = p.side ? matR : matL;
        need(!mat.conductor, `${rpath}: saturating kinetics is for dissolved substrates, not a metal's ${species[p.i].name}`);
        need(rdef.K[species[p.i].name] !== undefined, `${rpath}.K.${species[p.i].name}: give a half-saturation concentration for every species the forward reaction consumes`);
      }
      return { part, fixedA, vmax, K: part.map((p) => (p.nu < 0 ? rdef.K[species[p.i].name] : 0)) };
    }
    if (rdef.srh !== undefined) {
      need(rdef.k0 === undefined && rdef.alpha === undefined, `${rpath}: give srh or k0 (and alpha), not both`);
      const srh = srhLaw(rdef.srh, `${rpath}.srh`, part, species, ['vn', 'vp']);
      for (const p of part) {
        const mat = p.side ? matR : matL;
        need(!mat.conductor, `${rpath}.srh: SRH kinetics is for carriers in semiconductors, not a metal's ${species[p.i].name}`);
      }
      return { part, fixedA, srh };
    }
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
    need(k === 0 || t >= v.t[k - 1], `${path}.t must not decrease`);
    need(k < 2 || t > v.t[k - 2], `${path}.t: at most two points at one time (a step)`);
  });
  v.values.forEach((x, k) => finite(x, `${path}.values[${k}]`));
  need(v.repeat === undefined || typeof v.repeat === 'boolean', `${path}.repeat must be true or false`);
  need(!v.repeat || v.t[v.t.length - 1] > v.t[0], `${path}: a repeating waveform needs a period (its last time after its first)`);
  return { t: Float64Array.from(v.t), values: Float64Array.from(v.values), repeat: v.repeat === true };
}

// Something must hold the device's overall level: a terminal held at a voltage, and if any is
// driven by a current, one held that passes current (one that passes none, a closed end, ties the
// level of nothing, and the driven terminal's voltage would float with φ).
function checkLevel(drives, owners) {
  need(drives.some((d) => d.kind === 'V'), "every terminal is driven by a current, so the device's overall level floats: hold at least one at a voltage V");
  const driven = drives.findIndex((d) => d.kind === 'I');
  if (driven < 0) return;
  const name = (k) => (k < 2 ? `contacts.${k === 0 ? 'left' : 'right'}` : `ports[${k - 2}]`);
  need(
    drives.some((d, k) => d.kind === 'V' && owners[k].passes),
    `${name(driven)} is driven by a current, but no terminal held at a voltage passes any, so nothing fixes the level its voltage is read against: hold one that does (the port of an electrode at V: 0, say: at the only electrode of a device, the same as I: 0)`,
  );
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
  model.ports.forEach((port, k) => need(port.passes || drives[2 + k].kind === 'V', `ports[${k}].I: this port passes no current (it exchanges only neutral species), so give it no drive`));
  checkLevel(drives, [model.contacts.left, model.contacts.right, ...model.ports]);
  return drives;
}

/**
 * A source's value at time t (s). At a step (two points at one time) it's the value after the
 * step; with `before`, the limit from earlier times, as a time step ending at t sees it.
 */
export function sourceAt(src, t, before = false) {
  if (src.value !== undefined) return src.value;
  const { t: ts, values: vs, repeat } = src, n = ts.length;
  if (before) t -= 8 * Number.EPSILON * Math.max(Math.abs(t), ts[n - 1] - ts[0], 1e-300); // just before: a few ulps
  if (repeat) {
    const T = ts[n - 1] - ts[0];
    t = ts[0] + ((((t - ts[0]) % T) + T) % T);
  }
  if (t < ts[0]) return vs[0];
  if (t >= ts[n - 1]) return vs[n - 1];
  let k = 0;
  while (ts[k] < t) k++;
  if (ts[k] === t) {
    while (k + 1 < n && ts[k + 1] === t) k++; // after a step
    return vs[k];
  }
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

// The device's cross-section A(x): planar (the default, A = 1 m², so currents read as densities),
// spherical or cylindrical about a centre r0 to the left of x = 0 (r = r0 + x), or a profile
// { x, values } of A in m². Every flux, current and amount is then a total through A.
function normalizeGeometry(g) {
  if (g === undefined || g === 'planar') return { type: 'planar' };
  if (g === 'spherical' || g === 'cylindrical') g = { type: g };
  need(isObject(g), "geometry must be 'planar', 'spherical', 'cylindrical', { type, r0 } or { area: { x, values } }");
  if (g.area !== undefined) {
    fields(g, 'geometry', ['area']);
    const p = profile(g.area, 'geometry.area', 'cross-sections (m²)', (a) => a >= 0, 'a cross-section ≥ 0 (m²)');
    need([...p.values].some((a) => a > 0), 'geometry.area: the cross-section is zero everywhere');
    return { type: 'profile', x: p.x, values: p.values };
  }
  fields(g, 'geometry', ['type', 'r0']);
  need(['planar', 'spherical', 'cylindrical'].includes(g.type), "geometry.type must be 'planar', 'spherical' or 'cylindrical'");
  if (g.type === 'planar') return { type: 'planar' };
  const r0 = g.r0 === undefined ? 0 : nonNegative(g.r0, 'geometry.r0 (m, the radius at x = 0)');
  return { type: g.type, r0 };
}

// An internal port: an outside phase with known levels (V_i = V + offset_i, or μ for neutral
// species, as at a contact), exchanging with every node in a window of one region.
// Links per species: 'equilibrium' (μ̄ held at the outside level throughout the window),
// { type: 'conductance', G } for charged species (G in S/m³: a source G (V_out − V_i)/(zF) per
// volume), { type: 'exchange', k, mu } for neutral ones (k in mol/(m³·s): a source
// k (μ_out − μ)/RT per volume), or 'blocked' (the default). And reactions with the port's
// terminal species, a metal's carrier at the port's level: an electrode surface spread through
// the window, `area` (m²/m³) of it per volume.
function normalizePort(pdef, path, regions, materials, species, speciesIndex, RT) {
  need(isObject(pdef), `${path} must be an object`);
  fields(pdef, path, ['name', 'region', 'from', 'to', 'V', 'I', 'R', 'terminal', 'species', 'reactions', 'area', 'surface', 'capacitance']);
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
  const surface = portSurface(pdef, path, speciesIndex);
  const reactions = portReactions(pdef, path, mat, terminal, species, speciesIndex, RT, surface);
  surface.forEach((sp, s) => need(reactions.some((rx) => rx.part.some((p) => p.side === 2 && p.s === s)), `${path}.surface.${sp.name}: no reaction of the port makes or uses it`));
  const capacitance = portCapacitance(pdef, path, mat);
  need(pdef.species === undefined ? reactions.length > 0 || capacitance !== null : isObject(pdef.species), `${path}.species must map species names to port links`);
  const links = species.map(() => ({ type: 'blocked' }));
  for (const [sname, raw] of Object.entries(pdef.species ?? {})) {
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
  need(links.some((l) => l.type !== 'blocked') || reactions.length > 0 || capacitance !== null, `${path}.species: the port exchanges no species`);
  // One that exchanges only neutral species (an O₂ supply) carries no current, so its voltage is
  // nobody's business: it can't be driven by one.
  const passes = reactions.length > 0 || capacitance !== null || links.some((l, i) => l.type !== 'blocked' && species[i].z !== 0);
  need(passes || drive.kind === 'V', `${path}.I: this port passes no current (it exchanges only neutral species), so give it no drive`);
  let area = null;
  if (reactions.length > 0 || capacitance !== null) {
    need(pdef.area !== undefined, `${path}.area: give the electrode's area per volume of the window (m²/m³, e.g. 1/h for a film of thickness h on it), a number or a profile`);
    area = typeof pdef.area === 'number' ? { value: positive(pdef.area, `${path}.area (m²/m³)`) } : profile(pdef.area, `${path}.area`, 'areas per volume (m²/m³)', (a) => a >= 0, 'an area per volume ≥ 0');
  } else need(pdef.area === undefined, `${path}.area: only an electrode (a port with reactions or a capacitance) has an area`);
  need(surface.length === 0 || reactions.length > 0, `${path}.surface: only an electrode (a port with reactions) has a surface`);
  return { name: pdef.name ?? path, region: r, from, to, span: reg.length, drive, terminal, species: links, passes, reactions, area, surface, capacitance };
}

// A capacitance spread through the window: a gate along a channel, an electrode's double layer.
// Per area of electrode it holds σ = C (V − zeroCharge − φ) on the port's side, and the window
// the opposite: the region's charge balance (Gauss's law, or neutrality at ε = 0) gains aσ.
function portCapacitance(pdef, path, mat) {
  if (pdef.capacitance === undefined) return null;
  const cpath = `${path}.capacitance`;
  need(isObject(pdef.capacitance), `${cpath} must be { C, zeroCharge }`);
  fields(pdef.capacitance, cpath, ['C', 'zeroCharge']);
  need(!mat.conductor && !mat.phiFree, `${cpath}: a capacitance couples to φ, which a ${mat.conductor ? 'conductor' : 'region with no free charge'} doesn't have`);
  return { C: positive(pdef.capacitance.C, `${cpath}.C (F/m² of electrode)`), zeroCharge: pdef.capacitance.zeroCharge === undefined ? 0 : finite(pdef.capacitance.zeroCharge, `${cpath}.zeroCharge (V)`) };
}

// An electrode port's surface species: coverages θ of its sites (Langmuir: each takes one site,
// μ = μ° + RT ln(θ/θ₀) with θ₀ = 1 − Σθ the bare fraction), Γ mol of sites per m² of electrode.
function portSurface(pdef, path, speciesIndex) {
  if (pdef.surface === undefined) return [];
  need(isObject(pdef.surface), `${path}.surface must map surface species to { mu0, capacity, theta0 }`);
  const list = Object.entries(pdef.surface).map(([name, d]) => {
    const spath = `${path}.surface.${name}`;
    need(!speciesIndex.has(name), `${spath}: '${name}' is a species of the device; a surface species is the port's own, with its own name`);
    need(isObject(d), `${spath} must be { mu0, capacity, theta0 }`);
    fields(d, spath, ['mu0', 'capacity', 'theta0']);
    const theta0 = d.theta0 === undefined ? 1e-6 : positive(d.theta0, `${spath}.theta0`);
    return { name, mu0: finite(d.mu0, `${spath}.mu0 (J/mol)`), capacity: positive(d.capacity, `${spath}.capacity (mol of sites per m² of electrode)`), theta0 };
  });
  need(list.length > 0, `${path}.surface: give at least one surface species, or leave it out`);
  const caps = new Set(list.map((sp) => sp.capacity));
  need(caps.size === 1, `${path}.surface: the species share the electrode's sites, so give them one capacity`);
  need(list.reduce((t, sp) => t + sp.theta0, 0) < 1, `${path}.surface: the starting coverages (theta0) must sum to less than 1`);
  return list;
}

// A port's reactions: Butler–Volmer, as at a face, between species of the port's region and the
// port's terminal species (a metal's carrier, activity 1, at the port's level), per area of the
// electrode: part side 0 for the region's species, side 1 for the carrier.
function portReactions(pdef, path, mat, terminal, species, speciesIndex, RT, surface = []) {
  need(pdef.reactions === undefined || Array.isArray(pdef.reactions), `${path}.reactions must be an array`);
  const list = pdef.reactions ?? [];
  if (list.length === 0) return [];
  need(!mat.conductor, `${path}.reactions: a port on a conductor is a wire; reactions belong on a face, or on a port in the solution`);
  // The metal's countercharge (each spot's double layer) is below the grid, so the solution
  // beside it must be neutral, as in porous-electrode theory; with ε > 0 the reactions would
  // leave a net charge with nothing to balance it.
  need(mat.epsr === 0 || mat.phiFree, `${path}.reactions: an electrode spread through a window needs its region strictly neutral (ε = 0): its double layers are below the grid`);
  need(terminal !== null && species[terminal].z !== 0, `${path}.terminal: a port with reactions needs its electrode's carrier as its terminal species (e.g. 'e-')`);
  need(!mat.present[terminal], `${path}.terminal: '${species[terminal].name}' is in the port's region too, so a reaction can't tell the electrode's from the region's`);
  const carrier = species[terminal].name;
  return list.map((rdef, k) => {
    const rpath = `${path}.reactions[${k}]`;
    need(isObject(rdef) && typeof rdef.equation === 'string', `${rpath} must be { equation, fixed, k0, alpha }`);
    fields(rdef, rpath, ['equation', 'fixed', 'k0', 'alpha', 'bare']);
    const holds = (name) => speciesIndex.has(name) && Boolean(mat.present[speciesIndex.get(name)]);
    const onSurface = (name) => surface.findIndex((sp) => sp.name === name);
    const sides = faceSides(parseEquation(rdef.equation, `${rpath}.equation`), { name: mat.name, holds, conductor: false }, { name: 'the electrode', holds: (name) => name === carrier || onSurface(name) >= 0, conductor: true }, rdef.fixed, `${rpath}.equation`);
    const part = [];
    let fixedA = 0, charge = 0;
    // The electrode's side: its carrier, and its surface species (side 2, neutral, by index s).
    const electrode = {};
    for (const [name, nu] of Object.entries(sides.right)) {
      const s = onSurface(name);
      if (s >= 0) part.push({ s, nu, side: 2 });
      else electrode[name] = nu;
    }
    for (const [map, side] of [[sides.left, 0], [electrode, 1]]) {
      const st = stoichiometry(map, `${rpath}.equation`, rdef.fixed, `${rpath}.fixed`, species, speciesIndex, RT);
      for (const { i, nu } of st.list) part.push({ i, nu, side });
      fixedA += st.fixedA;
      charge += st.charge;
    }
    need(part.some((p) => p.side === 1), `${rpath}: no '${carrier}' from the electrode takes part`);
    need(part.some((p) => p.side === 0 || p.side === 2), `${rpath}: no species of the port's region or its surface takes part`);
    need(charge === 0, `${rpath}: charge is not balanced (Σ ν z = ${charge})`);
    need(rdef.bare === undefined || typeof rdef.bare === 'boolean', `${rpath}.bare must be true or false`);
    return { part, fixedA, k0: positive(rdef.k0, `${rpath}.k0 (mol/(m²·s))`), alpha: transferCoefficient(rdef.alpha, `${rpath}.alpha`), bare: rdef.bare === true };
  });
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
  for (const port of ports) if (port.reactions.length > 0 || port.capacitance || port.species.some((l, i) => l.type !== 'blocked' && species[i].z !== 0)) anchored.add(find(port.region));
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
  fields(idef, where, ['phi', 'dipole', 'step', 'sheetCharge', 'species', 'reactions', 'gates']);
  const bothNeutral = (matL.epsr === 0 && matR.epsr === 0) || matL.phiFree || matR.phiFree;
  const rawPhi = idef.phi ?? (bothNeutral ? 'neutral' : 'pinned');
  const phi = typeof rawPhi === 'string' ? { type: rawPhi } : { ...rawPhi };
  need(isObject(phi), `${where}.phi must be a law name or { type, C }`);
  fields(phi, `${where}.phi`, phi.type === 'capacitive' ? ['type', 'C'] : ['type']);
  need(['pinned', 'neutral', 'capacitive'].includes(phi.type), `${where}.phi must be 'pinned', 'neutral' or { type: 'capacitive', C }`);
  if (phi.type === 'capacitive') positive(phi.C, `${where}.phi.C`);
  // A pinned face is a capacitor of infinite C. With no field on either side (ε = 0), its
  // charge would sit as free excess ions in the edge cells at no cost, so nothing determines it.
  need(
    !(phi.type === 'pinned' && matL.epsr === 0 && matR.epsr === 0),
    `${where}.phi: 'pinned' between two ε = 0 materials leaves the face's charge undetermined (no field on either side to hold it); ` +
      "use 'neutral' (the default: a free Donnan jump), or { type: 'capacitive', C } for a charged layer such as a lipid bilayer",
  );
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
  // blocked, an ohmic interface conductance G (S/m²), or a permeability P (m/s): electrodiffusion
  // through a thin membrane in a constant field, Goldman–Hodgkin–Katz.
  const links = defaultInterfaceLinks(matL, matR, species);
  const gates = normalizeGates(idef.gates, where, matL, matR);
  if (idef.species !== undefined) {
    need(isObject(idef.species), `${where}.species must map species names to interface links`);
    for (const [sname, raw] of Object.entries(idef.species)) {
      const lpath = `${where}.species.${sname}`;
      need(speciesIndex.has(sname), `${lpath}: unknown species '${sname}'${known(speciesIndex)}`);
      const i = speciesIndex.get(sname);
      const link = typeof raw === 'string' ? { type: raw } : raw;
      need(isObject(link) && INTERFACE_LINK_TYPES.has(link.type), `${lpath}.type must be one of ${[...INTERFACE_LINK_TYPES].join(', ')}`);
      fields(link, lpath, link.type === 'conductance' ? ['type', 'G'] : link.type === 'permeability' ? ['type', 'P', 'gates'] : ['type']);
      if (link.type !== 'blocked') need(matL.present[i] && matR.present[i], `${lpath}: '${sname}' must be present on both sides`);
      if (link.type === 'conductance') {
        need(species[i].z !== 0, `${lpath}: a conductance link needs a charged species`);
        positive(link.G, `${lpath}.G`);
      }
      if (link.type === 'permeability') {
        positive(link.P, `${lpath}.P (m/s)`);
        need(matL.modelOf[i] < 0 && matR.modelOf[i] < 0, `${lpath}: a permeability link needs ideal (dilute) statistics for '${sname}' on both sides`);
      }
      // Gated: P × Π x_q^p over the face's gates (a channel open with probability m³h, say).
      let gated = [];
      if (link.gates !== undefined) {
        need(isObject(link.gates), `${lpath}.gates must map the face's gate names to exponents, e.g. { m: 3, h: 1 }`);
        gated = Object.entries(link.gates).map(([gname, p]) => {
          const q = gates.findIndex((gt) => gt.name === gname);
          need(q >= 0, `${lpath}.gates.${gname}: no such gate on this face${gates.length ? ` (it has ${gates.map((gt) => gt.name).join(', ')})` : ` (define it in ${where}.gates)`}`);
          need(Number.isInteger(p) && p >= 1 && p <= 8, `${lpath}.gates.${gname}: the exponent must be a whole number from 1 to 8, got ${JSON.stringify(p)}`);
          return [q, p];
        });
      }
      links[i] = { ...link, gates: gated };
    }
  }

  const reactions = normalizeFaceReactions(idef, where, matL, matR, species, speciesIndex, RT);
  checkReactingLinks(reactions, idef, where, matL, matR, species);
  return { phi, dipole, sheetCharge, links, reactions, gates, conductor: null };
}

// A face's gates (Hodgkin–Huxley): each a fraction x in [0, 1] with first-order kinetics in the
// voltage across the face, V = φ_right − φ_left,
//   dx/dt = α(V) (1 − x) − β(V) x,
// α and β each one of NeuroML's three standard forms, with rate in 1/s and midpoint and scale in
// volts: 'exp', rate·e^((V − midpoint)/scale); 'sigmoid', rate/(1 + e^((midpoint − V)/scale));
// 'expLinear', rate·y/(1 − e^(−y)) with y = (V − midpoint)/scale. They scale permeabilities
// (a link's gates), and nothing else: a gate holds no charge and exchanges no free energy.
const GATE_RATE_TYPES = ['exp', 'sigmoid', 'expLinear'];
function normalizeGates(raw, where, matL, matR) {
  if (raw === undefined) return [];
  need(isObject(raw), `${where}.gates must map gate names to { alpha, beta }`);
  need(!matL.phiFree && !matR.phiFree, `${where}.gates: a gate follows the voltage across the face, which needs φ on both sides`);
  return Object.entries(raw).map(([name, gdef]) => {
    const gpath = `${where}.gates.${name}`;
    need(isObject(gdef), `${gpath} must be { alpha, beta }`);
    fields(gdef, gpath, ['alpha', 'beta']);
    const rate = (r, rpath) => {
      need(isObject(r) && GATE_RATE_TYPES.includes(r.type), `${rpath} must be { type, rate, midpoint, scale } with type one of ${GATE_RATE_TYPES.join(', ')}`);
      fields(r, rpath, ['type', 'rate', 'midpoint', 'scale']);
      const scale = finite(r.scale, `${rpath}.scale (V)`);
      need(scale !== 0, `${rpath}.scale (V) must not be zero`);
      return { type: r.type, rate: positive(r.rate, `${rpath}.rate (1/s)`), midpoint: finite(r.midpoint, `${rpath}.midpoint (V)`), scale };
    };
    return { name, alpha: rate(gdef.alpha, `${gpath}.alpha`), beta: rate(gdef.beta, `${gpath}.beta`) };
  });
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
      need(link.type !== 'permeability', `${lpath}: a permeability link is for a membrane between two solutions, not at a conductor`);
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
  return { phi, dipole: 0, zeroCharge, sheetCharge: 0, links, reactions, gates: [], conductor: { side, i: metal.i } };
}

function defaultInterfaceLinks(matL, matR, species) {
  return species.map((_, i) => ({ type: matL.present[i] && matR.present[i] ? 'equilibrium' : 'blocked' }));
}

function transferCoefficient(v, path) {
  if (v === undefined) return 0.5;
  need(isFiniteNumber(v) && v >= 0 && v <= 1, `${path} must be between 0 and 1, got ${JSON.stringify(v)}`);
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
      phi = { ...raw, zeroCharge: finite(cdef.zeroCharge, `${path}.zeroCharge (V − φ_edge at zero charge: a pzc, a barrier, or a gate's flat-band voltage less the bulk φ)`) };
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
