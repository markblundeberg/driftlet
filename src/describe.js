// A readable summary of a device definition, for driftlet/kit, with warnings for numbers that
// look like unit slips (a D in cm²/s, a μ° in eV, a capacitance in µF/cm²) and for double layers
// the grid won't resolve. Heuristics only: a warning is a prompt to check, never an error.

import { normalizeDevice } from './device.js';
import { EPS0, FARADAY, GAS_CONSTANT, AVOGADRO } from './constants.js';

const isObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v);

// A length, time or other quantity with an SI prefix: 1.5e-8 → '15 n'.
function si(v, unit) {
  if (v === 0 || !Number.isFinite(v)) return `${v} ${unit}`;
  const prefixes = [[1e9, 'G'], [1e6, 'M'], [1e3, 'k'], [1, ''], [1e-3, 'm'], [1e-6, 'µ'], [1e-9, 'n'], [1e-12, 'p'], [1e-15, 'f']];
  const a = Math.abs(v);
  const [f, p] = prefixes.find(([f]) => a >= f * 0.9995) ?? prefixes[prefixes.length - 1];
  return `${+(v / f).toPrecision(3)} ${p}${unit}`;
}
const num = (v) => (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-2 && v !== 0) ? v.toExponential(2) : String(+v.toPrecision(3)));
const perCm3 = (c) => `${(c * AVOGADRO * 1e-6).toExponential(2)} cm⁻³`;

// An equation from signed stoichiometries, for display.
function equation(...maps) {
  const lhs = [], rhs = [];
  for (const m of maps) {
    for (const [name, nu] of Object.entries(m ?? {})) (nu < 0 ? lhs : rhs).push(`${Math.abs(nu) === 1 ? '' : `${Math.abs(nu)} `}${name}`);
  }
  return `${lhs.join(' + ') || '0'} = ${rhs.join(' + ') || '0'}`;
}

function drive(c) {
  const v = (x) => (isObject(x) ? 'a waveform' : `${num(x)}`);
  if (c?.I !== undefined) return `driven at I = ${v(c.I)} A/m²`;
  if (c?.R !== undefined) return `a source of ${v(c.V ?? 0)} V behind R = ${num(c.R)} Ω·m²`;
  return `held at ${v(c?.V ?? 0)} V`;
}

/**
 * Warnings for numbers that look like unit slips: SI is m, mol/m³, J/mol, m²/s, F/m², C/m³.
 * @param {object} def device definition
 * @returns {string[]}
 */
export function unitWarnings(def) {
  const out = [];
  const warn = (path, msg) => out.push(`${path}: ${msg}`);
  if (def.T !== undefined && def.T < 60) warn('T', `${def.T} K is very cold; T is in kelvin (if that's °C, add 273.15)`);
  if (def.T > 5000) warn('T', `${def.T} K is hotter than any device; T is in kelvin`);
  const conc = (path, c) => {
    if (c > 1e6) warn(path, `${num(c)} mol/m³ is over 1000 M; concentrations are mol/m³ (units.molar, and units.perCm3 for carriers)`);
  };
  (def.species ?? []).forEach((sp, i) => sp?.cRef !== undefined && conc(`species[${i}].cRef`, sp.cRef));
  for (const [mname, mat] of Object.entries(def.materials ?? {})) {
    if (!isObject(mat)) continue;
    const path = `materials.${mname}`;
    if (mat.epsr > 0 && mat.epsr < 1) warn(`${path}.epsr`, `${mat.epsr} is below 1; epsr is relative to vacuum`);
    if (mat.conductor?.conductivity > 1e9) warn(`${path}.conductor.conductivity`, `${num(mat.conductor.conductivity)} S/m is beyond any metal (silver is 6.3e7 S/m)`);
    for (const [sname, sp] of Object.entries(mat.species ?? {})) {
      if (!isObject(sp)) continue;
      const sp_ = `${path}.species.${sname}`;
      const z = (def.species ?? []).find((x) => x?.name === sname)?.z;
      const ion = z !== undefined && z !== 0 && sname !== 'e-' && sname !== 'h+';
      if (sp.D > 1) warn(`${sp_}.D`, `${num(sp.D)} m²/s is beyond any real diffusivity; D is m²/s (units.cm2PerS)`);
      else if (ion && sp.D > 1e-7) warn(`${sp_}.D`, `${num(sp.D)} m²/s is fast for an ion (H⁺ in water is 9.3e-9 m²/s); D is m²/s, so cm²/s needs units.cm2PerS`);
      if (sp.mu0 !== 0 && Math.abs(sp.mu0) < 2000) warn(`${sp_}.mu0`, `${sp.mu0} J/mol is small for a standard potential; mu0 is J/mol (kJ/mol × 1000, units.eV for per-particle energies, × F for volts)`);
      if (Math.abs(sp.mu0) > 1e7) warn(`${sp_}.mu0`, `${num(sp.mu0)} J/mol is over 100 eV per particle; mu0 is J/mol`);
      if (sp.cRef !== undefined) conc(`${sp_}.cRef`, sp.cRef);
    }
  }
  const regionX0 = [0];
  for (const reg of def.regions ?? []) regionX0.push(regionX0.at(-1) + (isObject(reg) && reg.length > 0 ? reg.length : 0));
  (def.regions ?? []).forEach((reg, r) => {
    if (!isObject(reg)) return;
    const path = `regions[${r}]`;
    if (reg.length > 1) warn(`${path}.length`, `${num(reg.length)} m is over a metre; lengths are m (units.um, units.nm)`);
    if (reg.length > 0 && reg.length < 1e-10) warn(`${path}.length`, `${num(reg.length)} m is less than an atom; lengths are m`);
    if (Math.abs(reg.fixedCharge) > 1e10) warn(`${path}.fixedCharge`, `${num(reg.fixedCharge)} C/m³ is over 6e22 cm⁻³ of charge; fixedCharge is C/m³ (units.perCm3(N) * FARADAY)`);
    for (const [sname, c] of Object.entries(reg.c0 ?? {})) {
      if (!isObject(c)) {
        conc(`${path}.c0.${sname}`, c);
        continue;
      }
      if (Array.isArray(c.c)) conc(`${path}.c0.${sname}.c`, Math.max(...c.c));
      // A profile is against the device's x, not the region's own.
      const x0 = regionX0[r], x1 = x0 + reg.length;
      if (Array.isArray(c.x) && c.x.length > 0 && (c.x.at(-1) < x0 || c.x[0] > x1)) {
        warn(`${path}.c0.${sname}.x`, `the profile (x from ${num(c.x[0])} to ${num(c.x.at(-1))} m) misses the region (${num(x0)} to ${num(x1)} m); x is the device's, in m`);
      }
    }
  });
  const capacitance = (path, phi) => {
    if (isObject(phi) && phi.C > 10) warn(`${path}.phi.C`, `${num(phi.C)} F/m² is over 1000 µF/cm²; C is F/m² (1 µF/cm² = 0.01 F/m²)`);
  };
  const kinetics = (path, reactions) => {
    (reactions ?? []).forEach((rx, k) => {
      if (rx?.k0 > 10) warn(`${path}.reactions[${k}].k0`, `${num(rx.k0)} mol/(m²·s) is an exchange current of ${num(rx.k0 * FARADAY)} A/m² per electron; k0 is mol/(m²·s)`);
    });
  };
  (def.interfaces ?? []).forEach((f, k) => {
    if (!isObject(f)) return;
    capacitance(`interfaces[${k}]`, f.phi);
    kinetics(`interfaces[${k}]`, f.reactions);
  });
  for (const side of ['left', 'right']) {
    const c = def.contacts?.[side];
    if (!isObject(c)) continue;
    capacitance(`contacts.${side}`, c.phi);
    for (const [sname, v] of Object.entries(c.bath?.c ?? {})) conc(`contacts.${side}.bath.c.${sname}`, v);
    if (typeof c.V === 'number' && Math.abs(c.V) > 50) warn(`contacts.${side}.V`, `${c.V} V is a very large bias; voltages are V`);
  }
  return out;
}

/**
 * A readable summary of a device: its regions (with doping, Debye lengths and the grid at their
 * ends), faces, contacts and reactions, its characteristic times, and warnings (unit slips, and
 * double layers coarser than the grid). It validates the definition first, so an invalid one
 * throws its DeviceError.
 * @param {object} def device definition
 * @returns {string}
 */
export function describe(def) {
  const model = normalizeDevice(def);
  const { species, materials, regions, grid } = model;
  const RT = GAS_CONSTANT * model.T;
  const lines = [];
  const warnings = unitWarnings(def);
  const total = grid.x[grid.x.length - 1] - grid.x[0];
  lines.push(`${regions.length} region${regions.length > 1 ? 's' : ''}, ${si(total, 'm')}, ${grid.x.length} nodes, T = ${model.T} K`);
  lines.push(`species: ${species.map((sp) => `${sp.name} (z = ${sp.z > 0 ? '+' : ''}${sp.z})`).join(', ')}`);

  lines.push('regions:');
  const times = [];
  regions.forEach((reg, r) => {
    const mat = materials[reg.material];
    const a = grid.regionStart[r], b = grid.regionEnd[r];
    const head = `  ${reg.name} [${mat.name}] ${si(reg.length, 'm')}`;
    if (mat.conductor) {
      lines.push(`${head}: conductor, ${species[mat.conductor.i].name} at ${num(mat.conductor.sigma)} S/m`);
      return;
    }
    const parts = [mat.epsr > 0 ? `εr ${mat.epsr}` : 'strictly neutral (ε = 0)'];
    if (reg.fixedCharge !== 0) parts.push(`fixed charge ${reg.fixedCharge > 0 ? 'donor-like' : 'acceptor-like'} ${perCm3(Math.abs(reg.fixedCharge) / FARADAY)}`);
    // The screening concentration: Σ z²c from the initial composition, else the fixed charge's.
    let zzc = 0;
    species.forEach((sp, i) => {
      if (!mat.present[i]) return;
      if (Number.isFinite(reg.c0[i])) zzc += sp.z * sp.z * reg.c0[i];
      else if (reg.c0Profile[i]) zzc += sp.z * sp.z * Math.max(...reg.c0Profile[i].c);
    });
    if (zzc === 0) zzc = Math.abs(reg.fixedCharge) / FARADAY;
    const Dmax = Math.max(0, ...species.map((sp, i) => (mat.present[i] && sp.z !== 0 ? mat.D[i] : 0)));
    const Dmin = Math.min(...species.map((sp, i) => (mat.present[i] && mat.D[i] > 0 ? mat.D[i] : Infinity)));
    if (mat.epsr > 0 && zzc > 0) {
      const lambda = Math.sqrt((mat.epsr * EPS0 * RT) / (FARADAY * FARADAY * zzc));
      // Double layers form at faces that aren't neutral, and at gate or pinned contacts.
      const pins = (side) => ['capacitive', 'pinned'].includes(model.contacts[side].phi.type);
      const atLeft = r === 0 ? pins('left') : model.interfaces[r - 1].phi.type !== 'neutral';
      const atRight = r === regions.length - 1 ? pins('right') : model.interfaces[r].phi.type !== 'neutral';
      const hEnd = Math.max(atLeft ? grid.segLength[a] : 0, atRight ? grid.segLength[b - 1] : 0);
      parts.push(`Debye length ${si(lambda, 'm')}${hEnd > 0 ? ` (end cells ${si(hEnd, 'm')})` : ''}`);
      if (lambda < 2e-10) {
        warnings.push(`${reg.name}: a Debye length of ${si(lambda, 'm')} is below atomic size, where a continuum double layer means little; consider epsr: 0 (strict neutrality)`);
      } else if (hEnd > lambda) {
        warnings.push(`${reg.name}: end cells of ${si(hEnd, 'm')} are coarser than the Debye length, ${si(lambda, 'm')}; a double layer there won't be resolved (grid.hmin)`);
      }
      if (Dmax > 0) times.push(`  dielectric relaxation in ${reg.name}: ~${si((lambda * lambda) / Dmax, 's')}`);
    }
    if (Number.isFinite(Dmin)) times.push(`  diffusion across ${reg.name}: L²/D up to ${si((reg.length * reg.length) / Dmin, 's')}`);
    lines.push(`${head}: ${parts.join(', ')}; ${b - a} cells`);
  });

  const idefs = def.interfaces ?? [];
  if (regions.length > 1) {
    lines.push('faces:');
    model.interfaces.forEach((itf, f) => {
      const L = regions[f], R = regions[f + 1];
      const idef = isObject(idefs[f]) ? idefs[f] : {};
      if (L.material === R.material && itf.phi.type === 'pinned' && itf.dipole === 0 && !(idef.reactions ?? []).length && idef.species === undefined) {
        lines.push(`  ${L.name} | ${R.name}: same material`);
        return;
      }
      const law = itf.phi.type === 'capacitive' ? `capacitive, C = ${num(itf.phi.C)} F/m²` : itf.phi.type;
      const align = itf.conductor ? (itf.phi.type === 'capacitive' ? `, zeroCharge ${num(itf.zeroCharge)} V` : '') : itf.phi.type === 'neutral' ? '' : `, dipole ${num(itf.dipole)} V`;
      const blocked = species.filter((sp, i) => itf.links[i].type === 'blocked' && materials[L.material].present[i] && materials[R.material].present[i]).map((sp) => sp.name);
      const rx = (idef.reactions ?? []).map((r) => r.equation ?? equation(r.left, r.right));
      lines.push(`  ${L.name} | ${R.name}: ${law}${align}${blocked.length ? `; blocked: ${blocked.join(', ')}` : ''}${rx.length ? `; reactions: ${rx.join('; ')}` : ''}`);
    });
  }

  lines.push('contacts:');
  for (const side of ['left', 'right']) {
    const cdef = def.contacts?.[side];
    const ct = model.contacts[side];
    const linked = species.filter((sp, i) => ct.species[i].type !== 'blocked').map((sp) => sp.name);
    const term = ct.terminal === null ? '' : ` (terminal ${species[ct.terminal].name})`;
    lines.push(`  ${side}: ${drive(cdef)}${term}; ${linked.length ? `exchanges ${linked.join(', ')}` : 'exchanges nothing'}; φ ${ct.phi.type}`);
  }
  (def.ports ?? []).forEach((p, k) => lines.push(`  port ${model.ports[k].name}: ${drive(p)}, in ${regions[model.ports[k].region].name}`));
  if ((def.bulkReactions ?? []).length) {
    lines.push('bulk reactions:');
    for (const rx of def.bulkReactions) lines.push(`  ${rx.equation ?? equation(rx.nu)} in ${Object.keys(rx.kf ?? {}).join(', ')}`);
  }
  if (times.length) lines.push('time scales:', ...times);
  if (warnings.length) lines.push('warnings:', ...warnings.map((w) => `  ${w}`));
  return lines.join('\n');
}
