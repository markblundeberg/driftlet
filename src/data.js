// A small data library for driftlet/kit, every entry with its source (see docs/data.md). Kept
// small on purpose: a wrong library number is worse than a wrong demo number. Bulk properties
// only. How levels line up at an interface, and a surface's work function, belong to that
// interface or surface and are never in here.

import { DeviceError, checkOptions, shown } from './errors.js';
import { FARADAY, GAS_CONSTANT, AVOGADRO } from './constants.js';
import { powr } from './pow.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const freeze = (table) => Object.freeze(Object.fromEntries(Object.entries(table).map(([k, v]) => [k, Object.freeze(v)])));

/**
 * Aqueous ions at 25 °C. D: diffusivity at infinite dilution, m²/s (CRC Handbook of Chemistry
 * and Physics, "Ionic conductivity and diffusion at infinite dilution"). mu0: standard Gibbs
 * energy of formation, J/mol, the usual table convention with ΔfG°(H⁺) = 0 (Wagman et al., "The
 * NBS tables of chemical thermodynamic properties", J. Phys. Chem. Ref. Data 11, Suppl. 2, 1982),
 * on the 1 mol/kg standard state, used here with c_ref = 1000 mol/m³ (1 M).
 * Names are the formula followed by the charge.
 */
export const IONS = freeze({
  'H+': { z: 1, D: 9.311e-9, mu0: 0 },
  'Li+': { z: 1, D: 1.029e-9, mu0: -293.31e3 },
  'Na+': { z: 1, D: 1.334e-9, mu0: -261.905e3 },
  'K+': { z: 1, D: 1.957e-9, mu0: -283.27e3 },
  'Ag+': { z: 1, D: 1.648e-9, mu0: 77.107e3 },
  'Cu2+': { z: 2, D: 0.714e-9, mu0: 65.49e3 },
  'Zn2+': { z: 2, D: 0.703e-9, mu0: -147.06e3 },
  'Fe2+': { z: 2, D: 0.719e-9, mu0: -78.9e3 },
  'Fe3+': { z: 3, D: 0.604e-9, mu0: -4.7e3 },
  'OH-': { z: -1, D: 5.273e-9, mu0: -157.244e3 },
  'Cl-': { z: -1, D: 2.032e-9, mu0: -131.228e3 },
  'NO3-': { z: -1, D: 1.902e-9, mu0: -111.25e3 },
  'SO42-': { z: -2, D: 1.065e-9, mu0: -744.53e3 },
});

/** μ° of liquid water at 25 °C, J/mol (Wagman et al. 1982), for reactions such as H₂O ⇌ H⁺ + OH⁻. */
export const H2O = -237.129e3;

/** Relative permittivity of water at 25 °C (Malmberg and Maryott, J. Res. NBS 56, 1, 1956). */
export const WATER_EPSR = 78.3;

/**
 * Semiconductors at 300 K (Sze, Physics of Semiconductor Devices, 2nd ed., 1981, appendix):
 * relative permittivity, band gap Eg (eV), electron affinity chi (eV), effective densities of
 * states Nc and Nv (cm⁻³), and drift mobilities of lightly doped material mun and mup
 * (cm²/(V·s)).
 */
export const SEMICONDUCTORS = freeze({
  Si: { epsr: 11.9, Eg: 1.12, chi: 4.05, Nc: 2.8e19, Nv: 1.04e19, mun: 1500, mup: 450 },
  Ge: { epsr: 16.0, Eg: 0.66, chi: 4.0, Nc: 1.04e19, Nv: 6.0e18, mun: 3900, mup: 1900 },
  GaAs: { epsr: 13.1, Eg: 1.424, chi: 4.07, Nc: 4.7e17, Nv: 7.0e18, mun: 8500, mup: 400 },
});

/**
 * Metals at 293 K: resistivity rho (10⁻⁸ Ω·m, CRC Handbook, "Electrical resistivity of pure
 * metals") and the conductivity 1/rho (S/m).
 */
export const METALS = freeze(
  Object.fromEntries(
    Object.entries({ Ag: 1.587, Cu: 1.678, Au: 2.214, Al: 2.65, Zn: 5.9, Li: 9.28, Pt: 10.5 }).map(([name, rho]) => [name, { rho, conductivity: 1e8 / rho }]),
  ),
);

/**
 * Water with the named ions, as a library piece for `build()`: each ion as a species
 * (c_ref = 1 M) and the material `water` holding them, with D and μ° from `IONS`. The data are
 * 25 °C's, and the material says so (T: 298.15): a device at another temperature warns.
 * @param {string[]} names ions from `IONS`
 * @param {{ epsr?: number, material?: string }} [opts] permittivity (default that of water;
 *   0 makes the solution strictly neutral, without double layers), and the material's name
 * @returns {{ species: object[], materials: Record<string, object> }}
 */
export function aqueous(names, { epsr = WATER_EPSR, material = 'water' } = {}) {
  checkOptions('aqueous', arguments[1], ['epsr', 'material']);
  if (!Array.isArray(names) || names.length === 0) fail('aqueous: give a list of ions, e.g. aqueous([\'Na+\', \'Cl-\'])');
  const species = [], entries = {};
  for (const name of names) {
    const ion = IONS[name];
    if (!ion) fail(`aqueous: no data for '${name}' (the library has ${Object.keys(IONS).join(', ')})`);
    species.push({ name, z: ion.z, cRef: 1000 });
    entries[name] = { D: ion.D, mu0: ion.mu0 };
  }
  // (25 °C data: a device at another temperature warns)
  return { species, materials: { [material]: { epsr, species: entries, T: 298.15 } } };
}

// How the semiconductors' data move with temperature, for semiconductor(name, { T }):
// - the gap by Varshni's form, Eg(T) = Eg0 − αT²/(T + β) (Y. P. Varshni, Physica 34, 149, 1967),
//   with Thurmond's parameters (C. D. Thurmond, J. Electrochem. Soc. 122, 1133, 1975, as
//   tabulated in Sze, 2nd ed.): Eg0 (eV), alpha (eV/K), beta (K). Only the shift from 300 K is
//   used, on the table's Eg, so T = 300 K is the table: Varshni's own Eg(300) differs from it by
//   4.5 meV (Si, 1.1245), 3.4 meV (Ge, 0.6634) and 1.5 meV (GaAs, 1.4225).
// - Nc and Nv ∝ T^1.5 (parabolic bands, effective masses held fixed).
// - mobilities ∝ T^−gn and T^−gp, the lattice-limited power laws for lightly doped material
//   (Sze, 2nd ed., ch. 1): approximate, and wrong where impurity scattering matters (heavy
//   doping, or cold). D follows from μ and T by Einstein's relation, D = μkT/q.
// - the electron affinity chi isn't moved.
const WITH_T = {
  Si: { Eg0: 1.17, alpha: 4.73e-4, beta: 636, gn: 2.42, gp: 2.2 },
  Ge: { Eg0: 0.7437, alpha: 4.774e-4, beta: 235, gn: 1.66, gp: 2.33 },
  GaAs: { Eg0: 1.519, alpha: 5.405e-4, beta: 204, gn: 1.0, gp: 2.1 },
};
const T_RANGE = [150, 500];

/**
 * A semiconductor from `SEMICONDUCTORS`, as a library piece: electrons `e-` and holes `h+`, with
 * c_ref = Nc and Nv, D from the mobilities by the Einstein relation at T, and μ°_e⁻ = 0 (φ
 * anchored at the conduction band) and μ°_h⁺ = Eg·F. The data are for 300 K; another T (150 to
 * 500 K) moves them: the gap by Varshni's shift from 300 K, Nc and Nv as T^1.5, the mobilities by
 * lattice-limited power laws (light doping). The material carries its `T`, and a device at
 * another temperature warns: give the device the same T.
 * Statistics are ideal (Boltzmann); add Fermi–Dirac to the material for degenerate doping.
 * @param {string} name 'Si', 'Ge' or 'GaAs'
 * @param {{ material?: string, T?: number }} [opts] the material's name, and the temperature
 *   (K, default 300, the table's)
 * @returns {{ species: object[], materials: Record<string, object> }}
 */
export function semiconductor(name, { material = name, T = 300 } = {}) {
  checkOptions('semiconductor', arguments[1], ['material', 'T']);
  const sc = SEMICONDUCTORS[name];
  if (!sc) fail(`semiconductor: no data for '${name}' (the library has ${Object.keys(SEMICONDUCTORS).join(', ')})`);
  if (typeof T !== 'number') fail(`semiconductor: T must be a number, the temperature in K, got ${shown(T)}`);
  if (!(Number.isFinite(T) && T >= T_RANGE[0] && T <= T_RANGE[1])) {
    fail(`semiconductor: T must be the temperature in K, from ${T_RANGE[0]} to ${T_RANGE[1]} (the range its temperature laws are good for), got ${T}`);
  }
  const { alpha, beta, gn, gp } = WITH_T[name];
  const t = T / 300, varshni = (T) => (alpha * T * T) / (T + beta);
  const Eg = sc.Eg - (varshni(T) - varshni(300));
  const VT = (GAS_CONSTANT * T) / FARADAY;
  const perCm3 = (n) => (n * 1e6) / AVOGADRO;
  return {
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: {
      [material]: {
        T, // the temperature the data are for: a device at another warns
        epsr: sc.epsr,
        species: {
          'e-': { D: sc.mun * powr(t, -gn) * 1e-4 * VT, mu0: 0, cRef: perCm3(sc.Nc * powr(t, 1.5)) },
          'h+': { D: sc.mup * powr(t, -gp) * 1e-4 * VT, mu0: Eg * FARADAY, cRef: perCm3(sc.Nv * powr(t, 1.5)) },
        },
      },
    },
  };
}

/**
 * A trap level in a semiconductor that holds charge: its two charge states as immobile species
 * (D = 0), filled and emptied by capturing carriers, e⁻ + X⁰ = X⁻ and X⁻ + h⁺ = X⁰ (acceptor-like)
 * or X⁺ + e⁻ = X⁰ and X⁰ + h⁺ = X⁺ (donor-like). In steady state it recombines at exactly the
 * `srh` law's rate with the same τ_n, τ_p and n₁, and its charge is in Poisson's equation, which
 * the `srh` law (the trap's occupancy eliminated) leaves out. Each node keeps its own total.
 * @param {{ species: object[], materials: Record<string, object> }} lib a library entry with e⁻
 *   and h⁺ (semiconductor()), returned with the trap's states added
 * @param {{ density: number, tauN: number, tauP: number, n1?: number, charge?: 'acceptor' | 'donor',
 *   name?: string, material?: string, T?: number }} opts density N_t (mol/m³); τ_n = 1/(c_n N_t)
 *   and τ_p = 1/(c_p N_t) (s) at that density; n₁, the electron density that half-fills it
 *   (mol/m³; default n_i, midgap); T (K, default the material's T, else 300)
 * @returns {{ library: object, bulkReactions: object[], c0: Record<string, number> }} the library
 *   for build(), the capture reactions, and c0 for each layer holding the trap (its total N_t)
 */
export function trap(lib, { density, tauN, tauP, n1, charge = 'acceptor', name = 'X', material, T } = {}) {
  checkOptions('trap', arguments[1], ['density', 'tauN', 'tauP', 'n1', 'charge', 'name', 'material', 'T']);
  if (!lib || typeof lib !== 'object' || typeof lib.materials !== 'object') fail('trap: give the library entry first, semiconductor(\'Si\') say, then { density, tauN, tauP }');
  const names = Object.keys(lib.materials);
  const m = material ?? (names.length === 1 ? names[0] : fail(`trap: the library has materials ${names.join(', ')}: say which, with material`));
  const mat = lib.materials[m];
  if (!mat) fail(`trap: no material '${m}' in the library (it has ${names.join(', ')})`);
  const el = mat.species?.['e-'], ho = mat.species?.['h+'];
  if (!el || !ho) fail(`trap: material '${m}' needs e- and h+ (a semiconductor)`);
  const pos = (v, what) => (Number.isFinite(v) && v > 0 ? v : fail(`trap: ${what}, got ${shown(v)}`));
  pos(density, 'density must be the trap density N_t > 0 (mol/m³; units.perCm3(1e15), say)');
  pos(tauN, 'tauN must be the electron lifetime τ_n = 1/(c_n N_t) > 0 (s)');
  pos(tauP, 'tauP must be the hole lifetime τ_p = 1/(c_p N_t) > 0 (s)');
  // (the material's own temperature, where it gives one: semiconductor()'s does)
  if (T === undefined) T = mat.T ?? 300;
  pos(T, 'T must be the temperature in K');
  if (mat.T !== undefined && Math.abs(T - mat.T) > 0.01) fail(`trap: T is ${T} K, but material '${m}' is for ${mat.T} K; leave T out (the trap takes the material's), or remake the material for ${T} K`);
  if (charge !== 'acceptor' && charge !== 'donor') fail(`trap: charge must be 'acceptor' (0/−) or 'donor' (+/0), got ${JSON.stringify(charge)}`);
  if (typeof name !== 'string' || !/^[A-Za-z][\w()]*$/.test(name)) fail(`trap: name must be a plain name for the trap (its states are name0 and name- or name+), got ${JSON.stringify(name)}`);
  const RT = GAS_CONSTANT * T, cRefE = el.cRef ?? lib.species.find((sp) => sp.name === 'e-')?.cRef;
  const cRefH = ho.cRef ?? lib.species.find((sp) => sp.name === 'h+')?.cRef;
  const ni = Math.sqrt(cRefE * cRefH * Math.exp(-(ho.mu0 - el.mu0) / RT));
  const half = n1 === undefined ? ni : pos(n1, 'n1 must be the electron density that half-fills the trap > 0 (mol/m³)');
  // The filled state's μ° sits an electron's μ at n₁ above the empty one's: half full at n = n₁.
  const [empty, full] = charge === 'acceptor' ? [`${name}0`, `${name}-`] : [`${name}+`, `${name}0`];
  const zEmpty = charge === 'acceptor' ? 0 : 1;
  const states = {
    [empty]: { D: 0, mu0: 0, cRef: density },
    [full]: { D: 0, mu0: el.mu0 + RT * Math.log(half / cRefE), cRef: density },
  };
  return {
    library: {
      species: [...lib.species, { name: empty, z: zEmpty }, { name: full, z: zEmpty - 1 }],
      materials: { ...lib.materials, [m]: { ...mat, species: { ...mat.species, ...states } } },
    },
    bulkReactions: [
      { equation: `e- + ${empty} = ${full}`, kf: { [m]: 1 / (tauN * density) } },
      { equation: `${full} + h+ = ${empty}`, kf: { [m]: 1 / (tauP * density) } },
    ],
    // (each node keeps the total, X⁰ + X⁻ say; the split is the solve's to find)
    c0: { [empty]: density / 2, [full]: density / 2 },
  };
}

/**
 * A metal from `METALS` as a conductor material, its carrier `e-`.
 * @param {string} name
 * @param {{ material?: string }} [opts]
 * @returns {{ species: object[], materials: Record<string, object> }}
 */
export function metal(name, { material = name } = {}) {
  checkOptions('metal', arguments[1], ['material']);
  const m = METALS[name];
  if (!m) fail(`metal: no data for '${name}' (the library has ${Object.keys(METALS).join(', ')})`);
  return { species: [{ name: 'e-', z: -1 }], materials: { [material]: { conductor: { species: 'e-', conductivity: m.conductivity } } } };
}
