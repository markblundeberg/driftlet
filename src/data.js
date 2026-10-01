// A small data library for driftlet/kit, every entry with its source (see docs/data.md). Kept
// small on purpose: a wrong library number is worse than a wrong demo number. Bulk properties
// only. How levels line up at an interface, and a surface's work function, belong to that
// interface or surface and are never in here.

import { DeviceError } from './device.js';
import { FARADAY, GAS_CONSTANT, AVOGADRO } from './constants.js';

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
export const WATER_EPSR = 78.36;

/**
 * Semiconductors at 300 K (Sze, Physics of Semiconductor Devices, 2nd ed., 1981, appendix):
 * relative permittivity, band gap Eg (eV), electron affinity chi (V), effective densities of
 * states Nc and Nv (cm⁻³), and drift mobilities of lightly doped material mun and mup
 * (cm²/(V·s)).
 */
export const SEMICONDUCTORS = freeze({
  Si: { epsr: 11.9, Eg: 1.12, chi: 4.05, Nc: 2.8e19, Nv: 1.04e19, mun: 1500, mup: 450 },
  Ge: { epsr: 16.0, Eg: 0.66, chi: 4.0, Nc: 1.04e19, Nv: 6.0e18, mun: 3900, mup: 1900 },
  GaAs: { epsr: 13.1, Eg: 1.424, chi: 4.07, Nc: 4.7e17, Nv: 7.0e18, mun: 8500, mup: 400 },
});

/** Electrical conductivity of metals at 20 °C, S/m (the inverse of the resistivities in the CRC Handbook). */
export const METALS = freeze({
  Ag: { conductivity: 6.29e7 },
  Cu: { conductivity: 5.96e7 },
  Au: { conductivity: 4.10e7 },
  Al: { conductivity: 3.77e7 },
  Zn: { conductivity: 1.69e7 },
  Li: { conductivity: 1.08e7 },
  Pt: { conductivity: 9.43e6 },
});

/**
 * Water with the named ions, as a library piece for `build()`: each ion as a species
 * (c_ref = 1 M) and the material `water` holding them, with D and μ° from `IONS`.
 * @param {string[]} names ions from `IONS`
 * @param {{ epsr?: number, material?: string }} [opts] permittivity (default that of water;
 *   0 makes the solution strictly neutral, without double layers), and the material's name
 * @returns {{ species: object[], materials: Record<string, object> }}
 */
export function aqueous(names, { epsr = WATER_EPSR, material = 'water' } = {}) {
  if (!Array.isArray(names) || names.length === 0) fail('aqueous: give a list of ions, e.g. aqueous([\'Na+\', \'Cl-\'])');
  const species = [], entries = {};
  for (const name of names) {
    const ion = IONS[name];
    if (!ion) fail(`aqueous: no data for '${name}' (the library has ${Object.keys(IONS).join(', ')})`);
    species.push({ name, z: ion.z, cRef: 1000 });
    entries[name] = { D: ion.D, mu0: ion.mu0 };
  }
  return { species, materials: { [material]: { epsr, species: entries } } };
}

/**
 * A semiconductor from `SEMICONDUCTORS`, as a library piece: electrons `e-` and holes `h+`, with
 * c_ref = Nc and Nv, D from the mobilities by the Einstein relation at 300 K, and μ°_e⁻ = 0 (φ
 * anchored at the conduction band) and μ°_h⁺ = Eg·F. The data are for 300 K, so use `T: 300`.
 * Statistics are ideal (Boltzmann); add Fermi–Dirac to the material for degenerate doping.
 * @param {string} name 'Si', 'Ge' or 'GaAs'
 * @param {{ material?: string }} [opts]
 * @returns {{ species: object[], materials: Record<string, object> }}
 */
export function semiconductor(name, { material = name } = {}) {
  const sc = SEMICONDUCTORS[name];
  if (!sc) fail(`semiconductor: no data for '${name}' (the library has ${Object.keys(SEMICONDUCTORS).join(', ')})`);
  const VT = (GAS_CONSTANT * 300) / FARADAY;
  const perCm3 = (n) => (n * 1e6) / AVOGADRO;
  return {
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: {
      [material]: {
        epsr: sc.epsr,
        species: {
          'e-': { D: sc.mun * 1e-4 * VT, mu0: 0, cRef: perCm3(sc.Nc) },
          'h+': { D: sc.mup * 1e-4 * VT, mu0: sc.Eg * FARADAY, cRef: perCm3(sc.Nv) },
        },
      },
    },
  };
}

/**
 * A metal from `METALS` as a conductor material, its carrier `e-`.
 * @param {string} name
 * @param {{ material?: string }} [opts]
 * @returns {{ species: object[], materials: Record<string, object> }}
 */
export function metal(name, { material = name } = {}) {
  const m = METALS[name];
  if (!m) fail(`metal: no data for '${name}' (the library has ${Object.keys(METALS).join(', ')})`);
  return { species: [{ name: 'e-', z: -1 }], materials: { [material]: { conductor: { species: 'e-', conductivity: m.conductivity } } } };
}
