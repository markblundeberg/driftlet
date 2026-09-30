// Physical constants (SI, CODATA exact values where defined) and unit helpers.

/** Elementary charge, C. */
export const ELEMENTARY_CHARGE = 1.602176634e-19;
/** Avogadro constant, 1/mol. */
export const AVOGADRO = 6.02214076e23;
/** Boltzmann constant, J/K. */
export const BOLTZMANN = 1.380649e-23;
/** Faraday constant, C/mol. */
export const FARADAY = ELEMENTARY_CHARGE * AVOGADRO;
/** Molar gas constant, J/(mol·K). */
export const GAS_CONSTANT = BOLTZMANN * AVOGADRO;
/** Vacuum permittivity, F/m (CODATA 2022). */
export const EPS0 = 8.8541878188e-12;

/** Unit conversions into driftlet's SI units. */
export const units = Object.freeze({
  /** mol/L → mol/m³ */
  molar: (x) => x * 1e3,
  /** mmol/L → mol/m³ */
  millimolar: (x) => x,
  /** particles per cm³ → mol/m³ */
  perCm3: (x) => (x * 1e6) / AVOGADRO,
  /** particles per m³ → mol/m³ */
  perM3: (x) => x / AVOGADRO,
  /** eV per particle → J/mol */
  eV: (x) => x * FARADAY,
  /** J/mol → eV per particle */
  toEV: (x) => x / FARADAY,
  /** nm → m */
  nm: (x) => x * 1e-9,
  /** µm → m */
  um: (x) => x * 1e-6,
  /** cm²/s → m²/s */
  cm2PerS: (x) => x * 1e-4,
});
