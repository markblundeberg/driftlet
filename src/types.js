// JSDoc types for device definitions and solutions (for editors; no runtime content).
// Reference: docs/device.md. Units are SI throughout.

/**
 * @typedef {object} SpeciesDef
 * @property {string} name
 * @property {number} z charge number (integer; 0 for neutral species)
 * @property {number} [cRef] default reference concentration, mol/m³
 */

/**
 * @typedef {object} MaterialSpeciesDef
 * @property {number} D coefficient in N = −(D c/RT)∇μ̄, m²/s
 * @property {number} mu0 standard chemical potential, J/mol
 * @property {number} [cRef] concentration that mu0 refers to, mol/m³
 */

/**
 * Statistics model: how the listed species' concentrations follow from their reduced
 * potentials ζ = (μ̄ − zFφ − μ°)/RT. Species not in any model are ideal, c = c_ref e^ζ.
 * @typedef {(
 *   { type: 'fermi-dirac', species: string[], order?: 0.5 | 0 } |
 *   { type: 'lattice', species: string[], cMax: number } |
 *   { type: 'redlich-kister', species: [string], cMax: number, A?: number[] } |
 *   { type: 'debye-huckel', species: string[], epsr?: number, a?: number } |
 *   { type: 'insertion', species: [string, string], cMax: number, A?: number[],
 *     ocv?: { x: number[], E: number[], muRef: number } } |
 *   { type: 'custom', species: string[],
 *     evaluate: (zeta: number[]) => { c: number[], dcdzeta: number[][] } }
 * )} StatisticsDef
 */

/**
 * @typedef {object} MaterialDef
 * @property {number} epsr relative permittivity; 0 makes the material strictly neutral
 * @property {Record<string, MaterialSpeciesDef>} species species present in this material
 * @property {StatisticsDef[]} [statistics] non-ideal statistics (default: all ideal)
 */

/**
 * @typedef {object} GridOptions
 * @property {number} [hmin] cell size at region ends, m
 * @property {number} [hmax] largest cell size, m
 * @property {number} [ratio] geometric growth ratio (default 1.2)
 * @property {number} [minCells] minimum cells per region (default 8)
 * @property {number} [hLeft] left-end cell size, overriding hmin
 * @property {number} [hRight] right-end cell size, overriding hmin
 */

/**
 * @typedef {object} RegionDef
 * @property {string} material material name
 * @property {number} length m
 * @property {string} [name]
 * @property {number} [fixedCharge] immobile charge density, C/m³
 * @property {Record<string, number>} [c0] initial concentrations, mol/m³ (required for species
 *   not connected to a contact: it fixes the amount they conserve)
 * @property {number} [velocity] imposed flow toward +x, m/s, carrying every mobile species
 * @property {number} [mixing] eddy (turbulent) mixing diffusivity, m²/s; current-free
 * @property {GridOptions} [grid]
 */

/** @typedef {'dipole' | 'neutral' | { type: 'capacitive', C: number }} InterfacePhiLaw */
/** @typedef {'equilibrium' | 'blocked' | { type: 'conductance', G: number }} InterfaceSpeciesLink */

/**
 * @typedef {object} TransferReactionDef
 * @property {Record<string, number>} transfer species crossing left → right (forward), with coefficients
 * @property {number} k0 standard rate constant, mol/(m²·s)
 * @property {number} [alpha] transfer coefficient, 0…1 (default 0.5)
 */

/**
 * One face between consecutive regions. For different materials under a 'dipole' or
 * 'capacitive' law, exactly one alignment (step or dipole) is required.
 * @typedef {object} InterfaceDef
 * @property {{ species: string, value: number }} [step] standard-level step right − left, J/mol
 * @property {number} [dipole] φ jump right − left, V
 * @property {InterfacePhiLaw} [phi] electrostatic law (default 'dipole'; 'neutral' between ε = 0 materials)
 * @property {number} [sheetCharge] C/m²
 * @property {Record<string, InterfaceSpeciesLink>} [species] per-species laws (default equilibrium)
 * @property {TransferReactionDef[]} [reactions] Butler–Volmer transfer across the face
 */

/**
 * @typedef {object} BulkReactionDef
 * @property {Record<string, number>} [reactants] participants and stoichiometric coefficients
 * @property {Record<string, number>} [products]
 * @property {Record<string, number>} [fixed] μ (J/mol) of fixed-activity participants that aren't species
 * @property {Record<string, number>} kf forward rate constant per material name
 */

/**
 * Contacts use the same laws as internal faces; the outside is a phase with known levels
 * V_i = V + offset_i (a rigid ladder that the circuit slides by V).
 * @typedef {'blocked' | 'equilibrium'
 *   | { type: 'equilibrium', offset?: number, mu?: number }
 *   | { type: 'conductance', G: number, offset?: number }} ContactSpeciesLink
 */
/**
 * 'bulk': the end node is plain bulk (locally neutral, no double layer); 'neutral': no charge at
 * the face (D = 0); capacitive: Stern layer or gate; dipole: pinned φ_edge = V − zeroCharge.
 * @typedef {'bulk' | 'neutral'
 *   | { type: 'capacitive', C: number, zeroCharge: number }
 *   | { type: 'dipole', zeroCharge: number }} ContactPhiLink
 */

/**
 * Electrode reaction, written as reduction when electrons > 0: Σν_R R + n e⁻(metal) ⇌ Σν_P P.
 * @typedef {object} ElectrodeReactionDef
 * @property {Record<string, number>} [reactants]
 * @property {Record<string, number>} [products]
 * @property {number} [electrons] electrons taken from the metal (integer)
 * @property {Record<string, number>} [fixed] μ (J/mol) of fixed-activity participants
 * @property {number} k0 standard rate constant, mol/(m²·s)
 * @property {number} [alpha] transfer coefficient, 0…1 (default 0.5)
 */

/**
 * @typedef {object} ContactDef
 * @property {number} [V] terminal voltage, V (default 0)
 * @property {string} [terminal] species whose voltage V is
 * @property {Record<string, ContactSpeciesLink>} [species] per-species links (default blocked)
 * @property {{ c: Record<string, number>, reference: string, offset?: number }} [bath] a neutral
 *   composition held at the contact, anchored through a charged reference species
 * @property {ContactPhiLink} [phi] required when any species or reaction connects (default 'neutral')
 * @property {ElectrodeReactionDef[]} [reactions]
 */

/**
 * @typedef {{ mode: 'voltage' }
 *   | { mode: 'current', I: number }
 *   | { mode: 'load', R: number, V?: number }} CircuitDef
 */

/**
 * A device definition: plain, serialisable data.
 * @typedef {object} DeviceDefinition
 * @property {number} [T] temperature, K (default 298.15)
 * @property {SpeciesDef[]} species
 * @property {Record<string, MaterialDef>} materials
 * @property {RegionDef[]} regions left to right
 * @property {(InterfaceDef | null | undefined)[]} [interfaces] one per face, interfaces[f] between regions f and f+1
 * @property {BulkReactionDef[]} [bulkReactions]
 * @property {{ left?: ContactDef, right?: ContactDef }} [contacts]
 * @property {CircuitDef} [circuit] default: each contact at its own V
 * @property {GridOptions} [grid]
 */

/**
 * @typedef {object} ContactResult
 * @property {number} V terminal voltage, V
 * @property {Record<string, number>} flux particle flux toward +x through this contact, mol/(m²·s)
 * @property {number} D displacement toward +x at this contact, C/m²
 * @property {number} current total current toward +x, A/m²
 */

/**
 * @typedef {object} Solution
 * @property {Float64Array} x node positions, m (interface positions appear twice)
 * @property {Int32Array} region region index per node
 * @property {Float64Array} phi bookkeeping φ, V (NaN where undefined)
 * @property {Record<string, Float64Array>} c concentrations, mol/m³ (NaN where absent)
 * @property {Record<string, Float64Array>} mu electrochemical potentials μ̄, J/mol
 * @property {Record<string, Float64Array>} muStd standard levels μ° + zFφ, J/mol
 * @property {Record<string, Float64Array>} V species voltages μ̄/(zF), V (charged species)
 * @property {Record<string, Float64Array>} Vstd standard levels as voltages, V
 * @property {number} current terminal current toward +x, A/m²
 * @property {number} terminalVoltage V_right − V_left, V
 * @property {{ left: ContactResult, right: ContactResult }} contacts
 * @property {{ left?: { V: number, D: number, charge: number }, right?: { V: number, D: number, charge: number } }} gates
 * @property {{ dipole: number, sheetCharge: number, D: number, N: Record<string, number> }[]} interfaces
 * @property {number} charge total charge in the device, C/m²
 * @property {{ species: string, regions: number[], spectator: boolean, connected: boolean, reactive: boolean,
 *   amount: number, reference: number, intake: number, drift: number }[]} conservation
 * @property {string[]} warnings
 * @property {boolean} converged
 * @property {number} iterations
 * @property {number} time s
 * @property {boolean} [done] advance(): whether tEnd was reached (false if the frame budget ran out)
 * @property {number} [rejected] advance(): steps rejected by error control or Newton failure
 * @property {{ t: number[], current: number[], voltage: number[] }} [trace] advance(): terminal
 *   current (A/m², toward +x) and voltage after each accepted step
 */

/**
 * @typedef {object} ImpedanceResult
 * @property {Float64Array} f frequencies, Hz
 * @property {{ re: Float64Array, im: Float64Array }} Z impedance −δV/δI, Ω·m²
 * @property {{ phi: ComplexProfile, mu: Record<string, ComplexProfile>, c: Record<string, ComplexProfile> }[]} [profiles]
 *   per frequency, per unit excitation (volt, or A/m² in current mode)
 */

/** @typedef {{ re: Float64Array, im: Float64Array }} ComplexProfile */

export {};
