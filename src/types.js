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
 *   { type: 'ideal', species: string[] } |
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
 * A material: either ordinary (epsr, species, optional statistics) or a conductor.
 * @typedef {object} MaterialDef
 * @property {number} [epsr] relative permittivity; 0 makes the material strictly neutral
 * @property {Record<string, MaterialSpeciesDef>} [species] species present in this material
 * @property {StatisticsDef[]} [statistics] non-ideal statistics (default: all ideal)
 * @property {{ species: string, conductivity: number }} [conductor] a metal or fast ion
 *   conductor: its one carrier and conductivity (S/m); its only level is the carrier's (a metal's
 *   Fermi level), and φ is undefined inside
 * @property {number} [T] the temperature its data are for (K): a device at another temperature
 *   warns (semiconductor() sets it)
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
 * A table against the device's x (m): piecewise linear between its points, constant beyond them.
 * @typedef {{ x: number[], values: number[] }} Profile
 */

/**
 * @typedef {object} RegionDef
 * @property {string} material material name
 * @property {number} length m
 * @property {string} [name]
 * @property {number} [fixedCharge] immobile charge density, C/m³
 * @property {Record<string, number | Profile>} [c0] initial concentrations,
 *   mol/m³ (required for species not connected to a contact: it fixes the amount they conserve);
 *   or a profile against the device's x (m), piecewise linear; either is the starting state
 * @property {number} [velocity] imposed flow toward +x, m/s, carrying every mobile species
 * @property {number} [mixing] eddy (turbulent) mixing diffusivity, m²/s; current-free
 * @property {GridOptions} [grid]
 */

/** @typedef {'pinned' | 'neutral' | { type: 'capacitive', C: number }} InterfacePhiLaw */
/**
 * A species' law across a face. conductance: J = G (V_L − V_R), G in S/m², charged species;
 * permeability: Goldman–Hodgkin–Katz electrodiffusion, P in m/s, ideal statistics on both sides,
 * not beside a conductor. gates: the face's gates that scale G or P, with their exponents
 * ({ m: 3, h: 1 } for m³h; not beside a conductor).
 * @typedef {'equilibrium' | 'blocked' | { type: 'equilibrium' | 'blocked' }
 *   | { type: 'conductance', G: number, gates?: Record<string, number> }
 *   | { type: 'permeability', P: number, gates?: Record<string, number> }} InterfaceSpeciesLink
 */

/**
 * One of a gate's rates (1/s) in the voltage V it follows, in NeuroML's standard forms (midpoint
 * and scale in V, scale ≠ 0): 'exp', rate·e^((V − midpoint)/scale); 'sigmoid',
 * rate/(1 + e^((midpoint − V)/scale)); 'expLinear', rate·y/(1 − e^(−y)), y = (V − midpoint)/scale.
 * @typedef {{ type: 'exp' | 'sigmoid' | 'expLinear', rate: number, midpoint: number, scale: number }} GateRate
 */
/**
 * A Hodgkin–Huxley gate: a fraction open x in [0, 1], dx/dt = α(V)(1 − x) − β(V)x.
 * @typedef {{ alpha: GateRate, beta: GateRate }} GateDef
 */

/**
 * A reaction at a face: participants on each side with signed stoichiometric coefficients
 * (ν < 0 consumed by the forward reaction). Butler–Volmer kinetics (k0, alpha), or instead SRH
 * (srh) or saturating (vmax and K), neither with a conductor's carrier.
 * @typedef {object} FaceReactionDef
 * @property {string} [equation] the reaction as an equation, e.g. 'Ag+ + e- = Ag(s)', each
 *   participant on the side that holds it ('Li+(left)' where both do); or give left and right
 * @property {Record<string, number>} [left] participants on the left side (species present there, or fixed neutrals)
 * @property {Record<string, number>} [right] participants on the right side
 * @property {Record<string, number>} [fixed] μ (J/mol) of fixed-activity participants that aren't species
 * @property {number} [k0] Butler–Volmer: standard rate constant, mol/(m²·s) (required unless srh or vmax)
 * @property {number} [alpha] with k0: transfer coefficient, 0…1 (default 0.5)
 * @property {{ vn: number, vp: number, n1?: number }} [srh] instead: trap-assisted recombination, for a
 *   reaction consuming one negative and one positive species and making none; capture velocities
 *   vn, vp in m/s, n1 in mol/m³ (default n_i, a midgap trap)
 * @property {number} [vmax] instead: saturating kinetics, the largest rate, mol/(m²·s)
 * @property {Record<string, number>} [K] with vmax: a half-saturation concentration (mol/m³) for
 *   every species the forward reaction consumes
 */
/**
 * A face's diffuse layer on one side, in a solution.
 * @typedef {object} DiffuseLayerResult
 * @property {number} psi the drop from the side's edge node to the plane, φ_plane − φ_edge (V)
 * @property {number} charge the layer's charge (C/m²)
 * @property {Record<string, number>} excess each ion's surface excess, ∫ (c − c_edge) dx (mol/m²)
 * @property {Record<string, number>} plane each ion's concentration at the plane, c e^{−zFψ/RT} (mol/m³; on a
 *   lattice, over its crowding; for what a face reaction takes or makes there, times its level's drop across the layer)
 */
/**
 * One face between consecutive regions. For different materials under a 'pinned' or
 * 'capacitive' law, exactly one alignment (step or dipole) is required.
 * @typedef {object} InterfaceDef
 * @property {{ species: string, value: number }} [step] standard-level step right − left, J/mol
 * @property {number} [dipole] φ jump right − left, V
 * @property {number} [zeroCharge] beside a conductor region: V_F − φ_edge at zero charge, V
 * @property {InterfacePhiLaw} [phi] electrostatic law (default 'pinned'; 'neutral' between ε = 0 materials;
 *   beside a conductor 'neutral' or capacitive, required next to an ε > 0 material)
 * @property {number} [sheetCharge] C/m²
 * @property {{ epsr: number }} [diffuse] a Gouy–Chapman diffuse layer below the grid on each strictly
 *   neutral (ε = 0) side of a charged face, of relative permittivity epsr: the side's share of the
 *   face's charge held as its ions' surface excesses, the face's reactions at its plane
 * @property {Record<string, InterfaceSpeciesLink>} [species] per-species laws (default equilibrium
 *   where present on both sides and mobile on one, blocked otherwise; a species reacting at the
 *   face and present on both sides has no default)
 * @property {FaceReactionDef[]} [reactions] reactions at the face
 * @property {Record<string, GateDef>} [gates] voltage-gated channels: gates following
 *   V = φ_right − φ_left, which its links name (not beside a conductor)
 */

/**
 * @typedef {object} BulkReactionDef
 * @property {Record<string, number>} [nu] signed stoichiometric coefficients (ν < 0 consumed)
 * @property {string} [equation] the reaction as an equation, e.g. 'e- + h+ = 0' (instead of nu)
 * @property {Record<string, number>} [fixed] μ (J/mol) of fixed-activity participants that aren't species
 * @property {Record<string, number | Profile>} [kf] mass action: forward rate constant per material
 *   name, or a profile of it against the device's x (an absorption profile, for photogeneration);
 *   give kf or srh
 * @property {Record<string, { tauN: number, tauP: number, n1?: number }>} [srh] instead: SRH
 *   (trap-assisted) recombination per material name, for a reaction consuming one negative and
 *   one positive species and making none; lifetimes in s, n1 in mol/m³ (default n_i, a midgap trap)
 */

/**
 * Contacts use the same laws as internal faces; the outside is a phase with known levels
 * V_i = V + offset_i (a rigid ladder that the circuit slides by V).
 * @typedef {'blocked' | 'equilibrium' | { type: 'blocked' }
 *   | { type: 'equilibrium', offset?: number, mu?: number }
 *   | { type: 'conductance', G: number, offset?: number }
 *   | { type: 'exchange', k: number, mu: number }
 *   | { type: 'velocity', v: number, offset?: number, mu?: number }} ContactSpeciesLink
 */
/**
 * 'bulk': the end node is plain bulk (locally neutral, no double layer); 'neutral': no charge at
 * the face (D = 0); capacitive: Stern layer or gate; 'pinned': φ_edge = V − zeroCharge. The last
 * two take the contact's zeroCharge.
 * @typedef {'bulk' | 'neutral' | 'pinned' | { type: 'capacitive', C: number }} ContactPhiLink
 */

/**
 * @typedef {object} ContactDef
 * @property {Source} [V] held terminal voltage, V (default 0, if no I)
 * @property {Source} [I] or: the current driven into the device, A/m² (the voltage floats)
 * @property {number} [R] with V: a series resistance, Ω·m² (a source behind a resistor)
 * @property {string} [terminal] species whose voltage V is
 * @property {Record<string, ContactSpeciesLink>} [species] per-species links (default blocked)
 * @property {{ c: Record<string, number>, reference?: string, offset?: number, material?: string }} [bath]
 *   instead of species: a neutral composition c (mol/m³) held at the contact, anchored through a
 *   charged reference species (offset: its V_i − V, V); with no reference, V is the bath's φ. It's
 *   the end material at that composition, or the material it names (a solution beside a charged region)
 * @property {ContactPhiLink} [phi] required when any species connects (default 'neutral'; 'bulk' with a bath)
 * @property {number} [zeroCharge] for a capacitive or pinned law: V − φ_edge at zero charge, V
 */

/**
 * A terminal's source: a constant, or a piecewise-linear waveform through (t, value), constant
 * beyond the points, or periodic with period t_last − t_0 if repeat is true.
 * @typedef {number | { t: number[], values: number[], repeat?: boolean }} Source
 */

/**
 * An internal port: an outside phase exchanging with a window of one region's nodes. On a
 * conductor region it attaches to the whole conductor (no window), and a conductance G is per area
 * (S/m²) instead of per volume (S/m³).
 * @typedef {object} PortDef
 * @property {string} [name]
 * @property {string | number} region region name or index
 * @property {number} [from] window start, m from the region's left end (default 0)
 * @property {number} [to] window end, m (default: the region's length)
 * @property {Source} [V] held terminal voltage, V (default 0, if no I)
 * @property {Source} [I] or: the current driven into the device, A/m² (e.g. 0 for a reference electrode)
 * @property {number} [R] with V: a series resistance, Ω·m²
 * @property {string} [terminal] species whose offset defaults to 0; with reactions, the electrode's carrier
 * @property {Record<string, 'blocked' | 'equilibrium' | { type: 'blocked' } | { type: 'equilibrium', offset?: number, mu?: number }
 *   | { type: 'conductance', G: number, offset?: number, gates?: Record<string, number> }
 *   | { type: 'exchange', k: number, mu: number }>} [species] per-species links, per volume (default
 *   blocked): G in S/m³ (S/m² on a conductor), k in mol/(m³·s); a conductance's gates name the
 *   port's gates that scale G, with their exponents
 * @property {{ c: Record<string, number> }} [bath] the outside's composition (mol/m³): V is then its φ,
 *   and each linked species' level there follows from it (unless its link gives an offset or mu)
 * @property {{ equation: string, fixed?: Record<string, number>, k0: number, alpha?: number, bare?: boolean }[]} [reactions] an
 *   electrode spread through the window: Butler–Volmer per area (k0 in mol/(m²·s), alpha default 0.5)
 *   between the region's species and the terminal species at the port's level
 * @property {Record<string, { mu0: number, capacity: number, theta0?: number }>} [surface] with reactions: species
 *   on the electrode's sites, Langmuir (μ = μ° + RT ln(θ/θ₀)), capacity in mol of sites per m² of electrode,
 *   shared, theta0 the starting coverage (default 1e-6); reactions run on Langmuir sites (what fills
 *   sites needs bare ones), and one with `bare: true` needs a bare site even if it fills none (θ₀ = 1 − Σθ)
 * @property {{ C: number, zeroCharge?: number }} [capacitance] a capacitance spread through the window (F/m² of
 *   electrode): σ = C (V − zeroCharge − φ) on the port's side, aσ in each node's charge balance; a gate
 *   along a channel, a double layer
 * @property {number | { x: number[], values: number[] }} [area] with reactions or a capacitance: electrode area per volume,
 *   m²/m³ (1/h for a film of thickness h), or a profile against the device's x
 * @property {Record<string, GateDef>} [gates] with a capacitance: voltage-gated channels at each node of
 *   the window, following the voltage across it, φ − (V − zeroCharge), which its conductance links name
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
 * @property {PortDef[]} [ports] internal ports
 * @property {'planar' | 'spherical' | 'cylindrical' | { type: 'planar' | 'spherical' | 'cylindrical', r0?: number }
 *   | { area: { x: number[], values: number[] } }} [geometry] the cross-section A(x): planar (A = 1 m², the
 *   default: currents read per m²), shells about a centre at r = r0 (m) left of x = 0 (cylindrical: per
 *   metre of length), or A in m² against x. Currents, fluxes, amounts and charges are then totals through A.
 * @property {GridOptions} [grid]
 */

/**
 * @typedef {object} ContactResult
 * @property {number} V terminal voltage, V
 * @property {Record<string, number>} flux particle flux toward +x through this contact, mol/(m²·s)
 * @property {number} D displacement toward +x at this contact, C/m²
 * @property {number} current total current toward +x, A/m²
 * @property {Record<string, string>} links how each species that crosses it does (its link type)
 */

/**
 * The device's totals through its cross-section (currents, fluxes through contacts and ports,
 * charges, amounts, reaction totals, impedance) are per m² of a planar device, as below, and
 * whole (A, mol/s, C, mol, Ω; per metre of a cylinder) for one with a `geometry`. Rates at a
 * face or an electrode are always per m² of it.
 * @typedef {object} Solution
 * @property {Float64Array} x node positions, m (interface positions appear twice)
 * @property {Int32Array} region region index per node
 * @property {{ name: string, z: number }[]} species the species, in order, with their charge numbers
 * @property {{ name: string, material: string, x0: number, x1: number }[]} regions each region's name, material and extent, m
 * @property {Float64Array} phi bookkeeping φ, V (NaN where undefined)
 * @property {Record<string, Float64Array>} c concentrations, mol/m³ (NaN where absent)
 * @property {Record<string, Float64Array>} mu electrochemical potentials μ̄, J/mol
 * @property {Record<string, Float64Array>} muStd standard levels μ° + zFφ, J/mol
 * @property {Record<string, Float64Array>} V species voltages μ̄/(zF), V (charged species)
 * @property {Record<string, Float64Array>} Vstd standard levels as voltages, V
 * @property {Record<string, Float64Array>} flux each species' flux toward +x across each segment, between
 *   nodes g and g + 1 (a face's across its pair of nodes), mol/(m²·s): the fluxes the balances add up
 * @property {Record<string, Float64Array>} D diffusivity at each node, m²/s (NaN where absent, or in a conductor)
 * @property {number} T temperature, K
 * @property {'planar' | 'cylindrical' | 'spherical' | 'profile'} geometry the device's (amounts and rates per area only if planar)
 * @property {number} current current toward +x through the device, A/m²
 * @property {number} terminalVoltage V_right − V_left, V
 * @property {Record<string, { V: number, current: number }>} terminals each terminal (left, right,
 *   then the ports by name): its voltage, and its current into the device (A/m²)
 * @property {{ left: ContactResult, right: ContactResult }} contacts
 * @property {{ left?: { V: number, D: number, charge: number }, right?: { V: number, D: number, charge: number } }} plates the charge on a gate or
 *   Stern plate at a capacitive contact
 * @property {{ name: string, V: number, flux: Record<string, number>, current: number, x?: Float64Array,
 *   area?: Float64Array, rates?: Float64Array[], coverage?: Record<string, Float64Array>, bare?: Float64Array,
 *   sigma?: Float64Array, charge?: number, gates?: Record<string, Float64Array>, Vm?: Float64Array }[]} ports
 *   what each internal port brings into the device (current: the charging current included); with
 *   reactions or a capacitance, also its window's nodes (x, m) and area per volume there (m²/m³); with
 *   reactions, each reaction's net rate per area at each node (mol/(m²·s), forward); with a surface, each
 *   species' coverage and the bare fraction; with a capacitance, σ per area (C/m²) and the charge its
 *   side holds; with gates, each gate's fraction open at each node and the voltage they follow (V)
 * @property {{ left: string, right: string, dipole: number, sheetCharge: number, D: number, N: Record<string, number>,
 *   rates: number[], links: Record<string, string>, oneWay?: Record<string, [number, number]>,
 *   gates?: Record<string, number>, V?: number,
 *   diffuse?: { left?: DiffuseLayerResult, right?: DiffuseLayerResult } }[]} interfaces each face: the regions it joins, its dipole
 *   (V) and sheet charge (C/m²), its displacement D (C/m²), each species' flux N across it and each face
 *   reaction's net rate (mol/(m²·s), toward +x and forward), how each species that crosses it does (its
 *   link type); for those crossing by permeability, the one-way fluxes [toward +x, toward −x] whose
 *   difference is N; with gates, each one's fraction open and the voltage they follow, φ_right − φ_left (V; between the diffuse layers' planes, where it has them);
 *   with a diffuse layer, each side's
 * @property {{ rate: Float64Array, forward: Float64Array, nu: Record<string, number>, regions: number[], total: number }[]} bulkReactions
 *   each bulk reaction's net rate per node (mol/(m³·s), forward minus backward, NaN where it doesn't
 *   run), its one-way forward rate (the backward is forward − rate), what one forward reaction makes
 *   of each species (ν, negative for what it consumes), and the net rate per region and in all (mol/(m²·s))
 * @property {number} charge total charge in the device, C/m²
 * @property {{ species: string, regions: number[], spectator: boolean, connected: boolean, reactive: boolean,
 *   amount: number, reference: number, intake: number, drift: number }[]} conservation
 * @property {string[]} warnings
 * @property {boolean} converged
 * @property {number} iterations Newton iterations
 * @property {number} [steps] solve(): the steps it took to the steady state; advance(): accepted steps
 * @property {number} [substeps] step(): how many pieces the interval was split into
 * @property {(number | object)[]} [history] solver bookkeeping: solve(), a record of its steps and
 *   stages; step(), Newton's largest update per iteration (thermal units)
 * @property {number} time s
 * @property {boolean} steady whether this is a converged steady state (solve()), not a transient's
 * @property {boolean} [done] advance(): whether tEnd was reached (false if the frame budget ran out)
 * @property {number} [rejected] advance(): steps rejected by error control or Newton failure
 * @property {{ t: number[], current: number[], voltage: number[], probes?: number[][] }} [trace] advance():
 *   the right contact's current (A/m², toward +x) and V_right − V_left after each accepted step, and
 *   each probe's reading: inside, interpolated between the nodes around its x (mol/m³ for quantity
 *   'c', V for 'V' or 'phi'); at a face, a species' flux through it (mol/(m²·s), toward +x) or a
 *   gate's fraction open
 */

/**
 * @typedef {object} ImpedanceResult
 * @property {Float64Array} f frequencies, Hz
 * @property {{ re: Float64Array, im: Float64Array }} Z impedance δV/δI at the terminal, I into the device, Ω·m²
 * @property {Float64Array} relError each Z's estimated error, relative to |Z| (on Re Z and Im Z alike): the current's round-off
 *   and the solves' residuals; cautious (it can exceed the true error many times). Over 1%, or over 10·|Re Z|/|Z|,
 *   a warning says so
 * @property {{ phi: ComplexProfile, mu: Record<string, ComplexProfile>, c: Record<string, ComplexProfile> }[]} [profiles]
 *   per frequency, per unit excitation (volt, or A/m² in current mode)
 * @property {string[]} warnings the steady state's warnings, which hold for its impedance too
 */

/** @typedef {{ re: Float64Array, im: Float64Array }} ComplexProfile */

export {};
