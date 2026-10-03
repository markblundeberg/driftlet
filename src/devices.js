// Curated devices, for driftlet/kit: whole definitions from a published parameter set (light and
// contacts included), with reference results from an independent code to check a run against.
// Each returns ordinary definition data, to print, edit and pass to new Device.

import { AVOGADRO, FARADAY, GAS_CONSTANT, units } from './constants.js';
import { photogeneration } from './sources.js';
import { drive } from './stack.js';

const m3 = (x) => x / AVOGADRO, eV = units.eV;

// IonMonger's default parameters (parameters_template.m, IonMonger 2.0): Courtier, Cave, Walker,
// Richardson & Foster, J. Comput. Electron. 18, 1435 (2019). Levels in eV on the vacuum scale,
// densities per m³, D in m²/s, lifetimes in s, recombination velocities in m/s.
const P = {
  T: 298, b: 400e-9, eps: 24.1, alpha: 1.3e7, Ec: -3.7, Ev: -5.4, Dn: 1.7e-4, Dp: 1.7e-4, gc: 8.1e24, gv: 5.8e24,
  N0: 1.6e25, DIinf: 6.5e-8, EAI: 0.58, // vacancy D = DIinf·exp(−EAI/kT), ≈ 1.0e-17 m²/s at 298 K
  dE: 1e24, gcE: 5e25, EcE: -4.0, bE: 100e-9, epsE: 10, DE: 1e-5,
  dH: 1e24, gvH: 5e25, EvH: -5.1, bH: 200e-9, epsH: 3, DH: 1e-6,
  tn: 3e-9, tp: 3e-7, vnE: 1e5, vpE: 10, vnH: 0.1, vpH: 1e5, Fph: 1.4e21,
};

/**
 * IonMonger's default planar perovskite solar cell: TiO₂ (100 nm) | MAPbI₃ (400 nm) with mobile
 * iodide vacancies over an equal immobile background | spiro-OMeTAD (200 nm), at 298 K, lit
 * through the TiO₂. Species e⁻, h⁺ and the vacancies 'V+'. Recombination is SRH at both faces
 * (IonMonger's interface law is driftlet's face SRH, with n₁ from the TiO₂'s electron density)
 * and, with `bulkSRH`, SRH through midgap traps in the perovskite (τₙ = 3 ns, τₚ = 300 ns).
 * Light is 1.4e21 photons/(m²·s) at α = 1.3e7 /m (IonMonger's one sun), times `light`. The left
 * contact is the TiO₂'s, at 0 V; the right is the spiro's, at `V`, so photocurrent is positive.
 *
 * `PEROVSKITE_SCANS` holds IonMonger's own results for this cell, scan by scan.
 * @param {{ V?: import('./types.js').Source | { V?: import('./types.js').Source, I?: import('./types.js').Source, R?: number }, light?: number, bulkSRH?: boolean }} [opts]
 * @returns {import('./types.js').DeviceDefinition}
 */
export function perovskiteCell({ V = 1.2, light = 1, bulkSRH = true } = {}) {
  const { T } = P, VT = (GAS_CONSTANT * T) / FARADAY;
  const ni = Math.sqrt(P.gc * P.gv) * Math.exp(-(P.Ec - P.Ev) / (2 * VT));
  const n0 = P.gc * Math.exp((P.EcE + VT * Math.log(P.dE / P.gcE) - P.Ec) / VT); // the perovskite's n at the TiO₂'s Fermi level
  const DI = P.DIinf * Math.exp(-P.EAI / VT);
  return {
    T,
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
      { name: 'V+', z: 1 }, // iodide vacancies
    ],
    materials: {
      // Levels on the vacuum scale: μ° is the band edge's depth (minus it, for electrons).
      'TiO₂': { epsr: P.epsE, species: { 'e-': { D: P.DE, mu0: eV(P.EcE), cRef: m3(P.gcE) } } },
      'MAPbI₃': {
        epsr: P.eps,
        species: {
          'e-': { D: P.Dn, mu0: eV(P.Ec), cRef: m3(P.gc) },
          'h+': { D: P.Dp, mu0: eV(-P.Ev), cRef: m3(P.gv) },
          'V+': { D: DI, mu0: 0, cRef: m3(P.N0) },
        },
      },
      'spiro-OMeTAD': { epsr: P.epsH, species: { 'h+': { D: P.DH, mu0: eV(-P.EvH), cRef: m3(P.gvH) } } },
    },
    regions: [
      { material: 'TiO₂', length: P.bE, fixedCharge: FARADAY * m3(P.dE) },
      { material: 'MAPbI₃', length: P.b, fixedCharge: -FARADAY * m3(P.N0), c0: { 'V+': m3(P.N0) } },
      { material: 'spiro-OMeTAD', length: P.bH, fixedCharge: -FARADAY * m3(P.dH) },
    ],
    interfaces: [
      // Recombination through interface states (SRH), saturating at one carrier's capture: holes
      // reaching the TiO₂ recombine at v_p however many electrons wait there.
      { dipole: 0, species: { 'e-': 'equilibrium' }, reactions: [{ equation: 'e-(left) + h+ = 0', srh: { vn: P.vnE, vp: P.vpE, n1: m3((P.dE / n0) * ni) } }] },
      { dipole: 0, species: { 'h+': 'equilibrium' }, reactions: [{ equation: 'e- + h+(right) = 0', srh: { vn: P.vnH, vp: P.vpH, n1: m3(ni) } }] },
    ],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { ...drive(V, 'perovskiteCell: V'), terminal: 'h+', species: { 'h+': 'equilibrium' }, phi: 'bulk' },
    },
    bulkReactions: [
      ...(bulkSRH ? [{ equation: 'e- + h+ = 0', srh: { 'MAPbI₃': { tauN: P.tn, tauP: P.tp } } }] : []),
      ...(light > 0 ? [photogeneration({ material: 'MAPbI₃', flux: m3(P.Fph * light), alpha: P.alpha, mu: eV(3), from: P.bE, to: P.bE + P.b })] : []),
    ],
    grid: { hmin: 0.05e-9, hmax: 5e-9, ratio: 1.15 },
  };
}

/**
 * IonMonger's J–V scans of `perovskiteCell()` (IonMonger 2.0, commit 1ae4f9e, run in GNU Octave):
 * steady state at 1.2 V under light, then 1.2 → 0 → 1.2 V at `rate` (V/s). `full` is the cell
 * as given (bulk SRH on); `noBulkSRH` is `perovskiteCell({ bulkSRH: false })`. Per scan: J0, the
 * current at 0 V (A/m²); each sweep's maximum power Pmax (W/m², the largest V·J among IonMonger's
 * 201 output points per loop) and Voc (V, interpolated); and the hysteresis index
 * hi = (Pmax,rev − Pmax,fwd)/Pmax,rev. `hysteresis()` reads the same from a run's trace.
 * driftlet reproduces these to 0.15 % in hi, 0.5 W/m² in Pmax and 2 mV in Voc (its test suite).
 */
export const PEROVSKITE_SCANS = Object.freeze({
  source: 'IonMonger 2.0 (github.com/PerovskiteSCModelling/IonMonger, commit 1ae4f9e), default parameters',
  full: [
    { rate: 0.001, J0: 220.66, rev: { Pmax: 201.95, Voc: 1.09296 }, fwd: { Pmax: 200.42, Voc: 1.08517 }, hi: 0.00756 },
    { rate: 0.01, J0: 220.89, rev: { Pmax: 210.97, Voc: 1.11157 }, fwd: { Pmax: 195.18, Voc: 1.06566 }, hi: 0.07483 },
    { rate: 0.1, J0: 222.28, rev: { Pmax: 216.46, Voc: 1.11698 }, fwd: { Pmax: 144.85, Voc: 1.01277 }, hi: 0.33083 },
    { rate: 1, J0: 223.01, rev: { Pmax: 217.02, Voc: 1.11712 }, fwd: { Pmax: 210.91, Voc: 1.08238 }, hi: 0.02815 },
    { rate: 10, J0: 223.01, rev: { Pmax: 217.04, Voc: 1.11708 }, fwd: { Pmax: 216.55, Voc: 1.11509 }, hi: 0.00225 },
    { rate: 100, J0: 223.06, rev: { Pmax: 217.09, Voc: 1.11708 }, fwd: { Pmax: 216.93, Voc: 1.11693 }, hi: 0.0007 },
    { rate: 1000, J0: 223.54, rev: { Pmax: 217.58, Voc: 1.11714 }, fwd: { Pmax: 216.48, Voc: 1.11699 }, hi: 0.00505 },
  ],
  noBulkSRH: [
    { rate: 0.001, J0: 220.96, rev: { Pmax: 203.06, Voc: 1.0957 }, fwd: { Pmax: 201.54, Voc: 1.08762 }, hi: 0.00751 },
    { rate: 0.01, J0: 221.17, rev: { Pmax: 214.35, Voc: 1.1168 }, fwd: { Pmax: 196.01, Voc: 1.06737 }, hi: 0.08556 },
    { rate: 0.1, J0: 222.43, rev: { Pmax: 221.67, Voc: 1.12199 }, fwd: { Pmax: 145.73, Voc: 1.01331 }, hi: 0.3426 },
    { rate: 1, J0: 223.02, rev: { Pmax: 222.39, Voc: 1.12257 }, fwd: { Pmax: 213.41, Voc: 1.08376 }, hi: 0.04038 },
    { rate: 10, J0: 223.02, rev: { Pmax: 222.54, Voc: 1.12263 }, fwd: { Pmax: 221.77, Voc: 1.11951 }, hi: 0.00348 },
    { rate: 100, J0: 223.07, rev: { Pmax: 222.61, Voc: 1.12264 }, fwd: { Pmax: 222.38, Voc: 1.1223 }, hi: 0.00103 },
    { rate: 1000, J0: 223.55, rev: { Pmax: 223.12, Voc: 1.12271 }, fwd: { Pmax: 221.98, Voc: 1.12252 }, hi: 0.00513 },
  ],
});

/**
 * A J–V loop's figures of merit, from a transient's trace (`{ t, voltage, current }`, as a
 * recorder or advance() gives it): the loop splits at time `turn` into the sweep before (`rev`)
 * and after (`fwd`). Each sweep's maximum power V·I (W/m², the largest among the trace's points)
 * and its open-circuit voltage (where I changes sign, interpolated linearly); and the
 * hysteresis index (Pmax,rev − Pmax,fwd)/Pmax,rev. Power out is positive, so a cell whose
 * current runs the other way should have it negated first.
 * @param {{ t: ArrayLike<number>, voltage: ArrayLike<number>, current: ArrayLike<number> }} trace
 * @param {number} turn s
 */
export function hysteresis(trace, turn) {
  const out = {};
  for (const [name, on] of [['rev', (t) => t <= turn], ['fwd', (t) => t >= turn]]) {
    const pts = [];
    for (let k = 0; k < trace.t.length; k++) if (on(trace.t[k]) && Number.isFinite(trace.current[k])) pts.push([trace.voltage[k], trace.current[k]]);
    pts.sort((a, b) => a[0] - b[0]);
    const k = pts.findIndex(([, I], j) => j > 0 && pts[j - 1][1] > 0 && I <= 0);
    const Voc = k > 0 ? pts[k - 1][0] + (pts[k - 1][1] * (pts[k][0] - pts[k - 1][0])) / (pts[k - 1][1] - pts[k][1]) : NaN;
    out[name] = { Pmax: Math.max(...pts.map(([V, I]) => V * I)), Voc };
  }
  return { ...out, hi: (out.rev.Pmax - out.fwd.Pmax) / out.rev.Pmax };
}
