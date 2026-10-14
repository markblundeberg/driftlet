// The battery film's simulation, apart from its pictures (so it runs in node too): a toy
// thin-film lithium-ion cell, and the story told with it, one continuous run.

import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer } from '../src/kit.js';

// The cell: a toy thin-film lithium-ion cell, graphite | liquid electrolyte | iron phosphate (LFP),
// each host a regular solution (an ideal lattice gas with interaction Ω) against Li metal, Li⁺
// crossing each surface through a linear charge-transfer resistance.
export const RT = GAS_CONSTANT * 298.15;
export const GR = { cMax: 30000, Om: 1.0 * RT, E0: 0.1, D: 1.5e-15 };
export const LFP = { cMax: 22800, Om: 1.9 * RT, E0: 3.42, D: 5e-14 };
export const La = 2.2e-6, Ls = 3e-6, Lc = 2.5e-6, L = La + Ls + Lc, SALT = 1000;
export const X_EMPTY = { a: 0.02, c: 0.97 };
export const Q = LFP.cMax * Lc * FARADAY; // the cell's capacity, C/m² (the cathode's)
export const mu0 = (h) => -h.E0 * FARADAY + h.Om;
export const E = (h, x) => -(mu0(h) + RT * Math.log(x / (1 - x)) - 2 * h.Om * x) / FARADAY;
const host = (h) => ({
  epsr: 0,
  species: { 'Li+': { D: h.D, mu0: mu0(h), cRef: h.cMax }, 'e-': { D: 1e-6, mu0: 0, cRef: h.cMax } },
  statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: h.cMax, A: [h.Om] }],
});
const surface = { species: { 'Li+': { type: 'conductance', G: 100 } } };
export const cell = () => build({
  species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }, { name: 'PF6-', z: -1 }],
  materials: {
    graphite: host(GR), lfp: host(LFP),
    electrolyte: { epsr: 0, species: { 'Li+': { D: 3e-11, mu0: 0, cRef: 1000 }, 'PF6-': { D: 5e-11, mu0: 0, cRef: 1000 } } },
  },
  stack: [
    { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, // copper
    layer('graphite', La, { c0: { 'Li+': X_EMPTY.a * GR.cMax, 'e-': X_EMPTY.a * GR.cMax } }),
    surface,
    layer('electrolyte', Ls, { c0: { 'Li+': SALT, 'PF6-': SALT } }),
    surface,
    layer('lfp', Lc, { c0: { 'Li+': X_EMPTY.c * LFP.cMax, 'e-': X_EMPTY.c * LFP.cMax } }),
    { V: E(LFP, X_EMPTY.c) - E(GR, X_EMPTY.a), terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, // aluminium
  ],
  grid: { hmin: 10e-9, hmax: 80e-9, ratio: 1.1 },
});
export const CUTOFF = [2.8, 3.65];

// The story, one continuous run: phases of constant current (C-rate; + charges), each to its end.
export const PHASES = [
  { key: 'intro', secs: 2.2, C: 0, sim: 60 },
  { key: 'slow', secs: 6.5, C: 0.5, until: 'full' },
  { key: 'use', secs: 4.5, C: -0.5, until: 'empty' },
  { key: 'fast', secs: 5.5, C: 4, until: 'plating' },
  { key: 'ease', secs: 5, C: 4, control: true, until: 'full' },
  { key: 'end', secs: 3, C: 0, sim: 120 },
];

// A seeded random source, so a film renders the same every time.
export const rng = (seed) => () => ((seed = (seed + 0x6d2b79f5) | 0), (((seed ^ (seed >>> 15)) * (1 | seed)) >>> 0) / 4294967296);

// The run: the device, the current, the phase, and what the story needs to know.
export function makeRun() {
  const dev = new Device(cell());
  dev.solve();
  dev.set({ contacts: { right: { V: undefined, I: 0 } } });
  return { dev, I: 0, phase: 0, phaseT: 0, plated: false, margin: Infinity, history: [] };
}
export const sol = (run) => run.dev.solution();
// The Li in each part (mol/m²), and the plating margin: the anode against Li/Li⁺ at its surface,
// V_e⁻ in the graphite less V_Li⁺ in the electrolyte beside it (lithium metal forms below zero).
export function inventory(S) {
  let g = 0, e = 0, c = 0;
  for (let k = 1; k < S.x.length; k++) {
    if (S.region[k] !== S.region[k - 1]) continue;
    const a = 0.5 * (S.c['Li+'][k] + S.c['Li+'][k - 1]) * (S.x[k] - S.x[k - 1]), r = S.region[k];
    if (r === 0) g += a; else if (r === 1) e += a; else c += a;
  }
  const ia = S.region.lastIndexOf(0), ic = S.region.indexOf(2);
  return { g, e, c, soc: 1 - c / (LFP.cMax * Lc), margin: S.V['e-'][ia] - S.V['Li+'][ia + 1], xa: S.c['Li+'][ia] / GR.cMax, xc: S.c['Li+'][ic] / LFP.cMax };
}

// Advance the run toward simulated time t, in steps short enough to catch the phase's end; returns
// true when the phase has ended.
export function advance(run, t) {
  const ph = PHASES[run.phase];
  while (run.dev.solver.time < t - 1e-6) {
    const t0 = run.dev.solver.time;
    const C = Math.abs(run.I) / (Q / 3600);
    const dtMax = C > 0 ? 3600 / C / 300 : Infinity;
    // A step that won't go (the steep last millivolts before a cutoff) ends the phase there.
    const s = run.dev.advance(t0 + Math.min(t - t0, dtMax, run.dtV ?? Infinity));
    if (!s.converged) return true;
    const inv = inventory(sol(run)), V = s.terminalVoltage;
    // Each host's surface filling, as ln(x/(1−x)), which its voltage follows and which diverges as
    // the surface fills or empties: the next step moves it by at most 0.25 (about 6 mV).
    const u = [inv.xa, inv.xc].map((x) => Math.log(x / (1 - x))), dt = Math.max(s.time - t0, 1e-12);
    const rate = run.u ? Math.max(...u.map((v, j) => Math.abs(v - run.u[j]))) / dt : 0;
    run.u = u;
    run.dtV = rate > 0 ? Math.max(1e-6, 0.25 / rate) : Infinity;
    run.margin = inv.margin;
    if (inv.margin < 0) run.plated = true;
    if (ph.control) {
      // Charging as fast as keeps the anode 15 mV above plating: a gentle proportional controller,
      // in short steps.
      const Imax = (ph.C * Q) / 3600, k = Math.exp(Math.max(-0.1, Math.min(0.05, (inv.margin - 0.015) / 0.02)));
      setI(run, Math.min(Imax, run.I * k));
      run.dtV = Math.min(run.dtV, 3600 / ph.C / 1500);
    }
    const cutoff = (run.I > 0 && V >= CUTOFF[1]) || (run.I < 0 && V <= CUTOFF[0]);
    const end = cutoff || (ph.control && run.I < (0.25 * Q) / 3600) || (ph.until === 'plating' && inv.margin < 0.015);
    if (end) return true;
  }
  return ph.sim !== undefined && run.dev.solver.time >= run.phaseT + ph.sim - 1e-9;
}
export function setI(run, I) {
  run.I = I;
  run.dev.set({ contacts: { right: { I } } });
}
export function enter(run, k) {
  run.phase = k;
  run.phaseT = run.dev.solver.time;
  run.dtV = undefined;
  run.u = undefined;
  setI(run, (PHASES[k].C * Q) / 3600);
}

