// Membranes, junctions and mixed conductors, their features combined at random: a family for
// bench/combos (see its README). Each case is one of seven kinds, and within a kind every feature
// is chosen on its own wherever the two make sense together:
//
// - cell: a membrane face (capacitive; GHK permeabilities, or Hodgkin–Huxley gates on permeability
//   or conductance links, from the kit's hodgkinHuxley()), a Na⁺/K⁺ pump, a surface charge, a
//   Gouy–Chapman layer, double layers resolved, flow and eddy mixing beside it, a buffer reacting
//   in the cytoplasm, Debye–Hückel activities; the cytoplasm ending in a bath held (a voltage
//   clamp), floating or driven by a current (a current clamp), a wall, a K⁺ electrode, or a second
//   membrane with the pump (an epithelium); a series R.
// - axon: a membrane port along the axoplasm (a capacitance, HH's gated conductances or a passive
//   leak, its outside as a bath), myelinated (nodes and internodes, a port each), a stimulus
//   injector, Debye–Hückel activities; the ends sealed, clamped by a bath, a K⁺ electrode, or a
//   floating pipette.
// - junction: liquid junctions stepped or graded, ion exchangers (Donnan), ions crowding on an
//   exchanger's sites (lattice), a bipolar membrane with water splitting, a dialysis membrane
//   (permeabilities at a face), a resolved charged pore, an exchanger bead (a sphere) in a shell of
//   solution; flow and mixing; Debye–Hückel; baths held (behind a series R), floating or driven by
//   a current, a bath right on an exchanger.
// - compartments: 2–4 solutions in a row, each face a wall blocking some species or a membrane,
//   complexes and an isomer reacting in the bulk, so that conserved combinations overlap; a bath,
//   a wall or a K⁺ electrode at the far end.
// - host: insertion hosts (OCV table or Redlich–Kister) with an electrolyte, Li⁺ crossing by a
//   transfer reaction or freely, a copper collector region, a second host (a full cell, held or at
//   open circuit), a Li⁺ bath or a lithium contact.
// - island: a silver island between two solutions (a bipolar electrode), its faces neutral or
//   capacitive beside the conductor (with a Gouy–Chapman layer, or the solution resolved); silver
//   electrodes or baths at the ends; a silver wire as a reference port at no current.
// - miec: a mixed ionic–electronic conductor (Ag⁺ and e⁻ on a fixed background) behind a silver
//   contact, ending at an electron-only contact (Hebb–Wagner) or at a solution, Ag⁺ crossing freely
//   or by a transfer reaction, read by a bath or a silver electrode.
//
// The plan: plan.V is absolute: each kind estimates its rest (GHK, chord conductances, Planck,
// Nernst, an OCV) and adds a small bias. plan.equilibrium is set only where every steady state of
// the definition at its drives is one (one bath and nothing driven, two alike, a closed cell
// without pump or flow, a host fed by one bath, an electrode held at its rest unrounded);
// plan.unique only where no amount is closed, no terminal floats or is driven, and every gated
// membrane's steady current–voltage curve falls everywhere (checked numerically below at the
// solutions' compositions: a Na⁺ channel's negative slope could hold two states). transientOnly
// marks a gated membrane not held to a voltage whose steady state is linearly unstable (it fires on
// its own: the transient is a limit cycle), by Routh's test on the space-clamped patch.
//
// Steered around (stated limits, each refused with its reason): a permeability link with
// Debye–Hückel ions, a Gouy–Chapman layer beside them, an electrode port in a resolved solution.

import { FARADAY, GAS_CONSTANT, EPS0 } from '../../../src/index.js';
import { aqueous, bath, layer, build, ohmic, metal, hodgkinHuxley, injector, pulse, planck, IONS } from '../../../src/kit.js';

const mV = 1e-3;
const VTof = (T) => (GAS_CONSTANT * T) / FARADAY;
const round = (x, d = 4) => +x.toPrecision(d);

// --- small physics, for the plans' estimates (never the solver's job)
const rate = (rt, V) => {
  const y = (V - rt.midpoint) / rt.scale;
  if (rt.type === 'exp') return rt.rate * Math.exp(y);
  if (rt.type === 'sigmoid') return rt.rate / (1 + Math.exp(-y));
  return Math.abs(y) < 1e-9 ? rt.rate : (rt.rate * y) / -Math.expm1(-y);
};
const xinf = (g, V) => rate(g.alpha, V) / (rate(g.alpha, V) + rate(g.beta, V));
const opened = (gates, which, V) => (which ? Object.entries(which).reduce((p, [k, n]) => p * xinf(gates[k], V) ** n, 1) : 1);
const slowestGate = (gates) => {
  let t = 0;
  for (const g of Object.values(gates ?? {})) for (let V = -0.12; V <= 0.06; V += 0.005) t = Math.max(t, 1 / (rate(g.alpha, V) + rate(g.beta, V)));
  return t;
};
// A species' level as a voltage, (μ° + RT ln(c/c_ref))/(zF), less φ: what a bath referenced to it reads.
const offset = (mu0, z, c, T, cRef = 1000) => (mu0 + GAS_CONSTANT * T * Math.log(c / cRef)) / (z * FARADAY);
const ionMu0 = (name) => IONS[name]?.mu0 ?? 0;
const zOf = (name) => IONS[name]?.z ?? { 'A-': -1, KA: 0 }[name];
// Steady current toward +x through a membrane's links at φ_R − φ_L = V (A/m²), the solutions'
// concentrations held: GHK for a permeability, G(V_L − V_R) for a conductance, gates at x∞(V).
function membraneCurrent(links, gates, cL, cR, V, T) {
  const VT = VTof(T);
  let I = 0;
  for (const [s, link] of Object.entries(links)) {
    if (typeof link !== 'object') continue;
    const z = zOf(s), o = opened(gates, link.gates, V);
    if (link.type === 'permeability') {
      const u = (z * V) / VT, P = link.P * o;
      const N = Math.abs(u) < 1e-9 ? P * (cL[s] - cR[s]) : (P * u * (cL[s] - cR[s] * Math.exp(u))) / Math.expm1(u);
      I += z * FARADAY * N;
    } else I += link.G * o * (-V + (VT / z) * Math.log(cL[s] / cR[s]));
  }
  return I;
}
// Where that current vanishes (a rest), by bisection from the most negative crossing; and whether
// the curve rises everywhere over a wide range (one steady state per drive).
function membraneRest(links, gates, cL, cR, T) {
  const f = (V) => membraneCurrent(links, gates, cL, cR, V, T);
  let a = -0.25, fa = f(a);
  for (let b = -0.245; b <= 0.25; b += 0.005) {
    const fb = f(b);
    if (Math.sign(fb) !== Math.sign(fa)) {
      let lo = b - 0.005, hi = b;
      for (let k = 0; k < 60; k++) {
        const m = (lo + hi) / 2;
        if (Math.sign(f(m)) === Math.sign(f(lo))) lo = m;
        else hi = m;
      }
      return (lo + hi) / 2;
    }
    [a, fa] = [b, fb];
  }
  return 0;
}
// Monotone: I toward +x falls as φ_R − φ_L rises (outward-positive current rises with V_m).
function monotone(links, gates, cL, cR, T) {
  let prev = membraneCurrent(links, gates, cL, cR, -0.25, T);
  for (let V = -0.249; V <= 0.25; V += 0.001) {
    const I = membraneCurrent(links, gates, cL, cR, V, T);
    if (!(I < prev)) return false;
    prev = I;
  }
  return true;
}
// Is the space-clamped membrane's steady state at V0 stable? Its state is (V, each gate), with
// C dV/dt = I(V, x) + I_ext(V) (I toward +x, inward, raising V = φ_R − φ_L) and HH kinetics for
// the gates, the solutions' concentrations held. Linearised there, the Jacobian's characteristic
// polynomial (Faddeev–LeVerrier) must pass Routh's test: a Hopf point past which HH fires on its
// own (a current clamp, a clamp behind a resistance) makes the transient a limit cycle, not the
// steady state. extSlope is dI_ext/dV (−1/R behind a resistance, 0 for a current).
function currentWith(links, x, cL, cR, V, T) {
  const VT = VTof(T);
  let I = 0;
  for (const [s, link] of Object.entries(links)) {
    if (typeof link !== 'object') continue;
    const z = zOf(s), o = link.gates ? Object.entries(link.gates).reduce((p, [k, n]) => p * x[k] ** n, 1) : 1;
    if (link.type === 'permeability') {
      const u = (z * V) / VT, P = link.P * o;
      I += z * FARADAY * (Math.abs(u) < 1e-9 ? P * (cL[s] - cR[s]) : (P * u * (cL[s] - cR[s] * Math.exp(u))) / Math.expm1(u));
    } else I += link.G * o * (-V + (VT / z) * Math.log(cL[s] / cR[s]));
  }
  return I;
}
function stableAt(links, gates, cL, cR, T, C, extSlope, V0) {
  if (!gates) return true;
  const names = Object.keys(gates), n = names.length + 1, h = 1e-7;
  const x0 = Object.fromEntries(names.map((k) => [k, xinf(gates[k], V0)]));
  const J = Array.from({ length: n }, () => new Array(n).fill(0));
  J[0][0] = ((currentWith(links, x0, cL, cR, V0 + h, T) - currentWith(links, x0, cL, cR, V0 - h, T)) / (2 * h) + extSlope) / C;
  names.forEach((k, j) => {
    const up = { ...x0, [k]: x0[k] + 1e-7 }, dn = { ...x0, [k]: x0[k] - 1e-7 };
    J[0][j + 1] = (currentWith(links, up, cL, cR, V0, T) - currentWith(links, dn, cL, cR, V0, T)) / 2e-7 / C;
    const g = gates[k], a = (V) => rate(g.alpha, V), b = (V) => rate(g.beta, V);
    J[j + 1][0] = ((a(V0 + h) - a(V0 - h)) / (2 * h)) * (1 - x0[k]) - ((b(V0 + h) - b(V0 - h)) / (2 * h)) * x0[k];
    J[j + 1][j + 1] = -(a(V0) + b(V0));
  });
  // det(sI − J) = s^n + c1 s^(n−1) + … + cn
  const mul = (A, B) => A.map((row) => B[0].map((_, j) => row.reduce((t, v, k) => t + v * B[k][j], 0)));
  let M = J.map((row) => row.map(() => 0));
  const coef = [1];
  for (let k = 1; k <= n; k++) {
    M = mul(J, M).map((row, i) => row.map((v, j) => v + (i === j ? coef[k - 1] : 0))); // M_k = J M_{k−1} + c_{k−1} I
    coef.push(-mul(J, M).reduce((t, row, i) => t + row[i], 0) / k); // c_k = −tr(J M_k)/k
  }
  // Routh: the first column positive throughout
  const rows = [coef.filter((_, i) => i % 2 === 0), coef.filter((_, i) => i % 2 === 1)];
  while (rows.at(-1).length && rows.length <= n) {
    const [a, b] = rows.slice(-2), next = [];
    if (!(b[0] > 0)) return false;
    for (let i = 0; i + 1 < a.length; i++) next.push((b[0] * a[i + 1] - a[0] * (b[i + 1] ?? 0)) / b[0]);
    rows.push(next);
  }
  return rows.every((rw) => !rw.length || rw[0] > 0);
}
// Where I(V) + I_ext(V) vanishes: every crossing over a wide range.
function crossings(f) {
  const out = [];
  let a = -0.25, fa = f(a);
  for (let b = -0.249; b <= 0.25; b += 0.001) {
    const fb = f(b);
    if (Math.sign(fb) !== Math.sign(fa)) {
      let lo = a, hi = b;
      for (let k = 0; k < 50; k++) {
        const m = (lo + hi) / 2;
        if (Math.sign(f(m)) === Math.sign(f(lo))) lo = m;
        else hi = m;
      }
      out.push((lo + hi) / 2);
    }
    [a, fa] = [b, fb];
  }
  return out;
}
const debye = (c, names, epsr, T) => Math.sqrt((epsr * EPS0 * GAS_CONSTANT * T) / (FARADAY * FARADAY * names.reduce((t, n) => t + zOf(n) ** 2 * (c[n] ?? 0), 0)));

// Water with the named ions and an impermeant anion A⁻ (and KA, a neutral complex, if asked).
function water(names, { epsr = 0, material = 'water', T, DA = 0.5e-9, complex = false } = {}) {
  const w = aqueous(names, { epsr, material, ...(T >= 273.15 && T <= 323.15 ? { T } : {}) });
  if (!w.species.some((s) => s.name === 'A-')) w.species.push({ name: 'A-', z: -1, cRef: 1000 });
  w.materials[material].species['A-'] = { D: DA, mu0: 0 };
  if (complex) {
    w.species.push({ name: 'KA', z: 0, cRef: 1000 });
    w.materials[material].species.KA = { D: 0.4e-9, mu0: complex };
  }
  return w;
}
// Debye–Hückel activity coefficients on a material's charged species (a 4 Å ion size)
const dhOn = (material, charged) => (material.statistics = [{ type: 'debye-huckel', species: [...charged], epsr: 78.3, a: 4e-10 }]);
const wall = (names) => ({ species: Object.fromEntries(names.map((s) => [s, 'blocked'])), phi: 'neutral' });

// ------------------------------------------------------------------------------------------------
// cell: outside | membrane | cytoplasm | (wall, bath, K⁺ electrode, or membrane | serosa | bath)
function cell(r) {
  const T = r.pick([279.45, 298.15, 310.15]), VT = VTof(T), RT = GAS_CONSTANT * T;
  const features = [];
  const epsr = r.chance(0.15) ? 74 : 0;
  if (epsr) features.push('resolved');
  const OUT = { 'Na+': r.u(100, 460), 'K+': r.u(3, 20), 'A-': r.u(1, 30) };
  OUT['Cl-'] = OUT['Na+'] + OUT['K+'] - OUT['A-'];
  const IN = { 'Na+': r.u(5, 50), 'K+': r.u(100, 400), 'Cl-': r.u(5, 60) };
  IN['A-'] = IN['Na+'] + IN['K+'] - IN['Cl-'];
  const right = r.pick(['clamp', 'clamp', 'float', 'current', 'closed', 'electrode', 'epithelium']);
  features.push(right);

  // The membrane's channels.
  const mode = r.pick(['ghk', 'ghk', 'hh', 'linear']);
  let links, gates;
  if (mode === 'ghk') {
    links = {};
    for (const [s, lo, hi] of [['K+', -9, -7.5], ['Na+', -10.5, -8.5], ['Cl-', -10, -8]]) links[s] = r.chance(0.15) ? 'blocked' : { type: 'permeability', P: r.log(lo, hi) };
    if (Object.values(links).every((l) => l === 'blocked')) links['K+'] = { type: 'permeability', P: 1e-8 };
  } else {
    features.push(mode === 'hh' ? 'HH-gates' : 'HH-linear');
    const g = { Na: 1200 * r.log(-0.5, 0.3), K: 360 * r.log(-0.5, 0.3), leak: 3 * r.log(-0.5, 0.5) };
    let hh;
    try {
      hh = mode === 'hh' ? hodgkinHuxley({ inside: IN, outside: OUT, T, g }) : hodgkinHuxley({ T, g, linear: true });
    } catch {
      hh = hodgkinHuxley({ inside: IN, outside: OUT, T, g, at: -60 * mV });
    }
    ({ gates, species: links } = hh);
  }
  const permeant = (s) => typeof links[s] === 'object';
  const pump = right !== 'epithelium' && permeant('Na+') && permeant('K+') && r.chance(0.35);
  // The Na⁺/K⁺ pump: Na⁺ out of the cell, K⁺ in. The cell is on the right of the outer membrane,
  // on the left of an epithelium's basolateral one.
  const pumpRx = (vmax, cell = 'right') => ({
    equation: cell === 'right' ? '3 Na+(right) + 2 K+(left) + ATP = 3 Na+(left) + 2 K+(right) + ADP' : '3 Na+(left) + 2 K+(right) + ATP = 3 Na+(right) + 2 K+(left) + ADP',
    fixed: { ATP: r.u(45e3, 55e3), ADP: 0 },
    vmax,
    K: { 'Na+': r.u(5, 20), 'K+': r.u(0.5, 3) },
  });
  const vmax = r.log(-8, -6);
  if (pump) features.push('pump');

  // The cytoplasm: plain, or with a buffer binding K⁺ (its own material, so a dipole at the face).
  const buffer = r.chance(0.25);
  const Kd = r.log(1.5, 2.5), kfBuf = r.log(-3, 0);
  if (buffer) {
    features.push('buffer');
    IN.KA = (IN['K+'] * IN['A-']) / Kd;
  }
  const muKA = ionMu0('K+') + 0 + RT * Math.log(Kd / 1000);
  const ions = ['Na+', 'K+', 'Cl-'];
  const library = [water(ions, { epsr, T }), ...(buffer ? [water(ions, { epsr, T, material: 'cyto', complex: muKA })] : [])];
  // Debye–Hückel activities: only where the membrane's links are conductances (a permeability link
  // needs ideal statistics on both sides, as its refusal says), and not with the K⁺ electrode,
  // whose rest is written from ideal levels.
  const dhCell = mode === 'linear' && right !== 'electrode' && right !== 'epithelium' && r.chance(0.5);
  if (dhCell) {
    for (const piece of library) for (const m of Object.values(piece.materials)) dhOn(m, [...ions, 'A-']);
    features.push('debye-huckel');
  }
  const cytoMat = buffer ? 'cyto' : 'water';

  const C = r.u(0.005, 0.02);
  const face = { phi: { type: 'capacitive', C }, species: { ...links, 'A-': 'blocked' } };
  if (gates) face.gates = gates;
  if (buffer) face.dipole = 0;
  if (pump) face.reactions = [pumpRx(vmax)];
  if (!epsr && !dhCell && r.chance(0.2)) {
    // (a diffuse layer takes ideal ions, or ions on one lattice: not Debye–Hückel's, docs/device.md)
    face.diffuse = { epsr: 78.4 };
    features.push('diffuse');
  }
  if (r.chance(0.2)) {
    face.sheetCharge = -r.log(-3, -1.7);
    features.push('sheetCharge');
  }

  // Lengths and the solution beside the membrane outside: flow, mixing.
  const Lo = epsr ? r.log(-6.7, -6) : r.log(-6, -4.7), Lc = epsr ? r.log(-6.7, -6) : r.log(-5.7, -4.7);
  const outside = { name: 'outside', c0: { ...OUT } };
  let v = 0;
  if (r.chance(0.2)) {
    v = (r.chance(0.5) ? 1 : -1) * r.log(-6, -4.5);
    outside.velocity = v;
    features.push('flow');
  }
  if (r.chance(0.15)) {
    outside.mixing = r.log(-9, -7.5);
    features.push('mixing');
  }

  // Terminals.
  const refL = r.pick(['Cl-', undefined]), refR = r.pick(['Cl-', undefined]);
  const offL = refL ? offset(ionMu0(refL), -1, OUT[refL], T) : 0;
  const offR = (c) => (refR ? offset(ionMu0(refR), -1, c[refR], T) : 0);
  const R = r.chance(0.2) ? r.log(-3, -1) : undefined;
  if (R) features.push('seriesR');
  const restPhi = membraneRest(links, gates, OUT, IN, T);
  let stack, plan;
  const inRight = { ...IN };
  const dV = r.u(0.002, 0.02);
  // Each permeant ion's permeability as its links hold it at rest (a conductance as the
  // permeability with its chord conductance at the cytoplasm's concentration, gates at x∞): a
  // compartment it alone feeds fills or empties in L/P.
  const effectiveP = (ls, V) => Object.entries(ls).filter(([, l]) => typeof l === 'object').map(([sp, l]) => (l.type === 'permeability' ? l.P : (l.G * RT) / (FARADAY * FARADAY * IN[sp])) * Math.max(opened(gates, l.gates, V), 1e-12));
  const slowExchange = (2 * Lc) / Math.min(...effectiveP(links, restPhi)); // (×2: coupled through V_m, it ran ~2× the bare L/P)
  const slowDiff = Math.max(Lo, Lc) ** 2 / 0.5e-9;
  const slows = [slowDiff, C / Math.max(1e-3, Math.abs(membraneCurrent(links, gates, OUT, IN, restPhi + 1e-3, T) - membraneCurrent(links, gates, OUT, IN, restPhi, T)) / 1e-3), slowestGate(gates)];
  const pumpSlow = pump ? (Lc * IN['Na+']) / vmax : 0; // where the pump alone sets a compartment
  if (buffer) slows.push(1 / (kfBuf * (IN['K+'] + IN['A-'])));
  if (v) slows.push(Lo / Math.abs(v));
  const mono = !gates || monotone(links, gates, OUT, IN, T);
  const drive = (V) => (R ? { V, R } : V);
  const lefts = (d) => bath({ ...OUT }, refL, d);
  const cyto = layer(cytoMat, Lc, { name: 'cytoplasm', c0: { ...IN } });
  const outL = layer('water', Lo, outside);
  const allNames = [...ions, 'A-', ...(buffer ? ['KA'] : [])];
  if (right === 'clamp' || right === 'float' || right === 'current') {
    // The cytoplasm's far side a bath: held (a voltage clamp), floating, or driven by a small
    // current (a current clamp, a few mV across the membrane's conductance).
    const rest = restPhi - offL + offR(inRight);
    const held = right === 'clamp';
    const I = round(r.u(-0.03, 0.03));
    const Vdef = round(rest + r.u(-0.03, 0.03));
    stack = [lefts(held ? 0 : drive(0)), outL, face, cyto, bath(inRight, refR, held ? drive(Vdef) : { I: right === 'current' ? I : 0 })];
    plan = { terminal: held ? 'right' : 'left', V: held ? round(rest + r.u(-0.05, 0.05)) : round(r.u(-0.05, 0.05)), equilibrium: false, unique: held && mono && !(gates && R) };
    // Gated channels not held to a voltage (behind a resistance, at a current, floating) may
    // fire on their own: where any steady state at the def's drive or at plan.V is unstable, the
    // transient is a limit cycle, so it isn't compared with the steady state.
    if (gates && !(held && !R)) {
      const at = (Vterm) => {
        const Vsrc = Vterm - offR(inRight) + offL; // the clamp's level as φ_in − φ_out
        const ext = held ? (V) => (Vsrc - V) / R : () => (right === 'current' ? I : 0);
        return crossings((V) => membraneCurrent(links, gates, OUT, IN, V, T) + ext(V)).every((V) => stableAt(links, gates, OUT, IN, T, C, held ? -1 / R : 0, V));
      };
      if (!(at(held ? Vdef : 0) && (!held || at(plan.V)))) {
        plan.transientOnly = true;
        features.push('unstable');
      }
    }
  } else if (right === 'closed') {
    stack = [lefts(drive(0)), outL, face, layer(cytoMat, Lc, { name: 'cytoplasm', c0: { ...IN } }), wall(allNames)];
    plan = { terminal: 'left', V: round(r.u(-0.05, 0.05)), equilibrium: !pump && !v, unique: false };
    slows.push(slowExchange, pumpSlow);
  } else if (right === 'electrode') {
    // An electrode reversible to K⁺ in the cytoplasm (a K⁺-selective wire): its level, held.
    // At rest (no current) every ion that crosses the membrane is level across it, so the
    // electrode reads the outside's K⁺ level against the left bath's reference.
    const rest = offset(ionMu0('K+'), 1, OUT['K+'], T) - offL;
    // (an equilibrium only if the def's drive is the rest itself, unrounded)
    const at = pump || v ? round(rest + r.u(-0.02, 0.02)) : rest;
    stack = [lefts(0), outL, face, cyto, { V: at, ...(R ? { R } : {}), terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' }];
    plan = { terminal: 'right', V: round(rest + r.u(-0.04, 0.04)), equilibrium: !pump && !v, unique: false };
    slows.push(slowExchange, pumpSlow);
  } else {
    // An epithelium: apical membrane (Na⁺ channels, some K⁺), the cell, a basolateral membrane
    // with K⁺ and Cl⁻ channels and the pump, then the serosal solution and its bath.
    const baso = {
      phi: { type: 'capacitive', C: r.u(0.005, 0.02) },
      species: { 'Na+': r.chance(0.5) ? 'blocked' : { type: 'permeability', P: r.log(-11, -9.5) }, 'K+': { type: 'permeability', P: r.log(-8.5, -7.5) }, 'Cl-': r.chance(0.5) ? 'blocked' : { type: 'permeability', P: r.log(-10, -8.5) }, 'A-': 'blocked' },
      reactions: [pumpRx(r.log(-7, -6), 'left')],
      ...(buffer ? { dipole: 0 } : {}),
    };
    if (buffer) baso.species.KA = 'blocked';
    face.species['Na+'] = { type: 'permeability', P: r.log(-10, -9) }; // ENaC: not voltage-gated
    const serosa = { ...OUT, 'Na+': OUT['Na+'] * r.u(0.8, 1.2) };
    serosa['Cl-'] = serosa['Na+'] + serosa['K+'] - serosa['A-'];
    const Ls = r.log(-6, -5);
    const held = r.chance(0.6);
    stack = [lefts(held ? 0 : drive(0)), outL, face, cyto, baso, layer('water', Ls, { name: 'serosa', c0: { ...serosa } }), bath(serosa, refR, held ? drive(round(offR(serosa) - offL)) : { I: 0 })];
    plan = { terminal: held ? 'right' : 'left', V: round((held ? offR(serosa) - offL : 0) + r.u(-0.02, 0.02)), equilibrium: false, unique: false };
    slows.push(slowExchange, Lc / Math.min(...effectiveP(baso.species, restPhi)), (Lc * IN['Na+']) / 1e-7);
    features.push('pump');
  }
  const Lmin = Math.min(Lo, Lc);
  const grid = epsr
    ? { hmin: Math.min(0.1e-9, debye(IN, ions, epsr, T) / 5), hmax: Math.min(20e-9, Lmin / 20), ratio: 1.15 }
    : { hmin: Math.max(20e-9, Lmin / r.u(50, 400)), hmax: Lmin / r.u(10, 30), ratio: r.u(1.1, 1.25) };
  const def = build({ T, library, stack, grid, ...(buffer ? { bulkReactions: [{ equation: 'K+ + A- = KA', kf: { cyto: kfBuf } }] } : {}) });
  return { def, params: { kind: 'cell', features, T, mode, right, OUT, IN, restPhi, Lo, Lc, C, refL, refR }, plan: { ...plan, dV, slow: Math.max(...slows) } };
}

// ------------------------------------------------------------------------------------------------
// axon: the axoplasm along x, its membrane a port (or a port per node and internode)
function axon(r) {
  const T = r.pick([279.45, 291.65, 298.15]), VT = VTof(T), RT = GAS_CONSTANT * T;
  const features = [];
  // Debye–Hückel in the axoplasm couples every ion's level to the others', so the port's bath then
  // lists them all, A⁻ too (docs/device.md, a membrane through a window), neutral.
  const dh = r.chance(0.2);
  const OUT = { 'Na+': r.u(400, 480), 'K+': r.u(10, 25), ...(dh ? { 'A-': r.u(1, 10) } : {}) };
  OUT['Cl-'] = OUT['Na+'] + OUT['K+'] - (OUT['A-'] ?? 0);
  const IN = { 'Na+': r.u(30, 70), 'K+': r.u(350, 420), 'Cl-': r.u(40, 120) };
  IN['A-'] = IN['Na+'] + IN['K+'] - IN['Cl-'];
  const ions = ['Na+', 'K+', 'Cl-'];
  const DA = r.log(-11, -10);
  const lib = aqueous(ions, { epsr: 0, material: 'axoplasm', T });
  lib.species.push({ name: 'A-', z: -1, cRef: 1000 });
  lib.materials.axoplasm.species['A-'] = { D: DA, mu0: 0 };
  if (dh) dhOn(lib.materials.axoplasm, [...ions, 'A-']);
  for (const s of ions) lib.materials.axoplasm.species[s].D *= 0.5;
  const a = r.log(-5.3, -3.4), am = 2 / a, Cm = r.u(0.008, 0.012);
  const passive = r.chance(0.2);
  const g = { Na: 1200 * r.log(-0.3, 0.2), K: 360 * r.log(-0.3, 0.2), leak: 3 * r.log(-0.5, 0.5) };
  const hh = hodgkinHuxley({ T, area: am, g });
  const species = passive
    ? { 'Na+': { type: 'conductance', G: 0.05 * am * r.log(-0.5, 0.5) }, 'K+': { type: 'conductance', G: 1 * am * r.log(-0.5, 0.5) }, 'Cl-': { type: 'conductance', G: 0.5 * am * r.log(-0.5, 0.5) } }
    : hh.species;
  const gates = passive ? undefined : hh.gates;
  features.push(passive ? 'passive' : 'HH', ...(dh ? ['debye-huckel'] : []));
  // The membrane as a face would see it, per area (G/area): for the rest and the slope.
  const perArea = Object.fromEntries(Object.entries(species).map(([s, l]) => [s, { ...l, G: l.G / am }]));
  // (as a face with the outside on its left: V = φ_in − φ_out, the gates' voltage)
  const restPhi = membraneRest(perArea, gates, OUT, IN, T);
  const mono = !gates || monotone(perArea, gates, OUT, IN, T);
  const gTot = Math.abs(membraneCurrent(perArea, gates, OUT, IN, restPhi + 1e-3, T) - membraneCurrent(perArea, gates, OUT, IN, restPhi, T)) / 1e-3;

  const myelin = r.chance(0.3);
  const regions = [], ports = [];
  let L, membraneFraction = 1; // of the axon's length, the part whose membrane passes ions
  const membranePort = (name, region) => ({ name, region, V: 0, area: am, capacitance: { C: Cm, zeroCharge: 0 }, bath: { c: { ...OUT } }, ...(gates ? { gates } : {}), species: JSON.parse(JSON.stringify(species)) });
  if (myelin) {
    features.push('myelin');
    const n = r.int(2, 4), Lint = r.u(100e-6, 500e-6), Lnode = r.u(1e-6, 3e-6);
    for (let k = 0; k < n; k++) {
      regions.push(layer('axoplasm', Lnode, { name: `node${k}`, c0: { ...IN }, grid: { minCells: 4 } }));
      ports.push(membranePort(`node${k}`, `node${k}`));
      if (k < n - 1) {
        regions.push(layer('axoplasm', Lint, { name: `internode${k}`, c0: { ...IN }, grid: { hmin: Lnode / 2, hmax: Lint / 20, ratio: 1.3 } }));
        ports.push({ name: `myelin${k}`, region: `internode${k}`, V: 0, area: am, capacitance: { C: Cm / r.u(50, 200), zeroCharge: 0 } });
      }
    }
    L = n * Lnode + (n - 1) * Lint;
    membraneFraction = (n * Lnode) / L;
  } else {
    L = r.log(-3.5, -2);
    regions.push(layer('axoplasm', L, { name: 'axon', c0: { ...IN }, grid: { hmin: L / r.u(40, 150), hmax: L / r.u(40, 150) } }));
    regions[0].grid.hmax = regions[0].grid.hmin;
    ports.push(membranePort('membrane', 'axon'));
  }
  const stim = r.chance(0.2);
  if (stim) {
    features.push('stimulus');
    const reg = regions[0].name, w = Math.min(regions[0].length, 1e-3);
    ports.push(injector({ name: 'stim', region: reg, from: 0, to: w, species: 'K+', I: pulse({ start: 1e-3, width: 0.2e-3, amplitude: r.u(5, 30) }) }));
  }
  const left = r.pick(['sealed', 'clamp', 'electrode']), right = r.pick(['sealed', 'sealed', 'pipette']);
  features.push(`left:${left}`, `right:${right}`);
  const sealed = wall([...ions, 'A-']);
  const offK = offset(ionMu0('K+'), 1, IN['K+'], T);
  const restAt = left === 'electrode' ? restPhi + offK : restPhi;
  const contacts = {
    left: left === 'sealed' ? sealed : left === 'clamp' ? bath({ ...IN }, undefined, round(restAt)) : { V: round(restAt), terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' },
    right: right === 'sealed' ? sealed : bath({ ...IN }, undefined, { I: 0 }),
  };
  const firstPort = ports[0].name;
  const def = { T, species: lib.species, materials: lib.materials, regions, contacts, ports };
  // How the composition relaxes, ion by ion: through the membrane (its amount per area of membrane,
  // a/2·c, over its chord conductance's F²/RT-scaled exchange, gates at x∞ at rest) and by
  // diffusion along the axon where an end is open (L²/D); in parallel where both act. A⁻ crosses
  // no membrane: only diffusion along, open end or not (a closed profile relaxes so too).
  const tauOf = (sp) => {
    const l = species[sp], D = lib.materials.axoplasm.species[sp].D;
    const g = (l.G / am) * Math.max(opened(gates, l.gates, restPhi), 1e-12);
    const tauEx = (2 * (a / 2) * IN[sp] * FARADAY * FARADAY) / (RT * g * membraneFraction), tauD = (L * L) / D;
    // (an end feeds it only as a bath: the K⁺ electrode feeds K⁺ alone)
    const fed = left === 'clamp' || right === 'pipette' || (left === 'electrode' && sp === 'K+');
    return fed ? 1 / (1 / tauEx + 1 / tauD) : tauEx;
  };
  const slow = Math.max(Cm / Math.max(gTot, 1e-3), slowestGate(gates), ...ions.map(tauOf), (L * L) / DA);
  // (a K⁺ electrode feeds no A⁻, so A⁻ is closed then; a pipette floats; a stimulus is driven)
  const unique = left === 'clamp' && right === 'sealed' && !stim && mono;
  const stable = crossings((V) => membraneCurrent(perArea, gates, OUT, IN, V, T)).every((V) => stableAt(perArea, gates, OUT, IN, T, Cm, 0, V));
  if (!stable) features.push('unstable');
  return {
    def,
    params: { kind: 'axon', features, T, a, L, IN, OUT, restPhi, left, right },
    plan: {
      terminal: left === 'sealed' ? firstPort : 'left',
      V: round(left === 'sealed' ? r.u(-0.03, 0.03) : restAt + r.u(-0.03, 0.03)),
      dV: r.u(0.001, 0.01),
      equilibrium: left === 'sealed' && right === 'sealed' && !stim,
      unique,
      slow,
      // a membrane patch that fires on its own at rest: the axon's transient is no steady state
      ...(stable ? {} : { transientOnly: true }),
    },
  };
}

// ------------------------------------------------------------------------------------------------
// junction: bath | solutions, exchangers, a dialysis membrane or a pore | bath
function junction(r) {
  const T = 298.15, features = [];
  const layout = r.pick(['free', 'free', 'exchanger', 'exchanger', 'bipolar', 'dialysis', 'pore', 'bead']);
  features.push(layout);
  const cats = layout === 'pore' ? ['K+'] : [...new Set([r.pick(['Na+', 'K+', 'Li+', 'H+']), ...(r.chance(0.4) ? [r.pick(['Na+', 'K+', 'Li+'])] : [])])];
  const ans = layout === 'pore' ? ['Cl-'] : [...new Set([r.pick(['Cl-', 'NO3-']), ...(r.chance(0.2) ? [r.pick(['Cl-', 'NO3-', 'SO42-'])] : [])])];
  const split = layout === 'bipolar' && r.chance(0.5);
  const names = [...cats, ...ans, ...(split ? ['H+', 'OH-'].filter((n) => !cats.includes(n)) : [])];
  const cw = split ? 1000 * Math.exp(-(-157.244e3 + 237.129e3) / (2 * GAS_CONSTANT * T)) : 0; // √K_w, mol/m³
  const side = (lo, hi) => {
    const c = {};
    for (const n of cats) c[n] = r.log(lo, hi);
    const plus = cats.reduce((t, n) => t + zOf(n) * c[n], 0), w = ans.map(() => r.u(0.2, 1)), wz = ans.reduce((t, n, k) => t + w[k] * -zOf(n), 0);
    ans.forEach((n, k) => (c[n] = (plus * w[k]) / wz));
    if (split) {
      // water at its own equilibrium, H⁺·OH⁻ = K_w (an acid's H⁺ sets OH⁻), still neutral
      if (c['H+'] === undefined) c['H+'] = c['OH-'] = cw;
      else {
        c['OH-'] = (cw * cw) / c['H+'];
        c[ans[0]] -= c['OH-'] / -zOf(ans[0]);
      }
    }
    return c;
  };
  const symmetric = layout !== 'bead' && r.chance(0.25);
  const cl = side(0, 3), cr = symmetric ? { ...cl } : side(0, 3);
  if (symmetric) features.push('symmetric');
  const epsr = layout === 'pore' ? 78.5 : 0;
  const lib = aqueous(names, { epsr });
  if (layout !== 'pore' && layout !== 'dialysis' && r.chance(0.15)) {
    dhOn(lib.materials.water, names);
    features.push('debye-huckel');
  }
  const libs = [lib];
  const Dscale = r.log(-1.3, -0.3);
  const X = r.log(1, 3.5);
  const crowded = r.chance(0.3);
  const exMat = (material) => {
    const m = aqueous(names, { epsr, material });
    for (const sp of Object.values(m.materials[material].species)) sp.D *= Dscale;
    // ions crowding on the exchanger's sites (a lattice), a few times what it holds
    if (crowded) m.materials[material].statistics = [{ type: 'lattice', species: [...names], cMax: r.u(1.5, 5) * (X + 2 * Math.max(...Object.values(cl), ...Object.values(cr))) }];
    libs.push(m);
  };
  const regions = [];
  const L = r.log(-5, -3.5);
  let v = 0;
  if (r.chance(0.25) && layout !== 'pore' && layout !== 'bead') {
    v = (r.chance(0.5) ? 1 : -1) * r.log(-7, -5);
    features.push('flow');
  }
  const flowAll = r.chance(0.5);
  const waterLayer = (len, c0, more = {}) => layer('water', len, { c0, ...(v ? { velocity: v } : {}), ...more });
  let bathMaterial = false;
  const faces = [];
  let geometry;
  if (layout === 'free') {
    if (r.chance(0.5)) {
      regions.push(waterLayer(L / 2, { ...cl }), waterLayer(L / 2, { ...cr }));
      features.push('stepped');
    } else {
      regions.push(waterLayer(L, Object.fromEntries(names.map((n) => [n, { x: [0, 0.4 * L, 0.6 * L, L], values: [cl[n], cl[n], cr[n], cr[n]] }]))));
      features.push('graded');
    }
  } else if (layout === 'exchanger' || layout === 'bipolar') {
    const Lm = r.log(-5.3, -4);
    const signs = layout === 'bipolar' ? (r.chance(0.5) ? [-1, 1] : [1, -1]) : [r.pick([-1, 1])];
    signs.forEach((s, k) => exMat(`ex${k}`));
    bathMaterial = layout === 'exchanger' && r.chance(0.3);
    if (!bathMaterial) regions.push(waterLayer(L / 2, { ...cl }));
    const mid = Object.fromEntries(names.map((n) => [n, (cl[n] + cr[n]) / 2])); // neutral, as both are
    signs.forEach((s, k) => {
      // c0 near its Donnan state: counter-ions raised by the fixed charge, co-ions as given.
      const c0 = { ...mid };
      const counter = s < 0 ? cats.find((n) => n !== 'H+') ?? cats[0] : ans[0];
      c0[counter] += X / Math.abs(zOf(counter));
      regions.push(layer(`ex${k}`, Lm, { c0, fixedCharge: s * X * FARADAY, ...(v && flowAll ? { velocity: v } : {}) }));
    });
    if (!bathMaterial) regions.push(waterLayer(L / 2, { ...cr }));
    if (layout === 'bipolar') features.push(split ? 'water-splitting' : 'no-splitting');
    if (crowded) features.push('crowded');
    if (bathMaterial) features.push('bath-on-exchanger');
  } else if (layout === 'dialysis') {
    regions.push(waterLayer(L / 2, { ...cl }), waterLayer(L / 2, { ...cr }));
    const capacitive = r.chance(0.5);
    faces[0] = { phi: capacitive ? { type: 'capacitive', C: r.log(-2, -0.5) } : 'neutral', species: Object.fromEntries(names.map((n) => [n, { type: 'permeability', P: r.log(-7, -5) }])) };
    if (capacitive) faces[0].dipole = 0;
    features.push(capacitive ? 'capacitive' : 'neutral');
  } else if (layout === 'bead') {
    // An ion-exchange bead (a sphere, sealed at its centre) in a shell of solution, one bath
    // outside: the solution starts at cl, the bath holds cr.
    exMat('ex0');
    if (crowded) features.push('crowded');
    const s = r.pick([-1, 1]), counter = s < 0 ? cats.find((n) => n !== 'H+') ?? cats[0] : ans[0];
    const c0 = { ...cl };
    c0[counter] += X / Math.abs(zOf(counter));
    regions.push(layer('ex0', r.log(-5.5, -4), { c0, fixedCharge: s * X * FARADAY }), waterLayer(L, { ...cl }));
    geometry = { type: 'spherical', r0: 0 };
  } else {
    // A charged pore, resolved (ε > 0): access layers, two halves (each its charge), access.
    const Lp = r.log(-8, -6.7);
    const pore = aqueous(['K+', 'Cl-'], { epsr, material: 'pore' });
    for (const sp of Object.values(pore.materials.pore.species)) sp.D *= 0.1;
    libs.push(pore);
    const X2 = r.chance(0.3) ? -X : r.chance(0.5) ? X : X * r.u(0.2, 1);
    const lD = debye({ 'K+': Math.min(cl['K+'], cr['K+']), 'Cl-': Math.min(cl['K+'], cr['K+']) }, ['K+', 'Cl-'], epsr, T);
    const access = Math.max(10e-9, 8 * lD);
    regions.push(layer('water', access, { c0: { ...cl } }), layer('pore', Lp / 2, { fixedCharge: -X * FARADAY }), layer('pore', Lp / 2, { fixedCharge: -X2 * FARADAY }), layer('water', access, { c0: { ...cr } }));
    faces[0] = { dipole: 0 };
    faces[2] = { dipole: 0 };
    if (X2 === -X) features.push('bipolar-pore');
  }
  if (r.chance(0.15) && layout !== 'pore' && layout !== 'bead' && !bathMaterial) {
    regions[0].mixing = r.log(-9, -7);
    features.push('mixing');
  }
  const ref = r.pick([ans[0], cats[0], undefined]);
  const off = (c) => (ref ? offset(ionMu0(ref), zOf(ref), c[ref], T) : 0);
  const ionTable = Object.fromEntries(names.map((n) => [n, { z: zOf(n), D: IONS[n].D }]));
  let restPhi;
  if (symmetric) restPhi = 0;
  else if (layout === 'exchanger' && X > 3 * Math.max(...Object.values(cl), ...Object.values(cr))) {
    // ideally selective: its counter-ion level across it
    const s = regions.find((g) => g.fixedCharge)?.fixedCharge < 0 ? cats[0] : ans[0];
    restPhi = -(VTof(T) / zOf(s)) * Math.log(cr[s] / cl[s]);
  } else restPhi = planck(cl, cr, ionTable, { T });
  const rest = restPhi + off(cr) - off(cl);
  // The right bath held (behind a series R, sometimes), floating, or driven by a current well
  // below the salt's limiting current, about F·D·c/L.
  const drive = layout === 'bead' ? 'held' : r.pick(['held', 'held', 'held', 'float', 'current']);
  const held = drive === 'held';
  if (!held) features.push(drive);
  const cSalt = Math.min(...cats.filter((n) => n !== 'H+' || !split).map((n) => Math.min(cl[n], cr[n])));
  const Ltot0 = regions.reduce((t, g) => t + g.length, 0);
  const iLim = (FARADAY * Math.min(...names.map((n) => IONS[n].D)) * cSalt) / Ltot0;
  const I = round((r.chance(0.5) ? 1 : -1) * r.u(0.05, 0.3) * iLim);
  const R = held && layout !== 'bead' && r.chance(0.2) ? r.log(-4, -2) : undefined;
  if (R) features.push('seriesR');
  const restV = symmetric ? 0 : round(rest);
  const leftBath = layout === 'bead' ? {} : bath({ ...cl }, ref, 0, bathMaterial ? { material: 'water' } : {});
  const rightBath = bath({ ...cr }, ref, held ? (R ? { V: restV, R } : restV) : { I: drive === 'current' ? I : 0 }, bathMaterial ? { material: 'water' } : {});
  const stack = [leftBath];
  regions.forEach((g, k) => {
    if (k > 0) stack.push(faces[k - 1] ?? {});
    stack.push(g);
  });
  stack.push(rightBath);
  const Ltot = regions.reduce((t, g) => t + g.length, 0);
  const cmin = Math.min(...Object.values(cl), ...Object.values(cr));
  const grid = layout === 'pore'
    ? { hmin: Math.min(0.02e-9, debye(cl, ['K+', 'Cl-'], epsr, T) / 8), hmax: 0.5e-9, ratio: 1.1 }
    : { hmin: Math.min(...regions.map((g) => g.length)) / r.u(200, 2000), hmax: Math.min(...regions.map((g) => g.length)) / r.u(15, 40), ratio: r.u(1.08, 1.2) };
  const def = build({
    T,
    library: libs,
    stack,
    grid,
    ...(geometry ? { geometry } : {}),
    ...(split ? { bulkReactions: [{ equation: 'H+ + OH- = H2O', fixed: { H2O: -237.129e3 }, kf: Object.fromEntries(['water', 'ex0', 'ex1'].map((m) => [m, 1.4e8])) }] } : {}),
  });
  const Dmin = Math.min(...names.map((n) => IONS[n].D)) * (layout === 'free' || layout === 'dialysis' ? 1 : layout === 'pore' ? 0.1 : Dscale);
  const slows = [(Ltot * Ltot) / Dmin];
  if (v) slows.push(Ltot / Math.abs(v));
  if (layout === 'dialysis') slows.push(L / Math.min(...Object.values(faces[0].species).map((l) => l.P)));
  return {
    def,
    params: { kind: 'junction', features, names, cl, cr, X, L, ref, restPhi, cmin, ...(drive === 'current' ? { I } : {}) },
    plan: {
      terminal: held ? 'right' : 'left',
      V: round(layout === 'bead' ? r.u(-0.1, 0.1) : held ? rest + r.u(-0.1, 0.1) : r.u(-0.05, 0.05)),
      dV: r.u(0.002, 0.02),
      // one bath (a bead), or two alike with nothing driven: every steady state is at rest
      equilibrium: layout === 'bead' || (symmetric && !v && drive !== 'current'),
      unique: held,
      slow: Math.max(...slows),
    },
  };
}

// ------------------------------------------------------------------------------------------------
// compartments: 2–4 solutions in a row, the faces between them each a wall blocking some species
// or a membrane (permeabilities, some blocked), complexes and an isomer reacting in the bulk, so
// conserved combinations overlap; a bath on the left, and on the right a bath, a wall, or a K⁺
// electrode at the bath's K⁺ level.
function compartments(r) {
  const T = 298.15, RT = GAS_CONSTANT * T, features = [];
  const all = [['K+ + A- = KA', 'K+', 'A-', 'KA'], ['Na+ + A- = NaA', 'Na+', 'A-', 'NaA'], ['B = C', 'B', 'C'], ['B + A- = BA-', 'B', 'A-', 'BA-']];
  const rxs = all.filter(() => r.chance(0.5));
  if (rxs.length === 0) rxs.push(all[0]);
  features.push(...rxs.map(([e]) => e));
  const zs = { 'K+': 1, 'Na+': 1, 'Cl-': -1, 'A-': -1, KA: 0, NaA: 0, B: 0, C: 0, 'BA-': -1 };
  const D = { 'K+': 1.957e-9, 'Na+': 1.334e-9, 'Cl-': 2.032e-9, 'A-': 0.7e-9, KA: 0.9e-9, NaA: 0.8e-9, B: 1e-9, C: 0.8e-9, 'BA-': 0.6e-9 };
  const mu0 = { 'K+': IONS['K+'].mu0, 'Na+': IONS['Na+'].mu0, 'Cl-': IONS['Cl-'].mu0, 'A-': -50e3, B: -5e3 };
  const K = {};
  for (const [, a, b, p] of rxs) {
    K[p ?? b] = p ? r.log(0.5, 2.5) : r.log(-1, 1); // a complex's K_d (mol/m³), the isomer's ratio
    mu0[p ?? b] = p ? mu0[a] + mu0[b] + RT * Math.log(K[p] / 1000) : mu0[a] - RT * Math.log(K[b]);
  }
  const names = [...new Set(['K+', 'Na+', 'Cl-', 'A-', ...rxs.flatMap(([, ...sp]) => sp)])];
  // A composition with every reaction at equilibrium (for the baths): the free species chosen,
  // the products from them, Cl⁻ for neutrality.
  const balanced = () => {
    const c = { 'K+': r.log(0.7, 2), 'Na+': r.log(0.7, 2), 'A-': r.log(0, 1.3) };
    if (names.includes('B')) c.B = r.log(0, 1.5);
    for (const [, a, b, p] of rxs) c[p ?? b] = p ? (c[a] * c[b]) / K[p] : c[a] * K[b];
    c['Cl-'] = names.reduce((t, n) => (n === 'Cl-' ? t : t + zs[n] * c[n]), 0);
    if (c['Cl-'] < 1) {
      c['K+'] += 1 - c['Cl-'];
      for (const [, a, b, p] of rxs) if (a === 'K+') c[p] = (c[a] * c[b]) / K[p];
      c['Cl-'] = names.reduce((t, n) => (n === 'Cl-' ? t : t + zs[n] * c[n]), 0);
    }
    return c;
  };
  const n = r.int(2, 4), regions = [], faces = [];
  for (let q = 0; q < n; q++) {
    const c0 = {};
    for (const s of names) if (s !== 'Cl-') c0[s] = r.log(0.5, 1.9);
    const q0 = names.reduce((t, s) => t + (s === 'Cl-' ? 0 : zs[s] * c0[s]), 0);
    if (q0 < 4) c0['K+'] += 4 - q0;
    c0['Cl-'] = Math.max(q0, 4);
    regions.push(layer('water', r.log(-5, -4.2), { name: `r${q}`, c0 }));
    if (q > 0) {
      if (r.chance(0.5)) {
        faces.push({ phi: 'neutral', species: Object.fromEntries(names.filter((s) => r.chance(0.35) && s !== 'Cl-').map((s) => [s, 'blocked'])) });
        features.push('wall');
      } else {
        faces.push({ phi: { type: 'capacitive', C: r.u(0.005, 0.05) }, species: Object.fromEntries(names.map((s) => [s, r.chance(0.3) ? 'blocked' : { type: 'permeability', P: r.log(-8, -6) }])) });
        features.push('membrane');
      }
    }
  }
  const left = balanced();
  const right = r.pick(['bath', 'closed', 'K+']);
  features.push(right);
  const offL = offset(mu0['Cl-'], -1, left['Cl-'], T);
  let rightContact, plan;
  if (right === 'bath') {
    const cr = balanced(), rest = offset(mu0['Cl-'], -1, cr['Cl-'], T) - offL;
    rightContact = bath(cr, 'Cl-', round(rest));
    // Each species' stretches open to a bath, so nothing closed: every face passes it, or it's
    // fed from both ends' stretches. A wall blocks a species where it says so; a membrane where
    // it's 'blocked'.
    const passes = (f, s) => f.species[s] !== 'blocked';
    const open = names.every((s) => {
      // the stretches between faces that block s: the first touches the left bath, the last the right
      let cut = 0;
      faces.forEach((f) => (cut += passes(f, s) ? 0 : 1));
      return cut <= 1;
    });
    plan = { terminal: 'right', V: round(rest + r.u(-0.03, 0.03)), equilibrium: false, unique: open };
  } else if (right === 'closed') {
    rightContact = wall(names);
    plan = { terminal: 'left', V: round(r.u(-0.05, 0.05)), equilibrium: true, unique: false };
  } else {
    // a K⁺ electrode, held at the left bath's K⁺ level: at rest, as a second bath alike would be
    const rest = offset(mu0['K+'], 1, left['K+'], T) - offL;
    rightContact = { V: rest, terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' };
    plan = { terminal: 'right', V: round(rest + r.u(-0.03, 0.03)), equilibrium: true, unique: false };
  }
  const stack = [bath(left, 'Cl-', 0)];
  regions.forEach((g, k) => {
    if (k > 0) stack.push(faces[k - 1]);
    stack.push(g);
  });
  stack.push(rightContact);
  const kf = r.log(-5, -2);
  const def = build({
    T,
    species: names.map((name) => ({ name, z: zs[name], cRef: 1000 })),
    materials: { water: { epsr: 0, species: Object.fromEntries(names.map((s) => [s, { D: D[s], mu0: mu0[s] }])) } },
    stack,
    bulkReactions: rxs.map(([equation]) => ({ equation, kf: { water: kf } })),
    grid: { hmin: 0.4e-6, hmax: 3e-6 },
  });
  const Ltot = regions.reduce((t, g) => t + g.length, 0);
  // The slowest relaxation: diffusion through the row; a reaction at its dilutest (k_f per
  // mol/m³, ~2 mol/m³; B = C is first order, k_f(1 + 1/K)); and each moiety (K, Na, A, B: its
  // forms) crossing each face by the forms that pass it, as the bath holds them: a compartment's
  // amount L·Σc over what crosses, Σ P·c at a membrane, Σ D·c/L at a wall. (A moiety carried out
  // only by a dilute form, A in BA⁻ through a membrane that blocks A⁻, KA and NaA, took 3e5 s
  // where L/P said 1e4.)
  const moieties = { K: ['K+', 'KA'], Na: ['Na+', 'NaA'], A: ['A-', 'KA', 'NaA', 'BA-'], B: ['B', 'C', 'BA-'], Cl: ['Cl-'] };
  const slows = [(Ltot * Ltot) / 0.6e-9, 1 / (kf * 2), K.C ? 1 / (kf * (1 + 1 / K.C)) : 0];
  faces.forEach((f, q) => {
    const Lq = Math.max(regions[q].length, regions[q + 1].length);
    for (const forms of Object.values(moieties)) {
      const here = forms.filter((s) => names.includes(s));
      if (!here.length) continue;
      const amount = Lq * here.reduce((t, s) => t + left[s], 0);
      const crossing = here.reduce((t, s) => {
        const l = f.species[s] ?? 'equilibrium';
        return t + (l === 'blocked' ? 0 : typeof l === 'object' ? l.P * left[s] : (D[s] * left[s]) / Lq);
      }, 0);
      if (crossing > 0) slows.push(amount / crossing);
    }
  });
  const slow = Math.max(...slows);
  return { def, params: { kind: 'compartments', features, names, K, n, right }, plan: { ...plan, dV: r.u(0.002, 0.02), slow } };
}

// ------------------------------------------------------------------------------------------------
// host: collector | (Cu |) insertion host | electrolyte | (bath, Li contact, or host | collector)
function host(r) {
  const T = 298.15, VT = VTof(T), RT = GAS_CONSTANT * T, features = [];
  const xs = Array.from({ length: 25 }, (_, k) => 0.02 + (0.96 * k) / 24);
  const makeHost = (E0) => {
    const cMax = r.u(2e4, 5e4), D = r.log(-15, -13);
    if (r.chance(0.5)) {
      const a = r.u(1, 2), b = r.u(0, 0.3);
      const E = (x) => E0 - a * VT * Math.log(x / (1 - x)) - b * x;
      return { cMax, D, E, stat: { type: 'insertion', species: ['Li+', 'e-'], cMax, ocv: { x: xs, E: xs.map(E), muRef: 0 } }, mu0: 0, how: 'ocv' };
    }
    const Om = r.u(-1.5, 1.5) * RT;
    const E = (x) => E0 - VT * Math.log(x / (1 - x)) + (2 * Om * x) / FARADAY;
    return { cMax, D, E, stat: { type: 'insertion', species: ['Li+', 'e-'], cMax, A: [Om] }, mu0: -FARADAY * E0, how: 'RK' };
  };
  const material = (h) => ({ epsr: 0, species: { 'Li+': { D: h.D, mu0: h.mu0, cRef: h.cMax }, 'e-': { D: 1e-4, mu0: 0, cRef: h.cMax } }, statistics: [h.stat] });
  const cathodeFirst = r.chance(0.5);
  const h1 = makeHost(cathodeFirst ? r.u(3.4, 4.2) : r.u(0.1, 0.6));
  features.push(`host:${h1.how}`);
  const x1 = r.u(0.15, 0.85);
  const c = r.u(300, 1500), Dli = r.log(-10, -9.5), DA = r.log(-9.7, -9.3);
  const materials = { h1: material(h1), elyte: { epsr: 0, species: { 'Li+': { D: Dli, mu0: 0 }, 'A-': { D: DA, mu0: 0 } } } };
  const species = [{ name: 'Li+', z: 1, cRef: 1000 }, { name: 'e-', z: -1 }, { name: 'A-', z: -1, cRef: 1000 }];
  const k0 = r.log(-6, -2);
  const transferFace = () => (r.chance(0.7) ? { species: { 'Li+': 'blocked' }, reactions: [{ equation: 'Li+(left) = Li+(right)', k0, alpha: r.u(0.3, 0.7) }] } : {});
  const f1 = transferFace();
  features.push(f1.reactions ? 'transfer' : 'free-Li');
  const Lh = r.log(-6.3, -5.5), Le = r.log(-5.3, -4.5);
  const stack = [];
  const copper = r.chance(0.3);
  if (copper) {
    stack.push(ohmic(0, ['e-']), layer('Cu', r.log(-6, -5)), {});
    features.push('copper');
  } else stack.push(ohmic(0, ['e-']));
  const graded = r.chance(0.3);
  stack.push(layer('h1', Lh, { name: 'host', c0: { 'Li+': x1 * h1.cMax } }), f1, layer('elyte', Le, { name: 'elyte', c0: graded ? { 'Li+': { x: [Lh, Lh + Le], values: [c * 0.7, c * 1.3] }, 'A-': { x: [Lh, Lh + Le], values: [c * 0.7, c * 1.3] } } : { 'Li+': c, 'A-': c } }));
  if (graded) features.push('graded');
  const right = r.pick(['bath', 'bath', 'lithium', 'full', 'full']);
  features.push(right);
  let plan;
  const restHalf = -h1.E(x1);
  const slows = [(Lh * Lh) / h1.D, (Le * Le) / Math.min(Dli, DA)];
  if (f1.reactions) slows.push((h1.cMax * x1 * (1 - x1) * Lh) / k0);
  if (right === 'bath') {
    stack.push(bath({ 'Li+': c, 'A-': c }, 'Li+', round(restHalf)));
    plan = { terminal: 'right', V: round(restHalf + r.u(-0.05, 0.05)), equilibrium: true, unique: true };
  } else if (right === 'lithium') {
    // a lithium metal contact, reversible to Li⁺ (its offset 0: μ_Li = 0 on the table convention)
    stack.push({ V: round(restHalf), terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' });
    plan = { terminal: 'right', V: round(restHalf + r.u(-0.05, 0.05)), equilibrium: true, unique: false };
  } else {
    const h2 = makeHost(cathodeFirst ? r.u(0.1, 0.6) : r.u(3.4, 4.2));
    materials.h2 = material(h2);
    const x2 = r.u(0.15, 0.85);
    const ocv = h2.E(x2) - h1.E(x1);
    const held = r.chance(0.5);
    const Lh2 = r.log(-6.3, -5.5);
    const f2 = transferFace();
    stack.push(f2, layer('h2', Lh2, { name: 'host2', c0: { 'Li+': x2 * h2.cMax } }), ohmic(held ? round(ocv) : { I: 0 }, ['e-']));
    plan = { terminal: held ? 'right' : 'left', V: round(held ? ocv + r.u(-0.03, 0.03) : r.u(-0.05, 0.05)), equilibrium: true, unique: false };
    slows.push((Lh2 * Lh2) / h2.D);
    if (f2.reactions) slows.push((h2.cMax * x2 * (1 - x2) * Lh2) / k0);
    features.push(held ? 'held' : 'open-circuit', `host2:${h2.how}`);
  }
  const lib = copper ? [metal('Cu')] : [];
  const def = build({ T, species, materials, library: lib, stack, grid: { hmin: r.log(-8.5, -7.5), hmax: r.log(-7.3, -6.5) } });
  return { def, params: { kind: 'host', features, x1, c, right, restHalf }, plan: { ...plan, dV: r.u(0.002, 0.01), slow: Math.max(...slows) } };
}

// ------------------------------------------------------------------------------------------------
// island: (Ag electrode | bath) | water | Ag island | water | (Ag electrode | bath), a reference wire
function island(r) {
  const T = 298.15, RT = GAS_CONSTANT * T, features = [];
  const support = r.chance(0.6);
  const names = ['Ag+', 'NO3-', ...(support ? ['K+'] : [])];
  const cAg = r.log(0, 2), cK = support ? r.log(1, 3) : 0;
  const c = { 'Ag+': cAg, ...(support ? { 'K+': cK } : {}), 'NO3-': cAg + cK };
  if (support) features.push('support');
  // The water strictly neutral (its faces neutral, capacitive, or capacitive with a Gouy–Chapman
  // layer below the grid), or its double layers resolved (capacitive faces, a fine grid).
  const resolved = r.chance(0.2);
  const lib = aqueous(names, { epsr: resolved ? 78.3 : 0 });
  if (resolved) features.push('resolved');
  const silver = { equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 } };
  const k0 = () => r.log(-5, -2);
  const rc = [];
  const face = () => {
    const kk = k0(), f = { reactions: [{ ...silver, k0: kk, alpha: r.u(0.3, 0.7) }] };
    if (resolved || r.chance(0.4)) {
      const C = r.u(0.1, 0.5);
      f.phi = { type: 'capacitive', C };
      f.zeroCharge = r.u(-0.3, 0.3);
      features.push('capacitive');
      if (!resolved && r.chance(0.4)) {
        f.diffuse = { epsr: 78.4 };
        features.push('diffuse');
      }
      // C·R_ct, R_ct = RT/(F² k0 (c/c_ref)^(1−α)) at its larger exponent's (the dilute side's)
      rc.push({ C, R: RT / (FARADAY * FARADAY * kk * (cAg / 1000) ** Math.max(f.reactions[0].alpha, 1 - f.reactions[0].alpha)) });
    }
    if (support) f.species = { 'K+': 'blocked' };
    return f;
  };
  const L1 = r.log(-5.5, -4), L2 = r.log(-5.5, -4);
  const ends = [r.pick(['bath', 'electrode']), r.pick(['bath', 'electrode'])];
  features.push(`ends:${ends.join('/')}`);
  const stack = [];
  if (ends[0] === 'bath') stack.push(bath({ ...c }, 'NO3-', 0));
  else stack.push(ohmic(0, ['e-']), layer('Ag', 1e-6), face());
  stack.push(layer('water', L1, { name: 'w1', c0: { ...c } }), face(), layer('Ag', r.log(-6, -4), { name: 'island' }), face(), layer('water', L2, { name: 'w2', c0: { ...c } }));
  // At rest both ends read the same thing (alike baths, or alike electrodes): 0 V. A bath against
  // an electrode reads NO₃⁻'s level against the metal's, where Ag⁺ is at the metal's level.
  const offAg = offset(IONS['Ag+'].mu0, 1, cAg, T), offN = offset(IONS['NO3-'].mu0, -1, cAg + cK, T);
  const rest = ends[0] === ends[1] ? 0 : ends[0] === 'electrode' ? offN - offAg : offAg - offN;
  if (ends[1] === 'bath') stack.push(bath({ ...c }, 'NO3-', rest));
  else stack.push(face(), layer('Ag', 1e-6), ohmic(rest, ['e-']));
  const wire = !resolved && r.chance(0.3); // (an electrode port needs a strictly neutral region: docs/device.md)
  const ports = wire ? [{ name: 'ref', region: 'w1', from: L1 * 0.4, to: L1 * 0.45, terminal: 'e-', I: 0, area: r.log(3, 5), reactions: [{ ...silver, k0: k0(), alpha: 0.5 }] }] : [];
  if (wire) features.push('reference-wire');
  const lD = debye(c, names, 78.3, T);
  const grid = resolved ? { hmin: lD / r.u(4, 10), hmax: Math.min(L1, L2) / r.u(15, 40), ratio: 1.12 } : { hmin: Math.min(L1, L2) / r.u(200, 1000), hmax: Math.min(L1, L2) / r.u(15, 40), ratio: 1.15 };
  const def = build({ T, library: [lib, metal('Ag')], stack, ports, grid });
  return {
    def,
    params: { kind: 'island', features, c, L1, L2 },
    // NO₃⁻ is blocked at every metal, so each solution's NO₃⁻ is closed unless a bath feeds it.
    plan: { terminal: 'right', V: round(rest + r.u(-0.1, 0.1)), dV: r.u(0.002, 0.02), equilibrium: true, unique: ends[0] === 'bath' && ends[1] === 'bath' && !wire, slow: Math.max((L1 + L2) ** 2 / IONS['Ag+'].D, rc.reduce((t, x) => t + x.C, 0) * rc.reduce((t, x) => t + x.R, 0)) }, // (Ag⁺ crosses the island: one mode spans both; the faces' charge through all their R_ct)
  };
}

// ------------------------------------------------------------------------------------------------
// miec: Ag contact | mixed conductor (Ag⁺, e⁻) | (electron-only contact, or solution | bath / Ag)
function miec(r) {
  const T = 298.15, RT = GAS_CONSTANT * T, features = [];
  const X = r.log(1.5, 3.5), ratio = r.log(-5, 0); // fixed charge (mol/m³); e⁻ against it
  const ce = X * ratio, cA = X + ce;
  const muA = -RT * Math.log((cA * ce) / 1e6); // Ag⁺ + e⁻ at Ag metal's μ (0): μ°_A + μ°_e = −RT ln(c_A c_e/c_ref²)
  const Lm = r.log(-6, -4.5);
  const DAg = r.log(-11, -9.5), De = r.log(-8, -6);
  const materials = { mixed: { epsr: 0, species: { 'Ag+': { D: DAg, mu0: muA, cRef: 1000 }, 'e-': { D: De, mu0: 0, cRef: 1000 } } } };
  const species = [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'e-', z: -1 }];
  const agContact = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'Ag+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
  const stack = [agContact(0), layer('mixed', Lm, { name: 'mixed', fixedCharge: -X * FARADAY, c0: { 'Ag+': cA, 'e-': ce } })];
  const right = r.pick(['electronic', 'solution', 'solution']);
  features.push(right);
  let rest = 0, unique = true;
  const slows = [(Lm * Lm) / DAg];
  const library = [];
  if (right === 'electronic') {
    stack.push(ohmic(0, ['e-']));
  } else {
    const c = r.log(0, 2), Lw = r.log(-5.5, -4.5);
    const w = aqueous(['Ag+', 'NO3-'], { epsr: 0 });
    library.push(w);
    const transfer = r.chance(0.5), kt = r.log(-5, -2);
    stack.push(transfer ? { species: { 'Ag+': 'blocked' }, reactions: [{ equation: 'Ag+(left) = Ag+(right)', k0: kt, alpha: 0.5 }] } : {});
    if (transfer) slows.push((Math.min(ce, cA) * Lm) / kt); // the conductor's chemical capacity through the transfer
    features.push(transfer ? 'transfer' : 'free-Ag');
    stack.push(layer('water', Lw, { c0: { 'Ag+': c, 'NO3-': c } }));
    const end = r.pick(['bath', 'silver']);
    features.push(end);
    if (end === 'bath') {
      rest = offset(IONS['NO3-'].mu0, -1, c, T) - offset(IONS['Ag+'].mu0, 1, c, T);
      stack.push(bath({ 'Ag+': c, 'NO3-': c }, 'NO3-', rest));
    } else {
      library.push(metal('Ag'));
      stack.push({ reactions: [{ equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: r.log(-5, -2), alpha: 0.5 }] }, layer('Ag', 1e-6), ohmic(0, ['e-']));
      unique = false; // NO₃⁻ is closed in the water
    }
    slows.push((Lw * Lw) / IONS['NO3-'].D);
  }
  const def = build({ T, species, materials, library, stack, grid: { hmin: Lm / r.u(200, 1000), hmax: Lm / r.u(15, 40), ratio: 1.15 } });
  return { def, params: { kind: 'miec', features, X, ratio, Lm, right }, plan: { terminal: 'right', V: round(rest + r.u(-0.1, 0.1)), dV: r.u(0.002, 0.02), equilibrium: true, unique, slow: Math.max(...slows) } };
}

const KINDS = { cell, axon, junction, compartments, host, island, miec };
const WEIGHTS = [['cell', 4], ['axon', 2], ['junction', 3], ['compartments', 2], ['host', 2], ['island', 1], ['miec', 1]];

// (for checking the stability test against Hodgkin and Huxley's own patch: rest at −65 mV stable,
// repetitive firing from ~0.1 A/m², depolarisation block past ~1.5 A/m²)
export const _test = { stableAt, crossings, membraneCurrent };

export default {
  name: 'membrane',
  make(r) {
    const total = WEIGHTS.reduce((t, [, w]) => t + w, 0);
    let u = r.u(0, total), kind = WEIGHTS[0][0];
    for (const [k, w] of WEIGHTS) {
      if (u < w) {
        kind = k;
        break;
      }
      u -= w;
    }
    const c = KINDS[kind](r);
    c.params.features = [kind, ...c.params.features];
    // Each kind's slow is a sum of physical times (L²/D, L/P, a moiety's amount over its crossing,
    // C·R, 1/k, the gates'); checked against advance() from plan.V to 10 × slow and a solve from
    // there (300 cases), the levels had settled to 1e-2 RT in all but those marked transientOnly,
    // the worst ~2× their estimate: hence the 2.
    c.plan.slow *= 2;
    c.def = structuredClone(c.def); // plain data, nothing shared (a kit helper's gates reused across ports)
    return c;
  },
};
