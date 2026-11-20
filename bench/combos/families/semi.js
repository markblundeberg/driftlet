// Semiconductor devices that combine features at random: the family's generator (see ../README.md).
//
// Three shapes of device, each drawing on one pool of features, every feature on or off on its own
// where the physics allows it:
//   stack  a contact, 1–4 layers, a contact: homo- and heterojunctions (Si, Ge, GaAs, and AlGaAs and
//          InGaAs made from GaAs as docs/data.md does), faces aligned explicitly (an e⁻ step from
//          electron affinities, or vacuumDipole), alloys and hand-graded gaps, doping profiles;
//          contacts ohmic, ohmic with a minority surface velocity, selective, Schottky (pinned or
//          capacitive, emission velocities or equilibrium links), or a metal region behind a
//          capacitive face; a series R, a current drive, a base port (a 1D bipolar transistor), an
//          injector port, a geometry.
//   mos    1–2 layers under a gate: a capacitive gate contact, or an oxide and a metal, or a floating
//          gate between two oxides (wired to a port, or not); oxide charge, interface SRH, a port
//          holding the minority carrier beside the oxide (source and drain's stand-in).
//   tft    a film along x between source and drain, its gate a port's capacitance through the
//          window (one or two gates), optional doped source and drain regions, ε = 0 or not.
// Shared by all three: the temperature, Fermi–Dirac statistics, bulk SRH, band (radiative) and
// Auger recombination, trap levels (kit trap(), acceptor- or donor-like, one or two), light.
//
// A case is plain data, and every alignment is given: band offsets from electron affinities
// (Anderson's rule) with a jitter, Schottky barriers from a work function clamped into the gap,
// gate flat bands from the metal's work function. Kinetics are the same in every material of a case
// (build() refuses a rate keyed to an alloy's `to` that differs from its `from`).
import { units, FARADAY, GAS_CONSTANT, EPS0 } from '../../../src/index.js';
import {
  build, layer, alloy, ohmic, combine, semiconductor, metal, trap, recombination, photogeneration, injector,
  vacuumZeroCharge, vacuumDipole, SEMICONDUCTORS,
} from '../../../src/kit.js';

const eV = units.eV, cm3 = units.perCm3;
const WORK = { Al: 4.28, Ag: 4.26, Cu: 4.65, Au: 5.1, Pt: 5.65 }; // V (CRC, polycrystalline)
const OXIDES = { SiO2: 3.9, Al2O3: 9, HfO2: 20 };

// --- materials: the library piece, its electron affinity χ (V) and gap (eV), at T
function material(name, T, x) {
  const VT = (GAS_CONSTANT * T) / FARADAY;
  if (SEMICONDUCTORS[name]) {
    const lib = semiconductor(name, { T }), m = lib.materials[name];
    return { name, lib, chi: SEMICONDUCTORS[name].chi, Eg: m.species['h+'].mu0 / FARADAY, m };
  }
  // AlₓGa₁₋ₓAs and InₓGa₁₋ₓAs from GaAs's entry, as docs/data.md shows: the gap moved by the
  // alloy's shift from GaAs's at 300 K, N_c by the electron mass, ΔE_c a fixed share of ΔE_g.
  const lib = semiconductor('GaAs', { material: name, T }), m = lib.materials[name], EgG = m.species['h+'].mu0 / FARADAY;
  let dEg, me, share, epsr, mun;
  if (name === 'AlGaAs') (dEg = 1.247 * x), (me = 0.067 + 0.083 * x), (share = 0.62), (epsr = 13.1 - 3 * x), (mun = 8000 - 22000 * x + 10000 * x * x); // (Adachi's fit, x < 0.45)
  else (dEg = -(1.07 * x - 0.477 * x * (1 - x))), (me = 0.067 - 0.044 * x), (share = 0.6), (epsr = 13.1 + 1.5 * x), (mun = 8500 + 10000 * x);
  m.species['h+'] = { ...m.species['h+'], mu0: eV(EgG + dEg) };
  m.species['e-'] = { ...m.species['e-'], D: mun * Math.pow(T / 300, -1) * 1e-4 * VT, cRef: m.species['e-'].cRef * Math.pow(me / 0.067, 1.5) };
  m.epsr = epsr;
  return { name, lib, chi: 4.07 - share * dEg, Eg: EgG + dEg, m };
}
// Pairs an alloy may join (the two alike in all but their numbers).
const ALLOYS = [['Si', 'Ge'], ['Ge', 'Si'], ['GaAs', 'AlGaAs'], ['AlGaAs', 'GaAs'], ['GaAs', 'InGaAs'], ['InGaAs', 'GaAs'], ['AlGaAs', 'InGaAs']];
const SETS = [['Si'], ['Si'], ['Ge'], ['GaAs'], ['GaAs'], ['Si', 'Ge'], ['GaAs', 'AlGaAs'], ['GaAs', 'InGaAs'], ['GaAs', 'AlGaAs', 'InGaAs'], ['Ge', 'GaAs'], ['Si', 'GaAs']];

export default {
  name: 'semi',
  make(r) {
    const features = new Set(), on = (f) => features.add(f);
    const T = r.chance(0.7) ? 300 : r.pick([250, 275, 325, 350, 400]);
    if (T !== 300) on('T');
    const RT = GAS_CONSTANT * T, VT = RT / FARADAY;
    const topology = r.pick(['stack', 'stack', 'stack', 'mos', 'mos', 'tft']);
    on(topology);
    const set = topology === 'tft' ? [r.pick(['Si', 'Si', 'GaAs', 'Ge', 'InGaAs'])] : r.pick(SETS);
    const xAl = r.u(0.1, 0.4), xIn = r.u(0.05, 0.2);
    const M = Object.fromEntries(set.map((n) => [n, material(n, T, n === 'AlGaAs' ? xAl : xIn)]));
    const fd = r.chance(0.3);
    if (fd) {
      on('fermi-dirac');
      const species = r.chance(0.8) ? ['e-', 'h+'] : [r.pick(['e-', 'h+'])];
      for (const n of set) M[n].m.statistics = [{ type: 'fermi-dirac', species }];
    }
    const epsr0 = topology === 'tft' && r.chance(0.4);
    if (epsr0) (M[set[0]].m.epsr = 0), on('epsr0');
    const Nc = (n) => M[n].m.species['e-'].cRef, Nv = (n) => M[n].m.species['h+'].cRef;
    const ni = (n) => Math.sqrt(Nc(n) * Nv(n) * Math.exp(-M[n].Eg / VT));
    const Dmin = Math.min(...set.flatMap((n) => [M[n].m.species['e-'].D, M[n].m.species['h+'].D]));
    const epsOf = (n) => (M[n].m.epsr || 11) * EPS0;
    const debye = (n, N) => Math.sqrt((epsOf(n) * RT) / (FARADAY * FARADAY * Math.max(N, ni(n), cm3(1e13))));
    // Doping: without Fermi–Dirac statistics, held below a fifth of the band's density of states
    // (a page would turn the statistics on past it; GaAs's N_c is 4.7e17 cm⁻³).
    const dope = (n, type, lo = 14, hi = fd ? 19.7 : 18.7) => {
      let N = cm3(r.log(lo, hi));
      if (!fd) N = Math.min(N, 0.2 * (type === 'n' ? Nc(n) : Nv(n)));
      return N;
    };

    // --- layers: the semiconductor part of every topology
    let x0 = 0;
    const layers = [];
    const genLayer = (k, mat, type, L, opts = {}) => {
      const N = type === 'i' ? 0 : opts.N ?? dope(mat, type);
      const spec = { name: opts.name ?? `L${k}`, mat, type, N, L, x0, alloy: null, grade: null };
      const more = { name: spec.name };
      // doping, perhaps compensated, perhaps a profile (exponential across the layer)
      if (type !== 'i') {
        const key = type === 'n' ? 'donors' : 'acceptors', other = type === 'n' ? 'acceptors' : 'donors';
        if (!opts.plain && r.chance(0.15)) {
          on('doping-profile');
          const s = r.u(-2, 2) * (r.chance(0.5) ? 1 : 0.5), xs = Array.from({ length: 9 }, (_, j) => x0 + (L * j) / 8);
          more[key] = { x: xs, values: xs.map((x) => N * Math.pow(10, (s * (x - x0)) / L - s / 2)) };
        } else more[key] = N;
        if (!opts.plain && r.chance(0.1)) (more[other] = N * r.u(0.05, 0.5)), on('compensated');
        spec.net = (type === 'n' ? 1 : -1) * (N - (more[other] ?? 0)); // donors − acceptors, at the layer's middle
      } else spec.net = 0;
      // an alloy to a partner in the set, or a hand-graded gap, or neither
      const partners = ALLOYS.filter(([a, b]) => a === mat && M[b]);
      if (!opts.plain && partners.length && r.chance(0.2)) {
        on('alloy');
        const [, to] = r.pick(partners), ends = () => r.pick([0, 1, r.u(0, 1)]);
        let fa = ends(), fb = ends();
        if (fa === fb) fb = 1 - fa;
        const prof = r.chance(0.25);
        const step = { species: 'e-', value: eV(M[mat].chi - M[to].chi) }; // E_c(to) − E_c(from), Anderson
        const bow = to !== 'Si' && mat !== 'Si' && r.chance(0.3) ? r.u(0, 0.4) : 0;
        spec.alloy = { to, fa, fb, step: step.value, bow };
        const o = { ...more, step, fraction: prof ? { x: [0, L * r.u(0.3, 0.7), L], values: [fa, r.u(0, 1), fb] } : [fa, fb] };
        if (bow) o.bow = { 'h+': { mu0: eV(bow) } }, on('bow');
        spec.item = alloy(mat, to, L, o);
      } else if (!opts.plain && r.chance(0.12)) {
        on('graded');
        const dEc = r.u(-0.15, 0.15), dEg = r.u(-0.15, 0.25), m = M[mat].m;
        spec.grade = { dEc, dEg };
        const lin = (a, b) => ({ x: [x0, x0 + L], values: [a, b] });
        const g = {
          'e-': { mu0: lin(m.species['e-'].mu0, m.species['e-'].mu0 + eV(dEc)) },
          'h+': { mu0: lin(m.species['h+'].mu0, m.species['h+'].mu0 + eV(dEg - dEc)) },
        };
        if (r.chance(0.3)) g['e-'].D = lin(m.species['e-'].D, m.species['e-'].D * r.u(0.4, 2.5));
        if (r.chance(0.3) && m.epsr > 0) g.epsr = lin(m.epsr, m.epsr * r.u(0.8, 1.2));
        spec.item = layer(mat, L, { ...more, grade: g });
      } else spec.item = layer(mat, L, more);
      layers.push(spec);
      x0 += L;
      return spec;
    };
    // A layer's edge: its composition there, as the band edges (eV, from vacuum) and μ°(e⁻)'s shift
    // from its material's (J/mol): what a face's alignment and a barrier's zeroCharge are read against.
    const edge = (s, side) => {
      const m = M[s.mat];
      if (s.alloy) {
        const f = side === 'L' ? s.alloy.fa : s.alloy.fb, t = M[s.alloy.to];
        const Ec = -m.chi + (f * s.alloy.step) / FARADAY, Eg = (1 - f) * m.Eg + f * t.Eg - s.alloy.bow * f * (1 - f);
        return { Ec, Eg, dmu: f * s.alloy.step, plain: f === 0 };
      }
      if (s.grade && side === 'R') return { Ec: -m.chi + s.grade.dEc, Eg: m.Eg + s.grade.dEg, dmu: eV(s.grade.dEc), plain: false };
      return { Ec: -m.chi, Eg: m.Eg, dmu: 0, plain: true };
    };
    const lib = () => combine(...set.map((n) => M[n].lib));
    let faceGen = 0; // generation at a gate's interface (mol/(m²·s) at n_i), for plan.slow
    let faceG = Infinity; // the smallest conductance link at a face (S/m²), for plan.slow
    let barrier = 0; // the largest band offset a carrier crosses (eV), for plan.slow
    // A face between two semiconductor layers: aligned explicitly where they differ.
    const face = (a, b) => {
      const L = edge(a, 'R'), R = edge(b, 'L'), f = {};
      const differ = a.mat !== b.mat || !L.plain || !R.plain || a.alloy || b.alloy || a.grade;
      let dEc = R.Ec - L.Ec;
      if (differ) {
        on('hetero');
        if (r.chance(0.4)) dEc += r.u(-0.05, 0.05); // (not Anderson's)
        const how = L.plain && R.plain && a.mat !== b.mat && r.chance(0.3) ? 'vacuum' : r.chance(0.3) ? 'dipole' : 'step';
        if (how === 'vacuum') {
          f.dipole = vacuumDipole(lib(), { material: a.mat, anchor: 'e-', offset: M[a.mat].chi }, { material: b.mat, anchor: 'e-', offset: M[b.mat].chi });
          dEc = R.Ec - L.Ec;
        } else if (how === 'dipole') f.dipole = (R.dmu - L.dmu - eV(dEc)) / FARADAY; // e⁻'s standard level steps by Δμ° − F·dipole
        else f.step = { species: 'e-', value: eV(dEc) };
        const dEv = dEc - (R.Eg - L.Eg);
        barrier = Math.max(barrier, Math.abs(dEc), Math.abs(dEv));
        if (r.chance(0.07)) (f.phi = { type: 'capacitive', C: r.log(-1.3, 0.5) }), on('face-capacitive');
      }
      if (r.chance(0.08)) (f.sheetCharge = r.pick([-1, 1]) * r.log(-4, -2.7)), on('sheet-charge');
      if (r.chance(differ ? 0.15 : 0.05)) {
        on('face-srh');
        f.species = { 'e-': 'equilibrium', 'h+': 'equilibrium' };
        f.reactions = [{ equation: r.pick(['e-(left) + h+(right) = 0', 'e-(right) + h+(left) = 0', 'e-(left) + h+(left) = 0']), srh: { vn: r.log(-1, 4), vp: r.log(-1, 4) } }];
      }
      if (r.chance(0.05)) {
        on('face-conductance');
        f.species = { ...(f.species ?? { 'h+': 'equilibrium' }), 'e-': { type: 'conductance', G: r.log(4, 8) } };
        faceG = Math.min(faceG, f.species['e-'].G);
      }
      return f;
    };
    const gridFor = (s, finer = 1) => {
      const ld = debye(s.mat, s.N || ni(s.mat)), hmin = Math.max(Math.min(ld, 3e-9) / r.u(1.5, 3) / finer, 5e-11);
      return { hmin: Math.min(hmin, s.L / 20), hmax: s.L / r.u(15, 40), ratio: r.u(1.08, 1.2) };
    };
    // The depletion width that a junction or barrier of about the gap, reverse biased by `extra`,
    // pushes into a layer of doping N.
    const wdep = (n, N, extra = 1.5) => Math.sqrt((2 * epsOf(n) * (M[n].Eg + extra)) / (FARADAY * Math.max(N, cm3(1e13))));

    // --- recombination: one set of kinetics for every material
    const bulkReactions = [], keyAll = (v) => Object.fromEntries(set.map((n) => [n, v]));
    let tauMax = 0, srhO = null;
    const srh = r.chance(0.5), band = r.chance(0.3), auger = r.chance(0.2);
    if (srh) {
      on('srh');
      const tauN = r.log(-9, -5), tauP = tauN * r.log(-0.5, 0.5), o = { tauN, tauP };
      if (r.chance(0.3)) (o.n1 = ni(set[0]) * Math.exp(r.u(-8, 8))), on('srh-n1');
      bulkReactions.push({ equation: 'e- + h+ = 0', srh: keyAll(o) });
      tauMax = Math.max(tauMax, tauN, tauP), (srhO = o);
    }
    let kfBand = 0;
    if (band) {
      on('band');
      const tau = r.log(-9, -6), majority = cm3(r.log(15, 18));
      const rx = recombination({ material: set[0], tau, majority });
      kfBand = rx.kf[set[0]];
      bulkReactions.push({ ...rx, kf: keyAll(kfBand) });
      tauMax = Math.max(tauMax, tau);
    }
    if (auger) {
      on('auger');
      bulkReactions.push({ equation: 'e- + h+ = 0', auger: keyAll({ Cn: units.cm6PerS(r.log(-31, -29.5)), Cp: units.cm6PerS(r.log(-31.5, -30)) }) });
    }
    // Traps: the same level in every material (an alloy's two ends must be alike), one or two of them.
    const traps = [];
    if (r.chance(0.25)) {
      on('traps');
      const count = r.chance(0.2) ? 2 : 1;
      for (let j = 0; j < count; j++) {
        const o = { density: cm3(r.log(12, 15.5)), tauN: r.log(-9, -6), tauP: r.log(-9, -6), charge: r.pick(['acceptor', 'donor']), name: j ? 'Y' : 'X' };
        if (r.chance(0.5)) o.n1 = ni(set[0]) * Math.exp(r.u(-6, 6));
        if (o.charge === 'donor') on('donor-trap');
        traps.push(o);
      }
    }
    let library = combine(...set.map((n) => M[n].lib)); // (the materials as now: statistics, ε)
    const trapMade = [];
    for (const o of traps) {
      const per = {};
      for (const n of set) {
        const t = trap(library, { ...o, material: n });
        library = t.library;
        per[n] = t;
      }
      // (one reaction per equation, keyed to every material: an alloy's two ends alike)
      const rx = per[set[0]].bulkReactions.map((b, k) => ({ equation: b.equation, kf: Object.fromEntries(set.map((n) => [n, per[n].bulkReactions[k].kf[n]])) }));
      bulkReactions.push(...rx);
      trapMade.push({ o, per });
    }
    const trapC0 = (s) => {
      const c0 = {}, mid = s.net >= 0 ? { donors: s.net } : { acceptors: -s.net };
      for (const { per } of trapMade) Object.assign(c0, per[s.mat].c0For(mid));
      return c0;
    };
    const light = r.chance(0.3);
    if (light) {
      on('light');
      if (!srh && !band && !auger && !traps.length) {
        bulkReactions.push({ equation: 'e- + h+ = 0', srh: keyAll({ tauN: 1e-7, tauP: 1e-7 }) });
        tauMax = Math.max(tauMax, 1e-7), (srhO = { tauN: 1e-7, tauP: 1e-7 });
      }
    }
    const photonMu = eV(Math.max(...set.map((n) => M[n].Eg)) + 0.6);

    // --- the topology
    let stack, ports = [], geometry, plan = {}, Cmax = 0, extraSlow = 0, unique = true, equilibrium = !light;
    const junctions = []; // each one's forward sign for the right terminal's V
    let Rspec = 0, gateC = 0;
    // what passes each carrier at the semiconductor chain's two ends, its ports, its gated layer
    let chain = { left: [], right: [], ports: [], gated: null };
    const typeSign = (t) => (t === 'n' ? 1 : t === 'p' ? -1 : 0);

    if (topology === 'stack') {
      // a bipolar transistor's three alternating layers now and then, else 1–4 of any type
      const bjt = r.chance(0.12) ? r.pick(['n', 'p']) : null;
      const n = bjt ? 3 : r.pick([1, 2, 2, 3, 3, 4]);
      const kinds = [0, 1].map(() => r.pick(['ohmic', 'ohmic', 'ohmic', 'velocity', 'selective', 'schottky', 'schottky', 'metal']));
      // (a metal region on the left moves every layer along: profiles are against the device's x)
      const metalL = kinds.map((k) => (k === 'metal' ? r.log(-8, -6.5) : 0));
      x0 = metalL[0];
      let prev = null;
      for (let k = 0; k < n; k++) {
        const mat = prev && set.length > 1 && r.chance(0.6) ? r.pick(set.filter((m) => m !== prev.mat)) : r.pick(set);
        const type = bjt ? (k === 1 ? (bjt === 'n' ? 'p' : 'n') : bjt) : r.pick(['n', 'p', 'n', 'p', 'i']);
        let L = r.log(-7.3, -5.5);
        const s = genLayer(k, mat, type, L);
        // An end layer holds the depletion a reverse bias pushes into it, mostly (cut off at a
        // 'bulk' contact, it's warned of and is what a page shouldn't do).
        if ((k === 0 || k === n - 1) && type !== 'i' && r.chance(0.85)) {
          const need = Math.min(1.3 * wdep(mat, s.N), 8e-6);
          if (need > L) (x0 += need - L), (s.L = L = need), (s.item = { ...s.item, length: need }), fixProfiles(s);
        }
        prev = s;
      }
      // Junctions, each with the sign of V at the right that forward-biases it (n|p: +1), across
      // any intrinsic layers between (p-i-n is one)
      const doped = layers.filter((s) => s.type !== 'i').map((s) => typeSign(s.type));
      for (let k = 0; k + 1 < doped.length; k++) if (doped[k] !== doped[k + 1]) junctions.push(doped[k] > 0 ? 1 : -1);

      const ends = ['left', 'right'].map((side, j) => ({ side, s: j ? layers[n - 1] : layers[0], kind: kinds[j], metalL: metalL[j] }));
      // Every carrier passes somewhere (a carrier blocked at both ends is a closed amount).
      // (holes meeting a metal's electrons by a face reaction don't reach a terminal: still unlinked)
      const passes = (e) => (e.kind === 'selective' ? [majority(e.s)] : e.kind === 'metal' ? ['e-'] : ['e-', 'h+']);
      for (const e of ends) if (e.kind === 'metal') e.hReact = r.chance(0.5);
      if (!['e-', 'h+'].every((c) => ends.some((e) => passes(e).includes(c)))) (ends[1].kind = 'ohmic'), (ends[1].metalL = 0);
      chain.left = passes(ends[0]);
      chain.right = passes(ends[1]);
      chain.barriers = [];
      for (const e of ends) {
        if (e.kind === 'schottky') chain[`${e.side}Barrier`] = true;
        if (e.kind === 'schottky' || e.kind === 'metal') chain.barriers.push({ name: e.s.name, C: epsOf(e.s.mat) / Math.min(e.s.L, wdep(e.s.mat, Math.max(e.s.N, ni(e.s.mat)), 0)) });
      }
      for (const e of ends) {
        if (e.kind === 'schottky' || e.kind === 'metal') {
          const t = typeSign(e.s.type);
          if (t) junctions.push(e.side === 'right' ? t : -t);
        }
      }
      // Drives: the right held, or behind R, or driven by a current (then the left is the one moved).
      const Idriven = r.chance(0.12);
      const R = !Idriven && r.chance(0.2) ? r.log(-6, -3) : 0;
      const Rleft = r.chance(0.08) ? r.log(-6, -3.5) : 0;
      const signs = new Set(junctions);
      const fwd = signs.size === 1 ? [...signs][0] : 0;
      let Iright = null;
      if (Idriven) {
        on('current-driven');
        unique = false;
        // A forward current needs a known forward direction, and somewhere for the carriers to meet
        const canFlow = srh || band || auger || traps.length || ends.some((e) => ['ohmic', 'velocity', 'schottky'].includes(e.kind) || e.hReact);
        // (and no intrinsic layer: an insulator passes only a space-charge-limited current, which a
        // 'bulk' contact, holding its end neutral, can't inject; docs/device.md, pinned contacts)
        const insulating = layers.some((s) => s.type === 'i');
        if (light || insulating || (junctions.length && (!fwd || !canFlow)) || r.chance(0.4)) Iright = 0;
        else Iright = (junctions.length ? fwd : r.pick([-1, 1])) * r.log(-1, 2);
        if (Iright !== 0) equilibrium = false;
      }
      // geometry, now and then (its areas turn R into Ω and I into A)
      let Aleft = 1, Aright = 1;
      const Ltot = x0 + ends[1].metalL;
      if (r.chance(0.08)) {
        on('geometry');
        const kind = r.pick(['cylindrical', 'spherical', 'profile']), r0 = r.log(-6.5, -5);
        if (kind === 'profile') {
          const A0 = r.log(-11, -8), A1 = A0 * r.log(-0.5, 0.5);
          geometry = { area: { x: [0, Ltot], values: [A0, A1] } };
          (Aleft = A0), (Aright = A1);
        } else {
          geometry = { type: kind, r0 };
          const A = (rr) => (kind === 'spherical' ? 4 * Math.PI * rr * rr : 2 * Math.PI * rr);
          (Aleft = A(r0)), (Aright = A(r0 + Ltot));
        }
      }
      const drives = {
        left: Rleft ? { V: 0, R: Rleft / Aleft } : 0,
        right: Iright !== null ? { I: Iright * Aright } : R ? { V: 0, R: R / Aright } : 0,
      };
      if (R || Rleft) on('series-R'), (Rspec = Math.max(R, Rleft));
      stack = [contact(ends[0], drives.left)];
      layers.forEach((s, k) => {
        if (k > 0) {
          const f = face(layers[k - 1], s);
          if (Object.keys(f).length) stack.push(f);
        }
        stack.push(s.item);
      });
      // metal-region contacts sit between the contact and the end layer
      if (ends[0].kind === 'metal') stack.splice(1, 0, ...metalSide(ends[0]));
      if (ends[1].kind === 'metal') stack.splice(stack.length, 0, ...metalSide(ends[1]).reverse());
      stack.push(contact(ends[1], drives.right));
      for (const e of ends) on(`contact-${e.kind}`);

      // Ports: a base contact on the middle of three alternating layers (a bipolar transistor), an
      // injector of the minority carrier.
      let terminal = Iright !== null ? 'left' : 'right';
      if (n === 3 && layers[1].type !== 'i' && layers[0].type === layers[2].type && layers[0].type !== layers[1].type && r.chance(bjt ? 0.8 : 0.5)) {
        on('base-port');
        const b = layers[1], c = majority(b);
        const from = b.L * r.u(0.25, 0.45), to = b.L * r.u(0.55, 0.75);
        const link = r.chance(0.7) ? 'equilibrium' : { type: 'conductance', G: r.log(5, 8) / (to - from) }; // (S/m³: 1e5–1e8 S/m² over the window)
        ports.push({ name: 'base', region: b.name, from, to, V: 0, terminal: c, species: { [c]: link } });
        chain.ports.push({ layer: b.name, c });
        if (Iright === null && r.chance(0.5)) terminal = 'base';
      }
      if (!geometry && r.chance(0.07)) {
        on('injector');
        const s = r.pick(layers), sp = s.type === 'n' ? 'h+' : s.type === 'p' ? 'e-' : r.pick(['e-', 'h+']);
        const J = r.log(-1, 1.5) * (sp === 'h+' ? 1 : -1);
        ports.push(injector({ region: s.name, from: s.L * 0.4, to: s.L * 0.6, species: sp, I: J }));
        unique = false, (equilibrium = false);
      }
      // pnpn: an S-shaped I–V, more than one steady state at one voltage
      if (junctions.length >= 3) unique = false;
      const EgMin = Math.min(...layers.map((s) => Math.min(edge(s, 'L').Eg, edge(s, 'R').Eg)));
      let V;
      if (terminal === 'base') V = r.u(-0.6, 0.6) * EgMin;
      else if (terminal === 'left') V = r.u(-0.4, 0.4);
      else if (!junctions.length) V = r.u(-0.15, 0.15);
      else if (!fwd) V = r.u(-1.5, 1.5);
      else V = fwd * r.u(-1.5, 0.75 * EgMin);
      if (R) V *= 1.3;
      plan = { terminal, V };
      Cmax = Math.max(...layers.map((s) => epsOf(s.mat) / debye(s.mat, s.N || ni(s.mat))));
      if (light) bulkReactions.push(...lightOn(layers, r.chance(0.7) ? '+x' : '-x'));
    } else if (topology === 'mos') {
      const n = r.chance(0.65) ? 1 : 2;
      const bodyType = r.pick(['n', 'p', 'p']);
      for (let k = 0; k < n; k++) {
        const mat = k > 0 && set.length > 1 && r.chance(0.6) ? r.pick(set.filter((m) => m !== layers[0].mat)) : r.pick(set);
        const type = k === 0 ? bodyType : r.pick([bodyType, bodyType, 'i']);
        const N0 = type === 'i' ? 0 : dope(mat, type, 14.5, fd ? 19 : 18);
        const L = k === 0 ? Math.min(Math.max(1.5 * wdep(mat, N0, 1), r.log(-6.5, -5.5)), 20e-6) : r.log(-8, -7);
        genLayer(k, mat, type, L, { N: N0, name: k === 0 ? 'body' : 'channel' });
      }
      const top = layers[n - 1], body = layers[0];
      let kind = r.pick(['contact', 'oxide', 'oxide', 'floating']);
      const oxide = r.pick(Object.keys(OXIDES)), tox = r.log(-8.5, -7.6), Cox = (OXIDES[oxide] * EPS0) / tox;
      const tox2 = r.log(-8.5, -7.6), Cox2 = (OXIDES[oxide] * EPS0) / tox2;
      const gm = r.pick(Object.keys(WORK)), W = WORK[gm] + r.u(-0.2, 0.2);
      const e = edge(top, 'R');
      // zeroCharge at a conductor across the oxide: W + (V_vac − φ) on the semiconductor's side
      // (the oxide's φ continuous with it), the edge's own μ°(e⁻) shift included.
      const zcAt = (Wm) => vacuumZeroCharge(lib(), Wm, { material: top.mat, anchor: 'e-', offset: M[top.mat].chi }) - e.dmu / FARADAY;
      const Vg0 = r.chance(0.3) ? r.u(-1, 1) : 0;
      const Rg = r.chance(0.15) ? r.log(-4, 0) : 0;
      if (Rg) on('series-R'), (Rspec = Rg);
      const gateDrive = Rg ? { V: Vg0, R: Rg } : Vg0;
      stack = [ohmic(0), ...layers.flatMap((s, k) => (k ? [face(layers[k - 1], s), s.item] : [s.item]))];
      const ox = { dipole: 0 };
      if (r.chance(0.3)) (ox.sheetCharge = r.pick([-1, 1]) * r.log(-4.8, -2.8)), on('oxide-charge');
      if (r.chance(0.25)) (ox.reactions = [{ equation: 'e- + h+ = 0', srh: { vn: r.log(-1, 3), vp: r.log(-1, 3) } }]), on('interface-srh');
      if (kind === 'contact' && (ox.sheetCharge || ox.reactions)) kind = 'oxide'; // (those need the oxide's face)
      on(`gate-${kind}`);
      const metalFace = (Wm) => ({ phi: { type: 'capacitive', C: r.u(2, 10) }, zeroCharge: zcAt(Wm) });
      const oxGrid = (t) => ({ hmin: t / 10, hmax: t / 4 });
      if (kind === 'contact') {
        stack.push({ ...(typeof gateDrive === 'number' ? { V: gateDrive } : gateDrive), phi: { type: 'capacitive', C: Cox }, zeroCharge: zcAt(W) });
        gateC = Cox;
      } else if (kind === 'oxide') {
        stack.push(ox, layer(oxide, tox, { grid: oxGrid(tox) }), metalFace(W), layer(gm, r.log(-8, -7)), ohmic(gateDrive, ['e-']));
        gateC = Cox;
      } else {
        // a floating gate between two oxides; the second oxide's φ anchored as the first's
        const fm = r.pick(Object.keys(WORK)), Wf = WORK[fm] + r.u(-0.2, 0.2);
        const zc2 = (Wm) => zcAt(Wm); // (the same constant: a speciesless oxide's φ is bookkeeping)
        stack.push(ox, layer(oxide, tox, { grid: oxGrid(tox) }), metalFace(Wf), layer(fm, r.log(-8, -7.3), { name: 'fg' }),
          { phi: { type: 'capacitive', C: r.u(2, 10) }, zeroCharge: zc2(Wf) }, layer(oxide, tox2, { grid: oxGrid(tox2) }), metalFace(W), layer(gm, r.log(-8, -7)), ohmic(gateDrive, ['e-']));
        gateC = 1 / (1 / Cox + 1 / Cox2);
        const wire = r.pick(['none', 'held', 'held', 'open']);
        if (wire !== 'none') {
          on(`fg-wire-${wire}`);
          const G = r.log(2, 6), Rw = r.chance(0.3) ? r.log(-4, -1) : 0;
          ports.push({ name: 'fgwire', region: 'fg', terminal: 'e-', species: { 'e-': { type: 'conductance', G } }, ...(wire === 'open' ? { I: 0 } : Rw ? { V: 0, R: Rw } : { V: 0 }) });
          extraSlow = Math.max(extraSlow, (1 / G + Rw) * (Cox + Cox2));
          if (wire === 'open') unique = false;
        } else unique = false; // (its charge is the start's: a closed amount)
      }
      // a port holding the minority carrier beside the oxide: source and drain's stand-in
      const minority = body.type === 'n' ? 'h+' : 'e-';
      let supplied = false;
      if (r.chance(0.4)) {
        on('minority-port');
        const w = Math.min(r.log(-8.3, -7.5), top.L / 3), link = r.chance(0.7) ? 'equilibrium' : { type: 'conductance', G: r.log(2, 6) / w };
        ports.push({ name: 'sd', region: top.name, from: top.L - w, to: top.L, V: 0, terminal: minority, species: { [minority]: link } });
        supplied = true;
        if (link !== 'equilibrium') extraSlow = Math.max(extraSlow, gateC / (link.G * w));
      }
      chain = { left: ['e-', 'h+'], right: [], ports: supplied ? [{ layer: top.name, c: minority }] : [], gated: { name: top.name, C: gateC } };
      if (ox.reactions) faceGen = ni(top.mat) * Math.min(ox.reactions[0].srh.vn, ox.reactions[0].srh.vp);
      Cmax = gateC;
      // Biased about the flat band, V_FB = zeroCharge + φ_bulk (the oxide charge aside): within
      // 0.4 V/nm of oxide, and under Boltzmann statistics no further past threshold than a surface
      // layer a fifth of the band's density of states (beyond, a page would turn Fermi–Dirac on).
      const VFB = zcAt(W) + flatPhi(top);
      const Qmax = 0.3 * FARADAY * Math.min(Nc(top.mat), Nv(top.mat)) * 5e-9;
      const span = Math.min(3, 0.4 * tox * 1e9, fd ? Infinity : M[top.mat].Eg + Qmax / gateC);
      plan = { terminal: 'right', V: VFB + r.u(-span, span) };
      if (light) bulkReactions.push(...lightOn(layers, '+x'));
    } else {
      // a thin-film transistor along x: the film between source and drain, the gate a port's
      // capacitance through it
      const mat = set[0], t = r.log(-8.3, -7.3), Lch = r.log(-6, -4.7);
      const type = r.pick(['n', 'p', 'i', 'n']);
      const sd = r.chance(0.4) ? r.pick(['n', 'p']) : null;
      if (sd) {
        on('source-drain');
        const Nsd = dope(mat, sd, 18, fd ? 19.5 : 18.7), Lsd = r.log(-6.5, -5.7);
        genLayer(0, mat, sd, Lsd, { N: Nsd, name: 'source', plain: true });
        genLayer(1, mat, type, Lch, { name: 'channel', N: type === 'i' ? 0 : dope(mat, type, 14, 17), plain: true });
        genLayer(2, mat, sd, Lsd, { N: Nsd, name: 'drain', plain: true });
        if (sd !== type && type !== 'i') junctions.push(1, -1);
      } else genLayer(0, mat, type, Lch, { name: 'channel', N: type === 'i' ? 0 : dope(mat, type, 14, 17), plain: true });
      const ch = layers.find((s) => s.name === 'channel');
      const Vg0 = r.chance(0.3) ? r.u(-1, 1.5) : 0, Cg = (OXIDES[r.pick(Object.keys(OXIDES))] * EPS0) / r.log(-8.3, -7);
      const inset = r.chance(0.3) ? Lch * r.u(0.02, 0.15) : 0;
      if (inset) on('gate-underlap');
      ports.push({ name: 'gate', region: 'channel', ...(inset ? { from: inset, to: Lch - inset } : {}), V: Vg0, area: 1 / t, capacitance: { C: Cg, zeroCharge: r.u(-0.8, 0.8) } });
      gateC = Cg;
      if (r.chance(0.15)) {
        on('dual-gate');
        const C2 = Cg * r.log(-0.5, 0.5);
        ports.push({ name: 'gate2', region: 'channel', V: r.u(-0.5, 0.5), area: 1 / t, capacitance: { C: C2, zeroCharge: r.u(-0.8, 0.8) } });
        gateC += C2;
      }
      const Rd = r.chance(0.15) ? r.log(-6, -3) : 0;
      if (Rd) on('series-R'), (Rspec = Rd);
      ch.item.grid = { hmin: Math.min(debye(mat, ch.N || ni(mat)), Lch / 400), hmax: Lch / r.u(40, 100), ratio: r.u(1.05, 1.15) };
      stack = [ohmic(0), ...layers.map((s) => s.item), ohmic(Rd ? { V: 0, R: Rd } : 0)];
      // The slowest charging: the gate through the film at its least conductive, both carriers at n_i.
      const sigmaSheet = (FARADAY * FARADAY * ni(mat) * (M[mat].m.species['e-'].D + M[mat].m.species['h+'].D) * t) / RT;
      extraSlow = Math.max(extraSlow, (Lch * Lch * gateC) / sigmaSheet);
      chain = { left: ['e-', 'h+'], right: ['e-', 'h+'], ports: [], gated: { name: 'channel', C: (gateC * Lch) / t }, hmin: ch.item.grid.hmin };
      Cmax = gateC;
      const terminal = r.chance(0.5) ? 'gate' : 'right';
      // the gate about its flat band, Boltzmann's span held as the MOS's is (a sheet in the film)
      const VFB = ports[0].capacitance.zeroCharge + flatPhi(ch);
      const Qmax = 0.3 * FARADAY * Math.min(Nc(mat), Nv(mat)) * t;
      const span = fd ? 2.5 : Math.min(2.5, M[mat].Eg + Qmax / Cg);
      plan = { terminal, V: terminal === 'gate' ? VFB + r.u(-span, span) : r.u(-1, 2) };
      if (light) bulkReactions.push({ equation: 'photon = e- + h+', fixed: { photon: photonMu }, kf: keyAll(r.log(-5, -2) / t) });
    }
    // trap totals (every layer holds its material's traps) and grids: on the layers' items, which
    // the stack holds
    for (const s of layers) {
      if (trapMade.length) s.item.c0 = trapC0(s);
      if (!s.item.grid) s.item.grid = gridFor(s, topology === 'mos' && s === layers.at(-1) ? 3 : 1);
    }

    const def = build({ T, library: [library, ...metals(stack)], materials: oxides(stack), stack, bulkReactions, ...(ports.length ? { ports } : {}), ...(geometry ? { geometry } : {}) });
    // c0For() splits a trap's total by the layer's net doping alone; where the bands or the doping
    // vary across a layer, or the traps' own charge matters to its neutrality, start its traps at
    // the split they hold at local, neutral equilibrium instead, their charge counted (from a split
    // a few times off at one end of a graded layer, a cold solve fails: see the report's suspects).
    if (trapMade.length) {
      let xr = 0;
      def.regions.forEach((reg) => {
        const x0r = xr;
        xr += reg.length;
        const mat = def.materials[reg.material];
        if (!mat?.species?.['e-']) return;
        const netAt = (x) => (reg.fixedCharge?.x ? profileAt(reg.fixedCharge, x) : reg.fixedCharge ?? 0) / FARADAY;
        const uniform = !reg.grade && !reg.fixedCharge?.x;
        if (uniform && traps.reduce((t, o) => t + o.density, 0) < 0.05 * Math.abs(netAt(x0r))) return; // (c0For's)
        const xs = uniform ? [x0r] : Array.from({ length: 65 }, (_, j) => x0r + (reg.length * j) / 64), c0 = {};
        const names = traps.map((o) => (o.charge === 'acceptor' ? [`${o.name}0`, `${o.name}-`] : [`${o.name}+`, `${o.name}0`]));
        names.forEach(([e, f]) => ((c0[e] = []), (c0[f] = [])));
        for (const x of xs) {
          const at = (sp, f) => (reg.grade?.[sp]?.[f] ? profileAt(reg.grade[sp][f], x) : mat.species[sp][f]);
          const ni2 = at('e-', 'cRef') * at('h+', 'cRef') * Math.exp(-(at('e-', 'mu0') + at('h+', 'mu0')) / RT), net = netAt(x);
          // each trap's filled/empty ratio per electron: from e⁻ + empty = full at equilibrium
          const per = traps.map((o, j) => (at(names[j][1], 'cRef') / at(names[j][0], 'cRef')) * Math.exp((at('e-', 'mu0') + at(names[j][0], 'mu0') - at(names[j][1], 'mu0')) / RT) / at('e-', 'cRef'));
          const split = (n) => traps.map((o, j) => [o.density / (1 + per[j] * n), o.density / (1 + 1 / (per[j] * n))]);
          // neutrality, net + p − n + the traps' charge = 0, by bisection in ln n (it falls as n rises)
          const charge = (n) => net + ni2 / n - n + split(n).reduce((t, [e, f], j) => t + (traps[j].charge === 'acceptor' ? -f : e), 0);
          let lo = Math.log(ni2) / 2 - 120, hi = Math.log(Math.max(Math.abs(net), Math.sqrt(ni2))) + 10;
          for (let it = 0; it < 200; it++) {
            const mid = (lo + hi) / 2;
            if (charge(Math.exp(mid)) > 0) lo = mid;
            else hi = mid;
          }
          split(Math.exp((lo + hi) / 2)).forEach(([e, f], j) => (c0[names[j][0]].push(e), c0[names[j][1]].push(f)));
        }
        for (const k of Object.keys(c0)) c0[k] = uniform ? c0[k][0] : { x: xs, values: c0[k] };
        reg.c0 = { ...reg.c0, ...c0 };
      });
    }

    // --- plan.slow: the largest of the device's own times, each an over- rather than an underestimate
    const Lsemi = layers.reduce((t, s) => t + s.L, 0);
    const times = { transport: ((Lsemi * Lsemi) / Dmin) * Math.exp(Math.min(barrier / VT, 15)) };
    for (const o of traps) {
      // a trap's slowest exchange, where both carriers are scarce: emission alone
      for (const n of set) {
        const n1 = o.n1 ?? ni(n), p1 = ni(n) ** 2 / n1;
        times[`trap ${o.name} ${n}`] = 1 / (n1 / (o.tauN * o.density) + p1 / (o.tauP * o.density));
      }
    }
    // a series R across the device's largest capacitance: depletion, gate, or a forward diode's
    // stored charge (1 kA/m² for its lifetime)
    if (Rspec) times.RC = Rspec * (Cmax + (1e3 * Math.min(tauMax || times.transport, times.transport)) / VT);
    if (faceG < Infinity) times.faceRC = Cmax / faceG;
    // a junction's (or a barrier's, or a gate's) charge, through the neutral layers' resistance
    const doped = layers.filter((s) => s.type !== 'i');
    if (doped.length && (junctions.length || gateC || stack.some((it) => it.phi?.type === 'capacitive' || it.phi === 'pinned'))) {
      const Rn = doped.reduce((t, s) => t + Math.max(s.L - wdep(s.mat, s.N, 0), 0.1 * s.L) / ((FARADAY * FARADAY * s.N * M[s.mat].m.species[majority(s)].D) / RT), 0);
      const Cj = Math.max(gateC, ...doped.map((s) => epsOf(s.mat) / Math.min(s.L, wdep(s.mat, s.N, 0))));
      times.junctionRC = Cj * Rn;
    }
    if (extraSlow) times.gate = extraSlow;
    Object.assign(times, isolationTimes());
    for (const s of layers) {
      const N = Math.max(s.N, ni(s.mat)), D = Math.max(M[s.mat].m.species['e-'].D, M[s.mat].m.species['h+'].D);
      if (s.L > 3 * debye(s.mat, N)) times[`dielectric ${s.name}`] = (epsOf(s.mat) * RT) / (FARADAY * FARADAY * N * D);
    }
    // plan.slow is ten times the largest. Checked against transients on the first 200 cases (a dV
    // step from rest, run to 30·slow, its terminal currents within 1e-3 of their first deviation
    // from the steady ones): the largest estimate alone fell short by up to 50×; ten times it, in
    // 7 of 200, of which 3 are lit traps whose totals the time steps don't keep (the report's
    // suspect) and 4 modes these estimates miss, still 7–11× short (a lit, lightly doped layer
    // behind a conductance face; a trap in a graded heterojunction; a capacitive barrier over
    // depleted layers; a TFT's underlapped gate). params.slowModes keeps each estimate, for triage.
    const slowest = Object.entries(times).sort((a, b) => b[1] - a[1])[0][0];
    const slow = 10 * times[slowest];

    // The rest is at equilibrium unless lit or driven by a current; every terminal starts at its
    // rest drive (a gate away from 0 V still passes nothing).
    return {
      def,
      params: {
        topology, T, materials: set, x: { AlGaAs: set.includes('AlGaAs') ? xAl : undefined, InGaAs: set.includes('InGaAs') ? xIn : undefined },
        layers: layers.map((s) => ({ name: s.name, mat: s.mat, type: s.type, N: s.N && units.toPerCm3(s.N), L: s.L, alloy: s.alloy?.to, grade: !!s.grade })),
        traps: traps.map((o) => ({ ...o, density: units.toPerCm3(o.density) })), junctions: junctions.length, barrier, slowest,
        slowModes: Object.fromEntries(Object.entries(times).map(([k, v]) => [k, +v.toPrecision(3)])),
        features: [...features],
      },
      plan: { terminal: plan.terminal, V: plan.V, dV: r.log(-2.3, -1.3) * (topology === 'stack' ? 1 : 2), equilibrium, unique, slow },
    };

    // --- helpers that read the case so far
    // Slow modes where a carrier reaches a terminal only through a region where it's a minority:
    // (1) a layer whose majority carrier has no path of its own to a terminal (a floating base, a
    // p layer behind a metal that blocks holes): its charge relaxes through the leak across its
    // neighbours, diffusion of the scarce carrier plus generation; (2) a gate's inversion carrier
    // fed only that way (a MOS without source and drain, a TFT whose contacts hold the channel's
    // majority). Each is the capacitance over that leak's conductance, linearised.
    function isolationTimes() {
      const out = {}, nL = layers.length, F2RT = (FARADAY * FARADAY) / RT;
      const dens = (s, c) => {
        const n2 = ni(s.mat) ** 2, net = s.net, root = Math.sqrt((net * net) / 4 + n2);
        const nE = net >= 0 ? net / 2 + root : n2 / (root - net / 2);
        return c === 'e-' ? nE : n2 / nE;
      };
      const W = (s) => Math.min(s.L, wdep(s.mat, Math.max(s.N, ni(s.mat)), 0));
      const gen = (s) => {
        const n2 = ni(s.mat) ** 2, n1 = srhO?.n1 ?? ni(s.mat);
        return (srhO ? (n2 * W(s)) / (srhO.tauP * n1 + srhO.tauN * (n2 / n1)) : 0) + (kfBand ? kfBand * n2 * W(s) : 0);
      };
      // conductance (S/m²) of carrier c through layer j, where it's scarce
      const leak = (j, c) => F2RT * ((M[layers[j].mat].m.species[c].D * dens(layers[j], c)) / layers[j].L + gen(layers[j]));
      const runOf = (k, ok) => {
        let a = k, b = k;
        while (a > 0 && ok(layers[a - 1])) a--;
        while (b < nL - 1 && ok(layers[b + 1])) b++;
        return [a, b];
      };
      // a boundary's slow leak for a carrier (holes into a metal by a face reaction)
      const edgeLeak = ([a, b], c) => (a === 0 ? chain.leftLeak?.[c] ?? 0 : 0) + (b === nL - 1 ? chain.rightLeak?.[c] ?? 0 : 0);
      const reaches = ([a, b], c) => (a === 0 && chain.left.includes(c)) || (b === nL - 1 && chain.right.includes(c)) || chain.ports.some((p) => p.c === c && layers.slice(a, b + 1).some((s) => s.name === p.layer));
      const gates = chain.gated ? [chain.gated] : chain.barriers ?? [];
      const gatedIn = ([a, b]) => gates.filter((q) => layers.slice(a, b + 1).some((s) => s.name === q.name)).reduce((t, q) => t + q.C, 0);
      layers.forEach((s, k) => {
        if (s.type === 'i') return;
        const c = majority(s), run = runOf(k, (q) => q.type === s.type || q.type === 'i');
        if (reaches(run, c)) return;
        const G = (run[0] > 0 ? leak(run[0] - 1, c) : 0) + (run[1] < nL - 1 ? leak(run[1] + 1, c) : 0) + edgeLeak(run, c);
        if (G > 0) out[`isolated ${s.name}`] = (epsOf(s.mat) / W(s) + gatedIn(run)) / G;
      });
      // (a barrier at an end, a Schottky contact or a metal, holds an inversion layer under it too)
      for (const gate of gates) {
        const g = layers.findIndex((s) => s.name === gate.name), sg = layers[g];
        for (const c of sg.type === 'i' ? ['e-', 'h+'] : [sg.type === 'n' ? 'h+' : 'e-']) {
          if (chain.ports.some((p) => p.layer === sg.name && p.c === c)) continue;
          let G = F2RT * (gen(sg) * (sg.L / W(sg)) + (chain.gated ? faceGen : 0)), fast = false;
          for (const dir of [-1, 1]) {
            let j = g + dir, len = 0, n = dens(sg, c);
            const maj = (q) => (q.type === 'n' ? 'e-' : q.type === 'p' ? 'h+' : null);
            while (j >= 0 && j < nL && maj(layers[j]) !== c) (len += layers[j].L), (n = Math.min(n, dens(layers[j], c))), (j += dir);
            const side = dir < 0 ? 'left' : 'right';
            if (j >= 0 && j < nL) {
              if (reaches(runOf(j, (q) => maj(q) === c || q.type === 'i'), c)) fast = true;
            } else if (chain[side].includes(c)) {
              // a barrier contact sets the carrier's density at its face itself: no shortage there
              if (chain[`${side}Barrier`]) fast = true;
              const Gd = (F2RT * M[sg.mat].m.species[c].D * n) / Math.max(len, chain.hmin ?? 1e-9), v = chain[`${side}Vel`]?.[c];
              G += v ? 1 / (1 / Gd + 1 / (F2RT * v * n)) : Gd; // (a surface velocity in series)
            } else G += chain[`${side}Leak`]?.[c] ?? 0;
          }
          if (!fast && G > 0) out[`inversion ${c} ${sg.name}`] = gate.C / G;
        }
      }
      return out;
    }
    function majority(s) {
      return s.type === 'p' ? 'h+' : 'e-';
    }
    // φ in a layer's neutral bulk at equilibrium with the electrons at 0 V (Boltzmann; its right edge)
    function flatPhi(s) {
      const e = edge(s, 'R'), n2 = ni(s.mat) ** 2, net = s.net, root = Math.sqrt((net * net) / 4 + n2);
      const nE = net >= 0 ? net / 2 + root : n2 / (root - net / 2);
      return e.dmu / FARADAY + VT * Math.log(nE / Nc(s.mat));
    }
    function fixProfiles(s) {
      // (a layer stretched after its profiles were laid out: lay them out again over its length)
      const scale = (p) => ({ x: p.x.map((x, j) => s.x0 + (s.L * j) / (p.x.length - 1)), values: p.values });
      for (const key of ['donors', 'acceptors']) if (s.item[key]?.x) s.item[key] = scale(s.item[key]);
      if (s.item.grade) for (const v of Object.values(s.item.grade)) for (const k of Object.keys(v ?? {})) if (v[k]?.x) v[k] = scale(v[k]);
      if (s.item.grade?.epsr?.x) s.item.grade.epsr = scale(s.item.grade.epsr);
      // (the last knot exactly the length: one an ulp short, alloy() wrote a zero-width piece; see the report)
      if (s.item.alloy?.fraction?.x) s.item.alloy.fraction = { x: s.item.alloy.fraction.x.map((x, j, a) => (j === a.length - 1 ? s.L : (s.L * x) / a.at(-1))), values: s.item.alloy.fraction.values };
    }
    // A Schottky barrier: the work function's barrier for electrons at that edge, clamped into the
    // gap (a pinned Fermi level, as real surfaces are), as zeroCharge = φ_B − Δμ°(e⁻)/F.
    function barrierZC(s, side, Wm, out = {}) {
      const e = edge(s, side), phiB = Math.min(Math.max(Wm + e.Ec, 0.2), e.Eg - 0.2);
      (out.phiB = phiB), (out.Eg = e.Eg);
      return vacuumZeroCharge(lib(), phiB + M[s.mat].chi, { material: s.mat, anchor: 'e-', offset: M[s.mat].chi }) - e.dmu / FARADAY;
    }
    function contact(e, d) {
      const drive = typeof d === 'number' ? { V: d } : d, s = e.s, maj = majority(s), min = maj === 'e-' ? 'h+' : 'e-';
      if (e.kind === 'ohmic') return ohmic(d);
      if (e.kind === 'selective') return ohmic(d, [maj]);
      if (e.kind === 'metal') return ohmic(d, ['e-']);
      if (e.kind === 'velocity') {
        const v = r.log(-1, 5);
        chain[`${e.side}Vel`] = { [min]: v };
        return { ...drive, terminal: maj, species: { [maj]: 'equilibrium', [min]: { type: 'velocity', v, offset: 0 } }, phi: 'bulk' };
      }
      // Schottky: thermionic-emission velocities v = A*T²/(F N), or equilibrium links (diffusion
      // theory); the barrier pinned, or behind an interlayer's capacitance
      const side = e.side === 'left' ? 'L' : 'R', Wm = r.u(4.2, 5.6);
      const vel = r.chance(0.6), Astar = r.log(4.5, 6.1);
      const ln = (c, N) => (vel ? { type: 'velocity', v: (Astar * T * T) / (FARADAY * N) } : 'equilibrium');
      const species = { 'e-': ln('e-', Nc(s.mat)), 'h+': vel ? { ...ln('h+', Nv(s.mat)), offset: 0 } : { type: 'equilibrium', offset: 0 } };
      const cap = r.chance(0.5);
      if (cap) on('schottky-capacitive');
      return { ...drive, terminal: 'e-', species, phi: cap ? { type: 'capacitive', C: r.log(-1, 1) } : 'pinned', zeroCharge: barrierZC(s, side, Wm) };
    }
    // A metal region and its capacitive face (left to right as for the left end; reversed for the right)
    function metalSide(e) {
      const side = e.side === 'left' ? 'L' : 'R', mname = r.pick(Object.keys(WORK)), Wm = WORK[mname] + r.u(-0.2, 0.2);
      const b = {}, f = { phi: { type: 'capacitive', C: r.log(-0.5, 1) }, zeroCharge: barrierZC(e.s, side, Wm, b) };
      if (e.hReact) {
        // holes meet the metal's electrons (minority exchange at the barrier), by Butler–Volmer:
        // a leak of exchange rate k0 (p/N_v)^½, p at the face behind the holes' barrier E_g − φ_B
        f.species = { 'e-': 'equilibrium' };
        const k0 = r.log(-6, -1);
        f.reactions = [{ equation: e.side === 'left' ? 'e-(left) + h+ = 0' : 'e-(right) + h+ = 0', k0 }];
        chain[e.side === 'left' ? 'leftLeak' : 'rightLeak'] = { 'h+': ((FARADAY * FARADAY) / RT) * k0 * Math.exp(-(b.Eg - b.phiB) / (2 * VT)) };
        on('metal-hole-exchange');
      }
      return [layer(mname, e.metalL), f];
    }
    function lightOn(ls, toward) {
      const order = toward === '+x' ? ls : [...ls].reverse(), flux = r.log(-5, -2), alpha = Object.fromEntries(set.map((n) => [n, r.log(4.5, 7)]));
      let left = flux;
      const out = [];
      for (const s of order) {
        const g = photogeneration({ material: s.mat, flux: left, alpha: alpha[s.mat], mu: photonMu, from: s.x0, to: s.x0 + s.L, toward });
        out.push({ ...g, kf: keyAll(g.kf[s.mat]) }); // (zero outside the layer, so every material may carry it)
        left *= Math.exp(-alpha[s.mat] * s.L);
      }
      return out;
    }
  },
};

// A profile { x, values } at x: piecewise linear, constant beyond its ends (as the device reads it).
function profileAt(p, x) {
  const { x: xs, values: v } = p;
  if (x <= xs[0]) return v[0];
  for (let j = 1; j < xs.length; j++) if (x <= xs[j]) return v[j - 1] + ((v[j] - v[j - 1]) * (x - xs[j - 1])) / (xs[j] - xs[j - 1]);
  return v.at(-1);
}

// The metals a stack names (as layers), as library pieces; the oxides as materials.
function metals(stack) {
  const names = new Set(stack.filter((it) => WORK[it.material] !== undefined).map((it) => it.material));
  return [...names].map((n) => metal(n));
}
function oxides(stack) {
  return Object.fromEntries(stack.filter((it) => OXIDES[it.material] !== undefined).map((it) => [it.material, { epsr: OXIDES[it.material], species: {} }]));
}
