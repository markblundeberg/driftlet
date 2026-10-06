// Stress test: random but plausible devices (semiconductor stacks with ohmic or Schottky
// contacts, MOS capacitors, electrolyte cells, electrodes with face reactions, liquid junctions), each put through what a user would do with it (a
// cold solve, a bias sweep, a transient after a step, an impedance), every result judged by
// check() and by invariants that hold whatever the device: flat levels at equilibrium, a warm
// sweep and a cold solve agreeing, a passive impedance about equilibrium. Every case comes from
// its family and index alone, so a failure reruns on its own.
//
//   npm run stress                       every family, 40 cases each
//   npm run stress -- --n 200            more cases
//   npm run stress -- semi 17            one case, with its definition and what went wrong
//   npm run stress -- semi 17 --def f    …and write its definition and plan to f
//   npm run stress -- --save             write the summary as bench/stress.json

import { readFileSync, writeFileSync } from 'node:fs';
import { Device, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, bath, aqueous, semiconductor, metal, check, photogeneration, half, vacuumZeroCharge, IONS } from '../src/kit.js';

const argv = process.argv.slice(2);
const flag = (name, dflt) => {
  const k = argv.indexOf(`--${name}`);
  return k >= 0 ? argv[k + 1] : dflt;
};
const N = +flag('n', 40), save = argv.includes('--save');
const positional = argv.filter((a, k) => !a.startsWith('--') && !(k > 0 && ['--n', '--def'].includes(argv[k - 1])));
const [onlyFamily, onlyCase] = [positional[0], positional[1] === undefined ? undefined : +positional[1]];

// --- a seeded generator per case
function rng(seed) {
  let a = seed >>> 0;
  const next = () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  return {
    u: (a0, b0) => a0 + (b0 - a0) * next(),
    log: (a0, b0) => 10 ** (a0 + (b0 - a0) * next()), // between 10^a and 10^b
    pick: (xs) => xs[Math.floor(next() * xs.length)],
    chance: (p) => next() < p,
  };
}
const hash = (s) => [...s].reduce((h, c) => Math.imul(h ^ c.charCodeAt(0), 16777619) >>> 0, 2166136261);

const T = 300, RT = GAS_CONSTANT * T, VT = RT / FARADAY;

// --- families: each returns { def, params, plan } where plan says what to do with it
const families = {
  // Semiconductor stacks of one material: 1–3 doped layers between ohmic contacts, with or without
  // recombination and light.
  semi(r) {
    const mat = r.pick(['Si', 'Ge', 'GaAs']), lib = semiconductor(mat);
    const sp = lib.materials[mat].species, Eg = (sp['h+'].mu0 + sp['e-'].mu0 * 0 - sp['e-'].mu0) / FARADAY;
    const n = r.pick([1, 2, 2, 3]), layers = [];
    let total = 0, majority = 0;
    for (let k = 0; k < n; k++) {
      const type = r.pick(['n', 'p', 'n', 'p', 'i']), N = units.perCm3(r.log(14, 19)), length = r.log(-7.3, -5);
      total += length;
      majority = Math.max(majority, type === 'i' ? 0 : N);
      layers.push(layer(mat, length, { name: `L${k}`, ...(type === 'n' ? { donors: N } : type === 'p' ? { acceptors: N } : {}) }));
    }
    const reactions = [];
    const rec = r.pick(['none', 'band', 'srh']), tau = r.log(-9, -5);
    if (rec === 'band') reactions.push({ equation: 'e- + h+ = 0', kf: { [mat]: 1 / (tau * Math.max(majority, units.perCm3(1e15))) } });
    if (rec === 'srh') reactions.push({ equation: 'e- + h+ = 0', srh: { [mat]: { tauN: tau, tauP: tau } } });
    const light = rec !== 'none' && r.chance(0.3);
    if (light) reactions.push(photogeneration({ material: mat, flux: r.log(-5, -2), alpha: r.log(4, 7), mu: FARADAY * (Eg + 0.3), to: total }));
    const def = build({
      T,
      library: [lib],
      stack: [ohmic(0), ...layers, ohmic(0)],
      bulkReactions: reactions,
      grid: { hmin: r.log(-10, -9), hmax: total / r.u(20, 100), ratio: r.u(1.08, 1.2) },
    });
    // A Schottky contact on the left now and then: both carriers at the metal's Fermi level (the
    // usual boundary condition, thermionic emission at a high velocity), the barrier from its
    // work function (Schottky–Mott).
    const schottky = r.chance(0.25);
    if (schottky) {
      const W = r.u(4.3, 5.3), chi = { Si: 4.05, Ge: 4.0, GaAs: 4.07 }[mat];
      def.contacts.left = { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: { type: 'capacitive', C: 10 }, zeroCharge: vacuumZeroCharge(def, W, { material: mat, anchor: 'e-', offset: chi }) };
    }
    return {
      def,
      params: { mat, layers: layers.map((l) => ({ length: l.length, donors: l.donors, acceptors: l.acceptors })), rec, tau, light, schottky },
      plan: { side: 'right', V: r.u(-1.5, 0.8), equilibrium: !light, step: r.u(0.05, 0.5), t: r.log(-9, -6), openCircuit: light },
    };
  },

  // MOS capacitors: a doped semiconductor, an oxide, a metal gate, with or without generation.
  mos(r) {
    const mat = r.pick(['Si', 'Ge', 'GaAs']), lib = semiconductor(mat), gateMetal = r.pick(['Au', 'Al', 'Pt']);
    const doping = units.perCm3(r.log(15, 18)), p = r.chance(0.5), L = r.log(-6.5, -5.5), tox = r.log(-9, -7.7);
    const tau = r.chance(0.5) ? r.log(-9, -4) : 0;
    const def = build({
      T,
      library: [lib, metal(gateMetal)],
      materials: { SiO2: { epsr: 3.9, species: {} } },
      stack: [
        ohmic(0),
        layer(mat, L, { name: 'bulk', ...(p ? { acceptors: doping } : { donors: doping }) }),
        { dipole: 0 },
        layer('SiO2', tox, { grid: { hmin: tox / 10, hmax: tox / 4 } }),
        { phi: { type: 'capacitive', C: 100 }, zeroCharge: r.u(-1, 1) },
        layer(gateMetal, 50e-9),
        ohmic(0, ['e-']),
      ],
      bulkReactions: tau ? [{ equation: 'e- + h+ = 0', srh: { [mat]: { tauN: tau, tauP: tau } } }] : [],
      grid: { hmin: r.log(-10, -9.3), hmax: L / r.u(20, 60), ratio: r.u(1.08, 1.15) },
    });
    return { def, params: { mat, gateMetal, doping, p, L, tox, tau }, plan: { side: 'right', V: r.u(-3, 3), equilibrium: true, step: r.u(0.01, 0.2), t: r.log(-6, 1), cap: (3.9 * EPS0) / tox } };
  },

  // An electrolyte cell: one to three salts (any of the data library's ions, balanced), strictly
  // neutral or with its double layers resolved, between baths, reversible metal electrodes or
  // blocking ones.
  cell(r) {
    const cations = ['H+', 'Li+', 'Na+', 'K+', 'Ag+', 'Cu2+', 'Zn2+', 'Fe2+'], anions = ['OH-', 'Cl-', 'NO3-', 'SO42-'];
    const cs = [...new Set([r.pick(cations), ...(r.chance(0.4) ? [r.pick(cations)] : [])])];
    const as = [...new Set([r.pick(anions), ...(r.chance(0.3) ? [r.pick(anions)] : [])])];
    const names = [...cs, ...as], z = (s) => IONS[s].z;
    const c = {};
    for (const s of cs) c[s] = r.log(-1, 3);
    const plus = cs.reduce((t, s) => t + z(s) * c[s], 0);
    const w = as.map(() => r.u(0.2, 1)), wz = as.reduce((t, s, k) => t + w[k] * -z(s), 0);
    as.forEach((s, k) => (c[s] = (plus * w[k]) / wz));
    const epsr = r.pick([0, 78.3]), L = r.log(-6, -3);
    const I = 0.5 * names.reduce((t, s) => t + z(s) * z(s) * c[s], 0), debye = Math.sqrt((epsr * EPS0 * RT) / (2 * FARADAY * FARADAY * I));
    const reversible = (s) => ({ 'Ag+': 'Ag', 'Cu2+': 'Cu', 'Zn2+': 'Zn', 'Li+': 'Li' })[s];
    const electrode = (V) => {
      const kinds = ['bath', 'bath'];
      const metalIon = cs.find(reversible);
      if (metalIon) kinds.push('reversible', 'reversible');
      if (epsr > 0) kinds.push('blocking');
      const kind = r.pick(kinds);
      if (kind === 'bath') return { kind, contact: bath({ ...c }, as[0], { V }) };
      if (kind === 'reversible') {
        return { kind, contact: { V, terminal: metalIon, species: Object.fromEntries(names.map((s) => [s, s === metalIon ? 'equilibrium' : 'blocked'])), phi: epsr > 0 ? 'bulk' : 'bulk' } };
      }
      return { kind, contact: { V, phi: { type: 'capacitive', C: r.log(-1.5, 0) }, zeroCharge: 0 } };
    };
    const left = electrode(0), right = electrode(0);
    const passes = (e) => e.kind !== 'blocking';
    if (!passes(left) && !passes(right)) left.contact = bath({ ...c }, as[0], { V: 0 }), (left.kind = 'bath');
    const def = build({
      T,
      library: [aqueous(names, { epsr })],
      stack: [left.contact, layer('water', L, { c0: { ...c } }), right.contact],
      grid: epsr > 0 ? { hmin: debye / r.u(3, 10), hmax: Math.max(L / r.u(20, 80), debye), ratio: r.u(1.08, 1.2) } : { hmin: L / r.u(200, 2000), hmax: L / r.u(20, 80), ratio: r.u(1.08, 1.2) },
    });
    const same = left.kind === right.kind && left.kind !== 'blocking';
    const blockingRight = right.kind === 'blocking';
    return {
      def,
      params: { names, c, epsr, L, debye, left: left.kind, right: right.kind },
      plan: { side: 'right', V: r.u(-0.3, 0.3), equilibrium: same, step: r.u(0.005, 0.1), t: r.log(-3, 1), noSteadyCurrent: blockingRight },
    };
  },

  // Electrodes with Butler–Volmer face reactions: a metal in its own ions (silver, copper, zinc) or
  // platinum in a redox couple, with or without supporting electrolyte, strictly neutral or with
  // the double layers resolved; the far side the same electrode or a bath.
  echem(r) {
    const kind = r.pick(['Ag', 'Cu', 'Zn', 'redox']);
    const anion = r.pick(['NO3-', 'Cl-', 'SO42-']), support = r.chance(0.5) ? r.pick(['K+', 'Na+']) : null;
    const active = kind === 'redox' ? ['Fe3+', 'Fe2+'] : [{ Ag: 'Ag+', Cu: 'Cu2+', Zn: 'Zn2+' }[kind]];
    const names = [...active, ...(support ? [support] : []), anion], z = (n) => IONS[n].z;
    const c = {};
    for (const n of active) c[n] = r.log(-1, 2);
    if (support) c[support] = r.log(1, 3);
    c[anion] = names.filter((n) => n !== anion).reduce((t, n) => t + z(n) * c[n], 0) / -z(anion);
    const epsr = r.pick([0, 78.3]), L = r.log(-6, -3.5);
    const I = 0.5 * names.reduce((t, n) => t + z(n) * z(n) * c[n], 0), debye = Math.sqrt((Math.max(epsr, 1) * EPS0 * RT) / (2 * FARADAY * FARADAY * I));
    const electrodeMetal = kind === 'redox' ? 'Pt' : kind;
    const rx = kind === 'redox' ? half('Fe3+ + e- = Fe2+') : half(`${active[0]} + ${-z(anion) > 0 && z(active[0]) === 2 ? '2 ' : ''}e- = ${kind}(s)`, { [`${kind}(s)`]: 0 });
    const reaction = { ...rx, k0: r.log(-5, -1), alpha: r.u(0.3, 0.7) };
    const face = epsr > 0 ? { phi: { type: 'capacitive', C: r.u(0.1, 0.5) }, zeroCharge: r.u(-0.3, 0.3), reactions: [reaction] } : { reactions: [reaction] };
    const farBath = r.chance(0.3);
    const def = build({
      T,
      library: [aqueous(names, { epsr }), metal(electrodeMetal)],
      stack: [
        ohmic(0, ['e-']),
        layer(electrodeMetal, 1e-6),
        face,
        layer('water', L, { c0: { ...c } }),
        ...(farBath ? [bath({ ...c }, anion)] : [face, layer(electrodeMetal, 1e-6), ohmic(0, ['e-'])]),
      ],
      grid: epsr > 0 ? { hmin: debye / r.u(3, 10), hmax: Math.max(L / r.u(20, 80), debye), ratio: r.u(1.08, 1.2) } : { hmin: L / r.u(200, 2000), hmax: L / r.u(20, 80), ratio: r.u(1.08, 1.2) },
    });
    return {
      def,
      params: { kind, names, c, epsr, L, debye, farBath, k0: reaction.k0 },
      plan: { side: 'right', V: r.u(-0.2, 0.2), equilibrium: !farBath, step: r.u(0.005, 0.1), t: r.log(-3, 1), Ifrac: r.u(-0.9, 0.9) },
    };
  },

  // Liquid junctions: two baths of the same salts in different proportions, the far one floating
  // (the junction's potential), strictly neutral or not; steady, and grown from a sharp boundary.
  junction(r) {
    const cations = ['H+', 'Li+', 'Na+', 'K+', 'Cu2+', 'Zn2+'], anions = ['Cl-', 'NO3-', 'SO42-'];
    const cs = [...new Set([r.pick(cations), r.pick(cations)])], as = [...new Set([r.pick(anions), ...(r.chance(0.3) ? [r.pick(anions)] : [])])];
    const names = [...cs, ...as], z = (n) => IONS[n].z;
    const side = () => {
      const c = {};
      for (const n of cs) c[n] = r.log(-3, 3);
      const plus = cs.reduce((t, n) => t + z(n) * c[n], 0), w = as.map(() => r.u(0.2, 1)), wz = as.reduce((t, n, k) => t + w[k] * -z(n), 0);
      as.forEach((n, k) => (c[n] = (plus * w[k]) / wz));
      return c;
    };
    const left = side(), right = side(), epsr = r.pick([0, 0, 78.3]), L = r.log(-5, -3);
    const minI = Math.min(...[left, right].map((c) => 0.5 * names.reduce((t, n) => t + z(n) * z(n) * c[n], 0)));
    const debye = Math.sqrt((Math.max(epsr, 1) * EPS0 * RT) / (2 * FARADAY * FARADAY * minI));
    const sharp = r.chance(0.5);
    const lib = aqueous(names, { epsr });
    const regions = sharp
      ? [{ material: 'water', length: L / 2, c0: { ...left } }, { material: 'water', length: L / 2, c0: { ...right } }]
      : [{ material: 'water', length: L, c0: Object.fromEntries(names.map((n) => [n, { x: [0, 0.45 * L, 0.55 * L, L], values: [left[n], left[n], right[n], right[n]] }])) }];
    const def = {
      T,
      species: lib.species,
      materials: lib.materials,
      regions,
      contacts: { left: { bath: { c: { ...left }, reference: as[0] }, V: 0 }, right: { bath: { c: { ...right }, reference: as[0] }, I: 0 } },
      grid: epsr > 0 ? { hmin: Math.min(debye / 3, L / 200), hmax: L / r.u(20, 80), ratio: 1.1 } : { hmin: L / r.u(400, 4000), hmax: L / r.u(20, 80), ratio: r.u(1.08, 1.2) },
    };
    return { def, params: { names, left, right, epsr, L, sharp }, plan: { side: 'left', V: 0, equilibrium: false, step: 0, t: r.log(-2, 2), transientOnly: true } };
  },
};

// --- what's done with each device, and how it's judged
const flatness = (sol) => {
  let worst = 0;
  for (const name of Object.keys(sol.mu)) {
    const mu = Array.from(sol.mu[name]).filter(Number.isFinite);
    if (mu.length) worst = Math.max(worst, (Math.max(...mu) - Math.min(...mu)) / RT);
  }
  return worst;
};
const judged = (dev, sol) => {
  if (!sol.converged) return 'not converged';
  const c = check(dev, sol, { refine: false });
  return c.ok ? '' : c.items.filter((it) => it.ok === false).map((it) => `${it.name}: ${it.summary}`).join('; ');
};
const drive = (side, V) => ({ contacts: { [side]: { V } } });
const sameCurrent = (a, b, scale) => Math.abs(a - b) <= 1e-6 * Math.max(Math.abs(a), Math.abs(b), scale) + 1e-12;

function runCase(family, k) {
  const r = rng(hash(family) + k), { def, params, plan } = families[family](r), out = { family, k, params, plan, scenarios: {} };
  if (flag('def')) writeFileSync(flag('def'), JSON.stringify({ def, plan }, null, 1));
  const note = (name, error) => (out.scenarios[name] = error || 'ok');
  const attempt = (name, fn) => {
    try {
      note(name, fn());
    } catch (e) {
      note(name, `threw: ${String(e.message).split('\n')[0].slice(0, 160)}`);
    }
  };
  let dev;
  try {
    dev = new Device(def);
  } catch (e) {
    note('build', `threw: ${String(e.message).split('\n')[0].slice(0, 160)}`);
    return out;
  }
  if (plan.transientOnly) {
    attempt('steady', () => {
      const d = new Device(def);
      return judged(d, d.solve());
    });
    attempt('transient', () => {
      const d = new Device(def);
      const s = d.advance(plan.t, { budgetMs: 4000, maxSteps: 20000 });
      if (!s.converged) return `failed at t = ${s.time.toPrecision(3)} s (${s.steps} steps)`;
      if (!s.done) return `too slow: t = ${s.time.toPrecision(3)} of ${plan.t.toPrecision(3)} s in ${s.steps} steps`;
      return judged(d, s);
    });
    return out;
  }
  // 1. Cold, at no bias: equilibrium where nothing drives it.
  attempt('cold at 0 V', () => {
    const sol = dev.solve();
    const bad = judged(dev, sol);
    if (bad) return bad;
    if (plan.equilibrium) {
      const f = flatness(sol);
      if (f > 1e-6) return `not equilibrium: a level varies by ${f.toExponential(1)} RT`;
      if (Math.abs(sol.current) > 1e-9 * (1 + Math.abs(sol.current))) return `a current at equilibrium: ${sol.current}`;
    }
    return '';
  });
  // 2. Cold, at a bias.
  let cold;
  attempt('cold at bias', () => {
    const d = new Device(def);
    d.set(drive(plan.side, plan.V));
    cold = d.solve();
    return judged(d, cold);
  });
  // 3. Warm, swept there in steps, and the same answer.
  attempt('warm sweep', () => {
    const d = new Device(def);
    let sol = d.solve();
    const n = Math.max(1, Math.ceil(Math.abs(plan.V) / 0.1));
    for (let j = 1; j <= n; j++) {
      d.set(drive(plan.side, (plan.V * j) / n));
      sol = d.solve();
      if (!sol.converged) return `not converged at ${((plan.V * j) / n).toFixed(3)} V`;
    }
    const bad = judged(d, sol);
    if (bad) return bad;
    if (cold?.converged) {
      const gross = Math.max(...Object.values(sol.terminals).map((t) => Math.abs(t.current)));
      if (!sameCurrent(sol.current, cold.current, gross * 1e-3)) return `warm ${sol.current} vs cold ${cold.current} A/m²`;
      const D = (s) => s.interfaces.at(-1)?.D ?? 0;
      if (plan.cap && Math.abs(D(sol) - D(cold)) > 1e-6 * plan.cap) return `gate charge warm ${D(sol)} vs cold ${D(cold)}`;
    }
    return '';
  });
  // 3b. Open circuit (a solar cell's V_oc), from cold.
  if (plan.openCircuit) {
    attempt('open circuit', () => {
      const d = new Device(def);
      d.set({ contacts: { [plan.side]: { I: 0 } } });
      return judged(d, d.solve());
    });
  }
  // 3c. Driven by a current: a fraction of what passes at ±1 V (kinetics or transport can hold
  // that to nanoamps), found by a warm sweep.
  if (plan.Ifrac !== undefined) {
    attempt('current-driven', () => {
      const d = new Device(def);
      d.solve();
      const V1 = Math.sign(plan.Ifrac);
      let sol;
      for (let j = 1; j <= 10; j++) {
        d.set(drive(plan.side, (V1 * j) / 10));
        sol = d.solve();
        if (!sol.converged) return `no steady state at ${((V1 * j) / 10).toFixed(1)} V to find the limit`;
      }
      const target = Math.abs(plan.Ifrac) * sol.terminals[plan.side].current;
      const e = new Device(def);
      e.solve();
      e.set({ contacts: { [plan.side]: { I: target } } });
      const s1 = e.solve();
      const bad = judged(e, s1);
      if (bad) return bad;
      const I = s1.terminals[plan.side].current;
      return Math.abs(I - target) <= 1e-6 * Math.abs(target) + 1e-12 ? '' : `drove ${target}, got ${I} A/m²`;
    });
  }
  // 4. A step, then time.
  attempt('transient', () => {
    const d = new Device(def);
    d.solve();
    d.set(drive(plan.side, plan.step));
    const s = d.advance(plan.t, { budgetMs: 4000, maxSteps: 20000 });
    if (!s.converged) return `failed at t = ${s.time.toPrecision(3)} s (${s.steps} steps)`;
    if (!s.done) return `too slow: t = ${s.time.toPrecision(3)} of ${plan.t.toPrecision(3)} s in ${s.steps} steps`;
    return judged(d, s);
  });
  // 5. The impedance about equilibrium: passive, and at the lowest frequency the conductance
  // of steady solves either side.
  if (plan.equilibrium) {
    attempt('impedance', () => {
      const d = new Device(def);
      const sol = d.solve();
      if (!sol.converged) return 'no equilibrium to linearise about';
      // (well below the slowest diffusion's corner, L²/D)
      const fDC = Math.min(1e-6, 1e-3 / (2 * Math.PI * d._solver.slowestTime()));
      const fs = [fDC, 1e-2, 1e2, 1e6, 1e9], Z = d.impedance(fs, { terminal: plan.side }).Z;
      for (let j = 0; j < fs.length; j++) {
        const re = Z.re[j], im = Z.im[j], m = Math.hypot(re, im);
        if (!Number.isFinite(m)) return `non-finite at ${fs[j]} Hz`;
        if (re < -1e-6 * m) return `not passive at ${fs[j]} Hz: Re Z = ${re.toExponential(2)} of |Z| ${m.toExponential(2)}`;
        if (plan.cap && -1 / (2 * Math.PI * fs[j] * im) > plan.cap * 1.001) return `C above C_ox at ${fs[j]} Hz`;
      }
      // (A MOS capacitor's G is zero, and its series resistance, generation's, can make the
      // corner hours long.)
      if (plan.cap) return '';
      const dV = 1e-4, I = (V) => {
        const e = new Device(def);
        e.set(drive(plan.side, V));
        const s = e.solve();
        return s.converged ? s.terminals[plan.side].current : NaN;
      };
      const G = (I(dV) - I(-dV)) / (2 * dV);
      if (!Number.isFinite(G)) return 'no steady solves either side';
      // …and below the corner G/C, where that's lower still (a junction leaking only minority
      // carriers: hours).
      const admittance = (re, im) => ({ re: re / (re * re + im * im), im: -im / (re * re + im * im) });
      let f = fs[0], Y = admittance(Z.re[0], Z.im[0]);
      const tau = Y.im / (2 * Math.PI * f) / Math.abs(G);
      if (2 * Math.PI * f * tau > 1e-3) {
        f = 1e-3 / (2 * Math.PI * tau);
        const z = d.impedance([f], { terminal: plan.side }).Z;
        Y = admittance(z.re[0], z.im[0]);
      }
      if (Math.abs(Y.re - G) > 1e-3 * Math.hypot(Y.re, Y.im, G) + 1e-12) return `at ${f.toPrecision(2)} Hz Re Y = ${Y.re.toExponential(3)}, steady dI/dV = ${G.toExponential(3)} S/m²`;
      return '';
    });
  }
  return out;
}

// --- run
const names = onlyFamily ? [onlyFamily] : Object.keys(families);
if (onlyFamily && !families[onlyFamily]) throw new Error(`no family '${onlyFamily}' (known: ${Object.keys(families).join(', ')})`);
const tally = {}, failures = [];
const t0 = performance.now();
for (const family of names) {
  const ks = onlyCase === undefined ? Array.from({ length: N }, (_, k) => k) : [onlyCase];
  for (const k of ks) {
    const res = runCase(family, k);
    for (const [s, v] of Object.entries(res.scenarios)) {
      const key = `${family}: ${s}`;
      tally[key] ??= { ok: 0, of: 0 };
      tally[key].of++;
      if (v === 'ok') tally[key].ok++;
      else failures.push(`${family} ${k} · ${s}: ${v}`);
    }
    if (onlyCase !== undefined) console.log(JSON.stringify(res, null, 1));
  }
}
const ms = performance.now() - t0;
const width = Math.max(...Object.keys(tally).map((s) => s.length));
for (const [key, { ok, of }] of Object.entries(tally)) console.log(`${key.padEnd(width)}  ${String(ok).padStart(4)} / ${of}`);
const total = Object.values(tally).reduce((t, v) => ({ ok: t.ok + v.ok, of: t.of + v.of }), { ok: 0, of: 0 });
console.log(`\n${total.ok} of ${total.of} passed (${(100 * total.ok / total.of).toFixed(1)}%), ${(ms / 1000).toFixed(1)} s`);
if (failures.length) console.log(`\nFailures (rerun one with: npm run stress -- <family> <case>):\n  ${failures.join('\n  ')}`);
if (save) {
  writeFileSync(new URL('./stress.json', import.meta.url), JSON.stringify({ n: N, tally, failures }, null, 2) + '\n');
  console.log('Saved bench/stress.json');
}
