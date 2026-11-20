// Electrochemical devices that combine features at random: aqueous electrolytes from the kit's ion
// data (1–3 salts, strictly neutral or with their double layers resolved, with or without a
// supporting salt, ideal, Debye–Hückel or crowded on a lattice), in one to three compartments
// whose faces block some species; at either end an electrode (a metal in its own ions, a couple on
// an inert metal, silver–silver chloride, hydrogen on platinum, or a blocking metal; Butler–Volmer
// or Marcus–Hush–Chidsey; a Stern layer, a diffuse layer below the grid, adsorbates or a tethered
// couple on its face; a second couple now and then, for a mixed potential), a bath read by a
// reference electrode, a species or nothing, an insertion host (OCV table or regular solution),
// or a closed wall; bulk complexations (a CE mechanism where the complexed ion is an electrode's)
// and an acetate buffer; internal ports (a reference electrode or a salt bridge at I: 0, a porous
// electrode through a window with its double layer and adsorbates, a capacitance through a
// window, an ion link held at a voltage), held at V, behind R, or floating; a spherical or
// cylindrical microelectrode; a series resistance; the far terminal driven by a current; a
// stirred compartment; a temperature other than 25 °C.
//
// Every held terminal is set at its rest level, computed from the compositions as the solver
// reads them (activities included), so a device whose rest is known exactly is at rest as
// defined. Where it isn't known exactly (a closed compartment whose double layers take ions at
// rest, two couples on one electrode, two baths of different compositions, …) it isn't claimed.
// The plan carries, beside the README's fields, `rest`: no current anywhere and every phase in
// equilibrium within itself, which a cell between unlike electrodes at its open-circuit voltage
// is; `equilibrium` is rest with every species' level flat through the whole device, the
// electrons' in every metal included (so not that cell). Concentrations are mol/m³ (1 mM).
// Checked (scripts kept outside the repo): over 400 cases no rest, equilibrium or unique claim
// is wrong. plan.slow against transients after a step (60 cases): 36 settle between 0.03 and
// 10 slow; 20 settle sooner (it's generous, most where a porous electrode or a held link is
// counted); 4 are 5–20× short (closed cells with slow kinetics behind blocking faces, a buffer).
import { EPS0, FARADAY, GAS_CONSTANT, AVOGADRO, units } from '../../../src/index.js';
import { build, layer, ohmic, bath, aqueous, metal, IONS } from '../../../src/kit.js';

const AgCl = -109.789e3; // ΔfG° of AgCl(s), as the kit's AgAgCl has it

// Acetate and acetic acid (ΔfG° of the aqueous species from the NBS tables, from memory: they give
// pKa 4.757, as measured; D of acetate from the CRC's limiting conductivity, HAc's a round value).
const EXTRA = { 'Ac-': { z: -1, D: 1.089e-9, mu0: -369.31e3 }, HAc: { z: 0, D: 1.2e-9, mu0: -396.46e3 } };

// The electrodes a face (or a port) can be: its metal, the solution species its couple needs, and
// the couple. `one`: a single electron, so Marcus–Hush–Chidsey applies.
const COUPLES = {
  Ag: { metals: ['Ag'], ions: ['Ag+'], equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, one: true },
  Cu: { metals: ['Cu'], ions: ['Cu2+'], equation: 'Cu2+ + 2 e- = Cu(s)', fixed: { 'Cu(s)': 0 }, one: false },
  Zn: { metals: ['Zn'], ions: ['Zn2+'], equation: 'Zn2+ + 2 e- = Zn(s)', fixed: { 'Zn(s)': 0 }, one: false },
  Fe: { metals: ['Pt', 'Au'], ions: ['Fe3+', 'Fe2+'], equation: 'Fe3+ + e- = Fe2+', fixed: {}, one: true },
  AgCl: { metals: ['Ag'], ions: ['Cl-'], equation: 'AgCl(s) + e- = Ag(s) + Cl-', fixed: { 'AgCl(s)': AgCl, 'Ag(s)': 0 }, one: true },
  H: { metals: ['Pt'], ions: ['H+'], equation: '2 H+ + 2 e- = H2', fixed: { H2: 0 }, one: false },
};
const LIGANDS = { 'Fe3+': ['Cl-', 'Ac-'], 'Cu2+': ['Cl-', 'Ac-'], 'Zn2+': ['Cl-', 'Ac-'], 'Fe2+': ['Cl-'] };

// Signed stoichiometry of an equation (products positive).
function nus(equation) {
  const out = {};
  const [l, r] = equation.split('=').map((s) => s.trim());
  for (const [side, sign] of [[l, -1], [r, 1]]) {
    if (side === '0') continue;
    for (const term of side.split(' + ')) {
      const m = term.trim().match(/^(\d+) (.+)$/);
      const name = m ? m[2] : term.trim();
      out[name] = (out[name] ?? 0) + sign * (m ? +m[1] : 1);
    }
  }
  return out;
}
const chargeTag = (z) => (z === 0 ? '' : `${Math.abs(z) === 1 ? '' : Math.abs(z)}${z > 0 ? '+' : '-'}`);
const fixedOf = (c) => (Object.keys(c.fixed).length ? { fixed: { ...c.fixed } } : {});

export default {
  name: 'echem',
  make(r) {
    const features = [];
    const on = (f) => {
      if (!features.includes(f)) features.push(f);
      return true;
    };
    // 25 °C, or another temperature now and then (the ions' D move with it, their μ° stay 25 °C's,
    // as aqueous() has them, and the device warns where it reads them)
    const T = r.chance(0.15) ? r.u(278.15, 318.15) : 298.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;
    if (T !== 298.15) on('another temperature');

    // --- what's at each end
    const leftEnd = r.chance(0.78) ? 'electrode' : 'bath';
    let rightEnd = r.pick(['electrode', 'electrode', 'electrode', 'bath', 'bath', 'bath', 'host', 'wall']);
    if (leftEnd === 'bath' && rightEnd === 'bath' && r.chance(0.5)) rightEnd = 'wall';
    const epsr = r.chance(0.35) ? 78.3 : 0;
    const neutral = epsr === 0;
    if (rightEnd === 'wall' && !neutral) rightEnd = 'bath'; // (a wall needs a porous electrode, which needs ε = 0)
    if (!neutral) on('double layers resolved');
    const geometry = leftEnd === 'electrode' && rightEnd === 'bath' && r.chance(0.3) ? r.pick(['spherical', 'cylindrical']) : 'planar';
    if (geometry !== 'planar') on(geometry);

    // --- electrodes at faces, and a porous electrode's couple (Ag⁺ and Cl⁻ never together: AgCl
    // would precipitate)
    const kinds = ['Ag', 'Cu', 'Zn', 'Fe', 'Fe', 'AgCl', 'H', 'blocking'];
    const fits = (k, used) => !used.some((u) => (k === 'Ag' && u === 'AgCl') || (k === 'AgCl' && u === 'Ag'));
    const pickKind = (used, notBlocking = false) => {
      for (;;) {
        const k = r.pick(kinds);
        if (fits(k, used) && !(notBlocking && k === 'blocking')) return k;
      }
    };
    const left = leftEnd === 'electrode' ? { kind: pickKind([]) } : null;
    const right = rightEnd === 'electrode' ? { kind: left && r.chance(0.4) ? left.kind : pickKind(left ? [left.kind] : []) } : null;
    if (left && right && left.kind === 'blocking' && right.kind === 'blocking') right.kind = pickKind([], true);
    if (left?.kind === 'blocking' && !right && rightEnd !== 'bath') left.kind = 'Fe'; // (something must fix the solution's level)
    const used = [left?.kind, right?.kind].filter((k) => k && k !== 'blocking');
    const portPlans = [];
    if (neutral && (rightEnd === 'wall' || r.chance(0.25))) {
      const k = pickKind(used, true);
      used.push(k);
      portPlans.push({ type: 'porous', kind: k });
    }

    // --- the electrolyte's species
    const ions = new Set();
    for (const k of used) COUPLES[k].ions.forEach((n) => ions.add(n));
    if (rightEnd === 'host') ions.add('Li+');
    const hasAg = ions.has('Ag+');
    const support = r.chance(0.55) || rightEnd === 'host' ? r.pick(['K+', 'Na+']) : null;
    if (support) ions.add(support);
    // bulk reactions
    const ligandsOf = (n) => (LIGANDS[n] ?? []).filter((l) => !(l === 'Cl-' && hasAg));
    const complexable = [...ions].filter((n) => ligandsOf(n).length);
    const complexation = complexable.length > 0 && r.chance(0.35) ? r.pick(complexable) : null;
    const ligand = complexation ? r.pick(ligandsOf(complexation)) : null;
    const buffer = r.chance(0.2);
    if (buffer || ligand === 'Ac-') ions.add('Ac-');
    if (buffer) ions.add('H+'), ions.add('HAc');
    if (ligand) ions.add(ligand);
    // anions: a balancing one, set by neutrality (not the ligand), and maybe another
    const spectators = ['NO3-', 'ClO4-', 'SO42-', ...(hasAg ? [] : ['Cl-'])].filter((a) => a !== ligand);
    const balancer = r.pick(spectators);
    ions.add(balancer);
    if (r.chance(0.3)) {
      const more = spectators.filter((a) => a !== balancer && !ions.has(a));
      if (more.length) ions.add(r.pick(more));
    }
    const complex = complexation
      ? { name: `${complexation.replace(/\d*\+$/, '')}${ligand === 'Cl-' ? 'Cl' : 'Ac'}${chargeTag(IONS[complexation].z - 1)}`, z: IONS[complexation].z - 1, M: complexation, L: ligand, logK: r.u(0.5, 2.5), D: r.u(0.5, 0.9) * 1e-9 }
      : null;
    if (complex) ions.add(complex.name), on('complexation');
    if (complex && used.some((k) => COUPLES[k].ions.includes(complex.M))) on('CE mechanism');
    if (buffer) on('acetate buffer');
    const z = (n) => (n === complex?.name ? complex.z : (IONS[n] ?? EXTRA[n]).z);
    const D = (n) => (n === complex?.name ? complex.D : (IONS[n] ?? EXTRA[n]).D);
    if (!support && (ligand || buffer || [...ions].some((n) => z(n) < 0 && n !== balancer) || ![...ions].some((n) => z(n) > 0))) ions.add('K+');
    const supportCation = support ?? (ions.has('K+') ? 'K+' : null);
    if (supportCation) on('supporting salt');
    const names = [...ions];
    const dataIons = names.filter((n) => IONS[n]);
    // (H⁺ moves by proton hopping, which aqueous() leaves to us at another T: its D from the 25 °C
    // one by about 1.4% per kelvin, a round figure from its conductivity's rise)
    const lib = aqueous(dataIons, { epsr: neutral ? 0 : undefined, ...(T !== 298.15 ? { T, D: dataIons.includes('H+') ? { 'H+': IONS['H+'].D * (1 + 0.014 * (T - 298.15)) } : {} } : {}) });
    const water = lib.materials.water;
    const mu0 = Object.fromEntries(dataIons.map((n) => [n, IONS[n].mu0]));
    for (const n of names) {
      if (IONS[n]) continue;
      lib.species.push({ name: n, z: z(n), cRef: 1000 });
      if (n !== complex?.name) mu0[n] = EXTRA[n].mu0;
    }
    if (complex) mu0[complex.name] = mu0[complex.M] + mu0[complex.L] - RT * Math.log(10 ** complex.logK);
    for (const n of names) if (!IONS[n]) water.species[n] = { D: D(n) * (T === 298.15 ? 1 : water.species[balancer].D / IONS[balancer].D), mu0: mu0[n] };
    const charged = names.filter((n) => z(n) !== 0);
    const anions1 = charged.filter((n) => z(n) === -1);

    // --- activity models, and the compositions: free concentrations drawn, the complex and HAc
    // from their equilibria, the balancing anion from neutrality, all in activities (so a bath and
    // every c0 sit at their reactions' equilibrium exactly)
    let stats = null;
    const statPick = r.u(0, 1);
    if (neutral && statPick < 0.2) stats = { type: 'debye-huckel', species: charged, epsr: 78.3, a: r.u(3, 6) * 1e-10 };
    else if (statPick > 0.85) stats = { type: 'lattice', species: charged, cMax: 1 / (AVOGADRO * (r.u(0.5, 0.66) * 1e-9) ** 3) };
    const excessOf = (c) => {
      const ex = Object.fromEntries(names.map((n) => [n, 0]));
      if (!stats) return ex;
      if (stats.type === 'debye-huckel') {
        const eps = stats.epsr * EPS0, lB = (FARADAY * FARADAY) / (4 * Math.PI * eps * AVOGADRO * RT);
        let q = 0;
        for (const n of stats.species) q += z(n) * z(n) * c[n];
        const kappa = Math.sqrt((FARADAY * FARADAY * q) / (eps * RT));
        for (const n of stats.species) ex[n] = (-z(n) * z(n) * lB * kappa) / (2 * (1 + kappa * stats.a));
      } else {
        let s = 0;
        for (const n of stats.species) s += c[n];
        for (const n of stats.species) ex[n] = Math.log(stats.cMax / (stats.cMax - s));
      }
      return ex;
    };
    const zeta = (c, ex, n) => Math.log(c[n] / 1000) + ex[n];
    const compose = () => {
      const c = {};
      for (const n of names) {
        if (n === balancer || n === complex?.name || n === 'HAc') continue;
        if (n === supportCation) c[n] = support ? r.log(1.7, 3) : r.log(0.5, 2);
        else if (n === 'H+') c[n] = buffer ? r.log(-2, 0.5) : r.log(0, 2);
        else if (n === 'Li+') c[n] = r.log(1, 2.5);
        else if (z(n) < 0) c[n] = r.log(0.3, 2);
        else c[n] = r.log(0, 2);
      }
      let ex = Object.fromEntries(names.map((n) => [n, 0]));
      for (let it = 0; it < 300; it++) {
        if (complex) c[complex.name] = 1000 * Math.exp((mu0[complex.M] + mu0[complex.L] - mu0[complex.name]) / RT + zeta(c, ex, complex.M) + zeta(c, ex, complex.L) - ex[complex.name]);
        if (buffer) c.HAc = 1000 * Math.exp((mu0['H+'] + mu0['Ac-'] - mu0.HAc) / RT + zeta(c, ex, 'H+') + zeta(c, ex, 'Ac-') - ex.HAc);
        let q = 0;
        for (const n of names) if (n !== balancer) q += z(n) * c[n];
        if (q < 1) {
          // (too little positive charge for the balancing anion: more of a cation)
          c[supportCation ?? names.find((n) => z(n) > 0 && n !== complex?.name)] += 2 - q;
          continue;
        }
        const cb = q / -z(balancer);
        const settled = it > 0 && Math.abs(cb - c[balancer]) <= 1e-15 * cb;
        c[balancer] = cb;
        const ex1 = excessOf(c);
        const moved = names.some((n) => Math.abs(ex1[n] - ex[n]) > 1e-15);
        ex = ex1;
        if (settled && !moved) break;
      }
      return c;
    };
    const ionic = (c) => 0.5 * names.reduce((t, n) => t + z(n) * z(n) * c[n], 0);
    const total = (c) => names.reduce((t, n) => t + c[n], 0);

    // --- compartments
    const nComp = r.pick([1, 1, 1, 2, 2, 3]);
    const Ltot = geometry !== 'planar' ? r.log(-5, -3.7) : neutral ? r.log(-5.3, -3.5) : r.log(-6, -4.3);
    const w = Array.from({ length: nComp }, () => r.u(0.5, 1)), ws = w.reduce((a, b) => a + b, 0);
    const lengths = w.map((x) => (Ltot * x) / ws);
    const sameComp = nComp === 1 || r.chance(0.5);
    let comps;
    for (let attempt = 0; ; attempt++) {
      comps = [compose()];
      for (let q = 1; q < nComp; q++) comps.push(sameComp ? { ...comps[0] } : compose());
      const tooStrong = stats?.type === 'debye-huckel' ? Math.max(...comps.map(ionic)) > 300 : stats?.type === 'lattice' ? Math.max(...comps.map(total)) > 0.4 * stats.cMax : false;
      if (!tooStrong) break;
      if (attempt >= 2) stats = null; // (dilute solutions only, for Debye–Hückel; a lattice not near full)
    }
    if (stats) {
      water.statistics = [stats.type === 'debye-huckel' ? { type: 'debye-huckel', species: stats.species, epsr: stats.epsr, a: stats.a } : { type: 'lattice', species: stats.species, cMax: stats.cMax }];
      on(stats.type);
    }
    const exs = comps.map(excessOf);
    // faces between compartments, each blocking some species (at least one ion still crosses)
    const walls = [];
    for (let q = 1; q < nComp; q++) {
      const blocked = r.chance(0.75) ? names.filter(() => r.chance(0.4)) : [];
      if (charged.every((n) => blocked.includes(n))) blocked.splice(blocked.indexOf(r.pick(charged)), 1);
      walls.push(blocked);
    }
    const blockingWalls = walls.some((b) => b.length);
    if (blockingWalls) on('faces blocking species');
    if (nComp > 1) on(sameComp ? 'compartments alike' : 'compartments differ');

    // φ in each compartment, φ₀ = 0 for now (every level moves with it, and the anchor shifts
    // them all at the end), across each face from an ion that crosses it. Exact where one ion alone
    // crosses (neutrality then holds each side's composition) or the compositions are alike.
    let exact = true;
    const phis = [0];
    for (let q = 1; q < nComp; q++) {
      const cross = charged.filter((n) => !walls[q - 1].includes(n));
      const s = cross[0];
      phis.push(phis[q - 1] + (RT * (zeta(comps[q - 1], exs[q - 1], s) - zeta(comps[q], exs[q], s))) / (z(s) * FARADAY));
      if (!sameComp) {
        const neutralsCross = names.filter((n) => z(n) === 0 && !walls[q - 1].includes(n));
        if (cross.length > 1 || neutralsCross.some((n) => Math.abs(comps[q][n] / comps[q - 1][n] - 1) > 1e-12)) exact = false;
        if (!neutral) exact = false; // (a pinned face between them: the double layers either side take ions)
      }
    }
    const mubar = (q, n) => mu0[n] + RT * zeta(comps[q], exs[q], n) + z(n) * FARADAY * phis[q];
    const levelIn = (q, couple) => {
      const nu = nus(couple.equation);
      let S = 0;
      for (const [n, v] of Object.entries(nu)) if (n !== 'e-') S += v * (couple.fixed[n] ?? mubar(q, n));
      return S / (nu['e-'] * FARADAY);
    };

    // --- kinetics
    const kinetics = (one, scale = 1) => {
      const k0 = scale * r.log(-5, -1.5);
      if (one && r.chance(0.35)) return on('Marcus–Hush–Chidsey'), { k0, lambda: units.eV(r.u(0.3, 1.2)) };
      return { k0, alpha: r.u(0.3, 0.7) };
    };
    const coupleOf = (kind) => {
      const c = COUPLES[kind];
      return { equation: c.equation, fixed: kind === 'H' ? { H2: RT * Math.log(r.log(-4, 0)) } : { ...c.fixed }, one: c.one, kind };
    };

    // --- electrodes at the ends (their faces are written once the voltages are known)
    const slow = [];
    const kappaOf = (c) => ((FARADAY * FARADAY) / RT) * names.reduce((t, n) => t + z(n) * z(n) * D(n) * c[n], 0);
    const solR = Ltot / Math.min(...comps.map(kappaOf));
    const metals = new Set();
    let mixed = false, closedSurface = false;
    const setupElectrode = (e, q, side) => {
      e.q = q;
      e.metal = e.kind === 'blocking' ? r.pick(['Pt', 'Au']) : r.pick(COUPLES[e.kind].metals);
      metals.add(e.metal);
      e.couple = e.kind === 'blocking' ? null : coupleOf(e.kind);
      e.level = e.couple ? levelIn(q, e.couple) : null;
      e.law = neutral ? (e.kind === 'blocking' ? r.pick(['capacitive', 'diffuse']) : r.pick(['neutral', 'capacitive', 'capacitive', 'diffuse'])) : 'capacitive';
      if (e.law === 'diffuse' && stats?.type === 'debye-huckel') e.law = 'capacitive'; // (a diffuse layer's ions are ideal or on a lattice)
      e.C = r.u(0.1, 0.5);
      e.zcOffset = r.chance(0.5) ? 0 : r.u(-0.25, 0.25);
      e.surface = e.kind === 'H' && r.chance(0.5) ? 'volmer' : r.chance(0.25) && anions1.length ? r.pick(['adsorbate', 'adsorbate', 'tethered']) : null;
      if (e.surface === 'tethered') closedSurface = true;
      if (e.law !== 'neutral') on(e.law === 'diffuse' ? 'diffuse layer' : 'Stern layer');
      if (e.surface) on(e.surface === 'volmer' ? 'Volmer–Tafel' : e.surface === 'adsorbate' ? 'adsorbate at a face' : 'tethered couple');
      e.rx = e.couple && e.surface !== 'volmer' ? kinetics(e.couple.one) : null;
      // a second couple on the same metal now and then: a mixed potential
      e.second = null;
      if (e.couple && r.chance(0.12)) {
        const others = ['Fe', 'H', ...(e.metal === 'Ag' ? ['AgCl'] : [])].filter((k) => k !== e.kind && COUPLES[k].ions.every((n) => ions.has(n)));
        if (others.length) {
          const couple = coupleOf(r.pick(others));
          e.second = { couple, rx: kinetics(couple.one, 1e-2) };
          mixed = true;
          on('mixed potential');
        }
      }
      on(`${side} electrode: ${e.kind}`);
    };
    if (left) setupElectrode(left, 0, 'left');
    if (right) setupElectrode(right, nComp - 1, 'right');

    // --- the host
    let host = null;
    if (rightEnd === 'host') {
      const cMax = r.u(2e4, 5e4), L = r.log(-6.7, -5.7), DLi = r.log(-14, -12.5);
      const type = r.chance(0.5) ? 'ocv' : 'rk';
      const E0 = r.u(2.6, 3.6); // V against Li metal, which puts the host within ~0.5 V of an aqueous Li⁺'s level
      let model, x0, Ex0, Om = 0;
      if (type === 'ocv') {
        // an ideal-solution curve made steeper (a) and sloped (b): its points, and the start at one of them
        const a = r.u(1, 1.6), b = r.u(0, 0.3), xs = Array.from({ length: 13 }, (_, k) => 0.04 + (0.92 * k) / 12);
        const E = xs.map((x) => E0 - a * VT * Math.log(x / (1 - x)) - b * (x - 0.5));
        const j = r.int(3, 9);
        x0 = xs[j], Ex0 = E[j];
        model = { type: 'insertion', species: ['Li+', 'e-'], cMax, ocv: { x: xs, E, muRef: 0 } };
        on('insertion host (OCV table)');
      } else {
        Om = r.u(-2, 1.5) * RT;
        x0 = r.u(0.15, 0.85);
        Ex0 = E0 - VT * Math.log(x0 / (1 - x0)) + (2 * Om * x0) / FARADAY;
        model = { type: 'insertion', species: ['Li+', 'e-'], cMax, A: [Om] };
        on('insertion host (regular solution)');
      }
      host = {
        type, cMax, L, DLi, x0,
        level: Ex0 + mubar(nComp - 1, 'Li+') / FARADAY, // μ_Li = μ̄_Li⁺ + μ̄_e⁻ = −F E(x)
        mat: { epsr: 0, species: { 'Li+': { D: DLi, mu0: type === 'ocv' ? 0 : -FARADAY * E0, cRef: cMax }, 'e-': { D: 1e-6, mu0: 0, cRef: cMax } }, statistics: [model] },
        k0: r.log(-6, -3),
      };
      const tf = type === 'rk' ? Math.max(0.1, 1 - (2 * Om * x0 * (1 - x0)) / RT) : 1;
      slow.push((L * L) / (DLi * tf));
      const Cchem = (FARADAY * FARADAY * cMax * x0 * (1 - x0) * L) / RT, Rct = RT / (FARADAY * FARADAY * host.k0 * Math.sqrt(comps[nComp - 1]['Li+'] / 1000));
      slow.push(Cchem * (Rct + solR));
    }

    // --- baths
    const AGAGCL = { equation: COUPLES.AgCl.equation, fixed: COUPLES.AgCl.fixed };
    const bathRef = () => {
      const opts = ['species', 'species', 'none'];
      if (names.includes('Cl-')) opts.push('AgAgCl', 'AgAgCl');
      const k = r.pick(opts);
      if (k === 'species') return { kind: k, species: r.pick(charged) };
      return { kind: k };
    };
    const bathLevel = (q, ref) => {
      // the bath's terminal voltage at its compartment's φ
      if (ref.kind === 'none') return phis[q];
      if (ref.kind === 'AgAgCl') return levelIn(q, AGAGCL);
      return mubar(q, ref.species) / (z(ref.species) * FARADAY);
    };
    const refArg = (ref) => (ref.kind === 'AgAgCl' ? AGAGCL : ref.kind === 'species' ? ref.species : undefined);
    const leftBath = leftEnd === 'bath' ? bathRef() : null;
    const rightBath = rightEnd === 'bath' ? bathRef() : null;
    if (leftBath) on(`bath read by ${leftBath.kind === 'species' ? 'a species' : leftBath.kind}`);
    if (rightBath) on(`bath read by ${rightBath.kind === 'species' ? 'a species' : rightBath.kind}`);
    if (leftBath && rightBath && !sameComp) exact = false;

    // --- ports
    const ports = [];
    const regionNames = lengths.map((_, q) => `w${q}`);
    if (r.chance(0.25)) portPlans.push({ type: 'ref' });
    if (r.chance(0.15)) portPlans.push({ type: 'cap' });
    if (r.chance(0.12)) portPlans.push({ type: 'link' });
    for (const p of portPlans) {
      const q = r.int(0, nComp - 1), Lq = lengths[q];
      p.q = q;
      if (p.type === 'porous') {
        // (a window with a surface keeps off a face that has one)
        const whole = r.chance(0.4);
        const from = whole ? 0 : r.u(0.1, 0.5) * Lq, to = whole ? Lq : from + r.u(0.2, 0.4) * Lq;
        const area = r.log(3, 5.5);
        const couple = coupleOf(p.kind);
        p.level = levelIn(q, couple);
        p.port = { name: 'porous', region: regionNames[q], from, to, terminal: 'e-', area, reactions: [{ equation: couple.equation, ...fixedOf(couple), ...kinetics(couple.one, 1e-3) }] };
        if (r.chance(0.5)) p.cap = { C: r.u(0.1, 0.4), offset: r.chance(0.5) ? 0 : r.u(-0.2, 0.2) }, on('porous double layer');
        if (r.chance(0.2) && !whole && anions1.length) p.ads = { X: r.pick(anions1), capacity: r.log(-6, -4.7), theta: r.u(0.05, 0.5) }, on('porous adsorbate');
        p.drive = geometry === 'planar' ? r.pick(['V', 'V', 'V', 'R', 'float']) : r.pick(['V', 'float']);
        on(`porous electrode (${p.drive === 'float' ? 'floating' : p.drive === 'R' ? 'behind R' : 'held'})`);
        const C = (p.cap?.C ?? 0) + (p.ads ? (FARADAY * FARADAY * p.ads.capacity * 0.25) / RT : 0);
        slow.push((area * Math.max(C, 0.01) * (to - from) ** 2) / kappaOf(comps[q]));
      } else if (p.type === 'ref') {
        // (a window of no width is the node there: at a face between compartments, where there is one)
        const point = q > 0 && r.chance(0.4);
        const x = point ? 0 : r.u(0.2, 0.8) * Lq, wd = point ? 0 : r.u(0.005, 0.03) * Lq;
        const avail = Object.keys(COUPLES).filter((k) => COUPLES[k].ions.every((n) => ions.has(n)));
        if (neutral && avail.length && r.chance(0.6)) {
          const couple = coupleOf(r.pick(avail));
          p.port = { name: 'ref', region: regionNames[q], from: x, to: x + wd, I: 0, terminal: 'e-', area: r.log(2, 4), reactions: [{ equation: couple.equation, ...fixedOf(couple), ...kinetics(couple.one) }] };
          on('reference electrode at I: 0');
        } else {
          const ion = r.pick(charged);
          p.port = { name: 'ref', region: regionNames[q], from: x, to: x + wd, I: 0, terminal: ion, species: { [ion]: { type: 'conductance', G: r.log(6, 10) } } };
          on('salt bridge at I: 0');
        }
      } else if (p.type === 'cap') {
        const from = r.u(0.1, 0.5) * Lq, to = from + r.u(0.1, 0.4) * Lq;
        p.cap = { C: r.u(0.05, 0.4), offset: r.chance(0.5) ? 0 : r.u(-0.2, 0.2) };
        p.port = { name: 'cap', region: regionNames[q], from, to, area: r.log(3, 5) };
        p.Vrel = r.u(-0.3, 0.3); // its voltage, from the solution's φ
        on('capacitance through a window');
        slow.push(p.port.area * p.cap.C * (to - from) * solR);
      } else {
        const from = r.u(0.1, 0.6) * Lq, to = from + r.u(0.05, 0.3) * Lq;
        p.ion = r.pick(charged);
        p.port = { name: 'link', region: regionNames[q], from, to, terminal: p.ion, species: { [p.ion]: { type: 'conductance', G: r.log(3, 7) } } };
        on('ion link held at V');
      }
      ports.push(p);
    }

    // --- the anchor: the left terminal at 0 V where it has a level, else the right; every other
    // held terminal at its rest level
    const leftLevel = left ? left.level : bathLevel(0, leftBath);
    const rightLevel = right ? right.level : rightBath ? bathLevel(nComp - 1, rightBath) : host ? host.level : null;
    const shift = leftLevel ?? rightLevel ?? ports.find((p) => p.type === 'porous')?.level ?? 0;
    const phi = (q) => phis[q] - shift;
    const VL = leftLevel === null ? 0 : leftLevel - shift;
    const VR = rightLevel === null ? (r.chance(0.5) ? VL : r.u(-0.3, 0.3)) : rightLevel - shift;
    // actual μ̄ of a solution species (in the anchored frame)
    const mu = (q, n) => mubar(q, n) - z(n) * FARADAY * shift;

    // drives: the right end behind R, or driven by a current, now and then
    const rightPasses = rightEnd !== 'wall' && !(right?.kind === 'blocking');
    const leftPasses = !(left?.kind === 'blocking');
    let rightDrive = { V: VR };
    if (rightPasses && r.chance(0.18)) rightDrive = { V: VR, R: r.log(-5, -2.5) }, on('series R');
    // (A steady current needs what one end uses to come from the other or a bath: not through a
    // face that blocks some species, not into a host that only stores it, not between unlike
    // couples in a closed cell, which would run one dry.)
    else if (rightPasses && leftPasses && !host && !blockingWalls && [left, right].every((e) => !e || e.rx) && (leftBath || rightBath || (left && right && left.kind === right.kind && !left.second && !right.second && !left.surface && !right.surface)) && r.chance(0.18)) {
      // a fraction of what each electrode can pass: its couple's most dilute ion's limiting current,
      // F D c / L, and its kinetics' (Marcus–Hush–Chidsey's rate constant levels off at about
      // k0 e^{λ/4RT}, Frumkin can starve it further: a few i0 where Butler–Volmer has no limit)
      let Imax = (FARADAY * 1e-9 * Math.min(...comps.map((c) => Math.min(...charged.map((n) => c[n]))))) / Ltot; // (between two baths: the most dilute ion's)
      for (const e of [left, right]) {
        if (!e) continue;
        const cs = COUPLES[e.kind].ions.map((n) => comps[e.q][n]), c = Math.min(...cs);
        const i0 = FARADAY * e.rx.k0 * (c / 1000);
        Imax = Math.min(Imax, (FARADAY * 0.7e-9 * c) / Ltot, e.rx.lambda ? 0.3 * i0 * Math.exp(e.rx.lambda / (4 * RT)) : 5 * i0);
      }
      rightDrive = { I: r.u(-0.2, 0.2) * Imax };
      on('driven by a current');
    }

    // --- faces at the electrodes
    const faceOf = (e, V) => {
      const f = {};
      if (e.law !== 'neutral') f.phi = { type: 'capacitive', C: e.C }, f.zeroCharge = V - phi(e.q) + e.zcOffset;
      if (e.law === 'diffuse') f.diffuse = { epsr: 78.4 };
      const reactions = [];
      if (e.rx) reactions.push({ equation: e.couple.equation, ...fixedOf(e.couple), ...e.rx });
      if (e.second) reactions.push({ equation: e.second.couple.equation, ...fixedOf(e.second.couple), ...e.second.rx });
      // each surface species starts at its coverage at rest, so no ions move to fill it
      const surface = {};
      if (e.surface === 'volmer') {
        const th = r.u(0.05, 0.5), capacity = r.log(-5.3, -4.7), k0 = r.log(-4, -1.5);
        surface.Hads = { mu0: mu(e.q, 'H+') - FARADAY * V - RT * Math.log(th / (1 - th)), capacity, theta0: th }; // μ = μ̄_H⁺ + μ̄_e⁻
        const volmer = kinetics(true);
        e.vk0 = volmer.k0;
        reactions.push({ equation: 'H+ + e- = Hads', ...volmer }, { equation: '2 Hads = H2', fixed: { H2: e.couple.fixed.H2 }, k0, alpha: 0.5 });
        slow.push(capacity / k0);
      }
      if (e.surface === 'adsorbate') {
        const X = r.pick(anions1), th = r.u(0.03, 0.4), capacity = r.log(-6, -4.7), k0 = r.log(-6, -3);
        surface.Xads = { mu0: mu(e.q, X) + FARADAY * V - RT * Math.log(th / (1 - th)), capacity, theta0: th }; // μ = μ̄_X⁻ − μ̄_e⁻
        reactions.push({ equation: `${X} = Xads + e-`, k0, alpha: 0.5 });
        slow.push(capacity / k0);
      }
      if (e.surface === 'tethered') {
        const X = r.pick(anions1), thR = r.u(0.05, 0.6), thO = r.u(0.05, 0.9 - thR), capacity = r.log(-6.3, -5), k0 = r.log(-7, -4);
        surface.Fc = { mu0: 0, capacity, theta0: thR };
        surface.FcX = { mu0: mu(e.q, X) + FARADAY * V + RT * Math.log(thR / thO), capacity, theta0: thO };
        reactions.push({ equation: `Fc + ${X} = FcX + e-`, k0, ...(r.chance(0.4) ? { lambda: units.eV(r.u(0.4, 1)) } : { alpha: 0.5 }) });
        slow.push(capacity / k0);
      }
      if (Object.keys(surface).length) f.surface = surface;
      if (reactions.length) f.reactions = reactions;
      // its time: the Stern layer charging through the solution and the reaction
      if (e.law !== 'neutral') slow.push(e.C * (solR + (e.rx ? RT / (FARADAY * FARADAY * e.rx.k0) : 0) + (rightDrive.R ?? 0)));
      return f;
    };

    // --- the stack (one compartment stirred now and then: eddy mixing, which carries no current)
    const stirred = r.chance(0.12) ? r.int(0, nComp - 1) : -1, stirD = r.log(-8.5, -6.5);
    if (stirred >= 0) on('stirred (mixing)');
    const stack = [];
    if (left) stack.push(ohmic(VL, ['e-']), layer(left.metal, 1e-6), faceOf(left, VL));
    else stack.push(bath({ ...comps[0] }, refArg(leftBath), VL));
    for (let q = 0; q < nComp; q++) {
      if (q > 0) stack.push(walls[q - 1].length ? { species: Object.fromEntries(walls[q - 1].map((n) => [n, 'blocked'])) } : {});
      stack.push(layer('water', lengths[q], { name: regionNames[q], c0: { ...comps[q] }, ...(q === stirred ? { mixing: stirD } : {}) }));
    }
    if (right) stack.push(faceOf(right, VR), layer(right.metal, 1e-6), ohmic(rightDrive, ['e-']));
    else if (rightBath) stack.push(bath({ ...comps[nComp - 1] }, refArg(rightBath), rightDrive));
    else if (host) {
      stack.push({ species: { 'Li+': 'blocked' }, reactions: [{ equation: 'Li+(water) = Li+(host)', k0: host.k0, alpha: 0.5 }] });
      stack.push(layer('host', host.L, { name: 'host', c0: { 'Li+': host.x0 * host.cMax, 'e-': host.x0 * host.cMax }, grid: { hmin: host.L / 100, hmax: host.L / 15, ratio: 1.15 } }));
      stack.push(ohmic(rightDrive, ['e-']));
    } else stack.push({ phi: 'neutral' });

    // ports, now that the voltages are known
    let needsElectrons = false;
    const defPorts = ports.map((p) => {
      const port = p.port;
      if (p.type === 'porous') {
        const V = p.level - shift;
        if (p.drive === 'V') port.V = V;
        else if (p.drive === 'R') (port.V = V), (port.R = r.log(-4, -1.5));
        else port.I = 0;
        if (p.cap) port.capacitance = { C: p.cap.C, zeroCharge: V - phi(p.q) + p.cap.offset };
        if (p.ads) {
          port.surface = { Xads: { mu0: mu(p.q, p.ads.X) + FARADAY * V - RT * Math.log(p.ads.theta / (1 - p.ads.theta)), capacity: p.ads.capacity, theta0: p.ads.theta } };
          port.reactions.push({ equation: `${p.ads.X} = Xads + e-`, k0: r.log(-8, -5), alpha: 0.5 });
        }
        needsElectrons = true;
      } else if (p.type === 'ref') needsElectrons ||= port.terminal === 'e-';
      else if (p.type === 'cap') {
        port.V = phi(p.q) + p.Vrel;
        port.capacitance = { C: p.cap.C, zeroCharge: p.Vrel + p.cap.offset };
      } else {
        port.V = mu(p.q, p.ion) / (z(p.ion) * FARADAY);
        if (geometry === 'planar' && r.chance(0.3)) port.R = r.log(-4, -2);
        // (it feeds its ion's amount in its compartment, a chemical capacitance, through G·w and R)
        slow.push(((FARADAY * FARADAY * comps[p.q][p.ion] * Ltot) / RT) * (1 / (port.species[p.ion].G * (port.to - port.from)) + (port.R ?? 0)));
      }
      return port;
    });

    // --- bulk reactions, and the grid: double layers resolved where ε > 0, reaction layers
    // (√(D/k), a quarter of it per cell) where something reacts in the bulk
    const bulkReactions = [];
    let reactionLayer = Infinity;
    if (complex) {
      const kf = r.log(0, 4), K = 10 ** complex.logK / 1000;
      bulkReactions.push({ equation: `${complex.M} + ${complex.L} = ${complex.name}`, kf: { water: kf } });
      for (const c of comps) {
        const k = kf * (c[complex.M] + c[complex.L]) + kf / K;
        reactionLayer = Math.min(reactionLayer, Math.sqrt(1e-9 / k));
        slow.push(1 / k);
      }
    }
    if (buffer) {
      const kf = r.log(0, 4);
      bulkReactions.push({ equation: 'HAc = H+ + Ac-', kf: { water: kf } });
      for (const c of comps) {
        const kb = (kf * c.HAc) / (c['H+'] * c['Ac-']), k = kf + kb * (c['H+'] + c['Ac-']);
        reactionLayer = Math.min(reactionLayer, Math.sqrt(1e-9 / k));
        slow.push(1 / k);
      }
    }
    const debye = Math.sqrt((78.3 * EPS0 * RT) / (2 * FARADAY * FARADAY * Math.max(...comps.map(ionic))));
    let hmin = neutral ? Ltot / r.u(200, 1000) : debye / r.u(3, 8);
    hmin = Math.min(hmin, reactionLayer / 4);
    if (geometry !== 'planar') hmin = Math.min(hmin, Ltot / 500);
    const hmax = Math.max(hmin * 2, Math.min(...lengths) / r.u(15, 40));
    const grid = { hmin, hmax, ratio: r.u(1.08, 1.2) };

    const library = [lib, ...[...metals].map((m) => metal(m))];
    if ((needsElectrons || host) && !metals.size) library.push({ species: [{ name: 'e-', z: -1 }], materials: {} });
    const r0 = geometry === 'planar' ? 0 : r.log(-6, -4.7);
    const def = build({
      T,
      library,
      ...(host ? { materials: { host: host.mat } } : {}),
      stack,
      bulkReactions,
      ports: defPorts,
      grid,
      ...(geometry === 'planar' ? {} : { geometry: { type: geometry, r0 } }),
    });
    // (with a cross-section, I and R are totals: scaled by the electrode's area)
    if (geometry !== 'planar') {
      const A = geometry === 'spherical' ? 4 * Math.PI * r0 * r0 : 2 * Math.PI * r0;
      const c = def.contacts.right;
      if (c.I !== undefined) c.I *= A;
      if (c.R !== undefined) c.R /= A;
    }

    // --- the plan
    slow.push((Ltot * Ltot) / Math.min(...names.map(D)));
    // A step moves the electrolyte's composition: its chemical capacitance, F² c L / RT for its most
    // dilute electroactive ion, charges through the electrodes' charge-transfer resistances
    // RT / (F² k₀ √(c/c_ref)) in series with the solution's (in a closed cell slower than diffusion
    // by far, where k₀ is small).
    // An ion's store counts what holds it in reserve (its complex, HAc for H⁺). Where faces block
    // some species, the solution's resistance is that of the ions that cross them all.
    const rct = (k0, kind, q, scale = 1) => RT / (FARADAY * FARADAY * k0 * scale * Math.sqrt(Math.min(...COUPLES[kind].ions.map((n) => comps[q][n])) / 1000));
    const store = (q, n) => comps[q][n] + (complex?.M === n || complex?.L === n ? comps[q][complex.name] : 0) + (buffer && (n === 'H+' || n === 'Ac-') ? comps[q].HAc : 0);
    const through = names.filter((n) => !walls.some((b) => b.includes(n)));
    const solRthrough = lengths.reduce((t, Lq, q) => t + Lq / (((FARADAY * FARADAY) / RT) * through.reduce((s, n) => s + z(n) * z(n) * D(n) * comps[q][n], 0)), 0);
    let Rtot = (blockingWalls && Number.isFinite(solRthrough) && solRthrough > 0 ? solRthrough : solR) + (rightDrive.R ?? 0), cAct = Infinity;
    const stores = (kind, q) => COUPLES[kind].ions.map((n) => store(q, n));
    for (const e of [left, right]) if (e?.rx || e?.vk0) (Rtot += rct(e.rx?.k0 ?? e.vk0, e.kind, e.q)), (cAct = Math.min(cAct, ...stores(e.kind, e.q)));
    for (const p of ports) if (p.type === 'porous' && p.drive !== 'float') (Rtot += rct(p.port.reactions[0].k0, p.kind, p.q, p.port.area * (p.port.to - p.port.from))), (cAct = Math.min(cAct, ...stores(p.kind, p.q)));
    if (Number.isFinite(cAct)) slow.push(((FARADAY * FARADAY * cAct * Ltot) / RT) * Rtot);
    // (a second couple on an electrode turns the composition over at its own, slower, rate)
    for (const e of [left, right]) if (e?.second) slow.push(((FARADAY * FARADAY * Math.min(...stores(e.second.couple.kind, e.q)) * Ltot) / RT) * rct(e.second.rx.k0, e.second.couple.kind, e.q));
    // the terminal to drive: the right end if it's held and does something, else a held port, else the left
    let terminal = 'right', V0 = VR;
    if (rightDrive.I !== undefined || rightEnd === 'wall' || (right?.kind === 'blocking' && r.chance(0.5))) {
      const held = defPorts.find((p) => p.V !== undefined);
      if (held) (terminal = held.name), (V0 = held.V);
      else (terminal = 'left'), (V0 = VL);
    }
    const anyBath = Boolean(leftBath || rightBath);
    const floating = rightDrive.I !== undefined || defPorts.some((p) => p.I !== undefined);
    // At rest, something holds charge: a Stern or diffuse layer off its zero charge, a port's
    // capacitance. In a closed stretch that takes ions from its own amounts, and moves the rest.
    const charging = [left, right].some((e) => e && e.law !== 'neutral' && e.zcOffset !== 0) || ports.some((p) => p.cap && p.cap.offset !== 0);
    // At rest: no current anywhere, every phase in equilibrium within itself (a cell between unlike
    // electrodes at its open-circuit voltage is at rest). At equilibrium, besides: every species'
    // level flat through the whole device, the electrons' in every metal and host included.
    const rest = exact && !mixed && rightDrive.I === undefined && !((!anyBath || blockingWalls) && charging);
    const eLevels = [left && VL, right && rightDrive.V, host && VR].filter((v) => v !== null && v !== undefined && v !== false);
    const flat = names.every((n) => comps.every((_, q) => Math.abs(mubar(q, n) - mubar(0, n)) <= 1e-9 * RT)) && eLevels.every((v) => Math.abs(v - eLevels[0]) <= 1e-9);
    const equilibrium = rest && flat;
    const unique = anyBath && !blockingWalls && !floating && !closedSurface;
    return {
      def,
      params: {
        features,
        leftEnd, rightEnd,
        left: left && { kind: left.kind, metal: left.metal, law: left.law, surface: left.surface, second: left.second?.couple.kind },
        right: right && { kind: right.kind, metal: right.metal, law: right.law, surface: right.surface, second: right.second?.couple.kind },
        epsr, names, comps, lengths, walls, stats: stats?.type ?? 'ideal', geometry, r0,
        ports: ports.map((p) => p.type + (p.kind ? `:${p.kind}` : '') + (p.drive ? `:${p.drive}` : '')),
        VL, VR, rightDrive, exact,
      },
      plan: { terminal, V: V0 + r.u(-0.12, 0.12), dV: r.u(0.005, 0.02), equilibrium, rest, unique, slow: Math.max(...slow.filter(Number.isFinite)) },
    };
  },
};
