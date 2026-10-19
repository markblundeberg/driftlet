import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, trap } from '../src/kit.js';

// SRH (trap-assisted) kinetics as a rate law: r = n p (1 − e^{−A/RT}) / (τ_p (n + n₁) + τ_n (p + p₁)),
// with p₁ = n p e^{−A/RT}/n₁ from the state, so equilibrium (A = 0) stays exact.

const T = 300, RT = GAS_CONSTANT * T;
const Si = semiconductor('Si');
const { 'e-': el, 'h+': ho } = Si.materials.Si.species;
const ni2 = el.cRef * ho.cRef * Math.exp(-(ho.mu0 + el.mu0) / RT);

// The root of f on [lo, hi] by bisection (f changes sign there).
const root = (f, lo, hi) => {
  for (let k = 0; k < 200; k++) {
    const m = 0.5 * (lo + hi);
    if (f(lo) * f(m) <= 0) hi = m;
    else lo = m;
  }
  return 0.5 * (lo + hi);
};

test('bulk SRH under uniform light: G = R_SRH(n, p) with n = N_D + p, for a midgap trap and a shallow one, low injection to high', () => {
  // An n-type slab lit uniformly, its holes blocked at both contacts: no gradients, so at every
  // node the generation G equals the SRH rate, and neutrality gives n = N_D + p.
  const ND = units.perCm3(1e15), tauN = 2e-6, tauP = 1e-6;
  for (const n1 of [undefined, units.perCm3(1e13)]) {
    for (const G of [units.perCm3(1e19), units.perCm3(1e22)]) {
      const s = new Device(
        build({
          T,
          library: [Si],
          stack: [ohmic(0, ['e-']), layer('Si', 10e-6, { donors: ND, c0: { 'h+': 1e-12 } }), ohmic(0, ['e-'])],
          bulkReactions: [
            { equation: 'e- + h+ = 0', srh: { Si: { tauN, tauP, ...(n1 === undefined ? {} : { n1 }) } } },
            { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: G } },
          ],
        }),
      ).solve();
      assert.ok(s.converged);
      const k1 = n1 ?? Math.sqrt(ni2), p1 = ni2 / k1;
      const R = (p) => {
        const n = ND + p;
        return (n * p - ni2) / (tauP * (n + k1) + tauN * (p + p1));
      };
      const p = root((q) => R(q) - G, 0, 1e3 * (ND + G * tauP));
      const mid = s.x.length >> 1;
      assert.ok(Math.abs(s.c['h+'][mid] / p - 1) < 1e-6, `n1 ${n1}, G ${G}: p ${s.c['h+'][mid]} vs ${p}`);
      // The solution reports each bulk reaction's rate: SRH matching G at every node, and in all.
      const [srh, light] = s.bulkReactions;
      assert.ok(Math.abs(srh.rate[mid] / G - 1) < 1e-6);
      assert.ok(Math.abs(srh.total / light.total - 1) < 1e-9, `R ${srh.total} vs G ${light.total} in all`);
      assert.ok(Math.abs(light.total / (G * 10e-6) - 1) < 1e-12);
      assert.ok(Math.abs(s.c['e-'][mid] / (ND + p) - 1) < 1e-6);
    }
  }
});

test('SRH at a face: zero in the dark, the law at the two edge nodes under light, and every pair either collected or recombined there', () => {
  // n⁺ | p silicon whose holes can't cross the face, only recombine there with the n⁺ side's
  // electrons: e⁻(left) + h⁺(right) = 0 with SRH velocities. Light makes pairs in the p side.
  const vn = 1e3, vp = 1e3, G = 10, Lp = 5e-6; // m/s, m/s, mol/(m³·s), m
  // (The p side is its own material, lit, so that no holes are made where they could neither
  // leave nor recombine.)
  const cell = (V, light) =>
    build({
      T,
      species: Si.species,
      materials: { n: Si.materials.Si, p: Si.materials.Si },
      stack: [
        ohmic(0, ['e-']),
        layer('n', 1e-6, { donors: units.perCm3(1e17), c0: { 'h+': 1e-12 } }),
        { dipole: 0, species: { 'e-': 'equilibrium', 'h+': 'blocked' }, reactions: [{ equation: 'e-(left) + h+(right) = 0', srh: { vn, vp } }] },
        layer('p', Lp, { acceptors: units.perCm3(1e16) }),
        ohmic(V, ['h+']),
      ],
      bulkReactions: light ? [{ equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { p: G } }] : [],
    });
  const dark = new Device(cell(0, false)).solve();
  assert.ok(dark.converged);
  assert.ok(Math.abs(dark.interfaces[0].rates[0]) < 1e-20, `dark rate ${dark.interfaces[0].rates[0]}`);
  for (const V of [0, 0.4, 0.5]) {
    const s = new Device(cell(V, true)).solve();
    assert.ok(s.converged);
    const gL = s.x.findIndex((x, k) => s.x[k + 1] === x); // the face's left node (it's doubled)
    const n = s.c['e-'][gL], p = s.c['h+'][gL + 1];
    const a = (s.mu['e-'][gL] + s.mu['h+'][gL + 1]) / RT;
    const n1 = Math.sqrt(n * p * Math.exp(-a)); // midgap, against the state
    const law = (n * p * -Math.expm1(-a)) / ((n + n1) / vp + (p + n1) / vn);
    const r = s.interfaces[0].rates[0];
    assert.ok(Math.abs(r / law - 1) < 1e-9, `V ${V}: rate ${r} vs the law ${law}`);
    // Every pair made in the p side is collected (electrons out by the left) or recombines at
    // the face; this is the current. (At 0.4 V the face takes 70% of them; at 0.5 V, what the
    // bias injects outweighs the light.)
    assert.ok(Math.abs(s.current / (FARADAY * (G * Lp - r)) - 1) < 1e-9, `V ${V}: ${s.current} vs ${FARADAY * (G * Lp - r)}`);
  }
});

test('SRH definitions are checked', () => {
  const def = (rx) => build({ T, library: [Si], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1 }), ohmic(0)], bulkReactions: [rx] });
  for (const [rx, message] of [
    [{ equation: 'e- + h+ = 0', kf: { Si: 1 }, srh: { Si: { tauN: 1, tauP: 1 } } }, /kf \(mass action\) or srh/],
    [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1 } } }, /tauP must be/],
    [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1, tauP: 1, vn: 1 } } }, /vn/],
    [{ equation: 'photon = e- + h+', fixed: { photon: 1e5 }, srh: { Si: { tauN: 1, tauP: 1 } } }, /consuming one negative and one positive species/],
  ]) {
    assert.throws(() => new Device(def(rx)), (e) => e instanceof DeviceError && message.test(e.message), JSON.stringify(rx));
  }
});

test('explicit immobile traps in steady state recombine exactly at the SRH rate, and keep their total node by node', () => {
  // The same lit n-type slab, recombining through trap states X⁰/X⁻ that don't move (D = 0):
  // e⁻ + X⁰ = X⁻ (capture coefficient c_n) and X⁻ + h⁺ = X⁰ (c_p). In steady state a single
  // trap level recombines at exactly R_SRH(n, p) with τ_n = 1/(c_n N_t), τ_p = 1/(c_p N_t), and
  // n₁ the electron density that half-fills it. The trap total X⁰ + X⁻ is conserved at each node
  // on its own, which the steady solve now holds as a row of that node's block (no time steps).
  const ND = units.perCm3(1e15), Nt = units.perCm3(1e12), tauN = 2e-6, tauP = 1e-6, G = units.perCm3(1e20);
  const n1 = units.perCm3(1e12);
  const cn = 1 / (tauN * Nt), cp = 1 / (tauP * Nt);
  const species = [...Si.species, { name: 'X0', z: 0, cRef: Nt }, { name: 'X-', z: -1, cRef: Nt }];
  const mat = { ...Si.materials.Si, species: { ...Si.materials.Si.species, X0: { D: 0, mu0: 0 }, 'X-': { D: 0, mu0: el.mu0 + RT * Math.log(n1 / el.cRef) } } };
  const dev = new Device(
    build({
      T,
      species,
      materials: { Si: mat },
      stack: [ohmic(0, ['e-']), layer('Si', 10e-6, { donors: ND, c0: { 'h+': 1e-12, X0: Nt / 2, 'X-': Nt / 2 } }), ohmic(0, ['e-'])],
      bulkReactions: [
        { equation: 'e- + X0 = X-', kf: { Si: cn } },
        { equation: 'X- + h+ = X0', kf: { Si: cp } },
        { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: G } },
      ],
    }),
  );
  const s = dev.solve();
  assert.ok(s.converged);
  assert.ok(s.history.every((h) => !(h.dt < Infinity)), 'solved directly, with no time steps');
  const mid = s.x.length >> 1, n = s.c['e-'][mid], p = s.c['h+'][mid];
  const R = (n * p - ni2) / (tauP * (n + n1) + tauN * (p + ni2 / n1));
  assert.ok(Math.abs(R / G - 1) < 1e-6, `R_SRH ${R} against G ${G}`);
  s.c.X0.forEach((x, g) => assert.ok(Math.abs((x + s.c['X-'][g]) / Nt - 1) < 1e-12, `trap total at node ${g}`));
});

test('trap(): the same trap written for you recombines at the srh law\'s rate, and its charge compensates the doping, n = N_D − N_t n/(n + n₁)', () => {
  // Lit, as above: in steady state the explicit trap recombines at R_SRH(n, p) exactly.
  const ND = units.perCm3(1e15), Nt = units.perCm3(1e12), tauN = 2e-6, tauP = 1e-6, G = units.perCm3(1e20), n1 = units.perCm3(1e12);
  const lit = trap(Si, { density: Nt, tauN, tauP, n1 });
  const s = new Device(build({
    T,
    library: [lit.library],
    stack: [ohmic(0, ['e-']), layer('Si', 10e-6, { donors: ND, c0: { 'h+': 1e-12, ...lit.c0 } }), ohmic(0, ['e-'])],
    bulkReactions: [...lit.bulkReactions, { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: G } }],
  })).solve();
  assert.ok(s.converged);
  const mid = s.x.length >> 1, n = s.c['e-'][mid], p = s.c['h+'][mid];
  const R = (n * p - ni2) / (tauP * (n + n1) + tauN * (p + ni2 / n1));
  assert.ok(Math.abs(R / G - 1) < 1e-6, `R_SRH ${R} against G ${G}`);
  // In the dark, deep traps as dense as the donors take a share of their electrons: neutrality
  // n + X⁻ = N_D (holes negligible), with the trap filled as f = n/(n + n₁). Acceptor-like
  // traps on n-type, donor-like on p-type (p + X⁺ = N_A, X⁺ the empty share).
  for (const charge of ['acceptor', 'donor']) {
    const N = units.perCm3(1e16), Ntd = units.perCm3(5e15), n1d = units.perCm3(1e15);
    const t = trap(Si, { density: Ntd, tauN: 1e-6, tauP: 1e-6, n1: n1d, charge });
    const dop = charge === 'acceptor' ? { donors: N } : { acceptors: N };
    const hold = charge === 'acceptor' ? ['e-'] : ['h+'];
    const sol = new Device(build({ T, library: [t.library], stack: [ohmic(0, hold), layer('Si', 2e-6, { ...dop, c0: { ...t.c0, [charge === 'acceptor' ? 'h+' : 'e-']: 1e-12 } }), ohmic(0, hold)], bulkReactions: t.bulkReactions })).solve();
    assert.ok(sol.converged, charge);
    const g = sol.x.length >> 1, ne = sol.c['e-'][g], ph = sol.c['h+'][g];
    if (charge === 'acceptor') {
      const want = root((x) => x + (Ntd * x) / (x + n1d) - N, 0, N);
      assert.ok(Math.abs(ne / want - 1) < 1e-6 && Math.abs(sol.c['X-'][g] / (Ntd * ne / (ne + n1d)) - 1) < 1e-6, `acceptor: n ${ne} vs ${want}`);
    } else {
      // p + X⁺ = N_A, with X⁺ = N_t (1 − f) and f = n/(n + n₁), n = n_i²/p
      const want = root((x) => x + Ntd * (1 - ni2 / x / (ni2 / x + n1d)) - N, 1e-30, N);
      assert.ok(Math.abs(ph / want - 1) < 1e-6, `donor: p ${ph} vs ${want}`);
    }
  }
});

// Gemini 4 Argon's review (its finding 4a). Under degeneracy, detailed balance scales each
// carrier's emission by its own γ = F₁/₂(ζ)/e^ζ = e^{−ex}: n₁ = γ_n n₁₀, p₁ = γ_p p₁₀. The product
// was always right (n₁p₁ = n p e^{−a}); the split put a degenerate n⁺ region's γ_n on its holes.
// In a device it hardly shows (where a carrier is degenerate, n ≫ n₁), so the rate is checked
// against the formula directly.
test('SRH under degeneracy: each carrier\'s trap density scales by its own γ, and the rate\'s derivatives match differences', () => {
  const dev = new Device(build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donors: units.perCm3(1e17) }), ohmic(0)], bulkReactions: [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1e-6, tauP: 2e-6 } } }] }));
  const S = dev.solver, tn = 1e-6, tp = 2e-6, ni0 = 1e-8, gn = 0.0349, gp = 1; // (mol/m³)
  const xn = -Math.log(gn), xp = -Math.log(gp);
  const cn = 5000, cp = 3e-15, a = 2.3; // so n p e^{−a} = γ_n γ_p n_i0²
  const shifted = (cn * cp * Math.exp(-a)) / (gn * gp);
  const rate = (law, nOk) => {
    const r = S._srhRate({ tn, tp, n1: law }, cn, cp, a, xn, xp).rate;
    const nPair = Math.sqrt(shifted), n1 = nOk ? gn * law : gn * nPair, p1 = nOk ? (gp * shifted) / law : gp * nPair;
    const want = (cn * cp * -Math.expm1(-a)) / (tp * (cn + n1) + tn * (cp + p1));
    return r / want - 1;
  };
  assert.ok(Math.abs(rate(NaN, false)) < 1e-12); // midgap: n₁ = γ_n n_i, p₁ = γ_p n_i
  assert.ok(Math.abs(rate(ni0 * 30, true)) < 1e-12); // a given n1: γ_n n1, and p₁ = γ_p n_i²/n1
  // Derivatives (Boltzmann carriers, s = 1): ∂/∂ln n, ∂/∂ln p, ∂/∂a against central differences.
  for (const n1 of [NaN, 1e-7]) {
    const law = { tn, tp, n1 }, at = (ln, lp, aa) => S._srhRate(law, Math.exp(ln), Math.exp(lp), aa).rate;
    const r0 = S._srhRate(law, 2e-3, 4e-6, -0.7), ln = Math.log(2e-3), lp = Math.log(4e-6), h = 1e-5;
    assert.ok(Math.abs((at(ln + h, lp, -0.7) - at(ln - h, lp, -0.7)) / (2 * h) / r0.Cn - 1) < 1e-5);
    assert.ok(Math.abs((at(ln, lp + h, -0.7) - at(ln, lp - h, -0.7)) / (2 * h) / r0.Cp - 1) < 1e-5);
    assert.ok(Math.abs((at(ln, lp, -0.7 + h) - at(ln, lp, -0.7 - h)) / (2 * h) / r0.Ca - 1) < 1e-5);
  }
  // Extraction far past e^709 (a quasi-Fermi split of −20 V): finite, the generation limit.
  const g = S._srhRate({ tn, tp, n1: NaN }, 1e-300, 1e-300, -800);
  assert.ok(Number.isFinite(g.rate) && Number.isFinite(g.Ca) && g.rate < 0);
});

test('SRH whose carriers share a statistics model (lattice, Debye–Hückel): its rows\' derivatives match differences, given n1 or midgap', () => {
  // Each carrier's excess potential depends on the other's concentration there: taken as a
  // function of its own alone (through c/K_nn), the e⁻ row's ∂/∂η_h was 8% off on a lattice and
  // 20% under Debye–Hückel, and a cold solve took 20 Newton iterations, not 12 (found by a cold
  // review of src/regions.js). Checked on the SRH term alone: in the full rows the fluxes dwarf it.
  for (const stats of [{ type: 'lattice', species: ['e-', 'h+'], cMax: 200 }, { type: 'debye-huckel', species: ['e-', 'h+'], epsr: 11.7, a: 1e-9 }]) {
    for (const n1 of [units.perCm3(3e19), undefined]) {
      const Si = semiconductor('Si');
      Si.materials.Si.statistics = [stats];
      const d = new Device(build({
        T: 300,
        library: [Si],
        stack: [ohmic(0, ['e-']), layer('Si', units.um(1), { donors: units.perCm3(1e19), c0: { 'h+': 1 } }), ohmic(0, ['e-'])],
        bulkReactions: [
          { equation: 'e- + h+ = 0', srh: { Si: { tauN: 1e-7, tauP: 2e-7, ...(n1 ? { n1 } : {}) } } },
          { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: units.perCm3(1e26) } },
        ],
        grid: { hmin: units.nm(20), hmax: units.nm(50), ratio: 1.3 },
      }));
      const sol = d.solve(), s = d.solver;
      assert.ok(sol.converged && sol.iterations <= 14, `${stats.type}, n1 ${n1}: ${sol.iterations} iterations`);
      const g = s.nNodes >> 1, b = s.blockOfNode[g], M = s.M, rxs = s.rxsIn[0].filter((rx) => rx.srh);
      const kernel = () => {
        s.sys.clear();
        s.res.fill(0);
        s.computeConcentrations();
        s._bulkReactions(g, rxs);
      };
      kernel();
      const m = s.sys.sizes[b], oB = s.sys.offB[b];
      const J = [1, 2].map((r) => [0, 1, 2].map((col) => s.sys.B[oB + s.loc[b * M + r] * m + s.loc[b * M + col]]));
      for (let col = 0; col < 3; col++) {
        const o = b * M + col, u0 = s.u[o], h = 1e-6;
        s.u[o] = u0 + h;
        kernel();
        const fp = [1, 2].map((r) => s.res[s.rix[b * M + r]]);
        s.u[o] = u0 - h;
        kernel();
        const fm = [1, 2].map((r) => s.res[s.rix[b * M + r]]);
        s.u[o] = u0;
        for (let r = 0; r < 2; r++) {
          const fd = (fp[r] - fm[r]) / (2 * h);
          assert.ok(Math.abs(J[r][col] - fd) <= 1e-6 * Math.abs(fd), `${stats.type}, n1 ${n1}, row ${r + 1}, column ${col}: ${J[r][col]} vs ${fd}`);
        }
      }
    }
  }
});

test('trap()\'s midgap default is the material\'s n_i whatever its φ reference: moving μ°(e⁻) by −F and μ°(h⁺) by +F changes no current', () => {
  // n_i² = N_c N_v e^{−(μ°_e + μ°_h)/RT}, a sum that a shift of the φ reference leaves alone. With
  // the difference, a diode's SRH current moved 1.1% under a 1 V shift.
  const diode = (shift) => {
    const S = semiconductor('Si'), m = S.materials.Si.species;
    m['e-'] = { ...m['e-'], mu0: m['e-'].mu0 - shift * FARADAY };
    m['h+'] = { ...m['h+'], mu0: m['h+'].mu0 + shift * FARADAY };
    const t = trap(S, { density: units.perCm3(1e15), tauN: 1e-6, tauP: 1e-6 });
    const s = new Device(build({ T, library: [t.library], stack: [ohmic(0), layer('Si', 2e-6, { donors: units.perCm3(1e17), name: 'n', c0: t.c0 }), layer('Si', 2e-6, { acceptors: units.perCm3(1e17), name: 'p', c0: t.c0 }), ohmic(0.5)], bulkReactions: t.bulkReactions })).solve();
    assert.ok(s.converged);
    return s.current;
  };
  const [I0, I1] = [diode(0), diode(1)];
  assert.ok(Math.abs(I1 / I0 - 1) < 1e-6, `${I0} vs ${I1} A/m²`);
});
