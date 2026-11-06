import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, bath, half, aqueous, metal, polarization, describe } from '../src/kit.js';
import { mhcIntegral } from '../src/marcus.js';

// A face's surface species (an electrode's, interfaces[f].surface): coverages on Langmuir sites,
// held at the node beside the face, their storage, their reactions with the metal and the solution.
// Against closed forms: a tethered couple's Nernstian coverages and its relaxation after a step
// (Butler–Volmer and Marcus–Hush–Chidsey, as in Chidsey's experiment), Laviron's voltammogram and
// adsorption impedance, and a Langmuir electrosorption isotherm.

const T = 298.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;
const G = 1e-6; // mol of sites per m²
// Ferrocene held at the end of a monolayer on gold, in 1 M HClO₄: Fc + ClO₄⁻ = FcClO₄ + e⁻.
const acid = { 'H+': 1000, 'ClO4-': 1000 }, L = 10e-6, D = { 'H+': 9.311e-9, 'ClO4-': 1.792e-9 };
const lib = {
  species: [{ name: 'H+', z: 1, cRef: 1000 }, { name: 'ClO4-', z: -1, cRef: 1000 }],
  materials: { acid: { epsr: 0, species: { 'H+': { D: D['H+'], mu0: 0 }, 'ClO4-': { D: D['ClO4-'], mu0: -8.52e3 } } } },
};
const mu0 = { Fc: 0, FcClO4: -10e3 };
const tethered = (V, kin = { alpha: 0.5 }, k0 = 1e-6) => build({
  T,
  library: [lib, metal('Au')],
  stack: [
    ohmic(V, ['e-']),
    layer('Au', 1e-6),
    {
      surface: { Fc: { mu0: mu0.Fc, capacity: G, theta0: 0.3 }, FcClO4: { mu0: mu0.FcClO4, capacity: G, theta0: 0.1 } },
      reactions: [{ equation: 'Fc + ClO4- = FcClO4 + e-', k0, ...kin }],
    },
    layer('acid', L, { c0: acid }),
    bath(acid, 'ClO4-'),
  ],
  grid: { hmin: 1e-8, hmax: 0.5e-6 },
});
const total = 0.4; // the couple's coverage, θ_R + θ_O, which it keeps
// E°′: equilibrium, μ°_R + RT ln θ_R + μ̄_ClO₄ = μ°_O + RT ln θ_O + μ̄_e, with θ_O = θ_R where μ̄_e = −F E°′
const formal = (sol, gR) => (mu0.FcClO4 - mu0.Fc - sol.mu['ClO4-'][gR]) / FARADAY;

test('a tethered couple: Nernstian coverages, its total kept, the bare sites reported', () => {
  const dev = new Device(tethered(0)), g = dev.model.grid, gL = g.regionEnd[0], gR = g.regionStart[1];
  for (const V of [-0.2, 0, 0.15, 0.4]) {
    dev.set({ contacts: { left: { V } } });
    const sol = dev.solve(), { coverage, bare } = sol.interfaces[0];
    assert.ok(sol.converged && Math.abs(sol.current) < 1e-15, `V ${V}`);
    // ln(θ_O/θ_R) = (μ°_R − μ°_O + μ̄_ClO₄ − μ̄_e)/RT
    const want = (mu0.Fc - mu0.FcClO4 + sol.mu['ClO4-'][gR] - sol.mu['e-'][gL]) / RT;
    assert.ok(Math.abs(Math.log(coverage.FcClO4 / coverage.Fc) - want) < 1e-9, `V ${V}: ${Math.log(coverage.FcClO4 / coverage.Fc)} vs ${want}`);
    assert.ok(Math.abs(coverage.Fc + coverage.FcClO4 - total) < 1e-12 && Math.abs(bare - (1 - total)) < 1e-12);
  }
});

// After a step to E, θ_O relaxes to its new equilibrium as e^{−kt}, k = (k0/Γ)(k_f + k_b) with the
// rate per site k0 (θ_R F(a°) − θ_O F(−a°)), F = e^{a/2} (Butler–Volmer, α = ½) or I(λ, a)/I(λ, 0),
// a° = F(E − E°′)/RT. Chidsey's experiment: ln k against E − E°′ bends over for Marcus–Hush–Chidsey.
test('a potential step: the coverage relaxes as e^{−kt}, k from Butler–Volmer or Marcus–Hush–Chidsey, as in Chidsey’s experiment', () => {
  const lam = units.eV(0.85) / RT, I0 = mhcIntegral(lam, 0).I;
  for (const [kin, F] of [[{ alpha: 0.5 }, (a) => Math.exp(a / 2)], [{ lambda: units.eV(0.85) }, (a) => mhcIntegral(lam, a).I / I0]]) {
    // (Butler–Volmer to 0.5 V: past it its k, 1e7 s⁻¹ by 0.9 V, outruns the acid's supply over nanoseconds)
    for (const eta of kin.alpha ? [0.1, 0.5] : [0.1, 0.5, 0.9]) {
      const dev = new Device(tethered(0, kin)), sol0 = dev.solve(), g = dev.model.grid, gR = g.regionStart[1];
      const E0 = formal(sol0, gR), V = E0 + eta;
      dev.set({ contacts: { left: { V } } });
      const a = eta / VT, kf = F(a), kb = F(-a), k = (1e-6 / G) * (kf + kb), eq = (total * kf) / (kf + kb);
      const th0 = sol0.interfaces[0].coverage.FcClO4;
      for (const t of [0.2 / k, 1 / k, 3 / k]) {
        const s = dev.advance(t, { tol: 1e-8 });
        const want = eq + (th0 - eq) * Math.exp(-k * t), th = s.interfaces[0].coverage.FcClO4;
        assert.ok(Math.abs(th - want) < 1e-4 * Math.abs(th0 - eq), `${kin.alpha ? 'BV' : 'MHC'} η ${eta}, t ${t}: ${th} vs ${want}`);
        // the current, oxidation into the solution: FΓ dθ_O/dt (through the gold, toward +x)
        const i = FARADAY * G * k * (eq - th);
        assert.ok(Math.abs(s.current / i - 1) < 1e-3, `current ${s.current} vs ${i}`);
      }
    }
  }
  // Chidsey's curvature: across 0.1 → 0.9 V, Butler–Volmer's k grows 5.7e6-fold (e^{0.4 V/V_T}, near enough),
  // Marcus–Hush–Chidsey's at λ = 0.85 eV 2500-fold
  const kOf = (F, eta) => F(eta / VT) + F(-eta / VT);
  const mhc = (a) => mhcIntegral(lam, a).I / I0, bv = (a) => Math.exp(a / 2);
  assert.ok(Math.abs(kOf(bv, 0.9) / kOf(bv, 0.1) / 5.66e6 - 1) < 0.01 && Math.abs(kOf(mhc, 0.9) / kOf(mhc, 0.1) / 2485 - 1) < 0.01);
});

// Laviron (1979): a fast surface couple swept slowly at ν, i = F Γ_c ν (F/RT) e^ξ/(1 + e^ξ)², ξ = F(E − E°′)/RT:
// a peak F² ν Γ_c/(4RT) at E°′, 3.53 RT/F (90.6 mV) wide at half height.
test('Laviron’s voltammogram: the peak current F²νΓ/4RT at E°′, 90.6 mV wide', () => {
  const dev = new Device(tethered(0, { alpha: 0.5 }, 1e-2)), sol0 = dev.solve(), E0 = formal(sol0, dev.model.grid.regionStart[1]);
  const nu = 0.01, span = 0.3, t1 = (2 * span) / nu;
  dev.set({ contacts: { left: { V: { t: [0, t1], values: [E0 - span, E0 + span] } } } });
  dev.solve(); // (from equilibrium at the start)
  const Es = [], is = [];
  for (let t = 0.5; t <= t1; t += 0.5) {
    const s = dev.advance(t, { tol: 1e-7 });
    Es.push(E0 - span + nu * t);
    is.push(s.current);
  }
  const top = is.reduce((b, v, k) => (v > is[b] ? k : b), 0), peak = (FARADAY * FARADAY * nu * G * total) / (4 * RT);
  assert.ok(Math.abs(is[top] / peak - 1) < 2e-3, `peak ${is[top]} vs ${peak}`);
  assert.ok(Math.abs(Es[top] - E0) < 0.006, `at ${Es[top] - E0} from E°′`);
  // the width at half height, by interpolation either side
  const cross = (k, dir) => {
    while (is[k] > is[top] / 2) k += dir;
    return Es[k] + ((Es[k - dir] - Es[k]) * (is[top] / 2 - is[k])) / (is[k - dir] - is[k]);
  };
  const width = cross(top, 1) - cross(top, -1);
  assert.ok(Math.abs(width / (3.5255 * VT) - 1) < 5e-3, `width ${width} vs ${3.5255 * VT}`);
});

// Laviron (1979): at equilibrium, Z = R_s + R_ct + 1/(iωC_ad), R_ct = RT/(F² r₀) with r₀ the exchange rate
// per area, C_ad = (F²Γ/RT) θ_O θ_R/(θ_O + θ_R), and R_s the acid's (Nernst–Einstein) over the layer.
test('the adsorption impedance: R_s + R_ct in series with C_ad, Laviron’s', () => {
  const dev = new Device(tethered(0.2, { alpha: 0.5 })), sol = dev.solve(), { Fc: thR, FcClO4: thO } = sol.interfaces[0].coverage;
  const r0 = 1e-6 * Math.sqrt(thR * thO * (acid['ClO4-'] / 1000)); // k0 θ_R^{1−α} θ_O^α (c/c_ref)^{1−α}, α = ½
  const Rct = RT / (FARADAY * FARADAY * r0), Cad = ((FARADAY * FARADAY * G) / RT) * ((thO * thR) / (thO + thR));
  const kappa = ((FARADAY * FARADAY) / RT) * (1000 * D['H+'] + 1000 * D['ClO4-']), Rs = L / kappa;
  const fs = [1e-4, 1e-2, 1, 1e2];
  const Z = dev.impedance(fs, { terminal: 'left' });
  fs.forEach((f, q) => {
    const w = 2 * Math.PI * f, re = Rs + Rct, im = -1 / (w * Cad);
    const err = Math.hypot(Z.Z.re[q] - re, Z.Z.im[q] - im) / Math.hypot(re, im);
    assert.ok(err < 1e-3, `f ${f}: ${Z.Z.re[q]} ${Z.Z.im[q]} vs ${re} ${im}`);
  });
});

// Cl⁻ + * = Cl(ads) + e⁻ on platinum: a Langmuir electrosorption isotherm, θ/(1 − θ) = e^x with
// x = (μ̄_Cl⁻ − μ̄_e − μ°_ads)/RT, its rate with the bare site it fills (the prefactor's θ₀).
const sorbing = (V, kin = { alpha: 0.5 }) => {
  const salt = { 'K+': 100, 'Cl-': 100 };
  return build({
    T,
    library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')],
    stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), { surface: { Clads: { mu0: 58e3, capacity: 2e-5, theta0: 0.1 } }, reactions: [{ equation: 'Cl- = Clads + e-', k0: 1e-3, ...kin }] }, layer('water', 10e-6, { c0: salt }), bath(salt, 'Cl-')],
    grid: { hmin: 1e-8, hmax: 0.5e-6 },
  });
};
test('electrosorption: a Langmuir isotherm against the metal’s level', () => {
  const dev = new Device(sorbing(0.6)), g = dev.model.grid, gL = g.regionEnd[0], gR = g.regionStart[1];
  const seen = [];
  for (const V of [0.45, 0.55, 0.6, 0.65, 0.75]) {
    dev.set({ contacts: { left: { V } } });
    const sol = dev.solve(), th = sol.interfaces[0].coverage.Clads;
    const x = (sol.mu['Cl-'][gR] - sol.mu['e-'][gL] - 58e3) / RT;
    assert.ok(Math.abs(th / (1 / (1 + Math.exp(-x))) - 1) < 1e-10 && Math.abs(sol.interfaces[0].bare - (1 - th)) < 1e-12, `V ${V}: θ ${th}`);
    seen.push(th);
  }
  assert.ok(seen[0] < 0.01 && seen.at(-1) > 0.99 && Math.abs(seen[2] - 0.5) < 0.02, `the isotherm, all of it: ${seen}`); // (the probe sees it)
});

test('polarization() at a face with a surface: the solved rates at the metal’s level, Butler–Volmer and Marcus–Hush–Chidsey', () => {
  for (const kin of [{ alpha: 0.4 }, { lambda: units.eV(0.6) }]) {
    for (const def of [tethered(0.3, kin), sorbing(0.6, kin)]) {
      // (mid-transient after a step, where the rates aren't zero)
      const dev = new Device(def);
      dev.solve();
      dev.set({ contacts: { left: { V: def.contacts.left.V + 0.15 } } });
      const sol = dev.advance(1e-4, { tol: 1e-7 }), Ve = sol.V['e-'][0];
      const r = polarization(dev, sol, { face: 0 }, Ve).reactions[0].rate[0];
      assert.ok(Math.abs(r / sol.interfaces[0].rates[0] - 1) < 1e-9, `${JSON.stringify(kin)}: ${r} vs ${sol.interfaces[0].rates[0]}`);
    }
  }
});

test('polarization() at a face with a chemical step beside the electron transfer: every rate the solver\'s, the chemical one passing no current (Argon\'s review of 1.6.1)', () => {
  // Volmer–Tafel: H⁺ + e⁻ = H(ads), 2 H(ads) = H₂. The Tafel step takes no electron from the metal;
  // polarization() once refused the whole face for it.
  const ions = {
    species: [{ name: 'H+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'H2', z: 0, cRef: 1000 }],
    materials: { el: { epsr: 0, species: { 'H+': { D: 9e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 }, H2: { D: 5e-9, mu0: 0 } } } },
  };
  const def = build({
    T, library: [ions, metal('Pt')],
    stack: [ohmic(0.05, ['e-']), layer('Pt', 1e-6),
      { surface: { 'H(ads)': { mu0: 0, capacity: 1e-5, theta0: 0.1 } }, reactions: [{ equation: 'H+ + e- = H(ads)', k0: 1e-3, alpha: 0.5 }, { equation: '2 H(ads) = H2', k0: 1e-3, alpha: 0.5 }] },
      layer('el', 10e-6), bath({ 'H+': 100, 'Cl-': 100, H2: 1 }, 'Cl-')],
  });
  const dev = new Device(def), sol = dev.solve();
  assert.ok(sol.converged);
  const curve = polarization(dev, sol, { face: 0 }, [sol.V['e-'][0], sol.V['e-'][0] + 0.1]);
  curve.reactions.forEach((r, k) => assert.ok(Math.abs(r.rate[0] / sol.interfaces[0].rates[k] - 1) < 1e-9, `${r.equation}: ${r.rate[0]} vs ${sol.interfaces[0].rates[k]}`));
  const tafel = curve.reactions[1];
  assert.ok(Number.isNaN(tafel.level) && tafel.current.every((c) => c === 0) && tafel.rate[1] === tafel.rate[0]); // (held coverages: the same rate at every V)
  assert.ok(Math.abs(curve.current[0] - curve.reactions[0].current[0]) < 1e-15);
  assert.ok(Math.abs(sol.interfaces[0].rates[0] / sol.interfaces[0].rates[1] - 2) < 1e-6); // (steady: two H made per H₂)
});

test("polarization() with the surface in equilibrium: a tethered couple keeps its total, and a species made from another surface species follows it, as steady solves at each V find them", () => {
  // Fc/FcClO₄ tethered (θ_R + θ_O = 0.4 kept, only its ratio set by V), H(ads) made from the acid, and
  // Y made only from H(ads). With nothing passing, each steady solve is the equilibrium at its V,
  // and the acid beside the gold is the same at every V, so one solution's curve must give each's
  // coverages. (Once the couple, and Y, were refused: no reaction made them from the solution alone.)
  const def = (V) => build({
    T,
    library: [lib, metal('Au')],
    stack: [
      ohmic(V, ['e-']), layer('Au', 1e-6),
      {
        surface: { Fc: { mu0: 0, capacity: G, theta0: 0.3 }, FcClO4: { mu0: -10e3, capacity: G, theta0: 0.1 }, Hads: { mu0: -5e3, capacity: G, theta0: 0.2 }, Y: { mu0: -2e3, capacity: G, theta0: 0.05 } },
        reactions: [
          { equation: 'Fc + ClO4- = FcClO4 + e-', k0: 1e-6, alpha: 0.5 },
          { equation: 'H+ + e- = Hads', k0: 1e-4, alpha: 0.5 },
          { equation: 'Hads + ClO4- = Y + e-', k0: 1e-4, alpha: 0.5 },
        ],
      },
      layer('acid', L, { c0: acid }), bath(acid, 'ClO4-'),
    ],
    grid: { hmin: 1e-8, hmax: 0.5e-6 },
  });
  const dev = new Device(def(0)), sol = dev.solve(), gL = dev.model.grid.regionEnd[0], Vs = [-0.3, 0, 0.1, 0.4];
  const curve = polarization(dev, sol, { face: 0 }, Vs.map((V) => sol.V['e-'][gL] + V), { surface: 'equilibrium' });
  Vs.forEach((V, j) => {
    const at = new Device(def(V)).solve();
    assert.ok(at.converged);
    for (const name of ['Fc', 'FcClO4', 'Hads', 'Y']) assert.ok(Math.abs(curve.coverage[name][j] / at.interfaces[0].coverage[name] - 1) < 1e-10, `V ${V}, ${name}: ${curve.coverage[name][j]} vs ${at.interfaces[0].coverage[name]}`);
    assert.ok(Math.abs(curve.bare[j] / at.interfaces[0].bare - 1) < 1e-10);
    assert.ok(Math.abs(curve.coverage.Fc[j] + curve.coverage.FcClO4[j] - 0.4) < 1e-12);
  });
});

test('a face’s surface is asked for clearly, and described', () => {
  const bad = (edit, re) => {
    const def = tethered(0);
    edit(def);
    assert.throws(() => new Device(def), (e) => e instanceof DeviceError && re.test(e.message), re);
  };
  bad((d) => (d.interfaces[0].surface.Unused = { mu0: 0, capacity: G }), /Unused: no reaction of the face makes or uses it/);
  bad((d) => (d.interfaces[0].surface['H+'] = { mu0: 0, capacity: G }), /a species of the device/);
  bad((d) => (d.interfaces[0].reactions[0] = { equation: 'Fc + ClO4- = FcClO4 + e-', vmax: 1, K: { 'ClO4-': 1 } }), /Butler–Volmer or Marcus–Hush–Chidsey/);
  // (only at a metal's face)
  const two = build({ T, library: [lib], stack: [bath(acid, 'ClO4-'), layer('acid', L, { c0: acid }), { surface: { S: { mu0: 0, capacity: G } } }, layer('acid', L, { c0: acid }), bath(acid, 'ClO4-')] });
  assert.throws(() => new Device(two), (e) => e instanceof DeviceError && /a metal on one side/.test(e.message));
  assert.match(describe(tethered(0)), /its surface: Fc, FcClO4 on 1.00e-6 mol\/m² of sites \(Langmuir\)/);
});

test("a face's surface species enter its cold start: a cold solve converges whichever reaction is listed first (round 17's fuzzer)", () => {
  // Pt with H(ads) (the Volmer step) beside Fe³⁺/Fe²⁺: the start balanced the first reaction at the face
  // without the surface species' μ° (−226 kJ/mol here), volts off, and Newton failed from it, saying the
  // Volmer step was out of equilibrium; listed the other way round, it converged
  const lib = aqueous(['K+', 'Cl-', 'H+', 'Fe3+', 'Fe2+'], { epsr: 0 });
  lib.materials.water.species['Fe2+'].mu0 = -194340;
  const c = { 'K+': 20.94, 'H+': 0.01695, 'Fe3+': 8.991, 'Fe2+': 2.5007 };
  c['Cl-'] = c['K+'] + c['H+'] + 3 * c['Fe3+'] + 2 * c['Fe2+'];
  const volmer = { equation: 'H+ + e- = Hads', k0: 1.374e-3, alpha: 0.5 }, iron = { equation: 'Fe3+ + e- = Fe2+', k0: 2.2045e-5, alpha: 0.6626 };
  const solve = (reactions, V) =>
    new Device(build({
      library: [lib, metal('Pt')],
      stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), { surface: { Hads: { mu0: -226342, capacity: 2.71e-6, theta0: 0.225 } }, reactions }, layer('water', 15.3e-6, { c0: c }), bath(c, 'Cl-')],
      grid: { hmin: 1.73e-8, hmax: 1.35e-6 },
    })).solve();
  for (const V of [0.3, 0.5048, 0.7]) {
    const a = solve([volmer, iron], V), b = solve([iron, volmer], V);
    assert.ok(a.converged && b.converged, `at ${V} V: ${a.warnings[0] ?? ''}`);
    assert.ok(Math.abs(a.terminals.left.current / b.terminals.left.current - 1) < 1e-9 && Math.abs(a.interfaces[0].coverage.Hads - b.interfaces[0].coverage.Hads) < 1e-9);
  }
  assert.ok(Math.abs(solve([volmer, iron], 0.5048).terminals.left.current + 0.0326621) < 1e-7);
  // The same electrode written on the right of its solution (both ends 1 V down), and as a floating
  // Pt between two baths: the start balanced neither (φ̂ left at 0; the floating metal's level
  // without the surface's μ°), and both failed. (Argon's review of 1.6.1, round 2.)
  const grid = { hmin: 1.73e-8, hmax: 1.35e-6 }, surface = { Hads: { mu0: -226342, capacity: 2.71e-6, theta0: 0.225 } }, low = { ...bath(c, 'Cl-'), V: -1 };
  const right = new Device(build({ library: [lib, metal('Pt')], stack: [low, layer('water', 15.3e-6, { c0: c }), { surface, reactions: [volmer, iron] }, layer('Pt', 1e-6), ohmic(0.5048 - 1, ['e-'])], grid })).solve();
  assert.ok(right.converged, right.warnings[0]);
  assert.ok(Math.abs(right.terminals.right.current + 0.0326621) < 1e-7, `${right.terminals.right.current}`);
  const floating = (reactions) => new Device(build({ library: [lib, metal('Pt')], stack: [bath(c, 'Cl-'), layer('water', 15.3e-6, { c0: c }), { surface, reactions }, layer('Pt', 1e-6), { surface, reactions }, layer('water', 15.3e-6, { c0: c }), bath(c, 'Cl-')], grid })).solve();
  const f1 = floating([volmer, iron]), f2 = floating([iron, volmer]);
  assert.ok(f1.converged && f2.converged, `${f1.warnings[0] ?? ''} ${f2.warnings[0] ?? ''}`);
  const Pt = (sol) => sol.V['e-'].find(Number.isFinite);
  assert.ok(Math.abs(Pt(f1) - Pt(f2)) < 1e-9 && Math.abs(f1.interfaces[0].coverage.Hads - f2.interfaces[0].coverage.Hads) < 1e-9);
});

test('bare: true on a face reaction that fills no site: its rate goes as the bare fraction, Butler–Volmer and Marcus–Hush–Chidsey (Argon\'s review of 1.6.1)', () => {
  // Ox⁺ + e⁻ = Red beside an adsorbate A(ads) in equilibrium with A in the solution (θ = θ₀ = ½):
  // on bare sites only, the reaction is the same one on clean metal at k0·θ₀, transport and all.
  // (Once `bare` was refused unless the reaction had surface participants of its own.)
  const ions = {
    species: ['Ox+', 'Red', 'Cl-', 'A'].map((name) => ({ name, z: name === 'Ox+' ? 1 : name === 'Cl-' ? -1 : 0, cRef: 1000 })),
    materials: { el: { epsr: 0, species: Object.fromEntries(['Ox+', 'Red', 'Cl-', 'A'].map((name) => [name, { D: 1e-9, mu0: 0 }])) } },
  };
  const cell = (kin, k0, surface) => build({
    T,
    library: [ions, metal('Au')],
    stack: [
      ohmic(-0.1, ['e-']), layer('Au', 1e-6),
      surface
        ? { surface: { 'A(ads)': { mu0: 0, capacity: G, theta0: 0.3 } }, reactions: [{ equation: 'Ox+ + e- = Red', k0, ...kin, bare: true }, { equation: 'A = A(ads)', k0: 1e-3, alpha: 0.5 }] }
        : { reactions: [{ equation: 'Ox+ + e- = Red', k0, ...kin }] },
      layer('el', 10e-6), bath({ 'Ox+': 10, 'Red': 10, 'Cl-': 10, A: 1000 }, 'Cl-'),
    ],
  });
  for (const kin of [{ alpha: 0.4 }, { lambda: units.eV(0.8) }]) {
    const a = new Device(cell(kin, 1e-5, true)).solve(), bare = a.interfaces[0].bare;
    const b = new Device(cell(kin, 1e-5 * bare, false)).solve();
    assert.ok(a.converged && b.converged);
    assert.ok(Math.abs(bare - 0.5) < 1e-6, `bare ${bare}`);
    assert.ok(Math.abs(a.current / b.current - 1) < 1e-9, `${JSON.stringify(kin)}: ${a.current} vs ${b.current}`);
  }
});

test('a face covered almost completely is warned of, and a tethered couple fed by a current is named as what it charges (Argon\'s review of 1.6.1)', () => {
  // H(ads) at μ° −60 kJ/mol held 0.5 V cathodic: bare 1e-17, which only a port's surface was
  // watched for. And a tethered couple at a port driven by a current: closed but for it, so no
  // steady state, which the warning put down to a capacitance the port doesn't have.
  const ions = (names) => ({ species: names.map(([name, z]) => ({ name, z, cRef: 1000 })), materials: { water: { epsr: 0, species: Object.fromEntries(names.map(([name]) => [name, { D: 1e-9, mu0: 0 }])) } } });
  const starved = new Device(build({
    library: [ions([['H+', 1], ['Cl-', -1]]), metal('Au')],
    stack: [{ V: -0.5, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, layer('Au', 1e-7),
      { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 }, surface: { Hads: { mu0: -60000, capacity: 1e-5, theta0: 0.99999999999 } }, reactions: [{ equation: 'H+ + e- = Hads', k0: 1e-2, alpha: 0.5 }] },
      layer('water', 1e-6, { c0: { 'H+': 100, 'Cl-': 100 } }), bath({ 'H+': 100, 'Cl-': 100 }, 'Cl-')],
  })).solve();
  assert.ok(starved.interfaces[0].bare < 1e-10);
  assert.ok(starved.warnings.some((w) => /interfaces\[0\]: its surface is covered to a bare fraction of/.test(w)), starved.warnings.join('\n'));
  const salt = { 'K+': 100, 'ClO4-': 100 };
  const fed = new Device(build({
    library: [ions([['K+', 1], ['ClO4-', -1]]), metal('Au')],
    stack: [bath(salt, 'ClO4-'), layer('water', 1e-6, { name: 'cell', c0: salt }), bath(salt, 'ClO4-')],
    ports: [{ name: 'work', region: 'cell', terminal: 'e-', area: 1e6, I: 1e-3, surface: { Fc: { mu0: 0, capacity: 1e-5, theta0: 0.45 }, FcClO4: { mu0: 0, capacity: 1e-5, theta0: 0.45 } }, reactions: [{ equation: 'Fc + ClO4- = FcClO4 + e-', k0: 1e-2, alpha: 0.5 }] }],
  })).solve();
  assert.ok(!fed.converged && fed.warnings.some((w) => /what it feeds \(Fc on ports work\) is closed but for it/.test(w)), fed.warnings.join('\n'));
});

test('a face\'s made covers its surface: the couple\'s turnover, G dθ/dt as it relaxes after a step', () => {
  const dev = new Device(tethered(0)), sol0 = dev.solve(), gR = dev.model.grid.regionStart[1];
  const eta = 0.1, V = formal(sol0, gR) + eta, a = eta / VT, kf = Math.exp(a / 2), kb = Math.exp(-a / 2), k = (1e-6 / G) * (kf + kb), eq = (total * kf) / (kf + kb);
  dev.set({ contacts: { left: { V } } });
  const s = dev.advance(1 / k, { tol: 1e-8 }), th = s.interfaces[0].coverage.FcClO4, made = s.interfaces[0].made.surface;
  const want = G * k * (eq - th); // (d/dt of G θ along e^{−kt})
  assert.ok(Math.abs(made.FcClO4 / want - 1) < 1e-3, `${made.FcClO4} against ${want}`);
  assert.ok(Math.abs(made.Fc / want + 1) < 1e-3, `${made.Fc}`);
});
