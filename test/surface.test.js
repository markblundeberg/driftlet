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
