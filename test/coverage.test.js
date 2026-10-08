import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer, bath, aqueous, polarization } from '../src/kit.js';

// An electrode port's surface: species on its sites, with coverages θ (Langmuir: one site each,
// μ = μ° + RT ln(θ/θ₀), θ₀ = 1 − Σθ the bare fraction) and Γ mol of sites per m² of electrode. A
// reaction marked `bare` runs only on the bare fraction: a passive film's blocking.

const T = 298.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;
const G = 1e-5, L = 10e-6, area = 1e4; // sites (mol/m²); a film 10 µm thick on 100 µm of electrode per … (area per volume, m²/m³)

// NaOH (with a little Fe²⁺) beside a bath that holds OH⁻ at V = 0 (so μ̄_OH⁻ = 0 there), on an
// electrode spread through the film, held at V.
const cell = (V, { reactions, surface, I }) => {
  const lib = aqueous(['Na+', 'OH-', 'Fe2+'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  return new Device(build({
    T,
    library: [lib],
    stack: [{}, layer('water', L, { name: 'film' }), bath({ 'Na+': 100, 'OH-': 102, 'Fe2+': 1 }, 'OH-')],
    ports: [{ name: 'metal', region: 'film', ...(I === undefined ? { V } : { I }), terminal: 'e-', area, surface, reactions }],
    grid: { hmin: L / 20, hmax: L / 20 },
  }));
};
const adsorb = { equation: 'OH- = OHads + e-', k0: 1e-4, alpha: 0.5, bare: true };

test('OH⁻ adsorbing on an electrode: the Langmuir isotherm θ/(1 − θ) = e^((μ̄_OH⁻ + FV − μ°)/RT)', () => {
  for (const V of [-0.1, -0.03, 0, 0.04, 0.1]) {
    const sol = cell(V, { surface: { OHads: { mu0: 0, capacity: G } }, reactions: [adsorb] }).solve();
    assert.ok(sol.converged);
    const p = sol.ports[0], g = sol.x.indexOf(p.x[3]), theta = p.coverage.OHads[3];
    const want = 1 / (1 + Math.exp(-(sol.mu['OH-'][g] + FARADAY * V) / RT));
    assert.ok(Math.abs(theta / want - 1) < 1e-9, `V = ${V}: θ ${theta} vs ${want}`);
    assert.ok(Math.abs(p.current) < 1e-12, 'nothing passes once the surface is in equilibrium');
  }
});

test('a surface covered almost completely is flagged', () => {
  const sol = cell(0.1, { surface: { OHads: { mu0: -60e3, capacity: G } }, reactions: [adsorb] }).solve();
  assert.ok(sol.converged && sol.ports[0].bare[3] < 1e-10);
  assert.ok(sol.warnings.some((w) => /covered to a bare fraction/.test(w)), sol.warnings.join('; '));
  assert.ok(!cell(0.1, { surface: { OHads: { mu0: 0, capacity: G } }, reactions: [adsorb] }).solve().warnings.some((w) => /bare fraction/.test(w)));
});

test('filling the surface passes F·Γ·(electrode area)·Δθ through the electrode, and it ends on the isotherm', () => {
  const dev = cell(0.05, { surface: { OHads: { mu0: 0, capacity: G, theta0: 1e-4 } }, reactions: [adsorb] });
  // The port's current (OH⁻ leaving the water for the surface), sampled on a log clock and integrated.
  let Q = 0, tPrev = 0, iPrev = dev.solution().ports[0].current, s;
  for (let j = 0; j <= 300; j++) {
    const t = 1e-4 * 10 ** (j / 60); // to 10 s, the surface's time constant ~0.1 s
    s = dev.advance(t, { tol: 1e-6 });
    assert.ok(s.converged);
    Q += ((s.ports[0].current + iPrev) / 2) * (t - tPrev);
    [tPrev, iPrev] = [t, s.ports[0].current];
  }
  const theta = s.ports[0].coverage.OHads[5], want = 1 / (1 + Math.exp(-0.05 / VT));
  assert.ok(Math.abs(theta / want - 1) < 1e-6, `${theta} vs ${want}: filled`);
  const charge = FARADAY * G * area * L * (want - 1e-4);
  assert.ok(Math.abs(Q / charge - 1) < 1e-3, `${Q} vs ${charge} C/m²`);
});

test('the coverage is under the time step\'s error control: Langmuir filling follows its exponential', () => {
  // A slow adsorption (k0 1e-6) barely disturbs the solution, so θ relaxes exponentially to the
  // isotherm: dθ/dt = (k_f θ₀ − k_b θ)/Γ, at the rate constants of the solution beside it.
  const dev = cell(0.05, { surface: { OHads: { mu0: 0, capacity: G, theta0: 1e-4 } }, reactions: [{ ...adsorb, k0: 1e-6 }] });
  const sol0 = dev.solution(), g = sol0.x.indexOf(sol0.ports[0].x[3]);
  const A = (sol0.mu['OH-'][g] + FARADAY * 0.05) / RT, pre = 1e-6 * Math.sqrt(sol0.c['OH-'][g] / 1000);
  const kf = pre * Math.exp(A / 2), kb = pre * Math.exp(-A / 2), thInf = kf / (kf + kb), tau = G / (kf + kb);
  for (const t of [0.3, 1, 3].map((x) => x * tau)) {
    const theta = dev.advance(t, { tol: 1e-5 }).ports[0].coverage.OHads[3];
    const exact = thInf + (1e-4 - thInf) * Math.exp(-t / tau);
    assert.ok(Math.abs(theta - exact) < 2e-4 * thInf, `t = ${t / tau} τ: θ ${theta} vs ${exact}`);
  }
});

test('a reaction on bare sites, on an electrode with no surface species: every site is bare', () => {
  const rate = (bare) => cell(0.05, { reactions: [{ ...adsorb, equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, bare }] }).solve().ports[0].rates[0][3];
  const [plain, onBare] = [rate(false), rate(true)];
  assert.ok(plain !== 0 && onBare === plain, `${onBare} vs ${plain}`);
});

test('a metal dissolving on bare sites only: the active–passive curve, the film\'s Langmuir blocking exact', () => {
  // Fe²⁺ + 2e⁻ ⇌ Fe on bare metal, and Fe + 2OH⁻ ⇌ Fe(OH)₂ + 2e⁻ as a passivating film on the
  // same sites. Dissolution climbs with V (Tafel, one V_T per e-fold at α = ½), the film's bare
  // fraction falls two V_T per e-fold once it's mostly covered, so the current peaks and falls
  // (Flade). Iron's couple sits at −0.15 V, the film half covers at 0.05 V.
  const probe = cell(0, { surface: { OHads: { mu0: 0, capacity: G } }, reactions: [adsorb] }).solve();
  const muFe2 = probe.mu['Fe2+'][3], muOH = probe.mu['OH-'][3];
  const muFe = muFe2 - 2 * FARADAY * -0.15; // Fe²⁺ + 2e⁻ ⇌ Fe in equilibrium at V = −0.15
  const mu0 = muFe + 2 * muOH + 2 * FARADAY * 0.05;
  const surface = { 'Fe(OH)2': { mu0, capacity: G } };
  const iron = { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': muFe }, k0: 1e-6, alpha: 0.5, bare: true };
  const film = { equation: 'Fe(s) + 2 OH- = Fe(OH)2 + 2 e-', fixed: { 'Fe(s)': muFe }, k0: 1e-4, alpha: 0.5, bare: true };
  const Vs = [-0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25];
  const dissolving = Vs.map((V) => {
    const dev = cell(V, { surface, reactions: [iron, film] }), sol = dev.solve();
    assert.ok(sol.converged);
    const p = sol.ports[0], w = 3, g = sol.x.indexOf(p.x[w]), theta = p.coverage['Fe(OH)2'][w];
    // The film in equilibrium with the solution beside it: θ/(1 − θ) = e^((μ_Fe + 2μ̄_OH⁻ + 2FV − μ°)/RT).
    const iso = 1 / (1 + Math.exp(-(muFe + 2 * sol.mu['OH-'][g] + 2 * FARADAY * V - mu0) / RT));
    assert.ok(Math.abs(theta / iso - 1) < 1e-6, `V = ${V}: θ ${theta} vs ${iso}`);
    // Dissolution on the bare fraction: Butler–Volmer times (1 − θ).
    const a = (sol.mu['Fe2+'][g] - 2 * FARADAY * V - muFe) / RT;
    const bv = 1e-6 * Math.sqrt(sol.c['Fe2+'][g] / 1000) * (Math.exp(a / 2) - Math.exp(-a / 2));
    assert.ok(Math.abs(p.rates[0][w] / (bv * (1 - theta)) - 1) < 1e-9, `${p.rates[0][w]} vs ${bv * (1 - theta)}`);
    // polarization() agrees at the metal's own level.
    const pol = polarization(dev, sol, { port: 'metal', x: p.x[w] }, V);
    assert.ok(Math.abs(pol.reactions[0].rate[0] / p.rates[0][w] - 1) < 1e-9);
    return -2 * FARADAY * p.rates[0][w];
  });
  // The steady-state curve from one solution: polarization() with the film re-equilibrated at
  // each potential, against the solves above (the solution beside the metal barely changes).
  const dev0 = cell(0, { surface, reactions: [iron, film] }), sol0 = dev0.solve();
  const curve = polarization(dev0, sol0, { port: 'metal', x: sol0.ports[0].x[3] }, Vs, { surface: 'equilibrium' });
  curve.reactions[0].rate.forEach((r, j) => assert.ok(Math.abs(-2 * FARADAY * r / dissolving[j] - 1) < 2e-2, `V = ${Vs[j]}: ${-2 * FARADAY * r} vs ${dissolving[j]}`));
  const peak = dissolving.indexOf(Math.max(...dissolving));
  assert.ok(peak > 0 && peak < Vs.length - 1, `peak inside: ${dissolving}`);
  assert.ok(dissolving.at(-1) < dissolving[peak] / 10, `passive: ${dissolving}`);
});

test('what passes through a surface is conserved with it: the steady state of A⁺ + e⁻ = S, S + e⁻ = B⁻ is where the transient ends', () => {
  // A and B are fed only through the surface species S, so A + B + (what S holds) is fixed; the
  // steady solve has to count the surface's share, or A and B look fed from outside.
  const species = [{ name: 'A+', z: 1, cRef: 1000 }, { name: 'B-', z: -1, cRef: 1000 }, { name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }];
  const D = { D: 1e-9, mu0: 0 }, Gs = 1e-5, a = 1e6;
  const dev = () =>
    new Device({
      T,
      species,
      materials: { water: { epsr: 0, species: { 'A+': D, 'B-': D, 'Na+': D, 'Cl-': D } } },
      regions: [{ material: 'water', length: 1e-6, c0: { 'A+': 10, 'B-': 10, 'Na+': 100, 'Cl-': 100 } }],
      contacts: { left: { V: 0, terminal: 'Na+', species: { 'Na+': 'equilibrium' }, phi: 'bulk' }, right: {} },
      ports: [{ name: 'el', region: 0, V: 0.03, terminal: 'e-', area: a, surface: { S: { mu0: 2000, capacity: Gs, theta0: 1e-3 } },
        reactions: [{ equation: 'A+ + e- = S', k0: 1e-4, alpha: 0.5 }, { equation: 'S + e- = B-', k0: 1e-4, alpha: 0.5 }] }],
      grid: { hmin: 50e-9, hmax: 50e-9 },
    });
  const total = (d, sol) => {
    const v = d.model.grid.vol, p = sol.ports[0];
    let t = 0;
    for (let g = 0; g < sol.x.length; g++) t += v[g] * (sol.c['A+'][g] + sol.c['B-'][g]);
    p.x.forEach((x, w) => (t += Gs * a * v[sol.x.indexOf(x)] * p.coverage.S[w]));
    return t;
  };
  const d = dev(), start = total(d, d.solution()), steady = d.solve();
  assert.ok(steady.converged && Math.abs(total(d, steady) / start - 1) < 1e-12, `${total(d, steady)} vs ${start}`);
  const e = dev(), late = e.advance(1e4);
  assert.ok(late.converged && Math.abs(total(e, late) / start - 1) < 1e-12);
  const [cs, cl] = [steady.c['A+'][5], late.c['A+'][5]], [ts, tl] = [steady.ports[0].coverage.S[3], late.ports[0].coverage.S[3]];
  assert.ok(Math.abs(cl / cs - 1) < 1e-9 && Math.abs(tl / ts - 1) < 1e-9 && Math.abs(cs - 10) > 1, `A⁺ ${cs} vs ${cl}, θ ${ts} vs ${tl}`);
});

test('a surface window widened by set() starts its new nodes at the starting coverage', () => {
  const def = cell(0.05, { surface: { OHads: { mu0: 0, capacity: G, theta0: 1e-4 } }, reactions: [adsorb] }).def;
  def.ports[0].to = L / 2;
  const dev = new Device(def);
  const half = dev.solve(), n0 = half.ports[0].x.length;
  assert.ok(half.ports[0].coverage.OHads[0] > 0.8, 'filled where the window was');
  dev.set({ ports: { metal: { to: L } } });
  const now = dev.solution().ports[0];
  assert.ok(now.x.length > n0 && Math.abs(now.coverage.OHads.at(-1) / 1e-4 - 1) < 1e-9, `${now.coverage.OHads.at(-1)}`);
  assert.ok(now.coverage.OHads[0] > 0.8, 'kept where it was');
});

test('surfaces are checked', () => {
  const r = (s) => [{ equation: 'OH- = OHads + e-', k0: 1, alpha: 0.5 }].map((x) => ({ ...x, ...s }));
  assert.throws(() => cell(0, { surface: { OHads: { mu0: 0, capacity: G } }, reactions: [{ equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 1, alpha: 0.5 }] }), (e) => e instanceof DeviceError && /no reaction of the port makes or uses it/.test(e.message));
  assert.throws(() => cell(0, { surface: { 'Na+': { mu0: 0, capacity: G } }, reactions: r({}) }), /a species of the device/);
  assert.throws(() => cell(0, { surface: { OHads: { mu0: 0, capacity: G, theta0: 0.6 }, X: { mu0: 0, capacity: G, theta0: 0.6 } }, reactions: r({}) }), /sum to less than 1/);
  assert.throws(() => cell(0, { surface: { OHads: { mu0: 0, capacity: G } }, reactions: r({ bare: 'yes' }) }), /bare must be true or false/);
});

test('adsorbates that only turn into each other keep their total on every site: θ_A + θ_B as it started, their ratio at equilibrium', () => {
  // A(ads) + OH⁻ = B(ads) + e⁻, both partners fed (the bath, the electrode): nothing changes how
  // many sites A and B hold together, node by node, so a steady solve keeps that (directly).
  const surface = { Aads: { mu0: 0, capacity: G, theta0: 0.2 }, Bads: { mu0: 3e3, capacity: G, theta0: 0.1 } };
  const reactions = [{ equation: 'Aads + OH- = Bads + e-', k0: 1e-4, alpha: 0.5 }];
  for (const V of [-0.05, 0.02]) {
    const sol = cell(V, { surface, reactions }).solve();
    assert.ok(sol.converged && sol.steps === 1);
    const p = sol.ports[0];
    for (let j = 0; j < p.x.length; j++) {
      const A = p.coverage.Aads[j], B = p.coverage.Bads[j], g = sol.x.indexOf(p.x[j]);
      const r = Math.exp((0 - 3e3 + sol.mu['OH-'][g] + FARADAY * V) / RT);
      assert.ok(Math.abs(A + B - 0.3) < 1e-12 && Math.abs(B / A / r - 1) < 1e-9, `V = ${V}, node ${j}: θ_A ${A}, θ_B ${B}, ratio ${B / A} vs ${r}`);
    }
  }
});

test('an adsorbed couple, O + e⁻ = R + Cl⁻, in a film given its composition: cold, it starts at the bath\'s levels and solves directly to θ_O/θ_R = e^((E − E°)/V_T)', () => {
  // Every ion's amount given (c0), nothing fixes the film's φ but the electrode's reaction: the
  // start balances it against the surface's starting coverages, where it had left them out and
  // started 1.4 V off, so a cold solve took the pseudo-transient ramp or failed, and advance()
  // stalled at its first step.
  const E0 = 0.3, Lf = 5e-6;
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const def = (E) => build({
    T,
    library: [lib],
    stack: [{ species: { 'K+': 'blocked', 'Cl-': 'blocked' }, phi: 'neutral' }, layer('water', Lf, { name: 'film', c0: { 'K+': 500, 'Cl-': 500 } }), bath({ 'K+': 500, 'Cl-': 500 }, 'Cl-', 0)],
    ports: [{
      name: 'we', region: 'film', terminal: 'e-', V: E, area: 1 / Lf,
      surface: { Oads: { mu0: 0, capacity: G, theta0: 0.4 }, Rads: { mu0: -FARADAY * E0, capacity: G, theta0: 0.4 } },
      reactions: [{ equation: 'Oads + e- = Rads + Cl-', k0: 1e-2, alpha: 0.5 }],
    }],
    grid: { hmin: 0.1e-6, hmax: 0.5e-6 },
  });
  for (const E of [0.3, 0.17, 0, -0.2]) {
    const sol = new Device(def(E)).solve();
    assert.ok(sol.converged && sol.steps === 1, `E = ${E}: ${sol.steps} steps`);
    const O = sol.ports[0].coverage.Oads[2], R = sol.ports[0].coverage.Rads[2];
    assert.ok(Math.abs(O + R - 0.8) < 1e-12 && Math.abs(O / R / Math.exp((E - E0) / VT) - 1) < 1e-9, `E = ${E}: θ_O ${O}, θ_R ${R}`);
  }
  const r = new Device(def(E0)).advance(10);
  assert.ok(r.done && r.converged && r.rejected === 0, `advance: ${r.steps} steps, ${r.rejected} rejected`);
});
