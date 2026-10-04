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
  const peak = dissolving.indexOf(Math.max(...dissolving));
  assert.ok(peak > 0 && peak < Vs.length - 1, `peak inside: ${dissolving}`);
  assert.ok(dissolving.at(-1) < dissolving[peak] / 10, `passive: ${dissolving}`);
});

test('surfaces are checked', () => {
  const r = (s) => [{ equation: 'OH- = OHads + e-', k0: 1, alpha: 0.5 }].map((x) => ({ ...x, ...s }));
  assert.throws(() => cell(0, { surface: { OHads: { mu0: 0, capacity: G } }, reactions: [{ equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 1, alpha: 0.5 }] }), (e) => e instanceof DeviceError && /no reaction of the port makes or uses it/.test(e.message));
  assert.throws(() => cell(0, { surface: { 'Na+': { mu0: 0, capacity: G } }, reactions: r({}) }), /a species of the device/);
  assert.throws(() => cell(0, { surface: { OHads: { mu0: 0, capacity: G, theta0: 0.6 }, X: { mu0: 0, capacity: G, theta0: 0.6 } }, reactions: r({}) }), /sum to less than 1/);
  assert.throws(() => cell(0, { surface: { OHads: { mu0: 0, capacity: G } }, reactions: r({ bare: 'yes' }) }), /bare must be true or false/);
});
