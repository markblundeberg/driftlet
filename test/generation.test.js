import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, photogeneration } from '../src/kit.js';

// Generation profiles: a rate constant tabulated against x, here light absorbed as it goes in
// (Beer–Lambert), against the short-circuit current that collection theory predicts for it.

test('photogeneration(): a Beer–Lambert table that generates exactly what the light loses', () => {
  for (const [alpha, to] of [[1e7, 60e-6], [1e5, 60e-6], [1e3, 60e-6], [1e5, 1e-6]]) {
    const rx = photogeneration({ material: 'Si', flux: 2e-4, alpha, mu: units.eV(3), from: 1e-6, to: 1e-6 + to });
    assert.equal(rx.equation, 'photon = e- + h+');
    const { x, values } = rx.kf.Si;
    let total = 0;
    for (let k = 1; k < x.length; k++) total += ((values[k] + values[k - 1]) / 2) * (x[k] - x[k - 1]);
    const absorbed = 2e-4 * -Math.expm1(-alpha * to);
    assert.ok(Math.abs(total / absorbed - 1) < 1e-6, `α = ${alpha}: ${total / absorbed}`);
    assert.ok(x[0] < 1e-6 && x[1] === 1e-6 && values[0] === 0, 'nothing before from');
    assert.ok(Math.abs(x.at(-2) - (1e-6 + to)) < 1e-18 && x.at(-1) > x.at(-2) && values.at(-1) === 0, 'nor after to');
  }
  // αL of exactly 50: the last point isn't repeated by round-off (x stays increasing).
  const { x } = photogeneration({ material: 'Si', flux: 1, alpha: 1e7, mu: 1e5, from: 0.5e-6, to: 5.5e-6 }).kf.Si;
  assert.ok(x.every((v, k) => k === 0 || v > x[k - 1]));
  assert.throws(() => photogeneration({ material: 'Si', flux: 1, alpha: 1e5, to: 1e-6 }), /mu must be the photons' μ/);
  assert.throws(() => photogeneration({ material: 'Si', flux: 1, alpha: 0, mu: 1e5, to: 1e-6 }), /absorption coefficient/);
  assert.throws(() => photogeneration({ material: 'Si', flux: 1, alpha: 1e5, mu: 1e5, from: 2e-6, to: 1e-6 }), /to > from/);
});

test('kf profiles: checked as data', () => {
  const def = (kf) => build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1 }), ohmic(0)], bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: kf } }] });
  assert.throws(() => new Device(def({ x: [0, 1e-6], values: [1, -1] })), (e) => e instanceof DeviceError && /values\[1\] must be a rate constant ≥ 0/.test(e.message));
  assert.throws(() => new Device(def({ x: [0, 1e-6], k: [1, 1] })), (e) => e instanceof DeviceError && /\bk\b/.test(e.message));
  assert.ok(new Device(def({ x: [0, 1e-6], values: [1e8, 2e8] })).solve().converged);
});

test('Beer–Lambert in an n⁺p cell: J_sc is qΦ∫αe^{−αx}η(x)dx, from blue light lost in the emitter to red collected deep in the base', () => {
  // A 0.5 µm n⁺ emitter on a 60 µm p base, ohmic (sinks) at both ends, band-to-band
  // recombination: L_p = 1.5 µm in the emitter, L_n = 88 µm in the base. A carrier generated at x
  // is collected with probability η(x): 1 in the depletion layer, sinh(x/L_p)/sinh(x_e/L_p)
  // across the emitter to its edge x_e, sinh((W − x)/L_n)/sinh((W − x_b)/L_n) across the base
  // from its edge x_b. (Gärtner's formula is the semi-infinite base without an emitter.) The
  // grid is the same coarse one for every colour (cells up to 1 µm, where blue light dies within
  // 0.1 µm): the profile is averaged over each node's box, so what's generated is exact.
  const Si = semiconductor('Si'), We = units.um(0.5), Wb = units.um(60), W = We + Wb;
  const ND = units.perCm3(1e19), NA = units.perCm3(1e16), kf = 1 / (2e-6 * NA), flux = 1e-4;
  const { 'e-': el, 'h+': ho } = Si.materials.Si.species, VT = (GAS_CONSTANT * 300) / FARADAY;
  const ni2 = el.cRef * ho.cRef * Math.exp(-(ho.mu0 - el.mu0) / (GAS_CONSTANT * 300));
  const Wd = Math.sqrt(((2 * Si.materials.Si.epsr * EPS0 * VT * Math.log((ND * NA) / ni2)) / FARADAY) * (1 / NA + 1 / ND));
  const xe = We - (Wd * NA) / (NA + ND), xb = We + (Wd * ND) / (NA + ND);
  const Ln = Math.sqrt(el.D / (kf * NA)), Lp = Math.sqrt(ho.D / (kf * ND));
  const eta = (x) => (x < xe ? Math.sinh(x / Lp) / Math.sinh(xe / Lp) : x <= xb ? 1 : Math.sinh((W - x) / Ln) / Math.sinh((W - xb) / Ln));
  for (const [alpha, tol] of [[1e7, 3e-3], [1e6, 5e-4], [1e5, 5e-4], [1e4, 5e-4], [1e3, 5e-4]]) {
    const s = new Device(
      build({
        T: 300,
        library: [Si],
        stack: [ohmic(0), layer('Si', We, { donors: ND }), layer('Si', Wb, { acceptors: NA }), ohmic(0)],
        bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: kf } }, photogeneration({ material: 'Si', flux, alpha, mu: units.eV(3), to: W })],
        grid: { hmin: units.nm(1), hmax: units.um(1) },
      }),
    ).solve();
    assert.ok(s.converged);
    // The integral by the midpoint rule, finely enough for e^{−αx} at the bluest.
    let collected = 0;
    const N = 200000;
    for (let k = 0; k < N; k++) {
      const x = ((k + 0.5) * W) / N;
      collected += alpha * Math.exp(-alpha * x) * eta(x) * (W / N);
    }
    const ratio = s.current / (FARADAY * flux * collected);
    assert.ok(Math.abs(ratio - 1) < tol, `α = ${alpha}/m: J_sc ${ratio} of theory`);
  }
});

test('photogeneration() in one layer of a material that runs through the device: nothing outside [from, to]', () => {
  // A p⁺-i-n⁺ diode all of one silicon, light absorbed only in the i layer, under reverse bias:
  // every pair is collected, J = qΦ(1 − e^{−αL}), none from the doped layers either side.
  const Si = semiconductor('Si'), Wd = units.um(0.5), L = units.um(5), N = units.perCm3(1e19), flux = 1e-3;
  for (const alpha of [1e5, 1e6]) {
    const s = new Device(
      build({
        T: 300,
        library: [Si],
        stack: [ohmic(0), layer('Si', Wd, { acceptors: N }), layer('Si', L), layer('Si', Wd, { donors: N }), ohmic(5)],
        bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e6 } }, photogeneration({ material: 'Si', flux, alpha, mu: units.eV(3), from: Wd, to: Wd + L })],
        grid: { hmin: units.nm(1), hmax: units.nm(50) },
      }),
    ).solve();
    const ratio = -s.current / (FARADAY * flux * -Math.expm1(-alpha * L));
    assert.ok(Math.abs(ratio - 1) < 1e-3, `α = ${alpha}/m: J ${ratio} of theory`);
  }
});
