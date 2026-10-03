import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Device, units, AVOGADRO, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { photogeneration } from '../src/kit.js';

// A benchmark against an independent code: J–V hysteresis in a planar perovskite solar cell,
// computed by IonMonger (finite elements in MATLAB/Octave, Courtier et al. 2019) and by
// driftlet, scan by scan from 1 mV/s to 1 kV/s. The reference scans and how they were made are
// in test/fixtures/ionmonger.
//
// TiO₂ | MAPbI₃ with mobile iodide vacancies | spiro-OMeTAD, IonMonger's default parameters,
// bulk SRH off. Three species: electrons, holes, vacancies (over an equal immobile background).
// IonMonger's interface recombination is SRH at each face, saturating at one carrier's capture;
// written with the ETL's own electron density it's exactly driftlet's face SRH law, with n₁ =
// (d_E/n₀)·n_i at the TiO₂ face (n₀ the perovskite's electron density at the ETL's Fermi level)
// and n₁ = n_i at the spiro face. Levels are on the vacuum scale, as IonMonger's are.

const T = 298, VT = (GAS_CONSTANT * T) / FARADAY, m3 = (x) => x / AVOGADRO, eV = units.eV;
const P = {
  b: 400e-9, eps: 24.1, alpha: 1.3e7, Ec: -3.7, Ev: -5.4, Dn: 1.7e-4, Dp: 1.7e-4, gc: 8.1e24, gv: 5.8e24,
  N0: 1.6e25, DI: 6.5e-8 * Math.exp(-0.58 / (8.61733035e-5 * T)), // IonMonger's Arrhenius form, ≈ 1.0e-17 m²/s
  dE: 1e24, gcE: 5e25, EcE: -4.0, bE: 100e-9, epsE: 10, DE: 1e-5,
  dH: 1e24, gvH: 5e25, EvH: -5.1, bH: 200e-9, epsH: 3, DH: 1e-6,
  tn: 3e-9, tp: 3e-7, vnE: 1e5, vpE: 10, vnH: 0.1, vpH: 1e5, Fph: 1.4e21,
};
const ni = Math.sqrt(P.gc * P.gv) * Math.exp(-(P.Ec - P.Ev) / (2 * VT));
const n0 = P.gc * Math.exp((P.EcE + VT * Math.log(P.dE / P.gcE) - P.Ec) / VT);

const cell = (V, bulk) => ({
  T,
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
    { name: 'V+', z: 1 }, // iodide vacancies
  ],
  materials: {
    TiO2: { epsr: P.epsE, species: { 'e-': { D: P.DE, mu0: eV(P.EcE), cRef: m3(P.gcE) } } },
    MAPI: {
      epsr: P.eps,
      species: {
        'e-': { D: P.Dn, mu0: eV(P.Ec), cRef: m3(P.gc) },
        'h+': { D: P.Dp, mu0: eV(-P.Ev), cRef: m3(P.gv) },
        'V+': { D: P.DI, mu0: 0, cRef: m3(P.N0) },
      },
    },
    spiro: { epsr: P.epsH, species: { 'h+': { D: P.DH, mu0: eV(-P.EvH), cRef: m3(P.gvH) } } },
  },
  regions: [
    { name: 'TiO₂', material: 'TiO2', length: P.bE, fixedCharge: FARADAY * m3(P.dE) },
    { name: 'MAPbI₃', material: 'MAPI', length: P.b, fixedCharge: -FARADAY * m3(P.N0), c0: { 'V+': m3(P.N0) } },
    { name: 'spiro-OMeTAD', material: 'spiro', length: P.bH, fixedCharge: -FARADAY * m3(P.dH) },
  ],
  interfaces: [
    { dipole: 0, species: { 'e-': 'equilibrium' }, reactions: [{ equation: 'e-(left) + h+ = 0', srh: { vn: P.vnE, vp: P.vpE, n1: m3((P.dE / n0) * ni) } }] },
    { dipole: 0, species: { 'h+': 'equilibrium' }, reactions: [{ equation: 'e- + h+(right) = 0', srh: { vn: P.vnH, vp: P.vpH, n1: m3(ni) } }] },
  ],
  contacts: {
    left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
    right: { V, terminal: 'h+', species: { 'h+': 'equilibrium' }, phi: 'bulk' },
  },
  bulkReactions: [
    // IonMonger's bulk SRH (midgap traps) is driftlet's bulk SRH law.
    ...(bulk ? [{ equation: 'e- + h+ = 0', srh: { MAPI: { tauN: P.tn, tauP: P.tp } } }] : []),
    photogeneration({ material: 'MAPI', flux: m3(P.Fph), alpha: P.alpha, mu: eV(3), from: P.bE, to: P.bE + P.b })],
  grid: { hmin: 0.05e-9, hmax: 5e-9, ratio: 1.15 },
});

// Preconditioned at 1.2 V (steady state, light on), then 1.2 → 0 → 1.2 V: [t, V, J (mA/cm²)].
function scan(rate, bulk) {
  const dev = new Device(cell(1.2, bulk));
  assert.ok(dev.solve().converged);
  const half = 1.2 / rate;
  dev.set({ contacts: { right: { V: { t: [0, half, 2 * half], values: [1.2, 0, 1.2] } } } });
  const r = dev.advance(2 * half, { tol: 1e-4, dtMax: (2 * half) / 480 });
  assert.ok(r.converged && r.done);
  return r.trace.t.map((t, k) => [t, r.trace.voltage[k], r.trace.current[k] / 10]);
}

// Each sweep's maximum power and open-circuit voltage, and the hysteresis index.
function loop(rows, rate) {
  const half = 1.2 / rate, out = {};
  for (const [name, on] of [['rev', (t) => t <= half], ['fwd', (t) => t >= half]]) {
    const pts = rows.filter(([t, , J]) => on(t) && Number.isFinite(J)).map(([, V, J]) => [V, J]).sort((a, b) => a[0] - b[0]);
    const k = pts.findIndex(([, J], j) => j > 0 && pts[j - 1][1] > 0 && J <= 0);
    const [[v0, j0], [v1, j1]] = [pts[k - 1], pts[k]];
    out[name] = { pmax: Math.max(...pts.map(([V, J]) => V * J)), voc: v0 + (j0 * (v1 - v0)) / (j0 - j1) };
  }
  return { ...out, hi: (out.rev.pmax - out.fwd.pmax) / out.rev.pmax };
}

const reference = (set, rate) =>
  readFileSync(new URL(`fixtures/ionmonger/${set}_${rate}.csv`, import.meta.url), 'utf8').trim().split('\n').map((line) => line.split(',').map(Number));

for (const [set, bulk, what] of [['scan', false, 'bulk SRH off'], ['full', true, "IonMonger's full defaults"]]) {
  test(`perovskite hysteresis against IonMonger (${what}), 1 mV/s to 1 kV/s: hysteresis index, maximum power, V_oc and the whole J–V loop`, () => {
    for (const rate of ['0.001', '0.01', '0.1', '1', '10', '100', '1000']) {
      const ours = scan(+rate, bulk), theirs = reference(set, rate);
      const a = loop(ours, +rate), b = loop(theirs, +rate);
      const at = `${rate} V/s`;
      assert.ok(Math.abs(a.hi - b.hi) < 1.5e-3, `${at}: hysteresis index ${a.hi} vs ${b.hi}`);
      for (const s of ['rev', 'fwd']) {
        assert.ok(Math.abs(a[s].pmax - b[s].pmax) < 0.05, `${at} ${s}: P_max ${a[s].pmax} vs ${b[s].pmax} mW/cm²`);
        assert.ok(Math.abs(a[s].voc - b[s].voc) < 2e-3, `${at} ${s}: V_oc ${a[s].voc} vs ${b[s].voc} V`);
      }
      // The loop itself, point by point in time (up to 1.1 V, short of the steep rise past V_oc).
      let worst = 0;
      for (const [t, V, J] of theirs) {
        if (!Number.isFinite(J) || V > 1.1 || t <= 0) continue;
        const k = ours.findIndex(([tk]) => tk >= t);
        if (k <= 0) continue;
        const [[t0, , J0], [t1, , J1]] = [ours[k - 1], ours[k]];
        worst = Math.max(worst, Math.abs(J0 + ((J1 - J0) * (t - t0)) / (t1 - t0) - J));
      }
      assert.ok(worst < 0.3, `${at}: J differs by up to ${worst} mA/cm²`);
    }
  });
}

test('the same cell in steady state: every pair made is collected or recombines (J = F(G − R_bulk − R_faces)), and describe() reads it right', async () => {
  const { describe } = await import('../src/kit.js');
  for (const V of [0, 0.9, 1.1]) {
    const s = new Device(cell(V, true)).solve();
    assert.ok(s.converged);
    const [srh, light] = s.bulkReactions, faces = s.interfaces[0].rates[0] + s.interfaces[1].rates[0];
    assert.ok(Math.abs(light.total / (m3(P.Fph) * -Math.expm1(-P.alpha * P.b)) - 1) < 1e-9, 'Beer–Lambert, absorbed in full');
    const J = FARADAY * (light.total - srh.total - faces);
    assert.ok(Math.abs(s.current / J - 1) < 1e-6, `${V} V: J ${s.current} vs F(G − R) ${J}`);
  }
  // The perovskite's Debye length is set by its vacancies alone (the background doesn't move),
  // and the SRH reaction is listed where it runs.
  const text = describe(cell(1.2, true));
  const lambda = Math.sqrt((P.eps * 8.8541878128e-12 * GAS_CONSTANT * T) / (FARADAY * FARADAY * m3(P.N0)));
  assert.match(text, new RegExp(`MAPbI₃ .*Debye length ${(lambda * 1e9).toPrecision(3)} nm`));
  assert.match(text, /e- \+ h\+ = 0 \(SRH\) in MAPI/);
});

