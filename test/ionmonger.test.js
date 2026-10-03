import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Device, AVOGADRO, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { perovskiteCell, PEROVSKITE_SCANS, hysteresis, recorder } from '../src/kit.js';

// A benchmark against an independent code: J–V hysteresis in a planar perovskite solar cell,
// computed by IonMonger (finite elements in MATLAB/Octave, Courtier et al. 2019) and by
// driftlet, scan by scan from 1 mV/s to 1 kV/s. The reference scans and how they were made are
// in test/fixtures/ionmonger; the cell is the kit's perovskiteCell().
//
// TiO₂ | MAPbI₃ with mobile iodide vacancies | spiro-OMeTAD, IonMonger's default parameters,
// with bulk SRH on and off. Three species: electrons, holes, vacancies (over an equal immobile
// background). IonMonger's interface recombination is SRH at each face, saturating at one
// carrier's capture; written with the ETL's own electron density it's exactly driftlet's face SRH
// law, with n₁ = (d_E/n₀)·n_i at the TiO₂ face (n₀ the perovskite's electron density at the
// ETL's Fermi level) and n₁ = n_i at the spiro face. Levels are on the vacuum scale, as
// IonMonger's are.

const T = 298, m3 = (x) => x / AVOGADRO;
const P = { b: 400e-9, eps: 24.1, alpha: 1.3e7, N0: 1.6e25, Fph: 1.4e21 };
const cell = (V, bulkSRH) => perovskiteCell({ V, bulkSRH });

// Preconditioned at 1.2 V (steady state, light on), then 1.2 → 0 → 1.2 V, recorded: the trace
// in SI, and rows [t, V, J (mA/cm²)].
function scan(rate, bulkSRH) {
  const dev = new Device(cell(1.2, bulkSRH));
  assert.ok(dev.solve().converged);
  const half = 1.2 / rate;
  dev.set({ contacts: { right: { V: { t: [0, half, 2 * half], values: [1.2, 0, 1.2] } } } });
  const rec = recorder(dev, { times: [half] });
  const r = rec.advance(2 * half, { tol: 1e-4, dtMax: (2 * half) / 480 });
  assert.ok(r.converged && r.done);
  assert.deepEqual(rec.frames.map((f) => f.time), [0, half]);
  const { trace } = rec;
  return { trace, rows: trace.t.map((t, k) => [t, trace.voltage[k], trace.current[k] / 10]) };
}

const reference = (set, rate) =>
  readFileSync(new URL(`fixtures/ionmonger/${set}_${rate}.csv`, import.meta.url), 'utf8').trim().split('\n').map((line) => line.split(',').map(Number));
const asTrace = (rows) => ({ t: rows.map((r) => r[0]), voltage: rows.map((r) => r[1]), current: rows.map((r) => r[2] * 10) });
const RATES = ['0.001', '0.01', '0.1', '1', '10', '100', '1000'];

test("PEROVSKITE_SCANS is IonMonger's runs, read by hysteresis()", () => {
  for (const [set, key] of [['full', 'full'], ['scan', 'noBulkSRH']]) {
    assert.deepEqual(PEROVSKITE_SCANS[key].map((r) => r.rate), RATES.map(Number));
    for (const ref of PEROVSKITE_SCANS[key]) {
      const rows = reference(set, ref.rate), half = 1.2 / ref.rate, m = hysteresis(asTrace(rows), half);
      const near = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${set} ${ref.rate} V/s ${what}: ${a} vs ${b}`);
      near(m.hi, ref.hi, 5e-6, 'hi');
      near(rows.find(([t]) => Math.abs(t - half) < 1e-9 * half)[2] * 10, ref.J0, 0.006, 'J0');
      for (const s of ['rev', 'fwd']) {
        near(m[s].Pmax, ref[s].Pmax, 0.006, `${s} Pmax`);
        near(m[s].Voc, ref[s].Voc, 5e-6, `${s} Voc`);
      }
    }
  }
});

for (const [set, key, bulkSRH, what] of [['scan', 'noBulkSRH', false, 'bulk SRH off'], ['full', 'full', true, "IonMonger's full defaults"]]) {
  test(`perovskite hysteresis against IonMonger (${what}), 1 mV/s to 1 kV/s: hysteresis index, maximum power, V_oc and the whole J–V loop (to 0.03 mA/cm² where gentle, 1 mV where steep)`, () => {
    for (const b of PEROVSKITE_SCANS[key]) {
      const rate = String(b.rate), { trace, rows: ours } = scan(b.rate, bulkSRH), theirs = reference(set, rate);
      const a = hysteresis(trace, 1.2 / b.rate);
      const at = `${rate} V/s`;
      assert.ok(Math.abs(a.hi - b.hi) < 1.5e-3, `${at}: hysteresis index ${a.hi} vs ${b.hi}`);
      for (const s of ['rev', 'fwd']) {
        assert.ok(Math.abs(a[s].Pmax - b[s].Pmax) < 0.5, `${at} ${s}: P_max ${a[s].Pmax} vs ${b[s].Pmax} W/m²`);
        assert.ok(Math.abs(a[s].Voc - b[s].Voc) < 2e-3, `${at} ${s}: V_oc ${a[s].Voc} vs ${b[s].Voc} V`);
      }
      // The loop itself, at IonMonger's own points (matched in time): as a current difference
      // where the curve is gentle, and as a voltage offset where it's steep (near V_oc it falls at
      // up to 1100 mA/cm² per volt, so a sub-millivolt shift would read as a large current gap).
      let gentle = 0, steep = 0;
      for (let q = 1; q < theirs.length - 1; q++) {
        const [t, V, J] = theirs[q];
        if (!Number.isFinite(J) || V > 1.15 || t <= 0) continue;
        const k = ours.findIndex(([tk]) => tk >= t);
        if (k <= 0) continue;
        const [[t0, , J0], [t1, , J1]] = [ours[k - 1], ours[k]];
        const gap = Math.abs(J0 + ((J1 - J0) * (t - t0)) / (t1 - t0) - J);
        const [, Va, Ja] = theirs[q - 1], [, Vb, Jb] = theirs[q + 1];
        const slope = Math.abs((Jb - Ja) / (Vb - Va)); // mA/cm² per V
        if (slope < 50) gentle = Math.max(gentle, gap);
        else steep = Math.max(steep, gap / slope);
      }
      assert.ok(gentle < 0.03, `${at}: J differs by up to ${gentle} mA/cm² where the curve is gentle`);
      assert.ok(steep < 1e-3, `${at}: the steep part is offset by up to ${steep * 1000} mV`);
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
  assert.match(text, new RegExp(`\\[MAPbI₃\\] .*Debye length ${(lambda * 1e9).toPrecision(3)} nm`));
  assert.match(text, /e- \+ h\+ = 0 \(SRH\) in MAPbI₃/);
});

