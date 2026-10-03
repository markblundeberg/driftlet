import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError } from '../src/index.js';
import { build, layer, bath, recorder } from '../src/kit.js';

// The recorder on a transient with a closed form: KCl in a strictly neutral electrolyte between
// a bath (c̄) and a blocking wall, laid down as c̄ + a·sin(πx/2L), its one diffusion mode, which
// decays as e^{−t/τ}, 1/τ = π²D/4L² with D = 2D₊D₋/(D₊ + D₋).
const L = 100e-6, cbar = 10, amp = 5, Dp = 1.96e-9, Dm = 2.03e-9;
const tau = (4 * L * L) / (Math.PI * Math.PI * ((2 * Dp * Dm) / (Dp + Dm)));
const mode = (x) => Math.sin((Math.PI * x) / (2 * L));
const salt = () => {
  const x = Array.from({ length: 401 }, (_, k) => (k * L) / 400), profile = { x, values: x.map((v) => cbar + amp * mode(v)) };
  return new Device(
    build({
      species: [
        { name: 'K+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'K+': { D: Dp, mu0: 0 }, 'Cl-': { D: Dm, mu0: 0 } } } },
      stack: [bath({ 'K+': cbar, 'Cl-': cbar }, 'Cl-'), layer('water', L, { c0: { 'K+': profile, 'Cl-': profile } }), { species: { 'K+': 'blocked', 'Cl-': 'blocked' }, phi: 'neutral' }],
      grid: { hmin: 2e-6, hmax: 2e-6 },
    }),
  );
};
// The mode's amplitude in a solution, (2/L)∫(c − c̄) sin(πx/2L) dx.
const amplitude = (s) => {
  let m = 0;
  for (let g = 1; g < s.x.length; g++) for (const k of [g - 1, g]) m += (s.c['K+'][k] - cbar) * mode(s.x[k]) * ((s.x[g] - s.x[g - 1]) / 2);
  return (2 / L) * m;
};
const probes = [{ x: L, species: 'K+' }];

test('recorder: a frame every τ/4, each the decaying mode at its own time; the trace and probes run unbroken from the start', () => {
  const dev = salt(), rec = recorder(dev, { every: tau / 4, probes });
  const s = rec.advance(2 * tau, { tol: 1e-7 });
  assert.ok(s.converged && s.done);
  assert.equal(rec.frames.length, 9);
  const a0 = amplitude(rec.frames[0]);
  rec.frames.forEach((f, k) => {
    assert.ok(Math.abs(f.time / ((k * tau) / 4) - 1) < 1e-12 || (k === 0 && f.time === 0), `frame ${k} at ${f.time}`);
    const ratio = amplitude(f) / (a0 * Math.exp(-f.time / tau));
    assert.ok(Math.abs(ratio - 1) < 3e-5, `frame ${k}: ${ratio} of theory`);
  });
  // Frames are separate snapshots, not views of one state.
  assert.notEqual(rec.frames[1].c['K+'][5], rec.frames[2].c['K+'][5]);
  // The trace starts at the start, never steps back, and lands on every frame; the probe at the
  // wall reads c̄ + a·e^{−t/τ} all along (the start included).
  const { t, probes: [wall] } = rec.trace;
  assert.equal(t[0], 0);
  for (let k = 1; k < t.length; k++) assert.ok(t[k] > t[k - 1]);
  for (const f of rec.frames) assert.ok(t.includes(f.time));
  assert.equal(wall.length, t.length);
  wall.forEach((c, k) => assert.ok(Math.abs((c - cbar) / (a0 * Math.exp(-t[k] / tau)) - 1) < 1e-3, `t = ${t[k]}: ${c}`));
  // Scrubbing: the frame at or before a time.
  assert.equal(rec.frame(0.3 * tau), rec.frames[1]);
  assert.equal(rec.frame(tau / 2), rec.frames[2]);
  assert.equal(rec.frame(-1), rec.frames[0]);
  assert.equal(rec.frame(10 * tau), rec.frames[8]);
});

test('recorder: a run in small time budgets (an animation) is the straight run exactly, and without a schedule each call ends in a frame', () => {
  const straight = recorder(salt(), { times: [tau / 3, tau, 0.5 * tau] });
  straight.advance(1.5 * tau, { tol: 1e-6 });
  const dev = salt(), rec = recorder(dev, { times: [tau / 3, tau, 0.5 * tau] });
  let calls = 0, s;
  do {
    s = rec.advance(1.5 * tau, { tol: 1e-6, budgetMs: 0.5 });
    calls++;
  } while (s.converged && !s.done && calls < 10000);
  assert.ok(s.done && calls > 1, `${calls} calls`);
  assert.deepEqual(rec.frames.map((f) => f.time), [0, tau / 3, 0.5 * tau, tau]);
  assert.deepEqual(rec.trace, straight.trace);
  assert.deepEqual(rec.frames.at(-1).c['K+'], straight.frames.at(-1).c['K+']);

  const free = recorder(salt());
  for (const t of [0.1, 0.2, 0.4]) free.advance(t * tau);
  free.advance(0.4 * tau); // already there: no new frame
  assert.deepEqual(free.frames.map((f) => f.time / tau), [0, 0.1, 0.2, 0.4]);
});

test('recorder options are checked', () => {
  const dev = salt();
  assert.throws(() => recorder(dev, { every: 1, times: [1] }), (e) => e instanceof DeviceError && /every or times/.test(e.message));
  assert.throws(() => recorder(dev, { every: 0 }), /every must be a time > 0/);
  assert.throws(() => recorder(dev, { times: [1, NaN] }), /times must be an array/);
  assert.throws(() => recorder(dev, { probes: [{ x: L, species: 'Na+' }] }), /no species "Na\+"/);
});
