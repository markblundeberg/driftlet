import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units } from '../src/index.js';
import { build, layer, ohmic, bath, half, traces, speciesRole, semiconductor, metal, aqueous } from '../src/kit.js';
import { bandDiagram, levelChart } from '../src/plot.js';

// Round 6's plot chase: pictures that stated something false (an axis titled in volts over
// energies, tick labels away from their ticks, a line leaving the frame), and ranges whose ticks
// never ended.

const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 });
const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 };
const cell = () =>
  new Device(
    build({
      species: [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }],
      materials: {
        water: { epsr: 78.5, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } },
        Ag: { conductor: { species: 'e-', conductivity: 6.3e7 } },
      },
      stack: [ohmic(0, ['e-']), layer('Ag', 1e-6), { ...stern, reactions: [{ ...silver, k0: 1e-3 }] },
        layer('water', 10e-6, { name: 'electrolyte', c0: { 'Ag+': 10, 'NO3-': 10 } }),
        { ...stern, reactions: [{ ...silver, k0: 1e-3 }] }, layer('Ag', 1e-6), ohmic(0.2, ['e-'])],
      grid: { hmin: 0.1e-9, hmax: 100e-9, ratio: 1.15 },
    }),
  ).solve();
const pn = (V) =>
  new Device(build({
    T: 300, library: [semiconductor('Si')],
    stack: [ohmic(0), layer('Si', 1e-6, { name: 'n', donors: units.perCm3(1e17) }), layer('Si', 1e-6, { name: 'p', acceptors: units.perCm3(1e16) }), ohmic(V)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
  })).solve();
const one = (x, y, extra = {}) => ({
  x: Float64Array.from(x), series: [{ id: 'a', label: 'a', kind: 'level', role: 'electron', slot: 0, y: Float64Array.from(y), ...extra }],
  regions: [{ name: 'r', material: 'm', x0: x[0], x1: x[x.length - 1] }], faces: [],
});

// What an SVG says: its y and x tick labels with their pixel positions, the frame, the paths.
const yTicks = (svg) => [...svg.matchAll(/<text x="50" y="([^"]*)" text-anchor="end"[^>]*>([^<]*)<\/text>/g)].map((m) => ({ px: +m[1] - 4, label: m[2] }));
const xTicks = (svg) => [...svg.matchAll(/<text x="([^"]*)" y="[^"]*" text-anchor="middle" fill="var\(--driftlet-muted\)">([-\d.e]+)<\/text>/g)].map((m) => ({ px: +m[1], label: m[2] }));
const frame = (svg) => {
  const m = svg.match(/<rect x="([^"]*)" y="([^"]*)" width="([^"]*)" height="([^"]*)" fill="none"/);
  return { x0: +m[1], y0: +m[2], x1: +m[1] + +m[3], y1: +m[2] + +m[4] };
};
const ylabelOf = (svg) => svg.match(/rotate\(-90\)[^>]*>(.*?)<\/text>/)[1].replace(/<[^>]+>/g, '');
const paths = (svg) => [...svg.matchAll(/<path d="([^"]*)" stroke="([^"]*)"([^>]*)>(?:<title>([^<]*)<\/title>)?/g)].map((m) => ({ d: m[1], stroke: m[2], rest: m[3], title: m[4] }));
// value → pixel on the y axis, from two of its tick labels
const yMap = (svg) => {
  const t = yTicks(svg);
  const [a, b] = [t[0], t[t.length - 1]];
  return (v) => a.px + ((v - +a.label) * (b.px - a.px)) / (+b.label - +a.label);
};

test('#2: traces say whether they are energies, and levelChart titles its axis from that', () => {
  const sol = pn(0.3);
  const tr = traces(sol, { energy: true });
  assert.equal(tr.energy, true);
  assert.equal(traces(sol).energy, false);
  assert.equal(ylabelOf(levelChart(tr)), 'electron energy (eV)');
  assert.equal(ylabelOf(levelChart(tr)), ylabelOf(bandDiagram(sol, { energy: true })));
  assert.equal(ylabelOf(levelChart(traces(sol))), 'voltage (V)');
  assert.equal(ylabelOf(levelChart(one([0, 1e-6], [0, 1]))), 'voltage (V)'); // hand-made traces
  assert.equal(ylabelOf(levelChart(tr, { ylabel: 'E (eV)' })), 'E (eV)');
});

test('#6: tick labels carry the decimals their step needs, so each names the position it marks', () => {
  // [1, 1.00004]: four significant digits wrote 1, 1, 1, 1
  const labels = yTicks(levelChart(one([0, 1e-6], [1, 1.00004]), { range: [1, 1.00004] })).map((t) => t.label);
  assert.equal(new Set(labels).size, labels.length, `distinct labels: ${labels}`);
  // An axis over [10.98, 11.00] µm, in nm: "10990" was drawn at 10985 nm.
  const xs = 10.98e-6, xe = 11e-6;
  const svg = levelChart(one([0, 11e-6], [0, 1]), { xlim: [xs, xe] });
  const f = frame(svg), ticks = xTicks(svg);
  assert.ok(ticks.length >= 4);
  for (const t of ticks) {
    const at = (xs + ((t.px - f.x0) / (f.x1 - f.x0)) * (xe - xs)) * 1e9;
    assert.ok(Math.abs(+t.label - at) < 0.01, `label ${t.label} drawn at ${at} nm`);
  }
  // and the y axis: each label where its value is, against the other labels
  const y = yTicks(levelChart(one([0, 1e-6], [0.3, 0.3000001]), { range: [0.3, 0.3000001] }));
  assert.equal(new Set(y.map((t) => t.label)).size, y.length);
  const step = (+y[1].label - +y[0].label) / (y[1].px - y[0].px);
  for (const t of y) assert.ok(Math.abs(+t.label - (+y[0].label + (t.px - y[0].px) * step)) < 1e-12, `${t.label}`);
});

test('#11: an energy-up diagram names a standard redox level as an energy, E°, not V°', () => {
  const sol = cell();
  const levels = [{ half: silver }, { half: silver, standard: true }];
  assert.deepEqual(traces(sol, { energy: true, levels }).series.slice(1).map((s) => s.label), ['Ag⁺ + e⁻ = Ag(s)', 'E° Ag⁺ + e⁻ = Ag(s)']);
  assert.deepEqual(traces(sol, { species: [], levels }).series.map((s) => s.label), ['Ag⁺ + e⁻ = Ag(s)', 'V° Ag⁺ + e⁻ = Ag(s)']);
});

test('#4: traces are copies, so editing one in place leaves the solution alone', () => {
  const sol = cell();
  const keep = { V: Float64Array.from(sol.V['Ag+']), Vstd: Float64Array.from(sol.Vstd['Ag+']), phi: Float64Array.from(sol.phi), x: Float64Array.from(sol.x) };
  const tr = traces(sol, { phi: true });
  for (const s of tr.series) for (let k = 0; k < s.y.length; k++) s.y[k] *= 1000; // to mV, in place
  for (let k = 0; k < tr.x.length; k++) tr.x[k] *= 1e6; // to µm, in place
  assert.deepEqual(sol.V['Ag+'], keep.V);
  assert.deepEqual(sol.Vstd['Ag+'], keep.Vstd);
  assert.deepEqual(sol.phi, keep.phi);
  assert.deepEqual(sol.x, keep.x);
  assert.notEqual(traces(sol).series[0].y, sol.V['Ag+']);
});

test('#1, #7: the range fits what is drawn: the lines at the window edges, and edited lines', () => {
  // A line rising steeply past the window's right edge: the fit to the nodes inside let it out.
  const x = [0, 1e-9, 2e-9, 3e-9], y = [0, 1, 2, 10];
  const tr = { ...one(x, y), range: [-0.6, 10.6] };
  for (const [xlim, inView] of [[[0.5e-9, 2.5e-9], [0.5, 6]], [[2.2e-9, 2.8e-9], [3.6, 8.4]]]) {
    const svg = levelChart(tr, { xlim });
    const py = yMap(svg), f = frame(svg);
    for (const v of inView) assert.ok(py(v) >= f.y0 && py(v) <= f.y1, `xlim ${xlim}: ${v} at ${py(v)} px, frame ${f.y0}–${f.y1}`);
    // and fitted: the frame not much taller than what's in view
    const span = Math.abs(py(inView[1]) - py(inView[0]));
    assert.ok(span > 0.8 * (f.y1 - f.y0), `xlim ${xlim}: the lines span ${span} px of ${f.y1 - f.y0}`);
  }
  // Edited traces with no range: their cached range was used, and the plot came out empty.
  const sol = cell();
  const ag = traces(sol, { species: ['Ag+'], standard: false });
  for (const s of ag.series) s.y = s.y.map((v) => 1000 * v);
  const svg = levelChart(ag, { ylabel: 'mV' }), f = frame(svg);
  const ys = [...paths(svg)[0].d.matchAll(/[ML][-\d.]+ ([-\d.]+)/g)].map((m) => +m[1]);
  assert.ok(Math.min(...ys) >= f.y0 && Math.max(...ys) <= f.y1, `line from ${Math.min(...ys)} to ${Math.max(...ys)} px, frame ${f.y0}–${f.y1}`);
  // Unedited, the same range as before: traces' suggestion, padded 6%.
  const fresh = traces(sol, { species: ['Ag+'], standard: false });
  assert.deepEqual(yTicks(levelChart(fresh)), yTicks(levelChart(fresh, { range: fresh.range })));
});

test('#12: a metal-only device has one Fermi level, E_F, whatever semiconductor its library carries', () => {
  for (const library of [[metal('Pt')], [metal('Pt'), semiconductor('Si')]]) {
    const sol = new Device(build({ library, stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), ohmic(0.1, ['e-'])] })).solve();
    assert.deepEqual(traces(sol, { energy: true }).series.map((s) => s.label), ['E_F']);
  }
  assert.ok(traces(pn(0), { energy: true }).series.some((s) => s.label === 'E_{Fn}'));
});

test('#5, #10: ranges an ulp wide or over ±∞ draw (or refuse) at once, never hang', () => {
  for (const lo of [1, 0.2, -2.5, 1e3]) {
    for (const ulps of [1, 2, 3]) {
      const hi = lo + ulps * Number.EPSILON * Math.abs(lo);
      const t0 = Date.now();
      const tr = one([0, 1e-6], [lo, hi]);
      const svg = levelChart(tr);
      levelChart(tr, { range: [lo, hi] });
      assert.ok(Date.now() - t0 < 1000);
      const labels = yTicks(svg).map((t) => t.label);
      assert.ok(labels.length >= 2 && new Set(labels).size === labels.length, `${lo} + ${ulps} ulp: ${labels}`);
      // traces' range reads it as flat too
      assert.ok(levelChart({ ...tr, range: [lo - 0.5, hi + 0.5] }).length > 0);
    }
  }
  const sol = cell();
  const tr = traces(sol, { species: ['Ag+'], standard: false });
  for (const range of [[-Infinity, Infinity], [0, Infinity], [NaN, 1], [1, 1], [1, 0]]) {
    assert.throws(() => levelChart(tr, { range }), /levelChart: range must be \[lo, hi\], finite, with lo < hi/);
  }
  assert.throws(() => levelChart(tr, { xlim: [0, Infinity] }), /levelChart: xlim must be/);
  // An infinite value in the lines: skipped by traces' range and by the fit.
  const g = sol.V['Ag+'].findIndex(Number.isFinite), keep = sol.V['Ag+'][g];
  sol.V['Ag+'][g] = Infinity;
  try {
    const inf = traces(sol, { species: ['Ag+'], standard: false });
    assert.ok(inf.range.every(Number.isFinite), `${inf.range}`);
    assert.ok(yTicks(levelChart(inf)).length >= 2);
  } finally {
    sol.V['Ag+'][g] = keep;
  }
});

test('#3: a level lying on another of its weight is drawn dashed over it, so both colours show', () => {
  const dashed = (svg, title) => / stroke-dasharray=/.test(paths(svg).find((p) => p.title === title).rest);
  const eq = bandDiagram(pn(0), { energy: true });
  assert.equal(dashed(eq, 'E_Fn'), false);
  assert.equal(dashed(eq, 'E_Fp'), true); // drawn last, over E_Fn, at equilibrium
  // Under bias they part, and both are solid; the band edges never coincide.
  const bias = bandDiagram(pn(0.3), { energy: true });
  for (const t of ['E_Fn', 'E_Fp', 'E_c', 'E_v']) assert.equal(dashed(bias, t), false, t);
  // The legend's swatches stay solid: the line is the same kind of line.
  assert.equal((eq.match(/<line [^>]*stroke-width="2.8" stroke-dasharray/g) ?? []).length, 0);
});

test('#8: a couple at two fixed μ is two levels, with two colours and labels that say which', () => {
  const sol = new Device(build({ library: [aqueous(['H+', 'Cl-'], { epsr: 0 })], stack: [bath({ 'H+': 100, 'Cl-': 100 }, 'Cl-', 0), layer('water', 1e-6), bath({ 'H+': 100, 'Cl-': 100 }, 'Cl-', 0)] })).solve();
  const a = half('2 H+ + 2 e- = H2', { H2: 0 }), b = half('2 H+ + 2 e- = H2', { H2: 1e4 });
  const tr = traces(sol, { species: [], levels: [{ half: a }, { half: b }, { half: b, standard: true }] });
  assert.deepEqual(tr.series.map((s) => s.slot), [0, 1, 1]); // a couple's two levels still share one
  assert.equal(new Set(tr.series.map((s) => s.label)).size, 3);
  assert.match(tr.series[1].label, /10 kJ\/mol/);
  // One equation alone keeps its plain label.
  assert.equal(traces(sol, { species: [], levels: [{ half: a }] }).series[0].label, '2 H⁺ + 2 e- = H2');
});

test('#9: a lone shifted line gives its offset in the note, with no legend to carry it', () => {
  const sol = cell();
  const lone = bandDiagram(sol, { species: ['Ag+'], standard: false, shifts: { 'Ag+': 1 } });
  assert.match(lone, />⌇ = per-species offset \(\+1 V\)</);
  const two = bandDiagram(sol, { species: ['Ag+'], shifts: { 'Ag+': -0.5 } });
  assert.match(two, />⌇ = per-species offset</);
  assert.doesNotMatch(two, /offset \(/);
});

test('#15: a species named constructor reads its own shift and role, not Object’s', () => {
  const sol = new Device(build({
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'constructor', z: -1, cRef: 1000 }],
    materials: { water: { epsr: 78.5, species: { 'K+': { D: 1.96e-9, mu0: -282.5e3 }, constructor: { D: 2e-9, mu0: -131.2e3 } } } },
    stack: [bath({ 'K+': 100, constructor: 100 }, 'constructor', 0), layer('water', 1e-6), bath({ 'K+': 10, constructor: 10 }, 'constructor', { I: 0 })],
  })).solve();
  const s = traces(sol, { standard: false }).series.find((t) => t.species === 'constructor');
  assert.equal(s.label, 'V_{constructor}');
  assert.equal(s.shift, undefined);
  assert.ok(s.y.every(Number.isFinite));
  assert.deepEqual(speciesRole(sol, 'constructor'), { role: 'anion', slot: 0 });
  assert.throws(() => speciesRole(cell(), 'constructor'), /no species 'constructor'/);
  assert.throws(() => speciesRole(cell(), 'toString'), /no species 'toString'/);
});
