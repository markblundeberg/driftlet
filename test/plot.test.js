import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { build, layer, ohmic, half, traces, typeset } from '../src/kit.js';
import { bandDiagram, levelChart } from '../src/plot.js';

// Plot-ready traces, and the SVG level diagram drawn from them.

const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 });
const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1 };
const cell = () =>
  new Device(
    build({
      species: [
        { name: 'Ag+', z: 1, cRef: 1000 },
        { name: 'NO3-', z: -1, cRef: 1000 },
        { name: 'e-', z: -1 },
      ],
      materials: {
        water: { epsr: 78.5, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } },
        Ag: { conductor: { species: 'e-', conductivity: 6.3e7 } },
      },
      stack: [
        ohmic(0, ['e-']),
        layer('Ag', 1e-6),
        { ...stern, reactions: [{ ...silver, k0: 1e-3 }] },
        layer('water', 10e-6, { name: 'electrolyte', c0: { 'Ag+': 10, 'NO3-': 10 } }),
        { ...stern, reactions: [{ ...silver, k0: 1e-3 }] },
        layer('Ag', 1e-6),
        ohmic(0.2, ['e-']),
      ],
      grid: { hmin: 0.1e-9, hmax: 100e-9, ratio: 1.15 },
    }),
  ).solve();

test('traces: a line per level, colour slots that follow the species, regions and a range', () => {
  const sol = cell();
  const tr = traces(sol, { phi: true, levels: [{ half: silver, label: 'Ag+/Ag' }] });
  assert.deepEqual(
    tr.series.map((s) => [s.id, s.kind, s.role, s.slot]),
    [
      ['V:Ag+', 'level', 'cation', 0],
      ['Vstd:Ag+', 'standard', 'cation', 0],
      ['V:NO3-', 'level', 'anion', 0],
      ['Vstd:NO3-', 'standard', 'anion', 0],
      ['V:e-', 'level', 'electron', 0], // a metal's Fermi level, with no standard level
      ['phi', 'phi', 'phi', 0],
      ['level:0', 'redox', 'redox', 0],
    ],
  );
  // Showing fewer species doesn't repaint the rest; a couple's two levels share a colour.
  assert.deepEqual(traces(sol, { species: ['NO3-'] }).series.map((s) => [s.role, s.slot]), [['anion', 0], ['anion', 0]]);
  const both = traces(sol, { species: [], levels: [{ half: silver }, { half: silver, standard: true }] }).series;
  assert.deepEqual(both.map((s) => [s.kind, s.slot]), [['redox', 0], ['redox-standard', 0]]);
  // A display offset moves all of a species' lines, and says so.
  const shifted = traces(sol, { species: ['NO3-'], shifts: { 'NO3-': -1.5 } }).series[0];
  const g = sol.x.length >> 1;
  assert.ok(Math.abs(shifted.y[g] - (sol.V['NO3-'][g] - 1.5)) < 1e-12 && shifted.shift === -1.5 && /⌇−1.5 V$/.test(shifted.label));
  assert.deepEqual(tr.regions.map((r) => [r.name, r.material]), [['region 0', 'Ag'], ['electrolyte', 'water'], ['region 2', 'Ag']]);
  assert.ok(Math.abs(tr.regions[1].x1 - 11e-6) < 1e-15 && tr.faces.length === 2 && tr.faces[0] === tr.regions[1].x0);
  for (const s of tr.series) for (const v of s.y) if (Number.isFinite(v)) assert.ok(v > tr.range[0] && v < tr.range[1]);
  assert.throws(() => traces(sol, { species: ['Cl-'] }), /no charged species 'Cl-'/);
});

test('the SVG level diagram: one path per line, broken where a level is undefined, a legend, a label', () => {
  const sol = cell();
  const svg = bandDiagram(sol, { title: 'silver cell' });
  assert.match(svg, /^<svg [^>]*role="img"[^>]*aria-label="silver cell: V Ag⁺, V° Ag⁺, V NO₃⁻, V° NO₃⁻, V e⁻ against x"/);
  assert.ok(svg.endsWith('</svg>'));
  // Tags balance (no parser in node, so count them).
  for (const tag of ['svg', 'g', 'path', 'text', 'style', 'clipPath', 'title']) {
    const open = svg.match(new RegExp(`<${tag}[ >]`, 'g'))?.length ?? 0, close = svg.match(new RegExp(`</${tag}>`, 'g'))?.length ?? 0;
    assert.equal(open, close, tag);
  }
  const paths = [...svg.matchAll(/<path d="([^"]*)"[^>]*><title>([^<]*)<\/title>/g)].map((m) => ({ d: m[1], label: m[2] }));
  // Thin standard levels first, so the thick species voltages sit on top.
  assert.deepEqual(paths.map((p) => p.label), ['V° Ag⁺', 'V° NO₃⁻', 'V Ag⁺', 'V NO₃⁻', 'V e⁻']);
  const byLabel = Object.fromEntries(paths.map((p) => [p.label, p]));
  // Ag⁺ lives in the solution only: one stretch. The Fermi level is in both metals: two.
  const moves = (p) => p.d.match(/M/g).length;
  assert.equal(moves(byLabel['V Ag⁺']), 1);
  assert.equal(moves(byLabel['V e⁻']), 2);
  assert.match(svg, />electrolyte \(water\)<\/text>/);
  // A single line needs no legend.
  const one = levelChart(traces(sol, { species: ['e-'] }));
  assert.equal((one.match(/<line [^>]*stroke-width=/g) ?? []).length, 0);
  // Species voltages thick, standard levels thin; a shifted species is marked and explained.
  assert.match(svg, /stroke-width="2.8"><title>V Ag⁺<\/title>/);
  assert.match(svg, /stroke-width="1.3"><title>V° Ag⁺<\/title>/);
  const marked = bandDiagram(sol, { shifts: { 'NO3-': -1 } });
  assert.match(marked, /⌇ = per-species offset/);
  assert.match(marked, /<path d="M[^"]*l6 4.5l-6 4.5[^"]*" stroke="var\(--driftlet-ink\)"[^>]*><title>V NO₃⁻ ⌇−1 V<\/title>/);
  // Two diagrams in one page don't share ids.
  assert.notEqual(svg.match(/id="([^"]+)"/)[1], one.match(/id="([^"]+)"/)[1]);
});

test('zooming with xlim: only the regions in view, ticks over the window, the range fitted to it', () => {
  const sol = cell();
  const full = bandDiagram(sol), zoom = bandDiagram(sol, { xlim: [0.9e-6, 1.1e-6] });
  assert.match(full, />electrolyte \(water\)<\/text>/);
  // The window covers the first face: metal and solution, and ticks in nm.
  assert.match(zoom, /x \(nm\)/);
  assert.equal((zoom.match(/stroke-dasharray="2 3"/g) ?? []).length, 1, 'one face in view');
  const visible = traces(sol).series.flatMap((s) => [...s.y].filter((v, g) => Number.isFinite(v) && sol.x[g] >= 0.9e-6 && sol.x[g] <= 1.1e-6));
  const ticks = [...zoom.matchAll(/text-anchor="end"[^>]*>([-\d.e]+)</g)].map((m) => +m[1]);
  assert.ok(Math.min(...ticks) >= Math.min(...visible) - 0.2 && Math.max(...ticks) <= Math.max(...visible) + 0.2);
});

test('ytick writes the y axis labels, escaped: say powers of ten for log profiles', () => {
  const sol = cell();
  const tr = { x: sol.x, regions: sol.regions, faces: [], series: [{ id: 'c', label: 'c', kind: 'level', role: 'cation', slot: 0, y: Array.from(sol.c['Ag+'], (c) => Math.log10(c)) }] };
  const svg = levelChart(tr, { range: [-3, 2], ytick: (v) => `<10^${v}>` });
  const labels = [...svg.matchAll(/text-anchor="end"[^>]*>([^<]*)</g)].map((m) => m[1]);
  assert.ok(labels.includes('&lt;10^-2&gt;') && labels.includes('&lt;10^1&gt;'), labels.join(' '));
});

test('typeset: a species name with its charge as a superscript and its counts as subscripts', () => {
  assert.equal(typeset('SO42-', -2), 'SO₄²⁻');
  assert.equal(typeset('Fe(CN)63-', -3), 'Fe(CN)₆³⁻');
  assert.equal(typeset('e-', -1), 'e⁻');
  assert.equal(typeset('H2PO4-', -1), 'H₂PO₄⁻');
  assert.equal(typeset('Li', 0), 'Li');
  assert.equal(typeset('odd', 1), 'odd', 'a name not ending in its charge is left alone');
  // Redox levels' default labels typeset the species they know.
  const tr = traces(cell(), { species: [], levels: [{ half: silver }] });
  assert.equal(tr.series[0].label, 'Ag⁺ + e⁻ = Ag(s)');
});
