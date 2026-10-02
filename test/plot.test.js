import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { build, layer, ohmic, half, traces } from '../src/kit.js';
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
    tr.series.map((s) => [s.id, s.kind, s.slot]),
    [
      ['V:Ag+', 'level', 0],
      ['Vstd:Ag+', 'standard', 0],
      ['V:NO3-', 'level', 1],
      ['Vstd:NO3-', 'standard', 1],
      ['V:e-', 'level', 2], // a metal's Fermi level, with no standard level
      ['phi', 'phi', 3],
      ['level:0', 'redox', 4],
    ],
  );
  // Showing fewer species doesn't repaint the rest.
  assert.deepEqual(traces(sol, { species: ['e-'] }).series.map((s) => s.slot), [2]);
  assert.deepEqual(tr.regions.map((r) => [r.name, r.material]), [['region 0', 'Ag'], ['electrolyte', 'water'], ['region 2', 'Ag']]);
  assert.ok(Math.abs(tr.regions[1].x1 - 11e-6) < 1e-15 && tr.faces.length === 2 && tr.faces[0] === tr.regions[1].x0);
  for (const s of tr.series) for (const v of s.y) if (Number.isFinite(v)) assert.ok(v > tr.range[0] && v < tr.range[1]);
  assert.throws(() => traces(sol, { species: ['Cl-'] }), /no charged species 'Cl-'/);
});

test('the SVG level diagram: one path per line, broken where a level is undefined, a legend, a label', () => {
  const sol = cell();
  const svg = bandDiagram(sol, { title: 'silver cell' });
  assert.match(svg, /^<svg [^>]*role="img"[^>]*aria-label="silver cell: V Ag\+, V° Ag\+, V NO3-, V° NO3-, V e- against x"/);
  assert.ok(svg.endsWith('</svg>'));
  // Tags balance (no parser in node, so count them).
  for (const tag of ['svg', 'g', 'path', 'text', 'style', 'clipPath', 'title']) {
    const open = svg.match(new RegExp(`<${tag}[ >]`, 'g'))?.length ?? 0, close = svg.match(new RegExp(`</${tag}>`, 'g'))?.length ?? 0;
    assert.equal(open, close, tag);
  }
  const paths = [...svg.matchAll(/<path d="([^"]*)"[^>]*><title>([^<]*)<\/title>/g)].map((m) => ({ d: m[1], label: m[2] }));
  assert.deepEqual(paths.map((p) => p.label), ['V Ag+', 'V° Ag+', 'V NO3-', 'V° NO3-', 'V e-']);
  // Ag⁺ lives in the solution only: one stretch. The Fermi level is in both metals: two.
  const moves = (p) => p.d.match(/M/g).length;
  assert.equal(moves(paths[0]), 1);
  assert.equal(moves(paths[4]), 2);
  assert.match(svg, />electrolyte \(water\)<\/text>/);
  // A single line needs no legend.
  const one = levelChart(traces(sol, { species: ['e-'] }));
  assert.equal((one.match(/<line [^>]*stroke-width="2"/g) ?? []).length, 0);
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
