// Shared by the demo pages: the page's own code shown as a worked example, and level diagrams.

import { bandDiagram, levelChart, THEME } from '../src/plot.js';
import { traces, speciesRole } from '../src/kit.js';

// driftlet/plot's theme colours as page-wide CSS variables (light and dark), so the canvas
// charts draw each species in the colour its level-diagram lines have.
{
  const vars = (mode) =>
    Object.entries(THEME[mode])
      .flatMap(([role, list]) => list.map((c, k) => `--driftlet-${role}${role === 'electron' ? '' : `-${k + 1}`}: ${c};`))
      .join(' ');
  const style = document.createElement('style');
  style.textContent =
    `:root { ${vars('light')} } @media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { ${vars('dark')} } } ` +
    `:root[data-theme="dark"] { ${vars('dark')} }`;
  document.head.append(style);
}

// A species' colour, as in the level diagrams (a CSS variable, which chart() resolves).
export function speciesColor(sol, name) {
  const { role, slot } = speciesRole(sol, name);
  return role === 'electron' ? 'var(--driftlet-electron)' : `var(--driftlet-${role}-${(slot % THEME.light[role].length) + 1})`;
}

// The page's module script, with its imports as they'd read from the package, in a panel.
export function showCode() {
  const script = document.querySelector('script[type="module"]');
  const code = script.textContent
    .replace(/^\n/, '')
    .replaceAll("'../src/index.js'", "'driftlet'")
    .replaceAll("'../src/kit.js'", "'driftlet/kit'")
    .replaceAll("'../src/plot.js'", "'driftlet/plot'");
  const details = document.createElement('details');
  details.className = 'code';
  details.innerHTML = '<summary>The code (the whole page script; <code>chart</code> and <code>common.js</code> are this site’s small helpers)</summary><pre><code></code></pre>';
  details.querySelector('pre code').textContent = code;
  if (!new URLSearchParams(location.search).has('shot')) document.querySelector('main').append(details);
}

// A level diagram from driftlet/plot into an element.
export function levels(el, sol, opts = {}) {
  el.innerHTML = bandDiagram(sol, { width: 720, height: 380, ...opts });
}

// Ions' species voltages sit volts apart (their μ° differ), and often only their changes matter:
// each drawn relative to its value at the left end, in mV. (traces() is plain data, edited here.)
export function relativeLevels(el, sol, { species, title, labels = {} }) {
  const tr = traces(sol, { species, standard: false, labels });
  for (const line of tr.series) {
    const ref = line.y.find(Number.isFinite);
    line.y = line.y.map((v) => 1000 * (v - ref));
  }
  const all = tr.series.flatMap((line) => [...line.y].filter(Number.isFinite));
  const lo = Math.min(...all), hi = Math.max(...all), pad = 0.08 * (hi - lo || 1);
  el.innerHTML = levelChart(tr, { width: 720, height: 360, range: [lo - pad, hi + pad], ylabel: 'V_i − V_i(left) (mV)', title });
}

// A slider's value from the URL (?V=0.4), for sharing a state or taking a screenshot.
export function fromQuery(input) {
  const v = new URLSearchParams(location.search).get(input.id);
  if (v !== null) input.value = v;
}

// ?shot: just the title and the figures (for the gallery's thumbnails and the README).
if (new URLSearchParams(location.search).has('shot')) {
  const style = document.createElement('style');
  style.textContent = 'main > p, .controls, details.code { display: none !important; } h1 { margin-bottom: 10px; }';
  document.head.append(style);
}
