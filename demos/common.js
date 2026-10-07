// Shared by the demo pages: the page's own code shown as a worked example, level diagrams, and
// profiles drawn to line up beneath them.

import { bandDiagram, levelChart } from '../src/plot.js';
import { speciesRole, typeset } from '../src/kit.js';

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

// Concentration profiles drawn like the level diagrams (same layout, x axis and species colours),
// so they line up beneath them: levelChart() takes any lines as traces. Each entry of `species` is
// a name, or { name, c, label, kind } to draw some other profile in that species' colour (kind
// 'standard' draws it thin), or { c, label, kind, role } for a profile of no species (role
// 'redox', or none for muted). `scale` converts units; `log` draws log₁₀, with powers of ten.
export function profiles(el, sol, { species, title, ylabel = 'concentration (mol/m³)', range, xlim, scale = 1, log = false, labels = {}, width = 720, height = 300 }) {
  const f = log ? (c) => (c > 0 ? Math.log10(c * scale) : NaN) : (c) => c * scale;
  const series = species.map((entry) => {
    const { name, c = sol.c[name], label = labels[name] ?? `c_{${typeset(name, sol.species.find((sp) => sp.name === name).z)}}`, kind = 'level', role } = typeof entry === 'string' ? { name: entry } : entry;
    return { id: `c:${name}:${label}`, label, kind, ...(name === undefined ? { role, slot: 0 } : speciesRole(sol, name)), y: Array.from(c, f) };
  });
  if (!range) {
    const [xs, xe] = xlim ?? [sol.x[0], sol.x.at(-1)];
    const shown = series.flatMap((s) => s.y.filter((v, g) => Number.isFinite(v) && sol.x[g] >= xs && sol.x[g] <= xe));
    const lo = Math.min(...shown), hi = Math.max(...shown);
    range = log ? [Math.floor(lo - 0.2), Math.ceil(hi + 0.2)] : [Math.min(0, lo), 1.08 * hi || 1];
  }
  const tr = { x: sol.x, series, regions: sol.regions, faces: sol.regions.slice(1).map((r) => r.x0), range };
  el.innerHTML = levelChart(tr, { width, height, title, ylabel, range, xlim, ytick: log ? powerOfTen : undefined });
}
const superscript = (n) => String(n).replace(/[-0-9]/g, (d) => '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'['-0123456789'.indexOf(d)]);
const powerOfTen = (v) => (Number.isInteger(v) ? `10${superscript(v)}` : '');

// A slider's value from the URL (?V=0.4), for sharing a state or taking a screenshot.
export function fromQuery(input) {
  const v = new URLSearchParams(location.search).get(input.id);
  if (v !== null) input.value = v;
}

// ?shot: just the title and the figures (for the gallery's thumbnails and the README).
if (new URLSearchParams(location.search).has('shot')) {
  const style = document.createElement('style');
  style.textContent = 'nav.site, main > p, main > ul, main > h2, .controls, details.code, .illustration { display: none !important; } h1 { margin-bottom: 10px; }';
  document.head.append(style);
}

// A live demo's pause: a button in its controls, and a pause of its own while the demo is
// scrolled out of view or the tab is hidden (back on when it returns, unless paused by hand), so
// a page doesn't keep a CPU busy unseen. The demo's frame loop asks `paused` before requesting its
// next frame; `resume` (called when it's unpaused) restarts the loop where it stopped. `watch` is
// the element whose visibility counts (default: the controls' section, or the page).
export function pauser(controls, resume, { watch } = {}) {
  let byHand = false, hidden = false, offscreen = false;
  const button = document.createElement('button');
  button.id = `${controls.id || 'controls'}-pause`;
  button.type = 'button';
  const paused = () => byHand || hidden || offscreen;
  const update = (was) => {
    button.textContent = byHand ? 'resume' : 'pause';
    if (was && !paused()) resume();
  };
  button.addEventListener('click', () => {
    const was = paused();
    byHand = !byHand;
    update(was);
  });
  // (Screenshots and ?t= states run headless on virtual time: never pause those.)
  if (!new URLSearchParams(location.search).has('shot')) {
    document.addEventListener('visibilitychange', () => {
      const was = paused();
      hidden = document.hidden;
      update(was);
    });
    const target = watch ?? controls.closest?.('section') ?? document.querySelector('main');
    // (Where there's no IntersectionObserver, as in a test's stub document, it just keeps running.)
    if (typeof IntersectionObserver !== 'undefined' && target) {
      new IntersectionObserver(([entry]) => {
        const was = paused();
        offscreen = !entry.isIntersecting;
        update(was);
      }).observe(target);
    }
  }
  update(false);
  controls.append(button);
  return { get paused() { return paused(); } };
}
