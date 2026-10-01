// driftlet/plot: a level diagram as an SVG string, with no dependencies and no DOM, so it
// works the same in a page, a Web Worker or Node. Put it in a page with
// `element.innerHTML = bandDiagram(sol)`. For your own plotting, `traces()` in driftlet/kit gives
// the same lines as data.
//
// Colour follows the species (its slot), solid lines are species voltages V_i, dashed ones
// standard levels V°_i, and regions are shaded bands labelled with their material. The colours
// are CSS custom properties (--driftlet-1 … --driftlet-8, --driftlet-ink, …) with light and dark
// defaults, so a page can restyle them.

import { traces } from './traces.js';

// Categorical hues in a fixed order (validated for colour-vision deficiency, light and dark).
const LIGHT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948'];
const DARK = ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767'];

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

function ticks(a, b, n) {
  const step0 = (b - a) / n;
  const mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= step0) ?? 10 * mag;
  const out = [];
  for (let t = Math.ceil(a / step) * step; t <= b + 1e-9 * step; t += step) out.push(Math.abs(t) < 1e-9 * step ? 0 : t);
  return out;
}
const fmt = (v) => String(+v.toPrecision(4));

// Length unit for the x axis, from the device's extent.
function lengthUnit(L) {
  if (L >= 1e-3) return { scale: 1e3, name: 'mm' };
  if (L >= 1e-6) return { scale: 1e6, name: 'µm' };
  return { scale: 1e9, name: 'nm' };
}

const style = (id) => {
  const vars = (cols) => cols.map((c, k) => `--driftlet-${k + 1}:${c};`).join('');
  return (
    `<style>` +
    `#${id}{${vars(LIGHT)}--driftlet-ink:#0b0b0b;--driftlet-muted:#52514e;--driftlet-grid:#e4e3df;--driftlet-band:#f1f0ec;--driftlet-face:#a8a7a1}` +
    `@media (prefers-color-scheme: dark){:root:where(:not([data-theme="light"])) #${id}{${vars(DARK)}--driftlet-ink:#ffffff;--driftlet-muted:#c3c2b7;--driftlet-grid:#3a3a37;--driftlet-band:#262624;--driftlet-face:#6b6a64}}` +
    `:root[data-theme="dark"] #${id}{${vars(DARK)}--driftlet-ink:#ffffff;--driftlet-muted:#c3c2b7;--driftlet-grid:#3a3a37;--driftlet-band:#262624;--driftlet-face:#6b6a64}` +
    `</style>`
  );
};

let counter = 0;

/**
 * Line chart of level-diagram traces as an SVG string.
 * @param {ReturnType<typeof traces>} tr from `traces()`
 * @param {{ width?: number, height?: number, ylabel?: string, range?: [number, number], title?: string }} [opts]
 * @returns {string}
 */
export function levelChart(tr, { width = 640, height = 360, ylabel = 'voltage (V)', range, title } = {}) {
  const id = `driftlet-plot-${++counter}`;
  const { x, series, regions } = tr;
  const [y0, y1] = range ?? tr.range;
  const xs = x[0], xe = x[x.length - 1];
  const unit = lengthUnit(xe - xs);

  // Legend entries, flowed into rows above the plot.
  const font = 12, charW = 0.6 * font;
  const rows = [[]];
  let rowW = 0;
  for (const s of series) {
    const w = 28 + s.label.length * charW + 14;
    if (rowW + w > width - 16 && rows[rows.length - 1].length) {
      rows.push([]);
      rowW = 0;
    }
    rows[rows.length - 1].push({ s, x: 8 + rowW });
    rowW += w;
  }
  const legendH = series.length > 1 ? rows.length * 18 + 8 : 0;
  const pad = { l: 56, r: 14, t: legendH + (title ? 22 : 8) + 16, b: 40 };
  const pw = width - pad.l - pad.r, ph = height - pad.t - pad.b;
  const px = (v) => pad.l + ((v - xs) / (xe - xs)) * pw;
  const py = (v) => pad.t + (1 - (v - y0) / (y1 - y0)) * ph;
  const r2 = (v) => Math.round(v * 100) / 100;

  const out = [];
  out.push(
    `<svg id="${id}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" ` +
      `aria-label="${esc(`${title ? `${title}: ` : ''}${series.map((s) => s.label).join(', ')} against x`)}" ` +
      `font-family="system-ui, sans-serif" font-size="${font}">`,
  );
  out.push(style(id));
  if (title) out.push(`<text x="8" y="16" fill="var(--driftlet-ink)" font-weight="600">${esc(title)}</text>`);

  // Region bands (every other one shaded) with their materials, and the faces between them.
  regions.forEach((r, k) => {
    const a = px(r.x0), b = px(r.x1);
    if (k % 2 === 1) out.push(`<rect x="${r2(a)}" y="${pad.t}" width="${r2(b - a)}" height="${ph}" fill="var(--driftlet-band)"/>`);
    if (b - a > r.material.length * charW + 4) out.push(`<text x="${r2((a + b) / 2)}" y="${pad.t - 5}" text-anchor="middle" fill="var(--driftlet-muted)">${esc(r.material)}</text>`);
  });
  // Grid and axes.
  for (const t of ticks(y0, y1, 5)) {
    out.push(`<line x1="${pad.l}" x2="${pad.l + pw}" y1="${r2(py(t))}" y2="${r2(py(t))}" stroke="var(--driftlet-grid)"/>`);
    out.push(`<text x="${pad.l - 6}" y="${r2(py(t) + 4)}" text-anchor="end" fill="var(--driftlet-muted)">${fmt(t)}</text>`);
  }
  for (const t of ticks(xs * unit.scale, xe * unit.scale, 6)) {
    out.push(`<text x="${r2(px(t / unit.scale))}" y="${pad.t + ph + 16}" text-anchor="middle" fill="var(--driftlet-muted)">${fmt(t)}</text>`);
  }
  for (const f of tr.faces) out.push(`<line x1="${r2(px(f))}" x2="${r2(px(f))}" y1="${pad.t}" y2="${pad.t + ph}" stroke="var(--driftlet-face)" stroke-dasharray="2 3"/>`);
  out.push(`<text x="${r2(pad.l + pw / 2)}" y="${height - 8}" text-anchor="middle" fill="var(--driftlet-ink)">x (${unit.name})</text>`);
  out.push(`<text transform="translate(14 ${r2(pad.t + ph / 2)}) rotate(-90)" text-anchor="middle" fill="var(--driftlet-ink)">${esc(ylabel)}</text>`);

  // Lines, clipped to the plot, broken where undefined.
  out.push(`<clipPath id="${id}-clip"><rect x="${pad.l}" y="${pad.t}" width="${pw}" height="${ph}"/></clipPath>`);
  out.push(`<g clip-path="url(#${id}-clip)" fill="none" stroke-width="2" stroke-linejoin="round">`);
  const dash = { level: '', standard: ' stroke-dasharray="6 4"', phi: ' stroke-dasharray="1 3"', redox: ' stroke-dasharray="10 3 2 3"' };
  const colour = (s) => `var(--driftlet-${(s.slot % 8) + 1})`;
  for (const s of series) {
    let d = '', pen = false;
    for (let g = 0; g < x.length; g++) {
      const v = s.y[g];
      if (!Number.isFinite(v)) {
        pen = false;
        continue;
      }
      d += `${pen ? 'L' : 'M'}${r2(px(x[g]))} ${r2(py(v))}`;
      pen = true;
    }
    if (d) out.push(`<path d="${d}" stroke="${colour(s)}"${dash[s.kind]}><title>${esc(s.label)}</title></path>`);
  }
  out.push(`</g>`);
  out.push(`<rect x="${pad.l}" y="${pad.t}" width="${pw}" height="${ph}" fill="none" stroke="var(--driftlet-face)"/>`);

  // Legend (two or more lines): a swatch in each line's colour and dash, its label in ink.
  if (series.length > 1) {
    rows.forEach((row, k) => {
      const yy = (title ? 22 : 0) + 14 + k * 18;
      for (const { s, x: xx } of row) {
        out.push(`<line x1="${xx}" x2="${xx + 22}" y1="${yy - 4}" y2="${yy - 4}" stroke="${colour(s)}" stroke-width="2"${dash[s.kind]}/>`);
        out.push(`<text x="${xx + 28}" y="${yy}" fill="var(--driftlet-ink)">${esc(s.label)}</text>`);
      }
    });
  }
  out.push(`</svg>`);
  return out.join('');
}

/**
 * A solution's level diagram as an SVG string: each charged species' voltage (solid) and
 * standard level (dashed), regions as bands. Options are those of `traces()` and `levelChart()`.
 * @param {import('./types.js').Solution} sol
 * @param {Parameters<typeof traces>[1] & Parameters<typeof levelChart>[1]} [opts]
 * @returns {string}
 */
export function bandDiagram(sol, opts = {}) {
  const { width, height, ylabel, range, title, ...pick } = opts;
  return levelChart(traces(sol, pick), { width, height, ylabel, range, title });
}
