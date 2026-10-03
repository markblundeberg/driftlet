// driftlet/plot: a level diagram as an SVG string, with no dependencies and no DOM, so it
// works the same in a page, a Web Worker or Node. Put it in a page with
// `element.innerHTML = bandDiagram(sol)`. For your own plotting, `traces()` in driftlet/kit gives
// the same lines as data.
//
// The theme (explained in docs/visualization.md): species voltages V_i are thick solid lines and standard levels
// V°_i thin solid ones, in the species' colour; redox levels are thick dashed lines and their
// standard levels thin dashed ones; φ is thin, dotted and grey. Electrons are steel blue, cations
// (holes too) warm colours, anions cool ones, redox levels blue-violets; each colour follows its
// species. A species drawn with a display offset carries a ⌇ mark. The colours are CSS custom
// properties (--driftlet-electron, --driftlet-cation-1 …, --driftlet-anion-1 …,
// --driftlet-redox-1 …, --driftlet-ink, …) with light and dark defaults of no specificity, so any
// rule on the page restyles them.

import { traces } from './traces.js';
import { powi } from './pow.js';

// Colours by role, light and dark (validated together for colour-vision deficiency; line weight,
// dash and the legend carry identity too). A role with more members than colours cycles.
export const THEME = Object.freeze({
  light: { electron: ['#3a74b4'], cation: ['#d1362f', '#c8650f', '#b8417a', '#8a5a12'], anion: ['#11968a', '#0a9bbd', '#2e8b57'], redox: ['#5a4bc6', '#9a3fae', '#3f5fa8'] },
  dark: { electron: ['#5b95d6'], cation: ['#e8584f', '#d0772a', '#d65d95', '#b88a3e'], anion: ['#219c8f', '#2aa5bb', '#3aa874'], redox: ['#8b7ff0', '#a86bc9', '#7a9ae0'] },
});
const ROLES = ['electron', 'cation', 'anion', 'redox'];

/**
 * The colour of a trace's role and slot (a CSS colour), for drawing other charts in the same
 * colours as the level diagrams.
 * @param {'electron' | 'cation' | 'anion' | 'redox'} role
 * @param {number} slot
 * @param {{ dark?: boolean }} [opts]
 */
export function themeColor(role, slot, { dark = false } = {}) {
  const list = THEME[dark ? 'dark' : 'light'][role];
  return list[slot % list.length];
}

// Line weight and dash by kind (thick: what carriers feel; thin: standard levels).
const STROKE = {
  level: { width: 2.8, dash: '' },
  standard: { width: 1.3, dash: '' },
  redox: { width: 2.8, dash: '8 5' },
  'redox-standard': { width: 1.3, dash: '6 4' },
  phi: { width: 1.3, dash: '1.5 3' },
};

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

/**
 * A label split into plain and subscript runs. Labels mark subscripts TeX-style: `V_{e⁻}`, or
 * `C_ox` for a run of letters, digits and charge signs. Everything else is literal.
 * @param {string} label
 * @returns {{ text: string, sub: boolean }[]}
 */
export function labelParts(label) {
  const out = [];
  const re = /_(?:\{([^}]*)\}|([\p{L}\p{N}⁺⁻]+))/gu;
  let at = 0;
  for (const m of String(label).matchAll(re)) {
    if (m.index > at) out.push({ text: label.slice(at, m.index), sub: false });
    out.push({ text: m[1] ?? m[2], sub: true });
    at = m.index + m[0].length;
  }
  if (at < label.length) out.push({ text: label.slice(at), sub: false });
  return out;
}
// The label as plain text (for titles and accessible names): braces dropped, `_` kept.
const plain = (label) => labelParts(label).map((p) => (p.sub ? `_${p.text}` : p.text)).join('');
// The label as SVG text content, subscripts lowered and shrunk (`dy`, which every browser
// supports, rather than baseline-shift).
const rich = (label, font) => {
  const dy = Math.round(0.3 * font * 10) / 10;
  let low = false;
  return labelParts(label)
    .map((p) => {
      const shift = p.sub === low ? '' : ` dy="${p.sub ? dy : -dy}"`;
      low = p.sub;
      return p.sub ? `<tspan${shift} font-size="${0.75 * font}">${esc(p.text)}</tspan>` : shift ? `<tspan${shift}>${esc(p.text)}</tspan>` : esc(p.text);
    })
    .join('');
};
// Its width in characters, for layout.
const span = (label) => labelParts(label).reduce((n, p) => n + (p.sub ? 0.75 : 1) * [...p.text].length, 0);

function ticks(a, b, n) {
  const step0 = (b - a) / n;
  const mag = powi(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= step0) ?? 10 * mag;
  const out = [];
  for (let t = Math.ceil(a / step) * step; t <= b + 1e-9 * step; t += step) out.push(Math.abs(t) < 1e-9 * step ? 0 : t);
  return out;
}
const fmt = (v) => String(+v.toPrecision(4));

// The range of the lines within [xs, xe], padded.
function fit(tr, xs, xe) {
  let lo = Infinity, hi = -Infinity;
  for (const s of tr.series) {
    for (let g = 0; g < tr.x.length; g++) {
      if (tr.x[g] < xs || tr.x[g] > xe || !Number.isFinite(s.y[g])) continue;
      lo = Math.min(lo, s.y[g]);
      hi = Math.max(hi, s.y[g]);
    }
  }
  if (!(hi > lo)) return tr.range;
  const pad = 0.06 * (hi - lo);
  return [lo - pad, hi + pad];
}

// Length unit for the x axis, from the device's extent.
function lengthUnit(L) {
  if (L >= 1e-3) return { scale: 1e3, name: 'mm' };
  if (L >= 1e-6) return { scale: 1e6, name: 'µm' };
  return { scale: 1e9, name: 'nm' };
}

const style = (id) => {
  const vars = (mode) => ROLES.flatMap((role) => THEME[mode][role].map((c, k) => `--driftlet-${role}${role === 'electron' ? '' : `-${k + 1}`}:${c};`)).join('');
  const light = `${vars('light')}--driftlet-ink:#0b0b0b;--driftlet-muted:#52514e;--driftlet-grid:#e4e3df;--driftlet-band:#f1f0ec;--driftlet-face:#a8a7a1`;
  const dark = `${vars('dark')}--driftlet-ink:#ffffff;--driftlet-muted:#c3c2b7;--driftlet-grid:#3a3a37;--driftlet-band:#262624;--driftlet-face:#6b6a64`;
  return (
    `<style>` +
    `:where(#${id}){${light}}` +
    `@media (prefers-color-scheme: dark){:where(:root:not([data-theme="light"]) #${id}){${dark}}}` +
    `:where(:root[data-theme="dark"] #${id}){${dark}}` +
    `</style>`
  );
};

let counter = 0;
// Ids unique on a page, even with two copies of this module loaded (a CDN one and a local one).
const instance = Math.random().toString(36).slice(2, 8);

/**
 * Line chart of level-diagram traces as an SVG string.
 * @param {ReturnType<typeof traces>} tr from `traces()`
 * @param {{ width?: number, height?: number, ylabel?: string, range?: [number, number], xlim?: [number, number], title?: string, ytick?: (v: number) => string }} [opts]
 *   `xlim` (m) zooms into part of the device; the range then fits what's shown. `ytick` writes
 *   the y axis's tick labels, say as powers of ten for lines that are log₁₀ of something
 * @returns {string}
 */
export function levelChart(tr, { width = 640, height = 360, ylabel = 'voltage (V)', range, xlim, title, ytick = fmt } = {}) {
  const id = `driftlet-plot-${instance}-${++counter}`;
  const { x, series, regions } = tr;
  const xs = xlim ? xlim[0] : x[0], xe = xlim ? xlim[1] : x[x.length - 1];
  const [y0, y1] = range ?? (xlim ? fit(tr, xs, xe) : tr.range);
  const unit = lengthUnit(xe - xs);

  // Legend entries, flowed into rows above the plot.
  const font = 12, charW = 0.6 * font;
  const rows = [[]];
  let rowW = 0;
  for (const s of series) {
    const w = 28 + span(s.label) * charW + 14;
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
    `<svg id="${id}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" style="max-width:100%;height:auto" role="img" ` +
      `aria-label="${esc(`${title ? `${plain(title)}: ` : ''}${series.map((s) => plain(s.label)).join(', ')} against x`)}" ` +
      `font-family="system-ui, sans-serif" font-size="${font}">`,
  );
  out.push(style(id));
  if (title) out.push(`<text x="8" y="16" fill="var(--driftlet-ink)" font-weight="600">${rich(title, font)}</text>`);

  // Region bands (every other one shaded) with their materials, and the faces between them.
  regions.forEach((r, k) => {
    if (r.x1 <= xs || r.x0 >= xe) return;
    const a = px(Math.max(r.x0, xs)), b = px(Math.min(r.x1, xe));
    if (k % 2 === 1) out.push(`<rect x="${r2(a)}" y="${pad.t}" width="${r2(b - a)}" height="${ph}" fill="var(--driftlet-band)"/>`);
    // A region's own name where it was given one (not the default 'region k'), with its material.
    const label = r.name && !/^region \d+$/.test(r.name) && r.name !== r.material ? `${r.name} (${r.material})` : r.material;
    const named = r.name && !/^region \d+$/.test(r.name) ? r.name : r.material;
    const shown = [label, named, r.material].find((t) => b - a > t.length * charW + 4) ?? null;
    if (shown) out.push(`<text x="${r2((a + b) / 2)}" y="${pad.t - 5}" text-anchor="middle" fill="var(--driftlet-muted)">${esc(shown)}</text>`);
  });
  // Grid and axes.
  for (const t of ticks(y0, y1, 5)) {
    out.push(`<line x1="${pad.l}" x2="${pad.l + pw}" y1="${r2(py(t))}" y2="${r2(py(t))}" stroke="var(--driftlet-grid)"/>`);
    out.push(`<text x="${pad.l - 6}" y="${r2(py(t) + 4)}" text-anchor="end" fill="var(--driftlet-muted)">${esc(ytick(t))}</text>`);
  }
  for (const t of ticks(xs * unit.scale, xe * unit.scale, 6)) {
    out.push(`<text x="${r2(px(t / unit.scale))}" y="${pad.t + ph + 16}" text-anchor="middle" fill="var(--driftlet-muted)">${fmt(t)}</text>`);
  }
  for (const f of tr.faces.filter((f) => f > xs && f < xe)) out.push(`<line x1="${r2(px(f))}" x2="${r2(px(f))}" y1="${pad.t}" y2="${pad.t + ph}" stroke="var(--driftlet-face)" stroke-dasharray="2 3"/>`);
  out.push(`<text x="${r2(pad.l + pw / 2)}" y="${height - 8}" text-anchor="middle" fill="var(--driftlet-ink)">x (${unit.name})</text>`);
  out.push(`<text transform="translate(14 ${r2(pad.t + ph / 2)}) rotate(-90)" text-anchor="middle" fill="var(--driftlet-ink)">${rich(ylabel, font)}</text>`);

  // Lines, clipped to the plot, broken where undefined.
  out.push(`<clipPath id="${id}-clip"><rect x="${pad.l}" y="${pad.t}" width="${pw}" height="${ph}"/></clipPath>`);
  out.push(`<g clip-path="url(#${id}-clip)" fill="none" stroke-linejoin="round">`);
  const colour = (s) => {
    if (s.role === 'phi' || !s.role) return 'var(--driftlet-muted)';
    if (s.role === 'electron') return 'var(--driftlet-electron)';
    return `var(--driftlet-${s.role}-${(s.slot % THEME.light[s.role].length) + 1})`;
  };
  const stroke = (s) => {
    const st = STROKE[s.kind] ?? STROKE.level;
    return `stroke="${colour(s)}" stroke-width="${st.width}"${st.dash ? ` stroke-dasharray="${st.dash}"` : ''}`;
  };
  // Thin lines first, so the thick ones the carriers feel sit on top.
  const order = [...series].sort((a, b) => (STROKE[a.kind]?.width ?? 2) - (STROKE[b.kind]?.width ?? 2));
  const marks = [];
  for (const s of order) {
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
    if (d) out.push(`<path d="${d}" ${stroke(s)}><title>${esc(plain(s.label))}</title></path>`);
    // A species drawn with a display offset: a small zigzag across its V_i, near the left.
    if (s.shift && s.kind === 'level') {
      const g = [...x.keys()].find((k) => x[k] >= xs + 0.04 * (xe - xs) && Number.isFinite(s.y[k]) && s.y[k] >= y0 && s.y[k] <= y1);
      if (g !== undefined) marks.push({ x: px(x[g]), y: py(s.y[g]), s });
    }
  }
  for (const m of marks) {
    const zz = `M${r2(m.x - 3)} ${r2(m.y - 9)}l6 4.5l-6 4.5l6 4.5l-6 4.5`;
    out.push(`<path d="${zz}" stroke="var(--driftlet-ink)" stroke-width="1.6" fill="none"><title>${esc(plain(m.s.label))}</title></path>`);
  }
  out.push(`</g>`);
  if (series.some((s) => s.shift)) {
    out.push(`<text x="${pad.l + pw - 6}" y="${pad.t + ph - 6}" text-anchor="end" fill="var(--driftlet-muted)">⌇ = per-species offset</text>`);
  }
  out.push(`<rect x="${pad.l}" y="${pad.t}" width="${pw}" height="${ph}" fill="none" stroke="var(--driftlet-face)"/>`);

  // Legend (two or more lines): a swatch in each line's colour and dash, its label in ink.
  if (series.length > 1) {
    rows.forEach((row, k) => {
      const yy = (title ? 22 : 0) + 14 + k * 18;
      for (const { s, x: xx } of row) {
        out.push(`<line x1="${xx}" x2="${xx + 22}" y1="${yy - 4}" y2="${yy - 4}" ${stroke(s)}/>`);
        out.push(`<text x="${xx + 28}" y="${yy}" fill="var(--driftlet-ink)">${rich(s.label, font)}</text>`);
      }
    });
  }
  out.push(`</svg>`);
  return out.join('');
}

/**
 * A solution's level diagram as an SVG string: each charged species' voltage (thick) and
 * standard level (thin), regions as bands; or with `energy: true`, the familiar energy-up band
 * diagram of electrons and holes. Options are those of `traces()` and `levelChart()`.
 * @param {import('./types.js').Solution} sol
 * @param {Parameters<typeof traces>[1] & Parameters<typeof levelChart>[1]} [opts]
 * @returns {string}
 */
export function bandDiagram(sol, opts = {}) {
  const { width, height, ylabel = opts.energy ? 'electron energy (eV)' : undefined, range, xlim, title, ytick, ...pick } = opts;
  return levelChart(traces(sol, pick), { width, height, ylabel, range, xlim, title, ytick });
}
