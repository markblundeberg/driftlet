// Minimal canvas line charts for the demos (no dependencies).
//
// chart(canvas, { series, xlabel, ylabel, xlim, ylim, logx, logy, xscale })
//   series: [{ x, y, color, label, dash, width }]. NaN values break a line; repeated x values
//   (doubled interface nodes) draw as vertical steps. Labels mark subscripts as driftlet/plot's
//   do: `C_ox`, `c_{Zn²⁺}`.

import { labelParts } from '../src/plot.js';

// The same categorical palette as driftlet/plot (style.css defines --c1 … --c8, light and dark).
export const colors = ['var(--c1)', 'var(--c2)', 'var(--c3)', 'var(--c4)', 'var(--c5)', 'var(--c6)', 'var(--c7)', 'var(--c8)'];

const css = (name, fallback) => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
const paint = (c) => (c.startsWith('var(') ? css(c.slice(4, -1), '#888') : c);

// Draws a label with its subscripts lowered and shrunk, aligned as given; returns its width.
function label(ctx, text, x, y, align = 'left') {
  const parts = labelParts(text), font = ctx.font;
  const small = font.replace(/^\d+px/, '9px');
  const widths = parts.map((p) => {
    ctx.font = p.sub ? small : font;
    return ctx.measureText(p.text).width;
  });
  const total = widths.reduce((a, b) => a + b, 0);
  let at = align === 'center' ? x - total / 2 : align === 'right' ? x - total : x;
  const was = ctx.textAlign;
  ctx.textAlign = 'left';
  parts.forEach((p, k) => {
    ctx.font = p.sub ? small : font;
    ctx.fillText(p.text, at, p.sub ? y + 3.6 : y);
    at += widths[k];
  });
  ctx.font = font;
  ctx.textAlign = was;
  return total;
}
const plain = (text) => labelParts(text).map((p) => (p.sub ? `_${p.text}` : p.text)).join('');

export function chart(canvas, { series, xlabel = '', ylabel = '', xlim, ylim, logx = false, logy = false, xscale = 1, legend = true }) {
  canvas.setAttribute('role', 'img');
  canvas.setAttribute('aria-label', `${plain(ylabel)} against ${plain(xlabel)}: ${series.map((s) => s.label).filter(Boolean).map(plain).join(', ')}`);
  const dpr = window.devicePixelRatio || 1;
  const W = canvas.clientWidth, H = canvas.clientHeight;
  if (canvas.width !== Math.round(W * dpr) || canvas.height !== Math.round(H * dpr)) {
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const ink = css('--ink', '#1f2937'), muted = css('--muted', '#6b7280'), grid = css('--grid', '#e5e7eb');
  ctx.clearRect(0, 0, W, H);
  const items = legend ? series.filter((s) => s.label) : [];
  const pad = { l: 64, r: 12, t: items.length ? 30 : 12, b: 42 };
  const fy = logy ? (v) => (v > 0 ? Math.log10(v) : NaN) : (v) => v;
  const fx = logx ? (v) => (v > 0 ? Math.log10(v * xscale) : NaN) : (v) => v * xscale;
  let [x0, x1] = xlim ? (logx ? xlim.map(Math.log10) : xlim) : extent(series.flatMap((s) => Array.from(s.x, fx)));
  let [y0, y1] = ylim ? ylim.map(fy) : extent(series.flatMap((s) => Array.from(s.y, fy)));
  if (!(x1 > x0)) [x0, x1] = [x0 - 1, x0 + 1];
  if (!(y1 > y0)) [y0, y1] = [y0 - 1, y0 + 1];
  if (!ylim) {
    const m = 0.05 * (y1 - y0);
    y0 -= m;
    y1 += m;
  }
  const px = (x) => pad.l + ((x - x0) / (x1 - x0)) * (W - pad.l - pad.r);
  const py = (y) => H - pad.b - ((y - y0) / (y1 - y0)) * (H - pad.t - pad.b);

  // grid and ticks
  ctx.font = '12px system-ui, sans-serif';
  ctx.fillStyle = muted;
  ctx.strokeStyle = grid;
  ctx.lineWidth = 1;
  for (const t of logx ? ticks(Math.ceil(x0), Math.floor(x1), 6).filter(Number.isInteger) : ticks(x0, x1, 6)) {
    ctx.beginPath();
    ctx.moveTo(px(t), pad.t);
    ctx.lineTo(px(t), H - pad.b);
    ctx.stroke();
    ctx.textAlign = 'center';
    ctx.fillText(logx ? `1e${t}` : fmt(t), px(t), H - pad.b + 16);
  }
  for (const t of logy ? ticks(Math.ceil(y0), Math.floor(y1), 6).filter(Number.isInteger) : ticks(y0, y1, 5)) {
    ctx.beginPath();
    ctx.moveTo(pad.l, py(t));
    ctx.lineTo(W - pad.r, py(t));
    ctx.stroke();
    ctx.textAlign = 'right';
    ctx.fillText(logy ? `1e${t}` : fmt(t), pad.l - 6, py(t) + 4);
  }
  ctx.fillStyle = ink;
  label(ctx, xlabel, (pad.l + W - pad.r) / 2, H - 6, 'center');
  ctx.save();
  ctx.translate(14, (pad.t + H - pad.b) / 2);
  ctx.rotate(-Math.PI / 2);
  label(ctx, ylabel, 0, 0, 'center');
  ctx.restore();

  // lines
  ctx.save();
  ctx.beginPath();
  ctx.rect(pad.l, pad.t, W - pad.l - pad.r, H - pad.t - pad.b);
  ctx.clip();
  series.forEach((s, k) => {
    ctx.strokeStyle = paint(s.color ?? colors[k % colors.length]);
    ctx.lineWidth = s.width ?? 2;
    ctx.setLineDash(s.dash ?? []);
    ctx.beginPath();
    let pen = false;
    for (let i = 0; i < s.x.length; i++) {
      const y = fy(s.y[i]);
      if (!Number.isFinite(y)) {
        pen = false;
        continue;
      }
      const X = px(fx(s.x[i])), Y = py(y);
      if (!Number.isFinite(X)) {
        pen = false;
        continue;
      }
      if (pen) ctx.lineTo(X, Y);
      else ctx.moveTo(X, Y);
      pen = true;
    }
    if (s.width !== 0) ctx.stroke();
    if (s.marker) {
      ctx.fillStyle = ctx.strokeStyle;
      for (let i = 0; i < s.x.length; i++) {
        const y = fy(s.y[i]);
        if (Number.isFinite(y)) {
          ctx.beginPath();
          ctx.arc(px(fx(s.x[i])), py(y), 3, 0, 2 * Math.PI);
          ctx.fill();
        }
      }
    }
  });
  ctx.restore();
  ctx.setLineDash([]);

  // legend: one row above the plot
  let xx = pad.l;
  ctx.textAlign = 'left';
  for (const s of items) {
    ctx.strokeStyle = paint(s.color ?? colors[series.indexOf(s) % colors.length]);
    ctx.lineWidth = 2;
    ctx.setLineDash(s.dash ?? []);
    ctx.beginPath();
    ctx.moveTo(xx, 14);
    ctx.lineTo(xx + 20, 14);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = ink;
    xx += 25 + label(ctx, s.label, xx + 25, 18) + 16;
  }
}

function extent(values) {
  let lo = Infinity, hi = -Infinity;
  for (const v of values) {
    if (!Number.isFinite(v)) continue;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return [lo, hi];
}

function ticks(a, b, n) {
  const step0 = (b - a) / n;
  const mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= step0) ?? 10 * mag;
  const out = [];
  for (let t = Math.ceil(a / step) * step; t <= b + 1e-9 * step; t += step) out.push(Math.abs(t) < 1e-12 * step ? 0 : t);
  return out;
}

function fmt(v) {
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v.toExponential(0);
  return String(+v.toPrecision(4));
}

// Initial control values from the URL (?V=0.4), for sharing a state or taking a screenshot.
export function fromQuery(input) {
  const v = new URLSearchParams(location.search).get(input.id);
  if (v !== null) input.value = v;
}
