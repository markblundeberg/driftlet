// A real-browser smoke test of the demos: headless Chrome, driven over the DevTools protocol
// with Node's own WebSocket (no dependencies). Each page loads from a local static server, runs
// for a while, then has every control moved to its extremes (sliders to min and max, every
// option of every select, checkboxes toggled, buttons pressed). It fails on any uncaught
// exception, console error or failed request, on a page without plots, or on NaN in a plot.
//
//   npm run smoke [-- pages…]       (CHROME=/path/to/chrome to choose the browser)

import { readdir, mkdtemp, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { root, serve, launch, connect } from './serve.js';

const settle = 2500; // ms a page runs after loading, and after each control is moved

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

// In the page: every control's extremes, as a list of moves to make one at a time.
const controlsScript = `(() => {
  const moves = [];
  document.querySelectorAll('input[type=range]').forEach((el) => moves.push([el.id, 'value', el.min], [el.id, 'value', el.max], [el.id, 'value', el.defaultValue]));
  document.querySelectorAll('select').forEach((el) => {
    for (const o of el.options) moves.push([el.id, 'value', o.value]);
    moves.push([el.id, 'value', el.value]);
  });
  document.querySelectorAll('input[type=checkbox]').forEach((el) => moves.push([el.id, 'checked', !el.checked], [el.id, 'checked', el.checked]));
  document.querySelectorAll('button').forEach((el) => moves.push([el.id, 'click']));
  return moves;
})()`;

const moveScript = ([id, prop, value]) => `(() => {
  const el = document.getElementById(${JSON.stringify(id)});
  if (${JSON.stringify(prop)} === 'click') return el.click();
  el[${JSON.stringify(prop)}] = ${JSON.stringify(value)};
  el.dispatchEvent(new Event('input', { bubbles: true }));
  el.dispatchEvent(new Event('change', { bubbles: true }));
})()`;

// In the page: its plots, and any NaN or Infinity drawn into them.
const plotsScript = `(() => {
  const svgs = [...document.querySelectorAll('svg')].filter((s) => !s.closest('a, nav'));
  const bad = svgs.filter((s) => /NaN|Infinity/.test(s.outerHTML)).map((s) => (s.querySelector('title')?.textContent ?? s.outerHTML.slice(0, 80)));
  return { plots: svgs.length, bad };
})()`;

async function check(cdp, base, page) {
  const problems = [];
  const { targetId } = await cdp.send('Target.createTarget', { url: 'about:blank' });
  const { sessionId } = await cdp.send('Target.attachToTarget', { targetId, flatten: true });
  const listen = (msg) => {
    if (msg.sessionId !== sessionId) return;
    const p = msg.params;
    if (msg.method === 'Runtime.exceptionThrown') {
      const d = p.exceptionDetails;
      problems.push(`exception: ${d.exception?.description ?? d.text}`);
    } else if (msg.method === 'Runtime.consoleAPICalled' && (p.type === 'error' || p.type === 'assert')) {
      problems.push(`console.${p.type}: ${p.args.map((a) => a.value ?? a.description).join(' ')}`);
    } else if (msg.method === 'Log.entryAdded' && p.entry.level === 'error') {
      problems.push(`${p.entry.source}: ${p.entry.text}${p.entry.url ? ` (${p.entry.url})` : ''}`);
    }
  };
  cdp.on(listen);
  const run = (expression) => cdp.send('Runtime.evaluate', { expression, returnByValue: true }, sessionId).then((r) => r.result.value);
  try {
    await cdp.send('Runtime.enable', {}, sessionId);
    await cdp.send('Log.enable', {}, sessionId);
    await cdp.send('Page.enable', {}, sessionId);
    await cdp.send('Emulation.setDeviceMetricsOverride', { width: 1100, height: 900, deviceScaleFactor: 1, mobile: false }, sessionId);
    const loaded = new Promise((ok) => {
      const f = (msg) => msg.sessionId === sessionId && msg.method === 'Page.loadEventFired' && (cdp.off(f), ok());
      cdp.on(f);
    });
    await cdp.send('Page.navigate', { url: `${base}/demos/${page}` }, sessionId);
    await loaded;
    await sleep(settle);
    const moves = await run(controlsScript);
    for (const move of moves) {
      await run(moveScript(move));
      await sleep(move[1] === 'click' ? settle : settle / 2);
    }
    const { plots, bad } = await run(plotsScript);
    if (page !== 'index.html' && plots === 0) problems.push('no plots on the page');
    for (const t of bad) problems.push(`NaN or Infinity drawn in: ${t}`);
    return { problems, moves: moves.length, plots };
  } finally {
    cdp.off(listen);
    await cdp.send('Target.closeTarget', { targetId });
  }
}

const args = process.argv.slice(2);
const pages = args.length > 0 ? args : (await readdir(join(root, 'demos'))).filter((f) => f.endsWith('.html')).sort();
const server = await serve();
const base = `http://127.0.0.1:${server.address().port}`;
const profile = await mkdtemp(join(tmpdir(), 'driftlet-smoke-'));
const { chrome, url } = await launch(profile);
let failed = 0;
try {
  const cdp = await connect(url);
  for (const page of pages) {
    const t0 = Date.now();
    const { problems, moves, plots } = await check(cdp, base, page);
    const secs = ((Date.now() - t0) / 1000).toFixed(1);
    if (problems.length === 0) console.log(`ok   ${page} (${plots} plots, ${moves} control moves, ${secs} s)`);
    else {
      failed++;
      console.log(`FAIL ${page} (${secs} s)`);
      for (const p of [...new Set(problems)]) console.log(`       ${p}`);
    }
  }
  cdp.close();
} finally {
  chrome.kill();
  server.close();
  await sleep(200);
  await rm(profile, { recursive: true, force: true });
}
if (failed > 0) {
  console.log(`${failed} of ${pages.length} pages failed`);
  process.exit(1);
}
