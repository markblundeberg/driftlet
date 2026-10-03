// A real-browser smoke test of the demos: headless Chrome, driven over the DevTools protocol
// with Node's own WebSocket (no dependencies). Each page loads from a local static server, runs
// for a while, then has every control moved to its extremes (sliders to min and max, every
// option of every select, checkboxes toggled, buttons pressed). It fails on any uncaught
// exception, console error or failed request, on a page without plots, or on NaN in a plot.
//
//   npm run smoke [-- pages…]       (CHROME=/path/to/chrome to choose the browser)

import { spawn } from 'node:child_process';
import { createServer } from 'node:http';
import { readFile, readdir, mkdtemp, rm } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, extname, resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const settle = 2500; // ms a page runs after loading, and after each control is moved
const types = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.svg': 'image/svg+xml', '.png': 'image/png', '.json': 'application/json' };

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

function findChrome() {
  const candidates = [process.env.CHROME, '/usr/bin/google-chrome', '/usr/bin/google-chrome-stable', '/usr/bin/chromium', '/usr/bin/chromium-browser',
    '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'];
  const found = candidates.find((p) => p && existsSync(p));
  if (!found) throw new Error('no Chrome found; set CHROME=/path/to/chrome');
  return found;
}

// The repo, served over HTTP (module scripts don't load from file://).
function serve() {
  const server = createServer(async (req, res) => {
    const path = join(root, decodeURIComponent(new URL(req.url, 'http://x').pathname));
    if (!path.startsWith(root)) return res.writeHead(403).end();
    try {
      const body = await readFile(path);
      res.writeHead(200, { 'content-type': types[extname(path)] ?? 'application/octet-stream' }).end(body);
    } catch {
      res.writeHead(404).end();
    }
  });
  return new Promise((ok) => server.listen(0, '127.0.0.1', () => ok(server)));
}

// Chrome with a debugging port; resolves to its browser-wide WebSocket URL.
async function launch(profile) {
  const chrome = spawn(findChrome(), ['--headless=new', '--remote-debugging-port=0', `--user-data-dir=${profile}`, '--no-first-run',
    '--no-default-browser-check', '--disable-gpu', '--disable-extensions', '--no-sandbox', 'about:blank'], { stdio: ['ignore', 'ignore', 'pipe'] });
  const url = await new Promise((ok, fail) => {
    let err = '';
    chrome.stderr.on('data', (d) => {
      err += d;
      const m = err.match(/DevTools listening on (ws:\/\/\S+)/);
      if (m) ok(m[1]);
    });
    chrome.on('exit', (code) => fail(new Error(`Chrome exited (${code}) before listening:\n${err}`)));
  });
  return { chrome, url };
}

// A minimal DevTools-protocol client: send(method, params, sessionId), and events by session.
async function connect(url) {
  const ws = new WebSocket(url);
  await new Promise((ok, fail) => {
    ws.onopen = ok;
    ws.onerror = () => fail(new Error(`can't connect to ${url}`));
  });
  let id = 0;
  const pending = new Map(), listeners = new Set();
  ws.onmessage = ({ data }) => {
    const msg = JSON.parse(data);
    if (msg.id !== undefined) {
      const p = pending.get(msg.id);
      pending.delete(msg.id);
      if (msg.error) p.fail(new Error(`${p.method}: ${msg.error.message}`));
      else p.ok(msg.result);
    } else for (const f of listeners) f(msg);
  };
  return {
    send: (method, params = {}, sessionId) =>
      new Promise((ok, fail) => {
        pending.set(++id, { ok, fail, method });
        ws.send(JSON.stringify({ id, method, params, sessionId }));
      }),
    on: (f) => listeners.add(f),
    off: (f) => listeners.delete(f),
    close: () => ws.close(),
  };
}

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
