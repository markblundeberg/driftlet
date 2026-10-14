// Shared by the browser tools: the repo over HTTP (module scripts don't load from file://), a
// Chrome to point at it, and a minimal client for its DevTools protocol.

import { spawn } from 'node:child_process';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { join, extname, resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

export const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.svg': 'image/svg+xml', '.png': 'image/png', '.woff2': 'font/woff2', '.txt': 'text/plain', '.json': 'application/json' };

export function findChrome() {
  const candidates = [process.env.CHROME, '/usr/bin/google-chrome', '/usr/bin/google-chrome-stable', '/usr/bin/chromium', '/usr/bin/chromium-browser',
    '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'];
  const found = candidates.find((p) => p && existsSync(p));
  if (!found) throw new Error('no Chrome found; set CHROME=/path/to/chrome');
  return found;
}

// The repo on a free local port; resolves to the server.
export function serve() {
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
export async function launch(profile) {
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
export async function connect(url) {
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
