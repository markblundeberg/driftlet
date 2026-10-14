// The films (films/*.html), rendered to video: each page in headless Chrome at 1080 × 1350,
// drawn one frame at a time on request (window.reel.frame(k)), so a film comes out the same
// every time, however long a frame takes to compute. Frames go straight into ffmpeg.
//
//   node tools/film.js battery [--from k] [--to k] [--png dir [--every n]]    (into films/out/)

import { spawn } from 'node:child_process';
import { mkdtemp, mkdir, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { root, serve, launch, connect } from './serve.js';

const args = process.argv.slice(2);
const name = (args[0] ?? 'battery').replace(/\.html$/, '');
const opt = (k, d) => { const i = args.indexOf(k); return i >= 0 ? args[i + 1] : d; };
const W = 1080, H = 1350;

const server = await serve();
const base = `http://127.0.0.1:${server.address().port}`;
const profile = await mkdtemp(join(tmpdir(), 'driftlet-film-'));
const { chrome, url } = await launch(profile);
const cdp = await connect(url);
const outDir = join(root, 'films', 'out');
await mkdir(outDir, { recursive: true });
try {
  const { targetId } = await cdp.send('Target.createTarget', { url: 'about:blank' });
  const { sessionId } = await cdp.send('Target.attachToTarget', { targetId, flatten: true });
  cdp.on((msg) => {
    if (msg.sessionId !== sessionId) return;
    if (msg.method === 'Runtime.exceptionThrown') console.error('page exception:', msg.params.exceptionDetails.exception?.description);
    if (msg.method === 'Runtime.consoleAPICalled') console.error('page:', msg.params.args.map((a) => a.value ?? a.description).join(' '));
  });
  await cdp.send('Runtime.enable', {}, sessionId);
  await cdp.send('Page.enable', {}, sessionId);
  await cdp.send('Emulation.setDeviceMetricsOverride', { width: W, height: H, deviceScaleFactor: 1, mobile: false }, sessionId);
  await cdp.send('Page.navigate', { url: `${base}/films/${name}.html?film` }, sessionId);
  const run = async (expression) => {
    const r = await cdp.send('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true }, sessionId);
    if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description ?? r.exceptionDetails.text);
    return r.result.value;
  };
  let info;
  for (let k = 0; k < 600 && !(info = await run('window.reel && { frames: reel.frames, fps: reel.fps }')); k++) await new Promise((r) => setTimeout(r, 100));
  if (!info) throw new Error('the page never set window.reel');
  console.error('film:', JSON.stringify(info));
  const from = +opt('--from', 0), to = Math.min(+opt('--to', info.frames - 1), info.frames - 1), pngs = opt('--png');
  const mp4 = join(outDir, `${name}.mp4`);
  const ff = pngs ? null : spawn('ffmpeg', ['-y', '-loglevel', 'error', '-f', 'image2pipe', '-framerate', String(info.fps), '-i', '-',
    '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '18', '-preset', 'slow', '-movflags', '+faststart', mp4], { stdio: ['pipe', 'inherit', 'inherit'] });
  if (pngs) await mkdir(pngs, { recursive: true });
  const t0 = Date.now();
  const every = +opt('--every', 1); // stills only: every nth frame
  for (let k = from; k <= to; k++) {
    await run(`reel.frame(${k})`);
    if (pngs && (k - from) % every) continue;
    const { data } = await cdp.send('Page.captureScreenshot', { format: 'png', clip: { x: 0, y: 0, width: W, height: H, scale: 1 } }, sessionId);
    const buf = Buffer.from(data, 'base64');
    if (pngs) await writeFile(join(pngs, `${String(k).padStart(4, '0')}.png`), buf);
    else if (!ff.stdin.write(buf)) await new Promise((r) => ff.stdin.once('drain', r));
    if (k % 30 === 0) process.stderr.write(`\rframe ${k}/${to} (${((Date.now() - t0) / 1000).toFixed(0)} s)`);
  }
  process.stderr.write('\n');
  if (ff) {
    ff.stdin.end();
    await new Promise((ok, fail) => ff.on('exit', (c) => (c === 0 ? ok() : fail(new Error(`ffmpeg exited ${c}`)))));
    console.log(mp4);
  } else console.log(pngs);
} finally {
  cdp.close();
  chrome.kill();
  server.close();
  await new Promise((r) => chrome.once('exit', r));
  await rm(profile, { recursive: true, force: true, maxRetries: 10, retryDelay: 200 });
}
