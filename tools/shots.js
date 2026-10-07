// The demos' screenshots (the gallery's thumbnails and the README's), retaken: each page with
// ?shot (title and figures only) and the state it's shown in, in headless Chrome on virtual
// time, so an animation is caught at the same moment every run.
//
//   npm run shots [-- pages…]       (into demos/screenshots/; CHROME=/path/to/chrome)

import { execFile } from 'node:child_process';
import { promisify } from 'node:util';
import { mkdtemp, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { findChrome, root, serve } from './serve.js';

// Each page's query (its state: a time to stop at, a slider) and how long it runs, virtually.
const SHOTS = {
  pn: { query: '', budget: 60000 },
  solar: { query: '', budget: 60000 },
  organic: { query: 'LX=25', budget: 60000 },
  perovskite: { query: 't=0.75', budget: 60000 },
  mos: { query: '', budget: 60000 },
  redox: { query: 't=15', budget: 60000 },
  daniell: { query: 't=3600', budget: 60000 },
  saturation: { query: 't=0.53', budget: 60000 },
  impedance: { query: 'f=1.25&phase=15', budget: 60000 },
  cell: { query: '', budget: 60000 },
  axon: { query: 't=72.2', budget: 60000 },
  channel: { query: '', budget: 60000 },
  junction: { query: 't=5', budget: 60000 },
  membrane: { query: '', budget: 60000 },
  'double-layer': { query: '', budget: 60000 },
  insertion: { query: 't=6000', budget: 60000 },
};

const args = process.argv.slice(2);
const names = args.length > 0 ? args.map((a) => a.replace(/\.html$/, '')) : Object.keys(SHOTS);
for (const n of names) if (!SHOTS[n]) throw new Error(`no screenshot settings for '${n}' (known: ${Object.keys(SHOTS).join(', ')})`);
const server = await serve();
const base = `http://127.0.0.1:${server.address().port}`;
const chrome = findChrome();
try {
  for (const n of names) {
    const { query, budget } = SHOTS[n];
    const profile = await mkdtemp(join(tmpdir(), 'driftlet-shot-'));
    const out = join(root, 'demos', 'screenshots', `${n}.png`);
    await promisify(execFile)(chrome, ['--headless', '--disable-gpu', '--hide-scrollbars', `--user-data-dir=${profile}`, `--virtual-time-budget=${budget}`,
      `--screenshot=${out}`, '--window-size=1000,640', `${base}/demos/${n}.html?${query ? `${query}&` : ''}shot`], { timeout: 300000 });
    await rm(profile, { recursive: true, force: true });
    console.log(`demos/screenshots/${n}.png`);
  }
} finally {
  server.close();
}
