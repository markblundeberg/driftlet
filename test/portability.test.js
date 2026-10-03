import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { Worker } from 'node:worker_threads';
import { Device, FARADAY, units } from '../src/index.js';

// driftlet runs unchanged in a browser, a Web Worker, or from a CDN: plain ES modules with
// relative imports, and nothing from Node or the DOM. A browser isn't available to node --test,
// so these check what makes that work.

const srcDir = new URL('../src/', import.meta.url);
const sources = readdirSync(srcDir).filter((f) => f.endsWith('.js'));
// Code only: comments and string literals stripped (roughly, which is enough here).
const code = (f) =>
  readFileSync(new URL(f, srcDir), 'utf8')
    .replace(/\/\*[\s\S]*?\*\//g, '')
    .replace(/(^|[^:])\/\/.*$/gm, '$1')
    .replace(/`(?:\\.|[^`\\])*`|'(?:\\.|[^'\\])*'/g, "''");

test('every import is a relative path to a .js file (what a browser and a CDN need)', () => {
  for (const f of sources) {
    const text = readFileSync(new URL(f, srcDir), 'utf8');
    for (const m of text.matchAll(/^\s*(?:import|export)\b[^;]*?\bfrom\s+'([^']+)'/gm)) {
      assert.match(m[1], /^\.\/[\w-]+\.js$/, `${f}: import from '${m[1]}'`);
    }
    assert.doesNotMatch(text, /\bimport\s*\(\s*'[^.]/, `${f}: dynamic import of a bare specifier`);
  }
});

test('no Node or DOM globals in the library', () => {
  for (const f of sources) {
    const c = code(f);
    for (const name of ['process', 'require', 'Buffer', '__dirname', 'window', 'document', 'self']) {
      assert.doesNotMatch(c, new RegExp(`\\b${name}(\\s*[.(]|\\[)`), `${f} uses ${name}`);
    }
  }
});

test('no ** or Math.pow in the library: their last bit differs between engine versions', () => {
  // V8's changed between Node 22 and 24, enough to flip adaptive time steps; src/pow.js has
  // versions that agree everywhere (repeated squaring, a square root, exp and log).
  for (const f of sources) {
    const found = code(f).match(/\*\*|Math\.pow/);
    assert.equal(found, null, `${f} uses ${found?.[0]}: use powi or powr from pow.js`);
  }
});

test('a definition posted to a worker solves there exactly as here', async () => {
  const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
  const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
  const def = {
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } } },
    regions: [
      { material: 'Si', length: 1e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'Si', length: 1e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
    ],
    contacts: { left: ohmic(0), right: ohmic(0.4) },
    bulkReactions: [{ nu: { 'e-': -1, 'h+': -1 }, kf: { Si: 1e-6 } }],
    grid: { hmin: 1e-9, hmax: 50e-9 },
  };
  const here = new Device(def).solve();
  const worker = new Worker(new URL('./fixtures/solve-worker.js', import.meta.url), { workerData: { def } });
  const there = await new Promise((resolve, reject) => {
    worker.once('message', resolve);
    worker.once('error', reject);
  });
  await worker.terminate();
  assert.ok(there.converged);
  assert.equal(there.current, here.current);
  assert.deepEqual(Array.from(there.phi), Array.from(here.phi));
});
