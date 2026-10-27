import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join } from 'node:path';

// Every demo page is a worked example, so it must keep running: its module script runs here
// against a stand-in DOM (inputs take the values written in the page; canvases accept anything),
// through a few animation frames, and must not throw. Pages with level diagrams must draw them.

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const demos = join(root, 'demos');
const pages = readdirSync(demos).filter((f) => f.endsWith('.html') && f !== 'index.html');

// Anything: every property is another stand-in, calls return one, and assignments are kept.
function anything(record) {
  const kept = {};
  const target = function () {};
  return new Proxy(target, {
    get(_, key) {
      if (key in kept) return kept[key];
      if (key === Symbol.toPrimitive) return (hint) => (hint === 'number' ? 0 : '');
      if (key === 'then') return undefined;
      return (kept[key] = anything(record));
    },
    set(_, key, value) {
      kept[key] = value;
      record?.(key, value);
      return true;
    },
    apply: () => anything(record),
  });
}

for (const page of pages) {
  test(`demos/${page} runs`, async () => {
    const html = readFileSync(join(demos, page), 'utf8');
    const script = /<script type="module">([\s\S]*?)<\/script>/.exec(html)[1];
    const svgs = [];
    const record = (key, value) => key === 'innerHTML' && typeof value === 'string' && value.startsWith('<svg') && svgs.push(value);
    const elements = new Map();
    const element = (id) => {
      if (!elements.has(id)) {
        const el = anything(record);
        const input = new RegExp(`<input id="${id}"[^>]*value="([^"]*)"`).exec(html);
        if (input) el.value = input[1];
        el.addEventListener = () => {};
        el.clientWidth = 600;
        el.clientHeight = 300;
        elements.set(id, el);
      }
      return elements.get(id);
    };
    let frames = 0;
    const errors = [];
    const saved = {};
    const globals = {
      document: Object.assign(anything(record), {
        getElementById: element,
        querySelector: (sel) => (sel === 'script[type="module"]' ? { textContent: script } : anything(record)),
        createElement: () => anything(record),
        head: anything(record),
      }),
      location: { search: '' },
      window: anything(record),
      devicePixelRatio: 1,
      getComputedStyle: () => ({ getPropertyValue: () => '' }),
      requestAnimationFrame: (fn) => {
        if (frames++ < 4) setTimeout(() => {
          try {
            fn();
          } catch (e) {
            errors.push(e);
          }
        }, 0);
      },
      addEventListener: () => {},
    };
    for (const [k, v] of Object.entries(globals)) {
      saved[k] = Object.getOwnPropertyDescriptor(globalThis, k);
      Object.defineProperty(globalThis, k, { value: v, configurable: true, writable: true });
    }
    try {
      const at = (path) => pathToFileURL(join(root, path)).href;
      const src = script
        .replaceAll("'../src/", `'${at('src/')}`)
        .replaceAll("'./plot.js'", `'${at('demos/plot.js')}'`)
        .replaceAll("'./common.js'", `'${at('demos/common.js')}'`);
      await import(`data:text/javascript;base64,${Buffer.from(src).toString('base64')}`);
      // live() and animations finish their work over a few tasks.
      for (let k = 0; k < 40 && (frames < 4 || svgs.length === 0); k++) await new Promise((r) => setTimeout(r, 25));
      await new Promise((r) => setTimeout(r, 50));
    } finally {
      for (const [k, d] of Object.entries(saved)) {
        if (d) Object.defineProperty(globalThis, k, d);
        else delete globalThis[k];
      }
    }
    assert.deepEqual(errors.map((e) => e.stack), []);
    if (/levels\(|profiles\(|levelChart\(/.test(script)) {
      assert.ok(svgs.length > 0, 'it draws its level diagram');
      for (const svg of svgs) assert.doesNotMatch(svg, /NaN/, 'no NaN in the drawing');
    }
  });
}

// The Frumkin page draws its simulation against its own theory and says how well they agree: here
// its script, up to where it needs a page, is run against the claims its text makes.
test("demos/frumkin.html: its theory is the simulation's, as the page says", async () => {
  const html = readFileSync(join(demos, 'frumkin.html'), 'utf8');
  const script = /<script type="module">([\s\S]*?)<\/script>/.exec(html)[1];
  const at = (path) => pathToFileURL(join(root, path)).href;
  const head = script
    .slice(0, script.indexOf('// (All above runs without a page'))
    .replace(/^import .*'\.\/(plot|common)\.js';$/gm, '')
    .replaceAll("'../src/", `'${at('src/')}`);
  assert.ok(head.length > 1000 && !/document\./.test(head), 'the page-free part of its script');
  const { Device, FARADAY } = await import('../src/index.js');
  const { electrode, theory, resistance, opened, ions, reactant, L, VT } = await import(
    `data:text/javascript;base64,${Buffer.from(`${head}\nexport { electrode, theory, resistance, opened, ions, reactant, L, VT };`).toString('base64')}`
  );
  const jlim = (2 * FARADAY * ions['S2O82-'].D * reactant) / L;
  for (const salt of [1, 10, 100]) {
    const dev = new Device(electrode(salt, 0.4));
    let dip = 0, plateau = 0, classic = 0;
    for (let k = 0; k <= 40; k++) {
      const E = 0.4 - 0.05 * k;
      dev.set({ contacts: { left: { V: E } } });
      const sol = dev.solve(), j = -sol.terminals.left.current, th = theory(salt, E);
      assert.ok(sol.converged && sol.warnings.length === 0, `${salt} mM, ${E} V: ${sol.warnings.join(' ')}`);
      if (j < 0.5 * jlim) dip = Math.max(dip, Math.abs(j / th - 1));
      plateau = Math.max(plateau, Math.abs(j - th) / jlim);
      classic = Math.max(classic, Math.abs(Math.log(j / theory(salt, E, 'classic'))));
    }
    // "in the dip they agree with the simulation to 0.2%", "the plateau by about 1%"
    assert.ok(dip < 2e-3 && plateau < 1.2e-2, `${salt} mM: ${dip} in the dip, ${plateau} of the limit`);
    // (and the check can tell: the correction without the layer's resistance is far off, in dilute salt)
    if (salt < 100) assert.ok(classic > Math.log(2), `${salt} mM: the classic correction within ${classic}`);
  }
  // Deep in the dip, the limiting current times L/(L + K); and the plane's level is what K says it is.
  const dev = new Device(electrode(1, 0.4));
  for (let E = 0.4; E > -1.2; E -= 0.05) {
    dev.set({ contacts: { left: { V: E } } });
    dev.solve();
  }
  dev.set({ contacts: { left: { V: -1.2 } } });
  const sol = dev.solve(), lay = opened(sol), K = resistance(lay.c, lay.d.psi / VT), j = -sol.terminals.left.current;
  assert.ok(Math.abs(j / jlim / (L / (L + K)) - 1) < 1e-3, `−1.2 V: ${j / jlim} of the limit, L/(L + K) ${L / (L + K)}`);
  // (the plane's level, as reported, and the page's own K: the edge's level, 1, from the far side)
  const far = lay.u[0] + (lay.N * K) / (ions['S2O82-'].D * lay.c['S2O82-']);
  assert.ok(lay.u[0] < 1e-3 && Math.abs(far - 1) < 1e-6, `the plane's level ${lay.u[0]}, and back out ${far}`);
  assert.ok(Math.abs(lay.u.at(-1) - 1) < 2e-2, `the drawn profile, out at ${lay.x.at(-1)} m: ${lay.u.at(-1)}`);
});
