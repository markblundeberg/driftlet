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
    if (/levels\(|relativeLevels\(/.test(script)) {
      assert.ok(svgs.length > 0, 'it draws its level diagram');
      for (const svg of svgs) assert.doesNotMatch(svg, /NaN/, 'no NaN in the drawing');
    }
  });
}
