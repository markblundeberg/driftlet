import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join } from 'node:path';

// Every ```js block in the README, the docs and the agent guide runs as a module against the
// source, so examples can't rot. Blocks marked ```js nocheck are skipped (sketches, fragments).
// An ```html page's module script runs against a stub document: it must draw, and redraw when
// its controls move.
const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const entry = (name) => pathToFileURL(join(root, `src/${name}.js`)).href;
const docs = existsSync(join(root, 'docs')) ? readdirSync(join(root, 'docs')) : [];
const files = ['README.md', 'llms.txt', ...docs.filter((f) => f.endsWith('.md')).map((f) => `docs/${f}`)];
const pkg = JSON.parse(readFileSync(join(root, 'package.json'), 'utf8'));

// Imports of the package or its CDN copy, pointed at the source.
const local = (code) =>
  code
    .replaceAll(/from 'driftlet(?:\/(\w+))?'/g, (_, sub) => `from '${entry(sub ?? 'index')}'`)
    .replaceAll(/'https:\/\/cdn\.jsdelivr\.net\/npm\/driftlet@[^/]+\/src\/(\w+)\.js'/g, (_, name) => `'${entry(name)}'`);
const run = (src) => import(`data:text/javascript;base64,${Buffer.from(src).toString('base64')}`);

// A document with the elements a page's script looks up: values, text, HTML, input listeners.
function stubDocument() {
  const elements = new Map();
  const element = (id) => {
    if (!elements.has(id)) {
      const listeners = [];
      elements.set(id, {
        value: '0',
        textContent: '',
        innerHTML: '',
        addEventListener: (type, fn) => listeners.push({ type, fn }),
        fire: (type) => listeners.filter((l) => l.type === type).forEach((l) => l.fn({ type })),
      });
    }
    return elements.get(id);
  };
  return { getElementById: element, element };
}
const until = async (cond, ms = 10000) => {
  const t0 = Date.now();
  while (!cond()) {
    if (Date.now() - t0 > ms) return false;
    await new Promise((r) => setTimeout(r, 5));
  }
  return true;
};

for (const file of files) {
  const text = readFileSync(join(root, file), 'utf8');
  const blocks = [...text.matchAll(/```js\n([\s\S]*?)```/g)].map((m) => m[1]);
  blocks.forEach((code, k) => {
    test(`${file}: example ${k + 1} runs`, async () => {
      const logs = [];
      globalThis.__log = (...a) => logs.push(a);
      await run(local(code).replaceAll('console.log(', '__log('));
      for (const line of logs) for (const v of line) if (typeof v === 'number') assert.ok(Number.isFinite(v), `${file}: ${line}`);
    });
  });
  const pages = [...text.matchAll(/```html\n([\s\S]*?)```/g)].map((m) => m[1]);
  pages.forEach((page, k) => {
    test(`${file}: page ${k + 1} draws, and redraws when its control moves`, async () => {
      const script = /<script type="module">([\s\S]*?)<\/script>/.exec(page)[1];
      const doc = stubDocument();
      globalThis.document = doc;
      try {
        await run(local(script));
        const ids = [...page.matchAll(/id="([^"]+)"/g)].map((m) => m[1]);
        const shown = () => ids.map((id) => doc.element(id).innerHTML).join('');
        assert.ok(await until(() => shown().includes('<svg')), 'it draws an SVG');
        const range = /<input id="([^"]+)" type="range"[^>]*value="([^"]+)"/.exec(page);
        assert.ok(range, 'a page with a range control');
        const el = doc.element(range[1]);
        const first = shown();
        el.value = String(+range[2] + 0.3);
        el.fire('input');
        assert.ok(await until(() => shown() !== first), 'it redraws');
      } finally {
        delete globalThis.document;
      }
    });
  });
}

test('the agent guide pins the CDN version being released', () => {
  const pinned = new Set([...readFileSync(join(root, 'llms.txt'), 'utf8').matchAll(/driftlet@([^/]+)\//g)].map((m) => m[1]));
  assert.equal(pinned.size, 1, [...pinned].join(', '));
  if (pkg.version !== '0.0.0') assert.deepEqual([...pinned], [pkg.version]);
});
