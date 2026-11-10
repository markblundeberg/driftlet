import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join } from 'node:path';
import { figureBlocks, drawFigure } from '../tools/figures.js';

// Every ```js block in the README, the docs and the agent guide runs as a module against the
// source, so examples can't rot. Blocks marked ```js nocheck are skipped (sketches, fragments).
// A logging line's comment that starts `// =>` is a claim about what it prints, and is held to
// it: each item listed (commas between; a bracketed remark, or anything after a colon or a
// semicolon, is commentary) must appear in that line's output. (Two of the templates' numbers
// had drifted from what they printed, 21.12 for 21.10 and 18.8 for 18.7, with nothing to notice.)
// An ```html page's module script runs against a stub document: it must draw, and redraw when
// its controls move. A block fenced ```js figure=NAME also draws the picture shown under it,
// docs/figures/NAME.svg, and the file is held to what the block draws (tools/figures.js).
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
  const blocks = [...text.matchAll(/```js(?: figure=[\w-]+)?\n([\s\S]*?)```/g)].map((m) => m[1]);
  blocks.forEach((code, k) => {
    test(`${file}: example ${k + 1} runs, and prints what its comments say`, async () => {
      const logs = [], printed = new Map();
      globalThis.__log = (i, ...a) => {
        logs.push(a);
        printed.set(i, `${printed.get(i) ?? ''}${a.map((v) => (typeof v === 'string' ? v : JSON.stringify(v))).join(' ')}\n`);
      };
      const lines = code.split('\n');
      await run(local(lines.map((l, i) => l.replaceAll('console.log(', `__log(${i}, `)).join('\n')));
      for (const line of logs) for (const v of line) if (typeof v === 'number') assert.ok(Number.isFinite(v), `${file}: ${line}`);
      lines.forEach((l, i) => {
        const claim = /\/\/ => ?(.*)$/.exec(l)?.[1];
        if (claim === undefined) return;
        assert.ok(l.includes('console.log('), `${file}: a "// =>" comment belongs on the line that logs: ${l.trim()}`);
        const items = claim.replace(/\([^)]*\)/g, '').split(/[:;]/)[0].split(',').map((x) => x.trim().replaceAll('−', '-')).filter(Boolean);
        assert.ok(items.length > 0, `${file}: "// =>" with nothing to check: ${l.trim()}`);
        const out = (printed.get(i) ?? '').replaceAll('−', '-');
        for (const item of items) assert.ok(out.includes(item), `${file}: the comment says ${item}, the line printed ${JSON.stringify(out.trim())}: ${l.trim()}`);
      });
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

// Two drawings of one figure: the same text, and every number within `tol` (the coordinates are
// rounded to 0.01 px, so a last-bit difference between two builds of Node can move one by a step).
function sameDrawing(a, b, tol = 0.011) {
  const number = /-?\d+(?:\.\d+)?/g;
  if (a.replace(number, '#') !== b.replace(number, '#')) return false;
  const na = a.match(number) ?? [], nb = b.match(number) ?? [];
  return na.every((v, k) => Math.abs(v - nb[k]) <= tol);
}
const figures = figureBlocks();
for (const block of figures) {
  test(`${block.file}: figure ${block.name} is what its example draws, and the page shows it`, async () => {
    const path = join(root, 'docs/figures', `${block.name}.svg`);
    assert.ok(existsSync(path), `no docs/figures/${block.name}.svg: npm run figures`);
    assert.ok(sameDrawing(readFileSync(path, 'utf8'), await drawFigure(block)), `docs/figures/${block.name}.svg isn't what its example in ${block.file} draws now: npm run figures, and look at what changed`);
    assert.ok(readFileSync(join(root, block.file), 'utf8').includes(`](figures/${block.name}.svg)`), `${block.file} doesn't show figures/${block.name}.svg`);
  });
}
test('every figure in docs/figures is drawn by one example', () => {
  const names = figures.map((b) => b.name);
  assert.deepEqual([...new Set(names)], names, 'two examples draw a figure of the same name');
  const files = existsSync(join(root, 'docs/figures')) ? readdirSync(join(root, 'docs/figures')) : [];
  assert.deepEqual(files.sort(), names.map((n) => `${n}.svg`).sort());
});

test('the agent guide pins the CDN version being released', () => {
  const pinned = new Set([...readFileSync(join(root, 'llms.txt'), 'utf8').matchAll(/driftlet@([^/]+)\//g)].map((m) => m[1]));
  assert.equal(pinned.size, 1, [...pinned].join(', '));
  if (pkg.version !== '0.0.0') assert.deepEqual([...pinned], [pkg.version]);
});

test("the docs' inline math is $`…`$, which markdown's emphasis rules can't break", () => {
  for (const file of [...docs.filter((f) => f.endsWith('.md')).map((f) => `docs/${f}`), 'ROADMAP.md', 'CHANGELOG.md']) {
    const prose = readFileSync(join(root, file), 'utf8').replace(/```[\s\S]*?```/g, '');
    prose.split('\n').forEach((line, k) => {
      const bare = line.replace(/\$`[^`]*`\$/g, '').replace(/`[^`]*`/g, '');
      assert.ok(!/(?<!\\)\$/.test(bare), `${file}:${k + 1}: a bare $ (write inline math as $\`…\`$): ${line.trim()}`);
    });
  }
});
