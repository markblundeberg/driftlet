// The docs' figures, drawn by the docs' own examples: a block fenced ```js figure=NAME is run,
// and the level diagram it draws (its one call of bandDiagram() or levelChart()) is the picture
// at docs/figures/NAME.svg. test/docs.test.js redraws each and holds the file to it, so a figure
// can't drift from the code beside it.
//
//   npm run figures       (writes docs/figures/; then look at what changed before committing)

import { readFileSync, readdirSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join, resolve } from 'node:path';

export const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const entry = (name) => pathToFileURL(join(root, `src/${name}.js`)).href;
const moduleUrl = (src) => `data:text/javascript;base64,${Buffer.from(src).toString('base64')}`;

// driftlet/plot with its two drawing calls recorded.
const recorder = moduleUrl(
  `export * from '${entry('plot')}';\n` +
    `import { bandDiagram as band, levelChart as chart } from '${entry('plot')}';\n` +
    `const keep = (draw) => (...args) => { const svg = draw(...args); globalThis.__figures.push(svg); return svg; };\n` +
    `export const bandDiagram = keep(band), levelChart = keep(chart);\n`,
);

/** The figure blocks of the docs: `{ file, name, code }`, in order. */
export function figureBlocks() {
  const out = [];
  for (const f of readdirSync(join(root, 'docs')).filter((f) => f.endsWith('.md')).sort()) {
    const text = readFileSync(join(root, 'docs', f), 'utf8');
    for (const m of text.matchAll(/```js figure=([\w-]+)\n([\s\S]*?)```/g)) out.push({ file: `docs/${f}`, name: m[1], code: m[2] });
  }
  return out;
}

// A page gives a diagram its background and its theme. A file on its own, shown as an image on
// any page, light or dark, needs its own: the demos' panel colour, under the light theme (the
// dark rules look for a page around the diagram, and a file has none). Its id is the figure's
// name, where a page's is random.
const standalone = (svg, name) =>
  svg
    .replaceAll(/driftlet-plot-[a-z0-9]+-\d+/g, `driftlet-figure-${name}`)
    .replace('</style>', '</style><rect width="100%" height="100%" rx="6" fill="#fdfbf7"/>');

let runs = 0;
/** Run a block and return the picture it draws, as the file's text. */
export async function drawFigure({ file, name, code }) {
  const src = code
    .replaceAll(/from 'driftlet\/plot'/g, `from '${recorder}'`)
    .replaceAll(/from 'driftlet(?:\/(\w+))?'/g, (_, sub) => `from '${entry(sub ?? 'index')}'`)
    .replaceAll('console.log(', '(() => {})(');
  globalThis.__figures = [];
  await import(moduleUrl(`${src}\n// run ${++runs}`)); // (a module's text is its identity: without this a second run would draw nothing)
  const drawn = globalThis.__figures;
  delete globalThis.__figures;
  if (drawn.length !== 1) throw new Error(`${file}: figure=${name} must draw one diagram (one bandDiagram() or levelChart() call); it drew ${drawn.length}`);
  return `${standalone(drawn[0], name)}\n`;
}

if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const dir = join(root, 'docs/figures');
  mkdirSync(dir, { recursive: true });
  const names = new Set();
  for (const block of figureBlocks()) {
    const path = join(dir, `${block.name}.svg`), svg = await drawFigure(block);
    const before = existsSync(path) ? readFileSync(path, 'utf8') : null;
    writeFileSync(path, svg);
    names.add(`${block.name}.svg`);
    console.log(`${before === null ? 'new      ' : before === svg ? 'unchanged' : 'CHANGED  '} docs/figures/${block.name}.svg  (${block.file}, ${(svg.length / 1024).toFixed(1)} kB)`);
  }
  for (const f of readdirSync(dir)) if (!names.has(f)) console.log(`no block draws docs/figures/${f}: delete it, or name it in a block's fence`);
}
