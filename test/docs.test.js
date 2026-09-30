import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join } from 'node:path';

// Every ```js block in the README and docs runs as a module against the source, so examples
// can't rot. Blocks marked ```js nocheck are skipped (sketches, fragments).
const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const entry = pathToFileURL(join(root, 'src/index.js')).href;
const docs = existsSync(join(root, 'docs')) ? readdirSync(join(root, 'docs')) : [];
const files = ['README.md', ...docs.filter((f) => f.endsWith('.md')).map((f) => `docs/${f}`)];

for (const file of files) {
  const text = readFileSync(join(root, file), 'utf8');
  const blocks = [...text.matchAll(/```js\n([\s\S]*?)```/g)].map((m) => m[1]);
  blocks.forEach((code, k) => {
    test(`${file}: example ${k + 1} runs`, async () => {
      const logs = [];
      const src = code
        .replaceAll("from 'driftlet'", `from '${entry}'`)
        .replaceAll('console.log(', '__log(');
      globalThis.__log = (...a) => logs.push(a);
      await import(`data:text/javascript;base64,${Buffer.from(src).toString('base64')}`);
      for (const line of logs) for (const v of line) if (typeof v === 'number') assert.ok(Number.isFinite(v), `${file}: ${line}`);
    });
  });
}
