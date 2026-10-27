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

// Pages that draw their simulation against their own theory say how well they agree: here each one's
// script, up to where it needs a page, is run against the claims its text makes.
async function pageTheory(page, names) {
  const html = readFileSync(join(demos, page), 'utf8');
  const script = /<script type="module">([\s\S]*?)<\/script>/.exec(html)[1];
  const at = (path) => pathToFileURL(join(root, path)).href;
  const head = script
    .slice(0, script.indexOf('// (All above runs without a page'))
    .replace(/^import .*'\.\/(plot|common)\.js';$/gm, '')
    .replaceAll("'../src/", `'${at('src/')}`);
  assert.ok(head.length > 1000 && !/document\./.test(head), 'the page-free part of its script');
  return import(`data:text/javascript;base64,${Buffer.from(`${head}\nexport { ${names.join(', ')} };`).toString('base64')}`);
}

test("demos/frumkin.html: its theory is the simulation's, as the page says", async () => {
  const { Device, FARADAY } = await import('../src/index.js');
  const { electrode, theory, resistance, opened, ions, reactant, L, VT } = await pageTheory('frumkin.html', ['electrode', 'theory', 'resistance', 'opened', 'ions', 'reactant', 'L', 'VT']);
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

test("demos/steps.html: Cottrell's and Sand's theory is the simulation's, as the page says", async () => {
  const { Device } = await import('../src/index.js');
  const { electrode, stepE, stepI, E0, Do, Dr, c } = await pageTheory('steps.html', ['electrode', 'stepE', 'stepI', 'E0', 'Do', 'Dr', 'c']);
  const rest = new Device(electrode({ I: 0 })).solve();
  assert.ok(Math.abs(rest.terminals.left.V - E0) < 1e-9, `at rest ${rest.terminals.left.V} V, E° ${E0}`);
  // A potential step: the current, with kinetics, from 10 µs to 1 s ("up to 1.5%" at first, the solution's
  // resistance; "0.15% from 1 ms on")
  for (const eta of [-0.5, -0.35, -0.2, -0.1, -0.05]) {
    const dev = new Device(electrode(E0));
    dev.solve();
    dev.set({ contacts: { left: { V: E0 + eta } } });
    const th = stepE(eta);
    let worst = 0, profile = 0;
    for (const t of [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1]) {
      const sol = dev.advance(t, { tol: 1e-4 });
      const bound = t < 1e-4 ? 1.6e-2 : t < 1e-3 ? 6e-3 : 1.6e-3;
      worst = Math.max(worst, Math.abs(-sol.terminals.left.current / th.current(t) - 1) / bound);
      // (the profiles until 0.1 s: by 1 s the semi-infinite theory is 0.4% from the bath's at its edge, 100 µm out)
      const x0 = sol.x[sol.region.indexOf(1)];
      sol.x.forEach((x, g) => {
        if (sol.region[g] !== 1 || t > 0.1) return;
        const off = Math.max(Math.abs(sol.c['Fe3+'][g] - th.profile(x - x0, t, Do, -1)), Math.abs(sol.c['Fe2+'][g] - th.profile(x - x0, t, Dr, 1)));
        profile = Math.max(profile, off / (2 * bound * c)); // (the same shift as the current's, near the surface)
      });
    }
    assert.ok(worst < 1 && profile < 1, `η ${eta}: the current at ${worst} of its bound, the profiles at ${profile} of theirs`);
  }
  // A current step: the potential until τ, and the transition at τ
  for (const i of [2, 5, 20]) {
    const dev = new Device(electrode({ I: 0 }));
    dev.solve();
    dev.set({ contacts: { left: { I: -i } } });
    const th = stepI(i);
    let worst = 0;
    for (const f of [0.1, 0.3, 0.5, 0.7, 0.9]) {
      const sol = dev.advance(f * th.tau, { tol: 1e-4 });
      worst = Math.max(worst, Math.abs(sol.terminals.left.V - th.potential(f * th.tau)));
    }
    const end = dev.advance(1.2 * th.tau, { tol: 1e-4 });
    assert.ok(worst < 1e-3, `${i} A/m²: E within ${worst} V`);
    assert.ok(!end.converged && Math.abs(end.time / th.tau - 1) < 5e-3, `${i} A/m²: the transition at ${end.time} s, τ ${th.tau}`);
  }
});

test('demos/rde.html: the Koutecký–Levich and Tafel analysis gives back k⁰ and α, as the page says', async () => {
  const { sweep, theory, rpms, Es } = await pageTheory('rde.html', ['sweep', 'theory', 'rpms', 'Es']);
  for (const k0 of [1e-4, 1e-2, 0.3]) {
    const run = sweep(k0);
    let worst = 0;
    rpms.forEach((rpm, r) => Es.forEach((E, k) => {
      const th = theory(E, rpm, k0);
      if (Math.abs(th) > 1e-2) worst = Math.max(worst, Math.abs(run.I[r][k] / th - 1));
    }));
    // "less migration, which the 1 M KCl keeps to 0.1%"; "k⁰ comes back within 0.02%"
    assert.ok(worst < 1.2e-3, `k⁰ ${k0}: the currents within ${worst} of the theory`);
    assert.ok(run.used.length >= 3 && Math.abs(run.back.k0 / k0 - 1) < 3e-4 && Math.abs(run.back.alpha - 0.5) < 5e-4, `k⁰ ${k0}: back ${JSON.stringify(run.back)} from ${run.used.length}`);
  }
  // "an electrode too fast for a disk to measure"
  assert.equal(sweep(1).used.length, 0);
});
