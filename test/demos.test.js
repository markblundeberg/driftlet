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
  const { electrode, theory, opened, planeRate, setSite, ions, reactant, L } = await pageTheory('frumkin.html', ['electrode', 'theory', 'opened', 'planeRate', 'setSite', 'ions', 'reactant', 'L']);
  const jlim = (2 * FARADAY * ions['S2O82-'].D * reactant) / L;
  const Es = Array.from({ length: 33 }, (_, k) => 0.4 - 0.05 * k);
  const runs = {};
  for (const a of [0.1, 0.7]) {
    setSite(a * 1e-9);
    for (const salt of [1, 10, 100, 316]) {
      const dev = new Device(electrode(salt, 0.4));
      let dip = 0, plateau = 0, kmax = 0, trim = 1, early = 0;
      const j = Es.map((E) => {
        dev.set({ contacts: { left: { V: E } } });
        const sol = dev.solve(), i = -sol.terminals.left.current, th = theory(salt, E), cl = theory(salt, E, 'classic');
        assert.ok(sol.converged && sol.warnings.length === 0, `${a} nm, ${salt} mM, ${E} V: ${sol.warnings.join(' ')}`);
        if (i < 0.5 * jlim) dip = Math.max(dip, Math.abs(i / th - 1));
        plateau = Math.max(plateau, Math.abs(i - th) / jlim);
        kmax = Math.max(kmax, planeRate(sol));
        trim = Math.min(trim, th / cl);
        if (planeRate(sol) < 0.1) early = Math.max(early, Math.abs(th / cl - 1));
        return i;
      });
      runs[`${a} ${salt}`] = j;
      // "about 5 m/s at −1.2 V in dilute salt, 10 to 20 in the strongest"
      if (salt === 1) assert.ok(planeRate(dev.solve()) > 4 && planeRate(dev.solve()) < 7, `${a} nm, 1 mM: k_f at −1.2 V ${planeRate(dev.solve())}`);
      if (salt === 316) assert.ok(kmax > 10, `${a} nm, 316 mM: k_f up to ${kmax}`);
      // "to 0.1% where the current is down", "the plateau by about 1%", "below 25 m/s in this window"
      assert.ok(dip < 1e-3 && plateau < 1.2e-2 && kmax < 25, `${a} nm, ${salt} mM: ${dip} in the dip, ${plateau} of the limit, k_f up to ${kmax} m/s`);
      // "costs little until ... k_f past about 1 m/s" (under 3% while it's below 0.1 m/s), "by up to a third at −1.2 V"
      assert.ok(early < 3e-2 && trim > 0.66, `${a} nm, ${salt} mM: the layer's resistance ${early} while k_f < 0.1 m/s, down to ${trim} of the classic`);
    }
  }
  // "with small sites the current turns and rises within the window, while on sites 0.7 nm across ... it keeps falling"
  const past = (j) => j.filter((_, k) => Es[k] < -0.3);
  const small = past(runs['0.1 1']), big = past(runs['0.7 1']);
  assert.ok(small.at(-1) > 1.5 * Math.min(...small), `0.1 nm: ${Math.min(...small)} at the least, ${small.at(-1)} at −1.2 V`);
  assert.ok(big.every((v, k) => k === 0 || v < big[k - 1]), `0.7 nm: ${big}`);
  // "the minimum moves the other way ... −0.79 V in 1 mM and −0.61 V in 0.1 M, where it's a dip of under 10%";
  // "it about doubles from the minimum to −1.2 V in dilute salt (the classic correction's triples)"
  setSite(0.1e-9);
  const fine = Array.from({ length: 101 }, (_, k) => -0.2 - 0.01 * k);
  const least = (salt, kind) => {
    const j = fine.map((E) => theory(salt, E, kind)), k = j.indexOf(Math.min(...j));
    return { E: fine[k], rise: j.at(-1) / j[k], depth: j[k] / theory(salt, 0.4, kind) };
  };
  const [m1, m100] = [least(1), least(100)], c1 = least(1, 'classic');
  assert.ok(Math.abs(m1.E + 0.79) < 0.015 && Math.abs(m100.E + 0.61) < 0.015 && m100.depth > 0.9, `${JSON.stringify(m1)}, ${JSON.stringify(m100)}`);
  assert.ok(Math.abs(m1.rise - 2.1) < 0.2 && Math.abs(c1.rise - 3) < 0.25, `${m1.rise}, ${c1.rise}`);
  // "However fast the reaction, that caps the current at the limiting current times L/(L + K) ... and the cap falls"
  for (const a of [0.1, 0.7]) {
    setSite(a * 1e-9);
    const cap = new Device(electrode(1, 0.4));
    let lastK = 0;
    for (const E of Es) {
      cap.set({ contacts: { left: { V: E } } });
      const sol = cap.solve();
      if (E > -0.3) continue;
      const { K } = opened(sol), j = -sol.terminals.left.current;
      assert.ok(K > lastK && j <= (1.02 * jlim * L) / (L + K), `${a} nm, ${E} V: K ${K} after ${lastK}, ${j / jlim} of the limit against ${L / (L + K)}`);
      lastK = K;
    }
  }
  // "Most of K builds up within an ångström or two of the plane"; and the plane's level is what K says it is.
  setSite(0.7e-9);
  const dev = new Device(electrode(1, 0.4));
  for (const E of Es) {
    dev.set({ contacts: { left: { V: E } } });
    const sol = dev.solve();
    if (E > -0.95) continue;
    const lay = opened(sol), u0 = lay.u[0], k2 = lay.x.findIndex((x) => x > 2e-10);
    const far = u0 + (lay.N * lay.K) / (ions['S2O82-'].D * lay.c['S2O82-']);
    assert.ok((lay.u[k2] - u0) / (lay.u.at(-1) - u0) > 0.5, `${E} V: ${(lay.u[k2] - u0) / (lay.u.at(-1) - u0)} of K within 2 Å`);
    assert.ok(Math.abs(far - 1) < 1e-6, `${E} V: the plane's level ${u0}, and back out ${far}`);
  }
});

test("demos/steps.html: Cottrell's and Sand's theory is the simulation's, as the page says", async () => {
  const { Device } = await import('../src/index.js');
  const { electrode, stepE, stepI, rates, E0, Do, Dr, c } = await pageTheory('steps.html', ['electrode', 'stepE', 'stepI', 'rates', 'E0', 'Do', 'Dr', 'c']);
  const rest = new Device(electrode({ I: 0 })).solve();
  assert.ok(Math.abs(rest.terminals.left.V - E0) < 1e-9, `at rest ${rest.terminals.left.V} V, E° ${E0}`);
  // "after about 1/H² ... (2 µs at −0.3 V, 5 ms at −0.1 V)"
  const settle = (eta) => {
    const { kf, kb } = rates(eta), H = kf / Math.sqrt(Do) + kb / Math.sqrt(Dr);
    return 1 / (H * H);
  };
  assert.ok(Math.abs(settle(-0.3) / 2e-6 - 1) < 0.05 && Math.abs(settle(-0.1) / 5e-3 - 1) < 0.06, `${settle(-0.3)}, ${settle(-0.1)} s`);
  // "within 5% of the line only after about ten times that"
  for (const eta of [-0.3, -0.1]) {
    const { kf, kb } = rates(eta), H = kf / Math.sqrt(Do) + kb / Math.sqrt(Dr), th = stepE(eta), t = 10 / (H * H);
    assert.ok(Math.abs(th.current(t) / th.cottrell(t) - 0.95) < 0.01, `${eta} V: ${th.current(t) / th.cottrell(t)} at 10/H²`);
  }
  // with a double layer, at η = −0.1 V: "1.8 times at 0.1 ms and 1.1 at 10 ms"
  {
    const dl = new Device(electrode(E0, true));
    dl.solve();
    dl.set({ contacts: { left: { V: E0 - 0.1 } } });
    const th = stepE(-0.1), r1 = -dl.advance(1e-4, { tol: 1e-4 }).terminals.left.current / th.current(1e-4), r2 = -dl.advance(1e-2, { tol: 1e-4 }).terminals.left.current / th.current(1e-2);
    assert.ok(Math.abs(r1 - 1.8) < 0.08 && Math.abs(r2 - 1.12) < 0.04, `${r1}, ${r2}`);
  }
  // A potential step: the current, with kinetics, from 10 µs to 1 s ("up to 1.2%" at first, the solution's
  // resistance; "0.2% from 1 ms on")
  for (const eta of [-0.5, -0.35, -0.2, -0.1, -0.05, -0.02]) {
    const dev = new Device(electrode(E0));
    dev.solve();
    dev.set({ contacts: { left: { V: E0 + eta } } });
    const th = stepE(eta);
    let worst = 0, profile = 0;
    for (const t of [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1]) {
      const sol = dev.advance(t, { tol: 1e-4 });
      const bound = t < 1e-4 ? 1.3e-2 : t < 1e-3 ? 6e-3 : 2e-3;
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
  // A current step: the potential until τ, and the transition at τ; with a double layer, read late
  // ("crossing E° − 0.25 V at 1.10τ at 0.2 mA/cm², 1.19τ at 0.5 mA/cm², and not by 1.3τ at 2 mA/cm²")
  for (const [i, late] of [[2, 1.1], [5, 1.19], [20, Infinity]]) {
    const th = stepI(i);
    const dev = new Device(electrode({ I: 0 }));
    dev.solve();
    dev.set({ contacts: { left: { I: -i } } });
    let worst = 0;
    for (const f of [0.1, 0.3, 0.5, 0.7, 0.9]) {
      const sol = dev.advance(f * th.tau, { tol: 1e-4 });
      worst = Math.max(worst, Math.abs(sol.terminals.left.V - th.potential(f * th.tau)));
    }
    const end = dev.advance(1.2 * th.tau, { tol: 1e-4 });
    assert.ok(worst < 1e-3, `${i} A/m²: E within ${worst} V`);
    assert.ok(!end.converged && Math.abs(end.time / th.tau - 1) < 5e-3, `${i} A/m²: the transition at ${end.time} s, τ ${th.tau}`);
    const dl = new Device(electrode({ I: 0 }, true));
    dl.solve();
    dl.set({ contacts: { left: { I: -i } } });
    let cross = Infinity;
    for (let k = 1; k <= 520 && cross === Infinity; k++) {
      const sol = dl.advance((k * 1.3 * th.tau) / 520, { tol: 1e-4 });
      assert.ok(sol.converged, `${i} A/m² with a double layer, at ${sol.time} s`);
      if (sol.terminals.left.V < E0 - 0.25) cross = sol.time / th.tau;
    }
    assert.ok(late === Infinity ? cross === Infinity : Math.abs(cross - late) < 0.01, `${i} A/m² with a double layer: past E° − 0.25 V at ${cross}τ`);
  }
});

test('demos/rde.html: the Koutecký–Levich and Tafel analysis gives back k⁰ and α, as the page says', async () => {
  const { sweep, theory, kinetic, rpms, Es, E0, Do, Dr, nu } = await pageTheory('rde.html', ['sweep', 'theory', 'kinetic', 'rpms', 'Es', 'E0', 'Do', 'Dr', 'nu']);
  for (const k0 of [1e-4, 1e-2, 0.3]) {
    const run = sweep(k0);
    let worst = 0;
    rpms.forEach((rpm, r) => Es.forEach((E, k) => {
      const th = theory(E, rpm, k0);
      if (Math.abs(th) > 1e-2) worst = Math.max(worst, Math.abs(run.I[r][k] / th - 1));
    }));
    // "less migration, which the supporting salt keeps to 0.1%"; "k⁰ and α come back within 0.05%"
    assert.ok(worst < 1e-3, `k⁰ ${k0}: the currents within ${worst} of the theory`);
    assert.ok(run.used.length >= 3 && Math.abs(run.back.k0 / k0 - 1) < 5e-4 && Math.abs(run.back.alpha / 0.5 - 1) < 5e-4, `k⁰ ${k0}: back ${JSON.stringify(run.back)} from ${run.used.length}`);
  }
  // "too fast for the disk, by the ×10 rule"; "the simulation's intercepts stay within 0.06% of the true kinetic
  // current however far past it they go"; with 0.5% noise, "about 4% at ten times the plateau and 30% at a hundred",
  // "k⁰ comes back within about 1% for a slow electrode, 3% at 0.01 cm/s and 12% at 0.03 cm/s"
  assert.equal(sweep(1).used.length, 0);
  const errs = (run, k0) => run.kl.filter((p) => p.E < E0 - 0.02).map((p) => [kinetic(p.E, k0) / run.plateau.at(-1), Math.abs(1 / p.a / kinetic(p.E, k0) - 1)]);
  const near = (rows, q0) => rows.reduce((best, r) => (Math.abs(Math.log(r[0] / q0)) < Math.abs(Math.log(best[0] / q0)) ? r : best))[1];
  const clean = errs(sweep(1), 1);
  assert.ok(Math.max(...clean.map((r) => r[0])) > 1000 && clean.every(([, e]) => e < 7e-4), JSON.stringify(clean));
  // over many draws of the noise: "an intercept scatters by about 10% at ten times the plateau, and by its own
  // size or more at a hundred"; "k⁰ comes back within about 1% for a slow electrode, 4% at 0.01 cm/s and 12% at
  // 0.03 cm/s (medians over draws), where a third of the draws leave fewer than the three potentials"
  const draws = (k0, n) => Array.from({ length: n }, (_, k) => sweep(k0, 0.005, k + 1));
  const median = (a) => [...a].sort((x, y) => x - y)[Math.floor(a.length / 2)];
  const rms = (a) => Math.sqrt(a.reduce((t, x) => t + x * x, 0) / a.length);
  const at01 = draws(0.1, 16), scatter = (lo, hi) => rms(at01.flatMap((run) => errs(run, 0.1).filter(([q]) => q > lo && q < hi).map(([, e]) => e)));
  assert.ok(scatter(8, 12) > 0.05 && scatter(8, 12) < 0.2 && scatter(80, 120) > 1, `${scatter(8, 12)}, ${scatter(80, 120)}`);
  for (const [k0, within] of [[1e-3, 0.012], [1e-2, 0.025], [0.1, 0.06]]) {
    const back = median((k0 === 0.1 ? at01 : draws(k0, 16)).map((run) => Math.abs(run.back.k0 / k0 - 1)));
    assert.ok(back < within && back > within / 4, `k⁰ ${k0} with noise: median ${back}`);
  }
  const fast = draws(0.316, 18), fits = fast.filter((run) => run.back);
  assert.ok(fits.length < 16 && fits.length > 6 && Math.abs(median(fits.map((run) => Math.abs(run.back.k0 / 0.316 - 1))) - 0.12) < 0.08, `${fits.length} fits`);
  // "about 3% above a real disk at this Schmidt number (Newman's correction)"; "1.19 times ... a real disk's would be 1.12"
  const sc = nu / Do, newman = 1 + 0.298 * Math.cbrt(1 / sc) + 0.14514 * Math.cbrt(1 / (sc * sc));
  assert.ok(Math.abs(newman - 1.03) < 0.005 && Math.abs(Dr / Do - 1.19) < 0.005 && Math.abs(Math.cbrt((Dr / Do) * (Dr / Do)) - 1.12) < 0.005, `${newman}, ${Dr / Do}`);
});

test('demos/randles.html: the Randles circuit, from the physics, is the simulated spectrum, as the page says', async () => {
  const { atRest, randles, fs } = await pageTheory('randles.html', ['atRest', 'randles', 'fs']);
  // "under 1% across these controls" (the corners of the sliders, and the middle)
  for (const p of [{ k0: 1, salt: 10, L: 3.2e-6 }, { k0: 1, salt: 10, L: 316e-6 }, { k0: 1e-4, salt: 1000, L: 316e-6 }, { k0: 1e-4, salt: 10, L: 3.2e-6 }, { k0: 1e-2, salt: 100, L: 17.8e-6 }]) {
    const dev = atRest(p), th = randles(p), { Z } = dev.impedance(fs);
    const worst = Math.max(...fs.map((f, k) => {
      const z = th.Z(f);
      return Math.hypot(Z.re[k] - z.re, Z.im[k] - z.im) / Math.hypot(z.re, z.im);
    }));
    assert.ok(worst < 1e-2, `${JSON.stringify(p)}: within ${worst}`);
  }
  // At the start: "the dip between the arcs sits about 11% beyond R_s + R_ct"; "the phase plot's tall peak ...
  // at 1/(2πC_dl√(R_s(R_s + R_ct))), near 1 kHz"; "the corner ... 37 Hz"
  const th = randles({ k0: 0.01, salt: 100, L: 10 ** 1.75 * 1e-6 }), ft = Array.from({ length: 4001 }, (_, k) => 10 ** (-2 + k / 500)), z = ft.map(th.Z);
  const fRC = 1 / (2 * Math.PI * th.Rct * th.Cdl);
  let dip = null;
  for (let k = 1; k < z.length - 1; k++) if (ft[k] < fRC && -z[k].im < -z[k - 1].im && -z[k].im < -z[k + 1].im) dip = z[k];
  const ph = z.map((w) => -Math.atan2(w.im, w.re)), peak = ft[ph.indexOf(Math.max(...ph.slice(1000)))];
  const formula = 1 / (2 * Math.PI * th.Cdl * Math.sqrt(th.Rs * (th.Rs + th.Rct)));
  assert.ok(Math.abs(dip.re / (th.Rs + th.Rct) - 1.11) < 0.02 && Math.abs(peak / formula - 1) < 0.02 && peak > 800 && peak < 1500 && Math.abs(fRC - 37) < 1, `${dip.re / (th.Rs + th.Rct)}, ${peak} against ${formula}, ${fRC}`);
});

test('demos/corrosion.html: the mixed potential and the Evans diagram against their closed forms, as the page says', async () => {
  const { Device, FARADAY, GAS_CONSTANT } = await import('../src/index.js');
  const { IONS } = await import('../src/kit.js');
  const { electrode, net, currents, corrosion } = await pageTheory('corrosion.html', ['electrode', 'net', 'currents', 'corrosion']);
  const B = (GAS_CONSTANT * 298.15) / FARADAY; // the anodic Tafel slope over 2.3, with α = ½ for two electrons
  for (const L of [10 ** -4.25, 10 ** -3.25, 1e-3]) {
    // "they agree to 0.1 mV and five digits"
    const rest = new Device(electrode({ I: 0 }, L)).solve(), th = corrosion(L), E = rest.terminals.left.V;
    const icorr = 2 * FARADAY * -rest.interfaces[0].rates[0], lim = currents(E, L).limit;
    assert.ok(rest.converged && Math.abs(E - th.E) < 1e-4 && Math.abs(icorr / th.i - 1) < 1e-5, `${L} m: E_corr ${E} against ${th.E}, i_corr ${icorr} against ${th.i}`);
    const dev = new Device(electrode(E, L));
    dev.solve();
    const at = (d) => (dev.set({ contacts: { left: { V: E + d } } }), dev.solve());
    if (L === 10 ** -3.25) {
      // the default: "in the range measured for steel in sea water (−0.35 to −0.45 V vs SHE), here at its
      // negative end"; "on that plateau"; "0.1 V below ... about 1% of its free rate, for a protection current
      // of about oxygen's limit"; Tafel extrapolation and Stern–Geary "within a few percent"
      const prot = currents(E - 0.1, L), tafel = 2 * FARADAY * 5e-6 * Math.exp((E - IONS['Fe2+'].mu0 / (2 * FARADAY)) / B); // the anodic Tafel line alone, at E_corr
      const Rp = 2e-3 / (at(1e-3).terminals.left.current - at(-1e-3).terminals.left.current);
      assert.ok(E < -0.42 && E > -0.45 && icorr / lim > 0.9, `${E} V, ${icorr / lim} of the limit`);
      assert.ok(prot.Fe / icorr < 0.02 && -prot.O2 / lim > 0.98, `${prot.Fe / icorr}, ${-prot.O2 / lim}`);
      assert.ok(Math.abs((Rp * icorr) / B - 1) < 0.035 && Math.abs(tafel / icorr - 1) < 0.05, `R_p·i_corr ${Rp * icorr} against ${B}; Tafel ${tafel / icorr}`);
      // fitted to the measured (net) points: "anywhere from half to 1.6 times i_corr depending on where you fit";
      // "add the plateau back and fit close in, 25 to 75 mV above E_corr, and it comes within 5%"
      const pts = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15].map((d) => [d, at(d).terminals.left.current]);
      const tafelFit = (lo, hi, back = 0) => {
        const p = pts.filter(([d]) => d > lo - 1e-9 && d < hi + 1e-9), n = p.length;
        const x = p.map((q) => q[0]), y = p.map((q) => Math.log(q[1] + back)), mx = x.reduce((a, v) => a + v, 0) / n, my = y.reduce((a, v) => a + v, 0) / n;
        const b = x.reduce((a, v, k) => a + (v - mx) * (y[k] - my), 0) / x.reduce((a, v) => a + (v - mx) * (v - mx), 0);
        return Math.exp(my - b * mx) / icorr;
      };
      const near = tafelFit(0.025, 0.075), far = tafelFit(0.1, 0.15), fixed = tafelFit(0.025, 0.075, lim);
      assert.ok(near > 0.45 && near < 0.65 && far > 1.45 && far < 1.7 && Math.abs(fixed - 1) < 0.06, `${near}, ${far}, ${fixed}`);
      // "about 0.15 V below, the Fermi level passes the iron's own redox level and iron would plate"
      assert.ok(currents(E - 0.13, L).Fe > 0 && currents(E - 0.17, L).Fe < 0, `${currents(E - 0.13, L).Fe}, ${currents(E - 0.17, L).Fe}`);
      // "9% by 0.1 V up, all of it the solution's potential at the face, 2.4 mV there: about half ohmic drop";
      // "with that taken out the two agree to 10⁻⁴"
      const up = at(0.1), g = up.region.indexOf(1), phi = up.phi[g], i = up.terminals.left.current;
      const kappa = ((FARADAY * FARADAY) / (GAS_CONSTANT * 298.15)) * (IONS['Na+'].D * 500 + IONS['Cl-'].D * 500), ohmic = (i * L) / kappa;
      const c = currents(E + 0.1 - phi, L);
      assert.ok(Math.abs(i / net(E + 0.1, L) - 0.91) < 0.01 && Math.abs(phi - 2.4e-3) < 1e-4 && ohmic / phi > 0.4 && ohmic / phi < 0.6, `${i / net(E + 0.1, L)}, φ ${phi}, ohmic ${ohmic}`);
      assert.ok(Math.abs(at(0.15).terminals.left.current / net(E + 0.15, L) - 0.68) < 0.01, '"by 0.15 V up it\'s 32% below"');
      assert.ok(Math.abs(i / (c.Fe + c.O2) - 1) < 2e-4, `with φ taken out: ${i / (c.Fe + c.O2)}`);
    }
    // "for the thinnest layers ... the kinetics share the control"
    if (L === 10 ** -4.25) assert.ok(icorr / lim < 0.6, `${icorr / lim}`);
    // "the cathodic side agrees everywhere"
    for (const d of [-0.25, -0.1, -0.05]) {
      const sol = at(d);
      assert.ok(sol.converged && sol.warnings.length === 0 && Math.abs(sol.terminals.left.current / net(E + d, L) - 1) < 1e-3, `${L} m, ${d} V`);
    }
  }
});

test('demos/lithium.html: Marcus–Hush–Chidsey plating against its closed form, and the numbers the page gives', async () => {
  const { FARADAY, GAS_CONSTANT } = await import('../src/index.js');
  const { mhcIntegral } = await import('../src/marcus.js');
  const { kinetic, steady, sweep, etas, limit, ceiling, states, mhcI, c, DLi } = await pageTheory('lithium.html', ['kinetic', 'steady', 'sweep', 'etas', 'limit', 'ceiling', 'states', 'mhcI', 'c', 'DLi']);
  const RT = GAS_CONSTANT * 298.15, VT = RT / FARADAY, kT = (lam) => (lam * FARADAY) / RT;
  // its own integral (Simpson) is the library's
  for (const [lam, a] of [[8.56, 0], [8.56, 9.7], [8.56, -9.7], [19.5, 30], [39, -20]]) assert.ok(Math.abs(mhcI(lam, a) / mhcIntegral(lam, a).I - 1) < 1e-12, `I(${lam}, ${a})`);
  // the simulated currents follow the full line: 3e-5 at the defaults, 3e-4 across the controls
  for (const [lam, r, tol] of [[0.22, 5e-6, 3e-5], [0.22, 1e-6, 3e-4], [0.22, 25e-6, 3e-4], [0.5, 5e-6, 3e-4], [0.5, 25e-6, 3e-4], [1, 1e-6, 3e-4]]) {
    const I = sweep(lam, r);
    I.forEach((i, k) => {
      const th = steady(etas[k], lam, r);
      assert.ok(Math.abs(i - th) < tol * Math.max(Math.abs(th), 1e-2), `λ ${lam} eV, r ${r} m, η ${etas[k]}: ${i} vs ${th}`);
    });
  }
  const mA = (j) => j / 10, r = 5e-6;
  // the kinetic ceiling, 295 mA/cm²; a third of a 5 µm sphere's limit; Butler–Volmer three times the current at −0.25 V; the plateau 225
  assert.equal(Math.round(mA(ceiling(0.22))), 295);
  assert.ok(Math.abs(ceiling(0.22) / limit(r) - 1 / 3) < 0.03, `${ceiling(0.22) / limit(r)}`);
  const bv = steady(-0.25, 0.22, r, 'bv') / steady(-0.25, 0.22, r);
  assert.ok(bv > 2.9 && bv < 3.1, `Butler–Volmer ×${bv}`);
  // past the data: Marcus–Hush–Chidsey levels off below its ceiling, Butler–Volmer runs on to the limit
  const plateau = steady(-0.6, 0.22, r);
  assert.ok(Math.round(mA(plateau)) === 225 && Math.abs(steady(-0.55, 0.22, r) / plateau - 1) < 0.01, `plateau ${plateau}`);
  assert.ok(steady(-0.6, 0.22, r, 'bv') / limit(r) > 0.98);
  // on 25 µm, transport caps the current below the ceiling; λ = 0.5 eV moves the plateau up
  assert.ok(limit(25e-6) < ceiling(0.22));
  assert.ok(steady(-0.6, 0.5, r) > 3 * plateau);
  // α = ½ at η = 0, and the local slope ½ − |η|F/4λ while |η| ≪ λ
  // (the forward rate's slope, d ln I/da at a = −η/kT, by the page's own integral)
  const slope = (eta, lam) => (Math.log(mhcI(kT(lam), -eta / VT + 1e-3)) - Math.log(mhcI(kT(lam), -eta / VT - 1e-3))) / 2e-3;
  assert.ok(Math.abs(slope(0, 0.22) - 0.5) < 1e-6, `${slope(0, 0.22)}`);
  assert.ok(Math.abs(kinetic(-1e-3, 0.22, c) / ((86 * 1e-3) / VT) - 1) < 1e-3, 'the exchange current, from the slope at η = 0');
  // (falling 30% slower than Marcus's ½ − F|η|/2λ at first)
  const fall = (slope(0.001, 0.22) - slope(-0.001, 0.22)) / 0.002 / (1 / (2 * 0.22)); // dα/dη in units of F/2λ
  assert.ok(Math.abs(fall - 0.7) < 0.01, `${fall}`);
  // Gerischer's overlaps are the forward and backward rates: their areas, against the integral
  for (const [eta, lam, cs] of [[-0.25, 0.22, 800], [0.1, 0.22, 1200], [-0.5, 0.6, 500]]) {
    const s = states(eta, lam, cs), dE = s.E[1] - s.E[0], area = (v) => v.reduce((t, x) => t + x, 0) * dE;
    const a0 = -eta / VT + Math.log(cs / c);
    const fwd = (cs / c) * mhcI(kT(lam), a0), back = mhcI(kT(lam), -a0);
    assert.ok(Math.abs(area(s.plating) / area(s.stripping) / (fwd / back) - 1) < 1e-3, `η ${eta}: ${area(s.plating) / area(s.stripping)} vs ${fwd / back}`);
    assert.ok(Math.abs((area(s.plating) * Math.sqrt(4 * Math.PI * kT(lam))) / fwd - 1) < 1e-3, `the area is the rate: ${(area(s.plating) * Math.sqrt(4 * Math.PI * kT(lam))) / fwd}`);
  }
  // Zeng et al.'s closed form at λ = 0.22 eV: 4% low at η = 0, 4.5% high at 0.25 V; the exact integral
  // puts the same curve at λ = 0.23 eV with an exchange current 3% higher, within about 1% rms
  const erfc = (x) => {
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const v = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? v : 2 - v;
  };
  const zeng = (l, a) => (Math.sqrt(Math.PI * l) / (1 + Math.exp(-a))) * erfc((l - Math.sqrt(1 + Math.sqrt(l) + a * a)) / (2 * Math.sqrt(l)));
  const L = kT(0.22), I = (l, a) => mhcIntegral(l, a).I;
  assert.ok(Math.abs(zeng(L, 0) / I(L, 0) - 1 + 0.04) < 0.005 && Math.abs(zeng(L, 0.25 / VT) / I(L, 0.25 / VT) - 1.045) < 0.005);
  const ev = Array.from({ length: 50 }, (_, k) => -0.25 + 0.01 * k + (k >= 25 ? 0.01 : 0));
  const net = (f, l, e) => (f(l, e / VT) - f(l, -e / VT)) / f(l, 0);
  let best = null;
  for (let lam = 0.2; lam <= 0.26; lam += 0.0005) {
    const d = ev.map((e) => Math.log(net(I, kT(lam), e) / net(zeng, L, e))), m = d.reduce((s, x) => s + x, 0) / d.length;
    const rms = Math.sqrt(d.reduce((s, x) => s + (x - m) * (x - m), 0) / d.length);
    if (!best || rms < best.rms) best = { lam, rms, i0: Math.exp(-m) };
  }
  assert.ok(Math.abs(best.lam - 0.23) < 0.005 && Math.abs(best.i0 - 1.03) < 0.01 && best.rms < 0.012, JSON.stringify(best));
  // the ideal conductivity, Nernst–Einstein's: 2.4 S/m
  const kappa = (FARADAY * FARADAY * c * (DLi + (DLi * 0.62) / 0.38)) / RT;
  assert.ok(Math.abs(kappa - 2.4) < 0.05, `${kappa}`);
});
