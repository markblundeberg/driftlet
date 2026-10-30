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
  const { kinetic, steady, sweep, etas, limit, ceiling, states, mhcI, c, DLi, marcusHush } = await pageTheory('lithium.html', ['kinetic', 'steady', 'sweep', 'etas', 'limit', 'ceiling', 'states', 'mhcI', 'c', 'DLi', 'marcusHush']);
  const RT = GAS_CONSTANT * 298.15, VT = RT / FARADAY, kT = (lam) => (lam * FARADAY) / RT;
  // its own integral (Simpson) is the library's
  for (const [lam, a] of [[8.56, 0], [8.56, 9.7], [8.56, -9.7], [19.5, 30], [39, -20]]) assert.ok(Math.abs(mhcI(lam, a) / mhcIntegral(lam, a).I - 1) < 1e-12, `I(${lam}, ${a})`);
  // the simulated currents follow the full line: 3e-5 at the defaults, 3e-4 across the controls
  for (const [lam, r, tol] of [[0.22, 5e-6, 3e-5], [0.22, 1e-6, 5e-4], [0.22, 25e-6, 5e-4], [0.5, 5e-6, 5e-4], [0.5, 25e-6, 5e-4], [1, 1e-6, 5e-4], [1, 25e-6, 5e-4], [0.1, 25e-6, 5e-4]]) {
    const I = sweep(lam, r);
    I.forEach((i, k) => {
      const th = steady(etas[k], lam, r);
      assert.ok(Math.abs(i - th) < tol * Math.max(Math.abs(th), 1e-2), `λ ${lam} eV, r ${r} m, η ${etas[k]}: ${i} vs ${th}`);
    });
  }
  const mA = (j) => j / 10, r = 5e-6;
  // the two fits, Marcus–Hush (0.34 eV, 8.8) and MHC (0.22 eV, 8.6), within 8% over the measured ±0.25 V;
  // past it MHC levels off and Marcus–Hush turns down
  for (let eta = -0.25; eta <= 0.25001; eta += 0.01) {
    if (Math.abs(eta) < 1e-9) continue;
    assert.ok(Math.abs(marcusHush(eta) / kinetic(eta, 0.22, c) - 1) < 0.08, `η ${eta}: ${marcusHush(eta) / kinetic(eta, 0.22, c)}`);
  }
  assert.ok(marcusHush(-0.6) < marcusHush(-0.45) && kinetic(-0.6, 0.22, c) > kinetic(-0.45, 0.22, c));
  // the ceiling: roughly 300 mA/cm², a quarter more or less for ±0.02 eV in λ, nearer 350 with the exact integral's pair
  assert.ok(Math.abs(ceiling(0.24) / ceiling(0.22) - 1.25) < 0.03 && Math.abs(ceiling(0.22) / ceiling(0.2) - 1.25) < 0.03, `${ceiling(0.24) / ceiling(0.22)}`);
  assert.ok(Math.abs((1.034 * ceiling(0.2326)) / 10 / 350 - 1) < 0.02, `${(1.034 * ceiling(0.2326)) / 10}`);
  // by 0.25 V, Butler–Volmer six times the kinetic current
  assert.ok(Math.abs(kinetic(-0.25, 0.22, c, 'bv') / kinetic(-0.25, 0.22, c) / 6.3 - 1) < 0.02);
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
  // puts the same curve at λ = 0.23 (to 0.24) eV with an exchange current 3% higher (within 5%)
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
  assert.ok(Math.abs(best.lam - 0.233) < 0.001 && Math.abs(best.i0 - 1.034) < 0.005 && best.rms < 0.012, JSON.stringify(best));
  // weighed on the currents themselves instead: 0.226 eV and 6% higher, its ceiling 335 mA/cm² (the page: 0.226 to 0.233, 3 to 6%, 335 to 350)
  let abs = null;
  for (let lam = 0.2; lam <= 0.26; lam += 0.0005) {
    const e = ev.map((x) => net(I, kT(lam), x)), z = ev.map((x) => net(zeng, L, x));
    const sc = z.reduce((t, v, k) => t + v * e[k], 0) / e.reduce((t, v) => t + v * v, 0), res = e.reduce((t, v, k) => t + (sc * v - z[k]) * (sc * v - z[k]), 0);
    if (!abs || res < abs.res) abs = { lam, res, i0: sc };
  }
  assert.ok(Math.abs(abs.lam - 0.2265) < 0.001 && Math.abs(abs.i0 - 1.056) < 0.006, JSON.stringify(abs));
  assert.ok(Math.abs((abs.i0 * ceiling(abs.lam)) / 10 - 335) < 3 && Math.abs((best.i0 * ceiling(best.lam)) / 10 - 350) < 4, `${(abs.i0 * ceiling(abs.lam)) / 10}`);
  // the ideal conductivity, Nernst–Einstein's: 2.4 S/m
  const kappa = (FARADAY * FARADAY * c * (DLi + (DLi * 0.62) / 0.38)) / RT;
  assert.ok(Math.abs(kappa - 2.4) < 0.05, `${kappa}`);
});

test("demos/chidsey.html: Chidsey's experiment, the fitted rates against the theory, and the numbers the page gives", async () => {
  const { mhcIntegral } = await import('../src/marcus.js');
  const { step, rate, etas, mhcI, kT, A, states, marcusFit } = await pageTheory('chidsey.html', ['step', 'rate', 'etas', 'mhcI', 'kT', 'A', 'states', 'marcusFit']);
  // one-level Marcus with its own λ over −1.0 to +0.8 V: 1.02 eV, within 9%; and it turns down past it
  const one = marcusFit(0.85, 298.15);
  assert.ok(Math.abs(one.lam - 1.02) < 0.005 && one.dev < 0.09 && one.dev > 0.08, JSON.stringify(one));
  assert.ok(rate(1.4, 0.85, 298.15, 'marcus1', one) / rate(1.4, 0.85, 298.15) < 0.2);
  for (const [lam, a] of [[33.1, 0], [33.1, 30], [33.1, -30], [19.5, 15]]) assert.ok(Math.abs(mhcI(lam, a) / mhcIntegral(lam, a).I - 1) < 1e-12);
  // the rate constants fitted to the simulated decays, against the theory: 2e-4, across η, T and λ
  for (const [T, lam] of [[298.15, 0.85], [274.15, 0.85], [320.15, 0.5], [298.15, 1.2]]) {
    for (const e of [-0.9, -0.3, 0.1, 0.5, 0.8]) {
      const s = step(e, lam, T);
      assert.ok(Math.abs(s.k / rate(e, lam, T) - 1) < 2e-4, `T ${T}, λ ${lam}, η ${e}: ${s.k} vs ${rate(e, lam, T)}`);
    }
  }
  assert.equal(etas.length, 18);
  // Chidsey's νρ (6.73e4, 7e4 in his abstract): at the formal potential 1.30 s⁻¹ each way at 25 °C, so 2.60 for
  // the step (4% over the 2.5 he set it by), 0.58 at 1 °C, 2.46 at 47 °C; the ceiling 3.5e4 s⁻¹
  assert.equal(A, 6.73e4);
  const k0 = (T) => rate(1e-9, 0.85, T) / 2;
  assert.ok(Math.abs(k0(298.15) - 1.30) < 0.005 && Math.abs(k0(274.15) - 0.58) < 0.005 && Math.abs(k0(320.15) - 2.46) < 0.005);
  assert.ok(Math.abs(2 * k0(298.15) / 2.5 - 1.04) < 0.005);
  assert.ok(Math.abs(A * kT(298.15) * Math.sqrt((4 * Math.PI * 0.85) / kT(298.15)) / 3.5e4 - 1) < 0.01);
  // by 0.8 V: Butler–Volmer 500 times faster; one-level Marcus at its peak, η = λ; the temperature 10% there, ×4 at E°′
  assert.ok(Math.abs(rate(0.8, 0.85, 298.15, 'bv') / rate(0.8, 0.85, 298.15) / 523 - 1) < 0.01);
  const mk = (e) => rate(e, 0.85, 298.15, 'marcus');
  assert.ok(mk(0.85) > mk(0.8) && mk(0.85) > mk(0.9));
  assert.ok(Math.abs(rate(0.8, 0.85, 320.15) / rate(0.8, 0.85, 274.15) - 1.10) < 0.01 && Math.abs(k0(320.15) / k0(274.15) - 4.2) < 0.05);
  // the ceiling's share: half for the oxidation at η = λ, 41% for the sum at 0.8 V, 99.9% by 1.5 V
  const ceil = (T) => A * kT(T) * Math.sqrt((4 * Math.PI * 0.85) / kT(T)), L25 = 0.85 / kT(298.15);
  assert.ok(Math.abs((A * kT(298.15) * mhcI(L25, L25)) / ceil(298.15) - 0.5) < 1e-9);
  assert.ok(Math.abs(rate(0.8, 0.85, 298.15) / ceil(298.15) - 0.41) < 0.005 && rate(1.5, 0.85, 298.15) / ceil(298.15) > 0.998);
  // the activation energy at the formal potential, 0.24 eV (λ/4 = 0.21); 2.60 s⁻¹ for the step's relaxation there
  const Ea = (Math.log(k0(320.15) / k0(274.15)) * (kT(298.15) / 298.15)) / (1 / 274.15 - 1 / 320.15); // (R/F = kT/T)
  assert.ok(Math.abs(Ea - 0.2365) < 0.001, `${Ea}`);
  assert.ok(Math.abs(rate(1e-9, 0.85, 298.15) - 2.60) < 0.005);
  // with the coupling held instead of A (A ∝ 1/√T): 2% at 0.8 V, 3.9-fold at the formal potential
  const fixed = (e, T) => rate(e, 0.85, T) * Math.sqrt(298.15 / T);
  assert.ok(Math.abs(fixed(0.8, 320.15) / fixed(0.8, 274.15) - 1.018) < 0.002 && Math.abs(fixed(1e-9, 320.15) / fixed(1e-9, 274.15) - 3.90) < 0.01);
  // Gerischer's overlaps are the oxidation and reduction rates: their areas against the integral
  for (const [e, lam, T] of [[0.5, 0.85, 298.15], [-0.3, 0.6, 320.15]]) {
    const E = Array.from({ length: 4001 }, (_, k) => -3 + 0.0015 * k), s = states(e, lam, T, E), area = (v) => v.reduce((t, x) => t + x, 0) * 0.0015;
    const L = lam / kT(T), want = mhcI(L, e / kT(T)) / Math.sqrt(4 * Math.PI * L);
    assert.ok(Math.abs(area(s.oxidation) / want - 1) < 1e-3 && Math.abs(area(s.reduction) / (mhcI(L, -e / kT(T)) / Math.sqrt(4 * Math.PI * L)) - 1) < 1e-3);
  }
});

test('demos/hydrogen.html: the three routes against their steady theory, the Tafel slopes and the volcano, as the page says', async () => {
  const { FARADAY } = await import('../src/index.js');
  const { steady, sweep, volcano, exchange, Es, dGs, MECHS } = await pageTheory('hydrogen.html', ['steady', 'sweep', 'volcano', 'exchange', 'Es', 'dGs', 'MECHS']);
  // the simulated currents and coverages against the theory: 4e-6 (the most at 4.6 A/cm²)
  for (const mech of ['volmer', 'heyrovsky', 'tafel']) {
    for (const dG of [-0.2, 0.1, 0.3]) {
      sweep(dG, mech).forEach((s, j) => {
        const th = steady(Es[j], dG, mech);
        assert.ok(Math.abs(s.i - th.i) < 4e-6 * Math.max(Math.abs(th.i), 1e-6) && Math.abs(s.theta / th.theta - 1) < 4e-6, `${mech} ${dG} ${Es[j]}: ${s.i} vs ${th.i}`);
      });
    }
  }
  // the volcano's points: 1e-7
  for (const mech of ['volmer', 'heyrovsky', 'tafel']) {
    volcano(mech).forEach((v, j) => {
      const th = exchange(dGs[j], mech);
      assert.ok(Math.abs(v / th - 1) < 1e-7, `${mech} ${dGs[j]}: ${v} vs ${th}`);
    });
  }
  // the Tafel slopes (mV per decade, from ±2 mV): a slow Volmer step 118 on a weak binder and 40 on a strong
  // one (θ near 1), the mirror image of a slow Heyrovsky step's 40 on a weak binder, which bends up through
  // 118 as the surface fills; a slow Tafel step 29.7, then its limit on a full surface, 2F k_T e^{ΔG/kT}
  const slope = (E, dG, mech) => 4 / Math.log10(steady(E - 0.002, dG, mech).i / steady(E + 0.002, dG, mech).i);
  assert.ok(Math.abs(slope(-0.25, 0.2, 'volmer') - 118.3) < 0.2);
  assert.ok(Math.abs(slope(-0.1, -0.2, 'volmer') - 40) < 1 && steady(-0.1, -0.2, 'volmer').theta > 0.9);
  assert.ok(Math.abs(slope(-0.075, 0.2, 'heyrovsky') - 39.5) < 0.3 && steady(-0.075, 0.2, 'heyrovsky').theta < 0.01);
  assert.ok(Math.abs(slope(-0.29, 0.2, 'heyrovsky') - 118) < 6);
  assert.ok(Math.abs(slope(-0.4, 0.2, 'heyrovsky') - 139) < 1 && Math.abs(steady(-0.4, 0.2, 'heyrovsky').x - 0.9) < 0.01); // (the acid's supply)
  // the mirror: the Volmer route at −ΔG the Heyrovsky route at +ΔG, the same current, covered and bare swapped
  for (const E of [0.05, -0.1, -0.2, -0.3, -0.4]) {
    const a = steady(E, -0.2, 'volmer'), b = steady(E, 0.2, 'heyrovsky');
    assert.ok(Math.abs(a.i / b.i - 1) < 1e-9 && Math.abs(a.theta + b.theta - 1) < 1e-9, `mirror at ${E}`);
  }
  // the Tafel route's limit, i₀(1 + e^{ΔG/kT})²: 1.06 i₀ at platinum's −0.09 eV (i₀ from ±1 mV, so 1e-3)
  for (const dG of [-0.09, 0, 0.1]) assert.ok(Math.abs(steady(-1, dG, 'tafel').i / exchange(dG, 'tafel') / (1 + Math.exp((dG * FARADAY) / (8.314462618 * 298.15))) ** 2 - 1) < 1e-3, `tafel limit ${dG}`);
  assert.ok(Math.abs(slope(-0.075, 0.2, 'tafel') - 29.7) < 0.3);
  const lim = 2 * FARADAY * MECHS.tafel.T * Math.exp((-0.2 * FARADAY) / (8.314462618 * 298.15));
  assert.ok(Math.abs(steady(-0.6, -0.2, 'tafel').i / lim - 1) < 1e-3 && steady(-0.6, -0.2, 'tafel').theta > 0.99);
  // the volcano: exactly symmetric about 0 for every route, its legs a decade per 118 meV (per 59 for the Tafel route)
  for (const mech of ['volmer', 'heyrovsky', 'tafel']) for (const d of [0.05, 0.15, 0.25]) assert.ok(Math.abs(exchange(d, mech) / exchange(-d, mech) - 1) < 1e-6, `${mech} ±${d}`);
  const leg = (mech) => 100 / Math.log10(exchange(0.2, mech) / exchange(0.3, mech)); // meV per decade
  assert.ok(Math.abs(leg('heyrovsky') - 118.3) < 1.5 && Math.abs(leg('volmer') - 118.3) < 1.5 && Math.abs(leg('tafel') - 59.2) < 1, `${leg('heyrovsky')} ${leg('tafel')}`);
});

test('demos/mott-schottky.html: the capacitance against Poisson–Boltzmann, both readings against the cell, and the fits, as the page says', async () => {
  const { Device, units } = await import('../src/index.js');
  const { theory, measure, reading, fit, cell } = await pageTheory('mott-schottky.html', ['theory', 'measure', 'reading', 'fit', 'cell']);
  const fd = Array.from({ length: 21 }, (_, k) => 10 ** (2 + 0.25 * k));
  for (const [nd, c, holes] of [[1e15, 1, 'none'], [1e16, 100, 'none'], [1e18, 1000, 'none'], [1e16, 100, 'equilibrium']]) {
    const p = { N: units.perCm3(nd), c, holes }, th = theory(p);
    // the flat band on the SHE scale (−0.186 V at 10¹⁶ cm⁻³): where the simulated face holds no charge
    const D = (E) => new Device(cell(E, p)).solve().interfaces[0].D;
    assert.ok(D(th.Vfb - 1e-4) < 0 && D(th.Vfb + 1e-4) > 0, `${nd}, ${c}: flat band ${th.Vfb}`);
    // C read in series at 100 Hz: within 0.5% from 0.2 V below the flat band up (short of inversion, with holes)
    for (const q of measure(p, [100])) {
      if (q.E < th.Vfb - 0.2 || (holes === 'equilibrium' && q.E > th.Einv - 0.15)) continue;
      assert.ok(Math.abs(reading(q.Z[0], 100, 'series') / th.C(q.E) - 1) < 5e-3, `${nd}, ${c}, ${holes} at ${q.E}: ${reading(q.Z[0], 100, 'series')} vs ${th.C(q.E)}`);
    }
    // both readings of Z, 100 Hz to 10 MHz, against the cell's circuit (the solution's ε/L across its R), in depletion
    for (const E of [th.Vfb + 0.15, th.Vfb + 0.45]) {
      const [h] = measure(p, fd, [E]);
      h.Z.forEach((z, k) => {
        for (const how of ['series', 'parallel']) assert.ok(Math.abs(reading(z, fd[k], how) / reading(th.Z(fd[k], E), fd[k], how) - 1) < 3e-3, `${nd}, ${c}, ${E}, ${fd[k]} Hz, ${how}`);
      });
    }
  }
  assert.ok(Math.abs(theory({ N: units.perCm3(1e16), c: 100 }).Vfb + 0.1863) < 1e-4);
  // the fits (1 kHz, in series): from 0.1 V above the flat band about 1% high and 2–4 mV low; from 0.3 V
  // above (the default, +0.1 to +0.5 V at 10¹⁶ cm⁻³), within 0.2% and 0.4 mV, as the theory's own points fit
  const p = { N: units.perCm3(1e16), c: 100 }, th = theory(p);
  const sim = measure(p, [1e3]).map((q) => ({ E: q.E, y: 1 / reading(q.Z[0], 1e3, 'series') ** 2 }));
  const near = fit(sim, -0.1, 0.3), dflt = fit(sim, 0.1, 0.5), own = fit(sim.map((q) => ({ E: q.E, y: 1 / th.C(q.E) ** 2 })), 0.1, 0.5);
  assert.ok(near.N / p.N - 1 > 0.005 && near.N / p.N - 1 < 0.012 && near.Vfb - th.Vfb < -0.002 && near.Vfb - th.Vfb > -0.004, `near: ${near.N / p.N - 1}, ${near.Vfb - th.Vfb}`);
  assert.ok(Math.abs(dflt.N / p.N - 1) < 2e-3 && Math.abs(dflt.Vfb - th.Vfb) < 4e-4 && Math.abs(dflt.N / own.N - 1) < 5e-4, `default: ${dflt.N / p.N - 1}, ${dflt.Vfb - th.Vfb}`);
  // across accumulation, far off
  assert.ok(fit(sim, -0.45, 0.1).N / p.N - 1 > 0.5);
  // read in parallel at 1 MHz: the slope nearly right, the flat band 31 mV low (1/C_p² = 1/C² + 2(ωR)² + …)
  const par = fit(measure(p, [1e6]).map((q) => ({ E: q.E, y: 1 / reading(q.Z[0], 1e6, 'parallel') ** 2 })), 0.1, 0.5);
  assert.ok(Math.abs(par.N / p.N - 1) < 2e-3 && Math.abs(par.Vfb - th.Vfb + 0.031) < 1.5e-3, `parallel: ${par.N / p.N - 1}, ${par.Vfb - th.Vfb}`);
  // heavily doped in 1 mM: 8% of E − E_fb in the solution at +0.3 V (3.5% of a change), the doping 2–3% high
  // and the flat band a few mV low; read in series at 100 kHz, C down by about a third; read in parallel there, the line turns over
  const hp = { N: units.perCm3(1e18), c: 1 }, hth = theory(hp);
  assert.ok(Math.abs((hth.split(0.3).yw * 0.025693) / (0.3 - hth.Vfb) - 0.083) < 0.003 && Math.abs(hth.C(0.3) / hth.Cw(0.3) - 0.035) < 0.002);
  const hs = measure(hp, [1e3]).map((q) => ({ E: q.E, y: 1 / reading(q.Z[0], 1e3, 'series') ** 2 }));
  for (const [a, b] of [[0.1, 0.5], [hth.Vfb + 0.1, hth.Vfb + 0.5]]) {
    const hf = fit(hs, Math.round(a * 20) / 20, Math.round(b * 20) / 20);
    assert.ok(hf.N / hp.N - 1 > 0.02 && hf.N / hp.N - 1 < 0.035 && hf.Vfb - hth.Vfb < -0.001 && hf.Vfb - hth.Vfb > -0.007, `1e18: ${hf.N / hp.N - 1}, ${hf.Vfb - hth.Vfb}`);
  }
  const [at] = measure(hp, [1e5], [0.3]), fall = reading(at.Z[0], 1e5, 'series') / hth.C(0.3);
  assert.ok(fall > 0.62 && fall < 0.72, `series at 100 kHz: ${fall}`);
  assert.ok(fit(measure(hp, [1e5]).map((q) => ({ E: q.E, y: 1 / reading(q.Z[0], 1e5, 'parallel') ** 2 })), 0.1, 0.5).slope < 0);
  // holes: none made, the line goes on straight past where inversion would begin (deep depletion); in
  // equilibrium, the simulated 1/C² goes flat (the holes can't follow), where the theory's falls toward 0
  const deep = measure(p, [1e3], [th.Einv + 0.1, th.Einv + 0.4]).map((q) => 1 / reading(q.Z[0], 1e3, 'series') ** 2);
  assert.ok(Math.abs(deep[0] / (1 / th.C(th.Einv + 0.1) ** 2) - 1) < 5e-3 && deep[1] > 1.3 * deep[0], `deep depletion: ${deep}`);
  const ep = { ...p, holes: 'equilibrium' }, eth = theory(ep);
  const inv = measure(ep, [1e3], [eth.Einv + 0.1, eth.Einv + 0.4]).map((q) => reading(q.Z[0], 1e3, 'series'));
  assert.ok(Math.abs(inv[1] / inv[0] - 1) < 0.02 && eth.C(eth.Einv + 0.4) > 10 * inv[1], `${inv} vs ${eth.C(eth.Einv + 0.4)}`);
});

test('demos/redox.html: each cycle against the same experiment worked out apart from the solver, as the page says', async () => {
  const { simulate, cvTheory, peaks, ironBulk } = await pageTheory('redox.html', ['simulate', 'cvTheory', 'peaks', 'ironBulk']);
  const gap = (v, k0, L, b, cycles = 3) => {
    const sim = simulate(v, k0, L, b, cycles), th = cvTheory({ v, k0, L, b, cycles });
    const at = (t) => {
      let k = th.t.findIndex((x) => x >= t);
      if (k < 1) k = t <= th.t[0] ? 1 : th.t.length - 1;
      return th.I[k - 1] + ((th.I[k] - th.I[k - 1]) * (t - th.t[k - 1])) / (th.t[k] - th.t[k - 1]);
    };
    const top = Math.max(...th.I.map(Math.abs));
    return { worst: Math.max(...sim.t.map((t, j) => Math.abs(sim.I[j] - at(t)) / top)), sim: peaks(sim), th: peaks(th) };
  };
  // within 1% of the peak current, at the controls' corners and the review's worst setting (the rest is migration)
  for (const [v, L, k0] of [[1e-3, 10 ** 0.5, 1e-6], [0.1, 1000, 10], [10, 10 ** 1.25, 1e-3], [10, 1000, 10], [0.1, 10 ** 1.25, 1e-3], [1e-2, 10 ** 2.25, 10]]) {
    const g = gap(v, k0, L * 1e-6);
    assert.ok(g.worst < 1e-2, `${v} V/s, ${L} µm, k0 ${k0}: ${g.worst}`);
  }
  // a thick layer with slow kinetics still drifts at the tenth cycle: compared cycle for cycle
  assert.ok(gap(1, 1e-5, 10 ** 2.25 * 1e-6, undefined, 10).worst < 1e-2);
  // with a tenth of the iron, a tenth of the gap at a slow scan: it's migration; at a fast one the
  // simulation's own resolution is left, up to 0.3%
  assert.ok(gap(0.1, 1e-3, 10 ** 1.25 * 1e-6, ironBulk(0.5)).worst < 1.5e-3);
  assert.ok(gap(10, 10, 10 ** 1.25 * 1e-6, ironBulk(0.5)).worst < 3.5e-3);
  // the separation by Nicholson's ψ: 57.5 mV reversible, 57.8, 59.9 and 80 at k⁰ 10, 1, 0.1 (0.1 V/s, no
  // ohmic drop), the same at 10 V/s with ten times the k⁰; 65 mV through 1 mm's ohmic drop at k⁰ 10
  const sep = (v, k0, iR = false) => {
    const p = peaks(cvTheory({ v, k0, L: 1e-3, iR }));
    return p.peaked ? 1000 * (p.Ea - p.Ec) : NaN;
  };
  for (const [k0, want] of [[1000, 57.5], [10, 57.8], [1, 59.9], [0.1, 80]]) assert.ok(Math.abs(sep(0.1, k0) - want) < 0.2, `k0 ${k0}: ${sep(0.1, k0)}`);
  assert.ok(Math.abs(sep(10, 1) - sep(0.1, 0.1)) < 0.05 && Math.abs(sep(0.1, 10, true) - 64.7) < 0.3);
  const still = gap(0.1, 10, 1e-3);
  assert.ok(Math.abs(1000 * (still.sim.Ea - still.sim.Ec) - 1000 * (still.th.Ea - still.th.Ec)) < 1.5);
  // the wave's shape: levelled off (a thin layer), still rising (slow kinetics), or peaked
  assert.equal(peaks(cvTheory({ v: 0.1, k0: 1e-3, L: 10 ** 1.25 * 1e-6 })).shape, 'level');
  assert.equal(peaks(cvTheory({ v: 0.1, k0: 1e-6, L: 10 ** 0.5 * 1e-6 })).shape, 'rising');
  assert.equal(still.th.shape, 'peaks');
});

test('demos/myelin.html: the spikes against the compartmental cables in its script, as the page says', async () => {
  const { Device } = await import('../src/index.js');
  const { cable, myelinated, bare, SETTLE } = await pageTheory('myelin.html', ['cable', 'myelinated', 'bare', 'SETTLE']);
  // when the spike crosses −30 mV at the cable's points, simulated from rest, against the cable
  const arrivals = (def, th) => {
    const d = new Device(def(false));
    d.advance(SETTLE, { tol: 1e-4 });
    d.set(def(true));
    const end = d.solution().x.at(-1), probes = th.x.filter((x) => x <= end).map((x) => ({ x, quantity: 'phi' }));
    const { t, probes: pr } = d.advance(SETTLE + 2e-3, { tol: 1e-5, probes }).trace;
    return probes.map((_, j) => {
      for (let i = 1; i < t.length; i++) if (pr[j][i - 1] < -0.03 && pr[j][i] >= -0.03) return t[i - 1] + ((t[i] - t[i - 1]) * (-0.03 - pr[j][i - 1])) / (pr[j][i] - pr[j][i - 1]) - SETTLE;
      return NaN;
    });
  };
  const last = (a, th) => {
    const k = a.findLastIndex(Number.isFinite);
    return { sim: a[k], cable: th.t[k] };
  };
  // within 1.5 µs at every point (0.477 ms to 1.9 cm at the defaults; 50 wraps, settling slowest), the bare axon within 4 µs
  const worst = (a, th) => Math.max(...a.map((v, j) => (Number.isFinite(v) ? Math.abs(v - th.t[j]) : 0)));
  for (const p of [{ wraps: 200, Lint: 1e-3, gNode: 5 }, { wraps: 50, Lint: 1e-3, gNode: 5 }]) {
    const th = cable(p), a = arrivals((on) => myelinated(p, on), th);
    assert.ok(worst(a, th) < 1.5e-6 && a.filter(Number.isFinite).length > 30, `${p.wraps} wraps: ${worst(a, th)} s`);
    if (p.wraps === 200) assert.ok(Math.abs(last(a, th).sim - 0.477e-3) < 1e-6);
  }
  const tb = cable({}, { bare: true }), b = arrivals(bare, tb);
  assert.ok(worst(b, tb) < 4e-6, `bare: ${worst(b, tb)} s`);
  // the speeds the text quotes (to the node nearest 1 cm, by the cable the simulation follows)
  const speed = (p) => {
    const step = 1e-6 + p.Lint, x = Math.round(0.01 / step) * step + 0.5e-6, th = cable(p, { tEnd: 4e-3 });
    return x / th.t[th.x.findIndex((xx) => Math.abs(xx - x) < 1e-12)];
  };
  const v = (wraps, Lint, gNode) => speed({ wraps, Lint: Lint * 1e-3, gNode });
  const near = (got, want, tol) => assert.ok(Math.abs(got - want) < tol, `${got} m/s for ${want}`);
  near(v(200, 1, 5), 35, 0.5);
  near(v(50, 1, 5), 11, 0.5);
  near(v(800, 1, 5), 79, 0.5);
  assert.ok(Number.isNaN(v(20, 1, 5)));
  near(v(200, 1, 1), 24, 0.5);
  near(v(200, 1, 20), 43, 0.5);
  for (const L of [0.5, 0.75, 1, 1.25]) near(v(200, L, 5), 35, 0.5);
  near(v(200, 0.25, 5), 33, 0.5);
  near(v(200, 3, 5), 32, 0.5);
});

test('demos/axon.html: the spike against the space-clamped Hodgkin–Huxley equations in its script, as the page says', async () => {
  const { Device } = await import('../src/index.js');
  const { axon, hhODE, stimulus, stimOn, SPACES } = await pageTheory('axon.html', ['axon', 'hhODE', 'stimulus', 'stimOn', 'SPACES']);
  // the open sea water: from each pulse, the peak within 1 mV and the upstroke (−30 mV) within 40 µs
  const space = SPACES.open, up = (ts, Vs) => {
    for (let i = 1; i < ts.length; i++) if (Vs[i - 1] < -0.03 && Vs[i] >= -0.03) return ts[i - 1] + ((ts[i] - ts[i - 1]) * (-0.03 - Vs[i - 1])) / (Vs[i] - Vs[i - 1]);
    return NaN;
  };
  for (const p of [{ I: 0.1, period: 25e-3, T: 279.45, na: 1, k: 1, space: 'open' }, { I: 0.1, period: 25e-3, T: 293.15, na: 1, k: 1, space: 'open' }, { I: 1, period: 25e-3, T: 279.45, na: 1, k: 0.2, space: 'open' }]) {
    const probes = [{ x: space.x, quantity: 'phi', region: space.outside }, { x: space.x, quantity: 'phi', region: 'axoplasm' }, ...['m', 'h', 'n'].map((g) => ({ interface: space.face, gate: g }))];
    const d = new Device(axon(p));
    d.advance(space.settle, { tol: 1e-4 });
    d.set({ ports: { stim: { I: stimulus(p) } } });
    const { t, probes: pr } = d.advance(space.settle + 2 * p.period - 1e-3, { tol: 1e-4, dtMax: 5e-5, probes }).trace;
    const k0 = t.findIndex((x, k) => k > 0 && stimOn(x, p) && !stimOn(t[k - 1], p)), V = t.map((_, k) => pr[1][k] - pr[0][k]);
    const ode = hhODE(p, { t0: t[k0], V: V[k0], m: pr[2][k0], h: pr[3][k0], n: pr[4][k0] }, t[k0] + 0.9 * p.period);
    const ts = t.slice(k0), vs = V.slice(k0), peak = Math.max(...vs), odePeak = Math.max(...ode.V);
    assert.ok(peak > -0.02 && Math.abs(peak - odePeak) < 1e-3, `${JSON.stringify(p)}: peak ${peak} vs ${odePeak}`); // (fired; at 20 °C it barely reaches −9 mV)
    assert.ok(Math.abs(up(ts, vs) - up(ode.t, ode.V)) < 40e-6, `${JSON.stringify(p)}: upstroke ${up(ts, vs)} vs ${up(ode.t, ode.V)}`);
  }
});
