import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import {
  AgAgCl, H2O, IONS, alloy, aqueous, bath, build, combine, half, hodgkinHuxley, injector, layer, metal, ohmic, photogeneration, pulse, ramp, semiconductor, trap,
} from '../src/kit.js';

// Faults in the kit's definition builders found by round 5's outside reading pass (OpenAI's Astra
// on the kit), each confirmed by a chaser's probe before it was fixed.

test('a fixed μ for the electrons is refused by half() and by a hand-made bath reference', () => {
  // Accepted, bath() summed it into its offset: 96485 J/mol moved the terminal by 1 V (0.512950
  // to −0.487046 V), while level() skipped it and a face reaction refused it.
  assert.throws(() => half('AgCl(s) + e- = Ag(s) + Cl-', { ...AgAgCl.fixed, 'e-': 96485 }), (e) => e instanceof DeviceError && /half: fixed\.e-: 'e-' is a species of the device/.test(e.message));
  const c = { 'K+': 500, 'Cl-': 500 };
  assert.throws(() => bath(c, { equation: AgAgCl.equation, fixed: { ...AgAgCl.fixed, 'e-': 96485 } }), /bath: the reference electrode's fixed\.e-/);
  assert.ok(Math.abs(bath(c, AgAgCl).bath.offset - 1.1378828) < 1e-6); // (the electrode itself unchanged)
});

test('bath(c, undefined, drive, opts) keeps its drive and its options', () => {
  const c = { 'Na+': 10, 'Cl-': 10 };
  assert.deepEqual(bath(c, undefined, 0, { material: 'water' }), bath(c, 0, { material: 'water' }));
  assert.equal(bath(c, undefined, 0.05, { material: 'water' }).V, 0.05);
  assert.deepEqual(bath(c, undefined, { I: 0 }, { material: 'water' }).bath.material, 'water');
  assert.equal(bath(c, undefined, { I: 0 }).I, 0);
  // The consequence: a membrane between two water baths. With the bath's own material the
  // membrane's end holds the partitioned 0.485 mol/m³; dropped, it held the bath's 10.
  const species = [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }];
  const materials = {
    water: { epsr: 78, species: { 'Na+': { D: 1.3e-9, cRef: 1000, mu0: 0 }, 'Cl-': { D: 2e-9, cRef: 1000, mu0: 0 } } },
    mem: { epsr: 10, species: { 'Na+': { D: 1e-11, cRef: 1000, mu0: 5e3 }, 'Cl-': { D: 1e-11, cRef: 1000, mu0: 1e4 } } },
  };
  const s = new Device(build({ species, materials, stack: [bath(c, undefined, 0, { material: 'water' }), layer('mem', 1e-7), bath({ 'Na+': 1, 'Cl-': 1 }, 0, { material: 'water' })], grid: { minCells: 40 } })).solve();
  assert.ok(s.converged);
  assert.ok(Math.abs(s.c['Na+'][0] / 0.48535 - 1) < 1e-3, `c(Na+) at the membrane's end ${s.c['Na+'][0]}`);
});

test("layer() and alloy() given another layer as a template take their own material and length", () => {
  const old = layer('Si', 1e-6, { name: 'abs', donors: 1 });
  assert.deepEqual(layer('GaAs', 2e-6, old), { name: 'abs', donors: 1, material: 'GaAs', length: 2e-6 });
  const a = alloy('GaAs', 'AlAs', 2e-6, old);
  assert.equal(a.material, 'GaAs');
  assert.equal(a.length, 2e-6);
  assert.equal(a.alloy.to, 'AlAs');
  assert.equal(alloy('GaAs', 'InAs', 3e-6, a).alloy.to, 'InAs');
});

test('two build()s from one library share no material objects', () => {
  const lib = semiconductor('Si'), D0 = lib.materials.Si.species['e-'].D;
  const def = () => build({ library: [lib], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1 }), ohmic(0)] });
  const ref = def(), cand = def();
  assert.notEqual(cand.materials.Si, lib.materials.Si);
  cand.materials.Si.species['e-'].D = 2 * D0;
  assert.equal(ref.materials.Si.species['e-'].D, D0);
  assert.equal(lib.materials.Si.species['e-'].D, D0);
  assert.deepEqual(combine(lib).materials, lib.materials); // (a copy, alike)
});

test('a material given twice must be alike in every defined key, whatever the order', () => {
  // A's mu0T: undefined against B's statistics: over A's keys alone, alike, and B's Fermi–Dirac
  // statistics were dropped when A came first.
  const base = semiconductor('Si');
  const A = { materials: { Si: { ...base.materials.Si, mu0T: undefined } } };
  const B = { materials: { Si: { ...base.materials.Si, statistics: [{ type: 'fermi-dirac', species: ['e-'] }] } } };
  for (const parts of [[A, B], [B, A]]) assert.throws(() => combine({ species: base.species }, ...parts), /material 'Si' is given twice, differently/);
  // and an undefined key is still no difference
  assert.deepEqual(Object.keys(combine({ species: base.species }, A, base).materials.Si).filter((k) => A.materials.Si[k] !== undefined), ['T', 'epsr', 'species']);
});

test("an alloy's end without its own cRef grades from the species' cRef", () => {
  const species = [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }];
  const materials = {
    A: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } },
    B: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0, cRef: 500 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } },
  };
  const def = build({ species, materials, stack: [bath({ 'Na+': 10, 'Cl-': 10 }), alloy('A', 'B', 1e-6, { dipole: 0, c0: { 'Na+': 10, 'Cl-': 10 } }), bath({ 'Na+': 10, 'Cl-': 10 })] });
  assert.deepEqual(def.regions[0].grade['Na+'].cRef.values, [1000, 500]); // (was NaN, NaN)
  assert.equal(def.regions[0].grade['Cl-'], undefined);
  assert.ok(new Device(def).solve().converged);
});

test('a second trap of the same name in a material is refused; in another material it shares the species', () => {
  // Merged, the second's density replaced the first's in one pool: p 1.0989 times the two traps'.
  const Si = semiconductor('Si');
  const t1 = trap(Si, { density: 1e-6, tauN: 2e-6, tauP: 1e-6 });
  assert.throws(() => trap(t1.library, { density: 1e-8, tauN: 1e-7, tauP: 1e-7 }), (e) => e instanceof DeviceError && /material 'Si' already has 'X0'.*give the second trap a name/.test(e.message));
  const t2 = trap(t1.library, { density: 1e-8, tauN: 1e-7, tauP: 1e-7, name: 'Y' });
  assert.deepEqual(t2.library.species.map((s) => s.name), ['e-', 'h+', 'X0', 'X-', 'Y0', 'Y-']);
  const two = combine(Si, semiconductor('Ge'));
  const a = trap(two, { density: 1e-6, tauN: 1e-6, tauP: 1e-6, material: 'Si' });
  const b = trap(a.library, { density: 2e-6, tauN: 1e-6, tauP: 1e-6, material: 'Ge' });
  assert.deepEqual(b.library.species.map((s) => s.name), ['e-', 'h+', 'X0', 'X-']);
  assert.equal(b.library.materials.Si.species.X0.cRef, 1e-6);
  assert.equal(b.library.materials.Ge.species.X0.cRef, 2e-6);
  // (a species of the name with another charge is refused)
  const odd = { species: [...Si.species, { name: 'X-', z: 1 }], materials: Si.materials };
  assert.throws(() => trap(odd, { density: 1e-6, tauN: 1e-6, tauP: 1e-6 }), /already has a species 'X-' of charge 1/);
});

test("trap c0For splits the total by the stable roots, the minority kept at any doping", () => {
  // n = net/2 + √(net²/4 + n_i²) lost the minority to cancellation: at 1e18 cm⁻³ acceptors X⁻ was
  // 0 (refused as the user's c0), at 1e17 0.6% low; and 1 − f, the empty fraction at heavy n-type.
  const T = 300, Si = semiconductor('Si'), { 'e-': el, 'h+': ho } = Si.materials.Si.species;
  const ni = Math.sqrt(el.cRef * ho.cRef * Math.exp(-(ho.mu0 + el.mu0) / (GAS_CONSTANT * T)));
  const Nt = 1e-6, t = trap(Si, { density: Nt, tauN: 1e-6, tauP: 1e-6 });
  for (const NAcm of [1e14, 1e16, 1e17, 1e18]) {
    const NA = units.perCm3(NAcm), n = (2 * ni * ni) / (NA + Math.sqrt(NA * NA + 4 * ni * ni)); // (n p = n_i², p = NA + n)
    const got = t.c0For({ acceptors: NA });
    assert.ok(Math.abs(got['X-'] / ((Nt * n) / (n + ni)) - 1) < 1e-12, `${NAcm}: X- ${got['X-']}`);
    assert.ok(Math.abs((got['X-'] + got.X0) / Nt - 1) < 1e-15);
  }
  const n1 = ni * Math.exp((-0.46 * FARADAY) / (GAS_CONSTANT * T)), deep = trap(Si, { density: 1e-3, tauN: 1e-6, tauP: 1e-6, n1 });
  for (const NDcm of [1e16, 1e18]) {
    const ND = units.perCm3(NDcm), n = ND / 2 + Math.sqrt((ND * ND) / 4 + ni * ni);
    assert.ok(Math.abs(deep.c0For({ donors: ND }).X0 / ((1e-3 * n1) / (n + n1)) - 1) < 1e-12); // (was 11% low at 1e18)
  }
  const NA = units.perCm3(1e18), ND = units.perCm3(1e16);
  const def = build({ T, library: [t.library], stack: [ohmic(0), layer('Si', 1e-6, { acceptors: NA, c0: t.c0For({ acceptors: NA }) }), layer('Si', 1e-6, { donors: ND, c0: t.c0For({ donors: ND }) }), ohmic(0)], bulkReactions: t.bulkReactions });
  assert.ok(new Device(def).solve().converged);
});

test("hodgkinHuxley refuses a leak carried by a gated channel's ion", () => {
  // The leak's link replaced the channel's: K⁺ at 3 S/m², ungated, where 360 n⁴ + 3 is right.
  const inside = { 'Na+': 50, 'K+': 400, 'Cl-': 40 }, outside = { 'Na+': 440, 'K+': 20, 'Cl-': 560 };
  for (const leak of ['K+', 'Na+']) {
    assert.throws(() => hodgkinHuxley({ T: 300, linear: true, leak, leakZ: 1 }), /can't be both gated and ungated/);
    assert.throws(() => hodgkinHuxley({ T: 279.45, inside, outside, leak, leakZ: 1 }), /can't be both gated and ungated/);
  }
  assert.throws(() => hodgkinHuxley({ T: 300, linear: true, names: { Na: 'K+' } }), /names\.Na and names\.K/);
  assert.deepEqual(Object.keys(hodgkinHuxley({ T: 300, linear: true }).species), ['Na+', 'K+', 'Cl-']);
});

test("a bulk rate keyed to an alloy's far material alone is refused, as it doesn't reach the layer", () => {
  // Generation keyed to B in a layer alloy('Si', 'B', …) at fraction 1: the layer is Si, and it
  // generated nothing there (total 1.66e-4 where 3.32e-4 is right), converged.
  const L = 1e-6, G = units.perCm3(1e20), ND = units.perCm3(1e16);
  const lib = [semiconductor('Si'), semiconductor('Si', { material: 'B' })];
  const mk = (kf) =>
    build({
      T: 300,
      library: lib,
      stack: [ohmic(0), layer('B', L, { donors: ND }), { step: { species: 'e-', value: 0 } }, alloy('Si', 'B', L, { fraction: [1, 1], step: { species: 'e-', value: 0 }, donors: ND }), ohmic(0)],
      bulkReactions: [{ equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf }, { equation: 'e- + h+ = 0', kf: { B: 1e9, Si: 1e9 } }],
      grid: { minCells: 40 },
    });
  assert.throws(() => mk({ B: G }), (e) => e instanceof DeviceError && /bulkReactions\[0\]\.kf\.B: stack\[3\] is alloy\('Si', 'B'\), a layer of 'Si' throughout/.test(e.message));
  assert.throws(() => mk({ B: G, Si: G / 2 }), /it runs at kf\.Si there/);
  const sol = new Device(mk({ B: G, Si: G })).solve();
  assert.ok(sol.converged);
  assert.ok(Math.abs(sol.bulkReactions[0].total / (2 * L * G) - 1) < 1e-6);
});

test('a ramp or a pulse edge below the clock\'s resolution at its start is refused', () => {
  // At t = 1e9 s (resolution ~1e-7 s) a 10 ns ramp was a jump, and a pulse's 1 µs width 0.95 µs.
  assert.throws(() => ramp({ from: 0, to: 1, duration: 1e-8, start: 1e9 }), (e) => e instanceof DeviceError && /ramp: the duration.*time's resolution.*device's time 0/.test(e.message));
  assert.throws(() => pulse({ start: 1e9, width: 1e-6, rise: 1e-8, amplitude: 1 }), /pulse: the width|pulse: the rise/);
  assert.throws(() => pulse({ start: 1e9, width: 1, rise: 1e-8, amplitude: 1 }), /pulse: the rise/);
  assert.deepEqual(ramp({ from: 0, to: 1, duration: 1e-8 }).t, [0, 1e-8]);
  assert.deepEqual(ramp({ from: 0, to: 1, duration: 1, start: 1e9 }).t, [1e9, 1e9 + 1]);
  assert.deepEqual(pulse({ start: 1e9, width: 10, amplitude: 1 }).t, [1e9, 1e9, 1e9 + 10, 1e9 + 10]);
});

test('step() over a breakpoint of a drive says it crosses it unlanded', () => {
  // One step of 3 ms over a 1 ms pulse of 1 A/m² passed no charge at all, without a word.
  const wire = () => new Device(build({ library: [metal('Cu')], stack: [ohmic(0, ['e-']), layer('Cu', 100e-6), ohmic({ I: pulse({ start: 1e-3, width: 1e-3, amplitude: 1 }) }, ['e-'])] }));
  let d = wire();
  d.solve();
  const s = d.step(3e-3);
  assert.ok(s.converged);
  assert.equal(d.solver.termCharge[1], 0); // (the stated limit: the pulse is missed)
  assert.equal(s.warnings.length, 1);
  assert.match(s.warnings[0], /step\(0\.003\): the drive of terminal 'right' has a breakpoint at t = 0\.001 s, inside this step.*advance\(tEnd\)/);
  d = wire();
  d.solve();
  for (let k = 0; k < 3; k++) assert.deepEqual(d.step(1e-3).warnings, []); // (landing on each: no word, and the charge passes)
  assert.ok(Math.abs(Math.abs(d.solver.termCharge[1]) - 1e-3) < 1e-15);
});

test("an injector's conductance counts the window's cross-section", () => {
  // With a geometry its I is a current (A): planar-calibrated, a 1e-4 m² device put the port
  // 1000 V above the window where 0.1 was asked.
  const L = 1e-6, A = 1e-4;
  const port = injector({ region: 'bar', from: 0, to: L, species: 'e-', I: 1e-8, headroom: 0.1, area: A });
  assert.equal(port.species['e-'].G, 1e-8 / (L * A * 0.1));
  const def = build({ T: 300, library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', L, { name: 'bar', donors: units.perCm3(1e16) }), ohmic(0)], ports: [port], grid: { minCells: 40 }, geometry: { area: { x: [0, L], values: [A, A] } } });
  const s = new Device(def).solve();
  assert.ok(s.converged);
  assert.ok(Math.abs(s.terminals.injector.V / 0.1 - 1) < 1e-6, `port at ${s.terminals.injector.V} V`);
  assert.throws(() => injector({ region: 'bar', from: 0, to: L, species: 'e-', I: 1, area: 0 }), /area must be the window's cross-section/);
});

test("photogeneration's guard outside its span is a negligible share at any absorption depth", () => {
  // A hair of 1e-10 of the span leaked G_entry·hair/2 = 5e-5 of the flux at 1e6 absorption lengths.
  const area = ({ x, values }, a, b) => {
    let s = 0;
    for (let k = 1; k < x.length; k++) {
      const lo = Math.max(a, x[k - 1]), hi = Math.min(b, x[k]);
      if (hi <= lo) continue;
      const f = (t) => values[k - 1] + ((values[k] - values[k - 1]) * (t - x[k - 1])) / (x[k] - x[k - 1]);
      s += ((hi - lo) * (f(lo) + f(hi))) / 2;
    }
    return s;
  };
  for (const [alpha, L, toward] of [[1e8, 1e-2, '+x'], [1e8, 1e-2, '-x'], [1e6, 1e-6, '+x'], [1e3, 1e-9, '+x']]) {
    const from = 1e-3, to = from + L, kf = photogeneration({ material: 'm', flux: 1, alpha, mu: 3e5, from, to, toward }).kf.m;
    const outside = area(kf, -Infinity, from) + area(kf, to, Infinity);
    assert.ok(outside < 1e-9 * (1 - Math.exp(-alpha * L)), `${alpha} ${L} ${toward}: ${outside}`);
    assert.ok(kf.x.every((x, k) => k === 0 || x > kf.x[k - 1]));
  }
});

test("the warning on μ° from another temperature says what the device's constants are", () => {
  // At 37 °C with 25 °C's μ°, Kw comes out e^(−ΔG°/RT) = 3.5e-14: neither 25 °C's nor 37 °C's.
  const T = 310.15, w = aqueous(['H+', 'OH-', 'Na+', 'Cl-'], { epsr: 0, T, D: { 'H+': 9.3e-9, 'OH-': 5.3e-9 } });
  const blocked = { 'H+': 'blocked', 'OH-': 'blocked', 'Na+': 'blocked', 'Cl-': 'blocked' };
  const dev = new Device({
    T, ...w,
    regions: [{ name: 'a', material: 'water', length: 1e-6, c0: { 'H+': 1e-3, 'OH-': 1e-3, 'Na+': 100, 'Cl-': 100 } }],
    bulkReactions: [{ equation: 'H2O = H+ + OH-', fixed: { H2O }, kf: { water: 1e-3 } }],
    contacts: { left: { phi: 'neutral', V: 0, terminal: 'Cl-', species: { ...blocked, 'Cl-': 'equilibrium' } }, right: { phi: 'neutral', species: blocked } },
    grid: { minCells: 8 },
  });
  const s = dev.solve();
  const [warning] = s.warnings.filter((x) => /its μ° are for 298\.15 K/.test(x));
  assert.match(warning, /e\^\(−ΔG°\/RT\) with 298\.15 K's ΔG° at 310\.15 K.*neither 298\.15 K's constant nor 310\.15 K's/);
  const dG = IONS['H+'].mu0 + IONS['OH-'].mu0 - H2O, Kw = (s.c['H+'][4] / 1000) * (s.c['OH-'][4] / 1000);
  assert.ok(Math.abs(Kw / Math.exp(-dG / (GAS_CONSTANT * T)) - 1) < 1e-6, `Kw ${Kw}`);
});
