// Round 5: check() (which state a solution is of, and its items), polarization() and the junction
// formulas, from Astra's fifth pass. Each test names its item.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer, ohmic, metal, aqueous, half, bath, check, polarization, henderson, planck, IONS } from '../src/kit.js';

const item = (r, name) => r.items.find((i) => i.name === name);
const F = FARADAY;

// A copper wire between two contacts: φ is NaN in the metal, everywhere.
const wire = (V = 0.05) => build({ library: [metal('Cu')], stack: [ohmic(0, ['e-']), layer('Cu', 100e-6), ohmic(V, ['e-'])] });
// A silver cell, Ag | AgNO₃ | Ag, at 50 mV.
const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 });
const agCell = (k0 = 1e-3, hmax = 1e-6) => {
  const face = { reactions: [{ ...silver, k0, alpha: 0.3 }] };
  return build({ library: [aqueous(['Ag+', 'NO3-'], { epsr: 0 }), metal('Ag')], stack: [ohmic(0, ['e-']), layer('Ag', 1e-6), face, layer('water', 20e-6, { c0: { 'Ag+': 10, 'NO3-': 10 } }), face, layer('Ag', 1e-6), ohmic(0.05, ['e-'])], grid: { hmin: 0.2e-6, hmax } });
};
// A Li⁺ conductor with a neutral tracer, held at both ends, and optionally a port.
const bar = ({ V = 0.01, G = 0, shift = 0 } = {}) => ({
  species: [{ name: 'Li+', z: 1, cRef: 1000 }, { name: 'X', z: 0 }],
  materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 }, X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * F, c0: { X: 1 } }],
  contacts: {
    left: { V: shift, terminal: 'Li+', species: { 'Li+': 'equilibrium', X: { type: 'equilibrium', mu: 0 } }, phi: 'neutral' },
    right: { V: shift + V, terminal: 'Li+', species: { 'Li+': 'equilibrium', X: { type: 'equilibrium', mu: 0 } }, phi: 'neutral' },
  },
  ports: G === 0 ? [] : [{ name: 'p', region: 'bar', from: 4e-6, to: 5e-6, V: shift + 0.01, terminal: 'Li+', species: { 'Li+': { type: 'conductance', G, offset: 0 } } }],
  grid: { hmax: 0.5e-6 },
});
// The first solve check() makes (the grid's, then the order check's two), with its result edited.
function editingSolve(which, edit, run) {
  const real = Device.prototype.solve;
  let n = 0;
  Device.prototype.solve = function (...a) {
    const r = real.apply(this, a);
    return ++n === which ? edit(r) : r;
  };
  try {
    return run();
  } finally {
    Device.prototype.solve = real;
  }
}

test('A3: a solution from before a set() of a material is not the present state, though its φ (NaN in a metal) and terminals match', () => {
  const dev = new Device(wire()), old = dev.solve();
  dev.set({ materials: { Cu: { conductor: { conductivity: 1e6 } } } });
  const now = dev.solve();
  assert.ok(Math.abs(old.current / now.current) > 50, `${old.current} vs ${now.current}`); // (60×)
  for (const opts of [{ refine: false }, {}]) {
    const r = check(dev, old, opts);
    assert.equal(r.ok, false, r.text);
    assert.match(item(r, 'solution').summary, /set\(\) to a new definition/);
    assert.ok(!item(r, 'balance') && !item(r, 'grid'), 'nothing else is compared');
  }
  assert.ok(check(dev, now, { refine: false }).ok);
});

test('A3: a copy of a solution (no mark) is compared by its numbers, NaN φ in the metals no longer masking the electrolyte\'s', () => {
  const dev = new Device(agCell()), old = dev.solve();
  dev.set({ materials: { water: { species: { 'Ag+': { D: 10 * dev.def.materials.water.species['Ag+'].D } } } } });
  dev.solve();
  const r = check(dev, { ...old }, { refine: false });
  assert.equal(r.ok, false, r.text);
  assert.match(item(r, 'solution').summary, /its φ differs from the device's by up to 2\.7\de-4 V/);
  assert.ok(check(dev, { ...dev.solution() }, { refine: false }).ok, 'a copy of the present state passes');
});

test('A10: check() refuses a state set() has changed the equations of, given or not, and a cold layout passed as a solution', () => {
  const d = new Device(wire(0.01));
  // (a cold layout, its snapshot given)
  const cold = check(d, d.solution(), { refine: false });
  assert.equal(cold.ok, false);
  assert.match(item(cold, 'converged').summary, /not solved: the device holds its starting layout/);
  d.advance(0.1);
  assert.ok(check(d).ok, check(d).text);
  d.set({ contacts: { right: { V: 1 } } });
  for (const r of [check(d), check(d, d.solution())]) {
    assert.equal(r.ok, false, r.text);
    assert.match(item(r, 'converged').summary, /not solved: the device's state hasn't been solved since a set\(\)/);
  }
  d.step(1e-6);
  assert.ok(check(d).ok, check(d).text);
});

test('A3: a rollback makes its checkpoint\'s solutions the present state again; another device\'s, of another definition, is refused', () => {
  const d = new Device(wire(0.01)), s = d.solve(), cp = d._checkpoint();
  d.set({ contacts: { right: { V: 0.02 } } });
  d.solve();
  assert.equal(check(d, s, { refine: false }).ok, false);
  d._rollback(cp);
  assert.ok(check(d, s, { refine: false }).ok);
  const other = new Device(wire(0.02)).solve();
  const r = check(d, other, { refine: false });
  assert.equal(r.ok, false);
  assert.match(item(r, 'solution').summary, /another device's, of another definition|terminal 'right'/);
  // (a twin, of the same definition, is the same state: compared by its numbers)
  assert.ok(check(d, new Device(wire(0.01)).solve(), { refine: false }).ok);
});

test('A20: a 2× error in a small terminal current beside a large one fails the grid check', () => {
  const d = new Device(bar({ G: 200 })), s = d.solve();
  assert.ok(Math.abs(s.terminals.p.current / s.terminals.right.current) < 1e-5);
  const r = editingSolve(1, (x) => ({ ...x, terminals: { ...x.terminals, p: { ...x.terminals.p, current: 2 * x.terminals.p.current } } }), () => check(d, s, { order: false }));
  assert.equal(item(r, 'grid').ok, false, r.text);
  assert.match(item(r, 'grid').summary, /largest change is p current/);
  assert.ok(check(d, s, { order: false }).ok, 'unedited, it passes');
});

test('A21: a NaN term fails the balance, and a NaN drift the conservation, instead of vanishing', () => {
  const d = new Device(bar()), s = d.solve();
  const nan = (x) => ({ ...x, flux: { ...x.flux, 'Li+': NaN } });
  const r = check(d, { ...s, contacts: { left: nan(s.contacts.left), right: nan(s.contacts.right) } }, { refine: false });
  assert.equal(item(r, 'balance').ok, false, r.text);
  assert.match(item(r, 'balance').summary, /not a number: Li\+'s left contact reads NaN/);
  d.set({ contacts: { right: { V: 0.02 } } });
  const t = d.advance(1e-3);
  assert.ok(check(d, t, { refine: false }).ok);
  const r2 = check(d, { ...t, conservation: t.conservation.map((c) => ({ ...c, drift: NaN })) }, { refine: false });
  assert.equal(item(r2, 'conservation').ok, false, r2.text);
  assert.match(item(r2, 'conservation').summary, /not a number: Li\+'s stretch \(drift NaN\)/);
});

test('A35: what balance can\'t tell from round-off is said beside the ledgers that balance', () => {
  const d = new Device(bar({ V: 0.01 })), s = d.solve();
  const bad = { ...s, contacts: { left: { ...s.contacts.left, flux: { ...s.contacts.left.flux, X: 1e-15 } }, right: { ...s.contacts.right, flux: { ...s.contacts.right.flux, X: 0 } } } };
  const b = item(check(d, bad, { refine: false }), 'balance');
  assert.match(b.summary, /^one ledger sums to zero .*; besides them, X at 1\.00e-15 mol\/\(m²·s\) is 1e-11 of what it could carry and doesn't balance/);
});

test('A36: a face\'s surface species has a ledger, and a surface out of balance fails it', () => {
  const dev = new Device(build({
    library: [aqueous(['H+', 'Cl-'], { epsr: 0 }), metal('Pt')],
    stack: [ohmic(-0.25, ['e-']), layer('Pt', 1e-6),
      { surface: { Hads: { mu0: -20e3, capacity: 1e-5 } }, reactions: [{ equation: 'H+ + e- = Hads', k0: 1e-4, alpha: 0.5 }, { equation: 'Hads + H+ + e- = H2', fixed: { H2: 0 }, k0: 1e-5, alpha: 0.5 }] },
      layer('water', 10e-6), bath({ 'H+': 100, 'Cl-': 100 }, 'Cl-')],
  }));
  const sol = dev.solve(), r = check(dev, sol, { refine: false, order: false });
  const b = item(r, 'balance');
  assert.ok(r.ok, r.text);
  assert.ok(b.details.ledgers.some((l) => l.title === 'Hads at face 0'), b.details.ledgers.map((l) => l.title).join(', '));
  const fake = structuredClone(sol), [r1, r2] = sol.interfaces[0].rates, dr = 0.3 * r2;
  fake.interfaces[0].rates = [r1 + dr, r2 - dr]; // (every volume ledger unchanged)
  assert.equal(item(check(dev, fake, { refine: false, order: false }), 'balance').ok, false);
});

test('A36: a port\'s surface species has a ledger too, over its window', () => {
  const G = 1e-5, lib = aqueous(['H+', 'Cl-'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const dev = new Device(build({ library: [lib], stack: [bath({ 'H+': 100, 'Cl-': 100 }, 'Cl-'), layer('water', 10e-6, { name: 'film' }), {}],
    ports: [{ name: 'pt', region: 'film', V: -0.2, terminal: 'e-', area: 1e3, surface: { Hads: { mu0: -20e3, capacity: G } },
      reactions: [{ equation: 'H+ + e- = Hads', k0: 1e-6, alpha: 0.5 }, { equation: 'Hads + H+ + e- = H2', fixed: { H2: 0 }, k0: 1e-7, alpha: 0.5 }] }],
    grid: { hmin: 1e-6, hmax: 1e-6 } }));
  const sol = dev.solve();
  assert.ok(sol.converged);
  const r = check(dev, sol, { refine: false, order: false }), b = item(r, 'balance');
  assert.ok(r.ok, r.text);
  const l = b.details.ledgers.find((x) => x.title === 'Hads at port pt');
  assert.ok(l && l.terms.length === 2 && Math.abs(l.residual) < 1e-9, JSON.stringify(l));
  const fake = structuredClone(sol);
  fake.ports[0].rates[0] = fake.ports[0].rates[0].map((x) => 1.3 * x);
  fake.ports[0].rates[1] = fake.ports[0].rates[1].map((x) => 0.7 * x);
  assert.equal(item(check(dev, fake, { refine: false, order: false }), 'balance').ok, false);
});

test('A19: a refinement that fails from a state with a history doesn\'t blame the grid from cold', () => {
  const def = {
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D: 1e-11, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: 10e-6, fixedCharge: -100 * F }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'a', region: 'bar', from: 1e-6, to: 3e-6, V: 0.05, area: 1e6, capacitance: { C: 0.2 } }],
    grid: { hmax: 0.25e-6 },
  };
  const d = new Device(def);
  d.solve();
  d.set({ ports: { a: { I: 0 } } });
  d.advance(1e-4);
  const sol = d.solve();
  assert.ok(check(d, sol, { order: false }).ok);
  const r = editingSolve(1, (x) => ({ ...x, converged: false }), () => check(d, sol, { order: false }));
  assert.equal(item(r, 'grid').ok, null, r.text);
  assert.match(item(r, 'grid').summary, /from cold it would hold other amounts/);
});

test('A25: the grid check counts a face\'s diffuse-layer charge, and compares each layer\'s', () => {
  const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 } };
  const d = new Device(build({ library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')], stack: [ohmic(0.5, ['e-']), layer('Pt', 1e-6), stern, layer('water', 10e-6, { c0: { 'K+': 100, 'Cl-': 100 } }), stern, layer('Pt', 1e-6), ohmic(0, ['e-'])], grid: { hmin: 2e-6, hmax: 2e-6 } }));
  const s = d.solve(), g = item(check(d, s, { order: false }), 'grid');
  assert.ok(g.ok, g.summary); // (was: "nothing to compare")
  const layer0 = g.details.diffs.find((x) => x.what === "face 0's diffuse layer charge (right)");
  assert.ok(layer0 && Math.abs(layer0.a / s.interfaces[0].diffuse.right.charge - 1) < 1e-12, JSON.stringify(g.details.diffs));
});

test('A28: in the order check a current change alone can fail', () => {
  const d = new Device(bar({ V: 1e-3 })), s = d.solve();
  // (check's solves: the grid's, then the order check's as given and in reverse)
  const r = editingSolve(3, (x) => ({ ...x, terminals: { ...x.terminals, right: { ...x.terminals.right, current: 1 }, left: { ...x.terminals.left, current: -1 } } }), () => check(d, s));
  assert.equal(item(r, 'order').ok, false, r.text);
  assert.ok(check(d, s).ok, 'unedited, it passes');
});

test('A29: a floating voltage\'s grid change reads the same whatever is 0 V', () => {
  const rel = (shift) => {
    const def = bar({ V: 0, G: 200, shift });
    def.ports[0] = { ...def.ports[0], V: undefined, I: 0 };
    def.regions[0].grid = { hmax: 2.5e-6 };
    def.contacts.right = { ...def.contacts.right, V: shift + 0.01 };
    const d = new Device(def), s = d.solve();
    return item(check(d, s, { order: false }), 'grid').details.diffs.find((x) => x.what === 'p voltage').rel;
  };
  const [a, b] = [rel(0), rel(1)];
  assert.ok(Math.abs(a / b - 1) < 1e-3, `${a} vs ${b}`); // (was 40× apart)
});

test('A22: polarization() refuses a solution from before a set() that rebuilt the device, and takes a later frame of the same definition', () => {
  const dev = new Device(agCell(1e-5)), sol = dev.solve(), Ve = sol.V['e-'][dev.model.grid.regionEnd[0]];
  assert.ok(Math.abs(polarization(dev, sol, { face: 0 }, Ve).reactions[0].rate[0] / sol.interfaces[0].rates[0] - 1) < 1e-12);
  const frame = dev.step(1e-6);
  assert.ok(Number.isFinite(polarization(dev, sol, { face: 0 }, Ve).current[0]), 'a solve or step since is fine');
  assert.ok(Number.isFinite(polarization(dev, frame, { face: 0 }, Ve).current[0]));
  const faster = { reactions: [{ ...silver, k0: 1e-4, alpha: 0.3 }] };
  dev.set({ interfaces: [faster, faster] });
  assert.throws(() => polarization(dev, sol, { face: 0 }, Ve), (e) => e instanceof DeviceError && /sol isn't of the device's present definition \(the device has been set\(\) to a new definition/.test(e.message));
  const regrid = new Device(agCell(1e-5)), s2 = regrid.solve();
  regrid.set({ grid: { hmax: 0.5e-6 } });
  assert.throws(() => polarization(regrid, s2, { face: 1 }, Ve), /present definition/);
  assert.throws(() => polarization(regrid, { ...s2 }, { face: 1 }, Ve), /it has 35 nodes, the device 50/);
});

test('A26: polarization() at a metal–metal face gives the far metal\'s carrier activity 1, as the solver does', () => {
  const d = new Device(build({ library: [metal('Au'), metal('Pt')], stack: [ohmic(0, ['e-']), layer('Au', 1e-6), { species: { 'e-': 'blocked' }, reactions: [{ equation: 'e-(left) = e-(right)', k0: 1e-4 }] }, layer('Pt', 1e-6), ohmic(0.01, ['e-'])] }));
  const s = d.solve(), p = polarization(d, s, { face: 0 }, s.V['e-'][d.model.grid.regionEnd[0]]);
  assert.ok(Math.abs(p.reactions[0].rate[0] / s.interfaces[0].rates[0] - 1) < 1e-9, `${p.reactions[0].rate[0]} vs ${s.interfaces[0].rates[0]}`); // (was NaN)
  assert.ok(Math.abs(p.current[0] / s.current - 1) < 1e-9);
});

test('A11: polarization() at the solved level is the solved rate, a fast couple\'s near-equilibrium remainder too', () => {
  for (const kFe of [1e-2, 1e2, 1e4, 1e5]) {
    const lib = aqueous(['Fe2+', 'H+', 'Cl-'], { epsr: 0 });
    lib.species.push({ name: 'e-', z: -1 });
    const dev = new Device(build({ library: [lib], stack: [bath({ 'Fe2+': 10, 'H+': 10, 'Cl-': 30 }, 'Cl-'), layer('water', 10e-6, { name: 'film' }), {}],
      ports: [{ name: 'iron', region: 'film', I: 0, terminal: 'e-', area: 1e3, reactions: [
        { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: kFe, alpha: 0.5 },
        { equation: '2 H+ + 2 e- = H2', fixed: { H2: 0 }, k0: 1e-9, alpha: 0.5 }] }],
      grid: { hmin: 1e-6, hmax: 1e-6 } }));
    const sol = dev.solve(), port = sol.ports[0], w = 4;
    const fe = polarization(dev, sol, { port: 'iron', x: port.x[w] }, port.V).reactions[0].rate[0];
    assert.ok(Math.abs(fe / port.rates[0][w] - 1) < 1e-9, `k0 ${kFe}: ${fe} vs ${port.rates[0][w]}`); // (was 0.56 off at 1e4, 0 at 1e5)
  }
});

test('A8: with the surface in equilibrium, the equilibrated step passes what the surface balance asks of it (Volmer–Tafel)', () => {
  const ions = {
    species: [{ name: 'H+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'H2', z: 0, cRef: 1000 }],
    materials: { el: { epsr: 0, species: { 'H+': { D: 9e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 }, H2: { D: 5e-9, mu0: 0 } } } },
  };
  // (Volmer at 1e-1, fast enough to be held at equilibrium: at 1e-3, as Tafel, it passed 0.44 of its
  // exchange rate, its coverage 31% under the solve's, and since Neon's R4 3.2 that's refused.)
  const volmer = { equation: 'H+ + e- = H(ads)', k0: 1e-1, alpha: 0.5 }, tafel = { equation: '2 H(ads) = H2', k0: 1e-3, alpha: 0.5 };
  const dev = new Device(build({ T: 298.15, library: [ions, metal('Pt')],
    stack: [ohmic(0.1, ['e-']), layer('Pt', 1e-6), { surface: { 'H(ads)': { mu0: 0, capacity: 1e-5, theta0: 0.1 } }, reactions: [volmer, tafel] }, layer('el', 10e-6), bath({ 'H+': 100, 'Cl-': 100, H2: 1 }, 'Cl-')] }));
  const sol = dev.solve(), c = polarization(dev, sol, { face: 0 }, sol.V['e-'][0], { surface: 'equilibrium' });
  const [v, t] = c.reactions;
  assert.equal(v.equilibrated, true);
  assert.equal(t.equilibrated, undefined);
  assert.ok(Math.abs(v.rate[0] / (2 * t.rate[0]) - 1) < 1e-12, `Volmer ${v.rate[0]}, Tafel ${t.rate[0]}`); // (Volmer was 0)
  // Volmer near equilibrium in the solve itself: the curve is the solved current to its coverage's difference.
  assert.ok(Math.abs(c.current[0] / sol.current - 1) < 1e-2, `${c.current[0]} vs ${sol.current}`);
});

// The binary salt's closed form, φ_R − φ_L = −V_T (D₊ − D₋)/(D₊ + D₋) ln(c_R/c_L).
const VT = (GAS_CONSTANT * 298.15) / F, Dp = IONS['Na+'].D, Dm = IONS['Cl-'].D;
const binary = (cL, cR) => -VT * ((Dp - Dm) / (Dp + Dm)) * Math.log(cR / cL);
const salt = (c) => ({ 'Na+': c, 'Cl-': c });

test('A1: planck() against the binary closed form to 1 nV at 1e2, 1e4, 1e6 and 1e12 to 1, both ways round', () => {
  for (const ratio of [1e2, 1e4, 1e6, 1e12]) {
    for (const [cL, cR] of [[100, 100 / ratio], [100 / ratio, 100]]) {
      const v = planck(salt(cL), salt(cR), IONS), exact = binary(cL, cR);
      assert.ok(Math.abs(v - exact) < 1e-9, `${cL} | ${cR}: ${v * 1e3} mV vs ${exact * 1e3}`); // (100 mM | 0.1 mM was −274 mV for −73.6)
    }
  }
});

test('A1: planck() says when its steps don\'t resolve the junction', () => {
  assert.throws(() => planck(salt(100), salt(1e-4), IONS, { steps: 8 }), (e) => e instanceof DeviceError && /isn't resolved: with 8 steps .* with 4 .*pass more steps/.test(e.message));
});

test('A46, A54, A53: a neutral solute is left out, T and steps are checked, and henderson holds at extreme ratios', () => {
  const ions = { 'A+': { z: 1, D: 2e-9 }, 'B-': { z: -1, D: 1e-9 }, N: { z: 0, D: 1e-9 } };
  const exact = -VT * (1 / 3) * Math.log(10 / 100);
  for (const f of [planck, henderson]) assert.ok(Math.abs(f({ 'A+': 100, 'B-': 100, N: 1 }, { 'A+': 10, 'B-': 10, N: 2 }, ions) - exact) < 1e-12, f.name);
  for (const T of [0, -300, NaN, Infinity]) {
    assert.throws(() => henderson(salt(100), salt(10), IONS, { T }), /T is the temperature in K, a positive number/);
    assert.throws(() => planck(salt(100), salt(10), IONS, { T }), /T is the temperature in K, a positive number/);
  }
  for (const steps of [0, 1.5, -5, NaN]) assert.throws(() => planck(salt(100), salt(10), IONS, { steps }), /steps is the number of integration steps .* a whole number ≥ 1/);
  for (const lo of [1e-12, 1e-15, 1e-16]) {
    assert.ok(Math.abs(henderson(salt(100), salt(lo), IONS) / binary(100, lo) - 1) < 1e-12, `100 | ${lo}`); // (−∞ at 1e-15)
    assert.ok(Math.abs(henderson(salt(lo), salt(100), IONS) / binary(lo, 100) - 1) < 1e-12);
  }
});
