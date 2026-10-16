import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, units } from '../src/index.js';
import { build, layer, ohmic, bath, semiconductor, photogeneration, perovskiteCell, check } from '../src/kit.js';

// An n⁺p silicon cell under light (Beer–Lambert), graded from 1 nm at each end to `hmax`.
const Si = semiconductor('Si'), W = units.um(60.5), flux = 3e-3, alpha = 1e5;
const solar = (V, { hmin = units.nm(1), hmax = units.um(1), ratio = 1.2 } = {}) =>
  build({
    T: 300,
    library: [Si],
    stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e19) }), layer('Si', units.um(60), { acceptors: units.perCm3(1e16) }), ohmic(V)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }, photogeneration({ material: 'Si', flux, alpha, mu: units.eV(3), to: W })],
    grid: { hmin, hmax, ratio },
  });
const item = (report, name) => report.items.find((it) => it.name === name);

test('check(): the ledger of a lit cell is J = F(G − R), term by term, and an equilibrium has none', () => {
  const dev = new Device(solar(0.5)), sol = dev.solve(), report = check(dev, sol, { refine: false });
  const { ledgers, worst } = item(report, 'balance').details;
  assert.ok(worst < 1e-12, `sums to ${worst}`);
  const e = ledgers.find((l) => l.species === 'e-'), term = (what) => e.terms.find((t) => t.what === what).rate;
  // Every photon absorbed makes a pair: Φ(1 − e^{−αW}).
  assert.ok(Math.abs(term('photon = e- + h+') / (flux * -Math.expm1(-alpha * W)) - 1) < 1e-9);
  // What leaves through the contacts is the current, and the rest recombined.
  const out = -term('left contact') - term('right contact'), hOut = ledgers.find((l) => l.species === 'h+').terms.filter((t) => /contact/.test(t.what)).reduce((a, t) => a - t.rate, 0);
  assert.ok(Math.abs(FARADAY * (-term('left contact') + 0) - FARADAY * -term('left contact')) === 0);
  assert.ok(Math.abs((FARADAY * (term('photon = e- + h+') + term('e- + h+ = 0') - out)) / sol.current) < 1e-9);
  assert.ok(out > 0 && hOut > 0);
  assert.match(report.text, /^ok {3}balance: 2 ledgers sum to zero/m);
  assert.match(report.text, /e-: photon = e- \+ h\+: \+2\.99e-3 \(289 A\/m²\)/);

  const dark = new Device(build({ T: 300, library: [Si], stack: [ohmic(0), layer('Si', 1e-6, { donors: units.perCm3(1e17) }), layer('Si', 1e-6, { acceptors: units.perCm3(1e16) }), ohmic(0)] }));
  const eq = check(dark, dark.solve());
  assert.ok(eq.ok, eq.text);
  assert.equal(item(eq, 'balance').summary, 'nothing flows or reacts (equilibrium)');
  assert.ok(item(eq, 'grid').ok, 'round-off currents are not compared');
});

test('check(): two reactions with one equation (radiative and SRH recombination) keep a ledger line each', () => {
  const def = solar(0.5);
  def.bulkReactions = [def.bulkReactions[0], { equation: 'e- + h+ = 0', srh: { Si: { tauN: 1e-6, tauP: 1e-6 } } }, def.bulkReactions[1]];
  const dev = new Device(def), sol = dev.solve(), report = check(dev, sol, { refine: false });
  const e = item(report, 'balance').details.ledgers.find((l) => l.species === 'e-');
  const term = (what) => e.terms.find((t) => t.what === what)?.rate;
  const total = (k) => sol.bulkReactions[k].regions.reduce((a, r) => a + r, 0);
  assert.ok(Math.abs(term('e- + h+ = 0 (bulkReactions[0])') / -total(0) - 1) < 1e-12);
  assert.ok(Math.abs(term('e- + h+ = 0 (bulkReactions[1])') / -total(1) - 1) < 1e-12);
  assert.equal(term('e- + h+ = 0'), undefined);
  assert.equal(term('photon = e- + h+') > 0, true, 'a reaction with its own equation keeps its plain name');
});

test('check(): the grid check finds a coarse grid, and its change estimates the error (second order: ¾ of it)', () => {
  const dev = new Device(solar(0.5)), report = check(dev, dev.solve());
  const grid = item(report, 'grid');
  assert.equal(grid.ok, false);
  assert.equal(report.ok, false);
  assert.match(grid.summary, /largest change is region 0 \(Si\) charge, .* over 0\.01: refine the grid/);
  // Against a grid 16 times as fine, as good as exact here: the coarse grid's error, and the change.
  const charge0 = (s) => grid.details.diffs.find((d) => /region 0/.test(d.what)) && s;
  const exact = new Device(solar(0.5, { hmin: units.nm(1) / 16, hmax: units.um(1) / 16, ratio: 1.2 ** (1 / 16) })).solve();
  const q = (s) => {
    let sum = 0;
    for (let g = 0; g < s.x.length; g++) {
      if (s.region[g] !== 0) continue;
      const h = (g > 0 && s.region[g - 1] === 0 ? s.x[g] - s.x[g - 1] : 0) + (s.region[g + 1] === 0 ? s.x[g + 1] - s.x[g] : 0);
      sum += (FARADAY * (units.perCm3(1e19) + s.c['h+'][g] - s.c['e-'][g]) * h) / 2;
    }
    return sum;
  };
  assert.ok(charge0(exact));
  const d = grid.details.diffs.find((x) => /region 0/.test(x.what)), error = d.a - q(exact), change = d.a - d.b;
  assert.ok(change / error > 0.65 && change / error < 0.85, `change ${change} against error ${error}`);
  // Finer, it passes; and a benchmark's tolerance can still ask for more.
  const finer = new Device(solar(0.5, { hmin: units.nm(0.25), hmax: units.um(0.25), ratio: 1.2 ** 0.25 }));
  const ok = check(finer, finer.solve());
  assert.ok(ok.ok, ok.text);
  assert.equal(check(finer, finer.solve(), { tol: 1e-4 }).ok, false);
});

test('check(): a floating terminal compares its voltage (V_oc); a strictly neutral junction has no charge to compare', () => {
  const dev = new Device(perovskiteCell({ V: { I: 0 } })), report = check(dev, dev.solve());
  assert.ok(report.ok, report.text);
  const voc = item(report, 'grid').details.diffs.find((d) => d.what === 'right voltage');
  assert.ok(voc && voc.a > 1.05 && Math.abs(voc.b - voc.a) < 1e-4, JSON.stringify(voc));
  const e = item(report, 'balance').details.ledgers.find((l) => l.species === 'e-');
  assert.ok(e.terms.every((t) => !/contact/.test(t.what) || Math.abs(t.rate) < 1e-12 * Math.max(...e.terms.map((u) => u.rate))), 'at open circuit nothing leaves');

  const junction = new Device(
    build({
      species: [
        { name: 'Na+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } } },
      stack: [bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-'), layer('water', 100e-6), bath({ 'Na+': 10, 'Cl-': 10 }, 'Cl-', { I: 0 })],
      grid: { hmin: 1e-6, hmax: 5e-6 },
    }),
  );
  const j = check(junction, junction.solve());
  assert.ok(j.ok, j.text);
  assert.ok(item(j, 'grid').details.diffs.every((d) => !/charge/.test(d.what)));
});

test('check(): a transient checks conservation; a failed solve and warnings are reported', () => {
  const dev = new Device(perovskiteCell({ V: 1.2 }));
  dev.solve();
  dev.set({ contacts: { right: { V: 1.0 } } });
  const s = dev.advance(1, { tol: 1e-4 }), report = check(dev, s);
  assert.ok(report.ok, report.text);
  assert.deepEqual(report.items.map((it) => it.name), ['converged', 'warnings', 'conservation', 'grid']);
  assert.match(item(report, 'conservation').summary, /^1 closed stretch .* kept to /);
  // (its grid isn't checked, and it says so rather than passing it)
  assert.equal(item(report, 'grid').ok, null);
  assert.match(item(report, 'grid').summary, /not checked in a transient/);

  const failed = check(dev, { ...s, converged: false });
  assert.equal(failed.ok, false);
  assert.deepEqual(failed.items.map((it) => it.name), ['converged', 'warnings']);
  const slip = solar(0);
  slip.materials.Si.species['e-'].D = 36; // cm²/s, not m²/s
  const slipped = new Device(slip), w = check(slipped, slipped.solve(), { refine: false });
  assert.equal(item(w, 'warnings').ok, null);
  assert.ok(w.ok, 'a warning asks you to look; it is not a failure');
  assert.match(w.text, /^\? {4}warnings: 1: materials\.Si\.species\.e-\.D/m);
});

// A contact holding one ion at a level far from its neighbours' (Na⁺ at −0.04 V beside a KCl
// bath) holds it at 1e37 mol/m³: denser than any material, and said so. (Its currents in a
// transient had read 1e28 A/m², with no warning.)
test('a species denser than any material is warned of', () => {
  const sol = new Device({
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'Na+', z: 1, cRef: 1000 }],
    materials: { w: { epsr: 0, species: { 'K+': { D: 1.96e-9, mu0: -283e3 }, 'Cl-': { D: 2.03e-9, mu0: -131e3 }, 'Na+': { D: 1.33e-9, mu0: -262e3 } } } },
    regions: [{ material: 'w', length: 50e-6, c0: { 'K+': 60, 'Cl-': 80, 'Na+': 20 } }],
    contacts: { left: { bath: { c: { 'K+': 80, 'Cl-': 80 }, reference: 'Cl-' }, V: 0 }, right: { terminal: 'Na+', species: { 'Na+': 'equilibrium' }, phi: 'bulk', V: -0.04 } },
    grid: { hmin: 0.3e-6, hmax: 3e-6 },
  }).solve();
  assert.ok(sol.warnings.some((w) => /Na\+ reaches .* denser than any material/.test(w)), sol.warnings.join('\n'));
});

// A steady solve that depends on the order the species are listed in has converged somewhere
// that isn't steady. Round 11's fuzzer found it by exactly this test (its F5: an electrode at open
// circuit with a fast couple and a slow one, fed through a slow transfer at a face, settling
// 20–80 mV apart by the listing, each converged without a warning); check() now makes it too.
test('check() solves the device again with its species listed in reverse, and fails where that differs', () => {
  // A solve made to land 1 mV off when the species come reversed stands in for one that isn't
  // steady (the fuzzer's devices that did it are fixed).
  const dev = new Device(solar(0.5)), first = dev.def.species[0].name, solve = Device.prototype.solve;
  Device.prototype.solve = function (...args) {
    const sol = solve.apply(this, args);
    if (this.def.species[0].name !== first) for (const t of Object.values(sol.terminals)) t.V += 1e-3;
    return sol;
  };
  let item;
  try {
    item = check(dev, dev.solve(), { tol: 1 }).items.find((it) => it.name === 'order');
  } finally {
    Device.prototype.solve = solve;
  }
  assert.ok(item && item.ok === false && /different steady state/.test(item.summary), item?.summary);
  // A device that solves the same either way passes it.
  const fine = new Device(solar(0.5));
  const ok = check(fine, fine.solve(), { tol: 1 }).items.find((it) => it.name === 'order');
  assert.ok(ok?.ok === true, ok?.summary);
});

// The fuzzer's F5 itself: 4 mV apart by the listing (0.274922 and 0.271016 V) while the contact
// and flat-stretch rows added a level's low part to its high part before subtracting the target,
// which rounds the low part away. They subtract first now, and it's the same either way.
test('an electrode at open circuit with a fast couple and a slow one, fed through a slow face: the same steady state in either listing', () => {
  const water = { epsr: 0, species: { 'B-': { D: 1e-9, mu0: 0 }, 'C2+': { D: 1e-9, mu0: 0 } } };
  const c0 = { 'B-': 474.7252163975241, 'C2+': 237.36260819876205 };
  const d = new Device({
    T: 298.15,
    species: [{ name: 'e-', z: -1 }, { name: 'B-', z: -1, cRef: 1000 }, { name: 'C2+', z: 2, cRef: 1000 }],
    materials: { water, water2: water, M: { conductor: { species: 'e-', conductivity: 1523338.8768275294 } } },
    regions: [{ name: 'M1', material: 'M', length: 1e-6 }, { name: 'water', material: 'water', length: 1e-6, c0 }, { name: 'water2', material: 'water2', length: 1e-6, c0 }],
    interfaces: [
      { reactions: [{ equation: 'B- = C2+ + 3 e-', k0: 1, fixed: { 'S00B(s)': 0 } }, { equation: 'B- = e- + S00B(s)', k0: 1e-20, fixed: { 'S00B(s)': -4108 } }] },
      { species: { 'B-': 'blocked' }, reactions: [{ equation: 'B-(left) = B-(right)', k0: 1e-17 }] },
    ],
    contacts: { left: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0.3, bath: { c: { 'B-': 21.62350860406672, 'C2+': 10.81175430203336 }, reference: 'B-' } } },
    grid: { hmin: 1e-7, hmax: 3e-7 },
  });
  const report = check(d, d.solve()), item = report.items.find((it) => it.name === 'order');
  assert.ok(item?.ok === true, item?.summary);
});

test('check(): a diode passing a small current has ledgers, not equilibrium (balanced flows count however small)', () => {
  // A pn diode's flows are tiny against D·c/L of its majority carriers (1e-19 of it at 1 µV), yet
  // real: linear in the bias, and balanced to round-off. Only 0 V is equilibrium.
  const dev = new Device(build({ T: 300, library: [Si], stack: [ohmic(0), layer('Si', 1.5e-6, { donors: units.perCm3(1e17) }), layer('Si', 3e-6, { acceptors: units.perCm3(1e16) }), ohmic(0)], bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }], grid: { hmin: 1e-9, hmax: 2e-8 } }));
  const balance = (V) => {
    dev.set({ contacts: { right: { V } } });
    return check(dev, dev.solve()).items.find((i) => i.name === 'balance');
  };
  assert.match(balance(0).summary, /equilibrium/);
  for (const V of [1e-6, 0.3, -1.5]) {
    const b = balance(V);
    assert.ok(b.ok && /(one ledger sums|2 ledgers sum) to zero/.test(b.summary), `${V} V: ${b.summary}`);
  }
});
