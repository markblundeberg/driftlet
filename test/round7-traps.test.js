// A trap's states under light, and an alloy's knots: what the semi family of random devices
// (bench/combos/families/semi.js) found, each settled by a probe (20 Nov 2026). The old numbers
// are in each test's comment.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units } from '../src/index.js';
import { alloy, build, layer, ohmic, semiconductor, trap } from '../src/kit.js';

// An undoped Si bar 10 µm long between ohmic contacts at 0 V, one trap of 1e15 cm⁻³ and light
// generating uniformly: symmetric, so it carries no current. With `neutral`, strictly neutral (ε =
// 0), an acceptor trap listed X⁻ first and donors filling it, so that at its neutral nodes X⁻ is
// the most abundant charged species: the charge rows' pivot is the species whose row the trap's
// total takes.
const litBar = ({ neutral = false } = {}) => {
  const si = semiconductor('Si');
  if (neutral) si.materials.Si.epsr = 0;
  const t = trap(si, { density: units.perCm3(neutral ? 1e17 : 1e15), tauN: 1e-7, tauP: 1e-7, charge: neutral ? 'acceptor' : 'donor' });
  if (neutral) {
    const sp = t.library.species;
    t.library.species = [sp[0], sp[1], sp[3], sp[2]];
  }
  const L = 10e-6, doping = neutral ? { donors: units.perCm3(1e16) } : {};
  return build({
    T: 300,
    library: [t.library],
    stack: [ohmic(0), layer('Si', L, { name: 'bar', ...doping, c0: t.c0For(doping), grid: { hmin: 5e-9, hmax: L / 60, ratio: 1.1 } }), ohmic(0)],
    bulkReactions: [...t.bulkReactions, { equation: 'photon = e- + h+', fixed: { photon: units.eV(2) }, kf: { Si: 1e3 } }],
  });
};
const totals = (sol) => Array.from(sol.x, (_, g) => sol.species.filter((s) => /^X/.test(s.name)).reduce((a, s) => a + sol.c[s.name][g], 0));

test("a trap's total stays what it was at every node through a long transient under light, and the symmetric bar carries no current", () => {
  // (each capture balance sums its storage against rates that nearly cancel, ~R per volume, and
  // past dt ~ X/(εR) the storage was below their round-off: the balances said only that the rates
  // balance, and the total floated. From a steady state, with nothing changed, the totals drifted
  // 3.5e-7 by 1e6 s and the bar came to carry −2.4e-6 A/m² (2.2e-9 and 3.2e-8 A/m² by 1e8 s in the
  // neutral bar, with its row moved off the pivot); now round-off, ~1e-14)
  for (const neutral of [false, true]) {
    const d = new Device(litBar({ neutral }));
    const S = d.solve();
    assert.ok(S.converged);
    const T0 = totals(S), scale = Math.abs(S.terminals.left.current) + 1e-9;
    for (const t of [1e6, 1e10]) {
      const a = d.advance(t, { tol: 1e-4 });
      assert.ok(a.converged, `advanced to ${t} s`);
      const drift = Math.max(...totals(a).map((v, g) => Math.abs(v / T0[g] - 1)));
      assert.ok(drift < 1e-12, `${neutral ? 'neutral' : 'Si'} bar at ${t} s: the trap's totals drifted ${drift.toExponential(2)}`);
      assert.ok(Math.abs(a.current) < scale, `${neutral ? 'neutral' : 'Si'} bar at ${t} s: it carries ${a.current.toExponential(3)} A/m²`);
    }
  }
});

test("an alloy's fraction profile whose last knot is an ulp short of the layer's end, with a bow, gives a grade new Device takes", () => {
  // (the last piece, an ulp long, was cut in 128 knots, all at the end: refused, "regions[0].grade.e-.mu0.x
  // must be finite and strictly increasing", a field the user never wrote; the semi family's case 205)
  const L = 7.568963116110493e-6;
  const GaAs = semiconductor('GaAs'), InGaAs = semiconductor('GaAs', { material: 'InGaAs' });
  InGaAs.materials.InGaAs.species['h+'] = { ...InGaAs.materials.InGaAs.species['h+'], mu0: GaAs.materials.GaAs.species['h+'].mu0 - units.eV(0.15) };
  for (const end of [L - L * Number.EPSILON, L]) {
    const def = build({
      T: 300,
      library: [GaAs, InGaAs],
      stack: [ohmic(0), alloy('GaAs', 'InGaAs', L, { step: { species: 'e-', value: units.eV(-0.09) }, fraction: { x: [0, L / 2, end], values: [0, 0.3, 1] }, bow: { 'h+': { mu0: units.eV(0.3) } }, donors: units.perCm3(1e15) }), ohmic(0)],
    });
    const { x } = def.regions[0].grade['e-'].mu0;
    assert.ok(x.every((v, k) => k === 0 || v > x[k - 1]), 'knots strictly increasing');
    assert.equal(x[x.length - 1], L);
    assert.ok(new Device(def).solve().converged);
  }
});
