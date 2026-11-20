// Round 4: the Kirchhoff readout (A20) and a failed solve's reason (A26), from Astra's fourth pass.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor } from '../src/kit.js';

// A20: a no-island n-GaAs resistor at 10 mV, its left contact's electron flux inflated after the
// bookkeeping by k times the true current, so the terminals miss Kirchhoff by k·I: warned at every k.
test('A20: a Kirchhoff miss beside another terminal\'s real current is not round-off', () => {
  const def = build({ T: 300, library: [semiconductor('GaAs')], stack: [ohmic(0), layer('GaAs', 1e-6, { donors: units.perCm3(1e17) }), ohmic(0.01)] });
  const d = new Device(def), s0 = d.solve(), sv = d.solver, I = s0.terminals.left.current, ie = def.species.findIndex((x) => x.name === 'e-');
  assert.equal(sv.islands.length, 0);
  assert.equal(s0.warnings.length, 0);
  const orig = sv._assembleBookkeeping.bind(sv);
  for (const k of [0.1, 1.01, 3, 1e6]) {
    sv._assembleBookkeeping = (dt) => { orig(dt); sv.contactFlux.left[ie] += (k * I) / (-1 * FARADAY); };
    assert.ok(d.solution().warnings.some((w) => w.includes("don't add up")), `k ${k}`);
  }
  sv._assembleBookkeeping = orig;
});

// A26: a solve failing at a singular border says Newton's reason.
test('A26: a failed solve\'s warning carries Newton\'s reason', () => {
  const d = new Device(build({ library: [semiconductor('Si')], stack: [ohmic(0), layer('Si', 1e-6, { donors: 1 }), ohmic(0.01)] }));
  d.solve();
  const f = d.solver._factor, why = 'the bordered system is singular: a conserved amount or a floating terminal is not determined';
  d.solver._factor = () => { throw new Error(why); };
  try {
    const s = d.solve({ maxSteps: 1, continuation: false });
    assert.ok(s.warnings.some((w) => w.includes(`(Newton's last try: ${why}), and nothing more specific showed why`)), s.warnings.join(' | '));
  } finally {
    d.solver._factor = f;
  }
});
