import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';

// Faults found by round 4's outside reading passes (OpenAI's Astra on the steady solve, Sol on
// impedance, the readout and the discretised equations), each confirmed by a probe against a
// second route before it was fixed.

test('eddy mixing in a graded region is the stated law on its segments, in either orientation', () => {
  // A neutral species in a 2 RT grade of μ°, mixed and diffusing on 8 cells, driven by 1 RT across.
  // Per segment the law is N = (D/h)[B(Δ)c_L − B(−Δ)c_R] − (D_mix/h) c̄ Δη, c̄ the log mean of the
  // two ends' c. Leaving the grade out of c̄ moved the contact flux by ±12% between the device and
  // its mirror (−6.649e-3 and −5.297e-3 mol/(m²·s)); the chaser's own Newton solve of the law on
  // this grid gives −5.928665e-3 both ways (Sol's fourth pass).
  const T = 298.15, RT = GAS_CONSTANT * T, L = 1e-6, cells = 8, D = 1e-9, Dm = 1e-8;
  const B = (x) => (Math.abs(x) < 1e-12 ? 1 - x / 2 : x / Math.expm1(x));
  const lm = (a, b) => (a === b ? a : (b - a) / Math.log(b / a));
  for (const mirror of [false, true]) {
    const g = mirror ? [2 * RT, 0] : [0, 2 * RT], mu = mirror ? [RT, 0] : [0, RT];
    const sol = new Device({
      T,
      species: [{ name: 'X', z: 0 }],
      materials: { f: { epsr: 0, species: { X: { D, cRef: 1, mu0: 0 } } } },
      regions: [{ material: 'f', length: L, c0: { X: 1 }, mixing: Dm, grade: { X: { mu0: { x: [0, L], values: g } } } }],
      contacts: {
        left: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: mu[0] } } },
        right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: mu[1] } } },
      },
      grid: { minCells: cells, hmax: L / cells },
    }).solve();
    assert.ok(sol.converged);
    const sgn = mirror ? -1 : 1, J = sgn * sol.contacts.left.flux.X;
    assert.ok(Math.abs(J / -5.928665e-3 - 1) < 1e-6, `${mirror ? 'mirror' : 'original'}: ${J}`);
    // and the law itself, evaluated on the solution's own c (η/RT = ln c + μ°/RT here)
    const x = sol.x, c = sol.c.X, m0 = (xx) => (g[0] + ((g[1] - g[0]) * xx) / L) / RT;
    for (let s = 0; s < x.length - 1; s++) {
      const h = x[s + 1] - x[s], d = m0(x[s + 1]) - m0(x[s]), de = Math.log(c[s + 1] / c[s]) + d;
      const N = (D / h) * (B(d) * c[s] - B(-d) * c[s + 1]) - (Dm / h) * lm(c[s], c[s + 1]) * de;
      assert.ok(Math.abs(N / sol.contacts.left.flux.X - 1) < 1e-6, `segment ${s}: ${N}`);
      assert.ok(Math.abs(sol.flux.X[s] / sol.contacts.left.flux.X - 1) < 1e-9, `readout ${s}`);
    }
  }
});

test("a metal linked into a semiconductor isn't summed as a metal island, which would drop the semiconductor's recombination", () => {
  // e⁻ from a metal contact into a semiconductor whose holes leave the other end, recombining in
  // the bulk and at the face. Driven at the current that −10 mV passes, the device settled at
  // −33.9 mV, converged and silent; its mirror, and a long transient of either, were right
  // (Sol's fourth pass).
  const T = 298.15, RT = GAS_CONSTANT * T;
  const def = (drive, mirror) => {
    const semi = { name: 'semi', material: 's', length: 10e-6, c0: { 'e-': 1, 'h+': 1 } }, met = { name: 'metal', material: 'm', length: 1e-6 };
    const mc = { ...drive, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium', 'h+': 'blocked' } };
    const hc = { V: 0, terminal: 'h+', phi: 'neutral', species: { 'h+': 'equilibrium', 'e-': 'blocked' } };
    const carrier = { D: 1e-9, cRef: 10, mu0: RT * Math.log(10) };
    return {
      T,
      species: [{ name: 'e-', z: -1 }, { name: 'h+', z: 1 }],
      materials: { m: { conductor: { species: 'e-', conductivity: 1e7 } }, s: { epsr: 0, species: { 'e-': carrier, 'h+': carrier } } },
      regions: mirror ? [semi, met] : [met, semi],
      interfaces: [{ phi: 'neutral', species: { 'e-': 'equilibrium' }, reactions: [{ equation: mirror ? 'e-(right) + h+(left) = 0' : 'e-(left) + h+(right) = 0', k0: 1e-5, alpha: 0.5 }] }],
      bulkReactions: [{ equation: 'e- + h+ = 0', kf: { s: 1 } }],
      contacts: mirror ? { left: hc, right: mc } : { left: mc, right: hc },
      grid: { hmin: 1e-6, hmax: 1e-6 },
    };
  };
  for (const mirror of [false, true]) {
    const side = mirror ? 'right' : 'left';
    const ref = new Device(def({ V: -0.01 }, mirror)).solve();
    const sol = new Device(def({ I: ref.terminals[side].current }, mirror)).solve();
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.terminals[side].V + 0.01) < 1e-9, `${side}: ${sol.terminals[side].V}`);
    assert.ok(Math.abs(sol.bulkReactions[0].total / ref.bulkReactions[0].total - 1) < 1e-6);
  }
});

test("a capacitive port charging in a floating neutral end box passes nothing through that end's I = 0 contact", () => {
  // A strictly neutral KCl end, floating with only K⁺ linked, and a gate's point port in its end
  // box stepped 0 → 10 mV. The box's ion storage is the gate's charging there; read after the
  // contact's current, 5.0 A/m² of K⁺ left through the I = 0 contact and the terminals summed to
  // +5. The same end closed (both ions blocked) is the answer (Sol's fourth pass).
  const def = (left) => ({
    T: 298.15,
    species: [{ name: 'K+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { w: { epsr: 0, species: { 'K+': { D: 1e-9, cRef: 100, mu0: 0 }, 'Cl-': { D: 1e-9, cRef: 100, mu0: 0 } } } },
    regions: [{ material: 'w', length: 10e-6, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left, right: { bath: { c: { 'K+': 100, 'Cl-': 100 }, reference: 'Cl-' } } },
    ports: [{ name: 'gate', region: 0, from: 0, to: 0, V: 0, area: 1e5, capacitance: { C: 1 } }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  const run = (left) => {
    const d = new Device(def(left));
    d.solve();
    d.set({ ports: { gate: { V: 0.01 } } });
    return d.step(100e-6);
  };
  const floating = run({ phi: 'neutral', I: 0, terminal: 'K+', species: { 'K+': 'equilibrium', 'Cl-': 'blocked' } });
  const closed = run({ phi: 'neutral', species: { 'K+': 'blocked', 'Cl-': 'blocked' } });
  assert.ok(floating.converged && closed.converged);
  const gate = closed.terminals.gate.current;
  assert.ok(gate > 4.9, `the gate charges: ${gate}`);
  assert.ok(Math.abs(FARADAY * floating.contacts.left.flux['K+']) < 1e-9 * gate, `through the contact: ${floating.contacts.left.flux['K+']}`);
  assert.ok(Math.abs(Object.values(floating.terminals).reduce((a, t) => a + t.current, 0)) < 1e-9 * gate);
  assert.ok(Math.abs(floating.terminals.gate.current / gate - 1) < 1e-9);
  assert.ok(Math.abs(floating.c['K+'][0] / closed.c['K+'][0] - 1) < 1e-9);
});

test("names like Object's own properties are the solution's own keys, and a reaction's ν of them a number", () => {
  // Assigned, a species named '__proto__' set each map's prototype and was missing from it, and a
  // port so named from terminals; ν of 'constructor' read Object's and came out NaN (Sol's
  // fourth pass).
  const own = (o, k) => Object.hasOwn(o, k);
  for (const name of ['__proto__', 'constructor', 'toString']) {
    const sol = new Device({
      species: [{ name: 'X', z: 0 }, { name, z: 0 }],
      materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 }, [name]: { D: 1e-9, cRef: 1, mu0: 0 } } } },
      regions: [{ material: 'f', length: 1e-6, c0: { X: 1, [name]: 1 } }],
      contacts: { left: { phi: 'neutral', species: { X: 'blocked', [name]: 'blocked' } }, right: { phi: 'neutral', species: { X: 'blocked', [name]: 'blocked' } } },
      bulkReactions: [{ equation: `${name} = X`, kf: { f: 1 } }],
      grid: { minCells: 8 },
    }).solve();
    assert.ok(sol.converged);
    for (const map of [sol.c, sol.mu, sol.muStd, sol.V, sol.Vstd, sol.flux, sol.D, sol.contacts.left.flux, sol.contacts.right.flux]) {
      assert.ok(own(map, name), name);
      assert.ok(ArrayBuffer.isView(map[name]) || typeof map[name] === 'number', name);
    }
    assert.ok(own(sol.bulkReactions[0].nu, name));
    assert.equal(sol.bulkReactions[0].nu[name], -1);
    assert.equal(sol.bulkReactions[0].nu.X, 1);
  }
  const F = FARADAY, c0 = 100;
  const sol = new Device({
    T: 300,
    species: [{ name: 'K+', z: 1 }],
    materials: { m: { epsr: 0, species: { 'K+': { D: 1e-9, cRef: c0, mu0: 0 } } } },
    regions: [{ material: 'm', length: 10e-6, c0: { 'K+': c0 }, fixedCharge: -F * c0 }],
    contacts: { left: { V: 0, terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' } },
    ports: [{ name: '__proto__', region: 0, from: 4e-6, to: 6e-6, V: 1e-3, terminal: 'K+', species: { 'K+': { type: 'conductance', G: 1e9 } } }],
    grid: { minCells: 20, hmax: 0.5e-6 },
  }).solve();
  assert.deepEqual(Object.keys(sol.terminals), ['left', 'right', '__proto__']);
  const sum = Object.values(sol.terminals).reduce((a, t) => a + t.current, 0);
  assert.ok(Math.abs(sum) < 1e-12 * Math.abs(sol.terminals.left.current), `${sum}`);
});
