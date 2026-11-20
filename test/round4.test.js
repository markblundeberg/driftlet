import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT, SolverError, units } from '../src/index.js';
import { complexSolve } from '../src/impedance.js';
import { aqueous, bath, build, half, layer, metal, ohmic, semiconductor } from '../src/kit.js';
import { addTerm } from '../src/sum.js';

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


// Faults found in the discretised equations by round 4's outside reading passes (a second Sol over
// the kernels), each confirmed by a probe against a second route before it was fixed.

const T = 298.15, RT = GAS_CONSTANT * T;

// B6: an affinity summed term by term lost its digits where the participants' μ°/RT are large.
test('addTerm: an affinity summed exactly, in any order (round 4, B6)', () => {
  // Exact sums of doubles as integers in units of 2^-200 (every term here is a multiple of that).
  const big = (x) => BigInt(x * 2 ** 200);
  // large terms that nearly cancel, as A's and B's μ°/RT do, the last leaving a ~ 1e-9
  const terms = [[1, 40.123456789], [-1, 39.98765432101], [3, 0.30000000000000004], [-2, 0.1], [0.5, 7.77]];
  terms.push([-1, terms.reduce((s, [nu, x]) => s + nu * x, 0) - 1e-9]);
  const want = terms.reduce((s, [nu, x]) => s + big(nu * 2) * big(x), 0n); // ν·2 keeps ½ an integer
  const scale = 2 ** 200, wantA = Number(want) / (2 * scale * scale);
  let plainSpread = 0;
  for (let k = 0; k < 720; k++) {
    // every order (6! of them)
    const order = [], rest = [...terms];
    for (let m = k, j = rest.length; j > 0; j--) order.push(rest.splice(m % j, 1)[0]), (m = Math.floor(m / j));
    const acc = new Float64Array(2);
    let plain = 0;
    for (const [nu, x] of order) addTerm(acc, nu, x, 0), (plain += nu * x);
    const a = acc[0] + acc[1];
    assert.ok(Math.abs(a - wantA) <= Math.abs(wantA) * Number.EPSILON + 1e-28, `order ${k}: ${a} vs ${wantA}`);
    plainSpread = Math.max(plainSpread, Math.abs(plain - wantA) / Math.abs(wantA));
  }
  assert.ok(plainSpread > 1e-8, `a plain sum loses a's digits here (${plainSpread})`); // the test's teeth
  // an ∞ passes through (a participant at c = 0)
  const acc = new Float64Array([2, 0]);
  addTerm(acc, -1, Infinity, 0);
  addTerm(acc, 3, 0.1, 0);
  assert.equal(acc[0] + acc[1], -Infinity);
});

// A + B + E = C + D in a 10 µm film, μ°/RT of 40, −40, 0.3, 0.3, 0: A, B, E, C held at the left at
// their μ°, D drawn off at the right half a thermal unit below. The same physics in any order of the
// equation; fast, it's local equilibrium (K = 1), solved by its reaction invariants: c_A = c_B = c_E
// = a, c_C = 2 − a, c_D = s − a, s linear from 2, (2 − a)(s − a) = a³, c_D = e^{−1/2} at the right.
test('a fast bulk reaction among large μ° comes out the same in any order of its equation (round 4, B6)', () => {
  const mu0 = { A: 40, B: -40, E: 0.3, C: 0.3, D: 0 };
  const film = (equation, kf) => new Device({
    T, species: Object.keys(mu0).map((name) => ({ name, z: 0 })),
    materials: { f: { epsr: 4, species: Object.fromEntries(Object.entries(mu0).map(([k, m]) => [k, { D: 1e-9, cRef: 1, mu0: m * RT }])) } },
    regions: [{ name: 'a', material: 'f', length: 1e-5, c0: { A: 1, B: 1, E: 1, C: 1, D: 1 } }],
    bulkReactions: [{ equation, kf: { f: kf } }],
    contacts: {
      left: { phi: 'neutral', species: { A: { type: 'equilibrium', mu: mu0.A * RT }, B: { type: 'equilibrium', mu: mu0.B * RT }, E: { type: 'equilibrium', mu: mu0.E * RT }, C: { type: 'equilibrium', mu: mu0.C * RT }, D: 'blocked' } },
      right: { phi: 'neutral', species: { A: 'blocked', B: 'blocked', E: 'blocked', C: 'blocked', D: { type: 'equilibrium', mu: (mu0.D - 0.5) * RT } } },
    },
    grid: { minCells: 32 },
  }).solve();
  const cD = Math.exp(-0.5);
  let lo = 0, hi = 2;
  for (let k = 0; k < 200; k++) {
    const a = (lo + hi) / 2;
    if (a * a * a - (2 - a) * cD > 0) hi = a;
    else lo = a;
  }
  const exact = (1e-9 * (2 - ((lo + hi) / 2 + cD))) / 1e-5;
  const N = ['A + B + E = C + D', 'A + E + B = C + D', 'E + A + B = C + D'].map((eq) => {
    const s = film(eq, 1e14);
    assert.ok(s.converged, eq);
    return s.contacts.right.flux.D;
  });
  // (summed term by term, these spread by 4.9e-4)
  assert.ok(Math.max(...N) / Math.min(...N) - 1 < 1e-12, `the orders' fluxes ${N}`);
  for (const n of N) assert.ok(Math.abs(n / exact - 1) < 1e-9, `${n} vs local equilibrium's ${exact}`);
});

// The same at a face: A + B + E on the left = C + D on the right, blocked across it otherwise.
test('a face reaction among large μ° comes out the same in any order of its equation (round 4, B6)', () => {
  const mu0 = { A: 40, B: -40, E: 0.3, C: 0.3, D: 0 };
  const cell = (equation) => new Device({
    T, species: Object.keys(mu0).map((name) => ({ name, z: 0 })),
    materials: { f: { epsr: 4, species: Object.fromEntries(Object.entries(mu0).map(([k, m]) => [k, { D: 1e-9, cRef: 1, mu0: m * RT }])) } },
    regions: ['a', 'b'].map((name) => ({ name, material: 'f', length: 1e-5, c0: { A: 1, B: 1, E: 1, C: 1, D: 1 } })),
    interfaces: [{ phi: 'neutral', species: { A: 'blocked', B: 'blocked', E: 'blocked', C: 'blocked', D: 'blocked' }, reactions: [{ equation, k0: 1e14 }] }],
    contacts: {
      left: { phi: 'neutral', species: { A: { type: 'equilibrium', mu: mu0.A * RT }, B: { type: 'equilibrium', mu: mu0.B * RT }, E: { type: 'equilibrium', mu: mu0.E * RT }, C: 'blocked', D: 'blocked' } },
      right: { phi: 'neutral', species: { A: 'blocked', B: 'blocked', E: 'blocked', C: { type: 'equilibrium', mu: mu0.C * RT }, D: { type: 'equilibrium', mu: (mu0.D - 0.5) * RT } } },
    },
    grid: { minCells: 32 },
  }).solve();
  const N = ['A(left) + B(left) + E(left) = C(right) + D(right)', 'A(left) + E(left) + B(left) = C(right) + D(right)', 'E(left) + A(left) + B(left) = C(right) + D(right)'].map((eq) => {
    const s = cell(eq);
    assert.ok(s.converged, eq);
    return s.contacts.right.flux.D;
  });
  // (summed term by term, these spread by 6.7e-15)
  assert.ok(Math.max(...N) / Math.min(...N) - 1 < 1e-15, `the orders' fluxes ${N}`);
});

// B13: the dominant adsorbate's storage diagonal, (Γθ_q/dt)(1 − θ_q), was 0 or noise where the bare
// fraction is below an ulp of 1; test/coverage.test.js's cell, its film 100 kJ/mol more stable.
test('a surface filling to a bare fraction of 4e-19 completes its transient, onto the isotherm (round 4, B13)', () => {
  const G = 1e-5, L = 10e-6, area = 1e4, mu0 = -100e3, V = 0.05;
  const lib = aqueous(['Na+', 'OH-', 'Fe2+'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const dev = new Device(build({
    T, library: [lib],
    stack: [{}, layer('water', L, { name: 'film' }), bath({ 'Na+': 100, 'OH-': 102, 'Fe2+': 1 }, 'OH-')],
    ports: [{ name: 'metal', region: 'film', V, terminal: 'e-', area, surface: { OHads: { mu0, capacity: G, theta0: 1e-4 } }, reactions: [{ equation: 'OH- = OHads + e-', k0: 1e-4, alpha: 0.5, bare: true }] }],
    grid: { hmin: L / 20, hmax: L / 20 },
  }));
  const s = dev.advance(1e-4, { tol: 1e-6 }); // (it stalled at 7.4e-9 s, its steps cut to 3e-19 s)
  assert.ok(s.converged, s.warnings.join('; '));
  const p = s.ports[0], g = s.x.indexOf(p.x[5]);
  const want = 1 / (1 + Math.exp((s.mu['OH-'][g] + FARADAY * V - mu0) / RT));
  assert.ok(want < 1e-18 && Math.abs(p.bare[5] / want - 1) < 1e-6, `bare ${p.bare[5]} vs the isotherm's ${want}`);
});

// B5: a velocity contact's flux v c (e^x − 1) took e^x − 1 as a difference.
test("a velocity contact's flux is exact for a small drive (round 4, B5)", () => {
  // A neutral tracer film, the velocity link on the left (outside μ = 0), equilibrium on the right at
  // μ = −δ RT: linear diffusion, so N = (1 − e^{−δ})/(1/v + L/D) exactly.
  const L = 1e-6, D = 1e-9, v = D / L, delta = 1e-8;
  const s = new Device({
    T, species: [{ name: 'X', z: 0 }], materials: { f: { epsr: 4, species: { X: { D, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: L, c0: { X: 1 } }],
    contacts: { left: { phi: 'neutral', species: { X: { type: 'velocity', v, mu: 0 } } }, right: { phi: 'neutral', species: { X: { type: 'equilibrium', mu: -delta * RT } } } },
    grid: { minCells: 8 },
  }).solve();
  const exact = -Math.expm1(-delta) / (1 / v + L / D);
  assert.ok(s.converged);
  assert.ok(Math.abs(s.contacts.left.flux.X / exact - 1) < 1e-12, `N ${s.contacts.left.flux.X} vs ${exact}`); // (4.3e-9 before)
});


// Faults found in the impedance and the readout by an outside reading pass (OpenAI's Sol, its
// fourth), each confirmed by a probe against a second route before it was fixed.

// A film with a contact passing K⁺ on the left, a membrane port `mem` passing K⁺ in its middle and
// a Cl⁻ reference port, both current-driven terminals sharing the K⁺ that the law conserves.
const sharedK = (left, mem) => {
  const film = { 'Na+': 100, 'K+': 10, 'Cl-': 110 }, K = (G) => ({ 'K+': { type: 'conductance', G, offset: 0 } });
  return build({
    T: 298.15,
    library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })],
    stack: [{ ...left, phi: 'neutral', species: K(1e6) }, layer('water', 1e-6, { name: 'film', c0: film }), { V: 0, phi: { type: 'capacitive', C: 0.5 }, zeroCharge: -0.05 }],
    ports: [{ name: 'mem', region: 'film', from: 0.4e-6, to: 0.6e-6, ...mem, species: K(1e13) }, { name: 'ref', region: 'film', from: 0.7e-6, to: 0.8e-6, V: 0, species: { 'Cl-': { type: 'equilibrium', offset: 0 } } }],
    grid: { hmin: 5e-8, hmax: 5e-8 },
  });
};
// A neutral film between a pinned contact and a capacitive plate, the plate driven as given.
const plate = (right) => new Device({
  species: [{ name: 'X', z: 0 }],
  materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
  contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } }, right: { ...right, phi: { type: 'capacitive', C: 1e-3 }, zeroCharge: 0, species: { X: 'blocked' } } },
  grid: { minCells: 8 },
});
const near = (a, b, q, tol) => Math.hypot(a.Z.re[q] - b.Z.re[q], a.Z.im[q] - b.Z.im[q]) / Math.hypot(b.Z.re[q], b.Z.im[q]) < tol;

test('a current-driven terminal measured where another shares its conserved law: the held twin\'s Z', () => {
  // (The left contact's row read by the law's change, iω F δS/λ = δI_left, held its δI to 0
  // where mem's own small signal fed the same K⁺: the left AC-shorted, mem's Im Z 37× too small
  // at 1 kHz and 3.6e7× at 1 Hz, relError 1e-13. A sine stepped through step() agrees with the
  // held twin.)
  const fs = [1e3, 1];
  const d = new Device(sharedK({ I: 0 }, { I: 0 })), V = d.solve().terminals.mem.V;
  const z = d.impedance(fs, { terminal: 'mem' });
  const twin = new Device(sharedK({ I: 0 }, { V })).impedance(fs, { terminal: 'mem' });
  fs.forEach((f, q) => assert.ok(near(z, twin, q, 1e-5), `${f} Hz: ${z.Z.re[q]} ${z.Z.im[q]}i, held twin ${twin.Z.re[q]} ${twin.Z.im[q]}i`));
  // (the left contact, its law its own, as before: its held twin's)
  const zl = d.impedance(fs, { terminal: 'left' });
  const twinL = new Device(sharedK({ V: d.solution().terminals.left.V }, { I: 0 })).impedance(fs, { terminal: 'left' });
  fs.forEach((f, q) => assert.ok(near(zl, twinL, q, 1e-5), `left, ${f} Hz: ${zl.Z.re[q]} ${zl.Z.im[q]}i, held twin ${twinL.Z.re[q]} ${twinL.Z.im[q]}i`));
});

test('a solution read after a set() of a floating terminal\'s drive alone: what it passed as solved, until it moves', () => {
  // (Read at the new drive, a capacitive plate set from I: 0 to 2 A/m² read 2 A/m² flowing at no
  // step, the plate uncharged and the terminals summing to −2, without a word; behind 10 Ω·m², a
  // source set from 0 to 0.5 V read 0.05 A/m². Held terminals were read as solved already.)
  const sum = (s) => Object.values(s.terminals).reduce((t, x) => t + x.current, 0);
  for (const [from, to] of [[{ I: 0 }, { I: 2 }], [{ V: 0, R: 10 }, { V: 0.5 }]]) {
    const d = plate(from);
    d.solve();
    d.set({ contacts: { right: to } });
    const s = d.solution();
    assert.ok(s.terminals.right.current === 0, `${JSON.stringify(to)}: ${s.terminals.right.current} A/m²`);
    assert.ok(sum(s) === 0, `${sum(s)}`);
    const t = d.step(1e-9); // (and the step from there passes the new drive)
    assert.ok(Math.abs(t.plates.right.charge - (to.I ?? 0.05) * 1e-9) < 1e-3 * (to.I ?? 0.05) * 1e-9, `${t.plates.right.charge}`);
  }
  // (a port's current likewise)
  const d = new Device(sharedK({ I: 0 }, { I: 0 }));
  d.solve();
  d.set({ ports: { mem: { I: 1e-3 } } });
  assert.ok(d.solution().terminals.mem.current === 0);
  assert.equal(d.step(1e-9).terminals.mem.current, 1e-3);
});

test('a refused impedance() leaves a pending step response as it was', () => {
  // (Each solved the steady state at the new drive before refusing: the next step read 0 where
  // the plate's charging current was 3420 A/m².)
  for (const bad of [[[0], {}], [[1], { terminal: 'nope' }], [[NaN], {}], [[1, -1], {}]]) {
    const d = plate({ V: 0 }), twin = plate({ V: 0 });
    for (const x of [d, twin]) {
      x.solve();
      x.set({ contacts: { right: { V: 0.1 } } });
    }
    assert.throws(() => d.impedance(...bad), /impedance: (frequencies|no terminal)/);
    const a = d.step(1e-9).terminals.right.current, b = twin.step(1e-9).terminals.right.current;
    assert.ok(b > 3e3 && a === b, `${JSON.stringify(bad)}: ${a} A/m², its twin ${b}`);
  }
});

test('a terminal no small-signal current returns through is refused, not read as NaN', () => {
  // (A dielectric up to a neutral contact passes nothing: Z came out NaN, with a warning that it
  // was round-off.)
  const d = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } }, right: {} },
    grid: { minCells: 8 },
  });
  for (const terminal of ['left', 'right']) assert.throws(() => d.impedance([1, 1e6], { terminal }), new RegExp(`no current returns through '${terminal}' at 1 Hz`));
});

test('a plate\'s impedance at an absurd frequency, and a singular small system refused', () => {
  // (Divided by c² + d², which overflowed: at 1e160 Hz Z read 0, and behind 0.4 Ω·m², NaN; a zero
  // pivot gave NaN for every floating terminal's δV.)
  const both = (R) => new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { V: 0, ...(R ? { R } : {}), phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } }, right: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } } },
    grid: { minCells: 8 },
  });
  const C = (4 * EPS0) / 1e-6;
  for (const f of [1e160, 1e300]) {
    for (const R of [0, 0.4]) {
      const z = both(R).impedance([f]), im = -1 / (2 * Math.PI * f * C);
      assert.ok(Math.abs(z.Z.re[0] - R) < 1e-12 && Math.abs(z.Z.im[0] / im - 1) < 1e-12, `${f} Hz, R ${R}: ${z.Z.re[0]} ${z.Z.im[0]}i, expected ${R} ${im}i`);
    }
  }
  const zero = () => [[Float64Array.from([0, 0, 1]), Float64Array.from([0, 0, 1])], [new Float64Array(3), new Float64Array(3)]];
  assert.throws(() => complexSolve(...zero(), 2, 1), (e) => e instanceof SolverError && /singular at 1 Hz/.test(e.message));
  const rank1 = [[Float64Array.from([1, 2, 1]), Float64Array.from([2, 4, 3])], [new Float64Array(3), new Float64Array(3)]];
  assert.throws(() => complexSolve(...rank1, 2, 1), SolverError);
});

test('a floating terminal\'s δV carries its readings\' error: a gate behind R at 1 nHz, a Pt island behind a bath\'s R', () => {
  // (Read from its laws and cuts with no error of its own, δV went into Kirchhoff's readout as
  // exact: a gate behind 1e4 Ω·m², measured at the back contact at 1 nHz, put Re Z − R at −5.4
  // for the held gate's 2869, relError 1e-9; a Pt island behind 1 Ω·m² read Z − R 0.3% off at
  // 3e-5–1e-4 Hz, relError 5e-12, its bath's cuts 2000× past their estimates. Now the floating
  // responses are split into their regions' shifts and the rest where a reading is unresolved.)
  const si = semiconductor('Si');
  si.materials.Si.species['e-'].mu0 -= FARADAY * 0.3;
  si.materials.Si.species['h+'].mu0 += FARADAY * 0.3;
  const mos = (gate) => build({
    T: 300, library: [si, metal('Au')], materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [ohmic(0), layer('Si', units.um(15), { acceptors: units.perCm3(1e16) }), { dipole: -0.3 }, layer('SiO2', units.nm(20), { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }), { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.7 }, layer('Au', units.nm(50)), ohmic(gate, ['e-'])],
    bulkReactions: [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1e-6, tauP: 1e-6 } } }],
    grid: { hmin: units.nm(0.1), hmax: units.nm(10), ratio: 1.1 },
  });
  const fs = [1e-9, 1e-8], ref = new Device(mos(0.6)).impedance(fs, { terminal: 'right' });
  for (const R of [1e3, 1e4]) {
    const z = new Device(mos({ V: 0.6, R })).impedance(fs, { terminal: 'left' });
    fs.forEach((f, q) => {
      const e = Math.abs(z.Z.re[q] - R - ref.Z.re[q]), said = z.relError[q] * Math.hypot(z.Z.re[q], z.Z.im[q]);
      assert.ok(e < 1e-3 * ref.Z.re[q] && e < said * 10 + 1e-3, `R ${R}, ${f} Hz: Re Z − R ${z.Z.re[q] - R} for ${ref.Z.re[q]}, its error said ${said}`);
    });
  }
  const c = { 'K+': 100, 'Cl-': 125, 'Fe3+': 5, 'Fe2+': 5 };
  const face = { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1e-3, alpha: 0.5 }], phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 };
  const island = (left) => build({ library: [aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+']), metal('Pt')], stack: [bath(c, 'Cl-', left), layer('water', 5e-6, { c0: c }), face, layer('Pt', 1e-6), face, layer('water', 5e-6, { c0: c }), bath(c, 'Cl-', 0)], grid: { hmin: 0.3e-9, hmax: 1e-6, ratio: 1.3 } });
  const gs = [3e-5, 1e-4], Z0 = new Device(island({ V: 0 })).impedance(gs);
  for (const R of [1, 1e3]) {
    const z = new Device(island({ V: 0, R })).impedance(gs);
    gs.forEach((f, q) => assert.ok(Math.hypot(z.Z.re[q] - R - Z0.Z.re[q], z.Z.im[q] - Z0.Z.im[q]) < 1e-5 * Math.hypot(Z0.Z.re[q], Z0.Z.im[q]), `R ${R}, ${f} Hz: Z − R ${z.Z.re[q] - R} ${z.Z.im[q]}i for ${Z0.Z.re[q]} ${Z0.Z.im[q]}i`));
    assert.deepEqual(z.warnings.filter((w) => /impedance:/.test(w)), []);
  }
});

// A Pt electrode in a Fe³⁺/Fe²⁺ bath at its open circuit: at 10 nHz its response is nearly a
// uniform shift, which the readout splits off (Sol's fourth pass, S5).
const redox = (right) => {
  const c = { 'Fe3+': 9.6, 'Fe2+': 0.536, 'Cl-': 29.872 };
  return build({
    T: 300, library: [aqueous(['Fe3+', 'Fe2+', 'Cl-'], { epsr: 78.3 }), metal('Pt')],
    stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.395 }, zeroCharge: -0.119, reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 0.0411, alpha: 0.377 }] }, layer('water', 9.83e-6, { c0: { ...c } }), bath({ ...c }, 'Cl-', right)],
    grid: { hmin: 2.6e-10, hmax: 3.7e-7, ratio: 1.136 },
  });
};

test('the profiles are those of the reading Z kept', () => {
  // (The split's solve overwrote them even where Z set it aside. Here its shift is made
  // nonsense, so that it is: the profiles moved with the nonsense, 0.1% at a node.)
  const run = (scale) => {
    const d = new Device(redox({ I: 0 }));
    d.solver._regionShift = (yr, yi) => [yr.map((v) => v * scale), yi.map((v) => v * scale)];
    return d.impedance([1e-8], { profiles: true });
  };
  const a = run(1e6), b = run(1e9);
  assert.ok(a.relError[0] > 1e-2, `${a.relError[0]}`); // (the split set aside: Z is the plain reading's, warned)
  assert.deepEqual(Array.from(a.profiles[0].phi.re), Array.from(b.profiles[0].phi.re));
  assert.deepEqual(Array.from(a.profiles[0].mu['Fe3+'].im), Array.from(b.profiles[0].mu['Fe3+'].im));
});

test('a GMRES left above its tolerance, where relError doesn\'t show it, is warned of', () => {
  // (No device in the tests leaves it so; the solver's report of where it did is passed on.)
  const d = new Device(redox({ I: 0 })), was = d.solver.impedance.bind(d.solver);
  d.solver.impedance = (fs, opts) => ({ ...was(fs, opts), unconverged: [1e-3] });
  const z = d.impedance([1e-3, 1]);
  assert.ok(!('unconverged' in z));
  assert.equal(z.warnings.filter((w) => /stopped short of its tolerance at 0\.001 Hz \(1 of 2/.test(w)).length, 1);
  assert.equal(new Device(redox({ I: 0 })).impedance([1e-3, 1]).warnings.filter((w) => /stopped short/.test(w)).length, 0);
});

test('a port that passes no charge, and eddy mixing, leave the current cuts as they were', () => {
  // (Each cut is the total current only where nothing takes current in or out: a neutral
  // tracer's exchange counted as a window anyway, and a mixing region's cuts were left out,
  // though mixing carries none. With mixing: 1e-30 every cut was gone, relError 50 at 10 nHz.)
  const withPort = (d, from, to) => {
    d.species.push({ name: 'N', z: 0, cRef: 1000 });
    d.materials.water.species.N = { D: 1e-9, mu0: 0 };
    const w = d.regions.find((r) => r.material === 'water');
    w.c0 = { ...w.c0, N: 1 };
    d.ports = [{ name: 'tracer', region: d.regions.indexOf(w), from, to, V: 0, species: { N: { type: 'exchange', k: 1e-3, mu: 0 } } }];
    return d;
  };
  const withMixing = (d) => {
    d.regions.find((r) => r.material === 'water').mixing = 1e-30;
    return d;
  };
  const fs = [1e-8, 1e-4];
  // (held at the voltage it has at no current, measured at either end)
  const held = (make) => new Device(make(redox({ V: new Device(make(redox({ I: 0 }))).solve().terminals.right.V })));
  for (const terminal of ['left', 'right']) {
    const plain = held((d) => d).impedance(fs, { terminal });
    for (const [tag, make] of [['tracer by the bath', (d) => withPort(d, 9.0e-6, 9.5e-6)], ['tracer by the Pt', (d) => withPort(d, 0.1e-6, 0.5e-6)], ['mixing', withMixing]]) {
      const z = held(make).impedance(fs, { terminal });
      fs.forEach((f, q) => {
        assert.ok(z.relError[q] < 3 * plain.relError[q], `${tag}, ${terminal}, ${f} Hz: relError ${z.relError[q]}, plain ${plain.relError[q]}`);
        assert.ok(Math.hypot(z.Z.re[q] - plain.Z.re[q], z.Z.im[q] - plain.Z.im[q]) < 1e-6 * Math.hypot(plain.Z.re[q], plain.Z.im[q]), `${tag}, ${terminal}, ${f} Hz`);
      });
    }
  }
});
