// The steady solve, against the pointers an outside re-read of it observed (Astra, round 8, A1–A11;
// settled with probes on f1f3094). Each test's comment gives the old numbers; the right ones come
// from a closed form, a transient run long enough to settle, or a steady state's own conditions
// (no flux through a blocked end, equilibrium where no reservoir reaches).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous, bath, build, check, layer } from '../src/kit.js';

const RT = GAS_CONSTANT * 300;
const entry = (D, cRef = 1) => ({ D, cRef, mu0: 0 });
const contact = () => ({ V: 0, terminal: 'e-', phi: 'bulk', species: { 'e-': 'equilibrium' } });
const tracer = (d, name, ...levels) => {
  d.species.push({ name, z: 0 });
  levels.forEach((level, j) => {
    d.materials[d.regions[j].material].species[name] = entry(1e-9);
    d.regions[j].c0[name] = level;
  });
};
// Y = Z in a neutral film; with a second region where Y can't move, the steady solve goes by giant
// steps (the partly immobile Y of test/round4-newton.test.js).
function film(oneRegion = false) {
  const def = {
    T: 300,
    species: [{ name: 'e-', z: -1 }, { name: 'Y', z: 0 }, { name: 'Z', z: 0 }],
    materials: {
      m1: { epsr: 0, species: { 'e-': entry(1e-4, 1000), Y: entry(1e-9, 1000), Z: entry(1e-9, 1000) } },
      m2: { epsr: 0, species: { 'e-': entry(1e-4, 1000), Y: entry(0, 1000) } },
    },
    regions: [{ name: 'A', material: 'm1', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 1, Z: 1 } }],
    contacts: { left: contact(), right: contact() },
    grid: { minCells: 8 },
  };
  if (!oneRegion) {
    def.regions.push({ name: 'B', material: 'm2', length: 1e-6, fixedCharge: FARADAY, c0: { Y: 1 } });
    def.bulkReactions = [{ equation: 'Y = Z', kf: { m1: 1 } }];
  }
  return def;
}
const named = (s, where) => s.warnings.some((w) => w.startsWith(`${where}:`));

// Y and Z held A RT apart at the left contact, blocked at the right, Y = Z at kf.
const heldApart = (A, kf, cells = 8) => {
  const d = film(true);
  d.regions[0].c0 = { Y: 1000, Z: 1000 };
  d.contacts.left.species.Y = { type: 'equilibrium', mu: A * RT };
  d.contacts.left.species.Z = { type: 'equilibrium', mu: 0 };
  d.bulkReactions = [{ equation: 'Y = Z', kf: { m1: kf } }];
  d.grid = { minCells: cells };
  return d;
};

test('a reacting group at one contact is held flat only where its affinity there is round-off, not under 1e-9 RT', () => {
  // (Held flat below 1e-9 RT, flat rows replaced every balance, so each box reacted at kf c A with
  // no flux between them: at kf 1e10, A = 5e-10 RT, the total was 5.0e-3 mol/(m²·s) against the
  // contact's 3.1e-4 (16×), and Y left through the blocked right end at 3.1e-4. At kf 1e3 on 64
  // cells, 5.0e-10 against the continuum's 3.1409e-10, the contact's intake 3.9e-12. Converged,
  // unwarned, check() passing it.)
  for (const A of [5e-10, 1e-12, 1e-15]) {
    const s = new Device(heldApart(A, 1e10)).solve();
    assert.ok(s.converged, s.warnings.join('; '));
    const total = s.bulkReactions[0].total, inY = s.contacts.left.flux.Y;
    assert.ok(Math.abs(total / inY - 1) < 1e-9, `A ${A}: total ${total}, the contact passes ${inY}`);
    assert.ok(Math.abs(s.contacts.right.flux.Y) < 1e-15 * inY, `A ${A}: ${s.contacts.right.flux.Y} through the blocked end`);
  }
  // Resolved, against the continuum: Y − Z decays from d₀ = c(e^A − 1) over ℓ = √(D/2kf), and the
  // total is kf d₀ ℓ tanh(L/ℓ), what the contact passes.
  for (const [kf, cells, tol] of [[1e3, 64, 3e-4], [1, 8, 1e-5]]) {
    const A = 5e-10, s = new Device(heldApart(A, kf, cells)).solve(), ell = Math.sqrt(1e-9 / (2 * kf));
    const want = kf * 1000 * Math.expm1(A) * ell * Math.tanh(1e-6 / ell);
    assert.ok(Math.abs(s.bulkReactions[0].total / want - 1) < tol, `kf ${kf}: total ${s.bulkReactions[0].total}, the continuum's ${want}`);
    assert.ok(Math.abs(s.contacts.left.flux.Y / want - 1) < tol, `kf ${kf}: the contact passes ${s.contacts.left.flux.Y}, the continuum's ${want}`);
  }
  // At an affinity that is zero, still held flat (a MOS capacitor's back contact needs it).
  assert.equal(new Device(heldApart(0, 1e10)).solver.flatStretches.length, 2);
});

// Two electrode ports in one film, each running Fe²⁺ = Fe³⁺ + e⁻, the Fe³⁺ closed between them.
const twoPorts = () => {
  const lib = aqueous(['Na+', 'Cl-', 'Fe2+', 'Fe3+'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  return build({
    T: 298.15,
    library: [lib],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'Cl-': 105, 'Fe2+': 1, 'Fe3+': 1 } }), bath({ 'Na+': 100, 'Cl-': 102, 'Fe2+': 1 }, 'Cl-')],
    ports: ['p1', 'p2'].map((name, j) => ({
      name, region: 'film', from: 6e-6 * j, to: 4e-6 + 6e-6 * j, I: 0, terminal: 'e-', area: 1e4,
      reactions: [{ equation: 'Fe2+ = Fe3+ + e-', k0: 1e-3, alpha: 0.5 }],
    })),
    grid: { hmin: 5e-7, hmax: 5e-7 },
  });
};

test("two electrode ports' opposed currents through one closed couple solve steady; a net feed is still refused", () => {
  // (Refused at zero steps, "what it feeds (Fe3+ in film) is closed but for it": the precheck
  // took any law with an electrode's charge column as charging, its net feed uncounted. The
  // transient settles at −0.64184280256907, −0.65531790550188 V.)
  const at = (I1, I2) => {
    const d = new Device(twoPorts());
    d.solve();
    return d.set({ ports: { p1: { I: I1 }, p2: { I: I2 } } });
  };
  const s = at(1e-3, -1e-3).solve();
  assert.ok(s.converged, s.warnings.join('; '));
  const t = at(1e-3, -1e-3);
  assert.ok(t.advance(1e6, { maxSteps: 20000 }).converged);
  const late = t.solution();
  for (const p of ['p1', 'p2']) assert.ok(Math.abs(s.terminals[p].V - late.terminals[p].V) < 1e-12, `${p}: ${s.terminals[p].V} V steady, ${late.terminals[p].V} V after 1e6 s`);
  // +1 in, −0.9 out: the couple charges for ever.
  const net = at(1e-3, -0.9e-3).solve();
  assert.equal(net.converged, false);
  assert.ok(net.warnings.some((w) => w.includes('no steady state') && w.includes("closed but for it and p2, whose currents into it don't balance")), net.warnings.join('; '));
});

test("a node-local inventory is kept through a cold retry and a nested held solve (the stages don't re-capture it)", () => {
  // (Q immobile, filled to 5 from a contact then blocked, a step taken; the first two stage
  // attempts refused. The cold retry laid out c0 = 1, the resistance's continuation re-captured
  // Q's references from it, and the solve returned Q = 1 everywhere, converged, unwarned.)
  const d = film(true);
  tracer(d, 'Q', 1);
  d.contacts.left.species.Q = { type: 'equilibrium', mu: RT * Math.log(5) };
  d.contacts.right.R = 0.1;
  const dev = new Device(d);
  dev.solve();
  dev.set({ materials: { m1: { species: { Q: { D: 0 } } } }, contacts: { left: { species: { Q: 'blocked' } } } });
  const start = dev.step(1e-3);
  assert.ok(start.c.Q.every((c) => Math.abs(c - 5) < 1e-9));
  const v = dev.solver, steps = v._steadySteps;
  let calls = 0;
  v._steadySteps = function (opts) {
    return ++calls <= 2 ? { converged: false, steps: 1, iterations: 1, history: [] } : steps.call(this, opts);
  };
  const s = dev.solve();
  assert.ok(calls > 2 && s.converged, s.warnings.join('; '));
  assert.ok(s.c.Q.every((c) => Math.abs(c - 5) < 1e-9), `Q: ${Array.from(s.c.Q).join(', ')}`);
});

test("two exchange contacts pulling one way and the other are balanced by their laws, solved directly or by giant steps", () => {
  // (W between reservoirs at 0.5 and 8 mol/m³, its level balanced at 2 where k is too small to
  // resolve. Directly at k = 1e-40, starts of 1 and 3 came back as 1 and 3: rateHeld read each
  // contact's intake from its end box, the transport's residual, not the exchange law. By giant
  // steps at k = 1e-32, the same, the two through-fluxes 3× apart; links pushing both ways
  // passed the one-sign test.)
  for (const [direct, k] of [[true, 1e-40], [false, 1e-32]]) {
    for (const start of [1, 3]) {
      const d = film(direct);
      tracer(d, 'W', ...d.regions.map(() => start));
      for (const [side, c] of [['left', 0.5], ['right', 8]]) d.contacts[side].species.W = { type: 'exchange', k, mu: RT * Math.log(c) };
      const s = new Device(d).solve();
      const silent = s.converged && s.warnings.length === 0;
      assert.ok(!silent || Math.abs(s.c.W[0] - 2) < 1e-6, `direct ${direct}, from ${start}: W at ${s.c.W[0]}, converged, unwarned`);
      if (!direct) assert.ok(named(s, "contacts.left.species['W'], contacts.right.species['W']"), s.warnings.join('; '));
    }
  }
  // Fast enough to resolve: no warning.
  const d = film();
  tracer(d, 'W', 1, 1);
  for (const [side, c] of [['left', 0.5], ['right', 8]]) d.contacts[side].species.W = { type: 'exchange', k: 1, mu: RT * Math.log(c) };
  const s = new Device(d).solve();
  assert.ok(s.converged && s.warnings.length === 0, s.warnings.join('; '));
});

test('where no reservoir reaches, a slow reaction left off equilibrium is warned of: a cycle, and one pool beside another', () => {
  // (U = V = W = U at kf 1e-32 from 1, 2, 9: any two span what the third changes, so no reaction
  // alone was needed, and the solve kept the start, unwarned; at kf 1 it's 4, 4, 4. P = Q in two
  // pools across a blocked face, one at (1, 7), the other at equilibrium (8, 8): one zero affinity
  // among the nodes defeated the one-sign test, unwarned.)
  for (const [k, slow] of [[1e-32, true], [1, false]]) {
    const d = film();
    for (const [name, c] of [['U', 1], ['V', 2], ['W', 9]]) tracer(d, name, c);
    for (const equation of ['U = V', 'V = W', 'W = U']) d.bulkReactions.push({ equation, kf: { m1: k } });
    const dev = new Device(d), s = dev.solve();
    assert.ok(s.converged);
    if (slow) {
      for (const x of [1, 2, 3]) assert.ok(named(s, `bulkReactions[${x}]`), s.warnings.join('; '));
      assert.equal(check(dev, s).items.find((it) => it.name === 'balance').ok, false);
    } else {
      assert.equal(s.warnings.length, 0, s.warnings.join('; '));
      for (const name of ['U', 'V', 'W']) assert.ok(Math.abs(s.c[name][0] - 4) < 1e-9);
    }
  }
  for (const [k, slow] of [[1e-32, true], [1, false]]) {
    const d = film();
    tracer(d, 'P', 1, 8);
    tracer(d, 'Q', 7, 8);
    d.interfaces = [{ phi: 'neutral', species: { P: 'blocked', Q: 'blocked' } }];
    d.bulkReactions.push({ equation: 'P = Q', kf: { m1: k, m2: k } });
    const s = new Device(d).solve();
    assert.ok(s.converged);
    assert.equal(named(s, 'bulkReactions[1]'), slow, s.warnings.join('; '));
    if (!slow) assert.ok(Math.abs(s.c.P[0] - 4) < 1e-9 && Math.abs(s.c.Q[0] - 4) < 1e-9);
  }
});

test("a slow face permeability, and a gate behind a huge R, aren't taken as steady from the first giant step", () => {
  // (W at 1 and 7 in two closed halves behind P = 1e-35 m/s: steady from the first giant step,
  // unwarned; any P > 0 gives 4 and 4. A gate behind R = 1e30 Ω·m², C = 1 µF/m², its source set
  // to 1 V: −0.179 V, converged, unwarned; a steady state passes no current, so V = 1.)
  for (const [P, slow] of [[1e-35, true], [1e-6, false]]) {
    const d = film();
    tracer(d, 'W', 1, 7);
    d.interfaces = [{ phi: 'neutral', species: { W: { type: 'permeability', P } } }];
    const s = new Device(d).solve();
    assert.ok(s.converged);
    assert.equal(named(s, "interfaces[0].species['W']"), slow, s.warnings.join('; '));
    if (!slow) assert.ok(Math.abs(s.c.W[0] - 4) < 1e-9 && Math.abs(s.c.W.at(-1) - 4) < 1e-9);
  }
  for (const [R, slow] of [[1e30, true], [1e20, false]]) {
    const d = film();
    d.ports = [{ name: 'gate', region: 'A', from: 0.2e-6, to: 0.8e-6, V: 0, R, area: 1, capacitance: { C: 1e-6, zeroCharge: 0 } }];
    const dev = new Device(d);
    dev.solve();
    const s = dev.set({ ports: { gate: { V: 1 } } }).solve();
    assert.ok(s.converged);
    assert.equal(named(s, 'gate'), slow, s.warnings.join('; '));
    if (!slow) assert.ok(Math.abs(s.terminals.gate.V - 1) < 1e-9);
  }
});

test("coupled spectators on a nearly full lattice are put back to their amounts in either species order", () => {
  // (P and Q at 4.999 each of a lattice's 10 sites, P's level raised 0.2 RT: swept P first both
  // came back exactly; Q first, the 50 passes left Q 7.0e-5 high.)
  for (const fill of [4.9, 4.999]) {
    for (const order of [['P', 'Q'], ['Q', 'P']]) {
      const d = film();
      for (const name of order) tracer(d, name, fill);
      d.materials.m1.statistics = [{ type: 'lattice', species: ['P', 'Q'], cMax: 10 }];
      const v = new Device(d).solver, i = v.model.species.findIndex((sp) => sp.name === 'P');
      for (let g = 0; g <= v.model.grid.regionEnd[0]; g++) v.u[v.blockOfNode[g] * v.M + 1 + i] += 0.2;
      v.computeConcentrations();
      v._renormalizeSpectators();
      for (const name of order) {
        const k = v.stretches.findIndex((st) => v.model.species[st.species].name === name);
        const off = v.amount(v.stretches[k]) / v.referenceAmounts[k] - 1;
        assert.ok(Math.abs(off) <= 2e-12, `fill ${fill}, ${order.join(' first, ')}: ${name} off by ${off}`);
      }
    }
  }
});
