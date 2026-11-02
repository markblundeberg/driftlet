import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY } from '../src/index.js';
import { build, layer, bath, aqueous, check, metal, ohmic } from '../src/kit.js';
import { choice, device } from '../tools/lawsweep-devices.js';

// What a steady solve keeps, at the edges of the law finder: devices from Gemini 4 Argon's review
// (round 2 of its breadth-first pass), each a combination the steady system once dropped or
// weighed wrongly, so that solve() converged where a long advance() from the same start doesn't go.
// The last five came from a cold review of src/laws.js against its contract alone.

const lib = (ions, more = []) => {
  const l = aqueous(ions, { epsr: 0 });
  l.species.push({ name: 'e-', z: -1 }, ...more);
  return l;
};
const near = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${what}: ${a} vs ${b}`);

test('two floating electrodes sharing a closed ion: the ion is kept, each electrode passing nothing', () => {
  // Fe²⁺ = Fe³⁺ + e⁻ on two electrodes at I = 0, Fe²⁺ from a bath, Fe³⁺ closed in the film: the
  // law, Fe³⁺ with both electrodes' charges, was dropped for having two. (With capacitances, it
  // moves with both voltages, in one's circuit row.)
  for (const cap of [{}, { capacitance: { C: 0.2, zeroCharge: 0 } }]) {
  const def = build({
    T: 298.15,
    library: [lib(['Na+', 'Cl-', 'Fe2+', 'Fe3+'])],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'Cl-': 105, 'Fe2+': 1, 'Fe3+': 1 } }), bath({ 'Na+': 100, 'Cl-': 102, 'Fe2+': 1 }, 'Cl-')],
    ports: ['p1', 'p2'].map((name, k) => ({ name, region: 'film', from: 6e-6 * k, to: 4e-6 + 6e-6 * k, I: 0, terminal: 'e-', area: 1e4, ...cap, reactions: [{ equation: 'Fe2+ = Fe3+ + e-', k0: 1e-3, alpha: 0.5 }] })),
    grid: { hmin: 5e-7, hmax: 5e-7 },
  });
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e5, { tol: 1e-8 });
  assert.ok(s.converged && (cap.capacitance || check(d, s).ok), check(d, s).text); // (with capacitances, its charge wants a finer grid)
  for (const p of ['p1', 'p2']) near(s.terminals[p].V, late.terminals[p].V, 1e-6, `${p}${cap.capacitance ? ', with capacitance' : ''}`);
  }
});

test('adsorbates turning into each other through a closed ion: their sites\' total is kept node by node, and the ion\'s with them', () => {
  // A(ads) + Cl⁻ = B(ads) + e⁻, Cl⁻ closed in the film: θ_A + θ_B stays what it was at each node,
  // and Cl⁻ + θ_B summed. Pivoted on Cl⁻ first, the surface's own law was mixed with it and lost.
  for (const [drive, extra, c0, out] of [
    [{ V: 0.05 }, [], { 'Na+': 100, 'OH-': 99, 'Cl-': 1 }, { 'Na+': 100, 'OH-': 100 }],
    [{ I: 0 }, ['Fe2+', 'Fe3+'], { 'Na+': 100, 'OH-': 104, 'Cl-': 1, 'Fe2+': 1, 'Fe3+': 1 }, { 'Na+': 100, 'OH-': 102, 'Fe2+': 1 }],
  ]) {
    const def = build({
      T: 298.15,
      library: [lib(['Na+', 'OH-', 'Cl-', ...extra])],
      stack: [{}, layer('water', 10e-6, { name: 'film', c0 }), bath(out, 'OH-')],
      ports: [{
        name: 'metal', region: 'film', ...drive, terminal: 'e-', area: 1e4,
        surface: { Aads: { mu0: 0, capacity: 1e-5, theta0: 0.2 }, Bads: { mu0: 3e3, capacity: 1e-5, theta0: 0.1 } },
        reactions: [{ equation: 'Aads + Cl- = Bads + e-', k0: 1e-4, alpha: 0.5 }, ...(extra.length ? [{ equation: 'Fe2+ = Fe3+ + e-', k0: 1e-4, alpha: 0.5 }] : [])],
      }],
      grid: { hmin: 5e-7, hmax: 5e-7 },
    });
    const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e6, { tol: 1e-9 }), p = s.ports[0];
    assert.ok(s.converged);
    p.x.forEach((_, j) => near(p.coverage.Aads[j] + p.coverage.Bads[j], 0.3, 1e-9, `${JSON.stringify(drive)}: θ_A + θ_B at node ${j}`));
    near(s.terminals.metal.V, late.terminals.metal.V, 1e-5, `${JSON.stringify(drive)}: V`);
    near(s.c['Cl-'][0], late.c['Cl-'][0], 1e-5, `${JSON.stringify(drive)}: Cl⁻`);
  }
});

// K⁺ between a held bath of Na⁺Cl⁻ and something at I = 0 that passes K⁺ and has a capacitance:
// what's kept is the K⁺ with that charge, n_K + Q/F, not the K⁺ alone.
const film = { 'Na+': 100, 'K+': 10, 'Cl-': 110 };
const capCell = (end) => build({
  T: 298.15,
  library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })],
  stack: [end.contact ?? {}, layer('water', 1e-6, { name: 'film', c0: film }), bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-')],
  ...(end.port ? { ports: [end.port] } : {}),
  grid: { hmin: 5e-8, hmax: 5e-8 },
});
const kept = (d, s, charge) => {
  const v = d.model.grid.vol;
  return s.c['K+'].reduce((t, c, g) => t + (Number.isFinite(c) ? v[g] * c : 0), 0) + charge(s) / FARADAY;
};

test('a capacitive port at no current that passes an ion: the ion and the port\'s charge are kept together', () => {
  const def = capCell({ port: { name: 'mem', region: 'film', I: 0, area: 1e6, capacitance: { C: 0.1, zeroCharge: -0.05 }, species: { 'K+': { type: 'conductance', G: 1e12, offset: 0 } } } });
  const Q = (s) => s.ports[0].charge;
  const d0 = new Device(def), start = kept(d0, d0.solution(), Q);
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e3, { tol: 1e-9 });
  assert.ok(s.converged);
  near(kept(d, s, Q), start, 1e-9, 'n_K + Q/F');
  near(s.terminals.mem.V, late.terminals.mem.V, 1e-6, 'V');
});

test('a capacitive contact at no current that passes an ion: the ion and the contact\'s charge are kept together', () => {
  const def = capCell({ contact: { I: 0, zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: { 'K+': { type: 'conductance', G: 1e6, offset: 0 } } } });
  const D = (s) => s.contacts.left.D;
  const d0 = new Device(def), start = kept(d0, d0.solution(), D);
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e3, { tol: 1e-9 });
  assert.ok(s.converged);
  near(kept(d, s, D), start, 1e-9, 'n_K + D/F');
  near(kept(d, late, D), start, 1e-9, 'n_K + D/F, the transient'); // (its long steps read the contact's current by its flux)
  near(s.terminals.left.V, late.terminals.left.V, 1e-6, 'V');
});

test('an immobile species that adsorbs: kept node by node with what the surface holds, Γ a θ (mol/m³), not θ', () => {
  const G = 1e-5, area = 1e4;
  const l = lib(['Na+', 'OH-'], [{ name: 'S0', z: 0 }]);
  l.materials.water.species.S0 = { D: 0, mu0: 0, cRef: 1 };
  const def = build({
    T: 298.15,
    library: [l],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'OH-': 100, S0: 1 } }), bath({ 'Na+': 100, 'OH-': 100 }, 'OH-')],
    ports: [{ name: 'metal', region: 'film', V: 0.1, terminal: 'e-', area, surface: { Aads: { mu0: 0, capacity: G, theta0: 0.2 } }, reactions: [{ equation: 'S0 + OH- = Aads + e-', k0: 1e-3, alpha: 0.5 }] }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  const s = new Device(def).solve(), late = new Device(def).advance(1e4, { tol: 1e-9 });
  assert.ok(s.converged);
  for (const sol of [s, late]) sol.ports[0].x.forEach((x, j) => near(sol.c.S0[sol.x.indexOf(x)] + G * area * sol.ports[0].coverage.Aads[j], 1 + G * area * 0.2, 1e-9, `c + Γaθ at node ${j}`));
  near(s.c.S0[0], late.c.S0[0], 1e-5, 'S0');
});

test('a face reaction that makes an immobile species, beside a floating electrode: every balance said once, and the electrode passes nothing', () => {
  const l = lib(['Na+', 'Cl-', 'Ag+'], [{ name: 'S0', z: 0 }]);
  l.materials.water.species.S0 = { D: 0, mu0: 0, cRef: 1 };
  const def = build({
    T: 298.15,
    library: [l],
    stack: [
      {},
      layer('water', 2e-6, { name: 'leftFilm', c0: { 'Na+': 100, 'Cl-': 101, 'Ag+': 1, S0: 1 } }),
      { species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium', 'Ag+': 'blocked', S0: 'blocked' }, reactions: [{ equation: 'Ag+(right) + Cl-(right) = S0(right)', k0: 1e-3, alpha: 0.5 }] },
      layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'Cl-': 101, 'Ag+': 1, S0: 0.1 } }),
      bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-'),
    ],
    ports: [{ name: 'metal', region: 'film', I: 0, terminal: 'e-', area: 1e4, reactions: [{ equation: 'Ag+ + e- = Ag', fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 }] }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e4, { tol: 1e-9 });
  assert.ok(s.converged && check(d, s).ok, check(d, s).text);
  near(s.terminals.metal.V, late.terminals.metal.V, 1e-6, 'V');
  near(s.c['Ag+'].at(-1), late.c['Ag+'].at(-1), 1e-6, 'Ag⁺');
});

const adsorbingOnFloating = () => {
  const G = 1e-5, area = 1e4;
  const l = lib(['Na+', 'OH-', 'Cl-'], [{ name: 'S0', z: 0 }]);
  l.materials.water.species.S0 = { D: 0, mu0: 0, cRef: 1 };
  const def = build({
    T: 298.15,
    library: [l],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'OH-': 1, 'Cl-': 99, S0: 1 } }), bath({ 'Na+': 100, 'OH-': 1, 'Cl-': 99 }, 'Cl-')],
    ports: [{ name: 'metal', region: 'film', from: 0, to: 4e-6, I: 0, terminal: 'e-', area, surface: { Cads: { mu0: 0, capacity: G, theta0: 0.2 } }, reactions: [{ equation: 'S0 + OH- = Cads + e-', k0: 1e-4, alpha: 0.5 }] }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  return { def, G, area };
};

test('an immobile species adsorbing on a floating electrode over part of its region: S + Γaθ kept where the electrode is, S alone elsewhere, and the electrons summed', () => {
  // S0 + OH⁻ = C(ads) + e⁻ at I = 0, the electrode over 4 of the film's 10 µm. Required to share
  // the stretch's nodes, the local law was lost, and the mean coverage drifted from 0.2 to 0.70.
  const { def, G, area } = adsorbingOnFloating();
  const s = new Device(def).solve(), late = new Device(def).advance(1e6, { tol: 1e-9 });
  assert.ok(s.converged && late.converged);
  for (const sol of [s, late]) {
    const p = sol.ports[0];
    near(p.coverage.Cads.reduce((t, x) => t + x, 0) / p.x.length, 0.2, 1e-9, 'the mean coverage');
    p.x.forEach((x, j) => near(sol.c.S0[sol.x.indexOf(x)] + G * area * p.coverage.Cads[j], 1 + G * area * 0.2, 1e-9, `S0 + Γaθ at node ${j}`));
  }
  near(s.terminals.metal.V, late.terminals.metal.V, 1e-6, 'V');
});

test('two contacts at no current passing K⁺, one or both capacitors: the ion kept with both their charges', () => {
  // (From Gemini 4 Argon's third round.) The law takes one contact's circuit row; the other's
  // charge was left out, and with the row's own contact plain, every capacitor's was.
  const K = { 'K+': { type: 'conductance', G: 1e6, offset: 0 } };
  const cap = { I: 0, zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: K };
  for (const left of [{ I: 0, phi: 'neutral', species: K }, cap]) {
    const def = build({
      T: 298.15,
      library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })],
      stack: [left, layer('water', 1e-6, { name: 'film', c0: film }), cap],
      ports: [{ name: 'ref', region: 'film', V: 0, species: { 'Cl-': { type: 'equilibrium', offset: 0 } } }],
      grid: { hmin: 5e-8, hmax: 5e-8 },
    });
    const D = (s) => s.contacts.left.D - s.contacts.right.D;
    const d0 = new Device(def), start = kept(d0, d0.solution(), D);
    const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e3, { tol: 1e-9 });
    assert.ok(s.converged && late.converged);
    near(kept(d, s, D), start, 1e-9, `n_K + Σ D/F, left ${left.phi.type ?? left.phi}`);
    near(s.terminals.right.V, late.terminals.right.V, 1e-6, 'V');
  }
});

test('a charged immobile species beside the face that makes it: where nothing touches it, its concentration stays put when φ moves', () => {
  // (From Gemini 4 Argon's third round.) Ag⁺(right) = S⁺(right) at a face, S⁺ immobile, a floating
  // Ag electrode: at the nodes away from the face, S⁺'s steady balance is 0 = 0, and its level
  // stayed where it was while φ moved with the bath, so S⁺ fell to e^{−0.2/V_T} there (and the face
  // node took the rest of its total). The electrode settled 19 mV off.
  const l = lib(['Na+', 'Cl-', 'Ag+'], [{ name: 'S+', z: 1 }]);
  l.materials.water.species['S+'] = { D: 0, mu0: 0, cRef: 1 };
  const def = (V) => build({
    T: 298.15,
    library: [l],
    stack: [
      {},
      layer('water', 2e-6, { name: 'leftFilm', c0: { 'Na+': 100, 'Cl-': 102, 'Ag+': 1, 'S+': 1 } }),
      { species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium', 'Ag+': 'blocked', 'S+': 'blocked' }, reactions: [{ equation: 'Ag+(right) = S+(right)', k0: 1e-3, alpha: 0.5 }] },
      layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'Cl-': 102, 'Ag+': 1, 'S+': 1 } }),
      bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-', V),
    ],
    ports: [{ name: 'metal', region: 'film', I: 0, terminal: 'e-', area: 1e4, reactions: [{ equation: 'Ag+ + e- = Ag', fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 }] }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  const V0 = new Device(def(0)).solve().terminals.metal.V;
  const s = new Device(def(0.2)).solve(), late = new Device(def(0.2)).advance(1e6, { tol: 1e-9 });
  assert.ok(s.converged && late.converged);
  near(s.terminals.metal.V, late.terminals.metal.V, 1e-6, 'V');
  near(s.terminals.metal.V - V0, 0.2, 1e-6, 'V moved with the bath');
  const g = s.x.length - 3; // (far from the face)
  near(s.c['S+'][g], 1, 1e-9, 'S⁺ away from the face');
});

test('an immobile species linked across a face: carried between the face\'s two nodes, so not kept node by node', () => {
  // D = 0 on both sides, the face passing S⁺ by a conductance: the pair evens out. Kept node by
  // node, solve() left them at 1 and 0.1.
  const l = aqueous(['Na+', 'Cl-'], { epsr: 0 });
  l.species.push({ name: 'S+', z: 1 });
  l.materials.water.species['S+'] = { D: 0, mu0: 0, cRef: 1 };
  const def = build({
    T: 298.15,
    library: [l],
    stack: [{}, layer('water', 2e-6, { name: 'left', c0: { 'Na+': 100, 'Cl-': 101, 'S+': 1 } }), { species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium', 'S+': { type: 'conductance', G: 1e-6 } } },
      layer('water', 4e-6, { name: 'right', c0: { 'Na+': 100, 'Cl-': 100.1, 'S+': 0.1 } }), bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-')],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  const d = new Device(def), g = d.model.grid, s = d.solve(), late = new Device(def).advance(1e8, { tol: 1e-9 });
  assert.ok(s.converged && late.converged);
  for (const at of [g.regionEnd[0], g.regionStart[1]]) near(s.c['S+'][at], late.c['S+'][at], 1e-6, `S⁺ at node ${at}`);
  near(s.c['S+'][g.regionEnd[0]], 0.4, 1e-6, 'S⁺ evened out');
});

test('an immobile species and two adsorbates on one electrode: each node-by-node law gets a row of its own', () => {
  // S⁺ + A(ads) + e⁻ = B(ads), the electrode held: A − S⁺ and B + S⁺ are kept at each node. Both
  // took S⁺'s row, one was lost, and θ_A went to 1 for the transient's 0.3.
  const l = lib(['Na+', 'OH-'], [{ name: 'S+', z: 1 }]);
  l.materials.water.species['S+'] = { D: 0, mu0: 0, cRef: 1 };
  const def = build({
    T: 298.15,
    library: [l],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'OH-': 101, 'S+': 1 } }), bath({ 'Na+': 100, 'OH-': 100 }, 'OH-')],
    ports: [{
      name: 'metal', region: 'film', V: 0.05, terminal: 'e-', area: 1e4,
      surface: { Aads: { mu0: 0, capacity: 1e-5, theta0: 0.2 }, Bads: { mu0: 3e3, capacity: 1e-5, theta0: 0.1 } },
      reactions: [{ equation: 'S+ + Aads + e- = Bads', k0: 1e-4, alpha: 0.5 }],
    }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e5, { tol: 1e-9 });
  assert.ok(s.converged && late.converged);
  s.ports[0].x.forEach((_, j) => near(s.ports[0].coverage.Aads[j], late.ports[0].coverage.Aads[j], 1e-6, `θ_A at node ${j}`));
});

test('an immobile ion that a port at no current links: kept in the port\'s circuit row, not node by node too', () => {
  // S⁺ immobile, its only charged species through a port at I = 0: kept both ways, the port's
  // I = 0 was lost, and solve() settled at 0 V for the transient's −1.42 V.
  const l = aqueous(['Na+', 'Cl-'], { epsr: 0 });
  l.species.push({ name: 'S+', z: 1 });
  l.materials.water.species['S+'] = { D: 0, mu0: 0, cRef: 1 };
  const def = build({
    T: 298.15,
    library: [l],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'Cl-': 101, 'S+': 1 } }), bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-')],
    ports: [{ name: 'p', region: 'film', from: 0, to: 6e-6, I: 0, species: { 'S+': { type: 'conductance', G: 1e3, offset: 0 } } }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e5, { tol: 1e-9 });
  assert.ok(s.converged && late.converged);
  near(s.terminals.p.V, late.terminals.p.V, 1e-6, 'V');
  assert.ok(Math.abs(s.ports[0].flux['S+']) < 1e-15, `the port passes ${s.ports[0].flux['S+']} mol/(m²·s) of S⁺ at I = 0`);
});

test('two laws through what a capacitive contact at no current passes: the contact\'s charge kept once, in its circuit row', () => {
  // K⁺ through a capacitor at I = 0, K⁺ + L⁻ = KL with L⁻ and KL closed: K⁺ − L⁻ and K⁺ + KL both
  // pass K⁺. The second took a balance row without the contact's charge, and solve() failed.
  const l = aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 });
  l.species.push({ name: 'L-', z: -1 }, { name: 'KL', z: 0 });
  l.materials.water.species['L-'] = { D: 1e-9, mu0: 0, cRef: 1 };
  l.materials.water.species.KL = { D: 1e-9, mu0: l.materials.water.species['K+'].mu0 - 3e3, cRef: 1 };
  const def = build({
    T: 298.15,
    library: [l],
    stack: [{ I: 0, zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: { 'K+': { type: 'conductance', G: 1e6, offset: 0 } } },
      layer('water', 1e-6, { name: 'film', c0: { 'Na+': 100, 'K+': 0.5, 'Cl-': 90.5, 'L-': 10, KL: 1 } }), bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-')],
    bulkReactions: [{ equation: 'K+ + L- = KL', kf: { water: 1e3 } }],
    grid: { hmin: 5e-8, hmax: 5e-8 },
  });
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1e3, { tol: 1e-9 });
  assert.ok(s.converged && late.converged);
  near(s.terminals.left.V, late.terminals.left.V, 1e-6, 'V');
  assert.equal(d.solver.constraints.filter((cs) => cs.circuit < 0 && cs.parts.some((p) => d.solver.stretches[p.stretch].species === 1)).length, 0, 'no balance-row law through K⁺');
});

test('a capacitive contact behind a resistance passing a trace of an ion: its current read by the flux, as its charge has no place in the change', () => {
  // Read by the change in K⁺ (a trace, so the change looked the cleaner), its current left out
  // dD/dt, and every step failed.
  const def = build({
    T: 298.15,
    library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })],
    stack: [{ V: 0.1, R: 1e-3, zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: { 'K+': { type: 'conductance', G: 1e6, offset: 0 } } },
      layer('water', 1e-6, { name: 'film', c0: { 'Na+': 100, 'K+': 1e-6, 'Cl-': 100 } }), bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-')],
    grid: { hmin: 5e-8, hmax: 5e-8 },
  });
  const d = new Device(def), s = d.advance(1e-2, { tol: 1e-7 });
  assert.equal(d.solver.timeLaws.length, 0);
  assert.ok(s.converged && s.terminals.left.V < -1, `V = ${s.terminals.left.V}`);
});

test('a node-by-node law in the window of a floating electrode: its row says nothing of the electrode\'s voltage', () => {
  // (From a cold review of src/lawRows.js.) The row S0 + Γaθ takes held the port reaction's
  // ∂/∂V from before, where the law has none: Newton took an extra iteration to its answer.
  const d = new Device(adsorbingOnFloating().def), sv = d.solver, assemble = sv.assemble;
  let checked = 0;
  sv.assemble = function (dt) {
    assemble.call(this, dt);
    if (!(this.steady && dt === Infinity)) return;
    for (const lc of this.localConstraints) {
      for (const g of lc.nodes) {
        const o = this.blockOfNode[g] * this.M + lc.row;
        if (this.loc[o] < 0) continue;
        for (const k of this.floating) assert.equal(this.termB[k][this.rix[o]], 0, `node ${g}: ∂/∂V of terminal ${k}`);
        checked++;
      }
    }
  };
  assert.ok(d.solve().converged);
  assert.ok(checked > 0);
});

test('the law sweep\'s device 19: a capacitive contact at no current beside a floating electrode, its charge held in the contact\'s own circuit row', () => {
  // (A mutation run: holding the charge in the first circuit row whose law left none to hold, the
  // tests all passed, and 85 of the sweep's devices came out wrong, this one at 4.4 mV for −0.09.)
  const def = device(choice(19)), s = new Device(def).solve(), late = new Device(def).advance(1e6, { tol: 1e-9 });
  assert.ok(s.converged && late.converged);
  for (const name of ['left', 'p1']) near(s.terminals[name].V, late.terminals[name].V, 1e-4, `V ${name}`);
});

// Two from Argon's review of 1.6.1, both through a face's surface species, which the law finder
// and the flat-stretch shortcut had read a face reaction without.
const surfaceIons = (names) => ({
  species: names.map(([name, z]) => ({ name, z, cRef: 1000 })),
  materials: { el: { epsr: 0, species: Object.fromEntries(names.map(([name]) => [name, { D: 1e-9, mu0: 0 }])) } },
});
const mechanism = (k0, extra = {}) => ({
  phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 },
  surface: { S: { mu0: 0, capacity: 1e-5, theta0: 0.5 } },
  reactions: [{ equation: 'Ox+ + e- = S', k0, alpha: 0.5 }, { equation: 'S = Red', k0, alpha: 0.5 }],
  ...extra,
});

test('a surface intermediate passing a closed couple on to a counter-electrode: Ox⁺ + Red + S is kept, and S = Red not held at equilibrium', () => {
  // Ox⁺ + e⁻ = S, S = Red at one electrode, Red = Ox⁺ + e⁻ at the other, the film closed but for
  // a reference port at no current. Over the stretches alone, S = Red changed Ox⁺ + Red + S, so it
  // took a row at A = 0 and that law was dropped: −10.98 A/m², where the transient settles at −1.82,
  // with 11× the couple the film started with.
  const def = build({
    library: [surfaceIons([['K+', 1], ['Cl-', -1], ['Ox+', 1], ['Red', 0]]), metal('Au')],
    ports: [{ name: 'ref', region: 'el', from: 4e-6, to: 6e-6, I: 0, terminal: 'Ox+', species: { 'Ox+': 'equilibrium' } }],
    stack: [
      ohmic(0, ['e-']), layer('Au', 1e-6), mechanism(1e-3),
      layer('el', 10e-6, { name: 'el', c0: { 'K+': 100, 'Cl-': 110, 'Ox+': 10, 'Red': 10 } }),
      { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 }, reactions: [{ equation: 'Red = Ox+ + e-', k0: 1e-3, alpha: 0.5 }] },
      layer('Au', 1e-6), ohmic(0.1, ['e-']),
    ],
  });
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1000);
  assert.ok(s.converged && s.warnings.length === 0, s.warnings.join('; '));
  near(s.current, late.current, 1e-9, 'current');
  near(s.interfaces[0].coverage.S, late.interfaces[0].coverage.S, 1e-9, 'coverage');
});

test('a bath species made only from a surface intermediate: not held flat at its bath\'s level', () => {
  // S = Red with Red at its bath's level reads A = 0 over the stretches alone, so Red stood in a
  // group of its own and was held flat at 1000 mol/m³ up to the face while S = Red carried the
  // current: 6582 A/m², where the transient settles at 4044.
  const def = build({
    library: [surfaceIons([['Ox+', 1], ['Cl-', -1], ['Red', 0]]), metal('Au')],
    stack: [
      { V: 0.1, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, layer('Au', 1e-6), mechanism(0.1),
      layer('el', 10e-6, { c0: { 'Ox+': 1000, 'Cl-': 1000, 'Red': 1000 } }), bath({ 'Ox+': 1000, 'Cl-': 1000, 'Red': 1000 }, 'Cl-'),
    ],
  });
  const d = new Device(def), s = d.solve(), late = new Device(def).advance(1000);
  assert.ok(s.converged && s.warnings.length === 0, s.warnings.join('; '));
  near(s.current, late.current, 1e-9, 'current');
  const j = s.c.Red.findIndex(Number.isFinite);
  near(s.c.Red[j], late.c.Red[j], 1e-9, 'Red at the face');
});

test('two steps through a surface intermediate at one face, both at equilibrium: one takes the surface\'s own row, and Ox⁺ + Red + S is kept', () => {
  // Ox⁺ + e⁻ = S (all but stuck) and S = Red at one electrode of a closed cell: a steady state has
  // both at A = 0. Pinned in Ox⁺'s and Red's rows, they left S's own row saying nothing new, and
  // the law had no row: dropped, the solve failed. (Argon's review of 1.6.1, round 2.)
  const surfaceAt = (k0) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, surface: { S: { mu0: 0, capacity: 1e-5, theta0: 0.5 } }, reactions: [{ equation: 'Ox+ + e- = S', k0, alpha: 0.5 }, { equation: 'S = Red', k0: 1e-3, alpha: 0.5 }] });
  const def = build({
    library: [surfaceIons([['K+', 1], ['Cl-', -1], ['Ox+', 1], ['Red', 0]]), metal('Au')],
    ports: [{ name: 'ref', region: 'el', from: 4e-6, to: 6e-6, I: 0, terminal: 'K+', species: { 'K+': 'equilibrium' } }],
    stack: [ohmic(0.15, ['e-']), layer('Au', 1e-6), surfaceAt(1e-18), layer('el', 10e-6, { name: 'el', c0: { 'K+': 100, 'Cl-': 110, 'Ox+': 10, 'Red': 10 } }), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }, layer('Au', 1e-6), ohmic(0, ['e-'])],
  });
  const total = (sol, dev) => {
    let a = 1e-5 * sol.interfaces[0].coverage.S;
    for (const name of ['Ox+', 'Red']) sol.c[name].forEach((c, g) => Number.isFinite(c) && (a += dev.model.grid.vol[g] * c));
    return a;
  };
  const d = new Device(def), s = d.solve();
  assert.ok(s.converged);
  near(total(s, d), 2.05e-4, 1e-12, 'Ox⁺ + Red + S');
  for (const r of s.interfaces[0].rates) assert.ok(Math.abs(r) < 1e-20, `rate ${r}`);
  const late = d.advance(1000); // (and steady: time moves nothing)
  near(late.interfaces[0].coverage.S, s.interfaces[0].coverage.S, 1e-9, 'θ_S');
});

test('a chain through two surface species, every step all but stuck: solved directly, at the equilibrium the same cell reaches when fast (Argon\'s review of 1.6.1, round 3)', () => {
  // Ox⁺ + e⁻ = Red, Red = S1, S1 = S2 at one face of a closed cell, k0 1e-18 to 1e-16. With the face
  // reactions read without their surface species (1.6.1), the solve fell back to time steps and
  // failed. The equilibrium doesn't depend on k0: the same cell 1e15 times faster, stepped there,
  // is the yardstick.
  const cell = (scale) => build({
    library: [surfaceIons([['K+', 1], ['Cl-', -1], ['Ox+', 1], ['Red', 0]]), metal('Au')],
    ports: [{ name: 'ref', region: 'el', from: 4e-6, to: 6e-6, I: 0, terminal: 'K+', species: { 'K+': 'equilibrium' } }],
    stack: [ohmic(0.15, ['e-']), layer('Au', 1e-6), {
      phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0,
      surface: { S1: { mu0: -5000, capacity: 1e-5, theta0: 0.25 }, S2: { mu0: -5000, capacity: 1e-5, theta0: 0.25 } },
      reactions: [{ equation: 'Ox+ + e- = Red', k0: 1e-18 * scale, alpha: 0.5 }, { equation: 'Red = S1', k0: 1e-17 * scale, alpha: 0.5 }, { equation: 'S1 = S2', k0: 1e-16 * scale, alpha: 0.5 }],
    }, layer('el', 10e-6, { name: 'el', c0: { 'K+': 100, 'Cl-': 110, 'Ox+': 10, 'Red': 10 } }), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 }, layer('Au', 1e-6), ohmic(0, ['e-'])],
  });
  const d = new Device(cell(1)), s = d.solve();
  assert.ok(s.converged && d.solver._directSteady(), s.warnings.join('; '));
  const fast = new Device(cell(1e15)).advance(1e4);
  near(s.interfaces[0].coverage.S1, fast.interfaces[0].coverage.S1, 1e-6, 'θ_S1');
  near(s.interfaces[0].coverage.S2, fast.interfaces[0].coverage.S2, 1e-6, 'θ_S2');
  near(s.c.Red.find(Number.isFinite), fast.c.Red.find(Number.isFinite), 1e-6, 'Red');
});
