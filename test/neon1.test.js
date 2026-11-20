// Neon's first review (20 Nov 2026), area 1: the cold start (src/start.js).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { equilibrateFaces } from '../src/start.js';
import { build, layer, metal, ohmic, check, aqueous, semiconductor } from '../src/kit.js';

// Two strictly neutral solutions meeting at a face that passes both ions, sealed, one or both
// under a held gate. The face's t = 0⁺ equilibrium (equilibrateFaces) left the gate's charge out
// of the nodes' neutrality: a gated left node's counter-ions were moved off it, and the first
// step, holding each node neutral, lost them, 5.2e-3 of the ions in a sealed cell, with no
// warning (check() alone said so).
const gated = (regions) => ({
  species: [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }],
  materials: { aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, cRef: 1000, mu0: 0 }, 'Cl-': { D: 2.03e-9, cRef: 1000, mu0: 0 } } } },
  regions: [
    { name: 'L', material: 'aq', length: 1e-6, c0: regions.includes('L') ? { 'Na+': 20, 'Cl-': 10 } : { 'Na+': 20, 'Cl-': 20 } },
    { name: 'R', material: 'aq', length: 1e-6, c0: regions.includes('R') ? { 'Na+': 100, 'Cl-': 110 } : { 'Na+': 100, 'Cl-': 100 } },
  ],
  interfaces: [{ phi: 'neutral', species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium' } }],
  contacts: { left: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } }, right: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } } },
  ports: regions.map((region, k) => ({ name: `g${k}`, region, V: -0.5 + 0.3 * k, area: 1e6, capacitance: { C: 1 } })),
  grid: { minCells: 8 },
});
for (const regions of [['L'], ['R'], ['L', 'R']]) {
  test(`a gated face's t = 0⁺ keeps every ion and the gates' charge (gates on ${regions.join(' and ')}; Neon 1.1)`, () => {
    const steady = new Device(gated(regions)).solve();
    const d = new Device(gated(regions));
    const t = d.advance(0.1);
    assert.ok(t.converged);
    for (let k = 0; k < regions.length; k++) {
      const Q = t.ports[k].charge, Qs = steady.ports[k].charge;
      assert.ok(Math.abs(Q - Qs) < 1e-9 * Math.abs(Qs), `gate ${k}: ${Q} after 0.1 s, ${Qs} C/m² steady`);
    }
    assert.match(check(d, t).text, /ok {3}conservation/);
    // and the face's t = 0⁺ (equilibrateFaces) is neutral at every node, ρ + aσ = 0, gates and all
    const e = new Device(gated(regions)), def = gated(regions);
    e.solution();
    const sv = e.solver, { grid } = sv.model;
    equilibrateFaces(sv);
    for (let g = 0; g < grid.nNodes; g++) {
      const r = grid.nodeRegion[g], b = sv.blockOfNode[g];
      let rho = FARADAY * (sv.c[2 * g] - sv.c[2 * g + 1]);
      def.ports.forEach((p) => {
        if (p.region === def.regions[r].name) rho += p.area * p.capacitance.C * (p.V - sv.VT * sv.u[b * sv.M]);
      });
      assert.ok(Math.abs(rho) < 1e-9 * FARADAY * (sv.c[2 * g] + sv.c[2 * g + 1]), `node ${g} (${def.regions[r].name}): ρ + aσ = ${rho} C/m³ at t = 0⁺`);
    }
  });
}

// Area 2: a group of stretches held flat at its one contact's levels as an equilibrium, where a
// reaction among surface species alone (S = F, a Tafel step) keeps it from one: the check left
// those reactions out, and the surface's terms in the others, so A = S, S = F with F held 25
// kJ/mol below A solved to no flux and A flat at its bath's 1000 mol/m³, with no warning, while
// each time step turned A into F at 9.8e-3 mol/(m²·s) (Neon 2.1).
test('a catalytic face, A = S then S = F to a held product: the steady solve carries what the transient settles to, and only an equilibrium is held flat (Neon 2.1)', () => {
  const lib = { species: [{ name: 'A', z: 0, cRef: 1000 }], materials: { aq: { epsr: 0, species: { A: { D: 1e-9, mu0: 0 } } } } };
  const cell = (dG) =>
    build({
      T: 298.15, library: [lib, metal('Pt')],
      stack: [
        { phi: 'neutral', species: { A: { type: 'equilibrium', mu: 0 } } },
        layer('aq', 1e-6),
        { surface: { S: { mu0: 0, capacity: 1e-5, theta0: 0.1 } }, reactions: [{ equation: 'A = S', k0: 1e-2 }, { equation: 'S = F', fixed: { F: dG }, k0: 1e-2 }] },
        layer('Pt', 1e-7),
        ohmic(0, ['e-']),
      ],
      grid: { minCells: 10 },
    });
  const d = new Device(cell(-25000)), s = d.solve(), t = new Device(cell(-25000)).advance(10);
  assert.ok(s.converged && t.converged);
  const [r, rt] = [s.interfaces[0].rates[1], t.interfaces[0].rates[1]];
  assert.ok(r > 9e-3 && Math.abs(r / rt - 1) < 1e-9, `steady ${r}, after 10 s ${rt} mol/(m²·s)`);
  assert.match(check(d, s).text, /ok {3}balance/);
  const eq = new Device(cell(0)); // (at F's own level the group is at equilibrium: still held flat)
  assert.deepEqual(eq.solver.flatStretches.map((st) => st.species), [0]);
  assert.equal(eq.solve().interfaces[0].rates[1], 0);
});

// Area 1's cold start, written another way round (Neon 1.3, 1.5): φ̂ ran from the left, set by the
// first thing there that set it, and a stack mirrored, or its ports listed in another order,
// started volts away and was refused where the other converged.
const fe = (() => {
  const lib = aqueous(['Fe3+', 'Fe2+', 'Cl-'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  return lib;
})();
const couple = { equation: 'Fe3+ + e- = Fe2+', k0: 1e-2, alpha: 0.5 };
const solA = layer('water', 1e-6, { name: 'A', c0: { 'Fe3+': 1, 'Fe2+': 1, 'Cl-': 5 } });
const solB = layer('water', 1e-6, { name: 'B', c0: { 'Fe3+': 2, 'Fe2+': 1, 'Cl-': 8 } });
const coldPhi = (def, region) => {
  const d = new Device(def), s = d.solution(), g = d.grid.regionStart[def.regions.findIndex((r) => r.name === region)];
  return s.phi[g];
};

test('a held electrode in one solution and a floating one beside it: the start is the held one\'s, whichever is listed first (Neon 1.3)', () => {
  const we = { name: 'we', region: 'A', from: 0, to: 0.4e-6, V: 5, terminal: 'e-', area: 1e4, reactions: [couple] };
  const sense = { name: 'sense', region: 'A', from: 0.6e-6, to: 1e-6, I: 0, terminal: 'e-', area: 1e4, reactions: [couple] };
  const cell = (ports) => build({ T: 298.15, library: [fe], stack: [{}, solA, {}], ports, grid: { minCells: 10 } });
  const [a, b] = [cell([we, sense]), cell([sense, we])];
  assert.ok(Math.abs(coldPhi(a, 'A') - coldPhi(b, 'A')) < 1e-9, `cold φ ${coldPhi(a, 'A')}, listed the other way ${coldPhi(b, 'A')} V`);
  const s = new Device(b).solve();
  assert.ok(s.converged && Math.abs(s.terminals.sense.V - 5) < 1e-6, `[sense, we]: ${s.converged}, sense at ${s.terminals.sense.V} V`);
});

test('an electrode held in the second of two solutions sets the first one\'s start too, as mirrored (Neon 1.3)', () => {
  const port = { name: 'we', region: 'B', V: 10, terminal: 'e-', area: 1e4, reactions: [couple] };
  const ab = build({ T: 298.15, library: [fe], stack: [{}, solA, solB, {}], ports: [port], grid: { minCells: 10 } });
  const ba = build({ T: 298.15, library: [fe], stack: [{}, solB, solA, {}], ports: [port], grid: { minCells: 10 } });
  assert.ok(Math.abs(coldPhi(ab, 'A') - coldPhi(ba, 'A')) < 1e-9, `cold φ in A ${coldPhi(ab, 'A')}, mirrored ${coldPhi(ba, 'A')} V`);
  assert.ok(new Device(ab).solve().converged);
});

test('a floating gated contact sets no level: the electrode across the cell sets the start, as mirrored (Neon 1.3)', () => {
  const cap = { I: 0, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 };
  const face = { reactions: [couple] };
  const pt = layer('Pt', 1e-7);
  const one = build({ T: 298.15, library: [aqueous(['Fe3+', 'Fe2+', 'Cl-'], { epsr: 0 }), metal('Pt')], stack: [cap, solA, solB, face, pt, ohmic(2, ['e-'])], grid: { minCells: 10 } });
  const two = build({ T: 298.15, library: [aqueous(['Fe3+', 'Fe2+', 'Cl-'], { epsr: 0 }), metal('Pt')], stack: [ohmic(2, ['e-']), pt, face, solB, solA, cap], grid: { minCells: 10 } });
  assert.ok(Math.abs(coldPhi(one, 'A') - coldPhi(two, 'A')) < 1e-9, `cold φ in A ${coldPhi(one, 'A')}, mirrored ${coldPhi(two, 'A')} V`);
});

test('a gate over part of a strictly neutral solution sets the whole region\'s start, either side of it (Neon 1.5)', () => {
  const cell = (from, to) => ({
    species: [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }],
    materials: { aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, cRef: 1000, mu0: 0 }, 'Cl-': { D: 2.03e-9, cRef: 1000, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'aq', length: 1e-6, c0: { 'Na+': 100, 'Cl-': 100 } }],
    contacts: { left: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } }, right: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } } },
    ports: [{ name: 'g', region: 'a', from, to, V: 2, area: 1e6, capacitance: { C: 1 } }],
    grid: { minCells: 10, hmax: 1e-7 },
  });
  for (const [from, to] of [[0.2e-6, 0.5e-6], [0.5e-6, 0.8e-6]]) {
    const phi = new Device(cell(from, to)).solution().phi;
    assert.ok(phi.every((p) => Math.abs(p - 2) < 1e-9), `gate on [${from}, ${to}]: cold φ ${[...phi].map((p) => p.toFixed(3))}`);
  }
});

test('a hole-only layer with 20 V across it starts with finite densities, and a contact whose own density is beyond any number is refused, saying why (Neon 1.2)', () => {
  // (the start's neutrality search clamps its exponents at ±700: the charge was flat there, read as
  // a fixed charge's floor, and the search stopped on its first step, starting at c = ∞)
  const layerOf = (contacts) => ({
    species: [{ name: 'h+', z: 1 }],
    materials: { d: { epsr: 10, species: { 'h+': { D: 1e-5, cRef: 10, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'd', length: 200e-9 }],
    contacts,
    grid: { minCells: 40 },
  });
  const fed = new Device(layerOf({ left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { 'h+': 'blocked' } }, right: { V: 20, terminal: 'h+', phi: 'pinned', zeroCharge: 0, species: { 'h+': 'equilibrium' } } }));
  assert.ok(fed.solution().c['h+'].every(Number.isFinite));
  assert.ok(fed.solve().converged);
  const absurd = new Device(layerOf({ left: { V: 20, terminal: 'h+', phi: 'pinned', zeroCharge: 20, species: { 'h+': 'equilibrium' } }, right: { phi: 'neutral', species: { 'h+': 'blocked' } } }));
  assert.throws(() => absurd.solution(), /contacts\.left: holds h\+ at cRef·e\^7\d\d at its end, beyond any density/);
});

test("a face between strictly neutral solutions finds its t = 0⁺ where an ion's partition is beyond a double's range (Neon 1.4)", () => {
  // (K = e^(Δμ°/RT + z D) overflowed, K·c_R was ∞·0 = NaN, and the face was left as c0 had it)
  const sp = (m) => ({ 'Na+': { D: 1.33e-9, cRef: 1000, mu0: m }, 'Cl-': { D: 2.03e-9, cRef: 1000, mu0: 0 } });
  const d = new Device({
    species: [{ name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }, { name: 'A-', z: -1 }],
    materials: { L: { epsr: 0, species: { ...sp(0), 'A-': { D: 0, cRef: 1000, mu0: 0 } } }, R: { epsr: 0, species: sp(2e6) } },
    regions: [{ name: 'L', material: 'L', length: 1e-6, c0: { 'Na+': 20, 'Cl-': 10, 'A-': 10 } }, { name: 'R', material: 'R', length: 1e-6, c0: { 'Na+': 10, 'Cl-': 10 } }],
    interfaces: [{ phi: 'neutral', species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium', 'A-': 'blocked' } }],
    contacts: { left: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' } }, right: { phi: 'neutral', species: { 'Na+': 'blocked', 'Cl-': 'blocked' } } },
    ports: [{ name: 'g', region: 'R', V: 0, area: 1e6, capacitance: { C: 1 } }],
    grid: { minCells: 8 },
  });
  d.solution();
  const s = d.solver, M = s.M, gL = d.grid.regionEnd[0], gR = d.grid.regionStart[1], bL = s.blockOfNode[gL] * M, bR = s.blockOfNode[gR] * M;
  const vol = d.grid.vol, before = [0, 1].map((i) => vol[gL] * s.c[gL * 3 + i] + vol[gR] * s.c[gR * 3 + i]);
  equilibrateFaces(s);
  assert.ok(s.u.every(Number.isFinite));
  for (const i of [0, 1]) {
    assert.ok(Math.abs(s.u[bL + 1 + i] - s.u[bR + 1 + i]) < 1e-6, `species ${i}: η ${s.u[bL + 1 + i]} | ${s.u[bR + 1 + i]}, not equal across the face`);
    const after = vol[gL] * s.c[gL * 3 + i] + vol[gR] * s.c[gR * 3 + i];
    assert.ok(Math.abs(after / before[i] - 1) < 1e-12, `species ${i}: amount ${before[i]} → ${after}`);
  }
});

// Area 4: what check() compares (src/checks.js).
test("check()'s grid item compares each species' flux, so a grid too coarse where the currents cancel fails it (Neon 4.2)", () => {
  // (A uniformly lit n⁺-Si slab between two ohmic contacts at 0 V: I = 0 by symmetry, each carrier
  // leaving both ends. On 0.5–2 µm cells the holes collected are 5.0% off the resolved 8.371e-4
  // mol/(m²·s), and the grid item compared the current alone, ~1e-10 of nothing: ok.)
  const slab = (grid) =>
    build({
      T: 300, library: [semiconductor('Si')],
      stack: [ohmic(0), layer('Si', 20e-6, { donors: units.perCm3(1e19) }), ohmic(0)],
      bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }, { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: 1e3 } }],
      grid,
    });
  const d = new Device(slab({ hmin: 5e-7, hmax: 2e-6 }));
  assert.match(check(d, d.solve()).text, /FAIL grid: .* h\+ in at the (left|right) contact/);
});

test("check()'s grid item compares a reactant's flux as it is, not against what it could carry (Neon 4.3)", () => {
  // (A against a flow, A = B to a held product: 12.7% off on 4 cells, a 9.2% change read against
  // 1e-3 of D c/L as 4.7e-3, ok; the floor is for spectators.)
  const RT = GAS_CONSTANT * 298.15;
  const d = new Device({
    species: [{ name: 'A', z: 0 }],
    materials: { w: { epsr: 0, species: { A: { D: 1e-9, mu0: 0, cRef: 1000 } } } },
    regions: [{ material: 'w', length: 10e-6, velocity: -5e-4, c0: { A: 1000 } }],
    contacts: { left: { species: { A: { type: 'equilibrium', mu: 0 } }, phi: 'neutral' }, right: { species: {}, phi: 'neutral' } },
    bulkReactions: [{ equation: 'A = B', fixed: { B: -6 * RT * Math.log(10) }, kf: { w: 2.5e-3 } }],
    grid: { minCells: 4 },
  });
  assert.match(check(d, d.solve()).text, /FAIL grid: .* A in at the left contact/);
});

test("check() knows a copy of a solution of neutral species alone from before a set() (Neon 4.4)", () => {
  // (φ is undefined throughout, and the stale-state test compared φ alone: a copy from before kf
  // doubled, its flux half the device's, passed with refine: false)
  const def = {
    species: [{ name: 'A', z: 0 }],
    materials: { w: { epsr: 0, species: { A: { D: 1e-9, mu0: 0, cRef: 1000 } } } },
    regions: [{ material: 'w', length: 10e-6, velocity: -5e-4, c0: { A: 1000 } }],
    contacts: { left: { species: { A: { type: 'equilibrium', mu: 0 } }, phi: 'neutral' }, right: { species: {}, phi: 'neutral' } },
    bulkReactions: [{ equation: 'A = B', fixed: { B: -34000 }, kf: { w: 2.5e-3 } }],
  };
  const d = new Device(def), old = structuredClone(d.solve());
  d.set({ bulkReactions: [{ ...def.bulkReactions[0], kf: { w: 5e-3 } }] });
  d.solve();
  assert.match(check(d, old, { refine: false }).text, /FAIL solution: not the device's present state \(its A differs/);
});

test("a closed trace beside an abundant species whose starting profile dips at the region's left end keeps its amount (Neon 2.4)", () => {
  // (the stretches' starting amounts, which order the conserved basis, read a profile at its
  // region's left end: B, 1000 mol/m³ but for a dip to 1e-12 at x = 0, counted as 1e-18 mol/m²,
  // so A − B and B + C were kept, and A + C, 2e-17 mol/m², came out 2.7% off with no warning;
  // at a 1e-8 trace the solve failed. The mirror image was exact.)
  const names = ['A', 'B', 'C'], L = 1e-6;
  const lib = { species: names.map((name) => ({ name, z: 0, cRef: 1000 })), materials: { film: { epsr: 0, species: Object.fromEntries(names.map((s) => [s, { D: 1e-9, mu0: 0 }])) } } };
  for (const trace of [1e-11, 1e-8]) {
    const B = { x: [0, L / 10, L], values: [1e-12, 1000, 1000] };
    const d = new Device(build({ library: [lib], stack: [{}, layer('film', L, { c0: { A: trace, C: trace, B } }), {}], bulkReactions: [{ equation: 'A + B = C', kf: { film: 1e3 } }] }));
    const s = d.solve(), amount = (sp) => s.c[sp].reduce((t, c, g) => t + d.grid.vol[g] * c, 0);
    assert.ok(s.converged, `trace ${trace}`);
    const AC = amount('A') + amount('C');
    assert.ok(Math.abs(AC / (2 * trace * L) - 1) < 1e-9, `trace ${trace}: A + C ${AC}, started ${2 * trace * L} mol/m²`);
  }
});

test('a law of several species on one lattice, beside diffuse layers, on a coarse grid: room for every layer each part reaches (Neon 2.2)', () => {
  // (each law's derivative had room for each layer once; lawAmount writes each layer once per
  // part that reaches it: Pt | four interconverting cations and Cl⁻ on one lattice | Pt, Stern and
  // diffuse layers, threw "needs 104 entries, its workspace holds 100" on 1 to 3 cells)
  const names = ['A+', 'B+', 'C+', 'D+', 'Cl-'];
  const el = {
    species: names.map((name) => ({ name, z: name === 'Cl-' ? -1 : 1, cRef: 1000 })),
    materials: { el: { epsr: 0, statistics: [{ type: 'lattice', species: names, cMax: 5000 }], species: Object.fromEntries(names.map((s) => [s, { D: 1e-9, mu0: 0 }])) } },
  };
  const c0 = { 'A+': 100, 'B+': 100, 'C+': 100, 'D+': 100, 'Cl-': 400 };
  const stern = { phi: { type: 'capacitive', C: 0.2 }, diffuse: { epsr: 78 }, zeroCharge: 0, species: Object.fromEntries(names.map((s) => [s, 'blocked'])) };
  for (const cells of [1, 2, 3]) {
    const def = build({
      library: [el, metal('Pt')],
      stack: [ohmic(0.1, ['e-']), layer('Pt', 1e-7), stern, layer('el', 1e-6, { c0 }), stern, layer('Pt', 1e-7), ohmic(0, ['e-'])],
      bulkReactions: [{ equation: 'A+ = B+', kf: { el: 1 } }, { equation: 'B+ = C+', kf: { el: 1 } }, { equation: 'C+ = D+', kf: { el: 1 } }],
      grid: { minCells: cells },
    });
    assert.ok(new Device(def).solve().converged, `${cells} cells`);
  }
});

test("check()'s balance counts a ledger as moving where a moving one shows its reaction's rate is real (Neon's reply on 4.1)", () => {
  // (A = S, S = F at a face at k0 1e-10: S's ledger balanced its two steps at 9.9e-11 mol/(m²·s), 1e-10
  // of its scale; a solve that dropped A's flux (2.1, before its fix) left A's ledger unbalanced by
  // all of it, below the 1e-8 at which a ledger moves, and balance said ok. Here the flux is zeroed
  // in a copy, as that solve had it.)
  const def = {
    T: 300,
    species: [{ name: 'A', z: 0 }, { name: 'e-', z: -1 }],
    materials: { aq: { epsr: 0, species: { A: { D: 1e-9, mu0: 0, cRef: 1000 } } }, pt: { conductor: { species: 'e-', conductivity: 1e7 } } },
    regions: [{ name: 'L', material: 'aq', length: 1e-6, c0: { A: 1000 }, grid: { hmax: 1e-7, minCells: 10 } }, { name: 'Pt', material: 'pt', length: 1e-7 }],
    interfaces: [{ phi: 'neutral', species: { A: 'blocked', 'e-': 'blocked' }, surface: { S: { mu0: 0, capacity: 1e-5, theta0: 0.5 } }, reactions: [{ equation: 'A = S', k0: 1e-10, alpha: 0.5 }, { equation: 'S = F', fixed: { F: -25000 }, k0: 1e-10, alpha: 0.5 }] }],
    contacts: { left: { phi: 'neutral', species: { A: { type: 'equilibrium', mu: 0 }, 'e-': 'blocked' } }, right: { V: 0, phi: 'bulk', species: { 'e-': { type: 'equilibrium', offset: 0 }, A: 'blocked' } } },
  };
  const d = new Device(def), s = d.solve();
  assert.match(check(d, s, { refine: false }).text, /ok {3}balance: 2 ledgers sum to zero/);
  const dropped = { ...s, contacts: { ...s.contacts, left: { ...s.contacts.left, flux: { ...s.contacts.left.flux, A: 0 } } } };
  assert.match(check(d, dropped, { refine: false }).text, /FAIL balance: .* not a steady state/);
});
