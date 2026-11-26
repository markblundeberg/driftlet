// Contract tests for the anchor-outward cold start (src/start.js):
//   - Direct mirror symmetry of the cold start state itself (φ, V, c, termV) across architectures
//   - Every Known item removed from src/start.js's contract (verified failing on c87b373)
//   - Interior metal right-face anchoring, neutral-family divalent salt beside a high-V gate,
//     round-17 mirrored face + reaction-order independence, and stretch source precedence.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { aqueous, bath, build, half, layer, metal, ohmic, semiconductor } from '../src/kit.js';

const T0 = 298.15;
const RT = GAS_CONSTANT * T0;
const VT = RT / FARADAY;

// Mirror a plain definition (regions, interfaces, contacts, ports) across x → L − x.
function mirrorDef(def) {
  const m = structuredClone(def);
  const L = m.regions.reduce((s, r) => s + r.length, 0);
  const oldStarts = [];
  let acc = 0;
  for (const r of m.regions) { oldStarts.push(acc); acc += r.length; }
  const R = m.regions.length;
  m.regions = m.regions.slice().reverse().map((r) => {
    const nr = { ...r };
    if (nr.grid?.hLeft !== undefined || nr.grid?.hRight !== undefined) {
      nr.grid = { ...nr.grid, hLeft: nr.grid.hRight, hRight: nr.grid.hLeft };
    }
    return nr;
  });
  if (m.interfaces) {
    m.interfaces = m.interfaces.slice().reverse().map((itf) => {
      if (!itf) return itf;
      const ni = { ...itf };
      if (ni.dipole) ni.dipole = -ni.dipole;
      if (ni.reactions) {
        ni.reactions = ni.reactions.map((rx) => ({
          ...rx,
          equation: rx.equation
            .replace(/\(left\)/g, '(__R__)')
            .replace(/\(right\)/g, '(left)')
            .replace(/\(__R__\)/g, '(right)'),
        }));
      }
      if (ni.surfaceSide) ni.surfaceSide = ni.surfaceSide === 'left' ? 'right' : 'left';
      return ni;
    });
  }
  m.contacts = { left: m.contacts?.right ?? {}, right: m.contacts?.left ?? {} };
  if (m.ports) {
    m.ports = m.ports.map((p) => {
      const np = { ...p };
      const rIdx = typeof p.region === 'number' ? p.region : def.regions.findIndex((r) => r.name === p.region);
      const rLen = def.regions[rIdx].length;
      if (typeof p.region === 'number') np.region = R - 1 - p.region;
      if (p.from !== undefined || p.to !== undefined) {
        const a = p.from ?? 0, b = p.to ?? rLen;
        np.from = rLen - b;
        np.to = rLen - a;
      }
      return np;
    });
  }
  return m;
}

// Assert that the cold start state of `def` and `mirrorDef(def)` are exact mirror images.
function assertColdMirrors(def, label, tol = 1e-12) {
  const a = new Device(structuredClone(def)).solution();
  const b = new Device(mirrorDef(def)).solution();
  const N = a.x.length;
  assert.equal(b.x.length, N, `${label}: grid node count differs`);
  for (let g = 0; g < N; g++) {
    const gm = N - 1 - g;
    if (Number.isFinite(a.phi[g]) || Number.isFinite(b.phi[gm])) {
      assert.ok(
        Math.abs(a.phi[g] - b.phi[gm]) <= tol * Math.max(1, Math.abs(a.phi[g]), Math.abs(b.phi[gm])),
        `${label}: cold φ[${g}] = ${a.phi[g]} vs mirror[${gm}] = ${b.phi[gm]}`,
      );
    }
    for (const sp of a.species) {
      const va = a.V[sp.name]?.[g], vb = b.V[sp.name]?.[gm];
      if (Number.isFinite(va) || Number.isFinite(vb)) {
        assert.ok(
          Math.abs(va - vb) <= tol * Math.max(1, Math.abs(va), Math.abs(vb)),
          `${label}: cold V.${sp.name}[${g}] = ${va} vs mirror[${gm}] = ${vb}`,
        );
      }
      const ca = a.c[sp.name]?.[g], cb = b.c[sp.name]?.[gm];
      if (Number.isFinite(ca) || Number.isFinite(cb)) {
        assert.ok(
          Math.abs(ca - cb) <= tol * Math.max(1e-12, Math.abs(ca), Math.abs(cb)),
          `${label}: cold c.${sp.name}[${g}] = ${ca} vs mirror[${gm}] = ${cb}`,
        );
      }
    }
  }
  assert.ok(Math.abs(a.terminals.left.V - b.terminals.right.V) <= tol * Math.max(1, Math.abs(a.terminals.left.V)), `${label}: left.V vs mirror right.V`);
  assert.ok(Math.abs(a.terminals.right.V - b.terminals.left.V) <= tol * Math.max(1, Math.abs(a.terminals.right.V)), `${label}: right.V vs mirror left.V`);
  for (const p of a.ports) {
    const pm = b.ports.find((q) => q.name === p.name);
    assert.ok(Math.abs(p.V - pm.V) <= tol * Math.max(1, Math.abs(p.V)), `${label}: port ${p.name}.V`);
  }
}

test('cold start is directly mirror-symmetric across architectures', () => {
  // 1. Closed left end, bath referenced to Cl⁻ on the right (old Known item: ~55 RT off).
  const libSalt = aqueous(['Na+', 'Cl-'], { epsr: 0 });
  const cSalt = { 'Na+': 100, 'Cl-': 100 };
  const closedLeftBathRight = build({
    library: [libSalt],
    stack: [{}, layer('water', 10e-6, { c0: cSalt }), bath(cSalt, 'Cl-')],
    grid: { minCells: 10 },
  });
  assertColdMirrors(closedLeftBathRight, 'closed | water | bath(Cl-)');

  // 2. Two held contacts at different V on a semiconductor stack with a dipole.
  const semiDiff = build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [ohmic(0), layer('Si', 400e-9, { acceptors: units.perCm3(1e17) }), { dipole: 0.25 }, layer('Si', 600e-9, { donors: units.perCm3(2e17) }), ohmic(0.4)],
    grid: { minCells: 12 },
  });
  assertColdMirrors(semiDiff, 'p-Si | dipole | n-Si at 0.4 V');

  // 3. Pure conductor stretch between two held contacts at different V (ohmic ramp).
  const wire = build({
    library: [metal('Cu')],
    stack: [ohmic(0, ['e-']), layer('Cu', 1e-3), ohmic(0.2, ['e-'])],
    grid: { minCells: 10 },
  });
  assertColdMirrors(wire, 'Cu wire 0 → 0.2 V');

  // 4. Two unlike reversible electrodes across an electrolyte (ohmic gap interpolation).
  const libCell = aqueous(['Cu2+', 'Zn2+', 'Na+', 'Cl-'], { epsr: 0 });
  const cCell = { 'Cu2+': 10, 'Zn2+': 10, 'Na+': 100, 'Cl-': 140 };
  const cell = build({
    library: [libCell, metal('Cu'), metal('Zn')],
    stack: [
      ohmic(0, ['e-']), layer('Cu', 1e-6),
      { reactions: [{ equation: 'Cu2+ + 2 e- = Cu(s)', fixed: { 'Cu(s)': 0 }, k0: 1e-3, alpha: 0.5 }] },
      layer('water', 20e-6, { c0: cCell }),
      { reactions: [{ equation: 'Zn2+ + 2 e- = Zn(s)', fixed: { 'Zn(s)': 0 }, k0: 1e-3, alpha: 0.5 }] },
      layer('Zn', 1e-6), ohmic(-1.1, ['e-']),
    ],
    grid: { minCells: 10 },
  });
  assertColdMirrors(cell, 'Cu | electrolyte | Zn');

  // 5. Interior metal whose only reaction is on its right face, behind a blocking face on its left.
  const libFe = aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 });
  const cFe = { 'K+': 100, 'Cl-': 125, 'Fe3+': 5, 'Fe2+': 5 };
  const interiorRight = build({
    library: [libFe, metal('Pt')],
    stack: [
      { phi: 'neutral' },
      layer('water', 10e-6, { name: 'far', c0: cFe }),
      { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 },
      layer('Pt', 2e-6, { name: 'mid' }),
      { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1e-3, alpha: 0.5 }] },
      layer('water', 10e-6, { name: 'near', c0: cFe }),
      bath(cFe, 'Cl-'),
    ],
    grid: { minCells: 10 },
  });
  assertColdMirrors(interiorRight, 'water | cap | Pt | Fe rx | water | bath');

  // 6. Strictly neutral gated region on right, ungated region on left.
  const gatedRight = build({
    library: [libSalt],
    stack: [
      {},
      layer('water', 10e-6, { name: 'leftW', c0: cSalt }),
      { phi: 'neutral', species: { 'Na+': 'equilibrium', 'Cl-': 'equilibrium' } },
      layer('water', 10e-6, { name: 'rightW', c0: cSalt }),
      {},
    ],
    ports: [{ name: 'g0', region: 'rightW', V: 3.5, area: 1e5, capacitance: { C: 0.2, zeroCharge: 0 } }],
    grid: { minCells: 8 },
  });
  assertColdMirrors(gatedRight, 'ungated | gated(3.5 V)');
});

test('gauge behind a closed left end with a Cl⁻-referenced bath on the right starts at the bath φ (not 0 V)', () => {
  // Removed Known item: previously started ~55 RT (~1.42 V) off at φ = 0.
  const lib = aqueous(['Na+', 'Cl-'], { epsr: 0 });
  const c = { 'Na+': 100, 'Cl-': 100 };
  const phiBath = (lib.materials.water.species['Cl-'].mu0 + RT * Math.log(100 / 1000)) / FARADAY;
  const d = new Device(build({
    library: [lib],
    stack: [{}, layer('water', 10e-6, { c0: c }), bath(c, 'Cl-')],
    grid: { minCells: 10 },
  }));
  const cold = d.solution();
  assert.ok(Math.abs(cold.phi[0] - phiBath) < 1e-9, `cold φ[0] = ${cold.phi[0]} V, expected ${phiBath} V`);
  assert.ok(Math.abs(cold.phi.at(-1) - phiBath) < 1e-9, `cold φ[end] = ${cold.phi.at(-1)} V, expected ${phiBath} V`);
  const s = d.solve();
  assert.ok(s.converged && s.iterations <= 2, `iters = ${s.iterations}`);
});

test('c0 given on a contact-fed carrier does not move the gauge (resistor at rest)', () => {
  // Removed Known item: a resistor at rest with c0 given on its majority carrier started ~0.15–0.2 V
  // off equilibrium and took ~17× as many steps.
  const Nd = units.perCm3(1e17);
  const lib = semiconductor('Si');
  const ni2 = lib.materials.Si.species['e-'].cRef * lib.materials.Si.species['h+'].cRef *
    Math.exp(-(lib.materials.Si.species['e-'].mu0 + lib.materials.Si.species['h+'].mu0) / (GAS_CONSTANT * 300));
  const pEq = ni2 / Nd;
  for (const [label, c0] of [['e- alone', { 'e-': Nd }], ['both carriers', { 'e-': Nd, 'h+': pEq }]]) {
    const d = new Device(build({
      T: 300,
      library: [lib],
      stack: [ohmic(0), layer('Si', 1e-6, { donors: Nd, c0 }), ohmic(0)],
      grid: { minCells: 12 },
    }));
    const cold = d.solution();
    const mid = cold.x.length >> 1;
    assert.ok(Math.abs(cold.V['e-'][mid]) < 1e-6, `${label}: cold V_e = ${cold.V['e-'][mid]} V`);
    const s = d.solve();
    assert.ok(s.converged && s.iterations <= 3, `${label}: iters = ${s.iterations}`);
  }
});

test('round-17 mirrored face (bath | water with c0 | Pt) and reaction-order independence', () => {
  // bath(Cl-) | water (with c0 given or dropped) | Pt with two reactions [volmer, iron] vs [iron, volmer]:
  // cold φ must not shift when c0 matches the bath, when mirrored, or when reactions are reordered.
  const lib = aqueous(['H+', 'Fe3+', 'Fe2+', 'K+', 'Cl-'], { epsr: 0 });
  const c = { 'H+': 10, 'Fe3+': 5, 'Fe2+': 5, 'K+': 100, 'Cl-': 135 };
  const volmer = { equation: '2 H+ + 2 e- = H2(g)', fixed: { 'H2(g)': 0 }, k0: 1e-4, alpha: 0.5 };
  const iron = { ...half('Fe3+ + e- = Fe2+'), k0: 1e-3, alpha: 0.5 };
  const phiBath = (lib.materials.water.species['Cl-'].mu0 + RT * Math.log(c['Cl-'] / 1000)) / FARADAY;
  const VFe = (lib.materials.water.species['Fe3+'].mu0 - lib.materials.water.species['Fe2+'].mu0) / FARADAY + phiBath;
  const make = (rxs, withC0, metalLeft) => {
    const w = layer('water', 10e-6, withC0 ? { c0: c } : {});
    const elec = [ohmic(VFe, ['e-']), layer('Pt', 1e-6), { reactions: rxs }];
    const stack = metalLeft ? [...elec, w, bath(c, 'Cl-')] : [bath(c, 'Cl-'), w, ...elec.slice().reverse()];
    return build({ library: [lib, metal('Pt')], stack, grid: { minCells: 12 } });
  };
  const cases = [
    ['[volmer, iron] c0 metalR', make([volmer, iron], true, false)],
    ['[iron, volmer] c0 metalR', make([iron, volmer], true, false)],
    ['[volmer, iron] no-c0 metalR', make([volmer, iron], false, false)],
    ['[volmer, iron] c0 metalL', make([volmer, iron], true, true)],
  ];
  for (const [label, def] of cases) {
    const d = new Device(def);
    const cold = d.solution();
    const wNode = cold.region.indexOf(label.includes('metalL') ? 1 : 0) + 2;
    assert.ok(Math.abs(cold.phi[wNode] - phiBath) < 1e-6, `${label}: cold φ = ${cold.phi[wNode]} V vs bath ${phiBath} V`);
    const s = d.solve();
    assert.ok(s.converged && s.iterations <= 20, `${label}: converged=${s.converged} iters=${s.iterations}`);
  }
});

test('held insertion hosts stay at their contacts (two-host battery) and unheld hosts shift their collector', () => {
  // Removed Known item: a two-host battery with both contacts held started 137 mV off; a host behind
  // a floating Cu collector at I = 0 started Cu at 0 V (4.09 V from OCV).
  const x = [0.05, 0.2, 0.5, 0.8, 0.95];
  const ocvOf = (E0) => ({ x, E: x.map((v) => E0 - VT * Math.log(v / (1 - v))), muRef: 0 });
  const hostMat = (E0) => ({
    epsr: 0,
    species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } },
    statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: ocvOf(E0) }],
  });
  const elMat = {
    epsr: 0,
    species: { 'Li+': { D: 1e-10, mu0: -293e3, cRef: 1000 }, 'PF6-': { D: 2e-10, mu0: 0, cRef: 1000 } },
  };
  const species = [{ name: 'Li+', z: 1 }, { name: 'PF6-', z: -1 }, { name: 'e-', z: -1 }];
  const passLi = { phi: 'neutral', species: { 'Li+': 'equilibrium' } };
  const passE = { phi: 'neutral', species: { 'e-': 'equilibrium' } };
  const E_an = 0.1, E_cat = 3.8091;
  const ocvCell = E_cat - E_an;
  // 1. Two-host battery with both contacts held at OCV.
  const bat = new Device({
    T: T0,
    species,
    materials: { anode: hostMat(E_an), el: elMat, cathode: hostMat(E_cat) },
    regions: [
      { material: 'anode', length: 5e-6, c0: { 'Li+': 15000 } },
      { material: 'el', length: 20e-6, c0: { 'Li+': 1000, 'PF6-': 1000 } },
      { material: 'cathode', length: 5e-6, c0: { 'Li+': 15000 } },
    ],
    interfaces: [passLi, passLi],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { V: ocvCell, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
    },
    grid: { minCells: 8 },
  });
  const bc = bat.solution();
  assert.ok(Math.abs(bc.V['e-'][0] - 0) < 1e-9, `anode V_e = ${bc.V['e-'][0]}`);
  assert.ok(Math.abs(bc.V['e-'].at(-1) - ocvCell) < 1e-9, `cathode V_e = ${bc.V['e-'].at(-1)}`);
  const bs = bat.solve();
  assert.ok(bs.converged && bs.iterations <= 10, `battery iters = ${bs.iterations}`);

  // 2. Cathode behind a floating Cu collector at I = 0: collector starts at cell OCV.
  const batFloat = new Device({
    T: T0,
    species,
    materials: { anode: hostMat(E_an), el: elMat, cathode: hostMat(E_cat), Cu: { conductor: { species: 'e-', conductivity: 6e7 } } },
    regions: [
      { material: 'anode', length: 5e-6, c0: { 'Li+': 15000 } },
      { material: 'el', length: 20e-6, c0: { 'Li+': 1000, 'PF6-': 1000 } },
      { material: 'cathode', length: 5e-6, c0: { 'Li+': 15000 } },
      { material: 'Cu', length: 2e-6 },
    ],
    interfaces: [passLi, passLi, passE],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
    },
    grid: { minCells: 8 },
  });
  const fc = batFloat.solution();
  assert.ok(Math.abs(fc.terminals.right.V - ocvCell) < 1e-6, `floating Cu cold V = ${fc.terminals.right.V} V vs OCV ${ocvCell} V`);
  const fs = batFloat.solve();
  assert.ok(fs.converged && Math.abs(fs.terminals.right.V - ocvCell) < 1e-6);
});

test('electrode reaction with a neutral level-held participant balances φ_cold exactly', () => {
  // Removed Known item: previously balancedAt returned NaN when a neutral participant was mode 1,
  // leaving the electrolyte 0.68 V off its Nernst potential.
  const dev = new Device({
    T: T0,
    species: [
      { name: 'Ox+', z: 1, cRef: 1000 },
      { name: 'Red', z: 0, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
      { name: 'e-', z: -1 },
    ],
    materials: {
      Pt: { conductor: { species: 'e-', conductivity: 9e6 } },
      el: {
        epsr: 0,
        species: {
          'Ox+': { D: 1e-9, mu0: 50e3 },
          Red: { D: 1e-9, mu0: -15e3 },
          'Cl-': { D: 2e-9, mu0: 0 },
        },
      },
    },
    regions: [
      { material: 'Pt', length: 1e-6 },
      { material: 'el', length: 20e-6, c0: { 'Ox+': 10, 'Cl-': 10 } },
    ],
    interfaces: [{ reactions: [{ equation: 'Ox+ + e- = Red', k0: 1e-3, alpha: 0.5 }] }],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
      right: { phi: 'neutral', species: { Red: { type: 'equilibrium', mu: -15e3 + RT * Math.log(100 / 1000) } } },
    },
    grid: { minCells: 10 },
  });
  // Equilibrium: μ̄(Ox+) + μ̄(e-) = μ̄(Red) => 50e3 + RT ln(10/1000) + F φ + 0 = -15e3 + RT ln(100/1000)
  const phiEq = (-65e3 + RT * Math.log(10)) / FARADAY;
  const cold = dev.solution();
  const gEl = cold.region.indexOf(1);
  assert.ok(Math.abs(cold.phi[gEl] - phiEq) < 1e-9, `cold φ = ${cold.phi[gEl]} V, expected ${phiEq} V`);
  const s = dev.solve();
  assert.ok(s.converged && s.iterations <= 3, `iters = ${s.iterations}`);
});

test('interior metal whose only reaction is on its right face anchors to that face symmetrically', () => {
  const lib = aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+'], { epsr: 0 });
  const c = { 'K+': 100, 'Cl-': 125, 'Fe3+': 5, 'Fe2+': 5 };
  const def = build({
    library: [lib, metal('Pt')],
    stack: [
      { phi: 'neutral' },
      layer('water', 10e-6, { name: 'far', c0: c }),
      { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 },
      layer('Pt', 2e-6, { name: 'mid' }),
      { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1e-3, alpha: 0.5 }] },
      layer('water', 10e-6, { name: 'near', c0: c }),
      bath(c, 'Cl-'),
    ],
    grid: { minCells: 10 },
  });
  const dA = new Device(def), dB = new Device(mirrorDef(def));
  const cA = dA.solution(), cB = dB.solution();
  const w = lib.materials.water.species;
  const phiBath = (w['Cl-'].mu0 + RT * Math.log(c['Cl-'] / 1000)) / FARADAY;
  const VFe = (w['Fe3+'].mu0 - w['Fe2+'].mu0) / FARADAY + phiBath;
  const gPtA = cA.region.indexOf(1), gPtB = cB.region.indexOf(1);
  assert.ok(Math.abs(cA.V['e-'][gPtA] - VFe) < 1e-6, `right-face Pt cold V_e = ${cA.V['e-'][gPtA]} V vs ${VFe} V`);
  assert.ok(Math.abs(cA.V['e-'][gPtA] - cB.V['e-'][gPtB]) < 1e-12, 'mirror V_e');
  assert.ok(Math.abs(cA.phi[0] - VFe) < 1e-6, `far water behind blocking face φ = ${cA.phi[0]} V vs ${VFe} V`);
  assert.ok(Math.abs(cA.phi[0] - cB.phi.at(-1)) < 1e-12, 'mirror far φ');
});

test('strictly neutral divalent salt left of a high-voltage gated region solves cold (neutral-family n14/n15)', () => {
  // Previously left-to-right sweep put the left ungated salt at φ = 0 instead of propagating the
  // right gate's ~10.9 V anchor back across the neutral face, failing to converge cold.
  const def = {
    T: T0,
    species: [
      { name: 'Ca2+', z: 2, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: {
      el: {
        epsr: 0,
        species: {
          'Ca2+': { D: 7.9e-10, mu0: -553.5e3 },
          'Cl-': { D: 2.0e-9, mu0: -131.2e3 },
        },
      },
    },
    regions: [
      { name: 'leftSalt', material: 'el', length: 20e-6, c0: { 'Ca2+': 50, 'Cl-': 100 } },
      { name: 'gated', material: 'el', length: 20e-6, c0: { 'Ca2+': 50, 'Cl-': 100 } },
    ],
    interfaces: [{ phi: 'neutral', species: { 'Ca2+': 'equilibrium', 'Cl-': 'equilibrium' } }],
    contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [
      { name: 'g0', region: 'gated', V: 10.9, area: 1e4, capacitance: { C: 0.2, zeroCharge: 0 } },
    ],
    grid: { minCells: 10 },
  };
  const d = new Device(def);
  const cold = d.solution();
  assert.ok(Math.abs(cold.phi[0] - 10.9) < 1e-6, `left ungated salt cold φ = ${cold.phi[0]} V, expected 10.9 V`);
  const s = d.solve();
  assert.ok(s.converged && Math.abs(s.phi[0] - 10.9) < 1e-6, `solved φ = ${s.phi[0]} V`);
  assertColdMirrors(def, 'divalent salt | 10.9 V gated');
});

test('stretch source precedence: held over floating, contact over port, bulk over non-bulk', () => {
  const lib = aqueous(['Na+', 'Cl-'], { epsr: 0 });
  const c = { 'Na+': 100, 'Cl-': 100 };
  const phiBathAt = (V) => V + (lib.materials.water.species['Cl-'].mu0 + RT * Math.log(100 / 1000)) / FARADAY;
  // 1. Floating left bath (I: 0) vs held right bath at V = 2.0 V: must anchor to held right.
  const dFloatHeld = new Device(build({
    library: [lib],
    stack: [bath(c, 'Cl-', { I: 0 }), layer('water', 10e-6), bath(c, 'Cl-', 2.0)],
    grid: { minCells: 10 },
  }));
  const cf = dFloatHeld.solution();
  assert.ok(Math.abs(cf.phi[0] - phiBathAt(2.0)) < 1e-9, `held right φ = ${cf.phi[0]} vs ${phiBathAt(2.0)}`);
  assert.ok(Math.abs(cf.terminals.left.V - 2.0) < 1e-9, `floating left V = ${cf.terminals.left.V} vs 2.0`);

  // 2. Non-bulk held left contact (conductance link) vs bulk held right bath at V = 2.0 V:
  //    'bulk' contact wins over non-'bulk' regardless of side (old anchorSide picked left!).
  const w = lib.materials.water.species;
  const offNa = (w['Na+'].mu0 + RT * Math.log(100 / 1000)) / FARADAY + phiBathAt(0);
  const offCl = -(w['Cl-'].mu0 + RT * Math.log(100 / 1000)) / FARADAY + phiBathAt(0);
  const dNonBulkBulk = new Device(build({
    library: [lib],
    stack: [
      { V: 0, phi: 'neutral', species: { 'Na+': { type: 'conductance', G: 1e-3, offset: offNa }, 'Cl-': { type: 'conductance', G: 1e-3, offset: offCl } } },
      layer('water', 10e-6),
      bath(c, 'Cl-', 2.0),
    ],
    grid: { minCells: 10 },
  }));
  const cnb = dNonBulkBulk.solution();
  assert.ok(Math.abs(cnb.phi[0] - phiBathAt(2.0)) < 1e-9, `bulk right φ = ${cnb.phi[0]} vs ${phiBathAt(2.0)}`);
});
