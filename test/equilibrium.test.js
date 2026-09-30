import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, GAS_CONSTANT } from '../src/index.js';

const T = 298.15;
const RT = GAS_CONSTANT * T;
const VT = RT / FARADAY;

const salt = () => [
  { name: 'Na+', z: 1, cRef: 1000 },
  { name: 'Cl-', z: -1, cRef: 1000 },
];
const water = { epsr: 78.5, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } };
const gate = (V, C = 0.2) => ({ phi: { type: 'capacitive', C, zeroCharge: 0, V } });

// A floating island of 1:1 salt between two gates. Every species is a conserved spectator.
function island({ c = 10, L = 300e-9, VL = 0, VR = 0, C = 0.2, hmin = 0.02e-9 } = {}) {
  return {
    species: salt(),
    materials: { water },
    regions: [{ material: 'water', length: L, c0: { 'Na+': c, 'Cl-': c } }],
    contacts: { left: gate(VL, C), right: gate(VR, C) },
    grid: { hmin, hmax: 5e-9, ratio: 1.1 },
  };
}

const debyeLength = (epsr, c) => Math.sqrt((epsr * EPS0 * RT) / (2 * FARADAY * FARADAY * c));
const spread = (a) => Math.max(...a) - Math.min(...a);
const mid = (a) => a[a.length >> 1];

function assertEquilibrium(sol, tolRT = 1e-9) {
  assert.ok(sol.converged, 'converged');
  for (const [name, mu] of Object.entries(sol.mu)) {
    const finite = mu.filter(Number.isFinite);
    assert.ok(spread(finite) / RT < tolRT, `${name}: μ̄ spread ${spread(finite) / RT} RT`);
  }
}

// Total charge (device + gates) vanishes, relative to the gross ionic charge present.
function assertGauss(sol) {
  const gates = (sol.gates.left?.charge ?? 0) + (sol.gates.right?.charge ?? 0);
  const gross = FARADAY * sol.conservation.reduce((s, st) => s + st.amount, 0);
  assert.ok(Math.abs(sol.charge + gates) < 1e-12 * gross, `Gauss: device ${sol.charge}, gates ${gates}`);
}

test('floating island: flat μ̄, Gauss, and exact conservation at equilibrium', () => {
  const sol = new Device(island({ VL: 0.1 })).solve();
  assertEquilibrium(sol);
  assertGauss(sol);
  for (const st of sol.conservation) {
    assert.ok(st.spectator);
    assert.ok(Math.abs(st.drift) < 1e-13, `${st.species} drift ${st.drift}`);
  }
});

test('floating island: amounts stay conserved through a warm-started gate sweep', () => {
  const dev = new Device(island());
  dev.solve();
  for (let k = 0; k <= 20; k++) {
    const V = -0.5 + k * 0.05;
    dev.set({ contacts: { left: { phi: { V } }, right: { phi: { V: -V } } } });
    const sol = dev.solve();
    assertEquilibrium(sol);
    assertGauss(sol);
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `V=${V} ${st.species} drift ${st.drift}`);
  }
});

test('floating island, linear limit: gate charge matches series C_gate + C_Debye', () => {
  const c = 10, C = 0.2, dV = 1e-4;
  const sol = new Device(island({ c, C, VL: dV })).solve();
  // Neutral island ⇒ its bulk sits midway; each face is C_gate in series with ε/λ_D.
  const Cd = (78.5 * EPS0) / debyeLength(78.5, c);
  const Cs = 1 / (1 / C + 1 / Cd);
  const expected = (Cs * dV) / 2;
  const got = sol.gates.left.charge;
  assert.ok(Math.abs(got / expected - 1) < 1e-3, `gate charge ${got} vs ${expected}`);
  assert.ok(Math.abs(mid(sol.phi) - dV / 2) < 1e-3 * dV);
});

test('floating island, nonlinear: Gouy–Chapman charge and profile', () => {
  // Symmetric gates at ±V: the island bulk stays at φ = 0 and each face is a Gouy–Chapman layer.
  const V = 0.4, epsr = 78.5, eps = epsr * EPS0;
  const sol = new Device(island({ c: 10, VL: V, VR: -V, C: 0.3 })).solve();
  assertEquilibrium(sol);
  const cb = mid(sol.c['Na+']);
  assert.ok(Math.abs(mid(sol.c['Cl-']) / cb - 1) < 1e-9, 'bulk neutral');
  const lambda = debyeLength(epsr, cb);
  const zs = (sol.phi[0] - mid(sol.phi)) / VT;
  assert.ok(zs > 4, `strongly nonlinear surface potential (${zs} thermal units)`);

  // Diffuse-layer charge: σ = −√(8εRTc) sinh(ζ/2), balanced by the left gate's charge.
  const sigma = -Math.sqrt(8 * eps * RT * cb) * Math.sinh(zs / 2);
  assert.ok(Math.abs(-sol.gates.left.charge / sigma - 1) < 2e-3, `σ ${-sol.gates.left.charge} vs ${sigma}`);

  // Profile: tanh(ζ(x)/4) = tanh(ζs/4)·exp(−x/λ).
  let worst = 0;
  for (let g = 0; g < sol.x.length && sol.x[g] < 6 * lambda; g++) {
    const z = (sol.phi[g] - mid(sol.phi)) / VT;
    const zExact = 4 * Math.atanh(Math.tanh(zs / 4) * Math.exp(-sol.x[g] / lambda));
    worst = Math.max(worst, Math.abs(z - zExact));
  }
  assert.ok(worst < 5e-3, `worst profile error ${worst} thermal units`);

  // Boltzmann populations against the local potential, with flat μ̄.
  const g = 3;
  const psi = (sol.phi[g] - mid(sol.phi)) / VT;
  assert.ok(Math.abs(sol.c['Cl-'][g] / (cb * Math.exp(psi)) - 1) < 1e-9);
  assert.ok(Math.abs(sol.c['Na+'][g] / (cb * Math.exp(-psi)) - 1) < 1e-9);
});

test('A | B | A stack with unequal alignments: each face is its own double layer', () => {
  // Two solvents; ions have different standard potentials in each (transfer energies), and
  // the two A|B faces carry different dipoles, as with different adsorbates on each.
  const dmuNa = 3 * RT, dmuCl = -1 * RT; // μ°(B) − μ°(A)
  const B = {
    epsr: 32.7,
    species: {
      'Na+': { D: 1e-9, mu0: water.species['Na+'].mu0 + dmuNa },
      'Cl-': { D: 1e-9, mu0: water.species['Cl-'].mu0 + dmuCl },
    },
  };
  const d1 = 0.05, d2 = -0.02; // interface dipoles, V (these need not cancel)
  const L = 400e-9;
  const def = {
    species: salt(),
    materials: { A: water, B },
    regions: [
      { material: 'A', length: L, c0: { 'Na+': 10, 'Cl-': 10 } },
      { material: 'B', length: L, c0: { 'Na+': 2, 'Cl-': 2 } },
      { material: 'A', length: L, c0: { 'Na+': 10, 'Cl-': 10 } },
    ],
    interfaces: [{ dipole: d1 }, { dipole: -d2 }],
    contacts: { left: gate(0), right: gate(0) },
    grid: { hmin: 0.02e-9, hmax: 5e-9, ratio: 1.1 },
  };
  const dev = new Device(def);
  const sol = dev.solve();
  assertEquilibrium(sol);
  assertGauss(sol);
  for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12);

  // Bulk values at the middle of each region.
  const g = dev.grid;
  const midNode = (r) => (g.regionStart[r] + g.regionEnd[r]) >> 1;
  const bulk = [0, 1, 2].map((r) => ({ phi: sol.phi[midNode(r)], c: sol.c['Na+'][midNode(r)] }));

  // The salt partition is gauge-free: c+c− ratio set by the transfer energies alone.
  const cB = bulk[1].c, cA = bulk[0].c;
  assert.ok(Math.abs(sol.c['Cl-'][midNode(1)] / cB - 1) < 1e-8, 'B bulk neutral');
  assert.ok(Math.abs((cB * cB) / (cA * cA) / Math.exp(-(dmuNa + dmuCl) / RT) - 1) < 1e-8);

  // Each face: A's diffuse drop a plus B's diffuse drop b plus the dipole adds up to the
  // bulk-to-bulk step, and the two diffuse charges cancel:
  //   √(ε_A c_A) sinh(a/2) = √(ε_B c_B) sinh(b/2)
  const faces = [
    { f: 0, rA: 0, rB: 1, sign: 1, dip: d1 },
    { f: 1, rA: 2, rB: 1, sign: -1, dip: d2 },
  ];
  for (const { f, rA, rB, sign, dip } of faces) {
    // Walk from A's bulk toward B's bulk; for the right face that is right-to-left.
    const nA = sign > 0 ? g.regionEnd[rA - 0] : g.regionStart[rA];
    const nB = sign > 0 ? g.regionStart[rB] : g.regionEnd[rB];
    const total = (bulk[rB].phi - bulk[rA].phi) / VT;
    const a = (sol.phi[nA] - bulk[rA].phi) / VT;
    const b = (bulk[rB].phi - sol.phi[nB]) / VT;
    const dipHat = (sol.phi[nB] - sol.phi[nA]) / VT;
    assert.ok(Math.abs(dipHat - dip / VT) < 1e-9, `face ${f} dipole`);
    assert.ok(Math.abs(a + dipHat + b - total) < 1e-9);

    // Analytic split by bisection on a.
    const kA = Math.sqrt(78.5 * cA), kB = Math.sqrt(32.7 * cB);
    const rest = total - dip / VT;
    let lo = -50, hi = 50;
    for (let it = 0; it < 200; it++) {
      const m = (lo + hi) / 2;
      const h = kA * Math.sinh(m / 2) - kB * Math.sinh((rest - m) / 2);
      if (h > 0) hi = m;
      else lo = m;
    }
    const aExact = (lo + hi) / 2;
    assert.ok(Math.abs(a - aExact) < 2e-3, `face ${f}: a ${a} vs ${aExact}`);
    const sigmaExact = -Math.sqrt(8 * 78.5 * EPS0 * RT * cA) * Math.sinh(aExact / 2);
    const sigma = sign * sol.interfaces[f].D; // charge in A's diffuse layer
    assert.ok(Math.abs(sigma / sigmaExact - 1) < 3e-3, `face ${f}: σ ${sigma} vs ${sigmaExact}`);
  }
});

test('ion-exchange membrane between two salt layers: Donnan equilibrium', () => {
  // The membrane carries fixed negative charge and excludes Cl⁻ entirely (absent species).
  const X = 50; // mol/m³ of fixed charge
  const M = { epsr: 20, species: { 'Na+': { D: 1e-10, mu0: water.species['Na+'].mu0 + 2 * RT } } };
  const L = 300e-9;
  const def = {
    species: salt(),
    materials: { water, M },
    regions: [
      { material: 'water', length: L, c0: { 'Na+': 10, 'Cl-': 10 } },
      { material: 'M', length: L, fixedCharge: -X * FARADAY, c0: { 'Na+': X } },
      { material: 'water', length: L, c0: { 'Na+': 10, 'Cl-': 10 } },
    ],
    interfaces: [{ dipole: 0 }, { dipole: 0 }],
    contacts: { left: gate(0), right: gate(0) },
    grid: { hmin: 0.02e-9, hmax: 5e-9, ratio: 1.1 },
  };
  const dev = new Device(def);
  const sol = dev.solve();
  assertEquilibrium(sol);
  assertGauss(sol);
  const g = dev.grid;
  const midNode = (r) => (g.regionStart[r] + g.regionEnd[r]) >> 1;
  // Cl⁻ is absent in the membrane: NaN there, and conserved separately on each side.
  assert.ok(sol.c['Cl-'].slice(g.regionStart[1], g.regionEnd[1] + 1).every(Number.isNaN));
  assert.equal(sol.conservation.filter((s) => s.species === 'Cl-').length, 2);
  // Membrane bulk is neutral: c_Na = X.
  assert.ok(Math.abs(sol.c['Na+'][midNode(1)] / X - 1) < 1e-8);
  // Donnan: c_Na(M)/c_Na(water) = exp(−Δμ°/RT − ΔφD/VT)
  const cW = sol.c['Na+'][midNode(0)];
  const dphi = (sol.phi[midNode(1)] - sol.phi[midNode(0)]) / VT;
  assert.ok(Math.abs(X / cW / Math.exp(-2 - dphi) - 1) < 1e-8);
});

test('equilibrium is reached from a far-off start and is independent of the path', () => {
  const a = new Device(island({ VL: 0.3, VR: -0.1 })).solve();
  const dev = new Device(island({ VL: -0.4, VR: 0.4 }));
  dev.solve();
  dev.set({ contacts: { left: { phi: { V: 0.3 } }, right: { phi: { V: -0.1 } } } });
  const b = dev.solve();
  for (let g = 0; g < a.phi.length; g++) assert.ok(Math.abs(a.phi[g] - b.phi[g]) < 1e-10);
});
