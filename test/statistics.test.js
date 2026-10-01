import { test } from 'node:test';
import assert from 'node:assert/strict';
import { AVOGADRO, Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { fermiHalf, fermiMinusHalf } from '../src/fermi.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// Independent Fermi–Dirac integrals 𝓕_j (normalised to e^x as x → −∞), by the trapezoid rule
// in u = √t, where the integrand is smooth and even.
const GAMMA = { 0.5: Math.sqrt(Math.PI) / 2, 1.5: (3 * Math.sqrt(Math.PI)) / 4, '-0.5': Math.sqrt(Math.PI) };
function fermi(j, x) {
  const h = 0.01, umax = Math.sqrt(Math.max(x, 0) + 45);
  const f = (u) => {
    const p = 2 * u ** (2 * j + 1), e = u * u - x;
    return e > 0 ? (p * Math.exp(-e)) / (1 + Math.exp(-e)) : p / (1 + Math.exp(e));
  };
  let s = f(0) / 2;
  for (let k = 1; k * h <= umax; k++) s += f(k * h);
  return (s * h) / GAMMA[j];
}

const throwsDevice = (def, pattern) =>
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && pattern.test(e.message), String(pattern));

test('Fermi–Dirac integrals 𝓕_{1/2} and 𝓕_{−1/2} to 1e-13 across every branch', () => {
  let worst = 0;
  for (let x = -40; x <= 90; x += 0.173) {
    worst = Math.max(worst, Math.abs(fermiHalf(x) / fermi(0.5, x) - 1), Math.abs(fermiMinusHalf(x) / fermi(-0.5, x) - 1));
  }
  assert.ok(worst < 1e-13, `worst relative error ${worst}`);
});

// Degenerately doped silicon with Fermi–Dirac electrons and holes.
const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19), Eg = units.eV(1.12);
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
const silicon = (ND, left) => ({
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: {
      epsr: 11.7,
      species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: Eg, cRef: Nv } },
      statistics: [{ type: 'fermi-dirac', species: ['e-', 'h+'] }],
    },
  },
  regions: [{ material: 'Si', length: 200e-9, fixedCharge: ND * FARADAY }],
  contacts: { left, right: ohmic(0) },
  grid: { hmin: 0.01e-9, hmax: 2e-9, ratio: 1.06 },
});

test('Fermi–Dirac: degenerate bulk n = N_c 𝓕_{1/2}(ζ), neutral, with flat μ̄', () => {
  const ND = units.perCm3(1e20);
  const sol = new Device(silicon(ND, ohmic(0))).solve();
  assert.ok(sol.converged);
  const m = sol.x.length >> 1;
  const zeta = (sol.mu['e-'][m] + FARADAY * sol.phi[m]) / RT; // (μ̄ − μ° − zFφ)/RT, z = −1
  assert.ok(zeta > 2, `degenerate (ζ = ${zeta})`);
  assert.ok(Math.abs(sol.c['e-'][m] / (Nc * fermi(0.5, zeta)) - 1) < 1e-12);
  // Two-dimensional density of states: 𝓕₀(ζ) = ln(1 + e^ζ).
  const bar = neutralBar([{ type: 'fermi-dirac', species: ['X'], order: 0 }], 1.3, 1.3, 10).solve();
  assert.ok(Math.abs(bar.c.X[4] / (1000 * Math.log1p(Math.exp(1.3))) - 1) < 1e-14);
  assert.ok(Math.abs((sol.c['e-'][m] - sol.c['h+'][m]) / ND - 1) < 1e-12);
  const mu = sol.mu['e-'];
  assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-12);
});

test('Fermi–Dirac: accumulation-layer charge from Gauss with 𝓕_{3/2} (degenerate surface)', () => {
  // Electrons in equilibrium with a metal whose pinned φ sets the surface ζ_s = −zeroCharge/V_T.
  const ND = units.perCm3(1e18);
  for (const zeroCharge of [-0.1, -0.25]) {
    const left = { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'pinned', zeroCharge };
    const dev = new Device(silicon(ND, left));
    const sol = dev.solve();
    assert.ok(sol.converged);
    const m = sol.x.length >> 1;
    let Q = 0;
    for (let k = 0; k < sol.x.length; k++) Q += dev.grid.vol[k] * FARADAY * (ND - sol.c['e-'][k] + sol.c['h+'][k]);
    const zb = (sol.mu['e-'][m] + FARADAY * sol.phi[m]) / RT;
    const zpb = (sol.mu['h+'][m] - Eg - FARADAY * sol.phi[m]) / RT;
    const psi = (sol.phi[0] - sol.phi[m]) / VT;
    assert.ok(Math.abs(zb + psi + zeroCharge / VT) < 1e-9, 'surface ζ');
    // ½ ε E_s² = −∫ρ dφ from bulk to surface: n and p integrate to N 𝓕_{3/2}.
    const I = ND * psi - Nc * (fermi(1.5, zb + psi) - fermi(1.5, zb)) - Nv * (fermi(1.5, zpb - psi) - fermi(1.5, zpb));
    const want = -Math.sign(psi) * Math.sqrt(2 * 11.7 * EPS0 * RT * Math.abs(I));
    assert.ok(Math.abs(Q / want - 1) < 1e-3, `zeroCharge ${zeroCharge}: ${Q} vs ${want}`);
  }
});

// A 1:1 electrolyte on a lattice of c_max sites shared by both ions (Bikerman).
const salt = [
  { name: 'Na+', z: 1, cRef: 1000 },
  { name: 'Cl-', z: -1, cRef: 1000 },
];
const brine = (statistics, epsr = 78.5) => ({
  epsr,
  species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } },
  statistics,
});
const bath = (c) => ({ bath: { c: { 'Na+': c, 'Cl-': c }, reference: 'Cl-' } });

const crowded = (statistics, c0, lam) => new Device(crowdedDef(statistics, c0, lam));
function crowdedDef(statistics, c0, lam) {
  return {
    species: salt,
    materials: { water: brine(statistics) },
    regions: [{ material: 'water', length: 40 * lam }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0 }, right: { V: 0, ...bath(c0) } },
    grid: { hmin: lam / 400, hmax: lam / 4, ratio: 1.05 },
  };
}

test('lattice gas: the Kilic–Bazant–Ajdari crowded double layer, charge ∝ √ψ at large ψ', () => {
  const c0 = 100, nu = 0.1, cMax = (2 * c0) / nu;
  const lam = Math.sqrt((78.5 * EPS0 * RT) / (2 * FARADAY * FARADAY * c0));
  const dev = crowded([{ type: 'lattice', species: ['Na+', 'Cl-'], cMax }], c0, lam);
  const phiBulk = dev.solve().phi.at(-1);
  for (const psiWant of [2, 10, 40, -20]) {
    // A blocking electrode pinned at φ_edge = V − zeroCharge, ψ thermal units above the bulk.
    dev.set({ contacts: { left: { V: phiBulk + psiWant * VT } } });
    const sol = dev.solve();
    assert.ok(sol.converged);
    let Q = 0;
    for (let k = 0; k < sol.x.length; k++) Q += dev.grid.vol[k] * FARADAY * (sol.c['Na+'][k] - sol.c['Cl-'][k]);
    const psi = (sol.phi[0] - sol.phi.at(-1)) / VT;
    assert.ok(Math.abs(psi - psiWant) < 1e-9);
    const want = -Math.sign(psi) * 2 * FARADAY * c0 * lam * Math.sqrt((2 / nu) * Math.log(1 + 2 * nu * Math.sinh(psi / 2) ** 2));
    assert.ok(Math.abs(Q / want - 1) < 5e-4, `ψ = ${psi}: ${Q} vs ${want}`);
    assert.ok(Math.max(...sol.c['Na+'], ...sol.c['Cl-']) <= cMax, 'never beyond the site density');
    assert.ok(Math.abs(sol.c['Na+'].at(-1) / c0 - 1) < 1e-12, 'bath composition reproduced');
  }
});

test('custom statistics: a user function reproduces the built-in lattice gas', () => {
  const c0 = 100, cMax = 2000;
  const lam = Math.sqrt((78.5 * EPS0 * RT) / (2 * FARADAY * FARADAY * c0));
  const lattice = (zeta) => {
    const w = zeta.map((z) => (1000 / cMax) * Math.exp(z));
    const S = 1 + w[0] + w[1];
    const c = w.map((wi) => (cMax * wi) / S);
    return { c, dcdzeta: c.map((ca, a) => c.map((cb, b) => (a === b ? ca : 0) - (ca * cb) / cMax)) };
  };
  const a = crowded([{ type: 'lattice', species: ['Na+', 'Cl-'], cMax }], c0, lam);
  const b = crowded([{ type: 'custom', species: ['Na+', 'Cl-'], evaluate: lattice }], c0, lam);
  for (const V of [-1.2, -1.5]) {
    a.set({ contacts: { left: { V } } });
    b.set({ contacts: { left: { V } } });
    const sa = a.solve(), sb = b.solve();
    assert.ok(sb.converged);
    assert.ok(Math.abs(sb.contacts.left.D / sa.contacts.left.D - 1) < 1e-10);
  }
  const bad = crowdedDef([{ type: 'custom', species: ['Na+', 'Cl-'], evaluate: () => ({ c: [1, 1], dcdzeta: [[1, 0.5], [0, 1]] }) }], c0, lam);
  throwsDevice(bad, /must be symmetric/);
  bad.materials.water.statistics[0].evaluate = () => ({ c: [1, 1], dcdzeta: [[1, 2], [2, 1]] });
  throwsDevice(bad, /positive definite/);
});

test('lattice gas: a floating island conserves its amounts exactly through a gate sweep', () => {
  const gate = (V) => ({ V, phi: { type: 'capacitive', C: 0.3 }, zeroCharge: 0 });
  const dev = new Device({
    species: salt,
    materials: { water: brine([{ type: 'lattice', species: ['Na+', 'Cl-'], cMax: 400 }]) },
    regions: [{ material: 'water', length: 100e-9, c0: { 'Na+': 100, 'Cl-': 100 } }],
    contacts: { left: gate(0), right: gate(0) },
    grid: { hmin: 0.02e-9, hmax: 2e-9, ratio: 1.1 },
  });
  for (const V of [0, 0.3, 1, -1]) {
    dev.set({ contacts: { left: { V }, right: { V: -V } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `V=${V}`);
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `V=${V} ${st.species}: ${st.drift}`);
    for (const mu of Object.values(sol.mu)) assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-9);
  }
});

// A neutral species diffusing between two reservoirs. With the flux −(D c/RT)∇μ, the steady
// flux is exactly N = −(D/L)[P(ζ_R) − P(ζ_L)], where P = ∫ c dζ is the grand potential density.
const cMax = 1000;
const neutralBar = (statistics, zL, zR, cells) => new Device(neutralBarDef(statistics, zL, zR, cells));
const neutralBarDef = (statistics, zL, zR, cells) => ({
    species: [{ name: 'X', z: 0, cRef: 1000 }],
    materials: { m: { epsr: 10, species: { X: { D: 1e-9, mu0: 0 } }, statistics } },
    regions: [{ material: 'm', length: 1e-6 }],
    contacts: {
      left: { species: { X: { type: 'equilibrium', mu: zL * RT } }, phi: 'neutral' },
      right: { species: { X: { type: 'equilibrium', mu: zR * RT } }, phi: 'neutral' },
    },
    grid: { hmin: 1e-6 / cells, hmax: 1e-6 / cells },
});
const simpson = (f, a, b, K = 20000) => {
  const h = (b - a) / K;
  let s = f(a) + f(b);
  for (let k = 1; k < K; k++) s += (k % 2 ? 4 : 2) * f(a + k * h);
  return (s * h) / 3;
};

test('steady transport through non-ideal statistics: N = −(D/L)ΔP (constant mobility)', () => {
  // Langmuir: P = c_max ln(1 + (c_ref/c_max) e^ζ). The excess is exactly linear here, so any
  // grid is exact.
  const P = (z) => cMax * Math.log1p(Math.exp(z));
  const sol = neutralBar([{ type: 'lattice', species: ['X'], cMax }], 3, -2, 10).solve();
  assert.ok(sol.converged);
  const want = -(1e-9 / 1e-6) * (P(-2) - P(3));
  assert.ok(Math.abs(sol.contacts.left.flux.X / want - 1) < 1e-12);

  // Regular solution (Redlich–Kister A₀ = Ω): second-order convergence to the exact flux.
  const Om = 1.6 * RT;
  const zeta = (x) => Math.log(x / (1 - x)) - (2 * Om * x) / RT;
  const xOf = (z) => {
    let lo = 1e-15, hi = 1 - 1e-15;
    for (let k = 0; k < 200; k++) (zeta(0.5 * (lo + hi)) < z ? (lo = 0.5 * (lo + hi)) : (hi = 0.5 * (lo + hi)));
    return 0.5 * (lo + hi);
  };
  // ΔP = c_max ∫ x ζ′(x) dx
  const exact = -(1e-9 / 1e-6) * cMax * simpson((x) => x * (1 / x + 1 / (1 - x) - (2 * Om) / RT), xOf(2), xOf(-2));
  const err = [40, 160].map((cells) => {
    const s = neutralBar([{ type: 'redlich-kister', species: ['X'], cMax, A: [Om] }], 2, -2, cells).solve();
    assert.ok(s.converged);
    assert.ok(Math.abs(s.c.X[0] / cMax - xOf(2)) < 1e-12, 'boundary occupancy from the isotherm');
    return Math.abs(s.contacts.left.flux.X / exact - 1);
  });
  assert.ok(err[1] < 1e-4 && err[0] / err[1] > 12, `errors ${err}`);
});

test('Redlich–Kister: a non-convex free energy (miscibility gap) is refused', () => {
  throwsDevice(neutralBarDef([{ type: 'redlich-kister', species: ['X'], cMax, A: [2.2 * RT] }], 0, 0, 10), /not convex/);
  assert.doesNotThrow(() => neutralBar([{ type: 'redlich-kister', species: ['X'], cMax, A: [1.9 * RT] }], 0, 0, 10)); // Ω < 2RT
});

test('Debye–Hückel: liquid-junction EMF 2t₊(RT/F) ln(a₁/a₂) with activities, and warnings', () => {
  const Dp = 1.33e-9, Dm = 2.03e-9, c1 = 100, c2 = 10, tp = Dp / (Dp + Dm), a = 0.4e-9;
  const eps = 78.5 * EPS0, lB = FARADAY ** 2 / (4 * Math.PI * eps * AVOGADRO * RT);
  const lnGamma = (c) => {
    const k = Math.sqrt((2 * FARADAY ** 2 * c) / (eps * RT));
    return (-lB * k) / (2 * (1 + k * a));
  };
  const def = (epsr) => ({
    species: salt,
    materials: { water: brine([{ type: 'debye-huckel', species: ['Na+', 'Cl-'], epsr: 78.5, a }], epsr) },
    regions: [{ material: 'water', length: 10e-6 }],
    contacts: { left: bath(c1), right: { ...bath(c2), I: 0 } },
    grid: { minCells: 1600 },
  });
  const sol = new Device(def(0)).solve();
  assert.ok(sol.converged);
  const emf = 2 * tp * VT * (Math.log(c1 / c2) + lnGamma(c1) - lnGamma(c2));
  assert.ok(Math.abs(sol.terminalVoltage / emf - 1) < 1e-6, `${sol.terminalVoltage} vs ${emf}`);
  // The bath's salt chemical potential includes the mean activity coefficient.
  const salt0 = (sol.mu['Na+'][0] + sol.mu['Cl-'][0] - (-261.9e3 - 131.2e3)) / RT;
  assert.ok(Math.abs(salt0 - 2 * (Math.log(c1 / 1000) + lnGamma(c1))) < 1e-12);
  assert.deepEqual(sol.warnings, []);
  const warned = new Device(def(78.5));
  assert.ok(warned.model.warnings.some((w) => /Debye–Hückel/.test(w)));
});

// An intercalation host (ε = 0) holding Li⁺ + e⁻, described by its OCV against Li metal.
const cHost = 20000, Om = 1.2 * RT, mu0Li = -293.3e3;
const E = (x) => -(mu0Li + RT * Math.log(x / (1 - x)) - 2 * Om * x) / FARADAY; // regular solution, μ_Li(metal) = 0
const host = (stat, left, right, cells = 50) => ({
  species: [
    { name: 'Li+', z: 1 },
    { name: 'e-', z: -1 },
  ],
  materials: {
    host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: mu0Li, cRef: cHost }, 'e-': { D: 1e-4, mu0: 0, cRef: cHost } }, statistics: [stat] },
  },
  regions: [{ material: 'host', length: 1e-6 }],
  contacts: { left, right },
  grid: { minCells: cells },
});
const liIon = (V) => ({ V, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' });
const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
const regular = { type: 'insertion', species: ['Li+', 'e-'], cMax: cHost, A: [Om] };

test('insertion host: composition follows the OCV of the neutral combination (table round trip)', () => {
  // Li⁺ from a Li reference on the left, e⁻ from a current collector at V on the right. At
  // equilibrium μ_Li = μ̄_Li⁺ + μ̄_e⁻ = −F V, so the host sits at the x where E(x) = V.
  const xs = [0.05, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9];
  const table = { type: 'insertion', species: ['Li+', 'e-'], cMax: cHost, ocv: { x: xs, E: xs.map(E), muRef: 0 } };
  const tab = new Device(host(table, liIon(0), collector(E(0.5))));
  const rk = new Device(host(regular, liIon(0), collector(E(0.5))));
  for (const x of [0.01, 0.1, 0.3, 0.52, 0.9, 0.99]) {
    tab.set({ contacts: { right: { V: E(x) } } });
    rk.set({ contacts: { right: { V: E(x) } } });
    const st = tab.solve(), sr = rk.solve();
    assert.ok(st.converged && sr.converged);
    assert.ok(Math.abs(sr.c['Li+'][10] / cHost - x) < 1e-12, `analytic isotherm at x=${x}`);
    // The table's non-ideal residual is interpolated, so a regular-solution table is exact
    // between its points and beyond its ends.
    assert.ok(Math.abs(st.c['Li+'][10] / cHost - x) < 1e-12, `table at x=${x}: ${st.c['Li+'][10] / cHost}`);
    assert.ok(Number.isNaN(st.phi[10]), 'φ is undefined inside the host');
    assert.ok(Math.abs(st.c['e-'][10] / st.c['Li+'][10] - 1) < 1e-14, 'neutral');
  }
});

test('insertion host: steady chemical diffusion N = −(D_Li/L)ΔP with fast electrons', () => {
  // Left: a Li reservoir (Li⁺ and e⁻) at x = 0.8. Right: Li⁺ only, at the level of x = 0.2.
  // Electrons can't leave on the right, so μ̄_e⁻ is flat and ∇μ̄_Li⁺ = ∇μ_Li.
  const xL = 0.8, xR = 0.2;
  const left = { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'Li+': { type: 'equilibrium', offset: -E(xL) } }, phi: 'bulk' };
  const exact = -(1e-14 / 1e-6) * cHost * simpson((x) => x * (1 / x + 1 / (1 - x) - (2 * Om) / RT), xL, xR);
  const err = [50, 200].map((cells) => {
    const sol = new Device(host(regular, left, liIon(-E(xR)), cells)).solve();
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.current / (FARADAY * sol.contacts.right.flux['Li+']) - 1) < 1e-12);
    return Math.abs(sol.contacts.right.flux['Li+'] / exact - 1);
  });
  assert.ok(err[1] < 1e-5 && err[0] / err[1] > 12, `errors ${err}`);
});

test('insertion host: transient uptake relaxes at π²D_chem/4L², D_chem = D_Li·x ζ′(x), conserving exactly', () => {
  // A small step in the collector voltage; Li⁺ enters on the left only. At late times the
  // deficit decays by 1/(1 + λΔt) per backward-Euler step, with λ the slowest diffusion mode.
  const x0 = 0.8, L = 1e-6;
  const dev = new Device(host(regular, liIon(0), collector(E(x0)), 100));
  dev.solve();
  dev.set({ contacts: { right: { V: E(x0) - 1e-6 } } });
  const final = new Device(dev.def).solve();
  const amount = (s) => s.conservation.find((st) => st.species === 'Li+').amount;
  const xf = final.c['Li+'][5] / cHost;
  const lambda = (Math.PI ** 2 * 1e-14 * xf * (1 / xf + 1 / (1 - xf) - (2 * Om) / RT)) / (4 * L * L);
  const dt = 0.05 / lambda;
  let prev, rate;
  for (let k = 0; k < 60; k++) {
    const sol = dev.step(dt);
    for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-13, `${st.species} drift ${st.drift}`);
    const deficit = amount(sol) - amount(final);
    if (prev !== undefined) rate = (prev / deficit - 1) / dt;
    prev = deficit;
  }
  assert.ok(Math.abs(rate / lambda - 1) < 1e-4, `rate ${rate} vs ${lambda}`);
});

test('insertion statistics are checked', () => {
  let def = host(regular, liIon(0), collector(0));
  def.materials.host.epsr = 5;
  throwsDevice(def, /needs epsr: 0/);
  def = host(regular, liIon(0), { ...collector(0), phi: 'pinned', zeroCharge: 0 });
  throwsDevice(def, /φ is undefined/);
  def = host({ ...regular, A: undefined, ocv: { x: [0.2, 0.5, 0.8], E: [3.5, 3.6, 3.4], muRef: 0 } }, liIon(0), collector(0));
  throwsDevice(def, /strictly decreasing/);
  def = host(regular, liIon(0), collector(0));
  def.regions[0].fixedCharge = -10 * FARADAY;
  throwsDevice(def, /opposite sign to that carrier/);
  def = host({ ...regular, species: ['Li+', 'Li+'] }, liIon(0), collector(0));
  throwsDevice(def, /already belongs/);
});
