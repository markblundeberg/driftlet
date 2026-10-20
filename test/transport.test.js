import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// A neutral species between two reservoirs at c_L and c_R.
const D = 1e-9, L = 1e-5, cL = 2, cR = 1;
const bar = (region) =>
  new Device({
    species: [{ name: 'X', z: 0, cRef: 1000 }],
    materials: { m: { epsr: 10, species: { X: { D, mu0: 0 } } } },
    regions: [{ material: 'm', length: L, ...region }],
    contacts: {
      left: { species: { X: { type: 'equilibrium', mu: RT * Math.log(cL / 1000) } }, phi: 'neutral' },
      right: { species: { X: { type: 'equilibrium', mu: RT * Math.log(cR / 1000) } }, phi: 'neutral' },
    },
    grid: { hmin: L / 10, hmax: L / 10 },
  });

test('advection: exact convection–diffusion profile and flux on a coarse grid, either direction', () => {
  // N = v c − D c′ is constant: c = c_L + (c_R − c_L)(e^{vx/D} − 1)/(e^{Pe} − 1), Pe = vL/D.
  for (const v of [1e-4, -3e-4, 1e-6]) {
    const sol = bar({ velocity: v }).solve();
    assert.ok(sol.converged);
    const Pe = (v * L) / D;
    const want = (v * (cL * Math.exp(Pe) - cR)) / Math.expm1(Pe);
    assert.ok(Math.abs(sol.contacts.left.flux.X / want - 1) < 1e-12, `v=${v}`);
    assert.ok(Math.abs(sol.contacts.right.flux.X / want - 1) < 1e-12);
    for (let g = 0; g < sol.x.length; g++) {
      const exact = cL + ((cR - cL) * Math.expm1((v * sol.x[g]) / D)) / Math.expm1(Pe);
      assert.ok(Math.abs(sol.c.X[g] / exact - 1) < 1e-12);
    }
  }
});

test('advection far past a cell Péclet number of 700, either direction: every flux v·c upstream', () => {
  // The flux's weights B(Δ) and expm1(Δη) underflow and overflow past Δ ≈ 710 (flow toward −x):
  // the product was 0·∞, every flux NaN, and the solve failed (found by a cold review of
  // src/regions.js). Toward +x it was exact. With lattice statistics too (the concentrated kernel).
  const lattice = { type: 'lattice', species: ['X'], cMax: 3000 };
  const mu = (c, stats) => RT * Math.log(c / 1000) + (stats ? RT * Math.log(3000 / (3000 - c)) : 0);
  for (const stats of [null, lattice]) {
    for (const v of [1, -1, 0.72, -0.72, -50]) {
      const def = {
        species: [{ name: 'X', z: 0, cRef: 1000 }],
        materials: { m: { epsr: 10, species: { X: { D, mu0: 0 } }, ...(stats ? { statistics: [stats] } : {}) } },
        regions: [{ material: 'm', length: L, velocity: v }],
        contacts: {
          left: { species: { X: { type: 'equilibrium', mu: mu(cL, stats) } }, phi: 'neutral' },
          right: { species: { X: { type: 'equilibrium', mu: mu(cR, stats) } }, phi: 'neutral' },
        },
        grid: { hmin: L / 10, hmax: L / 10 },
      };
      const sol = new Device(def).solve(), want = v * (v > 0 ? cL : cR), tag = `${stats ? 'lattice' : 'ideal'}, v = ${v}`;
      assert.ok(sol.converged, tag);
      assert.ok(Math.abs(sol.contacts.left.flux.X / want - 1) < 1e-12, `${tag}: ${sol.contacts.left.flux.X}`);
      assert.ok(sol.flux.X.every((N) => Math.abs(N / want - 1) < 1e-12), `${tag}: ${sol.flux.X}`);
    }
  }
});

test('eddy mixing: a neutral species diffuses with D + D_mix, exactly', () => {
  for (const mixing of [1e-9, 1e-7]) {
    const sol = bar({ mixing }).solve();
    assert.ok(sol.converged && sol.iterations <= 8);
    assert.ok(Math.abs(sol.contacts.left.flux.X / (((D + mixing) * (cL - cR)) / L) - 1) < 1e-12);
  }
});

test('eddy mixing carries no current: junction EMF unchanged, salt flux raised by D_mix', () => {
  const Dp = 1.33e-9, Dm = 2.03e-9, tp = Dp / (Dp + Dm), c1 = 100, c2 = 10, Lj = 10e-6;
  const bath = (c) => ({ bath: { c: { 'Na+': c, 'Cl-': c }, reference: 'Cl-' } });
  const junction = (mixing, statistics) =>
    new Device({
      species: [
        { name: 'Na+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'Na+': { D: Dp, mu0: -261.9e3 }, 'Cl-': { D: Dm, mu0: -131.2e3 } }, ...(statistics ? { statistics } : {}) } },
      regions: [{ material: 'water', length: Lj, mixing }],
      contacts: { left: bath(c1), right: { ...bath(c2), I: 0 } },
      grid: { minCells: 400 },
    }).solve();
  const Ds = (2 * Dp * Dm) / (Dp + Dm), emf = 2 * tp * VT * Math.log(c1 / c2);
  for (const mixing of [0, 1e-8]) {
    const sol = junction(mixing);
    assert.ok(sol.converged);
    assert.ok(Math.abs(sol.terminalVoltage / emf - 1) < 1e-5, `EMF with D_mix = ${mixing}`);
    const flux = ((Ds + mixing) * (c1 - c2)) / Lj;
    assert.ok(Math.abs(sol.contacts.left.flux['Na+'] / flux - 1) < 1e-6, `salt flux with D_mix = ${mixing}`);
  }
  // Crowded, the EMF is mixing's to the same 1e-5 (mixing carries no current whatever K̄ is).
  const lattice = [{ type: 'lattice', species: ['Na+', 'Cl-'], cMax: 400 }];
  const [still, stirred] = [0, 1e-8].map((mixing) => junction(mixing, lattice));
  assert.ok(stirred.converged && Math.abs(stirred.terminalVoltage / still.terminalVoltage - 1) < 1e-5, `${stirred.terminalVoltage} vs ${still.terminalVoltage}`);
});

test('eddy mixing leaves equilibrium double layers alone', () => {
  // A gated island with strong mixing: still flat μ̄ and the same gate charge as without it.
  const gate = (V) => ({ V, phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 });
  // Crowded too (a lattice of 300 mM, nearly full at the gates), where K̄ isn't diagonal.
  const island = (mixing, crowded) =>
    new Device({
      species: [
        { name: 'Na+', z: 1, cRef: 1000 },
        { name: 'Cl-', z: -1, cRef: 1000 },
      ],
      materials: {
        water: {
          epsr: 78.5,
          species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } },
          ...(crowded ? { statistics: [{ type: 'lattice', species: ['Na+', 'Cl-'], cMax: 300 }] } : {}),
        },
      },
      regions: [{ material: 'water', length: 200e-9, c0: { 'Na+': 10, 'Cl-': 10 }, mixing }],
      contacts: { left: gate(0.3), right: gate(-0.3) },
      grid: { hmin: 0.05e-9, hmax: 5e-9, ratio: 1.1 },
    }).solve();
  for (const crowded of [false, true]) {
    const a = island(0, crowded), b = island(1e-6, crowded);
    assert.ok(b.converged);
    assert.ok(Math.abs(b.plates.left.charge / a.plates.left.charge - 1) < 1e-10);
    for (const mu of Object.values(b.mu)) assert.ok((Math.max(...mu) - Math.min(...mu)) / RT < 1e-9);
    if (crowded) assert.ok(Math.max(...b.c['Na+'], ...b.c['Cl-']) > 200, 'crowded at a gate');
  }
});

test('eddy mixing in a crowded solution mixes concentration, −D_mix Δc, not activity', () => {
  // Eddies carry parcels of fluid, so they mix c. A neutral X on a lattice of cMax = 3, between
  // c = 2 and 1: its molecular flux integrates to D·cMax·ln((cMax − c_R)/(cMax − c_L)), exactly on
  // any grid, and its mixing to D_mix·Δc, at second order in the cells. (Mixing by activity, as
  // before 1.2.8, ran 2.08× too fast here; found by a cold review of src/regions.js.)
  const cMax = 3, Dx = 1e-10, Dm = 1e-9;
  const mu = (c) => RT * (Math.log(c / 1000) + Math.log(cMax / (cMax - c)));
  const want = (Dx * cMax * Math.log((cMax - cR) / (cMax - cL)) + Dm * (cL - cR)) / L;
  const err = (cells) => {
    const sol = new Device({
      species: [{ name: 'X', z: 0, cRef: 1000 }],
      materials: { m: { epsr: 10, species: { X: { D: Dx, mu0: 0 } }, statistics: [{ type: 'lattice', species: ['X'], cMax }] } },
      regions: [{ material: 'm', length: L, mixing: Dm }],
      contacts: {
        left: { species: { X: { type: 'equilibrium', mu: mu(cL) } }, phi: 'neutral' },
        right: { species: { X: { type: 'equilibrium', mu: mu(cR) } }, phi: 'neutral' },
      },
      grid: { hmin: L / cells, hmax: L / cells },
    }).solve();
    assert.ok(sol.converged && sol.iterations <= 8, `${cells} cells`);
    const N = sol.contacts.left.flux.X;
    assert.ok(sol.flux.X.every((x) => Math.abs(x / N - 1) < 1e-12), `${cells} cells: the readout`);
    return N / want - 1;
  };
  const e = [12, 24, 48].map(err);
  assert.ok(Math.abs(e[0]) < 2e-5 && Math.abs(e[2]) < 1e-6, `${e}`);
  for (let k = 1; k < 3; k++) assert.ok(Math.abs(e[k - 1] / e[k] - 4) < 0.05, `second order: ${e}`);
  // Ions in a neutral crowded solution (ε = 0), lattice and Debye–Hückel, where mixing far
  // outruns diffusion: every ion's flux is D_mix·Δc/L, segment by segment (within 1e-4 on 20
  // cells; mixing by activity was 19% off).
  const salt = (c) => ({ bath: { c, reference: 'Cl-' } });
  const left = { 'Na+': 200, 'K+': 20, 'Cl-': 220 }, right = { 'Na+': 50, 'K+': 100, 'Cl-': 150 };
  for (const statistics of [{ type: 'lattice', cMax: 1000 }, { type: 'debye-huckel', epsr: 78.5, a: 4e-10 }]) {
    const ion = (mu0) => ({ D: 1e-16, mu0 });
    const sol = new Device({
      species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
      materials: { water: { epsr: 0, species: { 'Na+': ion(-261.9e3), 'K+': ion(-282.5e3), 'Cl-': ion(-131.2e3) }, statistics: [{ ...statistics, species: ['Na+', 'K+', 'Cl-'] }] } },
      regions: [{ material: 'water', length: L, mixing: 1e-9 }],
      contacts: { left: salt(left), right: { ...salt(right), I: 0 } },
      grid: { hmin: L / 20, hmax: L / 20 },
    }).solve();
    assert.ok(sol.converged, statistics.type);
    for (const sp of ['Na+', 'K+', 'Cl-']) {
      const want = (1e-9 * (left[sp] - right[sp])) / L;
      assert.ok(Math.abs(sol.contacts.left.flux[sp] / want - 1) < 1e-4, `${statistics.type} ${sp}: ${sol.contacts.left.flux[sp]} vs ${want}`);
      sol.flux[sp].forEach((N) => assert.ok(Math.abs(N / want - 1) < 1e-4, `${statistics.type} ${sp}: segment ${N} vs ${want}`));
    }
  }
});

test('flow and mixing are checked', () => {
  const def = bar({}).def;
  def.regions[0].mixing = -1;
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /mixing must be a non-negative number/.test(e.message));
  def.regions[0].mixing = 0;
  def.regions[0].velocity = 'fast';
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /velocity must be a finite number/.test(e.message));
});

test('exchange link at a contact: N = k (μ_out − μ)/RT in series with diffusion', () => {
  // The bar's left end exchanges X with an outside phase at c_out through a rate constant k; the
  // right end holds c_R. Steady: N = D (c_0 − c_R)/L = k ln(c_out/c_0), solved for c_0.
  const k = 1e-4, cOut = 2;
  const dev = bar();
  dev.set({ contacts: { left: { species: { X: { type: 'exchange', k, mu: RT * Math.log(cOut / 1000) } } } } });
  const sol = dev.solve();
  assert.ok(sol.converged);
  let lo = cR, hi = cOut;
  for (let it = 0; it < 200; it++) {
    const c0 = 0.5 * (lo + hi);
    if ((D * (c0 - cR)) / L > k * Math.log(cOut / c0)) hi = c0;
    else lo = c0;
  }
  const want = (D * (lo - cR)) / L;
  assert.ok(Math.abs(sol.contacts.left.flux.X / want - 1) < 1e-10, `${sol.contacts.left.flux.X} vs ${want}`);
  assert.ok(Math.abs(sol.c.X[0] / lo - 1) < 1e-10);
});

test('a neutral species diffusing through two regions of a dielectric: no charge anywhere, so no φ to solve for', () => {
  const RT = GAS_CONSTANT * 298.15;
  const bath = (c) => ({ species: { O2: { type: 'equilibrium', mu: RT * Math.log(c / 1000) } }, phi: 'neutral' });
  for (const epsr of [3, 0]) {
    for (const grid of [{ minCells: 8 }, { hmin: 1e-8, hmax: 1e-6 }]) {
      const sol = new Device({
        species: [{ name: 'O2', z: 0, cRef: 1000 }],
        materials: { polymer: { epsr, species: { O2: { D: 1e-11, mu0: 0 } } } },
        regions: [
          { material: 'polymer', length: 5e-6 },
          { material: 'polymer', length: 5e-6 },
        ],
        contacts: { left: bath(1), right: bath(0.1) },
        grid,
      }).solve();
      assert.ok(sol.converged, `epsr ${epsr}, ${JSON.stringify(grid)}`);
      assert.ok(Math.abs(sol.contacts.left.flux.O2 / ((1e-11 * 0.9) / 10e-6) - 1) < 1e-12);
    }
  }
});
