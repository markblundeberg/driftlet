import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor } from '../src/kit.js';

// A capacitance spread through a port's window: σ = C (V − zeroCharge − φ) per area of electrode
// on the port's side, and the window's charge balance gains aσ (a the electrode's area per
// volume). A gate along a thin channel (the gradual-channel approximation), or the double layer of
// a porous electrode.

const T = 298.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;

// Complex helpers on [re, im] pairs.
const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const div = (a, b) => {
  const q = b[0] * b[0] + b[1] * b[1];
  return [(a[0] * b[0] + a[1] * b[1]) / q, (a[1] * b[0] - a[0] * b[1]) / q];
};
const sqrt = ([x, y]) => {
  const r = Math.hypot(x, y);
  return [Math.sqrt((r + x) / 2), Math.sign(y || 1) * Math.sqrt((r - x) / 2)];
};
const coth = ([x, y]) => {
  const d = Math.cosh(2 * x) - Math.cos(2 * y);
  return [Math.sinh(2 * x) / d, -Math.sin(2 * y) / d];
};

// A thin-film transistor: an undoped film t thick (no body, no holes), its electrons from the
// source and drain contacts, a gate along its length through C per area. Strictly neutral across
// the film's section: F n = a σ, a = 1/t. That's the charge-sheet model, whose current is exact
// from below threshold to saturation:
//   I/W = (μ/L) [(σ_s² − σ_d²)/(2C) + V_T (σ_s − σ_d)],
// with σ at each end set by neutrality against that contact's electron level.
const t = 10e-9, L = 10e-6, Cg = 1e-3, mu = 1e-2, Nc = 46; // m, m, F/m², m²/(V·s), mol/m³ (2.8e19 cm⁻³)
const tft = (Vg, Vd) =>
  new Device({
    T,
    species: [{ name: 'e-', z: -1, cRef: Nc }],
    materials: { film: { epsr: 0, species: { 'e-': { D: mu * VT, mu0: 0 } } } },
    regions: [{ name: 'channel', material: 'film', length: L }],
    contacts: {
      left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' },
      right: { V: Vd, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'neutral' },
    },
    ports: [{ name: 'gate', region: 'channel', V: Vg, area: 1 / t, capacitance: { C: Cg, zeroCharge: 0.3 } }],
    grid: { hmin: L / 4000, hmax: L / 200, ratio: 1.05 }, // (second order: 1e-3 at twice these, deep in saturation)
  });
// σ at a contact holding the electrons at level V (V_e⁻ = V): F n(φ) = σ/t, n = N_c e^((φ − V)/V_T).
const sigmaAt = (Vg, V) => {
  let lo = -5, hi = 5;
  for (let k = 0; k < 200; k++) {
    const phi = (lo + hi) / 2, s = Cg * (Vg - 0.3 - phi);
    if (FARADAY * Nc * Math.exp((phi - V) / VT) * t > s) hi = phi;
    else lo = phi;
  }
  return Cg * (Vg - 0.3 - (lo + hi) / 2);
};

test('a thin-film transistor with its gate along the channel: the charge-sheet current, linear to saturated, above and below threshold', () => {
  for (const Vg of [0.2, 1, 2]) {
    const I = [0.05, 0.5, 2.5, 3].map((Vd) => {
      const sol = tft(Vg, Vd).solve();
      assert.ok(sol.converged);
      const ss = sigmaAt(Vg, 0), sd = sigmaAt(Vg, Vd);
      const want = ((mu / L) * ((ss * ss - sd * sd) / (2 * Cg) + VT * (ss - sd))) / t; // A per m² of the film's section
      assert.ok(Math.abs(-sol.current / want - 1) < 5e-4, `Vg ${Vg}, Vd ${Vd}: ${-sol.current} vs ${want} A/m²`);
      // The gate's charge is the channel's, end to end.
      const sigma = sol.ports[0].sigma;
      const near = (x, y) => Math.abs(x - y) <= 1e-6 * Math.abs(y) + 1e-12; // (deep in saturation σ_d is 0 to the last bit)
      assert.ok(near(sigma[0], ss) && near(sigma.at(-1), sd), `σ ${sigma[0]} vs ${ss}, ${sigma.at(-1)} vs ${sd}`);
      return -sol.current;
    });
    // Saturation: past V_D ≈ V_G − V_T the current hardly grows.
    if (Vg === 2) assert.ok(I[3] / I[2] < 1.01 && I[2] / I[0] > 10, `${I}`);
  }
});

test('a porous electrode\'s double layer on an ionic conductor: the de Levie impedance √(r/y) coth(L√(ry))', () => {
  // Li⁺ on a fixed background (strictly neutral), fed at its left end, with a capacitance spread
  // along it to a port at 0 V: a resistance r = 1/σ per length and a shunt y = iω C_v per volume,
  // C_v the double layer's aC in series with Li⁺'s own chemical capacitance F²c/RT.
  const X = 100, D = 1e-11, Lb = 10e-6, a = 1e6, C = 0.2;
  const dev = new Device({
    T,
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: Lb, fixedCharge: -X * FARADAY }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'dl', region: 'bar', V: 0, area: a, capacitance: { C } }],
    grid: { hmin: Lb / 400, hmax: Lb / 400 },
  });
  const sol = dev.solve();
  const c = sol.c['Li+'][200], sigma = (FARADAY * FARADAY * D * c) / RT;
  const Cv = 1 / (1 / (a * C) + RT / (FARADAY * FARADAY * c));
  const fs = [0.01, 1, 100, 1e4];
  const Z = dev.impedance(fs, { terminal: 'left' });
  fs.forEach((f, q) => {
    const y = [0, 2 * Math.PI * f * Cv], r = [1 / sigma, 0];
    const exact = mul(sqrt(div(r, y)), coth(mul(sqrt(mul(r, y)), [Lb, 0])));
    const err = Math.hypot(Z.Z.re[q] - exact[0], Z.Z.im[q] - exact[1]) / Math.hypot(...exact);
    assert.ok(err < 1e-3, `f = ${f}: ${Z.Z.re[q]} ${Z.Z.im[q]} vs ${exact}, ${err}`);
  });
});

test('the gate\'s impedance at low frequency is its charge capacitance, dQ/dV_G from two steady states', () => {
  // The gate charges the channel through the channel's own resistance: at low frequency that's a
  // pure capacitance, the quasi-static one (the gate in series with the channel's chemical
  // capacitance, integrated along it).
  const Q = (Vg) => tft(Vg, 0.5).solve().ports[0].charge, dV = 1e-4;
  const Cqs = (Q(1 + dV) - Q(1 - dV)) / (2 * dV);
  const fs = [1, 100];
  const Z = tft(1, 0.5).impedance(fs, { terminal: 'gate' });
  fs.forEach((f, q) => {
    const want = -1 / (2 * Math.PI * f * Cqs);
    assert.ok(Math.abs(Z.Z.im[q] / want - 1) < 1e-6, `f = ${f}: ${Z.Z.im[q]} vs ${want} Ω·m²`);
    assert.ok(Math.abs(Z.Z.re[q]) < 1e-3 * Math.abs(want), 'nearly lossless this slowly');
  });
});

test('a double layer charging in a closed, strictly neutral electrolyte keeps every ion it holds', () => {
  // The ions supply the double layer's countercharge, so their net charge isn't zero but −aσ;
  // the solver's charge-continuity row has to count its change, or the most abundant ion leaks.
  const ions = [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }];
  const dev = new Device({
    T,
    species: ions,
    materials: { water: { epsr: 0, species: { 'Na+': { D: 1.3e-9, mu0: -261.9e3 }, 'Cl-': { D: 2e-9, mu0: -131.2e3 } } } },
    regions: [{ material: 'water', length: 1e-3, c0: { 'Na+': 500, 'Cl-': 500 } }],
    contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ name: 'dl', region: 0, from: 0.5e-3, V: { t: [0, 1e-3], values: [0, 0.2] }, area: 1e4, capacitance: { C: 0.2 } }],
    grid: { hmin: 10e-6, hmax: 50e-6 },
  });
  for (const time of [1e-3, 1, 100]) {
    const s = dev.advance(time);
    assert.ok(s.converged);
    for (const st of s.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `${st.species} at ${time} s: ${st.drift}`);
  }
});

test('two capacitances on one window add: two gates of C/2 hold what one of C does, ions conserved', () => {
  const ions = [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }];
  const gate = (name, C) => ({ name, region: 0, from: 0.5e-3, V: { t: [0, 1e-3], values: [0, 0.2] }, area: 1e4, capacitance: { C } });
  const run = (ports) => {
    const dev = new Device({
      T,
      species: ions,
      materials: { water: { epsr: 0, species: { 'Na+': { D: 1.3e-9, mu0: -261.9e3 }, 'Cl-': { D: 2e-9, mu0: -131.2e3 } } } },
      regions: [{ material: 'water', length: 1e-3, c0: { 'Na+': 500, 'Cl-': 500 } }],
      contacts: { left: { phi: 'neutral' }, right: { phi: 'neutral' } },
      ports,
      grid: { hmin: 10e-6, hmax: 50e-6 },
    });
    return [1e-3, 1].map((time) => {
      const s = dev.advance(time);
      assert.ok(s.converged);
      for (const st of s.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `${st.species} at ${time} s: ${st.drift}`);
      return s.ports.reduce((t, p) => t + p.charge, 0);
    });
  };
  const one = run([gate('a', 0.2)]), two = run([gate('a', 0.1), gate('b', 0.1)]);
  one.forEach((q, j) => assert.ok(Math.abs(two[j] / q - 1) < 1e-6, `${two[j]} vs ${q}`));
});

test('set() patches ports by name: a gate\'s voltage is a drive, changed in place', () => {
  const dev = tft(1, 0.5), solver = dev.solver, before = dev.solve().current;
  dev.set({ ports: { gate: { V: 2 } } });
  assert.equal(dev.solver, solver, 'a drive alone: the same solver, its state kept');
  const after = dev.solve().current;
  assert.ok(Math.abs(after / tft(2, 0.5).solve().current - 1) < 1e-9 && Math.abs(after) > 2 * Math.abs(before));
  assert.equal(dev.def.ports[0].region, 'channel', 'the rest of the port kept');
  dev.set({ ports: { gate: { I: 0 } } }); // a new kind of drive drops the old
  assert.equal(dev.def.ports[0].V, undefined);
  assert.throws(() => dev.set({ ports: { gat: { V: 1 } } }), (e) => e instanceof DeviceError && /no port named 'gat' \(the ports: gate\)/.test(e.message));
  assert.throws(() => dev.set({ ports: 3 }), (e) => e instanceof DeviceError && /ports must be an array/.test(e.message));
});

test('a capacitance is checked', () => {
  const def = tft(1, 0.1).def;
  def.ports[0].capacitance = { C: -1 };
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /capacitance\.C/.test(e.message));
  def.ports[0].capacitance = { C: 1e-3 };
  delete def.ports[0].area;
  assert.throws(() => new Device(def), /area: give the electrode's area/);
});

test('a floating electrode with a double layer starts uncharged and charges to its rest potential', () => {
  // Silver spread through AgNO₃ at I: 0, against a silver contact: the double layer starts empty
  // (as the neutral starting composition has it), the plating current charges it, and the
  // electrode ends level with the contact's silver, its charge C(V − zeroCharge − φ) per area.
  const ions = [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }];
  const def = (port) => ({
    T,
    species: ions,
    materials: { water: { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
    regions: [{ material: 'water', length: 10e-6, c0: { 'NO3-': 10, 'Ag+': 10 } }],
    contacts: { left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium' }, phi: 'neutral' }, right: {} },
    ports: [{ name: 'el', region: 0, from: 5e-6, terminal: 'e-', area: 1e5, capacitance: { C: 0.2, zeroCharge: 0.1 },
      reactions: [{ equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-4, alpha: 0.5 }], ...port }],
    grid: { hmin: 0.5e-6, hmax: 0.5e-6 },
  });
  const dev = new Device(def({ I: 0 }));
  assert.ok(Math.abs(dev.solution().ports[0].charge) < 1e-15, 'starts uncharged');
  const s = dev.advance(100);
  assert.ok(s.converged);
  const p = s.ports[0];
  assert.ok(Math.abs(p.V - s.V['Ag+'][0]) < 1e-6, `rests level with the silver: ${p.V} vs ${s.V['Ag+'][0]}`);
  assert.ok(Math.abs(p.charge) > 1e-3, `charged: ${p.charge}`);
  // With every held terminal passing nothing, a floating one's level is nobody's: an error.
  const closed = def({ I: 0 });
  closed.contacts.left = {};
  assert.throws(() => new Device(closed), (e) => e instanceof DeviceError && /no terminal held at a voltage passes any/.test(e.message));
});

test('a capacitance driven by a current only charges: no steady state, but a transient at I·t, and at no current it keeps its charge', () => {
  // A double layer in NaCl, fed 1 mA/m² from an uncharged start, the current returning through
  // an electrode that exchanges Cl⁻ (which also sets the potentials' level).
  const dev = () =>
    new Device({
      T,
      species: [{ name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
      materials: { water: { epsr: 0, species: { 'Na+': { D: 1.3e-9, mu0: -261.9e3 }, 'Cl-': { D: 2e-9, mu0: -131.2e3 } } } },
      regions: [{ material: 'water', length: 1e-3, c0: { 'Na+': 500, 'Cl-': 500 } }],
      contacts: { left: { V: 0, terminal: 'Cl-', species: { 'Cl-': 'equilibrium' }, phi: 'neutral' }, right: { phi: 'neutral' } },
      ports: [{ name: 'dl', region: 0, from: 0.5e-3, I: 1e-3, area: 1e4, capacitance: { C: 0.2 } }],
      grid: { hmin: 10e-6, hmax: 50e-6 },
    });
  const none = dev().solve();
  assert.ok(!none.converged && none.warnings.some((w) => w.startsWith('dl: driven at 0.00100 A/m², it has no steady state: it passes current only by charging its capacitance')), none.warnings.join('\n'));
  const d = dev();
  assert.ok(Math.abs(d.solution().ports[0].charge) < 1e-12, 'starts uncharged');
  for (const t of [1e-3, 1]) {
    const s = d.advance(t);
    assert.ok(s.converged && Math.abs(s.ports[0].charge / (1e-3 * t) - 1) < 1e-6, `at ${t} s: ${s.ports[0].charge} C/m²`);
  }
  // At no current it keeps the charge it holds: its steady state.
  const kept = d.set({ ports: { dl: { I: 0 } } }).solve();
  assert.ok(kept.converged && Math.abs(kept.ports[0].charge / 1e-3 - 1) < 1e-9, `${kept.ports[0].charge} C/m²`);
});

test('a gate left floating (I = 0) keeps its charge: moving the back contact moves it alike', () => {
  // p-Si under a gate's capacitance, held at 0.5 V and solved, then opened: its charge stays, so
  // where the back contact goes, the whole device and the gate go too.
  const Si = semiconductor('Si');
  const dev = new Device(build({
    T: 300,
    library: [Si],
    stack: [ohmic(0), layer('Si', 1e-6, { acceptors: units.perCm3(1e17) }), { V: 0.5, phi: { type: 'capacitive', C: 3.45e-3 }, zeroCharge: 0 }],
    grid: { hmin: 1e-9, hmax: 2e-8 },
  }));
  const held = dev.solve();
  const Q = held.terminals.right.charge ?? held.contacts.right.D;
  dev.set({ contacts: { right: { I: 0 } } });
  const open = dev.solve();
  assert.ok(open.converged && Math.abs(open.terminals.right.V - 0.5) < 1e-9, `${open.terminals.right.V} V`);
  const moved = dev.set({ contacts: { left: { V: 0.1 } } }).solve();
  assert.ok(moved.converged && Math.abs(moved.terminals.right.V - 0.6) < 1e-9, `${moved.terminals.right.V} V`);
  assert.ok(Math.abs((moved.terminals.right.charge ?? moved.contacts.right.D) / Q - 1) < 1e-9);
});
