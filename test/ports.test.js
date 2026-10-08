import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, metal, aqueous, bath, check } from '../src/kit.js';

const RT = GAS_CONSTANT * 298.15;

// A single-ion conductor (Li⁺ on a fixed background, strictly neutral, so exactly ohmic) fed at
// its left end and drained all along its length by a port at 0 V: a transmission line.
const X = 100, D = 1e-11, L = 10e-6;
const sigma = (FARADAY * FARADAY * D * X) / RT;
const line = (G, cells, V1 = 0.01) =>
  new Device({
    species: [{ name: 'Li+', z: 1, cRef: 1000 }],
    materials: { solid: { epsr: 0, species: { 'Li+': { D, mu0: 0 } } } },
    regions: [{ name: 'bar', material: 'solid', length: L, fixedCharge: -X * FARADAY }],
    contacts: { left: { V: V1, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { phi: 'neutral' } },
    ports: [{ region: 'bar', V: 0, terminal: 'Li+', species: { 'Li+': { type: 'conductance', G } } }],
    grid: { minCells: cells },
  });

test('distributed conductance to a port: transmission-line current σV₁ tanh(L/λ)/λ, second order', () => {
  for (const G of [1e6, 1e9]) {
    const lam = Math.sqrt(sigma / G);
    const want = ((sigma * 0.01) / lam) * Math.tanh(L / lam);
    const err = [100, 400].map((cells) => {
      const sol = line(G, cells).solve();
      assert.ok(sol.converged);
      // Current in = current out through the port (the right end is blocked).
      assert.ok(Math.abs(sol.contacts.left.current + sol.ports[0].current - sol.contacts.right.current) < 1e-12 * want);
      for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12);
      return Math.abs(sol.contacts.left.current / want - 1);
    });
    assert.ok(err[1] < 1e-4 && err[0] / err[1] > 12, `G=${G}: errors ${err}`);
  }
});

test('transient through a port: the bookkeeping includes what the port brings in', () => {
  const dev = line(1e9, 100, 0);
  dev.solve();
  dev.set({ contacts: { left: { V: 0.02 } } });
  const sol = dev.advance(1e-3);
  for (const st of sol.conservation) assert.ok(Math.abs(st.drift) < 1e-12, `drift ${st.drift}`);
});

test('a held port level: a neutral species pinned in a window, fed by diffusion through the rest', () => {
  // X diffuses from a reservoir at μ_L on the left; a port holds μ_P in the right half. In
  // steady state the left half carries the linear profile between the two levels.
  const cL = 2, cP = 0.5, Dx = 1e-9, Lx = 2e-6;
  const sol = new Device({
    species: [{ name: 'X', z: 0, cRef: 1000 }],
    materials: { m: { epsr: 10, species: { X: { D: Dx, mu0: 0 } } } },
    regions: [{ name: 'm', material: 'm', length: Lx }],
    contacts: { left: { species: { X: { type: 'equilibrium', mu: RT * Math.log(cL / 1000) } }, phi: 'neutral' }, right: { phi: 'neutral' } },
    ports: [{ region: 'm', from: Lx / 2, V: 0, species: { X: { type: 'equilibrium', mu: RT * Math.log(cP / 1000) } } }],
    grid: { hmin: Lx / 100, hmax: Lx / 100 },
  }).solve();
  assert.ok(sol.converged);
  const N = (Dx * (cL - cP)) / (Lx / 2);
  assert.ok(Math.abs(sol.contacts.left.flux.X / N - 1) < 1e-9);
  assert.ok(Math.abs(sol.ports[0].flux.X / -N - 1) < 1e-9, 'the port absorbs what enters');
});

// An electrode spread through a window: a port with reactions, its carrier at the port's level.
// In linear response each node's reaction is a conductance to the port, so the single-ion bar
// with lithium plating along it is the transmission line above, G = F² a k₀ (c/c_ref)^(1−α)/RT.
test('an electrode spread along a bar (reactions on a port): in linear response, the transmission line', () => {
  for (const G of [1e6, 1e9]) {
    const lam = Math.sqrt(sigma / G), V1 = 1e-5;
    const want = ((sigma * V1) / lam) * Math.tanh(L / lam);
    const area = 1e6, k0 = (G * RT) / (FARADAY * FARADAY * area * Math.sqrt(X / 1000));
    const err = [100, 400].map((cells) => {
      const def = line(G, cells, V1).def;
      def.species.push({ name: 'e-', z: -1 });
      def.ports = [{ region: 'bar', V: 0, terminal: 'e-', area, reactions: [{ equation: 'Li+ + e- = Li(s)', fixed: { 'Li(s)': 0 }, k0, alpha: 0.5 }] }];
      const sol = new Device(def).solve();
      assert.ok(sol.converged);
      assert.ok(Math.abs(sol.contacts.left.current + sol.ports[0].current) < 1e-9 * want);
      // Each node's rate, per area of electrode, follows the line's cosh profile.
      const r = sol.ports[0].rates[0], x = sol.ports[0].x;
      const at = (xx) => (G * V1 * Math.cosh((L - xx) / lam)) / Math.cosh(L / lam) / (FARADAY * area);
      assert.ok(Math.abs(r[0] / at(x[0]) - 1) < 2e-3 && Math.abs(r.at(-1) / at(x.at(-1)) - 1) < 2e-3, `${r[0]} vs ${at(x[0])}`);
      return Math.abs(sol.contacts.left.current / want - 1);
    });
    // (A line much longer than λ's scale, G = 1e6, is nearly uniform: already at the floor.)
    assert.ok(err[1] < 1e-4 && (err[0] / err[1] > 12 || err[0] < 1e-6), `G=${G}: errors ${err}`);
  }
});

test('a floating electrode spread through a solution corrodes at the Wagner–Traud mixed potential, between the two couples\' levels', () => {
  // Iron in acid: Fe²⁺ + 2e⁻ ⇌ Fe and 2H⁺ + 2e⁻ ⇌ H₂, each Butler–Volmer with α = ½, on iron
  // spread thinly through 10 µm of solution beside a bath (so the composition is the bath's).
  // At I = 0 the iron's electron level V sits where the two rates cancel. With α = ½ that's
  // exact: e^{2V/V_T} = (k₁e^{v₁} + k₂e^{v₂})/(k₁e^{−v₁} + k₂e^{−v₂}), where v = V_i/V_T is each
  // couple's level (Fe²⁺'s V_i for iron, H⁺'s for hydrogen) and k the prefactor.
  const T = 298.15, VT = (GAS_CONSTANT * T) / FARADAY;
  const lib = aqueous(['Fe2+', 'H+', 'Cl-'], { epsr: 0 });
  lib.species.push({ name: 'e-', z: -1 });
  const c = { 'Fe2+': 10, 'H+': 10, 'Cl-': 30 }, k1 = 1e-7, k2 = 1e-9;
  let V0;
  for (const I of [0, 0.05]) {
    const dev = new Device(build({
      T,
      library: [lib],
      stack: [bath(c, 'Cl-'), layer('water', 10e-6, { name: 'film' }), {}],
      ports: [{ name: 'iron', region: 'film', I, terminal: 'e-', area: 1e3, reactions: [
        { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: k1, alpha: 0.5 },
        { equation: '2 H+ + 2 e- = H2', fixed: { H2: 0 }, k0: k2, alpha: 0.5 },
      ] }],
      grid: { hmin: 0.5e-6, hmax: 0.5e-6 },
    }));
    const sol = dev.solve();
    assert.ok(sol.converged);
    const p = sol.ports[0], g = 0;
    // Electrons the forward reactions take from the iron, per second, over the window: −I/F.
    const len = (w) => ((w < p.x.length - 1 ? p.x[w + 1] - p.x[w] : 0) + (w > 0 ? p.x[w] - p.x[w - 1] : 0)) / 2;
    const taken = p.x.reduce((sum, _, w) => sum + 2 * (p.rates[0][w] + p.rates[1][w]) * p.area[w] * len(w), 0);
    assert.ok(Math.abs(-FARADAY * taken - I) < 1e-9 * Math.max(I, FARADAY * Math.abs(p.rates[1][0]) * 1e3 * 10e-6), `${-FARADAY * taken} vs ${I} A/m²`);
    if (I === 0) {
      const v1 = sol.V['Fe2+'][g] / VT, v2 = sol.V['H+'][g] / VT;
      const K1 = k1 * Math.sqrt(sol.c['Fe2+'][g] / 1000), K2 = k2 * (sol.c['H+'][g] / 1000);
      const V = (VT / 2) * Math.log((K1 * Math.exp(v1) + K2 * Math.exp(v2)) / (K1 * Math.exp(-v1) + K2 * Math.exp(-v2)));
      assert.ok(Math.abs(p.V - V) < 1e-6, `${p.V} vs ${V} V`);
      assert.ok(v1 * VT < p.V && p.V < v2 * VT);
      assert.ok(p.rates[0][0] < 0 && p.rates[1][0] > 0, 'iron dissolves (backward) as hydrogen comes off (forward)');
    } else assert.ok(p.V > V0 + 0.01, 'driven anodically, the iron sits higher');
    V0 = p.V;
  }
});

test('an exchange port feeding a face reaction in a closed, strictly neutral electrolyte: the transient runs from the start, the iron corroding at the mixed potential', () => {
  // A drop of salt water on iron (two iron faces, one piece of metal), O₂ supplied from the air
  // by a port over part of it and reduced at the faces while the iron dissolves: the water's
  // potential is set only by the face reactions. The port's window was once left out of the
  // better-conditioned terms a strictly neutral transient is solved in, and no step converged.
  const water = aqueous(['Na+', 'Cl-', 'Fe2+', 'OH-'], { epsr: 0 });
  water.species.push({ name: 'O2', z: 0, cRef: 1000 });
  water.materials.water.species.O2 = { D: 2e-9, mu0: 16.4e3 };
  const iron = { species: [{ name: 'e-', z: -1 }], materials: { Fe: { conductor: { species: 'e-', conductivity: 1e7 } } } };
  const face = {
    reactions: [
      { equation: 'Fe2+ + 2 e- = Fe(s)', fixed: { 'Fe(s)': 0 }, k0: 5e-6, alpha: 0.5 },
      { equation: 'O2 + 2 H2O + 4 e- = 4 OH-', fixed: { H2O: -237.129e3 }, k0: 3e-9, alpha: 0.125 },
    ],
  };
  const Lw = 1e-3, cO2 = 0.26;
  const def = (drive = {}) =>
    build({
      library: [water, iron],
      stack: [ohmic(0, ['e-']), layer('Fe', 1e-4), face, layer('water', Lw, { name: 'drop', c0: { 'Na+': 500, 'Cl-': 500, 'Fe2+': 1e-6, 'OH-': 2e-6, O2: cO2 } }), face, layer('Fe', 1e-4), ohmic(0, ['e-'])],
      ports: [{ name: 'air', region: 'drop', from: Lw / 2, species: { O2: { type: 'exchange', k: 0.05, mu: 16.4e3 + RT * Math.log(cO2 / 1000) } }, ...drive }],
      grid: { hmin: 5e-6, hmax: 50e-6 },
    });
  const dev = new Device(def());
  for (const t of [1, 1000]) {
    const s = dev.advance(t);
    assert.ok(s.converged, `to ${t} s: ${s.steps} steps, ${s.rejected} rejected`);
    assert.ok(check(dev, s, { refine: false }).ok, check(dev, s, { refine: false }).text);
    // No current leaves the drop: what the iron loses, the O₂ takes, at whatever potential that needs.
    const [a, b] = s.interfaces.map((f) => f.rates), iFe = -2 * FARADAY * (a[0] + b[0]), iO2 = 4 * FARADAY * (a[1] + b[1]);
    assert.ok(iFe > 0.1 && Math.abs(iO2 / iFe - 1) < 1e-9, `${iFe} vs ${iO2} A/m²`);
  }
  // A port that exchanges only neutral species carries no current, so can't be driven by one.
  assert.throws(() => new Device(def({ I: 0 })), (e) => e instanceof DeviceError && /passes no current/.test(e.message));
});

test('MOS: a port grounding the channel gives the low-frequency C–V, inversion included', () => {
  // p-Si under a gate (right). Without generation, inversion electrons could only arrive by
  // minority-carrier diffusion from the back contact (minutes, through a bulk with ~1e3 of them
  // per cm³), and the inversion layer's Fermi level is held only weakly (see the next tests). A
  // port holding the
  // electrons at ground beside the oxide (as source and drain would, in 2D) anchors the channel,
  // and the small-signal capacitance then follows the quasi-static C–V.
  const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19), NA = units.perCm3(1e17);
  const Cg = (3.9 * EPS0) / 5e-9, Lsi = 0.5e-6;
  const mos = (port) =>
    new Device({
      species: [
        { name: 'e-', z: -1 },
        { name: 'h+', z: 1 },
      ],
      materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } } },
      regions: [{ name: 'Si', material: 'Si', length: Lsi, fixedCharge: -NA * FARADAY }],
      contacts: {
        left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' },
        right: { V: -3, phi: { type: 'capacitive', C: Cg }, zeroCharge: -0.9 },
      },
      ports: port ? [{ region: 'Si', from: Lsi - 5e-9, V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' } }] : undefined,
      grid: { hmin: 0.1e-9, hmax: 10e-9, ratio: 1.1 },
    });
  const w = 2 * Math.PI * 1e3;
  const C1k = (d) => -1 / (w * d.impedance([1e3]).Z.im[0]);
  const grounded = mos(true), floating = mos(false);
  const seen = [];
  for (const V of [-3, -1.5, -1, 0, 1]) {
    grounded.set({ contacts: { right: { V } } });
    const C = C1k(grounded);
    const q = (v) => (grounded.set({ contacts: { right: { V: v } } }), grounded.solve().gates.right.charge);
    const Cqs = Math.abs(q(V + 1e-4) - q(V - 1e-4)) / 2e-4;
    grounded.set({ contacts: { right: { V } } });
    assert.ok(Math.abs(C / Cqs - 1) < 2e-3, `V=${V}: ${C / Cg} vs quasi-static ${Cqs / Cg} (× C_ox)`);
    seen.push(C / Cg);
    if (V < -1.2) {
      // Accumulation and depletion: no channel yet, so the port changes nothing.
      floating.set({ contacts: { right: { V } } });
      assert.ok(Math.abs(C1k(floating) / C - 1) < 1e-6, `V=${V}`);
    }
  }
  // The full low-frequency curve: accumulation near C_ox, a depletion dip, inversion back up.
  assert.ok(seen[0] > 0.9 && Math.min(...seen) < 0.2 && seen.at(-1) > 0.9, `C/C_ox: ${seen.map((c) => c.toFixed(3))}`);
});

test('MOS without a port: its steady states are equilibrium, the same as with the channel held, on any grid', () => {
  // Without generation, the inversion electrons reach only the back contact (the oxide blocks
  // them), so in a steady state they carry nothing and their level is the contact's, flat:
  // equilibrium, which a port holding the channel at that level gives too. The steady solve pins
  // such a level outright; found through the bulk's ~1e3 cm⁻³ of electrons instead, it failed at
  // some gate voltages and grids and not others.
  const Lsi = units.um(0.5), tox = units.nm(5), Cox = (3.9 * EPS0) / tox;
  const def = (hmin, port) => build({
    T: 300,
    library: [semiconductor('Si'), metal('Au')],
    materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [
      ohmic(0),
      layer('Si', Lsi, { name: 'p-Si', acceptors: units.perCm3(1e17) }),
      { dipole: 0 },
      layer('SiO2', tox, { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }),
      { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.1 },
      layer('Au', units.nm(50)),
      ohmic(0, ['e-']),
    ],
    ports: port ? [{ region: 'p-Si', from: Lsi - units.nm(5), V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' } }] : undefined,
    grid: { hmin, hmax: units.nm(10), ratio: 1.1 },
  });
  const charge = (d, V) => {
    d.set({ contacts: { right: { V } } });
    const s = d.solve();
    assert.ok(s.converged, `V = ${V}`);
    return s.interfaces.at(-1).D; // the gate's charge, C/m²
  };
  for (const hmin of [units.nm(0.1), units.nm(0.13), units.nm(0.2)]) {
    const bare = new Device(def(hmin)), held = new Device(def(hmin, true)), C = [];
    for (const V of [-1.5, -0.5, 0, 0.3, 0.6, 1, 1.5]) {
      const q = charge(bare, V), want = charge(held, V);
      assert.ok(Math.abs(q - want) < 1e-6 * Math.abs(Cox * 1.5), `hmin ${hmin}, V = ${V}: ${q} vs ${want} C/m²`);
      C.push(q);
    }
    // Accumulation, depletion, then inversion: its electrons add charge fast above threshold.
    assert.ok(C[0] > 0 && C.at(-1) < C[3] && (C.at(-1) - C.at(-2)) / 0.5 < 0.5 * -Cox, `${C}`);
    // Cold too, at each voltage
    for (const V of [0.3, 0.4, 0.9]) assert.ok(new Device(def(hmin)).set({ contacts: { right: { V } } }).solve().converged, `cold at ${V} V, hmin ${hmin}`);
  }
});

test('MOS without a port: after a gate step, the inversion layer fills over minutes to the low-frequency charge', () => {
  // The electrons come by minority diffusion from the back contact, so the gate's charge creeps
  // from the high-frequency value to the low-frequency one, dQ/dV from two steady states. The
  // inversion layer's 0.1 nm cells couple their electrons ~1e18 times more tightly than its
  // storage over a second: Newton's factorised solves lose the layer's level to round-off, and
  // refine it by GMRES (it once took 80,000 steps to get here, 1.5% off).
  const Lsi = units.um(0.5), tox = units.nm(5), Cox = (3.9 * EPS0) / tox;
  const def = build({
    T: 300,
    library: [semiconductor('Si'), metal('Au')],
    materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [
      ohmic(0),
      layer('Si', Lsi, { name: 'p-Si', acceptors: units.perCm3(1e17) }),
      { dipole: 0 },
      layer('SiO2', tox, { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }),
      { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.1 },
      layer('Au', units.nm(50)),
      ohmic(0, ['e-']),
    ],
    grid: { hmin: units.nm(0.1), hmax: units.nm(10), ratio: 1.1 },
  });
  const at = (V) => new Device(def).set({ contacts: { right: { V } } });
  const D = (V) => at(V).solve().interfaces.at(-1).D;
  const low = (D(0.61) - D(0.6)) / 0.01;
  const dev = at(0.6), D0 = dev.solve().interfaces.at(-1).D;
  dev.set({ contacts: { right: { V: 0.61 } } });
  const s = dev.advance(1e5);
  assert.ok(s.converged && s.done && s.steps < 1000, `${s.steps} steps`);
  const C = (s.interfaces.at(-1).D - D0) / 0.01;
  assert.ok(Math.abs(C / low - 1) < 1e-3 && low / -Cox > 0.85, `${C / -Cox} vs ${low / -Cox} × C_ox`);
  assert.ok(Math.abs(s.current) < 1e-12, `still ${s.current} A/m² through the oxide`);
});

test('MOS without a port: the impedance gives the high-frequency C–V at 1 Hz, as a time-domain run does', () => {
  // The inversion layer's electrons can only come by minority diffusion from the back contact,
  // through a bulk with ~1e3 cm⁻³ of them: minutes. The assembled Jacobian alone (two huge entries
  // per inversion-layer flux, nearly cancelling) gave the layer an exchange path that followed
  // the gate at 1 Hz; the impedance's GMRES now takes J·v with each flux on its η difference.
  const Lsi = units.um(0.5), tox = units.nm(5), Cox = (3.9 * EPS0) / tox;
  const def = build({
    T: 300,
    library: [semiconductor('Si'), metal('Au')],
    materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [
      ohmic(0),
      layer('Si', Lsi, { name: 'p-Si', acceptors: units.perCm3(1e17) }),
      { dipole: 0 },
      layer('SiO2', tox, { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }),
      { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.1 },
      layer('Au', units.nm(50)),
      ohmic(0, ['e-']),
    ],
    grid: { hmin: units.nm(0.1), hmax: units.nm(10), ratio: 1.1 },
  });
  const at = (V) => {
    const d = new Device(def);
    d.set({ contacts: { right: { V } } });
    assert.ok(d.solve().converged);
    return d;
  };
  const dev = at(0.6);
  const C = (f) => -1 / (2 * Math.PI * f * dev.impedance([f]).Z.im[0]) / Cox;
  const [low, high] = [C(1e-6), C(1)];
  assert.ok(low > 0.85 && high < 0.15, `C/C_ox ${low} at 1 µHz, ${high} at 1 Hz`);
  // A 10 mV gate step, stepped in time: the charge at 1 s is the high-frequency value.
  const run = at(0.6);
  run.set({ contacts: { right: { V: 0.61 } } });
  let q = 0, t0 = 0;
  const s = run.advance(1);
  s.trace.t.forEach((t, k) => ((q += s.trace.current[k] * (t - t0)), (t0 = t)));
  assert.ok(Math.abs(Math.abs(q) / (Cox * 0.01) / high - 1) < 0.1, `step: ${Math.abs(q) / (Cox * 0.01)} vs ${high}`);
});

test('a level held by one port where another\'s electrode reacts: the same in either order, and the current is conserved', () => {
  // Ag⁺ held at V = 0 through a window where silver plates at 10 mV: the hold supplies what the
  // electrode takes, so their currents cancel. A held level replaces its node's balance, and a
  // reaction once added into the replaced row when its port came after the hold's.
  const ions = [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }];
  const water = { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } };
  const window = { region: 0, from: 0.4e-6, to: 0.6e-6 };
  const hold = { name: 'hold', ...window, V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium' } };
  const plate = { name: 'el', ...window, V: 0.01, terminal: 'e-', area: 1e6, reactions: [{ equation: 'Ag+ + e- = Ag(s)', fixed: { 'Ag(s)': 0 }, k0: 1e-3, alpha: 0.4 }] };
  const currents = [[hold, plate], [plate, hold]].map((ports) => {
    const sol = new Device({
      species: ions,
      materials: { water },
      regions: [{ material: 'water', length: 1e-6, c0: { 'NO3-': 10, 'Ag+': 10 } }],
      contacts: { left: { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' }, right: {} },
      ports,
      grid: { hmin: 2e-9, hmax: 100e-9, ratio: 1.5, minCells: 4 },
    }).solve();
    assert.ok(sol.converged);
    const total = Object.values(sol.terminals).reduce((t, x) => t + x.current, 0);
    assert.ok(Math.abs(total) < 1e-9, `the terminals' currents sum to ${total}`);
    return Object.fromEntries(sol.ports.map((p) => [p.name, p.current]));
  });
  assert.ok(currents[0].el > 0.1 && Math.abs(currents[0].el + currents[0].hold) < 1e-9 * currents[0].el, JSON.stringify(currents[0]));
  assert.ok(Math.abs(currents[1].el / currents[0].el - 1) < 1e-12, JSON.stringify(currents));
});

test('ports are checked', () => {
  const def = line(1e6, 10).def;
  def.ports[0].species['Li+'] = { type: 'equilibrium', mu: 0 };
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /not by mu/.test(e.message));
  def.ports[0].species['Li+'] = { type: 'conductance', G: 1, offset: 0, mu: 0 };
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /Li\+\.mu: not a field here/.test(e.message));
  def.ports[0].species['Li+'] = { type: 'exchange', k: 1 };
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /conductance/.test(e.message));
  def.ports[0].species['Li+'] = { type: 'conductance', G: 1 };
  def.ports[0].region = 'nowhere';
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /region's name or index/.test(e.message));
  def.ports[0].region = 0;
  def.ports[0].from = 2.5e-6;
  def.ports[0].to = 2.6e-6; // between nodes (1 µm cells)
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /holds no grid node/.test(e.message));
  // An electrode port: its reactions need its carrier and its area per volume.
  const el = line(1e6, 10).def;
  el.species.push({ name: 'e-', z: -1 });
  const plate = { equation: 'Li+ + e- = Li(s)', fixed: { 'Li(s)': 0 }, k0: 1e-6, alpha: 0.5 };
  el.ports = [{ region: 'bar', reactions: [plate], area: 1e6 }];
  assert.throws(() => new Device(el), (e) => e instanceof DeviceError && /a reacting port's terminal/.test(e.message));
  el.ports = [{ region: 'bar', terminal: 'e-', reactions: [plate] }];
  assert.throws(() => new Device(el), (e) => e instanceof DeviceError && /area: give the electrode's area per volume/.test(e.message));
  el.ports = [{ region: 'bar', terminal: 'e-', reactions: [{ ...plate, equation: 'Li+ = Li(s)' }], area: 1e6 }];
  assert.throws(() => new Device(el), (e) => e instanceof DeviceError && /no 'e-' from the electrode takes part/.test(e.message));
  el.ports = [{ region: 'bar', terminal: 'e-', reactions: [plate], area: 1e6 }];
  el.materials.solid.epsr = 10;
  assert.throws(() => new Device(el), (e) => e instanceof DeviceError && /needs its region strictly neutral/.test(e.message));
});

test('a port driven by a current pulse: it delivers I·t_p, and switching off is quick however far its voltage falls', () => {
  // Holes injected for 0.5 µs into n-Si under a field, with no recombination: the port's
  // voltage floats hundreds of volts above the window at this current (small G), then collapses
  // when the current stops. Holes are conserved until they reach the right contact.
  const ND = units.perCm3(1e15), Lb = 3e-3, E = 3000, tp = 0.5e-6, I = (ND * FARADAY * 50e-6) / tp; // A/m²: about n₀ over the window
  const dev = new Device(
    build({
      T: 300,
      library: [semiconductor('Si')],
      stack: [ohmic(E * Lb), layer('Si', Lb, { name: 'bar', donors: ND }), ohmic(0)],
      ports: [{ name: 'emitter', region: 'bar', from: 0.5e-3, to: 0.55e-3, terminal: 'h+', species: { 'h+': { type: 'conductance', G: 1e6 } },
        I: { t: [0, 1e-9, tp, tp + 1e-9], values: [0, I, I, 0] } }],
      grid: { hmin: 10e-6, hmax: 10e-6 },
    }),
  );
  const rest = dev.solve();
  const holes = (s) => s.x.reduce((a, x, g) => (g ? a + ((s.c['h+'][g] - rest.c['h+'][g] + s.c['h+'][g - 1] - rest.c['h+'][g - 1]) / 2) * (x - s.x[g - 1]) : 0), 0);
  const on = dev.advance(tp);
  assert.ok(on.ports[0].V > 100, `the port floats far above the window while driven (${on.ports[0].V} V)`);
  const off = dev.advance(tp + 50e-9);
  assert.ok(off.converged && off.steps + off.rejected < 40, `switching off took ${off.steps} steps and ${off.rejected} rejections`);
  assert.ok(Math.abs(off.ports[0].V) < 20, `the port's voltage collapses once the current stops (${off.ports[0].V} V)`);
  // Everything injected is still in the bar (the pulse is far from the right contact).
  assert.ok(Math.abs(holes(off) / ((I * tp) / FARADAY) - 1) < 1e-3, `${holes(off)} mol/m² of holes`);
});

test('a capacitance-only port given a terminal species that no material holds is told it needs none', () => {
  const def = {
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }],
    materials: { water: { epsr: 0, species: { 'K+': { D: 2e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 } } } },
    regions: [{ name: 'pore', material: 'water', length: 1e-4, c0: { 'K+': 100, 'Cl-': 100 } }],
    contacts: { left: { phi: 'neutral' }, right: { V: 0, bath: { c: { 'K+': 100, 'Cl-': 100 } } } },
    ports: [{ name: 'mat', region: 'pore', V: 0, terminal: 'e-', area: 1e6, capacitance: { C: 0.2, zeroCharge: 0 } }],
  };
  assert.throws(() => new Device(def), /terminal of port "mat".*only with reactions.*needs no terminal/);
});

test('a port holding a species at a device end node that the contact there blocks: a window of that node alone passes current', () => {
  // A particle charged through its surface: Li⁺ from the electrolyte (the right contact), e⁻ from
  // a wire (a port) on the surface node alone. Every charge driven in arrives as inserted Li.
  const cMax = 30000, a = 5e-6, I = -1e-12; // A, a total in a sphere
  const dev = new Device({
    species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: {
      host: {
        epsr: 0,
        species: { 'Li+': { D: 1e-14, mu0: 0, cRef: cMax }, 'e-': { D: 1e-4, mu0: 0, cRef: cMax } },
        statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax, ocv: { x: [0.01, 0.5, 0.99], E: [4, 3.9, 3.8], muRef: 0 } }],
      },
    },
    regions: [{ name: 'p', material: 'host', length: a, c0: { 'Li+': 0.4 * cMax, 'e-': 0.4 * cMax } }],
    geometry: { type: 'spherical', r0: 0 },
    contacts: { left: { phi: 'neutral' }, right: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' } },
    ports: [{ name: 'wire', region: 'p', from: a - 2e-9, to: a, V: 3.9, terminal: 'e-', species: { 'e-': 'equilibrium' } }],
    grid: { hmin: 1e-8, hmax: 1e-7 },
  });
  assert.equal(dev.model.ports[0].nodes.length, 1, 'the surface node alone');
  const before = dev.solve().conservation.find((c) => c.species === 'Li+').amount;
  dev.set({ ports: { wire: { I } } });
  const run = dev.advance(1, { dtMax: 0.2 });
  assert.ok(run.converged && run.done);
  const inserted = run.conservation.find((c) => c.species === 'Li+').amount - before;
  assert.ok(Math.abs(inserted / (-I / FARADAY) - 1) < 1e-6, `${inserted} mol inserted, ${-I / FARADAY} driven`);
});
