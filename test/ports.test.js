import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, semiconductor, metal } from '../src/kit.js';

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

test('MOS: a port grounding the channel gives the low-frequency C–V, inversion included', () => {
  // p-Si under a gate (right). Without generation, inversion electrons could only arrive by
  // minority-carrier diffusion from the back contact (minutes, through a bulk with ~1e3 of them
  // per cm³), and the inversion layer's Fermi level is held so weakly that the steady solve is
  // fragile (see the next test, and the roadmap). A port holding the
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

test('MOS without a port: the impedance gives the high-frequency C–V at 1 Hz, as a time-domain run does', () => {
  // The inversion layer's electrons can only come by minority diffusion from the back contact,
  // through a bulk with ~1e3 cm⁻³ of them: minutes. The assembled Jacobian alone (two huge entries
  // per inversion-layer flux, nearly cancelling) gave the layer an exchange path that followed
  // the gate at 1 Hz; the impedance now takes J·v from the residual, by GMRES.
  // (Its steady solve is itself fragile, the roadmap's weakly held minority: it converges on
  // this grid, built from units' values exactly, and not on some within round-off of it.)
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
