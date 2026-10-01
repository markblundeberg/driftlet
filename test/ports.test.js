import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT, units } from '../src/index.js';

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
  // minority-carrier diffusion from the back contact (weeks), and the inversion layer's Fermi
  // level is floating to within round-off: no steady state can be computed. A port holding the
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
        right: { V: -3, phi: { type: 'capacitive', C: Cg, zeroCharge: -0.9 } },
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

test('ports are checked', () => {
  const def = line(1e6, 10).def;
  def.ports[0].species['Li+'] = { type: 'conductance', G: 1, offset: 0, mu: 0 };
  assert.throws(() => new Device(def), (e) => e instanceof DeviceError && /not by mu/.test(e.message));
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
