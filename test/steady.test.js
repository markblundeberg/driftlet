import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, AVOGADRO, units } from '../src/index.js';
import { build, layer, ohmic, bath, half, aqueous, metal, semiconductor, photogeneration } from '../src/kit.js';

const VT = (GAS_CONSTANT * 298.15) / FARADAY;

// Ag | AgNO₃ | Ag, strictly neutral, with the nitrate blocked at both electrodes. In steady state
// the nitrate carries no flux, so c ∝ e^{φ̂} and the silver flux N = −2D∇c makes c linear about
// its mean c₀, which only the conserved nitrate amount fixes. The electrode voltage is then
// V = 2V_T ln(c_L/c_R), so I = (4 F D c₀ / L) tanh(V / 4V_T).
test('steady state with a conserved spectator, solved directly: I = (4FDc₀/L) tanh(V/4V_T)', () => {
  const D = 1.65e-9, c0 = 10, L = 20e-6;
  const electrode = (V) => ({ V, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' });
  const cell = (cells) =>
    new Device({
      species: [
        { name: 'Ag+', z: 1, cRef: 1000 },
        { name: 'NO3-', z: -1, cRef: 1000 },
      ],
      materials: { water: { epsr: 0, species: { 'Ag+': { D, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
      regions: [{ material: 'water', length: L, c0: { 'NO3-': c0 } }],
      contacts: { left: electrode(0), right: electrode(0) },
      grid: { minCells: cells },
    });
  const error = (dev, V) => {
    dev.set({ contacts: { right: { V } } });
    const sol = dev.solve();
    assert.ok(sol.converged);
    assert.equal(sol.steps, 1, 'one steady (dt = ∞) solve, with the nitrate amount as a constraint');
    for (const st of sol.conservation) if (st.spectator) assert.ok(Math.abs(st.drift) < 1e-13, `drift ${st.drift}`);
    return Math.abs(sol.current) / (((4 * FARADAY * D * c0) / L) * Math.tanh(V / (4 * VT))) - 1;
  };
  const coarse = cell(400), fine = cell(800);
  for (const V of [0.01, 0.1]) assert.ok(Math.abs(error(coarse, V)) < 1e-5, `V=${V}`);
  // Near depletion (c_R ≈ 0.006 c₀ at 0.3 V) the error is larger, and second order.
  const e1 = error(coarse, 0.3), e2 = error(fine, 0.3);
  assert.ok(Math.abs(e1) < 5e-4 && e1 / e2 > 3.6 && e1 / e2 < 4.4, `${e1}, ${e2}`);
});

test('set() that closes a stretch conserves what it holds then, not what it held at the start', () => {
  // O₂ diffuses into a film for a while; then its contact is blocked. The film keeps that O₂.
  const RT = 8.314462618 * 298.15, L = 10e-6;
  const open = (c) => ({ species: { O2: { type: 'equilibrium', mu: RT * Math.log(c / 1000) } }, phi: 'neutral' });
  const def = (left) => ({
    species: [{ name: 'O2', z: 0, cRef: 1000 }],
    materials: { polymer: { epsr: 3, species: { O2: { D: 1e-11, mu0: 0 } } } },
    regions: [{ material: 'polymer', length: L, c0: { O2: 1 } }],
    contacts: { left, right: { phi: 'neutral' } },
    grid: { minCells: 50 },
  });
  const dev = new Device(def(open(1)));
  dev.solve();
  dev.set({ contacts: { left: open(2) } });
  const held = dev.advance(1).conservation[0].amount; // partly filled toward c = 2
  assert.ok(held > 1.3e-5 && held < 1.4e-5);
  dev.set({ contacts: { left: { species: { O2: 'blocked' } } } });
  const closed = dev.solve();
  assert.ok(closed.converged);
  assert.ok(Math.abs(closed.conservation[0].amount / held - 1) < 1e-12, `${closed.conservation[0].amount} vs ${held}`);
  assert.ok(Math.abs(closed.conservation[0].drift) < 1e-12);
});

// Slow ions sharing a layer with charged, reacting, immobile traps (SRH through explicit trap
// states, X⁰ + e⁻ = X⁻, X⁻ + h⁺ = X⁰): a cold start needs steps far below anything the ions'
// time scale suggests. (From a perovskite cell: iodide vacancies, D = 1e-17 m²/s, in MAPbI₃.)
test('a cold start with slow ions and immobile traps: solve() and advance() both reach equilibrium', () => {
  const F = FARADAY, T = 298, RT = GAS_CONSTANT * T, m3 = (n) => n / AVOGADRO;
  const gc = 8.1e24, gv = 5.8e24, Eg = 1.7, N0 = 1.6e25, Nt = 1e20;
  const ni = Math.sqrt(gc * gv) * Math.exp(-(Eg * F) / (2 * RT));
  const def = build({
    T,
    species: [{ name: 'e-', z: -1 }, { name: 'h+', z: 1 }, { name: 'X0', z: 0 }, { name: 'X-', z: -1 }, { name: 'V+', z: 1 }],
    materials: {
      MAPI: {
        epsr: 24.1,
        species: {
          'e-': { D: 1.7e-4, mu0: 0, cRef: m3(gc) },
          'h+': { D: 1.7e-4, mu0: units.eV(Eg), cRef: m3(gv) },
          X0: { D: 0, mu0: 0, cRef: m3(Nt) },
          'X-': { D: 0, mu0: RT * Math.log(ni / gc), cRef: m3(Nt) },
          'V+': { D: 1e-17, mu0: 0, cRef: m3(N0) },
        },
      },
    },
    stack: [ohmic(0), layer('MAPI', 400e-9, { fixedCharge: -F * m3(N0), c0: { X0: m3(Nt) / 2, 'X-': m3(Nt) / 2, 'V+': m3(N0) } }), ohmic(0)],
    bulkReactions: [
      { equation: 'e- + X0 = X-', kf: { MAPI: 1 / (3e-9 * m3(Nt)) } },
      { equation: 'X- + h+ = X0', kf: { MAPI: 1 / (3e-7 * m3(Nt)) } },
    ],
    grid: { hmin: 0.05e-9, hmax: 10e-9, ratio: 1.15 },
  });
  const flat = (s) => {
    for (const name of ['e-', 'h+', 'X-']) {
      const mu = s.mu[name];
      if ((Math.max(...mu) - Math.min(...mu)) / RT > 1e-6) return `${name} not flat`;
    }
    return '';
  };
  const s = new Device(def).solve();
  assert.ok(s.converged, 'cold steady solve');
  assert.equal(flat(s), '');
  const dev = new Device(def);
  const r = dev.advance(1e5, { tol: 1e-4 }); // from the cold start, with the default first step
  assert.ok(r.done && r.converged, `advance: done ${r.done}, ${r.steps} steps`);
  assert.equal(flat(dev.solution()), '');
});

test('a level flat across a face: holes held by one contact through two regions, blocked at the other', () => {
  // A Schottky contact that takes only electrons, on intrinsic GaAs, then heavily doped p-GaAs to
  // an ohmic contact: the holes reach that contact alone, so a steady state pins their level flat
  // there, through the face between the two layers, whose hole flux is then zero (an unknown
  // that the edge balances alone had set, and the solve was singular).
  const def = build({
    T: 300,
    library: [semiconductor('GaAs')],
    stack: [ohmic(0), layer('GaAs', 1.2e-6, { name: 'i' }), layer('GaAs', 1.1e-7, { name: 'p', acceptors: units.perCm3(7e18) }), ohmic(0)],
    grid: { hmin: 2e-10, hmax: 1.8e-8, ratio: 1.16 },
  });
  def.contacts.left = { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: { type: 'capacitive', C: 10 }, zeroCharge: 0.25 };
  const s = new Device(def).solve();
  assert.ok(s.converged);
  for (const name of ['e-', 'h+']) {
    const mu = Array.from(s.mu[name]).filter(Number.isFinite);
    assert.ok((Math.max(...mu) - Math.min(...mu)) / (GAS_CONSTANT * 300) < 1e-9, `${name} flat`);
  }
});

test('cold steady solves that need continuation: open circuit, a lit cell at forward bias, a current far from the start', () => {
  // Open circuit: the floating terminal held at a voltage, marched to where its current crosses
  // zero, and floated from there; the same V_oc on two grids.
  const Si = semiconductor('Si'), NA = units.perCm3(1e16), W = units.um(60.5);
  const cell = (hmin) =>
    build({
      T: 300,
      library: [Si],
      stack: [ohmic(0), layer('Si', units.um(0.5), { donors: units.perCm3(1e19) }), layer('Si', units.um(60), { acceptors: NA }), { ...ohmic(0), V: undefined, I: 0 }],
      bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1 / (60e-6 * NA) } }, photogeneration({ material: 'Si', flux: 3e-5, alpha: 1630, mu: units.eV(1.18), to: W })],
      grid: { hmin, hmax: units.um(1) },
    });
  const voc = [units.nm(1), units.nm(0.25)].map((h) => {
    const s = new Device(cell(h)).solve();
    assert.ok(s.converged, `hmin ${h}`);
    return s.terminals.right.V;
  });
  assert.ok(voc[0] > 0.3 && Math.abs(voc[1] - voc[0]) < 1e-5, `V_oc ${voc}`);

  // A lit GaAs diode with a Schottky contact, at 0.78 V forward from cold: the voltage
  // continuation starts with the terminals level, which lit needs the light ramped up too.
  const diode = () => {
    const def = build({
      T: 300,
      library: [semiconductor('GaAs')],
      stack: [ohmic(0), layer('GaAs', 1.73e-7, { acceptors: units.perCm3(3.3e18) }), layer('GaAs', 1.98e-7, { donors: units.perCm3(1.6e17) }), ohmic(0)],
      bulkReactions: [
        { equation: 'e- + h+ = 0', srh: { GaAs: { tauN: 3.2e-7, tauP: 3.2e-7 } } },
        photogeneration({ material: 'GaAs', flux: 2e-4, alpha: 1e6, mu: units.eV(1.72), to: 3.71e-7 }),
      ],
      grid: { hmin: 3e-10, hmax: 7.7e-9, ratio: 1.13 },
    });
    def.contacts.left = { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: { type: 'capacitive', C: 10 }, zeroCharge: 0.67 };
    return new Device(def);
  };
  const cold = diode().set({ contacts: { right: { V: 0.78 } } }).solve();
  const warm = diode();
  warm.solve();
  for (let V = 0.1; V < 0.785; V += 0.1) warm.set({ contacts: { right: { V: Math.min(V, 0.78) } } }).solve();
  const w = warm.set({ contacts: { right: { V: 0.78 } } }).solve();
  assert.ok(cold.converged && w.converged && Math.abs(cold.current - w.current) <= 1e-6 * Math.abs(w.current), `${cold.current} vs ${w.current}`);

  // Fe³⁺/Fe²⁺ on platinum against a bath, driven at a reducing current that needs about 3 V:
  // the march goes up, the way that raises the current into a passive device.
  const iron = new Device(build({
    library: [aqueous(['Fe3+', 'Fe2+', 'SO42-'], { epsr: 0 }), metal('Pt')],
    stack: [
      ohmic(0, ['e-']),
      layer('Pt', 1e-6),
      { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1.4e-5, alpha: 0.33 }] },
      layer('water', 1.97e-5, { c0: { 'Fe3+': 27.66, 'Fe2+': 0.435, 'SO42-': 41.925 } }),
      bath({ 'Fe3+': 27.66, 'Fe2+': 0.435, 'SO42-': 41.925 }, 'SO42-', { I: 0.89 }),
    ],
    grid: { hmin: 1.3e-8, hmax: 4.2e-7, ratio: 1.17 },
  }));
  const r = iron.solve();
  assert.ok(r.converged && Math.abs(r.terminals.right.current / 0.89 - 1) < 1e-6 && r.terminals.right.V > 2, `${r.terminals.right.current} A/m² at ${r.terminals.right.V} V`);
  // Oxidizing beyond the scarce Fe²⁺'s limiting current (~1.5 A/m²): no voltage gets there, and
  // the warning says so.
  const over = iron.set({ contacts: { right: { I: -20 } } }).solve();
  assert.ok(!over.converged && over.warnings.some((w) => w.startsWith('right: no steady state passes the driven current -20.0 A/m²')), over.warnings.join('\n'));
  // (and as the opposite sign passes, it says how a terminal's current is signed)
  assert.ok(over.warnings.some((w) => w.endsWith("(It passes 20.0 A/m²: if that was meant, a terminal's current is into the device.)")), over.warnings.join('\n'));
});

test('a closed Fe³⁺/Fe²⁺ cell driven by a current near its limit: the held voltage that passes it', () => {
  // Pt | Fe³⁺, Fe²⁺, Cl⁻ | Pt with resolved double layers, closed, at 70% of what passes at −1 V:
  // Fe³⁺ falls ~20 orders below Fe²⁺ at the cathode, and floated, the system loses ~33 digits. A
  // current-driven steady state is the held one at the voltage that passes the current, so that
  // is what the solve finds where the floated one can't.
  const face = { phi: { type: 'capacitive', C: 0.35 }, zeroCharge: -0.3, reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 0.07, alpha: 0.4 }] };
  const def = (right) => ({
    T: 300,
    species: [{ name: 'Fe3+', z: 3, cRef: 1000 }, { name: 'Fe2+', z: 2, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }],
    materials: {
      water: { epsr: 78.3, species: { 'Fe3+': { D: 6.04e-10, mu0: -4700 }, 'Fe2+': { D: 7.19e-10, mu0: -78900 }, 'Cl-': { D: 2.032e-9, mu0: -131228 } } },
      Pt: { conductor: { species: 'e-', conductivity: 9.5e6 } },
    },
    regions: [{ material: 'Pt', length: 1e-6 }, { material: 'water', length: 1.37e-6, c0: { 'Fe3+': 0.2, 'Fe2+': 15.5, 'Cl-': 31.6 } }, { material: 'Pt', length: 1e-6 }],
    interfaces: [face, face],
    contacts: { left: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk', ...right } },
    grid: { hmin: 4e-10, hmax: 2e-8, ratio: 1.15 },
  });
  const swept = new Device(def({ V: 0 }));
  swept.solve();
  let s;
  for (let j = 1; j <= 10; j++) s = swept.set({ contacts: { right: { V: -j / 10 } } }).solve();
  const I = 0.7 * s.terminals.right.current;
  const driven = new Device(def({ V: 0 }));
  driven.solve();
  const d = driven.set({ contacts: { right: { I } } }).solve();
  assert.ok(d.converged && Math.abs(d.terminals.right.current / I - 1) < 1e-9, `${d.converged}: ${d.terminals.right.current} vs ${I}`);
  const held = new Device(def({ V: d.terminals.right.V })).solve();
  assert.ok(held.converged && Math.abs(held.terminals.right.current / I - 1) < 1e-6, `held at ${d.terminals.right.V} V: ${held.terminals.right.current} vs ${I}`);
});

test('past the limiting current: a closed zinc cell, its blocked sulfate held flat, as a long transient finds it', () => {
  // Zn | ZnSO₄ | Zn, resolved double layers, at 1.1 i_lim (8FD₊c/L for a closed cell of a 2:2
  // salt): extended space charge at the cathode excludes the sulfate to ~1e-25 mol/m³. Its level
  // is flat at steady state, so steady solves hold it flat (the amount fixing where) rather than
  // find it through conductances that small, which had lost the system 14 digits.
  const face = { phi: { type: 'capacitive', C: 0.26 }, zeroCharge: -0.12, reactions: [{ ...half('Zn2+ + 2 e- = Zn(s)', { 'Zn(s)': 0 }), k0: 4e-3, alpha: 0.5 }] };
  const c = { 'Zn2+': 0.12, 'SO42-': 0.12 };
  const def = (V) =>
    build({
      T: 300,
      library: [aqueous(['Zn2+', 'SO42-'], { epsr: 78.3 }), metal('Zn')],
      stack: [ohmic(0, ['e-']), layer('Zn', 1e-6), face, layer('water', 4.9e-5, { c0: c }), face, layer('Zn', 1e-6), ohmic(V, ['e-'])],
      grid: { hmin: 1.8e-9, hmax: 1.4e-6, ratio: 1.15 },
    });
  const cold = new Device(def(0.7)).solve();
  const d = new Device(def(0));
  d.solve();
  const t = d.set({ contacts: { right: { V: 0.7 } } }).advance(100);
  const iLim = (8 * FARADAY * 7.03e-10 * 0.12) / 4.9e-5;
  assert.ok(cold.converged && t.done && Math.abs(cold.current / t.current - 1) < 1e-9, `${cold.current} vs ${t.current} A/m²`);
  assert.ok(Math.abs(cold.current) > 1.08 * iLim, `${cold.current} vs i_lim ${iLim}`);
});
