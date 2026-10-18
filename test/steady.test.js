import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, AVOGADRO, units } from '../src/index.js';
import { build, layer, ohmic, bath, half, aqueous, metal, semiconductor, photogeneration, check } from '../src/kit.js';

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
  // Driven by a waveform, the steady state is at the waveform's value now: the pseudo-transient
  // ramp's steps don't move the time it's read at (they had: the solve came out at 1 V).
  const ramped = new Device({ ...def, contacts: { ...def.contacts, right: { ...def.contacts.right, V: { t: [0, 1], values: [0, 1] } } } }).solve();
  assert.ok(ramped.converged && ramped.steps > 1 && ramped.terminals.right.V === 0, `ramped: V ${ramped.terminals.right.V}, ${ramped.steps} steps`);
  assert.equal(flat(ramped), '');
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

test('a closed battery at open circuit keeps its charge: a cold solve at the OCV of its c0, and after a charge, where the transient relaxes', () => {
  // Two insertion hosts and an electrolyte, the right collector at I = 0: its electrons change
  // only through that terminal, so their amount is conserved, and a steady solve keeps it (one
  // steady state per state of charge). The start puts each host's Li⁺ level with the
  // electrolyte's (its φ is a gauge), so the terminal starts at the OCV too.
  const Ea = (x) => 0.09 + 0.55 * Math.exp(-x / 0.12) + 0.12 * (1 - x) - 0.012 * Math.log(x / (1 - x));
  const Ec = (x) => 4.25 - 0.55 * x - 0.06 * Math.log(x / (1 - x));
  const xs = Array.from({ length: 201 }, (_, i) => 0.005 + (0.99 * i) / 200);
  const host = (cMax, E) => ({
    epsr: 0,
    species: { 'Li+': { D: 5e-14, mu0: 0, cRef: cMax }, 'e-': { D: 1e-4, mu0: 0, cRef: cMax } },
    statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax, ocv: { x: xs, E: xs.map(E), muRef: 0 } }],
  });
  const face = { species: { 'Li+': 'blocked' }, reactions: [{ equation: 'Li+(left) = Li+(right)', k0: 1e-4 }] };
  const collector = (drive) => ({ ...drive, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const cell = (right) =>
    new Device({
      species: [{ name: 'Li+', z: 1, cRef: 1000 }, { name: 'e-', z: -1, cRef: 1000 }, { name: 'A-', z: -1, cRef: 1000 }],
      materials: { anode: host(30000, Ea), cathode: host(25000, Ec), elyte: { epsr: 0, species: { 'Li+': { D: 1.5e-10, mu0: 0 }, 'A-': { D: 2.5e-10, mu0: 0 } } } },
      regions: [
        { name: 'anode', material: 'anode', length: 1.5e-6, c0: { 'Li+': 0.6 * 30000 } },
        { name: 'sep', material: 'elyte', length: 40e-6, c0: { 'Li+': 1000, 'A-': 1000 } },
        { name: 'cathode', material: 'cathode', length: 2e-6, c0: { 'Li+': 0.35 * 25000 } },
      ],
      interfaces: [face, face],
      contacts: { left: collector({ V: 0 }), right: collector(right) },
      grid: { hmin: 5e-9, hmax: 50e-9 },
    });
  const ocv = Ec(0.35) - Ea(0.6);
  const dev = cell({ I: 0 });
  assert.ok(Math.abs(dev.solver.termV[1] - ocv) < 1e-7, `starts at ${dev.solver.termV[1]} V, OCV ${ocv} V`);
  const cold = dev.solve();
  assert.ok(cold.converged && Math.abs(cold.terminals.right.V - ocv) < 1e-7, `${cold.terminals.right.V} V, OCV ${ocv} V`);
  // Charged at a current for a minute, then left open: the transient relaxes to where the steady
  // solve goes.
  dev.set({ contacts: { right: { I: 5 } } });
  assert.ok(dev.advance(60).converged);
  dev.set({ contacts: { right: { I: 0 } } });
  const relaxed = dev.advance(3600);
  const steady = dev.solve();
  assert.ok(steady.converged && Math.abs(steady.terminals.right.V - relaxed.terminals.right.V) < 1e-6, `${steady.terminals.right.V} vs ${relaxed.terminals.right.V} V`);
  assert.ok(steady.terminals.right.V > ocv + 0.05, 'charged');
  // Driven at a current, it only stores what comes in: no steady state, and the solve says so.
  const driven = cell({ I: 2 }).solve();
  assert.ok(!driven.converged && driven.warnings.some((w) => w.startsWith('right: driven at 2.00 A/m², it has no steady state: what it feeds (e- in cathode)')), driven.warnings.join('\n'));
});

test('a closed population with immobile traps: their electrons conserved with it, each node\'s traps on their own, solved directly as a transient ends', () => {
  // An insulating film between gates, its electrons mobile and closed, trap states X⁰/X⁻ fixed in
  // place (e⁻ + X⁰ = X⁻): e⁻ + X⁻ is conserved over the film, X⁰ + X⁻ at every node.
  const def = (V) => ({
    species: [{ name: 'e-', z: -1 }, { name: 'X0', z: 0 }, { name: 'X-', z: -1 }],
    materials: { film: { epsr: 4, species: { 'e-': { D: 1e-6, mu0: 0, cRef: 1 }, X0: { D: 0, mu0: 0, cRef: 1 }, 'X-': { D: 0, mu0: -3000, cRef: 1 } } } },
    regions: [{ name: 'film', material: 'film', length: 50e-9, fixedCharge: 2 * FARADAY, c0: { 'e-': 1, X0: 1, 'X-': 1 } }],
    bulkReactions: [{ equation: 'e- + X0 = X-', kf: { film: 1e3 } }],
    contacts: { left: { V: 0, phi: { type: 'capacitive', C: 0.1 }, zeroCharge: 0 }, right: { V, phi: { type: 'capacitive', C: 0.1 }, zeroCharge: 0 } },
    grid: { minCells: 60 },
  });
  const steady = new Device(def(0.05)).solve();
  assert.ok(steady.converged && steady.steps === 1, `steps ${steady.steps}`);
  const amount = (s, name) => s.conservation.find((c) => c.species === name).amount;
  assert.ok(Math.abs((amount(steady, 'e-') + amount(steady, 'X-')) / 1e-7 - 1) < 1e-12);
  for (let g = 0; g < steady.x.length; g++) assert.ok(Math.abs(steady.c.X0[g] + steady.c['X-'][g] - 2) < 1e-12, `node ${g}`);
  const run = new Device(def(0.05)).advance(1);
  for (const name of ['e-', 'X-']) assert.ok(Math.abs(amount(run, name) / amount(steady, name) - 1) < 1e-6, `${name}: ${amount(run, name)} vs ${amount(steady, name)}`);
});

// A closed solution with two complexations sharing A (A + B = AB, A + C = AC): three amounts
// conserved together (B + AB, C + AC, A + AB + AC), each in place of one balance row, the three
// independent rows whatever order the species are listed in. At steady state it's homogeneous
// chemical equilibrium: c_AB = K₁ c_A c_B / c_ref, c_AC = K₂ c_A c_C / c_ref.
test('two complexations sharing a species: the amounts in place of independent rows, in any species order', () => {
  const D = { A: 1e-9, B: 1.3e-9, C: 0.9e-9, AB: 0.8e-9, AC: 0.7e-9, 'Na+': 1.3e-9, 'Cl-': 2e-9 };
  const mu0 = { A: 0, B: 0, C: 0, AB: -5e3, AC: -4e3, 'Na+': -262e3, 'Cl-': -131e3 };
  const c0 = { A: 20, B: 5, C: 7, AB: 30, AC: 12, 'Na+': 30, 'Cl-': 30 };
  const solve = (order) =>
    new Device({
      species: order.map((name) => ({ name, z: name === 'Na+' ? 1 : name === 'Cl-' ? -1 : 0, cRef: 1000 })),
      materials: { w: { epsr: 0, species: Object.fromEntries(order.map((n) => [n, { D: D[n], mu0: mu0[n] }])) } },
      regions: [{ material: 'w', length: 40e-6, c0 }],
      bulkReactions: [{ equation: 'A + B = AB', kf: { w: 1e-5 } }, { equation: 'A + C = AC', kf: { w: 1e-5 } }],
      contacts: { left: { bath: { c: { 'Na+': 30, 'Cl-': 30 }, reference: 'Cl-' }, V: 0 }, right: { phi: 'neutral' } },
      grid: { hmin: 0.5e-6, hmax: 3e-6 },
    }).solve();
  const RT = GAS_CONSTANT * 298.15, K1 = Math.exp(5e3 / RT), K2 = Math.exp(4e3 / RT);
  const [B, C, A] = [c0.B + c0.AB, c0.C + c0.AC, c0.A + c0.AB + c0.AC];
  let lo = 0, hi = A;
  for (let k = 0; k < 200; k++) {
    const a = 0.5 * (lo + hi);
    if (a + (B * K1 * a) / (1000 + K1 * a) + (C * K2 * a) / (1000 + K2 * a) > A) hi = a;
    else lo = a;
  }
  const cA = 0.5 * (lo + hi);
  const permutations = (a) => (a.length <= 1 ? [a] : a.flatMap((x, j) => permutations([...a.slice(0, j), ...a.slice(j + 1)]).map((p) => [x, ...p])));
  for (const order of permutations(['A', 'B', 'C', 'AB', 'AC'])) {
    const sol = solve([...order, 'Na+', 'Cl-']);
    assert.ok(sol.converged && sol.steps === 1, order.join());
    for (const g of [0, sol.x.length - 1]) assert.ok(Math.abs(sol.c.A[g] / cA - 1) < 1e-9, `${order}: c_A ${sol.c.A[g]} vs ${cA}`);
  }
});

test('electrons held flat behind an open circuit: their amount stands in for the circuit row, and the solve goes direct', () => {
  // n-Si between a gate and a contact that passes only electrons, at I = 0: the electrons carry
  // nothing at steady state (flat at the contact's level) and nothing but that contact feeds
  // them, so their amount is what fixes the level. (Pinned as a balance row beside the flat rows,
  // the circuit row was left all zeros, and the solve failed.)
  const Si = semiconductor('Si');
  const def = {
    T: 300,
    species: Si.species,
    materials: { Si: Si.materials.Si },
    regions: [{ material: 'Si', length: 1e-6, fixedCharge: units.perCm3(1e16) * FARADAY, c0: { 'h+': units.perCm3(1e4) } }],
    contacts: { left: { V: 0, phi: { type: 'capacitive', C: 0.05 }, zeroCharge: 0 }, right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
    grid: { hmin: 1e-9, hmax: 2e-8 },
  };
  const steady = new Device(def).solve();
  const run = new Device(def).advance(1);
  assert.ok(steady.converged && steady.steps === 1, `steps ${steady.steps}`);
  assert.ok(Math.abs(steady.terminals.right.V - run.terminals.right.V) < 1e-9, `${steady.terminals.right.V} vs ${run.terminals.right.V} V`);
});

test('a MOS capacitor\'s back contact moved under an inverted gate: the electrons\' level follows it, flat', () => {
  // Without recombination the inversion layer's electrons reach the back contact only through a
  // bulk with ~1e3 of them per cm³, so found through their own conduction the level can't move
  // (each solve failed through ~90 pseudo-transient steps); a steady state has them flat at the
  // contact's level, and steady solves hold them there.
  const Si = semiconductor('Si');
  const dev = new Device(build({
    T: 300,
    library: [Si],
    stack: [ohmic(0), layer('Si', 2e-6, { acceptors: units.perCm3(1e17) }), { V: 1.5, phi: { type: 'capacitive', C: 3.45e-3 }, zeroCharge: 0 }],
    grid: { hmin: 1e-10, hmax: 5e-8 },
  }));
  assert.ok(dev.solve().converged);
  for (const Vb of [0.1, -0.3]) {
    const s = dev.set({ contacts: { left: { V: Vb } } }).solve();
    assert.ok(s.converged && s.iterations < 20, `${Vb} V: ${s.iterations} iterations`);
    const mu = Array.from(s.mu['e-']);
    assert.ok(Math.abs(Math.max(...mu) - Math.min(...mu)) < 1e-9 * GAS_CONSTANT * 300 && Math.abs(mu[0] + FARADAY * Vb) < 1e-6, `${Vb} V`);
  }
});

test('a level held flat across a face: the zinc cell past its limit with its electrolyte in two regions', () => {
  // As above, but the sulfate crosses a face (a separator of the same water). Found through its
  // own conduction there, it took twice the iterations short of the limit and failed past it.
  const face = { phi: { type: 'capacitive', C: 0.26 }, zeroCharge: -0.12, reactions: [{ ...half('Zn2+ + 2 e- = Zn(s)', { 'Zn(s)': 0 }), k0: 4e-3, alpha: 0.5 }] };
  const c = { 'Zn2+': 0.12, 'SO42-': 0.12 };
  const cell = (split) =>
    new Device(build({
      T: 300,
      library: [aqueous(['Zn2+', 'SO42-'], { epsr: 78.3 }), metal('Zn')],
      stack: [ohmic(0, ['e-']), layer('Zn', 1e-6), face, ...(split ? [layer('water', 2.45e-5, { c0: c }), { dipole: 0 }, layer('water', 2.45e-5, { c0: c })] : [layer('water', 4.9e-5, { c0: c })]), face, layer('Zn', 1e-6), ohmic(0.9, ['e-'])],
      grid: { hmin: 1.8e-9, hmax: 1.4e-6, ratio: 1.15 },
    })).solve();
  const one = cell(false), two = cell(true);
  assert.ok(one.converged && two.converged && Math.abs(two.current / one.current - 1) < 1e-3, `${two.current} vs ${one.current} A/m²`);
});

test('a MOS capacitor whose carriers recombine, its back contact moved: still at equilibrium, its gate charge the gate\'s moved the other way', () => {
  // GaAs (~1e-4 minority electrons per cm³ in the bulk) with slow SRH recombination, under a gate
  // in inversion. The electrons and holes react only with each other and reach only the back
  // contact, so a steady state is equilibrium: both flat at its levels. Found through generation
  // and their own conduction, moving that contact by 70 mV failed to converge, or converged to
  // a gate charge 2× off.
  const def = (Vb, Vg) => build({
    T: 300,
    library: [semiconductor('GaAs'), metal('Pt')],
    materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [
      ohmic(Vb),
      layer('GaAs', 9.4e-7, { name: 'bulk', acceptors: units.perCm3(2.4e16) }),
      { dipole: 0 },
      layer('SiO2', 1e-9, { grid: { hmin: 1e-10, hmax: 2.5e-10 } }),
      { phi: { type: 'capacitive', C: 100 }, zeroCharge: -0.51 },
      layer('Pt', 5e-8),
      ohmic(Vg, ['e-']),
    ],
    bulkReactions: [{ equation: 'e- + h+ = 0', srh: { GaAs: { tauN: 8e-9, tauP: 8e-9 } } }],
    grid: { hmin: 2.8e-10, hmax: 1.7e-8, ratio: 1.1 },
  });
  const dev = new Device(def(0, 1.19));
  assert.ok(dev.solve().converged);
  const moved = dev.set({ contacts: { left: { V: 0.07 } } }).solve();
  const cold = new Device(def(0, 1.12)).solve();
  assert.ok(moved.converged && cold.converged, `${moved.converged}, ${cold.converged}`);
  const D = (s) => s.interfaces.at(-1).D;
  assert.ok(Math.abs(D(moved) / D(cold) - 1) < 1e-6, `${D(moved)} vs ${D(cold)} C/m²`);
  const bulk = (s, name) => Array.from(s.mu[name]).slice(0, dev.grid.regionEnd[0] + 1);
  for (const name of ['e-', 'h+']) {
    const mu = bulk(moved, name);
    assert.ok(Math.max(...mu) - Math.min(...mu) < 1e-9 * GAS_CONSTANT * 300, `${name} flat`);
  }
});

test('a floating gate, its back contact jumped far: the holes held flat carry exactly nothing, and the terminal currents add up', () => {
  // The p bulk's holes reach the back contact alone, so a steady solve holds their level flat at
  // it; flat only to round-off after a jump to −6 V, they had read as 3e-4 A/m² at that contact,
  // through a bulk conducting ~1e10 S/m², against the 3.5e-7 A/m² the inversion layer leaks.
  const oxide = { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } };
  const dev = new Device(build({
    T: 300,
    library: [semiconductor('Si'), metal('Al'), oxide],
    stack: [ohmic(1, ['e-']), layer('Al', 20e-9), { phi: { type: 'capacitive', C: 10 }, zeroCharge: -0.95 }, layer('SiO2', 10e-9), { dipole: 0 },
      layer('Si', 1e-6, { name: 'Si', acceptors: units.perCm3(1e17) }), ohmic(0)],
    ports: [{ name: 'channel', region: 'Si', from: 0, to: 1e-9, V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' } }],
    grid: { hmin: 0.25e-9, hmax: 50e-9, ratio: 1.15 },
  }));
  dev.solve();
  dev.set({ contacts: { left: { I: 0 } } });
  const Q = dev.solve().interfaces[0].D;
  dev.set({ contacts: { right: { V: -6 } } });
  const s = dev.solve(), { left, right, channel } = s.terminals;
  assert.ok(s.converged && Math.abs(s.interfaces[0].D / Q - 1) < 1e-12, `gate charge ${s.interfaces[0].D} vs ${Q}`);
  assert.ok(channel.current > 1e-8 && Math.abs(left.current + right.current + channel.current) < 1e-9 * channel.current, `currents ${left.current}, ${right.current}, ${channel.current}`);
});

// An electrode port at no current reducing M⁺ + e⁻ = M(s), M⁺ closed in a film whose K⁺ and X⁻
// a bath holds. What the reactions take, the electrode's charge gives: M⁺ + Q/F is conserved
// (with no capacitance, Q = 0, the M⁺ alone), which a steady solve keeps in place of the port's
// circuit row. Steady, M⁺ is uniform and at Nernst with the electrode, X⁻ − K⁺ = c(M⁺) with
// K⁺X⁻ = c_b²: V = (μ°(M⁺) + μ°(X⁻) + RT ln(c(M⁺) c(X⁻)/c_ref²))/F, the bath's X⁻ at 0 V.
// Silver and iron couples on one platinum electrode at open circuit, in a closed solution: at
// steady state both are at equilibrium with the metal's electrons, so with each other,
// Fe³⁺ + Ag(s) = Fe²⁺ + Ag⁺ at K = exp(−ΔG°/RT), however slow either is. The silver couple's rate
// was the only row holding Ag⁺'s level, and at k0 = 1e-20 it was lost to diffusion's round-off:
// the solve settled where silver wasn't (−0.4661 V), and 1e-15 to 1e-18 failed.
test('a reaction 1e27 times slower than the rest still reaches its equilibrium in a steady solve: Ag⁺/Ag beside Fe³⁺/Fe²⁺ at open circuit', () => {
  const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 }), iron = half('Fe3+ + e- = Fe2+');
  const ions = ['K+', 'NO3-', 'Ag+', 'Fe3+', 'Fe2+'], L = 50e-6;
  const cell = (kAg, order) => new Device(build({
    library: [aqueous(order, { epsr: 0 }), metal('Pt')],
    stack: [
      bath({ 'K+': 200, 'NO3-': 200 }, 'NO3-', 0),
      layer('water', L, { c0: { 'K+': 200, 'NO3-': 285, 'Ag+': 10, 'Fe3+': 5, 'Fe2+': 30 } }),
      { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [{ ...silver, k0: kAg }, { ...iron, k0: 1e-3 }] },
      layer('Pt', 1e-6),
      ohmic({ I: 0 }, ['e-']),
    ],
    grid: { hmin: 0.5e-9, hmax: 3e-6, ratio: 1.2 },
  }));
  const K = Math.exp(-(-78.9e3 + 77.107e3 + 4.7e3) / (VT * FARADAY)); // (μ° of Fe²⁺ and Ag⁺ less Fe³⁺'s; Ag(s) at 0)
  let V0;
  for (const kAg of [1e-3, 1e-15, 1e-30]) {
    for (const order of [ions, [...ions].reverse()]) {
      const s = cell(kAg, order).solve(), at = `k0 ${kAg}, ${order}`;
      assert.ok(s.converged, at);
      const g = s.x.findIndex((x) => x > L / 2), c = (sp) => s.c[sp][g];
      assert.ok(Math.abs((c('Fe2+') * c('Ag+')) / (c('Fe3+') * 1000) / K - 1) < 1e-9, `${at}: K ${(c('Fe2+') * c('Ag+')) / (c('Fe3+') * 1000)} for ${K}`);
      const fe = s.conservation.filter((x) => x.species === 'Fe3+' || x.species === 'Fe2+').reduce((t, x) => t + x.amount, 0);
      assert.ok(Math.abs(fe / (35 * L) - 1) < 1e-9, `${at}: iron ${fe}`);
      V0 ??= s.terminals.right.V;
      assert.ok(Math.abs(s.terminals.right.V - V0) < 1e-9, `${at}: V ${s.terminals.right.V} for ${V0}`);
    }
  }
  // Ag⁺ complexed too (Ag⁺ + NO₃⁻ = AgNO₃, fast): Ag⁺ is then in two reactions, and it's the
  // silver, Ag⁺ + AgNO₃, that the couple alone changes.
  const K2 = Math.exp(3e3 / (VT * FARADAY));
  let V1;
  for (const kAg of [1e-3, 1e-25]) {
    const d = build({
      library: [aqueous(ions, { epsr: 0 }), metal('Pt')],
      stack: [
        bath({ 'K+': 200, 'NO3-': 200 }, 'NO3-', 0),
        layer('water', L, { c0: { 'K+': 200, 'NO3-': 285, 'Ag+': 10, 'Fe3+': 5, 'Fe2+': 30 } }),
        { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [{ ...silver, k0: kAg }, { ...iron, k0: 1e-3 }] },
        layer('Pt', 1e-6),
        ohmic({ I: 0 }, ['e-']),
      ],
      grid: { hmin: 0.5e-9, hmax: 3e-6, ratio: 1.2 },
    });
    const water = d.materials.water.species;
    d.species.push({ name: 'AgNO3', z: 0, cRef: 1000 });
    water.AgNO3 = { D: 1e-9, mu0: water['Ag+'].mu0 + water['NO3-'].mu0 - 3e3 };
    d.regions.find((r) => r.c0?.['Ag+']).c0.AgNO3 = 1;
    d.bulkReactions = [{ equation: 'Ag+ + NO3- = AgNO3', kf: { water: 1e3 } }];
    const s = new Device(d).solve();
    assert.ok(s.converged, `complexed, k0 ${kAg}`);
    const g = s.x.findIndex((x) => x > L / 2), c = (sp) => s.c[sp][g];
    assert.ok(Math.abs((c('Fe2+') * c('Ag+')) / (c('Fe3+') * 1000) / K - 1) < 1e-9, `complexed, k0 ${kAg}: K ${(c('Fe2+') * c('Ag+')) / (c('Fe3+') * 1000)}`);
    assert.ok(Math.abs((c('AgNO3') * 1000) / (c('Ag+') * c('NO3-')) / K2 - 1) < 1e-9, `complexed, k0 ${kAg}: K2`);
    V1 ??= s.terminals.right.V;
    assert.ok(Math.abs(s.terminals.right.V - V1) < 1e-9, `complexed, k0 ${kAg}: V ${s.terminals.right.V} for ${V1}`);
  }
});

// An electrode at open circuit between a bath's Fe³⁺/Fe²⁺ and a closed pool of Ag⁺ (in a film the
// bath doesn't reach for silver): the pool has no sink but the electrode, so the silver couple
// settles at equilibrium, and with no current, the iron one too, at the bath's level, however
// slow it is. The current-driven electrode kept the silver couple off its affinity, and from an
// iron k0 of 1e-15 the solve converged 78 mV off, its iron not stationary.
test('an open-circuit electrode feeding a closed pool: both couples at equilibrium, however slow the other', () => {
  const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 }), iron = half('Fe3+ + e- = Fe2+');
  const bathC = { 'K+': 5000, 'Fe3+': 10, 'Fe2+': 10, 'NO3-': 5050 };
  const cell = (kFe, order) => new Device(build({
    library: [aqueous(order, { epsr: 0 }), metal('Pt')],
    stack: [
      ohmic({ I: 0 }, ['e-']),
      layer('Pt', 1e-6),
      { reactions: [{ ...iron, k0: kFe }, { ...silver, k0: 1e-3 }] },
      layer('water', 20e-6, { c0: { ...bathC, 'Ag+': 1, 'NO3-': 5051 } }),
      bath(bathC, 'NO3-'),
    ],
    grid: { hmin: 10e-9, hmax: 0.5e-6 },
  }));
  const K = Math.exp(-(-78.9e3 + 77.107e3 + 4.7e3) / (VT * FARADAY)); // Fe³⁺ + Ag(s) = Fe²⁺ + Ag⁺
  const ions = ['K+', 'NO3-', 'Fe3+', 'Fe2+', 'Ag+'];
  let V0;
  for (const kFe of [1e-3, 1e-20]) {
    for (const order of [ions, [...ions].reverse()]) {
      const s = cell(kFe, order).solve(), at = `k0 ${kFe}, ${order[0]} first`;
      assert.ok(s.converged && s.warnings.length === 0, `${at}: ${s.warnings.join(' ')}`);
      const g = s.c['Ag+'].findIndex((x) => x > 0), c = (sp) => s.c[sp][g]; // (the water's first node, at the face)
      assert.ok(Math.abs((c('Fe2+') * c('Ag+')) / (c('Fe3+') * 1000) / K - 1) < 1e-9, `${at}: K ${(c('Fe2+') * c('Ag+')) / (c('Fe3+') * 1000)} for ${K}`);
      V0 ??= s.terminals.left.V;
      assert.ok(Math.abs(s.terminals.left.V - V0) < 1e-9, `${at}: V ${s.terminals.left.V} for ${V0}`);
    }
  }
});

// A reaction a steady state needs at equilibrium, too slow for the solve to see, said so: a slow
// complexation in a cell whose only reservoir is a K⁺ contact kept its starting composition
// (K = 0.1 where it's 7.5) and converged without a word.
test('a reaction too slow for a steady solve to resolve is warned of, not passed off as steady', () => {
  const cell = (kf) => new Device({
    species: [{ name: 'K+', z: 1, cRef: 1000 }, { name: 'A-', z: -1, cRef: 1000 }, { name: 'KA', z: 0, cRef: 1000 }],
    materials: { w: { epsr: 0, species: { 'K+': { D: 1e-9, mu0: 0 }, 'A-': { D: 1e-9, mu0: 0 }, KA: { D: 1e-9, mu0: -5e3 } } } },
    regions: [{ material: 'w', length: 1e-4, c0: { 'K+': 100, 'A-': 100, KA: 1 } }],
    bulkReactions: [{ equation: 'K+ + A- = KA', kf: { w: kf } }],
    contacts: { left: { V: 0, terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, species: {}, phi: 'neutral' } },
  });
  const fast = cell(1).solve(), K = Math.exp(5e3 / (VT * FARADAY));
  assert.ok(fast.converged && fast.warnings.length === 0, fast.warnings.join(' '));
  assert.ok(Math.abs((fast.c.KA[5] * 1000) / (fast.c['K+'][5] * fast.c['A-'][5]) / K - 1) < 1e-9);
  const d = cell(1e-30), slow = d.solve();
  assert.ok(slow.warnings.some((w) => /bulkReactions\[0\]: out of equilibrium/.test(w)), slow.warnings.join(' '));
  // and check() doesn't call it equilibrium
  const bal = check(d, slow, { refine: false }).items.find((it) => it.name === 'balance');
  assert.ok(bal.ok === false && /not a steady state/.test(bal.summary), bal.summary);
});

// Found by round 11's fuzzer. Two metals about a closed film, the right one at open circuit with
// one couple written twice (once each way, k0 8.8e-5 and 504). A pivot cancelled to exactly
// nothing, the perturbed solve left one spelling's rate row 41.7 mol/(m²·s) off, and the small
// update passed for convergence: V 0.7293, where the couple's equilibrium is 0.6906 V.
test('a pivot perturbed in a steady solve: the residual is checked, not the update alone', () => {
  const water = { epsr: 0, species: { 'A+': { D: 1e-9, mu0: -20000 }, 'B-': { D: 1e-9, mu0: -12000 }, 'D2-': { D: 1e-9, mu0: 14000 }, 'F-': { D: 1e-9, mu0: 1000 } } };
  const cell = (couple, order) => {
    const species = [{ name: 'A+', z: 1, cRef: 1000 }, { name: 'B-', z: -1, cRef: 1000 }, { name: 'D2-', z: -2, cRef: 1000 }, { name: 'F-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }];
    return new Device({
      T: 298.15,
      species: order === 'reversed' ? species.reverse() : species,
      materials: { water, M: { conductor: { species: 'e-', conductivity: 47994048.15234138 } } },
      regions: [
        { name: 'M1', material: 'M', length: 1e-6 },
        { name: 'water', material: 'water', length: 4.811448493222383e-6, c0: { 'A+': 298.8472666344369, 'B-': 93.32160813739429, 'D2-': 88.61516488336346, 'F-': 28.295328730315656 } },
        { name: 'M2', material: 'M', length: 1e-6 },
      ],
      interfaces: [{ reactions: [{ equation: 'B- = A+ + 2 e-', k0: 360.6517926504172 }] }, { reactions: couple }],
      contacts: { left: { V: 0.8281985982321203, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
      grid: { hmin: 4.8114484932223834e-8, hmax: 3.207632328814922e-7 }, // (the fuzzer's own numbers: the pivot cancels exactly only at them)
    }).solve();
  };
  const slow = { equation: 'F- + e- = D2-', k0: 8.82494895894919e-5 }, fast = { equation: 'D2- = F- + e-', k0: 504.44001895513924 };
  const V = cell([fast]).terminals.right.V;
  for (const [couple, order] of [[[slow], 'as written'], [[slow, fast], 'as written'], [[slow, fast], 'reversed'], [[fast, slow], 'as written']]) {
    const s = cell(couple, order);
    assert.ok(s.converged && Math.abs(s.terminals.right.V - V) < 1e-9, `${couple.length} spellings, ${order}: ${s.converged}, V ${s.terminals.right.V} for ${V}`);
  }
});

test('a level only the face after it holds: the same steady state in any listing of the species', () => {
  // Strictly neutral water between two electrodes: its level is set only by the fast couple at its
  // right face (C²⁺ + 2e⁻ = S₂(s), the right electrode held), so eliminated from the left, its
  // last node's pivot is round-off alone, 7e-32 or 3e-16 as the species happen to be listed. At
  // 7e-32 it scaled the face's block by 1e31, which then lost its pivot outright, and the solve
  // failed. The left electrode, at no current, then sits where its own couple is in equilibrium.
  const water = { epsr: 0, species: { 'F-': { D: 1e-9, mu0: -18755 }, 'C2+': { D: 1e-9, mu0: -12698 } } };
  const species = [{ name: 'F-', z: -1, cRef: 1000 }, { name: 'C2+', z: 2, cRef: 1000 }, { name: 'e-', z: -1 }];
  const VR = 0.48157, RT = GAS_CONSTANT * 298.15;
  const phi = (-23475 + 12698 - RT * Math.log(100 / 1000) + 2 * FARADAY * VR) / (2 * FARADAY);
  const VL = (18230 + 18755 - RT * Math.log(200 / 1000) + FARADAY * phi) / FARADAY;
  for (const order of ['as written', 'reversed']) {
    const s = new Device({
      species: order === 'reversed' ? [...species].reverse() : species,
      materials: { water, M: { conductor: { species: 'e-', conductivity: 4e5 } } },
      regions: [{ name: 'M1', material: 'M', length: 1e-6 }, { name: 'water', material: 'water', length: 6.4e-6, c0: { 'F-': 200, 'C2+': 100 } }, { name: 'M2', material: 'M', length: 1e-6 }],
      interfaces: [
        { reactions: [{ equation: 'F- = e- + S(s)', k0: 2.3e-15, alpha: 0.7, fixed: { 'S(s)': 18230 } }] },
        { reactions: [{ equation: 'C2+ + 2 e- = S1(s)', k0: 6.1e-16, alpha: 1, fixed: { 'S1(s)': 355 } }, { equation: 'C2+ + 2 e- = S2(s)', k0: 59, alpha: 0.5, fixed: { 'S2(s)': -23475 } }] },
      ],
      contacts: { left: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: VR, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
      grid: { hmin: 6.4e-8, hmax: 4.3e-7 }, // (the fuzzer's own numbers)
    }).solve();
    assert.ok(s.converged && Math.abs(s.terminals.left.V - VL) < 1e-8, `${order}: ${s.converged}, V ${s.terminals.left.V} for ${VL}`);
  }
});

test('the equilibria of reactions take the same rows in any listing: two electrodes with solid couples, at equilibrium', () => {
  // Each couple with a fixed solid sets its ion's level against its electrode's, so at
  // equilibrium the voltages follow from the solids alone. Picked by where the species were
  // listed, the rows that the couples' affinities stand in for differed, and the solve converged
  // in one listing and failed in another.
  const solve = (def, order) => {
    if (order === 'species reversed') def.species.reverse();
    if (order === 'reactions reversed') for (const itf of def.interfaces) itf.reactions.reverse();
    return new Device(def).solve();
  };
  const metal = { conductor: { species: 'e-', conductivity: 1085645.496591289 } };
  const silver = (order) =>
    solve(
      {
        species: [['A+', 1], ['F-', -1], ['G3+', 3], ['N', 0], ['AG4+', 4], ['e-', -1]].map(([name, z]) => (name === 'e-' ? { name, z } : { name, z, cRef: 1000 })),
        materials: {
          water: { epsr: 0, species: { 'A+': { D: 2.3243172943985234e-10, mu0: -21607 }, 'F-': { D: 1.6362096300932947e-9, mu0: 12229 }, 'G3+': { D: 3.0234488685206066e-10, mu0: 3130 }, N: { D: 7.592792804803207e-11, mu0: 22977 }, 'AG4+': { D: 1.4210738829864827e-10, mu0: -42893 } } },
          M: metal,
        },
        regions: [
          { name: 'M1', material: 'M', length: 1e-6 },
          { name: 'water', material: 'water', length: 1.9969045428519623e-6, c0: { 'A+': 14.458449156935107, 'F-': 469.4061672549616, 'G3+': 151.1841389745432, N: 0.13570880266728685, 'AG4+': 0.3488252935992165 } },
          { name: 'M2', material: 'M', length: 1e-6 },
        ],
        interfaces: [
          { reactions: [{ equation: 'S00A(s) = A+ + e-', k0: 32.947969652281095, alpha: 0.5, fixed: { 'S00A(s)': 25802 } }, { equation: 'AG4+ + 3 e- = A+', k0: 3.973379615067484e-15, alpha: 0 }] },
          { reactions: [{ equation: 'AG4+ + 4 e- = N', k0: 0.013451793412656647, alpha: 0.7 }, { equation: 'S11A(s) = A+ + e-', k0: 2.2354030130687728e-6, alpha: 0, fixed: { 'S11A(s)': 22363 } }] },
        ],
        bulkReactions: [{ equation: 'A+ + G3+ = AG4+', kf: { water: 0.29414660725796493 } }],
        contacts: { left: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { V: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
        grid: { hmin: 1.9969045428519622e-8, hmax: 1.3312696952346416e-7 }, // (the fuzzer's own numbers)
      },
      order,
    );
  // A⁺ level with both solids: μ̄ = μ(S00A) + F V_L = μ(S11A) + F V_R.
  for (const order of ['as written', 'species reversed', 'reactions reversed']) {
    const s = silver(order), VL = (22363 - 25802) / FARADAY;
    assert.ok(s.converged && Math.abs(s.terminals.left.V - VL) < 1e-9, `${order}: ${s.converged}, V ${s.terminals.left.V} for ${VL}`);
  }
  const complexed = (order) =>
    solve(
      {
        species: [['C2+', 2], ['F-', -1], ['E+', 1], ['CF+', 1], ['e-', -1]].map(([name, z]) => (name === 'e-' ? { name, z } : { name, z, cRef: 1000 })),
        materials: {
          water: { epsr: 0, species: { 'C2+': { D: 5.183246578681108e-10, mu0: -22707 }, 'F-': { D: 1.7039236094838368e-9, mu0: -2371 }, 'E+': { D: 7.345979739182845e-10, mu0: 16778 }, 'CF+': { D: 4.0102842475617254e-11, mu0: -61263 } } },
          M: { conductor: { species: 'e-', conductivity: 13503964.589985639 } },
        },
        regions: [
          { name: 'M1', material: 'M', length: 1e-6 },
          { name: 'water', material: 'water', length: 4.5574857648705386e-6, c0: { 'C2+': 295.6958039221225, 'F-': 797.4761452984887, 'E+': 206.0481658485044, 'CF+': 0.03637160573936982 } },
          { name: 'M2', material: 'M', length: 1e-6 },
        ],
        interfaces: [
          { reactions: [{ equation: 'C2+ + 2 e- = S00C(s)', k0: 295.1343507546302, alpha: 0.7, fixed: { 'S00C(s)': -14333 } }, { equation: 'CF+ + e- = S01CF(s)', k0: 0.6036981923097351, alpha: 1, fixed: { 'S01CF(s)': -27421 } }] },
          { reactions: [{ equation: 'C2+ + 3 e- = F-', k0: 0.20554335828811326, alpha: 0 }] },
        ],
        bulkReactions: [{ equation: 'C2+ + F- = CF+', kf: { water: 0.006209375707347837 } }],
        contacts: { left: { V: -0.11393122468143702, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' }, right: { I: 0, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
        grid: { hmin: 4.5574857648705384e-8, hmax: 3.0383238432470257e-7 },
      },
      order,
    );
  // F⁻ is CF⁺ less C²⁺, each level with its solid, and the right electrode level with all three:
  // V_R = V_L + (2μ(S00C) − μ(S01CF))/3F.
  for (const order of ['as written', 'species reversed']) {
    const s = complexed(order), VR = -0.11393122468143702 + (2 * -14333 + 27421) / (3 * FARADAY);
    assert.ok(s.converged && Math.abs(s.terminals.right.V - VR) < 1e-9, `${order}: ${s.converged}, V ${s.terminals.right.V} for ${VR}`);
  }
});

test('an electrode at no current that reduces a closed ion: its charge and the ion conserved together, in any species order', () => {
  const RT = GAS_CONSTANT * 298.15, cb = 60.42185, L = 36e-6;
  const def = (order, capacitance) => ({
    species: order.map((name) => ({ name, z: name === 'K+' || name === 'M+' ? 1 : -1, cRef: 1000 })),
    materials: { film: { epsr: 0, species: { 'K+': { D: 2e-9, mu0: -283e3 }, 'X-': { D: 2e-9, mu0: -131e3 }, 'M+': { D: 1e-9, mu0: 70e3 } } } },
    regions: [{ name: 'film', material: 'film', length: L, c0: { 'K+': 60.4, 'X-': 70.4, 'M+': 10 } }],
    contacts: { left: { bath: { c: { 'K+': cb, 'X-': cb }, reference: 'X-' }, V: 0 }, right: { phi: 'neutral' } },
    ports: [{ name: 'el', region: 'film', terminal: 'e-', I: 0, area: 1e5, ...(capacitance ? { capacitance: { C: 0.2, zeroCharge: 0 } } : {}), reactions: [{ equation: 'M+ + e- = M(s)', fixed: { 'M(s)': 0 }, k0: 4e-7 }] }],
    grid: { hmin: 1e-6, hmax: 6e-6 },
  });
  const nernst = (cM) => (70e3 - 131e3 + RT * Math.log((cM * (cM + Math.sqrt(cM * cM + 4 * cb * cb))) / 2 / 1e6)) / FARADAY;
  const amount = (s) => s.conservation.find((c) => c.species === 'M+').amount;
  for (const order of [['K+', 'X-', 'M+', 'e-'], ['M+', 'e-', 'K+', 'X-']]) {
    const s = new Device(def(order, false)).solve();
    assert.ok(s.converged && s.steps === 1, `${order}: steps ${s.steps}`);
    for (const g of [0, s.x.length - 1]) assert.ok(Math.abs(s.c['M+'][g] / 10 - 1) < 1e-9, `${order}: c(M+) ${s.c['M+'][g]}`);
    assert.ok(Math.abs(s.terminals.el.V - nernst(10)) < 1e-9, `${order}: V ${s.terminals.el.V} vs ${nernst(10)}`);
  }
  // With a capacitance, the reaction's charge is left on the electrode as it settles, from the
  // start's V (level with the film's φ): as a long transient finds it.
  const run = new Device(def(['K+', 'X-', 'M+', 'e-'], true));
  const s0 = run.solution(), start = amount(s0) + s0.ports[0].charge / FARADAY;
  for (const order of [['K+', 'X-', 'M+', 'e-'], ['M+', 'e-', 'K+', 'X-']]) {
    const s = new Device(def(order, true)).solve(), cM = s.c['M+'][0];
    assert.ok(s.converged && s.steps === 1, `${order}: steps ${s.steps}`);
    assert.ok(Math.abs((amount(s) + s.ports[0].charge / FARADAY) / start - 1) < 1e-12 && Math.abs(cM / 10 - 1) > 1e-3, `${order}: M+ + Q/F ${amount(s) + s.ports[0].charge / FARADAY} vs ${start}`);
    assert.ok(Math.abs(s.terminals.el.V - nernst(cM)) < 1e-9, `${order}: V ${s.terminals.el.V} vs ${nernst(cM)}`);
  }
  const a = run.advance(1e4);
  assert.ok(Math.abs(a.terminals.el.V - new Device(def(['K+', 'X-', 'M+', 'e-'], true)).solve().terminals.el.V) < 1e-7, `transient V ${a.terminals.el.V}`);
});

// A lit MOS capacitor, its gate jumped to 1.5 V: a nanosecond on, its surface is in deep
// depletion, minority electrons ten orders below where they end, and Newton from there fails
// where it converges from a cold layout for the drives. (A continuation from the last steady
// solve's voltages no longer applies once a transient has moved the state.)
test('a lit MOS capacitor solved a moment into a transient: the steady state its warm solve finds', () => {
  const def = build({
    T: 300,
    library: [semiconductor('Si'), metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }],
    stack: [ohmic(0, ['e-']), layer('Al', units.nm(20)), { phi: { type: 'capacitive', C: 10 }, zeroCharge: 0.05 }, layer('SiO2', units.nm(10)), { dipole: 0 }, layer('Si', units.um(1), { name: 'Si', acceptors: units.perCm3(1e17) }), ohmic(0)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }, { equation: 'photon = e- + h+', fixed: { photon: units.eV(3) }, kf: { Si: 1e-3 } }],
    grid: { hmin: units.nm(0.25), hmax: units.nm(50), ratio: 1.15 },
  });
  const start = () => {
    const d = new Device(def);
    d.solve();
    d.set({ contacts: { left: { V: 1.5 } } });
    return d;
  };
  const warm = start().solve();
  assert.ok(warm.converged);
  for (const t of [1e-9, 1e-3]) {
    const d = start();
    d.advance(t, { tol: 1e-4 });
    const s = d.solve();
    assert.ok(s.converged && s.iterations < 100, `after ${t} s: ${s.iterations} iterations`);
    const n = s.c['e-'], n0 = warm.c['e-'];
    for (let g = 0; g < n.length; g++) if (Number.isFinite(n0[g])) assert.ok(Math.abs(n[g] / n0[g] - 1) < 1e-9, `after ${t} s, node ${g}: ${n[g]} vs ${n0[g]}`);
  }
});

// A floated gate keeps its charge through a regrid: set({ grid }) carries the state over,
// interpolated onto the new nodes. (It had started cold, and the gate fell to flat band, its
// charge lost.)
test('a floating gate keeps its charge through a change of grid', () => {
  const def = build({
    T: 300,
    library: [semiconductor('Si'), metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }],
    stack: [ohmic(1, ['e-']), layer('Al', units.nm(20)), { phi: { type: 'capacitive', C: 10 }, zeroCharge: -0.95 }, layer('SiO2', units.nm(10)), { dipole: 0 }, layer('Si', units.um(1), { name: 'Si', acceptors: units.perCm3(1e17) }), ohmic(0)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
    grid: { hmin: units.nm(0.25), hmax: units.nm(50), ratio: 1.15 },
  });
  const d = new Device(def), Q = d.solve().interfaces[0].D;
  d.set({ contacts: { left: { I: 0 } } });
  d.set({ grid: { hmin: units.nm(0.125), hmax: units.nm(25), ratio: 1.07 } });
  const s = d.solve();
  assert.ok(s.converged && Math.abs(s.interfaces[0].D / Q - 1) < 1e-12, `gate charge ${s.interfaces[0].D} vs ${Q}`);
  assert.ok(Math.abs(s.terminals.left.V - 1) < 2e-3, `V ${s.terminals.left.V}`);
});

// Trap states (D = 0) in two regions of different density: nothing carries them across the
// face, so by default they're blocked there, and each node keeps its own total. (Crossing freely
// by default, the two edge boxes' levels were held together, whatever flux that took: the steady
// solve failed, and a transient mixed the regions' traps, totals off by 1.8 at the face.)
test('immobile trap states in two regions keep their own densities across the face', () => {
  const def = {
    species: [{ name: 'e-', z: -1 }, { name: 'X0', z: 0 }, { name: 'X-', z: -1 }],
    materials: { film: { epsr: 4, species: { 'e-': { D: 1e-6, mu0: 0, cRef: 1 }, X0: { D: 0, mu0: 0, cRef: 1 }, 'X-': { D: 0, mu0: -3000, cRef: 1 } } } },
    regions: [
      { name: 'a', material: 'film', length: 30e-9, fixedCharge: 2 * FARADAY, c0: { 'e-': 1, X0: 1, 'X-': 1 } },
      { name: 'b', material: 'film', length: 20e-9, fixedCharge: 4 * FARADAY, c0: { 'e-': 1, X0: 2, 'X-': 3 } },
    ],
    bulkReactions: [{ equation: 'e- + X0 = X-', kf: { film: 1e3 } }],
    contacts: { left: { V: 0, phi: { type: 'capacitive', C: 0.1 }, zeroCharge: 0 }, right: { V: 0.05, phi: { type: 'capacitive', C: 0.1 }, zeroCharge: 0 } },
    grid: { minCells: 30 },
  };
  for (const s of [new Device(def).solve(), new Device(def).advance(1)]) {
    assert.ok(s.converged);
    s.x.forEach((_, g) => assert.ok(Math.abs(s.c.X0[g] + s.c['X-'][g] - (s.region[g] === 0 ? 2 : 5)) < 1e-12, `node ${g}`));
  }
});

// A diode driven at a current of the wrong sign (into reverse, past what it can pass): the solve
// fails at ±26 V, where depletion reaches both contacts. That state isn't a solution, so its
// accuracy warnings (the space charge cut off at the contacts) are left out: only the reason it
// failed is said, which names the sign slip.
test('a failed steady solve says why it failed, without accuracy warnings about the state it stopped in', () => {
  const def = build({
    T: 300,
    library: [semiconductor('Si')],
    stack: [ohmic(0), layer('Si', 1e-6, { acceptors: units.perCm3(1e16) }), layer('Si', 1e-6, { donors: units.perCm3(1e17) }), ohmic({ I: 100 })],
    grid: { hmin: 1e-9, hmax: 20e-9 },
  });
  const sol = new Device(def).solve();
  assert.ok(!sol.converged);
  assert.equal(sol.warnings.length, 1, sol.warnings.join('\n'));
  assert.match(sol.warnings[0], /^right: no steady state passes the driven current 100 A\/m².*a terminal's current is into the device/);
});

// Gemini 4 Argon's review (its finding 3): the conserved combinations came from eliminating the
// stretches in listing order, so B listed first gave B + C and A − B, with A + C (1e-14 beside
// B's 1000) lost to round-off: no convergence, or converged with A + C 65% off.
test('A + B ⇌ C with B at 1 M and A, C in trace: the same conserved amounts, kept exactly, in either listing', () => {
  for (const trace of [1e-11, 1e-14]) {
    for (const order of [['A', 'B', 'C'], ['B', 'A', 'C'], ['C', 'B', 'A']]) {
      const dev = new Device({
        species: [...order.map((name) => ({ name, z: 0 })), { name: 'Na+', z: 1 }, { name: 'Cl-', z: -1 }],
        materials: {
          water: {
            epsr: 78.3,
            species: { A: { D: 1e-9, cRef: 1000, mu0: 0 }, B: { D: 1e-9, cRef: 1000, mu0: 0 }, C: { D: 1e-9, cRef: 1000, mu0: -15e3 }, 'Na+': { D: 1e-9, cRef: 1000, mu0: 0 }, 'Cl-': { D: 1e-9, cRef: 1000, mu0: 0 } },
          },
        },
        regions: [{ name: 'bath', material: 'water', length: 1e-4, c0: { A: trace, B: 1000, C: trace, 'Na+': 1000, 'Cl-': 1000 } }],
        bulkReactions: [{ equation: 'A + B = C', kf: { water: 1e2 } }],
        contacts: {
          left: { V: 0, phi: 'neutral', terminal: 'Na+', species: { 'Na+': 'equilibrium' } },
          right: { V: 0, phi: 'neutral', terminal: 'Na+', species: { 'Na+': 'equilibrium' } },
        },
      });
      const sol = dev.solve(), ac = sol.c.A[0] + sol.c.C[0];
      assert.ok(sol.converged, `${trace}, [${order}]: ${sol.warnings.join('; ')}`);
      assert.ok(Math.abs(ac / (2 * trace) - 1) < 1e-9, `${trace}, [${order}]: A + C = ${ac}`);
      // At equilibrium C/A = (c_B/c_ref) e^(15 kJ/mol / RT) = 425.
      assert.ok(Math.abs(sol.c.C[0] / sol.c.A[0] / Math.exp(15e3 / (8.314462618 * dev.model.T)) - 1) < 1e-6);
    }
  }
});

test('a face reaction with a fixed participant, between a level a contact holds flat and one it doesn\'t: the steady state is where the transient ends', () => {
  // A(left) + C = B(left), C fixed at −1000 J/mol, the left contact holding μ_A = 0 and μ_B = 1000:
  // the face reaction's affinity, A/RT = fixed − Σ ν η, decides whether A and B are flat to the
  // contact. With + ν it held A at the contact's level and ran the face 3× too fast.
  const def = {
    T: 298.15,
    species: [{ name: 'A', z: 0 }, { name: 'B', z: 0 }, { name: 'S', z: 0 }],
    materials: { w: { epsr: 0, species: { A: { D: 1e-12, mu0: 0, cRef: 1000 }, B: { D: 1e-12, mu0: 0, cRef: 1000 }, S: { D: 1e-9, mu0: 0, cRef: 1000 } } } },
    regions: [{ material: 'w', length: 1e-6, c0: { S: 1 } }, { material: 'w', length: 1e-6, c0: { S: 1, A: 1, B: 1 } }],
    interfaces: [{ species: { A: 'blocked', B: 'blocked', S: 'equilibrium' }, reactions: [{ equation: 'A(left) + C = B(left)', fixed: { C: -1000 }, k0: 1e-3 }] }],
    contacts: { left: { V: 0, phi: 'neutral', species: { A: { type: 'equilibrium', mu: 0 }, B: { type: 'equilibrium', mu: 1000 }, S: 'blocked' } }, right: { V: 0, species: { A: 'blocked', B: 'blocked', S: 'blocked' } } },
  };
  const d = new Device(def), steady = d.solve(), late = new Device(def).advance(1e3, { tol: 1e-8 });
  const r = (s) => s.interfaces[0].rates[0];
  assert.ok(steady.converged && late.converged && Math.abs(r(steady) / r(late) - 1) < 1e-6, `face rate ${r(steady)} vs the transient's ${r(late)}`);
  assert.ok(check(d, steady).ok, check(d, steady).text);
});
