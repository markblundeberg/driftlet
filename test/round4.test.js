import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, EPS0, FARADAY, SolverError, units } from '../src/index.js';
import { complexSolve } from '../src/impedance.js';
import { aqueous, bath, build, half, layer, metal, ohmic, semiconductor } from '../src/kit.js';

// Faults found in the impedance and the readout by an outside reading pass (OpenAI's Sol, its
// fourth), each confirmed by a probe against a second route before it was fixed.

// A film with a contact passing K⁺ on the left, a membrane port `mem` passing K⁺ in its middle and
// a Cl⁻ reference port, both current-driven terminals sharing the K⁺ that the law conserves.
const sharedK = (left, mem) => {
  const film = { 'Na+': 100, 'K+': 10, 'Cl-': 110 }, K = (G) => ({ 'K+': { type: 'conductance', G, offset: 0 } });
  return build({
    T: 298.15,
    library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })],
    stack: [{ ...left, phi: 'neutral', species: K(1e6) }, layer('water', 1e-6, { name: 'film', c0: film }), { V: 0, phi: { type: 'capacitive', C: 0.5 }, zeroCharge: -0.05 }],
    ports: [{ name: 'mem', region: 'film', from: 0.4e-6, to: 0.6e-6, ...mem, species: K(1e13) }, { name: 'ref', region: 'film', from: 0.7e-6, to: 0.8e-6, V: 0, species: { 'Cl-': { type: 'equilibrium', offset: 0 } } }],
    grid: { hmin: 5e-8, hmax: 5e-8 },
  });
};
// A neutral film between a pinned contact and a capacitive plate, the plate driven as given.
const plate = (right) => new Device({
  species: [{ name: 'X', z: 0 }],
  materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
  regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
  contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } }, right: { ...right, phi: { type: 'capacitive', C: 1e-3 }, zeroCharge: 0, species: { X: 'blocked' } } },
  grid: { minCells: 8 },
});
const near = (a, b, q, tol) => Math.hypot(a.Z.re[q] - b.Z.re[q], a.Z.im[q] - b.Z.im[q]) / Math.hypot(b.Z.re[q], b.Z.im[q]) < tol;

test('a current-driven terminal measured where another shares its conserved law: the held twin\'s Z', () => {
  // (The left contact's row read by the law's change, iω F δS/λ = δI_left, held its δI to 0
  // where mem's own small signal fed the same K⁺: the left AC-shorted, mem's Im Z 37× too small
  // at 1 kHz and 3.6e7× at 1 Hz, relError 1e-13. A sine stepped through step() agrees with the
  // held twin.)
  const fs = [1e3, 1];
  const d = new Device(sharedK({ I: 0 }, { I: 0 })), V = d.solve().terminals.mem.V;
  const z = d.impedance(fs, { terminal: 'mem' });
  const twin = new Device(sharedK({ I: 0 }, { V })).impedance(fs, { terminal: 'mem' });
  fs.forEach((f, q) => assert.ok(near(z, twin, q, 1e-5), `${f} Hz: ${z.Z.re[q]} ${z.Z.im[q]}i, held twin ${twin.Z.re[q]} ${twin.Z.im[q]}i`));
  // (the left contact, its law its own, as before: its held twin's)
  const zl = d.impedance(fs, { terminal: 'left' });
  const twinL = new Device(sharedK({ V: d.solution().terminals.left.V }, { I: 0 })).impedance(fs, { terminal: 'left' });
  fs.forEach((f, q) => assert.ok(near(zl, twinL, q, 1e-5), `left, ${f} Hz: ${zl.Z.re[q]} ${zl.Z.im[q]}i, held twin ${twinL.Z.re[q]} ${twinL.Z.im[q]}i`));
});

test('a solution read after a set() of a floating terminal\'s drive alone: what it passed as solved, until it moves', () => {
  // (Read at the new drive, a capacitive plate set from I: 0 to 2 A/m² read 2 A/m² flowing at no
  // step, the plate uncharged and the terminals summing to −2, without a word; behind 10 Ω·m², a
  // source set from 0 to 0.5 V read 0.05 A/m². Held terminals were read as solved already.)
  const sum = (s) => Object.values(s.terminals).reduce((t, x) => t + x.current, 0);
  for (const [from, to] of [[{ I: 0 }, { I: 2 }], [{ V: 0, R: 10 }, { V: 0.5 }]]) {
    const d = plate(from);
    d.solve();
    d.set({ contacts: { right: to } });
    const s = d.solution();
    assert.ok(s.terminals.right.current === 0, `${JSON.stringify(to)}: ${s.terminals.right.current} A/m²`);
    assert.ok(sum(s) === 0, `${sum(s)}`);
    const t = d.step(1e-9); // (and the step from there passes the new drive)
    assert.ok(Math.abs(t.plates.right.charge - (to.I ?? 0.05) * 1e-9) < 1e-3 * (to.I ?? 0.05) * 1e-9, `${t.plates.right.charge}`);
  }
  // (a port's current likewise)
  const d = new Device(sharedK({ I: 0 }, { I: 0 }));
  d.solve();
  d.set({ ports: { mem: { I: 1e-3 } } });
  assert.ok(d.solution().terminals.mem.current === 0);
  assert.equal(d.step(1e-9).terminals.mem.current, 1e-3);
});

test('a refused impedance() leaves a pending step response as it was', () => {
  // (Each solved the steady state at the new drive before refusing: the next step read 0 where
  // the plate's charging current was 3420 A/m².)
  for (const bad of [[[0], {}], [[1], { terminal: 'nope' }], [[NaN], {}], [[1, -1], {}]]) {
    const d = plate({ V: 0 }), twin = plate({ V: 0 });
    for (const x of [d, twin]) {
      x.solve();
      x.set({ contacts: { right: { V: 0.1 } } });
    }
    assert.throws(() => d.impedance(...bad), /impedance: (frequencies|no terminal)/);
    const a = d.step(1e-9).terminals.right.current, b = twin.step(1e-9).terminals.right.current;
    assert.ok(b > 3e3 && a === b, `${JSON.stringify(bad)}: ${a} A/m², its twin ${b}`);
  }
});

test('a terminal no small-signal current returns through is refused, not read as NaN', () => {
  // (A dielectric up to a neutral contact passes nothing: Z came out NaN, with a warning that it
  // was round-off.)
  const d = new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } }, right: {} },
    grid: { minCells: 8 },
  });
  for (const terminal of ['left', 'right']) assert.throws(() => d.impedance([1, 1e6], { terminal }), new RegExp(`no current returns through '${terminal}' at 1 Hz`));
});

test('a plate\'s impedance at an absurd frequency, and a singular small system refused', () => {
  // (Divided by c² + d², which overflowed: at 1e160 Hz Z read 0, and behind 0.4 Ω·m², NaN; a zero
  // pivot gave NaN for every floating terminal's δV.)
  const both = (R) => new Device({
    species: [{ name: 'X', z: 0 }],
    materials: { f: { epsr: 4, species: { X: { D: 1e-9, cRef: 1, mu0: 0 } } } },
    regions: [{ name: 'a', material: 'f', length: 1e-6, c0: { X: 1 } }],
    contacts: { left: { V: 0, ...(R ? { R } : {}), phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } }, right: { V: 0, phi: 'pinned', zeroCharge: 0, species: { X: 'blocked' } } },
    grid: { minCells: 8 },
  });
  const C = (4 * EPS0) / 1e-6;
  for (const f of [1e160, 1e300]) {
    for (const R of [0, 0.4]) {
      const z = both(R).impedance([f]), im = -1 / (2 * Math.PI * f * C);
      assert.ok(Math.abs(z.Z.re[0] - R) < 1e-12 && Math.abs(z.Z.im[0] / im - 1) < 1e-12, `${f} Hz, R ${R}: ${z.Z.re[0]} ${z.Z.im[0]}i, expected ${R} ${im}i`);
    }
  }
  const zero = () => [[Float64Array.from([0, 0, 1]), Float64Array.from([0, 0, 1])], [new Float64Array(3), new Float64Array(3)]];
  assert.throws(() => complexSolve(...zero(), 2, 1), (e) => e instanceof SolverError && /singular at 1 Hz/.test(e.message));
  const rank1 = [[Float64Array.from([1, 2, 1]), Float64Array.from([2, 4, 3])], [new Float64Array(3), new Float64Array(3)]];
  assert.throws(() => complexSolve(...rank1, 2, 1), SolverError);
});

test('a floating terminal\'s δV carries its readings\' error: a gate behind R at 1 nHz, a Pt island behind a bath\'s R', () => {
  // (Read from its laws and cuts with no error of its own, δV went into Kirchhoff's readout as
  // exact: a gate behind 1e4 Ω·m², measured at the back contact at 1 nHz, put Re Z − R at −5.4
  // for the held gate's 2869, relError 1e-9; a Pt island behind 1 Ω·m² read Z − R 0.3% off at
  // 3e-5–1e-4 Hz, relError 5e-12, its bath's cuts 2000× past their estimates. Now the floating
  // responses are split into their regions' shifts and the rest where a reading is unresolved.)
  const si = semiconductor('Si');
  si.materials.Si.species['e-'].mu0 -= FARADAY * 0.3;
  si.materials.Si.species['h+'].mu0 += FARADAY * 0.3;
  const mos = (gate) => build({
    T: 300, library: [si, metal('Au')], materials: { SiO2: { epsr: 3.9, species: {} } },
    stack: [ohmic(0), layer('Si', units.um(15), { acceptors: units.perCm3(1e16) }), { dipole: -0.3 }, layer('SiO2', units.nm(20), { grid: { hmin: units.nm(0.5), hmax: units.nm(1) } }), { phi: { type: 'capacitive', C: 100 }, zeroCharge: 0.7 }, layer('Au', units.nm(50)), ohmic(gate, ['e-'])],
    bulkReactions: [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: 1e-6, tauP: 1e-6 } } }],
    grid: { hmin: units.nm(0.1), hmax: units.nm(10), ratio: 1.1 },
  });
  const fs = [1e-9, 1e-8], ref = new Device(mos(0.6)).impedance(fs, { terminal: 'right' });
  for (const R of [1e3, 1e4]) {
    const z = new Device(mos({ V: 0.6, R })).impedance(fs, { terminal: 'left' });
    fs.forEach((f, q) => {
      const e = Math.abs(z.Z.re[q] - R - ref.Z.re[q]), said = z.relError[q] * Math.hypot(z.Z.re[q], z.Z.im[q]);
      assert.ok(e < 1e-3 * ref.Z.re[q] && e < said * 10 + 1e-3, `R ${R}, ${f} Hz: Re Z − R ${z.Z.re[q] - R} for ${ref.Z.re[q]}, its error said ${said}`);
    });
  }
  const c = { 'K+': 100, 'Cl-': 125, 'Fe3+': 5, 'Fe2+': 5 };
  const face = { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1e-3, alpha: 0.5 }], phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0 };
  const island = (left) => build({ library: [aqueous(['K+', 'Cl-', 'Fe3+', 'Fe2+']), metal('Pt')], stack: [bath(c, 'Cl-', left), layer('water', 5e-6, { c0: c }), face, layer('Pt', 1e-6), face, layer('water', 5e-6, { c0: c }), bath(c, 'Cl-', 0)], grid: { hmin: 0.3e-9, hmax: 1e-6, ratio: 1.3 } });
  const gs = [3e-5, 1e-4], Z0 = new Device(island({ V: 0 })).impedance(gs);
  for (const R of [1, 1e3]) {
    const z = new Device(island({ V: 0, R })).impedance(gs);
    gs.forEach((f, q) => assert.ok(Math.hypot(z.Z.re[q] - R - Z0.Z.re[q], z.Z.im[q] - Z0.Z.im[q]) < 1e-5 * Math.hypot(Z0.Z.re[q], Z0.Z.im[q]), `R ${R}, ${f} Hz: Z − R ${z.Z.re[q] - R} ${z.Z.im[q]}i for ${Z0.Z.re[q]} ${Z0.Z.im[q]}i`));
    assert.deepEqual(z.warnings.filter((w) => /impedance:/.test(w)), []);
  }
});

// A Pt electrode in a Fe³⁺/Fe²⁺ bath at its open circuit: at 10 nHz its response is nearly a
// uniform shift, which the readout splits off (Sol's fourth pass, S5).
const redox = (right) => {
  const c = { 'Fe3+': 9.6, 'Fe2+': 0.536, 'Cl-': 29.872 };
  return build({
    T: 300, library: [aqueous(['Fe3+', 'Fe2+', 'Cl-'], { epsr: 78.3 }), metal('Pt')],
    stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.395 }, zeroCharge: -0.119, reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 0.0411, alpha: 0.377 }] }, layer('water', 9.83e-6, { c0: { ...c } }), bath({ ...c }, 'Cl-', right)],
    grid: { hmin: 2.6e-10, hmax: 3.7e-7, ratio: 1.136 },
  });
};

test('the profiles are those of the reading Z kept', () => {
  // (The split's solve overwrote them even where Z set it aside. Here its shift is made
  // nonsense, so that it is: the profiles moved with the nonsense, 0.1% at a node.)
  const run = (scale) => {
    const d = new Device(redox({ I: 0 }));
    d.solver._regionShift = (yr, yi) => [yr.map((v) => v * scale), yi.map((v) => v * scale)];
    return d.impedance([1e-8], { profiles: true });
  };
  const a = run(1e6), b = run(1e9);
  assert.ok(a.relError[0] > 1e-2, `${a.relError[0]}`); // (the split set aside: Z is the plain reading's, warned)
  assert.deepEqual(Array.from(a.profiles[0].phi.re), Array.from(b.profiles[0].phi.re));
  assert.deepEqual(Array.from(a.profiles[0].mu['Fe3+'].im), Array.from(b.profiles[0].mu['Fe3+'].im));
});

test('a GMRES left above its tolerance, where relError doesn\'t show it, is warned of', () => {
  // (No device in the tests leaves it so; the solver's report of where it did is passed on.)
  const d = new Device(redox({ I: 0 })), was = d.solver.impedance.bind(d.solver);
  d.solver.impedance = (fs, opts) => ({ ...was(fs, opts), unconverged: [1e-3] });
  const z = d.impedance([1e-3, 1]);
  assert.ok(!('unconverged' in z));
  assert.equal(z.warnings.filter((w) => /stopped short of its tolerance at 0\.001 Hz \(1 of 2/.test(w)).length, 1);
  assert.equal(new Device(redox({ I: 0 })).impedance([1e-3, 1]).warnings.filter((w) => /stopped short/.test(w)).length, 0);
});
