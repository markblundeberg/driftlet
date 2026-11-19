import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device } from '../src/index.js';
import { aqueous, build, layer } from '../src/kit.js';

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
