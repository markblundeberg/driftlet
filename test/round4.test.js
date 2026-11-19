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
