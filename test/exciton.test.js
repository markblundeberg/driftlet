import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { photogeneration } from '../src/kit.js';

// An organic bilayer solar cell, where light makes neutral excitons rather than free carriers: a
// third species, beside electrons and holes, that diffuses regardless of the field, decays, and
// splits only at the donor/acceptor interface. Illustrative numbers of the order of CuPc/C60.
//
//   anode (W 5.0 eV) | donor: holes, excitons | acceptor: electrons | cathode (W 4.2 eV)
//
// The donor's HOMO is 5.2 eV below vacuum, the acceptor's LUMO 4.0 eV (so the interface's
// charge-transfer gap is E_CT = 1.2 eV), the exciton's energy 1.7 eV; the levels are aligned on
// the vacuum level (dipole 0), the contacts pinned at their work functions. At the interface,
// X → h⁺(donor) + e⁻(acceptor) splits excitons, and e⁻ + h⁺ → 0 is charge-transfer recombination.

const T = 300, VT = (GAS_CONSTANT * T) / FARADAY, N = units.perCm3(1e21), eV = units.eV;
const EX = 1.7, HOMO = 5.2, LUMO = 4.0, tauX = 1e-9, alpha = 1e7, flux = 1e-3, kSplit = 1e3, kCT = 1e-2;

const cell = ({ d, LX, right = { V: 0 }, light = 1, h = 0.5e-9 }) => ({
  T,
  species: [
    { name: 'e-', z: -1, cRef: N },
    { name: 'h+', z: 1, cRef: N },
    { name: 'X', z: 0, cRef: N },
  ],
  materials: {
    donor: { epsr: 3.5, species: { 'h+': { D: 1e-7 * VT, mu0: eV(HOMO) }, X: { D: (LX * LX) / tauX, mu0: eV(EX) } } },
    acceptor: { epsr: 3.5, species: { 'e-': { D: 1e-6 * VT, mu0: eV(-LUMO) } } },
  },
  regions: [
    { material: 'donor', length: d, c0: { X: 1e-12 } },
    { material: 'acceptor', length: 40e-9 },
  ],
  interfaces: [{ dipole: 0, reactions: [{ equation: 'X = h+ + e-', k0: kSplit }, { equation: 'e- + h+ = 0', k0: kCT }] }],
  contacts: {
    left: { V: 0, terminal: 'h+', species: { 'h+': 'equilibrium' }, phi: 'pinned', zeroCharge: 5.0 },
    right: { ...right, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'pinned', zeroCharge: 4.2 },
  },
  bulkReactions: [
    { equation: 'X = 0', kf: { donor: 1 / tauX } },
    ...(light > 0 ? [photogeneration({ material: 'donor', flux: flux * light, alpha, mu: eV(2.0), to: d, makes: 'X' })] : []),
  ],
  grid: { hmin: 0.05e-9, hmax: h },
});

test('excitons diffuse to the interface and split: the current is qΦ∫αe^{−αx}η(x)dx, with η from exciton diffusion', () => {
  // Light enters through the anode, which blocks excitons (X′ = 0); the interface takes them at
  // the velocity v its splitting reaction gives, k₀ e^{α(E_X − E_CT)/kT}/c_ref, so −D X′ = v X
  // there. Then η(x) = cosh(x/L)/(cosh(d/L) + (D/vL) sinh(d/L)) with L = √(Dτ), and the integral
  // is closed-form. Every split exciton is one electron of current, less what recombines at the
  // interface: J = F(r_split − r_CT), exactly.
  const v = (kSplit * Math.exp((0.5 * (EX - (HOMO - LUMO))) / VT)) / N;
  for (const [d, LX] of [[20e-9, 10e-9], [60e-9, 10e-9], [20e-9, 30e-9], [100e-9, 30e-9]]) {
    const s = new Device(cell({ d, LX })).solve();
    assert.ok(s.converged);
    const D = (LX * LX) / tauX, k = 1 / LX;
    const integral = 0.5 * (-Math.expm1(-(alpha - k) * d) / (alpha - k) + -Math.expm1(-(alpha + k) * d) / (alpha + k));
    const split = (flux * alpha * integral) / (Math.cosh(d / LX) + (D / (v * LX)) * Math.sinh(d / LX));
    const [rSplit, rCT] = s.interfaces[0].rates;
    assert.deepEqual([s.interfaces[0].left, s.interfaces[0].right], ['region 0', 'region 1']); // the regions it joins
    assert.ok(Math.abs(rSplit / split - 1) < 5e-4, `d ${d}, L ${LX}: split ${rSplit / split} of theory`);
    assert.ok(Math.abs(-s.current / (FARADAY * (rSplit - rCT)) - 1) < 1e-12); // holes leave by the anode: J < 0
    assert.ok(rCT < 1e-2 * rSplit, `little recombines at short circuit: ${rCT / rSplit}`);
  }
});

test('the open-circuit voltage is set by charge-transfer recombination: V_oc = V_T ln(J_ph/J₀ + 1), J₀ = F k₀ e^{−(1−α)E_CT/kT}', () => {
  // At open circuit no current flows in either layer, so the quasi-Fermi levels are flat from
  // each contact to the interface, and there μ̄_e + μ̄_h = qV_oc. The interface's recombination
  // e⁻ + h⁺ → 0 is then k₀(c_e c_h/c_ref²)^{1−α}(e^{αa} − e^{−(1−α)a}) with a = (qV_oc − E_CT)/kT
  // + ln(c_e c_h/c_ref²)... which collapses to J₀(e^{V_oc/V_T} − 1): detailed balance across the
  // charge-transfer gap, whatever the layers' transport. (Recombination back through the exciton,
  // 0.5 eV higher, is e^{−20} smaller.)
  const J0 = FARADAY * kCT * Math.exp((-0.5 * (HOMO - LUMO)) / VT);
  for (const light of [0.01, 1, 3]) {
    const s = new Device(cell({ d: 20e-9, LX: 10e-9, right: { I: 0 }, light })).solve();
    assert.ok(s.converged);
    const Voc = s.terminals.left.V - s.terminals.right.V; // the anode is the positive one
    const theory = VT * Math.log((FARADAY * s.interfaces[0].rates[0]) / J0 + 1);
    assert.ok(Math.abs(Voc - theory) < 1e-6, `light ${light}×: V_oc ${Voc}, theory ${theory}`);
  }
});
