import { test } from 'node:test';
import assert from 'node:assert/strict';
import { DeviceError, units } from '../src/index.js';
import { build, layer, ohmic, describe, unitWarnings } from '../src/kit.js';

// describe(): a readable summary, with warnings for likely unit slips and unresolved double layers.

const silicon = (eSpecies = {}) => ({
  species: [
    { name: 'e-', z: -1 },
    { name: 'h+', z: 1 },
  ],
  materials: {
    Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: units.perCm3(2.8e19), ...eSpecies }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: units.perCm3(1.04e19) } } },
  },
});
const pn = (lib, grid = { hmin: 1e-9, hmax: 20e-9 }) =>
  build({
    library: [lib],
    stack: [ohmic(0), layer('Si', 1e-6, { name: 'n', donors: units.perCm3(1e17) }), layer('Si', 1e-6, { name: 'p', acceptors: units.perCm3(1e16) }), ohmic(0.5)],
    bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e13 } }],
    grid,
  });

test('a clean pn diode: regions with doping and Debye lengths, faces, contacts, time scales, no warnings', () => {
  const text = describe(pn(silicon()));
  assert.match(text, /^2 regions, 2 µm, \d+ nodes, T = 298.15 K/);
  assert.match(text, /n \[Si\] 1 µm: εr 11.7, fixed charge donor-like 1.00e\+17 cm⁻³, Debye length 12.9 nm/);
  assert.match(text, /p \[Si\] 1 µm: .*acceptor-like 1.00e\+16 cm⁻³, Debye length 40.8 nm/);
  assert.match(text, /n \| p: same material/);
  assert.match(text, /right: held at 0.5 V \(terminal e-\); exchanges e-, h\+; φ bulk/);
  assert.match(text, /e- \+ h\+ = 0 in Si/);
  assert.match(text, /diffusion across n: L²\/D up to 833 ps/);
  assert.doesNotMatch(text, /warnings/);
});

test('unit slips are flagged where they are, and a coarse grid against the Debye length', () => {
  const slips = unitWarnings(pn(silicon({ D: 36, mu0: 4.05 })));
  assert.equal(slips.length, 2);
  assert.match(slips[0], /^materials\.Si\.species\.e-\.D: .*cm2PerS/);
  assert.match(slips[1], /^materials\.Si\.species\.e-\.mu0: 4.05 J\/mol is tiny/);
  const more = unitWarnings({
    T: 25,
    species: [{ name: 'Na+', z: 1, cRef: 1e8 }],
    regions: [{ material: 'x', length: 5, fixedCharge: 1e17, c0: { 'Na+': 2e6 } }],
    interfaces: [{ phi: { type: 'capacitive', C: 20 }, reactions: [{ k0: 100 }] }],
    contacts: { left: { V: 300 } },
  });
  for (const path of ['T', 'species[0].cRef', 'regions[0].length', 'regions[0].fixedCharge', 'regions[0].c0.Na+', 'interfaces[0].phi.C', 'interfaces[0].reactions[0].k0', 'contacts.left.V']) {
    assert.ok(more.some((w) => w.startsWith(`${path}:`)), path);
  }
  const coarse = describe(pn(silicon(), { hmin: 30e-9, hmax: 50e-9 }));
  assert.match(coarse, /warnings:\n {2}n: end cells of [\d.]+ nm are coarser than the Debye length, 12.9 nm/);
  assert.throws(() => describe({ species: [] }), DeviceError);
});
