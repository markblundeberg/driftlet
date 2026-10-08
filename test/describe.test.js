import { test } from 'node:test';
import assert from 'node:assert/strict';
import { DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, bath, aqueous, metal, IONS, describe, unitWarnings } from '../src/kit.js';

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

test("the Debye length counts the doping's carriers, even when c0 lists a trace of something else", () => {
  const def = pn(silicon());
  def.species.push({ name: 'T', z: 1, cRef: 1 });
  def.materials.Si.species.T = { D: 0, mu0: 0 };
  def.regions[0].c0 = { T: 1e-12 };
  assert.match(describe(def), /n \[Si\] 1 µm: .*Debye length 12.9 nm/);
});

test('unit slips are flagged where they are, and a coarse grid against the Debye length', () => {
  const slips = unitWarnings(pn(silicon({ D: 36, mu0: 4.05 })));
  assert.equal(slips.length, 2);
  assert.match(slips[0], /^materials\.Si\.species\.e-\.D: .*cm2PerS/);
  assert.match(slips[1], /^materials\.Si\.species\.e-\.mu0: 4.05 J\/mol is small/);
  const more = unitWarnings({
    T: 25,
    species: [{ name: 'Na+', z: 1, cRef: 1e8 }],
    regions: [{ material: 'x', length: 5, fixedCharge: 1e17, c0: { 'Na+': 2e6 } }],
    interfaces: [{ phi: { type: 'capacitive', C: 20 }, reactions: [{ k0: 100 }] }],
    contacts: { left: { V: 300 } },
    ports: [{ capacitance: { C: 20 }, reactions: [{ k0: 100 }] }],
  });
  for (const path of ['T', 'species[0].cRef', 'regions[0].length', 'regions[0].fixedCharge', 'regions[0].c0.Na+', 'interfaces[0].phi.C', 'interfaces[0].reactions[0].k0', 'contacts.left.V', 'ports[0].capacitance.C', 'ports[0].reactions[0].k0']) {
    assert.ok(more.some((w) => w.startsWith(`${path}:`)), path);
  }
  // Ions: D in cm²/s and μ° in kJ/mol, the commonest slips in electrochemistry.
  const ions = {
    species: [
      { name: 'Na+', z: 1, cRef: 1000 },
      { name: 'Cl-', z: -1, cRef: 1000 },
    ],
    materials: { water: { epsr: 78.4, species: { 'Na+': { D: 1.334e-5, mu0: -261.905 }, 'Cl-': { D: 2.032e-9, mu0: -131.228e3 } } } },
    contacts: { left: { bath: { c: { 'Na+': 1e20, 'Cl-': 1e20 }, reference: 'Cl-' } } },
  };
  const ionSlips = unitWarnings(ions);
  for (const path of ['materials.water.species.Na+.D', 'materials.water.species.Na+.mu0', 'contacts.left.bath.c.Na+']) {
    assert.ok(ionSlips.some((w) => w.startsWith(`${path}:`)), path);
  }
  assert.ok(!ionSlips.some((w) => w.startsWith('materials.water.species.Cl-')));
  assert.deepEqual(unitWarnings({ T: 77 }), [], 'liquid nitrogen is a real temperature');
  const coarse = describe(pn(silicon(), { hmin: 30e-9, hmax: 50e-9 }));
  assert.match(coarse, /warnings:\n {2}n: end cells of [\d.]+ nm are coarser than the Debye length, 12.9 nm/);
  assert.throws(() => describe({ species: [] }), DeviceError);
});

test('a bath read through a reference ion says where the SHE level sits against its terminal', () => {
  // On the usual tables' μ°, SHE is φ in a solution, and a Cl⁻ terminal reads φ − (μ°_Cl⁻ + RT ln a)/(−F).
  const cell = (c, offset) =>
    build({
      library: [aqueous(Object.keys(c), { epsr: 0 }), metal('Pt')],
      stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), {}, layer('water', 1e-4), bath(c, 'Cl-', 0, offset === undefined ? {} : { offset })],
    });
  const RT = GAS_CONSTANT * 298.15, gap = (c) => -(IONS['Cl-'].mu0 + RT * Math.log(c / 1000)) / FARADAY;
  const kcl = { 'K+': 500, 'Cl-': 500 };
  assert.match(describe(cell(kcl)), new RegExp(`right: .*the bath's φ is at V − ${gap(500).toFixed(3)} V \\(on the usual tables' μ°`));
  assert.match(describe(cell(kcl, gap(500))), /the bath's φ is at V \(/, 'an offset of that much reads SHE');
  // With H⁺ in the bath, its own μ° says where SHE is.
  assert.match(describe(cell({ 'H+': 100, 'Cl-': 100 })), new RegExp(`the bath's SHE level is at V − ${(gap(100) - IONS['H+'].mu0 / FARADAY).toFixed(3)} V`));
  const salt = cell(kcl);
  assert.doesNotMatch(describe({ ...salt, contacts: { ...salt.contacts, right: { V: 0, bath: { c: kcl } } } }), /SHE/, 'a bath read by its φ says nothing');
});

test('an outside out of equilibrium with a reaction at the contact is warned of; ohmic contacts with recombination are not', () => {
  // H⁺ + OH⁻ ⇌ H₂O with the tables' μ°: a bath at 0.1 mM each is a little short of K_w.
  const water = (c) =>
    build({
      library: [aqueous(['Na+', 'Cl-', 'H+', 'OH-'], { epsr: 0 })],
      stack: [bath({ 'Na+': 100, 'Cl-': 100, 'H+': c, 'OH-': c }, 'Cl-'), layer('water', 1e-5), bath({ 'Na+': 100, 'Cl-': 100, 'H+': c, 'OH-': c }, 'Cl-')],
      bulkReactions: [{ equation: 'H+ + OH- = H2O', fixed: { H2O: -237.1e3 }, kf: { water: 1.4e8 } }],
    });
  const off = describe(water(1e-4));
  const a = +off.match(/contacts\.left: the outside isn't in equilibrium with bulkReactions\[0\] \(H\+ \+ OH- = H2O; A = (-?[\d.]+) RT/)[1];
  assert.match(off, /contacts\.right: the outside isn't in equilibrium/);
  // at c e^{−a/2}, each μ up by −a/2 RT, it's in equilibrium
  assert.doesNotMatch(describe(water(1e-4 * Math.exp(-a / 2))), /isn't in equilibrium/);
  assert.doesNotMatch(describe(pn(silicon())), /isn't in equilibrium/, 'ohmic contacts hold e⁻ and h⁺ at one level');
});
