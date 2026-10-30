import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT, units } from '../src/index.js';
import { build, layer, ohmic, bath, half, aqueous, metal, semiconductor, polarization, unitWarnings } from '../src/kit.js';
import { mhcIntegral } from '../src/marcus.js';

// Marcus–Hush–Chidsey electron transfer (src/marcus.js): the integral against an independent
// high-precision quadrature, its exact identities and limits, and an electrode's steady current
// against the textbook closed form for kinetics in series with diffusion.

const RT = GAS_CONSTANT * 298.15;
const VT = RT / FARADAY;

// [λ/RT, a, I(λ, a), ∂I/∂a] from mpmath at 30 digits (mp.quad over (−∞, ∞) with breakpoints at
// the Gaussian's centre, ±5√λ and ±20√λ about it, 0, ±5 and ±40; quoted to 17 digits; where a < −30,
// as e^{a} I(λ, −a)). Two entries at a = 200 have ∂I/∂a below mpmath's noise, ~1e−38.
const REFERENCE = [
  [0.3, -200, 2.6870093566078806e-87, 2.6870093566078806e-87], // (far enough back that only the reflection gets it)
  [0.3, -100, 7.2229959117179891e-44, 7.2229959117179891e-44],
  [0.3, -30, 1.8169003236223378e-13, 1.816900323622028e-13],
  [0.3, -5, 0.012925409681114518, 0.01277157511566042],
  [0.3, -0.5, 0.64065848317406594, 0.38272870761817921],
  [0.3, 0, 0.84300853940048788, 0.42150426970024394],
  [0.3, 0.5, 1.0562672684635627, 0.42525430730590042],
  [0.3, 5, 1.9183008834692328, 0.022831073838214698],
  [0.3, 20, 1.9416259052635994, 7.2920998118664276e-9],
  [0.3, 60, 1.9416259125556993, 3.0979423689766243e-26],
  [0.3, 200, 1.9416259125556993, 1.248959218037643e-38],
  [4, -30, 6.6343819447264319e-13, 6.6343819428757953e-13],
  [4, -5, 0.029531361099831875, 0.024056619515506422],
  [4, -0.5, 0.6323398445254518, 0.33786717171670416],
  [4, 0, 0.8189382741758503, 0.40946913708792515],
  [4, 0.5, 1.0425521519803244, 0.48550335929970149],
  [4, 5, 4.3828425934249896, 0.81252369379998572],
  [4, 20, 7.0897722412174859, 4.2859225042250704e-5],
  [4, 60, 7.0898154036220641, 1.8506396755882244e-22],
  [4, 200, 7.0898154036220641, -3.4897525286559567e-38],
  [20, -100, 5.8975514659726814e-43, 5.8975514659726814e-43],
  [20, -30, 1.3882051240132755e-12, 1.3599086449314131e-12],
  [20, -5, 0.0012123721163845125, 0.00073095229733032094],
  [20, -0.5, 0.014847323396541013, 0.0075770970188823888],
  [20, 0, 0.019113659027958031, 0.0095568295139790153],
  [20, 0.5, 0.024479097896840842, 0.011986576871650917],
  [20, 5, 0.17993197580050218, 0.071449036200423379],
  [20, 20, 7.926654595212022, 0.9620331907272771],
  [20, 60, 15.853309174085991, 1.4355154636766848e-8],
  [20, 200, 15.853309190424044, 2.204048191412357e-38],
  [40, -30, 2.8651052426826969e-13, 2.362053572042752e-13],
  [40, -5, 9.6243954454005861e-6, 5.3527015340713729e-6],
  [40, -0.5, 0.00010493627571341358, 5.3057937900829006e-5],
  [40, 0, 0.00013493030856446541, 6.7465154282232705e-5],
  [40, 0.5, 0.00017301066983675824, 8.5532919040174948e-5],
  [40, 5, 0.0014283869325043508, 0.00063397558809961027],
  [40, 20, 0.31880938104461426, 0.088776397409528904],
  [40, 60, 22.101155484547101, 0.088776397409528904],
  [40, 200, 22.419964865591715, -1.7714355396190859e-37],
  [150, -30, 1.1448598552437833e-23, 6.8327530578176712e-24], // (as e^{−30} I(150, 30), and by a second substitution)
  [150, -5, 1.2615713544795335e-17, 6.5116115738021901e-18],
  [150, -0.5, 1.2457613303692481e-16, 6.2489272715628498e-17],
  [150, 0, 1.6002352264931065e-16, 8.0011761324655326e-17],
  [150, 0.5, 2.0539132035954689e-16, 1.0236392724270902e-16],
  [150, 5, 1.8723379015162407e-15, 9.0592905699935702e-16],
  [150, 20, 1.8475981641538238e-12, 8.0448199441141574e-13],
  [150, 60, 5.159972355883463e-6, 1.5838875705608551e-6],
  [150, 200, 43.32724101584692, 0.016134461812171273],
];

test('the Marcus–Hush–Chidsey integral and its slope against mpmath at 30 digits', () => {
  for (const [lam, a, I, dI] of REFERENCE) {
    const r = mhcIntegral(lam, a);
    assert.ok(Math.abs(r.I / I - 1) < 1e-13, `λ ${lam}, a ${a}: ${r.I} vs ${I}`);
    // (relative to itself, unless it's below the integral's round-off)
    assert.ok(Math.abs(r.dI - dI) < 1e-13 * Math.max(Math.abs(dI), 1e-2 * I), `λ ${lam}, a ${a}: ∂I/∂a ${r.dI} vs ${dI}`);
  }
});

test('its exact identities: the mirror about a = λ, slope ½ at a = 0, saturation at √(4πλ)', () => {
  for (const lam of [0.5, 3, 19.5, 40, 300]) {
    const full = Math.sqrt(4 * Math.PI * lam), I0 = mhcIntegral(lam, 0);
    // I(λ, λ + δ) + I(λ, λ − δ) = √(4πλ): the occupied and empty levels share the Gaussian
    for (const d of [0, 1, 0.5 * lam, lam, 1.5 * lam + 3]) {
      const sum = mhcIntegral(lam, lam + d).I + mhcIntegral(lam, lam - d).I;
      assert.ok(Math.abs(sum / full - 1) < 1e-13, `λ ${lam}, δ ${d}: ${sum} vs ${full}`);
    }
    // d ln I/da = ½ at a = 0 (from I(λ, −a) = e^{−a} I(λ, a)), the symmetric Butler–Volmer α
    assert.ok(Math.abs(I0.dI / I0.I - 0.5) < 1e-13, `λ ${lam}: slope ${I0.dI / I0.I}`);
    assert.ok(Math.abs(mhcIntegral(lam, 4 * lam + 200).I / full - 1) < 1e-13, `λ ${lam}: saturation`);
  }
});

test('its transfer coefficient d ln I/da: ½ at a = 0, falling as 1/(2λ) for λ ≫ 1, and slower at moderate λ', () => {
  const al = (lam, a) => mhcIntegral(lam, a).dI / mhcIntegral(lam, a).I;
  const fall = (lam) => ((al(lam, -0.01) - al(lam, 0.01)) / 0.02) * 2 * lam; // in units of Marcus's 1/(2λ)
  assert.ok(Math.abs(fall(400) - 1) < 0.02 && Math.abs(fall(1000) - 1) < 0.01, `${fall(400)}, ${fall(1000)}`);
  assert.ok(Math.abs(fall(0.22 * 38.92) - 0.70) < 0.01, `λ = 0.22 eV: ${fall(0.22 * 38.92)}`); // 30% slower
});

test("its limits: Marcus's parabola for λ ≫ |a|, a single level for λ → 0, and Zeng et al.'s closed form", () => {
  // λ ≫ |a|: I(λ, a)/I(λ, 0) → exp(a/2 − a²/(4λ)), the Fermi sea's spread a correction ~1/λ²
  for (const lam of [400, 1000]) {
    const I0 = mhcIntegral(lam, 0).I;
    for (const a of [-4, -1, 1, 4]) {
      const dev = Math.log(mhcIntegral(lam, a).I / I0) - (a / 2 - (a * a) / (4 * lam));
      assert.ok(Math.abs(dev) < 50 / (lam * lam), `λ ${lam}, a ${a}: off the parabola by ${dev}`);
    }
  }
  // λ → 0: the Gaussian is narrower than the Fermi edge, I → √(4πλ)/(1 + e^{λ − a})
  for (const a of [-3, 0, 2]) {
    const lam = 1e-4, want = Math.sqrt(4 * Math.PI * lam) / (1 + Math.exp(lam - a));
    assert.ok(Math.abs(mhcIntegral(lam, a).I / want - 1) < 1e-4, `a ${a}`);
  }
  // Zeng, Smith, Bai and Bazant, J. Electroanal. Chem. 735, 77 (2014): √(πλ)/(1 + e^{−a})
  // erfc((λ − √(1 + √λ + a²))/(2√λ)), an approximation within 12% here for λ up to 19.5 RT (0.5 eV);
  // past that it runs high at a = 0 (18% for λ = 30 RT, 26% at 1 eV, each also by mpmath's own quad).
  const erfc = (x) => {
    // (Numerical Recipes' erfcc, fractional error < 1.2e−7)
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  };
  const zeng = (lam, a) => (Math.sqrt(Math.PI * lam) / (1 + Math.exp(-a))) * erfc((lam - Math.sqrt(1 + Math.sqrt(lam) + a * a)) / (2 * Math.sqrt(lam)));
  for (const lam of [0.5, 1, 2, 4, 8, 10, 19.46]) {
    for (let a = -40; a <= 80; a += 0.5) assert.ok(Math.abs(zeng(lam, a) / mhcIntegral(lam, a).I - 1) < 0.12, `λ ${lam}, a ${a}`);
  }
  assert.ok(Math.abs(zeng(38.92, 0) / mhcIntegral(38.92, 0).I - 1.26) < 0.005, 'at 1 eV'); // (λ = 1 eV, 38.9 RT)
});

// O⁺ + e⁻ ⇌ R at a platinum electrode behind a 20 µm unstirred layer, in 1 M KCl (as in
// test/kinetics.test.js), its kinetics Marcus–Hush–Chidsey: k0 small enough that the rate
// saturates well below the limiting current, so the plateau is the kinetics' own.
const L = 20e-6, cO = 1, cR = 1, DO = 1e-9, DR = 0.8e-9, lambda = units.eV(0.5);
const mu0 = { 'O+': 10e3, R: -86.485e3, 'K+': -283.3e3, 'Cl-': -131.2e3 };
const electrode = (V, rx) => ({
  species: [
    { name: 'O+', z: 1, cRef: 1000 },
    { name: 'R', z: 0, cRef: 1000 },
    { name: 'K+', z: 1, cRef: 1000 },
    { name: 'Cl-', z: -1, cRef: 1000 },
    { name: 'e-', z: -1 },
  ],
  materials: {
    water: {
      epsr: 78.5,
      species: { 'O+': { D: DO, mu0: mu0['O+'] }, R: { D: DR, mu0: mu0.R }, 'K+': { D: 1.96e-9, mu0: mu0['K+'] }, 'Cl-': { D: 2.03e-9, mu0: mu0['Cl-'] } },
    },
    Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } },
  },
  regions: [{ material: 'Pt', length: 1e-6 }, { material: 'water', length: L }],
  interfaces: [{ phi: 'neutral', reactions: [{ left: { 'e-': -1 }, right: { 'O+': -1, R: 1 }, k0: 1e-5, lambda, ...rx }] }],
  contacts: {
    left: { V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' },
    right: { bath: { c: { 'O+': cO, R: cR, 'K+': 1000, 'Cl-': 1001 }, reference: 'Cl-' } },
  },
  grid: { hmin: 0.05e-9, hmax: 200e-9, ratio: 1.15 },
});
const phiB = (mu0['Cl-'] + RT * Math.log(1001 / 1000)) / FARADAY;
const E0p = (mu0['O+'] - mu0.R) / FARADAY + phiB, Veq = E0p + VT * Math.log(cO / cR);
const kp = 1e-5 / 1000; // k0/c_ref, m/s

test('a Marcus–Hush–Chidsey electrode: none at equilibrium, mixed control against the closed form, the rate saturating below the limiting current', () => {
  const dev = new Device(electrode(Veq));
  assert.ok(Math.abs(dev.solve().current) < 1e-12, 'no current at the Nernst potential');
  // Textbook: r = k′(c_O,s F(x) − c_R,s F(−x)), F = I(λ, x)/I(λ, 0), x = (E°′ − V)/V_T, with
  // linear diffusion to the bath, in closed form (F(−x) = e^{−x} F(x) is detailed balance).
  const lam = lambda / RT, I0 = mhcIntegral(lam, 0).I, F = (x) => mhcIntegral(lam, x).I / I0;
  const theory = (V) => {
    const x = (E0p - V) / VT, f = F(x), b = F(-x);
    return (-FARADAY * kp * (cO * f - cR * b)) / (1 + (kp * L * f) / DO + (kp * L * b) / DR);
  };
  const limit = (-FARADAY * DO * cO) / L, sat = -FARADAY * kp * cO * F(Infinity);
  const rows = [];
  for (const eta of [0.02, 0.1, -0.05, -0.2, -0.4, -0.6, -0.8, -1.0, -1.4]) {
    dev.set({ contacts: { left: { V: Veq + eta } } });
    const sol = dev.solve();
    assert.ok(sol.converged, `η=${eta}`);
    const want = theory(Veq + eta);
    assert.ok(Math.abs(sol.current / want - 1) < 5e-4, `η=${eta}: ${sol.current} vs ${want}`);
    rows.push(sol.current);
  }
  // Far cathodic, the current levels off at the saturated rate k′c F(∞), a quarter of the limiting
  // current here, where Butler–Volmer's would have run on to the limit.
  assert.ok(sat / limit > 0.13 && sat / limit < 0.15, `saturated ${sat} vs limit ${limit}`);
  const plateau = rows.at(-1) / theory(-50);
  assert.ok(Math.abs(rows.at(-1) / rows.at(-2) - 1) < 0.01 && Math.abs(plateau - 1) < 0.01, `plateau ${rows.at(-1)} ${rows.at(-2)}`);
  const bv = new Device(electrode(Veq - 1.4, { lambda: undefined, alpha: 0.5 })).solve().current;
  assert.ok(bv / limit > 0.99, `Butler–Volmer at the same k0: ${bv} vs limit ${limit}`);
});

test('a Marcus–Hush–Chidsey polarization curve: the solved rate at the metal\'s level, and a Tafel slope that bends from α = ½', () => {
  const dev = new Device(electrode(Veq - 0.3));
  const sol = dev.solve();
  const Ve = sol.V['e-'][0];
  const p = polarization(dev, sol, { face: 0 }, [Ve, Veq - 0.4, Veq - 0.401]);
  const r = p.reactions[0];
  assert.ok(Math.abs(r.rate[0] / sol.interfaces[0].rates[0] - 1) < 1e-9, `${r.rate[0]} vs ${sol.interfaces[0].rates[0]}`);
  // Near the couple's level, half the potential's slope (α = ½, as for any λ); 0.4 V down, less
  // (the local α = d ln k/dx, which Marcus's parabola puts at ½ − x/(2λ) for λ ≫ 1, here from the integral itself)
  const lam = lambda / RT, ln = (j, k) => Math.log(Math.abs(r.rate[k] / r.rate[j])) / (Math.abs(p.V[k] - p.V[j]) / VT);
  const surf = (V) => sol.c['O+'][dev.model.grid.regionStart[1]] * mhcIntegral(lam, (E0p - V) / VT).I;
  assert.ok(ln(1, 2) < 0.4 && Math.abs(ln(1, 2) / (Math.log(surf(p.V[2]) / surf(p.V[1])) / (0.001 / VT)) - 1) < 1e-3, `local α ${ln(1, 2)}`);
});

test('Marcus–Hush–Chidsey kinetics is asked for clearly', () => {
  const bad = (rx, re) => assert.throws(() => new Device(electrode(Veq, rx)), (e) => e instanceof DeviceError && re.test(e.message), re);
  bad({ alpha: 0.5 }, /lambda.*or alpha.*not both/);
  bad({ lambda: -1 }, /lambda/);
  bad({ lambda: 1e9 }, /1000 RT/);
  bad({ left: { 'e-': -2 }, right: { 'O+': -2, R: 2 } }, /one electron/);
  const port = electrode(Veq);
  port.materials.water.epsr = 0; // (an electrode port's region is strictly neutral)
  port.ports = [{ name: 'wire', region: 1, from: 0, to: 1e-6, V: 0, terminal: 'e-', area: 1e4, reactions: [{ equation: 'O+ + e- = R', k0: 1e-5, lambda }] }];
  assert.throws(() => new Device(port), (e) => e instanceof DeviceError && /for a face's reactions/.test(e.message));
  // not between a metal's band and a semiconductor's (Marcus's picture is a molecule's or an ion's)
  const band = build({
    T: 298.15, library: [semiconductor('Si', { T: 298.15 }), metal('Pt')],
    stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.3, species: { 'e-': 'blocked' }, reactions: [{ left: { 'e-': -1 }, right: { 'e-': 1 }, k0: 1e-3, lambda }] }, layer('Si', 2e-6, { donors: units.perCm3(1e16) }), ohmic(0.2)],
    grid: { hmin: 1e-9, hmax: 50e-9, ratio: 1.2 },
  });
  assert.throws(() => new Device(band), (e) => e instanceof DeviceError && /not between a metal's band and a semiconductor's/.test(e.message));
  // λ in eV where J/mol was meant: warned of
  const slip = electrode(Veq, { lambda: 0.5 });
  assert.ok(unitWarnings(slip).some((w) => /lambda: 0.5 J\/mol/.test(w)), unitWarnings(slip).join('; '));
  assert.ok(!unitWarnings(electrode(Veq)).some((w) => /lambda/.test(w)));
});

// Beside a diffuse layer (src/diffuse.js), the couple reacts at the plane: the textbook Frumkin-corrected
// rate, from what the solution reports there (each ion at the plane, the drop ψ, the metal's μ̄), and
// polarization() at the metal's own level, for Butler–Volmer too. Dilute and crowded (Bikerman).
test('beside a diffuse layer: the Frumkin-corrected rate at the plane, and polarization() there', () => {
  const fe = { 'Fe3+': 1, 'Fe2+': 5, 'Cl-': 13 };
  const cell = (kin, V, crowd) => {
    const def = build({
      T: 298.15,
      library: [aqueous(['Fe3+', 'Fe2+', 'Cl-'], { epsr: 0 }), metal('Pt')],
      stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.2, diffuse: { epsr: 78.4 }, reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1e-2, ...kin }] }, layer('water', 20e-6, { c0: fe }), bath(fe, 'Cl-')],
      grid: { hmin: 1e-7, hmax: 2e-6 },
    });
    if (crowd) {
      def.species.push({ name: 'S', z: 0, cRef: 1000 });
      def.materials.water.species.S = { D: 1e-9, mu0: 0 };
      for (const r of def.regions) if (r.material === 'water') r.c0 = { ...r.c0, S: 500 };
      def.contacts.right.bath.c.S = 500;
      def.materials.water.statistics = [{ type: 'lattice', species: ['Fe3+', 'Fe2+', 'Cl-', 'S'], cMax: 13000 }];
    }
    return def;
  };
  const lam = units.eV(0.5) / RT, I0 = mhcIntegral(lam, 0).I, F = (x) => mhcIntegral(lam, x).I / I0;
  for (const [V, crowd] of [[0.9, false], [-0.3, false], [0.2, false], [1.4, true]]) {
    for (const kin of [{ lambda: units.eV(0.5) }, { alpha: 0.5 }]) {
      const dev = new Device(cell(kin, V, crowd)), sol = dev.solve();
      assert.ok(sol.converged);
      const gL = dev.model.grid.regionEnd[0], gR = dev.model.grid.regionStart[1], dl = sol.interfaces[0].diffuse.right;
      if (kin.lambda) {
        // a° = (μ°_O + 3Fφ_plane + μ̄_e − μ°_R − 2Fφ_plane)/RT, with μ° + zFφ at the plane its edge's plus zFψ
        const a0 = (sol.muStd['Fe3+'][gR] + 3 * FARADAY * dl.psi + sol.mu['e-'][gL] - sol.muStd['Fe2+'][gR] - 2 * FARADAY * dl.psi) / RT;
        const want = 1e-2 * ((dl.plane['Fe3+'] / 1000) * F(a0) - (dl.plane['Fe2+'] / 1000) * F(-a0));
        assert.ok(Math.abs(sol.interfaces[0].rates[0] / want - 1) < 1e-10, `V ${V}${crowd ? ', crowded' : ''}: ${sol.interfaces[0].rates[0]} vs ${want}`);
      }
      const p = polarization(dev, sol, { face: 0 }, sol.V['e-'][gL]).reactions[0].rate[0];
      assert.ok(Math.abs(p / sol.interfaces[0].rates[0] - 1) < 1e-9, `${kin.lambda ? 'MHC' : 'BV'}, V ${V}${crowd ? ', crowded' : ''}: polarization ${p} vs ${sol.interfaces[0].rates[0]}`);
    }
  }
});
