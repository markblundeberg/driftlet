// Doped Si in a resolved electrolyte. A soluble one-electron couple can exchange
// conduction electrons at the face. Its concentrations, formal potential and band
// alignment are model inputs, not material data for a particular redox couple.
import { EPS0, AVOGADRO, FARADAY as F, GAS_CONSTANT } from '../../../src/index.js';
import { semiconductor, ohmic, bath, build, layer } from '../../../src/kit.js';

export default {
  name: 'wetjunction',
  make(r) {
    const T = 298.15, RT = GAS_CONSTANT * T, VT = RT / F;
    const si = semiconductor('Si', { T }), bands = si.materials.Si.species;
    const N = r.log(-1.5, 0), c = r.log(-0.3, 1);
    const L = r.u(60, 150) * 1e-9, W = r.u(30, 90) * 1e-9;
    const model = r.pick(['ideal', 'lattice', 'debye-huckel']);
    const reactive = r.chance(0.72), compact = r.chance(0.5);
    const cMax = r.u(2500, 6000), a = r.u(0.35, 0.55) * 1e-9;
    const composition = { 'Ox+': c, 'A-': c, ...(reactive ? { Red: c } : {}) };
    const ions = Object.keys(composition);
    const species = Object.fromEntries(ions.map((s) => [s, { D: r.log(-9.3, -8.85), mu0: 0, cRef: c }]));
    const statistics = model === 'ideal' ? [] : [model === 'lattice' ? { type: model, species: ions, cMax } : { type: model, species: ['Ox+', 'A-'], epsr: 78.3, a }];
    // At equal concentrations the lattice excess cancels between Ox+ and Red.
    // In Debye-Hückel only the charged reactant contributes an excess to E_rest.
    const eps = 78.3 * EPS0, kappa = Math.sqrt(2 * c * F * F / (eps * RT));
    const lB = F * F / (4 * Math.PI * eps * AVOGADRO * RT);
    const rest = reactive && model === 'debye-huckel' ? -VT * lB * kappa / (2 * (1 + a * kappa)) : 0;
    const ni2 = bands['e-'].cRef * bands['h+'].cRef * Math.exp(-bands['h+'].mu0 / RT);
    const n = N / 2 + Math.sqrt(N * N / 4 + ni2);
    const dipole = -rest - VT * Math.log(n / bands['e-'].cRef) + r.u(-0.035, 0.035);
    const C = r.u(0.08, 0.35), k0 = r.log(-5.5, -3.5), alpha = r.u(0.35, 0.65);
    const face = { dipole, ...(compact ? { phi: { type: 'capacitive', C } } : {}),
      ...(reactive ? { reactions: [{ equation: 'Ox+ + e- = Red', k0, alpha }] } : {}) };
    const def = build({
      T, library: [si], species: ions.map((name) => ({ name, z: name === 'Ox+' ? 1 : name === 'A-' ? -1 : 0 })),
      materials: { liquid: { epsr: 78.3, species, statistics } },
      stack: [ohmic(rest), layer('Si', L, { donors: N }), face,
        layer('liquid', W, { c0: composition }), bath(composition, 0)],
      grid: { hmax: Math.min(L, W) / 7, hmin: 1e-9, ratio: 1.45 },
    });
    const features = ['n-Si', model, reactive ? 'band redox' : 'blocking junction', compact ? 'Stern capacitance' : 'continuous displacement'];
    return { def, params: { features, N, L, W, c, dipole, rest, ...(compact ? { C } : {}), ...(reactive ? { k0, alpha } : {}),
      ...(model === 'lattice' ? { cMax } : model === 'debye-huckel' ? { a } : {}), D: Object.fromEntries(ions.map((s) => [s, species[s].D])) },
      plan: { terminal: 'left', V: rest + r.u(-0.05, 0.05), dV: r.u(0.008, 0.022), equilibrium: true, unique: true,
        slow: Math.max(W * W / Math.min(...ions.map((s) => species[s].D)), L * L / bands['h+'].D) } };
  },
};
