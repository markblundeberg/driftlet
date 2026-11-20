// Electrochemical gating of a redox-active film: electrons supplied by the
// back contact reduce a mobile dopant Ox+ to its neutral form Red. The gate is
// insulated. Ox+ and Red can share sites and have different diffusivities.
import { FARADAY as F, GAS_CONSTANT } from '../../../src/index.js';
import { photogeneration } from '../../../src/kit.js';

export default {
  name: 'redoxgate',
  make(r) {
    const T = 300, RT = T * GAS_CONSTANT;
    const L = r.u(50, 150) * 1e-9, n = r.log(-1, 0.5), total = r.log(-0.3, 1);
    const fraction = r.u(0.2, 0.8), ox = fraction * total, red = total - ox;
    const lattice = r.chance(0.5), open = r.chance(0.35), compact = r.chance(0.6);
    const cMax = total * r.u(2, 10), C = r.u(0.02, 0.15);
    const Dox = r.log(-13, -10), Dred = Dox * r.u(0.3, 3), Dn = r.log(-7, -5);
    const epsr = r.u(20, 50), kf = r.log(4, 8);
    // At phi = 0 all three levels are zero at the concentrations given.
    // A common lattice excess cancels between the two redox states.
    const mu0 = lattice ? RT * Math.log1p(-total / cMax) : 0;
    const species = {
      'e-': { D: Dn, mu0: RT * Math.log(50), cRef: 50 * n },
      'Ox+': { D: Dox, mu0, cRef: ox },
      Red: { D: Dred, mu0, cRef: red },
    };
    const def = { T,
      species: [{ name: 'e-', z: -1 }, { name: 'Ox+', z: 1 }, { name: 'Red', z: 0 }],
      materials: { film: { epsr, species,
        statistics: lattice ? [{ type: 'lattice', species: ['Ox+', 'Red'], cMax }] : [] } },
      regions: [{ material: 'film', length: L, fixedCharge: F * (n - ox), c0: { 'e-': n, 'Ox+': ox, Red: red } }],
      bulkReactions: [{ equation: 'Ox+ + e- = Red', kf: { film: kf } }],
      contacts: { left: { V: 0, phi: 'pinned', zeroCharge: 0, terminal: 'e-',
        species: { 'e-': 'equilibrium', ...(open ? { Red: { type: 'equilibrium', mu: 0 } } : {}) } },
        right: { V: 0, phi: compact ? { type: 'capacitive', C } : 'pinned', zeroCharge: 0 } },
      grid: { hmax: L / 12 },
    };
    const features = ['bulk electron transfer', open ? 'buffered Red' : 'closed dopant',
      lattice ? 'shared lattice' : 'dilute dopant', compact ? 'capacitive gate' : 'pinned gate'];
    const plan = { terminal: 'right', V: r.u(-0.06, 0.06), dV: r.u(0.01, 0.025),
      equilibrium: true, unique: open, slow: Math.max(L * L / Math.min(Dox, Dred), 1 / (kf * n)) };
    const params = { features, L, n, total, fraction, Dox, Dred, Dn, epsr, kf,
      ...(lattice ? { cMax } : {}), ...(compact ? { C } : {}) };
    // In a closed film, light can cycle the finite dopant pool. Ionization has
    // first-order dependence on Red. Calibrate the profile to the initial Red
    // density, so the initial source integrates to the absorbed photon flux.
    if (!open && r.chance(0.6)) {
      const flux = r.log(-5, -2.8), alpha = r.u(0.4, 2) / L;
      const light = photogeneration({ material: 'film', flux: flux / red, alpha, to: L, mu: 2 * F });
      light.equation = 'photon + Red = Ox+ + e-';
      def.bulkReactions.push(light);
      features.push('dopant photoionization');
      Object.assign(params, { flux, alpha });
      plan.equilibrium = false;
    }
    return { def, params, plan };
  },
};
