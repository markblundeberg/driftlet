// Two plain devices, to try the properties on before the feature families: a silicon pn diode
// with SRH recombination, and a silver nitrate cell between silver electrodes.
import { units } from '../../../src/index.js';
import { build, layer, ohmic, semiconductor } from '../../../src/kit.js';

export default {
  name: 'basic',
  make(r) {
    if (r.chance(0.5)) {
      const lib = semiconductor('Si'), L = r.log(-6.5, -5.5), tau = r.log(-8, -6);
      const N = units.perCm3(r.log(15, 17));
      const def = build({
        T: 300,
        library: [lib],
        stack: [ohmic(0), layer('Si', L, { name: 'p', acceptors: N }), layer('Si', L, { name: 'n', donors: N }), ohmic(0)],
        bulkReactions: [{ equation: 'e- + h+ = 0', srh: { Si: { tauN: tau, tauP: tau } } }],
        grid: { hmin: 1e-9, hmax: L / 40, ratio: 1.15 },
      });
      return { def, params: { kind: 'diode', L, tau, N, features: ['srh'] }, plan: { terminal: 'right', V: r.u(-0.6, 0.3), dV: 0.05, equilibrium: true, unique: true, slow: Math.max(tau, (4 * L * L) / 3.5e-3) } };
    }
    const c = r.log(0, 2), L = r.log(-5, -4);
    const electrode = { V: 0, terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
    const def = {
      T: 298.15,
      species: [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }],
      materials: { water: { epsr: 0, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } } },
      regions: [{ name: 'w', material: 'water', length: L, c0: { 'NO3-': c } }],
      contacts: { left: { ...electrode, species: { ...electrode.species } }, right: { ...electrode, species: { ...electrode.species } } },
      grid: { hmin: L / 2000, hmax: L / 50 },
    };
    // (NO3⁻ is blocked at both ends: its amount is closed, so the state isn't unique)
    const plan = { terminal: 'right', V: r.u(-0.05, 0.05), dV: 0.01, equilibrium: true, unique: false, slow: (L * L) / 1.65e-9 }, features = ['blocked'];
    if (r.chance(0.4)) {
      // Two ports pumping NO₃⁻ in and out at ±J, its amount closed between them: a law whose
      // change two driven terminals share. (Neon's second review, R2 1.1: impedance() left such
      // a port with no admittance, Re Z −19.3 Ω·m² where held it was 0.0036, unwarned.)
      const sigma = (96485.33212 ** 2 * c * (1.65e-9 + 1.9e-9)) / (8.314462618 * 298.15), area = r.log(2, 4), C = 0.2;
      const J = r.u(0.1, 1) * ((0.2 * c * 96485.33212 * 1.9e-9) / L);
      def.ports = [0, 1].map((k) => ({ name: `p${k + 1}`, region: 'w', from: (k ? 0.6 : 0.1) * L, to: (k ? 0.9 : 0.4) * L, area, capacitance: { C, zeroCharge: 0 }, I: k ? J : -J, species: { 'NO3-': { type: 'conductance', G: 1e8, offset: 0 } } }));
      plan.slow = Math.max(plan.slow, (C * area * L * L) / sigma);
      plan.equilibrium = false;
      features.push('pumps');
    }
    return { def, params: { kind: 'cell', c, L, features }, plan };
  },
};
