// A mixed ionic/electronic film with a mobile positive defect on a compensating
// background. One electrode can exchange the defect, as at a reversible electrode,
// or both can block it. The two electronic contacts may be carrier selective.
import { FARADAY, GAS_CONSTANT } from '../../../src/index.js';
import { photogeneration } from '../../../src/kit.js';

export default {
  name: 'mixedfilm',
  make(r) {
    const T = 300, RT = GAS_CONSTANT * T;
    const L = r.u(60, 150) * 1e-9, D = r.log(-12, -10);
    const c = r.log(-1, 0.7), cMax = c * r.u(3, 20);
    const lattice = r.chance(0.55), fermi = r.chance(0.4), open = r.chance(0.6);
    const selective = r.chance(0.5), lit = r.chance(0.5);
    const n = r.log(-2, -0.5), band = r.u(0.09, 0.17) * FARADAY;
    const tau = r.log(-8, -6), flux = r.log(-6, -4), alpha = r.u(0.4, 2) / L;
    const toward = r.pick(['+x', '-x']);
    const statistics = [];
    if (fermi) statistics.push({ type: 'fermi-dirac', species: ['e-', 'h+'] });
    if (lattice) statistics.push({ type: 'lattice', species: ['V+'], cMax });
    // At phi=0 and equal terminal voltages both bands have equal occupations.
    // Invert the vacancy lattice at c so that the optional ion contact is at rest.
    const species = {
      'e-': { D: r.log(-7, -5), mu0: band, cRef: n * Math.exp(band / RT) },
      'h+': { D: r.log(-7, -5), mu0: band, cRef: n * Math.exp(band / RT) },
      'V+': { D, mu0: lattice ? RT * Math.log1p(-c / cMax) : 0, cRef: c },
    };
    const contact = (names) => ({
      V: 0, phi: 'pinned', zeroCharge: 0, terminal: names[0],
      species: Object.fromEntries(names.map((s, k) => [s, k ? { type: 'equilibrium', offset: 0 } : 'equilibrium'])),
    });
    const bulkReactions = [{ equation: 'e- + h+ = 0', kf: { film: 1 / (tau * n) } }];
    if (lit) bulkReactions.push(photogeneration({ material: 'film', flux, alpha, toward, to: L, mu: 2 * FARADAY }));
    const def = {
      T, species: Object.keys(species).map((name) => ({ name, z: name === 'e-' ? -1 : 1 })),
      materials: { film: { epsr: r.u(15, 35), species, statistics } },
      regions: [{ material: 'film', length: L, fixedCharge: -FARADAY * c, c0: { 'V+': c } }],
      contacts: { left: contact(['e-', ...(selective ? [] : ['h+']), ...(open ? ['V+'] : [])]), right: contact(selective ? ['h+'] : ['e-', 'h+']) },
      grid: { hmax: L / 12 }, bulkReactions,
    };
    const features = [open ? 'exchangeable vacancy' : 'closed vacancy', selective ? 'selective contacts' : 'ambipolar contacts',
      lattice ? 'vacancy lattice' : 'dilute vacancy', fermi ? 'Fermi bands' : 'Boltzmann bands', lit ? 'Beer-Lambert light' : 'dark'];
    return { def, params: { features, L, D, c, cMax, n, band, tau, ...(lit ? { flux, alpha, toward } : {}),
      epsr: def.materials.film.epsr, Dn: species['e-'].D, Dp: species['h+'].D },
      plan: { terminal: 'right', V: r.u(-0.06, 0.06), dV: r.u(0.008, 0.025), equilibrium: !lit, unique: open, slow: Math.max(L * L / D, tau) } };
  },
};
