// Strictly neutral charge as a step will start from it: the neutralG property's two readings
// (bench/combos/run.js), kept apart so a script can read a device with them too.
import { FARADAY } from '../../src/constants.js';

// --- strictly neutral charge, as a step will start from it
// Each node of a strictly neutral (ε = 0) region where φ is defined, but at a contact that sets φ:
// its Gauss row, assembled at the state (steady off, so no law takes a row's place), against the
// row's own terms. That row is G(u) = −v(ρ_fixed + F Σ z c) − v·(the ports' σ) + what the box
// holds of a face's or contact's D. Where the contact holds the node neutral by a bulk law, the
// ions alone. The worst |G|/scale, and where.
export function gaussOff(d) {
  const s = d.solver, { n, M, z, c, present } = s, grid = s.model.grid, R = s.rix;
  const nodes = [];
  for (let g = 0; g < s.nNodes; g++) {
    const mat = s.model.materials[s.nodeMaterial[g]];
    if (mat.conductor || mat.epsr !== 0 || s.phiUndefined[g] || s.nodeConductor[g] >= 0) continue;
    const side = g === 0 ? 'left' : g === s.nNodes - 1 ? 'right' : null;
    if (side && s.model.contacts[side].phi.type === 'pinned') continue;
    nodes.push([g, side && s.model.contacts[side].phi.type === 'bulk']);
  }
  if (!nodes.length) return null;
  const steady = s.steady;
  s.steady = false;
  try {
    s.assemble(Infinity);
  } finally {
    s.steady = steady;
  }
  let worst = 0, at = '';
  for (const [g, bulk] of nodes) {
    const v = grid.vol[g], b = s.blockOfNode[g];
    let ions = s.rhoFixed[g], gross = Math.abs(s.rhoFixed[g]);
    for (let i = 0; i < n; i++) if (present[g * n + i]) (ions += FARADAY * z[i] * c[g * n + i]), (gross += FARADAY * Math.abs(z[i] * c[g * n + i]));
    const G = bulk ? v * ions : s.loc[b * M] >= 0 ? s.res[R[b * M]] : NaN;
    const scale = v * (gross + Math.abs(s.qg[g])) + Math.abs(G + v * (ions + s.qg[g]));
    const rel = Math.abs(G) / Math.max(scale, 1e-300);
    if (!(rel <= worst)) (worst = rel), (at = `x = ${grid.x[g].toPrecision(4)} m (node ${g}${bulk ? ', bulk contact' : ''}): G ${G.toExponential(2)} of ${scale.toExponential(2)}`);
  }
  return { worst, at };
}

// What the next step takes as each gated stretch's and each capacitive contact's charge at its
// start (portQEnd, contactDEnd, where a step or set() left them) against what the state holds
// there (its ports' charges as read now; the contact's D): summed over each group of capacitive
// ports sharing strictly neutral boxes (src/regions.js, _gateGroupsOf), where a gate's charge can
// move to another gate on the same boxes, but not out of them. Relative to the charges involved.
export function startCharges(d) {
  const s = d.solver;
  if (!s.portQEnd && !s.contactDEnd) return { worst: 0, at: '' };
  s._assembleBookkeeping(Infinity);
  let worst = 0, at = '';
  const see = (rel, what) => !(rel <= worst) && ((worst = rel), (at = what));
  if (s.portQEnd) for (const { ports, nodes } of s.gateGroups) {
    let start = 0, now = 0, gross = 0;
    for (const k of ports) (start += s.portQEnd[k]), (now += s.portQ[k]), (gross += Math.abs(s.portQEnd[k]) + Math.abs(s.portQ[k]));
    let ions = 0;
    for (const g of nodes) {
      let q = Math.abs(s.rhoFixed[g]);
      for (let i = 0; i < s.n; i++) if (s.present[g * s.n + i]) q += FARADAY * Math.abs(s.z[i] * s.c[g * s.n + i]);
      ions += s.model.grid.vol[g] * q;
    }
    see(Math.abs(start - now) / Math.max(gross + 1e-4 * ions, 1e-300), `gates ${ports.map((k) => s.model.ports[k].name).join(', ')}: start ${start.toExponential(3)}, the state's ${now.toExponential(3)}`);
  }
  return { worst, at };
}

