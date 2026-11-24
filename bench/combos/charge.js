// Strictly neutral charge as a step will start from it: the neutralG property's two readings
// (bench/combos/run.js), kept apart so a script can read a device with them too.
import { FARADAY } from '../../src/constants.js';

// The signed charge a first step from a solved state will read as current through two bulk
// baths. Summing the charged balances and using Gauss at the end leaves the start's Gauss
// residual, divided by dt. Ungated charge rows and the bath end boxes already cancel their
// own start charge. The other boxes still carry their rounded concentrations: 2.6e-15 C/m² at
// a capacitive face and −2.0e-15 C/m² in a resolved membrane read 3.0e-8 and −1.17e-7 A/m²
// (membrane combos 18 and 20, kirchhoff; Sol's chase, 23 Nov).
// Called on a second device so the check reads the same start without touching the step's
// state. Restrict the reading to bulk baths with only volume storage, where the sum above
// applies. A real charge defect is a failure: cap both the local error (forming concentrations
// from levels can lose several bits) and its sum relative to the stored ions' charge.
export function bathStartRoundoff(d) {
  if (!d) return 0;
  const s = d.solver, { model, n, z } = s;
  if (model.ports.length || model.surfaces.length || s.diffList.length || s.nodeConductor.some((i) => i >= 0)) return 0;
  if (Object.values(model.contacts).some((ct) => ct.phi.type !== 'bulk' || ct.species.some((link) => link.type !== 'equilibrium'))) return 0;
  d._settleRegrid();
  // A drive-only set() keeps the solved state, including its storage, at the previous drive.
  if (!s.atSteady && !s.readV) return 0;
  s.assemble(Infinity);
  let defect = 0, gross = 0;
  for (let g = 1; g < s.nNodes - 1; g++) {
    if (s.phiUndefined[g] || s.chargeNode[g]) continue;
    const v = model.grid.vol[g], G = s.res[s.rix[s.blockOfNode[g] * s.M]];
    let ions = 0;
    for (let i = 0; i < n; i++) ions += Math.abs(z[i] * s.c[g * n + i]);
    const scale = v * FARADAY * ions;
    if (!(Math.abs(G) <= 128 * Number.EPSILON * scale)) return 0;
    defect += G;
    gross += scale;
  }
  return Math.abs(defect) <= 16 * Number.EPSILON * gross ? defect : 0;
}

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

// What the next step takes as each capacitor's charge at its start (portQEnd, contactDEnd, where a
// step or set() left them) against what the state holds there (its ports' charges as read now; a
// contact's D), summed over each group of strictly neutral boxes that capacitors join: a port's
// window (the boxes it holds charge in), a capacitive face between two such boxes. Within a group a
// charge can pass from one capacitor to another through the solution, but not out of it, so only
// the sums must agree (a gate whose hardware changed keeps its own charge while φ̂ jumps under its
// neighbours). Relative to the charges involved. (Built here from the definition, apart from the
// solver's own groups: the gates over a capacitive face's box were dropped from those, and from
// this reader with them; Neon's fourth review, R4 1.2, and its contacts weren't read at all, R4 1.1.)
export function startCharges(d) {
  const s = d.solver, grid = s.model.grid, N = s.nNodes;
  if (!s.portQEnd && !s.contactDEnd) return { worst: 0, at: '' };
  s._assembleBookkeeping(Infinity);
  const side = (g) => (g === 0 ? 'left' : g === N - 1 ? 'right' : null);
  const neutral = new Uint8Array(N);
  for (let g = 0; g < N; g++) {
    const mat = s.model.materials[s.nodeMaterial[g]];
    if (mat.conductor || mat.epsr !== 0 || s.phiUndefined[g] || s.nodeConductor[g] >= 0) continue;
    const ct = side(g) && s.model.contacts[side(g)];
    if (ct && (ct.phi.type === 'pinned' || ct.phi.type === 'bulk')) continue;
    neutral[g] = 1;
  }
  const root = Int32Array.from({ length: N }, (_, g) => g);
  const find = (g) => (root[g] === g ? g : (root[g] = find(root[g])));
  const join = (a, b) => (root[find(a)] = find(b));
  const held = s.model.ports.map((p) => (p.capacitance ? p.nodes.filter((g) => neutral[g] && s._portHolds(g)) : []));
  for (const ns of held) for (const g of ns) join(g, ns[0]);
  s.model.interfaces.forEach((itf, f) => {
    const gL = grid.regionEnd[f];
    if (itf.phi.type === 'capacitive' && neutral[gL] && neutral[gL + 1]) join(gL, gL + 1);
  });
  const groups = new Map(), at = (g) => groups.get(find(g)) ?? groups.set(find(g), { start: 0, now: 0, gross: 0, names: [], nodes: new Set() }).get(find(g));
  held.forEach((ns, k) => {
    if (ns.length === 0 || !s.portQEnd) return;
    const gr = at(ns[0]);
    gr.start += s.portQEnd[k];
    gr.now += s.portQ[k];
    gr.gross += Math.abs(s.portQEnd[k]) + Math.abs(s.portQ[k]);
    gr.names.push(s.model.ports[k].name);
  });
  for (const g of [0, N - 1]) {
    const sd = side(g);
    if (!neutral[g] || !s.contactDEnd || s.model.contacts[sd].phi.type !== 'capacitive') continue;
    // (its plate's charge, as a port's σ is its electrode's: contactD is the displacement toward +x)
    const gr = at(g), sg = sd === 'left' ? 1 : -1;
    gr.start += sg * s.contactDEnd[sd];
    gr.now += sg * s.contactD[sd];
    gr.gross += Math.abs(s.contactDEnd[sd]) + Math.abs(s.contactD[sd]);
    gr.names.push(`contact ${sd}`);
  }
  // (a face to a metal whose sheet set() left to start from a D of its own, sheetDEnd: its share,
  // −A·D on the box to its left, +A·D on the one to its right)
  s.model.interfaces.forEach((itf, f) => {
    const gL = grid.regionEnd[f], D0 = s.sheetDEnd?.[f];
    if (itf.phi.type !== 'capacitive' || !Number.isFinite(D0) || neutral[gL] === neutral[gL + 1]) return;
    const g = neutral[gL] ? gL : gL + 1, sg = neutral[gL] ? -1 : 1, A = grid.area[gL], D = s.u[s.blockOfFace[f] * s.M];
    const gr = at(g);
    gr.start += sg * A * D0;
    gr.now += sg * A * D;
    gr.gross += A * (Math.abs(D0) + Math.abs(D));
    gr.names.push(`face ${f}`);
  });
  for (let g = 0; g < N; g++) if (neutral[g] && groups.has(find(g))) groups.get(find(g)).nodes.add(g);
  let worst = 0, where = '';
  for (const gr of groups.values()) {
    let ions = 0;
    for (const g of gr.nodes) {
      let q = Math.abs(s.rhoFixed[g]);
      for (let i = 0; i < s.n; i++) if (s.present[g * s.n + i]) q += FARADAY * Math.abs(s.z[i] * s.c[g * s.n + i]);
      ions += grid.vol[g] * q;
    }
    const rel = Math.abs(gr.start - gr.now) / Math.max(gr.gross + 1e-4 * ions, 1e-300);
    if (!(rel <= worst)) (worst = rel), (where = `${gr.names.join(', ')}: start ${gr.start.toExponential(3)}, the state's ${gr.now.toExponential(3)}`);
  }
  return { worst, at: where };
}
