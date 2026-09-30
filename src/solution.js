// Plotting-ready snapshot of a solver state.
//
// Arrays are indexed by grid node; doubled interface nodes appear twice at the same x, so
// steps plot as vertical lines. Quantities of a species are NaN where it is absent, and
// the voltage views are NaN for neutral species.

import { FARADAY } from './constants.js';

export function makeSolution(solver, result = {}) {
  const { model, n, M, u, c, VT } = solver;
  const { grid, species, contacts, interfaces } = model;
  const nNodes = grid.nNodes;
  const RT = model.RT, F = FARADAY;

  const phi = new Float64Array(nNodes);
  for (let g = 0; g < nNodes; g++) phi[g] = VT * u[solver.blockOfNode[g] * M];

  const sol = {
    x: Float64Array.from(grid.x),
    region: Int32Array.from(grid.nodeRegion),
    phi,
    c: {},
    mu: {},
    muStd: {},
    V: {},
    Vstd: {},
    time: solver.time,
    converged: result.converged ?? true,
    iterations: result.iterations ?? 0,
    steps: result.steps,
    history: result.history,
  };

  for (let i = 0; i < n; i++) {
    const { name, z } = species[i];
    const cc = new Float64Array(nNodes), mu = new Float64Array(nNodes), mus = new Float64Array(nNodes);
    const V = new Float64Array(nNodes), Vs = new Float64Array(nNodes);
    for (let g = 0; g < nNodes; g++) {
      const k = g * n + i;
      if (!solver.present[k]) {
        cc[g] = mu[g] = mus[g] = V[g] = Vs[g] = NaN;
        continue;
      }
      cc[g] = c[k];
      const o = solver.blockOfNode[g] * M + 1 + i;
      mu[g] = RT * (u[o] + solver.uLo[o]);
      mus[g] = RT * solver.mu0hat[k] + z * F * phi[g];
      V[g] = z === 0 ? NaN : mu[g] / (z * F);
      Vs[g] = z === 0 ? NaN : mus[g] / (z * F);
    }
    sol.c[name] = cc;
    sol.mu[name] = mu;
    sol.muStd[name] = mus;
    sol.V[name] = V;
    sol.Vstd[name] = Vs;
  }

  // Contacts: fluxes and displacement at the final state (re-evaluated from the balance rows).
  solver.assemble(solver.lastDt);
  sol.contacts = {};
  sol.gates = {};
  for (const side of ['left', 'right']) {
    const ct = contacts[side];
    const flux = {};
    let conduction = 0;
    for (let i = 0; i < n; i++) {
      flux[species[i].name] = solver.contactFlux[side][i];
      conduction += F * species[i].z * solver.contactFlux[side][i];
    }
    const D = solver.contactD[side];
    const dt = solver.lastDt;
    const displacement = Number.isFinite(dt) ? (D - solver.contactDOld[side]) / dt : 0;
    // Current toward +x through this contact; in steady state both contacts agree.
    sol.contacts[side] = { V: ct.V, flux, D, current: conduction + displacement };
    if (ct.phi.type === 'capacitive') {
      // Charge per area on the gate (or metal) plate: +D at the left, −D at the right.
      sol.gates[side] = { V: ct.V, D, charge: side === 'left' ? D : -D };
    }
  }
  sol.current = sol.contacts.right.current;
  sol.terminalVoltage = contacts.right.V - contacts.left.V;

  // Interfaces: dipole and the fluxes carried by each flux node.
  sol.interfaces = interfaces.map((itf, f) => {
    const b = solver.blockOfFace[f];
    const N = {};
    for (let i = 0; i < n; i++) N[species[i].name] = u[b * M + 1 + i];
    return { dipole: itf.dipole, sheetCharge: itf.sheetCharge, D: u[b * M], N };
  });

  // Total charge per area in the device (space charge plus sheet charges).
  let q = 0;
  for (let g = 0; g < nNodes; g++) {
    let rho = solver.rhoFixed[g];
    for (let i = 0; i < n; i++) rho += F * species[i].z * c[g * n + i];
    q += grid.vol[g] * rho;
  }
  for (const itf of interfaces) q += itf.sheetCharge;
  sol.charge = q;

  // Conservation bookkeeping for each species stretch.
  sol.conservation = solver.stretches.map((st, k) => {
    const amount = solver.amount(st);
    const reference = solver.referenceAmounts[k];
    return {
      species: species[st.species].name,
      regions: st.regions.slice(),
      spectator: st.spectator,
      amount,
      reference,
      drift: st.spectator ? (amount - reference) / reference : NaN,
    };
  });

  return sol;
}
