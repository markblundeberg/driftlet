// Plotting-ready snapshot of a solver state.
//
// Arrays are indexed by grid node; doubled interface nodes appear twice at the same x, so
// steps plot as vertical lines. Quantities of a species are NaN where it is absent, and
// the voltage views are NaN for neutral species.

import { EPS0, FARADAY } from './constants.js';

export function makeSolution(solver, result = {}) {
  const { model, n, M, u, c, VT } = solver;
  const { grid, species, contacts, interfaces } = model;
  const nNodes = grid.nNodes;
  const RT = model.RT, F = FARADAY;

  const phi = new Float64Array(nNodes);
  for (let g = 0; g < nNodes; g++) phi[g] = solver.phiUndefined[g] ? NaN : VT * u[solver.blockOfNode[g] * M];

  const sol = {
    x: Float64Array.from(grid.x),
    region: Int32Array.from(grid.nodeRegion),
    species: species.map(({ name, z }) => ({ name, z })),
    regions: model.regions.map((r, k) => ({ name: r.name, material: model.materials[r.material].name, x0: grid.x[grid.regionStart[k]], x1: grid.x[grid.regionEnd[k]] })),
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
    substeps: result.substeps,
    history: result.history,
  };
  if (result.done !== undefined) {
    sol.done = result.done;
    sol.rejected = result.rejected;
    sol.trace = result.trace;
  }

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
      const o = solver.blockOfNode[g] * M + 1 + i;
      mu[g] = RT * (u[o] + solver.uLo[o]);
      if (solver.nodeConductor[g] >= 0) {
        // A metal has a Fermi level but no concentration or standard level of its own.
        cc[g] = mus[g] = Vs[g] = NaN;
        V[g] = mu[g] / (z * F);
        continue;
      }
      cc[g] = c[k];
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

  // Contacts: fluxes and displacement at the final state (re-evaluated from the balance rows of
  // the boxes they're read from).
  solver.computeConcentrations();
  solver._assembleBookkeeping(solver.lastDt);
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
    const V = solver.termV[side === 'left' ? 0 : 1];
    sol.contacts[side] = { V, flux, D, current: conduction + displacement };
    if (ct.phi.type === 'capacitive' || ct.phi.type === 'pinned') {
      // Charge per area on the gate (or metal) plate: +D at the left, −D at the right.
      sol.gates[side] = { V, D, charge: side === 'left' ? D : -D };
    }
  }
  sol.current = sol.contacts.right.current;
  // Internal ports: what each brings into the device (in steady state, the contact currents
  // differ by the ports' total).
  sol.ports = model.ports.map((port, k) => {
    const flux = {};
    let current = 0;
    for (let i = 0; i < n; i++) {
      flux[species[i].name] = solver.portFlux[k][i];
      current += F * species[i].z * solver.portFlux[k][i];
    }
    return { name: port.name, V: solver.termV[2 + k], flux, current };
  });
  sol.terminalVoltage = sol.contacts.right.V - sol.contacts.left.V;
  // Every terminal (the contacts, then the ports by name): its voltage and its current into the
  // device, conduction plus displacement. They sum to zero in steady state.
  sol.terminals = {};
  model.terminals.forEach((t, k) => {
    const current = t.kind === 'port' ? sol.ports[t.index].current : t.side === 'left' ? sol.contacts.left.current : -sol.contacts.right.current;
    sol.terminals[t.name] = { V: solver.termV[k], current };
  });

  // Interfaces: dipole and the fluxes carried by each flux node.
  sol.interfaces = interfaces.map((itf, f) => {
    const b = solver.blockOfFace[f];
    const N = {};
    for (let i = 0; i < n; i++) N[species[i].name] = u[b * M + 1 + i];
    const rates = itf.reactions.map((_, k) => u[b * M + 1 + n + k]); // mol/(m²·s), forward
    return { left: model.regions[f].name, right: model.regions[f + 1].name, dipole: itf.dipole, sheetCharge: itf.sheetCharge, D: u[b * M], N, rates };
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

  // Resolution warnings: where the model resolves a double layer (dipole or capacitive faces,
  // gate contacts), check the local Debye length against the adjacent cell.
  sol.warnings = [...(model.warnings ?? [])];
  const cond = solver.conditioning;
  if (result.converged === false && cond && cond.digits > 12) {
    const digits = Number.isFinite(cond.digits) ? `lost ${cond.digits.toFixed(0)} of its ~16 digits` : 'was exactly singular';
    sol.warnings.push(
      `the steady system ${digits} near x = ${cond.x.toExponential(3)} m (${cond.where}): part of the device is held ` +
        'only weakly, e.g. a floating region coupled to the rest through tiny conductances or rates, or a stiff ' +
        'chain whose level nothing pins. Strengthen that coupling, or anchor the region (a port, a contact).',
    );
  }
  const debye = (g) => {
    const mat = model.materials[model.regions[grid.nodeRegion[g]].material];
    const s2 = solver.screening(g); // zᵀKz: Σ z² c for ideal statistics
    return mat.epsr > 0 && s2 > 0 ? Math.sqrt((mat.epsr * EPS0 * RT) / (F * F * s2)) : NaN;
  };
  const check = (g, h, where) => {
    const lam = debye(g);
    if (lam < h) {
      sol.warnings.push(
        `${where}: double layer unresolved (cell ${h.toExponential(2)} m vs Debye length ${lam.toExponential(2)} m); ` +
          "its charge will depend on the mesh. Refine the grid there, or use phi: 'neutral' for a macroscopic model.",
      );
    }
  };
  interfaces.forEach((itf, f) => {
    if (itf.phi.type === 'neutral') return;
    const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
    check(gL, grid.segLength[gL - 1], `interfaces[${f}] (left side)`);
    check(gR, grid.segLength[gR], `interfaces[${f}] (right side)`);
  });
  const resolves = (ct) => ct.phi.type === 'capacitive' || ct.phi.type === 'pinned';
  if (resolves(contacts.left)) check(0, grid.segLength[0], 'contacts.left');
  if (resolves(contacts.right)) check(nNodes - 1, grid.segLength[nNodes - 2], 'contacts.right');

  // Steep profiles carrying current: Scharfetter–Gummel takes the field as uniform across each
  // cell, which fails where it isn't, as at an electrode where a species is nearly depleted (there
  // φ̂ ≈ ln c, and on a coarse cell the flux comes out too large: by (1+r)ln(1/r)/(2(1−r)) for
  // a concentration ratio r across it). Flag a region's end cell whose field is well above its
  // neighbour's while its concentration changes steeply, for species carrying the current. (A
  // pn depletion region's field changes smoothly from cell to cell, which SG handles.)
  if (Number.isFinite(sol.current)) {
    const share = species.map((sp, i) =>
      Math.max(Math.abs(sp.z * (sol.contacts.left.flux[sp.name] ?? 0)), Math.abs(sp.z * (sol.contacts.right.flux[sp.name] ?? 0))),
    );
    const total = share.reduce((a, b) => a + b, 0);
    species.forEach((sp, i) => {
      if (sp.z === 0 || !(share[i] > 0.1 * total)) return;
      const cc = sol.c[sp.name];
      let lo = Infinity, hi = -Infinity;
      for (const v of cc) if (Number.isFinite(v)) [lo, hi] = [Math.min(lo, v), Math.max(hi, v)];
      const span = hi - lo;
      if (!(span > 0)) return;
      // At each region end, the two outermost cells, each with its inner neighbour.
      const ends = [];
      model.regions.forEach((reg, r) => {
        const a = grid.regionStart[r], b = grid.regionEnd[r];
        if (model.materials[reg.material].conductor || b - a < 3) return;
        ends.push([[a, a + 1], [a + 1, a + 2]], [[b - 1, b - 2], [b - 2, b - 3]]);
      });
      const field = (j) => (sol.phi[j + 1] - sol.phi[j]) / grid.segLength[j];
      let worst = 0, at = -1, hAt = 0;
      for (const cells of ends) {
        // The bulk fixes the flux, so a cell whose flux comes out too large by e takes too small
        // a share of the concentration drop; the missing drop shifts the whole profile, and the
        // current with it, by about e times the cell's share of the drop.
        let err = 0;
        for (const [s, k] of cells) {
          const ratio = Math.min(cc[s], cc[s + 1]) / Math.max(cc[s], cc[s + 1]);
          if (!(ratio < 0.5) || !(Math.abs(field(s)) > 2 * Math.abs(field(k))) || Math.abs(sol.phi[s + 1] - sol.phi[s]) < VT) continue;
          const e = ((1 + ratio) * Math.log(1 / ratio)) / (2 * (1 - ratio)) - 1;
          err += (e * Math.abs(cc[s + 1] - cc[s])) / span;
        }
        if (err > worst) [worst, at, hAt] = [err, cells[0][0], grid.segLength[cells[0][0]]];
      }
      if (worst > 0.003) {
        sol.warnings.push(
          `${sp.name}: a steep profile near x = ${grid.x[at].toExponential(3)} m on a ${hAt.toExponential(2)} m cell, where the field ` +
            `isn't uniform across the cell; the current may come out too large (an estimate, likely low: ${(100 * worst).toFixed(worst < 0.1 ? 1 : 0)}%). ` +
            'Refine the grid toward that end (a smaller hmin).',
        );
      }
    });
  }

  // Conservation bookkeeping for each species stretch.
  // Amount now vs the reference amount plus what came in through the contacts. The drift is
  // relative to the larger of the two; NaN where a reaction also makes or consumes it.
  sol.conservation = solver.stretches.map((st, k) => {
    const amount = solver.amount(st);
    const reference = solver.referenceAmounts[k];
    const intake = solver.boundaryIntake[k];
    const expected = reference + intake;
    const scale = Math.max(Math.abs(expected), Math.abs(amount), Math.abs(intake), 1e-300);
    return {
      species: species[st.species].name,
      regions: st.regions.slice(),
      spectator: st.spectator,
      connected: st.connected,
      reactive: st.reactive,
      amount,
      reference,
      intake,
      drift: st.reactive ? NaN : (amount - expected) / scale,
    };
  });

  return sol;
}
