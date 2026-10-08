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
    T: model.T,
    time: solver.time,
    steady: solver.atSteady,
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
  // Each species' flux toward +x across each segment (between nodes g and g + 1; a face's across
  // its pair of nodes), and its diffusivity at each node (NaN where it's absent or in a metal).
  sol.flux = {};
  sol.D = {};
  solver.segmentFluxes().forEach((J, i) => (sol.flux[species[i].name] = J));
  for (let i = 0; i < n; i++) {
    sol.D[species[i].name] = Float64Array.from(grid.nodeRegion, (r, g) => {
      const mat = model.materials[model.regions[r].material];
      return mat.conductor || !solver.present[g * n + i] ? NaN : mat.D[i];
    });
  }
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
    // How each species crosses it (a link type), for those that do.
    sol.contacts[side].links = Object.fromEntries(ct.species.flatMap((l, i) => (l.type === 'blocked' ? [] : [[species[i].name, l.type]])));
    if (ct.phi.type === 'capacitive' || ct.phi.type === 'pinned') {
      // Charge on the gate (or metal) plate, per area if the device is planar: +D at the left, −D
      // at the right.
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
    const out = { name: port.name, V: solver.termV[2 + k], flux, current };
    if (port.capacitance) {
      // A capacitance's charging current, and its charge (the port's side), like a gate's; and
      // its charge per area of electrode at each node, σ = C (V − zeroCharge − φ) (none at an
      // end whose φ its contact sets: the contact holds that charge).
      const dt = solver.lastDt;
      out.current += Number.isFinite(dt) ? (solver.portQ[k] - solver.portQOld[k]) / dt : 0;
      out.charge = solver.portQ[k];
      const { C: Cs, zeroCharge } = port.capacitance, last = model.grid.nNodes - 1;
      const sets = (side) => ['pinned', 'bulk'].includes(model.contacts[side].phi.type);
      out.sigma = Float64Array.from(port.nodes, (g) => ((g === 0 && sets('left')) || (g === last && sets('right')) ? 0 : Cs * (solver.termV[2 + k] - zeroCharge - phi[g])));
    }
    if (port.gates.length > 0) {
      // A membrane's gates at each node, and the voltage they follow, φ − (V − zeroCharge).
      out.gates = Object.fromEntries(port.gates.map((gate, q) => [gate.name, Float64Array.from(port.nodes, (g) => u[solver.blockOfNode[g] * M + 1 + n + port.surface.length + q])]));
      out.Vm = Float64Array.from(port.nodes, (g) => solver._portVoltage(k, g));
    }
    if (port.area !== null) {
      out.x = Float64Array.from(port.nodes, (g) => model.grid.x[g]);
      out.area = Float64Array.from(solver.portArea[k]);
    }
    if (port.reactions.length > 0) {
      // An electrode spread through the window: each reaction's rate at its nodes.
      out.rates = solver.portRates[k].map((r) => Float64Array.from(r)); // mol/(m²·s), forward
      // its surface: each species' coverage θ at each node of the window
      if (port.surface.length > 0) {
        out.coverage = Object.fromEntries(port.surface.map((sp, q) => [sp.name, Float64Array.from(port.nodes, (g) => solver.th[g * solver.nSurf + q])]));
        out.bare = Float64Array.from(port.nodes, (g) => solver.th0[g]); // θ₀, to full precision however covered
      }
    }
    return out;
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
    const out = { left: model.regions[f].name, right: model.regions[f + 1].name, dipole: itf.dipole, sheetCharge: itf.sheetCharge, D: u[b * M], N, rates };
    // How each species crosses it (a link type), for those that do; and for those crossing by
    // permeability, the one-way fluxes [toward +x, toward −x], per area.
    out.links = Object.fromEntries(itf.links.flatMap((l, i) => (l.type === 'blocked' ? [] : [[species[i].name, l.type]])));
    const oneWay = solver.faceOneWay(f);
    if (oneWay.some((w) => w)) out.oneWay = Object.fromEntries(oneWay.flatMap((w, i) => (w ? [[species[i].name, w]] : [])));
    if (itf.gates.length > 0) {
      // Each gate's fraction open, and the voltage across the face that drives them.
      out.gates = Object.fromEntries(itf.gates.map((gate, q) => [gate.name, u[solver._gateSlot(f, q)]]));
      out.V = solver._faceVoltage(f);
    }
    return out;
  });

  // Bulk reactions: each one's forward rate at every node (mol/(m³·s), NaN where it doesn't run),
  // and integrated over each region and the device (mol/(m²·s)), by the nodes' control volumes,
  // just as the balances count it.
  sol.bulkReactions = solver.rxs.map((rx) => {
    const rate = new Float64Array(nNodes).fill(NaN), forward = new Float64Array(nNodes).fill(NaN), regions = model.regions.map(() => 0);
    let total = 0;
    for (let g = 0; g < nNodes; g++) {
      const m = solver.nodeMaterial[g];
      if (!(rx.kf[m] > 0)) continue;
      rate[g] = solver.bulkRate(rx, g);
      forward[g] = solver.bulkForward(rx, g);
      const amount = grid.vol[g] * rate[g];
      regions[grid.nodeRegion[g]] += amount;
      total += amount;
    }
    // Its two one-way rates are forward and forward − rate; what each makes of each species, per
    // forward reaction: ν, negative for what it consumes.
    const nu = {};
    rx.sp.forEach((i, p) => (nu[species[i].name] = (nu[species[i].name] ?? 0) - rx.nu[p]));
    return { rate, forward, nu, regions, total };
  });

  // Total charge in the device (space charge plus sheet charges), per area if it's planar.
  let q = 0;
  for (let g = 0; g < nNodes; g++) {
    let rho = solver.rhoFixed[g];
    for (let i = 0; i < n; i++) rho += F * species[i].z * c[g * n + i];
    q += grid.vol[g] * rho;
  }
  interfaces.forEach((itf, f) => (q += itf.sheetCharge * grid.area[grid.regionEnd[f]]));
  sol.charge = q;

  // Resolution warnings: where the model resolves a double layer (dipole or capacitive faces,
  // gate contacts), check the local Debye length against the adjacent cell.
  sol.warnings = [...(model.warnings ?? []), ...(result.stopped ? [result.stopped] : [])];
  // A surface covered almost completely: its coverages' η barely move anything below ~1e-10 bare,
  // and Newton loses them (a transient can stall there).
  sol.ports.forEach((p) => {
    if (!p.bare) return;
    const min = Math.min(...p.bare);
    if (min < 1e-10) {
      sol.warnings.push(
        `ports ${p.name}: its surface is covered to a bare fraction of ${min.toExponential(1)} at x = ${p.x[p.bare.indexOf(min)].toExponential(3)} m. ` +
          'Below ~1e-10 the coverage barely affects anything and the solver struggles with it; a passive metal blocks to 1e-3–1e-6, so make the film less stable (raise its μ°).',
      );
    }
  });
  const cond = solver.conditioning, un = solver.unreached, ch = solver.charging;
  if (result.converged === false && ch) {
    const A = model.geometry.type === 'planar' ? 'A/m²' : 'A';
    sol.warnings.push(
      `${ch.terminal}: driven at ${ch.current.toPrecision(3)} ${A}, it has no steady state: what it feeds (${ch.what.join(', ')}) is closed but for it, ` +
        'so it only stores what comes in (a host filling, a capacitor charging). advance() in time instead; at no current, a steady solve keeps the charge it holds.',
    );
  } else if (result.converged === false && un) {
    const [[V0, I0], [V1, I1]] = [un.low, un.high], g = (x) => x.toPrecision(3), A = model.geometry.type === 'planar' ? 'A/m²' : 'A';
    // The opposite sign within reach: likely a sign slip (a drive is the current into the device).
    const flipped = Math.min(I0, I1) <= -un.target && -un.target <= Math.max(I0, I1);
    sol.warnings.push(
      `${un.terminal}: no steady state passes the driven current ${g(un.target)} ${A}; held from ${g(V0)} V to ${g(V1)} V ` +
        `it passed ${g(I0)} to ${g(I1)} ${A} (a limiting current, or kinetics too slow). Drive less, or hold a voltage.` +
        (flipped ? ` (It passes ${g(-un.target)} ${A}: if that was meant, a terminal's current is into the device.)` : '') +
        // Nothing passed at any voltage: the device only stores what comes in (a host filling, a
        // capacitor charging), so it has no steady state at a current.
        (Math.max(Math.abs(I0), Math.abs(I1)) < 1e-6 * Math.abs(un.target)
          ? ' It passes almost nothing steadily at any voltage, so it may only store what comes in (a host filling, a capacitor charging), which has no steady state at a current: advance() in time instead.'
          : ''),
    );
  } else if (result.converged === false && cond && cond.digits > 12) {
    const digits = Number.isFinite(cond.digits) ? `lost ${cond.digits.toFixed(0)} of its ~16 digits` : 'was exactly singular';
    // An ion there so scarce beside the others (a minority swept out of a junction by a large
    // bias, excluded at a sharp neutral face) that its level can't be held in double precision.
    let scarce = null;
    for (const g of cond.nodes ?? []) {
      let all = 0;
      for (let i = 0; i < model.species.length; i++) if (model.species[i].z !== 0 && solver.c[g * model.species.length + i] > 0) all += Math.abs(model.species[i].z) * solver.c[g * model.species.length + i];
      for (let i = 0; i < model.species.length; i++) {
        const c = solver.c[g * model.species.length + i], r = c / all;
        if (model.species[i].z !== 0 && c > 0 && r < 1e-10 && !(r >= scarce?.r)) scarce = { r, name: model.species[i].name };
      }
    }
    sol.warnings.push(
      `the steady system ${digits} near x = ${cond.x.toExponential(3)} m (${cond.where}): ` +
        (scarce?.r < 1e-24
          ? `there ${scarce.name} is ${scarce.r.toExponential(0)} of the ions around it, too scarce for its level to be held in double precision ` +
            '(a minority swept out of a junction by a large bias, or excluded at a sharp neutral face). Drive less; or, at a face, ' +
            'resolve its double layer (ε > 0) rather than make it a sharp neutral step.'
          : 'part of the device is held only weakly, e.g. a floating region coupled to the rest through tiny conductances or rates, or a stiff ' +
            'chain whose level nothing pins. Strengthen that coupling, or anchor the region (a port, a contact).' +
            (scarce
              ? ` Or a species there too scarce to hold: ${scarce.name} is ${scarce.r.toExponential(0)} of the ions around it (a layer depleted past a limiting current, say), which a smaller drive avoids.`
              : '')),
    );
  }
  if (result.converged === false && !result.stopped && !un && !ch && !(cond && cond.digits > 12)) {
    sol.warnings.push(
      'the solve did not converge, and nothing specific showed why. ' +
        'Check the definition with describe() (from driftlet/kit); approach this state in smaller steps (a sweep, or advance() in time); or refine the grid.',
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
  // Beside a strictly neutral (ε = 0) material, a face's charge is held in the edge node's box,
  // which then acts as a diffuse layer as wide as the box: a capacitance F²Σz²c·w/RT in series
  // with the face's (a Helmholtz C, or for a pinned face whatever holds the charge on the other
  // side: a conductor's sheet, or a resolved diffuse layer, ε/λ_D), which shrinks as the grid is
  // refined. Flag it where it lowers that capacitance by over 2%. (Where the box is the smaller
  // part, as an ε = 0 metal's against a semiconductor's depletion, it's a fine stand-in for the
  // metal's own screening.)
  const boxed = (g, h, phi, other, where) => {
    const mat = model.materials[model.regions[grid.nodeRegion[g]].material];
    if (mat.conductor || mat.epsr !== 0 || solver.phiUndefined?.[g]) return;
    const s2 = solver.screening(g), Cbox = (F * F * s2 * h) / (2 * RT);
    if (!(s2 > 0)) return;
    let C = phi.C;
    if (phi.type === 'pinned') {
      const om = other < 0 ? null : model.materials[model.regions[grid.nodeRegion[other]].material];
      C = !om || om.conductor ? Infinity : (om.epsr * EPS0) / debye(other);
    }
    if (!(C / (C + Cbox) > 0.02)) return;
    const lowers = Number.isFinite(C) ? `lowering the double-layer capacitance by ${((100 * C) / (C + Cbox)).toFixed(0)}%` : 'which alone sets its charge';
    sol.warnings.push(
      `${where}: the face's charge sits in the strictly neutral cell beside it (${h.toExponential(2)} m), which acts as a diffuse layer that wide ` +
        `(${Cbox.toPrecision(2)} F/m² in series), ${lowers}, so it depends on the mesh. ` +
        (Number.isFinite(C)
          ? `Coarsen the grid there: a cell of ${((2 * 49 * C * RT) / (F * F * s2)).toExponential(1)} m or more keeps it within 2%.`
          : "Give the face a Helmholtz capacitance (phi: { type: 'capacitive', C }), or make it neutral."),
    );
  };
  interfaces.forEach((itf, f) => {
    if (itf.phi.type === 'neutral') return;
    const gL = grid.regionEnd[f], gR = grid.regionStart[f + 1];
    check(gL, grid.segLength[gL - 1], `interfaces[${f}] (left side)`);
    check(gR, grid.segLength[gR], `interfaces[${f}] (right side)`);
    boxed(gL, grid.segLength[gL - 1], itf.phi, gR, `interfaces[${f}] (left side)`);
    boxed(gR, grid.segLength[gR], itf.phi, gL, `interfaces[${f}] (right side)`);
  });
  const resolves = (ct) => ct.phi.type === 'capacitive' || ct.phi.type === 'pinned';
  if (resolves(contacts.left)) {
    check(0, grid.segLength[0], 'contacts.left');
    boxed(0, grid.segLength[0], contacts.left.phi, -1, 'contacts.left');
  }
  if (resolves(contacts.right)) {
    check(nNodes - 1, grid.segLength[nNodes - 2], 'contacts.right');
    boxed(nNodes - 1, grid.segLength[nNodes - 2], contacts.right.phi, -1, 'contacts.right');
  }
  // Electrons or holes above their band's effective density of states (the material's cRef)
  // under Boltzmann statistics: degenerate, where those statistics no longer hold (a surface
  // driven past what a reaction can take away piles holes up to 1e9 mol/m³ and still converges).
  model.regions.forEach((reg, r) => {
    const mat = model.materials[reg.material];
    if (mat.conductor) return;
    model.species.forEach((sp, i) => {
      if ((sp.name !== 'e-' && sp.name !== 'h+') || !mat.present[i] || mat.modelOf?.[i] >= 0) return;
      let top = 0, at = NaN;
      for (let g = grid.regionStart[r]; g <= grid.regionEnd[r]; g++) {
        const c = solver.c[g * model.species.length + i];
        if (c > top) [top, at] = [c, grid.x[g]];
      }
      if (top > mat.cRef[i]) {
        sol.warnings.push(
          `${reg.name}: ${sp.name} reaches ${top.toExponential(2)} mol/m³ at x = ${at.toExponential(3)} m, above its cRef (${mat.cRef[i].toExponential(2)}, the band's effective density of states), where Boltzmann statistics no longer hold: ` +
            'give the material Fermi–Dirac statistics (docs/statistics.md), or ask whether this state is physical (a surface driven past what can carry its carriers away).',
        );
      }
    });
  });
  // A contact that holds its end neutral ('bulk': a bath, an ohmic contact) cuts off any space
  // charge reaching it: a double layer longer than its region. The charge it cuts off is about
  // ρ λ_D at the node beside it; against the region's own space charge (half Σ|ρ| over it, the
  // two signs of a double layer), warn past 1e-3.
  for (const side of ['left', 'right']) {
    if (contacts[side].phi.type !== 'bulk') continue;
    const r = side === 'left' ? 0 : model.regions.length - 1, a = grid.regionStart[r], b = grid.regionEnd[r];
    const mat = model.materials[model.regions[r].material];
    if (!(mat.epsr > 0) || mat.conductor || b - a < 3) continue;
    const n = model.species.length, fixed = model.regions[r].fixedCharge ?? 0;
    const rho = (g) => {
      let q = fixed;
      for (let i = 0; i < n; i++) if (solver.c[g * n + i] > 0) q += F * model.species[i].z * solver.c[g * n + i];
      return q;
    };
    // (and beside the contact the material is itself out of neutrality, by more than round-off:
    // a quasi-neutral layer split off a junction's holds almost no space charge, so the ratio
    // below alone would read its round-off as a cut-off layer)
    const gross = (g) => {
      let q = Math.abs(fixed);
      for (let i = 0; i < n; i++) if (solver.c[g * n + i] > 0) q += F * Math.abs(model.species[i].z) * solver.c[g * n + i];
      return q;
    };
    let total = 0;
    for (let g = a; g <= b; g++) total += (Math.abs(rho(g)) * grid.vol[g]) / 2;
    const g1 = side === 'left' ? a + 1 : b - 1, lam = debye(g1), w = (grid.x[g1 + 1] - grid.x[g1 - 1]) / 2;
    const cut = (Math.abs(rho(g1)) * lam * grid.vol[g1]) / w;
    if (Number.isFinite(cut) && total > 0 && cut > 1e-3 * total && Math.abs(rho(g1)) > 1e-3 * gross(g1)) {
      sol.warnings.push(
        `contacts.${side}: the space charge in ${model.regions[r].name} reaches the contact, which holds its end neutral and so cuts it off (a double layer or depletion longer than the region): ` +
          `the charge beside the contact is ${(cut / total).toPrecision(2)} of the space charge in the region, and the error in that charge can be several times more. ` +
          'Lengthen the region (a crowded double layer grows as the root of the voltage).',
      );
    }
  }

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
