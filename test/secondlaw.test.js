import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT, bernoulli, units } from '../src/index.js';

// The second law, end to end. In a steady state, the free energy brought in through the terminals,
// Σ N_in μ̄_out over every species at every terminal, equals the total dissipation inside: each
// flux down its own μ̄ drop, N·(μ̄_L − μ̄_R), and each reaction's rate times its affinity, r·A. Every
// one of those terms is non-negative. The identity is exact for the discrete equations (sum each
// box's balance times its μ̄), so it checks the bookkeeping of fluxes, readouts and reactions
// against an independent recomputation of every dissipation term from the solution.
//
// V·I is the electrical part of the input; terminals exchanging several species also trade
// chemical free energy through their offsets (a concentration cell at open circuit runs on that
// alone).

const RT = GAS_CONSTANT * 298.15;

function freeEnergyBalance(dev, sol) {
  const { model } = dev, { grid, species } = model, n = species.length;
  const mu = species.map((sp) => sol.mu[sp.name]);
  const VT = RT / FARADAY;
  // Input: through each contact and port, at the outside phase's levels.
  const outside = (link, V, i) => (species[i].z === 0 ? link.mu : species[i].z * FARADAY * (V + link.offset));
  let input = 0, electrical = 0;
  for (const side of ['left', 'right']) {
    const ct = model.contacts[side], r = sol.contacts[side], V = sol.terminals[side].V;
    species.forEach((sp, i) => {
      const link = ct.species[i];
      if (link.type === 'blocked') return;
      const nIn = side === 'left' ? r.flux[sp.name] : -r.flux[sp.name];
      input += nIn * outside(link, V, i);
    });
    electrical += sol.terminals[side].current * V;
  }
  model.ports.forEach((port, k) => {
    const V = sol.terminals[port.name].V;
    species.forEach((sp, i) => {
      if (port.species[i].type !== 'blocked') input += sol.ports[k].flux[sp.name] * outside(port.species[i], V, i);
    });
    electrical += sol.terminals[port.name].current * V;
  });

  // Dissipation, term by term.
  const terms = { transport: 0, conduction: 0, bulk: 0, interface: 0, links: 0 };
  let worst = 0; // the most negative single term
  const add = (key, v) => {
    terms[key] += v;
    worst = Math.min(worst, v);
  };
  for (let s = 0; s < grid.nNodes - 1; s++) {
    const r = grid.segRegion[s];
    if (r < 0) continue;
    const mat = model.materials[model.regions[r].material], h = grid.segLength[s];
    if (mat.conductor) {
      // Ohm's law: J = −g Δη, g = σRT/(z²F²h), dissipating g RT Δη².
      const i = mat.conductor.i, z = species[i].z;
      const g = (mat.conductor.sigma * RT) / (z * z * FARADAY * FARADAY * h), deta = (mu[i][s + 1] - mu[i][s]) / RT;
      add('conduction', g * RT * deta * deta);
      continue;
    }
    for (let i = 0; i < n; i++) {
      if (!mat.present[i] || mat.D[i] === 0) continue;
      // Scharfetter–Gummel: N = −(D/h) B(Δ) c_L expm1(Δη), Δ = zΔφ/V_T.
      const d = (species[i].z * (sol.phi[s + 1] - sol.phi[s])) / VT, deta = (mu[i][s + 1] - mu[i][s]) / RT;
      const N = -(mat.D[i] / h) * bernoulli(d) * sol.c[species[i].name][s] * Math.expm1(deta);
      add('transport', N * (mu[i][s] - mu[i][s + 1]));
    }
  }
  for (let g = 0; g < grid.nNodes; g++) {
    const m = model.regions[grid.nodeRegion[g]].material;
    for (const rx of model.reactions) {
      if (!(rx.kf[m] > 0)) continue;
      let P = rx.kf[m], a = rx.fixedA;
      for (const { i, nu } of rx.reactants) {
        P *= sol.c[species[i].name][g] ** nu;
        a += (nu * mu[i][g]) / RT;
      }
      for (const { i, nu } of rx.products) a -= (nu * mu[i][g]) / RT;
      add('bulk', grid.vol[g] * P * -Math.expm1(-a) * a * RT);
    }
  }
  model.interfaces.forEach((itf, f) => {
    const gL = grid.regionEnd[f], gR = gL + 1;
    species.forEach((sp, i) => {
      if (itf.links[i].type === 'conductance') add('interface', sol.interfaces[f].N[sp.name] * (mu[i][gL] - mu[i][gR]));
    });
    itf.reactions.forEach((rx, k) => {
      let a = rx.fixedA;
      for (const p of rx.part) a -= (p.nu * mu[p.i][p.side ? gR : gL]) / RT;
      add('interface', sol.interfaces[f].rates[k] * a * RT);
    });
  });
  // Links to outside phases that aren't in equilibrium (conductance, exchange): each flux in
  // drops from the outside level to the node's, at a contact or across a port's window.
  model.ports.forEach((port) => {
    const V = sol.terminals[port.name].V;
    species.forEach((sp, i) => {
      const link = port.species[i];
      if (link.type !== 'conductance' && link.type !== 'exchange') return;
      const muOut = outside(link, V, i);
      for (const g of port.nodes) {
        const k = link.type === 'conductance' ? (link.G * VT) / (sp.z * sp.z * FARADAY) : link.k;
        const nIn = grid.vol[g] * k * ((muOut - mu[i][g]) / RT);
        add('links', nIn * (muOut - mu[i][g]));
      }
    });
  });
  for (const side of ['left', 'right']) {
    const ct = model.contacts[side], V = sol.terminals[side].V, g = side === 'left' ? 0 : grid.nNodes - 1;
    species.forEach((sp, i) => {
      const link = ct.species[i];
      if (link.type !== 'conductance' && link.type !== 'exchange') return;
      const nIn = side === 'left' ? sol.contacts[side].flux[sp.name] : -sol.contacts[side].flux[sp.name];
      add('links', nIn * (outside(link, V, i) - mu[i][g]));
    });
  }
  const dissipation = terms.transport + terms.conduction + terms.bulk + terms.interface + terms.links;
  return { input, electrical, dissipation, terms, worst };
}

const check = (label, dev) => {
  const sol = dev.solve();
  assert.ok(sol.converged, label);
  const b = freeEnergyBalance(dev, sol);
  assert.ok(b.input > 0, `${label}: free energy must come in (${b.input})`);
  assert.ok(Math.abs(b.dissipation / b.input - 1) < 1e-10, `${label}: dissipation ${b.dissipation} vs input ${b.input}`);
  assert.ok(b.worst > -1e-12 * b.input, `${label}: a negative dissipation term (${b.worst})`);
  return b;
};

const Nc = units.perCm3(2.8e19), Nv = units.perCm3(1.04e19);
const ohmic = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium', 'h+': { type: 'equilibrium', offset: 0 } }, phi: 'bulk' });
const ions = [
  { name: 'Ag+', z: 1, cRef: 1000 },
  { name: 'NO3-', z: -1, cRef: 1000 },
];
const water = (epsr) => ({ epsr, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } });

test('second law: a forward-biased pn diode turns all its electrical input into heat', () => {
  const dev = new Device({
    species: [
      { name: 'e-', z: -1 },
      { name: 'h+', z: 1 },
    ],
    materials: { Si: { epsr: 11.7, species: { 'e-': { D: 36e-4, mu0: 0, cRef: Nc }, 'h+': { D: 12e-4, mu0: units.eV(1.12), cRef: Nv } } } },
    regions: [
      { material: 'Si', length: 1e-6, fixedCharge: units.perCm3(1e17) * FARADAY },
      { material: 'Si', length: 1e-6, fixedCharge: -units.perCm3(1e16) * FARADAY },
    ],
    contacts: { left: ohmic(0), right: ohmic(0.5) },
    bulkReactions: [{ nu: { 'e-': -1, 'h+': -1 }, kf: { Si: 1e13 } }], // diffusion length ~0.15 µm
    grid: { hmin: 1e-9, hmax: 20e-9 },
  });
  const b = check('pn diode', dev);
  assert.ok(b.terms.bulk > 0.01 * b.input, 'recombination takes a real share');
  // Electrons and holes both at the metal's level: all the input is electrical, V·I.
  assert.ok(Math.abs(b.electrical / b.input - 1) < 1e-9);
  assert.ok(b.terms.transport > 0 && b.terms.bulk > 0);
});

test('second law: a liquid junction at open circuit runs on chemical free energy alone', () => {
  const salt = [
    { name: 'Na+', z: 1, cRef: 1000 },
    { name: 'Cl-', z: -1, cRef: 1000 },
  ];
  const bath = (c) => ({ bath: { c: { 'Na+': c, 'Cl-': c }, reference: 'Cl-' } });
  const dev = new Device({
    species: salt,
    materials: { water: { epsr: 78.5, species: { 'Na+': { D: 1.33e-9, mu0: -261.9e3 }, 'Cl-': { D: 2.03e-9, mu0: -131.2e3 } } } },
    regions: [{ material: 'water', length: 10e-6 }],
    contacts: { left: bath(100), right: { ...bath(10), I: 0 } },
    grid: { hmin: 0.2e-9, hmax: 50e-9, ratio: 1.15 },
  });
  const b = check('liquid junction', dev);
  assert.ok(Math.abs(b.electrical) < 1e-12 * b.input, 'no electrical work at open circuit');
});

test('second law: silver electrodes and a bipolar plate, with resistance, reactions and an interface conductance', () => {
  const plating = (metal) => ({ [metal]: { 'e-': -1, Ag: 1 }, [metal === 'left' ? 'right' : 'left']: { 'Ag+': -1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.3 });
  const face = (metal) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating(metal)] });
  const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
  const salt = { 'NO3-': 10, 'Ag+': 10 };
  const dev = new Device({
    species: [...ions, { name: 'e-', z: -1 }],
    materials: { water: water(78.5), gel: water(78.5), Ag: { conductor: { species: 'e-', conductivity: 1e3 } } },
    regions: [
      { material: 'Ag', length: 1e-6 },
      { material: 'water', length: 5e-6, c0: salt },
      { material: 'gel', length: 5e-6, c0: salt },
      { material: 'Ag', length: 1e-6 },
      { material: 'water', length: 5e-6, c0: salt },
      { material: 'Ag', length: 1e-6 },
    ],
    interfaces: [face('left'), { dipole: 0, species: { 'Ag+': { type: 'conductance', G: 20 } } }, face('right'), face('left'), face('right')],
    contacts: { left: collector(0), right: collector(0.4) },
    grid: { hmin: 0.1e-9, hmax: 100e-9, ratio: 1.15 },
  });
  const b = check('silver cell', dev);
  assert.ok(Math.abs(b.electrical / b.input - 1) < 1e-9, 'only electrons cross the terminals');
  assert.ok(b.terms.transport > 0 && b.terms.conduction > 0 && b.terms.interface > 0);
});

test('second law: silver nitrate with a port behind a conductance and a contact through one', () => {
  const links = { terminal: 'Ag+', species: { 'Ag+': 'equilibrium', 'NO3-': 'blocked' }, phi: 'bulk' };
  const dev = new Device({
    species: ions,
    materials: { water: water(0) },
    regions: [{ name: 'cell', material: 'water', length: 20e-6, c0: { 'NO3-': 10 } }],
    contacts: { left: { V: 0, ...links }, right: { V: 0.05, ...links, species: { 'Ag+': { type: 'conductance', G: 50 }, 'NO3-': 'blocked' } } },
    ports: [{ name: 'p', region: 'cell', from: 8e-6, to: 12e-6, V: 0.02, terminal: 'Ag+', species: { 'Ag+': { type: 'conductance', G: 1e7 } } }],
    grid: { hmin: 10e-9, hmax: 0.5e-6 },
  });
  const b = check('port and conductance links', dev);
  assert.ok(b.terms.links > 0.01 * b.input, 'the links take a real share');
});
