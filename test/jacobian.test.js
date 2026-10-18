import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, units } from '../src/index.js';
import { build, layer, ohmic, bath, half, aqueous, metal, semiconductor } from '../src/kit.js';

// The Newton system against finite differences of its own residual: every row the solver writes
// (the kernels', the laws' in place of the rows they take, the pins held aside, the floating
// terminals' circuit rows) must have the derivatives of what its residual says. A wrong one
// still converges to the right answer, the residual being right, so no other test sees it: a
// node-by-node law kept a port reaction's ∂/∂V from the row it replaced, and Newton took an
// extra iteration (found by a cold review of src/lawRows.js, whose harness this is).
//
// Checked inside a real solve(), at its second steady assembly (the state on its way, not at a
// solution), and inside a real advance(), at its second time step's assembly, as Newton solves
// it there: at a strictly neutral node, its rows combined into charge continuity and neutrality,
// and φ̂ moving with every charged species' η (η_i = η'_i + z_i φ̂', see Solver._chargeRows).

const steadyAt = (s, dt) => s.steady && dt === Infinity;

// The full residual as Newton solves it: block rows, each pin's law in its row, then the
// floating terminals' circuit rows.
function residual(s, dt) {
  s.assemble(dt);
  const N = s.sys.size, f = new Float64Array(N + s.floating.length);
  for (let j = 0; j < N; j++) f[j] = s.res[j];
  for (const p of s.pins) f[s.rix[p.row]] = p.res;
  s.floating.forEach((k, a) => (f[N + a] = s.termRes[k]));
  return f;
}

// The Jacobian the solver's system stands for, dense: the blocks with each floating terminal's
// column (termB), a pin's law in place of its row (a flat stretch's pin moving its whole
// stretch), and the circuit rows.
function jacobian(s) {
  const { sys, nB } = s, N = sys.size, fl = s.floating, K = fl.length, R = s.rix;
  const { A, B, C, sizes, offA, offB, offC, offX } = sys;
  const J = Array.from({ length: N + K }, () => new Float64Array(N + K));
  for (let b = 0; b < nB; b++) {
    const m = sizes[b];
    for (let l = 0; l < m; l++) {
      const r = offX[b] + l, row = J[r];
      for (let k = 0; k < m; k++) row[offX[b] + k] += B[offB[b] + l * m + k];
      if (b > 0) for (let k = 0, mp = sizes[b - 1]; k < mp; k++) row[offX[b - 1] + k] += A[offA[b] + l * mp + k];
      if (b < nB - 1) for (let k = 0, mn = sizes[b + 1]; k < mn; k++) row[offX[b + 1] + k] += C[offC[b] + l * mn + k];
      fl.forEach((k, a) => (row[N + a] += s.termB[k][r]));
    }
  }
  for (const p of s.pins) {
    const r = R[p.row];
    J[r].fill(0);
    for (let j = 0; j < p.len; j++) J[r][p.idx[j]] += p.w[j];
    if (p.col) for (const o of p.col) if (o !== p.row && s.loc[o] >= 0) J[R[o]][r] -= 1;
  }
  fl.forEach((k, a) => {
    const row = J[N + a];
    for (let j = 0; j < N; j++) row[j] = s.termC[k][j];
    row[N + a] = s.termDI[k];
    fl.forEach((k2, b) => b !== a && (row[N + b] = s.termDX[k][k2]));
  });
  return J;
}

// Each column of the Jacobian against central differences of the residual, where they differ
// by more than 1e-5 relative and 1e-7 of the row's largest entry. A flat stretch's nodes all
// follow its first (the assembly snaps them there), so its first node's column is the sum of
// theirs and the others' are empty. In a time step's unknowns, a neutral node's φ̂ moves its
// charged species' η with it.
function differences(s, dt) {
  const N = s.sys.size, fl = s.floating, K = fl.length, M = s.M;
  const comp = new Int32Array(N).fill(-1);
  for (let o = 0; o < s.nB * M; o++) if (s.loc[o] >= 0) comp[s.rix[o]] = o;
  const u0 = Float64Array.from(s.u), l0 = Float64Array.from(s.uLo), v0 = Float64Array.from(s.termV);
  const restore = () => {
    s.u.set(u0);
    s.uLo.set(l0);
    s.termV.set(v0);
  };
  residual(s, dt);
  const J = jacobian(s), rowMax = J.map((row) => row.reduce((m, x) => Math.max(m, Math.abs(x)), 0));
  const together = new Map(), moves = new Map();
  for (const g of s.transformed ? s.chargeNodes : []) {
    const lb = s.blockOfNode[g] * M;
    if (s.loc[lb] < 0) continue;
    const shift = [[lb, 1]];
    for (let i = 0; i < s.n; i++) if (s.present[g * s.n + i] && s.z[i] !== 0) shift.push([lb + 1 + i, s.z[i]]);
    moves.set(s.rix[lb], shift);
  }
  for (const st of dt === Infinity && s.steady ? s.flatStretches : []) {
    const first = s.rix[s.blockOfNode[st.nodes[0]] * M + 1 + st.species], members = [];
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
      const o = s.blockOfNode[g] * M + 1 + st.species;
      if (s.loc[o] >= 0) members.push(s.rix[o]);
    }
    for (const r of members) together.set(r, r === first ? members : []);
  }
  const bad = [];
  for (let j = 0; j < N + K; j++) {
    const nudge = (h) => {
      if (j >= N) s.termV[fl[j - N]] += h;
      else for (const [o, w] of moves.get(j) ?? [[comp[j], 1]]) s.u[o] += w * h;
    };
    const h = j < N ? 1e-6 * Math.max(1, Math.abs(u0[comp[j]])) : 1e-6;
    nudge(h);
    const fp = residual(s, dt);
    restore();
    nudge(-h);
    const fm = residual(s, dt);
    restore();
    const cols = together.get(j) ?? [j];
    for (let i = 0; i < N + K; i++) {
      let code = 0;
      for (const c of cols) code += J[i][c];
      const fd = (fp[i] - fm[i]) / (2 * h), diff = Math.abs(code - fd), floor = 1e-7 * rowMax[i] + 1e-300;
      if (diff > floor && diff / (Math.abs(code) + Math.abs(fd) + floor) > 1e-5) bad.push(`row ${i} col ${j}: ${code.toExponential(4)} written, ${fd.toExponential(4)} by differences`);
    }
  }
  restore();
  residual(s, dt); // (the system as it was, for the solve to go on with)
  return bad;
}

// Run `go` on the device, checking the system at the nth assembly `when` picks.
function checkInside(def, go, when, nth = 2) {
  const d = new Device(def), s = d.solver, assemble = s.assemble;
  let seen = 0, bad = null, busy = false;
  s.assemble = function (dt) {
    assemble.call(this, dt);
    if (busy || bad || !when(this, dt) || ++seen < nth) return;
    busy = true;
    try {
      bad = differences(this, dt);
    } finally {
      busy = false;
    }
  };
  go(d);
  delete s.assemble;
  return bad;
}

// The devices: each law and kernel the steady system writes, from the other tests.
const film = { 'Na+': 100, 'K+': 10, 'Cl-': 110 };
const capCell = (end) =>
  build({
    T: 298.15,
    library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 80 })],
    stack: [end.contact ?? {}, layer('water', 1e-6, { name: 'film', c0: film }), bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-')],
    ...(end.port ? { ports: [end.port] } : {}),
    grid: { hmin: 5e-8, hmax: 5e-8 },
  });
const Kc = { 'K+': { type: 'conductance', G: 1e6, offset: 0 } };
const lattice = (def, species, cMax) => {
  for (const m of Object.values(def.materials)) if (m.species && !m.conductor) m.statistics = [{ type: 'lattice', species: species.filter((x) => m.species[x]), cMax }];
  return def;
};
const adsorbing = () => {
  const l = aqueous(['Na+', 'OH-', 'Cl-'], { epsr: 0 });
  l.species.push({ name: 'e-', z: -1 }, { name: 'S0', z: 0 });
  l.materials.water.species.S0 = { D: 0, mu0: 0, cRef: 1 };
  return build({
    T: 298.15,
    library: [l],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: { 'Na+': 100, 'OH-': 1, 'Cl-': 99, S0: 1 } }), bath({ 'Na+': 100, 'OH-': 1, 'Cl-': 99 }, 'Cl-')],
    ports: [{ name: 'metal', region: 'film', from: 0, to: 4e-6, I: 0, terminal: 'e-', area: 1e4, surface: { Cads: { mu0: 0, capacity: 1e-5, theta0: 0.2 } }, reactions: [{ equation: 'S0 + OH- = Cads + e-', k0: 1e-4, alpha: 0.5 }] }],
    grid: { hmin: 1e-6, hmax: 1e-6 },
  });
};
const devices = {
  // An immobile species adsorbing on a floating electrode: a node-by-node law with a surface, a
  // circuit-row law (and the same with crowding, its ∂c through the statistics).
  'adsorbing on a floating electrode': adsorbing,
  'adsorbing, crowded': () => lattice(adsorbing(), ['Na+', 'Cl-', 'S0'], 2000),
  // A contact at I = 0 passing K⁺, its charge in the law: capacitive, bulk; a membrane port's.
  'a capacitive contact passing K⁺': () => capCell({ contact: { I: 0, zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: Kc } }),
  'a capacitive contact, crowded': () => lattice(capCell({ contact: { I: 0, zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: Kc } }), ['Na+', 'K+', 'Cl-'], 3000),
  'a bulk contact passing K⁺': () => capCell({ contact: { I: 0, phi: 'bulk', species: Kc } }),
  'a membrane capacitor passing K⁺': () => capCell({ port: { name: 'mem', region: 'film', I: 0, area: 1e6, capacitance: { C: 0.1, zeroCharge: -0.05 }, species: { 'K+': { type: 'conductance', G: 1e12, offset: 0 } } } }),
  // Face reactions at a floating metal, two couples (test/steady.test.js's pool).
  'a floating Pt with two couples': () => {
    const silver = half('Ag+ + e- = Ag(s)', { 'Ag(s)': 0 }), iron = half('Fe3+ + e- = Fe2+');
    const bathC = { 'K+': 5000, 'Fe3+': 10, 'Fe2+': 10, 'NO3-': 5050 };
    return build({
      library: [aqueous(['K+', 'NO3-', 'Fe3+', 'Fe2+', 'Ag+'], { epsr: 0 }), metal('Pt')],
      stack: [ohmic({ I: 0 }, ['e-']), layer('Pt', 1e-6), { reactions: [{ ...iron, k0: 1e-3 }, { ...silver, k0: 1e-3 }] }, layer('water', 20e-6, { c0: { ...bathC, 'Ag+': 1, 'NO3-': 5051 } }), bath(bathC, 'NO3-')],
      grid: { hmin: 10e-9, hmax: 0.5e-6 },
    });
  },
  // A metal floating between capacitive faces (test/metal.test.js).
  'a floating gate': () => ({
    species: [{ name: 'e-', z: -1 }],
    materials: { ox: { epsr: 3.9, species: {} }, Au: { conductor: { species: 'e-', conductivity: 4e7 } } },
    regions: [{ material: 'ox', length: 5e-9 }, { material: 'Au', length: 20e-9 }, { material: 'ox', length: 10e-9 }],
    interfaces: [{ phi: { type: 'capacitive', C: 10 }, zeroCharge: 0.3 }, { phi: { type: 'capacitive', C: 10 }, zeroCharge: 0.3 }],
    contacts: { left: { V: 0, phi: { type: 'capacitive', C: 1 }, zeroCharge: 0 }, right: { V: 1, phi: { type: 'capacitive', C: 1 }, zeroCharge: 0.2 } },
    grid: { minCells: 4 },
  }),
  // A MOS floating gate over silicon: a semiconductor's flat stretches and a metal's charge.
  'a MOS floating gate': () =>
    build({
      T: 300,
      library: [semiconductor('Si'), metal('Al'), { species: [], materials: { SiO2: { epsr: 3.9, species: {} } } }],
      stack: [ohmic({ I: 0 }, ['e-']), layer('Al', units.nm(20)), { phi: { type: 'capacitive', C: 10 }, zeroCharge: -0.95 }, layer('SiO2', units.nm(10)), { dipole: 0 }, layer('Si', units.um(1), { name: 'Si', acceptors: units.perCm3(1e17) }), ohmic(0)],
      bulkReactions: [{ equation: 'e- + h+ = 0', kf: { Si: 1e8 } }],
      grid: { hmin: units.nm(0.5), hmax: units.nm(100), ratio: 1.3 },
    }),
  // A bipolar silver plate between two plating faces, held at a bias.
  'a bipolar plate': () => {
    const water = { epsr: 78.5, species: { 'Ag+': { D: 1.65e-9, mu0: 77.1e3 }, 'NO3-': { D: 1.9e-9, mu0: -111.3e3 } } };
    const plating = (m) => ({ [m]: { 'e-': -1, Ag: 1 }, [m === 'left' ? 'right' : 'left']: { 'Ag+': -1 }, fixed: { Ag: 0 }, k0: 1e-3, alpha: 0.5 });
    const face = (m) => ({ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0.1, reactions: [plating(m)] });
    const salt = { 'NO3-': 10, 'Ag+': 10 }, Ag = { material: 'Ag', length: 1e-6 };
    const collector = (V) => ({ V, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' });
    return {
      species: [{ name: 'Ag+', z: 1, cRef: 1000 }, { name: 'NO3-', z: -1, cRef: 1000 }, { name: 'e-', z: -1 }],
      materials: { water, Ag: { conductor: { species: 'e-', conductivity: 6e7 } } },
      regions: [Ag, { material: 'water', length: 10e-6, c0: salt }, { material: 'Ag', length: 2e-6 }, { material: 'water', length: 10e-6, c0: salt }, Ag],
      interfaces: [face('left'), face('right'), face('left'), face('right')],
      contacts: { left: collector(0), right: collector(0.2) },
      grid: { hmin: 2e-9, hmax: 500e-9, ratio: 1.3 },
    };
  },
};

for (const [name, def] of Object.entries(devices)) {
  test(`the Newton system's derivatives are its residual's: ${name}`, () => {
    const steady = checkInside(def(), (d) => d.solve(), steadyAt);
    assert.ok(steady, 'no steady assembly was checked');
    assert.deepEqual(steady.slice(0, 5), [], 'steady');
    const step = checkInside(def(), (d) => d.advance(1e-3, { tol: 1e-6 }), (s, dt) => Number.isFinite(dt));
    assert.ok(step, 'no time step was checked');
    assert.deepEqual(step.slice(0, 5), [], 'a time step');
  });
}
