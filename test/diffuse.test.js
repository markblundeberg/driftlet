import { test } from 'node:test';
import assert from 'node:assert/strict';
import { AVOGADRO, Device, DeviceError, EPS0, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer, ohmic, bath, half, metal, aqueous, describe } from '../src/kit.js';
import { diffuseLayer, diffuseScratch, diffuseUnits } from '../src/diffuse.js';

// A sub-grid Gouy–Chapman layer at a face beside a strictly neutral electrolyte (src/diffuse.js,
// interfaces[f].diffuse). Its answer sheet was written first: the Poisson–Boltzmann problem
// solved in x (scipy's solve_bvp, within 2e-9), mpmath quadrature of the first integral (the
// table below) and the z:z closed forms (within 1e-16) agree; the kernel is held to the table.

// [name, c (mol/m³), z, s = Fψ/RT, q (C/m²), Γ (mol/m²)], εr 78.4 at 298.15 K, by the answer
// sheet's constants (F 96485.33212, R 8.314462618, ε0 8.8541878128e-12).
const TABLE = [
  ['KCl 0.1 M', [100, 100], [1, -1], -20, 408.6266729451325, [0.004234924762003149, -1.9226528674435084e-7]],
  ['KCl 0.1 M', [100, 100], [1, -1], -5, 0.22448221660722878, [2.1501030223978235e-6, -1.76491203634766e-7]],
  ['KCl 0.1 M', [100, 100], [1, -1], -1, 0.019334326568207803, [1.2473224396343482e-7, -7.565393021857927e-8]],
  ['KCl 0.1 M', [100, 100], [1, -1], 1, -0.019334326568207803, [-7.565393021857927e-8, 1.2473224396343482e-7]],
  ['KCl 0.1 M', [100, 100], [1, -1], 5, -0.22448221660722878, [-1.76491203634766e-7, 2.1501030223978235e-6]],
  ['KCl 0.1 M', [100, 100], [1, -1], 20, -408.6266729451325, [-1.9226528674435084e-7, 0.004234924762003149]],
  ['MgSO4 10 mM', [10, 10], [2, -2], 4, -0.32019467646172134, [-2.9844374025794124e-8, 1.6294476107055805e-6]],
  ['ferri/ferro in KCl', [5, 5, 1035, 1000], [-3, -4, 1, -1], -2, 0.1405280480067089, [-3.3937660898466594e-9, -3.8138986839007775e-9, 1.0565352496651903e-6, -3.7449844118931936e-7]],
  ['ferri/ferro in KCl', [5, 5, 1035, 1000], [-3, -4, 1, -1], 2, -0.27769282610550294, [7.673074596303917e-8, 3.724453797555037e-7, -3.4043536217789106e-7, 8.176742123766643e-7]],
  ['ferri/ferro in KCl', [5, 5, 1035, 1000], [-3, -4, 1, -1], 6, -675.987901377099, [8.535055497436696e-6, 0.0017447071139870293, -3.73869860675932e-7, 1.313402182261649e-6]],
  ['CaCl2 10 mM', [10, 20], [2, -1], -3, 0.11740822809594156, [5.876477257381608e-7, -4.155508441569852e-8]],
  ['CaCl2 10 mM', [10, 20], [2, -1], 3, -0.035768351712957586, [-4.29445290988018e-8, 2.848237841282335e-7]],
];
const sheetUnits = diffuseUnits(78.4 * 8.8541878128e-12, 8.314462618 * 298.15, 96485.33212);

test("a diffuse layer's charge and surface excesses match the Poisson–Boltzmann answer sheet", () => {
  for (const [name, c, z, s, q, gam] of TABLE) {
    const n = c.length, out = diffuseLayer(Float64Array.from(c), z, n, s, sheetUnits, diffuseScratch(n));
    assert.ok(Math.abs(out.q - q) <= 1e-13 * Math.abs(q), `${name} at s = ${s}: q ${out.q} for ${q}`);
    gam.forEach((g, i) => assert.ok(Math.abs(out.gam[i] - g) <= 1e-13 * Math.abs(g), `${name} at s = ${s}: Γ[${i}] ${out.gam[i]} for ${g}`));
    // Gauss's law across the layer: F Σ z Γ = q
    const sum = 96485.33212 * z.reduce((a, zi, i) => a + zi * out.gam[i], 0);
    assert.ok(Math.abs(sum - q) <= 1e-14 * Math.abs(q), `${name} at s = ${s}: F Σ z Γ ${sum} for q ${q}`);
  }
});

test("a diffuse layer's derivatives are its functions' own", () => {
  for (const [name, c, z, s] of TABLE.filter((r) => Math.abs(r[3]) <= 5)) {
    const n = c.length, at = (cc, ss) => diffuseLayer(Float64Array.from(cc), z, n, ss, sheetUnits, diffuseScratch(n));
    const o = at(c, s), h = 1e-5, p = at(c, s + h), m = at(c, s - h);
    const near = (fd, a, what) => assert.ok(Math.abs(fd - a) <= 1e-7 * Math.abs(a) + 1e-30, `${name} at s = ${s}: ${what} ${a} written, ${fd} by differences`);
    near((p.q - m.q) / (2 * h), o.dqds, '∂q/∂s');
    for (let i = 0; i < n; i++) near((p.gam[i] - m.gam[i]) / (2 * h), o.dgds[i], `∂Γ[${i}]/∂s`);
    for (let j = 0; j < n; j++) {
      const up = [...c], dn = [...c], dc = 1e-5 * c[j];
      up[j] += dc;
      dn[j] -= dc;
      const a = at(up, s), b = at(dn, s);
      near((a.q - b.q) / (2 * dc), o.dqdc[j], `∂q/∂c[${j}]`);
      // (the co-ions' excesses are small beside the counter-ions': their derivatives to the row's scale)
      const row = Math.max(...Array.from({ length: n }, (_, i) => Math.abs(o.dgdc[i * n + j])));
      for (let i = 0; i < n; i++) assert.ok(Math.abs((a.gam[i] - b.gam[i]) / (2 * dc) - o.dgdc[i * n + j]) <= 1e-7 * row, `${name} at s = ${s}: ∂Γ[${i}]/∂c[${j}]`);
    }
  }
});

// Crowded (Bikerman): the ions on a lattice of cMax sites, as the bulk's lattice statistics have
// them, c_i(s) = c_i e^{−z_i s}/D(s), D = 1 + Σ θ_j (e^{−z_j s} − 1), θ_j = c_j/cMax. Its answer
// sheet (bikerman.py): the BVP in x (within 5e-10), mpmath quadrature of G = cMax ln D + s Σ c z
// (below), and Kilic–Bazant–Ajdari's closed form for q of a z:z salt (to round-off).
// [name, c, z, on the lattice, cMax, s, q, Γ]
const CROWDED = [
  ['KCl 0.1 M, 0.71 nm sites', [100, 100], [1, -1], [1, 1], 4639.530685488421, -20, 0.5080158268921463, [4.996315542685478e-6, -2.6889747705435907e-7]],
  ['KCl 0.1 M, 0.71 nm sites', [100, 100], [1, -1], [1, 1], 4639.530685488421, 1, -0.019222580254097645, [-7.6331925817044e-8, 1.2289607943239935e-7]],
  ['KCl 0.1 M, 0.71 nm sites', [100, 100], [1, -1], [1, 1], 4639.530685488421, 40, -0.7598886436245672, [-3.2516346564739623e-7, 7.550527345879109e-6]],
  ['CuSO4 0.1 M, 0.71 nm', [100, 100], [2, -2], [1, 1], 4639.530685488421, -8, 0.4406926642529582, [2.1567997517142786e-6, -1.2692905229120115e-7]],
  ['ferri/ferro in KCl, 0.71 nm', [5, 5, 1035, 1000], [-3, -4, 1, -1], [1, 1, 1, 1], 4639.530685488421, -6, 0.2682340850523549, [-5.485619530516571e-9, -5.8955162932153116e-9, 1.9321940013056307e-6, -8.078174635086997e-7]],
  ['ferri/ferro in KCl, 0.71 nm', [5, 5, 1035, 1000], [-3, -4, 1, -1], [1, 1, 1, 1], 4639.530685488421, 10, -0.7277341684393017, [7.59152333056317e-8, 1.5916334031666858e-6, -7.422436661649903e-7, 2.0591018129076168e-7]],
  ['KCl 0.5 M + 1 M neutral, 0.6 nm', [500, 500, 1000], [1, -1, 0], [1, 1, 1], 7687.680866545586, -10, 0.4385127922421983, [3.964424972638353e-6, -5.804398540361775e-7, -5.060027812526466e-7]],
];

test("a crowded diffuse layer's charge and excesses match the Bikerman answer sheet, its neutral species displaced", () => {
  for (const [name, c, z, on, cMax, s, q, gam] of CROWDED) {
    const n = c.length, lat = { cMax, on: Uint8Array.from(on) };
    const at = (cc, ss) => diffuseLayer(Float64Array.from(cc), z, n, ss, sheetUnits, diffuseScratch(n), lat), out = at(c, s);
    assert.ok(Math.abs(out.q - q) <= 1e-13 * Math.abs(q), `${name} at s = ${s}: q ${out.q} for ${q}`);
    gam.forEach((g, i) => assert.ok(Math.abs(out.gam[i] - g) <= 1e-13 * Math.abs(g), `${name} at s = ${s}: Γ[${i}] ${out.gam[i]} for ${g}`));
    // its derivatives, by differences
    const h = 1e-5 * Math.max(1, Math.abs(s)), p = at(c, s + h), m = at(c, s - h);
    assert.ok(Math.abs((p.q - m.q) / (2 * h) - out.dqds) < 1e-7 * Math.abs(out.dqds), `${name} at s = ${s}: ∂q/∂s`);
    assert.ok(Math.abs((p.D - m.D) / (2 * h) - out.dDds) < 1e-6 * Math.abs(out.dDds), `${name} at s = ${s}: ∂D/∂s`);
    for (let j = 0; j < n; j++) {
      const up = [...c], dn = [...c], dc = 1e-5 * c[j];
      up[j] += dc;
      dn[j] -= dc;
      const a = at(up, s), b = at(dn, s), row = Math.max(...Array.from({ length: n }, (_, i) => Math.abs(out.dgdc[i * n + j])));
      assert.ok(Math.abs((a.q - b.q) / (2 * dc) - out.dqdc[j]) <= 1e-7 * Math.abs(out.dqdc[j]) + 1e-30, `${name} at s = ${s}: ∂q/∂c[${j}]`);
      for (let i = 0; i < n; i++) assert.ok(Math.abs((a.gam[i] - b.gam[i]) / (2 * dc) - out.dgdc[i * n + j]) <= 1e-7 * row, `${name} at s = ${s}: ∂Γ[${i}]/∂c[${j}]`);
    }
  }
});

// A closed blocking cell, Pt | 0.1 M KCl, 10 µm | Pt, a Stern layer of 0.2 F/m² (zeroCharge 0) and
// a diffuse layer on each face: by symmetry the bulk sits at V/2 and the left plane at V/2 + ψ,
// with C (V/2 − ψ) = √(8εRTc) sinh(Fψ/2RT), and the salt the layers take comes out of the bulk:
// c = c0 − (Γ₊(ψ) + Γ₊(−ψ))/L. Solved with mpmath (the answer sheet's cell.py).
const CELL = [
  { V: 0.05, sigma: 0.003917043109234853, psi: 0.005414784453825734, cb: 99.99978629783573, gK: -1.9230134698271056e-8, gCl: 2.1367156341042578e-8 },
  { V: 0.5, sigma: 0.04032406553579317, psi: 0.04837967232103414, cb: 99.98166086158574, gK: -1.1726905235359616e-7, gCl: 3.0066043649615234e-7 },
  { V: 1.0, sigma: 0.08400739467109519, psi: 0.07996302664452405, cb: 99.94326675549475, gK: -1.5167142321938999e-7, gCl: 7.190038682718826e-7 },
];
const stern = { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 } };
const cell = (V) => build({ library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')], stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), stern, layer('water', 10e-6, { c0: { 'K+': 100, 'Cl-': 100 } }), stern, layer('Pt', 1e-6), ohmic(0, ['e-'])], grid: { hmin: 0.5e-6, hmax: 0.5e-6 } });
// Each ion's amount (mol/m²): the boxes' and what the diffuse layers hold.
const amounts = (dev, s) => Object.fromEntries(['K+', 'Cl-'].map((name) => {
  let a = 0;
  s.c[name].forEach((c, g) => Number.isFinite(c) && (a += dev.model.grid.vol[g] * c));
  for (const itf of s.interfaces) for (const d of Object.values(itf.diffuse ?? {})) a += d.excess[name];
  return [name, a];
}));

test('a closed blocking cell with diffuse layers: the charge, drop and excesses of Stern and Gouy–Chapman in series, its salt kept', () => {
  const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
  for (const t of CELL) {
    const dev = new Device(cell(t.V)), s = dev.solve(), left = s.interfaces[0].diffuse.right, right = s.interfaces[1].diffuse.left;
    assert.ok(s.converged && s.warnings.length === 0, s.warnings.join('\n'));
    assert.ok(rel(s.interfaces[0].D, t.sigma) < 1e-9, `${t.V} V: σ ${s.interfaces[0].D} for ${t.sigma}`);
    assert.ok(rel(left.psi, t.psi) < 1e-9 && rel(right.psi, -t.psi) < 1e-9, `${t.V} V: ψ ${left.psi}, ${right.psi} for ±${t.psi}`);
    assert.ok(rel(left.excess['K+'], t.gK) < 1e-9 && rel(left.excess['Cl-'], t.gCl) < 1e-9, `${t.V} V: Γ ${JSON.stringify(left.excess)}`);
    assert.ok(rel(right.excess['K+'], t.gCl) < 1e-9 && rel(right.excess['Cl-'], t.gK) < 1e-9, `${t.V} V: Γ ${JSON.stringify(right.excess)} on the right`);
    assert.ok(rel(s.c['K+'][Math.floor(s.x.length / 2)], t.cb) < 1e-12, `${t.V} V: the bulk's ${s.c['K+'][Math.floor(s.x.length / 2)]} for ${t.cb}`);
    const a = amounts(dev, s);
    assert.ok(rel(a['K+'], 1e-3) < 1e-12 && rel(a['Cl-'], 1e-3) < 1e-12, `${t.V} V: amounts ${JSON.stringify(a)} for 1e-3 each`);
    assert.ok(Math.abs(s.charge) < 1e-9 * t.sigma, `${t.V} V: the device's charge ${s.charge} (its ions' ~100 C/m², summed)`);
  }
  // In time: from 0 V, stepped to 0.5 V, the layers charge through the bulk and draw salt from it,
  // which diffuses back over L²/D: every step keeps both amounts, and the end is the steady state.
  const dev = new Device(cell(0));
  dev.solve();
  dev.set({ contacts: { left: { V: 0.5 } } });
  let s;
  for (const t of [1e-8, 1e-6, 1e-4, 1e-2, 1]) {
    s = dev.advance(t);
    const a = amounts(dev, s);
    assert.ok(s.converged && rel(a['K+'], 1e-3) < 1e-12 && rel(a['Cl-'], 1e-3) < 1e-12, `at ${t} s: amounts ${JSON.stringify(a)}`);
  }
  assert.ok(rel(s.interfaces[0].D, CELL[1].sigma) < 1e-9, `σ after 1 s: ${s.interfaces[0].D} for ${CELL[1].sigma}`);
});

test('a diffuse layer charged by a current of ions: what comes in is held, read by the change of what it feeds', () => {
  // K⁺ driven in at the left end of a closed KCl film, against a Pt face that blocks it: the only
  // way in for anything, so the left terminal's current is read by the K⁺ amount's change, the
  // layer's excess included (src/lawRows.js, circuitLaws); the metal's charge moves by −I t.
  const c = { 'K+': 100, 'Cl-': 100 }, I = 0.5;
  const dev = new Device(build({
    library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')],
    stack: [{ I: 0, terminal: 'K+', species: { 'K+': 'equilibrium' }, phi: 'bulk' }, layer('water', 10e-6, { c0: c }), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 } }, layer('Pt', 1e-6), ohmic(0, ['e-'])],
    grid: { hmin: 0.5e-6, hmax: 0.5e-6 },
  }));
  const s0 = dev.solve(), K0 = amounts(dev, s0)['K+'], D0 = s0.interfaces[0].D;
  dev.set({ contacts: { left: { I } } });
  for (const t of [1e-4, 1e-3, 1e-2]) {
    const s = dev.advance(t), dK = amounts(dev, s)['K+'] - K0;
    assert.ok(s.converged && Math.abs(dK / ((I * t) / FARADAY) - 1) < 1e-8, `at ${t} s: K⁺ in ${dK} for ${(I * t) / FARADAY}`);
    assert.ok(Math.abs(s.interfaces[0].D - D0 - I * t) < 1e-8 * I * t, `at ${t} s: the metal's charge ${-s.interfaces[0].D} from ${-D0}`);
  }
});

test('a diffuse layer charged over very long steps, read by the change of what it holds', () => {
  // The same through a port passing K⁺ at 1 pA/m², by steps of 1e10 s: read by its flux, the first
  // step diverges, and it's taken again read by the change of the K⁺ amount (src/lawRows.js),
  // whose start must count what the layer held then.
  const c = { 'K+': 100, 'Cl-': 100 }, I = 1e-12;
  const dev = new Device(build({
    library: [aqueous(['K+', 'Cl-'], { epsr: 0 }), metal('Pt')],
    stack: [{}, layer('water', 10e-6, { name: 'film', c0: c }), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, diffuse: { epsr: 78.4 } }, layer('Pt', 1e-6), ohmic(0, ['e-'])],
    ports: [{ name: 'src', region: 'film', from: 0, to: 1e-6, I: 0, species: { 'K+': { type: 'conductance', G: 1e3, offset: 0 } } }],
    grid: { hmin: 0.5e-6, hmax: 0.5e-6 },
  }));
  const K0 = amounts(dev, dev.solve())['K+'];
  dev.set({ ports: { src: { I } } });
  for (const t of [1e10, 2e10]) {
    const s = dev.step(1e10), dK = amounts(dev, s)['K+'] - K0;
    assert.ok(s.converged && Math.abs(dK / ((I * t) / FARADAY) - 1) < 1e-8, `at ${t} s: K⁺ in ${dK} for ${(I * t) / FARADAY}`);
  }
});

test('a crowded diffuse layer is the limit of the resolved one, and keeps what a closed cell holds', () => {
  // Pt | 0.1 M KCl on a lattice of 0.71 nm sites | a bath, a stiff Stern layer (1 F/m²) held 1 V
  // below the pzc, the plane full of K⁺: against the same solution resolved (ε = 78.4, its lattice
  // in the bulk too), which converges onto it at second order (1.4e-4, 3.6e-5, 8.7e-6 at cells
  // from 0.02, 0.01, 0.005 nm growing by 5%, 2.5%, 1.25%).
  const salt = { 'K+': 100, 'Cl-': 100 }, cMax = 1 / (AVOGADRO * 0.71e-9 ** 3);
  const cell = (V, resolved) => {
    const l = aqueous(['K+', 'Cl-'], { epsr: resolved ? 78.4 : 0 });
    l.materials.water.statistics = [{ type: 'lattice', species: ['K+', 'Cl-'], cMax }];
    return build({ library: [l, metal('Pt')], stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 1 }, zeroCharge: 0, ...(resolved ? {} : { diffuse: { epsr: 78.4 } }) }, layer('water', 2e-6, { c0: salt }), bath(salt, 'Cl-')], grid: resolved ? { hmin: 0.01e-9, hmax: 0.2e-6, ratio: 1.025 } : { hmin: 1e-7, hmax: 0.2e-6, ratio: 1.15 } });
  };
  const rest = new Device(cell(0, false)).solve(), pzc = rest.phi[rest.region.indexOf(1)];
  const s = new Device(cell(pzc - 1, false)).solve(), r = new Device(cell(pzc - 1, true)).solve(), d = s.interfaces[0].diffuse.right;
  assert.ok(s.converged && s.warnings.length === 0, s.warnings.join('\n'));
  assert.ok(Math.abs(d.plane['K+'] / cMax - 1) < 1e-6, `the plane's K⁺ ${d.plane['K+']} of ${cMax} sites`);
  const gap = (r.interfaces[0].D - s.interfaces[0].D) / s.interfaces[0].D;
  assert.ok(gap > 2e-5 && gap < 6e-5, `resolved σ ${r.interfaces[0].D} against ${s.interfaces[0].D}: ${gap}`);
  // A closed cell (Pt | 0.5 M KCl and 1 M of a neutral species, all on 0.6 nm sites | Pt) stepped to
  // 1.5 V: the planes fill, the neutral species is pushed out of them, and every amount is kept.
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 });
  lib.species.push({ name: 'S', z: 0, cRef: 1000 });
  lib.materials.water.species.S = { D: 1e-9, mu0: 0 };
  lib.materials.water.statistics = [{ type: 'lattice', species: ['K+', 'Cl-', 'S'], cMax: 1 / (AVOGADRO * 0.6e-9 ** 3) }];
  const stern1 = { phi: { type: 'capacitive', C: 1 }, zeroCharge: 0, diffuse: { epsr: 78.4 } }, c0 = { 'K+': 500, 'Cl-': 500, S: 1000 };
  const closed = (V) => build({ library: [lib, metal('Pt')], stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), stern1, layer('water', 10e-6, { c0 }), stern1, layer('Pt', 1e-6), ohmic(0, ['e-'])], grid: { hmin: 0.5e-6, hmax: 0.5e-6 } });
  const all = (dev, sol) => Object.fromEntries(['K+', 'Cl-', 'S'].map((name) => {
    let a = 0;
    sol.c[name].forEach((cc, g) => Number.isFinite(cc) && (a += dev.model.grid.vol[g] * cc));
    for (const itf of sol.interfaces) for (const dd of Object.values(itf.diffuse ?? {})) a += dd.excess[name];
    return [name, a];
  }));
  const dev = new Device(closed(0)), a0 = all(dev, dev.solve());
  dev.set({ contacts: { left: { V: 1.5 } } });
  let sol;
  for (const t of [1e-6, 1e-2, 1]) {
    sol = dev.advance(t);
    const a = all(dev, sol);
    for (const k of ['K+', 'Cl-', 'S']) assert.ok(sol.converged && Math.abs(a[k] / a0[k] - 1) < 1e-12, `at ${t} s: ${k} ${a[k]} from ${a0[k]}`);
  }
  assert.ok(sol.interfaces[0].diffuse.right.excess.S < -5e-7, `S displaced from the planes: ${sol.interfaces[0].diffuse.right.excess.S}`);
  const st = new Device(closed(1.5)).solve();
  assert.ok(Math.abs(sol.interfaces[0].D / st.interfaces[0].D - 1) < 1e-8, `σ after 1 s ${sol.interfaces[0].D}, steady ${st.interfaces[0].D}`);
});

test("a reaction at a crowded plane sees what the resolved lattice's does, where the layer's own transport doesn't count", () => {
  // The ferri/ferrocyanide electrode on a lattice of 0.71 nm sites, a Stern layer of 1 F/m² at
  // zeroCharge 0.25 V, its plane half full of Cl⁻: R_dc as the resolved electrode's (ε = 78.4,
  // its lattice in the bulk) to 1.3e-5. (The layer is taken in equilibrium across itself, its own
  // transport left out: about λ_D c/(L c_plane) of the cell's, which grows where a reacting species
  // is driven out of the plane, 3.6e-4 at 6% of the bulk, or a packed plane leaves no room to
  // move, 6% at 99.7% full.)
  const c = { FeIII: 5, FeII: 5, 'Cl-': 1000, 'K+': 1035 }, cMax = 1 / (AVOGADRO * 0.71e-9 ** 3);
  const def = (grid, resolved) => {
    const d = ferri(30e-9, { diffuse: resolved ? null : { epsr: 78.4 }, epsr: resolved ? 78.4 : 0, zeroCharge: 0.25 });
    d.materials.water.statistics = [{ type: 'lattice', species: ['FeIII', 'FeII', 'K+', 'Cl-'], cMax }];
    d.interfaces[0].phi.C = 1;
    d.grid = grid;
    return d;
  };
  const g = rdc(def({ hmin: 30e-9, hmax: 2e-6, ratio: 1.2 }, false)), r = rdc(def({ hmin: 0.01e-9, hmax: 2e-6, ratio: 1.025 }, true));
  assert.ok(Math.abs(r - g) < 3e-5 * g, `R_dc ${g} Ω·cm², resolved ${r}`);
  const d = new Device(def({ hmin: 30e-9, hmax: 2e-6, ratio: 1.2 }, false)).solve().interfaces[0].diffuse.right;
  assert.ok(d.plane['Cl-'] > 0.4 * cMax, `the plane's Cl⁻ ${d.plane['Cl-']} of ${cMax}`);
});

test("the cell's low-frequency capacitance is Stern's and Gouy–Chapman's in series", () => {
  const RT = GAS_CONSTANT * 298.15, kappa = Math.sqrt((2 * FARADAY ** 2 * 100) / (78.4 * EPS0 * RT));
  const Cdl = 1 / (1 / 0.2 + 1 / (78.4 * EPS0 * kappa)); // each face at its zero charge: C_d = εκ
  const dev = new Device(cell(0));
  dev.solve();
  const fs = [1e-2, 1e2], { Z } = dev.impedance(fs);
  fs.forEach((f, k) => {
    const C = -1 / (2 * Math.PI * f * Z.im[k]);
    assert.ok(Math.abs(C - Cdl / 2) < 1e-6 * Cdl, `${f} Hz: C ${C} for ${Cdl / 2} (two faces in series)`);
  });
});

// Pt | 5 mM ferri/ferrocyanide in 1 M KCl, the Stern layer charged at rest (zeroCharge 0, the
// couple at +0.36 V: σ ≈ 0.07 C/m²). Without a diffuse layer the strictly neutral edge cell holds
// that charge, and the reaction sees its composition: the DC resistance moved from 56.6 to 47.4
// Ω·cm² between cells of 300 and 1 nm. With one it doesn't depend on the grid, and it's the limit
// of the same electrode with its double layer resolved (ε = 78.4, cells to 0.0125 nm). The real
// Frumkin effect is −18% here: the face draws the anions in, and only uncharged at rest does the
// electrode have the resistance without it (56.69).
const ferri = (hmin, { diffuse = { epsr: 78.4 }, epsr = 0, zeroCharge = 0 } = {}) => {
  const c = { FeIII: 5, FeII: 5, 'Cl-': 1000, 'K+': 1035 };
  return build({
    library: [{ species: [{ name: 'FeIII', z: -3, cRef: 1000 }, { name: 'FeII', z: -4, cRef: 1000 }, { name: 'K+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
      materials: { water: { epsr, species: { FeIII: { D: 0.72e-9, mu0: 0 }, FeII: { D: 0.67e-9, mu0: -FARADAY * 0.36 }, 'K+': { D: 1.957e-9, mu0: 0 }, 'Cl-': { D: 2.032e-9, mu0: 0 } } } } }, metal('Pt')],
    stack: [ohmic({ I: 0 }, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge, ...(diffuse ? { diffuse } : {}), reactions: [{ ...half('FeIII + e- = FeII'), k0: 5e-2, alpha: 0.5 }] }, layer('water', 30e-6, { c0: c }), bath(c, 'Cl-')],
    grid: { hmin, hmax: 2e-6, ratio: 1.2 },
  });
};
// The DC resistance (Ω·cm²) by ±1 mA/m²
const rdc = (def) => {
  const V = (I) => {
    const d = structuredClone(def);
    d.contacts.left.I = I;
    const s = new Device(d).solve();
    assert.ok(s.converged && s.warnings.length === 0, s.warnings.join('\n'));
    return s.terminals.left.V;
  };
  return ((V(1e-3) - V(-1e-3)) / 2e-3) * 1e4;
};

test("an electrode's DC resistance beside a diffuse layer doesn't depend on the grid, and is the resolved double layer's", () => {
  const R = [300e-9, 30e-9, 1e-9].map((h) => rdc(ferri(h)));
  assert.ok(Math.abs(R[1] - R[0]) < 1e-9 * R[0] && Math.abs(R[2] - R[0]) < 1e-9 * R[0], `R_dc ${R.join(', ')} Ω·cm²`);
  assert.ok(Math.abs(R[0] - 46.399) < 1e-3, `R_dc ${R[0]}`);
  // the resolved electrode converges onto it (8.9e-4, 2.9e-4, 1.2e-4, 6.1e-5 at 0.1 nm to 0.0125 nm)
  const resolved = rdc(ferri(0.0125e-9, { diffuse: null, epsr: 78.4 }));
  assert.ok(resolved > R[0] && resolved - R[0] < 1e-4 * R[0], `resolved ${resolved} for ${R[0]}`);
  // uncharged at rest, the face draws nothing in: the resistance without a Frumkin effect, on any grid, either way
  const rest = new Device(ferri(30e-9)).solve(), g = rest.region.indexOf(1), zc = rest.terminals.left.V - rest.phi[g];
  const bare = rdc(ferri(30e-9, { zeroCharge: zc })), plain = rdc(ferri(300e-9, { diffuse: null, zeroCharge: zc }));
  assert.ok(Math.abs(bare - 56.688) < 1e-3 && Math.abs(bare - plain) < 1e-6 * bare, `uncharged ${bare}, without a diffuse layer ${plain}`);
});

test('the reaction at a diffuse layer sees the plane, which the solution reports', () => {
  const s = new Device(ferri(30e-9)).solve(), d = s.interfaces[0].diffuse.right, g = s.region.indexOf(1);
  assert.ok(d.psi > 0 && d.charge < 0, `ψ ${d.psi} V, charge ${d.charge} C/m²: a positive metal's layer is negative`);
  const VT = (GAS_CONSTANT * 298.15) / FARADAY;
  for (const [name, z] of [['FeIII', -3], ['FeII', -4], ['K+', 1], ['Cl-', -1]]) {
    assert.ok(Math.abs(d.plane[name] - s.c[name][g] * Math.exp((-z * d.psi) / VT)) < 1e-12 * d.plane[name], `${name} at the plane`);
  }
  assert.ok(Math.abs(d.charge + s.interfaces[0].D) < 1e-9 * Math.abs(d.charge), `the layer holds the metal's charge: ${d.charge} against ${s.interfaces[0].D}`);
  assert.match(describe(ferri(30e-9)), /a diffuse layer below the grid \(Gouy–Chapman, εr 78\.4\) on the region 1 side/);
});

test('a diffuse layer is refused where it would be wrong or meaningless', () => {
  const lib = aqueous(['K+', 'Cl-'], { epsr: 0 }), c = { 'K+': 100, 'Cl-': 100 };
  const make = (face, { epsr = 0, mod } = {}) => {
    const l = aqueous(['K+', 'Cl-'], { epsr });
    mod?.(l);
    return build({ library: [l, metal('Pt')], stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, ...face }, layer('water', 1e-6, { c0: c }), bath(c, 'Cl-')], grid: { hmin: 1e-7, hmax: 1e-7 } });
  };
  const refused = (def, pattern) => assert.throws(() => new Device(def), (e) => e instanceof DeviceError && pattern.test(e.message), pattern.source);
  refused(make({ diffuse: {} }), /diffuse\.epsr must be a positive number/);
  refused(make({ diffuse: { epsr: 78, a: 1 } }), /diffuse\.a: not a field here/);
  refused(make({ diffuse: { epsr: 78 }, phi: 'neutral', zeroCharge: undefined }), /a neutral face holds no charge, so it has no diffuse layer/);
  refused(make({ diffuse: { epsr: 78 } }, { epsr: 78 }), /neither side is a strictly neutral.*the grid resolves the double layer itself/);
  refused(make({ diffuse: { epsr: 78 } }, { mod: (l) => (l.materials.water.statistics = [{ type: 'debye-huckel', species: ['K+', 'Cl-'], epsr: 78.4 }]) }), /'K\+' in 'water' has statistics a diffuse layer can't take: its ions must all be ideal \(a dilute layer\) or all on one lattice model/);
  refused(make({ diffuse: { epsr: 78 } }, { mod: (l) => (l.materials.water.statistics = [{ type: 'lattice', species: ['K+'], cMax: 5000 }]) }), /'Cl-' in 'water' has statistics a diffuse layer can't take/);
  refused(make({ diffuse: { epsr: 78 } }, { mod: (l) => (l.materials.water.species['K+'].D = 0) }), /'K\+' is charged and can't move in 'water' \(D = 0\)/);
});

// A charged membrane (C 0.01 F/m², a sheet of −0.02 C/m² on its right face) between two strictly
// neutral solutions, a diffuse layer on each side; K⁺ crosses by a permeability gated by n
// (sigmoid α, β about −30 mV), the rest are blocked. Its answer sheet (membrane.py, mpmath): the
// drops from Gauss and the Stern law alone, q_L(ψ_L) = D, q_R(ψ_R) = −(D + σ), D = −C V_m with
// V_m = Δφ + ψ_R − ψ_L between the planes; the gate at α/(α + β) of V_m; Goldman–Hodgkin–Katz
// between the planes' concentrations c e^{−zFψ/RT}.
const MEMBRANE = [
  { V: -0.1, D: 0.00110145342394159, psiL: -0.0014884579314256, psiR: -0.0207976995439835, Vm: -0.110145342394159, n: 0.000330512920061092, NK: 1.40028794604074e-12 },
  { V: -0.06, D: 0.000710804212490204, psiL: -0.000960629355138855, psiR: -0.0212049494436698, Vm: -0.0710804212490204, n: 0.0161740304336446, NK: -7.11608873808973e-10 },
  { V: 0, D: 0.000124810284894188, psiL: -0.000168686448995856, psiR: -0.0218135019386403, Vm: -0.0124810284894188, n: 0.852191928423054, NK: -2.11175716567407e-7 },
  { V: 0.05, D: -0.000363537632357396, psiL: 0.000491328861962143, psiR: -0.0223184772918532, Vm: 0.0363537632357396, n: 0.99868863742237, NK: -6.20685552170249e-7 },
];

test('a charged membrane between diffuse layers: its permeability and gates work between the planes', () => {
  const cL = { 'Na+': 100, 'K+': 5, 'Cl-': 105 }, cR = { 'Na+': 10, 'K+': 140, 'Cl-': 150 };
  const gate = { alpha: { type: 'sigmoid', rate: 100, midpoint: -0.03, scale: 0.01 }, beta: { type: 'sigmoid', rate: 100, midpoint: -0.03, scale: -0.01 } };
  for (const t of MEMBRANE) {
    const def = build({
      library: [aqueous(['Na+', 'K+', 'Cl-'], { epsr: 0 })],
      stack: [bath(cL, 'Cl-'), layer('water', 5e-6, { c0: cL }), { phi: { type: 'capacitive', C: 0.01 }, sheetCharge: -0.02, diffuse: { epsr: 78.4 }, gates: { n: gate }, species: { 'K+': { type: 'permeability', P: 1e-9, gates: { n: 1 } }, 'Na+': 'blocked', 'Cl-': 'blocked' } }, layer('water', 5e-6, { c0: cR }), bath(cR, 'Cl-', t.V)],
      grid: { hmin: 0.5e-6, hmax: 1e-6 },
    });
    const s = new Device(def).solve(), f = s.interfaces[0];
    assert.ok(s.converged && s.warnings.length === 0, s.warnings.join('\n'));
    for (const [what, got, want] of [['D', f.D, t.D], ['ψ left', f.diffuse.left.psi, t.psiL], ['ψ right', f.diffuse.right.psi, t.psiR], ['V between the planes', f.V, t.Vm], ['n', f.gates.n, t.n], ['K⁺ flux', f.N['K+'], t.NK]]) {
      assert.ok(Math.abs(got - want) < 1e-9 * Math.abs(want), `${t.V} V: ${what} ${got} for ${want}`);
    }
  }
});
