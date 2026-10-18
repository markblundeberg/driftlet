// A sweep for steady states that keep the wrong amounts: small devices combining the features at
// the edges of the steady system's conservation laws (floating and capacitive electrodes, two of
// them, surfaces, immobile species, face reactions, capacitive contacts), each solved directly
// and run long as a transient from the same start. Where the transient has settled (it moves no
// further over 100× as long), the two must agree.
//   node tools/lawsweep.js [--from k] [--to k] [--show k] [--out file.json]
// Prints one line per device that disagrees (or fails), and a tally. --show k prints device k's
// definition and both answers.
import { writeFileSync } from 'node:fs';
import { Device } from '../src/index.js';
import { build, layer, bath, aqueous } from '../src/kit.js';

const args = process.argv.slice(2), opt = (name, d) => (args.includes(name) ? args[args.indexOf(name) + 1] : d);

// The features, each a few options; a device is one choice of each.
const AXES = {
  left: ['closed', 'capK', 'heldK'], // the left end: blocked, a capacitor at I = 0 passing K⁺, or held passing K⁺
  port: ['none', 'redox', 'redoxSurf', 'surfOnly', 'S0surf'], // an electrode in the film
  drive: ['V', 'I0', 'R'], // its drive
  cap: [false, true], // and a capacitance on it
  second: ['none', 'redoxI0', 'membrane'], // a second port: another floating electrode, or a membrane capacitor passing K⁺
  face: ['none', 'open', 'makesP0'], // a face splitting the film: Fe³⁺ blocked there, or a reaction making an immobile P0 at it
};
// A second family, after the first's devices: an immobile ion and what carries it (a face that
// passes it, a port at I = 0 that links it, an electrode's surface it adsorbs on), beside a
// capacitive contact passing K⁺ (at I = 0, or behind a resistance) and K⁺'s complex, each where
// a cold review of src/laws.js found a law kept wrongly.
const AXES2 = {
  left: ['closed', 'capK', 'RcapK'], // the left end: blocked, a capacitor passing K⁺ at I = 0, or behind a resistance
  immobile: ['none', 'plain', 'faceLinked', 'portLinked', 'adsorbs'], // an immobile S⁺, and what carries it
  complex: [false, true], // K⁺ + L⁻ = KL, L⁻ and KL closed
  K: ['trace', 'some'], // K⁺ in the film: 1e-6 or 1 mol/m³ (which law comes first, and how a step reads a current)
};
const families = [AXES, AXES2].map((axes) => ({ axes, keys: Object.keys(axes), size: Object.values(axes).reduce((n, a) => n * a.length, 1) }));
const total = families.reduce((n, f) => n + f.size, 0);
const choice = (idx) => {
  let fam = 0;
  while (idx >= families[fam].size) idx -= families[fam++].size;
  const { axes, keys } = families[fam], c = { family: fam };
  for (const k of keys) {
    c[k] = axes[k][idx % axes[k].length];
    idx = Math.floor(idx / axes[k].length);
  }
  return c;
};

// The film is strictly neutral water (ε = 0), as an electrode spread through it needs; a
// capacitor (the left one, or the membrane) sits in a thin layer of polar water (ε = 80) at the
// left, on a grid fine enough for its double layer.
const G = 1e-5, L = 10e-6, Lp = 1e-6;
function deviceA(c) {
  const polar = c.left === 'capK' || c.second === 'membrane';
  const lib = aqueous(['Na+', 'K+', 'Cl-', 'OH-', 'Fe2+', 'Fe3+'], { epsr: 0 });
  lib.species.push({ name: 'S0', z: 0 }, { name: 'P0', z: 0 });
  if (c.port !== 'none' || c.second === 'redoxI0') lib.species.push({ name: 'e-', z: -1 }); // (an electrode's carrier, with one)
  lib.materials.water.species.S0 = { D: 0, mu0: 0, cRef: 1 };
  // Every reaction starts at its equilibrium (far from one, a first step is violent and the
  // transient fails, leaving nothing to compare): P0, the face's immobile product, and the
  // adsorbates' μ° from the film's levels with the electrode at the iron couple's.
  const RT0 = 8.314462618 * 298.15, w0 = lib.materials.water.species, a = { 'Fe3+': 1e-3, 'Cl-': 0.112, 'OH-': 1e-3 };
  lib.materials.water.species.P0 = { D: 0, mu0: w0['Fe3+'].mu0 + 3 * w0['Cl-'].mu0 + RT0 * Math.log(a['Fe3+'] * a['Cl-'] ** 3), cRef: 1 };
  const eMu = w0['Fe2+'].mu0 - w0['Fe3+'].mu0; // (μ̄_e at a 1:1 couple, less −Fφ, which cancels below)
  const mu0B = (x) => w0[x].mu0 + RT0 * Math.log(a[x]) - eMu - RT0 * Math.log(0.1 / 0.2); // A(ads) + X = B(ads) + e⁻, θ_A 0.2, θ_B 0.1
  const mu0C = w0['OH-'].mu0 + RT0 * Math.log(a['OH-']) - eMu - RT0 * Math.log(0.2 / 0.8); // S0 + OH⁻ = C(ads) + e⁻, θ_C 0.2 of 1
  lib.materials.polar = { ...lib.materials.water, epsr: 80 };
  const film = { 'Na+': 100, 'K+': 10, 'OH-': 1, 'Fe2+': 1, 'Fe3+': 1, S0: 1, P0: 1 };
  film['Cl-'] = film['Na+'] + film['K+'] + 2 * film['Fe2+'] + 3 * film['Fe3+'] - film['OH-'];
  // (K⁺ passed at the film's own level, 10 mM: what a start near its steady state wants)
  // (its species voltage there: μ/F plus φ, which the bath's Cl⁻ at V = 0 sets, (μ°_Cl + RT ln a)/F)
  const RT = 8.314462618 * 298.15, sp = lib.materials.water.species;
  const VK = (sp['K+'].mu0 + RT * Math.log(film['K+'] / 1000) + sp['Cl-'].mu0 + RT * Math.log(101 / 1000)) / 96485.33212;
  const K = { 'K+': { type: 'conductance', G: 1e6, offset: VK } };
  const left = { closed: {}, capK: { I: 0, zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: K }, heldK: { V: 0, phi: 'bulk', species: K } }[c.left];
  const all = (law) => Object.fromEntries(['Na+', 'K+', 'Cl-', 'OH-', 'Fe2+'].map((s) => [s, law]));
  const stack = [left];
  if (polar) stack.push(layer('polar', Lp, { name: 'skin', c0: film, grid: { hmin: 2e-10, hmax: 5e-8, ratio: 1.25 } }), { phi: 'neutral', species: { ...all('equilibrium'), 'Fe3+': 'equilibrium', S0: 'blocked', P0: 'blocked' } });
  if (c.face === 'none') stack.push(layer('water', L, { name: 'film', c0: film }));
  else stack.push(
    layer('water', L / 2, { name: 'film', c0: film }),
    { species: { ...all('equilibrium'), 'Fe3+': 'blocked', S0: 'blocked', P0: 'blocked' }, ...(c.face === 'makesP0' ? { reactions: [{ equation: 'Fe3+(right) + 3 Cl-(right) = P0(right)', k0: 1e-3, alpha: 0.5 }] } : {}) },
    layer('water', L / 2, { name: 'far', c0: film }),
  );
  stack.push(bath({ 'Na+': 100, 'Cl-': 101, 'OH-': 1, 'Fe2+': 1 }, 'Cl-'));
  const ports = [];
  const w = c.face === 'none' ? L : L / 2;
  if (c.port !== 'none') {
    const surface = c.port === 'redoxSurf' || c.port === 'surfOnly' ? { Aads: { mu0: 0, capacity: G, theta0: 0.2 }, Bads: { mu0: mu0B(c.port === 'redoxSurf' ? 'Cl-' : 'OH-'), capacity: G, theta0: 0.1 } } : c.port === 'S0surf' ? { Cads: { mu0: mu0C, capacity: G, theta0: 0.2 } } : undefined;
    const reactions = [
      ...(c.port.startsWith('redox') ? [{ equation: 'Fe2+ = Fe3+ + e-', k0: 1e-4, alpha: 0.5 }] : []),
      ...(c.port === 'redoxSurf' ? [{ equation: 'Aads + Cl- = Bads + e-', k0: 1e-4, alpha: 0.5 }] : []),
      ...(c.port === 'surfOnly' ? [{ equation: 'Aads + OH- = Bads + e-', k0: 1e-4, alpha: 0.5 }] : []),
      ...(c.port === 'S0surf' ? [{ equation: 'S0 + OH- = Cads + e-', k0: 1e-4, alpha: 0.5 }] : []),
    ];
    // (held near the iron couple's level, where the film's 1:1 iron sits: E° = (μ°_Fe³⁺ − μ°_Fe²⁺)/F,
    // on the bath's φ scale; far from it, all the iron goes one way and the steady state is extreme)
    const E = (sp['Fe3+'].mu0 - sp['Fe2+'].mu0 + sp['Cl-'].mu0 + RT * Math.log(101 / 1000)) / 96485.33212;
    const drive = { V: { V: E + 0.02 }, I0: { I: 0 }, R: { V: E + 0.02, R: 10 } }[c.drive];
    ports.push({ name: 'p1', region: 'film', from: 0, to: 0.4 * w, ...drive, terminal: 'e-', area: 1e4, ...(c.cap ? { capacitance: { C: 0.2, zeroCharge: 0 } } : {}), ...(surface ? { surface } : {}), reactions });
  } else if (c.drive !== 'V' || c.cap) return null; // (drives and capacitances of no port: the same device again)
  if (c.second === 'redoxI0') ports.push({ name: 'p2', region: 'film', from: 0.6 * w, to: w, I: 0, terminal: 'e-', area: 1e4, reactions: [{ equation: 'Fe2+ = Fe3+ + e-', k0: 1e-4, alpha: 0.5 }] });
  if (c.second === 'membrane') ports.push({ name: 'mem', region: 'skin', I: 0, area: 1e6, capacitance: { C: 0.1, zeroCharge: -0.05 }, species: { 'K+': { type: 'conductance', G: 1e12, offset: VK } } });
  return build({ T: 298.15, library: [lib], stack, ports, grid: { hmin: 1e-6, hmax: 1e-6 } });
}

// The second family's: a film of strictly neutral water (ε = 0) cut by a face at its middle,
// with an immobile S⁺ (D = 0) where it has one, after a thin polar skin where the left end is a
// capacitor.
function device2(c) {
  const ions = ['Na+', 'K+', 'Cl-'], lib = aqueous(ions, { epsr: 0 });
  const RT = 8.314462618 * 298.15, w = lib.materials.water.species;
  const cK = c.K === 'trace' ? 1e-6 : 1, cL = 10, cKL = 1;
  if (c.immobile !== 'none') {
    lib.species.push({ name: 'S+', z: 1 });
    w['S+'] = { D: 0, mu0: 0, cRef: 1 };
  }
  if (c.immobile === 'adsorbs') lib.species.push({ name: 'e-', z: -1 }); // (the electrode's carrier)
  if (c.complex) {
    lib.species.push({ name: 'L-', z: -1 }, { name: 'KL', z: 0 });
    w['L-'] = { D: 1e-9, mu0: 0, cRef: 1 };
    // (KL at its equilibrium with the film's K⁺ and L⁻)
    w.KL = { D: 1e-9, mu0: w['K+'].mu0 + w['L-'].mu0 + RT * Math.log((cK / lib.species.find((sp) => sp.name === 'K+').cRef) * cL / cKL), cRef: 1 };
  }
  const has = (name) => lib.species.some((sp) => sp.name === name), only = (o) => Object.fromEntries(Object.entries(o).filter(([k]) => has(k)));
  lib.materials.polar = { ...lib.materials.water, epsr: 80 };
  const film = (S) => {
    const f = { 'Na+': 100, 'K+': cK, ...(c.complex ? { 'L-': cL, KL: cKL } : {}), ...(c.immobile !== 'none' ? { 'S+': S } : {}) };
    f['Cl-'] = f['Na+'] + f['K+'] + (f['S+'] ?? 0) - (f['L-'] ?? 0);
    return f;
  };
  // (K⁺ passed at the film's own level, as in the first family: far from it, the contact pumps K⁺
  // to an extreme state that the steady solve can't reach)
  const VK = (w['K+'].mu0 + RT * Math.log(cK / 1000) + w['Cl-'].mu0 + RT * Math.log(100 / 1000)) / 96485.33212;
  const K = { 'K+': { type: 'conductance', G: 1e6, offset: VK } };
  const cap = { zeroCharge: -0.05, phi: { type: 'capacitive', C: 0.5 }, species: K };
  const left = { closed: {}, capK: { I: 0, ...cap }, RcapK: { V: 0, R: 1e-3, ...cap } }[c.left];
  const pass = only({ 'Na+': 'equilibrium', 'K+': 'equilibrium', 'Cl-': 'equilibrium', 'L-': 'blocked', KL: 'blocked', 'S+': 'blocked' });
  const stack = [left];
  if (c.left !== 'closed') stack.push(layer('polar', 1e-6, { name: 'skin', c0: film(1), grid: { hmin: 2e-10, hmax: 5e-8, ratio: 1.25 } }), { phi: 'neutral', species: pass });
  stack.push(
    layer('water', 5e-6, { name: 'film', c0: film(1) }),
    { species: { ...pass, ...(c.immobile === 'faceLinked' ? { 'S+': { type: 'conductance', G: 1e-6 } } : {}) } },
    layer('water', 5e-6, { name: 'far', c0: film(0.1) }),
    bath({ 'Na+': 100, 'Cl-': 100 }, 'Cl-'),
  );
  const ports = [];
  if (c.immobile === 'portLinked') ports.push({ name: 'p', region: 'film', from: 0, to: 3e-6, I: 0, species: { 'S+': { type: 'conductance', G: 1e3, offset: 0 } } });
  if (c.immobile === 'adsorbs') {
    ports.push({ name: 'metal', region: 'film', V: 0.05, terminal: 'e-', area: 1e4,
      surface: { Aads: { mu0: 0, capacity: 1e-5, theta0: 0.2 }, Bads: { mu0: 3e3, capacity: 1e-5, theta0: 0.1 } },
      reactions: [{ equation: 'S+ + Aads + e- = Bads', k0: 1e-4, alpha: 0.5 }] });
  }
  return build({ T: 298.15, library: [lib], stack, ports, ...(c.complex ? { bulkReactions: [{ equation: 'K+ + L- = KL', kf: { water: 1e3, polar: 1e3 } }] } : {}), grid: { hmin: 1e-6, hmax: 1e-6 } });
}
const device = (c) => (c.family === 0 ? deviceA(c) : device2(c));

// What a device keeps: each species' amount, each surface species' sites, and each floating
// terminal's voltage.
function measures(d, s) {
  const v = d.model.grid.vol, m = {};
  for (const sp of d.model.species) {
    const c = s.c[sp.name];
    if (!c) continue;
    m[`n ${sp.name}`] = c.reduce((t, x, g) => t + (Number.isFinite(x) ? v[g] * x : 0), 0);
  }
  s.ports.forEach((p) => {
    for (const [name, th] of Object.entries(p.coverage ?? {})) m[`θ ${name}`] = p.x.reduce((t, x, j) => t + v[s.x.indexOf(x)] * p.area[j] * th[j], 0);
  });
  for (const side of ['left', 'right']) if (s.contacts[side].D) m[`D ${side}`] = s.contacts[side].D;
  for (const k of d.solver.floating) {
    const name = k < 2 ? ['left', 'right'][k] : d.model.ports[k - 2].name;
    m[`V ${name}`] = s.terminals[name].V;
  }
  return m;
}
// (Voltages to a mV; charges to a µC/m²; amounts against themselves, floored at 1e-10 of the
// device's whole inventory, below which they're round-off.)
const differ = (a, b) => {
  let worst = { what: null, rel: 0 };
  const inventory = Object.keys(a).filter((k) => k.startsWith('n ')).reduce((t, k) => t + Math.abs(a[k]), 0);
  for (const k of Object.keys(a)) {
    const scale = k.startsWith('V ') ? 1e-3 : k.startsWith('D ') ? Math.max(Math.abs(a[k]), Math.abs(b[k]), 1e-6) : Math.max(Math.abs(a[k]), Math.abs(b[k]), 1e-10 * inventory);
    const rel = Math.abs(a[k] - b[k]) / scale;
    if (rel > worst.rel) worst = { what: k, rel, a: a[k], b: b[k] };
  }
  return worst;
};

const from = +opt('--from', 0), to = Math.min(+opt('--to', total), total), show = opt('--show');
const tally = { ok: 0, skipped: 0, invalid: 0, transientFailed: 0, unsettled: 0, solveFailed: 0, differ: 0 }, found = [];
for (let idx = show !== undefined ? +show : from; idx < (show !== undefined ? +show + 1 : to); idx++) {
  const c = choice(idx), def = device(c);
  if (!def) {
    tally.skipped++;
    continue;
  }
  const tag = `${idx} ${families[c.family].keys.map((k) => c[k]).join(' ')}`;
  try {
    new Device(def);
  } catch (e) {
    tally.invalid++;
    if (show !== undefined) console.log(tag, 'invalid:', e.message);
    continue;
  }
  try {
    const d = new Device(def), s = d.solve();
    const e = new Device(def), late = e.advance(1e6, { tol: 1e-9 });
    const ms = measures(d, s), ml = measures(e, late);
    if (show !== undefined) console.log(JSON.stringify(def, null, 1), '\nsolve', ms, '\nadvance', ml);
    if (!late.converged || !late.done) {
      tally.transientFailed++; // (nothing to compare with)
      continue;
    }
    if (!s.converged) {
      tally.solveFailed++;
      console.log(`${tag}: solve() did not converge`);
      found.push({ idx, c, kind: 'solveFailed' });
      continue;
    }
    const w = differ(ms, ml);
    if (w.rel < 1e-5) {
      tally.ok++;
      continue;
    }
    // Settled? 100× as long again moves it no further.
    const last = e.advance(1e8, { tol: 1e-9 }), later = measures(e, last), drift = differ(ml, later);
    if (show !== undefined) console.log('transient from 1e6 to 1e8 s:', last.converged, drift);
    if (!last.converged || !last.done || drift.rel > 1e-6) {
      tally.unsettled++;
      continue;
    }
    tally.differ++;
    console.log(`${tag}: ${w.what} ${w.a.toPrecision(6)} solved, ${w.b.toPrecision(6)} transient (${w.rel.toExponential(1)})`);
    found.push({ idx, c, kind: 'differ', worst: w });
  } catch (err) {
    tally.solveFailed++;
    console.log(`${tag}: threw ${err.message.slice(0, 120)}`);
    found.push({ idx, c, kind: 'threw', message: err.message });
  }
}
console.log(`of ${to - from} (${total} in all):`, JSON.stringify(tally));
if (opt('--out')) writeFileSync(opt('--out'), JSON.stringify(found, null, 1));
