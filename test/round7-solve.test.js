// The steady solve and its start, against four suspects the combos fuzzer's membrane and echem
// families turned up (20 Nov 2026): a converged root that wasn't one, a start that depended on the
// gauge, a floating port that overrode its solution's level, a law with a NaN weight, and a cell
// behind a resistance laid out from the wrong voltage. Each test's comment gives the old numbers;
// the right ones come from equilibrium, worked by hand from μ° and the compositions.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { aqueous, bath, build, layer, metal, ohmic } from '../src/kit.js';

const RT = GAS_CONSTANT * 298.15;
const levels = (lib, c) => (name) => lib.materials.water.species[name].mu0 + RT * Math.log(c[name] / 1000);

test("a myelinated axon's Na⁺, fed only through its nodes' gated channels, solves cold to where the channels balance, not to a level they leak from", () => {
  // (A K⁺ electrode on the left, the right end sealed, three HH nodes and two myelin internodes;
  // the membrane combos' case 158. The cold solve converged at φ = +0.424 V, Na⁺ and Cl⁻ at
  // 2.07e10 mol/m³, the h gates shut to 1.7e-12: Na⁺ lost 2.7e-15 mol/(m²·s) through them with
  // nothing coming in, 35 thermal units from balance, each row solved to 1e-23 of its largest
  // term. Listed in reverse, it solved to −8.2 mV; check() failed it on 'order'.)
  const def = JSON.parse(readFileSync(new URL('./fixtures/myelin-k-electrode.json', import.meta.url), 'utf8'));
  const d = new Device(def), v = d.solver, refused = [];
  const newton = v.newton;
  v.newton = function (dt, opts) {
    const r = newton.call(this, dt, opts);
    if (!r.converged && r.error) refused.push(r.error);
    return r;
  };
  const s = d.solve();
  assert.ok(s.converged, s.warnings.join('; '));
  // Steady, Na⁺ and Cl⁻ have nowhere to go but the channels, so at each node each is at its
  // bath's level: V_i inside = (μ° + RT ln(c_out/c_ref))/(zF) + V_port (V = 0 there).
  const rt = GAS_CONSTANT * def.T, mat = def.materials.axoplasm.species;
  for (const port of def.ports.filter((p) => p.species?.['Na+'])) {
    const g = s.x.findIndex((x, j) => s.region[j] === def.regions.findIndex((r) => r.name === port.region) && x > 0);
    for (const ion of ['Na+', 'Cl-']) {
      const z = ion === 'Na+' ? 1 : -1, out = (mat[ion].mu0 + rt * Math.log(port.bath.c[ion] / 1000)) / (z * FARADAY);
      assert.ok(Math.abs(s.V[ion][g] - out) < 1e-6, `${port.name}: ${ion} at ${s.V[ion][g]} V, its bath's ${out}`);
    }
  }
  assert.ok(s.c['Na+'][s.x.length >> 1] < 1e3, `Na⁺ at ${s.c['Na+'][s.x.length >> 1]} mol/m³`);
  // (Newton reaches the false root first, from cold, and it's refused, saying why; the retry
  // with smaller updates finds the right one.)
  assert.ok(refused.some((e) => /^Na\+ in node0–node2 loses .* thermal units from where they do$/.test(e)), `the false root refused, saying why: ${refused[0]}`);
});

test('two solutions in a row behind a bath, or one behind a blocking metal, start at their electrode\'s level in any gauge', () => {
  // (Bath | CuSO₄ + NaCl in two regions | Cu at rest, every voltage moved by φ: from φ = 2 V
  // (1.5 V converged) the cold start put both solutions at φ = 0, 83 RT from the copper's
  // level, and the solve failed, where one region, or the mirror image, converged at any φ;
  // so with the bath read by Na⁺, φ = 2.75 V. And Pt (blocking, a Stern layer) | KNO₃ +
  // Zn(NO₃)₂ | Zn at rest, ε = 78.3: started at φ = 0, 0.83 V from the zinc's level, it failed;
  // moved to φ = 0, converged. The echem combos' cases 212 and 341.)
  const c = { 'Cu2+': 40, 'Na+': 250, 'Cl-': 80, 'SO42-': 125 };
  const lib = aqueous(['Cu2+', 'Na+', 'Cl-', 'SO42-'], { epsr: 0 }), mu = levels(lib, c);
  const copper = { reactions: [{ equation: 'Cu2+ + 2 e- = Cu(s)', fixed: { 'Cu(s)': 0 }, k0: 3e-5, alpha: 0.5 }] };
  for (const phi of [2, 3]) {
    const VCu = mu('Cu2+') / (2 * FARADAY) + phi;
    const water = [layer('water', 0.5e-6, { name: 'w0', c0: c }), {}, layer('water', 0.5e-6, { name: 'w1', c0: c })];
    const s = new Device(build({ library: [lib, metal('Cu')], stack: [bath(c, phi), ...water, copper, layer('Cu', 1e-6), ohmic(VCu, ['e-'])], grid: { hmin: 5e-9, hmax: 2e-8, ratio: 1.15 } })).solve();
    assert.ok(s.converged, `φ = ${phi} V: ${s.warnings[0]}`);
    for (const [g, x] of s.x.entries()) if (s.region[g] < 2 && Number.isFinite(s.phi[g])) assert.ok(Math.abs(s.phi[g] - phi) < 1e-9, `φ = ${phi} V: ${s.phi[g]} at x = ${x}`);
  }
  const cz = { 'Zn2+': 6.3, 'K+': 41, 'NO3-': 53.6 };
  const libz = aqueous(['Zn2+', 'K+', 'NO3-'], { epsr: 78.3 }), phiZ = -levels(libz, cz)('Zn2+') / (2 * FARADAY);
  const zinc = { phi: { type: 'capacitive', C: 0.29 }, zeroCharge: -phiZ, reactions: [{ equation: 'Zn2+ + 2 e- = Zn(s)', fixed: { 'Zn(s)': 0 }, k0: 5e-3, alpha: 0.65 }] };
  const s = new Device(build({
    library: [libz, metal('Pt'), metal('Zn')],
    stack: [ohmic(0, ['e-']), layer('Pt', 1e-6), { phi: { type: 'capacitive', C: 0.25 }, zeroCharge: -phiZ }, layer('water', 9.6e-6, { c0: cz }), zinc, layer('Zn', 1e-6), ohmic(0, ['e-'])],
    grid: { hmin: 3e-10, hmax: 1.4e-7, ratio: 1.1 },
  })).solve();
  assert.ok(s.converged, s.warnings[0]);
  const mid = s.x.length >> 1;
  assert.ok(Math.abs(s.phi[mid] - phiZ) < 1e-9, `φ ${s.phi[mid]} V, the zinc's level puts it at ${phiZ}`);
});

test("a floating porous electrode beside an unlike one starts at its own couple's level, where the solution's electrode put the solution", () => {
  // (Zn | Zn²⁺, Fe³⁺, Fe²⁺, KClO₄ with a porous Fe³⁺/Fe²⁺ port at I = 0 | bath, at rest: the
  // port's V as it stood, 0 V, set the solution's φ (−0.725 V, the couple at 0 V) over the
  // zinc's 0.827, its mixed potential then found there; the solve failed, where the port held
  // at its couple's 1.5518 V converged. The echem combos' cases 207 and 291.)
  const c = { 'Zn2+': 6.3, 'Fe3+': 2.2, 'Fe2+': 12.4, 'K+': 500, 'ClO4-': 500 + 2 * 6.3 + 3 * 2.2 + 2 * 12.4 };
  const lib = aqueous(['Zn2+', 'Fe3+', 'Fe2+', 'K+', 'ClO4-'], { epsr: 0 }), mu = levels(lib, c);
  const phi = -mu('Zn2+') / (2 * FARADAY), VFe = (mu('Fe3+') - mu('Fe2+')) / FARADAY + phi, L = 44e-6;
  for (const kin of [{ alpha: 0.5 }, { lambda: 47.3e3 }]) {
    const s = new Device(build({
      library: [lib, metal('Zn')],
      stack: [ohmic(0, ['e-']), layer('Zn', 1e-6), { reactions: [{ equation: 'Zn2+ + 2 e- = Zn(s)', fixed: { 'Zn(s)': 0 }, k0: 4.8e-4, alpha: 0.48 }] }, layer('water', L, { name: 'w', c0: c }), bath(c, phi)],
      ports: [{ name: 'porous', region: 'w', from: 0, to: L, terminal: 'e-', area: 6.2e3, I: 0, reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 3.2e-6, ...kin }] }],
      grid: { hmin: L / 500, hmax: L / 30, ratio: 1.15 },
    })).solve();
    assert.ok(s.converged, `${JSON.stringify(kin)}: ${s.warnings[0]}`);
    assert.ok(Math.abs(s.ports[0].V - VFe) < 1e-9, `${JSON.stringify(kin)}: the port at ${s.ports[0].V} V, its couple's level ${VFe}`);
  }
});

test("a floating electrode's double layer stays out of a law its terminal passes none of (the total iron, in its circuit row)", () => {
  // (Ag/AgCl | Fe³⁺, Fe²⁺, KClO₄, KCl with a porous Fe³⁺/Fe²⁺ port, a double layer, at I = 0 | Pt
  // at the couple's level, all at rest. The total iron took the port's circuit row, and the
  // port's charge went in with it at λ·z, λ undefined: NaN, and every steady solve failed with a
  // non-finite update, at rest from the start. Without the double layer, or held, it converged.
  // The echem combos' case 65.)
  const c = { 'Fe3+': 4.27, 'Fe2+': 2.14, 'K+': 160.5, 'Cl-': 20.8 };
  c['ClO4-'] = 160.5 + 3 * 4.27 + 2 * 2.14 - 20.8;
  const lib = aqueous(['Fe3+', 'Fe2+', 'K+', 'ClO4-', 'Cl-'], { epsr: 0 }), mu = levels(lib, c);
  lib.species.push({ name: 'e-', z: -1 });
  const phi = (mu('Cl-') + 109789) / FARADAY, VFe = (mu('Fe3+') - mu('Fe2+')) / FARADAY + phi, L = 200e-6;
  const def = build({
    library: [lib, metal('Pt'), metal('Ag')],
    stack: [
      ohmic(0, ['e-']), layer('Ag', 1e-6), { reactions: [{ equation: 'AgCl(s) + e- = Ag(s) + Cl-', fixed: { 'AgCl(s)': -109789, 'Ag(s)': 0 }, k0: 3.5e-5, alpha: 0.5 }] },
      layer('water', L, { name: 'w', c0: c }), { reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 3.6e-4, alpha: 0.5 }] }, layer('Pt', 1e-6), ohmic(VFe, ['e-']),
    ],
    ports: [{ name: 'porous', region: 'w', terminal: 'e-', area: 8.4e4, I: 0, capacitance: { C: 0.2, zeroCharge: VFe - phi }, reactions: [{ equation: 'Fe3+ + e- = Fe2+', k0: 1.2e-7, alpha: 0.43 }] }],
    grid: { hmin: 2.6e-7, hmax: 1.1e-5, ratio: 1.12 },
  });
  const d = new Device(def);
  for (const cs of d.solver.constraints) for (const p of cs.carriers) assert.ok(Number.isFinite(p.w), `a carrier at weight ${p.w} in ${cs.key}`);
  const s = d.solve();
  assert.ok(s.converged, s.warnings[0]);
  assert.ok(Math.abs(s.ports[0].V - VFe) < 1e-9, `the port at ${s.ports[0].V} V, its couple's level ${VFe}`);
});

test('a closed unlike cell with Stern layers, held at its open circuit behind a resistance, solves cold at rest', () => {
  // (Cu | CuCl₂ + ZnCl₂ + NaCl | Zn, closed, at its OCV, −1.1258 V, behind R: laid out with the
  // zinc at 0 V (level with the copper), the held stage solved from that layout and failed, for
  // R from 1e-6 to 1e-2 Ω·m², where the same cell at R = 0, or without the Stern layers,
  // converged cold. The echem combos' case 38.)
  const c = { 'Cu2+': 10, 'Zn2+': 1.5, 'Na+': 100, 'Cl-': 123 };
  const lib = aqueous(['Cu2+', 'Zn2+', 'Na+', 'Cl-'], { epsr: 0 }), mu = levels(lib, c);
  const ocv = (mu('Zn2+') - mu('Cu2+')) / (2 * FARADAY), phi = -mu('Cu2+') / (2 * FARADAY);
  for (const [R, off] of [[1e-4, 0], [1e-2, 0.23]]) {
    const s = new Device(build({
      library: [lib, metal('Cu'), metal('Zn')],
      stack: [
        ohmic(0, ['e-']), layer('Cu', 1e-6),
        { phi: { type: 'capacitive', C: 0.45 }, zeroCharge: -phi + off, reactions: [{ equation: 'Cu2+ + 2 e- = Cu(s)', fixed: { 'Cu(s)': 0 }, k0: 1e-3, alpha: 0.4 }] },
        layer('water', 12.6e-6, { c0: c }),
        { phi: { type: 'capacitive', C: 0.22 }, zeroCharge: ocv - phi, reactions: [{ equation: 'Zn2+ + 2 e- = Zn(s)', fixed: { 'Zn(s)': 0 }, k0: 7.6e-5, alpha: 0.49 }] },
        layer('Zn', 1e-6), ohmic({ V: ocv, R }, ['e-']),
      ],
      grid: { hmin: 1.35e-8, hmax: 1.3e-7, ratio: 1.1 },
    })).solve();
    assert.ok(s.converged, `R ${R}: ${s.warnings[0]}`);
    assert.ok(Math.abs(s.terminals.right.V - ocv) < 1e-9 && Math.abs(s.current) < 1e-12, `R ${R}: ${s.terminals.right.V} V (its OCV ${ocv}), ${s.current} A/m²`);
  }
});
