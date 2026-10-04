import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { build, layer, bath, aqueous, recorder, check, describe } from '../src/kit.js';

// A cell membrane as a face: a capacitor (1 µF/cm²) that K⁺, Na⁺ and Cl⁻ cross by electrodiffusion
// in a constant field (a permeability each, Goldman–Hodgkin–Katz), with the impermeant anions of
// the cytoplasm (A⁻) and the Na⁺/K⁺ pump as a saturating, ATP-driven reaction on the same face.
// Outside (left) | membrane | cytoplasm (right), at 37 °C, typical mammalian concentrations.

const T = 310.15, RT = GAS_CONSTANT * T, VT = RT / FARADAY;
const OUT = { 'Na+': 145, 'K+': 5, 'Cl-': 110, 'A-': 40 }, IN = { 'Na+': 12, 'K+': 140, 'Cl-': 10, 'A-': 142 };
const REST = { 'K+': 1e-8, 'Na+': 0.04e-8, 'Cl-': 0.45e-8 }; // Hodgkin and Katz's 1 : 0.04 : 0.45
const IONS = ['Na+', 'K+', 'Cl-'];

const water = (epsr) => {
  const w = aqueous(IONS, { epsr });
  w.species.push({ name: 'A-', z: -1, cRef: 1000 });
  w.materials.water.species['A-'] = { D: 0.5e-9, mu0: 0 };
  return w;
};
const membrane = (P, reactions = []) => ({
  phi: { type: 'capacitive', C: 0.01 },
  species: { ...Object.fromEntries(IONS.map((s) => [s, P[s] ? { type: 'permeability', P: P[s] } : 'blocked'])), 'A-': 'blocked' },
  reactions,
});
// Between two baths (the inside one floating: no current, as for a cell), or with a closed
// cytoplasm (a wall that nothing crosses), so the cell's own contents set its concentrations.
const cell = ({ P = REST, epsr = 0, inside = IN, closed = false, reactions, right } = {}) =>
  build({
    T,
    library: [water(epsr)],
    stack: [
      bath(OUT, 'Cl-'),
      layer('water', 1e-6, { name: 'outside', c0: OUT }),
      membrane(P, reactions),
      layer('water', closed ? 5e-6 : 1e-6, { name: 'cytoplasm', c0: inside }),
      closed ? { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' } : bath(inside, 'Cl-', right ?? { I: 0 }),
    ],
    grid: epsr ? { hmin: 0.05e-9, hmax: 50e-9 } : { hmin: 50e-9, hmax: 200e-9 },
  });
const Vm = (s) => s.phi.at(-1) - s.phi[0]; // φ_in − φ_out
const ghk = (P, o, i) => VT * Math.log((P['K+'] * o['K+'] + P['Na+'] * o['Na+'] + P['Cl-'] * i['Cl-']) / (P['K+'] * i['K+'] + P['Na+'] * i['Na+'] + P['Cl-'] * o['Cl-']));
// GHK's current into the cell (A/m²) at membrane voltage V, from the two solutions' concentrations.
const ghkCurrent = (P, o, i, V) =>
  IONS.reduce((sum, s) => {
    const z = s === 'Cl-' ? -1 : 1, u = (z * V) / VT;
    return sum + z * FARADAY * P[s] * (Math.abs(u) < 1e-9 ? o[s] - i[s] : (u * (o[s] - i[s] * Math.exp(u))) / Math.expm1(u));
  }, 0);

test('a membrane face: the resting potential is Goldman–Hodgkin–Katz, with the solutions strictly neutral or their diffuse layers resolved', () => {
  for (const P of [REST, { 'K+': 1e-8, 'Na+': 20e-8, 'Cl-': 0.45e-8 }, { 'K+': 1e-8, 'Na+': 0.003e-8, 'Cl-': 0.03e-8 }]) {
    for (const epsr of [0, 74]) {
      const dev = new Device(cell({ P, epsr })), s = dev.solve();
      assert.ok(s.converged);
      // To a few µV: the leaks polarise the 1 µm of unstirred water beside the membrane (most with
      // Na⁺ channels open, ~1e-4 of the Na⁺ there).
      const expected = ghk(P, OUT, IN);
      assert.ok(Math.abs(Vm(s) - expected) < 1e-5, `ε ${epsr}, P_Na ${P['Na+']}: ${Vm(s)} vs GHK ${expected}`);
      // No current, but each ion leaks: K⁺ out, Na⁺ in, and their currents cancel with Cl⁻'s.
      const N = s.interfaces[0].N;
      assert.ok(N['K+'] < 0 && N['Na+'] > 0);
      assert.ok(Math.abs(N['Na+'] + N['K+'] - N['Cl-']) < 1e-9 * Math.abs(N['K+']));
      // The membrane holds C·V, on the cytoplasm's side as on the outside's.
      const g = s.x.findIndex((_, k) => s.region[k] === 0 && s.region[k + 1] === 1);
      assert.ok(Math.abs(s.interfaces[0].D / (-0.01 * (s.phi[g + 1] - s.phi[g])) - 1) < 1e-9);
      assert.ok(check(dev, s).ok);
    }
  }
});

test('a membrane face: the current–voltage curve is GHK at every voltage, and so is each ion\'s flux at the face', () => {
  const dev = new Device(cell({ right: 0 }));
  // The cytoplasm's bath held at a voltage: its Cl⁻ level, so Δφ = V + V_T ln(Cl_in/Cl_out).
  for (const V of [-0.15, -0.09, -0.06, 0, 0.03, 0.08]) {
    dev.set({ contacts: { right: { V: V - VT * Math.log(IN['Cl-'] / OUT['Cl-']) } } });
    const s = dev.solve(), dphi = Vm(s);
    assert.ok(Math.abs(dphi - V) < 1e-4, `Δφ ${dphi} at ${V}`);
    // Against the solutions' own concentrations, to the shift the membrane's charge makes beside
    // it (its box at the face holds C·V: here a few tenths of a percent at −150 mV).
    const I = s.current; // into the cell: toward +x
    const expected = ghkCurrent(REST, OUT, IN, dphi);
    assert.ok(Math.abs(I - expected) < 3e-3 * Math.abs(expected) + 1e-9, `${V} V: ${I} vs GHK ${expected} A/m²`);
    // Each flux is the law at the face's own state: its edge concentrations and its φ jump.
    const itf = s.interfaces[0], gL = s.x.findIndex((x, g) => s.region[g] === 0 && s.region[g + 1] === 1), gR = gL + 1;
    const u = (s.phi[gR] - s.phi[gL]) / VT;
    for (const sp of IONS) {
      const z = sp === 'Cl-' ? -1 : 1, cL = s.c[sp][gL], cR = s.c[sp][gR];
      const law = Math.abs(u) < 1e-9 ? REST[sp] * (cL - cR) : (REST[sp] * z * u * (cL - cR * Math.exp(z * u))) / (Math.exp(z * u) - 1);
      assert.ok(Math.abs(itf.N[sp] - law) < 1e-10 * Math.abs(law), `${sp} at ${V} V: ${itf.N[sp]} vs ${law}`);
    }
  }
});

test('a membrane face is the thin, ion-poor limit of a resolved lipid layer (the constant field)', () => {
  // The same membrane as a 5 nm region of lipid (εr 5) that ions dissolve in only sparingly
  // (25 kJ/mol above water), with D chosen for the same permeability, P = D·β/L.
  const L = 5e-9, dG = 25e3, beta = Math.exp(-dG / RT);
  for (const P of [REST, { 'K+': 1e-8, 'Na+': 20e-8, 'Cl-': 0.45e-8 }]) {
    const w = water(74);
    w.materials.lipid = { epsr: 5, species: Object.fromEntries(IONS.map((s) => [s, { D: (P[s] * L) / beta, mu0: w.materials.water.species[s].mu0 + dG }])) };
    const lipid = new Device(
      build({
        T,
        library: [w],
        stack: [bath(OUT, 'Cl-'), layer('water', 20e-9, { c0: OUT }), { dipole: 0 }, layer('lipid', L, { c0: Object.fromEntries(IONS.map((s) => [s, beta * Math.sqrt(OUT[s] * IN[s])])) }), { dipole: 0 }, layer('water', 20e-9, { c0: IN }), bath(IN, 'Cl-', { I: 0 })],
        grid: { hmin: 0.02e-9, hmax: 1e-9 },
      }),
    ).solve();
    const face = new Device(cell({ P, epsr: 74 })).solve();
    assert.ok(Math.abs(Vm(lipid) - Vm(face)) < 2e-5, `${Vm(lipid)} vs ${Vm(face)}`);
  }
});

test('opening K⁺ channels: the membrane charges at the GHK current over C, and a trace of φ starts at zero', () => {
  // Both solutions as laid down (no voltage); then the membrane potential falls toward E_K.
  const dev = new Device(cell({ P: { 'K+': 1e-8 } }));
  const rec = recorder(dev, { probes: [{ x: 0, quantity: 'phi' }, { x: 2e-6, quantity: 'phi' }] });
  rec.advance(2e-3, { tol: 1e-5 });
  const { t, voltage, probes: [out, inside] } = rec.trace;
  const dphi = t.map((_, k) => inside[k] - out[k]), offset = VT * Math.log(IN['Cl-'] / OUT['Cl-']);
  // The terminal reads the cytoplasm's Cl⁻ level against the outside's: Δφ less Cl⁻'s Nernst term, from the start.
  t.forEach((_, k) => assert.ok(Math.abs(voltage[k] + offset - dphi[k]) < 1e-9, `t = ${t[k]}: ${voltage[k] + offset} vs ${dphi[k]}`));
  assert.equal(dphi[0], 0);
  // dV/dt = I/C at first, I = F·P_K·(c_out − c_in) at V = 0.
  const k = t.findIndex((x) => x >= 2e-6), slope = dphi[k] / t[k], expected = (FARADAY * 1e-8 * (OUT['K+'] - IN['K+'])) / 0.01;
  assert.ok(Math.abs(slope / expected - 1) < 0.01, `${slope} vs ${expected} V/s`);
  assert.ok(dphi.at(-1) < -0.015 && dphi.at(-1) > VT * Math.log(OUT['K+'] / IN['K+']));
});

const pump = (dG = 50e3, vmax = 1e-7) => ({ equation: '3 Na+(right) + 2 K+(left) + ATP = 3 Na+(left) + 2 K+(right) + ADP', fixed: { ATP: dG, ADP: 0 }, vmax, K: { 'Na+': 10, 'K+': 1.5 } });
const inside = (s) => Object.fromEntries([...IONS, 'A-'].map((sp) => [sp, s.c[sp].at(-1)]));

test('the Na⁺/K⁺ pump in a closed cell: the pump–leak steady state is Mullins–Noda, and the pump runs at its saturating law', () => {
  // Steady state with no net flux of anything: Na⁺ leaks in as fast as the pump puts it out (3
  // per cycle), K⁺ leaks out as the pump brings it in (2 per cycle), so the leaks are in the
  // ratio 3:2 whatever the pump's kinetics, and V_m = V_T ln((r P_K K_o + P_Na Na_o)/(r P_K K_i + P_Na Na_i)), r = 3/2.
  for (const vmax of [1e-7, 1e-6]) {
    const dev = new Device(cell({ closed: true, reactions: [pump(50e3, vmax)] })), s = dev.solve();
    assert.ok(s.converged);
    const i = inside(s), r = 1.5, P = REST;
    const mn = VT * Math.log((r * P['K+'] * OUT['K+'] + P['Na+'] * OUT['Na+']) / (r * P['K+'] * i['K+'] + P['Na+'] * i['Na+']));
    assert.ok(Math.abs(Vm(s) - mn) < 1e-9, `${vmax}: V_m ${Vm(s)} vs Mullins–Noda ${mn}`);
    // Cl⁻ has no pump, so it sits at its Nernst potential; the impermeant anions stay.
    assert.ok(Math.abs(Vm(s) + VT * Math.log(OUT['Cl-'] / i['Cl-'])) < 1e-9);
    assert.ok(Math.abs(i['A-'] / IN['A-'] - 1) < 1e-3);
    // The rate at the face's state: vmax (Na_in/(Na_in + K_Na))³ (K_out/(K_out + K_K))² (1 − e^{−A/RT}).
    const g = s.x.findIndex((_, k) => s.region[k] === 0 && s.region[k + 1] === 1);
    const nai = s.c['Na+'][g + 1], ko = s.c['K+'][g];
    const A = 50e3 + 3 * (s.mu['Na+'][g + 1] - s.mu['Na+'][g]) + 2 * (s.mu['K+'][g] - s.mu['K+'][g + 1]);
    const law = vmax * (nai / (nai + 10)) ** 3 * (ko / (ko + 1.5)) ** 2 * -Math.expm1(-A / RT);
    assert.ok(Math.abs(s.interfaces[0].rates[0] / law - 1) < 1e-9, `${s.interfaces[0].rates[0]} vs ${law}`);
    // check() keeps a ledger per compartment: in the cytoplasm, Na⁺ leaks in across the membrane
    // exactly as fast as the pump puts it out, three per cycle.
    const report = check(dev, s), { ledgers } = report.items.find((it) => it.name === 'balance').details;
    assert.ok(report.ok, report.text);
    const na = ledgers.find((l) => l.species === 'Na+' && l.compartment === 'cytoplasm');
    const leak = na.terms.find((t) => /across face 0 \(permeability\)/.test(t.what)).rate, pumped = na.terms.find((t) => /ATP/.test(t.what)).rate;
    assert.ok(leak > 0 && Math.abs(pumped / leak + 1) < 1e-12 && Math.abs(-pumped / (3 * s.interfaces[0].rates[0]) - 1) < 1e-12);
    assert.match(report.text, /Na\+ in cytoplasm: /);
  }
  // A stronger pump keeps the cell further from equilibrium: more K⁺ in, less Na⁺, more negative.
  const weak = new Device(cell({ closed: true, reactions: [pump(50e3, 1e-7)] })).solve(), strong = new Device(cell({ closed: true, reactions: [pump(50e3, 1e-6)] })).solve();
  assert.ok(inside(strong)['K+'] > inside(weak)['K+'] && Vm(strong) < Vm(weak));
});

test('a pump with nothing leaking back stalls at its static head (3Δμ̄_Na − 2Δμ̄_K = ΔG_ATP); with the pump off, the cell sits at Donnan equilibrium', () => {
  for (const dG of [40e3, 50e3]) {
    const s = new Device(cell({ closed: true, P: {}, reactions: [pump(dG)] })).solve(), i = inside(s);
    const dmu = (sp) => RT * Math.log(i[sp] / OUT[sp]) + FARADAY * Vm(s); // μ̄_in − μ̄_out
    assert.ok(Math.abs((3 * -dmu('Na+') + 2 * dmu('K+')) / dG - 1) < 1e-9);
    assert.ok(Math.abs(s.interfaces[0].rates[0]) < 1e-20);
  }
  // Pump off: every permeant ion at equilibrium across the membrane, so K⁺, Na⁺ and Cl⁻ share one
  // Donnan ratio, set by the cytoplasm's impermeant anions, and V_m = −V_T ln r.
  const s = new Device(cell({ closed: true })).solve(), i = inside(s);
  const r = i['K+'] / OUT['K+'];
  assert.ok(r > 1);
  assert.ok(Math.abs(i['Na+'] / OUT['Na+'] / r - 1) < 1e-9 && Math.abs(OUT['Cl-'] / i['Cl-'] / r - 1) < 1e-9);
  assert.ok(Math.abs(Vm(s) + VT * Math.log(r)) < 1e-9);
});

test('membrane definitions are checked', () => {
  const def = (link, rx) => {
    const d = cell();
    d.interfaces[0].species['K+'] = link;
    if (rx) d.interfaces[0].reactions = [rx];
    return () => new Device(d);
  };
  const err = (re) => (e) => e instanceof DeviceError && re.test(e.message);
  assert.throws(def({ type: 'permeability' }), err(/species\.K\+\.P \(m\/s\)/));
  assert.throws(def({ type: 'permeability', P: 1e-8, G: 1 }), err(/\bG\b/));
  assert.throws(def({ type: 'permeability', P: 1e-8 }, { ...pump(), K: { 'Na+': 10 } }), err(/K\.K\+: give a half-saturation concentration for every species the forward reaction consumes/));
  assert.throws(def({ type: 'permeability', P: 1e-8 }, { ...pump(), K: { 'Na+': 10, 'K+': 1, 'Cl-': 1 } }), err(/K\.Cl-: not a species the forward reaction consumes/));
  assert.throws(def({ type: 'permeability', P: 1e-8 }, { ...pump(), k0: 1 }), err(/vmax and K \(saturating\), srh, or k0 and alpha/));
  // describe() lists what crosses by a law, and the kinetics.
  const d = cell({ reactions: [pump()] });
  assert.match(describe(d), /outside \| cytoplasm: capacitive, C = 0\.01 F\/m², dipole 0 V; crossing by a law: Na\+ 4\.00e-10 m\/s, K\+ 1\.00e-8 m\/s, Cl- 4\.50e-9 m\/s; blocked: A-; reactions: .* \(saturating, vmax 1\.00e-7 mol\/\(m²·s\)\)/);
});
