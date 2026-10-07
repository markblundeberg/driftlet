import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, FARADAY, GAS_CONSTANT } from '../src/index.js';
import { hodgkinHuxley, injector, pulse } from '../src/kit.js';

// A propagating action potential. Along a squid giant axon (x runs along it), the axoplasm is a
// strictly neutral solution whose ions drift in the field along x, and the membrane is a port
// through the whole region: a capacitance (1 µF/cm², 2/a of it per volume of a radius-a axon) and
// Hodgkin and Huxley's channels as conductance links to the sea water outside, gated by m³h and n⁴.
// Nothing of cable theory is put in; with the outside held (no external resistance), it should
// come out. Hodgkin and Huxley's own case (1952): a = 238 µm, R_i = 35.4 Ω·cm, 18.5 °C, where they
// computed 18.8 m/s.

const T = 273.15 + 18.5, RT = GAS_CONSTANT * T, VT = RT / FARADAY, mV = 1e-3;
const a = 238e-6, am = 2 / a, Cm = 0.01, L = 0.06, h = 2e-4;
const z = { 'Na+': 1, 'K+': 1, 'Cl-': -1, 'A-': -1 };
// Their reversal potentials (modern convention: rest −65 mV), the leak's carried by Cl⁻.
const E = { 'Na+': 50 * mV, 'K+': -77 * mV, 'Cl-': -54.387 * mV };
const IN = { 'Na+': 50, 'K+': 400, 'Cl-': 52.5, 'A-': 397.5 };
const OUT = Object.fromEntries(['Na+', 'K+', 'Cl-'].map((s) => [s, IN[s] * Math.exp((z[s] * E[s]) / VT)]));
const G = { 'Na+': 1200, 'K+': 360, 'Cl-': 3 }; // S/m² of membrane (120, 36, 0.3 mS/cm²)
// The small ions share one D, set so the axoplasm's resistivity is theirs, 35.4 Ω·cm (the
// impermeant anions barely move).
const DA = 1e-11;
const D = (1 / 0.354 / ((FARADAY * FARADAY) / RT) - DA * IN['A-']) / (IN['Na+'] + IN['K+'] + IN['Cl-']);
const gates = hodgkinHuxley({ inside: IN, outside: { ...OUT, 'A-': 1 }, T }).gates;

function axon() {
  const species = (Dv) => ({ D: Dv, mu0: 0, cRef: 1000 });
  // Each ion's level outside, against φ there (the port's V, held at 0): its chemical part.
  const offset = (s) => (VT / z[s]) * Math.log(OUT[s] / 1000);
  const sealed = { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' };
  return new Device({
    T,
    species: Object.entries(z).map(([name, zz]) => ({ name, z: zz })),
    materials: { axoplasm: { epsr: 0, species: { 'Na+': species(D), 'K+': species(D), 'Cl-': species(D), 'A-': species(DA) } } },
    regions: [{ name: 'axon', material: 'axoplasm', length: L, c0: IN }],
    contacts: { left: sealed, right: sealed },
    ports: [
      {
        name: 'membrane', region: 'axon', V: 0, terminal: 'K+', area: am, capacitance: { C: Cm, zeroCharge: 0 }, gates,
        species: {
          'Na+': { type: 'conductance', G: G['Na+'] * am, offset: offset('Na+'), gates: { m: 3, h: 1 } },
          'K+': { type: 'conductance', G: G['K+'] * am, offset: offset('K+'), gates: { n: 4 } },
          'Cl-': { type: 'conductance', G: G['Cl-'] * am, offset: offset('Cl-') },
        },
      },
      // the stimulus: K⁺ into the first millimetre, 0.2 ms at 30 A/m² of cross-section
      injector({ name: 'stim', region: 'axon', from: 0, to: 1e-3, species: 'K+', I: pulse({ start: 30e-3, width: 0.2e-3, amplitude: 30 }) }),
    ],
    grid: { hmin: h, hmax: h },
  });
}

const rate = (r, V) => {
  const y = (V - r.midpoint) / r.scale;
  if (r.type === 'exp') return r.rate * Math.exp(y);
  if (r.type === 'sigmoid') return r.rate / (1 + Math.exp(-y));
  return (r.rate * y) / -Math.expm1(-y);
};
const inf = (g, V) => rate(g.alpha, V) / (rate(g.alpha, V) + rate(g.beta, V));
const Iion = (V, m, hh, n) => G['Na+'] * m ** 3 * hh * (V - E['Na+']) + G['K+'] * n ** 4 * (V - E['K+']) + G['Cl-'] * (V - E['Cl-']);

// The cable equation itself, C V_t = (σ/a_m) V_xx − I_ion + injected, by the method of lines
// (explicit, second order) on the same grid, sealed ends: the time its spike crosses 0 mV at x.
function cable(rest) {
  const sigma = ((FARADAY * FARADAY) / RT) * (D * (IN['Na+'] + IN['K+'] + IN['Cl-']) + DA * IN['A-']);
  const Dv = sigma / (am * Cm), Nx = Math.round(L / h) + 1, dt = (0.4 * h * h) / Dv;
  let y = [new Float64Array(Nx).fill(rest), ...['m', 'h', 'n'].map((k) => new Float64Array(Nx).fill(inf(gates[k], rest)))];
  const f = (t, [V, m, hh, n]) => {
    const out = [0, 0, 0, 0].map(() => new Float64Array(Nx));
    for (let i = 0; i < Nx; i++) {
      const l = V[i > 0 ? i - 1 : 1], r = V[i < Nx - 1 ? i + 1 : Nx - 2];
      const inj = i * h <= 1e-3 && t >= 0.1e-3 && t < 0.3e-3 ? 30 / 1e-3 / am : 0;
      out[0][i] = (Dv * (l - 2 * V[i] + r)) / (h * h) + (inj - Iion(V[i], m[i], hh[i], n[i])) / Cm;
      ['m', 'h', 'n'].forEach((k, j) => (out[1 + j][i] = rate(gates[k].alpha, V[i]) * (1 - [m, hh, n][j][i]) - rate(gates[k].beta, V[i]) * [m, hh, n][j][i]));
    }
    return out;
  };
  const at = { 0.02: NaN, 0.04: NaN };
  for (let t = 0; t < 4e-3; t += dt) {
    const k1 = f(t, y), y1 = y.map((v, j) => v.map((x, i) => x + dt * k1[j][i])), k2 = f(t + dt, y1);
    const next = y.map((v, j) => v.map((x, i) => x + (dt / 2) * (k1[j][i] + k2[j][i])));
    for (const x of [0.02, 0.04]) {
      const i = Math.round(x / h), v0 = y[0][i], v1 = next[0][i];
      if (Number.isNaN(at[x]) && v0 < 0 && v1 >= 0) at[x] = t + (dt * -v0) / (v1 - v0);
    }
    y = next;
  }
  return 0.02 / (at[0.04] - at[0.02]);
}

test('a propagating action potential: the cable equation, out of ions drifting along the axon, at Hodgkin and Huxley\'s 18.8 m/s', () => {
  const d = axon();
  // From the concentrations given, the membrane settles to rest, in a damped oscillation over
  // ~20 ms (no pump: it has no steady state).
  const rest = d.advance(30e-3, { tol: 1e-4 }).ports[0], j = rest.Vm.length >> 1;
  const Vr = rest.Vm[j];
  // Rest is where their ionic current vanishes, the gates at their limits: −65.0 mV.
  assert.ok(Math.abs(Iion(Vr, inf(gates.m, Vr), inf(gates.h, Vr), inf(gates.n, Vr))) < 1e-3 && Math.abs(Vr + 65 * mV) < 0.01 * mV, `rest ${Vr}`);
  for (const k of ["m", "h", "n"]) assert.ok(Math.abs(rest.gates[k][j] - inf(gates[k], Vr)) < 1e-4, `${k} at rest: ${rest.gates[k][j]} vs ${inf(gates[k], Vr)}`);
  // The spike crosses 0 mV at 2 cm and then at 4 cm.
  const sol = d.advance(34e-3, { tol: 1e-4, dtMax: 5e-6, probes: [{ x: 0.02, quantity: 'phi' }, { x: 0.04, quantity: 'phi' }] });
  assert.ok(sol.converged);
  const { t, probes } = sol.trace;
  const cross = (v) => {
    for (let i = 1; i < t.length; i++) if (v[i - 1] < 0 && v[i] >= 0) return t[i - 1] + ((t[i] - t[i - 1]) * -v[i - 1]) / (v[i] - v[i - 1]);
    return NaN;
  };
  const v = 0.02 / (cross(probes[1]) - cross(probes[0]));
  const ref = cable(Vr);
  assert.ok(Math.abs(v / ref - 1) < 2e-3, `${v} m/s against the cable equation's ${ref}`);
  assert.ok(Math.abs(v / 18.8 - 1) < 0.01, `${v} m/s against Hodgkin and Huxley's 18.8`);
  // It overshoots, and the ions that carry it barely change the axoplasm (well under 0.1%).
  assert.ok(Math.max(...probes[0]) > 20 * mV);
  const c = d.solution().c;
  for (const s of ['Na+', 'K+']) assert.ok(c[s].every((x) => Math.abs(x / IN[s] - 1) < 1e-3), s);
});

test('saltatory conduction: a myelinated axon\'s spike jumps node to node, as a compartmental cable says', () => {
  // A thin axon (a = 5 µm) myelinated: nodes of Ranvier 1 µm long, every millimetre, each with
  // Hodgkin and Huxley's channels at five times their density (nodes are packed with them), and
  // internodes whose membrane is myelin, a capacitance 200 times smaller and no channels. Each node
  // and internode is a region of the same axoplasm, with a port of its own.
  const a5 = 5e-6, am5 = 2 / a5, nodes = 12, Ln = 1e-6, Li = 1e-3, Cmy = Cm / 200, gN = 5;
  const offset = (s) => (VT / z[s]) * Math.log(OUT[s] / 1000);
  const species = (Dv) => ({ D: Dv, mu0: 0, cRef: 1000 });
  const sealed = { species: { 'Na+': 'blocked', 'K+': 'blocked', 'Cl-': 'blocked', 'A-': 'blocked' }, phi: 'neutral' };
  const regions = [], ports = [];
  for (let k = 0; k < nodes; k++) {
    regions.push({ name: `node${k}`, material: 'axoplasm', length: Ln, c0: IN });
    ports.push({
      name: `node${k}`, region: `node${k}`, V: 0, terminal: 'K+', area: am5, capacitance: { C: Cm, zeroCharge: 0 }, gates,
      species: {
        'Na+': { type: 'conductance', G: gN * G['Na+'] * am5, offset: offset('Na+'), gates: { m: 3, h: 1 } },
        'K+': { type: 'conductance', G: gN * G['K+'] * am5, offset: offset('K+'), gates: { n: 4 } },
        'Cl-': { type: 'conductance', G: G['Cl-'] * am5, offset: offset('Cl-') },
      },
    });
    if (k === nodes - 1) break;
    regions.push({ name: `internode${k}`, material: 'axoplasm', length: Li, c0: IN });
    ports.push({ name: `myelin${k}`, region: `internode${k}`, V: 0, area: am5, capacitance: { C: Cmy, zeroCharge: 0 } });
  }
  ports.push(injector({ name: 'stim', region: 'node0', from: 0, to: Ln, species: 'K+', I: pulse({ start: 30.1e-3, width: 0.1e-3, amplitude: 100 }) }));
  const d = new Device({
    T,
    species: Object.entries(z).map(([name, zz]) => ({ name, z: zz })),
    materials: { axoplasm: { epsr: 0, species: { 'Na+': species(D), 'K+': species(D), 'Cl-': species(D), 'A-': species(DA) } } },
    regions, contacts: { left: sealed, right: sealed }, ports,
    grid: { hmin: 0.25e-6, hmax: 25e-6, ratio: 1.3 },
  });
  const rest = d.advance(30e-3, { tol: 1e-4 }).ports[0].Vm[2];

  // The compartmental cable: each node isopotential, each internode 50 passive compartments, backward
  // Euler for V (tridiagonal), the gates exponentially (Rush–Larsen), from rest.
  const len = [], isNode = [], nodeAt = [];
  for (let k = 0; k < nodes; k++) {
    nodeAt.push(len.length);
    len.push(Ln);
    isNode.push(true);
    if (k < nodes - 1) for (let j = 0; j < 50; j++) (len.push(Li / 50), isNode.push(false));
  }
  const N = len.length, sigma = 1 / 0.354, dt = 1e-7;
  const cap = len.map((l, i) => (isNode[i] ? Cm : Cmy) * am5 * l), gAx = len.slice(1).map((l, i) => sigma / (0.5 * (l + len[i])));
  const Inode = (v, m, hh, n) => gN * (G['Na+'] * m ** 3 * hh * (v - E['Na+']) + G['K+'] * n ** 4 * (v - E['K+'])) + G['Cl-'] * (v - E['Cl-']);
  let Vr = -0.065;
  for (let it = 0; it < 60; it++) Vr -= Inode(Vr, inf(gates.m, Vr), inf(gates.h, Vr), inf(gates.n, Vr)) / ((Inode(Vr + 1e-7, inf(gates.m, Vr + 1e-7), inf(gates.h, Vr + 1e-7), inf(gates.n, Vr + 1e-7)) - Inode(Vr - 1e-7, inf(gates.m, Vr - 1e-7), inf(gates.h, Vr - 1e-7), inf(gates.n, Vr - 1e-7))) / 2e-7);
  assert.ok(Math.abs(rest - Vr) < 0.01 * mV, `rest ${rest} vs ${Vr}`);
  const V = new Float64Array(N).fill(Vr), x = ['m', 'h', 'n'].map((k) => new Float64Array(N).fill(inf(gates[k], Vr)));
  const lo = new Float64Array(N), di = new Float64Array(N), up = new Float64Array(N), rhs = new Float64Array(N);
  const arrive = new Float64Array(nodes).fill(NaN);
  for (let t = 0; t < 0.45e-3; t += dt) {
    for (let i = 0; i < N; i++) {
      if (!isNode[i]) continue;
      ['m', 'h', 'n'].forEach((k, j) => {
        const al = rate(gates[k].alpha, V[i]), be = rate(gates[k].beta, V[i]), xi = al / (al + be);
        x[j][i] = xi + (x[j][i] - xi) * Math.exp(-(al + be) * dt);
      });
    }
    for (let i = 0; i < N; i++) {
      const A = am5 * len[i];
      let g = 0, Ie = 0;
      if (isNode[i]) {
        const gNa = gN * G['Na+'] * x[0][i] ** 3 * x[1][i], gK = gN * G['K+'] * x[2][i] ** 4, gL = G['Cl-'];
        g = A * (gNa + gK + gL);
        Ie = A * (gNa * E['Na+'] + gK * E['K+'] + gL * E['Cl-']);
      }
      lo[i] = i > 0 ? -gAx[i - 1] : 0;
      up[i] = i < N - 1 ? -gAx[i] : 0;
      di[i] = cap[i] / dt + g - lo[i] - up[i];
      rhs[i] = (cap[i] / dt) * V[i] + Ie + (i === 0 && t >= 0.1e-3 && t < 0.2e-3 ? 100 : 0);
    }
    for (let i = 1; i < N; i++) {
      const w = lo[i] / di[i - 1];
      di[i] -= w * up[i - 1];
      rhs[i] -= w * rhs[i - 1];
    }
    const prev = nodeAt.map((i) => V[i]);
    V[N - 1] = rhs[N - 1] / di[N - 1];
    for (let i = N - 2; i >= 0; i--) V[i] = (rhs[i] - up[i] * V[i + 1]) / di[i];
    nodeAt.forEach((i, k) => {
      if (Number.isNaN(arrive[k]) && prev[k] < -0.03 && V[i] >= -0.03) arrive[k] = t + (dt * (-0.03 - prev[k])) / (V[i] - prev[k]);
    });
  }

  // driftlet, sampled every microsecond: when the spike crosses −30 mV at each node's middle.
  const sites = [2, 4, 6, 8, 10], trace = sites.map(() => []), ts = [];
  for (let k = 1; k <= 450; k++) {
    const s = d.advance(30e-3 + k * 1e-6, { tol: 1e-5 });
    ts.push(s.time - 30e-3);
    sites.forEach((n, j) => trace[j].push(s.phi[s.x.findIndex((xx) => xx >= n * (Ln + Li) + Ln / 2)]));
  }
  const at = (v) => {
    for (let i = 1; i < ts.length; i++) if (v[i - 1] < -0.03 && v[i] >= -0.03) return ts[i - 1] + ((ts[i] - ts[i - 1]) * (-0.03 - v[i - 1])) / (v[i] - v[i - 1]);
    return NaN;
  };
  sites.forEach((n, j) => assert.ok(Math.abs(at(trace[j]) - arrive[n]) < 1e-6, `node ${n}: ${at(trace[j])} s vs ${arrive[n]}`));
  // 8 mm in about 0.2 ms: some 40 m/s, against 2.7 m/s for the same axon bare (18.73 √(5/238)).
  const v = (8 * (Ln + Li)) / (at(trace[4]) - at(trace[0]));
  assert.ok(v > 10 * 18.73 * Math.sqrt(a5 / a), `${v} m/s`);
});

test("a membrane port from its bath and the kit's hodgkinHuxley({ area }): the same axon as the explicit one", () => {
  // The bath gives each linked ion its level outside, against φ there (the port's V): its chemical
  // part, (μ° + RT ln(c/c_ref))/(zF), here the offsets written out above.
  const hh = hodgkinHuxley({ T, area: am });
  const explicit = axon();
  const def = structuredClone(explicit.def);
  const { name, region, V, area: ar, capacitance } = def.ports[0];
  def.ports[0] = { name, region, V, area: ar, capacitance, gates: hh.gates, species: hh.species, bath: { c: { 'Na+': OUT['Na+'], 'K+': OUT['K+'], 'Cl-': OUT['Cl-'] } } };
  const viaBath = new Device(def);
  const links = (d) => d.model.ports[0].species.map((l) => (l.type === 'conductance' ? [l.G, l.offset, l.gates] : null));
  const [a, b] = [links(explicit), links(viaBath)];
  a.forEach((l, i) => l && assert.ok(Math.abs(l[0] - b[i][0]) <= 1e-12 * l[0] && Math.abs(l[1] - b[i][1]) < 1e-12 && JSON.stringify(l[2]) === JSON.stringify(b[i][2]), `link ${i}: ${l} vs ${b[i]}`));
  // and so through a spike, step for step
  for (const d of [explicit, viaBath]) d.advance(30e-3, { tol: 1e-4 });
  const sa = explicit.advance(31.5e-3, { tol: 1e-4 }), sb = viaBath.advance(31.5e-3, { tol: 1e-4 });
  assert.ok(sa.phi.every((v, g) => Math.abs(v - sb.phi[g]) < 1e-12), 'the two axons part');
  // A bath's ion must be in the region; and a charged link with neither an offset nor a bath entry is refused.
  assert.throws(() => new Device({ ...def, ports: [{ ...def.ports[0], bath: { c: { 'Na+': 1, 'Ca2+': 1 } } }, def.ports[1]] }), /bath\.c\.Ca2\+: unknown species/);
  assert.throws(() => new Device({ ...def, ports: [{ ...def.ports[0], bath: { c: { 'Na+': OUT['Na+'] } } }, def.ports[1]] }), /give 'K\+' in the bath/);
});
