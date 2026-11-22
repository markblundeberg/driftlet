// Neon's fourth review (21 Nov 2026), its own pick, the kit's polarization(): the surface held in
// equilibrium at each V (R4 3.1, 3.2), and what fuzzing the fix found next to them.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Device, DeviceError, FARADAY } from '../src/index.js';
import { build, layer, ohmic, bath, metal, polarization } from '../src/kit.js';

const salt = { 'Na+': 100, 'Cl-': 100 };
// A Pt face in 100 mM NaCl with tethered couples Ox_k + Na⁺ + e⁻ = Red_k (Neon's R4 3.1 device,
// one couple or two), the metal at V.
const tethered = (V, couples) => ({
  species: [{ name: 'e-', z: -1, cRef: 1000 }, { name: 'Na+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
  materials: {
    Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } },
    aq: { epsr: 0, species: { 'Na+': { D: 1.33e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } },
  },
  regions: [{ name: 'm', material: 'Pt', length: 1e-7 }, { name: 'el', material: 'aq', length: 1e-6, c0: salt }],
  interfaces: [{
    phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0,
    surface: Object.fromEntries(couples.flatMap(({ mu0, theta0 }, k) => [[`Ox${k}`, { capacity: 1e-5, mu0, theta0 }], [`Red${k}`, { capacity: 1e-5, mu0, theta0 }]])),
    reactions: couples.map((_, k) => ({ equation: `Ox${k} + Na+ + e- = Red${k}`, k0: 1e-2, alpha: 0.5 })),
  }],
  contacts: { left: { V, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, right: { V: 0, bath: { c: salt } } },
  grid: { minCells: 8 },
});

// Each couple's coverages against the closed form, nothing passing: its total kept from sol, its
// ratio Nernst's, ln(θ_Red/θ_Ox) falling by (V − V₀)/V_T; and the bare fraction sol's, at every V.
const againstNernst = (d, sol, Vs, n) => {
  const curve = polarization(d, sol, { face: 0 }, Vs, { surface: 'equilibrium' });
  const f = sol.interfaces[0], V0 = sol.V['e-'][d.model.grid.regionEnd[0]], VT = d.model.RT / FARADAY;
  let worst = 0;
  Vs.forEach((V, j) => {
    for (let k = 0; k < n; k++) {
      const ox = f.coverage[`Ox${k}`], red = f.coverage[`Red${k}`], L = Math.log(red / ox) - (V - V0) / VT;
      worst = Math.max(worst, Math.abs(curve.coverage[`Ox${k}`][j] / ((ox + red) / (1 + Math.exp(L))) - 1), Math.abs(curve.coverage[`Red${k}`][j] / ((ox + red) / (1 + Math.exp(-L))) - 1));
    }
    worst = Math.max(worst, Math.abs(curve.bare[j] / f.bare - 1));
  });
  return worst;
};

// R4 3.1: a tethered couple 99.99% covered, solved at 0.1 V and swept to −0.5, threw "the
// surface's equilibrium at V = -0.5 didn't converge", and so did 1e-8 and 1e-15 bare. Newton
// reached the root, but there its step was the gradient's rounding (n·θ − n·θ_held, two numbers
// near 1) over a Hessian of 5e-5, 1.1e-11, above the 1e-12 it stopped at, and the line search on
// F's value couldn't move it: 100 iterations. Made right, at 1e-15 the free direction from
// Gram–Schmidt (n_Ox and n_Red an ulp apart) still let the bare fraction drift to half of sol's.
// Neon's fourth review; now within 2e-14 of the closed form, the bare fraction kept.
test('polarization(surface: equilibrium): a tethered couple 99.99%, 99.999999% and 99.9999999999999% covered keeps its total and follows Nernst over ±1 V (Neon R4 3.1)', () => {
  const RT = 8.314462618 * 298.15;
  for (const [mu0, bare] of [[0, 1e-4], [0, 1e-8], [-45 * RT, 1e-15]]) {
    const d = new Device(tethered(0.1, [{ mu0, theta0: (1 - bare) / 2 }])), sol = d.solve();
    assert.ok(sol.converged);
    const worst = againstNernst(d, sol, [-0.5, 0, 0.5, -1, 1], 1);
    assert.ok(worst < 1e-12, `bare ${bare}: ${worst.toExponential(2)} off the closed form`);
  }
});

// Next to R4 3.1, two tethered couples at a bare fraction near 1e-15, swept over ±2 V: before, it
// threw at −2 V ("didn't converge"). Fuzzing the fix found two more ways it broke there: with the
// Hessian's entries 0.24 and its least eigenvalue the bare fraction (along the couples' common
// total), eliminating it gave a zero pivot; and in the echelon basis the step along that total was
// the other direction's rounding, 2e-8, and never shrank. Now the step is from the Hessian's
// square root, by Householder, with the filled total a direction of its own; within 3e-13.
test('polarization(surface: equilibrium): two tethered couples, nearly full, keep their totals and the bare fraction over ±2 V', () => {
  const b = 2e-15, d = new Device(tethered(0, [{ mu0: 0, theta0: (0.7 * (1 - b)) / 2 }, { mu0: -10e3, theta0: (0.3 * (1 - b)) / 2 }])), sol = d.solve();
  assert.ok(sol.converged && sol.interfaces[0].bare < 1e-14, `bare ${sol.interfaces[0].bare}`);
  const worst = againstNernst(d, sol, [-2, -1.1, -0.3, 0.4, 1.2, 2], 2);
  assert.ok(worst < 1e-11, `${worst.toExponential(2)} off the closed form`);
});

// Volmer–Tafel hydrogen evolution in 100 mM HCl at −0.2 V (Neon's R4 3.2 device): Volmer fast
// (k0 1e2), Tafel slow (1e-3), listed either way round.
const her = (V, order, kV = 1e2, kT = 1e-3) => {
  const rxs = { volmer: { equation: 'H+ + e- = H', k0: kV, alpha: 0.5 }, tafel: { equation: '2 H = H2', k0: kT, alpha: 0.5, fixed: { H2: 0 } } };
  return new Device({
    species: [{ name: 'e-', z: -1, cRef: 1000 }, { name: 'H+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }],
    materials: {
      Pt: { conductor: { species: 'e-', conductivity: 9.4e6 } },
      aq: { epsr: 0, species: { 'H+': { D: 9.3e-9, mu0: 0 }, 'Cl-': { D: 2.03e-9, mu0: 0 } } },
    },
    regions: [{ name: 'm', material: 'Pt', length: 1e-7 }, { name: 'el', material: 'aq', length: 1e-6, c0: { 'H+': 100, 'Cl-': 100 } }],
    interfaces: [{ phi: { type: 'capacitive', C: 0.2 }, zeroCharge: 0, surface: { H: { capacity: 1e-5, mu0: 0, theta0: 0.1 } }, reactions: order.map((r) => rxs[r]) }],
    contacts: { left: { V, terminal: 'e-', phi: 'neutral', species: { 'e-': 'equilibrium' } }, right: { V: 0, bath: { c: { 'H+': 100, 'Cl-': 100 } } } },
    grid: { minCells: 8 },
  });
};

// R4 3.2: which step was held at equilibrium went by listing order (single-surface-species steps
// first, then the rest, each group as listed). Listed [Volmer, Tafel] the curve was the solve's
// −191.37 A/m² to 8e-6; listed [Tafel, Volmer], Tafel was held at equilibrium, θ_H stuck at ½,
// and the curve read −2.39e7 A/m², 1.2e5 times the solved current, unwarned (Neon's fourth
// review). Now the faster step, by exchange rate, is held, whichever is listed first.
test('polarization(surface: equilibrium): Volmer–Tafel holds the fast step at equilibrium whichever is listed first (Neon R4 3.2)', () => {
  const curves = [['volmer', 'tafel'], ['tafel', 'volmer']].map((order) => {
    const d = her(-0.2, order), sol = d.solve();
    assert.ok(sol.converged);
    const c = polarization(d, sol, { face: 0 }, [sol.V['e-'][d.model.grid.regionEnd[0]]], { surface: 'equilibrium' });
    const volmer = c.reactions.find((r) => r.equation === 'H+ + e- = H'), tafel = c.reactions.find((r) => r.equation === '2 H = H2');
    assert.equal(volmer.equilibrated, true, order.join());
    assert.equal(tafel.equilibrated, undefined, order.join());
    assert.ok(Math.abs(c.current[0] / sol.current - 1) < 1e-4, `${order}: ${c.current[0]} against the solve's ${sol.current}`);
    return c;
  });
  assert.equal(curves[1].current[0], curves[0].current[0]);
  assert.equal(curves[1].coverage.H[0], curves[0].coverage.H[0]);
});

// What Volmer before Tafel had hidden: a step held at equilibrium is a quasi-equilibrium only
// while it's fast for what it passes. With Volmer and Tafel both at k0 1e-3 (round 5's A8 device),
// at −0.1 V the curve read −20.8 A/m² against the solve's −8.0, unwarned, and at 0.1 V its
// coverage was 31% under the solve's (Volmer passing 0.44 of its exchange rate). With Volmer at
// 1e2 against 1e-3 it holds to −0.4 V (passing 0.05 of it, the bare fraction 5% off the solve's)
// and at −0.5 V would pass 0.33, the bare fraction a quarter off. Now refused past a tenth, naming
// the step.
test('polarization(surface: equilibrium) refuses a step too slow for what it passes, and says which', () => {
  const ions = {
    species: [{ name: 'H+', z: 1, cRef: 1000 }, { name: 'Cl-', z: -1, cRef: 1000 }, { name: 'H2', z: 0, cRef: 1000 }],
    materials: { el: { epsr: 0, species: { 'H+': { D: 9e-9, mu0: 0 }, 'Cl-': { D: 2e-9, mu0: 0 }, H2: { D: 5e-9, mu0: 0 } } } },
  };
  const a8 = (V) => new Device(build({ T: 298.15, library: [ions, metal('Pt')],
    stack: [ohmic(V, ['e-']), layer('Pt', 1e-6), { surface: { 'H(ads)': { mu0: 0, capacity: 1e-5, theta0: 0.1 } }, reactions: [{ equation: 'H+ + e- = H(ads)', k0: 1e-3, alpha: 0.5 }, { equation: '2 H(ads) = H2', k0: 1e-3, alpha: 0.5 }] }, layer('el', 10e-6), bath({ 'H+': 100, 'Cl-': 100, H2: 1 }, 'Cl-')] }));
  for (const V of [0.1, -0.1]) {
    const d = a8(V), sol = d.solve();
    assert.ok(sol.converged);
    assert.throws(() => polarization(d, sol, { face: 0 }, sol.V['e-'][0], { surface: 'equilibrium' }), (e) => e instanceof DeviceError && /isn't near equilibrium: H\+ \+ e- = H\(ads\), held at equilibrium as the fastest step .* would pass .* times its exchange rate/.test(e.message));
  }
  // (Volmer fast: the curve to −0.4 V from one solution, each point's bare fraction within 5% of a solve's there, and −0.5 V refused)
  const d = her(-0.2, ['volmer', 'tafel']), sol = d.solve(), shift = sol.V['e-'][d.model.grid.regionEnd[0]] + 0.2;
  const Vs = [0, -0.2, -0.4], curve = polarization(d, sol, { face: 0 }, Vs.map((V) => V + shift), { surface: 'equilibrium' });
  Vs.forEach((V, j) => {
    const at = her(V, ['volmer', 'tafel']).solve();
    assert.ok(Math.abs(curve.bare[j] / at.interfaces[0].bare - 1) < 0.05 && Math.abs(curve.current[j] / at.current - 1) < 2e-3, `V ${V}: bare ${curve.bare[j]} against ${at.interfaces[0].bare}, I ${curve.current[j]} against ${at.current}`);
  });
  assert.throws(() => polarization(d, sol, { face: 0 }, [-0.5 + shift], { surface: 'equilibrium' }), /isn't near equilibrium: H\+ \+ e- = H, .* would pass 0\.00200 mol\/\(m²·s\), 0\.33\d times its exchange rate/);
});
