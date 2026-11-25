import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Device } from '../src/index.js';
import { aqueous, bath, build, half, layer, metal, ohmic } from '../src/kit.js';
import { LEARNED } from '../src/solver.js';

const pastRange = () => new Device(JSON.parse(readFileSync(new URL('./fixtures/debye-huckel-past-range.json', import.meta.url), 'utf8')));

// A passive Fe³⁺/Fe²⁺ electrode driven at a current whose steady state needs over 2 V.
const iron = () => new Device(build({
  library: [aqueous(['Fe3+', 'Fe2+', 'SO42-'], { epsr: 0 }), metal('Pt')],
  stack: [
    ohmic(0, ['e-']),
    layer('Pt', 1e-6),
    { reactions: [{ ...half('Fe3+ + e- = Fe2+'), k0: 1.4e-5, alpha: 0.33 }] },
    layer('water', 1.97e-5, { c0: { 'Fe3+': 27.66, 'Fe2+': 0.435, 'SO42-': 41.925 } }),
    bath({ 'Fe3+': 27.66, 'Fe2+': 0.435, 'SO42-': 41.925 }, 'SO42-', { I: 0.89 }),
  ],
  grid: { hmin: 1.3e-8, hmax: 4.2e-7, ratio: 1.17 },
}));

test('the Debye–Hückel solve stops on its wall-time budget before reporting a refusal', (t) => {
  const d = pastRange(), s = d.solver, newton = s.newton.bind(s);
  let elapsed = 0, rounds = 0, iterations = 0;
  t.mock.method(performance, 'now', () => elapsed);
  t.mock.method(s, 'newton', (...args) => {
    assert.ok(++rounds <= 1, 'the clock expired before the second Newton round');
    const result = newton(...args);
    iterations += result.iterations;
    elapsed = 201;
    return result;
  });
  const sol = d.solve({ budgetMs: 200 });
  assert.equal(rounds, 1);
  assert.equal(sol.iterations, iterations);
  assert.equal(sol.converged, false);
  assert.equal(sol.budgetSpent, true);
  assert.match(sol.warnings.join('\n'), /budget spent after 201 ms/);
  assert.doesNotMatch(sol.warnings.join('\n'), /refused by a statistics model/);
});

test('a zero steady budget allows one complete Newton round', (t) => {
  const d = pastRange(), s = d.solver, newton = s.newton.bind(s);
  let rounds = 0;
  t.mock.method(performance, 'now', () => 0);
  t.mock.method(s, 'newton', (...args) => {
    assert.ok(++rounds <= 1, 'a zero budget must stop before a retry even on an unchanged clock');
    return newton(...args);
  });
  const sol = d.solve({ budgetMs: 0 });
  assert.equal(rounds, 1);
  assert.ok(sol.iterations > 0);
  assert.equal(sol.converged, false);
  assert.equal(sol.budgetSpent, true);
});

test('an unlimited steady solve never reads the clock and matches a solve with time to finish', (t) => {
  let reads = 0;
  t.mock.method(performance, 'now', () => (reads++, 0));
  const plain = iron().solve();
  assert.ok(plain.converged);
  assert.equal(reads, 0);
  assert.deepEqual(iron().solve({ budgetMs: Infinity }), plain);
  assert.equal(reads, 0);
  assert.deepEqual(iron().solve({ budgetMs: 1 }), plain);
  assert.ok(reads > 0);
});

test('a steady budget spent during a held solve restores the state and its next answer', (t) => {
  const d = iron(), s = d.solver, newton = s.newton.bind(s);
  const steps = s._steadySteps.bind(s);
  const before = s._snapshot(), learned = LEARNED.map((k) => s[k]);
  const drives = s.terms.map((term) => term.drive);
  let elapsed = 0, first = true;
  t.mock.method(performance, 'now', () => elapsed);
  // Send the solve into continuation even if the floated attempt now converges outright.
  t.mock.method(s, '_steadySteps', (opts) => {
    const result = steps(opts);
    if (!first) return result;
    first = false;
    return { ...result, converged: false };
  });
  t.mock.method(s, 'newton', (...args) => {
    assert.equal(elapsed, 0, 'the expired call must start no more rounds');
    const result = newton(...args);
    if (s.terms[1].drive.kind === 'V') elapsed = 2;
    return result;
  });
  const stop = d.solve({ budgetMs: 1 });
  assert.equal(elapsed, 2);
  assert.equal(stop.converged, false);
  assert.equal(stop.budgetSpent, true);
  assert.deepEqual(s._snapshot(), before);
  assert.deepEqual(s.terms.map((term) => term.drive), drives);
  LEARNED.forEach((k, i) => assert.equal(s[k], learned[i], k));
  t.mock.restoreAll();
  const retry = d.solve(), fresh = iron().solve();
  assert.ok(fresh.converged);
  assert.deepEqual(retry, fresh);
});

test('a solve that converges inside its budget keeps its answer: the polish after it is not cut short', (t) => {
  // (A converged solve whose budget ran out in its last Newton round went on to the Kirchhoff
  // polish, which checked the clock and threw the converged answer away as budgetSpent. The polish
  // is one step, and only refines what's already an answer; the maintainer's review of Sol's budget,
  // 25 Nov.) The insertion electrode at rest of test/newton.test.js, which the polish takes from
  // Newton's tolerance to round-off (solved at 0.4 V, then at 0.3 V from there: from cold, either
  // lands exactly, with nothing to polish).
  const RT = 8.314462618 * 298.15, F = 96485.33212, x = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];
  const E = x.map((v) => 0.4 - (RT / F) * Math.log(v / (1 - v)));
  const host = () => new Device({
    species: [{ name: 'Li+', z: 1 }, { name: 'e-', z: -1 }],
    materials: { host: { epsr: 0, species: { 'Li+': { D: 1e-14, mu0: 0, cRef: 30000 }, 'e-': { D: 1e-4, mu0: 0, cRef: 30000 } }, statistics: [{ type: 'insertion', species: ['Li+', 'e-'], cMax: 30000, ocv: { x, E, muRef: 0 } }] } },
    regions: [{ material: 'host', length: 1e-6 }],
    contacts: { left: { V: 0, terminal: 'Li+', species: { 'Li+': 'equilibrium' }, phi: 'bulk' }, right: { V: 0.4, terminal: 'e-', species: { 'e-': 'equilibrium' }, phi: 'bulk' } },
  });
  const warm = () => {
    const dev = host();
    dev.solve();
    dev.set({ contacts: { right: { V: 0.3 } } });
    return dev;
  };
  const plain = warm().solve();
  const d = warm(), s = d.solver, newton = s.newton.bind(s);
  let elapsed = 0, polished = 0;
  const misses = s._missesKirchhoff.bind(s);
  t.mock.method(performance, 'now', () => elapsed);
  t.mock.method(s, 'newton', (...args) => {
    const result = newton(...args);
    if (result.converged) elapsed = 1e9; // (the budget runs out in the round that converges)
    return result;
  });
  t.mock.method(s, '_missesKirchhoff', () => { const m = misses(); polished += m ? 1 : 0; return m; });
  const sol = d.solve({ budgetMs: 1 });
  assert.equal(polished, 1);
  assert.ok(sol.converged && !sol.budgetSpent, sol.warnings.join('\n'));
  assert.ok(Math.abs(sol.terminals.right.current) < 1e-20, `${sol.terminals.right.current} A/m²`);
  assert.deepEqual(sol, plain);
});
