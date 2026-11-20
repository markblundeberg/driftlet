// A transient recorder, for driftlet/kit: advance() in pieces (a run to completion, or a little
// each animation frame), with the trace joined up across calls and whole solutions kept as
// frames, for scrubbing back through a run or plotting snapshots side by side.

import { DeviceError, checkOptions } from './errors.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

/**
 * Record a device's transient. Frames are whole solutions: at the start, then at each of
 * `times` (or every `every` seconds from the start), which `advance()` lands on exactly; with
 * neither, at the end of each `advance()` call. The trace (`t`, `current`, `voltage`,
 * `terminals` with each one's V, current and charge passed, and `probes` if given) runs from the start through every accepted step of every call.
 * @param {import('./index.js').Device} device
 * @param {{ every?: number, times?: number[], probes?: ({ x: number, species?: string, quantity?: 'c' | 'V' | 'phi', region?: string | number } | { interface: number, species?: string, gate?: string })[] }} [opts]
 */
export function recorder(device, opts = {}) {
  checkOptions('recorder', arguments[1], ['every', 'times', 'probes']);
  const { every, times, probes } = opts;
  if (every !== undefined && times !== undefined) fail('recorder: give every or times, not both');
  if (every !== undefined && !(every > 0 && Number.isFinite(every))) fail('recorder: every must be a time > 0 (s)');
  if (times !== undefined && !(Array.isArray(times) && times.every(Number.isFinite))) fail('recorder: times must be an array of times (s)');
  const first = device.solution();
  // (and each terminal's charge passed, joined up as advance()'s traces give it, from the running
  // total the device carries: summed back from the currents, an RC's charge came out 22% short;
  // the pieces' charges were in the frames' own traces only; Sol's third review)
  const startCharge = (nm) => device.solver.termCharge[device.solver.terms.findIndex((t) => t.name === nm)];
  const t0 = first.time;
  const schedule = times ? [...times].filter((t) => t > t0).sort((a, b) => a - b) : null;
  const scheduled = schedule !== null || every !== undefined;
  let next = every !== undefined ? 1 : 0; // the next frame due: schedule[next], or t0 + next·every
  const due = () => (schedule ? schedule[next] : every !== undefined ? t0 + next * every : undefined);
  // Two times the same but for round-off: a frame due at t0 + 3·0.1 = 0.30000000000000004 is the
  // one an advance(0.3) ends on. (Taken as later, it wasn't kept, and frame(0.3) gave the state
  // at 0.2; Astra's fifth pass.)
  const same = (a, b) => Math.abs(a - b) <= 4 * Number.EPSILON * Math.max(Math.abs(a), Math.abs(b));

  const trace = { t: [first.time], current: [first.current], voltage: [first.terminals.right.V - first.terminals.left.V],
    terminals: Object.fromEntries(Object.entries(first.terminals).map(([nm, r]) => [nm, { V: [r.V], current: [r.current], charge: [startCharge(nm)] }])) };
  if (probes) trace.probes = device.solver._probes(probes).read().map((v) => [v]); // the start's readings
  const frames = [first];

  /**
   * Advance to tEnd (s), stopping at each frame time on the way to keep a frame. Takes
   * advance()'s options; with `budgetMs` it returns when that much wall time is spent, `done`
   * saying whether tEnd was reached. Returns the latest solution (its own `trace` is only the
   * last piece: the recorder's `trace` is the whole run; its `steps` and `rejected` are the call's).
   * `dt0` is the call's first step and `maxSteps` the call's limit, as for one advance(); probes
   * are given to recorder().
   * @param {number} tEnd
   * @param {object} [aopts] as for Device.advance()
   */
  function advance(tEnd, aopts = {}) {
    // (refused: the pieces' probes went into no trace, and with probes given to recorder() too,
    // the solution's trace held those instead of the ones asked for; Astra's fifth pass)
    if (aopts?.probes !== undefined) fail('recorder: advance(): probes are given to recorder(device, { probes }), which reads them into its trace across calls; leave them out of advance()');
    // (refused before a step is taken: a port added since was left out of the trace, and one
    // removed threw a TypeError with the device already advanced, the trace ragged; Astra's fifth pass)
    const names = device.solver.terms.map((t) => t.name), kept = Object.keys(trace.terminals);
    if (names.length !== kept.length || !names.every((nm) => kept.includes(nm))) {
      fail(`recorder: advance(): the device's terminals changed since recorder() (recorded: ${kept.join(', ')}; now: ${names.join(', ')}), so its trace can't run on. Make a new recorder(device) for the device as it is now`);
    }
    // Frames due at or before the device's present time (advanced directly, say) are passed: not
    // kept as copies of now (frames at 0, 2.5, 2.5, 3 where 0, 1, 2, 3 were due; Astra's fifth pass).
    const tNow = device.solver.time;
    if (every !== undefined) next = Math.max(next, Math.floor((tNow - t0) / every));
    while (due() !== undefined && (due() <= tNow || same(due(), tNow))) next++;
    // The call's dt0 is its first step, and its maxSteps counts the steps of every piece (each
    // piece given them, maxSteps 1 ran a step per frame to the end, and dt0 restarted the steps
    // after every frame: 43 steps where 18 did; Astra's fifth pass).
    const { dt0, maxSteps, ...rest } = aopts;
    const start = now(), budget = aopts.budgetMs ?? Infinity;
    let sol, steps = 0, rejected = 0;
    for (;;) {
      const stop = due();
      const end = stop !== undefined && same(stop, tEnd); // (the frame due where the call ends)
      const target = stop !== undefined && stop < tEnd && !end ? stop : tEnd;
      const left = budget - (now() - start);
      const limit = maxSteps === undefined || !(maxSteps >= 0) ? maxSteps : Math.max(0, maxSteps - steps - rejected);
      sol = device.advance(target, {
        ...rest,
        ...(steps + rejected === 0 && dt0 !== undefined ? { dt0 } : {}),
        ...(maxSteps !== undefined ? { maxSteps: limit } : {}),
        ...(probes ? { probes } : {}),
        budgetMs: left > 0 ? left : 0,
      });
      steps += sol.steps ?? 0;
      rejected += sol.rejected ?? 0;
      const tr = sol.trace;
      for (let k = 0; k < tr.t.length; k++) {
        trace.t.push(tr.t[k]);
        trace.current.push(tr.current[k]);
        trace.voltage.push(tr.voltage[k]);
        for (const nm in trace.terminals) {
          trace.terminals[nm].V.push(tr.terminals[nm].V[k]);
          trace.terminals[nm].current.push(tr.terminals[nm].current[k]);
          trace.terminals[nm].charge.push(tr.terminals[nm].charge[k]);
        }
        if (probes) tr.probes.forEach((p, j) => trace.probes[j].push(p[k]));
      }
      if (!sol.converged || !sol.done) break;
      if (target === stop || end) {
        frames.push(sol);
        next++;
      }
      if (target === tEnd) break;
      if (now() - start > budget) {
        sol.done = false;
        break;
      }
      if (maxSteps !== undefined && steps + rejected >= maxSteps) {
        sol.done = false;
        break;
      }
    }
    if (maxSteps !== undefined && sol.converged && !sol.done && steps + rejected >= maxSteps) {
      // (said as one advance() says it, of the call: the piece's own said its share of the steps
      // and its own end)
      const said = `advance stopped at t = ${sol.time.toPrecision(6)} s of ${tEnd.toPrecision(6)}: maxSteps (${maxSteps}) used up (${rejected} of them rejected). Call again to go on, or loosen tol.`;
      sol.warnings = [...sol.warnings.filter((w) => !/^advance stopped at .*: maxSteps \(\d+\) used up/.test(w)), said];
    }
    if (!scheduled && sol.converged && sol.time > frames.at(-1).time) frames.push(sol);
    sol.steps = steps;
    sol.rejected = rejected;
    return sol;
  }

  /**
   * The frame at or just before time t (the first frame before the start).
   * @param {number} t
   */
  function frame(t) {
    let lo = 0, hi = frames.length - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (frames[mid].time <= t || same(frames[mid].time, t)) lo = mid; // (at, but for round-off: every 0.1's frame at 0.30000000000000004 is frame(0.3))
      else hi = mid - 1;
    }
    return frames[lo];
  }

  return { advance, frame, frames, trace };
}
