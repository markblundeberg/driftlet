// A transient recorder, for driftlet/kit: advance() in pieces (a run to completion, or a little
// each animation frame), with the trace joined up across calls and whole solutions kept as
// frames, for scrubbing back through a run or plotting snapshots side by side.

import { DeviceError } from './errors.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

/**
 * Record a device's transient. Frames are whole solutions: at the start, then at each of
 * `times` (or every `every` seconds from the start), which `advance()` lands on exactly; with
 * neither, at the end of each `advance()` call. The trace (`t`, `current`, `voltage`, and
 * `probes` if given) runs from the start through every accepted step of every call.
 * @param {import('./index.js').Device} device
 * @param {{ every?: number, times?: number[], probes?: { x: number, species: string, quantity?: 'c' | 'V', region?: string | number }[] }} [opts]
 */
export function recorder(device, opts = {}) {
  const { every, times, probes } = opts;
  if (every !== undefined && times !== undefined) fail('recorder: give every or times, not both');
  if (every !== undefined && !(every > 0 && Number.isFinite(every))) fail('recorder: every must be a time > 0 (s)');
  if (times !== undefined && !(Array.isArray(times) && times.every(Number.isFinite))) fail('recorder: times must be an array of times (s)');
  const first = device.solution();
  const t0 = first.time;
  const schedule = times ? [...times].filter((t) => t > t0).sort((a, b) => a - b) : null;
  const scheduled = schedule !== null || every !== undefined;
  let next = every !== undefined ? 1 : 0; // the next frame due: schedule[next], or t0 + next·every
  const due = () => (schedule ? schedule[next] : every !== undefined ? t0 + next * every : undefined);

  const trace = { t: [first.time], current: [first.current], voltage: [first.terminals.right.V - first.terminals.left.V] };
  if (probes) trace.probes = device.solver._probes(probes).read().map((v) => [v]); // the start's readings
  const frames = [first];

  /**
   * Advance to tEnd (s), stopping at each frame time on the way to keep a frame. Takes
   * advance()'s options; with `budgetMs` it returns when that much wall time is spent, `done`
   * saying whether tEnd was reached. Returns the latest solution (its own `trace` is only the
   * last piece: the recorder's `trace` is the whole run).
   * @param {number} tEnd
   * @param {object} [aopts] as for Device.advance()
   */
  function advance(tEnd, aopts = {}) {
    const start = now(), budget = aopts.budgetMs ?? Infinity;
    let sol;
    for (;;) {
      const stop = due();
      const target = stop !== undefined && stop < tEnd ? stop : tEnd;
      const left = budget - (now() - start);
      sol = device.advance(target, { ...aopts, ...(probes ? { probes } : {}), budgetMs: left > 0 ? left : 0 });
      const tr = sol.trace;
      for (let k = 0; k < tr.t.length; k++) {
        trace.t.push(tr.t[k]);
        trace.current.push(tr.current[k]);
        trace.voltage.push(tr.voltage[k]);
        if (probes) tr.probes.forEach((p, j) => trace.probes[j].push(p[k]));
      }
      if (!sol.converged || !sol.done) break;
      if (target === stop) {
        frames.push(sol);
        next++;
      }
      if (target === tEnd) break;
      if (now() - start > budget) {
        sol.done = false;
        break;
      }
    }
    if (!scheduled && sol.converged && sol.time > frames.at(-1).time) frames.push(sol);
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
      if (frames[mid].time <= t) lo = mid;
      else hi = mid - 1;
    }
    return frames[lo];
  }

  return { advance, frame, frames, trace };
}
