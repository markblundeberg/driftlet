// A live device for sliders, for driftlet/kit. Changes arrive faster than a device can solve
// them, so they're merged while a solve runs and only the latest is solved next; each solve
// starts warm from the last; a change that fails from there is ramped to from the last good
// state, then tried from scratch; and if it still fails, the last good solution stays. Optionally the solving runs in a
// Web Worker, keeping the page responsive.

import { Device } from './index.js';
import { merge } from './merge.js';

const isPlainObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v) && !ArrayBuffer.isView(v);
const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

// A point on the way from definition a to b (0 ≤ s ≤ 1): every finite number that differs is
// interpolated, everything else is b's.
export function between(a, b, s) {
  if (typeof a === 'number' && typeof b === 'number' && Number.isFinite(a) && Number.isFinite(b)) return a + s * (b - a);
  if (Array.isArray(a) && Array.isArray(b) && a.length === b.length) return b.map((v, k) => between(a[k], v, s));
  if (isPlainObject(a) && isPlainObject(b)) {
    const out = {};
    for (const k of Object.keys(b)) out[k] = k in a ? between(a[k], b[k], s) : b[k];
    return out;
  }
  return b;
}

/**
 * What a live device does with each change, on whichever thread solves it: apply, solve warm,
 * ramp from the last good state if that fails, and fall back to the last good solution.
 */
export class Session {
  /** @param {object} def device definition */
  constructor(def, { maxRamp = 32 } = {}) {
    this.device = new Device(def);
    this.maxRamp = maxRamp;
    this.good = null; // { cp, sol }: the last converged state
  }

  /**
   * @param {object} patch merged into the definition, as by Device.set()
   * @returns {{ solution: object | null, info: { ms: number, ramp: number, failed: boolean, error?: string, warnings?: string[] } }}
   */
  update(patch) {
    const t0 = now();
    const dev = this.device;
    let target;
    try {
      target = merge(dev.def, patch);
      dev.set(patch);
    } catch (e) {
      // An invalid change leaves the device as it was.
      return { solution: this.good?.sol ?? null, info: { ms: now() - t0, ramp: 0, failed: true, error: e.message } };
    }
    let sol = dev.solve();
    let ramp = 0;
    // Ramp from the last good state in more and more steps.
    for (let n = 2; !sol.converged && this.good && n <= this.maxRamp; n *= 2) {
      ramp = n;
      dev._rollback(this.good.cp);
      const from = this.good.cp.def;
      for (let k = 1; k <= n; k++) {
        dev.set(between(from, target, k / n));
        sol = dev.solve();
        if (!sol.converged) break;
      }
    }
    if (!sol.converged) {
      // Last, from scratch (with the solver's own continuation in bias).
      const fresh = new Device(target);
      const cold = fresh.solve();
      if (cold.converged) {
        this.device = fresh;
        sol = cold;
      }
    }
    if (sol.converged) {
      this.good = { cp: this.device._checkpoint(), sol };
      return { solution: sol, info: { ms: now() - t0, ramp, failed: false } };
    }
    if (this.good) this.device._rollback(this.good.cp);
    return { solution: this.good?.sol ?? sol, info: { ms: now() - t0, ramp, failed: true, warnings: sol.warnings } };
  }
}

/**
 * A live device.
 * @typedef {object} Live
 * @property {(patch: object) => Promise<{ solution: object | null, info: object }>} set merge a change
 *   (as Device.set()) and solve it, after the solve under way; changes made meanwhile are merged
 *   and solved once, and all their promises resolve with that solution
 * @property {Promise<{ solution: object | null, info: object }>} ready the first solve
 * @property {object | null} solution the last good solution
 * @property {() => void} close stop the worker, if any
 */

/**
 * A device for live demos: `live.set(patch)` as often as a slider moves, and draw what comes
 * back (or in `onsolution`). Solves run one at a time, each warm from the last, letting the page
 * paint in between; changes made during a solve are merged and solved together.
 *
 * On failure, the solution is the last good one, with `info.failed` (and `info.error` for an
 * invalid change, `info.warnings` for one that didn't converge). `info.ramp` says how many steps
 * a ramp from the last good state took (0: none needed), and `info.ms` the time.
 *
 * With `worker: true`, solving happens in a Web Worker (module workers are needed, as in every
 * current browser). A function returning a Worker-like object can be given instead.
 * @param {object} def device definition
 * @param {{ worker?: boolean | (() => any), onsolution?: (solution: object, info: object) => void, maxRamp?: number }} [opts]
 * @returns {Live}
 */
export function live(def, { worker = false, onsolution, maxRamp = 32 } = {}) {
  const engine = worker ? workerEngine(def, worker, maxRamp) : localEngine(def, maxRamp);
  let pending = null, waiters = [], busy = false;
  const api = {
    solution: null,
    set(patch) {
      pending = pending ? merge(pending, patch) : patch;
      const p = new Promise((resolve) => waiters.push(resolve));
      if (!busy) pump();
      return p;
    },
    close: () => engine.close(),
  };
  function pump() {
    if (pending === null) return;
    busy = true;
    const patch = pending, ws = waiters;
    pending = null;
    waiters = [];
    engine.run(patch).then((res) => {
      if (res.solution && !res.info.failed) api.solution = res.solution;
      if (onsolution && res.solution) onsolution(res.solution, res.info);
      for (const w of ws) w(res);
      busy = false;
      pump();
    });
  }
  api.ready = api.set({});
  return api;
}

// Solving on this thread, one macrotask per solve so the page can paint and take input.
function localEngine(def, maxRamp) {
  let session = null, error = null;
  try {
    session = new Session(def, { maxRamp });
  } catch (e) {
    error = e;
  }
  return {
    run: (patch) =>
      new Promise((resolve) =>
        setTimeout(() => resolve(session ? session.update(patch) : { solution: null, info: { ms: 0, ramp: 0, failed: true, error: error.message } }), 0),
      ),
    close() {},
  };
}

// Solving in a worker running worker.js; replies come back in order.
function workerEngine(def, worker, maxRamp) {
  const w = typeof worker === 'function' ? worker() : startWorker();
  const queue = [];
  w.onmessage = (e) => queue.shift()(e.data);
  w.postMessage({ type: 'init', def, maxRamp });
  return {
    run(patch) {
      return new Promise((resolve) => {
        queue.push(resolve);
        w.postMessage({ type: 'update', patch });
      });
    },
    close: () => w.terminate(),
  };
}

function startWorker() {
  const url = new URL('./worker.js', import.meta.url);
  try {
    return new Worker(url, { type: 'module' });
  } catch {
    // Cross-origin (the library loaded from a CDN): a same-origin stub that imports it.
    const blob = new Blob([`import ${JSON.stringify(url.href)};`], { type: 'text/javascript' });
    return new Worker(URL.createObjectURL(blob), { type: 'module' });
  }
}
