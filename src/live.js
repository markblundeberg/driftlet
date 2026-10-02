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
   * Apply changes in order, then solve once. Each change is applied on its own, so an invalid
   * one is refused (its message in `errors`) without losing the others.
   * @param {object[]} patches each merged into the definition, as by Device.set()
   * @returns {{ solution: object | null, errors: (string | undefined)[], info: { ms: number, ramp: number, failed: boolean, error?: string, warnings?: string[] } }}
   */
  update(patches) {
    const t0 = now();
    const errors = patches.map(() => undefined);
    let ramp = 0;
    const done = (sol, info) => ({ solution: sol, errors, info: { ms: now() - t0, ramp, ...info } });
    patches.forEach((patch, k) => {
      try {
        this.device.set(patch); // a change that fails leaves the device as it was
      } catch (e) {
        errors[k] = e.message;
      }
    });
    if (this.good && errors.every(Boolean)) return done(this.good.sol, { failed: true }); // nothing changed
    const target = this.device.def;
    try {
      let sol = this.device.solve();
      // Ramp from the last good state in more and more steps.
      for (let n = 2; !sol.converged && this.good && n <= this.maxRamp; n *= 2) {
        ramp = n;
        this.device._rollback(this.good.cp);
        const from = this.good.cp.def;
        for (let k = 1; k <= n; k++) {
          this.device.set(between(from, target, k / n));
          sol = this.device.solve();
          if (!sol.converged) break;
        }
      }
      if (!sol.converged) {
        // Last, from scratch (with the solver's own continuations).
        const fresh = new Device(target);
        const cold = fresh.solve();
        if (cold.converged) {
          this.device = fresh;
          sol = cold;
        }
      }
      if (sol.converged) {
        this.good = { cp: this.device._checkpoint(), sol };
        return done(sol, { failed: false });
      }
      if (this.good) this.device._rollback(this.good.cp);
      return done(this.good?.sol ?? sol, { failed: true, warnings: sol.warnings });
    } catch (e) {
      // Something the definition's checks can't see, found while building or solving.
      if (this.good) this.device._rollback(this.good.cp);
      return done(this.good?.sol ?? null, { failed: true, error: e.message });
    }
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
  let pending = [], busy = false;
  const api = {
    solution: null,
    set(patch) {
      const p = new Promise((resolve) => pending.push({ patch, resolve }));
      if (!busy) pump();
      return p;
    },
    close: () => engine.close(),
  };
  function pump() {
    if (pending.length === 0) return;
    busy = true;
    const batch = pending;
    pending = [];
    engine
      .run(batch.map((b) => b.patch))
      .catch((e) => ({ solution: api.solution, errors: [], info: { ms: 0, ramp: 0, failed: true, error: String(e?.message ?? e) } }))
      .then((res) => {
        if (res.solution && !res.info.failed) api.solution = res.solution;
        batch.forEach((b, k) => {
          const own = res.errors?.[k];
          b.resolve({ solution: res.solution, info: own ? { ...res.info, failed: true, error: own } : res.info });
        });
        busy = false;
        pump();
        if (onsolution && res.solution) {
          try {
            onsolution(res.solution, res.info);
          } catch (e) {
            // The page's own error: reported (as an uncaught error would be, in a browser) without
            // stopping the device.
            if (typeof globalThis.reportError === 'function') globalThis.reportError(e);
            else console.error(e);
          }
        }
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
    run: (patches) =>
      new Promise((resolve) =>
        setTimeout(() => resolve(session ? session.update(patches) : { solution: null, errors: [], info: { ms: 0, ramp: 0, failed: true, error: error.message } }), 0),
      ),
    close() {},
  };
}

// Solving in a worker running worker.js; replies come back in order.
function workerEngine(def, worker, maxRamp) {
  const w = typeof worker === 'function' ? worker() : startWorker();
  const queue = [];
  const failAll = (message) => {
    while (queue.length) queue.shift()({ solution: null, errors: [], info: { ms: 0, ramp: 0, failed: true, error: message } });
  };
  w.onmessage = (e) => queue.shift()?.(e.data);
  w.onerror = (e) => failAll(`the worker failed: ${e?.message ?? e}`);
  w.postMessage({ type: 'init', def, maxRamp });
  return {
    run(patches) {
      return new Promise((resolve) => {
        queue.push(resolve);
        try {
          w.postMessage({ type: 'update', patches });
        } catch (e) {
          // A change that can't be posted (a function, such as custom statistics).
          queue.pop();
          resolve({ solution: null, errors: patches.map(() => `can't be sent to a worker: ${e.message}`), info: { ms: 0, ramp: 0, failed: true } });
        }
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
