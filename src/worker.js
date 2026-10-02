// The Web Worker behind live(def, { worker: true }): it holds the device and solves each change
// it's sent, replying in order with { solution, info }.

import { Session } from './live.js';

let session = null, error = null;
globalThis.addEventListener('message', (e) => {
  const m = e.data;
  if (m.type === 'init') {
    try {
      session = new Session(m.def, { maxRamp: m.maxRamp });
    } catch (err) {
      error = err.message;
    }
  } else if (m.type === 'update') {
    let reply;
    try {
      reply = session ? session.update(m.patches) : { solution: null, errors: [], info: { ms: 0, ramp: 0, failed: true, error } };
    } catch (err) {
      reply = { solution: null, errors: [], info: { ms: 0, ramp: 0, failed: true, error: err.message } };
    }
    globalThis.postMessage(reply);
  }
});
