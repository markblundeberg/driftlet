// driftlet's Web Worker (src/worker.js) on a node worker thread, through the two globals a
// browser worker has (see live.test.js). Run on its own, it does nothing.
import { parentPort } from 'node:worker_threads';

if (parentPort) {
  globalThis.postMessage = (m) => parentPort.postMessage(m);
  globalThis.addEventListener = (type, fn) => parentPort.on(type, (data) => fn({ data }));
  await import('../../src/worker.js');
}
