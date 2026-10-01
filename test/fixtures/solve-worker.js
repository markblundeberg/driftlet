// Solves a posted device definition off the main thread (see portability.test.js). Run on its
// own (node --test picks up every file under test/), it does nothing.
import { parentPort, workerData } from 'node:worker_threads';
import { Device } from '../../src/index.js';

if (parentPort) {
  const sol = new Device(workerData.def).solve();
  parentPort.postMessage({ converged: sol.converged, current: sol.current, phi: sol.phi });
}
