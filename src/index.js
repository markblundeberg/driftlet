// driftlet: 1D drift–diffusion–reaction solver.

import { normalizeDevice } from './device.js';

export { DeviceError, normalizeDevice } from './device.js';
export { buildGrid, gradedCells } from './grid.js';
export { BlockTridiagonal } from './blockTridiagonal.js';
export { bernoulli, bernoulliDerivative } from './bernoulli.js';
export * from './constants.js';

export class Device {
  /** @param {object} def plain, serialisable device definition */
  constructor(def) {
    this.def = def;
    this.model = normalizeDevice(def);
  }

  get grid() {
    return this.model.grid;
  }
}
