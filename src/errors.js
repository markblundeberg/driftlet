// The error a definition's validation throws, naming the offending path.

export class DeviceError extends Error {
  constructor(message) {
    super(message);
    this.name = 'DeviceError';
  }
}
