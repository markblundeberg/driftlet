// Deep merge of plain objects, as Device.set() applies a patch: plain objects merge key by key,
// and arrays, typed arrays and every other value are replaced wholesale, as is an object marked
// WHOLE (a contact the kit wrote: ohmic(), bath()), which says all of what it is.

// (a symbol, not enumerable, so a definition's JSON and a comparison of definitions don't see
// it; a spread drops it, and the copy merges again)
export const WHOLE = Symbol.for('driftlet.whole');
export const markWhole = (o) => Object.defineProperty(o, WHOLE, { value: true });

const isPlainObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v) && !ArrayBuffer.isView(v);
// A key set as the object's own: assigned, '__proto__' (a species' name, a JSON key) set the
// prototype instead, and the key was gone from the copy (Astra's third pass).
const put = (o, k, v) => (k === '__proto__' ? Object.defineProperty(o, k, { value: v, enumerable: true, writable: true, configurable: true }) : (o[k] = v));

// The WHOLE marks in a value, as paths of keys, and the marks put back: a structured clone
// (postMessage to a worker) drops a symbol, and live's worker merged the kit's contact it should
// have replaced, a bath swapped from Ag/AgCl to bare Cl⁻ keeping its 1.14 V offset (a fresh-eyes
// sweep). So live's worker engine sends the paths beside what it posts.
export function wholePaths(v, at = [], out = []) {
  if (v === null || typeof v !== 'object' || ArrayBuffer.isView(v)) return out;
  if (v[WHOLE]) out.push(at);
  for (const [k, x] of Object.entries(v)) wholePaths(x, [...at, k], out);
  return out;
}
export function markPaths(v, paths = []) {
  for (const p of paths) {
    const o = p.reduce((o, k) => (o !== null && typeof o === 'object' ? o[k] : undefined), v);
    if (o !== null && typeof o === 'object' && !o[WHOLE]) markWhole(o);
  }
  return v;
}

export function merge(base, patch) {
  if (!isPlainObject(base) || !isPlainObject(patch)) return patch;
  const out = { ...base };
  for (const [k, v] of Object.entries(patch)) put(out, k, isPlainObject(v) && Object.hasOwn(base, k) && isPlainObject(base[k]) && !v[WHOLE] ? merge(base[k], v) : v);
  return out;
}

// Whether two definitions' values are the same: plain objects and arrays by their contents (typed
// arrays too), anything else (numbers, strings, functions) by identity.
export function deepEqual(x, y) {
  if (x === y) return true;
  if (typeof x !== 'object' || typeof y !== 'object' || x === null || y === null) return false;
  if (ArrayBuffer.isView(x) || ArrayBuffer.isView(y)) return ArrayBuffer.isView(x) && ArrayBuffer.isView(y) && x.length === y.length && x.every((v, k) => v === y[k]);
  if (Array.isArray(x) !== Array.isArray(y)) return false;
  if (Array.isArray(x) && x.length !== y.length) return false; // (holes and a trailing undefined count: Astra's second pass)
  const keys = new Set([...Object.keys(x), ...Object.keys(y)]);
  for (const k of keys) if (!deepEqual(x[k], y[k])) return false;
  return true;
}

// A definition's own copy: plain objects and arrays copied all the way down (a WHOLE mark kept),
// anything else (functions, typed arrays, numbers) as it is. Kept by reference, a caller who edited
// the object it had passed and set() it again was told nothing changed, and the edit leaked in at
// the next rebuild (Astra's second pass).
export function ownCopy(v) {
  if (Array.isArray(v)) return v.map(ownCopy);
  if (v === null || typeof v !== 'object' || ArrayBuffer.isView(v) || Object.getPrototypeOf(v) !== Object.prototype) return v;
  const out = {};
  for (const [k, x] of Object.entries(v)) put(out, k, ownCopy(x));
  return v[WHOLE] ? markWhole(out) : out;
}
