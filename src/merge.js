// Deep merge of plain objects, as Device.set() applies a patch: plain objects merge key by key,
// and arrays, typed arrays and every other value are replaced wholesale, as is an object marked
// WHOLE (a contact the kit wrote: ohmic(), bath()), which says all of what it is.

// (a symbol, not enumerable, so a definition's JSON and a comparison of definitions don't see
// it; a spread drops it, and the copy merges again)
export const WHOLE = Symbol.for('driftlet.whole');
export const markWhole = (o) => Object.defineProperty(o, WHOLE, { value: true });

const isPlainObject = (v) => v !== null && typeof v === 'object' && !Array.isArray(v) && !ArrayBuffer.isView(v);

export function merge(base, patch) {
  if (!isPlainObject(base) || !isPlainObject(patch)) return patch;
  const out = { ...base };
  for (const [k, v] of Object.entries(patch)) out[k] = isPlainObject(v) && isPlainObject(base[k]) && !v[WHOLE] ? merge(base[k], v) : v;
  return out;
}
