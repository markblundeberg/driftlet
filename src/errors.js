// The error a definition's validation throws, naming the offending path.

export class DeviceError extends Error {
  constructor(message) {
    super(message);
    this.name = 'DeviceError';
  }
}

// The nearest of `names` to a mistyped one: the same but for case, '_' or '-', or one or two
// letters off (dtmax for dtMax, engery for energy).
function nearest(name, names) {
  const plain = (s) => s.toLowerCase().replace(/[_-]/g, '');
  const edits = (a, b) => {
    let row = Array.from({ length: b.length + 1 }, (_, j) => j);
    for (let i = 1; i <= a.length; i++) {
      const next = [i];
      for (let j = 1; j <= b.length; j++) next[j] = Math.min(row[j] + 1, next[j - 1] + 1, row[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
      row = next;
    }
    return row[b.length];
  };
  return names.find((n) => plain(n) === plain(name)) ?? names.find((n) => n.length > 3 && edits(plain(n), plain(name)) <= 2);
}

/**
 * A function's options: a name it doesn't take (a typo) is an error, not silently ignored.
 * @param {string} where the function, for the message
 * @param {unknown} opts
 * @param {string[]} allowed
 * @param {new (message: string) => Error} [Kind]
 */
export function checkOptions(where, opts, allowed, Kind = DeviceError) {
  if (opts === undefined || opts === null) return;
  if (typeof opts !== 'object' || Array.isArray(opts)) throw new Kind(`${where}: options must be an object (${allowed.join(', ')})`);
  for (const k of Object.keys(opts)) {
    if (opts[k] === undefined || allowed.includes(k)) continue;
    const near = nearest(k, allowed);
    throw new Kind(`${where}: no option '${k}'${near ? ` (did you mean '${near}'?)` : ''}; the options are ${allowed.join(', ')}`);
  }
}

/**
 * A value as an error message shows it: a number as itself (NaN and Infinity too, which JSON
 * would print as null), a string quoted and called one ("350", a string, where a number was
 * wanted).
 * @param {unknown} v
 * @returns {string}
 */
export function shown(v) {
  if (typeof v === 'number' || typeof v === 'bigint' || v === undefined) return String(v);
  if (typeof v === 'string') return `${JSON.stringify(v)} (a string)`;
  if (typeof v === 'function') return 'a function';
  try {
    return JSON.stringify(v, (k, x) => (typeof x === 'number' && !Number.isFinite(x) ? String(x) : x)) ?? String(v);
  } catch {
    return String(v);
  }
}
