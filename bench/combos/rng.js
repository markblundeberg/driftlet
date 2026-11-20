// A seeded generator, so a case is known by its family and index alone (as in bench/stress.js).
export function rng(seed) {
  let a = seed >>> 0;
  const next = () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  return {
    next,
    u: (a0, b0) => a0 + (b0 - a0) * next(),
    log: (a0, b0) => 10 ** (a0 + (b0 - a0) * next()), // between 10^a and 10^b
    pick: (xs) => xs[Math.floor(next() * xs.length)],
    chance: (p) => next() < p,
    int: (lo, hi) => lo + Math.floor(next() * (hi - lo + 1)), // lo..hi inclusive
  };
}
export const hash = (s) => [...s].reduce((h, c) => Math.imul(h ^ c.charCodeAt(0), 16777619) >>> 0, 2166136261);
export const caseRng = (family, k) => rng(hash(family) + k);
