// A preload that checks, as a solver runs, the invariants of src/laws.js and src/lawRows.js that
// are cheap to check: load it into anything (the tests, the law sweep, the fuzzers) with
//   node --import ./tools/check-hook.js …
// and a breach throws (an Error, not a SolverError, so nothing catches it as a failed solve).
//   - findLaws leaves its structure as it found it (laws.js, invariant 7).
//   - After each steady assembly, no two laws hold the same row, and each holds an active one
//     (laws.js 2, lawRows.js 1). A flat stretch's rows are its own, but for its first node's:
//     the pin of its own amount (one with `col`) takes that.
//   - A node-by-node law's row has no terminal's ∂/∂V in it, and a flat stretch's row only its
//     contact's (lawRows.js 1).
import { createHash } from 'node:crypto';

const { Solver } = await import(new URL('../src/solver.js', import.meta.url).href);

const fail = (what) => {
  throw new Error(`invariant: ${what}`);
};

// A hash of everything findLaws reads (typed arrays by their bytes, the rest as JSON).
const digest = (s) => {
  const h = createHash('sha1'), seen = new WeakSet();
  const walk = (v) => {
    if (ArrayBuffer.isView(v)) return h.update(new Uint8Array(v.buffer, v.byteOffset, v.byteLength));
    if (v && typeof v === 'object') {
      if (seen.has(v)) return h.update('↺');
      seen.add(v);
      for (const k of Object.keys(v).sort()) {
        if (k === 'conserved') continue; // (a stretch's, written by Solver._laws from what findLaws returns)
        h.update(k);
        walk(v[k]);
      }
      return;
    }
    h.update(String(v));
  };
  walk(s._structure());
  return h.digest('hex');
};

const laws = Solver.prototype._laws;
Solver.prototype._laws = function (...a) {
  const before = digest(this);
  const r = laws.apply(this, a);
  if (digest(this) !== before) fail('findLaws wrote its structure');
  return r;
};

const assemble = Solver.prototype.assemble;
Solver.prototype.assemble = function (dt) {
  assemble.call(this, dt);
  if (!(this.steady && dt === Infinity)) return;
  const { M, loc, rix: R } = this, held = new Map();
  const take = (o, who) => {
    if (loc[o] < 0) return;
    if (held.has(o)) fail(`row ${o} (block ${Math.floor(o / M)}, slot ${o % M}) taken by ${held.get(o)} and by ${who}`);
    held.set(o, who);
  };
  const firstOf = new Map(); // a flat stretch's first node's row → the stretch
  for (const st of this.flatStretches) {
    const first = this.blockOfNode[st.nodes[0]] * M + 1 + st.species;
    firstOf.set(first, st);
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) {
      const o = this.blockOfNode[g] * M + 1 + st.species;
      if (o !== first) take(o, `a flat stretch of species ${st.species}`);
    }
  }
  const pin = (p, who) => {
    if (loc[p.row] < 0) fail(`${who} pinned at an inactive row ${p.row}`);
    if (firstOf.has(p.row) && !p.col) fail(`${who} took the first row of a flat stretch it doesn't move (no col)`);
    take(p.row, who);
  };
  for (const [o, st] of firstOf) if (!this.constraints.some((cs) => cs.row === o && cs.col)) take(o, `a flat stretch of species ${st.species}`);
  this.constraints.forEach((cs, j) => cs.circuit < 0 && pin(cs, `constraint ${j}`));
  this.affinities.forEach((af, j) => pin(af, `affinity ${j}`));
  this.pinnedIslands.forEach((isl, j) => pin(isl, `island ${j}`));
  this.localConstraints.forEach((lc, j) => lc.nodes.forEach((g) => take(this.blockOfNode[g] * M + lc.row, `local constraint ${j}`)));
  // What a terminal's ∂/∂V may sit in.
  const own = new Map();
  for (const st of this.flatStretches) {
    const side = st.leftOpen ? 0 : st.rightOpen ? 1 : -1;
    for (let g = st.nodes[0]; g <= st.nodes[1]; g++) own.set(this.blockOfNode[g] * M + 1 + st.species, side);
  }
  for (const [o, who] of held) {
    if (!(own.has(o) || who.startsWith('local'))) continue;
    if (who.startsWith('constraint')) continue; // (a pin: its row is free, see lawRows.js)
    for (const k of this.floating) if (this.termB[k][R[o]] !== 0 && k !== own.get(o)) fail(`${who}'s row ${o} holds terminal ${k}'s ∂/∂V`);
  }
};
