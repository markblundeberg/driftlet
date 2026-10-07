// driftlet/kit: voltage-gated channels as plain definition data. Hodgkin and Huxley's squid axon:
// its gates, and permeabilities for its Na⁺, K⁺ and leak channels, for a membrane face.

import { DeviceError } from './errors.js';
import { FARADAY, GAS_CONSTANT } from './constants.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const need = (ok, message) => ok || fail(message);
const num = (v) => typeof v === 'number' && Number.isFinite(v);

// Hodgkin and Huxley (1952), J. Physiol. 117, 500, at their 6.3 °C, in the modern convention
// (V = φ_in − φ_out, rest at −65 mV), in SI: rates in 1/s, voltages in V.
const mV = 1e-3;
const HH_GATES = {
  m: { alpha: { type: 'expLinear', rate: 1000, midpoint: -40 * mV, scale: 10 * mV }, beta: { type: 'exp', rate: 4000, midpoint: -65 * mV, scale: -18 * mV } },
  h: { alpha: { type: 'exp', rate: 70, midpoint: -65 * mV, scale: -20 * mV }, beta: { type: 'sigmoid', rate: 1000, midpoint: -35 * mV, scale: 10 * mV } },
  n: { alpha: { type: 'expLinear', rate: 100, midpoint: -55 * mV, scale: 10 * mV }, beta: { type: 'exp', rate: 125, midpoint: -65 * mV, scale: -80 * mV } },
};
const HH_T = 273.15 + 6.3;

// The GHK current (outward, A/m²) of an ion at unit permeability, at voltage V (inside − outside).
export function ghkCurrent(z, inside, outside, V, T) {
  const u = (z * V * FARADAY) / (GAS_CONSTANT * T);
  if (Math.abs(u) < 1e-9) return z * FARADAY * (inside - outside);
  return (z * FARADAY * u * (inside - outside * Math.exp(-u))) / -Math.expm1(-u);
}

/**
 * Hodgkin and Huxley's squid-axon channels: for a membrane face with the outside on its left
 * (so the gates' voltage, φ_right − φ_left, is the membrane potential): the face's `gates`, and
 * `species` links for Na⁺ (P·m³h), K⁺ (P·n⁴) and the leak (ungated), to spread into the face
 * beside its φ law, `{ phi: { type: 'capacitive', C: 0.01 }, ...hh }`, adding links (or blocks)
 * for any other species.
 *
 * The currents are GHK's, which rectify, so a permeability matches a conductance at one voltage
 * only: each is the one whose chord conductance at `at` (rest, −65 mV) is Hodgkin and Huxley's
 * (`g`: 120, 36 and 0.3 mS/cm²), at the concentrations given. The leak is carried by one ion
 * (`leak`, of charge `leakZ`: Cl⁻ by default, whose Nernst level in squid is near their leak
 * reversal, −54 mV). The
 * rates are theirs at 6.3 °C, scaled to `T` by Q₁₀ = 3.
 *
 * Or, given `area` (membrane per volume of the region, 2/a for an axon of radius a), for a membrane
 * port along a region: links that are their linear conductances, G = g × area per volume, gated the
 * same, each ion's reversal potential following from the port's bath and the concentrations inside
 * (so `inside`, `outside` and `at` aren't needed).
 * @param {{ inside: object, outside: object, T?: number, g?: { Na?: number, K?: number, leak?: number },
 *   leak?: string, leakZ?: number, at?: number, names?: { Na?: string, K?: string }, area?: number }} opts concentrations in
 *   mol/m³ by species name, K, S/m², V
 * @returns {{ gates: object, species: object }}
 */
export function hodgkinHuxley({ inside, outside, T = HH_T, g = {}, leak = 'Cl-', leakZ = -1, at = -65 * mV, names = {}, area } = {}) {
  need(num(T) && T > 0, `hodgkinHuxley: T must be a temperature in K, got ${T}`);
  need(num(at), `hodgkinHuxley: at must be a voltage (V), got ${at}`);
  const G = { Na: 1200, K: 360, leak: 3, ...g };
  const Na = names.Na ?? 'Na+', K = names.K ?? 'K+';
  const q10 = Math.exp(Math.log(3) * ((T - HH_T) / 10)); // Q₁₀ = 3
  const scaled = (r) => ({ ...r, rate: r.rate * q10 });
  const gates = Object.fromEntries(Object.entries(HH_GATES).map(([k, gt]) => [k, { alpha: scaled(gt.alpha), beta: scaled(gt.beta) }]));
  if (area !== undefined) {
    need(num(area) && area > 0, `hodgkinHuxley: area must be the membrane per volume (m²/m³; 2/a for an axon of radius a), got ${area}`);
    for (const [k, v] of Object.entries(G)) need(num(v) && v > 0, `hodgkinHuxley: g.${k} must be a conductance > 0 (S/m²), got ${v}`);
    return {
      gates,
      species: {
        [Na]: { type: 'conductance', G: G.Na * area, gates: { m: 3, h: 1 } },
        [K]: { type: 'conductance', G: G.K * area, gates: { n: 4 } },
        [leak]: { type: 'conductance', G: G.leak * area },
      },
    };
  }
  need(inside && outside && typeof inside === 'object' && typeof outside === 'object', 'hodgkinHuxley: give inside and outside concentrations (mol/m³) by species name (or, for a membrane port, its area)');
  const P = (name, z, gval, what) => {
    need(num(inside[name]) && inside[name] > 0 && num(outside[name]) && outside[name] > 0, `hodgkinHuxley: give ${name}'s concentration inside and outside (mol/m³), for its ${what}`);
    need(num(gval) && gval > 0, `hodgkinHuxley: g.${what} must be a conductance > 0 (S/m²), got ${gval}`);
    const E = ((GAS_CONSTANT * T) / (z * FARADAY)) * Math.log(outside[name] / inside[name]);
    need(Math.abs(at - E) > 1e-4, `hodgkinHuxley: ${name}'s Nernst level, ${(E * 1e3).toFixed(1)} mV, is at the calibration voltage; choose another 'at'`);
    return (gval * (at - E)) / ghkCurrent(z, inside[name], outside[name], at, T);
  };
  need(Number.isInteger(leakZ) && leakZ !== 0, `hodgkinHuxley: leakZ must be the leak ion's charge number, got ${leakZ}`);
  return {
    gates,
    species: {
      [Na]: { type: 'permeability', P: P(Na, 1, G.Na, 'Na channels'), gates: { m: 3, h: 1 } },
      [K]: { type: 'permeability', P: P(K, 1, G.K, 'K channels'), gates: { n: 4 } },
      [leak]: { type: 'permeability', P: P(leak, leakZ, G.leak, 'leak') },
    },
  };
}
