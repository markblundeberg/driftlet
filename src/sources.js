// driftlet/kit: sources as plain definition data. Waveforms for a terminal's V or I, an injector
// port, and band-to-band recombination from a lifetime. Each returns an ordinary spec value.

import { DeviceError } from './errors.js';

const fail = (message) => {
  throw new DeviceError(message);
};
const need = (ok, message) => ok || fail(message);
const num = (v) => typeof v === 'number' && Number.isFinite(v);

/**
 * A rectangular pulse, as a waveform for a terminal's `V` or `I`: `base` until `start`, then
 * `base + amplitude` for `width`, then `base` again. With a `rise` > 0 its edges are ramps of
 * that length, placed so the area stays `amplitude · width`; without one they're jumps, which
 * the time stepping handles exactly (so a pulse with no rise must start after t = 0, the time a
 * device starts at and solve() reads: it's off then).
 * @param {{ start?: number, width: number, amplitude: number, rise?: number, base?: number }} opts s, s, V or A/m², s, V or A/m²
 * @returns {{ t: number[], values: number[] }}
 */
export function pulse({ start = 0, width, amplitude, rise = 0, base = 0 } = {}) {
  need(num(width) && width > 0, `pulse: width must be > 0 (s), got ${width}`);
  need(num(amplitude), `pulse: amplitude must be a finite number (V or A/m²), got ${amplitude}`);
  need(num(start) && start >= 0, `pulse: start must be ≥ 0 (s), got ${start}`);
  need(num(rise) && rise >= 0 && rise <= width, `pulse: rise must be between 0 and the width (s), got ${rise}`);
  need(num(base), `pulse: base must be a finite number, got ${base}`);
  need(rise > 0 || start > 0, 'pulse: with no rise, start after t = 0, so the device starts (and solve() reads it) before the pulse');
  const top = base + amplitude;
  return { t: [start, start + rise, start + width, start + width + rise], values: [base, top, top, base] };
}

/**
 * A square wave, repeating: `high` for `duty · period`, then `low`, from t = 0.
 * @param {{ high: number, low: number, period: number, duty?: number }} opts
 * @returns {{ t: number[], values: number[], repeat: true }}
 */
export function square({ high, low, period, duty = 0.5 } = {}) {
  need(num(high) && num(low), `square: high and low must be finite numbers, got ${high} and ${low}`);
  need(num(period) && period > 0, `square: period must be > 0 (s), got ${period}`);
  need(num(duty) && duty > 0 && duty < 1, `square: duty must be between 0 and 1, got ${duty}`);
  const on = duty * period;
  return { t: [0, on, on, period], values: [high, high, low, low], repeat: true };
}

/**
 * A triangle wave, repeating: from `from` to `to` and back in each `period`, from t = 0 (a
 * cyclic voltammogram's sweep, from its starting potential).
 * @param {{ from: number, to: number, period: number }} opts
 * @returns {{ t: number[], values: number[], repeat: true }}
 */
export function triangle({ from, to, period } = {}) {
  need(num(from) && num(to), `triangle: from and to must be finite numbers, got ${from} and ${to}`);
  need(num(period) && period > 0, `triangle: period must be > 0 (s), got ${period}`);
  return { t: [0, period / 2, period], values: [from, to, from], repeat: true };
}

/**
 * A ramp from `from` to `to` over `duration`, starting at `start`, constant on either side.
 * @param {{ from: number, to: number, duration: number, start?: number }} opts
 * @returns {{ t: number[], values: number[] }}
 */
export function ramp({ from, to, duration, start = 0 } = {}) {
  need(num(from) && num(to), `ramp: from and to must be finite numbers, got ${from} and ${to}`);
  need(num(duration) && duration > 0, `ramp: duration must be > 0 (s), got ${duration}`);
  need(num(start) && start >= 0, `ramp: start must be ≥ 0 (s), got ${start}`);
  return { t: [start, start + duration], values: [from, to] };
}

/**
 * A port that injects a species into a window of a region at a driven current `I` (A/m², a
 * number or a waveform), spread evenly over the window. It's an ordinary port: driven by `I`,
 * its terminal the species, linked by a conductance small enough that at the peak current the
 * port sits `headroom` volts above the window's level of that species, so every node of the
 * window takes the same share (like a current source). At I = 0 it still links the window,
 * weakly.
 * @param {{ name?: string, region: string | number, from: number, to: number, species: string,
 *   I: number | { t: number[], values: number[], repeat?: boolean }, headroom?: number }} opts
 * @returns {object} a port, for the definition's `ports`
 */
export function injector({ name = 'injector', region, from, to, species, I, headroom = 2.5 } = {}) {
  need(region !== undefined, 'injector: give the region (its name or index)');
  need(num(from) && num(to) && to > from, `injector: from and to must bound the window, in m from the region's left end (got ${from}, ${to})`);
  need(typeof species === 'string', 'injector: give the species it injects, e.g. species: \'h+\'');
  need(num(headroom) && headroom > 0, `injector: headroom must be > 0 (V), got ${headroom}`);
  const peak = num(I) ? Math.abs(I) : Array.isArray(I?.values) && I.values.every(num) ? Math.max(...I.values.map(Math.abs)) : NaN;
  need(Number.isFinite(peak), 'injector: I must be a current (A/m²) or a waveform { t, values }');
  const G = peak > 0 ? peak / ((to - from) * headroom) : 1; // S/m³
  return { name, region, from, to, terminal: species, species: { [species]: { type: 'conductance', G } }, I };
}

/**
 * Band-to-band recombination e⁻ + h⁺ ⇌ 0 in a material, from the minority-carrier lifetime it
 * should give at a majority density: kf = 1/(τ · majority). (The rate is kf (np − n_i²), so τ is
 * the low-injection lifetime of whichever carrier is the minority.)
 * @param {{ material: string, tau: number, majority: number }} opts material name, s, mol/m³
 * @returns {{ equation: string, kf: Record<string, number> }} a bulk reaction
 */
export function recombination({ material, tau, majority } = {}) {
  need(typeof material === 'string', 'recombination: give the material, by name');
  need(num(tau) && tau > 0, `recombination: tau must be > 0 (s), got ${tau}`);
  need(num(majority) && majority > 0, `recombination: majority must be the majority carrier density > 0 (mol/m³; units.perCm3 converts), got ${majority}`);
  return { equation: 'e- + h+ = 0', kf: { [material]: 1 / (tau * majority) } };
}

/**
 * Photogeneration by light absorbed in a material (Beer–Lambert): a photon flux entering at x =
 * `from` and travelling toward +x, generating pairs at G(x) = flux · α · e^{−α(x − from)} up to
 * x = `to`. Written as photon ⇌ e⁻ + h⁺ (or whatever `makes` says, such as an exciton 'X') from a
 * photon reservoir at μ = `mu` (J/mol; well above what it makes, so the rate is G itself), with
 * kf as a profile against the device's x. The table is
 * exact in what it generates in all (each segment's mean is the exponential's), and the solver
 * averages it over each node's box, so the total is right on any grid.
 * @param {{ material: string, flux: number, alpha: number, mu: number, from?: number, to: number, makes?: string }} opts
 *   material name; photons absorbed, mol/(m²·s); absorption coefficient, 1/m; the photons' μ,
 *   J/mol; where the light enters and where the material ends (m, the device's x); the
 *   reaction's right-hand side (default 'e- + h+')
 * @returns {{ equation: string, fixed: Record<string, number>, kf: Record<string, { x: number[], values: number[] }> }} a bulk reaction
 */
export function photogeneration({ material, flux, alpha, mu, from = 0, to, makes = 'e- + h+' } = {}) {
  need(typeof material === 'string', 'photogeneration: give the material, by name');
  need(num(flux) && flux >= 0, `photogeneration: flux must be the photon flux ≥ 0 entering the material (mol/(m²·s)), got ${flux}`);
  need(num(alpha) && alpha > 0, `photogeneration: alpha must be the absorption coefficient > 0 (1/m; 1/cm × 100), got ${alpha}`);
  need(num(mu), `photogeneration: mu must be the photons' μ (J/mol), well above the gap: units.eV(3), say`);
  need(typeof makes === 'string' && makes.trim() !== '', `photogeneration: makes must be what the light makes, as an equation's side ('e- + h+', 'X'), got ${makes}`);
  need(num(from) && num(to) && to > from, `photogeneration: from and to must be positions with to > from (m, the device's x), got ${from} and ${to}`);
  // Evenly spaced points, at most 0.05 absorption lengths apart, out to 50 of them (beyond,
  // nothing is left), then `to`.
  const depth = alpha * (to - from), reach = Math.min(depth, 50), n = Math.ceil(reach / 0.05), step = reach / n;
  const taus = Array.from({ length: n + 1 }, (_, k) => k * step);
  if (depth > 50) taus.push(depth);
  // A straight line between two points of e^{−τ} encloses (Δτ/2)·coth(Δτ/2) times too much: each
  // value is scaled back by that, so the total is the exponential's own (to ~Δτ⁴).
  const fix = (step / 2) / Math.tanh(step / 2);
  return {
    equation: `photon = ${makes}`,
    fixed: { photon: mu },
    kf: { [material]: { x: taus.map((t) => from + t / alpha), values: taus.map((t) => (flux * alpha * Math.exp(-t)) / fix) } },
  };
}
