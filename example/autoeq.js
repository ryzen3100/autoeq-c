/*
 * Copyright (C) 2026 PEQdB Inc.
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/* Build with: zig build (produces zig-out/bin/autoeq-wasm.wasm) */

const N = 384,
      X0 = 20,
      X1 = 20_000,
      LX0 = Math.log(X0),
      LX1 = Math.log(X1),
      LXR = LX1 - LX0;

export const LX = Array.from({ length: N }, (_, i) => LXR/(N - 1)*i + LX0);
export const X  = LX.map(Math.exp);

/**
 * @typedef {'PK'|'LSC'|'HSC'} FilterType
 */

/**
 * @typedef {{
 *   type: FilterType,
 *   f0: number,
 *   gain: number,
 *   q: number,
 * }} Filter
 */

/**
 * @readonly
 * @enum {number}
 */
const Type = Object.freeze({
  PK: 0,
  LSC: 1,
  HSC: 2,
});

/** @type {ReadonlyArray<FilterType>} */
const TYPE_NAMES = Object.freeze(['PK', 'LSC', 'HSC']);

/**
 * @readonly
 * @enum {number}
 */
export const Smooth = Object.freeze({
  NONE: 0,
  IE: 1,
  OE: 2,
});

/**
 * @typedef {[number, number]} Lim
 */

/**
 * @typedef {{
 *   type: number,
 *   f0: Lim,
 *   gain: Lim,
 *   q: Lim,
 * }} Spec
 */

/** @typedef {{ m: WebAssembly.Instance }} Inst */

/** @param {string|URL|ArrayBuffer} [wasmSource] @returns {Promise<Inst>} */
export async function make(wasmSource) {
  let bytes;
  if (wasmSource instanceof ArrayBuffer) {
    bytes = wasmSource;
  } else if (typeof wasmSource === 'string' || wasmSource instanceof URL) {
    const url = typeof wasmSource === 'string' ? new URL(wasmSource, import.meta.url) : wasmSource;
    if (url.protocol === 'file:') {
      const fs = await import('node:fs/promises');
      bytes = await fs.readFile(url);
    } else {
      const resp = await fetch(url);
      bytes = await resp.arrayBuffer();
    }
  } else {
    const url = new URL('./autoeq-wasm.wasm', import.meta.url);
    if (url.protocol === 'file:') {
      const fs = await import('node:fs/promises');
      bytes = await fs.readFile(url);
    } else {
      const resp = await fetch(url);
      bytes = await resp.arrayBuffer();
    }
  }
  const { instance } = await WebAssembly.instantiate(bytes);
  return { m: instance };
}

/**
 * @typedef {{
 *   specs: Spec[],
 *   smooth: boolean,
 *   demean: boolean,
 * }} Config
 */

/**
 * @type {{
 *   STANDARD: (n?: number) => Config,
 *   PRECISE: (n?: number) => Config,
 * }}
 */
export const CONFIGS = {
  STANDARD: (n = 10) => ({
    specs: [
      { type: Type.LSC, f0: [20, 16_000], gain: [-16, 16], q: [.4, 3.] },
      { type: Type.HSC, f0: [20, 16_000], gain: [-16, 16], q: [.4, 3.] },
      ...Array(n - 2).fill(
        { type: Type.PK , f0: [20, 16_000], gain: [-16, 16], q: [.4, 4.] })
    ],
    smooth: true,
    demean: true,
  }),
  PRECISE: (n = 10) => ({
    specs: [
      { type: Type.LSC, f0: [20, 16_000], gain: [-16, 16], q: [.4, 3.] },
      { type: Type.HSC, f0: [20, 16_000], gain: [-16, 16], q: [.4, 3.] },
      ...Array(n - 2).fill(
        { type: Type.PK , f0: [20, 16_000], gain: [-16, 16], q: [.4, 4.] })
    ],
    smooth: false,
    demean: true,
  }),
};

/**
 * @param {Inst} s
 * @param {ArrayLike<number>} dst
 * @param {ArrayLike<number>} src
 * @param {Config} config
 * @param {number} [smooth]
 * @param {number} [steps]
 * @param {number} [fs]
 * @returns {{ filters: Filter[], time: number, loss: number, amp: number } | null}
 */
export function run(s, dst, src, config, smooth = Smooth.NONE, steps = 3000, fs = 48_000) {
  if (dst.length !== N || src.length !== N) {
    console.error(`dst and src must have length ${N}`);
    return null;
  }

  const types = Int32Array.from(config.specs.map(s => s.type));

  const exp = s.m.exports;
  const mem = exp.memory.buffer;

  let allocs = [];

  function alloc(sz) {
    const p = exp.malloc(sz);
    allocs.push(p);
    return p;
  }

  const k = N,
        n = types.length;

  const pF   = alloc(k * 4),
        pDst = alloc(k * 4),
        pSrc = alloc(k * 4),
        pR   = alloc(k * 4);

  new Float32Array(mem).set(X, pF >> 2);
  new Float32Array(mem).set(dst, pDst >> 2);
  new Float32Array(mem).set(src, pSrc >> 2);

  const pTypes = alloc(n * 4),
        pF0    = alloc(n * 4),
        pGain  = alloc(n * 4),
        pQ     = alloc(n * 4),
        pAmp   = alloc(4);

  new Int32Array(mem).set(types, pTypes >> 2);

  const pF0Lim   = alloc(n * 8),
        pGainLim = alloc(n * 8),
        pQLim    = alloc(n * 8);

  const heap = new Float32Array(mem);
  config.specs.forEach((spec, i) => {
    let p = (pF0Lim >> 2) + 2*i;
    heap[p + 0] = spec.f0[0];
    heap[p + 1] = spec.f0[1];

    p = (pGainLim >> 2) + 2*i;
    heap[p + 0] = spec.gain[0];
    heap[p + 1] = spec.gain[1];

    p = (pQLim >> 2) + 2*i;
    heap[p + 0] = spec.q[0];
    heap[p + 1] = spec.q[1];
  });

  let pSmooth = !config.smooth ? 0
    : smooth == Smooth.IE ? exp.get_ie_smooth()
    : smooth == Smooth.OE ? exp.get_oe_smooth()
    : 0;

  const t0 = performance.now();

  const mean = exp.preprocess_export(pF, pDst, pSrc, pR, pSmooth, config.demean);
  const loss = exp.autoeq_export(
    steps, pTypes, pF0, pGain, pQ, config.demean ? pAmp : 0,
    pF0Lim, pGainLim, pQLim,
    n, pF, pR, fs);

  const t1 = performance.now();

  const f0   = new Float32Array(mem, pF0, n),
        gain = new Float32Array(mem, pGain, n),
        q    = new Float32Array(mem, pQ, n);

  const amp = config.demean ? new Float32Array(mem)[pAmp >> 2] : 0;

  /** @type {Filter[]} */
  const filters = [];

  for (let i = 0; i < n; i++) {
    filters.push({
      type: TYPE_NAMES[types[i]] ?? /** @type {any} */ (types[i]),
      f0: f0[i],
      gain: gain[i],
      q: q[i],
    });
  }

  allocs.forEach(p => exp.free(p));

  return {
    filters,
    time: t1 - t0,
    loss,
    amp: mean + amp,
  };
}

/**
 * @param {number[]} x
 * @param {number[]} y
 * @returns {number[]}
 */
export function interp(x, y) {
  if (x.length !== y.length) {
    console.error('x and y must have the same length');
    return [];
  }

  if (x.length < 2) {
    console.error('x and y must have length >= 2');
    return [];
  }

  const n = x.length;

  const lx = x.map(Math.log);
  const out = Array(X.length);

  let i = 0;
  for (let j = 0; j < X.length; ++j) {
    const t = Math.log(X[j]);

    if (t <= lx[0]) {
      out[j] = y[0];
      continue;
    }

    if (t >= lx[n - 1]) {
      out[j] = y[n - 1];
      continue;
    }

    while (i + 1 < n - 1 && lx[i + 1] < t)
      ++i;

    const x0 = lx[i],
          x1 = lx[i + 1];

    const den = x1 - x0;
    const u = den === 0 ? 0 : (t - x0) / den;

    out[j] = y[i] + u * (y[i + 1] - y[i]);
  }

  return out;
}
