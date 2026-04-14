# autoeq-c

Blazing fast and accurate AutoEQ engine that produces high-quality, rig-aware results.

Originally written in C, now implemented in Zig (0.15.2). The original C source in `src/` is kept for reference.

This repo produces two build targets via a single `build.zig`:

- **WebAssembly module** — consumed from JS/TS (raw WebAssembly API, no Emscripten)
- **Native CLI binary** — consumed via stdin/stdout (float32 binary input, text output)

## Build

Requirements:

- [Zig](https://ziglang.org/) 0.15.2

```sh
zig build -Doptimize=ReleaseFast
```

This produces:
- `zig-out/bin/autoeq` — native CLI binary
- `zig-out/bin/autoeq.wasm` — WebAssembly module

Run tests:

```sh
zig build test
```

The WASM module is consumed by the JavaScript wrapper in [example/autoeq.js](example/autoeq.js).

Three functions are provided:

- `make()`
- `run(inst, dst, src, config, smooth, steps, fs)`
- `interp(x, y)`

Notes:

- The frequency axis is fixed to 384 log-spaced points `X` from 20 Hz → 20 kHz. If your data isn’t already on the internal `X` grid, use the provided linear interpolation function `interp(x, y)` from `example/autoeq.js`.
- `interp()` is linear interpolation. If you resample very dense data (e.g. raw FFT / measurement output) directly down to 384 points you will introduce aliasing artifacts. Prefer to bin/average onto a log-spaced grid and/or apply some smoothing before interpolating.
- The number of filters must be <= 32. Default configs are provided as `autoeq.CONFIGS.STANDARD(N)` and `autoeq.CONFIGS.PRECISE(N)`.

### TypeScript example ([example/test.ts](example/test.ts))

```ts
import * as autoeq from './autoeq';

const inst = await autoeq.make();

// suppose you have arbitrary measurement points
// resample them onto the internal 384-point grid `autoeq.X` before running

const x = [20, 50, 100, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000];

// these are just example curves
const dstRaw = x.map(x => 6*Math.sin(x)); // target curve in dB at x
const srcRaw = x.map(x => 6*Math.cos(x)); // source curve in dB at x

const dst = autoeq.interp(x, dstRaw);
const src = autoeq.interp(x, srcRaw);

// autoeq src to dst with 8 filters and IE smoothing

const config = autoeq.CONFIGS.STANDARD(8);
const res = autoeq.run(inst, dst, src, config, autoeq.Smooth.IE);

if (res) {
    console.log(res.filters);
    console.log({ amp: res.amp, loss: res.loss, time: res.time });
}
```

Example output:

```
[
  {
    type: "LSC",
    f0: 23.711910247802734,
    gain: 7.024210453033447,
    q: 0.565746009349823,
  },
  /* ... */
  {
    type: "PK",
    f0: 9249.556640625,
    gain: 3.049382209777832,
    q: 1.3290164470672607,
  }
]
{
  amp: 1.4850640296936035,
  loss: 0.027285713702440262,
  time: 30.329625,
}
```

## Native CLI usage

The native binary is built by the same `zig build -Doptimize=ReleaseFast` command above (outputs to `zig-out/bin/autoeq`).

It is a low-level tool meant to be called from a script or part of a larger program.

- `./autoeq [iox] <N> <?steps>`

Where:

- `iox`
    - `i` = In-ear tuned smoothing preset
    - `o` = Over-ear smoothing preset
    - `x` = No smoothing
- `N` = Number of filters to fit (<= 32)
- `steps` (optional) = Number of optimizer steps

Input (stdin):

- 384 `float32` values for `dst` followed by 384 `float32` values for `src`

Output (stdout):

- First line: `<amp> <max>`
  - `amp` is the fitted overall gain offset (includes mean removal)
  - `max` is the max dB of the fitted filter stack response
- Then `N` lines, one per filter:
  - `<type> <f0> <gain> <Q>`

There is a minimal Python example in [example/test.py](example/test.py) with `NumPy` and `subprocess` that demonstrates how to call the native binary.

## License

LGPL-3.0-or-later. See [LICENSE](LICENSE).
