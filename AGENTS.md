You are an experienced, pragmatic software engineering AI agent. Do not over-engineer a solution when a simple one is possible. Keep edits minimal. If you want an exception to ANY rule, you MUST stop and get permission first.

# Project Overview

**autoeq-c** is a fast AutoEQ (automatic equalization) engine by PEQdB Inc. It fits parametric EQ filters (Peaking, Low Shelf, High Shelf) to match a target frequency response given a source measurement. The core algorithm uses gradient-based optimization (AdaBelief) with analytical derivatives to minimize MSE between the residual and zero across a fixed 384-point log-spaced frequency grid (20 Hz–20 kHz).

The project produces two build targets via a single `build.zig`:
1. **WebAssembly module** (native Zig cross-compilation, no Emscripten) — consumed from JS/TS
2. **Native CLI binary** — consumed via stdin/stdout (float32 binary input, text output)

Licensed under LGPL-3.0-or-later.

## Technology Stack

| Layer | Technology |
|-------|-----------|
| Core language | Zig 0.15.2 |
| Build system | `build.zig` (Zig's built-in build system) |
| WASM build | `zig build` cross-compilation to `wasm32-freestanding` |
| JS wrapper | Vanilla JavaScript (ES modules, JSDoc types, raw WebAssembly API) |
| Tests | `zig build test` (Zig's built-in test runner) |
| No package manager, no npm/pip, no CI |

# Reference

## Directory Structure

```
src-zig/
  types.zig    — Types (Biquad, Filter, Lim, Type enum), constants (MAX_N=32, K=384)
  biquad.zig   — Biquad coefficient computation with analytical derivatives, spectrum()
  init.zig     — Filter initialization via peak-finding (adapted from scipy/jaakkopasanen/AutoEq)
  smooth.zig   — Smoothing presets (in-ear, over-ear)
  optimize.zig — AdaBelief optimizer, gradient computation, fitting
  root.zig     — Preprocessing, top-level autoeq orchestration
  cli.zig      — Native CLI main (stdin/stdout)
  wasm.zig     — WebAssembly exports
example/
  autoeq.js       — JS/TS wrapper for WASM module (raw WebAssembly API, memory management, configs, Smooth enum)
  autoeq-wasm.d.ts — Minimal TS declaration for WASM module
  test.ts         — TypeScript usage example
  test.py         — Python wrapper for native binary (subprocess + float32 stdin)
build.zig         — Zig build file (native + WASM targets, test step)
```

## Architecture

```
dst (target) + src (measurement)
  → preprocess(): smooth src, compute residual = dst - smoothed_src, demean, treble rolloff
  → autoeq(): for each filter, init via peak-finding, subtract spectrum, then fit all jointly
  → fit(): AdaBelief optimizer minimizing MSE over 3N+1 parameters
  → Output: {f0, gain, Q} per filter + overall amp
```

Key design decisions:
- **Reparameterization**: Optimizes in log-frequency and bandwidth space for better conditioning
- **Analytical gradients**: Hand-derived partial derivatives for all biquad coefficients — core performance trick
- **Box constraints via projection**: Parameters clipped after each step; optimizer running mean zeroed for constrained params
- **Fixed frequency grid**: 384 log-spaced points, all computation on this grid
- **Greedy sequential initialization**: Find largest peak → fit filter → subtract → repeat

## Zig Conventions

- `@setFloatMode(.optimized)` on hot paths for fast-math semantics (equivalent to C's `-ffast-math`)
- `@import` for module dependencies; `pub` for exported symbols
- Slices and arrays indexed with `[n]` syntax; `for` loops over slices
- `comptime` where possible for compile-time evaluation
- `std.mem.copy` / `std.mem.set` instead of `memcpy`/`memset`
- Exported WASM functions use `export fn` with ABI-compatible types
- Native CLI uses Zig's `std.process` for argument parsing and I/O
- Tabs for indentation (per `.editorconfig`)

## Lessons Learned

- **WASM build**: Zig cross-compiles to `wasm32-freestanding` natively — no Emscripten needed. The WASM module exports flat C-ABI functions.
- **Float mode**: `@setFloatMode(.optimized)` must be set on hot-path functions (biquad, optimizer) to get performance comparable to C's `-ffast-math`.
- **Performance (overview)**: Zig is ~14% slower than C (`-O3 -ffast-math`) on this workload. The gap is **entirely due to glibc's libmvec** (vectorized math). C auto-vectorizes transcendental calls (sin, cos, log, exp, sinh, asinh) into `_ZGVbN4v_*` SIMD implementations; Zig uses scalar `compiler-rt` math. `@Vector(4, f32)` helps with arithmetic (dot products, polynomial eval) but ~90% of runtime is in transcendental functions. `exe.linkSystemLibrary("m")` does NOT help — Zig's built-in math takes precedence. Possible fixes: call libmvec directly via `extern fn` + inline asm (Linux/glibc only), wait for Zig veclib support, override compiler-rt per-function, or use polynomial approximations. Native binary ~1.1MB vs C's ~39KB.
- **JS wrapper**: Updated from Emscripten module to raw WebAssembly API (`WebAssembly.instantiateStreaming`). Same public interface preserved.
- **Binary sizes**: Native ~1.1MB, WASM ~641KB.

# Essential Commands

```sh
# Build everything (native + WASM)
zig build -Doptimize=ReleaseFast

# Run tests
zig build test

# Native CLI usage (binary is at zig-out/bin/autoeq)
# i=in-ear smooth, o=over-ear smooth, x=no smooth
# N=number of filters (1-32), steps=optimizer iterations (default 3000)
./zig-out/bin/autoeq [iox] <N> <?steps>
# Input: 768 float32 values via stdin (384 dst + 384 src)
```

No linter or formatter is configured beyond `.editorconfig`.

# Commit and Pull Request Guidelines

- Commit messages: lowercase, brief imperative (e.g., `add python example`, `readme changes`, `ts -> js`).
- Validate with `zig build` and `zig build test` before committing.
- No PR template or issue templates exist.
# Notes

You have MemPalace agents. Run mempalace_list_agents to see them.

