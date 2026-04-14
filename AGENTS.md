You are an experienced, pragmatic software engineering AI agent. Do not over-engineer a solution when a simple one is possible. Keep edits minimal. If you want an exception to ANY rule, you MUST stop and get permission first.

# Project Overview

**autoeq-c** is a fast AutoEQ (automatic equalization) engine written in C by PEQdB Inc. It fits parametric EQ filters (Peaking, Low Shelf, High Shelf) to match a target frequency response given a source measurement. The core algorithm uses gradient-based optimization (AdaBelief) with analytical derivatives to minimize MSE between the residual and zero across a fixed 384-point log-spaced frequency grid (20 Hz–20 kHz).

The project produces two build targets:
1. **WebAssembly module** (via Emscripten) — consumed from JS/TS
2. **Native CLI binary** — consumed via stdin/stdout (float32 binary input, text output)

Licensed under LGPL-3.0-or-later.

## Technology Stack

| Layer | Technology |
|-------|-----------|
| Core language | C (C11: `_Static_assert`, `bool`, compound literals, `restrict`) |
| WASM build | Emscripten (`emcc`) with `-O3 -ffast-math -msimd128` |
| JS wrapper | Vanilla JavaScript (ES modules, JSDoc types) |
| Build tooling | Single `./build` shell script; `cc` for native |
| No package manager, no npm/pip, no Makefile, no CI |

# Reference

## Directory Structure

```
src/
  lib.h        — Types (Biquad, Filter, Lim, Type enum), constants (MAX_N=32, K=384), inline helpers
  biquad.c     — Biquad coefficient computation with analytical derivatives, spectrum()
  init.c       — Filter initialization via peak-finding (adapted from scipy/jaakkopasanen/AutoEq)
  autoeq.c     — Core: gradient, AdaBelief optimizer, preprocessing, CLI main()
example/
  autoeq.js    — JS/TS wrapper for WASM module (memory management, configs, Smooth enum)
  autoeq-wasm.d.ts — Minimal TS declaration for WASM module
  test.ts      — TypeScript usage example
  test.py      — Python wrapper for native binary (subprocess + float32 stdin)
web/
  pre.js       — Emscripten --pre-js: routes stderr to console.debug
build          — Shell script: emcc build for WASM
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

## C Conventions

- Custom typedefs: `i32`, `u32`, `f32`, `f64`
- Macros: `FOR_K()`, `FOR_N()`, `FOR_W()` for grid/filter iteration
- `restrict` qualifiers throughout for aliasing optimization
- `__builtin_expect` for branch hints
- Designated initializers for arrays/structs
- `#ifdef WASM` for conditional compilation (WASM vs native)
- `INFO()` macro → `fprintf(stderr, ...)` for logging
- Tabs for indentation (per `.editorconfig`)

# Essential Commands

```sh
# WASM build (requires emcc in PATH)
./build

# Native build
cc src/*.c -lm -o autoeq

# Native CLI usage
# i=in-ear smooth, o=over-ear smooth, x=no smooth
# N=number of filters (1-32), steps=optimizer iterations (default 3000)
./autoeq [iox] <N> <?steps>
# Input: 768 float32 values via stdin (384 dst + 384 src)
```

No test framework, linter, or formatter is configured.

# Commit and Pull Request Guidelines

- Commit messages: lowercase, brief imperative (e.g., `add python example`, `readme changes`, `ts -> js`).
- Validate with `./build` (WASM) and/or `cc src/*.c -lm -o autoeq` (native) before committing.
- No PR template or issue templates exist.
# Notes

You have MemPalace agents. Run mempalace_list_agents to see them.

