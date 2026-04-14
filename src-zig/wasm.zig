// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later

const std = @import("std");
const root = @import("autoeq");

const K = root.K;
const Type = root.Type;
const Lim = root.Lim;
const Smooth = root.Smooth;
const IE_SMOOTH = root.IE_SMOOTH;
const OE_SMOOTH = root.OE_SMOOTH;
const preprocess = root.preprocess;
const autoeq = root.autoeq;

// Simple bump allocator for WASM linear memory
var heap_offset: usize = 0x10000;

export fn malloc(size: usize) [*]u8 {
    const aligned = (size + 15) & ~@as(usize, 15);
    const ptr = heap_offset;
    heap_offset += aligned;
    const needed = heap_offset;
    const current = @wasmMemorySize(0) * 65536;
    if (needed > current) {
        const pages = (needed - current + 65535) / 65536;
        _ = @wasmMemoryGrow(0, pages);
    }
    return @ptrFromInt(ptr);
}

export fn free(ptr: [*]u8) void {
    _ = ptr;
}

export fn preprocess_export(
    f: [*]const f32,
    dst: [*]const f32,
    src: [*]const f32,
    r: [*]f32,
    smooth: ?*const Smooth,
    demean: bool,
) f32 {
    return preprocess(f, dst, src, r, smooth, demean);
}

export fn autoeq_export(
    steps: i32,
    types: [*]const Type,
    f0: [*]f32,
    gain: [*]f32,
    Q: [*]f32,
    amp: ?*f32,
    f0_lim: [*]const Lim,
    gain_lim: [*]const Lim,
    Q_lim: [*]const Lim,
    N: usize,
    f: [*]const f32,
    r: [*]const f32,
    fs: f32,
) f32 {
    return autoeq(steps, types, f0, gain, Q, amp, f0_lim, gain_lim, Q_lim, N, f, r, fs);
}

export fn get_ie_smooth() *const Smooth {
    return &IE_SMOOTH;
}

export fn get_oe_smooth() *const Smooth {
    return &OE_SMOOTH;
}

pub fn panic(msg: []const u8, _: ?*std.builtin.StackTrace, _: ?usize) noreturn {
    _ = msg;
    while (true) {}
}
