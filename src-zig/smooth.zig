// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later

const std = @import("std");
const types = @import("types.zig");

const K = types.K;
const Lim = types.Lim;
const sq = types.sq;

pub const Smooth = extern struct {
    smooth_f0: f32, smooth_f1: f32,

    smooth_lo: f32,
    smooth_hi: f32,

    bias_f0: f32, bias_f1: f32,
    bias_f2: f32, bias_f3: f32,

    bias_lo: f32,
    bias_md: f32,
    bias_hi: f32,

    clip_f: f32,
};

pub const IE_SMOOTH: Smooth = .{
    .smooth_lo = 0.3,
    .smooth_hi = 0.03,

    .smooth_f0 = 3000.0, .smooth_f1 = 12000.0,

    .bias_lo = 0.0,
    .bias_md = 0.15,
    .bias_hi = 0.03,

    .bias_f0 = 10000.0, .bias_f1 = 13000.0,
    .bias_f2 = 14000.0, .bias_f3 = 20000.0,

    .clip_f = 18500.0,
};

pub const OE_SMOOTH: Smooth = .{
    .smooth_lo = 0.3,
    .smooth_hi = 0.03,

    .smooth_f0 = 5000.0, .smooth_f1 = 15000.0,

    .bias_lo = 0.0,
    .bias_md = 0.3,
    .bias_hi = 0.2,

    .bias_f0 = 6000.0, .bias_f1 = 9000.0,
    .bias_f2 = 9000.0, .bias_f3 = 20000.0,

    .clip_f = 17000.0,
};

/// Binary search — x[] is monotonic
pub fn search(x: [*]const f32, n: usize, v: f32) i32 {
    var lo: i32 = 0;
    var hi: i32 = @intCast(n - 1);
    while (lo < hi) {
        const mid = lo + @divTrunc(hi - lo, 2);
        if (x[@intCast(mid)] < v)
            lo = mid + 1
        else
            hi = mid;
    }
    if (lo > 0 and @abs(x[@intCast(lo - 1)] - v) < @abs(x[@intCast(lo)] - v))
        lo -= 1;
    return lo;
}

/// Sigmoid-based step function
fn sgm(x: f32, x0: f32, x1: f32) f32 {
    const SMOOTH: f32 = 4.0;
    const k = SMOOTH / (x1 - x0);
    const m = 0.5 * (x0 + x1);
    const y = k * (x - m);
    return 0.5 * std.math.tanh(0.5 * y) + 0.5;
}

/// Progressively smooth data from 1kHz onwards using bilateral-ish filtering.
pub fn adaptive_smooth(s: *const Smooth, f: [*]const f32, r: [*]f32) void {
    comptime {
        std.debug.assert(K == 384);
    }
    const H: i32 = 48;

    const smooth_l0 = @log(s.smooth_f0);
    const smooth_l1 = @log(s.smooth_f1);

    const bias_l0 = @log(s.bias_f0);
    const bias_l1 = @log(s.bias_f1);
    const bias_l2 = @log(s.bias_f2);
    const bias_l3 = @log(s.bias_f3);

    var x: [K]f32 = undefined;
    @memcpy(&x, r[0..K]);

    const clip_idx = search(f, K, s.clip_f);

    for (0..K) |k_| {
        const f_k = f[k_];
        const l = @log(f_k);
        const x_k = x[k_];

        const sigma = s.smooth_lo + (s.smooth_hi - s.smooth_lo) * sgm(l, smooth_l0, smooth_l1);
        const bias = s.bias_lo + (s.bias_md - s.bias_lo) * sgm(l, bias_l0, bias_l1)
            + (s.bias_hi - s.bias_md) * sgm(l, bias_l2, bias_l3);

        var a: f32 = 0.0;
        var c: f32 = 0.0;

        var j: i32 = -H;
        while (j <= H) : (j += 1) {
            var jj: i32 = @as(i32, @intCast(k_)) + j;
            if (jj < 0) jj = 0;
            if (jj > clip_idx) jj = clip_idx;

            const x_jj = x[@intCast(jj)];
            const d_spatial = sq(@as(f32, @floatFromInt(j)) * sigma);
            const d_range = bias * (x_jj - x_k);

            const w = @exp(-0.5 * d_spatial + d_range);

            a += w * x[@intCast(jj)];
            c += w;
        }

        r[k_] = a / c;
    }
}

pub fn center_mean(x: [*]f32, n: usize) f32 {
    var sum: f32 = 0.0;
    for (0..n) |i| sum += x[i];
    const mean = sum / @as(f32, @floatFromInt(n));
    for (0..n) |i| x[i] -= mean;
    return mean;
}

pub fn treble_rolloff(f: [*]const f32, r: [*]f32, f_treble: f32) void {
    const treble_idx: usize = @intCast(search(f, K, f_treble));
    const n_treble: usize = K - treble_idx;
    const inv = 1.0 / @as(f32, @floatFromInt(n_treble - 1));

    for (0..n_treble) |i| {
        const t: f32 = @floatFromInt(i);
        const w = @cos(0.5 * std.math.pi * t * inv);
        r[treble_idx + i] *= w;
    }
}

pub fn preprocess(
    f: [*]const f32,
    dst: [*]const f32,
    src: [*]const f32,
    r: [*]f32,
    smooth: ?*const Smooth,
    demean: bool,
) f32 {
    const F_TREBLE_SMOOTH: f32 = 16000.0;
    const F_TREBLE_UNSMOOTH: f32 = 18500.0;

    var b: [K]f32 = undefined;
    @memcpy(&b, src[0..K]);

    if (smooth) |s| adaptive_smooth(s, f, &b);

    for (0..K) |k_| r[k_] = dst[k_] - b[k_];

    var mean: f32 = 0.0;
    if (demean)
        mean = center_mean(r, K);

    treble_rolloff(f, r, if (smooth != null) F_TREBLE_SMOOTH else F_TREBLE_UNSMOOTH);

    return mean;
}
