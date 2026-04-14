// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later

const std = @import("std");
const types = @import("types.zig");
const biquad = @import("biquad.zig");

const K = types.K;
const Filter = types.Filter;
const Lim = types.Lim;
const Type = types.Type;
const limit = types.limit;
const spectrum = biquad.spectrum;

const Peak = struct {
    width: f32,
    height: f32,
    idx: i32,
};

/// Function adapted from scipy
fn largest_peak(x: [*]const f32, f: [*]const f32, lim: Lim) Peak {
    const H: usize = K / 2;

    var peaks: [H]i32 = undefined;
    var n: usize = 0;
    const i_max: i32 = @as(i32, @intCast(K)) - 1;

    var i: i32 = 1;
    while (i < i_max) : (i += 1) {
        if (f[@intCast(i)] < lim.lo or f[@intCast(i)] > lim.hi)
            continue;

        if (x[@intCast(i - 1)] >= x[@intCast(i)])
            continue;

        var i_ahead: i32 = i + 1;
        while (i_ahead < i_max and x[@intCast(i_ahead)] == x[@intCast(i)])
            i_ahead += 1;

        if (x[@intCast(i_ahead)] < x[@intCast(i)]) {
            const left_edge = i;
            const right_edge = i_ahead - 1;
            peaks[n] = @divTrunc(left_edge + right_edge, 2);
            n += 1;
            i = i_ahead;
        }
    }

    var prominences: [H]f32 = undefined;
    var left_bases: [H]i32 = undefined;
    var right_bases: [H]i32 = undefined;

    for (0..n) |p| {
        const peak = peaks[p];

        const x_peak = x[@intCast(peak)];

        left_bases[p] = peak;
        var left_min = x_peak;
        var ii: i32 = peak;
        while (ii >= 0 and x[@intCast(ii)] <= x_peak) : (ii -= 1) {
            if (x[@intCast(ii)] < left_min) {
                left_min = x[@intCast(ii)];
                left_bases[p] = ii;
            }
        }

        right_bases[p] = peak;
        var right_min = x_peak;
        ii = peak;
        while (ii <= i_max and x[@intCast(ii)] <= x_peak) : (ii += 1) {
            if (x[@intCast(ii)] < right_min) {
                right_min = x[@intCast(ii)];
                right_bases[p] = ii;
            }
        }

        prominences[p] = x_peak - @max(left_min, right_min);
    }

    var largest = Peak{ .width = 0, .height = 0, .idx = -1 };
    var largest_size: f32 = 0.0;

    for (0..n) |p| {
        const lo = left_bases[p];
        const hi = right_bases[p];
        const peak = peaks[p];

        const x_peak = x[@intCast(peak)];
        const height = x_peak - 0.5 * prominences[p];

        var ii: i32 = peak;
        while (lo < ii and height < x[@intCast(ii)])
            ii -= 1;

        var left_ip: f32 = @floatFromInt(ii);
        if (x[@intCast(ii)] < height)
            left_ip += (height - x[@intCast(ii)]) / (x[@intCast(ii + 1)] - x[@intCast(ii)]);

        ii = peak;
        while (ii < hi and height < x[@intCast(ii)])
            ii += 1;

        var right_ip: f32 = @floatFromInt(ii);
        if (x[@intCast(ii)] < height)
            right_ip -= (height - x[@intCast(ii)]) / (x[@intCast(ii - 1)] - x[@intCast(ii)]);

        const width = right_ip - left_ip;
        const size = width * x_peak;

        if (size > largest_size) {
            largest = .{ .idx = peak, .width = width, .height = x_peak };
            largest_size = size;
        }
    }

    return largest;
}

/// Initialization strategies adapted from jaakkopasanen/AutoEq
pub fn init_pk(y: [*]f32, f: [*]const f32, fs: f32, lim_f0: Lim, lim_gain: Lim, lim_q: Lim) Filter {
    _ = fs;

    var rect: [K]f32 = undefined;

    for (0..K) |k_| rect[k_] = @max(y[k_], 0.0);
    const peak = largest_peak(&rect, f, lim_f0);

    for (0..K) |k_| rect[k_] = @max(-y[k_], 0.0);
    const dip = largest_peak(&rect, f, lim_f0);

    const p = if (peak.width * peak.height > dip.width * dip.height) peak else dip;

    const f0 = f[@intCast(p.idx)];
    var gain: f32 = if (p.idx == peak.idx) peak.height else -dip.height;
    const bw = p.width * @log2(f[1] / f[0]);
    const bw_exp2 = @exp2(bw);
    var Q: f32 = @sqrt(bw_exp2) / (bw_exp2 - 1.0);

    _ = limit(&gain, lim_gain);
    _ = limit(&Q, lim_q);

    return .{ .f0 = f0, .gain = gain, .Q = Q };
}

fn init_shelf(filter_type: Type, y: [*]f32, f: [*]const f32, fs: f32, lim_f0: Lim, lim_gain: Lim, lim_q: Lim, reverse: bool) Filter {
    

    var lf0 = lim_f0;
    lf0.lo = @max(lf0.lo, 40.0);
    lf0.hi = @min(lf0.hi, 10000.0);

    var best: f32 = 0.0;
    var best_idx: usize = 0;

    var a: f32 = 0.0;
    for (0..K) |k_| {
        const idx: usize = if (reverse) K - 1 - k_ else k_;
        a += y[idx];

        const avg = @abs(a / @as(f32, @floatFromInt(k_ + 1)));
        if (avg > best) {
            best = avg;
            best_idx = idx;
        }
    }

    var f0 = f[best_idx];
    var Q: f32 = std.math.sqrt1_2;

    _ = limit(&f0, lf0);
    _ = limit(&Q, lim_q);

    var w: [K]f32 = @splat(0.0);
    spectrum(filter_type, f0, 1.0, Q, fs, f, &w);

    var p: f32 = 0.0;
    var c: f32 = 0.0;
    for (0..K) |k_| {
        p += w[k_] * y[k_];
        c += w[k_];
    }

    var gain = p / c;
    _ = limit(&gain, lim_gain);

    return .{ .f0 = f0, .gain = gain, .Q = Q };
}

pub fn init_lsc(y: [*]f32, f: [*]const f32, fs: f32, lim_f0: Lim, lim_gain: Lim, lim_q: Lim) Filter {
    return init_shelf(.lsc, y, f, fs, lim_f0, lim_gain, lim_q, false);
}

pub fn init_hsc(y: [*]f32, f: [*]const f32, fs: f32, lim_f0: Lim, lim_gain: Lim, lim_q: Lim) Filter {
    return init_shelf(.hsc, y, f, fs, lim_f0, lim_gain, lim_q, true);
}

pub fn init_fn(filter_type: Type) *const fn ([*]f32, [*]const f32, f32, Lim, Lim, Lim) Filter {
    return switch (filter_type) {
        .pk => &init_pk,
        .lsc => &init_lsc,
        .hsc => &init_hsc,
        else => &init_pk,
    };
}

test "largest_peak finds obvious peak" {
    // Create K-sized arrays with an obvious peak at index 50
    var data: [K]f32 = @splat(-1.0);
    var freqs: [K]f32 = undefined;
    for (0..K) |i| {
        freqs[i] = 20.0 * @exp(@as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(K - 1)) * @log(1000.0));
        data[i] = -@abs(@as(f32, @floatFromInt(i)) - 50.0);
    }
    data[50] = 10.0;

    const result = largest_peak(&data, &freqs, .{ .lo = freqs[0], .hi = freqs[K - 1] });
    try std.testing.expect(result.idx == 50);
    try std.testing.expect(result.height > 5.0);
    try std.testing.expect(result.width > 0);
}
