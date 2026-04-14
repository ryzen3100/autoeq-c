// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later

const std = @import("std");
const types = @import("types.zig");

pub const Type = types.Type;
pub const Biquad = types.Biquad;
pub const K = types.K;
pub const sq = types.sq;
pub const exp10f_inline = types.exp10f_inline;

pub fn pk(A: f32, cos_w: f32, alpha: f32) Biquad {
    const rA = 1.0 / A;
    return .{
        .b0 = A * alpha + 1.0,
        .db0_dA = alpha,
        .db0_dalpha = A,
        .db0_dcos = 0.0,

        .b1 = -2.0 * cos_w,
        .db1_dA = 0.0,
        .db1_dcos = -2.0,

        .b2 = -A * alpha + 1.0,
        .db2_dA = -alpha,
        .db2_dalpha = -A,
        .db2_dcos = 0.0,

        .a0 = (A + alpha) * rA,
        .da0_dA = -alpha * sq(rA),
        .da0_dalpha = rA,
        .da0_dcos = 0.0,

        .a1 = -2.0 * cos_w,
        .da1_dA = 0.0,
        .da1_dcos = -2.0,

        .a2 = (A - alpha) * rA,
        .da2_dA = alpha * sq(rA),
        .da2_dalpha = -rA,
        .da2_dcos = 0.0,
    };
}

pub fn lsc(A: f32, cos_w: f32, alpha: f32) Biquad {
    const p1 = A + 1.0;
    const m1 = A - 1.0;
    const sqrt_A = @sqrt(A);
    const k = 2 * sqrt_A * alpha;
    const dk_dA = alpha / sqrt_A;
    const dk_dalpha = 2.0 * sqrt_A;

    return .{
        .b0 = A * (-cos_w * m1 + k + p1),
        .db0_dA = A * dk_dA - A * cos_w + A - cos_w * m1 + k + p1,
        .db0_dalpha = A * dk_dalpha,
        .db0_dcos = -A * m1,

        .b1 = 2.0 * A * (-cos_w * p1 + m1),
        .db1_dA = -2.0 * A * cos_w + 2.0 * A - 2.0 * cos_w * p1 + 2.0 * m1,
        .db1_dcos = -2.0 * A * p1,

        .b2 = A * (-cos_w * m1 - k + p1),
        .db2_dA = -A * dk_dA - A * cos_w + A - cos_w * m1 - k + p1,
        .db2_dalpha = -A * dk_dalpha,
        .db2_dcos = -A * m1,

        .a0 = cos_w * m1 + k + p1,
        .da0_dA = dk_dA + cos_w + 1.0,
        .da0_dalpha = dk_dalpha,
        .da0_dcos = m1,

        .a1 = -2.0 * cos_w * p1 - 2.0 * m1,
        .da1_dA = -2.0 * cos_w - 2.0,
        .da1_dcos = -2.0 * p1,

        .a2 = cos_w * m1 - k + p1,
        .da2_dA = -dk_dA + cos_w + 1.0,
        .da2_dalpha = -dk_dalpha,
        .da2_dcos = m1,
    };
}

pub fn hsc(A: f32, cos_w: f32, alpha: f32) Biquad {
    const p1 = A + 1.0;
    const m1 = A - 1.0;
    const sqrt_A = @sqrt(A);
    const k = 2 * sqrt_A * alpha;
    const dk_dA = alpha / sqrt_A;
    const dk_dalpha = 2.0 * sqrt_A;

    return .{
        .b0 = A * (cos_w * m1 + k + p1),
        .db0_dA = A * dk_dA + A * cos_w + A + cos_w * m1 + k + p1,
        .db0_dalpha = A * dk_dalpha,
        .db0_dcos = A * m1,

        .b1 = -2.0 * A * (cos_w * p1 + m1),
        .db1_dA = -2.0 * A * cos_w - 2.0 * A - 2.0 * cos_w * p1 - 2.0 * m1,
        .db1_dcos = -2.0 * A * p1,

        .b2 = A * (cos_w * m1 - k + p1),
        .db2_dA = -A * dk_dA + A * cos_w + A + cos_w * m1 - k + p1,
        .db2_dalpha = -A * dk_dalpha,
        .db2_dcos = A * m1,

        .a0 = -cos_w * m1 + k + p1,
        .da0_dA = dk_dA - cos_w + 1.0,
        .da0_dalpha = dk_dalpha,
        .da0_dcos = -m1,

        .a1 = -2.0 * cos_w * p1 + 2.0 * m1,
        .da1_dA = 2.0 - 2.0 * cos_w,
        .da1_dcos = -2.0 * p1,

        .a2 = -cos_w * m1 - k + p1,
        .da2_dA = -dk_dA - cos_w + 1.0,
        .da2_dalpha = -dk_dalpha,
        .da2_dcos = -m1,
    };
}

pub fn biquad_fn(filter_type: Type) *const fn (f32, f32, f32) Biquad {
    return switch (filter_type) {
        .pk => &pk,
        .lsc => &lsc,
        .hsc => &hsc,
        else => &pk,
    };
}

pub fn spectrum(filter_type: Type, f0: f32, gain: f32, Q: f32, fs: f32, f: [*]const f32, y: [*]f32) void {
    @setFloatMode(.optimized);
    const A = exp10f_inline(gain / 40.0);
    const w0 = 2.0 * std.math.pi / fs * f0;
    const cos_w = @cos(w0);
    const sin_w = @sin(w0);
    const alpha = sin_w * 0.5 / Q;

    const s = biquad_fn(filter_type)(A, cos_w, alpha);

    const b_x0 = sq(s.b0 + s.b1 + s.b2);
    const b_x1 = -4.0 * (s.b0 * s.b1 + 4.0 * s.b0 * s.b2 + s.b1 * s.b2);
    const b_x2 = 16.0 * s.b0 * s.b2;
    const a_x0 = sq(s.a0 + s.a1 + s.a2);
    const a_x1 = -4.0 * (s.a0 * s.a1 + 4.0 * s.a0 * s.a2 + s.a1 * s.a2);
    const a_x2 = 16.0 * s.a0 * s.a2;

    for (0..K) |k_| {
        const phi = sq(@sin(std.math.pi / fs * f[k_]));

        const b_poly = b_x0 + phi * (b_x1 + phi * b_x2);
        const a_poly = a_x0 + phi * (a_x1 + phi * a_x2);

        y[k_] += (10.0 / std.math.ln10) * (@log(b_poly) - @log(a_poly));
    }
}

test "pk returns valid biquad for positive gain" {
    const A: f32 = 2.0;
    const cos_w: f32 = 0.5;
    const alpha: f32 = 0.1;
    const b = pk(A, cos_w, alpha);
    try std.testing.expect(!std.math.isNan(b.b0));
    try std.testing.expect(!std.math.isNan(b.b1));
    try std.testing.expect(!std.math.isNan(b.b2));
    try std.testing.expect(!std.math.isNan(b.a1));
    try std.testing.expect(!std.math.isNan(b.a2));
    try std.testing.expect(!std.math.isInf(b.b0));
    // Positive gain peaking filter: b0 should be positive
    try std.testing.expect(b.b0 > 0);
}

test "lsc returns valid biquad" {
    const A: f32 = 1.5;
    const cos_w: f32 = 0.707;
    const alpha: f32 = 0.05;
    const b = lsc(A, cos_w, alpha);
    try std.testing.expect(!std.math.isNan(b.b0));
    try std.testing.expect(!std.math.isNan(b.a1));
    try std.testing.expect(!std.math.isInf(b.b0));
}

test "hsc returns valid biquad" {
    const A: f32 = 0.8;
    const cos_w: f32 = 0.707;
    const alpha: f32 = 0.05;
    const b = hsc(A, cos_w, alpha);
    try std.testing.expect(!std.math.isNan(b.b0));
    try std.testing.expect(!std.math.isNan(b.a1));
    try std.testing.expect(!std.math.isInf(b.b0));
}

test "spectrum pk at center frequency matches gain" {
    const fs: f32 = 48000.0;
    const f0: f32 = 1000.0;
    const gain: f32 = 6.0;
    const Q: f32 = 1.0;

    var freqs: [K]f32 = undefined;
    for (0..K) |i| {
        freqs[i] = 20.0 * @exp(@as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(K - 1)) * @log(1000.0));
    }
    var y: [K]f32 = @splat(0.0);

    spectrum(.pk, f0, gain, Q, fs, &freqs, &y);

    // Response at center frequency should be close to the gain (in dB)
    // Find the closest grid point to f0
    var best_i: usize = 0;
    var best_diff: f32 = 1e10;
    for (0..K) |i| {
        const diff = @abs(freqs[i] - f0);
        if (diff < best_diff) {
            best_diff = diff;
            best_i = i;
        }
    }
    // At center frequency, response should be approximately gain dB
    try std.testing.expectApproxEqAbs(@as(f32, gain), y[best_i], 0.5);

    // Response far from center should be close to 0 dB
    try std.testing.expect(@abs(y[0]) < 2.0);
}
