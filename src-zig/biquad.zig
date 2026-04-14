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
