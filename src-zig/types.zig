// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// Shared types, constants, and helpers.
// This file has NO imports from sibling files to avoid circular deps.

const std = @import("std");

pub const max_n: usize = 32;
pub const K: usize = 384;

pub const Type = enum(u32) {
    pk = 0,
    lsc = 1,
    hsc = 2,
    _,
};

pub const type_names = [_][]const u8{ "PK", "LSC", "HSC" };

pub const Biquad = extern struct {
    b0: f32, b1: f32, b2: f32, a0: f32, a1: f32, a2: f32,

    db0_dA: f32, db0_dalpha: f32, db0_dcos: f32,
    db1_dA: f32, db1_dcos: f32,
    db2_dA: f32, db2_dalpha: f32, db2_dcos: f32,
    da0_dA: f32, da0_dalpha: f32, da0_dcos: f32,
    da1_dA: f32, da1_dcos: f32,
    da2_dA: f32, da2_dalpha: f32, da2_dcos: f32,
};

pub const Filter = extern struct {
    f0: f32,
    gain: f32,
    Q: f32,
};

pub const Lim = extern struct {
    lo: f32,
    hi: f32,
};

pub fn clip(x: f32, lo: f32, hi: f32) f32 {
    return @max(@min(x, hi), lo);
}

pub fn sq(x: f32) f32 {
    return x * x;
}

pub fn exp10f_inline(x: f32) f32 {
    return @exp(std.math.ln10 * x);
}

pub fn limit(x: *f32, lim: Lim) bool {
    const orig = x.*;
    x.* = clip(x.*, lim.lo, lim.hi);
    return x.* != orig;
}
test "clip within range" {
    try std.testing.expectEqual(@as(f32, 5.0), clip(5.0, 0.0, 10.0));
    try std.testing.expectEqual(@as(f32, -3.0), clip(-3.0, -5.0, 5.0));
}

test "clip at boundaries" {
    try std.testing.expectEqual(@as(f32, 0.0), clip(0.0, 0.0, 10.0));
    try std.testing.expectEqual(@as(f32, 10.0), clip(10.0, 0.0, 10.0));
}

test "clip outside range" {
    try std.testing.expectEqual(@as(f32, 0.0), clip(-5.0, 0.0, 10.0));
    try std.testing.expectEqual(@as(f32, 10.0), clip(15.0, 0.0, 10.0));
}

test "sq" {
    try std.testing.expectEqual(@as(f32, 0.0), sq(0.0));
    try std.testing.expectEqual(@as(f32, 4.0), sq(2.0));
    try std.testing.expectEqual(@as(f32, 9.0), sq(-3.0));
}

test "exp10f_inline" {
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), exp10f_inline(0.0), 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 10.0), exp10f_inline(1.0), 1e-4);
    try std.testing.expectApproxEqAbs(@as(f32, 100.0), exp10f_inline(2.0), 1e-3);
    try std.testing.expectApproxEqAbs(@as(f32, 0.1), exp10f_inline(-1.0), 1e-5);
}

test "limit clamps and reports change" {
    var val1: f32 = 5.0;
    const changed1 = limit(&val1, .{ .lo = 0.0, .hi = 10.0 });
    try std.testing.expectEqual(@as(f32, 5.0), val1);
    try std.testing.expect(!changed1);

    var val2: f32 = -5.0;
    const changed2 = limit(&val2, .{ .lo = 0.0, .hi = 10.0 });
    try std.testing.expectEqual(@as(f32, 0.0), val2);
    try std.testing.expect(changed2);

    var val3: f32 = 15.0;
    const changed3 = limit(&val3, .{ .lo = 0.0, .hi = 10.0 });
    try std.testing.expectEqual(@as(f32, 10.0), val3);
    try std.testing.expect(changed3);
}

