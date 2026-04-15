// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later

const std = @import("std");
const builtin = @import("builtin");
const tp = @import("types.zig");
const biquad_mod = @import("biquad.zig");
const init_mod = @import("init.zig");
const smooth_mod = @import("smooth.zig");

const K = tp.K;
const max_n = tp.max_n;
const Type = tp.Type;
const Biquad = tp.Biquad;
const Filter = tp.Filter;
const Lim = tp.Lim;
const sq = tp.sq;
const limit = tp.limit;
const exp10f_inline = tp.exp10f_inline;
const spectrum = biquad_mod.spectrum;
const biquad_fn = biquad_mod.biquad_fn;
const init_fn = init_mod.init_fn;
const Smooth = smooth_mod.Smooth;
const preprocess = smooth_mod.preprocess;

pub const max_w = 3 * max_n + 1;

/// 4-wide SIMD vector type for K=384 vectorized reductions (384/4=96 iterations).
const Vec4 = @Vector(4, f32);
const v4s: Vec4 = @splat(0.0);

pub const Wrt = struct {
    v: [max_w]f32 = @splat(0.0),
};

inline fn lf_at(wrt: *Wrt, N: usize, i: usize) *f32 {
    return &wrt.v[0 * N + i];
}

inline fn gain_at(wrt: *Wrt, N: usize, i: usize) *f32 {
    return &wrt.v[1 * N + i];
}

inline fn bw_at(wrt: *Wrt, N: usize, i: usize) *f32 {
    return &wrt.v[2 * N + i];
}

inline fn amp_at(wrt: *Wrt, N: usize) *f32 {
    return &wrt.v[3 * N];
}

inline fn lf_at_const(wrt: *const Wrt, N: usize, i: usize) f32 {
    return wrt.v[0 * N + i];
}

inline fn gain_at_const(wrt: *const Wrt, N: usize, i: usize) f32 {
    return wrt.v[1 * N + i];
}

inline fn bw_at_const(wrt: *const Wrt, N: usize, i: usize) f32 {
    return wrt.v[2 * N + i];
}

inline fn amp_at_const(wrt: *const Wrt, N: usize) f32 {
    return wrt.v[3 * N];
}

fn q_to_bw(Q: f32) f32 {
    return 2.0 / std.math.ln2 * std.math.asinh(0.5 / Q);
}

fn bw_to_q(bw: f32) f32 {
    return 0.5 / std.math.sinh(0.5 * std.math.ln2 * bw);
}

const Consts = struct {
    filter_types: [*]const Type,
    phi: [*]const f32,
    r: [*]const f32,
    fs: f32,
    N: usize,
    opt_amp: bool,
};

fn grad(c: *const Consts, x: *Wrt, g: *Wrt) f32 {
    @setFloatMode(.optimized);
    const N = c.N;
    const rK: f32 = 1.0 / @as(f32, @floatFromInt(K));

    const opt_amp = c.opt_amp;

    var dy_dw0: [max_n][K]f32 = undefined;
    var dy_dgain: [max_n][K]f32 = undefined;
    var dy_dbw: [max_n][K]f32 = undefined;

    var w0_v: [max_n]f32 = undefined;

    var pred: [K]f32 = undefined;
    const pred_init: f32 = if (opt_amp) exp10f_inline(amp_at_const(x, N) / 10.0) else 1.0;
    for (0..K) |k_| pred[k_] = pred_init;

    for (0..N) |n_| {
        const f0 = @exp(lf_at_const(x, N, n_));
        const gain = gain_at_const(x, N, n_);
        const bw = bw_at_const(x, N, n_);

        const A = exp10f_inline(gain / 40.0);
        const w0 = 2.0 * std.math.pi / c.fs * f0;
        const cos_w = @cos(w0);
        const sin_w = @sin(w0);
        const kQ = std.math.sinh(0.5 * std.math.ln2 * bw);
        const alpha = sin_w * kQ;

        w0_v[n_] = w0;

        const s: Biquad = biquad_fn(c.filter_types[n_])(A, cos_w, alpha);

        const dA_dgain = A * std.math.ln10 / 40.0;
        const dalpha_dw0 = cos_w * kQ;
        const dalpha_dbw = sin_w * std.math.cosh(0.5 * std.math.ln2 * bw) * 0.5 * std.math.ln2;
        const dcos_dw0 = -sin_w;

        const b_x0 = sq(s.b0 + s.b1 + s.b2);
        const b_x1 = -4.0 * (s.b0 * s.b1 + 4.0 * s.b0 * s.b2 + s.b1 * s.b2);
        const b_x2 = 16.0 * s.b0 * s.b2;
        const a_x0 = sq(s.a0 + s.a1 + s.a2);
        const a_x1 = -4.0 * (s.a0 * s.a1 + 4.0 * s.a0 * s.a2 + s.a1 * s.a2);
        const a_x2 = 16.0 * s.a0 * s.a2;

        const ba = s.b0 + s.b1 + s.b2;
        const aa = s.a0 + s.a1 + s.a2;

        // Precompute scalar constants used inside the k-loop to avoid
        // redundant work when we vectorise 4 elements at a time.
        const bm_c: f32 = 20.0 / std.math.ln10;
        const am_c: f32 = -20.0 / std.math.ln10;
        const b1_4b2: f32 = s.b1 + 4.0 * s.b2;
        const b0_b2: f32 = s.b0 + s.b2;
        const _4b0_b1: f32 = 4.0 * s.b0 + s.b1;
        const a1_4a2: f32 = s.a1 + 4.0 * s.a2;
        const a0_a2: f32 = s.a0 + s.a2;
        const _4a0_a1: f32 = 4.0 * s.a0 + s.a1;

        // --- Vectorised inner loop: 4 elements at a time ---
        const K4 = K / 4;
        var ki: usize = 0;
        while (ki < K4) : (ki += 1) {
            const base = ki * 4;

            // Load 4 phi values
            const phi_v: Vec4 = .{ c.phi[base], c.phi[base + 1], c.phi[base + 2], c.phi[base + 3] };

            // Evaluate polynomials  (Horner form via vector ops)
            const b_x0_v: Vec4 = @splat(b_x0);
            const b_x1_v: Vec4 = @splat(b_x1);
            const b_x2_v: Vec4 = @splat(b_x2);
            const a_x0_v: Vec4 = @splat(a_x0);
            const a_x1_v: Vec4 = @splat(a_x1);
            const a_x2_v: Vec4 = @splat(a_x2);
            const b_poly_v: Vec4 = b_x0_v + phi_v * (b_x1_v + phi_v * b_x2_v);
            const a_poly_v: Vec4 = a_x0_v + phi_v * (a_x1_v + phi_v * a_x2_v);

            // Update prediction
            const ratio_v: Vec4 = b_poly_v / a_poly_v;
            const pred_in: Vec4 = .{ pred[base], pred[base + 1], pred[base + 2], pred[base + 3] };
            const pred_v: Vec4 = pred_in * ratio_v;
            pred[base] = pred_v[0];
            pred[base + 1] = pred_v[1];
            pred[base + 2] = pred_v[2];
            pred[base + 3] = pred_v[3];

            // Gradient computation (scalar — complex chain rule per-element)
            var j: usize = 0;
            while (j < 4) : (j += 1) {
                const k_ = base + j;
                const phi_k = c.phi[k_];
                const b_poly = b_poly_v[j];
                const a_poly = a_poly_v[j];
                const inv_b_poly = 1.0 / b_poly;
                const inv_a_poly = 1.0 / a_poly;

                const bm = bm_c * inv_b_poly;
                const am = am_c * inv_a_poly;

                const _8phi2 = 8.0 * phi_k * phi_k;
                const _2phi = 2.0 * phi_k;

                const dy_db0 = bm * (ba - _2phi * b1_4b2 + _8phi2 * s.b2);
                const dy_db1 = bm * (ba - _2phi * b0_b2);
                const dy_db2 = bm * (ba - _2phi * _4b0_b1 + _8phi2 * s.b0);

                const dy_da0 = am * (aa - _2phi * a1_4a2 + _8phi2 * s.a2);
                const dy_da1 = am * (aa - _2phi * a0_a2);
                const dy_da2 = am * (aa - _2phi * _4a0_a1 + _8phi2 * s.a0);

                const dy_dA =
                    dy_db0 * s.db0_dA +
                    dy_db1 * s.db1_dA +
                    dy_db2 * s.db2_dA +
                    dy_da0 * s.da0_dA +
                    dy_da1 * s.da1_dA +
                    dy_da2 * s.da2_dA;
                const dy_dalpha =
                    dy_db0 * s.db0_dalpha +
                    dy_db2 * s.db2_dalpha +
                    dy_da0 * s.da0_dalpha +
                    dy_da2 * s.da2_dalpha;
                const dy_dcos =
                    dy_db0 * s.db0_dcos +
                    dy_db1 * s.db1_dcos +
                    dy_db2 * s.db2_dcos +
                    dy_da0 * s.da0_dcos +
                    dy_da1 * s.da1_dcos +
                    dy_da2 * s.da2_dcos;

                dy_dw0[n_][k_] = dy_dalpha * dalpha_dw0 + dy_dcos * dcos_dw0;
                dy_dgain[n_][k_] = dy_dA * dA_dgain;
                dy_dbw[n_][k_] = dy_dalpha * dalpha_dbw;
            }
        }
    }

    var L: f32 = 0.0;
    var dL_dy: [K]f32 = undefined;
    var dL_dy_sum: f32 = 0.0;

    // Vectorised loss/dL_dy computation — 4 elements at a time.
    const log_coeff: f32 = 10.0 / std.math.ln10;
    const log_coeff_v: Vec4 = @splat(log_coeff);
    const two_v: Vec4 = @splat(2.0);
    var L_v: Vec4 = v4s;
    var dL_dy_sum_v: Vec4 = v4s;
    const K4 = K / 4;
    var ki: usize = 0;
    while (ki < K4) : (ki += 1) {
        const base = ki * 4;
        const pred_v: Vec4 = .{ pred[base], pred[base + 1], pred[base + 2], pred[base + 3] };
        const r_v: Vec4 = .{ c.r[base], c.r[base + 1], c.r[base + 2], c.r[base + 3] };
        const d_v: Vec4 = log_coeff_v * @log(pred_v) - r_v;
        L_v += d_v * d_v;
        const dL_v: Vec4 = two_v * d_v;
        dL_dy[base] = dL_v[0];
        dL_dy[base + 1] = dL_v[1];
        dL_dy[base + 2] = dL_v[2];
        dL_dy[base + 3] = dL_v[3];
        dL_dy_sum_v += dL_v;
    }
    L = @reduce(.Add, L_v) * rK;
    dL_dy_sum = @reduce(.Add, dL_dy_sum_v);

    amp_at(g, N).* = if (opt_amp) dL_dy_sum * rK else 0.0;

    // Vectorised dot-product reductions — 4 elements at a time.
    for (0..N) |n_| {
        var glf_v: Vec4 = v4s;
        var ggain_v: Vec4 = v4s;
        var gbw_v: Vec4 = v4s;

        var ki2: usize = 0;
        while (ki2 < K4) : (ki2 += 1) {
            const base = ki2 * 4;
            const dy: Vec4 = .{ dL_dy[base], dL_dy[base + 1], dL_dy[base + 2], dL_dy[base + 3] };
            const dw0: Vec4 = .{ dy_dw0[n_][base], dy_dw0[n_][base + 1], dy_dw0[n_][base + 2], dy_dw0[n_][base + 3] };
            const dg: Vec4 = .{ dy_dgain[n_][base], dy_dgain[n_][base + 1], dy_dgain[n_][base + 2], dy_dgain[n_][base + 3] };
            const db: Vec4 = .{ dy_dbw[n_][base], dy_dbw[n_][base + 1], dy_dbw[n_][base + 2], dy_dbw[n_][base + 3] };
            glf_v += dy * dw0;
            ggain_v += dy * dg;
            gbw_v += dy * db;
        }

        lf_at(g, N, n_).* = @reduce(.Add, glf_v) * rK * w0_v[n_];
        gain_at(g, N, n_).* = @reduce(.Add, ggain_v) * rK;
        bw_at(g, N, n_).* = @reduce(.Add, gbw_v) * rK;
    }

    return L;
}

/// AdaBelief optimizer — https://arxiv.org/abs/2010.07468
const AdaBelief = struct {
    m: Wrt,
    s: Wrt,
    b1: f32, b2: f32,
    b1t: f32, b2t: f32,
    eps: f32, eps_root: f32,
    lr: f32,
    N: usize,
    step_count: i32,

    fn init(N: usize) AdaBelief {
        var self: AdaBelief = .{
            .m = .{},
            .s = .{},
            .b1 = 0.9,
            .b1t = 0.9,
            .b2 = 0.99,
            .b2t = 0.99,
            .eps = 1e-12,
            .eps_root = 1e-12,
            .lr = 0.03,
            .N = N,
            .step_count = 0,
        };
        const w_count = 3 * N + 1;
        @memset(self.m.v[0..w_count], 0.0);
        @memset(self.s.v[0..w_count], 0.0);
        return self;
    }

    fn do_step(self: *AdaBelief, x: *Wrt, g: *Wrt) void {
        const N = self.N;
        const w_count = 3 * N + 1;

        for (0..w_count) |w| {
            self.m.v[w] = self.b1 * self.m.v[w] + (1.0 - self.b1) * g.v[w];
            self.s.v[w] = self.b2 * self.s.v[w] + (1.0 - self.b2) * sq(g.v[w] - self.m.v[w]);

            const m_hat = self.m.v[w] / (1.0 - self.b1t);
            const s_hat = self.s.v[w] / (1.0 - self.b2t);

            const den = @sqrt(s_hat + self.eps_root) + self.eps;

            x.v[w] -= self.lr * m_hat / den;
        }

        self.b1t *= self.b1;
        self.b2t *= self.b2;
        self.step_count += 1;
    }
};

fn info(comptime fmt: [:0]const u8, args: anytype) void {
    if (builtin.target.os.tag == .freestanding) return;
    @call(.never_inline, info_impl, .{ fmt, args });
}

fn info_impl(comptime fmt: [:0]const u8, args: anytype) void {
    std.debug.print(fmt, args);
}

pub fn fit(
    steps: i32,
    filter_types: [*]const Type,
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
    @setFloatMode(.optimized);
    const LOG: i32 = 250;

    var lf_lim: [max_n]Lim = undefined;
    var bw_lim: [max_n]Lim = undefined;

    for (0..N) |n_| {
        lf_lim[n_] = .{
            .lo = @log(f0_lim[n_].lo),
            .hi = @log(f0_lim[n_].hi),
        };
        bw_lim[n_] = .{
            .lo = q_to_bw(Q_lim[n_].hi),
            .hi = q_to_bw(Q_lim[n_].lo),
        };
    }

    var phi: [K]f32 = undefined;
    for (0..K) |k_| phi[k_] = sq(@sin(std.math.pi / fs * f[k_]));

    var x: Wrt = .{};

    for (0..N) |n_| {
        lf_at(&x, N, n_).* = @log(f0[n_]);
        gain_at(&x, N, n_).* = gain[n_];
        bw_at(&x, N, n_).* = q_to_bw(Q[n_]);
    }
    amp_at(&x, N).* = if (amp) |a| a.* else 0.0;

    info("\n", .{});
    if (amp) |a|
        info("{s:>12}: {d: >34.2} dB\n", .{ "P", @as(f64, a.*) });
    for (0..N) |n_| {
        info("{d:>12}: {s:>8} {d:>12.1} Hz {d:>9.2} dB {d:>10.3} Q\n", .{
            n_ + 1,
            tp.type_names[@intFromEnum(filter_types[n_])],
            @as(f64, f0[n_]),
            @as(f64, gain[n_]),
            @as(f64, Q[n_]),
        });
    }
    info("\n", .{});

    var g: Wrt = .{};
    var best: Wrt = .{};
    var L: f32 = 0.0;
    var best_L: f32 = 1e9;

    var t0: i128 = if (builtin.target.os.tag != .freestanding) std.time.nanoTimestamp() else 0;

    const c = Consts{
        .filter_types = filter_types,
        .phi = &phi,
        .r = r,
        .fs = fs,
        .N = N,
        .opt_amp = amp != null,
    };

    var opt = AdaBelief.init(N);

    const w_count = 3 * N + 1;
    var step_idx: i32 = 0;
    while (step_idx < steps) : (step_idx += 1) {
        L = grad(&c, &x, &g);

        opt.do_step(&x, &g);

        for (0..N) |n_| {
            if (limit(lf_at(&x, N, n_), lf_lim[n_]))
                lf_at(&opt.m, N, n_).* = 0.0;

            if (limit(gain_at(&x, N, n_), gain_lim[n_]))
                gain_at(&opt.m, N, n_).* = 0.0;

            if (limit(bw_at(&x, N, n_), bw_lim[n_]))
                bw_at(&opt.m, N, n_).* = 0.0;
        }

        if (L < best_L) {
            best_L = L;
            @memcpy(best.v[0..w_count], x.v[0..w_count]);
        }

        if (@rem(step_idx, LOG) == 0 and step_idx != 0) {
            if (builtin.target.os.tag != .freestanding) {
                const t1 = std.time.nanoTimestamp();
                const elapsed = @as(f64, @floatFromInt(t1 - t0)) / 1e9;
                info("L: {d:.5} | {d:>6}/{d} [{d:>6.0} it/s]\n", .{
                    L,
                    step_idx,
                    steps,
                    @as(f64, @floatFromInt(LOG)) / elapsed,
                });
                t0 = t1;
            }
        }
    }

    for (0..N) |n_| {
        f0[n_] = @exp(lf_at_const(&best, N, n_));
        gain[n_] = gain_at_const(&best, N, n_);
        Q[n_] = bw_to_q(bw_at_const(&best, N, n_));
    }

    if (amp) |a|
        a.* = amp_at_const(&best, N);

    info("L: {d:.5}\n\n", .{best_L});
    if (amp) |a|
        info("{s:>12}: {d: >34.2} dB\n", .{ "P", @as(f64, a.*) });
    for (0..N) |n_| {
        info("{d:>12}: {s:>8} {d:>12.1} Hz {d:>9.2} dB {d:>10.3} Q\n", .{
            n_ + 1,
            tp.type_names[@intFromEnum(filter_types[n_])],
            @as(f64, f0[n_]),
            @as(f64, gain[n_]),
            @as(f64, Q[n_]),
        });
    }
    info("\n", .{});

    return best_L;
}

pub fn autoeq(
    steps: i32,
    filter_types: [*]const Type,
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
    if (steps <= 0) return 0.0;

    var r_init: [K]f32 = undefined;
    @memcpy(&r_init, r[0..K]);

    for (0..N) |n_| {
        const filter_type = filter_types[n_];
        const p: Filter = init_fn(filter_type)(&r_init, f, fs, f0_lim[n_], gain_lim[n_], Q_lim[n_]);
        spectrum(filter_type, p.f0, -p.gain, p.Q, fs, f, &r_init);

        f0[n_] = p.f0;
        gain[n_] = p.gain;
        Q[n_] = p.Q;
    }

    if (amp) |a|
        a.* = 0.0;

    return fit(steps, filter_types, f0, gain, Q, amp, f0_lim, gain_lim, Q_lim, N, f, r, fs);
}

test "grad returns positive loss and finite gradients" {
    const N: usize = 2;
    var freqs: [K]f32 = undefined;
    for (0..K) |i| {
        freqs[i] = 20.0 * @exp(@as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(K - 1)) * @log(1000.0));
    }
    var r: [K]f32 = @splat(1.0);

    var types: [max_n]Type = @splat(.pk);
    types[0] = .lsc;
    types[1] = .hsc;

    // phi must be sin²(π f / fs) — same as computed in fit()
    var phi: [K]f32 = undefined;
    for (0..K) |i| {
        phi[i] = sq(@sin(std.math.pi / 48000.0 * freqs[i]));
    }

    const consts = Consts{
        .filter_types = &types,
        .phi = &phi,
        .r = &r,
        .fs = 48000.0,
        .N = @intCast(N),
        .opt_amp = true,
    };

    var x = Wrt{};
    // Set some initial values
    x.v[0] = 7.0; // lf0 for filter 0
    x.v[1] = 6.0; // lf0 for filter 1
    x.v[2] = 0.0; // gain for filter 0
    x.v[3] = 0.0; // gain for filter 1
    x.v[4] = 0.5; // bw for filter 0
    x.v[5] = 0.5; // bw for filter 1
    x.v[6] = 0.0; // amp

    var g = Wrt{};
    const loss = grad(&consts, &x, &g);
    try std.testing.expect(loss > 0);
    try std.testing.expect(std.math.isFinite(loss));

    const w = 3 * N + 1;
    for (0..w) |i| {
        try std.testing.expect(std.math.isFinite(g.v[i]));
    }
}

test "autoeq end-to-end with minimal config" {
    var freqs: [K]f32 = undefined;
    for (0..K) |i| {
        freqs[i] = 20.0 * @exp(@as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(K - 1)) * @log(1000.0));
    }
    var r: [K]f32 = @splat(1.0);
    var f0: [max_n]f32 = @splat(1000.0);
    var gain: [max_n]f32 = @splat(0.0);
    var Q: [max_n]f32 = @splat(1.0);
    var amp: f32 = 0.0;

    var types: [max_n]Type = @splat(.pk);
    types[0] = .lsc;
    types[1] = .hsc;

    var f0_lim: [max_n]Lim = @splat(.{ .lo = 20.0, .hi = 16000.0 });
    var gain_lim: [max_n]Lim = @splat(.{ .lo = -16.0, .hi = 16.0 });
    var q_lim: [max_n]Lim = @splat(.{ .lo = 0.4, .hi = 4.0 });

    const loss = autoeq(
        50, &types, &f0, &gain, &Q, &amp,
        &f0_lim, &gain_lim, &q_lim,
        2, &freqs, &r, 48000.0,
    );

    try std.testing.expect(std.math.isFinite(loss));
    try std.testing.expect(loss >= 0);
    // Filter parameters should be within bounds
    for (0..2) |i| {
        try std.testing.expect(f0[i] >= 20.0);
        try std.testing.expect(f0[i] <= 16000.0);
        try std.testing.expect(gain[i] >= -16.0);
        try std.testing.expect(gain[i] <= 16.0);
        try std.testing.expect(Q[i] >= 0.4);
        try std.testing.expect(Q[i] <= 4.0);
    }
}
