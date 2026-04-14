// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later

const std = @import("std");
const root = @import("autoeq");

const K = root.K;
const max_n = root.max_n;
const Type = root.Type;
const Lim = root.Lim;
const spectrum = root.spectrum;
const Smooth = root.Smooth;
const IE_SMOOTH = root.IE_SMOOTH;
const OE_SMOOTH = root.OE_SMOOTH;
const preprocess = root.preprocess;
const autoeq = root.autoeq;
const type_names = root.type_names;

const Spec = struct {
    type: Type,
    f0: Lim,
    gain: Lim,
    q: Lim,
};

fn make_spec(comptime t: Type, comptime f: [2]f32, comptime g: [2]f32, comptime q: [2]f32) Spec {
    return .{
        .type = t,
        .f0 = .{ .lo = f[0], .hi = f[1] },
        .gain = .{ .lo = g[0], .hi = g[1] },
        .q = .{ .lo = q[0], .hi = q[1] },
    };
}

const SPECS: [max_n]Spec = blk: {
    var specs: [max_n]Spec = undefined;
    specs[0] = make_spec(.lsc, .{ 20.0, 16000.0 }, .{ -16.0, 16.0 }, .{ 0.4, 3.0 });
    specs[1] = make_spec(.hsc, .{ 20.0, 16000.0 }, .{ -16.0, 16.0 }, .{ 0.4, 3.0 });
    for (2..max_n) |i| {
        specs[i] = make_spec(.pk, .{ 20.0, 16000.0 }, .{ -16.0, 16.0 }, .{ 0.4, 4.0 });
    }
    break :blk specs;
};

comptime {
    std.debug.assert(SPECS.len == max_n);
}

const TYPES: [max_n]Type = blk: {
    var types: [max_n]Type = undefined;
    for (0..max_n) |i| types[i] = SPECS[i].type;
    break :blk types;
};

const F0_LIM: [max_n]Lim = blk: {
    var lims: [max_n]Lim = undefined;
    for (0..max_n) |i| lims[i] = SPECS[i].f0;
    break :blk lims;
};

const GAIN_LIM: [max_n]Lim = blk: {
    var lims: [max_n]Lim = undefined;
    for (0..max_n) |i| lims[i] = SPECS[i].gain;
    break :blk lims;
};

const Q_LIM: [max_n]Lim = blk: {
    var lims: [max_n]Lim = undefined;
    for (0..max_n) |i| lims[i] = SPECS[i].q;
    break :blk lims;
};

pub fn main() !void {
    const t0 = std.time.nanoTimestamp();

    const stdout = std.fs.File.stdout();
    const stdin = std.fs.File.stdin();

    const args = try std.process.argsAlloc(std.heap.page_allocator);
    defer std.process.argsFree(std.heap.page_allocator, args);

    if (args.len < 3) {
        std.debug.print("usage: autoeq [iox] <N> <?steps>\n", .{});
        std.process.exit(1);
    }

    const smooth_opt: ?*const Smooth = switch (args[1][0]) {
        'i' => &IE_SMOOTH,
        'o' => &OE_SMOOTH,
        'x' => null,
        else => {
            std.debug.print("invalid type\n", .{});
            std.process.exit(1);
        },
    };

    const N_parsed = std.fmt.parseInt(i32, args[2], 10) catch {
        std.debug.print("invalid count\n", .{});
        std.process.exit(1);
    };
    if (N_parsed <= 0 or N_parsed > max_n) {
        std.debug.print("invalid count\n", .{});
        std.process.exit(1);
    }
    const N_u: usize = @intCast(N_parsed);

    var steps: i32 = 3000;
    if (args.len == 4) {
        steps = std.fmt.parseInt(i32, args[3], 10) catch {
            std.debug.print("invalid steps\n", .{});
            std.process.exit(1);
        };
        if (steps <= 0) {
            std.debug.print("invalid steps\n", .{});
            std.process.exit(1);
        }
    }

    const FS: f32 = 48000.0;
    const F0: f32 = 20.0;
    const F1: f32 = 20000.0;
    const L0 = @log(F0);
    const L1 = @log(F1);
    const LR = L1 - L0;

    var freq: [K]f32 = undefined;
    for (0..K) |k_| freq[k_] = @exp(L0 + LR / @as(f32, @floatFromInt(K - 1)) * @as(f32, @floatFromInt(k_)));

    var dst: [K]f32 = undefined;
    var src: [K]f32 = undefined;
    read_stdin(&dst, stdin) catch {
        std.debug.print("error reading dst from stdin\n", .{});
        std.process.exit(1);
    };
    read_stdin(&src, stdin) catch {
        std.debug.print("error reading src from stdin\n", .{});
        std.process.exit(1);
    };

    var gain: [max_n]f32 = undefined;
    var f0: [max_n]f32 = undefined;
    var q: [max_n]f32 = undefined;
    var amp: f32 = 0.0;

    var r: [K]f32 = undefined;
    _ = preprocess(&freq, &dst, &src, &r, smooth_opt, true);

    _ = autoeq(steps, &TYPES, &f0, &gain, &q, &amp, &F0_LIM, &GAIN_LIM, &Q_LIM, N_u, &freq, &r, FS);

    var y: [K]f32 = @splat(0.0);
    for (0..N_u) |n_|
        spectrum(TYPES[n_], f0[n_], gain[n_], q[n_], FS, &freq, &y);

    var max_val: f32 = -1e9;
    for (0..K) |k_| max_val = @max(max_val, y[k_]);

    // Output — match C format exactly
    var buf: [1024]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();

    w.print("{e:.7} {e:.7}\n", .{ amp, max_val }) catch unreachable;
    for (0..N_u) |n_| {
        w.print("{s} {e:.7} {e:.7} {e:.7}\n", .{
            type_names[@intFromEnum(TYPES[n_])],
            f0[n_],
            gain[n_],
            q[n_],
        }) catch unreachable;
    }
    stdout.writeAll(fbs.getWritten()) catch unreachable;

    const t1 = std.time.nanoTimestamp();
    const elapsed_ms = @as(f64, @floatFromInt(t1 - t0)) / 1e6;
    std.debug.print("\ntime: {d:.0} ms\n\n", .{elapsed_ms});
}

fn read_stdin(out: *[K]f32, file: std.fs.File) !void {
    const bytes = std.mem.sliceAsBytes(out);
    var offset: usize = 0;
    while (offset < bytes.len) {
        const n = try file.read(bytes[offset..]);
        if (n == 0) return error.UnexpectedEof;
        offset += n;
    }
}
