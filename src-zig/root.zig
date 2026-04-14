// Copyright (C) 2026 PEQdB Inc.
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// Aggregation module — re-exports everything from sub-modules.

const std = @import("std");

// Re-export types
pub const max_n: usize = 32;
pub const K: usize = 384;
pub const Type = @import("types.zig").Type;
pub const type_names = @import("types.zig").type_names;
pub const Biquad = @import("types.zig").Biquad;
pub const Filter = @import("types.zig").Filter;
pub const Lim = @import("types.zig").Lim;
pub const clip = @import("types.zig").clip;
pub const sq = @import("types.zig").sq;
pub const exp10f_inline = @import("types.zig").exp10f_inline;
pub const limit = @import("types.zig").limit;

// Re-export sub-modules
pub const biquad_mod = @import("biquad.zig");
pub const init_mod = @import("init.zig");
pub const smooth_mod = @import("smooth.zig");
pub const optimize_mod = @import("optimize.zig");

// Convenience re-exports
pub const spectrum = biquad_mod.spectrum;
pub const preprocess = smooth_mod.preprocess;
pub const autoeq = optimize_mod.autoeq;
pub const Smooth = smooth_mod.Smooth;
pub const IE_SMOOTH = smooth_mod.IE_SMOOTH;
pub const OE_SMOOTH = smooth_mod.OE_SMOOTH;

test {
    _ = @import("types.zig");
    _ = biquad_mod;
    _ = init_mod;
    _ = smooth_mod;
    _ = optimize_mod;
}
