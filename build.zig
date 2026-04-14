const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Native CLI binary — root module is cli.zig, with "root" import = root.zig
    const cli_mod = b.createModule(.{
        .root_source_file = b.path("src-zig/cli.zig"),
        .target = target,
        .optimize = optimize,
    });

    // The "root" import chain: root.zig -> types.zig, biquad.zig, init.zig, smooth.zig, optimize.zig
    // Since root.zig uses @import("types.zig") etc., these resolve as relative file imports
    // within the same module. We need to register root.zig as a named import.
    const root_mod = b.createModule(.{
        .root_source_file = b.path("src-zig/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    cli_mod.addImport("autoeq", root_mod);

    const exe = b.addExecutable(.{
        .name = "autoeq",
        .root_module = cli_mod,
    });
    b.installArtifact(exe);

    // WASM module
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });
    const wasm_root_mod = b.createModule(.{
        .root_source_file = b.path("src-zig/root.zig"),
        .target = wasm_target,
        .optimize = .ReleaseFast,
    });
    const wasm_mod = b.createModule(.{
        .root_source_file = b.path("src-zig/wasm.zig"),
        .target = wasm_target,
        .optimize = .ReleaseFast,
    });
    wasm_mod.addImport("autoeq", wasm_root_mod);

    const wasm_exe = b.addExecutable(.{
        .name = "autoeq-wasm",
        .root_module = wasm_mod,
    });
    wasm_exe.entry = .disabled;
    wasm_exe.rdynamic = true;
    b.installArtifact(wasm_exe);

    // Tests
    const test_mod = b.createModule(.{
        .root_source_file = b.path("src-zig/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    const lib_tests = b.addTest(.{
        .root_module = test_mod,
    });
    const run_tests = b.addRunArtifact(lib_tests);
    b.step("test", "Run tests").dependOn(&run_tests.step);
}
