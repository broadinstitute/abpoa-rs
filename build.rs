use std::path::PathBuf;
use std::env;

use glob::glob;


fn main() {
    let mut compiler_args = vec![
        "-Wall",
        "-Wno-unused-function",
        "-Wno-misleading-indentation",
        "-DUSE_SIMDE",
        "-DSIMDE_ENABLE_NATIVE_ALIASES",
    ];

    if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
        compiler_args.push("-D__AVX2__");
        compiler_args.push("-march=armv8-a+simd");
    }

    let mut builder = cc::Build::new();

    builder
        .files(
            glob("abPOA/src/*.c")
                .expect("Could not find abPOA source files")
                .map(|x| x.unwrap())
                .filter(|x| !x.ends_with("abpoa.c")) // exclude the CLI main file (which also includes a main() function)
        )
        .include("abPOA/include");

    for flag in &compiler_args {
        builder.flag(flag);
    }

    builder.compile("abpoa");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("abPOA/include/abpoa.h")
        .header("abPOA/src/abpoa_seq.h")
        .clang_args(&compiler_args)
        .opaque_type("abpoa_simd_matrix_t")
        .blocklist_type(r".*\d+x\d+_t")
        .blocklist_type(r".*\d+x\d+x\d+_t")
        .blocklist_type(r"neon_.*")
        .blocklist_type(r"simde_.*")
        .blocklist_type(r"__m\d+[dfi]?")
        .blocklist_type(r"SIMD[if]")
        .blocklist_item(r"F[EP]_.*")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap()).join("bindings.rs");
    bindings
        .write_to_file(out_path)
        .expect("Couldn't write bindings!");
}
