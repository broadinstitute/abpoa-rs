use std::path::PathBuf;
use std::env;

use glob::glob;


fn main() {
    let mut builder = cc::Build::new();
        
    builder
        .files(
            glob("abPOA/src/*.c")
            .expect("Could not find abPOA source files")
            .map(|x| x.unwrap())
        )
        .include("abPOA/include")
        .define("USE_SIMDE", None)
        .define("SIMDE_ENABLE_NATIVE_ALIASES", None)
        .flag("-Wall")
        .flag("-Wno-unused-function")
        .flag("-Wno-misleading-indentation");
        
    if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
        builder.define("__AVX2__", None);
        builder.flag("-march=armv8-a+simd");
    }
    
    builder.compile("abpoa");
    
    let compiler_args = [
        "-Wall",
        "-Wno-unused-function",
        "-Wno-misleading-indentation",
        "-DUSE_SIMDE",
        "-DSIMDE_ENABLE_NATIVE_ALIASES",
    ];
    
    
    
    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("abPOA/include/abpoa.h")
        .clang_args(&compiler_args)
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