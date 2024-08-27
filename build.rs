use std::path::PathBuf;
use std::env;
use std::process::Command;


fn main() {
    
    let make_result = Command::new("make")
        .arg("libabpoa")
        .current_dir("abPOA")
        .status()
        .expect("Could not build abPOA! GNU make installed?");
    
    if !make_result.success() {
        panic!("Build of abPOA unsuccessful!");
    }
    
    println!("cargo:rustc-link-search=native={}/abPOA/lib", env::current_dir().unwrap().display());
    println!("cargo:rustc-link-lib=abpoa");
    
    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("abPOA/include/abpoa.h")
        .clang_args([
            "-Wall",
            "-Wno-unused-function",
            "-Wno-misleading-indentation",
            "-DUSE_SIMDE",
            "-DSIMDE_ENABLE_NATIVE_ALIASES",
        ])
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