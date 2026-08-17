fn main() {
    println!("cargo:rerun-if-changed=../src/gpu/knearest.cpp");
    println!("cargo:rerun-if-changed=../src/gpu/capi.cpp");
    println!("cargo:rerun-if-changed=../include/linkcell_gpu.h");
    println!("cargo:rerun-if-changed=../include/linkcell_gpu.hpp");
    println!("cargo:rerun-if-changed=../include/linkcell.hpp");
    println!("cargo:rerun-if-changed=../include/linkcell.h");
    println!("cargo:rerun-if-changed=../third_party/gpulite/gpulite/gpulite.hpp");

    let manifest = std::path::PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let root = manifest.parent().unwrap();

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++17")
        .define("LINKCELL_HAS_GPULITE", "1")
        .include(root.join("include"))
        .include(root.join("third_party/gpulite"))
        .file(root.join("src/gpu/knearest.cpp"))
        .file(root.join("src/gpu/capi.cpp"))
        .warnings(false);
    if cfg!(target_os = "linux") || cfg!(target_os = "macos") {
        build.flag_if_supported("-fPIC");
    }
    build.compile("linkcell_gpu");

    println!("cargo:rustc-link-lib=dylib=dl");
    println!("cargo:rustc-link-lib=dylib=pthread");
    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=dylib=c++");
    } else {
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }
}
