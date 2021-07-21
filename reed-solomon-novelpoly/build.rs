#![allow(unused)]
#![feature(const_generics)]
#![feature(const_evaluatable_checked)]
#![feature(associated_type_defaults)]
#![feature(destructuring_assignment)]

use std::io;

include!("src/util.rs");
include!("src/field/traits.rs");
include!("src/field/macros.rs");
include!("src/field/field_util.rs");
include!("gen_field_tables.rs");

mod f2e16 {
    include!("src/field/f2e16.rs");
}

mod f256 {
    include!("src/field/f256.rs");
}

#[cfg(feature = "with-alt-cxx-impl")]
fn gen_ffi_novel_poly_basis_lib() {
    cc::Build::new().file("cxx/RSErasureCode.c").include("cxx").compile("novelpolycxxffi");
}

#[cfg(feature = "with-alt-cxx-impl")]
fn gen_ffi_novel_poly_basis_bindgen() {
    use std::env;
    use std::path::PathBuf;

    println!("cargo:rustc-link-lib=novelpolycxxffi");

    println!("cargo:rerun-if-changed=wrapper.h");

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings.write_to_file(out_path.join("bindings.rs")).expect("Couldn't write bindings!");
}

fn main() -> io::Result<()> {    
    gen_field_tables::<f2e16::F2e16>()?;
    gen_field_tables::<f256::F256>()?;
    // to avoid a circular loop, we need to import a dummy
    // table, such that we do not depend on the thing we are
    // about to spawn
    println!("cargo:rustc-cfg=table_bootstrap_complete");

    #[cfg(feature = "with-alt-cxx-impl")]
    {
        gen_ffi_novel_poly_basis_lib();
        gen_ffi_novel_poly_basis_bindgen();
    }

    Ok(())
}
