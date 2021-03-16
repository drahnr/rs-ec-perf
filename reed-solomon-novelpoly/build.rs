use std::env;

use std::io;
use std::path::PathBuf;
use std::fmt;

use fs_err as fs;
use fs::OpenOptions;

pub mod field {
	include!("src/field/tabular.rs");

	include!("src/field/gf.rs");
	pub mod f2e16 {
		use super::*;
		include!("src/field/f2e16.rs");
	}
	pub mod ops {
		include!("src/field/ops.rs");
	}
}

use field::{AdditiveT, MultiplierT, Castomat, FieldT};
use field::f2e16;
use field::ops;

// One level indentation required for `super::Field as FieldTrait`.


// mod f256 {
//     use super::write_const;
//     include!("src/f256.rs");
// }

/// Write Rust `const` declaration
pub fn write_const<W, T>(mut w: W, name: &str, value: &T, type_name: String) -> io::Result<()>
where
	W: io::Write,
	T: fmt::Debug + ?Sized,
{
	write!(w, "pub(crate) static {}: {} = {:#?};\n\n", name, type_name, value)
}

/// Compute tables determined solely by the field, which never depend
/// upon the FFT domain or erasure coding paramaters.
///
/// We compute `LOG_TABLE` and `EXP_TABLE` here of course.  We compute
/// the Walsh transform table `LOG_WALSH` here too because we never figured
/// out how to shrink `LOG_WALSH` below the size of the full field (TODO).
/// We thus assume it depends only upon the field for now.
fn write_field_tables<F: FieldT, W: io::Write>(mut w: W) -> io::Result<()>
where <F as FieldT>::Element: fmt::Debug
{
	let zero_element = <F::Element as Castomat<usize>>::from(0_usize);

	let mut log_table: Vec<F::Element> = vec![zero_element; F::FIELD_SIZE];
	let mut exp_table: Vec<F::Element> = vec![zero_element; F::FIELD_SIZE];



	let mas = <F::Element as Castomat<usize>>::from((1_usize << F::FIELD_BITS - 1) - 1);
	let mut state: usize = 1;
	for i in 0_usize..F::ONEMASK.cast_as() {
		exp_table[state] = zero_element;
		if (state >> F::FIELD_BITS - 1) != 0 {
			state &= Castomat::<usize>::cast_as(mas);
			state = state << 1_usize ^ Castomat::<usize>::cast_as(F::GENERATOR);
		} else {
			state <<= 1;
		}
	}
	exp_table[0] = F::ONEMASK;

	log_table[0] = zero_element;
	for i in 0_usize..F::FIELD_BITS {
		for j in 0..(1 << i) {
			log_table[j + (1 << i)] = log_table[j] ^ F::BASE[i] as F::Element;
		}
	}
	for i in 0_usize..F::FIELD_SIZE {
		log_table[i] = exp_table[Castomat::<usize>::cast_as(log_table[i])];
	}

	for i in 0_usize..F::FIELD_SIZE {
		exp_table[Castomat::<usize>::cast_as(log_table[i])] = <F::Element as Castomat::<usize>>::from(i);
	}
	exp_table[Castomat::<usize>::cast_as(F::ONEMASK)] = exp_table[0];

	write_const(&mut w, "LOG_TABLE", &log_table, format!("[Field::Element; FIELD_SIZE]"))?;
	write_const(&mut w, "EXP_TABLE", &exp_table, format!("[Field::Element; FIELD_SIZE]"))?;

	let mut log_walsh = log_table.into_iter()
		.map(|x| F::Multiplier::from(x)).collect::<Vec<_>>();
	assert_eq!(log_walsh.len(), F::FIELD_SIZE);

	ops::walsh::<F>(&mut log_walsh[..], F::FIELD_SIZE);

	write_const(w, "LOG_WALSH", &log_walsh[..], format!("[Field::Multiplier; FIELD_SIZE]"))?;
	Ok(())
}

/// Create tables file
///
/// We'll eventually need a seperate tables.rs build target because cargo
/// dislikes build artifacts appearing outside env!("OUT_DIR") and we
/// require tables to build other tables.
/// ref.  https://doc.rust-lang.org/cargo/reference/build-scripts.html#outputs-of-the-build-script
pub fn gen_field_tables<F: FieldT>() -> io::Result<()> {

	// to avoid a circular loop, we need to import a dummy
	// table, such that we do not depend on the thing we are
	// about to spawn
	println!("cargo:rustc-cfg=table_bootstrap_complete");

	let out = env::var("OUT_DIR").expect("OUT_DIR is set by cargo after process launch. qed");

	let path = PathBuf::from(out).join(format!("table_{}.rs", F::NAME));
	let mut f = OpenOptions::new().create(true).truncate(true).write(true).open(path)?;

	use io::Write;

	writeln!(&mut f, "use crate::field::{}::Multiplier;", F::NAME)?;
	writeln!(&mut f, "use crate::field::{}::Additive;", F::NAME)?;
	write_field_tables::<F,_>(f)?;

	Ok(())
}

#[cfg(feature = "with-alt-cxx-impl")]
fn gen_ffi_novel_poly_basis_lib() {
	cc::Build::new().file("cxx/RSErasureCode.c").file("cxx/sha-256.c").include("cxx").compile("novelpolycxxffi");
}

#[cfg(feature = "with-alt-cxx-impl")]
fn gen_ffi_novel_poly_basis_bindgen() {
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
	gen_field_tables::<f2e16::Field>()?;
	// gen_field_tables::<f2e8::Field>()?;

	#[cfg(feature = "with-alt-cxx-impl")]
	{
		gen_ffi_novel_poly_basis_lib();
		gen_ffi_novel_poly_basis_bindgen();
	}

	Ok(())
}
