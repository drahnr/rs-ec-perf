use std::io;
use std::env;
use std::fmt;
use std::path::PathBuf;
use fs_err as fs;
use core::convert::{TryFrom,TryInto};

/// Write Rust `const` declaration
pub fn write_static<W, T>(mut w: W, name: &str, value: &T, type_name: &str) -> io::Result<()>
where
	W: io::Write,
	T: fmt::Debug,
{
	writeln!(w, r###"#[allow(unused)]
     const {}: {} = {:#?};"###, name, type_name, value)

//	writeln!(w, r###"#[allow(unused)]
//pub(crate) static {}: {} = {:#?};"###, name, type_name, value)
}


/// Compute tables determined solely by the field, which never depend
/// upon the FFT domain or erasure coding paramaters.
///
/// We compute `LOG_TABLE` and `EXP_TABLE` here of course.  We compute
/// the Walsh transform table `LOG_WALSH` here too because we never figured
/// out how to shrink `LOG_WALSH` below the size of the full field (TODO).
/// We thus assume it depends only upon the field for now.
#[allow(unused)]

fn write_field_tables<W: io::Write, F: FieldAdd>(mut w: W) -> io::Result<()>
where
    <F::Wide as TryInto<F::Element>>::Error: core::fmt::Debug,
<F::Element as TryFrom<usize>>::Error: core::fmt::Debug,
[F::Element; F::FIELD_BITS] : Sized,
[F::Element; F::FIELD_SIZE] : Sized,
{
    let mut base: [F::Element; F::FIELD_BITS] = F::generate_cantor_basis(F::BASE_FINAL).unwrap(
    );
    // Find new Cantor basis if desired, but requires correct multiplication
    /*
    unwrap_or_else(|| {
        for bf in (0..=255).rev() {
            if let Some(base) = generate_cantor_basis(bf) { return base; }
        }
        panic!();
    });
    */

	let mut log_table: [F::Element; F::FIELD_SIZE] = [F::ZERO_ELEMENT; F::FIELD_SIZE];
	let mut exp_table: [F::Element; F::FIELD_SIZE] = [F::ZERO_ELEMENT; F::FIELD_SIZE];

	let mas: F::Element = (F::ONE_ELEMENT << F::FIELD_BITS - 1) - F::ONE_ELEMENT;
	let mut state: usize = 1;
	for i in 0_usize..(F::ONEMASK.into()) {
		exp_table[state] = i.try_into().unwrap();
		if (state >> F::FIELD_BITS - 1) != 0 {
            let mas_usize : usize = mas.into();
			state  &= mas_usize;
            let gen_usize : usize = F::GENERATOR.into();
			state = state << 1_usize ^ gen_usize;
		} else {
			state <<= 1;
		}
	}
	exp_table[0] = F::ONEMASK;

	log_table[0] = F::ZERO_ELEMENT;
	for i in 0..F::FIELD_BITS {
		for j in 0..(1 << i) {
			log_table[j + (1 << i)] = log_table[j] ^ base[i];
		}
	}
	for i in 0..F::FIELD_SIZE {
		log_table[i] = exp_table[Into::<usize>::into(log_table[i])];
	}

	for i in 0..F::FIELD_SIZE {
		exp_table[Into::<usize>::into(log_table[i])] = i.try_into().unwrap();
	}
	exp_table[Into::<usize>::into(F::ONEMASK)] = exp_table[0];

	// write_static(&mut w, "BASE", &base, ["[",F::FIELD_NAME,"::Element; ",F::FIELD_NAME,"::FIELD_BITS]"].concat().as_str()) ?;

	// write_static(&mut w, "LOG_TABLE", &log_table, ["[",F::FIELD_NAME,"::Element; ",F::FIELD_NAME,"::FIELD_SIZE]"].concat().as_str()) ?;
	// write_static(&mut w, "EXP_TABLE", &exp_table, ["[",F::FIELD_NAME,"::Element; ",F::FIELD_NAME,"::FIELD_SIZE]"].concat().as_str()) ?;
	write_static(&mut w, "BASE", &base, "[Self::Element; Self::FIELD_BITS]") ?;

	write_static(&mut w, "LOG_TABLE", &log_table, "[Self::Element; Self::FIELD_SIZE]") ?;
	write_static(&mut w, "EXP_TABLE", &exp_table, "[Self::Element; Self::FIELD_SIZE]") ?;

	// mem_cpy(&mut log_walsh[..], &log_table[..]);
	let mut log_walsh = log_table.clone().iter_mut().map(|element| Logarithm::<F>(*element)).collect::<Vec<Logarithm<F>>>();
	
	log_walsh[0] = Logarithm(F::ZERO_ELEMENT);
	walsh(&mut log_walsh[..], F::FIELD_SIZE);

	write_static(w, "LOG_WALSH", &log_walsh, "[Logarithm<Self>; Self::FIELD_SIZE]")?;
	// write_static(w, "LOG_WALSH", &log_walsh, ["[Logarithm::<",F::FIELD_NAME,">; ",F::FIELD_NAME,"::FIELD_SIZE]"].concat().as_str())?;
	Ok(())
}

/// Create tables file
///
/// We'll eventually need a seperate tables.rs build target because cargo
/// dislikes build artifacts appearing outside env!("OUT_DIR") and we
/// require tables to build other tables.
/// ref.  https://doc.rust-lang.org/cargo/reference/build-scripts.html#outputs-of-the-build-script
    pub fn gen_field_tables<F: FieldAdd>() -> io::Result<()> where
                           <F::Wide as TryInto<F::Element>>::Error: core::fmt::Debug,
                           <F::Element as TryFrom<usize>>::Error: core::fmt::Debug,
                           [F::Element; F::FIELD_BITS] : Sized,
                           [F::Element; F::FIELD_SIZE] : Sized

                            {
	let out = env::var("OUT_DIR").expect("OUT_DIR is set by cargo after process launch. qed");
	let path = PathBuf::from(out).join(format!("table_{}.rs", F::FIELD_NAME));
	let f = fs::OpenOptions::new().create(true).truncate(true).write(true).open(path)?;
	write_field_tables::<_, F>(f) 
}
