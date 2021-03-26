use std::io;
use std::env;
use std::fmt;
use std::path::PathBuf;
use fs_err as fs;

/// Write Rust `const` declaration
pub fn write_static<W, T>(mut w: W, name: &str, value: &T, type_name: &str) -> io::Result<()>
where
	W: io::Write,
	T: fmt::Debug,
{
	writeln!(w, r###"#[allow(unused)]
pub(crate) static {}: {} = {:#?};"###, name, type_name, value)
}


/// Compute tables determined solely by the field, which never depend
/// upon the FFT domain or erasure coding paramaters.
///
/// We compute `LOG_TABLE` and `EXP_TABLE` here of course.  We compute
/// the Walsh transform table `LOG_WALSH` here too because we never figured
/// out how to shrink `LOG_WALSH` below the size of the full field (TODO).
/// We thus assume it depends only upon the field for now.
#[allow(unused)]
fn write_field_tables<W: io::Write>(mut w: W) -> io::Result<()> {
    let mut base: [Elt; FIELD_BITS] = generate_cantor_basis(BASE_FINAL).unwrap(
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

	let mut log_table: [Elt; FIELD_SIZE] = [0; FIELD_SIZE];
	let mut exp_table: [Elt; FIELD_SIZE] = [0; FIELD_SIZE];

	let mas: Elt = (1 << FIELD_BITS - 1) - 1;
	let mut state: usize = 1;
	for i in 0_usize..(ONEMASK as usize) {
		exp_table[state] = i as Elt;
		if (state >> FIELD_BITS - 1) != 0 {
			state &= mas as usize;
			state = state << 1_usize ^ GENERATOR as usize;
		} else {
			state <<= 1;
		}
	}
	exp_table[0] = ONEMASK;

	log_table[0] = 0;
	for i in 0..FIELD_BITS {
		for j in 0..(1 << i) {
			log_table[j + (1 << i)] = log_table[j] ^ base[i];
		}
	}
	for i in 0..FIELD_SIZE {
		log_table[i] = exp_table[log_table[i] as usize];
	}

	for i in 0..FIELD_SIZE {
		exp_table[log_table[i] as usize] = i as Elt;
	}
	exp_table[ONEMASK as usize] = exp_table[0];

	write_static(&mut w, "BASE", &base, "[Elt; FIELD_BITS]") ?;

	write_static(&mut w, "LOG_TABLE", &log_table, "[Elt; FIELD_SIZE]") ?;
	write_static(&mut w, "EXP_TABLE", &exp_table, "[Elt; FIELD_SIZE]") ?;

	// mem_cpy(&mut log_walsh[..], &log_table[..]);
	let log_walsh = log_table.clone();
	let mut log_walsh = unsafe { core::mem::transmute::<_, [Logarithm; FIELD_SIZE]>(log_walsh) };
	log_walsh[0] = Logarithm(0);
	walsh(&mut log_walsh[..], FIELD_SIZE);

	write_static(w, "LOG_WALSH", &log_walsh, "[Logarithm; FIELD_SIZE]")?;
	Ok(())
}

/// Create tables file
///
/// We'll eventually need a seperate tables.rs build target because cargo
/// dislikes build artifacts appearing outside env!("OUT_DIR") and we
/// require tables to build other tables.
/// ref.  https://doc.rust-lang.org/cargo/reference/build-scripts.html#outputs-of-the-build-script
pub fn gen_field_tables() -> io::Result<()> {
	let out = env::var("OUT_DIR").expect("OUT_DIR is set by cargo after process launch. qed");
	let path = PathBuf::from(out).join(format!("table_{}.rs", FIELD_NAME));
	let f = fs::OpenOptions::new().create(true).truncate(true).write(true).open(path)?;
	write_field_tables(f)
}
