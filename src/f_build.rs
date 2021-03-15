
/// Compute tables determined solely by the field, which never depend
/// upon the FFT domain or erasure coding paramaters.
///
/// We compute `LOG_TABLE` and `EXP_TABLE` here of course.  We compute
/// the Walsh transform table `LOG_WALSH` here too because we never figured
/// out how to shrink `LOG_WALSH` below the size of the full field (TODO).
/// We thus assume it depends only upon the field for now.
#[allow(dead_code)]
pub(crate) fn write_field_tables(out: std::path::PathBuf) -> std::io::Result<()> {
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
			log_table[j + (1 << i)] = log_table[j] ^ BASE[i];
		}
	}
	for i in 0..FIELD_SIZE {
		log_table[i] = exp_table[log_table[i] as usize];
	}

	for i in 0..FIELD_SIZE {
		exp_table[log_table[i] as usize] = i as GFSymbol;
	}
	exp_table[ONEMASK as usize] = exp_table[0];

	let path = out.join(format!("table_{}.rs", FIELD_NAME));
	let mut w = fs_err::OpenOptions::new().create(true).truncate(true).write(true).open(path) ?;

	write_const(&mut w, "LOG_TABLE", &log_table, "[Elt; FIELD_SIZE]") ?;
	write_const(&mut w, "EXP_TABLE", &exp_table, "[Elt; FIELD_SIZE]") ?;

	// mem_cpy(&mut log_walsh[..], &log_table[..]);
	let log_walsh = log_table.clone();
	let mut log_walsh = unsafe { core::mem::transmute::<_, [Multiplier; FIELD_SIZE]>(log_walsh) };
	log_walsh[0] = Multiplier(0);
	walsh(&mut log_walsh[..], FIELD_SIZE);

	write_const(w, "LOG_WALSH", &log_walsh, "[Multiplier; FIELD_SIZE]") ?;
	Ok(())
}
