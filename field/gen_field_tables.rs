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

	write_static(&mut w, "BASE", &base, ["[<",F::FIELD_NAME," as FieldAdd>::Element; <",F::FIELD_NAME," as FieldAdd>::FIELD_BITS]"].concat().as_str()) ?;

	write_static(&mut w, "LOG_TABLE", &log_table, ["[<",F::FIELD_NAME," as FieldAdd>::Element; <",F::FIELD_NAME," as FieldAdd>::FIELD_SIZE]"].concat().as_str()) ?;
	write_static(&mut w, "EXP_TABLE", &exp_table, ["[<",F::FIELD_NAME," as FieldAdd>::Element; <",F::FIELD_NAME," as FieldAdd>::FIELD_SIZE]"].concat().as_str()) ?;

	// write_static(&mut w, "BASE", &base, "[Self::Element; Self::FIELD_BITS]") ?;

	// write_static(&mut w, "LOG_TABLE", &log_table, "[Self::Element; Self::FIELD_SIZE]") ?;
	// write_static(&mut w, "EXP_TABLE", &exp_table, "[Self::Element; Self::FIELD_SIZE]") ?;

	// mem_cpy(&mut log_walsh[..], &log_table[..]);
	let mut log_walsh = log_table.clone().iter_mut().map(|element| Logarithm::<F>(*element)).collect::<Vec<Logarithm<F>>>();
	
	log_walsh[0] = Logarithm(F::ZERO_ELEMENT);
	walsh(&mut log_walsh[..], F::FIELD_SIZE);

	write_static(&mut w, "LOG_WALSH", &log_walsh, ["[Logarithm<",F::FIELD_NAME,">; <",F::FIELD_NAME," as FieldAdd>::FIELD_SIZE]"].concat().as_str())?;
	                                       // write_static(w, "LOG_WALSH", &log_walsh, ["[Logarithm::<<",F::FIELD_NAME," as FieldAdd>>; <",F::FIELD_NAME," as FieldAdd>::FIELD_SIZE]"].concat().as_str())?;

   let (afft_skew_table, b) = initialize_afft_table::<F>(&log_table, &exp_table);
	write_static(&mut w, "AFFT_SKEW_TABLE", &afft_skew_table, ["[Logarithm<",F::FIELD_NAME,">; <",F::FIELD_NAME," as FieldAdd>::ONEMASK_USIZE]"].concat().as_str())?;
	                                       // write_static(w, "LOG_WALSH", &log_walsh, ["[Logarithm::<<",F::FIELD_NAME," as FieldAdd>>; <",F::FIELD_NAME," as FieldAdd>::F
//                                           AFFT_SKEW_TABLE                                           
	Ok(())
}

fn initialize_afft_table<F: FieldAdd>(log_table: &[F::Element], exp_table: &[F::Element]) -> (Vec::<Logarithm<F>>, Vec::<Logarithm<F>>) 
where  <F::Wide as TryInto<F::Element>>::Error: core::fmt::Debug,
{                        
    // We cannot yet identify if base has an additive or multiplicative
    // representation, or mybe something else entirely.  (TODO)
    let mut base: Vec<F::Element> = vec![F::ZERO_ELEMENT; F::FIELD_BITS - 1];

        let mut skews_additive: Vec<Additive<F>> = vec![Additive(F::ZERO_ELEMENT); F::ONEMASK_USIZE];
    
    	for i in 1..F::FIELD_BITS {
//            base[i - 1] = F::ONE_ELEMENT << i;
        }

    	// We construct SKEW_FACTOR in additive form to be \bar{s}_j(omega)
    	// from page 6285 for all omega in the field.
    	for m in 0..(F::FIELD_BITS - 1) {
    		let step = 1 << (m + 1);
    		skews_additive[(1 << m) - 1] = Additive::<F>(F::ZERO_ELEMENT);
    		for i in m..(F::FIELD_BITS - 1) {
    			let s = 1 << (i + 1);

    			let mut j = (1 << m) - 1;
    			while j < s {
    				// Justified by (5) page 6285, except..
    				// we expect SKEW_FACTOR[j ^ field_base[i]] or similar
    				skews_additive[j + s] = skews_additive[j] ^ Additive::<F>(base[i]);
    				j += step;
    			}
    		}

    		// Compute base[m] = ONEMASK - base[m] * EXP[LOG[base[m] ^ 1]]
    		// = ONEMASK - base[m] * (base[m] ^ 1)
    		// TODO: But why?
    		//
    		// let idx = mul_table(base[m], LOG_TABLE[(base[m] ^ 1_u16) as usize]);
            let log_value_rhs = log_table[Into::<usize>::into(base[m] ^ F::ONE_ELEMENT)];
    		let idx = log_table[Into::<usize>::into(base[m])] + log_value_rhs;
            // WTF?!?
    		// base[m] = ONEMASK - LOG_TABLE[idx as usize];
            base[m] = F::ONEMASK - idx;

    		// Compute base[i] = base[i] * EXP[b % ONEMASK]
    		// where b = base[m] + LOG[base[i] ^ 1_u16].
    		// As ONEMASK is the order of the multiplicative grou,
    		// base[i] = base[i] * EXP[base[m]] * (base[i] ^ 1)
    		// TODO: But why?
    		for i in (m + 1)..(F::FIELD_BITS - 1) {
    			// WTF?!?
    			// let b = LOG_TABLE[(base[i] as u16 ^ 1_u16) as usize] as u32 + base[m] as u32;
                let rhs: F::Wide = (base[m].into());
                let lhs: F::Wide = From::<F::Element>::from(log_table[Into::<usize>::into(base[i] ^ F::ONE_ELEMENT)]);
    			let b : F::Wide = lhs + rhs;
    			let b : F::Wide = b % (F::ONEMASK_WIDE);
                let b : F::Element  = TruncateTo::<F>::truncate(b);
    			// base[i] = mul_table(base[i], b as u16);
    			base[i] = exp_table[Into::<usize>::into(log_table[Into::<usize>::into(base[i])] + b)];
    		}
    	}

    	// Convert skew factors from Additive to Logarithm form
        let mut skews_multiplier = vec![Logarithm(F::ZERO_ELEMENT); F::ONEMASK_USIZE];
    	for i in 0..(F::ONEMASK_USIZE) {
    		// SKEW_FACTOR[i] = LOG_TABLE[SKEW_FACTOR[i] as usize];
    		skews_multiplier[i] = Logarithm::<F>(log_table[Into::<usize>::into(skews_additive[i].0)]);
    	}

        let mut b = vec![Logarithm(F::ZERO_ELEMENT); F::FIELD_SIZE >> 1];

    	// TODO: How does this alter base?
    	base[0] = F::ONEMASK - base[0];
    	for i in 1..(F::FIELD_BITS - 1) {
    		base[i] = TruncateTo::<F>::truncate( (
                (F::ONEMASK_WIDE) - (base[i].into()) + (base[i - 1].into())
            ) % (F::ONEMASK_WIDE) );
    	}

    	// TODO: What is b anyways?
    	b[0] = Logarithm::<F>(F::ZERO_ELEMENT);
    	for i in 0..(F::FIELD_BITS - 1) {
    		let depart = 1 << i;
    		for j in 0..depart {
                let wide_base : F::Wide = (base[i].into());
                let wide_exponent : F::Wide = (b[j].to_wide() + wide_base) % F::ONEMASK_WIDE;
                let exponent : F::Element = TruncateTo::<F>::truncate(wide_exponent);
    			b[j + depart] = Logarithm::<F>(exponent);
    		}
    	}

        (skews_multiplier, b)
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

                           
