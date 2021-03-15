use static_init::dynamic;

#[dynamic(0)]
pub(crate) static AFFT: AdditiveFFT<super::field::f2e16::Field> = AdditiveFFT::<F>::initalize();


/// Additive FFT and inverse in the "novel polynomial basis"
pub(crate) struct AdditiveFFT<F> {
    /// Multiplier form of twisted factors used in AdditiveFFT
    pub skews: [Multiplier<F>; F::ONEMASK.cast_as()], // skew_multiplier
}


/// Formal derivative of polynomial in the new?? basis
pub(crate) fn formal_derivative<F>(cos: &mut [Additive<F>], size: usize) {
	for i in 1..size {
		let length = ((i ^ i - 1) + 1) >> 1;
		for j in (i - length)..i {
			cos[j] ^= cos.get(j + length).copied().unwrap_or_default();
		}
	}
	let mut i = size;
	while i < F::FIELD_SIZE && i < cos.len() {
		for j in 0..size {
			cos[j] ^= cos.get(j + i).copied().unwrap_or_default();
		}
		i <<= 1;
	}
}


// We want the low rate scheme given in
// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf
// and https://github.com/catid/leopard/blob/master/docs/LowRateDecoder.pdf
// but this code resembles https://github.com/catid/leopard which
// implements the high rate decoder in
// https://github.com/catid/leopard/blob/master/docs/HighRateDecoder.pdf
// We're hunting for the differences and trying to undersrtand the algorithm.

/// Inverse additive FFT in the "novel polynomial basis"
pub(crate) fn inverse_afft<F: FieldT>(data: &mut [Additive<F>], size: usize, index: usize) {
    unsafe { &AFFT }.inverse_afft(data,size,index)
}

/// Additive FFT in the "novel polynomial basis"
pub(crate) fn afft<F: FieldT>(data: &mut [Additive<F>], size: usize, index: usize) {
    unsafe { &AFFT }.afft(data,size,index)
}


impl<F: Field> AdditiveFFT<F> {

    /// Inverse additive FFT in the "novel polynomial basis"
    pub(crate) fn inverse_afft(&self, data: &mut [Additive<F>], size: usize, index: usize) {
    	// All line references to Algorithm 2 page 6288 of
    	// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf

    	// Depth of the recursion on line 7 and 8 is given by depart_no
    	// aka 1 << ((k of Algorithm 2) - (i of Algorithm 2)) where
    	// k of Algorithm 1 is read as FIELD_BITS here.
    	// Recusion base layer implicitly imports d_r aka ala line 1.
    	// After this, we start at depth (i of Algorithm 2) = (k of Algorithm 2) - 1
    	// and progress through FIELD_BITS-1 steps, obtaining \Psi_\beta(0,0).
    	let mut depart_no = 1_usize;
    	while depart_no < size {
    		// Agrees with for loop (j of Algorithm 2) in (0..2^{k-i-1}) from line 3,
    		// except we've j in (depart_no..size).step_by(2*depart_no), meaning
    		// the doubled step compensated for the halve size exponent, and
    		// somehow this j captures the subscript on \omega_{j 2^{i+1}}.	 (TODO)
    		let mut j = depart_no;
    		while j < size {
    			// At this point loops over i in (j - depart_no)..j give a bredth
    			// first loop across the recursion branches from lines 7 and 8,
    			// so the i loop corresponds to r in Algorithm 2.  In fact,
    			// data[i] and data[i + depart_no] together cover everything,
    			// thanks to the outer j loop.

    			// Loop on line 3, so i corresponds to j in Algorithm 2
    			for i in (j - depart_no)..j {
    				// Line 4, justified by (34) page 6288, but
    				// adding depart_no acts like the r+2^i superscript.
    				data[i + depart_no] ^= data[i];
    			}

    			// Algorithm 2 indexs the skew factor in line 5 page 6288
    			// by i and \omega_{j 2^{i+1}}, but not by r explicitly.
    			// We further explore this confusion below. (TODO)
    			let skew = self.skews[j + index - 1];
    			// It's reasonale to skip the loop if skew is zero, but doing so with
    			// all bits set requires justification.	 (TODO)
    			if skew.0 != F::ONEMASK {
    				// Again loop on line 3, except skew should depend upon i aka j in Algorithm 2 (TODO)
    				for i in (j - depart_no)..j {
    					// Line 5, justified by (35) page 6288, but
    					// adding depart_no acts like the r+2^i superscript.
    					data[i] ^= data[i + depart_no].mul(skew);
    				}
    			}

    			// Increment by double depart_no in agreement with
    			// our updating 2*depart_no elements at this depth.
    			j += depart_no << 1;
    		}
    		depart_no <<= 1;
    	}
    }

    /// Additive FFT in the "novel polynomial basis"
    pub(crate) fn afft(&self, data: &mut [Additive<F>], size: usize, index: usize) {
    	// All line references to Algorithm 1 page 6287 of
    	// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf

    	// Depth of the recursion on line 3 and 4 is given by depart_no
    	// aka 1 << ((k of Algorithm 1) - (i of Algorithm 1)) where
    	// k of Algorithm 1 is read as FIELD_BITS here.
    	// Recusion base layer implicitly imports d_r aka ala line 1.
    	// After this, we start at depth (i of Algorithm 1) = (k of Algorithm 1) - 1
    	// and progress through FIELD_BITS-1 steps, obtaining \Psi_\beta(0,0).
    	let mut depart_no = size >> 1_usize;
    	while depart_no > 0 {
    		// Agrees with for loop (j of Algorithm 1) in (0..2^{k-i-1}) from line 5,
    		// except we've j in (depart_no..size).step_by(2*depart_no), meaning
    		// the doubled step compensated for the halve size exponent, and
    		// somehow this j captures the subscript on \omega_{j 2^{i+1}}.	 (TODO)
    		let mut j = depart_no;
    		while j < size {
    			// At this point loops over i in (j - depart_no)..j give a bredth
    			// first loop across the recursion branches from lines 3 and 4,
    			// so the i loop corresponds to r in Algorithm 1.  In fact,
    			// data[i] and data[i + depart_no] together cover everything,
    			// thanks to the outer j loop.

    			// Algorithm 1 indexs the skew factor in line 6 aka (28) page 6287
    			// by i and \omega_{j 2^{i+1}}, but not by r explicitly.
    			// We doubt the lack of explicit dependence upon r justifies
    			// extracting the skew factor outside the loop here.
    			// As indexing by \omega_{j 2^{i+1}} appears absolute elsewhere,
    			// we think r actually appears but the skew factor repeats itself
    			// like in (19) in the proof of Lemma 4.  (TODO)
    			// We should understand the rest of this basis story, like (8) too.	 (TODO)
    			let skew = self.skews[j + index - 1];
    			// It's reasonale to skip the loop if skew is zero, but doing so with
    			// all bits set requires justification.	 (TODO)
    			if skew.0 != ONEMASK {
    				// Loop on line 5, except skew should depend upon i aka j in Algorithm 1 (TODO)
    				for i in (j - depart_no)..j {
    					// Line 6, explained by (28) page 6287, but
    					// adding depart_no acts like the r+2^i superscript.
    					data[i] ^= data[i + depart_no].mul(skew);
    				}
    			}

    			// Again loop on line 5, so i corresponds to j in Algorithm 1
    			for i in (j - depart_no)..j {
    				// Line 7, explained by (31) page 6287, but
    				// adding depart_no acts like the r+2^i superscript.
    				data[i + depart_no] ^= data[i];
    			}

    			// Increment by double depart_no in agreement with
    			// our updating 2*depart_no elements at this depth.
    			j += depart_no << 1;
    		}
    		depart_no >>= 1;
    	}
    }


    // initialize SKEW_FACTOR and B
    fn initalize() -> AdditiveFFT<F> {
        // We cannot yet identify if base has an additive or multiplicative
        // representation, or mybe something else entirely.  (TODO)
    	let mut base: [F::Element; F::FIELD_BITS - 1] = Default::default();

        let mut skews_additive = [Additive<F>(0); F::ONEMASK.cast_as()];

    	for i in 1..F::FIELD_BITS {
    		base[i - 1] = 1 << i;
    	}

    	// We construct SKEW_FACTOR in additive form to be \bar{s}_j(omega)
    	// from page 6285 for all omega in the field.
    	for m in 0..(F::FIELD_BITS - 1) {
    		let step = 1 << (m + 1);
    		skews_additive[(1 << m) - 1] = Additive<F>(0);
    		for i in m..(F::FIELD_BITS - 1) {
    			let s = 1 << (i + 1);

    			let mut j = (1 << m) - 1;
    			while j < s {
    				// Justified by (5) page 6285, except..
    				// we expect SKEW_FACTOR[j ^ field_base[i]] or similar
    				skews_additive[j + s] = skews_additive[j] ^ Additive<F>(base[i]);
    				j += step;
    			}
    		}

    		// Compute base[m] = ONEMASK - base[m] * EXP[LOG[base[m] ^ 1]]
    		// = ONEMASK - base[m] * (base[m] ^ 1)
    		// TODO: But why?
    		//
    		// let idx = mul_table(base[m], LOG_TABLE[(base[m] ^ 1_u16).cast_as()]);
    		let idx = Additive<F>(base[m]).mul( Additive<F>(base[m] ^ 1).to_multiplier() );
            // WTF?!?
    		// base[m] = ONEMASK - LOG_TABLE[idx.cast_as()];
            base[m] = F::ONEMASK - idx.to_multiplier().0;

    		// Compute base[i] = base[i] * EXP[b % ONEMASK]
    		// where b = base[m] + LOG[base[i] ^ 1_u16].
    		// As ONEMASK is the order of the multiplicative grou,
    		// base[i] = base[i] * EXP[base[m]] * (base[i] ^ 1)
    		// TODO: But why?
    		for i in (m + 1)..(FIELD_BITS - 1) {
    			// WTF?!?
    			// let b = LOG_TABLE[(base[i] as u16 ^ 1_u16).cast_as()] as u32 + base[m] as u32;
    			let b = Additive<F>(base[i] ^ 1).to_multiplier().to_wide() + (base[m] as Wide);
    			let b = b % (F::ONEMASK as Wide);
    			// base[i] = mul_table(base[i], b as u16);
    			base[i] = Additive<F>(base[i]).mul(Multiplier<F>(b as Elt)).0;
    		}
    	}

    	// Convert skew factors from Additive to Multiplier form
        let mut skews_multiplier = [Multiplier<F>(0); F::ONEMASK.cast_as()];
    	for i in 0..(ONEMASK.cast_as()) {
    		// SKEW_FACTOR[i] = LOG_TABLE[SKEW_FACTOR[i].cast_as()];
    		skews_multiplier[i] = skews_additive[i].to_multiplier();
    	}

        Self {
            // skews_additive,
            skews: skews_multiplier,
        }
    }

}
