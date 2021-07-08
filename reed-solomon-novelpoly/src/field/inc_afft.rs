

use static_init::{dynamic};

// use core::ops::Mul;

use super::afft::*;


#[dynamic(0)]
pub static AFFT_TABLES: AfftTables = AfftTables::initalize();

/// Additive FFT and IFFT skew table and formal derivative transform table
/// for computing Reed-Solomon in the "novel polynomial basis".
#[allow(non_snake_case)]
pub struct AfftTables {
    /// Logarithm form of twisted factors used in our additive FFT
    pub skews: [Logarithm; ONEMASK as usize], // skew_multiplier
    /// Factors used in formal derivative, actually all zero if field was constructed correctly.
    pub B: [Logarithm; FIELD_SIZE >> 1],
}


/// Formal derivative of polynomial in tweaked?? basis
#[allow(non_snake_case)]
pub fn tweaked_formal_derivative(codeword: &mut [Additive], n: usize) {
    #[cfg(b_is_not_one)]
    let B = unsafe { &AFFT_TABLES.B };

    // We change nothing when multiplying by b from B.
	#[cfg(b_is_not_one)]
	for i in (0..n).into_iter().step_by(2) {
		let b = Logarithm(ONEMASK) - B[i >> 1];
		codeword[i] = codeword[i].mul(b);
		codeword[i + 1] = codeword[i + 1].mul(b);
	}

	formal_derivative(codeword, n);

	// Again changes nothing by multiplying by b although b differs here.
	#[cfg(b_is_not_one)]
	for i in (0..n).into_iter().step_by(2) {
		let b = B[i >> 1];
		codeword[i] = codeword[i].mul(b);
		codeword[i + 1] = codeword[i + 1].mul(b);
	}
}

/// This test ensure that b can safely be bypassed in tweaked_formal_derivative
// TODO: This is wrong right now, use negation!!
#[cfg(not(b_is_not_one))]
#[allow(non_snake_case)]
#[test]
fn b_is_one() {
    let B = unsafe { &AFFT_TABLES.B };
    fn test_b(b: Logarithm) {
        for x in 0..FIELD_SIZE {
            let x = Additive(x as Elt);
        	assert_eq!(x, x.mul(b));
        }
    }
    let mut old_b = None;
	for i in (0..FIELD_SIZE).into_iter().step_by(256) {
        let b = B[i >> 1];
        if old_b != Some(b) {
            test_b( Logarithm(ONEMASK) - b );
            test_b( b );
            old_b = Some(b);
        }
    }
}


impl AfftField for Additive {
    type Multiplier = Logarithm;

    #[inline(always)]
    fn compute_skew(_depart_no: usize, j: usize, index: usize) -> Option<Self::Multiplier> {
        // let i = (j+index) >> depart_no.trailing_zeros();
        // debug_assert_eq!(i, (j+index)/depart_no);
        // let new_skew = Additive((i - 1) as Elt).to_multiplier();
        let old_skew = unsafe { &AFFT_TABLES }.skews[j + index - 1];
        // Actually this does not yet work, indicating a mistake
        /*
        use core::any::TypeId;
        debug_assert!( TypeId::of::<Elt>() != TypeId::of::<u8>() || new_skew == old_skew );
        */

		// It's reasonale to skip the loop if skew is zero, but doing so with
		// all bits set requires justification.	 (TODO)
		if old_skew.0 != ONEMASK { Some(old_skew) } else { None }
    }
}


impl AfftTables {

    /// Initialize SKEW_FACTOR and B
    #[allow(non_snake_case)]
    fn initalize() -> AfftTables {
        // We cannot yet identify if base has an additive or multiplicative
        // representation, or mybe something else entirely.  (TODO)
    	let mut base: [Elt; FIELD_BITS - 1] = Default::default();

        let mut skews_additive = [Additive(0); ONEMASK as usize];

    	for i in 1..FIELD_BITS {
    		base[i - 1] = 1 << i;
    	}

    	// We construct SKEW_FACTOR in additive form to be \bar{s}_j(omega)
    	// from page 6285 for all omega in the field.
    	for m in 0..(FIELD_BITS - 1) {
    		let step = 1 << (m + 1);
    		skews_additive[(1 << m) - 1] = Additive(0);
    		for i in m..(FIELD_BITS - 1) {
    			let s = 1 << (i + 1);

    			let mut j = (1 << m) - 1;
    			while j < s {
    				// Justified by (5) page 6285, except..
    				// we expect SKEW_FACTOR[j ^ field_base[i]] or similar
    				skews_additive[j + s] = skews_additive[j] ^ Additive(base[i]);
    				j += step;
    			}
    		}

    		// Compute base[m] = ONEMASK - base[m] * EXP[LOG[base[m] ^ 1]]
    		// = ONEMASK - base[m] * (base[m] ^ 1)
    		// TODO: But why?
    		//
    		// let idx = mul_table(base[m], LOG_TABLE[(base[m] ^ 1_u16) as usize]);
    		let idx = Additive(base[m]).mul( Additive(base[m] ^ 1).to_multiplier() );
            // WTF?!?
    		// base[m] = ONEMASK - LOG_TABLE[idx as usize];
            base[m] = ONEMASK - idx.to_multiplier().0;

    		// Compute base[i] = base[i] * EXP[b % ONEMASK]
    		// where b = base[m] + LOG[base[i] ^ 1_u16].
    		// As ONEMASK is the order of the multiplicative grou,
    		// base[i] = base[i] * EXP[base[m]] * (base[i] ^ 1)
    		// TODO: But why?
    		for i in (m + 1)..(FIELD_BITS - 1) {
    			// WTF?!?
    			// let b = LOG_TABLE[(base[i] as u16 ^ 1_u16) as usize] as u32 + base[m] as u32;
    			let b = Additive(base[i] ^ 1).to_multiplier().to_wide() + (base[m] as Wide);
    			let b = b % (ONEMASK as Wide);
    			// base[i] = mul_table(base[i], b as u16);
    			base[i] = Additive(base[i]).mul(Logarithm(b as Elt)).0;
    		}
    	}

    	// Convert skew factors from Additive to Logarithm form
        let mut skews_multiplier = [Logarithm(0); ONEMASK as usize];
    	for i in 0..(ONEMASK as usize) {
    		// SKEW_FACTOR[i] = LOG_TABLE[SKEW_FACTOR[i] as usize];
    		skews_multiplier[i] = skews_additive[i].to_multiplier();
    	}

        let mut B = [Logarithm(0); FIELD_SIZE >> 1];

    	// TODO: How does this alter base?
    	base[0] = ONEMASK - base[0];
    	for i in 1..(FIELD_BITS - 1) {
    		base[i] = ( (
                (ONEMASK as Wide) - (base[i] as Wide) + (base[i - 1] as Wide)
            ) % (ONEMASK as Wide) ) as Elt;
    	}

    	// TODO: What is B anyways?
    	B[0] = Logarithm(0);
    	for i in 0..(FIELD_BITS - 1) {
    		let depart = 1 << i;
    		for j in 0..depart {
    			B[j + depart] = Logarithm( ((
                    B[j].to_wide() + (base[i] as Wide)
                ) % (ONEMASK as Wide)) as Elt);
    		}
    	}

        AfftTables {
            // skews_additive,
            skews: skews_multiplier,
			B,
        }
    }

} // impl AdditiveFFT

#[test]
fn flt_back_and_forth() {
    use rand::prelude::*;
    let mut rng = thread_rng();
    let uni = rand::distributions::Uniform::<Elt>::new_inclusive(0, ONEMASK);

    const N: usize = 128;
    for (k,n) in &[(N/4,N)] {
        let mut data = (0..*n).into_iter().map( |_| Additive(uni.sample(&mut rng)) ).collect::<Vec<Additive>>();
        let expected = data.clone();

        afft(&mut data, *n, *k);

        // make sure something is done
        assert!(data.iter().zip(expected.iter()).filter(|(a, b)| { a != b }).count() > 0);

        
        inverse_afft(&mut data, *n, *k);

        assert_eq!(data, expected);
    }
}

