use std::convert::TryInto;

use core::ops::Mul;
use crate::{FieldAdd, Logarithm, Additive, TruncateTo, FieldMul};


pub trait AfftField : FieldAdd where [(); Self::FIELD_BYTES]: Sized,
    <Self::Wide as TryInto<<Self as FieldAdd>::Element>>::Error : core::fmt::Debug,
{
    type Multiplier = Logarithm<Self>;
    //#[dynamic(0)]

    fn get_skew(i: usize) -> Logarithm<Self> {
        <Self as FieldAdd>::AFFT_SKEW_TABLE[i]
    }

    fn initialize() -> (Vec::<Logarithm<Self>>, Vec::<Logarithm<Self>>) {
        // We cannot yet identify if base has an additive or multiplicative
        // representation, or mybe something else entirely.  (TODO)
    	let mut base: Vec<Self::Element> = vec![Self::ZERO_ELEMENT; Self::FIELD_BITS - 1];

        let mut skews_additive: Vec<Additive<Self>> = vec![Additive(Self::ZERO_ELEMENT); Self::ONEMASK_USIZE];
    
    	for i in 1..Self::FIELD_BITS {
//            base[i - 1] = Self::ONE_ELEMENT << i;
        }

    	// We construct SKEW_FACTOR in additive form to be \bar{s}_j(omega)
    	// from page 6285 for all omega in the field.
    	for m in 0..(Self::FIELD_BITS - 1) {
    		let step = 1 << (m + 1);
    		skews_additive[(1 << m) - 1] = Additive::<Self>(Self::ZERO_ELEMENT);
    		for i in m..(Self::FIELD_BITS - 1) {
    			let s = 1 << (i + 1);

    			let mut j = (1 << m) - 1;
    			while j < s {
    				// Justified by (5) page 6285, except..
    				// we expect SKEW_FACTOR[j ^ field_base[i]] or similar
    				skews_additive[j + s] = skews_additive[j] ^ Additive::<Self>(base[i]);
    				j += step;
    			}
    		}

    		// Compute base[m] = ONEMASK - base[m] * EXP[LOG[base[m] ^ 1]]
    		// = ONEMASK - base[m] * (base[m] ^ 1)
    		// TODO: But why?
    		//
    		// let idx = mul_table(base[m], LOG_TABLE[(base[m] ^ 1_u16) as usize]);
            let lhs = Additive::<Self>(base[m]);
            let log_value_rhs = Logarithm::<Self>(base[m] ^ Self::ONE_ELEMENT);
    		let idx : Additive::<Self> =  lhs * log_value_rhs;
            // WTF?!?
    		// base[m] = ONEMASK - LOG_TABLE[idx as usize];
            base[m] = Self::ONEMASK - FieldMul::<Self, Logarithm<Self>>::to_multiplier(idx).0;

    		// Compute base[i] = base[i] * EXP[b % ONEMASK]
    		// where b = base[m] + LOG[base[i] ^ 1_u16].
    		// As ONEMASK is the order of the multiplicative grou,
    		// base[i] = base[i] * EXP[base[m]] * (base[i] ^ 1)
    		// TODO: But why?
    		for i in (m + 1)..(Self::FIELD_BITS - 1) {
    			// WTF?!?
    			// let b = LOG_TABLE[(base[i] as u16 ^ 1_u16) as usize] as u32 + base[m] as u32;
                let rhs: Self::Wide = (base[m].into());
                let lhs: Self::Wide = Additive::<Self>(base[i] ^ Self::ONE_ELEMENT).to_multiplier().to_wide();
    			let b : Self::Wide = lhs + rhs;
    			let b : Self::Wide = b % (Self::ONEMASK_WIDE);
                let b : Self::Element  = TruncateTo::<Self>::truncate(b);
    			// base[i] = mul_table(base[i], b as u16);
    			base[i] = Mul::<Logarithm<Self>>::mul(Additive::<Self>(base[i]),Logarithm::<Self>(b)).0;
    		}
    	}

    	// Convert skew factors from Additive to Logarithm form
        let mut skews_multiplier = vec![Logarithm(Self::ZERO_ELEMENT); Self::ONEMASK_USIZE];
    	for i in 0..(Self::ONEMASK_USIZE) {
    		// SKEW_FACTOR[i] = LOG_TABLE[SKEW_FACTOR[i] as usize];
    		skews_multiplier[i] = skews_additive[i].to_multiplier();
    	}

        let mut B = vec![Logarithm(Self::ZERO_ELEMENT); Self::FIELD_SIZE >> 1];

    	// TODO: How does this alter base?
    	base[0] = Self::ONEMASK - base[0];
    	for i in 1..(Self::FIELD_BITS - 1) {
    		base[i] = TruncateTo::<Self>::truncate( (
                (Self::ONEMASK_WIDE) - (base[i].into()) + (base[i - 1].into())
            ) % (Self::ONEMASK_WIDE) );
    	}

    	// TODO: What is B anyways?
    	B[0] = Logarithm::<Self>(Self::ZERO_ELEMENT);
    	for i in 0..(Self::FIELD_BITS - 1) {
    		let depart = 1 << i;
    		for j in 0..depart {
                let wide_base : Self::Wide = (base[i].into());
                let wide_exponent : Self::Wide = (B[j].to_wide() + wide_base) % Self::ONEMASK_WIDE;
                let exponent : Self::Element = TruncateTo::<Self>::truncate(wide_exponent);
    			B[j + depart] = Logarithm::<Self>(exponent);
    		}
    	}

        (skews_multiplier, B)
    }

    #[inline(always)]
    fn compute_skew(_depart_no: usize, j: usize, index: usize) -> Option<Logarithm<Self>> {
        // let i = (j+index) >> depart_no.trailing_zeros();
        // debug_assert_eq!(i, (j+index)/depart_no);
        // let new_skew = Additive((i - 1) as Elt).to_multiplier();
        let old_skew = Self::get_skew(j + index - 1);
        // Actually this does not yet work, indicating a mistake
        /*
        use core::any::TypeId;
        debug_assert!( TypeId::of::<Elt>() != TypeId::of::<u8>() || new_skew == old_skew );
        */

		// It's reasonale to skip the loop if skew is zero, but doing so with
		// all bits set requires justification.	 (TODO)
		if old_skew.0 != Self::ONEMASK { Some(old_skew) } else { None }
    }

    /// TODO: Why this is not part of AfftField trait?
    /// Inverse additive FFT in the "novel polynomial basis"
    fn inverse_afft(data: &mut [Additive<Self>], size: usize, index: usize)
    where <<Self as FieldAdd>::Wide as TryInto<<Self as FieldAdd>::Element>>::Error: core::fmt::Debug,
    [(); Self::FIELD_BYTES]: Sized,
    [(); Self::ONEMASK_USIZE]: Sized,
    [(); Self::FIELD_SIZE >> 1]: Sized
    {
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
                if let Some(skew) = Self::compute_skew(depart_no, j, index) {
                    // Again loop on line 3, except skew should depend upon i aka j in Algorithm 2 (TODO)
                    for i in (j - depart_no)..j {
                        // Line 5, justified by (35) page 6288, but
                        // adding depart_no acts like the r+2^i superscript.
                        data[i] ^= data[i + depart_no] * skew;
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
    fn afft(data: &mut [Additive<Self>], size: usize, index: usize)
    where [u8; Self::FIELD_BYTES]: Sized,
    <Self::Wide as TryInto<<Self as FieldAdd>::Element>>::Error : core::fmt::Debug,
    //[(); F::ONEMASK_USIZE]: Sized,
    //[(); F::FIELD_SIZE >> 1]: Sized
    {
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
                if let Some(skew) = Self::compute_skew(depart_no, j, index) {
                    // Loop on line 5, except skew should depend upon i aka j in Algorithm 1 (TODO)
                    for i in (j - depart_no)..j {
                        // Line 6, explained by (28) page 6287, but
                        // adding depart_no acts like the r+2^i superscript.
                        data[i] ^= data[i + depart_no] * skew;
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
    
}


#[test]
fn flt_back_and_forth() {
    use rand::prelude::*;
    let mut rng = thread_rng();
    let uni = rand::distributions::Uniform::<Elt>::new_inclusive(0, F::ONEMASK);

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

// We want the low rate scheme given in
// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf
// and https://github.com/catid/leopard/blob/master/docs/LowRateDecoder.pdf
// but this code resembles https://github.com/catid/leopard which
// implements the high rate decoder in
// https://github.com/catid/leopard/blob/master/docs/HighRateDecoder.pdf
// We're hunting for the differences and trying to undersrtand the algorithm.


