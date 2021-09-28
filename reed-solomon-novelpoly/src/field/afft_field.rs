use std::convert::{TryInto};
use crate::{FieldAdd, Logarithm, Additive};


pub trait AfftField : FieldAdd where [(); Self::FIELD_BYTES]: Sized,
    <Self::Wide as TryInto<<Self as FieldAdd>::Element>>::Error : core::fmt::Debug,
{
    type Multiplier = Logarithm<Self>;
    //#[dynamic(0)]

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

use crate::f256::{F256};
use crate::f2e16::F2e16;

#[cfg(test)]
mod tests {

   use super::*;
    fn flt_back_and_forth<F: AfftField>()
    where
        [u8; F::FIELD_BYTES]: Sized,
    [(); F::FIELD_BYTES]: Sized,
    [(); F::ONEMASK_USIZE]: Sized,
    [(); F::FIELD_SIZE >> 1]: Sized,
    <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug,
    <F::Element as TryFrom<usize>>::Error : core::fmt::Debug

    {
        use rand::prelude::*;
        let mut rng = thread_rng();
        let uni = rand::distributions::Uniform::<usize>::new_inclusive(0, F::ONEMASK_USIZE);
        
        const N: usize = 128;
        for (k,n) in &[(N/4,N)] {
            let mut data = (0..*n).into_iter().map( |_| Additive::<F>(<F::Element as TryFrom::<usize>>::try_from(uni.sample(&mut rng)).unwrap()) ).collect::<Vec<Additive::<F>>>();
            let expected = data.clone();
            
 F::afft(&mut data, *n, *k);
            
            // make sure something is done
            assert!(data.iter().zip(expected.iter()).filter(|(a, b)| { a != b }).count() > 0);
            
            
            F::inverse_afft(&mut data, *n, *k);
            
            assert_eq!(data, expected);
        }
    }
    
    
    test_all_fields_for!(flt_back_and_forth);


// // We want the low rate scheme given in
// // https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf
// // and https://github.com/catid/leopard/blob/master/docs/LowRateDecoder.pdf
// // but this code resembles https://github.com/catid/leopard which
// // implements the high rate decoder in
// // https://github.com/catid/leopard/blob/master/docs/HighRateDecoder.pdf
// // We're hunting for the differences and trying to undersrtand the algorithm.


}
