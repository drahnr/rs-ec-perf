use super::{FieldAdd, FieldMul};

/// Formal derivative of polynomial in the new?? basis
pub fn formal_derivative<F: AfftField>(cos: &mut [F], size: usize) {
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

/// ..
pub trait AfftField: FieldMul<Self::Multiplier> {
    type Multiplier: Clone + Copy;
    fn compute_skew(depart_no: usize, j: usize, index: usize) -> Option<Self::Multiplier>;
}

// We want the low rate scheme given in
// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf
// and https://github.com/catid/leopard/blob/master/docs/LowRateDecoder.pdf
// but this code resembles https://github.com/catid/leopard which
// implements the high rate decoder in
// https://github.com/catid/leopard/blob/master/docs/HighRateDecoder.pdf
// We're hunting for the differences and trying to undersrtand the algorithm.

/// Inverse additive FFT in the "novel polynomial basis"
pub fn inverse_afft<Additive: AfftField>(data: &mut [Additive], size: usize, index: usize) {
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
            if let Some(skew) = Additive::compute_skew(depart_no, j, index) {
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
pub fn afft<Additive: AfftField>(data: &mut [Additive], size: usize, index: usize) {
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
            if let Some(skew) = Additive::compute_skew(depart_no, j, index) {
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
