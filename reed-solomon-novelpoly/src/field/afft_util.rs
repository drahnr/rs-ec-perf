use super::{FieldAdd, FieldMul, Logarithm, Additive};
use core::convert::{TryFrom,TryInto, Into};

use static_init::{dynamic};

/// Formal derivative of polynomial in tweaked?? basis
#[allow(non_snake_case)]
pub fn tweaked_formal_derivative<F: FieldAdd>(codeword: &mut [Additive<F>], n: usize) {
    #[cfg(b_is_not_one)]
    let B = unsafe { &AFFT_TABLES.B };

    // We change nothing when multiplying by b from B.
	#[cfg(b_is_not_one)]
	for i in (0..n).into_iter().step_by(2) {
		let b = Logarithm::<F>(F::ONEMASK) - B[i >> 1];
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
// TODO: Provide B again and activate test
////////////////////////////////////////////////////////////////
// #[test]                                                    //
// fn b_is_one() {                                            //
//     let B = unsafe { &AFFT_TABLES.B };                     //
//     fn test_b(b: Logarithm) {                              //
//         for x in 0..FIELD_SIZE {                           //
//             let x = Additive(x as Elt);                    //
//         	assert_eq!(x, x.mul(b));                          //
//         }                                                  //
//     }                                                      //
//     let mut old_b = None;                                  //
// 	for i in (0..FIELD_SIZE).into_iter().step_by(256) {       //
//         let b = B[i >> 1];                                 //
//         if old_b != Some(b) {                              //
//             test_b( Logarithm(ONEMASK) - b );              //
//             test_b( b );                                   //
//             old_b = Some(b);                               //
//         }                                                  //
//     }                                                      //
// }                                                          //
////////////////////////////////////////////////////////////////


/// Formal derivative of polynomial in the new?? basis
pub fn formal_derivative<F: FieldAdd>(cos: &mut [Additive<F>], size: usize) {
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


