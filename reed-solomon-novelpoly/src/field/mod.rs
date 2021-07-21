#![feature(const_generics)]
#![feature(const_evaluatable_checked)]
#![feature(destructuring_assignment)]

use crate::errors::*;
use crate::util::*;

#[macro_use]
pub mod macros;
mod traits;
mod field_util;
mod afft_util;
mod afft_field;

pub use traits::{FieldAdd, FieldMul, TruncateTo, Logarithm, Additive};
pub use field_util::{walsh, gf_mul_bitpoly_reduced};
pub use afft_util::{tweaked_formal_derivative};
pub use afft_field::{AfftField};

pub mod f256;
pub mod f2e16;

#[test]
#[ignore]
fn agreement_f2e16_with_f256() {
    use core::ops::Mul;

    for i in 1..=255 {
        let i_f256 = f256::Additive(i).to_multiplier();
        let i_f2e16 = f2e16::Additive(i as u16).to_multiplier();
        for j in 0..=255 {
            let j_f256 = f256::Additive(j).mul(i_f256);
            let j_f2e16 = f2e16::Additive(j as u16).mul(i_f2e16);
            assert_eq!(j_f256.0 as u16, j_f2e16.0);
        }
    }
}
