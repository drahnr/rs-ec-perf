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

use f256::F256;
use f2e16::F2e16;

//pub use all_field_registery;

//pub static mut all_field_registery : Vec<Box< dyn FieldAdd>> = vec![];

#[test]
#[ignore]
fn agreement_f2e16_with_f256() {
     use core::ops::Mul;

     for i in 1..=255 {
         let i_f256 = Additive::<F256>(i).to_multiplier();
         let i_f2e16 = Additive::<F2e16>(i as u16).to_multiplier();
         for j in 0..=255 {
             let j_f256 = Additive::<F256>(j).mul(i_f256);
             let j_f2e16 = Additive::<F2e16>(j as u16).mul(i_f2e16);
             assert_eq!(j_f256.0 as u16, j_f2e16.0);
         }
     }
}
