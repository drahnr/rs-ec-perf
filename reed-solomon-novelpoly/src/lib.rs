#![forbid(unused_crate_dependencies)]
#![feature(const_generics)]
#![feature(const_evaluatable_checked)]
#![feature(destructuring_assignment)]
#![feature(associated_type_defaults)]
#![feature(array_windows)]

pub mod errors;
pub use errors::*;

pub mod util;
pub use util::*;

pub mod field;
pub use self::field::f256;
pub use self::field::f2e16;
pub use self::field::f256::F256;
pub use self::field::f2e16::F2e16;
pub use self::field::{FieldAdd, FieldMul, TruncateTo, Logarithm, Additive, AfftField, walsh, gf_mul_bitpoly_reduced};

#[macro_use]
pub use self::field::macros;

mod novel_poly_basis;
pub use self::novel_poly_basis::*;
pub use self::novel_poly_basis::availability_util::*;

pub mod shard;
pub use self::shard::{Shard};

pub mod wrapped_shard;
pub use self::wrapped_shard::WrappedShard;

#[cfg(feature = "with-alt-cxx-impl")]
pub mod cxx;

#[cfg(test)]
mod test {
    use super::*;
    use std::convert::TryInto;
    use reed_solomon_tester::{roundtrip, BYTES, N_SHARDS};

    #[cfg(feature = "naive")]
    fn status_quo_roundtrip<F: FieldAdd>() -> Result<()> {
        roundtrip(status_quo::encode::<F, WrappedShard>, status_quo::reconstruct::<F, WrappedShard>, &BYTES[..1337], N_SHARDS)
    }

    fn novel_poly_basis_roundtrip<F: AfftField>() -> Result<()>
        where
     [u8; F::FIELD_BYTES]: Sized,
    [(); F::FIELD_BYTES]: Sized,
    [(); F::ONEMASK_USIZE]: Sized,
    [(); F::FIELD_SIZE >> 1]: Sized,
    <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug
    {
        roundtrip(
            novel_poly_basis::encode::<F, WrappedShard>,
            novel_poly_basis::reconstruct::<F, WrappedShard>,
            &BYTES[..1337],
            N_SHARDS,
        )
    }
    test_all_fields_for!(novel_poly_basis_roundtrip);
    
}
