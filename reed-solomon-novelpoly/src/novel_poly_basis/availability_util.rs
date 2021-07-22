/// Contains interface for collator and validator to generate pieces
/// to distribute in order to ensure availability properties

use core::convert::{TryInto};

use crate::{Shard, ReedSolomon};
use crate::{AfftField};
use crate::errors::*;

use crate::util::recoverablity_subset_size;

/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct<F : AfftField, S: Shard<F>>(
    received_shards: Vec<Option<S>>,
    validator_count: usize,
) -> Result<Vec<u8>>
where <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug,
[u8; F::FIELD_BYTES]: Sized,
[(); F::FIELD_BYTES]: Sized,
[(); F::ONEMASK_USIZE]: Sized,
[(); F::FIELD_SIZE >> 1]: Sized,
{
    let rs = ReedSolomon::<F>::new(validator_count, recoverablity_subset_size(validator_count))?;

    rs.reconstruct(received_shards)
}

pub fn encode<F: AfftField, S: Shard<F>>(bytes: &[u8], validator_count: usize) -> Result<Vec<S>>
where <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug,
[u8; F::FIELD_BYTES]: Sized,
[(); F::FIELD_BYTES]: Sized,
[(); F::ONEMASK_USIZE]: Sized,
[(); F::FIELD_SIZE >> 1]: Sized,
{
    
    let rs = ReedSolomon::<F>::new(validator_count, recoverablity_subset_size(validator_count))?;

    rs.encode::<S>(bytes)
}
