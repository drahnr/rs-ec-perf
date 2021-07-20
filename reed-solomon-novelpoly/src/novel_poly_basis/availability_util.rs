/// Contains interface for collator and validator to generate pieces
/// to distribute in order to ensure availability properties

use crate::{Shard, ReedSolomon};
use crate::field::afft::*;
use crate::field::FieldAdd;
use crate::field::TruncateTo;

//use crate::shard::ShardHold;

pub use super::util::*;

use super::field::f2e16;

/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct<S: Shard<F: AfftField>>(
    received_shards: Vec<Option<S>>,
    validator_count: usize,
) -> Result<Vec<u8>> {
    let rs = ReedSolomon::<F>::new(validator_count, recoverablity_subset_size(validator_count))?;

    rs.reconstruct(received_shards)
}

pub fn encode<S: Shard<F: AfftField>>(bytes: &[u8], validator_count: usize) -> Result<Vec<S>> {
    let rs = ReedSolomon::<F>::new(validator_count, recoverablity_subset_size(validator_count))?;

    rs.encode::<S>(bytes)
}
