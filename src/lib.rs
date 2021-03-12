use rand::prelude::*;
use rand::seq::index::IndexVec;


use reed_solomon_tester::BYTES;

#[cfg(feature = "novelpoly")]
pub use reed_solomon_novelpoly as novelpoly;

#[cfg(feature = "naive")]
pub mod naive;

use color_eyre::eyre::Result;

const N_SHARDS: usize = 123;
const TEST_DATA_CHUNK_SIZE: usize = 1337;

#[cfg(test)]
mod test {
	use super::*;

	#[cfg(feature = "naive")]
	#[test]
	fn naive_roundtrip() -> Result<()> {
		roundtrip(naive::encode, naive::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)
	}

	#[cfg(feature = "novelpoly")]

	#[test]
	fn novelpoly_roundtrip() -> Result<()> {
		roundtrip(novelpoly::encode, novelpoly::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)
	}
}
