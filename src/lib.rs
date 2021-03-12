#[cfg(feature = "novelpoly")]
pub use reed_solomon_novelpoly as novelpoly;

#[cfg(feature = "naive")]
pub mod naive;

pub const N_SHARDS: usize = 123;
pub const TEST_DATA_CHUNK_SIZE: usize = 1337;
pub use reed_solomon_tester::BYTES;

#[cfg(test)]
mod test {

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
