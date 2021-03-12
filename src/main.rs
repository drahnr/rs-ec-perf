use reed_solomon_performance::*;
use reed_solomon_tester::*;

use color_eyre::{Result, eyre, install};

fn main() -> Result<()> {
	color_eyre::install()?;

	#[cfg(feature = "novelpoly")]
	{
		roundtrip(novelpoly::encode, novelpoly::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)?;
	}

	#[cfg(feature = "novelpoly-with-alt-cxx-impl")]
	{
		roundtrip(novelpoly::cxx::encode, novelpoly::cxx::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)?;
	}


	#[cfg(feature = "naive")]
	{
		roundtrip(naive::encode, naive::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)?;
	}

	Ok(())
}
