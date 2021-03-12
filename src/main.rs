use color_eyre::Result;

fn main() -> Result<()> {
	color_eyre::install()?;

	#[allow(unused)]
	use reed_solomon_performance::{BYTES, TEST_DATA_CHUNK_SIZE, N_SHARDS};

	#[cfg(feature = "novelpoly")]
	{
		use reed_solomon_performance::novelpoly;
		reed_solomon_tester::roundtrip(novelpoly::encode, novelpoly::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)?;
	}

	#[cfg(feature = "novelpoly-with-alt-cxx-impl")]
	{
		use reed_solomon_performance::novelpoly;
		reed_solomon_tester::roundtrip(novelpoly::cxx::encode, novelpoly::cxx::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)?;
	}

	#[cfg(feature = "naive")]
	{
		use reed_solomon_performance::naive;
		reed_solomon_tester::roundtrip(naive::encode, naive::reconstruct, &BYTES[..TEST_DATA_CHUNK_SIZE], N_SHARDS)?;
	}

	Ok(())
}
