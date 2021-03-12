pub mod errors;
pub use errors::*;

pub mod f2e16;

mod novel_poly_basis;
pub use novel_poly_basis::*;

#[cfg(feature = "with-alt-cxx-impl")]
pub mod cxx;

#[cfg(test)]
mod test {
	use super::*;
	use reed_solomon_tester::{roundtrip, BYTES};

	const N_SHARDS: usize = 2000;

	#[cfg(feature = "naive")]
	#[test]
	fn status_quo_roundtrip() -> Result<()> {
		roundtrip(status_quo::encode, status_quo::reconstruct, &BYTES[..1337], N_SHARDS)
	}

	#[test]
	fn novel_poly_basis_roundtrip() -> Result<()> {
		roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &BYTES[..1337], N_SHARDS)
	}
}
