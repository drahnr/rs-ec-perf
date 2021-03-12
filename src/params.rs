
pub const fn log2(mut x: usize) -> usize {
	let mut o: usize = 0;
	while x > 1 {
		x >>= 1;
		o += 1;
	}
	o
}

pub const fn is_nonzero_power_of_2(x: usize) -> bool {
	x > 0_usize && x & (x - 1) == 0
}

use itertools::Itertools;

pub const fn next_higher_power_of_2(k: usize) -> usize {
	if is_nonzero_power_of_2(k) {
		k
	} else {
		1 << (log2(k) + 1)
	}
}

pub const fn next_lower_power_of_2(k: usize) -> usize {
	if is_nonzero_power_of_2(k) {
		k
	} else {
		1 << log2(k)
	}
}

/// Params for the encoder / decoder
/// derived from a target validator count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CodeParams {
	/// total number of message symbols to send
	/// Invariant is a power of base 2
	pub n: usize,
	/// number of information containing chunks
	/// Invariant is a power of base 2, `k < n`
	pub k: usize,
	/// Avoid copying unnecessary chunks.
	pub real_n: usize,
}

impl CodeParams {
	/// Create a new reed solomon erasure encoding wrapper
	pub fn derive_from_third_plus_epsilon(real_n: usize) -> Result<Self> {
		if real_n < 2 {
			return Err(Error::ValidatorCountTooLow(real_n));
		}
		// we need to be able to reconstruct from 1/3 - eps
		let k = std::cmp::max(real_n / 3, 1); // for the odd case of 2 validators
		let k = next_lower_power_of_2(k);
		let n = next_higher_power_of_2(real_n);
		if n > FIELD_SIZE as usize {
            // TODO: Not the correct error here
			return Err(Error::ValidatorCountTooLow(real_n));
		}
		Ok(Self { n, k, real_n })
	}

	pub fn check(&self) -> Result<CheckedParams> {
        let n = self.n;
        let k = self.k;
		if !is_nonzero_power_of_2(n) || !is_nonzero_power_of_2(k) {
			Err(Error::ParamterMustBePowerOf2 { n, k })
		} else {
			Ok(CheckedParams { params: *self })
		}
	}

}


pub struct CheckedParams {
    params: CodeParams,
}
impl core::ops::Deref for CheckedParams {
    type Target = CodeParams;
    fn deref(&self) -> &CodeParams { &self.params }
}

