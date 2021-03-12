// Encoding/erasure decoding for Reed-Solomon codes over binary extension fields
//
// Derived impl of `RSAErasureCode.c`.
//
// Lin, Han and Chung, "Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes," FOCS14.
// (http://arxiv.org/abs/1404.3458)

#![allow(dead_code)]

use itertools::Itertools;

mod util;
mod wrapped_shard;

use self::wrapped_shard::*;

use crate::errors::*;

pub use self::util::*;

use super::f2e16::*;

//-----Used in decoding procedure-------
//twisted factors used in FFT
static mut SKEW_FACTOR: [Elt; ONEMASK as usize] = [0_u16; ONEMASK as usize];

//factors used in formal derivative
static mut B: [Multiplier; FIELD_SIZE >> 1] = [Multiplier(0_u16); FIELD_SIZE >> 1];

//formal derivative of polynomial in the new basis
pub fn formal_derivative(cos: &mut [Additive], size: usize) {
	for i in 1..size {
		let length = ((i ^ i - 1) + 1) >> 1;
		for j in (i - length)..i {
			cos[j] ^= cos.get(j + length).copied().unwrap_or_default();
		}
	}
	let mut i = size;
	while i < FIELD_SIZE && i < cos.len() {
		for j in 0..size {
			cos[j] ^= cos.get(j + i).copied().unwrap_or_default();
		}
		i <<= 1;
	}
}

// We want the low rate scheme given in
// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf
// and https://github.com/catid/leopard/blob/master/docs/LowRateDecoder.pdf
// but this code resembles https://github.com/catid/leopard which
// implements the high rate decoder in
// https://github.com/catid/leopard/blob/master/docs/HighRateDecoder.pdf
// We're hunting for the differences and trying to undersrtand the algorithm.

//IFFT in the proposed basis
pub fn inverse_afft_in_novel_poly_basis(data: &mut [Additive], size: usize, index: usize) {
	// All line references to Algorithm 2 page 6288 of
	// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf

	// Depth of the recursion on line 7 and 8 is given by depart_no
	// aka 1 << ((k of Algorithm 2) - (i of Algorithm 2)) where
	// k of Algorithm 1 is read as FIELD_BITS here.
	// Recusion base layer implicitly imports d_r aka ala line 1.
	// After this, we start at depth (i of Algorithm 2) = (k of Algorithm 2) - 1
	// and progress through FIELD_BITS-1 steps, obtaining \Psi_\beta(0,0).
	let mut depart_no = 1_usize;
	while depart_no < size {
		// Agrees with for loop (j of Algorithm 2) in (0..2^{k-i-1}) from line 3,
		// except we've j in (depart_no..size).step_by(2*depart_no), meaning
		// the doubled step compensated for the halve size exponent, and
		// somehow this j captures the subscript on \omega_{j 2^{i+1}}.	 (TODO)
		let mut j = depart_no;
		while j < size {
			// At this point loops over i in (j - depart_no)..j give a bredth
			// first loop across the recursion branches from lines 7 and 8,
			// so the i loop corresponds to r in Algorithm 2.  In fact,
			// data[i] and data[i + depart_no] together cover everything,
			// thanks to the outer j loop.

			// Loop on line 3, so i corresponds to j in Algorithm 2
			for i in (j - depart_no)..j {
				// Line 4, justified by (34) page 6288, but
				// adding depart_no acts like the r+2^i superscript.
				data[i + depart_no] ^= data[i];
			}

			// Algorithm 2 indexs the skew factor in line 5 page 6288
			// by i and \omega_{j 2^{i+1}}, but not by r explicitly.
			// We further explore this confusion below. (TODO)
			let skew = Multiplier(unsafe { SKEW_FACTOR[j + index - 1] });
			// It's reasonale to skip the loop if skew is zero, but doing so with
			// all bits set requires justification.	 (TODO)
			if skew.0 != ONEMASK {
				// Again loop on line 3, except skew should depend upon i aka j in Algorithm 2 (TODO)
				for i in (j - depart_no)..j {
					// Line 5, justified by (35) page 6288, but
					// adding depart_no acts like the r+2^i superscript.
					data[i] ^= data[i + depart_no].mul(skew);
				}
			}

			// Increment by double depart_no in agreement with
			// our updating 2*depart_no elements at this depth.
			j += depart_no << 1;
		}
		depart_no <<= 1;
	}
}

//FFT in the proposed basis
pub fn afft_in_novel_poly_basis(data: &mut [Additive], size: usize, index: usize) {
	// All line references to Algorithm 1 page 6287 of
	// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf

	// Depth of the recursion on line 3 and 4 is given by depart_no
	// aka 1 << ((k of Algorithm 1) - (i of Algorithm 1)) where
	// k of Algorithm 1 is read as FIELD_BITS here.
	// Recusion base layer implicitly imports d_r aka ala line 1.
	// After this, we start at depth (i of Algorithm 1) = (k of Algorithm 1) - 1
	// and progress through FIELD_BITS-1 steps, obtaining \Psi_\beta(0,0).
	let mut depart_no = size >> 1_usize;
	while depart_no > 0 {
		// Agrees with for loop (j of Algorithm 1) in (0..2^{k-i-1}) from line 5,
		// except we've j in (depart_no..size).step_by(2*depart_no), meaning
		// the doubled step compensated for the halve size exponent, and
		// somehow this j captures the subscript on \omega_{j 2^{i+1}}.	 (TODO)
		let mut j = depart_no;
		while j < size {
			// At this point loops over i in (j - depart_no)..j give a bredth
			// first loop across the recursion branches from lines 3 and 4,
			// so the i loop corresponds to r in Algorithm 1.  In fact,
			// data[i] and data[i + depart_no] together cover everything,
			// thanks to the outer j loop.

			// Algorithm 1 indexs the skew factor in line 6 aka (28) page 6287
			// by i and \omega_{j 2^{i+1}}, but not by r explicitly.
			// We doubt the lack of explicit dependence upon r justifies
			// extracting the skew factor outside the loop here.
			// As indexing by \omega_{j 2^{i+1}} appears absolute elsewhere,
			// we think r actually appears but the skew factor repeats itself
			// like in (19) in the proof of Lemma 4.  (TODO)
			// We should understand the rest of this basis story, like (8) too.	 (TODO)
			let skew = Multiplier(unsafe { SKEW_FACTOR[j + index - 1] });
			// It's reasonale to skip the loop if skew is zero, but doing so with
			// all bits set requires justification.	 (TODO)
			if skew.0 != ONEMASK {
				// Loop on line 5, except skew should depend upon i aka j in Algorithm 1 (TODO)
				for i in (j - depart_no)..j {
					// Line 6, explained by (28) page 6287, but
					// adding depart_no acts like the r+2^i superscript.
					data[i] ^= data[i + depart_no].mul(skew);
				}
			}

			// Again loop on line 5, so i corresponds to j in Algorithm 1
			for i in (j - depart_no)..j {
				// Line 7, explained by (31) page 6287, but
				// adding depart_no acts like the r+2^i superscript.
				data[i + depart_no] ^= data[i];
			}

			// Increment by double depart_no in agreement with
			// our updating 2*depart_no elements at this depth.
			j += depart_no << 1;
		}
		depart_no >>= 1;
	}
}

//initialize SKEW_FACTOR and B
unsafe fn init_dec() {
	// We cannot yet identify if base has an additive or multiplicative
	// representation, or mybe something else entirely.  (TODO)
	let mut base: [Elt; FIELD_BITS - 1] = Default::default();

	let mut skew_factor: [Additive; ONEMASK as usize] = [Additive(0_u16); ONEMASK as usize];

	for i in 1..FIELD_BITS {
		base[i - 1] = 1 << i;
	}

	// We construct SKEW_FACTOR in additive form to be \bar{s}_j(omega)
	// from page 6285 for all omega in the field.
	for m in 0..(FIELD_BITS - 1) {
		let step = 1 << (m + 1);
		skew_factor[(1 << m) - 1] = Additive(0);
		for i in m..(FIELD_BITS - 1) {
			let s = 1 << (i + 1);

			let mut j = (1 << m) - 1;
			while j < s {
				// Justified by (5) page 6285, except..
				// we expect SKEW_FACTOR[j ^ field_base[i]] or similar
				skew_factor[j + s] = skew_factor[j] ^ Additive(base[i]);
				j += step;
			}
		}

		// Compute base[m] = ONEMASK - base[m] * EXP[LOG[base[m] ^ 1]]
		// = ONEMASK - base[m] * (base[m] ^ 1)
		// TODO: But why?
		//
		// let idx = mul_table(base[m], LOG_TABLE[(base[m] ^ 1_u16) as usize]);
		let idx = Additive(base[m]).mul(Additive(base[m] ^ 1_u16).to_multiplier());
		// WTF?!?
		// base[m] = ONEMASK - LOG_TABLE[idx as usize];
		base[m] = ONEMASK - idx.to_multiplier().0;

		// Compute base[i] = base[i] * EXP[b % ONEMASK]
		// where b = base[m] + LOG[base[i] ^ 1_u16].
		// As ONEMASK is the order of the multiplicative grou,
		// base[i] = base[i] * EXP[base[m]] * (base[m] ^ 1)
		// TODO: But why?
		for i in (m + 1)..(FIELD_BITS - 1) {
			// WTF?!?
			// let b = LOG_TABLE[(base[i] as u16 ^ 1_u16) as usize] as u32 + base[m] as u32;
			let b = Additive(base[i] ^ 1_u16).to_multiplier().to_wide() + (base[m] as Wide);
			let b = b % (ONEMASK as Wide);
			// base[i] = mul_table(base[i], b as u16);
			base[i] = Additive(base[i]).mul(Multiplier(b as Elt)).0;
		}
	}
	// Convert SKEW_FACTOR from Additive to Multiplier form
	for i in 0..(ONEMASK as usize) {
		// SKEW_FACTOR[i] = LOG_TABLE[SKEW_FACTOR[i] as usize];
		SKEW_FACTOR[i] = skew_factor[i].to_multiplier().0;
	}

	#[cfg(test)]
	{
		// TODO: How does this alter base?
		base[0] = ONEMASK - base[0];
		for i in 1..(FIELD_BITS - 1) {
			base[i] = (((ONEMASK as Wide) - (base[i] as Wide) + (base[i - 1] as Wide)) % (ONEMASK as Wide)) as Elt;
		}

		// TODO: What is B anyways?
		B[0] = Multiplier(0);
		for i in 0..(FIELD_BITS - 1) {
			let depart = 1 << i;
			for j in 0..depart {
				B[j + depart] = Multiplier(((B[j].to_wide() + (base[i] as Wide)) % (ONEMASK as Wide)) as Elt);
			}
		}
	} // cfg(test)
}

/// Setup both decoder and encoder.
pub fn setup() {
	use std::sync::Once;

	static SETUP: Once = Once::new();

	SETUP.call_once(|| unsafe {
		init_dec();
	});
}

// Encoding alg for k/n < 0.5: message is a power of two
pub fn encode_low(data: &[Additive], k: usize, codeword: &mut [Additive], n: usize) {
	assert!(k + k <= n);
	assert_eq!(codeword.len(), n);
	assert_eq!(data.len(), n);

	assert!(is_power_of_2(n));
	assert!(is_power_of_2(k));

	// k | n is guaranteed
	assert_eq!((n / k) * k, n);

	// move the data to the codeword
	mem_cpy(&mut codeword[0..], &data[0..]);

	// split after the first k
	let (codeword_first_k, codeword_skip_first_k) = codeword.split_at_mut(k);

	inverse_afft_in_novel_poly_basis(codeword_first_k, k, 0);

	// the first codeword is now the basis for the remaining transforms
	// denoted `M_topdash`

	for shift in (k..n).into_iter().step_by(k) {
		let codeword_at_shift = &mut codeword_skip_first_k[(shift - k)..shift];
		// copy `M_topdash` to the position we are currently at, the n transform
		mem_cpy(codeword_at_shift, codeword_first_k);
		afft_in_novel_poly_basis(codeword_at_shift, k, shift);
	}

	// restore `M` from the derived ones
	mem_cpy(&mut codeword[0..k], &data[0..k]);
}

fn mem_cpy(dest: &mut [Additive], src: &[Additive]) {
	let sl = src.len();
	debug_assert_eq!(dest.len(), sl);
	for i in 0..sl {
		dest[i] = src[i];
	}
}

//data: message array. parity: parity array. mem: buffer(size>= n-k)
//Encoding alg for k/n>0.5: parity is a power of two.
pub fn encode_high(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	let t: usize = n - k;

	// mem_zero(&mut parity[0..t]);
	for i in 0..t {
		parity[i] = Additive(0);
	}

	let mut i = t;
	while i < n {
		mem_cpy(&mut mem[..t], &data[(i - t)..t]);

		inverse_afft_in_novel_poly_basis(mem, t, i);
		for j in 0..t {
			parity[j] ^= mem[j];
		}
		i += t;
	}
	afft_in_novel_poly_basis(parity, t, 0);
}

// Compute the evaluations of the error locator polynomial
// `fn decode_init`
// since this has only to be called once per reconstruction
pub fn eval_error_polynomial(erasure: &[bool], log_walsh2: &mut [Multiplier], n: usize) {
	let z = std::cmp::min(n, erasure.len());
	for i in 0..z {
		log_walsh2[i] = Multiplier(erasure[i] as Elt);
	}
	for i in z..n {
		log_walsh2[i] = Multiplier(0);
	}
	walsh(log_walsh2, FIELD_SIZE);
	for i in 0..n {
		let tmp = log_walsh2[i].to_wide() * LOG_WALSH[i].to_wide();
		log_walsh2[i] = Multiplier((tmp % ONEMASK as Wide) as Elt);
	}
	walsh(log_walsh2, FIELD_SIZE);
	for i in 0..z {
		if erasure[i] {
			log_walsh2[i] = Multiplier(ONEMASK) - log_walsh2[i];
		}
	}
}

/// recover determines how many shards to recover (starting from 0)
// technically we only need to recover
// the first `k` instead of all `n` which
// would include parity chunks.
fn decode_main(codeword: &mut [Additive], recover_up_to: usize, erasure: &[bool], log_walsh2: &[Multiplier], n: usize) {
	assert_eq!(codeword.len(), n);
	assert!(n >= recover_up_to);
	assert_eq!(erasure.len(), n);

	for i in 0..n {
		codeword[i] = if erasure[i] { Additive(0) } else { codeword[i].mul(log_walsh2[i]) };
	}

	inverse_afft_in_novel_poly_basis(codeword, n, 0);

	// formal derivative

	// We should change nothing when multiplying by b from B.
	for i in (0..n).into_iter().step_by(2) {
		let b = Multiplier(ONEMASK) - unsafe { B[i >> 1] };
		#[cfg(test)]
		let x: [_; 2] = [codeword[i], codeword[i + 1]];
		codeword[i] = codeword[i].mul(b);
		codeword[i + 1] = codeword[i + 1].mul(b);
		#[cfg(test)]
		assert_eq!(x, [codeword[i], codeword[i + 1]]);
	}

	formal_derivative(codeword, n);

	// Again changes nothing by multiplying by b although b differs here.
	for i in (0..n).into_iter().step_by(2) {
		#[cfg(test)]
		let x: [_; 2] = [codeword[i], codeword[i + 1]];
		let b = unsafe { B[i >> 1] };
		codeword[i] = codeword[i].mul(b);
		codeword[i + 1] = codeword[i + 1].mul(b);
		#[cfg(test)]
		assert_eq!(x, [codeword[i], codeword[i + 1]]);
	}

	afft_in_novel_poly_basis(codeword, n, 0);

	for i in 0..recover_up_to {
		codeword[i] = if erasure[i] { codeword[i].mul(log_walsh2[i]) } else { Additive(0) };
	}
}

/// Params for the encoder / decoder
/// derived from a target validator count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CodeParams {
	/// total number of message symbols to send
	/// Invariant is a power of base 2
	n: usize,
	/// number of information containing chunks
	/// Invariant is a power of base 2, `k < n`
	k: usize,
	/// Avoid copying unnecessary chunks.
	validator_count: usize,
}

impl CodeParams {
	/// Create a new reed solomon erasure encoding wrapper
	pub fn derive_from_validator_count(validator_count: usize) -> Result<Self> {
		if validator_count < 2 {
			return Err(Error::ValidatorCountTooLow(validator_count));
		}
		// we need to be able to reconstruct from 1/3 - eps
		let k = std::cmp::max(validator_count / 3, 1); // for the odd case of 2 validators
		let k = next_lower_power_of_2(k);
		let n = next_higher_power_of_2(validator_count);
		if n > FIELD_SIZE as usize {
			return Err(Error::ValidatorCountTooLow(validator_count));
		}
		Ok(Self { n, k, validator_count })
	}

	// make a reed-solomon instance.
	pub fn make_encoder(&self) -> ReedSolomon {
		ReedSolomon::new(self.n, self.k, self.validator_count)
			.expect("this struct is not created with invalid shard number; qed")
	}
}

pub fn encode(bytes: &[u8], validator_count: usize) -> Result<Vec<WrappedShard>> {
	let params = CodeParams::derive_from_validator_count(validator_count)?;

	let rs = params.make_encoder();
	rs.encode(bytes)
}

/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct<'a>(received_shards: Vec<Option<WrappedShard>>, validator_count: usize) -> Result<Vec<u8>> {
	let params = CodeParams::derive_from_validator_count(validator_count)?;

	let rs = params.make_encoder();
	rs.reconstruct(received_shards)
}

pub struct ReedSolomon {
	n: usize,
	k: usize,
	validator_count: usize,
}

impl ReedSolomon {
	/// Returns the size per shard in bytes
	pub fn shard_len(&self, payload_size: usize) -> usize {
		let payload_symbols = (payload_size + 1) / 2;
		let shard_symbols_ceil = (payload_symbols + self.k - 1) / self.k;
		let shard_bytes = shard_symbols_ceil * 2;
		shard_bytes
	}

	pub fn new(n: usize, k: usize, validator_count: usize) -> Result<Self> {
		setup();
		if !is_power_of_2(n) && !is_power_of_2(k) {
			Err(Error::ParamterMustBePowerOf2 { n, k })
		} else {
			Ok(Self { validator_count, n, k })
		}
	}

	pub fn encode(&self, bytes: &[u8]) -> Result<Vec<WrappedShard>> {
		if bytes.is_empty() {
			return Err(Error::PayloadSizeIsZero);
		}

		// setup the shards, n is likely _larger_, so use the truely required number of shards

		// required shard length in bytes, rounded to full symbols
		let shard_len = self.shard_len(bytes.len());
		assert!(shard_len > 0);
		// collect all sub encoding runs

		let validator_count = self.validator_count;
		let k2 = self.k * 2;
		// prepare one wrapped shard per validator
		let mut shards = vec![
			WrappedShard::new({
				let mut v = Vec::<u8>::with_capacity(shard_len);
				unsafe { v.set_len(shard_len) }
				v
			});
			validator_count
		];

		for (chunk_idx, i) in (0..bytes.len()).into_iter().step_by(k2).enumerate() {
			let end = std::cmp::min(i + k2, bytes.len());
			assert_ne!(i, end);
			let data_piece = &bytes[i..end];
			assert!(!data_piece.is_empty());
			assert!(data_piece.len() <= k2);
			let encoding_run = encode_sub(data_piece, self.n, self.k)?;
			for val_idx in 0..validator_count {
				AsMut::<[[u8; 2]]>::as_mut(&mut shards[val_idx])[chunk_idx] = encoding_run[val_idx].0.to_be_bytes();
			}
		}

		Ok(shards)
	}

	/// each shard contains one symbol of one run of erasure coding
	pub fn reconstruct(&self, received_shards: Vec<Option<WrappedShard>>) -> Result<Vec<u8>> {
		let gap = self.n.saturating_sub(received_shards.len());

		// obtain a sample of a shard length and assume that is the truth
		// XXX make sure all shards have equal length
		let shard_len_in_syms = received_shards
			.iter()
			.find_map(|x| {
				x.as_ref().map(|x| {
					let x = AsRef::<[[u8; 2]]>::as_ref(x);
					x.len()
				})
			})
			.unwrap();

		let received_shards =
			received_shards.into_iter().take(self.n).chain(std::iter::repeat(None).take(gap)).collect::<Vec<_>>();

		assert_eq!(received_shards.len(), self.n);

		// must be collected after expanding `received_shards` to the anticipated size
		let mut existential_count = 0_usize;
		let erasures = received_shards
			.iter()
			.map(|x| x.is_none())
			.inspect(|erased| existential_count += !*erased as usize)
			.collect::<Vec<bool>>();

		if existential_count < self.k {
			return Err(Error::NeedMoreShards { have: existential_count, min: self.k, all: self.n });
		}

		// Evaluate error locator polynomial only once
		let mut error_poly_in_log = [Multiplier(0); FIELD_SIZE];
		eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

		let mut acc = Vec::<u8>::with_capacity(shard_len_in_syms * 2 * self.k);
		for i in 0..shard_len_in_syms {
			// take the i-th element of all shards and try to recover
			let decoding_run = received_shards
				.iter()
				.map(|x| {
					x.as_ref().map(|x| {
						let z = AsRef::<[[u8; 2]]>::as_ref(&x)[i];
						Additive(u16::from_be_bytes(z))
					})
				})
				.collect::<Vec<Option<Additive>>>();

			assert_eq!(decoding_run.len(), self.n);

			// reconstruct from one set of symbols which was spread over all erasure chunks
			let piece = reconstruct_sub(&decoding_run[..], &erasures, self.n, self.k, &error_poly_in_log).unwrap();
			acc.extend_from_slice(&piece[..]);
		}

		Ok(acc)
	}
}

/// Bytes shall only contain payload data
pub fn encode_sub(bytes: &[u8], n: usize, k: usize) -> Result<Vec<Additive>> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert!(bytes.len() <= k << 1);
	assert!(k <= n / 2);

	// must be power of 2
	let dl = bytes.len();
	let l = if is_power_of_2(dl) {
		dl
	} else {
		let l = log2(dl);
		let l = 1 << l;
		let l = if l >= dl { l } else { l << 1 };
		l
	};
	assert!(is_power_of_2(l));
	assert!(l >= dl);

	// pad the incoming bytes with trailing 0s
	// so we get a buffer of size `N` in `GF` symbols
	let zero_bytes_to_add = n * 2 - dl;
	let data: Vec<Additive> = bytes
		.into_iter()
		.copied()
		.chain(std::iter::repeat(0u8).take(zero_bytes_to_add))
		.tuple_windows()
		.step_by(2)
		.map(|(a, b)| Additive(u16::from_be_bytes([a, b])))
		.collect::<Vec<Additive>>();

	// update new data bytes with zero padded bytes
	// `l` is now `GF(2^16)` symbols
	let l = data.len();
	assert_eq!(l, n);

	let mut codeword = data.clone();
	assert_eq!(codeword.len(), n);

	encode_low(&data[..], k, &mut codeword[..], n);

	Ok(codeword)
}

pub fn reconstruct_sub(
	codewords: &[Option<Additive>],
	erasures: &[bool],
	n: usize,
	k: usize,
	error_poly: &[Multiplier; FIELD_SIZE],
) -> Result<Vec<u8>> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert_eq!(codewords.len(), n);
	assert!(k <= n / 2);

	// the first k suffice for the original k message codewords
	let recover_up_to = k; // n;

	// The recovered _payload_ chunks AND parity chunks
	let mut recovered = vec![Additive(0); recover_up_to];

	// get rid of all `None`s
	let mut codeword = codewords
		.into_iter()
		.enumerate()
		.map(|(idx, sym)| {
			// fill the gaps with `0_u16` codewords
			if let Some(sym) = sym {
				(idx, *sym)
			} else {
				(idx, Additive(0))
			}
		})
		.map(|(idx, codeword)| {
			if idx < recovered.len() {
				recovered[idx] = codeword;
			}
			codeword
		})
		.collect::<Vec<Additive>>();

	// filled up the remaining spots with 0s
	assert_eq!(codeword.len(), n);

	//---------Erasure decoding----------------

	decode_main(&mut codeword[..], recover_up_to, &erasures[..], &error_poly[..], n);

	for idx in 0..recover_up_to {
		if erasures[idx] {
			recovered[idx] = codeword[idx];
		};
	}

	let mut recovered_bytes = Vec::with_capacity(recover_up_to * 2);
	recovered.into_iter().take(k).for_each(|x| recovered_bytes.extend_from_slice(&x.0.to_be_bytes()[..]));
	Ok(recovered_bytes)
}

#[cfg(test)]
mod tests;
