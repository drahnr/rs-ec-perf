// Encoding/erasure decoding for Reed-Solomon codes over binary extension fields
//
// Derived impl of `RSAErasureCode.c`.
//
// Lin, Han and Chung, "Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes," FOCS14.
// (http://arxiv.org/abs/1404.3458)

use std::marker::PhantomData;

use crate::errors::*;
use crate::f2e16::*;
use crate::Shard;
use crate::field::afft::*;

pub use super::util::*;

use super::field::f2e16;

/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct<'a, S: Shard>(received_shards: Vec<Option<S>>, validator_count: usize) -> Result<Vec<u8>> {
	let rs : ReedSolomon<f2e16::Additive> = ErasureCodeFactory::make_encoder(validator_count, recoverablity_subset_size(validator_count))?;

	rs.reconstruct(received_shards)
}

pub fn encode<S: Shard>(bytes: &[u8], validator_count: usize) -> Result<Vec<S>> {
	let rs : ReedSolomon<f2e16::Additive> = ErasureCodeFactory::make_encoder(validator_count, recoverablity_subset_size(validator_count))?;

	rs.encode::<S>(bytes)
}

/// Params for the encoder / decoder
/// derived from a target validator count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ErasureCodeFactory<F: AfftField> {
    _marker: PhantomData<F>,
}

impl <F: AfftField> ErasureCodeFactory<F> {
	/// Create a new reed solomon erasure encoding wrapper
	/// `k` the intended number of data shards needed to recover.
	/// `n` the intended number of resulting shards.
	///
	/// Assures that the derived paramters retain at most the given coding
	/// rate, and as such assure recoverability with at least an equiv fraction
	/// as provided by the input `n`, and `k` parameterset.
	pub fn make_encoder(n: usize, k: usize) -> Result<ReedSolomon<F>> {
		if n < 2 {
			return Err(Error::WantedShardCountTooLow(n));
		}
		if k < 1 {
			return Err(Error::WantedPayloadShardCountTooLow(k));
		}
		let k_po2 = next_lower_power_of_2(k);
		let n_po2 = next_higher_power_of_2(n);
		// If the coding rate of the power of 2 variants, is higher,
		// we would have to lower k by one order of magnitude base 2
		// which is true by definition
		assert!(n * k_po2 <= n_po2 * k);

		if n_po2 > FIELD_SIZE as usize {
			return Err(Error::WantedShardCountTooHigh(n));
		}

	    // make a reed-solomon instance.
		Ok(ReedSolomon::new(n_po2, k_po2, n)
			.expect("this struct is not created with invalid shard number; qed"))		
	}
}

pub struct ReedSolomon<F: AfftField> {
	/// Avoid copying unnecessary chunks.
	wanted_n: usize,
	/// total number of message symbols to send
	/// Invariant is a power of base 2
	n: usize,
	/// number of information containing chunks
	/// Invariant is a power of base 2, `k < n`
	k: usize,
    _marker: PhantomData<F>,
}

impl <F:AfftField> ReedSolomon<F> {
	/// Returns the size per shard in bytes
	pub fn shard_len(&self, payload_size: usize) -> usize {
		let payload_symbols = (payload_size + 1) / 2;
		let shard_symbols_ceil = (payload_symbols + self.k - 1) / self.k;
		let shard_bytes = shard_symbols_ceil * 2;
		shard_bytes
	}

	pub(crate) fn new(n: usize, k: usize, wanted_n: usize) -> Result<Self> {
		if !is_power_of_2(n) && !is_power_of_2(k) {
			Err(Error::ParamterMustBePowerOf2 { n, k })
		} else {
			Ok(Self { wanted_n, n, k, _marker: PhantomData })
		}
	}

	pub fn encode<S: Shard>(&self, bytes: &[u8]) -> Result<Vec<S>> {
		if bytes.is_empty() {
			return Err(Error::PayloadSizeIsZero);
		}

		// setup the shards, n is likely _larger_, so use the truely required number of shards

		// required shard length in bytes, rounded to full symbols
		let shard_len = self.shard_len(bytes.len());
		assert!(shard_len > 0);
		// collect all sub encoding runs

		let validator_count = self.wanted_n;
		let k2 = self.k * 2;
		// prepare one wrapped shard per validator
		let mut shards = vec![
			<S as From<Vec<u8>>>::from({
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
			let encoding_run = Self::encode_sub(data_piece, self.n, self.k)?;
			for val_idx in 0..validator_count {
				AsMut::<[[u8; 2]]>::as_mut(&mut shards[val_idx])[chunk_idx] = encoding_run[val_idx].0.to_be_bytes();
			}
		}

		Ok(shards)
	}

	/// each shard contains one symbol of one run of erasure coding
	pub fn reconstruct<S: Shard>(&self, received_shards: Vec<Option<S>>) -> Result<Vec<u8>> {
		let gap = self.n.saturating_sub(received_shards.len());

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

		// obtain a sample of a shard length and assume that is the truth
		let shard_len_in_syms = {
			let (first_shard_idx, first_shard_len) = received_shards
				.iter()
				.enumerate()
				.find_map(|(idx, shard)| {
					shard.as_ref().map(|shard| {
						let shard = AsRef::<[[u8; 2]]>::as_ref(shard);
						(idx, shard.len())
					})
				})
				.expect("Existential shard count is at least k shards. qed");

			// make sure all shards have the same length as the first one
			if let Some(other_shard_len) = received_shards[(first_shard_idx + 1)..].iter().find_map(|shard| {
				shard.as_ref().and_then(|shard| {
					let shard = AsRef::<[[u8; 2]]>::as_ref(shard);
					if first_shard_len != shard.len() {
						Some(shard.len())
					} else {
						None
					}
				})
			}) {
				return Err(Error::InconsistentShardLengths { first: first_shard_len, other: other_shard_len });
			}

			first_shard_len
		};

		// Evaluate error locator polynomial only once
		let mut error_poly_in_log = [Logarithm(0); FIELD_SIZE];
		Self::eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

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
			let piece =
				Self::reconstruct_sub(&decoding_run[..], &erasures, self.n, self.k, &error_poly_in_log).unwrap();
			acc.extend_from_slice(&piece[..]);
		}

		Ok(acc)
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
        codeword.copy_from_slice(data);

	    // split after the first k
	    let (codeword_first_k, codeword_skip_first_k) = codeword.split_at_mut(k);

        inverse_afft(codeword_first_k, k, 0);

	    // the first codeword is now the basis for the remaining transforms
	    // denoted `M_topdash`

	    for shift in (k..n).into_iter().step_by(k) {
		    let codeword_at_shift = &mut codeword_skip_first_k[(shift - k)..shift];
		    // copy `M_topdash` to the position we are currently at, the n transform
		    codeword_at_shift.copy_from_slice(codeword_first_k);
		    afft(codeword_at_shift, k, shift);
	    }

	    // restore `M` from the derived ones
	    (&mut codeword[0..k]).copy_from_slice(&data[0..k]);
    }


    // TODO: Make encode_high actually work again!  Add tests!

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
		    (&mut mem[..t]).copy_from_slice(&data[(i - t)..t]);

		    inverse_afft(mem, t, i);
		    for j in 0..t {
			    parity[j] ^= mem[j];
		    }
		    i += t;
	    }
	    afft(parity, t, 0);
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

        // tuple_windows are only used here
        use itertools::Itertools;

	    // pad the incoming bytes with trailing 0s
	    // so we get a buffer of size `N` in `GF` symbols
	    let zero_bytes_to_add = n * 2 - dl;
	    let data: Vec<Additive> = bytes
		    .into_iter()
		    .copied()
		    .chain(std::iter::repeat(0u8).take(zero_bytes_to_add))
		    .tuple_windows()
		    .step_by(2)
		    .map(|(a, b)| Additive(Elt::from_be_bytes([a, b])))
		    .collect::<Vec<Additive>>();

	    // update new data bytes with zero padded bytes
	    // `l` is now `GF(2^16)` symbols
	    let l = data.len();
	    assert_eq!(l, n);

	    let mut codeword = data.clone();
	    assert_eq!(codeword.len(), n);

	    Self::encode_low(&data[..], k, &mut codeword[..], n);

	    Ok(codeword)
    }

    
    pub fn reconstruct_sub(
	    codewords: &[Option<Additive>],
	    erasures: &[bool],
	    n: usize,
	    k: usize,
	    error_poly: &[Logarithm; FIELD_SIZE],
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

	    Self::decode_main(&mut codeword[..], recover_up_to, &erasures[..], &error_poly[..], n);

	    for idx in 0..recover_up_to {
		if erasures[idx] {
			recovered[idx] = codeword[idx];
		};
	    }

	    let mut recovered_bytes = Vec::with_capacity(recover_up_to * 2);
	    recovered.into_iter().take(k).for_each(|x| recovered_bytes.extend_from_slice(&x.0.to_be_bytes()[..]));
	    Ok(recovered_bytes)
    }
    
    
    
    /// recover determines how many shards to recover (starting from 0)
    // technically we only need to recover
    // the first `k` instead of all `n` which
    // would include parity chunks.
    pub(crate) fn decode_main(codeword: &mut [Additive], recover_up_to: usize, erasure: &[bool], log_walsh2: &[Logarithm], n: usize) {
	    assert_eq!(codeword.len(), n);
	assert!(n >= recover_up_to);
	    assert_eq!(erasure.len(), n);

	    for i in 0..n {
		codeword[i] = if erasure[i] { Additive(0) } else { codeword[i].mul(log_walsh2[i]) };
	    }

	    inverse_afft(codeword, n, 0);

	    tweaked_formal_derivative(codeword, n);

	    afft(codeword, n, 0);

	    for i in 0..recover_up_to {
		codeword[i] = if erasure[i] { codeword[i].mul(log_walsh2[i]) } else { Additive(0) };
	    }
}
    

    // Compute the evaluations of the error locator polynomial
    // `fn decode_init`
    // since this has only to be called once per reconstruction
    pub fn eval_error_polynomial(erasure: &[bool], log_walsh2: &mut [Logarithm], n: usize) {
	let z = std::cmp::min(n, erasure.len());
	    for i in 0..z {
		log_walsh2[i] = Logarithm(erasure[i] as Elt);
	    }
	for i in z..n {
		log_walsh2[i] = Logarithm(0);
	}
	    walsh(log_walsh2, FIELD_SIZE);
	for i in 0..n {
		let tmp = log_walsh2[i].to_wide() * LOG_WALSH[i].to_wide();
		log_walsh2[i] = Logarithm((tmp % ONEMASK as Wide) as Elt);
	}
	walsh(log_walsh2, FIELD_SIZE);
	    for i in 0..z {
		if erasure[i] {
			log_walsh2[i] = Logarithm(ONEMASK) - log_walsh2[i];
		}
	    }
    }
    
    
}

#[cfg(test)]
mod tests;
