// Encoding/erasure decoding for Reed-Solomon codes over binary extension fields
//
// Derived impl of `RSAErasureCode.c`.
//
// Lin, Han and Chung, "Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes," FOCS14.
// (http://arxiv.org/abs/1404.3458)

#![allow(dead_code)]

use super::*;

use super::f2e16::*;



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

include!("encode.rs");
include!("decode.rs");


use itertools::Itertools;

pub const fn next_higher_power_of_2(k: usize) -> usize {
	if !is_nonzero_power_of_2(k) {
		1 << (log2(k) + 1)
	} else {
		k
	}
}

pub const fn next_lower_power_of_2(k: usize) -> usize {
	if !is_nonzero_power_of_2(k) {
		1 << log2(k)
	} else {
		k
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
	real_n: usize,
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
			return Err(Error::ValidatorCountTooLow(real_n));
		}
		Ok(Self { n, k, real_n })
	}

	// make a reed-solomon instance.
	pub fn make_encoder(&self) -> ReedSolomon {
		ReedSolomon::new(self.n, self.k, self.real_n)
			.expect("this struct is not created with invalid shard number; qed")
	}
}

pub fn encode(bytes: &[u8], real_n: usize) -> Result<Vec<WrappedShard>> {
	let params = CodeParams::derive_from_third_plus_epsilon(real_n)?;

	let rs = params.make_encoder();
	rs.encode(bytes)
}

/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct(received_shards: Vec<Option<WrappedShard>>, real_n: usize) -> Result<Vec<u8>> {
	let params = CodeParams::derive_from_third_plus_epsilon(real_n)?;

	let rs = params.make_encoder();
	rs.reconstruct(received_shards)
}

pub struct ReedSolomon {
	n: usize,
	k: usize,
	real_n: usize,
}

impl ReedSolomon {
	/// Returns the size per shard in bytes
	pub fn shard_len(&self, payload_size: usize) -> usize {
		let payload_symbols = (payload_size + 1) / 2;
		let shard_symbols_ceil = (payload_symbols + self.k - 1) / self.k;
		let shard_bytes = shard_symbols_ceil * 2;
		shard_bytes
	}

	pub fn new(n: usize, k: usize, real_n: usize) -> Result<Self> {
		// setup();
		if !is_nonzero_power_of_2(n) && !is_nonzero_power_of_2(k) {
			Err(Error::ParamterMustBePowerOf2 { n, k })
		} else {
			Ok(Self { real_n, n, k })
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

		let real_n = self.real_n;
		let k2 = self.k * 2;
		// prepare one wrapped shard per validator
		let mut shards = vec![
			WrappedShard::new({
				let mut v = Vec::<u8>::with_capacity(shard_len);
				unsafe { v.set_len(shard_len) }
				v
			});
			real_n
		];

		for (chunk_idx, i) in (0..bytes.len()).into_iter().step_by(k2).enumerate() {
			let end = std::cmp::min(i + k2, bytes.len());
			assert_ne!(i, end);
			let data_piece = &bytes[i..end];
			assert!(!data_piece.is_empty());
			assert!(data_piece.len() <= k2);
			let encoding_run = encode_sub(data_piece, self.n, self.k)?;
			for val_idx in 0..real_n {
				AsMut::<[[u8; 2]]>::as_mut(&mut shards[val_idx])[chunk_idx] = encoding_run[val_idx].0.to_be_bytes();
			}
		}

		Ok(shards)
	}

	/// each shard contains one symbol of one run of erasure coding
	pub fn reconstruct(&self, received_shards: Vec<Option<WrappedShard>>) -> Result<Vec<u8>> {
		let gap = self.n.saturating_sub(received_shards.len());
		let received_shards =
			received_shards.into_iter().take(self.n).chain(std::iter::repeat(None).take(gap)).collect::<Vec<_>>();

		assert_eq!(received_shards.len(), self.n);

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

		// TODO check shard length is what we'd expect

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
	assert!(is_nonzero_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_nonzero_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert!(bytes.len() <= k << 1);
	assert!(k <= n / 2);

	// must be power of 2
	let dl = bytes.len();
	let l = if is_nonzero_power_of_2(dl) {
		dl
	} else {
		let l = log2(dl);
		let l = 1 << l;
		let l = if l >= dl { l } else { l << 1 };
		l
	};
	assert!(is_nonzero_power_of_2(l));
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
	assert!(is_nonzero_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_nonzero_power_of_2(k), "Algorithm only works for 2^i sizes for K");
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

	// the first k would suffice for the original k message codewords
	let recover_up_to = k;

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
mod test {
	use assert_matches::assert_matches;
	use rand::distributions::Uniform;
	use rand::seq::index::IndexVec;

	use super::*;

	/*
	// If this passes then you do not require the b_not_one feature
	fn b_is_one() {
		let n = FIELD_SIZE >> 1;
		// Everything
		// for i in (0..n) {
		// Just like in decode_main
		for i in (0..n).into_iter().step_by(2) {
			let b = Multiplier(ONEMASK) - unsafe { B[i >> 1] };
			assert_eq!(b, ONEMASK);
		}
	}
	*/

	fn print_sha256(txt: &'static str, data: &[Additive]) {
		use sha2::Digest;
		let data = unsafe { ::std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 2) };

		let mut digest = sha2::Sha256::new();
		digest.update(data);
		println!("sha256(rs|{}):", txt);
		for byte in digest.finalize().into_iter() {
			print!("{:02x}", byte);
		}
		println!("")
	}

	/// Generate a random index
	fn rand_gf_element() -> Additive {
		let mut rng = thread_rng();
		let uni = Uniform::<Elt>::new_inclusive(0, ONEMASK);
		Additive(uni.sample(&mut rng))
	}

	#[test]
	fn base_2_powers_of_2() {
		assert!(!is_nonzero_power_of_2(0));
		for i in 0..20 {
			assert!(is_nonzero_power_of_2(1 << i));
		}
		for i in 0..20 {
			assert!(!is_nonzero_power_of_2(7 << i));
		}
		let mut f = 3;
		for _i in 0..20 {
			f *= 7;
			assert!(!is_nonzero_power_of_2(f));
		}
		assert_eq!(is_nonzero_power_of_2(3), false);
	}

	#[test]
	fn base_2_upper_bound() {
		for i in 1_usize..=1024 {
			let upper = next_higher_power_of_2(i);
			if is_nonzero_power_of_2(i) {
				assert_eq!(upper, i);
			} else {
				assert!(upper > i);
			}
		}
	}

	#[test]
	fn k_n_construction() {
		// skip the two, it's a special case
		for real_n in 3_usize..=8200 {
			let CodeParams { n, k, .. } = CodeParams::derive_from_third_plus_epsilon(real_n).unwrap();
			assert!(real_n <= n, "vc={} <= n={} violated", real_n, n);
			assert!(real_n / 3 >= k, "vc={} / 3 >= k={} violated", real_n, k);
			assert!(real_n >= k * 3, "vc={} <= k={} *3  violated", real_n, k);
		}
	}

	#[test]
	fn flt_back_and_forth() {
		const N: usize = 128;

		let mut data = (0..N).into_iter().map(|_x| rand_gf_element()).collect::<Vec<Additive>>();
		let expected = data.clone();

		afft(&mut data, N, N / 4);

		// make sure something is done
		assert!(data.iter().zip(expected.iter()).filter(|(a, b)| { a != b }).count() > 0);

		inverse_afft(&mut data, N, N / 4);

		itertools::assert_equal(data, expected);
	}

	#[test]
	fn sub_encode_decode() -> Result<()> {
		// setup();
		let mut rng = rand::thread_rng();

		const N: usize = 32;
		const K: usize = 4;

		const K2: usize = K * 2;
		let mut data = [0u8; K2];
		rng.fill_bytes(&mut data[..]);

		let codewords = encode_sub(&data, N, K)?;
		let mut codewords = codewords.into_iter().map(|x| Some(x)).collect::<Vec<_>>();
		assert_eq!(codewords.len(), N);
		codewords[0] = None;
		codewords[1] = None;
		codewords[2] = None;
		codewords[N - 3] = None;
		codewords[N - 2] = None;
		codewords[N - 1] = None;

		let erasures = codewords.iter().map(|x| x.is_none()).collect::<Vec<bool>>();

		// Evaluate error locator polynomial only once
		let mut error_poly_in_log = [Multiplier(0); FIELD_SIZE];
		eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

		let reconstructed = reconstruct_sub(&codewords[..], &erasures[..], N, K, &error_poly_in_log)?;
		itertools::assert_equal(data.iter(), reconstructed.iter().take(K2));
		Ok(())
	}

	fn deterministic_drop_shards<T: Sized, G: rand::SeedableRng + rand::Rng>(
		codewords: &mut [Option<T>],
		n: usize,
		k: usize,
		_rng: &mut G,
	) -> IndexVec {
		let l = codewords.len();
		let mut v = Vec::with_capacity(n - k);
		// k is a power of 2
		let half = (n - k) >> 1;
		for i in 0..half {
			codewords[i] = None;
			v.push(i);
		}
		// if the codewords is shorter than n
		// the remaining ones were
		// already dropped implicitly
		for i in n - half..n {
			if i < l {
				codewords[i] = None;
				v.push(i);
			}
		}
		IndexVec::from(v)
	}

	fn deterministic_drop_shards_clone<T: Sized + Clone>(
		codewords: &[T],
		n: usize,
		k: usize,
	) -> (Vec<Option<T>>, IndexVec) {
		let mut rng = SmallRng::from_seed(crate::SMALL_RNG_SEED);
		let mut codewords = codewords.into_iter().map(|x| Some(x.clone())).collect::<Vec<Option<T>>>();
		let idx = deterministic_drop_shards::<T, SmallRng>(&mut codewords, n, k, &mut rng);
		assert!(idx.len() <= n - k);
		(codewords, idx)
	}

	// for shards of length 1
	fn wrapped_shard_len1_as_gf_sym(w: &WrappedShard) -> Additive {
		let val = AsRef::<[[u8; 2]]>::as_ref(w)[0];
		Additive(u16::from_be_bytes(val))
	}

	#[test]
	fn sub_eq_big_for_small_messages() {
		const N_VALIDATORS: usize = 128;
		const N: usize = N_VALIDATORS;
		const K: usize = 32;

		// setup();

		const K2: usize = K * 2;

		// assure the derived sizes match
		let rs = CodeParams::derive_from_third_plus_epsilon(N_VALIDATORS).unwrap();
		assert_eq!(rs.n, N);
		assert_eq!(rs.k, K);

		// create random predictable bytes
		// and create a message that results in 1 GF element symbols
		// per validator
		let data = {
			let mut rng = SmallRng::from_seed(crate::SMALL_RNG_SEED);
			let mut data = [0u8; K2];
			rng.fill_bytes(&mut data[..]);
			data
		};

		let mut codewords = encode(&data, rs.n).unwrap();
		let mut codewords_sub = encode_sub(&data, N, K).unwrap();

		itertools::assert_equal(codewords.iter().map(wrapped_shard_len1_as_gf_sym), codewords_sub.iter().copied());

		let (codewords, _) = deterministic_drop_shards_clone(&mut codewords, N, K);
		let (codewords_sub, _) = deterministic_drop_shards_clone(&mut codewords_sub, N, K);

		itertools::assert_equal(
			codewords.iter().map(|w| w.as_ref().map(wrapped_shard_len1_as_gf_sym)),
			codewords_sub.iter().copied(),
		);

		let erasures = codewords.iter().map(|x| x.is_none()).collect::<Vec<bool>>();

		// Evaluate error locator polynomial only once
		let mut error_poly_in_log = [Multiplier(0); FIELD_SIZE];
		eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

		let reconstructed_sub = reconstruct_sub(&codewords_sub[..], &erasures[..], N, K, &error_poly_in_log).unwrap();
		let reconstructed = reconstruct(codewords, rs.n).unwrap();
		itertools::assert_equal(reconstructed.iter().take(K2), reconstructed_sub.iter().take(K2));
		itertools::assert_equal(reconstructed.iter().take(K2), data.iter());
		itertools::assert_equal(reconstructed_sub.iter().take(K2), data.iter());
	}

	#[test]
	fn roundtrip_for_large_messages() -> Result<()> {
		const N_VALIDATORS: usize = 2000;
		const N: usize = 2048;
		const K: usize = 512;

		// setup();

		const K2: usize = K * 2;

		// assure the derived sizes match
		let rs = CodeParams::derive_from_third_plus_epsilon(N_VALIDATORS).unwrap();
		assert_eq!(rs.n, N);
		assert_eq!(rs.k, K);

		// make sure each shard is more than one byte to
		// test the shard size
		// in GF symbols
		let shard_length: usize = 23;

		let payload = &crate::BYTES[0..K2 * shard_length];
		// let payload = &crate::BYTES[..];

		let mut shards = encode(payload, N_VALIDATORS).unwrap();

		// for (idx, shard) in shards.iter().enumerate() {
		//	let sl = AsRef::<[[u8; 2]]>::as_ref(&shard).len();
		//	assert_eq!(shard_length, sl, "Shard #{} has an unxpected length {} (expected: {})", idx, sl, shard_length);
		// }
		let (received_shards, dropped_indices) = deterministic_drop_shards_clone(&mut shards, rs.n, rs.k);

		let reconstructed_payload = reconstruct(received_shards, N_VALIDATORS).unwrap();

		assert_recovery(payload, &reconstructed_payload, dropped_indices);

		// verify integrity with criterion tests
		roundtrip_w_drop_closure::<_, _, _, SmallRng>(
			encode,
			reconstruct,
			payload,
			N_VALIDATORS,
			deterministic_drop_shards::<WrappedShard, SmallRng>,
		)?;

		roundtrip_w_drop_closure::<_, _, _, SmallRng>(
			encode,
			reconstruct,
			payload,
			N_VALIDATORS,
			crate::drop_random_max,
		)?;

		Ok(())
	}

	macro_rules! simplicissimus {
		($name:ident: validators: $real_n:literal, payload: $payload_size:literal; $matchmaker:pat) => {
			simplicissimus!($name: validators: $real_n, payload: $payload_size; $matchmaker => {});
		};
		($name:ident: validators: $real_n:literal, payload: $payload_size:literal) => {
			simplicissimus!($name: validators: $real_n, payload: $payload_size; Ok(x) => { let _ = x; });
		};
		($name:ident: validators: $real_n:literal, payload: $payload_size:literal; $matchmaker:pat => $assertive:expr) => {
			#[test]
			fn $name () {
				let res = roundtrip_w_drop_closure::<_,_,_,SmallRng>(
					encode,
					reconstruct,
					&BYTES[0..$payload_size], $real_n,
					 deterministic_drop_shards::<WrappedShard, SmallRng>);
				assert_matches::assert_matches!(res, $matchmaker => {
					$assertive
				});
			}
		};
	}

	simplicissimus!(case_0: validators: 2003, payload: 0; Err(Error::PayloadSizeIsZero));

	// Roughly one GFSymbol per validator payload
	simplicissimus!(case_1: validators: 10, payload: 16);

	// Unit payload, but mayn validators
	simplicissimus!(case_2: validators: 100, payload: 1);

	// Common case of way ore payload than validators
	simplicissimus!(case_3: validators: 4, payload: 100);

	// Way more validators than payload bytes
	simplicissimus!(case_4: validators: 2003, payload: 17);

	#[test]
	fn flt_roundtrip_small() {
		const N: usize = 16;
		const EXPECTED: [Additive; N] =
			unsafe { std::mem::transmute([1_u16, 2, 3, 5, 8, 13, 21, 44, 65, 0, 0xFFFF, 2, 3, 5, 7, 11]) };

		let mut data = EXPECTED.clone();

		afft(&mut data, N, N / 4);

		println!("novel basis(rust):");
		data.iter().for_each(|sym| {
			print!(" {:04X}", sym.0);
		});
		println!("");

		inverse_afft(&mut data, N, N / 4);
		itertools::assert_equal(data.iter(), EXPECTED.iter());
	}

	#[test]
	fn ported_c_test() {
		const N: usize = 256;
		const K: usize = 8;

		// setup();

		//-----------Generating message----------
		//message array
		let mut data = [Additive(0); N];

		for i in 0..K {
			//filled with random numbers
			data[i] = Additive((i * i % ONEMASK as usize) as u16);
			// data[i] = rand_gf_element();
		}

		assert_eq!(data.len(), N);

		println!("Message(Last n-k are zeros): ");
		for i in 0..K {
			print!("{:04x} ", data[i].0);
		}
		println!("");
		print_sha256("data", &data[..]);

		//---------encoding----------
		let mut codeword = [Additive(0); N];

		if K + K > N && false {
			let (data_till_t, data_skip_t) = data.split_at_mut(N - K);
			encode_high(data_skip_t, K, data_till_t, &mut codeword[..], N);
		} else {
			encode_low(&data[..], K, &mut codeword[..], N);
		}

		// println!("Codeword:");
		// for i in K..(K+100) {
		// print!("{:04x} ", codeword[i]);
		// }
		// println!("");

		print_sha256("encoded", &codeword);

		//--------erasure simulation---------

		// Array indicating erasures
		let mut erasure = [false; N];

		let erasures_iv = if false {
			// erase random `(N-K)` codewords
			let mut rng = rand::thread_rng();
			let erasures_iv: IndexVec = rand::seq::index::sample(&mut rng, N, N - K);

			erasures_iv
		} else {
			IndexVec::from((0..(N - K)).into_iter().collect::<Vec<usize>>())
		};
		assert_eq!(erasures_iv.len(), N - K);

		for i in erasures_iv {
			//erasure codeword symbols
			erasure[i] = true;
			codeword[i] = Additive(0);
		}

		print_sha256("erased", &codeword);

		//---------Erasure decoding----------------
		let mut log_walsh2: [Multiplier; FIELD_SIZE] = [Multiplier(0); FIELD_SIZE];

		eval_error_polynomial(&erasure[..], &mut log_walsh2[..], FIELD_SIZE);

		// TODO: Make print_sha256 polymorphic
		// print_sha256("log_walsh2", &log_walsh2);

		decode_main(&mut codeword[..], K, &erasure[..], &log_walsh2[..], N);

		print_sha256("decoded", &codeword[0..K]);

		println!("Decoded result:");
		for i in 0..N {
			// the data word plus a few more
			print!("{:04x} ", codeword[i].0);
		}
		println!("");

		for i in 0..K {
			//Check the correctness of the result
			if data[i] != codeword[i] {
				println!("🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍");
				panic!("Decoding ERROR! value at [{}] should={:04x} vs is={:04x}", i, data[i].0, codeword[i].0);
			}
		}
		println!(
			r#">>>>>>>>> 🎉🎉🎉🎉
>>>>>>>>> > Decoding is **SUCCESS** ful! 🎈
>>>>>>>>>"#
		);
	}

	#[test]
	fn test_code_params() {
		assert_matches!(CodeParams::derive_from_third_plus_epsilon(0), Err(_));

		assert_matches!(CodeParams::derive_from_third_plus_epsilon(1), Err(_));

		assert_eq!(CodeParams::derive_from_third_plus_epsilon(2), Ok(CodeParams { n: 2, k: 1, real_n: 2 }));

		assert_eq!(CodeParams::derive_from_third_plus_epsilon(3), Ok(CodeParams { n: 4, k: 1, real_n: 3 }));

		assert_eq!(CodeParams::derive_from_third_plus_epsilon(4), Ok(CodeParams { n: 4, k: 1, real_n: 4 }));

		assert_eq!(
			CodeParams::derive_from_third_plus_epsilon(100),
			Ok(CodeParams { n: 128, k: 32, real_n: 100 })
		);
	}

	#[test]
	fn shard_len_is_reasonable() {
		let rs = CodeParams { n: 16, k: 4, real_n: 5 }.make_encoder();

		// since n must be a power of 2
		// the chunk sizes becomes slightly larger
		// than strictly necessary
		assert_eq!(rs.shard_len(100), 26);
		assert_eq!(rs.shard_len(99), 26);

		// see if it rounds up to 2.
		assert_eq!(rs.shard_len(95), 24);
		assert_eq!(rs.shard_len(94), 24);

		assert_eq!(rs.shard_len(90), 24);

		// needs 3 bytes to fit, rounded up to next even number.
		assert_eq!(rs.shard_len(19), 6);
	}
}
