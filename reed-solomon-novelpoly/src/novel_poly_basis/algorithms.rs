#![allow(dead_code)]

use super::*;

// Encoding alg for k/n < 0.5: message is a power of two
pub(crate) fn encode_low(data: &[Additive], k: usize, codeword: &mut [Additive], n: usize) {
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

//data: message array. parity: parity array. mem: buffer(size>= n-k)
//Encoding alg for k/n>0.5: parity is a power of two.
pub(crate) fn encode_high(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	let t: usize = n - k;

	for i in 0..t {
		parity[i] = Additive(0);
	}

	let mut i = t;
	while i < n {
		(&mut mem[..t]).copy_from_sliec(&data[(i - t)..t]);

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
pub(crate) fn eval_error_polynomial(erasure: &[bool], log_walsh2: &mut [Multiplier], n: usize) {
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
pub(crate) fn decode_main(codeword: &mut [Additive], recover_up_to: usize, erasure: &[bool], log_walsh2: &[Multiplier], n: usize) {
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
