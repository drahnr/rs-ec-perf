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

	if t % 8 == 0 && n % 8 == 0 && k % 8 == 0 {
		let mut mem8 = vec![Additive8x::zero(); mem.len()/8];
		let data8 = Vec::from_iter(data.iter().step_by(8).enumerate().map(|(piece_idx, _offset)| {
			Additive8x::load(&data[(piece_idx*8)..(piece_idx+1)*8])
		}));
		let mut parity8 = Vec::from_iter(data.iter().step_by(8).enumerate().map(|(piece_idx, _offset)| {
			Additive8x::load(&*parity[(piece_idx*8)..(piece_idx+1)*8])
		}));
		
		encode_high_faster8(&data8, &mut parity8, &mut mem8, n);
		
		parity8.copy_to(&mut parity);
		mem8.copy_to(&mut mem);
		return;
	}
	
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Additive8x(pub simd::u16x8);

pub fn encode_high_faster8(data: &[Additive8x], k: usize, parity: &mut [Additive8x], mem: &mut [Additive8x], n: usize) {
	let t: usize = n - k;
	assert!(t >= 8);
	assert_eq!(t % 8, 0);

	let t8s = (t >> 3);
	for i in 0..t8s {
		parity[i] = Additive8x::zero();
	}

	let mut i = t8s;
	while i < n {
		(&mut mem[..t8s]).copy_from_slice(&data[(i - t8s)..t]);

		inverse_afft_faster8(mem, t8s, i);
		for j in 0..t8s {
			parity[j] ^= mem[j];
		}
		i += t8s;
	}
	afft_faster8(parity, t8s, 0);
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

	encode_low(&data[..], k, &mut codeword[..], n);

	Ok(codeword)
}
