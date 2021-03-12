pub const fn log2(mut x: usize) -> usize {
	let mut o: usize = 0;
	while x > 1 {
		x >>= 1;
		o += 1;
	}
	o
}

pub const fn is_power_of_2(x: usize) -> bool {
	x > 0_usize && x & (x - 1) == 0
}

pub const fn next_higher_power_of_2(k: usize) -> usize {
	if !is_power_of_2(k) {
		1 << (log2(k) + 1)
	} else {
		k
	}
}

pub const fn next_lower_power_of_2(k: usize) -> usize {
	if !is_power_of_2(k) {
		1 << log2(k)
	} else {
		k
	}
}
