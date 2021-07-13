
/// Fast Walsh–Hadamard transform over modulo ONEMASK
pub fn walsh<F: FieldAdd>(data: &mut [Logarithm], size: usize) {
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		let depart_no_next = depart_no << 1;
		while j < size {
			for i in j..(depart_no + j) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				let mask = F::ONEMASK as F::Wide;
				let tmp2: F::Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: F::Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				data[i] = Logarithm(((tmp1 & mask) + (tmp1 >> F::FIELD_BITS)) as F::Element);
				data[i + depart_no] = Logarithm(((tmp2 & mask) + (tmp2 >> F::FIELD_BITS)) as F::Element);
			}
			j += depart_no_next;
		}
		depart_no = depart_no_next;
	}
}


#[allow(unused)]
fn bitpoly_mul<F: FieldAdd>(mut a: F::Wide, mut b: F::Wide) -> F::Wide {
    let mut r: F::Wide =0;
    for i in 0..F::FIELD_BITS {
        if (b>>i) & 1 != 0 {
			r ^= a<<i;
		}
    }
    r
}

#[allow(unused)]
fn gf_mul_bitpoly_reduced<F: FieldAdd>(a: F::Element, b: F::Element) -> F::Element {
    use core::convert::TryInto;
    let len = F::FIELD_BITS;
    let mut r: F::Wide = bitpoly_mul(a as F::Wide,b as F::Wide);
    let red : F::Wide = (1 << F::FIELD_BITS) + (F::GENERATOR as F::Wide);
    for i in (len..=(len*2-1)).rev() {
        if r & (1<<i) != 0 {
			r ^= red<<(i-len);
		}
    }
    r.try_into().unwrap()
}

#[test]
fn is_cantor_basis() {
    for w in BASE.windows(2) {
        let b = w[1];
        let square = gf_mul_bitpoly_reduced(b,b);
        let a = w[0];
        // let eq = if a == (square ^ b) { "==" } else { "!=" };
        // println!("{:#b} {} {:#b}\n", a, eq, square ^ b);
        assert_eq!(a, square ^ b);
    }
}

fn generate_cantor_basis<F: FieldAdd>(mut next: F::Element) -> Option<[F::Element; F::FIELD_BITS]> {
    let mut base: [F::Element; F::FIELD_BITS] = [0; F::FIELD_BITS];
    base[0] = 1;
    for b in base.iter_mut().rev() {
        if next == 0 || (next == 1 && *b != 1) { return None; }
        *b = next;
        let square = gf_mul_bitpoly_reduced(next,next);
        next ^= square;
    }
    if next == 0 { Some(base) } else { None }
}

