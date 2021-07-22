#[cfg(table_bootstrap_complete)]
use core::convert::{TryFrom, TryInto, Into};
#[cfg(table_bootstrap_complete)]
use super::*;
    
/// Fast Walsh–Hadamard transform over modulo ONEMASK
pub fn walsh<F: FieldAdd>(data: &mut [Logarithm<F>], size: usize) where <<F as FieldAdd>::Wide as TryInto<<F as FieldAdd>::Element>>::Error: core::fmt::Debug{
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		let depart_no_next = depart_no << 1;
		while j < size {
			for i in j..(depart_no + j) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				let mask : F::Wide = F::ONEMASK.into();
				let tmp2: F::Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: F::Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				data[i] = Logarithm(((tmp1 & mask) + (tmp1 >> F::FIELD_BITS)).try_into().ok().unwrap());
				data[i + depart_no] = Logarithm(((tmp2 & mask) + (tmp2 >> F::FIELD_BITS)).try_into().unwrap());
			}
			j += depart_no_next;
		}
		depart_no = depart_no_next;
	}
}


#[allow(unused)]
fn bitpoly_mul<F: FieldAdd>(mut a: F::Wide, mut b: F::Wide) -> F::Wide {
    let mut r: F::Wide = F::ZERO_ELEMENT.into();
    for i in 0..F::FIELD_BITS {
        if (b>>i) & F::ONE_ELEMENT.into() != F::ZERO_ELEMENT.into() {
			r ^= a << i;
		}
    }
    r
}

#[allow(unused)]
pub fn gf_mul_bitpoly_reduced<F: FieldAdd>(a: F::Element, b: F::Element) -> F::Element where <F::Wide as TryInto<<F as FieldAdd>::Element>>::Error : core::fmt::Debug {
    use core::convert::TryInto;
    let len = F::FIELD_BITS;
    let mut r: F::Wide = bitpoly_mul::<F>(a.into(),b.into());
    let red : F::Wide = (F::ONE_ELEMENT_WIDE << F::FIELD_BITS) + F::GENERATOR.into();
    for i in (len..=(len*2-1)).rev() {
        if r & (F::ONE_ELEMENT_WIDE <<i) != F::ZERO_ELEMENT_WIDE {
			r ^= red<<(i-len);
		}
    }
    r.try_into().unwrap()
}

#[cfg(test)]        
fn is_cantor_basis<F: FieldAdd>() {
    for w in F::BASE.windows(2) {
        let b = w[1];
        let square = gf_mul_bitpoly_reduced(b,b);
        let a = w[0];
        // let eq = if a == (square ^ b) { "==" } else { "!=" };
        // println!("{:#b} {} {:#b}\n", a, eq, square ^ b);
        assert_eq!(a, square ^ b);
    }
}

#[cfg(test)]        
test_all_fields_for!(is_cantor_basis);
