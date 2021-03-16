
use crate::field::{Castomat, AdditiveT, MultiplierT, FieldT};

/// Fast Walsh–Hadamard transform over modulo ONEMASK
pub fn walsh<F: FieldT>(data: &mut [F::Multiplier], size: usize) {
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		let depart_no_next = depart_no << 1;
		while j < size {
			for i in j..(depart_no + j) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				let mask: F::Wide = <F::Wide as Castomat<F::Element>>::from(F::ONEMASK);
				let tmp2: F::Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: F::Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				let field_bits = <F::Wide as Castomat<usize>>::from(F::FIELD_BITS);
				data[i] = F::Multiplier::from((tmp1 & mask) + (tmp1 >> field_bits));
				data[i + depart_no] = F::Multiplier::from(((tmp2 & mask) + (tmp2 >> field_bits)));
			}
			j += depart_no_next;
		}
		depart_no = depart_no_next;
	}
}
