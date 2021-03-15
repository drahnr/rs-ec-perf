
use super::{FieldT, Castomat};

use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign};

/// Additive via XOR form of f2e16
#[derive(Clone, Copy, Debug, Default, BitXor, BitXorAssign, PartialEq, Eq)] // PartialOrd,Ord
pub struct Additive<F: FieldT>(pub F::Element);

impl<F: FieldT> Additive<F> {
    #[inline(always)]
	pub fn to_wide(self) -> F::Wide {
		F::Wide::from(self.0)
	}

    #[inline(always)]
	pub fn from_wide(x: F::Wide) -> Additive<F> {
		Additive(F::Element::from(x))
	}

	const ZERO: Additive<F> = Additive(0.cast_as() as F::Element);
}

#[cfg(table_bootstrap_complete)]
impl<F: FieldT> Additive<F> {
	/// Return multiplier prepared form
    #[inline(always)]
	pub fn to_multiplier(self) -> Multiplier<F> {
		Multiplier(LOG_TABLE[self.0.cast_as()])
	}

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
	pub fn mul(self, other: Multiplier<F>) -> Additive<F> {
		if self == Self::ZERO {
			return Self::ZERO;
		}
		let log = F::Wide::cast_from(LOG_TABLE[self.0.cast_as()]) + F::Wide::cast_from(other.0);
		let offset = F::Wide::cast_from(log & ONEMASK) + (log >> FIELD_BITS);
		Additive(EXP_TABLE[offset.cast_as()])
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	pub fn mul_assign_slice(selfy: &mut [Self], other: Multiplier<F>) {
		for s in selfy {
			*s = s.mul(other);
		}
	}
}


/// Multiplicaiton friendly LOG form of f2e16
#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, PartialEq, Eq)] // Default, PartialOrd,Ord
pub struct Multiplier<F: FieldT>(pub F::Element);

impl<F: FieldT> From<usize> for Multiplier<F> {
	fn from(src: usize) -> Self {
		Self(F::Element::from(src))
	}
}

impl<F: FieldT> Multiplier<F> {
    #[inline(always)]
	pub fn to_wide(self) -> F::Wide {
		self.0 as F::Wide
	}

	#[inline(always)]
	pub fn from_wide(x: F::Wide) -> Multiplier<F> {
		Multiplier(x as F::Element)
	}

	pub fn as_usize(&self) -> Self {
		usize::from(self.0)
	}
}


/// Fast Walsh–Hadamard transform over modulo ONEMASK
pub fn walsh<F: FieldT>(data: &mut [Multiplier<F>], size: usize) {
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		let depart_no_next = depart_no << 1;
		while j < size {
			for i in j..(depart_no + j) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				let mask: F::Wide = F::ONEMASK.cast_as() as F::Wide;
				let tmp2: F::Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: F::Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				data[i] = Multiplier(((tmp1 & mask) + (tmp1 >> F::FIELD_BITS)).cast_as() as F::Element);
				data[i + depart_no] = Multiplier(((tmp2 & mask) + (tmp2 >> F::FIELD_BITS)).cast_as() as F::Element);
			}
			j += depart_no_next;
		}
		depart_no = depart_no_next;
	}
}
