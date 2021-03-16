
// Do not change. Must stay as is for buid script.
use super::FieldT;
use super::AdditiveT;
use super::MultiplierT;
use super::ops;

#[derive(Debug, Clone, Copy)]
pub struct Field;

impl FieldT for Field {
	const NAME: &'static str = "f2e16";
	type Element = u16;
	type Wide = u32;
	type Additive = Additive;
	type Multiplier = Multiplier;
	const GENERATOR: Self::Element = 0x2D;
	const FIELD_BITS: usize = 16;

	const FIELD_SIZE: usize = 1_usize << Self::FIELD_BITS;
    const ONEMASK: Self::Element = (Self::FIELD_SIZE - 1) as Self::Element;

	// cantor base
	const BASE: &'static [Self::Element] =
	&[1 as Self::Element, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198];
}



use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign};

/// Additive via XOR form of f2e16
#[derive(Clone, Copy, Debug, Default, BitXor, BitXorAssign, PartialEq, Eq)] // PartialOrd,Ord
pub struct Additive(pub u16);

impl From<u32> for Additive {
	fn from(src: u32) -> Self {
		Self(src as u16)
	}
}

impl From<u16> for Additive {
	fn from(src: u16) -> Self {
		Self(src)
	}
}

impl From<usize> for Additive {
	fn from(src: usize) -> Self {
		Self(src as u16)
	}
}

impl AdditiveT<Field> for Additive {
    #[inline(always)]
	fn to_wide(self) -> u32 {
		self.0 as u32
	}

    #[inline(always)]
	fn from_wide(x: u32) -> Additive {
		Additive(x as u16)
	}

	const ZERO: Additive = Additive(0 as u16);

	/// Return multiplier prepared form
    #[inline(always)]
	#[cfg(table_bootstrap_complete)]
	fn to_multiplier(self) -> Multiplier {
		Multiplier(LOG_TABLE[self.0 as usize])
	}

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
	#[cfg(table_bootstrap_complete)]
	fn mul(self, other: Multiplier) -> Additive {
		if self == Self::ZERO {
			return Self::ZERO;
		}
		let log = LOG_TABLE[self.0 as usize] as u32 + other.0;
		let offset = (log & ONEMASK) as u32 + (log >> FIELD_BITS);
		Additive(EXP_TABLE[offset as usize])
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	#[cfg(table_bootstrap_complete)]
	fn mul_assign_slice(selfy: &mut [Self], other: Multiplier) {
		for s in selfy {
			*s = s.mul(other);
		}
	}
}


/// Multiplicaiton friendly LOG form of f2e16
#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, PartialEq, Eq)] // Default, PartialOrd,Ord
pub struct Multiplier(pub u16);

impl From<usize> for Multiplier {
	fn from(src: usize) -> Self {
		Self(src as u16)
	}
}

impl From<u32> for Multiplier {
	fn from(src: u32) -> Self {
		Self(src as u16)
	}
}

impl From<u16> for Multiplier {
	fn from(src: u16) -> Self {
		Self(src)
	}
}

impl MultiplierT<Field> for Multiplier {
    #[inline(always)]
	fn to_wide(self) -> u32 {
		self.0 as u32
	}

	#[inline(always)]
	fn from_wide(x: u32) -> Multiplier {
		Multiplier(x as u16)
	}

	fn as_usize(&self) -> usize {
		self.0 as usize
	}
}
