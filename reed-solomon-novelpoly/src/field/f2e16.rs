
// Do not change. Must stay as is for buid script.
use super::FieldT;
use super::ops;

#[derive(Debug, Clone, Copy)]
pub struct Field;

impl FieldT for Field {
	const NAME: &'static str = "f2e16";
	type Element = u16;
	type Wide = u32;
	// type Additive = ops::Additive<Self>;
	// type Multiplier = ops::Multiplier<Self>;
	const GENERATOR: Self::Element = 0x2D;
	const FIELD_BITS: usize = 16;

	const FIELD_SIZE: usize = 1_usize << Self::FIELD_BITS;
    const ONEMASK: Self::Element = (Self::FIELD_SIZE - 1) as Self::Element;

	// cantor base
	const BASE: &'static [Self::Element] =
	&[1 as Self::Element, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198];
}
