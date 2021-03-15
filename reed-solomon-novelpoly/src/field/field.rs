use std::ops::*;
use std::fmt::Debug;

pub trait Primitive: 'static + num::ToPrimitive + num::FromPrimitive + Copy + Clone + AddAssign + Add + Sub + SubAssign + Shl + Shr + ShlAssign + ShrAssign + BitXor + BitXorAssign + BitAnd + BitAndAssign {
}

impl<T> Primitive for T where T: 'static + num::ToPrimitive + num::FromPrimitive + Debug + Copy + Clone + AddAssign + Add + Sub + SubAssign + Shl + Shr + ShlAssign + ShrAssign + BitXor + BitXorAssign + BitAnd + BitAndAssign {}

// trait AdditiveT<F: FieldT>: Clone + Copy + Debug + Default + BitXor + BitXorAssign + PartialEq + Eq {
// 	fn to_wide(self) -> F::Wide;
// 	fn from_wide(x: F::Wide) -> Additive<F>;

// 	const ZERO: Additive<F><F>;
// }

// trait MultiplierT<F: FieldT>: Clone + Copy + Debug + Default + BitXor + BitXorAssign + PartialEq + Eq {
// 	fn to_multiplier(self) -> Multiplier<F>;
// 	fn mul(self, other: Multiplier<F>) -> Additive<F>;
// 	fn mul_assign_slice(selfy: &mut [Self], other: Multiplier<F>);
// }

pub trait Castomat<Y> {
	fn cast_as(self) -> Y;
	fn cast_from(y: Y) -> Self;
}

macro_rules! castomat_impl {
	($x:ty as $y:ty as ..) => {
		castomat_impl!($x as $y);
		castomat_impl!($y as $x);
	};
	($x:ty as $y:ty) => {
		impl Castomat<$y> for $x {
			#[inline(always)]
			fn cast_as(self) -> $y {
				self as $y
			}

			#[inline(always)]
			fn cast_from(y: $y) -> $x {
				y as $x
			}
		}
	};
}

castomat_impl!(u16 as u32 as ..);
castomat_impl!(u16 as u64 as ..);
castomat_impl!(u32 as u64 as ..);
castomat_impl!(u16 as usize as ..);
castomat_impl!(u32 as usize as ..);
castomat_impl!(u64 as usize as ..);

pub trait FieldT: Debug + Sized {
	const NAME: &'static str;

	type Element: Castomat<Self::Wide> + Primitive;
	type Wide: Castomat<Self::Element> + Primitive;
	// type Additive: From<Self::Element> + From<Self::Wide> + AdditiveT<Self>;
	// type Multiplier: From<Self::Element> + From<Self::Wide> + MultiplierT<Self>;

    const GENERATOR: Self::Element;

	const FIELD_BITS: usize;
	const FIELD_SIZE: usize;
	const ONEMASK: Self::Element;
	const BASE: &'static [Self::Element];
}



#[cfg(test)]
fn verify_cantor_basis<F: Field>() {
	unimplemented!("TODO {}", F::NAME);
	/*
    for b in &BASE {
        let b = Additive(*b);
        let square = b.mul(b.to_multiplier());
        assert!(BASE.contains(& (square ^ b).0 ));
    }
    */

    // for (previous, current) in CANTOR_BASE.tuple_windows(2) {
	// 	if (b>>i) & != 0 {
	// 		r ^= a << i;
	// 	}
    //     let b = Additive(w[1]);
    //     let square = b.mul(b.to_multiplier());
    //     assert_eq!(a, (square ^ b).0);
    // }

	// assert_eq!()
}
