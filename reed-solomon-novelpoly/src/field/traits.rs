use core::ops::{BitXor, BitXorAssign, Mul, MulAssign};
use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign, AsRef};

/// Additive field representation
pub trait FieldAdd:
{
    const FIELD_BITS: usize;
    const FIELD_BYTES: usize = Self::FIELD_BITS / 8;
    const FIELD_SIZE: usize = 1_usize << Self::FIELD_BITS;
    const FIELD_NAME: &'static str;

    type ElementAsBytes;
    type Element : Sized;
    type Wide;
    // const ONE: Self;

    const ONEMASK: Self::Element = (Self::FIELD_SIZE - 1) as Self::Element;

    /// Quotient ideal generator given by tail of irreducible polynomial
    const GENERATOR: Self::Element;

    /// Cantor basis' final element
    const BASE_FINAL: Self::Element;

    //fn get_internals(&self) -> Elm;
    fn get_field_name() -> str;

    fn generate_cantor_basis(mut next: Self::Element) -> Option<[Self::Element; Self::FIELD_BITS]> {
        let mut base: [Self::Element; Self::FIELD_BITS] = [0; Self::FIELD_BITS];
        base[0] = 1;
        for b in base.iter_mut().rev() {
            if next == 0 || (next == 1 && *b != 1) { return None; }
            *b = next;
            let square = gf_mul_bitpoly_reduced(next,next);
            next ^= square;
        }
        if next == 0 { Some(base) } else { None }
    }
    const ZERO: Self;

    fn to_be_bytes(&self) -> [u8; Self::FIELD_BYTES];
    fn from_be_bytes(serialized_element: [u8; Self::FIELD_BYTES]) -> Self;
    fn to_wide(self) -> Self::Wide;
    fn from_wide(x: Self::Wide) -> Self;

}

trait FieldAddElement<F: FieldAdd> : Clone + Copy + core::fmt::Debug + Default + PartialEq<Self> + Eq + BitXor<Self, Output = Self> + BitXorAssign<Self> {
}

/// Declare field Element
///
/// Additive via XOR form
///    
#[derive(Clone, Copy, Debug, Default, BitXor, BitXorAssign, PartialEq, Eq)] // PartialOrd,Ord
pub struct Additive<Elt> (pub Elt);

impl<Elt> FieldAdd for Additive<Elt> {

    const ZERO: Self = Additive::<Elt>(0);
    type ElementAsBytes = [u8; { Self::FIELD_BYTES }];

    // const ONE: Additive = Additive(1);
    fn to_be_bytes(&self) -> [u8; Self::FIELD_BYTES] {
        self.0.to_be_bytes()
    }

    fn from_be_bytes(serialized_element: [u8; Self::FIELD_BYTES]) -> Self {
        Self(Self::Element::from_be_bytes(serialized_element))
    }

    // impl AsRef<[u8; <Self as FieldAdd>::FIELD_BYTES]> for Additive {
    //     fn as_ref(&self) -> &[u8; <Self as FieldAdd>::FIELD_BYTES] {
    //         unimplemented!()
    //     }
    // }

    #[inline(always)]
    fn to_wide(self) -> Self::Wide {
        self.0 as Self::Wide
    }

    #[inline(always)]
    fn from_wide(x: Self::Wide) -> Self {
        Additive(x as Self::Element)
    }
}

/// Paramaterized field multiplier representation
pub trait FieldMul<F, Multiplier>: Mul<Multiplier, Output = Self> + MulAssign<Multiplier>
where
    F: FieldAdd,
    [u8; F::FIELD_BYTES] : Sized,
    Multiplier: Clone + Copy,
{
    fn to_multiplier(self) -> Multiplier;

    /// Multiply field elements by a single multiplier, using SIMD if available
    ///
    /// We avoid using `where [Self]: MulAssign<Multiplier>` here because
    /// providing a default trait impl requires specialization.
    #[inline(always)]
    fn mul_assign_slice(selfy: &mut [Self], other: Multiplier) where Self: Sized{
        for s in selfy {
            *s = (*s) * (other);
        }
    }
}


/// Fast Walsh–Hadamard transform over modulo ONEMASK
pub fn walsh<F: FieldAdd>(data: &mut [Logarithm<F>], size: usize) {
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

/// Multiplicaiton friendly LOG form of f2e16
#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, PartialEq, Eq)] // Default, PartialOrd,Ord
pub struct Logarithm<F: FieldAdd>(pub F::Element);

impl <F>  Logarithm<F> where F: FieldAdd {
    #[inline(always)]
	pub fn to_wide(self) -> F::Wide {
		self.0 as F::Wide
	}
    #[inline(always)]
	pub fn from_wide(x: F::Wide) -> Logarithm<F> {
		Logarithm(x as F::Element)
	}
}

impl <F> Mul<Logarithm<F>> for Logarithm<F> where F: FieldAdd, [u8; F::FIELD_BYTES] : Sized
 {
    type Output = Additive<F>;

	/// TODO:  Leaky abstraction!  Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
    #[cfg(table_bootstrap_complete)]
    fn mul(self, other: Logarithm) -> Additive {
		let log = self.0 as Wide + other.0 as Wide;
        // Compute sum of logarithms modulo 2^FIELD_BITS-1 perhaps? 
		let offset = (log & ONEMASK as Wide) + (log >> FIELD_BITS);
		Additive(EXP_TABLE[offset as usize])
    }

    #[cfg(not(table_bootstrap_complete))]
    fn mul(self, other: Logarithm<F>) -> Additive<F> { panic!(); }
}

impl <F: FieldAdd> Mul<Logarithm<F>> for Additive<F> where [u8; F::FIELD_BYTES]: Sized {
    type Output = Additive<F>;

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
    #[cfg(table_bootstrap_complete)]
    fn mul(self, other: Logarithm<F>) -> Additive<F> {
		if self == Self::ZERO {
			return Self::ZERO;
		}
        self.to_multiplier() * other
    }

    #[cfg(not(table_bootstrap_complete))]
    fn mul(self, other: Logarithm<F>) -> Additive<F> { panic!(); }
}

impl <F: FieldAdd> MulAssign<Logarithm<F>> for Additive<F> where [u8; F::FIELD_BYTES]: Sized {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Logarithm<F>) {
        *self = (*self) * rhs;
    }
}

#[cfg(table_bootstrap_complete)]
impl FieldMul<Logarithm> for Additive {
	/// Return multiplier prepared form
    #[inline(always)]
	fn to_multiplier(self) -> Logarithm {
		Logarithm(LOG_TABLE[self.0 as usize])
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	fn mul_assign_slice(selfy: &mut [Self], other: Logarithm) {
		for s in selfy {
			*s *= other;
		}
	}
}



#[cfg(table_bootstrap_complete)]
include!(concat!(env!("OUT_DIR"), "/table_", $name, ".rs"));
