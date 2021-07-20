use core::ops::{BitXor, BitXorAssign, Mul, MulAssign, Add, Shl, Shr, Index};
use std::ops::{BitAnd, Sub};
use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign, BitAnd};
use num::Integer;
use core::convert::{TryFrom,TryInto, Into};

/// Additive field representation
pub trait FieldAdd : Clone + Copy + core::fmt::Debug + Default + PartialEq<Self> + Eq where
{    
    const FIELD_BITS: usize;
    const FIELD_BYTES: usize = Self::FIELD_BITS / 8;
    const FIELD_SIZE: usize = 1_usize << Self::FIELD_BITS;
    const FIELD_NAME: &'static str;
    
    type Element : Sized + Clone + Copy + core::fmt::Debug + Default + PartialEq<Self::Element> + Eq + BitXor<Self::Element, Output = Self::Element> + BitXorAssign<Self::Element> + TryFrom<Self::Wide> + Integer + core::ops::AddAssign +  core::ops::SubAssign + Shl<usize, Output = Self::Element> + Shr<usize, Output = Self::Element> +  Into<usize> + Into<Self::Wide> + TryFrom<usize> + Into<usize> + Default;
    type Wide: Copy + Integer + BitAnd<Self::Wide, Output = Self::Wide>  + Shl<usize, Output = Self::Wide> + Shr<usize, Output = Self::Wide> + BitXor<Self::Wide, Output = Self::Wide> + BitXorAssign<Self::Wide> + TryInto<Self::Element> + From<Self::Element>;
    // const ONE: Self;
    type BaseArray;
    type FieldTableArray : Index<usize, Output=Self::Element>; //TODO: just index over the field elemennts arrays are not that sexy.
    type LogWalshTable;
    type AfftSkewTable : Index<usize, Output=Logarithm<Self>>;

    /// Quotient ideal generator given by tail of irreducible polynomial
    const GENERATOR: Self::Element;

    /// Cantor basis' final element
    const BASE_FINAL: Self::Element;

    //fn get_internals(&self) -> Elm;
    //fn get_field_name() -> str;

    const ZERO_ELEMENT: Self::Element;
    const ONE_ELEMENT: Self::Element;
    const ONEMASK: Self::Element; // should be  = (F::FIELD_SIZE - 1);
    const ONEMASK_USIZE: usize; // should be  = (F::FIELD_SIZE - 1);
    

    const ZERO_ELEMENT_WIDE: Self::Wide;
    const ONE_ELEMENT_WIDE: Self::Wide;
    const ONEMASK_WIDE: Self::Wide; // should be  = (F::FIELD_SIZE - 1);

    const BASE: Self::BaseArray;
    const LOG_TABLE: Self::FieldTableArray;
    const EXP_TABLE: Self::FieldTableArray;
    const LOG_WALSH: Self::LogWalshTable;
    /// Additive FFT and IFFT skew table and formal derivative transform table
    /// for computing Reed-Solomon in the "novel polynomial basis".
    const AFFT_SKEW_TABLE: Self::AfftSkewTable;

    fn generate_cantor_basis(mut next: Self::Element) -> Option<[Self::Element; Self::FIELD_BITS]> where
        <<Self as FieldAdd>::Wide as TryInto<<Self as FieldAdd>::Element>>::Error: core::fmt::Debug,
    [Self::Element; Self::FIELD_BITS]: Sized,
    
    {
        let mut base: [Self::Element; Self::FIELD_BITS] = [Self::ZERO_ELEMENT; Self::FIELD_BITS];
        base[0] = Self::ONE_ELEMENT;
        for b in base.iter_mut().rev() {
            if next == Self::ZERO_ELEMENT || (next == Self::ONE_ELEMENT && *b != Self::ONE_ELEMENT) { return None; }
            *b = next;
            let square : Self::Element = gf_mul_bitpoly_reduced::<Self>(next,next) ;
            next = next ^ square as Self::Element;
        }
        if next == Self::ZERO_ELEMENT { Some(base) } else { None }
    }

    fn from_be_bytes_to_element(bytes: [u8; Self::FIELD_BYTES]) -> Self::Element;
    fn from_element_to_be_bytes(element: Self::Element) -> [u8; Self::FIELD_BYTES];

}

pub trait TruncateTo<F: FieldAdd> : From<F::Element> + TryInto<F::Element> + BitAnd<Output = Self>  {
    const TRUNCATE_MASK: Self;
    fn truncate(self) -> F::Element where <Self as TryInto<F::Element>>::Error: core::fmt::Debug {
        let truncated_wide = Self::TRUNCATE_MASK & self;
        truncated_wide.try_into().expect("already truncated. q.e.d")
    }

}

impl<F: FieldAdd> TruncateTo<F> for F::Wide  {
    const TRUNCATE_MASK: Self = F::ONEMASK_WIDE;
}


trait FieldAddElement<F: FieldAdd> : Clone + Copy + core::fmt::Debug + Default + PartialEq<Self> + Eq + BitXor<Self, Output = Self> + BitXorAssign<Self>  where
    <F::Wide as TryInto<<F as FieldAdd>::Element>>::Error : core::fmt::Debug,
{
    
    const ZERO: Self;
    //type ElementAsBytes: Sized;
}

/// Declare field Element
///
/// Additive via XOR form
///    
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)] // PartialOrd,Ord
pub struct Additive<F: FieldAdd> (pub F::Element);

impl <F: FieldAdd> BitXor for  Additive<F> {
    type Output = Self;
    fn bitxor(self, _: Self) -> <Self as BitXor<Self>>::Output { todo!() }
}

impl <F: FieldAdd> BitXorAssign for  Additive<F> {
    fn bitxor_assign(&mut self, _: Self) { todo!() }
}

impl<F: FieldAdd> FieldAddElement<F> for Additive<F> where
    F: FieldAdd,
<F::Wide as TryInto<<F as FieldAdd>::Element>>::Error : core::fmt::Debug,
{

    const ZERO: Self = Additive::<F>(F::ZERO_ELEMENT);
}

impl<F: FieldAdd> Additive<F> where
    F: FieldAdd,
<F::Wide as TryInto<<F as FieldAdd>::Element>>::Error : core::fmt::Debug,
{

    pub fn to_be_bytes(&self) -> [u8; F::FIELD_BYTES] {
        F::from_element_to_be_bytes(self.0)
    }

    pub fn from_be_bytes(serialized_element: [u8; F::FIELD_BYTES]) -> Self {
        Self(F::from_be_bytes_to_element(serialized_element))
    }

    #[inline(always)]
    pub fn to_wide(self) -> F::Wide {
        self.0.into()
    }

    #[inline(always)]
    pub fn from_wide(x: F::Wide) -> Additive<F> {
        Additive(TruncateTo::<F>::truncate(x))
    }
}

/// Paramaterized field multiplier representation
pub trait FieldMul<F, Multiplier>: Clone + Copy + Mul<Multiplier, Output = Self> + MulAssign<Multiplier>
where
    F: FieldAdd,
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
fn gf_mul_bitpoly_reduced<F: FieldAdd>(a: F::Element, b: F::Element) -> F::Element where <F::Wide as TryInto<<F as FieldAdd>::Element>>::Error : core::fmt::Debug {
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
		self.0.into()
	}
    #[inline(always)]
	pub fn from_wide(x: F::Wide) -> Logarithm<F> where <<F as FieldAdd>::Wide as TryInto<<F as FieldAdd>::Element>>::Error: core::fmt::Debug {
		Logarithm(TruncateTo::<F>::truncate(x))
	}
}

impl <F> Mul<Logarithm<F>> for Logarithm<F> where F: FieldAdd, [u8; F::FIELD_BYTES] : Sized, [(); F::FIELD_SIZE]: Sized, <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug,
 {
    type Output = Additive<F>;

	/// TODO:  Leaky abstraction!  Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
    #[cfg(table_bootstrap_complete)]
    fn mul(self, other: Logarithm<F>) -> Additive<F> where [F::Element; F::FIELD_SIZE]: Sized,
    [(); F::FIELD_SIZE]: Sized,
    <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug
    {
        let lhs: F::Wide =self.0.into();
        let rhs : F::Wide = other.0.into();
		let log : F::Wide = lhs + rhs;
        // Compute sum of logarithms modulo 2^FIELD_BITS-1 perhaps?
        let lhs: F::Element = TruncateTo::<F>::truncate(log);
        let rhs: F::Element = TruncateTo::<F>::truncate(log >> F::FIELD_BITS);
		let offset : F::Element = lhs + rhs;
		Additive(F::EXP_TABLE[offset.into()])
    }

    #[cfg(not(table_bootstrap_complete))]
    fn mul(self, other: Logarithm<F>) -> Additive<F> { panic!(); }
}

impl <F: FieldAdd> Mul<Logarithm<F>> for Additive<F>
        where  <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug,            
{
    type Output = Additive<F>;

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
    #[cfg(table_bootstrap_complete)]
    fn mul(self, other: Logarithm<F>) -> Additive<F> {
		if self == Self::ZERO {
			return Self::ZERO;
		}
         Additive::<F>(F::EXP_TABLE[(self.to_multiplier().0 + other.0).into()])
    }

    #[cfg(not(table_bootstrap_complete))]
    fn mul(self, other: Logarithm<F>) -> Additive<F> { panic!(); }
}

impl <F: FieldAdd> MulAssign<Logarithm<F>> for Additive<F>
where  <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug
{
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Logarithm<F>) {
        *self = (*self) * rhs;
    }
}

impl<F: FieldAdd> FieldMul<F, Logarithm<F>> for Additive<F>  where
    <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug
{
    //type Output=Additive<F>;
	/// Return multiplier prepared form
    #[inline(always)]
	fn to_multiplier(self) -> Logarithm<F> {
		Logarithm(F::LOG_TABLE[self.0.into()])
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	fn mul_assign_slice(selfy: &mut [Self], other: Logarithm<F>) {
		for s in selfy {
			*s *= other;
		}
	}
}


