use core::ops::{BitXor, BitXorAssign, Mul, MulAssign, Add, Shl, Shr, Index};
use std::ops::{BitAnd, Sub};
use derive_more::{Add, AddAssign, Sub, SubAssign, Display};
use num::{PrimInt, NumCast};
use core::convert::{TryFrom,TryInto, Into};


#[cfg(table_bootstrap_complete)]
use super::gf_mul_bitpoly_reduced;

/// Additive finite field representation
pub trait FieldAdd : Clone + Copy + core::fmt::Debug + Default + PartialEq<Self> + Eq  + 'static
{    
    const FIELD_BITS: usize;
    const FIELD_BYTES: usize = Self::FIELD_BITS / 8;
    const FIELD_SIZE: usize = 1_usize << Self::FIELD_BITS;
    const FIELD_NAME: &'static str;
    
    type Element : PrimInt + Default + core::fmt::Debug + BitXorAssign<Self::Element> + TryFrom<Self::Wide> +  core::ops::AddAssign +  core::ops::SubAssign + Into<usize>  + Into<Self::Wide> + TryFrom<usize> +  From<bool> + NumCast;

    type Wide: PrimInt + core::fmt::Debug + BitXorAssign<Self::Wide> +  From<Self::Element>;

    // const ONE: Self;
    /// Quotient ideal generator given by tail of irreducible polynomial
    const GENERATOR: Self::Element;

    /// Cantor basis' final element
    const BASE_FINAL: Self::Element;

    //fn get_internals(&self) -> Elm;
    //fn get_field_name() -> str;

    const ZERO_ELEMENT: Self::Element;
    const ONE_ELEMENT: Self::Element;
    const ONEMASK: Self::Element; // should be  = (F::FIELD_SIZE - 1);
    const ONEMASK_USIZE: usize = Self::FIELD_SIZE - 1;
    
    const ZERO_ELEMENT_WIDE: Self::Wide;
    const ONE_ELEMENT_WIDE: Self::Wide;
    const ONEMASK_WIDE: Self::Wide; // should be  = (F::FIELD_SIZE - 1);

    const BASE: &'static [Self::Element];
    const LOG_TABLE: &'static [Self::Element];
    const EXP_TABLE: &'static [Self::Element];
    const LOG_WALSH: &'static [Logarithm<Self>];
    /// Additive FFT and IFFT skew table and formal derivative transform table
    /// for computing Reed-Solomon in the "novel polynomial basis".
    const AFFT_SKEW_TABLE: &'static [Logarithm<Self>];

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

    //because accessing trait tables carshes the compiler we are going to have these
    //auxilary functions
    fn get_base_table(index: usize) -> Self::Element;
    fn get_exp_table(index: usize) -> Self::Element;
    fn get_log_table(index: usize) -> Self::Element;
    fn get_log_walsh(index: usize) -> Logarithm<Self>;
    fn get_skew(i: usize) -> Logarithm<Self>;


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

/// Declare field Element
///
/// Additive via XOR form
///    
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Display)] // PartialOrd,Ord
pub struct Additive<F: FieldAdd> (pub F::Element);


impl <F: FieldAdd> BitXor for  Additive<F> {
    type Output = Self;
    fn bitxor(self, rhs: Self) -> <Self as BitXor<Self>>::Output { Additive::<F>(self.0 ^ rhs.0)  }
}

impl <F: FieldAdd> BitXorAssign for  Additive<F> {
    fn bitxor_assign(&mut self, rhs: Self) { self.0 ^= rhs.0; }
}

impl<F: FieldAdd> Additive<F> where
    F: FieldAdd,
<F::Wide as TryInto<<F as FieldAdd>::Element>>::Error : core::fmt::Debug,
{

    #[inline(always)]
    pub fn to_be_bytes(&self) -> [u8; F::FIELD_BYTES] {
        F::from_element_to_be_bytes(self.0)
    }

    #[inline(always)]
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

    #[inline(always)]
    pub fn zero() -> Self {
        Additive::<F>(F::ZERO_ELEMENT)
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

impl<F: FieldAdd> FieldMul<F, Logarithm<F>> for Additive<F>  where
    <F::Wide as TryInto<F::Element>>::Error : core::fmt::Debug
{
    //type Output=Additive<F>;
	/// Return multiplier prepared form
    #[inline(always)]
	fn to_multiplier(self) -> Logarithm<F> {
	    //Logarithm(F::LOG_TABLE[<F::Element as Into<usize>>::into(self.0)])
	    Logarithm(F::get_log_table(<F::Element as Into<usize>>::into(self.0)))
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	fn mul_assign_slice(selfy: &mut [Self], other: Logarithm<F>) {
		for s in selfy {
			*s *= other;
		}
	}
}

///
/// Multiplicaiton friendly LOG form
///
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
		Additive(F::EXP_TABLE[<F::Element as Into<usize>>::into(offset)])
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
		if self == Self::zero() {
			return Self::zero();
		}
        //TODO: Why don't we check for self.to_multiplier().0 + other.0 being too big
	//TODO: F::FIELD_SIZE - 1 a constant
         Additive::<F>(F::get_exp_table(<F::Element as Into<usize>>::into(self.to_multiplier().0) + <F::Element as Into<usize>>::into(other.0)))
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



