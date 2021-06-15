
pub use core::ops::{Mul, MulAssign};

include!("inc_cantor_basis.rs");

/// Multiplicaiton friendly LOG form of f2e16
#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, PartialEq, Eq)] // Default, PartialOrd,Ord
pub struct Logarithm(pub Elt);

impl Logarithm {
    #[inline(always)]
	pub fn to_wide(self) -> Wide {
		self.0 as Wide
	}
    #[inline(always)]
	pub fn from_wide(x: Wide) -> Logarithm {
		Logarithm(x as Elt)
	}
}

impl Mul<Logarithm> for Logarithm {
    type Output = Additive;

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
    fn mul(self, other: Logarithm) -> Additive { panic!(); }
}

impl Mul<Logarithm> for Additive {
    type Output = Additive;

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
    #[cfg(table_bootstrap_complete)]
    fn mul(self, other: Logarithm) -> Additive {
		if self == Self::ZERO {
			return Self::ZERO;
		}
        self.to_multiplier() * other
    }

    #[cfg(not(table_bootstrap_complete))]
    fn mul(self, other: Logarithm) -> Additive { panic!(); }
}

impl MulAssign<Logarithm> for Additive {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Logarithm) {
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



/*
Actually our to_multiplier abstraction is leaky

#[test]
fn multiply_by_zero() {
    let zero_mul = Additive::ZERO.to_multiplier();
    for i in 0..FIELD_SIZE {
        let i = Additive(i as Elt);
        // assert_eq!(Additive::ZERO, Additive::ZERO.mul(i.to_multiplier()) );
        assert_eq!(Additive::ZERO, i.mul(zero_mul) );
    }
}
*/

#[test]
fn embedded_gf16() {
    // We've a leaky to_multiplier abstraction that fucks up zero, so start at 1.
    let mask: Elt = !0xF;
    for i in 1..16 {
        let i = Additive(i as Elt).to_multiplier();
        for j in 0..16 {
            let j = Additive(j as Elt);
            assert!(j.mul(i).0 & mask == 0);
        }
    }
}


/*
#[test]
fn print_gf256() {
    use std::io::Write;
	let mut w = std::fs::OpenOptions::new().create(true).truncate(true).write(true).open(Additive::FIELD_NAME).unwrap();

    write!(w, "\n\n\n{} :\n", Additive::FIELD_NAME);
    for i in 1..=255 {
        write!(w, "{:#b} * .. = ", i);
        let i = Additive(i).to_multiplier();
        for j in 0..=255 {
            let j = Additive(j);
            write!(w, "{:#b} ", j.mul(i).0);
        }
        write!(w, "\n");
    }    
}
*/
