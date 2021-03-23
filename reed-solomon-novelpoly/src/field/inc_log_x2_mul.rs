


/// Multiplicaiton friendly double LOG form for degree two extension field
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SplitLogarithm {
    /// Log form of coefficient of X^0 and X^1
    pub xs: [HalfLogarithm; 2],
    /// Log form of mixed term
    pub beta_x1_plus_x0: HalfLogarithm,
}

impl Mul<Logarithm> for Additive {
    type Output = Additive;

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
    #[cfg(table_bootstrap_complete)]
    fn mul(self, other: Logarithm) -> Additive {

        mul_prepared( (x T + y), (v_log,u_log, u_plus_v_log) )  -> (z, xu_plus_xv)
        where
           x_log = LOG[x]
           y_log = LOG[y]
           xu_plus_xv = EXP[x_log + u_log] ^ EXP[y_log + u_log]
           z = EXP[w_plus_v_log + x]


		if self == Self::ZERO {
			return Self::ZERO;
		}
        let s = self.split();
        let s = [ s[0].to_multiplier(), s[1].to_multiplier() ];
        
		let log = (LOG_TABLE[self.0 as usize] as Wide) + other.0 as Wide;
		let offset = (log & ONEMASK as Wide) + (log >> FIELD_BITS);
		Additive(EXP_TABLE[offset as usize])
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
        let xs = self.split();
        let beta_x1_plus_x0 = (xs[1].mul_beta() ^ xa[0]);
        let xs = [ xs[0].to_multiplier(), xs[1].to_multiplier() ];
		SplitLogarithm { xs, beta_x1_plus_x0 }
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







use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign};


/// Additive via XOR form of f2e16
#[derive(Clone, Copy, Debug, Default, BitXor, BitXorAssign, PartialEq, Eq)] // PartialOrd,Ord
pub struct Additive {
    pub adds: [Half::Additive; 2],
}

impl Additive {
    /*
    #[inline(always)]
	pub fn to_wide(self) -> Wide {
		self.0 as Wide
	}
    #[inline(always)]
	pub fn from_wide(x: Wide) -> Additive {
		Additive(x as Half)
	}
    */

	pub const ZERO: Additive = Additive { adds: [Half::Additive::ZERO; 2] };
}

#[cfg(table_bootstrap_complete)]
impl Additive {
	/// Return multiplier prepared form
    #[inline(always)]
	pub fn to_multiplier(self) -> Logarithm {
        let muls = [ self.adds[0].to_multiplier(), self.adds[0].to_multiplier() ];
        let mul_of_2nd_b_xor_1st = (muls[1].mul(Half::OMEGA) ^ muls[0]).to_multiplier();
		DoubleLogarithm { muls, mul_of_2nd_b_xor_1st }
	}

    #[inline(always)]
	pub fn mul(self, other: DoubleLogarithm) -> Additive {
		if self == Self::ZERO {
			return Self::ZERO;
		}
        let muls = [ self.adds[0].to_multiplier(), self.adds[0].to_multiplier() ];
        Additive { adds: [
            muls[0].mul(other.muls[0]) ^ muls[1].mul(other.muls[1]),
            mul_of_2nd_b_xor_1st.mul(other.muls[1]),
        ] }
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	pub fn mul_assign_slice(selfy: &mut [Self], other: DoubleLogarithm) {
		// TODO: SIMD
		for s in selfy {
			*s = s.mul(other);
		}
	}
}


/// Multiplicaiton friendly LOG form of f2e16
#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, PartialEq, Eq)] // Default, PartialOrd,Ord
pub struct DoubleLogarithm {
    pub muls: [Half::Logarithm; 2],
    pub mul_of_2nd_b_xor_1st: Logarithm,
}

// impl DoubleLogarithm { }
    