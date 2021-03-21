use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign};

/// Additive via XOR form of f2e16
#[derive(Clone, Copy, Debug, Default, BitXor, BitXorAssign, PartialEq, Eq)] // PartialOrd,Ord
pub struct Additive(pub Elt);

impl Additive {
    #[inline(always)]
	pub fn to_wide(self) -> Wide {
		self.0 as Wide
	}
    #[inline(always)]
	pub fn from_wide(x: Wide) -> Additive {
		Additive(x as Elt)
	}

	pub const ZERO: Additive = Additive(0);
}

#[cfg(table_bootstrap_complete)]
impl Additive {
	/// Return multiplier prepared form
    #[inline(always)]
	pub fn to_multiplier(self) -> Multiplier {
		Multiplier(LOG_TABLE[self.0 as usize])
	}

	/// Return a*EXP_TABLE[b] over GF(2^r)
    #[inline(always)]
	pub fn mul(self, other: Multiplier) -> Additive {
		if self == Self::ZERO {
			return Self::ZERO;
		}
		let log = (LOG_TABLE[self.0 as usize] as Wide) + other.0 as Wide;
		let offset = (log & ONEMASK as Wide) + (log >> FIELD_BITS);
		Additive(EXP_TABLE[offset as usize])
	}

	/// Multiply field elements by a single multiplier, using SIMD if available
    #[inline(always)]
	pub fn mul_assign_slice(selfy: &mut [Self], other: Multiplier) {
		for s in selfy {
			*s = s.mul(other);
		}
	}
}


/// Multiplicaiton friendly LOG form of f2e16
#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, PartialEq, Eq)] // Default, PartialOrd,Ord
pub struct Multiplier(pub Elt);

impl Multiplier {
    #[inline(always)]
	pub fn to_wide(self) -> Wide {
		self.0 as Wide
	}
    #[inline(always)]
	pub fn from_wide(x: Wide) -> Multiplier {
		Multiplier(x as Elt)
	}
}


/// Fast Walsh–Hadamard transform over modulo ONEMASK
pub fn walsh(data: &mut [Multiplier], size: usize) {
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		let depart_no_next = depart_no << 1;
		while j < size {
			for i in j..(depart_no + j) {
				// We deal with data in log form here, but field form looks like:
				//			 data[i] := data[i] / data[i+depart_no]
				// data[i+depart_no] := data[i] * data[i+depart_no]
				let mask = ONEMASK as Wide;
				let tmp2: Wide = data[i].to_wide() + mask - data[i + depart_no].to_wide();
				let tmp1: Wide = data[i].to_wide() + data[i + depart_no].to_wide();
				data[i] = Multiplier(((tmp1 & mask) + (tmp1 >> FIELD_BITS)) as Elt);
				data[i + depart_no] = Multiplier(((tmp2 & mask) + (tmp2 >> FIELD_BITS)) as Elt);
			}
			j += depart_no_next;
		}
		depart_no = depart_no_next;
	}
}

#[allow(unused)]
fn bitpoly_mul(mut a: Wide, mut b: Wide) -> Wide {
    let mut r: Wide =0;
    for i in 0..FIELD_BITS {
        if (b>>i) & 1 != 0 {
			r ^= a<<i;
		}
    }
    r
}

#[allow(unused)]
fn gf_mul_bitpoly_reduced(a: Elt, b: Elt) -> Elt {
    use core::convert::TryInto;
    let len = FIELD_BITS;
    let mut r: Wide = bitpoly_mul(a as Wide,b as Wide);
    let red : Wide = (1 << FIELD_BITS) + (GENERATOR as Wide);
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

fn generate_cantor_basis(mut next: Elt) -> Option<[Elt; FIELD_BITS]> {
    let mut base: [Elt; FIELD_BITS] = [0; FIELD_BITS];
    base[0] = 1;
    for b in base.iter_mut().rev() {
        if next == 0 || (next == 1 && *b != 1) { return None; }
        *b = next;
        let square = gf_mul_bitpoly_reduced(next,next);
        next ^= square;
    }
    if next == 0 { Some(base) } else { None }
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
	let mut w = std::fs::OpenOptions::new().create(true).truncate(true).write(true).open(FIELD_NAME).unwrap();

    write!(w, "\n\n\n{} :\n", FIELD_NAME);
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

