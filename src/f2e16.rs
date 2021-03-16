
pub const FIELD_NAME: &'static str = "f2e16";

pub type Elt = u16;
pub type Wide = u32;

pub const FIELD_BITS: usize = 16;
pub const FIELD_SIZE: usize = 1_usize << FIELD_BITS;

include!("f_log_mul.rs");

#[cfg(table_bootstrap_complete)]
include!("afft.rs");

// must be placed in a separate file, such that the preproc never tries to eval OUT_DIR
// in env which does not exist in the build.rs case
#[cfg(table_bootstrap_complete)]
include!(concat!(env!("OUT_DIR"), "/table_", "f2e16", ".rs"));


/* Needs Cleanup  */

pub type GFSymbol = Elt;
pub const ONEMASK: Elt = (FIELD_SIZE - 1) as Elt;

/// Quotient ideal generator given by tail of irreducible polynomial
pub const GENERATOR: Elt = 0x2D; // x^16 + x^5 + x^3 + x^2 + 1

// impl Additive {
//     pub const ONE: Additive = Additive(???);
// }

// Cantor basis
pub const BASE: [Elt; FIELD_BITS] =
    [1_u16, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198];



// TODO: Abstract this test over 0x1002d and move it into f_log_mut.rs

fn bitpoly_mul16(mut a: Wide, mut b: Wide) -> Wide {
    let mut r: Wide =0;
    for i in 0..FIELD_BITS {
        if (b>>i) & 1 != 0 { r ^= (a<<i); }
    }
    r
}

fn gf_mul_0x1002d(mut a: Elt, mut b: Elt) -> Elt {
    use core::convert::TryInto;
    let len = FIELD_BITS;
    let mut r: Wide = bitpoly_mul16(a as Wide,b as Wide);
    let red: Wide = 0x1002d;
    for i in (len..=(len*2-1)).rev() {
        if r & (1<<i) != 0 { r ^= (red<<(i-len)); }
    }
    r.try_into().unwrap()
}

#[test]
fn cantor_basis() {
    for w in BASE.windows(2) {
        let b = w[1];
        let square = gf_mul_0x1002d(b,b);
        let a = w[0];
        // let eq = if a == (square ^ b) { "==" } else { "!=" };
        // println!("{:#b} {} {:#b}\n", a, eq, square ^ b);
        assert_eq!(a, square ^ b);
    }
}
    