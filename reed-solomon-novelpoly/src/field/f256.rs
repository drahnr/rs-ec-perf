#[cfg(table_bootstrap_complete)]
use super::*;
use core::convert::{TryInto};
use crate::{FieldAdd, TruncateTo, Logarithm, walsh};

decl_field_additive!(F256, bits = 8, generator = 0x1D, elt = u8, wide = u16, cantor_base_final_elt = 230);

#[cfg(table_bootstrap_complete)]
use crate::AfftField;
#[cfg(table_bootstrap_complete)]
impl AfftField for F256 {}

//    const AFFT_TABLES: AfftTables<F> = AfftTables::<F>::initialize();

// Valid Cantor basis, but fails embedded_gf16
//pub const GENERATOR: Elt = 0x1D; //GF(2^8): x^8 + x^4 + x^3 + x^2 + 1
// Valid Cantor basis, but fails both embedded_gf16 and b_is_one.
// pub const GENERATOR: Elt = 0x1B; //GF(2^8): x^8 + x^4 + x^3 + x + 1
// pub const GENERATOR: Elt = 0x3F; //GF(2^8): x^8 + x^5 + x^4 + x^3 + x^2 + x + 1
// pub const GENERATOR: Elt = 0x39; //GF(2^8): x^8 + x^5 + x^4 + x^3 + 1
// pub const GENERATOR: Elt = 0x77; //GF(2^8): z^8 + z^6 + z^5 + z^4 + z^3 + 1

// Is this Chen's suggested tower?  Does not yield a Cantor basis.
// pub const GENERATOR: Elt = 0x7B; //GF(2^8): z^8 + z^6 + z^5 + z^4 + z^3 + z + 1

// Select for GFNI compatability, but lacks an embedded GF(16).
// pub const GENERATOR: Elt = 0x1B; //GF(2^8): x^8 + x^4 + x^3 + x + 1

// impl Additive {
//     pub const ONE: Additive = Additive(???);
// }

// pub const BASE_FINAL: Elt = 238;

// /// Cantor basis
// pub const BASE: [Elt; FIELD_BITS] = [1, 214, 152, 146, 86, 200, 88, 230];

//include!("inc_logarithm.rs");


//#[cfg(table_bootstrap_complete)]
//include!("inc_afft.rs");
