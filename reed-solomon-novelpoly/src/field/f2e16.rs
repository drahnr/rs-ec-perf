use crate::{FieldAdd, Logarithm};

decl_field_additive!(F2e16, bits = 16, generator = 0x2D, elt = u16, wide = u32, cantor_base_final_elt = 39198);

#[cfg(table_bootstrap_complete)]
use crate::AfftField;
#[cfg(table_bootstrap_complete)]
impl AfftField for F2e16 {}


// #[cfg(table_bootstrap_complete)]
// include!("inc_afft.rs");

#[test]
fn embedded_gf256() {
    // We've a leaky to_multiplier abstraction that fucks up zero, so start at 1.
    //println!("{:?}", <F2e16 as FieldAdd>::LOG_TABLE);
    let mask: <F2e16 as FieldAdd>::Element = !0xFF;
     for i in 1..=255 {
	#[cfg(table_bootstrap_complete)]
     	 let i = Additive::<F2e16>(i as <F2e16 as FieldAdd>::Element);
         let log_i = <Additive<F2e16> as FieldMul<F2e16, Logarithm<F2e16>>>::to_multiplier(i);
         for j in 0..256 {
             let j = Additive::<F2e16>(j as <F2e16 as FieldAdd>::Element);
             assert!(<Additive<F2e16> as Mul<Logarithm<F2e16>>>::mul(j,log_i).0 != 257);
          }
    }
}

#[test]
fn flt_roundtrip_small() {
    const N: usize = 16;
    const EXPECTED: [Additive<F2e16>; N] =
        unsafe { std::mem::transmute([1_u16, 2, 3, 5, 8, 13, 21, 44, 65, 0, 0xFFFF, 2, 3, 5, 7, 11]) };

    let mut data = EXPECTED.clone();

    AfftField::afft(&mut data, N, N / 4);

    /*
    println!("novel basis(rust):");
    data.iter().for_each(|sym| {
        print!(" {:04X}", sym.0);
    });
    println!("");
    */
    
    AfftField::inverse_afft(&mut data, N, N / 4);

    assert_eq!(data, EXPECTED);

}
