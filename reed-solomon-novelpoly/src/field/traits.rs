use core::ops::{BitXor, BitXorAssign, Mul, MulAssign};

/// Additive field representation 
pub trait FieldAdd :
    Clone + Copy + core::fmt::Debug + Default
    + PartialEq<Self> + Eq
    + BitXor<Self, Output=Self> + BitXorAssign<Self>

{
    const FIELD_BITS: usize;
    const FIELD_BYTES: usize = Self::FIELD_BITS / 8;
    const FIELD_SIZE: usize = 1_usize << Self::FIELD_BITS;
    const ZERO: Self;

    // const ONE: Self;
    type ElementAsBytes; // = [u8; Self::FIELD_BYTES];

    //fn get_internals(&self) -> Elm;
    //fn to_be_bytes(&self) -> Self::ElementAsBytes;
}

/// Paramaterized field multiplier representation
pub trait FieldMul<Multiplier> :
    FieldAdd + Mul<Multiplier, Output=Self> + MulAssign<Multiplier>
where
    Multiplier: Clone+Copy
{
	fn to_multiplier(self) -> Multiplier;

	/// Multiply field elements by a single multiplier, using SIMD if available
    ///
    /// We avoid using `where [Self]: MulAssign<Multiplier>` here because
    /// providing a default trait impl requires specialization.
    #[inline(always)]
	fn mul_assign_slice(selfy: &mut [Self], other: Multiplier) {
		for s in selfy {
			*s = (*s) * (other);
		}
	}
}


/// Declare field and include its tables
///
/// Requires Elt and Wide be defined previosuly.
macro_rules! decl_field_additive {
	($name:literal, bits = $fbits:literal) => {

        use derive_more::{Add, AddAssign, BitXor, BitXorAssign, Sub, SubAssign, AsRef};
        use super::{FieldAdd,FieldMul};

        /// Additive via XOR form
        #[derive(Clone, Copy, Debug, Default, BitXor, BitXorAssign, PartialEq, Eq)] // PartialOrd,Ord
        pub struct Additive(pub Elt);

        impl FieldAdd for Additive {
            type ElementAsBytes = [u8; {Self::FIELD_BYTES}];
            
        	const FIELD_BITS: usize = $fbits;
        	const ZERO: Additive = Additive(0);
        	// const ONE: Additive = Additive(1);

            // fn get_internals(&self) -> Elt{
            //     self.0;
            // }
        }

        impl AsRef<[u8; <Self as FieldAdd>::FIELD_BYTES]> for Additive {
            fn as_ref(&self) -> &[u8; <Self as FieldAdd>::FIELD_BYTES] {
                unimplemented!()
            }

        }
        
        impl Additive {
            #[inline(always)]
        	pub fn to_wide(self) -> Wide {
        		self.0 as Wide
        	}
            #[inline(always)]
        	pub fn from_wide(x: Wide) -> Additive {
        		Additive(x as Elt)
        	}

	    pub fn to_be_bytes(&self) -> [u8; <Self as FieldAdd>::FIELD_BYTES] {
                 self.0.to_be_bytes() 
            }

	    pub fn from_be_bytes(serialized_element: [u8; <Self as FieldAdd>::FIELD_BYTES]) -> Self {
		Self(Elt::from_be_bytes(serialized_element))
	    }

        }

        pub const FIELD_BITS: usize = Additive::FIELD_BITS;
        pub const FIELD_SIZE: usize = Additive::FIELD_SIZE;

        pub const ONEMASK: Elt = (Additive::FIELD_SIZE - 1) as Elt;

		pub const FIELD_NAME: &'static str = $name;

		#[cfg(table_bootstrap_complete)]
		include!(concat!(env!("OUT_DIR"), "/table_", $name, ".rs"));

	};
} // macro_rules! decl_field_additive


