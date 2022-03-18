/// Declare field and include its tables
///
macro_rules! decl_field_additive {
 	($name: tt, bits = $fbits:literal, generator = $generator: literal, elt = $elt:tt, wide = $wide:tt, cantor_base_final_elt = $cantor_base_final_elt: literal) => {

        #[cfg(table_bootstrap_complete)]
        include!(concat!(env!("OUT_DIR"), "/table_", stringify!($name), ".rs"));

        #[derive(Clone, Copy,Debug, Default, PartialEq, Eq)]
        pub struct $name;
        impl FieldAdd for $name where
        {
            type Element = $elt;
            type Wide = $wide;
            const FIELD_BITS : usize = $fbits;

            /// Quotient ideal generator given by tail of irreducible polynomial
            const GENERATOR: Self::Element = $generator;

            const FIELD_NAME: &'static str = stringify!($name);

            const ZERO_ELEMENT: Self::Element = 0 as $elt;
            const ONE_ELEMENT: Self::Element = 1 as $elt;
            const ZERO_ELEMENT_WIDE: Self::Wide  = 0 as $wide;
            const ONE_ELEMENT_WIDE: Self::Wide = 1 as $wide;
            const ONEMASK: Self::Element = (Self::FIELD_SIZE - 1) as $elt;
            const ONEMASK_USIZE: usize = (Self::FIELD_SIZE - 1) as usize;
            const ONEMASK_WIDE: Self::Wide = (Self::FIELD_SIZE - 1) as $wide;

            const BASE_FINAL: Self::Element = $cantor_base_final_elt;

            fn from_be_bytes_to_element(bytes: [u8; Self::FIELD_BYTES]) -> Self::Element {
                Self::Element::from_be_bytes(bytes)
            }

            fn from_element_to_be_bytes(element: Self::Element) -> [u8; Self::FIELD_BYTES] {
                element.to_be_bytes()
            }

            #[cfg(table_bootstrap_complete)]
            const BASE: &'static [Self::Element] = &BASE;
            #[cfg(not(table_bootstrap_complete))]            
            const BASE: &'static [Self::Element] = &[Self::ZERO_ELEMENT; Self::FIELD_BITS];

            #[cfg(table_bootstrap_complete)]
            const LOG_TABLE: &'static [Self::Element] = &LOG_TABLE;
            #[cfg(not(table_bootstrap_complete))]            
            const LOG_TABLE: &'static [Self::Element] = &[Self::ZERO_ELEMENT; Self::FIELD_SIZE];

            #[cfg(table_bootstrap_complete)]
            const EXP_TABLE: &'static [Self::Element] = &EXP_TABLE;
            #[cfg(not(table_bootstrap_complete))]            
            const EXP_TABLE: &'static [Self::Element] = &[Self::ZERO_ELEMENT; Self::FIELD_SIZE];

            #[cfg(table_bootstrap_complete)]
            const LOG_WALSH: &'static [Logarithm<Self>] = &LOG_WALSH;
            #[cfg(not(table_bootstrap_complete))]            
            const LOG_WALSH: &'static [Logarithm<Self>] = &[Logarithm(Self::ZERO_ELEMENT); Self::FIELD_SIZE];          

            #[cfg(table_bootstrap_complete)]
            const AFFT_SKEW_TABLE: &'static [Logarithm<Self>] = &AFFT_SKEW_TABLE;
            #[cfg(not(table_bootstrap_complete))]
            /// Logarithm form of twisted factors used in our additive FFT
            const AFFT_SKEW_TABLE: &'static [Logarithm<Self>] = &[Logarithm(Self::ZERO_ELEMENT); Self::ONEMASK_USIZE] ;            

	    fn get_base_table(index: usize) -> Self::Element {
		#[cfg(table_bootstrap_complete)]
		return BASE[index];
			
		#[cfg(not(table_bootstrap_complete))]
		return Self::ZERO_ELEMENT;
	    }
	    
	    fn get_log_table(index: usize) -> Self::Element {
		#[cfg(table_bootstrap_complete)]
		return LOG_TABLE[index];
			
		#[cfg(not(table_bootstrap_complete))]
		return Self::ZERO_ELEMENT;
		
	    }

	    fn get_exp_table(index: usize) -> Self::Element {
		#[cfg(table_bootstrap_complete)]
		return EXP_TABLE[index  % (Self::FIELD_SIZE - 1)];
			
		#[cfg(not(table_bootstrap_complete))]
		return Self::ZERO_ELEMENT;
		
	    }

	    fn get_log_walsh(index: usize) -> Logarithm<Self> {
		#[cfg(table_bootstrap_complete)]
		return LOG_WALSH[index];

		#[cfg(not(table_bootstrap_complete))]
		return Logarithm(Self::ZERO_ELEMENT);

	    }

	    fn get_skew(i: usize) -> Logarithm<Self> {
		#[cfg(table_bootstrap_complete)]
		return AFFT_SKEW_TABLE[i];
			
		#[cfg(not(table_bootstrap_complete))]
		return Logarithm(Self::ZERO_ELEMENT);
		
	    }


        }


 	};

        
 } // macro_rules! decl_field_additive

#[cfg(test)]
#[macro_export]
macro_rules! test_all_fields_for {
        // Arguments are module name and function name of function to test bench
        ($fn_name:ident) => {
            // The macro will expand into the contents of this block.
            paste::item! {
                #[test]
                fn [< test_ $fn_name F256>] () {
                    $fn_name::<F256>();
                }

                #[test]
                fn [< test_ $fn_name F2e16>] () {
                    $fn_name::<F2e16>();
                }
            }
        };
 }//macro_rules! tess_all_fields_for
