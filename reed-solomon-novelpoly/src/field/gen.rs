/// Declare field and include its tables
macro_rules! decl_field {
	($name:literal, bits = $fbits:literal) => {
		pub const FIELD_NAME: &'static str = $name;

		pub const FIELD_BITS: usize = $fbits;
		pub const FIELD_SIZE: usize = 1_usize << FIELD_BITS;

		#[cfg(table_bootstrap_complete)]
		include!(concat!(env!("OUT_DIR"), "/table_", $name, ".rs"));
	};
}



