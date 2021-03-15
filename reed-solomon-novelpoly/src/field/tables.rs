
macro_rules! inc_table {
    ($fieldname:tt) => {

        pub(crate) mod $fieldname {
            use super::super::$fieldname::$ft;

            #[cfg(not(table_bootstrap_complete))]
            pub(crate) const LOG_TABLE: [$ft::Element; $ft::FIELD_SIZE] = [0 as $ft::Element; $ft::FIELD_SIZE];

            // #[cfg(not(table_bootstrap_complete))]
            // pub(crate) const LOG_WALSH: [u16; FIELD_SIZE] = [0; FIELD_SIZE];

            #[cfg(not(table_bootstrap_complete))]
            pub(crate) const EXP_TABLE: [$ft::Element; $ft::FIELD_SIZE] = [0 as $ft::Element; $ft::FIELD_SIZE];

            #[cfg(table_bootstrap_complete)]
            include!(concat!(env!("OUT_DIR"), "/table_", stringify!($ft::NAME), ".rs"));
        }
    }
}


// inc_table!(f256);
inc_table!(f2e16);
