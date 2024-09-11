pub mod ffi {
    #![allow(non_upper_case_globals)]
    #![allow(non_camel_case_types)]
    #![allow(non_snake_case)]
    #![allow(dead_code)]

    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

    // Define pointer to static char arrays defined in abPOA (align_seq.c)
    // Use a trick found at: https://github.com/rust-lang/rust/issues/54450
    // See https://play.rust-lang.org/?version=stable&mode=debug&edition=2015&gist=f67a2476f7743d87225e69ae6efd2910
    pub struct ExternPtr<T>(*const T);
    unsafe impl<T> Sync for ExternPtr<T> {}

    impl<T> ExternPtr<T> {
        pub fn as_ptr(&self) -> *const T {
            self.0
        }
    }

    /// Transforms AaCcGgTtNn to 0-4
    pub static AB_NT4_TABLE: ExternPtr<u8> = {
        extern "C" {
            #[link_name = "ab_nt4_table"]
            static inner: u8;
        }
        ExternPtr(unsafe { &inner as *const u8 })
    };

    /// Transforms 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
    pub static AB_NT256_TABLE: ExternPtr<u8> = {
        extern "C" {
            #[link_name = "ab_nt256_table"]
            static inner: u8;
        }
        ExternPtr(unsafe { &inner as *const u8 })
    };

    /// Transforms amino acids to 0-25
    pub static AB_AA26_TABLE: ExternPtr<u8> = {
        extern "C" {
            #[link_name = "ab_aa26_table"]
            static inner: u8;
        }
        ExternPtr(unsafe { &inner as *const u8 })
    };

    /// Transforms 0-25 to amino acids
    pub static AB_AA256_TABLE: ExternPtr<u8> = {
        extern "C" {
            #[link_name = "ab_aa26_table"]
            static inner: u8;
        }
        ExternPtr(unsafe { &inner as *const u8 })
    };

    /// Pre-computed ilog2 for [0-2^16)
    pub static AB_LOG_TABLE_65536: ExternPtr<u8> = {
        extern "C" {
            #[link_name = "ab_LogTable65536"]
            static inner: u8;
        }
        ExternPtr(unsafe { &inner as *const u8 })
    };

    /// Pre-computed bit flag table (16 bit)
    pub static AB_BIT_TABLE_16: ExternPtr<u8> = {
        extern "C" {
            #[link_name = "ab_bit_table16"]
            static inner: u8;
        }
        ExternPtr(unsafe { &inner as *const u8 })
    };
}
