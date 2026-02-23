//! Fixed point values
mod binary_ops;
mod fixed_point_i32_i16;
mod num_traits_fp;
mod raw;
pub use fixed_point_i32_i16::FixedPoint_i32_16;

pub(crate) use raw::{SignedRaw3232, UnsignedRaw3232};

pub(crate) trait FPRaw: num_traits::Num + Copy + std::fmt::Debug {}

pub(crate) trait FPInt: num_traits::Num + Copy + std::fmt::Debug {}
pub(crate) trait FPFrac: num_traits::Num + Copy + std::fmt::Debug {}

impl<T> FPInt for T where T: num_traits::Num + Copy + std::fmt::Debug {}
impl<T> FPFrac for T where T: num_traits::Num + Copy + std::fmt::Debug {}

// Trait implemented by fixed point numbers
pub(crate) trait FixedPoint: Sized {
    /// Number of bits in the fractional part
    const FRAC_BITS: usize;
    /// True if signed, False if unsigned
    const SIGNED: bool;
    /// Number of bits in integer part (including sign bit)
    const INT_BITS: usize;
    const FRAC_MASK: u64 = (1_64 << Self::FRAC_BITS) - 1;
    const INT_MASK: u64 = !Self::FRAC_MASK;

    /// Type for pure integer part, signed if SIGNED is true
    ///
    /// Used for intermediate operations; must have at least 2 more bits than the actual integer
    type INT: FPInt;
    /// Unsigned type for pure fractional part; must have at least 2 more bits than the actual fraction
    type FRAC: FPFrac;
    /// Type used to store the value
    type RAW: FPInt;
    /// Intermediate operation integer type with at least double the number of integer bits
    type DBLINT: FPInt;
    /// Intermediate operation unsigned fractional type with at least double the number of fraction bits
    type DBLFRAC: FPFrac;

    /// Return the raw value
    fn as_raw(&self) -> Self::RAW;

    /// Set from a raw value
    fn of_raw(value: &Self::RAW) -> Self;

    /// Map to an intermediate value, integer part
    fn inter_int(&self) -> Self::DBLINT;

    /// Map to an intermediate value, fractional part
    fn inter_frac(&self) -> Self::DBLFRAC;

    /// Map to an intermediate value
    ///
    /// None if there is an overflow
    fn of_inter(value: &(Self::DBLINT, Self::DBLFRAC)) -> Option<Self>;

    // #[inline(always)]
    // fn integer(&self) -> Self::INT {
    //    (self.as_raw() >> Self::FRACT_BITS) as Self::INT
    // }
}
