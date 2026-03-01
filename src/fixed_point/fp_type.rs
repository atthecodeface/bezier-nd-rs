use super::Int;

use num_traits::{ConstOne, ConstZero};

/// A trait for how 'FpType<T,N>` can be a fixed point value with base type `T` and `N` fractional bits
///
/// The fixed point type must have one bit for sign, at least one bit for
/// integer (at least two for proper logs, three for proper trig constants), and
/// at least one bit for fractional part (otherwise a pure integer can be used)
pub trait HowIsFixedPoint<T: Int> {
    /// Number of bits in T
    const NB: usize = 0;
    /// Number of bits in the fractional part
    const NB_FRAC: usize = 0;
    /// Number of bits in the integer part (including the sign bit)
    ///
    /// This *MUST* conform to NB = NB_SIGN_AND_INT + NB_FRAC
    const NB_SIGN_AND_INT: usize = Self::NB - Self::NB_FRAC;
    /// Number of bits in the pure integer part (excluding the sign bit)
    ///
    /// This *MUST* conform to NB = 1 + NB_INT + NB_FRAC
    const NB_INT: usize = Self::NB - Self::NB_FRAC;
    /// Number of fractional bits in the double, when treated as our associated fixed point value
    ///
    /// The double has twice the number of integer bits and twice the number of fractional bits
    const NB_DBLFRAC: usize = Self::NB_FRAC * 2;
    /// Number of integer bits in the double, when treated as our associated fixed point value
    ///
    /// The double has twice the number of integer bits and twice the number of fractional bits
    const NB_DBLINT: usize = Self::NB_SIGN_AND_INT * 2;
    /// Bit mask of T of all zeros except the top bit (the sign bit)
    ///
    /// Must equal T::ONE << (Self::NB - 1)
    const SIGN_MASK: T;
    /// Bit mask of T of zeros in the fractional part, ones elsewhwere
    ///
    /// Must equal  (-T::ONE) << Self::NB_FRAC
    const SIGNED_INT_MASK: T;
    /// Bit mask of T of ones in the fractional part, zeros elsewhwere
    ///
    /// Must equal  !Self::SIGNED_INT_MASK
    const FRAC_MASK: T;
    /// Bit mask of a T with lower NB_INT-1 bits set
    ///
    /// If a value of T contains an integer where any of these bits are set but
    /// not all then the value cannot be represented by the fixed point type
    const MAX_INT_MASK: T;
    /// The value of one as a fixed point number
    ///
    /// Must equal T::ONE << Self::NB_FRAC;
    const ONE: T;
}

impl HowIsFixedPoint<i8> for FPType<i8, 4> {
    const NB: usize = 8;
    const NB_FRAC: usize = 4;
    const SIGN_MASK: i8 = i8::ONE << (Self::NB - 1);
    const SIGNED_INT_MASK: i8 = (-i8::ONE) << Self::NB_FRAC;
    const FRAC_MASK: i8 = !Self::SIGNED_INT_MASK;
    const MAX_INT_MASK: i8 = (!i8::ZERO) >> (Self::NB_FRAC + 1);
    const ONE: i8 = i8::ONE << Self::NB_FRAC;
}

/// Type for which `HowIsFixedPoint<T>` must be implemented for `Fixed<T,N>` to be valid
pub struct FPType<T, const N: usize>([T; N]);
