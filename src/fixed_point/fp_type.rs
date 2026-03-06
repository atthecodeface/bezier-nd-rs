use super::Int;

/// A trait for how 'FpType<T,N>` can be a fixed point value with base type `T` and `N` fractional bits
///
/// The fixed point type must have one bit for sign, at least one bit for
/// integer (at least two for proper logs, three for proper trig constants), and
/// at least one bit for fractional part (otherwise a pure integer can be used)
pub trait HowIsFixedPoint<T: Int> {
    /// Number of bits in T
    const NB: usize = 0;
    /// Number of bits in the fractional part
    ///
    /// This is given by 'N' in `Fixed<T,N>`; it is here because it is needed for the other constants
    const NB_FRAC: usize = 0;
    /// Number of bits in the integer part (including the sign bit)
    ///
    /// This *MUST* conform to NB = NB_SIGN_AND_INT + NB_FRAC
    const NB_SIGN_AND_INT: usize = Self::NB - Self::NB_FRAC;
    /// Number of bits in the pure integer part (excluding the sign bit)
    ///
    /// This *MUST* conform to NB = 1 + NB_INT + NB_FRAC for two's complement values
    ///
    /// This *MUST* conform to NB = NB_INT + NB_FRAC for values with a separate sign
    const NB_INT: usize = Self::NB - Self::NB_FRAC - 1;
    /// Number of fractional bits in the double, when treated as our associated fixed point value
    ///
    /// The double has twice the number of integer bits and twice the number of fractional bits
    const NB_DBLFRAC: usize = Self::NB_FRAC * 2;
    /// Number of integer bits in the double, when treated as our associated fixed point value
    ///
    /// The double has twice the number of integer bits and twice the number of fractional bits
    const NB_DBLINT: usize = Self::NB_SIGN_AND_INT * 2;
    /// The value of one as a fixed point number
    ///
    /// Must equal T::ONE << Self::NB_FRAC;
    const ONE: T;
    /// True if there is a dedicated sign bit and the value is held as an unsigned value; false if the value is held using two's complement
    const DEDICATED_SIGN: bool = false;
}

/// Type for which `HowIsFixedPoint<T>` must be implemented for `Fixed<T,N>` to be valid
pub struct FPType<T, const N: usize>([T; N]);
