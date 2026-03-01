use super::{FPType, Fixed, HowIsFixedPoint, UsefulInt};

use num_traits::{ConstOne, ConstZero, FromPrimitive};

impl<T: UsefulInt, const N: usize> crate::NumOps for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    /// Create a value from a fraction with signed numerator and unsigned denominator
    fn frac(n: i32, u: u32) -> Self {
        Self::from_i32(n).unwrap() / Self::from_u32(u).unwrap()
    }

    /// Create a value from an [i32]
    fn of_i32(n: i32) -> Self {
        Self::from_i32(n).unwrap()
    }

    /// Convert an f64 to this value
    fn of_f64(_v: f64) -> Self {
        Self::default()
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn reciprocal(self) -> Self {
        Self::ONE / self
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn is_unreliable_divisor(self) -> bool {
        false
    }

    fn powi(self, p: i32) -> Self {
        (0..p).fold(Self::ONE, |acc, _| acc * self)
    }

    /// Return an estimate of the square root to within a precision
    fn sqrt_est(self) -> Self {
        // sqrt with a u32_24 input produces a u32_28 output
        // ; a u64_0 (pure integer) produces a u64_32 output
        let nb = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        if (nb & 1) == 0 {
            let x_uns_half_frac = super::functions::sqrt(self.value.unsigned());
            Self {
                value: T::of_unsigned(x_uns_half_frac) >> (nb / 2),
            }
        } else {
            let x_uns_half_frac = super::functions::sqrt(self.value.unsigned() >> 1_u8);
            Self {
                value: T::of_unsigned(x_uns_half_frac) >> (nb / 2),
            }
        }
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn cbrt_est(self) -> Self {
        crate::utils::cbrt_est::<_, 5>(self)
    }
    /// Return true if the value is negative
    fn is_sign_negative(self) -> bool {
        self < Self::ZERO
    }
}
