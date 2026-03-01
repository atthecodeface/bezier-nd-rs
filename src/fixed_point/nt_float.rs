use super::{FPType, Fixed, HowIsFixedPoint, UsefulInt};

use num_traits::{ConstOne, ConstZero, Float, Zero};

impl<T: UsefulInt, const N: usize> Float for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn nan() -> Self {
        panic!("Fixed point numbers do not have a NaN");
    }
    fn infinity() -> Self {
        panic!("Fixed point numbers do not have infinity");
    }
    fn neg_infinity() -> Self {
        panic!("Fixed point numbers do not have infinity");
    }
    fn neg_zero() -> Self {
        Self::ZERO
    }
    fn min_value() -> Self {
        Self { value: T::ONE }
    }
    fn min_positive_value() -> Self {
        Self { value: T::ONE }
    }
    fn max_value() -> Self {
        let sign_bit = <FPType<T, N> as HowIsFixedPoint<T>>::SIGN_MASK;
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        Self {
            value: int_mask & !sign_bit,
        }
    }
    fn is_nan(self) -> bool {
        false
    }
    fn is_infinite(self) -> bool {
        false
    }
    fn is_finite(self) -> bool {
        true
    }
    fn is_normal(self) -> bool {
        !self.is_zero()
    }
    fn classify(self) -> std::num::FpCategory {
        if self.is_zero() {
            std::num::FpCategory::Zero
        } else {
            std::num::FpCategory::Normal
        }
    }
    fn floor(self) -> Self {
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        Self {
            value: self.value & int_mask,
        }
    }
    fn ceil(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        let int_one = <FPType<T, N> as HowIsFixedPoint<T>>::ONE;
        if self.value & frac_mask == T::ZERO {
            self
        } else {
            Self {
                // This should panic on overflow
                value: (self.value & int_mask) + int_one,
            }
        }
    }
    fn round(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        let half_frac_mask = frac_mask >> 1_u8;
        if self.value > T::ZERO {
            Self {
                value: (self.value + (half_frac_mask + T::ONE)) & frac_mask,
            }
        } else {
            Self {
                value: (self.value + half_frac_mask) & frac_mask,
            }
        }
    }
    fn trunc(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        // remove the fractional part; for +ve numbers this is clear the bits, for -ve numbers add the fractional mask the clear the bits
        if self.value < T::ZERO {
            Self {
                value: (self.value + frac_mask) & int_mask,
            }
        } else {
            Self {
                value: self.value & int_mask,
            }
        }
    }
    fn fract(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        Self {
            value: self.value & frac_mask,
        }
    }
    fn abs(self) -> Self {
        if self.value < T::ZERO {
            Self { value: -self.value }
        } else {
            Self { value: -self.value }
        }
    }
    fn signum(self) -> Self {
        if self.value < T::ZERO {
            -Self::ONE
        } else {
            Self::ONE
        }
    }
    fn is_sign_positive(self) -> bool {
        self.value >= T::ZERO
    }
    fn is_sign_negative(self) -> bool {
        self.value < T::ZERO
    }
    fn mul_add(self, _a: Self, _b: Self) -> Self {
        todo!();
    }
    fn recip(self) -> Self {
        todo!();
    }
    fn powi(self, _n: i32) -> Self {
        todo!();
    }
    fn powf(self, _n: Self) -> Self {
        todo!();
    }
    fn sqrt(self) -> Self {
        todo!();
    }
    fn exp(self) -> Self {
        todo!();
    }
    fn exp2(self) -> Self {
        todo!();
    }
    fn ln(self) -> Self {
        todo!();
    }
    fn log(self, _base: Self) -> Self {
        todo!();
    }
    fn log2(self) -> Self {
        todo!();
    }
    fn log10(self) -> Self {
        todo!();
    }
    fn max(self, other: Self) -> Self {
        if self.value > other.value {
            self
        } else {
            other
        }
    }
    fn min(self, other: Self) -> Self {
        if self.value < other.value {
            self
        } else {
            other
        }
    }
    fn abs_sub(self, other: Self) -> Self {
        if self.value > other.value {
            self - other
        } else {
            other - self
        }
    }
    fn cbrt(self) -> Self {
        todo!();
    }
    fn hypot(self, _other: Self) -> Self {
        todo!();
    }
    fn sin(self) -> Self {
        todo!();
    }
    fn cos(self) -> Self {
        todo!();
    }
    fn tan(self) -> Self {
        todo!();
    }
    fn asin(self) -> Self {
        todo!();
    }
    fn acos(self) -> Self {
        todo!();
    }
    fn atan(self) -> Self {
        todo!();
    }
    fn atan2(self, _other: Self) -> Self {
        todo!();
    }
    fn sin_cos(self) -> (Self, Self) {
        todo!();
    }
    fn exp_m1(self) -> Self {
        todo!();
    }
    fn ln_1p(self) -> Self {
        todo!();
    }
    fn sinh(self) -> Self {
        todo!();
    }
    fn cosh(self) -> Self {
        todo!();
    }
    fn tanh(self) -> Self {
        todo!();
    }
    fn asinh(self) -> Self {
        todo!();
    }
    fn acosh(self) -> Self {
        todo!();
    }
    fn atanh(self) -> Self {
        todo!();
    }
    fn integer_decode(self) -> (u64, i16, i8) {
        todo!();
    }

    fn epsilon() -> Self {
        Self { value: T::ONE }
    }
    fn is_subnormal(self) -> bool {
        false
    }
    fn to_degrees(self) -> Self {
        todo!();
    }
    fn to_radians(self) -> Self {
        todo!();
    }
    fn clamp(self, min: Self, max: Self) -> Self {
        if self < min {
            min
        } else if self > max {
            max
        } else {
            self
        }
    }
    fn copysign(self, sign: Self) -> Self {
        if (self.value < T::ZERO) == (sign.value < T::ZERO) {
            self
        } else {
            Self { value: -self.value }
        }
    }
}
