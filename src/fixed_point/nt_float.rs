use super::{FPType, Fixed, HowIsFixedPoint, UsefulConsts, UsefulInt};

use num_traits::{ConstOne, ConstZero, Float, Zero};

impl<T: UsefulInt + UsefulConsts, const N: usize> Float for Fixed<T, N>
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
    fn min_positive_value() -> Self {
        Self { value: T::ONE }
    }
    #[inline]
    fn min_value() -> Self {
        <Self as num_traits::Bounded>::min_value()
    }
    #[inline]
    fn max_value() -> Self {
        <Self as num_traits::Bounded>::max_value()
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

    fn trunc(self) -> Self {
        let dedicated_sign = <FPType<T, N> as HowIsFixedPoint<T>>::DEDICATED_SIGN;
        let one = <FPType<T, N> as HowIsFixedPoint<T>>::ONE;
        let int_mask = (!T::ZERO) << N;
        if dedicated_sign || self.value >= T::ZERO {
            Self {
                value: self.value & int_mask,
            }
        } else if (self.value & int_mask) == self.value {
            self
        } else {
            // negative two's complement
            Self {
                value: (self.value & int_mask) + one,
            }
        }
    }

    fn floor(self) -> Self {
        let dedicated_sign = <FPType<T, N> as HowIsFixedPoint<T>>::DEDICATED_SIGN;
        let one = <FPType<T, N> as HowIsFixedPoint<T>>::ONE;
        let int_mask = (!T::ZERO) << N;
        if !dedicated_sign || self.value >= T::ZERO {
            Self {
                value: self.value & int_mask,
            }
        } else if (self.value & int_mask) == self.value {
            self
        } else {
            // negative non-integer with a dedicated sign
            Self {
                value: (self.value & int_mask) - one,
            }
        }
    }

    fn ceil(self) -> Self {
        let dedicated_sign = <FPType<T, N> as HowIsFixedPoint<T>>::DEDICATED_SIGN;
        let one = <FPType<T, N> as HowIsFixedPoint<T>>::ONE;
        let int_mask = (!T::ZERO) << N;
        if !dedicated_sign || self.value <= T::ZERO {
            Self {
                value: self.value & int_mask,
            }
        } else if (self.value & int_mask) == self.value {
            self
        } else {
            // positive non-integer with a dedicated sign
            Self {
                value: (self.value & int_mask) + one,
            }
        }
    }

    // This is trunc for values with a sign bit, floor for two's complement values
    fn round(self) -> Self {
        todo!();
    }

    fn fract(self) -> Self {
        let dedicated_sign = <FPType<T, N> as HowIsFixedPoint<T>>::DEDICATED_SIGN;
        let one = <FPType<T, N> as HowIsFixedPoint<T>>::ONE;
        let frac_mask = !((!T::ZERO) << N);
        if dedicated_sign || self.value >= T::ZERO {
            Self {
                value: self.value & frac_mask,
            }
        } else {
            // negative two's complement: fraction - 1
            Self {
                value: (self.value - one) & frac_mask,
            }
        }
    }

    #[inline]
    fn abs(self) -> Self {
        <Self as num_traits::Signed>::abs(&self)
    }

    #[inline]
    fn signum(self) -> Self {
        <Self as num_traits::Signed>::signum(&self)
    }

    #[inline]
    fn is_sign_positive(self) -> bool {
        !<Self as num_traits::Signed>::is_negative(&self)
    }

    #[inline]
    fn is_sign_negative(self) -> bool {
        <Self as num_traits::Signed>::is_negative(&self)
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
        <Self as num_traits::Signed>::abs_sub(&self, &other)
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }

    fn recip(self) -> Self {
        Self::ONE / self
    }

    fn sqrt(self) -> Self {
        // Invoke the fixed point function sqrt() on self (possibly shifted right by one) and shift left by frac_bits/2
        todo!();
    }

    fn cbrt(self) -> Self {
        todo!();
    }

    fn powi(self, n: i32) -> Self {
        <Self as num_traits::pow::Pow<_>>::pow(self, n)
    }

    fn powf(self, _n: Self) -> Self {
        // y^x = 2^(x/log2(y))
        todo!();
    }
    fn exp(self) -> Self {
        todo!();
    }
    fn exp2(self) -> Self {
        // 2^x = (2^x.frac())<<x.floor()
        todo!();
    }
    fn ln(self) -> Self {
        // ln(x) = ln(2)*log2(x)
        todo!();
    }
    fn log(self, _base: Self) -> Self {
        // log2(e^x) = x*log2(e)
        // e^x = 2^(log2(e^x)) = 2^(x.log2(e))
        // logb(x) = ln(x) * ln(b)
        todo!();
    }
    fn log2(self) -> Self {
        todo!();
    }
    fn log10(self) -> Self {
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

    fn integer_decode(self) -> (u64, i16, i8) {
        if self.is_zero() {
            (0, 0, 1)
        } else {
            let num_bits = <FPType<T, N> as HowIsFixedPoint<T>>::NB;
            let num_frac_bits = N as i16;
            let (sign, s) = if self.is_sign_negative() {
                (-1, -self.value)
            } else {
                (1, self.value)
            };
            if let Some(value) = s.to_u64() {
                return (value, -num_frac_bits, sign);
            }
            // Can optimize this with a binary search
            let mut top_bit: usize = 63;
            let mut x = s >> top_bit;
            for _ in 0..num_bits {
                if x <= T::ONE {
                    break;
                }
                x = x >> 1_usize;
                top_bit += 1;
            }

            // If the top bit is bit 79 then we need to shift that down to bit 63, ie. by 16
            //
            // This *must* (or else things are pear-shaped) be convertible to a u64
            let s = s >> (top_bit - 63);
            (
                s.to_u64().unwrap(),
                (top_bit as i16) - num_frac_bits - 63,
                sign,
            )
        }
    }

    fn epsilon() -> Self {
        Self { value: T::ONE }
    }
    fn is_subnormal(self) -> bool {
        false
    }
    fn to_degrees(self) -> Self {
        self * Self::cnst(T::TO_DEGREES, 6)
    }
    fn to_radians(self) -> Self {
        self * Self::cnst(T::TO_RADIANS, -5)
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
