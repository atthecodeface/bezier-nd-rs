use super::{FPType, Fixed, HowIsFixedPoint, UsefulConsts, UsefulInt};

use num_traits::FloatConst;

impl<T: UsefulInt + UsefulConsts, const N: usize> Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn cnst(c: T, int_bits: isize) -> Self {
        let shr =
            ((<FPType<T, N> as HowIsFixedPoint<T>>::NB_SIGN_AND_INT as isize) - int_bits) as usize;
        Self { value: c >> shr }
    }
}

impl<T: UsefulInt + UsefulConsts, const N: usize> FloatConst for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn SQRT_2() -> Self {
        // sqrt(2) = 1.41, so one integer bits
        Self::cnst(T::SQRT_2, 1)
    }
    fn FRAC_1_SQRT_2() -> Self {
        // 1/sqrt(2) = 0.707, so zero integer bits
        Self::cnst(T::SQRT_2, 0)
    }

    fn TAU() -> Self {
        assert!(<FPType<T, N> as HowIsFixedPoint<T>>::NB_INT > 1);
        // TAU is 6.28, so three integer bits
        Self::cnst(T::PI, 3)
    }
    fn PI() -> Self {
        // PI is 3.14, so two integer bits
        Self::cnst(T::PI, 2)
    }
    fn FRAC_PI_2() -> Self {
        Self::cnst(T::PI, 1)
    }
    fn FRAC_PI_4() -> Self {
        Self::cnst(T::PI, 0)
    }
    fn FRAC_PI_8() -> Self {
        Self::cnst(T::PI, -1)
    }

    fn FRAC_2_PI() -> Self {
        // 2/pi = 0.6366..., so zero integer bits
        Self::cnst(T::FRAC_2_PI, 0)
    }
    fn FRAC_1_PI() -> Self {
        // 1/pi = 0.3183..., so minus one integer bits
        Self::cnst(T::FRAC_2_PI, -1)
    }

    fn FRAC_2_SQRT_PI() -> Self {
        // 2/pi.sqrt() = 1.1283..., so one integer bits
        Self::cnst(T::FRAC_1_SQRT_PI, 1)
    }
    fn FRAC_PI_3() -> Self {
        // pi/3 = 1.047..., so one integer bits
        Self::cnst(T::FRAC_PI_3, 1)
    }
    fn FRAC_PI_6() -> Self {
        // pi/6 = 0.5235..., so zero integer bits
        Self::cnst(T::FRAC_PI_3, 0)
    }

    fn E() -> Self {
        // e is 2.7, so two integer bits
        Self::cnst(T::E, 2)
    }

    fn LN_10() -> Self {
        // ln(10) is 2.3025..., so two integer bits
        Self::cnst(T::LN_10, 2)
    }
    fn LN_2() -> Self {
        // ln(2) is 0.693..., so zero integer bits
        Self::cnst(T::LN_2, 0)
    }
    fn LOG10_E() -> Self {
        // log10(e) is 0.4342..., so minus one integer bits
        Self::cnst(T::LOG10_E, -1)
    }
    fn LOG2_E() -> Self {
        // log2(e) is 1,442..., so one integer bit
        Self::cnst(T::LOG2_E, 1)
    }
    fn LOG10_2() -> Self {
        // log10(2) is 0.3010..., so minus one integer bit
        Self::cnst(T::LOG10_2, -1)
    }
    fn LOG2_10() -> Self {
        // log2(10) is 3.321..., so two integer bits
        Self::cnst(T::LOG2_10, 2)
    }
}
