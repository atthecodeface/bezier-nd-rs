use super::{FPType, Fixed, HowIsFixedPoint, UsefulInt};

use num_traits::{Bounded, ConstOne, ConstZero, Float, Signed, Zero};

impl<T: UsefulInt, const N: usize> Signed for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn abs(&self) -> Self {
        if self.value < T::ZERO {
            -*self
        } else {
            *self
        }
    }
    fn signum(&self) -> Self {
        if self.value < T::ZERO {
            Self { value: -T::ONE }
        } else {
            Self { value: T::ONE }
        }
    }

    fn is_positive(&self) -> bool {
        self.value > T::ZERO
    }
    fn is_negative(&self) -> bool {
        self.value < T::ZERO
    }
    fn abs_sub(&self, other: &Self) -> Self {
        if self.value > other.value {
            *self - *other
        } else {
            *other - *self
        }
    }
}

impl<T: UsefulInt, const N: usize> Bounded for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn min_value() -> Self {
        Self {
            value: T::SIGN_MASK,
        }
    }
    fn max_value() -> Self {
        Self {
            value: !T::SIGN_MASK,
        }
    }
}
