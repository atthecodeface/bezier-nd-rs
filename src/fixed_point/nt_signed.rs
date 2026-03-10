use super::{FPType, Fixed, HowIsFixedPoint, UsefulInt};

use num_traits::{Bounded, Signed};

impl<T: UsefulInt, const N: usize> Signed for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    #[inline]
    fn abs(&self) -> Self {
        if self.value < T::ZERO {
            -*self
        } else {
            *self
        }
    }

    #[inline]
    fn signum(&self) -> Self {
        if self.value < T::ZERO {
            Self { value: -T::ONE }
        } else {
            Self { value: T::ONE }
        }
    }

    #[inline]
    fn is_positive(&self) -> bool {
        self.value > T::ZERO
    }

    #[inline]
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
    #[inline]
    fn min_value() -> Self {
        Self {
            value: T::min_value(),
        }
    }
    #[inline]
    fn max_value() -> Self {
        Self {
            value: T::max_value(),
        }
    }
}
