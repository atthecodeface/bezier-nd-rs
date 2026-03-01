use super::{FPType, Fixed, HowIsFixedPoint, UsefulInt};

use num_traits::{ConstOne, ConstZero, One, Zero};

impl<T: UsefulInt, const N: usize> Zero for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn is_zero(&self) -> bool {
        self == &Self::ZERO
    }
    fn zero() -> Self {
        Self::ZERO
    }
}

impl<T: UsefulInt, const N: usize> ConstZero for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    const ZERO: Self = Self { value: T::ZERO };
}

impl<T: UsefulInt, const N: usize> One for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn is_one(&self) -> bool {
        self == &Self::ONE
    }
    fn one() -> Self {
        Self::ONE
    }
}

impl<T: UsefulInt, const N: usize> ConstOne for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    const ONE: Self = Self {
        value: <FPType<T, N> as HowIsFixedPoint<T>>::ONE,
    };
}
