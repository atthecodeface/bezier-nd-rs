use super::UsefulInt;

use super::{ArithCode, FPType, Fixed, HowIsFixedPoint};

use num_traits::Bounded;

impl<T: UsefulInt, const N: usize> num_traits::CheckedAdd for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn checked_add(&self, v: &Self) -> Option<Self> {
        let mut result = *self;
        let ac = result.do_add(v);
        ac.checked(result)
    }
}

impl<T: UsefulInt, const N: usize> num_traits::CheckedSub for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn checked_sub(&self, v: &Self) -> Option<Self> {
        let mut result = *self;
        let ac = result.do_sub(v);
        ac.checked(result)
    }
}

impl<T: UsefulInt, const N: usize> num_traits::CheckedMul for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn checked_mul(&self, v: &Self) -> Option<Self> {
        let mut result = *self;
        let ac = result.do_mul(v);
        ac.checked(result)
    }
}

impl<T: UsefulInt, const N: usize> num_traits::CheckedDiv for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn checked_div(&self, v: &Self) -> Option<Self> {
        let mut result = *self;
        let ac = result.do_div(v);
        ac.checked(result)
    }
}

impl<T: UsefulInt, const N: usize> num_traits::CheckedRem for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn checked_rem(&self, v: &Self) -> Option<Self> {
        let mut result = *self;
        let ac = result.do_rem(v);
        ac.checked(result)
    }
}
