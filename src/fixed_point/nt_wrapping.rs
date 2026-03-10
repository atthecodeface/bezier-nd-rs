use super::UsefulInt;

use super::{FPType, Fixed, HowIsFixedPoint};

impl<T: UsefulInt, const N: usize> num_traits::WrappingAdd for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn wrapping_add(&self, v: &Self) -> Self {
        let mut result = *self;
        let _ = result.do_add(v);
        result
    }
}

impl<T: UsefulInt, const N: usize> num_traits::WrappingSub for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn wrapping_sub(&self, v: &Self) -> Self {
        let mut result = *self;
        let _ = result.do_sub(v);
        result
    }
}

impl<T: UsefulInt, const N: usize> num_traits::WrappingMul for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn wrapping_mul(&self, v: &Self) -> Self {
        let mut result = *self;
        let _ = result.do_mul(v);
        result
    }
}
