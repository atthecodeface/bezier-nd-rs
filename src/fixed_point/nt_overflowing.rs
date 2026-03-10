use super::UsefulInt;

use super::{FPType, Fixed, HowIsFixedPoint};

impl<T: UsefulInt, const N: usize> num_traits::ops::overflowing::OverflowingAdd for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn overflowing_add(&self, v: &Self) -> (Self, bool) {
        let mut result = *self;
        let ac = result.do_add(v);
        (result, ac.overflowing("addition"))
    }
}

impl<T: UsefulInt, const N: usize> num_traits::ops::overflowing::OverflowingSub for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn overflowing_sub(&self, v: &Self) -> (Self, bool) {
        let mut result = *self;
        let ac = result.do_sub(v);
        (result, ac.overflowing("subtraction"))
    }
}

impl<T: UsefulInt, const N: usize> num_traits::ops::overflowing::OverflowingMul for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn overflowing_mul(&self, v: &Self) -> (Self, bool) {
        let mut result = *self;
        let ac = result.do_mul(v);
        (result, ac.overflowing("multiplication"))
    }
}
