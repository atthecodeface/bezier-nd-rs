use super::UsefulInt;

use super::{ArithCode, FPType, Fixed, HowIsFixedPoint};

use num_traits::Bounded;

impl<T: UsefulInt, const N: usize> num_traits::SaturatingAdd for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn saturating_add(&self, v: &Self) -> Self {
        let mut result = *self;
        match result.do_add(v) {
            ArithCode::OverflowMax => Self::max_value(),
            ArithCode::OverflowMin => Self::min_value(),
            _ => result,
        }
    }
}

impl<T: UsefulInt, const N: usize> num_traits::SaturatingSub for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn saturating_sub(&self, v: &Self) -> Self {
        let mut result = *self;
        match result.do_sub(v) {
            ArithCode::OverflowMax => Self::max_value(),
            ArithCode::OverflowMin => Self::min_value(),
            _ => result,
        }
    }
}

impl<T: UsefulInt, const N: usize> num_traits::SaturatingMul for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn saturating_mul(&self, v: &Self) -> Self {
        let mut result = *self;
        match result.do_mul(v) {
            ArithCode::OverflowMax => Self::max_value(),
            ArithCode::OverflowMin => Self::min_value(),
            _ => result,
        }
    }
}

impl<T: UsefulInt, const N: usize> num_traits::Saturating for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn saturating_add(self, v: Self) -> Self {
        let mut result = self;
        match result.do_add(&v) {
            ArithCode::OverflowMax => Self::max_value(),
            ArithCode::OverflowMin => Self::min_value(),
            _ => result,
        }
    }
    fn saturating_sub(self, v: Self) -> Self {
        let mut result = self;
        match result.do_sub(&v) {
            ArithCode::OverflowMax => Self::max_value(),
            ArithCode::OverflowMin => Self::min_value(),
            _ => result,
        }
    }
}
