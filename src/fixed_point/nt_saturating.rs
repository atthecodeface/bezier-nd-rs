use super::{UsefulInt, UsefulUInt};

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
