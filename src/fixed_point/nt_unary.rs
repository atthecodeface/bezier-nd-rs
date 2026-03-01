use super::{FPType, Fixed, HowIsFixedPoint, UsefulConsts, UsefulInt};
use std::ops::*;

use num_traits::FloatConst;

impl<T: UsefulInt, const N: usize> Neg for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self { value: -self.value }
    }
}
