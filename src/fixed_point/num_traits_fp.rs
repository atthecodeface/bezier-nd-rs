use num_traits::{ConstOne, ConstZero, One, Zero};

use super::FixedPoint;
use super::FixedPoint_i32_16;

impl Zero for FixedPoint_i32_16 {
    fn is_zero(&self) -> bool {
        self == &Self::ZERO
    }
    fn zero() -> Self {
        Self::ZERO
    }
}
impl ConstZero for FixedPoint_i32_16 {
    const ZERO: Self = Self::of_i8(0);
}

/*
impl One for FixedPoint_i32_16 {
    fn is_one(&self) -> bool {
        self == Self::ONE
    }
    fn one() -> Self {
        Self::ONE
    }
}
impl ConstOne for FixedPoint_i32_16 {
    const ONE: Self = Self::of_i8(1);
}
*/
