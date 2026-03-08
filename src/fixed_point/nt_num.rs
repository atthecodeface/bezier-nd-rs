use super::UsefulInt;

use super::{FPType, Fixed, HowIsFixedPoint};

// num_traits::PrimInt ?
// num_traits::Signed (not Unsigned)
// num_traits::Pow ?
// num_traits::Inv ?
// num_traits::Wrapping  Neg/Shl/Shr
// num_traits::Saturating
// num_traits::Checked Neg/Shl/Shr
// num_traits::ops::overflowing::*

impl<T: UsefulInt, const N: usize> num_traits::Num for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    type FromStrRadixErr = ();
    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Ok(Self::default())
    }
}
