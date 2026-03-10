use super::UsefulInt;
use super::{FPType, Fixed, HowIsFixedPoint};

use num_traits::ConstOne;

// num_traits::PrimInt - no, this is not an integer
// num_traits::Pow ?
// num_traits::Inv ?
// num_traits::Wrapping  Neg/Shl/Shr
// num_traits::Checked Neg/Shl/Shr

impl<T: UsefulInt, const N: usize> num_traits::Num for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    type FromStrRadixErr = ();
    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Ok(Self::default())
    }
}

macro_rules! pow_of_unsigned {
    {$t:ty} => {
        impl<T: UsefulInt, const N: usize> num_traits::pow::Pow<$t> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Self;
            fn pow(self, rhs: $t) -> Self::Output {
                num_traits::pow(self, rhs as usize)
            }
        }
    };
}
macro_rules! pow_of_signed {
    {$t:ty} => {
        impl<T: UsefulInt, const N: usize> num_traits::pow::Pow<$t> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Self;
            fn pow(self, rhs: $t) -> Self::Output {
                if rhs<0 {
                    Self::ONE / num_traits::pow(self, (-rhs) as usize)
                }else {
                num_traits::pow(self, rhs as usize)}
            }
        }
    };
}

pow_of_unsigned!(u8);
pow_of_unsigned!(u16);
pow_of_unsigned!(u32);
pow_of_unsigned!(u64);
pow_of_unsigned!(usize);

pow_of_signed!(i8);
pow_of_signed!(i16);
pow_of_signed!(i32);
pow_of_signed!(i64);
pow_of_signed!(isize);
