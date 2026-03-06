use super::{BorrowingSub, CarryingAdd, OverflowingMul};
use num_traits::{ConstOne, ConstZero};

/// Internal trait that bundles up traits required for UsefulInt and UsefulUInt
pub(crate) trait Int:
    Copy
    + 'static
    + Default
    + PartialOrd
    + Ord
    + Send
    + Sync
    + std::fmt::Debug
    + std::fmt::Display
    + std::fmt::LowerHex
    + num_traits::Num
    + num_traits::NumCast
    + num_traits::ToPrimitive
    + num_traits::FromPrimitive
    // + num_traits::PrimInt
    + num_traits::ConstOne
    + num_traits::ConstZero
    + std::ops::AddAssign<Self>
    + std::ops::SubAssign<Self>
    + std::ops::MulAssign<Self>
    + std::ops::BitAndAssign<Self>
    + std::ops::BitOrAssign<Self>
    + std::ops::BitXorAssign<Self>
    + std::ops::MulAssign<Self>
    + std::ops::BitAnd<Output = Self>
    + std::ops::BitOr<Output = Self>
    + std::ops::Not<Output = Self>
    + std::ops::Shl<usize, Output = Self>
    + std::ops::Shr<usize, Output = Self>
    + std::ops::Shr<u8, Output = Self>
{
}
impl<T> Int for T where
    T: Copy
        + 'static
        + Default
        + PartialOrd
        + Ord
        + Send
        + Sync
        + std::fmt::Debug
        + std::fmt::Display
        + std::fmt::LowerHex
        + num_traits::Num
        + num_traits::NumCast
        //        + num_traits::PrimInt
        + num_traits::FromPrimitive
        + ConstOne
        + ConstZero
        + std::ops::AddAssign<Self>
        + std::ops::SubAssign<Self>
        + std::ops::MulAssign<Self>
        + std::ops::BitAndAssign<Self>
        + std::ops::BitOrAssign<Self>
        + std::ops::BitXorAssign<Self>
        + std::ops::BitAnd<Output = Self>
        + std::ops::BitOr<Output = Self>
        + std::ops::Not<Output = Self>
        + std::ops::Shl<usize, Output = Self>
        + std::ops::Shr<usize, Output = Self>
        + std::ops::Shr<u8, Output = Self>
{
}

pub(crate) trait SignedInt:
    Int
    + std::ops::Neg<Output = Self>
    + num_traits::ops::overflowing::OverflowingSub
    + num_traits::ops::overflowing::OverflowingAdd // + CarryingAdd + BorrowingSub + OverflowingMul
{
}
impl<T> SignedInt for T where
    T: Int
        + std::ops::Neg<Output = Self>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub // + CarryingAdd + BorrowingSub + OverflowingMul
{
}

pub(crate) trait UsefulInt: SignedInt {
    type Unsigned: UsefulUInt;
    type Dbl: Int;
    type DblUnsigned: Int;
    const NB: usize;
    const NB_DBL: usize = Self::NB * 2;
    fn as_dbl(self) -> Self::Dbl;
    fn as_dbl_upper(self) -> Self::Dbl;
    fn as_dbl_unsigned(self) -> Self::DblUnsigned;
    fn unsigned(self) -> Self::Unsigned;
    fn of_unsigned(v: Self::Unsigned) -> Self;
    fn of_dbl(dbl: Self::Dbl) -> (Self, bool);
    // fn unsigned_of_dbl(dbl: Self::Dbl) -> (Self::Unsigned, Self::Unsigned);
    fn dbl_mult(self, other: &Self) -> Self::Dbl;
    fn min_value() -> Self;
    fn max_value() -> Self;
}

pub trait UsefulUInt: Int {
    type Dbl: Int;
    const NB: usize;
    fn as_dbl(self) -> Self::Dbl;
}

macro_rules! make_useful {
    {$t:ty, $uns:ty, $dbl:ty, $dbl_uns:ty, $nb:expr} => {
        impl UsefulInt for $t {
            type Unsigned = $uns;
            type Dbl = $dbl;
            type DblUnsigned = $dbl_uns;
            const NB: usize = $nb;
            fn as_dbl(self) -> Self::Dbl {
                self as $dbl
            }
            fn as_dbl_upper(self) -> Self::Dbl {
                self.as_dbl() << Self::NB
            }
            fn as_dbl_unsigned(self) -> Self::DblUnsigned {
                self as $dbl_uns
            }
            fn unsigned(self) -> Self::Unsigned {
                self as $uns
            }
            fn of_unsigned(v: Self::Unsigned) -> Self {
                v as Self
            }
            fn dbl_mult(self, other: &Self) -> Self::Dbl {
                (self as $dbl) * (*other as $dbl)
            }
            // Get a value from the least significant half of the Dbl
            #[inline]
            fn of_dbl(dbl: Self::Dbl) -> (Self, bool) {
                // Consider the upper $nb+1 bits - if they are all one, or all zero, then there is no overflow
                let upper_bits = dbl >> ($nb-1);
                let overflow = !(upper_bits==0) || (upper_bits+1==0);
                (Self::of_unsigned(dbl as $uns), overflow)
            }
            fn min_value() -> Self { ((1 as $uns)<<($nb-1)) as $t}
            fn max_value() -> Self { (((1 as $uns)<<($nb-1))-1) as $t }
            // fn unsigned_of_dbl(dbl:Self::Dbl) -> (Self::Unsigned, Self::Unsigned) {
            //     ((dbl >> $nb) as $uns, dbl as $uns)
            // }
        }
        impl UsefulUInt for $uns {
            type Dbl = $dbl_uns;
            const NB: usize = $nb;
            fn as_dbl(self) -> Self::Dbl {
                self as $dbl_uns
            }
        }
    };
}

make_useful!(i8, u8, i16, u16, 8);
make_useful!(i16, u16, i32, u32, 16);
make_useful!(i32, u32, i64, u64, 32);
make_useful!(i64, u64, i128, u128, 64);
