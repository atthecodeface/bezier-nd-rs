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
    + num_traits::PrimInt
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
        + num_traits::PrimInt
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

pub trait SignedInt:
    Int + std::ops::Neg<Output = Self> + CarryingAdd + BorrowingSub + OverflowingMul
{
}
impl<T> SignedInt for T where
    T: Int + std::ops::Neg<Output = Self> + CarryingAdd + BorrowingSub + OverflowingMul
{
}

pub trait UsefulInt: SignedInt {
    type Unsigned: UsefulUInt;
    type Dbl: Int;
    const NB: usize;
    const NB_DBL: usize = Self::NB * 2;
    const SIGN_MASK: Self;
    const DBL_LOWER_MASK: Self::Dbl;
    const DBL_UPPER_MASK: Self::Dbl;
    const DBL_SIGN_MASK: Self::Dbl;
    fn as_dbl(self) -> Self::Dbl;
    fn as_dbl_upper(self) -> Self::Dbl;
    fn unsigned(self) -> Self::Unsigned;
    fn of_unsigned(v: Self::Unsigned) -> Self;
    #[inline]
    fn of_dbl_unchecked(dbl: Self::Dbl) -> Self {
        Self::of_unsigned(Self::unsigned_of_dbl(dbl).0)
    }
    fn unsigned_of_dbl(dbl: Self::Dbl) -> (Self::Unsigned, Self::Unsigned);
    fn dbl_mult(self, other: &Self) -> Self::Dbl;
}

pub trait UsefulUInt: Int {
    type Dbl: Int;
    const NB: usize;
    fn as_dbl(self) -> Self::Dbl;
    fn as_dbl_upper(self) -> Self::Dbl {
        self.as_dbl() << Self::NB
    }
}

macro_rules! make_useful {
    {$t:ty, $uns:ty, $dbl:ty, $dbl_uns:ty, $nb:expr} => {
    impl UsefulInt for $t {
        type Unsigned = $uns;
        type Dbl = $dbl;
        const NB: usize = $nb;
        const SIGN_MASK : Self = Self::ONE << ($nb-1);
        const DBL_LOWER_MASK : Self::Dbl = !Self::DBL_UPPER_MASK;
        const DBL_UPPER_MASK : Self::Dbl = (-Self::Dbl::ONE) << $nb;
        const DBL_SIGN_MASK : Self::Dbl = Self::Dbl::ONE << ($nb*2-1);
        fn as_dbl(self) -> Self::Dbl {
            self as $dbl
        }
        fn as_dbl_upper(self) -> Self::Dbl {
            self.as_dbl() << Self::NB
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
        fn unsigned_of_dbl(dbl:Self::Dbl) -> (Self::Unsigned, Self::Unsigned) {
            ((dbl >> $nb) as $uns, dbl as $uns)
        }
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
