use num_traits::{ConstOne, ConstZero, Num};
pub trait Int:
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
    + Num
    + ConstOne
    + ConstZero
    + std::ops::AddAssign<Self>
    + std::ops::SubAssign<Self>
    + std::ops::MulAssign<Self>
    + std::ops::BitAndAssign<Self>
    + std::ops::BitOrAssign<Self>
    + std::ops::BitXorAssign<Self>
    + std::ops::MulAssign<Self>
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
        + Num
        + ConstOne
        + ConstZero
        + std::ops::AddAssign<Self>
        + std::ops::SubAssign<Self>
        + std::ops::MulAssign<Self>
        + std::ops::BitAndAssign<Self>
        + std::ops::BitOrAssign<Self>
        + std::ops::BitXorAssign<Self>
        + std::ops::Not<Output = Self>
        + std::ops::Shl<usize, Output = Self>
        + std::ops::Shr<usize, Output = Self>
        + std::ops::Shr<u8, Output = Self>
{
}

pub trait UsefulInt: Int + std::ops::Neg<Output = Self> {
    type Unsigned: UsefulUInt;
    type Dbl: Int;
    const NB: usize;
    const SIGN_MASK: Self;
    const DBL_LOWER_MASK: Self::Dbl;
    const DBL_UPPER_MASK: Self::Dbl;
    const DBL_SIGN_MASK: Self::Dbl;
    fn as_dbl(self) -> Self::Dbl;
    fn as_dbl_upper(self) -> Self::Dbl;
    fn unsigned(self) -> Self::Unsigned;
    fn of_unsigned(v: Self::Unsigned) -> Self;
    fn unsigned_of_dbl(dbl: Self::Dbl) -> (Self::Unsigned, Self::Unsigned);
    fn of_dbl(dbl: Self::Dbl) -> Option<Self>;
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
        const DBL_LOWER_MASK : Self::Dbl = (Self::Dbl::ONE << $nb) -Self::Dbl::ONE;
        const DBL_UPPER_MASK : Self::Dbl = !Self::DBL_LOWER_MASK;
        const DBL_SIGN_MASK : Self::Dbl = Self::DBL_UPPER_MASK | (Self::Dbl::ONE << ($nb-1));
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
        fn of_dbl(dbl: Self::Dbl) -> Option<Self> {
            let mut sgn = dbl;
            sgn &= Self::DBL_SIGN_MASK;
            match sgn {
                Self::DBL_SIGN_MASK => Some(dbl as $t),
                0 => Some(dbl as $t),
                _ => None
            }
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
