use num_traits::{ConstOne, ConstZero, Num};
pub trait OverflowingMul: Sized {
    fn overflowing_mul(self, rhs: Self) -> (Self, bool);
}
pub trait CarryingAdd: Sized {
    /// Calculate the subtraction of rhs from self, and a further subtraction of the borrow if set, returning the wrapped value and a boolean indication of overflow
    ///
    /// Since the bool in the result indicates *overflow* not carry out, it must not be used as a carry in to other operations
    fn carrying_add(self, rhs: Self, carry: bool) -> (Self, bool);
}
pub trait BorrowingSub: Sized {
    /// Calculate the sum of self and rhs and a carry in, returning the wrapped value and a boolean indication of overflow
    ///
    /// Since the bool in the result indicates *overflow* not carry out, it must not be used as a carry in to other operations
    fn borrowing_sub(self, rhs: Self, carry: bool) -> (Self, bool);
}

macro_rules! make_signed_int_deps {
    {$t:ty} => {
    impl OverflowingMul for $t {
        fn overflowing_mul(self, rhs: Self) -> (Self, bool) {
            self.overflowing_mul(rhs)
        }
    }
    impl CarryingAdd for $t {
        fn carrying_add(self, rhs: Self, carry: bool) -> (Self, bool) {
            let (r, o) = self.overflowing_add(rhs);
            if !carry {
                (r, o)
            } else {
                let (r2, o2) = r.overflowing_add(1);
                (r2, o | o2)
            }
        }
    }
    impl BorrowingSub for $t {
        fn borrowing_sub(self, rhs: Self, borrow: bool) -> (Self, bool) {
            let (r, o) = self.overflowing_sub(rhs);
            if !borrow {
                (r, o)
            } else {
                let (r2, o2) = r.overflowing_sub(1);
                (r2, o | o2)
            }
        }
    }
    }
}

make_signed_int_deps!(i8);
make_signed_int_deps!(i16);
make_signed_int_deps!(i32);
make_signed_int_deps!(i64);
make_signed_int_deps!(isize);
make_signed_int_deps!(i128);

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
        + Num
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

/// ```ignore
/// const PI: Self = Self::TAU >> 1;
/// const FRAC_PI_2: Self = Self::PI >> 1;
/// const FRAC_PI_4: Self = Self::FRAC_PI_2 >> 1;
/// const FRAC_PI_6: Self = Self::FRAC_PI_3 >> 1;
/// const FRAC_PI_8: Self = Self::FRAC_PI_4 >> 1;
/// const FRAC_1_PI: Self = Self::FRAC_2_PI >> 1;
/// ```
pub trait UsefulConsts {
    const E: Self;
    const LN_2: Self;
    const LN_10: Self;
    const LOG2_E: Self;
    const LOG2_10: Self;
    const LOG10_2: Self;
    const LOG10_E: Self;
    const SQRT_2: Self;
    const TAU: Self;

    const FRAC_2_PI: Self;
    const FRAC_PI_3: Self;
    const FRAC_1_SQRT_2: Self;
    const FRAC_2_SQRT_PI: Self;

    const PI: Self;
    const FRAC_PI_2: Self;
    const FRAC_PI_4: Self;
    const FRAC_PI_6: Self;
    const FRAC_PI_8: Self;

    const FRAC_1_PI: Self;
}
