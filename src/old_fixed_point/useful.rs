use num_traits::{ConstOne, ConstZero};

/// A trait that bundles up traits required for UsefulInt and UsefulUInt
///
/// We want UsefulInt and UsefulUInt to probably be visible externally
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

pub trait SignedInt:
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

/// The trait supplying the methods and associated types required for a value to be used in fixed point calculations
pub trait UsefulInt: SignedInt {
    /// The associated unsigned value type
    type Unsigned: UsefulUInt;
    /// The associated signed value type with double the number of bits
    type Dbl: Int;
    /// The associated unsigned value type with double the number of bits
    type DblUnsigned: Int;
    /// The number of bits (including a sign bit if negative numbers are twos complement)
    const NB: usize;
    /// The number of bits in the double type
    const NB_DBL: usize = Self::NB * 2;
    /// Return this sign extended into a value with double the number of bits
    fn as_dbl(self) -> Self::Dbl;
    /// Return this as a value with double the number of bits, with the least significant bits set to zero
    fn as_dbl_upper(self) -> Self::Dbl;
    /// Return this zero-extended into a value with double the number of bits
    fn as_dbl_unsigned(self) -> Self::DblUnsigned;
    /// Reutrn this as an unsigned value using the associated type
    ///
    /// Does this panic if the value is negative?
    fn unsigned(self) -> Self::Unsigned;
    /// Return a value given an unsigned value
    ///
    /// Does this panic if the value has the top bit set and this is twos complement?
    fn of_unsigned(v: Self::Unsigned) -> Self;
    /// Return the value of the least significant bits of the double as this
    /// type, and indicate whether the value overflowed (i.e. the value returned
    /// does not match the value in the original double)
    fn of_dbl(dbl: Self::Dbl) -> (Self, bool);
    /// Multiply by another value producing a value that is double the number of bits; this loses no precision
    fn dbl_mult(self, other: &Self) -> Self::Dbl;
    /// The minimum value representable; for two's complement this is the value with only the top bit set
    fn min_value() -> Self;
    /// The maximum value representable; for two's complement this is the value with only the top bit clear
    fn max_value() -> Self;
}
/// The trait supplying the methods and associated types required for a value to be used in fixed point calculations
pub trait UsefulIntTrig: UsefulInt {
    /// Given self in the range -1/2 to +1/2 reflecting an angle of -PI/4 to +PI/4, calculate sin and cos
    fn sincos_first_quad(self) -> (Self, Self);
    fn atan2(self, x: Self) -> Self {
        Self::ZERO
    }
    fn asin(self) -> Self {
        Self::ZERO
    }
    fn acos(self) -> Self {
        Self::ZERO
    }
}

/// The trait that must be implemented for the Fixed value to be able to use stuff
///
/// Can we lose this one if we have the methods in UsefulInt?
pub trait UsefulUInt: Int {
    /// A type with double the number of bits
    type Dbl: Int;
    /// The number of bits
    const NB: usize;
    /// Convert Self to a double version with the most significant bits all zeros
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

impl UsefulIntTrig for i32 {
    fn sincos_first_quad(self) -> (Self, Self) {
        use super::functions::i32_28::{
            ATAN_ANGLES_I32_28, COS_SCALE_U32_28, HALF_PI_U32_28, NEG_POW2_I32_28,
        };
        let a = self as i64 * (HALF_PI_U32_28 as i64);
        let (c, s) = super::functions::apply_rotation_table::<_, 1>(
            (a >> 32) as i32,
            (COS_SCALE_U32_28 as i32, 0),
            ATAN_ANGLES_I32_28,
            NEG_POW2_I32_28,
        )
        .1;
        (s, c)
    }
}
