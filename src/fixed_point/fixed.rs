use geo_nd::vector::add;

use super::{Int, UsefulConsts, UsefulInt, UsefulUInt};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ArithCode {
    /// Indicates a perfectly valid result was obtained
    Ok,
    /// Indicates result exceeded the range of the most positive value; the actual result provided is ideally the Wrapping value
    ///
    /// If saturating, then set to the most positive value; if checked, return None; if wrapping, then ignore; if plain, then panic
    ///
    /// For values that support inifinites, then the result may be replaced with +Infinity
    OverflowMax,
    /// Indicates result exceeded the range of the most negative value; the actual result provided is ideally the Wrapping value
    ///
    /// If saturating, then set to the most negative value; if checked, return None; if wrapping, then ignore; if plain, then panic
    ///
    /// For values that support inifinites, then the result may be replaced with +Infinity
    OverflowMin,
    /// A division by zero
    ///
    /// If saturating, then set to the most negative value; if checked, return None; if wrapping, then ignore; if plain, then panic
    DivideByZero,
    /// Not a number
    ///
    /// This may be returned for example for a square root of a negative number, or asin of 2, or log of a negative number, etc
    NotANumber,
}

pub trait HowIsFixedPoint<T: UsefulInt> {
    /// Number of bits in T
    const NB: usize = 0;
    /// Number of bits in the fractional part
    const NB_FRAC: usize = 0;
    /// Number of bits in the integer part (including the sign bit)
    ///
    /// This *MUST* conform to NB = NB_INT + NB_FRAC
    const NB_INT: usize = Self::NB - Self::NB_FRAC;
    /// Number of fractional bits in the double, when treated as our associated fixed point value
    ///
    /// The double has twice the number of integer bits and twice the number of fractional bits
    const NB_DBLFRAC: usize = Self::NB_FRAC * 2;
    /// Number of integer bits in the double, when treated as our associated fixed point value
    ///
    /// The double has twice the number of integer bits and twice the number of fractional bits
    const NB_DBLINT: usize = Self::NB_INT * 2;
    /// Bit mask of T of all zeros except the top bit (the sign bit)
    ///
    /// Must equal T::ONE << (Self::NB - 1)
    const SIGN_MASK: T;
    /// Bit mask of T of zeros in the fractional part, ones elsewhwere
    ///
    /// Must equal  (-T::ONE) << Self::NB_FRAC
    const SIGNED_INT_MASK: T;
    /// Bit mask of T of ones in the fractional part, zeros elsewhwere
    ///
    /// Must equal  !Self::SIGNED_INT_MASK
    const FRAC_MASK: T;
    /// The value of one as a fixed point number
    ///
    /// Must equal T::ONE << Self::NB_FRAC;
    const ONE: T;
}

impl HowIsFixedPoint<i8> for FPType<i8, 4> {
    const NB: usize = 8;
    const NB_FRAC: usize = 4;
    const SIGN_MASK: i8 = i8::ONE << (Self::NB - 1);
    const SIGNED_INT_MASK: i8 = (-i8::ONE) << Self::NB_FRAC;
    const FRAC_MASK: i8 = !Self::SIGNED_INT_MASK;
    const ONE: i8 = i8::ONE << Self::NB_FRAC;
}

pub struct FPType<T, const N: usize>([T; N]);

#[repr(transparent)]
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, PartialOrd)]
pub struct Fixed<T: UsefulInt, const N: usize>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    value: T,
}

impl<T: UsefulInt, const N: usize> std::fmt::Display for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <Self as std::fmt::Debug>::fmt(self, f)
    }
}

impl<T: UsefulInt, const N: usize> Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    /// Take a double and shift it down, using the dropped bits to help round appropriately; return True if it does not overflow,
    /// but always set the value
    ///
    /// Round to nearest, with ties broken *away* from zero (same as f32 or f64.round)
    ///
    /// Away from zero means: for +ve, add 1 if fraction >=1/2, for -ve subtract 1 if fraction >=1/2
    ///
    /// Consider an i16_12, that is 0xd800, so a shift of 8 produces an i8_4, then the i8_4 value should be 0xd8, or -2.5
    ///
    /// 0x2800 represents 2.5
    /// 0x287f represents 2.531005859375 - must round to 2.5 (i.e. no add after shift)
    /// 0x2880 represents 2.53125- must round to 2.5625 (i.e. add 1 after shift)
    /// 0x2880 represents 2.531494140625
    /// 0x2900 represents 2.5625
    /// 0xd800 represents -2.5
    /// 0xd880 represents -2.46875 - must round to -2.5 (i.e. no add after shift)
    /// 0xd881 represents -2.468505859375 - must round to -2.4375 (i.e. add 1 after shift)
    /// 0xd900 represents -2.4375
    ///
    /// Hence if the dropped bits are 0x8000... then add 1 after shift
    fn reduce_double(&mut self, dbl: &T::Dbl, shr: usize) -> bool {
        // Shift left by double size minus the shift right gives an (unsigned) bits we will drop
        // If this is (signed) negative then we are dropping >=1/2; if -itself==itself then it is zero or half
        let bits_to_drop = *dbl << (T::NB_DBL - shr);
        let mut add_one_after_shift = !(bits_to_drop & T::DBL_SIGN_MASK).is_zero();
        if *dbl < T::Dbl::ZERO && bits_to_drop == T::DBL_SIGN_MASK {
            add_one_after_shift = false;
        }
        let mut dbl = (*dbl >> shr);
        if add_one_after_shift {
            dbl += T::Dbl::ONE;
        }
        self.value = T::of_dbl_unchecked(dbl);
        let all_set_or_clr = T::DBL_SIGN_MASK | T::DBL_UPPER_MASK;
        dbl &= all_set_or_clr;
        (dbl == all_set_or_clr) || dbl.is_zero()
    }

    #[inline(always)]
    pub(crate) fn do_add(&mut self, other: &Self) -> ArithCode {
        let (r, o) = self.value.carrying_add(other.value, false);
        self.value = r;
        if !o {
            ArithCode::Ok
        } else if self.value < T::ZERO {
            ArithCode::OverflowMin
        } else {
            ArithCode::OverflowMax
        }
    }
    #[inline(always)]
    pub(crate) fn do_sub(&mut self, other: &Self) -> ArithCode {
        let (r, o) = self.value.borrowing_sub(other.value, false);
        self.value = r;
        if !o {
            ArithCode::Ok
        } else if self.value < T::ZERO {
            ArithCode::OverflowMin
        } else {
            ArithCode::OverflowMax
        }
    }
    #[inline(always)]
    pub(crate) fn do_bit_and(&mut self, other: &Self) {
        self.value &= other.value;
    }
    #[inline(always)]
    pub(crate) fn do_bit_or(&mut self, other: &Self) {
        self.value |= other.value;
    }
    #[inline(always)]
    pub(crate) fn do_bit_xor(&mut self, other: &Self) {
        self.value ^= other.value;
    }
    // Should probably change to checked mul returning a bool for overflow (and possibly underflow, where result is 0 but neither input is)
    #[inline(always)]
    pub(crate) fn do_mul(&mut self, other: &Self) {
        // If T = x*2^NB_FRAC, then dbl = x*y*2^(NB_FRAC*2), and the result should be x*y*2^NB_FRAC, so shift right by NB_FRAC
        let shr = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        let dbl = self.value.dbl_mult(&other.value);
        if !self.reduce_double(&dbl, shr) {
            panic!("Overflow in multiply");
        };
    }
    // Should probably change to checked div returning a bool for underflow, overflow, div by zero
    #[inline(always)]
    pub(crate) fn do_div(&mut self, other: &Self) {
        // If T = x*2^NB_FRAC, then s=x*2^(NB_FRAC+NB), o=y*2^NB_FRAC, r=x/y*2^NB, and the result should be x/y*2^NB_FRAC, so shift right by NB-NB_FRAC = NB_INT
        let shr = <FPType<T, N> as HowIsFixedPoint<T>>::NB_INT;
        let s = self.value.as_dbl_upper();
        let o = other.value.as_dbl();
        // r = self/other * 2^(T::NB);
        let r = s / o;
        if !self.reduce_double(&r, shr) {
            panic!("Overflow in divide");
        };
    }
    // x remainder y = x - (x / y).trunc() * y
    //
    // If x=64 (2^14 in i16_8) and y = 1/256 (1 in i16_8) then x/y = 2^14
    //
    // If x,y are each X or Y*2^N, then x/y = is an integer (i.e. already truncated or rounded in some fashion)
    #[inline(always)]
    pub(crate) fn do_rem(&mut self, other: &Self) {
        // Not conviced this works for negative values
        //
        // And possibly it can just be self.value % other.value...
        let x_div_y_trunc = self.value / other.value;
        self.value -= x_div_y_trunc * other.value;
    }
}

use num_traits::{ConstOne, ConstZero, One, Zero};
impl<T: UsefulInt, const N: usize> Zero for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn is_zero(&self) -> bool {
        self == &Self::ZERO
    }
    fn zero() -> Self {
        Self::ZERO
    }
}

impl<T: UsefulInt, const N: usize> ConstZero for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    const ZERO: Self = Self { value: T::ZERO };
}

impl<T: UsefulInt, const N: usize> One for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn is_one(&self) -> bool {
        self == &Self::ONE
    }
    fn one() -> Self {
        Self::ONE
    }
}

impl<T: UsefulInt, const N: usize> ConstOne for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    const ONE: Self = Self {
        value: <FPType<T, N> as HowIsFixedPoint<T>>::ONE,
    };
}

#[test]
fn test_thing() {
    let x = Fixed::<i8, 4>::ONE;
}

macro_rules! binary_op {
    {$trait:ident, $fn:ident, $f:ident} => {
        impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (mut self, other:Fixed<T, N>) -> Fixed<T, N> {
                self.$f(&other); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (mut self, other:&Fixed<T, N>) -> Fixed<T, N> {
                self.$f(other); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for &Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (self, other:Fixed<T, N>) -> Fixed<T, N> {
                let mut s = *self; s.$f(&other); s
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for &Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (self, other:&Fixed<T, N>) -> Fixed<T, N> {
                let mut s = *self; s.$f(&other); s
            }
        }
    }
}

macro_rules! binary_assign_op {
    {$trait:ident, $fn:ident, $f:ident} => {
    impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for Fixed<T, N>
    where
        FPType<T, N>: HowIsFixedPoint<T>,
    {
            fn $fn (&mut self, other:Fixed<T, N>) {
                self.$f(&other);
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            FPType<T, N>: HowIsFixedPoint<T>,
        {
            fn $fn (&mut self, other:&Fixed<T, N>) {
                self.$f(other);
            }
        }
    }
}
use std::ops::*;

impl<T: UsefulInt, const N: usize> Neg for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self { value: -self.value }
    }
}

binary_op! {Add, add, do_add}
binary_assign_op! {AddAssign, add_assign, do_add}
binary_op! {Sub, sub, do_sub}
binary_assign_op! {SubAssign, sub_assign, do_sub}
binary_op! {Mul, mul, do_mul}
binary_assign_op! {MulAssign, mul_assign, do_mul}
binary_op! {Div, div, do_div}
binary_assign_op! {DivAssign, div_assign, do_div}
binary_op! {Rem, rem, do_rem}
binary_assign_op! {RemAssign, rem_assign, do_rem}
binary_op! {BitAnd, bitand, do_bit_and}
binary_assign_op! {BitAndAssign, bitand_assign, do_bit_and}
binary_op! {BitOr, bitor, do_bit_or}
binary_assign_op! {BitOrAssign, bitor_assign, do_bit_or}
binary_op! {BitXor, bitxor, do_bit_xor}
binary_assign_op! {BitXorAssign, bitxor_assign, do_bit_xor}

// num_traits::PrimInt ?
// num_traits::Signed (not Unsigned)
// num_traits::Pow ?
// num_traits::Inv ?
// num_traits::Wrapping*
// num_traits::Saturating*
// num_traits::Checked*

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

impl<T: UsefulInt, const N: usize> num_traits::Signed for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn abs(&self) -> Self {
        if self.value < T::ZERO {
            -*self
        } else {
            *self
        }
    }
    fn signum(&self) -> Self {
        if self.value < T::ZERO {
            Self { value: -T::ONE }
        } else {
            Self { value: T::ONE }
        }
    }

    fn is_positive(&self) -> bool {
        self.value > T::ZERO
    }
    fn is_negative(&self) -> bool {
        self.value < T::ZERO
    }
    fn abs_sub(&self, other: &Self) -> Self {
        if self.value > other.value {
            *self - *other
        } else {
            *other - *self
        }
    }
}

impl<T: UsefulInt, const N: usize> num_traits::Num for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    type FromStrRadixErr = ();
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Ok(Self::default())
    }
}
impl<T: UsefulInt, const N: usize> num_traits::Bounded for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn min_value() -> Self {
        Self {
            value: T::SIGN_MASK,
        }
    }
    fn max_value() -> Self {
        Self {
            value: !T::SIGN_MASK,
        }
    }
}

impl<T: UsefulInt, const N: usize> num_traits::FromPrimitive for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn from_u64(x: u64) -> Option<Self> {
        todo!();
        None
    }
    fn from_i64(n: i64) -> Option<Self> {
        todo!();
        None
    }
}
impl<T: UsefulInt, const N: usize> num_traits::ToPrimitive for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn to_i64(&self) -> Option<i64> {
        todo!();
        None
    }
    fn to_u64(&self) -> Option<u64> {
        todo!();
        None
    }
}

impl<T: UsefulInt, const N: usize> num_traits::NumCast for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn from<V: num::ToPrimitive>(n: V) -> Option<Self> {
        n.to_u64().map(Self::from_u64).flatten()
    }
}

impl<T: UsefulInt + UsefulConsts, const N: usize> num_traits::FloatConst for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn E() -> Self {
        Self { value: T::E }
    }
    fn FRAC_1_PI() -> Self {
        Self {
            value: T::FRAC_1_PI,
        }
    }
    fn FRAC_1_SQRT_2() -> Self {
        Self {
            value: T::FRAC_1_SQRT_2,
        }
    }
    fn FRAC_2_PI() -> Self {
        Self {
            value: T::FRAC_2_PI,
        }
    }
    fn FRAC_2_SQRT_PI() -> Self {
        Self {
            value: T::FRAC_2_SQRT_PI,
        }
    }
    fn FRAC_PI_2() -> Self {
        Self {
            value: T::FRAC_PI_2,
        }
    }
    fn FRAC_PI_3() -> Self {
        Self {
            value: T::FRAC_PI_3,
        }
    }
    fn FRAC_PI_4() -> Self {
        Self {
            value: T::FRAC_PI_4,
        }
    }
    fn FRAC_PI_6() -> Self {
        Self {
            value: T::FRAC_PI_6,
        }
    }
    fn FRAC_PI_8() -> Self {
        Self {
            value: T::FRAC_PI_8,
        }
    }
    fn LN_10() -> Self {
        Self { value: T::LN_10 }
    }
    fn LN_2() -> Self {
        Self { value: T::LN_2 }
    }
    fn LOG10_E() -> Self {
        Self { value: T::LOG10_E }
    }
    fn LOG2_E() -> Self {
        Self { value: T::LOG2_E }
    }
    fn PI() -> Self {
        Self { value: T::PI }
    }
    fn SQRT_2() -> Self {
        Self { value: T::SQRT_2 }
    }

    fn TAU() -> Self {
        Self { value: T::TAU }
    }

    fn LOG10_2() -> Self {
        Self { value: T::LOG10_2 }
    }

    fn LOG2_10() -> Self {
        Self { value: T::LOG2_10 }
    }
}

impl<T: UsefulInt, const N: usize> num_traits::Float for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    fn nan() -> Self {
        panic!("Fixed point numbers do not have a NaN");
    }
    fn infinity() -> Self {
        panic!("Fixed point numbers do not have infinity");
    }
    fn neg_infinity() -> Self {
        panic!("Fixed point numbers do not have infinity");
    }
    fn neg_zero() -> Self {
        Self::ZERO
    }
    fn min_value() -> Self {
        Self { value: T::ONE }
    }
    fn min_positive_value() -> Self {
        Self { value: T::ONE }
    }
    fn max_value() -> Self {
        let sign_bit = <FPType<T, N> as HowIsFixedPoint<T>>::SIGN_MASK;
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        Self {
            value: int_mask & !sign_bit,
        }
    }
    fn is_nan(self) -> bool {
        false
    }
    fn is_infinite(self) -> bool {
        false
    }
    fn is_finite(self) -> bool {
        true
    }
    fn is_normal(self) -> bool {
        !self.is_zero()
    }
    fn classify(self) -> std::num::FpCategory {
        if self.is_zero() {
            std::num::FpCategory::Zero
        } else {
            std::num::FpCategory::Normal
        }
    }
    fn floor(self) -> Self {
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        Self {
            value: self.value & int_mask,
        }
    }
    fn ceil(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        let int_one = <FPType<T, N> as HowIsFixedPoint<T>>::ONE;
        if self.value & frac_mask == T::ZERO {
            self
        } else {
            Self {
                // This should panic on overflow
                value: (self.value & int_mask) + int_one,
            }
        }
    }
    fn round(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        let half_frac_mask = frac_mask >> 1_u8;
        if self.value > T::ZERO {
            Self {
                value: (self.value + (half_frac_mask + T::ONE)) & frac_mask,
            }
        } else {
            Self {
                value: (self.value + half_frac_mask) & frac_mask,
            }
        }
    }
    fn trunc(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        let int_mask = <FPType<T, N> as HowIsFixedPoint<T>>::SIGNED_INT_MASK;
        // remove the fractional part; for +ve numbers this is clear the bits, for -ve numbers add the fractional mask the clear the bits
        if self.value < T::ZERO {
            Self {
                value: (self.value + frac_mask) & int_mask,
            }
        } else {
            Self {
                value: self.value & int_mask,
            }
        }
    }
    fn fract(self) -> Self {
        let frac_mask = <FPType<T, N> as HowIsFixedPoint<T>>::FRAC_MASK;
        Self {
            value: self.value & frac_mask,
        }
    }
    fn abs(self) -> Self {
        if self.value < T::ZERO {
            Self { value: -self.value }
        } else {
            Self { value: -self.value }
        }
    }
    fn signum(self) -> Self {
        if self.value < T::ZERO {
            -Self::ONE
        } else {
            Self::ONE
        }
    }
    fn is_sign_positive(self) -> bool {
        self.value >= T::ZERO
    }
    fn is_sign_negative(self) -> bool {
        self.value < T::ZERO
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        todo!();
    }
    fn recip(self) -> Self {
        todo!();
    }
    fn powi(self, n: i32) -> Self {
        todo!();
    }
    fn powf(self, n: Self) -> Self {
        todo!();
    }
    fn sqrt(self) -> Self {
        todo!();
    }
    fn exp(self) -> Self {
        todo!();
    }
    fn exp2(self) -> Self {
        todo!();
    }
    fn ln(self) -> Self {
        todo!();
    }
    fn log(self, base: Self) -> Self {
        todo!();
    }
    fn log2(self) -> Self {
        todo!();
    }
    fn log10(self) -> Self {
        todo!();
    }
    fn max(self, other: Self) -> Self {
        if self.value > other.value {
            self
        } else {
            other
        }
    }
    fn min(self, other: Self) -> Self {
        if self.value < other.value {
            self
        } else {
            other
        }
    }
    fn abs_sub(self, other: Self) -> Self {
        if self.value > other.value {
            self - other
        } else {
            other - self
        }
    }
    fn cbrt(self) -> Self {
        todo!();
    }
    fn hypot(self, other: Self) -> Self {
        todo!();
    }
    fn sin(self) -> Self {
        todo!();
    }
    fn cos(self) -> Self {
        todo!();
    }
    fn tan(self) -> Self {
        todo!();
    }
    fn asin(self) -> Self {
        todo!();
    }
    fn acos(self) -> Self {
        todo!();
    }
    fn atan(self) -> Self {
        todo!();
    }
    fn atan2(self, other: Self) -> Self {
        todo!();
    }
    fn sin_cos(self) -> (Self, Self) {
        todo!();
    }
    fn exp_m1(self) -> Self {
        todo!();
    }
    fn ln_1p(self) -> Self {
        todo!();
    }
    fn sinh(self) -> Self {
        todo!();
    }
    fn cosh(self) -> Self {
        todo!();
    }
    fn tanh(self) -> Self {
        todo!();
    }
    fn asinh(self) -> Self {
        todo!();
    }
    fn acosh(self) -> Self {
        todo!();
    }
    fn atanh(self) -> Self {
        todo!();
    }
    fn integer_decode(self) -> (u64, i16, i8) {
        todo!();
    }

    fn epsilon() -> Self {
        Self { value: T::ONE }
    }
    fn is_subnormal(self) -> bool {
        false
    }
    fn to_degrees(self) -> Self {
        todo!();
    }
    fn to_radians(self) -> Self {
        todo!();
    }
    fn clamp(self, min: Self, max: Self) -> Self {
        if self < min {
            min
        } else if self > max {
            max
        } else {
            self
        }
    }
    fn copysign(self, sign: Self) -> Self {
        if (self.value < T::ZERO) == (sign.value < T::ZERO) {
            self
        } else {
            Self { value: -self.value }
        }
    }
}

use num_traits::FromPrimitive;

impl<T: UsefulInt, const N: usize> crate::NumOps for Fixed<T, N>
where
    FPType<T, N>: HowIsFixedPoint<T>,
{
    /// Create a value from a fraction with signed numerator and unsigned denominator
    fn frac(n: i32, u: u32) -> Self {
        Self::from_i32(n).unwrap() / Self::from_u32(u).unwrap()
    }

    /// Create a value from an [i32]
    fn of_i32(n: i32) -> Self {
        Self::from_i32(n).unwrap()
    }

    /// Convert an f64 to this value
    fn of_f64(v: f64) -> Self {
        Self::default()
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn reciprocal(self) -> Self {
        Self::ONE / self
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn is_unreliable_divisor(self) -> bool {
        false
    }

    fn powi(self, p: i32) -> Self {
        (0..p).fold(Self::ONE, |acc, _| acc * self)
    }

    /// Return an estimate of the square root to within a precision
    fn sqrt_est(self) -> Self {
        // sqrt with a u32_24 input produces a u32_28 output
        // ; a u64_0 (pure integer) produces a u64_32 output
        let nb = <FPType<T, N> as HowIsFixedPoint<T>>::NB_FRAC;
        if (nb & 1) == 0 {
            let x_uns_half_frac = super::functions::sqrt(self.value.unsigned());
            Self {
                value: T::of_unsigned(x_uns_half_frac) >> (nb / 2),
            }
        } else {
            let x_uns_half_frac = super::functions::sqrt(self.value.unsigned() >> 1_u8);
            Self {
                value: T::of_unsigned(x_uns_half_frac) >> (nb / 2),
            }
        }
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn cbrt_est(self) -> Self {
        crate::utils::cbrt_est::<_, 5>(self)
    }
    /// Return true if the value is negative
    fn is_sign_negative(self) -> bool {
        self < Self::ZERO
    }
}
