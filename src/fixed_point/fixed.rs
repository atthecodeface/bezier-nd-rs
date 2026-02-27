use super::{Int, UsefulInt, UsefulUInt};

pub trait IsFixedPointThing<T: Int> {
    const ONE: T = T::ONE;
    const FRAC: usize = 0;
    const DBLFRAC: usize = 0;
}

impl IsFixedPointThing<i8> for Blob<i8, 4> {
    const ONE: i8 = 1 << 4;
    const FRAC: usize = 4;
    const DBLFRAC: usize = 8;
}

pub struct Blob<T, const N: usize>([T; N]);

#[repr(transparent)]
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, PartialOrd)]
pub struct Fixed<T: UsefulInt, const N: usize>
where
    Blob<T, N>: IsFixedPointThing<T>,
{
    value: T,
}

impl<T: UsefulInt, const N: usize> std::fmt::Display for Fixed<T, N>
where
    Blob<T, N>: IsFixedPointThing<T>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <Self as std::fmt::Debug>::fmt(self, f)
    }
}

impl<T: UsefulInt, const N: usize> Fixed<T, N>
where
    Blob<T, N>: IsFixedPointThing<T>,
{
    #[inline(always)]
    pub(crate) fn do_add(&mut self, other: &Self) {
        self.value += other.value;
    }
    #[inline(always)]
    pub(crate) fn do_sub(&mut self, other: &Self) {
        self.value -= other.value;
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
    #[inline(always)]
    pub(crate) fn do_mul(&mut self, other: &Self) {
        let Some(result) = T::of_dbl(
            self.value.dbl_mult(&other.value) >> <Blob<T, N> as IsFixedPointThing<T>>::DBLFRAC,
        ) else {
            panic!("Overflow in multiply");
        };
        self.value = result;
    }
    #[inline(always)]
    pub(crate) fn do_div(&mut self, other: &Self) {
        // let s: SignedRaw3232 = (*self).into();
        // let o: SignedRaw3232 = (*other).into();
        // *self = (s / o).try_into().unwrap();
    }
    #[inline(always)]
    pub(crate) fn do_rem(&mut self, other: &Self) {
        // let s: SignedRaw3232 = (*self).into();
        // let o: SignedRaw3232 = (*other).into();
        // *self = (s % o).try_into().unwrap();
    }
}

use num_traits::{ConstOne, ConstZero, One, Zero};
impl<T: UsefulInt, const N: usize> Zero for Fixed<T, N>
where
    Blob<T, N>: IsFixedPointThing<T>,
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
    Blob<T, N>: IsFixedPointThing<T>,
{
    const ZERO: Self = Self { value: T::ZERO };
}

impl<T: UsefulInt, const N: usize> One for Fixed<T, N>
where
    Blob<T, N>: IsFixedPointThing<T>,
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
    Blob<T, N>: IsFixedPointThing<T>,
{
    const ONE: Self = Self {
        value: <Blob<T, N> as IsFixedPointThing<T>>::ONE,
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
            Blob<T, N>: IsFixedPointThing<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (mut self, other:Fixed<T, N>) -> Fixed<T, N> {
                self.$f(&other); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            Blob<T, N>: IsFixedPointThing<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (mut self, other:&Fixed<T, N>) -> Fixed<T, N> {
                self.$f(other); self
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<Fixed<T,N>> for &Fixed<T, N>
        where
            Blob<T, N>: IsFixedPointThing<T>,
        {
            type Output = Fixed<T, N>;
            fn $fn (self, other:Fixed<T, N>) -> Fixed<T, N> {
                let mut s = *self; s.$f(&other); s
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for &Fixed<T, N>
        where
            Blob<T, N>: IsFixedPointThing<T>,
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
        Blob<T, N>: IsFixedPointThing<T>,
    {
            fn $fn (&mut self, other:Fixed<T, N>) {
                self.$f(&other);
            }
        }
        impl<T: UsefulInt, const N: usize> $trait<&Fixed<T,N>> for Fixed<T, N>
        where
            Blob<T, N>: IsFixedPointThing<T>,
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
    Blob<T, N>: IsFixedPointThing<T>,
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

impl<T: UsefulInt, const N: usize> num_traits::Num for Fixed<T, N>
where
    Blob<T, N>: IsFixedPointThing<T>,
{
    type FromStrRadixErr = ();
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Ok(Self::default())
    }
}

impl<T: UsefulInt, const N: usize> num_traits::FromPrimitive for Fixed<T, N>
where
    Blob<T, N>: IsFixedPointThing<T>,
{
    fn from_u64(x: u64) -> Option<Self> {
        None
    }
    fn from_i64(n: i64) -> Option<Self> {
        None
    }
}

use num_traits::FromPrimitive;

impl<T: UsefulInt, const N: usize> crate::NumOps for Fixed<T, N>
where
    Blob<T, N>: IsFixedPointThing<T>,
{
    /// Create a value from a fraction with signed numerator and unsigned denominator
    fn frac(n: i32, u: u32) -> Self {
        Self::default()
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
        crate::utils::sqrt_est::<_, 5>(self, true)
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
