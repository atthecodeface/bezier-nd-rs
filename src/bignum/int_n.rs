use super::UIntN;
use num_traits::Zero;
use num_traits::{ops::overflowing::OverflowingMul, FromPrimitive};
use std::ops::*;

/// A signed integer, of +- 64*N bits (stored with sign and magnitude),
/// supporting copy
///
/// This is not optimized for performance, but to enable use of algorithms
/// that require num_traits::Num
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IntN<const N: usize> {
    /// Sign of the number
    pub(crate) is_neg: bool,
    /// Magnitude of the number
    pub(crate) value: UIntN<N>,
}

impl<const N: usize> std::convert::From<u64> for IntN<N> {
    fn from(value: u64) -> Self {
        Self {
            is_neg: false,
            value: value.into(),
        }
    }
}

impl<const N: usize> std::convert::From<i64> for IntN<N> {
    fn from(value: i64) -> Self {
        if value < 0 {
            Self {
                is_neg: true,
                value: ((-value) as u64).into(),
            }
        } else {
            (value as u64).into()
        }
    }
}

impl<const N: usize> std::convert::From<UIntN<N>> for IntN<N> {
    fn from(value: UIntN<N>) -> Self {
        Self {
            is_neg: false,
            value,
        }
    }
}

impl<const N: usize> std::convert::From<(bool, UIntN<N>)> for IntN<N> {
    fn from((is_neg, value): (bool, UIntN<N>)) -> Self {
        Self { is_neg, value }
    }
}

impl<const N: usize> std::cmp::Ord for IntN<N> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self.is_neg, other.is_neg) {
            (true, true) => other.value.cmp(&self.value),
            (false, false) => self.value.cmp(&other.value),
            (false, true) => std::cmp::Ordering::Greater,
            (true, false) => std::cmp::Ordering::Less,
        }
    }
}

impl<const N: usize> std::cmp::PartialOrd for IntN<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> std::fmt::Display for IntN<N> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let s: String = self.value.as_digits(10).collect();
        fmt.pad_integral(!self.is_neg, "", &s)
    }
}

impl<const N: usize> std::fmt::Binary for IntN<N> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let value_string = format!("{:b}", self.value);
        fmt.pad_integral(!self.is_neg, "0b", &value_string)
    }
}

impl<const N: usize> std::fmt::LowerHex for IntN<N> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let value_string = format!("{:x}", self.value);
        fmt.pad_integral(!self.is_neg, "0x", &value_string)
    }
}

impl<const N: usize> std::fmt::UpperHex for IntN<N> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let value_string = format!("{:X}", self.value);
        fmt.pad_integral(!self.is_neg, "0x", &value_string)
    }
}

impl<const N: usize> std::ops::Not for IntN<N> {
    type Output = Self;
    fn not(mut self) -> Self {
        self.value = !self.value;
        self
    }
}

impl<const N: usize> std::ops::Neg for IntN<N> {
    type Output = Self;
    fn neg(mut self) -> Self {
        if !self.value_is_zero() {
            self.is_neg = !self.is_neg;
        }
        self
    }
}

macro_rules! binary_op {
    {$trait:ident, $fn:ident, $f:ident} => {
        impl<const N: usize> $trait<IntN<N>> for IntN<N> {

            type Output = IntN<N>;
            fn $fn (mut self, other:IntN<N>) -> IntN<N> {
                self.$f(&other); self
            }
        }
        impl<const N: usize> $trait<&IntN<N>> for IntN<N> {

            type Output = IntN<N>;
            fn $fn (mut self, other:&IntN<N>) -> IntN<N> {
                self.$f(other); self
            }
        }
        impl<const N: usize> $trait<IntN<N>> for &IntN<N> {

            type Output = IntN<N>;
            fn $fn (self, other:IntN<N>) -> IntN<N> {
                let mut s = *self; s.$f(&other); s
            }
        }
        impl<const N: usize> $trait<&IntN<N>> for &IntN<N> {

            type Output = IntN<N>;
            fn $fn (self, other:&IntN<N>) -> IntN<N> {
                let mut s = *self; s.$f(&other); s
            }
        }
    }
}

macro_rules! binary_assign_op {
    {$trait:ident, $fn:ident, $f:ident} => {
        impl<const N: usize> $trait<IntN<N>> for IntN<N> {

            fn $fn (&mut self, other:IntN< N>) {
                self.$f(&other);
            }
        }
        impl<const N: usize> $trait<&IntN<N>> for IntN<N> {
            fn $fn (&mut self, other:&IntN<N>) {
                self.$f(other);
            }
        }
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

macro_rules! extra_binary_ops {
    {overflowing, $trait:ident, $fn:ident, $f:ident} => {
        impl<const N: usize> num_traits::ops::overflowing::$trait for IntN<N> {
            fn $fn(&self, v: &Self) -> (Self, bool) {
                let mut s = *self;
                let overflow = s.$f(v);
                (s, overflow)
            }
        }
    };
{wrapping, $trait:ident, $fn:ident, $f:ident} => {
impl<const N: usize> num_traits::ops::wrapping::$trait for IntN<N> {
    fn $fn(&self, v: &Self) -> Self {
        let mut s = *self;
        s.$f(v);
        s
    }
}
};
{checked, $trait:ident, $fn:ident, $f:ident} => {
impl<const N: usize> num_traits::ops::checked::$trait for IntN<N> {
    fn $fn(&self, v: &Self) -> Option<Self> {
        let mut s = *self;
        s.$f(v).then_some(s)
    }
}
}
}
extra_binary_ops! {overflowing, OverflowingAdd, overflowing_add, do_add_overflow}
extra_binary_ops! {wrapping, WrappingAdd, wrapping_add, do_add_overflow}
extra_binary_ops! {checked, CheckedAdd, checked_add, do_add_overflow}
extra_binary_ops! {overflowing, OverflowingSub, overflowing_sub, do_sub_overflow}
extra_binary_ops! {wrapping, WrappingSub, wrapping_sub, do_sub_overflow}
extra_binary_ops! {checked, CheckedSub, checked_sub, do_sub_overflow}
extra_binary_ops! {overflowing, OverflowingMul, overflowing_mul, do_mul_overflow}
extra_binary_ops! {wrapping, WrappingMul, wrapping_mul, do_mul_overflow}
extra_binary_ops! {checked, CheckedMul, checked_mul, do_mul_overflow}

macro_rules! shf_op {
    {$trait:ident, $t:ty, $fn:ident, $f:ident} => {
        impl<const N: usize> $trait<$t> for IntN<N> {
            type Output = IntN<N>;
            fn $fn(mut self, rhs: $t) ->IntN<N> {
                let rhs = rhs as u32;
                self. $f (rhs);
                self
            }
        }
        impl<const N: usize> $trait<&$t> for IntN<N> {
            type Output = IntN<N>;
            fn $fn(mut self, rhs: &$t) ->IntN<N> {
                let rhs = *rhs as u32;
                self. $f (rhs);
                self
            }
        }
        impl<const N: usize> $trait<$t> for &IntN<N> {
            type Output = IntN<N>;
            fn $fn(self, rhs: $t) -> IntN<N> {
                let mut r = *self;
                let rhs = rhs as u32;
                r. $f (rhs);
                r
            }
        }
        impl<const N: usize> $trait<&$t> for &IntN<N> {
            type Output = IntN<N>;
            fn $fn(self, rhs: &$t) -> IntN<N> {
                let mut r = *self;
                let rhs = *rhs as u32;
                r. $f (rhs);
                r
            }
        }
    }
}
macro_rules! shf_assign_op {
    {$trait:ident, $t:ty, $fn:ident, $f:ident} => {
        impl<const N: usize> $trait<$t> for IntN<N> {
            fn $fn(&mut self, rhs: $t) {
                let rhs = rhs as u32;
                self. $f (rhs);
            }
        }
        impl<const N: usize> $trait<&$t> for IntN<N> {
            fn $fn(&mut self, rhs: &$t) {
                let rhs = *rhs as u32;
                self. $f (rhs);
            }
        }
    }
}

shf_op! {Shr, u8, shr, shift_right}
shf_op! {Shr, u16, shr, shift_right}
shf_op! {Shr, u32, shr, shift_right}
shf_op! {Shr, u64, shr, shift_right}
shf_op! {Shr, usize, shr, shift_right}
shf_op! {Shl, u8, shl, shift_left}
shf_op! {Shl, u16, shl, shift_left}
shf_op! {Shl, u32, shl, shift_left}
shf_op! {Shl, u64, shl, shift_left}
shf_op! {Shl, usize, shl, shift_left}

shf_assign_op! {ShrAssign, u8, shr_assign, shift_right}
shf_assign_op! {ShrAssign, u16, shr_assign, shift_right}
shf_assign_op! {ShrAssign, u32, shr_assign, shift_right}
shf_assign_op! {ShrAssign, u64, shr_assign, shift_right}
shf_assign_op! {ShrAssign, usize, shr_assign, shift_right}
shf_assign_op! {ShlAssign, u8, shl_assign, shift_left}
shf_assign_op! {ShlAssign, u16, shl_assign, shift_left}
shf_assign_op! {ShlAssign, u32, shl_assign, shift_left}
shf_assign_op! {ShlAssign, u64, shl_assign, shift_left}
shf_assign_op! {ShlAssign, usize, shl_assign, shift_left}

impl<const N: usize> std::ops::Mul<&UIntN<N>> for IntN<N> {
    type Output = Self;

    fn mul(mut self, other: &UIntN<N>) -> Self {
        self.value *= *other;
        self
    }
}
impl<const N: usize> std::ops::Mul<UIntN<N>> for IntN<N> {
    type Output = Self;

    fn mul(mut self, other: UIntN<N>) -> Self {
        self.value *= other;
        self
    }
}

impl<const N: usize> std::ops::Div<UIntN<N>> for IntN<N> {
    type Output = Self;

    #[track_caller]
    fn div(mut self, other: UIntN<N>) -> Self {
        self.value /= other;
        self
    }
}

impl<const N: usize> IntN<N> {
    /// Create a new value
    pub const fn new(is_neg: bool, value: &[u64; N]) -> Self {
        Self {
            is_neg,
            value: UIntN::new(value),
        }
    }

    /// Return the value with just the specified bit set
    pub const fn with_bit_set(bit: u32) -> Self {
        Self {
            is_neg: false,
            value: UIntN::with_bit_set(bit),
        }
    }

    /// Get the magnitude of the integer
    pub fn magnitude(&self) -> &UIntN<N> {
        &self.value
    }

    /// Get the sign of the integer
    pub fn is_neg(&self) -> bool {
        self.is_neg
    }

    fn value_is_zero(&self) -> bool {
        use num_traits::Zero;
        self.value.is_zero()
    }

    /// 20 / 7 = 2 remainder 6
    /// 20 / -7 = -2 remainder 6
    /// -20 / 7 = -2 remainder -6
    /// -20 / -7 = 2 remainder -6
    fn do_div_rem(&self, other: &Self) -> Option<(Self, Self)> {
        if other.value_is_zero() {
            None
        } else {
            let Some((div, rem)) = self.value.do_div_rem(&other.value) else {
                panic!("Divide by zero");
            };
            let div = Self {
                is_neg: self.is_neg != other.is_neg,
                value: div,
            };
            let rem = Self {
                is_neg: self.is_neg,
                value: rem,
            };
            Some((div, rem))
        }
    }

    #[track_caller]
    fn do_bit_and(&mut self, other: &Self) {
        self.value &= other.value;
        self.is_neg = self.is_neg && other.is_neg;
    }

    #[track_caller]
    fn do_bit_or(&mut self, other: &Self) {
        self.value |= other.value;
        self.is_neg = self.is_neg || other.is_neg;
    }

    #[track_caller]
    fn do_bit_xor(&mut self, other: &Self) {
        self.value ^= other.value;
        self.is_neg = self.is_neg != other.is_neg;
    }

    #[track_caller]
    fn do_add_overflow(&mut self, other: &Self) -> bool {
        self.do_add_sub(&other, false)
    }
    #[track_caller]
    fn do_sub_overflow(&mut self, other: &Self) -> bool {
        self.do_add_sub(&other, true)
    }
    #[track_caller]
    fn do_mul_overflow(&mut self, other: &Self) -> bool {
        self.is_neg = self.is_neg != other.is_neg;
        let (a, b) = self.value.overflowing_mul(&other.value);
        self.value = a;
        b
    }

    #[track_caller]
    fn do_add(&mut self, other: &Self) {
        assert!(!self.do_add_overflow(other), "Addition overflowed");
    }

    #[track_caller]
    fn do_sub(&mut self, other: &Self) {
        assert!(!self.do_sub_overflow(other), "Subtratcion underflowed");
    }

    #[track_caller]
    pub fn shift_left(&mut self, amount: u32) -> bool {
        self.value.shift_left(amount)
    }
    #[track_caller]
    pub fn shift_right(&mut self, amount: u32) {
        self.value.shift_right(amount)
    }

    #[track_caller]
    fn do_mul(&mut self, other: &Self) {
        self.is_neg = self.is_neg != other.is_neg;
        self.value *= other.value;
    }

    #[track_caller]
    fn do_div(&mut self, other: &Self) {
        let Some((d, _r)) = self.do_div_rem(other) else {
            panic!("Division by zero");
        };
        *self = d;
    }

    #[track_caller]
    fn do_rem(&mut self, other: &Self) {
        let Some((_d, r)) = self.do_div_rem(other) else {
            panic!("Division by zero");
        };
        *self = r;
    }

    #[track_caller]
    fn do_add_sub(&mut self, other: &Self, negate_other: bool) -> bool {
        let mut other_neg = other.is_neg;
        if negate_other {
            other_neg = !other_neg;
        }
        if self.is_neg != other_neg {
            // if self is +1 and other is -3 then subtract yields a borrow
            // and self becomes [u64::MAX, .. u64::MAX, , u64::MAX-1]
            //
            // we need to replace with U64::MAX-v for each v (effectively +1 is
            // we really store things as ones complement)
            use num_traits::ops::overflowing::OverflowingSub;
            let (value, borrow) = self.value.overflowing_sub(&other.value);
            self.value = value;
            if borrow {
                self.is_neg = !self.is_neg;
                self.value.twos_complement();
            }
            if self.value.is_zero() {
                self.is_neg = false;
            }
            false
        } else {
            use num_traits::ops::overflowing::OverflowingAdd;
            let (value, carry) = self.value.overflowing_add(&other.value);
            self.value = value;
            if self.value.is_zero() {
                self.is_neg = false;
            }
            carry
        }
    }
}

impl<const N: usize> num_traits::ToPrimitive for IntN<N> {
    fn to_i64(&self) -> Option<i64> {
        self.value
            .to_i64()
            .map(|s| if self.is_neg { -s } else { s })
    }
    fn to_u64(&self) -> Option<u64> {
        if self.is_neg {
            None
        } else {
            self.value.to_u64()
        }
    }
}

impl<const N: usize> num_traits::FromPrimitive for IntN<N> {
    fn from_i64(n: i64) -> Option<Self> {
        Some(n.into())
    }
    fn from_u64(n: u64) -> Option<Self> {
        Some(n.into())
    }
}

impl<const N: usize> num_traits::NumCast for IntN<N> {
    fn from<V: num_traits::ToPrimitive>(n: V) -> Option<Self> {
        n.to_u64().map(Self::from_u64).flatten()
    }
}

impl<const N: usize> num_traits::Num for IntN<N> {
    type FromStrRadixErr = std::num::ParseIntError;
    fn from_str_radix(src: &str, radix: u32) -> Result<Self, std::num::ParseIntError> {
        use num_traits::Zero;
        let mut is_neg = src.starts_with('-');
        let mut c = src.chars();
        if is_neg {
            let _ = c.next();
        }
        let value = UIntN::from_str_radix(c.as_str(), radix)?;
        if value.is_zero() {
            is_neg = false;
        }
        Ok(Self { is_neg, value })
    }
}

impl<const N: usize> num_traits::identities::One for IntN<N> {
    fn one() -> Self {
        Self {
            value: UIntN::one(),
            ..Default::default()
        }
    }
}

impl<const N: usize> num_traits::identities::Zero for IntN<N> {
    fn zero() -> Self {
        Self::default()
    }
    fn is_zero(&self) -> bool {
        self.value_is_zero()
    }
}

impl<const N: usize> num_traits::identities::ConstZero for IntN<N> {
    const ZERO: Self = Self {
        is_neg: false,
        value: UIntN::ZERO,
    };
}

impl<const N: usize> num_traits::identities::ConstOne for IntN<N> {
    const ONE: Self = Self {
        is_neg: false,
        value: <UIntN<_> as num_traits::identities::ConstOne>::ONE,
    };
}
