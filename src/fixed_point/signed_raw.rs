use super::UnsignedRaw3232;
use num_traits::{
    ConstOne, ConstZero, FromPrimitive, Num, NumCast, One, Signed, ToPrimitive, Zero,
};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct SignedRaw3232 {
    is_neg: bool,
    raw: UnsignedRaw3232,
}

impl std::convert::From<u32> for SignedRaw3232 {
    fn from(i: u32) -> Self {
        Self {
            is_neg: false,
            raw: i.into(),
        }
    }
}

impl std::convert::From<i32> for SignedRaw3232 {
    fn from(i: i32) -> Self {
        Self {
            is_neg: i < 0,
            raw: ((i & i32::MAX) as u32).into(),
        }
    }
}

impl std::convert::From<(i32, u32)> for SignedRaw3232 {
    fn from((i, f): (i32, u32)) -> Self {
        Self {
            is_neg: i < 0,
            raw: (((i & i32::MAX) as u32), f).into(),
        }
    }
}

impl std::convert::TryFrom<i64> for SignedRaw3232 {
    type Error = ();
    #[inline]
    fn try_from(v: i64) -> Result<Self, ()> {
        if v >= (u32::MIN as i64) && v < ((u32::MAX) as i64) {
            Ok(Self {
                is_neg: v < 0,
                raw: (v as u32).into(),
            })
        } else {
            Err(())
        }
    }
}

impl std::convert::TryFrom<u64> for SignedRaw3232 {
    type Error = ();
    #[inline]
    fn try_from(v: u64) -> Result<Self, ()> {
        if v < 1 << 32 {
            Ok((v as u32).into())
        } else {
            Err(())
        }
    }
}

impl std::convert::TryFrom<f32> for SignedRaw3232 {
    type Error = ();
    #[inline]
    fn try_from(f: f32) -> Result<Self, ()> {
        Self::of_f32_bits(f.to_bits()).ok_or(())
    }
}

impl std::convert::TryFrom<f64> for SignedRaw3232 {
    type Error = ();
    #[inline]
    fn try_from(f: f64) -> Result<Self, ()> {
        Self::of_f64_bits(f.to_bits()).ok_or(())
    }
}

impl std::convert::From<SignedRaw3232> for f64 {
    #[inline]
    fn from(r: SignedRaw3232) -> f64 {
        r.as_f64_bits().map(f64::from_bits).unwrap_or(f64::NAN)
    }
}

impl std::convert::From<SignedRaw3232> for f32 {
    #[inline]
    fn from(r: SignedRaw3232) -> f32 {
        r.as_f32_bits().map(f32::from_bits).unwrap_or(f32::NAN)
    }
}

impl std::cmp::Ord for SignedRaw3232 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self.is_neg, other.is_neg) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            _ => self.raw.cmp(&other.raw),
        }
    }
}

impl std::cmp::PartialOrd for SignedRaw3232 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl SignedRaw3232 {
    pub(crate) fn is_neg(&self) -> bool {
        self.is_neg
    }
    pub(crate) fn int(&self) -> u32 {
        self.raw.int()
    }
    pub(crate) fn frac(&self) -> u32 {
        self.raw.frac()
    }
    /// Take a mantissa (bit 63 set if not zero) and exponent (0=> divide by 2^63)
    /// to generate a (normalized) Raw number
    pub(crate) fn of_sign_mantissa_exp(is_neg: bool, mantissa: u64, exp: i32) -> Option<Self> {
        UnsignedRaw3232::of_mantissa_exp(mantissa, exp).map(|raw| Self {
            is_neg: is_neg && !raw.is_zero(),
            raw,
        })
    }

    /// Get a [UnsignedRaw3232] from f64 bits, if it fits in the range
    /// for the numerator/denominator, by using powers of two
    /// for the denominator
    pub fn of_f64_bits(bits: u64) -> Option<Self> {
        let is_neg = (bits & (1 << 63)) != 0;
        UnsignedRaw3232::of_f64_bits(bits & !(1 << 63)).map(|raw| Self { is_neg, raw })
    }

    /// Get an [UnsignedRaw3232] from f32 bits, if it fits in the range
    /// for the numerator/denominator, by using powers of two
    /// for the denominator
    pub fn of_f32_bits(bits: u32) -> Option<Self> {
        let is_neg = (bits & (1 << 31)) != 0;
        UnsignedRaw3232::of_f32_bits(bits & !(1 << 31)).map(|raw| Self { is_neg, raw })
    }

    /// Get the bit-value of this RationalN as an f64
    ///
    /// This can be converted to an f64 with f64::from_bits()
    pub(crate) fn as_f64_bits(&self) -> Option<u64> {
        self.raw
            .as_f64_bits()
            .map(|r| if self.is_neg { r | (1 << 63) } else { r })
    }

    /// Get the bit-value of this RationalN as an f32
    ///
    /// This can be converted to an f32 with f32::from_bits()
    pub(crate) fn as_f32_bits(&self) -> Option<u32> {
        self.raw
            .as_f32_bits()
            .map(|r| if self.is_neg { r | (1 << 31) } else { r })
    }

    /// Calculate the mantissa and exponent. Ignores sign.
    ///
    /// For zero, return (0,0); for 1 return (1<<63,0); 2 returns (1<<63,1), half (1<<63,-1)
    ///
    /// Shift numerator up until top bit is set, and denominator up until top bit is set
    /// if denominator is bigger than numerator then shift it right by one
    /// shift denominator right by 64
    ///
    /// Divide shifted numerator by shifted numerator; the dividend is the mantissa (could do some rounding...)
    pub fn to_mantissa64_exp(&self) -> (u64, i32) {
        self.raw.to_mantissa64_exp()
    }

    fn do_add_sub(&mut self, other: &Self, is_sub: bool) {
        if (self.is_neg == other.is_neg) == is_sub {
            let (r, b) = self.raw.overflowing_sub(&other.raw);
            if b {
                self.is_neg = !self.is_neg;
                todo!();
            }
            self.raw = r;
        } else {
            self.raw += other.raw;
        }
    }

    #[inline(always)]
    pub(crate) fn do_add(&mut self, other: &Self) {
        self.do_add_sub(other, false);
    }
    #[inline(always)]
    pub(crate) fn do_sub(&mut self, other: &Self) {
        self.do_add_sub(other, true);
    }
    #[inline(always)]
    pub(crate) fn do_bit_and(&mut self, other: &Self) {
        self.is_neg &= other.is_neg;
        self.raw &= other.raw;
    }
    #[inline(always)]
    pub(crate) fn do_bit_or(&mut self, other: &Self) {
        self.is_neg |= other.is_neg;
        self.raw |= other.raw;
    }
    #[inline(always)]
    pub(crate) fn do_bit_xor(&mut self, other: &Self) {
        self.is_neg ^= other.is_neg;
        self.raw ^= other.raw;
    }

    #[inline(always)]
    pub(crate) fn do_mul(&mut self, other: &Self) {
        self.raw.do_mul(&other.raw);
        self.is_neg ^= other.is_neg;
    }

    #[inline(always)]
    pub(crate) fn do_div(&mut self, other: &Self) {
        self.raw.do_div(&other.raw);
        self.is_neg ^= other.is_neg;
    }

    #[inline(always)]
    pub(crate) fn do_rem(&mut self, _other: &Self) {
        todo!();
    }
}
impl std::ops::Neg for SignedRaw3232 {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        self.is_neg = !self.is_neg;
        self
    }
}

impl Num for SignedRaw3232 {
    type FromStrRadixErr = ();
    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err(())
    }
}

impl Zero for SignedRaw3232 {
    fn is_zero(&self) -> bool {
        self == &Self::ZERO
    }
    fn zero() -> Self {
        Self::ZERO
    }
}
impl ConstZero for SignedRaw3232 {
    const ZERO: Self = Self {
        is_neg: false,
        raw: UnsignedRaw3232::ZERO,
    };
}

impl One for SignedRaw3232 {
    fn is_one(&self) -> bool {
        self == &Self::ONE
    }
    fn one() -> Self {
        Self::ONE
    }
}
impl ConstOne for SignedRaw3232 {
    const ONE: Self = Self {
        is_neg: false,
        raw: UnsignedRaw3232::ONE,
    };
}

impl FromPrimitive for SignedRaw3232 {
    fn from_i64(v: i64) -> Option<Self> {
        v.try_into().ok()
    }
    fn from_u64(v: u64) -> Option<Self> {
        v.try_into().ok()
    }
    fn from_f64(v: f64) -> Option<Self> {
        v.try_into().ok()
    }
    fn from_f32(v: f32) -> Option<Self> {
        v.try_into().ok()
    }
}

impl ToPrimitive for SignedRaw3232 {
    fn to_i64(&self) -> Option<i64> {
        Some(self.int() as i64)
    }
    fn to_u64(&self) -> Option<u64> {
        Some(self.int() as u64)
    }
    fn to_f64(&self) -> Option<f64> {
        (*self).try_into().ok()
    }
    fn to_f32(&self) -> Option<f32> {
        (*self).try_into().ok()
    }
}
impl NumCast for SignedRaw3232 {
    fn from<N: ToPrimitive>(n: N) -> Option<Self> {
        n.to_f64().map(|n| n.try_into().ok()).flatten()
    }
}

impl Signed for SignedRaw3232 {
    fn abs(&self) -> Self {
        let mut s = *self;
        s.is_neg = false;
        s
    }
    fn abs_sub(&self, other: &Self) -> Self {
        (*self - *other).abs()
    }
    fn signum(&self) -> Self {
        if self.is_neg() {
            -Self::ONE
        } else {
            Self::ONE
        }
    }
    fn is_positive(&self) -> bool {
        !self.is_neg()
    }
    fn is_negative(&self) -> bool {
        self.is_neg()
    }
}
