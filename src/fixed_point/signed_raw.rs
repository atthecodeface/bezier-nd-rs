use super::FPRaw;
use super::UnsignedRaw3232;
use num_traits::{ConstOne, ConstZero, FromPrimitive, Num, NumCast, One, ToPrimitive, Zero};

// Value is (i + f/2^32) * sign(is_neg)
#[derive(Debug, Clone, Copy, Default, PartialEq, Hash)]
pub struct SignedRaw3232 {
    is_neg: bool,
    raw: UnsignedRaw3232,
}

#[test]
fn test_raw() -> Result<(), Box<dyn std::error::Error>> {
    let mut r = UnsignedRaw3232::ZERO;
    let mut s = UnsignedRaw3232::ONE;
    assert!(r.is_zero());
    assert!(!r.is_one());
    assert!(!s.is_zero());
    assert!(s.is_one());

    r += UnsignedRaw3232::ONE;
    assert!(!r.is_zero());
    assert!(r.is_one());
    r += r;
    assert!(!r.is_zero());
    assert!(!r.is_one());

    r *= r;
    r -= r;
    assert!(r.is_zero());
    assert!(!r.is_one());

    r = 1_f64.try_into().unwrap();
    assert!(!r.is_zero());
    assert!(r.is_one());

    assert_eq!(1.0_f64, r.into(), "From f64 value of 1 must be 1");
    r = 4.0_f64.try_into().unwrap();
    assert_eq!(4.0_f64, r.into(), "From f64 value of 4 must be 4");
    s = 0.25_f64.try_into().unwrap();
    assert_eq!(0.25_f64, s.into(), "From f64 value of 1/4 must be 1/4");

    r = r * s;
    assert!(!r.is_zero());
    assert!(r.is_one(), "1/4 times 4 is 1");

    s = 1.0_f64.try_into().unwrap();
    r = 4.0_f64.try_into().unwrap();
    r /= s;
    assert_eq!(4.0_f64, r.into(), "4/1 must be 4");

    s = s / r;
    assert_eq!(0.25_f64, s.into(), "1/4 must be 1/4");

    s *= s;
    assert_eq!(0.0625_f64, s.into(), "1/4/4 must be 1/16");

    r = r / s;
    assert_eq!(64.0_f64, r.into(), "4 / 1/16 must be 64");
    assert_eq!(r * UnsignedRaw3232::from_u32(5).unwrap(), 320_u32.into());

    assert_eq!(
        r / UnsignedRaw3232::from_u32(5).unwrap(),
        (64.0_f64 / 5.0).try_into().unwrap()
    );

    assert_eq!(
        r % UnsignedRaw3232::from_u32(5).unwrap(),
        (64.0_f64 % 5.0).try_into().unwrap()
    );
    Ok(())
}

impl SignedRaw3232 {
    /// Take a mantissa (bit 63 set if not zero) and exponent (0=> divide by 2^63)
    /// to generate a (normalized) Raw number
    pub(crate) fn of_sign_mantissa_exp(is_neg: bool, mantissa: u64, exp: i32) -> Option<Self> {
        UnsignedRaw3232::of_mantissa_exp(mantissa, exp).map(|raw| Self {
            is_neg: is_neg && !raw.is_zero(),
            raw,
        })
    }

    // (i0+f0*2^-32) * (i1+f1*2^-32) =
    // (i0*i1) + (f0*i1 + f1*i0)*2^-32 + (f0*f1)*2^-64
    #[inline(always)]
    pub(crate) fn do_mul(&mut self, other: &Self) {
        self.raw.do_mul(&other.raw);
        self.is_neg ^= other.is_neg;
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

    #[inline(always)]
    pub(crate) fn do_add(&mut self, other: &Self) {}
    #[inline(always)]
    pub(crate) fn do_sub(&mut self, other: &Self) {}
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
    pub(crate) fn do_div(&mut self, other: &Self) {}
    #[inline(always)]
    pub(crate) fn do_rem(&mut self, other: &Self) {}
}
