use num_traits::{ConstOne, ConstZero, FromPrimitive, Num, NumCast, One, ToPrimitive, Zero};

// Value is (i + f/2^32)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct UnsignedRaw3232 {
    i: u32,
    f: u32,
}

impl std::convert::From<u32> for UnsignedRaw3232 {
    fn from(i: u32) -> Self {
        Self { i, f: 0 }
    }
}

impl std::convert::From<(u32, u32)> for UnsignedRaw3232 {
    fn from((i, f): (u32, u32)) -> Self {
        Self { i, f }
    }
}
impl std::convert::TryFrom<i64> for UnsignedRaw3232 {
    type Error = ();
    #[inline]
    fn try_from(v: i64) -> Result<Self, ()> {
        if v >= 0 && v < 1 << 32 {
            Ok((v as u32).into())
        } else {
            Err(())
        }
    }
}

impl std::convert::TryFrom<u64> for UnsignedRaw3232 {
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

impl std::convert::TryFrom<f32> for UnsignedRaw3232 {
    type Error = ();
    #[inline]
    fn try_from(f: f32) -> Result<Self, ()> {
        Self::of_f32_bits(f.to_bits()).ok_or(())
    }
}

impl std::convert::TryFrom<f64> for UnsignedRaw3232 {
    type Error = ();
    #[inline]
    fn try_from(f: f64) -> Result<Self, ()> {
        Self::of_f64_bits(f.to_bits()).ok_or(())
    }
}

impl std::convert::From<UnsignedRaw3232> for f64 {
    #[inline]
    fn from(r: UnsignedRaw3232) -> f64 {
        r.as_f64_bits().map(f64::from_bits).unwrap_or(f64::NAN)
    }
}

impl std::convert::From<UnsignedRaw3232> for f32 {
    #[inline]
    fn from(r: UnsignedRaw3232) -> f32 {
        r.as_f32_bits().map(f32::from_bits).unwrap_or(f32::NAN)
    }
}

impl std::cmp::Ord for UnsignedRaw3232 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.i.cmp(&other.i) {
            std::cmp::Ordering::Equal => self.f.cmp(&other.f),
            ord => ord,
        }
    }
}

impl std::cmp::PartialOrd for UnsignedRaw3232 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl UnsignedRaw3232 {
    pub(crate) fn int(&self) -> u32 {
        self.i
    }
    pub(crate) fn frac(&self) -> u32 {
        self.f
    }

    /// Take a mantissa (bit 63 set if not zero) and exponent (0=> divide by 2^63)
    /// to generate a (normalized) Raw number
    ///
    /// The top integer bit set in the result will be is exponent
    ///
    /// If exponent > 31 then it is not possible to store the value
    ///
    /// If exponent is 31 then the result is integer=mantissa>>32, fraction = (mantissa>>0) & 32 bits
    ///
    /// If exponent is zero then the result is integer=mantissa>>63, fraction = (mantissa>>31) & 32 bits
    ///
    /// If exponent is -1 then the result is integer=0, fraction = (mantissa>>32) & 32 bits
    ///
    /// If exponent is -32 then the result is integer=0, fraction = (mantissa>>63) & 32 bits
    ///
    /// If exponent is less than -33 then the result is zero
    pub(crate) fn of_mantissa_exp(mantissa: u64, exp: i32) -> Option<Self> {
        if mantissa == 0 || exp < -33 {
            Some(Self::default())
        } else if exp > 31 {
            None
        } else if exp >= 0 {
            let i = mantissa >> (63 - exp);
            let f = mantissa >> (31 - exp);
            Some(Self {
                i: i as u32,
                f: f as u32,
            })
        } else {
            let f = mantissa >> (31 - exp);
            Some(Self { i: 0, f: f as u32 })
        }
    }

    /// Get a [UnsignedRaw3232] from f64 bits, if it fits in the range
    /// for the numerator/denominator, by using powers of two
    /// for the denominator
    pub fn of_f64_bits(bits: u64) -> Option<Self> {
        if bits & (u64::MAX >> 1) == 0 {
            Some(Self::default())
        } else {
            let is_neg = (bits & (1 << 63)) != 0;
            if is_neg {
                None
            } else {
                let mantissa = (1 << 52) | (bits & ((1 << 52) - 1));
                let exp = (bits >> 52) & 0x7ff;
                let exp = (exp as i32) - 1023;
                Self::of_mantissa_exp(mantissa << 11, exp)
            }
        }
    }

    /// Get an [UnsignedRaw3232] from f32 bits, if it fits in the range
    /// for the numerator/denominator, by using powers of two
    /// for the denominator
    pub fn of_f32_bits(bits: u32) -> Option<Self> {
        // eprintln!("of_f32_bits: {bits:x} {}", f32::from_bits(bits));
        if bits & (u32::MAX >> 1) == 0 {
            Some(Self::default())
        } else {
            let is_neg = (bits & (1 << 31)) != 0;
            if is_neg {
                None
            } else {
                let mantissa = (1 << 23) | (bits & ((1 << 23) - 1));
                let exp = (bits >> 23) & 0xff;
                let exp = (exp as i32) - 127;
                Self::of_mantissa_exp((mantissa as u64) << 40, exp)
            }
        }
    }

    /// Get the bit-value of this RationalN as an f64
    ///
    /// This can be converted to an f64 with f64::from_bits()
    pub(crate) fn as_f64_bits(&self) -> Option<u64> {
        let mut result = 0;
        let (mantissa, exp) = self.to_mantissa64_exp();
        if mantissa == 0 {
            Some(0)
        } else if !(-1022..=1023).contains(&exp) {
            None
        } else {
            let mut rounded_mantissa = mantissa >> 11;
            if (mantissa & ((1 << 11) - 1)) >= (1 << 10) {
                rounded_mantissa += 1;
            }
            result |= (((exp + 1023) & 0x7ff) as u64) << 52;
            result |= rounded_mantissa & ((1 << 52) - 1);
            Some(result)
        }
    }

    /// Get the bit-value of this RationalN as an f32
    ///
    /// This can be converted to an f32 with f32::from_bits()
    pub(crate) fn as_f32_bits(&self) -> Option<u32> {
        let mut result = 0;
        let (mantissa, exp) = self.to_mantissa64_exp();
        if mantissa == 0 {
            Some(0)
        } else if !(-126..=127).contains(&exp) {
            None
        } else {
            let mut rounded_mantissa = mantissa >> 40;
            if (mantissa & ((1 << 40) - 1)) >= (1 << 39) {
                rounded_mantissa += 1;
            }
            result |= (((exp + 127) & 0xff) as u32) << 23;
            result |= (rounded_mantissa & ((1 << 23) - 1)) as u32;
            Some(result)
        }
    }

    /// The exponent is the top bit set in 'i', or if the integer part is zero,
    /// the top bit set in 'f' minus 32. Hence the output is in the range -32 to 31.
    ///
    /// For the smallest number (i=0, f=1) the exponent is -32; for the largest number the exponent is 31.
    ///
    /// Hence shifting left by 31-exponent puts the top bit into bit 31 of the integer part
    fn exponent(&self) -> i32 {
        if self.i != 0 {
            (0..32_usize)
                .rev()
                .find(|i| ((self.i >> i) & 1) != 0)
                .unwrap() as i32
        } else if self.f != 0 {
            ((0..32).rev().find(|i| ((self.f >> i) & 1) != 0).unwrap() as i32) - 32
        } else {
            0
        }
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
        if self.i == 0 && self.f == 0 {
            return (0, 0);
        };
        let mantissa: u64 = ((self.i as u64) << 32) | (self.f as u64);
        let exponent = self.exponent();
        let mantissa = mantissa << ((31 - exponent) as usize);
        (mantissa, exponent)
    }

    pub(crate) fn overflowing_sub(mut self, other: &Self) -> (Self, bool) {
        let (f, b) = self.f.overflowing_sub(other.f);
        let (i, b) = self.i.borrowing_sub(other.i, b);
        self.f = f;
        self.i = i;
        (self, b)
    }
    pub(crate) fn overflowing_add(mut self, other: &Self) -> (Self, bool) {
        let (f, c) = self.f.overflowing_add(other.f);
        let (i, c) = self.i.carrying_add(other.i, c);
        self.f = f;
        self.i = i;
        (self, c)
    }
    #[inline(always)]
    pub(crate) fn do_add(&mut self, other: &Self) {
        let (s, c) = self.overflowing_add(other);
        assert!(!c, "Overflow in add");
        *self = s;
    }
    #[inline(always)]
    pub(crate) fn do_sub(&mut self, other: &Self) {
        let (s, b) = self.overflowing_sub(other);
        assert!(!b, "Underflow in sub");
        *self = s;
    }
    #[inline(always)]
    pub(crate) fn do_bit_and(&mut self, other: &Self) {
        self.i &= other.i;
        self.f &= other.f;
    }
    #[inline(always)]
    pub(crate) fn do_bit_or(&mut self, other: &Self) {
        self.i |= other.i;
        self.f |= other.f;
    }
    #[inline(always)]
    pub(crate) fn do_bit_xor(&mut self, other: &Self) {
        self.i ^= other.i;
        self.f ^= other.f;
    }
    // (i0+f0*2^-32) * (i1+f1*2^-32) =
    // (i0*i1) + (f0*i1 + f1*i0)*2^-32 + (f0*f1)*2^-64
    #[inline(always)]
    pub(crate) fn do_mul(&mut self, other: &Self) {
        let i0 = self.i as u64;
        let i1 = other.i as u64;
        let f0 = self.f as u64;
        let f1 = other.f as u64;
        let mut i_i = i0 * i1;
        let f_i_0 = i0 * f1;
        let f_i_1 = i1 * f0;
        let f_f = f0 * f1;
        let (f_i, c) = f_i_0.carrying_add(f_i_1, false);
        let (f_i, c2) = f_i.carrying_add(f_f >> 32, false);
        // eprintln!("{i0}.{f0} {i1}.{f1} {i_i} {f_i_0} {f_i_1}");
        if c {
            i_i += 1;
        }
        if c2 {
            i_i += 1;
        }
        i_i += f_i >> 32 ;
        // eprintln!("{i_i} {f_i}");
        self.i = i_i as u32;
        self.f = f_i as u32;
    }

    #[inline(always)]
    pub(crate) fn do_div(&mut self, other: &Self) {
        // self = sm * 2^se; other = om * 2^oe; result = sm/om * 2^(se-oe)
        //
        // se, oe are in the range -32 to +31; so se-oe is in the range -63 to +63
        let (sm, se) = self.to_mantissa64_exp();
        let (om, oe) = other.to_mantissa64_exp();
        assert_ne!(om, 0, "Division by zero");
        let sm = (sm as u128) << 64;
        let om = om as u128 ;
        // rm = sm/om * 2^64; this is in the range of just over (1/2)<<64 to just under (2)<<64
        //
        // Thos requires rm*2^32 so that the fractional part is the bottom 32 bits, and the integer part is above
        //
        // Must multiply by 2^(se-oe-64+32); se-oe-32 is in the range -95 to +31
        let mut rm = sm / om;
        let mut re = se - oe;
        if rm >= 1 << 64 {
            rm >>= 1;
            re += 1;
        }
        // rm is 1/2<<64 to u64::MAX; i.e. it has bit 63 set
        //
        // re == 0 implies rm = result * 2^64, so i is 0, f is rm >> 32
        //
        // re == 1 implies rm = result * 2^63, so i is rm >> 63, f is rm >> 31
        //
        // re == 32 implies rm = result * 2^32, so i is rm >> 32, f is rm >> 0
        //
        // re > 32 implies there is an integer overflow
        assert!(re <= 32, "Overflow in division");
        if re > 0 {
            rm <<= re;
        } else if re < -32 {
            rm = 0;
        } else if re < 0 {
            rm >>= (-re) as usize;
        }
        let i = (rm >> 64) as u32;
        let f = (rm >> 32) as u32;
        *self = Self { i, f };
    }

    // Rem x % y = x - (x / y).trunc() * y
    #[inline(always)]
    pub(crate) fn do_rem(&mut self, other: &Self) {
        let mut s = *self;
        s /= other;
        s.f = 0;
        s *= other;
        *self -= s;
    }
}

impl Num for UnsignedRaw3232 {
    type FromStrRadixErr = ();
    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err(())
    }
}

impl Zero for UnsignedRaw3232 {
    fn is_zero(&self) -> bool {
        self == &Self::ZERO
    }
    fn zero() -> Self {
        Self::ZERO
    }
}
impl ConstZero for UnsignedRaw3232 {
    const ZERO: Self = Self { i: 0, f: 0 };
}
impl One for UnsignedRaw3232 {
    fn is_one(&self) -> bool {
        self == &Self::ONE
    }
    fn one() -> Self {
        Self::ONE
    }
}
impl ConstOne for UnsignedRaw3232 {
    const ONE: Self = Self { i: 1, f: 0 };
}

impl FromPrimitive for UnsignedRaw3232 {
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

impl ToPrimitive for UnsignedRaw3232 {
    fn to_i64(&self) -> Option<i64> {
        Some(self.i as i64)
    }
    fn to_u64(&self) -> Option<u64> {
        Some(self.i as u64)
    }
    fn to_f64(&self) -> Option<f64> {
        (*self).try_into().ok()
    }
    fn to_f32(&self) -> Option<f32> {
        (*self).try_into().ok()
    }
}
impl NumCast for UnsignedRaw3232 {
    fn from<N: ToPrimitive>(n: N) -> Option<Self> {
        n.to_f64().map(|n| n.try_into().ok()).flatten()
    }
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
