use super::{IntN, UIntN};
use crate::utils;
use num_traits::{ConstOne, One, Zero};

/// A Rational number, with numerator and denominator as N*64-bit integers
/// The type is Copy
///
/// This is not highly optimized for performance; it is designed to permit
/// precise operations of algorithms that require num_traits::Num
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RationalN<const N: usize> {
    /// Numerator - signed integer
    numer: IntN<N>,
    /// Denominator - unsigned integer
    denom: UIntN<N>,
}

impl<const N: usize> std::default::Default for RationalN<N> {
    fn default() -> Self {
        Self {
            numer: IntN::default(),
            denom: UIntN::ONE,
        }
    }
}

impl<const N: usize> std::convert::From<(IntN<N>, UIntN<N>)> for RationalN<N> {
    fn from((numer, denom): (IntN<N>, UIntN<N>)) -> Self {
        Self { numer, denom }
    }
}

impl<const N: usize> std::convert::From<u64> for RationalN<N> {
    fn from(value: u64) -> Self {
        Self {
            numer: value.into(),
            denom: UIntN::ONE,
        }
    }
}

impl<const N: usize> std::convert::From<i64> for RationalN<N> {
    fn from(value: i64) -> Self {
        Self {
            numer: value.into(),
            denom: UIntN::ONE,
        }
    }
}

impl<const N: usize> std::convert::From<(i64, u64)> for RationalN<N> {
    fn from((numer, denom): (i64, u64)) -> Self {
        Self {
            numer: numer.into(),
            denom: denom.into(),
        }
    }
}

impl<const N: usize> std::convert::From<f32> for RationalN<N> {
    #[inline]
    fn from(f: f32) -> Self {
        Self::of_f32_bits(f.to_bits()).unwrap()
    }
}

impl<const N: usize> std::convert::From<f64> for RationalN<N> {
    #[inline]
    fn from(f: f64) -> Self {
        Self::of_f64_bits(f.to_bits()).unwrap()
    }
}

impl<const N: usize> std::convert::From<RationalN<N>> for f32 {
    #[inline]
    fn from(r: RationalN<N>) -> f32 {
        (&r).into()
    }
}

impl<const N: usize> std::convert::From<RationalN<N>> for f64 {
    #[inline]
    fn from(r: RationalN<N>) -> f64 {
        (&r).into()
    }
}

impl<const N: usize> std::convert::From<&RationalN<N>> for f32 {
    fn from(r: &RationalN<N>) -> f32 {
        r.as_f32_bits().map(f32::from_bits).unwrap_or(f32::NAN)
    }
}

impl<const N: usize> std::convert::From<&RationalN<N>> for f64 {
    fn from(r: &RationalN<N>) -> f64 {
        r.as_f64_bits().map(f64::from_bits).unwrap_or(f64::NAN)
    }
}

impl<const N: usize> std::cmp::Ord for RationalN<N> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self.is_neg(), other.is_neg()) {
            (true, true) => {
                let sn_od = other.denom * self.numer.magnitude();
                let on_sd = self.denom * other.numer.magnitude();
                on_sd.cmp(&sn_od)
            }
            (false, false) => {
                let sn_od = other.denom * self.numer.magnitude();
                let on_sd = self.denom * other.numer.magnitude();
                sn_od.cmp(&on_sd)
            }
            (false, true) => std::cmp::Ordering::Greater,
            (true, false) => std::cmp::Ordering::Less,
        }
    }
}

impl<const N: usize> std::cmp::PartialOrd for RationalN<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> std::fmt::Display for RationalN<N> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        if self.denom().is_one() {
            self.numer.fmt(fmt)
        } else {
            let s = format!("{}/{}", self.numer, self.denom);
            fmt.pad(&s)
        }
    }
}

impl<const N: usize> std::ops::Neg for RationalN<N> {
    type Output = Self;
    fn neg(mut self) -> Self {
        self.numer = -self.numer;
        self
    }
}

impl<const N: usize> std::ops::Add for RationalN<N> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.do_add_sub(&other, false)
    }
}

impl<const N: usize> std::ops::AddAssign for RationalN<N> {
    fn add_assign(&mut self, other: Self) {
        *self = self.do_add_sub(&other, false)
    }
}

impl<const N: usize> std::ops::Sub for RationalN<N> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.do_add_sub(&other, true)
    }
}

impl<const N: usize> std::ops::SubAssign for RationalN<N> {
    fn sub_assign(&mut self, other: Self) {
        *self = self.do_add_sub(&other, true);
    }
}

impl<const N: usize> std::ops::Mul for RationalN<N> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.do_multiply(&other)
    }
}

impl<const N: usize> std::ops::MulAssign for RationalN<N> {
    fn mul_assign(&mut self, other: Self) {
        *self = self.do_multiply(&other);
    }
}

impl<const N: usize> std::ops::Div for RationalN<N> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self.do_divide(&other)
    }
}

impl<const N: usize> std::ops::DivAssign for RationalN<N> {
    fn div_assign(&mut self, other: Self) {
        *self = self.do_divide(&other);
    }
}

impl<const N: usize> std::ops::Rem for RationalN<N> {
    type Output = Self;

    fn rem(self, _other: Self) -> Self {
        todo!();
    }
}

impl<const N: usize> std::ops::RemAssign for RationalN<N> {
    fn rem_assign(&mut self, _other: Self) {
        todo!();
    }
}

impl<const N: usize> num_traits::FromPrimitive for RationalN<N> {
    fn from_i64(n: i64) -> Option<Self> {
        IntN::from_i64(n).map(|s| Self {
            numer: s,
            denom: num_traits::identities::ConstOne::ONE,
        })
    }
    fn from_u64(n: u64) -> Option<Self> {
        IntN::from_u64(n).map(|s| Self {
            numer: s,
            denom: num_traits::identities::ConstOne::ONE,
        })
    }
}

impl<const N: usize> num_traits::identities::One for RationalN<N> {
    fn one() -> Self {
        Self {
            numer: num_traits::identities::ConstOne::ONE,
            denom: num_traits::identities::ConstOne::ONE,
        }
    }
}

impl<const N: usize> num_traits::identities::Zero for RationalN<N> {
    fn zero() -> Self {
        Self::default()
    }
    fn is_zero(&self) -> bool {
        self.numer.is_zero()
    }
}

impl<const N: usize> num_traits::identities::ConstZero for RationalN<N> {
    const ZERO: Self = Self {
        numer: IntN::ZERO,
        denom: UIntN::ONE,
    };
}

impl<const N: usize> num_traits::identities::ConstOne for RationalN<N> {
    const ONE: Self = Self {
        numer: IntN::ONE,
        denom: UIntN::ONE,
    };
}

impl<const N: usize> RationalN<N> {
    /// Return a reference to the magnitude of the numerator of the rational
    pub fn numer(&self) -> &UIntN<N> {
        self.numer.magnitude()
    }

    /// Return a reference to the denominator of the rational
    pub fn denom(&self) -> &UIntN<N> {
        &self.denom
    }

    /// Return the sign of the rational
    pub fn is_neg(&self) -> bool {
        self.numer.is_neg()
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
        let (numer, numer_top_bit) = self.numer().most_significant_u128();
        if numer == 0 {
            return (0, 0);
        };
        let (mut denom, denom_top_bit) = self.denom().most_significant_u128();
        let mut exponent = (numer_top_bit as i32) - (denom_top_bit as i32);
        // eprintln!("denom , numer {} {denom} {numer} ", denom > numer);
        if denom > numer {
            exponent -= 1;
            denom >>= 1;
        }
        denom >>= 63;
        // eprintln!("{numer} / {denom}, {exponent}");
        ((numer / denom) as u64, exponent)
    }

    /// Normalize a RationalN by reducing numer/denom by GCD, setting zero
    /// correctly, etc
    ///
    /// If the RationalN were more optimized (such as having flags for zero,
    /// integer, and exponent, etc) then this would set those appropriately
    pub fn normalize(&mut self) -> bool {
        let gcd = self.numer().gcd(&self.denom);
        if gcd.is_zero() {
            if !self.denom.is_one() {
                self.denom = UIntN::<N>::ONE;
                true
            } else {
                false
            }
        } else if gcd.is_one() {
            false
        } else {
            self.numer = self.numer / gcd;
            self.denom /= gcd;
            true
        }
    }

    /// Take a mantissa (bit 63 set if not zero) and exponent (0=> divide by 2^63)
    /// to generate a (normalized) rational
    pub fn of_sign_mantissa_exp(is_neg: bool, mantissa: u64, exp: i32) -> Option<Self> {
        // eprintln!("of_sign_me: {is_neg} {mantissa:x} {exp}");
        if mantissa == 0 {
            Some(Self::default())
        } else {
            // bottom_bit is the right-most bit of mantissa that is set
            // num_bits is the number of bits set in mantissa (top bit downwards)
            //
            // for 1 (mantissa = 1<<63, exp=0); num_bits=1, bottom_bit=63
            // for 3 (mantissa = 3<<62, exp=1); num_bits=2, bottom_bit=62
            // for int N of num_bits mantissa = N<<(63-num_bits), exp=num_bits-1
            let bottom_bit = (0..64).find(|i| (mantissa >> (*i) & 1) != 0).unwrap();
            let num_bits = 64 - bottom_bit;
            let mantissa = mantissa >> (bottom_bit as u64);
            let left_shift = exp + 1 - num_bits;
            if left_shift == 0 {
                Some(mantissa.into())
            } else if left_shift + num_bits >= 64 * (N as i32) {
                None
            } else if left_shift > 0 {
                let mut numer: UIntN<N> = mantissa.into();
                numer.shift_left(left_shift as u32);
                let denom = UIntN::ONE;
                Some(Self {
                    numer: (is_neg, numer).into(),
                    denom,
                })
            } else if (-left_shift) >= 64 * (N as i32) {
                None
            } else {
                let numer: UIntN<N> = mantissa.into();
                let mut denom = UIntN::ONE;
                denom.shift_left((-left_shift) as u32);
                Some(Self {
                    numer: (is_neg, numer).into(),
                    denom,
                })
            }
        }
    }

    /// Get a RationalN from f64 bits, if it fits in the range
    /// for the numerator/denominator, by using powers of two
    /// for the denominator
    pub fn of_f64_bits(bits: u64) -> Option<Self> {
        // eprintln!("of_f64_bits: {bits:x} {}", f64::from_bits(bits));
        if bits & (u64::MAX >> 1) == 0 {
            Some(Self::default())
        } else {
            let is_neg = (bits & (1 << 63)) != 0;
            let mantissa = (1 << 52) | (bits & ((1 << 52) - 1));
            let exp = (bits >> 52) & 0x7ff;
            let exp = (exp as i32) - 1023;
            Self::of_sign_mantissa_exp(is_neg, mantissa << 11, exp)
        }
    }

    /// Get a RationalN from f32 bits, if it fits in the range
    /// for the numerator/denominator, by using powers of two
    /// for the denominator
    pub fn of_f32_bits(bits: u32) -> Option<Self> {
        // eprintln!("of_f32_bits: {bits:x} {}", f32::from_bits(bits));
        if bits & (u32::MAX >> 1) == 0 {
            Some(Self::default())
        } else {
            let is_neg = (bits & (1 << 31)) != 0;
            let mantissa = (1 << 23) | (bits & ((1 << 23) - 1));
            let exp = (bits >> 23) & 0xff;
            let exp = (exp as i32) - 127;
            Self::of_sign_mantissa_exp(is_neg, (mantissa as u64) << 40, exp)
        }
    }

    /// Get the bit-value of this RationalN as an f64
    ///
    /// This can be converted to an f64 with f64::from_bits()
    pub fn as_f64_bits(&self) -> Option<u64> {
        let mut result = 0;
        if self.is_neg() {
            result |= 1 << 63;
        }
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
    pub fn as_f32_bits(&self) -> Option<u32> {
        let mut result = 0;
        if self.is_neg() {
            result |= 1 << 31;
        }
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

    fn do_add_sub(&self, other: &Self, negate: bool) -> Self {
        let denom_gcd = self.denom.gcd(&other.denom);
        let self_denom_no_gcd: IntN<_> = (negate, self.denom / denom_gcd).into();
        let other_denom_no_gcd = other.denom / denom_gcd;
        let mut numer = self.numer * other_denom_no_gcd;
        // eprintln!(
        //     "{self} +-? {other} {numer} {} {}",
        //     other.numer, self_denom_no_gcd
        // );
        numer += other.numer * self_denom_no_gcd;
        let mut denom = self.denom * other_denom_no_gcd;
        if numer.is_zero() {
            Self::default()
        } else {
            let gcd = numer.magnitude().gcd(&denom);
            numer /= gcd.into();
            denom /= gcd;
            Self { numer, denom }
        }
    }

    /// Multiply A/B by C/D
    ///
    /// Find g0 = gcd(A,D); g1 = gcd(B,C); A = a.g0; B = b.g1; C = c.g1 ; D = d.g0
    ///
    /// A/B * C/D = AC / BD = a.g0.c.g1 / b.g1.d.g0 = a.c / b.d
    fn do_multiply(&self, other: &Self) -> Self {
        if self.numer.is_zero() || other.numer.is_zero() {
            return Self::default();
        }
        let g0 = self.numer.magnitude().gcd(&other.denom);
        let g1 = other.numer.magnitude().gcd(&self.denom);
        let a = self.numer / g0;
        let b = self.denom / g1;
        let c = other.numer / g1;
        let d = other.denom / g0;
        let numer = *a.magnitude() * c.magnitude();
        let denom = b * d;
        Self {
            numer: (self.numer.is_neg() != other.numer.is_neg(), numer).into(),
            denom,
        }
    }

    /// Divide A/B by D/C = multiply A/B by C/D
    ///
    /// Find g0 = gcd(A,D); g1 = gcd(B,C); A = a.g0; B = b.g1; C = c.g1 ; D = d.g0
    ///
    /// A/B * C/D = AC / BD = a.g0.c.g1 / b.g1.d.g0 = a.c / b.d
    fn do_divide(&self, other: &Self) -> Self {
        assert!(!other.numer.is_zero(), "Division by zero");
        if self.numer.is_zero() {
            return *self;
        }

        let g0 = self.numer.magnitude().gcd(other.numer.magnitude());
        let g1 = other.denom.gcd(&self.denom);
        let a = self.numer / g0;
        let b = self.denom / g1;
        let c = other.denom / g1;
        let d = other.numer / g0;
        let numer = c * a.magnitude();
        let denom = b * d.magnitude();
        Self {
            numer: (self.numer.is_neg() != other.numer.is_neg(), numer).into(),
            denom,
        }
    }

    /// Generate a string for an array, and return the LCM
    ///
    /// Should only include this if alloc
    pub fn with_common_denom<'a, I: Iterator<Item = &'a Self> + Clone + 'a>(
        iter: I,
    ) -> (impl Iterator<Item = IntN<N>> + 'a, UIntN<N>) {
        let lcm = iter
            .clone()
            .fold(UIntN::ONE, |lcm, v| lcm.lcm(&v.denom).unwrap());
        (
            iter.map(move |v| {
                let m = lcm / v.denom;
                v.numer * m
            }),
            lcm,
        )
    }
}

impl<const N: usize> num_traits::Num for RationalN<N> {
    type FromStrRadixErr = std::num::ParseIntError;
    fn from_str_radix(src: &str, radix: u32) -> Result<Self, std::num::ParseIntError> {
        if let Some((n, d)) = src.split_once('/') {
            let numer = IntN::from_str_radix(n, radix)?;
            let denom = UIntN::from_str_radix(d, radix)?;
            Ok(Self { numer, denom })
        } else {
            let numer = IntN::from_str_radix(src, radix)?;
            Ok(Self {
                numer,
                denom: UIntN::ONE,
            })
        }
    }
}

impl<const N: usize> crate::NumOps for RationalN<N> {
    fn of_f64(v: f64) -> Self {
        v.into()
    }

    fn is_unreliable_divisor(self) -> bool {
        self.numer.magnitude().find_top_bit_set() + 32 < self.denom.find_top_bit_set()
    }

    fn is_sign_negative(self) -> bool {
        self.numer.is_neg()
    }

    /// Return an estimate of the square root to within a precision
    fn sqrt_est(self) -> Self {
        utils::sqrt_est::<_, 5>(self, true)
    }

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn cbrt_est(self) -> Self {
        utils::cbrt_est::<_, 5>(self)
    }
}
