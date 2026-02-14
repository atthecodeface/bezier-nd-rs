use super::UIntN;

/// A signed integer, of +- 64*N bits (stored with sign and magnitude),
/// supporting copy
///
/// This is not optimized for performance, but to enable use of algorithms
/// that require num_traits::Num
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IntN<const N: usize> {
    /// Sign of the number
    is_neg: bool,
    /// Magnitude of the number
    value: UIntN<N>,
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

impl<const N: usize> std::ops::Neg for IntN<N> {
    type Output = Self;
    fn neg(mut self) -> Self {
        if !self.value_is_zero() {
            self.is_neg = !self.is_neg;
        }
        self
    }
}

impl<const N: usize> std::ops::Add for IntN<N> {
    type Output = Self;

    #[track_caller]
    fn add(self, other: Self) -> Self {
        self.do_add_sub(&other, false)
    }
}

impl<const N: usize> std::ops::AddAssign for IntN<N> {
    #[track_caller]
    fn add_assign(&mut self, other: Self) {
        *self = self.do_add_sub(&other, false);
    }
}

impl<const N: usize> std::ops::Sub for IntN<N> {
    type Output = Self;

    #[track_caller]
    fn sub(self, other: Self) -> Self {
        self.do_add_sub(&other, true)
    }
}

impl<const N: usize> std::ops::SubAssign for IntN<N> {
    #[track_caller]
    fn sub_assign(&mut self, other: Self) {
        *self = self.do_add_sub(&other, true);
    }
}

impl<const N: usize> std::ops::Mul for IntN<N> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.do_multiply(&other)
    }
}

impl<const N: usize> std::ops::Mul<&IntN<N>> for IntN<N> {
    type Output = Self;

    fn mul(self, other: &IntN<N>) -> Self {
        self.do_multiply(other)
    }
}

impl<const N: usize> std::ops::MulAssign for IntN<N> {
    fn mul_assign(&mut self, other: Self) {
        *self = self.do_multiply(&other);
    }
}

impl<const N: usize> std::ops::Div for IntN<N> {
    type Output = Self;

    #[track_caller]
    fn div(self, other: Self) -> Self {
        let Some((d, _r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        d
    }
}

impl<const N: usize> std::ops::DivAssign for IntN<N> {
    fn div_assign(&mut self, other: Self) {
        let Some((d, _r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        *self = d;
    }
}

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

impl<const N: usize> std::ops::Rem for IntN<N> {
    type Output = Self;
    #[track_caller]
    fn rem(self, other: Self) -> Self {
        let Some((_d, r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        r
    }
}

impl<const N: usize> std::ops::RemAssign for IntN<N> {
    #[track_caller]
    fn rem_assign(&mut self, other: Self) {
        let Some((_d, r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        *self = r;
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

impl<const N: usize> IntN<N> {
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
    fn do_add_sub(mut self, other: &Self, negate_other: bool) -> Self {
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
            if self.value.subtract(&other.value) {
                self.is_neg = !self.is_neg;
                self.value.twos_complement();
            }
            if self.value_is_zero() {
                self.is_neg = false;
            }
        } else {
            self.value += other.value;
        }
        self
    }

    fn do_multiply(&self, other: &Self) -> Self {
        Self {
            value: self.value * other.value,
            is_neg: self.is_neg != other.is_neg,
        }
    }

    /// Calculate the mantissa and exponent. Ignores sign.
    ///
    /// For zero, return (0,0); for 1 return (1<<63,0); 2 returns (1<<63,1)
    pub fn to_mantissa64_exp(&self) -> (u64, u32) {
        // (1,0) for 1; (0,0) for 0; (u128::MAX, 127) for u128::MAX
        let (numer, exponent) = self.value.most_significant_u128();
        if numer == 0 {
            return (0, 0);
        };
        ((numer >> 64) as u64, exponent)
    }

    /// Get the bit-value of this IntN as an f64
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
        } else if exp > 1023 {
            None
        } else {
            result |= (((exp + 1023) & 0x7ff) as u64) << 52;
            result |= (mantissa >> 11) & ((1 << 52) - 1);
            Some(result)
        }
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
