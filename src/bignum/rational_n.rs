use super::{IntN, UIntN};
use num_traits::{ConstOne, ConstZero, One, Zero};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RationalN<const N: usize> {
    numer: IntN<N>,
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
            write!(fmt, "{}/{}", self.numer, self.denom)
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
    pub fn numer(&self) -> &UIntN<N> {
        self.numer.magnitude()
    }
    pub fn denom(&self) -> &UIntN<N> {
        &self.denom
    }
    pub fn is_neg(&self) -> bool {
        self.numer.is_neg()
    }

    fn do_add_sub(&self, other: &Self, negate: bool) -> Self {
        let denom_gcd = self.denom.gcd(&other.denom);
        let self_denom_no_gcd: IntN<_> = (negate, self.denom / denom_gcd).into();
        let other_denom_no_gcd = other.denom / denom_gcd;
        let mut numer = self.numer * other_denom_no_gcd;
        eprintln!(
            "{self} +-? {other} {numer} {} {}",
            other.numer, self_denom_no_gcd
        );
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

        let g0 = self.numer.magnitude().gcd(&other.numer.magnitude());
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
