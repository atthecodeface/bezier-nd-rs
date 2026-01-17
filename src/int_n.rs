#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntN<const N: usize> {
    is_neg: bool,
    // Value is stored with most significant at [0] so that PartialOrd operates correctly
    value: [u64; N],
}

impl<const N: usize> std::default::Default for IntN<N> {
    fn default() -> Self {
        Self {
            is_neg: false,
            value: [0; N],
        }
    }
}

impl<const N: usize> std::cmp::PartialOrd for IntN<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self.is_neg, other.is_neg) {
            (true, true) => other.value.partial_cmp(&self.value),
            (false, false) => self.value.partial_cmp(&other.value),
            (false, true) => Some(std::cmp::Ordering::Greater),
            (true, false) => Some(std::cmp::Ordering::Less),
        }
    }
}

impl<const N: usize> std::fmt::Display for IntN<N> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        std::fmt::Debug::fmt(self, fmt)
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

    fn add(self, other: Self) -> Self {
        self.do_add(&other)
    }
}

impl<const N: usize> std::ops::AddAssign for IntN<N> {
    fn add_assign(&mut self, other: Self) {
        *self = self.do_add(&other);
    }
}

impl<const N: usize> std::ops::Sub for IntN<N> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.do_subtract(&other)
    }
}

impl<const N: usize> std::ops::SubAssign for IntN<N> {
    fn sub_assign(&mut self, other: Self) {
        *self = self.do_subtract(&other);
    }
}

impl<const N: usize> std::ops::Mul for IntN<N> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.do_multiply(&other)
    }
}

impl<const N: usize> std::ops::MulAssign for IntN<N> {
    fn mul_assign(&mut self, other: Self) {
        *self = self.do_multiply(&other);
    }
}

impl<const N: usize> std::ops::Div for IntN<N> {
    type Output = Self;

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

impl<const N: usize> std::ops::Rem for IntN<N> {
    type Output = Self;

    fn rem(self, other: Self) -> Self {
        let Some((_d, r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        r
    }
}

impl<const N: usize> std::ops::RemAssign for IntN<N> {
    fn rem_assign(&mut self, other: Self) {
        let Some((_d, r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        *self = r;
    }
}

impl<const N: usize> num_traits::FromPrimitive for IntN<N> {
    fn from_i64(n: i64) -> Option<Self> {
        let mut s = Self::default();
        if n < 0 {
            s.value[N - 1] = (-n) as u64;
        } else {
            s.value[N - 1] = n as u64;
        }
        Some(s)
    }
    fn from_u64(n: u64) -> Option<Self> {
        let mut s = Self::default();
        s.value[N - 1] = n;
        Some(s)
    }
}

impl<const N: usize> num_traits::identities::One for IntN<N> {
    fn one() -> Self {
        let mut s = Self::default();
        s.value[N - 1] = 1;
        s
    }
}

impl<const N: usize> num_traits::identities::Zero for IntN<N> {
    fn zero() -> Self {
        Self {
            is_neg: false,
            value: [0; N],
        }
    }
    fn is_zero(&self) -> bool {
        self.value_is_zero()
    }
}

impl<const N: usize> num_traits::identities::ConstZero for IntN<N> {
    const ZERO: Self = Self {
        is_neg: false,
        value: [0; N],
    };
}

impl<const N: usize> num_traits::identities::ConstOne for IntN<N> {
    const ONE: Self = {
        let mut s = Self {
            is_neg: false,
            value: [0; N],
        };
        s.value[N - 1] = 1;
        s
    };
}

impl<const N: usize> IntN<N> {
    fn value_is_zero(&self) -> bool {
        self.value.iter().all(|s| *s == 0)
    }

    /// Find top bit of u64
    ///
    /// Return false if v is zero
    fn find_top_bit(v: u64) -> Option<u32> {
        (0..64).rev().find(|i| ((v >> i) & 1) != 0)
    }

    /// Return n where self.value[n] is nonzero
    ///
    /// Returns None if self is zero
    fn most_significant_n(v: &[u64; N]) -> Option<u32> {
        v.iter()
            .enumerate()
            .find(|(_, v)| (**v) != 0)
            .map(|(i, _)| i as u32)
    }

    /// Amount to shift left v by (in bits) to make its top bit line up with others top bit
    ///
    /// Returns None if either value is zero, or if v's top bit is already to the left of others
    fn left_shift_amount(v: &[u64; N], other: &[u64; N]) -> Option<u32> {
        let Some(v_i) = Self::most_significant_n(v) else {
            return None;
        };
        let Some(o_i) = Self::most_significant_n(other) else {
            return None;
        };
        if o_i > v_i {
            return None;
        } else {
            let o_tb = Self::find_top_bit(other[o_i as usize]).unwrap();
            let v_tb = Self::find_top_bit(v[o_i as usize]).unwrap();
            if o_i == v_i && v_tb < o_tb {
                None
            } else if o_i == v_i {
                Some(v_tb - o_tb)
            } else {
                Some((v_i - o_i) * 64 + v_tb - o_tb)
            }
        }
    }

    /// Shift left by n bits
    fn shift_left(v: &mut [u64; N], amount: u32) -> bool {
        if amount == 0 {
            return false;
        }
        let bit_shift = amount & 63;
        let byte_shift = (amount >> 6) as usize;
        if bit_shift == 0 {
            let overflow = v[0..byte_shift].iter().any(|v| *v != 0);
            v.copy_within(byte_shift..N, 0);
            v[N - byte_shift..N].fill(0);
            overflow
        } else {
            let mut shift_in = 0;
            let right_shift = 64 - bit_shift;
            for i in 0..N {
                let v_i = v[i];
                let shift_out = v_i.wrapping_shl(bit_shift);
                let next_shift_in = v_i.wrapping_shr(right_shift);
                if i >= byte_shift {
                    v[i - byte_shift] = shift_out | shift_in;
                } else if shift_out != 0 || shift_in != 9 {
                    return true;
                }
                shift_in = next_shift_in;
            }
            false
        }
    }

    /// Shift right by one bit
    fn shift_right(v: &mut [u64; N]) {
        let mut shift_in = false;
        for i in 0..N {
            let v_i = v[i];
            let shift_out = (v_i & 1) != 0;
            v[i] = v_i.wrapping_shr(1) | {
                if shift_in {
                    0x8000_0000_0000_0000
                } else {
                    0
                }
            };
            shift_in = shift_out;
        }
    }

    /// Set bit n
    fn set_bit(v: &mut [u64; N], n: u32) {
        let bit = 1 << (n & 31);
        v[(n >> 5) as usize] |= bit;
    }

    /// Divide value by divisor which MUST NOT BE ZERO
    ///
    /// Put the remainder in rem
    fn div_rem(value: &[u64; N], divisor: &[u64; N], div: &mut [u64; N], rem: &mut [u64; N]) {
        *rem = *value;
        *div = [0; N];
        let mut thing = *divisor;
        let Some(mut bit_num) = Self::left_shift_amount(&thing, rem) else {
            return;
        };
        Self::shift_left(&mut thing, bit_num);
        while bit_num < (64 * N as u32) {
            if *rem > thing {
                assert!(
                    !Self::subtract_value(rem, &thing),
                    "Subtract cannot underflow as rem > thing"
                );
                Self::set_bit(div, bit_num);
            }
            Self::shift_right(&mut thing);
            bit_num += 1;
        }
    }

    /// 20 / 7 = 2 remainder 6
    /// 20 / -7 = -2 remainder 6
    /// -20 / 7 = -2 remainder -6
    /// -20 / -7 = 2 remainder -6
    fn do_div_rem(&self, other: &Self) -> Option<(Self, Self)> {
        if other.value_is_zero() {
            None
        } else {
            let mut div = Self::default();
            let mut rem = Self::default();
            Self::div_rem(&self.value, &other.value, &mut div.value, &mut rem.value);
            div.is_neg = self.is_neg != other.is_neg;
            rem.is_neg = self.is_neg;
            Some((div, rem))
        }
    }

    // Add two values and return true if overflow
    fn add_value(v0: &mut [u64; N], v1: &[u64; N]) -> bool {
        let mut carry = false;
        for (p0, p1) in v0.iter_mut().rev().zip(v1.iter().rev()) {
            let (r, carry_out) = (*p0).carrying_add(*p1, carry);
            *p0 = r;
            carry = carry_out;
        }
        carry
    }

    // Subtract two values and return true if borrow
    fn subtract_value(v0: &mut [u64; N], v1: &[u64; N]) -> bool {
        let mut borrow = false;
        for (p0, p1) in v0.iter_mut().rev().zip(v1.iter().rev()) {
            let (r, borrow_out) = (*p0).borrowing_sub(*p1, borrow);
            *p0 = r;
            borrow = borrow_out;
        }
        borrow
    }

    // Multiply v0*v and accumulate into r[..i] indicating overflow
    fn mul_acc_value_step(r: &mut [u64; N], v0: &[u64; N], v: u64, i: usize) -> bool {
        let mut overflow = 0;
        for (j, p) in v0.iter().enumerate().rev() {
            let beyond_n = (i + j) >= N;
            let r_in = {
                if beyond_n {
                    0
                } else {
                    r[i + j]
                }
            };
            let (overflow_out, mul_acc_out) = (*p).carrying_mul_add(v, r_in, overflow);
            if beyond_n {
                if mul_acc_out != 0 {
                    return true;
                }
            } else {
                r[i + j] = mul_acc_out;
            }
            overflow = overflow_out;
        }
        overflow != 0
    }

    // Multiply two values and return true if borrow
    fn multiply_value(r: &mut [u64; N], v0: &[u64; N], v1: &[u64; N]) -> bool {
        for (i, v) in v1.iter().rev().enumerate() {
            if Self::mul_acc_value_step(r, v0, *v, i) {
                return true;
            }
        }
        false
    }

    fn do_add(mut self, other: &Self) -> Self {
        if self.is_neg != other.is_neg {
            // if self is +1 and other is -3 then subtract yields a borrow
            // and self becomes [u64::MAX, .. u64::MAX, , u64::MAX-1]
            //
            // we need to replace with U64::MAX-v for each v
            if Self::subtract_value(&mut self.value, &other.value) {
                self.is_neg = !self.is_neg;
                for v in self.value.iter_mut() {
                    *v = u64::MAX - (*v);
                }
            }
            if self.value_is_zero() {
                self.is_neg = false;
            }
        } else {
            let overflow = Self::add_value(&mut self.value, &other.value);
            assert!(overflow, "Addition overflowed");
        }
        self
    }

    fn do_subtract(mut self, other: &Self) -> Self {
        if self.is_neg == other.is_neg {
            // if self is +1 and other is +3 then subtract yields a borrow
            // and self becomes [u64::MAX, .. u64::MAX, , u64::MAX-1]
            //
            // we need to replace with U64::MAX-v for each v
            if Self::subtract_value(&mut self.value, &other.value) {
                self.is_neg = !self.is_neg;
                for v in self.value.iter_mut() {
                    *v = u64::MAX - (*v);
                }
            }
            if self.value_is_zero() {
                self.is_neg = false;
            }
        } else {
            let overflow = Self::add_value(&mut self.value, &other.value);
            assert!(overflow, "Subtraction overflowed");
        }
        self
    }
    fn do_multiply(&self, other: &Self) -> Self {
        let mut r = Self::default();
        r.is_neg = self.is_neg != other.is_neg;
        let overflow = Self::multiply_value(&mut r.value, &self.value, &other.value);
        assert!(overflow, "Multiplication overflowed");
        r
    }
}

impl<const N: usize> num_traits::Num for IntN<N> {
    type FromStrRadixErr = std::num::ParseIntError;
    fn from_str_radix(src: &str, radix: u32) -> Result<Self, std::num::ParseIntError> {
        u64::from_str_radix(src, radix)
            .map(|s| <Self as num_traits::FromPrimitive>::from_u64(s).unwrap())
    }
}
