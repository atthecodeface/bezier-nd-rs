use num_traits::{ConstOne, ConstZero, Zero};
use std::fmt::Write;


/// Unsigned integer of N*64 bits, supporting copy
///
/// This is not heavily optimized for performance, but for space; an optimization
/// might include flags for 'is zero', 'is one', and possibly the number of u64 values
/// that are non-zero; instead this is *purely* the array of u64 that make up the value
///
/// The purpose of this 'bignum' is to enable larger numbers simply for algorithms that
/// require num_traits::Num
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UIntN<const N: usize> {
    /// Value is stored with most significant at index 0 so that PartialOrd operates correctly
    value: [u64; N],
}

impl<const N: usize> std::default::Default for UIntN<N> {
    fn default() -> Self {
        Self { value: [0; N] }
    }
}

impl<const N: usize> std::ops::Deref for UIntN<N> {
    type Target = [u64; N];
    fn deref(&self) -> &Self::Target {
        &self.value
    }
}

impl<const N: usize> std::convert::TryFrom<i64> for UIntN<N> {
    type Error = std::num::TryFromIntError;
    fn try_from(value: i64) -> Result<Self, std::num::TryFromIntError> {
        if value < 0 {
            Err(<u64 as std::convert::TryFrom<i64>>::try_from(-1).unwrap_err())
        } else {
            Ok((value as u64).into())
        }
    }
}

impl<const N: usize> std::convert::From<u64> for UIntN<N> {
    fn from(value: u64) -> Self {
        let mut s = Self::default();
        s.value[N - 1] = value;
        s
    }
}

impl<const N: usize> std::cmp::PartialOrd for UIntN<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<const N: usize> std::cmp::Ord for UIntN<N> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.value.cmp(&other.value)
    }
}

pub struct UIntNDigitIter<const N: usize> {
    remaining: UIntN<N>,
    scan: UIntN<N>,
    radix: u32,
    radix_s: UIntN<N>,
    is_last: bool,
    is_complete: bool,
}

impl<const N: usize> UIntNDigitIter<N> {
    pub fn new(value: UIntN<N>, radix: u32) -> Self {
        if value.is_zero() {
            Self {
                is_last: true,
                is_complete: false,
                radix,
                remaining: UIntN::ZERO,
                scan: UIntN::ZERO,
                radix_s: UIntN::ZERO,
            }
        } else {
            let radix_s: UIntN<_> = (radix as u64).into();
            let mut scan = radix_s;
            while scan < value {
                let (overflow, scan_times_radix) = scan.multiply_by_u64(radix as u64);
                if overflow {
                    break;
                }
                scan = scan_times_radix;
            }
            while scan > value {
                scan /= radix_s;
            }
            Self {
                is_last: false,
                is_complete: false,
                radix,
                remaining: value,
                scan,
                radix_s,
            }
        }
    }
}

impl<const N: usize> std::iter::Iterator for UIntNDigitIter<N> {
    type Item = u32;
    fn next(&mut self) -> Option<u32> {
        if self.is_complete {
            None
        } else if self.is_last {
            self.is_complete = true;
            Some(self.remaining.value[N - 1] as u32)
        } else if self.scan.cmp_u64(self.radix as u64) == std::cmp::Ordering::Less {
            self.is_last = true;
            self.next()
        } else {
            let (div, rem) = self.remaining.do_div_rem(&self.scan).unwrap();
            let digit = div.value[N - 1] as u32;
            self.scan /= self.radix_s;
            self.remaining = rem;
            Some(digit)
        }
    }
}

impl<const N: usize> std::fmt::Display for UIntN<N> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        for c in self.as_digits(10) {
            fmt.write_char(c)?;
        }
        Ok(())
    }
}

impl<const N: usize> std::ops::Neg for UIntN<N> {
    type Output = Self;
    fn neg(self) -> Self {
        panic!("UIntN does not support negation");
    }
}

impl<const N: usize> std::ops::Add for UIntN<N> {
    type Output = Self;

    #[track_caller]
    fn add(self, other: Self) -> Self {
        self.do_add(&other)
    }
}

impl<const N: usize> std::ops::AddAssign for UIntN<N> {
    #[track_caller]
    fn add_assign(&mut self, other: Self) {
        *self = self.do_add(&other);
    }
}

impl<const N: usize> std::ops::Sub for UIntN<N> {
    type Output = Self;

    #[track_caller]
    fn sub(self, other: Self) -> Self {
        self.do_subtract(&other)
    }
}

impl<const N: usize> std::ops::SubAssign for UIntN<N> {
    #[track_caller]
    fn sub_assign(&mut self, other: Self) {
        *self = self.do_subtract(&other);
    }
}

impl<const N: usize> std::ops::Mul<&UIntN<N>> for UIntN<N> {
    type Output = UIntN<N>;

    #[track_caller]
    fn mul(self, other: &UIntN<N>) -> UIntN<N> {
        self.do_multiply(other)
    }
}

impl<const N: usize> std::ops::Mul for UIntN<N> {
    type Output = Self;

    #[track_caller]
    fn mul(self, other: Self) -> Self {
        self.do_multiply(&other)
    }
}

impl<const N: usize> std::ops::MulAssign for UIntN<N> {
    #[track_caller]
    fn mul_assign(&mut self, other: Self) {
        *self = self.do_multiply(&other);
    }
}

impl<const N: usize> std::ops::Div for UIntN<N> {
    type Output = Self;

    #[track_caller]
    fn div(self, other: Self) -> Self {
        let Some((d, _r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        d
    }
}

impl<const N: usize> std::ops::DivAssign for UIntN<N> {
    #[track_caller]
    fn div_assign(&mut self, other: Self) {
        // eprintln!("Divide {:?} {:?}", self.value, other.value);
        let Some((d, _r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        // eprintln!("Result {:?} {:?}", d, _r);
        *self = d;
    }
}

impl<const N: usize> std::ops::Rem for UIntN<N> {
    type Output = Self;

    #[track_caller]
    fn rem(self, other: Self) -> Self {
        let Some((_d, r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        r
    }
}

impl<const N: usize> std::ops::RemAssign for UIntN<N> {
    #[track_caller]
    fn rem_assign(&mut self, other: Self) {
        let Some((_d, r)) = self.do_div_rem(&other) else {
            panic!("Division by zero");
        };
        *self = r;
    }
}

impl<const N: usize> num_traits::FromPrimitive for UIntN<N> {
    fn from_i64(n: i64) -> Option<Self> {
        n.try_into().ok()
    }
    fn from_u64(n: u64) -> Option<Self> {
        Some(n.into())
    }
}

impl<const N: usize> num_traits::identities::One for UIntN<N> {
    fn one() -> Self {
        Self::ONE
    }
}

impl<const N: usize> num_traits::identities::Zero for UIntN<N> {
    fn zero() -> Self {
        Self::ZERO
    }
    fn is_zero(&self) -> bool {
        self.value_is_zero()
    }
}

impl<const N: usize> num_traits::identities::ConstZero for UIntN<N> {
    const ZERO: Self = Self { value: [0; N] };
}

impl<const N: usize> num_traits::identities::ConstOne for UIntN<N> {
    const ONE: Self = {
        let mut s = Self { value: [0; N] };
        s.value[N - 1] = 1;
        s
    };
}

impl<const N: usize> UIntN<N> {
    /// Generate an iterator of char for the digits of the number given a radix
    pub fn as_digits(&self, radix: u32) -> impl Iterator<Item = char> {
        UIntNDigitIter::new(*self, radix).map(move |c| char::from_digit(c, radix).unwrap())
    }

    fn value_is_zero(&self) -> bool {
        self.value.iter().all(|s| *s == 0)
    }

    /// Find top bit set + 1
    pub(crate) fn find_top_bit_set(&self) -> u32 {
        if let Some(n) = self.most_significant_n() {
            ((N as u32) - n - 1) * 64 + Self::find_top_bit(self.value[n as usize]).unwrap()
        } else {
            0
        }
    }

    /// Find top bit of u64
    ///
    /// Return false if v is zero
    fn find_top_bit(v: u64) -> Option<u32> {
        (0..64).rev().find(|i| ((v >> i) & 1) != 0)
    }

    /// Return most significant n where `self.value[n]` is nonzero
    ///
    /// Returns None if self is zero
    fn most_significant_n(&self) -> Option<u32> {
        self.value
            .iter()
            .enumerate()
            .find(|(_, v)| (**v) != 0)
            .map(|(i, _)| i as u32)
    }

    /// Amount to shift left self by (in bits) to make its top bit line up with others top bit
    ///
    /// Returns None if either value is zero, or if self's top bit is already to the left of others
    fn left_shift_amount(&self, other: &Self) -> Option<u32> {
        // eprintln!("Find left shift amount to get {self:?} to {other:?}");
        let Some(v_i) = self.most_significant_n() else {
            // eprintln!("No most significant for {self:?}");
            return None;
        };
        let Some(o_i) = other.most_significant_n() else {
            // eprintln!("No most significant for {other:?}");
            return None;
        };
        if o_i > v_i {
            // eprintln!("{o_i} > {v_i} - no shift needed");
            None
        } else {
            // Self.value[v_i] is the most significant (lowest v_i) that is nonzero
            // other.value[o_i] is the most significant (lowest o_i) that is nonzero
            //
            // v_i can be bigger than o_i but not less than
            let o_tb = Self::find_top_bit(other.value[o_i as usize]).unwrap();
            let v_tb = Self::find_top_bit(self.value[v_i as usize]).unwrap();
            // eprintln!("Other i/tb, self i/tb : {o_i} . {o_tb} , {v_i} . {v_tb}");
            if o_i == v_i && v_tb > o_tb {
                None
            } else if o_i == v_i {
                Some(o_tb - v_tb)
            } else {
                Some((v_i - o_i) * 64 + o_tb - v_tb)
            }
        }
    }

    /// Find the most significant 128 bits and the shift
    /// to apply
    ///
    /// The result is:
    /// * (0,0) for zero
    /// * (0x80000000_00000000_00000000_00000000,0) for one
    /// * (0x80000000_00000000_00000000_00000000,64) for 1<<64
    /// * (0xc0000000_00000000_00000000_00000000,1) for 3
    pub fn most_significant_u128(&self) -> (u128, u32) {
        let Some(i) = self.most_significant_n() else {
            return (0, 0);
        };
        let i = i as usize;
        let bit = Self::find_top_bit(self.value[i]).unwrap() as usize;
        // If self==1 then i=N-1 and bit=0
        let top_bit = 64 * (N - 1 - i) + bit;
        let mut top_bits = self.value[i] as u128;
        // top_bits is all the bits in value[i], throwing away the top '127-bit' zeros
        //
        // top_bits will have bits[127-bit..127] valid
        top_bits <<= 127 - bit;
        if i + 1 < N {
            // If there are more bits, shift them up and or them in to make [63-bit..127] valid
            top_bits |= (self.value[i + 1] as u128) << (63 - bit);
        }
        if i + 2 < N && bit < 63 {
            // If there are even more bits, shift them up by -1-bit (i.e. right by bit+1) to make [0..127] valid
            top_bits |= (self.value[i + 2] as u128) >> (bit + 1);
        }
        (top_bits, top_bit as u32)
    }

    /// Shift left by n bits
    ///
    /// Return true if it overflows
    #[track_caller]
    pub fn shift_left(&mut self, amount: u32) -> bool {
        if amount == 0 {
            return false;
        }
        let bit_shift = amount & 63;
        let byte_shift = (amount >> 6) as usize;
        if bit_shift == 0 {
            let overflow = self.value[0..byte_shift].iter().any(|v| *v != 0);
            self.value.copy_within(byte_shift..N, 0);
            self.value[N - byte_shift..N].fill(0);
            overflow
        } else {
            let right_shift = 64 - bit_shift;
            let mut v_i = self.value[byte_shift];
            for i in 0..N {
                let v_ip1 = {
                    if i + byte_shift + 1 < N {
                        self.value[i + byte_shift + 1]
                    } else {
                        0
                    }
                };
                let shift_in = v_ip1.wrapping_shr(right_shift);
                self.value[i] = v_i.wrapping_shl(bit_shift) | shift_in;
                v_i = v_ip1;
            }
            false
        }
    }

    /// Shift right by 64 bits
    pub fn shift_right_64(&mut self) {
        self.value.copy_within(0..N - 1, 1);
        self.value[0] = 0;
    }

    /// Shift right by one bit
    pub fn shift_right_one(&mut self) {
        let mut shift_in = false;
        for i in 0..N {
            let v_i = self.value[i];
            let shift_out = (v_i & 1) != 0;
            self.value[i] = v_i.wrapping_shr(1) | {
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
    fn set_bit(&mut self, n: u32) {
        let bit = 1 << (n & 63);
        self.value[N - 1 - (n >> 6) as usize] |= bit;
    }

    /// Divide value by divisor which MUST NOT BE ZERO
    ///
    /// Put the remainder in rem
    pub fn do_div_rem(&self, divisor: &Self) -> Option<(Self, Self)> {
        // eprintln!("Do_div_rem {self:?} {divisor:?}");
        if divisor.value_is_zero() {
            return None;
        }
        let mut div = Self::default();
        let mut rem = *self;
        let mut thing = *divisor;
        let Some(mut bit_num) = thing.left_shift_amount(&rem) else {
            // eprintln!("{thing:?} {rem:?} no bit_num ");
            return Some((div, rem));
        };
        // eprintln!("need to shift left {thing:?} by {bit_num}");
        thing.shift_left(bit_num);
        // eprintln!("thing now {thing:?} rem {rem:?}");
        loop {
            // eprintln!("loop : {thing:?} rem {rem:?} bit_num:{bit_num}");
            if rem >= thing {
                rem -= thing;
                div.set_bit(bit_num);
            }
            thing.shift_right_one();
            if bit_num == 0 {
                break;
            }
            bit_num -= 1;
        }
        Some((div, rem))
    }

    /// Multiply v0*v and accumulate into self.value[..i] indicating overflow
    fn mul_acc_value_step(&mut self, v0: &[u64; N], v: u64, i: usize) -> bool {
        let mut overflow = 0;
        for (j, p) in v0.iter().rev().enumerate() {
            let beyond_n = (i + j) >= N;
            let r_in = {
                if beyond_n {
                    0
                } else {
                    self.value[N - 1 - (i + j)]
                }
            };
            let (mul_acc_out, overflow_out) = (*p).carrying_mul_add(v, r_in, overflow);
            if beyond_n {
                if mul_acc_out != 0 {
                    return true;
                }
            } else {
                self.value[N - 1 - (i + j)] = mul_acc_out;
            }
            overflow = overflow_out;
        }
        overflow != 0
    }

    // Multiply two values and return true if overflow
    fn multiply_value(&self, other: &Self) -> (bool, Self) {
        let mut result = Self::default();
        for (i, v) in other.value.iter().rev().enumerate() {
            if result.mul_acc_value_step(&self.value, *v, i) {
                return (true, result);
            }
        }
        (false, result)
    }

    // Multiply by a u64 and return true if overflow
    fn multiply_by_u64(&self, value: u64) -> (bool, Self) {
        let mut result = Self::default();
        if result.mul_acc_value_step(&self.value, value, 0) {
            return (true, result);
        }
        (false, result)
    }

    #[track_caller]
    fn do_add(mut self, other: &Self) -> Self {
        let mut carry = false;
        for (p0, p1) in self.value.iter_mut().rev().zip(other.value.iter().rev()) {
            let (r, carry_out) = (*p0).carrying_add(*p1, carry);
            *p0 = r;
            carry = carry_out;
        }
        assert!(!carry, "Addition overflowed");
        self
    }

    /// Negate the value using twos complement
    pub fn twos_complement(&mut self) {
        let mut carry = true;
        for v in self.value.iter_mut().rev() {
            let (r, carry_out) = (!(*v)).carrying_add(0, carry);
            *v = r;
            carry = carry_out;
        }
    }

    /// Subtract other value from self
    pub fn subtract(&mut self, other: &Self) -> bool {
        let mut borrow = false;
        for (p0, p1) in self.value.iter_mut().rev().zip(other.value.iter().rev()) {
            let (r, borrow_out) = (*p0).borrowing_sub(*p1, borrow);
            *p0 = r;
            borrow = borrow_out;
        }
        borrow
    }

    fn do_subtract(mut self, other: &Self) -> Self {
        let mut borrow = false;
        for (p0, p1) in self.value.iter_mut().rev().zip(other.value.iter().rev()) {
            let (r, borrow_out) = (*p0).borrowing_sub(*p1, borrow);
            *p0 = r;
            borrow = borrow_out;
        }
        assert!(!borrow, "Subtraction underflowed");
        self
    }

    #[track_caller]
    fn do_multiply(&self, other: &Self) -> Self {
        let (overflow, r) = self.multiply_value(other);
        assert!(!overflow, "Multiplication overflowed");
        r
    }

    /// Calculate the GCD of two UIntN
    pub fn gcd(&self, other: &Self) -> Self {
        if self.value_is_zero() || other.value_is_zero() {
            Self::default()
        } else {
            let (mut a, mut b) = match self.cmp(other) {
                std::cmp::Ordering::Equal => {
                    return *self;
                }
                std::cmp::Ordering::Less => (*other, *self),
                std::cmp::Ordering::Greater => (*self, *other),
            };
            loop {
                let rem = a % b;
                if rem.value_is_zero() {
                    return b;
                }
                let rem2 = b % rem;
                if rem2.value_is_zero() {
                    return rem;
                }
                a = rem;
                b = rem2;
            }
        }
    }

    /// Calculate the LCMof two UIntN
    pub fn lcm(&self, other: &Self) -> Option<Self> {
        if self.value_is_zero() || other.value_is_zero() {
            None
        } else {
            let gcd = self.gcd(other);
            let result = *self / gcd;
            Some(result * other)
        }
    }

    /// Compare this UIntN with a pure u64 value
    pub fn cmp_u64(&self, value: u64) -> std::cmp::Ordering {
        if self.value[0..N - 1].iter().all(|s| *s == 0) {
            self.value[N - 1].cmp(&value)
        } else {
            std::cmp::Ordering::Greater
        }
    }
}

impl<const N: usize> num_traits::Num for UIntN<N> {
    type FromStrRadixErr = std::num::ParseIntError;
    fn from_str_radix(src: &str, radix: u32) -> Result<Self, std::num::ParseIntError> {
        let mut result = Self::default();
        let radix_s: Self = (radix as u64).into();
        let mut has_chars = false;
        for c in src.chars() {
            let Some(d) = c.to_digit(radix) else {
                return Err(i32::from_str_radix("a", 10).unwrap_err());
            };
            let d: Self = (d as u64).into();
            result = result * radix_s + d;
            has_chars = true;
        }
        if has_chars {
            Ok(result)
        } else {
            Err("".parse::<i32>().unwrap_err())
        }
    }
}

#[test]
fn test_most_significant_u128() {
    let x = UIntN::<2>::from(u64::MAX);
    assert_eq!(UIntN::<2>::from(0_u64).most_significant_u128(), (0, 0));
    assert_eq!(
        UIntN::<2>::from(1_u64).most_significant_u128(),
        (0x80000000_00000000_00000000_00000000_u128, 0)
    );
    assert_eq!(
        UIntN::<2>::from(2_u64).most_significant_u128(),
        (0x80000000_00000000_00000000_00000000_u128, 1)
    );
    assert_eq!(
        UIntN::<2>::from(3_u64).most_significant_u128(),
        (0xc0000000_00000000_00000000_00000000_u128, 1)
    );
    assert_eq!(
        UIntN::<2>::from(u64::MAX).most_significant_u128(),
        (0xffffffff_ffffffff_00000000_00000000_u128, 63)
    );
    assert_eq!(
        (x + x).most_significant_u128(),
        (0xffffffff_ffffffff_00000000_00000000_u128, 64)
    );
}
