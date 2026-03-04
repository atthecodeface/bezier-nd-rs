use num_traits::{sign, Zero};

use crate::bignum::UIntN;

/// Calculate PI as a fixed point value with *3* bits of integer value
///
/// e.g. UIntN<2>, with 128 bits, has the top three bits set as '011 = 3' from '3.1415...'
///
/// This uses the BBP digit extraction algorithm from Wikipedia (Simon Plouffe, 1995)
pub fn pi_scaled_by<const N: usize>() -> UIntN<N> {
    let mut sum = UIntN::<N>::default();
    let num_bits = UIntN::<N>::num_bits();
    let four_scale = UIntN::<N>::with_bit_set(num_bits - 1);
    let two_scale = UIntN::<N>::with_bit_set(num_bits - 2);
    let one_scale = UIntN::<N>::with_bit_set(num_bits - 3);
    for k in 0..(num_bits as u64 / 4) {
        let s0 = four_scale / (8 * k + 1).into();
        let s1 = two_scale / (8 * k + 4).into();
        let s2 = one_scale / (8 * k + 5).into();
        let s3 = one_scale / (8 * k + 6).into();
        let mut s = s0 - s1 - s2 - s3;
        s.shift_right(4 * (k as u32));
        sum += s;
    }
    sum
}

pub fn e_scaled_by<const N: usize>() -> UIntN<N> {
    let mut sum = UIntN::<N>::default();
    let num_bits = UIntN::<N>::num_bits();
    let one_scale = UIntN::<N>::with_bit_set(num_bits - 2);
    let mut r_factorial = one_scale;
    for k in 0.. {
        sum += r_factorial;
        r_factorial /= (k + 1).into();
        if r_factorial.is_zero() {
            break;
        }
    }
    sum
}

/// Calculate the natural log of two
///
/// The calculation is of -ln(1 + (-1/2))
///
/// Uses the Taylor expansion of ln(1+x) = -2^(1) - 2^(-2)/2 - 2^(-3)/3 - 2^(-4)/4 ...
pub fn ln_two_scaled_by<const N: usize>() -> UIntN<N> {
    let num_bits = UIntN::<N>::num_bits();
    let one_scale = UIntN::<N>::with_bit_set(num_bits - 1);

    let mut sum = UIntN::<N>::default();

    for k in 1.. {
        let mut s = one_scale;
        s = s / k.into();
        s.shift_right((k - 1) as u32);
        if s.is_zero() {
            break;
        }
        sum += s;
    }
    sum
}

/// Calculate log(1+2^-power) * 2^scale
///
/// Uses the Taylor expansion, with a specified scaling for the base.
///
/// For natural logarithm, base can be zero
///
/// Use the Taylor expansion at 1 with x=2^-power yields 2^-power - 2^(-2*poewr)/2 + 2^(-3*power)/3 - 2^(-4*power)/4 ...
///
pub fn ln_one_plus_two_neg_power<const N: usize>(power: u32, base: UIntN<N>) -> UIntN<N> {
    let num_bits = UIntN::<N>::num_bits();
    let one_scale = UIntN::<N>::with_bit_set(num_bits - 1);

    let mut sum = UIntN::<N>::default();
    let mut value = {
        if base.is_zero() {
            one_scale
        } else {
            base
        }
    };
    value.shift_right(power);
    let mut sign_positive = true;
    let mut k = 1;
    loop {
        if sign_positive {
            sum += value / k.into();
        } else {
            sum -= value / k.into();
        }
        value.shift_right(power);
        if value.is_zero() {
            break;
        }
        k += 1;
        sign_positive = !sign_positive;
    }
    sum
}

/// Calculate atan(2^(-power))
///
/// Uses the Taylor expansion
///
/// arctan(x) = x - x^3/3 + x^5/5 - x^7/7 etc
///
pub fn atan_two_neg_power<const N: usize>(power: u32) -> UIntN<N> {
    assert!(power!=0, "atan(1) does not converge fast enough; this is PI/2, so use a scaled version of PI instead");
    let num_bits = UIntN::<N>::num_bits();
    let one_scale = UIntN::<N>::with_bit_set(num_bits - 1);

    let mut sum = UIntN::<N>::default();
    let mut value = { one_scale };
    value.shift_right(power);
    let mut sign_positive = true;
    let mut k = 1;
    loop {
        if sign_positive {
            sum += value / k.into();
        } else {
            sum -= value / k.into();
        }
        value.shift_right(power * 2);
        if value.is_zero() {
            break;
        }
        k += 2;
        sign_positive = !sign_positive;
    }
    sum
}
