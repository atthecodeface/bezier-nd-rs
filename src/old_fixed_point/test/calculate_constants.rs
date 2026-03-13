use num_traits::{FromPrimitive, Zero};

use crate::bignum::UIntN;

/// Calculate PI as a fixed point value with *3* bits of integer value
///
/// e.g. UIntN<2>, with 128 bits, has the top three bits set as '011 = 3' from '3.1415...'
///
/// This uses the BBP digit extraction algorithm from Wikipedia (Simon Plouffe, 1995)
pub fn pi_scaled_by<const N: usize>(num_frac_bits: u32) -> UIntN<N> {
    assert!(
        num_frac_bits + 3 <= UIntN::<N>::num_bits(),
        "Request for PI to more fractional bits than the value can hold"
    );
    let mut sum = UIntN::<N>::default();
    let four_scale = UIntN::<N>::with_bit_set(num_frac_bits + 2);
    let two_scale = UIntN::<N>::with_bit_set(num_frac_bits + 1);
    let one_scale = UIntN::<N>::with_bit_set(num_frac_bits);
    for k in 0..=(num_frac_bits as u64 / 4) {
        let s0 = four_scale / UIntN::<N>::from_u64(8 * k + 1).unwrap();
        let s1 = two_scale / UIntN::<N>::from_u64(8 * k + 4).unwrap();
        let s2 = one_scale / UIntN::<N>::from_u64(8 * k + 5).unwrap();
        let s3 = one_scale / UIntN::<N>::from_u64(8 * k + 6).unwrap();
        let mut s = s0 - s1 - s2 - s3;
        s >>= 4 * k;
        sum += s;
    }
    sum
}

/// Calculate E as a fixed point value with *2* bits of integer value
///
/// e.g. UIntN<2>, with 128 bits, has the top two bits set as '10 = 2' from '2.718...'
///
/// This uses the (for x=1) the series 1 + x + x^2/2! + x^3/3! ...
pub fn e_scaled_by<const N: usize>(num_frac_bits: u32) -> UIntN<N> {
    assert!(
        num_frac_bits + 2 <= UIntN::<N>::num_bits(),
        "Request for E to more fractional bits than the value can hold"
    );
    let mut sum = UIntN::<N>::default();
    let one_scale = UIntN::<N>::with_bit_set(num_frac_bits);
    let mut r_factorial = one_scale;
    for k in 0.. {
        sum += r_factorial;
        r_factorial /= UIntN::<N>::from_u32(k + 1).unwrap();
        if r_factorial.is_zero() {
            break;
        }
    }
    sum
}

/// Calculate the natural log of two; this is a pure fraction
///
/// The calculation is of -ln(1 + (-1/2))
///
/// Uses the Taylor expansion of ln(1+x) = -2^(1) - 2^(-2)/2 - 2^(-3)/3 - 2^(-4)/4 ...
pub fn ln_two_scaled_by<const N: usize>(num_frac_bits: u32) -> UIntN<N> {
    assert!(
        num_frac_bits <= UIntN::<N>::num_bits(),
        "Request for ln2(e) to more fractional bits than the value can hold"
    );
    let mut sum = UIntN::<N>::default();
    for k in 1.. {
        let s = UIntN::<N>::with_bit_set(num_frac_bits - k) / UIntN::<N>::from_u32(k).unwrap();
        if s.is_zero() {
            break;
        }
        sum += s;
    }
    sum
}

/// Calculate log(1+2^-power) either given a base with the appropriate number of
/// fractional bits, or as a natural logarithm with the specified number of
/// fractional bits
///
/// Uses the Taylor expansion, with a specified scaling for the base.
///
/// For natural logarithm `num_frac_bits` is used, and base should be None
///
/// Use the Taylor expansion at 1 with x=2^-power yields 2^-power - 2^(-2*poewr)/2 + 2^(-3*power)/3 - 2^(-4*power)/4 ...
///
pub fn ln_one_plus_two_neg_power<const N: usize>(
    num_frac_bits: u32,
    power: u32,
    base: Option<UIntN<N>>,
) -> UIntN<N> {
    let mut sum = UIntN::<N>::default();
    let mut value = {
        if let Some(base) = base {
            base >> power
        } else {
            assert!(
                num_frac_bits <= UIntN::<N>::num_bits(),
                "Request for ln2(1+2^-power) to more fractional bits than the value can hold"
            );
            UIntN::<N>::with_bit_set(num_frac_bits - power)
        }
    };

    let mut sign_positive = true;
    let mut k: u32 = 1;
    loop {
        if sign_positive {
            sum += value / UIntN::<N>::from_u32(k).unwrap();
        } else {
            sum -= value / UIntN::<N>::from_u32(k).unwrap();
        }
        value >>= power;
        if value.is_zero() {
            break;
        }
        k += 1;
        sign_positive = !sign_positive;
    }
    sum
}

/// Calculate atan(2^(-power)) with the specified number of fractional bits
///
/// Uses the Taylor expansion
///
/// arctan(x) = x - x^3/3 + x^5/5 - x^7/7 etc
///
pub fn atan_two_neg_power<const N: usize>(num_frac_bits: u32, power: u32) -> UIntN<N> {
    assert!(power!=0, "atan(1) as the series expansion used does not converge fast enough; this is PI/2, so use a scaled version of PI instead");
    let mut sum = UIntN::<N>::default();
    let mut value = UIntN::<N>::with_bit_set(num_frac_bits - power);
    let mut sign_positive = true;
    let mut k = 1;
    loop {
        if sign_positive {
            sum += value / UIntN::<N>::from_u32(k).unwrap();
        } else {
            sum -= value / UIntN::<N>::from_u32(k).unwrap();
        }
        value >>= power * 2;
        if value.is_zero() {
            break;
        }
        k += 2;
        sign_positive = !sign_positive;
    }
    sum
}
