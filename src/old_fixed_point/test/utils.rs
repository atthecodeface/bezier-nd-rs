use super::{BackingKind, TestKind};

pub fn int_bits<T, I, const N: usize>(_t: &T) -> usize
where
    T: TestKind<I, N>,
    I: BackingKind,
{
    std::mem::size_of::<I>() * 8 - N - 1
}

/// Compare the most significant `num_bits` of the value with the `top_bits`
///
/// The integer_decode yields the most significant bits in 'm'; this can be
/// shifted until the top bit that is set is 'num_bits-1`
pub fn matches_float<T, I, const N: usize>(value: T, f: f64) -> bool
where
    T: TestKind<I, N> + num_traits::Float,
    I: BackingKind,
{
    let (m, e, s) = value.integer_decode();
    let value_f = (m as f64) * (s as f64) * 2.0_f64.powi(e as i32);
    let diff = (value_f - f).abs();
    let epsilon = 2.0_f64.powi(-(N as i32));
    diff <= epsilon
}

/// Compare values to be nearly equal - to within `frac_value*2^-N`
///
/// The value from integer_decode is m*2^e, so:
///
/// * if e==-N then compare m and frac_value
///
/// * if e+k==-N then compare m and frac_value<<k
///
/// * if e==k-N then compare m<<k and frac_value
pub fn nearly_equal<T, I, const N: usize>(t: T, e: T, frac_value: u64) -> bool
where
    T: TestKind<I, N> + num_traits::Float,
    I: BackingKind,
{
    let diff = (t - e).abs();
    let (m, e, s) = diff.integer_decode();
    assert_eq!(s, 1, "Sign of absolute value must be positive");
    let e_diff = e + N as i16;
    dbg!(m, e, s, N, frac_value, e_diff);
    if m == 0 {
        true
    } else if e_diff >= 0 {
        if m << e_diff == 0 {
            false
        } else {
            m << e_diff <= frac_value
        }
    } else {
        m <= frac_value << -e_diff
    }
}
