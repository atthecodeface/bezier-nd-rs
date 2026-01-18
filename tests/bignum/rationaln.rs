use bezier_nd::bignum::RationalN;
use geo_nd::Num;
use num_traits::{one, zero, ConstZero};

fn prime_25519<N: Num + From<u64>>() -> N {
    let a: N = 2.into();
    let a2_4 = a * a * a * a;
    let a2_16 = a2_4 * a2_4 * a2_4 * a2_4;
    let a2_64 = a2_16 * a2_16 * a2_16 * a2_16;
    let a = a2_64;
    let a = a * a;
    let a = a * a;
    // a = 2^256
    let a = a / 2.into();
    let a = a - 19.into();
    a
}

fn test_add_subtract_n<const N: usize>(a: RationalN<N>) {
    let b = a + a;
    let c = b + a;

    assert_eq!(c - b, a);
    assert_eq!(c - a, b);
    assert_eq!(b - a, a);
    assert_eq!(a - a, 0_i64.into());
    assert_eq!(a + a, b);

    assert_eq!(b - c, -a);
    assert_eq!(a - c, -b);
    assert_eq!(a - b, -a);

    let m_a = -a;
    let m_b = -b;
    let m_c = m_a - b;

    assert_eq!(b + m_c, m_a);
    assert_eq!(a + m_c, m_b);
    assert_eq!(a + m_b, m_a);

    assert_eq!(m_c + b, m_a);
    assert_eq!(m_c + a, m_b);
    assert_eq!(m_b + a, m_a);

    assert_eq!(m_a + a, 0_u64.into());

    assert!(b > a);
    assert!(a < b);

    assert!(c > a);
    assert!(a < c);

    assert!(c > b);
    assert!(b < c);

    assert!(m_b < m_a);
    assert!(m_a > m_b);

    assert!(m_c < m_a);
    assert!(m_a > m_c);

    assert!(m_c < m_b);
    assert!(m_b > m_c);

    assert!(m_b < a);
    assert!(a > m_b);

    use num_traits::Num;
    assert_eq!(RationalN::from_str_radix(&a.to_string(), 10).unwrap(), a);
    assert_eq!(
        RationalN::from_str_radix(&m_a.to_string(), 10).unwrap(),
        m_a
    );
}

fn test_multiply_small_n<const N: usize>(a: RationalN<N>) {
    let zero: RationalN<N> = zero();
    let one: RationalN<N> = one();
    let two: RationalN<N> = 2_u64.into();
    let three: RationalN<N> = 3_u64.into();
    assert_eq!(a + a, a * two);
    assert_eq!(a + a + a, a * three);
    assert_eq!(a + a, two * a);
    assert_eq!(a + a + a, three * a);
    assert_eq!(three * a - a, two * a);

    assert!(two * a > a);
    assert!(a < three * a);

    assert_eq!(zero / a, zero);
    assert_eq!(zero * a, zero);
    assert_eq!(one * a, a);

    assert_eq!(a.numer().gcd(a.denom()), 1.into());
}

fn test_number_n<const N: usize>(a: RationalN<N>) {
    test_add_subtract_n(a);
    test_multiply_small_n(a);
    use num_traits::Num;
    assert!(RationalN::<N>::from_str_radix("", 10).is_err());
    assert!(RationalN::<N>::from_str_radix("x", 10).is_err());
}

#[test]
fn test_rational_tiny_int() {
    test_number_n::<1>(1_u64.into());
    test_number_n::<1>((u64::MAX >> 2).into());

    test_number_n::<8>(1_u64.into());
    test_number_n::<8>((u64::MAX >> 2).into());
}

#[test]
fn test_rational_tiny_fract() {
    use num_traits::FromPrimitive;
    let half = RationalN::<1>::from_u64(1).unwrap() / RationalN::<1>::from_u64(2).unwrap();
    let big: RationalN<1> = (u64::MAX - 3).into();
    let big = big * half * half;
    test_number_n(half);
    test_number_n(big);

    let half = RationalN::<8>::from_u64(1).unwrap() / RationalN::<8>::from_u64(2).unwrap();
    let big: RationalN<8> = u64::MAX.into();
    test_number_n(half);
    test_number_n(big);
}

#[test]
fn test_int_big() {
    let a = prime_25519::<RationalN<8>>();
    test_number_n(a);

    let a = prime_25519::<RationalN<5>>();
    test_number_n(a);
}

#[test]
#[should_panic]
fn test_int_add_overflow() {
    let a: RationalN<1> = u64::MAX.into();
    let one: RationalN<_> = 1_i64.into();
    let _b = a + one;
}

#[test]
#[should_panic]
fn test_int_mul_overflow_from_last() {
    let a: RationalN<1> = u64::MAX.into();
    let two: RationalN<_> = 2_i64.into();
    let _b = two * a;
}

#[test]
#[should_panic]
fn test_int_mul_overflow_from_middle() {
    let a: RationalN<4> = u64::MAX.into();
    let a2 = a * a;
    let a3 = a * a2;
    let _b = a2 * a3;
}

#[test]
#[should_panic]
fn test_int_div_zero() {
    let a = prime_25519::<RationalN<8>>();
    let _b = a / RationalN::ZERO;
}

#[test]
#[should_panic]
fn test_int_div_assign_zero() {
    let mut a = prime_25519::<RationalN<8>>();
    a /= RationalN::ZERO;
}

#[test]
#[should_panic]
fn test_int_rem_zero() {
    let a = prime_25519::<RationalN<8>>();
    let _b = a % RationalN::ZERO;
}

#[test]
#[should_panic]
fn test_int_rem_assign_zero() {
    let mut a = prime_25519::<RationalN<8>>();
    a %= RationalN::ZERO;
}
