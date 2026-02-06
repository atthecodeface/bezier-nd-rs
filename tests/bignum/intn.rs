use bezier_nd::bignum::IntN;
use geo_nd::Num;
use num_traits::{ConstOne, ConstZero};

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

fn test_add_subtract_n<N: Num + From<u64>>(a: N) {
    let b = a + a;
    let c = b + a;

    assert_eq!(c - b, a);
    assert_eq!(c - a, b);
    assert_eq!(b - a, a);
    assert_eq!(a - a, 0.into());
    assert_eq!(a + a, b);

    assert_eq!(b - c, -a);
    assert_eq!(a - c, -b);
    assert_eq!(a - b, -a);

    let mut c = b;
    c += a;
    assert_eq!(b - c, -a);
    assert_eq!(a - c, -b);
    c -= b;
    assert_eq!(c, a);
    c += b;

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
}

fn test_multiply_small_n<N: Num + From<u64>>(a: N) {
    let zero: N = N::zero();
    let one: N = N::one();
    let two: N = 2.into();
    let three: N = 3.into();
    assert_eq!(a + a, a * two);
    assert_eq!(a + a + a, a * three);
    assert_eq!(a + a, two * a);
    assert_eq!(a + a + a, three * a);
    assert_eq!(three * a - a, two * a);

    assert!(two * a > a);
    assert!(a < three * a);

    assert_eq!(zero / a, zero);
    assert_eq!(zero * a, zero);
    assert_eq!((one - one / a == zero), (a == one));
    assert_eq!(one * a, a);
}

fn test_number_n<N: Num + From<u64>>(a: N) {
    test_add_subtract_n(a);
    test_multiply_small_n(a);
    assert!(<N as num_traits::Num>::from_str_radix("", 10).is_err());
    assert!(<N as num_traits::Num>::from_str_radix("x", 10).is_err());
}

#[test]
fn test_int_tiny() {
    use num_traits::Num;
    assert_eq!(IntN::<1>::ZERO.to_string(), "0");

    test_number_n::<IntN<1>>(1_u64.into());
    test_number_n::<IntN<1>>((u64::MAX >> 2).into());

    test_number_n::<IntN<8>>(1_u64.into());
    test_number_n::<IntN<8>>((u64::MAX >> 2).into());

    let v: IntN<1> = (5_i64).into();
    let vn: IntN<1> = (-1_i64).into();
    let v = v + vn;
    let v = v - vn;
    assert_eq!(v.to_string(), "5");
    assert_eq!(vn.to_string(), "-1");
    let v = v * &v;
    let mut v = v * v;
    v *= v;
    assert_eq!(v.to_string(), "390625");
    test_number_n::<IntN<1>>(v);

    let v: IntN<1> = u64::MAX.into();
    assert_eq!(v.to_string(), "18446744073709551615");
    assert_eq!(v, IntN::from_str_radix(&v.to_string(), 10).unwrap());
    let v = -v;
    assert_eq!(v.to_string(), "-18446744073709551615");
    assert_eq!(v, IntN::from_str_radix(&v.to_string(), 10).unwrap());
    let v = -v;
    assert_eq!(v.to_string(), "18446744073709551615");
    assert_eq!(v, IntN::from_str_radix(&v.to_string(), 10).unwrap());

    assert_eq!(
        <IntN::<4> as num_traits::Num>::from_str_radix("0", 10).unwrap(),
        0_i64.into()
    );
    assert_eq!(
        <IntN::<4> as num_traits::Num>::from_str_radix("-0", 10).unwrap(),
        0_i64.into()
    );
}

#[test]
fn test_int_big() {
    let a = prime_25519::<IntN<8>>();
    assert_eq!(
        a.to_string(),
        "57896044618658097711785492504343953926634992332820282019728792003956564819949"
    );
    test_number_n(a);
    assert_eq!(IntN::<8>::ONE / a, IntN::<8>::ZERO);

    let a = prime_25519::<IntN<5>>();
    assert_eq!(
        a.to_string(),
        "57896044618658097711785492504343953926634992332820282019728792003956564819949"
    );
    test_number_n(a);
    let a_as_string = a.to_string();
    let a_from_string = <IntN<5> as num_traits::Num>::from_str_radix(&a_as_string, 10).unwrap();
    assert_eq!(a, a_from_string);
}

#[test]
#[should_panic]
fn test_int_add_overflow() {
    let a: IntN<1> = u64::MAX.into();
    let one: IntN<_> = 1_i64.into();
    let _b = a + one;
}

#[test]
#[should_panic]
fn test_int_mul_overflow_from_last() {
    let a: IntN<1> = u64::MAX.into();
    let two: IntN<_> = 2_i64.into();
    let _b = two * a;
}

#[test]
#[should_panic]
fn test_int_mul_overflow_from_middle() {
    let a: IntN<4> = u64::MAX.into();
    let a2 = a * a;
    let a3 = a * a2;
    let _b = a2 * a3;
}

#[test]
#[should_panic]
fn test_int_div_zero() {
    let a = prime_25519::<IntN<8>>();
    let _b = a / IntN::ZERO;
}

#[test]
#[should_panic]
fn test_int_div_assign_zero() {
    let mut a = prime_25519::<IntN<8>>();
    a /= IntN::ZERO;
}

#[test]
#[should_panic]
fn test_int_rem_zero() {
    let a = prime_25519::<IntN<8>>();
    let _b = a % IntN::ZERO;
}

#[test]
#[should_panic]
fn test_int_rem_assign_zero() {
    let mut a = prime_25519::<IntN<8>>();
    a %= IntN::ZERO;
}
