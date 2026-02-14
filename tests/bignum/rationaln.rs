use bezier_nd::bignum::RationalN;
use num_traits::{one, zero, ConstZero, Num, Zero};

fn prime_25519<N: Num + Copy + From<u64>>() -> N {
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

fn assert_bits_almost_eq<const N: usize>(a: RationalN<N>, v: f64) {
    let a_bits = a.as_f64_bits().unwrap();
    let v_bits = v.to_bits();

    assert_eq!(
        a_bits & !1,
        v_bits & !1,
        "Mismatch for {a} {a_bits} to {v_bits}",
    );

    let a_bits = a.as_f32_bits().unwrap();
    let v_bits = (v as f32).to_bits();
    assert_eq!(
        a_bits & !1,
        v_bits & !1,
        "Mismatch for {a} {v} : {a_bits:x} to {v_bits:x}",
    );

    let a_f64: f64 = a.into();
    let a_f32: f32 = a.into();
    assert!((v - a_f64).abs() < 1E-12);
    assert!(((v as f32) - a_f32).abs() < 1E-6);
}

fn test_mantissa<const N: usize>(a: RationalN<N>, m: u64, exp: i32, v: f64, precise: bool) {
    assert_eq!(a.to_mantissa64_exp(), (m, exp));
    assert_eq!((-a).to_mantissa64_exp(), (m, exp));
    assert_bits_almost_eq(a, v);

    if precise {
        assert_eq!(a, v.try_into().unwrap());
        assert_eq!(a, (v as f32).try_into().unwrap());
    }

    // Note that -0.0 !+ 0.0 but 0/1 == -0/1
    if !a.is_zero() {
        assert_bits_almost_eq(-a, -v);
        if precise {
            assert_eq!(-v, (-a).into());
            assert_eq!((-v) as f32, (-a).into());
        }
    }
}

fn mantissa_exp<const N: usize>() {
    test_mantissa::<N>(0_u64.into(), 0, 0, 0.0, true);
    test_mantissa::<N>(1_i64.into(), 1 << 63, 0, 1.0, true);
    test_mantissa::<N>(2_i64.into(), 1 << 63, 1, 2.0, true);
    test_mantissa::<N>((1_i64, 2_u64).into(), 1 << 63, -1, 0.5, true);
    test_mantissa::<N>((1_i64, 4_u64).into(), 1 << 63, -2, 0.25, true);
    test_mantissa::<N>((3_i64, 4_u64).into(), 3 << 62, -1, 0.75, true);
    test_mantissa::<N>((3_i64, 1_u64).into(), 3 << 62, 1, 3.0, true);
    test_mantissa::<N>(
        (1_i64, 3_u64).into(),
        0xaaaaaaaaaaaaaaaa,
        -2,
        1.0 / 3.0,
        false,
    );
}

#[test]
fn test_mantissa_exp() {
    mantissa_exp::<1>();
    mantissa_exp::<4>();
    mantissa_exp::<16>();

    let f = 1.0_f64 / 12.0;
    let value: RationalN<8> = f.try_into().unwrap();
    eprintln!("Value {value}");
    let value_f: f64 = value.into();
    assert!((f - value_f).abs() < 1E-8);
}

fn assert_approx_eq<const N: usize>(a: &RationalN<N>, scale: &RationalN<N>, value: f64) {
    let a_f64: f64 = a.into();
    let scale_f64: f64 = scale.into();
    let delta = (a_f64 * scale_f64 - value).abs();
    assert!(
        delta < 1E-8,
        "{a}*{scale} = {a_f64}*{scale_f64} should equal {value}",
    );
}

fn test_arithmetic<const N: usize>(value: f64, scale: f64) {
    let mut v_r: RationalN<N> = value.try_into().unwrap();
    let s_r: RationalN<N> = scale.try_into().unwrap();
    v_r /= s_r;
    assert_approx_eq(&v_r, &s_r, value);
    assert_approx_eq(&(v_r + v_r), &s_r, value * 2.0);
    //    assert_approx_eq(&(v_r * v_r), &s_r, value * value / scale);
    //    assert_approx_eq(&((v_r + v_r) / v_r), &s_r, 2.0 * scale);
}

#[test]
fn test_arithmetic_one() {
    for scale in [1.0, 2., 3., 4., 8.] {
        test_arithmetic::<1>(0.25, scale);
        test_arithmetic::<1>(0.5, scale);
        test_arithmetic::<1>(1000., scale);
    }
}

#[test]
fn test_arithmetic_two() {
    for scale in [1.0, 2., 3., 4., 8.] {
        test_arithmetic::<2>(0.25, scale);
        test_arithmetic::<2>(0.5, scale);
        test_arithmetic::<2>(1000., scale);
        test_arithmetic::<2>(0.1, scale);
        test_arithmetic::<2>(0.01, scale);
    }
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
