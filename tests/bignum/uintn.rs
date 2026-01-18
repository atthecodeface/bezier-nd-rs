use bezier_nd::bignum::{IntN, UIntN};
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

    assert!(b > a);
    assert!(a < b);

    assert!(c > a);
    assert!(a < c);

    assert!(c > b);
    assert!(b < c);
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

fn test_gcd_n<const N: usize>(a: UIntN<N>) {
    assert_eq!(a.gcd(&a), a);
    assert_eq!(a.gcd(&(a + 1.into())), 1.into());
}

fn test_number_n<N: Num + From<u64>>(a: N) {
    test_add_subtract_n(a);
    test_multiply_small_n(a);
    assert!(<N as num_traits::Num>::from_str_radix("", 10).is_err());
    assert!(<N as num_traits::Num>::from_str_radix("x", 10).is_err());
}

#[test]
fn test_uint_tiny() {
    assert_eq!(UIntN::<1>::ZERO.to_string(), "0");

    test_number_n::<UIntN<1>>(1.into());
    test_number_n::<UIntN<1>>((u64::MAX >> 2).into());

    test_number_n::<UIntN<8>>(1.into());
    test_number_n::<UIntN<8>>((u64::MAX >> 2).into());

    test_gcd_n::<1>(1.into());
    let a: UIntN<8> = 369.into();
    let b: UIntN<8> = 594.into();
    let zero = UIntN::<8>::ZERO;
    assert_eq!(a.gcd(&b), 9.into());
    assert_eq!(zero.gcd(&a), zero);
    assert_eq!(a.gcd(&zero), zero);

    let v: UIntN<1> = u64::MAX.into();
    assert_eq!(v.to_string(), "18446744073709551615");
}

#[test]
fn test_uint_big() {
    let a = prime_25519::<UIntN<8>>();
    assert_eq!(
        a.to_string(),
        "57896044618658097711785492504343953926634992332820282019728792003956564819949"
    );
    test_number_n(a);
    assert_eq!(UIntN::<8>::ONE / a, UIntN::<8>::ZERO);

    let a = prime_25519::<UIntN<5>>();
    assert_eq!(
        a.to_string(),
        "57896044618658097711785492504343953926634992332820282019728792003956564819949"
    );
    test_number_n(a);
    let a_as_string = a.to_string();
    let a_from_string = <UIntN<5> as num_traits::Num>::from_str_radix(&a_as_string, 10).unwrap();
    assert_eq!(a, a_from_string);

    assert_eq!(a.gcd(&a), a);
    assert_eq!(a.gcd(&((u64::MAX).into())), 1_u64.into());
}

#[test]
#[should_panic]
fn test_uint_neg() {
    let a = prime_25519::<UIntN<8>>();
    let _b = -a;
}

#[test]
#[should_panic]
fn test_uint_add_overflow() {
    let a: UIntN<1> = u64::MAX.into();
    let one: UIntN<_> = 1.into();
    let _b = a + one;
}

#[test]
#[should_panic]
fn test_uint_mul_overflow_from_last() {
    let a: UIntN<1> = u64::MAX.into();
    let two: UIntN<_> = 2.into();
    let _b = two * a;
}

#[test]
#[should_panic]
fn test_uint_mul_overflow_from_middle() {
    let a: UIntN<4> = u64::MAX.into();
    let a2 = a * a;
    let a3 = a * a2;
    let _b = a2 * a3;
}

#[test]
#[should_panic]
fn test_uint_div_zero() {
    let a = prime_25519::<UIntN<8>>();
    let _b = a / UIntN::ZERO;
}

#[test]
#[should_panic]
fn test_uint_div_assign_zero() {
    let mut a = prime_25519::<UIntN<8>>();
    a /= UIntN::ZERO;
}

#[test]
#[should_panic]
fn test_uint_rem_zero() {
    let a = prime_25519::<UIntN<8>>();
    let _b = a % UIntN::ZERO;
}

#[test]
#[should_panic]
fn test_uint_rem_assign_zero() {
    let mut a = prime_25519::<UIntN<8>>();
    a %= UIntN::ZERO;
}
