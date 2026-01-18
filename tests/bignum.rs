use bezier_nd::bignum::{IntN, RationalN, UIntN};
use geo_nd::Num;

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
}

fn test_multiply_small_n<N: Num + From<u64>>(a: N) {
    let two: N = 2.into();
    let three: N = 3.into();
    assert_eq!(a + a, a * two);
    assert_eq!(a + a + a, a * three);
    assert_eq!(a + a, two * a);
    assert_eq!(a + a + a, three * a);
    assert_eq!(three * a - a, two * a);
}

fn test_number_n<N: Num + From<u64>>(a: N) {
    test_add_subtract_n(a);
    test_multiply_small_n(a);
}

#[test]
fn test_uint_tiny() {
    test_number_n::<UIntN<1>>(1.into());
    test_number_n::<UIntN<1>>((u64::MAX >> 2).into());

    test_number_n::<UIntN<8>>(1.into());
    test_number_n::<UIntN<8>>((u64::MAX >> 2).into());
}

#[test]
fn test_int_tiny() {
    test_number_n::<IntN<1>>(1_u64.into());
    test_number_n::<IntN<1>>((u64::MAX >> 2).into());

    test_number_n::<IntN<8>>(1_u64.into());
    test_number_n::<IntN<8>>((u64::MAX >> 2).into());
}

#[test]
fn test_uint_big() {
    let a = prime_25519::<UIntN<8>>();
    assert_eq!(
        a.to_string(),
        "57896044618658097711785492504343953926634992332820282019728792003956564819949"
    );
    test_number_n(a);

    let a = prime_25519::<UIntN<5>>();
    assert_eq!(
        a.to_string(),
        "57896044618658097711785492504343953926634992332820282019728792003956564819949"
    );
    test_number_n(a);
}
