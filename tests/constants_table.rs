use bezier_nd::Num;
use bezier_nd::CONSTANTS_TABLE;

fn assert_eq_value<F: Num>(index: usize, value: F) -> impl Fn(&[F]) -> bool {
    move |table| {
        assert_eq!(
            table[index], value,
            "Expected table index {index} to match value: {table:?}"
        );
        true
    }
}
#[test]
fn test_table_f64() {
    assert!(CONSTANTS_TABLE.use_constants_table(assert_eq_value(0, 5.0), &[5., 6., 7., 8.],));
}

#[test]
fn test_table_f32() {
    assert!(CONSTANTS_TABLE.use_constants_table(assert_eq_value(0, 5.0), &[5., 6., 7., 8.],));
}

/*
#[test]
fn test_table_rational() {
    assert!(CONSTANTS_TABLE.use_constants_table(
        assert_eq_value::<crate::bignum::RationalN<2>>(0, (5, 1).into()),
        &[5., 6., 7., 8.],
    ));
}
*/
