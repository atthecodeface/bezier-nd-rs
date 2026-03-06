use super::super::Fixed;
use super::constants;
use crate::fixed_point::IntN;
#[test]
fn basic() {
    constants::simple_constants::<Fixed<i8, 4>, i8, 4>();
    constants::simple_constants::<Fixed<i8, 2>, i8, 2>();
    constants::simple_constants::<Fixed<i8, 6>, i8, 6>();

    constants::float_constants::<Fixed<i8, 4>, i8, 4>();

    constants::float_constants::<Fixed<i64, 60>, i64, 60>();
    constants::float_constants::<Fixed<IntN<4>, 240>, IntN<4>, 240>();

    //    assert!(false, "Forced failure");
}
