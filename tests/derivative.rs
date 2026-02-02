//a Imports
use bezier_nd::{BezierEval, BezierND};
mod utils;
use utils::assert_near_equal;

#[test]
fn derivative() {
    let p0 = [0., 0.];
    let p1 = [10., -1.];
    let p2 = [6., 1.];
    let p3 = [20., 5.];

    let b = BezierND::<f64, 4, 2>::new(&[p0, p1, p2, p3]);
    assert_eq!(b.degree(), 3);

    let (db, db_n) = b.nth_derivative(1);
    assert_eq!(db.degree(), 2);
    assert_eq!(db_n, 3.0);
    assert_near_equal(&db.pts()[0], &[10., -1.]);
    assert_near_equal(&db.pts()[1], &[-4., 2.]);
    assert_near_equal(&db.pts()[2], &[14., 4.]);

    let (d2b, d2b_n) = b.nth_derivative(2);
    assert_eq!(d2b.degree(), 1);
    assert_eq!(d2b_n, 6.0);
    assert_near_equal(&d2b.pts()[0], &[-14., 3.]);
    assert_near_equal(&d2b.pts()[1], &[18., 2.]);

    let (d3b, d3b_n) = b.nth_derivative(3);
    assert_eq!(d3b.degree(), 0);
    assert_eq!(d3b_n, 6.0);
    assert_near_equal(&d3b.pts()[0], &[32., -1.]);
}
