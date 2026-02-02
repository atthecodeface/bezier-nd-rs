mod utils;
use bezier_nd::{bernstein_fns, Bezier};
use geo_nd::FArray;
use utils::test_beziers_approx_eq;

#[test]
fn elevate() {
    let p0: FArray<f64, 2> = [0., 0.].into();
    let p1: FArray<f64, 2> = [10., 0.].into();
    let p2: FArray<f64, 2> = [6., 1.].into();
    let p3: FArray<f64, 2> = [20., 5.].into();

    let b = Bezier::quadratic(&p0, &p1, &p2);
    let b2 = b.clone().elevate();
    eprintln!("Elevated {b} to {b2}");
    test_beziers_approx_eq(&b, &b2);

    let b = Bezier::quadratic(&p3, &p1, &p2);
    let b2 = b.clone().elevate();
    eprintln!("Elevated {b} to {b2}");
    test_beziers_approx_eq(&b, &b2);

    let b = Bezier::line(&p3, &p1);
    let b2 = b.clone().elevate();
    eprintln!("Elevated {b} to {b2}");
    test_beziers_approx_eq(&b, &b2);
}

#[test]
fn elevate_matrix() {
    let mut m = [0.0f64; 100];

    // Check point to linear matrix
    assert_eq!(
        bernstein_fns::generate_elevate_by_one_matrix(&mut m, 0),
        1.0
    );
    assert_eq!(&m[0..2], &[1., 1.]);

    // Check linear - quadratic matrix
    assert_eq!(
        bernstein_fns::generate_elevate_by_one_matrix(&mut m, 1),
        2.0
    );
    assert_eq!(&m[0..6], &[2., 0., 1.0, 1.0, 0., 2.]);

    // Check cubic to order 4 matrix
    assert_eq!(
        bernstein_fns::generate_elevate_by_one_matrix(&mut m, 3),
        4.0
    );
    assert_eq!(&m[0..4], &[4., 0., 0., 0.]);
    assert_eq!(&m[4..8], &[1.0, 3., 0., 0.]);
    assert_eq!(&m[8..12], &[0., 2., 2., 0.]);
    assert_eq!(&m[12..16], &[0., 0., 3., 1.]);
    assert_eq!(&m[16..20], &[0., 0., 0., 4.]);
}
