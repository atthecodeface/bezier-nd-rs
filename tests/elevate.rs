mod utils;
use bezier_nd::{bernstein_fns, BezierElevate, BezierND};
use utils::test_beziers_approx_eq;

#[test]
fn elevate_bezier() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [6., 1.];
    let p3 = [20., 5.];

    let b = BezierND::<_, 4, _>::new(&[p0, p1, p2]);
    let b2 = b.elevate_by_one().unwrap();
    eprintln!("Elevated {b} to {b2}");
    test_beziers_approx_eq(&b, &b2);

    let b = BezierND::<_, 4, _>::new(&[p3, p1, p2]);
    let b2 = b.elevate_by_one().unwrap();
    eprintln!("Elevated {b} to {b2}");
    test_beziers_approx_eq(&b, &b2);

    let b = BezierND::<_, 4, _>::new(&[p3, p1]);
    let b2 = b.elevate_by_one().unwrap();
    eprintln!("Elevated {b} to {b2}");
    test_beziers_approx_eq(&b, &b2);
}

#[test]
fn elevate_farray() {
    let p0 = [0.0_f32, 0.];
    let p1 = [10., 0.];
    let p2 = [6., 1.];
    let p3 = [20., 5.];

    let b = [p0, p1, p2];
    let b2 = b.elevate_by_one().unwrap();
    eprintln!("Elevated {b:?} to {b2:?}");
    test_beziers_approx_eq(&b, &b2);

    let b = [p3, p1, p2];
    let b2 = b.elevate_by_one().unwrap();
    eprintln!("Elevated {b:?} to {b2:?}");
    test_beziers_approx_eq(&b, &b2);

    let b = [p3, p1];
    let b2 = b.elevate_by_one().unwrap();
    eprintln!("Elevated {b:?} to {b2:?}");
    test_beziers_approx_eq(&b, &b2);
}

#[test]
fn elevate_fvec() {
    let seed = "banana4";
    let mut rng = utils::make_random(seed);
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    let b = utils::new_random_point_vec::<_, f32, _, 1>(&mut rng, &distribution, 2);
    let b1 = b.elevate_by_one().unwrap();
    let b2 = b1.elevate_by_one().unwrap();
    let b2_v2 = b.elevate_by(2).unwrap();

    eprintln!("Compare first elevate-by-1");
    test_beziers_approx_eq(&b, &b1);
    eprintln!("Compare second elevate-by-1");
    test_beziers_approx_eq(&b, &b2);
    eprintln!("Compare elevate-by-2 with two elevate-by-1");
    test_beziers_approx_eq(&b, &b2_v2);
}

#[test]
fn elevate_matrix() {
    // Check point to linear matrix
    let (f, m) = bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix::<f32>(0);
    assert_eq!(f, 1.0);
    assert_eq!(&m[0..2], &[1., 1.]);

    // Check linear - quadratic matrix
    let (f, m) = bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix::<f32>(1);
    assert_eq!(f, 2.0);
    assert_eq!(&m[0..6], &[2., 0., 1.0, 1.0, 0., 2.]);

    // Check cubic to order 4 matrix
    let (f, m) = bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix::<f32>(3);
    assert_eq!(f, 4.0);
    assert_eq!(&m[0..4], &[4., 0., 0., 0.]);
    assert_eq!(&m[4..8], &[1.0, 3., 0., 0.]);
    assert_eq!(&m[8..12], &[0., 2., 2., 0.]);
    assert_eq!(&m[12..16], &[0., 0., 3., 1.]);
    assert_eq!(&m[16..20], &[0., 0., 0., 4.]);
}
