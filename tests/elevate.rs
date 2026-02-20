mod utils;
use bezier_nd::{bernstein_fns, constants, BezierElevate, BezierMap, BezierND};
use geo_nd::matrix;
use utils::test_beziers_approx_eq;

use crate::utils::assert_near_identity;

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

    // Generate a linear Bezier and elevate it twice
    let b = utils::new_random_point_vec::<_, f32, _, 1>(&mut rng, &distribution, 2);
    let b1 = b.elevate_by_one().unwrap();
    let b2 = b1.elevate_by_one().unwrap();

    eprintln!("Compare first elevate-by-1");
    test_beziers_approx_eq(&b, &b1);
    eprintln!("Compare second elevate-by-1");
    test_beziers_approx_eq(&b, &b2);

    // e_deg_3 is 5 by 4 matrix; etc
    //
    // If they are multiplied togehther it is 7 by 4 matrix
    let e_deg_3 = constants::elevate_table_of_match(3).unwrap();
    let e_deg_4 = constants::elevate_table_of_match(4).unwrap();
    let e_deg_5 = constants::elevate_table_of_match(5).unwrap();

    let r_deg_4 =
        constants::reduce_table_of_match(bezier_nd::BezierReduction::LeastSquares, 4).unwrap();
    let r_deg_5 =
        constants::reduce_table_of_match(bezier_nd::BezierReduction::LeastSquares, 5).unwrap();
    let r_deg_6 =
        constants::reduce_table_of_match(bezier_nd::BezierReduction::LeastSquares, 6).unwrap();

    let mut e_deg_35 = vec![0_f64; 24];
    matrix::multiply_dyn(6, 5, 4, e_deg_4, e_deg_3, &mut e_deg_35);
    let mut e_deg_36 = vec![0_f64; 28];
    matrix::multiply_dyn(7, 6, 4, e_deg_5, &e_deg_35, &mut e_deg_36);

    let mut r_deg_64 = vec![0_f64; 35];
    matrix::multiply_dyn(5, 6, 7, r_deg_5, r_deg_6, &mut r_deg_64);
    let mut r_deg_63 = vec![0_f64; 28];
    matrix::multiply_dyn(4, 5, 7, r_deg_4, &r_deg_64, &mut r_deg_63);

    utils::display::eprintln_matrix(&e_deg_36, 4);
    utils::display::eprintln_matrix(&r_deg_63, 7);

    let mut test = vec![0_f64; 16];
    matrix::multiply_dyn(4, 7, 4, &r_deg_63, &e_deg_36, &mut test);
    utils::display::eprintln_matrix(&test, 4);
    assert_near_identity(4, &test);

    let b3 = utils::new_random_point_vec::<_, f64, _, 1>(&mut rng, &distribution, 4);
    let b6 = b3
        .mapped_to_degree(6, &e_deg_36)
        .expect("Must be able to map fvec bezier");
    let b3_v2 = b6
        .mapped_to_degree(3, &r_deg_63)
        .expect("Must be able to map fvec bezier");
    eprintln!("{b3:?}");
    eprintln!("{b6:?}");
    eprintln!("{b3_v2:?}");
    utils::test_beziers_approx_eq(&b3, &b6);
    utils::test_beziers_approx_eq(&b3_v2, &b6);
    utils::test_beziers_approx_eq(&b3_v2, &b3);
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
