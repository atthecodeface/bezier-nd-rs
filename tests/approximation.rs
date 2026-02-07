//! Test the Approximation type
//!
//! The Approximation type is a list of points/line segments
//! and associated values of the 't' parameter such that all
//! the points on the approximation are within a closeness_sq of
//! the point on the Bezier with the same 't' parameter value
//!
//! Approximation can be made from any Bezier that has a point
//! of type [F;D] and a parameter of type F, and that provides
//! BezierEval and BezierSplit.

//a Imports
use bezier_nd::Approximation;
use bezier_nd::Bezier;
use bezier_nd::BezierEval;
use bezier_nd::BezierND;
use bezier_nd::BezierSplit;
mod utils;
use geo_nd::vector;
use utils::test_beziers_approx_eq;

fn test_approximation_closeness_sq<
    B: BezierEval<f32, [f32; D]> + BezierSplit + Clone + std::fmt::Debug,
    const D: usize,
>(
    bezier: &B,
    closeness_sq: f32,
) {
    eprintln!("Testing bexier approximation {bezier:?} {closeness_sq}");
    let a = Approximation::new(bezier, closeness_sq);

    assert_eq!(a.num_control_points(), bezier.num_control_points());
    assert_eq!(a.degree(), bezier.degree());
    for i in 0..a.num_control_points() {
        assert_eq!(a.control_point(i), bezier.control_point(i));
    }
    assert!(a.closeness_sq_to_line() > bezier.closeness_sq_to_line());
    assert!(
        a.closeness_sq_to_line()
            <= (bezier.closeness_sq_to_line().sqrt() + closeness_sq.sqrt()).powi(2)
    );
    assert!(a.dc_sq_from_line() > bezier.dc_sq_from_line());
    assert!(a.dc_sq_from_line() <= (bezier.dc_sq_from_line().sqrt() + closeness_sq.sqrt()).powi(2));

    eprintln!("{a:?}");
    let max_dist_sq = utils::max_distance_sq(bezier, &a, 1000);
    assert!(
        max_dist_sq <= closeness_sq,
        "Approximation not close enough to original Bezier (max {max_dist_sq}, closeness_sq {closeness_sq})"
    );

    for (i, (t, p)) in a.iter_t_pts().enumerate() {
        assert!(
            vector::distance_sq(&bezier.point_at(t), &p) < 1E-6,
            "Approximation should have same points for values of t"
        );
        assert_eq!(a.points()[i], p);
        assert_eq!(a.ts()[i], t);
    }

    for (i, (p0, p1)) in a.iter_lines().enumerate() {
        assert_eq!(a.points()[i], *p0);
        assert_eq!(a.points()[i + 1], *p1);
    }

    let mut a_pt = *a.endpoints().0;
    let mut b_pt = a_pt;
    let steps = 1000;
    let dt = 1.0 / ((steps + 1) as f32);
    for t in utils::float_iter(steps) {
        let a_dt = a.derivative_at(t);
        let b_dt = bezier.derivative_at(t);
        a_pt = vector::add(a_pt, &a_dt.1, a_dt.0 * dt);
        b_pt = vector::add(b_pt, &b_dt.1, b_dt.0 * dt);
    }
    let d_sq = vector::distance_sq(&a_pt, &b_pt);
    assert!(d_sq < (steps as f32)*closeness_sq, "Expected Sum(derivatives) for closeness {closeness_sq} to be about equal but was {d_sq} distance sq apart {a_pt:?} <> {b_pt:?}");
}

#[test]
fn test_f_array() {
    let mut rng = utils::make_random("banana");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    // Test with a quadratic Bezier with 10 diumensional points
    let bezier = utils::new_random_point_array::<_, _, _, 3, 10>(&mut rng, &distribution);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation_closeness_sq(&bezier, closeness_sq);
    }

    // Test with a cubic Bezier with 6 diumensional points
    let bezier = utils::new_random_point_array::<_, _, _, 4, 6>(&mut rng, &distribution);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation_closeness_sq(&bezier, closeness_sq);
    }
}

#[test]
fn test_fvec() {
    let mut rng = utils::make_random("banana_split");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    for degree in 1..39 {
        let mut bezier = vec![[0.0_f32; 2]; degree + 1];
        utils::set_random_point_array(&mut rng, &distribution, &mut bezier);
        for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
            test_approximation_closeness_sq(&bezier, closeness_sq);
        }
    }
}

#[test]
fn test_bezier() {
    let mut rng = utils::make_random("banana");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    // Test with a linear Bezier with 10 dimensional points
    let bezier = utils::new_random_point_array::<_, _, _, 2, 10>(&mut rng, &distribution);
    let bezier_nd: Bezier<_, _> = bezier.into();
    let tf_bezier: [[f32; 10]; 2] = (&bezier_nd).try_into().unwrap();
    let opt_tf: Result<[[f32; 10]; 3], _> = (&bezier_nd).try_into();
    assert!(opt_tf.is_err());
    test_beziers_approx_eq(&bezier_nd, &tf_bezier);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation_closeness_sq(&bezier_nd, closeness_sq);
    }

    // Test with a quadratic Bezier with 10 dimensional points
    let bezier = utils::new_random_point_array::<_, _, _, 3, 10>(&mut rng, &distribution);
    let bezier_nd: Bezier<_, _> = bezier.into();
    let tf_bezier: [[f32; 10]; 3] = (&bezier_nd).try_into().unwrap();
    let opt_tf: Result<[[f32; 10]; 2], _> = (&bezier_nd).try_into();
    assert!(opt_tf.is_err());
    test_beziers_approx_eq(&bezier_nd, &tf_bezier);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation_closeness_sq(&bezier_nd, closeness_sq);
    }

    // Test with a degree 7 Bezier with 10 dimensional points
    let bezier = utils::new_random_point_array::<_, _, _, 4, 10>(&mut rng, &distribution);
    let bezier_nd: Bezier<_, _> = bezier.into();
    let tf_bezier: [[f32; 10]; 4] = (&bezier_nd).try_into().unwrap();
    let opt_tf: Result<[[f32; 10]; 2], _> = (&bezier_nd).try_into();
    assert!(opt_tf.is_err());
    test_beziers_approx_eq(&bezier_nd, &tf_bezier);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation_closeness_sq(&bezier_nd, closeness_sq);
    }
}

#[test]
fn test_bezier_nd() {
    let mut rng = utils::make_random("banana");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    // Test with a quadratic Bezier with 10 diumensional points
    let bezier = utils::new_random_point_array::<_, _, _, 3, 10>(&mut rng, &distribution);
    let bezier_nd = BezierND::<_, 6, _>::new(&bezier);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation_closeness_sq(&bezier_nd, closeness_sq);
    }

    // Test with a degree 7 Bezier with 10 diumensional points
    let bezier = utils::new_random_point_array::<_, _, _, 8, 10>(&mut rng, &distribution);
    let bezier_nd = BezierND::<_, 9, _>::new(&bezier);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation_closeness_sq(&bezier_nd, closeness_sq);
    }
}
