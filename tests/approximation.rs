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
use bezier_nd::BezierEval;
use bezier_nd::BezierND;
use bezier_nd::BezierSplit;
mod utils;
use geo_nd::vector;

fn test_approximation<
    B: BezierEval<f32, [f32; D]> + BezierSplit + Clone + std::fmt::Debug,
    const D: usize,
>(
    bezier: &B,
    closeness_sq: f32,
) {
    eprintln!("Testing bexier approximation {bezier:?} {closeness_sq}");
    let a = Approximation::new(bezier, closeness_sq);
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
}

#[test]
fn test_f_array() {
    let mut rng = utils::make_random("banana");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    // Test with a quadratic Bezier with 10 diumensional points
    let bezier = utils::new_random_point_array::<_, _, _, 3, 10>(&mut rng, &distribution);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation(&bezier, closeness_sq);
    }

    // Test with a cubic Bezier with 6 diumensional points
    let bezier = utils::new_random_point_array::<_, _, _, 4, 6>(&mut rng, &distribution);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation(&bezier, closeness_sq);
    }
}

#[test]
fn test_fvec() {
    let mut rng = utils::make_random("banana_split");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    for degree in 1..10 {
        let mut bezier = vec![[0.0_f32; 2]; (degree + 1)];
        utils::set_random_point_array(&mut rng, &distribution, &mut bezier);
        for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
            test_approximation(&bezier, closeness_sq);
        }
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
        test_approximation(&bezier_nd, closeness_sq);
    }

    // Test with a degree 7 Bezier with 10 diumensional points
    let bezier = utils::new_random_point_array::<_, _, _, 8, 10>(&mut rng, &distribution);
    let bezier_nd = BezierND::<_, 9, _>::new(&bezier);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation(&bezier_nd, closeness_sq);
    }
}
