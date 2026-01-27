//! Test the Approximation type
//!
//! As part of the iterators, this also exercises the splitting into
//! line segments for a Bezier within a closeness_sq
//!
//! Specifically this tests types:
//!
//! * [[f32;N];2] - Linear Bezier
//! * [[f32;N];3] - Quadratic Bezier
//! * [[f32;N];4] - Cubic Bezier

//a Imports
use bezier_nd::Approximation;
use bezier_nd::BezierEval;
use bezier_nd::BezierSplit;
mod utils;

fn test_approximation<
    B: BezierEval<f32, [f32; D]> + BezierSplit + Clone + std::fmt::Debug,
    const D: usize,
>(
    bezier: &B,
    closeness_sq: f32,
) {
    let a = Approximation::new(bezier, closeness_sq);
    let max_dist_sq = utils::max_distance_sq(bezier, &a, 1000);
    assert!(
        max_dist_sq <= closeness_sq,
        "Approximation not close enough to original Bezier (max {max_dist_sq}, closeness_sq {closeness_sq})"
    );
}

#[test]
fn t0() {
    let mut rng = utils::make_random("banana");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();
    let mut bezier = [[0.0_f32; 1]; 3];
    utils::random_point_array(&mut rng, &distribution, &mut bezier);
    for closeness_sq in [1.0, 0.3, 0.1, 0.01, 0.001, 0.0001] {
        test_approximation(&bezier, closeness_sq);
    }
}
