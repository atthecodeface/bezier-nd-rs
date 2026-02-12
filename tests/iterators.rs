//! Test the PointIter, PointTIter, LineIter, LineTIter for Beziers
//!
//! As part of the iterators, this also exercises the splitting into
//! line segments for a Bezier within a closeness_sq
//!
//! Specifically this tests types:
//!
//! * [[f32;N];2] - Linear Bezier
//! * [[f32;N];3] - Quadratic Bezier
//! * [[f32;N];4] - Cubic Bezier

use bezier_nd::BezierEval;
use bezier_nd::BezierFlatIterator;
use bezier_nd::BezierSplit;
use geo_nd::vector;
mod utils;

/// Test the PointIter, PointTIter, LineIter, LineTIter for Beziers
///
/// This will also test BezierSplit, and 'point_at' of BezierEval
fn test_points<
    B: BezierEval<f32, [f32; D]> + BezierSplit + Clone + std::fmt::Debug,
    const D: usize,
>(
    bezier: &B,
    closeness_sq: f32,
) {
    eprintln!(
        "test_points for Bezier of degree {} and dimension {D} with closeness_sq {closeness_sq}",
        bezier.degree()
    );
    eprintln!("Bezier {bezier:?}");
    for (t, p) in bezier.as_t_points(closeness_sq) {
        eprintln!("Next t/p {t} {p:?}");
        let bp = bezier.point_at(t);
        let d_sq = vector::distance_sq(&bp, &p);
        assert!(
            d_sq < 1E-10,
            "Distance between point at {t} {p:?} and actual bezier at t {bp:?} is too big {d_sq}"
        );
    }

    // Create a set of points for 0<=t<=1 on the Bezier
    let bezier_all_pts = utils::BezierPtSet::of_point_at(bezier, 100);
    let max_pt_sep_sq = bezier_all_pts.max_pt_separation_sq();
    eprintln!(
        "Bezier point set with {} pts has max separation_sq of {max_pt_sep_sq}",
        bezier_all_pts.len()
    );

    // Create a set of points where the whole Bezier is within closeness_sq, and check the points are on the Bezier
    let bezier_pts_of_closeness = utils::BezierPtSet::of_point_iter(bezier.as_points(closeness_sq));
    let max_d_sq = bezier_pts_of_closeness.max_distance_sq_of_all_pts(&bezier_all_pts);
    assert!(max_d_sq < max_pt_sep_sq + 1E-10, "Distance squared between bezier as points is too far ({max_d_sq}) from the Bezier (or rather the closest of {} equidistant points on the Bezier", bezier_all_pts.len());

    // Create a set of points that is the approximation as lines interpolated by 100 for each line segment
    let bezier_lines_pts =
        utils::BezierPtSet::of_line_iter_interpolated(bezier_pts_of_closeness.as_lines(), 100);
    let max_interpolated_pt_sep_sq = bezier_lines_pts.max_pt_separation_sq();
    let max_d_sq = bezier_all_pts.max_distance_sq_of_all_pts(&bezier_lines_pts);
    assert!(max_d_sq < closeness_sq + max_interpolated_pt_sep_sq, "Distance squared between bezier and the approximation as lines is too far ({max_d_sq}) from the approximation max {} = closenss_sq {closeness_sq} + error bar {max_interpolated_pt_sep_sq}",
        closeness_sq + max_interpolated_pt_sep_sq
    )
}

fn test_point_iter<const D: usize, const N: usize>(seed: &str)
where
    [[f32; D]; N]: BezierEval<f32, [f32; D]> + BezierSplit,
{
    let mut rng = utils::make_random(seed);
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();
    let bezier = utils::new_random_point_array(&mut rng, &distribution);
    for closeness_sq in [1.0, 0.1, 0.0001] {
        test_points(&bezier, closeness_sq);
    }
}

#[test]
fn test_point_iterators() {
    // Test straight bezier
    test_point_iter::<1, 2>("straight bezier dim 1");
    test_point_iter::<1, 3>("quadratic bezier dim 1");
    test_point_iter::<1, 4>("cubic bezier dim 1");

    test_point_iter::<10, 2>("straight bezier dim 10");
    test_point_iter::<10, 3>("quadratic bezier dim 10");
    test_point_iter::<10, 4>("cubic bezier dim 10");
}
