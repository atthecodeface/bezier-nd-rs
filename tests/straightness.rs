//a Imports
use bezier_nd::Bezier;
use bezier_nd::BezierEval;
use bezier_nd::Float;
use bezier_nd::{BezierIntoIterator, BezierOps, BezierSplit};
mod utils;
use geo_nd::{vector, FArray, Vector};

#[track_caller]
fn bezier_straight_as<F: Float>(bezier: &Bezier<F, 2>, straightness: F)
where
    FArray<F, 2>: Vector<F, 2>,
{
    for i in 0..30 {
        let s: F = (1.4_f32).powi(i - 15).into();
        //eprintln!(
        //    "{bezier} is_straight {straightness} {s} {}",
        //    bezier.is_straight(s)
        //);

        if straightness < s {
            assert!(
                bezier.is_straight(s),
                "Bezier {bezier} should be straighter than {straightness}",
            );
        } else {
            assert!(
                !bezier.is_straight(s),
                "Bezier {bezier} should not be straighter than {straightness}",
            );
        }
    }
}

/// Create NPTS points on the bezier, and store them in a Vec
///
/// Split the bezier into lines given the straightness
///
/// For each line segment find NSEG points along the line and store all of these in a Vec
///
/// Find the largest closest distance between points on the bezier to segmented points
///
/// Find the largest closest distance between segmented points to points on the bezier
///
/// Both of these should be less than 'straightness'
fn bezier_lines_within_straightness<F: Float, const D: usize, B>(bezier: &B, straightness_sq: F)
where
    B: BezierEval<F, [F; D]> + BezierSplit + Clone + BezierIntoIterator<F, [F; D]>,
{
    const NPTS: isize = 1000;
    const NSEG_PTS: usize = 1000;
    let bezier_pts: Vec<[F; D]> = (0..NPTS)
        .map(|i: isize| {
            let t: F = ((i as f32) / (NPTS as f32)).into();
            bezier.point_at(t)
        })
        .collect();
    let mut segment_pts = Vec::<[F; D]>::new();
    for (p0, p1) in bezier.as_lines(straightness_sq) {
        for t in utils::float_iter(NSEG_PTS) {
            segment_pts.push(vector::sum_scaled(&[p0, p1], &[F::ONE - t, t]));
        }
    }
    utils::assert_pts_all_within_closeness_sq_of_pts(
        &segment_pts,
        bezier_pts.iter(),
        straightness_sq,
    );
    utils::assert_pts_all_within_closeness_sq_of_pts(
        &bezier_pts,
        segment_pts.iter(),
        straightness_sq,
    );
}

//fi test_straight_as
fn test_straight_as() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [10., 1.];
    let p3 = [20., 0.];
    let p4 = [20., 1.];

    let sp0 = vector::scale(p0, 10.0);
    let sp1 = vector::scale(p1, 10.0);
    let sp2 = vector::scale(p2, 10.0);
    let sp3 = vector::scale(p3, 10.0);
    let sp4 = vector::scale(p4, 10.0);

    let mut b: Bezier<_, _> = [p0, p1].into();
    let sb: Bezier<_, _> = [sp0, sp1].into();
    b.scale(10.);
    utils::bezier_eq(&b, sb.control_points());
    utils::bezier_eq(&sb, b.control_points());

    bezier_straight_as(&[p0, p1].into(), 1E-10);
    bezier_straight_as(&[p0, p2].into(), 1E-10);
    bezier_straight_as(&[p0, p3].into(), 1E-10);
    bezier_straight_as(&[p0, p4].into(), 1E-10);
    bezier_straight_as(&[sp0, sp1].into(), 1E-10);
    bezier_straight_as(&[sp0, sp2].into(), 1E-10);
    bezier_straight_as(&[sp0, sp3].into(), 1E-10);
    bezier_straight_as(&[sp0, sp4].into(), 1E-10);

    // P0, P1, P3 are in a line so should be perfectly straight
    bezier_straight_as(&[p0, p1, p3].into(), 1E-10);
    bezier_straight_as(&[sp0, sp1, sp3].into(), 1E-10);

    // P0 -> P3 with P2 as a control; P2 is a small amount above the centre of P3
    // This basically tests scaling of straightness is linear, and it matches
    // these know good values of straightness
    //
    // (The values for straightness here are determined by hand...)
    bezier_straight_as(&[p0, p2, p3].into(), 0.8);
    bezier_straight_as(&[sp0, sp2, sp3].into(), 8.0);

    bezier_straight_as(&[p0, p1, p4].into(), 0.5);
    bezier_straight_as(&[sp0, sp1, sp4].into(), 5.0);

    bezier_straight_as(&[p0, p1, p2, p3].into(), 0.8);
    bezier_straight_as(&[sp0, sp1, sp2, sp3].into(), 8.0);

    let mut b = [p0, p1, p2, p4];
    let sb = [sp0, sp1, sp2, sp4];

    b.scale(10.);
    utils::bezier_eq(&b, sb.control_points());
    utils::bezier_eq(&sb, b.control_points());

    bezier_straight_as(&[p0, p1, p2, p4].into(), 0.6);
    bezier_straight_as(&[sp0, sp1, sp2, sp4].into(), 6.0);
}

//fi test_within_straightness
fn test_within_straightness() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [10., 1.];
    let p3 = [20., 0.];

    for straightness_sq in [0.1, 0.01, 0.001] {
        let b: Bezier<_, _> = [p0, p2].into();
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p2], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p2], straightness_sq);

        let b: Bezier<_, _> = [p0, p3, p1].into();
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p1, p2], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p1, p2], straightness_sq);

        let b: Bezier<_, _> = [p0, p1, p2, p3].into();
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p1, p2, p3], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p1, p2, p3], straightness_sq);

        let b = [p0, p2, p1, p3];
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p2, p1, p3], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p2, p1, p3], straightness_sq);

        let b = vec![p3, p2, p1, p0];
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p3, p2, p1, p0], straightness_sq);
        bezier_lines_within_straightness(&vec![p3, p2, p1, p0], straightness_sq);

        let b: Bezier<_, _> = [p2, p3, p0, p1].into();
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p2, p3, p0, p1], straightness_sq);
        bezier_lines_within_straightness(&vec![p2, p3, p0, p1], straightness_sq);
    }
}

//a Tests
#[test]
fn test_f32_within_straightness() {
    test_within_straightness();
}

#[test]
fn test_f32_straight_as() {
    test_straight_as();
}
