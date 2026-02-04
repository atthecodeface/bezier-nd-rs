//a Imports
use bezier_nd::Bezier;
use bezier_nd::BezierEval;
use bezier_nd::{BezierIntoIterator, BezierSplit};
use bezier_nd::Float;
mod utils;
use geo_nd::{vector, FArray, Vector};

fn bezier_straight_as<F: Float>(bezier: &Bezier<F, 2>, straightness: F)
where
    FArray<F, 2>: Vector<F, 2>,
{
    for i in 0..30 {
        let s: F = (1.4_f32).powi(i - 15).into();
        eprintln!("{} {} {}", straightness, s, bezier.is_straight(s));
        assert_eq!(
            straightness < s,
            bezier.is_straight(s),
            "Bezier {bezier} .is_straight({s}) mismatched with straightness<s for {straightness}",
        );
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
    B: BezierEval<F, [F; D]> + BezierSplit + Clone + BezierIntoIterator<F, B, [F; D]>,
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
        for t in utils::float_iter(F::ZERO, F::ONE, NSEG_PTS) {
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
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [10., 1.].into();
    let p3: FArray<f32, 2> = [20., 0.].into();
    let p4: FArray<f32, 2> = [20., 1.].into();

    let sp0 = p0 * 10.;
    let sp1 = p1 * 10.;
    let sp2 = p2 * 10.;
    let sp3 = p3 * 10.;
    let sp4 = p4 * 10.;

    let mut b = Bezier::line(&p0, &p1);
    let sb = Bezier::line(&sp0, &sp1);
    b.scale(10.);
    utils::vec_eq(b.control_point(0), sb.control_point(0));
    utils::vec_eq(b.control_point(1), sb.control_point(1));

    bezier_straight_as(&Bezier::line(&p0, &p1), 1E-10);
    bezier_straight_as(&Bezier::line(&p0, &p2), 1E-10);
    bezier_straight_as(&Bezier::line(&p0, &p3), 1E-10);
    bezier_straight_as(&Bezier::line(&p0, &p4), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp1), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp2), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp3), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp4), 1E-10);

    // P0, P1, P3 are in a line so should be perfectly straight
    bezier_straight_as(&Bezier::quadratic(&p0, &p1, &p3), 1E-10);
    bezier_straight_as(&Bezier::quadratic(&sp0, &sp1, &sp3), 1E-10);

    // P0 -> P3 with P2 as a control; P2 is a small amount above the centre of P3
    // This basically tests scaling of straightness is linear, and it matches
    // these know good values of straightness
    //
    // (The values for straightness here are determined by hand...)
    bezier_straight_as(&Bezier::quadratic(&p0, &p2, &p3), 0.8);
    bezier_straight_as(&Bezier::quadratic(&sp0, &sp2, &sp3), 8.0);

    bezier_straight_as(&Bezier::quadratic(&p0, &p1, &p4), 0.5);
    bezier_straight_as(&Bezier::quadratic(&sp0, &sp1, &sp4), 5.0);

    bezier_straight_as(&Bezier::cubic(&p0, &p1, &p2, &p3), 0.8);
    bezier_straight_as(&Bezier::cubic(&sp0, &sp1, &sp2, &sp3), 8.0);

    let mut b = Bezier::cubic(&p0, &p1, &p2, &p4);
    let sb = Bezier::cubic(&sp0, &sp1, &sp2, &sp4);
    b.scale(10.);
    utils::vec_eq(b.control_point(0), sb.control_point(0));
    utils::vec_eq(b.control_point(1), sb.control_point(1));
    utils::vec_eq(b.control_point(2), sb.control_point(2));
    utils::vec_eq(b.control_point(3), sb.control_point(3));

    bezier_straight_as(&Bezier::cubic(&p0, &p1, &p2, &p4), 0.6);
    bezier_straight_as(&Bezier::cubic(&sp0, &sp1, &sp2, &sp4), 6.0);
}

//fi test_within_straightness
fn test_within_straightness() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [10., 1.];
    let p3 = [20., 0.];

    for straightness_sq in [0.1, 0.01, 0.001] {
        let b = Bezier::line(&p0, &p2);
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p2], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p2], straightness_sq);

        let b = Bezier::quadratic(&p0, &p3, &p1);
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p1, p2], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p1, p2], straightness_sq);

        let b = Bezier::cubic(&p0, &p1, &p2, &p3);
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p1, p2, p3], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p1, p2, p3], straightness_sq);

        let b = Bezier::cubic(&p0, &p2, &p1, &p3);
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p0, p2, p1, p3], straightness_sq);
        bezier_lines_within_straightness(&vec![p0, p2, p1, p3], straightness_sq);

        let b = Bezier::cubic(&p3, &p2, &p1, &p0);
        bezier_lines_within_straightness(&b, straightness_sq);
        bezier_lines_within_straightness(&[p3, p2, p1, p0], straightness_sq);
        bezier_lines_within_straightness(&vec![p3, p2, p1, p0], straightness_sq);

        let b = Bezier::cubic(&p2, &p3, &p0, &p1);
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
