use bezier_nd::{bernstein_fns, BezierEval, BezierSplit};

use bezier_nd::Num;
use geo_nd::{matrix, vector};

#[track_caller]
pub fn bezier_eq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bez: &B, v: &[[F; D]]) {
    let c = bez.control_points();
    assert_eq!(
        c.len(),
        v.len(),
        "Expected same number of control points for beziers to be equal"
    );
    for (p, v) in c.iter().zip(v.iter()) {
        super::vec_eq(p, v);
    }
}

#[track_caller]
/// Find the maximum distance between `b0` for t0<=t<=t1 and `b1` for 0<=t<=1 in 'steps' intervals
pub fn test_subsection<
    F: Num + From<f32>,
    const D: usize,
    B0: BezierEval<F, [F; D]> + std::fmt::Debug,
    B1: BezierEval<F, [F; D]> + std::fmt::Debug,
>(
    b0: &B0,
    b1: &B1,
    t0: f32,
    t1: f32,
    steps: usize,
) {
    eprintln!("Testing {b0:?} between {t0} and {t1} is {b1:?}");
    let d_sq = super::max_distance_sq_subsection(b0, b1, t0.into(), t1.into(), steps);
    assert!(
        d_sq < 0.001_f32.into(),
        "Max distance squared (actually {d_sq}) between subsections exceeeded"
    );
}

/// Test that Beziers are approximately equal
pub fn test_beziers_approx_eq<
    B0: BezierEval<F, [F; D]> + BezierSplit<F> + std::fmt::Debug,
    B1: BezierEval<F, [F; D]> + BezierSplit<F> + std::fmt::Debug,
    F: Num + From<f32>,
    const D: usize,
>(
    b0: &B0,
    b1: &B1,
) {
    eprintln!("Comparing whole Bezier {b0:?} {b1:?}");
    test_subsection(b0, b1, 0.0, 1.0, 100);
    eprintln!("Comparing subsection from 0.0 -> 1.0");
    test_subsection(
        b0,
        &b1.section(0.0_f32.into(), 1.0_f32.into()),
        0.0,
        1.0,
        100,
    );
    eprintln!("Comparing subsection from 0.1 -> 0.4");
    test_subsection(
        b0,
        &b1.section(0.1_f32.into(), 0.4_f32.into()),
        0.1,
        0.4,
        100,
    );
    eprintln!("Comparing first half from 0.0 -> 0.5");
    test_subsection(b0, &b1.split().0, 0.0, 0.5, 100);
    eprintln!("Comparing second half from 0.5 -> 1.0");
    test_subsection(b0, &b1.split().1, 0.5, 1.0, 100);

    eprintln!("Comparing first half of split_at 0.2");
    test_subsection(b0, &b1.split_at(0.2_f32.into()).0, 0.0, 0.2, 100);
    eprintln!("Comparing second half of split_at 0.2");
    test_subsection(b0, &b1.split_at(0.2_f32.into()).1, 0.2, 1.0, 100);
}
