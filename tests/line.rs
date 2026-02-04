//a Imports
use bezier_nd::BezierEval;
use bezier_nd::BezierSplit;
use bezier_nd::Num;
use bezier_nd::{BasicBezier, Bezier};
mod utils;
use bezier_nd::BezierMinMax;

/// Test a Bezier that is a line between two points
fn test_line_between<F: Num, B: BasicBezier<F, [F; 2]> + std::fmt::Debug>(
    bezier: &B,
    p0: &[F; 2],
    p1: &[F; 2],
) {
    let pm = [
        (p0[0] + p1[0]) / 2.0_f32.into(),
        (p0[1] + p1[1]) / 2.0_f32.into(),
    ];
    let dp = [(p1[0] - p0[0]), (p1[1] - p0[1])];

    assert_eq!(bezier.degree(), 1);
    assert_eq!(bezier.num_control_points(), 2);

    utils::vec_eq(&bezier.endpoints().0, p0);
    utils::vec_eq(&bezier.endpoints().1, p1);

    utils::vec_eq(&bezier.point_at(F::ZERO), p0);
    utils::vec_eq(&bezier.point_at(F::ONE), p1);
    utils::vec_eq(&bezier.point_at(0.5_f32.into()), &pm);

    utils::vec_eq(&bezier.split().0.point_at(F::ZERO), p0);
    utils::vec_eq(&bezier.split().0.point_at(F::ONE), &pm);
    utils::vec_eq(&bezier.split().1.point_at(F::ZERO), &pm);
    utils::vec_eq(&bezier.split().1.point_at(F::ONE), p1);

    for (t0, t1) in [(0., 1.), (0.1, 0.3), (0.3, 0.7), (0.8, 1.0)] {
        utils::test_subsection(bezier, &bezier.section(t0.into(), t1.into()), t0, t1, 100);
    }

    let (b0, b1) = bezier.split();
    utils::test_subsection(bezier, &b0, 0., 0.5, 100);
    utils::test_subsection(bezier, &b1, 0.5, 1.0, 100);

    let scale = bezier.derivative_at(F::ZERO).0;
    eprintln!("Scale {scale}");
    utils::assert_near_equal_scale(&bezier.derivative_at(F::ZERO).1, &dp, scale);
    utils::assert_near_equal_scale(&bezier.derivative_at(0.5_f32.into()).1, &dp, scale);
    utils::assert_near_equal_scale(&bezier.derivative_at(F::ONE).1, &dp, scale);

    let mut v = Vec::new();
    v.clear();
    for (a, _b) in bezier.as_lines(0.1_f32.into()) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        1,
        "We know that at any straightness there must be 1 line segments"
    );

    utils::assert_min_max_coords(bezier);
}

fn test_line() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [10., 1.];
    let b01: Bezier<_, _> = [p0, p1].as_ref().try_into().unwrap();
    let b02: Bezier<_, _> = [p0, p2].into();

    test_line_between(&b01, &p0, &p1);
    test_line_between(&[p0, p1], &p0, &p1);
    //    test_line_between(&vec![p0, p1], &p0, &p1);

    test_line_between(&b02, &p0, &p2);
    test_line_between(&[p0, p2], &p0, &p2);
    //    test_line_between(&vec![p0, p2], &p0, &p2);
}

#[test]
fn test_f32_line() {
    test_line();
}
