//a Imports
use bezier_nd::Num;
use bezier_nd::{BasicBezier, BezierIterationType, BezierND};
mod utils;
use rand::prelude::*;

/// Test a Bezier that is a line between two points
fn test_line_between<F: Num + From<f32>, B: BasicBezier<F, [F; 2]> + std::fmt::Debug>(
    bezier: &B,
    pts: [[F; 2]; 2],
) {
    let p0 = &pts[0];
    let p1 = &pts[1];
    let pm = [
        (p0[0] + p1[0]) / 2.0_f32.into(),
        (p0[1] + p1[1]) / 2.0_f32.into(),
    ];
    let dp = [(p1[0] - p0[0]), (p1[1] - p0[1])];

    assert_eq!(bezier.degree(), 1);
    assert_eq!(bezier.num_control_points(), 2);
    utils::bezier_eq(bezier, &pts);

    assert_eq!(bezier.closeness_sq_to_line(), F::ZERO);
    assert_eq!(bezier.dc_sq_from_line(), F::ZERO);

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
    for (_, a, _, _b) in bezier.as_t_lines(BezierIterationType::ClosenessSq(0.1_f32.into())) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        1,
        "We know that at any straightness there must be 1 line segments"
    );

    utils::assert_min_max_coords(bezier);
}

fn test_line<
    F: Num + From<f32>,
    R: Rng,
    D: Distribution<f32>,
    B: BasicBezier<F, [F; 2]> + std::fmt::Debug,
    M: Fn([[F; 2]; 2]) -> B,
>(
    rng: &mut R,
    distribution: &D,

    create_bezier: M,
) {
    for _ in 0..10 {
        let pts: [[F; 2]; 2] = utils::new_random_point_array(rng, distribution);
        let bezier = create_bezier(pts);
        test_line_between(&bezier, pts);
    }
}

#[test]
fn test_lines() {
    let mut rng = utils::make_random("test_lines_seed");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    test_line::<f32, _, _, _, _>(&mut rng, &distribution, |pts| {
        let b = BezierND::<_, 5, _>::new(&pts);
        b
    });
    test_line::<f32, _, _, _, _>(&mut rng, &distribution, |pts| {
        let b: Vec<_> = pts.into();
        b
    });
    test_line::<f32, _, _, _, _>(&mut rng, &distribution, |pts| pts);
}
