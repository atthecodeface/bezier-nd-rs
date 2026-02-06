//a Imports
use bezier_nd::{Approximation, BasicBezier, Bezier, BezierEval};
use bezier_nd::{Float, Num};
mod utils;
use geo_nd::vector;
use rand::prelude::*;

fn test_cubic_with<F: Num, B: BasicBezier<F, [F; 2]> + std::fmt::Debug>(
    bezier: &B,
    pts: [[F; 2]; 4],
) {
    let p0 = &pts[0];
    let p1 = &pts[3];
    let pm = vector::sum_scaled(
        &pts,
        &[
            0.125_f32.into(),
            0.375_f32.into(),
            0.375_f32.into(),
            0.125_f32.into(),
        ],
    );

    assert_eq!(bezier.degree(), 3);
    assert_eq!(bezier.num_control_points(), 4);
    for i in 0..bezier.num_control_points() {
        assert_eq!(bezier.control_point(i), &pts[i]);
    }

    utils::vec_eq(&bezier.endpoints().0, p0);
    utils::vec_eq(&bezier.endpoints().1, p1);

    assert!(bezier.dc_sq_from_line() >= bezier.closeness_sq_to_line());

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

    // Check that the tangent at t=0 is the direction to the control point
    // Check that the tangent at t=1 is the direction to the control point
    let scale = bezier.derivative_at(F::ZERO).0;
    eprintln!("Scale {scale}");
    utils::assert_near_equal_scale(
        &bezier.derivative_at(F::ZERO).1,
        &vector::sum_scaled(&pts, &[-scale, scale, F::ZERO, F::ZERO]),
        scale,
    );
    utils::assert_near_equal_scale(
        &bezier.derivative_at(F::ONE).1,
        &vector::sum_scaled(&pts, &[F::ZERO, F::ZERO, -scale, scale]),
        scale,
    );
    for t in utils::float_iter(10) {
        utils::assert_near_equal_scale(
            &bezier.derivative_at(t).1,
            &vector::scale(
                [
                    vector::sub(pts[1], &pts[0], F::ONE),
                    vector::sub(pts[2], &pts[1], F::ONE),
                    vector::sub(pts[3], &pts[2], F::ONE),
                ]
                .point_at(t),
                (3.0_f32).into(),
            ),
            bezier.derivative_at(t).0,
        );
    }
    utils::assert_min_max_coords(bezier)
}

fn test_cubic<
    F: Float,
    R: Rng,
    D: Distribution<f32>,
    B: BasicBezier<F, [F; 2]> + std::fmt::Debug,
    M: Fn([[F; 2]; 4]) -> B,
>(
    rng: &mut R,
    distribution: &D,

    create_bezier: M,
) {
    for _ in 0..10 {
        let pts: [[F; 2]; 4] = utils::new_random_point_array(rng, distribution);
        let bezier = create_bezier(pts);
        test_cubic_with(&bezier, pts);
    }
}
fn test_arc<F: Float>() {
    let bezier = bezier_nd::bernstein_fns::arc::arc(
        F::PI() / (2.0_f32.into()),
        F::ONE,
        &[F::ZERO, F::ZERO],
        &[F::ONE, F::ZERO],
        &[F::ZERO, F::ONE],
        F::ZERO,
    );

    let approx = Approximation::new(&bezier, 0.0001_f32.into());

    eprintln!("{:?}", approx);
    utils::approx_eq(
        F::PI() * 0.5_f32.into(),
        approx.section_length(F::ZERO, F::ONE),
        0.001_f32,
        "Length of 90-degree arc of circle radius 1 should be PI/2",
    );

    utils::approx_eq(
        0.5_f32.into(),
        approx.t_of_distance(F::PI() / 4.0_f32.into()).unwrap(),
        0.001,
        "t of half-way round 90-degree arc of circle radius 1",
    );
    utils::approx_eq(
        0.245_f32.into(),
        approx.t_of_distance(F::PI() / 8.0_f32.into()).unwrap(),
        0.001,
        "t of quarter-way round 90-degree arc of circle radius 1",
    );
    utils::approx_eq(
        0.755_f32.into(),
        approx
            .t_of_distance(F::PI() * (3.0_f32 / 8.0).into())
            .unwrap(),
        0.001,
        "t of three-quarters-way round 90-degree arc of circle radius 1",
    );
}

#[test]
fn test_cubics() {
    let mut rng = utils::make_random("test_cubics_seed");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    eprintln!("*************************************************\nTesting 'Bezier<>' type (from 'try_into())");
    test_cubic::<f32, _, _, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.as_ref().try_into().unwrap();
        b
    });

    eprintln!("*************************************************\nTesting 'Bezier<>' type again (from 'into')");
    test_cubic::<f32, _, _, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Vec' type");
    test_cubic::<f32, _, _, _, _>(&mut rng, &distribution, |pts| {
        let b: Vec<_> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Array' type");
    test_cubic::<f32, _, _, _, _>(&mut rng, &distribution, |pts| pts);
}

#[test]
fn test_cubic_arc() {
    eprintln!("*************************************************\nTesting 'Bezier<>' type (from 'try_into())");
    test_arc::<f32>();
}
