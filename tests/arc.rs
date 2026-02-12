//a Imports
use bezier_nd::BezierEval;
use bezier_nd::{bernstein_fns, metrics, BezierFlatIterator};
mod utils;
use geo_nd::{vector, FArray};

fn test_cubic_arc() {
    let p0 = [0., 0.];

    let x = bernstein_fns::arc::arc(
        (90.0f32).to_radians(),
        1.,
        &p0,
        &[1., 0.].into(),
        &[0., 1.].into(),
        0.,
    );
    eprintln!("{:?}", x);
    use std::f32::consts::PI;
    utils::approx_eq(
        0.5,
        metrics::length_of_lines(x.as_lines(1E-6)) / PI,
        0.001,
        "Length of 90-degree arc of circle radius 1 should be PI/2",
    );

    /*
    utils::approx_eq(
        0.5,
        x.t__of_distance(0.001, PI / 4.).0,
        0.001,
        "t of half-way round 90-degree arc of circle radius 1",
    );
    utils::approx_eq(
        0.245,
        x.t_of_distance(0.001, PI / 8.).0,
        0.001,
        "t of quarter-way round 90-degree arc of circle radius 1",
    );
    utils::approx_eq(
        0.755,
        x.t_of_distance(0.001, PI * 3. / 8.).0,
        0.001,
        "t of three-quarters-way round 90-degree arc of circle radius 1",
    );*/
}

//fi arc_ave_square_error
fn arc_ave_square_error<B: BezierEval<f32, [f32; 2]>>(
    arc: &B,
    center: &[f32; 2],
    radius: f32,
    t0: f32,
    t1: f32,
    n: usize,
) -> f32 {
    let delta_t = (t1 - t0) / (n - 1) as f32;
    let mut error2 = 0.;
    for i in 0..n {
        let p = arc.point_at((i as f32) * delta_t + t0);
        let e = vector::distance(&p, center) - radius;
        error2 += e * e;
    }
    error2 / (n as f32)
}

//fi test_arc
fn test_arc() {
    let x_axis = [1.0, 0.0];
    let y_axis = [0.0, 1.0];

    for (angle, radius, cx, cy, rotate) in [
        (90.0f32, 1., 0., 0., 0.0f32),
        (90.0f32, 1., 0., 0., 10.),
        (90.0f32, 1., 0., 0., 20.),
        (90.0f32, 1., 0., 0., 30.),
        (10.0f32, 1., 0., 0., 0.),
        (20.0f32, 1., 0., 0., 0.),
        (30.0f32, 1., 0., 0., 0.),
        (40.0f32, 1., 0., 0., 0.),
        (50.0f32, 1., 0., 0., 0.),
        (60.0f32, 1., 0., 0., 0.),
        (60.0f32, 1., 0., 0., 0.),
        (70.0f32, 1., 0., 0., 0.),
        (80.0f32, 1., 0., 0., 0.),
        (100.0f32, 1., 0., 0., 0.),
        (110.0f32, 1., 0., 0., 0.),
        (120.0f32, 1., 0., 0., 0.),
        (130.0f32, 1., 0., 0., 0.),
        (140.0f32, 1., 0., 0., 0.),
        (150.0f32, 1., 0., 0., 0.),
        (160.0f32, 1., 0., 0., 0.),
        (170.0f32, 1., 0., 0., 0.),
        (10.0f32, 2., 0., 0., 0.),
        (30.0f32, 2., 0., 0., 0.),
        (60.0f32, 2., 0., 0., 0.),
        (80.0f32, 2., 0., 0., 0.),
        (90.0f32, 2., 0., 0., 30.),
        (100.0f32, 2., 0., 0., 0.),
        (120.0f32, 2., 0., 0., 0.),
        (150.0f32, 2., 0., 0., 0.),
        (170.0f32, 2., 0., 0., 0.),
        (10.0f32, 2., 1., 9., 0.),
        (30.0f32, 2., 2., 8., 0.),
        (60.0f32, 2., 3., 7., 0.),
        (80.0f32, 2., 4., 6., 0.),
        (90.0f32, 2., 5., 5., 30.),
        (100.0f32, 2., 6., 4., 0.),
        (120.0f32, 2., 7., 3., 0.),
        (150.0f32, 2., 8., 2., 0.),
        (170.0f32, 2., 9., 1., 0.),
    ] {
        let center: FArray<f32, 2> = [cx, cy].into();
        let x = bernstein_fns::arc::arc(
            angle.to_radians(),
            radius,
            &center,
            &x_axis,
            &y_axis,
            rotate.to_radians(),
        );

        let (c, r) = bernstein_fns::arc::center_radius_of_bezier_arc(&x);
        let e2 = arc_ave_square_error(&x, &c, r, 0., 1., 10);
        println!("c {c:?} r {r} arc_ave_square_error {e2}");
        utils::approx_eq(radius, r, 0.000015, "Radius of arc should be as requested");
        assert!(
            vector::distance(&c, &center) < 0.0001,
            "Center {c:?} should match {center:?}",
        );
        assert!(e2 < 0.001, "Eccentricity of arc should be < 0.001");
    }
}

//fi test_round
fn test_round() {
    let x_axis: FArray<f32, 2> = [1.0, 0.0].into();
    let y_axis: FArray<f32, 2> = [0.0, 1.0].into();

    let sqrt2 = 2.0_f32.sqrt();
    let r_sqrt2 = 1.0 / sqrt2;

    let x = bernstein_fns::arc::of_round_corner(
        &(x_axis + y_axis),
        &(y_axis * 3.),
        &(x_axis * 0.5),
        1.,
    );
    println!("of_round_corner(&[1.,1.], &[0.,3.], &[0.5,0.], 1.) : {x:?}",);
    let (c, r) = bernstein_fns::arc::center_radius_of_bezier_arc(&x);
    let e2 = arc_ave_square_error(&x, &c, r, 0., 1., 10);
    println!("arc_ave_square_error {}", e2);
    assert!(
        e2 < 0.01,
        "Eccentricity of 90deg rounded corner should be <1%"
    );
    utils::pt_eq(&c, 0., 0.);
    utils::pt_eq(&x.point_at(0.), 1., 0.);
    utils::pt_eq(&x.point_at(1.), 0., 1.);

    let x = bernstein_fns::arc::of_round_corner(
        &(x_axis + y_axis),
        &(x_axis * 0.5),
        &(y_axis * 3.),
        1.,
    );
    println!("of_round_corner(&[1.,1.], &[0.5,0.], &[0.,3.],  1.) : {x:?}",);
    let (c, r) = bernstein_fns::arc::center_radius_of_bezier_arc(&x);
    let e2 = arc_ave_square_error(&x, &c, r, 0., 1., 10);
    println!("arc_ave_square_error {}", e2);
    assert!(
        e2 < 0.01,
        "Eccentricity of 90deg rounded corner should be <1%"
    );
    utils::pt_eq(&c, 0., 0.);
    utils::pt_eq(&x.point_at(0.), 0., 1.);
    utils::pt_eq(&x.point_at(1.), 1., 0.);

    let x = bernstein_fns::arc::of_round_corner(
        &(x_axis * sqrt2),
        &(x_axis + y_axis),
        &(x_axis - y_axis),
        1.,
    );
    println!("of_round_corner(&[sqrt2,0.], &[1.,1.], &[1.,-1.], 1.) : {x:?}",);
    let (c, r) = bernstein_fns::arc::center_radius_of_bezier_arc(&x);
    let e2 = arc_ave_square_error(&x, &c, r, 0., 1., 10);
    println!("arc_ave_square_error {}", e2);
    assert!(
        e2 < 0.01,
        "Eccentricity of 90deg rounded corner should be <1%"
    );
    utils::pt_eq(&c, 0., 0.);
    utils::pt_eq(&x.point_at(0.), r_sqrt2, -r_sqrt2);
    utils::pt_eq(&x.point_at(1.), r_sqrt2, r_sqrt2);

    let radius = 0.5; // F::frac(1,2);
                      // let radius = 1.0;// F::frac(1,2);
    for (dx, dy, ase) in [
        (1., 1., 0.001),
        (1., 30., 0.001),
        (1., 15., 0.001),
        (1., 11., 0.001),
        (1., 9., 0.001),
        (1., 7., 0.001),
        (1., 5., 0.001),
        (1., 3., 0.001),
        (1., -3., 0.001),
        (1., 2.5, 0.001),
        (1., 2., 0.001),
        (1., 1.5, 0.001),
        (1., 1., 0.001),
        (1.5, 1., 0.001),
        (2.0, 1., 0.001),
        (2.5, 1., 0.001),
        (3., 1., 0.001),
        (5., 1., 0.001),
        (7., 1., 0.001),
        (11., 1., 0.001),
        (15., 1., 0.001),
        (30., 1., 0.001),
    ] {
        let x = bernstein_fns::arc::of_round_corner(
            &x_axis,
            &(x_axis * dx + y_axis * dy),
            &(x_axis * dx - y_axis * dy),
            radius,
        );
        let (c, r) = bernstein_fns::arc::center_radius_of_bezier_arc(&x);
        let e2 = arc_ave_square_error(&x, &c, radius, 0., 1., 10);
        eprintln!("arc_ave_square_error for dx {dx} dy {dy} is {e2} (radius found {r})",);
        assert!(
            e2 < ase,
            "Eccentricity of rounded corner {dx} {dy} should be < {ase}",
        );
    }
}

#[test]
fn test_f32_cubic_arc() {
    test_cubic_arc();
}

#[test]
fn test_f32_arc() {
    test_arc();
}

#[test]
fn test_f32_round() {
    test_round();
}
