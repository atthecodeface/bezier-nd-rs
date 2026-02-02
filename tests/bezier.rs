//a Imports
use bezier_nd::Bezier;
use bezier_nd::BezierEval;
use bezier_nd::BezierSplit;
use bezier_nd::Float;
mod utils;
use bezier_nd::BezierMinMax;
use geo_nd::{vector, FArray, Vector};

//fi does_bisect
/// Check that bisecting into two beziers yields the same points
fn does_bisect<F: Float, const D: usize>(bezier: &Bezier<F, D>)
where
    FArray<F, D>: Vector<F, D>,
{
    let (b0, b1) = bezier.split();
    println!("Test bisection of {} into {}, {}", bezier, b0, b1);
    for i in 0..21 {
        let t: F = ((i as f32) / 20.0).into();
        let half = 0.5_f32.into();
        let p0 = bezier.point_at(t * half);
        let p1 = bezier.point_at(t * half + half);
        //println!("t {} : {:?} : {:?}", t, p0, p1);
        utils::vec_eq(&b0.point_at(t), &p0);
        utils::vec_eq(&b1.point_at(t), &p1);
    }
}

//fi does_split
/// Test that between the sub-bezier between t of t0 and t1 for bezier has
///  the same points as that bezier itself
fn does_split<F: Float, const D: usize>(bezier: &Bezier<F, D>, t0: F, t1: F)
where
    FArray<F, D>: Vector<F, D>,
{
    #[allow(non_snake_case)]
    let EPSILON: F = (1E-5_f32).into();
    let b = bezier.bezier_between(t0, t1);
    for i in 0..21 {
        let bt: F = ((i as f32) / 20.0).into();
        let t = t0 + (t1 - t0) * bt;
        let p = bezier.point_at(t);
        let pb = b.point_at(bt);
        // println!("t {} : {:?} : {:?}", t, p, pb);
        let close_enough = EPSILON;
        utils::approx_eq(
            p[0],
            pb[0],
            close_enough,
            &format!("Bezier split x {} {} {} : {} : {}", t, t0, t1, bezier, b),
        );
        utils::approx_eq(
            p[1],
            pb[1],
            close_enough,
            &format!("Bezier split y {} {} {} : {} : {}", t, t0, t1, bezier, b),
        );
    }
}

fn min_max_coords<
    F: bezier_nd::Num,
    const D: usize,
    B: BezierMinMax<F> + BezierEval<F, [F; D]> + std::fmt::Debug,
>(
    bezier: &B,
) {
    eprintln!("Min/max of {bezier:?}");
    let pts = utils::BezierPtSet::of_point_at(bezier, 1000);
    let bbox = pts.bbox();
    eprintln!("Bbox of pts set = {bbox:?}");
    for d in 0..D {
        let (t_min_d, f_min_d) = bezier.t_coord_at_min_max(false, d).unwrap();
        let (t_max_d, f_max_d) = bezier.t_coord_at_min_max(true, d).unwrap();
        eprintln!("BBox coord {d} min {f_min_d} <> max {f_max_d} at {t_min_d} {t_max_d}");
        utils::approx_eq(f_min_d, bezier.point_at(t_min_d)[d], 1E-6, "Pt at t_min_d");
        utils::approx_eq(f_max_d, bezier.point_at(t_max_d)[d], 1E-6, "Pt at t_min_d");
        utils::approx_eq(f_min_d, bbox.0[d], 1E-4, "Bbox min coord d");
        utils::approx_eq(f_max_d, bbox.1[d], 1E-4, "Bbox min coord d");
    }
}

//a Test
//fi test_line
fn test_line() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [10., 1.];
    let b01 = Bezier::line(&p0, &p1);
    let b02 = Bezier::line(&p0, &p2);

    utils::pt_eq(&b01.point_at(0.), p0[0], p0[1]);
    utils::pt_eq(
        &b01.point_at(0.5),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    utils::pt_eq(&b01.point_at(1.), p1[0], p1[1]);
    utils::pt_eq(&b02.point_at(0.), p0[0], p0[1]);
    utils::pt_eq(
        &b02.point_at(0.5),
        (p0[0] + p2[0]) / 2.,
        (p0[1] + p2[1]) / 2.,
    );
    utils::pt_eq(&b02.point_at(1.), p2[0], p2[1]);

    utils::pt_eq(&b01.split().0.point_at(0.), p0[0], p0[1]);
    utils::pt_eq(
        &b01.split().0.point_at(1.),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    utils::pt_eq(
        &b01.split().1.point_at(0.),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    utils::pt_eq(&b01.split().1.point_at(1.), p1[0], p1[1]);

    does_split(&b01, 0., 1.);
    does_split(&b01, 0.1, 0.3);
    does_split(&b01, 0.3, 0.7);
    does_split(&b01, 0.7, 1.0);

    does_split(&b02, 0., 1.);
    does_split(&b02, 0.1, 0.3);
    does_split(&b02, 0.3, 0.7);
    does_split(&b02, 0.7, 1.0);

    does_bisect(&b01);
    does_bisect(&b02);

    utils::pt_eq(&b01.derivative_at(0.).1, p1[0] - p0[0], p1[1] - p0[1]);
    utils::pt_eq(&b01.derivative_at(0.5).1, p1[0] - p0[0], p1[1] - p0[1]);
    utils::pt_eq(&b01.derivative_at(1.0).1, p1[0] - p0[0], p1[1] - p0[1]);
    utils::pt_eq(&b02.derivative_at(0.).1, p2[0] - p0[0], p2[1] - p0[1]);
    utils::pt_eq(&b02.derivative_at(0.5).1, p2[0] - p0[0], p2[1] - p0[1]);
    utils::pt_eq(&b02.derivative_at(1.0).1, p2[0] - p0[0], p2[1] - p0[1]);

    let mut v = Vec::new();
    v.clear();
    for (a, _b) in b01.as_lines(0.1) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        1,
        "We know that at any straightness there must be 1 line segments"
    );

    min_max_coords(&[p0, p1]);
    min_max_coords(&[p1, p0]);
    min_max_coords(&[p0, p2]);
    min_max_coords(&[p2, p0]);
}

//fi test_quadratic
fn test_quadratic() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [10., 1.];
    let p3 = [20., 0.];

    // b is a quad with endpoints (0,0).and (10,1), with the middle control point of (10.0)
    let b = Bezier::quadratic(&p0, &p1, &p2);

    // These are generically true...
    utils::pt_eq(&b.point_at(0.), p0[0], p0[1]);
    utils::pt_eq(
        &b.point_at(0.5),
        (p0[0] + p2[0]) / 4. + p1[0] / 2.,
        (p0[1] + p2[1]) / 4. + p1[1] / 2.,
    );
    utils::pt_eq(&b.point_at(1.), p2[0], p2[1]);

    // Check that bisecting into two beziers yields the same points
    does_bisect(&b);

    // Check that splitting the bezier into a smaller section has the same points as this
    does_split(&b, 0., 1.);
    does_split(&b, 0.1, 0.3);
    does_split(&b, 0.3, 0.7);
    does_split(&b, 0.7, 1.0);

    // Check that the tangent at t=0 is the direction to the control point
    utils::pt_eq(
        &b.derivative_at(0.).1,
        1. * (p1[0] - p0[0]),
        1. * (p1[1] - p0[1]),
    );
    // pt_eq( &b.derivative_at(0.5), p1[0]-p0[0], p1[1]-p0[1] );
    // Check that the tangent at t=1 is the direction to the control point
    utils::pt_eq(
        &b.derivative_at(1.0).1,
        1. * (p2[0] - p1[0]),
        1. * (p2[1] - p1[1]),
    );

    // Convert to lines
    //
    // Specific to this curve which is 0,0 -> 10,1.0 with control 10.0,0.0 so is always within 1.0
    let mut v = Vec::new();
    for (p0, p1) in b.as_lines(1.0) {
        // eprintln!("added {p0:?} {p1:?} to bezier {b}");
        v.push((p0, p1));
    }
    eprintln!("{v:?}");
    assert_eq!(
        v.len(),
        1,
        "We know that at straightness 0.1 there must be 1 line segments"
    );

    let mut v = Vec::new();
    for (a, _b) in b.as_lines(0.001) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        10,
        "We know that at straightness_sq 0.001 there must be 10 line segments"
    );

    // Use a Bezier with a control point outside the endpoints
    //
    // Convert to lines
    //
    // Specific to this curve which is 0,0 -> 10,0.0 with control 20.0,0.0
    //
    // This is effectively a line from (0,0) to (13.332,0) then back to (10,0)
    let b = Bezier::quadratic(&p0, &p3, &p1);
    let mut v = Vec::new();
    for (a, b) in b.as_lines(1.0) {
        v.push((a, b));
    }
    assert_eq!(
        v.len(),
        3,
        "We know that at straightness 0.1 there must be 3 line segments"
    );

    let mut v = Vec::new();
    for (a, b) in b.as_lines(0.001) {
        v.push((a, b));
    }
    assert_eq!(
        v.len(),
        6,
        "We know that at straightness 0.001 there must be 6 line segments"
    );

    min_max_coords(&[p0, p1, p2]);
    min_max_coords(&[p1, p0, p2]);
    min_max_coords(&[p0, p2, p1]);
    min_max_coords(&[p2, p0, p1]);
}

//fi test_cubic
fn test_cubic() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [6., 1.].into();
    let p3: FArray<f32, 2> = [20., 5.].into();
    let b = Bezier::cubic(&p0, &p1, &p2, &p3);

    utils::pt_eq(&b.point_at(0.), p0[0], p0[1]);
    utils::pt_eq(&b.point_at(1.), p3[0], p3[1]);

    utils::pt_eq(&b.derivative_at(0.).1, p1[0] - p0[0], p1[1] - p0[1]);
    utils::pt_eq(&b.derivative_at(1.0).1, p3[0] - p2[0], p3[1] - p2[1]);

    does_bisect(&b);

    does_split(&b, 0., 1.);
    does_split(&b, 0.1, 0.3);
    does_split(&b, 0.3, 0.7);
    does_split(&b, 0.7, 1.0);

    let x = Bezier::arc(
        (90.0f32).to_radians(),
        1.,
        &p0,
        &[1., 0.].into(),
        &[0., 1.].into(),
        0.,
    );
    println!("{}", x);
    use std::f32::consts::PI;
    utils::approx_eq(
        0.5,
        x.length(0.001) / PI,
        0.001,
        "Length of 90-degree arc of circle radius 1 should be PI/2",
    );

    utils::approx_eq(
        0.5,
        x.t_of_distance(0.001, PI / 4.).0,
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
    );

    let mut v = Vec::new();
    v.clear();
    for (a, _b) in b.as_lines(0.5) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        6, // was 3 before changing to traits
        "We know that at straightness_sq of 0.5 there should be 6 line segments"
    );

    v.clear();
    for (a, _b) in b.as_lines(0.01) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        13,
        "We know that at straightness_sq of 0.01 there should be 13 line segments"
    );

    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [6., 1.];
    let p3 = [20., 5.];

    min_max_coords(&[p0, p1, p2, p3]);
    min_max_coords(&[p1, p0, p2, p3]);
    min_max_coords(&[p0, p3, p2, p1]);
    min_max_coords(&[p2, p3, p0, p1]);
}

//fi arc_ave_square_error
fn arc_ave_square_error(
    arc: &Bezier<f32, 2>,
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
        let x = Bezier::arc(
            angle.to_radians(),
            radius,
            &center,
            &x_axis,
            &y_axis,
            rotate.to_radians(),
        );
        let (c, r) = x.center_radius_of_bezier_arc();
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

    let x = Bezier::of_round_corner(&(x_axis + y_axis), &(y_axis * 3.), &(x_axis * 0.5), 1.);
    println!("of_round_corner(&[1.,1.], &[0.,3.], &[0.5,0.], 1.) : {}", x);
    let (c, r) = x.center_radius_of_bezier_arc();
    let e2 = arc_ave_square_error(&x, &c, r, 0., 1., 10);
    println!("arc_ave_square_error {}", e2);
    assert!(
        e2 < 0.01,
        "Eccentricity of 90deg rounded corner should be <1%"
    );
    utils::pt_eq(&c, 0., 0.);
    utils::pt_eq(&x.point_at(0.), 1., 0.);
    utils::pt_eq(&x.point_at(1.), 0., 1.);

    let x = Bezier::of_round_corner(&(x_axis + y_axis), &(x_axis * 0.5), &(y_axis * 3.), 1.);
    println!(
        "of_round_corner(&[1.,1.], &[0.5,0.], &[0.,3.],  1.) : {}",
        x
    );
    let (c, r) = x.center_radius_of_bezier_arc();
    let e2 = arc_ave_square_error(&x, &c, r, 0., 1., 10);
    println!("arc_ave_square_error {}", e2);
    assert!(
        e2 < 0.01,
        "Eccentricity of 90deg rounded corner should be <1%"
    );
    utils::pt_eq(&c, 0., 0.);
    utils::pt_eq(&x.point_at(0.), 0., 1.);
    utils::pt_eq(&x.point_at(1.), 1., 0.);

    let x = Bezier::of_round_corner(
        &(x_axis * sqrt2),
        &(x_axis + y_axis),
        &(x_axis - y_axis),
        1.,
    );
    println!(
        "of_round_corner(&[sqrt2,0.], &[1.,1.], &[1.,-1.], 1.) : {}",
        x
    );
    let (c, r) = x.center_radius_of_bezier_arc();
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
        let x = Bezier::of_round_corner(
            &x_axis,
            &(x_axis * dx + y_axis * dy),
            &(x_axis * dx - y_axis * dy),
            radius,
        );
        let (c, r) = x.center_radius_of_bezier_arc();
        let e2 = arc_ave_square_error(&x, &c, radius, 0., 1., 10);
        println!(
            "arc_ave_square_error for dx {} dy {} is {} (radius found {})",
            dx, dy, e2, r
        );
        assert!(
            e2 < ase,
            "Eccentricity of rounded corner {} {} should be < {}",
            dx,
            dy,
            ase
        );
    }
}

#[test]
fn test_f32_line() {
    test_line();
}

#[test]
fn test_f32_quadratic() {
    test_quadratic();
}

#[test]
fn test_f32_cubic() {
    test_cubic();
}

#[test]
fn test_f32_arc() {
    test_arc();
}

#[test]
fn test_f32_round() {
    test_round();
}
