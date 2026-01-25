//a Imports
use bezier_nd::Bezier;
use bezier_nd::BezierEval;
use bezier_nd::BezierSplit;
use geo_nd::{vector, FArray, Float, Vector};

// type Point<F:Float> = FArray<F,2>;
// type B<F:Float> = Bezier<F, Point<F>, 2>;

// type Point<const D:usize> = geometry::simd::F32x2Vec2;
// type B<const D: usize> = Bezier<f32, Point<D>, D>;

///a 'equality' tests
//fi vec_eq
pub fn vec_eq<F: Float, const D: usize>(v0: &[F; D], v1: &[F; D]) {
    #[allow(non_snake_case)]
    let EPSILON: F = (1E-5_f32).into();
    let d = vector::distance(v0, v1);
    assert!(d < EPSILON, "mismatch in {:?} {:?}", v0, v1);
}

//fi pt_eq
pub fn pt_eq<F: Float, const D: usize>(v: &[F; D], x: F, y: F) {
    #[allow(non_snake_case)]
    let EPSILON: F = (1E-5_f32).into();
    assert!(
        (v[0] - x).abs() < EPSILON,
        "mismatch in x {:?} {:?} {:?}",
        v,
        x,
        y
    );
    assert!(
        (v[1] - y).abs() < EPSILON,
        "mismatch in y {:?} {:?} {:?}",
        v,
        x,
        y
    );
}

//fi approx_eq
pub fn approx_eq<F: Float>(a: F, b: F, tolerance: F, msg: &str) {
    assert!((a - b).abs() < tolerance, "{} {:?} {:?}", msg, a, b);
}

//fi bezier_eq
pub fn bezier_eq<F: Float, const D: usize>(bez: &Bezier<F, D>, v: Vec<[F; D]>) {
    assert_eq!(bez.degree(), 4, "bezier_eq works only for cubics");
    vec_eq(bez.control_point(0), &v[0].into());
    vec_eq(bez.control_point(1), &v[1].into());
    vec_eq(bez.control_point(2), &v[2].into());
    vec_eq(bez.control_point(3), &v[3].into());
}

//a Bezier test subfns
//fi bezier_straight_as
// Check that the bezier 'is_straight' is true for >=straightness, false for < straightness
fn bezier_straight_as<F: Float>(bezier: &Bezier<F, 2>, straightness: F)
where
    FArray<F, 2>: Vector<F, 2>,
{
    for i in 0..30 {
        let s: F = (1.4_f32).powi(i - 15).into();
        println!("{} {} {}", s, straightness, bezier.is_straight(s));
        assert_eq!(
            straightness < s,
            bezier.is_straight(s),
            "Bezier {bezier} .is_straight({s}) failed for {straightness}",
        );
    }
}

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
        println!("t {} : {:?} : {:?}", t, p0, p1);
        vec_eq(&b0.point_at(t), &p0);
        vec_eq(&b1.point_at(t), &p1);
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
        println!("t {} : {:?} : {:?}", t, p, pb);
        let close_enough = EPSILON;
        approx_eq(
            p[0],
            pb[0],
            close_enough,
            &format!("Bezier split x {} {} {} : {} : {}", t, t0, t1, bezier, b),
        );
        approx_eq(
            p[1],
            pb[1],
            close_enough,
            &format!("Bezier split y {} {} {} : {} : {}", t, t0, t1, bezier, b),
        );
    }
}

//fi bezier_lines_within_straightness
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
fn bezier_lines_within_straightness<F: Float>(bezier: &Bezier<F, 2>, straightness: F)
where
    FArray<F, 2>: Vector<F, 2>,
{
    const NPTS: isize = 1000;
    const NSEG_PTS: isize = 100;
    let bezier_pts: Vec<FArray<F, 2>> = (0..NPTS)
        .map(|i: isize| {
            let t: F = ((i as f32) / (NPTS as f32)).into();
            bezier.point_at(t).into()
        })
        .collect();
    let mut segment_pts = Vec::<FArray<F, 2>>::new();
    for (p0, p1) in bezier.as_lines(straightness) {
        let p0: FArray<F, 2> = p0.into();
        let p1: FArray<F, 2> = p1.into();
        let dp = p1 - p0;
        for i in 0..NSEG_PTS {
            let t: F = ((i as f32) / 100.0).into();
            segment_pts.push(p0 + dp * t);
        }
    }
    let mut max_excursion = F::zero();
    for p in &bezier_pts {
        let min_d = segment_pts
            .iter()
            .fold(1.0E10_f32.into(), |md: F, q| md.min(p.distance(q)));
        max_excursion = max_excursion.max(min_d);
    }
    dbg!(straightness, max_excursion);
    assert!(
        max_excursion < straightness,
        "Worst case Min dist from bezier to line segments should be less than straightness {} {}",
        max_excursion,
        straightness.sqrt()
    );

    let mut max_excursion = F::zero();
    for p in &segment_pts {
        let min_d = bezier_pts
            .iter()
            .fold(1.0E10_f32.into(), |md: F, q| md.min(p.distance(q)));
        max_excursion = max_excursion.max(min_d);
    }
    dbg!(straightness, max_excursion);
    assert!(
        max_excursion < straightness,
        "Worst case Min dist from line segments to bezier should be less than straightness {} {}",
        max_excursion,
        straightness.sqrt()
    );
}

//a Test
//fi test_line
fn test_line() {
    let p0 = [0., 0.];
    let p1 = [10., 0.];
    let p2 = [10., 1.];
    let b01 = Bezier::line(&p0, &p1);
    let b02 = Bezier::line(&p0, &p2);

    pt_eq(&b01.point_at(0.), p0[0], p0[1]);
    pt_eq(
        &b01.point_at(0.5),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    pt_eq(&b01.point_at(1.), p1[0], p1[1]);
    pt_eq(&b02.point_at(0.), p0[0], p0[1]);
    pt_eq(
        &b02.point_at(0.5),
        (p0[0] + p2[0]) / 2.,
        (p0[1] + p2[1]) / 2.,
    );
    pt_eq(&b02.point_at(1.), p2[0], p2[1]);

    pt_eq(&b01.split().0.point_at(0.), p0[0], p0[1]);
    pt_eq(
        &b01.split().0.point_at(1.),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    pt_eq(
        &b01.split().1.point_at(0.),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    pt_eq(&b01.split().1.point_at(1.), p1[0], p1[1]);

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

    pt_eq(&b01.derivative_at(0.).1, p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b01.derivative_at(0.5).1, p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b01.derivative_at(1.0).1, p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b02.derivative_at(0.).1, p2[0] - p0[0], p2[1] - p0[1]);
    pt_eq(&b02.derivative_at(0.5).1, p2[0] - p0[0], p2[1] - p0[1]);
    pt_eq(&b02.derivative_at(1.0).1, p2[0] - p0[0], p2[1] - p0[1]);

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
    pt_eq(&b.point_at(0.), p0[0], p0[1]);
    pt_eq(
        &b.point_at(0.5),
        (p0[0] + p2[0]) / 4. + p1[0] / 2.,
        (p0[1] + p2[1]) / 4. + p1[1] / 2.,
    );
    pt_eq(&b.point_at(1.), p2[0], p2[1]);

    // Check that bisecting into two beziers yields the same points
    does_bisect(&b);

    // Check that splitting the bezier into a smaller section has the same points as this
    does_split(&b, 0., 1.);
    does_split(&b, 0.1, 0.3);
    does_split(&b, 0.3, 0.7);
    does_split(&b, 0.7, 1.0);

    // Check that the tangent at t=0 is the direction to the control point
    pt_eq(
        &b.derivative_at(0.).1,
        1. * (p1[0] - p0[0]),
        1. * (p1[1] - p0[1]),
    );
    // pt_eq( &b.derivative_at(0.5), p1[0]-p0[0], p1[1]-p0[1] );
    // Check that the tangent at t=1 is the direction to the control point
    pt_eq(
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
    assert_eq!(
        v.len(),
        4, // was 1 before traits
        "We know that at straightness 0.1 there must be 1 line segments"
    );

    let mut v = Vec::new();
    for (a, _b) in b.as_lines(0.001) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        128, // was 53 before traits
        "We know that at straightness 0.001  there must be 53 line segments"
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
        4, // was 3 before traits
        "We know that at straightness 0.1 there must be 3 line segments"
    );

    let mut v = Vec::new();
    for (a, b) in b.as_lines(0.001) {
        v.push((a, b));
    }
    assert_eq!(
        v.len(),
        128, // was 8 before traits
        "We know that at straightness 0.001  there must be 8 line segments"
    );
}

//fi test_cubic
fn test_cubic() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [6., 1.].into();
    let p3: FArray<f32, 2> = [20., 5.].into();
    let b = Bezier::cubic(&p0, &p1, &p2, &p3);

    pt_eq(&b.point_at(0.), p0[0], p0[1]);
    pt_eq(&b.point_at(1.), p3[0], p3[1]);

    pt_eq(&b.derivative_at(0.).1, p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b.derivative_at(1.0).1, p3[0] - p2[0], p3[1] - p2[1]);

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
    approx_eq(
        0.5,
        x.length(0.001) / PI,
        0.001,
        "Length of 90-degree arc of circle radius 1 should be PI/2",
    );

    approx_eq(
        0.5,
        x.t_of_distance(0.001, PI / 4.).0,
        0.001,
        "t of half-way round 90-degree arc of circle radius 1",
    );
    approx_eq(
        0.245,
        x.t_of_distance(0.001, PI / 8.).0,
        0.001,
        "t of quarter-way round 90-degree arc of circle radius 1",
    );
    approx_eq(
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
        22, // was 3 before changing to traits
        "We know that at straightness 0.5 there should be 3 line segments"
    );

    v.clear();
    for (a, _b) in b.as_lines(0.01) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        980, // was 23 before changing to traits
        "We know that at straightness 0.01 there should be 23 line segments"
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
    vec_eq(b.control_point(0), sb.control_point(0));
    vec_eq(b.control_point(1), sb.control_point(1));

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
    vec_eq(b.control_point(0), sb.control_point(0));
    vec_eq(b.control_point(1), sb.control_point(1));
    vec_eq(b.control_point(2), sb.control_point(2));
    vec_eq(b.control_point(3), sb.control_point(3));

    bezier_straight_as(&Bezier::cubic(&p0, &p1, &p2, &p4), 0.6);
    bezier_straight_as(&Bezier::cubic(&sp0, &sp1, &sp2, &sp4), 6.0);
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
        approx_eq(radius, r, 0.000015, "Radius of arc should be as requested");
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
    pt_eq(&c, 0., 0.);
    pt_eq(&x.point_at(0.), 1., 0.);
    pt_eq(&x.point_at(1.), 0., 1.);

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
    pt_eq(&c, 0., 0.);
    pt_eq(&x.point_at(0.), 0., 1.);
    pt_eq(&x.point_at(1.), 1., 0.);

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
    pt_eq(&c, 0., 0.);
    pt_eq(&x.point_at(0.), r_sqrt2, -r_sqrt2);
    pt_eq(&x.point_at(1.), r_sqrt2, r_sqrt2);

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

//fi test_within_straightness
fn test_within_straightness() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [10., 1.].into();
    let p3: FArray<f32, 2> = [20., 0.].into();
    // let p4: FArray<f32, 2> = [20., 1.].into();

    let b = Bezier::quadratic(&p0, &p3, &p1);
    bezier_lines_within_straightness(&b, 0.1);

    let b = Bezier::cubic(&p0, &p1, &p2, &p3);
    bezier_lines_within_straightness(&b, 0.1);

    let b = Bezier::cubic(&p0, &p2, &p1, &p3);
    bezier_lines_within_straightness(&b, 0.1);

    let b = Bezier::cubic(&p3, &p2, &p1, &p0);
    bezier_lines_within_straightness(&b, 0.1);

    let b = Bezier::cubic(&p2, &p3, &p0, &p1);
    bezier_lines_within_straightness(&b, 0.1);

    let b = Bezier::cubic(&p0, &p3, &p3, &p1);
    bezier_lines_within_straightness(&b, 0.1);
}

//a Tests
#[test]
fn test_f32_within_straightness() {
    test_within_straightness();
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
fn test_f32_straight_as() {
    test_straight_as();
}

#[test]
fn test_f32_arc() {
    test_arc();
}

#[test]
fn test_f32_round() {
    test_round();
}
