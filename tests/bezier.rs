//a Imports
use bezier_nd::Bezier;
use geo_nd::{FArray, Float, Vector};

// type Point<F:Float> = FArray<F,2>;
// type B<F:Float> = Bezier<F, Point<F>, 2>;

// type Point<const D:usize> = geometry::simd::F32x2Vec2;
// type B<const D: usize> = Bezier<f32, Point<D>, D>;

///a 'equality' tests
//fi vec_eq
pub fn vec_eq<F: Float, V: Vector<F, D>, const D: usize>(v0: &V, v1: &V) {
    #[allow(non_snake_case)]
    let EPSILON: F = (1E-5_f32).into();
    let d = v0.distance(v1);
    assert!(d < EPSILON, "mismatch in {:?} {:?}", v0, v1);
}

//fi pt_eq
pub fn pt_eq<F: Float, V: Vector<F, D>, const D: usize>(v: &V, x: F, y: F) {
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
pub fn bezier_eq<F: Float, V: Vector<F, D>, const D: usize>(bez: &Bezier<F, V, D>, v: Vec<[F; D]>) {
    assert_eq!(bez.degree(), 4, "bezier_eq works only for cubics");
    vec_eq(bez.borrow_pt(0), &v[0].into());
    vec_eq(bez.borrow_pt(2), &v[1].into());
    vec_eq(bez.borrow_pt(3), &v[2].into());
    vec_eq(bez.borrow_pt(1), &v[3].into());
}

//a Bezier test subfns
//fi bezier_straight_as
fn bezier_straight_as<F: Float, V: Vector<F, 2>>(bezier: &Bezier<F, V, 2>, straightness: F) {
    for i in 0..30 {
        let s: F = (1.4_f32).powi(i - 15).into();
        println!("{} {} {}", s, straightness, bezier.is_straight(s));
        assert_eq!(
            straightness < s,
            bezier.is_straight(s),
            "Bezier {} .is_straight({}) failed for {}",
            bezier,
            s,
            straightness
        );
    }
}

//fi does_bisect
fn does_bisect<F: Float, V: Vector<F, D>, const D: usize>(bezier: &Bezier<F, V, D>) {
    let (b0, b1) = bezier.bisect();
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
fn does_split<F: Float, V: Vector<F, D>, const D: usize>(bezier: &Bezier<F, V, D>, t0: F, t1: F) {
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

//a Test
//fi test_line
fn test_line<Point>()
where
    Point: Vector<f32, 2>,
{
    let p0: Point = [0., 0.].into();
    let p1: Point = [10., 0.].into();
    let p2: Point = [10., 1.].into();
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

    pt_eq(&b01.bisect().0.point_at(0.), p0[0], p0[1]);
    pt_eq(
        &b01.bisect().0.point_at(1.),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    pt_eq(
        &b01.bisect().1.point_at(0.),
        (p0[0] + p1[0]) / 2.,
        (p0[1] + p1[1]) / 2.,
    );
    pt_eq(&b01.bisect().1.point_at(1.), p1[0], p1[1]);

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

    pt_eq(&b01.tangent_at(0.), p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b01.tangent_at(0.5), p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b01.tangent_at(1.0), p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b02.tangent_at(0.), p2[0] - p0[0], p2[1] - p0[1]);
    pt_eq(&b02.tangent_at(0.5), p2[0] - p0[0], p2[1] - p0[1]);
    pt_eq(&b02.tangent_at(1.0), p2[0] - p0[0], p2[1] - p0[1]);

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
fn test_quadratic<Point>()
where
    Point: Vector<f32, 2>,
{
    let p0: Point = [0., 0.].into();
    let p1: Point = [10., 0.].into();
    let p2: Point = [10., 1.].into();
    let b = Bezier::quadratic(&p0, &p1, &p2);

    pt_eq(&b.point_at(0.), p0[0], p0[1]);
    pt_eq(
        &b.point_at(0.5),
        (p0[0] + p2[0]) / 4. + p1[0] / 2.,
        (p0[1] + p2[1]) / 4. + p1[1] / 2.,
    );
    pt_eq(&b.point_at(1.), p2[0], p2[1]);

    does_bisect(&b);

    does_split(&b, 0., 1.);
    does_split(&b, 0.1, 0.3);
    does_split(&b, 0.3, 0.7);
    does_split(&b, 0.7, 1.0);

    pt_eq(
        &b.tangent_at(0.),
        1. * (p1[0] - p0[0]),
        1. * (p1[1] - p0[1]),
    );
    // pt_eq( &b.tangent_at(0.5), p1[0]-p0[0], p1[1]-p0[1] );
    pt_eq(
        &b.tangent_at(1.0),
        1. * (p2[0] - p1[0]),
        1. * (p2[1] - p1[1]),
    );

    let mut v = Vec::new();
    v.clear();
    for (a, _b) in b.as_lines(0.1) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        1,
        "We know that at straightness 0.1 there must be 1 line segments"
    );

    let mut v = Vec::new();
    v.clear();
    for (a, _b) in b.as_lines(0.01) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        52,
        "We know that at straightness 0.01  there must be 52 line segments"
    );
}

//fi test_cubic
fn test_cubic<Point>()
where
    Point: Vector<f32, 2>,
{
    let p0: Point = [0., 0.].into();
    let p1: Point = [10., 0.].into();
    let p2: Point = [6., 1.].into();
    let p3: Point = [20., 5.].into();
    let b = Bezier::cubic(&p0, &p1, &p2, &p3);

    pt_eq(&b.point_at(0.), p0[0], p0[1]);
    pt_eq(&b.point_at(1.), p3[0], p3[1]);

    pt_eq(&b.tangent_at(0.), p1[0] - p0[0], p1[1] - p0[1]);
    pt_eq(&b.tangent_at(1.0), p3[0] - p2[0], p3[1] - p2[1]);

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
    for (a, _b) in b.as_lines(0.1) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        3,
        "We know that at straightness 0.1 there should be 3 line segments"
    );

    v.clear();
    for (a, _b) in b.as_lines(0.01) {
        v.push(a);
    }
    assert_eq!(
        v.len(),
        24,
        "We know that at straightness 0.01 there should be 24 line segments"
    );
}

//fi test_straight
fn test_straight<Point>()
where
    Point: Vector<f32, 2>,
{
    let p0: Point = [0., 0.].into();
    let p1: Point = [10., 0.].into();
    let p2: Point = [10., 1.].into();
    let p3: Point = [20., 0.].into();
    let p4: Point = [20., 1.].into();

    let sp0 = p0 * 10.;
    let sp1 = p1 * 10.;
    let sp2 = p2 * 10.;
    let sp3 = p3 * 10.;
    let sp4 = p4 * 10.;

    let mut b = Bezier::line(&p0, &p1);
    let sb = Bezier::line(&sp0, &sp1);
    b.scale(10.);
    vec_eq(b.borrow_pt(0), sb.borrow_pt(0));
    vec_eq(b.borrow_pt(1), sb.borrow_pt(1));

    bezier_straight_as(&Bezier::line(&p0, &p1), 1E-10);
    bezier_straight_as(&Bezier::line(&p0, &p2), 1E-10);
    bezier_straight_as(&Bezier::line(&p0, &p3), 1E-10);
    bezier_straight_as(&Bezier::line(&p0, &p4), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp1), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp2), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp3), 1E-10);
    bezier_straight_as(&Bezier::line(&sp0, &sp4), 1E-10);

    bezier_straight_as(&Bezier::quadratic(&p0, &p1, &p3), 1E-10);
    bezier_straight_as(&Bezier::quadratic(&sp0, &sp1, &sp3), 1E-10);

    bezier_straight_as(&Bezier::quadratic(&p0, &p2, &p3), 0.05);
    bezier_straight_as(&Bezier::quadratic(&sp0, &sp2, &sp3), 0.05);

    bezier_straight_as(&Bezier::quadratic(&p0, &p1, &p4), 0.03);
    bezier_straight_as(&Bezier::quadratic(&sp0, &sp1, &sp4), 0.03);

    bezier_straight_as(&Bezier::cubic(&p0, &p1, &p2, &p3), 0.05);
    bezier_straight_as(&Bezier::cubic(&sp0, &sp1, &sp2, &sp3), 0.05);

    let mut b = Bezier::cubic(&p0, &p1, &p2, &p4);
    let sb = Bezier::cubic(&sp0, &sp1, &sp2, &sp4);
    b.scale(10.);
    vec_eq(b.borrow_pt(0), sb.borrow_pt(0));
    vec_eq(b.borrow_pt(1), sb.borrow_pt(1));
    vec_eq(b.borrow_pt(2), sb.borrow_pt(2));
    vec_eq(b.borrow_pt(3), sb.borrow_pt(3));

    bezier_straight_as(&Bezier::cubic(&p0, &p1, &p2, &p4), 0.065);
    bezier_straight_as(&Bezier::cubic(&sp0, &sp1, &sp2, &sp4), 0.065);
}

//fi arc_ave_square_error
fn arc_ave_square_error<Point>(
    arc: &Bezier<f32, Point, 2>,
    center: &Point,
    radius: f32,
    t0: f32,
    t1: f32,
    n: usize,
) -> f32
where
    Point: Vector<f32, 2>,
{
    let delta_t = (t1 - t0) / (n - 1) as f32;
    let mut error2 = 0.;
    for i in 0..n {
        let e = arc.point_at((i as f32) * delta_t + t0).distance(center) - radius;
        error2 += e * e;
    }
    error2 / (n as f32)
}

//fi test_arc
fn test_arc<Point>()
where
    Point: Vector<f32, 2>,
{
    let x_axis = [1.0, 0.0].into();
    let y_axis = [0.0, 1.0].into();

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
        let center: Point = [cx, cy].into();
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
        println!("c {} r {} arc_ave_square_error {}", c, r, e2);
        approx_eq(radius, r, 0.000015, "Radius of arc should be as requested");
        assert!(
            c.distance(&center) < 0.0001,
            "Center {} should match {}",
            c,
            center
        );
        assert!(e2 < 0.001, "Eccentricity of arc should be < 0.001");
    }
}

//fi test_round
fn test_round<Point>()
where
    Point: Vector<f32, 2>,
{
    let x_axis: Point = [1.0, 0.0].into();
    let y_axis: Point = [0.0, 1.0].into();

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

//a Tests
#[test]
fn test_f32_line() {
    test_line::<FArray<f32, 2>>();
}

#[test]
fn test_f32_quadratic() {
    test_quadratic::<FArray<f32, 2>>();
}

#[test]
fn test_f32_cubic() {
    test_cubic::<FArray<f32, 2>>();
}

#[test]
fn test_f32_straight() {
    test_straight::<FArray<f32, 2>>();
}

#[test]
fn test_f32_arc() {
    test_arc::<FArray<f32, 2>>();
}

#[test]
fn test_f32_round() {
    test_round::<FArray<f32, 2>>();
}
