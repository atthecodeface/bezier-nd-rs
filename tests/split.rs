//a Imports
use bezier_nd::Bezier;
use geo_nd::{vector, FArray};

fn float_iter(t0: f32, t1: f32, n: usize) -> impl Iterator<Item = f32> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= 0.0, "Float range must be positive");
    (0..n).map(move |i: usize| (i as f32) / (r * ((n - 1) as f32)))
}

fn test_subsection(b: &Bezier<f32, 2>, sub: &Bezier<f32, 2>, t0: f32, t1: f32) {
    eprintln!("Testing subsections of beziers {b} {sub} {t0} {t1}");
    for sub_t in float_iter(0.0, 1.0, 100) {
        let t = t0 + (t1 - t0) * sub_t;
        let p = b.point_at(t);
        let sub_p = sub.point_at(sub_t);
        let d2 = vector::distance_sq(&p, &sub_p);
        assert!(d2<1E-4,
        "Points at bezier {t} : {p:?} and subbezier {sub_t} : {sub_p:?} should be roughly the same but have distance {d2}");
    }
}

fn test_beziers_approx_eq(b0: &Bezier<f32, 2>, b1: &Bezier<f32, 2>) {
    test_subsection(b0, b1, 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.0, 1.0), 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.1, 0.4), 0.1, 0.4);
    test_subsection(b0, &b1.bisect().0, 0.0, 0.5);
    test_subsection(b0, &b1.bisect().1, 0.5, 1.0);
}

#[test]
fn bisect() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [6., 1.].into();
    let p3: FArray<f32, 2> = [20., 5.].into();

    test_beziers_approx_eq(
        &Bezier::cubic(&p0, &p1, &p2, &p3),
        &Bezier::cubic(&p0, &p1, &p2, &p3),
    );
    test_beziers_approx_eq(
        &Bezier::cubic(&p2, &p1, &p3, &p0),
        &Bezier::cubic(&p2, &p1, &p3, &p0),
    );
    test_beziers_approx_eq(
        &Bezier::quadratic(&p2, &p1, &p3),
        &Bezier::quadratic(&p2, &p1, &p3),
    );
    test_beziers_approx_eq(
        &Bezier::quadratic(&p0, &p1, &p2),
        &Bezier::quadratic(&p0, &p1, &p2),
    );
    test_beziers_approx_eq(&Bezier::line(&p0, &p2), &Bezier::line(&p0, &p2));
}

#[test]
fn elevate() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let p2: FArray<f32, 2> = [6., 1.].into();
    let p3: FArray<f32, 2> = [20., 5.].into();

    let b = Bezier::quadratic(&p0, &p1, &p2);
    let b2 = b.clone().elevate();
    test_beziers_approx_eq(&b, &b2);

    let b = Bezier::quadratic(&p3, &p1, &p2);
    let b2 = b.clone().elevate();
    test_beziers_approx_eq(&b, &b2);

    let b = Bezier::line(&p3, &p1);
    let b2 = b.clone().elevate();
    test_beziers_approx_eq(&b, &b2);
}
