use bezier_nd::Bezier;
use geo_nd::Float;
use geo_nd::{
    matrix,
    vector::{self, reduce},
    FArray,
};

pub fn float_iter<F: Float>(t0: F, t1: F, n: usize) -> impl Iterator<Item = F> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= F::zero(), "Float range must be positive");
    (0..n).map(move |i: usize| (<f32 as Into<F>>::into((i as f32))) / (r * ((n - 1) as f32).into()))
}

#[track_caller]
pub fn assert_near_identity<F: Float>(n: usize, m: &[F]) {
    let eps = (1E-4_f32).into();
    for i in 0..n * n {
        if i % (n + 1) == 0 {
            assert!((m[i] - F::one()).abs() < eps);
        } else {
            assert!(m[i].abs() < eps);
        }
    }
}

#[track_caller]
pub fn assert_near_equal<F: Float>(m0: &[F], m1: &[F]) {
    let eps = (1E-4_f32).into();

    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        assert!(
            (*v0 - *v1).abs() < eps,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

#[track_caller]
pub fn assert_near_equal_scale<F: Float>(m0: &[F], m1: &[F], scale: F) {
    let eps = (1E-4_f32).into();

    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        assert!(
            ((*v0) * scale - *v1).abs() < eps,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

pub fn test_subsection(b: &Bezier<f64, 2>, sub: &Bezier<f64, 2>, t0: f64, t1: f64) {
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

pub fn test_beziers_approx_eq(b0: &Bezier<f64, 2>, b1: &Bezier<f64, 2>) {
    test_subsection(b0, b1, 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.0, 1.0), 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.1, 0.4), 0.1, 0.4);
    test_subsection(b0, &b1.bisect().0, 0.0, 0.5);
    test_subsection(b0, &b1.bisect().1, 0.5, 1.0);
}
