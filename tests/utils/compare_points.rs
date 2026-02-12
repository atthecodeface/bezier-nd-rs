use bezier_nd::{bernstein_fns, BezierEval, BezierSplit};

use bezier_nd::Float;
use bezier_nd::Num;
use geo_nd::{matrix, vector};

#[track_caller]
pub fn vec_eq<F: Num, const D: usize>(v0: &[F; D], v1: &[F; D]) {
    #[allow(non_snake_case)]
    let EPSILON: F = (1E-8_f32).into();
    let d = vector::distance_sq(v0, v1);
    assert!(d < EPSILON, "mismatch in {:?} {:?}", v0, v1);
}

//fi pt_eq
pub fn pt_eq<F: Num, const D: usize>(v: &[F; D], x: F, y: F) {
    #[allow(non_snake_case)]
    let EPSILON: F = (1E-5_f32).into();
    assert!(
        super::abs(v[0] - x) < EPSILON,
        "mismatch in x {:?} {:?} {:?}",
        v,
        x,
        y
    );
    assert!(
        super::abs(v[1] - y) < EPSILON,
        "mismatch in y {:?} {:?} {:?}",
        v,
        x,
        y
    );
}

#[track_caller]
pub fn approx_eq<F: Num, I: Into<F>>(a: F, b: F, tolerance: I, msg: &str) {
    let tolerance = tolerance.into();
    let diff = a - b;
    let diff = if diff < F::ZERO { -diff } else { diff };
    assert!(
        diff < tolerance,
        "{msg}: a {a:?} != b {b:?}: difference {diff} tolerance {tolerance}"
    );
}

/// Assert that a matrix is near the identity matrix
#[track_caller]
pub fn assert_near_identity<N: geo_nd::Num>(n: usize, m: &[N]) {
    let eps = N::from_usize(1).unwrap() / N::from_usize(10000).unwrap();
    for i in 0..n * n {
        let exp = {
            if i % (n + 1) == 0 {
                N::one()
            } else {
                N::zero()
            }
        };
        assert!(m[i] - exp > -eps && m[i] - exp < eps);
    }
}

/// Assert that a matrix is nearly equal to another
#[track_caller]
pub fn assert_near_equal<N: Num>(m0: &[N], m1: &[N]) {
    let eps = N::from_usize(1).unwrap() / N::from_usize(10000).unwrap();

    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        let dv = *v0 - *v1;
        assert!(
            dv >= -eps && dv <= eps,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

/// Assert that a scaled matrix is nearly equal to another unscaled matrix
#[track_caller]
pub fn assert_near_equal_scale<N: Num, I: Into<N> + Copy>(m0: &[N], m1: &[N], scale: I) {
    let eps = N::from_usize(1).unwrap() / N::from_usize(10000).unwrap();
    let scale = scale.into();
    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        let dv = (*v0) * scale - *v1;
        assert!(
            dv >= -eps && dv <= eps,
            "Data {i} is mismatch in {m0:?} * {scale} <> {m1:?}"
        );
    }
}

/// Assert that a scaled matrix is nearly equal to another unscaled matrix
#[track_caller]
pub fn assert_near_equal_sorted_scale<N: Num, const D: usize, I: Into<N> + Copy>(
    m0: &[N; D],
    m1: &[N; D],
    scale: I,
) {
    let mut m0 = *m0;
    let mut m1 = *m1;
    m0.sort_by(|a, b| a.partial_cmp(b).unwrap());
    m1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_near_equal_scale(&m0, &m1, scale);
}
