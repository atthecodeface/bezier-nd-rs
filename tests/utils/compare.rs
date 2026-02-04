use bezier_nd::{bernstein_fns, BezierEval, BezierMinMax, BezierSection, BezierSplit};

use bezier_nd::Float;
use bezier_nd::{Bezier, Num};
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
        "approx_eq mismatch: {} diff {diff} tolerance {tolerance} {:?}<>{:?}",
        msg,
        a,
        b
    );
}

#[track_caller]
pub fn bezier_eq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bez: &B, v: &[[F; D]]) {
    for (i, v) in v.iter().enumerate() {
        vec_eq(bez.control_point(i), v);
    }
}
