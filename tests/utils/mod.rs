#![allow(dead_code)]
#![allow(unused_imports)]
use bezier_nd::BezierEval;
use bezier_nd::BezierSplit;

use bezier_nd::Float;
use bezier_nd::{Bezier, Num};
use geo_nd::{matrix, vector};

mod point_set;
mod random;
pub use point_set::BezierPtSet;
pub use random::*;
mod measures;
pub use measures::*;

pub fn vec_eq<F: Num, const D: usize>(v0: &[F; D], v1: &[F; D]) {
    #[allow(non_snake_case)]
    let EPSILON: F = (1E-8_f32).into();
    let d = vector::distance_sq(v0, v1);
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
pub fn approx_eq<F: Num, I: Into<F>>(a: F, b: F, tolerance: I, msg: &str) {
    let tolerance = tolerance.into();
    let diff = a - b;
    let diff = if diff < F::ZERO { -diff } else { diff };
    assert!(diff < tolerance, "{} {:?} {:?}", msg, a, b);
}

//fi bezier_eq
pub fn bezier_eq<F: Float, const D: usize>(bez: &Bezier<F, D>, v: Vec<[F; D]>) {
    assert_eq!(bez.degree(), 4, "bezier_eq works only for cubics");
    vec_eq(bez.control_point(0), &v[0].into());
    vec_eq(bez.control_point(1), &v[1].into());
    vec_eq(bez.control_point(2), &v[2].into());
    vec_eq(bez.control_point(3), &v[3].into());
}

/// Iterate in 'n' steps from t0 to t1 inclusive
pub fn float_iter<N: Num>(t0: N, t1: N, n: usize) -> impl Iterator<Item = N> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= N::zero(), "Float range must be positive");
    (0..n).map(move |i: usize| r * N::from_usize(i).unwrap() / N::from_usize(n - 1).unwrap())
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
pub fn assert_near_equal_scale<N: Num>(m0: &[N], m1: &[N], scale: N) {
    let eps = N::from_usize(1).unwrap() / N::from_usize(10000).unwrap();

    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        let dv = (*v0) * scale - *v1;
        assert!(
            dv >= -eps && dv <= eps,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}

/// Test that a subsection of a Bezier has the same points as another between t0 and t1
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

/// Test that Beziers are approximately equal
pub fn test_beziers_approx_eq(b0: &Bezier<f64, 2>, b1: &Bezier<f64, 2>) {
    test_subsection(b0, b1, 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.0, 1.0), 0.0, 1.0);
    test_subsection(b0, &b1.bezier_between(0.1, 0.4), 0.1, 0.4);
    test_subsection(b0, &b1.split().0, 0.0, 0.5);
    test_subsection(b0, &b1.split().1, 0.5, 1.0);
}

/// Generate Bernstein matrices for a given degree and values of t
#[track_caller]
pub fn generate_bernstein_matrix<F: Float>(matrix: &mut [F], degree: usize, ts: &[F]) {
    assert_eq!(
        matrix.len(),
        (degree + 1) * ts.len(),
        "Trying to generate an {} by {} matrix, so must be given a matrix of that size",
        ts.len(),
        degree + 1
    );
    for (r, t) in ts.iter().enumerate() {
        for c in 0..(degree + 1) {
            matrix[r * (degree + 1) + c] =
                bezier_nd::bernstein::bezier_fns::basis_coeff(degree, c, *t);
        }
    }
}

/// Generate Bernstein matrices for a given degree and values of t
pub fn bernstein_basis_coeff_br<N: geo_nd::Num + From<i64>>(degree: usize, i: usize, t: N) -> N {
    let u = N::ONE - t;
    let mut result = N::ONE;
    for c in 0..degree {
        if c < i {
            result *= t;
        } else {
            result *= u;
        }
    }
    // Multiply by n! / (n-i)! / i!
    //
    // This is multiply by (n-j) for j = 0..i-1 inclusive and divide by j for j = 1..i inclusive
    for j in 0..=i {
        if j < i {
            let f: N = ((degree - j) as i64).into();
            result *= f;
        }

        if j >= 1 {
            let f: N = (j as i64).into();
            result /= f;
        }
    }
    result
}

/// Generate a Bernstein matrix for a given degree and values of t
///
/// This generates an N * (degree+1) matrix which can be applied to 'degree' control points
/// to generate the positions on the Bezier (of that degree with those control points) at
/// N specified values of 't' (given by ts)
#[track_caller]
pub fn generate_bernstein_matrix_br<N: geo_nd::Num + From<i64>>(
    matrix: &mut [N],
    degree: usize,
    ts: &[N],
) {
    assert_eq!(
        matrix.len(),
        (degree + 1) * ts.len(),
        "Trying to generate an {} by {} matrix, so must be given a matrix of that size",
        ts.len(),
        degree + 1
    );
    for (r, t) in ts.iter().enumerate() {
        for c in 0..(degree + 1) {
            matrix[r * (degree + 1) + c] = bernstein_basis_coeff_br(degree, c, *t);
        }
    }
}

/// Generate an elevate-by-N matrix for a given degree
///
/// This is a (degree+1+N)*(degree+1) matrix that should be applied to the control points
/// to generate a new array of control points of the elevated Bezier
#[track_caller]
#[must_use]
pub fn generate_elevate_by_n_matrix<N: bezier_nd::Num>(degree: usize, by_n: usize) -> Vec<N> {
    use bezier_nd::bernstein::bezier_fns::generate_elevate_by_one_matrix;

    assert!(by_n >= 1);
    let mut result = vec![N::ZERO; (degree + 1) * (degree + 1 + by_n)];
    if by_n == 1 {
        let scale = generate_elevate_by_one_matrix(&mut result, degree);
        for e in result.iter_mut() {
            *e /= scale;
        }
    } else {
        // m_a is degree+2 * degree+1
        // m_a is degree+1+by_n * degree+2
        let m_a = generate_elevate_by_n_matrix(degree, 1);
        let m_b = generate_elevate_by_n_matrix(degree + 1, by_n - 1);
        matrix::multiply_dyn(
            degree + 1 + by_n,
            degree + 2,
            degree + 1,
            &m_b,
            &m_a,
            &mut result,
        );
    }
    result
}
