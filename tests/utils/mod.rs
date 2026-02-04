#![allow(dead_code)]
#![allow(unused_imports)]
use bezier_nd::{bernstein_fns, BezierEval, BezierMinMax, BezierSection, BezierSplit};

use bezier_nd::Float;
use bezier_nd::{Bezier, Num};
use geo_nd::{matrix, vector};

mod point_set;
mod random;
pub use point_set::BezierPtSet;
pub use random::*;
mod measures;
pub use measures::*;
mod compare;
pub use compare::*;

pub fn abs<F: Num>(f: F) -> F {
    if f < F::ZERO {
        -f
    } else {
        f
    }
}
/// Iterate in 'n' steps from t0 to t1 inclusive
pub fn float_iter<N: Num>(t0: N, t1: N, n: usize) -> impl Iterator<Item = N> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= N::zero(), "Float range must be positive");
    (0..n).map(move |i: usize| t0 + r * N::from_usize(i).unwrap() / N::from_usize(n - 1).unwrap())
}

#[track_caller]
/// Find the maximum distance between `b0` for t0<=t<=t1 and `b1` for 0<=t<=1 in 'steps' intervals
pub fn test_subsection<
    F: Num,
    const D: usize,
    B0: BezierEval<F, [F; D]> + std::fmt::Debug,
    B1: BezierEval<F, [F; D]> + std::fmt::Debug,
>(
    b0: &B0,
    b1: &B1,
    t0: f32,
    t1: f32,
    steps: usize,
) {
    eprintln!("Testing {b0:?} between {t0} and {t1} is {b1:?}");
    let d_sq = max_distance_sq_subsection(b0, b1, t0.into(), t1.into(), steps);
    assert!(
        d_sq < 0.001_f32.into(),
        "Max distance squared (actually {d_sq}) between subsections exceeeded"
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

/// Test that Beziers are approximately equal
pub fn test_beziers_approx_eq<
    B0: BezierEval<F, [F; D]> + BezierSplit + BezierSection<F> + std::fmt::Debug,
    B1: BezierEval<F, [F; D]> + BezierSplit + BezierSection<F> + std::fmt::Debug,
    F: Num,
    const D: usize,
>(
    b0: &B0,
    b1: &B1,
) {
    eprintln!("Comparing whole Bezier {b0:?} {b1:?}");
    test_subsection(b0, b1, 0.0, 1.0, 100);
    eprintln!("Comparing subsection from 0.0 -> 1.0");
    test_subsection(
        b0,
        &b1.section(0.0_f32.into(), 1.0_f32.into()),
        0.0,
        1.0,
        100,
    );
    eprintln!("Comparing subsection from 0.1 -> 0.4");
    test_subsection(
        b0,
        &b1.section(0.1_f32.into(), 0.4_f32.into()),
        0.1,
        0.4,
        100,
    );
    eprintln!("Comparing first half from 0.0 -> 0.5");
    test_subsection(b0, &b1.split().0, 0.0, 0.5, 100);
    eprintln!("Comparing second half from 0.5 -> 1.0");
    test_subsection(b0, &b1.split().1, 0.5, 1.0, 100);
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
            matrix[r * (degree + 1) + c] = bernstein_fns::basis_coeff(degree, c, *t);
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
    use bernstein_fns::generate_elevate_by_one_matrix;

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

pub fn assert_min_max_coords<
    F: bezier_nd::Num,
    const D: usize,
    B: BezierMinMax<F> + BezierEval<F, [F; D]> + std::fmt::Debug,
>(
    bezier: &B,
) {
    eprintln!("Min/max of {bezier:?}");
    let pts = BezierPtSet::of_point_at(bezier, 1000);
    let bbox = pts.bbox();
    eprintln!("Bbox of pts set = {bbox:?}");
    for d in 0..D {
        let (t_min_d, f_min_d) = bezier.t_coord_at_min_max(false, d).unwrap();
        let (t_max_d, f_max_d) = bezier.t_coord_at_min_max(true, d).unwrap();
        eprintln!("BBox coord {d} min {f_min_d} <> max {f_max_d} at {t_min_d} {t_max_d}");
        approx_eq(f_min_d, bezier.point_at(t_min_d)[d], 1E-6, "Pt at t_min_d");
        approx_eq(f_max_d, bezier.point_at(t_max_d)[d], 1E-6, "Pt at t_min_d");
        approx_eq(f_min_d, bbox.0[d], 1E-4, "Bbox min coord d");
        approx_eq(f_max_d, bbox.1[d], 1E-4, "Bbox min coord d");
    }
}
