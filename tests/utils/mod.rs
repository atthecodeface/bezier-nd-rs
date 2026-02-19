#![allow(dead_code)]
#![allow(unused_imports)]
use bezier_nd::{bernstein_fns, BezierEval, BezierSplit};

use bezier_nd::Num;
use geo_nd::{matrix, vector};

mod compare_beziers;
mod compare_points;
mod measures;
mod minmax_bbox;
mod point_set;
mod random;

pub use compare_beziers::*;
pub use compare_points::*;
pub use measures::*;
pub use minmax_bbox::*;
pub use point_set::BezierPtSet;
pub use random::*;
pub mod display;

pub fn abs<F: Num>(f: F) -> F {
    if f < F::ZERO {
        -f
    } else {
        f
    }
}

/// Iterate in 'n' steps from t0 to t1 inclusive
pub fn float_iter_between<N: Num>(t0: N, t1: N, n: usize) -> impl Iterator<Item = N> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= N::zero(), "Float range must be positive");
    (0..n).map(move |i: usize| t0 + r * N::from_usize(i).unwrap() / N::from_usize(n - 1).unwrap())
}

/// Iterate in 'n' steps from 0 to 1 inclusive
pub fn float_iter<F: Num>(n: usize) -> impl Iterator<Item = F> {
    float_iter_between(F::ZERO, F::ONE, n)
}

/// Generate Bernstein matrices for a given degree and values of t
#[track_caller]
pub fn generate_bernstein_matrix<F: Num>(matrix: &mut [F], degree: usize, ts: &[F]) {
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
pub fn bernstein_basis_coeff_br<N: Num>(degree: usize, i: usize, t: N) -> N {
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
            let f = N::of_usize(degree - j);
            result *= f;
        }

        if j >= 1 {
            let f = N::of_usize(j);
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
pub fn generate_bernstein_matrix_br<N: Num>(matrix: &mut [N], degree: usize, ts: &[N]) {
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
pub fn generate_elevate_by_n_matrix<N: Num>(degree: usize, by_n: usize) -> Vec<N> {
    use bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix;
    assert!(by_n >= 1);
    if by_n == 1 {
        let (scale, mut result) = generate_elevate_by_one_matrix(degree);
        for e in result.iter_mut() {
            *e /= scale;
        }
        result
    } else {
        // m_a is degree+2 * degree+1
        // m_a is degree+1+by_n * degree+2
        let mut result = vec![N::ZERO; (degree + 1) * (degree + 1 + by_n)];
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
        result
    }
}
