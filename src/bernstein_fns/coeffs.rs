//! A collection of functions for use with Bernstein polynomial Beziers

use crate::constants::BINOMIALS_U;
use crate::{Float, Num};
use geo_nd::vector;

/// Calculate the ith Bernstein polynomial coefficient at 't' for a given degree.
///
/// This is (1-t)^(degree-i) * t^i * (degree! / (i!.(degree-i)!) )
///
/// This requires F:Float as it uses 'powi', which is not available if only F:Num
#[inline]
#[track_caller]
pub fn basis_coeff_enum<F: Float>(degree: usize, t: F) -> impl Iterator<Item = (usize, F)> {
    assert!(
        degree < BINOMIALS_U.len(),
        "Maximum degree of Bezier supported is {}",
        BINOMIALS_U.len()
    );
    let u = F::ONE - t;
    let coeffs = BINOMIALS_U[degree];
    (0..=degree).map(move |i| {
        (
            i,
            t.powi(i as i32) * u.powi((degree - i) as i32) * F::from_usize(coeffs[1 + i]).unwrap(),
        )
    })
}

/// Calculate the coefficient for the first derivative of the ith Bernstein polynomial at 't' for a given degree.
///
/// This is (i^2)/n.B[i-1,n-1](t) - (n-i)^2/n.B[i,n-1](t)
///
/// This requires F:Float as it uses 'powi', which is not available if only F:Num
#[inline]
#[track_caller]
pub fn basis_dt_coeff_enum<F: Float>(degree: usize, t: F) -> (F, impl Iterator<Item = (usize, F)>) {
    let n = degree + 1;
    let _n_f = F::from_usize(n).unwrap();
    (
        F::from_usize(n - 1).unwrap(),
        std::iter::once((0, F::ZERO))
            .chain(basis_coeff_enum(degree - 1, t))
            .zip(basis_coeff_enum(degree - 1, t).chain(std::iter::once((0, F::ZERO))))
            .enumerate()
            .map(move |(i, ((_, c0), (_i1, c1)))| {
                (
                    i,
                    (F::from_usize(i * i).unwrap() * c0
                        - F::from_usize((n - i - 1) * (n - i - 1)).unwrap() * c1),
                )
            }),
    )
}

/// Calculate the coefficient for the first derivative of the ith Bernstein polynomial at 't' for a given degree.
///
/// This is (i^2)/n.B[i-1,n-1](t) - (n-i)^2/n.B[i,n-1](t)
///
/// This requires F:Num
#[inline]
#[track_caller]
pub fn basis_dt_coeff_enum_num<F: Num>(
    degree: usize,
    t: F,
) -> (F, impl Iterator<Item = (usize, F)>) {
    let n = degree + 1;
    let _n_f = F::from_usize(n).unwrap();
    (
        F::from_usize(n - 1).unwrap(),
        std::iter::once((0, F::ZERO))
            .chain(basis_coeff_enum_num(degree - 1, t))
            .zip(basis_coeff_enum_num(degree - 1, t).chain(std::iter::once((0, F::ZERO))))
            .enumerate()
            .map(move |(i, ((_, c0), (_i1, c1)))| {
                (
                    i,
                    (F::from_usize(i * i).unwrap() * c0
                        - F::from_usize((n - i - 1) * (n - i - 1)).unwrap() * c1),
                )
            }),
    )
}

/// Calculate the ith Bernstein polynomial coefficient at 't' for a given degree.
///
/// This is (1-t)^(degree-i) * t^i * (degree! / (i!.(degree-i)!) )
///
/// This requires F:Float as it uses 'powi', which is not available if only F:Num
#[track_caller]
#[inline]
pub fn basis_coeff_enum_num<F: Num>(degree: usize, t: F) -> impl Iterator<Item = (usize, F)> {
    assert!(
        degree < BINOMIALS_U.len(),
        "Maximum degree of Bezier supported is {}",
        BINOMIALS_U.len()
    );
    let u = F::ONE - t;
    let coeffs = BINOMIALS_U[degree];
    (0..=degree).map(move |i| {
        let mut x = F::from_usize(coeffs[1 + i]).unwrap();
        for _ in 0..i {
            x *= t;
        }
        for _ in i..degree {
            x *= u;
        }
        (i, x)
    })
}

/// Calculate the ith Bernstein polynomial coefficient at 't' for a given degree.
///
/// This is (1-t)^(degree-i) * t^i * (degree! / (i!.(degree-i)!) )
///
/// This requires F:Float as it uses 'powi', which is not available if only F:Num
#[inline]
#[track_caller]
pub fn basis_coeff<F: Float>(degree: usize, i: usize, t: F) -> F {
    assert!(
        degree < BINOMIALS_U.len(),
        "Maximum degree of Bezier supported is {}",
        BINOMIALS_U.len()
    );

    let u = F::one() - t;
    let coeffs = BINOMIALS_U[degree];
    t.powi(i as i32) * u.powi((degree - i) as i32) * F::from_usize(coeffs[1 + i]).unwrap()
}

/// Calculate the control points of the nth derivative of a Bernstein Bezier
///
/// Updates the provided d_pts slice, and returns the scaling which should be applied (multiply by the result)
#[track_caller]
pub fn nth_derivative<F: Num, const D: usize>(pts: &[[F; D]], n: usize, d_pts: &mut [[F; D]]) -> F {
    assert!(
        pts.len() > n,
        "Bezier of degree {}-1 must actually be of degree {n} or higher to have an {n}th derivative",
        pts.len()
    );
    assert!(
        d_pts.len() >= pts.len()-n, // note pts.len() > n already
        "Provided storage for derivative points must be at least {} to provide {n}th derivative of degree {} Bezier",
        pts.len()-n, pts.len()-1
    );
    let degree = pts.len() - 1;
    let mut scale = F::one();
    for i in 0..n {
        scale *= F::from_usize(degree - i).unwrap();
    }
    for (i, sp) in d_pts.iter_mut().take(degree + 1 - n).enumerate() {
        let mut m1_n_positive = (n & 1) == 0;
        for (c, p) in BINOMIALS_U[n][1..].iter().zip(pts[i..].iter()) {
            let mut c = F::from_usize(*c).unwrap();
            if !m1_n_positive {
                c = -c;
            }
            *sp = vector::add(*sp, p, c);
            m1_n_positive = !m1_n_positive;
        }
    }
    scale
}

/// Elevate a Bezier by one degree, that should be reduced by 'F'
///
/// This generates and applies the elevate-by-one matrix
pub fn elevate_by_one<F: Num, const D: usize>(pts: &mut [[F; D]], ele: &mut [[F; D]]) -> F {
    assert!(
        ele.len() > pts.len(),
        "At least {} points required to elevate, but a slice with only {} was provided",
        pts.len() + 1,
        ele.len()
    );
    // n = number of points, i.e. self.pts[n-1] is the last valid point
    let n = pts.len();
    ele[0] = pts[0];
    ele[n + 1] = pts[n];
    for (j, e) in ele.iter_mut().skip(1).take(n - 1).enumerate() {
        *e = vector::add(
            vector::scale(pts[j], (j as f32).into()),
            &pts[j + 1],
            ((n + 1 - j) as f32).into(),
        );
    }
    (1.0 / ((n + 1) as f32)).into()
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
        for (c, bc) in basis_coeff_enum_num(degree, *t) {
            matrix[r * (degree + 1) + c] = bc;
        }
    }
}
