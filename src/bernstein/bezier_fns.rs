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
pub fn basis_coeff_iter<F: Float>(degree: usize, t: F) -> impl Iterator<Item = F> {
    let u = F::ONE - t;
    let coeffs = BINOMIALS_U[degree];
    (0..=degree).map(move |i| {
        t.powi(i as i32) * u.powi((degree - i) as i32) * F::from_usize(coeffs[1 + i]).unwrap()
    })
}

/// Calculate the ith Bernstein polynomial coefficient at 't' for a given degree.
///
/// This is (1-t)^(degree-i) * t^i * (degree! / (i!.(degree-i)!) )
///
/// This requires F:Float as it uses 'powi', which is not available if only F:Num
#[inline]
pub fn basis_coeff_iter_num<F: Num>(degree: usize, t: F) -> impl Iterator<Item = F> {
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
        x
    })
}

/// Calculate the ith Bernstein polynomial coefficient at 't' for a given degree.
///
/// This is (1-t)^(degree-i) * t^i * (degree! / (i!.(degree-i)!) )
///
/// This requires F:Float as it uses 'powi', which is not available if only F:Num
#[inline]
pub fn basis_coeff<F: Float>(degree: usize, i: usize, t: F) -> F {
    let u = F::one() - t;
    let coeffs = BINOMIALS_U[degree];
    t.powi(i as i32) * u.powi((degree - i) as i32) * F::from_usize(coeffs[1 + i]).unwrap()
}

/// Calculate the Mij element of the elevation matrix to elevate by one degree
/// given the current degree. The elevation matrix is (degree+2)x(degree+1)
///
/// This is given by `P'[j] = Sum((n i).(1 j-i) / (n+1 j).P[i])`
/// hence `Mij = (n i).(1 j-i) / (n+1 j)`
///
/// * Mij = 1 for i=j=0 and for i=degree, j=degree+1
/// * Mij = j/n+1 for i=j-1
/// * Mij = (n+1-j)/n+1 for i=j
/// * Mij = 0 otherwise
///
#[inline]
pub fn elevation_by_one_matrix_ele<N: Num>(degree: usize, i: usize, j: usize) -> (N, N) {
    let scale = N::from_usize(degree + 1).unwrap();
    if (j == 0 && i == 0) || (j == degree + 1 && i == degree) {
        (scale, scale)
    } else if i + 1 == j {
        // note j!=n+1 as if j == n+1 then i = n, and previous case is used
        // note j != 0
        (N::from_usize(j).unwrap(), scale)
    } else if i == j {
        // note j!=0 as if j == 0 then i=0 and previous case is used
        // note j != n+1 as i<=n
        (N::from_usize(degree + 1 - j).unwrap(), scale)
    } else {
        (N::zero(), scale)
    }
}

/// Generate an 'elevate by one degree' matrix given a specific degree, subject to
/// a scaling down by the result
///
/// This is a (degree+2)*(degree+1) matrix that should be applied to the control points
/// to generate a new array of control points of the elevated Bezier
#[track_caller]
#[must_use]
pub fn generate_elevate_by_one_matrix<N: Num>(matrix: &mut [N], degree: usize) -> N {
    assert!(
        matrix.len() >= (degree + 1) * (degree + 2),
        "Must have enough room in matrix for coeffs for degree -> degreee+1"
    );
    for ((i, j), m) in (0..(degree + 2))
        .flat_map(|j| (0..degree + 1).map(move |i| (i, j)))
        .zip(matrix.iter_mut())
    {
        let (n, _) = elevation_by_one_matrix_ele::<N>(degree, i, j);
        *m = n;
    }
    N::from_usize(degree + 1).unwrap()
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

//mp split_at_de_cast
/// Use de Casteljau's algorithm to split a Bernstein Bezier control
/// points set into two other Bernstein Bezier control point sets
/// at a given parameter t
///
/// The first Bezier returned has parameter t0 where 0<=t0<=1 maps to 0<=t*t0<=t
///
/// The second Bezier returned has parameter t1 where 0<=t1<=1 maps to t<=t+(1-t)*t1<=1
///
/// This destroys the provided points
pub fn split_at_de_cast<F: Num, const D: usize>(
    pts: &mut [[F; D]],
    t: F,
    b0: &mut [[F; D]],
    b1: &mut [[F; D]],
) {
    // eprintln!("Split {pts:?} at {t}");
    // Beta[0][i] = pts[i]
    // for j = 1..=n
    //  for i = 0..=n-j
    //   Beta[j][i] = (1-t)*Beta[j-1][i] + t*Beta[j-1][i+1]
    let n = pts.len();
    assert!(b0.len()>= n, "Splitting Bezier with {n} control points requires target Bezier 0 to have the same number of control points");
    assert!(b1.len()>= n, "Splitting Bezier with {n} control points requires target Bezier 1 to have the same number of control points");
    b0[0] = pts[0];
    b1[n - 1] = pts[n - 1];
    let u = F::one() - t;
    for j in 1..n {
        // For j=1, p[0] = u*p[0]+t*p[1], p[1] = u*p[1]+t*p[2], ... p[n-2] = u*p[n-2]+t*p[n-1]
        // For j=n-1, p[0] = u*p[0] + t*p[1]
        for i in 0..(n - j) {
            pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
        }
        b0[j] = pts[0];
        b1[n - 1 - j] = pts[n - 1 - j];
    }
    // eprintln!("Split into {b0:?} {b1:?}");
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
        for (c, bc) in basis_coeff_iter_num(degree, *t).enumerate() {
            matrix[r * (degree + 1) + c] = bc;
        }
    }
}
