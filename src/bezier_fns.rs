use crate::constants::BINOMIALS;
use geo_nd::{vector, Float};

/// Calculate the ith Bernstein polynomial coefficient at 't' for a given degree.
///
/// This is (1-t)^(degree-i) * t^i * (degree! / (i!.(degree-i)!) )
#[inline]
pub fn bernstein_basis_coeff<F: Float>(degree: usize, i: usize, t: F) -> F {
    let u = F::one() - t;
    let coeffs = BINOMIALS[degree];
    t.powi(i as i32) * u.powi((degree - i) as i32) * (coeffs[1 + i]).into()
}

/// Calculate the Mij element of the elevation matrix to elevate by one degree
/// given the current degree. The elevation matrix is (degree+2)x(degree+1)
///
/// This is given by P'[j] = Sum((n i).(1 j-i) / (n+1 j).P[i])
/// hence Mij = (n i).(1 j-i) / (n+1 j)
///
/// * Mij = 1 for i=j=0 and for i=degree, j=degree+1
/// * Mij = j/n+1 for i=j-1
/// * Mij = (n+1-j)/n+1 for i=j
/// * Mij = 0 otherwise
///
#[inline]
pub fn elevation_by_one_matrix_ele<F: Float>(n: usize, i: usize, j: usize) -> F {
    if (j == 0 && i == 0) || (j == n + 1 && i == n) {
        F::one()
    } else if i + 1 == j {
        // note j!=n+1 as if j == n+1 then i = n, and previous case is used
        // note j != 0
        ((j as f32) / (n + 1) as f32).into()
    } else if i == j {
        // note j!=0 as if j == 0 then i=0 and previous case is used
        // note j != n+1 as i<=n
        (((n + 1 - j) as f32) / ((n + 1) as f32)).into()
    } else {
        F::zero()
    }
}

/// Calculate the derivative Bernstein points of the nth derivative of a Bernstein Bezier
///
/// Updates the provided d_pts slice, and returns the scaling which should be applied (multiply by F)
#[track_caller]
pub fn nth_bernstein_derivative<F: Float, const D: usize>(
    pts: &[[F; D]],
    n: usize,
    d_pts: &mut [[F; D]],
) -> F {
    assert!(
        pts.len() >= n + 1,
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
        scale *= ((degree - i) as f32).into();
    }
    for (i, sp) in d_pts.iter_mut().take(degree + 1 - n).enumerate() {
        let mut m1_n_positive = (n & 1) == 0;
        for (c, p) in BINOMIALS[n][1..].iter().zip(pts[i..].iter()) {
            if m1_n_positive {
                *sp = vector::add(*sp, p, (*c).into());
            } else {
                *sp = vector::sub(*sp, p, (*c).into());
            }
            m1_n_positive = !m1_n_positive;
        }
    }
    scale
}
