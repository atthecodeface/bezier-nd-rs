use crate::constants::BINOMIALS;
use geo_nd::Float;

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
