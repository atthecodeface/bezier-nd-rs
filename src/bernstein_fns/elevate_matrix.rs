use crate::constants::BINOMIALS_U;
use crate::{Float, Num};
use geo_nd::vector;

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
