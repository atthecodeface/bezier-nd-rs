use geo_nd::matrix;

use crate::{utils, Num};

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
pub fn generate_elevate_by_one_matrix<N: Num>(degree: usize) -> (N, Vec<N>) {
    let scale = N::from_usize(degree + 1).unwrap();
    let matrix: Vec<_> = (0..(degree + 2))
        .flat_map(|j| (0..degree + 1).map(move |i| (i, j)))
        .map(|(i, j)| elevation_by_one_matrix_ele::<N>(degree, i, j).0)
        .collect();
    (scale, matrix)
}

/// Generate reduce-by-one matrix with minimum squared distance of all points
#[track_caller]
#[must_use]
pub fn reduce_uniform_matrix<F: Num>(reduce_from_degree: usize, degree: usize) -> Vec<F> {
    let r = reduce_from_degree + 1;
    let n = degree + 1;

    let ts: Vec<F> = utils::float_iter(n).collect();

    // bern_n is the matrix that can be applied to a Bezier of 'degree' to
    // generate the points at the values of t given in 'ts'
    //
    // bern_n is (degree+1) * (degree+1)
    let mut bern_n = vec![F::ZERO; n * n];
    crate::bernstein_fns::coeffs::generate_bernstein_matrix(&mut bern_n, degree, &ts);

    // bern_np1 is the matrix that can be applied to a Bezier of 'reduce_from_degree' to
    // generate the points at the values of t given in 'ts'
    //
    // bern_np1 is (degree+1) * (reduce_from_degree+1)
    let mut bern_np1 = vec![F::ZERO; r * n];
    crate::bernstein_fns::coeffs::generate_bernstein_matrix(&mut bern_np1, reduce_from_degree, &ts);

    // bern_n_inverse is such that bern_n_inverse * bern_n == Identity
    //
    // bern_n_inverse is (degree+1) * (degree+1)
    let mut bern_n_inverse = bern_n.clone();
    let mut lu = vec![F::ZERO; n * n];
    let mut pivot = vec![0; n];
    assert_ne!(
        matrix::lup_decompose(n, &bern_n, &mut lu, &mut pivot),
        F::ZERO,
        "Matrix is invertible"
    );
    assert!(
        matrix::lup_invert(
            n,
            &lu,
            &pivot,
            &mut bern_n_inverse,
            &mut vec![F::ZERO; n],
            &mut vec![F::ZERO; n],
        ),
        "Matrix is invertible"
    );

    // reduce = bern_n_inverse * bern_np1
    //
    // (degree+1) * (degree+1) times (degree+1) * (reduce_from_degree+1) yields...
    //
    // reduce is (degree+1) *(reduce_from_degree+1)
    let mut reduce = vec![F::ZERO; n * r];
    matrix::multiply_dyn(n, n, r, &bern_n_inverse, &bern_np1, &mut reduce);
    reduce
}

/// Generate reduce-by-one matrix with minimum squared distance of all points
pub fn reduce_l2_min_by_one_matrix<F: Num>(degree: usize) -> Vec<F> {
    let n2 = degree + 2;
    let n1 = degree + 1;

    // e is n2 rows of n1 columns
    let (scale, mut e) = generate_elevate_by_one_matrix(degree);
    for e in e.iter_mut() {
        *e /= scale;
    }

    let mut e_t = e.clone();
    geo_nd::matrix::transpose_dyn(n2, n1, &e, &mut e_t);

    // e transpose is n1 rows of n2 columns
    let mut e_e_t = vec![F::ZERO; n1 * n1];
    geo_nd::matrix::multiply_dyn(n1, n2, n1, &e_t, &e, &mut e_e_t);

    let mut e_e_t_inv = e_e_t.clone();
    let mut lu = e_e_t.clone();
    let mut pivot = vec![0; n1];
    let mut tr1 = vec![F::ZERO; n1];
    let mut tr2 = vec![F::ZERO; n1];
    let _det = geo_nd::matrix::lup_decompose(n1, &e_e_t, &mut lu, &mut pivot);

    assert!(geo_nd::matrix::lup_invert(
        n1,
        &lu,
        &pivot,
        &mut e_e_t_inv,
        &mut tr1,
        &mut tr2
    ));

    // The result reduces from n2 points to n1 points
    //
    // result is n1 rows of n2 columns
    let mut r = vec![F::ZERO; n1 * n2];
    geo_nd::matrix::multiply_dyn(n1, n1, n2, &e_e_t_inv, &e_t, &mut r);

    r
}
