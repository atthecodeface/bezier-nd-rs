use crate::constants::BINOMIALS;
use geo_nd::{cast, vector, Float, Num};

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
pub fn elevation_by_one_matrix_ele<N: Num + num_traits::FromPrimitive>(
    degree: usize,
    i: usize,
    j: usize,
) -> (N, N) {
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

#[track_caller]
#[must_use]
pub fn generate_elevate_by_one_matrix<N: Num + num_traits::FromPrimitive>(
    matrix: &mut [N],
    degree: usize,
) -> N {
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

//mp bernstein_split_at_de_cast
/// Use de Casteljau's algorithm to split a Bernstein Bezier control
/// points set into two other Bernstein Bezier control point sets
/// at a given parameter t
///
/// The first Bezier returned has parameter t0 where 0<=t0<=1 maps to 0<=t*t0<=t
///
/// The second Bezier returned has parameter t1 where 0<=t1<=1 maps to t<=t+(1-t)*t1<=1
///
/// This destroys the provided points
pub fn bernstein_split_at_de_cast<F: Float, const D: usize>(
    pts: &mut [[F; D]],
    t: F,
    b0: &mut [[F; D]],
    b1: &mut [[F; D]],
) {
    // Beta[0][i] = pts[i]
    // for j = 1..=n
    //  for i = 0..=n-j
    // Beta[j][i] = (1-t)*Beta[j-1][i] + t*Beta[j-1][i+1]
    let n = pts.len();
    assert!(b0.len()>= n, "Splitting Bezier with {n} control points requires target Bezier 0 to have the same number of control points");
    assert!(b1.len()>= n, "Splitting Bezier with {n} control points requires target Bezier 1 to have the same number of control points");
    let u = F::one() - t;
    for j in 1..n {
        for i in 0..(n - j) {
            pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
        }
        b0[j] = pts[0];
        b1[n - 1 - j] = pts[n - 1 - j];
    }
}

/// Elevate a Bezier by one degree, that should be reduced by 'F'
///
/// This generates and applies the elevate-by-one matrix
pub fn elevate_by_one<F: Float, const D: usize>(pts: &mut [[F; D]], ele: &mut [[F; D]]) -> F {
    assert!(
        ele.len() >= pts.len() + 1,
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
