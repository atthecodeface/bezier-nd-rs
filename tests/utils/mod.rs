#![allow(dead_code)]
use bezier_nd::BezierEval;
use bezier_nd::BezierSplit;

use bezier_nd::Float;
use bezier_nd::{Bezier, Num};
use geo_nd::{matrix, vector};

mod random;
pub use random::*;

/// Find the maximum distance between the Bezier for 0<=t<=1 in 'steps' intervals
pub fn max_distance_sq<
    F: Num,
    const D: usize,
    B0: BezierEval<F, [F; D]>,
    B1: BezierEval<F, [F; D]>,
>(
    b0: &B0,
    b1: &B1,
    steps: usize,
) -> F {
    let mut max_d_sq = F::ZERO;
    for t in float_iter(F::ZERO, F::ONE, steps) {
        let p0 = b0.point_at(t);
        let p1 = b1.point_at(t);
        let d_sq = vector::distance_sq(&p0, &p1);
        eprintln!("t: {t} pts: {p0:?} {p1:?} d_sq:{d_sq}");
        if d_sq > max_d_sq {
            max_d_sq = d_sq;
        }
    }
    max_d_sq
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

/// Create NPTS points on the bezier, and store them in a Vec
///
/// Split the bezier into lines given the straightness
///
/// For each line segment find NSEG points along the line and store all of these in a Vec
///
/// Find the largest closest distance between points on the bezier to segmented points
///
/// Find the largest closest distance between segmented points to points on the bezier
///
/// Both of these should be less than 'straightness'
pub struct BezierPtSet<F: Num, const D: usize> {
    pts: Vec<[F; D]>,
}

impl<F: Num, const D: usize> BezierPtSet<F, D> {
    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }
    pub fn len(&self) -> usize {
        self.pts.len()
    }
    pub fn of_point_at<B: BezierEval<F, [F; D]>>(bezier: &B, steps: usize) -> Self {
        let pts = float_iter(F::ZERO, F::ONE, steps)
            .map(|t| bezier.point_at(t))
            .collect();
        Self { pts }
    }

    pub fn of_point_iter<I: Iterator<Item = [F; D]>>(iter: I) -> Self {
        let pts = iter.collect();
        Self { pts }
    }

    pub fn of_line_iter<I: Iterator<Item = ([F; D], [F; D])>>(mut iter: I) -> Self {
        let l0 = iter.next().unwrap();
        let mut pts = vec![l0.0, l0.1];
        for (_, p) in iter {
            pts.push(p);
        }
        Self { pts }
    }

    pub fn of_line_iter_interpolated<I: Iterator<Item = ([F; D], [F; D])>>(
        iter: I,
        steps: usize,
    ) -> Self {
        let pts = iter
            .flat_map(|(p0, p1)| {
                float_iter(F::ZERO, F::ONE, steps).map(move |t| vector::mix(&p0, &p1, t))
            })
            .collect();
        Self { pts }
    }

    pub fn as_points(&self) -> impl Iterator<Item = [F; D]> + '_ {
        self.pts.iter().copied()
    }

    pub fn as_lines(&self) -> impl Iterator<Item = ([F; D], [F; D])> + '_ {
        self.pts
            .iter()
            .copied()
            .zip(self.pts.iter().skip(1).copied())
    }

    /// Calculate the maximum distance between successive points
    ///
    /// This helps to provide an error bound for a point supposedly 'on' the same Bezier
    pub fn max_pt_separation_sq(&self) -> F {
        let p0_p1_iter = self.pts.iter().zip(self.pts.iter().skip(1));
        p0_p1_iter.fold(F::ZERO, |max_d, (p0, p1)| {
            max(max_d, vector::distance_sq(p0, p1))
        })
    }

    #[track_caller]
    pub fn min_distance_sq_to_pt(&self, pt: &[F; D]) -> F {
        let mut min_d_sq = 1E10_f32.into();
        for p in &self.pts {
            let d_sq = vector::distance_sq(p, pt);
            min_d_sq = min(min_d_sq, d_sq);
        }
        min_d_sq
    }

    /// For BezierPtSet that are of the same length and corresponding points, find
    /// the maximum distance between any of the two corresponding points
    #[track_caller]
    pub fn max_distance_sq_between(&self, other: &Self) -> F {
        let mut max_d_sq = F::ZERO;
        assert_eq!(self.pts.len(), other.pts.len());
        for (p0, p1) in self.pts.iter().zip(other.pts.iter()) {
            let d_sq = vector::distance_sq(p0, p1);
            eprintln!("pts: {p0:?} {p1:?} d_sq:{d_sq}");
            if d_sq > max_d_sq {
                max_d_sq = d_sq;
            }
        }
        max_d_sq
    }

    /// For BezierPtSet for what should be (approx) the same Bezier's but with
    /// different t's for the points, find the max d_sq between *self*'s points and the closest of *other*'s points.
    ///
    /// This should be called only if other contains a full Bezier trace
    pub fn max_distance_sq_of_all_pts(&self, other: &Self) -> F {
        let mut max_excursion = F::zero();
        for p in self.pts.iter() {
            let min_d = other.pts.iter().fold(1.0E10_f32.into(), |md: F, q| {
                min(md, vector::distance_sq(p, q))
            });
            max_excursion = max(max_excursion, min_d);
        }
        max_excursion
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
