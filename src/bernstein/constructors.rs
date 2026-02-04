use crate::Num;
use crate::{bernstein_fns, BezierBuilder, BezierConstruct};
use geo_nd::{matrix, vector};

use super::Bezier;

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
    /// Create a new Bezier given all the control points
    pub fn new(pts: &[[F; D]]) -> Self {
        let n = pts.len();
        assert!(n >= 1, "Beziers have at least 1 point");
        assert!(
            n <= N,
            "Attempt to create a Bezier of max {N} pts with actually {n} pts"
        );
        let mut s = Self {
            degree: n - 1,
            ..Default::default()
        };
        s.pts.split_at_mut(n).0.copy_from_slice(pts);
        s
    }

    /// Create a new Bezier that is the nth derivative of this
    /// subject to a scaling factor (i.e. the actual Bezier should be scaled up by F)
    ///
    /// As this type uses Bernstein polynomials, this is the Bezier of degree N-n
    /// that has points control points:
    ///
    /// `Pn[i] = (N n) . Sum((-1)^j * P[i+n-j] * (n j))`
    ///
    /// (n j) is kept in `bezier_nd::constants::BINOMIALS[n][j+1]`
    pub fn nth_derivative(&self, n: usize) -> (Self, F) {
        let mut s = Self {
            degree: self.degree - n,
            ..Default::default()
        };
        s.degree = self.degree - n;
        let scale = bernstein_fns::nth_derivative(&self.pts[0..self.degree + 1], n, &mut s.pts);
        (s, scale)
    }

    /// Apply a (new_degree+1) by (degree+1) matrix to the points to generate a new Bezier
    /// of a new degree
    pub fn apply_matrix(&self, matrix: &[F], new_degree: usize) -> Self {
        assert!(
            new_degree < N + 1,
            "Cannot create a Bezier<{N}> of {new_degree}",
        );
        assert_eq!(
            matrix.len(),
            (new_degree + 1) * (self.degree + 1),
            "Matrix to apply to Bezier of degree {} to degree {new_degree} must have {} elements",
            self.degree,
            (new_degree + 1) * (self.degree + 1)
        );
        let _nr = new_degree + 1;
        let nc = self.degree + 1;
        let mut s = Self {
            degree: new_degree,
            ..Default::default()
        };
        for (m, sp) in matrix.chunks_exact(nc).zip(s.pts.iter_mut()) {
            let mut sum = [F::zero(); D];
            for (coeff, p) in m.iter().zip(self.pts.iter()) {
                sum = vector::add(sum, p, *coeff);
            }
            *sp = sum;
        }
        s
    }

    /// Elevate a Bezier by one degree
    ///
    /// This generates and applies the elevate-by-one matrix
    pub fn elevate_by_one(&self) -> Self {
        assert!(
            self.degree < N + 2,
            "Cannot elevate Bezier<{N}> which already has {} pts",
            self.degree + 1
        );
        // n = number of points, i.e. self.pts[n-1] is the last valid point
        let n = self.degree + 1;
        let mut s = Self {
            degree: n,
            ..Default::default()
        };
        s.pts[0] = self.pts[0];
        s.pts[self.degree] = self.pts[self.degree];
        let scale: F = (1.0 / ((n + 1) as f32)).into();
        for j in 1..self.degree {
            s.pts[j] = vector::scale(
                vector::add(
                    vector::scale(self.pts[j], (j as f32).into()),
                    &self.pts[j + 1],
                    ((n + 1 - j) as f32).into(),
                ),
                scale,
            );
        }
        s
    }
}

impl<F: Num, const N: usize, const D: usize> BezierConstruct<F, D> for Bezier<F, N, D> {
    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, ()> {
        let (mut matrix, pts) = builder.get_matrix_pts()?;
        if pts.len() > N {
            return Err(());
        }
        let degree = pts.len() - 1;
        let bezier = Self::new(&pts);

        let mut lu = matrix.clone();
        let mut pivot = vec![0; degree + 1];
        let mut tr0 = vec![F::zero(); degree + 1];
        let mut tr1 = vec![F::zero(); degree + 1];

        if matrix::lup_decompose(degree + 1, &matrix, &mut lu, &mut pivot) == F::ZERO {
            return Err(());
        }

        assert!(
            matrix::lup_invert(degree + 1, &lu, &pivot, &mut matrix, &mut tr0, &mut tr1),
            "Matrix must be invertible"
        );
        Ok(bezier.apply_matrix(&matrix, degree))
    }
}
