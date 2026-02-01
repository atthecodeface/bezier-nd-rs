use crate::BezierBuilder;
use crate::Num;
use geo_nd::{matrix, vector};

use super::{bezier_fns, Bezier};

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
        let scale = bezier_fns::nth_derivative(&self.pts[0..self.degree + 1], n, &mut s.pts);
        (s, scale)
    }

    /// Use de Casteljau's algorithm to split
    pub fn split_at_de_cast(mut self, t: F) -> (Self, Self) {
        let mut s0 = self;
        let mut s1 = self;
        bezier_fns::split_at_de_cast(
            &mut self.pts[0..self.degree + 1],
            t,
            &mut s0.pts,
            &mut s1.pts,
        );
        (s0, s1)
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

    /// Returns two Bezier's that split the curve at parameter t=0.5
    pub fn bisect(&self) -> (Self, Self) {
        self.split_at_de_cast(0.5_f32.into())
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

    //mp bezier_between
    /// Returns the Bezier that is a subset of this Bezier between two parameters 0 <= t0 < t1 <= 1
    pub fn bezier_between(&self, t0: F, t1: F) -> Self {
        let dt = t1 - t0;
        assert!(t0 < F::one(), "Must select a t0 that is less than 1.0");
        assert!(dt > F::zero(), "Must select a t range that is > 0");
        if t0 == F::zero() {
            self.split_at_de_cast(dt).0
        } else if t1 == F::one() {
            self.split_at_de_cast(dt).1
        } else {
            // b is the subset from t0 to 1.0
            let (_, b) = self.split_at_de_cast(t0);
            b.split_at_de_cast(dt / (F::one() - t0)).0
        }
    }

    /// Construct a [Bezier] from a builder, of the minimum degree
    ///
    /// Return Err if the builder has a degree larger than this type
    /// permits
    pub fn of_builder(builder: BezierBuilder<F, D>) -> Result<Self, ()> {
        let degree = builder.bezier_min_degree();
        if degree > Self::max_degree() {
            return Err(());
        }
        let n2 = (degree + 1) * (degree + 1);
        let mut bern_n = [F::zero(); 100];

        let mut basis = vec![];
        let mut pts = vec![];
        for c in builder.iter() {
            let t = c.at();
            let pt = c.posn();
            pts.push(*pt);
            bezier_fns::generate_bernstein_matrix(&mut bern_n[0..(degree + 1)], degree, &[t]);
            basis.extend(bern_n.iter().take(degree + 1));
        }

        if pts.len() != degree + 1 {
            return Err(());
        }

        let bezier = Self::new(&pts);

        let mut basis_inverse = basis.clone();
        let mut lu = basis.clone();
        let mut pivot = vec![0; degree + 1];
        if matrix::lup_decompose(degree + 1, &basis[0..n2], &mut lu[0..n2], &mut pivot) == F::ZERO {
            return Err(());
        }

        let mut tr0 = vec![F::zero(); degree + 1];
        let mut tr1 = vec![F::zero(); degree + 1];
        assert!(
            matrix::lup_invert(
                degree + 1,
                &lu,
                &pivot,
                &mut basis_inverse,
                &mut tr0,
                &mut tr1
            ),
            "Matrix is invertible"
        );
        Ok(bezier.apply_matrix(&basis_inverse, degree))
    }
}
