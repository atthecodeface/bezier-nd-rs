use crate::{utils, Num};
use geo_nd::vector;

use super::Bezier;
use crate::BezierEval;

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
    /// The L2 norm between two beziers given a step dt
    ///
    /// This is an approximation meant for testing/analysis - do not use if
    /// performance is required!
    ///
    /// L2 norm is the integral of the distance squared between each point
    /// at the same t. This is a simple estiamte using summation.
    pub fn metric_l2_est(&self, other: &Self, num_steps: usize) -> F {
        let mut total_d2 = F::zero();
        let ns: F = (num_steps as f32).into();
        for i in 0..num_steps {
            let t: F = (i as f32).into();
            let t = t / ns;
            let d2 = vector::distance_sq(&self.point_at(t), &other.point_at(t));
            total_d2 += d2;
        }
        total_d2 / ns
    }

    /// The maximum difference between two beziers given a step dt
    ///
    /// This is an approximation meant for testing/analysis - do not use if
    /// performance is required!
    ///
    /// Note that dm_sq >= l2, since l2, as can be intuited by it being 1/n.Sum(n)(d_sq(i/n)), for
    /// n approaching infinity, and dm_sq is Max(d_sq) = 1/n.Sum(n)(Max(d_sq)), and d_sq() at some t
    /// is always at at most Max(d_sq()).
    pub fn metric_dm_sq_est(&self, other: &Self, num_steps: usize) -> F {
        let mut max_d2 = F::zero();
        let ns: F = (num_steps as f32).into();
        for i in 0..num_steps {
            let t: F = (i as f32).into();
            let t = t / ns;
            let d2 = vector::distance_sq(&self.point_at(t), &other.point_at(t));
            if max_d2 < d2 {
                max_d2 = d2;
            }
        }
        max_d2
    }

    /// Maximum of the squared difference betwen the control points of two Beziers
    ///
    /// This metric is `Max(P0[j]-P1[j])^2 * (Sum(B[i](t)))^2`, which is
    /// the same as `Sum(B[i](t) * Max(P0[i]-P1[i]))^2`; this is more than
    /// `Sum(Max(B[i](t) * (P0[i]-P1[i])))^2`
    ///
    /// Hence this value is guaranteed to be greater than or equal to the `dm_sq` metric
    ///
    pub fn metric_dc_sq(&self, other: &Self) -> F {
        assert!(
            self.degree == other.degree,
            "Degrees of Beziers must match for a metric"
        );
        let mut max_d2 = F::zero();
        for (s, o) in self.pts.iter().zip(other.pts.iter()).take(self.degree + 1) {
            let d2 = vector::distance_sq(s, o);
            if max_d2 < d2 {
                max_d2 = d2;
            }
        }
        max_d2
    }

    /// Total of the squared difference betwen the control points of two Beziers
    ///
    /// This is clearly more than the *maximum* of the squared difference between the
    /// control points, as that is only one of the `N` control points
    ///
    /// This, though, is clearly less than or equal to `N` times the *maximum* of the squared
    /// difference between the control points.
    pub fn metric_df_sq(&self, other: &Self) -> F {
        assert!(
            self.degree == other.degree,
            "Degrees of Beziers must match for a metric"
        );
        let mut d2 = F::zero();
        for (s, o) in self.pts.iter().zip(other.pts.iter()).take(self.degree + 1) {
            d2 += vector::distance_sq(s, o);
        }
        d2
    }

    /// Apply a (degree+1) by (degree+1) matrix (should be elevate-of-reduce) to the points
    /// and calculate the new dc squared
    pub fn dc2_of_ele_red(&self, matrix: &[F]) -> F {
        assert_eq!(
            matrix.len(),
            (self.degree + 1) * (self.degree + 1),
            "Matrix to apply to Bezier of degree {} must have {} elements",
            self.degree,
            (self.degree + 1) * (self.degree + 1)
        );
        let mut max_d2 = F::zero();
        let nc = self.degree + 1;
        for (m, s) in matrix.chunks_exact(nc).zip(self.pts.iter()) {
            let mut sum = [F::zero(); D];
            for (coeff, p) in m.iter().zip(self.pts.iter()) {
                sum = vector::add(sum, p, *coeff);
            }
            let d2 = vector::distance_sq(s, &sum);
            if max_d2 < d2 {
                max_d2 = d2;
            }
        }
        max_d2
    }
}
