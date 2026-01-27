use crate::Num;
use geo_nd::vector;

use super::Bezier;
use crate::BezierEval;

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
    //mp metric_dm_est
    /// The maximum difference between two beziers given a step dt
    ///
    /// This is an approximation meant for testing/analysis - do not use if
    /// performance is required!
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

    //mp metric_df
    /// The square-root of the sum of the distance-squred between the control points of two Beziers
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

    //mp metric_dc
    /// Maximum of the difference betwen the control points of two Beziers
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

    //mp dc_of_ele_red
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
