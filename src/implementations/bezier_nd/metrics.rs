use crate::Num;
use geo_nd::vector;

use super::Bezier;
use crate::BezierEval;

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
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
