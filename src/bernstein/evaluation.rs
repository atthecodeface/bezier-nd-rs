use crate::{BezierEval, Num};
use geo_nd::vector;

use super::{bezier_fns, Bezier};

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
    /// Returns the point at parameter 't' along the Bezier using de Casteljau's algorithm
    ///
    /// This does not require powi, but it is an O(n^2) operations and O(n) space operation
    pub fn point_at_de_cast(&self, t: F) -> [F; D] {
        // Beta[0][i] = pt[i]
        let mut pts = self.pts;
        // for j = 1..=n
        //  for i = 0..=n-j
        // Beta[j][i] = (1-t)*Beta[j-1][i] + t*Beta[j-1][i+1]
        let u = F::one() - t;
        for j in 1..(self.degree + 1) {
            for i in 0..(self.degree + 1 - j) {
                pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
            }
        }
        pts[0]
    }

    /// Returns the value of the nth deriviative at parameter 't' along the Bezier
    ///
    /// This could be optimized to not store the points, but that seems ultimately to be
    /// unnecessary
    pub fn nth_derivative_value_at(&self, n: usize, t: F) -> [F; D] {
        if n > self.degree {
            [F::zero(); D]
        } else {
            let (dn, f) = self.nth_derivative(n);
            vector::scale(dn.point_at_de_cast(t), f)
        }
    }
}

impl<F, const N: usize, const D: usize> BezierEval<F, [F; D]> for Bezier<F, N, D>
where
    F: Num,
{
    fn point_at(&self, t: F) -> [F; D] {
        let mut r = [F::zero(); D];
        for (c, pt) in bezier_fns::basis_coeff_iter_num(self.degree, t).zip(self.pts.iter()) {
            r = vector::add(r, pt, c);
        }
        r
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        (F::ONE, self.nth_derivative_value_at(1, t))
    }

    /// Borrow the endpoints of the Bezier
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self.pts[0], &self.pts[self.degree])
    }

    fn closeness_sq_to_line(&self) -> F {
        self.dc_sq_from_line()
    }
    fn dc_sq_from_line(&self) -> F {
        if self.degree < 2 {
            F::ZERO
        } else {
            let n = (self.degree + 1) as f32;
            let mut max_dc_sq = F::ZERO;
            let p01 = [self.pts[0], self.pts[self.degree]];
            for (i, p) in self.pts.iter().take(self.degree).skip(1).enumerate() {
                let t = ((i as f32) / n).into();
                let pt = vector::sum_scaled(&p01, &[F::ONE - t, t]);
                let dc_sq = vector::distance_sq(&pt, p);
                if dc_sq > max_dc_sq {
                    max_dc_sq = dc_sq;
                }
            }
            max_dc_sq
        }
    }

    fn num_control_points(&self) -> usize {
        self.degree + 1
    }

    fn control_point(&self, n: usize) -> &[F; D] {
        &self.pts[n]
    }

    fn degree(&self) -> usize {
        self.degree
    }
}
