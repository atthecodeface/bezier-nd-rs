use geo_nd::{vector, Float};

use super::{bezier_fns, Bezier};

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
    //mi vector_of
    /// Returns a vector of a combination of the vectors of the bezier
    #[inline]
    pub fn vector_of(&self, sc: &[F], reduce: F) -> [F; D] {
        let mut r = [F::zero(); D];
        for (pt, sc) in self.pts.iter().zip(sc.iter()) {
            for (j, rj) in r.iter_mut().enumerate() {
                *rj += *sc * (*pt)[j];
            }
        }
        vector::reduce(r, reduce)
    }

    //mp point_at_de_cast
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

    //mp point_at
    /// Returns the point at parameter 't' along the Bezier
    pub fn point_at(&self, t: F) -> [F; D] {
        let mut r = [F::zero(); D];
        for (c, pt) in bezier_fns::basis_coeff_iter(self.degree, t).zip(self.pts.iter()) {
            r = vector::add(r, pt, c);
        }
        r
    }

    //mp nth_derivative_value_at
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
