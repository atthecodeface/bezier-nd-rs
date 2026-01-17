use geo_nd::Float;

use super::Bezier;

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
    //ap degree
    /// Return degree of the Bezier (e.g. 3 for a cubic)
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Return a slice of the control points
    pub fn pts(&self) -> &[[F; D]] {
        &self.pts[0..self.degree + 1]
    }
}
