use crate::Num;

use super::Bezier;

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
    /// Return a slice of the control points
    pub fn pts(&self) -> &[[F; D]] {
        &self.pts[0..self.degree + 1]
    }
}
