use crate::{BezierOps, Num};
use geo_nd::vector;

use super::Bezier;

impl<F, const N: usize, const D: usize> BezierOps<F, [F; D]> for Bezier<F, N, D>
where
    F: Num,
{
    fn add(&mut self, other: &Self) -> bool {
        if self.degree != other.degree {
            false
        } else {
            for (s, o) in self.pts.iter_mut().zip(other.pts.iter()) {
                *s = vector::add(*s, o, F::ONE)
            }
            true
        }
    }

    /// Subtract another Bezier from this
    fn sub(&mut self, other: &Self) -> bool {
        if self.degree != other.degree {
            false
        } else {
            for (s, o) in self.pts.iter_mut().zip(other.pts.iter()) {
                *s = vector::sub(*s, o, F::ONE)
            }
            true
        }
    }

    /// Scale the points
    fn scale(&mut self, scale: F) {
        self.map_pts(&|_, p| vector::scale(*p, scale));
    }

    /// Map the points
    fn map_pts(&mut self, map: &dyn Fn(usize, &[F; D]) -> [F; D]) {
        for (i, p) in self.pts.iter_mut().take(self.degree + 1).enumerate() {
            *p = map(i, p);
        }
    }
}
