use crate::Num;
use geo_nd::vector;

use super::Bezier;

impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
    //mp map_pts
    /// Apply a function to all of the points in the Bezier
    pub fn map_pts<Map: Fn([F; D]) -> [F; D]>(&mut self, map: Map) {
        for p in self.pts.iter_mut().take(self.degree + 1) {
            *p = map(*p);
        }
    }

    //mp scale
    /// Scale the Bezier by applying the scale factor to all of the points
    ///
    /// This is an example of the [Bezier::map_pts] method
    pub fn scale(&mut self, s: F) {
        self.map_pts(|p| vector::scale(p, s));
    }
}
