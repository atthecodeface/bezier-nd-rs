use crate::{bernstein_fns, BezierEval, BezierSection, BezierSplit, Num};
use geo_nd::vector;

use super::Bezier;

impl<F, const N: usize, const D: usize> BezierSplit for Bezier<F, N, D>
where
    F: crate::Num,
{
    fn split(&self) -> (Self, Self) {
        self.split_at_de_cast(0.5_f32.into())
    }
}

impl<F, const N: usize, const D: usize> BezierSection<F> for Bezier<F, N, D>
where
    F: crate::Num,
{
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut temp = *self;
        let mut s0 = *self;
        let mut s1 = *self;
        bernstein_fns::split::into_two_at_de_cast(
            &mut temp.pts[0..self.degree + 1],
            t,
            &mut s0.pts,
            &mut s1.pts,
        );
        (s0, s1)
    }
    fn section(&self, t0: F, t: F) -> Self {
        todo!();
    }
}
