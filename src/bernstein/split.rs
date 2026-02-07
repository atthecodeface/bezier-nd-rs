use crate::{bernstein_fns, BezierSection, BezierSplit, Num};

use super::Bezier;

impl<F, const N: usize, const D: usize> BezierSplit for Bezier<F, N, D>
where
    F: Num,
{
    fn split(&self) -> (Self, Self) {
        self.split_at(0.5_f32.into())
    }
}

impl<F, const N: usize, const D: usize> BezierSection<F> for Bezier<F, N, D>
where
    F: Num,
{
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut first = *self;
        let mut latter = *self;
        bernstein_fns::de_casteljau::split_at(
            &mut latter.pts[0..self.degree + 1],
            t,
            &mut first.pts,
        );
        (first, latter)
    }
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = *self;
        if t0 > F::ZERO {
            bernstein_fns::de_casteljau::bezier_from(&mut to_split.pts, t0);
        }
        if t1 < F::ONE {
            let t10 = (t1 - t0) / (F::ONE - t0);
            bernstein_fns::de_casteljau::bezier_to(&mut to_split.pts, t10);
        }
        to_split
    }
}
