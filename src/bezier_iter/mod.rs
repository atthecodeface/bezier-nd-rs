mod beziers;
mod lines;
mod points;
mod split;

pub use beziers::BezierQuadIter;
pub use lines::{BezierLineIter, BezierLineTIter};
pub use points::{BezierPointIter, BezierPointTIter};
pub use split::{BezierSplitIter, BezierSplitTIter};

pub trait BezierIntoIterator<F, B, P>
where
    F: crate::Num,
    B: crate::BezierSplit + crate::BezierEval<F, P> + Clone,
    P: Clone,
{
    fn as_lines(&self, closeness_sq: F) -> impl Iterator<Item = (P, P)>;
    fn as_points(&self, closeness_sq: F) -> impl Iterator<Item = P> {
        BezierPointIter::new(self.as_lines(closeness_sq))
    }
    fn as_t_lines(&self, closeness_sq: F) -> impl Iterator<Item = (F, P, F, P)>;
    fn as_t_points(&self, closeness_sq: F) -> impl Iterator<Item = (F, P)> {
        BezierPointTIter::new(self.as_t_lines(closeness_sq))
    }
}
impl<F, B, P> BezierIntoIterator<F, B, P> for B
where
    F: crate::Num,
    B: crate::BezierSplit + crate::BezierEval<F, P> + Clone,
    P: Clone,
{
    fn as_lines(&self, closeness_sq: F) -> impl Iterator<Item = (P, P)> {
        BezierLineIter::new(self, closeness_sq)
    }
    fn as_t_lines(&self, closeness_sq: F) -> impl Iterator<Item = (F, P, F, P)> {
        BezierLineTIter::new(self, closeness_sq)
    }
}
