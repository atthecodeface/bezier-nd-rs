use crate::BezierSplitIter;
use crate::Num;
use crate::{BezierEval, BezierReduce, BezierSplit};
use std::marker::PhantomData;

/// An iterator with Item = (P, P) of straight lines that form a single Bezier
#[derive(Clone)]
pub struct BezierQuadIter<
    F: Num,
    B: BezierSplit + BezierEval<F, P> + Clone + BezierReduce<F, P>,
    P: Clone,
> {
    /// Bezier iterator
    split_iter: BezierSplitIter<B>,
    /// Maximum curviness of the line segments returned
    closeness_sq: F,
    phantom: PhantomData<P>,
}

impl<F, B, P> BezierQuadIter<F, B, P>
where
    F: Num,
    B: BezierSplit + BezierEval<F, P> + Clone + BezierReduce<F, P>,
    P: Clone,
{
    //fp new
    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn new(bezier: &B, closeness_sq: F) -> Self {
        let split_iter = BezierSplitIter::new(bezier);
        Self {
            closeness_sq,
            split_iter,
            phantom: PhantomData,
        }
    }

    //zz All done
}

//ip Iterator for BezierLineIter
impl<F, B, P> std::iter::Iterator for BezierQuadIter<F, B, P>
where
    F: Num,
    B: BezierSplit + BezierEval<F, P> + Clone + BezierReduce<F, P>,
    P: Clone,
{
    /// Item is a pair of points that make a straight line
    type Item = <B as BezierReduce<F, P>>::Quadratic;
    /// next - return None or Some(pa,pb)
    ///
    /// It pops the first Bezier from the stack: this is (pa,px); if
    /// this is straight enough then return it, else split it in two
    /// (pa,pm), (pm,px) and push them in reverse order, then recurse.
    ///
    /// This forces the segment returned (eventually!) to be (pa,pb)
    /// and to leave the top of the stack starting with pb.
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(bezier) = self.split_iter.next() {
            if bezier.closeness_sq_to_quadratic() < self.closeness_sq {
                return bezier.reduced_to_quadratic();
            } else {
                self.split_iter.add_split(&bezier);
            }
        }
        None
    }

    //zz All done
}
