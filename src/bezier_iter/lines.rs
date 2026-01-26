use crate::Num;
use crate::{BezierEval, BezierReduce, BezierSplit, BezierSplitIter, BezierSplitTIter};
use std::marker::PhantomData;

/// An iterator with Item = (P, P) of straight lines that form a single Bezier
///
/// An iteration will provide (Pa, Pb) pairs of points, with
/// the next iteration providing (Pb, Pc), then (Pc, Pd), etc;
/// sharing the end/start points.
///
/// This iterates over a Bezier by repeated splitting (at t=0.5) until an individual
/// segment is within the specified 'straightness'
#[derive(Clone)]
pub struct BezierLineIter<F: Num, B: BezierSplit + BezierEval<F, P> + Clone, P: Clone> {
    /// Bezier iterator
    split_iter: BezierSplitIter<B>,
    /// Maximum curviness of the line segments returned
    straightness_sq: F,
    phantom: PhantomData<P>,
}

//pi BezierLineIter
impl<F, B, P> BezierLineIter<F, B, P>
where
    F: Num,
    B: BezierSplit + BezierEval<F, P> + Clone,
    P: Clone,
{
    //fp new
    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn new(bezier: &B, straightness_sq: F) -> Self {
        let split_iter = BezierSplitIter::new(bezier);
        Self {
            straightness_sq,
            split_iter,
            phantom: PhantomData,
        }
    }

    //mp restart
    /// Clear the Bezier line iterator and restart with a new Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn restart(&mut self, bezier: &B, straightness_sq: F) {
        self.straightness_sq = straightness_sq;
        self.split_iter.restart(bezier);
    }

    //zz All done
}

//ip Iterator for BezierLineIter
impl<F, B, P> std::iter::Iterator for BezierLineIter<F, B, P>
where
    F: Num,
    B: BezierSplit + BezierEval<F, P> + Clone,
    P: Clone,
{
    /// Item is a pair of points that make a straight line
    type Item = (P, P);
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
            if bezier.closeness_sq_to_line() < self.straightness_sq {
                let (ep0, ep1) = bezier.endpoints();
                return Some((ep0.clone(), ep1.clone()));
            } else {
                self.split_iter.add_split(&bezier);
            }
        }
        None
    }

    //zz All done
}

/// An iterator with Item = (P, P) of straight lines that form a single Bezier
///
/// An iteration will provide (Pa, Pb) pairs of points, with
/// the next iteration providing (Pb, Pc), then (Pc, Pd), etc;
/// sharing the end/start points.
///
/// This iterates over a Bezier by repeated splitting (at t=0.5) until an individual
/// segment is within the specified 'straightness'
#[derive(Clone)]
pub struct BezierLineTIter<F: Num, B: BezierSplit + BezierEval<F, P> + Clone, P: Clone> {
    /// Bezier iterator
    split_iter: BezierSplitTIter<B, F>,
    /// Maximum curviness of the line segments returned
    straightness_sq: F,
    phantom: PhantomData<P>,
}

//pi BezierLineTIter
impl<F, B, P> BezierLineTIter<F, B, P>
where
    F: Num,
    B: BezierSplit + BezierEval<F, P> + Clone,
    P: Clone,
{
    //fp new
    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn new(bezier: &B, straightness_sq: F) -> Self {
        let split_iter = BezierSplitTIter::new(bezier);
        Self {
            straightness_sq,
            split_iter,
            phantom: PhantomData,
        }
    }
}

//ip Iterator for BezierLineTIter
impl<F, B, P> std::iter::Iterator for BezierLineTIter<F, B, P>
where
    F: Num,
    B: BezierSplit + BezierEval<F, P> + Clone,
    P: Clone,
{
    /// Item is a pair of points that make a straight line with the 't' values
    type Item = (F, P, F, P);

    /// next - return None or Some(ta,pa,tb,pb)
    ///
    /// It pops the first Bezier from the stack: this is (pa,px); if
    /// this is straight enough then return it, else split it in two
    /// (pa,pm), (pm,px) and push them in reverse order, then recurse.
    ///
    /// This forces the segment returned (eventually!) to be (pa,pb)
    /// and to leave the top of the stack starting with pb.
    fn next(&mut self) -> Option<Self::Item> {
        while let Some((t0, t1, bezier)) = self.split_iter.next() {
            if bezier.closeness_sq_to_line() < self.straightness_sq {
                let (ep0, ep1) = bezier.endpoints();
                return Some((t0, ep0.clone(), t1, ep1.clone()));
            } else {
                self.split_iter.add_split(&(t0, t1, bezier));
            }
        }
        None
    }

    //zz All done
}
