use crate::Num;
use crate::{BezierEval, BezierIterationType, BezierSplit, BezierSplitTIter};
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
pub struct BezierLineTIter<F: Num, B: BezierSplit<F> + BezierEval<F, P> + Clone, P: Clone> {
    /// Bezier iterator
    split_iter: BezierSplitTIter<B, F>,
    /// Maximum curviness of the line segments returned
    straightness_sq: F,
    metric: for<'a> fn(&'a B) -> F,
    phantom: PhantomData<P>,
}

//pi BezierLineTIter
impl<F, B, P> BezierLineTIter<F, B, P>
where
    F: Num,
    B: BezierSplit<F> + BezierEval<F, P> + Clone,
    P: Clone,
{
    /// Create a [BezierLineTIter] from a Bezier, using a specified [BezierIterationType]
    ///
    /// This might mean creating an iterator that uses uniform `t` values
    /// between 0 and 1, or by splitting the Bezier curve into halves
    /// recursively until a specific 'closeness' to a line is achieved
    pub fn of_iter_type(bezier: &B, iter_type: BezierIterationType<F>) -> Self {
        match iter_type {
            BezierIterationType::ClosenessSq(f) => Self::split_bezier(bezier, f, false),
            BezierIterationType::DcClosenessSq(f) => Self::split_bezier(bezier, f, true),
            BezierIterationType::Uniform(n) => Self::uniform(bezier, F::ZERO, F::ONE, n),
        }
    }

    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn split_bezier(bezier: &B, straightness_sq: F, use_dc: bool) -> Self {
        let split_iter = BezierSplitTIter::split_bezier(bezier);
        let metric = if use_dc {
            B::dc_sq_from_line
        } else {
            B::closeness_sq_to_line
        };

        Self {
            straightness_sq,
            split_iter,
            metric,
            phantom: PhantomData,
        }
    }

    /// Create a uniform line iterator with `t0<=t<=t1` with the specified number of steps for the Bezier
    pub fn uniform(bezier: &B, t0: F, t1: F, num_steps: usize) -> Self {
        let split_iter = BezierSplitTIter::uniform(bezier, t0, t1, num_steps);
        Self {
            straightness_sq: F::ZERO,
            split_iter,
            metric: B::closeness_sq_to_line,
            phantom: PhantomData,
        }
    }
}

//ip Iterator for BezierLineTIter
impl<F, B, P> std::iter::Iterator for BezierLineTIter<F, B, P>
where
    F: Num,
    B: BezierSplit<F> + BezierEval<F, P> + Clone,
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
        if self.straightness_sq == F::ZERO {
            self.split_iter
                .next()
                .map(|(t0, t1, b)| (t0, b.point_at(t0), t1, b.point_at(t1)))
        } else {
            while let Some((t0, t1, bezier)) = self.split_iter.next() {
                let metric = (self.metric)(&bezier);
                if metric < self.straightness_sq {
                    let (ep0, ep1) = bezier.endpoints();
                    return Some((t0, ep0.clone(), t1, ep1.clone()));
                } else {
                    self.split_iter.add_split(&(t0, t1, bezier));
                }
            }
            None
        }
    }
}
