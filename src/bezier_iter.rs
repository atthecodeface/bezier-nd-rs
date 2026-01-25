use crate::Num;
use crate::{BezierEval, BezierSplit};
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
pub struct BezierSplitIter<B: BezierSplit + Clone> {
    /// A stack of future beziers to examine
    /// The top of the stack is p0->p1; below that is p1->p2, etc
    /// These beziers may need to be split to achieve some critertion
    stack: Vec<B>,
}

//pi BezierLineIter
impl<B> BezierSplitIter<B>
where
    B: BezierSplit + Clone,
{
    //fp new
    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn new(bezier: &B) -> Self {
        let stack = vec![bezier.clone()];
        Self { stack }
    }

    /// Clear the Bezier line iterator and restart with a new Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn restart(&mut self, bezier: &B) {
        self.stack.clear();
        self.stack.push(bezier.clone());
    }

    /// Add a Bezier to the iterator
    ///
    /// This clones the Bezier.
    pub fn add(&mut self, bezier: &B) {
        self.stack.push(bezier.clone());
    }

    /// Splits and adds both Bezier to the iterator; the first half
    /// of the split Bezier will be returned next by the iterator
    ///
    /// This clones the Bezier.
    pub fn add_split(&mut self, bezier: &B) {
        let (b0, b1) = bezier.split();
        self.stack.push(b1);
        self.stack.push(b0);
    }
}

//ip Iterator for BezierLineIter
impl<B> std::iter::Iterator for BezierSplitIter<B>
where
    B: BezierSplit + Clone,
{
    /// Item is a pair of points that make a straight line
    type Item = B;

    /// next - return the next Bezier
    ///
    /// It pops the first Bezier from the stack; a client may then
    /// consume this Bezier, or it might choose to add_split the Bezier, or whatever.
    fn next(&mut self) -> Option<Self::Item> {
        self.stack.pop()
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

/// An iterator with Item = P of points that form a single Bezier curve where the
/// steps between points would be lines that are 'straight enough'.
///
/// This iterator returns the points that BezierLineIter uses, in the
/// same order (pa, pb, ...).
///
/// The first point returned will be the starting point of the Bezier
/// (control point `p0`); the last point returned will be the end
/// point of the Bezier (control point `p1`).
pub struct BezierPointIter<P, I>
where
    P: Clone,
    I: Iterator<Item = (P, P)>,
{
    /// A line iterator that returns the next line segment required;
    /// usually the first point of this segment that this iterator
    /// provides is returned as the next point.
    ///
    /// When this returns none, the end-point of the previous
    /// iteration needs to be returned as the last point.
    lines: I,

    /// The last point to be returned - if this is valid then the line
    /// iterator has finished, and just the last point on the Bezier
    /// needs to be returned.
    last_point: Option<P>,
}

//ip BezierPointIter
impl<P, I> BezierPointIter<P, I>
where
    P: Clone,
    I: Iterator<Item = (P, P)>,
{
    //fp new
    /// Create a new point iterator from a line iterator
    pub fn new(lines: I) -> Self {
        Self {
            lines,
            last_point: None,
        }
    }

    //zz All done
}

//ii BezierPointIter
impl<P, I> std::iter::Iterator for BezierPointIter<P, I>
where
    P: Clone,
    I: Iterator<Item = (P, P)>,
{
    /// Iterator returns Point's
    type Item = P;

    /// Return the first point of any line segment provided by the
    /// line iterator, but record the endpoint of that segment first;
    /// if the line iterator has finished then return any recorded
    /// endpoint, deleting it first.
    fn next(&mut self) -> Option<Self::Item> {
        if let Some((p0, p1)) = self.lines.next() {
            self.last_point = Some(p1);
            Some(p0)
        } else {
            std::mem::take(&mut self.last_point)
        }
    }

    //zz All done
}
