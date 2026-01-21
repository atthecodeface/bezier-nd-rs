use crate::BezierSplit;

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
