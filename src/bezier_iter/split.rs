use crate::BezierSplit;
use crate::Num;

/// An iterator with Item = Bezier, allowing for splitting dynamically during
/// iteration
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

/// An iterator with Item = (F, F, Bezier), allowing for splitting dynamically during
/// iteration
#[derive(Clone)]
pub struct BezierSplitTIter<B: BezierSplit + Clone, F: Num> {
    /// A stack of future beziers to examine
    /// The top of the stack is p0->p1; below that is p1->p2, etc
    /// These beziers may need to be split to achieve some critertion
    stack: Vec<(F, F, B)>,
}

//pi BezierLineIter
impl<B, F> BezierSplitTIter<B, F>
where
    B: BezierSplit + Clone,
    F: Num,
{
    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn new(bezier: &B) -> Self {
        let stack = vec![(F::ZERO, F::ONE, bezier.clone())];
        Self { stack }
    }

    /// Add a Bezier to the iterator
    ///
    /// This clones the Bezier.
    pub fn add(&mut self, ts: (F, F), bezier: &B) {
        self.stack.push((ts.0, ts.1, bezier.clone()));
    }

    /// Splits and adds both Bezier to the iterator; the first half
    /// of the split Bezier will be returned next by the iterator
    pub fn add_split(&mut self, (t0, t1, bezier): &(F, F, B)) {
        let tm = (*t0 + *t1) / (2.0_f32).into();
        let (b0, b1) = bezier.split();
        self.stack.push((tm, *t1, b1));
        self.stack.push((*t0, tm, b0));
    }
}

//ip Iterator for BezierLineIter
impl<B, F> std::iter::Iterator for BezierSplitTIter<B, F>
where
    B: BezierSplit + Clone,
    F: Num,
{
    /// Item is a pair of points that make a straight line
    type Item = (F, F, B);

    /// next - return the next Bezier
    ///
    /// It pops the first Bezier from the stack; a client may then
    /// consume this Bezier, or it might choose to add_split the Bezier, or whatever.
    fn next(&mut self) -> Option<Self::Item> {
        self.stack.pop()
    }
}
