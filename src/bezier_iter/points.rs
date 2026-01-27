
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

impl<P, I> std::convert::From<I> for BezierPointIter<P, I>
where
    P: Clone,
    I: Iterator<Item = (P, P)>,
{
    fn from(lines: I) -> Self {
        Self::new(lines)
    }
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

/// An iterator with Item = P of points that form a single Bezier curve where the
/// steps between points would be lines that are 'straight enough'.
///
/// This iterator returns the points that BezierLineIter uses, in the
/// same order (pa, pb, ...).
///
/// The first point returned will be the starting point of the Bezier
/// (control point `p0`); the last point returned will be the end
/// point of the Bezier (control point `p1`).
pub struct BezierPointTIter<F, P, I>
where
    P: Clone,
    I: Iterator<Item = (F, P, F, P)>,
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
    last_point: Option<(F, P)>,
}

//ip BezierPointIter
impl<F, P, I> BezierPointTIter<F, P, I>
where
    P: Clone,
    I: Iterator<Item = (F, P, F, P)>,
{
    //fp new
    /// Create a new point iterator from a line iterator
    pub fn new(lines: I) -> Self {
        Self {
            lines,
            last_point: None,
        }
    }
}

impl<F, P, I> std::convert::From<I> for BezierPointTIter<F, P, I>
where
    P: Clone,
    I: Iterator<Item = (F, P, F, P)>,
{
    fn from(lines: I) -> Self {
        Self::new(lines)
    }
}

//ii BezierPointIter
impl<F, P, I> std::iter::Iterator for BezierPointTIter<F, P, I>
where
    P: Clone,
    I: Iterator<Item = (F, P, F, P)>,
{
    /// Iterator returns Point's
    type Item = (F, P);

    /// Return the first point of any line segment provided by the
    /// line iterator, but record the endpoint of that segment first;
    /// if the line iterator has finished then return any recorded
    /// endpoint, deleting it first.
    fn next(&mut self) -> Option<Self::Item> {
        if let Some((t0, p0, t1, p1)) = self.lines.next() {
            self.last_point = Some((t1, p1));
            Some((t0, p0))
        } else {
            std::mem::take(&mut self.last_point)
        }
    }

    //zz All done
}
