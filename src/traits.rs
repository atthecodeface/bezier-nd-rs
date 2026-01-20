use geo_nd::{vector, Num};

/// A trait of a Bezier that has a parameter of type 'F' and points of type 'P'
///
/// This provides for the evaluation of the Bezier, and determination of how
/// straight (or how close to a quadratic/cubic Bezier) it is.
///
/// It also provides access to the points
///
/// This is dyn-compatible
pub trait BezierEval<F: Num, P: Clone> {
    /// Evaluate the Bezier at parameter 't' and return the point P
    fn point_at(&self, t: F) -> P;

    /// Evaluate the deriative of the Bezier at parameter 't' and return the point P
    fn derivative_at(&self, t: F) -> P;

    /// Borrow the endpoints of the Bezier
    fn endpoints(&self) -> (&P, &P);

    /// Return true if the Bezier is within 'straightness' of a straight line
    fn is_straight(&self, straightness: F) -> bool;

    /// Find the closess of the Bezier to a quadratic with the same endpoints
    fn closeness_to_quad(&self) -> F;

    /// Find the closess of the Bezier to a cubic with the same endpoints
    fn closeness_to_cubic(&self) -> F;

    /// Return the number of control points in the Bezier
    fn num_control_points(&self) -> usize;

    /// Borrow the nth control point (in 0..num_control_points())
    fn control_point(&self, n: usize) -> &P;

    /// Returns degree of the Bezier
    ///
    /// Cubic beziers return 3
    /// Quadratic beziers return 2
    /// Linear beziers (lines...) return 1
    #[inline]
    fn degree(&self) -> usize {
        self.num_control_points() - 1
    }
}

pub trait BoxedBezier<F: Num, P: Clone>: BezierEval<F, P> {
    /// Optionally reduce the bezier by one degree in some manner
    fn boxed_reduce(&self) -> Option<Box<dyn BoxedBezier<F, P>>> {
        None
    }

    /// Optionally split the bezier into two at t=0.5
    fn boxed_split(&self) -> Option<(Box<dyn BoxedBezier<F, P>>, Box<dyn BoxedBezier<F, P>>)> {
        None
    }

    /// Optionally reduce the Bezier down to a quadratic Bezier
    ///
    /// If bezier.closeness_to_quad() returns 'straightness' then the
    /// Bezier returned by this call should be within 'straightness' for all
    /// parameter values t in 0 <= t <= 1.0
    ///
    /// A Bezier of degree 2 or lower (a quadratic or linear Bezier) should return None
    fn boxed_reduce_to_quad(&self) -> Option<Box<dyn BoxedBezier<F, P>>> {
        None
    }

    /// Optionally reduce the Bezier down to a cubic Bezier
    ///
    /// If bezier.closeness_to_cubic() returns 'straightness' then the
    /// Bezier returned by this call should be within 'straightness' for all
    /// parameter values t in 0 <= t <= 1.0
    ///
    /// A Bezier of degree 3 or lower (a cubic, quadratic or linear Bezier) should return None
    fn boxed_reduce_to_cubic(&self) -> Option<Box<dyn BoxedBezier<F, P>>> {
        None
    }
}

/// A trait provided by a Bezier to allow it to be split
pub trait BezierSplit: Sized {
    /// Bisect the Bezier into two of the same degree
    fn split(&self) -> (Self, Self);
}

/// A trait provided by a Bezier to allow it to be reduced by one degree
pub trait BezierReduce<F: Num, P: Clone> {
    /// The Bezier that this reduces to
    type Reduced: BezierEval<F, P>;

    /// Return a Bezier reduced by one degree
    fn reduce(&self) -> Self::Reduced;
}
