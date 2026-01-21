use geo_nd::Num;

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

    /// Find how close the Bezier is to a quadratic with the same endpoints
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    fn closeness_sq_to_quad(&self) -> F;

    /// Find how close the Bezier to a cubic with the same endpoints
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    fn closeness_sq_to_cubic(&self) -> F;

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

/// A trait that provides for measurement of distance from a point to the Bezier, and closest
/// point on a Bezier to an actual point
///
/// This is not analytically feasible for higher degree Beziers, and so has methods that return
/// Option. A method is included to help determine if a point is potentially closer than a specific
/// distance, which allows that Bezier to be refined (and split, reduced etc) *only* if the point
/// *might* b closer than a specific distance
pub trait BezierDistance<F: Num, P> {
    /// The value of t (in range 0<=t<=1) of the Bezier whose point is a minimum distance from p
    ///
    /// If there is not a minimum in the distance function between the Bezier and the point
    /// for t in the range 0<=t<=1 then None should be returned.
    ///
    /// If there is only one minimum in the distance function between the Bezier and the point
    /// for t in the range 0<=t<=1 then it should be returned
    ///
    /// If there are two minima in the distance function between the Bezier and the
    /// point for t in the range 0<=t<=1 then the value with the smaller distance should be returned
    ///
    /// If there are two equal minima in the distance function between the Bezier and the
    /// point for t in the range 0<=t<=1 then either value should be returned
    ///
    /// Note that the distance between P and an endpoint (t=0 or 1) *may* be less than the Bezier point
    /// at the value of t provided by this function
    fn t_at_min_distance(&self, p: &P) -> Option<F>;

    /// The minimum distance squared from the Bezier to the point if
    /// the Bezier has a minimum in its distance function to the point provided
    /// with 0<=t<=1 (if minimum is outside 0<=t<=1 then return None)
    ///
    /// If Some(t) = self.t_at_min_distance(p) then
    /// this = distance squared between P and bezier.position_at(t)
    ///
    /// If this returns Some(f) then f must be greater than or equal to the value returned by est_min_distance_sq_to()
    ///
    /// The true minimum distance squared to the Bezier including endpoints
    /// is the minimum of (d0_sq, d1_sq, dp_sq) (d0_sq = distance squared to endpoint 0,
    /// d1_sq is distance squared to endpoint 1)
    fn min_distance_sq_to(&self, p: &P) -> Option<F>;

    /// An estimate of the minimum distance squared from the Bezier to the point
    /// for 0<=t<=1. This must ALWAYS be less than or equal to the true minimum distance.
    /// If a Bezier does not support this estimate then it *can* return ZERO.
    ///
    /// This can be used to determine if a point *cannot* be closer than distance D from the
    /// Bezier, and so the Bezier can be ignored if the point is closer to a different Bezier in a set
    /// and the 'closest distance to set of Beziers' is being calculated.
    ///
    /// The distance squared from the point to the line between the two end-points can be termed 'd_line_sq'.
    ///
    /// For a degree-1 Bezier (straight line) the value returned by this method should be d_line_sq
    ///
    /// For a Bazier that is not a line then the whole Bezier is within a distance dc_sq of
    /// the line between the two endpoints, where dc_sq is the maximum distance squared between
    /// each control point and where it *would* be on a straight Bezier of the same degree.
    /// So the value returned by this method can be (d_line_sq.sqrt() - dc_sq.sqrt())^2, or 0 if d_line_sq<dc_sq
    ///
    /// If d_line = dc+D, D>0, then:
    ///
    /// (d_line_sq.sqrt() - dc_sq.sqrt())^2 = d_line_sq + dc_sq - 2.d_line.dc
    /// = d_line_sq - dc_sq - 2D.dc < d_sq - d_line_sq
    ///
    /// So a value of d_line_sq - dc_sq can be returned, or 0. if this is negative.
    fn est_min_distance_sq_to(&self, p: &P) -> F {
        F::ZERO
    }
}

/// A dyn-compatible trait supported by a Bezier that allows it to be reduced/split
/// into potentially *different* Bezier types that themselves can be reduced/split
///
/// Implementing this requires 'alloc', but permits a Bezier to be split into
/// lines, quadratic Beziers, or cubic Beziers, to within a certain straightness
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
