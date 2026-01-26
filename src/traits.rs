/// A trait that is required for most of the 'F' parameters
///
/// This is a superset of geo_nd::Num as its requirements are needed
pub trait Num:
    Copy
    + PartialEq
    + PartialOrd
    + std::fmt::Display
    + std::fmt::Debug
    + std::ops::Neg<Output = Self>
    + num_traits::Num
    + num_traits::NumAssignOps
    + num_traits::ConstOne
    + num_traits::ConstZero
    + num_traits::FromPrimitive
    + From<f32>
    + 'static
{
    fn recip_f32(f: f32) -> Self;
    fn is_unreliable_divisor(self) -> bool;
}
use geo_nd::vector;

pub trait Float: Num + num_traits::Float + num_traits::FloatConst {}

impl<T> Float for T where T: Num + num_traits::Float + num_traits::FloatConst {}

impl<T> Num for T
where
    T: Copy
        + PartialEq
        + PartialOrd
        + std::fmt::Display
        + std::fmt::Debug
        + std::ops::Neg<Output = Self>
        + num_traits::Num
        + num_traits::NumAssignOps
        + num_traits::ConstOne
        + num_traits::ConstZero
        + num_traits::FromPrimitive
        + From<f32>
        + 'static,
{
    fn recip_f32(f: f32) -> Self {
        (1.0 / f).into()
    }
    fn is_unreliable_divisor(self) -> bool {
        self <= f32::EPSILON.into() && self >= (-f32::EPSILON).into()
    }
}

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

    /// Evaluate the deriative of the Bezier at parameter 't' and return the
    /// value and a scale factor (i.e. the actual deriviative is F*P)
    fn derivative_at(&self, t: F) -> (F, P);

    /// Borrow the endpoints of the Bezier
    fn endpoints(&self) -> (&P, &P);

    /// Find how close the Bezier is to a line segment with the same endpoints
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    ///
    /// All points on the Bezier with 0<=t<=1 must line within this closeness of the
    /// line segment between the endpoints
    fn closeness_sq_to_line(&self) -> F;

    /// Find how close the Bezier is to a Bezier of degree 2 with the same endpoints
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    fn closeness_sq_to_quadratic(&self) -> F;

    /// Find how close the Bezier to a Bezier of degree 3 with the same endpoints
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

pub trait BezierFloatArray<F: Float, const D: usize> {
    fn closeness_to_line(&self) -> F;
    fn closeness_to_quadratic(&self) -> F;
    fn closeness_to_cubic(&self) -> F;
    fn bbox_estimate(&self) -> ([F; D], [F; D]);
}

impl<B, F: Float, const D: usize> BezierFloatArray<F, D> for B
where
    B: BezierEval<F, [F; D]>,
{
    fn closeness_to_line(&self) -> F {
        self.closeness_sq_to_line().sqrt()
    }
    fn closeness_to_quadratic(&self) -> F {
        self.closeness_sq_to_quadratic().sqrt()
    }
    fn closeness_to_cubic(&self) -> F {
        self.closeness_sq_to_cubic().sqrt()
    }
    fn bbox_estimate(&self) -> ([F; D], [F; D]) {
        let pts = self.endpoints();
        let mut min = *pts.0;
        let mut max = *pts.0;
        let _ = vector::update_bbox(&[*pts.1], &mut min, &mut max);
        (min, max)
    }
}

/// A trait that provides for measurement of distance from a point to the Bezier, and closest
/// point on a Bezier to an actual point.
///
/// This is dyn-compatible.
///
/// Precise calculation is not analytically feasible for higher degree Beziers, and so has methods that return
/// Option. A method is included to help determine if a point is potentially closer than a specific
/// distance, which allows that Bezier to be refined (and split, reduced etc) *only* if the point
/// *might* b closer than a specific distance
pub trait BezierDistance<F: Num, P> {
    /// A value of t, 0<=t<=1, for which the distance between the point P and the Bezier at that parameter t
    /// (point Q) is D, and where D is the minimum distance between the point B and any point on
    /// the Bezier with 0<=t<=1.
    ///
    /// Also returns the distance squared to the Bezier with parameter value t.
    ///
    /// If this is too complicated to calculate (without iteration, for example) then None should be returned
    ///
    /// If this value returns Some((t,d_sq)), then no other point on the Bezier with t in that range
    /// has a smaller distance to the point P
    ///
    /// In other words:
    ///
    /// * if there is at least one minimum in the distance function between the Bezier and the point
    ///   P for t in the range 0<=t<=1 then this returns t at the minimum with the smallest distance
    ///   (or any one of them, if more than one are equidistant from P)
    ///
    /// * if there are no minima in the distance function between the Bezier and the point
    ///   for t in the range 0<=t<=1 then this returns t=0 or t=1 depending on which endpoint is closer to P,
    ///   (i.e. all points on the Bezier with 0<=t<=1 are at least as far as the endpoint).
    ///
    fn t_dsq_closest_to_pt(&self, pt: &P) -> Option<(F, F)>;

    /// An estimate of the minimum distance squared from the Bezier to the point
    /// for 0<=t<=1. This must ALWAYS be less than or equal to the true minimum distance.
    /// If a Bezier does not support this estimate then it *can* return ZERO.
    ///
    /// Another way to term it is the point cannot be any closer to the Bezier than this distance squared.
    ///
    /// This value *could* be determined by the distance from the point to the convex hull of a
    /// Bernstein Bezier's control points.
    ///
    /// This method can be used to determine if a point *cannot* be closer than distance D from the
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
    fn est_min_distance_sq_to(&self, _p: &P) -> F {
        F::ZERO
    }
}

/// A trait that provides for estimating and calcuating the bounding box for a Bezier.
///
/// This is dyn-compatible.
///
/// Precise calculation is not analytically feasible for higher degree Beziers, and so the method returns
/// Option.
///
/// An *estimate* of the bounding box can always be obtained by finding the endpoints of the Bezier, and finding
/// the closeness of the Bezier to a line segment; find the BBox of the endpoints, and expand by the closeness.
pub trait BezierMinMax<F: Num> {
    /// Find a value of t, 0<=t<=1 and the associated point coordinate value where the point coordinate
    /// value is the minimum (or maximum if 'use_max' is true).
    ///
    /// This enables finding the bounding box of a Bezier by requesting the minimum and maximum for each
    /// index in the dimensions of the point coordinate space.
    ///
    /// The value may not be analytically calculatable, in which case None can be returned.
    fn t_coord_at_min_max(&self, use_max: bool, pt_index: usize) -> Option<(F, F)>;
}

/// A trait provided by a Bezier to allow it to be split into two
///
/// This trait could be enhanced with a 'split_at' method, but that then
/// requires it to have a 'F:Num' generic, which means that any generic
/// type that just uses split-in-half would need F:Num, which is onerous
///
/// Hence BezierSplitAt is a separate trait
pub trait BezierSplit: Sized {
    /// Bisect the Bezier into two of the same degree
    fn split(&self) -> (Self, Self);
}

/// A trait provided by a Bezier to allow it to be split in a manner more
/// complex than just bisection
///
/// This is separate from bisection, as to split at an arbitrary point requires
/// a generic 'F'
pub trait BezierSplitExtra<F: Num>: Sized {
    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self);

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=1
    fn section(&self, t0: F, t1: F) -> Self;
}

/// A trait provided by a Bezier to allow it to be reduced
///
/// This provides methods to generate reduced Beziers; at a
/// minimum the type must provide a way to produce a Bezier
/// reduced by at least one degree.
///
///
///
///  by one degree
pub trait BezierReduce<F: Num, P: Clone>: BezierEval<F, P> {
    /// The Bezier that this reduces to
    ///
    /// If the Bezier does not readily reduce (because it is already
    /// a linear Bezier, for example) then this can be Self
    type Reduced: BezierEval<F, P>;

    /// The quadratic Bezier that this reduces to
    ///
    /// If the Bezier does not have a mechanism to reduce it
    /// to a quadratic then this can be 'Self', and 'reduced_to_quadratic' should return None
    type Quadratic: BezierEval<F, P>;

    /// The cubic Bezier that this reduces to
    ///
    /// If the Bezier does not have a mechanism to reduce it
    /// to a cubic then this can be 'Self', and 'reduced_to_cubic' should return None
    type Cubic: BezierEval<F, P>;

    /// Determine if reduction actually reduces
    ///
    /// If this returns false then the 'reduce' method should not be invoked, and the
    /// 'closeness_sq_to_reduction' method will return a value that must be ignored
    fn can_reduce() -> bool;

    /// Find how close the Bezier is to the (potenital) reduction of the Bezier
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    ///
    /// If this returns None then 'reduce' need not return a sensible reult.
    fn closeness_sq_to_reduction(&self) -> Option<F>;

    /// Return the Bezier reduced by at least one degree
    ///
    /// The reduced Bezier *can* be reduced by more than one degree (such as
    /// a cubic Bezier from a Bezier of degree 10); the reduction *must* have a lower
    /// degree than the original.
    fn reduce(&self) -> Self::Reduced;

    /// Return a Bezier of this reduced to a quadratic Bezier
    ///
    /// Invoking self.closeness_sq_to_quadratic will give a measure of how
    /// close the quadratic would be to 'self'
    ///
    /// If the Bezier does not have a mechanism to reduce it
    /// to a quadratic then this should return None
    fn reduced_to_quadratic(&self) -> Option<Self::Quadratic> {
        None
    }

    /// Return a Bezier of this reduced to a cubic Bezier
    ///
    /// Invoking self.closeness_sq_to_cubic will give a measure of how
    /// close the cubic would be to 'self'
    ///
    /// If the Bezier does not have a mechanism to reduce it
    /// to a cubic then this should return None
    fn reduced_to_cubic(&self) -> Option<Self::Cubic> {
        None
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

    /// Find how close the Bezier is to the (potenital) reduction of the Bezier
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    ///
    /// If this returns None then 'reduce' need not return a sensible reult.
    fn closeness_sq_to_reduction(&self) -> Option<F>;

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
