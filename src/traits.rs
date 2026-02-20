use crate::{BezierBuilder, BezierError, BezierPointTIter};

/// This enumeration indicates the kind of straightening iterator over a Bezier
/// (for points or lines) that is desired; a Bezier may be straightened purely
/// by selecting a uniform `N` points along it (with parameter `0<=t<=1`), or by
/// splitting the Bezier until the difference between the line segments and a
/// straight line is within a certain closeness (squared)
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum BezierIterationType<F: Num> {
    /// Split the line until its 'closeness_sq_from_line` is within the
    /// specified value
    ClosenessSq(F),

    /// Split the line until its 'dc_sq_from_line` is within the specified value
    DcClosenessSq(F),

    /// Split the line into a specified number of sections using a uniform
    /// splitting of the parameter t between 0 and 1
    Uniform(usize),
}

/// How to reduce a Bezier curve
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum BezierReduction {
    /// Reduce using LeastSquares
    LeastSquares,

    /// Reduce using a mapping preserving points for paramater `t` uniformly from 0 to 1
    UniformT,

    /// Reduce preserving some depth (if possible, not yet implemented...)
    Preserving(u8),
}

/// Which metric to return
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum BezierMetric {
    /// The maximum of the distance squared between each point
    /// at the same t. This is a simple estimate using summation over `n` points.
    ///
    /// This is an approximation meant for testing/analysis - do not use if
    /// performance is required!
    /// Use maximim of distance squared for N uniform `t` values
    MaxDistanceSquared(usize),

    /// The integral of the distance squared between each point
    /// at the same t. This is a simple estimate using summation over `n` points.
    ///
    /// This is an approximation meant for testing/analysis - do not use if
    /// performance is required!
    SumDistanceSquared(usize),

    /// Maximum squared control point distance
    MaxControlSquared,

    /// Sum of squared control point distance
    SumControlSquared,
}

/// Trait that must be supplied to enable a value to be used as a parameter or point in a Bezier
pub trait NumOps:
    Sized + 'static + std::ops::Neg<Output = Self> + std::cmp::PartialOrd + num_traits::ConstZero + Copy
{
    /// Create a value from a fraction with signed numerator and unsigned denominator
    fn frac(n: i32, u: u32) -> Self;

    /// Create a value from an [i32]
    fn of_i32(n: i32) -> Self {
        Self::frac(n, 1)
    }

    /// Create a value from a [usize]
    fn of_usize(n: usize) -> Self {
        Self::frac(n as i32, 1)
    }

    /// Convert an f64 to this value
    fn of_f64(v: f64) -> Self;

    /// Return true if the divisor is too close to zero
    /// to provide a useful result
    fn is_unreliable_divisor(self) -> bool;

    /// Raise to the i'th power
    fn powi(self, p: i32) -> Self;

    /// Return an estimate of the square root
    fn sqrt_est(self) -> Self;

    /// Return an estimate of the cube root
    fn cbrt_est(self) -> Self;

    /// Return true if the value is negative
    fn is_sign_negative(self) -> bool {
        self < Self::ZERO
    }
}

impl NumOps for f32 {
    fn frac(n: i32, u: u32) -> Self {
        (n as f32) / (u as f32)
    }

    /// Convert an f64 to this value
    fn of_f64(v: f64) -> Self {
        v as f32
    }

    fn is_unreliable_divisor(self) -> bool {
        self.abs() <= f32::EPSILON
    }

    fn powi(self, p: i32) -> Self {
        <Self as num_traits::Float>::powi(self, p)
    }

    fn sqrt_est(self) -> Self {
        self.sqrt()
    }

    fn cbrt_est(self) -> Self {
        self.cbrt()
    }

    fn is_sign_negative(self) -> bool {
        <Self as num_traits::Float>::is_sign_negative(self)
    }
}

impl NumOps for f64 {
    fn frac(n: i32, u: u32) -> Self {
        (n as f64) / (u as f64)
    }

    fn of_f64(v: f64) -> Self {
        v
    }

    fn is_unreliable_divisor(self) -> bool {
        self.abs() <= f64::EPSILON
    }

    fn powi(self, p: i32) -> Self {
        <Self as num_traits::Float>::powi(self, p)
    }

    fn sqrt_est(self) -> Self {
        self.sqrt()
    }

    fn cbrt_est(self) -> Self {
        self.cbrt()
    }
    fn is_sign_negative(self) -> bool {
        <Self as num_traits::Float>::is_sign_negative(self)
    }
}

/// A trait that is the basic requirement for the 't' parameter for Beziers; it
/// requires basic arithmetic operations, and conversion from f32.
///
/// A blanket implementation is provided for any type that provides the required trait implementations, and hence is implemented for example by [f32] and [f64].
///
/// This is also readily supported by rational number types (with some good approximation of From f32).
pub trait Num:
    Copy
    + std::any::Any
    + PartialEq
    + PartialOrd
    + Send
    + Sync
    + std::fmt::Display
    + std::fmt::Debug
    + std::ops::Neg<Output = Self>
    + num_traits::Num
    + num_traits::NumAssignOps
    + num_traits::ConstOne
    + num_traits::ConstZero
    + num_traits::FromPrimitive
    + 'static
    + NumOps
{
    /// Return the absolute value (not all [Num] types support the [num_traits::Float] trait)
    fn nabs(self) -> Self {
        if self.is_sign_negative() {
            -self
        } else {
            self
        }
    }

    /// Return the minimum of two values
    fn min(self, other: Self) -> Self {
        if self <= other {
            self
        } else {
            other
        }
    }

    /// Return the maximum of two values
    fn max(self, other: Self) -> Self {
        if other <= self {
            self
        } else {
            other
        }
    }
}

impl<T> Num for T where
    T: Copy
        + std::any::Any
        + PartialEq
        + PartialOrd
        + Send
        + Sync
        + std::fmt::Display
        + std::fmt::Debug
        + std::ops::Neg<Output = Self>
        + num_traits::Num
        + num_traits::NumAssignOps
        + num_traits::ConstOne
        + num_traits::ConstZero
        + num_traits::FromPrimitive
        + 'static
        + NumOps
{
}

pub trait BasicBezier<F: Num, P: Clone>:
    BezierEval<F, P> + BezierOps<F, P> + BezierSplit<F> + BezierFlatIterator<F, P> + Clone
{
}

impl<F: Num, P: Clone, T: BezierEval<F, P>> BasicBezier<F, P> for T where
    T: BezierEval<F, P> + BezierOps<F, P> + BezierSplit<F> + BezierFlatIterator<F, P> + Clone
{
}

/// A trait of a Bezier that has a parameter of type 'F' and points of type 'P'
///
/// This provides for the evaluation of the Bezier, and determination of how
/// straight (or how close to a quadratic/cubic Bezier) it is.
///
/// It also provides access to the points and first derivative at a parameter 't' of the Bezier,
/// given 0<=t<=1
///
/// Precise calculation of distances and parameters of `t` relative to points are not always analytically feasible for higher degree Beziers, and so some methods return
/// Option for these precise calculations. Some additional methods are then provided for cheaper estimations, allowing (for example) a Bezier to be refined (and split, reduced etc) *only* if the point
/// *might* be closer than a specific distance
///
/// An *estimate* of the bounding box can always be obtained by: finding the endpoints of the Bezier, and finding
/// the closeness of the Bezier to a line segment; find the BBox of the endpoints, and expand by the closeness.
///
/// This is explicitly dyn-compatible, and supported by `Approximation`
pub trait BezierEval<F: Num, P: Clone> {
    /// Distance between
    ///
    /// Method provided to permit length-of-a-vector
    fn distance_sq_between(&self, p0: &P, p1: &P) -> F;

    /// Evaluate the point on the Bezier at parameter 't'
    fn point_at(&self, t: F) -> P;

    /// Evaluate the first deriative of the Bezier at parameter 't' and return the
    /// value and a scale factor (i.e. the actual deriviative is F*P)
    fn derivative_at(&self, t: F) -> (F, P);

    /// Get the endpoints of the Bezier
    #[track_caller]
    fn endpoints(&self) -> (P, P) {
        (
            self.control_points()
                .first()
                .expect("Bezier with no control points has no endpoints")
                .clone(),
            self.control_points().last().unwrap().clone(),
        )
    }

    /// Borrow the control points
    ///
    /// The endpoints are at n=0 and n=self.degree()
    fn control_points(&self) -> &[P];

    /// Return the number of control points in the Bezier
    ///
    /// The minimum number of control points is 2 (the end points); further control
    /// points increase the degree of the Bezier.
    fn num_control_points(&self) -> usize {
        self.control_points().len()
    }

    /// Returns degree of the Bezier.
    ///
    /// This is one fewer than the number of control points on the Bezier, and refers
    /// to the largest power of the parameter 't' in the Bezier evaluation.
    ///
    /// Hence cubic Beziers return 3, quadratic 2, and lines return 1
    #[inline]
    fn degree(&self) -> usize {
        self.num_control_points() - 1
    }

    /// Find how close the Bezier is to a line segment with the same endpoints
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    ///
    /// All points on the Bezier with 0<=t<=1 must lie within this closeness of the
    /// line segment between the endpoints
    ///
    /// If a Bezier is split into Beziers such that all those Beziers have a
    /// value less than `s_sq` for this closeness, then using just the endpoints of those
    /// Beziers as endpoints of straight lines will approximate the original Bezier such
    /// that the distance squared between the Bezier and the line segments is no more than
    /// the value `s_sq`
    ///
    /// This can thus be used to generate an Approximation that provides the same point locus
    /// as the original Bezier to within the closeness_sq, while losing the precise mapping from
    /// `t` to a point as the original Bezier.
    ///
    /// This value must never be greater than that returned by dc_sq_from_line
    fn closeness_sq_to_line(&self) -> F;

    /// Get the 'dc_sq' difference metric between the Bezier and a line segment with the same endpoints
    ///
    /// This is the maximum of the squared distance between the `i`th
    /// control point (i between 1 and n-1 for degree n) and the point `i/n` along
    /// the line between the endpoints
    ///
    /// If a Bezier is split into Beziers (remembering 't' values of the sub-Beziers)
    /// such that all those Beziers have a value less than `s_sq` for this closeness, then
    /// treating those Beziers as *linear* Bezier curves
    /// will provide (using the remembered t values to select the correct Bezier etc)
    /// an approximation A where the distance squared between `P(t)` on the original Bezier
    /// and the respective `A(t)` - using t linearly along the sub-Beziers - is no more than
    /// the value `s_sq`.
    ///
    /// This can thus be used to generate an Approximation that maintains the same mapping from
    /// `t` to a point as the original Bezier, to within the closeness squared
    fn dc_sq_from_line(&self) -> F;

    /// Calculate the chosen metric of this Bezier from another given by its control points
    ///
    /// If `other` is None then calculate the metric relative to a linear Bezier with the same endpoints
    ///
    /// Must return Some if `other` has the same degree as this Bezier; should return None if not
    fn metric_from(&self, other: Option<&[P]>, metric: BezierMetric) -> Option<F>;

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
    fn est_min_distance_sq_to(&self, pt: &P) -> F;

    /// Find the values of t, 0<=t<=1 and the associated point coordinate values where the point coordinate
    /// values are at a minimum and maximum (return None for min if give_min is false, and None for max if give_max is false)
    ///
    /// If *only* the minimum is required then `give_max` should be false; this might enable an optimized implementation which
    /// does not bother to calculate the maximum value. Similarly if the maximum only is required.
    ///
    /// This enables finding the bounding box of a Bezier by requesting the minimum and maximum for each
    /// index in the dimensions of the point coordinate space.
    ///
    /// The value may not be analytically calculatable, in which case None can be returned; the bounding box might then be evaluated either
    /// by iterating through points on the Bezier, or estimating with the control points perhaps expanded by the dc_sq_from_line metric
    fn t_coords_at_min_max(
        &self,
        pt_index: usize,
        give_min: bool,
        give_max: bool,
    ) -> (Option<(F, F)>, Option<(F, F)>);
}

/// A set of traits that modify the control points of a Bezier
///
/// Not supported by `Approximation`
pub trait BezierOps<F: Num, P: Clone> {
    /// Add this Bezier to another
    fn add(&mut self, other: &Self) -> bool;

    /// Subtract another Bezier from this
    fn sub(&mut self, other: &Self) -> bool;

    /// Scale the points
    fn scale(&mut self, scale: F);

    /// Map the points one at a time
    fn map_pts(&mut self, map: &dyn Fn(usize, &P) -> P);

    /// Map the points
    ///
    /// The mapping function takes a mutable slice of all the control points,
    /// and returns `true` if a modification was made
    ///
    /// Returns the value returned by `map`
    ///
    /// This cannot change the degree of a Bezier - see the [BezierMap] trait to do that
    fn map_all_pts<'a>(&'a mut self, map: &'a mut dyn FnMut(&'a mut [P]) -> bool) -> bool;
}

/// A trait provided by a Bezier to allow it to be split into two
///
/// This trait could be enhanced with a 'split_at' method, but that then
/// requires it to have a 'F:Num' generic, which means that any generic
/// type that just uses split-in-half would need F:Num, which is onerous
///
/// Hence BezierSplitAt is a separate trait
///
/// The [BezierIntoIterator] trait provides methods to iterate over a Bezier
/// curve as points or lines; it is implemented for any type that provides
/// *this* trait [BezierSplit] and [BezierEval] (plus [Clone])
///
/// This is not dyn-compatible, and not supported by `Approximation`
///
/// A trait provided by a Bezier to allow it to be split in a manner more
/// complex than just bisection
///
/// This is separate from bisection, as to split at an arbitrary point requires
/// a generic 'F'
///
/// This is not dyn-compatible, and not supported by `Approximation`
pub trait BezierSplit<F: Num>: Sized {
    /// Bisect the Bezier into two of the same degree
    fn split(&self) -> (Self, Self);
    /// Split the Bezier into two (one covering parameter values 0..t, the other t..1)
    fn split_at(&self, t: F) -> (Self, Self);

    /// Generate the Bezier with parameter t' where 0<=t'<=1 provides the same points
    /// as the original Bezier with t0<=t<=t1
    ///
    /// Note that t0 and t1 can be outside of the range `[0,1]`
    fn section(&self, t0: F, t1: F) -> Self;
}

/// Implementing BezierIntoIterator provides a Bezier with the ability
/// to be iterated as a list of points or line segments, where all the points
/// on the returned (potentially implied) lines are within a 'closeness_sq' of
/// the Bezier
///
/// This trait is frequently used to, for example, render a Bezier as lines which
/// can then be drawn on a canvas; if the closeness_sq is set to be half the width
/// of a pixel then the Bezier is rendered to within a pixel (which is the best possible approximation).
///
/// There is a *blanket* implementation of this trait for any Bezier that supports
/// the required traits of [BezierSplit], [BezierEval] and [Clone]
///
/// This is not dyn-compatible
///
/// This is supported by `Approximation`
pub trait BezierFlatIterator<F, P>: BezierEval<F, P>
where
    F: Num,
    P: Clone,
{
    /// Return an iterator of N points along the Bezier at even steps of `t`
    fn as_t_lines(&self, iter_type: BezierIterationType<F>) -> impl Iterator<Item = (F, P, F, P)>;

    /// Return an iterator of N points along the Bezier at even steps of `t`
    fn as_t_points(&self, iter_type: BezierIterationType<F>) -> impl Iterator<Item = (F, P)> {
        BezierPointTIter::new(self.as_t_lines(iter_type))
    }

    /// Calculates the length of an iterator
    fn iter_length(&self, iter_type: BezierIterationType<F>) -> F {
        self.as_t_lines(iter_type)
            .fold(F::ZERO, |acc, (_t0, p0, _t1, p1)| {
                acc + self.distance_sq_between(&p0, &p1).sqrt_est()
            })
    }

    /// Calculates the length of a setion of the Bezier
    fn iter_t_of_distance(&self, iter_type: BezierIterationType<F>, mut distance: F) -> Option<F> {
        for (t0, p0, t1, p1) in self.as_t_lines(iter_type) {
            let dp = self.distance_sq_between(&p0, &p1).sqrt_est();

            if distance <= dp {
                if dp.is_unreliable_divisor() {
                    return Some(t0);
                } else {
                    return Some(t0 + (t1 - t0) * distance / dp);
                }
            }
            distance -= dp;
        }
        None
    }
}

/// A trait provided by a Bezier to allow it to be elevated with unchanged locus
///
/// This provides methods to generate elevated Beziers by one degree, or by many;
/// the types returned for each may be different.
///
/// If an implementor of this trait does not support elevation by many degrees then
/// the associated type should be 'Self', and the elevate_by method should return None
///
/// Implementing this trait *ought* to imply support for single degree elevation, but this it not an absolute requirement.
pub trait BezierElevate<F: Num, P: Clone>: BezierEval<F, P> {
    /// Type of Bezier returned when elevated by one, or Self if that always fails
    type ElevatedByOne: BezierEval<F, P>;

    /// Elevate the Bezier by one degree; return None only if this is *never* supported
    fn elevate_by_one(&self) -> Option<Self::ElevatedByOne>;
}

/// A trait provided by a Bezier to allow it to be reduced
///
/// This provides methods to generate reduced Beziers; at a
/// minimum the type must provide a way to produce a Bezier
/// reduced by at least one degree.
pub trait BezierReduce<F: Num, P: Clone>: BezierEval<F, P> {
    /// The Bezier that this reduces to
    ///
    /// If the Bezier does not readily reduce (because it is already
    /// a linear Bezier, for example) then this can be Self
    type Reduced: BezierEval<F, P>;

    /// Determine if reduction actually reduces
    ///
    /// If this returns false then the 'reduce' method should not be invoked, and the
    /// 'closeness_sq_to_reduction' method will return a value that must be ignored
    fn can_reduce(&self, method: BezierReduction) -> bool;

    /// Find how close the Bezier is to the (potenital) reduction of the Bezier
    ///
    /// The metric should be quadratic with respect to the units of length of P,
    /// i.e. if P were in reality measured in meters, then this metric should be
    /// in units of square meteres
    ///
    /// If this returns None then 'reduce' need not return a sensible reult.
    fn dc_sq_from_reduction(&self, method: BezierReduction) -> F;

    /// Return the Bezier reduced by at least one degree
    ///
    /// The reduced Bezier *can* be reduced by more than one degree (such as
    /// a cubic Bezier from a Bezier of degree 10); the reduction *must* have a lower
    /// degree than the original.
    fn reduce(&self, method: BezierReduction) -> Option<Self::Reduced>;
}

/// A trait provided by a Bezier to allow it to be mapped to a Bezier of a different degree
///
/// This provides methods to generate a Bezier of a potenitally different degree by applying a mapping to its control points.
/// Normally this is not used to map to a Bezier of the same degree - that can be achieved by simply mapping the control points themselves.
///
/// A type may provide this, but only support mapping to one specific (other) degree of Bezier
pub trait BezierMap<F: Num, P: Clone>: BezierEval<F, P> {
    /// The Bezier that this maps to, if a mapping can succeed

    type Mapped: BezierEval<F, P>;

    /// Return the Bezier mapped by a mapping matrix to a Bezier (of a different degree, possibly)
    ///
    /// If the mapping is not supported then return None
    fn mapped_to_degree(&self, to_degree: usize, matrix: &[F]) -> Option<Self::Mapped>;

    /// Return the c_sq metric for the Bezier mapped by a mapping matrix
    ///
    /// The mapping matrix will often be
    ///
    /// If the mapping is not supported then return None
    fn dc_sq_of_mapped_from_line(&self, to_degree: usize, matrix: &[F]) -> Option<F>;
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

/// A trait that allows a type to provide building from a [BezierBuilder]
///
/// Implementing this permits a [BezierBuilder] to have its `construct` method
/// used to attempt to create a Bezier from its constraints.
pub trait BezierConstruct<F, const D: usize>: Sized
where
    F: Num,
{
    /// Create the curve with the required control points if possible
    ///
    /// The builder may be overconstrained for the type - for example, a type that permits
    /// at most a cubic Bezier cannot be built from a Builder that has five constraints
    ///
    /// The builder may be underconstrained for the type - for example a type that only
    /// provides Cubic Beziers may not permit building from a builder that defines only
    /// the endpoints of a line.
    ///
    /// A builder *can* be inconsistent - it might specify two different points for
    /// the Bezier to pass through for the same value of parameter `t`
    ///
    /// In all such cases the type can return an error; it might be that an underconstrained
    /// Bezier is built to its required degree, and then elevated to the degree of the
    /// actual type (if that is fixed), but this is not required.
    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, BezierError>;
}
