use std::num;

//a Imports
use geo_nd::vector;
use geo_nd::Float;

use crate::{BezierLineIter, BezierPointIter};

const BINOMIALS: &[&[f32]] = &[
    &[1., 1.],
    &[2., 1., 1.],
    &[4., 1., 2., 1.],
    &[8., 1., 3., 3., 1.],
    &[16., 1., 4., 6., 4., 1.],
    &[32., 1., 5., 10., 10., 5., 1.],
    &[64., 1., 6., 15., 20., 15., 6., 1.],
    &[128., 1., 7., 21., 35., 35., 21., 7., 1.],
    &[256., 1., 8., 28., 56., 70., 56., 28., 8., 1.],
    &[512., 1., 9., 36., 84., 126., 126., 84., 36., 9., 1.],
];

const BINOMIAL_DIFFS: &[&[f32]] = &[
    &[1., 0.],
    &[1., -1., 1.],
    &[1., 1., -1., 1.],
    &[8., 1., 3., 3., 1.],
    &[16., 1., 4., 6., 4., 1.],
    &[32., 1., 5., 10., 10., 5., 1.],
    &[64., 1., 6., 15., 20., 15., 6., 1.],
    &[128., 1., 7., 21., 35., 35., 21., 7., 1.],
    &[256., 1., 8., 28., 56., 70., 56., 28., 8., 1.],
    &[512., 1., 9., 36., 84., 126., 126., 84., 36., 9., 1.],
];

pub const ELEVATE_BY_ONE_MATRICES: &[&[f32]] = &[
    &[1., 1.],
    &[1., 0., 0.5, 0.5, 0., 1.],
    &[
        1.0, 0.0, 0.0, 0.33333334, 0.6666667, 0.0, 0.0, 0.6666667, 0.33333334, 0.0, 0.0, 1.0,
    ],
    &[
        1., 0., 0., 0., 0.25, 0.75, 0., 0., 0., 0.5, 0.5, 0., 0., 0., 0.75, 0.25, 0., 0., 0., 1.,
    ],
    &[
        1.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.8, 0.0, 0.0, 0.0, 0.0, 0.4, 0.6, 0.0, 0.0, 0.0, 0.0, 0.6,
        0.4, 0.0, 0.0, 0.0, 0.0, 0.8, 0.2, 0.0, 0.0, 0.0, 0.0, 1.0,
    ],
    &[
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.16666667, 0.8333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.33333334,
        0.6666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6666667,
        0.33333334, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8333333, 0.16666667, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ],
    &[
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.14285715, 0.85714287, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.2857143, 0.71428573, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42857143, 0.5714286, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.5714286, 0.42857143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.71428573, 0.2857143,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.85714287, 0.14285715, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ],
    &[
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.25, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.625, 0.375, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.75, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.875, 0.125, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ],
    &[
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.11111111, 0.8888889, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.22222222, 0.7777778, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.33333334, 0.6666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.44444445, 0.5555556, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5555556, 0.44444445, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.6666667, 0.33333334, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7777778,
        0.22222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8888889, 0.11111111, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ],
];

//a Bezier
//tp Bezier
/// A [Bezier] is an implementation of an arbitrary bezier
///
/// The control points are stored to be used with Bernstein basis polynomials
///
/// ``` ignore
///    b[i,n](t) = n!/(i!(n-i)!) . u^i . t^(n-i) where u=1-t, n=Bezier degree, 0<=i <=n
/// ```
/// A point at parameter t (0 <= t <= 1) is then Sum(b[i,n](t).pts[i])
///
/// If n were two then this is u^2.pts[0] + 2.u.t.pts[1] + t^2.pts[2], which can be rewritten as:
/// ``` ignore
///   (1-2t+t^2).pts[0] + 2.(1-t).t.pts[1] + t^2.pts[2]
///   = t^2.(pts[0] - 2.pts[1] + pts[2]) + 2t.(pts[1] - pts[0]) + pts[0]
/// ```
///
/// i.e. a set of different point values where a monomial basis (1, t, t^2, ...) is used
///
/// To convert from Bernstein basis to monomial basis for degree 4 (five points) the pts would map
/// using a 5x5 matrix:
///
/// ``` ignore
///     |  1   0   0   0   0 |
///     | -4   4   0   0   0 |
///     |  6 -12   6   0   0 |
///     | -4  12 -12   4   0 |
///     |  1  -4   6  -4   1 |
///```
///
/// Actual value Mij for degree N is (-1)^(j-i) . (n C j) . (j c i) for i<=j, 0 otherwise (by inspection...)
///
/// Note that this matrix is guaranteed to be invertible, hence all Bernstein representations are equivalent
/// to monomial representations, and vice versa
///
/// Note that P[N] is only used in conjunction with monomial t^N, and hence if this were to be used
/// to derive points etc then it would be unstable near t=1. However, it has the potentially useful
/// property (for t<1) that 1 > t > t^2 > t^3 > 0.
///
/// Another possible basis would be for the basis 1, (t-T), (t-T)^2, etc for a fixed T (effectively
/// monomial is this with T=0). This is termed the Taylor basis vectors in https://arxiv.org/pdf/2201.07834
/// This has the property that it is locally (at t=0) approximatable (as one can ignore high powers of (t-T)).
///
/// Metrics ('distance' between two Beziers)
///
/// Various metrics (real number) can be used for 'measuring' the difference between two Beziers. Note that a metric
/// is a special kind of measurement (such as the distance between two points), where:
///
/// * The metric is always >= 0
///
/// * If the metric is zero then the two items are the same
///
/// * The metric between A and B is the same as that between B and A
///
/// * A metric between A and B, plus the metric between B and C, must be less than or equal to that between A and C
///
/// One metric (maximum parameterwise distance, call it dM) would be the maximum distance between two Beziers
/// for the same parameter t (i.e. max[t] |B(t)-C(t)|).
///
/// One metric is the square-root of the sum of the distance-squred between the control points of two Beziers
/// (this is termed here dF)
///
/// Another is the maximum of the difference betwen the control points of two Beziers (termed dC). This
/// has a neat property that both Beziers with a metric will lie within (1) the first Bezier plus a
/// ball (of D dimensions) of radius dC, and (2) within the second Bezier plus a similar ball. Another way to put
/// this is that the maximum distance between any two points on the two Beziers for the same parameter.
///
/// Note that dF is guaranteed to be >= dC (since it is only *one* control point), but df <= sqrt(N)*dC as there are
/// N points.
///
/// As more points are added to a Bezier, N increases, and dF can increase; dM, howwver, is independent of N.
///
/// # Least-squares reduction of degree N+1 to degree N
///
/// The reduction matrix for degree n+1 to n is given by
///
/// RL2 = E.transpose() * (E*E.transpose()).inverse()
///
/// Note that E * RL2 = E * E.transpose * (E*E.transpose()).inverse() = X * X.inverse() = identity
///
/// # Basis reduction using points on the Bezier
///
/// The polynomial stored in the structure used the Bernstein polynomials as a basis.
///
/// That is, point(t) = BernsteinMatrix(t) * Polynomials
///
/// If we select 'N' different values of t, ti, with (0<=ti<=1), then we find 'N' points pi
/// where pi = Polynomials * B(ti)
/// i.e. [pi] = Polynomials * [B(ti)]
/// or Polynomials = [pi] * [B(ti]].inverse()
///
/// So if we want to determine Polynomials we can do so from N points pi at N values ti. These are points on the
/// Bezier.
///
/// To reduce from degree N+1 to N we will require N+1 points (as a degree N curve requires N+1 points). We can
/// use *any* N+1 ti in the range 0<=ti<=1, as long as they are all different values.
///
/// We must in some sense find pi for these N+1 points given these ti for order N+1, and then map them back to
/// the Bernstein basis for order N.
///
/// We will get
///
///   Polynomials[N] = [pi] * [B[N](ti)].inverse()
///   Polynomials[N] = Polynomials[N+1] * [B[N+1](ti)] * B[N](ti)].inverse()
///
/// i.e. we can use a reduction matrix
///
///   R[t0..tN] = [B[N+1](ti)] * [B[N](ti)].inverse()
///
/// Note that elevating-by-one of B[N](t) produces B[N+1](t), so
///   E[N] * R[t0..tN] = E[N] * [B[N+1](ti)] * [B[N](ti)].inverse()
///                    = [E[N] * B[N+1](ti)] * [B[N](ti)].inverse()
///                    = [B[N](ti)] * [B[N](ti)].inverse()
///                    = Identity
///
/// i.e. the reduction chosen has the property that elevate(reduce) is the identity, which we need for a reduction
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Bezier<F, const N: usize, const D: usize>
where
    F: Float,
{
    /// Degree, which is one less than number of valid control points
    ///
    /// i.e. for a cubic this is three
    ///
    /// A Bezier of degree 0 is a fixed point
    degree: usize,
    /// Control points - 0..=degree are valid
    pts: [[F; D]; N],
}

//ti Default for Bezier
impl<F, const N: usize, const D: usize> std::default::Default for Bezier<F, N, D>
where
    F: Float,
{
    fn default() -> Self {
        let pts = [[F::zero(); D]; N];
        Self { degree: 0, pts }
    }
}

//ti Display for Bezier
impl<F, const N: usize, const D: usize> std::fmt::Display for Bezier<F, N, D>
where
    F: Float,
{
    //mp fmt - format a `Bezier` for display
    /// Display the `Bezier' as sets of points
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        let mut needs_comma = false;
        for p in self.pts.iter().take(self.degree + 1) {
            if needs_comma {
                write!(f, ", ")?;
            }
            needs_comma = true;
            vector::fmt(f, p);
        }
        Ok(())
    }

    //zz All done
}

//ip Bezier
impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
    //fp new
    /// Create a new Bezier given all the control points
    pub fn new(pts: &[[F; D]]) -> Self {
        let n = pts.len();
        assert!(n >= 1, "Beziers have at least 1 point");
        assert!(
            n < N,
            "Attempt to create a Bezier of max {N} pts with actually {n} pts"
        );
        let mut s = Self::default();
        s.degree = n - 1;
        s.pts.split_at_mut(n).0.copy_from_slice(pts);
        s
    }
    //mp derivative
    /// Create a new Bezier that is the derivative of this
    ///
    /// This is the Bezier of degre N-1 that has points `P'[i]` where
    ///
    ///   `P'[i] = N * (P[i+1] - P[i])`
    pub fn derivative(&self) -> Self {
        assert!(
            self.degree >= 1,
            "Bezier must be degree one or higher to have a derivative"
        );
        let mut s = Self::default();
        s.degree = self.degree - 1;
        let n: F = (self.degree as f32).into();
        for i in 0..self.degree {
            s.pts[i] = vector::sub(self.pts[i + 1], &self.pts[i], n);
        }
        s
    }

    //mp bernstein_basis_coeff
    #[inline]
    pub fn bernstein_basis_coeff(degree: usize, i: usize, t: F) -> F {
        let u = F::one() - t;
        let coeffs = BINOMIALS[degree];
        t.powi(i as i32) * u.powi((degree - i) as i32) * (coeffs[1 + i]).into()
    }

    //mp elevate
    /// Elevate a Bezier by one degree
    ///
    /// This is the Bezier of degree N+1 that has points `P'[i]` where
    ///
    /// ``` ignore
    ///   P'[0] = P[0]
    ///   P'[j] = Sum((n i).(1 j-i) / (n+1 j).P[i])
    ///   P'[j] = Sum( n! / (n-i)!i! * 1!(1-j+i)!/(1-j+i)! * j!(n+1-j)!/(n+1)! . P[i]) (0<=i<=j)
    ///   P'[j] = Sum( n! / (n-i)!i! * j!(n+1-j)!/(n+1)! . P[i]) (i==j-1 or i==j)
    ///   P'[j] = Sum( n!/(n+1)! * j!/i! * (n+1-j)! / (n-i)! . P[i]) (i==j-1 or i==j)
    ///   P'[j] = (j * P[i-1] + (n+1-j) * P[i]) / (n+1)
    ///   P'[N+1] = P[N]
    /// ```
    ///
    /// Elevation matrix element (col i, row j) for degree n to n+1
    /// 0 <= i <= n
    /// 0 <= j <= n+1
    #[inline]
    pub fn elevation_by_one_matrix_ele(n: usize, i: usize, j: usize) -> F {
        if j == 0 && i == 0 {
            F::one()
        } else if j == n + 1 && i == n {
            F::one()
        } else if i + 1 == j {
            // note j!=n+1 as if j == n+1 then i = n, and previous case is used
            // note j != 0
            ((j as f32) / (n + 1) as f32).into()
        } else if i == j {
            // note j!=0 as if j == 0 then i=0 and previous case is used
            // note j != n+1 as i<=n
            (((n + 1 - j) as f32) / ((n + 1) as f32)).into()
        } else {
            F::zero()
        }
    }

    pub fn elevate_by_one(&self) -> Self {
        assert!(
            self.degree < N + 2,
            "Cannot elevate Bezier<{N}> which already has {} pts",
            self.degree + 1
        );
        // n = number of points, i.e. self.pts[n-1] is the last valid point
        let n = self.degree + 1;
        let mut s = Self::default();
        s.degree = n;
        s.pts[0] = self.pts[0];
        s.pts[self.degree] = self.pts[self.degree];
        let scale: F = (1.0 / ((n + 1) as f32)).into();
        for j in 1..self.degree {
            s.pts[j] = vector::scale(
                vector::add(
                    vector::scale(self.pts[j], (j as f32).into()),
                    &self.pts[j + 1],
                    ((n + 1 - j) as f32).into(),
                ),
                scale,
            );
        }
        s
    }

    //mp apply_matrix
    pub fn apply_matrix(&mut self, data: &[F], new_degree: usize) {
        // geo_nd::matrix::inverse4(m)
        todo!();
    }

    //mp degree
    /// Return degree of the Bezier (e.g. 3 for a cubic)
    pub fn degree(&self) -> usize {
        self.degree
    }

    //mp map_pts
    /// Apply a function to all of the points in the Bezier
    pub fn map_pts<Map: Fn([F; D]) -> [F; D]>(&mut self, map: Map) {
        for p in self.pts.iter_mut().take(self.degree + 1) {
            *p = map(*p);
        }
    }

    //mp metric_dm_est
    /// The maximum difference between two beziers given a step dt
    ///
    /// This is an approximation meant for testing/analysis - do not use if
    /// performance is required!
    pub fn metric_dm_est(&self, other: &Self, num_steps: usize) -> F {
        let mut d2 = F::zero();
        for i in 0..num_steps {
            let t: F = (i as f32).into();
            d2 = d2.max(vector::distance_sq(&self.point_at(t), &other.point_at(t)));
        }
        d2.sqrt()
    }

    //mp metric_df
    /// Rhe square-root of the sum of the distance-squred between the control points of two Beziers
    pub fn metric_df(&self, other: &Self) -> F {
        assert!(
            self.degree == other.degree,
            "Degrees of Beziers must match for a metric"
        );
        let mut d2 = F::zero();
        for (s, o) in self.pts.iter().zip(other.pts.iter()).take(self.degree + 1) {
            d2 += vector::distance_sq(s, o);
        }
        d2.sqrt()
    }

    /// Maximum of the difference betwen the control points of two Beziers
    pub fn metric_dc(&self, other: &Self) -> F {
        assert!(
            self.degree == other.degree,
            "Degrees of Beziers must match for a metric"
        );
        let mut d2 = F::zero();
        for (s, o) in self.pts.iter().zip(other.pts.iter()).take(self.degree + 1) {
            d2 = d2.max(vector::distance_sq(s, o));
        }
        d2.sqrt()
    }

    //mp scale
    /// Scale the Bezier by applying the scale factor to all of the points
    ///
    /// This is an example of the [Bezier::map_pts] method
    pub fn scale(&mut self, s: F) {
        self.map_pts(|p| vector::scale(p, s));
    }

    //mi vector_of
    /// Returns a vector of a combination of the vectors of the bezier
    #[inline]
    pub fn vector_of(&self, sc: &[F], reduce: F) -> [F; D] {
        let mut r = [F::zero(); D];
        for (pt, sc) in self.pts.iter().zip(sc.iter()) {
            for (j, rj) in r.iter_mut().enumerate() {
                *rj += *sc * (*pt)[j];
            }
        }
        vector::reduce(r, reduce)
    }

    //mp point_at
    /// Returns the point at parameter 't' along the Bezier
    pub fn point_at(&self, t: F) -> [F; D] {
        let mut r = [F::zero(); D];
        let u = F::one() - t;
        let n = self.degree;
        let coeffs = BINOMIALS[n];
        for (i, (c, pt)) in coeffs[1..].iter().zip(self.pts.iter()).enumerate() {
            let scale = t.powi(i as i32) * u.powi((n - i) as i32) * (*c).into();
            r = vector::add(r, pt, scale);
        }
        r
    }

    //mp tangent_at
    /// Returns the tangent vector at parameter 't' along the Bezier
    ///
    /// Note that this is not necessarily a unit vector
    pub fn tangent_at(&self, t: F) -> [F; D] {
        todo!();
    }

    //mp bisect
    /// Returns two Bezier's that split the curve at parameter t=0.5
    ///
    /// For quadratics the midpoint is 1/4(p0 + 2*c + p1)
    pub fn bisect(&self) -> (Self, Self) {
        todo!();
    }

    //mp bezier_between
    /// Returns the Bezier that is a subset of this Bezier between two parameters 0 <= t0 < t1 <= 1
    pub fn bezier_between(&self, t0: F, t1: F) -> Self {
        todo!();
    }

    //mp as_lines
    /// Return a [BezierLineIter] iterator that provides line segments
    /// when the Bezier is broken down into 'straight' enough through
    /// bisection.
    pub fn as_lines(&self, straightness: F) -> BezierLineIter<F, D> {
        todo!();
    }

    //mp as_points
    /// Return a [BezierPointIter] iterator that provides points along
    /// the curve when the Bezier is broken down into 'straight'
    /// enough through bisection.
    pub fn as_points(&self, straightness: F) -> BezierPointIter<F, D> {
        todo!();
    }

    //mp is_straight
    /// Returns true if the Bezier is straighter than a 'straightness' measure
    ///
    /// A linear bezier is always straight.
    ///
    /// A straightness measure for a quadratic bezier (one control
    /// point) can be thought of as the ratio between the area of the
    /// triangle formed by the two endpoints and the control point
    /// (three points must form a triangle on a plane) in relation to
    /// the distance between the endpoints (the curve will be entirely
    /// within the triangle.
    ///
    /// A straightness measure for a cubic bezier (two control points)
    /// can be though of similarly, except that the curve now must fit
    /// within a volume given by the two control points and the
    /// endpoints; hence the straightness is measured in some way by
    /// the volume in relation to the distance between the endpoints,
    /// but also should be no straighter than the area of any one
    /// control point in relation to the disnance between the
    /// endpoints (the Bezier may be a planar curve that is quite
    /// unstraight but with a volume of zero).
    ///
    /// Hence the straightness here is defined as the sum of (the
    /// ratio between (the distance of each control point from the
    /// straight line between the two endpoints) and (the distance
    /// between the two endpoints))
    ///
    /// `straightness` is thus independent of the length of the Bezier
    pub fn is_straight(&self, straightness: F) -> bool {
        todo!();
    }

    //mp length
    /// Calculates the length of the Bezier when it is rendered down
    /// to the given a straightness
    ///
    /// `straightness` is independent of the length of the Bezier
    pub fn length(&self, straightness: F) -> F {
        todo!();
    }

    //fp arc
    /// Create a Cubic Bezier that approximates closely a circular arc
    ///
    /// The arc has a center C, a radius R, and is of an angle (should be <= PI/2).
    ///
    /// The arc sweeps through points a distance R from C, in a circle
    /// using a pair of the planar unit vectors in the vector space for the
    /// points.
    ///
    /// The arc will be between an angle A1 and A2, where A2-A1 == angle, and A1==rotate
    ///
    pub fn arc(
        angle: F,
        radius: F,
        center: &[F; D],
        unit: &[F; D],
        normal: &[F; D],
        rotate: F,
    ) -> Self {
        todo!();
    }

    //zz All done
}
