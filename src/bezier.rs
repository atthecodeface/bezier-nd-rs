//a Imports
use geo_nd::vector;
use geo_nd::Float;

use crate::constants::BINOMIALS;

//a Bezier
//tp Bezier
/// A [Bezier] is an implementation of an arbitrary bezier up to a maximum degree
/// given in the type.
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
            vector::fmt(f, p)?;
        }
        write!(f, "]")
    }

    //zz All done
}

//ip Bezier accessors
impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
    //ap degree
    /// Return degree of the Bezier (e.g. 3 for a cubic)
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Return a slice of the control points
    pub fn pts(&self) -> &[[F; D]] {
        &self.pts[0..self.degree + 1]
    }
}

//ip Bezier constructors / splitters etc
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
            n <= N,
            "Attempt to create a Bezier of max {N} pts with actually {n} pts"
        );
        let mut s = Self {
            degree: n - 1,
            ..Default::default()
        };
        s.pts.split_at_mut(n).0.copy_from_slice(pts);
        s
    }

    //mp nth_derivative
    /// Create a new Bezier that is the nth derivative of this
    /// subject to a scaling factor (i.e. the actual Bezier should be scaled up by F)
    ///
    /// As this type uses Bernstein polynomials, this is the Bezier of degree N-n
    /// that has points control points:
    ///
    /// Pn[i] = (N n) . Sum((-1)^j * P[i+n-j] * (n j))
    ///
    /// (n j) is kept in BINOMIALS[n][j+1]
    pub fn nth_derivative(&self, n: usize) -> (Self, F) {
        assert!(
            self.degree >= n,
            "Bezier of degree {N} must actually be of degree {n} or higher to have an nth derivative"
        );
        let mut s = Self {
            degree: self.degree - n,
            ..Default::default()
        };
        s.degree = self.degree - n;
        let mut scale = F::one();
        for i in 0..n {
            scale *= ((self.degree - i) as f32).into();
        }
        for (i, sp) in s.pts.iter_mut().take(self.degree + 1 - n).enumerate() {
            let mut m1_n_positive = (n & 1) == 0;
            for (c, p) in BINOMIALS[n][1..].iter().zip(self.pts[i..].iter()) {
                if m1_n_positive {
                    *sp = vector::add(*sp, p, (*c).into());
                } else {
                    *sp = vector::sub(*sp, p, (*c).into());
                }
                m1_n_positive = !m1_n_positive;
            }
        }
        (s, scale)
    }

    //mp split_at_de_cast
    /// Use de Casteljau's algorithm to split
    pub fn split_at_de_cast(&self, t: F) -> (Self, Self) {
        let mut s0 = *self;
        let mut s1 = *self;
        // Beta[0][i] = pt[i]
        let mut pts = self.pts;
        // for j = 1..=n
        //  for i = 0..=n-j
        // Beta[j][i] = (1-t)*Beta[j-1][i] + t*Beta[j-1][i+1]
        let u = F::one() - t;
        for j in 1..(self.degree + 1) {
            for i in 0..(self.degree + 1 - j) {
                pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
            }
            s0.pts[j] = pts[0];
            s1.pts[self.degree - j] = pts[self.degree - j];
        }
        (s0, s1)
    }

    //mp apply_matrix
    /// Apply a (new_degree+1) by (degree+1) matrix to the points to generate a new Bezier
    /// of a new degree
    pub fn apply_matrix(&self, matrix: &[F], new_degree: usize) -> Self {
        assert!(
            new_degree < N + 1,
            "Cannot create a Bezier<{N}> of {new_degree}",
        );
        assert_eq!(
            matrix.len(),
            (new_degree + 1) * (self.degree + 1),
            "Matrix to apply to Bezier of degree {} to degree {new_degree} must have {} elements",
            self.degree,
            (new_degree + 1) * (self.degree + 1)
        );
        let _nr = new_degree + 1;
        let nc = self.degree + 1;
        let mut s = Self {
            degree: new_degree,
            ..Default::default()
        };
        for (m, sp) in matrix.chunks_exact(nc).zip(s.pts.iter_mut()) {
            let mut sum = [F::zero(); D];
            for (coeff, p) in m.iter().zip(self.pts.iter()) {
                sum = vector::add(sum, p, *coeff);
            }
            *sp = sum;
        }
        s
    }

    //mp bisect
    /// Returns two Bezier's that split the curve at parameter t=0.5
    ///
    /// For quadratics the midpoint is 1/4(p0 + 2*c + p1)
    pub fn bisect(&self) -> (Self, Self) {
        self.split_at_de_cast(0.5_f32.into())
    }

    //mp elevate
    /// Elevate a Bezier by one degree
    ///
    /// This generates and applies the elevate-by-one matrix
    pub fn elevate_by_one(&self) -> Self {
        assert!(
            self.degree < N + 2,
            "Cannot elevate Bezier<{N}> which already has {} pts",
            self.degree + 1
        );
        // n = number of points, i.e. self.pts[n-1] is the last valid point
        let n = self.degree + 1;
        let mut s = Self {
            degree: n,
            ..Default::default()
        };
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

    //mp bezier_between
    /// Returns the Bezier that is a subset of this Bezier between two parameters 0 <= t0 < t1 <= 1
    pub fn bezier_between(&self, t0: F, t1: F) -> Self {
        let dt = t1 - t0;
        assert!(t0 < F::one(), "Must select a t0 that is less than 1.0");
        assert!(dt > F::zero(), "Must select a t range that is > 0");
        if t0 == F::zero() {
            self.split_at_de_cast(dt).0
        } else if t1 == F::one() {
            self.split_at_de_cast(dt).1
        } else {
            // b is the subset from t0 to 1.0
            let (_, b) = self.split_at_de_cast(t0);
            b.split_at_de_cast(dt / (F::one() - t0)).0
        }
    }
}

//ip Bezier metrics
impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
    //mp metric_dm_est
    /// The maximum difference between two beziers given a step dt
    ///
    /// This is an approximation meant for testing/analysis - do not use if
    /// performance is required!
    pub fn metric_dm_est(&self, other: &Self, num_steps: usize) -> F {
        let mut d2 = F::zero();
        let ns: F = (num_steps as f32).into();
        for i in 0..num_steps {
            let t: F = (i as f32).into();
            let t = t / ns;
            d2 = d2.max(vector::distance_sq(&self.point_at(t), &other.point_at(t)));
        }
        d2.sqrt()
    }

    //mp metric_df
    /// The square-root of the sum of the distance-squred between the control points of two Beziers
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

    //mp metric_dc
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

    //mp dc_of_ele_red
    /// Apply a (degree+1) by (degree+1) matrix (should be elevate-of-reduce) to the points
    /// and calculate the new dc squared
    pub fn dc2_of_ele_red(&self, matrix: &[F]) -> F {
        assert_eq!(
            matrix.len(),
            (self.degree + 1) * (self.degree + 1),
            "Matrix to apply to Bezier of degree {} must have {} elements",
            self.degree,
            (self.degree + 1) * (self.degree + 1)
        );
        let mut d2 = F::zero();
        let nc = self.degree + 1;
        for (m, s) in matrix.chunks_exact(nc).zip(self.pts.iter()) {
            let mut sum = [F::zero(); D];
            for (coeff, p) in m.iter().zip(self.pts.iter()) {
                sum = vector::add(sum, p, *coeff);
            }
            d2 = d2.max(vector::distance_sq(s, &sum));
        }
        d2
    }
}
//ip Bezier manipulation
impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
    //mp map_pts
    /// Apply a function to all of the points in the Bezier
    pub fn map_pts<Map: Fn([F; D]) -> [F; D]>(&mut self, map: Map) {
        for p in self.pts.iter_mut().take(self.degree + 1) {
            *p = map(*p);
        }
    }

    //mp scale
    /// Scale the Bezier by applying the scale factor to all of the points
    ///
    /// This is an example of the [Bezier::map_pts] method
    pub fn scale(&mut self, s: F) {
        self.map_pts(|p| vector::scale(p, s));
    }
}

//ip Bezier evaluation
impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
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

    //mp point_at_de_cast
    /// Returns the point at parameter 't' along the Bezier using de Casteljau's algorithm
    pub fn point_at_de_cast(&self, t: F) -> [F; D] {
        // Beta[0][i] = pt[i]
        let mut pts = self.pts;
        // for j = 1..=n
        //  for i = 0..=n-j
        // Beta[j][i] = (1-t)*Beta[j-1][i] + t*Beta[j-1][i+1]
        let u = F::one() - t;
        for j in 1..(self.degree + 1) {
            for i in 0..(self.degree + 1 - j) {
                pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
            }
        }
        pts[0]
    }

    //mp point_at
    /// Returns the point at parameter 't' along the Bezier
    pub fn point_at(&self, t: F) -> [F; D] {
        let mut r = [F::zero(); D];
        let u = F::one() - t;
        let n = self.degree;
        for (i, (c, pt)) in BINOMIALS[n][1..].iter().zip(self.pts.iter()).enumerate() {
            let scale = t.powi(i as i32) * u.powi((n - i) as i32) * (*c).into();
            r = vector::add(r, pt, scale);
        }
        r
    }

    //mp nth_derivative_value_at
    /// Returns the value of the nth deriviative at parameter 't' along the Bezier
    ///
    /// This could be optimized to not store the points, but that seems ultimately to be
    /// unnecessary
    pub fn nth_derivative_value_at(&self, n: usize, t: F) -> [F; D] {
        if n > self.degree {
            [F::zero(); D]
        } else {
            let (dn, f) = self.nth_derivative(n);
            vector::scale(dn.point_at_de_cast(t), f)
        }
    }
}

//ip Bezier iterators
impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Float,
{
    //mp Reduce-and-split iterator
    /// Apply a (new_degree+1) by (degree+1) matrix to the points to generate a new Bezier
    /// of a new degree
    pub fn reduce_and_split_iter<'a>(
        &'a self,
        reduce_matrix: &'a [F],
        elev_reduce_matrix: &'a [F],
        reduce_degree: usize,
        max_dc_sq: F,
    ) -> BezierReduceIter<'a, F, N, D> {
        BezierReduceIter {
            reduce_matrix,
            elev_reduce_matrix,
            reduce_degree,
            max_dc_sq,
            stack: vec![],
        }
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
        let s2 = straightness * straightness;
        self.pts
            .iter()
            .skip(1)
            .take(self.degree)
            .all(|m| vector::length_sq(m) <= s2)
    }

    //zz All done
}

//tp BezierReduceIter
/// A type that provided an iterator implementaion for splitting a Bezier
/// into reduced Beziers given a maximum 'dc' metric
#[derive(Clone, PartialEq, Debug)]
pub struct BezierReduceIter<'a, F, const N: usize, const D: usize>
where
    F: Float,
{
    reduce_matrix: &'a [F],
    elev_reduce_matrix: &'a [F],
    reduce_degree: usize,
    max_dc_sq: F,
    stack: Vec<(usize, Bezier<F, N, D>)>,
}

//ip Iterator for BezierReduceIter
impl<'a, F, const N: usize, const D: usize> std::iter::Iterator for BezierReduceIter<'a, F, N, D>
where
    F: Float,
{
    type Item = (usize, Bezier<F, N, D>);
    fn next(&mut self) -> Option<(usize, Bezier<F, N, D>)> {
        while let Some((n, b)) = self.stack.pop() {
            let dc2 = b.dc2_of_ele_red(self.elev_reduce_matrix);
            if dc2 > self.max_dc_sq {
                let (b0, b1) = b.bisect();
                self.stack.push((n + 1, b1));
                self.stack.push((n + 1, b0));
            } else {
                return Some((n, b.apply_matrix(self.reduce_matrix, self.reduce_degree)));
            }
        }
        None
    }
}
