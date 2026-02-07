//!  A [Bezier] is an implementation of an arbitrary bezier up to a maximum degree
//! given in the type.
//!
//! The control points are stored to be used with Bernstein basis polynomials
//!
//! ``` ignore
//!    b[i,n](t) = n!/(i!(n-i)!) . u^i . t^(n-i) where u=1-t, n=Bezier degree, 0<=i <=n
//! ```
//! A point at parameter t (0 <= t <= 1) is then `Sum(b[i,n](t).pts[i])`
//!
//! If n were two then this is `u^2.pts[0] + 2.u.t.pts[1] + t^2.pts[2]`, which can be rewritten as:
//! ``` ignore
//!   (1-2t+t^2).pts[0] + 2.(1-t).t.pts[1] + t^2.pts[2]
//!   = t^2.(pts[0] - 2.pts[1] + pts[2]) + 2t.(pts[1] - pts[0]) + pts[0]
//! ```
//!
//! i.e. a set of different point values where a monomial basis (1, t, t^2, ...) is used
//!
//! To convert from Bernstein basis to monomial basis for degree 4 (five points) the pts would map
//! using a 5x5 matrix:
//!
//! ``` ignore
//!     |  1   0   0   0   0 |
//!     | -4   4   0   0   0 |
//!     |  6 -12   6   0   0 |
//!     | -4  12 -12   4   0 |
//!     |  1  -4   6  -4   1 |
//!```
//!
//! Actual value Mij for degree N is (-1)^(j-i) . (n C j) . (j c i) for i<=j, 0 otherwise (by inspection...)
//!
//! Note that this matrix is guaranteed to be invertible, hence all Bernstein representations are equivalent
//! to monomial representations, and vice versa
//!
//! Note that `P[N]` is only used in conjunction with monomial t^N, and hence if this were to be used
//! to derive points etc then it would be unstable near t=1. However, it has the potentially useful
//! property (for t<1) that 1 > t > t^2 > t^3 > 0.
//!
//! Another possible basis would be for the basis 1, (t-T), (t-T)^2, etc for a fixed T (effectively
//! monomial is this with T=0). This is termed the Taylor basis vectors in <https://arxiv.org/pdf/2201.07834>
//! This has the property that it is locally (at t=0) approximatable (as one can ignore high powers of (t-T)).
//!
//! Metrics ('distance' between two Beziers)
//!
//! Various metrics (real number) can be used for 'measuring' the difference between two Beziers. Note that a metric
//! is a special kind of measurement (such as the distance between two points), where:
//!
//! * The metric is always >= 0
//!
//! * If the metric is zero then the two items are the same
//!
//! * The metric between A and B is the same as that between B and A
//!
//! * A metric between A and B, plus the metric between B and C, must be less than or equal to that between A and C
//!
//! One metric (maximum parameterwise distance, call it dM) would be the maximum distance between two Beziers
//! for the same parameter t (i.e. `max[t] |B(t)-C(t)|`).
//!
//! One metric is the square-root of the sum of the distance-squred between the control points of two Beziers
//! (this is termed here dF)
//!
//! Another is the maximum of the difference betwen the control points of two Beziers (termed dC). This
//! has a neat property that both Beziers with a metric will lie within (1) the first Bezier plus a
//! ball (of D dimensions) of radius dC, and (2) within the second Bezier plus a similar ball. Another way to put
//! this is that the maximum distance between any two points on the two Beziers for the same parameter.
//!
//! Note that dF is guaranteed to be >= dC (since it is only *one* control point), but df <= sqrt(N)*dC as there are
//! N points.
//!
//! As more points are added to a Bezier, N increases, and dF can increase; dM, howwver, is independent of N.
//!
//! # Least-squares reduction of degree N+1 to degree N
//!
//! The reduction matrix for degree n+1 to n is given by
//!
//! RL2 = E.transpose() * (E*E.transpose()).inverse()
//!
//! Note that E * RL2 = E * E.transpose * (E*E.transpose()).inverse() = X * X.inverse() = identity
//!
//! # Basis reduction using points on the Bezier
//!
//! The polynomial stored in the structure used the Bernstein polynomials as a basis.
//!
//! That is, point(t) = BernsteinMatrix(t) * Polynomials
//!
//! If we select 'N' different values of t, ti, with (0<=ti<=1), then we find 'N' points pi
//! where `pi = Polynomials * B(ti)`
//! i.e. `[pi] = Polynomials * [B(ti)]`
//! or `Polynomials = [pi] * [B(ti]].inverse()`
//!
//! So if we want to determine Polynomials we can do so from N points pi at N values ti. These are points on the
//! Bezier.
//!
//! To reduce from degree N+1 to N we will require N+1 points (as a degree N curve requires N+1 points). We can
//! use *any* N+1 ti in the range 0<=ti<=1, as long as they are all different values.
//!
//! We must in some sense find pi for these N+1 points given these ti for order N+1, and then map them back to
//! the Bernstein basis for order N.
//!
//! We will get
//!
//! ```text
//!   Polynomials[N] = [pi] * [B[N](ti)].inverse()
//!   Polynomials[N] = Polynomials[N+1] * [B[N+1](ti)] * B[N](ti)].inverse()
//! ```
//!
//! i.e. we can use a reduction matrix
//!
//! ```text
//!   R[t0..tN] = [B[N+1](ti)] * [B[N](ti)].inverse()
//! ```
//!
//! Note that elevating-by-one of `B[N](t)` produces `B[N+1](t)`, so
//!
//! ```text
//!   E[N] * R[t0..tN] = E[N] * [B[N+1](ti)] * [B[N](ti)].inverse()
//!                    = [E[N] * B[N+1](ti)] * [B[N](ti)].inverse()
//!                    = [B[N](ti)] * [B[N](ti)].inverse()
//!                    = Identity
//! ```
//!
//! i.e. the reduction chosen has the property that elevate(reduce) is the identity, which we need for a reduction

use crate::BezierEval;
use crate::{bernstein_fns, BezierBuilder, BezierConstruct};
use crate::{BezierOps, BezierSection, BezierSplit, Num};

use geo_nd::matrix;
use geo_nd::vector;

/// This type stores a Bernstein Bezier of up to N control points each of dimension D
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct BezierND<F, const N: usize, const D: usize>
where
    F: Num,
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
impl<F, const N: usize, const D: usize> std::default::Default for BezierND<F, N, D>
where
    F: Num,
{
    fn default() -> Self {
        let pts = [[F::ZERO; D]; N];
        Self { degree: 0, pts }
    }
}

//ti Display for Bezier
impl<F, const N: usize, const D: usize> std::fmt::Display for BezierND<F, N, D>
where
    F: Num,
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

//ip Bezier iterators
impl<F, const N: usize, const D: usize> BezierND<F, N, D>
where
    F: Num,
{
    /// Find the maximum degree this type can encode
    pub fn max_degree() -> usize {
        N - 1
    }

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
}

//tp BezierReduceIter
/// A type that provided an iterator implementaion for splitting a Bezier
/// into reduced Beziers given a maximum 'dc' metric
#[derive(Clone, PartialEq, Debug)]
pub struct BezierReduceIter<'a, F, const N: usize, const D: usize>
where
    F: Num,
{
    reduce_matrix: &'a [F],
    elev_reduce_matrix: &'a [F],
    reduce_degree: usize,
    max_dc_sq: F,
    stack: Vec<(usize, BezierND<F, N, D>)>,
}

impl<F, const N: usize, const D: usize> BezierND<F, N, D>
where
    F: Num,
{
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

    /// Create a new Bezier that is the nth derivative of this
    /// subject to a scaling factor (i.e. the actual Bezier should be scaled up by F)
    ///
    /// As this type uses Bernstein polynomials, this is the Bezier of degree N-n
    /// that has points control points:
    ///
    /// `Pn[i] = (N n) . Sum((-1)^j * P[i+n-j] * (n j))`
    ///
    /// (n j) is kept in `bezier_nd::constants::BINOMIALS[n][j+1]`
    pub fn nth_derivative(&self, n: usize) -> (Self, F) {
        let mut s = Self {
            degree: self.degree - n,
            ..Default::default()
        };
        s.degree = self.degree - n;
        let scale = bernstein_fns::nth_derivative(&self.pts[0..self.degree + 1], n, &mut s.pts);
        (s, scale)
    }

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
}

impl<F, const N: usize, const D: usize> BezierND<F, N, D>
where
    F: Num,
{
    /// Returns the point at parameter 't' along the Bezier using de Casteljau's algorithm
    ///
    /// This does not require powi, but it is an O(n^2) operations and O(n) space operation
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

impl<F, const N: usize, const D: usize> BezierEval<F, [F; D]> for BezierND<F, N, D>
where
    F: Num,
{
    fn point_at(&self, t: F) -> [F; D] {
        let mut r = [F::zero(); D];
        for ((_i, c), pt) in
            bernstein_fns::basis_coeff_enum_num(self.degree, t).zip(self.pts.iter())
        {
            r = vector::add(r, pt, c);
        }
        r
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        (F::ONE, self.nth_derivative_value_at(1, t))
    }

    /// Borrow the endpoints of the Bezier
    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self.pts[0], &self.pts[self.degree])
    }

    fn closeness_sq_to_line(&self) -> F {
        self.dc_sq_from_line()
    }
    fn dc_sq_from_line(&self) -> F {
        if self.degree < 2 {
            F::ZERO
        } else {
            let n = (self.degree + 1) as f32;
            let mut max_dc_sq = F::ZERO;
            let p01 = [self.pts[0], self.pts[self.degree]];
            for (i, p) in self.pts.iter().take(self.degree).skip(1).enumerate() {
                let t = ((i as f32) / n).into();
                let pt = vector::sum_scaled(&p01, &[F::ONE - t, t]);
                let dc_sq = vector::distance_sq(&pt, p);
                if dc_sq > max_dc_sq {
                    max_dc_sq = dc_sq;
                }
            }
            max_dc_sq
        }
    }

    fn num_control_points(&self) -> usize {
        self.degree + 1
    }

    fn control_point(&self, n: usize) -> &[F; D] {
        &self.pts[n]
    }

    fn degree(&self) -> usize {
        self.degree
    }
    fn for_each_control_point(&self, map: &mut dyn FnMut(usize, &[F; D])) {
        self.pts
            .iter()
            .take(self.degree + 1)
            .enumerate()
            .for_each(|(i, pt)| map(i, pt))
    }
}

impl<F: Num, const N: usize, const D: usize> BezierConstruct<F, D> for BezierND<F, N, D> {
    fn of_builder(builder: &BezierBuilder<F, D>) -> Result<Self, ()> {
        let (mut matrix, pts) = builder.get_matrix_pts()?;
        if pts.len() > N {
            return Err(());
        }
        let degree = pts.len() - 1;
        let bezier = Self::new(&pts);

        let mut lu = matrix.clone();
        let mut pivot = vec![0; degree + 1];
        let mut tr0 = vec![F::zero(); degree + 1];
        let mut tr1 = vec![F::zero(); degree + 1];

        if matrix::lup_decompose(degree + 1, &matrix, &mut lu, &mut pivot) == F::ZERO {
            return Err(());
        }

        assert!(
            matrix::lup_invert(degree + 1, &lu, &pivot, &mut matrix, &mut tr0, &mut tr1),
            "Matrix must be invertible"
        );
        Ok(bezier.apply_matrix(&matrix, degree))
    }
}

//ip Iterator for BezierReduceIter
impl<'a, F, const N: usize, const D: usize> std::iter::Iterator for BezierReduceIter<'a, F, N, D>
where
    F: Num,
{
    type Item = (usize, BezierND<F, N, D>);
    fn next(&mut self) -> Option<(usize, BezierND<F, N, D>)> {
        while let Some((n, b)) = self.stack.pop() {
            let dc2 = b.dc2_of_ele_red(self.elev_reduce_matrix);
            if dc2 > self.max_dc_sq {
                let (b0, b1) = b.split();
                self.stack.push((n + 1, b1));
                self.stack.push((n + 1, b0));
            } else {
                return Some((n, b.apply_matrix(self.reduce_matrix, self.reduce_degree)));
            }
        }
        None
    }
}

impl<F, const N: usize, const D: usize> BezierND<F, N, D>
where
    F: Num,
{
    /// Return a slice of the control points
    pub fn pts(&self) -> &[[F; D]] {
        &self.pts[0..self.degree + 1]
    }
}

impl<F, const N: usize, const D: usize> crate::BezierReduce<F, [F; D]> for BezierND<F, N, D>
where
    F: crate::Num,
{
    type Reduced = Self;
    type Quadratic = Self;
    type Cubic = Self;

    fn reduce(&self) -> Self::Reduced {
        todo!()
    }
    fn can_reduce(&self) -> bool {
        self.degree > 1
    }
    fn closeness_sq_to_reduction(&self) -> Option<F> {
        None
    }

    fn closeness_sq_to_quadratic(&self) -> F {
        todo!()
    }

    fn closeness_sq_to_cubic(&self) -> F {
        todo!()
    }
    fn reduced_to_quadratic(&self) -> Option<Self::Quadratic> {
        None
    }
    fn reduced_to_cubic(&self) -> Option<Self::Cubic> {
        None
    }
}

impl<F, const N: usize, const D: usize> BezierSplit for BezierND<F, N, D>
where
    F: Num,
{
    fn split(&self) -> (Self, Self) {
        self.split_at(0.5_f32.into())
    }
}

impl<F, const N: usize, const D: usize> BezierSection<F> for BezierND<F, N, D>
where
    F: Num,
{
    fn split_at(&self, t: F) -> (Self, Self) {
        let mut first = *self;
        let mut latter = *self;
        bernstein_fns::de_casteljau::split_at(
            &mut latter.pts[0..self.degree + 1],
            t,
            &mut first.pts,
        );
        (first, latter)
    }
    fn section(&self, t0: F, t1: F) -> Self {
        let mut to_split = *self;
        if t0 > F::ZERO {
            bernstein_fns::de_casteljau::bezier_from(&mut to_split.pts, t0);
        }
        if t1 < F::ONE {
            let t10 = (t1 - t0) / (F::ONE - t0);
            bernstein_fns::de_casteljau::bezier_to(&mut to_split.pts, t10);
        }
        to_split
    }
}

impl<F, const N: usize, const D: usize> BezierND<F, N, D>
where
    F: Num,
{
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
        let mut max_d2 = F::zero();
        let nc = self.degree + 1;
        for (m, s) in matrix.chunks_exact(nc).zip(self.pts.iter()) {
            let mut sum = [F::zero(); D];
            for (coeff, p) in m.iter().zip(self.pts.iter()) {
                sum = vector::add(sum, p, *coeff);
            }
            let d2 = vector::distance_sq(s, &sum);
            if max_d2 < d2 {
                max_d2 = d2;
            }
        }
        max_d2
    }
}

impl<F, const N: usize, const D: usize> BezierOps<F, [F; D]> for BezierND<F, N, D>
where
    F: Num,
{
    fn add(&mut self, other: &Self) -> bool {
        if self.degree != other.degree {
            false
        } else {
            for (s, o) in self.pts.iter_mut().zip(other.pts.iter()) {
                *s = vector::add(*s, o, F::ONE)
            }
            true
        }
    }

    /// Subtract another Bezier from this
    fn sub(&mut self, other: &Self) -> bool {
        if self.degree != other.degree {
            false
        } else {
            for (s, o) in self.pts.iter_mut().zip(other.pts.iter()) {
                *s = vector::sub(*s, o, F::ONE)
            }
            true
        }
    }

    /// Scale the points
    fn scale(&mut self, scale: F) {
        self.map_pts(&|_, p| vector::scale(*p, scale));
    }

    /// Map the points
    fn map_pts(&mut self, map: &dyn Fn(usize, &[F; D]) -> [F; D]) {
        for (i, p) in self.pts.iter_mut().take(self.degree + 1).enumerate() {
            *p = map(i, p);
        }
    }
}
