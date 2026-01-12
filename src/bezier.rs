//a Imports
use geo_nd::vector;
use geo_nd::Float;

use crate::{BezierLineIter, BezierPointIter};

const BINOMIALS: &[&[f32]] = &[
    &[1., 1.],
    &[1., 2., 1.],
    &[1., 3., 3., 1.],
    &[1., 4., 6., 4., 1.],
    &[1., 5., 10., 10., 5., 1.],
    &[1., 6., 15., 20., 15., 6., 1.],
    &[1., 7., 21., 35., 35., 21., 7., 1.],
    &[1., 8., 28., 56., 70., 56., 28., 8., 1.],
    &[1., 9., 36., 84., 126., 126., 84., 36., 9., 1.],
];

//a Bezier
//tp Bezier
/// A [Bezier] is an implementation of an arbitrary bezier
///
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Bezier<F, const N: usize, const D: usize>
where
    F: Float,
{
    /// Number of valid control points
    num: usize,
    /// Control points - 0..(num-1) are valid
    pts: [[F; D]; N],
}

//ti Default for Bezier
impl<F, const N: usize, const D: usize> std::default::Default for Bezier<F, N, D>
where
    F: Float,
{
    fn default() -> Self {
        let pts = [[F::zero(); D]; N];
        Self { num: 2, pts }
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
        for p in self.pts.iter().take(self.num) {
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
        assert!(n > 1, "Beziers have at least 2 points");
        assert!(
            n < N,
            "Attempt to create a Bezier of max {N} pts with actually {n} pts"
        );
        let mut s = Self::default();
        s.num = n;
        s.pts.split_at_mut(n).0.copy_from_slice(pts);
        s
    }

    //mp elevate
    /// Elevate a Bezier by one degree (cannot elevate a Cubic)
    pub fn elevate(&mut self) {
        assert!(
            self.num < N,
            "Cannot elevate Bezier<{N}> which already has {} pts",
            self.num
        );
    }

    //mp degree
    /// Returns number of points used for the Bezier (2 to 4)
    ///
    /// Cubic beziers return 3
    /// Quadratic beziers return 2
    /// Linear beziers (lines...) return 1
    pub fn degree(&self) -> usize {
        self.num - 1
    }

    //mp map_pts
    /// Apply a function to all of the points in the Bezier
    pub fn map_pts<Map: Fn([F; D]) -> [F; D]>(&mut self, map: Map) {
        for p in self.pts.iter_mut().take(self.num) {
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

    //mi vector_of
    /// Returns a vector of a combination of the vectors of the bezier
    #[inline]
    fn vector_of(&self, sc: &[F], reduce: F) -> [F; D] {
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
        let n = self.num;
        for (i, pt) in self.pts.iter().take(n).enumerate() {
            r = vector::add(r, pt, t.powi(i as i32) * u.powi((n - i) as i32));
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
