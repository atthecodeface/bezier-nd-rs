//a Imports
use crate::{bernstein_fns, Float};
use geo_nd::vector;

use crate::{
    BezierElevate, BezierEval, BezierIntoIterator, BezierMinMax, BezierOps, BezierReduce,
    BezierSection, BezierSplit,
};

//a Bezier
//tp Bezier
/// A [Bezier] is an implementation of a linear, quadratic or cubic Bezier
/// curve using a parameter which has the [Float] trait, and consists
/// of points of `[F; D]`.
///
/// To split a quadratic bezier at t is simple: the split point is p(t),
/// and the two control points (cl, cr) are:
///
///   cl(t) = u.p0 + t.c ; cr = u.c + t.p1
///
/// Hence the Quadratic Bezier between t0 and t1 can be calculated
/// by splitting to get the right-hand Bezier of t0->1, and splitting
/// this to get the left-hand Bezier at (t1-t0)/u0 = (t2,u2)
///
///    Note `t2 = (t1-t0)/u0; u2=1-t2 = (u0+t0-t1)/u0 = (1-t1)/u0 = u1/u0`
///
/// ```text
///    cl(t0) = u0.p0 + t0.c
///    cr(t0) = u0.c  + t1.p1
///     p(t0) = u0.cl(t0)  + t0.cr(t0)
/// ```
///
///    Bezier t0->1 : p(t0), cr(t0), p1
///
/// ```text
///  c(t0,t1)  = u2.p(t0)  + t2.cr(t0)
///            = u2.u0.cl(t0) + u2.t0.cr(t0) + t2.cr(t0)
///            = u2.u0.cl(t0) + (u2.t0+t2).cr(t0)
///  But u2.u0    = u1
///  And u2.t0+t2 = u1/u0.t0+(t1-t0)/u0
///               = (t0.u1+t1-t0)/u0
///               = (t0 - t1.t0 + t1 - t0) / u0
///               = (t1 - t1.t0) / u0
///               = t1(1-t0) / (1-t0)
///               = t1
/// ```
///  Hence
/// ```text
///  c(t0,t1)  = u1.cl(t0) + t1.cr(t0)
///            = u0.u1.p0 + u1.t0.c + u0.t1.c + t0.t1.p1
///            = u0.u1.p0 + (u1.t0+u0.t1).c + t0.t1.p1
/// ```
///  And the points are:
/// ```text
///      p(t0) = u0.u0.p0 + 2(u0.t0).c + t0.t0.p1
///      p(t1) = u1.u1.p0 + 2(u1.t1).c + t1.t1.p1
/// ```
///
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Bezier<F, const D: usize>
where
    F: Float,
{
    /// Number of valid control points (2-4)
    num: usize,
    /// Control points - endpoints are always 0 and [num-1]
    pts: [[F; D]; 4],
}

//ti Display for Bezier
impl<F, const D: usize> std::fmt::Display for Bezier<F, D>
where
    F: Float,
{
    //mp fmt - format a `Bezier` for display
    /// Display the `Bezier' as sets of points
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        vector::fmt(f, &self.pts[0])?;
        write!(f, "<-")?;
        if self.num > 2 {
            vector::fmt(f, &self.pts[1])?;
        }
        if self.num > 3 {
            write!(f, ":")?;
            vector::fmt(f, &self.pts[2])?;
        }
        write!(f, "->")?;
        vector::fmt(f, &self.pts[self.num - 1])
    }

    //zz All done
}

impl<F: Float, const D: usize> std::convert::From<[[F; D]; 2]> for Bezier<F, D> {
    fn from(value: [[F; D]; 2]) -> Self {
        Self {
            num: 2,
            pts: [value[0], value[1], [F::zero(); D], [F::zero(); D]],
        }
    }
}

impl<F: Float, const D: usize> std::convert::From<[[F; D]; 3]> for Bezier<F, D> {
    fn from(value: [[F; D]; 3]) -> Self {
        Self {
            num: 3,
            pts: [value[0], value[1], value[2], [F::zero(); D]],
        }
    }
}

impl<F: Float, const D: usize> std::convert::From<[[F; D]; 4]> for Bezier<F, D> {
    fn from(value: [[F; D]; 4]) -> Self {
        Self { num: 4, pts: value }
    }
}

impl<F, const D: usize> std::convert::TryFrom<&[[F; D]]> for Bezier<F, D>
where
    F: Float,
{
    type Error = ();
    fn try_from(pts: &[[F; D]]) -> Result<Self, ()> {
        match pts.len() {
            2 => Ok([pts[0], pts[1]].into()),
            3 => Ok([pts[0], pts[1], pts[2]].into()),
            4 => Ok([pts[0], pts[1], pts[2], pts[3]].into()),
            _ => Err(()),
        }
    }
}

impl<F, const D: usize> std::convert::TryFrom<&Bezier<F, D>> for [[F; D]; 2]
where
    F: Float,
{
    type Error = ();
    fn try_from(b: &Bezier<F, D>) -> std::result::Result<[[F; D]; 2], ()> {
        if b.num == 2 {
            Ok([b.pts[0], b.pts[1]])
        } else {
            Err(())
        }
    }
}

impl<F, const D: usize> std::convert::TryFrom<&Bezier<F, D>> for [[F; D]; 3]
where
    F: Float,
{
    type Error = ();
    fn try_from(b: &Bezier<F, D>) -> std::result::Result<[[F; D]; 3], ()> {
        if b.num == 3 {
            Ok([b.pts[0], b.pts[1], b.pts[2]])
        } else {
            Err(())
        }
    }
}

impl<F, const D: usize> std::convert::TryFrom<&Bezier<F, D>> for [[F; D]; 4]
where
    F: Float,
{
    type Error = ();
    fn try_from(b: &Bezier<F, D>) -> std::result::Result<[[F; D]; 4], ()> {
        if b.num == 4 {
            Ok([b.pts[0], b.pts[1], b.pts[2], b.pts[3]])
        } else {
            Err(())
        }
    }
}

impl<F, const D: usize> BezierEval<F, [F; D]> for Bezier<F, D>
where
    F: Float,
{
    fn point_at(&self, t: F) -> [F; D] {
        match self.num {
            2 => self.as_array_2().point_at(t),
            3 => self.as_array_3().point_at(t),
            _ => self.as_array_4().point_at(t),
        }
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        match self.num {
            2 => self.as_array_2().derivative_at(t),
            3 => self.as_array_3().derivative_at(t),
            _ => self.as_array_4().derivative_at(t),
        }
    }

    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        match self.num {
            2 => (&self.pts[0], &self.pts[1]),
            3 => (&self.pts[0], &self.pts[2]),
            _ => (&self.pts[0], &self.pts[3]),
        }
    }

    fn closeness_sq_to_line(&self) -> F {
        match self.num {
            2 => F::ZERO,
            3 => self.as_array_3().closeness_sq_to_line(),
            _ => self.as_array_4().closeness_sq_to_line(),
        }
    }
    fn dc_sq_from_line(&self) -> F {
        match self.num {
            2 => F::ZERO,
            3 => self.as_array_3().dc_sq_from_line(),
            _ => self.as_array_4().dc_sq_from_line(),
        }
    }
    fn num_control_points(&self) -> usize {
        self.num
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self.pts[n]
    }
    fn for_each_control_point(&self, map: &mut dyn FnMut(usize, &[F; D])) {
        self.pts
            .iter()
            .take(self.num)
            .enumerate()
            .for_each(|(i, pt)| map(i, pt))
    }
}

impl<F, const D: usize> BezierElevate<F, [F; D]> for Bezier<F, D>
where
    F: Float,
{
    type ElevatedByOne = Self;
    type Elevated = Self;
    fn elevate_by_one(&self) -> Option<Self> {
        match self.num {
            2 => {
                let c = vector::sum_scaled(&self.pts[0..2], &[(0.5_f32).into(), (0.5_f32).into()]);
                Some(Self::quadratic(&self.pts[0], &c, &self.pts[1]))
            }
            3 => {
                let c0 = vector::sum_scaled(
                    &self.pts[0..2],
                    &[(0.33333333_f32).into(), (0.66666667_f32).into()],
                );
                let c1 = vector::sum_scaled(
                    &self.pts[1..3],
                    &[(0.66666667_f32).into(), (0.33333333_f32).into()],
                );
                Some(Self::cubic(&self.pts[0], &c0, &c1, &self.pts[2]))
            }
            _ => None,
        }
    }
    fn elevate_by(&self, degree: usize) -> Option<Self> {
        match degree {
            0 => Some(self.clone()),
            1 => self.elevate_by_one(),
            2 => self.elevate_by_one().map(|s| s.elevate_by_one()).flatten(),
            _ => None,
        }
    }
}

impl<F, const D: usize> BezierReduce<F, [F; D]> for Bezier<F, D>
where
    F: Float,
{
    type Reduced = Self;
    type Quadratic = Self;
    type Cubic = Self;
    fn reduce(&self) -> Self::Reduced {
        match self.num {
            3 => [self.pts[0], self.pts[2]].into(),
            4 => self.as_array_4().reduced_to_quadratic().unwrap().into(),
            _ => *self,
        }
    }
    fn can_reduce(&self) -> bool {
        self.num >= 3
    }
    fn closeness_sq_to_reduction(&self) -> Option<F> {
        match self.num {
            3 => Some(self.closeness_sq_to_line()),
            4 => Some(self.closeness_sq_to_quadratic()),
            _ => None,
        }
    }

    fn closeness_sq_to_quadratic(&self) -> F {
        if self.num == 4 {
            self.as_array_4().closeness_sq_to_quadratic()
        } else {
            F::ZERO
        }
    }
    fn closeness_sq_to_cubic(&self) -> F {
        F::ZERO
    }
    fn reduced_to_quadratic(&self) -> Option<Self::Quadratic> {
        if self.num < 4 {
            None
        } else {
            Some(self.reduce())
        }
    }
    fn reduced_to_cubic(&self) -> Option<Self::Cubic> {
        None
    }
}

impl<F, const D: usize> BezierSplit for Bezier<F, D>
where
    F: Float,
{
    /// Returns two Bezier's that split the curve at parameter t=0.5
    ///
    /// For quadratics the midpoint is 1/4(p0 + 2*c + p1)
    fn split(&self) -> (Self, Self) {
        match self.num {
            2 => {
                let (b0, b1) = self.as_array_2().split();
                (b0.into(), b1.into())
            }
            3 => {
                let (b0, b1) = self.as_array_3().split();
                (b0.into(), b1.into())
            }
            _ => {
                let (b0, b1) = self.as_array_4().split();
                (b0.into(), b1.into())
            }
        }
    }
}

impl<F, const D: usize> BezierSection<F> for Bezier<F, D>
where
    F: Float,
{
    fn split_at(&self, t: F) -> (Self, Self) {
        (self.section(F::ZERO, t), self.section(t, F::ONE))
    }
    fn section(&self, t0: F, t1: F) -> Self {
        match self.num {
            2 => self.as_array_2().section(t0, t1).into(),
            3 => self.as_array_3().section(t0, t1).into(),
            _ => self.as_array_4().section(t0, t1).into(),
        }
    }
}

impl<F: Float, const D: usize> BezierOps<F, [F; D]> for Bezier<F, D> {
    fn add(&mut self, other: &Self) -> bool {
        for (s, o) in self.pts.iter_mut().zip(other.pts.iter()) {
            *s = vector::add(*s, o, F::ONE)
        }
        true
    }
    fn sub(&mut self, other: &Self) -> bool {
        for (s, o) in self.pts.iter_mut().zip(other.pts.iter()) {
            *s = vector::sub(*s, o, F::ONE)
        }
        true
    }
    fn scale(&mut self, scale: F) {
        for s in self.pts.iter_mut() {
            *s = vector::scale(*s, scale);
        }
    }
    fn map_pts(&mut self, map: &dyn Fn(usize, &[F; D]) -> [F; D]) {
        for (i, p) in self.pts.iter_mut().take(self.num).enumerate() {
            *p = map(i, p);
        }
    }
}

impl<F: Float, const D: usize> BezierMinMax<F> for Bezier<F, D> {
    fn t_coord_at_min_max(&self, use_max: bool, pt_index: usize) -> Option<(F, F)> {
        match self.num {
            2 => self.as_array_2().t_coord_at_min_max(use_max, pt_index),
            3 => self.as_array_3().t_coord_at_min_max(use_max, pt_index),
            _ => self.as_array_4().t_coord_at_min_max(use_max, pt_index),
        }
    }
}

impl<F, const D: usize> Bezier<F, D>
where
    F: Float,
{
    fn as_array_2(&self) -> [[F; D]; 2] {
        assert!(self.num == 2);
        self.pts[0..2].try_into().unwrap()
    }

    fn as_array_3(&self) -> [[F; D]; 3] {
        assert!(self.num == 3);
        self.pts[0..3].try_into().unwrap()
    }

    fn as_array_4(&self) -> [[F; D]; 4] {
        assert!(self.num == 4);
        self.pts[0..4].try_into().unwrap()
    }

    //mp scale
    /// Scale the Bezier by applying the scale factor to all of the points
    ///
    /// This is an example of the [Bezier::map_pts] method
    pub fn scale(&mut self, s: F) {
        self.map_pts(|p| vector::scale(p, s));
    }

    //mp map_pts
    /// Apply a function to all of the points in the Bezier
    pub fn map_pts<Map: Fn([F; D]) -> [F; D]>(&mut self, map: Map) {
        for p in self.pts.iter_mut() {
            *p = map(*p);
        }
    }

    //mp length
    /// Calculates the length of the Bezier when it is rendered down
    /// to the given a straightness
    ///
    /// `straightness` is independent of the length of the Bezier
    pub fn length(&self, straightness: F) -> F {
        self.as_lines(straightness * straightness)
            .fold(F::ZERO, |acc, (p0, p1)| {
                acc + geo_nd::vector::distance(&p0, &p1)
            })
    }

    //fi t_of_distance_rec
    /// Internal function used to find the distance recursively
    fn t_of_distance_rec(
        &self,
        straightness: F,
        distance: F,
        t_start: F,
        t_scale: F,
        acc_length: F,
    ) -> (Option<F>, F) {
        let zero = F::zero();
        let two = (2.0_f32).into();
        if distance <= acc_length {
            (Some(t_start), zero)
        } else if self.is_straight(straightness) {
            let (ep0, ep1) = self.endpoints();
            let d = vector::distance(ep0, ep1);
            if distance > acc_length + d {
                (None, acc_length + d)
            } else if d < F::epsilon() {
                (Some(t_start + t_scale), acc_length + d)
            } else {
                let rel_d = distance - acc_length;
                (Some(t_start + t_scale * rel_d / d), acc_length + d)
            }
        } else {
            let t_subscale = t_scale / two;
            let (b0, b1) = self.split();
            match b0.t_of_distance_rec(straightness, distance, t_start, t_subscale, acc_length) {
                (None, length) => b1.t_of_distance_rec(
                    straightness,
                    distance,
                    t_start + t_subscale,
                    t_subscale,
                    length,
                ),
                r => r,
            }
        }
    }

    //mp t_of_distance
    /// Calculates the parameter 't' at a certain distance along the Bezier given a straightness
    ///
    /// `straightness` is independent of the length of the Bezier
    ///
    /// Returns t,true if the distance is along the Bezier
    ///
    /// Returns 0.,false if the distance is before the start of the Bezier
    ///
    /// Returns 1.,false if the distance is beyond the end of the Bezier
    pub fn t_of_distance(&self, straightness: F, distance: F) -> (F, bool) {
        let zero = F::zero();
        let one = F::one();
        if distance < zero {
            (zero, false)
        } else {
            match self.t_of_distance_rec(straightness, distance, zero, one, zero) {
                (None, _) => (one, false),
                (Some(t), _) => (t, true),
            }
        }
    }

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
        bernstein_fns::arc::arc(angle, radius, center, unit, normal, rotate).into()
    }

    //fp of_round_corner
    /// Create a Cubic Bezier that is a circular arc focused on the corner point,
    /// with v0 and v1 are vectors IN to the point (P)
    ///
    /// As it is a circular arc we have a kite P, P+k.v0, C, P+k.v1, where
    ///
    /// ```text
    /// |P+k.v0 - C| = |P+k.v1 - C| = r; |P-C| = d (i.e. side lengths are r, r, k, k)
    /// ```
    ///
    /// with two corners being right-angles. (and d is the length of
    /// the kite diagonal opposite these right-angles).
    ///
    /// The kite is formed from two d, r, k right-angled triangles; it
    /// has two other angles, alpha and 2*theta, (alpha = angle
    /// between v0 and v1). Hence alpha = 180 - 2*theta, theta = 90-(alpha/2)
    ///
    /// Hence d^2 = r^2 + k^2; r/d = cos(theta), k/d=sin(theta)
    ///
    /// We know cos(alpha) = v0.v1 (assuming unit vectors).
    ///
    /// ```text
    /// cos(alpha) = cos(180-2*theta)
    ///            = -cos(2*theta)
    ///            = -(2cos^2(theta) - 1)
    ///            = 1 - 2cos^2(theta)
    ///
    /// cos^2(theta) = (1 - cos(alpha)) / 2 = r^2/d^2
    ///
    /// sin^2(theta) = (1 + cos(alpha)) / 2
    ///
    /// => d^2 = 2*r^2  / (1 - cos(alpha))
    /// ```
    ///
    /// Hence also k^2, and hence d and k.
    ///
    /// Then we require an arc given the angle of the arc is 2*theta
    pub fn of_round_corner(corner: &[F; D], v0: &[F; D], v1: &[F; D], radius: F) -> Self {
        bernstein_fns::arc::of_round_corner(corner, v0, v1, radius).into()
    }

    //mp center_radius_of_bezier_arc
    /// Find the center and radius of a Bezier if it is assumed to
    /// be a circular arc
    ///
    /// what is the center of the circle
    /// given point p0 and unit tangent t0
    /// and point p1 and unit tangent t1
    ///
    /// ```text
    /// |p0-c|^2 = |p1-c|^2 = r^2
    /// (p0-c) . t0 = 0
    /// (p1-c) . t1 = 0
    /// ```
    ///
    /// Consider c = k0.t0 + k1.t1
    ///
    /// (given t0.t0 == 1 and t1.t1==1)
    ///
    /// ```text
    /// (p0-c) . t0 = (p0 - k0.t0 - k1.t1).t0 = 0
    ///       p0.t0 = k0 + k1(t1.t0)
    /// ```
    /// similarly
    /// ```text
    ///       p1.t1 = k1 + k0(t1.t0)
    /// ```
    ///
    /// hence
    /// ```text
    ///  (t1.t0) * (p1.t1)         = k0.(t1.t0)^2 + k1(t1.t0)
    ///  p0.t0 - (t1.t0) * (p1.t1) = k0 ( 1 - (t1.t0)^2)
    ///  k0 = (p0.t0 - p1.t1 * t1.t0) / ( 1 - (t1.t0)^2)
    ///  k1 = (p1.t1 - p0.t0 * t1.t0) / ( 1 - (t1.t0)^2)
    /// ```
    pub fn center_radius_of_bezier_arc(&self) -> ([F; D], F) {
        bernstein_fns::arc::center_radius_of_bezier_arc(self)
    }

    //mp pt_distance_sq_from
    /// Calculate the distance of a point from a line segment of this kind of Bezier
    ///
    /// Also return whether the point is beyond the ends of the line (i.e. the distance
    /// is measured from the point to one of the endpoints)
    pub fn pt_distance_sq_from(pt: &[F; D], l0: &[F; D], l1: &[F; D]) -> (F, bool) {
        // Make everything relative to line.0
        //
        // so point is P and line is L (i.e. on the line = k*L for some 0 <= k <= 1)
        //
        // Note P = kL + N for some k and normal to the line N (N.L = 0)
        //
        // P.L = k(L.L) + N.L = k|L|^2
        let pt_m_l0 = vector::sub(*pt, l0, F::one());
        // len_pt_m_l0_sq = |P|^2
        let len_pt_m_l0_sq = vector::length_sq(&pt_m_l0);
        if len_pt_m_l0_sq < F::epsilon() {
            return (len_pt_m_l0_sq, false);
        }
        let l1_m_l0 = vector::sub(*l1, l0, F::one());
        // len_l1_m_l0_sq = |L|^2
        let len_l1_m_l0_sq = vector::length_sq(&l1_m_l0);
        if len_l1_m_l0_sq < F::epsilon() {
            return (len_pt_m_l0_sq, true);
        }
        // pt_along_line = P . L = k|L|^2
        let pt_along_line = vector::dot(&pt_m_l0, &l1_m_l0);

        // / len_pt_m_l0_sq.sqrt();
        // If k < 0 then distance is from line 0
        if pt_along_line < F::zero() {
            return (len_pt_m_l0_sq, true);
        }
        // If k > 1 then distance is from line 1
        if pt_along_line > len_l1_m_l0_sq {
            return (vector::distance_sq(pt, l1), true);
        }
        let pt_along_line_rel = pt_along_line / len_l1_m_l0_sq;
        let pt_projected = vector::scale(l1_m_l0, pt_along_line_rel);
        (vector::distance_sq(&pt_m_l0, &pt_projected), false)
    }
    //zz All done
}

// Deprecated methods
impl<F, const D: usize> Bezier<F, D>
where
    F: Float,
{
    //fp line
    /// Deprecated - use .into() with `From<[[F;D];2]>` trait
    ///
    /// Create a new Bezier that is a line between two points
    pub fn line(p0: &[F; D], p1: &[F; D]) -> Self {
        [*p0, *p1].into()
    }

    //fp quadratic
    /// Deprecated - use .into() with `From<[[F;D];3]>` trait
    ///
    /// Create a new Quadratic Bezier that is a line between two points
    /// with one absolute control points
    pub fn quadratic(p0: &[F; D], c: &[F; D], p1: &[F; D]) -> Self {
        [*p0, *c, *p1].into()
    }

    //fp cubic
    /// Deprecated - use .into() with `From<[[F;D];4]>` trait
    ///
    /// Create a new Cubic Bezier that is a line between two points
    /// with two absolute control points
    pub fn cubic(p0: &[F; D], c0: &[F; D], c1: &[F; D], p1: &[F; D]) -> Self {
        [*p0, *c0, *c1, *p1].into()
    }

    /// Deprecated - use BezierSection trait
    ///
    /// Returns the Bezier that is a subset of this Bezier between two parameters 0 <= t0 < t1 <= 1
    ///
    pub fn bezier_between(&self, t0: F, t1: F) -> Self {
        self.section(t0, t1)
    }

    //mp as_lines
    /// Deprecated - use BezierIntoIterator trait
    ///
    /// Return a [BezierLineIter] iterator that provides line segments
    /// when the Bezier is broken down into 'straight' enough through
    /// bisection.
    ///
    pub fn as_lines(&self, straightness_sq: F) -> impl Iterator<Item = ([F; D], [F; D])> + '_ {
        <Self as BezierIntoIterator<F, [F; D]>>::as_lines(self, straightness_sq)
    }

    //mp as_points
    /// Deprecated - use BezierIntoIterator trait
    ///
    /// Return a [BezierPointIter] iterator that provides points along
    /// the curve when the Bezier is broken down into 'straight'
    /// enough through bisection.
    ///
    pub fn as_points(&self, straightness_sq: F) -> impl Iterator<Item = [F; D]> + '_ {
        <Self as BezierIntoIterator<F, [F; D]>>::as_points(self, straightness_sq)
    }

    //mp is_straight
    /// Deprecated - use BezierEval and closeness_sq_to_line < straightness
    ///
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
        // p is the vector between two endpoints
        // lp2 is the length squared of p
        // c is a relative control point (i.e. mid-control relative to an endpoint)
        //
        // Return (effectively) (|c||p|sin(angle between p and c))^2 and |p|^2
        //
        // These two form a ration that should be used to reflect ( |c|sin(angle) )^2
        //
        // If |c| is tiny then return |c|sin(angle) of 0
        //
        // If |p| is tiny then return (|c|sin(90))^2
        fn straightness2_of_control<F, const D: usize>(p: &[F; D], lp2: F, c: &[F; D]) -> (F, F)
        where
            F: Float,
        {
            let lc2 = vector::length_sq(c);
            if lc2 < F::epsilon() {
                (F::zero(), lp2)
            } else if lp2 < F::epsilon() {
                (lc2, F::one())
            } else {
                // cdp is |c| |p| cos(angle between)
                let cdp = vector::dot(c, p);
                // c_p_s = |c|^2 |p|^2 * (1 - cos^2(angle between))
                //       = |c|^2 |p|^2 * sin^2(angle between)
                let c_p_s = lp2 * lc2 - cdp * cdp;
                (c_p_s, lp2)
            }
        }
        match self.num {
            2 => true,
            3 => {
                // Now make everything relative to p0 and test perpendicular distance of c from the line
                let c = vector::sub(self.pts[1], &self.pts[0], F::one());
                let p = vector::sub(self.pts[2], &self.pts[0], F::one());
                let lp2 = vector::length_sq(&p);

                // Need to test if (c-p0).(p1-p0) < -straightness.|p1-p0| or (c-p1).(p1-p0)>straightness.|p1-p0|
                let c_p1 = vector::sub(self.pts[1], &self.pts[2], F::one());
                if vector::dot(&c, &p) < -straightness * lp2.sqrt() {
                    return false;
                }
                if vector::dot(&c_p1, &p) > straightness * lp2.sqrt() {
                    return false;
                }
                // get |c||p|sin(angle) ^2 and |p|^2
                let (c_p_s_sq, p_sq) = straightness2_of_control(&p, lp2, &c);
                // return true if (|c|sin(angle))^2 *|p|^2 <= straightness^2 *|p|^2
                // i.e. |c|.|sin(angle)| < straightness
                c_p_s_sq <= straightness * straightness * p_sq
            }
            _ => {
                let p = vector::sub(self.pts[3], &self.pts[0], F::one());
                let c0 = vector::sub(self.pts[1], &self.pts[0], F::one());
                let c1 = vector::sub(self.pts[2], &self.pts[0], F::one());
                let lp2 = vector::length_sq(&p);

                // Need to test if (c-p0).(p1-p0) < -straightness.|p1-p0| or (c-p1).(p1-p0)>straightness.|p1-p0|
                let c0_p1 = vector::sub(self.pts[1], &self.pts[3], F::one());
                if vector::dot(&c0, &p) < -straightness * lp2.sqrt() {
                    return false;
                }
                if vector::dot(&c0_p1, &p) > straightness * lp2.sqrt() {
                    return false;
                }
                let c1_p1 = vector::sub(self.pts[2], &self.pts[3], F::one());
                if vector::dot(&c1, &p) < -straightness * lp2.sqrt() {
                    return false;
                }
                if vector::dot(&c1_p1, &p) > straightness * lp2.sqrt() {
                    return false;
                }

                // get |c||p|sin(angle) ^2 and |p|^2 for each control point
                let (c0_p_s_sq, p_sq_0) = straightness2_of_control(&p, lp2, &c0);
                let (c1_p_s_sq, p_sq_1) = straightness2_of_control(&p, lp2, &c1);
                // return true if Sum( (|c|sin(angle))^2*|p|^2  ) <= straightness *|p|^2
                (c0_p_s_sq + c1_p_s_sq) <= straightness * straightness * F::max(p_sq_0, p_sq_1)
            }
        }
    }
}

//ip Bezier<F, 2>
impl<F> Bezier<F, 2>
where
    F: Float,
{
    /// Return the normal at a paramenter 't' (for 2D beziers only)
    pub fn normal_at(&self, t: F) -> [F; 2] {
        let (sc, x) = self.derivative_at(t);
        [x[1] * sc, x[0] * sc]
    }

    //mp min_distance_sq_from
    /// Calculate an estimate for the minimum distance of a point from this Bezier
    ///
    /// If the point is 'within' the Bezier control points then this is zero
    pub fn min_distance_sq_from(&self, pt: &[F; 2]) -> F {
        // return (Pt-P0) x (P1-P0), which is -ve if Pt-P0 is anticlockwise of P1-P0
        fn side_of_line<F: Float>(pt: &[F; 2], p0: &[F; 2], p1: &[F; 2]) -> F {
            (pt[0] - p0[0]) * (p1[1] - p0[1]) - (pt[1] - p0[1]) * (p1[0] - p0[0])
        }
        let distance_sq_from_p0_p1 = Bezier::pt_distance_sq_from(pt, &self.pts[0], &self.pts[1]).0;
        match self.num {
            2 => distance_sq_from_p0_p1,
            3 => {
                let ctl_side_of_line = side_of_line(&self.pts[2], &self.pts[0], &self.pts[1]);
                let p_side_of_p0_p2 = side_of_line(pt, &self.pts[0], &self.pts[2]);
                let p_side_of_p2_p1 = side_of_line(pt, &self.pts[2], &self.pts[1]);
                let p_side_of_p1_p0 = side_of_line(pt, &self.pts[1], &self.pts[0]);

                if ctl_side_of_line == F::zero() {
                    distance_sq_from_p0_p1
                } else if ctl_side_of_line < F::zero()
                    && p_side_of_p0_p2 >= F::zero()
                    && p_side_of_p2_p1 >= F::zero()
                    && p_side_of_p1_p0 >= F::zero()
                {
                    F::zero()
                } else if ctl_side_of_line > F::zero()
                    && p_side_of_p0_p2 <= F::zero()
                    && p_side_of_p2_p1 <= F::zero()
                    && p_side_of_p1_p0 <= F::zero()
                {
                    F::zero()
                } else {
                    // Outside the triangle of the bezier and its control points
                    //
                    // Distance is at least the distance to the triangle
                    let distance_sq_from_p0_c =
                        Bezier::pt_distance_sq_from(pt, &self.pts[0], &self.pts[2]).0;
                    let distance_sq_from_p1_c =
                        Bezier::pt_distance_sq_from(pt, &self.pts[1], &self.pts[2]).0;
                    distance_sq_from_p0_c.min(distance_sq_from_p1_c.min(distance_sq_from_p0_p1))
                }
            }
            _ => F::zero(),
        }
    }
}
