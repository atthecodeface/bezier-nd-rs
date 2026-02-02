//a Imports
use crate::Float;
use geo_nd::vector;

use crate::{
    BezierEval, BezierIntoIterator, BezierLineIter, BezierPointIter, BezierReduce, BezierSection,
    BezierSplit,
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
    /// Control points - endpoints are always 0 and 1
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
            vector::fmt(f, &self.pts[2])?;
        }
        if self.num > 3 {
            write!(f, ":")?;
            vector::fmt(f, &self.pts[3])?;
        }
        write!(f, "->")?;
        vector::fmt(f, &self.pts[1])
    }

    //zz All done
}

impl<F, const D: usize> BezierEval<F, [F; D]> for Bezier<F, D>
where
    F: Float,
{
    fn point_at(&self, t: F) -> [F; D] {
        let two: F = (2.0_f32).into();
        let three: F = (3.0_f32).into();
        let omt = F::one() - t;
        match self.num {
            2 => self.vector_of(&[omt, t], F::one()),
            3 => {
                let p0_sc = omt * omt;
                let c_sc = two * omt * t;
                let p1_sc = t * t;
                self.vector_of(&[p0_sc, p1_sc, c_sc], F::one())
            }
            _ => {
                let p0_sc = omt * omt * omt;
                let c0_sc = three * omt * omt * t;
                let c1_sc = three * omt * t * t;
                let p1_sc = t * t * t;
                self.vector_of(&[p0_sc, p1_sc, c0_sc, c1_sc], F::one())
            }
        }
    }
    fn derivative_at(&self, t: F) -> (F, [F; D]) {
        let one = F::one();
        let two: F = (2.0_f32).into();
        let three: F = (3.0_f32).into();
        let four: F = (4.0_f32).into();
        match self.num {
            2 => (one, self.vector_of(&[-one, one], one)),
            3 => {
                let p0_sc = t - one; // d/dt (1-t)^2
                let c_sc = one - two * t; // d/dt 2t(1-t)
                let p1_sc = t; // d/dt t^2
                (one, self.vector_of(&[p0_sc, p1_sc, c_sc], one))
            }
            _ => {
                let p0_sc = two * t - t * t - one; // d/dt (1-t)^3
                let c0_sc = three * t * t - four * t + one; // d/dt 3t(1-t)^2
                let c1_sc = two * t - three * t * t; // d/dt 3t^2(1-t)
                let p1_sc = t * t; // d/dt t^3
                (one, self.vector_of(&[p0_sc, p1_sc, c0_sc, c1_sc], one))
            }
        }
    }

    fn endpoints(&self) -> (&[F; D], &[F; D]) {
        (&self.pts[0], &self.pts[1])
    }

    fn closeness_sq_to_line(&self) -> F {
        match self.num {
            2 => F::ZERO,
            3 => [self.pts[0], self.pts[2], self.pts[1]].closeness_sq_to_line(),
            _ => [self.pts[0], self.pts[2], self.pts[3], self.pts[1]].closeness_sq_to_line(),
        }
    }
    fn dc_sq_from_line(&self) -> F {
        match self.num {
            2 => F::ZERO,
            3 => [self.pts[0], self.pts[2], self.pts[1]].dc_sq_from_line(),
            _ => [self.pts[0], self.pts[2], self.pts[3], self.pts[1]].dc_sq_from_line(),
        }
    }
    fn num_control_points(&self) -> usize {
        self.num
    }
    fn control_point(&self, n: usize) -> &[F; D] {
        &self.pts[n]
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
            3 => Self::line(&self.pts[0], &self.pts[2]),
            4 => {
                let [p0, c, p1] = [self.pts[0], self.pts[2], self.pts[3], self.pts[1]]
                    .reduced_to_quadratic()
                    .unwrap();
                Self::quadratic(&p0, &c, &p1)
            }
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
        if self.num <= 3 {
            F::ZERO
        } else {
            [self.pts[0], self.pts[2], self.pts[3], self.pts[1]].closeness_sq_to_quadratic()
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
    //mp bisect
    /// Returns two Bezier's that split the curve at parameter t=0.5
    ///
    /// For quadratics the midpoint is 1/4(p0 + 2*c + p1)
    fn split(&self) -> (Self, Self) {
        let zero = F::zero();
        let one = F::one();
        let two = (2.0_f32).into();
        let three: F = (3.0_f32).into();
        let four: F = (4.0_f32).into();
        let eight: F = (8.0_f32).into();
        match self.num {
            2 => {
                let pm = self.vector_of(&[one, one], two);
                (Self::line(&self.pts[0], &pm), Self::line(&pm, &self.pts[1]))
            }
            3 => {
                let c0 = self.vector_of(&[one, zero, one], two);
                let c1 = self.vector_of(&[zero, one, one], two);
                let pm = vector::add(c0, &c1, F::one());
                let pm = vector::reduce(pm, 2.0_f32.into());
                (
                    Self::quadratic(&self.pts[0], &c0, &pm),
                    Self::quadratic(&pm, &c1, &self.pts[1]),
                )
            }
            _ => {
                let pm = self.vector_of(&[one, one, three, three], eight);
                let c00 = self.vector_of(&[one, zero, one], two);
                let c01 = self.vector_of(&[one, zero, two, one], four);
                let c10 = self.vector_of(&[zero, one, one, two], four);
                let c11 = self.vector_of(&[zero, one, zero, one], two);
                (
                    Self::cubic(&self.pts[0], &c00, &c01, &pm),
                    Self::cubic(&pm, &c10, &c11, &self.pts[1]),
                )
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
            2 => {
                let u0 = F::one() - t0;
                let u1 = F::one() - t1;
                let r0 = self.vector_of(&[u0, t0], F::one());
                let r1 = self.vector_of(&[u1, t1], F::one());
                Self::line(&r0, &r1)
            }
            3 => {
                let two: F = (2.0_f32).into();
                let u0 = F::one() - t0;
                let u1 = F::one() - t1;
                let rp0 = self.vector_of(&[u0 * u0, t0 * t0, two * u0 * t0], F::one());
                let rp1 = self.vector_of(&[u1 * u1, t1 * t1, two * u1 * t1], F::one());
                let rc0 = self.vector_of(&[u0 * u1, t1 * t0, u0 * t1 + u1 * t0], F::one());
                Self::quadratic(&rp0, &rc0, &rp1)
            }
            _ => {
                // simply: c0 = p0 + tangent(0)
                // and if we scale the curve to t1-t0 in size, tangents scale the same
                let rp0 = self.point_at(t0);
                let rp1 = self.point_at(t1);
                let (_sc0, rt0) = self.derivative_at(t0);
                let (_sc1, rt1) = self.derivative_at(t1);

                let t1_m_t0 = t1 - t0;

                let mut rc0 = [F::zero(); D];
                let mut rc1 = [F::zero(); D];
                for i in 0..D {
                    rc0[i] = rp0[i] + t1_m_t0 * rt0[i]; // *sc0*sc1?
                    rc1[i] = rp1[i] - t1_m_t0 * rt1[i];
                }

                Self::cubic(&rp0, &rc0, &rc1, &rp1)
            }
        }
    }
}

impl<F, const D: usize> Bezier<F, D>
where
    F: Float,
{
    //fp line
    /// Create a new Bezier that is a line between two points
    pub fn line(p0: &[F; D], p1: &[F; D]) -> Self {
        Self {
            num: 2,
            pts: [*p0, *p1, [F::zero(); D], [F::zero(); D]],
        }
    }

    //fp quadratic
    /// Create a new Quadratic Bezier that is a line between two points
    /// with one absolute control points
    pub fn quadratic(p0: &[F; D], c: &[F; D], p1: &[F; D]) -> Self {
        Self {
            num: 3,
            pts: [*p0, *p1, *c, [F::zero(); D]],
        }
    }

    //fp cubic
    /// Create a new Cubic Bezier that is a line between two points
    /// with two absolute control points
    pub fn cubic(p0: &[F; D], c0: &[F; D], c1: &[F; D], p1: &[F; D]) -> Self {
        Self {
            num: 4,
            pts: [*p0, *p1, *c0, *c1],
        }
    }

    //mp elevate
    /// Elevate a Bezier by one degree (cannot elevate a Cubic)
    pub fn elevate(mut self) -> Self {
        match self.num {
            2 => {
                self.num = 3;
                self.pts[2] = vector::reduce(
                    vector::add(self.pts[0], &self.pts[1], F::one()),
                    (2.0_f32).into(),
                );
            }
            3 => {
                // Set pts[3] using pts[0..2] before pts[2] is changed...
                self.num = 4;
                self.pts[3] = vector::reduce(
                    vector::add(self.pts[1], &self.pts[2], (2.0_f32).into()),
                    (3.0_f32).into(),
                );
                self.pts[2] = vector::reduce(
                    vector::add(self.pts[0], &self.pts[2], (2.0_f32).into()),
                    (3.0_f32).into(),
                );
            }
            _ => {
                panic!("Cannot elevate a cubic Bezier at this point");
            }
        }
        self
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

    //mi vector_of
    /// Returns a vector of a combination of the vectors of the bezier
    #[inline]
    fn vector_of(&self, sc: &[F], reduce: F) -> [F; D] {
        let mut r = [F::zero(); D];
        for (i, sc) in sc.iter().enumerate() {
            for (j, rj) in r.iter_mut().enumerate() {
                *rj += *sc * self.pts[i][j];
            }
        }
        vector::reduce(r, reduce)
    }

    //mp bezier_between
    /// Returns the Bezier that is a subset of this Bezier between two parameters 0 <= t0 < t1 <= 1
    pub fn bezier_between(&self, t0: F, t1: F) -> Self {
        self.section(t0, t1)
    }

    //mp as_lines
    /// Return a [BezierLineIter] iterator that provides line segments
    /// when the Bezier is broken down into 'straight' enough through
    /// bisection.
    pub fn as_lines(&self, straightness_sq: F) -> impl Iterator<Item = ([F; D], [F; D])> + '_ {
        <Self as BezierIntoIterator<F, _, [F; D]>>::as_lines(self, straightness_sq)
    }

    //mp as_points
    /// Return a [BezierPointIter] iterator that provides points along
    /// the curve when the Bezier is broken down into 'straight'
    /// enough through bisection.
    pub fn as_points(&self, straightness_sq: F) -> impl Iterator<Item = [F; D]> + '_ {
        <Self as BezierIntoIterator<F, _, [F; D]>>::as_points(self, straightness_sq)
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
                let c = vector::sub(self.pts[2], &self.pts[0], F::one());
                let p = vector::sub(self.pts[1], &self.pts[0], F::one());
                let lp2 = vector::length_sq(&p);

                // Need to test if (c-p0).(p1-p0) < -straightness.|p1-p0| or (c-p1).(p1-p0)>straightness.|p1-p0|
                let c_p1 = vector::sub(self.pts[2], &self.pts[1], F::one());
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
                let p = vector::sub(self.pts[1], &self.pts[0], F::one());
                let c0 = vector::sub(self.pts[2], &self.pts[0], F::one());
                let c1 = vector::sub(self.pts[3], &self.pts[0], F::one());
                let lp2 = vector::length_sq(&p);

                // Need to test if (c-p0).(p1-p0) < -straightness.|p1-p0| or (c-p1).(p1-p0)>straightness.|p1-p0|
                let c0_p1 = vector::sub(self.pts[2], &self.pts[1], F::one());
                if vector::dot(&c0, &p) < -straightness * lp2.sqrt() {
                    return false;
                }
                if vector::dot(&c0_p1, &p) > straightness * lp2.sqrt() {
                    return false;
                }
                let c1_p1 = vector::sub(self.pts[3], &self.pts[1], F::one());
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

    //mp length
    /// Calculates the length of the Bezier when it is rendered down
    /// to the given a straightness
    ///
    /// `straightness` is independent of the length of the Bezier
    pub fn length(&self, straightness: F) -> F {
        if self.is_straight(straightness) {
            let (ep0, ep1) = self.endpoints();
            vector::distance(ep0, ep1)
        } else {
            let (b0, b1) = self.split();
            b0.length(straightness) + b1.length(straightness)
        }
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

    //fi lambda_of_k_d
    fn lambda_of_k_d(k: F, d: F) -> F {
        // There are numerous versions of calculating
        // the lambda for the arc from the angle of the arc
        //
        // For a 90 degree arc the *best* values is 0.2652165 apparently
        //
        // One equation that provides this  is
        //   lambda = four_thirds * radius * (d/k - one);
        //
        // Another is:
        //   theta = (k/d).asin() / F::int(4);
        //   lambda = four_thirds * theta.tan();
        //
        // This table is captures the values for this second
        //
        // Actually attempting a better approximation leads to the following for (k/d)^2 -> lambda
        //
        // 0.0011099165 0.009397572
        // 0.004424813 0.04415609
        // 0.008196682 0.06044403
        // 0.012195113 0.07385845
        // 0.019999988 0.09472578
        // 0.038461603 0.13201918
        // 0.100000046 0.21637033
        // 0.100000046 0.21637033
        // 0.137931 0.2567711
        // 0.20000009 0.314736
        // 0.3076923 0.40363038
        // 0.5 0.5519717
        // 0.6923078 0.71254206
        // 0.8000001 0.822074
        // 0.862069 0.89976513
        // 0.8999999 0.9571549
        // 0.9615385 1.0864261
        // 0.98 1.1479391
        // 0.99180335 1.2072284
        // 0.9955752 1.2359663
        // 0.9988901 1.2764238
        //
        // With a quintic polynomial of coeffs (x^0 + x^1 +... + x^5):
        // 3.1603235091816735e-002
        // 2.7950542994656820e+000
        // -1.1486743224812313e+001
        // 2.8975368657401102e+001
        // -3.2845222512637491e+001
        // 1.3779429574112177e+001
        //
        // Or for r^2/d^2 -> lambda
        // 0.9988901 0.009397572
        // 0.9955752 0.04415609
        // 0.99180335 0.06044403
        // 0.9878049 0.07385845
        // 0.98 0.09472578
        // 0.9615384 0.13201918
        // 0.9 0.21637033
        // 0.9 0.21637033
        // 0.862069 0.2567711
        // 0.79999995 0.314736
        // 0.6923077 0.40363038
        // 0.5 0.5519717
        // 0.3076923 0.71254206
        // 0.20000002 0.822074
        // 0.13793105 0.89976513
        // 0.10000005 0.9571549
        // 0.038461536 1.0864261
        // 0.02000001 1.1479391
        // 0.008196682 1.2072284
        // 0.0044247806 1.2359663
        // 0.0011098981 1.2764238
        //
        // 1.2494900596889080e+000
        // -4.2639321404424191e+000
        // 1.6162330324360198e+001
        // -3.5388797293367219e+001
        // 3.6051953254575963e+001
        // -1.3779440945693199e+001

        let k_d = k / d;
        // let four_thirds  = F::frac(4,3);
        // let  theta = (k/d).asin() / F::int(4);
        // let  lambda = four_thirds * theta.tan();
        // lambda
        let k_d = k_d * k_d;
        let a0: F = 3.160_323_6e-2_f32.into();
        let a1: F = 2.795_054_2_f32.into();
        let a2: F = (-1.148_674_3e1_f32).into();
        let a3: F = 2.897_536_8e1_f32.into();
        let a4: F = (-3.284_522_2e1_f32).into();
        let a5: F = 1.377_942_9e1_f32.into();
        a0 + a1 * k_d
            + a2 * k_d * k_d
            + a3 * k_d * k_d * k_d
            + a4 * k_d * k_d * k_d * k_d
            + a5 * k_d * k_d * k_d * k_d * k_d
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
        let two = (2.0_f32).into();
        let half_angle = angle / two;
        let s = half_angle.sin();
        let lambda = radius * Self::lambda_of_k_d(s, F::one());

        let d0a = rotate;
        let (d0s, d0c) = d0a.sin_cos();
        let d1a = rotate + angle;
        let (d1s, d1c) = d1a.sin_cos();

        let mut p0 = [F::zero(); D];
        let mut p1 = [F::zero(); D];
        let mut c0 = [F::zero(); D];
        let mut c1 = [F::zero(); D];
        for i in 0..D {
            p0[i] = center[i] + unit[i] * (d0c * radius) + normal[i] * (d0s * radius);
            p1[i] = center[i] + unit[i] * (d1c * radius) + normal[i] * (d1s * radius);

            c0[i] = p0[i] - unit[i] * (d0s * lambda) + normal[i] * (d0c * lambda);
            c1[i] = p1[i] + unit[i] * (d1s * lambda) - normal[i] * (d1c * lambda);
        }

        Self::cubic(&p0, &c0, &c1, &p1)
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
        let nearly_one = (0.999_999_f32).into();
        let one = F::one();
        let two: F = (2.0_f32).into();
        let v0 = vector::normalize(*v0);
        let v1 = vector::normalize(*v1);
        let cos_alpha = vector::dot(&v0, &v1);
        if cos_alpha.abs() >= nearly_one {
            // v0 and v1 point in the same direction
            let mut p0 = [F::zero(); D];
            let mut p1 = [F::zero(); D];
            for i in 0..D {
                p0[i] = corner[i] - radius * v0[i];
                p1[i] = corner[i] - radius * v1[i];
            }
            Self::quadratic(&p0, corner, &p1)
        } else {
            let r2 = radius * radius;
            let d2 = two * r2 / (one - cos_alpha);
            let k2 = d2 - r2;
            let d = d2.sqrt();
            let k = k2.sqrt();

            let lambda = radius * Self::lambda_of_k_d(k, d);

            let mut p0 = [F::zero(); D];
            let mut p1 = [F::zero(); D];
            let mut c0 = [F::zero(); D];
            let mut c1 = [F::zero(); D];
            for i in 0..D {
                p0[i] = corner[i] - k * v0[i];
                p1[i] = corner[i] - k * v1[i];
                c0[i] = p0[i] + lambda * v0[i];
                c1[i] = p1[i] + lambda * v1[i];
            }
            Self::cubic(&p0, &c0, &c1, &p1)
        }
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
        let zero = F::zero();
        let one = F::one();
        let p0 = self.point_at(zero);
        let p1 = self.point_at(one);
        let (_sc0, t0) = self.derivative_at(zero);
        let (_sc1, t1) = self.derivative_at(one);
        let t0 = vector::normalize(t0);
        let t1 = vector::normalize(t1);
        let t1_d_t0 = vector::dot(&t1, &t0);
        let p0_d_t0 = vector::dot(&p0, &t0);
        let p1_d_t1 = vector::dot(&p1, &t1);
        let k0 = (p0_d_t0 - p1_d_t1 * t1_d_t0) / (one - t1_d_t0 * t1_d_t0);
        let k1 = (p1_d_t1 - p0_d_t0 * t1_d_t0) / (one - t1_d_t0 * t1_d_t0);

        let mut c = [F::zero(); D];
        for i in 0..D {
            c[i] = t0[i] * k0 + t1[i] * k1;
        }

        let r = (vector::distance(&c, &p0) + vector::distance(&c, &p1)) / (2.0_f32).into();
        (c, r)
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
