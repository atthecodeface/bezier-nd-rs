/*a Copyright

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@file    curve.rs
@brief   Bezier curves of N-dimensions
 */

//a Imports
use geo_nd::vector;
use geo_nd::{Float, Vector};
use crate::BezierLineIter;
use crate::BezierPointIter;

//a Bezier
//tp Bezier
/// This library supports Bezier curves of up to order 3 - i.e. up to
/// Cubic; these have two control poits.
///
/// Note: in this section we use u=1-t
///
/// A linear Bezier has two points, p0 and p1, and provides points
/// along the line as:
///    p(t,u=1-t) = u*p0 + t*p1
///
/// A linear Bezier may be split at t into (p0, u*p0+t*p1); (u*p0+t*p1, p1).
///
/// A quadratic Bezier has three points, p0, c and p1, and provides
/// points along the curve as:
///
///    p(t,u=1-t) = u^2.p0 + 2.t.u.c + t^2.p1
///
/// or, viewing it is a linear Bezier between two linear beziers:
///
///    p(t) = u(u.p0 + t.c) + t(u.c + t.p1)
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
///    Note t2 = (t1-t0)/u0; u2=1-t2 = (u0+t0-t1)/u0 = (1-t1)/u0 = u1/u0
///
///    cl(t0) = u0.p0 + t0.c
///    cr(t0) = u0.c  + t1.p1
///     p(t0) = u0.cl(t0)  + t0.cr(t0)
///
///    Bezier t0->1 : p(t0), cr(t0), p1
///
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
///  Hence
///  c(t0,t1)  = u1.cl(t0) + t1.cr(t0)
///            = u0.u1.p0 + u1.t0.c + u0.t1.c + t0.t1.p1
///            = u0.u1.p0 + (u1.t0+u0.t1).c + t0.t1.p1
///  And the points are:
///      p(t0) = u0.u0.p0 + 2(u0.t0).c + t0.t0.p1
///      p(t1) = u1.u1.p0 + 2(u1.t1).c + t1.t1.p1
///
use std::marker::PhantomData;
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Bezier<F, V, const D:usize>
where F:Float, V:Vector<F,D> {
    /// Number of valid control points (2-4)
    num : usize,
    /// Control points - endpoints are always 0 and 1
    pts : [V;4],
    f:PhantomData<F>,
}

//ti Display for Bezier
impl <F, V, const D:usize> std::fmt::Display for Bezier<F, V, D>
where F:Float, V:Vector<F,D> {

    //mp fmt - format a `Bezier` for display
    /// Display the `Bezier' as sets of points
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[" )?;
        vector::fmt(f,self.pts[0].as_ref())?;
        write!(f, "<-" )?;
        if self.num > 2 {vector::fmt(f,self.pts[2].as_ref())?;}
        if self.num > 3 {
            write!(f, ":" )?;
            vector::fmt(f,self.pts[3].as_ref())?;
        }
        write!(f, "->" )?;
        vector::fmt(f,self.pts[1].as_ref())
    }

    //zz All done
}

//ip Bezier
impl <F, V, const D:usize> Bezier<F, V, D>
where F:Float, V:Vector<F,D> {
    //mp borrow_pt
    /// Get the start or end point of the Bezier - index 0 gives the
    /// start point, index 1 the end point
    pub fn borrow_pt(&self, index:usize) -> &V {
        &self.pts[index]
    }

    //mp endpoints
    /// Deconstruct and get the endpoints
    pub fn endpoints(self) -> (V, V) {
        (self.pts[0], self.pts[1])
    }

    //mp get_distance
    /// Get the distance between the start and end points
    pub fn get_distance(&self) -> F {
        self.pts[0].distance(&self.pts[1])
    }

    //fp line
    /// Create a new Bezier that is a line between two points
    pub fn line(p0:&V, p1:&V) -> Self {
        Self { num:2, pts:[*p0, *p1, V::zero(), V::zero()], f:PhantomData }
    }

    //fp quadratic
    /// Create a new Quadratic Bezier that is a line between two points
    /// with one absolute control points
    pub fn quadratic(p0:&V, c:&V, p1:&V) -> Self {
        Self { num:3, pts:[*p0, *p1, *c, V::zero()], f:PhantomData }
    }

    //fp cubic
    /// Create a new Cubic Bezier that is a line between two points
    /// with two absolute control points
    pub fn cubic(p0:&V, c0:&V, c1:&V, p1:&V) -> Self {
        Self { num:4, pts:[*p0, *p1, *c0, *c1], f:PhantomData }
    }

    //mp is_line
    /// Returns true if the Bezier is a straight line
    pub fn is_line(&self) -> bool { self.num == 2 }

    //mp is_quadratic
    /// Returns true if the Bezier is a quadratic
    pub fn is_quadratic(&self) -> bool { self.num == 3 }

    //mp is_cubic
    /// Returns true if the Bezier is a cubic
    pub fn is_cubic(&self) -> bool { self.num == 4 }

    //mp scale
    /// Consume the Bezier and return a new Bezier scaled separately in X and Y by two scaling parameters
    pub fn scale(&mut self, s:F)  {
        //for i in 0..self.pts.len() {
        //self.pts[i] *= s;
        //}
        self.map_pts(|p| p*s);
    }

    //mp map_pts
    /// Apply a function to all of the points in the Bezier
    pub fn map_pts<Map:Fn(V)->V> (&mut self, map:Map) {
        for p in self.pts.iter_mut() {
            *p = map(*p);
        }
    }

    //mp vector_of
    /// Returns a vector of a combination of the vectors of the bezier
    #[inline]
    pub fn vector_of(&self, sc:&[F], reduce:F) -> V {
        let mut r = self.pts[0] * sc[0];
        for i in 1..sc.len() {
            r += self.pts[i] * sc[i];
        }
        r / reduce
    }

    //mp point_at
    /// Returns the point at parameter 't' along the Bezier
    pub fn point_at(&self, t:F) -> V {
        let omt = F::int(1) - t;
        match self.num {
            2 => {
                self.vector_of(&[omt, t], F::int(1))
            },
            3 => {
                let p0_sc =             omt*omt;
                let c_sc  = F::int(2) * omt*t;
                let p1_sc =             t*t;
                self.vector_of(&[p0_sc, p1_sc, c_sc], F::int(1))
            },
            _  => {
                let p0_sc =             omt*omt*omt;
                let c0_sc = F::int(3) * omt*omt*t;
                let c1_sc = F::int(3) * omt*t*t;
                let p1_sc =             t*t*t;
                self.vector_of(&[p0_sc, p1_sc, c0_sc, c1_sc], F::int(1))
            },
        }
    }

    //mp tangent_at
    /// Returns the tangent vector at parameter 't' along the Bezier
    ///
    /// This is not necessarily a unit vector
    pub fn tangent_at(&self, t:F) -> V {
        let one   = F::int(1);
        let two   = F::int(2);
        let three = F::int(3);
        let four  = F::int(4);
        match self.num {
            2 => {
                self.vector_of(&[-one, one], one)
            },
            3 => {
                let p0_sc = t - one;    // d/dt (1-t)^2
                let c_sc  = one - two*t; // d/dt 2t(1-t)
                let p1_sc = t;         // d/dt t^2
                self.vector_of(&[p0_sc, p1_sc, c_sc], one)
            },
            _ => {
                let p0_sc = two*t     - t*t - one ; // d/dt (1-t)^3
                let c0_sc = three*t*t - four*t + one; // d/dt 3t(1-t)^2
                let c1_sc = two*t     - three*t*t   ; // d/dt 3t^2(1-t)
                let p1_sc = t*t               ; // d/dt t^3
                self.vector_of(&[p0_sc, p1_sc, c0_sc, c1_sc], one)
            },
        }
    }

    //mp bisect
    /// Returns two Bezier's that split the curve at parameter t=0.5
    ///
    /// For quadratics the midpoint is 1/4(p0 + 2*c + p1)
    pub fn bisect(&self) -> (Self, Self) {
        let zero   = F::zero();
        let one    = F::int(1);
        let two    = F::int(2);
        let three  = F::int(3);
        let four   = F::int(4);
        let eight  = F::int(8);
        match self.num {
            2 => {
                let pm = self.vector_of(&[one,one],two);
                (Self::line(&self.pts[0], &pm), Self::line(&pm, &self.pts[1]))
            },
            3 => {
                let c0 = self.vector_of(&[one,zero,one],two);
                let c1 = self.vector_of(&[zero,one,one],two);
                let pm = (c0 + c1) / two;
                (Self::quadratic(&self.pts[0], &c0, &pm), Self::quadratic(&pm, &c1, &self.pts[1]))
            },
            _ => {
                let pm  = self.vector_of(&[one,one,three,three],eight);
                let c00 = self.vector_of(&[one,zero,one],two);
                let c01 = self.vector_of(&[one,zero,two,one],four);
                let c10 = self.vector_of(&[zero,one,one,two],four);
                let c11 = self.vector_of(&[zero,one,zero,one],two);
                (Self::cubic(&self.pts[0],&c00,&c01,&pm), Self::cubic(&pm,&c10,&c11,&self.pts[1]))
            },
        }
    }

    //mp bezier_between
    /// Returns the Bezier between two parameters 0 <= t0 < t1 <= 1
    pub fn bezier_between(&self, t0:F, t1:F) -> Self {
        let two = F::int(2);
        let p0 = &self.pts[0];
        let p1 = &self.pts[1];
        match self.num {
            2 => {
                let u0 = F::one() - t0;
                let u1 = F::one() - t1;
                let r0 = self.pts[0]*u0 + self.pts[1] * t0;
                let r1 = self.pts[0]*u1 + self.pts[1] * t1;
                Self::line(&r0, &r1)
            },
            3 => {
                let c = &self.pts[2];
                let u0 = F::one() - t0;
                let u1 = F::one() - t1;
                let rp0 = *p0 * (u0*u0) + *c * (two*u0*t0)     + *p1 * (t0*t0);
                let rp1 = *p0 * (u1*u1) + *c * (two*u1*t1)     + *p1 * (t1*t1);
                let rc0 = *p0 * (u0*u1) + *c * (u0*t1 + u1*t0) + *p1 * (t1*t0);
                Self::quadratic(&rp0, &rc0, &rp1)
            },
            _ => {
                // simply: c0 = p0 + tangent(0)
                // and if we scale the curve to t1-t0 in size, tangents scale the same
                let rp0 = self.point_at(t0);
                let rt0 = self.tangent_at(t0);
                let rt1 = self.tangent_at(t1);
                let rp1 = self.point_at(t1);
                let t1_m_t0 = t1 - t0;
                let rc0 = rp0 + rt0 * t1_m_t0;
                let rc1 = rp1 - rt1 * t1_m_t0;
                Self::cubic(&rp0, &rc0, &rc1, &rp1)
            },
        }
    }

    //mp as_lines
    /// Iterate over line segments that are 'straight' enough
    pub fn as_lines(&self, straightness:F) -> BezierLineIter<F, V, D> {
        BezierLineIter::new(self, straightness)
    }

    //mp as_points
    /// Iterate over points that make 'straight' enough lines
    pub fn as_points(&self, straightness:F) -> BezierPointIter<F, V, D> {
        BezierPointIter::new(BezierLineIter::new(self, straightness))
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
    pub fn is_straight(&self, straightness:F) -> bool {
        fn straightness_of_control<F, V, const D:usize>(p:&V, lp2:F, c:&V) -> (F,F)
        where F:Float, V:Vector<F,D> {
            let lc2 = c.length_sq();
            if lc2 < F::epsilon() {
                (F::zero(), lp2)
            } else if lp2 < F::epsilon() {
                (lc2, F::one())
            } else {
                let cdp = c.dot(p);
                let c_s = F::sqrt(lp2*lc2 - cdp*cdp);
                (c_s, lp2)
            }
        }
        match self.num {
            2 => true,
            3 => {
                let p = self.pts[1] - self.pts[0];
                let lp2 = p.length_sq();
                let c = self.pts[2] - self.pts[0];
                let (c_s, sc) = straightness_of_control(&p, lp2, &c);
                c_s <= straightness * sc
            },
            _ => {
                let p = self.pts[1] - self.pts[0];
                let lp2 = p.length_sq();
                let c0 = self.pts[2] - self.pts[0];
                let (c0_s, sc0) = straightness_of_control(&p, lp2, &c0);
                let c1 = self.pts[3] - self.pts[0];
                let (c1_s, sc1) = straightness_of_control(&p, lp2, &c1);
                (c0_s + c1_s) <= straightness * F::max(sc0, sc1)
            },
        }
    }

    //mp length
    /// Calculates the length given a straightness
    ///
    /// `straightness` is independent of the length of the Bezier
    pub fn length(&self, straightness:F) -> F {
        if self.is_straight(straightness) {
            self.get_distance()
        } else {
            let (b0, b1) = self.bisect();
            b0.length(straightness) + b1.length(straightness)
        }
    }

    //mp t_of_distance
    /// Calculates the parameter 't' at a certain distance along the Bezier given a straightness
    ///
    /// `straightness` is independent of the length of the Bezier
    ///
    /// Returns t,true if the distance is along the Bezier
    /// Returns 0.,false if the distance is before the start of the Bezier
    /// Returns 1.,false if the distance is beyond the end of the Bezier
    fn t_of_distance_rec(&self, straightness:F, distance:F, t_start:F, t_scale:F, acc_length:F) -> (Option<F>, F) {
        let zero   = F::zero();
        let two    = F::int(2);
        if distance <= acc_length {
            (Some(t_start), zero)
        } else if self.is_straight(straightness) {
            let d     = self.get_distance();
            if distance > acc_length+d {
                (None, acc_length+d)
            } else if d < F::epsilon() {
                (Some(t_start + t_scale), acc_length+d)
            } else {
                let rel_d = distance - acc_length;
                (Some(t_start + t_scale * rel_d / d), acc_length+d)
            }
        } else {
            let t_subscale = t_scale / two;
            let (b0, b1) = self.bisect();
            match b0.t_of_distance_rec(straightness, distance, t_start, t_subscale, acc_length) {
                (None, length) => {
                    b1.t_of_distance_rec( straightness, distance, t_start + t_subscale, t_subscale, length )
                }
                r => r
            }
        }
    }
    pub fn t_of_distance(&self, straightness:F, distance:F) -> (F, bool) {
        let zero   = F::zero();
        let one    = F::int(1);
        if distance < zero {
            (zero, false)
        } else {
            match self.t_of_distance_rec(straightness, distance, zero, one, zero) {
                (None, _)    => (one, false),
                (Some(t), _) => (t, true),
            }
        }
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
    pub fn arc(angle:F, radius:F, center:&V, unit:&V, normal:&V, rotate:F) -> Self {
        let one    = F::int(1);
        let two    = F::int(2);
        let four_thirds  = F::frac(4,3);
        let half_angle = angle / two;
        let s = half_angle.sin();
        let lambda = radius * four_thirds * (one/s - one);

        let d0a = rotate;
        let (d0s,d0c) = d0a.sin_cos();
        let d1a = rotate+angle;
        let (d1s,d1c) = d1a.sin_cos();

        let p0 = *center + (*unit) * (d0c*radius) + *normal * (d0s*radius);
        let p1 = *center + (*unit) * (d1c*radius) + *normal * (d1s*radius);

        let c0 = p0 - (*unit) * (d0s*lambda) + *normal * (d0c*lambda);
        let c1 = p1 + (*unit) * (d1s*lambda) - *normal * (d1c*lambda);

        Self::cubic(&p0, &c0, &c1, &p1)
    }

    //fp of_round_corner
    /// Create a Cubic Bezier that is a circular arc focused on the corner point,
    /// with v0 and v1 are vectors IN to the point (P)
    ///
    /// As it is a circular arc we have a kite P, P+k.v0, C, P+k.v1, where
    ///
    /// |P+k.v0 - C| = |P+k.v1 - C| = r; |P-C| = d (i.e. side lengths are r, r, k, k)
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
    ///
    /// Hence also k^2, and hence d and k.
    ///
    /// Then we require an arc given the angle of the arc is 2*theta, which requires a lambda of
    /// 4/3 * r * (1/sin(theta)-1) = 4/3 * r * (d/k - 1)
    ///
    /// Note though that d^2/k^2 = 1/sin^2(theta) = 2/(1+cos(alpha))
    ///
    /// hence d/k = sqrt(2/(1+cos(alpha)))
    ///
    /// hence lambda = 4/3 * r * (sqrt(2/(1+cos(alpha))) - 1)
    pub fn of_round_corner(corner:&V, v0:&V, v1:&V, radius:F) -> Self {
        let nearly_one = F::frac(99_999, 100_000);
        let one    = F::int(1);
        let two    = F::int(2);
        let four_thirds  = F::frac(4,3);
        let mut v0    = v0.clone();
        let mut v1    = v1.clone();
        v0.normalize();
        v1.normalize();
        println!("{}, {}",v0,v1);
        // let reverse = v0[0]*v1[1] - v1[0]*v0[1] > F::zero();
        // println!("Round rev {}",reverse);
        let cos_alpha = v0.dot(&v1);
        if cos_alpha >= nearly_one {
            // v0 and v1 point in the same direction
            let p0 = *corner - (v0 * radius);
            let p1 = *corner - (v1 * radius);
            Self::quadratic(&p0, corner, &p1)
        } else if cos_alpha <= -nearly_one {
            // basically 180 degress apart
            let p0 = *corner - (v0 * radius);
            let p1 = *corner - (v1 * radius);
            Self::quadratic(&p0, corner, &p1)
        } else {
            let r2 = radius * radius;
            let d2 = two * r2 / (one - cos_alpha);
            let k2 = d2 - r2;
            let d = d2.sqrt();
            let k = k2.sqrt();
            /*
            println!("{} r^2 {} d^2 {} {} {} {}", cos_alpha, r2, d2, k2, d, k);
            println!("alpha {} theta {}, sin(theta) {} : {}",
                     cos_alpha.acos().to_degrees(),
                     (F::int(180) - cos_alpha.acos().to_degrees())/F::int(2),
                     ((F::int(180) - cos_alpha.acos().to_degrees())/F::int(2)).to_radians().sin(),
                     k/d,
            );
             */
            // let v0_plus_v1_u = vector::normalize(vector::add(v0.clone(), &v1, one));
            // let center = vector::add(corner.clone(), &v0_plus_v1_u, d);
            // let lambda = four/three * radius * ((two / (one + cos_alpha)).sqrt() - one);

            // let lambda = four_thirds * radius * (d/k - one);
            let theta = (k/d).asin() / F::int(4);
            let lambda = four_thirds * theta.tan();
            // println!("lambda {} : {} ", lambda, theta);

            let p0 = *corner - (v0 * k);
            let p1 = *corner - (v1 * k);
            // let k = radius / (n1u.dot(&v0u));
            // let center = Point::new( corner.x-k*(v0u.x+v1u.x), corner.y-k*(v0u.y+v1u.y) );
            // let normal_diff = Point::new(n0u.x-n1u.x, n0u.y-n1u.y);
            // let vector_sum  = Point::new(v0u.x+v1u.x, v0u.y+v1u.y);
            // let l2 = vector_sum.x*vector_sum.x + vector_sum.y*vector_sum.y;
            // let l = l2.sqrt();
            // let lambda = 4.0*radius/(3.*l2) * (2.*l + (normal_diff.x*vector_sum.x + normal_diff.y*vector_sum.y));
            // let c0 = Point::new(p0.x + lambda * v0u.x, p0.y + lambda * v0u.y);
            // let c1 = Point::new(p1.x + lambda * v1u.x, p1.y + lambda * v1u.y);
            let c0 = p0 + (v0 * lambda);
            let c1 = p1 + (v1 * lambda);

            let b = Self::cubic(&p0, &c0, &c1, &p1);
            /*
            let mut v0_p_v1 = (v0 + v1); v0_p_v1.normalize();
            let center = *corner - v0_p_v1 * d;
            println!("center {} p0 {} p1 {} p0_2_c {} p1_2_c {}",
                     center, p0, p1, p0.distance(&center), p1.distance(&center) );

            println!("b(0.5) {} b(0.5)_2_c {}\n\n", b.point_at(F::frac(1,2)), b.point_at(F::frac(1,2)).distance(&center));
            */
            b
        }
    }

    //zz All done
}

