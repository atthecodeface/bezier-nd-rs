use crate::polynomial::{PolyNewtonRaphson, Polynomial};
use crate::utils;
use crate::{BezierEval, Num};

use geo_nd::vector;

/// A value of t and distance squared to it, 0<=t<=1, for which the distance between the point
/// P and the quadratic Bezier between the two points is a minimum.
///
/// The distance squared between a quadratic Bezier and a point P is a quartic equation in t (with
/// a positive coefficient of t^4)
/// and the derivative is a cubic which has at least one real root (as there is a least one minimum
/// in the distance squared); indeed, if there are three real roots then they can be ordered t0<t1<t2,
/// and the minima are at t0 and t2.
///
/// This calculates the roots, the distances to the endpoints and the minima roots, and returns the smallest
/// distance squared and t value.
// Note p(t) = u^2.p0 + 2ut.p1 + t^2.p2
// or p0 + 2(p1 - p0)t + (p0 - 2p1 + p2).t^2
// and
// p'(t) = 2(p1-p0) + 2(p0-2p1+p2)t
//
// Define p012 = p0-2*p1+p2, p01=p1-p0
//
// dsq = distance_sq to point = |p(t)-p|^2 = (p(t)-p) . (p(t)-p)
// Define p0p = p0-p
//
// d/dt (dsq) = 2 (p(t)-p) . d/dt(p(t)-p) = 2 (p(t)-p) . p'(t)
// = 0 at a minimum distance for the curve to p
//
// (p(t)-p) . p'(t) = (p0 + 2(p1 - p0)t + (p0 - 2p1 + p2).t^2 - p) . (2(p1-p0) + 2(p0-2p1+p2)t)
// = (p0p + 2p01.t + p012.t^2) . 2 . (p01 + p012.t)
// = 2(p012.p012.t^3 + (2p01.p012 + p012.p01)t^2 + (2.p01.p01 + p0p.p012)t + p0p.p01)
// = 2(p012.p012.t^3 + 3.p01.p012.t^2 + (2.p01.p01 + p0p.p012)t + p0p.p01)
// = 0 at minimum or maximum distance squared for varying t
//
// This cubic can be solved for t, leading to one or three roots
//
// d_sq is a quartic that is always positive; it will have a positive coefficient of t^4. It
// s derivative will have either one or three real roots; if there are three then they can be
// sorted (into t0, t1, t2; t0 <= t1 <= t2 and the minima of d_sq will be at t0 and t2 (but
// the closest points on 0<=t<=1 may be at t=0 or t=1)
//
pub fn bezier_quad_t_dsq_closest_to_pt<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(
    bezier: &B,
    pt: &[F; D],
) -> Option<(F, F)> {
    assert_eq!(bezier.degree(), 2);

    let (p0, p1) = bezier.endpoints();
    let pts = [*p0, *bezier.control_point(1), *p1];

    // First find the closest endpoint in case we have to bail (rather than panic)
    let d0_sq = vector::distance_sq(pt, &pts[0]);
    let d1_sq = vector::distance_sq(pt, &pts[2]);
    let mut t_min = F::ZERO;
    let mut dsq_min = d0_sq;
    if d1_sq < dsq_min {
        t_min = F::ONE;
        dsq_min = d1_sq;
    }

    // eprintln!("Cloeset endpoint {t_min}:{dsq_min}");
    let p012 = vector::sum_scaled(&pts, &[F::ONE, (-2.0_f32).into(), F::ONE]);
    let p01 = vector::sum_scaled(&pts, &[-F::ONE, F::ONE]);
    let p0p = vector::add(pts[0], pt, -F::ONE);
    let a = vector::dot(&p012, &p012);
    let d_dt_dsq_poly = [
        vector::dot(&p01, &p0p),
        (vector::dot(&p01, &p01) * (2.0_f32).into() + vector::dot(&p012, &p0p)),
        (vector::dot(&p01, &p012) * (3.0_f32).into()),
        a,
    ];
    let (opt_t0, opt_t1, opt_t2) = utils::find_root_cubic_num(d_dt_dsq_poly);
    let Some(t0) = opt_t0 else {
        panic!("Failed to find root for polynomial {d_dt_dsq_poly:?}");
        return Some((t_min, dsq_min));
    };

    let t1 = opt_t1.unwrap_or(t0);
    let t2 = opt_t2.unwrap_or(t0);
    eprintln!("Poly at t1 {t1} - {}", d_dt_dsq_poly.calc(t1));
    eprintln!("Poly at t2 {t2} - {}", d_dt_dsq_poly.calc(t2));
    // Make t0 < t1; t1 and t2 will be unordered, and t0 and t2 will be unordered
    let (t0, t1) = {
        if t0 < t1 {
            (t0, t1)
        } else {
            (t1, t0)
        }
    };
    // Make t0 < t2 ordered; t0 < t1 and t0<t2; ordering of t1 and t2 uncertain
    let (t0, t2) = {
        if t0 < t2 {
            (t0, t2)
        } else {
            (t2, t0)
        }
    };
    // Make t0 < t2 ordered; t0 < t1 and t0<t2; ordering of t1 and t2 uncertain
    let (t1, t2) = {
        if t1 < t2 {
            (t1, t2)
        } else {
            (t2, t1)
        }
    };
    // t1 is not a minimum, so we don't use that any further
    let _t1 = t1;

    // At this point t0 and t2 are minima, and t0<t2
    //
    // There are six cases:
    //          0.             1.           Minimum
    //  --------------------------------------------
    // 1: t0 t2  |              |          | at t=0
    // 2: t0     |  t2          |          | at t=0 or t2
    // 3: t0     |              |  t2      | at t=0 or t=1
    // 4:        |  t0 t2       |          | at t0 or t2
    // 5:        |  t0          |  t2      | at t0 or t=1
    // 6:        |              |  t0 t2.  } at t=1
    // eprintln!("{t0} {t1} {t2}");
    let t0_lt_zero = t0 <= F::ZERO;
    let t2_lt_zero = t2 <= F::ZERO;
    let t0_gt_one = t0 >= F::ONE;
    let t2_gt_one = t2 >= F::ONE;
    // Handle #1, #3, #6
    if t2_lt_zero || t0_gt_one || (t0_lt_zero && t2_gt_one) {
        return Some((t_min, dsq_min));
    }
    let dsq_at_t0 = vector::distance_sq(pt, &bezier.point_at(t0));
    let dsq_at_t2 = vector::distance_sq(pt, &bezier.point_at(t2));
    if !t0_lt_zero && dsq_at_t0 < dsq_min {
        dsq_min = dsq_at_t0;
        t_min = t0;
    }
    if !t2_gt_one && dsq_at_t2 < dsq_min {
        dsq_min = dsq_at_t2;
        t_min = t2;
    }
    Some((t_min, dsq_min))
}
