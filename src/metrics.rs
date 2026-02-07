//! This library provides functions that return *Metrics* for Bezier curves or the difference between
//! them and other curves or lines.
//!
use crate::{utils, BezierEval, Num};
use geo_nd::vector;

/// Maximum of the squared length of the control points
///
/// This provides a bound on the maximum distance (squared) from the origin of every point on the Bezier curve
pub fn c_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    bezier.for_each_control_point(&mut |_i, pt| {
        result = utils::max(result, vector::length_sq(pt));
    });
    result
}

/// Total of the squared length of the control points
///
/// This is guaranteed to be at least as large as c_sq
pub fn f_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    bezier.for_each_control_point(&mut |_i, pt| {
        result += vector::length_sq(pt);
    });
    result
}

/// The L2 norm between two beziers given a step dt
///
/// This is an approximation meant for testing/analysis - do not use if
/// performance is required!
///
/// L2 norm is the integral of the distance squared between each point
/// at the same t. This is a simple estiamte using summation.
pub fn dl_sq_est<F: Num, const D: usize, B: BezierEval<F, [F; D]>, B2: BezierEval<F, [F; D]>>(
    bezier: &B,
    other: &B2,
    num_steps: usize,
) -> Option<F> {
    if bezier.num_control_points() != other.num_control_points() {
        None
    } else {
        let mut total_d2 = F::zero();
        let ns: F = (num_steps as f32).into();
        for i in 0..num_steps {
            let t: F = (i as f32).into();
            let t = t / ns;
            let d2 = vector::distance_sq(&bezier.point_at(t), &other.point_at(t));
            total_d2 += d2;
        }
        Some(total_d2 / ns)
    }
}

/// The maximum difference between two beziers given a step dt
///
/// This is an approximation meant for testing/analysis - do not use if
/// performance is required!
///
/// Note that dm_sq >= l2, since l2, as can be intuited by it being 1/n.Sum(n)(d_sq(i/n)), for
/// n approaching infinity, and dm_sq is Max(d_sq) = 1/n.Sum(n)(Max(d_sq)), and d_sq() at some t
/// is always at at most Max(d_sq()).
pub fn dm_sq_est<F: Num, const D: usize, B: BezierEval<F, [F; D]>, B2: BezierEval<F, [F; D]>>(
    bezier: &B,
    other: &B2,
    num_steps: usize,
) -> Option<F> {
    if bezier.num_control_points() != other.num_control_points() {
        None
    } else {
        let mut max_d2 = F::zero();
        let ns: F = (num_steps as f32).into();
        for i in 0..num_steps {
            let t: F = (i as f32).into();
            let t = t / ns;
            let d2 = vector::distance_sq(&bezier.point_at(t), &other.point_at(t));
            if max_d2 < d2 {
                max_d2 = d2;
            }
        }
        Some(max_d2)
    }
}

/// Maximum of the squared length of the difference in control points between two Bezier curves
///
/// This requires the two Bezier curves to be of the same degree, else `None` is returned
///
/// This metric is `Max(P0[j]-P1[j])^2 * (Sum(B[i](t)))^2`, which is
/// the same as `Sum(B[i](t) * Max(P0[i]-P1[i]))^2`; this is more than
/// `Sum(Max(B[i](t) * (P0[i]-P1[i])))^2`
///
/// Hence this value is guaranteed to be greater than or equal to the `dm_sq` metric
///
/// This provides a simply computable bound on the maximum distance (squared) between the respective points on the two Bezier
/// curves for the same parameter `t`.
///
pub fn dc_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>, B2: BezierEval<F, [F; D]>>(
    bezier: &B,
    other: &B2,
) -> Option<F> {
    if bezier.num_control_points() != other.num_control_points() {
        None
    } else {
        let mut result = F::ZERO;
        bezier.for_each_control_point(&mut |i, pt| {
            result = utils::max(result, vector::distance_sq(pt, other.control_point(i)));
        });
        Some(result)
    }
}

/// Total of the squared length of the difference in control points between two Bezier curves
///
/// This requires the two Bezier curves to be of the same degree, else `None` is returned
///
/// This is clearly more than the `dc_sq`, the *maximum* of the squared difference between the
/// control points, as that is only one of the `N` control points
///
/// This, though, is clearly less than or equal to `N` times `dc_sq`, the *maximum* of the squared
/// difference between the control points.
pub fn df_sq<F: Num, const D: usize, B: BezierEval<F, [F; D]>, B2: BezierEval<F, [F; D]>>(
    bezier: &B,
    other: &B2,
) -> Option<F> {
    if bezier.num_control_points() != other.num_control_points() {
        None
    } else {
        let mut result = F::ZERO;
        bezier.for_each_control_point(&mut |i, pt| {
            result += vector::distance_sq(pt, other.control_point(i));
        });
        Some(result)
    }
}

/// Maximum of the squared length of the difference in control points between a bezier and the straight line between its endpoints
///
/// This provides a metric of how far from a straight line the curve is; all points on the
/// Bezier curve will be no further from the straight line than this distance squared
pub fn dc_sq_from_line<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    let mut iter = utils::float_iter(bezier.num_control_points());
    let (l0, l1) = bezier.endpoints();
    bezier.for_each_control_point(&mut |_i, pt| {
        let l = vector::mix(l0, l1, iter.next().unwrap());
        result = utils::max(result, vector::distance_sq(pt, &l));
    });
    result
}

/// Total of the squared length of the difference in control points between a bezier and the straight line between its endpoints
///
/// This is guaranteed to be at least as large as `dc_sq_from_line`
pub fn df_sq_from_line<F: Num, const D: usize, B: BezierEval<F, [F; D]>>(bezier: &B) -> F {
    let mut result = F::ZERO;
    let mut iter = utils::float_iter(bezier.num_control_points());
    let (l0, l1) = bezier.endpoints();
    bezier.for_each_control_point(&mut |_i, pt| {
        let l = vector::mix(l0, l1, iter.next().unwrap());
        result += vector::distance_sq(pt, &l);
    });
    result
}
